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

#include "ThreeCenterElectronRepulsionPrimRecHSI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsi,
                                 size_t idx_eri_0_fsi,
                                 size_t idx_eri_1_fsi,
                                 size_t idx_eri_1_gsh,
                                 size_t idx_eri_1_gsi,
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

    /// Set up components of auxilary buffer : FSI

    auto g_xxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi);

    auto g_xxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 1);

    auto g_xxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 2);

    auto g_xxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 3);

    auto g_xxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 4);

    auto g_xxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 5);

    auto g_xxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 6);

    auto g_xxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 7);

    auto g_xxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 8);

    auto g_xxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 9);

    auto g_xxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 10);

    auto g_xxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 11);

    auto g_xxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 12);

    auto g_xxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 13);

    auto g_xxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 14);

    auto g_xxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 15);

    auto g_xxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 16);

    auto g_xxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 17);

    auto g_xxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 18);

    auto g_xxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 19);

    auto g_xxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 20);

    auto g_xxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 21);

    auto g_xxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 22);

    auto g_xxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 23);

    auto g_xxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 24);

    auto g_xxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 25);

    auto g_xxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 26);

    auto g_xxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 27);

    auto g_xxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 28);

    auto g_xxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 30);

    auto g_xxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 33);

    auto g_xxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 37);

    auto g_xxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 42);

    auto g_xxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 48);

    auto g_xxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 56);

    auto g_xxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 57);

    auto g_xxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 59);

    auto g_xxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 62);

    auto g_xxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 66);

    auto g_xxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 71);

    auto g_xyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 85);

    auto g_xyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 87);

    auto g_xyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 88);

    auto g_xyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 90);

    auto g_xyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 91);

    auto g_xyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 92);

    auto g_xyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 94);

    auto g_xyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 95);

    auto g_xyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 96);

    auto g_xyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 97);

    auto g_xyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 99);

    auto g_xyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 100);

    auto g_xyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 101);

    auto g_xyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 102);

    auto g_xyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 103);

    auto g_xyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 105);

    auto g_xyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 106);

    auto g_xyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 107);

    auto g_xyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 108);

    auto g_xyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 109);

    auto g_xyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 110);

    auto g_xyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 111);

    auto g_xzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 142);

    auto g_xzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 144);

    auto g_xzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 145);

    auto g_xzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 147);

    auto g_xzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 148);

    auto g_xzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 149);

    auto g_xzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 151);

    auto g_xzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 152);

    auto g_xzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 153);

    auto g_xzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 154);

    auto g_xzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 156);

    auto g_xzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 157);

    auto g_xzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 158);

    auto g_xzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 159);

    auto g_xzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 160);

    auto g_xzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 161);

    auto g_xzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 162);

    auto g_xzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 163);

    auto g_xzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 164);

    auto g_xzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 165);

    auto g_xzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 166);

    auto g_xzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 167);

    auto g_yyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 168);

    auto g_yyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 169);

    auto g_yyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 170);

    auto g_yyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 171);

    auto g_yyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 172);

    auto g_yyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 173);

    auto g_yyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 174);

    auto g_yyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 175);

    auto g_yyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 176);

    auto g_yyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 177);

    auto g_yyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 178);

    auto g_yyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 179);

    auto g_yyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 180);

    auto g_yyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 181);

    auto g_yyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 182);

    auto g_yyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 183);

    auto g_yyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 184);

    auto g_yyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 185);

    auto g_yyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 186);

    auto g_yyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 187);

    auto g_yyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 188);

    auto g_yyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 189);

    auto g_yyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 190);

    auto g_yyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 191);

    auto g_yyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 192);

    auto g_yyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 193);

    auto g_yyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 194);

    auto g_yyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 195);

    auto g_yyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 197);

    auto g_yyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 199);

    auto g_yyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 202);

    auto g_yyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 206);

    auto g_yyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 211);

    auto g_yyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 217);

    auto g_yzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 224);

    auto g_yzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 226);

    auto g_yzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 228);

    auto g_yzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 229);

    auto g_yzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 231);

    auto g_yzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 232);

    auto g_yzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 233);

    auto g_yzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 235);

    auto g_yzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 236);

    auto g_yzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 237);

    auto g_yzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 238);

    auto g_yzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 240);

    auto g_yzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 241);

    auto g_yzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 242);

    auto g_yzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 243);

    auto g_yzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 244);

    auto g_yzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 246);

    auto g_yzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 247);

    auto g_yzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 248);

    auto g_yzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 249);

    auto g_yzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 250);

    auto g_yzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 251);

    auto g_zzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 252);

    auto g_zzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 253);

    auto g_zzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 254);

    auto g_zzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 255);

    auto g_zzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 256);

    auto g_zzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 257);

    auto g_zzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 258);

    auto g_zzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 259);

    auto g_zzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 260);

    auto g_zzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 261);

    auto g_zzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 262);

    auto g_zzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 263);

    auto g_zzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 264);

    auto g_zzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 265);

    auto g_zzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 266);

    auto g_zzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 267);

    auto g_zzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 268);

    auto g_zzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 269);

    auto g_zzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 270);

    auto g_zzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 271);

    auto g_zzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 272);

    auto g_zzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 273);

    auto g_zzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 274);

    auto g_zzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 275);

    auto g_zzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 276);

    auto g_zzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 277);

    auto g_zzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 278);

    auto g_zzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 279);

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

    auto g_xxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 28);

    auto g_xxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 30);

    auto g_xxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 33);

    auto g_xxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 37);

    auto g_xxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 42);

    auto g_xxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 48);

    auto g_xxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 56);

    auto g_xxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 57);

    auto g_xxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 59);

    auto g_xxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 62);

    auto g_xxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 66);

    auto g_xxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 71);

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

    auto g_xyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 111);

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

    auto g_xzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 161);

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

    auto g_yyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 197);

    auto g_yyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 199);

    auto g_yyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 202);

    auto g_yyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 206);

    auto g_yyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 211);

    auto g_yyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 217);

    auto g_yzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 224);

    auto g_yzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 226);

    auto g_yzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 228);

    auto g_yzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 229);

    auto g_yzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 231);

    auto g_yzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 232);

    auto g_yzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 233);

    auto g_yzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 235);

    auto g_yzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 236);

    auto g_yzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 237);

    auto g_yzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 238);

    auto g_yzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 240);

    auto g_yzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 241);

    auto g_yzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 242);

    auto g_yzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 243);

    auto g_yzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 244);

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

    /// Set up components of auxilary buffer : GSH

    auto g_xxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh);

    auto g_xxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 1);

    auto g_xxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 2);

    auto g_xxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 3);

    auto g_xxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 4);

    auto g_xxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 5);

    auto g_xxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 6);

    auto g_xxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 7);

    auto g_xxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 8);

    auto g_xxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 9);

    auto g_xxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 10);

    auto g_xxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 11);

    auto g_xxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 12);

    auto g_xxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 13);

    auto g_xxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 14);

    auto g_xxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 15);

    auto g_xxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 16);

    auto g_xxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 17);

    auto g_xxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 18);

    auto g_xxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 19);

    auto g_xxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 20);

    auto g_xxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 44);

    auto g_xxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 46);

    auto g_xxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 47);

    auto g_xxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 49);

    auto g_xxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 50);

    auto g_xxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 51);

    auto g_xxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 53);

    auto g_xxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 54);

    auto g_xxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 55);

    auto g_xxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 56);

    auto g_xxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 58);

    auto g_xxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 59);

    auto g_xxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 60);

    auto g_xxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 61);

    auto g_xxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 62);

    auto g_xxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 63);

    auto g_xxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 64);

    auto g_xxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 65);

    auto g_xxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 66);

    auto g_xxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 67);

    auto g_xxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 68);

    auto g_xxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 69);

    auto g_xxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 70);

    auto g_xxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 71);

    auto g_xxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 72);

    auto g_xxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 73);

    auto g_xxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 74);

    auto g_xxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 75);

    auto g_xxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 76);

    auto g_xxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 77);

    auto g_xxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 78);

    auto g_xxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 79);

    auto g_xxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 80);

    auto g_xxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 81);

    auto g_xxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 82);

    auto g_xxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 83);

    auto g_xxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 105);

    auto g_xxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 106);

    auto g_xxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 107);

    auto g_xxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 108);

    auto g_xxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 109);

    auto g_xxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 110);

    auto g_xxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 111);

    auto g_xxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 112);

    auto g_xxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 113);

    auto g_xxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 114);

    auto g_xxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 115);

    auto g_xxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 116);

    auto g_xxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 117);

    auto g_xxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 118);

    auto g_xxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 119);

    auto g_xxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 120);

    auto g_xxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 121);

    auto g_xxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 122);

    auto g_xxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 123);

    auto g_xxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 124);

    auto g_xxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 125);

    auto g_xyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 127);

    auto g_xyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 129);

    auto g_xyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 130);

    auto g_xyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 132);

    auto g_xyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 133);

    auto g_xyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 134);

    auto g_xyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 136);

    auto g_xyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 137);

    auto g_xyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 138);

    auto g_xyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 139);

    auto g_xyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 141);

    auto g_xyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 142);

    auto g_xyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 143);

    auto g_xyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 144);

    auto g_xyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 145);

    auto g_xzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 191);

    auto g_xzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 193);

    auto g_xzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 194);

    auto g_xzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 196);

    auto g_xzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 197);

    auto g_xzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 198);

    auto g_xzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 200);

    auto g_xzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 201);

    auto g_xzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 202);

    auto g_xzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 203);

    auto g_xzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 205);

    auto g_xzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 206);

    auto g_xzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 207);

    auto g_xzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 208);

    auto g_xzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 209);

    auto g_yyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 210);

    auto g_yyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 211);

    auto g_yyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 212);

    auto g_yyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 213);

    auto g_yyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 214);

    auto g_yyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 215);

    auto g_yyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 216);

    auto g_yyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 217);

    auto g_yyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 218);

    auto g_yyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 219);

    auto g_yyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 220);

    auto g_yyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 221);

    auto g_yyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 222);

    auto g_yyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 223);

    auto g_yyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 224);

    auto g_yyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 225);

    auto g_yyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 226);

    auto g_yyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 227);

    auto g_yyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 228);

    auto g_yyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 229);

    auto g_yyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 230);

    auto g_yyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 233);

    auto g_yyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 235);

    auto g_yyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 236);

    auto g_yyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 238);

    auto g_yyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 239);

    auto g_yyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 240);

    auto g_yyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 242);

    auto g_yyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 243);

    auto g_yyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 244);

    auto g_yyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 245);

    auto g_yyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 247);

    auto g_yyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 248);

    auto g_yyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 249);

    auto g_yyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 250);

    auto g_yyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 251);

    auto g_yyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 252);

    auto g_yyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 253);

    auto g_yyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 254);

    auto g_yyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 255);

    auto g_yyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 256);

    auto g_yyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 257);

    auto g_yyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 258);

    auto g_yyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 259);

    auto g_yyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 260);

    auto g_yyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 261);

    auto g_yyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 262);

    auto g_yyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 263);

    auto g_yyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 264);

    auto g_yyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 265);

    auto g_yyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 266);

    auto g_yyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 267);

    auto g_yyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 268);

    auto g_yyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 269);

    auto g_yyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 270);

    auto g_yyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 271);

    auto g_yyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 272);

    auto g_yzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 274);

    auto g_yzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 275);

    auto g_yzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 276);

    auto g_yzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 277);

    auto g_yzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 278);

    auto g_yzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 279);

    auto g_yzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 280);

    auto g_yzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 281);

    auto g_yzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 282);

    auto g_yzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 283);

    auto g_yzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 284);

    auto g_yzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 285);

    auto g_yzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 286);

    auto g_yzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 287);

    auto g_yzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 288);

    auto g_yzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 289);

    auto g_yzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 290);

    auto g_yzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 291);

    auto g_yzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 292);

    auto g_yzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 293);

    auto g_zzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 294);

    auto g_zzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 295);

    auto g_zzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 296);

    auto g_zzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 297);

    auto g_zzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 298);

    auto g_zzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 299);

    auto g_zzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 300);

    auto g_zzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 301);

    auto g_zzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 302);

    auto g_zzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 303);

    auto g_zzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 304);

    auto g_zzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 305);

    auto g_zzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 306);

    auto g_zzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 307);

    auto g_zzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 308);

    auto g_zzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 309);

    auto g_zzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 310);

    auto g_zzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 311);

    auto g_zzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 312);

    auto g_zzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 313);

    auto g_zzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 314);

    /// Set up components of auxilary buffer : GSI

    auto g_xxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi);

    auto g_xxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 1);

    auto g_xxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 2);

    auto g_xxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 3);

    auto g_xxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 4);

    auto g_xxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 5);

    auto g_xxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 6);

    auto g_xxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 7);

    auto g_xxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 8);

    auto g_xxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 9);

    auto g_xxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 10);

    auto g_xxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 11);

    auto g_xxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 12);

    auto g_xxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 13);

    auto g_xxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 14);

    auto g_xxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 15);

    auto g_xxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 16);

    auto g_xxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 17);

    auto g_xxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 18);

    auto g_xxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 19);

    auto g_xxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 20);

    auto g_xxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 21);

    auto g_xxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 22);

    auto g_xxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 23);

    auto g_xxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 24);

    auto g_xxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 25);

    auto g_xxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 26);

    auto g_xxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 27);

    auto g_xxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 28);

    auto g_xxxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 29);

    auto g_xxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 30);

    auto g_xxxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 31);

    auto g_xxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 33);

    auto g_xxxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 34);

    auto g_xxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 37);

    auto g_xxxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 38);

    auto g_xxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 42);

    auto g_xxxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 43);

    auto g_xxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 48);

    auto g_xxxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 49);

    auto g_xxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 56);

    auto g_xxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 57);

    auto g_xxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 58);

    auto g_xxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 59);

    auto g_xxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 60);

    auto g_xxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 61);

    auto g_xxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 62);

    auto g_xxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 63);

    auto g_xxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 64);

    auto g_xxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 65);

    auto g_xxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 66);

    auto g_xxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 67);

    auto g_xxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 68);

    auto g_xxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 69);

    auto g_xxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 70);

    auto g_xxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 71);

    auto g_xxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 72);

    auto g_xxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 73);

    auto g_xxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 74);

    auto g_xxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 75);

    auto g_xxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 76);

    auto g_xxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 78);

    auto g_xxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 79);

    auto g_xxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 80);

    auto g_xxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 81);

    auto g_xxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 82);

    auto g_xxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 83);

    auto g_xxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 84);

    auto g_xxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 85);

    auto g_xxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 86);

    auto g_xxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 87);

    auto g_xxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 88);

    auto g_xxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 89);

    auto g_xxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 90);

    auto g_xxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 91);

    auto g_xxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 92);

    auto g_xxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 93);

    auto g_xxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 94);

    auto g_xxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 95);

    auto g_xxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 96);

    auto g_xxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 97);

    auto g_xxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 98);

    auto g_xxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 99);

    auto g_xxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 100);

    auto g_xxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 101);

    auto g_xxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 102);

    auto g_xxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 103);

    auto g_xxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 104);

    auto g_xxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 105);

    auto g_xxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 106);

    auto g_xxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 107);

    auto g_xxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 108);

    auto g_xxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 109);

    auto g_xxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 110);

    auto g_xxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 111);

    auto g_xxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 140);

    auto g_xxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 141);

    auto g_xxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 142);

    auto g_xxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 143);

    auto g_xxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 144);

    auto g_xxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 145);

    auto g_xxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 146);

    auto g_xxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 147);

    auto g_xxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 148);

    auto g_xxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 149);

    auto g_xxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 150);

    auto g_xxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 151);

    auto g_xxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 152);

    auto g_xxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 153);

    auto g_xxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 154);

    auto g_xxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 155);

    auto g_xxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 156);

    auto g_xxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 157);

    auto g_xxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 158);

    auto g_xxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 159);

    auto g_xxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 160);

    auto g_xxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 161);

    auto g_xxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 162);

    auto g_xxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 163);

    auto g_xxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 164);

    auto g_xxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 165);

    auto g_xxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 166);

    auto g_xxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 167);

    auto g_xyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 168);

    auto g_xyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 169);

    auto g_xyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 171);

    auto g_xyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 172);

    auto g_xyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 174);

    auto g_xyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 175);

    auto g_xyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 176);

    auto g_xyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 178);

    auto g_xyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 179);

    auto g_xyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 180);

    auto g_xyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 181);

    auto g_xyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 183);

    auto g_xyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 184);

    auto g_xyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 185);

    auto g_xyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 186);

    auto g_xyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 187);

    auto g_xyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 189);

    auto g_xyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 190);

    auto g_xyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 191);

    auto g_xyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 192);

    auto g_xyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 193);

    auto g_xyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 194);

    auto g_xyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 195);

    auto g_xzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 252);

    auto g_xzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 254);

    auto g_xzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 256);

    auto g_xzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 257);

    auto g_xzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 259);

    auto g_xzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 260);

    auto g_xzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 261);

    auto g_xzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 263);

    auto g_xzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 264);

    auto g_xzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 265);

    auto g_xzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 266);

    auto g_xzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 268);

    auto g_xzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 269);

    auto g_xzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 270);

    auto g_xzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 271);

    auto g_xzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 272);

    auto g_xzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 273);

    auto g_xzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 274);

    auto g_xzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 275);

    auto g_xzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 276);

    auto g_xzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 277);

    auto g_xzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 278);

    auto g_xzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 279);

    auto g_yyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 280);

    auto g_yyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 281);

    auto g_yyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 282);

    auto g_yyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 283);

    auto g_yyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 284);

    auto g_yyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 285);

    auto g_yyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 286);

    auto g_yyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 287);

    auto g_yyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 288);

    auto g_yyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 289);

    auto g_yyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 290);

    auto g_yyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 291);

    auto g_yyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 292);

    auto g_yyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 293);

    auto g_yyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 294);

    auto g_yyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 295);

    auto g_yyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 296);

    auto g_yyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 297);

    auto g_yyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 298);

    auto g_yyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 299);

    auto g_yyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 300);

    auto g_yyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 301);

    auto g_yyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 302);

    auto g_yyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 303);

    auto g_yyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 304);

    auto g_yyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 305);

    auto g_yyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 306);

    auto g_yyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 307);

    auto g_yyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 309);

    auto g_yyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 310);

    auto g_yyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 311);

    auto g_yyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 312);

    auto g_yyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 313);

    auto g_yyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 314);

    auto g_yyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 315);

    auto g_yyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 316);

    auto g_yyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 317);

    auto g_yyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 318);

    auto g_yyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 319);

    auto g_yyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 320);

    auto g_yyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 321);

    auto g_yyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 322);

    auto g_yyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 323);

    auto g_yyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 324);

    auto g_yyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 325);

    auto g_yyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 326);

    auto g_yyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 327);

    auto g_yyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 328);

    auto g_yyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 329);

    auto g_yyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 330);

    auto g_yyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 331);

    auto g_yyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 332);

    auto g_yyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 333);

    auto g_yyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 334);

    auto g_yyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 335);

    auto g_yyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 336);

    auto g_yyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 337);

    auto g_yyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 338);

    auto g_yyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 339);

    auto g_yyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 340);

    auto g_yyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 341);

    auto g_yyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 342);

    auto g_yyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 343);

    auto g_yyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 344);

    auto g_yyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 345);

    auto g_yyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 346);

    auto g_yyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 347);

    auto g_yyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 348);

    auto g_yyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 349);

    auto g_yyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 350);

    auto g_yyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 351);

    auto g_yyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 352);

    auto g_yyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 353);

    auto g_yyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 354);

    auto g_yyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 355);

    auto g_yyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 356);

    auto g_yyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 357);

    auto g_yyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 358);

    auto g_yyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 359);

    auto g_yyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 360);

    auto g_yyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 361);

    auto g_yyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 362);

    auto g_yyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 363);

    auto g_yzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 364);

    auto g_yzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 365);

    auto g_yzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 366);

    auto g_yzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 367);

    auto g_yzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 368);

    auto g_yzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 369);

    auto g_yzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 370);

    auto g_yzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 371);

    auto g_yzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 372);

    auto g_yzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 373);

    auto g_yzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 374);

    auto g_yzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 375);

    auto g_yzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 376);

    auto g_yzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 377);

    auto g_yzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 378);

    auto g_yzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 379);

    auto g_yzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 380);

    auto g_yzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 381);

    auto g_yzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 382);

    auto g_yzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 383);

    auto g_yzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 384);

    auto g_yzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 385);

    auto g_yzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 386);

    auto g_yzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 387);

    auto g_yzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 388);

    auto g_yzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 389);

    auto g_yzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 390);

    auto g_yzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 391);

    auto g_zzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 392);

    auto g_zzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 393);

    auto g_zzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 394);

    auto g_zzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 395);

    auto g_zzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 396);

    auto g_zzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 397);

    auto g_zzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 398);

    auto g_zzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 399);

    auto g_zzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 400);

    auto g_zzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 401);

    auto g_zzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 402);

    auto g_zzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 403);

    auto g_zzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 404);

    auto g_zzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 405);

    auto g_zzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 406);

    auto g_zzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 407);

    auto g_zzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 408);

    auto g_zzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 409);

    auto g_zzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 410);

    auto g_zzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 411);

    auto g_zzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 412);

    auto g_zzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 413);

    auto g_zzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 414);

    auto g_zzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 415);

    auto g_zzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 416);

    auto g_zzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 417);

    auto g_zzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 418);

    auto g_zzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 419);

    /// Set up 0-28 components of targeted buffer : HSI

    auto g_xxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi);

    auto g_xxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 1);

    auto g_xxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 2);

    auto g_xxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 3);

    auto g_xxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 4);

    auto g_xxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 5);

    auto g_xxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 6);

    auto g_xxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 7);

    auto g_xxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 8);

    auto g_xxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 9);

    auto g_xxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 10);

    auto g_xxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 11);

    auto g_xxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 12);

    auto g_xxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 13);

    auto g_xxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 14);

    auto g_xxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 15);

    auto g_xxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 16);

    auto g_xxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 17);

    auto g_xxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 18);

    auto g_xxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 19);

    auto g_xxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 20);

    auto g_xxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 21);

    auto g_xxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 22);

    auto g_xxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 23);

    auto g_xxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 24);

    auto g_xxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 25);

    auto g_xxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 26);

    auto g_xxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 27);

    #pragma omp simd aligned(g_xxx_0_xxxxxx_0, g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxy_0, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxxz_0, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxyy_0, g_xxx_0_xxxxyy_1, g_xxx_0_xxxxyz_0, g_xxx_0_xxxxyz_1, g_xxx_0_xxxxzz_0, g_xxx_0_xxxxzz_1, g_xxx_0_xxxyyy_0, g_xxx_0_xxxyyy_1, g_xxx_0_xxxyyz_0, g_xxx_0_xxxyyz_1, g_xxx_0_xxxyzz_0, g_xxx_0_xxxyzz_1, g_xxx_0_xxxzzz_0, g_xxx_0_xxxzzz_1, g_xxx_0_xxyyyy_0, g_xxx_0_xxyyyy_1, g_xxx_0_xxyyyz_0, g_xxx_0_xxyyyz_1, g_xxx_0_xxyyzz_0, g_xxx_0_xxyyzz_1, g_xxx_0_xxyzzz_0, g_xxx_0_xxyzzz_1, g_xxx_0_xxzzzz_0, g_xxx_0_xxzzzz_1, g_xxx_0_xyyyyy_0, g_xxx_0_xyyyyy_1, g_xxx_0_xyyyyz_0, g_xxx_0_xyyyyz_1, g_xxx_0_xyyyzz_0, g_xxx_0_xyyyzz_1, g_xxx_0_xyyzzz_0, g_xxx_0_xyyzzz_1, g_xxx_0_xyzzzz_0, g_xxx_0_xyzzzz_1, g_xxx_0_xzzzzz_0, g_xxx_0_xzzzzz_1, g_xxx_0_yyyyyy_0, g_xxx_0_yyyyyy_1, g_xxx_0_yyyyyz_0, g_xxx_0_yyyyyz_1, g_xxx_0_yyyyzz_0, g_xxx_0_yyyyzz_1, g_xxx_0_yyyzzz_0, g_xxx_0_yyyzzz_1, g_xxx_0_yyzzzz_0, g_xxx_0_yyzzzz_1, g_xxx_0_yzzzzz_0, g_xxx_0_yzzzzz_1, g_xxx_0_zzzzzz_0, g_xxx_0_zzzzzz_1, g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxxyz_1, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxyy_1, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxxyyz_1, g_xxxx_0_xxxyz_1, g_xxxx_0_xxxyzz_1, g_xxxx_0_xxxzz_1, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxyyy_1, g_xxxx_0_xxyyyy_1, g_xxxx_0_xxyyyz_1, g_xxxx_0_xxyyz_1, g_xxxx_0_xxyyzz_1, g_xxxx_0_xxyzz_1, g_xxxx_0_xxyzzz_1, g_xxxx_0_xxzzz_1, g_xxxx_0_xxzzzz_1, g_xxxx_0_xyyyy_1, g_xxxx_0_xyyyyy_1, g_xxxx_0_xyyyyz_1, g_xxxx_0_xyyyz_1, g_xxxx_0_xyyyzz_1, g_xxxx_0_xyyzz_1, g_xxxx_0_xyyzzz_1, g_xxxx_0_xyzzz_1, g_xxxx_0_xyzzzz_1, g_xxxx_0_xzzzz_1, g_xxxx_0_xzzzzz_1, g_xxxx_0_yyyyy_1, g_xxxx_0_yyyyyy_1, g_xxxx_0_yyyyyz_1, g_xxxx_0_yyyyz_1, g_xxxx_0_yyyyzz_1, g_xxxx_0_yyyzz_1, g_xxxx_0_yyyzzz_1, g_xxxx_0_yyzzz_1, g_xxxx_0_yyzzzz_1, g_xxxx_0_yzzzz_1, g_xxxx_0_yzzzzz_1, g_xxxx_0_zzzzz_1, g_xxxx_0_zzzzzz_1, g_xxxxx_0_xxxxxx_0, g_xxxxx_0_xxxxxy_0, g_xxxxx_0_xxxxxz_0, g_xxxxx_0_xxxxyy_0, g_xxxxx_0_xxxxyz_0, g_xxxxx_0_xxxxzz_0, g_xxxxx_0_xxxyyy_0, g_xxxxx_0_xxxyyz_0, g_xxxxx_0_xxxyzz_0, g_xxxxx_0_xxxzzz_0, g_xxxxx_0_xxyyyy_0, g_xxxxx_0_xxyyyz_0, g_xxxxx_0_xxyyzz_0, g_xxxxx_0_xxyzzz_0, g_xxxxx_0_xxzzzz_0, g_xxxxx_0_xyyyyy_0, g_xxxxx_0_xyyyyz_0, g_xxxxx_0_xyyyzz_0, g_xxxxx_0_xyyzzz_0, g_xxxxx_0_xyzzzz_0, g_xxxxx_0_xzzzzz_0, g_xxxxx_0_yyyyyy_0, g_xxxxx_0_yyyyyz_0, g_xxxxx_0_yyyyzz_0, g_xxxxx_0_yyyzzz_0, g_xxxxx_0_yyzzzz_0, g_xxxxx_0_yzzzzz_0, g_xxxxx_0_zzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_xxxxxx_0[i] = 4.0 * g_xxx_0_xxxxxx_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxx_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxy_0[i] = 4.0 * g_xxx_0_xxxxxy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxz_0[i] = 4.0 * g_xxx_0_xxxxxz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyy_0[i] = 4.0 * g_xxx_0_xxxxyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyz_0[i] = 4.0 * g_xxx_0_xxxxyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxzz_0[i] = 4.0 * g_xxx_0_xxxxzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyy_0[i] = 4.0 * g_xxx_0_xxxyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyz_0[i] = 4.0 * g_xxx_0_xxxyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyzz_0[i] = 4.0 * g_xxx_0_xxxyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxzzz_0[i] = 4.0 * g_xxx_0_xxxzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyy_0[i] = 4.0 * g_xxx_0_xxyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyz_0[i] = 4.0 * g_xxx_0_xxyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyzz_0[i] = 4.0 * g_xxx_0_xxyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyzzz_0[i] = 4.0 * g_xxx_0_xxyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxzzzz_0[i] = 4.0 * g_xxx_0_xxzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyy_0[i] = 4.0 * g_xxx_0_xyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyz_0[i] = 4.0 * g_xxx_0_xyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyzz_0[i] = 4.0 * g_xxx_0_xyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyzz_1[i] * fz_be_0 + g_xxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyzzz_0[i] = 4.0 * g_xxx_0_xyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyzzz_1[i] * fz_be_0 + g_xxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyzzzz_0[i] = 4.0 * g_xxx_0_xyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyzzzz_1[i] * fz_be_0 + g_xxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xzzzzz_0[i] = 4.0 * g_xxx_0_xzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xzzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyy_0[i] = 4.0 * g_xxx_0_yyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyz_0[i] = 4.0 * g_xxx_0_yyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyzz_0[i] = 4.0 * g_xxx_0_yyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyzz_1[i] * fz_be_0 + g_xxxx_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyzzz_0[i] = 4.0 * g_xxx_0_yyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyzzz_1[i] * fz_be_0 + g_xxxx_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyzzzz_0[i] = 4.0 * g_xxx_0_yyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyzzzz_1[i] * fz_be_0 + g_xxxx_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yzzzzz_0[i] = 4.0 * g_xxx_0_yzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yzzzzz_1[i] * fz_be_0 + g_xxxx_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_zzzzzz_0[i] = 4.0 * g_xxx_0_zzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_zzzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 28-56 components of targeted buffer : HSI

    auto g_xxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 28);

    auto g_xxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 29);

    auto g_xxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 30);

    auto g_xxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 31);

    auto g_xxxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 32);

    auto g_xxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 33);

    auto g_xxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 34);

    auto g_xxxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 35);

    auto g_xxxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 36);

    auto g_xxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 37);

    auto g_xxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 38);

    auto g_xxxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 39);

    auto g_xxxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 40);

    auto g_xxxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 41);

    auto g_xxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 42);

    auto g_xxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 43);

    auto g_xxxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 44);

    auto g_xxxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 45);

    auto g_xxxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 46);

    auto g_xxxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 47);

    auto g_xxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 48);

    auto g_xxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 49);

    auto g_xxxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 50);

    auto g_xxxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 51);

    auto g_xxxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 52);

    auto g_xxxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 53);

    auto g_xxxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 54);

    auto g_xxxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 55);

    #pragma omp simd aligned(g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxxyz_1, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxyy_1, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxxyyz_1, g_xxxx_0_xxxyz_1, g_xxxx_0_xxxyzz_1, g_xxxx_0_xxxzz_1, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxyyy_1, g_xxxx_0_xxyyyy_1, g_xxxx_0_xxyyyz_1, g_xxxx_0_xxyyz_1, g_xxxx_0_xxyyzz_1, g_xxxx_0_xxyzz_1, g_xxxx_0_xxyzzz_1, g_xxxx_0_xxzzz_1, g_xxxx_0_xxzzzz_1, g_xxxx_0_xyyyy_1, g_xxxx_0_xyyyyy_1, g_xxxx_0_xyyyyz_1, g_xxxx_0_xyyyz_1, g_xxxx_0_xyyyzz_1, g_xxxx_0_xyyzz_1, g_xxxx_0_xyyzzz_1, g_xxxx_0_xyzzz_1, g_xxxx_0_xyzzzz_1, g_xxxx_0_xzzzz_1, g_xxxx_0_xzzzzz_1, g_xxxx_0_yyyyy_1, g_xxxx_0_yyyyyy_1, g_xxxx_0_yyyyyz_1, g_xxxx_0_yyyyz_1, g_xxxx_0_yyyyzz_1, g_xxxx_0_yyyzz_1, g_xxxx_0_yyyzzz_1, g_xxxx_0_yyzzz_1, g_xxxx_0_yyzzzz_1, g_xxxx_0_yzzzz_1, g_xxxx_0_yzzzzz_1, g_xxxx_0_zzzzz_1, g_xxxx_0_zzzzzz_1, g_xxxxy_0_xxxxxx_0, g_xxxxy_0_xxxxxy_0, g_xxxxy_0_xxxxxz_0, g_xxxxy_0_xxxxyy_0, g_xxxxy_0_xxxxyz_0, g_xxxxy_0_xxxxzz_0, g_xxxxy_0_xxxyyy_0, g_xxxxy_0_xxxyyz_0, g_xxxxy_0_xxxyzz_0, g_xxxxy_0_xxxzzz_0, g_xxxxy_0_xxyyyy_0, g_xxxxy_0_xxyyyz_0, g_xxxxy_0_xxyyzz_0, g_xxxxy_0_xxyzzz_0, g_xxxxy_0_xxzzzz_0, g_xxxxy_0_xyyyyy_0, g_xxxxy_0_xyyyyz_0, g_xxxxy_0_xyyyzz_0, g_xxxxy_0_xyyzzz_0, g_xxxxy_0_xyzzzz_0, g_xxxxy_0_xzzzzz_0, g_xxxxy_0_yyyyyy_0, g_xxxxy_0_yyyyyz_0, g_xxxxy_0_yyyyzz_0, g_xxxxy_0_yyyzzz_0, g_xxxxy_0_yyzzzz_0, g_xxxxy_0_yzzzzz_0, g_xxxxy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_xxxxxx_0[i] = g_xxxx_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxy_0[i] = g_xxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxz_0[i] = g_xxxx_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyy_0[i] = 2.0 * g_xxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyz_0[i] = g_xxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxzz_0[i] = g_xxxx_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyy_0[i] = 3.0 * g_xxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyz_0[i] = 2.0 * g_xxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyzz_0[i] = g_xxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxzzz_0[i] = g_xxxx_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyy_0[i] = 4.0 * g_xxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyz_0[i] = 3.0 * g_xxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyzz_0[i] = 2.0 * g_xxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyzzz_0[i] = g_xxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxzzzz_0[i] = g_xxxx_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyy_0[i] = 5.0 * g_xxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyz_0[i] = 4.0 * g_xxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyzz_0[i] = 3.0 * g_xxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyzzz_0[i] = 2.0 * g_xxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyzzzz_0[i] = g_xxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xzzzzz_0[i] = g_xxxx_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyy_0[i] = 6.0 * g_xxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyz_0[i] = 5.0 * g_xxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyzz_0[i] = 4.0 * g_xxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyzzz_0[i] = 3.0 * g_xxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyzzzz_0[i] = 2.0 * g_xxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yzzzzz_0[i] = g_xxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_zzzzzz_0[i] = g_xxxx_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 56-84 components of targeted buffer : HSI

    auto g_xxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 56);

    auto g_xxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 57);

    auto g_xxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 58);

    auto g_xxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 59);

    auto g_xxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 60);

    auto g_xxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 61);

    auto g_xxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 62);

    auto g_xxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 63);

    auto g_xxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 64);

    auto g_xxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 65);

    auto g_xxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 66);

    auto g_xxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 67);

    auto g_xxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 68);

    auto g_xxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 69);

    auto g_xxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 70);

    auto g_xxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 71);

    auto g_xxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 72);

    auto g_xxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 73);

    auto g_xxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 74);

    auto g_xxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 75);

    auto g_xxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 76);

    auto g_xxxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 77);

    auto g_xxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 78);

    auto g_xxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 79);

    auto g_xxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 80);

    auto g_xxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 81);

    auto g_xxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 82);

    auto g_xxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 83);

    #pragma omp simd aligned(g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxxyz_1, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxyy_1, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxxyyz_1, g_xxxx_0_xxxyz_1, g_xxxx_0_xxxyzz_1, g_xxxx_0_xxxzz_1, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxyyy_1, g_xxxx_0_xxyyyy_1, g_xxxx_0_xxyyyz_1, g_xxxx_0_xxyyz_1, g_xxxx_0_xxyyzz_1, g_xxxx_0_xxyzz_1, g_xxxx_0_xxyzzz_1, g_xxxx_0_xxzzz_1, g_xxxx_0_xxzzzz_1, g_xxxx_0_xyyyy_1, g_xxxx_0_xyyyyy_1, g_xxxx_0_xyyyyz_1, g_xxxx_0_xyyyz_1, g_xxxx_0_xyyyzz_1, g_xxxx_0_xyyzz_1, g_xxxx_0_xyyzzz_1, g_xxxx_0_xyzzz_1, g_xxxx_0_xyzzzz_1, g_xxxx_0_xzzzz_1, g_xxxx_0_xzzzzz_1, g_xxxx_0_yyyyy_1, g_xxxx_0_yyyyyy_1, g_xxxx_0_yyyyyz_1, g_xxxx_0_yyyyz_1, g_xxxx_0_yyyyzz_1, g_xxxx_0_yyyzz_1, g_xxxx_0_yyyzzz_1, g_xxxx_0_yyzzz_1, g_xxxx_0_yyzzzz_1, g_xxxx_0_yzzzz_1, g_xxxx_0_yzzzzz_1, g_xxxx_0_zzzzz_1, g_xxxx_0_zzzzzz_1, g_xxxxz_0_xxxxxx_0, g_xxxxz_0_xxxxxy_0, g_xxxxz_0_xxxxxz_0, g_xxxxz_0_xxxxyy_0, g_xxxxz_0_xxxxyz_0, g_xxxxz_0_xxxxzz_0, g_xxxxz_0_xxxyyy_0, g_xxxxz_0_xxxyyz_0, g_xxxxz_0_xxxyzz_0, g_xxxxz_0_xxxzzz_0, g_xxxxz_0_xxyyyy_0, g_xxxxz_0_xxyyyz_0, g_xxxxz_0_xxyyzz_0, g_xxxxz_0_xxyzzz_0, g_xxxxz_0_xxzzzz_0, g_xxxxz_0_xyyyyy_0, g_xxxxz_0_xyyyyz_0, g_xxxxz_0_xyyyzz_0, g_xxxxz_0_xyyzzz_0, g_xxxxz_0_xyzzzz_0, g_xxxxz_0_xzzzzz_0, g_xxxxz_0_yyyyyy_0, g_xxxxz_0_yyyyyz_0, g_xxxxz_0_yyyyzz_0, g_xxxxz_0_yyyzzz_0, g_xxxxz_0_yyzzzz_0, g_xxxxz_0_yzzzzz_0, g_xxxxz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_xxxxxx_0[i] = g_xxxx_0_xxxxxx_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxy_0[i] = g_xxxx_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxz_0[i] = g_xxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyy_0[i] = g_xxxx_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyz_0[i] = g_xxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxzz_0[i] = 2.0 * g_xxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyy_0[i] = g_xxxx_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyz_0[i] = g_xxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyzz_0[i] = 2.0 * g_xxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxzzz_0[i] = 3.0 * g_xxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyy_0[i] = g_xxxx_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyz_0[i] = g_xxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyzz_0[i] = 2.0 * g_xxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyzzz_0[i] = 3.0 * g_xxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxzzzz_0[i] = 4.0 * g_xxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyy_0[i] = g_xxxx_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyz_0[i] = g_xxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyzz_0[i] = 2.0 * g_xxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyzzz_0[i] = 3.0 * g_xxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyzzzz_0[i] = 4.0 * g_xxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xzzzzz_0[i] = 5.0 * g_xxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyy_0[i] = g_xxxx_0_yyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyz_0[i] = g_xxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyzz_0[i] = 2.0 * g_xxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyzzz_0[i] = 3.0 * g_xxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxx_0_yyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyzzzz_0[i] = 4.0 * g_xxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yzzzzz_0[i] = 5.0 * g_xxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_zzzzzz_0[i] = 6.0 * g_xxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxx_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 84-112 components of targeted buffer : HSI

    auto g_xxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 84);

    auto g_xxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 85);

    auto g_xxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 86);

    auto g_xxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 87);

    auto g_xxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 88);

    auto g_xxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 89);

    auto g_xxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 90);

    auto g_xxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 91);

    auto g_xxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 92);

    auto g_xxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 93);

    auto g_xxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 94);

    auto g_xxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 95);

    auto g_xxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 96);

    auto g_xxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 97);

    auto g_xxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 98);

    auto g_xxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 99);

    auto g_xxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 100);

    auto g_xxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 101);

    auto g_xxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 102);

    auto g_xxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 103);

    auto g_xxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 104);

    auto g_xxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 105);

    auto g_xxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 106);

    auto g_xxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 107);

    auto g_xxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 108);

    auto g_xxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 109);

    auto g_xxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 110);

    auto g_xxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 111);

    #pragma omp simd aligned(g_xxx_0_xxxxxx_0, g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxz_0, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxzz_0, g_xxx_0_xxxxzz_1, g_xxx_0_xxxzzz_0, g_xxx_0_xxxzzz_1, g_xxx_0_xxzzzz_0, g_xxx_0_xxzzzz_1, g_xxx_0_xzzzzz_0, g_xxx_0_xzzzzz_1, g_xxxy_0_xxxxxx_1, g_xxxy_0_xxxxxz_1, g_xxxy_0_xxxxzz_1, g_xxxy_0_xxxzzz_1, g_xxxy_0_xxzzzz_1, g_xxxy_0_xzzzzz_1, g_xxxyy_0_xxxxxx_0, g_xxxyy_0_xxxxxy_0, g_xxxyy_0_xxxxxz_0, g_xxxyy_0_xxxxyy_0, g_xxxyy_0_xxxxyz_0, g_xxxyy_0_xxxxzz_0, g_xxxyy_0_xxxyyy_0, g_xxxyy_0_xxxyyz_0, g_xxxyy_0_xxxyzz_0, g_xxxyy_0_xxxzzz_0, g_xxxyy_0_xxyyyy_0, g_xxxyy_0_xxyyyz_0, g_xxxyy_0_xxyyzz_0, g_xxxyy_0_xxyzzz_0, g_xxxyy_0_xxzzzz_0, g_xxxyy_0_xyyyyy_0, g_xxxyy_0_xyyyyz_0, g_xxxyy_0_xyyyzz_0, g_xxxyy_0_xyyzzz_0, g_xxxyy_0_xyzzzz_0, g_xxxyy_0_xzzzzz_0, g_xxxyy_0_yyyyyy_0, g_xxxyy_0_yyyyyz_0, g_xxxyy_0_yyyyzz_0, g_xxxyy_0_yyyzzz_0, g_xxxyy_0_yyzzzz_0, g_xxxyy_0_yzzzzz_0, g_xxxyy_0_zzzzzz_0, g_xxyy_0_xxxxxy_1, g_xxyy_0_xxxxy_1, g_xxyy_0_xxxxyy_1, g_xxyy_0_xxxxyz_1, g_xxyy_0_xxxyy_1, g_xxyy_0_xxxyyy_1, g_xxyy_0_xxxyyz_1, g_xxyy_0_xxxyz_1, g_xxyy_0_xxxyzz_1, g_xxyy_0_xxyyy_1, g_xxyy_0_xxyyyy_1, g_xxyy_0_xxyyyz_1, g_xxyy_0_xxyyz_1, g_xxyy_0_xxyyzz_1, g_xxyy_0_xxyzz_1, g_xxyy_0_xxyzzz_1, g_xxyy_0_xyyyy_1, g_xxyy_0_xyyyyy_1, g_xxyy_0_xyyyyz_1, g_xxyy_0_xyyyz_1, g_xxyy_0_xyyyzz_1, g_xxyy_0_xyyzz_1, g_xxyy_0_xyyzzz_1, g_xxyy_0_xyzzz_1, g_xxyy_0_xyzzzz_1, g_xxyy_0_yyyyy_1, g_xxyy_0_yyyyyy_1, g_xxyy_0_yyyyyz_1, g_xxyy_0_yyyyz_1, g_xxyy_0_yyyyzz_1, g_xxyy_0_yyyzz_1, g_xxyy_0_yyyzzz_1, g_xxyy_0_yyzzz_1, g_xxyy_0_yyzzzz_1, g_xxyy_0_yzzzz_1, g_xxyy_0_yzzzzz_1, g_xxyy_0_zzzzzz_1, g_xyy_0_xxxxxy_0, g_xyy_0_xxxxxy_1, g_xyy_0_xxxxyy_0, g_xyy_0_xxxxyy_1, g_xyy_0_xxxxyz_0, g_xyy_0_xxxxyz_1, g_xyy_0_xxxyyy_0, g_xyy_0_xxxyyy_1, g_xyy_0_xxxyyz_0, g_xyy_0_xxxyyz_1, g_xyy_0_xxxyzz_0, g_xyy_0_xxxyzz_1, g_xyy_0_xxyyyy_0, g_xyy_0_xxyyyy_1, g_xyy_0_xxyyyz_0, g_xyy_0_xxyyyz_1, g_xyy_0_xxyyzz_0, g_xyy_0_xxyyzz_1, g_xyy_0_xxyzzz_0, g_xyy_0_xxyzzz_1, g_xyy_0_xyyyyy_0, g_xyy_0_xyyyyy_1, g_xyy_0_xyyyyz_0, g_xyy_0_xyyyyz_1, g_xyy_0_xyyyzz_0, g_xyy_0_xyyyzz_1, g_xyy_0_xyyzzz_0, g_xyy_0_xyyzzz_1, g_xyy_0_xyzzzz_0, g_xyy_0_xyzzzz_1, g_xyy_0_yyyyyy_0, g_xyy_0_yyyyyy_1, g_xyy_0_yyyyyz_0, g_xyy_0_yyyyyz_1, g_xyy_0_yyyyzz_0, g_xyy_0_yyyyzz_1, g_xyy_0_yyyzzz_0, g_xyy_0_yyyzzz_1, g_xyy_0_yyzzzz_0, g_xyy_0_yyzzzz_1, g_xyy_0_yzzzzz_0, g_xyy_0_yzzzzz_1, g_xyy_0_zzzzzz_0, g_xyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyy_0_xxxxxx_0[i] = g_xxx_0_xxxxxx_0[i] * fbe_0 - g_xxx_0_xxxxxx_1[i] * fz_be_0 + g_xxxy_0_xxxxxx_1[i] * wa_y[i];

        g_xxxyy_0_xxxxxy_0[i] = 2.0 * g_xyy_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxyy_0_xxxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxz_0[i] = g_xxx_0_xxxxxz_0[i] * fbe_0 - g_xxx_0_xxxxxz_1[i] * fz_be_0 + g_xxxy_0_xxxxxz_1[i] * wa_y[i];

        g_xxxyy_0_xxxxyy_0[i] = 2.0 * g_xyy_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxyz_0[i] = 2.0 * g_xyy_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxzz_0[i] = g_xxx_0_xxxxzz_0[i] * fbe_0 - g_xxx_0_xxxxzz_1[i] * fz_be_0 + g_xxxy_0_xxxxzz_1[i] * wa_y[i];

        g_xxxyy_0_xxxyyy_0[i] = 2.0 * g_xyy_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxyyz_0[i] = 2.0 * g_xyy_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxyzz_0[i] = 2.0 * g_xyy_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxzzz_0[i] = g_xxx_0_xxxzzz_0[i] * fbe_0 - g_xxx_0_xxxzzz_1[i] * fz_be_0 + g_xxxy_0_xxxzzz_1[i] * wa_y[i];

        g_xxxyy_0_xxyyyy_0[i] = 2.0 * g_xyy_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxyyyz_0[i] = 2.0 * g_xyy_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxyyzz_0[i] = 2.0 * g_xyy_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxyzzz_0[i] = 2.0 * g_xyy_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxzzzz_0[i] = g_xxx_0_xxzzzz_0[i] * fbe_0 - g_xxx_0_xxzzzz_1[i] * fz_be_0 + g_xxxy_0_xxzzzz_1[i] * wa_y[i];

        g_xxxyy_0_xyyyyy_0[i] = 2.0 * g_xyy_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xyyyyz_0[i] = 2.0 * g_xyy_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xyyyzz_0[i] = 2.0 * g_xyy_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyzz_1[i] * fz_be_0 + g_xxyy_0_yyyzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xyyzzz_0[i] = 2.0 * g_xyy_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyzzz_1[i] * fz_be_0 + g_xxyy_0_yyzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xyzzzz_0[i] = 2.0 * g_xyy_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyzzzz_1[i] * fz_be_0 + g_xxyy_0_yzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xzzzzz_0[i] = g_xxx_0_xzzzzz_0[i] * fbe_0 - g_xxx_0_xzzzzz_1[i] * fz_be_0 + g_xxxy_0_xzzzzz_1[i] * wa_y[i];

        g_xxxyy_0_yyyyyy_0[i] = 2.0 * g_xyy_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_yyyyyz_0[i] = 2.0 * g_xyy_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_yyyyzz_0[i] = 2.0 * g_xyy_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyzz_1[i] * fz_be_0 + g_xxyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_yyyzzz_0[i] = 2.0 * g_xyy_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyzzz_1[i] * fz_be_0 + g_xxyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_yyzzzz_0[i] = 2.0 * g_xyy_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyzzzz_1[i] * fz_be_0 + g_xxyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_yzzzzz_0[i] = 2.0 * g_xyy_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yzzzzz_1[i] * fz_be_0 + g_xxyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_zzzzzz_0[i] = 2.0 * g_xyy_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_zzzzzz_1[i] * fz_be_0 + g_xxyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 112-140 components of targeted buffer : HSI

    auto g_xxxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 112);

    auto g_xxxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 113);

    auto g_xxxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 114);

    auto g_xxxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 115);

    auto g_xxxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 116);

    auto g_xxxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 117);

    auto g_xxxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 118);

    auto g_xxxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 119);

    auto g_xxxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 120);

    auto g_xxxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 121);

    auto g_xxxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 122);

    auto g_xxxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 123);

    auto g_xxxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 124);

    auto g_xxxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 125);

    auto g_xxxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 126);

    auto g_xxxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 127);

    auto g_xxxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 128);

    auto g_xxxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 129);

    auto g_xxxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 130);

    auto g_xxxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 131);

    auto g_xxxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 132);

    auto g_xxxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 133);

    auto g_xxxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 134);

    auto g_xxxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 135);

    auto g_xxxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 136);

    auto g_xxxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 137);

    auto g_xxxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 138);

    auto g_xxxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 139);

    #pragma omp simd aligned(g_xxxy_0_xxxxxy_1, g_xxxy_0_xxxxyy_1, g_xxxy_0_xxxyyy_1, g_xxxy_0_xxyyyy_1, g_xxxy_0_xyyyyy_1, g_xxxy_0_yyyyyy_1, g_xxxyz_0_xxxxxx_0, g_xxxyz_0_xxxxxy_0, g_xxxyz_0_xxxxxz_0, g_xxxyz_0_xxxxyy_0, g_xxxyz_0_xxxxyz_0, g_xxxyz_0_xxxxzz_0, g_xxxyz_0_xxxyyy_0, g_xxxyz_0_xxxyyz_0, g_xxxyz_0_xxxyzz_0, g_xxxyz_0_xxxzzz_0, g_xxxyz_0_xxyyyy_0, g_xxxyz_0_xxyyyz_0, g_xxxyz_0_xxyyzz_0, g_xxxyz_0_xxyzzz_0, g_xxxyz_0_xxzzzz_0, g_xxxyz_0_xyyyyy_0, g_xxxyz_0_xyyyyz_0, g_xxxyz_0_xyyyzz_0, g_xxxyz_0_xyyzzz_0, g_xxxyz_0_xyzzzz_0, g_xxxyz_0_xzzzzz_0, g_xxxyz_0_yyyyyy_0, g_xxxyz_0_yyyyyz_0, g_xxxyz_0_yyyyzz_0, g_xxxyz_0_yyyzzz_0, g_xxxyz_0_yyzzzz_0, g_xxxyz_0_yzzzzz_0, g_xxxyz_0_zzzzzz_0, g_xxxz_0_xxxxxx_1, g_xxxz_0_xxxxxz_1, g_xxxz_0_xxxxyz_1, g_xxxz_0_xxxxz_1, g_xxxz_0_xxxxzz_1, g_xxxz_0_xxxyyz_1, g_xxxz_0_xxxyz_1, g_xxxz_0_xxxyzz_1, g_xxxz_0_xxxzz_1, g_xxxz_0_xxxzzz_1, g_xxxz_0_xxyyyz_1, g_xxxz_0_xxyyz_1, g_xxxz_0_xxyyzz_1, g_xxxz_0_xxyzz_1, g_xxxz_0_xxyzzz_1, g_xxxz_0_xxzzz_1, g_xxxz_0_xxzzzz_1, g_xxxz_0_xyyyyz_1, g_xxxz_0_xyyyz_1, g_xxxz_0_xyyyzz_1, g_xxxz_0_xyyzz_1, g_xxxz_0_xyyzzz_1, g_xxxz_0_xyzzz_1, g_xxxz_0_xyzzzz_1, g_xxxz_0_xzzzz_1, g_xxxz_0_xzzzzz_1, g_xxxz_0_yyyyyz_1, g_xxxz_0_yyyyz_1, g_xxxz_0_yyyyzz_1, g_xxxz_0_yyyzz_1, g_xxxz_0_yyyzzz_1, g_xxxz_0_yyzzz_1, g_xxxz_0_yyzzzz_1, g_xxxz_0_yzzzz_1, g_xxxz_0_yzzzzz_1, g_xxxz_0_zzzzz_1, g_xxxz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyz_0_xxxxxx_0[i] = g_xxxz_0_xxxxxx_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxy_0[i] = g_xxxy_0_xxxxxy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxxz_0[i] = g_xxxz_0_xxxxxz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxyy_0[i] = g_xxxy_0_xxxxyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxyz_0[i] = g_xxxz_0_xxxxz_1[i] * fi_acd_0 + g_xxxz_0_xxxxyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxzz_0[i] = g_xxxz_0_xxxxzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyyy_0[i] = g_xxxy_0_xxxyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxyyz_0[i] = 2.0 * g_xxxz_0_xxxyz_1[i] * fi_acd_0 + g_xxxz_0_xxxyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyzz_0[i] = g_xxxz_0_xxxzz_1[i] * fi_acd_0 + g_xxxz_0_xxxyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxzzz_0[i] = g_xxxz_0_xxxzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyyy_0[i] = g_xxxy_0_xxyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxyyyz_0[i] = 3.0 * g_xxxz_0_xxyyz_1[i] * fi_acd_0 + g_xxxz_0_xxyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyzz_0[i] = 2.0 * g_xxxz_0_xxyzz_1[i] * fi_acd_0 + g_xxxz_0_xxyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyzzz_0[i] = g_xxxz_0_xxzzz_1[i] * fi_acd_0 + g_xxxz_0_xxyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxzzzz_0[i] = g_xxxz_0_xxzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyyy_0[i] = g_xxxy_0_xyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xyyyyz_0[i] = 4.0 * g_xxxz_0_xyyyz_1[i] * fi_acd_0 + g_xxxz_0_xyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyzz_0[i] = 3.0 * g_xxxz_0_xyyzz_1[i] * fi_acd_0 + g_xxxz_0_xyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyzzz_0[i] = 2.0 * g_xxxz_0_xyzzz_1[i] * fi_acd_0 + g_xxxz_0_xyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyzzzz_0[i] = g_xxxz_0_xzzzz_1[i] * fi_acd_0 + g_xxxz_0_xyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xzzzzz_0[i] = g_xxxz_0_xzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyyy_0[i] = g_xxxy_0_yyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_yyyyyz_0[i] = 5.0 * g_xxxz_0_yyyyz_1[i] * fi_acd_0 + g_xxxz_0_yyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyzz_0[i] = 4.0 * g_xxxz_0_yyyzz_1[i] * fi_acd_0 + g_xxxz_0_yyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyzzz_0[i] = 3.0 * g_xxxz_0_yyzzz_1[i] * fi_acd_0 + g_xxxz_0_yyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyzzzz_0[i] = 2.0 * g_xxxz_0_yzzzz_1[i] * fi_acd_0 + g_xxxz_0_yyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yzzzzz_0[i] = g_xxxz_0_zzzzz_1[i] * fi_acd_0 + g_xxxz_0_yzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_zzzzzz_0[i] = g_xxxz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 140-168 components of targeted buffer : HSI

    auto g_xxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 140);

    auto g_xxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 141);

    auto g_xxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 142);

    auto g_xxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 143);

    auto g_xxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 144);

    auto g_xxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 145);

    auto g_xxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 146);

    auto g_xxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 147);

    auto g_xxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 148);

    auto g_xxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 149);

    auto g_xxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 150);

    auto g_xxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 151);

    auto g_xxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 152);

    auto g_xxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 153);

    auto g_xxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 154);

    auto g_xxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 155);

    auto g_xxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 156);

    auto g_xxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 157);

    auto g_xxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 158);

    auto g_xxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 159);

    auto g_xxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 160);

    auto g_xxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 161);

    auto g_xxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 162);

    auto g_xxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 163);

    auto g_xxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 164);

    auto g_xxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 165);

    auto g_xxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 166);

    auto g_xxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 167);

    #pragma omp simd aligned(g_xxx_0_xxxxxx_0, g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxy_0, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxyy_0, g_xxx_0_xxxxyy_1, g_xxx_0_xxxyyy_0, g_xxx_0_xxxyyy_1, g_xxx_0_xxyyyy_0, g_xxx_0_xxyyyy_1, g_xxx_0_xyyyyy_0, g_xxx_0_xyyyyy_1, g_xxxz_0_xxxxxx_1, g_xxxz_0_xxxxxy_1, g_xxxz_0_xxxxyy_1, g_xxxz_0_xxxyyy_1, g_xxxz_0_xxyyyy_1, g_xxxz_0_xyyyyy_1, g_xxxzz_0_xxxxxx_0, g_xxxzz_0_xxxxxy_0, g_xxxzz_0_xxxxxz_0, g_xxxzz_0_xxxxyy_0, g_xxxzz_0_xxxxyz_0, g_xxxzz_0_xxxxzz_0, g_xxxzz_0_xxxyyy_0, g_xxxzz_0_xxxyyz_0, g_xxxzz_0_xxxyzz_0, g_xxxzz_0_xxxzzz_0, g_xxxzz_0_xxyyyy_0, g_xxxzz_0_xxyyyz_0, g_xxxzz_0_xxyyzz_0, g_xxxzz_0_xxyzzz_0, g_xxxzz_0_xxzzzz_0, g_xxxzz_0_xyyyyy_0, g_xxxzz_0_xyyyyz_0, g_xxxzz_0_xyyyzz_0, g_xxxzz_0_xyyzzz_0, g_xxxzz_0_xyzzzz_0, g_xxxzz_0_xzzzzz_0, g_xxxzz_0_yyyyyy_0, g_xxxzz_0_yyyyyz_0, g_xxxzz_0_yyyyzz_0, g_xxxzz_0_yyyzzz_0, g_xxxzz_0_yyzzzz_0, g_xxxzz_0_yzzzzz_0, g_xxxzz_0_zzzzzz_0, g_xxzz_0_xxxxxz_1, g_xxzz_0_xxxxyz_1, g_xxzz_0_xxxxz_1, g_xxzz_0_xxxxzz_1, g_xxzz_0_xxxyyz_1, g_xxzz_0_xxxyz_1, g_xxzz_0_xxxyzz_1, g_xxzz_0_xxxzz_1, g_xxzz_0_xxxzzz_1, g_xxzz_0_xxyyyz_1, g_xxzz_0_xxyyz_1, g_xxzz_0_xxyyzz_1, g_xxzz_0_xxyzz_1, g_xxzz_0_xxyzzz_1, g_xxzz_0_xxzzz_1, g_xxzz_0_xxzzzz_1, g_xxzz_0_xyyyyz_1, g_xxzz_0_xyyyz_1, g_xxzz_0_xyyyzz_1, g_xxzz_0_xyyzz_1, g_xxzz_0_xyyzzz_1, g_xxzz_0_xyzzz_1, g_xxzz_0_xyzzzz_1, g_xxzz_0_xzzzz_1, g_xxzz_0_xzzzzz_1, g_xxzz_0_yyyyyy_1, g_xxzz_0_yyyyyz_1, g_xxzz_0_yyyyz_1, g_xxzz_0_yyyyzz_1, g_xxzz_0_yyyzz_1, g_xxzz_0_yyyzzz_1, g_xxzz_0_yyzzz_1, g_xxzz_0_yyzzzz_1, g_xxzz_0_yzzzz_1, g_xxzz_0_yzzzzz_1, g_xxzz_0_zzzzz_1, g_xxzz_0_zzzzzz_1, g_xzz_0_xxxxxz_0, g_xzz_0_xxxxxz_1, g_xzz_0_xxxxyz_0, g_xzz_0_xxxxyz_1, g_xzz_0_xxxxzz_0, g_xzz_0_xxxxzz_1, g_xzz_0_xxxyyz_0, g_xzz_0_xxxyyz_1, g_xzz_0_xxxyzz_0, g_xzz_0_xxxyzz_1, g_xzz_0_xxxzzz_0, g_xzz_0_xxxzzz_1, g_xzz_0_xxyyyz_0, g_xzz_0_xxyyyz_1, g_xzz_0_xxyyzz_0, g_xzz_0_xxyyzz_1, g_xzz_0_xxyzzz_0, g_xzz_0_xxyzzz_1, g_xzz_0_xxzzzz_0, g_xzz_0_xxzzzz_1, g_xzz_0_xyyyyz_0, g_xzz_0_xyyyyz_1, g_xzz_0_xyyyzz_0, g_xzz_0_xyyyzz_1, g_xzz_0_xyyzzz_0, g_xzz_0_xyyzzz_1, g_xzz_0_xyzzzz_0, g_xzz_0_xyzzzz_1, g_xzz_0_xzzzzz_0, g_xzz_0_xzzzzz_1, g_xzz_0_yyyyyy_0, g_xzz_0_yyyyyy_1, g_xzz_0_yyyyyz_0, g_xzz_0_yyyyyz_1, g_xzz_0_yyyyzz_0, g_xzz_0_yyyyzz_1, g_xzz_0_yyyzzz_0, g_xzz_0_yyyzzz_1, g_xzz_0_yyzzzz_0, g_xzz_0_yyzzzz_1, g_xzz_0_yzzzzz_0, g_xzz_0_yzzzzz_1, g_xzz_0_zzzzzz_0, g_xzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzz_0_xxxxxx_0[i] = g_xxx_0_xxxxxx_0[i] * fbe_0 - g_xxx_0_xxxxxx_1[i] * fz_be_0 + g_xxxz_0_xxxxxx_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxy_0[i] = g_xxx_0_xxxxxy_0[i] * fbe_0 - g_xxx_0_xxxxxy_1[i] * fz_be_0 + g_xxxz_0_xxxxxy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxz_0[i] = 2.0 * g_xzz_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxzz_0_xxxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxyy_0[i] = g_xxx_0_xxxxyy_0[i] * fbe_0 - g_xxx_0_xxxxyy_1[i] * fz_be_0 + g_xxxz_0_xxxxyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxyz_0[i] = 2.0 * g_xzz_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxzz_0[i] = 2.0 * g_xzz_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyyy_0[i] = g_xxx_0_xxxyyy_0[i] * fbe_0 - g_xxx_0_xxxyyy_1[i] * fz_be_0 + g_xxxz_0_xxxyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxyyz_0[i] = 2.0 * g_xzz_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyzz_0[i] = 2.0 * g_xzz_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxzzz_0[i] = 2.0 * g_xzz_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyyy_0[i] = g_xxx_0_xxyyyy_0[i] * fbe_0 - g_xxx_0_xxyyyy_1[i] * fz_be_0 + g_xxxz_0_xxyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxyyyz_0[i] = 2.0 * g_xzz_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyzz_0[i] = 2.0 * g_xzz_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyzzz_0[i] = 2.0 * g_xzz_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxzzzz_0[i] = 2.0 * g_xzz_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyyy_0[i] = g_xxx_0_xyyyyy_0[i] * fbe_0 - g_xxx_0_xyyyyy_1[i] * fz_be_0 + g_xxxz_0_xyyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xyyyyz_0[i] = 2.0 * g_xzz_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyzz_0[i] = 2.0 * g_xzz_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyzz_1[i] * fz_be_0 + g_xxzz_0_yyyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyzzz_0[i] = 2.0 * g_xzz_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyzzz_1[i] * fz_be_0 + g_xxzz_0_yyzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyzzzz_0[i] = 2.0 * g_xzz_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyzzzz_1[i] * fz_be_0 + g_xxzz_0_yzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xzzzzz_0[i] = 2.0 * g_xzz_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xzzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzzz_1[i] * fi_acd_0 + g_xxzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyy_0[i] = 2.0 * g_xzz_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyy_1[i] * fz_be_0 + g_xxzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyz_0[i] = 2.0 * g_xzz_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyzz_0[i] = 2.0 * g_xzz_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyzz_1[i] * fz_be_0 + g_xxzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyzzz_0[i] = 2.0 * g_xzz_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyzzz_1[i] * fz_be_0 + g_xxzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyzzzz_0[i] = 2.0 * g_xzz_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyzzzz_1[i] * fz_be_0 + g_xxzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yzzzzz_0[i] = 2.0 * g_xzz_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yzzzzz_1[i] * fz_be_0 + g_xxzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_zzzzzz_0[i] = 2.0 * g_xzz_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_zzzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 168-196 components of targeted buffer : HSI

    auto g_xxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 168);

    auto g_xxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 169);

    auto g_xxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 170);

    auto g_xxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 171);

    auto g_xxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 172);

    auto g_xxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 173);

    auto g_xxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 174);

    auto g_xxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 175);

    auto g_xxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 176);

    auto g_xxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 177);

    auto g_xxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 178);

    auto g_xxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 179);

    auto g_xxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 180);

    auto g_xxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 181);

    auto g_xxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 182);

    auto g_xxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 183);

    auto g_xxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 184);

    auto g_xxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 185);

    auto g_xxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 186);

    auto g_xxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 187);

    auto g_xxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 188);

    auto g_xxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 189);

    auto g_xxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 190);

    auto g_xxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 191);

    auto g_xxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 192);

    auto g_xxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 193);

    auto g_xxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 194);

    auto g_xxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 195);

    #pragma omp simd aligned(g_xxy_0_xxxxxx_0, g_xxy_0_xxxxxx_1, g_xxy_0_xxxxxz_0, g_xxy_0_xxxxxz_1, g_xxy_0_xxxxzz_0, g_xxy_0_xxxxzz_1, g_xxy_0_xxxzzz_0, g_xxy_0_xxxzzz_1, g_xxy_0_xxzzzz_0, g_xxy_0_xxzzzz_1, g_xxy_0_xzzzzz_0, g_xxy_0_xzzzzz_1, g_xxyy_0_xxxxxx_1, g_xxyy_0_xxxxxz_1, g_xxyy_0_xxxxzz_1, g_xxyy_0_xxxzzz_1, g_xxyy_0_xxzzzz_1, g_xxyy_0_xzzzzz_1, g_xxyyy_0_xxxxxx_0, g_xxyyy_0_xxxxxy_0, g_xxyyy_0_xxxxxz_0, g_xxyyy_0_xxxxyy_0, g_xxyyy_0_xxxxyz_0, g_xxyyy_0_xxxxzz_0, g_xxyyy_0_xxxyyy_0, g_xxyyy_0_xxxyyz_0, g_xxyyy_0_xxxyzz_0, g_xxyyy_0_xxxzzz_0, g_xxyyy_0_xxyyyy_0, g_xxyyy_0_xxyyyz_0, g_xxyyy_0_xxyyzz_0, g_xxyyy_0_xxyzzz_0, g_xxyyy_0_xxzzzz_0, g_xxyyy_0_xyyyyy_0, g_xxyyy_0_xyyyyz_0, g_xxyyy_0_xyyyzz_0, g_xxyyy_0_xyyzzz_0, g_xxyyy_0_xyzzzz_0, g_xxyyy_0_xzzzzz_0, g_xxyyy_0_yyyyyy_0, g_xxyyy_0_yyyyyz_0, g_xxyyy_0_yyyyzz_0, g_xxyyy_0_yyyzzz_0, g_xxyyy_0_yyzzzz_0, g_xxyyy_0_yzzzzz_0, g_xxyyy_0_zzzzzz_0, g_xyyy_0_xxxxxy_1, g_xyyy_0_xxxxy_1, g_xyyy_0_xxxxyy_1, g_xyyy_0_xxxxyz_1, g_xyyy_0_xxxyy_1, g_xyyy_0_xxxyyy_1, g_xyyy_0_xxxyyz_1, g_xyyy_0_xxxyz_1, g_xyyy_0_xxxyzz_1, g_xyyy_0_xxyyy_1, g_xyyy_0_xxyyyy_1, g_xyyy_0_xxyyyz_1, g_xyyy_0_xxyyz_1, g_xyyy_0_xxyyzz_1, g_xyyy_0_xxyzz_1, g_xyyy_0_xxyzzz_1, g_xyyy_0_xyyyy_1, g_xyyy_0_xyyyyy_1, g_xyyy_0_xyyyyz_1, g_xyyy_0_xyyyz_1, g_xyyy_0_xyyyzz_1, g_xyyy_0_xyyzz_1, g_xyyy_0_xyyzzz_1, g_xyyy_0_xyzzz_1, g_xyyy_0_xyzzzz_1, g_xyyy_0_yyyyy_1, g_xyyy_0_yyyyyy_1, g_xyyy_0_yyyyyz_1, g_xyyy_0_yyyyz_1, g_xyyy_0_yyyyzz_1, g_xyyy_0_yyyzz_1, g_xyyy_0_yyyzzz_1, g_xyyy_0_yyzzz_1, g_xyyy_0_yyzzzz_1, g_xyyy_0_yzzzz_1, g_xyyy_0_yzzzzz_1, g_xyyy_0_zzzzzz_1, g_yyy_0_xxxxxy_0, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxyy_0, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyz_0, g_yyy_0_xxxxyz_1, g_yyy_0_xxxyyy_0, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyz_0, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyzz_0, g_yyy_0_xxxyzz_1, g_yyy_0_xxyyyy_0, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyz_0, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyzz_0, g_yyy_0_xxyyzz_1, g_yyy_0_xxyzzz_0, g_yyy_0_xxyzzz_1, g_yyy_0_xyyyyy_0, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyz_0, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyzz_0, g_yyy_0_xyyyzz_1, g_yyy_0_xyyzzz_0, g_yyy_0_xyyzzz_1, g_yyy_0_xyzzzz_0, g_yyy_0_xyzzzz_1, g_yyy_0_yyyyyy_0, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyz_0, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyzz_0, g_yyy_0_yyyyzz_1, g_yyy_0_yyyzzz_0, g_yyy_0_yyyzzz_1, g_yyy_0_yyzzzz_0, g_yyy_0_yyzzzz_1, g_yyy_0_yzzzzz_0, g_yyy_0_yzzzzz_1, g_yyy_0_zzzzzz_0, g_yyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyy_0_xxxxxx_0[i] = 2.0 * g_xxy_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxx_1[i] * fz_be_0 + g_xxyy_0_xxxxxx_1[i] * wa_y[i];

        g_xxyyy_0_xxxxxy_0[i] = g_yyy_0_xxxxxy_0[i] * fbe_0 - g_yyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyyy_0_xxxxy_1[i] * fi_acd_0 + g_xyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxz_0[i] = 2.0 * g_xxy_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxz_1[i] * fz_be_0 + g_xxyy_0_xxxxxz_1[i] * wa_y[i];

        g_xxyyy_0_xxxxyy_0[i] = g_yyy_0_xxxxyy_0[i] * fbe_0 - g_yyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyy_1[i] * fi_acd_0 + g_xyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxyz_0[i] = g_yyy_0_xxxxyz_0[i] * fbe_0 - g_yyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyz_1[i] * fi_acd_0 + g_xyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxzz_0[i] = 2.0 * g_xxy_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxzz_1[i] * fz_be_0 + g_xxyy_0_xxxxzz_1[i] * wa_y[i];

        g_xxyyy_0_xxxyyy_0[i] = g_yyy_0_xxxyyy_0[i] * fbe_0 - g_yyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyy_1[i] * fi_acd_0 + g_xyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxyyz_0[i] = g_yyy_0_xxxyyz_0[i] * fbe_0 - g_yyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyz_1[i] * fi_acd_0 + g_xyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxyzz_0[i] = g_yyy_0_xxxyzz_0[i] * fbe_0 - g_yyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyzz_1[i] * fi_acd_0 + g_xyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxzzz_0[i] = 2.0 * g_xxy_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxzzz_1[i] * fz_be_0 + g_xxyy_0_xxxzzz_1[i] * wa_y[i];

        g_xxyyy_0_xxyyyy_0[i] = g_yyy_0_xxyyyy_0[i] * fbe_0 - g_yyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyy_1[i] * fi_acd_0 + g_xyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxyyyz_0[i] = g_yyy_0_xxyyyz_0[i] * fbe_0 - g_yyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyz_1[i] * fi_acd_0 + g_xyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxyyzz_0[i] = g_yyy_0_xxyyzz_0[i] * fbe_0 - g_yyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyzz_1[i] * fi_acd_0 + g_xyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxyzzz_0[i] = g_yyy_0_xxyzzz_0[i] * fbe_0 - g_yyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyzzz_1[i] * fi_acd_0 + g_xyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxzzzz_0[i] = 2.0 * g_xxy_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxzzzz_1[i] * fz_be_0 + g_xxyy_0_xxzzzz_1[i] * wa_y[i];

        g_xxyyy_0_xyyyyy_0[i] = g_yyy_0_xyyyyy_0[i] * fbe_0 - g_yyy_0_xyyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyyy_1[i] * fi_acd_0 + g_xyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xyyyyz_0[i] = g_yyy_0_xyyyyz_0[i] * fbe_0 - g_yyy_0_xyyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyyz_1[i] * fi_acd_0 + g_xyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xyyyzz_0[i] = g_yyy_0_xyyyzz_0[i] * fbe_0 - g_yyy_0_xyyyzz_1[i] * fz_be_0 + g_xyyy_0_yyyzz_1[i] * fi_acd_0 + g_xyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xyyzzz_0[i] = g_yyy_0_xyyzzz_0[i] * fbe_0 - g_yyy_0_xyyzzz_1[i] * fz_be_0 + g_xyyy_0_yyzzz_1[i] * fi_acd_0 + g_xyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xyzzzz_0[i] = g_yyy_0_xyzzzz_0[i] * fbe_0 - g_yyy_0_xyzzzz_1[i] * fz_be_0 + g_xyyy_0_yzzzz_1[i] * fi_acd_0 + g_xyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xzzzzz_0[i] = 2.0 * g_xxy_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xzzzzz_1[i] * fz_be_0 + g_xxyy_0_xzzzzz_1[i] * wa_y[i];

        g_xxyyy_0_yyyyyy_0[i] = g_yyy_0_yyyyyy_0[i] * fbe_0 - g_yyy_0_yyyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_yyyyyz_0[i] = g_yyy_0_yyyyyz_0[i] * fbe_0 - g_yyy_0_yyyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_yyyyzz_0[i] = g_yyy_0_yyyyzz_0[i] * fbe_0 - g_yyy_0_yyyyzz_1[i] * fz_be_0 + g_xyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_yyyzzz_0[i] = g_yyy_0_yyyzzz_0[i] * fbe_0 - g_yyy_0_yyyzzz_1[i] * fz_be_0 + g_xyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_yyzzzz_0[i] = g_yyy_0_yyzzzz_0[i] * fbe_0 - g_yyy_0_yyzzzz_1[i] * fz_be_0 + g_xyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_yzzzzz_0[i] = g_yyy_0_yzzzzz_0[i] * fbe_0 - g_yyy_0_yzzzzz_1[i] * fz_be_0 + g_xyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_zzzzzz_0[i] = g_yyy_0_zzzzzz_0[i] * fbe_0 - g_yyy_0_zzzzzz_1[i] * fz_be_0 + g_xyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 196-224 components of targeted buffer : HSI

    auto g_xxyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 196);

    auto g_xxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 197);

    auto g_xxyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 198);

    auto g_xxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 199);

    auto g_xxyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 200);

    auto g_xxyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 201);

    auto g_xxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 202);

    auto g_xxyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 203);

    auto g_xxyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 204);

    auto g_xxyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 205);

    auto g_xxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 206);

    auto g_xxyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 207);

    auto g_xxyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 208);

    auto g_xxyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 209);

    auto g_xxyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 210);

    auto g_xxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 211);

    auto g_xxyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 212);

    auto g_xxyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 213);

    auto g_xxyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 214);

    auto g_xxyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 215);

    auto g_xxyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 216);

    auto g_xxyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 217);

    auto g_xxyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 218);

    auto g_xxyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 219);

    auto g_xxyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 220);

    auto g_xxyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 221);

    auto g_xxyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 222);

    auto g_xxyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 223);

    #pragma omp simd aligned(g_xxyy_0_xxxxx_1, g_xxyy_0_xxxxxx_1, g_xxyy_0_xxxxxy_1, g_xxyy_0_xxxxxz_1, g_xxyy_0_xxxxy_1, g_xxyy_0_xxxxyy_1, g_xxyy_0_xxxxyz_1, g_xxyy_0_xxxxz_1, g_xxyy_0_xxxxzz_1, g_xxyy_0_xxxyy_1, g_xxyy_0_xxxyyy_1, g_xxyy_0_xxxyyz_1, g_xxyy_0_xxxyz_1, g_xxyy_0_xxxyzz_1, g_xxyy_0_xxxzz_1, g_xxyy_0_xxxzzz_1, g_xxyy_0_xxyyy_1, g_xxyy_0_xxyyyy_1, g_xxyy_0_xxyyyz_1, g_xxyy_0_xxyyz_1, g_xxyy_0_xxyyzz_1, g_xxyy_0_xxyzz_1, g_xxyy_0_xxyzzz_1, g_xxyy_0_xxzzz_1, g_xxyy_0_xxzzzz_1, g_xxyy_0_xyyyy_1, g_xxyy_0_xyyyyy_1, g_xxyy_0_xyyyyz_1, g_xxyy_0_xyyyz_1, g_xxyy_0_xyyyzz_1, g_xxyy_0_xyyzz_1, g_xxyy_0_xyyzzz_1, g_xxyy_0_xyzzz_1, g_xxyy_0_xyzzzz_1, g_xxyy_0_xzzzz_1, g_xxyy_0_xzzzzz_1, g_xxyy_0_yyyyy_1, g_xxyy_0_yyyyyy_1, g_xxyy_0_yyyyyz_1, g_xxyy_0_yyyyz_1, g_xxyy_0_yyyyzz_1, g_xxyy_0_yyyzz_1, g_xxyy_0_yyyzzz_1, g_xxyy_0_yyzzz_1, g_xxyy_0_yyzzzz_1, g_xxyy_0_yzzzz_1, g_xxyy_0_yzzzzz_1, g_xxyy_0_zzzzz_1, g_xxyy_0_zzzzzz_1, g_xxyyz_0_xxxxxx_0, g_xxyyz_0_xxxxxy_0, g_xxyyz_0_xxxxxz_0, g_xxyyz_0_xxxxyy_0, g_xxyyz_0_xxxxyz_0, g_xxyyz_0_xxxxzz_0, g_xxyyz_0_xxxyyy_0, g_xxyyz_0_xxxyyz_0, g_xxyyz_0_xxxyzz_0, g_xxyyz_0_xxxzzz_0, g_xxyyz_0_xxyyyy_0, g_xxyyz_0_xxyyyz_0, g_xxyyz_0_xxyyzz_0, g_xxyyz_0_xxyzzz_0, g_xxyyz_0_xxzzzz_0, g_xxyyz_0_xyyyyy_0, g_xxyyz_0_xyyyyz_0, g_xxyyz_0_xyyyzz_0, g_xxyyz_0_xyyzzz_0, g_xxyyz_0_xyzzzz_0, g_xxyyz_0_xzzzzz_0, g_xxyyz_0_yyyyyy_0, g_xxyyz_0_yyyyyz_0, g_xxyyz_0_yyyyzz_0, g_xxyyz_0_yyyzzz_0, g_xxyyz_0_yyzzzz_0, g_xxyyz_0_yzzzzz_0, g_xxyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_xxxxxx_0[i] = g_xxyy_0_xxxxxx_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxy_0[i] = g_xxyy_0_xxxxxy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxz_0[i] = g_xxyy_0_xxxxx_1[i] * fi_acd_0 + g_xxyy_0_xxxxxz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyy_0[i] = g_xxyy_0_xxxxyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyz_0[i] = g_xxyy_0_xxxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxxyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxzz_0[i] = 2.0 * g_xxyy_0_xxxxz_1[i] * fi_acd_0 + g_xxyy_0_xxxxzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyy_0[i] = g_xxyy_0_xxxyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyz_0[i] = g_xxyy_0_xxxyy_1[i] * fi_acd_0 + g_xxyy_0_xxxyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyzz_0[i] = 2.0 * g_xxyy_0_xxxyz_1[i] * fi_acd_0 + g_xxyy_0_xxxyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxzzz_0[i] = 3.0 * g_xxyy_0_xxxzz_1[i] * fi_acd_0 + g_xxyy_0_xxxzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyy_0[i] = g_xxyy_0_xxyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyz_0[i] = g_xxyy_0_xxyyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyzz_0[i] = 2.0 * g_xxyy_0_xxyyz_1[i] * fi_acd_0 + g_xxyy_0_xxyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyzzz_0[i] = 3.0 * g_xxyy_0_xxyzz_1[i] * fi_acd_0 + g_xxyy_0_xxyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxzzzz_0[i] = 4.0 * g_xxyy_0_xxzzz_1[i] * fi_acd_0 + g_xxyy_0_xxzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyy_0[i] = g_xxyy_0_xyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyz_0[i] = g_xxyy_0_xyyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyzz_0[i] = 2.0 * g_xxyy_0_xyyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyzzz_0[i] = 3.0 * g_xxyy_0_xyyzz_1[i] * fi_acd_0 + g_xxyy_0_xyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyzzzz_0[i] = 4.0 * g_xxyy_0_xyzzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xzzzzz_0[i] = 5.0 * g_xxyy_0_xzzzz_1[i] * fi_acd_0 + g_xxyy_0_xzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyy_0[i] = g_xxyy_0_yyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyz_0[i] = g_xxyy_0_yyyyy_1[i] * fi_acd_0 + g_xxyy_0_yyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyzz_0[i] = 2.0 * g_xxyy_0_yyyyz_1[i] * fi_acd_0 + g_xxyy_0_yyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyzzz_0[i] = 3.0 * g_xxyy_0_yyyzz_1[i] * fi_acd_0 + g_xxyy_0_yyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyzzzz_0[i] = 4.0 * g_xxyy_0_yyzzz_1[i] * fi_acd_0 + g_xxyy_0_yyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yzzzzz_0[i] = 5.0 * g_xxyy_0_yzzzz_1[i] * fi_acd_0 + g_xxyy_0_yzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_zzzzzz_0[i] = 6.0 * g_xxyy_0_zzzzz_1[i] * fi_acd_0 + g_xxyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 224-252 components of targeted buffer : HSI

    auto g_xxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 224);

    auto g_xxyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 225);

    auto g_xxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 226);

    auto g_xxyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 227);

    auto g_xxyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 228);

    auto g_xxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 229);

    auto g_xxyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 230);

    auto g_xxyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 231);

    auto g_xxyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 232);

    auto g_xxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 233);

    auto g_xxyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 234);

    auto g_xxyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 235);

    auto g_xxyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 236);

    auto g_xxyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 237);

    auto g_xxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 238);

    auto g_xxyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 239);

    auto g_xxyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 240);

    auto g_xxyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 241);

    auto g_xxyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 242);

    auto g_xxyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 243);

    auto g_xxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 244);

    auto g_xxyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 245);

    auto g_xxyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 246);

    auto g_xxyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 247);

    auto g_xxyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 248);

    auto g_xxyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 249);

    auto g_xxyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 250);

    auto g_xxyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 251);

    #pragma omp simd aligned(g_xxyzz_0_xxxxxx_0, g_xxyzz_0_xxxxxy_0, g_xxyzz_0_xxxxxz_0, g_xxyzz_0_xxxxyy_0, g_xxyzz_0_xxxxyz_0, g_xxyzz_0_xxxxzz_0, g_xxyzz_0_xxxyyy_0, g_xxyzz_0_xxxyyz_0, g_xxyzz_0_xxxyzz_0, g_xxyzz_0_xxxzzz_0, g_xxyzz_0_xxyyyy_0, g_xxyzz_0_xxyyyz_0, g_xxyzz_0_xxyyzz_0, g_xxyzz_0_xxyzzz_0, g_xxyzz_0_xxzzzz_0, g_xxyzz_0_xyyyyy_0, g_xxyzz_0_xyyyyz_0, g_xxyzz_0_xyyyzz_0, g_xxyzz_0_xyyzzz_0, g_xxyzz_0_xyzzzz_0, g_xxyzz_0_xzzzzz_0, g_xxyzz_0_yyyyyy_0, g_xxyzz_0_yyyyyz_0, g_xxyzz_0_yyyyzz_0, g_xxyzz_0_yyyzzz_0, g_xxyzz_0_yyzzzz_0, g_xxyzz_0_yzzzzz_0, g_xxyzz_0_zzzzzz_0, g_xxzz_0_xxxxx_1, g_xxzz_0_xxxxxx_1, g_xxzz_0_xxxxxy_1, g_xxzz_0_xxxxxz_1, g_xxzz_0_xxxxy_1, g_xxzz_0_xxxxyy_1, g_xxzz_0_xxxxyz_1, g_xxzz_0_xxxxz_1, g_xxzz_0_xxxxzz_1, g_xxzz_0_xxxyy_1, g_xxzz_0_xxxyyy_1, g_xxzz_0_xxxyyz_1, g_xxzz_0_xxxyz_1, g_xxzz_0_xxxyzz_1, g_xxzz_0_xxxzz_1, g_xxzz_0_xxxzzz_1, g_xxzz_0_xxyyy_1, g_xxzz_0_xxyyyy_1, g_xxzz_0_xxyyyz_1, g_xxzz_0_xxyyz_1, g_xxzz_0_xxyyzz_1, g_xxzz_0_xxyzz_1, g_xxzz_0_xxyzzz_1, g_xxzz_0_xxzzz_1, g_xxzz_0_xxzzzz_1, g_xxzz_0_xyyyy_1, g_xxzz_0_xyyyyy_1, g_xxzz_0_xyyyyz_1, g_xxzz_0_xyyyz_1, g_xxzz_0_xyyyzz_1, g_xxzz_0_xyyzz_1, g_xxzz_0_xyyzzz_1, g_xxzz_0_xyzzz_1, g_xxzz_0_xyzzzz_1, g_xxzz_0_xzzzz_1, g_xxzz_0_xzzzzz_1, g_xxzz_0_yyyyy_1, g_xxzz_0_yyyyyy_1, g_xxzz_0_yyyyyz_1, g_xxzz_0_yyyyz_1, g_xxzz_0_yyyyzz_1, g_xxzz_0_yyyzz_1, g_xxzz_0_yyyzzz_1, g_xxzz_0_yyzzz_1, g_xxzz_0_yyzzzz_1, g_xxzz_0_yzzzz_1, g_xxzz_0_yzzzzz_1, g_xxzz_0_zzzzz_1, g_xxzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_xxxxxx_0[i] = g_xxzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxy_0[i] = g_xxzz_0_xxxxx_1[i] * fi_acd_0 + g_xxzz_0_xxxxxy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxz_0[i] = g_xxzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyy_0[i] = 2.0 * g_xxzz_0_xxxxy_1[i] * fi_acd_0 + g_xxzz_0_xxxxyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyz_0[i] = g_xxzz_0_xxxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxzz_0[i] = g_xxzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyy_0[i] = 3.0 * g_xxzz_0_xxxyy_1[i] * fi_acd_0 + g_xxzz_0_xxxyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyz_0[i] = 2.0 * g_xxzz_0_xxxyz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyzz_0[i] = g_xxzz_0_xxxzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxzzz_0[i] = g_xxzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyy_0[i] = 4.0 * g_xxzz_0_xxyyy_1[i] * fi_acd_0 + g_xxzz_0_xxyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyz_0[i] = 3.0 * g_xxzz_0_xxyyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyzz_0[i] = 2.0 * g_xxzz_0_xxyzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyzzz_0[i] = g_xxzz_0_xxzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxzzzz_0[i] = g_xxzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyy_0[i] = 5.0 * g_xxzz_0_xyyyy_1[i] * fi_acd_0 + g_xxzz_0_xyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyz_0[i] = 4.0 * g_xxzz_0_xyyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyzz_0[i] = 3.0 * g_xxzz_0_xyyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyzzz_0[i] = 2.0 * g_xxzz_0_xyzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyzzzz_0[i] = g_xxzz_0_xzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xzzzzz_0[i] = g_xxzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyy_0[i] = 6.0 * g_xxzz_0_yyyyy_1[i] * fi_acd_0 + g_xxzz_0_yyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyz_0[i] = 5.0 * g_xxzz_0_yyyyz_1[i] * fi_acd_0 + g_xxzz_0_yyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyzz_0[i] = 4.0 * g_xxzz_0_yyyzz_1[i] * fi_acd_0 + g_xxzz_0_yyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyzzz_0[i] = 3.0 * g_xxzz_0_yyzzz_1[i] * fi_acd_0 + g_xxzz_0_yyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyzzzz_0[i] = 2.0 * g_xxzz_0_yzzzz_1[i] * fi_acd_0 + g_xxzz_0_yyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yzzzzz_0[i] = g_xxzz_0_zzzzz_1[i] * fi_acd_0 + g_xxzz_0_yzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_zzzzzz_0[i] = g_xxzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 252-280 components of targeted buffer : HSI

    auto g_xxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 252);

    auto g_xxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 253);

    auto g_xxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 254);

    auto g_xxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 255);

    auto g_xxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 256);

    auto g_xxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 257);

    auto g_xxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 258);

    auto g_xxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 259);

    auto g_xxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 260);

    auto g_xxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 261);

    auto g_xxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 262);

    auto g_xxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 263);

    auto g_xxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 264);

    auto g_xxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 265);

    auto g_xxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 266);

    auto g_xxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 267);

    auto g_xxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 268);

    auto g_xxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 269);

    auto g_xxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 270);

    auto g_xxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 271);

    auto g_xxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 272);

    auto g_xxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 273);

    auto g_xxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 274);

    auto g_xxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 275);

    auto g_xxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 276);

    auto g_xxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 277);

    auto g_xxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 278);

    auto g_xxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 279);

    #pragma omp simd aligned(g_xxz_0_xxxxxx_0, g_xxz_0_xxxxxx_1, g_xxz_0_xxxxxy_0, g_xxz_0_xxxxxy_1, g_xxz_0_xxxxyy_0, g_xxz_0_xxxxyy_1, g_xxz_0_xxxyyy_0, g_xxz_0_xxxyyy_1, g_xxz_0_xxyyyy_0, g_xxz_0_xxyyyy_1, g_xxz_0_xyyyyy_0, g_xxz_0_xyyyyy_1, g_xxzz_0_xxxxxx_1, g_xxzz_0_xxxxxy_1, g_xxzz_0_xxxxyy_1, g_xxzz_0_xxxyyy_1, g_xxzz_0_xxyyyy_1, g_xxzz_0_xyyyyy_1, g_xxzzz_0_xxxxxx_0, g_xxzzz_0_xxxxxy_0, g_xxzzz_0_xxxxxz_0, g_xxzzz_0_xxxxyy_0, g_xxzzz_0_xxxxyz_0, g_xxzzz_0_xxxxzz_0, g_xxzzz_0_xxxyyy_0, g_xxzzz_0_xxxyyz_0, g_xxzzz_0_xxxyzz_0, g_xxzzz_0_xxxzzz_0, g_xxzzz_0_xxyyyy_0, g_xxzzz_0_xxyyyz_0, g_xxzzz_0_xxyyzz_0, g_xxzzz_0_xxyzzz_0, g_xxzzz_0_xxzzzz_0, g_xxzzz_0_xyyyyy_0, g_xxzzz_0_xyyyyz_0, g_xxzzz_0_xyyyzz_0, g_xxzzz_0_xyyzzz_0, g_xxzzz_0_xyzzzz_0, g_xxzzz_0_xzzzzz_0, g_xxzzz_0_yyyyyy_0, g_xxzzz_0_yyyyyz_0, g_xxzzz_0_yyyyzz_0, g_xxzzz_0_yyyzzz_0, g_xxzzz_0_yyzzzz_0, g_xxzzz_0_yzzzzz_0, g_xxzzz_0_zzzzzz_0, g_xzzz_0_xxxxxz_1, g_xzzz_0_xxxxyz_1, g_xzzz_0_xxxxz_1, g_xzzz_0_xxxxzz_1, g_xzzz_0_xxxyyz_1, g_xzzz_0_xxxyz_1, g_xzzz_0_xxxyzz_1, g_xzzz_0_xxxzz_1, g_xzzz_0_xxxzzz_1, g_xzzz_0_xxyyyz_1, g_xzzz_0_xxyyz_1, g_xzzz_0_xxyyzz_1, g_xzzz_0_xxyzz_1, g_xzzz_0_xxyzzz_1, g_xzzz_0_xxzzz_1, g_xzzz_0_xxzzzz_1, g_xzzz_0_xyyyyz_1, g_xzzz_0_xyyyz_1, g_xzzz_0_xyyyzz_1, g_xzzz_0_xyyzz_1, g_xzzz_0_xyyzzz_1, g_xzzz_0_xyzzz_1, g_xzzz_0_xyzzzz_1, g_xzzz_0_xzzzz_1, g_xzzz_0_xzzzzz_1, g_xzzz_0_yyyyyy_1, g_xzzz_0_yyyyyz_1, g_xzzz_0_yyyyz_1, g_xzzz_0_yyyyzz_1, g_xzzz_0_yyyzz_1, g_xzzz_0_yyyzzz_1, g_xzzz_0_yyzzz_1, g_xzzz_0_yyzzzz_1, g_xzzz_0_yzzzz_1, g_xzzz_0_yzzzzz_1, g_xzzz_0_zzzzz_1, g_xzzz_0_zzzzzz_1, g_zzz_0_xxxxxz_0, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxyz_0, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxzz_0, g_zzz_0_xxxxzz_1, g_zzz_0_xxxyyz_0, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyzz_0, g_zzz_0_xxxyzz_1, g_zzz_0_xxxzzz_0, g_zzz_0_xxxzzz_1, g_zzz_0_xxyyyz_0, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyzz_0, g_zzz_0_xxyyzz_1, g_zzz_0_xxyzzz_0, g_zzz_0_xxyzzz_1, g_zzz_0_xxzzzz_0, g_zzz_0_xxzzzz_1, g_zzz_0_xyyyyz_0, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyzz_0, g_zzz_0_xyyyzz_1, g_zzz_0_xyyzzz_0, g_zzz_0_xyyzzz_1, g_zzz_0_xyzzzz_0, g_zzz_0_xyzzzz_1, g_zzz_0_xzzzzz_0, g_zzz_0_xzzzzz_1, g_zzz_0_yyyyyy_0, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyz_0, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyzz_0, g_zzz_0_yyyyzz_1, g_zzz_0_yyyzzz_0, g_zzz_0_yyyzzz_1, g_zzz_0_yyzzzz_0, g_zzz_0_yyzzzz_1, g_zzz_0_yzzzzz_0, g_zzz_0_yzzzzz_1, g_zzz_0_zzzzzz_0, g_zzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzz_0_xxxxxx_0[i] = 2.0 * g_xxz_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxx_1[i] * fz_be_0 + g_xxzz_0_xxxxxx_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxy_0[i] = 2.0 * g_xxz_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxy_1[i] * fz_be_0 + g_xxzz_0_xxxxxy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxz_0[i] = g_zzz_0_xxxxxz_0[i] * fbe_0 - g_zzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzzz_0_xxxxz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxyy_0[i] = 2.0 * g_xxz_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxyy_1[i] * fz_be_0 + g_xxzz_0_xxxxyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxyz_0[i] = g_zzz_0_xxxxyz_0[i] * fbe_0 - g_zzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxyz_1[i] * fi_acd_0 + g_xzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxzz_0[i] = g_zzz_0_xxxxzz_0[i] * fbe_0 - g_zzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyyy_0[i] = 2.0 * g_xxz_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxyyy_1[i] * fz_be_0 + g_xxzz_0_xxxyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxyyz_0[i] = g_zzz_0_xxxyyz_0[i] * fbe_0 - g_zzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyyz_1[i] * fi_acd_0 + g_xzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyzz_0[i] = g_zzz_0_xxxyzz_0[i] * fbe_0 - g_zzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyzz_1[i] * fi_acd_0 + g_xzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxzzz_0[i] = g_zzz_0_xxxzzz_0[i] * fbe_0 - g_zzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyyy_0[i] = 2.0 * g_xxz_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxyyyy_1[i] * fz_be_0 + g_xxzz_0_xxyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxyyyz_0[i] = g_zzz_0_xxyyyz_0[i] * fbe_0 - g_zzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyyz_1[i] * fi_acd_0 + g_xzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyzz_0[i] = g_zzz_0_xxyyzz_0[i] * fbe_0 - g_zzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyzz_1[i] * fi_acd_0 + g_xzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyzzz_0[i] = g_zzz_0_xxyzzz_0[i] * fbe_0 - g_zzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyzzz_1[i] * fi_acd_0 + g_xzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxzzzz_0[i] = g_zzz_0_xxzzzz_0[i] * fbe_0 - g_zzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyyy_0[i] = 2.0 * g_xxz_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xyyyyy_1[i] * fz_be_0 + g_xxzz_0_xyyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xyyyyz_0[i] = g_zzz_0_xyyyyz_0[i] * fbe_0 - g_zzz_0_xyyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyyz_1[i] * fi_acd_0 + g_xzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyzz_0[i] = g_zzz_0_xyyyzz_0[i] * fbe_0 - g_zzz_0_xyyyzz_1[i] * fz_be_0 + g_xzzz_0_yyyzz_1[i] * fi_acd_0 + g_xzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyzzz_0[i] = g_zzz_0_xyyzzz_0[i] * fbe_0 - g_zzz_0_xyyzzz_1[i] * fz_be_0 + g_xzzz_0_yyzzz_1[i] * fi_acd_0 + g_xzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyzzzz_0[i] = g_zzz_0_xyzzzz_0[i] * fbe_0 - g_zzz_0_xyzzzz_1[i] * fz_be_0 + g_xzzz_0_yzzzz_1[i] * fi_acd_0 + g_xzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xzzzzz_0[i] = g_zzz_0_xzzzzz_0[i] * fbe_0 - g_zzz_0_xzzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzzz_1[i] * fi_acd_0 + g_xzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyy_0[i] = g_zzz_0_yyyyyy_0[i] * fbe_0 - g_zzz_0_yyyyyy_1[i] * fz_be_0 + g_xzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyz_0[i] = g_zzz_0_yyyyyz_0[i] * fbe_0 - g_zzz_0_yyyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyzz_0[i] = g_zzz_0_yyyyzz_0[i] * fbe_0 - g_zzz_0_yyyyzz_1[i] * fz_be_0 + g_xzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyzzz_0[i] = g_zzz_0_yyyzzz_0[i] * fbe_0 - g_zzz_0_yyyzzz_1[i] * fz_be_0 + g_xzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyzzzz_0[i] = g_zzz_0_yyzzzz_0[i] * fbe_0 - g_zzz_0_yyzzzz_1[i] * fz_be_0 + g_xzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yzzzzz_0[i] = g_zzz_0_yzzzzz_0[i] * fbe_0 - g_zzz_0_yzzzzz_1[i] * fz_be_0 + g_xzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_zzzzzz_0[i] = g_zzz_0_zzzzzz_0[i] * fbe_0 - g_zzz_0_zzzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 280-308 components of targeted buffer : HSI

    auto g_xyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 280);

    auto g_xyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 281);

    auto g_xyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 282);

    auto g_xyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 283);

    auto g_xyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 284);

    auto g_xyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 285);

    auto g_xyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 286);

    auto g_xyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 287);

    auto g_xyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 288);

    auto g_xyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 289);

    auto g_xyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 290);

    auto g_xyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 291);

    auto g_xyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 292);

    auto g_xyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 293);

    auto g_xyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 294);

    auto g_xyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 295);

    auto g_xyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 296);

    auto g_xyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 297);

    auto g_xyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 298);

    auto g_xyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 299);

    auto g_xyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 300);

    auto g_xyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 301);

    auto g_xyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 302);

    auto g_xyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 303);

    auto g_xyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 304);

    auto g_xyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 305);

    auto g_xyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 306);

    auto g_xyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 307);

    #pragma omp simd aligned(g_xyyyy_0_xxxxxx_0, g_xyyyy_0_xxxxxy_0, g_xyyyy_0_xxxxxz_0, g_xyyyy_0_xxxxyy_0, g_xyyyy_0_xxxxyz_0, g_xyyyy_0_xxxxzz_0, g_xyyyy_0_xxxyyy_0, g_xyyyy_0_xxxyyz_0, g_xyyyy_0_xxxyzz_0, g_xyyyy_0_xxxzzz_0, g_xyyyy_0_xxyyyy_0, g_xyyyy_0_xxyyyz_0, g_xyyyy_0_xxyyzz_0, g_xyyyy_0_xxyzzz_0, g_xyyyy_0_xxzzzz_0, g_xyyyy_0_xyyyyy_0, g_xyyyy_0_xyyyyz_0, g_xyyyy_0_xyyyzz_0, g_xyyyy_0_xyyzzz_0, g_xyyyy_0_xyzzzz_0, g_xyyyy_0_xzzzzz_0, g_xyyyy_0_yyyyyy_0, g_xyyyy_0_yyyyyz_0, g_xyyyy_0_yyyyzz_0, g_xyyyy_0_yyyzzz_0, g_xyyyy_0_yyzzzz_0, g_xyyyy_0_yzzzzz_0, g_xyyyy_0_zzzzzz_0, g_yyyy_0_xxxxx_1, g_yyyy_0_xxxxxx_1, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxxz_1, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxxz_1, g_yyyy_0_xxxxzz_1, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyz_1, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxxzz_1, g_yyyy_0_xxxzzz_1, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyzz_1, g_yyyy_0_xxyzzz_1, g_yyyy_0_xxzzz_1, g_yyyy_0_xxzzzz_1, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyzz_1, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyzzz_1, g_yyyy_0_xyzzzz_1, g_yyyy_0_xzzzz_1, g_yyyy_0_xzzzzz_1, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyzz_1, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyzzz_1, g_yyyy_0_yyzzzz_1, g_yyyy_0_yzzzz_1, g_yyyy_0_yzzzzz_1, g_yyyy_0_zzzzz_1, g_yyyy_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_xxxxxx_0[i] = 6.0 * g_yyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxx_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxy_0[i] = 5.0 * g_yyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxz_0[i] = 5.0 * g_yyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyy_0[i] = 4.0 * g_yyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyz_0[i] = 4.0 * g_yyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxzz_0[i] = 4.0 * g_yyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyy_0[i] = 3.0 * g_yyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyz_0[i] = 3.0 * g_yyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyzz_0[i] = 3.0 * g_yyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxzzz_0[i] = 3.0 * g_yyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyy_0[i] = 2.0 * g_yyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyz_0[i] = 2.0 * g_yyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyzz_0[i] = 2.0 * g_yyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyzzz_0[i] = 2.0 * g_yyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxzzzz_0[i] = 2.0 * g_yyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyy_0[i] = g_yyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyz_0[i] = g_yyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyzz_0[i] = g_yyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyzzz_0[i] = g_yyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyzzzz_0[i] = g_yyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xzzzzz_0[i] = g_yyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyy_0[i] = g_yyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyz_0[i] = g_yyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyzz_0[i] = g_yyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyzzz_0[i] = g_yyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyzzzz_0[i] = g_yyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yzzzzz_0[i] = g_yyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_zzzzzz_0[i] = g_yyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 308-336 components of targeted buffer : HSI

    auto g_xyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 308);

    auto g_xyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 309);

    auto g_xyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 310);

    auto g_xyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 311);

    auto g_xyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 312);

    auto g_xyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 313);

    auto g_xyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 314);

    auto g_xyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 315);

    auto g_xyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 316);

    auto g_xyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 317);

    auto g_xyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 318);

    auto g_xyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 319);

    auto g_xyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 320);

    auto g_xyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 321);

    auto g_xyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 322);

    auto g_xyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 323);

    auto g_xyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 324);

    auto g_xyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 325);

    auto g_xyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 326);

    auto g_xyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 327);

    auto g_xyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 328);

    auto g_xyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 329);

    auto g_xyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 330);

    auto g_xyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 331);

    auto g_xyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 332);

    auto g_xyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 333);

    auto g_xyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 334);

    auto g_xyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 335);

    #pragma omp simd aligned(g_xyyy_0_xxxxxx_1, g_xyyy_0_xxxxxy_1, g_xyyy_0_xxxxyy_1, g_xyyy_0_xxxyyy_1, g_xyyy_0_xxyyyy_1, g_xyyy_0_xyyyyy_1, g_xyyyz_0_xxxxxx_0, g_xyyyz_0_xxxxxy_0, g_xyyyz_0_xxxxxz_0, g_xyyyz_0_xxxxyy_0, g_xyyyz_0_xxxxyz_0, g_xyyyz_0_xxxxzz_0, g_xyyyz_0_xxxyyy_0, g_xyyyz_0_xxxyyz_0, g_xyyyz_0_xxxyzz_0, g_xyyyz_0_xxxzzz_0, g_xyyyz_0_xxyyyy_0, g_xyyyz_0_xxyyyz_0, g_xyyyz_0_xxyyzz_0, g_xyyyz_0_xxyzzz_0, g_xyyyz_0_xxzzzz_0, g_xyyyz_0_xyyyyy_0, g_xyyyz_0_xyyyyz_0, g_xyyyz_0_xyyyzz_0, g_xyyyz_0_xyyzzz_0, g_xyyyz_0_xyzzzz_0, g_xyyyz_0_xzzzzz_0, g_xyyyz_0_yyyyyy_0, g_xyyyz_0_yyyyyz_0, g_xyyyz_0_yyyyzz_0, g_xyyyz_0_yyyzzz_0, g_xyyyz_0_yyzzzz_0, g_xyyyz_0_yzzzzz_0, g_xyyyz_0_zzzzzz_0, g_yyyz_0_xxxxxz_1, g_yyyz_0_xxxxyz_1, g_yyyz_0_xxxxz_1, g_yyyz_0_xxxxzz_1, g_yyyz_0_xxxyyz_1, g_yyyz_0_xxxyz_1, g_yyyz_0_xxxyzz_1, g_yyyz_0_xxxzz_1, g_yyyz_0_xxxzzz_1, g_yyyz_0_xxyyyz_1, g_yyyz_0_xxyyz_1, g_yyyz_0_xxyyzz_1, g_yyyz_0_xxyzz_1, g_yyyz_0_xxyzzz_1, g_yyyz_0_xxzzz_1, g_yyyz_0_xxzzzz_1, g_yyyz_0_xyyyyz_1, g_yyyz_0_xyyyz_1, g_yyyz_0_xyyyzz_1, g_yyyz_0_xyyzz_1, g_yyyz_0_xyyzzz_1, g_yyyz_0_xyzzz_1, g_yyyz_0_xyzzzz_1, g_yyyz_0_xzzzz_1, g_yyyz_0_xzzzzz_1, g_yyyz_0_yyyyyy_1, g_yyyz_0_yyyyyz_1, g_yyyz_0_yyyyz_1, g_yyyz_0_yyyyzz_1, g_yyyz_0_yyyzz_1, g_yyyz_0_yyyzzz_1, g_yyyz_0_yyzzz_1, g_yyyz_0_yyzzzz_1, g_yyyz_0_yzzzz_1, g_yyyz_0_yzzzzz_1, g_yyyz_0_zzzzz_1, g_yyyz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyz_0_xxxxxx_0[i] = g_xyyy_0_xxxxxx_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxy_0[i] = g_xyyy_0_xxxxxy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxz_0[i] = 5.0 * g_yyyz_0_xxxxz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxyy_0[i] = g_xyyy_0_xxxxyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxyz_0[i] = 4.0 * g_yyyz_0_xxxyz_1[i] * fi_acd_0 + g_yyyz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxzz_0[i] = 4.0 * g_yyyz_0_xxxzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyyy_0[i] = g_xyyy_0_xxxyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxyyz_0[i] = 3.0 * g_yyyz_0_xxyyz_1[i] * fi_acd_0 + g_yyyz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyzz_0[i] = 3.0 * g_yyyz_0_xxyzz_1[i] * fi_acd_0 + g_yyyz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxzzz_0[i] = 3.0 * g_yyyz_0_xxzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyyy_0[i] = g_xyyy_0_xxyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxyyyz_0[i] = 2.0 * g_yyyz_0_xyyyz_1[i] * fi_acd_0 + g_yyyz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyzz_0[i] = 2.0 * g_yyyz_0_xyyzz_1[i] * fi_acd_0 + g_yyyz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyzzz_0[i] = 2.0 * g_yyyz_0_xyzzz_1[i] * fi_acd_0 + g_yyyz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxzzzz_0[i] = 2.0 * g_yyyz_0_xzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyyy_0[i] = g_xyyy_0_xyyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xyyyyz_0[i] = g_yyyz_0_yyyyz_1[i] * fi_acd_0 + g_yyyz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyzz_0[i] = g_yyyz_0_yyyzz_1[i] * fi_acd_0 + g_yyyz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyzzz_0[i] = g_yyyz_0_yyzzz_1[i] * fi_acd_0 + g_yyyz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyzzzz_0[i] = g_yyyz_0_yzzzz_1[i] * fi_acd_0 + g_yyyz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xzzzzz_0[i] = g_yyyz_0_zzzzz_1[i] * fi_acd_0 + g_yyyz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyy_0[i] = g_yyyz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyz_0[i] = g_yyyz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyzz_0[i] = g_yyyz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyzzz_0[i] = g_yyyz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyzzzz_0[i] = g_yyyz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yzzzzz_0[i] = g_yyyz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_zzzzzz_0[i] = g_yyyz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 336-364 components of targeted buffer : HSI

    auto g_xyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 336);

    auto g_xyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 337);

    auto g_xyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 338);

    auto g_xyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 339);

    auto g_xyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 340);

    auto g_xyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 341);

    auto g_xyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 342);

    auto g_xyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 343);

    auto g_xyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 344);

    auto g_xyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 345);

    auto g_xyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 346);

    auto g_xyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 347);

    auto g_xyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 348);

    auto g_xyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 349);

    auto g_xyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 350);

    auto g_xyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 351);

    auto g_xyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 352);

    auto g_xyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 353);

    auto g_xyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 354);

    auto g_xyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 355);

    auto g_xyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 356);

    auto g_xyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 357);

    auto g_xyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 358);

    auto g_xyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 359);

    auto g_xyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 360);

    auto g_xyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 361);

    auto g_xyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 362);

    auto g_xyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 363);

    #pragma omp simd aligned(g_xyyzz_0_xxxxxx_0, g_xyyzz_0_xxxxxy_0, g_xyyzz_0_xxxxxz_0, g_xyyzz_0_xxxxyy_0, g_xyyzz_0_xxxxyz_0, g_xyyzz_0_xxxxzz_0, g_xyyzz_0_xxxyyy_0, g_xyyzz_0_xxxyyz_0, g_xyyzz_0_xxxyzz_0, g_xyyzz_0_xxxzzz_0, g_xyyzz_0_xxyyyy_0, g_xyyzz_0_xxyyyz_0, g_xyyzz_0_xxyyzz_0, g_xyyzz_0_xxyzzz_0, g_xyyzz_0_xxzzzz_0, g_xyyzz_0_xyyyyy_0, g_xyyzz_0_xyyyyz_0, g_xyyzz_0_xyyyzz_0, g_xyyzz_0_xyyzzz_0, g_xyyzz_0_xyzzzz_0, g_xyyzz_0_xzzzzz_0, g_xyyzz_0_yyyyyy_0, g_xyyzz_0_yyyyyz_0, g_xyyzz_0_yyyyzz_0, g_xyyzz_0_yyyzzz_0, g_xyyzz_0_yyzzzz_0, g_xyyzz_0_yzzzzz_0, g_xyyzz_0_zzzzzz_0, g_yyzz_0_xxxxx_1, g_yyzz_0_xxxxxx_1, g_yyzz_0_xxxxxy_1, g_yyzz_0_xxxxxz_1, g_yyzz_0_xxxxy_1, g_yyzz_0_xxxxyy_1, g_yyzz_0_xxxxyz_1, g_yyzz_0_xxxxz_1, g_yyzz_0_xxxxzz_1, g_yyzz_0_xxxyy_1, g_yyzz_0_xxxyyy_1, g_yyzz_0_xxxyyz_1, g_yyzz_0_xxxyz_1, g_yyzz_0_xxxyzz_1, g_yyzz_0_xxxzz_1, g_yyzz_0_xxxzzz_1, g_yyzz_0_xxyyy_1, g_yyzz_0_xxyyyy_1, g_yyzz_0_xxyyyz_1, g_yyzz_0_xxyyz_1, g_yyzz_0_xxyyzz_1, g_yyzz_0_xxyzz_1, g_yyzz_0_xxyzzz_1, g_yyzz_0_xxzzz_1, g_yyzz_0_xxzzzz_1, g_yyzz_0_xyyyy_1, g_yyzz_0_xyyyyy_1, g_yyzz_0_xyyyyz_1, g_yyzz_0_xyyyz_1, g_yyzz_0_xyyyzz_1, g_yyzz_0_xyyzz_1, g_yyzz_0_xyyzzz_1, g_yyzz_0_xyzzz_1, g_yyzz_0_xyzzzz_1, g_yyzz_0_xzzzz_1, g_yyzz_0_xzzzzz_1, g_yyzz_0_yyyyy_1, g_yyzz_0_yyyyyy_1, g_yyzz_0_yyyyyz_1, g_yyzz_0_yyyyz_1, g_yyzz_0_yyyyzz_1, g_yyzz_0_yyyzz_1, g_yyzz_0_yyyzzz_1, g_yyzz_0_yyzzz_1, g_yyzz_0_yyzzzz_1, g_yyzz_0_yzzzz_1, g_yyzz_0_yzzzzz_1, g_yyzz_0_zzzzz_1, g_yyzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_xxxxxx_0[i] = 6.0 * g_yyzz_0_xxxxx_1[i] * fi_acd_0 + g_yyzz_0_xxxxxx_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxy_0[i] = 5.0 * g_yyzz_0_xxxxy_1[i] * fi_acd_0 + g_yyzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxz_0[i] = 5.0 * g_yyzz_0_xxxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyy_0[i] = 4.0 * g_yyzz_0_xxxyy_1[i] * fi_acd_0 + g_yyzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyz_0[i] = 4.0 * g_yyzz_0_xxxyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxzz_0[i] = 4.0 * g_yyzz_0_xxxzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyy_0[i] = 3.0 * g_yyzz_0_xxyyy_1[i] * fi_acd_0 + g_yyzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyz_0[i] = 3.0 * g_yyzz_0_xxyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyzz_0[i] = 3.0 * g_yyzz_0_xxyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxzzz_0[i] = 3.0 * g_yyzz_0_xxzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyy_0[i] = 2.0 * g_yyzz_0_xyyyy_1[i] * fi_acd_0 + g_yyzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyz_0[i] = 2.0 * g_yyzz_0_xyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyzz_0[i] = 2.0 * g_yyzz_0_xyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyzzz_0[i] = 2.0 * g_yyzz_0_xyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxzzzz_0[i] = 2.0 * g_yyzz_0_xzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyy_0[i] = g_yyzz_0_yyyyy_1[i] * fi_acd_0 + g_yyzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyz_0[i] = g_yyzz_0_yyyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyzz_0[i] = g_yyzz_0_yyyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyzzz_0[i] = g_yyzz_0_yyzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyzzzz_0[i] = g_yyzz_0_yzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xzzzzz_0[i] = g_yyzz_0_zzzzz_1[i] * fi_acd_0 + g_yyzz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyy_0[i] = g_yyzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyz_0[i] = g_yyzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyzz_0[i] = g_yyzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyzzz_0[i] = g_yyzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyzzzz_0[i] = g_yyzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yzzzzz_0[i] = g_yyzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_zzzzzz_0[i] = g_yyzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 364-392 components of targeted buffer : HSI

    auto g_xyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 364);

    auto g_xyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 365);

    auto g_xyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 366);

    auto g_xyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 367);

    auto g_xyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 368);

    auto g_xyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 369);

    auto g_xyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 370);

    auto g_xyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 371);

    auto g_xyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 372);

    auto g_xyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 373);

    auto g_xyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 374);

    auto g_xyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 375);

    auto g_xyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 376);

    auto g_xyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 377);

    auto g_xyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 378);

    auto g_xyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 379);

    auto g_xyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 380);

    auto g_xyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 381);

    auto g_xyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 382);

    auto g_xyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 383);

    auto g_xyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 384);

    auto g_xyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 385);

    auto g_xyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 386);

    auto g_xyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 387);

    auto g_xyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 388);

    auto g_xyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 389);

    auto g_xyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 390);

    auto g_xyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 391);

    #pragma omp simd aligned(g_xyzzz_0_xxxxxx_0, g_xyzzz_0_xxxxxy_0, g_xyzzz_0_xxxxxz_0, g_xyzzz_0_xxxxyy_0, g_xyzzz_0_xxxxyz_0, g_xyzzz_0_xxxxzz_0, g_xyzzz_0_xxxyyy_0, g_xyzzz_0_xxxyyz_0, g_xyzzz_0_xxxyzz_0, g_xyzzz_0_xxxzzz_0, g_xyzzz_0_xxyyyy_0, g_xyzzz_0_xxyyyz_0, g_xyzzz_0_xxyyzz_0, g_xyzzz_0_xxyzzz_0, g_xyzzz_0_xxzzzz_0, g_xyzzz_0_xyyyyy_0, g_xyzzz_0_xyyyyz_0, g_xyzzz_0_xyyyzz_0, g_xyzzz_0_xyyzzz_0, g_xyzzz_0_xyzzzz_0, g_xyzzz_0_xzzzzz_0, g_xyzzz_0_yyyyyy_0, g_xyzzz_0_yyyyyz_0, g_xyzzz_0_yyyyzz_0, g_xyzzz_0_yyyzzz_0, g_xyzzz_0_yyzzzz_0, g_xyzzz_0_yzzzzz_0, g_xyzzz_0_zzzzzz_0, g_xzzz_0_xxxxxx_1, g_xzzz_0_xxxxxz_1, g_xzzz_0_xxxxzz_1, g_xzzz_0_xxxzzz_1, g_xzzz_0_xxzzzz_1, g_xzzz_0_xzzzzz_1, g_yzzz_0_xxxxxy_1, g_yzzz_0_xxxxy_1, g_yzzz_0_xxxxyy_1, g_yzzz_0_xxxxyz_1, g_yzzz_0_xxxyy_1, g_yzzz_0_xxxyyy_1, g_yzzz_0_xxxyyz_1, g_yzzz_0_xxxyz_1, g_yzzz_0_xxxyzz_1, g_yzzz_0_xxyyy_1, g_yzzz_0_xxyyyy_1, g_yzzz_0_xxyyyz_1, g_yzzz_0_xxyyz_1, g_yzzz_0_xxyyzz_1, g_yzzz_0_xxyzz_1, g_yzzz_0_xxyzzz_1, g_yzzz_0_xyyyy_1, g_yzzz_0_xyyyyy_1, g_yzzz_0_xyyyyz_1, g_yzzz_0_xyyyz_1, g_yzzz_0_xyyyzz_1, g_yzzz_0_xyyzz_1, g_yzzz_0_xyyzzz_1, g_yzzz_0_xyzzz_1, g_yzzz_0_xyzzzz_1, g_yzzz_0_yyyyy_1, g_yzzz_0_yyyyyy_1, g_yzzz_0_yyyyyz_1, g_yzzz_0_yyyyz_1, g_yzzz_0_yyyyzz_1, g_yzzz_0_yyyzz_1, g_yzzz_0_yyyzzz_1, g_yzzz_0_yyzzz_1, g_yzzz_0_yyzzzz_1, g_yzzz_0_yzzzz_1, g_yzzz_0_yzzzzz_1, g_yzzz_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzz_0_xxxxxx_0[i] = g_xzzz_0_xxxxxx_1[i] * wa_y[i];

        g_xyzzz_0_xxxxxy_0[i] = 5.0 * g_yzzz_0_xxxxy_1[i] * fi_acd_0 + g_yzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxz_0[i] = g_xzzz_0_xxxxxz_1[i] * wa_y[i];

        g_xyzzz_0_xxxxyy_0[i] = 4.0 * g_yzzz_0_xxxyy_1[i] * fi_acd_0 + g_yzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxyz_0[i] = 4.0 * g_yzzz_0_xxxyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxzz_0[i] = g_xzzz_0_xxxxzz_1[i] * wa_y[i];

        g_xyzzz_0_xxxyyy_0[i] = 3.0 * g_yzzz_0_xxyyy_1[i] * fi_acd_0 + g_yzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxyyz_0[i] = 3.0 * g_yzzz_0_xxyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxyzz_0[i] = 3.0 * g_yzzz_0_xxyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxzzz_0[i] = g_xzzz_0_xxxzzz_1[i] * wa_y[i];

        g_xyzzz_0_xxyyyy_0[i] = 2.0 * g_yzzz_0_xyyyy_1[i] * fi_acd_0 + g_yzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxyyyz_0[i] = 2.0 * g_yzzz_0_xyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxyyzz_0[i] = 2.0 * g_yzzz_0_xyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxyzzz_0[i] = 2.0 * g_yzzz_0_xyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxzzzz_0[i] = g_xzzz_0_xxzzzz_1[i] * wa_y[i];

        g_xyzzz_0_xyyyyy_0[i] = g_yzzz_0_yyyyy_1[i] * fi_acd_0 + g_yzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xyyyyz_0[i] = g_yzzz_0_yyyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xyyyzz_0[i] = g_yzzz_0_yyyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xyyzzz_0[i] = g_yzzz_0_yyzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xyzzzz_0[i] = g_yzzz_0_yzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xzzzzz_0[i] = g_xzzz_0_xzzzzz_1[i] * wa_y[i];

        g_xyzzz_0_yyyyyy_0[i] = g_yzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_yyyyyz_0[i] = g_yzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_yyyyzz_0[i] = g_yzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_yyyzzz_0[i] = g_yzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_yyzzzz_0[i] = g_yzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_yzzzzz_0[i] = g_yzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_zzzzzz_0[i] = g_yzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 392-420 components of targeted buffer : HSI

    auto g_xzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 392);

    auto g_xzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 393);

    auto g_xzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 394);

    auto g_xzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 395);

    auto g_xzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 396);

    auto g_xzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 397);

    auto g_xzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 398);

    auto g_xzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 399);

    auto g_xzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 400);

    auto g_xzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 401);

    auto g_xzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 402);

    auto g_xzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 403);

    auto g_xzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 404);

    auto g_xzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 405);

    auto g_xzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 406);

    auto g_xzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 407);

    auto g_xzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 408);

    auto g_xzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 409);

    auto g_xzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 410);

    auto g_xzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 411);

    auto g_xzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 412);

    auto g_xzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 413);

    auto g_xzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 414);

    auto g_xzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 415);

    auto g_xzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 416);

    auto g_xzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 417);

    auto g_xzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 418);

    auto g_xzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 419);

    #pragma omp simd aligned(g_xzzzz_0_xxxxxx_0, g_xzzzz_0_xxxxxy_0, g_xzzzz_0_xxxxxz_0, g_xzzzz_0_xxxxyy_0, g_xzzzz_0_xxxxyz_0, g_xzzzz_0_xxxxzz_0, g_xzzzz_0_xxxyyy_0, g_xzzzz_0_xxxyyz_0, g_xzzzz_0_xxxyzz_0, g_xzzzz_0_xxxzzz_0, g_xzzzz_0_xxyyyy_0, g_xzzzz_0_xxyyyz_0, g_xzzzz_0_xxyyzz_0, g_xzzzz_0_xxyzzz_0, g_xzzzz_0_xxzzzz_0, g_xzzzz_0_xyyyyy_0, g_xzzzz_0_xyyyyz_0, g_xzzzz_0_xyyyzz_0, g_xzzzz_0_xyyzzz_0, g_xzzzz_0_xyzzzz_0, g_xzzzz_0_xzzzzz_0, g_xzzzz_0_yyyyyy_0, g_xzzzz_0_yyyyyz_0, g_xzzzz_0_yyyyzz_0, g_xzzzz_0_yyyzzz_0, g_xzzzz_0_yyzzzz_0, g_xzzzz_0_yzzzzz_0, g_xzzzz_0_zzzzzz_0, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxy_1, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxy_1, g_zzzz_0_xxxxyy_1, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxyy_1, g_zzzz_0_xxxyyy_1, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxzz_1, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxyyy_1, g_zzzz_0_xxyyyy_1, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyzz_1, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxzzz_1, g_zzzz_0_xxzzzz_1, g_zzzz_0_xyyyy_1, g_zzzz_0_xyyyyy_1, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyzz_1, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyzzz_1, g_zzzz_0_xyzzzz_1, g_zzzz_0_xzzzz_1, g_zzzz_0_xzzzzz_1, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyzz_1, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyzzz_1, g_zzzz_0_yyzzzz_1, g_zzzz_0_yzzzz_1, g_zzzz_0_yzzzzz_1, g_zzzz_0_zzzzz_1, g_zzzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_xxxxxx_0[i] = 6.0 * g_zzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxx_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxy_0[i] = 5.0 * g_zzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxz_0[i] = 5.0 * g_zzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyy_0[i] = 4.0 * g_zzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyz_0[i] = 4.0 * g_zzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxzz_0[i] = 4.0 * g_zzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyy_0[i] = 3.0 * g_zzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyz_0[i] = 3.0 * g_zzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyzz_0[i] = 3.0 * g_zzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxzzz_0[i] = 3.0 * g_zzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyy_0[i] = 2.0 * g_zzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyz_0[i] = 2.0 * g_zzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyzz_0[i] = 2.0 * g_zzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyzzz_0[i] = 2.0 * g_zzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxzzzz_0[i] = 2.0 * g_zzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyy_0[i] = g_zzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyz_0[i] = g_zzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyzz_0[i] = g_zzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyzzz_0[i] = g_zzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyzzzz_0[i] = g_zzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xzzzzz_0[i] = g_zzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyy_0[i] = g_zzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyz_0[i] = g_zzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyzz_0[i] = g_zzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyzzz_0[i] = g_zzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyzzzz_0[i] = g_zzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yzzzzz_0[i] = g_zzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_zzzzzz_0[i] = g_zzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 420-448 components of targeted buffer : HSI

    auto g_yyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 420);

    auto g_yyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 421);

    auto g_yyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 422);

    auto g_yyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 423);

    auto g_yyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 424);

    auto g_yyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 425);

    auto g_yyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 426);

    auto g_yyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 427);

    auto g_yyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 428);

    auto g_yyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 429);

    auto g_yyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 430);

    auto g_yyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 431);

    auto g_yyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 432);

    auto g_yyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 433);

    auto g_yyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 434);

    auto g_yyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 435);

    auto g_yyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 436);

    auto g_yyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 437);

    auto g_yyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 438);

    auto g_yyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 439);

    auto g_yyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 440);

    auto g_yyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 441);

    auto g_yyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 442);

    auto g_yyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 443);

    auto g_yyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 444);

    auto g_yyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 445);

    auto g_yyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 446);

    auto g_yyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 447);

    #pragma omp simd aligned(g_yyy_0_xxxxxx_0, g_yyy_0_xxxxxx_1, g_yyy_0_xxxxxy_0, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxxz_0, g_yyy_0_xxxxxz_1, g_yyy_0_xxxxyy_0, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyz_0, g_yyy_0_xxxxyz_1, g_yyy_0_xxxxzz_0, g_yyy_0_xxxxzz_1, g_yyy_0_xxxyyy_0, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyz_0, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyzz_0, g_yyy_0_xxxyzz_1, g_yyy_0_xxxzzz_0, g_yyy_0_xxxzzz_1, g_yyy_0_xxyyyy_0, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyz_0, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyzz_0, g_yyy_0_xxyyzz_1, g_yyy_0_xxyzzz_0, g_yyy_0_xxyzzz_1, g_yyy_0_xxzzzz_0, g_yyy_0_xxzzzz_1, g_yyy_0_xyyyyy_0, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyz_0, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyzz_0, g_yyy_0_xyyyzz_1, g_yyy_0_xyyzzz_0, g_yyy_0_xyyzzz_1, g_yyy_0_xyzzzz_0, g_yyy_0_xyzzzz_1, g_yyy_0_xzzzzz_0, g_yyy_0_xzzzzz_1, g_yyy_0_yyyyyy_0, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyz_0, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyzz_0, g_yyy_0_yyyyzz_1, g_yyy_0_yyyzzz_0, g_yyy_0_yyyzzz_1, g_yyy_0_yyzzzz_0, g_yyy_0_yyzzzz_1, g_yyy_0_yzzzzz_0, g_yyy_0_yzzzzz_1, g_yyy_0_zzzzzz_0, g_yyy_0_zzzzzz_1, g_yyyy_0_xxxxx_1, g_yyyy_0_xxxxxx_1, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxxz_1, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxxz_1, g_yyyy_0_xxxxzz_1, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyz_1, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxxzz_1, g_yyyy_0_xxxzzz_1, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyzz_1, g_yyyy_0_xxyzzz_1, g_yyyy_0_xxzzz_1, g_yyyy_0_xxzzzz_1, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyzz_1, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyzzz_1, g_yyyy_0_xyzzzz_1, g_yyyy_0_xzzzz_1, g_yyyy_0_xzzzzz_1, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyzz_1, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyzzz_1, g_yyyy_0_yyzzzz_1, g_yyyy_0_yzzzz_1, g_yyyy_0_yzzzzz_1, g_yyyy_0_zzzzz_1, g_yyyy_0_zzzzzz_1, g_yyyyy_0_xxxxxx_0, g_yyyyy_0_xxxxxy_0, g_yyyyy_0_xxxxxz_0, g_yyyyy_0_xxxxyy_0, g_yyyyy_0_xxxxyz_0, g_yyyyy_0_xxxxzz_0, g_yyyyy_0_xxxyyy_0, g_yyyyy_0_xxxyyz_0, g_yyyyy_0_xxxyzz_0, g_yyyyy_0_xxxzzz_0, g_yyyyy_0_xxyyyy_0, g_yyyyy_0_xxyyyz_0, g_yyyyy_0_xxyyzz_0, g_yyyyy_0_xxyzzz_0, g_yyyyy_0_xxzzzz_0, g_yyyyy_0_xyyyyy_0, g_yyyyy_0_xyyyyz_0, g_yyyyy_0_xyyyzz_0, g_yyyyy_0_xyyzzz_0, g_yyyyy_0_xyzzzz_0, g_yyyyy_0_xzzzzz_0, g_yyyyy_0_yyyyyy_0, g_yyyyy_0_yyyyyz_0, g_yyyyy_0_yyyyzz_0, g_yyyyy_0_yyyzzz_0, g_yyyyy_0_yyzzzz_0, g_yyyyy_0_yzzzzz_0, g_yyyyy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_xxxxxx_0[i] = 4.0 * g_yyy_0_xxxxxx_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxx_1[i] * fz_be_0 + g_yyyy_0_xxxxxx_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxy_0[i] = 4.0 * g_yyy_0_xxxxxy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxy_1[i] * fz_be_0 + g_yyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxz_0[i] = 4.0 * g_yyy_0_xxxxxz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxz_1[i] * fz_be_0 + g_yyyy_0_xxxxxz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyy_0[i] = 4.0 * g_yyy_0_xxxxyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyz_0[i] = 4.0 * g_yyy_0_xxxxyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyz_1[i] * fz_be_0 + g_yyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxzz_0[i] = 4.0 * g_yyy_0_xxxxzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxzz_1[i] * fz_be_0 + g_yyyy_0_xxxxzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyy_0[i] = 4.0 * g_yyy_0_xxxyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyz_0[i] = 4.0 * g_yyy_0_xxxyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyzz_0[i] = 4.0 * g_yyy_0_xxxyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyzz_1[i] * fz_be_0 + g_yyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxzzz_0[i] = 4.0 * g_yyy_0_xxxzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxzzz_1[i] * fz_be_0 + g_yyyy_0_xxxzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyy_0[i] = 4.0 * g_yyy_0_xxyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyz_0[i] = 4.0 * g_yyy_0_xxyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyzz_0[i] = 4.0 * g_yyy_0_xxyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyzzz_0[i] = 4.0 * g_yyy_0_xxyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyzzz_1[i] * fz_be_0 + g_yyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxzzzz_0[i] = 4.0 * g_yyy_0_xxzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxzzzz_1[i] * fz_be_0 + g_yyyy_0_xxzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyy_0[i] = 4.0 * g_yyy_0_xyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyz_0[i] = 4.0 * g_yyy_0_xyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyzz_0[i] = 4.0 * g_yyy_0_xyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyzzz_0[i] = 4.0 * g_yyy_0_xyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyzzzz_0[i] = 4.0 * g_yyy_0_xyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyzzzz_1[i] * fz_be_0 + g_yyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xzzzzz_0[i] = 4.0 * g_yyy_0_xzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xzzzzz_1[i] * fz_be_0 + g_yyyy_0_xzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyy_0[i] = 4.0 * g_yyy_0_yyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyz_0[i] = 4.0 * g_yyy_0_yyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyzz_0[i] = 4.0 * g_yyy_0_yyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyzzz_0[i] = 4.0 * g_yyy_0_yyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyzzzz_0[i] = 4.0 * g_yyy_0_yyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yzzzzz_0[i] = 4.0 * g_yyy_0_yzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yzzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_zzzzzz_0[i] = 4.0 * g_yyy_0_zzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_zzzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 448-476 components of targeted buffer : HSI

    auto g_yyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 448);

    auto g_yyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 449);

    auto g_yyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 450);

    auto g_yyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 451);

    auto g_yyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 452);

    auto g_yyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 453);

    auto g_yyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 454);

    auto g_yyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 455);

    auto g_yyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 456);

    auto g_yyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 457);

    auto g_yyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 458);

    auto g_yyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 459);

    auto g_yyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 460);

    auto g_yyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 461);

    auto g_yyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 462);

    auto g_yyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 463);

    auto g_yyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 464);

    auto g_yyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 465);

    auto g_yyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 466);

    auto g_yyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 467);

    auto g_yyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 468);

    auto g_yyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 469);

    auto g_yyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 470);

    auto g_yyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 471);

    auto g_yyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 472);

    auto g_yyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 473);

    auto g_yyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 474);

    auto g_yyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 475);

    #pragma omp simd aligned(g_yyyy_0_xxxxx_1, g_yyyy_0_xxxxxx_1, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxxz_1, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxxz_1, g_yyyy_0_xxxxzz_1, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyz_1, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxxzz_1, g_yyyy_0_xxxzzz_1, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyzz_1, g_yyyy_0_xxyzzz_1, g_yyyy_0_xxzzz_1, g_yyyy_0_xxzzzz_1, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyzz_1, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyzzz_1, g_yyyy_0_xyzzzz_1, g_yyyy_0_xzzzz_1, g_yyyy_0_xzzzzz_1, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyzz_1, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyzzz_1, g_yyyy_0_yyzzzz_1, g_yyyy_0_yzzzz_1, g_yyyy_0_yzzzzz_1, g_yyyy_0_zzzzz_1, g_yyyy_0_zzzzzz_1, g_yyyyz_0_xxxxxx_0, g_yyyyz_0_xxxxxy_0, g_yyyyz_0_xxxxxz_0, g_yyyyz_0_xxxxyy_0, g_yyyyz_0_xxxxyz_0, g_yyyyz_0_xxxxzz_0, g_yyyyz_0_xxxyyy_0, g_yyyyz_0_xxxyyz_0, g_yyyyz_0_xxxyzz_0, g_yyyyz_0_xxxzzz_0, g_yyyyz_0_xxyyyy_0, g_yyyyz_0_xxyyyz_0, g_yyyyz_0_xxyyzz_0, g_yyyyz_0_xxyzzz_0, g_yyyyz_0_xxzzzz_0, g_yyyyz_0_xyyyyy_0, g_yyyyz_0_xyyyyz_0, g_yyyyz_0_xyyyzz_0, g_yyyyz_0_xyyzzz_0, g_yyyyz_0_xyzzzz_0, g_yyyyz_0_xzzzzz_0, g_yyyyz_0_yyyyyy_0, g_yyyyz_0_yyyyyz_0, g_yyyyz_0_yyyyzz_0, g_yyyyz_0_yyyzzz_0, g_yyyyz_0_yyzzzz_0, g_yyyyz_0_yzzzzz_0, g_yyyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_xxxxxx_0[i] = g_yyyy_0_xxxxxx_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxy_0[i] = g_yyyy_0_xxxxxy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxz_0[i] = g_yyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyy_0[i] = g_yyyy_0_xxxxyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyz_0[i] = g_yyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxzz_0[i] = 2.0 * g_yyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyy_0[i] = g_yyyy_0_xxxyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyz_0[i] = g_yyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyzz_0[i] = 2.0 * g_yyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxzzz_0[i] = 3.0 * g_yyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyy_0[i] = g_yyyy_0_xxyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyz_0[i] = g_yyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyzz_0[i] = 2.0 * g_yyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyzzz_0[i] = 3.0 * g_yyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxzzzz_0[i] = 4.0 * g_yyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyy_0[i] = g_yyyy_0_xyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyz_0[i] = g_yyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyzz_0[i] = 2.0 * g_yyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyzzz_0[i] = 3.0 * g_yyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyzzzz_0[i] = 4.0 * g_yyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xzzzzz_0[i] = 5.0 * g_yyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyy_0[i] = g_yyyy_0_yyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyz_0[i] = g_yyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyzz_0[i] = 2.0 * g_yyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyzzz_0[i] = 3.0 * g_yyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyy_0_yyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyzzzz_0[i] = 4.0 * g_yyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yzzzzz_0[i] = 5.0 * g_yyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_zzzzzz_0[i] = 6.0 * g_yyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 476-504 components of targeted buffer : HSI

    auto g_yyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 476);

    auto g_yyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 477);

    auto g_yyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 478);

    auto g_yyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 479);

    auto g_yyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 480);

    auto g_yyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 481);

    auto g_yyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 482);

    auto g_yyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 483);

    auto g_yyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 484);

    auto g_yyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 485);

    auto g_yyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 486);

    auto g_yyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 487);

    auto g_yyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 488);

    auto g_yyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 489);

    auto g_yyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 490);

    auto g_yyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 491);

    auto g_yyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 492);

    auto g_yyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 493);

    auto g_yyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 494);

    auto g_yyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 495);

    auto g_yyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 496);

    auto g_yyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 497);

    auto g_yyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 498);

    auto g_yyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 499);

    auto g_yyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 500);

    auto g_yyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 501);

    auto g_yyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 502);

    auto g_yyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 503);

    #pragma omp simd aligned(g_yyy_0_xxxxxy_0, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxyy_0, g_yyy_0_xxxxyy_1, g_yyy_0_xxxyyy_0, g_yyy_0_xxxyyy_1, g_yyy_0_xxyyyy_0, g_yyy_0_xxyyyy_1, g_yyy_0_xyyyyy_0, g_yyy_0_xyyyyy_1, g_yyy_0_yyyyyy_0, g_yyy_0_yyyyyy_1, g_yyyz_0_xxxxxy_1, g_yyyz_0_xxxxyy_1, g_yyyz_0_xxxyyy_1, g_yyyz_0_xxyyyy_1, g_yyyz_0_xyyyyy_1, g_yyyz_0_yyyyyy_1, g_yyyzz_0_xxxxxx_0, g_yyyzz_0_xxxxxy_0, g_yyyzz_0_xxxxxz_0, g_yyyzz_0_xxxxyy_0, g_yyyzz_0_xxxxyz_0, g_yyyzz_0_xxxxzz_0, g_yyyzz_0_xxxyyy_0, g_yyyzz_0_xxxyyz_0, g_yyyzz_0_xxxyzz_0, g_yyyzz_0_xxxzzz_0, g_yyyzz_0_xxyyyy_0, g_yyyzz_0_xxyyyz_0, g_yyyzz_0_xxyyzz_0, g_yyyzz_0_xxyzzz_0, g_yyyzz_0_xxzzzz_0, g_yyyzz_0_xyyyyy_0, g_yyyzz_0_xyyyyz_0, g_yyyzz_0_xyyyzz_0, g_yyyzz_0_xyyzzz_0, g_yyyzz_0_xyzzzz_0, g_yyyzz_0_xzzzzz_0, g_yyyzz_0_yyyyyy_0, g_yyyzz_0_yyyyyz_0, g_yyyzz_0_yyyyzz_0, g_yyyzz_0_yyyzzz_0, g_yyyzz_0_yyzzzz_0, g_yyyzz_0_yzzzzz_0, g_yyyzz_0_zzzzzz_0, g_yyzz_0_xxxxxx_1, g_yyzz_0_xxxxxz_1, g_yyzz_0_xxxxyz_1, g_yyzz_0_xxxxz_1, g_yyzz_0_xxxxzz_1, g_yyzz_0_xxxyyz_1, g_yyzz_0_xxxyz_1, g_yyzz_0_xxxyzz_1, g_yyzz_0_xxxzz_1, g_yyzz_0_xxxzzz_1, g_yyzz_0_xxyyyz_1, g_yyzz_0_xxyyz_1, g_yyzz_0_xxyyzz_1, g_yyzz_0_xxyzz_1, g_yyzz_0_xxyzzz_1, g_yyzz_0_xxzzz_1, g_yyzz_0_xxzzzz_1, g_yyzz_0_xyyyyz_1, g_yyzz_0_xyyyz_1, g_yyzz_0_xyyyzz_1, g_yyzz_0_xyyzz_1, g_yyzz_0_xyyzzz_1, g_yyzz_0_xyzzz_1, g_yyzz_0_xyzzzz_1, g_yyzz_0_xzzzz_1, g_yyzz_0_xzzzzz_1, g_yyzz_0_yyyyyz_1, g_yyzz_0_yyyyz_1, g_yyzz_0_yyyyzz_1, g_yyzz_0_yyyzz_1, g_yyzz_0_yyyzzz_1, g_yyzz_0_yyzzz_1, g_yyzz_0_yyzzzz_1, g_yyzz_0_yzzzz_1, g_yyzz_0_yzzzzz_1, g_yyzz_0_zzzzz_1, g_yyzz_0_zzzzzz_1, g_yzz_0_xxxxxx_0, g_yzz_0_xxxxxx_1, g_yzz_0_xxxxxz_0, g_yzz_0_xxxxxz_1, g_yzz_0_xxxxyz_0, g_yzz_0_xxxxyz_1, g_yzz_0_xxxxzz_0, g_yzz_0_xxxxzz_1, g_yzz_0_xxxyyz_0, g_yzz_0_xxxyyz_1, g_yzz_0_xxxyzz_0, g_yzz_0_xxxyzz_1, g_yzz_0_xxxzzz_0, g_yzz_0_xxxzzz_1, g_yzz_0_xxyyyz_0, g_yzz_0_xxyyyz_1, g_yzz_0_xxyyzz_0, g_yzz_0_xxyyzz_1, g_yzz_0_xxyzzz_0, g_yzz_0_xxyzzz_1, g_yzz_0_xxzzzz_0, g_yzz_0_xxzzzz_1, g_yzz_0_xyyyyz_0, g_yzz_0_xyyyyz_1, g_yzz_0_xyyyzz_0, g_yzz_0_xyyyzz_1, g_yzz_0_xyyzzz_0, g_yzz_0_xyyzzz_1, g_yzz_0_xyzzzz_0, g_yzz_0_xyzzzz_1, g_yzz_0_xzzzzz_0, g_yzz_0_xzzzzz_1, g_yzz_0_yyyyyz_0, g_yzz_0_yyyyyz_1, g_yzz_0_yyyyzz_0, g_yzz_0_yyyyzz_1, g_yzz_0_yyyzzz_0, g_yzz_0_yyyzzz_1, g_yzz_0_yyzzzz_0, g_yzz_0_yyzzzz_1, g_yzz_0_yzzzzz_0, g_yzz_0_yzzzzz_1, g_yzz_0_zzzzzz_0, g_yzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzz_0_xxxxxx_0[i] = 2.0 * g_yzz_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxx_1[i] * fz_be_0 + g_yyzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxy_0[i] = g_yyy_0_xxxxxy_0[i] * fbe_0 - g_yyy_0_xxxxxy_1[i] * fz_be_0 + g_yyyz_0_xxxxxy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxxz_0[i] = 2.0 * g_yzz_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxz_1[i] * fz_be_0 + g_yyzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxyy_0[i] = g_yyy_0_xxxxyy_0[i] * fbe_0 - g_yyy_0_xxxxyy_1[i] * fz_be_0 + g_yyyz_0_xxxxyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxyz_0[i] = 2.0 * g_yzz_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxyz_1[i] * fz_be_0 + g_yyzz_0_xxxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxzz_0[i] = 2.0 * g_yzz_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxzz_1[i] * fz_be_0 + g_yyzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyyy_0[i] = g_yyy_0_xxxyyy_0[i] * fbe_0 - g_yyy_0_xxxyyy_1[i] * fz_be_0 + g_yyyz_0_xxxyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxyyz_0[i] = 2.0 * g_yzz_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxxyz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyzz_0[i] = 2.0 * g_yzz_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyzz_1[i] * fz_be_0 + g_yyzz_0_xxxzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxzzz_0[i] = 2.0 * g_yzz_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxzzz_1[i] * fz_be_0 + g_yyzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyyy_0[i] = g_yyy_0_xxyyyy_0[i] * fbe_0 - g_yyy_0_xxyyyy_1[i] * fz_be_0 + g_yyyz_0_xxyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxyyyz_0[i] = 2.0 * g_yzz_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xxyyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyzz_0[i] = 2.0 * g_yzz_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxyzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyzzz_0[i] = 2.0 * g_yzz_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyzzz_1[i] * fz_be_0 + g_yyzz_0_xxzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxzzzz_0[i] = 2.0 * g_yzz_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxzzzz_1[i] * fz_be_0 + g_yyzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyyy_0[i] = g_yyy_0_xyyyyy_0[i] * fbe_0 - g_yyy_0_xyyyyy_1[i] * fz_be_0 + g_yyyz_0_xyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xyyyyz_0[i] = 2.0 * g_yzz_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_xyyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyzz_0[i] = 2.0 * g_yzz_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xyyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyzzz_0[i] = 2.0 * g_yzz_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xyzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyzzzz_0[i] = 2.0 * g_yzz_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyzzzz_1[i] * fz_be_0 + g_yyzz_0_xzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xzzzzz_0[i] = 2.0 * g_yzz_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xzzzzz_1[i] * fz_be_0 + g_yyzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyyy_0[i] = g_yyy_0_yyyyyy_0[i] * fbe_0 - g_yyy_0_yyyyyy_1[i] * fz_be_0 + g_yyyz_0_yyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_yyyyyz_0[i] = 2.0 * g_yzz_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzz_0_yyyyz_1[i] * fi_acd_0 + g_yyzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyzz_0[i] = 2.0 * g_yzz_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_yyyzz_1[i] * fi_acd_0 + g_yyzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyzzz_0[i] = 2.0 * g_yzz_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_yyzzz_1[i] * fi_acd_0 + g_yyzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyzzzz_0[i] = 2.0 * g_yzz_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_yzzzz_1[i] * fi_acd_0 + g_yyzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yzzzzz_0[i] = 2.0 * g_yzz_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yzzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzzz_1[i] * fi_acd_0 + g_yyzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_zzzzzz_0[i] = 2.0 * g_yzz_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_zzzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 504-532 components of targeted buffer : HSI

    auto g_yyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 504);

    auto g_yyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 505);

    auto g_yyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 506);

    auto g_yyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 507);

    auto g_yyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 508);

    auto g_yyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 509);

    auto g_yyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 510);

    auto g_yyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 511);

    auto g_yyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 512);

    auto g_yyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 513);

    auto g_yyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 514);

    auto g_yyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 515);

    auto g_yyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 516);

    auto g_yyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 517);

    auto g_yyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 518);

    auto g_yyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 519);

    auto g_yyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 520);

    auto g_yyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 521);

    auto g_yyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 522);

    auto g_yyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 523);

    auto g_yyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 524);

    auto g_yyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 525);

    auto g_yyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 526);

    auto g_yyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 527);

    auto g_yyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 528);

    auto g_yyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 529);

    auto g_yyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 530);

    auto g_yyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 531);

    #pragma omp simd aligned(g_yyz_0_xxxxxy_0, g_yyz_0_xxxxxy_1, g_yyz_0_xxxxyy_0, g_yyz_0_xxxxyy_1, g_yyz_0_xxxyyy_0, g_yyz_0_xxxyyy_1, g_yyz_0_xxyyyy_0, g_yyz_0_xxyyyy_1, g_yyz_0_xyyyyy_0, g_yyz_0_xyyyyy_1, g_yyz_0_yyyyyy_0, g_yyz_0_yyyyyy_1, g_yyzz_0_xxxxxy_1, g_yyzz_0_xxxxyy_1, g_yyzz_0_xxxyyy_1, g_yyzz_0_xxyyyy_1, g_yyzz_0_xyyyyy_1, g_yyzz_0_yyyyyy_1, g_yyzzz_0_xxxxxx_0, g_yyzzz_0_xxxxxy_0, g_yyzzz_0_xxxxxz_0, g_yyzzz_0_xxxxyy_0, g_yyzzz_0_xxxxyz_0, g_yyzzz_0_xxxxzz_0, g_yyzzz_0_xxxyyy_0, g_yyzzz_0_xxxyyz_0, g_yyzzz_0_xxxyzz_0, g_yyzzz_0_xxxzzz_0, g_yyzzz_0_xxyyyy_0, g_yyzzz_0_xxyyyz_0, g_yyzzz_0_xxyyzz_0, g_yyzzz_0_xxyzzz_0, g_yyzzz_0_xxzzzz_0, g_yyzzz_0_xyyyyy_0, g_yyzzz_0_xyyyyz_0, g_yyzzz_0_xyyyzz_0, g_yyzzz_0_xyyzzz_0, g_yyzzz_0_xyzzzz_0, g_yyzzz_0_xzzzzz_0, g_yyzzz_0_yyyyyy_0, g_yyzzz_0_yyyyyz_0, g_yyzzz_0_yyyyzz_0, g_yyzzz_0_yyyzzz_0, g_yyzzz_0_yyzzzz_0, g_yyzzz_0_yzzzzz_0, g_yyzzz_0_zzzzzz_0, g_yzzz_0_xxxxxx_1, g_yzzz_0_xxxxxz_1, g_yzzz_0_xxxxyz_1, g_yzzz_0_xxxxz_1, g_yzzz_0_xxxxzz_1, g_yzzz_0_xxxyyz_1, g_yzzz_0_xxxyz_1, g_yzzz_0_xxxyzz_1, g_yzzz_0_xxxzz_1, g_yzzz_0_xxxzzz_1, g_yzzz_0_xxyyyz_1, g_yzzz_0_xxyyz_1, g_yzzz_0_xxyyzz_1, g_yzzz_0_xxyzz_1, g_yzzz_0_xxyzzz_1, g_yzzz_0_xxzzz_1, g_yzzz_0_xxzzzz_1, g_yzzz_0_xyyyyz_1, g_yzzz_0_xyyyz_1, g_yzzz_0_xyyyzz_1, g_yzzz_0_xyyzz_1, g_yzzz_0_xyyzzz_1, g_yzzz_0_xyzzz_1, g_yzzz_0_xyzzzz_1, g_yzzz_0_xzzzz_1, g_yzzz_0_xzzzzz_1, g_yzzz_0_yyyyyz_1, g_yzzz_0_yyyyz_1, g_yzzz_0_yyyyzz_1, g_yzzz_0_yyyzz_1, g_yzzz_0_yyyzzz_1, g_yzzz_0_yyzzz_1, g_yzzz_0_yyzzzz_1, g_yzzz_0_yzzzz_1, g_yzzz_0_yzzzzz_1, g_yzzz_0_zzzzz_1, g_yzzz_0_zzzzzz_1, g_zzz_0_xxxxxx_0, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxz_0, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxyz_0, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxzz_0, g_zzz_0_xxxxzz_1, g_zzz_0_xxxyyz_0, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyzz_0, g_zzz_0_xxxyzz_1, g_zzz_0_xxxzzz_0, g_zzz_0_xxxzzz_1, g_zzz_0_xxyyyz_0, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyzz_0, g_zzz_0_xxyyzz_1, g_zzz_0_xxyzzz_0, g_zzz_0_xxyzzz_1, g_zzz_0_xxzzzz_0, g_zzz_0_xxzzzz_1, g_zzz_0_xyyyyz_0, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyzz_0, g_zzz_0_xyyyzz_1, g_zzz_0_xyyzzz_0, g_zzz_0_xyyzzz_1, g_zzz_0_xyzzzz_0, g_zzz_0_xyzzzz_1, g_zzz_0_xzzzzz_0, g_zzz_0_xzzzzz_1, g_zzz_0_yyyyyz_0, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyzz_0, g_zzz_0_yyyyzz_1, g_zzz_0_yyyzzz_0, g_zzz_0_yyyzzz_1, g_zzz_0_yyzzzz_0, g_zzz_0_yyzzzz_1, g_zzz_0_yzzzzz_0, g_zzz_0_yzzzzz_1, g_zzz_0_zzzzzz_0, g_zzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzz_0_xxxxxx_0[i] = g_zzz_0_xxxxxx_0[i] * fbe_0 - g_zzz_0_xxxxxx_1[i] * fz_be_0 + g_yzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxy_0[i] = 2.0 * g_yyz_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxxy_1[i] * fz_be_0 + g_yyzz_0_xxxxxy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxxz_0[i] = g_zzz_0_xxxxxz_0[i] * fbe_0 - g_zzz_0_xxxxxz_1[i] * fz_be_0 + g_yzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxyy_0[i] = 2.0 * g_yyz_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxyy_1[i] * fz_be_0 + g_yyzz_0_xxxxyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxyz_0[i] = g_zzz_0_xxxxyz_0[i] * fbe_0 - g_zzz_0_xxxxyz_1[i] * fz_be_0 + g_yzzz_0_xxxxz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxzz_0[i] = g_zzz_0_xxxxzz_0[i] * fbe_0 - g_zzz_0_xxxxzz_1[i] * fz_be_0 + g_yzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyyy_0[i] = 2.0 * g_yyz_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxyyy_1[i] * fz_be_0 + g_yyzz_0_xxxyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxyyz_0[i] = g_zzz_0_xxxyyz_0[i] * fbe_0 - g_zzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxxyz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyzz_0[i] = g_zzz_0_xxxyzz_0[i] * fbe_0 - g_zzz_0_xxxyzz_1[i] * fz_be_0 + g_yzzz_0_xxxzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxzzz_0[i] = g_zzz_0_xxxzzz_0[i] * fbe_0 - g_zzz_0_xxxzzz_1[i] * fz_be_0 + g_yzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyyy_0[i] = 2.0 * g_yyz_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxyyyy_1[i] * fz_be_0 + g_yyzz_0_xxyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxyyyz_0[i] = g_zzz_0_xxyyyz_0[i] * fbe_0 - g_zzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xxyyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyzz_0[i] = g_zzz_0_xxyyzz_0[i] * fbe_0 - g_zzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxyzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyzzz_0[i] = g_zzz_0_xxyzzz_0[i] * fbe_0 - g_zzz_0_xxyzzz_1[i] * fz_be_0 + g_yzzz_0_xxzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxzzzz_0[i] = g_zzz_0_xxzzzz_0[i] * fbe_0 - g_zzz_0_xxzzzz_1[i] * fz_be_0 + g_yzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyyy_0[i] = 2.0 * g_yyz_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xyyyyy_1[i] * fz_be_0 + g_yyzz_0_xyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xyyyyz_0[i] = g_zzz_0_xyyyyz_0[i] * fbe_0 - g_zzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_xyyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyzz_0[i] = g_zzz_0_xyyyzz_0[i] * fbe_0 - g_zzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xyyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyzzz_0[i] = g_zzz_0_xyyzzz_0[i] * fbe_0 - g_zzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xyzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyzzzz_0[i] = g_zzz_0_xyzzzz_0[i] * fbe_0 - g_zzz_0_xyzzzz_1[i] * fz_be_0 + g_yzzz_0_xzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xzzzzz_0[i] = g_zzz_0_xzzzzz_0[i] * fbe_0 - g_zzz_0_xzzzzz_1[i] * fz_be_0 + g_yzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyyy_0[i] = 2.0 * g_yyz_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_yyyyyy_1[i] * fz_be_0 + g_yyzz_0_yyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_yyyyyz_0[i] = g_zzz_0_yyyyyz_0[i] * fbe_0 - g_zzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzz_0_yyyyz_1[i] * fi_acd_0 + g_yzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyzz_0[i] = g_zzz_0_yyyyzz_0[i] * fbe_0 - g_zzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_yyyzz_1[i] * fi_acd_0 + g_yzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyzzz_0[i] = g_zzz_0_yyyzzz_0[i] * fbe_0 - g_zzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_yyzzz_1[i] * fi_acd_0 + g_yzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyzzzz_0[i] = g_zzz_0_yyzzzz_0[i] * fbe_0 - g_zzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_yzzzz_1[i] * fi_acd_0 + g_yzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yzzzzz_0[i] = g_zzz_0_yzzzzz_0[i] * fbe_0 - g_zzz_0_yzzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzzz_1[i] * fi_acd_0 + g_yzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_zzzzzz_0[i] = g_zzz_0_zzzzzz_0[i] * fbe_0 - g_zzz_0_zzzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 532-560 components of targeted buffer : HSI

    auto g_yzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 532);

    auto g_yzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 533);

    auto g_yzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 534);

    auto g_yzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 535);

    auto g_yzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 536);

    auto g_yzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 537);

    auto g_yzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 538);

    auto g_yzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 539);

    auto g_yzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 540);

    auto g_yzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 541);

    auto g_yzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 542);

    auto g_yzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 543);

    auto g_yzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 544);

    auto g_yzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 545);

    auto g_yzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 546);

    auto g_yzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 547);

    auto g_yzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 548);

    auto g_yzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 549);

    auto g_yzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 550);

    auto g_yzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 551);

    auto g_yzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 552);

    auto g_yzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 553);

    auto g_yzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 554);

    auto g_yzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 555);

    auto g_yzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 556);

    auto g_yzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 557);

    auto g_yzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 558);

    auto g_yzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 559);

    #pragma omp simd aligned(g_yzzzz_0_xxxxxx_0, g_yzzzz_0_xxxxxy_0, g_yzzzz_0_xxxxxz_0, g_yzzzz_0_xxxxyy_0, g_yzzzz_0_xxxxyz_0, g_yzzzz_0_xxxxzz_0, g_yzzzz_0_xxxyyy_0, g_yzzzz_0_xxxyyz_0, g_yzzzz_0_xxxyzz_0, g_yzzzz_0_xxxzzz_0, g_yzzzz_0_xxyyyy_0, g_yzzzz_0_xxyyyz_0, g_yzzzz_0_xxyyzz_0, g_yzzzz_0_xxyzzz_0, g_yzzzz_0_xxzzzz_0, g_yzzzz_0_xyyyyy_0, g_yzzzz_0_xyyyyz_0, g_yzzzz_0_xyyyzz_0, g_yzzzz_0_xyyzzz_0, g_yzzzz_0_xyzzzz_0, g_yzzzz_0_xzzzzz_0, g_yzzzz_0_yyyyyy_0, g_yzzzz_0_yyyyyz_0, g_yzzzz_0_yyyyzz_0, g_yzzzz_0_yyyzzz_0, g_yzzzz_0_yyzzzz_0, g_yzzzz_0_yzzzzz_0, g_yzzzz_0_zzzzzz_0, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxy_1, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxy_1, g_zzzz_0_xxxxyy_1, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxyy_1, g_zzzz_0_xxxyyy_1, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxzz_1, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxyyy_1, g_zzzz_0_xxyyyy_1, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyzz_1, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxzzz_1, g_zzzz_0_xxzzzz_1, g_zzzz_0_xyyyy_1, g_zzzz_0_xyyyyy_1, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyzz_1, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyzzz_1, g_zzzz_0_xyzzzz_1, g_zzzz_0_xzzzz_1, g_zzzz_0_xzzzzz_1, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyzz_1, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyzzz_1, g_zzzz_0_yyzzzz_1, g_zzzz_0_yzzzz_1, g_zzzz_0_yzzzzz_1, g_zzzz_0_zzzzz_1, g_zzzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_xxxxxx_0[i] = g_zzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxy_0[i] = g_zzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxz_0[i] = g_zzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyy_0[i] = 2.0 * g_zzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyz_0[i] = g_zzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxzz_0[i] = g_zzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyy_0[i] = 3.0 * g_zzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyz_0[i] = 2.0 * g_zzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyzz_0[i] = g_zzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxzzz_0[i] = g_zzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyy_0[i] = 4.0 * g_zzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyz_0[i] = 3.0 * g_zzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyzz_0[i] = 2.0 * g_zzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyzzz_0[i] = g_zzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxzzzz_0[i] = g_zzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyy_0[i] = 5.0 * g_zzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyz_0[i] = 4.0 * g_zzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyzz_0[i] = 3.0 * g_zzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyzzz_0[i] = 2.0 * g_zzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyzzzz_0[i] = g_zzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xzzzzz_0[i] = g_zzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyy_0[i] = 6.0 * g_zzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyz_0[i] = 5.0 * g_zzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyzz_0[i] = 4.0 * g_zzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyzzz_0[i] = 3.0 * g_zzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyzzzz_0[i] = 2.0 * g_zzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yzzzzz_0[i] = g_zzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_zzzzzz_0[i] = g_zzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 560-588 components of targeted buffer : HSI

    auto g_zzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 560);

    auto g_zzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 561);

    auto g_zzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 562);

    auto g_zzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 563);

    auto g_zzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 564);

    auto g_zzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 565);

    auto g_zzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 566);

    auto g_zzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 567);

    auto g_zzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 568);

    auto g_zzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 569);

    auto g_zzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 570);

    auto g_zzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 571);

    auto g_zzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 572);

    auto g_zzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 573);

    auto g_zzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 574);

    auto g_zzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 575);

    auto g_zzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 576);

    auto g_zzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 577);

    auto g_zzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 578);

    auto g_zzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 579);

    auto g_zzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 580);

    auto g_zzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 581);

    auto g_zzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 582);

    auto g_zzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 583);

    auto g_zzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 584);

    auto g_zzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 585);

    auto g_zzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 586);

    auto g_zzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 587);

    #pragma omp simd aligned(g_zzz_0_xxxxxx_0, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxy_0, g_zzz_0_xxxxxy_1, g_zzz_0_xxxxxz_0, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxyy_0, g_zzz_0_xxxxyy_1, g_zzz_0_xxxxyz_0, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxzz_0, g_zzz_0_xxxxzz_1, g_zzz_0_xxxyyy_0, g_zzz_0_xxxyyy_1, g_zzz_0_xxxyyz_0, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyzz_0, g_zzz_0_xxxyzz_1, g_zzz_0_xxxzzz_0, g_zzz_0_xxxzzz_1, g_zzz_0_xxyyyy_0, g_zzz_0_xxyyyy_1, g_zzz_0_xxyyyz_0, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyzz_0, g_zzz_0_xxyyzz_1, g_zzz_0_xxyzzz_0, g_zzz_0_xxyzzz_1, g_zzz_0_xxzzzz_0, g_zzz_0_xxzzzz_1, g_zzz_0_xyyyyy_0, g_zzz_0_xyyyyy_1, g_zzz_0_xyyyyz_0, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyzz_0, g_zzz_0_xyyyzz_1, g_zzz_0_xyyzzz_0, g_zzz_0_xyyzzz_1, g_zzz_0_xyzzzz_0, g_zzz_0_xyzzzz_1, g_zzz_0_xzzzzz_0, g_zzz_0_xzzzzz_1, g_zzz_0_yyyyyy_0, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyz_0, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyzz_0, g_zzz_0_yyyyzz_1, g_zzz_0_yyyzzz_0, g_zzz_0_yyyzzz_1, g_zzz_0_yyzzzz_0, g_zzz_0_yyzzzz_1, g_zzz_0_yzzzzz_0, g_zzz_0_yzzzzz_1, g_zzz_0_zzzzzz_0, g_zzz_0_zzzzzz_1, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxy_1, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxy_1, g_zzzz_0_xxxxyy_1, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxyy_1, g_zzzz_0_xxxyyy_1, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxzz_1, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxyyy_1, g_zzzz_0_xxyyyy_1, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyzz_1, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxzzz_1, g_zzzz_0_xxzzzz_1, g_zzzz_0_xyyyy_1, g_zzzz_0_xyyyyy_1, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyzz_1, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyzzz_1, g_zzzz_0_xyzzzz_1, g_zzzz_0_xzzzz_1, g_zzzz_0_xzzzzz_1, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyzz_1, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyzzz_1, g_zzzz_0_yyzzzz_1, g_zzzz_0_yzzzz_1, g_zzzz_0_yzzzzz_1, g_zzzz_0_zzzzz_1, g_zzzz_0_zzzzzz_1, g_zzzzz_0_xxxxxx_0, g_zzzzz_0_xxxxxy_0, g_zzzzz_0_xxxxxz_0, g_zzzzz_0_xxxxyy_0, g_zzzzz_0_xxxxyz_0, g_zzzzz_0_xxxxzz_0, g_zzzzz_0_xxxyyy_0, g_zzzzz_0_xxxyyz_0, g_zzzzz_0_xxxyzz_0, g_zzzzz_0_xxxzzz_0, g_zzzzz_0_xxyyyy_0, g_zzzzz_0_xxyyyz_0, g_zzzzz_0_xxyyzz_0, g_zzzzz_0_xxyzzz_0, g_zzzzz_0_xxzzzz_0, g_zzzzz_0_xyyyyy_0, g_zzzzz_0_xyyyyz_0, g_zzzzz_0_xyyyzz_0, g_zzzzz_0_xyyzzz_0, g_zzzzz_0_xyzzzz_0, g_zzzzz_0_xzzzzz_0, g_zzzzz_0_yyyyyy_0, g_zzzzz_0_yyyyyz_0, g_zzzzz_0_yyyyzz_0, g_zzzzz_0_yyyzzz_0, g_zzzzz_0_yyzzzz_0, g_zzzzz_0_yzzzzz_0, g_zzzzz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_xxxxxx_0[i] = 4.0 * g_zzz_0_xxxxxx_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxx_1[i] * fz_be_0 + g_zzzz_0_xxxxxx_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxy_0[i] = 4.0 * g_zzz_0_xxxxxy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxy_1[i] * fz_be_0 + g_zzzz_0_xxxxxy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxz_0[i] = 4.0 * g_zzz_0_xxxxxz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxz_1[i] * fz_be_0 + g_zzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyy_0[i] = 4.0 * g_zzz_0_xxxxyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyy_1[i] * fz_be_0 + g_zzzz_0_xxxxyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyz_0[i] = 4.0 * g_zzz_0_xxxxyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyz_1[i] * fz_be_0 + g_zzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxzz_0[i] = 4.0 * g_zzz_0_xxxxzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyy_0[i] = 4.0 * g_zzz_0_xxxyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyy_1[i] * fz_be_0 + g_zzzz_0_xxxyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyz_0[i] = 4.0 * g_zzz_0_xxxyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyz_1[i] * fz_be_0 + g_zzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyzz_0[i] = 4.0 * g_zzz_0_xxxyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxzzz_0[i] = 4.0 * g_zzz_0_xxxzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyy_0[i] = 4.0 * g_zzz_0_xxyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyy_1[i] * fz_be_0 + g_zzzz_0_xxyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyz_0[i] = 4.0 * g_zzz_0_xxyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyz_1[i] * fz_be_0 + g_zzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyzz_0[i] = 4.0 * g_zzz_0_xxyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyzzz_0[i] = 4.0 * g_zzz_0_xxyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxzzzz_0[i] = 4.0 * g_zzz_0_xxzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyy_0[i] = 4.0 * g_zzz_0_xyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyy_1[i] * fz_be_0 + g_zzzz_0_xyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyz_0[i] = 4.0 * g_zzz_0_xyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyz_1[i] * fz_be_0 + g_zzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyzz_0[i] = 4.0 * g_zzz_0_xyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyzzz_0[i] = 4.0 * g_zzz_0_xyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyzzzz_0[i] = 4.0 * g_zzz_0_xyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xzzzzz_0[i] = 4.0 * g_zzz_0_xzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyy_0[i] = 4.0 * g_zzz_0_yyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyy_1[i] * fz_be_0 + g_zzzz_0_yyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyz_0[i] = 4.0 * g_zzz_0_yyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyz_1[i] * fz_be_0 + g_zzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyzz_0[i] = 4.0 * g_zzz_0_yyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyzzz_0[i] = 4.0 * g_zzz_0_yyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzz_0_yyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyzzzz_0[i] = 4.0 * g_zzz_0_yyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yzzzzz_0[i] = 4.0 * g_zzz_0_yzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_zzzzzz_0[i] = 4.0 * g_zzz_0_zzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzz_0_zzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

