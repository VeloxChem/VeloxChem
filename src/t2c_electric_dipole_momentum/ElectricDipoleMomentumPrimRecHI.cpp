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

#include "ElectricDipoleMomentumPrimRecHI.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_hi(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_hi,
                                      const size_t              idx_dip_fi,
                                      const size_t              idx_dip_gh,
                                      const size_t              idx_ovl_gi,
                                      const size_t              idx_dip_gi,
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

    auto tr_x_xxx_xxxxxx = pbuffer.data(idx_dip_fi);

    auto tr_x_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 1);

    auto tr_x_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 2);

    auto tr_x_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 3);

    auto tr_x_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 4);

    auto tr_x_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 5);

    auto tr_x_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 6);

    auto tr_x_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 7);

    auto tr_x_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 8);

    auto tr_x_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 9);

    auto tr_x_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 10);

    auto tr_x_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 11);

    auto tr_x_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 12);

    auto tr_x_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 13);

    auto tr_x_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 14);

    auto tr_x_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 15);

    auto tr_x_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 16);

    auto tr_x_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 17);

    auto tr_x_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 18);

    auto tr_x_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 19);

    auto tr_x_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 20);

    auto tr_x_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 21);

    auto tr_x_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 22);

    auto tr_x_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 23);

    auto tr_x_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 24);

    auto tr_x_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 25);

    auto tr_x_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 26);

    auto tr_x_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 27);

    auto tr_x_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 28);

    auto tr_x_xxy_xxxxxy = pbuffer.data(idx_dip_fi + 29);

    auto tr_x_xxy_xxxxxz = pbuffer.data(idx_dip_fi + 30);

    auto tr_x_xxy_xxxxyy = pbuffer.data(idx_dip_fi + 31);

    auto tr_x_xxy_xxxxyz = pbuffer.data(idx_dip_fi + 32);

    auto tr_x_xxy_xxxxzz = pbuffer.data(idx_dip_fi + 33);

    auto tr_x_xxy_xxxyyy = pbuffer.data(idx_dip_fi + 34);

    auto tr_x_xxy_xxxyyz = pbuffer.data(idx_dip_fi + 35);

    auto tr_x_xxy_xxxyzz = pbuffer.data(idx_dip_fi + 36);

    auto tr_x_xxy_xxxzzz = pbuffer.data(idx_dip_fi + 37);

    auto tr_x_xxy_xxyyyy = pbuffer.data(idx_dip_fi + 38);

    auto tr_x_xxy_xxyyyz = pbuffer.data(idx_dip_fi + 39);

    auto tr_x_xxy_xxyyzz = pbuffer.data(idx_dip_fi + 40);

    auto tr_x_xxy_xxyzzz = pbuffer.data(idx_dip_fi + 41);

    auto tr_x_xxy_xxzzzz = pbuffer.data(idx_dip_fi + 42);

    auto tr_x_xxy_xyyyyy = pbuffer.data(idx_dip_fi + 43);

    auto tr_x_xxy_xyyyyz = pbuffer.data(idx_dip_fi + 44);

    auto tr_x_xxy_xyyyzz = pbuffer.data(idx_dip_fi + 45);

    auto tr_x_xxy_xyyzzz = pbuffer.data(idx_dip_fi + 46);

    auto tr_x_xxy_xyzzzz = pbuffer.data(idx_dip_fi + 47);

    auto tr_x_xxy_xzzzzz = pbuffer.data(idx_dip_fi + 48);

    auto tr_x_xxy_zzzzzz = pbuffer.data(idx_dip_fi + 55);

    auto tr_x_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 56);

    auto tr_x_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 57);

    auto tr_x_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 58);

    auto tr_x_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 59);

    auto tr_x_xxz_xxxxyz = pbuffer.data(idx_dip_fi + 60);

    auto tr_x_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 61);

    auto tr_x_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 62);

    auto tr_x_xxz_xxxyyz = pbuffer.data(idx_dip_fi + 63);

    auto tr_x_xxz_xxxyzz = pbuffer.data(idx_dip_fi + 64);

    auto tr_x_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 65);

    auto tr_x_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 66);

    auto tr_x_xxz_xxyyyz = pbuffer.data(idx_dip_fi + 67);

    auto tr_x_xxz_xxyyzz = pbuffer.data(idx_dip_fi + 68);

    auto tr_x_xxz_xxyzzz = pbuffer.data(idx_dip_fi + 69);

    auto tr_x_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 70);

    auto tr_x_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 71);

    auto tr_x_xxz_xyyyyz = pbuffer.data(idx_dip_fi + 72);

    auto tr_x_xxz_xyyyzz = pbuffer.data(idx_dip_fi + 73);

    auto tr_x_xxz_xyyzzz = pbuffer.data(idx_dip_fi + 74);

    auto tr_x_xxz_xyzzzz = pbuffer.data(idx_dip_fi + 75);

    auto tr_x_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 76);

    auto tr_x_xxz_yyyyyy = pbuffer.data(idx_dip_fi + 77);

    auto tr_x_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 83);

    auto tr_x_xyy_xxxxxx = pbuffer.data(idx_dip_fi + 84);

    auto tr_x_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 85);

    auto tr_x_xyy_xxxxxz = pbuffer.data(idx_dip_fi + 86);

    auto tr_x_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 87);

    auto tr_x_xyy_xxxxzz = pbuffer.data(idx_dip_fi + 89);

    auto tr_x_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 90);

    auto tr_x_xyy_xxxzzz = pbuffer.data(idx_dip_fi + 93);

    auto tr_x_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 94);

    auto tr_x_xyy_xxzzzz = pbuffer.data(idx_dip_fi + 98);

    auto tr_x_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 99);

    auto tr_x_xyy_xzzzzz = pbuffer.data(idx_dip_fi + 104);

    auto tr_x_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 105);

    auto tr_x_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 106);

    auto tr_x_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 107);

    auto tr_x_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 108);

    auto tr_x_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 109);

    auto tr_x_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 110);

    auto tr_x_xyz_xxxxxz = pbuffer.data(idx_dip_fi + 114);

    auto tr_x_xyz_xxxxzz = pbuffer.data(idx_dip_fi + 117);

    auto tr_x_xyz_xxxzzz = pbuffer.data(idx_dip_fi + 121);

    auto tr_x_xyz_xxzzzz = pbuffer.data(idx_dip_fi + 126);

    auto tr_x_xyz_xzzzzz = pbuffer.data(idx_dip_fi + 132);

    auto tr_x_xzz_xxxxxx = pbuffer.data(idx_dip_fi + 140);

    auto tr_x_xzz_xxxxxy = pbuffer.data(idx_dip_fi + 141);

    auto tr_x_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 142);

    auto tr_x_xzz_xxxxyy = pbuffer.data(idx_dip_fi + 143);

    auto tr_x_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 145);

    auto tr_x_xzz_xxxyyy = pbuffer.data(idx_dip_fi + 146);

    auto tr_x_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 149);

    auto tr_x_xzz_xxyyyy = pbuffer.data(idx_dip_fi + 150);

    auto tr_x_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 154);

    auto tr_x_xzz_xyyyyy = pbuffer.data(idx_dip_fi + 155);

    auto tr_x_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 160);

    auto tr_x_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 162);

    auto tr_x_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 163);

    auto tr_x_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 164);

    auto tr_x_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 165);

    auto tr_x_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 166);

    auto tr_x_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 167);

    auto tr_x_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 168);

    auto tr_x_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 169);

    auto tr_x_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 170);

    auto tr_x_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 171);

    auto tr_x_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 172);

    auto tr_x_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 173);

    auto tr_x_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 174);

    auto tr_x_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 175);

    auto tr_x_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 176);

    auto tr_x_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 177);

    auto tr_x_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 178);

    auto tr_x_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 179);

    auto tr_x_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 180);

    auto tr_x_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 181);

    auto tr_x_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 182);

    auto tr_x_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 183);

    auto tr_x_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 184);

    auto tr_x_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 185);

    auto tr_x_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 186);

    auto tr_x_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 187);

    auto tr_x_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 188);

    auto tr_x_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 189);

    auto tr_x_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 190);

    auto tr_x_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 191);

    auto tr_x_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 192);

    auto tr_x_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 193);

    auto tr_x_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 194);

    auto tr_x_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 195);

    auto tr_x_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 197);

    auto tr_x_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 198);

    auto tr_x_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 199);

    auto tr_x_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 201);

    auto tr_x_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 202);

    auto tr_x_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 205);

    auto tr_x_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 206);

    auto tr_x_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 210);

    auto tr_x_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 211);

    auto tr_x_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 216);

    auto tr_x_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 217);

    auto tr_x_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 223);

    auto tr_x_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 224);

    auto tr_x_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 226);

    auto tr_x_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 228);

    auto tr_x_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 229);

    auto tr_x_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 231);

    auto tr_x_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 232);

    auto tr_x_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 233);

    auto tr_x_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 235);

    auto tr_x_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 236);

    auto tr_x_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 237);

    auto tr_x_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 238);

    auto tr_x_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 240);

    auto tr_x_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 241);

    auto tr_x_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 242);

    auto tr_x_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 243);

    auto tr_x_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 244);

    auto tr_x_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 246);

    auto tr_x_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 247);

    auto tr_x_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 248);

    auto tr_x_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 249);

    auto tr_x_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 250);

    auto tr_x_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 251);

    auto tr_x_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 252);

    auto tr_x_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 253);

    auto tr_x_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 254);

    auto tr_x_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 255);

    auto tr_x_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 256);

    auto tr_x_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 257);

    auto tr_x_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 258);

    auto tr_x_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 259);

    auto tr_x_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 260);

    auto tr_x_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 261);

    auto tr_x_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 262);

    auto tr_x_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 263);

    auto tr_x_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 264);

    auto tr_x_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 265);

    auto tr_x_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 266);

    auto tr_x_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 267);

    auto tr_x_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 268);

    auto tr_x_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 269);

    auto tr_x_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 270);

    auto tr_x_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 271);

    auto tr_x_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 272);

    auto tr_x_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 273);

    auto tr_x_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 274);

    auto tr_x_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 275);

    auto tr_x_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 276);

    auto tr_x_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 277);

    auto tr_x_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 278);

    auto tr_x_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 279);

    auto tr_y_xxx_xxxxxx = pbuffer.data(idx_dip_fi + 280);

    auto tr_y_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 281);

    auto tr_y_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 282);

    auto tr_y_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 283);

    auto tr_y_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 284);

    auto tr_y_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 285);

    auto tr_y_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 286);

    auto tr_y_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 287);

    auto tr_y_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 288);

    auto tr_y_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 289);

    auto tr_y_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 290);

    auto tr_y_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 291);

    auto tr_y_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 292);

    auto tr_y_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 293);

    auto tr_y_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 294);

    auto tr_y_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 295);

    auto tr_y_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 296);

    auto tr_y_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 297);

    auto tr_y_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 298);

    auto tr_y_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 299);

    auto tr_y_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 300);

    auto tr_y_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 301);

    auto tr_y_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 302);

    auto tr_y_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 303);

    auto tr_y_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 304);

    auto tr_y_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 305);

    auto tr_y_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 306);

    auto tr_y_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 307);

    auto tr_y_xxy_xxxxxy = pbuffer.data(idx_dip_fi + 309);

    auto tr_y_xxy_xxxxyy = pbuffer.data(idx_dip_fi + 311);

    auto tr_y_xxy_xxxxyz = pbuffer.data(idx_dip_fi + 312);

    auto tr_y_xxy_xxxyyy = pbuffer.data(idx_dip_fi + 314);

    auto tr_y_xxy_xxxyyz = pbuffer.data(idx_dip_fi + 315);

    auto tr_y_xxy_xxxyzz = pbuffer.data(idx_dip_fi + 316);

    auto tr_y_xxy_xxyyyy = pbuffer.data(idx_dip_fi + 318);

    auto tr_y_xxy_xxyyyz = pbuffer.data(idx_dip_fi + 319);

    auto tr_y_xxy_xxyyzz = pbuffer.data(idx_dip_fi + 320);

    auto tr_y_xxy_xxyzzz = pbuffer.data(idx_dip_fi + 321);

    auto tr_y_xxy_xyyyyy = pbuffer.data(idx_dip_fi + 323);

    auto tr_y_xxy_xyyyyz = pbuffer.data(idx_dip_fi + 324);

    auto tr_y_xxy_xyyyzz = pbuffer.data(idx_dip_fi + 325);

    auto tr_y_xxy_xyyzzz = pbuffer.data(idx_dip_fi + 326);

    auto tr_y_xxy_xyzzzz = pbuffer.data(idx_dip_fi + 327);

    auto tr_y_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 329);

    auto tr_y_xxy_yyyyyz = pbuffer.data(idx_dip_fi + 330);

    auto tr_y_xxy_yyyyzz = pbuffer.data(idx_dip_fi + 331);

    auto tr_y_xxy_yyyzzz = pbuffer.data(idx_dip_fi + 332);

    auto tr_y_xxy_yyzzzz = pbuffer.data(idx_dip_fi + 333);

    auto tr_y_xxy_yzzzzz = pbuffer.data(idx_dip_fi + 334);

    auto tr_y_xxy_zzzzzz = pbuffer.data(idx_dip_fi + 335);

    auto tr_y_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 336);

    auto tr_y_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 337);

    auto tr_y_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 339);

    auto tr_y_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 342);

    auto tr_y_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 346);

    auto tr_y_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 351);

    auto tr_y_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 358);

    auto tr_y_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 359);

    auto tr_y_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 360);

    auto tr_y_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 361);

    auto tr_y_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 362);

    auto tr_y_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 363);

    auto tr_y_xyy_xxxxxx = pbuffer.data(idx_dip_fi + 364);

    auto tr_y_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 365);

    auto tr_y_xyy_xxxxxz = pbuffer.data(idx_dip_fi + 366);

    auto tr_y_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 367);

    auto tr_y_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 368);

    auto tr_y_xyy_xxxxzz = pbuffer.data(idx_dip_fi + 369);

    auto tr_y_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 370);

    auto tr_y_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 371);

    auto tr_y_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 372);

    auto tr_y_xyy_xxxzzz = pbuffer.data(idx_dip_fi + 373);

    auto tr_y_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 374);

    auto tr_y_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 375);

    auto tr_y_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 376);

    auto tr_y_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 377);

    auto tr_y_xyy_xxzzzz = pbuffer.data(idx_dip_fi + 378);

    auto tr_y_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 379);

    auto tr_y_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 380);

    auto tr_y_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 381);

    auto tr_y_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 382);

    auto tr_y_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 383);

    auto tr_y_xyy_xzzzzz = pbuffer.data(idx_dip_fi + 384);

    auto tr_y_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 385);

    auto tr_y_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 386);

    auto tr_y_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 387);

    auto tr_y_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 388);

    auto tr_y_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 389);

    auto tr_y_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 390);

    auto tr_y_xyy_zzzzzz = pbuffer.data(idx_dip_fi + 391);

    auto tr_y_xyz_yyyyyz = pbuffer.data(idx_dip_fi + 414);

    auto tr_y_xyz_yyyyzz = pbuffer.data(idx_dip_fi + 415);

    auto tr_y_xyz_yyyzzz = pbuffer.data(idx_dip_fi + 416);

    auto tr_y_xyz_yyzzzz = pbuffer.data(idx_dip_fi + 417);

    auto tr_y_xyz_yzzzzz = pbuffer.data(idx_dip_fi + 418);

    auto tr_y_xyz_zzzzzz = pbuffer.data(idx_dip_fi + 419);

    auto tr_y_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 422);

    auto tr_y_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 424);

    auto tr_y_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 425);

    auto tr_y_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 427);

    auto tr_y_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 428);

    auto tr_y_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 429);

    auto tr_y_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 431);

    auto tr_y_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 432);

    auto tr_y_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 433);

    auto tr_y_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 434);

    auto tr_y_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 436);

    auto tr_y_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 437);

    auto tr_y_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 438);

    auto tr_y_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 439);

    auto tr_y_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 440);

    auto tr_y_xzz_yyyyyy = pbuffer.data(idx_dip_fi + 441);

    auto tr_y_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 442);

    auto tr_y_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 443);

    auto tr_y_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 444);

    auto tr_y_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 445);

    auto tr_y_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 446);

    auto tr_y_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 447);

    auto tr_y_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 448);

    auto tr_y_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 449);

    auto tr_y_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 450);

    auto tr_y_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 451);

    auto tr_y_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 452);

    auto tr_y_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 453);

    auto tr_y_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 454);

    auto tr_y_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 455);

    auto tr_y_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 456);

    auto tr_y_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 457);

    auto tr_y_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 458);

    auto tr_y_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 459);

    auto tr_y_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 460);

    auto tr_y_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 461);

    auto tr_y_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 462);

    auto tr_y_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 463);

    auto tr_y_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 464);

    auto tr_y_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 465);

    auto tr_y_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 466);

    auto tr_y_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 467);

    auto tr_y_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 468);

    auto tr_y_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 469);

    auto tr_y_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 470);

    auto tr_y_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 471);

    auto tr_y_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 472);

    auto tr_y_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 473);

    auto tr_y_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 474);

    auto tr_y_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 475);

    auto tr_y_yyz_xxxxxx = pbuffer.data(idx_dip_fi + 476);

    auto tr_y_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 477);

    auto tr_y_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 479);

    auto tr_y_yyz_xxxxyz = pbuffer.data(idx_dip_fi + 480);

    auto tr_y_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 482);

    auto tr_y_yyz_xxxyyz = pbuffer.data(idx_dip_fi + 483);

    auto tr_y_yyz_xxxyzz = pbuffer.data(idx_dip_fi + 484);

    auto tr_y_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 486);

    auto tr_y_yyz_xxyyyz = pbuffer.data(idx_dip_fi + 487);

    auto tr_y_yyz_xxyyzz = pbuffer.data(idx_dip_fi + 488);

    auto tr_y_yyz_xxyzzz = pbuffer.data(idx_dip_fi + 489);

    auto tr_y_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 491);

    auto tr_y_yyz_xyyyyz = pbuffer.data(idx_dip_fi + 492);

    auto tr_y_yyz_xyyyzz = pbuffer.data(idx_dip_fi + 493);

    auto tr_y_yyz_xyyzzz = pbuffer.data(idx_dip_fi + 494);

    auto tr_y_yyz_xyzzzz = pbuffer.data(idx_dip_fi + 495);

    auto tr_y_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 497);

    auto tr_y_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 498);

    auto tr_y_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 499);

    auto tr_y_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 500);

    auto tr_y_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 501);

    auto tr_y_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 502);

    auto tr_y_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 503);

    auto tr_y_yzz_xxxxxy = pbuffer.data(idx_dip_fi + 505);

    auto tr_y_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 506);

    auto tr_y_yzz_xxxxyy = pbuffer.data(idx_dip_fi + 507);

    auto tr_y_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 508);

    auto tr_y_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 509);

    auto tr_y_yzz_xxxyyy = pbuffer.data(idx_dip_fi + 510);

    auto tr_y_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 511);

    auto tr_y_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 512);

    auto tr_y_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 513);

    auto tr_y_yzz_xxyyyy = pbuffer.data(idx_dip_fi + 514);

    auto tr_y_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 515);

    auto tr_y_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 516);

    auto tr_y_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 517);

    auto tr_y_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 518);

    auto tr_y_yzz_xyyyyy = pbuffer.data(idx_dip_fi + 519);

    auto tr_y_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 520);

    auto tr_y_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 521);

    auto tr_y_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 522);

    auto tr_y_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 523);

    auto tr_y_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 524);

    auto tr_y_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 525);

    auto tr_y_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 526);

    auto tr_y_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 527);

    auto tr_y_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 528);

    auto tr_y_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 529);

    auto tr_y_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 530);

    auto tr_y_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 531);

    auto tr_y_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 532);

    auto tr_y_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 533);

    auto tr_y_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 534);

    auto tr_y_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 535);

    auto tr_y_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 536);

    auto tr_y_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 537);

    auto tr_y_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 538);

    auto tr_y_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 539);

    auto tr_y_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 540);

    auto tr_y_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 541);

    auto tr_y_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 542);

    auto tr_y_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 543);

    auto tr_y_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 544);

    auto tr_y_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 545);

    auto tr_y_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 546);

    auto tr_y_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 547);

    auto tr_y_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 548);

    auto tr_y_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 549);

    auto tr_y_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 550);

    auto tr_y_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 551);

    auto tr_y_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 552);

    auto tr_y_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 553);

    auto tr_y_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 554);

    auto tr_y_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 555);

    auto tr_y_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 556);

    auto tr_y_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 557);

    auto tr_y_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 558);

    auto tr_y_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 559);

    auto tr_z_xxx_xxxxxx = pbuffer.data(idx_dip_fi + 560);

    auto tr_z_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 561);

    auto tr_z_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 562);

    auto tr_z_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 563);

    auto tr_z_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 564);

    auto tr_z_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 565);

    auto tr_z_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 566);

    auto tr_z_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 567);

    auto tr_z_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 568);

    auto tr_z_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 569);

    auto tr_z_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 570);

    auto tr_z_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 571);

    auto tr_z_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 572);

    auto tr_z_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 573);

    auto tr_z_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 574);

    auto tr_z_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 575);

    auto tr_z_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 576);

    auto tr_z_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 577);

    auto tr_z_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 578);

    auto tr_z_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 579);

    auto tr_z_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 580);

    auto tr_z_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 581);

    auto tr_z_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 582);

    auto tr_z_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 583);

    auto tr_z_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 584);

    auto tr_z_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 585);

    auto tr_z_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 586);

    auto tr_z_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 587);

    auto tr_z_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 588);

    auto tr_z_xxy_xxxxxz = pbuffer.data(idx_dip_fi + 590);

    auto tr_z_xxy_xxxxzz = pbuffer.data(idx_dip_fi + 593);

    auto tr_z_xxy_xxxzzz = pbuffer.data(idx_dip_fi + 597);

    auto tr_z_xxy_xxzzzz = pbuffer.data(idx_dip_fi + 602);

    auto tr_z_xxy_xzzzzz = pbuffer.data(idx_dip_fi + 608);

    auto tr_z_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 609);

    auto tr_z_xxy_yyyyyz = pbuffer.data(idx_dip_fi + 610);

    auto tr_z_xxy_yyyyzz = pbuffer.data(idx_dip_fi + 611);

    auto tr_z_xxy_yyyzzz = pbuffer.data(idx_dip_fi + 612);

    auto tr_z_xxy_yyzzzz = pbuffer.data(idx_dip_fi + 613);

    auto tr_z_xxy_yzzzzz = pbuffer.data(idx_dip_fi + 614);

    auto tr_z_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 616);

    auto tr_z_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 618);

    auto tr_z_xxz_xxxxyz = pbuffer.data(idx_dip_fi + 620);

    auto tr_z_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 621);

    auto tr_z_xxz_xxxyyz = pbuffer.data(idx_dip_fi + 623);

    auto tr_z_xxz_xxxyzz = pbuffer.data(idx_dip_fi + 624);

    auto tr_z_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 625);

    auto tr_z_xxz_xxyyyz = pbuffer.data(idx_dip_fi + 627);

    auto tr_z_xxz_xxyyzz = pbuffer.data(idx_dip_fi + 628);

    auto tr_z_xxz_xxyzzz = pbuffer.data(idx_dip_fi + 629);

    auto tr_z_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 630);

    auto tr_z_xxz_xyyyyz = pbuffer.data(idx_dip_fi + 632);

    auto tr_z_xxz_xyyyzz = pbuffer.data(idx_dip_fi + 633);

    auto tr_z_xxz_xyyzzz = pbuffer.data(idx_dip_fi + 634);

    auto tr_z_xxz_xyzzzz = pbuffer.data(idx_dip_fi + 635);

    auto tr_z_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 636);

    auto tr_z_xxz_yyyyyy = pbuffer.data(idx_dip_fi + 637);

    auto tr_z_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 638);

    auto tr_z_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 639);

    auto tr_z_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 640);

    auto tr_z_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 641);

    auto tr_z_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 642);

    auto tr_z_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 643);

    auto tr_z_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 645);

    auto tr_z_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 647);

    auto tr_z_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 648);

    auto tr_z_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 650);

    auto tr_z_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 651);

    auto tr_z_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 652);

    auto tr_z_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 654);

    auto tr_z_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 655);

    auto tr_z_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 656);

    auto tr_z_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 657);

    auto tr_z_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 659);

    auto tr_z_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 660);

    auto tr_z_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 661);

    auto tr_z_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 662);

    auto tr_z_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 663);

    auto tr_z_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 665);

    auto tr_z_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 666);

    auto tr_z_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 667);

    auto tr_z_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 668);

    auto tr_z_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 669);

    auto tr_z_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 670);

    auto tr_z_xyy_zzzzzz = pbuffer.data(idx_dip_fi + 671);

    auto tr_z_xyz_yyyyyy = pbuffer.data(idx_dip_fi + 693);

    auto tr_z_xyz_yyyyyz = pbuffer.data(idx_dip_fi + 694);

    auto tr_z_xyz_yyyyzz = pbuffer.data(idx_dip_fi + 695);

    auto tr_z_xyz_yyyzzz = pbuffer.data(idx_dip_fi + 696);

    auto tr_z_xyz_yyzzzz = pbuffer.data(idx_dip_fi + 697);

    auto tr_z_xyz_yzzzzz = pbuffer.data(idx_dip_fi + 698);

    auto tr_z_xzz_xxxxxx = pbuffer.data(idx_dip_fi + 700);

    auto tr_z_xzz_xxxxxy = pbuffer.data(idx_dip_fi + 701);

    auto tr_z_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 702);

    auto tr_z_xzz_xxxxyy = pbuffer.data(idx_dip_fi + 703);

    auto tr_z_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 704);

    auto tr_z_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 705);

    auto tr_z_xzz_xxxyyy = pbuffer.data(idx_dip_fi + 706);

    auto tr_z_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 707);

    auto tr_z_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 708);

    auto tr_z_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 709);

    auto tr_z_xzz_xxyyyy = pbuffer.data(idx_dip_fi + 710);

    auto tr_z_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 711);

    auto tr_z_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 712);

    auto tr_z_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 713);

    auto tr_z_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 714);

    auto tr_z_xzz_xyyyyy = pbuffer.data(idx_dip_fi + 715);

    auto tr_z_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 716);

    auto tr_z_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 717);

    auto tr_z_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 718);

    auto tr_z_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 719);

    auto tr_z_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 720);

    auto tr_z_xzz_yyyyyy = pbuffer.data(idx_dip_fi + 721);

    auto tr_z_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 722);

    auto tr_z_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 723);

    auto tr_z_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 724);

    auto tr_z_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 725);

    auto tr_z_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 726);

    auto tr_z_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 727);

    auto tr_z_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 728);

    auto tr_z_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 729);

    auto tr_z_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 730);

    auto tr_z_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 731);

    auto tr_z_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 732);

    auto tr_z_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 733);

    auto tr_z_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 734);

    auto tr_z_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 735);

    auto tr_z_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 736);

    auto tr_z_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 737);

    auto tr_z_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 738);

    auto tr_z_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 739);

    auto tr_z_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 740);

    auto tr_z_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 741);

    auto tr_z_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 742);

    auto tr_z_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 743);

    auto tr_z_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 744);

    auto tr_z_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 745);

    auto tr_z_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 746);

    auto tr_z_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 747);

    auto tr_z_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 748);

    auto tr_z_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 749);

    auto tr_z_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 750);

    auto tr_z_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 751);

    auto tr_z_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 752);

    auto tr_z_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 753);

    auto tr_z_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 754);

    auto tr_z_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 755);

    auto tr_z_yyz_xxxxxx = pbuffer.data(idx_dip_fi + 756);

    auto tr_z_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 758);

    auto tr_z_yyz_xxxxyz = pbuffer.data(idx_dip_fi + 760);

    auto tr_z_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 761);

    auto tr_z_yyz_xxxyyz = pbuffer.data(idx_dip_fi + 763);

    auto tr_z_yyz_xxxyzz = pbuffer.data(idx_dip_fi + 764);

    auto tr_z_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 765);

    auto tr_z_yyz_xxyyyz = pbuffer.data(idx_dip_fi + 767);

    auto tr_z_yyz_xxyyzz = pbuffer.data(idx_dip_fi + 768);

    auto tr_z_yyz_xxyzzz = pbuffer.data(idx_dip_fi + 769);

    auto tr_z_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 770);

    auto tr_z_yyz_xyyyyz = pbuffer.data(idx_dip_fi + 772);

    auto tr_z_yyz_xyyyzz = pbuffer.data(idx_dip_fi + 773);

    auto tr_z_yyz_xyyzzz = pbuffer.data(idx_dip_fi + 774);

    auto tr_z_yyz_xyzzzz = pbuffer.data(idx_dip_fi + 775);

    auto tr_z_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 776);

    auto tr_z_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 777);

    auto tr_z_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 778);

    auto tr_z_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 779);

    auto tr_z_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 780);

    auto tr_z_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 781);

    auto tr_z_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 782);

    auto tr_z_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 783);

    auto tr_z_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 784);

    auto tr_z_yzz_xxxxxy = pbuffer.data(idx_dip_fi + 785);

    auto tr_z_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 786);

    auto tr_z_yzz_xxxxyy = pbuffer.data(idx_dip_fi + 787);

    auto tr_z_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 788);

    auto tr_z_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 789);

    auto tr_z_yzz_xxxyyy = pbuffer.data(idx_dip_fi + 790);

    auto tr_z_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 791);

    auto tr_z_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 792);

    auto tr_z_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 793);

    auto tr_z_yzz_xxyyyy = pbuffer.data(idx_dip_fi + 794);

    auto tr_z_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 795);

    auto tr_z_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 796);

    auto tr_z_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 797);

    auto tr_z_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 798);

    auto tr_z_yzz_xyyyyy = pbuffer.data(idx_dip_fi + 799);

    auto tr_z_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 800);

    auto tr_z_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 801);

    auto tr_z_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 802);

    auto tr_z_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 803);

    auto tr_z_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 804);

    auto tr_z_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 805);

    auto tr_z_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 806);

    auto tr_z_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 807);

    auto tr_z_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 808);

    auto tr_z_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 809);

    auto tr_z_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 810);

    auto tr_z_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 811);

    auto tr_z_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 812);

    auto tr_z_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 813);

    auto tr_z_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 814);

    auto tr_z_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 815);

    auto tr_z_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 816);

    auto tr_z_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 817);

    auto tr_z_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 818);

    auto tr_z_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 819);

    auto tr_z_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 820);

    auto tr_z_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 821);

    auto tr_z_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 822);

    auto tr_z_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 823);

    auto tr_z_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 824);

    auto tr_z_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 825);

    auto tr_z_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 826);

    auto tr_z_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 827);

    auto tr_z_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 828);

    auto tr_z_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 829);

    auto tr_z_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 830);

    auto tr_z_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 831);

    auto tr_z_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 832);

    auto tr_z_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 833);

    auto tr_z_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 834);

    auto tr_z_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 835);

    auto tr_z_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 836);

    auto tr_z_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 837);

    auto tr_z_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 838);

    auto tr_z_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 839);

    // Set up components of auxiliary buffer : GH

    auto tr_x_xxxx_xxxxx = pbuffer.data(idx_dip_gh);

    auto tr_x_xxxx_xxxxy = pbuffer.data(idx_dip_gh + 1);

    auto tr_x_xxxx_xxxxz = pbuffer.data(idx_dip_gh + 2);

    auto tr_x_xxxx_xxxyy = pbuffer.data(idx_dip_gh + 3);

    auto tr_x_xxxx_xxxyz = pbuffer.data(idx_dip_gh + 4);

    auto tr_x_xxxx_xxxzz = pbuffer.data(idx_dip_gh + 5);

    auto tr_x_xxxx_xxyyy = pbuffer.data(idx_dip_gh + 6);

    auto tr_x_xxxx_xxyyz = pbuffer.data(idx_dip_gh + 7);

    auto tr_x_xxxx_xxyzz = pbuffer.data(idx_dip_gh + 8);

    auto tr_x_xxxx_xxzzz = pbuffer.data(idx_dip_gh + 9);

    auto tr_x_xxxx_xyyyy = pbuffer.data(idx_dip_gh + 10);

    auto tr_x_xxxx_xyyyz = pbuffer.data(idx_dip_gh + 11);

    auto tr_x_xxxx_xyyzz = pbuffer.data(idx_dip_gh + 12);

    auto tr_x_xxxx_xyzzz = pbuffer.data(idx_dip_gh + 13);

    auto tr_x_xxxx_xzzzz = pbuffer.data(idx_dip_gh + 14);

    auto tr_x_xxxx_yyyyy = pbuffer.data(idx_dip_gh + 15);

    auto tr_x_xxxx_yyyyz = pbuffer.data(idx_dip_gh + 16);

    auto tr_x_xxxx_yyyzz = pbuffer.data(idx_dip_gh + 17);

    auto tr_x_xxxx_yyzzz = pbuffer.data(idx_dip_gh + 18);

    auto tr_x_xxxx_yzzzz = pbuffer.data(idx_dip_gh + 19);

    auto tr_x_xxxx_zzzzz = pbuffer.data(idx_dip_gh + 20);

    auto tr_x_xxxy_xxxxx = pbuffer.data(idx_dip_gh + 21);

    auto tr_x_xxxy_xxxxy = pbuffer.data(idx_dip_gh + 22);

    auto tr_x_xxxy_xxxxz = pbuffer.data(idx_dip_gh + 23);

    auto tr_x_xxxy_xxxyy = pbuffer.data(idx_dip_gh + 24);

    auto tr_x_xxxy_xxxyz = pbuffer.data(idx_dip_gh + 25);

    auto tr_x_xxxy_xxxzz = pbuffer.data(idx_dip_gh + 26);

    auto tr_x_xxxy_xxyyy = pbuffer.data(idx_dip_gh + 27);

    auto tr_x_xxxy_xxyyz = pbuffer.data(idx_dip_gh + 28);

    auto tr_x_xxxy_xxyzz = pbuffer.data(idx_dip_gh + 29);

    auto tr_x_xxxy_xxzzz = pbuffer.data(idx_dip_gh + 30);

    auto tr_x_xxxy_xyyyy = pbuffer.data(idx_dip_gh + 31);

    auto tr_x_xxxy_xyyyz = pbuffer.data(idx_dip_gh + 32);

    auto tr_x_xxxy_xyyzz = pbuffer.data(idx_dip_gh + 33);

    auto tr_x_xxxy_xyzzz = pbuffer.data(idx_dip_gh + 34);

    auto tr_x_xxxy_xzzzz = pbuffer.data(idx_dip_gh + 35);

    auto tr_x_xxxz_xxxxx = pbuffer.data(idx_dip_gh + 42);

    auto tr_x_xxxz_xxxxy = pbuffer.data(idx_dip_gh + 43);

    auto tr_x_xxxz_xxxxz = pbuffer.data(idx_dip_gh + 44);

    auto tr_x_xxxz_xxxyy = pbuffer.data(idx_dip_gh + 45);

    auto tr_x_xxxz_xxxyz = pbuffer.data(idx_dip_gh + 46);

    auto tr_x_xxxz_xxxzz = pbuffer.data(idx_dip_gh + 47);

    auto tr_x_xxxz_xxyyy = pbuffer.data(idx_dip_gh + 48);

    auto tr_x_xxxz_xxyyz = pbuffer.data(idx_dip_gh + 49);

    auto tr_x_xxxz_xxyzz = pbuffer.data(idx_dip_gh + 50);

    auto tr_x_xxxz_xxzzz = pbuffer.data(idx_dip_gh + 51);

    auto tr_x_xxxz_xyyyy = pbuffer.data(idx_dip_gh + 52);

    auto tr_x_xxxz_xyyyz = pbuffer.data(idx_dip_gh + 53);

    auto tr_x_xxxz_xyyzz = pbuffer.data(idx_dip_gh + 54);

    auto tr_x_xxxz_xyzzz = pbuffer.data(idx_dip_gh + 55);

    auto tr_x_xxxz_xzzzz = pbuffer.data(idx_dip_gh + 56);

    auto tr_x_xxxz_yyyyz = pbuffer.data(idx_dip_gh + 58);

    auto tr_x_xxxz_yyyzz = pbuffer.data(idx_dip_gh + 59);

    auto tr_x_xxxz_yyzzz = pbuffer.data(idx_dip_gh + 60);

    auto tr_x_xxxz_yzzzz = pbuffer.data(idx_dip_gh + 61);

    auto tr_x_xxxz_zzzzz = pbuffer.data(idx_dip_gh + 62);

    auto tr_x_xxyy_xxxxx = pbuffer.data(idx_dip_gh + 63);

    auto tr_x_xxyy_xxxxy = pbuffer.data(idx_dip_gh + 64);

    auto tr_x_xxyy_xxxxz = pbuffer.data(idx_dip_gh + 65);

    auto tr_x_xxyy_xxxyy = pbuffer.data(idx_dip_gh + 66);

    auto tr_x_xxyy_xxxyz = pbuffer.data(idx_dip_gh + 67);

    auto tr_x_xxyy_xxxzz = pbuffer.data(idx_dip_gh + 68);

    auto tr_x_xxyy_xxyyy = pbuffer.data(idx_dip_gh + 69);

    auto tr_x_xxyy_xxyyz = pbuffer.data(idx_dip_gh + 70);

    auto tr_x_xxyy_xxyzz = pbuffer.data(idx_dip_gh + 71);

    auto tr_x_xxyy_xxzzz = pbuffer.data(idx_dip_gh + 72);

    auto tr_x_xxyy_xyyyy = pbuffer.data(idx_dip_gh + 73);

    auto tr_x_xxyy_xyyyz = pbuffer.data(idx_dip_gh + 74);

    auto tr_x_xxyy_xyyzz = pbuffer.data(idx_dip_gh + 75);

    auto tr_x_xxyy_xyzzz = pbuffer.data(idx_dip_gh + 76);

    auto tr_x_xxyy_xzzzz = pbuffer.data(idx_dip_gh + 77);

    auto tr_x_xxyy_yyyyy = pbuffer.data(idx_dip_gh + 78);

    auto tr_x_xxyy_yyyyz = pbuffer.data(idx_dip_gh + 79);

    auto tr_x_xxyy_yyyzz = pbuffer.data(idx_dip_gh + 80);

    auto tr_x_xxyy_yyzzz = pbuffer.data(idx_dip_gh + 81);

    auto tr_x_xxyy_yzzzz = pbuffer.data(idx_dip_gh + 82);

    auto tr_x_xxzz_xxxxx = pbuffer.data(idx_dip_gh + 105);

    auto tr_x_xxzz_xxxxy = pbuffer.data(idx_dip_gh + 106);

    auto tr_x_xxzz_xxxxz = pbuffer.data(idx_dip_gh + 107);

    auto tr_x_xxzz_xxxyy = pbuffer.data(idx_dip_gh + 108);

    auto tr_x_xxzz_xxxyz = pbuffer.data(idx_dip_gh + 109);

    auto tr_x_xxzz_xxxzz = pbuffer.data(idx_dip_gh + 110);

    auto tr_x_xxzz_xxyyy = pbuffer.data(idx_dip_gh + 111);

    auto tr_x_xxzz_xxyyz = pbuffer.data(idx_dip_gh + 112);

    auto tr_x_xxzz_xxyzz = pbuffer.data(idx_dip_gh + 113);

    auto tr_x_xxzz_xxzzz = pbuffer.data(idx_dip_gh + 114);

    auto tr_x_xxzz_xyyyy = pbuffer.data(idx_dip_gh + 115);

    auto tr_x_xxzz_xyyyz = pbuffer.data(idx_dip_gh + 116);

    auto tr_x_xxzz_xyyzz = pbuffer.data(idx_dip_gh + 117);

    auto tr_x_xxzz_xyzzz = pbuffer.data(idx_dip_gh + 118);

    auto tr_x_xxzz_xzzzz = pbuffer.data(idx_dip_gh + 119);

    auto tr_x_xxzz_yyyyy = pbuffer.data(idx_dip_gh + 120);

    auto tr_x_xxzz_yyyyz = pbuffer.data(idx_dip_gh + 121);

    auto tr_x_xxzz_yyyzz = pbuffer.data(idx_dip_gh + 122);

    auto tr_x_xxzz_yyzzz = pbuffer.data(idx_dip_gh + 123);

    auto tr_x_xxzz_yzzzz = pbuffer.data(idx_dip_gh + 124);

    auto tr_x_xxzz_zzzzz = pbuffer.data(idx_dip_gh + 125);

    auto tr_x_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 127);

    auto tr_x_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 129);

    auto tr_x_xyyy_xxxyz = pbuffer.data(idx_dip_gh + 130);

    auto tr_x_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 132);

    auto tr_x_xyyy_xxyyz = pbuffer.data(idx_dip_gh + 133);

    auto tr_x_xyyy_xxyzz = pbuffer.data(idx_dip_gh + 134);

    auto tr_x_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 136);

    auto tr_x_xyyy_xyyyz = pbuffer.data(idx_dip_gh + 137);

    auto tr_x_xyyy_xyyzz = pbuffer.data(idx_dip_gh + 138);

    auto tr_x_xyyy_xyzzz = pbuffer.data(idx_dip_gh + 139);

    auto tr_x_xzzz_xxxxx = pbuffer.data(idx_dip_gh + 189);

    auto tr_x_xzzz_xxxxy = pbuffer.data(idx_dip_gh + 190);

    auto tr_x_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 191);

    auto tr_x_xzzz_xxxyy = pbuffer.data(idx_dip_gh + 192);

    auto tr_x_xzzz_xxxyz = pbuffer.data(idx_dip_gh + 193);

    auto tr_x_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 194);

    auto tr_x_xzzz_xxyyy = pbuffer.data(idx_dip_gh + 195);

    auto tr_x_xzzz_xxyyz = pbuffer.data(idx_dip_gh + 196);

    auto tr_x_xzzz_xxyzz = pbuffer.data(idx_dip_gh + 197);

    auto tr_x_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 198);

    auto tr_x_xzzz_xyyyy = pbuffer.data(idx_dip_gh + 199);

    auto tr_x_xzzz_xyyyz = pbuffer.data(idx_dip_gh + 200);

    auto tr_x_xzzz_xyyzz = pbuffer.data(idx_dip_gh + 201);

    auto tr_x_xzzz_xyzzz = pbuffer.data(idx_dip_gh + 202);

    auto tr_x_xzzz_xzzzz = pbuffer.data(idx_dip_gh + 203);

    auto tr_x_yyyy_xxxxx = pbuffer.data(idx_dip_gh + 210);

    auto tr_x_yyyy_xxxxy = pbuffer.data(idx_dip_gh + 211);

    auto tr_x_yyyy_xxxxz = pbuffer.data(idx_dip_gh + 212);

    auto tr_x_yyyy_xxxyy = pbuffer.data(idx_dip_gh + 213);

    auto tr_x_yyyy_xxxyz = pbuffer.data(idx_dip_gh + 214);

    auto tr_x_yyyy_xxxzz = pbuffer.data(idx_dip_gh + 215);

    auto tr_x_yyyy_xxyyy = pbuffer.data(idx_dip_gh + 216);

    auto tr_x_yyyy_xxyyz = pbuffer.data(idx_dip_gh + 217);

    auto tr_x_yyyy_xxyzz = pbuffer.data(idx_dip_gh + 218);

    auto tr_x_yyyy_xxzzz = pbuffer.data(idx_dip_gh + 219);

    auto tr_x_yyyy_xyyyy = pbuffer.data(idx_dip_gh + 220);

    auto tr_x_yyyy_xyyyz = pbuffer.data(idx_dip_gh + 221);

    auto tr_x_yyyy_xyyzz = pbuffer.data(idx_dip_gh + 222);

    auto tr_x_yyyy_xyzzz = pbuffer.data(idx_dip_gh + 223);

    auto tr_x_yyyy_xzzzz = pbuffer.data(idx_dip_gh + 224);

    auto tr_x_yyyy_yyyyy = pbuffer.data(idx_dip_gh + 225);

    auto tr_x_yyyy_yyyyz = pbuffer.data(idx_dip_gh + 226);

    auto tr_x_yyyy_yyyzz = pbuffer.data(idx_dip_gh + 227);

    auto tr_x_yyyy_yyzzz = pbuffer.data(idx_dip_gh + 228);

    auto tr_x_yyyy_yzzzz = pbuffer.data(idx_dip_gh + 229);

    auto tr_x_yyyy_zzzzz = pbuffer.data(idx_dip_gh + 230);

    auto tr_x_yyzz_xxxxz = pbuffer.data(idx_dip_gh + 254);

    auto tr_x_yyzz_xxxyz = pbuffer.data(idx_dip_gh + 256);

    auto tr_x_yyzz_xxxzz = pbuffer.data(idx_dip_gh + 257);

    auto tr_x_yyzz_xxyyz = pbuffer.data(idx_dip_gh + 259);

    auto tr_x_yyzz_xxyzz = pbuffer.data(idx_dip_gh + 260);

    auto tr_x_yyzz_xxzzz = pbuffer.data(idx_dip_gh + 261);

    auto tr_x_yyzz_xyyyz = pbuffer.data(idx_dip_gh + 263);

    auto tr_x_yyzz_xyyzz = pbuffer.data(idx_dip_gh + 264);

    auto tr_x_yyzz_xyzzz = pbuffer.data(idx_dip_gh + 265);

    auto tr_x_yyzz_xzzzz = pbuffer.data(idx_dip_gh + 266);

    auto tr_x_yyzz_yyyyz = pbuffer.data(idx_dip_gh + 268);

    auto tr_x_yyzz_yyyzz = pbuffer.data(idx_dip_gh + 269);

    auto tr_x_yyzz_yyzzz = pbuffer.data(idx_dip_gh + 270);

    auto tr_x_yyzz_yzzzz = pbuffer.data(idx_dip_gh + 271);

    auto tr_x_yyzz_zzzzz = pbuffer.data(idx_dip_gh + 272);

    auto tr_x_yzzz_xxxxz = pbuffer.data(idx_dip_gh + 275);

    auto tr_x_yzzz_xxxyz = pbuffer.data(idx_dip_gh + 277);

    auto tr_x_yzzz_xxxzz = pbuffer.data(idx_dip_gh + 278);

    auto tr_x_yzzz_xxyyz = pbuffer.data(idx_dip_gh + 280);

    auto tr_x_yzzz_xxyzz = pbuffer.data(idx_dip_gh + 281);

    auto tr_x_yzzz_xxzzz = pbuffer.data(idx_dip_gh + 282);

    auto tr_x_yzzz_xyyyz = pbuffer.data(idx_dip_gh + 284);

    auto tr_x_yzzz_xyyzz = pbuffer.data(idx_dip_gh + 285);

    auto tr_x_yzzz_xyzzz = pbuffer.data(idx_dip_gh + 286);

    auto tr_x_yzzz_xzzzz = pbuffer.data(idx_dip_gh + 287);

    auto tr_x_yzzz_yyyyz = pbuffer.data(idx_dip_gh + 289);

    auto tr_x_yzzz_yyyzz = pbuffer.data(idx_dip_gh + 290);

    auto tr_x_yzzz_yyzzz = pbuffer.data(idx_dip_gh + 291);

    auto tr_x_yzzz_yzzzz = pbuffer.data(idx_dip_gh + 292);

    auto tr_x_yzzz_zzzzz = pbuffer.data(idx_dip_gh + 293);

    auto tr_x_zzzz_xxxxx = pbuffer.data(idx_dip_gh + 294);

    auto tr_x_zzzz_xxxxy = pbuffer.data(idx_dip_gh + 295);

    auto tr_x_zzzz_xxxxz = pbuffer.data(idx_dip_gh + 296);

    auto tr_x_zzzz_xxxyy = pbuffer.data(idx_dip_gh + 297);

    auto tr_x_zzzz_xxxyz = pbuffer.data(idx_dip_gh + 298);

    auto tr_x_zzzz_xxxzz = pbuffer.data(idx_dip_gh + 299);

    auto tr_x_zzzz_xxyyy = pbuffer.data(idx_dip_gh + 300);

    auto tr_x_zzzz_xxyyz = pbuffer.data(idx_dip_gh + 301);

    auto tr_x_zzzz_xxyzz = pbuffer.data(idx_dip_gh + 302);

    auto tr_x_zzzz_xxzzz = pbuffer.data(idx_dip_gh + 303);

    auto tr_x_zzzz_xyyyy = pbuffer.data(idx_dip_gh + 304);

    auto tr_x_zzzz_xyyyz = pbuffer.data(idx_dip_gh + 305);

    auto tr_x_zzzz_xyyzz = pbuffer.data(idx_dip_gh + 306);

    auto tr_x_zzzz_xyzzz = pbuffer.data(idx_dip_gh + 307);

    auto tr_x_zzzz_xzzzz = pbuffer.data(idx_dip_gh + 308);

    auto tr_x_zzzz_yyyyy = pbuffer.data(idx_dip_gh + 309);

    auto tr_x_zzzz_yyyyz = pbuffer.data(idx_dip_gh + 310);

    auto tr_x_zzzz_yyyzz = pbuffer.data(idx_dip_gh + 311);

    auto tr_x_zzzz_yyzzz = pbuffer.data(idx_dip_gh + 312);

    auto tr_x_zzzz_yzzzz = pbuffer.data(idx_dip_gh + 313);

    auto tr_x_zzzz_zzzzz = pbuffer.data(idx_dip_gh + 314);

    auto tr_y_xxxx_xxxxx = pbuffer.data(idx_dip_gh + 315);

    auto tr_y_xxxx_xxxxy = pbuffer.data(idx_dip_gh + 316);

    auto tr_y_xxxx_xxxxz = pbuffer.data(idx_dip_gh + 317);

    auto tr_y_xxxx_xxxyy = pbuffer.data(idx_dip_gh + 318);

    auto tr_y_xxxx_xxxyz = pbuffer.data(idx_dip_gh + 319);

    auto tr_y_xxxx_xxxzz = pbuffer.data(idx_dip_gh + 320);

    auto tr_y_xxxx_xxyyy = pbuffer.data(idx_dip_gh + 321);

    auto tr_y_xxxx_xxyyz = pbuffer.data(idx_dip_gh + 322);

    auto tr_y_xxxx_xxyzz = pbuffer.data(idx_dip_gh + 323);

    auto tr_y_xxxx_xxzzz = pbuffer.data(idx_dip_gh + 324);

    auto tr_y_xxxx_xyyyy = pbuffer.data(idx_dip_gh + 325);

    auto tr_y_xxxx_xyyyz = pbuffer.data(idx_dip_gh + 326);

    auto tr_y_xxxx_xyyzz = pbuffer.data(idx_dip_gh + 327);

    auto tr_y_xxxx_xyzzz = pbuffer.data(idx_dip_gh + 328);

    auto tr_y_xxxx_xzzzz = pbuffer.data(idx_dip_gh + 329);

    auto tr_y_xxxx_yyyyy = pbuffer.data(idx_dip_gh + 330);

    auto tr_y_xxxx_yyyyz = pbuffer.data(idx_dip_gh + 331);

    auto tr_y_xxxx_yyyzz = pbuffer.data(idx_dip_gh + 332);

    auto tr_y_xxxx_yyzzz = pbuffer.data(idx_dip_gh + 333);

    auto tr_y_xxxx_yzzzz = pbuffer.data(idx_dip_gh + 334);

    auto tr_y_xxxx_zzzzz = pbuffer.data(idx_dip_gh + 335);

    auto tr_y_xxxy_xxxxy = pbuffer.data(idx_dip_gh + 337);

    auto tr_y_xxxy_xxxyy = pbuffer.data(idx_dip_gh + 339);

    auto tr_y_xxxy_xxxyz = pbuffer.data(idx_dip_gh + 340);

    auto tr_y_xxxy_xxyyy = pbuffer.data(idx_dip_gh + 342);

    auto tr_y_xxxy_xxyyz = pbuffer.data(idx_dip_gh + 343);

    auto tr_y_xxxy_xxyzz = pbuffer.data(idx_dip_gh + 344);

    auto tr_y_xxxy_xyyyy = pbuffer.data(idx_dip_gh + 346);

    auto tr_y_xxxy_xyyyz = pbuffer.data(idx_dip_gh + 347);

    auto tr_y_xxxy_xyyzz = pbuffer.data(idx_dip_gh + 348);

    auto tr_y_xxxy_xyzzz = pbuffer.data(idx_dip_gh + 349);

    auto tr_y_xxxy_yyyyy = pbuffer.data(idx_dip_gh + 351);

    auto tr_y_xxxy_yyyyz = pbuffer.data(idx_dip_gh + 352);

    auto tr_y_xxxy_yyyzz = pbuffer.data(idx_dip_gh + 353);

    auto tr_y_xxxy_yyzzz = pbuffer.data(idx_dip_gh + 354);

    auto tr_y_xxxy_yzzzz = pbuffer.data(idx_dip_gh + 355);

    auto tr_y_xxyy_xxxxx = pbuffer.data(idx_dip_gh + 378);

    auto tr_y_xxyy_xxxxy = pbuffer.data(idx_dip_gh + 379);

    auto tr_y_xxyy_xxxxz = pbuffer.data(idx_dip_gh + 380);

    auto tr_y_xxyy_xxxyy = pbuffer.data(idx_dip_gh + 381);

    auto tr_y_xxyy_xxxyz = pbuffer.data(idx_dip_gh + 382);

    auto tr_y_xxyy_xxxzz = pbuffer.data(idx_dip_gh + 383);

    auto tr_y_xxyy_xxyyy = pbuffer.data(idx_dip_gh + 384);

    auto tr_y_xxyy_xxyyz = pbuffer.data(idx_dip_gh + 385);

    auto tr_y_xxyy_xxyzz = pbuffer.data(idx_dip_gh + 386);

    auto tr_y_xxyy_xxzzz = pbuffer.data(idx_dip_gh + 387);

    auto tr_y_xxyy_xyyyy = pbuffer.data(idx_dip_gh + 388);

    auto tr_y_xxyy_xyyyz = pbuffer.data(idx_dip_gh + 389);

    auto tr_y_xxyy_xyyzz = pbuffer.data(idx_dip_gh + 390);

    auto tr_y_xxyy_xyzzz = pbuffer.data(idx_dip_gh + 391);

    auto tr_y_xxyy_xzzzz = pbuffer.data(idx_dip_gh + 392);

    auto tr_y_xxyy_yyyyy = pbuffer.data(idx_dip_gh + 393);

    auto tr_y_xxyy_yyyyz = pbuffer.data(idx_dip_gh + 394);

    auto tr_y_xxyy_yyyzz = pbuffer.data(idx_dip_gh + 395);

    auto tr_y_xxyy_yyzzz = pbuffer.data(idx_dip_gh + 396);

    auto tr_y_xxyy_yzzzz = pbuffer.data(idx_dip_gh + 397);

    auto tr_y_xxyy_zzzzz = pbuffer.data(idx_dip_gh + 398);

    auto tr_y_xxzz_xxxxz = pbuffer.data(idx_dip_gh + 422);

    auto tr_y_xxzz_xxxyz = pbuffer.data(idx_dip_gh + 424);

    auto tr_y_xxzz_xxxzz = pbuffer.data(idx_dip_gh + 425);

    auto tr_y_xxzz_xxyyz = pbuffer.data(idx_dip_gh + 427);

    auto tr_y_xxzz_xxyzz = pbuffer.data(idx_dip_gh + 428);

    auto tr_y_xxzz_xxzzz = pbuffer.data(idx_dip_gh + 429);

    auto tr_y_xxzz_xyyyz = pbuffer.data(idx_dip_gh + 431);

    auto tr_y_xxzz_xyyzz = pbuffer.data(idx_dip_gh + 432);

    auto tr_y_xxzz_xyzzz = pbuffer.data(idx_dip_gh + 433);

    auto tr_y_xxzz_xzzzz = pbuffer.data(idx_dip_gh + 434);

    auto tr_y_xxzz_yyyyz = pbuffer.data(idx_dip_gh + 436);

    auto tr_y_xxzz_yyyzz = pbuffer.data(idx_dip_gh + 437);

    auto tr_y_xxzz_yyzzz = pbuffer.data(idx_dip_gh + 438);

    auto tr_y_xxzz_yzzzz = pbuffer.data(idx_dip_gh + 439);

    auto tr_y_xxzz_zzzzz = pbuffer.data(idx_dip_gh + 440);

    auto tr_y_xyyy_xxxxx = pbuffer.data(idx_dip_gh + 441);

    auto tr_y_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 442);

    auto tr_y_xyyy_xxxxz = pbuffer.data(idx_dip_gh + 443);

    auto tr_y_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 444);

    auto tr_y_xyyy_xxxyz = pbuffer.data(idx_dip_gh + 445);

    auto tr_y_xyyy_xxxzz = pbuffer.data(idx_dip_gh + 446);

    auto tr_y_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 447);

    auto tr_y_xyyy_xxyyz = pbuffer.data(idx_dip_gh + 448);

    auto tr_y_xyyy_xxyzz = pbuffer.data(idx_dip_gh + 449);

    auto tr_y_xyyy_xxzzz = pbuffer.data(idx_dip_gh + 450);

    auto tr_y_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 451);

    auto tr_y_xyyy_xyyyz = pbuffer.data(idx_dip_gh + 452);

    auto tr_y_xyyy_xyyzz = pbuffer.data(idx_dip_gh + 453);

    auto tr_y_xyyy_xyzzz = pbuffer.data(idx_dip_gh + 454);

    auto tr_y_xyyy_xzzzz = pbuffer.data(idx_dip_gh + 455);

    auto tr_y_xyyy_yyyyy = pbuffer.data(idx_dip_gh + 456);

    auto tr_y_xyyy_yyyyz = pbuffer.data(idx_dip_gh + 457);

    auto tr_y_xyyy_yyyzz = pbuffer.data(idx_dip_gh + 458);

    auto tr_y_xyyy_yyzzz = pbuffer.data(idx_dip_gh + 459);

    auto tr_y_xyyy_yzzzz = pbuffer.data(idx_dip_gh + 460);

    auto tr_y_xyyy_zzzzz = pbuffer.data(idx_dip_gh + 461);

    auto tr_y_xyzz_xxxyz = pbuffer.data(idx_dip_gh + 487);

    auto tr_y_xyzz_xxyyz = pbuffer.data(idx_dip_gh + 490);

    auto tr_y_xyzz_xxyzz = pbuffer.data(idx_dip_gh + 491);

    auto tr_y_xyzz_xyyyz = pbuffer.data(idx_dip_gh + 494);

    auto tr_y_xyzz_xyyzz = pbuffer.data(idx_dip_gh + 495);

    auto tr_y_xyzz_xyzzz = pbuffer.data(idx_dip_gh + 496);

    auto tr_y_xyzz_yyyyz = pbuffer.data(idx_dip_gh + 499);

    auto tr_y_xyzz_yyyzz = pbuffer.data(idx_dip_gh + 500);

    auto tr_y_xyzz_yyzzz = pbuffer.data(idx_dip_gh + 501);

    auto tr_y_xyzz_yzzzz = pbuffer.data(idx_dip_gh + 502);

    auto tr_y_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 506);

    auto tr_y_xzzz_xxxyz = pbuffer.data(idx_dip_gh + 508);

    auto tr_y_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 509);

    auto tr_y_xzzz_xxyyz = pbuffer.data(idx_dip_gh + 511);

    auto tr_y_xzzz_xxyzz = pbuffer.data(idx_dip_gh + 512);

    auto tr_y_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 513);

    auto tr_y_xzzz_xyyyz = pbuffer.data(idx_dip_gh + 515);

    auto tr_y_xzzz_xyyzz = pbuffer.data(idx_dip_gh + 516);

    auto tr_y_xzzz_xyzzz = pbuffer.data(idx_dip_gh + 517);

    auto tr_y_xzzz_xzzzz = pbuffer.data(idx_dip_gh + 518);

    auto tr_y_xzzz_yyyyz = pbuffer.data(idx_dip_gh + 520);

    auto tr_y_xzzz_yyyzz = pbuffer.data(idx_dip_gh + 521);

    auto tr_y_xzzz_yyzzz = pbuffer.data(idx_dip_gh + 522);

    auto tr_y_xzzz_yzzzz = pbuffer.data(idx_dip_gh + 523);

    auto tr_y_xzzz_zzzzz = pbuffer.data(idx_dip_gh + 524);

    auto tr_y_yyyy_xxxxx = pbuffer.data(idx_dip_gh + 525);

    auto tr_y_yyyy_xxxxy = pbuffer.data(idx_dip_gh + 526);

    auto tr_y_yyyy_xxxxz = pbuffer.data(idx_dip_gh + 527);

    auto tr_y_yyyy_xxxyy = pbuffer.data(idx_dip_gh + 528);

    auto tr_y_yyyy_xxxyz = pbuffer.data(idx_dip_gh + 529);

    auto tr_y_yyyy_xxxzz = pbuffer.data(idx_dip_gh + 530);

    auto tr_y_yyyy_xxyyy = pbuffer.data(idx_dip_gh + 531);

    auto tr_y_yyyy_xxyyz = pbuffer.data(idx_dip_gh + 532);

    auto tr_y_yyyy_xxyzz = pbuffer.data(idx_dip_gh + 533);

    auto tr_y_yyyy_xxzzz = pbuffer.data(idx_dip_gh + 534);

    auto tr_y_yyyy_xyyyy = pbuffer.data(idx_dip_gh + 535);

    auto tr_y_yyyy_xyyyz = pbuffer.data(idx_dip_gh + 536);

    auto tr_y_yyyy_xyyzz = pbuffer.data(idx_dip_gh + 537);

    auto tr_y_yyyy_xyzzz = pbuffer.data(idx_dip_gh + 538);

    auto tr_y_yyyy_xzzzz = pbuffer.data(idx_dip_gh + 539);

    auto tr_y_yyyy_yyyyy = pbuffer.data(idx_dip_gh + 540);

    auto tr_y_yyyy_yyyyz = pbuffer.data(idx_dip_gh + 541);

    auto tr_y_yyyy_yyyzz = pbuffer.data(idx_dip_gh + 542);

    auto tr_y_yyyy_yyzzz = pbuffer.data(idx_dip_gh + 543);

    auto tr_y_yyyy_yzzzz = pbuffer.data(idx_dip_gh + 544);

    auto tr_y_yyyy_zzzzz = pbuffer.data(idx_dip_gh + 545);

    auto tr_y_yyyz_xxxxy = pbuffer.data(idx_dip_gh + 547);

    auto tr_y_yyyz_xxxxz = pbuffer.data(idx_dip_gh + 548);

    auto tr_y_yyyz_xxxyy = pbuffer.data(idx_dip_gh + 549);

    auto tr_y_yyyz_xxxyz = pbuffer.data(idx_dip_gh + 550);

    auto tr_y_yyyz_xxxzz = pbuffer.data(idx_dip_gh + 551);

    auto tr_y_yyyz_xxyyy = pbuffer.data(idx_dip_gh + 552);

    auto tr_y_yyyz_xxyyz = pbuffer.data(idx_dip_gh + 553);

    auto tr_y_yyyz_xxyzz = pbuffer.data(idx_dip_gh + 554);

    auto tr_y_yyyz_xxzzz = pbuffer.data(idx_dip_gh + 555);

    auto tr_y_yyyz_xyyyy = pbuffer.data(idx_dip_gh + 556);

    auto tr_y_yyyz_xyyyz = pbuffer.data(idx_dip_gh + 557);

    auto tr_y_yyyz_xyyzz = pbuffer.data(idx_dip_gh + 558);

    auto tr_y_yyyz_xyzzz = pbuffer.data(idx_dip_gh + 559);

    auto tr_y_yyyz_xzzzz = pbuffer.data(idx_dip_gh + 560);

    auto tr_y_yyyz_yyyyy = pbuffer.data(idx_dip_gh + 561);

    auto tr_y_yyyz_yyyyz = pbuffer.data(idx_dip_gh + 562);

    auto tr_y_yyyz_yyyzz = pbuffer.data(idx_dip_gh + 563);

    auto tr_y_yyyz_yyzzz = pbuffer.data(idx_dip_gh + 564);

    auto tr_y_yyyz_yzzzz = pbuffer.data(idx_dip_gh + 565);

    auto tr_y_yyyz_zzzzz = pbuffer.data(idx_dip_gh + 566);

    auto tr_y_yyzz_xxxxx = pbuffer.data(idx_dip_gh + 567);

    auto tr_y_yyzz_xxxxy = pbuffer.data(idx_dip_gh + 568);

    auto tr_y_yyzz_xxxxz = pbuffer.data(idx_dip_gh + 569);

    auto tr_y_yyzz_xxxyy = pbuffer.data(idx_dip_gh + 570);

    auto tr_y_yyzz_xxxyz = pbuffer.data(idx_dip_gh + 571);

    auto tr_y_yyzz_xxxzz = pbuffer.data(idx_dip_gh + 572);

    auto tr_y_yyzz_xxyyy = pbuffer.data(idx_dip_gh + 573);

    auto tr_y_yyzz_xxyyz = pbuffer.data(idx_dip_gh + 574);

    auto tr_y_yyzz_xxyzz = pbuffer.data(idx_dip_gh + 575);

    auto tr_y_yyzz_xxzzz = pbuffer.data(idx_dip_gh + 576);

    auto tr_y_yyzz_xyyyy = pbuffer.data(idx_dip_gh + 577);

    auto tr_y_yyzz_xyyyz = pbuffer.data(idx_dip_gh + 578);

    auto tr_y_yyzz_xyyzz = pbuffer.data(idx_dip_gh + 579);

    auto tr_y_yyzz_xyzzz = pbuffer.data(idx_dip_gh + 580);

    auto tr_y_yyzz_xzzzz = pbuffer.data(idx_dip_gh + 581);

    auto tr_y_yyzz_yyyyy = pbuffer.data(idx_dip_gh + 582);

    auto tr_y_yyzz_yyyyz = pbuffer.data(idx_dip_gh + 583);

    auto tr_y_yyzz_yyyzz = pbuffer.data(idx_dip_gh + 584);

    auto tr_y_yyzz_yyzzz = pbuffer.data(idx_dip_gh + 585);

    auto tr_y_yyzz_yzzzz = pbuffer.data(idx_dip_gh + 586);

    auto tr_y_yyzz_zzzzz = pbuffer.data(idx_dip_gh + 587);

    auto tr_y_yzzz_xxxxx = pbuffer.data(idx_dip_gh + 588);

    auto tr_y_yzzz_xxxxy = pbuffer.data(idx_dip_gh + 589);

    auto tr_y_yzzz_xxxxz = pbuffer.data(idx_dip_gh + 590);

    auto tr_y_yzzz_xxxyy = pbuffer.data(idx_dip_gh + 591);

    auto tr_y_yzzz_xxxyz = pbuffer.data(idx_dip_gh + 592);

    auto tr_y_yzzz_xxxzz = pbuffer.data(idx_dip_gh + 593);

    auto tr_y_yzzz_xxyyy = pbuffer.data(idx_dip_gh + 594);

    auto tr_y_yzzz_xxyyz = pbuffer.data(idx_dip_gh + 595);

    auto tr_y_yzzz_xxyzz = pbuffer.data(idx_dip_gh + 596);

    auto tr_y_yzzz_xxzzz = pbuffer.data(idx_dip_gh + 597);

    auto tr_y_yzzz_xyyyy = pbuffer.data(idx_dip_gh + 598);

    auto tr_y_yzzz_xyyyz = pbuffer.data(idx_dip_gh + 599);

    auto tr_y_yzzz_xyyzz = pbuffer.data(idx_dip_gh + 600);

    auto tr_y_yzzz_xyzzz = pbuffer.data(idx_dip_gh + 601);

    auto tr_y_yzzz_xzzzz = pbuffer.data(idx_dip_gh + 602);

    auto tr_y_yzzz_yyyyy = pbuffer.data(idx_dip_gh + 603);

    auto tr_y_yzzz_yyyyz = pbuffer.data(idx_dip_gh + 604);

    auto tr_y_yzzz_yyyzz = pbuffer.data(idx_dip_gh + 605);

    auto tr_y_yzzz_yyzzz = pbuffer.data(idx_dip_gh + 606);

    auto tr_y_yzzz_yzzzz = pbuffer.data(idx_dip_gh + 607);

    auto tr_y_yzzz_zzzzz = pbuffer.data(idx_dip_gh + 608);

    auto tr_y_zzzz_xxxxx = pbuffer.data(idx_dip_gh + 609);

    auto tr_y_zzzz_xxxxy = pbuffer.data(idx_dip_gh + 610);

    auto tr_y_zzzz_xxxxz = pbuffer.data(idx_dip_gh + 611);

    auto tr_y_zzzz_xxxyy = pbuffer.data(idx_dip_gh + 612);

    auto tr_y_zzzz_xxxyz = pbuffer.data(idx_dip_gh + 613);

    auto tr_y_zzzz_xxxzz = pbuffer.data(idx_dip_gh + 614);

    auto tr_y_zzzz_xxyyy = pbuffer.data(idx_dip_gh + 615);

    auto tr_y_zzzz_xxyyz = pbuffer.data(idx_dip_gh + 616);

    auto tr_y_zzzz_xxyzz = pbuffer.data(idx_dip_gh + 617);

    auto tr_y_zzzz_xxzzz = pbuffer.data(idx_dip_gh + 618);

    auto tr_y_zzzz_xyyyy = pbuffer.data(idx_dip_gh + 619);

    auto tr_y_zzzz_xyyyz = pbuffer.data(idx_dip_gh + 620);

    auto tr_y_zzzz_xyyzz = pbuffer.data(idx_dip_gh + 621);

    auto tr_y_zzzz_xyzzz = pbuffer.data(idx_dip_gh + 622);

    auto tr_y_zzzz_xzzzz = pbuffer.data(idx_dip_gh + 623);

    auto tr_y_zzzz_yyyyy = pbuffer.data(idx_dip_gh + 624);

    auto tr_y_zzzz_yyyyz = pbuffer.data(idx_dip_gh + 625);

    auto tr_y_zzzz_yyyzz = pbuffer.data(idx_dip_gh + 626);

    auto tr_y_zzzz_yyzzz = pbuffer.data(idx_dip_gh + 627);

    auto tr_y_zzzz_yzzzz = pbuffer.data(idx_dip_gh + 628);

    auto tr_y_zzzz_zzzzz = pbuffer.data(idx_dip_gh + 629);

    auto tr_z_xxxx_xxxxx = pbuffer.data(idx_dip_gh + 630);

    auto tr_z_xxxx_xxxxy = pbuffer.data(idx_dip_gh + 631);

    auto tr_z_xxxx_xxxxz = pbuffer.data(idx_dip_gh + 632);

    auto tr_z_xxxx_xxxyy = pbuffer.data(idx_dip_gh + 633);

    auto tr_z_xxxx_xxxyz = pbuffer.data(idx_dip_gh + 634);

    auto tr_z_xxxx_xxxzz = pbuffer.data(idx_dip_gh + 635);

    auto tr_z_xxxx_xxyyy = pbuffer.data(idx_dip_gh + 636);

    auto tr_z_xxxx_xxyyz = pbuffer.data(idx_dip_gh + 637);

    auto tr_z_xxxx_xxyzz = pbuffer.data(idx_dip_gh + 638);

    auto tr_z_xxxx_xxzzz = pbuffer.data(idx_dip_gh + 639);

    auto tr_z_xxxx_xyyyy = pbuffer.data(idx_dip_gh + 640);

    auto tr_z_xxxx_xyyyz = pbuffer.data(idx_dip_gh + 641);

    auto tr_z_xxxx_xyyzz = pbuffer.data(idx_dip_gh + 642);

    auto tr_z_xxxx_xyzzz = pbuffer.data(idx_dip_gh + 643);

    auto tr_z_xxxx_xzzzz = pbuffer.data(idx_dip_gh + 644);

    auto tr_z_xxxx_yyyyy = pbuffer.data(idx_dip_gh + 645);

    auto tr_z_xxxx_yyyyz = pbuffer.data(idx_dip_gh + 646);

    auto tr_z_xxxx_yyyzz = pbuffer.data(idx_dip_gh + 647);

    auto tr_z_xxxx_yyzzz = pbuffer.data(idx_dip_gh + 648);

    auto tr_z_xxxx_yzzzz = pbuffer.data(idx_dip_gh + 649);

    auto tr_z_xxxx_zzzzz = pbuffer.data(idx_dip_gh + 650);

    auto tr_z_xxxz_xxxxx = pbuffer.data(idx_dip_gh + 672);

    auto tr_z_xxxz_xxxxy = pbuffer.data(idx_dip_gh + 673);

    auto tr_z_xxxz_xxxxz = pbuffer.data(idx_dip_gh + 674);

    auto tr_z_xxxz_xxxyy = pbuffer.data(idx_dip_gh + 675);

    auto tr_z_xxxz_xxxyz = pbuffer.data(idx_dip_gh + 676);

    auto tr_z_xxxz_xxxzz = pbuffer.data(idx_dip_gh + 677);

    auto tr_z_xxxz_xxyyy = pbuffer.data(idx_dip_gh + 678);

    auto tr_z_xxxz_xxyyz = pbuffer.data(idx_dip_gh + 679);

    auto tr_z_xxxz_xxyzz = pbuffer.data(idx_dip_gh + 680);

    auto tr_z_xxxz_xxzzz = pbuffer.data(idx_dip_gh + 681);

    auto tr_z_xxxz_xyyyy = pbuffer.data(idx_dip_gh + 682);

    auto tr_z_xxxz_xyyyz = pbuffer.data(idx_dip_gh + 683);

    auto tr_z_xxxz_xyyzz = pbuffer.data(idx_dip_gh + 684);

    auto tr_z_xxxz_xyzzz = pbuffer.data(idx_dip_gh + 685);

    auto tr_z_xxxz_xzzzz = pbuffer.data(idx_dip_gh + 686);

    auto tr_z_xxxz_yyyyz = pbuffer.data(idx_dip_gh + 688);

    auto tr_z_xxxz_yyyzz = pbuffer.data(idx_dip_gh + 689);

    auto tr_z_xxxz_yyzzz = pbuffer.data(idx_dip_gh + 690);

    auto tr_z_xxxz_yzzzz = pbuffer.data(idx_dip_gh + 691);

    auto tr_z_xxxz_zzzzz = pbuffer.data(idx_dip_gh + 692);

    auto tr_z_xxyy_xxxxy = pbuffer.data(idx_dip_gh + 694);

    auto tr_z_xxyy_xxxyy = pbuffer.data(idx_dip_gh + 696);

    auto tr_z_xxyy_xxxyz = pbuffer.data(idx_dip_gh + 697);

    auto tr_z_xxyy_xxyyy = pbuffer.data(idx_dip_gh + 699);

    auto tr_z_xxyy_xxyyz = pbuffer.data(idx_dip_gh + 700);

    auto tr_z_xxyy_xxyzz = pbuffer.data(idx_dip_gh + 701);

    auto tr_z_xxyy_xyyyy = pbuffer.data(idx_dip_gh + 703);

    auto tr_z_xxyy_xyyyz = pbuffer.data(idx_dip_gh + 704);

    auto tr_z_xxyy_xyyzz = pbuffer.data(idx_dip_gh + 705);

    auto tr_z_xxyy_xyzzz = pbuffer.data(idx_dip_gh + 706);

    auto tr_z_xxyy_yyyyy = pbuffer.data(idx_dip_gh + 708);

    auto tr_z_xxyy_yyyyz = pbuffer.data(idx_dip_gh + 709);

    auto tr_z_xxyy_yyyzz = pbuffer.data(idx_dip_gh + 710);

    auto tr_z_xxyy_yyzzz = pbuffer.data(idx_dip_gh + 711);

    auto tr_z_xxyy_yzzzz = pbuffer.data(idx_dip_gh + 712);

    auto tr_z_xxzz_xxxxx = pbuffer.data(idx_dip_gh + 735);

    auto tr_z_xxzz_xxxxy = pbuffer.data(idx_dip_gh + 736);

    auto tr_z_xxzz_xxxxz = pbuffer.data(idx_dip_gh + 737);

    auto tr_z_xxzz_xxxyy = pbuffer.data(idx_dip_gh + 738);

    auto tr_z_xxzz_xxxyz = pbuffer.data(idx_dip_gh + 739);

    auto tr_z_xxzz_xxxzz = pbuffer.data(idx_dip_gh + 740);

    auto tr_z_xxzz_xxyyy = pbuffer.data(idx_dip_gh + 741);

    auto tr_z_xxzz_xxyyz = pbuffer.data(idx_dip_gh + 742);

    auto tr_z_xxzz_xxyzz = pbuffer.data(idx_dip_gh + 743);

    auto tr_z_xxzz_xxzzz = pbuffer.data(idx_dip_gh + 744);

    auto tr_z_xxzz_xyyyy = pbuffer.data(idx_dip_gh + 745);

    auto tr_z_xxzz_xyyyz = pbuffer.data(idx_dip_gh + 746);

    auto tr_z_xxzz_xyyzz = pbuffer.data(idx_dip_gh + 747);

    auto tr_z_xxzz_xyzzz = pbuffer.data(idx_dip_gh + 748);

    auto tr_z_xxzz_xzzzz = pbuffer.data(idx_dip_gh + 749);

    auto tr_z_xxzz_yyyyy = pbuffer.data(idx_dip_gh + 750);

    auto tr_z_xxzz_yyyyz = pbuffer.data(idx_dip_gh + 751);

    auto tr_z_xxzz_yyyzz = pbuffer.data(idx_dip_gh + 752);

    auto tr_z_xxzz_yyzzz = pbuffer.data(idx_dip_gh + 753);

    auto tr_z_xxzz_yzzzz = pbuffer.data(idx_dip_gh + 754);

    auto tr_z_xxzz_zzzzz = pbuffer.data(idx_dip_gh + 755);

    auto tr_z_xyyy_xxxxy = pbuffer.data(idx_dip_gh + 757);

    auto tr_z_xyyy_xxxyy = pbuffer.data(idx_dip_gh + 759);

    auto tr_z_xyyy_xxxyz = pbuffer.data(idx_dip_gh + 760);

    auto tr_z_xyyy_xxyyy = pbuffer.data(idx_dip_gh + 762);

    auto tr_z_xyyy_xxyyz = pbuffer.data(idx_dip_gh + 763);

    auto tr_z_xyyy_xxyzz = pbuffer.data(idx_dip_gh + 764);

    auto tr_z_xyyy_xyyyy = pbuffer.data(idx_dip_gh + 766);

    auto tr_z_xyyy_xyyyz = pbuffer.data(idx_dip_gh + 767);

    auto tr_z_xyyy_xyyzz = pbuffer.data(idx_dip_gh + 768);

    auto tr_z_xyyy_xyzzz = pbuffer.data(idx_dip_gh + 769);

    auto tr_z_xyyy_yyyyy = pbuffer.data(idx_dip_gh + 771);

    auto tr_z_xyyy_yyyyz = pbuffer.data(idx_dip_gh + 772);

    auto tr_z_xyyy_yyyzz = pbuffer.data(idx_dip_gh + 773);

    auto tr_z_xyyy_yyzzz = pbuffer.data(idx_dip_gh + 774);

    auto tr_z_xyyy_yzzzz = pbuffer.data(idx_dip_gh + 775);

    auto tr_z_xyyz_xxxyz = pbuffer.data(idx_dip_gh + 781);

    auto tr_z_xyyz_xxyyz = pbuffer.data(idx_dip_gh + 784);

    auto tr_z_xyyz_xxyzz = pbuffer.data(idx_dip_gh + 785);

    auto tr_z_xyyz_xyyyz = pbuffer.data(idx_dip_gh + 788);

    auto tr_z_xyyz_xyyzz = pbuffer.data(idx_dip_gh + 789);

    auto tr_z_xyyz_xyzzz = pbuffer.data(idx_dip_gh + 790);

    auto tr_z_xyyz_yyyyz = pbuffer.data(idx_dip_gh + 793);

    auto tr_z_xyyz_yyyzz = pbuffer.data(idx_dip_gh + 794);

    auto tr_z_xyyz_yyzzz = pbuffer.data(idx_dip_gh + 795);

    auto tr_z_xyyz_yzzzz = pbuffer.data(idx_dip_gh + 796);

    auto tr_z_xzzz_xxxxx = pbuffer.data(idx_dip_gh + 819);

    auto tr_z_xzzz_xxxxy = pbuffer.data(idx_dip_gh + 820);

    auto tr_z_xzzz_xxxxz = pbuffer.data(idx_dip_gh + 821);

    auto tr_z_xzzz_xxxyy = pbuffer.data(idx_dip_gh + 822);

    auto tr_z_xzzz_xxxyz = pbuffer.data(idx_dip_gh + 823);

    auto tr_z_xzzz_xxxzz = pbuffer.data(idx_dip_gh + 824);

    auto tr_z_xzzz_xxyyy = pbuffer.data(idx_dip_gh + 825);

    auto tr_z_xzzz_xxyyz = pbuffer.data(idx_dip_gh + 826);

    auto tr_z_xzzz_xxyzz = pbuffer.data(idx_dip_gh + 827);

    auto tr_z_xzzz_xxzzz = pbuffer.data(idx_dip_gh + 828);

    auto tr_z_xzzz_xyyyy = pbuffer.data(idx_dip_gh + 829);

    auto tr_z_xzzz_xyyyz = pbuffer.data(idx_dip_gh + 830);

    auto tr_z_xzzz_xyyzz = pbuffer.data(idx_dip_gh + 831);

    auto tr_z_xzzz_xyzzz = pbuffer.data(idx_dip_gh + 832);

    auto tr_z_xzzz_xzzzz = pbuffer.data(idx_dip_gh + 833);

    auto tr_z_xzzz_yyyyy = pbuffer.data(idx_dip_gh + 834);

    auto tr_z_xzzz_yyyyz = pbuffer.data(idx_dip_gh + 835);

    auto tr_z_xzzz_yyyzz = pbuffer.data(idx_dip_gh + 836);

    auto tr_z_xzzz_yyzzz = pbuffer.data(idx_dip_gh + 837);

    auto tr_z_xzzz_yzzzz = pbuffer.data(idx_dip_gh + 838);

    auto tr_z_xzzz_zzzzz = pbuffer.data(idx_dip_gh + 839);

    auto tr_z_yyyy_xxxxx = pbuffer.data(idx_dip_gh + 840);

    auto tr_z_yyyy_xxxxy = pbuffer.data(idx_dip_gh + 841);

    auto tr_z_yyyy_xxxxz = pbuffer.data(idx_dip_gh + 842);

    auto tr_z_yyyy_xxxyy = pbuffer.data(idx_dip_gh + 843);

    auto tr_z_yyyy_xxxyz = pbuffer.data(idx_dip_gh + 844);

    auto tr_z_yyyy_xxxzz = pbuffer.data(idx_dip_gh + 845);

    auto tr_z_yyyy_xxyyy = pbuffer.data(idx_dip_gh + 846);

    auto tr_z_yyyy_xxyyz = pbuffer.data(idx_dip_gh + 847);

    auto tr_z_yyyy_xxyzz = pbuffer.data(idx_dip_gh + 848);

    auto tr_z_yyyy_xxzzz = pbuffer.data(idx_dip_gh + 849);

    auto tr_z_yyyy_xyyyy = pbuffer.data(idx_dip_gh + 850);

    auto tr_z_yyyy_xyyyz = pbuffer.data(idx_dip_gh + 851);

    auto tr_z_yyyy_xyyzz = pbuffer.data(idx_dip_gh + 852);

    auto tr_z_yyyy_xyzzz = pbuffer.data(idx_dip_gh + 853);

    auto tr_z_yyyy_xzzzz = pbuffer.data(idx_dip_gh + 854);

    auto tr_z_yyyy_yyyyy = pbuffer.data(idx_dip_gh + 855);

    auto tr_z_yyyy_yyyyz = pbuffer.data(idx_dip_gh + 856);

    auto tr_z_yyyy_yyyzz = pbuffer.data(idx_dip_gh + 857);

    auto tr_z_yyyy_yyzzz = pbuffer.data(idx_dip_gh + 858);

    auto tr_z_yyyy_yzzzz = pbuffer.data(idx_dip_gh + 859);

    auto tr_z_yyyy_zzzzz = pbuffer.data(idx_dip_gh + 860);

    auto tr_z_yyyz_xxxxx = pbuffer.data(idx_dip_gh + 861);

    auto tr_z_yyyz_xxxxy = pbuffer.data(idx_dip_gh + 862);

    auto tr_z_yyyz_xxxxz = pbuffer.data(idx_dip_gh + 863);

    auto tr_z_yyyz_xxxyy = pbuffer.data(idx_dip_gh + 864);

    auto tr_z_yyyz_xxxyz = pbuffer.data(idx_dip_gh + 865);

    auto tr_z_yyyz_xxxzz = pbuffer.data(idx_dip_gh + 866);

    auto tr_z_yyyz_xxyyy = pbuffer.data(idx_dip_gh + 867);

    auto tr_z_yyyz_xxyyz = pbuffer.data(idx_dip_gh + 868);

    auto tr_z_yyyz_xxyzz = pbuffer.data(idx_dip_gh + 869);

    auto tr_z_yyyz_xxzzz = pbuffer.data(idx_dip_gh + 870);

    auto tr_z_yyyz_xyyyy = pbuffer.data(idx_dip_gh + 871);

    auto tr_z_yyyz_xyyyz = pbuffer.data(idx_dip_gh + 872);

    auto tr_z_yyyz_xyyzz = pbuffer.data(idx_dip_gh + 873);

    auto tr_z_yyyz_xyzzz = pbuffer.data(idx_dip_gh + 874);

    auto tr_z_yyyz_xzzzz = pbuffer.data(idx_dip_gh + 875);

    auto tr_z_yyyz_yyyyy = pbuffer.data(idx_dip_gh + 876);

    auto tr_z_yyyz_yyyyz = pbuffer.data(idx_dip_gh + 877);

    auto tr_z_yyyz_yyyzz = pbuffer.data(idx_dip_gh + 878);

    auto tr_z_yyyz_yyzzz = pbuffer.data(idx_dip_gh + 879);

    auto tr_z_yyyz_yzzzz = pbuffer.data(idx_dip_gh + 880);

    auto tr_z_yyyz_zzzzz = pbuffer.data(idx_dip_gh + 881);

    auto tr_z_yyzz_xxxxx = pbuffer.data(idx_dip_gh + 882);

    auto tr_z_yyzz_xxxxy = pbuffer.data(idx_dip_gh + 883);

    auto tr_z_yyzz_xxxxz = pbuffer.data(idx_dip_gh + 884);

    auto tr_z_yyzz_xxxyy = pbuffer.data(idx_dip_gh + 885);

    auto tr_z_yyzz_xxxyz = pbuffer.data(idx_dip_gh + 886);

    auto tr_z_yyzz_xxxzz = pbuffer.data(idx_dip_gh + 887);

    auto tr_z_yyzz_xxyyy = pbuffer.data(idx_dip_gh + 888);

    auto tr_z_yyzz_xxyyz = pbuffer.data(idx_dip_gh + 889);

    auto tr_z_yyzz_xxyzz = pbuffer.data(idx_dip_gh + 890);

    auto tr_z_yyzz_xxzzz = pbuffer.data(idx_dip_gh + 891);

    auto tr_z_yyzz_xyyyy = pbuffer.data(idx_dip_gh + 892);

    auto tr_z_yyzz_xyyyz = pbuffer.data(idx_dip_gh + 893);

    auto tr_z_yyzz_xyyzz = pbuffer.data(idx_dip_gh + 894);

    auto tr_z_yyzz_xyzzz = pbuffer.data(idx_dip_gh + 895);

    auto tr_z_yyzz_xzzzz = pbuffer.data(idx_dip_gh + 896);

    auto tr_z_yyzz_yyyyy = pbuffer.data(idx_dip_gh + 897);

    auto tr_z_yyzz_yyyyz = pbuffer.data(idx_dip_gh + 898);

    auto tr_z_yyzz_yyyzz = pbuffer.data(idx_dip_gh + 899);

    auto tr_z_yyzz_yyzzz = pbuffer.data(idx_dip_gh + 900);

    auto tr_z_yyzz_yzzzz = pbuffer.data(idx_dip_gh + 901);

    auto tr_z_yyzz_zzzzz = pbuffer.data(idx_dip_gh + 902);

    auto tr_z_yzzz_xxxxx = pbuffer.data(idx_dip_gh + 903);

    auto tr_z_yzzz_xxxxy = pbuffer.data(idx_dip_gh + 904);

    auto tr_z_yzzz_xxxxz = pbuffer.data(idx_dip_gh + 905);

    auto tr_z_yzzz_xxxyy = pbuffer.data(idx_dip_gh + 906);

    auto tr_z_yzzz_xxxyz = pbuffer.data(idx_dip_gh + 907);

    auto tr_z_yzzz_xxxzz = pbuffer.data(idx_dip_gh + 908);

    auto tr_z_yzzz_xxyyy = pbuffer.data(idx_dip_gh + 909);

    auto tr_z_yzzz_xxyyz = pbuffer.data(idx_dip_gh + 910);

    auto tr_z_yzzz_xxyzz = pbuffer.data(idx_dip_gh + 911);

    auto tr_z_yzzz_xxzzz = pbuffer.data(idx_dip_gh + 912);

    auto tr_z_yzzz_xyyyy = pbuffer.data(idx_dip_gh + 913);

    auto tr_z_yzzz_xyyyz = pbuffer.data(idx_dip_gh + 914);

    auto tr_z_yzzz_xyyzz = pbuffer.data(idx_dip_gh + 915);

    auto tr_z_yzzz_xyzzz = pbuffer.data(idx_dip_gh + 916);

    auto tr_z_yzzz_xzzzz = pbuffer.data(idx_dip_gh + 917);

    auto tr_z_yzzz_yyyyy = pbuffer.data(idx_dip_gh + 918);

    auto tr_z_yzzz_yyyyz = pbuffer.data(idx_dip_gh + 919);

    auto tr_z_yzzz_yyyzz = pbuffer.data(idx_dip_gh + 920);

    auto tr_z_yzzz_yyzzz = pbuffer.data(idx_dip_gh + 921);

    auto tr_z_yzzz_yzzzz = pbuffer.data(idx_dip_gh + 922);

    auto tr_z_yzzz_zzzzz = pbuffer.data(idx_dip_gh + 923);

    auto tr_z_zzzz_xxxxx = pbuffer.data(idx_dip_gh + 924);

    auto tr_z_zzzz_xxxxy = pbuffer.data(idx_dip_gh + 925);

    auto tr_z_zzzz_xxxxz = pbuffer.data(idx_dip_gh + 926);

    auto tr_z_zzzz_xxxyy = pbuffer.data(idx_dip_gh + 927);

    auto tr_z_zzzz_xxxyz = pbuffer.data(idx_dip_gh + 928);

    auto tr_z_zzzz_xxxzz = pbuffer.data(idx_dip_gh + 929);

    auto tr_z_zzzz_xxyyy = pbuffer.data(idx_dip_gh + 930);

    auto tr_z_zzzz_xxyyz = pbuffer.data(idx_dip_gh + 931);

    auto tr_z_zzzz_xxyzz = pbuffer.data(idx_dip_gh + 932);

    auto tr_z_zzzz_xxzzz = pbuffer.data(idx_dip_gh + 933);

    auto tr_z_zzzz_xyyyy = pbuffer.data(idx_dip_gh + 934);

    auto tr_z_zzzz_xyyyz = pbuffer.data(idx_dip_gh + 935);

    auto tr_z_zzzz_xyyzz = pbuffer.data(idx_dip_gh + 936);

    auto tr_z_zzzz_xyzzz = pbuffer.data(idx_dip_gh + 937);

    auto tr_z_zzzz_xzzzz = pbuffer.data(idx_dip_gh + 938);

    auto tr_z_zzzz_yyyyy = pbuffer.data(idx_dip_gh + 939);

    auto tr_z_zzzz_yyyyz = pbuffer.data(idx_dip_gh + 940);

    auto tr_z_zzzz_yyyzz = pbuffer.data(idx_dip_gh + 941);

    auto tr_z_zzzz_yyzzz = pbuffer.data(idx_dip_gh + 942);

    auto tr_z_zzzz_yzzzz = pbuffer.data(idx_dip_gh + 943);

    auto tr_z_zzzz_zzzzz = pbuffer.data(idx_dip_gh + 944);

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

    auto ts_xxxz_xxxxxz = pbuffer.data(idx_ovl_gi + 58);

    auto ts_xxxz_xxxxzz = pbuffer.data(idx_ovl_gi + 61);

    auto ts_xxxz_xxxzzz = pbuffer.data(idx_ovl_gi + 65);

    auto ts_xxxz_xxzzzz = pbuffer.data(idx_ovl_gi + 70);

    auto ts_xxxz_xzzzzz = pbuffer.data(idx_ovl_gi + 76);

    auto ts_xxyy_xxxxxy = pbuffer.data(idx_ovl_gi + 85);

    auto ts_xxyy_xxxxyy = pbuffer.data(idx_ovl_gi + 87);

    auto ts_xxyy_xxxyyy = pbuffer.data(idx_ovl_gi + 90);

    auto ts_xxyy_xxyyyy = pbuffer.data(idx_ovl_gi + 94);

    auto ts_xxyy_xyyyyy = pbuffer.data(idx_ovl_gi + 99);

    auto ts_xxyy_yyyyyy = pbuffer.data(idx_ovl_gi + 105);

    auto ts_xxyy_yyyyyz = pbuffer.data(idx_ovl_gi + 106);

    auto ts_xxyy_yyyyzz = pbuffer.data(idx_ovl_gi + 107);

    auto ts_xxyy_yyyzzz = pbuffer.data(idx_ovl_gi + 108);

    auto ts_xxyy_yyzzzz = pbuffer.data(idx_ovl_gi + 109);

    auto ts_xxyy_yzzzzz = pbuffer.data(idx_ovl_gi + 110);

    auto ts_xxzz_xxxxxx = pbuffer.data(idx_ovl_gi + 140);

    auto ts_xxzz_xxxxxz = pbuffer.data(idx_ovl_gi + 142);

    auto ts_xxzz_xxxxzz = pbuffer.data(idx_ovl_gi + 145);

    auto ts_xxzz_xxxzzz = pbuffer.data(idx_ovl_gi + 149);

    auto ts_xxzz_xxzzzz = pbuffer.data(idx_ovl_gi + 154);

    auto ts_xxzz_xzzzzz = pbuffer.data(idx_ovl_gi + 160);

    auto ts_xxzz_yyyyyz = pbuffer.data(idx_ovl_gi + 162);

    auto ts_xxzz_yyyyzz = pbuffer.data(idx_ovl_gi + 163);

    auto ts_xxzz_yyyzzz = pbuffer.data(idx_ovl_gi + 164);

    auto ts_xxzz_yyzzzz = pbuffer.data(idx_ovl_gi + 165);

    auto ts_xxzz_yzzzzz = pbuffer.data(idx_ovl_gi + 166);

    auto ts_xxzz_zzzzzz = pbuffer.data(idx_ovl_gi + 167);

    auto ts_xyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 189);

    auto ts_xyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 190);

    auto ts_xyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 191);

    auto ts_xyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 192);

    auto ts_xyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 193);

    auto ts_xyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 194);

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

    auto ts_yyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 330);

    auto ts_yyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 331);

    auto ts_yyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 332);

    auto ts_yyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 333);

    auto ts_yyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 334);

    auto ts_yyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 335);

    auto ts_yyzz_xxxxxz = pbuffer.data(idx_ovl_gi + 338);

    auto ts_yyzz_xxxxyz = pbuffer.data(idx_ovl_gi + 340);

    auto ts_yyzz_xxxxzz = pbuffer.data(idx_ovl_gi + 341);

    auto ts_yyzz_xxxyyz = pbuffer.data(idx_ovl_gi + 343);

    auto ts_yyzz_xxxyzz = pbuffer.data(idx_ovl_gi + 344);

    auto ts_yyzz_xxxzzz = pbuffer.data(idx_ovl_gi + 345);

    auto ts_yyzz_xxyyyz = pbuffer.data(idx_ovl_gi + 347);

    auto ts_yyzz_xxyyzz = pbuffer.data(idx_ovl_gi + 348);

    auto ts_yyzz_xxyzzz = pbuffer.data(idx_ovl_gi + 349);

    auto ts_yyzz_xxzzzz = pbuffer.data(idx_ovl_gi + 350);

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

    auto ts_yzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 366);

    auto ts_yzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 369);

    auto ts_yzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 373);

    auto ts_yzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 378);

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

    // Set up components of auxiliary buffer : GI

    auto tr_x_xxxx_xxxxxx = pbuffer.data(idx_dip_gi);

    auto tr_x_xxxx_xxxxxy = pbuffer.data(idx_dip_gi + 1);

    auto tr_x_xxxx_xxxxxz = pbuffer.data(idx_dip_gi + 2);

    auto tr_x_xxxx_xxxxyy = pbuffer.data(idx_dip_gi + 3);

    auto tr_x_xxxx_xxxxyz = pbuffer.data(idx_dip_gi + 4);

    auto tr_x_xxxx_xxxxzz = pbuffer.data(idx_dip_gi + 5);

    auto tr_x_xxxx_xxxyyy = pbuffer.data(idx_dip_gi + 6);

    auto tr_x_xxxx_xxxyyz = pbuffer.data(idx_dip_gi + 7);

    auto tr_x_xxxx_xxxyzz = pbuffer.data(idx_dip_gi + 8);

    auto tr_x_xxxx_xxxzzz = pbuffer.data(idx_dip_gi + 9);

    auto tr_x_xxxx_xxyyyy = pbuffer.data(idx_dip_gi + 10);

    auto tr_x_xxxx_xxyyyz = pbuffer.data(idx_dip_gi + 11);

    auto tr_x_xxxx_xxyyzz = pbuffer.data(idx_dip_gi + 12);

    auto tr_x_xxxx_xxyzzz = pbuffer.data(idx_dip_gi + 13);

    auto tr_x_xxxx_xxzzzz = pbuffer.data(idx_dip_gi + 14);

    auto tr_x_xxxx_xyyyyy = pbuffer.data(idx_dip_gi + 15);

    auto tr_x_xxxx_xyyyyz = pbuffer.data(idx_dip_gi + 16);

    auto tr_x_xxxx_xyyyzz = pbuffer.data(idx_dip_gi + 17);

    auto tr_x_xxxx_xyyzzz = pbuffer.data(idx_dip_gi + 18);

    auto tr_x_xxxx_xyzzzz = pbuffer.data(idx_dip_gi + 19);

    auto tr_x_xxxx_xzzzzz = pbuffer.data(idx_dip_gi + 20);

    auto tr_x_xxxx_yyyyyy = pbuffer.data(idx_dip_gi + 21);

    auto tr_x_xxxx_yyyyyz = pbuffer.data(idx_dip_gi + 22);

    auto tr_x_xxxx_yyyyzz = pbuffer.data(idx_dip_gi + 23);

    auto tr_x_xxxx_yyyzzz = pbuffer.data(idx_dip_gi + 24);

    auto tr_x_xxxx_yyzzzz = pbuffer.data(idx_dip_gi + 25);

    auto tr_x_xxxx_yzzzzz = pbuffer.data(idx_dip_gi + 26);

    auto tr_x_xxxx_zzzzzz = pbuffer.data(idx_dip_gi + 27);

    auto tr_x_xxxy_xxxxxx = pbuffer.data(idx_dip_gi + 28);

    auto tr_x_xxxy_xxxxxy = pbuffer.data(idx_dip_gi + 29);

    auto tr_x_xxxy_xxxxxz = pbuffer.data(idx_dip_gi + 30);

    auto tr_x_xxxy_xxxxyy = pbuffer.data(idx_dip_gi + 31);

    auto tr_x_xxxy_xxxxyz = pbuffer.data(idx_dip_gi + 32);

    auto tr_x_xxxy_xxxxzz = pbuffer.data(idx_dip_gi + 33);

    auto tr_x_xxxy_xxxyyy = pbuffer.data(idx_dip_gi + 34);

    auto tr_x_xxxy_xxxyyz = pbuffer.data(idx_dip_gi + 35);

    auto tr_x_xxxy_xxxyzz = pbuffer.data(idx_dip_gi + 36);

    auto tr_x_xxxy_xxxzzz = pbuffer.data(idx_dip_gi + 37);

    auto tr_x_xxxy_xxyyyy = pbuffer.data(idx_dip_gi + 38);

    auto tr_x_xxxy_xxyyyz = pbuffer.data(idx_dip_gi + 39);

    auto tr_x_xxxy_xxyyzz = pbuffer.data(idx_dip_gi + 40);

    auto tr_x_xxxy_xxyzzz = pbuffer.data(idx_dip_gi + 41);

    auto tr_x_xxxy_xxzzzz = pbuffer.data(idx_dip_gi + 42);

    auto tr_x_xxxy_xyyyyy = pbuffer.data(idx_dip_gi + 43);

    auto tr_x_xxxy_xyyyyz = pbuffer.data(idx_dip_gi + 44);

    auto tr_x_xxxy_xyyyzz = pbuffer.data(idx_dip_gi + 45);

    auto tr_x_xxxy_xyyzzz = pbuffer.data(idx_dip_gi + 46);

    auto tr_x_xxxy_xyzzzz = pbuffer.data(idx_dip_gi + 47);

    auto tr_x_xxxy_xzzzzz = pbuffer.data(idx_dip_gi + 48);

    auto tr_x_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 49);

    auto tr_x_xxxy_zzzzzz = pbuffer.data(idx_dip_gi + 55);

    auto tr_x_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 56);

    auto tr_x_xxxz_xxxxxy = pbuffer.data(idx_dip_gi + 57);

    auto tr_x_xxxz_xxxxxz = pbuffer.data(idx_dip_gi + 58);

    auto tr_x_xxxz_xxxxyy = pbuffer.data(idx_dip_gi + 59);

    auto tr_x_xxxz_xxxxyz = pbuffer.data(idx_dip_gi + 60);

    auto tr_x_xxxz_xxxxzz = pbuffer.data(idx_dip_gi + 61);

    auto tr_x_xxxz_xxxyyy = pbuffer.data(idx_dip_gi + 62);

    auto tr_x_xxxz_xxxyyz = pbuffer.data(idx_dip_gi + 63);

    auto tr_x_xxxz_xxxyzz = pbuffer.data(idx_dip_gi + 64);

    auto tr_x_xxxz_xxxzzz = pbuffer.data(idx_dip_gi + 65);

    auto tr_x_xxxz_xxyyyy = pbuffer.data(idx_dip_gi + 66);

    auto tr_x_xxxz_xxyyyz = pbuffer.data(idx_dip_gi + 67);

    auto tr_x_xxxz_xxyyzz = pbuffer.data(idx_dip_gi + 68);

    auto tr_x_xxxz_xxyzzz = pbuffer.data(idx_dip_gi + 69);

    auto tr_x_xxxz_xxzzzz = pbuffer.data(idx_dip_gi + 70);

    auto tr_x_xxxz_xyyyyy = pbuffer.data(idx_dip_gi + 71);

    auto tr_x_xxxz_xyyyyz = pbuffer.data(idx_dip_gi + 72);

    auto tr_x_xxxz_xyyyzz = pbuffer.data(idx_dip_gi + 73);

    auto tr_x_xxxz_xyyzzz = pbuffer.data(idx_dip_gi + 74);

    auto tr_x_xxxz_xyzzzz = pbuffer.data(idx_dip_gi + 75);

    auto tr_x_xxxz_xzzzzz = pbuffer.data(idx_dip_gi + 76);

    auto tr_x_xxxz_yyyyyy = pbuffer.data(idx_dip_gi + 77);

    auto tr_x_xxxz_yyyyyz = pbuffer.data(idx_dip_gi + 78);

    auto tr_x_xxxz_yyyyzz = pbuffer.data(idx_dip_gi + 79);

    auto tr_x_xxxz_yyyzzz = pbuffer.data(idx_dip_gi + 80);

    auto tr_x_xxxz_yyzzzz = pbuffer.data(idx_dip_gi + 81);

    auto tr_x_xxxz_yzzzzz = pbuffer.data(idx_dip_gi + 82);

    auto tr_x_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 83);

    auto tr_x_xxyy_xxxxxx = pbuffer.data(idx_dip_gi + 84);

    auto tr_x_xxyy_xxxxxy = pbuffer.data(idx_dip_gi + 85);

    auto tr_x_xxyy_xxxxxz = pbuffer.data(idx_dip_gi + 86);

    auto tr_x_xxyy_xxxxyy = pbuffer.data(idx_dip_gi + 87);

    auto tr_x_xxyy_xxxxyz = pbuffer.data(idx_dip_gi + 88);

    auto tr_x_xxyy_xxxxzz = pbuffer.data(idx_dip_gi + 89);

    auto tr_x_xxyy_xxxyyy = pbuffer.data(idx_dip_gi + 90);

    auto tr_x_xxyy_xxxyyz = pbuffer.data(idx_dip_gi + 91);

    auto tr_x_xxyy_xxxyzz = pbuffer.data(idx_dip_gi + 92);

    auto tr_x_xxyy_xxxzzz = pbuffer.data(idx_dip_gi + 93);

    auto tr_x_xxyy_xxyyyy = pbuffer.data(idx_dip_gi + 94);

    auto tr_x_xxyy_xxyyyz = pbuffer.data(idx_dip_gi + 95);

    auto tr_x_xxyy_xxyyzz = pbuffer.data(idx_dip_gi + 96);

    auto tr_x_xxyy_xxyzzz = pbuffer.data(idx_dip_gi + 97);

    auto tr_x_xxyy_xxzzzz = pbuffer.data(idx_dip_gi + 98);

    auto tr_x_xxyy_xyyyyy = pbuffer.data(idx_dip_gi + 99);

    auto tr_x_xxyy_xyyyyz = pbuffer.data(idx_dip_gi + 100);

    auto tr_x_xxyy_xyyyzz = pbuffer.data(idx_dip_gi + 101);

    auto tr_x_xxyy_xyyzzz = pbuffer.data(idx_dip_gi + 102);

    auto tr_x_xxyy_xyzzzz = pbuffer.data(idx_dip_gi + 103);

    auto tr_x_xxyy_xzzzzz = pbuffer.data(idx_dip_gi + 104);

    auto tr_x_xxyy_yyyyyy = pbuffer.data(idx_dip_gi + 105);

    auto tr_x_xxyy_yyyyyz = pbuffer.data(idx_dip_gi + 106);

    auto tr_x_xxyy_yyyyzz = pbuffer.data(idx_dip_gi + 107);

    auto tr_x_xxyy_yyyzzz = pbuffer.data(idx_dip_gi + 108);

    auto tr_x_xxyy_yyzzzz = pbuffer.data(idx_dip_gi + 109);

    auto tr_x_xxyy_yzzzzz = pbuffer.data(idx_dip_gi + 110);

    auto tr_x_xxyy_zzzzzz = pbuffer.data(idx_dip_gi + 111);

    auto tr_x_xxyz_xxxxxz = pbuffer.data(idx_dip_gi + 114);

    auto tr_x_xxyz_xxxxzz = pbuffer.data(idx_dip_gi + 117);

    auto tr_x_xxyz_xxxzzz = pbuffer.data(idx_dip_gi + 121);

    auto tr_x_xxyz_xxzzzz = pbuffer.data(idx_dip_gi + 126);

    auto tr_x_xxyz_xzzzzz = pbuffer.data(idx_dip_gi + 132);

    auto tr_x_xxyz_zzzzzz = pbuffer.data(idx_dip_gi + 139);

    auto tr_x_xxzz_xxxxxx = pbuffer.data(idx_dip_gi + 140);

    auto tr_x_xxzz_xxxxxy = pbuffer.data(idx_dip_gi + 141);

    auto tr_x_xxzz_xxxxxz = pbuffer.data(idx_dip_gi + 142);

    auto tr_x_xxzz_xxxxyy = pbuffer.data(idx_dip_gi + 143);

    auto tr_x_xxzz_xxxxyz = pbuffer.data(idx_dip_gi + 144);

    auto tr_x_xxzz_xxxxzz = pbuffer.data(idx_dip_gi + 145);

    auto tr_x_xxzz_xxxyyy = pbuffer.data(idx_dip_gi + 146);

    auto tr_x_xxzz_xxxyyz = pbuffer.data(idx_dip_gi + 147);

    auto tr_x_xxzz_xxxyzz = pbuffer.data(idx_dip_gi + 148);

    auto tr_x_xxzz_xxxzzz = pbuffer.data(idx_dip_gi + 149);

    auto tr_x_xxzz_xxyyyy = pbuffer.data(idx_dip_gi + 150);

    auto tr_x_xxzz_xxyyyz = pbuffer.data(idx_dip_gi + 151);

    auto tr_x_xxzz_xxyyzz = pbuffer.data(idx_dip_gi + 152);

    auto tr_x_xxzz_xxyzzz = pbuffer.data(idx_dip_gi + 153);

    auto tr_x_xxzz_xxzzzz = pbuffer.data(idx_dip_gi + 154);

    auto tr_x_xxzz_xyyyyy = pbuffer.data(idx_dip_gi + 155);

    auto tr_x_xxzz_xyyyyz = pbuffer.data(idx_dip_gi + 156);

    auto tr_x_xxzz_xyyyzz = pbuffer.data(idx_dip_gi + 157);

    auto tr_x_xxzz_xyyzzz = pbuffer.data(idx_dip_gi + 158);

    auto tr_x_xxzz_xyzzzz = pbuffer.data(idx_dip_gi + 159);

    auto tr_x_xxzz_xzzzzz = pbuffer.data(idx_dip_gi + 160);

    auto tr_x_xxzz_yyyyyy = pbuffer.data(idx_dip_gi + 161);

    auto tr_x_xxzz_yyyyyz = pbuffer.data(idx_dip_gi + 162);

    auto tr_x_xxzz_yyyyzz = pbuffer.data(idx_dip_gi + 163);

    auto tr_x_xxzz_yyyzzz = pbuffer.data(idx_dip_gi + 164);

    auto tr_x_xxzz_yyzzzz = pbuffer.data(idx_dip_gi + 165);

    auto tr_x_xxzz_yzzzzz = pbuffer.data(idx_dip_gi + 166);

    auto tr_x_xxzz_zzzzzz = pbuffer.data(idx_dip_gi + 167);

    auto tr_x_xyyy_xxxxxx = pbuffer.data(idx_dip_gi + 168);

    auto tr_x_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 169);

    auto tr_x_xyyy_xxxxxz = pbuffer.data(idx_dip_gi + 170);

    auto tr_x_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 171);

    auto tr_x_xyyy_xxxxyz = pbuffer.data(idx_dip_gi + 172);

    auto tr_x_xyyy_xxxxzz = pbuffer.data(idx_dip_gi + 173);

    auto tr_x_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 174);

    auto tr_x_xyyy_xxxyyz = pbuffer.data(idx_dip_gi + 175);

    auto tr_x_xyyy_xxxyzz = pbuffer.data(idx_dip_gi + 176);

    auto tr_x_xyyy_xxxzzz = pbuffer.data(idx_dip_gi + 177);

    auto tr_x_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 178);

    auto tr_x_xyyy_xxyyyz = pbuffer.data(idx_dip_gi + 179);

    auto tr_x_xyyy_xxyyzz = pbuffer.data(idx_dip_gi + 180);

    auto tr_x_xyyy_xxyzzz = pbuffer.data(idx_dip_gi + 181);

    auto tr_x_xyyy_xxzzzz = pbuffer.data(idx_dip_gi + 182);

    auto tr_x_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 183);

    auto tr_x_xyyy_xyyyyz = pbuffer.data(idx_dip_gi + 184);

    auto tr_x_xyyy_xyyyzz = pbuffer.data(idx_dip_gi + 185);

    auto tr_x_xyyy_xyyzzz = pbuffer.data(idx_dip_gi + 186);

    auto tr_x_xyyy_xyzzzz = pbuffer.data(idx_dip_gi + 187);

    auto tr_x_xyyy_xzzzzz = pbuffer.data(idx_dip_gi + 188);

    auto tr_x_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 189);

    auto tr_x_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 190);

    auto tr_x_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 191);

    auto tr_x_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 192);

    auto tr_x_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 193);

    auto tr_x_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 194);

    auto tr_x_xyyz_xxxxxy = pbuffer.data(idx_dip_gi + 197);

    auto tr_x_xyyz_xxxxxz = pbuffer.data(idx_dip_gi + 198);

    auto tr_x_xyyz_xxxxyy = pbuffer.data(idx_dip_gi + 199);

    auto tr_x_xyyz_xxxxzz = pbuffer.data(idx_dip_gi + 201);

    auto tr_x_xyyz_xxxyyy = pbuffer.data(idx_dip_gi + 202);

    auto tr_x_xyyz_xxxzzz = pbuffer.data(idx_dip_gi + 205);

    auto tr_x_xyyz_xxyyyy = pbuffer.data(idx_dip_gi + 206);

    auto tr_x_xyyz_xxzzzz = pbuffer.data(idx_dip_gi + 210);

    auto tr_x_xyyz_xyyyyy = pbuffer.data(idx_dip_gi + 211);

    auto tr_x_xyyz_xzzzzz = pbuffer.data(idx_dip_gi + 216);

    auto tr_x_xyzz_xxxxxx = pbuffer.data(idx_dip_gi + 224);

    auto tr_x_xyzz_xxxxxz = pbuffer.data(idx_dip_gi + 226);

    auto tr_x_xyzz_xxxxzz = pbuffer.data(idx_dip_gi + 229);

    auto tr_x_xyzz_xxxzzz = pbuffer.data(idx_dip_gi + 233);

    auto tr_x_xyzz_xxzzzz = pbuffer.data(idx_dip_gi + 238);

    auto tr_x_xyzz_xzzzzz = pbuffer.data(idx_dip_gi + 244);

    auto tr_x_xzzz_xxxxxx = pbuffer.data(idx_dip_gi + 252);

    auto tr_x_xzzz_xxxxxy = pbuffer.data(idx_dip_gi + 253);

    auto tr_x_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 254);

    auto tr_x_xzzz_xxxxyy = pbuffer.data(idx_dip_gi + 255);

    auto tr_x_xzzz_xxxxyz = pbuffer.data(idx_dip_gi + 256);

    auto tr_x_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 257);

    auto tr_x_xzzz_xxxyyy = pbuffer.data(idx_dip_gi + 258);

    auto tr_x_xzzz_xxxyyz = pbuffer.data(idx_dip_gi + 259);

    auto tr_x_xzzz_xxxyzz = pbuffer.data(idx_dip_gi + 260);

    auto tr_x_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 261);

    auto tr_x_xzzz_xxyyyy = pbuffer.data(idx_dip_gi + 262);

    auto tr_x_xzzz_xxyyyz = pbuffer.data(idx_dip_gi + 263);

    auto tr_x_xzzz_xxyyzz = pbuffer.data(idx_dip_gi + 264);

    auto tr_x_xzzz_xxyzzz = pbuffer.data(idx_dip_gi + 265);

    auto tr_x_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 266);

    auto tr_x_xzzz_xyyyyy = pbuffer.data(idx_dip_gi + 267);

    auto tr_x_xzzz_xyyyyz = pbuffer.data(idx_dip_gi + 268);

    auto tr_x_xzzz_xyyyzz = pbuffer.data(idx_dip_gi + 269);

    auto tr_x_xzzz_xyyzzz = pbuffer.data(idx_dip_gi + 270);

    auto tr_x_xzzz_xyzzzz = pbuffer.data(idx_dip_gi + 271);

    auto tr_x_xzzz_xzzzzz = pbuffer.data(idx_dip_gi + 272);

    auto tr_x_xzzz_yyyyyz = pbuffer.data(idx_dip_gi + 274);

    auto tr_x_xzzz_yyyyzz = pbuffer.data(idx_dip_gi + 275);

    auto tr_x_xzzz_yyyzzz = pbuffer.data(idx_dip_gi + 276);

    auto tr_x_xzzz_yyzzzz = pbuffer.data(idx_dip_gi + 277);

    auto tr_x_xzzz_yzzzzz = pbuffer.data(idx_dip_gi + 278);

    auto tr_x_xzzz_zzzzzz = pbuffer.data(idx_dip_gi + 279);

    auto tr_x_yyyy_xxxxxx = pbuffer.data(idx_dip_gi + 280);

    auto tr_x_yyyy_xxxxxy = pbuffer.data(idx_dip_gi + 281);

    auto tr_x_yyyy_xxxxxz = pbuffer.data(idx_dip_gi + 282);

    auto tr_x_yyyy_xxxxyy = pbuffer.data(idx_dip_gi + 283);

    auto tr_x_yyyy_xxxxyz = pbuffer.data(idx_dip_gi + 284);

    auto tr_x_yyyy_xxxxzz = pbuffer.data(idx_dip_gi + 285);

    auto tr_x_yyyy_xxxyyy = pbuffer.data(idx_dip_gi + 286);

    auto tr_x_yyyy_xxxyyz = pbuffer.data(idx_dip_gi + 287);

    auto tr_x_yyyy_xxxyzz = pbuffer.data(idx_dip_gi + 288);

    auto tr_x_yyyy_xxxzzz = pbuffer.data(idx_dip_gi + 289);

    auto tr_x_yyyy_xxyyyy = pbuffer.data(idx_dip_gi + 290);

    auto tr_x_yyyy_xxyyyz = pbuffer.data(idx_dip_gi + 291);

    auto tr_x_yyyy_xxyyzz = pbuffer.data(idx_dip_gi + 292);

    auto tr_x_yyyy_xxyzzz = pbuffer.data(idx_dip_gi + 293);

    auto tr_x_yyyy_xxzzzz = pbuffer.data(idx_dip_gi + 294);

    auto tr_x_yyyy_xyyyyy = pbuffer.data(idx_dip_gi + 295);

    auto tr_x_yyyy_xyyyyz = pbuffer.data(idx_dip_gi + 296);

    auto tr_x_yyyy_xyyyzz = pbuffer.data(idx_dip_gi + 297);

    auto tr_x_yyyy_xyyzzz = pbuffer.data(idx_dip_gi + 298);

    auto tr_x_yyyy_xyzzzz = pbuffer.data(idx_dip_gi + 299);

    auto tr_x_yyyy_xzzzzz = pbuffer.data(idx_dip_gi + 300);

    auto tr_x_yyyy_yyyyyy = pbuffer.data(idx_dip_gi + 301);

    auto tr_x_yyyy_yyyyyz = pbuffer.data(idx_dip_gi + 302);

    auto tr_x_yyyy_yyyyzz = pbuffer.data(idx_dip_gi + 303);

    auto tr_x_yyyy_yyyzzz = pbuffer.data(idx_dip_gi + 304);

    auto tr_x_yyyy_yyzzzz = pbuffer.data(idx_dip_gi + 305);

    auto tr_x_yyyy_yzzzzz = pbuffer.data(idx_dip_gi + 306);

    auto tr_x_yyyy_zzzzzz = pbuffer.data(idx_dip_gi + 307);

    auto tr_x_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 309);

    auto tr_x_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 310);

    auto tr_x_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 311);

    auto tr_x_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 313);

    auto tr_x_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 314);

    auto tr_x_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 317);

    auto tr_x_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 318);

    auto tr_x_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 322);

    auto tr_x_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 323);

    auto tr_x_yyyz_xzzzzz = pbuffer.data(idx_dip_gi + 328);

    auto tr_x_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 329);

    auto tr_x_yyyz_yyyyyz = pbuffer.data(idx_dip_gi + 330);

    auto tr_x_yyyz_yyyyzz = pbuffer.data(idx_dip_gi + 331);

    auto tr_x_yyyz_yyyzzz = pbuffer.data(idx_dip_gi + 332);

    auto tr_x_yyyz_yyzzzz = pbuffer.data(idx_dip_gi + 333);

    auto tr_x_yyyz_yzzzzz = pbuffer.data(idx_dip_gi + 334);

    auto tr_x_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 335);

    auto tr_x_yyzz_xxxxxx = pbuffer.data(idx_dip_gi + 336);

    auto tr_x_yyzz_xxxxxy = pbuffer.data(idx_dip_gi + 337);

    auto tr_x_yyzz_xxxxxz = pbuffer.data(idx_dip_gi + 338);

    auto tr_x_yyzz_xxxxyy = pbuffer.data(idx_dip_gi + 339);

    auto tr_x_yyzz_xxxxyz = pbuffer.data(idx_dip_gi + 340);

    auto tr_x_yyzz_xxxxzz = pbuffer.data(idx_dip_gi + 341);

    auto tr_x_yyzz_xxxyyy = pbuffer.data(idx_dip_gi + 342);

    auto tr_x_yyzz_xxxyyz = pbuffer.data(idx_dip_gi + 343);

    auto tr_x_yyzz_xxxyzz = pbuffer.data(idx_dip_gi + 344);

    auto tr_x_yyzz_xxxzzz = pbuffer.data(idx_dip_gi + 345);

    auto tr_x_yyzz_xxyyyy = pbuffer.data(idx_dip_gi + 346);

    auto tr_x_yyzz_xxyyyz = pbuffer.data(idx_dip_gi + 347);

    auto tr_x_yyzz_xxyyzz = pbuffer.data(idx_dip_gi + 348);

    auto tr_x_yyzz_xxyzzz = pbuffer.data(idx_dip_gi + 349);

    auto tr_x_yyzz_xxzzzz = pbuffer.data(idx_dip_gi + 350);

    auto tr_x_yyzz_xyyyyy = pbuffer.data(idx_dip_gi + 351);

    auto tr_x_yyzz_xyyyyz = pbuffer.data(idx_dip_gi + 352);

    auto tr_x_yyzz_xyyyzz = pbuffer.data(idx_dip_gi + 353);

    auto tr_x_yyzz_xyyzzz = pbuffer.data(idx_dip_gi + 354);

    auto tr_x_yyzz_xyzzzz = pbuffer.data(idx_dip_gi + 355);

    auto tr_x_yyzz_xzzzzz = pbuffer.data(idx_dip_gi + 356);

    auto tr_x_yyzz_yyyyyy = pbuffer.data(idx_dip_gi + 357);

    auto tr_x_yyzz_yyyyyz = pbuffer.data(idx_dip_gi + 358);

    auto tr_x_yyzz_yyyyzz = pbuffer.data(idx_dip_gi + 359);

    auto tr_x_yyzz_yyyzzz = pbuffer.data(idx_dip_gi + 360);

    auto tr_x_yyzz_yyzzzz = pbuffer.data(idx_dip_gi + 361);

    auto tr_x_yyzz_yzzzzz = pbuffer.data(idx_dip_gi + 362);

    auto tr_x_yyzz_zzzzzz = pbuffer.data(idx_dip_gi + 363);

    auto tr_x_yzzz_xxxxxx = pbuffer.data(idx_dip_gi + 364);

    auto tr_x_yzzz_xxxxxz = pbuffer.data(idx_dip_gi + 366);

    auto tr_x_yzzz_xxxxyz = pbuffer.data(idx_dip_gi + 368);

    auto tr_x_yzzz_xxxxzz = pbuffer.data(idx_dip_gi + 369);

    auto tr_x_yzzz_xxxyyz = pbuffer.data(idx_dip_gi + 371);

    auto tr_x_yzzz_xxxyzz = pbuffer.data(idx_dip_gi + 372);

    auto tr_x_yzzz_xxxzzz = pbuffer.data(idx_dip_gi + 373);

    auto tr_x_yzzz_xxyyyz = pbuffer.data(idx_dip_gi + 375);

    auto tr_x_yzzz_xxyyzz = pbuffer.data(idx_dip_gi + 376);

    auto tr_x_yzzz_xxyzzz = pbuffer.data(idx_dip_gi + 377);

    auto tr_x_yzzz_xxzzzz = pbuffer.data(idx_dip_gi + 378);

    auto tr_x_yzzz_xyyyyz = pbuffer.data(idx_dip_gi + 380);

    auto tr_x_yzzz_xyyyzz = pbuffer.data(idx_dip_gi + 381);

    auto tr_x_yzzz_xyyzzz = pbuffer.data(idx_dip_gi + 382);

    auto tr_x_yzzz_xyzzzz = pbuffer.data(idx_dip_gi + 383);

    auto tr_x_yzzz_xzzzzz = pbuffer.data(idx_dip_gi + 384);

    auto tr_x_yzzz_yyyyyy = pbuffer.data(idx_dip_gi + 385);

    auto tr_x_yzzz_yyyyyz = pbuffer.data(idx_dip_gi + 386);

    auto tr_x_yzzz_yyyyzz = pbuffer.data(idx_dip_gi + 387);

    auto tr_x_yzzz_yyyzzz = pbuffer.data(idx_dip_gi + 388);

    auto tr_x_yzzz_yyzzzz = pbuffer.data(idx_dip_gi + 389);

    auto tr_x_yzzz_yzzzzz = pbuffer.data(idx_dip_gi + 390);

    auto tr_x_yzzz_zzzzzz = pbuffer.data(idx_dip_gi + 391);

    auto tr_x_zzzz_xxxxxx = pbuffer.data(idx_dip_gi + 392);

    auto tr_x_zzzz_xxxxxy = pbuffer.data(idx_dip_gi + 393);

    auto tr_x_zzzz_xxxxxz = pbuffer.data(idx_dip_gi + 394);

    auto tr_x_zzzz_xxxxyy = pbuffer.data(idx_dip_gi + 395);

    auto tr_x_zzzz_xxxxyz = pbuffer.data(idx_dip_gi + 396);

    auto tr_x_zzzz_xxxxzz = pbuffer.data(idx_dip_gi + 397);

    auto tr_x_zzzz_xxxyyy = pbuffer.data(idx_dip_gi + 398);

    auto tr_x_zzzz_xxxyyz = pbuffer.data(idx_dip_gi + 399);

    auto tr_x_zzzz_xxxyzz = pbuffer.data(idx_dip_gi + 400);

    auto tr_x_zzzz_xxxzzz = pbuffer.data(idx_dip_gi + 401);

    auto tr_x_zzzz_xxyyyy = pbuffer.data(idx_dip_gi + 402);

    auto tr_x_zzzz_xxyyyz = pbuffer.data(idx_dip_gi + 403);

    auto tr_x_zzzz_xxyyzz = pbuffer.data(idx_dip_gi + 404);

    auto tr_x_zzzz_xxyzzz = pbuffer.data(idx_dip_gi + 405);

    auto tr_x_zzzz_xxzzzz = pbuffer.data(idx_dip_gi + 406);

    auto tr_x_zzzz_xyyyyy = pbuffer.data(idx_dip_gi + 407);

    auto tr_x_zzzz_xyyyyz = pbuffer.data(idx_dip_gi + 408);

    auto tr_x_zzzz_xyyyzz = pbuffer.data(idx_dip_gi + 409);

    auto tr_x_zzzz_xyyzzz = pbuffer.data(idx_dip_gi + 410);

    auto tr_x_zzzz_xyzzzz = pbuffer.data(idx_dip_gi + 411);

    auto tr_x_zzzz_xzzzzz = pbuffer.data(idx_dip_gi + 412);

    auto tr_x_zzzz_yyyyyy = pbuffer.data(idx_dip_gi + 413);

    auto tr_x_zzzz_yyyyyz = pbuffer.data(idx_dip_gi + 414);

    auto tr_x_zzzz_yyyyzz = pbuffer.data(idx_dip_gi + 415);

    auto tr_x_zzzz_yyyzzz = pbuffer.data(idx_dip_gi + 416);

    auto tr_x_zzzz_yyzzzz = pbuffer.data(idx_dip_gi + 417);

    auto tr_x_zzzz_yzzzzz = pbuffer.data(idx_dip_gi + 418);

    auto tr_x_zzzz_zzzzzz = pbuffer.data(idx_dip_gi + 419);

    auto tr_y_xxxx_xxxxxx = pbuffer.data(idx_dip_gi + 420);

    auto tr_y_xxxx_xxxxxy = pbuffer.data(idx_dip_gi + 421);

    auto tr_y_xxxx_xxxxxz = pbuffer.data(idx_dip_gi + 422);

    auto tr_y_xxxx_xxxxyy = pbuffer.data(idx_dip_gi + 423);

    auto tr_y_xxxx_xxxxyz = pbuffer.data(idx_dip_gi + 424);

    auto tr_y_xxxx_xxxxzz = pbuffer.data(idx_dip_gi + 425);

    auto tr_y_xxxx_xxxyyy = pbuffer.data(idx_dip_gi + 426);

    auto tr_y_xxxx_xxxyyz = pbuffer.data(idx_dip_gi + 427);

    auto tr_y_xxxx_xxxyzz = pbuffer.data(idx_dip_gi + 428);

    auto tr_y_xxxx_xxxzzz = pbuffer.data(idx_dip_gi + 429);

    auto tr_y_xxxx_xxyyyy = pbuffer.data(idx_dip_gi + 430);

    auto tr_y_xxxx_xxyyyz = pbuffer.data(idx_dip_gi + 431);

    auto tr_y_xxxx_xxyyzz = pbuffer.data(idx_dip_gi + 432);

    auto tr_y_xxxx_xxyzzz = pbuffer.data(idx_dip_gi + 433);

    auto tr_y_xxxx_xxzzzz = pbuffer.data(idx_dip_gi + 434);

    auto tr_y_xxxx_xyyyyy = pbuffer.data(idx_dip_gi + 435);

    auto tr_y_xxxx_xyyyyz = pbuffer.data(idx_dip_gi + 436);

    auto tr_y_xxxx_xyyyzz = pbuffer.data(idx_dip_gi + 437);

    auto tr_y_xxxx_xyyzzz = pbuffer.data(idx_dip_gi + 438);

    auto tr_y_xxxx_xyzzzz = pbuffer.data(idx_dip_gi + 439);

    auto tr_y_xxxx_xzzzzz = pbuffer.data(idx_dip_gi + 440);

    auto tr_y_xxxx_yyyyyy = pbuffer.data(idx_dip_gi + 441);

    auto tr_y_xxxx_yyyyyz = pbuffer.data(idx_dip_gi + 442);

    auto tr_y_xxxx_yyyyzz = pbuffer.data(idx_dip_gi + 443);

    auto tr_y_xxxx_yyyzzz = pbuffer.data(idx_dip_gi + 444);

    auto tr_y_xxxx_yyzzzz = pbuffer.data(idx_dip_gi + 445);

    auto tr_y_xxxx_yzzzzz = pbuffer.data(idx_dip_gi + 446);

    auto tr_y_xxxx_zzzzzz = pbuffer.data(idx_dip_gi + 447);

    auto tr_y_xxxy_xxxxxx = pbuffer.data(idx_dip_gi + 448);

    auto tr_y_xxxy_xxxxxy = pbuffer.data(idx_dip_gi + 449);

    auto tr_y_xxxy_xxxxyy = pbuffer.data(idx_dip_gi + 451);

    auto tr_y_xxxy_xxxxyz = pbuffer.data(idx_dip_gi + 452);

    auto tr_y_xxxy_xxxyyy = pbuffer.data(idx_dip_gi + 454);

    auto tr_y_xxxy_xxxyyz = pbuffer.data(idx_dip_gi + 455);

    auto tr_y_xxxy_xxxyzz = pbuffer.data(idx_dip_gi + 456);

    auto tr_y_xxxy_xxyyyy = pbuffer.data(idx_dip_gi + 458);

    auto tr_y_xxxy_xxyyyz = pbuffer.data(idx_dip_gi + 459);

    auto tr_y_xxxy_xxyyzz = pbuffer.data(idx_dip_gi + 460);

    auto tr_y_xxxy_xxyzzz = pbuffer.data(idx_dip_gi + 461);

    auto tr_y_xxxy_xyyyyy = pbuffer.data(idx_dip_gi + 463);

    auto tr_y_xxxy_xyyyyz = pbuffer.data(idx_dip_gi + 464);

    auto tr_y_xxxy_xyyyzz = pbuffer.data(idx_dip_gi + 465);

    auto tr_y_xxxy_xyyzzz = pbuffer.data(idx_dip_gi + 466);

    auto tr_y_xxxy_xyzzzz = pbuffer.data(idx_dip_gi + 467);

    auto tr_y_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 469);

    auto tr_y_xxxy_yyyyyz = pbuffer.data(idx_dip_gi + 470);

    auto tr_y_xxxy_yyyyzz = pbuffer.data(idx_dip_gi + 471);

    auto tr_y_xxxy_yyyzzz = pbuffer.data(idx_dip_gi + 472);

    auto tr_y_xxxy_yyzzzz = pbuffer.data(idx_dip_gi + 473);

    auto tr_y_xxxy_yzzzzz = pbuffer.data(idx_dip_gi + 474);

    auto tr_y_xxxy_zzzzzz = pbuffer.data(idx_dip_gi + 475);

    auto tr_y_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 476);

    auto tr_y_xxxz_xxxxxy = pbuffer.data(idx_dip_gi + 477);

    auto tr_y_xxxz_xxxxxz = pbuffer.data(idx_dip_gi + 478);

    auto tr_y_xxxz_xxxxyy = pbuffer.data(idx_dip_gi + 479);

    auto tr_y_xxxz_xxxxzz = pbuffer.data(idx_dip_gi + 481);

    auto tr_y_xxxz_xxxyyy = pbuffer.data(idx_dip_gi + 482);

    auto tr_y_xxxz_xxxzzz = pbuffer.data(idx_dip_gi + 485);

    auto tr_y_xxxz_xxyyyy = pbuffer.data(idx_dip_gi + 486);

    auto tr_y_xxxz_xxzzzz = pbuffer.data(idx_dip_gi + 490);

    auto tr_y_xxxz_xyyyyy = pbuffer.data(idx_dip_gi + 491);

    auto tr_y_xxxz_xzzzzz = pbuffer.data(idx_dip_gi + 496);

    auto tr_y_xxxz_yyyyyz = pbuffer.data(idx_dip_gi + 498);

    auto tr_y_xxxz_yyyyzz = pbuffer.data(idx_dip_gi + 499);

    auto tr_y_xxxz_yyyzzz = pbuffer.data(idx_dip_gi + 500);

    auto tr_y_xxxz_yyzzzz = pbuffer.data(idx_dip_gi + 501);

    auto tr_y_xxxz_yzzzzz = pbuffer.data(idx_dip_gi + 502);

    auto tr_y_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 503);

    auto tr_y_xxyy_xxxxxx = pbuffer.data(idx_dip_gi + 504);

    auto tr_y_xxyy_xxxxxy = pbuffer.data(idx_dip_gi + 505);

    auto tr_y_xxyy_xxxxxz = pbuffer.data(idx_dip_gi + 506);

    auto tr_y_xxyy_xxxxyy = pbuffer.data(idx_dip_gi + 507);

    auto tr_y_xxyy_xxxxyz = pbuffer.data(idx_dip_gi + 508);

    auto tr_y_xxyy_xxxxzz = pbuffer.data(idx_dip_gi + 509);

    auto tr_y_xxyy_xxxyyy = pbuffer.data(idx_dip_gi + 510);

    auto tr_y_xxyy_xxxyyz = pbuffer.data(idx_dip_gi + 511);

    auto tr_y_xxyy_xxxyzz = pbuffer.data(idx_dip_gi + 512);

    auto tr_y_xxyy_xxxzzz = pbuffer.data(idx_dip_gi + 513);

    auto tr_y_xxyy_xxyyyy = pbuffer.data(idx_dip_gi + 514);

    auto tr_y_xxyy_xxyyyz = pbuffer.data(idx_dip_gi + 515);

    auto tr_y_xxyy_xxyyzz = pbuffer.data(idx_dip_gi + 516);

    auto tr_y_xxyy_xxyzzz = pbuffer.data(idx_dip_gi + 517);

    auto tr_y_xxyy_xxzzzz = pbuffer.data(idx_dip_gi + 518);

    auto tr_y_xxyy_xyyyyy = pbuffer.data(idx_dip_gi + 519);

    auto tr_y_xxyy_xyyyyz = pbuffer.data(idx_dip_gi + 520);

    auto tr_y_xxyy_xyyyzz = pbuffer.data(idx_dip_gi + 521);

    auto tr_y_xxyy_xyyzzz = pbuffer.data(idx_dip_gi + 522);

    auto tr_y_xxyy_xyzzzz = pbuffer.data(idx_dip_gi + 523);

    auto tr_y_xxyy_xzzzzz = pbuffer.data(idx_dip_gi + 524);

    auto tr_y_xxyy_yyyyyy = pbuffer.data(idx_dip_gi + 525);

    auto tr_y_xxyy_yyyyyz = pbuffer.data(idx_dip_gi + 526);

    auto tr_y_xxyy_yyyyzz = pbuffer.data(idx_dip_gi + 527);

    auto tr_y_xxyy_yyyzzz = pbuffer.data(idx_dip_gi + 528);

    auto tr_y_xxyy_yyzzzz = pbuffer.data(idx_dip_gi + 529);

    auto tr_y_xxyy_yzzzzz = pbuffer.data(idx_dip_gi + 530);

    auto tr_y_xxyy_zzzzzz = pbuffer.data(idx_dip_gi + 531);

    auto tr_y_xxyz_xxxxxy = pbuffer.data(idx_dip_gi + 533);

    auto tr_y_xxyz_xxxxyy = pbuffer.data(idx_dip_gi + 535);

    auto tr_y_xxyz_xxxyyy = pbuffer.data(idx_dip_gi + 538);

    auto tr_y_xxyz_xxyyyy = pbuffer.data(idx_dip_gi + 542);

    auto tr_y_xxyz_xyyyyy = pbuffer.data(idx_dip_gi + 547);

    auto tr_y_xxyz_yyyyyz = pbuffer.data(idx_dip_gi + 554);

    auto tr_y_xxyz_yyyyzz = pbuffer.data(idx_dip_gi + 555);

    auto tr_y_xxyz_yyyzzz = pbuffer.data(idx_dip_gi + 556);

    auto tr_y_xxyz_yyzzzz = pbuffer.data(idx_dip_gi + 557);

    auto tr_y_xxyz_yzzzzz = pbuffer.data(idx_dip_gi + 558);

    auto tr_y_xxyz_zzzzzz = pbuffer.data(idx_dip_gi + 559);

    auto tr_y_xxzz_xxxxxx = pbuffer.data(idx_dip_gi + 560);

    auto tr_y_xxzz_xxxxxy = pbuffer.data(idx_dip_gi + 561);

    auto tr_y_xxzz_xxxxxz = pbuffer.data(idx_dip_gi + 562);

    auto tr_y_xxzz_xxxxyy = pbuffer.data(idx_dip_gi + 563);

    auto tr_y_xxzz_xxxxyz = pbuffer.data(idx_dip_gi + 564);

    auto tr_y_xxzz_xxxxzz = pbuffer.data(idx_dip_gi + 565);

    auto tr_y_xxzz_xxxyyy = pbuffer.data(idx_dip_gi + 566);

    auto tr_y_xxzz_xxxyyz = pbuffer.data(idx_dip_gi + 567);

    auto tr_y_xxzz_xxxyzz = pbuffer.data(idx_dip_gi + 568);

    auto tr_y_xxzz_xxxzzz = pbuffer.data(idx_dip_gi + 569);

    auto tr_y_xxzz_xxyyyy = pbuffer.data(idx_dip_gi + 570);

    auto tr_y_xxzz_xxyyyz = pbuffer.data(idx_dip_gi + 571);

    auto tr_y_xxzz_xxyyzz = pbuffer.data(idx_dip_gi + 572);

    auto tr_y_xxzz_xxyzzz = pbuffer.data(idx_dip_gi + 573);

    auto tr_y_xxzz_xxzzzz = pbuffer.data(idx_dip_gi + 574);

    auto tr_y_xxzz_xyyyyy = pbuffer.data(idx_dip_gi + 575);

    auto tr_y_xxzz_xyyyyz = pbuffer.data(idx_dip_gi + 576);

    auto tr_y_xxzz_xyyyzz = pbuffer.data(idx_dip_gi + 577);

    auto tr_y_xxzz_xyyzzz = pbuffer.data(idx_dip_gi + 578);

    auto tr_y_xxzz_xyzzzz = pbuffer.data(idx_dip_gi + 579);

    auto tr_y_xxzz_xzzzzz = pbuffer.data(idx_dip_gi + 580);

    auto tr_y_xxzz_yyyyyy = pbuffer.data(idx_dip_gi + 581);

    auto tr_y_xxzz_yyyyyz = pbuffer.data(idx_dip_gi + 582);

    auto tr_y_xxzz_yyyyzz = pbuffer.data(idx_dip_gi + 583);

    auto tr_y_xxzz_yyyzzz = pbuffer.data(idx_dip_gi + 584);

    auto tr_y_xxzz_yyzzzz = pbuffer.data(idx_dip_gi + 585);

    auto tr_y_xxzz_yzzzzz = pbuffer.data(idx_dip_gi + 586);

    auto tr_y_xxzz_zzzzzz = pbuffer.data(idx_dip_gi + 587);

    auto tr_y_xyyy_xxxxxx = pbuffer.data(idx_dip_gi + 588);

    auto tr_y_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 589);

    auto tr_y_xyyy_xxxxxz = pbuffer.data(idx_dip_gi + 590);

    auto tr_y_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 591);

    auto tr_y_xyyy_xxxxyz = pbuffer.data(idx_dip_gi + 592);

    auto tr_y_xyyy_xxxxzz = pbuffer.data(idx_dip_gi + 593);

    auto tr_y_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 594);

    auto tr_y_xyyy_xxxyyz = pbuffer.data(idx_dip_gi + 595);

    auto tr_y_xyyy_xxxyzz = pbuffer.data(idx_dip_gi + 596);

    auto tr_y_xyyy_xxxzzz = pbuffer.data(idx_dip_gi + 597);

    auto tr_y_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 598);

    auto tr_y_xyyy_xxyyyz = pbuffer.data(idx_dip_gi + 599);

    auto tr_y_xyyy_xxyyzz = pbuffer.data(idx_dip_gi + 600);

    auto tr_y_xyyy_xxyzzz = pbuffer.data(idx_dip_gi + 601);

    auto tr_y_xyyy_xxzzzz = pbuffer.data(idx_dip_gi + 602);

    auto tr_y_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 603);

    auto tr_y_xyyy_xyyyyz = pbuffer.data(idx_dip_gi + 604);

    auto tr_y_xyyy_xyyyzz = pbuffer.data(idx_dip_gi + 605);

    auto tr_y_xyyy_xyyzzz = pbuffer.data(idx_dip_gi + 606);

    auto tr_y_xyyy_xyzzzz = pbuffer.data(idx_dip_gi + 607);

    auto tr_y_xyyy_xzzzzz = pbuffer.data(idx_dip_gi + 608);

    auto tr_y_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 609);

    auto tr_y_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 610);

    auto tr_y_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 611);

    auto tr_y_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 612);

    auto tr_y_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 613);

    auto tr_y_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 614);

    auto tr_y_xyyy_zzzzzz = pbuffer.data(idx_dip_gi + 615);

    auto tr_y_xyyz_yyyyyz = pbuffer.data(idx_dip_gi + 638);

    auto tr_y_xyyz_yyyyzz = pbuffer.data(idx_dip_gi + 639);

    auto tr_y_xyyz_yyyzzz = pbuffer.data(idx_dip_gi + 640);

    auto tr_y_xyyz_yyzzzz = pbuffer.data(idx_dip_gi + 641);

    auto tr_y_xyyz_yzzzzz = pbuffer.data(idx_dip_gi + 642);

    auto tr_y_xyyz_zzzzzz = pbuffer.data(idx_dip_gi + 643);

    auto tr_y_xyzz_xxxxyz = pbuffer.data(idx_dip_gi + 648);

    auto tr_y_xyzz_xxxyyz = pbuffer.data(idx_dip_gi + 651);

    auto tr_y_xyzz_xxxyzz = pbuffer.data(idx_dip_gi + 652);

    auto tr_y_xyzz_xxyyyz = pbuffer.data(idx_dip_gi + 655);

    auto tr_y_xyzz_xxyyzz = pbuffer.data(idx_dip_gi + 656);

    auto tr_y_xyzz_xxyzzz = pbuffer.data(idx_dip_gi + 657);

    auto tr_y_xyzz_xyyyyz = pbuffer.data(idx_dip_gi + 660);

    auto tr_y_xyzz_xyyyzz = pbuffer.data(idx_dip_gi + 661);

    auto tr_y_xyzz_xyyzzz = pbuffer.data(idx_dip_gi + 662);

    auto tr_y_xyzz_xyzzzz = pbuffer.data(idx_dip_gi + 663);

    auto tr_y_xyzz_yyyyyy = pbuffer.data(idx_dip_gi + 665);

    auto tr_y_xyzz_yyyyyz = pbuffer.data(idx_dip_gi + 666);

    auto tr_y_xyzz_yyyyzz = pbuffer.data(idx_dip_gi + 667);

    auto tr_y_xyzz_yyyzzz = pbuffer.data(idx_dip_gi + 668);

    auto tr_y_xyzz_yyzzzz = pbuffer.data(idx_dip_gi + 669);

    auto tr_y_xyzz_yzzzzz = pbuffer.data(idx_dip_gi + 670);

    auto tr_y_xyzz_zzzzzz = pbuffer.data(idx_dip_gi + 671);

    auto tr_y_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 674);

    auto tr_y_xzzz_xxxxyz = pbuffer.data(idx_dip_gi + 676);

    auto tr_y_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 677);

    auto tr_y_xzzz_xxxyyz = pbuffer.data(idx_dip_gi + 679);

    auto tr_y_xzzz_xxxyzz = pbuffer.data(idx_dip_gi + 680);

    auto tr_y_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 681);

    auto tr_y_xzzz_xxyyyz = pbuffer.data(idx_dip_gi + 683);

    auto tr_y_xzzz_xxyyzz = pbuffer.data(idx_dip_gi + 684);

    auto tr_y_xzzz_xxyzzz = pbuffer.data(idx_dip_gi + 685);

    auto tr_y_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 686);

    auto tr_y_xzzz_xyyyyz = pbuffer.data(idx_dip_gi + 688);

    auto tr_y_xzzz_xyyyzz = pbuffer.data(idx_dip_gi + 689);

    auto tr_y_xzzz_xyyzzz = pbuffer.data(idx_dip_gi + 690);

    auto tr_y_xzzz_xyzzzz = pbuffer.data(idx_dip_gi + 691);

    auto tr_y_xzzz_xzzzzz = pbuffer.data(idx_dip_gi + 692);

    auto tr_y_xzzz_yyyyyy = pbuffer.data(idx_dip_gi + 693);

    auto tr_y_xzzz_yyyyyz = pbuffer.data(idx_dip_gi + 694);

    auto tr_y_xzzz_yyyyzz = pbuffer.data(idx_dip_gi + 695);

    auto tr_y_xzzz_yyyzzz = pbuffer.data(idx_dip_gi + 696);

    auto tr_y_xzzz_yyzzzz = pbuffer.data(idx_dip_gi + 697);

    auto tr_y_xzzz_yzzzzz = pbuffer.data(idx_dip_gi + 698);

    auto tr_y_xzzz_zzzzzz = pbuffer.data(idx_dip_gi + 699);

    auto tr_y_yyyy_xxxxxx = pbuffer.data(idx_dip_gi + 700);

    auto tr_y_yyyy_xxxxxy = pbuffer.data(idx_dip_gi + 701);

    auto tr_y_yyyy_xxxxxz = pbuffer.data(idx_dip_gi + 702);

    auto tr_y_yyyy_xxxxyy = pbuffer.data(idx_dip_gi + 703);

    auto tr_y_yyyy_xxxxyz = pbuffer.data(idx_dip_gi + 704);

    auto tr_y_yyyy_xxxxzz = pbuffer.data(idx_dip_gi + 705);

    auto tr_y_yyyy_xxxyyy = pbuffer.data(idx_dip_gi + 706);

    auto tr_y_yyyy_xxxyyz = pbuffer.data(idx_dip_gi + 707);

    auto tr_y_yyyy_xxxyzz = pbuffer.data(idx_dip_gi + 708);

    auto tr_y_yyyy_xxxzzz = pbuffer.data(idx_dip_gi + 709);

    auto tr_y_yyyy_xxyyyy = pbuffer.data(idx_dip_gi + 710);

    auto tr_y_yyyy_xxyyyz = pbuffer.data(idx_dip_gi + 711);

    auto tr_y_yyyy_xxyyzz = pbuffer.data(idx_dip_gi + 712);

    auto tr_y_yyyy_xxyzzz = pbuffer.data(idx_dip_gi + 713);

    auto tr_y_yyyy_xxzzzz = pbuffer.data(idx_dip_gi + 714);

    auto tr_y_yyyy_xyyyyy = pbuffer.data(idx_dip_gi + 715);

    auto tr_y_yyyy_xyyyyz = pbuffer.data(idx_dip_gi + 716);

    auto tr_y_yyyy_xyyyzz = pbuffer.data(idx_dip_gi + 717);

    auto tr_y_yyyy_xyyzzz = pbuffer.data(idx_dip_gi + 718);

    auto tr_y_yyyy_xyzzzz = pbuffer.data(idx_dip_gi + 719);

    auto tr_y_yyyy_xzzzzz = pbuffer.data(idx_dip_gi + 720);

    auto tr_y_yyyy_yyyyyy = pbuffer.data(idx_dip_gi + 721);

    auto tr_y_yyyy_yyyyyz = pbuffer.data(idx_dip_gi + 722);

    auto tr_y_yyyy_yyyyzz = pbuffer.data(idx_dip_gi + 723);

    auto tr_y_yyyy_yyyzzz = pbuffer.data(idx_dip_gi + 724);

    auto tr_y_yyyy_yyzzzz = pbuffer.data(idx_dip_gi + 725);

    auto tr_y_yyyy_yzzzzz = pbuffer.data(idx_dip_gi + 726);

    auto tr_y_yyyy_zzzzzz = pbuffer.data(idx_dip_gi + 727);

    auto tr_y_yyyz_xxxxxx = pbuffer.data(idx_dip_gi + 728);

    auto tr_y_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 729);

    auto tr_y_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 730);

    auto tr_y_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 731);

    auto tr_y_yyyz_xxxxyz = pbuffer.data(idx_dip_gi + 732);

    auto tr_y_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 733);

    auto tr_y_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 734);

    auto tr_y_yyyz_xxxyyz = pbuffer.data(idx_dip_gi + 735);

    auto tr_y_yyyz_xxxyzz = pbuffer.data(idx_dip_gi + 736);

    auto tr_y_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 737);

    auto tr_y_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 738);

    auto tr_y_yyyz_xxyyyz = pbuffer.data(idx_dip_gi + 739);

    auto tr_y_yyyz_xxyyzz = pbuffer.data(idx_dip_gi + 740);

    auto tr_y_yyyz_xxyzzz = pbuffer.data(idx_dip_gi + 741);

    auto tr_y_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 742);

    auto tr_y_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 743);

    auto tr_y_yyyz_xyyyyz = pbuffer.data(idx_dip_gi + 744);

    auto tr_y_yyyz_xyyyzz = pbuffer.data(idx_dip_gi + 745);

    auto tr_y_yyyz_xyyzzz = pbuffer.data(idx_dip_gi + 746);

    auto tr_y_yyyz_xyzzzz = pbuffer.data(idx_dip_gi + 747);

    auto tr_y_yyyz_xzzzzz = pbuffer.data(idx_dip_gi + 748);

    auto tr_y_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 749);

    auto tr_y_yyyz_yyyyyz = pbuffer.data(idx_dip_gi + 750);

    auto tr_y_yyyz_yyyyzz = pbuffer.data(idx_dip_gi + 751);

    auto tr_y_yyyz_yyyzzz = pbuffer.data(idx_dip_gi + 752);

    auto tr_y_yyyz_yyzzzz = pbuffer.data(idx_dip_gi + 753);

    auto tr_y_yyyz_yzzzzz = pbuffer.data(idx_dip_gi + 754);

    auto tr_y_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 755);

    auto tr_y_yyzz_xxxxxx = pbuffer.data(idx_dip_gi + 756);

    auto tr_y_yyzz_xxxxxy = pbuffer.data(idx_dip_gi + 757);

    auto tr_y_yyzz_xxxxxz = pbuffer.data(idx_dip_gi + 758);

    auto tr_y_yyzz_xxxxyy = pbuffer.data(idx_dip_gi + 759);

    auto tr_y_yyzz_xxxxyz = pbuffer.data(idx_dip_gi + 760);

    auto tr_y_yyzz_xxxxzz = pbuffer.data(idx_dip_gi + 761);

    auto tr_y_yyzz_xxxyyy = pbuffer.data(idx_dip_gi + 762);

    auto tr_y_yyzz_xxxyyz = pbuffer.data(idx_dip_gi + 763);

    auto tr_y_yyzz_xxxyzz = pbuffer.data(idx_dip_gi + 764);

    auto tr_y_yyzz_xxxzzz = pbuffer.data(idx_dip_gi + 765);

    auto tr_y_yyzz_xxyyyy = pbuffer.data(idx_dip_gi + 766);

    auto tr_y_yyzz_xxyyyz = pbuffer.data(idx_dip_gi + 767);

    auto tr_y_yyzz_xxyyzz = pbuffer.data(idx_dip_gi + 768);

    auto tr_y_yyzz_xxyzzz = pbuffer.data(idx_dip_gi + 769);

    auto tr_y_yyzz_xxzzzz = pbuffer.data(idx_dip_gi + 770);

    auto tr_y_yyzz_xyyyyy = pbuffer.data(idx_dip_gi + 771);

    auto tr_y_yyzz_xyyyyz = pbuffer.data(idx_dip_gi + 772);

    auto tr_y_yyzz_xyyyzz = pbuffer.data(idx_dip_gi + 773);

    auto tr_y_yyzz_xyyzzz = pbuffer.data(idx_dip_gi + 774);

    auto tr_y_yyzz_xyzzzz = pbuffer.data(idx_dip_gi + 775);

    auto tr_y_yyzz_xzzzzz = pbuffer.data(idx_dip_gi + 776);

    auto tr_y_yyzz_yyyyyy = pbuffer.data(idx_dip_gi + 777);

    auto tr_y_yyzz_yyyyyz = pbuffer.data(idx_dip_gi + 778);

    auto tr_y_yyzz_yyyyzz = pbuffer.data(idx_dip_gi + 779);

    auto tr_y_yyzz_yyyzzz = pbuffer.data(idx_dip_gi + 780);

    auto tr_y_yyzz_yyzzzz = pbuffer.data(idx_dip_gi + 781);

    auto tr_y_yyzz_yzzzzz = pbuffer.data(idx_dip_gi + 782);

    auto tr_y_yyzz_zzzzzz = pbuffer.data(idx_dip_gi + 783);

    auto tr_y_yzzz_xxxxxx = pbuffer.data(idx_dip_gi + 784);

    auto tr_y_yzzz_xxxxxy = pbuffer.data(idx_dip_gi + 785);

    auto tr_y_yzzz_xxxxxz = pbuffer.data(idx_dip_gi + 786);

    auto tr_y_yzzz_xxxxyy = pbuffer.data(idx_dip_gi + 787);

    auto tr_y_yzzz_xxxxyz = pbuffer.data(idx_dip_gi + 788);

    auto tr_y_yzzz_xxxxzz = pbuffer.data(idx_dip_gi + 789);

    auto tr_y_yzzz_xxxyyy = pbuffer.data(idx_dip_gi + 790);

    auto tr_y_yzzz_xxxyyz = pbuffer.data(idx_dip_gi + 791);

    auto tr_y_yzzz_xxxyzz = pbuffer.data(idx_dip_gi + 792);

    auto tr_y_yzzz_xxxzzz = pbuffer.data(idx_dip_gi + 793);

    auto tr_y_yzzz_xxyyyy = pbuffer.data(idx_dip_gi + 794);

    auto tr_y_yzzz_xxyyyz = pbuffer.data(idx_dip_gi + 795);

    auto tr_y_yzzz_xxyyzz = pbuffer.data(idx_dip_gi + 796);

    auto tr_y_yzzz_xxyzzz = pbuffer.data(idx_dip_gi + 797);

    auto tr_y_yzzz_xxzzzz = pbuffer.data(idx_dip_gi + 798);

    auto tr_y_yzzz_xyyyyy = pbuffer.data(idx_dip_gi + 799);

    auto tr_y_yzzz_xyyyyz = pbuffer.data(idx_dip_gi + 800);

    auto tr_y_yzzz_xyyyzz = pbuffer.data(idx_dip_gi + 801);

    auto tr_y_yzzz_xyyzzz = pbuffer.data(idx_dip_gi + 802);

    auto tr_y_yzzz_xyzzzz = pbuffer.data(idx_dip_gi + 803);

    auto tr_y_yzzz_xzzzzz = pbuffer.data(idx_dip_gi + 804);

    auto tr_y_yzzz_yyyyyy = pbuffer.data(idx_dip_gi + 805);

    auto tr_y_yzzz_yyyyyz = pbuffer.data(idx_dip_gi + 806);

    auto tr_y_yzzz_yyyyzz = pbuffer.data(idx_dip_gi + 807);

    auto tr_y_yzzz_yyyzzz = pbuffer.data(idx_dip_gi + 808);

    auto tr_y_yzzz_yyzzzz = pbuffer.data(idx_dip_gi + 809);

    auto tr_y_yzzz_yzzzzz = pbuffer.data(idx_dip_gi + 810);

    auto tr_y_yzzz_zzzzzz = pbuffer.data(idx_dip_gi + 811);

    auto tr_y_zzzz_xxxxxx = pbuffer.data(idx_dip_gi + 812);

    auto tr_y_zzzz_xxxxxy = pbuffer.data(idx_dip_gi + 813);

    auto tr_y_zzzz_xxxxxz = pbuffer.data(idx_dip_gi + 814);

    auto tr_y_zzzz_xxxxyy = pbuffer.data(idx_dip_gi + 815);

    auto tr_y_zzzz_xxxxyz = pbuffer.data(idx_dip_gi + 816);

    auto tr_y_zzzz_xxxxzz = pbuffer.data(idx_dip_gi + 817);

    auto tr_y_zzzz_xxxyyy = pbuffer.data(idx_dip_gi + 818);

    auto tr_y_zzzz_xxxyyz = pbuffer.data(idx_dip_gi + 819);

    auto tr_y_zzzz_xxxyzz = pbuffer.data(idx_dip_gi + 820);

    auto tr_y_zzzz_xxxzzz = pbuffer.data(idx_dip_gi + 821);

    auto tr_y_zzzz_xxyyyy = pbuffer.data(idx_dip_gi + 822);

    auto tr_y_zzzz_xxyyyz = pbuffer.data(idx_dip_gi + 823);

    auto tr_y_zzzz_xxyyzz = pbuffer.data(idx_dip_gi + 824);

    auto tr_y_zzzz_xxyzzz = pbuffer.data(idx_dip_gi + 825);

    auto tr_y_zzzz_xxzzzz = pbuffer.data(idx_dip_gi + 826);

    auto tr_y_zzzz_xyyyyy = pbuffer.data(idx_dip_gi + 827);

    auto tr_y_zzzz_xyyyyz = pbuffer.data(idx_dip_gi + 828);

    auto tr_y_zzzz_xyyyzz = pbuffer.data(idx_dip_gi + 829);

    auto tr_y_zzzz_xyyzzz = pbuffer.data(idx_dip_gi + 830);

    auto tr_y_zzzz_xyzzzz = pbuffer.data(idx_dip_gi + 831);

    auto tr_y_zzzz_xzzzzz = pbuffer.data(idx_dip_gi + 832);

    auto tr_y_zzzz_yyyyyy = pbuffer.data(idx_dip_gi + 833);

    auto tr_y_zzzz_yyyyyz = pbuffer.data(idx_dip_gi + 834);

    auto tr_y_zzzz_yyyyzz = pbuffer.data(idx_dip_gi + 835);

    auto tr_y_zzzz_yyyzzz = pbuffer.data(idx_dip_gi + 836);

    auto tr_y_zzzz_yyzzzz = pbuffer.data(idx_dip_gi + 837);

    auto tr_y_zzzz_yzzzzz = pbuffer.data(idx_dip_gi + 838);

    auto tr_y_zzzz_zzzzzz = pbuffer.data(idx_dip_gi + 839);

    auto tr_z_xxxx_xxxxxx = pbuffer.data(idx_dip_gi + 840);

    auto tr_z_xxxx_xxxxxy = pbuffer.data(idx_dip_gi + 841);

    auto tr_z_xxxx_xxxxxz = pbuffer.data(idx_dip_gi + 842);

    auto tr_z_xxxx_xxxxyy = pbuffer.data(idx_dip_gi + 843);

    auto tr_z_xxxx_xxxxyz = pbuffer.data(idx_dip_gi + 844);

    auto tr_z_xxxx_xxxxzz = pbuffer.data(idx_dip_gi + 845);

    auto tr_z_xxxx_xxxyyy = pbuffer.data(idx_dip_gi + 846);

    auto tr_z_xxxx_xxxyyz = pbuffer.data(idx_dip_gi + 847);

    auto tr_z_xxxx_xxxyzz = pbuffer.data(idx_dip_gi + 848);

    auto tr_z_xxxx_xxxzzz = pbuffer.data(idx_dip_gi + 849);

    auto tr_z_xxxx_xxyyyy = pbuffer.data(idx_dip_gi + 850);

    auto tr_z_xxxx_xxyyyz = pbuffer.data(idx_dip_gi + 851);

    auto tr_z_xxxx_xxyyzz = pbuffer.data(idx_dip_gi + 852);

    auto tr_z_xxxx_xxyzzz = pbuffer.data(idx_dip_gi + 853);

    auto tr_z_xxxx_xxzzzz = pbuffer.data(idx_dip_gi + 854);

    auto tr_z_xxxx_xyyyyy = pbuffer.data(idx_dip_gi + 855);

    auto tr_z_xxxx_xyyyyz = pbuffer.data(idx_dip_gi + 856);

    auto tr_z_xxxx_xyyyzz = pbuffer.data(idx_dip_gi + 857);

    auto tr_z_xxxx_xyyzzz = pbuffer.data(idx_dip_gi + 858);

    auto tr_z_xxxx_xyzzzz = pbuffer.data(idx_dip_gi + 859);

    auto tr_z_xxxx_xzzzzz = pbuffer.data(idx_dip_gi + 860);

    auto tr_z_xxxx_yyyyyy = pbuffer.data(idx_dip_gi + 861);

    auto tr_z_xxxx_yyyyyz = pbuffer.data(idx_dip_gi + 862);

    auto tr_z_xxxx_yyyyzz = pbuffer.data(idx_dip_gi + 863);

    auto tr_z_xxxx_yyyzzz = pbuffer.data(idx_dip_gi + 864);

    auto tr_z_xxxx_yyzzzz = pbuffer.data(idx_dip_gi + 865);

    auto tr_z_xxxx_yzzzzz = pbuffer.data(idx_dip_gi + 866);

    auto tr_z_xxxx_zzzzzz = pbuffer.data(idx_dip_gi + 867);

    auto tr_z_xxxy_xxxxxx = pbuffer.data(idx_dip_gi + 868);

    auto tr_z_xxxy_xxxxxz = pbuffer.data(idx_dip_gi + 870);

    auto tr_z_xxxy_xxxxzz = pbuffer.data(idx_dip_gi + 873);

    auto tr_z_xxxy_xxxzzz = pbuffer.data(idx_dip_gi + 877);

    auto tr_z_xxxy_xxzzzz = pbuffer.data(idx_dip_gi + 882);

    auto tr_z_xxxy_xzzzzz = pbuffer.data(idx_dip_gi + 888);

    auto tr_z_xxxy_yyyyyy = pbuffer.data(idx_dip_gi + 889);

    auto tr_z_xxxy_yyyyyz = pbuffer.data(idx_dip_gi + 890);

    auto tr_z_xxxy_yyyyzz = pbuffer.data(idx_dip_gi + 891);

    auto tr_z_xxxy_yyyzzz = pbuffer.data(idx_dip_gi + 892);

    auto tr_z_xxxy_yyzzzz = pbuffer.data(idx_dip_gi + 893);

    auto tr_z_xxxy_yzzzzz = pbuffer.data(idx_dip_gi + 894);

    auto tr_z_xxxz_xxxxxx = pbuffer.data(idx_dip_gi + 896);

    auto tr_z_xxxz_xxxxxy = pbuffer.data(idx_dip_gi + 897);

    auto tr_z_xxxz_xxxxxz = pbuffer.data(idx_dip_gi + 898);

    auto tr_z_xxxz_xxxxyy = pbuffer.data(idx_dip_gi + 899);

    auto tr_z_xxxz_xxxxyz = pbuffer.data(idx_dip_gi + 900);

    auto tr_z_xxxz_xxxxzz = pbuffer.data(idx_dip_gi + 901);

    auto tr_z_xxxz_xxxyyy = pbuffer.data(idx_dip_gi + 902);

    auto tr_z_xxxz_xxxyyz = pbuffer.data(idx_dip_gi + 903);

    auto tr_z_xxxz_xxxyzz = pbuffer.data(idx_dip_gi + 904);

    auto tr_z_xxxz_xxxzzz = pbuffer.data(idx_dip_gi + 905);

    auto tr_z_xxxz_xxyyyy = pbuffer.data(idx_dip_gi + 906);

    auto tr_z_xxxz_xxyyyz = pbuffer.data(idx_dip_gi + 907);

    auto tr_z_xxxz_xxyyzz = pbuffer.data(idx_dip_gi + 908);

    auto tr_z_xxxz_xxyzzz = pbuffer.data(idx_dip_gi + 909);

    auto tr_z_xxxz_xxzzzz = pbuffer.data(idx_dip_gi + 910);

    auto tr_z_xxxz_xyyyyy = pbuffer.data(idx_dip_gi + 911);

    auto tr_z_xxxz_xyyyyz = pbuffer.data(idx_dip_gi + 912);

    auto tr_z_xxxz_xyyyzz = pbuffer.data(idx_dip_gi + 913);

    auto tr_z_xxxz_xyyzzz = pbuffer.data(idx_dip_gi + 914);

    auto tr_z_xxxz_xyzzzz = pbuffer.data(idx_dip_gi + 915);

    auto tr_z_xxxz_xzzzzz = pbuffer.data(idx_dip_gi + 916);

    auto tr_z_xxxz_yyyyyy = pbuffer.data(idx_dip_gi + 917);

    auto tr_z_xxxz_yyyyyz = pbuffer.data(idx_dip_gi + 918);

    auto tr_z_xxxz_yyyyzz = pbuffer.data(idx_dip_gi + 919);

    auto tr_z_xxxz_yyyzzz = pbuffer.data(idx_dip_gi + 920);

    auto tr_z_xxxz_yyzzzz = pbuffer.data(idx_dip_gi + 921);

    auto tr_z_xxxz_yzzzzz = pbuffer.data(idx_dip_gi + 922);

    auto tr_z_xxxz_zzzzzz = pbuffer.data(idx_dip_gi + 923);

    auto tr_z_xxyy_xxxxxx = pbuffer.data(idx_dip_gi + 924);

    auto tr_z_xxyy_xxxxxy = pbuffer.data(idx_dip_gi + 925);

    auto tr_z_xxyy_xxxxxz = pbuffer.data(idx_dip_gi + 926);

    auto tr_z_xxyy_xxxxyy = pbuffer.data(idx_dip_gi + 927);

    auto tr_z_xxyy_xxxxyz = pbuffer.data(idx_dip_gi + 928);

    auto tr_z_xxyy_xxxxzz = pbuffer.data(idx_dip_gi + 929);

    auto tr_z_xxyy_xxxyyy = pbuffer.data(idx_dip_gi + 930);

    auto tr_z_xxyy_xxxyyz = pbuffer.data(idx_dip_gi + 931);

    auto tr_z_xxyy_xxxyzz = pbuffer.data(idx_dip_gi + 932);

    auto tr_z_xxyy_xxxzzz = pbuffer.data(idx_dip_gi + 933);

    auto tr_z_xxyy_xxyyyy = pbuffer.data(idx_dip_gi + 934);

    auto tr_z_xxyy_xxyyyz = pbuffer.data(idx_dip_gi + 935);

    auto tr_z_xxyy_xxyyzz = pbuffer.data(idx_dip_gi + 936);

    auto tr_z_xxyy_xxyzzz = pbuffer.data(idx_dip_gi + 937);

    auto tr_z_xxyy_xxzzzz = pbuffer.data(idx_dip_gi + 938);

    auto tr_z_xxyy_xyyyyy = pbuffer.data(idx_dip_gi + 939);

    auto tr_z_xxyy_xyyyyz = pbuffer.data(idx_dip_gi + 940);

    auto tr_z_xxyy_xyyyzz = pbuffer.data(idx_dip_gi + 941);

    auto tr_z_xxyy_xyyzzz = pbuffer.data(idx_dip_gi + 942);

    auto tr_z_xxyy_xyzzzz = pbuffer.data(idx_dip_gi + 943);

    auto tr_z_xxyy_xzzzzz = pbuffer.data(idx_dip_gi + 944);

    auto tr_z_xxyy_yyyyyy = pbuffer.data(idx_dip_gi + 945);

    auto tr_z_xxyy_yyyyyz = pbuffer.data(idx_dip_gi + 946);

    auto tr_z_xxyy_yyyyzz = pbuffer.data(idx_dip_gi + 947);

    auto tr_z_xxyy_yyyzzz = pbuffer.data(idx_dip_gi + 948);

    auto tr_z_xxyy_yyzzzz = pbuffer.data(idx_dip_gi + 949);

    auto tr_z_xxyy_yzzzzz = pbuffer.data(idx_dip_gi + 950);

    auto tr_z_xxyy_zzzzzz = pbuffer.data(idx_dip_gi + 951);

    auto tr_z_xxyz_xxxxxx = pbuffer.data(idx_dip_gi + 952);

    auto tr_z_xxyz_xxxxxz = pbuffer.data(idx_dip_gi + 954);

    auto tr_z_xxyz_xxxxzz = pbuffer.data(idx_dip_gi + 957);

    auto tr_z_xxyz_xxxzzz = pbuffer.data(idx_dip_gi + 961);

    auto tr_z_xxyz_xxzzzz = pbuffer.data(idx_dip_gi + 966);

    auto tr_z_xxyz_xzzzzz = pbuffer.data(idx_dip_gi + 972);

    auto tr_z_xxyz_yyyyyy = pbuffer.data(idx_dip_gi + 973);

    auto tr_z_xxyz_yyyyyz = pbuffer.data(idx_dip_gi + 974);

    auto tr_z_xxyz_yyyyzz = pbuffer.data(idx_dip_gi + 975);

    auto tr_z_xxyz_yyyzzz = pbuffer.data(idx_dip_gi + 976);

    auto tr_z_xxyz_yyzzzz = pbuffer.data(idx_dip_gi + 977);

    auto tr_z_xxyz_yzzzzz = pbuffer.data(idx_dip_gi + 978);

    auto tr_z_xxzz_xxxxxx = pbuffer.data(idx_dip_gi + 980);

    auto tr_z_xxzz_xxxxxy = pbuffer.data(idx_dip_gi + 981);

    auto tr_z_xxzz_xxxxxz = pbuffer.data(idx_dip_gi + 982);

    auto tr_z_xxzz_xxxxyy = pbuffer.data(idx_dip_gi + 983);

    auto tr_z_xxzz_xxxxyz = pbuffer.data(idx_dip_gi + 984);

    auto tr_z_xxzz_xxxxzz = pbuffer.data(idx_dip_gi + 985);

    auto tr_z_xxzz_xxxyyy = pbuffer.data(idx_dip_gi + 986);

    auto tr_z_xxzz_xxxyyz = pbuffer.data(idx_dip_gi + 987);

    auto tr_z_xxzz_xxxyzz = pbuffer.data(idx_dip_gi + 988);

    auto tr_z_xxzz_xxxzzz = pbuffer.data(idx_dip_gi + 989);

    auto tr_z_xxzz_xxyyyy = pbuffer.data(idx_dip_gi + 990);

    auto tr_z_xxzz_xxyyyz = pbuffer.data(idx_dip_gi + 991);

    auto tr_z_xxzz_xxyyzz = pbuffer.data(idx_dip_gi + 992);

    auto tr_z_xxzz_xxyzzz = pbuffer.data(idx_dip_gi + 993);

    auto tr_z_xxzz_xxzzzz = pbuffer.data(idx_dip_gi + 994);

    auto tr_z_xxzz_xyyyyy = pbuffer.data(idx_dip_gi + 995);

    auto tr_z_xxzz_xyyyyz = pbuffer.data(idx_dip_gi + 996);

    auto tr_z_xxzz_xyyyzz = pbuffer.data(idx_dip_gi + 997);

    auto tr_z_xxzz_xyyzzz = pbuffer.data(idx_dip_gi + 998);

    auto tr_z_xxzz_xyzzzz = pbuffer.data(idx_dip_gi + 999);

    auto tr_z_xxzz_xzzzzz = pbuffer.data(idx_dip_gi + 1000);

    auto tr_z_xxzz_yyyyyy = pbuffer.data(idx_dip_gi + 1001);

    auto tr_z_xxzz_yyyyyz = pbuffer.data(idx_dip_gi + 1002);

    auto tr_z_xxzz_yyyyzz = pbuffer.data(idx_dip_gi + 1003);

    auto tr_z_xxzz_yyyzzz = pbuffer.data(idx_dip_gi + 1004);

    auto tr_z_xxzz_yyzzzz = pbuffer.data(idx_dip_gi + 1005);

    auto tr_z_xxzz_yzzzzz = pbuffer.data(idx_dip_gi + 1006);

    auto tr_z_xxzz_zzzzzz = pbuffer.data(idx_dip_gi + 1007);

    auto tr_z_xyyy_xxxxxy = pbuffer.data(idx_dip_gi + 1009);

    auto tr_z_xyyy_xxxxyy = pbuffer.data(idx_dip_gi + 1011);

    auto tr_z_xyyy_xxxxyz = pbuffer.data(idx_dip_gi + 1012);

    auto tr_z_xyyy_xxxyyy = pbuffer.data(idx_dip_gi + 1014);

    auto tr_z_xyyy_xxxyyz = pbuffer.data(idx_dip_gi + 1015);

    auto tr_z_xyyy_xxxyzz = pbuffer.data(idx_dip_gi + 1016);

    auto tr_z_xyyy_xxyyyy = pbuffer.data(idx_dip_gi + 1018);

    auto tr_z_xyyy_xxyyyz = pbuffer.data(idx_dip_gi + 1019);

    auto tr_z_xyyy_xxyyzz = pbuffer.data(idx_dip_gi + 1020);

    auto tr_z_xyyy_xxyzzz = pbuffer.data(idx_dip_gi + 1021);

    auto tr_z_xyyy_xyyyyy = pbuffer.data(idx_dip_gi + 1023);

    auto tr_z_xyyy_xyyyyz = pbuffer.data(idx_dip_gi + 1024);

    auto tr_z_xyyy_xyyyzz = pbuffer.data(idx_dip_gi + 1025);

    auto tr_z_xyyy_xyyzzz = pbuffer.data(idx_dip_gi + 1026);

    auto tr_z_xyyy_xyzzzz = pbuffer.data(idx_dip_gi + 1027);

    auto tr_z_xyyy_yyyyyy = pbuffer.data(idx_dip_gi + 1029);

    auto tr_z_xyyy_yyyyyz = pbuffer.data(idx_dip_gi + 1030);

    auto tr_z_xyyy_yyyyzz = pbuffer.data(idx_dip_gi + 1031);

    auto tr_z_xyyy_yyyzzz = pbuffer.data(idx_dip_gi + 1032);

    auto tr_z_xyyy_yyzzzz = pbuffer.data(idx_dip_gi + 1033);

    auto tr_z_xyyy_yzzzzz = pbuffer.data(idx_dip_gi + 1034);

    auto tr_z_xyyy_zzzzzz = pbuffer.data(idx_dip_gi + 1035);

    auto tr_z_xyyz_xxxxyz = pbuffer.data(idx_dip_gi + 1040);

    auto tr_z_xyyz_xxxyyz = pbuffer.data(idx_dip_gi + 1043);

    auto tr_z_xyyz_xxxyzz = pbuffer.data(idx_dip_gi + 1044);

    auto tr_z_xyyz_xxyyyz = pbuffer.data(idx_dip_gi + 1047);

    auto tr_z_xyyz_xxyyzz = pbuffer.data(idx_dip_gi + 1048);

    auto tr_z_xyyz_xxyzzz = pbuffer.data(idx_dip_gi + 1049);

    auto tr_z_xyyz_xyyyyz = pbuffer.data(idx_dip_gi + 1052);

    auto tr_z_xyyz_xyyyzz = pbuffer.data(idx_dip_gi + 1053);

    auto tr_z_xyyz_xyyzzz = pbuffer.data(idx_dip_gi + 1054);

    auto tr_z_xyyz_xyzzzz = pbuffer.data(idx_dip_gi + 1055);

    auto tr_z_xyyz_yyyyyy = pbuffer.data(idx_dip_gi + 1057);

    auto tr_z_xyyz_yyyyyz = pbuffer.data(idx_dip_gi + 1058);

    auto tr_z_xyyz_yyyyzz = pbuffer.data(idx_dip_gi + 1059);

    auto tr_z_xyyz_yyyzzz = pbuffer.data(idx_dip_gi + 1060);

    auto tr_z_xyyz_yyzzzz = pbuffer.data(idx_dip_gi + 1061);

    auto tr_z_xyyz_yzzzzz = pbuffer.data(idx_dip_gi + 1062);

    auto tr_z_xyyz_zzzzzz = pbuffer.data(idx_dip_gi + 1063);

    auto tr_z_xyzz_yyyyyy = pbuffer.data(idx_dip_gi + 1085);

    auto tr_z_xyzz_yyyyyz = pbuffer.data(idx_dip_gi + 1086);

    auto tr_z_xyzz_yyyyzz = pbuffer.data(idx_dip_gi + 1087);

    auto tr_z_xyzz_yyyzzz = pbuffer.data(idx_dip_gi + 1088);

    auto tr_z_xyzz_yyzzzz = pbuffer.data(idx_dip_gi + 1089);

    auto tr_z_xyzz_yzzzzz = pbuffer.data(idx_dip_gi + 1090);

    auto tr_z_xzzz_xxxxxx = pbuffer.data(idx_dip_gi + 1092);

    auto tr_z_xzzz_xxxxxy = pbuffer.data(idx_dip_gi + 1093);

    auto tr_z_xzzz_xxxxxz = pbuffer.data(idx_dip_gi + 1094);

    auto tr_z_xzzz_xxxxyy = pbuffer.data(idx_dip_gi + 1095);

    auto tr_z_xzzz_xxxxyz = pbuffer.data(idx_dip_gi + 1096);

    auto tr_z_xzzz_xxxxzz = pbuffer.data(idx_dip_gi + 1097);

    auto tr_z_xzzz_xxxyyy = pbuffer.data(idx_dip_gi + 1098);

    auto tr_z_xzzz_xxxyyz = pbuffer.data(idx_dip_gi + 1099);

    auto tr_z_xzzz_xxxyzz = pbuffer.data(idx_dip_gi + 1100);

    auto tr_z_xzzz_xxxzzz = pbuffer.data(idx_dip_gi + 1101);

    auto tr_z_xzzz_xxyyyy = pbuffer.data(idx_dip_gi + 1102);

    auto tr_z_xzzz_xxyyyz = pbuffer.data(idx_dip_gi + 1103);

    auto tr_z_xzzz_xxyyzz = pbuffer.data(idx_dip_gi + 1104);

    auto tr_z_xzzz_xxyzzz = pbuffer.data(idx_dip_gi + 1105);

    auto tr_z_xzzz_xxzzzz = pbuffer.data(idx_dip_gi + 1106);

    auto tr_z_xzzz_xyyyyy = pbuffer.data(idx_dip_gi + 1107);

    auto tr_z_xzzz_xyyyyz = pbuffer.data(idx_dip_gi + 1108);

    auto tr_z_xzzz_xyyyzz = pbuffer.data(idx_dip_gi + 1109);

    auto tr_z_xzzz_xyyzzz = pbuffer.data(idx_dip_gi + 1110);

    auto tr_z_xzzz_xyzzzz = pbuffer.data(idx_dip_gi + 1111);

    auto tr_z_xzzz_xzzzzz = pbuffer.data(idx_dip_gi + 1112);

    auto tr_z_xzzz_yyyyyy = pbuffer.data(idx_dip_gi + 1113);

    auto tr_z_xzzz_yyyyyz = pbuffer.data(idx_dip_gi + 1114);

    auto tr_z_xzzz_yyyyzz = pbuffer.data(idx_dip_gi + 1115);

    auto tr_z_xzzz_yyyzzz = pbuffer.data(idx_dip_gi + 1116);

    auto tr_z_xzzz_yyzzzz = pbuffer.data(idx_dip_gi + 1117);

    auto tr_z_xzzz_yzzzzz = pbuffer.data(idx_dip_gi + 1118);

    auto tr_z_xzzz_zzzzzz = pbuffer.data(idx_dip_gi + 1119);

    auto tr_z_yyyy_xxxxxx = pbuffer.data(idx_dip_gi + 1120);

    auto tr_z_yyyy_xxxxxy = pbuffer.data(idx_dip_gi + 1121);

    auto tr_z_yyyy_xxxxxz = pbuffer.data(idx_dip_gi + 1122);

    auto tr_z_yyyy_xxxxyy = pbuffer.data(idx_dip_gi + 1123);

    auto tr_z_yyyy_xxxxyz = pbuffer.data(idx_dip_gi + 1124);

    auto tr_z_yyyy_xxxxzz = pbuffer.data(idx_dip_gi + 1125);

    auto tr_z_yyyy_xxxyyy = pbuffer.data(idx_dip_gi + 1126);

    auto tr_z_yyyy_xxxyyz = pbuffer.data(idx_dip_gi + 1127);

    auto tr_z_yyyy_xxxyzz = pbuffer.data(idx_dip_gi + 1128);

    auto tr_z_yyyy_xxxzzz = pbuffer.data(idx_dip_gi + 1129);

    auto tr_z_yyyy_xxyyyy = pbuffer.data(idx_dip_gi + 1130);

    auto tr_z_yyyy_xxyyyz = pbuffer.data(idx_dip_gi + 1131);

    auto tr_z_yyyy_xxyyzz = pbuffer.data(idx_dip_gi + 1132);

    auto tr_z_yyyy_xxyzzz = pbuffer.data(idx_dip_gi + 1133);

    auto tr_z_yyyy_xxzzzz = pbuffer.data(idx_dip_gi + 1134);

    auto tr_z_yyyy_xyyyyy = pbuffer.data(idx_dip_gi + 1135);

    auto tr_z_yyyy_xyyyyz = pbuffer.data(idx_dip_gi + 1136);

    auto tr_z_yyyy_xyyyzz = pbuffer.data(idx_dip_gi + 1137);

    auto tr_z_yyyy_xyyzzz = pbuffer.data(idx_dip_gi + 1138);

    auto tr_z_yyyy_xyzzzz = pbuffer.data(idx_dip_gi + 1139);

    auto tr_z_yyyy_xzzzzz = pbuffer.data(idx_dip_gi + 1140);

    auto tr_z_yyyy_yyyyyy = pbuffer.data(idx_dip_gi + 1141);

    auto tr_z_yyyy_yyyyyz = pbuffer.data(idx_dip_gi + 1142);

    auto tr_z_yyyy_yyyyzz = pbuffer.data(idx_dip_gi + 1143);

    auto tr_z_yyyy_yyyzzz = pbuffer.data(idx_dip_gi + 1144);

    auto tr_z_yyyy_yyzzzz = pbuffer.data(idx_dip_gi + 1145);

    auto tr_z_yyyy_yzzzzz = pbuffer.data(idx_dip_gi + 1146);

    auto tr_z_yyyy_zzzzzz = pbuffer.data(idx_dip_gi + 1147);

    auto tr_z_yyyz_xxxxxx = pbuffer.data(idx_dip_gi + 1148);

    auto tr_z_yyyz_xxxxxy = pbuffer.data(idx_dip_gi + 1149);

    auto tr_z_yyyz_xxxxxz = pbuffer.data(idx_dip_gi + 1150);

    auto tr_z_yyyz_xxxxyy = pbuffer.data(idx_dip_gi + 1151);

    auto tr_z_yyyz_xxxxyz = pbuffer.data(idx_dip_gi + 1152);

    auto tr_z_yyyz_xxxxzz = pbuffer.data(idx_dip_gi + 1153);

    auto tr_z_yyyz_xxxyyy = pbuffer.data(idx_dip_gi + 1154);

    auto tr_z_yyyz_xxxyyz = pbuffer.data(idx_dip_gi + 1155);

    auto tr_z_yyyz_xxxyzz = pbuffer.data(idx_dip_gi + 1156);

    auto tr_z_yyyz_xxxzzz = pbuffer.data(idx_dip_gi + 1157);

    auto tr_z_yyyz_xxyyyy = pbuffer.data(idx_dip_gi + 1158);

    auto tr_z_yyyz_xxyyyz = pbuffer.data(idx_dip_gi + 1159);

    auto tr_z_yyyz_xxyyzz = pbuffer.data(idx_dip_gi + 1160);

    auto tr_z_yyyz_xxyzzz = pbuffer.data(idx_dip_gi + 1161);

    auto tr_z_yyyz_xxzzzz = pbuffer.data(idx_dip_gi + 1162);

    auto tr_z_yyyz_xyyyyy = pbuffer.data(idx_dip_gi + 1163);

    auto tr_z_yyyz_xyyyyz = pbuffer.data(idx_dip_gi + 1164);

    auto tr_z_yyyz_xyyyzz = pbuffer.data(idx_dip_gi + 1165);

    auto tr_z_yyyz_xyyzzz = pbuffer.data(idx_dip_gi + 1166);

    auto tr_z_yyyz_xyzzzz = pbuffer.data(idx_dip_gi + 1167);

    auto tr_z_yyyz_xzzzzz = pbuffer.data(idx_dip_gi + 1168);

    auto tr_z_yyyz_yyyyyy = pbuffer.data(idx_dip_gi + 1169);

    auto tr_z_yyyz_yyyyyz = pbuffer.data(idx_dip_gi + 1170);

    auto tr_z_yyyz_yyyyzz = pbuffer.data(idx_dip_gi + 1171);

    auto tr_z_yyyz_yyyzzz = pbuffer.data(idx_dip_gi + 1172);

    auto tr_z_yyyz_yyzzzz = pbuffer.data(idx_dip_gi + 1173);

    auto tr_z_yyyz_yzzzzz = pbuffer.data(idx_dip_gi + 1174);

    auto tr_z_yyyz_zzzzzz = pbuffer.data(idx_dip_gi + 1175);

    auto tr_z_yyzz_xxxxxx = pbuffer.data(idx_dip_gi + 1176);

    auto tr_z_yyzz_xxxxxy = pbuffer.data(idx_dip_gi + 1177);

    auto tr_z_yyzz_xxxxxz = pbuffer.data(idx_dip_gi + 1178);

    auto tr_z_yyzz_xxxxyy = pbuffer.data(idx_dip_gi + 1179);

    auto tr_z_yyzz_xxxxyz = pbuffer.data(idx_dip_gi + 1180);

    auto tr_z_yyzz_xxxxzz = pbuffer.data(idx_dip_gi + 1181);

    auto tr_z_yyzz_xxxyyy = pbuffer.data(idx_dip_gi + 1182);

    auto tr_z_yyzz_xxxyyz = pbuffer.data(idx_dip_gi + 1183);

    auto tr_z_yyzz_xxxyzz = pbuffer.data(idx_dip_gi + 1184);

    auto tr_z_yyzz_xxxzzz = pbuffer.data(idx_dip_gi + 1185);

    auto tr_z_yyzz_xxyyyy = pbuffer.data(idx_dip_gi + 1186);

    auto tr_z_yyzz_xxyyyz = pbuffer.data(idx_dip_gi + 1187);

    auto tr_z_yyzz_xxyyzz = pbuffer.data(idx_dip_gi + 1188);

    auto tr_z_yyzz_xxyzzz = pbuffer.data(idx_dip_gi + 1189);

    auto tr_z_yyzz_xxzzzz = pbuffer.data(idx_dip_gi + 1190);

    auto tr_z_yyzz_xyyyyy = pbuffer.data(idx_dip_gi + 1191);

    auto tr_z_yyzz_xyyyyz = pbuffer.data(idx_dip_gi + 1192);

    auto tr_z_yyzz_xyyyzz = pbuffer.data(idx_dip_gi + 1193);

    auto tr_z_yyzz_xyyzzz = pbuffer.data(idx_dip_gi + 1194);

    auto tr_z_yyzz_xyzzzz = pbuffer.data(idx_dip_gi + 1195);

    auto tr_z_yyzz_xzzzzz = pbuffer.data(idx_dip_gi + 1196);

    auto tr_z_yyzz_yyyyyy = pbuffer.data(idx_dip_gi + 1197);

    auto tr_z_yyzz_yyyyyz = pbuffer.data(idx_dip_gi + 1198);

    auto tr_z_yyzz_yyyyzz = pbuffer.data(idx_dip_gi + 1199);

    auto tr_z_yyzz_yyyzzz = pbuffer.data(idx_dip_gi + 1200);

    auto tr_z_yyzz_yyzzzz = pbuffer.data(idx_dip_gi + 1201);

    auto tr_z_yyzz_yzzzzz = pbuffer.data(idx_dip_gi + 1202);

    auto tr_z_yyzz_zzzzzz = pbuffer.data(idx_dip_gi + 1203);

    auto tr_z_yzzz_xxxxxx = pbuffer.data(idx_dip_gi + 1204);

    auto tr_z_yzzz_xxxxxy = pbuffer.data(idx_dip_gi + 1205);

    auto tr_z_yzzz_xxxxxz = pbuffer.data(idx_dip_gi + 1206);

    auto tr_z_yzzz_xxxxyy = pbuffer.data(idx_dip_gi + 1207);

    auto tr_z_yzzz_xxxxyz = pbuffer.data(idx_dip_gi + 1208);

    auto tr_z_yzzz_xxxxzz = pbuffer.data(idx_dip_gi + 1209);

    auto tr_z_yzzz_xxxyyy = pbuffer.data(idx_dip_gi + 1210);

    auto tr_z_yzzz_xxxyyz = pbuffer.data(idx_dip_gi + 1211);

    auto tr_z_yzzz_xxxyzz = pbuffer.data(idx_dip_gi + 1212);

    auto tr_z_yzzz_xxxzzz = pbuffer.data(idx_dip_gi + 1213);

    auto tr_z_yzzz_xxyyyy = pbuffer.data(idx_dip_gi + 1214);

    auto tr_z_yzzz_xxyyyz = pbuffer.data(idx_dip_gi + 1215);

    auto tr_z_yzzz_xxyyzz = pbuffer.data(idx_dip_gi + 1216);

    auto tr_z_yzzz_xxyzzz = pbuffer.data(idx_dip_gi + 1217);

    auto tr_z_yzzz_xxzzzz = pbuffer.data(idx_dip_gi + 1218);

    auto tr_z_yzzz_xyyyyy = pbuffer.data(idx_dip_gi + 1219);

    auto tr_z_yzzz_xyyyyz = pbuffer.data(idx_dip_gi + 1220);

    auto tr_z_yzzz_xyyyzz = pbuffer.data(idx_dip_gi + 1221);

    auto tr_z_yzzz_xyyzzz = pbuffer.data(idx_dip_gi + 1222);

    auto tr_z_yzzz_xyzzzz = pbuffer.data(idx_dip_gi + 1223);

    auto tr_z_yzzz_xzzzzz = pbuffer.data(idx_dip_gi + 1224);

    auto tr_z_yzzz_yyyyyy = pbuffer.data(idx_dip_gi + 1225);

    auto tr_z_yzzz_yyyyyz = pbuffer.data(idx_dip_gi + 1226);

    auto tr_z_yzzz_yyyyzz = pbuffer.data(idx_dip_gi + 1227);

    auto tr_z_yzzz_yyyzzz = pbuffer.data(idx_dip_gi + 1228);

    auto tr_z_yzzz_yyzzzz = pbuffer.data(idx_dip_gi + 1229);

    auto tr_z_yzzz_yzzzzz = pbuffer.data(idx_dip_gi + 1230);

    auto tr_z_yzzz_zzzzzz = pbuffer.data(idx_dip_gi + 1231);

    auto tr_z_zzzz_xxxxxx = pbuffer.data(idx_dip_gi + 1232);

    auto tr_z_zzzz_xxxxxy = pbuffer.data(idx_dip_gi + 1233);

    auto tr_z_zzzz_xxxxxz = pbuffer.data(idx_dip_gi + 1234);

    auto tr_z_zzzz_xxxxyy = pbuffer.data(idx_dip_gi + 1235);

    auto tr_z_zzzz_xxxxyz = pbuffer.data(idx_dip_gi + 1236);

    auto tr_z_zzzz_xxxxzz = pbuffer.data(idx_dip_gi + 1237);

    auto tr_z_zzzz_xxxyyy = pbuffer.data(idx_dip_gi + 1238);

    auto tr_z_zzzz_xxxyyz = pbuffer.data(idx_dip_gi + 1239);

    auto tr_z_zzzz_xxxyzz = pbuffer.data(idx_dip_gi + 1240);

    auto tr_z_zzzz_xxxzzz = pbuffer.data(idx_dip_gi + 1241);

    auto tr_z_zzzz_xxyyyy = pbuffer.data(idx_dip_gi + 1242);

    auto tr_z_zzzz_xxyyyz = pbuffer.data(idx_dip_gi + 1243);

    auto tr_z_zzzz_xxyyzz = pbuffer.data(idx_dip_gi + 1244);

    auto tr_z_zzzz_xxyzzz = pbuffer.data(idx_dip_gi + 1245);

    auto tr_z_zzzz_xxzzzz = pbuffer.data(idx_dip_gi + 1246);

    auto tr_z_zzzz_xyyyyy = pbuffer.data(idx_dip_gi + 1247);

    auto tr_z_zzzz_xyyyyz = pbuffer.data(idx_dip_gi + 1248);

    auto tr_z_zzzz_xyyyzz = pbuffer.data(idx_dip_gi + 1249);

    auto tr_z_zzzz_xyyzzz = pbuffer.data(idx_dip_gi + 1250);

    auto tr_z_zzzz_xyzzzz = pbuffer.data(idx_dip_gi + 1251);

    auto tr_z_zzzz_xzzzzz = pbuffer.data(idx_dip_gi + 1252);

    auto tr_z_zzzz_yyyyyy = pbuffer.data(idx_dip_gi + 1253);

    auto tr_z_zzzz_yyyyyz = pbuffer.data(idx_dip_gi + 1254);

    auto tr_z_zzzz_yyyyzz = pbuffer.data(idx_dip_gi + 1255);

    auto tr_z_zzzz_yyyzzz = pbuffer.data(idx_dip_gi + 1256);

    auto tr_z_zzzz_yyzzzz = pbuffer.data(idx_dip_gi + 1257);

    auto tr_z_zzzz_yzzzzz = pbuffer.data(idx_dip_gi + 1258);

    auto tr_z_zzzz_zzzzzz = pbuffer.data(idx_dip_gi + 1259);

    // Set up 0-28 components of targeted buffer : HI

    auto tr_x_xxxxx_xxxxxx = pbuffer.data(idx_dip_hi);

    auto tr_x_xxxxx_xxxxxy = pbuffer.data(idx_dip_hi + 1);

    auto tr_x_xxxxx_xxxxxz = pbuffer.data(idx_dip_hi + 2);

    auto tr_x_xxxxx_xxxxyy = pbuffer.data(idx_dip_hi + 3);

    auto tr_x_xxxxx_xxxxyz = pbuffer.data(idx_dip_hi + 4);

    auto tr_x_xxxxx_xxxxzz = pbuffer.data(idx_dip_hi + 5);

    auto tr_x_xxxxx_xxxyyy = pbuffer.data(idx_dip_hi + 6);

    auto tr_x_xxxxx_xxxyyz = pbuffer.data(idx_dip_hi + 7);

    auto tr_x_xxxxx_xxxyzz = pbuffer.data(idx_dip_hi + 8);

    auto tr_x_xxxxx_xxxzzz = pbuffer.data(idx_dip_hi + 9);

    auto tr_x_xxxxx_xxyyyy = pbuffer.data(idx_dip_hi + 10);

    auto tr_x_xxxxx_xxyyyz = pbuffer.data(idx_dip_hi + 11);

    auto tr_x_xxxxx_xxyyzz = pbuffer.data(idx_dip_hi + 12);

    auto tr_x_xxxxx_xxyzzz = pbuffer.data(idx_dip_hi + 13);

    auto tr_x_xxxxx_xxzzzz = pbuffer.data(idx_dip_hi + 14);

    auto tr_x_xxxxx_xyyyyy = pbuffer.data(idx_dip_hi + 15);

    auto tr_x_xxxxx_xyyyyz = pbuffer.data(idx_dip_hi + 16);

    auto tr_x_xxxxx_xyyyzz = pbuffer.data(idx_dip_hi + 17);

    auto tr_x_xxxxx_xyyzzz = pbuffer.data(idx_dip_hi + 18);

    auto tr_x_xxxxx_xyzzzz = pbuffer.data(idx_dip_hi + 19);

    auto tr_x_xxxxx_xzzzzz = pbuffer.data(idx_dip_hi + 20);

    auto tr_x_xxxxx_yyyyyy = pbuffer.data(idx_dip_hi + 21);

    auto tr_x_xxxxx_yyyyyz = pbuffer.data(idx_dip_hi + 22);

    auto tr_x_xxxxx_yyyyzz = pbuffer.data(idx_dip_hi + 23);

    auto tr_x_xxxxx_yyyzzz = pbuffer.data(idx_dip_hi + 24);

    auto tr_x_xxxxx_yyzzzz = pbuffer.data(idx_dip_hi + 25);

    auto tr_x_xxxxx_yzzzzz = pbuffer.data(idx_dip_hi + 26);

    auto tr_x_xxxxx_zzzzzz = pbuffer.data(idx_dip_hi + 27);

#pragma omp simd aligned(pa_x,                  \
                             tr_x_xxx_xxxxxx,   \
                             tr_x_xxx_xxxxxy,   \
                             tr_x_xxx_xxxxxz,   \
                             tr_x_xxx_xxxxyy,   \
                             tr_x_xxx_xxxxyz,   \
                             tr_x_xxx_xxxxzz,   \
                             tr_x_xxx_xxxyyy,   \
                             tr_x_xxx_xxxyyz,   \
                             tr_x_xxx_xxxyzz,   \
                             tr_x_xxx_xxxzzz,   \
                             tr_x_xxx_xxyyyy,   \
                             tr_x_xxx_xxyyyz,   \
                             tr_x_xxx_xxyyzz,   \
                             tr_x_xxx_xxyzzz,   \
                             tr_x_xxx_xxzzzz,   \
                             tr_x_xxx_xyyyyy,   \
                             tr_x_xxx_xyyyyz,   \
                             tr_x_xxx_xyyyzz,   \
                             tr_x_xxx_xyyzzz,   \
                             tr_x_xxx_xyzzzz,   \
                             tr_x_xxx_xzzzzz,   \
                             tr_x_xxx_yyyyyy,   \
                             tr_x_xxx_yyyyyz,   \
                             tr_x_xxx_yyyyzz,   \
                             tr_x_xxx_yyyzzz,   \
                             tr_x_xxx_yyzzzz,   \
                             tr_x_xxx_yzzzzz,   \
                             tr_x_xxx_zzzzzz,   \
                             tr_x_xxxx_xxxxx,   \
                             tr_x_xxxx_xxxxxx,  \
                             tr_x_xxxx_xxxxxy,  \
                             tr_x_xxxx_xxxxxz,  \
                             tr_x_xxxx_xxxxy,   \
                             tr_x_xxxx_xxxxyy,  \
                             tr_x_xxxx_xxxxyz,  \
                             tr_x_xxxx_xxxxz,   \
                             tr_x_xxxx_xxxxzz,  \
                             tr_x_xxxx_xxxyy,   \
                             tr_x_xxxx_xxxyyy,  \
                             tr_x_xxxx_xxxyyz,  \
                             tr_x_xxxx_xxxyz,   \
                             tr_x_xxxx_xxxyzz,  \
                             tr_x_xxxx_xxxzz,   \
                             tr_x_xxxx_xxxzzz,  \
                             tr_x_xxxx_xxyyy,   \
                             tr_x_xxxx_xxyyyy,  \
                             tr_x_xxxx_xxyyyz,  \
                             tr_x_xxxx_xxyyz,   \
                             tr_x_xxxx_xxyyzz,  \
                             tr_x_xxxx_xxyzz,   \
                             tr_x_xxxx_xxyzzz,  \
                             tr_x_xxxx_xxzzz,   \
                             tr_x_xxxx_xxzzzz,  \
                             tr_x_xxxx_xyyyy,   \
                             tr_x_xxxx_xyyyyy,  \
                             tr_x_xxxx_xyyyyz,  \
                             tr_x_xxxx_xyyyz,   \
                             tr_x_xxxx_xyyyzz,  \
                             tr_x_xxxx_xyyzz,   \
                             tr_x_xxxx_xyyzzz,  \
                             tr_x_xxxx_xyzzz,   \
                             tr_x_xxxx_xyzzzz,  \
                             tr_x_xxxx_xzzzz,   \
                             tr_x_xxxx_xzzzzz,  \
                             tr_x_xxxx_yyyyy,   \
                             tr_x_xxxx_yyyyyy,  \
                             tr_x_xxxx_yyyyyz,  \
                             tr_x_xxxx_yyyyz,   \
                             tr_x_xxxx_yyyyzz,  \
                             tr_x_xxxx_yyyzz,   \
                             tr_x_xxxx_yyyzzz,  \
                             tr_x_xxxx_yyzzz,   \
                             tr_x_xxxx_yyzzzz,  \
                             tr_x_xxxx_yzzzz,   \
                             tr_x_xxxx_yzzzzz,  \
                             tr_x_xxxx_zzzzz,   \
                             tr_x_xxxx_zzzzzz,  \
                             tr_x_xxxxx_xxxxxx, \
                             tr_x_xxxxx_xxxxxy, \
                             tr_x_xxxxx_xxxxxz, \
                             tr_x_xxxxx_xxxxyy, \
                             tr_x_xxxxx_xxxxyz, \
                             tr_x_xxxxx_xxxxzz, \
                             tr_x_xxxxx_xxxyyy, \
                             tr_x_xxxxx_xxxyyz, \
                             tr_x_xxxxx_xxxyzz, \
                             tr_x_xxxxx_xxxzzz, \
                             tr_x_xxxxx_xxyyyy, \
                             tr_x_xxxxx_xxyyyz, \
                             tr_x_xxxxx_xxyyzz, \
                             tr_x_xxxxx_xxyzzz, \
                             tr_x_xxxxx_xxzzzz, \
                             tr_x_xxxxx_xyyyyy, \
                             tr_x_xxxxx_xyyyyz, \
                             tr_x_xxxxx_xyyyzz, \
                             tr_x_xxxxx_xyyzzz, \
                             tr_x_xxxxx_xyzzzz, \
                             tr_x_xxxxx_xzzzzz, \
                             tr_x_xxxxx_yyyyyy, \
                             tr_x_xxxxx_yyyyyz, \
                             tr_x_xxxxx_yyyyzz, \
                             tr_x_xxxxx_yyyzzz, \
                             tr_x_xxxxx_yyzzzz, \
                             tr_x_xxxxx_yzzzzz, \
                             tr_x_xxxxx_zzzzzz, \
                             ts_xxxx_xxxxxx,    \
                             ts_xxxx_xxxxxy,    \
                             ts_xxxx_xxxxxz,    \
                             ts_xxxx_xxxxyy,    \
                             ts_xxxx_xxxxyz,    \
                             ts_xxxx_xxxxzz,    \
                             ts_xxxx_xxxyyy,    \
                             ts_xxxx_xxxyyz,    \
                             ts_xxxx_xxxyzz,    \
                             ts_xxxx_xxxzzz,    \
                             ts_xxxx_xxyyyy,    \
                             ts_xxxx_xxyyyz,    \
                             ts_xxxx_xxyyzz,    \
                             ts_xxxx_xxyzzz,    \
                             ts_xxxx_xxzzzz,    \
                             ts_xxxx_xyyyyy,    \
                             ts_xxxx_xyyyyz,    \
                             ts_xxxx_xyyyzz,    \
                             ts_xxxx_xyyzzz,    \
                             ts_xxxx_xyzzzz,    \
                             ts_xxxx_xzzzzz,    \
                             ts_xxxx_yyyyyy,    \
                             ts_xxxx_yyyyyz,    \
                             ts_xxxx_yyyyzz,    \
                             ts_xxxx_yyyzzz,    \
                             ts_xxxx_yyzzzz,    \
                             ts_xxxx_yzzzzz,    \
                             ts_xxxx_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxx_xxxxxx[i] =
            4.0 * tr_x_xxx_xxxxxx[i] * fe_0 + 6.0 * tr_x_xxxx_xxxxx[i] * fe_0 + ts_xxxx_xxxxxx[i] * fe_0 + tr_x_xxxx_xxxxxx[i] * pa_x[i];

        tr_x_xxxxx_xxxxxy[i] =
            4.0 * tr_x_xxx_xxxxxy[i] * fe_0 + 5.0 * tr_x_xxxx_xxxxy[i] * fe_0 + ts_xxxx_xxxxxy[i] * fe_0 + tr_x_xxxx_xxxxxy[i] * pa_x[i];

        tr_x_xxxxx_xxxxxz[i] =
            4.0 * tr_x_xxx_xxxxxz[i] * fe_0 + 5.0 * tr_x_xxxx_xxxxz[i] * fe_0 + ts_xxxx_xxxxxz[i] * fe_0 + tr_x_xxxx_xxxxxz[i] * pa_x[i];

        tr_x_xxxxx_xxxxyy[i] =
            4.0 * tr_x_xxx_xxxxyy[i] * fe_0 + 4.0 * tr_x_xxxx_xxxyy[i] * fe_0 + ts_xxxx_xxxxyy[i] * fe_0 + tr_x_xxxx_xxxxyy[i] * pa_x[i];

        tr_x_xxxxx_xxxxyz[i] =
            4.0 * tr_x_xxx_xxxxyz[i] * fe_0 + 4.0 * tr_x_xxxx_xxxyz[i] * fe_0 + ts_xxxx_xxxxyz[i] * fe_0 + tr_x_xxxx_xxxxyz[i] * pa_x[i];

        tr_x_xxxxx_xxxxzz[i] =
            4.0 * tr_x_xxx_xxxxzz[i] * fe_0 + 4.0 * tr_x_xxxx_xxxzz[i] * fe_0 + ts_xxxx_xxxxzz[i] * fe_0 + tr_x_xxxx_xxxxzz[i] * pa_x[i];

        tr_x_xxxxx_xxxyyy[i] =
            4.0 * tr_x_xxx_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxxx_xxyyy[i] * fe_0 + ts_xxxx_xxxyyy[i] * fe_0 + tr_x_xxxx_xxxyyy[i] * pa_x[i];

        tr_x_xxxxx_xxxyyz[i] =
            4.0 * tr_x_xxx_xxxyyz[i] * fe_0 + 3.0 * tr_x_xxxx_xxyyz[i] * fe_0 + ts_xxxx_xxxyyz[i] * fe_0 + tr_x_xxxx_xxxyyz[i] * pa_x[i];

        tr_x_xxxxx_xxxyzz[i] =
            4.0 * tr_x_xxx_xxxyzz[i] * fe_0 + 3.0 * tr_x_xxxx_xxyzz[i] * fe_0 + ts_xxxx_xxxyzz[i] * fe_0 + tr_x_xxxx_xxxyzz[i] * pa_x[i];

        tr_x_xxxxx_xxxzzz[i] =
            4.0 * tr_x_xxx_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxxx_xxzzz[i] * fe_0 + ts_xxxx_xxxzzz[i] * fe_0 + tr_x_xxxx_xxxzzz[i] * pa_x[i];

        tr_x_xxxxx_xxyyyy[i] =
            4.0 * tr_x_xxx_xxyyyy[i] * fe_0 + 2.0 * tr_x_xxxx_xyyyy[i] * fe_0 + ts_xxxx_xxyyyy[i] * fe_0 + tr_x_xxxx_xxyyyy[i] * pa_x[i];

        tr_x_xxxxx_xxyyyz[i] =
            4.0 * tr_x_xxx_xxyyyz[i] * fe_0 + 2.0 * tr_x_xxxx_xyyyz[i] * fe_0 + ts_xxxx_xxyyyz[i] * fe_0 + tr_x_xxxx_xxyyyz[i] * pa_x[i];

        tr_x_xxxxx_xxyyzz[i] =
            4.0 * tr_x_xxx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxx_xyyzz[i] * fe_0 + ts_xxxx_xxyyzz[i] * fe_0 + tr_x_xxxx_xxyyzz[i] * pa_x[i];

        tr_x_xxxxx_xxyzzz[i] =
            4.0 * tr_x_xxx_xxyzzz[i] * fe_0 + 2.0 * tr_x_xxxx_xyzzz[i] * fe_0 + ts_xxxx_xxyzzz[i] * fe_0 + tr_x_xxxx_xxyzzz[i] * pa_x[i];

        tr_x_xxxxx_xxzzzz[i] =
            4.0 * tr_x_xxx_xxzzzz[i] * fe_0 + 2.0 * tr_x_xxxx_xzzzz[i] * fe_0 + ts_xxxx_xxzzzz[i] * fe_0 + tr_x_xxxx_xxzzzz[i] * pa_x[i];

        tr_x_xxxxx_xyyyyy[i] = 4.0 * tr_x_xxx_xyyyyy[i] * fe_0 + tr_x_xxxx_yyyyy[i] * fe_0 + ts_xxxx_xyyyyy[i] * fe_0 + tr_x_xxxx_xyyyyy[i] * pa_x[i];

        tr_x_xxxxx_xyyyyz[i] = 4.0 * tr_x_xxx_xyyyyz[i] * fe_0 + tr_x_xxxx_yyyyz[i] * fe_0 + ts_xxxx_xyyyyz[i] * fe_0 + tr_x_xxxx_xyyyyz[i] * pa_x[i];

        tr_x_xxxxx_xyyyzz[i] = 4.0 * tr_x_xxx_xyyyzz[i] * fe_0 + tr_x_xxxx_yyyzz[i] * fe_0 + ts_xxxx_xyyyzz[i] * fe_0 + tr_x_xxxx_xyyyzz[i] * pa_x[i];

        tr_x_xxxxx_xyyzzz[i] = 4.0 * tr_x_xxx_xyyzzz[i] * fe_0 + tr_x_xxxx_yyzzz[i] * fe_0 + ts_xxxx_xyyzzz[i] * fe_0 + tr_x_xxxx_xyyzzz[i] * pa_x[i];

        tr_x_xxxxx_xyzzzz[i] = 4.0 * tr_x_xxx_xyzzzz[i] * fe_0 + tr_x_xxxx_yzzzz[i] * fe_0 + ts_xxxx_xyzzzz[i] * fe_0 + tr_x_xxxx_xyzzzz[i] * pa_x[i];

        tr_x_xxxxx_xzzzzz[i] = 4.0 * tr_x_xxx_xzzzzz[i] * fe_0 + tr_x_xxxx_zzzzz[i] * fe_0 + ts_xxxx_xzzzzz[i] * fe_0 + tr_x_xxxx_xzzzzz[i] * pa_x[i];

        tr_x_xxxxx_yyyyyy[i] = 4.0 * tr_x_xxx_yyyyyy[i] * fe_0 + ts_xxxx_yyyyyy[i] * fe_0 + tr_x_xxxx_yyyyyy[i] * pa_x[i];

        tr_x_xxxxx_yyyyyz[i] = 4.0 * tr_x_xxx_yyyyyz[i] * fe_0 + ts_xxxx_yyyyyz[i] * fe_0 + tr_x_xxxx_yyyyyz[i] * pa_x[i];

        tr_x_xxxxx_yyyyzz[i] = 4.0 * tr_x_xxx_yyyyzz[i] * fe_0 + ts_xxxx_yyyyzz[i] * fe_0 + tr_x_xxxx_yyyyzz[i] * pa_x[i];

        tr_x_xxxxx_yyyzzz[i] = 4.0 * tr_x_xxx_yyyzzz[i] * fe_0 + ts_xxxx_yyyzzz[i] * fe_0 + tr_x_xxxx_yyyzzz[i] * pa_x[i];

        tr_x_xxxxx_yyzzzz[i] = 4.0 * tr_x_xxx_yyzzzz[i] * fe_0 + ts_xxxx_yyzzzz[i] * fe_0 + tr_x_xxxx_yyzzzz[i] * pa_x[i];

        tr_x_xxxxx_yzzzzz[i] = 4.0 * tr_x_xxx_yzzzzz[i] * fe_0 + ts_xxxx_yzzzzz[i] * fe_0 + tr_x_xxxx_yzzzzz[i] * pa_x[i];

        tr_x_xxxxx_zzzzzz[i] = 4.0 * tr_x_xxx_zzzzzz[i] * fe_0 + ts_xxxx_zzzzzz[i] * fe_0 + tr_x_xxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : HI

    auto tr_x_xxxxy_xxxxxx = pbuffer.data(idx_dip_hi + 28);

    auto tr_x_xxxxy_xxxxxy = pbuffer.data(idx_dip_hi + 29);

    auto tr_x_xxxxy_xxxxxz = pbuffer.data(idx_dip_hi + 30);

    auto tr_x_xxxxy_xxxxyy = pbuffer.data(idx_dip_hi + 31);

    auto tr_x_xxxxy_xxxxyz = pbuffer.data(idx_dip_hi + 32);

    auto tr_x_xxxxy_xxxxzz = pbuffer.data(idx_dip_hi + 33);

    auto tr_x_xxxxy_xxxyyy = pbuffer.data(idx_dip_hi + 34);

    auto tr_x_xxxxy_xxxyyz = pbuffer.data(idx_dip_hi + 35);

    auto tr_x_xxxxy_xxxyzz = pbuffer.data(idx_dip_hi + 36);

    auto tr_x_xxxxy_xxxzzz = pbuffer.data(idx_dip_hi + 37);

    auto tr_x_xxxxy_xxyyyy = pbuffer.data(idx_dip_hi + 38);

    auto tr_x_xxxxy_xxyyyz = pbuffer.data(idx_dip_hi + 39);

    auto tr_x_xxxxy_xxyyzz = pbuffer.data(idx_dip_hi + 40);

    auto tr_x_xxxxy_xxyzzz = pbuffer.data(idx_dip_hi + 41);

    auto tr_x_xxxxy_xxzzzz = pbuffer.data(idx_dip_hi + 42);

    auto tr_x_xxxxy_xyyyyy = pbuffer.data(idx_dip_hi + 43);

    auto tr_x_xxxxy_xyyyyz = pbuffer.data(idx_dip_hi + 44);

    auto tr_x_xxxxy_xyyyzz = pbuffer.data(idx_dip_hi + 45);

    auto tr_x_xxxxy_xyyzzz = pbuffer.data(idx_dip_hi + 46);

    auto tr_x_xxxxy_xyzzzz = pbuffer.data(idx_dip_hi + 47);

    auto tr_x_xxxxy_xzzzzz = pbuffer.data(idx_dip_hi + 48);

    auto tr_x_xxxxy_yyyyyy = pbuffer.data(idx_dip_hi + 49);

    auto tr_x_xxxxy_yyyyyz = pbuffer.data(idx_dip_hi + 50);

    auto tr_x_xxxxy_yyyyzz = pbuffer.data(idx_dip_hi + 51);

    auto tr_x_xxxxy_yyyzzz = pbuffer.data(idx_dip_hi + 52);

    auto tr_x_xxxxy_yyzzzz = pbuffer.data(idx_dip_hi + 53);

    auto tr_x_xxxxy_yzzzzz = pbuffer.data(idx_dip_hi + 54);

    auto tr_x_xxxxy_zzzzzz = pbuffer.data(idx_dip_hi + 55);

#pragma omp simd aligned(pa_y,                  \
                             tr_x_xxxx_xxxxx,   \
                             tr_x_xxxx_xxxxxx,  \
                             tr_x_xxxx_xxxxxy,  \
                             tr_x_xxxx_xxxxxz,  \
                             tr_x_xxxx_xxxxy,   \
                             tr_x_xxxx_xxxxyy,  \
                             tr_x_xxxx_xxxxyz,  \
                             tr_x_xxxx_xxxxz,   \
                             tr_x_xxxx_xxxxzz,  \
                             tr_x_xxxx_xxxyy,   \
                             tr_x_xxxx_xxxyyy,  \
                             tr_x_xxxx_xxxyyz,  \
                             tr_x_xxxx_xxxyz,   \
                             tr_x_xxxx_xxxyzz,  \
                             tr_x_xxxx_xxxzz,   \
                             tr_x_xxxx_xxxzzz,  \
                             tr_x_xxxx_xxyyy,   \
                             tr_x_xxxx_xxyyyy,  \
                             tr_x_xxxx_xxyyyz,  \
                             tr_x_xxxx_xxyyz,   \
                             tr_x_xxxx_xxyyzz,  \
                             tr_x_xxxx_xxyzz,   \
                             tr_x_xxxx_xxyzzz,  \
                             tr_x_xxxx_xxzzz,   \
                             tr_x_xxxx_xxzzzz,  \
                             tr_x_xxxx_xyyyy,   \
                             tr_x_xxxx_xyyyyy,  \
                             tr_x_xxxx_xyyyyz,  \
                             tr_x_xxxx_xyyyz,   \
                             tr_x_xxxx_xyyyzz,  \
                             tr_x_xxxx_xyyzz,   \
                             tr_x_xxxx_xyyzzz,  \
                             tr_x_xxxx_xyzzz,   \
                             tr_x_xxxx_xyzzzz,  \
                             tr_x_xxxx_xzzzz,   \
                             tr_x_xxxx_xzzzzz,  \
                             tr_x_xxxx_yyyyy,   \
                             tr_x_xxxx_yyyyyy,  \
                             tr_x_xxxx_yyyyyz,  \
                             tr_x_xxxx_yyyyz,   \
                             tr_x_xxxx_yyyyzz,  \
                             tr_x_xxxx_yyyzz,   \
                             tr_x_xxxx_yyyzzz,  \
                             tr_x_xxxx_yyzzz,   \
                             tr_x_xxxx_yyzzzz,  \
                             tr_x_xxxx_yzzzz,   \
                             tr_x_xxxx_yzzzzz,  \
                             tr_x_xxxx_zzzzz,   \
                             tr_x_xxxx_zzzzzz,  \
                             tr_x_xxxxy_xxxxxx, \
                             tr_x_xxxxy_xxxxxy, \
                             tr_x_xxxxy_xxxxxz, \
                             tr_x_xxxxy_xxxxyy, \
                             tr_x_xxxxy_xxxxyz, \
                             tr_x_xxxxy_xxxxzz, \
                             tr_x_xxxxy_xxxyyy, \
                             tr_x_xxxxy_xxxyyz, \
                             tr_x_xxxxy_xxxyzz, \
                             tr_x_xxxxy_xxxzzz, \
                             tr_x_xxxxy_xxyyyy, \
                             tr_x_xxxxy_xxyyyz, \
                             tr_x_xxxxy_xxyyzz, \
                             tr_x_xxxxy_xxyzzz, \
                             tr_x_xxxxy_xxzzzz, \
                             tr_x_xxxxy_xyyyyy, \
                             tr_x_xxxxy_xyyyyz, \
                             tr_x_xxxxy_xyyyzz, \
                             tr_x_xxxxy_xyyzzz, \
                             tr_x_xxxxy_xyzzzz, \
                             tr_x_xxxxy_xzzzzz, \
                             tr_x_xxxxy_yyyyyy, \
                             tr_x_xxxxy_yyyyyz, \
                             tr_x_xxxxy_yyyyzz, \
                             tr_x_xxxxy_yyyzzz, \
                             tr_x_xxxxy_yyzzzz, \
                             tr_x_xxxxy_yzzzzz, \
                             tr_x_xxxxy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxy_xxxxxx[i] = tr_x_xxxx_xxxxxx[i] * pa_y[i];

        tr_x_xxxxy_xxxxxy[i] = tr_x_xxxx_xxxxx[i] * fe_0 + tr_x_xxxx_xxxxxy[i] * pa_y[i];

        tr_x_xxxxy_xxxxxz[i] = tr_x_xxxx_xxxxxz[i] * pa_y[i];

        tr_x_xxxxy_xxxxyy[i] = 2.0 * tr_x_xxxx_xxxxy[i] * fe_0 + tr_x_xxxx_xxxxyy[i] * pa_y[i];

        tr_x_xxxxy_xxxxyz[i] = tr_x_xxxx_xxxxz[i] * fe_0 + tr_x_xxxx_xxxxyz[i] * pa_y[i];

        tr_x_xxxxy_xxxxzz[i] = tr_x_xxxx_xxxxzz[i] * pa_y[i];

        tr_x_xxxxy_xxxyyy[i] = 3.0 * tr_x_xxxx_xxxyy[i] * fe_0 + tr_x_xxxx_xxxyyy[i] * pa_y[i];

        tr_x_xxxxy_xxxyyz[i] = 2.0 * tr_x_xxxx_xxxyz[i] * fe_0 + tr_x_xxxx_xxxyyz[i] * pa_y[i];

        tr_x_xxxxy_xxxyzz[i] = tr_x_xxxx_xxxzz[i] * fe_0 + tr_x_xxxx_xxxyzz[i] * pa_y[i];

        tr_x_xxxxy_xxxzzz[i] = tr_x_xxxx_xxxzzz[i] * pa_y[i];

        tr_x_xxxxy_xxyyyy[i] = 4.0 * tr_x_xxxx_xxyyy[i] * fe_0 + tr_x_xxxx_xxyyyy[i] * pa_y[i];

        tr_x_xxxxy_xxyyyz[i] = 3.0 * tr_x_xxxx_xxyyz[i] * fe_0 + tr_x_xxxx_xxyyyz[i] * pa_y[i];

        tr_x_xxxxy_xxyyzz[i] = 2.0 * tr_x_xxxx_xxyzz[i] * fe_0 + tr_x_xxxx_xxyyzz[i] * pa_y[i];

        tr_x_xxxxy_xxyzzz[i] = tr_x_xxxx_xxzzz[i] * fe_0 + tr_x_xxxx_xxyzzz[i] * pa_y[i];

        tr_x_xxxxy_xxzzzz[i] = tr_x_xxxx_xxzzzz[i] * pa_y[i];

        tr_x_xxxxy_xyyyyy[i] = 5.0 * tr_x_xxxx_xyyyy[i] * fe_0 + tr_x_xxxx_xyyyyy[i] * pa_y[i];

        tr_x_xxxxy_xyyyyz[i] = 4.0 * tr_x_xxxx_xyyyz[i] * fe_0 + tr_x_xxxx_xyyyyz[i] * pa_y[i];

        tr_x_xxxxy_xyyyzz[i] = 3.0 * tr_x_xxxx_xyyzz[i] * fe_0 + tr_x_xxxx_xyyyzz[i] * pa_y[i];

        tr_x_xxxxy_xyyzzz[i] = 2.0 * tr_x_xxxx_xyzzz[i] * fe_0 + tr_x_xxxx_xyyzzz[i] * pa_y[i];

        tr_x_xxxxy_xyzzzz[i] = tr_x_xxxx_xzzzz[i] * fe_0 + tr_x_xxxx_xyzzzz[i] * pa_y[i];

        tr_x_xxxxy_xzzzzz[i] = tr_x_xxxx_xzzzzz[i] * pa_y[i];

        tr_x_xxxxy_yyyyyy[i] = 6.0 * tr_x_xxxx_yyyyy[i] * fe_0 + tr_x_xxxx_yyyyyy[i] * pa_y[i];

        tr_x_xxxxy_yyyyyz[i] = 5.0 * tr_x_xxxx_yyyyz[i] * fe_0 + tr_x_xxxx_yyyyyz[i] * pa_y[i];

        tr_x_xxxxy_yyyyzz[i] = 4.0 * tr_x_xxxx_yyyzz[i] * fe_0 + tr_x_xxxx_yyyyzz[i] * pa_y[i];

        tr_x_xxxxy_yyyzzz[i] = 3.0 * tr_x_xxxx_yyzzz[i] * fe_0 + tr_x_xxxx_yyyzzz[i] * pa_y[i];

        tr_x_xxxxy_yyzzzz[i] = 2.0 * tr_x_xxxx_yzzzz[i] * fe_0 + tr_x_xxxx_yyzzzz[i] * pa_y[i];

        tr_x_xxxxy_yzzzzz[i] = tr_x_xxxx_zzzzz[i] * fe_0 + tr_x_xxxx_yzzzzz[i] * pa_y[i];

        tr_x_xxxxy_zzzzzz[i] = tr_x_xxxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : HI

    auto tr_x_xxxxz_xxxxxx = pbuffer.data(idx_dip_hi + 56);

    auto tr_x_xxxxz_xxxxxy = pbuffer.data(idx_dip_hi + 57);

    auto tr_x_xxxxz_xxxxxz = pbuffer.data(idx_dip_hi + 58);

    auto tr_x_xxxxz_xxxxyy = pbuffer.data(idx_dip_hi + 59);

    auto tr_x_xxxxz_xxxxyz = pbuffer.data(idx_dip_hi + 60);

    auto tr_x_xxxxz_xxxxzz = pbuffer.data(idx_dip_hi + 61);

    auto tr_x_xxxxz_xxxyyy = pbuffer.data(idx_dip_hi + 62);

    auto tr_x_xxxxz_xxxyyz = pbuffer.data(idx_dip_hi + 63);

    auto tr_x_xxxxz_xxxyzz = pbuffer.data(idx_dip_hi + 64);

    auto tr_x_xxxxz_xxxzzz = pbuffer.data(idx_dip_hi + 65);

    auto tr_x_xxxxz_xxyyyy = pbuffer.data(idx_dip_hi + 66);

    auto tr_x_xxxxz_xxyyyz = pbuffer.data(idx_dip_hi + 67);

    auto tr_x_xxxxz_xxyyzz = pbuffer.data(idx_dip_hi + 68);

    auto tr_x_xxxxz_xxyzzz = pbuffer.data(idx_dip_hi + 69);

    auto tr_x_xxxxz_xxzzzz = pbuffer.data(idx_dip_hi + 70);

    auto tr_x_xxxxz_xyyyyy = pbuffer.data(idx_dip_hi + 71);

    auto tr_x_xxxxz_xyyyyz = pbuffer.data(idx_dip_hi + 72);

    auto tr_x_xxxxz_xyyyzz = pbuffer.data(idx_dip_hi + 73);

    auto tr_x_xxxxz_xyyzzz = pbuffer.data(idx_dip_hi + 74);

    auto tr_x_xxxxz_xyzzzz = pbuffer.data(idx_dip_hi + 75);

    auto tr_x_xxxxz_xzzzzz = pbuffer.data(idx_dip_hi + 76);

    auto tr_x_xxxxz_yyyyyy = pbuffer.data(idx_dip_hi + 77);

    auto tr_x_xxxxz_yyyyyz = pbuffer.data(idx_dip_hi + 78);

    auto tr_x_xxxxz_yyyyzz = pbuffer.data(idx_dip_hi + 79);

    auto tr_x_xxxxz_yyyzzz = pbuffer.data(idx_dip_hi + 80);

    auto tr_x_xxxxz_yyzzzz = pbuffer.data(idx_dip_hi + 81);

    auto tr_x_xxxxz_yzzzzz = pbuffer.data(idx_dip_hi + 82);

    auto tr_x_xxxxz_zzzzzz = pbuffer.data(idx_dip_hi + 83);

#pragma omp simd aligned(pa_z,                  \
                             tr_x_xxxx_xxxxx,   \
                             tr_x_xxxx_xxxxxx,  \
                             tr_x_xxxx_xxxxxy,  \
                             tr_x_xxxx_xxxxxz,  \
                             tr_x_xxxx_xxxxy,   \
                             tr_x_xxxx_xxxxyy,  \
                             tr_x_xxxx_xxxxyz,  \
                             tr_x_xxxx_xxxxz,   \
                             tr_x_xxxx_xxxxzz,  \
                             tr_x_xxxx_xxxyy,   \
                             tr_x_xxxx_xxxyyy,  \
                             tr_x_xxxx_xxxyyz,  \
                             tr_x_xxxx_xxxyz,   \
                             tr_x_xxxx_xxxyzz,  \
                             tr_x_xxxx_xxxzz,   \
                             tr_x_xxxx_xxxzzz,  \
                             tr_x_xxxx_xxyyy,   \
                             tr_x_xxxx_xxyyyy,  \
                             tr_x_xxxx_xxyyyz,  \
                             tr_x_xxxx_xxyyz,   \
                             tr_x_xxxx_xxyyzz,  \
                             tr_x_xxxx_xxyzz,   \
                             tr_x_xxxx_xxyzzz,  \
                             tr_x_xxxx_xxzzz,   \
                             tr_x_xxxx_xxzzzz,  \
                             tr_x_xxxx_xyyyy,   \
                             tr_x_xxxx_xyyyyy,  \
                             tr_x_xxxx_xyyyyz,  \
                             tr_x_xxxx_xyyyz,   \
                             tr_x_xxxx_xyyyzz,  \
                             tr_x_xxxx_xyyzz,   \
                             tr_x_xxxx_xyyzzz,  \
                             tr_x_xxxx_xyzzz,   \
                             tr_x_xxxx_xyzzzz,  \
                             tr_x_xxxx_xzzzz,   \
                             tr_x_xxxx_xzzzzz,  \
                             tr_x_xxxx_yyyyy,   \
                             tr_x_xxxx_yyyyyy,  \
                             tr_x_xxxx_yyyyyz,  \
                             tr_x_xxxx_yyyyz,   \
                             tr_x_xxxx_yyyyzz,  \
                             tr_x_xxxx_yyyzz,   \
                             tr_x_xxxx_yyyzzz,  \
                             tr_x_xxxx_yyzzz,   \
                             tr_x_xxxx_yyzzzz,  \
                             tr_x_xxxx_yzzzz,   \
                             tr_x_xxxx_yzzzzz,  \
                             tr_x_xxxx_zzzzz,   \
                             tr_x_xxxx_zzzzzz,  \
                             tr_x_xxxxz_xxxxxx, \
                             tr_x_xxxxz_xxxxxy, \
                             tr_x_xxxxz_xxxxxz, \
                             tr_x_xxxxz_xxxxyy, \
                             tr_x_xxxxz_xxxxyz, \
                             tr_x_xxxxz_xxxxzz, \
                             tr_x_xxxxz_xxxyyy, \
                             tr_x_xxxxz_xxxyyz, \
                             tr_x_xxxxz_xxxyzz, \
                             tr_x_xxxxz_xxxzzz, \
                             tr_x_xxxxz_xxyyyy, \
                             tr_x_xxxxz_xxyyyz, \
                             tr_x_xxxxz_xxyyzz, \
                             tr_x_xxxxz_xxyzzz, \
                             tr_x_xxxxz_xxzzzz, \
                             tr_x_xxxxz_xyyyyy, \
                             tr_x_xxxxz_xyyyyz, \
                             tr_x_xxxxz_xyyyzz, \
                             tr_x_xxxxz_xyyzzz, \
                             tr_x_xxxxz_xyzzzz, \
                             tr_x_xxxxz_xzzzzz, \
                             tr_x_xxxxz_yyyyyy, \
                             tr_x_xxxxz_yyyyyz, \
                             tr_x_xxxxz_yyyyzz, \
                             tr_x_xxxxz_yyyzzz, \
                             tr_x_xxxxz_yyzzzz, \
                             tr_x_xxxxz_yzzzzz, \
                             tr_x_xxxxz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxz_xxxxxx[i] = tr_x_xxxx_xxxxxx[i] * pa_z[i];

        tr_x_xxxxz_xxxxxy[i] = tr_x_xxxx_xxxxxy[i] * pa_z[i];

        tr_x_xxxxz_xxxxxz[i] = tr_x_xxxx_xxxxx[i] * fe_0 + tr_x_xxxx_xxxxxz[i] * pa_z[i];

        tr_x_xxxxz_xxxxyy[i] = tr_x_xxxx_xxxxyy[i] * pa_z[i];

        tr_x_xxxxz_xxxxyz[i] = tr_x_xxxx_xxxxy[i] * fe_0 + tr_x_xxxx_xxxxyz[i] * pa_z[i];

        tr_x_xxxxz_xxxxzz[i] = 2.0 * tr_x_xxxx_xxxxz[i] * fe_0 + tr_x_xxxx_xxxxzz[i] * pa_z[i];

        tr_x_xxxxz_xxxyyy[i] = tr_x_xxxx_xxxyyy[i] * pa_z[i];

        tr_x_xxxxz_xxxyyz[i] = tr_x_xxxx_xxxyy[i] * fe_0 + tr_x_xxxx_xxxyyz[i] * pa_z[i];

        tr_x_xxxxz_xxxyzz[i] = 2.0 * tr_x_xxxx_xxxyz[i] * fe_0 + tr_x_xxxx_xxxyzz[i] * pa_z[i];

        tr_x_xxxxz_xxxzzz[i] = 3.0 * tr_x_xxxx_xxxzz[i] * fe_0 + tr_x_xxxx_xxxzzz[i] * pa_z[i];

        tr_x_xxxxz_xxyyyy[i] = tr_x_xxxx_xxyyyy[i] * pa_z[i];

        tr_x_xxxxz_xxyyyz[i] = tr_x_xxxx_xxyyy[i] * fe_0 + tr_x_xxxx_xxyyyz[i] * pa_z[i];

        tr_x_xxxxz_xxyyzz[i] = 2.0 * tr_x_xxxx_xxyyz[i] * fe_0 + tr_x_xxxx_xxyyzz[i] * pa_z[i];

        tr_x_xxxxz_xxyzzz[i] = 3.0 * tr_x_xxxx_xxyzz[i] * fe_0 + tr_x_xxxx_xxyzzz[i] * pa_z[i];

        tr_x_xxxxz_xxzzzz[i] = 4.0 * tr_x_xxxx_xxzzz[i] * fe_0 + tr_x_xxxx_xxzzzz[i] * pa_z[i];

        tr_x_xxxxz_xyyyyy[i] = tr_x_xxxx_xyyyyy[i] * pa_z[i];

        tr_x_xxxxz_xyyyyz[i] = tr_x_xxxx_xyyyy[i] * fe_0 + tr_x_xxxx_xyyyyz[i] * pa_z[i];

        tr_x_xxxxz_xyyyzz[i] = 2.0 * tr_x_xxxx_xyyyz[i] * fe_0 + tr_x_xxxx_xyyyzz[i] * pa_z[i];

        tr_x_xxxxz_xyyzzz[i] = 3.0 * tr_x_xxxx_xyyzz[i] * fe_0 + tr_x_xxxx_xyyzzz[i] * pa_z[i];

        tr_x_xxxxz_xyzzzz[i] = 4.0 * tr_x_xxxx_xyzzz[i] * fe_0 + tr_x_xxxx_xyzzzz[i] * pa_z[i];

        tr_x_xxxxz_xzzzzz[i] = 5.0 * tr_x_xxxx_xzzzz[i] * fe_0 + tr_x_xxxx_xzzzzz[i] * pa_z[i];

        tr_x_xxxxz_yyyyyy[i] = tr_x_xxxx_yyyyyy[i] * pa_z[i];

        tr_x_xxxxz_yyyyyz[i] = tr_x_xxxx_yyyyy[i] * fe_0 + tr_x_xxxx_yyyyyz[i] * pa_z[i];

        tr_x_xxxxz_yyyyzz[i] = 2.0 * tr_x_xxxx_yyyyz[i] * fe_0 + tr_x_xxxx_yyyyzz[i] * pa_z[i];

        tr_x_xxxxz_yyyzzz[i] = 3.0 * tr_x_xxxx_yyyzz[i] * fe_0 + tr_x_xxxx_yyyzzz[i] * pa_z[i];

        tr_x_xxxxz_yyzzzz[i] = 4.0 * tr_x_xxxx_yyzzz[i] * fe_0 + tr_x_xxxx_yyzzzz[i] * pa_z[i];

        tr_x_xxxxz_yzzzzz[i] = 5.0 * tr_x_xxxx_yzzzz[i] * fe_0 + tr_x_xxxx_yzzzzz[i] * pa_z[i];

        tr_x_xxxxz_zzzzzz[i] = 6.0 * tr_x_xxxx_zzzzz[i] * fe_0 + tr_x_xxxx_zzzzzz[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : HI

    auto tr_x_xxxyy_xxxxxx = pbuffer.data(idx_dip_hi + 84);

    auto tr_x_xxxyy_xxxxxy = pbuffer.data(idx_dip_hi + 85);

    auto tr_x_xxxyy_xxxxxz = pbuffer.data(idx_dip_hi + 86);

    auto tr_x_xxxyy_xxxxyy = pbuffer.data(idx_dip_hi + 87);

    auto tr_x_xxxyy_xxxxyz = pbuffer.data(idx_dip_hi + 88);

    auto tr_x_xxxyy_xxxxzz = pbuffer.data(idx_dip_hi + 89);

    auto tr_x_xxxyy_xxxyyy = pbuffer.data(idx_dip_hi + 90);

    auto tr_x_xxxyy_xxxyyz = pbuffer.data(idx_dip_hi + 91);

    auto tr_x_xxxyy_xxxyzz = pbuffer.data(idx_dip_hi + 92);

    auto tr_x_xxxyy_xxxzzz = pbuffer.data(idx_dip_hi + 93);

    auto tr_x_xxxyy_xxyyyy = pbuffer.data(idx_dip_hi + 94);

    auto tr_x_xxxyy_xxyyyz = pbuffer.data(idx_dip_hi + 95);

    auto tr_x_xxxyy_xxyyzz = pbuffer.data(idx_dip_hi + 96);

    auto tr_x_xxxyy_xxyzzz = pbuffer.data(idx_dip_hi + 97);

    auto tr_x_xxxyy_xxzzzz = pbuffer.data(idx_dip_hi + 98);

    auto tr_x_xxxyy_xyyyyy = pbuffer.data(idx_dip_hi + 99);

    auto tr_x_xxxyy_xyyyyz = pbuffer.data(idx_dip_hi + 100);

    auto tr_x_xxxyy_xyyyzz = pbuffer.data(idx_dip_hi + 101);

    auto tr_x_xxxyy_xyyzzz = pbuffer.data(idx_dip_hi + 102);

    auto tr_x_xxxyy_xyzzzz = pbuffer.data(idx_dip_hi + 103);

    auto tr_x_xxxyy_xzzzzz = pbuffer.data(idx_dip_hi + 104);

    auto tr_x_xxxyy_yyyyyy = pbuffer.data(idx_dip_hi + 105);

    auto tr_x_xxxyy_yyyyyz = pbuffer.data(idx_dip_hi + 106);

    auto tr_x_xxxyy_yyyyzz = pbuffer.data(idx_dip_hi + 107);

    auto tr_x_xxxyy_yyyzzz = pbuffer.data(idx_dip_hi + 108);

    auto tr_x_xxxyy_yyzzzz = pbuffer.data(idx_dip_hi + 109);

    auto tr_x_xxxyy_yzzzzz = pbuffer.data(idx_dip_hi + 110);

    auto tr_x_xxxyy_zzzzzz = pbuffer.data(idx_dip_hi + 111);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_x_xxx_xxxxxx,   \
                             tr_x_xxx_xxxxxy,   \
                             tr_x_xxx_xxxxxz,   \
                             tr_x_xxx_xxxxyy,   \
                             tr_x_xxx_xxxxyz,   \
                             tr_x_xxx_xxxxzz,   \
                             tr_x_xxx_xxxyyy,   \
                             tr_x_xxx_xxxyyz,   \
                             tr_x_xxx_xxxyzz,   \
                             tr_x_xxx_xxxzzz,   \
                             tr_x_xxx_xxyyyy,   \
                             tr_x_xxx_xxyyyz,   \
                             tr_x_xxx_xxyyzz,   \
                             tr_x_xxx_xxyzzz,   \
                             tr_x_xxx_xxzzzz,   \
                             tr_x_xxx_xyyyyy,   \
                             tr_x_xxx_xyyyyz,   \
                             tr_x_xxx_xyyyzz,   \
                             tr_x_xxx_xyyzzz,   \
                             tr_x_xxx_xyzzzz,   \
                             tr_x_xxx_xzzzzz,   \
                             tr_x_xxx_zzzzzz,   \
                             tr_x_xxxy_xxxxx,   \
                             tr_x_xxxy_xxxxxx,  \
                             tr_x_xxxy_xxxxxy,  \
                             tr_x_xxxy_xxxxxz,  \
                             tr_x_xxxy_xxxxy,   \
                             tr_x_xxxy_xxxxyy,  \
                             tr_x_xxxy_xxxxyz,  \
                             tr_x_xxxy_xxxxz,   \
                             tr_x_xxxy_xxxxzz,  \
                             tr_x_xxxy_xxxyy,   \
                             tr_x_xxxy_xxxyyy,  \
                             tr_x_xxxy_xxxyyz,  \
                             tr_x_xxxy_xxxyz,   \
                             tr_x_xxxy_xxxyzz,  \
                             tr_x_xxxy_xxxzz,   \
                             tr_x_xxxy_xxxzzz,  \
                             tr_x_xxxy_xxyyy,   \
                             tr_x_xxxy_xxyyyy,  \
                             tr_x_xxxy_xxyyyz,  \
                             tr_x_xxxy_xxyyz,   \
                             tr_x_xxxy_xxyyzz,  \
                             tr_x_xxxy_xxyzz,   \
                             tr_x_xxxy_xxyzzz,  \
                             tr_x_xxxy_xxzzz,   \
                             tr_x_xxxy_xxzzzz,  \
                             tr_x_xxxy_xyyyy,   \
                             tr_x_xxxy_xyyyyy,  \
                             tr_x_xxxy_xyyyyz,  \
                             tr_x_xxxy_xyyyz,   \
                             tr_x_xxxy_xyyyzz,  \
                             tr_x_xxxy_xyyzz,   \
                             tr_x_xxxy_xyyzzz,  \
                             tr_x_xxxy_xyzzz,   \
                             tr_x_xxxy_xyzzzz,  \
                             tr_x_xxxy_xzzzz,   \
                             tr_x_xxxy_xzzzzz,  \
                             tr_x_xxxy_zzzzzz,  \
                             tr_x_xxxyy_xxxxxx, \
                             tr_x_xxxyy_xxxxxy, \
                             tr_x_xxxyy_xxxxxz, \
                             tr_x_xxxyy_xxxxyy, \
                             tr_x_xxxyy_xxxxyz, \
                             tr_x_xxxyy_xxxxzz, \
                             tr_x_xxxyy_xxxyyy, \
                             tr_x_xxxyy_xxxyyz, \
                             tr_x_xxxyy_xxxyzz, \
                             tr_x_xxxyy_xxxzzz, \
                             tr_x_xxxyy_xxyyyy, \
                             tr_x_xxxyy_xxyyyz, \
                             tr_x_xxxyy_xxyyzz, \
                             tr_x_xxxyy_xxyzzz, \
                             tr_x_xxxyy_xxzzzz, \
                             tr_x_xxxyy_xyyyyy, \
                             tr_x_xxxyy_xyyyyz, \
                             tr_x_xxxyy_xyyyzz, \
                             tr_x_xxxyy_xyyzzz, \
                             tr_x_xxxyy_xyzzzz, \
                             tr_x_xxxyy_xzzzzz, \
                             tr_x_xxxyy_yyyyyy, \
                             tr_x_xxxyy_yyyyyz, \
                             tr_x_xxxyy_yyyyzz, \
                             tr_x_xxxyy_yyyzzz, \
                             tr_x_xxxyy_yyzzzz, \
                             tr_x_xxxyy_yzzzzz, \
                             tr_x_xxxyy_zzzzzz, \
                             tr_x_xxyy_yyyyyy,  \
                             tr_x_xxyy_yyyyyz,  \
                             tr_x_xxyy_yyyyzz,  \
                             tr_x_xxyy_yyyzzz,  \
                             tr_x_xxyy_yyzzzz,  \
                             tr_x_xxyy_yzzzzz,  \
                             tr_x_xyy_yyyyyy,   \
                             tr_x_xyy_yyyyyz,   \
                             tr_x_xyy_yyyyzz,   \
                             tr_x_xyy_yyyzzz,   \
                             tr_x_xyy_yyzzzz,   \
                             tr_x_xyy_yzzzzz,   \
                             ts_xxyy_yyyyyy,    \
                             ts_xxyy_yyyyyz,    \
                             ts_xxyy_yyyyzz,    \
                             ts_xxyy_yyyzzz,    \
                             ts_xxyy_yyzzzz,    \
                             ts_xxyy_yzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyy_xxxxxx[i] = tr_x_xxx_xxxxxx[i] * fe_0 + tr_x_xxxy_xxxxxx[i] * pa_y[i];

        tr_x_xxxyy_xxxxxy[i] = tr_x_xxx_xxxxxy[i] * fe_0 + tr_x_xxxy_xxxxx[i] * fe_0 + tr_x_xxxy_xxxxxy[i] * pa_y[i];

        tr_x_xxxyy_xxxxxz[i] = tr_x_xxx_xxxxxz[i] * fe_0 + tr_x_xxxy_xxxxxz[i] * pa_y[i];

        tr_x_xxxyy_xxxxyy[i] = tr_x_xxx_xxxxyy[i] * fe_0 + 2.0 * tr_x_xxxy_xxxxy[i] * fe_0 + tr_x_xxxy_xxxxyy[i] * pa_y[i];

        tr_x_xxxyy_xxxxyz[i] = tr_x_xxx_xxxxyz[i] * fe_0 + tr_x_xxxy_xxxxz[i] * fe_0 + tr_x_xxxy_xxxxyz[i] * pa_y[i];

        tr_x_xxxyy_xxxxzz[i] = tr_x_xxx_xxxxzz[i] * fe_0 + tr_x_xxxy_xxxxzz[i] * pa_y[i];

        tr_x_xxxyy_xxxyyy[i] = tr_x_xxx_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxxy_xxxyy[i] * fe_0 + tr_x_xxxy_xxxyyy[i] * pa_y[i];

        tr_x_xxxyy_xxxyyz[i] = tr_x_xxx_xxxyyz[i] * fe_0 + 2.0 * tr_x_xxxy_xxxyz[i] * fe_0 + tr_x_xxxy_xxxyyz[i] * pa_y[i];

        tr_x_xxxyy_xxxyzz[i] = tr_x_xxx_xxxyzz[i] * fe_0 + tr_x_xxxy_xxxzz[i] * fe_0 + tr_x_xxxy_xxxyzz[i] * pa_y[i];

        tr_x_xxxyy_xxxzzz[i] = tr_x_xxx_xxxzzz[i] * fe_0 + tr_x_xxxy_xxxzzz[i] * pa_y[i];

        tr_x_xxxyy_xxyyyy[i] = tr_x_xxx_xxyyyy[i] * fe_0 + 4.0 * tr_x_xxxy_xxyyy[i] * fe_0 + tr_x_xxxy_xxyyyy[i] * pa_y[i];

        tr_x_xxxyy_xxyyyz[i] = tr_x_xxx_xxyyyz[i] * fe_0 + 3.0 * tr_x_xxxy_xxyyz[i] * fe_0 + tr_x_xxxy_xxyyyz[i] * pa_y[i];

        tr_x_xxxyy_xxyyzz[i] = tr_x_xxx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxy_xxyzz[i] * fe_0 + tr_x_xxxy_xxyyzz[i] * pa_y[i];

        tr_x_xxxyy_xxyzzz[i] = tr_x_xxx_xxyzzz[i] * fe_0 + tr_x_xxxy_xxzzz[i] * fe_0 + tr_x_xxxy_xxyzzz[i] * pa_y[i];

        tr_x_xxxyy_xxzzzz[i] = tr_x_xxx_xxzzzz[i] * fe_0 + tr_x_xxxy_xxzzzz[i] * pa_y[i];

        tr_x_xxxyy_xyyyyy[i] = tr_x_xxx_xyyyyy[i] * fe_0 + 5.0 * tr_x_xxxy_xyyyy[i] * fe_0 + tr_x_xxxy_xyyyyy[i] * pa_y[i];

        tr_x_xxxyy_xyyyyz[i] = tr_x_xxx_xyyyyz[i] * fe_0 + 4.0 * tr_x_xxxy_xyyyz[i] * fe_0 + tr_x_xxxy_xyyyyz[i] * pa_y[i];

        tr_x_xxxyy_xyyyzz[i] = tr_x_xxx_xyyyzz[i] * fe_0 + 3.0 * tr_x_xxxy_xyyzz[i] * fe_0 + tr_x_xxxy_xyyyzz[i] * pa_y[i];

        tr_x_xxxyy_xyyzzz[i] = tr_x_xxx_xyyzzz[i] * fe_0 + 2.0 * tr_x_xxxy_xyzzz[i] * fe_0 + tr_x_xxxy_xyyzzz[i] * pa_y[i];

        tr_x_xxxyy_xyzzzz[i] = tr_x_xxx_xyzzzz[i] * fe_0 + tr_x_xxxy_xzzzz[i] * fe_0 + tr_x_xxxy_xyzzzz[i] * pa_y[i];

        tr_x_xxxyy_xzzzzz[i] = tr_x_xxx_xzzzzz[i] * fe_0 + tr_x_xxxy_xzzzzz[i] * pa_y[i];

        tr_x_xxxyy_yyyyyy[i] = 2.0 * tr_x_xyy_yyyyyy[i] * fe_0 + ts_xxyy_yyyyyy[i] * fe_0 + tr_x_xxyy_yyyyyy[i] * pa_x[i];

        tr_x_xxxyy_yyyyyz[i] = 2.0 * tr_x_xyy_yyyyyz[i] * fe_0 + ts_xxyy_yyyyyz[i] * fe_0 + tr_x_xxyy_yyyyyz[i] * pa_x[i];

        tr_x_xxxyy_yyyyzz[i] = 2.0 * tr_x_xyy_yyyyzz[i] * fe_0 + ts_xxyy_yyyyzz[i] * fe_0 + tr_x_xxyy_yyyyzz[i] * pa_x[i];

        tr_x_xxxyy_yyyzzz[i] = 2.0 * tr_x_xyy_yyyzzz[i] * fe_0 + ts_xxyy_yyyzzz[i] * fe_0 + tr_x_xxyy_yyyzzz[i] * pa_x[i];

        tr_x_xxxyy_yyzzzz[i] = 2.0 * tr_x_xyy_yyzzzz[i] * fe_0 + ts_xxyy_yyzzzz[i] * fe_0 + tr_x_xxyy_yyzzzz[i] * pa_x[i];

        tr_x_xxxyy_yzzzzz[i] = 2.0 * tr_x_xyy_yzzzzz[i] * fe_0 + ts_xxyy_yzzzzz[i] * fe_0 + tr_x_xxyy_yzzzzz[i] * pa_x[i];

        tr_x_xxxyy_zzzzzz[i] = tr_x_xxx_zzzzzz[i] * fe_0 + tr_x_xxxy_zzzzzz[i] * pa_y[i];
    }

    // Set up 112-140 components of targeted buffer : HI

    auto tr_x_xxxyz_xxxxxx = pbuffer.data(idx_dip_hi + 112);

    auto tr_x_xxxyz_xxxxxy = pbuffer.data(idx_dip_hi + 113);

    auto tr_x_xxxyz_xxxxxz = pbuffer.data(idx_dip_hi + 114);

    auto tr_x_xxxyz_xxxxyy = pbuffer.data(idx_dip_hi + 115);

    auto tr_x_xxxyz_xxxxyz = pbuffer.data(idx_dip_hi + 116);

    auto tr_x_xxxyz_xxxxzz = pbuffer.data(idx_dip_hi + 117);

    auto tr_x_xxxyz_xxxyyy = pbuffer.data(idx_dip_hi + 118);

    auto tr_x_xxxyz_xxxyyz = pbuffer.data(idx_dip_hi + 119);

    auto tr_x_xxxyz_xxxyzz = pbuffer.data(idx_dip_hi + 120);

    auto tr_x_xxxyz_xxxzzz = pbuffer.data(idx_dip_hi + 121);

    auto tr_x_xxxyz_xxyyyy = pbuffer.data(idx_dip_hi + 122);

    auto tr_x_xxxyz_xxyyyz = pbuffer.data(idx_dip_hi + 123);

    auto tr_x_xxxyz_xxyyzz = pbuffer.data(idx_dip_hi + 124);

    auto tr_x_xxxyz_xxyzzz = pbuffer.data(idx_dip_hi + 125);

    auto tr_x_xxxyz_xxzzzz = pbuffer.data(idx_dip_hi + 126);

    auto tr_x_xxxyz_xyyyyy = pbuffer.data(idx_dip_hi + 127);

    auto tr_x_xxxyz_xyyyyz = pbuffer.data(idx_dip_hi + 128);

    auto tr_x_xxxyz_xyyyzz = pbuffer.data(idx_dip_hi + 129);

    auto tr_x_xxxyz_xyyzzz = pbuffer.data(idx_dip_hi + 130);

    auto tr_x_xxxyz_xyzzzz = pbuffer.data(idx_dip_hi + 131);

    auto tr_x_xxxyz_xzzzzz = pbuffer.data(idx_dip_hi + 132);

    auto tr_x_xxxyz_yyyyyy = pbuffer.data(idx_dip_hi + 133);

    auto tr_x_xxxyz_yyyyyz = pbuffer.data(idx_dip_hi + 134);

    auto tr_x_xxxyz_yyyyzz = pbuffer.data(idx_dip_hi + 135);

    auto tr_x_xxxyz_yyyzzz = pbuffer.data(idx_dip_hi + 136);

    auto tr_x_xxxyz_yyzzzz = pbuffer.data(idx_dip_hi + 137);

    auto tr_x_xxxyz_yzzzzz = pbuffer.data(idx_dip_hi + 138);

    auto tr_x_xxxyz_zzzzzz = pbuffer.data(idx_dip_hi + 139);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_x_xxxy_xxxxxy,  \
                             tr_x_xxxy_xxxxyy,  \
                             tr_x_xxxy_xxxyyy,  \
                             tr_x_xxxy_xxyyyy,  \
                             tr_x_xxxy_xyyyyy,  \
                             tr_x_xxxy_yyyyyy,  \
                             tr_x_xxxyz_xxxxxx, \
                             tr_x_xxxyz_xxxxxy, \
                             tr_x_xxxyz_xxxxxz, \
                             tr_x_xxxyz_xxxxyy, \
                             tr_x_xxxyz_xxxxyz, \
                             tr_x_xxxyz_xxxxzz, \
                             tr_x_xxxyz_xxxyyy, \
                             tr_x_xxxyz_xxxyyz, \
                             tr_x_xxxyz_xxxyzz, \
                             tr_x_xxxyz_xxxzzz, \
                             tr_x_xxxyz_xxyyyy, \
                             tr_x_xxxyz_xxyyyz, \
                             tr_x_xxxyz_xxyyzz, \
                             tr_x_xxxyz_xxyzzz, \
                             tr_x_xxxyz_xxzzzz, \
                             tr_x_xxxyz_xyyyyy, \
                             tr_x_xxxyz_xyyyyz, \
                             tr_x_xxxyz_xyyyzz, \
                             tr_x_xxxyz_xyyzzz, \
                             tr_x_xxxyz_xyzzzz, \
                             tr_x_xxxyz_xzzzzz, \
                             tr_x_xxxyz_yyyyyy, \
                             tr_x_xxxyz_yyyyyz, \
                             tr_x_xxxyz_yyyyzz, \
                             tr_x_xxxyz_yyyzzz, \
                             tr_x_xxxyz_yyzzzz, \
                             tr_x_xxxyz_yzzzzz, \
                             tr_x_xxxyz_zzzzzz, \
                             tr_x_xxxz_xxxxxx,  \
                             tr_x_xxxz_xxxxxz,  \
                             tr_x_xxxz_xxxxyz,  \
                             tr_x_xxxz_xxxxz,   \
                             tr_x_xxxz_xxxxzz,  \
                             tr_x_xxxz_xxxyyz,  \
                             tr_x_xxxz_xxxyz,   \
                             tr_x_xxxz_xxxyzz,  \
                             tr_x_xxxz_xxxzz,   \
                             tr_x_xxxz_xxxzzz,  \
                             tr_x_xxxz_xxyyyz,  \
                             tr_x_xxxz_xxyyz,   \
                             tr_x_xxxz_xxyyzz,  \
                             tr_x_xxxz_xxyzz,   \
                             tr_x_xxxz_xxyzzz,  \
                             tr_x_xxxz_xxzzz,   \
                             tr_x_xxxz_xxzzzz,  \
                             tr_x_xxxz_xyyyyz,  \
                             tr_x_xxxz_xyyyz,   \
                             tr_x_xxxz_xyyyzz,  \
                             tr_x_xxxz_xyyzz,   \
                             tr_x_xxxz_xyyzzz,  \
                             tr_x_xxxz_xyzzz,   \
                             tr_x_xxxz_xyzzzz,  \
                             tr_x_xxxz_xzzzz,   \
                             tr_x_xxxz_xzzzzz,  \
                             tr_x_xxxz_yyyyyz,  \
                             tr_x_xxxz_yyyyz,   \
                             tr_x_xxxz_yyyyzz,  \
                             tr_x_xxxz_yyyzz,   \
                             tr_x_xxxz_yyyzzz,  \
                             tr_x_xxxz_yyzzz,   \
                             tr_x_xxxz_yyzzzz,  \
                             tr_x_xxxz_yzzzz,   \
                             tr_x_xxxz_yzzzzz,  \
                             tr_x_xxxz_zzzzz,   \
                             tr_x_xxxz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyz_xxxxxx[i] = tr_x_xxxz_xxxxxx[i] * pa_y[i];

        tr_x_xxxyz_xxxxxy[i] = tr_x_xxxy_xxxxxy[i] * pa_z[i];

        tr_x_xxxyz_xxxxxz[i] = tr_x_xxxz_xxxxxz[i] * pa_y[i];

        tr_x_xxxyz_xxxxyy[i] = tr_x_xxxy_xxxxyy[i] * pa_z[i];

        tr_x_xxxyz_xxxxyz[i] = tr_x_xxxz_xxxxz[i] * fe_0 + tr_x_xxxz_xxxxyz[i] * pa_y[i];

        tr_x_xxxyz_xxxxzz[i] = tr_x_xxxz_xxxxzz[i] * pa_y[i];

        tr_x_xxxyz_xxxyyy[i] = tr_x_xxxy_xxxyyy[i] * pa_z[i];

        tr_x_xxxyz_xxxyyz[i] = 2.0 * tr_x_xxxz_xxxyz[i] * fe_0 + tr_x_xxxz_xxxyyz[i] * pa_y[i];

        tr_x_xxxyz_xxxyzz[i] = tr_x_xxxz_xxxzz[i] * fe_0 + tr_x_xxxz_xxxyzz[i] * pa_y[i];

        tr_x_xxxyz_xxxzzz[i] = tr_x_xxxz_xxxzzz[i] * pa_y[i];

        tr_x_xxxyz_xxyyyy[i] = tr_x_xxxy_xxyyyy[i] * pa_z[i];

        tr_x_xxxyz_xxyyyz[i] = 3.0 * tr_x_xxxz_xxyyz[i] * fe_0 + tr_x_xxxz_xxyyyz[i] * pa_y[i];

        tr_x_xxxyz_xxyyzz[i] = 2.0 * tr_x_xxxz_xxyzz[i] * fe_0 + tr_x_xxxz_xxyyzz[i] * pa_y[i];

        tr_x_xxxyz_xxyzzz[i] = tr_x_xxxz_xxzzz[i] * fe_0 + tr_x_xxxz_xxyzzz[i] * pa_y[i];

        tr_x_xxxyz_xxzzzz[i] = tr_x_xxxz_xxzzzz[i] * pa_y[i];

        tr_x_xxxyz_xyyyyy[i] = tr_x_xxxy_xyyyyy[i] * pa_z[i];

        tr_x_xxxyz_xyyyyz[i] = 4.0 * tr_x_xxxz_xyyyz[i] * fe_0 + tr_x_xxxz_xyyyyz[i] * pa_y[i];

        tr_x_xxxyz_xyyyzz[i] = 3.0 * tr_x_xxxz_xyyzz[i] * fe_0 + tr_x_xxxz_xyyyzz[i] * pa_y[i];

        tr_x_xxxyz_xyyzzz[i] = 2.0 * tr_x_xxxz_xyzzz[i] * fe_0 + tr_x_xxxz_xyyzzz[i] * pa_y[i];

        tr_x_xxxyz_xyzzzz[i] = tr_x_xxxz_xzzzz[i] * fe_0 + tr_x_xxxz_xyzzzz[i] * pa_y[i];

        tr_x_xxxyz_xzzzzz[i] = tr_x_xxxz_xzzzzz[i] * pa_y[i];

        tr_x_xxxyz_yyyyyy[i] = tr_x_xxxy_yyyyyy[i] * pa_z[i];

        tr_x_xxxyz_yyyyyz[i] = 5.0 * tr_x_xxxz_yyyyz[i] * fe_0 + tr_x_xxxz_yyyyyz[i] * pa_y[i];

        tr_x_xxxyz_yyyyzz[i] = 4.0 * tr_x_xxxz_yyyzz[i] * fe_0 + tr_x_xxxz_yyyyzz[i] * pa_y[i];

        tr_x_xxxyz_yyyzzz[i] = 3.0 * tr_x_xxxz_yyzzz[i] * fe_0 + tr_x_xxxz_yyyzzz[i] * pa_y[i];

        tr_x_xxxyz_yyzzzz[i] = 2.0 * tr_x_xxxz_yzzzz[i] * fe_0 + tr_x_xxxz_yyzzzz[i] * pa_y[i];

        tr_x_xxxyz_yzzzzz[i] = tr_x_xxxz_zzzzz[i] * fe_0 + tr_x_xxxz_yzzzzz[i] * pa_y[i];

        tr_x_xxxyz_zzzzzz[i] = tr_x_xxxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : HI

    auto tr_x_xxxzz_xxxxxx = pbuffer.data(idx_dip_hi + 140);

    auto tr_x_xxxzz_xxxxxy = pbuffer.data(idx_dip_hi + 141);

    auto tr_x_xxxzz_xxxxxz = pbuffer.data(idx_dip_hi + 142);

    auto tr_x_xxxzz_xxxxyy = pbuffer.data(idx_dip_hi + 143);

    auto tr_x_xxxzz_xxxxyz = pbuffer.data(idx_dip_hi + 144);

    auto tr_x_xxxzz_xxxxzz = pbuffer.data(idx_dip_hi + 145);

    auto tr_x_xxxzz_xxxyyy = pbuffer.data(idx_dip_hi + 146);

    auto tr_x_xxxzz_xxxyyz = pbuffer.data(idx_dip_hi + 147);

    auto tr_x_xxxzz_xxxyzz = pbuffer.data(idx_dip_hi + 148);

    auto tr_x_xxxzz_xxxzzz = pbuffer.data(idx_dip_hi + 149);

    auto tr_x_xxxzz_xxyyyy = pbuffer.data(idx_dip_hi + 150);

    auto tr_x_xxxzz_xxyyyz = pbuffer.data(idx_dip_hi + 151);

    auto tr_x_xxxzz_xxyyzz = pbuffer.data(idx_dip_hi + 152);

    auto tr_x_xxxzz_xxyzzz = pbuffer.data(idx_dip_hi + 153);

    auto tr_x_xxxzz_xxzzzz = pbuffer.data(idx_dip_hi + 154);

    auto tr_x_xxxzz_xyyyyy = pbuffer.data(idx_dip_hi + 155);

    auto tr_x_xxxzz_xyyyyz = pbuffer.data(idx_dip_hi + 156);

    auto tr_x_xxxzz_xyyyzz = pbuffer.data(idx_dip_hi + 157);

    auto tr_x_xxxzz_xyyzzz = pbuffer.data(idx_dip_hi + 158);

    auto tr_x_xxxzz_xyzzzz = pbuffer.data(idx_dip_hi + 159);

    auto tr_x_xxxzz_xzzzzz = pbuffer.data(idx_dip_hi + 160);

    auto tr_x_xxxzz_yyyyyy = pbuffer.data(idx_dip_hi + 161);

    auto tr_x_xxxzz_yyyyyz = pbuffer.data(idx_dip_hi + 162);

    auto tr_x_xxxzz_yyyyzz = pbuffer.data(idx_dip_hi + 163);

    auto tr_x_xxxzz_yyyzzz = pbuffer.data(idx_dip_hi + 164);

    auto tr_x_xxxzz_yyzzzz = pbuffer.data(idx_dip_hi + 165);

    auto tr_x_xxxzz_yzzzzz = pbuffer.data(idx_dip_hi + 166);

    auto tr_x_xxxzz_zzzzzz = pbuffer.data(idx_dip_hi + 167);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_x_xxx_xxxxxx,   \
                             tr_x_xxx_xxxxxy,   \
                             tr_x_xxx_xxxxxz,   \
                             tr_x_xxx_xxxxyy,   \
                             tr_x_xxx_xxxxyz,   \
                             tr_x_xxx_xxxxzz,   \
                             tr_x_xxx_xxxyyy,   \
                             tr_x_xxx_xxxyyz,   \
                             tr_x_xxx_xxxyzz,   \
                             tr_x_xxx_xxxzzz,   \
                             tr_x_xxx_xxyyyy,   \
                             tr_x_xxx_xxyyyz,   \
                             tr_x_xxx_xxyyzz,   \
                             tr_x_xxx_xxyzzz,   \
                             tr_x_xxx_xxzzzz,   \
                             tr_x_xxx_xyyyyy,   \
                             tr_x_xxx_xyyyyz,   \
                             tr_x_xxx_xyyyzz,   \
                             tr_x_xxx_xyyzzz,   \
                             tr_x_xxx_xyzzzz,   \
                             tr_x_xxx_xzzzzz,   \
                             tr_x_xxx_yyyyyy,   \
                             tr_x_xxxz_xxxxx,   \
                             tr_x_xxxz_xxxxxx,  \
                             tr_x_xxxz_xxxxxy,  \
                             tr_x_xxxz_xxxxxz,  \
                             tr_x_xxxz_xxxxy,   \
                             tr_x_xxxz_xxxxyy,  \
                             tr_x_xxxz_xxxxyz,  \
                             tr_x_xxxz_xxxxz,   \
                             tr_x_xxxz_xxxxzz,  \
                             tr_x_xxxz_xxxyy,   \
                             tr_x_xxxz_xxxyyy,  \
                             tr_x_xxxz_xxxyyz,  \
                             tr_x_xxxz_xxxyz,   \
                             tr_x_xxxz_xxxyzz,  \
                             tr_x_xxxz_xxxzz,   \
                             tr_x_xxxz_xxxzzz,  \
                             tr_x_xxxz_xxyyy,   \
                             tr_x_xxxz_xxyyyy,  \
                             tr_x_xxxz_xxyyyz,  \
                             tr_x_xxxz_xxyyz,   \
                             tr_x_xxxz_xxyyzz,  \
                             tr_x_xxxz_xxyzz,   \
                             tr_x_xxxz_xxyzzz,  \
                             tr_x_xxxz_xxzzz,   \
                             tr_x_xxxz_xxzzzz,  \
                             tr_x_xxxz_xyyyy,   \
                             tr_x_xxxz_xyyyyy,  \
                             tr_x_xxxz_xyyyyz,  \
                             tr_x_xxxz_xyyyz,   \
                             tr_x_xxxz_xyyyzz,  \
                             tr_x_xxxz_xyyzz,   \
                             tr_x_xxxz_xyyzzz,  \
                             tr_x_xxxz_xyzzz,   \
                             tr_x_xxxz_xyzzzz,  \
                             tr_x_xxxz_xzzzz,   \
                             tr_x_xxxz_xzzzzz,  \
                             tr_x_xxxz_yyyyyy,  \
                             tr_x_xxxzz_xxxxxx, \
                             tr_x_xxxzz_xxxxxy, \
                             tr_x_xxxzz_xxxxxz, \
                             tr_x_xxxzz_xxxxyy, \
                             tr_x_xxxzz_xxxxyz, \
                             tr_x_xxxzz_xxxxzz, \
                             tr_x_xxxzz_xxxyyy, \
                             tr_x_xxxzz_xxxyyz, \
                             tr_x_xxxzz_xxxyzz, \
                             tr_x_xxxzz_xxxzzz, \
                             tr_x_xxxzz_xxyyyy, \
                             tr_x_xxxzz_xxyyyz, \
                             tr_x_xxxzz_xxyyzz, \
                             tr_x_xxxzz_xxyzzz, \
                             tr_x_xxxzz_xxzzzz, \
                             tr_x_xxxzz_xyyyyy, \
                             tr_x_xxxzz_xyyyyz, \
                             tr_x_xxxzz_xyyyzz, \
                             tr_x_xxxzz_xyyzzz, \
                             tr_x_xxxzz_xyzzzz, \
                             tr_x_xxxzz_xzzzzz, \
                             tr_x_xxxzz_yyyyyy, \
                             tr_x_xxxzz_yyyyyz, \
                             tr_x_xxxzz_yyyyzz, \
                             tr_x_xxxzz_yyyzzz, \
                             tr_x_xxxzz_yyzzzz, \
                             tr_x_xxxzz_yzzzzz, \
                             tr_x_xxxzz_zzzzzz, \
                             tr_x_xxzz_yyyyyz,  \
                             tr_x_xxzz_yyyyzz,  \
                             tr_x_xxzz_yyyzzz,  \
                             tr_x_xxzz_yyzzzz,  \
                             tr_x_xxzz_yzzzzz,  \
                             tr_x_xxzz_zzzzzz,  \
                             tr_x_xzz_yyyyyz,   \
                             tr_x_xzz_yyyyzz,   \
                             tr_x_xzz_yyyzzz,   \
                             tr_x_xzz_yyzzzz,   \
                             tr_x_xzz_yzzzzz,   \
                             tr_x_xzz_zzzzzz,   \
                             ts_xxzz_yyyyyz,    \
                             ts_xxzz_yyyyzz,    \
                             ts_xxzz_yyyzzz,    \
                             ts_xxzz_yyzzzz,    \
                             ts_xxzz_yzzzzz,    \
                             ts_xxzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzz_xxxxxx[i] = tr_x_xxx_xxxxxx[i] * fe_0 + tr_x_xxxz_xxxxxx[i] * pa_z[i];

        tr_x_xxxzz_xxxxxy[i] = tr_x_xxx_xxxxxy[i] * fe_0 + tr_x_xxxz_xxxxxy[i] * pa_z[i];

        tr_x_xxxzz_xxxxxz[i] = tr_x_xxx_xxxxxz[i] * fe_0 + tr_x_xxxz_xxxxx[i] * fe_0 + tr_x_xxxz_xxxxxz[i] * pa_z[i];

        tr_x_xxxzz_xxxxyy[i] = tr_x_xxx_xxxxyy[i] * fe_0 + tr_x_xxxz_xxxxyy[i] * pa_z[i];

        tr_x_xxxzz_xxxxyz[i] = tr_x_xxx_xxxxyz[i] * fe_0 + tr_x_xxxz_xxxxy[i] * fe_0 + tr_x_xxxz_xxxxyz[i] * pa_z[i];

        tr_x_xxxzz_xxxxzz[i] = tr_x_xxx_xxxxzz[i] * fe_0 + 2.0 * tr_x_xxxz_xxxxz[i] * fe_0 + tr_x_xxxz_xxxxzz[i] * pa_z[i];

        tr_x_xxxzz_xxxyyy[i] = tr_x_xxx_xxxyyy[i] * fe_0 + tr_x_xxxz_xxxyyy[i] * pa_z[i];

        tr_x_xxxzz_xxxyyz[i] = tr_x_xxx_xxxyyz[i] * fe_0 + tr_x_xxxz_xxxyy[i] * fe_0 + tr_x_xxxz_xxxyyz[i] * pa_z[i];

        tr_x_xxxzz_xxxyzz[i] = tr_x_xxx_xxxyzz[i] * fe_0 + 2.0 * tr_x_xxxz_xxxyz[i] * fe_0 + tr_x_xxxz_xxxyzz[i] * pa_z[i];

        tr_x_xxxzz_xxxzzz[i] = tr_x_xxx_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxxz_xxxzz[i] * fe_0 + tr_x_xxxz_xxxzzz[i] * pa_z[i];

        tr_x_xxxzz_xxyyyy[i] = tr_x_xxx_xxyyyy[i] * fe_0 + tr_x_xxxz_xxyyyy[i] * pa_z[i];

        tr_x_xxxzz_xxyyyz[i] = tr_x_xxx_xxyyyz[i] * fe_0 + tr_x_xxxz_xxyyy[i] * fe_0 + tr_x_xxxz_xxyyyz[i] * pa_z[i];

        tr_x_xxxzz_xxyyzz[i] = tr_x_xxx_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxxz_xxyyz[i] * fe_0 + tr_x_xxxz_xxyyzz[i] * pa_z[i];

        tr_x_xxxzz_xxyzzz[i] = tr_x_xxx_xxyzzz[i] * fe_0 + 3.0 * tr_x_xxxz_xxyzz[i] * fe_0 + tr_x_xxxz_xxyzzz[i] * pa_z[i];

        tr_x_xxxzz_xxzzzz[i] = tr_x_xxx_xxzzzz[i] * fe_0 + 4.0 * tr_x_xxxz_xxzzz[i] * fe_0 + tr_x_xxxz_xxzzzz[i] * pa_z[i];

        tr_x_xxxzz_xyyyyy[i] = tr_x_xxx_xyyyyy[i] * fe_0 + tr_x_xxxz_xyyyyy[i] * pa_z[i];

        tr_x_xxxzz_xyyyyz[i] = tr_x_xxx_xyyyyz[i] * fe_0 + tr_x_xxxz_xyyyy[i] * fe_0 + tr_x_xxxz_xyyyyz[i] * pa_z[i];

        tr_x_xxxzz_xyyyzz[i] = tr_x_xxx_xyyyzz[i] * fe_0 + 2.0 * tr_x_xxxz_xyyyz[i] * fe_0 + tr_x_xxxz_xyyyzz[i] * pa_z[i];

        tr_x_xxxzz_xyyzzz[i] = tr_x_xxx_xyyzzz[i] * fe_0 + 3.0 * tr_x_xxxz_xyyzz[i] * fe_0 + tr_x_xxxz_xyyzzz[i] * pa_z[i];

        tr_x_xxxzz_xyzzzz[i] = tr_x_xxx_xyzzzz[i] * fe_0 + 4.0 * tr_x_xxxz_xyzzz[i] * fe_0 + tr_x_xxxz_xyzzzz[i] * pa_z[i];

        tr_x_xxxzz_xzzzzz[i] = tr_x_xxx_xzzzzz[i] * fe_0 + 5.0 * tr_x_xxxz_xzzzz[i] * fe_0 + tr_x_xxxz_xzzzzz[i] * pa_z[i];

        tr_x_xxxzz_yyyyyy[i] = tr_x_xxx_yyyyyy[i] * fe_0 + tr_x_xxxz_yyyyyy[i] * pa_z[i];

        tr_x_xxxzz_yyyyyz[i] = 2.0 * tr_x_xzz_yyyyyz[i] * fe_0 + ts_xxzz_yyyyyz[i] * fe_0 + tr_x_xxzz_yyyyyz[i] * pa_x[i];

        tr_x_xxxzz_yyyyzz[i] = 2.0 * tr_x_xzz_yyyyzz[i] * fe_0 + ts_xxzz_yyyyzz[i] * fe_0 + tr_x_xxzz_yyyyzz[i] * pa_x[i];

        tr_x_xxxzz_yyyzzz[i] = 2.0 * tr_x_xzz_yyyzzz[i] * fe_0 + ts_xxzz_yyyzzz[i] * fe_0 + tr_x_xxzz_yyyzzz[i] * pa_x[i];

        tr_x_xxxzz_yyzzzz[i] = 2.0 * tr_x_xzz_yyzzzz[i] * fe_0 + ts_xxzz_yyzzzz[i] * fe_0 + tr_x_xxzz_yyzzzz[i] * pa_x[i];

        tr_x_xxxzz_yzzzzz[i] = 2.0 * tr_x_xzz_yzzzzz[i] * fe_0 + ts_xxzz_yzzzzz[i] * fe_0 + tr_x_xxzz_yzzzzz[i] * pa_x[i];

        tr_x_xxxzz_zzzzzz[i] = 2.0 * tr_x_xzz_zzzzzz[i] * fe_0 + ts_xxzz_zzzzzz[i] * fe_0 + tr_x_xxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : HI

    auto tr_x_xxyyy_xxxxxx = pbuffer.data(idx_dip_hi + 168);

    auto tr_x_xxyyy_xxxxxy = pbuffer.data(idx_dip_hi + 169);

    auto tr_x_xxyyy_xxxxxz = pbuffer.data(idx_dip_hi + 170);

    auto tr_x_xxyyy_xxxxyy = pbuffer.data(idx_dip_hi + 171);

    auto tr_x_xxyyy_xxxxyz = pbuffer.data(idx_dip_hi + 172);

    auto tr_x_xxyyy_xxxxzz = pbuffer.data(idx_dip_hi + 173);

    auto tr_x_xxyyy_xxxyyy = pbuffer.data(idx_dip_hi + 174);

    auto tr_x_xxyyy_xxxyyz = pbuffer.data(idx_dip_hi + 175);

    auto tr_x_xxyyy_xxxyzz = pbuffer.data(idx_dip_hi + 176);

    auto tr_x_xxyyy_xxxzzz = pbuffer.data(idx_dip_hi + 177);

    auto tr_x_xxyyy_xxyyyy = pbuffer.data(idx_dip_hi + 178);

    auto tr_x_xxyyy_xxyyyz = pbuffer.data(idx_dip_hi + 179);

    auto tr_x_xxyyy_xxyyzz = pbuffer.data(idx_dip_hi + 180);

    auto tr_x_xxyyy_xxyzzz = pbuffer.data(idx_dip_hi + 181);

    auto tr_x_xxyyy_xxzzzz = pbuffer.data(idx_dip_hi + 182);

    auto tr_x_xxyyy_xyyyyy = pbuffer.data(idx_dip_hi + 183);

    auto tr_x_xxyyy_xyyyyz = pbuffer.data(idx_dip_hi + 184);

    auto tr_x_xxyyy_xyyyzz = pbuffer.data(idx_dip_hi + 185);

    auto tr_x_xxyyy_xyyzzz = pbuffer.data(idx_dip_hi + 186);

    auto tr_x_xxyyy_xyzzzz = pbuffer.data(idx_dip_hi + 187);

    auto tr_x_xxyyy_xzzzzz = pbuffer.data(idx_dip_hi + 188);

    auto tr_x_xxyyy_yyyyyy = pbuffer.data(idx_dip_hi + 189);

    auto tr_x_xxyyy_yyyyyz = pbuffer.data(idx_dip_hi + 190);

    auto tr_x_xxyyy_yyyyzz = pbuffer.data(idx_dip_hi + 191);

    auto tr_x_xxyyy_yyyzzz = pbuffer.data(idx_dip_hi + 192);

    auto tr_x_xxyyy_yyzzzz = pbuffer.data(idx_dip_hi + 193);

    auto tr_x_xxyyy_yzzzzz = pbuffer.data(idx_dip_hi + 194);

    auto tr_x_xxyyy_zzzzzz = pbuffer.data(idx_dip_hi + 195);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_x_xxy_xxxxxx,   \
                             tr_x_xxy_xxxxxy,   \
                             tr_x_xxy_xxxxxz,   \
                             tr_x_xxy_xxxxyy,   \
                             tr_x_xxy_xxxxyz,   \
                             tr_x_xxy_xxxxzz,   \
                             tr_x_xxy_xxxyyy,   \
                             tr_x_xxy_xxxyyz,   \
                             tr_x_xxy_xxxyzz,   \
                             tr_x_xxy_xxxzzz,   \
                             tr_x_xxy_xxyyyy,   \
                             tr_x_xxy_xxyyyz,   \
                             tr_x_xxy_xxyyzz,   \
                             tr_x_xxy_xxyzzz,   \
                             tr_x_xxy_xxzzzz,   \
                             tr_x_xxy_xyyyyy,   \
                             tr_x_xxy_xyyyyz,   \
                             tr_x_xxy_xyyyzz,   \
                             tr_x_xxy_xyyzzz,   \
                             tr_x_xxy_xyzzzz,   \
                             tr_x_xxy_xzzzzz,   \
                             tr_x_xxy_zzzzzz,   \
                             tr_x_xxyy_xxxxx,   \
                             tr_x_xxyy_xxxxxx,  \
                             tr_x_xxyy_xxxxxy,  \
                             tr_x_xxyy_xxxxxz,  \
                             tr_x_xxyy_xxxxy,   \
                             tr_x_xxyy_xxxxyy,  \
                             tr_x_xxyy_xxxxyz,  \
                             tr_x_xxyy_xxxxz,   \
                             tr_x_xxyy_xxxxzz,  \
                             tr_x_xxyy_xxxyy,   \
                             tr_x_xxyy_xxxyyy,  \
                             tr_x_xxyy_xxxyyz,  \
                             tr_x_xxyy_xxxyz,   \
                             tr_x_xxyy_xxxyzz,  \
                             tr_x_xxyy_xxxzz,   \
                             tr_x_xxyy_xxxzzz,  \
                             tr_x_xxyy_xxyyy,   \
                             tr_x_xxyy_xxyyyy,  \
                             tr_x_xxyy_xxyyyz,  \
                             tr_x_xxyy_xxyyz,   \
                             tr_x_xxyy_xxyyzz,  \
                             tr_x_xxyy_xxyzz,   \
                             tr_x_xxyy_xxyzzz,  \
                             tr_x_xxyy_xxzzz,   \
                             tr_x_xxyy_xxzzzz,  \
                             tr_x_xxyy_xyyyy,   \
                             tr_x_xxyy_xyyyyy,  \
                             tr_x_xxyy_xyyyyz,  \
                             tr_x_xxyy_xyyyz,   \
                             tr_x_xxyy_xyyyzz,  \
                             tr_x_xxyy_xyyzz,   \
                             tr_x_xxyy_xyyzzz,  \
                             tr_x_xxyy_xyzzz,   \
                             tr_x_xxyy_xyzzzz,  \
                             tr_x_xxyy_xzzzz,   \
                             tr_x_xxyy_xzzzzz,  \
                             tr_x_xxyy_zzzzzz,  \
                             tr_x_xxyyy_xxxxxx, \
                             tr_x_xxyyy_xxxxxy, \
                             tr_x_xxyyy_xxxxxz, \
                             tr_x_xxyyy_xxxxyy, \
                             tr_x_xxyyy_xxxxyz, \
                             tr_x_xxyyy_xxxxzz, \
                             tr_x_xxyyy_xxxyyy, \
                             tr_x_xxyyy_xxxyyz, \
                             tr_x_xxyyy_xxxyzz, \
                             tr_x_xxyyy_xxxzzz, \
                             tr_x_xxyyy_xxyyyy, \
                             tr_x_xxyyy_xxyyyz, \
                             tr_x_xxyyy_xxyyzz, \
                             tr_x_xxyyy_xxyzzz, \
                             tr_x_xxyyy_xxzzzz, \
                             tr_x_xxyyy_xyyyyy, \
                             tr_x_xxyyy_xyyyyz, \
                             tr_x_xxyyy_xyyyzz, \
                             tr_x_xxyyy_xyyzzz, \
                             tr_x_xxyyy_xyzzzz, \
                             tr_x_xxyyy_xzzzzz, \
                             tr_x_xxyyy_yyyyyy, \
                             tr_x_xxyyy_yyyyyz, \
                             tr_x_xxyyy_yyyyzz, \
                             tr_x_xxyyy_yyyzzz, \
                             tr_x_xxyyy_yyzzzz, \
                             tr_x_xxyyy_yzzzzz, \
                             tr_x_xxyyy_zzzzzz, \
                             tr_x_xyyy_yyyyyy,  \
                             tr_x_xyyy_yyyyyz,  \
                             tr_x_xyyy_yyyyzz,  \
                             tr_x_xyyy_yyyzzz,  \
                             tr_x_xyyy_yyzzzz,  \
                             tr_x_xyyy_yzzzzz,  \
                             tr_x_yyy_yyyyyy,   \
                             tr_x_yyy_yyyyyz,   \
                             tr_x_yyy_yyyyzz,   \
                             tr_x_yyy_yyyzzz,   \
                             tr_x_yyy_yyzzzz,   \
                             tr_x_yyy_yzzzzz,   \
                             ts_xyyy_yyyyyy,    \
                             ts_xyyy_yyyyyz,    \
                             ts_xyyy_yyyyzz,    \
                             ts_xyyy_yyyzzz,    \
                             ts_xyyy_yyzzzz,    \
                             ts_xyyy_yzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyy_xxxxxx[i] = 2.0 * tr_x_xxy_xxxxxx[i] * fe_0 + tr_x_xxyy_xxxxxx[i] * pa_y[i];

        tr_x_xxyyy_xxxxxy[i] = 2.0 * tr_x_xxy_xxxxxy[i] * fe_0 + tr_x_xxyy_xxxxx[i] * fe_0 + tr_x_xxyy_xxxxxy[i] * pa_y[i];

        tr_x_xxyyy_xxxxxz[i] = 2.0 * tr_x_xxy_xxxxxz[i] * fe_0 + tr_x_xxyy_xxxxxz[i] * pa_y[i];

        tr_x_xxyyy_xxxxyy[i] = 2.0 * tr_x_xxy_xxxxyy[i] * fe_0 + 2.0 * tr_x_xxyy_xxxxy[i] * fe_0 + tr_x_xxyy_xxxxyy[i] * pa_y[i];

        tr_x_xxyyy_xxxxyz[i] = 2.0 * tr_x_xxy_xxxxyz[i] * fe_0 + tr_x_xxyy_xxxxz[i] * fe_0 + tr_x_xxyy_xxxxyz[i] * pa_y[i];

        tr_x_xxyyy_xxxxzz[i] = 2.0 * tr_x_xxy_xxxxzz[i] * fe_0 + tr_x_xxyy_xxxxzz[i] * pa_y[i];

        tr_x_xxyyy_xxxyyy[i] = 2.0 * tr_x_xxy_xxxyyy[i] * fe_0 + 3.0 * tr_x_xxyy_xxxyy[i] * fe_0 + tr_x_xxyy_xxxyyy[i] * pa_y[i];

        tr_x_xxyyy_xxxyyz[i] = 2.0 * tr_x_xxy_xxxyyz[i] * fe_0 + 2.0 * tr_x_xxyy_xxxyz[i] * fe_0 + tr_x_xxyy_xxxyyz[i] * pa_y[i];

        tr_x_xxyyy_xxxyzz[i] = 2.0 * tr_x_xxy_xxxyzz[i] * fe_0 + tr_x_xxyy_xxxzz[i] * fe_0 + tr_x_xxyy_xxxyzz[i] * pa_y[i];

        tr_x_xxyyy_xxxzzz[i] = 2.0 * tr_x_xxy_xxxzzz[i] * fe_0 + tr_x_xxyy_xxxzzz[i] * pa_y[i];

        tr_x_xxyyy_xxyyyy[i] = 2.0 * tr_x_xxy_xxyyyy[i] * fe_0 + 4.0 * tr_x_xxyy_xxyyy[i] * fe_0 + tr_x_xxyy_xxyyyy[i] * pa_y[i];

        tr_x_xxyyy_xxyyyz[i] = 2.0 * tr_x_xxy_xxyyyz[i] * fe_0 + 3.0 * tr_x_xxyy_xxyyz[i] * fe_0 + tr_x_xxyy_xxyyyz[i] * pa_y[i];

        tr_x_xxyyy_xxyyzz[i] = 2.0 * tr_x_xxy_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxyy_xxyzz[i] * fe_0 + tr_x_xxyy_xxyyzz[i] * pa_y[i];

        tr_x_xxyyy_xxyzzz[i] = 2.0 * tr_x_xxy_xxyzzz[i] * fe_0 + tr_x_xxyy_xxzzz[i] * fe_0 + tr_x_xxyy_xxyzzz[i] * pa_y[i];

        tr_x_xxyyy_xxzzzz[i] = 2.0 * tr_x_xxy_xxzzzz[i] * fe_0 + tr_x_xxyy_xxzzzz[i] * pa_y[i];

        tr_x_xxyyy_xyyyyy[i] = 2.0 * tr_x_xxy_xyyyyy[i] * fe_0 + 5.0 * tr_x_xxyy_xyyyy[i] * fe_0 + tr_x_xxyy_xyyyyy[i] * pa_y[i];

        tr_x_xxyyy_xyyyyz[i] = 2.0 * tr_x_xxy_xyyyyz[i] * fe_0 + 4.0 * tr_x_xxyy_xyyyz[i] * fe_0 + tr_x_xxyy_xyyyyz[i] * pa_y[i];

        tr_x_xxyyy_xyyyzz[i] = 2.0 * tr_x_xxy_xyyyzz[i] * fe_0 + 3.0 * tr_x_xxyy_xyyzz[i] * fe_0 + tr_x_xxyy_xyyyzz[i] * pa_y[i];

        tr_x_xxyyy_xyyzzz[i] = 2.0 * tr_x_xxy_xyyzzz[i] * fe_0 + 2.0 * tr_x_xxyy_xyzzz[i] * fe_0 + tr_x_xxyy_xyyzzz[i] * pa_y[i];

        tr_x_xxyyy_xyzzzz[i] = 2.0 * tr_x_xxy_xyzzzz[i] * fe_0 + tr_x_xxyy_xzzzz[i] * fe_0 + tr_x_xxyy_xyzzzz[i] * pa_y[i];

        tr_x_xxyyy_xzzzzz[i] = 2.0 * tr_x_xxy_xzzzzz[i] * fe_0 + tr_x_xxyy_xzzzzz[i] * pa_y[i];

        tr_x_xxyyy_yyyyyy[i] = tr_x_yyy_yyyyyy[i] * fe_0 + ts_xyyy_yyyyyy[i] * fe_0 + tr_x_xyyy_yyyyyy[i] * pa_x[i];

        tr_x_xxyyy_yyyyyz[i] = tr_x_yyy_yyyyyz[i] * fe_0 + ts_xyyy_yyyyyz[i] * fe_0 + tr_x_xyyy_yyyyyz[i] * pa_x[i];

        tr_x_xxyyy_yyyyzz[i] = tr_x_yyy_yyyyzz[i] * fe_0 + ts_xyyy_yyyyzz[i] * fe_0 + tr_x_xyyy_yyyyzz[i] * pa_x[i];

        tr_x_xxyyy_yyyzzz[i] = tr_x_yyy_yyyzzz[i] * fe_0 + ts_xyyy_yyyzzz[i] * fe_0 + tr_x_xyyy_yyyzzz[i] * pa_x[i];

        tr_x_xxyyy_yyzzzz[i] = tr_x_yyy_yyzzzz[i] * fe_0 + ts_xyyy_yyzzzz[i] * fe_0 + tr_x_xyyy_yyzzzz[i] * pa_x[i];

        tr_x_xxyyy_yzzzzz[i] = tr_x_yyy_yzzzzz[i] * fe_0 + ts_xyyy_yzzzzz[i] * fe_0 + tr_x_xyyy_yzzzzz[i] * pa_x[i];

        tr_x_xxyyy_zzzzzz[i] = 2.0 * tr_x_xxy_zzzzzz[i] * fe_0 + tr_x_xxyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 196-224 components of targeted buffer : HI

    auto tr_x_xxyyz_xxxxxx = pbuffer.data(idx_dip_hi + 196);

    auto tr_x_xxyyz_xxxxxy = pbuffer.data(idx_dip_hi + 197);

    auto tr_x_xxyyz_xxxxxz = pbuffer.data(idx_dip_hi + 198);

    auto tr_x_xxyyz_xxxxyy = pbuffer.data(idx_dip_hi + 199);

    auto tr_x_xxyyz_xxxxyz = pbuffer.data(idx_dip_hi + 200);

    auto tr_x_xxyyz_xxxxzz = pbuffer.data(idx_dip_hi + 201);

    auto tr_x_xxyyz_xxxyyy = pbuffer.data(idx_dip_hi + 202);

    auto tr_x_xxyyz_xxxyyz = pbuffer.data(idx_dip_hi + 203);

    auto tr_x_xxyyz_xxxyzz = pbuffer.data(idx_dip_hi + 204);

    auto tr_x_xxyyz_xxxzzz = pbuffer.data(idx_dip_hi + 205);

    auto tr_x_xxyyz_xxyyyy = pbuffer.data(idx_dip_hi + 206);

    auto tr_x_xxyyz_xxyyyz = pbuffer.data(idx_dip_hi + 207);

    auto tr_x_xxyyz_xxyyzz = pbuffer.data(idx_dip_hi + 208);

    auto tr_x_xxyyz_xxyzzz = pbuffer.data(idx_dip_hi + 209);

    auto tr_x_xxyyz_xxzzzz = pbuffer.data(idx_dip_hi + 210);

    auto tr_x_xxyyz_xyyyyy = pbuffer.data(idx_dip_hi + 211);

    auto tr_x_xxyyz_xyyyyz = pbuffer.data(idx_dip_hi + 212);

    auto tr_x_xxyyz_xyyyzz = pbuffer.data(idx_dip_hi + 213);

    auto tr_x_xxyyz_xyyzzz = pbuffer.data(idx_dip_hi + 214);

    auto tr_x_xxyyz_xyzzzz = pbuffer.data(idx_dip_hi + 215);

    auto tr_x_xxyyz_xzzzzz = pbuffer.data(idx_dip_hi + 216);

    auto tr_x_xxyyz_yyyyyy = pbuffer.data(idx_dip_hi + 217);

    auto tr_x_xxyyz_yyyyyz = pbuffer.data(idx_dip_hi + 218);

    auto tr_x_xxyyz_yyyyzz = pbuffer.data(idx_dip_hi + 219);

    auto tr_x_xxyyz_yyyzzz = pbuffer.data(idx_dip_hi + 220);

    auto tr_x_xxyyz_yyzzzz = pbuffer.data(idx_dip_hi + 221);

    auto tr_x_xxyyz_yzzzzz = pbuffer.data(idx_dip_hi + 222);

    auto tr_x_xxyyz_zzzzzz = pbuffer.data(idx_dip_hi + 223);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_x_xxyy_xxxxxx,  \
                             tr_x_xxyy_xxxxxy,  \
                             tr_x_xxyy_xxxxy,   \
                             tr_x_xxyy_xxxxyy,  \
                             tr_x_xxyy_xxxxyz,  \
                             tr_x_xxyy_xxxyy,   \
                             tr_x_xxyy_xxxyyy,  \
                             tr_x_xxyy_xxxyyz,  \
                             tr_x_xxyy_xxxyz,   \
                             tr_x_xxyy_xxxyzz,  \
                             tr_x_xxyy_xxyyy,   \
                             tr_x_xxyy_xxyyyy,  \
                             tr_x_xxyy_xxyyyz,  \
                             tr_x_xxyy_xxyyz,   \
                             tr_x_xxyy_xxyyzz,  \
                             tr_x_xxyy_xxyzz,   \
                             tr_x_xxyy_xxyzzz,  \
                             tr_x_xxyy_xyyyy,   \
                             tr_x_xxyy_xyyyyy,  \
                             tr_x_xxyy_xyyyyz,  \
                             tr_x_xxyy_xyyyz,   \
                             tr_x_xxyy_xyyyzz,  \
                             tr_x_xxyy_xyyzz,   \
                             tr_x_xxyy_xyyzzz,  \
                             tr_x_xxyy_xyzzz,   \
                             tr_x_xxyy_xyzzzz,  \
                             tr_x_xxyy_yyyyy,   \
                             tr_x_xxyy_yyyyyy,  \
                             tr_x_xxyy_yyyyyz,  \
                             tr_x_xxyy_yyyyz,   \
                             tr_x_xxyy_yyyyzz,  \
                             tr_x_xxyy_yyyzz,   \
                             tr_x_xxyy_yyyzzz,  \
                             tr_x_xxyy_yyzzz,   \
                             tr_x_xxyy_yyzzzz,  \
                             tr_x_xxyy_yzzzz,   \
                             tr_x_xxyy_yzzzzz,  \
                             tr_x_xxyyz_xxxxxx, \
                             tr_x_xxyyz_xxxxxy, \
                             tr_x_xxyyz_xxxxxz, \
                             tr_x_xxyyz_xxxxyy, \
                             tr_x_xxyyz_xxxxyz, \
                             tr_x_xxyyz_xxxxzz, \
                             tr_x_xxyyz_xxxyyy, \
                             tr_x_xxyyz_xxxyyz, \
                             tr_x_xxyyz_xxxyzz, \
                             tr_x_xxyyz_xxxzzz, \
                             tr_x_xxyyz_xxyyyy, \
                             tr_x_xxyyz_xxyyyz, \
                             tr_x_xxyyz_xxyyzz, \
                             tr_x_xxyyz_xxyzzz, \
                             tr_x_xxyyz_xxzzzz, \
                             tr_x_xxyyz_xyyyyy, \
                             tr_x_xxyyz_xyyyyz, \
                             tr_x_xxyyz_xyyyzz, \
                             tr_x_xxyyz_xyyzzz, \
                             tr_x_xxyyz_xyzzzz, \
                             tr_x_xxyyz_xzzzzz, \
                             tr_x_xxyyz_yyyyyy, \
                             tr_x_xxyyz_yyyyyz, \
                             tr_x_xxyyz_yyyyzz, \
                             tr_x_xxyyz_yyyzzz, \
                             tr_x_xxyyz_yyzzzz, \
                             tr_x_xxyyz_yzzzzz, \
                             tr_x_xxyyz_zzzzzz, \
                             tr_x_xxyz_xxxxxz,  \
                             tr_x_xxyz_xxxxzz,  \
                             tr_x_xxyz_xxxzzz,  \
                             tr_x_xxyz_xxzzzz,  \
                             tr_x_xxyz_xzzzzz,  \
                             tr_x_xxyz_zzzzzz,  \
                             tr_x_xxz_xxxxxz,   \
                             tr_x_xxz_xxxxzz,   \
                             tr_x_xxz_xxxzzz,   \
                             tr_x_xxz_xxzzzz,   \
                             tr_x_xxz_xzzzzz,   \
                             tr_x_xxz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyz_xxxxxx[i] = tr_x_xxyy_xxxxxx[i] * pa_z[i];

        tr_x_xxyyz_xxxxxy[i] = tr_x_xxyy_xxxxxy[i] * pa_z[i];

        tr_x_xxyyz_xxxxxz[i] = tr_x_xxz_xxxxxz[i] * fe_0 + tr_x_xxyz_xxxxxz[i] * pa_y[i];

        tr_x_xxyyz_xxxxyy[i] = tr_x_xxyy_xxxxyy[i] * pa_z[i];

        tr_x_xxyyz_xxxxyz[i] = tr_x_xxyy_xxxxy[i] * fe_0 + tr_x_xxyy_xxxxyz[i] * pa_z[i];

        tr_x_xxyyz_xxxxzz[i] = tr_x_xxz_xxxxzz[i] * fe_0 + tr_x_xxyz_xxxxzz[i] * pa_y[i];

        tr_x_xxyyz_xxxyyy[i] = tr_x_xxyy_xxxyyy[i] * pa_z[i];

        tr_x_xxyyz_xxxyyz[i] = tr_x_xxyy_xxxyy[i] * fe_0 + tr_x_xxyy_xxxyyz[i] * pa_z[i];

        tr_x_xxyyz_xxxyzz[i] = 2.0 * tr_x_xxyy_xxxyz[i] * fe_0 + tr_x_xxyy_xxxyzz[i] * pa_z[i];

        tr_x_xxyyz_xxxzzz[i] = tr_x_xxz_xxxzzz[i] * fe_0 + tr_x_xxyz_xxxzzz[i] * pa_y[i];

        tr_x_xxyyz_xxyyyy[i] = tr_x_xxyy_xxyyyy[i] * pa_z[i];

        tr_x_xxyyz_xxyyyz[i] = tr_x_xxyy_xxyyy[i] * fe_0 + tr_x_xxyy_xxyyyz[i] * pa_z[i];

        tr_x_xxyyz_xxyyzz[i] = 2.0 * tr_x_xxyy_xxyyz[i] * fe_0 + tr_x_xxyy_xxyyzz[i] * pa_z[i];

        tr_x_xxyyz_xxyzzz[i] = 3.0 * tr_x_xxyy_xxyzz[i] * fe_0 + tr_x_xxyy_xxyzzz[i] * pa_z[i];

        tr_x_xxyyz_xxzzzz[i] = tr_x_xxz_xxzzzz[i] * fe_0 + tr_x_xxyz_xxzzzz[i] * pa_y[i];

        tr_x_xxyyz_xyyyyy[i] = tr_x_xxyy_xyyyyy[i] * pa_z[i];

        tr_x_xxyyz_xyyyyz[i] = tr_x_xxyy_xyyyy[i] * fe_0 + tr_x_xxyy_xyyyyz[i] * pa_z[i];

        tr_x_xxyyz_xyyyzz[i] = 2.0 * tr_x_xxyy_xyyyz[i] * fe_0 + tr_x_xxyy_xyyyzz[i] * pa_z[i];

        tr_x_xxyyz_xyyzzz[i] = 3.0 * tr_x_xxyy_xyyzz[i] * fe_0 + tr_x_xxyy_xyyzzz[i] * pa_z[i];

        tr_x_xxyyz_xyzzzz[i] = 4.0 * tr_x_xxyy_xyzzz[i] * fe_0 + tr_x_xxyy_xyzzzz[i] * pa_z[i];

        tr_x_xxyyz_xzzzzz[i] = tr_x_xxz_xzzzzz[i] * fe_0 + tr_x_xxyz_xzzzzz[i] * pa_y[i];

        tr_x_xxyyz_yyyyyy[i] = tr_x_xxyy_yyyyyy[i] * pa_z[i];

        tr_x_xxyyz_yyyyyz[i] = tr_x_xxyy_yyyyy[i] * fe_0 + tr_x_xxyy_yyyyyz[i] * pa_z[i];

        tr_x_xxyyz_yyyyzz[i] = 2.0 * tr_x_xxyy_yyyyz[i] * fe_0 + tr_x_xxyy_yyyyzz[i] * pa_z[i];

        tr_x_xxyyz_yyyzzz[i] = 3.0 * tr_x_xxyy_yyyzz[i] * fe_0 + tr_x_xxyy_yyyzzz[i] * pa_z[i];

        tr_x_xxyyz_yyzzzz[i] = 4.0 * tr_x_xxyy_yyzzz[i] * fe_0 + tr_x_xxyy_yyzzzz[i] * pa_z[i];

        tr_x_xxyyz_yzzzzz[i] = 5.0 * tr_x_xxyy_yzzzz[i] * fe_0 + tr_x_xxyy_yzzzzz[i] * pa_z[i];

        tr_x_xxyyz_zzzzzz[i] = tr_x_xxz_zzzzzz[i] * fe_0 + tr_x_xxyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 224-252 components of targeted buffer : HI

    auto tr_x_xxyzz_xxxxxx = pbuffer.data(idx_dip_hi + 224);

    auto tr_x_xxyzz_xxxxxy = pbuffer.data(idx_dip_hi + 225);

    auto tr_x_xxyzz_xxxxxz = pbuffer.data(idx_dip_hi + 226);

    auto tr_x_xxyzz_xxxxyy = pbuffer.data(idx_dip_hi + 227);

    auto tr_x_xxyzz_xxxxyz = pbuffer.data(idx_dip_hi + 228);

    auto tr_x_xxyzz_xxxxzz = pbuffer.data(idx_dip_hi + 229);

    auto tr_x_xxyzz_xxxyyy = pbuffer.data(idx_dip_hi + 230);

    auto tr_x_xxyzz_xxxyyz = pbuffer.data(idx_dip_hi + 231);

    auto tr_x_xxyzz_xxxyzz = pbuffer.data(idx_dip_hi + 232);

    auto tr_x_xxyzz_xxxzzz = pbuffer.data(idx_dip_hi + 233);

    auto tr_x_xxyzz_xxyyyy = pbuffer.data(idx_dip_hi + 234);

    auto tr_x_xxyzz_xxyyyz = pbuffer.data(idx_dip_hi + 235);

    auto tr_x_xxyzz_xxyyzz = pbuffer.data(idx_dip_hi + 236);

    auto tr_x_xxyzz_xxyzzz = pbuffer.data(idx_dip_hi + 237);

    auto tr_x_xxyzz_xxzzzz = pbuffer.data(idx_dip_hi + 238);

    auto tr_x_xxyzz_xyyyyy = pbuffer.data(idx_dip_hi + 239);

    auto tr_x_xxyzz_xyyyyz = pbuffer.data(idx_dip_hi + 240);

    auto tr_x_xxyzz_xyyyzz = pbuffer.data(idx_dip_hi + 241);

    auto tr_x_xxyzz_xyyzzz = pbuffer.data(idx_dip_hi + 242);

    auto tr_x_xxyzz_xyzzzz = pbuffer.data(idx_dip_hi + 243);

    auto tr_x_xxyzz_xzzzzz = pbuffer.data(idx_dip_hi + 244);

    auto tr_x_xxyzz_yyyyyy = pbuffer.data(idx_dip_hi + 245);

    auto tr_x_xxyzz_yyyyyz = pbuffer.data(idx_dip_hi + 246);

    auto tr_x_xxyzz_yyyyzz = pbuffer.data(idx_dip_hi + 247);

    auto tr_x_xxyzz_yyyzzz = pbuffer.data(idx_dip_hi + 248);

    auto tr_x_xxyzz_yyzzzz = pbuffer.data(idx_dip_hi + 249);

    auto tr_x_xxyzz_yzzzzz = pbuffer.data(idx_dip_hi + 250);

    auto tr_x_xxyzz_zzzzzz = pbuffer.data(idx_dip_hi + 251);

#pragma omp simd aligned(pa_y,                  \
                             tr_x_xxyzz_xxxxxx, \
                             tr_x_xxyzz_xxxxxy, \
                             tr_x_xxyzz_xxxxxz, \
                             tr_x_xxyzz_xxxxyy, \
                             tr_x_xxyzz_xxxxyz, \
                             tr_x_xxyzz_xxxxzz, \
                             tr_x_xxyzz_xxxyyy, \
                             tr_x_xxyzz_xxxyyz, \
                             tr_x_xxyzz_xxxyzz, \
                             tr_x_xxyzz_xxxzzz, \
                             tr_x_xxyzz_xxyyyy, \
                             tr_x_xxyzz_xxyyyz, \
                             tr_x_xxyzz_xxyyzz, \
                             tr_x_xxyzz_xxyzzz, \
                             tr_x_xxyzz_xxzzzz, \
                             tr_x_xxyzz_xyyyyy, \
                             tr_x_xxyzz_xyyyyz, \
                             tr_x_xxyzz_xyyyzz, \
                             tr_x_xxyzz_xyyzzz, \
                             tr_x_xxyzz_xyzzzz, \
                             tr_x_xxyzz_xzzzzz, \
                             tr_x_xxyzz_yyyyyy, \
                             tr_x_xxyzz_yyyyyz, \
                             tr_x_xxyzz_yyyyzz, \
                             tr_x_xxyzz_yyyzzz, \
                             tr_x_xxyzz_yyzzzz, \
                             tr_x_xxyzz_yzzzzz, \
                             tr_x_xxyzz_zzzzzz, \
                             tr_x_xxzz_xxxxx,   \
                             tr_x_xxzz_xxxxxx,  \
                             tr_x_xxzz_xxxxxy,  \
                             tr_x_xxzz_xxxxxz,  \
                             tr_x_xxzz_xxxxy,   \
                             tr_x_xxzz_xxxxyy,  \
                             tr_x_xxzz_xxxxyz,  \
                             tr_x_xxzz_xxxxz,   \
                             tr_x_xxzz_xxxxzz,  \
                             tr_x_xxzz_xxxyy,   \
                             tr_x_xxzz_xxxyyy,  \
                             tr_x_xxzz_xxxyyz,  \
                             tr_x_xxzz_xxxyz,   \
                             tr_x_xxzz_xxxyzz,  \
                             tr_x_xxzz_xxxzz,   \
                             tr_x_xxzz_xxxzzz,  \
                             tr_x_xxzz_xxyyy,   \
                             tr_x_xxzz_xxyyyy,  \
                             tr_x_xxzz_xxyyyz,  \
                             tr_x_xxzz_xxyyz,   \
                             tr_x_xxzz_xxyyzz,  \
                             tr_x_xxzz_xxyzz,   \
                             tr_x_xxzz_xxyzzz,  \
                             tr_x_xxzz_xxzzz,   \
                             tr_x_xxzz_xxzzzz,  \
                             tr_x_xxzz_xyyyy,   \
                             tr_x_xxzz_xyyyyy,  \
                             tr_x_xxzz_xyyyyz,  \
                             tr_x_xxzz_xyyyz,   \
                             tr_x_xxzz_xyyyzz,  \
                             tr_x_xxzz_xyyzz,   \
                             tr_x_xxzz_xyyzzz,  \
                             tr_x_xxzz_xyzzz,   \
                             tr_x_xxzz_xyzzzz,  \
                             tr_x_xxzz_xzzzz,   \
                             tr_x_xxzz_xzzzzz,  \
                             tr_x_xxzz_yyyyy,   \
                             tr_x_xxzz_yyyyyy,  \
                             tr_x_xxzz_yyyyyz,  \
                             tr_x_xxzz_yyyyz,   \
                             tr_x_xxzz_yyyyzz,  \
                             tr_x_xxzz_yyyzz,   \
                             tr_x_xxzz_yyyzzz,  \
                             tr_x_xxzz_yyzzz,   \
                             tr_x_xxzz_yyzzzz,  \
                             tr_x_xxzz_yzzzz,   \
                             tr_x_xxzz_yzzzzz,  \
                             tr_x_xxzz_zzzzz,   \
                             tr_x_xxzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzz_xxxxxx[i] = tr_x_xxzz_xxxxxx[i] * pa_y[i];

        tr_x_xxyzz_xxxxxy[i] = tr_x_xxzz_xxxxx[i] * fe_0 + tr_x_xxzz_xxxxxy[i] * pa_y[i];

        tr_x_xxyzz_xxxxxz[i] = tr_x_xxzz_xxxxxz[i] * pa_y[i];

        tr_x_xxyzz_xxxxyy[i] = 2.0 * tr_x_xxzz_xxxxy[i] * fe_0 + tr_x_xxzz_xxxxyy[i] * pa_y[i];

        tr_x_xxyzz_xxxxyz[i] = tr_x_xxzz_xxxxz[i] * fe_0 + tr_x_xxzz_xxxxyz[i] * pa_y[i];

        tr_x_xxyzz_xxxxzz[i] = tr_x_xxzz_xxxxzz[i] * pa_y[i];

        tr_x_xxyzz_xxxyyy[i] = 3.0 * tr_x_xxzz_xxxyy[i] * fe_0 + tr_x_xxzz_xxxyyy[i] * pa_y[i];

        tr_x_xxyzz_xxxyyz[i] = 2.0 * tr_x_xxzz_xxxyz[i] * fe_0 + tr_x_xxzz_xxxyyz[i] * pa_y[i];

        tr_x_xxyzz_xxxyzz[i] = tr_x_xxzz_xxxzz[i] * fe_0 + tr_x_xxzz_xxxyzz[i] * pa_y[i];

        tr_x_xxyzz_xxxzzz[i] = tr_x_xxzz_xxxzzz[i] * pa_y[i];

        tr_x_xxyzz_xxyyyy[i] = 4.0 * tr_x_xxzz_xxyyy[i] * fe_0 + tr_x_xxzz_xxyyyy[i] * pa_y[i];

        tr_x_xxyzz_xxyyyz[i] = 3.0 * tr_x_xxzz_xxyyz[i] * fe_0 + tr_x_xxzz_xxyyyz[i] * pa_y[i];

        tr_x_xxyzz_xxyyzz[i] = 2.0 * tr_x_xxzz_xxyzz[i] * fe_0 + tr_x_xxzz_xxyyzz[i] * pa_y[i];

        tr_x_xxyzz_xxyzzz[i] = tr_x_xxzz_xxzzz[i] * fe_0 + tr_x_xxzz_xxyzzz[i] * pa_y[i];

        tr_x_xxyzz_xxzzzz[i] = tr_x_xxzz_xxzzzz[i] * pa_y[i];

        tr_x_xxyzz_xyyyyy[i] = 5.0 * tr_x_xxzz_xyyyy[i] * fe_0 + tr_x_xxzz_xyyyyy[i] * pa_y[i];

        tr_x_xxyzz_xyyyyz[i] = 4.0 * tr_x_xxzz_xyyyz[i] * fe_0 + tr_x_xxzz_xyyyyz[i] * pa_y[i];

        tr_x_xxyzz_xyyyzz[i] = 3.0 * tr_x_xxzz_xyyzz[i] * fe_0 + tr_x_xxzz_xyyyzz[i] * pa_y[i];

        tr_x_xxyzz_xyyzzz[i] = 2.0 * tr_x_xxzz_xyzzz[i] * fe_0 + tr_x_xxzz_xyyzzz[i] * pa_y[i];

        tr_x_xxyzz_xyzzzz[i] = tr_x_xxzz_xzzzz[i] * fe_0 + tr_x_xxzz_xyzzzz[i] * pa_y[i];

        tr_x_xxyzz_xzzzzz[i] = tr_x_xxzz_xzzzzz[i] * pa_y[i];

        tr_x_xxyzz_yyyyyy[i] = 6.0 * tr_x_xxzz_yyyyy[i] * fe_0 + tr_x_xxzz_yyyyyy[i] * pa_y[i];

        tr_x_xxyzz_yyyyyz[i] = 5.0 * tr_x_xxzz_yyyyz[i] * fe_0 + tr_x_xxzz_yyyyyz[i] * pa_y[i];

        tr_x_xxyzz_yyyyzz[i] = 4.0 * tr_x_xxzz_yyyzz[i] * fe_0 + tr_x_xxzz_yyyyzz[i] * pa_y[i];

        tr_x_xxyzz_yyyzzz[i] = 3.0 * tr_x_xxzz_yyzzz[i] * fe_0 + tr_x_xxzz_yyyzzz[i] * pa_y[i];

        tr_x_xxyzz_yyzzzz[i] = 2.0 * tr_x_xxzz_yzzzz[i] * fe_0 + tr_x_xxzz_yyzzzz[i] * pa_y[i];

        tr_x_xxyzz_yzzzzz[i] = tr_x_xxzz_zzzzz[i] * fe_0 + tr_x_xxzz_yzzzzz[i] * pa_y[i];

        tr_x_xxyzz_zzzzzz[i] = tr_x_xxzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : HI

    auto tr_x_xxzzz_xxxxxx = pbuffer.data(idx_dip_hi + 252);

    auto tr_x_xxzzz_xxxxxy = pbuffer.data(idx_dip_hi + 253);

    auto tr_x_xxzzz_xxxxxz = pbuffer.data(idx_dip_hi + 254);

    auto tr_x_xxzzz_xxxxyy = pbuffer.data(idx_dip_hi + 255);

    auto tr_x_xxzzz_xxxxyz = pbuffer.data(idx_dip_hi + 256);

    auto tr_x_xxzzz_xxxxzz = pbuffer.data(idx_dip_hi + 257);

    auto tr_x_xxzzz_xxxyyy = pbuffer.data(idx_dip_hi + 258);

    auto tr_x_xxzzz_xxxyyz = pbuffer.data(idx_dip_hi + 259);

    auto tr_x_xxzzz_xxxyzz = pbuffer.data(idx_dip_hi + 260);

    auto tr_x_xxzzz_xxxzzz = pbuffer.data(idx_dip_hi + 261);

    auto tr_x_xxzzz_xxyyyy = pbuffer.data(idx_dip_hi + 262);

    auto tr_x_xxzzz_xxyyyz = pbuffer.data(idx_dip_hi + 263);

    auto tr_x_xxzzz_xxyyzz = pbuffer.data(idx_dip_hi + 264);

    auto tr_x_xxzzz_xxyzzz = pbuffer.data(idx_dip_hi + 265);

    auto tr_x_xxzzz_xxzzzz = pbuffer.data(idx_dip_hi + 266);

    auto tr_x_xxzzz_xyyyyy = pbuffer.data(idx_dip_hi + 267);

    auto tr_x_xxzzz_xyyyyz = pbuffer.data(idx_dip_hi + 268);

    auto tr_x_xxzzz_xyyyzz = pbuffer.data(idx_dip_hi + 269);

    auto tr_x_xxzzz_xyyzzz = pbuffer.data(idx_dip_hi + 270);

    auto tr_x_xxzzz_xyzzzz = pbuffer.data(idx_dip_hi + 271);

    auto tr_x_xxzzz_xzzzzz = pbuffer.data(idx_dip_hi + 272);

    auto tr_x_xxzzz_yyyyyy = pbuffer.data(idx_dip_hi + 273);

    auto tr_x_xxzzz_yyyyyz = pbuffer.data(idx_dip_hi + 274);

    auto tr_x_xxzzz_yyyyzz = pbuffer.data(idx_dip_hi + 275);

    auto tr_x_xxzzz_yyyzzz = pbuffer.data(idx_dip_hi + 276);

    auto tr_x_xxzzz_yyzzzz = pbuffer.data(idx_dip_hi + 277);

    auto tr_x_xxzzz_yzzzzz = pbuffer.data(idx_dip_hi + 278);

    auto tr_x_xxzzz_zzzzzz = pbuffer.data(idx_dip_hi + 279);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_x_xxz_xxxxxx,   \
                             tr_x_xxz_xxxxxy,   \
                             tr_x_xxz_xxxxxz,   \
                             tr_x_xxz_xxxxyy,   \
                             tr_x_xxz_xxxxyz,   \
                             tr_x_xxz_xxxxzz,   \
                             tr_x_xxz_xxxyyy,   \
                             tr_x_xxz_xxxyyz,   \
                             tr_x_xxz_xxxyzz,   \
                             tr_x_xxz_xxxzzz,   \
                             tr_x_xxz_xxyyyy,   \
                             tr_x_xxz_xxyyyz,   \
                             tr_x_xxz_xxyyzz,   \
                             tr_x_xxz_xxyzzz,   \
                             tr_x_xxz_xxzzzz,   \
                             tr_x_xxz_xyyyyy,   \
                             tr_x_xxz_xyyyyz,   \
                             tr_x_xxz_xyyyzz,   \
                             tr_x_xxz_xyyzzz,   \
                             tr_x_xxz_xyzzzz,   \
                             tr_x_xxz_xzzzzz,   \
                             tr_x_xxz_yyyyyy,   \
                             tr_x_xxzz_xxxxx,   \
                             tr_x_xxzz_xxxxxx,  \
                             tr_x_xxzz_xxxxxy,  \
                             tr_x_xxzz_xxxxxz,  \
                             tr_x_xxzz_xxxxy,   \
                             tr_x_xxzz_xxxxyy,  \
                             tr_x_xxzz_xxxxyz,  \
                             tr_x_xxzz_xxxxz,   \
                             tr_x_xxzz_xxxxzz,  \
                             tr_x_xxzz_xxxyy,   \
                             tr_x_xxzz_xxxyyy,  \
                             tr_x_xxzz_xxxyyz,  \
                             tr_x_xxzz_xxxyz,   \
                             tr_x_xxzz_xxxyzz,  \
                             tr_x_xxzz_xxxzz,   \
                             tr_x_xxzz_xxxzzz,  \
                             tr_x_xxzz_xxyyy,   \
                             tr_x_xxzz_xxyyyy,  \
                             tr_x_xxzz_xxyyyz,  \
                             tr_x_xxzz_xxyyz,   \
                             tr_x_xxzz_xxyyzz,  \
                             tr_x_xxzz_xxyzz,   \
                             tr_x_xxzz_xxyzzz,  \
                             tr_x_xxzz_xxzzz,   \
                             tr_x_xxzz_xxzzzz,  \
                             tr_x_xxzz_xyyyy,   \
                             tr_x_xxzz_xyyyyy,  \
                             tr_x_xxzz_xyyyyz,  \
                             tr_x_xxzz_xyyyz,   \
                             tr_x_xxzz_xyyyzz,  \
                             tr_x_xxzz_xyyzz,   \
                             tr_x_xxzz_xyyzzz,  \
                             tr_x_xxzz_xyzzz,   \
                             tr_x_xxzz_xyzzzz,  \
                             tr_x_xxzz_xzzzz,   \
                             tr_x_xxzz_xzzzzz,  \
                             tr_x_xxzz_yyyyyy,  \
                             tr_x_xxzzz_xxxxxx, \
                             tr_x_xxzzz_xxxxxy, \
                             tr_x_xxzzz_xxxxxz, \
                             tr_x_xxzzz_xxxxyy, \
                             tr_x_xxzzz_xxxxyz, \
                             tr_x_xxzzz_xxxxzz, \
                             tr_x_xxzzz_xxxyyy, \
                             tr_x_xxzzz_xxxyyz, \
                             tr_x_xxzzz_xxxyzz, \
                             tr_x_xxzzz_xxxzzz, \
                             tr_x_xxzzz_xxyyyy, \
                             tr_x_xxzzz_xxyyyz, \
                             tr_x_xxzzz_xxyyzz, \
                             tr_x_xxzzz_xxyzzz, \
                             tr_x_xxzzz_xxzzzz, \
                             tr_x_xxzzz_xyyyyy, \
                             tr_x_xxzzz_xyyyyz, \
                             tr_x_xxzzz_xyyyzz, \
                             tr_x_xxzzz_xyyzzz, \
                             tr_x_xxzzz_xyzzzz, \
                             tr_x_xxzzz_xzzzzz, \
                             tr_x_xxzzz_yyyyyy, \
                             tr_x_xxzzz_yyyyyz, \
                             tr_x_xxzzz_yyyyzz, \
                             tr_x_xxzzz_yyyzzz, \
                             tr_x_xxzzz_yyzzzz, \
                             tr_x_xxzzz_yzzzzz, \
                             tr_x_xxzzz_zzzzzz, \
                             tr_x_xzzz_yyyyyz,  \
                             tr_x_xzzz_yyyyzz,  \
                             tr_x_xzzz_yyyzzz,  \
                             tr_x_xzzz_yyzzzz,  \
                             tr_x_xzzz_yzzzzz,  \
                             tr_x_xzzz_zzzzzz,  \
                             tr_x_zzz_yyyyyz,   \
                             tr_x_zzz_yyyyzz,   \
                             tr_x_zzz_yyyzzz,   \
                             tr_x_zzz_yyzzzz,   \
                             tr_x_zzz_yzzzzz,   \
                             tr_x_zzz_zzzzzz,   \
                             ts_xzzz_yyyyyz,    \
                             ts_xzzz_yyyyzz,    \
                             ts_xzzz_yyyzzz,    \
                             ts_xzzz_yyzzzz,    \
                             ts_xzzz_yzzzzz,    \
                             ts_xzzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzz_xxxxxx[i] = 2.0 * tr_x_xxz_xxxxxx[i] * fe_0 + tr_x_xxzz_xxxxxx[i] * pa_z[i];

        tr_x_xxzzz_xxxxxy[i] = 2.0 * tr_x_xxz_xxxxxy[i] * fe_0 + tr_x_xxzz_xxxxxy[i] * pa_z[i];

        tr_x_xxzzz_xxxxxz[i] = 2.0 * tr_x_xxz_xxxxxz[i] * fe_0 + tr_x_xxzz_xxxxx[i] * fe_0 + tr_x_xxzz_xxxxxz[i] * pa_z[i];

        tr_x_xxzzz_xxxxyy[i] = 2.0 * tr_x_xxz_xxxxyy[i] * fe_0 + tr_x_xxzz_xxxxyy[i] * pa_z[i];

        tr_x_xxzzz_xxxxyz[i] = 2.0 * tr_x_xxz_xxxxyz[i] * fe_0 + tr_x_xxzz_xxxxy[i] * fe_0 + tr_x_xxzz_xxxxyz[i] * pa_z[i];

        tr_x_xxzzz_xxxxzz[i] = 2.0 * tr_x_xxz_xxxxzz[i] * fe_0 + 2.0 * tr_x_xxzz_xxxxz[i] * fe_0 + tr_x_xxzz_xxxxzz[i] * pa_z[i];

        tr_x_xxzzz_xxxyyy[i] = 2.0 * tr_x_xxz_xxxyyy[i] * fe_0 + tr_x_xxzz_xxxyyy[i] * pa_z[i];

        tr_x_xxzzz_xxxyyz[i] = 2.0 * tr_x_xxz_xxxyyz[i] * fe_0 + tr_x_xxzz_xxxyy[i] * fe_0 + tr_x_xxzz_xxxyyz[i] * pa_z[i];

        tr_x_xxzzz_xxxyzz[i] = 2.0 * tr_x_xxz_xxxyzz[i] * fe_0 + 2.0 * tr_x_xxzz_xxxyz[i] * fe_0 + tr_x_xxzz_xxxyzz[i] * pa_z[i];

        tr_x_xxzzz_xxxzzz[i] = 2.0 * tr_x_xxz_xxxzzz[i] * fe_0 + 3.0 * tr_x_xxzz_xxxzz[i] * fe_0 + tr_x_xxzz_xxxzzz[i] * pa_z[i];

        tr_x_xxzzz_xxyyyy[i] = 2.0 * tr_x_xxz_xxyyyy[i] * fe_0 + tr_x_xxzz_xxyyyy[i] * pa_z[i];

        tr_x_xxzzz_xxyyyz[i] = 2.0 * tr_x_xxz_xxyyyz[i] * fe_0 + tr_x_xxzz_xxyyy[i] * fe_0 + tr_x_xxzz_xxyyyz[i] * pa_z[i];

        tr_x_xxzzz_xxyyzz[i] = 2.0 * tr_x_xxz_xxyyzz[i] * fe_0 + 2.0 * tr_x_xxzz_xxyyz[i] * fe_0 + tr_x_xxzz_xxyyzz[i] * pa_z[i];

        tr_x_xxzzz_xxyzzz[i] = 2.0 * tr_x_xxz_xxyzzz[i] * fe_0 + 3.0 * tr_x_xxzz_xxyzz[i] * fe_0 + tr_x_xxzz_xxyzzz[i] * pa_z[i];

        tr_x_xxzzz_xxzzzz[i] = 2.0 * tr_x_xxz_xxzzzz[i] * fe_0 + 4.0 * tr_x_xxzz_xxzzz[i] * fe_0 + tr_x_xxzz_xxzzzz[i] * pa_z[i];

        tr_x_xxzzz_xyyyyy[i] = 2.0 * tr_x_xxz_xyyyyy[i] * fe_0 + tr_x_xxzz_xyyyyy[i] * pa_z[i];

        tr_x_xxzzz_xyyyyz[i] = 2.0 * tr_x_xxz_xyyyyz[i] * fe_0 + tr_x_xxzz_xyyyy[i] * fe_0 + tr_x_xxzz_xyyyyz[i] * pa_z[i];

        tr_x_xxzzz_xyyyzz[i] = 2.0 * tr_x_xxz_xyyyzz[i] * fe_0 + 2.0 * tr_x_xxzz_xyyyz[i] * fe_0 + tr_x_xxzz_xyyyzz[i] * pa_z[i];

        tr_x_xxzzz_xyyzzz[i] = 2.0 * tr_x_xxz_xyyzzz[i] * fe_0 + 3.0 * tr_x_xxzz_xyyzz[i] * fe_0 + tr_x_xxzz_xyyzzz[i] * pa_z[i];

        tr_x_xxzzz_xyzzzz[i] = 2.0 * tr_x_xxz_xyzzzz[i] * fe_0 + 4.0 * tr_x_xxzz_xyzzz[i] * fe_0 + tr_x_xxzz_xyzzzz[i] * pa_z[i];

        tr_x_xxzzz_xzzzzz[i] = 2.0 * tr_x_xxz_xzzzzz[i] * fe_0 + 5.0 * tr_x_xxzz_xzzzz[i] * fe_0 + tr_x_xxzz_xzzzzz[i] * pa_z[i];

        tr_x_xxzzz_yyyyyy[i] = 2.0 * tr_x_xxz_yyyyyy[i] * fe_0 + tr_x_xxzz_yyyyyy[i] * pa_z[i];

        tr_x_xxzzz_yyyyyz[i] = tr_x_zzz_yyyyyz[i] * fe_0 + ts_xzzz_yyyyyz[i] * fe_0 + tr_x_xzzz_yyyyyz[i] * pa_x[i];

        tr_x_xxzzz_yyyyzz[i] = tr_x_zzz_yyyyzz[i] * fe_0 + ts_xzzz_yyyyzz[i] * fe_0 + tr_x_xzzz_yyyyzz[i] * pa_x[i];

        tr_x_xxzzz_yyyzzz[i] = tr_x_zzz_yyyzzz[i] * fe_0 + ts_xzzz_yyyzzz[i] * fe_0 + tr_x_xzzz_yyyzzz[i] * pa_x[i];

        tr_x_xxzzz_yyzzzz[i] = tr_x_zzz_yyzzzz[i] * fe_0 + ts_xzzz_yyzzzz[i] * fe_0 + tr_x_xzzz_yyzzzz[i] * pa_x[i];

        tr_x_xxzzz_yzzzzz[i] = tr_x_zzz_yzzzzz[i] * fe_0 + ts_xzzz_yzzzzz[i] * fe_0 + tr_x_xzzz_yzzzzz[i] * pa_x[i];

        tr_x_xxzzz_zzzzzz[i] = tr_x_zzz_zzzzzz[i] * fe_0 + ts_xzzz_zzzzzz[i] * fe_0 + tr_x_xzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : HI

    auto tr_x_xyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 280);

    auto tr_x_xyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 281);

    auto tr_x_xyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 282);

    auto tr_x_xyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 283);

    auto tr_x_xyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 284);

    auto tr_x_xyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 285);

    auto tr_x_xyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 286);

    auto tr_x_xyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 287);

    auto tr_x_xyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 288);

    auto tr_x_xyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 289);

    auto tr_x_xyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 290);

    auto tr_x_xyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 291);

    auto tr_x_xyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 292);

    auto tr_x_xyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 293);

    auto tr_x_xyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 294);

    auto tr_x_xyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 295);

    auto tr_x_xyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 296);

    auto tr_x_xyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 297);

    auto tr_x_xyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 298);

    auto tr_x_xyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 299);

    auto tr_x_xyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 300);

    auto tr_x_xyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 301);

    auto tr_x_xyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 302);

    auto tr_x_xyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 303);

    auto tr_x_xyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 304);

    auto tr_x_xyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 305);

    auto tr_x_xyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 306);

    auto tr_x_xyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 307);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_x_xyy_xxxxxx,   \
                             tr_x_xyy_xxxxxz,   \
                             tr_x_xyy_xxxxzz,   \
                             tr_x_xyy_xxxzzz,   \
                             tr_x_xyy_xxzzzz,   \
                             tr_x_xyy_xzzzzz,   \
                             tr_x_xyyy_xxxxxx,  \
                             tr_x_xyyy_xxxxxz,  \
                             tr_x_xyyy_xxxxzz,  \
                             tr_x_xyyy_xxxzzz,  \
                             tr_x_xyyy_xxzzzz,  \
                             tr_x_xyyy_xzzzzz,  \
                             tr_x_xyyyy_xxxxxx, \
                             tr_x_xyyyy_xxxxxy, \
                             tr_x_xyyyy_xxxxxz, \
                             tr_x_xyyyy_xxxxyy, \
                             tr_x_xyyyy_xxxxyz, \
                             tr_x_xyyyy_xxxxzz, \
                             tr_x_xyyyy_xxxyyy, \
                             tr_x_xyyyy_xxxyyz, \
                             tr_x_xyyyy_xxxyzz, \
                             tr_x_xyyyy_xxxzzz, \
                             tr_x_xyyyy_xxyyyy, \
                             tr_x_xyyyy_xxyyyz, \
                             tr_x_xyyyy_xxyyzz, \
                             tr_x_xyyyy_xxyzzz, \
                             tr_x_xyyyy_xxzzzz, \
                             tr_x_xyyyy_xyyyyy, \
                             tr_x_xyyyy_xyyyyz, \
                             tr_x_xyyyy_xyyyzz, \
                             tr_x_xyyyy_xyyzzz, \
                             tr_x_xyyyy_xyzzzz, \
                             tr_x_xyyyy_xzzzzz, \
                             tr_x_xyyyy_yyyyyy, \
                             tr_x_xyyyy_yyyyyz, \
                             tr_x_xyyyy_yyyyzz, \
                             tr_x_xyyyy_yyyzzz, \
                             tr_x_xyyyy_yyzzzz, \
                             tr_x_xyyyy_yzzzzz, \
                             tr_x_xyyyy_zzzzzz, \
                             tr_x_yyyy_xxxxxy,  \
                             tr_x_yyyy_xxxxy,   \
                             tr_x_yyyy_xxxxyy,  \
                             tr_x_yyyy_xxxxyz,  \
                             tr_x_yyyy_xxxyy,   \
                             tr_x_yyyy_xxxyyy,  \
                             tr_x_yyyy_xxxyyz,  \
                             tr_x_yyyy_xxxyz,   \
                             tr_x_yyyy_xxxyzz,  \
                             tr_x_yyyy_xxyyy,   \
                             tr_x_yyyy_xxyyyy,  \
                             tr_x_yyyy_xxyyyz,  \
                             tr_x_yyyy_xxyyz,   \
                             tr_x_yyyy_xxyyzz,  \
                             tr_x_yyyy_xxyzz,   \
                             tr_x_yyyy_xxyzzz,  \
                             tr_x_yyyy_xyyyy,   \
                             tr_x_yyyy_xyyyyy,  \
                             tr_x_yyyy_xyyyyz,  \
                             tr_x_yyyy_xyyyz,   \
                             tr_x_yyyy_xyyyzz,  \
                             tr_x_yyyy_xyyzz,   \
                             tr_x_yyyy_xyyzzz,  \
                             tr_x_yyyy_xyzzz,   \
                             tr_x_yyyy_xyzzzz,  \
                             tr_x_yyyy_yyyyy,   \
                             tr_x_yyyy_yyyyyy,  \
                             tr_x_yyyy_yyyyyz,  \
                             tr_x_yyyy_yyyyz,   \
                             tr_x_yyyy_yyyyzz,  \
                             tr_x_yyyy_yyyzz,   \
                             tr_x_yyyy_yyyzzz,  \
                             tr_x_yyyy_yyzzz,   \
                             tr_x_yyyy_yyzzzz,  \
                             tr_x_yyyy_yzzzz,   \
                             tr_x_yyyy_yzzzzz,  \
                             tr_x_yyyy_zzzzzz,  \
                             ts_yyyy_xxxxxy,    \
                             ts_yyyy_xxxxyy,    \
                             ts_yyyy_xxxxyz,    \
                             ts_yyyy_xxxyyy,    \
                             ts_yyyy_xxxyyz,    \
                             ts_yyyy_xxxyzz,    \
                             ts_yyyy_xxyyyy,    \
                             ts_yyyy_xxyyyz,    \
                             ts_yyyy_xxyyzz,    \
                             ts_yyyy_xxyzzz,    \
                             ts_yyyy_xyyyyy,    \
                             ts_yyyy_xyyyyz,    \
                             ts_yyyy_xyyyzz,    \
                             ts_yyyy_xyyzzz,    \
                             ts_yyyy_xyzzzz,    \
                             ts_yyyy_yyyyyy,    \
                             ts_yyyy_yyyyyz,    \
                             ts_yyyy_yyyyzz,    \
                             ts_yyyy_yyyzzz,    \
                             ts_yyyy_yyzzzz,    \
                             ts_yyyy_yzzzzz,    \
                             ts_yyyy_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyy_xxxxxx[i] = 3.0 * tr_x_xyy_xxxxxx[i] * fe_0 + tr_x_xyyy_xxxxxx[i] * pa_y[i];

        tr_x_xyyyy_xxxxxy[i] = 5.0 * tr_x_yyyy_xxxxy[i] * fe_0 + ts_yyyy_xxxxxy[i] * fe_0 + tr_x_yyyy_xxxxxy[i] * pa_x[i];

        tr_x_xyyyy_xxxxxz[i] = 3.0 * tr_x_xyy_xxxxxz[i] * fe_0 + tr_x_xyyy_xxxxxz[i] * pa_y[i];

        tr_x_xyyyy_xxxxyy[i] = 4.0 * tr_x_yyyy_xxxyy[i] * fe_0 + ts_yyyy_xxxxyy[i] * fe_0 + tr_x_yyyy_xxxxyy[i] * pa_x[i];

        tr_x_xyyyy_xxxxyz[i] = 4.0 * tr_x_yyyy_xxxyz[i] * fe_0 + ts_yyyy_xxxxyz[i] * fe_0 + tr_x_yyyy_xxxxyz[i] * pa_x[i];

        tr_x_xyyyy_xxxxzz[i] = 3.0 * tr_x_xyy_xxxxzz[i] * fe_0 + tr_x_xyyy_xxxxzz[i] * pa_y[i];

        tr_x_xyyyy_xxxyyy[i] = 3.0 * tr_x_yyyy_xxyyy[i] * fe_0 + ts_yyyy_xxxyyy[i] * fe_0 + tr_x_yyyy_xxxyyy[i] * pa_x[i];

        tr_x_xyyyy_xxxyyz[i] = 3.0 * tr_x_yyyy_xxyyz[i] * fe_0 + ts_yyyy_xxxyyz[i] * fe_0 + tr_x_yyyy_xxxyyz[i] * pa_x[i];

        tr_x_xyyyy_xxxyzz[i] = 3.0 * tr_x_yyyy_xxyzz[i] * fe_0 + ts_yyyy_xxxyzz[i] * fe_0 + tr_x_yyyy_xxxyzz[i] * pa_x[i];

        tr_x_xyyyy_xxxzzz[i] = 3.0 * tr_x_xyy_xxxzzz[i] * fe_0 + tr_x_xyyy_xxxzzz[i] * pa_y[i];

        tr_x_xyyyy_xxyyyy[i] = 2.0 * tr_x_yyyy_xyyyy[i] * fe_0 + ts_yyyy_xxyyyy[i] * fe_0 + tr_x_yyyy_xxyyyy[i] * pa_x[i];

        tr_x_xyyyy_xxyyyz[i] = 2.0 * tr_x_yyyy_xyyyz[i] * fe_0 + ts_yyyy_xxyyyz[i] * fe_0 + tr_x_yyyy_xxyyyz[i] * pa_x[i];

        tr_x_xyyyy_xxyyzz[i] = 2.0 * tr_x_yyyy_xyyzz[i] * fe_0 + ts_yyyy_xxyyzz[i] * fe_0 + tr_x_yyyy_xxyyzz[i] * pa_x[i];

        tr_x_xyyyy_xxyzzz[i] = 2.0 * tr_x_yyyy_xyzzz[i] * fe_0 + ts_yyyy_xxyzzz[i] * fe_0 + tr_x_yyyy_xxyzzz[i] * pa_x[i];

        tr_x_xyyyy_xxzzzz[i] = 3.0 * tr_x_xyy_xxzzzz[i] * fe_0 + tr_x_xyyy_xxzzzz[i] * pa_y[i];

        tr_x_xyyyy_xyyyyy[i] = tr_x_yyyy_yyyyy[i] * fe_0 + ts_yyyy_xyyyyy[i] * fe_0 + tr_x_yyyy_xyyyyy[i] * pa_x[i];

        tr_x_xyyyy_xyyyyz[i] = tr_x_yyyy_yyyyz[i] * fe_0 + ts_yyyy_xyyyyz[i] * fe_0 + tr_x_yyyy_xyyyyz[i] * pa_x[i];

        tr_x_xyyyy_xyyyzz[i] = tr_x_yyyy_yyyzz[i] * fe_0 + ts_yyyy_xyyyzz[i] * fe_0 + tr_x_yyyy_xyyyzz[i] * pa_x[i];

        tr_x_xyyyy_xyyzzz[i] = tr_x_yyyy_yyzzz[i] * fe_0 + ts_yyyy_xyyzzz[i] * fe_0 + tr_x_yyyy_xyyzzz[i] * pa_x[i];

        tr_x_xyyyy_xyzzzz[i] = tr_x_yyyy_yzzzz[i] * fe_0 + ts_yyyy_xyzzzz[i] * fe_0 + tr_x_yyyy_xyzzzz[i] * pa_x[i];

        tr_x_xyyyy_xzzzzz[i] = 3.0 * tr_x_xyy_xzzzzz[i] * fe_0 + tr_x_xyyy_xzzzzz[i] * pa_y[i];

        tr_x_xyyyy_yyyyyy[i] = ts_yyyy_yyyyyy[i] * fe_0 + tr_x_yyyy_yyyyyy[i] * pa_x[i];

        tr_x_xyyyy_yyyyyz[i] = ts_yyyy_yyyyyz[i] * fe_0 + tr_x_yyyy_yyyyyz[i] * pa_x[i];

        tr_x_xyyyy_yyyyzz[i] = ts_yyyy_yyyyzz[i] * fe_0 + tr_x_yyyy_yyyyzz[i] * pa_x[i];

        tr_x_xyyyy_yyyzzz[i] = ts_yyyy_yyyzzz[i] * fe_0 + tr_x_yyyy_yyyzzz[i] * pa_x[i];

        tr_x_xyyyy_yyzzzz[i] = ts_yyyy_yyzzzz[i] * fe_0 + tr_x_yyyy_yyzzzz[i] * pa_x[i];

        tr_x_xyyyy_yzzzzz[i] = ts_yyyy_yzzzzz[i] * fe_0 + tr_x_yyyy_yzzzzz[i] * pa_x[i];

        tr_x_xyyyy_zzzzzz[i] = ts_yyyy_zzzzzz[i] * fe_0 + tr_x_yyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 308-336 components of targeted buffer : HI

    auto tr_x_xyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 308);

    auto tr_x_xyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 309);

    auto tr_x_xyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 310);

    auto tr_x_xyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 311);

    auto tr_x_xyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 312);

    auto tr_x_xyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 313);

    auto tr_x_xyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 314);

    auto tr_x_xyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 315);

    auto tr_x_xyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 316);

    auto tr_x_xyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 317);

    auto tr_x_xyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 318);

    auto tr_x_xyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 319);

    auto tr_x_xyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 320);

    auto tr_x_xyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 321);

    auto tr_x_xyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 322);

    auto tr_x_xyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 323);

    auto tr_x_xyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 324);

    auto tr_x_xyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 325);

    auto tr_x_xyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 326);

    auto tr_x_xyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 327);

    auto tr_x_xyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 328);

    auto tr_x_xyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 329);

    auto tr_x_xyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 330);

    auto tr_x_xyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 331);

    auto tr_x_xyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 332);

    auto tr_x_xyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 333);

    auto tr_x_xyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 334);

    auto tr_x_xyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 335);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             tr_x_xyyy_xxxxxx,  \
                             tr_x_xyyy_xxxxxy,  \
                             tr_x_xyyy_xxxxy,   \
                             tr_x_xyyy_xxxxyy,  \
                             tr_x_xyyy_xxxxyz,  \
                             tr_x_xyyy_xxxyy,   \
                             tr_x_xyyy_xxxyyy,  \
                             tr_x_xyyy_xxxyyz,  \
                             tr_x_xyyy_xxxyz,   \
                             tr_x_xyyy_xxxyzz,  \
                             tr_x_xyyy_xxyyy,   \
                             tr_x_xyyy_xxyyyy,  \
                             tr_x_xyyy_xxyyyz,  \
                             tr_x_xyyy_xxyyz,   \
                             tr_x_xyyy_xxyyzz,  \
                             tr_x_xyyy_xxyzz,   \
                             tr_x_xyyy_xxyzzz,  \
                             tr_x_xyyy_xyyyy,   \
                             tr_x_xyyy_xyyyyy,  \
                             tr_x_xyyy_xyyyyz,  \
                             tr_x_xyyy_xyyyz,   \
                             tr_x_xyyy_xyyyzz,  \
                             tr_x_xyyy_xyyzz,   \
                             tr_x_xyyy_xyyzzz,  \
                             tr_x_xyyy_xyzzz,   \
                             tr_x_xyyy_xyzzzz,  \
                             tr_x_xyyy_yyyyyy,  \
                             tr_x_xyyyz_xxxxxx, \
                             tr_x_xyyyz_xxxxxy, \
                             tr_x_xyyyz_xxxxxz, \
                             tr_x_xyyyz_xxxxyy, \
                             tr_x_xyyyz_xxxxyz, \
                             tr_x_xyyyz_xxxxzz, \
                             tr_x_xyyyz_xxxyyy, \
                             tr_x_xyyyz_xxxyyz, \
                             tr_x_xyyyz_xxxyzz, \
                             tr_x_xyyyz_xxxzzz, \
                             tr_x_xyyyz_xxyyyy, \
                             tr_x_xyyyz_xxyyyz, \
                             tr_x_xyyyz_xxyyzz, \
                             tr_x_xyyyz_xxyzzz, \
                             tr_x_xyyyz_xxzzzz, \
                             tr_x_xyyyz_xyyyyy, \
                             tr_x_xyyyz_xyyyyz, \
                             tr_x_xyyyz_xyyyzz, \
                             tr_x_xyyyz_xyyzzz, \
                             tr_x_xyyyz_xyzzzz, \
                             tr_x_xyyyz_xzzzzz, \
                             tr_x_xyyyz_yyyyyy, \
                             tr_x_xyyyz_yyyyyz, \
                             tr_x_xyyyz_yyyyzz, \
                             tr_x_xyyyz_yyyzzz, \
                             tr_x_xyyyz_yyzzzz, \
                             tr_x_xyyyz_yzzzzz, \
                             tr_x_xyyyz_zzzzzz, \
                             tr_x_xyyz_xxxxxz,  \
                             tr_x_xyyz_xxxxzz,  \
                             tr_x_xyyz_xxxzzz,  \
                             tr_x_xyyz_xxzzzz,  \
                             tr_x_xyyz_xzzzzz,  \
                             tr_x_xyz_xxxxxz,   \
                             tr_x_xyz_xxxxzz,   \
                             tr_x_xyz_xxxzzz,   \
                             tr_x_xyz_xxzzzz,   \
                             tr_x_xyz_xzzzzz,   \
                             tr_x_yyyz_yyyyyz,  \
                             tr_x_yyyz_yyyyzz,  \
                             tr_x_yyyz_yyyzzz,  \
                             tr_x_yyyz_yyzzzz,  \
                             tr_x_yyyz_yzzzzz,  \
                             tr_x_yyyz_zzzzzz,  \
                             ts_yyyz_yyyyyz,    \
                             ts_yyyz_yyyyzz,    \
                             ts_yyyz_yyyzzz,    \
                             ts_yyyz_yyzzzz,    \
                             ts_yyyz_yzzzzz,    \
                             ts_yyyz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyz_xxxxxx[i] = tr_x_xyyy_xxxxxx[i] * pa_z[i];

        tr_x_xyyyz_xxxxxy[i] = tr_x_xyyy_xxxxxy[i] * pa_z[i];

        tr_x_xyyyz_xxxxxz[i] = 2.0 * tr_x_xyz_xxxxxz[i] * fe_0 + tr_x_xyyz_xxxxxz[i] * pa_y[i];

        tr_x_xyyyz_xxxxyy[i] = tr_x_xyyy_xxxxyy[i] * pa_z[i];

        tr_x_xyyyz_xxxxyz[i] = tr_x_xyyy_xxxxy[i] * fe_0 + tr_x_xyyy_xxxxyz[i] * pa_z[i];

        tr_x_xyyyz_xxxxzz[i] = 2.0 * tr_x_xyz_xxxxzz[i] * fe_0 + tr_x_xyyz_xxxxzz[i] * pa_y[i];

        tr_x_xyyyz_xxxyyy[i] = tr_x_xyyy_xxxyyy[i] * pa_z[i];

        tr_x_xyyyz_xxxyyz[i] = tr_x_xyyy_xxxyy[i] * fe_0 + tr_x_xyyy_xxxyyz[i] * pa_z[i];

        tr_x_xyyyz_xxxyzz[i] = 2.0 * tr_x_xyyy_xxxyz[i] * fe_0 + tr_x_xyyy_xxxyzz[i] * pa_z[i];

        tr_x_xyyyz_xxxzzz[i] = 2.0 * tr_x_xyz_xxxzzz[i] * fe_0 + tr_x_xyyz_xxxzzz[i] * pa_y[i];

        tr_x_xyyyz_xxyyyy[i] = tr_x_xyyy_xxyyyy[i] * pa_z[i];

        tr_x_xyyyz_xxyyyz[i] = tr_x_xyyy_xxyyy[i] * fe_0 + tr_x_xyyy_xxyyyz[i] * pa_z[i];

        tr_x_xyyyz_xxyyzz[i] = 2.0 * tr_x_xyyy_xxyyz[i] * fe_0 + tr_x_xyyy_xxyyzz[i] * pa_z[i];

        tr_x_xyyyz_xxyzzz[i] = 3.0 * tr_x_xyyy_xxyzz[i] * fe_0 + tr_x_xyyy_xxyzzz[i] * pa_z[i];

        tr_x_xyyyz_xxzzzz[i] = 2.0 * tr_x_xyz_xxzzzz[i] * fe_0 + tr_x_xyyz_xxzzzz[i] * pa_y[i];

        tr_x_xyyyz_xyyyyy[i] = tr_x_xyyy_xyyyyy[i] * pa_z[i];

        tr_x_xyyyz_xyyyyz[i] = tr_x_xyyy_xyyyy[i] * fe_0 + tr_x_xyyy_xyyyyz[i] * pa_z[i];

        tr_x_xyyyz_xyyyzz[i] = 2.0 * tr_x_xyyy_xyyyz[i] * fe_0 + tr_x_xyyy_xyyyzz[i] * pa_z[i];

        tr_x_xyyyz_xyyzzz[i] = 3.0 * tr_x_xyyy_xyyzz[i] * fe_0 + tr_x_xyyy_xyyzzz[i] * pa_z[i];

        tr_x_xyyyz_xyzzzz[i] = 4.0 * tr_x_xyyy_xyzzz[i] * fe_0 + tr_x_xyyy_xyzzzz[i] * pa_z[i];

        tr_x_xyyyz_xzzzzz[i] = 2.0 * tr_x_xyz_xzzzzz[i] * fe_0 + tr_x_xyyz_xzzzzz[i] * pa_y[i];

        tr_x_xyyyz_yyyyyy[i] = tr_x_xyyy_yyyyyy[i] * pa_z[i];

        tr_x_xyyyz_yyyyyz[i] = ts_yyyz_yyyyyz[i] * fe_0 + tr_x_yyyz_yyyyyz[i] * pa_x[i];

        tr_x_xyyyz_yyyyzz[i] = ts_yyyz_yyyyzz[i] * fe_0 + tr_x_yyyz_yyyyzz[i] * pa_x[i];

        tr_x_xyyyz_yyyzzz[i] = ts_yyyz_yyyzzz[i] * fe_0 + tr_x_yyyz_yyyzzz[i] * pa_x[i];

        tr_x_xyyyz_yyzzzz[i] = ts_yyyz_yyzzzz[i] * fe_0 + tr_x_yyyz_yyzzzz[i] * pa_x[i];

        tr_x_xyyyz_yzzzzz[i] = ts_yyyz_yzzzzz[i] * fe_0 + tr_x_yyyz_yzzzzz[i] * pa_x[i];

        tr_x_xyyyz_zzzzzz[i] = ts_yyyz_zzzzzz[i] * fe_0 + tr_x_yyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 336-364 components of targeted buffer : HI

    auto tr_x_xyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 336);

    auto tr_x_xyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 337);

    auto tr_x_xyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 338);

    auto tr_x_xyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 339);

    auto tr_x_xyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 340);

    auto tr_x_xyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 341);

    auto tr_x_xyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 342);

    auto tr_x_xyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 343);

    auto tr_x_xyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 344);

    auto tr_x_xyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 345);

    auto tr_x_xyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 346);

    auto tr_x_xyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 347);

    auto tr_x_xyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 348);

    auto tr_x_xyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 349);

    auto tr_x_xyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 350);

    auto tr_x_xyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 351);

    auto tr_x_xyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 352);

    auto tr_x_xyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 353);

    auto tr_x_xyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 354);

    auto tr_x_xyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 355);

    auto tr_x_xyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 356);

    auto tr_x_xyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 357);

    auto tr_x_xyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 358);

    auto tr_x_xyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 359);

    auto tr_x_xyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 360);

    auto tr_x_xyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 361);

    auto tr_x_xyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 362);

    auto tr_x_xyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 363);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             tr_x_xyy_xxxxxy,   \
                             tr_x_xyy_xxxxyy,   \
                             tr_x_xyy_xxxyyy,   \
                             tr_x_xyy_xxyyyy,   \
                             tr_x_xyy_xyyyyy,   \
                             tr_x_xyyz_xxxxxy,  \
                             tr_x_xyyz_xxxxyy,  \
                             tr_x_xyyz_xxxyyy,  \
                             tr_x_xyyz_xxyyyy,  \
                             tr_x_xyyz_xyyyyy,  \
                             tr_x_xyyzz_xxxxxx, \
                             tr_x_xyyzz_xxxxxy, \
                             tr_x_xyyzz_xxxxxz, \
                             tr_x_xyyzz_xxxxyy, \
                             tr_x_xyyzz_xxxxyz, \
                             tr_x_xyyzz_xxxxzz, \
                             tr_x_xyyzz_xxxyyy, \
                             tr_x_xyyzz_xxxyyz, \
                             tr_x_xyyzz_xxxyzz, \
                             tr_x_xyyzz_xxxzzz, \
                             tr_x_xyyzz_xxyyyy, \
                             tr_x_xyyzz_xxyyyz, \
                             tr_x_xyyzz_xxyyzz, \
                             tr_x_xyyzz_xxyzzz, \
                             tr_x_xyyzz_xxzzzz, \
                             tr_x_xyyzz_xyyyyy, \
                             tr_x_xyyzz_xyyyyz, \
                             tr_x_xyyzz_xyyyzz, \
                             tr_x_xyyzz_xyyzzz, \
                             tr_x_xyyzz_xyzzzz, \
                             tr_x_xyyzz_xzzzzz, \
                             tr_x_xyyzz_yyyyyy, \
                             tr_x_xyyzz_yyyyyz, \
                             tr_x_xyyzz_yyyyzz, \
                             tr_x_xyyzz_yyyzzz, \
                             tr_x_xyyzz_yyzzzz, \
                             tr_x_xyyzz_yzzzzz, \
                             tr_x_xyyzz_zzzzzz, \
                             tr_x_xyzz_xxxxxx,  \
                             tr_x_xyzz_xxxxxz,  \
                             tr_x_xyzz_xxxxzz,  \
                             tr_x_xyzz_xxxzzz,  \
                             tr_x_xyzz_xxzzzz,  \
                             tr_x_xyzz_xzzzzz,  \
                             tr_x_xzz_xxxxxx,   \
                             tr_x_xzz_xxxxxz,   \
                             tr_x_xzz_xxxxzz,   \
                             tr_x_xzz_xxxzzz,   \
                             tr_x_xzz_xxzzzz,   \
                             tr_x_xzz_xzzzzz,   \
                             tr_x_yyzz_xxxxyz,  \
                             tr_x_yyzz_xxxyyz,  \
                             tr_x_yyzz_xxxyz,   \
                             tr_x_yyzz_xxxyzz,  \
                             tr_x_yyzz_xxyyyz,  \
                             tr_x_yyzz_xxyyz,   \
                             tr_x_yyzz_xxyyzz,  \
                             tr_x_yyzz_xxyzz,   \
                             tr_x_yyzz_xxyzzz,  \
                             tr_x_yyzz_xyyyyz,  \
                             tr_x_yyzz_xyyyz,   \
                             tr_x_yyzz_xyyyzz,  \
                             tr_x_yyzz_xyyzz,   \
                             tr_x_yyzz_xyyzzz,  \
                             tr_x_yyzz_xyzzz,   \
                             tr_x_yyzz_xyzzzz,  \
                             tr_x_yyzz_yyyyyy,  \
                             tr_x_yyzz_yyyyyz,  \
                             tr_x_yyzz_yyyyz,   \
                             tr_x_yyzz_yyyyzz,  \
                             tr_x_yyzz_yyyzz,   \
                             tr_x_yyzz_yyyzzz,  \
                             tr_x_yyzz_yyzzz,   \
                             tr_x_yyzz_yyzzzz,  \
                             tr_x_yyzz_yzzzz,   \
                             tr_x_yyzz_yzzzzz,  \
                             tr_x_yyzz_zzzzzz,  \
                             ts_yyzz_xxxxyz,    \
                             ts_yyzz_xxxyyz,    \
                             ts_yyzz_xxxyzz,    \
                             ts_yyzz_xxyyyz,    \
                             ts_yyzz_xxyyzz,    \
                             ts_yyzz_xxyzzz,    \
                             ts_yyzz_xyyyyz,    \
                             ts_yyzz_xyyyzz,    \
                             ts_yyzz_xyyzzz,    \
                             ts_yyzz_xyzzzz,    \
                             ts_yyzz_yyyyyy,    \
                             ts_yyzz_yyyyyz,    \
                             ts_yyzz_yyyyzz,    \
                             ts_yyzz_yyyzzz,    \
                             ts_yyzz_yyzzzz,    \
                             ts_yyzz_yzzzzz,    \
                             ts_yyzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzz_xxxxxx[i] = tr_x_xzz_xxxxxx[i] * fe_0 + tr_x_xyzz_xxxxxx[i] * pa_y[i];

        tr_x_xyyzz_xxxxxy[i] = tr_x_xyy_xxxxxy[i] * fe_0 + tr_x_xyyz_xxxxxy[i] * pa_z[i];

        tr_x_xyyzz_xxxxxz[i] = tr_x_xzz_xxxxxz[i] * fe_0 + tr_x_xyzz_xxxxxz[i] * pa_y[i];

        tr_x_xyyzz_xxxxyy[i] = tr_x_xyy_xxxxyy[i] * fe_0 + tr_x_xyyz_xxxxyy[i] * pa_z[i];

        tr_x_xyyzz_xxxxyz[i] = 4.0 * tr_x_yyzz_xxxyz[i] * fe_0 + ts_yyzz_xxxxyz[i] * fe_0 + tr_x_yyzz_xxxxyz[i] * pa_x[i];

        tr_x_xyyzz_xxxxzz[i] = tr_x_xzz_xxxxzz[i] * fe_0 + tr_x_xyzz_xxxxzz[i] * pa_y[i];

        tr_x_xyyzz_xxxyyy[i] = tr_x_xyy_xxxyyy[i] * fe_0 + tr_x_xyyz_xxxyyy[i] * pa_z[i];

        tr_x_xyyzz_xxxyyz[i] = 3.0 * tr_x_yyzz_xxyyz[i] * fe_0 + ts_yyzz_xxxyyz[i] * fe_0 + tr_x_yyzz_xxxyyz[i] * pa_x[i];

        tr_x_xyyzz_xxxyzz[i] = 3.0 * tr_x_yyzz_xxyzz[i] * fe_0 + ts_yyzz_xxxyzz[i] * fe_0 + tr_x_yyzz_xxxyzz[i] * pa_x[i];

        tr_x_xyyzz_xxxzzz[i] = tr_x_xzz_xxxzzz[i] * fe_0 + tr_x_xyzz_xxxzzz[i] * pa_y[i];

        tr_x_xyyzz_xxyyyy[i] = tr_x_xyy_xxyyyy[i] * fe_0 + tr_x_xyyz_xxyyyy[i] * pa_z[i];

        tr_x_xyyzz_xxyyyz[i] = 2.0 * tr_x_yyzz_xyyyz[i] * fe_0 + ts_yyzz_xxyyyz[i] * fe_0 + tr_x_yyzz_xxyyyz[i] * pa_x[i];

        tr_x_xyyzz_xxyyzz[i] = 2.0 * tr_x_yyzz_xyyzz[i] * fe_0 + ts_yyzz_xxyyzz[i] * fe_0 + tr_x_yyzz_xxyyzz[i] * pa_x[i];

        tr_x_xyyzz_xxyzzz[i] = 2.0 * tr_x_yyzz_xyzzz[i] * fe_0 + ts_yyzz_xxyzzz[i] * fe_0 + tr_x_yyzz_xxyzzz[i] * pa_x[i];

        tr_x_xyyzz_xxzzzz[i] = tr_x_xzz_xxzzzz[i] * fe_0 + tr_x_xyzz_xxzzzz[i] * pa_y[i];

        tr_x_xyyzz_xyyyyy[i] = tr_x_xyy_xyyyyy[i] * fe_0 + tr_x_xyyz_xyyyyy[i] * pa_z[i];

        tr_x_xyyzz_xyyyyz[i] = tr_x_yyzz_yyyyz[i] * fe_0 + ts_yyzz_xyyyyz[i] * fe_0 + tr_x_yyzz_xyyyyz[i] * pa_x[i];

        tr_x_xyyzz_xyyyzz[i] = tr_x_yyzz_yyyzz[i] * fe_0 + ts_yyzz_xyyyzz[i] * fe_0 + tr_x_yyzz_xyyyzz[i] * pa_x[i];

        tr_x_xyyzz_xyyzzz[i] = tr_x_yyzz_yyzzz[i] * fe_0 + ts_yyzz_xyyzzz[i] * fe_0 + tr_x_yyzz_xyyzzz[i] * pa_x[i];

        tr_x_xyyzz_xyzzzz[i] = tr_x_yyzz_yzzzz[i] * fe_0 + ts_yyzz_xyzzzz[i] * fe_0 + tr_x_yyzz_xyzzzz[i] * pa_x[i];

        tr_x_xyyzz_xzzzzz[i] = tr_x_xzz_xzzzzz[i] * fe_0 + tr_x_xyzz_xzzzzz[i] * pa_y[i];

        tr_x_xyyzz_yyyyyy[i] = ts_yyzz_yyyyyy[i] * fe_0 + tr_x_yyzz_yyyyyy[i] * pa_x[i];

        tr_x_xyyzz_yyyyyz[i] = ts_yyzz_yyyyyz[i] * fe_0 + tr_x_yyzz_yyyyyz[i] * pa_x[i];

        tr_x_xyyzz_yyyyzz[i] = ts_yyzz_yyyyzz[i] * fe_0 + tr_x_yyzz_yyyyzz[i] * pa_x[i];

        tr_x_xyyzz_yyyzzz[i] = ts_yyzz_yyyzzz[i] * fe_0 + tr_x_yyzz_yyyzzz[i] * pa_x[i];

        tr_x_xyyzz_yyzzzz[i] = ts_yyzz_yyzzzz[i] * fe_0 + tr_x_yyzz_yyzzzz[i] * pa_x[i];

        tr_x_xyyzz_yzzzzz[i] = ts_yyzz_yzzzzz[i] * fe_0 + tr_x_yyzz_yzzzzz[i] * pa_x[i];

        tr_x_xyyzz_zzzzzz[i] = ts_yyzz_zzzzzz[i] * fe_0 + tr_x_yyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : HI

    auto tr_x_xyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 364);

    auto tr_x_xyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 365);

    auto tr_x_xyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 366);

    auto tr_x_xyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 367);

    auto tr_x_xyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 368);

    auto tr_x_xyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 369);

    auto tr_x_xyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 370);

    auto tr_x_xyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 371);

    auto tr_x_xyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 372);

    auto tr_x_xyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 373);

    auto tr_x_xyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 374);

    auto tr_x_xyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 375);

    auto tr_x_xyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 376);

    auto tr_x_xyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 377);

    auto tr_x_xyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 378);

    auto tr_x_xyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 379);

    auto tr_x_xyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 380);

    auto tr_x_xyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 381);

    auto tr_x_xyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 382);

    auto tr_x_xyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 383);

    auto tr_x_xyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 384);

    auto tr_x_xyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 385);

    auto tr_x_xyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 386);

    auto tr_x_xyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 387);

    auto tr_x_xyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 388);

    auto tr_x_xyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 389);

    auto tr_x_xyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 390);

    auto tr_x_xyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 391);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_x_xyzzz_xxxxxx, \
                             tr_x_xyzzz_xxxxxy, \
                             tr_x_xyzzz_xxxxxz, \
                             tr_x_xyzzz_xxxxyy, \
                             tr_x_xyzzz_xxxxyz, \
                             tr_x_xyzzz_xxxxzz, \
                             tr_x_xyzzz_xxxyyy, \
                             tr_x_xyzzz_xxxyyz, \
                             tr_x_xyzzz_xxxyzz, \
                             tr_x_xyzzz_xxxzzz, \
                             tr_x_xyzzz_xxyyyy, \
                             tr_x_xyzzz_xxyyyz, \
                             tr_x_xyzzz_xxyyzz, \
                             tr_x_xyzzz_xxyzzz, \
                             tr_x_xyzzz_xxzzzz, \
                             tr_x_xyzzz_xyyyyy, \
                             tr_x_xyzzz_xyyyyz, \
                             tr_x_xyzzz_xyyyzz, \
                             tr_x_xyzzz_xyyzzz, \
                             tr_x_xyzzz_xyzzzz, \
                             tr_x_xyzzz_xzzzzz, \
                             tr_x_xyzzz_yyyyyy, \
                             tr_x_xyzzz_yyyyyz, \
                             tr_x_xyzzz_yyyyzz, \
                             tr_x_xyzzz_yyyzzz, \
                             tr_x_xyzzz_yyzzzz, \
                             tr_x_xyzzz_yzzzzz, \
                             tr_x_xyzzz_zzzzzz, \
                             tr_x_xzzz_xxxxx,   \
                             tr_x_xzzz_xxxxxx,  \
                             tr_x_xzzz_xxxxxy,  \
                             tr_x_xzzz_xxxxxz,  \
                             tr_x_xzzz_xxxxy,   \
                             tr_x_xzzz_xxxxyy,  \
                             tr_x_xzzz_xxxxyz,  \
                             tr_x_xzzz_xxxxz,   \
                             tr_x_xzzz_xxxxzz,  \
                             tr_x_xzzz_xxxyy,   \
                             tr_x_xzzz_xxxyyy,  \
                             tr_x_xzzz_xxxyyz,  \
                             tr_x_xzzz_xxxyz,   \
                             tr_x_xzzz_xxxyzz,  \
                             tr_x_xzzz_xxxzz,   \
                             tr_x_xzzz_xxxzzz,  \
                             tr_x_xzzz_xxyyy,   \
                             tr_x_xzzz_xxyyyy,  \
                             tr_x_xzzz_xxyyyz,  \
                             tr_x_xzzz_xxyyz,   \
                             tr_x_xzzz_xxyyzz,  \
                             tr_x_xzzz_xxyzz,   \
                             tr_x_xzzz_xxyzzz,  \
                             tr_x_xzzz_xxzzz,   \
                             tr_x_xzzz_xxzzzz,  \
                             tr_x_xzzz_xyyyy,   \
                             tr_x_xzzz_xyyyyy,  \
                             tr_x_xzzz_xyyyyz,  \
                             tr_x_xzzz_xyyyz,   \
                             tr_x_xzzz_xyyyzz,  \
                             tr_x_xzzz_xyyzz,   \
                             tr_x_xzzz_xyyzzz,  \
                             tr_x_xzzz_xyzzz,   \
                             tr_x_xzzz_xyzzzz,  \
                             tr_x_xzzz_xzzzz,   \
                             tr_x_xzzz_xzzzzz,  \
                             tr_x_xzzz_zzzzzz,  \
                             tr_x_yzzz_yyyyyy,  \
                             tr_x_yzzz_yyyyyz,  \
                             tr_x_yzzz_yyyyzz,  \
                             tr_x_yzzz_yyyzzz,  \
                             tr_x_yzzz_yyzzzz,  \
                             tr_x_yzzz_yzzzzz,  \
                             ts_yzzz_yyyyyy,    \
                             ts_yzzz_yyyyyz,    \
                             ts_yzzz_yyyyzz,    \
                             ts_yzzz_yyyzzz,    \
                             ts_yzzz_yyzzzz,    \
                             ts_yzzz_yzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzz_xxxxxx[i] = tr_x_xzzz_xxxxxx[i] * pa_y[i];

        tr_x_xyzzz_xxxxxy[i] = tr_x_xzzz_xxxxx[i] * fe_0 + tr_x_xzzz_xxxxxy[i] * pa_y[i];

        tr_x_xyzzz_xxxxxz[i] = tr_x_xzzz_xxxxxz[i] * pa_y[i];

        tr_x_xyzzz_xxxxyy[i] = 2.0 * tr_x_xzzz_xxxxy[i] * fe_0 + tr_x_xzzz_xxxxyy[i] * pa_y[i];

        tr_x_xyzzz_xxxxyz[i] = tr_x_xzzz_xxxxz[i] * fe_0 + tr_x_xzzz_xxxxyz[i] * pa_y[i];

        tr_x_xyzzz_xxxxzz[i] = tr_x_xzzz_xxxxzz[i] * pa_y[i];

        tr_x_xyzzz_xxxyyy[i] = 3.0 * tr_x_xzzz_xxxyy[i] * fe_0 + tr_x_xzzz_xxxyyy[i] * pa_y[i];

        tr_x_xyzzz_xxxyyz[i] = 2.0 * tr_x_xzzz_xxxyz[i] * fe_0 + tr_x_xzzz_xxxyyz[i] * pa_y[i];

        tr_x_xyzzz_xxxyzz[i] = tr_x_xzzz_xxxzz[i] * fe_0 + tr_x_xzzz_xxxyzz[i] * pa_y[i];

        tr_x_xyzzz_xxxzzz[i] = tr_x_xzzz_xxxzzz[i] * pa_y[i];

        tr_x_xyzzz_xxyyyy[i] = 4.0 * tr_x_xzzz_xxyyy[i] * fe_0 + tr_x_xzzz_xxyyyy[i] * pa_y[i];

        tr_x_xyzzz_xxyyyz[i] = 3.0 * tr_x_xzzz_xxyyz[i] * fe_0 + tr_x_xzzz_xxyyyz[i] * pa_y[i];

        tr_x_xyzzz_xxyyzz[i] = 2.0 * tr_x_xzzz_xxyzz[i] * fe_0 + tr_x_xzzz_xxyyzz[i] * pa_y[i];

        tr_x_xyzzz_xxyzzz[i] = tr_x_xzzz_xxzzz[i] * fe_0 + tr_x_xzzz_xxyzzz[i] * pa_y[i];

        tr_x_xyzzz_xxzzzz[i] = tr_x_xzzz_xxzzzz[i] * pa_y[i];

        tr_x_xyzzz_xyyyyy[i] = 5.0 * tr_x_xzzz_xyyyy[i] * fe_0 + tr_x_xzzz_xyyyyy[i] * pa_y[i];

        tr_x_xyzzz_xyyyyz[i] = 4.0 * tr_x_xzzz_xyyyz[i] * fe_0 + tr_x_xzzz_xyyyyz[i] * pa_y[i];

        tr_x_xyzzz_xyyyzz[i] = 3.0 * tr_x_xzzz_xyyzz[i] * fe_0 + tr_x_xzzz_xyyyzz[i] * pa_y[i];

        tr_x_xyzzz_xyyzzz[i] = 2.0 * tr_x_xzzz_xyzzz[i] * fe_0 + tr_x_xzzz_xyyzzz[i] * pa_y[i];

        tr_x_xyzzz_xyzzzz[i] = tr_x_xzzz_xzzzz[i] * fe_0 + tr_x_xzzz_xyzzzz[i] * pa_y[i];

        tr_x_xyzzz_xzzzzz[i] = tr_x_xzzz_xzzzzz[i] * pa_y[i];

        tr_x_xyzzz_yyyyyy[i] = ts_yzzz_yyyyyy[i] * fe_0 + tr_x_yzzz_yyyyyy[i] * pa_x[i];

        tr_x_xyzzz_yyyyyz[i] = ts_yzzz_yyyyyz[i] * fe_0 + tr_x_yzzz_yyyyyz[i] * pa_x[i];

        tr_x_xyzzz_yyyyzz[i] = ts_yzzz_yyyyzz[i] * fe_0 + tr_x_yzzz_yyyyzz[i] * pa_x[i];

        tr_x_xyzzz_yyyzzz[i] = ts_yzzz_yyyzzz[i] * fe_0 + tr_x_yzzz_yyyzzz[i] * pa_x[i];

        tr_x_xyzzz_yyzzzz[i] = ts_yzzz_yyzzzz[i] * fe_0 + tr_x_yzzz_yyzzzz[i] * pa_x[i];

        tr_x_xyzzz_yzzzzz[i] = ts_yzzz_yzzzzz[i] * fe_0 + tr_x_yzzz_yzzzzz[i] * pa_x[i];

        tr_x_xyzzz_zzzzzz[i] = tr_x_xzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : HI

    auto tr_x_xzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 392);

    auto tr_x_xzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 393);

    auto tr_x_xzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 394);

    auto tr_x_xzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 395);

    auto tr_x_xzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 396);

    auto tr_x_xzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 397);

    auto tr_x_xzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 398);

    auto tr_x_xzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 399);

    auto tr_x_xzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 400);

    auto tr_x_xzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 401);

    auto tr_x_xzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 402);

    auto tr_x_xzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 403);

    auto tr_x_xzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 404);

    auto tr_x_xzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 405);

    auto tr_x_xzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 406);

    auto tr_x_xzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 407);

    auto tr_x_xzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 408);

    auto tr_x_xzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 409);

    auto tr_x_xzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 410);

    auto tr_x_xzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 411);

    auto tr_x_xzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 412);

    auto tr_x_xzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 413);

    auto tr_x_xzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 414);

    auto tr_x_xzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 415);

    auto tr_x_xzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 416);

    auto tr_x_xzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 417);

    auto tr_x_xzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 418);

    auto tr_x_xzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 419);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_x_xzz_xxxxxx,   \
                             tr_x_xzz_xxxxxy,   \
                             tr_x_xzz_xxxxyy,   \
                             tr_x_xzz_xxxyyy,   \
                             tr_x_xzz_xxyyyy,   \
                             tr_x_xzz_xyyyyy,   \
                             tr_x_xzzz_xxxxxx,  \
                             tr_x_xzzz_xxxxxy,  \
                             tr_x_xzzz_xxxxyy,  \
                             tr_x_xzzz_xxxyyy,  \
                             tr_x_xzzz_xxyyyy,  \
                             tr_x_xzzz_xyyyyy,  \
                             tr_x_xzzzz_xxxxxx, \
                             tr_x_xzzzz_xxxxxy, \
                             tr_x_xzzzz_xxxxxz, \
                             tr_x_xzzzz_xxxxyy, \
                             tr_x_xzzzz_xxxxyz, \
                             tr_x_xzzzz_xxxxzz, \
                             tr_x_xzzzz_xxxyyy, \
                             tr_x_xzzzz_xxxyyz, \
                             tr_x_xzzzz_xxxyzz, \
                             tr_x_xzzzz_xxxzzz, \
                             tr_x_xzzzz_xxyyyy, \
                             tr_x_xzzzz_xxyyyz, \
                             tr_x_xzzzz_xxyyzz, \
                             tr_x_xzzzz_xxyzzz, \
                             tr_x_xzzzz_xxzzzz, \
                             tr_x_xzzzz_xyyyyy, \
                             tr_x_xzzzz_xyyyyz, \
                             tr_x_xzzzz_xyyyzz, \
                             tr_x_xzzzz_xyyzzz, \
                             tr_x_xzzzz_xyzzzz, \
                             tr_x_xzzzz_xzzzzz, \
                             tr_x_xzzzz_yyyyyy, \
                             tr_x_xzzzz_yyyyyz, \
                             tr_x_xzzzz_yyyyzz, \
                             tr_x_xzzzz_yyyzzz, \
                             tr_x_xzzzz_yyzzzz, \
                             tr_x_xzzzz_yzzzzz, \
                             tr_x_xzzzz_zzzzzz, \
                             tr_x_zzzz_xxxxxz,  \
                             tr_x_zzzz_xxxxyz,  \
                             tr_x_zzzz_xxxxz,   \
                             tr_x_zzzz_xxxxzz,  \
                             tr_x_zzzz_xxxyyz,  \
                             tr_x_zzzz_xxxyz,   \
                             tr_x_zzzz_xxxyzz,  \
                             tr_x_zzzz_xxxzz,   \
                             tr_x_zzzz_xxxzzz,  \
                             tr_x_zzzz_xxyyyz,  \
                             tr_x_zzzz_xxyyz,   \
                             tr_x_zzzz_xxyyzz,  \
                             tr_x_zzzz_xxyzz,   \
                             tr_x_zzzz_xxyzzz,  \
                             tr_x_zzzz_xxzzz,   \
                             tr_x_zzzz_xxzzzz,  \
                             tr_x_zzzz_xyyyyz,  \
                             tr_x_zzzz_xyyyz,   \
                             tr_x_zzzz_xyyyzz,  \
                             tr_x_zzzz_xyyzz,   \
                             tr_x_zzzz_xyyzzz,  \
                             tr_x_zzzz_xyzzz,   \
                             tr_x_zzzz_xyzzzz,  \
                             tr_x_zzzz_xzzzz,   \
                             tr_x_zzzz_xzzzzz,  \
                             tr_x_zzzz_yyyyyy,  \
                             tr_x_zzzz_yyyyyz,  \
                             tr_x_zzzz_yyyyz,   \
                             tr_x_zzzz_yyyyzz,  \
                             tr_x_zzzz_yyyzz,   \
                             tr_x_zzzz_yyyzzz,  \
                             tr_x_zzzz_yyzzz,   \
                             tr_x_zzzz_yyzzzz,  \
                             tr_x_zzzz_yzzzz,   \
                             tr_x_zzzz_yzzzzz,  \
                             tr_x_zzzz_zzzzz,   \
                             tr_x_zzzz_zzzzzz,  \
                             ts_zzzz_xxxxxz,    \
                             ts_zzzz_xxxxyz,    \
                             ts_zzzz_xxxxzz,    \
                             ts_zzzz_xxxyyz,    \
                             ts_zzzz_xxxyzz,    \
                             ts_zzzz_xxxzzz,    \
                             ts_zzzz_xxyyyz,    \
                             ts_zzzz_xxyyzz,    \
                             ts_zzzz_xxyzzz,    \
                             ts_zzzz_xxzzzz,    \
                             ts_zzzz_xyyyyz,    \
                             ts_zzzz_xyyyzz,    \
                             ts_zzzz_xyyzzz,    \
                             ts_zzzz_xyzzzz,    \
                             ts_zzzz_xzzzzz,    \
                             ts_zzzz_yyyyyy,    \
                             ts_zzzz_yyyyyz,    \
                             ts_zzzz_yyyyzz,    \
                             ts_zzzz_yyyzzz,    \
                             ts_zzzz_yyzzzz,    \
                             ts_zzzz_yzzzzz,    \
                             ts_zzzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzz_xxxxxx[i] = 3.0 * tr_x_xzz_xxxxxx[i] * fe_0 + tr_x_xzzz_xxxxxx[i] * pa_z[i];

        tr_x_xzzzz_xxxxxy[i] = 3.0 * tr_x_xzz_xxxxxy[i] * fe_0 + tr_x_xzzz_xxxxxy[i] * pa_z[i];

        tr_x_xzzzz_xxxxxz[i] = 5.0 * tr_x_zzzz_xxxxz[i] * fe_0 + ts_zzzz_xxxxxz[i] * fe_0 + tr_x_zzzz_xxxxxz[i] * pa_x[i];

        tr_x_xzzzz_xxxxyy[i] = 3.0 * tr_x_xzz_xxxxyy[i] * fe_0 + tr_x_xzzz_xxxxyy[i] * pa_z[i];

        tr_x_xzzzz_xxxxyz[i] = 4.0 * tr_x_zzzz_xxxyz[i] * fe_0 + ts_zzzz_xxxxyz[i] * fe_0 + tr_x_zzzz_xxxxyz[i] * pa_x[i];

        tr_x_xzzzz_xxxxzz[i] = 4.0 * tr_x_zzzz_xxxzz[i] * fe_0 + ts_zzzz_xxxxzz[i] * fe_0 + tr_x_zzzz_xxxxzz[i] * pa_x[i];

        tr_x_xzzzz_xxxyyy[i] = 3.0 * tr_x_xzz_xxxyyy[i] * fe_0 + tr_x_xzzz_xxxyyy[i] * pa_z[i];

        tr_x_xzzzz_xxxyyz[i] = 3.0 * tr_x_zzzz_xxyyz[i] * fe_0 + ts_zzzz_xxxyyz[i] * fe_0 + tr_x_zzzz_xxxyyz[i] * pa_x[i];

        tr_x_xzzzz_xxxyzz[i] = 3.0 * tr_x_zzzz_xxyzz[i] * fe_0 + ts_zzzz_xxxyzz[i] * fe_0 + tr_x_zzzz_xxxyzz[i] * pa_x[i];

        tr_x_xzzzz_xxxzzz[i] = 3.0 * tr_x_zzzz_xxzzz[i] * fe_0 + ts_zzzz_xxxzzz[i] * fe_0 + tr_x_zzzz_xxxzzz[i] * pa_x[i];

        tr_x_xzzzz_xxyyyy[i] = 3.0 * tr_x_xzz_xxyyyy[i] * fe_0 + tr_x_xzzz_xxyyyy[i] * pa_z[i];

        tr_x_xzzzz_xxyyyz[i] = 2.0 * tr_x_zzzz_xyyyz[i] * fe_0 + ts_zzzz_xxyyyz[i] * fe_0 + tr_x_zzzz_xxyyyz[i] * pa_x[i];

        tr_x_xzzzz_xxyyzz[i] = 2.0 * tr_x_zzzz_xyyzz[i] * fe_0 + ts_zzzz_xxyyzz[i] * fe_0 + tr_x_zzzz_xxyyzz[i] * pa_x[i];

        tr_x_xzzzz_xxyzzz[i] = 2.0 * tr_x_zzzz_xyzzz[i] * fe_0 + ts_zzzz_xxyzzz[i] * fe_0 + tr_x_zzzz_xxyzzz[i] * pa_x[i];

        tr_x_xzzzz_xxzzzz[i] = 2.0 * tr_x_zzzz_xzzzz[i] * fe_0 + ts_zzzz_xxzzzz[i] * fe_0 + tr_x_zzzz_xxzzzz[i] * pa_x[i];

        tr_x_xzzzz_xyyyyy[i] = 3.0 * tr_x_xzz_xyyyyy[i] * fe_0 + tr_x_xzzz_xyyyyy[i] * pa_z[i];

        tr_x_xzzzz_xyyyyz[i] = tr_x_zzzz_yyyyz[i] * fe_0 + ts_zzzz_xyyyyz[i] * fe_0 + tr_x_zzzz_xyyyyz[i] * pa_x[i];

        tr_x_xzzzz_xyyyzz[i] = tr_x_zzzz_yyyzz[i] * fe_0 + ts_zzzz_xyyyzz[i] * fe_0 + tr_x_zzzz_xyyyzz[i] * pa_x[i];

        tr_x_xzzzz_xyyzzz[i] = tr_x_zzzz_yyzzz[i] * fe_0 + ts_zzzz_xyyzzz[i] * fe_0 + tr_x_zzzz_xyyzzz[i] * pa_x[i];

        tr_x_xzzzz_xyzzzz[i] = tr_x_zzzz_yzzzz[i] * fe_0 + ts_zzzz_xyzzzz[i] * fe_0 + tr_x_zzzz_xyzzzz[i] * pa_x[i];

        tr_x_xzzzz_xzzzzz[i] = tr_x_zzzz_zzzzz[i] * fe_0 + ts_zzzz_xzzzzz[i] * fe_0 + tr_x_zzzz_xzzzzz[i] * pa_x[i];

        tr_x_xzzzz_yyyyyy[i] = ts_zzzz_yyyyyy[i] * fe_0 + tr_x_zzzz_yyyyyy[i] * pa_x[i];

        tr_x_xzzzz_yyyyyz[i] = ts_zzzz_yyyyyz[i] * fe_0 + tr_x_zzzz_yyyyyz[i] * pa_x[i];

        tr_x_xzzzz_yyyyzz[i] = ts_zzzz_yyyyzz[i] * fe_0 + tr_x_zzzz_yyyyzz[i] * pa_x[i];

        tr_x_xzzzz_yyyzzz[i] = ts_zzzz_yyyzzz[i] * fe_0 + tr_x_zzzz_yyyzzz[i] * pa_x[i];

        tr_x_xzzzz_yyzzzz[i] = ts_zzzz_yyzzzz[i] * fe_0 + tr_x_zzzz_yyzzzz[i] * pa_x[i];

        tr_x_xzzzz_yzzzzz[i] = ts_zzzz_yzzzzz[i] * fe_0 + tr_x_zzzz_yzzzzz[i] * pa_x[i];

        tr_x_xzzzz_zzzzzz[i] = ts_zzzz_zzzzzz[i] * fe_0 + tr_x_zzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : HI

    auto tr_x_yyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 420);

    auto tr_x_yyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 421);

    auto tr_x_yyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 422);

    auto tr_x_yyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 423);

    auto tr_x_yyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 424);

    auto tr_x_yyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 425);

    auto tr_x_yyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 426);

    auto tr_x_yyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 427);

    auto tr_x_yyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 428);

    auto tr_x_yyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 429);

    auto tr_x_yyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 430);

    auto tr_x_yyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 431);

    auto tr_x_yyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 432);

    auto tr_x_yyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 433);

    auto tr_x_yyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 434);

    auto tr_x_yyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 435);

    auto tr_x_yyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 436);

    auto tr_x_yyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 437);

    auto tr_x_yyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 438);

    auto tr_x_yyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 439);

    auto tr_x_yyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 440);

    auto tr_x_yyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 441);

    auto tr_x_yyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 442);

    auto tr_x_yyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 443);

    auto tr_x_yyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 444);

    auto tr_x_yyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 445);

    auto tr_x_yyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 446);

    auto tr_x_yyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 447);

#pragma omp simd aligned(pa_y,                  \
                             tr_x_yyy_xxxxxx,   \
                             tr_x_yyy_xxxxxy,   \
                             tr_x_yyy_xxxxxz,   \
                             tr_x_yyy_xxxxyy,   \
                             tr_x_yyy_xxxxyz,   \
                             tr_x_yyy_xxxxzz,   \
                             tr_x_yyy_xxxyyy,   \
                             tr_x_yyy_xxxyyz,   \
                             tr_x_yyy_xxxyzz,   \
                             tr_x_yyy_xxxzzz,   \
                             tr_x_yyy_xxyyyy,   \
                             tr_x_yyy_xxyyyz,   \
                             tr_x_yyy_xxyyzz,   \
                             tr_x_yyy_xxyzzz,   \
                             tr_x_yyy_xxzzzz,   \
                             tr_x_yyy_xyyyyy,   \
                             tr_x_yyy_xyyyyz,   \
                             tr_x_yyy_xyyyzz,   \
                             tr_x_yyy_xyyzzz,   \
                             tr_x_yyy_xyzzzz,   \
                             tr_x_yyy_xzzzzz,   \
                             tr_x_yyy_yyyyyy,   \
                             tr_x_yyy_yyyyyz,   \
                             tr_x_yyy_yyyyzz,   \
                             tr_x_yyy_yyyzzz,   \
                             tr_x_yyy_yyzzzz,   \
                             tr_x_yyy_yzzzzz,   \
                             tr_x_yyy_zzzzzz,   \
                             tr_x_yyyy_xxxxx,   \
                             tr_x_yyyy_xxxxxx,  \
                             tr_x_yyyy_xxxxxy,  \
                             tr_x_yyyy_xxxxxz,  \
                             tr_x_yyyy_xxxxy,   \
                             tr_x_yyyy_xxxxyy,  \
                             tr_x_yyyy_xxxxyz,  \
                             tr_x_yyyy_xxxxz,   \
                             tr_x_yyyy_xxxxzz,  \
                             tr_x_yyyy_xxxyy,   \
                             tr_x_yyyy_xxxyyy,  \
                             tr_x_yyyy_xxxyyz,  \
                             tr_x_yyyy_xxxyz,   \
                             tr_x_yyyy_xxxyzz,  \
                             tr_x_yyyy_xxxzz,   \
                             tr_x_yyyy_xxxzzz,  \
                             tr_x_yyyy_xxyyy,   \
                             tr_x_yyyy_xxyyyy,  \
                             tr_x_yyyy_xxyyyz,  \
                             tr_x_yyyy_xxyyz,   \
                             tr_x_yyyy_xxyyzz,  \
                             tr_x_yyyy_xxyzz,   \
                             tr_x_yyyy_xxyzzz,  \
                             tr_x_yyyy_xxzzz,   \
                             tr_x_yyyy_xxzzzz,  \
                             tr_x_yyyy_xyyyy,   \
                             tr_x_yyyy_xyyyyy,  \
                             tr_x_yyyy_xyyyyz,  \
                             tr_x_yyyy_xyyyz,   \
                             tr_x_yyyy_xyyyzz,  \
                             tr_x_yyyy_xyyzz,   \
                             tr_x_yyyy_xyyzzz,  \
                             tr_x_yyyy_xyzzz,   \
                             tr_x_yyyy_xyzzzz,  \
                             tr_x_yyyy_xzzzz,   \
                             tr_x_yyyy_xzzzzz,  \
                             tr_x_yyyy_yyyyy,   \
                             tr_x_yyyy_yyyyyy,  \
                             tr_x_yyyy_yyyyyz,  \
                             tr_x_yyyy_yyyyz,   \
                             tr_x_yyyy_yyyyzz,  \
                             tr_x_yyyy_yyyzz,   \
                             tr_x_yyyy_yyyzzz,  \
                             tr_x_yyyy_yyzzz,   \
                             tr_x_yyyy_yyzzzz,  \
                             tr_x_yyyy_yzzzz,   \
                             tr_x_yyyy_yzzzzz,  \
                             tr_x_yyyy_zzzzz,   \
                             tr_x_yyyy_zzzzzz,  \
                             tr_x_yyyyy_xxxxxx, \
                             tr_x_yyyyy_xxxxxy, \
                             tr_x_yyyyy_xxxxxz, \
                             tr_x_yyyyy_xxxxyy, \
                             tr_x_yyyyy_xxxxyz, \
                             tr_x_yyyyy_xxxxzz, \
                             tr_x_yyyyy_xxxyyy, \
                             tr_x_yyyyy_xxxyyz, \
                             tr_x_yyyyy_xxxyzz, \
                             tr_x_yyyyy_xxxzzz, \
                             tr_x_yyyyy_xxyyyy, \
                             tr_x_yyyyy_xxyyyz, \
                             tr_x_yyyyy_xxyyzz, \
                             tr_x_yyyyy_xxyzzz, \
                             tr_x_yyyyy_xxzzzz, \
                             tr_x_yyyyy_xyyyyy, \
                             tr_x_yyyyy_xyyyyz, \
                             tr_x_yyyyy_xyyyzz, \
                             tr_x_yyyyy_xyyzzz, \
                             tr_x_yyyyy_xyzzzz, \
                             tr_x_yyyyy_xzzzzz, \
                             tr_x_yyyyy_yyyyyy, \
                             tr_x_yyyyy_yyyyyz, \
                             tr_x_yyyyy_yyyyzz, \
                             tr_x_yyyyy_yyyzzz, \
                             tr_x_yyyyy_yyzzzz, \
                             tr_x_yyyyy_yzzzzz, \
                             tr_x_yyyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyy_xxxxxx[i] = 4.0 * tr_x_yyy_xxxxxx[i] * fe_0 + tr_x_yyyy_xxxxxx[i] * pa_y[i];

        tr_x_yyyyy_xxxxxy[i] = 4.0 * tr_x_yyy_xxxxxy[i] * fe_0 + tr_x_yyyy_xxxxx[i] * fe_0 + tr_x_yyyy_xxxxxy[i] * pa_y[i];

        tr_x_yyyyy_xxxxxz[i] = 4.0 * tr_x_yyy_xxxxxz[i] * fe_0 + tr_x_yyyy_xxxxxz[i] * pa_y[i];

        tr_x_yyyyy_xxxxyy[i] = 4.0 * tr_x_yyy_xxxxyy[i] * fe_0 + 2.0 * tr_x_yyyy_xxxxy[i] * fe_0 + tr_x_yyyy_xxxxyy[i] * pa_y[i];

        tr_x_yyyyy_xxxxyz[i] = 4.0 * tr_x_yyy_xxxxyz[i] * fe_0 + tr_x_yyyy_xxxxz[i] * fe_0 + tr_x_yyyy_xxxxyz[i] * pa_y[i];

        tr_x_yyyyy_xxxxzz[i] = 4.0 * tr_x_yyy_xxxxzz[i] * fe_0 + tr_x_yyyy_xxxxzz[i] * pa_y[i];

        tr_x_yyyyy_xxxyyy[i] = 4.0 * tr_x_yyy_xxxyyy[i] * fe_0 + 3.0 * tr_x_yyyy_xxxyy[i] * fe_0 + tr_x_yyyy_xxxyyy[i] * pa_y[i];

        tr_x_yyyyy_xxxyyz[i] = 4.0 * tr_x_yyy_xxxyyz[i] * fe_0 + 2.0 * tr_x_yyyy_xxxyz[i] * fe_0 + tr_x_yyyy_xxxyyz[i] * pa_y[i];

        tr_x_yyyyy_xxxyzz[i] = 4.0 * tr_x_yyy_xxxyzz[i] * fe_0 + tr_x_yyyy_xxxzz[i] * fe_0 + tr_x_yyyy_xxxyzz[i] * pa_y[i];

        tr_x_yyyyy_xxxzzz[i] = 4.0 * tr_x_yyy_xxxzzz[i] * fe_0 + tr_x_yyyy_xxxzzz[i] * pa_y[i];

        tr_x_yyyyy_xxyyyy[i] = 4.0 * tr_x_yyy_xxyyyy[i] * fe_0 + 4.0 * tr_x_yyyy_xxyyy[i] * fe_0 + tr_x_yyyy_xxyyyy[i] * pa_y[i];

        tr_x_yyyyy_xxyyyz[i] = 4.0 * tr_x_yyy_xxyyyz[i] * fe_0 + 3.0 * tr_x_yyyy_xxyyz[i] * fe_0 + tr_x_yyyy_xxyyyz[i] * pa_y[i];

        tr_x_yyyyy_xxyyzz[i] = 4.0 * tr_x_yyy_xxyyzz[i] * fe_0 + 2.0 * tr_x_yyyy_xxyzz[i] * fe_0 + tr_x_yyyy_xxyyzz[i] * pa_y[i];

        tr_x_yyyyy_xxyzzz[i] = 4.0 * tr_x_yyy_xxyzzz[i] * fe_0 + tr_x_yyyy_xxzzz[i] * fe_0 + tr_x_yyyy_xxyzzz[i] * pa_y[i];

        tr_x_yyyyy_xxzzzz[i] = 4.0 * tr_x_yyy_xxzzzz[i] * fe_0 + tr_x_yyyy_xxzzzz[i] * pa_y[i];

        tr_x_yyyyy_xyyyyy[i] = 4.0 * tr_x_yyy_xyyyyy[i] * fe_0 + 5.0 * tr_x_yyyy_xyyyy[i] * fe_0 + tr_x_yyyy_xyyyyy[i] * pa_y[i];

        tr_x_yyyyy_xyyyyz[i] = 4.0 * tr_x_yyy_xyyyyz[i] * fe_0 + 4.0 * tr_x_yyyy_xyyyz[i] * fe_0 + tr_x_yyyy_xyyyyz[i] * pa_y[i];

        tr_x_yyyyy_xyyyzz[i] = 4.0 * tr_x_yyy_xyyyzz[i] * fe_0 + 3.0 * tr_x_yyyy_xyyzz[i] * fe_0 + tr_x_yyyy_xyyyzz[i] * pa_y[i];

        tr_x_yyyyy_xyyzzz[i] = 4.0 * tr_x_yyy_xyyzzz[i] * fe_0 + 2.0 * tr_x_yyyy_xyzzz[i] * fe_0 + tr_x_yyyy_xyyzzz[i] * pa_y[i];

        tr_x_yyyyy_xyzzzz[i] = 4.0 * tr_x_yyy_xyzzzz[i] * fe_0 + tr_x_yyyy_xzzzz[i] * fe_0 + tr_x_yyyy_xyzzzz[i] * pa_y[i];

        tr_x_yyyyy_xzzzzz[i] = 4.0 * tr_x_yyy_xzzzzz[i] * fe_0 + tr_x_yyyy_xzzzzz[i] * pa_y[i];

        tr_x_yyyyy_yyyyyy[i] = 4.0 * tr_x_yyy_yyyyyy[i] * fe_0 + 6.0 * tr_x_yyyy_yyyyy[i] * fe_0 + tr_x_yyyy_yyyyyy[i] * pa_y[i];

        tr_x_yyyyy_yyyyyz[i] = 4.0 * tr_x_yyy_yyyyyz[i] * fe_0 + 5.0 * tr_x_yyyy_yyyyz[i] * fe_0 + tr_x_yyyy_yyyyyz[i] * pa_y[i];

        tr_x_yyyyy_yyyyzz[i] = 4.0 * tr_x_yyy_yyyyzz[i] * fe_0 + 4.0 * tr_x_yyyy_yyyzz[i] * fe_0 + tr_x_yyyy_yyyyzz[i] * pa_y[i];

        tr_x_yyyyy_yyyzzz[i] = 4.0 * tr_x_yyy_yyyzzz[i] * fe_0 + 3.0 * tr_x_yyyy_yyzzz[i] * fe_0 + tr_x_yyyy_yyyzzz[i] * pa_y[i];

        tr_x_yyyyy_yyzzzz[i] = 4.0 * tr_x_yyy_yyzzzz[i] * fe_0 + 2.0 * tr_x_yyyy_yzzzz[i] * fe_0 + tr_x_yyyy_yyzzzz[i] * pa_y[i];

        tr_x_yyyyy_yzzzzz[i] = 4.0 * tr_x_yyy_yzzzzz[i] * fe_0 + tr_x_yyyy_zzzzz[i] * fe_0 + tr_x_yyyy_yzzzzz[i] * pa_y[i];

        tr_x_yyyyy_zzzzzz[i] = 4.0 * tr_x_yyy_zzzzzz[i] * fe_0 + tr_x_yyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 448-476 components of targeted buffer : HI

    auto tr_x_yyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 448);

    auto tr_x_yyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 449);

    auto tr_x_yyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 450);

    auto tr_x_yyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 451);

    auto tr_x_yyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 452);

    auto tr_x_yyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 453);

    auto tr_x_yyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 454);

    auto tr_x_yyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 455);

    auto tr_x_yyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 456);

    auto tr_x_yyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 457);

    auto tr_x_yyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 458);

    auto tr_x_yyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 459);

    auto tr_x_yyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 460);

    auto tr_x_yyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 461);

    auto tr_x_yyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 462);

    auto tr_x_yyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 463);

    auto tr_x_yyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 464);

    auto tr_x_yyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 465);

    auto tr_x_yyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 466);

    auto tr_x_yyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 467);

    auto tr_x_yyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 468);

    auto tr_x_yyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 469);

    auto tr_x_yyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 470);

    auto tr_x_yyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 471);

    auto tr_x_yyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 472);

    auto tr_x_yyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 473);

    auto tr_x_yyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 474);

    auto tr_x_yyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 475);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_x_yyyy_xxxxxx,  \
                             tr_x_yyyy_xxxxxy,  \
                             tr_x_yyyy_xxxxy,   \
                             tr_x_yyyy_xxxxyy,  \
                             tr_x_yyyy_xxxxyz,  \
                             tr_x_yyyy_xxxyy,   \
                             tr_x_yyyy_xxxyyy,  \
                             tr_x_yyyy_xxxyyz,  \
                             tr_x_yyyy_xxxyz,   \
                             tr_x_yyyy_xxxyzz,  \
                             tr_x_yyyy_xxyyy,   \
                             tr_x_yyyy_xxyyyy,  \
                             tr_x_yyyy_xxyyyz,  \
                             tr_x_yyyy_xxyyz,   \
                             tr_x_yyyy_xxyyzz,  \
                             tr_x_yyyy_xxyzz,   \
                             tr_x_yyyy_xxyzzz,  \
                             tr_x_yyyy_xyyyy,   \
                             tr_x_yyyy_xyyyyy,  \
                             tr_x_yyyy_xyyyyz,  \
                             tr_x_yyyy_xyyyz,   \
                             tr_x_yyyy_xyyyzz,  \
                             tr_x_yyyy_xyyzz,   \
                             tr_x_yyyy_xyyzzz,  \
                             tr_x_yyyy_xyzzz,   \
                             tr_x_yyyy_xyzzzz,  \
                             tr_x_yyyy_yyyyy,   \
                             tr_x_yyyy_yyyyyy,  \
                             tr_x_yyyy_yyyyyz,  \
                             tr_x_yyyy_yyyyz,   \
                             tr_x_yyyy_yyyyzz,  \
                             tr_x_yyyy_yyyzz,   \
                             tr_x_yyyy_yyyzzz,  \
                             tr_x_yyyy_yyzzz,   \
                             tr_x_yyyy_yyzzzz,  \
                             tr_x_yyyy_yzzzz,   \
                             tr_x_yyyy_yzzzzz,  \
                             tr_x_yyyyz_xxxxxx, \
                             tr_x_yyyyz_xxxxxy, \
                             tr_x_yyyyz_xxxxxz, \
                             tr_x_yyyyz_xxxxyy, \
                             tr_x_yyyyz_xxxxyz, \
                             tr_x_yyyyz_xxxxzz, \
                             tr_x_yyyyz_xxxyyy, \
                             tr_x_yyyyz_xxxyyz, \
                             tr_x_yyyyz_xxxyzz, \
                             tr_x_yyyyz_xxxzzz, \
                             tr_x_yyyyz_xxyyyy, \
                             tr_x_yyyyz_xxyyyz, \
                             tr_x_yyyyz_xxyyzz, \
                             tr_x_yyyyz_xxyzzz, \
                             tr_x_yyyyz_xxzzzz, \
                             tr_x_yyyyz_xyyyyy, \
                             tr_x_yyyyz_xyyyyz, \
                             tr_x_yyyyz_xyyyzz, \
                             tr_x_yyyyz_xyyzzz, \
                             tr_x_yyyyz_xyzzzz, \
                             tr_x_yyyyz_xzzzzz, \
                             tr_x_yyyyz_yyyyyy, \
                             tr_x_yyyyz_yyyyyz, \
                             tr_x_yyyyz_yyyyzz, \
                             tr_x_yyyyz_yyyzzz, \
                             tr_x_yyyyz_yyzzzz, \
                             tr_x_yyyyz_yzzzzz, \
                             tr_x_yyyyz_zzzzzz, \
                             tr_x_yyyz_xxxxxz,  \
                             tr_x_yyyz_xxxxzz,  \
                             tr_x_yyyz_xxxzzz,  \
                             tr_x_yyyz_xxzzzz,  \
                             tr_x_yyyz_xzzzzz,  \
                             tr_x_yyyz_zzzzzz,  \
                             tr_x_yyz_xxxxxz,   \
                             tr_x_yyz_xxxxzz,   \
                             tr_x_yyz_xxxzzz,   \
                             tr_x_yyz_xxzzzz,   \
                             tr_x_yyz_xzzzzz,   \
                             tr_x_yyz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyz_xxxxxx[i] = tr_x_yyyy_xxxxxx[i] * pa_z[i];

        tr_x_yyyyz_xxxxxy[i] = tr_x_yyyy_xxxxxy[i] * pa_z[i];

        tr_x_yyyyz_xxxxxz[i] = 3.0 * tr_x_yyz_xxxxxz[i] * fe_0 + tr_x_yyyz_xxxxxz[i] * pa_y[i];

        tr_x_yyyyz_xxxxyy[i] = tr_x_yyyy_xxxxyy[i] * pa_z[i];

        tr_x_yyyyz_xxxxyz[i] = tr_x_yyyy_xxxxy[i] * fe_0 + tr_x_yyyy_xxxxyz[i] * pa_z[i];

        tr_x_yyyyz_xxxxzz[i] = 3.0 * tr_x_yyz_xxxxzz[i] * fe_0 + tr_x_yyyz_xxxxzz[i] * pa_y[i];

        tr_x_yyyyz_xxxyyy[i] = tr_x_yyyy_xxxyyy[i] * pa_z[i];

        tr_x_yyyyz_xxxyyz[i] = tr_x_yyyy_xxxyy[i] * fe_0 + tr_x_yyyy_xxxyyz[i] * pa_z[i];

        tr_x_yyyyz_xxxyzz[i] = 2.0 * tr_x_yyyy_xxxyz[i] * fe_0 + tr_x_yyyy_xxxyzz[i] * pa_z[i];

        tr_x_yyyyz_xxxzzz[i] = 3.0 * tr_x_yyz_xxxzzz[i] * fe_0 + tr_x_yyyz_xxxzzz[i] * pa_y[i];

        tr_x_yyyyz_xxyyyy[i] = tr_x_yyyy_xxyyyy[i] * pa_z[i];

        tr_x_yyyyz_xxyyyz[i] = tr_x_yyyy_xxyyy[i] * fe_0 + tr_x_yyyy_xxyyyz[i] * pa_z[i];

        tr_x_yyyyz_xxyyzz[i] = 2.0 * tr_x_yyyy_xxyyz[i] * fe_0 + tr_x_yyyy_xxyyzz[i] * pa_z[i];

        tr_x_yyyyz_xxyzzz[i] = 3.0 * tr_x_yyyy_xxyzz[i] * fe_0 + tr_x_yyyy_xxyzzz[i] * pa_z[i];

        tr_x_yyyyz_xxzzzz[i] = 3.0 * tr_x_yyz_xxzzzz[i] * fe_0 + tr_x_yyyz_xxzzzz[i] * pa_y[i];

        tr_x_yyyyz_xyyyyy[i] = tr_x_yyyy_xyyyyy[i] * pa_z[i];

        tr_x_yyyyz_xyyyyz[i] = tr_x_yyyy_xyyyy[i] * fe_0 + tr_x_yyyy_xyyyyz[i] * pa_z[i];

        tr_x_yyyyz_xyyyzz[i] = 2.0 * tr_x_yyyy_xyyyz[i] * fe_0 + tr_x_yyyy_xyyyzz[i] * pa_z[i];

        tr_x_yyyyz_xyyzzz[i] = 3.0 * tr_x_yyyy_xyyzz[i] * fe_0 + tr_x_yyyy_xyyzzz[i] * pa_z[i];

        tr_x_yyyyz_xyzzzz[i] = 4.0 * tr_x_yyyy_xyzzz[i] * fe_0 + tr_x_yyyy_xyzzzz[i] * pa_z[i];

        tr_x_yyyyz_xzzzzz[i] = 3.0 * tr_x_yyz_xzzzzz[i] * fe_0 + tr_x_yyyz_xzzzzz[i] * pa_y[i];

        tr_x_yyyyz_yyyyyy[i] = tr_x_yyyy_yyyyyy[i] * pa_z[i];

        tr_x_yyyyz_yyyyyz[i] = tr_x_yyyy_yyyyy[i] * fe_0 + tr_x_yyyy_yyyyyz[i] * pa_z[i];

        tr_x_yyyyz_yyyyzz[i] = 2.0 * tr_x_yyyy_yyyyz[i] * fe_0 + tr_x_yyyy_yyyyzz[i] * pa_z[i];

        tr_x_yyyyz_yyyzzz[i] = 3.0 * tr_x_yyyy_yyyzz[i] * fe_0 + tr_x_yyyy_yyyzzz[i] * pa_z[i];

        tr_x_yyyyz_yyzzzz[i] = 4.0 * tr_x_yyyy_yyzzz[i] * fe_0 + tr_x_yyyy_yyzzzz[i] * pa_z[i];

        tr_x_yyyyz_yzzzzz[i] = 5.0 * tr_x_yyyy_yzzzz[i] * fe_0 + tr_x_yyyy_yzzzzz[i] * pa_z[i];

        tr_x_yyyyz_zzzzzz[i] = 3.0 * tr_x_yyz_zzzzzz[i] * fe_0 + tr_x_yyyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 476-504 components of targeted buffer : HI

    auto tr_x_yyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 476);

    auto tr_x_yyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 477);

    auto tr_x_yyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 478);

    auto tr_x_yyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 479);

    auto tr_x_yyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 480);

    auto tr_x_yyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 481);

    auto tr_x_yyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 482);

    auto tr_x_yyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 483);

    auto tr_x_yyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 484);

    auto tr_x_yyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 485);

    auto tr_x_yyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 486);

    auto tr_x_yyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 487);

    auto tr_x_yyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 488);

    auto tr_x_yyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 489);

    auto tr_x_yyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 490);

    auto tr_x_yyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 491);

    auto tr_x_yyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 492);

    auto tr_x_yyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 493);

    auto tr_x_yyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 494);

    auto tr_x_yyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 495);

    auto tr_x_yyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 496);

    auto tr_x_yyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 497);

    auto tr_x_yyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 498);

    auto tr_x_yyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 499);

    auto tr_x_yyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 500);

    auto tr_x_yyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 501);

    auto tr_x_yyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 502);

    auto tr_x_yyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 503);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_x_yyy_xxxxxy,   \
                             tr_x_yyy_xxxxyy,   \
                             tr_x_yyy_xxxyyy,   \
                             tr_x_yyy_xxyyyy,   \
                             tr_x_yyy_xyyyyy,   \
                             tr_x_yyy_yyyyyy,   \
                             tr_x_yyyz_xxxxxy,  \
                             tr_x_yyyz_xxxxyy,  \
                             tr_x_yyyz_xxxyyy,  \
                             tr_x_yyyz_xxyyyy,  \
                             tr_x_yyyz_xyyyyy,  \
                             tr_x_yyyz_yyyyyy,  \
                             tr_x_yyyzz_xxxxxx, \
                             tr_x_yyyzz_xxxxxy, \
                             tr_x_yyyzz_xxxxxz, \
                             tr_x_yyyzz_xxxxyy, \
                             tr_x_yyyzz_xxxxyz, \
                             tr_x_yyyzz_xxxxzz, \
                             tr_x_yyyzz_xxxyyy, \
                             tr_x_yyyzz_xxxyyz, \
                             tr_x_yyyzz_xxxyzz, \
                             tr_x_yyyzz_xxxzzz, \
                             tr_x_yyyzz_xxyyyy, \
                             tr_x_yyyzz_xxyyyz, \
                             tr_x_yyyzz_xxyyzz, \
                             tr_x_yyyzz_xxyzzz, \
                             tr_x_yyyzz_xxzzzz, \
                             tr_x_yyyzz_xyyyyy, \
                             tr_x_yyyzz_xyyyyz, \
                             tr_x_yyyzz_xyyyzz, \
                             tr_x_yyyzz_xyyzzz, \
                             tr_x_yyyzz_xyzzzz, \
                             tr_x_yyyzz_xzzzzz, \
                             tr_x_yyyzz_yyyyyy, \
                             tr_x_yyyzz_yyyyyz, \
                             tr_x_yyyzz_yyyyzz, \
                             tr_x_yyyzz_yyyzzz, \
                             tr_x_yyyzz_yyzzzz, \
                             tr_x_yyyzz_yzzzzz, \
                             tr_x_yyyzz_zzzzzz, \
                             tr_x_yyzz_xxxxxx,  \
                             tr_x_yyzz_xxxxxz,  \
                             tr_x_yyzz_xxxxyz,  \
                             tr_x_yyzz_xxxxz,   \
                             tr_x_yyzz_xxxxzz,  \
                             tr_x_yyzz_xxxyyz,  \
                             tr_x_yyzz_xxxyz,   \
                             tr_x_yyzz_xxxyzz,  \
                             tr_x_yyzz_xxxzz,   \
                             tr_x_yyzz_xxxzzz,  \
                             tr_x_yyzz_xxyyyz,  \
                             tr_x_yyzz_xxyyz,   \
                             tr_x_yyzz_xxyyzz,  \
                             tr_x_yyzz_xxyzz,   \
                             tr_x_yyzz_xxyzzz,  \
                             tr_x_yyzz_xxzzz,   \
                             tr_x_yyzz_xxzzzz,  \
                             tr_x_yyzz_xyyyyz,  \
                             tr_x_yyzz_xyyyz,   \
                             tr_x_yyzz_xyyyzz,  \
                             tr_x_yyzz_xyyzz,   \
                             tr_x_yyzz_xyyzzz,  \
                             tr_x_yyzz_xyzzz,   \
                             tr_x_yyzz_xyzzzz,  \
                             tr_x_yyzz_xzzzz,   \
                             tr_x_yyzz_xzzzzz,  \
                             tr_x_yyzz_yyyyyz,  \
                             tr_x_yyzz_yyyyz,   \
                             tr_x_yyzz_yyyyzz,  \
                             tr_x_yyzz_yyyzz,   \
                             tr_x_yyzz_yyyzzz,  \
                             tr_x_yyzz_yyzzz,   \
                             tr_x_yyzz_yyzzzz,  \
                             tr_x_yyzz_yzzzz,   \
                             tr_x_yyzz_yzzzzz,  \
                             tr_x_yyzz_zzzzz,   \
                             tr_x_yyzz_zzzzzz,  \
                             tr_x_yzz_xxxxxx,   \
                             tr_x_yzz_xxxxxz,   \
                             tr_x_yzz_xxxxyz,   \
                             tr_x_yzz_xxxxzz,   \
                             tr_x_yzz_xxxyyz,   \
                             tr_x_yzz_xxxyzz,   \
                             tr_x_yzz_xxxzzz,   \
                             tr_x_yzz_xxyyyz,   \
                             tr_x_yzz_xxyyzz,   \
                             tr_x_yzz_xxyzzz,   \
                             tr_x_yzz_xxzzzz,   \
                             tr_x_yzz_xyyyyz,   \
                             tr_x_yzz_xyyyzz,   \
                             tr_x_yzz_xyyzzz,   \
                             tr_x_yzz_xyzzzz,   \
                             tr_x_yzz_xzzzzz,   \
                             tr_x_yzz_yyyyyz,   \
                             tr_x_yzz_yyyyzz,   \
                             tr_x_yzz_yyyzzz,   \
                             tr_x_yzz_yyzzzz,   \
                             tr_x_yzz_yzzzzz,   \
                             tr_x_yzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzz_xxxxxx[i] = 2.0 * tr_x_yzz_xxxxxx[i] * fe_0 + tr_x_yyzz_xxxxxx[i] * pa_y[i];

        tr_x_yyyzz_xxxxxy[i] = tr_x_yyy_xxxxxy[i] * fe_0 + tr_x_yyyz_xxxxxy[i] * pa_z[i];

        tr_x_yyyzz_xxxxxz[i] = 2.0 * tr_x_yzz_xxxxxz[i] * fe_0 + tr_x_yyzz_xxxxxz[i] * pa_y[i];

        tr_x_yyyzz_xxxxyy[i] = tr_x_yyy_xxxxyy[i] * fe_0 + tr_x_yyyz_xxxxyy[i] * pa_z[i];

        tr_x_yyyzz_xxxxyz[i] = 2.0 * tr_x_yzz_xxxxyz[i] * fe_0 + tr_x_yyzz_xxxxz[i] * fe_0 + tr_x_yyzz_xxxxyz[i] * pa_y[i];

        tr_x_yyyzz_xxxxzz[i] = 2.0 * tr_x_yzz_xxxxzz[i] * fe_0 + tr_x_yyzz_xxxxzz[i] * pa_y[i];

        tr_x_yyyzz_xxxyyy[i] = tr_x_yyy_xxxyyy[i] * fe_0 + tr_x_yyyz_xxxyyy[i] * pa_z[i];

        tr_x_yyyzz_xxxyyz[i] = 2.0 * tr_x_yzz_xxxyyz[i] * fe_0 + 2.0 * tr_x_yyzz_xxxyz[i] * fe_0 + tr_x_yyzz_xxxyyz[i] * pa_y[i];

        tr_x_yyyzz_xxxyzz[i] = 2.0 * tr_x_yzz_xxxyzz[i] * fe_0 + tr_x_yyzz_xxxzz[i] * fe_0 + tr_x_yyzz_xxxyzz[i] * pa_y[i];

        tr_x_yyyzz_xxxzzz[i] = 2.0 * tr_x_yzz_xxxzzz[i] * fe_0 + tr_x_yyzz_xxxzzz[i] * pa_y[i];

        tr_x_yyyzz_xxyyyy[i] = tr_x_yyy_xxyyyy[i] * fe_0 + tr_x_yyyz_xxyyyy[i] * pa_z[i];

        tr_x_yyyzz_xxyyyz[i] = 2.0 * tr_x_yzz_xxyyyz[i] * fe_0 + 3.0 * tr_x_yyzz_xxyyz[i] * fe_0 + tr_x_yyzz_xxyyyz[i] * pa_y[i];

        tr_x_yyyzz_xxyyzz[i] = 2.0 * tr_x_yzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_yyzz_xxyzz[i] * fe_0 + tr_x_yyzz_xxyyzz[i] * pa_y[i];

        tr_x_yyyzz_xxyzzz[i] = 2.0 * tr_x_yzz_xxyzzz[i] * fe_0 + tr_x_yyzz_xxzzz[i] * fe_0 + tr_x_yyzz_xxyzzz[i] * pa_y[i];

        tr_x_yyyzz_xxzzzz[i] = 2.0 * tr_x_yzz_xxzzzz[i] * fe_0 + tr_x_yyzz_xxzzzz[i] * pa_y[i];

        tr_x_yyyzz_xyyyyy[i] = tr_x_yyy_xyyyyy[i] * fe_0 + tr_x_yyyz_xyyyyy[i] * pa_z[i];

        tr_x_yyyzz_xyyyyz[i] = 2.0 * tr_x_yzz_xyyyyz[i] * fe_0 + 4.0 * tr_x_yyzz_xyyyz[i] * fe_0 + tr_x_yyzz_xyyyyz[i] * pa_y[i];

        tr_x_yyyzz_xyyyzz[i] = 2.0 * tr_x_yzz_xyyyzz[i] * fe_0 + 3.0 * tr_x_yyzz_xyyzz[i] * fe_0 + tr_x_yyzz_xyyyzz[i] * pa_y[i];

        tr_x_yyyzz_xyyzzz[i] = 2.0 * tr_x_yzz_xyyzzz[i] * fe_0 + 2.0 * tr_x_yyzz_xyzzz[i] * fe_0 + tr_x_yyzz_xyyzzz[i] * pa_y[i];

        tr_x_yyyzz_xyzzzz[i] = 2.0 * tr_x_yzz_xyzzzz[i] * fe_0 + tr_x_yyzz_xzzzz[i] * fe_0 + tr_x_yyzz_xyzzzz[i] * pa_y[i];

        tr_x_yyyzz_xzzzzz[i] = 2.0 * tr_x_yzz_xzzzzz[i] * fe_0 + tr_x_yyzz_xzzzzz[i] * pa_y[i];

        tr_x_yyyzz_yyyyyy[i] = tr_x_yyy_yyyyyy[i] * fe_0 + tr_x_yyyz_yyyyyy[i] * pa_z[i];

        tr_x_yyyzz_yyyyyz[i] = 2.0 * tr_x_yzz_yyyyyz[i] * fe_0 + 5.0 * tr_x_yyzz_yyyyz[i] * fe_0 + tr_x_yyzz_yyyyyz[i] * pa_y[i];

        tr_x_yyyzz_yyyyzz[i] = 2.0 * tr_x_yzz_yyyyzz[i] * fe_0 + 4.0 * tr_x_yyzz_yyyzz[i] * fe_0 + tr_x_yyzz_yyyyzz[i] * pa_y[i];

        tr_x_yyyzz_yyyzzz[i] = 2.0 * tr_x_yzz_yyyzzz[i] * fe_0 + 3.0 * tr_x_yyzz_yyzzz[i] * fe_0 + tr_x_yyzz_yyyzzz[i] * pa_y[i];

        tr_x_yyyzz_yyzzzz[i] = 2.0 * tr_x_yzz_yyzzzz[i] * fe_0 + 2.0 * tr_x_yyzz_yzzzz[i] * fe_0 + tr_x_yyzz_yyzzzz[i] * pa_y[i];

        tr_x_yyyzz_yzzzzz[i] = 2.0 * tr_x_yzz_yzzzzz[i] * fe_0 + tr_x_yyzz_zzzzz[i] * fe_0 + tr_x_yyzz_yzzzzz[i] * pa_y[i];

        tr_x_yyyzz_zzzzzz[i] = 2.0 * tr_x_yzz_zzzzzz[i] * fe_0 + tr_x_yyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 504-532 components of targeted buffer : HI

    auto tr_x_yyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 504);

    auto tr_x_yyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 505);

    auto tr_x_yyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 506);

    auto tr_x_yyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 507);

    auto tr_x_yyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 508);

    auto tr_x_yyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 509);

    auto tr_x_yyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 510);

    auto tr_x_yyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 511);

    auto tr_x_yyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 512);

    auto tr_x_yyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 513);

    auto tr_x_yyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 514);

    auto tr_x_yyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 515);

    auto tr_x_yyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 516);

    auto tr_x_yyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 517);

    auto tr_x_yyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 518);

    auto tr_x_yyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 519);

    auto tr_x_yyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 520);

    auto tr_x_yyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 521);

    auto tr_x_yyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 522);

    auto tr_x_yyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 523);

    auto tr_x_yyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 524);

    auto tr_x_yyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 525);

    auto tr_x_yyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 526);

    auto tr_x_yyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 527);

    auto tr_x_yyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 528);

    auto tr_x_yyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 529);

    auto tr_x_yyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 530);

    auto tr_x_yyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 531);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_x_yyz_xxxxxy,   \
                             tr_x_yyz_xxxxyy,   \
                             tr_x_yyz_xxxyyy,   \
                             tr_x_yyz_xxyyyy,   \
                             tr_x_yyz_xyyyyy,   \
                             tr_x_yyz_yyyyyy,   \
                             tr_x_yyzz_xxxxxy,  \
                             tr_x_yyzz_xxxxyy,  \
                             tr_x_yyzz_xxxyyy,  \
                             tr_x_yyzz_xxyyyy,  \
                             tr_x_yyzz_xyyyyy,  \
                             tr_x_yyzz_yyyyyy,  \
                             tr_x_yyzzz_xxxxxx, \
                             tr_x_yyzzz_xxxxxy, \
                             tr_x_yyzzz_xxxxxz, \
                             tr_x_yyzzz_xxxxyy, \
                             tr_x_yyzzz_xxxxyz, \
                             tr_x_yyzzz_xxxxzz, \
                             tr_x_yyzzz_xxxyyy, \
                             tr_x_yyzzz_xxxyyz, \
                             tr_x_yyzzz_xxxyzz, \
                             tr_x_yyzzz_xxxzzz, \
                             tr_x_yyzzz_xxyyyy, \
                             tr_x_yyzzz_xxyyyz, \
                             tr_x_yyzzz_xxyyzz, \
                             tr_x_yyzzz_xxyzzz, \
                             tr_x_yyzzz_xxzzzz, \
                             tr_x_yyzzz_xyyyyy, \
                             tr_x_yyzzz_xyyyyz, \
                             tr_x_yyzzz_xyyyzz, \
                             tr_x_yyzzz_xyyzzz, \
                             tr_x_yyzzz_xyzzzz, \
                             tr_x_yyzzz_xzzzzz, \
                             tr_x_yyzzz_yyyyyy, \
                             tr_x_yyzzz_yyyyyz, \
                             tr_x_yyzzz_yyyyzz, \
                             tr_x_yyzzz_yyyzzz, \
                             tr_x_yyzzz_yyzzzz, \
                             tr_x_yyzzz_yzzzzz, \
                             tr_x_yyzzz_zzzzzz, \
                             tr_x_yzzz_xxxxxx,  \
                             tr_x_yzzz_xxxxxz,  \
                             tr_x_yzzz_xxxxyz,  \
                             tr_x_yzzz_xxxxz,   \
                             tr_x_yzzz_xxxxzz,  \
                             tr_x_yzzz_xxxyyz,  \
                             tr_x_yzzz_xxxyz,   \
                             tr_x_yzzz_xxxyzz,  \
                             tr_x_yzzz_xxxzz,   \
                             tr_x_yzzz_xxxzzz,  \
                             tr_x_yzzz_xxyyyz,  \
                             tr_x_yzzz_xxyyz,   \
                             tr_x_yzzz_xxyyzz,  \
                             tr_x_yzzz_xxyzz,   \
                             tr_x_yzzz_xxyzzz,  \
                             tr_x_yzzz_xxzzz,   \
                             tr_x_yzzz_xxzzzz,  \
                             tr_x_yzzz_xyyyyz,  \
                             tr_x_yzzz_xyyyz,   \
                             tr_x_yzzz_xyyyzz,  \
                             tr_x_yzzz_xyyzz,   \
                             tr_x_yzzz_xyyzzz,  \
                             tr_x_yzzz_xyzzz,   \
                             tr_x_yzzz_xyzzzz,  \
                             tr_x_yzzz_xzzzz,   \
                             tr_x_yzzz_xzzzzz,  \
                             tr_x_yzzz_yyyyyz,  \
                             tr_x_yzzz_yyyyz,   \
                             tr_x_yzzz_yyyyzz,  \
                             tr_x_yzzz_yyyzz,   \
                             tr_x_yzzz_yyyzzz,  \
                             tr_x_yzzz_yyzzz,   \
                             tr_x_yzzz_yyzzzz,  \
                             tr_x_yzzz_yzzzz,   \
                             tr_x_yzzz_yzzzzz,  \
                             tr_x_yzzz_zzzzz,   \
                             tr_x_yzzz_zzzzzz,  \
                             tr_x_zzz_xxxxxx,   \
                             tr_x_zzz_xxxxxz,   \
                             tr_x_zzz_xxxxyz,   \
                             tr_x_zzz_xxxxzz,   \
                             tr_x_zzz_xxxyyz,   \
                             tr_x_zzz_xxxyzz,   \
                             tr_x_zzz_xxxzzz,   \
                             tr_x_zzz_xxyyyz,   \
                             tr_x_zzz_xxyyzz,   \
                             tr_x_zzz_xxyzzz,   \
                             tr_x_zzz_xxzzzz,   \
                             tr_x_zzz_xyyyyz,   \
                             tr_x_zzz_xyyyzz,   \
                             tr_x_zzz_xyyzzz,   \
                             tr_x_zzz_xyzzzz,   \
                             tr_x_zzz_xzzzzz,   \
                             tr_x_zzz_yyyyyz,   \
                             tr_x_zzz_yyyyzz,   \
                             tr_x_zzz_yyyzzz,   \
                             tr_x_zzz_yyzzzz,   \
                             tr_x_zzz_yzzzzz,   \
                             tr_x_zzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzz_xxxxxx[i] = tr_x_zzz_xxxxxx[i] * fe_0 + tr_x_yzzz_xxxxxx[i] * pa_y[i];

        tr_x_yyzzz_xxxxxy[i] = 2.0 * tr_x_yyz_xxxxxy[i] * fe_0 + tr_x_yyzz_xxxxxy[i] * pa_z[i];

        tr_x_yyzzz_xxxxxz[i] = tr_x_zzz_xxxxxz[i] * fe_0 + tr_x_yzzz_xxxxxz[i] * pa_y[i];

        tr_x_yyzzz_xxxxyy[i] = 2.0 * tr_x_yyz_xxxxyy[i] * fe_0 + tr_x_yyzz_xxxxyy[i] * pa_z[i];

        tr_x_yyzzz_xxxxyz[i] = tr_x_zzz_xxxxyz[i] * fe_0 + tr_x_yzzz_xxxxz[i] * fe_0 + tr_x_yzzz_xxxxyz[i] * pa_y[i];

        tr_x_yyzzz_xxxxzz[i] = tr_x_zzz_xxxxzz[i] * fe_0 + tr_x_yzzz_xxxxzz[i] * pa_y[i];

        tr_x_yyzzz_xxxyyy[i] = 2.0 * tr_x_yyz_xxxyyy[i] * fe_0 + tr_x_yyzz_xxxyyy[i] * pa_z[i];

        tr_x_yyzzz_xxxyyz[i] = tr_x_zzz_xxxyyz[i] * fe_0 + 2.0 * tr_x_yzzz_xxxyz[i] * fe_0 + tr_x_yzzz_xxxyyz[i] * pa_y[i];

        tr_x_yyzzz_xxxyzz[i] = tr_x_zzz_xxxyzz[i] * fe_0 + tr_x_yzzz_xxxzz[i] * fe_0 + tr_x_yzzz_xxxyzz[i] * pa_y[i];

        tr_x_yyzzz_xxxzzz[i] = tr_x_zzz_xxxzzz[i] * fe_0 + tr_x_yzzz_xxxzzz[i] * pa_y[i];

        tr_x_yyzzz_xxyyyy[i] = 2.0 * tr_x_yyz_xxyyyy[i] * fe_0 + tr_x_yyzz_xxyyyy[i] * pa_z[i];

        tr_x_yyzzz_xxyyyz[i] = tr_x_zzz_xxyyyz[i] * fe_0 + 3.0 * tr_x_yzzz_xxyyz[i] * fe_0 + tr_x_yzzz_xxyyyz[i] * pa_y[i];

        tr_x_yyzzz_xxyyzz[i] = tr_x_zzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_yzzz_xxyzz[i] * fe_0 + tr_x_yzzz_xxyyzz[i] * pa_y[i];

        tr_x_yyzzz_xxyzzz[i] = tr_x_zzz_xxyzzz[i] * fe_0 + tr_x_yzzz_xxzzz[i] * fe_0 + tr_x_yzzz_xxyzzz[i] * pa_y[i];

        tr_x_yyzzz_xxzzzz[i] = tr_x_zzz_xxzzzz[i] * fe_0 + tr_x_yzzz_xxzzzz[i] * pa_y[i];

        tr_x_yyzzz_xyyyyy[i] = 2.0 * tr_x_yyz_xyyyyy[i] * fe_0 + tr_x_yyzz_xyyyyy[i] * pa_z[i];

        tr_x_yyzzz_xyyyyz[i] = tr_x_zzz_xyyyyz[i] * fe_0 + 4.0 * tr_x_yzzz_xyyyz[i] * fe_0 + tr_x_yzzz_xyyyyz[i] * pa_y[i];

        tr_x_yyzzz_xyyyzz[i] = tr_x_zzz_xyyyzz[i] * fe_0 + 3.0 * tr_x_yzzz_xyyzz[i] * fe_0 + tr_x_yzzz_xyyyzz[i] * pa_y[i];

        tr_x_yyzzz_xyyzzz[i] = tr_x_zzz_xyyzzz[i] * fe_0 + 2.0 * tr_x_yzzz_xyzzz[i] * fe_0 + tr_x_yzzz_xyyzzz[i] * pa_y[i];

        tr_x_yyzzz_xyzzzz[i] = tr_x_zzz_xyzzzz[i] * fe_0 + tr_x_yzzz_xzzzz[i] * fe_0 + tr_x_yzzz_xyzzzz[i] * pa_y[i];

        tr_x_yyzzz_xzzzzz[i] = tr_x_zzz_xzzzzz[i] * fe_0 + tr_x_yzzz_xzzzzz[i] * pa_y[i];

        tr_x_yyzzz_yyyyyy[i] = 2.0 * tr_x_yyz_yyyyyy[i] * fe_0 + tr_x_yyzz_yyyyyy[i] * pa_z[i];

        tr_x_yyzzz_yyyyyz[i] = tr_x_zzz_yyyyyz[i] * fe_0 + 5.0 * tr_x_yzzz_yyyyz[i] * fe_0 + tr_x_yzzz_yyyyyz[i] * pa_y[i];

        tr_x_yyzzz_yyyyzz[i] = tr_x_zzz_yyyyzz[i] * fe_0 + 4.0 * tr_x_yzzz_yyyzz[i] * fe_0 + tr_x_yzzz_yyyyzz[i] * pa_y[i];

        tr_x_yyzzz_yyyzzz[i] = tr_x_zzz_yyyzzz[i] * fe_0 + 3.0 * tr_x_yzzz_yyzzz[i] * fe_0 + tr_x_yzzz_yyyzzz[i] * pa_y[i];

        tr_x_yyzzz_yyzzzz[i] = tr_x_zzz_yyzzzz[i] * fe_0 + 2.0 * tr_x_yzzz_yzzzz[i] * fe_0 + tr_x_yzzz_yyzzzz[i] * pa_y[i];

        tr_x_yyzzz_yzzzzz[i] = tr_x_zzz_yzzzzz[i] * fe_0 + tr_x_yzzz_zzzzz[i] * fe_0 + tr_x_yzzz_yzzzzz[i] * pa_y[i];

        tr_x_yyzzz_zzzzzz[i] = tr_x_zzz_zzzzzz[i] * fe_0 + tr_x_yzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 532-560 components of targeted buffer : HI

    auto tr_x_yzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 532);

    auto tr_x_yzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 533);

    auto tr_x_yzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 534);

    auto tr_x_yzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 535);

    auto tr_x_yzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 536);

    auto tr_x_yzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 537);

    auto tr_x_yzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 538);

    auto tr_x_yzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 539);

    auto tr_x_yzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 540);

    auto tr_x_yzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 541);

    auto tr_x_yzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 542);

    auto tr_x_yzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 543);

    auto tr_x_yzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 544);

    auto tr_x_yzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 545);

    auto tr_x_yzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 546);

    auto tr_x_yzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 547);

    auto tr_x_yzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 548);

    auto tr_x_yzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 549);

    auto tr_x_yzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 550);

    auto tr_x_yzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 551);

    auto tr_x_yzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 552);

    auto tr_x_yzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 553);

    auto tr_x_yzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 554);

    auto tr_x_yzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 555);

    auto tr_x_yzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 556);

    auto tr_x_yzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 557);

    auto tr_x_yzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 558);

    auto tr_x_yzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 559);

#pragma omp simd aligned(pa_y,                  \
                             tr_x_yzzzz_xxxxxx, \
                             tr_x_yzzzz_xxxxxy, \
                             tr_x_yzzzz_xxxxxz, \
                             tr_x_yzzzz_xxxxyy, \
                             tr_x_yzzzz_xxxxyz, \
                             tr_x_yzzzz_xxxxzz, \
                             tr_x_yzzzz_xxxyyy, \
                             tr_x_yzzzz_xxxyyz, \
                             tr_x_yzzzz_xxxyzz, \
                             tr_x_yzzzz_xxxzzz, \
                             tr_x_yzzzz_xxyyyy, \
                             tr_x_yzzzz_xxyyyz, \
                             tr_x_yzzzz_xxyyzz, \
                             tr_x_yzzzz_xxyzzz, \
                             tr_x_yzzzz_xxzzzz, \
                             tr_x_yzzzz_xyyyyy, \
                             tr_x_yzzzz_xyyyyz, \
                             tr_x_yzzzz_xyyyzz, \
                             tr_x_yzzzz_xyyzzz, \
                             tr_x_yzzzz_xyzzzz, \
                             tr_x_yzzzz_xzzzzz, \
                             tr_x_yzzzz_yyyyyy, \
                             tr_x_yzzzz_yyyyyz, \
                             tr_x_yzzzz_yyyyzz, \
                             tr_x_yzzzz_yyyzzz, \
                             tr_x_yzzzz_yyzzzz, \
                             tr_x_yzzzz_yzzzzz, \
                             tr_x_yzzzz_zzzzzz, \
                             tr_x_zzzz_xxxxx,   \
                             tr_x_zzzz_xxxxxx,  \
                             tr_x_zzzz_xxxxxy,  \
                             tr_x_zzzz_xxxxxz,  \
                             tr_x_zzzz_xxxxy,   \
                             tr_x_zzzz_xxxxyy,  \
                             tr_x_zzzz_xxxxyz,  \
                             tr_x_zzzz_xxxxz,   \
                             tr_x_zzzz_xxxxzz,  \
                             tr_x_zzzz_xxxyy,   \
                             tr_x_zzzz_xxxyyy,  \
                             tr_x_zzzz_xxxyyz,  \
                             tr_x_zzzz_xxxyz,   \
                             tr_x_zzzz_xxxyzz,  \
                             tr_x_zzzz_xxxzz,   \
                             tr_x_zzzz_xxxzzz,  \
                             tr_x_zzzz_xxyyy,   \
                             tr_x_zzzz_xxyyyy,  \
                             tr_x_zzzz_xxyyyz,  \
                             tr_x_zzzz_xxyyz,   \
                             tr_x_zzzz_xxyyzz,  \
                             tr_x_zzzz_xxyzz,   \
                             tr_x_zzzz_xxyzzz,  \
                             tr_x_zzzz_xxzzz,   \
                             tr_x_zzzz_xxzzzz,  \
                             tr_x_zzzz_xyyyy,   \
                             tr_x_zzzz_xyyyyy,  \
                             tr_x_zzzz_xyyyyz,  \
                             tr_x_zzzz_xyyyz,   \
                             tr_x_zzzz_xyyyzz,  \
                             tr_x_zzzz_xyyzz,   \
                             tr_x_zzzz_xyyzzz,  \
                             tr_x_zzzz_xyzzz,   \
                             tr_x_zzzz_xyzzzz,  \
                             tr_x_zzzz_xzzzz,   \
                             tr_x_zzzz_xzzzzz,  \
                             tr_x_zzzz_yyyyy,   \
                             tr_x_zzzz_yyyyyy,  \
                             tr_x_zzzz_yyyyyz,  \
                             tr_x_zzzz_yyyyz,   \
                             tr_x_zzzz_yyyyzz,  \
                             tr_x_zzzz_yyyzz,   \
                             tr_x_zzzz_yyyzzz,  \
                             tr_x_zzzz_yyzzz,   \
                             tr_x_zzzz_yyzzzz,  \
                             tr_x_zzzz_yzzzz,   \
                             tr_x_zzzz_yzzzzz,  \
                             tr_x_zzzz_zzzzz,   \
                             tr_x_zzzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzz_xxxxxx[i] = tr_x_zzzz_xxxxxx[i] * pa_y[i];

        tr_x_yzzzz_xxxxxy[i] = tr_x_zzzz_xxxxx[i] * fe_0 + tr_x_zzzz_xxxxxy[i] * pa_y[i];

        tr_x_yzzzz_xxxxxz[i] = tr_x_zzzz_xxxxxz[i] * pa_y[i];

        tr_x_yzzzz_xxxxyy[i] = 2.0 * tr_x_zzzz_xxxxy[i] * fe_0 + tr_x_zzzz_xxxxyy[i] * pa_y[i];

        tr_x_yzzzz_xxxxyz[i] = tr_x_zzzz_xxxxz[i] * fe_0 + tr_x_zzzz_xxxxyz[i] * pa_y[i];

        tr_x_yzzzz_xxxxzz[i] = tr_x_zzzz_xxxxzz[i] * pa_y[i];

        tr_x_yzzzz_xxxyyy[i] = 3.0 * tr_x_zzzz_xxxyy[i] * fe_0 + tr_x_zzzz_xxxyyy[i] * pa_y[i];

        tr_x_yzzzz_xxxyyz[i] = 2.0 * tr_x_zzzz_xxxyz[i] * fe_0 + tr_x_zzzz_xxxyyz[i] * pa_y[i];

        tr_x_yzzzz_xxxyzz[i] = tr_x_zzzz_xxxzz[i] * fe_0 + tr_x_zzzz_xxxyzz[i] * pa_y[i];

        tr_x_yzzzz_xxxzzz[i] = tr_x_zzzz_xxxzzz[i] * pa_y[i];

        tr_x_yzzzz_xxyyyy[i] = 4.0 * tr_x_zzzz_xxyyy[i] * fe_0 + tr_x_zzzz_xxyyyy[i] * pa_y[i];

        tr_x_yzzzz_xxyyyz[i] = 3.0 * tr_x_zzzz_xxyyz[i] * fe_0 + tr_x_zzzz_xxyyyz[i] * pa_y[i];

        tr_x_yzzzz_xxyyzz[i] = 2.0 * tr_x_zzzz_xxyzz[i] * fe_0 + tr_x_zzzz_xxyyzz[i] * pa_y[i];

        tr_x_yzzzz_xxyzzz[i] = tr_x_zzzz_xxzzz[i] * fe_0 + tr_x_zzzz_xxyzzz[i] * pa_y[i];

        tr_x_yzzzz_xxzzzz[i] = tr_x_zzzz_xxzzzz[i] * pa_y[i];

        tr_x_yzzzz_xyyyyy[i] = 5.0 * tr_x_zzzz_xyyyy[i] * fe_0 + tr_x_zzzz_xyyyyy[i] * pa_y[i];

        tr_x_yzzzz_xyyyyz[i] = 4.0 * tr_x_zzzz_xyyyz[i] * fe_0 + tr_x_zzzz_xyyyyz[i] * pa_y[i];

        tr_x_yzzzz_xyyyzz[i] = 3.0 * tr_x_zzzz_xyyzz[i] * fe_0 + tr_x_zzzz_xyyyzz[i] * pa_y[i];

        tr_x_yzzzz_xyyzzz[i] = 2.0 * tr_x_zzzz_xyzzz[i] * fe_0 + tr_x_zzzz_xyyzzz[i] * pa_y[i];

        tr_x_yzzzz_xyzzzz[i] = tr_x_zzzz_xzzzz[i] * fe_0 + tr_x_zzzz_xyzzzz[i] * pa_y[i];

        tr_x_yzzzz_xzzzzz[i] = tr_x_zzzz_xzzzzz[i] * pa_y[i];

        tr_x_yzzzz_yyyyyy[i] = 6.0 * tr_x_zzzz_yyyyy[i] * fe_0 + tr_x_zzzz_yyyyyy[i] * pa_y[i];

        tr_x_yzzzz_yyyyyz[i] = 5.0 * tr_x_zzzz_yyyyz[i] * fe_0 + tr_x_zzzz_yyyyyz[i] * pa_y[i];

        tr_x_yzzzz_yyyyzz[i] = 4.0 * tr_x_zzzz_yyyzz[i] * fe_0 + tr_x_zzzz_yyyyzz[i] * pa_y[i];

        tr_x_yzzzz_yyyzzz[i] = 3.0 * tr_x_zzzz_yyzzz[i] * fe_0 + tr_x_zzzz_yyyzzz[i] * pa_y[i];

        tr_x_yzzzz_yyzzzz[i] = 2.0 * tr_x_zzzz_yzzzz[i] * fe_0 + tr_x_zzzz_yyzzzz[i] * pa_y[i];

        tr_x_yzzzz_yzzzzz[i] = tr_x_zzzz_zzzzz[i] * fe_0 + tr_x_zzzz_yzzzzz[i] * pa_y[i];

        tr_x_yzzzz_zzzzzz[i] = tr_x_zzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 560-588 components of targeted buffer : HI

    auto tr_x_zzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 560);

    auto tr_x_zzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 561);

    auto tr_x_zzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 562);

    auto tr_x_zzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 563);

    auto tr_x_zzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 564);

    auto tr_x_zzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 565);

    auto tr_x_zzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 566);

    auto tr_x_zzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 567);

    auto tr_x_zzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 568);

    auto tr_x_zzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 569);

    auto tr_x_zzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 570);

    auto tr_x_zzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 571);

    auto tr_x_zzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 572);

    auto tr_x_zzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 573);

    auto tr_x_zzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 574);

    auto tr_x_zzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 575);

    auto tr_x_zzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 576);

    auto tr_x_zzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 577);

    auto tr_x_zzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 578);

    auto tr_x_zzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 579);

    auto tr_x_zzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 580);

    auto tr_x_zzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 581);

    auto tr_x_zzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 582);

    auto tr_x_zzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 583);

    auto tr_x_zzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 584);

    auto tr_x_zzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 585);

    auto tr_x_zzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 586);

    auto tr_x_zzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 587);

#pragma omp simd aligned(pa_z,                  \
                             tr_x_zzz_xxxxxx,   \
                             tr_x_zzz_xxxxxy,   \
                             tr_x_zzz_xxxxxz,   \
                             tr_x_zzz_xxxxyy,   \
                             tr_x_zzz_xxxxyz,   \
                             tr_x_zzz_xxxxzz,   \
                             tr_x_zzz_xxxyyy,   \
                             tr_x_zzz_xxxyyz,   \
                             tr_x_zzz_xxxyzz,   \
                             tr_x_zzz_xxxzzz,   \
                             tr_x_zzz_xxyyyy,   \
                             tr_x_zzz_xxyyyz,   \
                             tr_x_zzz_xxyyzz,   \
                             tr_x_zzz_xxyzzz,   \
                             tr_x_zzz_xxzzzz,   \
                             tr_x_zzz_xyyyyy,   \
                             tr_x_zzz_xyyyyz,   \
                             tr_x_zzz_xyyyzz,   \
                             tr_x_zzz_xyyzzz,   \
                             tr_x_zzz_xyzzzz,   \
                             tr_x_zzz_xzzzzz,   \
                             tr_x_zzz_yyyyyy,   \
                             tr_x_zzz_yyyyyz,   \
                             tr_x_zzz_yyyyzz,   \
                             tr_x_zzz_yyyzzz,   \
                             tr_x_zzz_yyzzzz,   \
                             tr_x_zzz_yzzzzz,   \
                             tr_x_zzz_zzzzzz,   \
                             tr_x_zzzz_xxxxx,   \
                             tr_x_zzzz_xxxxxx,  \
                             tr_x_zzzz_xxxxxy,  \
                             tr_x_zzzz_xxxxxz,  \
                             tr_x_zzzz_xxxxy,   \
                             tr_x_zzzz_xxxxyy,  \
                             tr_x_zzzz_xxxxyz,  \
                             tr_x_zzzz_xxxxz,   \
                             tr_x_zzzz_xxxxzz,  \
                             tr_x_zzzz_xxxyy,   \
                             tr_x_zzzz_xxxyyy,  \
                             tr_x_zzzz_xxxyyz,  \
                             tr_x_zzzz_xxxyz,   \
                             tr_x_zzzz_xxxyzz,  \
                             tr_x_zzzz_xxxzz,   \
                             tr_x_zzzz_xxxzzz,  \
                             tr_x_zzzz_xxyyy,   \
                             tr_x_zzzz_xxyyyy,  \
                             tr_x_zzzz_xxyyyz,  \
                             tr_x_zzzz_xxyyz,   \
                             tr_x_zzzz_xxyyzz,  \
                             tr_x_zzzz_xxyzz,   \
                             tr_x_zzzz_xxyzzz,  \
                             tr_x_zzzz_xxzzz,   \
                             tr_x_zzzz_xxzzzz,  \
                             tr_x_zzzz_xyyyy,   \
                             tr_x_zzzz_xyyyyy,  \
                             tr_x_zzzz_xyyyyz,  \
                             tr_x_zzzz_xyyyz,   \
                             tr_x_zzzz_xyyyzz,  \
                             tr_x_zzzz_xyyzz,   \
                             tr_x_zzzz_xyyzzz,  \
                             tr_x_zzzz_xyzzz,   \
                             tr_x_zzzz_xyzzzz,  \
                             tr_x_zzzz_xzzzz,   \
                             tr_x_zzzz_xzzzzz,  \
                             tr_x_zzzz_yyyyy,   \
                             tr_x_zzzz_yyyyyy,  \
                             tr_x_zzzz_yyyyyz,  \
                             tr_x_zzzz_yyyyz,   \
                             tr_x_zzzz_yyyyzz,  \
                             tr_x_zzzz_yyyzz,   \
                             tr_x_zzzz_yyyzzz,  \
                             tr_x_zzzz_yyzzz,   \
                             tr_x_zzzz_yyzzzz,  \
                             tr_x_zzzz_yzzzz,   \
                             tr_x_zzzz_yzzzzz,  \
                             tr_x_zzzz_zzzzz,   \
                             tr_x_zzzz_zzzzzz,  \
                             tr_x_zzzzz_xxxxxx, \
                             tr_x_zzzzz_xxxxxy, \
                             tr_x_zzzzz_xxxxxz, \
                             tr_x_zzzzz_xxxxyy, \
                             tr_x_zzzzz_xxxxyz, \
                             tr_x_zzzzz_xxxxzz, \
                             tr_x_zzzzz_xxxyyy, \
                             tr_x_zzzzz_xxxyyz, \
                             tr_x_zzzzz_xxxyzz, \
                             tr_x_zzzzz_xxxzzz, \
                             tr_x_zzzzz_xxyyyy, \
                             tr_x_zzzzz_xxyyyz, \
                             tr_x_zzzzz_xxyyzz, \
                             tr_x_zzzzz_xxyzzz, \
                             tr_x_zzzzz_xxzzzz, \
                             tr_x_zzzzz_xyyyyy, \
                             tr_x_zzzzz_xyyyyz, \
                             tr_x_zzzzz_xyyyzz, \
                             tr_x_zzzzz_xyyzzz, \
                             tr_x_zzzzz_xyzzzz, \
                             tr_x_zzzzz_xzzzzz, \
                             tr_x_zzzzz_yyyyyy, \
                             tr_x_zzzzz_yyyyyz, \
                             tr_x_zzzzz_yyyyzz, \
                             tr_x_zzzzz_yyyzzz, \
                             tr_x_zzzzz_yyzzzz, \
                             tr_x_zzzzz_yzzzzz, \
                             tr_x_zzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzz_xxxxxx[i] = 4.0 * tr_x_zzz_xxxxxx[i] * fe_0 + tr_x_zzzz_xxxxxx[i] * pa_z[i];

        tr_x_zzzzz_xxxxxy[i] = 4.0 * tr_x_zzz_xxxxxy[i] * fe_0 + tr_x_zzzz_xxxxxy[i] * pa_z[i];

        tr_x_zzzzz_xxxxxz[i] = 4.0 * tr_x_zzz_xxxxxz[i] * fe_0 + tr_x_zzzz_xxxxx[i] * fe_0 + tr_x_zzzz_xxxxxz[i] * pa_z[i];

        tr_x_zzzzz_xxxxyy[i] = 4.0 * tr_x_zzz_xxxxyy[i] * fe_0 + tr_x_zzzz_xxxxyy[i] * pa_z[i];

        tr_x_zzzzz_xxxxyz[i] = 4.0 * tr_x_zzz_xxxxyz[i] * fe_0 + tr_x_zzzz_xxxxy[i] * fe_0 + tr_x_zzzz_xxxxyz[i] * pa_z[i];

        tr_x_zzzzz_xxxxzz[i] = 4.0 * tr_x_zzz_xxxxzz[i] * fe_0 + 2.0 * tr_x_zzzz_xxxxz[i] * fe_0 + tr_x_zzzz_xxxxzz[i] * pa_z[i];

        tr_x_zzzzz_xxxyyy[i] = 4.0 * tr_x_zzz_xxxyyy[i] * fe_0 + tr_x_zzzz_xxxyyy[i] * pa_z[i];

        tr_x_zzzzz_xxxyyz[i] = 4.0 * tr_x_zzz_xxxyyz[i] * fe_0 + tr_x_zzzz_xxxyy[i] * fe_0 + tr_x_zzzz_xxxyyz[i] * pa_z[i];

        tr_x_zzzzz_xxxyzz[i] = 4.0 * tr_x_zzz_xxxyzz[i] * fe_0 + 2.0 * tr_x_zzzz_xxxyz[i] * fe_0 + tr_x_zzzz_xxxyzz[i] * pa_z[i];

        tr_x_zzzzz_xxxzzz[i] = 4.0 * tr_x_zzz_xxxzzz[i] * fe_0 + 3.0 * tr_x_zzzz_xxxzz[i] * fe_0 + tr_x_zzzz_xxxzzz[i] * pa_z[i];

        tr_x_zzzzz_xxyyyy[i] = 4.0 * tr_x_zzz_xxyyyy[i] * fe_0 + tr_x_zzzz_xxyyyy[i] * pa_z[i];

        tr_x_zzzzz_xxyyyz[i] = 4.0 * tr_x_zzz_xxyyyz[i] * fe_0 + tr_x_zzzz_xxyyy[i] * fe_0 + tr_x_zzzz_xxyyyz[i] * pa_z[i];

        tr_x_zzzzz_xxyyzz[i] = 4.0 * tr_x_zzz_xxyyzz[i] * fe_0 + 2.0 * tr_x_zzzz_xxyyz[i] * fe_0 + tr_x_zzzz_xxyyzz[i] * pa_z[i];

        tr_x_zzzzz_xxyzzz[i] = 4.0 * tr_x_zzz_xxyzzz[i] * fe_0 + 3.0 * tr_x_zzzz_xxyzz[i] * fe_0 + tr_x_zzzz_xxyzzz[i] * pa_z[i];

        tr_x_zzzzz_xxzzzz[i] = 4.0 * tr_x_zzz_xxzzzz[i] * fe_0 + 4.0 * tr_x_zzzz_xxzzz[i] * fe_0 + tr_x_zzzz_xxzzzz[i] * pa_z[i];

        tr_x_zzzzz_xyyyyy[i] = 4.0 * tr_x_zzz_xyyyyy[i] * fe_0 + tr_x_zzzz_xyyyyy[i] * pa_z[i];

        tr_x_zzzzz_xyyyyz[i] = 4.0 * tr_x_zzz_xyyyyz[i] * fe_0 + tr_x_zzzz_xyyyy[i] * fe_0 + tr_x_zzzz_xyyyyz[i] * pa_z[i];

        tr_x_zzzzz_xyyyzz[i] = 4.0 * tr_x_zzz_xyyyzz[i] * fe_0 + 2.0 * tr_x_zzzz_xyyyz[i] * fe_0 + tr_x_zzzz_xyyyzz[i] * pa_z[i];

        tr_x_zzzzz_xyyzzz[i] = 4.0 * tr_x_zzz_xyyzzz[i] * fe_0 + 3.0 * tr_x_zzzz_xyyzz[i] * fe_0 + tr_x_zzzz_xyyzzz[i] * pa_z[i];

        tr_x_zzzzz_xyzzzz[i] = 4.0 * tr_x_zzz_xyzzzz[i] * fe_0 + 4.0 * tr_x_zzzz_xyzzz[i] * fe_0 + tr_x_zzzz_xyzzzz[i] * pa_z[i];

        tr_x_zzzzz_xzzzzz[i] = 4.0 * tr_x_zzz_xzzzzz[i] * fe_0 + 5.0 * tr_x_zzzz_xzzzz[i] * fe_0 + tr_x_zzzz_xzzzzz[i] * pa_z[i];

        tr_x_zzzzz_yyyyyy[i] = 4.0 * tr_x_zzz_yyyyyy[i] * fe_0 + tr_x_zzzz_yyyyyy[i] * pa_z[i];

        tr_x_zzzzz_yyyyyz[i] = 4.0 * tr_x_zzz_yyyyyz[i] * fe_0 + tr_x_zzzz_yyyyy[i] * fe_0 + tr_x_zzzz_yyyyyz[i] * pa_z[i];

        tr_x_zzzzz_yyyyzz[i] = 4.0 * tr_x_zzz_yyyyzz[i] * fe_0 + 2.0 * tr_x_zzzz_yyyyz[i] * fe_0 + tr_x_zzzz_yyyyzz[i] * pa_z[i];

        tr_x_zzzzz_yyyzzz[i] = 4.0 * tr_x_zzz_yyyzzz[i] * fe_0 + 3.0 * tr_x_zzzz_yyyzz[i] * fe_0 + tr_x_zzzz_yyyzzz[i] * pa_z[i];

        tr_x_zzzzz_yyzzzz[i] = 4.0 * tr_x_zzz_yyzzzz[i] * fe_0 + 4.0 * tr_x_zzzz_yyzzz[i] * fe_0 + tr_x_zzzz_yyzzzz[i] * pa_z[i];

        tr_x_zzzzz_yzzzzz[i] = 4.0 * tr_x_zzz_yzzzzz[i] * fe_0 + 5.0 * tr_x_zzzz_yzzzz[i] * fe_0 + tr_x_zzzz_yzzzzz[i] * pa_z[i];

        tr_x_zzzzz_zzzzzz[i] = 4.0 * tr_x_zzz_zzzzzz[i] * fe_0 + 6.0 * tr_x_zzzz_zzzzz[i] * fe_0 + tr_x_zzzz_zzzzzz[i] * pa_z[i];
    }

    // Set up 588-616 components of targeted buffer : HI

    auto tr_y_xxxxx_xxxxxx = pbuffer.data(idx_dip_hi + 588);

    auto tr_y_xxxxx_xxxxxy = pbuffer.data(idx_dip_hi + 589);

    auto tr_y_xxxxx_xxxxxz = pbuffer.data(idx_dip_hi + 590);

    auto tr_y_xxxxx_xxxxyy = pbuffer.data(idx_dip_hi + 591);

    auto tr_y_xxxxx_xxxxyz = pbuffer.data(idx_dip_hi + 592);

    auto tr_y_xxxxx_xxxxzz = pbuffer.data(idx_dip_hi + 593);

    auto tr_y_xxxxx_xxxyyy = pbuffer.data(idx_dip_hi + 594);

    auto tr_y_xxxxx_xxxyyz = pbuffer.data(idx_dip_hi + 595);

    auto tr_y_xxxxx_xxxyzz = pbuffer.data(idx_dip_hi + 596);

    auto tr_y_xxxxx_xxxzzz = pbuffer.data(idx_dip_hi + 597);

    auto tr_y_xxxxx_xxyyyy = pbuffer.data(idx_dip_hi + 598);

    auto tr_y_xxxxx_xxyyyz = pbuffer.data(idx_dip_hi + 599);

    auto tr_y_xxxxx_xxyyzz = pbuffer.data(idx_dip_hi + 600);

    auto tr_y_xxxxx_xxyzzz = pbuffer.data(idx_dip_hi + 601);

    auto tr_y_xxxxx_xxzzzz = pbuffer.data(idx_dip_hi + 602);

    auto tr_y_xxxxx_xyyyyy = pbuffer.data(idx_dip_hi + 603);

    auto tr_y_xxxxx_xyyyyz = pbuffer.data(idx_dip_hi + 604);

    auto tr_y_xxxxx_xyyyzz = pbuffer.data(idx_dip_hi + 605);

    auto tr_y_xxxxx_xyyzzz = pbuffer.data(idx_dip_hi + 606);

    auto tr_y_xxxxx_xyzzzz = pbuffer.data(idx_dip_hi + 607);

    auto tr_y_xxxxx_xzzzzz = pbuffer.data(idx_dip_hi + 608);

    auto tr_y_xxxxx_yyyyyy = pbuffer.data(idx_dip_hi + 609);

    auto tr_y_xxxxx_yyyyyz = pbuffer.data(idx_dip_hi + 610);

    auto tr_y_xxxxx_yyyyzz = pbuffer.data(idx_dip_hi + 611);

    auto tr_y_xxxxx_yyyzzz = pbuffer.data(idx_dip_hi + 612);

    auto tr_y_xxxxx_yyzzzz = pbuffer.data(idx_dip_hi + 613);

    auto tr_y_xxxxx_yzzzzz = pbuffer.data(idx_dip_hi + 614);

    auto tr_y_xxxxx_zzzzzz = pbuffer.data(idx_dip_hi + 615);

#pragma omp simd aligned(pa_x,                  \
                             tr_y_xxx_xxxxxx,   \
                             tr_y_xxx_xxxxxy,   \
                             tr_y_xxx_xxxxxz,   \
                             tr_y_xxx_xxxxyy,   \
                             tr_y_xxx_xxxxyz,   \
                             tr_y_xxx_xxxxzz,   \
                             tr_y_xxx_xxxyyy,   \
                             tr_y_xxx_xxxyyz,   \
                             tr_y_xxx_xxxyzz,   \
                             tr_y_xxx_xxxzzz,   \
                             tr_y_xxx_xxyyyy,   \
                             tr_y_xxx_xxyyyz,   \
                             tr_y_xxx_xxyyzz,   \
                             tr_y_xxx_xxyzzz,   \
                             tr_y_xxx_xxzzzz,   \
                             tr_y_xxx_xyyyyy,   \
                             tr_y_xxx_xyyyyz,   \
                             tr_y_xxx_xyyyzz,   \
                             tr_y_xxx_xyyzzz,   \
                             tr_y_xxx_xyzzzz,   \
                             tr_y_xxx_xzzzzz,   \
                             tr_y_xxx_yyyyyy,   \
                             tr_y_xxx_yyyyyz,   \
                             tr_y_xxx_yyyyzz,   \
                             tr_y_xxx_yyyzzz,   \
                             tr_y_xxx_yyzzzz,   \
                             tr_y_xxx_yzzzzz,   \
                             tr_y_xxx_zzzzzz,   \
                             tr_y_xxxx_xxxxx,   \
                             tr_y_xxxx_xxxxxx,  \
                             tr_y_xxxx_xxxxxy,  \
                             tr_y_xxxx_xxxxxz,  \
                             tr_y_xxxx_xxxxy,   \
                             tr_y_xxxx_xxxxyy,  \
                             tr_y_xxxx_xxxxyz,  \
                             tr_y_xxxx_xxxxz,   \
                             tr_y_xxxx_xxxxzz,  \
                             tr_y_xxxx_xxxyy,   \
                             tr_y_xxxx_xxxyyy,  \
                             tr_y_xxxx_xxxyyz,  \
                             tr_y_xxxx_xxxyz,   \
                             tr_y_xxxx_xxxyzz,  \
                             tr_y_xxxx_xxxzz,   \
                             tr_y_xxxx_xxxzzz,  \
                             tr_y_xxxx_xxyyy,   \
                             tr_y_xxxx_xxyyyy,  \
                             tr_y_xxxx_xxyyyz,  \
                             tr_y_xxxx_xxyyz,   \
                             tr_y_xxxx_xxyyzz,  \
                             tr_y_xxxx_xxyzz,   \
                             tr_y_xxxx_xxyzzz,  \
                             tr_y_xxxx_xxzzz,   \
                             tr_y_xxxx_xxzzzz,  \
                             tr_y_xxxx_xyyyy,   \
                             tr_y_xxxx_xyyyyy,  \
                             tr_y_xxxx_xyyyyz,  \
                             tr_y_xxxx_xyyyz,   \
                             tr_y_xxxx_xyyyzz,  \
                             tr_y_xxxx_xyyzz,   \
                             tr_y_xxxx_xyyzzz,  \
                             tr_y_xxxx_xyzzz,   \
                             tr_y_xxxx_xyzzzz,  \
                             tr_y_xxxx_xzzzz,   \
                             tr_y_xxxx_xzzzzz,  \
                             tr_y_xxxx_yyyyy,   \
                             tr_y_xxxx_yyyyyy,  \
                             tr_y_xxxx_yyyyyz,  \
                             tr_y_xxxx_yyyyz,   \
                             tr_y_xxxx_yyyyzz,  \
                             tr_y_xxxx_yyyzz,   \
                             tr_y_xxxx_yyyzzz,  \
                             tr_y_xxxx_yyzzz,   \
                             tr_y_xxxx_yyzzzz,  \
                             tr_y_xxxx_yzzzz,   \
                             tr_y_xxxx_yzzzzz,  \
                             tr_y_xxxx_zzzzz,   \
                             tr_y_xxxx_zzzzzz,  \
                             tr_y_xxxxx_xxxxxx, \
                             tr_y_xxxxx_xxxxxy, \
                             tr_y_xxxxx_xxxxxz, \
                             tr_y_xxxxx_xxxxyy, \
                             tr_y_xxxxx_xxxxyz, \
                             tr_y_xxxxx_xxxxzz, \
                             tr_y_xxxxx_xxxyyy, \
                             tr_y_xxxxx_xxxyyz, \
                             tr_y_xxxxx_xxxyzz, \
                             tr_y_xxxxx_xxxzzz, \
                             tr_y_xxxxx_xxyyyy, \
                             tr_y_xxxxx_xxyyyz, \
                             tr_y_xxxxx_xxyyzz, \
                             tr_y_xxxxx_xxyzzz, \
                             tr_y_xxxxx_xxzzzz, \
                             tr_y_xxxxx_xyyyyy, \
                             tr_y_xxxxx_xyyyyz, \
                             tr_y_xxxxx_xyyyzz, \
                             tr_y_xxxxx_xyyzzz, \
                             tr_y_xxxxx_xyzzzz, \
                             tr_y_xxxxx_xzzzzz, \
                             tr_y_xxxxx_yyyyyy, \
                             tr_y_xxxxx_yyyyyz, \
                             tr_y_xxxxx_yyyyzz, \
                             tr_y_xxxxx_yyyzzz, \
                             tr_y_xxxxx_yyzzzz, \
                             tr_y_xxxxx_yzzzzz, \
                             tr_y_xxxxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxx_xxxxxx[i] = 4.0 * tr_y_xxx_xxxxxx[i] * fe_0 + 6.0 * tr_y_xxxx_xxxxx[i] * fe_0 + tr_y_xxxx_xxxxxx[i] * pa_x[i];

        tr_y_xxxxx_xxxxxy[i] = 4.0 * tr_y_xxx_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxxx_xxxxy[i] * fe_0 + tr_y_xxxx_xxxxxy[i] * pa_x[i];

        tr_y_xxxxx_xxxxxz[i] = 4.0 * tr_y_xxx_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxxx_xxxxz[i] * fe_0 + tr_y_xxxx_xxxxxz[i] * pa_x[i];

        tr_y_xxxxx_xxxxyy[i] = 4.0 * tr_y_xxx_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxxx_xxxyy[i] * fe_0 + tr_y_xxxx_xxxxyy[i] * pa_x[i];

        tr_y_xxxxx_xxxxyz[i] = 4.0 * tr_y_xxx_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxxx_xxxyz[i] * fe_0 + tr_y_xxxx_xxxxyz[i] * pa_x[i];

        tr_y_xxxxx_xxxxzz[i] = 4.0 * tr_y_xxx_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxxx_xxxzz[i] * fe_0 + tr_y_xxxx_xxxxzz[i] * pa_x[i];

        tr_y_xxxxx_xxxyyy[i] = 4.0 * tr_y_xxx_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxxx_xxyyy[i] * fe_0 + tr_y_xxxx_xxxyyy[i] * pa_x[i];

        tr_y_xxxxx_xxxyyz[i] = 4.0 * tr_y_xxx_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxxx_xxyyz[i] * fe_0 + tr_y_xxxx_xxxyyz[i] * pa_x[i];

        tr_y_xxxxx_xxxyzz[i] = 4.0 * tr_y_xxx_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxxx_xxyzz[i] * fe_0 + tr_y_xxxx_xxxyzz[i] * pa_x[i];

        tr_y_xxxxx_xxxzzz[i] = 4.0 * tr_y_xxx_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxxx_xxzzz[i] * fe_0 + tr_y_xxxx_xxxzzz[i] * pa_x[i];

        tr_y_xxxxx_xxyyyy[i] = 4.0 * tr_y_xxx_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxxx_xyyyy[i] * fe_0 + tr_y_xxxx_xxyyyy[i] * pa_x[i];

        tr_y_xxxxx_xxyyyz[i] = 4.0 * tr_y_xxx_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxxx_xyyyz[i] * fe_0 + tr_y_xxxx_xxyyyz[i] * pa_x[i];

        tr_y_xxxxx_xxyyzz[i] = 4.0 * tr_y_xxx_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxxx_xyyzz[i] * fe_0 + tr_y_xxxx_xxyyzz[i] * pa_x[i];

        tr_y_xxxxx_xxyzzz[i] = 4.0 * tr_y_xxx_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxxx_xyzzz[i] * fe_0 + tr_y_xxxx_xxyzzz[i] * pa_x[i];

        tr_y_xxxxx_xxzzzz[i] = 4.0 * tr_y_xxx_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxxx_xzzzz[i] * fe_0 + tr_y_xxxx_xxzzzz[i] * pa_x[i];

        tr_y_xxxxx_xyyyyy[i] = 4.0 * tr_y_xxx_xyyyyy[i] * fe_0 + tr_y_xxxx_yyyyy[i] * fe_0 + tr_y_xxxx_xyyyyy[i] * pa_x[i];

        tr_y_xxxxx_xyyyyz[i] = 4.0 * tr_y_xxx_xyyyyz[i] * fe_0 + tr_y_xxxx_yyyyz[i] * fe_0 + tr_y_xxxx_xyyyyz[i] * pa_x[i];

        tr_y_xxxxx_xyyyzz[i] = 4.0 * tr_y_xxx_xyyyzz[i] * fe_0 + tr_y_xxxx_yyyzz[i] * fe_0 + tr_y_xxxx_xyyyzz[i] * pa_x[i];

        tr_y_xxxxx_xyyzzz[i] = 4.0 * tr_y_xxx_xyyzzz[i] * fe_0 + tr_y_xxxx_yyzzz[i] * fe_0 + tr_y_xxxx_xyyzzz[i] * pa_x[i];

        tr_y_xxxxx_xyzzzz[i] = 4.0 * tr_y_xxx_xyzzzz[i] * fe_0 + tr_y_xxxx_yzzzz[i] * fe_0 + tr_y_xxxx_xyzzzz[i] * pa_x[i];

        tr_y_xxxxx_xzzzzz[i] = 4.0 * tr_y_xxx_xzzzzz[i] * fe_0 + tr_y_xxxx_zzzzz[i] * fe_0 + tr_y_xxxx_xzzzzz[i] * pa_x[i];

        tr_y_xxxxx_yyyyyy[i] = 4.0 * tr_y_xxx_yyyyyy[i] * fe_0 + tr_y_xxxx_yyyyyy[i] * pa_x[i];

        tr_y_xxxxx_yyyyyz[i] = 4.0 * tr_y_xxx_yyyyyz[i] * fe_0 + tr_y_xxxx_yyyyyz[i] * pa_x[i];

        tr_y_xxxxx_yyyyzz[i] = 4.0 * tr_y_xxx_yyyyzz[i] * fe_0 + tr_y_xxxx_yyyyzz[i] * pa_x[i];

        tr_y_xxxxx_yyyzzz[i] = 4.0 * tr_y_xxx_yyyzzz[i] * fe_0 + tr_y_xxxx_yyyzzz[i] * pa_x[i];

        tr_y_xxxxx_yyzzzz[i] = 4.0 * tr_y_xxx_yyzzzz[i] * fe_0 + tr_y_xxxx_yyzzzz[i] * pa_x[i];

        tr_y_xxxxx_yzzzzz[i] = 4.0 * tr_y_xxx_yzzzzz[i] * fe_0 + tr_y_xxxx_yzzzzz[i] * pa_x[i];

        tr_y_xxxxx_zzzzzz[i] = 4.0 * tr_y_xxx_zzzzzz[i] * fe_0 + tr_y_xxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 616-644 components of targeted buffer : HI

    auto tr_y_xxxxy_xxxxxx = pbuffer.data(idx_dip_hi + 616);

    auto tr_y_xxxxy_xxxxxy = pbuffer.data(idx_dip_hi + 617);

    auto tr_y_xxxxy_xxxxxz = pbuffer.data(idx_dip_hi + 618);

    auto tr_y_xxxxy_xxxxyy = pbuffer.data(idx_dip_hi + 619);

    auto tr_y_xxxxy_xxxxyz = pbuffer.data(idx_dip_hi + 620);

    auto tr_y_xxxxy_xxxxzz = pbuffer.data(idx_dip_hi + 621);

    auto tr_y_xxxxy_xxxyyy = pbuffer.data(idx_dip_hi + 622);

    auto tr_y_xxxxy_xxxyyz = pbuffer.data(idx_dip_hi + 623);

    auto tr_y_xxxxy_xxxyzz = pbuffer.data(idx_dip_hi + 624);

    auto tr_y_xxxxy_xxxzzz = pbuffer.data(idx_dip_hi + 625);

    auto tr_y_xxxxy_xxyyyy = pbuffer.data(idx_dip_hi + 626);

    auto tr_y_xxxxy_xxyyyz = pbuffer.data(idx_dip_hi + 627);

    auto tr_y_xxxxy_xxyyzz = pbuffer.data(idx_dip_hi + 628);

    auto tr_y_xxxxy_xxyzzz = pbuffer.data(idx_dip_hi + 629);

    auto tr_y_xxxxy_xxzzzz = pbuffer.data(idx_dip_hi + 630);

    auto tr_y_xxxxy_xyyyyy = pbuffer.data(idx_dip_hi + 631);

    auto tr_y_xxxxy_xyyyyz = pbuffer.data(idx_dip_hi + 632);

    auto tr_y_xxxxy_xyyyzz = pbuffer.data(idx_dip_hi + 633);

    auto tr_y_xxxxy_xyyzzz = pbuffer.data(idx_dip_hi + 634);

    auto tr_y_xxxxy_xyzzzz = pbuffer.data(idx_dip_hi + 635);

    auto tr_y_xxxxy_xzzzzz = pbuffer.data(idx_dip_hi + 636);

    auto tr_y_xxxxy_yyyyyy = pbuffer.data(idx_dip_hi + 637);

    auto tr_y_xxxxy_yyyyyz = pbuffer.data(idx_dip_hi + 638);

    auto tr_y_xxxxy_yyyyzz = pbuffer.data(idx_dip_hi + 639);

    auto tr_y_xxxxy_yyyzzz = pbuffer.data(idx_dip_hi + 640);

    auto tr_y_xxxxy_yyzzzz = pbuffer.data(idx_dip_hi + 641);

    auto tr_y_xxxxy_yzzzzz = pbuffer.data(idx_dip_hi + 642);

    auto tr_y_xxxxy_zzzzzz = pbuffer.data(idx_dip_hi + 643);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_y_xxxx_xxxxxx,  \
                             tr_y_xxxx_xxxxxz,  \
                             tr_y_xxxx_xxxxzz,  \
                             tr_y_xxxx_xxxzzz,  \
                             tr_y_xxxx_xxzzzz,  \
                             tr_y_xxxx_xzzzzz,  \
                             tr_y_xxxxy_xxxxxx, \
                             tr_y_xxxxy_xxxxxy, \
                             tr_y_xxxxy_xxxxxz, \
                             tr_y_xxxxy_xxxxyy, \
                             tr_y_xxxxy_xxxxyz, \
                             tr_y_xxxxy_xxxxzz, \
                             tr_y_xxxxy_xxxyyy, \
                             tr_y_xxxxy_xxxyyz, \
                             tr_y_xxxxy_xxxyzz, \
                             tr_y_xxxxy_xxxzzz, \
                             tr_y_xxxxy_xxyyyy, \
                             tr_y_xxxxy_xxyyyz, \
                             tr_y_xxxxy_xxyyzz, \
                             tr_y_xxxxy_xxyzzz, \
                             tr_y_xxxxy_xxzzzz, \
                             tr_y_xxxxy_xyyyyy, \
                             tr_y_xxxxy_xyyyyz, \
                             tr_y_xxxxy_xyyyzz, \
                             tr_y_xxxxy_xyyzzz, \
                             tr_y_xxxxy_xyzzzz, \
                             tr_y_xxxxy_xzzzzz, \
                             tr_y_xxxxy_yyyyyy, \
                             tr_y_xxxxy_yyyyyz, \
                             tr_y_xxxxy_yyyyzz, \
                             tr_y_xxxxy_yyyzzz, \
                             tr_y_xxxxy_yyzzzz, \
                             tr_y_xxxxy_yzzzzz, \
                             tr_y_xxxxy_zzzzzz, \
                             tr_y_xxxy_xxxxxy,  \
                             tr_y_xxxy_xxxxy,   \
                             tr_y_xxxy_xxxxyy,  \
                             tr_y_xxxy_xxxxyz,  \
                             tr_y_xxxy_xxxyy,   \
                             tr_y_xxxy_xxxyyy,  \
                             tr_y_xxxy_xxxyyz,  \
                             tr_y_xxxy_xxxyz,   \
                             tr_y_xxxy_xxxyzz,  \
                             tr_y_xxxy_xxyyy,   \
                             tr_y_xxxy_xxyyyy,  \
                             tr_y_xxxy_xxyyyz,  \
                             tr_y_xxxy_xxyyz,   \
                             tr_y_xxxy_xxyyzz,  \
                             tr_y_xxxy_xxyzz,   \
                             tr_y_xxxy_xxyzzz,  \
                             tr_y_xxxy_xyyyy,   \
                             tr_y_xxxy_xyyyyy,  \
                             tr_y_xxxy_xyyyyz,  \
                             tr_y_xxxy_xyyyz,   \
                             tr_y_xxxy_xyyyzz,  \
                             tr_y_xxxy_xyyzz,   \
                             tr_y_xxxy_xyyzzz,  \
                             tr_y_xxxy_xyzzz,   \
                             tr_y_xxxy_xyzzzz,  \
                             tr_y_xxxy_yyyyy,   \
                             tr_y_xxxy_yyyyyy,  \
                             tr_y_xxxy_yyyyyz,  \
                             tr_y_xxxy_yyyyz,   \
                             tr_y_xxxy_yyyyzz,  \
                             tr_y_xxxy_yyyzz,   \
                             tr_y_xxxy_yyyzzz,  \
                             tr_y_xxxy_yyzzz,   \
                             tr_y_xxxy_yyzzzz,  \
                             tr_y_xxxy_yzzzz,   \
                             tr_y_xxxy_yzzzzz,  \
                             tr_y_xxxy_zzzzzz,  \
                             tr_y_xxy_xxxxxy,   \
                             tr_y_xxy_xxxxyy,   \
                             tr_y_xxy_xxxxyz,   \
                             tr_y_xxy_xxxyyy,   \
                             tr_y_xxy_xxxyyz,   \
                             tr_y_xxy_xxxyzz,   \
                             tr_y_xxy_xxyyyy,   \
                             tr_y_xxy_xxyyyz,   \
                             tr_y_xxy_xxyyzz,   \
                             tr_y_xxy_xxyzzz,   \
                             tr_y_xxy_xyyyyy,   \
                             tr_y_xxy_xyyyyz,   \
                             tr_y_xxy_xyyyzz,   \
                             tr_y_xxy_xyyzzz,   \
                             tr_y_xxy_xyzzzz,   \
                             tr_y_xxy_yyyyyy,   \
                             tr_y_xxy_yyyyyz,   \
                             tr_y_xxy_yyyyzz,   \
                             tr_y_xxy_yyyzzz,   \
                             tr_y_xxy_yyzzzz,   \
                             tr_y_xxy_yzzzzz,   \
                             tr_y_xxy_zzzzzz,   \
                             ts_xxxx_xxxxxx,    \
                             ts_xxxx_xxxxxz,    \
                             ts_xxxx_xxxxzz,    \
                             ts_xxxx_xxxzzz,    \
                             ts_xxxx_xxzzzz,    \
                             ts_xxxx_xzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxy_xxxxxx[i] = ts_xxxx_xxxxxx[i] * fe_0 + tr_y_xxxx_xxxxxx[i] * pa_y[i];

        tr_y_xxxxy_xxxxxy[i] = 3.0 * tr_y_xxy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxxy_xxxxy[i] * fe_0 + tr_y_xxxy_xxxxxy[i] * pa_x[i];

        tr_y_xxxxy_xxxxxz[i] = ts_xxxx_xxxxxz[i] * fe_0 + tr_y_xxxx_xxxxxz[i] * pa_y[i];

        tr_y_xxxxy_xxxxyy[i] = 3.0 * tr_y_xxy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxxy_xxxyy[i] * fe_0 + tr_y_xxxy_xxxxyy[i] * pa_x[i];

        tr_y_xxxxy_xxxxyz[i] = 3.0 * tr_y_xxy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxxy_xxxyz[i] * fe_0 + tr_y_xxxy_xxxxyz[i] * pa_x[i];

        tr_y_xxxxy_xxxxzz[i] = ts_xxxx_xxxxzz[i] * fe_0 + tr_y_xxxx_xxxxzz[i] * pa_y[i];

        tr_y_xxxxy_xxxyyy[i] = 3.0 * tr_y_xxy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxxy_xxyyy[i] * fe_0 + tr_y_xxxy_xxxyyy[i] * pa_x[i];

        tr_y_xxxxy_xxxyyz[i] = 3.0 * tr_y_xxy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxxy_xxyyz[i] * fe_0 + tr_y_xxxy_xxxyyz[i] * pa_x[i];

        tr_y_xxxxy_xxxyzz[i] = 3.0 * tr_y_xxy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxxy_xxyzz[i] * fe_0 + tr_y_xxxy_xxxyzz[i] * pa_x[i];

        tr_y_xxxxy_xxxzzz[i] = ts_xxxx_xxxzzz[i] * fe_0 + tr_y_xxxx_xxxzzz[i] * pa_y[i];

        tr_y_xxxxy_xxyyyy[i] = 3.0 * tr_y_xxy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxxy_xyyyy[i] * fe_0 + tr_y_xxxy_xxyyyy[i] * pa_x[i];

        tr_y_xxxxy_xxyyyz[i] = 3.0 * tr_y_xxy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxxy_xyyyz[i] * fe_0 + tr_y_xxxy_xxyyyz[i] * pa_x[i];

        tr_y_xxxxy_xxyyzz[i] = 3.0 * tr_y_xxy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxxy_xyyzz[i] * fe_0 + tr_y_xxxy_xxyyzz[i] * pa_x[i];

        tr_y_xxxxy_xxyzzz[i] = 3.0 * tr_y_xxy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxxy_xyzzz[i] * fe_0 + tr_y_xxxy_xxyzzz[i] * pa_x[i];

        tr_y_xxxxy_xxzzzz[i] = ts_xxxx_xxzzzz[i] * fe_0 + tr_y_xxxx_xxzzzz[i] * pa_y[i];

        tr_y_xxxxy_xyyyyy[i] = 3.0 * tr_y_xxy_xyyyyy[i] * fe_0 + tr_y_xxxy_yyyyy[i] * fe_0 + tr_y_xxxy_xyyyyy[i] * pa_x[i];

        tr_y_xxxxy_xyyyyz[i] = 3.0 * tr_y_xxy_xyyyyz[i] * fe_0 + tr_y_xxxy_yyyyz[i] * fe_0 + tr_y_xxxy_xyyyyz[i] * pa_x[i];

        tr_y_xxxxy_xyyyzz[i] = 3.0 * tr_y_xxy_xyyyzz[i] * fe_0 + tr_y_xxxy_yyyzz[i] * fe_0 + tr_y_xxxy_xyyyzz[i] * pa_x[i];

        tr_y_xxxxy_xyyzzz[i] = 3.0 * tr_y_xxy_xyyzzz[i] * fe_0 + tr_y_xxxy_yyzzz[i] * fe_0 + tr_y_xxxy_xyyzzz[i] * pa_x[i];

        tr_y_xxxxy_xyzzzz[i] = 3.0 * tr_y_xxy_xyzzzz[i] * fe_0 + tr_y_xxxy_yzzzz[i] * fe_0 + tr_y_xxxy_xyzzzz[i] * pa_x[i];

        tr_y_xxxxy_xzzzzz[i] = ts_xxxx_xzzzzz[i] * fe_0 + tr_y_xxxx_xzzzzz[i] * pa_y[i];

        tr_y_xxxxy_yyyyyy[i] = 3.0 * tr_y_xxy_yyyyyy[i] * fe_0 + tr_y_xxxy_yyyyyy[i] * pa_x[i];

        tr_y_xxxxy_yyyyyz[i] = 3.0 * tr_y_xxy_yyyyyz[i] * fe_0 + tr_y_xxxy_yyyyyz[i] * pa_x[i];

        tr_y_xxxxy_yyyyzz[i] = 3.0 * tr_y_xxy_yyyyzz[i] * fe_0 + tr_y_xxxy_yyyyzz[i] * pa_x[i];

        tr_y_xxxxy_yyyzzz[i] = 3.0 * tr_y_xxy_yyyzzz[i] * fe_0 + tr_y_xxxy_yyyzzz[i] * pa_x[i];

        tr_y_xxxxy_yyzzzz[i] = 3.0 * tr_y_xxy_yyzzzz[i] * fe_0 + tr_y_xxxy_yyzzzz[i] * pa_x[i];

        tr_y_xxxxy_yzzzzz[i] = 3.0 * tr_y_xxy_yzzzzz[i] * fe_0 + tr_y_xxxy_yzzzzz[i] * pa_x[i];

        tr_y_xxxxy_zzzzzz[i] = 3.0 * tr_y_xxy_zzzzzz[i] * fe_0 + tr_y_xxxy_zzzzzz[i] * pa_x[i];
    }

    // Set up 644-672 components of targeted buffer : HI

    auto tr_y_xxxxz_xxxxxx = pbuffer.data(idx_dip_hi + 644);

    auto tr_y_xxxxz_xxxxxy = pbuffer.data(idx_dip_hi + 645);

    auto tr_y_xxxxz_xxxxxz = pbuffer.data(idx_dip_hi + 646);

    auto tr_y_xxxxz_xxxxyy = pbuffer.data(idx_dip_hi + 647);

    auto tr_y_xxxxz_xxxxyz = pbuffer.data(idx_dip_hi + 648);

    auto tr_y_xxxxz_xxxxzz = pbuffer.data(idx_dip_hi + 649);

    auto tr_y_xxxxz_xxxyyy = pbuffer.data(idx_dip_hi + 650);

    auto tr_y_xxxxz_xxxyyz = pbuffer.data(idx_dip_hi + 651);

    auto tr_y_xxxxz_xxxyzz = pbuffer.data(idx_dip_hi + 652);

    auto tr_y_xxxxz_xxxzzz = pbuffer.data(idx_dip_hi + 653);

    auto tr_y_xxxxz_xxyyyy = pbuffer.data(idx_dip_hi + 654);

    auto tr_y_xxxxz_xxyyyz = pbuffer.data(idx_dip_hi + 655);

    auto tr_y_xxxxz_xxyyzz = pbuffer.data(idx_dip_hi + 656);

    auto tr_y_xxxxz_xxyzzz = pbuffer.data(idx_dip_hi + 657);

    auto tr_y_xxxxz_xxzzzz = pbuffer.data(idx_dip_hi + 658);

    auto tr_y_xxxxz_xyyyyy = pbuffer.data(idx_dip_hi + 659);

    auto tr_y_xxxxz_xyyyyz = pbuffer.data(idx_dip_hi + 660);

    auto tr_y_xxxxz_xyyyzz = pbuffer.data(idx_dip_hi + 661);

    auto tr_y_xxxxz_xyyzzz = pbuffer.data(idx_dip_hi + 662);

    auto tr_y_xxxxz_xyzzzz = pbuffer.data(idx_dip_hi + 663);

    auto tr_y_xxxxz_xzzzzz = pbuffer.data(idx_dip_hi + 664);

    auto tr_y_xxxxz_yyyyyy = pbuffer.data(idx_dip_hi + 665);

    auto tr_y_xxxxz_yyyyyz = pbuffer.data(idx_dip_hi + 666);

    auto tr_y_xxxxz_yyyyzz = pbuffer.data(idx_dip_hi + 667);

    auto tr_y_xxxxz_yyyzzz = pbuffer.data(idx_dip_hi + 668);

    auto tr_y_xxxxz_yyzzzz = pbuffer.data(idx_dip_hi + 669);

    auto tr_y_xxxxz_yzzzzz = pbuffer.data(idx_dip_hi + 670);

    auto tr_y_xxxxz_zzzzzz = pbuffer.data(idx_dip_hi + 671);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_y_xxxx_xxxxx,   \
                             tr_y_xxxx_xxxxxx,  \
                             tr_y_xxxx_xxxxxy,  \
                             tr_y_xxxx_xxxxxz,  \
                             tr_y_xxxx_xxxxy,   \
                             tr_y_xxxx_xxxxyy,  \
                             tr_y_xxxx_xxxxyz,  \
                             tr_y_xxxx_xxxxz,   \
                             tr_y_xxxx_xxxxzz,  \
                             tr_y_xxxx_xxxyy,   \
                             tr_y_xxxx_xxxyyy,  \
                             tr_y_xxxx_xxxyyz,  \
                             tr_y_xxxx_xxxyz,   \
                             tr_y_xxxx_xxxyzz,  \
                             tr_y_xxxx_xxxzz,   \
                             tr_y_xxxx_xxxzzz,  \
                             tr_y_xxxx_xxyyy,   \
                             tr_y_xxxx_xxyyyy,  \
                             tr_y_xxxx_xxyyyz,  \
                             tr_y_xxxx_xxyyz,   \
                             tr_y_xxxx_xxyyzz,  \
                             tr_y_xxxx_xxyzz,   \
                             tr_y_xxxx_xxyzzz,  \
                             tr_y_xxxx_xxzzz,   \
                             tr_y_xxxx_xxzzzz,  \
                             tr_y_xxxx_xyyyy,   \
                             tr_y_xxxx_xyyyyy,  \
                             tr_y_xxxx_xyyyyz,  \
                             tr_y_xxxx_xyyyz,   \
                             tr_y_xxxx_xyyyzz,  \
                             tr_y_xxxx_xyyzz,   \
                             tr_y_xxxx_xyyzzz,  \
                             tr_y_xxxx_xyzzz,   \
                             tr_y_xxxx_xyzzzz,  \
                             tr_y_xxxx_xzzzz,   \
                             tr_y_xxxx_xzzzzz,  \
                             tr_y_xxxx_yyyyyy,  \
                             tr_y_xxxxz_xxxxxx, \
                             tr_y_xxxxz_xxxxxy, \
                             tr_y_xxxxz_xxxxxz, \
                             tr_y_xxxxz_xxxxyy, \
                             tr_y_xxxxz_xxxxyz, \
                             tr_y_xxxxz_xxxxzz, \
                             tr_y_xxxxz_xxxyyy, \
                             tr_y_xxxxz_xxxyyz, \
                             tr_y_xxxxz_xxxyzz, \
                             tr_y_xxxxz_xxxzzz, \
                             tr_y_xxxxz_xxyyyy, \
                             tr_y_xxxxz_xxyyyz, \
                             tr_y_xxxxz_xxyyzz, \
                             tr_y_xxxxz_xxyzzz, \
                             tr_y_xxxxz_xxzzzz, \
                             tr_y_xxxxz_xyyyyy, \
                             tr_y_xxxxz_xyyyyz, \
                             tr_y_xxxxz_xyyyzz, \
                             tr_y_xxxxz_xyyzzz, \
                             tr_y_xxxxz_xyzzzz, \
                             tr_y_xxxxz_xzzzzz, \
                             tr_y_xxxxz_yyyyyy, \
                             tr_y_xxxxz_yyyyyz, \
                             tr_y_xxxxz_yyyyzz, \
                             tr_y_xxxxz_yyyzzz, \
                             tr_y_xxxxz_yyzzzz, \
                             tr_y_xxxxz_yzzzzz, \
                             tr_y_xxxxz_zzzzzz, \
                             tr_y_xxxz_yyyyyz,  \
                             tr_y_xxxz_yyyyzz,  \
                             tr_y_xxxz_yyyzzz,  \
                             tr_y_xxxz_yyzzzz,  \
                             tr_y_xxxz_yzzzzz,  \
                             tr_y_xxxz_zzzzzz,  \
                             tr_y_xxz_yyyyyz,   \
                             tr_y_xxz_yyyyzz,   \
                             tr_y_xxz_yyyzzz,   \
                             tr_y_xxz_yyzzzz,   \
                             tr_y_xxz_yzzzzz,   \
                             tr_y_xxz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxz_xxxxxx[i] = tr_y_xxxx_xxxxxx[i] * pa_z[i];

        tr_y_xxxxz_xxxxxy[i] = tr_y_xxxx_xxxxxy[i] * pa_z[i];

        tr_y_xxxxz_xxxxxz[i] = tr_y_xxxx_xxxxx[i] * fe_0 + tr_y_xxxx_xxxxxz[i] * pa_z[i];

        tr_y_xxxxz_xxxxyy[i] = tr_y_xxxx_xxxxyy[i] * pa_z[i];

        tr_y_xxxxz_xxxxyz[i] = tr_y_xxxx_xxxxy[i] * fe_0 + tr_y_xxxx_xxxxyz[i] * pa_z[i];

        tr_y_xxxxz_xxxxzz[i] = 2.0 * tr_y_xxxx_xxxxz[i] * fe_0 + tr_y_xxxx_xxxxzz[i] * pa_z[i];

        tr_y_xxxxz_xxxyyy[i] = tr_y_xxxx_xxxyyy[i] * pa_z[i];

        tr_y_xxxxz_xxxyyz[i] = tr_y_xxxx_xxxyy[i] * fe_0 + tr_y_xxxx_xxxyyz[i] * pa_z[i];

        tr_y_xxxxz_xxxyzz[i] = 2.0 * tr_y_xxxx_xxxyz[i] * fe_0 + tr_y_xxxx_xxxyzz[i] * pa_z[i];

        tr_y_xxxxz_xxxzzz[i] = 3.0 * tr_y_xxxx_xxxzz[i] * fe_0 + tr_y_xxxx_xxxzzz[i] * pa_z[i];

        tr_y_xxxxz_xxyyyy[i] = tr_y_xxxx_xxyyyy[i] * pa_z[i];

        tr_y_xxxxz_xxyyyz[i] = tr_y_xxxx_xxyyy[i] * fe_0 + tr_y_xxxx_xxyyyz[i] * pa_z[i];

        tr_y_xxxxz_xxyyzz[i] = 2.0 * tr_y_xxxx_xxyyz[i] * fe_0 + tr_y_xxxx_xxyyzz[i] * pa_z[i];

        tr_y_xxxxz_xxyzzz[i] = 3.0 * tr_y_xxxx_xxyzz[i] * fe_0 + tr_y_xxxx_xxyzzz[i] * pa_z[i];

        tr_y_xxxxz_xxzzzz[i] = 4.0 * tr_y_xxxx_xxzzz[i] * fe_0 + tr_y_xxxx_xxzzzz[i] * pa_z[i];

        tr_y_xxxxz_xyyyyy[i] = tr_y_xxxx_xyyyyy[i] * pa_z[i];

        tr_y_xxxxz_xyyyyz[i] = tr_y_xxxx_xyyyy[i] * fe_0 + tr_y_xxxx_xyyyyz[i] * pa_z[i];

        tr_y_xxxxz_xyyyzz[i] = 2.0 * tr_y_xxxx_xyyyz[i] * fe_0 + tr_y_xxxx_xyyyzz[i] * pa_z[i];

        tr_y_xxxxz_xyyzzz[i] = 3.0 * tr_y_xxxx_xyyzz[i] * fe_0 + tr_y_xxxx_xyyzzz[i] * pa_z[i];

        tr_y_xxxxz_xyzzzz[i] = 4.0 * tr_y_xxxx_xyzzz[i] * fe_0 + tr_y_xxxx_xyzzzz[i] * pa_z[i];

        tr_y_xxxxz_xzzzzz[i] = 5.0 * tr_y_xxxx_xzzzz[i] * fe_0 + tr_y_xxxx_xzzzzz[i] * pa_z[i];

        tr_y_xxxxz_yyyyyy[i] = tr_y_xxxx_yyyyyy[i] * pa_z[i];

        tr_y_xxxxz_yyyyyz[i] = 3.0 * tr_y_xxz_yyyyyz[i] * fe_0 + tr_y_xxxz_yyyyyz[i] * pa_x[i];

        tr_y_xxxxz_yyyyzz[i] = 3.0 * tr_y_xxz_yyyyzz[i] * fe_0 + tr_y_xxxz_yyyyzz[i] * pa_x[i];

        tr_y_xxxxz_yyyzzz[i] = 3.0 * tr_y_xxz_yyyzzz[i] * fe_0 + tr_y_xxxz_yyyzzz[i] * pa_x[i];

        tr_y_xxxxz_yyzzzz[i] = 3.0 * tr_y_xxz_yyzzzz[i] * fe_0 + tr_y_xxxz_yyzzzz[i] * pa_x[i];

        tr_y_xxxxz_yzzzzz[i] = 3.0 * tr_y_xxz_yzzzzz[i] * fe_0 + tr_y_xxxz_yzzzzz[i] * pa_x[i];

        tr_y_xxxxz_zzzzzz[i] = 3.0 * tr_y_xxz_zzzzzz[i] * fe_0 + tr_y_xxxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 672-700 components of targeted buffer : HI

    auto tr_y_xxxyy_xxxxxx = pbuffer.data(idx_dip_hi + 672);

    auto tr_y_xxxyy_xxxxxy = pbuffer.data(idx_dip_hi + 673);

    auto tr_y_xxxyy_xxxxxz = pbuffer.data(idx_dip_hi + 674);

    auto tr_y_xxxyy_xxxxyy = pbuffer.data(idx_dip_hi + 675);

    auto tr_y_xxxyy_xxxxyz = pbuffer.data(idx_dip_hi + 676);

    auto tr_y_xxxyy_xxxxzz = pbuffer.data(idx_dip_hi + 677);

    auto tr_y_xxxyy_xxxyyy = pbuffer.data(idx_dip_hi + 678);

    auto tr_y_xxxyy_xxxyyz = pbuffer.data(idx_dip_hi + 679);

    auto tr_y_xxxyy_xxxyzz = pbuffer.data(idx_dip_hi + 680);

    auto tr_y_xxxyy_xxxzzz = pbuffer.data(idx_dip_hi + 681);

    auto tr_y_xxxyy_xxyyyy = pbuffer.data(idx_dip_hi + 682);

    auto tr_y_xxxyy_xxyyyz = pbuffer.data(idx_dip_hi + 683);

    auto tr_y_xxxyy_xxyyzz = pbuffer.data(idx_dip_hi + 684);

    auto tr_y_xxxyy_xxyzzz = pbuffer.data(idx_dip_hi + 685);

    auto tr_y_xxxyy_xxzzzz = pbuffer.data(idx_dip_hi + 686);

    auto tr_y_xxxyy_xyyyyy = pbuffer.data(idx_dip_hi + 687);

    auto tr_y_xxxyy_xyyyyz = pbuffer.data(idx_dip_hi + 688);

    auto tr_y_xxxyy_xyyyzz = pbuffer.data(idx_dip_hi + 689);

    auto tr_y_xxxyy_xyyzzz = pbuffer.data(idx_dip_hi + 690);

    auto tr_y_xxxyy_xyzzzz = pbuffer.data(idx_dip_hi + 691);

    auto tr_y_xxxyy_xzzzzz = pbuffer.data(idx_dip_hi + 692);

    auto tr_y_xxxyy_yyyyyy = pbuffer.data(idx_dip_hi + 693);

    auto tr_y_xxxyy_yyyyyz = pbuffer.data(idx_dip_hi + 694);

    auto tr_y_xxxyy_yyyyzz = pbuffer.data(idx_dip_hi + 695);

    auto tr_y_xxxyy_yyyzzz = pbuffer.data(idx_dip_hi + 696);

    auto tr_y_xxxyy_yyzzzz = pbuffer.data(idx_dip_hi + 697);

    auto tr_y_xxxyy_yzzzzz = pbuffer.data(idx_dip_hi + 698);

    auto tr_y_xxxyy_zzzzzz = pbuffer.data(idx_dip_hi + 699);

#pragma omp simd aligned(pa_x,                  \
                             tr_y_xxxyy_xxxxxx, \
                             tr_y_xxxyy_xxxxxy, \
                             tr_y_xxxyy_xxxxxz, \
                             tr_y_xxxyy_xxxxyy, \
                             tr_y_xxxyy_xxxxyz, \
                             tr_y_xxxyy_xxxxzz, \
                             tr_y_xxxyy_xxxyyy, \
                             tr_y_xxxyy_xxxyyz, \
                             tr_y_xxxyy_xxxyzz, \
                             tr_y_xxxyy_xxxzzz, \
                             tr_y_xxxyy_xxyyyy, \
                             tr_y_xxxyy_xxyyyz, \
                             tr_y_xxxyy_xxyyzz, \
                             tr_y_xxxyy_xxyzzz, \
                             tr_y_xxxyy_xxzzzz, \
                             tr_y_xxxyy_xyyyyy, \
                             tr_y_xxxyy_xyyyyz, \
                             tr_y_xxxyy_xyyyzz, \
                             tr_y_xxxyy_xyyzzz, \
                             tr_y_xxxyy_xyzzzz, \
                             tr_y_xxxyy_xzzzzz, \
                             tr_y_xxxyy_yyyyyy, \
                             tr_y_xxxyy_yyyyyz, \
                             tr_y_xxxyy_yyyyzz, \
                             tr_y_xxxyy_yyyzzz, \
                             tr_y_xxxyy_yyzzzz, \
                             tr_y_xxxyy_yzzzzz, \
                             tr_y_xxxyy_zzzzzz, \
                             tr_y_xxyy_xxxxx,   \
                             tr_y_xxyy_xxxxxx,  \
                             tr_y_xxyy_xxxxxy,  \
                             tr_y_xxyy_xxxxxz,  \
                             tr_y_xxyy_xxxxy,   \
                             tr_y_xxyy_xxxxyy,  \
                             tr_y_xxyy_xxxxyz,  \
                             tr_y_xxyy_xxxxz,   \
                             tr_y_xxyy_xxxxzz,  \
                             tr_y_xxyy_xxxyy,   \
                             tr_y_xxyy_xxxyyy,  \
                             tr_y_xxyy_xxxyyz,  \
                             tr_y_xxyy_xxxyz,   \
                             tr_y_xxyy_xxxyzz,  \
                             tr_y_xxyy_xxxzz,   \
                             tr_y_xxyy_xxxzzz,  \
                             tr_y_xxyy_xxyyy,   \
                             tr_y_xxyy_xxyyyy,  \
                             tr_y_xxyy_xxyyyz,  \
                             tr_y_xxyy_xxyyz,   \
                             tr_y_xxyy_xxyyzz,  \
                             tr_y_xxyy_xxyzz,   \
                             tr_y_xxyy_xxyzzz,  \
                             tr_y_xxyy_xxzzz,   \
                             tr_y_xxyy_xxzzzz,  \
                             tr_y_xxyy_xyyyy,   \
                             tr_y_xxyy_xyyyyy,  \
                             tr_y_xxyy_xyyyyz,  \
                             tr_y_xxyy_xyyyz,   \
                             tr_y_xxyy_xyyyzz,  \
                             tr_y_xxyy_xyyzz,   \
                             tr_y_xxyy_xyyzzz,  \
                             tr_y_xxyy_xyzzz,   \
                             tr_y_xxyy_xyzzzz,  \
                             tr_y_xxyy_xzzzz,   \
                             tr_y_xxyy_xzzzzz,  \
                             tr_y_xxyy_yyyyy,   \
                             tr_y_xxyy_yyyyyy,  \
                             tr_y_xxyy_yyyyyz,  \
                             tr_y_xxyy_yyyyz,   \
                             tr_y_xxyy_yyyyzz,  \
                             tr_y_xxyy_yyyzz,   \
                             tr_y_xxyy_yyyzzz,  \
                             tr_y_xxyy_yyzzz,   \
                             tr_y_xxyy_yyzzzz,  \
                             tr_y_xxyy_yzzzz,   \
                             tr_y_xxyy_yzzzzz,  \
                             tr_y_xxyy_zzzzz,   \
                             tr_y_xxyy_zzzzzz,  \
                             tr_y_xyy_xxxxxx,   \
                             tr_y_xyy_xxxxxy,   \
                             tr_y_xyy_xxxxxz,   \
                             tr_y_xyy_xxxxyy,   \
                             tr_y_xyy_xxxxyz,   \
                             tr_y_xyy_xxxxzz,   \
                             tr_y_xyy_xxxyyy,   \
                             tr_y_xyy_xxxyyz,   \
                             tr_y_xyy_xxxyzz,   \
                             tr_y_xyy_xxxzzz,   \
                             tr_y_xyy_xxyyyy,   \
                             tr_y_xyy_xxyyyz,   \
                             tr_y_xyy_xxyyzz,   \
                             tr_y_xyy_xxyzzz,   \
                             tr_y_xyy_xxzzzz,   \
                             tr_y_xyy_xyyyyy,   \
                             tr_y_xyy_xyyyyz,   \
                             tr_y_xyy_xyyyzz,   \
                             tr_y_xyy_xyyzzz,   \
                             tr_y_xyy_xyzzzz,   \
                             tr_y_xyy_xzzzzz,   \
                             tr_y_xyy_yyyyyy,   \
                             tr_y_xyy_yyyyyz,   \
                             tr_y_xyy_yyyyzz,   \
                             tr_y_xyy_yyyzzz,   \
                             tr_y_xyy_yyzzzz,   \
                             tr_y_xyy_yzzzzz,   \
                             tr_y_xyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyy_xxxxxx[i] = 2.0 * tr_y_xyy_xxxxxx[i] * fe_0 + 6.0 * tr_y_xxyy_xxxxx[i] * fe_0 + tr_y_xxyy_xxxxxx[i] * pa_x[i];

        tr_y_xxxyy_xxxxxy[i] = 2.0 * tr_y_xyy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xxyy_xxxxy[i] * fe_0 + tr_y_xxyy_xxxxxy[i] * pa_x[i];

        tr_y_xxxyy_xxxxxz[i] = 2.0 * tr_y_xyy_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxyy_xxxxz[i] * fe_0 + tr_y_xxyy_xxxxxz[i] * pa_x[i];

        tr_y_xxxyy_xxxxyy[i] = 2.0 * tr_y_xyy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xxyy_xxxyy[i] * fe_0 + tr_y_xxyy_xxxxyy[i] * pa_x[i];

        tr_y_xxxyy_xxxxyz[i] = 2.0 * tr_y_xyy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxyy_xxxyz[i] * fe_0 + tr_y_xxyy_xxxxyz[i] * pa_x[i];

        tr_y_xxxyy_xxxxzz[i] = 2.0 * tr_y_xyy_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxyy_xxxzz[i] * fe_0 + tr_y_xxyy_xxxxzz[i] * pa_x[i];

        tr_y_xxxyy_xxxyyy[i] = 2.0 * tr_y_xyy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xxyy_xxyyy[i] * fe_0 + tr_y_xxyy_xxxyyy[i] * pa_x[i];

        tr_y_xxxyy_xxxyyz[i] = 2.0 * tr_y_xyy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxyy_xxyyz[i] * fe_0 + tr_y_xxyy_xxxyyz[i] * pa_x[i];

        tr_y_xxxyy_xxxyzz[i] = 2.0 * tr_y_xyy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxyy_xxyzz[i] * fe_0 + tr_y_xxyy_xxxyzz[i] * pa_x[i];

        tr_y_xxxyy_xxxzzz[i] = 2.0 * tr_y_xyy_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxyy_xxzzz[i] * fe_0 + tr_y_xxyy_xxxzzz[i] * pa_x[i];

        tr_y_xxxyy_xxyyyy[i] = 2.0 * tr_y_xyy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xxyy_xyyyy[i] * fe_0 + tr_y_xxyy_xxyyyy[i] * pa_x[i];

        tr_y_xxxyy_xxyyyz[i] = 2.0 * tr_y_xyy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxyy_xyyyz[i] * fe_0 + tr_y_xxyy_xxyyyz[i] * pa_x[i];

        tr_y_xxxyy_xxyyzz[i] = 2.0 * tr_y_xyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxyy_xyyzz[i] * fe_0 + tr_y_xxyy_xxyyzz[i] * pa_x[i];

        tr_y_xxxyy_xxyzzz[i] = 2.0 * tr_y_xyy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxyy_xyzzz[i] * fe_0 + tr_y_xxyy_xxyzzz[i] * pa_x[i];

        tr_y_xxxyy_xxzzzz[i] = 2.0 * tr_y_xyy_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxyy_xzzzz[i] * fe_0 + tr_y_xxyy_xxzzzz[i] * pa_x[i];

        tr_y_xxxyy_xyyyyy[i] = 2.0 * tr_y_xyy_xyyyyy[i] * fe_0 + tr_y_xxyy_yyyyy[i] * fe_0 + tr_y_xxyy_xyyyyy[i] * pa_x[i];

        tr_y_xxxyy_xyyyyz[i] = 2.0 * tr_y_xyy_xyyyyz[i] * fe_0 + tr_y_xxyy_yyyyz[i] * fe_0 + tr_y_xxyy_xyyyyz[i] * pa_x[i];

        tr_y_xxxyy_xyyyzz[i] = 2.0 * tr_y_xyy_xyyyzz[i] * fe_0 + tr_y_xxyy_yyyzz[i] * fe_0 + tr_y_xxyy_xyyyzz[i] * pa_x[i];

        tr_y_xxxyy_xyyzzz[i] = 2.0 * tr_y_xyy_xyyzzz[i] * fe_0 + tr_y_xxyy_yyzzz[i] * fe_0 + tr_y_xxyy_xyyzzz[i] * pa_x[i];

        tr_y_xxxyy_xyzzzz[i] = 2.0 * tr_y_xyy_xyzzzz[i] * fe_0 + tr_y_xxyy_yzzzz[i] * fe_0 + tr_y_xxyy_xyzzzz[i] * pa_x[i];

        tr_y_xxxyy_xzzzzz[i] = 2.0 * tr_y_xyy_xzzzzz[i] * fe_0 + tr_y_xxyy_zzzzz[i] * fe_0 + tr_y_xxyy_xzzzzz[i] * pa_x[i];

        tr_y_xxxyy_yyyyyy[i] = 2.0 * tr_y_xyy_yyyyyy[i] * fe_0 + tr_y_xxyy_yyyyyy[i] * pa_x[i];

        tr_y_xxxyy_yyyyyz[i] = 2.0 * tr_y_xyy_yyyyyz[i] * fe_0 + tr_y_xxyy_yyyyyz[i] * pa_x[i];

        tr_y_xxxyy_yyyyzz[i] = 2.0 * tr_y_xyy_yyyyzz[i] * fe_0 + tr_y_xxyy_yyyyzz[i] * pa_x[i];

        tr_y_xxxyy_yyyzzz[i] = 2.0 * tr_y_xyy_yyyzzz[i] * fe_0 + tr_y_xxyy_yyyzzz[i] * pa_x[i];

        tr_y_xxxyy_yyzzzz[i] = 2.0 * tr_y_xyy_yyzzzz[i] * fe_0 + tr_y_xxyy_yyzzzz[i] * pa_x[i];

        tr_y_xxxyy_yzzzzz[i] = 2.0 * tr_y_xyy_yzzzzz[i] * fe_0 + tr_y_xxyy_yzzzzz[i] * pa_x[i];

        tr_y_xxxyy_zzzzzz[i] = 2.0 * tr_y_xyy_zzzzzz[i] * fe_0 + tr_y_xxyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 700-728 components of targeted buffer : HI

    auto tr_y_xxxyz_xxxxxx = pbuffer.data(idx_dip_hi + 700);

    auto tr_y_xxxyz_xxxxxy = pbuffer.data(idx_dip_hi + 701);

    auto tr_y_xxxyz_xxxxxz = pbuffer.data(idx_dip_hi + 702);

    auto tr_y_xxxyz_xxxxyy = pbuffer.data(idx_dip_hi + 703);

    auto tr_y_xxxyz_xxxxyz = pbuffer.data(idx_dip_hi + 704);

    auto tr_y_xxxyz_xxxxzz = pbuffer.data(idx_dip_hi + 705);

    auto tr_y_xxxyz_xxxyyy = pbuffer.data(idx_dip_hi + 706);

    auto tr_y_xxxyz_xxxyyz = pbuffer.data(idx_dip_hi + 707);

    auto tr_y_xxxyz_xxxyzz = pbuffer.data(idx_dip_hi + 708);

    auto tr_y_xxxyz_xxxzzz = pbuffer.data(idx_dip_hi + 709);

    auto tr_y_xxxyz_xxyyyy = pbuffer.data(idx_dip_hi + 710);

    auto tr_y_xxxyz_xxyyyz = pbuffer.data(idx_dip_hi + 711);

    auto tr_y_xxxyz_xxyyzz = pbuffer.data(idx_dip_hi + 712);

    auto tr_y_xxxyz_xxyzzz = pbuffer.data(idx_dip_hi + 713);

    auto tr_y_xxxyz_xxzzzz = pbuffer.data(idx_dip_hi + 714);

    auto tr_y_xxxyz_xyyyyy = pbuffer.data(idx_dip_hi + 715);

    auto tr_y_xxxyz_xyyyyz = pbuffer.data(idx_dip_hi + 716);

    auto tr_y_xxxyz_xyyyzz = pbuffer.data(idx_dip_hi + 717);

    auto tr_y_xxxyz_xyyzzz = pbuffer.data(idx_dip_hi + 718);

    auto tr_y_xxxyz_xyzzzz = pbuffer.data(idx_dip_hi + 719);

    auto tr_y_xxxyz_xzzzzz = pbuffer.data(idx_dip_hi + 720);

    auto tr_y_xxxyz_yyyyyy = pbuffer.data(idx_dip_hi + 721);

    auto tr_y_xxxyz_yyyyyz = pbuffer.data(idx_dip_hi + 722);

    auto tr_y_xxxyz_yyyyzz = pbuffer.data(idx_dip_hi + 723);

    auto tr_y_xxxyz_yyyzzz = pbuffer.data(idx_dip_hi + 724);

    auto tr_y_xxxyz_yyzzzz = pbuffer.data(idx_dip_hi + 725);

    auto tr_y_xxxyz_yzzzzz = pbuffer.data(idx_dip_hi + 726);

    auto tr_y_xxxyz_zzzzzz = pbuffer.data(idx_dip_hi + 727);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             tr_y_xxxy_xxxxxx,  \
                             tr_y_xxxy_xxxxxy,  \
                             tr_y_xxxy_xxxxy,   \
                             tr_y_xxxy_xxxxyy,  \
                             tr_y_xxxy_xxxxyz,  \
                             tr_y_xxxy_xxxyy,   \
                             tr_y_xxxy_xxxyyy,  \
                             tr_y_xxxy_xxxyyz,  \
                             tr_y_xxxy_xxxyz,   \
                             tr_y_xxxy_xxxyzz,  \
                             tr_y_xxxy_xxyyy,   \
                             tr_y_xxxy_xxyyyy,  \
                             tr_y_xxxy_xxyyyz,  \
                             tr_y_xxxy_xxyyz,   \
                             tr_y_xxxy_xxyyzz,  \
                             tr_y_xxxy_xxyzz,   \
                             tr_y_xxxy_xxyzzz,  \
                             tr_y_xxxy_xyyyy,   \
                             tr_y_xxxy_xyyyyy,  \
                             tr_y_xxxy_xyyyyz,  \
                             tr_y_xxxy_xyyyz,   \
                             tr_y_xxxy_xyyyzz,  \
                             tr_y_xxxy_xyyzz,   \
                             tr_y_xxxy_xyyzzz,  \
                             tr_y_xxxy_xyzzz,   \
                             tr_y_xxxy_xyzzzz,  \
                             tr_y_xxxy_yyyyyy,  \
                             tr_y_xxxyz_xxxxxx, \
                             tr_y_xxxyz_xxxxxy, \
                             tr_y_xxxyz_xxxxxz, \
                             tr_y_xxxyz_xxxxyy, \
                             tr_y_xxxyz_xxxxyz, \
                             tr_y_xxxyz_xxxxzz, \
                             tr_y_xxxyz_xxxyyy, \
                             tr_y_xxxyz_xxxyyz, \
                             tr_y_xxxyz_xxxyzz, \
                             tr_y_xxxyz_xxxzzz, \
                             tr_y_xxxyz_xxyyyy, \
                             tr_y_xxxyz_xxyyyz, \
                             tr_y_xxxyz_xxyyzz, \
                             tr_y_xxxyz_xxyzzz, \
                             tr_y_xxxyz_xxzzzz, \
                             tr_y_xxxyz_xyyyyy, \
                             tr_y_xxxyz_xyyyyz, \
                             tr_y_xxxyz_xyyyzz, \
                             tr_y_xxxyz_xyyzzz, \
                             tr_y_xxxyz_xyzzzz, \
                             tr_y_xxxyz_xzzzzz, \
                             tr_y_xxxyz_yyyyyy, \
                             tr_y_xxxyz_yyyyyz, \
                             tr_y_xxxyz_yyyyzz, \
                             tr_y_xxxyz_yyyzzz, \
                             tr_y_xxxyz_yyzzzz, \
                             tr_y_xxxyz_yzzzzz, \
                             tr_y_xxxyz_zzzzzz, \
                             tr_y_xxxz_xxxxxz,  \
                             tr_y_xxxz_xxxxzz,  \
                             tr_y_xxxz_xxxzzz,  \
                             tr_y_xxxz_xxzzzz,  \
                             tr_y_xxxz_xzzzzz,  \
                             tr_y_xxyz_yyyyyz,  \
                             tr_y_xxyz_yyyyzz,  \
                             tr_y_xxyz_yyyzzz,  \
                             tr_y_xxyz_yyzzzz,  \
                             tr_y_xxyz_yzzzzz,  \
                             tr_y_xxyz_zzzzzz,  \
                             tr_y_xyz_yyyyyz,   \
                             tr_y_xyz_yyyyzz,   \
                             tr_y_xyz_yyyzzz,   \
                             tr_y_xyz_yyzzzz,   \
                             tr_y_xyz_yzzzzz,   \
                             tr_y_xyz_zzzzzz,   \
                             ts_xxxz_xxxxxz,    \
                             ts_xxxz_xxxxzz,    \
                             ts_xxxz_xxxzzz,    \
                             ts_xxxz_xxzzzz,    \
                             ts_xxxz_xzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyz_xxxxxx[i] = tr_y_xxxy_xxxxxx[i] * pa_z[i];

        tr_y_xxxyz_xxxxxy[i] = tr_y_xxxy_xxxxxy[i] * pa_z[i];

        tr_y_xxxyz_xxxxxz[i] = ts_xxxz_xxxxxz[i] * fe_0 + tr_y_xxxz_xxxxxz[i] * pa_y[i];

        tr_y_xxxyz_xxxxyy[i] = tr_y_xxxy_xxxxyy[i] * pa_z[i];

        tr_y_xxxyz_xxxxyz[i] = tr_y_xxxy_xxxxy[i] * fe_0 + tr_y_xxxy_xxxxyz[i] * pa_z[i];

        tr_y_xxxyz_xxxxzz[i] = ts_xxxz_xxxxzz[i] * fe_0 + tr_y_xxxz_xxxxzz[i] * pa_y[i];

        tr_y_xxxyz_xxxyyy[i] = tr_y_xxxy_xxxyyy[i] * pa_z[i];

        tr_y_xxxyz_xxxyyz[i] = tr_y_xxxy_xxxyy[i] * fe_0 + tr_y_xxxy_xxxyyz[i] * pa_z[i];

        tr_y_xxxyz_xxxyzz[i] = 2.0 * tr_y_xxxy_xxxyz[i] * fe_0 + tr_y_xxxy_xxxyzz[i] * pa_z[i];

        tr_y_xxxyz_xxxzzz[i] = ts_xxxz_xxxzzz[i] * fe_0 + tr_y_xxxz_xxxzzz[i] * pa_y[i];

        tr_y_xxxyz_xxyyyy[i] = tr_y_xxxy_xxyyyy[i] * pa_z[i];

        tr_y_xxxyz_xxyyyz[i] = tr_y_xxxy_xxyyy[i] * fe_0 + tr_y_xxxy_xxyyyz[i] * pa_z[i];

        tr_y_xxxyz_xxyyzz[i] = 2.0 * tr_y_xxxy_xxyyz[i] * fe_0 + tr_y_xxxy_xxyyzz[i] * pa_z[i];

        tr_y_xxxyz_xxyzzz[i] = 3.0 * tr_y_xxxy_xxyzz[i] * fe_0 + tr_y_xxxy_xxyzzz[i] * pa_z[i];

        tr_y_xxxyz_xxzzzz[i] = ts_xxxz_xxzzzz[i] * fe_0 + tr_y_xxxz_xxzzzz[i] * pa_y[i];

        tr_y_xxxyz_xyyyyy[i] = tr_y_xxxy_xyyyyy[i] * pa_z[i];

        tr_y_xxxyz_xyyyyz[i] = tr_y_xxxy_xyyyy[i] * fe_0 + tr_y_xxxy_xyyyyz[i] * pa_z[i];

        tr_y_xxxyz_xyyyzz[i] = 2.0 * tr_y_xxxy_xyyyz[i] * fe_0 + tr_y_xxxy_xyyyzz[i] * pa_z[i];

        tr_y_xxxyz_xyyzzz[i] = 3.0 * tr_y_xxxy_xyyzz[i] * fe_0 + tr_y_xxxy_xyyzzz[i] * pa_z[i];

        tr_y_xxxyz_xyzzzz[i] = 4.0 * tr_y_xxxy_xyzzz[i] * fe_0 + tr_y_xxxy_xyzzzz[i] * pa_z[i];

        tr_y_xxxyz_xzzzzz[i] = ts_xxxz_xzzzzz[i] * fe_0 + tr_y_xxxz_xzzzzz[i] * pa_y[i];

        tr_y_xxxyz_yyyyyy[i] = tr_y_xxxy_yyyyyy[i] * pa_z[i];

        tr_y_xxxyz_yyyyyz[i] = 2.0 * tr_y_xyz_yyyyyz[i] * fe_0 + tr_y_xxyz_yyyyyz[i] * pa_x[i];

        tr_y_xxxyz_yyyyzz[i] = 2.0 * tr_y_xyz_yyyyzz[i] * fe_0 + tr_y_xxyz_yyyyzz[i] * pa_x[i];

        tr_y_xxxyz_yyyzzz[i] = 2.0 * tr_y_xyz_yyyzzz[i] * fe_0 + tr_y_xxyz_yyyzzz[i] * pa_x[i];

        tr_y_xxxyz_yyzzzz[i] = 2.0 * tr_y_xyz_yyzzzz[i] * fe_0 + tr_y_xxyz_yyzzzz[i] * pa_x[i];

        tr_y_xxxyz_yzzzzz[i] = 2.0 * tr_y_xyz_yzzzzz[i] * fe_0 + tr_y_xxyz_yzzzzz[i] * pa_x[i];

        tr_y_xxxyz_zzzzzz[i] = 2.0 * tr_y_xyz_zzzzzz[i] * fe_0 + tr_y_xxyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 728-756 components of targeted buffer : HI

    auto tr_y_xxxzz_xxxxxx = pbuffer.data(idx_dip_hi + 728);

    auto tr_y_xxxzz_xxxxxy = pbuffer.data(idx_dip_hi + 729);

    auto tr_y_xxxzz_xxxxxz = pbuffer.data(idx_dip_hi + 730);

    auto tr_y_xxxzz_xxxxyy = pbuffer.data(idx_dip_hi + 731);

    auto tr_y_xxxzz_xxxxyz = pbuffer.data(idx_dip_hi + 732);

    auto tr_y_xxxzz_xxxxzz = pbuffer.data(idx_dip_hi + 733);

    auto tr_y_xxxzz_xxxyyy = pbuffer.data(idx_dip_hi + 734);

    auto tr_y_xxxzz_xxxyyz = pbuffer.data(idx_dip_hi + 735);

    auto tr_y_xxxzz_xxxyzz = pbuffer.data(idx_dip_hi + 736);

    auto tr_y_xxxzz_xxxzzz = pbuffer.data(idx_dip_hi + 737);

    auto tr_y_xxxzz_xxyyyy = pbuffer.data(idx_dip_hi + 738);

    auto tr_y_xxxzz_xxyyyz = pbuffer.data(idx_dip_hi + 739);

    auto tr_y_xxxzz_xxyyzz = pbuffer.data(idx_dip_hi + 740);

    auto tr_y_xxxzz_xxyzzz = pbuffer.data(idx_dip_hi + 741);

    auto tr_y_xxxzz_xxzzzz = pbuffer.data(idx_dip_hi + 742);

    auto tr_y_xxxzz_xyyyyy = pbuffer.data(idx_dip_hi + 743);

    auto tr_y_xxxzz_xyyyyz = pbuffer.data(idx_dip_hi + 744);

    auto tr_y_xxxzz_xyyyzz = pbuffer.data(idx_dip_hi + 745);

    auto tr_y_xxxzz_xyyzzz = pbuffer.data(idx_dip_hi + 746);

    auto tr_y_xxxzz_xyzzzz = pbuffer.data(idx_dip_hi + 747);

    auto tr_y_xxxzz_xzzzzz = pbuffer.data(idx_dip_hi + 748);

    auto tr_y_xxxzz_yyyyyy = pbuffer.data(idx_dip_hi + 749);

    auto tr_y_xxxzz_yyyyyz = pbuffer.data(idx_dip_hi + 750);

    auto tr_y_xxxzz_yyyyzz = pbuffer.data(idx_dip_hi + 751);

    auto tr_y_xxxzz_yyyzzz = pbuffer.data(idx_dip_hi + 752);

    auto tr_y_xxxzz_yyzzzz = pbuffer.data(idx_dip_hi + 753);

    auto tr_y_xxxzz_yzzzzz = pbuffer.data(idx_dip_hi + 754);

    auto tr_y_xxxzz_zzzzzz = pbuffer.data(idx_dip_hi + 755);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_y_xxx_xxxxxx,   \
                             tr_y_xxx_xxxxxy,   \
                             tr_y_xxx_xxxxyy,   \
                             tr_y_xxx_xxxyyy,   \
                             tr_y_xxx_xxyyyy,   \
                             tr_y_xxx_xyyyyy,   \
                             tr_y_xxxz_xxxxxx,  \
                             tr_y_xxxz_xxxxxy,  \
                             tr_y_xxxz_xxxxyy,  \
                             tr_y_xxxz_xxxyyy,  \
                             tr_y_xxxz_xxyyyy,  \
                             tr_y_xxxz_xyyyyy,  \
                             tr_y_xxxzz_xxxxxx, \
                             tr_y_xxxzz_xxxxxy, \
                             tr_y_xxxzz_xxxxxz, \
                             tr_y_xxxzz_xxxxyy, \
                             tr_y_xxxzz_xxxxyz, \
                             tr_y_xxxzz_xxxxzz, \
                             tr_y_xxxzz_xxxyyy, \
                             tr_y_xxxzz_xxxyyz, \
                             tr_y_xxxzz_xxxyzz, \
                             tr_y_xxxzz_xxxzzz, \
                             tr_y_xxxzz_xxyyyy, \
                             tr_y_xxxzz_xxyyyz, \
                             tr_y_xxxzz_xxyyzz, \
                             tr_y_xxxzz_xxyzzz, \
                             tr_y_xxxzz_xxzzzz, \
                             tr_y_xxxzz_xyyyyy, \
                             tr_y_xxxzz_xyyyyz, \
                             tr_y_xxxzz_xyyyzz, \
                             tr_y_xxxzz_xyyzzz, \
                             tr_y_xxxzz_xyzzzz, \
                             tr_y_xxxzz_xzzzzz, \
                             tr_y_xxxzz_yyyyyy, \
                             tr_y_xxxzz_yyyyyz, \
                             tr_y_xxxzz_yyyyzz, \
                             tr_y_xxxzz_yyyzzz, \
                             tr_y_xxxzz_yyzzzz, \
                             tr_y_xxxzz_yzzzzz, \
                             tr_y_xxxzz_zzzzzz, \
                             tr_y_xxzz_xxxxxz,  \
                             tr_y_xxzz_xxxxyz,  \
                             tr_y_xxzz_xxxxz,   \
                             tr_y_xxzz_xxxxzz,  \
                             tr_y_xxzz_xxxyyz,  \
                             tr_y_xxzz_xxxyz,   \
                             tr_y_xxzz_xxxyzz,  \
                             tr_y_xxzz_xxxzz,   \
                             tr_y_xxzz_xxxzzz,  \
                             tr_y_xxzz_xxyyyz,  \
                             tr_y_xxzz_xxyyz,   \
                             tr_y_xxzz_xxyyzz,  \
                             tr_y_xxzz_xxyzz,   \
                             tr_y_xxzz_xxyzzz,  \
                             tr_y_xxzz_xxzzz,   \
                             tr_y_xxzz_xxzzzz,  \
                             tr_y_xxzz_xyyyyz,  \
                             tr_y_xxzz_xyyyz,   \
                             tr_y_xxzz_xyyyzz,  \
                             tr_y_xxzz_xyyzz,   \
                             tr_y_xxzz_xyyzzz,  \
                             tr_y_xxzz_xyzzz,   \
                             tr_y_xxzz_xyzzzz,  \
                             tr_y_xxzz_xzzzz,   \
                             tr_y_xxzz_xzzzzz,  \
                             tr_y_xxzz_yyyyyy,  \
                             tr_y_xxzz_yyyyyz,  \
                             tr_y_xxzz_yyyyz,   \
                             tr_y_xxzz_yyyyzz,  \
                             tr_y_xxzz_yyyzz,   \
                             tr_y_xxzz_yyyzzz,  \
                             tr_y_xxzz_yyzzz,   \
                             tr_y_xxzz_yyzzzz,  \
                             tr_y_xxzz_yzzzz,   \
                             tr_y_xxzz_yzzzzz,  \
                             tr_y_xxzz_zzzzz,   \
                             tr_y_xxzz_zzzzzz,  \
                             tr_y_xzz_xxxxxz,   \
                             tr_y_xzz_xxxxyz,   \
                             tr_y_xzz_xxxxzz,   \
                             tr_y_xzz_xxxyyz,   \
                             tr_y_xzz_xxxyzz,   \
                             tr_y_xzz_xxxzzz,   \
                             tr_y_xzz_xxyyyz,   \
                             tr_y_xzz_xxyyzz,   \
                             tr_y_xzz_xxyzzz,   \
                             tr_y_xzz_xxzzzz,   \
                             tr_y_xzz_xyyyyz,   \
                             tr_y_xzz_xyyyzz,   \
                             tr_y_xzz_xyyzzz,   \
                             tr_y_xzz_xyzzzz,   \
                             tr_y_xzz_xzzzzz,   \
                             tr_y_xzz_yyyyyy,   \
                             tr_y_xzz_yyyyyz,   \
                             tr_y_xzz_yyyyzz,   \
                             tr_y_xzz_yyyzzz,   \
                             tr_y_xzz_yyzzzz,   \
                             tr_y_xzz_yzzzzz,   \
                             tr_y_xzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzz_xxxxxx[i] = tr_y_xxx_xxxxxx[i] * fe_0 + tr_y_xxxz_xxxxxx[i] * pa_z[i];

        tr_y_xxxzz_xxxxxy[i] = tr_y_xxx_xxxxxy[i] * fe_0 + tr_y_xxxz_xxxxxy[i] * pa_z[i];

        tr_y_xxxzz_xxxxxz[i] = 2.0 * tr_y_xzz_xxxxxz[i] * fe_0 + 5.0 * tr_y_xxzz_xxxxz[i] * fe_0 + tr_y_xxzz_xxxxxz[i] * pa_x[i];

        tr_y_xxxzz_xxxxyy[i] = tr_y_xxx_xxxxyy[i] * fe_0 + tr_y_xxxz_xxxxyy[i] * pa_z[i];

        tr_y_xxxzz_xxxxyz[i] = 2.0 * tr_y_xzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xxzz_xxxyz[i] * fe_0 + tr_y_xxzz_xxxxyz[i] * pa_x[i];

        tr_y_xxxzz_xxxxzz[i] = 2.0 * tr_y_xzz_xxxxzz[i] * fe_0 + 4.0 * tr_y_xxzz_xxxzz[i] * fe_0 + tr_y_xxzz_xxxxzz[i] * pa_x[i];

        tr_y_xxxzz_xxxyyy[i] = tr_y_xxx_xxxyyy[i] * fe_0 + tr_y_xxxz_xxxyyy[i] * pa_z[i];

        tr_y_xxxzz_xxxyyz[i] = 2.0 * tr_y_xzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xxzz_xxyyz[i] * fe_0 + tr_y_xxzz_xxxyyz[i] * pa_x[i];

        tr_y_xxxzz_xxxyzz[i] = 2.0 * tr_y_xzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xxzz_xxyzz[i] * fe_0 + tr_y_xxzz_xxxyzz[i] * pa_x[i];

        tr_y_xxxzz_xxxzzz[i] = 2.0 * tr_y_xzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_xxzz_xxzzz[i] * fe_0 + tr_y_xxzz_xxxzzz[i] * pa_x[i];

        tr_y_xxxzz_xxyyyy[i] = tr_y_xxx_xxyyyy[i] * fe_0 + tr_y_xxxz_xxyyyy[i] * pa_z[i];

        tr_y_xxxzz_xxyyyz[i] = 2.0 * tr_y_xzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xxzz_xyyyz[i] * fe_0 + tr_y_xxzz_xxyyyz[i] * pa_x[i];

        tr_y_xxxzz_xxyyzz[i] = 2.0 * tr_y_xzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xxzz_xyyzz[i] * fe_0 + tr_y_xxzz_xxyyzz[i] * pa_x[i];

        tr_y_xxxzz_xxyzzz[i] = 2.0 * tr_y_xzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xxzz_xyzzz[i] * fe_0 + tr_y_xxzz_xxyzzz[i] * pa_x[i];

        tr_y_xxxzz_xxzzzz[i] = 2.0 * tr_y_xzz_xxzzzz[i] * fe_0 + 2.0 * tr_y_xxzz_xzzzz[i] * fe_0 + tr_y_xxzz_xxzzzz[i] * pa_x[i];

        tr_y_xxxzz_xyyyyy[i] = tr_y_xxx_xyyyyy[i] * fe_0 + tr_y_xxxz_xyyyyy[i] * pa_z[i];

        tr_y_xxxzz_xyyyyz[i] = 2.0 * tr_y_xzz_xyyyyz[i] * fe_0 + tr_y_xxzz_yyyyz[i] * fe_0 + tr_y_xxzz_xyyyyz[i] * pa_x[i];

        tr_y_xxxzz_xyyyzz[i] = 2.0 * tr_y_xzz_xyyyzz[i] * fe_0 + tr_y_xxzz_yyyzz[i] * fe_0 + tr_y_xxzz_xyyyzz[i] * pa_x[i];

        tr_y_xxxzz_xyyzzz[i] = 2.0 * tr_y_xzz_xyyzzz[i] * fe_0 + tr_y_xxzz_yyzzz[i] * fe_0 + tr_y_xxzz_xyyzzz[i] * pa_x[i];

        tr_y_xxxzz_xyzzzz[i] = 2.0 * tr_y_xzz_xyzzzz[i] * fe_0 + tr_y_xxzz_yzzzz[i] * fe_0 + tr_y_xxzz_xyzzzz[i] * pa_x[i];

        tr_y_xxxzz_xzzzzz[i] = 2.0 * tr_y_xzz_xzzzzz[i] * fe_0 + tr_y_xxzz_zzzzz[i] * fe_0 + tr_y_xxzz_xzzzzz[i] * pa_x[i];

        tr_y_xxxzz_yyyyyy[i] = 2.0 * tr_y_xzz_yyyyyy[i] * fe_0 + tr_y_xxzz_yyyyyy[i] * pa_x[i];

        tr_y_xxxzz_yyyyyz[i] = 2.0 * tr_y_xzz_yyyyyz[i] * fe_0 + tr_y_xxzz_yyyyyz[i] * pa_x[i];

        tr_y_xxxzz_yyyyzz[i] = 2.0 * tr_y_xzz_yyyyzz[i] * fe_0 + tr_y_xxzz_yyyyzz[i] * pa_x[i];

        tr_y_xxxzz_yyyzzz[i] = 2.0 * tr_y_xzz_yyyzzz[i] * fe_0 + tr_y_xxzz_yyyzzz[i] * pa_x[i];

        tr_y_xxxzz_yyzzzz[i] = 2.0 * tr_y_xzz_yyzzzz[i] * fe_0 + tr_y_xxzz_yyzzzz[i] * pa_x[i];

        tr_y_xxxzz_yzzzzz[i] = 2.0 * tr_y_xzz_yzzzzz[i] * fe_0 + tr_y_xxzz_yzzzzz[i] * pa_x[i];

        tr_y_xxxzz_zzzzzz[i] = 2.0 * tr_y_xzz_zzzzzz[i] * fe_0 + tr_y_xxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 756-784 components of targeted buffer : HI

    auto tr_y_xxyyy_xxxxxx = pbuffer.data(idx_dip_hi + 756);

    auto tr_y_xxyyy_xxxxxy = pbuffer.data(idx_dip_hi + 757);

    auto tr_y_xxyyy_xxxxxz = pbuffer.data(idx_dip_hi + 758);

    auto tr_y_xxyyy_xxxxyy = pbuffer.data(idx_dip_hi + 759);

    auto tr_y_xxyyy_xxxxyz = pbuffer.data(idx_dip_hi + 760);

    auto tr_y_xxyyy_xxxxzz = pbuffer.data(idx_dip_hi + 761);

    auto tr_y_xxyyy_xxxyyy = pbuffer.data(idx_dip_hi + 762);

    auto tr_y_xxyyy_xxxyyz = pbuffer.data(idx_dip_hi + 763);

    auto tr_y_xxyyy_xxxyzz = pbuffer.data(idx_dip_hi + 764);

    auto tr_y_xxyyy_xxxzzz = pbuffer.data(idx_dip_hi + 765);

    auto tr_y_xxyyy_xxyyyy = pbuffer.data(idx_dip_hi + 766);

    auto tr_y_xxyyy_xxyyyz = pbuffer.data(idx_dip_hi + 767);

    auto tr_y_xxyyy_xxyyzz = pbuffer.data(idx_dip_hi + 768);

    auto tr_y_xxyyy_xxyzzz = pbuffer.data(idx_dip_hi + 769);

    auto tr_y_xxyyy_xxzzzz = pbuffer.data(idx_dip_hi + 770);

    auto tr_y_xxyyy_xyyyyy = pbuffer.data(idx_dip_hi + 771);

    auto tr_y_xxyyy_xyyyyz = pbuffer.data(idx_dip_hi + 772);

    auto tr_y_xxyyy_xyyyzz = pbuffer.data(idx_dip_hi + 773);

    auto tr_y_xxyyy_xyyzzz = pbuffer.data(idx_dip_hi + 774);

    auto tr_y_xxyyy_xyzzzz = pbuffer.data(idx_dip_hi + 775);

    auto tr_y_xxyyy_xzzzzz = pbuffer.data(idx_dip_hi + 776);

    auto tr_y_xxyyy_yyyyyy = pbuffer.data(idx_dip_hi + 777);

    auto tr_y_xxyyy_yyyyyz = pbuffer.data(idx_dip_hi + 778);

    auto tr_y_xxyyy_yyyyzz = pbuffer.data(idx_dip_hi + 779);

    auto tr_y_xxyyy_yyyzzz = pbuffer.data(idx_dip_hi + 780);

    auto tr_y_xxyyy_yyzzzz = pbuffer.data(idx_dip_hi + 781);

    auto tr_y_xxyyy_yzzzzz = pbuffer.data(idx_dip_hi + 782);

    auto tr_y_xxyyy_zzzzzz = pbuffer.data(idx_dip_hi + 783);

#pragma omp simd aligned(pa_x,                  \
                             tr_y_xxyyy_xxxxxx, \
                             tr_y_xxyyy_xxxxxy, \
                             tr_y_xxyyy_xxxxxz, \
                             tr_y_xxyyy_xxxxyy, \
                             tr_y_xxyyy_xxxxyz, \
                             tr_y_xxyyy_xxxxzz, \
                             tr_y_xxyyy_xxxyyy, \
                             tr_y_xxyyy_xxxyyz, \
                             tr_y_xxyyy_xxxyzz, \
                             tr_y_xxyyy_xxxzzz, \
                             tr_y_xxyyy_xxyyyy, \
                             tr_y_xxyyy_xxyyyz, \
                             tr_y_xxyyy_xxyyzz, \
                             tr_y_xxyyy_xxyzzz, \
                             tr_y_xxyyy_xxzzzz, \
                             tr_y_xxyyy_xyyyyy, \
                             tr_y_xxyyy_xyyyyz, \
                             tr_y_xxyyy_xyyyzz, \
                             tr_y_xxyyy_xyyzzz, \
                             tr_y_xxyyy_xyzzzz, \
                             tr_y_xxyyy_xzzzzz, \
                             tr_y_xxyyy_yyyyyy, \
                             tr_y_xxyyy_yyyyyz, \
                             tr_y_xxyyy_yyyyzz, \
                             tr_y_xxyyy_yyyzzz, \
                             tr_y_xxyyy_yyzzzz, \
                             tr_y_xxyyy_yzzzzz, \
                             tr_y_xxyyy_zzzzzz, \
                             tr_y_xyyy_xxxxx,   \
                             tr_y_xyyy_xxxxxx,  \
                             tr_y_xyyy_xxxxxy,  \
                             tr_y_xyyy_xxxxxz,  \
                             tr_y_xyyy_xxxxy,   \
                             tr_y_xyyy_xxxxyy,  \
                             tr_y_xyyy_xxxxyz,  \
                             tr_y_xyyy_xxxxz,   \
                             tr_y_xyyy_xxxxzz,  \
                             tr_y_xyyy_xxxyy,   \
                             tr_y_xyyy_xxxyyy,  \
                             tr_y_xyyy_xxxyyz,  \
                             tr_y_xyyy_xxxyz,   \
                             tr_y_xyyy_xxxyzz,  \
                             tr_y_xyyy_xxxzz,   \
                             tr_y_xyyy_xxxzzz,  \
                             tr_y_xyyy_xxyyy,   \
                             tr_y_xyyy_xxyyyy,  \
                             tr_y_xyyy_xxyyyz,  \
                             tr_y_xyyy_xxyyz,   \
                             tr_y_xyyy_xxyyzz,  \
                             tr_y_xyyy_xxyzz,   \
                             tr_y_xyyy_xxyzzz,  \
                             tr_y_xyyy_xxzzz,   \
                             tr_y_xyyy_xxzzzz,  \
                             tr_y_xyyy_xyyyy,   \
                             tr_y_xyyy_xyyyyy,  \
                             tr_y_xyyy_xyyyyz,  \
                             tr_y_xyyy_xyyyz,   \
                             tr_y_xyyy_xyyyzz,  \
                             tr_y_xyyy_xyyzz,   \
                             tr_y_xyyy_xyyzzz,  \
                             tr_y_xyyy_xyzzz,   \
                             tr_y_xyyy_xyzzzz,  \
                             tr_y_xyyy_xzzzz,   \
                             tr_y_xyyy_xzzzzz,  \
                             tr_y_xyyy_yyyyy,   \
                             tr_y_xyyy_yyyyyy,  \
                             tr_y_xyyy_yyyyyz,  \
                             tr_y_xyyy_yyyyz,   \
                             tr_y_xyyy_yyyyzz,  \
                             tr_y_xyyy_yyyzz,   \
                             tr_y_xyyy_yyyzzz,  \
                             tr_y_xyyy_yyzzz,   \
                             tr_y_xyyy_yyzzzz,  \
                             tr_y_xyyy_yzzzz,   \
                             tr_y_xyyy_yzzzzz,  \
                             tr_y_xyyy_zzzzz,   \
                             tr_y_xyyy_zzzzzz,  \
                             tr_y_yyy_xxxxxx,   \
                             tr_y_yyy_xxxxxy,   \
                             tr_y_yyy_xxxxxz,   \
                             tr_y_yyy_xxxxyy,   \
                             tr_y_yyy_xxxxyz,   \
                             tr_y_yyy_xxxxzz,   \
                             tr_y_yyy_xxxyyy,   \
                             tr_y_yyy_xxxyyz,   \
                             tr_y_yyy_xxxyzz,   \
                             tr_y_yyy_xxxzzz,   \
                             tr_y_yyy_xxyyyy,   \
                             tr_y_yyy_xxyyyz,   \
                             tr_y_yyy_xxyyzz,   \
                             tr_y_yyy_xxyzzz,   \
                             tr_y_yyy_xxzzzz,   \
                             tr_y_yyy_xyyyyy,   \
                             tr_y_yyy_xyyyyz,   \
                             tr_y_yyy_xyyyzz,   \
                             tr_y_yyy_xyyzzz,   \
                             tr_y_yyy_xyzzzz,   \
                             tr_y_yyy_xzzzzz,   \
                             tr_y_yyy_yyyyyy,   \
                             tr_y_yyy_yyyyyz,   \
                             tr_y_yyy_yyyyzz,   \
                             tr_y_yyy_yyyzzz,   \
                             tr_y_yyy_yyzzzz,   \
                             tr_y_yyy_yzzzzz,   \
                             tr_y_yyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyy_xxxxxx[i] = tr_y_yyy_xxxxxx[i] * fe_0 + 6.0 * tr_y_xyyy_xxxxx[i] * fe_0 + tr_y_xyyy_xxxxxx[i] * pa_x[i];

        tr_y_xxyyy_xxxxxy[i] = tr_y_yyy_xxxxxy[i] * fe_0 + 5.0 * tr_y_xyyy_xxxxy[i] * fe_0 + tr_y_xyyy_xxxxxy[i] * pa_x[i];

        tr_y_xxyyy_xxxxxz[i] = tr_y_yyy_xxxxxz[i] * fe_0 + 5.0 * tr_y_xyyy_xxxxz[i] * fe_0 + tr_y_xyyy_xxxxxz[i] * pa_x[i];

        tr_y_xxyyy_xxxxyy[i] = tr_y_yyy_xxxxyy[i] * fe_0 + 4.0 * tr_y_xyyy_xxxyy[i] * fe_0 + tr_y_xyyy_xxxxyy[i] * pa_x[i];

        tr_y_xxyyy_xxxxyz[i] = tr_y_yyy_xxxxyz[i] * fe_0 + 4.0 * tr_y_xyyy_xxxyz[i] * fe_0 + tr_y_xyyy_xxxxyz[i] * pa_x[i];

        tr_y_xxyyy_xxxxzz[i] = tr_y_yyy_xxxxzz[i] * fe_0 + 4.0 * tr_y_xyyy_xxxzz[i] * fe_0 + tr_y_xyyy_xxxxzz[i] * pa_x[i];

        tr_y_xxyyy_xxxyyy[i] = tr_y_yyy_xxxyyy[i] * fe_0 + 3.0 * tr_y_xyyy_xxyyy[i] * fe_0 + tr_y_xyyy_xxxyyy[i] * pa_x[i];

        tr_y_xxyyy_xxxyyz[i] = tr_y_yyy_xxxyyz[i] * fe_0 + 3.0 * tr_y_xyyy_xxyyz[i] * fe_0 + tr_y_xyyy_xxxyyz[i] * pa_x[i];

        tr_y_xxyyy_xxxyzz[i] = tr_y_yyy_xxxyzz[i] * fe_0 + 3.0 * tr_y_xyyy_xxyzz[i] * fe_0 + tr_y_xyyy_xxxyzz[i] * pa_x[i];

        tr_y_xxyyy_xxxzzz[i] = tr_y_yyy_xxxzzz[i] * fe_0 + 3.0 * tr_y_xyyy_xxzzz[i] * fe_0 + tr_y_xyyy_xxxzzz[i] * pa_x[i];

        tr_y_xxyyy_xxyyyy[i] = tr_y_yyy_xxyyyy[i] * fe_0 + 2.0 * tr_y_xyyy_xyyyy[i] * fe_0 + tr_y_xyyy_xxyyyy[i] * pa_x[i];

        tr_y_xxyyy_xxyyyz[i] = tr_y_yyy_xxyyyz[i] * fe_0 + 2.0 * tr_y_xyyy_xyyyz[i] * fe_0 + tr_y_xyyy_xxyyyz[i] * pa_x[i];

        tr_y_xxyyy_xxyyzz[i] = tr_y_yyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_xyyy_xyyzz[i] * fe_0 + tr_y_xyyy_xxyyzz[i] * pa_x[i];

        tr_y_xxyyy_xxyzzz[i] = tr_y_yyy_xxyzzz[i] * fe_0 + 2.0 * tr_y_xyyy_xyzzz[i] * fe_0 + tr_y_xyyy_xxyzzz[i] * pa_x[i];

        tr_y_xxyyy_xxzzzz[i] = tr_y_yyy_xxzzzz[i] * fe_0 + 2.0 * tr_y_xyyy_xzzzz[i] * fe_0 + tr_y_xyyy_xxzzzz[i] * pa_x[i];

        tr_y_xxyyy_xyyyyy[i] = tr_y_yyy_xyyyyy[i] * fe_0 + tr_y_xyyy_yyyyy[i] * fe_0 + tr_y_xyyy_xyyyyy[i] * pa_x[i];

        tr_y_xxyyy_xyyyyz[i] = tr_y_yyy_xyyyyz[i] * fe_0 + tr_y_xyyy_yyyyz[i] * fe_0 + tr_y_xyyy_xyyyyz[i] * pa_x[i];

        tr_y_xxyyy_xyyyzz[i] = tr_y_yyy_xyyyzz[i] * fe_0 + tr_y_xyyy_yyyzz[i] * fe_0 + tr_y_xyyy_xyyyzz[i] * pa_x[i];

        tr_y_xxyyy_xyyzzz[i] = tr_y_yyy_xyyzzz[i] * fe_0 + tr_y_xyyy_yyzzz[i] * fe_0 + tr_y_xyyy_xyyzzz[i] * pa_x[i];

        tr_y_xxyyy_xyzzzz[i] = tr_y_yyy_xyzzzz[i] * fe_0 + tr_y_xyyy_yzzzz[i] * fe_0 + tr_y_xyyy_xyzzzz[i] * pa_x[i];

        tr_y_xxyyy_xzzzzz[i] = tr_y_yyy_xzzzzz[i] * fe_0 + tr_y_xyyy_zzzzz[i] * fe_0 + tr_y_xyyy_xzzzzz[i] * pa_x[i];

        tr_y_xxyyy_yyyyyy[i] = tr_y_yyy_yyyyyy[i] * fe_0 + tr_y_xyyy_yyyyyy[i] * pa_x[i];

        tr_y_xxyyy_yyyyyz[i] = tr_y_yyy_yyyyyz[i] * fe_0 + tr_y_xyyy_yyyyyz[i] * pa_x[i];

        tr_y_xxyyy_yyyyzz[i] = tr_y_yyy_yyyyzz[i] * fe_0 + tr_y_xyyy_yyyyzz[i] * pa_x[i];

        tr_y_xxyyy_yyyzzz[i] = tr_y_yyy_yyyzzz[i] * fe_0 + tr_y_xyyy_yyyzzz[i] * pa_x[i];

        tr_y_xxyyy_yyzzzz[i] = tr_y_yyy_yyzzzz[i] * fe_0 + tr_y_xyyy_yyzzzz[i] * pa_x[i];

        tr_y_xxyyy_yzzzzz[i] = tr_y_yyy_yzzzzz[i] * fe_0 + tr_y_xyyy_yzzzzz[i] * pa_x[i];

        tr_y_xxyyy_zzzzzz[i] = tr_y_yyy_zzzzzz[i] * fe_0 + tr_y_xyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 784-812 components of targeted buffer : HI

    auto tr_y_xxyyz_xxxxxx = pbuffer.data(idx_dip_hi + 784);

    auto tr_y_xxyyz_xxxxxy = pbuffer.data(idx_dip_hi + 785);

    auto tr_y_xxyyz_xxxxxz = pbuffer.data(idx_dip_hi + 786);

    auto tr_y_xxyyz_xxxxyy = pbuffer.data(idx_dip_hi + 787);

    auto tr_y_xxyyz_xxxxyz = pbuffer.data(idx_dip_hi + 788);

    auto tr_y_xxyyz_xxxxzz = pbuffer.data(idx_dip_hi + 789);

    auto tr_y_xxyyz_xxxyyy = pbuffer.data(idx_dip_hi + 790);

    auto tr_y_xxyyz_xxxyyz = pbuffer.data(idx_dip_hi + 791);

    auto tr_y_xxyyz_xxxyzz = pbuffer.data(idx_dip_hi + 792);

    auto tr_y_xxyyz_xxxzzz = pbuffer.data(idx_dip_hi + 793);

    auto tr_y_xxyyz_xxyyyy = pbuffer.data(idx_dip_hi + 794);

    auto tr_y_xxyyz_xxyyyz = pbuffer.data(idx_dip_hi + 795);

    auto tr_y_xxyyz_xxyyzz = pbuffer.data(idx_dip_hi + 796);

    auto tr_y_xxyyz_xxyzzz = pbuffer.data(idx_dip_hi + 797);

    auto tr_y_xxyyz_xxzzzz = pbuffer.data(idx_dip_hi + 798);

    auto tr_y_xxyyz_xyyyyy = pbuffer.data(idx_dip_hi + 799);

    auto tr_y_xxyyz_xyyyyz = pbuffer.data(idx_dip_hi + 800);

    auto tr_y_xxyyz_xyyyzz = pbuffer.data(idx_dip_hi + 801);

    auto tr_y_xxyyz_xyyzzz = pbuffer.data(idx_dip_hi + 802);

    auto tr_y_xxyyz_xyzzzz = pbuffer.data(idx_dip_hi + 803);

    auto tr_y_xxyyz_xzzzzz = pbuffer.data(idx_dip_hi + 804);

    auto tr_y_xxyyz_yyyyyy = pbuffer.data(idx_dip_hi + 805);

    auto tr_y_xxyyz_yyyyyz = pbuffer.data(idx_dip_hi + 806);

    auto tr_y_xxyyz_yyyyzz = pbuffer.data(idx_dip_hi + 807);

    auto tr_y_xxyyz_yyyzzz = pbuffer.data(idx_dip_hi + 808);

    auto tr_y_xxyyz_yyzzzz = pbuffer.data(idx_dip_hi + 809);

    auto tr_y_xxyyz_yzzzzz = pbuffer.data(idx_dip_hi + 810);

    auto tr_y_xxyyz_zzzzzz = pbuffer.data(idx_dip_hi + 811);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_y_xxyy_xxxxx,   \
                             tr_y_xxyy_xxxxxx,  \
                             tr_y_xxyy_xxxxxy,  \
                             tr_y_xxyy_xxxxxz,  \
                             tr_y_xxyy_xxxxy,   \
                             tr_y_xxyy_xxxxyy,  \
                             tr_y_xxyy_xxxxyz,  \
                             tr_y_xxyy_xxxxz,   \
                             tr_y_xxyy_xxxxzz,  \
                             tr_y_xxyy_xxxyy,   \
                             tr_y_xxyy_xxxyyy,  \
                             tr_y_xxyy_xxxyyz,  \
                             tr_y_xxyy_xxxyz,   \
                             tr_y_xxyy_xxxyzz,  \
                             tr_y_xxyy_xxxzz,   \
                             tr_y_xxyy_xxxzzz,  \
                             tr_y_xxyy_xxyyy,   \
                             tr_y_xxyy_xxyyyy,  \
                             tr_y_xxyy_xxyyyz,  \
                             tr_y_xxyy_xxyyz,   \
                             tr_y_xxyy_xxyyzz,  \
                             tr_y_xxyy_xxyzz,   \
                             tr_y_xxyy_xxyzzz,  \
                             tr_y_xxyy_xxzzz,   \
                             tr_y_xxyy_xxzzzz,  \
                             tr_y_xxyy_xyyyy,   \
                             tr_y_xxyy_xyyyyy,  \
                             tr_y_xxyy_xyyyyz,  \
                             tr_y_xxyy_xyyyz,   \
                             tr_y_xxyy_xyyyzz,  \
                             tr_y_xxyy_xyyzz,   \
                             tr_y_xxyy_xyyzzz,  \
                             tr_y_xxyy_xyzzz,   \
                             tr_y_xxyy_xyzzzz,  \
                             tr_y_xxyy_xzzzz,   \
                             tr_y_xxyy_xzzzzz,  \
                             tr_y_xxyy_yyyyyy,  \
                             tr_y_xxyyz_xxxxxx, \
                             tr_y_xxyyz_xxxxxy, \
                             tr_y_xxyyz_xxxxxz, \
                             tr_y_xxyyz_xxxxyy, \
                             tr_y_xxyyz_xxxxyz, \
                             tr_y_xxyyz_xxxxzz, \
                             tr_y_xxyyz_xxxyyy, \
                             tr_y_xxyyz_xxxyyz, \
                             tr_y_xxyyz_xxxyzz, \
                             tr_y_xxyyz_xxxzzz, \
                             tr_y_xxyyz_xxyyyy, \
                             tr_y_xxyyz_xxyyyz, \
                             tr_y_xxyyz_xxyyzz, \
                             tr_y_xxyyz_xxyzzz, \
                             tr_y_xxyyz_xxzzzz, \
                             tr_y_xxyyz_xyyyyy, \
                             tr_y_xxyyz_xyyyyz, \
                             tr_y_xxyyz_xyyyzz, \
                             tr_y_xxyyz_xyyzzz, \
                             tr_y_xxyyz_xyzzzz, \
                             tr_y_xxyyz_xzzzzz, \
                             tr_y_xxyyz_yyyyyy, \
                             tr_y_xxyyz_yyyyyz, \
                             tr_y_xxyyz_yyyyzz, \
                             tr_y_xxyyz_yyyzzz, \
                             tr_y_xxyyz_yyzzzz, \
                             tr_y_xxyyz_yzzzzz, \
                             tr_y_xxyyz_zzzzzz, \
                             tr_y_xyyz_yyyyyz,  \
                             tr_y_xyyz_yyyyzz,  \
                             tr_y_xyyz_yyyzzz,  \
                             tr_y_xyyz_yyzzzz,  \
                             tr_y_xyyz_yzzzzz,  \
                             tr_y_xyyz_zzzzzz,  \
                             tr_y_yyz_yyyyyz,   \
                             tr_y_yyz_yyyyzz,   \
                             tr_y_yyz_yyyzzz,   \
                             tr_y_yyz_yyzzzz,   \
                             tr_y_yyz_yzzzzz,   \
                             tr_y_yyz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyz_xxxxxx[i] = tr_y_xxyy_xxxxxx[i] * pa_z[i];

        tr_y_xxyyz_xxxxxy[i] = tr_y_xxyy_xxxxxy[i] * pa_z[i];

        tr_y_xxyyz_xxxxxz[i] = tr_y_xxyy_xxxxx[i] * fe_0 + tr_y_xxyy_xxxxxz[i] * pa_z[i];

        tr_y_xxyyz_xxxxyy[i] = tr_y_xxyy_xxxxyy[i] * pa_z[i];

        tr_y_xxyyz_xxxxyz[i] = tr_y_xxyy_xxxxy[i] * fe_0 + tr_y_xxyy_xxxxyz[i] * pa_z[i];

        tr_y_xxyyz_xxxxzz[i] = 2.0 * tr_y_xxyy_xxxxz[i] * fe_0 + tr_y_xxyy_xxxxzz[i] * pa_z[i];

        tr_y_xxyyz_xxxyyy[i] = tr_y_xxyy_xxxyyy[i] * pa_z[i];

        tr_y_xxyyz_xxxyyz[i] = tr_y_xxyy_xxxyy[i] * fe_0 + tr_y_xxyy_xxxyyz[i] * pa_z[i];

        tr_y_xxyyz_xxxyzz[i] = 2.0 * tr_y_xxyy_xxxyz[i] * fe_0 + tr_y_xxyy_xxxyzz[i] * pa_z[i];

        tr_y_xxyyz_xxxzzz[i] = 3.0 * tr_y_xxyy_xxxzz[i] * fe_0 + tr_y_xxyy_xxxzzz[i] * pa_z[i];

        tr_y_xxyyz_xxyyyy[i] = tr_y_xxyy_xxyyyy[i] * pa_z[i];

        tr_y_xxyyz_xxyyyz[i] = tr_y_xxyy_xxyyy[i] * fe_0 + tr_y_xxyy_xxyyyz[i] * pa_z[i];

        tr_y_xxyyz_xxyyzz[i] = 2.0 * tr_y_xxyy_xxyyz[i] * fe_0 + tr_y_xxyy_xxyyzz[i] * pa_z[i];

        tr_y_xxyyz_xxyzzz[i] = 3.0 * tr_y_xxyy_xxyzz[i] * fe_0 + tr_y_xxyy_xxyzzz[i] * pa_z[i];

        tr_y_xxyyz_xxzzzz[i] = 4.0 * tr_y_xxyy_xxzzz[i] * fe_0 + tr_y_xxyy_xxzzzz[i] * pa_z[i];

        tr_y_xxyyz_xyyyyy[i] = tr_y_xxyy_xyyyyy[i] * pa_z[i];

        tr_y_xxyyz_xyyyyz[i] = tr_y_xxyy_xyyyy[i] * fe_0 + tr_y_xxyy_xyyyyz[i] * pa_z[i];

        tr_y_xxyyz_xyyyzz[i] = 2.0 * tr_y_xxyy_xyyyz[i] * fe_0 + tr_y_xxyy_xyyyzz[i] * pa_z[i];

        tr_y_xxyyz_xyyzzz[i] = 3.0 * tr_y_xxyy_xyyzz[i] * fe_0 + tr_y_xxyy_xyyzzz[i] * pa_z[i];

        tr_y_xxyyz_xyzzzz[i] = 4.0 * tr_y_xxyy_xyzzz[i] * fe_0 + tr_y_xxyy_xyzzzz[i] * pa_z[i];

        tr_y_xxyyz_xzzzzz[i] = 5.0 * tr_y_xxyy_xzzzz[i] * fe_0 + tr_y_xxyy_xzzzzz[i] * pa_z[i];

        tr_y_xxyyz_yyyyyy[i] = tr_y_xxyy_yyyyyy[i] * pa_z[i];

        tr_y_xxyyz_yyyyyz[i] = tr_y_yyz_yyyyyz[i] * fe_0 + tr_y_xyyz_yyyyyz[i] * pa_x[i];

        tr_y_xxyyz_yyyyzz[i] = tr_y_yyz_yyyyzz[i] * fe_0 + tr_y_xyyz_yyyyzz[i] * pa_x[i];

        tr_y_xxyyz_yyyzzz[i] = tr_y_yyz_yyyzzz[i] * fe_0 + tr_y_xyyz_yyyzzz[i] * pa_x[i];

        tr_y_xxyyz_yyzzzz[i] = tr_y_yyz_yyzzzz[i] * fe_0 + tr_y_xyyz_yyzzzz[i] * pa_x[i];

        tr_y_xxyyz_yzzzzz[i] = tr_y_yyz_yzzzzz[i] * fe_0 + tr_y_xyyz_yzzzzz[i] * pa_x[i];

        tr_y_xxyyz_zzzzzz[i] = tr_y_yyz_zzzzzz[i] * fe_0 + tr_y_xyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 812-840 components of targeted buffer : HI

    auto tr_y_xxyzz_xxxxxx = pbuffer.data(idx_dip_hi + 812);

    auto tr_y_xxyzz_xxxxxy = pbuffer.data(idx_dip_hi + 813);

    auto tr_y_xxyzz_xxxxxz = pbuffer.data(idx_dip_hi + 814);

    auto tr_y_xxyzz_xxxxyy = pbuffer.data(idx_dip_hi + 815);

    auto tr_y_xxyzz_xxxxyz = pbuffer.data(idx_dip_hi + 816);

    auto tr_y_xxyzz_xxxxzz = pbuffer.data(idx_dip_hi + 817);

    auto tr_y_xxyzz_xxxyyy = pbuffer.data(idx_dip_hi + 818);

    auto tr_y_xxyzz_xxxyyz = pbuffer.data(idx_dip_hi + 819);

    auto tr_y_xxyzz_xxxyzz = pbuffer.data(idx_dip_hi + 820);

    auto tr_y_xxyzz_xxxzzz = pbuffer.data(idx_dip_hi + 821);

    auto tr_y_xxyzz_xxyyyy = pbuffer.data(idx_dip_hi + 822);

    auto tr_y_xxyzz_xxyyyz = pbuffer.data(idx_dip_hi + 823);

    auto tr_y_xxyzz_xxyyzz = pbuffer.data(idx_dip_hi + 824);

    auto tr_y_xxyzz_xxyzzz = pbuffer.data(idx_dip_hi + 825);

    auto tr_y_xxyzz_xxzzzz = pbuffer.data(idx_dip_hi + 826);

    auto tr_y_xxyzz_xyyyyy = pbuffer.data(idx_dip_hi + 827);

    auto tr_y_xxyzz_xyyyyz = pbuffer.data(idx_dip_hi + 828);

    auto tr_y_xxyzz_xyyyzz = pbuffer.data(idx_dip_hi + 829);

    auto tr_y_xxyzz_xyyzzz = pbuffer.data(idx_dip_hi + 830);

    auto tr_y_xxyzz_xyzzzz = pbuffer.data(idx_dip_hi + 831);

    auto tr_y_xxyzz_xzzzzz = pbuffer.data(idx_dip_hi + 832);

    auto tr_y_xxyzz_yyyyyy = pbuffer.data(idx_dip_hi + 833);

    auto tr_y_xxyzz_yyyyyz = pbuffer.data(idx_dip_hi + 834);

    auto tr_y_xxyzz_yyyyzz = pbuffer.data(idx_dip_hi + 835);

    auto tr_y_xxyzz_yyyzzz = pbuffer.data(idx_dip_hi + 836);

    auto tr_y_xxyzz_yyzzzz = pbuffer.data(idx_dip_hi + 837);

    auto tr_y_xxyzz_yzzzzz = pbuffer.data(idx_dip_hi + 838);

    auto tr_y_xxyzz_zzzzzz = pbuffer.data(idx_dip_hi + 839);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             tr_y_xxy_xxxxxy,   \
                             tr_y_xxy_xxxxyy,   \
                             tr_y_xxy_xxxyyy,   \
                             tr_y_xxy_xxyyyy,   \
                             tr_y_xxy_xyyyyy,   \
                             tr_y_xxyz_xxxxxy,  \
                             tr_y_xxyz_xxxxyy,  \
                             tr_y_xxyz_xxxyyy,  \
                             tr_y_xxyz_xxyyyy,  \
                             tr_y_xxyz_xyyyyy,  \
                             tr_y_xxyzz_xxxxxx, \
                             tr_y_xxyzz_xxxxxy, \
                             tr_y_xxyzz_xxxxxz, \
                             tr_y_xxyzz_xxxxyy, \
                             tr_y_xxyzz_xxxxyz, \
                             tr_y_xxyzz_xxxxzz, \
                             tr_y_xxyzz_xxxyyy, \
                             tr_y_xxyzz_xxxyyz, \
                             tr_y_xxyzz_xxxyzz, \
                             tr_y_xxyzz_xxxzzz, \
                             tr_y_xxyzz_xxyyyy, \
                             tr_y_xxyzz_xxyyyz, \
                             tr_y_xxyzz_xxyyzz, \
                             tr_y_xxyzz_xxyzzz, \
                             tr_y_xxyzz_xxzzzz, \
                             tr_y_xxyzz_xyyyyy, \
                             tr_y_xxyzz_xyyyyz, \
                             tr_y_xxyzz_xyyyzz, \
                             tr_y_xxyzz_xyyzzz, \
                             tr_y_xxyzz_xyzzzz, \
                             tr_y_xxyzz_xzzzzz, \
                             tr_y_xxyzz_yyyyyy, \
                             tr_y_xxyzz_yyyyyz, \
                             tr_y_xxyzz_yyyyzz, \
                             tr_y_xxyzz_yyyzzz, \
                             tr_y_xxyzz_yyzzzz, \
                             tr_y_xxyzz_yzzzzz, \
                             tr_y_xxyzz_zzzzzz, \
                             tr_y_xxzz_xxxxxx,  \
                             tr_y_xxzz_xxxxxz,  \
                             tr_y_xxzz_xxxxzz,  \
                             tr_y_xxzz_xxxzzz,  \
                             tr_y_xxzz_xxzzzz,  \
                             tr_y_xxzz_xzzzzz,  \
                             tr_y_xyzz_xxxxyz,  \
                             tr_y_xyzz_xxxyyz,  \
                             tr_y_xyzz_xxxyz,   \
                             tr_y_xyzz_xxxyzz,  \
                             tr_y_xyzz_xxyyyz,  \
                             tr_y_xyzz_xxyyz,   \
                             tr_y_xyzz_xxyyzz,  \
                             tr_y_xyzz_xxyzz,   \
                             tr_y_xyzz_xxyzzz,  \
                             tr_y_xyzz_xyyyyz,  \
                             tr_y_xyzz_xyyyz,   \
                             tr_y_xyzz_xyyyzz,  \
                             tr_y_xyzz_xyyzz,   \
                             tr_y_xyzz_xyyzzz,  \
                             tr_y_xyzz_xyzzz,   \
                             tr_y_xyzz_xyzzzz,  \
                             tr_y_xyzz_yyyyyy,  \
                             tr_y_xyzz_yyyyyz,  \
                             tr_y_xyzz_yyyyz,   \
                             tr_y_xyzz_yyyyzz,  \
                             tr_y_xyzz_yyyzz,   \
                             tr_y_xyzz_yyyzzz,  \
                             tr_y_xyzz_yyzzz,   \
                             tr_y_xyzz_yyzzzz,  \
                             tr_y_xyzz_yzzzz,   \
                             tr_y_xyzz_yzzzzz,  \
                             tr_y_xyzz_zzzzzz,  \
                             tr_y_yzz_xxxxyz,   \
                             tr_y_yzz_xxxyyz,   \
                             tr_y_yzz_xxxyzz,   \
                             tr_y_yzz_xxyyyz,   \
                             tr_y_yzz_xxyyzz,   \
                             tr_y_yzz_xxyzzz,   \
                             tr_y_yzz_xyyyyz,   \
                             tr_y_yzz_xyyyzz,   \
                             tr_y_yzz_xyyzzz,   \
                             tr_y_yzz_xyzzzz,   \
                             tr_y_yzz_yyyyyy,   \
                             tr_y_yzz_yyyyyz,   \
                             tr_y_yzz_yyyyzz,   \
                             tr_y_yzz_yyyzzz,   \
                             tr_y_yzz_yyzzzz,   \
                             tr_y_yzz_yzzzzz,   \
                             tr_y_yzz_zzzzzz,   \
                             ts_xxzz_xxxxxx,    \
                             ts_xxzz_xxxxxz,    \
                             ts_xxzz_xxxxzz,    \
                             ts_xxzz_xxxzzz,    \
                             ts_xxzz_xxzzzz,    \
                             ts_xxzz_xzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzz_xxxxxx[i] = ts_xxzz_xxxxxx[i] * fe_0 + tr_y_xxzz_xxxxxx[i] * pa_y[i];

        tr_y_xxyzz_xxxxxy[i] = tr_y_xxy_xxxxxy[i] * fe_0 + tr_y_xxyz_xxxxxy[i] * pa_z[i];

        tr_y_xxyzz_xxxxxz[i] = ts_xxzz_xxxxxz[i] * fe_0 + tr_y_xxzz_xxxxxz[i] * pa_y[i];

        tr_y_xxyzz_xxxxyy[i] = tr_y_xxy_xxxxyy[i] * fe_0 + tr_y_xxyz_xxxxyy[i] * pa_z[i];

        tr_y_xxyzz_xxxxyz[i] = tr_y_yzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xyzz_xxxyz[i] * fe_0 + tr_y_xyzz_xxxxyz[i] * pa_x[i];

        tr_y_xxyzz_xxxxzz[i] = ts_xxzz_xxxxzz[i] * fe_0 + tr_y_xxzz_xxxxzz[i] * pa_y[i];

        tr_y_xxyzz_xxxyyy[i] = tr_y_xxy_xxxyyy[i] * fe_0 + tr_y_xxyz_xxxyyy[i] * pa_z[i];

        tr_y_xxyzz_xxxyyz[i] = tr_y_yzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xyzz_xxyyz[i] * fe_0 + tr_y_xyzz_xxxyyz[i] * pa_x[i];

        tr_y_xxyzz_xxxyzz[i] = tr_y_yzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xyzz_xxyzz[i] * fe_0 + tr_y_xyzz_xxxyzz[i] * pa_x[i];

        tr_y_xxyzz_xxxzzz[i] = ts_xxzz_xxxzzz[i] * fe_0 + tr_y_xxzz_xxxzzz[i] * pa_y[i];

        tr_y_xxyzz_xxyyyy[i] = tr_y_xxy_xxyyyy[i] * fe_0 + tr_y_xxyz_xxyyyy[i] * pa_z[i];

        tr_y_xxyzz_xxyyyz[i] = tr_y_yzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xyzz_xyyyz[i] * fe_0 + tr_y_xyzz_xxyyyz[i] * pa_x[i];

        tr_y_xxyzz_xxyyzz[i] = tr_y_yzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xyzz_xyyzz[i] * fe_0 + tr_y_xyzz_xxyyzz[i] * pa_x[i];

        tr_y_xxyzz_xxyzzz[i] = tr_y_yzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xyzz_xyzzz[i] * fe_0 + tr_y_xyzz_xxyzzz[i] * pa_x[i];

        tr_y_xxyzz_xxzzzz[i] = ts_xxzz_xxzzzz[i] * fe_0 + tr_y_xxzz_xxzzzz[i] * pa_y[i];

        tr_y_xxyzz_xyyyyy[i] = tr_y_xxy_xyyyyy[i] * fe_0 + tr_y_xxyz_xyyyyy[i] * pa_z[i];

        tr_y_xxyzz_xyyyyz[i] = tr_y_yzz_xyyyyz[i] * fe_0 + tr_y_xyzz_yyyyz[i] * fe_0 + tr_y_xyzz_xyyyyz[i] * pa_x[i];

        tr_y_xxyzz_xyyyzz[i] = tr_y_yzz_xyyyzz[i] * fe_0 + tr_y_xyzz_yyyzz[i] * fe_0 + tr_y_xyzz_xyyyzz[i] * pa_x[i];

        tr_y_xxyzz_xyyzzz[i] = tr_y_yzz_xyyzzz[i] * fe_0 + tr_y_xyzz_yyzzz[i] * fe_0 + tr_y_xyzz_xyyzzz[i] * pa_x[i];

        tr_y_xxyzz_xyzzzz[i] = tr_y_yzz_xyzzzz[i] * fe_0 + tr_y_xyzz_yzzzz[i] * fe_0 + tr_y_xyzz_xyzzzz[i] * pa_x[i];

        tr_y_xxyzz_xzzzzz[i] = ts_xxzz_xzzzzz[i] * fe_0 + tr_y_xxzz_xzzzzz[i] * pa_y[i];

        tr_y_xxyzz_yyyyyy[i] = tr_y_yzz_yyyyyy[i] * fe_0 + tr_y_xyzz_yyyyyy[i] * pa_x[i];

        tr_y_xxyzz_yyyyyz[i] = tr_y_yzz_yyyyyz[i] * fe_0 + tr_y_xyzz_yyyyyz[i] * pa_x[i];

        tr_y_xxyzz_yyyyzz[i] = tr_y_yzz_yyyyzz[i] * fe_0 + tr_y_xyzz_yyyyzz[i] * pa_x[i];

        tr_y_xxyzz_yyyzzz[i] = tr_y_yzz_yyyzzz[i] * fe_0 + tr_y_xyzz_yyyzzz[i] * pa_x[i];

        tr_y_xxyzz_yyzzzz[i] = tr_y_yzz_yyzzzz[i] * fe_0 + tr_y_xyzz_yyzzzz[i] * pa_x[i];

        tr_y_xxyzz_yzzzzz[i] = tr_y_yzz_yzzzzz[i] * fe_0 + tr_y_xyzz_yzzzzz[i] * pa_x[i];

        tr_y_xxyzz_zzzzzz[i] = tr_y_yzz_zzzzzz[i] * fe_0 + tr_y_xyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 840-868 components of targeted buffer : HI

    auto tr_y_xxzzz_xxxxxx = pbuffer.data(idx_dip_hi + 840);

    auto tr_y_xxzzz_xxxxxy = pbuffer.data(idx_dip_hi + 841);

    auto tr_y_xxzzz_xxxxxz = pbuffer.data(idx_dip_hi + 842);

    auto tr_y_xxzzz_xxxxyy = pbuffer.data(idx_dip_hi + 843);

    auto tr_y_xxzzz_xxxxyz = pbuffer.data(idx_dip_hi + 844);

    auto tr_y_xxzzz_xxxxzz = pbuffer.data(idx_dip_hi + 845);

    auto tr_y_xxzzz_xxxyyy = pbuffer.data(idx_dip_hi + 846);

    auto tr_y_xxzzz_xxxyyz = pbuffer.data(idx_dip_hi + 847);

    auto tr_y_xxzzz_xxxyzz = pbuffer.data(idx_dip_hi + 848);

    auto tr_y_xxzzz_xxxzzz = pbuffer.data(idx_dip_hi + 849);

    auto tr_y_xxzzz_xxyyyy = pbuffer.data(idx_dip_hi + 850);

    auto tr_y_xxzzz_xxyyyz = pbuffer.data(idx_dip_hi + 851);

    auto tr_y_xxzzz_xxyyzz = pbuffer.data(idx_dip_hi + 852);

    auto tr_y_xxzzz_xxyzzz = pbuffer.data(idx_dip_hi + 853);

    auto tr_y_xxzzz_xxzzzz = pbuffer.data(idx_dip_hi + 854);

    auto tr_y_xxzzz_xyyyyy = pbuffer.data(idx_dip_hi + 855);

    auto tr_y_xxzzz_xyyyyz = pbuffer.data(idx_dip_hi + 856);

    auto tr_y_xxzzz_xyyyzz = pbuffer.data(idx_dip_hi + 857);

    auto tr_y_xxzzz_xyyzzz = pbuffer.data(idx_dip_hi + 858);

    auto tr_y_xxzzz_xyzzzz = pbuffer.data(idx_dip_hi + 859);

    auto tr_y_xxzzz_xzzzzz = pbuffer.data(idx_dip_hi + 860);

    auto tr_y_xxzzz_yyyyyy = pbuffer.data(idx_dip_hi + 861);

    auto tr_y_xxzzz_yyyyyz = pbuffer.data(idx_dip_hi + 862);

    auto tr_y_xxzzz_yyyyzz = pbuffer.data(idx_dip_hi + 863);

    auto tr_y_xxzzz_yyyzzz = pbuffer.data(idx_dip_hi + 864);

    auto tr_y_xxzzz_yyzzzz = pbuffer.data(idx_dip_hi + 865);

    auto tr_y_xxzzz_yzzzzz = pbuffer.data(idx_dip_hi + 866);

    auto tr_y_xxzzz_zzzzzz = pbuffer.data(idx_dip_hi + 867);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_y_xxz_xxxxxx,   \
                             tr_y_xxz_xxxxxy,   \
                             tr_y_xxz_xxxxyy,   \
                             tr_y_xxz_xxxyyy,   \
                             tr_y_xxz_xxyyyy,   \
                             tr_y_xxz_xyyyyy,   \
                             tr_y_xxzz_xxxxxx,  \
                             tr_y_xxzz_xxxxxy,  \
                             tr_y_xxzz_xxxxyy,  \
                             tr_y_xxzz_xxxyyy,  \
                             tr_y_xxzz_xxyyyy,  \
                             tr_y_xxzz_xyyyyy,  \
                             tr_y_xxzzz_xxxxxx, \
                             tr_y_xxzzz_xxxxxy, \
                             tr_y_xxzzz_xxxxxz, \
                             tr_y_xxzzz_xxxxyy, \
                             tr_y_xxzzz_xxxxyz, \
                             tr_y_xxzzz_xxxxzz, \
                             tr_y_xxzzz_xxxyyy, \
                             tr_y_xxzzz_xxxyyz, \
                             tr_y_xxzzz_xxxyzz, \
                             tr_y_xxzzz_xxxzzz, \
                             tr_y_xxzzz_xxyyyy, \
                             tr_y_xxzzz_xxyyyz, \
                             tr_y_xxzzz_xxyyzz, \
                             tr_y_xxzzz_xxyzzz, \
                             tr_y_xxzzz_xxzzzz, \
                             tr_y_xxzzz_xyyyyy, \
                             tr_y_xxzzz_xyyyyz, \
                             tr_y_xxzzz_xyyyzz, \
                             tr_y_xxzzz_xyyzzz, \
                             tr_y_xxzzz_xyzzzz, \
                             tr_y_xxzzz_xzzzzz, \
                             tr_y_xxzzz_yyyyyy, \
                             tr_y_xxzzz_yyyyyz, \
                             tr_y_xxzzz_yyyyzz, \
                             tr_y_xxzzz_yyyzzz, \
                             tr_y_xxzzz_yyzzzz, \
                             tr_y_xxzzz_yzzzzz, \
                             tr_y_xxzzz_zzzzzz, \
                             tr_y_xzzz_xxxxxz,  \
                             tr_y_xzzz_xxxxyz,  \
                             tr_y_xzzz_xxxxz,   \
                             tr_y_xzzz_xxxxzz,  \
                             tr_y_xzzz_xxxyyz,  \
                             tr_y_xzzz_xxxyz,   \
                             tr_y_xzzz_xxxyzz,  \
                             tr_y_xzzz_xxxzz,   \
                             tr_y_xzzz_xxxzzz,  \
                             tr_y_xzzz_xxyyyz,  \
                             tr_y_xzzz_xxyyz,   \
                             tr_y_xzzz_xxyyzz,  \
                             tr_y_xzzz_xxyzz,   \
                             tr_y_xzzz_xxyzzz,  \
                             tr_y_xzzz_xxzzz,   \
                             tr_y_xzzz_xxzzzz,  \
                             tr_y_xzzz_xyyyyz,  \
                             tr_y_xzzz_xyyyz,   \
                             tr_y_xzzz_xyyyzz,  \
                             tr_y_xzzz_xyyzz,   \
                             tr_y_xzzz_xyyzzz,  \
                             tr_y_xzzz_xyzzz,   \
                             tr_y_xzzz_xyzzzz,  \
                             tr_y_xzzz_xzzzz,   \
                             tr_y_xzzz_xzzzzz,  \
                             tr_y_xzzz_yyyyyy,  \
                             tr_y_xzzz_yyyyyz,  \
                             tr_y_xzzz_yyyyz,   \
                             tr_y_xzzz_yyyyzz,  \
                             tr_y_xzzz_yyyzz,   \
                             tr_y_xzzz_yyyzzz,  \
                             tr_y_xzzz_yyzzz,   \
                             tr_y_xzzz_yyzzzz,  \
                             tr_y_xzzz_yzzzz,   \
                             tr_y_xzzz_yzzzzz,  \
                             tr_y_xzzz_zzzzz,   \
                             tr_y_xzzz_zzzzzz,  \
                             tr_y_zzz_xxxxxz,   \
                             tr_y_zzz_xxxxyz,   \
                             tr_y_zzz_xxxxzz,   \
                             tr_y_zzz_xxxyyz,   \
                             tr_y_zzz_xxxyzz,   \
                             tr_y_zzz_xxxzzz,   \
                             tr_y_zzz_xxyyyz,   \
                             tr_y_zzz_xxyyzz,   \
                             tr_y_zzz_xxyzzz,   \
                             tr_y_zzz_xxzzzz,   \
                             tr_y_zzz_xyyyyz,   \
                             tr_y_zzz_xyyyzz,   \
                             tr_y_zzz_xyyzzz,   \
                             tr_y_zzz_xyzzzz,   \
                             tr_y_zzz_xzzzzz,   \
                             tr_y_zzz_yyyyyy,   \
                             tr_y_zzz_yyyyyz,   \
                             tr_y_zzz_yyyyzz,   \
                             tr_y_zzz_yyyzzz,   \
                             tr_y_zzz_yyzzzz,   \
                             tr_y_zzz_yzzzzz,   \
                             tr_y_zzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzz_xxxxxx[i] = 2.0 * tr_y_xxz_xxxxxx[i] * fe_0 + tr_y_xxzz_xxxxxx[i] * pa_z[i];

        tr_y_xxzzz_xxxxxy[i] = 2.0 * tr_y_xxz_xxxxxy[i] * fe_0 + tr_y_xxzz_xxxxxy[i] * pa_z[i];

        tr_y_xxzzz_xxxxxz[i] = tr_y_zzz_xxxxxz[i] * fe_0 + 5.0 * tr_y_xzzz_xxxxz[i] * fe_0 + tr_y_xzzz_xxxxxz[i] * pa_x[i];

        tr_y_xxzzz_xxxxyy[i] = 2.0 * tr_y_xxz_xxxxyy[i] * fe_0 + tr_y_xxzz_xxxxyy[i] * pa_z[i];

        tr_y_xxzzz_xxxxyz[i] = tr_y_zzz_xxxxyz[i] * fe_0 + 4.0 * tr_y_xzzz_xxxyz[i] * fe_0 + tr_y_xzzz_xxxxyz[i] * pa_x[i];

        tr_y_xxzzz_xxxxzz[i] = tr_y_zzz_xxxxzz[i] * fe_0 + 4.0 * tr_y_xzzz_xxxzz[i] * fe_0 + tr_y_xzzz_xxxxzz[i] * pa_x[i];

        tr_y_xxzzz_xxxyyy[i] = 2.0 * tr_y_xxz_xxxyyy[i] * fe_0 + tr_y_xxzz_xxxyyy[i] * pa_z[i];

        tr_y_xxzzz_xxxyyz[i] = tr_y_zzz_xxxyyz[i] * fe_0 + 3.0 * tr_y_xzzz_xxyyz[i] * fe_0 + tr_y_xzzz_xxxyyz[i] * pa_x[i];

        tr_y_xxzzz_xxxyzz[i] = tr_y_zzz_xxxyzz[i] * fe_0 + 3.0 * tr_y_xzzz_xxyzz[i] * fe_0 + tr_y_xzzz_xxxyzz[i] * pa_x[i];

        tr_y_xxzzz_xxxzzz[i] = tr_y_zzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_xzzz_xxzzz[i] * fe_0 + tr_y_xzzz_xxxzzz[i] * pa_x[i];

        tr_y_xxzzz_xxyyyy[i] = 2.0 * tr_y_xxz_xxyyyy[i] * fe_0 + tr_y_xxzz_xxyyyy[i] * pa_z[i];

        tr_y_xxzzz_xxyyyz[i] = tr_y_zzz_xxyyyz[i] * fe_0 + 2.0 * tr_y_xzzz_xyyyz[i] * fe_0 + tr_y_xzzz_xxyyyz[i] * pa_x[i];

        tr_y_xxzzz_xxyyzz[i] = tr_y_zzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_xzzz_xyyzz[i] * fe_0 + tr_y_xzzz_xxyyzz[i] * pa_x[i];

        tr_y_xxzzz_xxyzzz[i] = tr_y_zzz_xxyzzz[i] * fe_0 + 2.0 * tr_y_xzzz_xyzzz[i] * fe_0 + tr_y_xzzz_xxyzzz[i] * pa_x[i];

        tr_y_xxzzz_xxzzzz[i] = tr_y_zzz_xxzzzz[i] * fe_0 + 2.0 * tr_y_xzzz_xzzzz[i] * fe_0 + tr_y_xzzz_xxzzzz[i] * pa_x[i];

        tr_y_xxzzz_xyyyyy[i] = 2.0 * tr_y_xxz_xyyyyy[i] * fe_0 + tr_y_xxzz_xyyyyy[i] * pa_z[i];

        tr_y_xxzzz_xyyyyz[i] = tr_y_zzz_xyyyyz[i] * fe_0 + tr_y_xzzz_yyyyz[i] * fe_0 + tr_y_xzzz_xyyyyz[i] * pa_x[i];

        tr_y_xxzzz_xyyyzz[i] = tr_y_zzz_xyyyzz[i] * fe_0 + tr_y_xzzz_yyyzz[i] * fe_0 + tr_y_xzzz_xyyyzz[i] * pa_x[i];

        tr_y_xxzzz_xyyzzz[i] = tr_y_zzz_xyyzzz[i] * fe_0 + tr_y_xzzz_yyzzz[i] * fe_0 + tr_y_xzzz_xyyzzz[i] * pa_x[i];

        tr_y_xxzzz_xyzzzz[i] = tr_y_zzz_xyzzzz[i] * fe_0 + tr_y_xzzz_yzzzz[i] * fe_0 + tr_y_xzzz_xyzzzz[i] * pa_x[i];

        tr_y_xxzzz_xzzzzz[i] = tr_y_zzz_xzzzzz[i] * fe_0 + tr_y_xzzz_zzzzz[i] * fe_0 + tr_y_xzzz_xzzzzz[i] * pa_x[i];

        tr_y_xxzzz_yyyyyy[i] = tr_y_zzz_yyyyyy[i] * fe_0 + tr_y_xzzz_yyyyyy[i] * pa_x[i];

        tr_y_xxzzz_yyyyyz[i] = tr_y_zzz_yyyyyz[i] * fe_0 + tr_y_xzzz_yyyyyz[i] * pa_x[i];

        tr_y_xxzzz_yyyyzz[i] = tr_y_zzz_yyyyzz[i] * fe_0 + tr_y_xzzz_yyyyzz[i] * pa_x[i];

        tr_y_xxzzz_yyyzzz[i] = tr_y_zzz_yyyzzz[i] * fe_0 + tr_y_xzzz_yyyzzz[i] * pa_x[i];

        tr_y_xxzzz_yyzzzz[i] = tr_y_zzz_yyzzzz[i] * fe_0 + tr_y_xzzz_yyzzzz[i] * pa_x[i];

        tr_y_xxzzz_yzzzzz[i] = tr_y_zzz_yzzzzz[i] * fe_0 + tr_y_xzzz_yzzzzz[i] * pa_x[i];

        tr_y_xxzzz_zzzzzz[i] = tr_y_zzz_zzzzzz[i] * fe_0 + tr_y_xzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 868-896 components of targeted buffer : HI

    auto tr_y_xyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 868);

    auto tr_y_xyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 869);

    auto tr_y_xyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 870);

    auto tr_y_xyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 871);

    auto tr_y_xyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 872);

    auto tr_y_xyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 873);

    auto tr_y_xyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 874);

    auto tr_y_xyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 875);

    auto tr_y_xyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 876);

    auto tr_y_xyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 877);

    auto tr_y_xyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 878);

    auto tr_y_xyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 879);

    auto tr_y_xyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 880);

    auto tr_y_xyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 881);

    auto tr_y_xyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 882);

    auto tr_y_xyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 883);

    auto tr_y_xyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 884);

    auto tr_y_xyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 885);

    auto tr_y_xyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 886);

    auto tr_y_xyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 887);

    auto tr_y_xyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 888);

    auto tr_y_xyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 889);

    auto tr_y_xyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 890);

    auto tr_y_xyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 891);

    auto tr_y_xyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 892);

    auto tr_y_xyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 893);

    auto tr_y_xyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 894);

    auto tr_y_xyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 895);

#pragma omp simd aligned(pa_x,                  \
                             tr_y_xyyyy_xxxxxx, \
                             tr_y_xyyyy_xxxxxy, \
                             tr_y_xyyyy_xxxxxz, \
                             tr_y_xyyyy_xxxxyy, \
                             tr_y_xyyyy_xxxxyz, \
                             tr_y_xyyyy_xxxxzz, \
                             tr_y_xyyyy_xxxyyy, \
                             tr_y_xyyyy_xxxyyz, \
                             tr_y_xyyyy_xxxyzz, \
                             tr_y_xyyyy_xxxzzz, \
                             tr_y_xyyyy_xxyyyy, \
                             tr_y_xyyyy_xxyyyz, \
                             tr_y_xyyyy_xxyyzz, \
                             tr_y_xyyyy_xxyzzz, \
                             tr_y_xyyyy_xxzzzz, \
                             tr_y_xyyyy_xyyyyy, \
                             tr_y_xyyyy_xyyyyz, \
                             tr_y_xyyyy_xyyyzz, \
                             tr_y_xyyyy_xyyzzz, \
                             tr_y_xyyyy_xyzzzz, \
                             tr_y_xyyyy_xzzzzz, \
                             tr_y_xyyyy_yyyyyy, \
                             tr_y_xyyyy_yyyyyz, \
                             tr_y_xyyyy_yyyyzz, \
                             tr_y_xyyyy_yyyzzz, \
                             tr_y_xyyyy_yyzzzz, \
                             tr_y_xyyyy_yzzzzz, \
                             tr_y_xyyyy_zzzzzz, \
                             tr_y_yyyy_xxxxx,   \
                             tr_y_yyyy_xxxxxx,  \
                             tr_y_yyyy_xxxxxy,  \
                             tr_y_yyyy_xxxxxz,  \
                             tr_y_yyyy_xxxxy,   \
                             tr_y_yyyy_xxxxyy,  \
                             tr_y_yyyy_xxxxyz,  \
                             tr_y_yyyy_xxxxz,   \
                             tr_y_yyyy_xxxxzz,  \
                             tr_y_yyyy_xxxyy,   \
                             tr_y_yyyy_xxxyyy,  \
                             tr_y_yyyy_xxxyyz,  \
                             tr_y_yyyy_xxxyz,   \
                             tr_y_yyyy_xxxyzz,  \
                             tr_y_yyyy_xxxzz,   \
                             tr_y_yyyy_xxxzzz,  \
                             tr_y_yyyy_xxyyy,   \
                             tr_y_yyyy_xxyyyy,  \
                             tr_y_yyyy_xxyyyz,  \
                             tr_y_yyyy_xxyyz,   \
                             tr_y_yyyy_xxyyzz,  \
                             tr_y_yyyy_xxyzz,   \
                             tr_y_yyyy_xxyzzz,  \
                             tr_y_yyyy_xxzzz,   \
                             tr_y_yyyy_xxzzzz,  \
                             tr_y_yyyy_xyyyy,   \
                             tr_y_yyyy_xyyyyy,  \
                             tr_y_yyyy_xyyyyz,  \
                             tr_y_yyyy_xyyyz,   \
                             tr_y_yyyy_xyyyzz,  \
                             tr_y_yyyy_xyyzz,   \
                             tr_y_yyyy_xyyzzz,  \
                             tr_y_yyyy_xyzzz,   \
                             tr_y_yyyy_xyzzzz,  \
                             tr_y_yyyy_xzzzz,   \
                             tr_y_yyyy_xzzzzz,  \
                             tr_y_yyyy_yyyyy,   \
                             tr_y_yyyy_yyyyyy,  \
                             tr_y_yyyy_yyyyyz,  \
                             tr_y_yyyy_yyyyz,   \
                             tr_y_yyyy_yyyyzz,  \
                             tr_y_yyyy_yyyzz,   \
                             tr_y_yyyy_yyyzzz,  \
                             tr_y_yyyy_yyzzz,   \
                             tr_y_yyyy_yyzzzz,  \
                             tr_y_yyyy_yzzzz,   \
                             tr_y_yyyy_yzzzzz,  \
                             tr_y_yyyy_zzzzz,   \
                             tr_y_yyyy_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyy_xxxxxx[i] = 6.0 * tr_y_yyyy_xxxxx[i] * fe_0 + tr_y_yyyy_xxxxxx[i] * pa_x[i];

        tr_y_xyyyy_xxxxxy[i] = 5.0 * tr_y_yyyy_xxxxy[i] * fe_0 + tr_y_yyyy_xxxxxy[i] * pa_x[i];

        tr_y_xyyyy_xxxxxz[i] = 5.0 * tr_y_yyyy_xxxxz[i] * fe_0 + tr_y_yyyy_xxxxxz[i] * pa_x[i];

        tr_y_xyyyy_xxxxyy[i] = 4.0 * tr_y_yyyy_xxxyy[i] * fe_0 + tr_y_yyyy_xxxxyy[i] * pa_x[i];

        tr_y_xyyyy_xxxxyz[i] = 4.0 * tr_y_yyyy_xxxyz[i] * fe_0 + tr_y_yyyy_xxxxyz[i] * pa_x[i];

        tr_y_xyyyy_xxxxzz[i] = 4.0 * tr_y_yyyy_xxxzz[i] * fe_0 + tr_y_yyyy_xxxxzz[i] * pa_x[i];

        tr_y_xyyyy_xxxyyy[i] = 3.0 * tr_y_yyyy_xxyyy[i] * fe_0 + tr_y_yyyy_xxxyyy[i] * pa_x[i];

        tr_y_xyyyy_xxxyyz[i] = 3.0 * tr_y_yyyy_xxyyz[i] * fe_0 + tr_y_yyyy_xxxyyz[i] * pa_x[i];

        tr_y_xyyyy_xxxyzz[i] = 3.0 * tr_y_yyyy_xxyzz[i] * fe_0 + tr_y_yyyy_xxxyzz[i] * pa_x[i];

        tr_y_xyyyy_xxxzzz[i] = 3.0 * tr_y_yyyy_xxzzz[i] * fe_0 + tr_y_yyyy_xxxzzz[i] * pa_x[i];

        tr_y_xyyyy_xxyyyy[i] = 2.0 * tr_y_yyyy_xyyyy[i] * fe_0 + tr_y_yyyy_xxyyyy[i] * pa_x[i];

        tr_y_xyyyy_xxyyyz[i] = 2.0 * tr_y_yyyy_xyyyz[i] * fe_0 + tr_y_yyyy_xxyyyz[i] * pa_x[i];

        tr_y_xyyyy_xxyyzz[i] = 2.0 * tr_y_yyyy_xyyzz[i] * fe_0 + tr_y_yyyy_xxyyzz[i] * pa_x[i];

        tr_y_xyyyy_xxyzzz[i] = 2.0 * tr_y_yyyy_xyzzz[i] * fe_0 + tr_y_yyyy_xxyzzz[i] * pa_x[i];

        tr_y_xyyyy_xxzzzz[i] = 2.0 * tr_y_yyyy_xzzzz[i] * fe_0 + tr_y_yyyy_xxzzzz[i] * pa_x[i];

        tr_y_xyyyy_xyyyyy[i] = tr_y_yyyy_yyyyy[i] * fe_0 + tr_y_yyyy_xyyyyy[i] * pa_x[i];

        tr_y_xyyyy_xyyyyz[i] = tr_y_yyyy_yyyyz[i] * fe_0 + tr_y_yyyy_xyyyyz[i] * pa_x[i];

        tr_y_xyyyy_xyyyzz[i] = tr_y_yyyy_yyyzz[i] * fe_0 + tr_y_yyyy_xyyyzz[i] * pa_x[i];

        tr_y_xyyyy_xyyzzz[i] = tr_y_yyyy_yyzzz[i] * fe_0 + tr_y_yyyy_xyyzzz[i] * pa_x[i];

        tr_y_xyyyy_xyzzzz[i] = tr_y_yyyy_yzzzz[i] * fe_0 + tr_y_yyyy_xyzzzz[i] * pa_x[i];

        tr_y_xyyyy_xzzzzz[i] = tr_y_yyyy_zzzzz[i] * fe_0 + tr_y_yyyy_xzzzzz[i] * pa_x[i];

        tr_y_xyyyy_yyyyyy[i] = tr_y_yyyy_yyyyyy[i] * pa_x[i];

        tr_y_xyyyy_yyyyyz[i] = tr_y_yyyy_yyyyyz[i] * pa_x[i];

        tr_y_xyyyy_yyyyzz[i] = tr_y_yyyy_yyyyzz[i] * pa_x[i];

        tr_y_xyyyy_yyyzzz[i] = tr_y_yyyy_yyyzzz[i] * pa_x[i];

        tr_y_xyyyy_yyzzzz[i] = tr_y_yyyy_yyzzzz[i] * pa_x[i];

        tr_y_xyyyy_yzzzzz[i] = tr_y_yyyy_yzzzzz[i] * pa_x[i];

        tr_y_xyyyy_zzzzzz[i] = tr_y_yyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 896-924 components of targeted buffer : HI

    auto tr_y_xyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 896);

    auto tr_y_xyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 897);

    auto tr_y_xyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 898);

    auto tr_y_xyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 899);

    auto tr_y_xyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 900);

    auto tr_y_xyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 901);

    auto tr_y_xyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 902);

    auto tr_y_xyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 903);

    auto tr_y_xyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 904);

    auto tr_y_xyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 905);

    auto tr_y_xyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 906);

    auto tr_y_xyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 907);

    auto tr_y_xyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 908);

    auto tr_y_xyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 909);

    auto tr_y_xyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 910);

    auto tr_y_xyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 911);

    auto tr_y_xyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 912);

    auto tr_y_xyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 913);

    auto tr_y_xyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 914);

    auto tr_y_xyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 915);

    auto tr_y_xyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 916);

    auto tr_y_xyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 917);

    auto tr_y_xyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 918);

    auto tr_y_xyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 919);

    auto tr_y_xyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 920);

    auto tr_y_xyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 921);

    auto tr_y_xyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 922);

    auto tr_y_xyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 923);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_y_xyyy_xxxxxx,  \
                             tr_y_xyyy_xxxxxy,  \
                             tr_y_xyyy_xxxxyy,  \
                             tr_y_xyyy_xxxyyy,  \
                             tr_y_xyyy_xxyyyy,  \
                             tr_y_xyyy_xyyyyy,  \
                             tr_y_xyyyz_xxxxxx, \
                             tr_y_xyyyz_xxxxxy, \
                             tr_y_xyyyz_xxxxxz, \
                             tr_y_xyyyz_xxxxyy, \
                             tr_y_xyyyz_xxxxyz, \
                             tr_y_xyyyz_xxxxzz, \
                             tr_y_xyyyz_xxxyyy, \
                             tr_y_xyyyz_xxxyyz, \
                             tr_y_xyyyz_xxxyzz, \
                             tr_y_xyyyz_xxxzzz, \
                             tr_y_xyyyz_xxyyyy, \
                             tr_y_xyyyz_xxyyyz, \
                             tr_y_xyyyz_xxyyzz, \
                             tr_y_xyyyz_xxyzzz, \
                             tr_y_xyyyz_xxzzzz, \
                             tr_y_xyyyz_xyyyyy, \
                             tr_y_xyyyz_xyyyyz, \
                             tr_y_xyyyz_xyyyzz, \
                             tr_y_xyyyz_xyyzzz, \
                             tr_y_xyyyz_xyzzzz, \
                             tr_y_xyyyz_xzzzzz, \
                             tr_y_xyyyz_yyyyyy, \
                             tr_y_xyyyz_yyyyyz, \
                             tr_y_xyyyz_yyyyzz, \
                             tr_y_xyyyz_yyyzzz, \
                             tr_y_xyyyz_yyzzzz, \
                             tr_y_xyyyz_yzzzzz, \
                             tr_y_xyyyz_zzzzzz, \
                             tr_y_yyyz_xxxxxz,  \
                             tr_y_yyyz_xxxxyz,  \
                             tr_y_yyyz_xxxxz,   \
                             tr_y_yyyz_xxxxzz,  \
                             tr_y_yyyz_xxxyyz,  \
                             tr_y_yyyz_xxxyz,   \
                             tr_y_yyyz_xxxyzz,  \
                             tr_y_yyyz_xxxzz,   \
                             tr_y_yyyz_xxxzzz,  \
                             tr_y_yyyz_xxyyyz,  \
                             tr_y_yyyz_xxyyz,   \
                             tr_y_yyyz_xxyyzz,  \
                             tr_y_yyyz_xxyzz,   \
                             tr_y_yyyz_xxyzzz,  \
                             tr_y_yyyz_xxzzz,   \
                             tr_y_yyyz_xxzzzz,  \
                             tr_y_yyyz_xyyyyz,  \
                             tr_y_yyyz_xyyyz,   \
                             tr_y_yyyz_xyyyzz,  \
                             tr_y_yyyz_xyyzz,   \
                             tr_y_yyyz_xyyzzz,  \
                             tr_y_yyyz_xyzzz,   \
                             tr_y_yyyz_xyzzzz,  \
                             tr_y_yyyz_xzzzz,   \
                             tr_y_yyyz_xzzzzz,  \
                             tr_y_yyyz_yyyyyy,  \
                             tr_y_yyyz_yyyyyz,  \
                             tr_y_yyyz_yyyyz,   \
                             tr_y_yyyz_yyyyzz,  \
                             tr_y_yyyz_yyyzz,   \
                             tr_y_yyyz_yyyzzz,  \
                             tr_y_yyyz_yyzzz,   \
                             tr_y_yyyz_yyzzzz,  \
                             tr_y_yyyz_yzzzz,   \
                             tr_y_yyyz_yzzzzz,  \
                             tr_y_yyyz_zzzzz,   \
                             tr_y_yyyz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyz_xxxxxx[i] = tr_y_xyyy_xxxxxx[i] * pa_z[i];

        tr_y_xyyyz_xxxxxy[i] = tr_y_xyyy_xxxxxy[i] * pa_z[i];

        tr_y_xyyyz_xxxxxz[i] = 5.0 * tr_y_yyyz_xxxxz[i] * fe_0 + tr_y_yyyz_xxxxxz[i] * pa_x[i];

        tr_y_xyyyz_xxxxyy[i] = tr_y_xyyy_xxxxyy[i] * pa_z[i];

        tr_y_xyyyz_xxxxyz[i] = 4.0 * tr_y_yyyz_xxxyz[i] * fe_0 + tr_y_yyyz_xxxxyz[i] * pa_x[i];

        tr_y_xyyyz_xxxxzz[i] = 4.0 * tr_y_yyyz_xxxzz[i] * fe_0 + tr_y_yyyz_xxxxzz[i] * pa_x[i];

        tr_y_xyyyz_xxxyyy[i] = tr_y_xyyy_xxxyyy[i] * pa_z[i];

        tr_y_xyyyz_xxxyyz[i] = 3.0 * tr_y_yyyz_xxyyz[i] * fe_0 + tr_y_yyyz_xxxyyz[i] * pa_x[i];

        tr_y_xyyyz_xxxyzz[i] = 3.0 * tr_y_yyyz_xxyzz[i] * fe_0 + tr_y_yyyz_xxxyzz[i] * pa_x[i];

        tr_y_xyyyz_xxxzzz[i] = 3.0 * tr_y_yyyz_xxzzz[i] * fe_0 + tr_y_yyyz_xxxzzz[i] * pa_x[i];

        tr_y_xyyyz_xxyyyy[i] = tr_y_xyyy_xxyyyy[i] * pa_z[i];

        tr_y_xyyyz_xxyyyz[i] = 2.0 * tr_y_yyyz_xyyyz[i] * fe_0 + tr_y_yyyz_xxyyyz[i] * pa_x[i];

        tr_y_xyyyz_xxyyzz[i] = 2.0 * tr_y_yyyz_xyyzz[i] * fe_0 + tr_y_yyyz_xxyyzz[i] * pa_x[i];

        tr_y_xyyyz_xxyzzz[i] = 2.0 * tr_y_yyyz_xyzzz[i] * fe_0 + tr_y_yyyz_xxyzzz[i] * pa_x[i];

        tr_y_xyyyz_xxzzzz[i] = 2.0 * tr_y_yyyz_xzzzz[i] * fe_0 + tr_y_yyyz_xxzzzz[i] * pa_x[i];

        tr_y_xyyyz_xyyyyy[i] = tr_y_xyyy_xyyyyy[i] * pa_z[i];

        tr_y_xyyyz_xyyyyz[i] = tr_y_yyyz_yyyyz[i] * fe_0 + tr_y_yyyz_xyyyyz[i] * pa_x[i];

        tr_y_xyyyz_xyyyzz[i] = tr_y_yyyz_yyyzz[i] * fe_0 + tr_y_yyyz_xyyyzz[i] * pa_x[i];

        tr_y_xyyyz_xyyzzz[i] = tr_y_yyyz_yyzzz[i] * fe_0 + tr_y_yyyz_xyyzzz[i] * pa_x[i];

        tr_y_xyyyz_xyzzzz[i] = tr_y_yyyz_yzzzz[i] * fe_0 + tr_y_yyyz_xyzzzz[i] * pa_x[i];

        tr_y_xyyyz_xzzzzz[i] = tr_y_yyyz_zzzzz[i] * fe_0 + tr_y_yyyz_xzzzzz[i] * pa_x[i];

        tr_y_xyyyz_yyyyyy[i] = tr_y_yyyz_yyyyyy[i] * pa_x[i];

        tr_y_xyyyz_yyyyyz[i] = tr_y_yyyz_yyyyyz[i] * pa_x[i];

        tr_y_xyyyz_yyyyzz[i] = tr_y_yyyz_yyyyzz[i] * pa_x[i];

        tr_y_xyyyz_yyyzzz[i] = tr_y_yyyz_yyyzzz[i] * pa_x[i];

        tr_y_xyyyz_yyzzzz[i] = tr_y_yyyz_yyzzzz[i] * pa_x[i];

        tr_y_xyyyz_yzzzzz[i] = tr_y_yyyz_yzzzzz[i] * pa_x[i];

        tr_y_xyyyz_zzzzzz[i] = tr_y_yyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 924-952 components of targeted buffer : HI

    auto tr_y_xyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 924);

    auto tr_y_xyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 925);

    auto tr_y_xyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 926);

    auto tr_y_xyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 927);

    auto tr_y_xyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 928);

    auto tr_y_xyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 929);

    auto tr_y_xyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 930);

    auto tr_y_xyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 931);

    auto tr_y_xyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 932);

    auto tr_y_xyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 933);

    auto tr_y_xyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 934);

    auto tr_y_xyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 935);

    auto tr_y_xyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 936);

    auto tr_y_xyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 937);

    auto tr_y_xyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 938);

    auto tr_y_xyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 939);

    auto tr_y_xyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 940);

    auto tr_y_xyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 941);

    auto tr_y_xyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 942);

    auto tr_y_xyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 943);

    auto tr_y_xyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 944);

    auto tr_y_xyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 945);

    auto tr_y_xyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 946);

    auto tr_y_xyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 947);

    auto tr_y_xyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 948);

    auto tr_y_xyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 949);

    auto tr_y_xyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 950);

    auto tr_y_xyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 951);

#pragma omp simd aligned(pa_x,                  \
                             tr_y_xyyzz_xxxxxx, \
                             tr_y_xyyzz_xxxxxy, \
                             tr_y_xyyzz_xxxxxz, \
                             tr_y_xyyzz_xxxxyy, \
                             tr_y_xyyzz_xxxxyz, \
                             tr_y_xyyzz_xxxxzz, \
                             tr_y_xyyzz_xxxyyy, \
                             tr_y_xyyzz_xxxyyz, \
                             tr_y_xyyzz_xxxyzz, \
                             tr_y_xyyzz_xxxzzz, \
                             tr_y_xyyzz_xxyyyy, \
                             tr_y_xyyzz_xxyyyz, \
                             tr_y_xyyzz_xxyyzz, \
                             tr_y_xyyzz_xxyzzz, \
                             tr_y_xyyzz_xxzzzz, \
                             tr_y_xyyzz_xyyyyy, \
                             tr_y_xyyzz_xyyyyz, \
                             tr_y_xyyzz_xyyyzz, \
                             tr_y_xyyzz_xyyzzz, \
                             tr_y_xyyzz_xyzzzz, \
                             tr_y_xyyzz_xzzzzz, \
                             tr_y_xyyzz_yyyyyy, \
                             tr_y_xyyzz_yyyyyz, \
                             tr_y_xyyzz_yyyyzz, \
                             tr_y_xyyzz_yyyzzz, \
                             tr_y_xyyzz_yyzzzz, \
                             tr_y_xyyzz_yzzzzz, \
                             tr_y_xyyzz_zzzzzz, \
                             tr_y_yyzz_xxxxx,   \
                             tr_y_yyzz_xxxxxx,  \
                             tr_y_yyzz_xxxxxy,  \
                             tr_y_yyzz_xxxxxz,  \
                             tr_y_yyzz_xxxxy,   \
                             tr_y_yyzz_xxxxyy,  \
                             tr_y_yyzz_xxxxyz,  \
                             tr_y_yyzz_xxxxz,   \
                             tr_y_yyzz_xxxxzz,  \
                             tr_y_yyzz_xxxyy,   \
                             tr_y_yyzz_xxxyyy,  \
                             tr_y_yyzz_xxxyyz,  \
                             tr_y_yyzz_xxxyz,   \
                             tr_y_yyzz_xxxyzz,  \
                             tr_y_yyzz_xxxzz,   \
                             tr_y_yyzz_xxxzzz,  \
                             tr_y_yyzz_xxyyy,   \
                             tr_y_yyzz_xxyyyy,  \
                             tr_y_yyzz_xxyyyz,  \
                             tr_y_yyzz_xxyyz,   \
                             tr_y_yyzz_xxyyzz,  \
                             tr_y_yyzz_xxyzz,   \
                             tr_y_yyzz_xxyzzz,  \
                             tr_y_yyzz_xxzzz,   \
                             tr_y_yyzz_xxzzzz,  \
                             tr_y_yyzz_xyyyy,   \
                             tr_y_yyzz_xyyyyy,  \
                             tr_y_yyzz_xyyyyz,  \
                             tr_y_yyzz_xyyyz,   \
                             tr_y_yyzz_xyyyzz,  \
                             tr_y_yyzz_xyyzz,   \
                             tr_y_yyzz_xyyzzz,  \
                             tr_y_yyzz_xyzzz,   \
                             tr_y_yyzz_xyzzzz,  \
                             tr_y_yyzz_xzzzz,   \
                             tr_y_yyzz_xzzzzz,  \
                             tr_y_yyzz_yyyyy,   \
                             tr_y_yyzz_yyyyyy,  \
                             tr_y_yyzz_yyyyyz,  \
                             tr_y_yyzz_yyyyz,   \
                             tr_y_yyzz_yyyyzz,  \
                             tr_y_yyzz_yyyzz,   \
                             tr_y_yyzz_yyyzzz,  \
                             tr_y_yyzz_yyzzz,   \
                             tr_y_yyzz_yyzzzz,  \
                             tr_y_yyzz_yzzzz,   \
                             tr_y_yyzz_yzzzzz,  \
                             tr_y_yyzz_zzzzz,   \
                             tr_y_yyzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzz_xxxxxx[i] = 6.0 * tr_y_yyzz_xxxxx[i] * fe_0 + tr_y_yyzz_xxxxxx[i] * pa_x[i];

        tr_y_xyyzz_xxxxxy[i] = 5.0 * tr_y_yyzz_xxxxy[i] * fe_0 + tr_y_yyzz_xxxxxy[i] * pa_x[i];

        tr_y_xyyzz_xxxxxz[i] = 5.0 * tr_y_yyzz_xxxxz[i] * fe_0 + tr_y_yyzz_xxxxxz[i] * pa_x[i];

        tr_y_xyyzz_xxxxyy[i] = 4.0 * tr_y_yyzz_xxxyy[i] * fe_0 + tr_y_yyzz_xxxxyy[i] * pa_x[i];

        tr_y_xyyzz_xxxxyz[i] = 4.0 * tr_y_yyzz_xxxyz[i] * fe_0 + tr_y_yyzz_xxxxyz[i] * pa_x[i];

        tr_y_xyyzz_xxxxzz[i] = 4.0 * tr_y_yyzz_xxxzz[i] * fe_0 + tr_y_yyzz_xxxxzz[i] * pa_x[i];

        tr_y_xyyzz_xxxyyy[i] = 3.0 * tr_y_yyzz_xxyyy[i] * fe_0 + tr_y_yyzz_xxxyyy[i] * pa_x[i];

        tr_y_xyyzz_xxxyyz[i] = 3.0 * tr_y_yyzz_xxyyz[i] * fe_0 + tr_y_yyzz_xxxyyz[i] * pa_x[i];

        tr_y_xyyzz_xxxyzz[i] = 3.0 * tr_y_yyzz_xxyzz[i] * fe_0 + tr_y_yyzz_xxxyzz[i] * pa_x[i];

        tr_y_xyyzz_xxxzzz[i] = 3.0 * tr_y_yyzz_xxzzz[i] * fe_0 + tr_y_yyzz_xxxzzz[i] * pa_x[i];

        tr_y_xyyzz_xxyyyy[i] = 2.0 * tr_y_yyzz_xyyyy[i] * fe_0 + tr_y_yyzz_xxyyyy[i] * pa_x[i];

        tr_y_xyyzz_xxyyyz[i] = 2.0 * tr_y_yyzz_xyyyz[i] * fe_0 + tr_y_yyzz_xxyyyz[i] * pa_x[i];

        tr_y_xyyzz_xxyyzz[i] = 2.0 * tr_y_yyzz_xyyzz[i] * fe_0 + tr_y_yyzz_xxyyzz[i] * pa_x[i];

        tr_y_xyyzz_xxyzzz[i] = 2.0 * tr_y_yyzz_xyzzz[i] * fe_0 + tr_y_yyzz_xxyzzz[i] * pa_x[i];

        tr_y_xyyzz_xxzzzz[i] = 2.0 * tr_y_yyzz_xzzzz[i] * fe_0 + tr_y_yyzz_xxzzzz[i] * pa_x[i];

        tr_y_xyyzz_xyyyyy[i] = tr_y_yyzz_yyyyy[i] * fe_0 + tr_y_yyzz_xyyyyy[i] * pa_x[i];

        tr_y_xyyzz_xyyyyz[i] = tr_y_yyzz_yyyyz[i] * fe_0 + tr_y_yyzz_xyyyyz[i] * pa_x[i];

        tr_y_xyyzz_xyyyzz[i] = tr_y_yyzz_yyyzz[i] * fe_0 + tr_y_yyzz_xyyyzz[i] * pa_x[i];

        tr_y_xyyzz_xyyzzz[i] = tr_y_yyzz_yyzzz[i] * fe_0 + tr_y_yyzz_xyyzzz[i] * pa_x[i];

        tr_y_xyyzz_xyzzzz[i] = tr_y_yyzz_yzzzz[i] * fe_0 + tr_y_yyzz_xyzzzz[i] * pa_x[i];

        tr_y_xyyzz_xzzzzz[i] = tr_y_yyzz_zzzzz[i] * fe_0 + tr_y_yyzz_xzzzzz[i] * pa_x[i];

        tr_y_xyyzz_yyyyyy[i] = tr_y_yyzz_yyyyyy[i] * pa_x[i];

        tr_y_xyyzz_yyyyyz[i] = tr_y_yyzz_yyyyyz[i] * pa_x[i];

        tr_y_xyyzz_yyyyzz[i] = tr_y_yyzz_yyyyzz[i] * pa_x[i];

        tr_y_xyyzz_yyyzzz[i] = tr_y_yyzz_yyyzzz[i] * pa_x[i];

        tr_y_xyyzz_yyzzzz[i] = tr_y_yyzz_yyzzzz[i] * pa_x[i];

        tr_y_xyyzz_yzzzzz[i] = tr_y_yyzz_yzzzzz[i] * pa_x[i];

        tr_y_xyyzz_zzzzzz[i] = tr_y_yyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 952-980 components of targeted buffer : HI

    auto tr_y_xyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 952);

    auto tr_y_xyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 953);

    auto tr_y_xyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 954);

    auto tr_y_xyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 955);

    auto tr_y_xyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 956);

    auto tr_y_xyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 957);

    auto tr_y_xyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 958);

    auto tr_y_xyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 959);

    auto tr_y_xyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 960);

    auto tr_y_xyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 961);

    auto tr_y_xyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 962);

    auto tr_y_xyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 963);

    auto tr_y_xyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 964);

    auto tr_y_xyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 965);

    auto tr_y_xyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 966);

    auto tr_y_xyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 967);

    auto tr_y_xyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 968);

    auto tr_y_xyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 969);

    auto tr_y_xyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 970);

    auto tr_y_xyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 971);

    auto tr_y_xyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 972);

    auto tr_y_xyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 973);

    auto tr_y_xyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 974);

    auto tr_y_xyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 975);

    auto tr_y_xyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 976);

    auto tr_y_xyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 977);

    auto tr_y_xyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 978);

    auto tr_y_xyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 979);

#pragma omp simd aligned(pa_x,                  \
                             tr_y_xyzzz_xxxxxx, \
                             tr_y_xyzzz_xxxxxy, \
                             tr_y_xyzzz_xxxxxz, \
                             tr_y_xyzzz_xxxxyy, \
                             tr_y_xyzzz_xxxxyz, \
                             tr_y_xyzzz_xxxxzz, \
                             tr_y_xyzzz_xxxyyy, \
                             tr_y_xyzzz_xxxyyz, \
                             tr_y_xyzzz_xxxyzz, \
                             tr_y_xyzzz_xxxzzz, \
                             tr_y_xyzzz_xxyyyy, \
                             tr_y_xyzzz_xxyyyz, \
                             tr_y_xyzzz_xxyyzz, \
                             tr_y_xyzzz_xxyzzz, \
                             tr_y_xyzzz_xxzzzz, \
                             tr_y_xyzzz_xyyyyy, \
                             tr_y_xyzzz_xyyyyz, \
                             tr_y_xyzzz_xyyyzz, \
                             tr_y_xyzzz_xyyzzz, \
                             tr_y_xyzzz_xyzzzz, \
                             tr_y_xyzzz_xzzzzz, \
                             tr_y_xyzzz_yyyyyy, \
                             tr_y_xyzzz_yyyyyz, \
                             tr_y_xyzzz_yyyyzz, \
                             tr_y_xyzzz_yyyzzz, \
                             tr_y_xyzzz_yyzzzz, \
                             tr_y_xyzzz_yzzzzz, \
                             tr_y_xyzzz_zzzzzz, \
                             tr_y_yzzz_xxxxx,   \
                             tr_y_yzzz_xxxxxx,  \
                             tr_y_yzzz_xxxxxy,  \
                             tr_y_yzzz_xxxxxz,  \
                             tr_y_yzzz_xxxxy,   \
                             tr_y_yzzz_xxxxyy,  \
                             tr_y_yzzz_xxxxyz,  \
                             tr_y_yzzz_xxxxz,   \
                             tr_y_yzzz_xxxxzz,  \
                             tr_y_yzzz_xxxyy,   \
                             tr_y_yzzz_xxxyyy,  \
                             tr_y_yzzz_xxxyyz,  \
                             tr_y_yzzz_xxxyz,   \
                             tr_y_yzzz_xxxyzz,  \
                             tr_y_yzzz_xxxzz,   \
                             tr_y_yzzz_xxxzzz,  \
                             tr_y_yzzz_xxyyy,   \
                             tr_y_yzzz_xxyyyy,  \
                             tr_y_yzzz_xxyyyz,  \
                             tr_y_yzzz_xxyyz,   \
                             tr_y_yzzz_xxyyzz,  \
                             tr_y_yzzz_xxyzz,   \
                             tr_y_yzzz_xxyzzz,  \
                             tr_y_yzzz_xxzzz,   \
                             tr_y_yzzz_xxzzzz,  \
                             tr_y_yzzz_xyyyy,   \
                             tr_y_yzzz_xyyyyy,  \
                             tr_y_yzzz_xyyyyz,  \
                             tr_y_yzzz_xyyyz,   \
                             tr_y_yzzz_xyyyzz,  \
                             tr_y_yzzz_xyyzz,   \
                             tr_y_yzzz_xyyzzz,  \
                             tr_y_yzzz_xyzzz,   \
                             tr_y_yzzz_xyzzzz,  \
                             tr_y_yzzz_xzzzz,   \
                             tr_y_yzzz_xzzzzz,  \
                             tr_y_yzzz_yyyyy,   \
                             tr_y_yzzz_yyyyyy,  \
                             tr_y_yzzz_yyyyyz,  \
                             tr_y_yzzz_yyyyz,   \
                             tr_y_yzzz_yyyyzz,  \
                             tr_y_yzzz_yyyzz,   \
                             tr_y_yzzz_yyyzzz,  \
                             tr_y_yzzz_yyzzz,   \
                             tr_y_yzzz_yyzzzz,  \
                             tr_y_yzzz_yzzzz,   \
                             tr_y_yzzz_yzzzzz,  \
                             tr_y_yzzz_zzzzz,   \
                             tr_y_yzzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzz_xxxxxx[i] = 6.0 * tr_y_yzzz_xxxxx[i] * fe_0 + tr_y_yzzz_xxxxxx[i] * pa_x[i];

        tr_y_xyzzz_xxxxxy[i] = 5.0 * tr_y_yzzz_xxxxy[i] * fe_0 + tr_y_yzzz_xxxxxy[i] * pa_x[i];

        tr_y_xyzzz_xxxxxz[i] = 5.0 * tr_y_yzzz_xxxxz[i] * fe_0 + tr_y_yzzz_xxxxxz[i] * pa_x[i];

        tr_y_xyzzz_xxxxyy[i] = 4.0 * tr_y_yzzz_xxxyy[i] * fe_0 + tr_y_yzzz_xxxxyy[i] * pa_x[i];

        tr_y_xyzzz_xxxxyz[i] = 4.0 * tr_y_yzzz_xxxyz[i] * fe_0 + tr_y_yzzz_xxxxyz[i] * pa_x[i];

        tr_y_xyzzz_xxxxzz[i] = 4.0 * tr_y_yzzz_xxxzz[i] * fe_0 + tr_y_yzzz_xxxxzz[i] * pa_x[i];

        tr_y_xyzzz_xxxyyy[i] = 3.0 * tr_y_yzzz_xxyyy[i] * fe_0 + tr_y_yzzz_xxxyyy[i] * pa_x[i];

        tr_y_xyzzz_xxxyyz[i] = 3.0 * tr_y_yzzz_xxyyz[i] * fe_0 + tr_y_yzzz_xxxyyz[i] * pa_x[i];

        tr_y_xyzzz_xxxyzz[i] = 3.0 * tr_y_yzzz_xxyzz[i] * fe_0 + tr_y_yzzz_xxxyzz[i] * pa_x[i];

        tr_y_xyzzz_xxxzzz[i] = 3.0 * tr_y_yzzz_xxzzz[i] * fe_0 + tr_y_yzzz_xxxzzz[i] * pa_x[i];

        tr_y_xyzzz_xxyyyy[i] = 2.0 * tr_y_yzzz_xyyyy[i] * fe_0 + tr_y_yzzz_xxyyyy[i] * pa_x[i];

        tr_y_xyzzz_xxyyyz[i] = 2.0 * tr_y_yzzz_xyyyz[i] * fe_0 + tr_y_yzzz_xxyyyz[i] * pa_x[i];

        tr_y_xyzzz_xxyyzz[i] = 2.0 * tr_y_yzzz_xyyzz[i] * fe_0 + tr_y_yzzz_xxyyzz[i] * pa_x[i];

        tr_y_xyzzz_xxyzzz[i] = 2.0 * tr_y_yzzz_xyzzz[i] * fe_0 + tr_y_yzzz_xxyzzz[i] * pa_x[i];

        tr_y_xyzzz_xxzzzz[i] = 2.0 * tr_y_yzzz_xzzzz[i] * fe_0 + tr_y_yzzz_xxzzzz[i] * pa_x[i];

        tr_y_xyzzz_xyyyyy[i] = tr_y_yzzz_yyyyy[i] * fe_0 + tr_y_yzzz_xyyyyy[i] * pa_x[i];

        tr_y_xyzzz_xyyyyz[i] = tr_y_yzzz_yyyyz[i] * fe_0 + tr_y_yzzz_xyyyyz[i] * pa_x[i];

        tr_y_xyzzz_xyyyzz[i] = tr_y_yzzz_yyyzz[i] * fe_0 + tr_y_yzzz_xyyyzz[i] * pa_x[i];

        tr_y_xyzzz_xyyzzz[i] = tr_y_yzzz_yyzzz[i] * fe_0 + tr_y_yzzz_xyyzzz[i] * pa_x[i];

        tr_y_xyzzz_xyzzzz[i] = tr_y_yzzz_yzzzz[i] * fe_0 + tr_y_yzzz_xyzzzz[i] * pa_x[i];

        tr_y_xyzzz_xzzzzz[i] = tr_y_yzzz_zzzzz[i] * fe_0 + tr_y_yzzz_xzzzzz[i] * pa_x[i];

        tr_y_xyzzz_yyyyyy[i] = tr_y_yzzz_yyyyyy[i] * pa_x[i];

        tr_y_xyzzz_yyyyyz[i] = tr_y_yzzz_yyyyyz[i] * pa_x[i];

        tr_y_xyzzz_yyyyzz[i] = tr_y_yzzz_yyyyzz[i] * pa_x[i];

        tr_y_xyzzz_yyyzzz[i] = tr_y_yzzz_yyyzzz[i] * pa_x[i];

        tr_y_xyzzz_yyzzzz[i] = tr_y_yzzz_yyzzzz[i] * pa_x[i];

        tr_y_xyzzz_yzzzzz[i] = tr_y_yzzz_yzzzzz[i] * pa_x[i];

        tr_y_xyzzz_zzzzzz[i] = tr_y_yzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 980-1008 components of targeted buffer : HI

    auto tr_y_xzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 980);

    auto tr_y_xzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 981);

    auto tr_y_xzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 982);

    auto tr_y_xzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 983);

    auto tr_y_xzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 984);

    auto tr_y_xzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 985);

    auto tr_y_xzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 986);

    auto tr_y_xzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 987);

    auto tr_y_xzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 988);

    auto tr_y_xzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 989);

    auto tr_y_xzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 990);

    auto tr_y_xzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 991);

    auto tr_y_xzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 992);

    auto tr_y_xzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 993);

    auto tr_y_xzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 994);

    auto tr_y_xzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 995);

    auto tr_y_xzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 996);

    auto tr_y_xzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 997);

    auto tr_y_xzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 998);

    auto tr_y_xzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 999);

    auto tr_y_xzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1000);

    auto tr_y_xzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1001);

    auto tr_y_xzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1002);

    auto tr_y_xzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1003);

    auto tr_y_xzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1004);

    auto tr_y_xzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1005);

    auto tr_y_xzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1006);

    auto tr_y_xzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1007);

#pragma omp simd aligned(pa_x,                  \
                             tr_y_xzzzz_xxxxxx, \
                             tr_y_xzzzz_xxxxxy, \
                             tr_y_xzzzz_xxxxxz, \
                             tr_y_xzzzz_xxxxyy, \
                             tr_y_xzzzz_xxxxyz, \
                             tr_y_xzzzz_xxxxzz, \
                             tr_y_xzzzz_xxxyyy, \
                             tr_y_xzzzz_xxxyyz, \
                             tr_y_xzzzz_xxxyzz, \
                             tr_y_xzzzz_xxxzzz, \
                             tr_y_xzzzz_xxyyyy, \
                             tr_y_xzzzz_xxyyyz, \
                             tr_y_xzzzz_xxyyzz, \
                             tr_y_xzzzz_xxyzzz, \
                             tr_y_xzzzz_xxzzzz, \
                             tr_y_xzzzz_xyyyyy, \
                             tr_y_xzzzz_xyyyyz, \
                             tr_y_xzzzz_xyyyzz, \
                             tr_y_xzzzz_xyyzzz, \
                             tr_y_xzzzz_xyzzzz, \
                             tr_y_xzzzz_xzzzzz, \
                             tr_y_xzzzz_yyyyyy, \
                             tr_y_xzzzz_yyyyyz, \
                             tr_y_xzzzz_yyyyzz, \
                             tr_y_xzzzz_yyyzzz, \
                             tr_y_xzzzz_yyzzzz, \
                             tr_y_xzzzz_yzzzzz, \
                             tr_y_xzzzz_zzzzzz, \
                             tr_y_zzzz_xxxxx,   \
                             tr_y_zzzz_xxxxxx,  \
                             tr_y_zzzz_xxxxxy,  \
                             tr_y_zzzz_xxxxxz,  \
                             tr_y_zzzz_xxxxy,   \
                             tr_y_zzzz_xxxxyy,  \
                             tr_y_zzzz_xxxxyz,  \
                             tr_y_zzzz_xxxxz,   \
                             tr_y_zzzz_xxxxzz,  \
                             tr_y_zzzz_xxxyy,   \
                             tr_y_zzzz_xxxyyy,  \
                             tr_y_zzzz_xxxyyz,  \
                             tr_y_zzzz_xxxyz,   \
                             tr_y_zzzz_xxxyzz,  \
                             tr_y_zzzz_xxxzz,   \
                             tr_y_zzzz_xxxzzz,  \
                             tr_y_zzzz_xxyyy,   \
                             tr_y_zzzz_xxyyyy,  \
                             tr_y_zzzz_xxyyyz,  \
                             tr_y_zzzz_xxyyz,   \
                             tr_y_zzzz_xxyyzz,  \
                             tr_y_zzzz_xxyzz,   \
                             tr_y_zzzz_xxyzzz,  \
                             tr_y_zzzz_xxzzz,   \
                             tr_y_zzzz_xxzzzz,  \
                             tr_y_zzzz_xyyyy,   \
                             tr_y_zzzz_xyyyyy,  \
                             tr_y_zzzz_xyyyyz,  \
                             tr_y_zzzz_xyyyz,   \
                             tr_y_zzzz_xyyyzz,  \
                             tr_y_zzzz_xyyzz,   \
                             tr_y_zzzz_xyyzzz,  \
                             tr_y_zzzz_xyzzz,   \
                             tr_y_zzzz_xyzzzz,  \
                             tr_y_zzzz_xzzzz,   \
                             tr_y_zzzz_xzzzzz,  \
                             tr_y_zzzz_yyyyy,   \
                             tr_y_zzzz_yyyyyy,  \
                             tr_y_zzzz_yyyyyz,  \
                             tr_y_zzzz_yyyyz,   \
                             tr_y_zzzz_yyyyzz,  \
                             tr_y_zzzz_yyyzz,   \
                             tr_y_zzzz_yyyzzz,  \
                             tr_y_zzzz_yyzzz,   \
                             tr_y_zzzz_yyzzzz,  \
                             tr_y_zzzz_yzzzz,   \
                             tr_y_zzzz_yzzzzz,  \
                             tr_y_zzzz_zzzzz,   \
                             tr_y_zzzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzz_xxxxxx[i] = 6.0 * tr_y_zzzz_xxxxx[i] * fe_0 + tr_y_zzzz_xxxxxx[i] * pa_x[i];

        tr_y_xzzzz_xxxxxy[i] = 5.0 * tr_y_zzzz_xxxxy[i] * fe_0 + tr_y_zzzz_xxxxxy[i] * pa_x[i];

        tr_y_xzzzz_xxxxxz[i] = 5.0 * tr_y_zzzz_xxxxz[i] * fe_0 + tr_y_zzzz_xxxxxz[i] * pa_x[i];

        tr_y_xzzzz_xxxxyy[i] = 4.0 * tr_y_zzzz_xxxyy[i] * fe_0 + tr_y_zzzz_xxxxyy[i] * pa_x[i];

        tr_y_xzzzz_xxxxyz[i] = 4.0 * tr_y_zzzz_xxxyz[i] * fe_0 + tr_y_zzzz_xxxxyz[i] * pa_x[i];

        tr_y_xzzzz_xxxxzz[i] = 4.0 * tr_y_zzzz_xxxzz[i] * fe_0 + tr_y_zzzz_xxxxzz[i] * pa_x[i];

        tr_y_xzzzz_xxxyyy[i] = 3.0 * tr_y_zzzz_xxyyy[i] * fe_0 + tr_y_zzzz_xxxyyy[i] * pa_x[i];

        tr_y_xzzzz_xxxyyz[i] = 3.0 * tr_y_zzzz_xxyyz[i] * fe_0 + tr_y_zzzz_xxxyyz[i] * pa_x[i];

        tr_y_xzzzz_xxxyzz[i] = 3.0 * tr_y_zzzz_xxyzz[i] * fe_0 + tr_y_zzzz_xxxyzz[i] * pa_x[i];

        tr_y_xzzzz_xxxzzz[i] = 3.0 * tr_y_zzzz_xxzzz[i] * fe_0 + tr_y_zzzz_xxxzzz[i] * pa_x[i];

        tr_y_xzzzz_xxyyyy[i] = 2.0 * tr_y_zzzz_xyyyy[i] * fe_0 + tr_y_zzzz_xxyyyy[i] * pa_x[i];

        tr_y_xzzzz_xxyyyz[i] = 2.0 * tr_y_zzzz_xyyyz[i] * fe_0 + tr_y_zzzz_xxyyyz[i] * pa_x[i];

        tr_y_xzzzz_xxyyzz[i] = 2.0 * tr_y_zzzz_xyyzz[i] * fe_0 + tr_y_zzzz_xxyyzz[i] * pa_x[i];

        tr_y_xzzzz_xxyzzz[i] = 2.0 * tr_y_zzzz_xyzzz[i] * fe_0 + tr_y_zzzz_xxyzzz[i] * pa_x[i];

        tr_y_xzzzz_xxzzzz[i] = 2.0 * tr_y_zzzz_xzzzz[i] * fe_0 + tr_y_zzzz_xxzzzz[i] * pa_x[i];

        tr_y_xzzzz_xyyyyy[i] = tr_y_zzzz_yyyyy[i] * fe_0 + tr_y_zzzz_xyyyyy[i] * pa_x[i];

        tr_y_xzzzz_xyyyyz[i] = tr_y_zzzz_yyyyz[i] * fe_0 + tr_y_zzzz_xyyyyz[i] * pa_x[i];

        tr_y_xzzzz_xyyyzz[i] = tr_y_zzzz_yyyzz[i] * fe_0 + tr_y_zzzz_xyyyzz[i] * pa_x[i];

        tr_y_xzzzz_xyyzzz[i] = tr_y_zzzz_yyzzz[i] * fe_0 + tr_y_zzzz_xyyzzz[i] * pa_x[i];

        tr_y_xzzzz_xyzzzz[i] = tr_y_zzzz_yzzzz[i] * fe_0 + tr_y_zzzz_xyzzzz[i] * pa_x[i];

        tr_y_xzzzz_xzzzzz[i] = tr_y_zzzz_zzzzz[i] * fe_0 + tr_y_zzzz_xzzzzz[i] * pa_x[i];

        tr_y_xzzzz_yyyyyy[i] = tr_y_zzzz_yyyyyy[i] * pa_x[i];

        tr_y_xzzzz_yyyyyz[i] = tr_y_zzzz_yyyyyz[i] * pa_x[i];

        tr_y_xzzzz_yyyyzz[i] = tr_y_zzzz_yyyyzz[i] * pa_x[i];

        tr_y_xzzzz_yyyzzz[i] = tr_y_zzzz_yyyzzz[i] * pa_x[i];

        tr_y_xzzzz_yyzzzz[i] = tr_y_zzzz_yyzzzz[i] * pa_x[i];

        tr_y_xzzzz_yzzzzz[i] = tr_y_zzzz_yzzzzz[i] * pa_x[i];

        tr_y_xzzzz_zzzzzz[i] = tr_y_zzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1008-1036 components of targeted buffer : HI

    auto tr_y_yyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 1008);

    auto tr_y_yyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1009);

    auto tr_y_yyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 1010);

    auto tr_y_yyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1011);

    auto tr_y_yyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1012);

    auto tr_y_yyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 1013);

    auto tr_y_yyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1014);

    auto tr_y_yyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1015);

    auto tr_y_yyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1016);

    auto tr_y_yyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 1017);

    auto tr_y_yyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1018);

    auto tr_y_yyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1019);

    auto tr_y_yyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1020);

    auto tr_y_yyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1021);

    auto tr_y_yyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 1022);

    auto tr_y_yyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1023);

    auto tr_y_yyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1024);

    auto tr_y_yyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1025);

    auto tr_y_yyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1026);

    auto tr_y_yyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1027);

    auto tr_y_yyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 1028);

    auto tr_y_yyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1029);

    auto tr_y_yyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1030);

    auto tr_y_yyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1031);

    auto tr_y_yyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1032);

    auto tr_y_yyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1033);

    auto tr_y_yyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1034);

    auto tr_y_yyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1035);

#pragma omp simd aligned(pa_y,                  \
                             tr_y_yyy_xxxxxx,   \
                             tr_y_yyy_xxxxxy,   \
                             tr_y_yyy_xxxxxz,   \
                             tr_y_yyy_xxxxyy,   \
                             tr_y_yyy_xxxxyz,   \
                             tr_y_yyy_xxxxzz,   \
                             tr_y_yyy_xxxyyy,   \
                             tr_y_yyy_xxxyyz,   \
                             tr_y_yyy_xxxyzz,   \
                             tr_y_yyy_xxxzzz,   \
                             tr_y_yyy_xxyyyy,   \
                             tr_y_yyy_xxyyyz,   \
                             tr_y_yyy_xxyyzz,   \
                             tr_y_yyy_xxyzzz,   \
                             tr_y_yyy_xxzzzz,   \
                             tr_y_yyy_xyyyyy,   \
                             tr_y_yyy_xyyyyz,   \
                             tr_y_yyy_xyyyzz,   \
                             tr_y_yyy_xyyzzz,   \
                             tr_y_yyy_xyzzzz,   \
                             tr_y_yyy_xzzzzz,   \
                             tr_y_yyy_yyyyyy,   \
                             tr_y_yyy_yyyyyz,   \
                             tr_y_yyy_yyyyzz,   \
                             tr_y_yyy_yyyzzz,   \
                             tr_y_yyy_yyzzzz,   \
                             tr_y_yyy_yzzzzz,   \
                             tr_y_yyy_zzzzzz,   \
                             tr_y_yyyy_xxxxx,   \
                             tr_y_yyyy_xxxxxx,  \
                             tr_y_yyyy_xxxxxy,  \
                             tr_y_yyyy_xxxxxz,  \
                             tr_y_yyyy_xxxxy,   \
                             tr_y_yyyy_xxxxyy,  \
                             tr_y_yyyy_xxxxyz,  \
                             tr_y_yyyy_xxxxz,   \
                             tr_y_yyyy_xxxxzz,  \
                             tr_y_yyyy_xxxyy,   \
                             tr_y_yyyy_xxxyyy,  \
                             tr_y_yyyy_xxxyyz,  \
                             tr_y_yyyy_xxxyz,   \
                             tr_y_yyyy_xxxyzz,  \
                             tr_y_yyyy_xxxzz,   \
                             tr_y_yyyy_xxxzzz,  \
                             tr_y_yyyy_xxyyy,   \
                             tr_y_yyyy_xxyyyy,  \
                             tr_y_yyyy_xxyyyz,  \
                             tr_y_yyyy_xxyyz,   \
                             tr_y_yyyy_xxyyzz,  \
                             tr_y_yyyy_xxyzz,   \
                             tr_y_yyyy_xxyzzz,  \
                             tr_y_yyyy_xxzzz,   \
                             tr_y_yyyy_xxzzzz,  \
                             tr_y_yyyy_xyyyy,   \
                             tr_y_yyyy_xyyyyy,  \
                             tr_y_yyyy_xyyyyz,  \
                             tr_y_yyyy_xyyyz,   \
                             tr_y_yyyy_xyyyzz,  \
                             tr_y_yyyy_xyyzz,   \
                             tr_y_yyyy_xyyzzz,  \
                             tr_y_yyyy_xyzzz,   \
                             tr_y_yyyy_xyzzzz,  \
                             tr_y_yyyy_xzzzz,   \
                             tr_y_yyyy_xzzzzz,  \
                             tr_y_yyyy_yyyyy,   \
                             tr_y_yyyy_yyyyyy,  \
                             tr_y_yyyy_yyyyyz,  \
                             tr_y_yyyy_yyyyz,   \
                             tr_y_yyyy_yyyyzz,  \
                             tr_y_yyyy_yyyzz,   \
                             tr_y_yyyy_yyyzzz,  \
                             tr_y_yyyy_yyzzz,   \
                             tr_y_yyyy_yyzzzz,  \
                             tr_y_yyyy_yzzzz,   \
                             tr_y_yyyy_yzzzzz,  \
                             tr_y_yyyy_zzzzz,   \
                             tr_y_yyyy_zzzzzz,  \
                             tr_y_yyyyy_xxxxxx, \
                             tr_y_yyyyy_xxxxxy, \
                             tr_y_yyyyy_xxxxxz, \
                             tr_y_yyyyy_xxxxyy, \
                             tr_y_yyyyy_xxxxyz, \
                             tr_y_yyyyy_xxxxzz, \
                             tr_y_yyyyy_xxxyyy, \
                             tr_y_yyyyy_xxxyyz, \
                             tr_y_yyyyy_xxxyzz, \
                             tr_y_yyyyy_xxxzzz, \
                             tr_y_yyyyy_xxyyyy, \
                             tr_y_yyyyy_xxyyyz, \
                             tr_y_yyyyy_xxyyzz, \
                             tr_y_yyyyy_xxyzzz, \
                             tr_y_yyyyy_xxzzzz, \
                             tr_y_yyyyy_xyyyyy, \
                             tr_y_yyyyy_xyyyyz, \
                             tr_y_yyyyy_xyyyzz, \
                             tr_y_yyyyy_xyyzzz, \
                             tr_y_yyyyy_xyzzzz, \
                             tr_y_yyyyy_xzzzzz, \
                             tr_y_yyyyy_yyyyyy, \
                             tr_y_yyyyy_yyyyyz, \
                             tr_y_yyyyy_yyyyzz, \
                             tr_y_yyyyy_yyyzzz, \
                             tr_y_yyyyy_yyzzzz, \
                             tr_y_yyyyy_yzzzzz, \
                             tr_y_yyyyy_zzzzzz, \
                             ts_yyyy_xxxxxx,    \
                             ts_yyyy_xxxxxy,    \
                             ts_yyyy_xxxxxz,    \
                             ts_yyyy_xxxxyy,    \
                             ts_yyyy_xxxxyz,    \
                             ts_yyyy_xxxxzz,    \
                             ts_yyyy_xxxyyy,    \
                             ts_yyyy_xxxyyz,    \
                             ts_yyyy_xxxyzz,    \
                             ts_yyyy_xxxzzz,    \
                             ts_yyyy_xxyyyy,    \
                             ts_yyyy_xxyyyz,    \
                             ts_yyyy_xxyyzz,    \
                             ts_yyyy_xxyzzz,    \
                             ts_yyyy_xxzzzz,    \
                             ts_yyyy_xyyyyy,    \
                             ts_yyyy_xyyyyz,    \
                             ts_yyyy_xyyyzz,    \
                             ts_yyyy_xyyzzz,    \
                             ts_yyyy_xyzzzz,    \
                             ts_yyyy_xzzzzz,    \
                             ts_yyyy_yyyyyy,    \
                             ts_yyyy_yyyyyz,    \
                             ts_yyyy_yyyyzz,    \
                             ts_yyyy_yyyzzz,    \
                             ts_yyyy_yyzzzz,    \
                             ts_yyyy_yzzzzz,    \
                             ts_yyyy_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyy_xxxxxx[i] = 4.0 * tr_y_yyy_xxxxxx[i] * fe_0 + ts_yyyy_xxxxxx[i] * fe_0 + tr_y_yyyy_xxxxxx[i] * pa_y[i];

        tr_y_yyyyy_xxxxxy[i] = 4.0 * tr_y_yyy_xxxxxy[i] * fe_0 + tr_y_yyyy_xxxxx[i] * fe_0 + ts_yyyy_xxxxxy[i] * fe_0 + tr_y_yyyy_xxxxxy[i] * pa_y[i];

        tr_y_yyyyy_xxxxxz[i] = 4.0 * tr_y_yyy_xxxxxz[i] * fe_0 + ts_yyyy_xxxxxz[i] * fe_0 + tr_y_yyyy_xxxxxz[i] * pa_y[i];

        tr_y_yyyyy_xxxxyy[i] =
            4.0 * tr_y_yyy_xxxxyy[i] * fe_0 + 2.0 * tr_y_yyyy_xxxxy[i] * fe_0 + ts_yyyy_xxxxyy[i] * fe_0 + tr_y_yyyy_xxxxyy[i] * pa_y[i];

        tr_y_yyyyy_xxxxyz[i] = 4.0 * tr_y_yyy_xxxxyz[i] * fe_0 + tr_y_yyyy_xxxxz[i] * fe_0 + ts_yyyy_xxxxyz[i] * fe_0 + tr_y_yyyy_xxxxyz[i] * pa_y[i];

        tr_y_yyyyy_xxxxzz[i] = 4.0 * tr_y_yyy_xxxxzz[i] * fe_0 + ts_yyyy_xxxxzz[i] * fe_0 + tr_y_yyyy_xxxxzz[i] * pa_y[i];

        tr_y_yyyyy_xxxyyy[i] =
            4.0 * tr_y_yyy_xxxyyy[i] * fe_0 + 3.0 * tr_y_yyyy_xxxyy[i] * fe_0 + ts_yyyy_xxxyyy[i] * fe_0 + tr_y_yyyy_xxxyyy[i] * pa_y[i];

        tr_y_yyyyy_xxxyyz[i] =
            4.0 * tr_y_yyy_xxxyyz[i] * fe_0 + 2.0 * tr_y_yyyy_xxxyz[i] * fe_0 + ts_yyyy_xxxyyz[i] * fe_0 + tr_y_yyyy_xxxyyz[i] * pa_y[i];

        tr_y_yyyyy_xxxyzz[i] = 4.0 * tr_y_yyy_xxxyzz[i] * fe_0 + tr_y_yyyy_xxxzz[i] * fe_0 + ts_yyyy_xxxyzz[i] * fe_0 + tr_y_yyyy_xxxyzz[i] * pa_y[i];

        tr_y_yyyyy_xxxzzz[i] = 4.0 * tr_y_yyy_xxxzzz[i] * fe_0 + ts_yyyy_xxxzzz[i] * fe_0 + tr_y_yyyy_xxxzzz[i] * pa_y[i];

        tr_y_yyyyy_xxyyyy[i] =
            4.0 * tr_y_yyy_xxyyyy[i] * fe_0 + 4.0 * tr_y_yyyy_xxyyy[i] * fe_0 + ts_yyyy_xxyyyy[i] * fe_0 + tr_y_yyyy_xxyyyy[i] * pa_y[i];

        tr_y_yyyyy_xxyyyz[i] =
            4.0 * tr_y_yyy_xxyyyz[i] * fe_0 + 3.0 * tr_y_yyyy_xxyyz[i] * fe_0 + ts_yyyy_xxyyyz[i] * fe_0 + tr_y_yyyy_xxyyyz[i] * pa_y[i];

        tr_y_yyyyy_xxyyzz[i] =
            4.0 * tr_y_yyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyyy_xxyzz[i] * fe_0 + ts_yyyy_xxyyzz[i] * fe_0 + tr_y_yyyy_xxyyzz[i] * pa_y[i];

        tr_y_yyyyy_xxyzzz[i] = 4.0 * tr_y_yyy_xxyzzz[i] * fe_0 + tr_y_yyyy_xxzzz[i] * fe_0 + ts_yyyy_xxyzzz[i] * fe_0 + tr_y_yyyy_xxyzzz[i] * pa_y[i];

        tr_y_yyyyy_xxzzzz[i] = 4.0 * tr_y_yyy_xxzzzz[i] * fe_0 + ts_yyyy_xxzzzz[i] * fe_0 + tr_y_yyyy_xxzzzz[i] * pa_y[i];

        tr_y_yyyyy_xyyyyy[i] =
            4.0 * tr_y_yyy_xyyyyy[i] * fe_0 + 5.0 * tr_y_yyyy_xyyyy[i] * fe_0 + ts_yyyy_xyyyyy[i] * fe_0 + tr_y_yyyy_xyyyyy[i] * pa_y[i];

        tr_y_yyyyy_xyyyyz[i] =
            4.0 * tr_y_yyy_xyyyyz[i] * fe_0 + 4.0 * tr_y_yyyy_xyyyz[i] * fe_0 + ts_yyyy_xyyyyz[i] * fe_0 + tr_y_yyyy_xyyyyz[i] * pa_y[i];

        tr_y_yyyyy_xyyyzz[i] =
            4.0 * tr_y_yyy_xyyyzz[i] * fe_0 + 3.0 * tr_y_yyyy_xyyzz[i] * fe_0 + ts_yyyy_xyyyzz[i] * fe_0 + tr_y_yyyy_xyyyzz[i] * pa_y[i];

        tr_y_yyyyy_xyyzzz[i] =
            4.0 * tr_y_yyy_xyyzzz[i] * fe_0 + 2.0 * tr_y_yyyy_xyzzz[i] * fe_0 + ts_yyyy_xyyzzz[i] * fe_0 + tr_y_yyyy_xyyzzz[i] * pa_y[i];

        tr_y_yyyyy_xyzzzz[i] = 4.0 * tr_y_yyy_xyzzzz[i] * fe_0 + tr_y_yyyy_xzzzz[i] * fe_0 + ts_yyyy_xyzzzz[i] * fe_0 + tr_y_yyyy_xyzzzz[i] * pa_y[i];

        tr_y_yyyyy_xzzzzz[i] = 4.0 * tr_y_yyy_xzzzzz[i] * fe_0 + ts_yyyy_xzzzzz[i] * fe_0 + tr_y_yyyy_xzzzzz[i] * pa_y[i];

        tr_y_yyyyy_yyyyyy[i] =
            4.0 * tr_y_yyy_yyyyyy[i] * fe_0 + 6.0 * tr_y_yyyy_yyyyy[i] * fe_0 + ts_yyyy_yyyyyy[i] * fe_0 + tr_y_yyyy_yyyyyy[i] * pa_y[i];

        tr_y_yyyyy_yyyyyz[i] =
            4.0 * tr_y_yyy_yyyyyz[i] * fe_0 + 5.0 * tr_y_yyyy_yyyyz[i] * fe_0 + ts_yyyy_yyyyyz[i] * fe_0 + tr_y_yyyy_yyyyyz[i] * pa_y[i];

        tr_y_yyyyy_yyyyzz[i] =
            4.0 * tr_y_yyy_yyyyzz[i] * fe_0 + 4.0 * tr_y_yyyy_yyyzz[i] * fe_0 + ts_yyyy_yyyyzz[i] * fe_0 + tr_y_yyyy_yyyyzz[i] * pa_y[i];

        tr_y_yyyyy_yyyzzz[i] =
            4.0 * tr_y_yyy_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyyy_yyzzz[i] * fe_0 + ts_yyyy_yyyzzz[i] * fe_0 + tr_y_yyyy_yyyzzz[i] * pa_y[i];

        tr_y_yyyyy_yyzzzz[i] =
            4.0 * tr_y_yyy_yyzzzz[i] * fe_0 + 2.0 * tr_y_yyyy_yzzzz[i] * fe_0 + ts_yyyy_yyzzzz[i] * fe_0 + tr_y_yyyy_yyzzzz[i] * pa_y[i];

        tr_y_yyyyy_yzzzzz[i] = 4.0 * tr_y_yyy_yzzzzz[i] * fe_0 + tr_y_yyyy_zzzzz[i] * fe_0 + ts_yyyy_yzzzzz[i] * fe_0 + tr_y_yyyy_yzzzzz[i] * pa_y[i];

        tr_y_yyyyy_zzzzzz[i] = 4.0 * tr_y_yyy_zzzzzz[i] * fe_0 + ts_yyyy_zzzzzz[i] * fe_0 + tr_y_yyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 1036-1064 components of targeted buffer : HI

    auto tr_y_yyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 1036);

    auto tr_y_yyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 1037);

    auto tr_y_yyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 1038);

    auto tr_y_yyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 1039);

    auto tr_y_yyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1040);

    auto tr_y_yyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 1041);

    auto tr_y_yyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 1042);

    auto tr_y_yyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1043);

    auto tr_y_yyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1044);

    auto tr_y_yyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 1045);

    auto tr_y_yyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 1046);

    auto tr_y_yyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1047);

    auto tr_y_yyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1048);

    auto tr_y_yyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1049);

    auto tr_y_yyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 1050);

    auto tr_y_yyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 1051);

    auto tr_y_yyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1052);

    auto tr_y_yyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1053);

    auto tr_y_yyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1054);

    auto tr_y_yyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1055);

    auto tr_y_yyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 1056);

    auto tr_y_yyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1057);

    auto tr_y_yyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1058);

    auto tr_y_yyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1059);

    auto tr_y_yyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1060);

    auto tr_y_yyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1061);

    auto tr_y_yyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1062);

    auto tr_y_yyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1063);

#pragma omp simd aligned(pa_z,                  \
                             tr_y_yyyy_xxxxx,   \
                             tr_y_yyyy_xxxxxx,  \
                             tr_y_yyyy_xxxxxy,  \
                             tr_y_yyyy_xxxxxz,  \
                             tr_y_yyyy_xxxxy,   \
                             tr_y_yyyy_xxxxyy,  \
                             tr_y_yyyy_xxxxyz,  \
                             tr_y_yyyy_xxxxz,   \
                             tr_y_yyyy_xxxxzz,  \
                             tr_y_yyyy_xxxyy,   \
                             tr_y_yyyy_xxxyyy,  \
                             tr_y_yyyy_xxxyyz,  \
                             tr_y_yyyy_xxxyz,   \
                             tr_y_yyyy_xxxyzz,  \
                             tr_y_yyyy_xxxzz,   \
                             tr_y_yyyy_xxxzzz,  \
                             tr_y_yyyy_xxyyy,   \
                             tr_y_yyyy_xxyyyy,  \
                             tr_y_yyyy_xxyyyz,  \
                             tr_y_yyyy_xxyyz,   \
                             tr_y_yyyy_xxyyzz,  \
                             tr_y_yyyy_xxyzz,   \
                             tr_y_yyyy_xxyzzz,  \
                             tr_y_yyyy_xxzzz,   \
                             tr_y_yyyy_xxzzzz,  \
                             tr_y_yyyy_xyyyy,   \
                             tr_y_yyyy_xyyyyy,  \
                             tr_y_yyyy_xyyyyz,  \
                             tr_y_yyyy_xyyyz,   \
                             tr_y_yyyy_xyyyzz,  \
                             tr_y_yyyy_xyyzz,   \
                             tr_y_yyyy_xyyzzz,  \
                             tr_y_yyyy_xyzzz,   \
                             tr_y_yyyy_xyzzzz,  \
                             tr_y_yyyy_xzzzz,   \
                             tr_y_yyyy_xzzzzz,  \
                             tr_y_yyyy_yyyyy,   \
                             tr_y_yyyy_yyyyyy,  \
                             tr_y_yyyy_yyyyyz,  \
                             tr_y_yyyy_yyyyz,   \
                             tr_y_yyyy_yyyyzz,  \
                             tr_y_yyyy_yyyzz,   \
                             tr_y_yyyy_yyyzzz,  \
                             tr_y_yyyy_yyzzz,   \
                             tr_y_yyyy_yyzzzz,  \
                             tr_y_yyyy_yzzzz,   \
                             tr_y_yyyy_yzzzzz,  \
                             tr_y_yyyy_zzzzz,   \
                             tr_y_yyyy_zzzzzz,  \
                             tr_y_yyyyz_xxxxxx, \
                             tr_y_yyyyz_xxxxxy, \
                             tr_y_yyyyz_xxxxxz, \
                             tr_y_yyyyz_xxxxyy, \
                             tr_y_yyyyz_xxxxyz, \
                             tr_y_yyyyz_xxxxzz, \
                             tr_y_yyyyz_xxxyyy, \
                             tr_y_yyyyz_xxxyyz, \
                             tr_y_yyyyz_xxxyzz, \
                             tr_y_yyyyz_xxxzzz, \
                             tr_y_yyyyz_xxyyyy, \
                             tr_y_yyyyz_xxyyyz, \
                             tr_y_yyyyz_xxyyzz, \
                             tr_y_yyyyz_xxyzzz, \
                             tr_y_yyyyz_xxzzzz, \
                             tr_y_yyyyz_xyyyyy, \
                             tr_y_yyyyz_xyyyyz, \
                             tr_y_yyyyz_xyyyzz, \
                             tr_y_yyyyz_xyyzzz, \
                             tr_y_yyyyz_xyzzzz, \
                             tr_y_yyyyz_xzzzzz, \
                             tr_y_yyyyz_yyyyyy, \
                             tr_y_yyyyz_yyyyyz, \
                             tr_y_yyyyz_yyyyzz, \
                             tr_y_yyyyz_yyyzzz, \
                             tr_y_yyyyz_yyzzzz, \
                             tr_y_yyyyz_yzzzzz, \
                             tr_y_yyyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyz_xxxxxx[i] = tr_y_yyyy_xxxxxx[i] * pa_z[i];

        tr_y_yyyyz_xxxxxy[i] = tr_y_yyyy_xxxxxy[i] * pa_z[i];

        tr_y_yyyyz_xxxxxz[i] = tr_y_yyyy_xxxxx[i] * fe_0 + tr_y_yyyy_xxxxxz[i] * pa_z[i];

        tr_y_yyyyz_xxxxyy[i] = tr_y_yyyy_xxxxyy[i] * pa_z[i];

        tr_y_yyyyz_xxxxyz[i] = tr_y_yyyy_xxxxy[i] * fe_0 + tr_y_yyyy_xxxxyz[i] * pa_z[i];

        tr_y_yyyyz_xxxxzz[i] = 2.0 * tr_y_yyyy_xxxxz[i] * fe_0 + tr_y_yyyy_xxxxzz[i] * pa_z[i];

        tr_y_yyyyz_xxxyyy[i] = tr_y_yyyy_xxxyyy[i] * pa_z[i];

        tr_y_yyyyz_xxxyyz[i] = tr_y_yyyy_xxxyy[i] * fe_0 + tr_y_yyyy_xxxyyz[i] * pa_z[i];

        tr_y_yyyyz_xxxyzz[i] = 2.0 * tr_y_yyyy_xxxyz[i] * fe_0 + tr_y_yyyy_xxxyzz[i] * pa_z[i];

        tr_y_yyyyz_xxxzzz[i] = 3.0 * tr_y_yyyy_xxxzz[i] * fe_0 + tr_y_yyyy_xxxzzz[i] * pa_z[i];

        tr_y_yyyyz_xxyyyy[i] = tr_y_yyyy_xxyyyy[i] * pa_z[i];

        tr_y_yyyyz_xxyyyz[i] = tr_y_yyyy_xxyyy[i] * fe_0 + tr_y_yyyy_xxyyyz[i] * pa_z[i];

        tr_y_yyyyz_xxyyzz[i] = 2.0 * tr_y_yyyy_xxyyz[i] * fe_0 + tr_y_yyyy_xxyyzz[i] * pa_z[i];

        tr_y_yyyyz_xxyzzz[i] = 3.0 * tr_y_yyyy_xxyzz[i] * fe_0 + tr_y_yyyy_xxyzzz[i] * pa_z[i];

        tr_y_yyyyz_xxzzzz[i] = 4.0 * tr_y_yyyy_xxzzz[i] * fe_0 + tr_y_yyyy_xxzzzz[i] * pa_z[i];

        tr_y_yyyyz_xyyyyy[i] = tr_y_yyyy_xyyyyy[i] * pa_z[i];

        tr_y_yyyyz_xyyyyz[i] = tr_y_yyyy_xyyyy[i] * fe_0 + tr_y_yyyy_xyyyyz[i] * pa_z[i];

        tr_y_yyyyz_xyyyzz[i] = 2.0 * tr_y_yyyy_xyyyz[i] * fe_0 + tr_y_yyyy_xyyyzz[i] * pa_z[i];

        tr_y_yyyyz_xyyzzz[i] = 3.0 * tr_y_yyyy_xyyzz[i] * fe_0 + tr_y_yyyy_xyyzzz[i] * pa_z[i];

        tr_y_yyyyz_xyzzzz[i] = 4.0 * tr_y_yyyy_xyzzz[i] * fe_0 + tr_y_yyyy_xyzzzz[i] * pa_z[i];

        tr_y_yyyyz_xzzzzz[i] = 5.0 * tr_y_yyyy_xzzzz[i] * fe_0 + tr_y_yyyy_xzzzzz[i] * pa_z[i];

        tr_y_yyyyz_yyyyyy[i] = tr_y_yyyy_yyyyyy[i] * pa_z[i];

        tr_y_yyyyz_yyyyyz[i] = tr_y_yyyy_yyyyy[i] * fe_0 + tr_y_yyyy_yyyyyz[i] * pa_z[i];

        tr_y_yyyyz_yyyyzz[i] = 2.0 * tr_y_yyyy_yyyyz[i] * fe_0 + tr_y_yyyy_yyyyzz[i] * pa_z[i];

        tr_y_yyyyz_yyyzzz[i] = 3.0 * tr_y_yyyy_yyyzz[i] * fe_0 + tr_y_yyyy_yyyzzz[i] * pa_z[i];

        tr_y_yyyyz_yyzzzz[i] = 4.0 * tr_y_yyyy_yyzzz[i] * fe_0 + tr_y_yyyy_yyzzzz[i] * pa_z[i];

        tr_y_yyyyz_yzzzzz[i] = 5.0 * tr_y_yyyy_yzzzz[i] * fe_0 + tr_y_yyyy_yzzzzz[i] * pa_z[i];

        tr_y_yyyyz_zzzzzz[i] = 6.0 * tr_y_yyyy_zzzzz[i] * fe_0 + tr_y_yyyy_zzzzzz[i] * pa_z[i];
    }

    // Set up 1064-1092 components of targeted buffer : HI

    auto tr_y_yyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 1064);

    auto tr_y_yyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 1065);

    auto tr_y_yyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 1066);

    auto tr_y_yyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 1067);

    auto tr_y_yyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 1068);

    auto tr_y_yyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 1069);

    auto tr_y_yyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 1070);

    auto tr_y_yyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 1071);

    auto tr_y_yyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 1072);

    auto tr_y_yyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 1073);

    auto tr_y_yyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 1074);

    auto tr_y_yyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 1075);

    auto tr_y_yyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 1076);

    auto tr_y_yyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 1077);

    auto tr_y_yyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 1078);

    auto tr_y_yyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 1079);

    auto tr_y_yyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 1080);

    auto tr_y_yyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 1081);

    auto tr_y_yyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 1082);

    auto tr_y_yyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 1083);

    auto tr_y_yyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 1084);

    auto tr_y_yyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1085);

    auto tr_y_yyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1086);

    auto tr_y_yyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1087);

    auto tr_y_yyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1088);

    auto tr_y_yyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1089);

    auto tr_y_yyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1090);

    auto tr_y_yyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 1091);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_y_yyy_xxxxxx,   \
                             tr_y_yyy_xxxxxy,   \
                             tr_y_yyy_xxxxyy,   \
                             tr_y_yyy_xxxxyz,   \
                             tr_y_yyy_xxxyyy,   \
                             tr_y_yyy_xxxyyz,   \
                             tr_y_yyy_xxxyzz,   \
                             tr_y_yyy_xxyyyy,   \
                             tr_y_yyy_xxyyyz,   \
                             tr_y_yyy_xxyyzz,   \
                             tr_y_yyy_xxyzzz,   \
                             tr_y_yyy_xyyyyy,   \
                             tr_y_yyy_xyyyyz,   \
                             tr_y_yyy_xyyyzz,   \
                             tr_y_yyy_xyyzzz,   \
                             tr_y_yyy_xyzzzz,   \
                             tr_y_yyy_yyyyyy,   \
                             tr_y_yyy_yyyyyz,   \
                             tr_y_yyy_yyyyzz,   \
                             tr_y_yyy_yyyzzz,   \
                             tr_y_yyy_yyzzzz,   \
                             tr_y_yyy_yzzzzz,   \
                             tr_y_yyyz_xxxxxx,  \
                             tr_y_yyyz_xxxxxy,  \
                             tr_y_yyyz_xxxxy,   \
                             tr_y_yyyz_xxxxyy,  \
                             tr_y_yyyz_xxxxyz,  \
                             tr_y_yyyz_xxxyy,   \
                             tr_y_yyyz_xxxyyy,  \
                             tr_y_yyyz_xxxyyz,  \
                             tr_y_yyyz_xxxyz,   \
                             tr_y_yyyz_xxxyzz,  \
                             tr_y_yyyz_xxyyy,   \
                             tr_y_yyyz_xxyyyy,  \
                             tr_y_yyyz_xxyyyz,  \
                             tr_y_yyyz_xxyyz,   \
                             tr_y_yyyz_xxyyzz,  \
                             tr_y_yyyz_xxyzz,   \
                             tr_y_yyyz_xxyzzz,  \
                             tr_y_yyyz_xyyyy,   \
                             tr_y_yyyz_xyyyyy,  \
                             tr_y_yyyz_xyyyyz,  \
                             tr_y_yyyz_xyyyz,   \
                             tr_y_yyyz_xyyyzz,  \
                             tr_y_yyyz_xyyzz,   \
                             tr_y_yyyz_xyyzzz,  \
                             tr_y_yyyz_xyzzz,   \
                             tr_y_yyyz_xyzzzz,  \
                             tr_y_yyyz_yyyyy,   \
                             tr_y_yyyz_yyyyyy,  \
                             tr_y_yyyz_yyyyyz,  \
                             tr_y_yyyz_yyyyz,   \
                             tr_y_yyyz_yyyyzz,  \
                             tr_y_yyyz_yyyzz,   \
                             tr_y_yyyz_yyyzzz,  \
                             tr_y_yyyz_yyzzz,   \
                             tr_y_yyyz_yyzzzz,  \
                             tr_y_yyyz_yzzzz,   \
                             tr_y_yyyz_yzzzzz,  \
                             tr_y_yyyzz_xxxxxx, \
                             tr_y_yyyzz_xxxxxy, \
                             tr_y_yyyzz_xxxxxz, \
                             tr_y_yyyzz_xxxxyy, \
                             tr_y_yyyzz_xxxxyz, \
                             tr_y_yyyzz_xxxxzz, \
                             tr_y_yyyzz_xxxyyy, \
                             tr_y_yyyzz_xxxyyz, \
                             tr_y_yyyzz_xxxyzz, \
                             tr_y_yyyzz_xxxzzz, \
                             tr_y_yyyzz_xxyyyy, \
                             tr_y_yyyzz_xxyyyz, \
                             tr_y_yyyzz_xxyyzz, \
                             tr_y_yyyzz_xxyzzz, \
                             tr_y_yyyzz_xxzzzz, \
                             tr_y_yyyzz_xyyyyy, \
                             tr_y_yyyzz_xyyyyz, \
                             tr_y_yyyzz_xyyyzz, \
                             tr_y_yyyzz_xyyzzz, \
                             tr_y_yyyzz_xyzzzz, \
                             tr_y_yyyzz_xzzzzz, \
                             tr_y_yyyzz_yyyyyy, \
                             tr_y_yyyzz_yyyyyz, \
                             tr_y_yyyzz_yyyyzz, \
                             tr_y_yyyzz_yyyzzz, \
                             tr_y_yyyzz_yyzzzz, \
                             tr_y_yyyzz_yzzzzz, \
                             tr_y_yyyzz_zzzzzz, \
                             tr_y_yyzz_xxxxxz,  \
                             tr_y_yyzz_xxxxzz,  \
                             tr_y_yyzz_xxxzzz,  \
                             tr_y_yyzz_xxzzzz,  \
                             tr_y_yyzz_xzzzzz,  \
                             tr_y_yyzz_zzzzzz,  \
                             tr_y_yzz_xxxxxz,   \
                             tr_y_yzz_xxxxzz,   \
                             tr_y_yzz_xxxzzz,   \
                             tr_y_yzz_xxzzzz,   \
                             tr_y_yzz_xzzzzz,   \
                             tr_y_yzz_zzzzzz,   \
                             ts_yyzz_xxxxxz,    \
                             ts_yyzz_xxxxzz,    \
                             ts_yyzz_xxxzzz,    \
                             ts_yyzz_xxzzzz,    \
                             ts_yyzz_xzzzzz,    \
                             ts_yyzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzz_xxxxxx[i] = tr_y_yyy_xxxxxx[i] * fe_0 + tr_y_yyyz_xxxxxx[i] * pa_z[i];

        tr_y_yyyzz_xxxxxy[i] = tr_y_yyy_xxxxxy[i] * fe_0 + tr_y_yyyz_xxxxxy[i] * pa_z[i];

        tr_y_yyyzz_xxxxxz[i] = 2.0 * tr_y_yzz_xxxxxz[i] * fe_0 + ts_yyzz_xxxxxz[i] * fe_0 + tr_y_yyzz_xxxxxz[i] * pa_y[i];

        tr_y_yyyzz_xxxxyy[i] = tr_y_yyy_xxxxyy[i] * fe_0 + tr_y_yyyz_xxxxyy[i] * pa_z[i];

        tr_y_yyyzz_xxxxyz[i] = tr_y_yyy_xxxxyz[i] * fe_0 + tr_y_yyyz_xxxxy[i] * fe_0 + tr_y_yyyz_xxxxyz[i] * pa_z[i];

        tr_y_yyyzz_xxxxzz[i] = 2.0 * tr_y_yzz_xxxxzz[i] * fe_0 + ts_yyzz_xxxxzz[i] * fe_0 + tr_y_yyzz_xxxxzz[i] * pa_y[i];

        tr_y_yyyzz_xxxyyy[i] = tr_y_yyy_xxxyyy[i] * fe_0 + tr_y_yyyz_xxxyyy[i] * pa_z[i];

        tr_y_yyyzz_xxxyyz[i] = tr_y_yyy_xxxyyz[i] * fe_0 + tr_y_yyyz_xxxyy[i] * fe_0 + tr_y_yyyz_xxxyyz[i] * pa_z[i];

        tr_y_yyyzz_xxxyzz[i] = tr_y_yyy_xxxyzz[i] * fe_0 + 2.0 * tr_y_yyyz_xxxyz[i] * fe_0 + tr_y_yyyz_xxxyzz[i] * pa_z[i];

        tr_y_yyyzz_xxxzzz[i] = 2.0 * tr_y_yzz_xxxzzz[i] * fe_0 + ts_yyzz_xxxzzz[i] * fe_0 + tr_y_yyzz_xxxzzz[i] * pa_y[i];

        tr_y_yyyzz_xxyyyy[i] = tr_y_yyy_xxyyyy[i] * fe_0 + tr_y_yyyz_xxyyyy[i] * pa_z[i];

        tr_y_yyyzz_xxyyyz[i] = tr_y_yyy_xxyyyz[i] * fe_0 + tr_y_yyyz_xxyyy[i] * fe_0 + tr_y_yyyz_xxyyyz[i] * pa_z[i];

        tr_y_yyyzz_xxyyzz[i] = tr_y_yyy_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyyz_xxyyz[i] * fe_0 + tr_y_yyyz_xxyyzz[i] * pa_z[i];

        tr_y_yyyzz_xxyzzz[i] = tr_y_yyy_xxyzzz[i] * fe_0 + 3.0 * tr_y_yyyz_xxyzz[i] * fe_0 + tr_y_yyyz_xxyzzz[i] * pa_z[i];

        tr_y_yyyzz_xxzzzz[i] = 2.0 * tr_y_yzz_xxzzzz[i] * fe_0 + ts_yyzz_xxzzzz[i] * fe_0 + tr_y_yyzz_xxzzzz[i] * pa_y[i];

        tr_y_yyyzz_xyyyyy[i] = tr_y_yyy_xyyyyy[i] * fe_0 + tr_y_yyyz_xyyyyy[i] * pa_z[i];

        tr_y_yyyzz_xyyyyz[i] = tr_y_yyy_xyyyyz[i] * fe_0 + tr_y_yyyz_xyyyy[i] * fe_0 + tr_y_yyyz_xyyyyz[i] * pa_z[i];

        tr_y_yyyzz_xyyyzz[i] = tr_y_yyy_xyyyzz[i] * fe_0 + 2.0 * tr_y_yyyz_xyyyz[i] * fe_0 + tr_y_yyyz_xyyyzz[i] * pa_z[i];

        tr_y_yyyzz_xyyzzz[i] = tr_y_yyy_xyyzzz[i] * fe_0 + 3.0 * tr_y_yyyz_xyyzz[i] * fe_0 + tr_y_yyyz_xyyzzz[i] * pa_z[i];

        tr_y_yyyzz_xyzzzz[i] = tr_y_yyy_xyzzzz[i] * fe_0 + 4.0 * tr_y_yyyz_xyzzz[i] * fe_0 + tr_y_yyyz_xyzzzz[i] * pa_z[i];

        tr_y_yyyzz_xzzzzz[i] = 2.0 * tr_y_yzz_xzzzzz[i] * fe_0 + ts_yyzz_xzzzzz[i] * fe_0 + tr_y_yyzz_xzzzzz[i] * pa_y[i];

        tr_y_yyyzz_yyyyyy[i] = tr_y_yyy_yyyyyy[i] * fe_0 + tr_y_yyyz_yyyyyy[i] * pa_z[i];

        tr_y_yyyzz_yyyyyz[i] = tr_y_yyy_yyyyyz[i] * fe_0 + tr_y_yyyz_yyyyy[i] * fe_0 + tr_y_yyyz_yyyyyz[i] * pa_z[i];

        tr_y_yyyzz_yyyyzz[i] = tr_y_yyy_yyyyzz[i] * fe_0 + 2.0 * tr_y_yyyz_yyyyz[i] * fe_0 + tr_y_yyyz_yyyyzz[i] * pa_z[i];

        tr_y_yyyzz_yyyzzz[i] = tr_y_yyy_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyyz_yyyzz[i] * fe_0 + tr_y_yyyz_yyyzzz[i] * pa_z[i];

        tr_y_yyyzz_yyzzzz[i] = tr_y_yyy_yyzzzz[i] * fe_0 + 4.0 * tr_y_yyyz_yyzzz[i] * fe_0 + tr_y_yyyz_yyzzzz[i] * pa_z[i];

        tr_y_yyyzz_yzzzzz[i] = tr_y_yyy_yzzzzz[i] * fe_0 + 5.0 * tr_y_yyyz_yzzzz[i] * fe_0 + tr_y_yyyz_yzzzzz[i] * pa_z[i];

        tr_y_yyyzz_zzzzzz[i] = 2.0 * tr_y_yzz_zzzzzz[i] * fe_0 + ts_yyzz_zzzzzz[i] * fe_0 + tr_y_yyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1092-1120 components of targeted buffer : HI

    auto tr_y_yyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1092);

    auto tr_y_yyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1093);

    auto tr_y_yyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1094);

    auto tr_y_yyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1095);

    auto tr_y_yyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1096);

    auto tr_y_yyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1097);

    auto tr_y_yyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1098);

    auto tr_y_yyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1099);

    auto tr_y_yyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1100);

    auto tr_y_yyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1101);

    auto tr_y_yyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1102);

    auto tr_y_yyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1103);

    auto tr_y_yyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1104);

    auto tr_y_yyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1105);

    auto tr_y_yyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1106);

    auto tr_y_yyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1107);

    auto tr_y_yyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1108);

    auto tr_y_yyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1109);

    auto tr_y_yyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1110);

    auto tr_y_yyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1111);

    auto tr_y_yyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1112);

    auto tr_y_yyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1113);

    auto tr_y_yyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1114);

    auto tr_y_yyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1115);

    auto tr_y_yyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1116);

    auto tr_y_yyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1117);

    auto tr_y_yyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1118);

    auto tr_y_yyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1119);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_y_yyz_xxxxxx,   \
                             tr_y_yyz_xxxxxy,   \
                             tr_y_yyz_xxxxyy,   \
                             tr_y_yyz_xxxxyz,   \
                             tr_y_yyz_xxxyyy,   \
                             tr_y_yyz_xxxyyz,   \
                             tr_y_yyz_xxxyzz,   \
                             tr_y_yyz_xxyyyy,   \
                             tr_y_yyz_xxyyyz,   \
                             tr_y_yyz_xxyyzz,   \
                             tr_y_yyz_xxyzzz,   \
                             tr_y_yyz_xyyyyy,   \
                             tr_y_yyz_xyyyyz,   \
                             tr_y_yyz_xyyyzz,   \
                             tr_y_yyz_xyyzzz,   \
                             tr_y_yyz_xyzzzz,   \
                             tr_y_yyz_yyyyyy,   \
                             tr_y_yyz_yyyyyz,   \
                             tr_y_yyz_yyyyzz,   \
                             tr_y_yyz_yyyzzz,   \
                             tr_y_yyz_yyzzzz,   \
                             tr_y_yyz_yzzzzz,   \
                             tr_y_yyzz_xxxxxx,  \
                             tr_y_yyzz_xxxxxy,  \
                             tr_y_yyzz_xxxxy,   \
                             tr_y_yyzz_xxxxyy,  \
                             tr_y_yyzz_xxxxyz,  \
                             tr_y_yyzz_xxxyy,   \
                             tr_y_yyzz_xxxyyy,  \
                             tr_y_yyzz_xxxyyz,  \
                             tr_y_yyzz_xxxyz,   \
                             tr_y_yyzz_xxxyzz,  \
                             tr_y_yyzz_xxyyy,   \
                             tr_y_yyzz_xxyyyy,  \
                             tr_y_yyzz_xxyyyz,  \
                             tr_y_yyzz_xxyyz,   \
                             tr_y_yyzz_xxyyzz,  \
                             tr_y_yyzz_xxyzz,   \
                             tr_y_yyzz_xxyzzz,  \
                             tr_y_yyzz_xyyyy,   \
                             tr_y_yyzz_xyyyyy,  \
                             tr_y_yyzz_xyyyyz,  \
                             tr_y_yyzz_xyyyz,   \
                             tr_y_yyzz_xyyyzz,  \
                             tr_y_yyzz_xyyzz,   \
                             tr_y_yyzz_xyyzzz,  \
                             tr_y_yyzz_xyzzz,   \
                             tr_y_yyzz_xyzzzz,  \
                             tr_y_yyzz_yyyyy,   \
                             tr_y_yyzz_yyyyyy,  \
                             tr_y_yyzz_yyyyyz,  \
                             tr_y_yyzz_yyyyz,   \
                             tr_y_yyzz_yyyyzz,  \
                             tr_y_yyzz_yyyzz,   \
                             tr_y_yyzz_yyyzzz,  \
                             tr_y_yyzz_yyzzz,   \
                             tr_y_yyzz_yyzzzz,  \
                             tr_y_yyzz_yzzzz,   \
                             tr_y_yyzz_yzzzzz,  \
                             tr_y_yyzzz_xxxxxx, \
                             tr_y_yyzzz_xxxxxy, \
                             tr_y_yyzzz_xxxxxz, \
                             tr_y_yyzzz_xxxxyy, \
                             tr_y_yyzzz_xxxxyz, \
                             tr_y_yyzzz_xxxxzz, \
                             tr_y_yyzzz_xxxyyy, \
                             tr_y_yyzzz_xxxyyz, \
                             tr_y_yyzzz_xxxyzz, \
                             tr_y_yyzzz_xxxzzz, \
                             tr_y_yyzzz_xxyyyy, \
                             tr_y_yyzzz_xxyyyz, \
                             tr_y_yyzzz_xxyyzz, \
                             tr_y_yyzzz_xxyzzz, \
                             tr_y_yyzzz_xxzzzz, \
                             tr_y_yyzzz_xyyyyy, \
                             tr_y_yyzzz_xyyyyz, \
                             tr_y_yyzzz_xyyyzz, \
                             tr_y_yyzzz_xyyzzz, \
                             tr_y_yyzzz_xyzzzz, \
                             tr_y_yyzzz_xzzzzz, \
                             tr_y_yyzzz_yyyyyy, \
                             tr_y_yyzzz_yyyyyz, \
                             tr_y_yyzzz_yyyyzz, \
                             tr_y_yyzzz_yyyzzz, \
                             tr_y_yyzzz_yyzzzz, \
                             tr_y_yyzzz_yzzzzz, \
                             tr_y_yyzzz_zzzzzz, \
                             tr_y_yzzz_xxxxxz,  \
                             tr_y_yzzz_xxxxzz,  \
                             tr_y_yzzz_xxxzzz,  \
                             tr_y_yzzz_xxzzzz,  \
                             tr_y_yzzz_xzzzzz,  \
                             tr_y_yzzz_zzzzzz,  \
                             tr_y_zzz_xxxxxz,   \
                             tr_y_zzz_xxxxzz,   \
                             tr_y_zzz_xxxzzz,   \
                             tr_y_zzz_xxzzzz,   \
                             tr_y_zzz_xzzzzz,   \
                             tr_y_zzz_zzzzzz,   \
                             ts_yzzz_xxxxxz,    \
                             ts_yzzz_xxxxzz,    \
                             ts_yzzz_xxxzzz,    \
                             ts_yzzz_xxzzzz,    \
                             ts_yzzz_xzzzzz,    \
                             ts_yzzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzz_xxxxxx[i] = 2.0 * tr_y_yyz_xxxxxx[i] * fe_0 + tr_y_yyzz_xxxxxx[i] * pa_z[i];

        tr_y_yyzzz_xxxxxy[i] = 2.0 * tr_y_yyz_xxxxxy[i] * fe_0 + tr_y_yyzz_xxxxxy[i] * pa_z[i];

        tr_y_yyzzz_xxxxxz[i] = tr_y_zzz_xxxxxz[i] * fe_0 + ts_yzzz_xxxxxz[i] * fe_0 + tr_y_yzzz_xxxxxz[i] * pa_y[i];

        tr_y_yyzzz_xxxxyy[i] = 2.0 * tr_y_yyz_xxxxyy[i] * fe_0 + tr_y_yyzz_xxxxyy[i] * pa_z[i];

        tr_y_yyzzz_xxxxyz[i] = 2.0 * tr_y_yyz_xxxxyz[i] * fe_0 + tr_y_yyzz_xxxxy[i] * fe_0 + tr_y_yyzz_xxxxyz[i] * pa_z[i];

        tr_y_yyzzz_xxxxzz[i] = tr_y_zzz_xxxxzz[i] * fe_0 + ts_yzzz_xxxxzz[i] * fe_0 + tr_y_yzzz_xxxxzz[i] * pa_y[i];

        tr_y_yyzzz_xxxyyy[i] = 2.0 * tr_y_yyz_xxxyyy[i] * fe_0 + tr_y_yyzz_xxxyyy[i] * pa_z[i];

        tr_y_yyzzz_xxxyyz[i] = 2.0 * tr_y_yyz_xxxyyz[i] * fe_0 + tr_y_yyzz_xxxyy[i] * fe_0 + tr_y_yyzz_xxxyyz[i] * pa_z[i];

        tr_y_yyzzz_xxxyzz[i] = 2.0 * tr_y_yyz_xxxyzz[i] * fe_0 + 2.0 * tr_y_yyzz_xxxyz[i] * fe_0 + tr_y_yyzz_xxxyzz[i] * pa_z[i];

        tr_y_yyzzz_xxxzzz[i] = tr_y_zzz_xxxzzz[i] * fe_0 + ts_yzzz_xxxzzz[i] * fe_0 + tr_y_yzzz_xxxzzz[i] * pa_y[i];

        tr_y_yyzzz_xxyyyy[i] = 2.0 * tr_y_yyz_xxyyyy[i] * fe_0 + tr_y_yyzz_xxyyyy[i] * pa_z[i];

        tr_y_yyzzz_xxyyyz[i] = 2.0 * tr_y_yyz_xxyyyz[i] * fe_0 + tr_y_yyzz_xxyyy[i] * fe_0 + tr_y_yyzz_xxyyyz[i] * pa_z[i];

        tr_y_yyzzz_xxyyzz[i] = 2.0 * tr_y_yyz_xxyyzz[i] * fe_0 + 2.0 * tr_y_yyzz_xxyyz[i] * fe_0 + tr_y_yyzz_xxyyzz[i] * pa_z[i];

        tr_y_yyzzz_xxyzzz[i] = 2.0 * tr_y_yyz_xxyzzz[i] * fe_0 + 3.0 * tr_y_yyzz_xxyzz[i] * fe_0 + tr_y_yyzz_xxyzzz[i] * pa_z[i];

        tr_y_yyzzz_xxzzzz[i] = tr_y_zzz_xxzzzz[i] * fe_0 + ts_yzzz_xxzzzz[i] * fe_0 + tr_y_yzzz_xxzzzz[i] * pa_y[i];

        tr_y_yyzzz_xyyyyy[i] = 2.0 * tr_y_yyz_xyyyyy[i] * fe_0 + tr_y_yyzz_xyyyyy[i] * pa_z[i];

        tr_y_yyzzz_xyyyyz[i] = 2.0 * tr_y_yyz_xyyyyz[i] * fe_0 + tr_y_yyzz_xyyyy[i] * fe_0 + tr_y_yyzz_xyyyyz[i] * pa_z[i];

        tr_y_yyzzz_xyyyzz[i] = 2.0 * tr_y_yyz_xyyyzz[i] * fe_0 + 2.0 * tr_y_yyzz_xyyyz[i] * fe_0 + tr_y_yyzz_xyyyzz[i] * pa_z[i];

        tr_y_yyzzz_xyyzzz[i] = 2.0 * tr_y_yyz_xyyzzz[i] * fe_0 + 3.0 * tr_y_yyzz_xyyzz[i] * fe_0 + tr_y_yyzz_xyyzzz[i] * pa_z[i];

        tr_y_yyzzz_xyzzzz[i] = 2.0 * tr_y_yyz_xyzzzz[i] * fe_0 + 4.0 * tr_y_yyzz_xyzzz[i] * fe_0 + tr_y_yyzz_xyzzzz[i] * pa_z[i];

        tr_y_yyzzz_xzzzzz[i] = tr_y_zzz_xzzzzz[i] * fe_0 + ts_yzzz_xzzzzz[i] * fe_0 + tr_y_yzzz_xzzzzz[i] * pa_y[i];

        tr_y_yyzzz_yyyyyy[i] = 2.0 * tr_y_yyz_yyyyyy[i] * fe_0 + tr_y_yyzz_yyyyyy[i] * pa_z[i];

        tr_y_yyzzz_yyyyyz[i] = 2.0 * tr_y_yyz_yyyyyz[i] * fe_0 + tr_y_yyzz_yyyyy[i] * fe_0 + tr_y_yyzz_yyyyyz[i] * pa_z[i];

        tr_y_yyzzz_yyyyzz[i] = 2.0 * tr_y_yyz_yyyyzz[i] * fe_0 + 2.0 * tr_y_yyzz_yyyyz[i] * fe_0 + tr_y_yyzz_yyyyzz[i] * pa_z[i];

        tr_y_yyzzz_yyyzzz[i] = 2.0 * tr_y_yyz_yyyzzz[i] * fe_0 + 3.0 * tr_y_yyzz_yyyzz[i] * fe_0 + tr_y_yyzz_yyyzzz[i] * pa_z[i];

        tr_y_yyzzz_yyzzzz[i] = 2.0 * tr_y_yyz_yyzzzz[i] * fe_0 + 4.0 * tr_y_yyzz_yyzzz[i] * fe_0 + tr_y_yyzz_yyzzzz[i] * pa_z[i];

        tr_y_yyzzz_yzzzzz[i] = 2.0 * tr_y_yyz_yzzzzz[i] * fe_0 + 5.0 * tr_y_yyzz_yzzzz[i] * fe_0 + tr_y_yyzz_yzzzzz[i] * pa_z[i];

        tr_y_yyzzz_zzzzzz[i] = tr_y_zzz_zzzzzz[i] * fe_0 + ts_yzzz_zzzzzz[i] * fe_0 + tr_y_yzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1120-1148 components of targeted buffer : HI

    auto tr_y_yzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1120);

    auto tr_y_yzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1121);

    auto tr_y_yzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1122);

    auto tr_y_yzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1123);

    auto tr_y_yzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1124);

    auto tr_y_yzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1125);

    auto tr_y_yzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1126);

    auto tr_y_yzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1127);

    auto tr_y_yzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1128);

    auto tr_y_yzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1129);

    auto tr_y_yzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1130);

    auto tr_y_yzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1131);

    auto tr_y_yzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1132);

    auto tr_y_yzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1133);

    auto tr_y_yzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1134);

    auto tr_y_yzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1135);

    auto tr_y_yzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1136);

    auto tr_y_yzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1137);

    auto tr_y_yzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1138);

    auto tr_y_yzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1139);

    auto tr_y_yzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1140);

    auto tr_y_yzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1141);

    auto tr_y_yzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1142);

    auto tr_y_yzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1143);

    auto tr_y_yzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1144);

    auto tr_y_yzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1145);

    auto tr_y_yzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1146);

    auto tr_y_yzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1147);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_y_yzz_xxxxxy,   \
                             tr_y_yzz_xxxxyy,   \
                             tr_y_yzz_xxxyyy,   \
                             tr_y_yzz_xxyyyy,   \
                             tr_y_yzz_xyyyyy,   \
                             tr_y_yzz_yyyyyy,   \
                             tr_y_yzzz_xxxxxy,  \
                             tr_y_yzzz_xxxxyy,  \
                             tr_y_yzzz_xxxyyy,  \
                             tr_y_yzzz_xxyyyy,  \
                             tr_y_yzzz_xyyyyy,  \
                             tr_y_yzzz_yyyyyy,  \
                             tr_y_yzzzz_xxxxxx, \
                             tr_y_yzzzz_xxxxxy, \
                             tr_y_yzzzz_xxxxxz, \
                             tr_y_yzzzz_xxxxyy, \
                             tr_y_yzzzz_xxxxyz, \
                             tr_y_yzzzz_xxxxzz, \
                             tr_y_yzzzz_xxxyyy, \
                             tr_y_yzzzz_xxxyyz, \
                             tr_y_yzzzz_xxxyzz, \
                             tr_y_yzzzz_xxxzzz, \
                             tr_y_yzzzz_xxyyyy, \
                             tr_y_yzzzz_xxyyyz, \
                             tr_y_yzzzz_xxyyzz, \
                             tr_y_yzzzz_xxyzzz, \
                             tr_y_yzzzz_xxzzzz, \
                             tr_y_yzzzz_xyyyyy, \
                             tr_y_yzzzz_xyyyyz, \
                             tr_y_yzzzz_xyyyzz, \
                             tr_y_yzzzz_xyyzzz, \
                             tr_y_yzzzz_xyzzzz, \
                             tr_y_yzzzz_xzzzzz, \
                             tr_y_yzzzz_yyyyyy, \
                             tr_y_yzzzz_yyyyyz, \
                             tr_y_yzzzz_yyyyzz, \
                             tr_y_yzzzz_yyyzzz, \
                             tr_y_yzzzz_yyzzzz, \
                             tr_y_yzzzz_yzzzzz, \
                             tr_y_yzzzz_zzzzzz, \
                             tr_y_zzzz_xxxxxx,  \
                             tr_y_zzzz_xxxxxz,  \
                             tr_y_zzzz_xxxxyz,  \
                             tr_y_zzzz_xxxxz,   \
                             tr_y_zzzz_xxxxzz,  \
                             tr_y_zzzz_xxxyyz,  \
                             tr_y_zzzz_xxxyz,   \
                             tr_y_zzzz_xxxyzz,  \
                             tr_y_zzzz_xxxzz,   \
                             tr_y_zzzz_xxxzzz,  \
                             tr_y_zzzz_xxyyyz,  \
                             tr_y_zzzz_xxyyz,   \
                             tr_y_zzzz_xxyyzz,  \
                             tr_y_zzzz_xxyzz,   \
                             tr_y_zzzz_xxyzzz,  \
                             tr_y_zzzz_xxzzz,   \
                             tr_y_zzzz_xxzzzz,  \
                             tr_y_zzzz_xyyyyz,  \
                             tr_y_zzzz_xyyyz,   \
                             tr_y_zzzz_xyyyzz,  \
                             tr_y_zzzz_xyyzz,   \
                             tr_y_zzzz_xyyzzz,  \
                             tr_y_zzzz_xyzzz,   \
                             tr_y_zzzz_xyzzzz,  \
                             tr_y_zzzz_xzzzz,   \
                             tr_y_zzzz_xzzzzz,  \
                             tr_y_zzzz_yyyyyz,  \
                             tr_y_zzzz_yyyyz,   \
                             tr_y_zzzz_yyyyzz,  \
                             tr_y_zzzz_yyyzz,   \
                             tr_y_zzzz_yyyzzz,  \
                             tr_y_zzzz_yyzzz,   \
                             tr_y_zzzz_yyzzzz,  \
                             tr_y_zzzz_yzzzz,   \
                             tr_y_zzzz_yzzzzz,  \
                             tr_y_zzzz_zzzzz,   \
                             tr_y_zzzz_zzzzzz,  \
                             ts_zzzz_xxxxxx,    \
                             ts_zzzz_xxxxxz,    \
                             ts_zzzz_xxxxyz,    \
                             ts_zzzz_xxxxzz,    \
                             ts_zzzz_xxxyyz,    \
                             ts_zzzz_xxxyzz,    \
                             ts_zzzz_xxxzzz,    \
                             ts_zzzz_xxyyyz,    \
                             ts_zzzz_xxyyzz,    \
                             ts_zzzz_xxyzzz,    \
                             ts_zzzz_xxzzzz,    \
                             ts_zzzz_xyyyyz,    \
                             ts_zzzz_xyyyzz,    \
                             ts_zzzz_xyyzzz,    \
                             ts_zzzz_xyzzzz,    \
                             ts_zzzz_xzzzzz,    \
                             ts_zzzz_yyyyyz,    \
                             ts_zzzz_yyyyzz,    \
                             ts_zzzz_yyyzzz,    \
                             ts_zzzz_yyzzzz,    \
                             ts_zzzz_yzzzzz,    \
                             ts_zzzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzz_xxxxxx[i] = ts_zzzz_xxxxxx[i] * fe_0 + tr_y_zzzz_xxxxxx[i] * pa_y[i];

        tr_y_yzzzz_xxxxxy[i] = 3.0 * tr_y_yzz_xxxxxy[i] * fe_0 + tr_y_yzzz_xxxxxy[i] * pa_z[i];

        tr_y_yzzzz_xxxxxz[i] = ts_zzzz_xxxxxz[i] * fe_0 + tr_y_zzzz_xxxxxz[i] * pa_y[i];

        tr_y_yzzzz_xxxxyy[i] = 3.0 * tr_y_yzz_xxxxyy[i] * fe_0 + tr_y_yzzz_xxxxyy[i] * pa_z[i];

        tr_y_yzzzz_xxxxyz[i] = tr_y_zzzz_xxxxz[i] * fe_0 + ts_zzzz_xxxxyz[i] * fe_0 + tr_y_zzzz_xxxxyz[i] * pa_y[i];

        tr_y_yzzzz_xxxxzz[i] = ts_zzzz_xxxxzz[i] * fe_0 + tr_y_zzzz_xxxxzz[i] * pa_y[i];

        tr_y_yzzzz_xxxyyy[i] = 3.0 * tr_y_yzz_xxxyyy[i] * fe_0 + tr_y_yzzz_xxxyyy[i] * pa_z[i];

        tr_y_yzzzz_xxxyyz[i] = 2.0 * tr_y_zzzz_xxxyz[i] * fe_0 + ts_zzzz_xxxyyz[i] * fe_0 + tr_y_zzzz_xxxyyz[i] * pa_y[i];

        tr_y_yzzzz_xxxyzz[i] = tr_y_zzzz_xxxzz[i] * fe_0 + ts_zzzz_xxxyzz[i] * fe_0 + tr_y_zzzz_xxxyzz[i] * pa_y[i];

        tr_y_yzzzz_xxxzzz[i] = ts_zzzz_xxxzzz[i] * fe_0 + tr_y_zzzz_xxxzzz[i] * pa_y[i];

        tr_y_yzzzz_xxyyyy[i] = 3.0 * tr_y_yzz_xxyyyy[i] * fe_0 + tr_y_yzzz_xxyyyy[i] * pa_z[i];

        tr_y_yzzzz_xxyyyz[i] = 3.0 * tr_y_zzzz_xxyyz[i] * fe_0 + ts_zzzz_xxyyyz[i] * fe_0 + tr_y_zzzz_xxyyyz[i] * pa_y[i];

        tr_y_yzzzz_xxyyzz[i] = 2.0 * tr_y_zzzz_xxyzz[i] * fe_0 + ts_zzzz_xxyyzz[i] * fe_0 + tr_y_zzzz_xxyyzz[i] * pa_y[i];

        tr_y_yzzzz_xxyzzz[i] = tr_y_zzzz_xxzzz[i] * fe_0 + ts_zzzz_xxyzzz[i] * fe_0 + tr_y_zzzz_xxyzzz[i] * pa_y[i];

        tr_y_yzzzz_xxzzzz[i] = ts_zzzz_xxzzzz[i] * fe_0 + tr_y_zzzz_xxzzzz[i] * pa_y[i];

        tr_y_yzzzz_xyyyyy[i] = 3.0 * tr_y_yzz_xyyyyy[i] * fe_0 + tr_y_yzzz_xyyyyy[i] * pa_z[i];

        tr_y_yzzzz_xyyyyz[i] = 4.0 * tr_y_zzzz_xyyyz[i] * fe_0 + ts_zzzz_xyyyyz[i] * fe_0 + tr_y_zzzz_xyyyyz[i] * pa_y[i];

        tr_y_yzzzz_xyyyzz[i] = 3.0 * tr_y_zzzz_xyyzz[i] * fe_0 + ts_zzzz_xyyyzz[i] * fe_0 + tr_y_zzzz_xyyyzz[i] * pa_y[i];

        tr_y_yzzzz_xyyzzz[i] = 2.0 * tr_y_zzzz_xyzzz[i] * fe_0 + ts_zzzz_xyyzzz[i] * fe_0 + tr_y_zzzz_xyyzzz[i] * pa_y[i];

        tr_y_yzzzz_xyzzzz[i] = tr_y_zzzz_xzzzz[i] * fe_0 + ts_zzzz_xyzzzz[i] * fe_0 + tr_y_zzzz_xyzzzz[i] * pa_y[i];

        tr_y_yzzzz_xzzzzz[i] = ts_zzzz_xzzzzz[i] * fe_0 + tr_y_zzzz_xzzzzz[i] * pa_y[i];

        tr_y_yzzzz_yyyyyy[i] = 3.0 * tr_y_yzz_yyyyyy[i] * fe_0 + tr_y_yzzz_yyyyyy[i] * pa_z[i];

        tr_y_yzzzz_yyyyyz[i] = 5.0 * tr_y_zzzz_yyyyz[i] * fe_0 + ts_zzzz_yyyyyz[i] * fe_0 + tr_y_zzzz_yyyyyz[i] * pa_y[i];

        tr_y_yzzzz_yyyyzz[i] = 4.0 * tr_y_zzzz_yyyzz[i] * fe_0 + ts_zzzz_yyyyzz[i] * fe_0 + tr_y_zzzz_yyyyzz[i] * pa_y[i];

        tr_y_yzzzz_yyyzzz[i] = 3.0 * tr_y_zzzz_yyzzz[i] * fe_0 + ts_zzzz_yyyzzz[i] * fe_0 + tr_y_zzzz_yyyzzz[i] * pa_y[i];

        tr_y_yzzzz_yyzzzz[i] = 2.0 * tr_y_zzzz_yzzzz[i] * fe_0 + ts_zzzz_yyzzzz[i] * fe_0 + tr_y_zzzz_yyzzzz[i] * pa_y[i];

        tr_y_yzzzz_yzzzzz[i] = tr_y_zzzz_zzzzz[i] * fe_0 + ts_zzzz_yzzzzz[i] * fe_0 + tr_y_zzzz_yzzzzz[i] * pa_y[i];

        tr_y_yzzzz_zzzzzz[i] = ts_zzzz_zzzzzz[i] * fe_0 + tr_y_zzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1148-1176 components of targeted buffer : HI

    auto tr_y_zzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1148);

    auto tr_y_zzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1149);

    auto tr_y_zzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1150);

    auto tr_y_zzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1151);

    auto tr_y_zzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1152);

    auto tr_y_zzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1153);

    auto tr_y_zzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1154);

    auto tr_y_zzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1155);

    auto tr_y_zzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1156);

    auto tr_y_zzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1157);

    auto tr_y_zzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1158);

    auto tr_y_zzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1159);

    auto tr_y_zzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1160);

    auto tr_y_zzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1161);

    auto tr_y_zzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1162);

    auto tr_y_zzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1163);

    auto tr_y_zzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1164);

    auto tr_y_zzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1165);

    auto tr_y_zzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1166);

    auto tr_y_zzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1167);

    auto tr_y_zzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1168);

    auto tr_y_zzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1169);

    auto tr_y_zzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1170);

    auto tr_y_zzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1171);

    auto tr_y_zzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1172);

    auto tr_y_zzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1173);

    auto tr_y_zzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1174);

    auto tr_y_zzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1175);

#pragma omp simd aligned(pa_z,                  \
                             tr_y_zzz_xxxxxx,   \
                             tr_y_zzz_xxxxxy,   \
                             tr_y_zzz_xxxxxz,   \
                             tr_y_zzz_xxxxyy,   \
                             tr_y_zzz_xxxxyz,   \
                             tr_y_zzz_xxxxzz,   \
                             tr_y_zzz_xxxyyy,   \
                             tr_y_zzz_xxxyyz,   \
                             tr_y_zzz_xxxyzz,   \
                             tr_y_zzz_xxxzzz,   \
                             tr_y_zzz_xxyyyy,   \
                             tr_y_zzz_xxyyyz,   \
                             tr_y_zzz_xxyyzz,   \
                             tr_y_zzz_xxyzzz,   \
                             tr_y_zzz_xxzzzz,   \
                             tr_y_zzz_xyyyyy,   \
                             tr_y_zzz_xyyyyz,   \
                             tr_y_zzz_xyyyzz,   \
                             tr_y_zzz_xyyzzz,   \
                             tr_y_zzz_xyzzzz,   \
                             tr_y_zzz_xzzzzz,   \
                             tr_y_zzz_yyyyyy,   \
                             tr_y_zzz_yyyyyz,   \
                             tr_y_zzz_yyyyzz,   \
                             tr_y_zzz_yyyzzz,   \
                             tr_y_zzz_yyzzzz,   \
                             tr_y_zzz_yzzzzz,   \
                             tr_y_zzz_zzzzzz,   \
                             tr_y_zzzz_xxxxx,   \
                             tr_y_zzzz_xxxxxx,  \
                             tr_y_zzzz_xxxxxy,  \
                             tr_y_zzzz_xxxxxz,  \
                             tr_y_zzzz_xxxxy,   \
                             tr_y_zzzz_xxxxyy,  \
                             tr_y_zzzz_xxxxyz,  \
                             tr_y_zzzz_xxxxz,   \
                             tr_y_zzzz_xxxxzz,  \
                             tr_y_zzzz_xxxyy,   \
                             tr_y_zzzz_xxxyyy,  \
                             tr_y_zzzz_xxxyyz,  \
                             tr_y_zzzz_xxxyz,   \
                             tr_y_zzzz_xxxyzz,  \
                             tr_y_zzzz_xxxzz,   \
                             tr_y_zzzz_xxxzzz,  \
                             tr_y_zzzz_xxyyy,   \
                             tr_y_zzzz_xxyyyy,  \
                             tr_y_zzzz_xxyyyz,  \
                             tr_y_zzzz_xxyyz,   \
                             tr_y_zzzz_xxyyzz,  \
                             tr_y_zzzz_xxyzz,   \
                             tr_y_zzzz_xxyzzz,  \
                             tr_y_zzzz_xxzzz,   \
                             tr_y_zzzz_xxzzzz,  \
                             tr_y_zzzz_xyyyy,   \
                             tr_y_zzzz_xyyyyy,  \
                             tr_y_zzzz_xyyyyz,  \
                             tr_y_zzzz_xyyyz,   \
                             tr_y_zzzz_xyyyzz,  \
                             tr_y_zzzz_xyyzz,   \
                             tr_y_zzzz_xyyzzz,  \
                             tr_y_zzzz_xyzzz,   \
                             tr_y_zzzz_xyzzzz,  \
                             tr_y_zzzz_xzzzz,   \
                             tr_y_zzzz_xzzzzz,  \
                             tr_y_zzzz_yyyyy,   \
                             tr_y_zzzz_yyyyyy,  \
                             tr_y_zzzz_yyyyyz,  \
                             tr_y_zzzz_yyyyz,   \
                             tr_y_zzzz_yyyyzz,  \
                             tr_y_zzzz_yyyzz,   \
                             tr_y_zzzz_yyyzzz,  \
                             tr_y_zzzz_yyzzz,   \
                             tr_y_zzzz_yyzzzz,  \
                             tr_y_zzzz_yzzzz,   \
                             tr_y_zzzz_yzzzzz,  \
                             tr_y_zzzz_zzzzz,   \
                             tr_y_zzzz_zzzzzz,  \
                             tr_y_zzzzz_xxxxxx, \
                             tr_y_zzzzz_xxxxxy, \
                             tr_y_zzzzz_xxxxxz, \
                             tr_y_zzzzz_xxxxyy, \
                             tr_y_zzzzz_xxxxyz, \
                             tr_y_zzzzz_xxxxzz, \
                             tr_y_zzzzz_xxxyyy, \
                             tr_y_zzzzz_xxxyyz, \
                             tr_y_zzzzz_xxxyzz, \
                             tr_y_zzzzz_xxxzzz, \
                             tr_y_zzzzz_xxyyyy, \
                             tr_y_zzzzz_xxyyyz, \
                             tr_y_zzzzz_xxyyzz, \
                             tr_y_zzzzz_xxyzzz, \
                             tr_y_zzzzz_xxzzzz, \
                             tr_y_zzzzz_xyyyyy, \
                             tr_y_zzzzz_xyyyyz, \
                             tr_y_zzzzz_xyyyzz, \
                             tr_y_zzzzz_xyyzzz, \
                             tr_y_zzzzz_xyzzzz, \
                             tr_y_zzzzz_xzzzzz, \
                             tr_y_zzzzz_yyyyyy, \
                             tr_y_zzzzz_yyyyyz, \
                             tr_y_zzzzz_yyyyzz, \
                             tr_y_zzzzz_yyyzzz, \
                             tr_y_zzzzz_yyzzzz, \
                             tr_y_zzzzz_yzzzzz, \
                             tr_y_zzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzz_xxxxxx[i] = 4.0 * tr_y_zzz_xxxxxx[i] * fe_0 + tr_y_zzzz_xxxxxx[i] * pa_z[i];

        tr_y_zzzzz_xxxxxy[i] = 4.0 * tr_y_zzz_xxxxxy[i] * fe_0 + tr_y_zzzz_xxxxxy[i] * pa_z[i];

        tr_y_zzzzz_xxxxxz[i] = 4.0 * tr_y_zzz_xxxxxz[i] * fe_0 + tr_y_zzzz_xxxxx[i] * fe_0 + tr_y_zzzz_xxxxxz[i] * pa_z[i];

        tr_y_zzzzz_xxxxyy[i] = 4.0 * tr_y_zzz_xxxxyy[i] * fe_0 + tr_y_zzzz_xxxxyy[i] * pa_z[i];

        tr_y_zzzzz_xxxxyz[i] = 4.0 * tr_y_zzz_xxxxyz[i] * fe_0 + tr_y_zzzz_xxxxy[i] * fe_0 + tr_y_zzzz_xxxxyz[i] * pa_z[i];

        tr_y_zzzzz_xxxxzz[i] = 4.0 * tr_y_zzz_xxxxzz[i] * fe_0 + 2.0 * tr_y_zzzz_xxxxz[i] * fe_0 + tr_y_zzzz_xxxxzz[i] * pa_z[i];

        tr_y_zzzzz_xxxyyy[i] = 4.0 * tr_y_zzz_xxxyyy[i] * fe_0 + tr_y_zzzz_xxxyyy[i] * pa_z[i];

        tr_y_zzzzz_xxxyyz[i] = 4.0 * tr_y_zzz_xxxyyz[i] * fe_0 + tr_y_zzzz_xxxyy[i] * fe_0 + tr_y_zzzz_xxxyyz[i] * pa_z[i];

        tr_y_zzzzz_xxxyzz[i] = 4.0 * tr_y_zzz_xxxyzz[i] * fe_0 + 2.0 * tr_y_zzzz_xxxyz[i] * fe_0 + tr_y_zzzz_xxxyzz[i] * pa_z[i];

        tr_y_zzzzz_xxxzzz[i] = 4.0 * tr_y_zzz_xxxzzz[i] * fe_0 + 3.0 * tr_y_zzzz_xxxzz[i] * fe_0 + tr_y_zzzz_xxxzzz[i] * pa_z[i];

        tr_y_zzzzz_xxyyyy[i] = 4.0 * tr_y_zzz_xxyyyy[i] * fe_0 + tr_y_zzzz_xxyyyy[i] * pa_z[i];

        tr_y_zzzzz_xxyyyz[i] = 4.0 * tr_y_zzz_xxyyyz[i] * fe_0 + tr_y_zzzz_xxyyy[i] * fe_0 + tr_y_zzzz_xxyyyz[i] * pa_z[i];

        tr_y_zzzzz_xxyyzz[i] = 4.0 * tr_y_zzz_xxyyzz[i] * fe_0 + 2.0 * tr_y_zzzz_xxyyz[i] * fe_0 + tr_y_zzzz_xxyyzz[i] * pa_z[i];

        tr_y_zzzzz_xxyzzz[i] = 4.0 * tr_y_zzz_xxyzzz[i] * fe_0 + 3.0 * tr_y_zzzz_xxyzz[i] * fe_0 + tr_y_zzzz_xxyzzz[i] * pa_z[i];

        tr_y_zzzzz_xxzzzz[i] = 4.0 * tr_y_zzz_xxzzzz[i] * fe_0 + 4.0 * tr_y_zzzz_xxzzz[i] * fe_0 + tr_y_zzzz_xxzzzz[i] * pa_z[i];

        tr_y_zzzzz_xyyyyy[i] = 4.0 * tr_y_zzz_xyyyyy[i] * fe_0 + tr_y_zzzz_xyyyyy[i] * pa_z[i];

        tr_y_zzzzz_xyyyyz[i] = 4.0 * tr_y_zzz_xyyyyz[i] * fe_0 + tr_y_zzzz_xyyyy[i] * fe_0 + tr_y_zzzz_xyyyyz[i] * pa_z[i];

        tr_y_zzzzz_xyyyzz[i] = 4.0 * tr_y_zzz_xyyyzz[i] * fe_0 + 2.0 * tr_y_zzzz_xyyyz[i] * fe_0 + tr_y_zzzz_xyyyzz[i] * pa_z[i];

        tr_y_zzzzz_xyyzzz[i] = 4.0 * tr_y_zzz_xyyzzz[i] * fe_0 + 3.0 * tr_y_zzzz_xyyzz[i] * fe_0 + tr_y_zzzz_xyyzzz[i] * pa_z[i];

        tr_y_zzzzz_xyzzzz[i] = 4.0 * tr_y_zzz_xyzzzz[i] * fe_0 + 4.0 * tr_y_zzzz_xyzzz[i] * fe_0 + tr_y_zzzz_xyzzzz[i] * pa_z[i];

        tr_y_zzzzz_xzzzzz[i] = 4.0 * tr_y_zzz_xzzzzz[i] * fe_0 + 5.0 * tr_y_zzzz_xzzzz[i] * fe_0 + tr_y_zzzz_xzzzzz[i] * pa_z[i];

        tr_y_zzzzz_yyyyyy[i] = 4.0 * tr_y_zzz_yyyyyy[i] * fe_0 + tr_y_zzzz_yyyyyy[i] * pa_z[i];

        tr_y_zzzzz_yyyyyz[i] = 4.0 * tr_y_zzz_yyyyyz[i] * fe_0 + tr_y_zzzz_yyyyy[i] * fe_0 + tr_y_zzzz_yyyyyz[i] * pa_z[i];

        tr_y_zzzzz_yyyyzz[i] = 4.0 * tr_y_zzz_yyyyzz[i] * fe_0 + 2.0 * tr_y_zzzz_yyyyz[i] * fe_0 + tr_y_zzzz_yyyyzz[i] * pa_z[i];

        tr_y_zzzzz_yyyzzz[i] = 4.0 * tr_y_zzz_yyyzzz[i] * fe_0 + 3.0 * tr_y_zzzz_yyyzz[i] * fe_0 + tr_y_zzzz_yyyzzz[i] * pa_z[i];

        tr_y_zzzzz_yyzzzz[i] = 4.0 * tr_y_zzz_yyzzzz[i] * fe_0 + 4.0 * tr_y_zzzz_yyzzz[i] * fe_0 + tr_y_zzzz_yyzzzz[i] * pa_z[i];

        tr_y_zzzzz_yzzzzz[i] = 4.0 * tr_y_zzz_yzzzzz[i] * fe_0 + 5.0 * tr_y_zzzz_yzzzz[i] * fe_0 + tr_y_zzzz_yzzzzz[i] * pa_z[i];

        tr_y_zzzzz_zzzzzz[i] = 4.0 * tr_y_zzz_zzzzzz[i] * fe_0 + 6.0 * tr_y_zzzz_zzzzz[i] * fe_0 + tr_y_zzzz_zzzzzz[i] * pa_z[i];
    }

    // Set up 1176-1204 components of targeted buffer : HI

    auto tr_z_xxxxx_xxxxxx = pbuffer.data(idx_dip_hi + 1176);

    auto tr_z_xxxxx_xxxxxy = pbuffer.data(idx_dip_hi + 1177);

    auto tr_z_xxxxx_xxxxxz = pbuffer.data(idx_dip_hi + 1178);

    auto tr_z_xxxxx_xxxxyy = pbuffer.data(idx_dip_hi + 1179);

    auto tr_z_xxxxx_xxxxyz = pbuffer.data(idx_dip_hi + 1180);

    auto tr_z_xxxxx_xxxxzz = pbuffer.data(idx_dip_hi + 1181);

    auto tr_z_xxxxx_xxxyyy = pbuffer.data(idx_dip_hi + 1182);

    auto tr_z_xxxxx_xxxyyz = pbuffer.data(idx_dip_hi + 1183);

    auto tr_z_xxxxx_xxxyzz = pbuffer.data(idx_dip_hi + 1184);

    auto tr_z_xxxxx_xxxzzz = pbuffer.data(idx_dip_hi + 1185);

    auto tr_z_xxxxx_xxyyyy = pbuffer.data(idx_dip_hi + 1186);

    auto tr_z_xxxxx_xxyyyz = pbuffer.data(idx_dip_hi + 1187);

    auto tr_z_xxxxx_xxyyzz = pbuffer.data(idx_dip_hi + 1188);

    auto tr_z_xxxxx_xxyzzz = pbuffer.data(idx_dip_hi + 1189);

    auto tr_z_xxxxx_xxzzzz = pbuffer.data(idx_dip_hi + 1190);

    auto tr_z_xxxxx_xyyyyy = pbuffer.data(idx_dip_hi + 1191);

    auto tr_z_xxxxx_xyyyyz = pbuffer.data(idx_dip_hi + 1192);

    auto tr_z_xxxxx_xyyyzz = pbuffer.data(idx_dip_hi + 1193);

    auto tr_z_xxxxx_xyyzzz = pbuffer.data(idx_dip_hi + 1194);

    auto tr_z_xxxxx_xyzzzz = pbuffer.data(idx_dip_hi + 1195);

    auto tr_z_xxxxx_xzzzzz = pbuffer.data(idx_dip_hi + 1196);

    auto tr_z_xxxxx_yyyyyy = pbuffer.data(idx_dip_hi + 1197);

    auto tr_z_xxxxx_yyyyyz = pbuffer.data(idx_dip_hi + 1198);

    auto tr_z_xxxxx_yyyyzz = pbuffer.data(idx_dip_hi + 1199);

    auto tr_z_xxxxx_yyyzzz = pbuffer.data(idx_dip_hi + 1200);

    auto tr_z_xxxxx_yyzzzz = pbuffer.data(idx_dip_hi + 1201);

    auto tr_z_xxxxx_yzzzzz = pbuffer.data(idx_dip_hi + 1202);

    auto tr_z_xxxxx_zzzzzz = pbuffer.data(idx_dip_hi + 1203);

#pragma omp simd aligned(pa_x,                  \
                             tr_z_xxx_xxxxxx,   \
                             tr_z_xxx_xxxxxy,   \
                             tr_z_xxx_xxxxxz,   \
                             tr_z_xxx_xxxxyy,   \
                             tr_z_xxx_xxxxyz,   \
                             tr_z_xxx_xxxxzz,   \
                             tr_z_xxx_xxxyyy,   \
                             tr_z_xxx_xxxyyz,   \
                             tr_z_xxx_xxxyzz,   \
                             tr_z_xxx_xxxzzz,   \
                             tr_z_xxx_xxyyyy,   \
                             tr_z_xxx_xxyyyz,   \
                             tr_z_xxx_xxyyzz,   \
                             tr_z_xxx_xxyzzz,   \
                             tr_z_xxx_xxzzzz,   \
                             tr_z_xxx_xyyyyy,   \
                             tr_z_xxx_xyyyyz,   \
                             tr_z_xxx_xyyyzz,   \
                             tr_z_xxx_xyyzzz,   \
                             tr_z_xxx_xyzzzz,   \
                             tr_z_xxx_xzzzzz,   \
                             tr_z_xxx_yyyyyy,   \
                             tr_z_xxx_yyyyyz,   \
                             tr_z_xxx_yyyyzz,   \
                             tr_z_xxx_yyyzzz,   \
                             tr_z_xxx_yyzzzz,   \
                             tr_z_xxx_yzzzzz,   \
                             tr_z_xxx_zzzzzz,   \
                             tr_z_xxxx_xxxxx,   \
                             tr_z_xxxx_xxxxxx,  \
                             tr_z_xxxx_xxxxxy,  \
                             tr_z_xxxx_xxxxxz,  \
                             tr_z_xxxx_xxxxy,   \
                             tr_z_xxxx_xxxxyy,  \
                             tr_z_xxxx_xxxxyz,  \
                             tr_z_xxxx_xxxxz,   \
                             tr_z_xxxx_xxxxzz,  \
                             tr_z_xxxx_xxxyy,   \
                             tr_z_xxxx_xxxyyy,  \
                             tr_z_xxxx_xxxyyz,  \
                             tr_z_xxxx_xxxyz,   \
                             tr_z_xxxx_xxxyzz,  \
                             tr_z_xxxx_xxxzz,   \
                             tr_z_xxxx_xxxzzz,  \
                             tr_z_xxxx_xxyyy,   \
                             tr_z_xxxx_xxyyyy,  \
                             tr_z_xxxx_xxyyyz,  \
                             tr_z_xxxx_xxyyz,   \
                             tr_z_xxxx_xxyyzz,  \
                             tr_z_xxxx_xxyzz,   \
                             tr_z_xxxx_xxyzzz,  \
                             tr_z_xxxx_xxzzz,   \
                             tr_z_xxxx_xxzzzz,  \
                             tr_z_xxxx_xyyyy,   \
                             tr_z_xxxx_xyyyyy,  \
                             tr_z_xxxx_xyyyyz,  \
                             tr_z_xxxx_xyyyz,   \
                             tr_z_xxxx_xyyyzz,  \
                             tr_z_xxxx_xyyzz,   \
                             tr_z_xxxx_xyyzzz,  \
                             tr_z_xxxx_xyzzz,   \
                             tr_z_xxxx_xyzzzz,  \
                             tr_z_xxxx_xzzzz,   \
                             tr_z_xxxx_xzzzzz,  \
                             tr_z_xxxx_yyyyy,   \
                             tr_z_xxxx_yyyyyy,  \
                             tr_z_xxxx_yyyyyz,  \
                             tr_z_xxxx_yyyyz,   \
                             tr_z_xxxx_yyyyzz,  \
                             tr_z_xxxx_yyyzz,   \
                             tr_z_xxxx_yyyzzz,  \
                             tr_z_xxxx_yyzzz,   \
                             tr_z_xxxx_yyzzzz,  \
                             tr_z_xxxx_yzzzz,   \
                             tr_z_xxxx_yzzzzz,  \
                             tr_z_xxxx_zzzzz,   \
                             tr_z_xxxx_zzzzzz,  \
                             tr_z_xxxxx_xxxxxx, \
                             tr_z_xxxxx_xxxxxy, \
                             tr_z_xxxxx_xxxxxz, \
                             tr_z_xxxxx_xxxxyy, \
                             tr_z_xxxxx_xxxxyz, \
                             tr_z_xxxxx_xxxxzz, \
                             tr_z_xxxxx_xxxyyy, \
                             tr_z_xxxxx_xxxyyz, \
                             tr_z_xxxxx_xxxyzz, \
                             tr_z_xxxxx_xxxzzz, \
                             tr_z_xxxxx_xxyyyy, \
                             tr_z_xxxxx_xxyyyz, \
                             tr_z_xxxxx_xxyyzz, \
                             tr_z_xxxxx_xxyzzz, \
                             tr_z_xxxxx_xxzzzz, \
                             tr_z_xxxxx_xyyyyy, \
                             tr_z_xxxxx_xyyyyz, \
                             tr_z_xxxxx_xyyyzz, \
                             tr_z_xxxxx_xyyzzz, \
                             tr_z_xxxxx_xyzzzz, \
                             tr_z_xxxxx_xzzzzz, \
                             tr_z_xxxxx_yyyyyy, \
                             tr_z_xxxxx_yyyyyz, \
                             tr_z_xxxxx_yyyyzz, \
                             tr_z_xxxxx_yyyzzz, \
                             tr_z_xxxxx_yyzzzz, \
                             tr_z_xxxxx_yzzzzz, \
                             tr_z_xxxxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxx_xxxxxx[i] = 4.0 * tr_z_xxx_xxxxxx[i] * fe_0 + 6.0 * tr_z_xxxx_xxxxx[i] * fe_0 + tr_z_xxxx_xxxxxx[i] * pa_x[i];

        tr_z_xxxxx_xxxxxy[i] = 4.0 * tr_z_xxx_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxxx_xxxxy[i] * fe_0 + tr_z_xxxx_xxxxxy[i] * pa_x[i];

        tr_z_xxxxx_xxxxxz[i] = 4.0 * tr_z_xxx_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxxx_xxxxz[i] * fe_0 + tr_z_xxxx_xxxxxz[i] * pa_x[i];

        tr_z_xxxxx_xxxxyy[i] = 4.0 * tr_z_xxx_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxxx_xxxyy[i] * fe_0 + tr_z_xxxx_xxxxyy[i] * pa_x[i];

        tr_z_xxxxx_xxxxyz[i] = 4.0 * tr_z_xxx_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxxx_xxxyz[i] * fe_0 + tr_z_xxxx_xxxxyz[i] * pa_x[i];

        tr_z_xxxxx_xxxxzz[i] = 4.0 * tr_z_xxx_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxxx_xxxzz[i] * fe_0 + tr_z_xxxx_xxxxzz[i] * pa_x[i];

        tr_z_xxxxx_xxxyyy[i] = 4.0 * tr_z_xxx_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxxx_xxyyy[i] * fe_0 + tr_z_xxxx_xxxyyy[i] * pa_x[i];

        tr_z_xxxxx_xxxyyz[i] = 4.0 * tr_z_xxx_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxxx_xxyyz[i] * fe_0 + tr_z_xxxx_xxxyyz[i] * pa_x[i];

        tr_z_xxxxx_xxxyzz[i] = 4.0 * tr_z_xxx_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxxx_xxyzz[i] * fe_0 + tr_z_xxxx_xxxyzz[i] * pa_x[i];

        tr_z_xxxxx_xxxzzz[i] = 4.0 * tr_z_xxx_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxxx_xxzzz[i] * fe_0 + tr_z_xxxx_xxxzzz[i] * pa_x[i];

        tr_z_xxxxx_xxyyyy[i] = 4.0 * tr_z_xxx_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxxx_xyyyy[i] * fe_0 + tr_z_xxxx_xxyyyy[i] * pa_x[i];

        tr_z_xxxxx_xxyyyz[i] = 4.0 * tr_z_xxx_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxxx_xyyyz[i] * fe_0 + tr_z_xxxx_xxyyyz[i] * pa_x[i];

        tr_z_xxxxx_xxyyzz[i] = 4.0 * tr_z_xxx_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxxx_xyyzz[i] * fe_0 + tr_z_xxxx_xxyyzz[i] * pa_x[i];

        tr_z_xxxxx_xxyzzz[i] = 4.0 * tr_z_xxx_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxxx_xyzzz[i] * fe_0 + tr_z_xxxx_xxyzzz[i] * pa_x[i];

        tr_z_xxxxx_xxzzzz[i] = 4.0 * tr_z_xxx_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxxx_xzzzz[i] * fe_0 + tr_z_xxxx_xxzzzz[i] * pa_x[i];

        tr_z_xxxxx_xyyyyy[i] = 4.0 * tr_z_xxx_xyyyyy[i] * fe_0 + tr_z_xxxx_yyyyy[i] * fe_0 + tr_z_xxxx_xyyyyy[i] * pa_x[i];

        tr_z_xxxxx_xyyyyz[i] = 4.0 * tr_z_xxx_xyyyyz[i] * fe_0 + tr_z_xxxx_yyyyz[i] * fe_0 + tr_z_xxxx_xyyyyz[i] * pa_x[i];

        tr_z_xxxxx_xyyyzz[i] = 4.0 * tr_z_xxx_xyyyzz[i] * fe_0 + tr_z_xxxx_yyyzz[i] * fe_0 + tr_z_xxxx_xyyyzz[i] * pa_x[i];

        tr_z_xxxxx_xyyzzz[i] = 4.0 * tr_z_xxx_xyyzzz[i] * fe_0 + tr_z_xxxx_yyzzz[i] * fe_0 + tr_z_xxxx_xyyzzz[i] * pa_x[i];

        tr_z_xxxxx_xyzzzz[i] = 4.0 * tr_z_xxx_xyzzzz[i] * fe_0 + tr_z_xxxx_yzzzz[i] * fe_0 + tr_z_xxxx_xyzzzz[i] * pa_x[i];

        tr_z_xxxxx_xzzzzz[i] = 4.0 * tr_z_xxx_xzzzzz[i] * fe_0 + tr_z_xxxx_zzzzz[i] * fe_0 + tr_z_xxxx_xzzzzz[i] * pa_x[i];

        tr_z_xxxxx_yyyyyy[i] = 4.0 * tr_z_xxx_yyyyyy[i] * fe_0 + tr_z_xxxx_yyyyyy[i] * pa_x[i];

        tr_z_xxxxx_yyyyyz[i] = 4.0 * tr_z_xxx_yyyyyz[i] * fe_0 + tr_z_xxxx_yyyyyz[i] * pa_x[i];

        tr_z_xxxxx_yyyyzz[i] = 4.0 * tr_z_xxx_yyyyzz[i] * fe_0 + tr_z_xxxx_yyyyzz[i] * pa_x[i];

        tr_z_xxxxx_yyyzzz[i] = 4.0 * tr_z_xxx_yyyzzz[i] * fe_0 + tr_z_xxxx_yyyzzz[i] * pa_x[i];

        tr_z_xxxxx_yyzzzz[i] = 4.0 * tr_z_xxx_yyzzzz[i] * fe_0 + tr_z_xxxx_yyzzzz[i] * pa_x[i];

        tr_z_xxxxx_yzzzzz[i] = 4.0 * tr_z_xxx_yzzzzz[i] * fe_0 + tr_z_xxxx_yzzzzz[i] * pa_x[i];

        tr_z_xxxxx_zzzzzz[i] = 4.0 * tr_z_xxx_zzzzzz[i] * fe_0 + tr_z_xxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 1204-1232 components of targeted buffer : HI

    auto tr_z_xxxxy_xxxxxx = pbuffer.data(idx_dip_hi + 1204);

    auto tr_z_xxxxy_xxxxxy = pbuffer.data(idx_dip_hi + 1205);

    auto tr_z_xxxxy_xxxxxz = pbuffer.data(idx_dip_hi + 1206);

    auto tr_z_xxxxy_xxxxyy = pbuffer.data(idx_dip_hi + 1207);

    auto tr_z_xxxxy_xxxxyz = pbuffer.data(idx_dip_hi + 1208);

    auto tr_z_xxxxy_xxxxzz = pbuffer.data(idx_dip_hi + 1209);

    auto tr_z_xxxxy_xxxyyy = pbuffer.data(idx_dip_hi + 1210);

    auto tr_z_xxxxy_xxxyyz = pbuffer.data(idx_dip_hi + 1211);

    auto tr_z_xxxxy_xxxyzz = pbuffer.data(idx_dip_hi + 1212);

    auto tr_z_xxxxy_xxxzzz = pbuffer.data(idx_dip_hi + 1213);

    auto tr_z_xxxxy_xxyyyy = pbuffer.data(idx_dip_hi + 1214);

    auto tr_z_xxxxy_xxyyyz = pbuffer.data(idx_dip_hi + 1215);

    auto tr_z_xxxxy_xxyyzz = pbuffer.data(idx_dip_hi + 1216);

    auto tr_z_xxxxy_xxyzzz = pbuffer.data(idx_dip_hi + 1217);

    auto tr_z_xxxxy_xxzzzz = pbuffer.data(idx_dip_hi + 1218);

    auto tr_z_xxxxy_xyyyyy = pbuffer.data(idx_dip_hi + 1219);

    auto tr_z_xxxxy_xyyyyz = pbuffer.data(idx_dip_hi + 1220);

    auto tr_z_xxxxy_xyyyzz = pbuffer.data(idx_dip_hi + 1221);

    auto tr_z_xxxxy_xyyzzz = pbuffer.data(idx_dip_hi + 1222);

    auto tr_z_xxxxy_xyzzzz = pbuffer.data(idx_dip_hi + 1223);

    auto tr_z_xxxxy_xzzzzz = pbuffer.data(idx_dip_hi + 1224);

    auto tr_z_xxxxy_yyyyyy = pbuffer.data(idx_dip_hi + 1225);

    auto tr_z_xxxxy_yyyyyz = pbuffer.data(idx_dip_hi + 1226);

    auto tr_z_xxxxy_yyyyzz = pbuffer.data(idx_dip_hi + 1227);

    auto tr_z_xxxxy_yyyzzz = pbuffer.data(idx_dip_hi + 1228);

    auto tr_z_xxxxy_yyzzzz = pbuffer.data(idx_dip_hi + 1229);

    auto tr_z_xxxxy_yzzzzz = pbuffer.data(idx_dip_hi + 1230);

    auto tr_z_xxxxy_zzzzzz = pbuffer.data(idx_dip_hi + 1231);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_z_xxxx_xxxxx,   \
                             tr_z_xxxx_xxxxxx,  \
                             tr_z_xxxx_xxxxxy,  \
                             tr_z_xxxx_xxxxxz,  \
                             tr_z_xxxx_xxxxy,   \
                             tr_z_xxxx_xxxxyy,  \
                             tr_z_xxxx_xxxxyz,  \
                             tr_z_xxxx_xxxxz,   \
                             tr_z_xxxx_xxxxzz,  \
                             tr_z_xxxx_xxxyy,   \
                             tr_z_xxxx_xxxyyy,  \
                             tr_z_xxxx_xxxyyz,  \
                             tr_z_xxxx_xxxyz,   \
                             tr_z_xxxx_xxxyzz,  \
                             tr_z_xxxx_xxxzz,   \
                             tr_z_xxxx_xxxzzz,  \
                             tr_z_xxxx_xxyyy,   \
                             tr_z_xxxx_xxyyyy,  \
                             tr_z_xxxx_xxyyyz,  \
                             tr_z_xxxx_xxyyz,   \
                             tr_z_xxxx_xxyyzz,  \
                             tr_z_xxxx_xxyzz,   \
                             tr_z_xxxx_xxyzzz,  \
                             tr_z_xxxx_xxzzz,   \
                             tr_z_xxxx_xxzzzz,  \
                             tr_z_xxxx_xyyyy,   \
                             tr_z_xxxx_xyyyyy,  \
                             tr_z_xxxx_xyyyyz,  \
                             tr_z_xxxx_xyyyz,   \
                             tr_z_xxxx_xyyyzz,  \
                             tr_z_xxxx_xyyzz,   \
                             tr_z_xxxx_xyyzzz,  \
                             tr_z_xxxx_xyzzz,   \
                             tr_z_xxxx_xyzzzz,  \
                             tr_z_xxxx_xzzzz,   \
                             tr_z_xxxx_xzzzzz,  \
                             tr_z_xxxx_zzzzzz,  \
                             tr_z_xxxxy_xxxxxx, \
                             tr_z_xxxxy_xxxxxy, \
                             tr_z_xxxxy_xxxxxz, \
                             tr_z_xxxxy_xxxxyy, \
                             tr_z_xxxxy_xxxxyz, \
                             tr_z_xxxxy_xxxxzz, \
                             tr_z_xxxxy_xxxyyy, \
                             tr_z_xxxxy_xxxyyz, \
                             tr_z_xxxxy_xxxyzz, \
                             tr_z_xxxxy_xxxzzz, \
                             tr_z_xxxxy_xxyyyy, \
                             tr_z_xxxxy_xxyyyz, \
                             tr_z_xxxxy_xxyyzz, \
                             tr_z_xxxxy_xxyzzz, \
                             tr_z_xxxxy_xxzzzz, \
                             tr_z_xxxxy_xyyyyy, \
                             tr_z_xxxxy_xyyyyz, \
                             tr_z_xxxxy_xyyyzz, \
                             tr_z_xxxxy_xyyzzz, \
                             tr_z_xxxxy_xyzzzz, \
                             tr_z_xxxxy_xzzzzz, \
                             tr_z_xxxxy_yyyyyy, \
                             tr_z_xxxxy_yyyyyz, \
                             tr_z_xxxxy_yyyyzz, \
                             tr_z_xxxxy_yyyzzz, \
                             tr_z_xxxxy_yyzzzz, \
                             tr_z_xxxxy_yzzzzz, \
                             tr_z_xxxxy_zzzzzz, \
                             tr_z_xxxy_yyyyyy,  \
                             tr_z_xxxy_yyyyyz,  \
                             tr_z_xxxy_yyyyzz,  \
                             tr_z_xxxy_yyyzzz,  \
                             tr_z_xxxy_yyzzzz,  \
                             tr_z_xxxy_yzzzzz,  \
                             tr_z_xxy_yyyyyy,   \
                             tr_z_xxy_yyyyyz,   \
                             tr_z_xxy_yyyyzz,   \
                             tr_z_xxy_yyyzzz,   \
                             tr_z_xxy_yyzzzz,   \
                             tr_z_xxy_yzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxy_xxxxxx[i] = tr_z_xxxx_xxxxxx[i] * pa_y[i];

        tr_z_xxxxy_xxxxxy[i] = tr_z_xxxx_xxxxx[i] * fe_0 + tr_z_xxxx_xxxxxy[i] * pa_y[i];

        tr_z_xxxxy_xxxxxz[i] = tr_z_xxxx_xxxxxz[i] * pa_y[i];

        tr_z_xxxxy_xxxxyy[i] = 2.0 * tr_z_xxxx_xxxxy[i] * fe_0 + tr_z_xxxx_xxxxyy[i] * pa_y[i];

        tr_z_xxxxy_xxxxyz[i] = tr_z_xxxx_xxxxz[i] * fe_0 + tr_z_xxxx_xxxxyz[i] * pa_y[i];

        tr_z_xxxxy_xxxxzz[i] = tr_z_xxxx_xxxxzz[i] * pa_y[i];

        tr_z_xxxxy_xxxyyy[i] = 3.0 * tr_z_xxxx_xxxyy[i] * fe_0 + tr_z_xxxx_xxxyyy[i] * pa_y[i];

        tr_z_xxxxy_xxxyyz[i] = 2.0 * tr_z_xxxx_xxxyz[i] * fe_0 + tr_z_xxxx_xxxyyz[i] * pa_y[i];

        tr_z_xxxxy_xxxyzz[i] = tr_z_xxxx_xxxzz[i] * fe_0 + tr_z_xxxx_xxxyzz[i] * pa_y[i];

        tr_z_xxxxy_xxxzzz[i] = tr_z_xxxx_xxxzzz[i] * pa_y[i];

        tr_z_xxxxy_xxyyyy[i] = 4.0 * tr_z_xxxx_xxyyy[i] * fe_0 + tr_z_xxxx_xxyyyy[i] * pa_y[i];

        tr_z_xxxxy_xxyyyz[i] = 3.0 * tr_z_xxxx_xxyyz[i] * fe_0 + tr_z_xxxx_xxyyyz[i] * pa_y[i];

        tr_z_xxxxy_xxyyzz[i] = 2.0 * tr_z_xxxx_xxyzz[i] * fe_0 + tr_z_xxxx_xxyyzz[i] * pa_y[i];

        tr_z_xxxxy_xxyzzz[i] = tr_z_xxxx_xxzzz[i] * fe_0 + tr_z_xxxx_xxyzzz[i] * pa_y[i];

        tr_z_xxxxy_xxzzzz[i] = tr_z_xxxx_xxzzzz[i] * pa_y[i];

        tr_z_xxxxy_xyyyyy[i] = 5.0 * tr_z_xxxx_xyyyy[i] * fe_0 + tr_z_xxxx_xyyyyy[i] * pa_y[i];

        tr_z_xxxxy_xyyyyz[i] = 4.0 * tr_z_xxxx_xyyyz[i] * fe_0 + tr_z_xxxx_xyyyyz[i] * pa_y[i];

        tr_z_xxxxy_xyyyzz[i] = 3.0 * tr_z_xxxx_xyyzz[i] * fe_0 + tr_z_xxxx_xyyyzz[i] * pa_y[i];

        tr_z_xxxxy_xyyzzz[i] = 2.0 * tr_z_xxxx_xyzzz[i] * fe_0 + tr_z_xxxx_xyyzzz[i] * pa_y[i];

        tr_z_xxxxy_xyzzzz[i] = tr_z_xxxx_xzzzz[i] * fe_0 + tr_z_xxxx_xyzzzz[i] * pa_y[i];

        tr_z_xxxxy_xzzzzz[i] = tr_z_xxxx_xzzzzz[i] * pa_y[i];

        tr_z_xxxxy_yyyyyy[i] = 3.0 * tr_z_xxy_yyyyyy[i] * fe_0 + tr_z_xxxy_yyyyyy[i] * pa_x[i];

        tr_z_xxxxy_yyyyyz[i] = 3.0 * tr_z_xxy_yyyyyz[i] * fe_0 + tr_z_xxxy_yyyyyz[i] * pa_x[i];

        tr_z_xxxxy_yyyyzz[i] = 3.0 * tr_z_xxy_yyyyzz[i] * fe_0 + tr_z_xxxy_yyyyzz[i] * pa_x[i];

        tr_z_xxxxy_yyyzzz[i] = 3.0 * tr_z_xxy_yyyzzz[i] * fe_0 + tr_z_xxxy_yyyzzz[i] * pa_x[i];

        tr_z_xxxxy_yyzzzz[i] = 3.0 * tr_z_xxy_yyzzzz[i] * fe_0 + tr_z_xxxy_yyzzzz[i] * pa_x[i];

        tr_z_xxxxy_yzzzzz[i] = 3.0 * tr_z_xxy_yzzzzz[i] * fe_0 + tr_z_xxxy_yzzzzz[i] * pa_x[i];

        tr_z_xxxxy_zzzzzz[i] = tr_z_xxxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 1232-1260 components of targeted buffer : HI

    auto tr_z_xxxxz_xxxxxx = pbuffer.data(idx_dip_hi + 1232);

    auto tr_z_xxxxz_xxxxxy = pbuffer.data(idx_dip_hi + 1233);

    auto tr_z_xxxxz_xxxxxz = pbuffer.data(idx_dip_hi + 1234);

    auto tr_z_xxxxz_xxxxyy = pbuffer.data(idx_dip_hi + 1235);

    auto tr_z_xxxxz_xxxxyz = pbuffer.data(idx_dip_hi + 1236);

    auto tr_z_xxxxz_xxxxzz = pbuffer.data(idx_dip_hi + 1237);

    auto tr_z_xxxxz_xxxyyy = pbuffer.data(idx_dip_hi + 1238);

    auto tr_z_xxxxz_xxxyyz = pbuffer.data(idx_dip_hi + 1239);

    auto tr_z_xxxxz_xxxyzz = pbuffer.data(idx_dip_hi + 1240);

    auto tr_z_xxxxz_xxxzzz = pbuffer.data(idx_dip_hi + 1241);

    auto tr_z_xxxxz_xxyyyy = pbuffer.data(idx_dip_hi + 1242);

    auto tr_z_xxxxz_xxyyyz = pbuffer.data(idx_dip_hi + 1243);

    auto tr_z_xxxxz_xxyyzz = pbuffer.data(idx_dip_hi + 1244);

    auto tr_z_xxxxz_xxyzzz = pbuffer.data(idx_dip_hi + 1245);

    auto tr_z_xxxxz_xxzzzz = pbuffer.data(idx_dip_hi + 1246);

    auto tr_z_xxxxz_xyyyyy = pbuffer.data(idx_dip_hi + 1247);

    auto tr_z_xxxxz_xyyyyz = pbuffer.data(idx_dip_hi + 1248);

    auto tr_z_xxxxz_xyyyzz = pbuffer.data(idx_dip_hi + 1249);

    auto tr_z_xxxxz_xyyzzz = pbuffer.data(idx_dip_hi + 1250);

    auto tr_z_xxxxz_xyzzzz = pbuffer.data(idx_dip_hi + 1251);

    auto tr_z_xxxxz_xzzzzz = pbuffer.data(idx_dip_hi + 1252);

    auto tr_z_xxxxz_yyyyyy = pbuffer.data(idx_dip_hi + 1253);

    auto tr_z_xxxxz_yyyyyz = pbuffer.data(idx_dip_hi + 1254);

    auto tr_z_xxxxz_yyyyzz = pbuffer.data(idx_dip_hi + 1255);

    auto tr_z_xxxxz_yyyzzz = pbuffer.data(idx_dip_hi + 1256);

    auto tr_z_xxxxz_yyzzzz = pbuffer.data(idx_dip_hi + 1257);

    auto tr_z_xxxxz_yzzzzz = pbuffer.data(idx_dip_hi + 1258);

    auto tr_z_xxxxz_zzzzzz = pbuffer.data(idx_dip_hi + 1259);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             tr_z_xxxx_xxxxxx,  \
                             tr_z_xxxx_xxxxxy,  \
                             tr_z_xxxx_xxxxyy,  \
                             tr_z_xxxx_xxxyyy,  \
                             tr_z_xxxx_xxyyyy,  \
                             tr_z_xxxx_xyyyyy,  \
                             tr_z_xxxxz_xxxxxx, \
                             tr_z_xxxxz_xxxxxy, \
                             tr_z_xxxxz_xxxxxz, \
                             tr_z_xxxxz_xxxxyy, \
                             tr_z_xxxxz_xxxxyz, \
                             tr_z_xxxxz_xxxxzz, \
                             tr_z_xxxxz_xxxyyy, \
                             tr_z_xxxxz_xxxyyz, \
                             tr_z_xxxxz_xxxyzz, \
                             tr_z_xxxxz_xxxzzz, \
                             tr_z_xxxxz_xxyyyy, \
                             tr_z_xxxxz_xxyyyz, \
                             tr_z_xxxxz_xxyyzz, \
                             tr_z_xxxxz_xxyzzz, \
                             tr_z_xxxxz_xxzzzz, \
                             tr_z_xxxxz_xyyyyy, \
                             tr_z_xxxxz_xyyyyz, \
                             tr_z_xxxxz_xyyyzz, \
                             tr_z_xxxxz_xyyzzz, \
                             tr_z_xxxxz_xyzzzz, \
                             tr_z_xxxxz_xzzzzz, \
                             tr_z_xxxxz_yyyyyy, \
                             tr_z_xxxxz_yyyyyz, \
                             tr_z_xxxxz_yyyyzz, \
                             tr_z_xxxxz_yyyzzz, \
                             tr_z_xxxxz_yyzzzz, \
                             tr_z_xxxxz_yzzzzz, \
                             tr_z_xxxxz_zzzzzz, \
                             tr_z_xxxz_xxxxxz,  \
                             tr_z_xxxz_xxxxyz,  \
                             tr_z_xxxz_xxxxz,   \
                             tr_z_xxxz_xxxxzz,  \
                             tr_z_xxxz_xxxyyz,  \
                             tr_z_xxxz_xxxyz,   \
                             tr_z_xxxz_xxxyzz,  \
                             tr_z_xxxz_xxxzz,   \
                             tr_z_xxxz_xxxzzz,  \
                             tr_z_xxxz_xxyyyz,  \
                             tr_z_xxxz_xxyyz,   \
                             tr_z_xxxz_xxyyzz,  \
                             tr_z_xxxz_xxyzz,   \
                             tr_z_xxxz_xxyzzz,  \
                             tr_z_xxxz_xxzzz,   \
                             tr_z_xxxz_xxzzzz,  \
                             tr_z_xxxz_xyyyyz,  \
                             tr_z_xxxz_xyyyz,   \
                             tr_z_xxxz_xyyyzz,  \
                             tr_z_xxxz_xyyzz,   \
                             tr_z_xxxz_xyyzzz,  \
                             tr_z_xxxz_xyzzz,   \
                             tr_z_xxxz_xyzzzz,  \
                             tr_z_xxxz_xzzzz,   \
                             tr_z_xxxz_xzzzzz,  \
                             tr_z_xxxz_yyyyyy,  \
                             tr_z_xxxz_yyyyyz,  \
                             tr_z_xxxz_yyyyz,   \
                             tr_z_xxxz_yyyyzz,  \
                             tr_z_xxxz_yyyzz,   \
                             tr_z_xxxz_yyyzzz,  \
                             tr_z_xxxz_yyzzz,   \
                             tr_z_xxxz_yyzzzz,  \
                             tr_z_xxxz_yzzzz,   \
                             tr_z_xxxz_yzzzzz,  \
                             tr_z_xxxz_zzzzz,   \
                             tr_z_xxxz_zzzzzz,  \
                             tr_z_xxz_xxxxxz,   \
                             tr_z_xxz_xxxxyz,   \
                             tr_z_xxz_xxxxzz,   \
                             tr_z_xxz_xxxyyz,   \
                             tr_z_xxz_xxxyzz,   \
                             tr_z_xxz_xxxzzz,   \
                             tr_z_xxz_xxyyyz,   \
                             tr_z_xxz_xxyyzz,   \
                             tr_z_xxz_xxyzzz,   \
                             tr_z_xxz_xxzzzz,   \
                             tr_z_xxz_xyyyyz,   \
                             tr_z_xxz_xyyyzz,   \
                             tr_z_xxz_xyyzzz,   \
                             tr_z_xxz_xyzzzz,   \
                             tr_z_xxz_xzzzzz,   \
                             tr_z_xxz_yyyyyy,   \
                             tr_z_xxz_yyyyyz,   \
                             tr_z_xxz_yyyyzz,   \
                             tr_z_xxz_yyyzzz,   \
                             tr_z_xxz_yyzzzz,   \
                             tr_z_xxz_yzzzzz,   \
                             tr_z_xxz_zzzzzz,   \
                             ts_xxxx_xxxxxx,    \
                             ts_xxxx_xxxxxy,    \
                             ts_xxxx_xxxxyy,    \
                             ts_xxxx_xxxyyy,    \
                             ts_xxxx_xxyyyy,    \
                             ts_xxxx_xyyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxz_xxxxxx[i] = ts_xxxx_xxxxxx[i] * fe_0 + tr_z_xxxx_xxxxxx[i] * pa_z[i];

        tr_z_xxxxz_xxxxxy[i] = ts_xxxx_xxxxxy[i] * fe_0 + tr_z_xxxx_xxxxxy[i] * pa_z[i];

        tr_z_xxxxz_xxxxxz[i] = 3.0 * tr_z_xxz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxxz_xxxxz[i] * fe_0 + tr_z_xxxz_xxxxxz[i] * pa_x[i];

        tr_z_xxxxz_xxxxyy[i] = ts_xxxx_xxxxyy[i] * fe_0 + tr_z_xxxx_xxxxyy[i] * pa_z[i];

        tr_z_xxxxz_xxxxyz[i] = 3.0 * tr_z_xxz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxxz_xxxyz[i] * fe_0 + tr_z_xxxz_xxxxyz[i] * pa_x[i];

        tr_z_xxxxz_xxxxzz[i] = 3.0 * tr_z_xxz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxxz_xxxzz[i] * fe_0 + tr_z_xxxz_xxxxzz[i] * pa_x[i];

        tr_z_xxxxz_xxxyyy[i] = ts_xxxx_xxxyyy[i] * fe_0 + tr_z_xxxx_xxxyyy[i] * pa_z[i];

        tr_z_xxxxz_xxxyyz[i] = 3.0 * tr_z_xxz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxxz_xxyyz[i] * fe_0 + tr_z_xxxz_xxxyyz[i] * pa_x[i];

        tr_z_xxxxz_xxxyzz[i] = 3.0 * tr_z_xxz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxxz_xxyzz[i] * fe_0 + tr_z_xxxz_xxxyzz[i] * pa_x[i];

        tr_z_xxxxz_xxxzzz[i] = 3.0 * tr_z_xxz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxxz_xxzzz[i] * fe_0 + tr_z_xxxz_xxxzzz[i] * pa_x[i];

        tr_z_xxxxz_xxyyyy[i] = ts_xxxx_xxyyyy[i] * fe_0 + tr_z_xxxx_xxyyyy[i] * pa_z[i];

        tr_z_xxxxz_xxyyyz[i] = 3.0 * tr_z_xxz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxxz_xyyyz[i] * fe_0 + tr_z_xxxz_xxyyyz[i] * pa_x[i];

        tr_z_xxxxz_xxyyzz[i] = 3.0 * tr_z_xxz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxxz_xyyzz[i] * fe_0 + tr_z_xxxz_xxyyzz[i] * pa_x[i];

        tr_z_xxxxz_xxyzzz[i] = 3.0 * tr_z_xxz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxxz_xyzzz[i] * fe_0 + tr_z_xxxz_xxyzzz[i] * pa_x[i];

        tr_z_xxxxz_xxzzzz[i] = 3.0 * tr_z_xxz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxxz_xzzzz[i] * fe_0 + tr_z_xxxz_xxzzzz[i] * pa_x[i];

        tr_z_xxxxz_xyyyyy[i] = ts_xxxx_xyyyyy[i] * fe_0 + tr_z_xxxx_xyyyyy[i] * pa_z[i];

        tr_z_xxxxz_xyyyyz[i] = 3.0 * tr_z_xxz_xyyyyz[i] * fe_0 + tr_z_xxxz_yyyyz[i] * fe_0 + tr_z_xxxz_xyyyyz[i] * pa_x[i];

        tr_z_xxxxz_xyyyzz[i] = 3.0 * tr_z_xxz_xyyyzz[i] * fe_0 + tr_z_xxxz_yyyzz[i] * fe_0 + tr_z_xxxz_xyyyzz[i] * pa_x[i];

        tr_z_xxxxz_xyyzzz[i] = 3.0 * tr_z_xxz_xyyzzz[i] * fe_0 + tr_z_xxxz_yyzzz[i] * fe_0 + tr_z_xxxz_xyyzzz[i] * pa_x[i];

        tr_z_xxxxz_xyzzzz[i] = 3.0 * tr_z_xxz_xyzzzz[i] * fe_0 + tr_z_xxxz_yzzzz[i] * fe_0 + tr_z_xxxz_xyzzzz[i] * pa_x[i];

        tr_z_xxxxz_xzzzzz[i] = 3.0 * tr_z_xxz_xzzzzz[i] * fe_0 + tr_z_xxxz_zzzzz[i] * fe_0 + tr_z_xxxz_xzzzzz[i] * pa_x[i];

        tr_z_xxxxz_yyyyyy[i] = 3.0 * tr_z_xxz_yyyyyy[i] * fe_0 + tr_z_xxxz_yyyyyy[i] * pa_x[i];

        tr_z_xxxxz_yyyyyz[i] = 3.0 * tr_z_xxz_yyyyyz[i] * fe_0 + tr_z_xxxz_yyyyyz[i] * pa_x[i];

        tr_z_xxxxz_yyyyzz[i] = 3.0 * tr_z_xxz_yyyyzz[i] * fe_0 + tr_z_xxxz_yyyyzz[i] * pa_x[i];

        tr_z_xxxxz_yyyzzz[i] = 3.0 * tr_z_xxz_yyyzzz[i] * fe_0 + tr_z_xxxz_yyyzzz[i] * pa_x[i];

        tr_z_xxxxz_yyzzzz[i] = 3.0 * tr_z_xxz_yyzzzz[i] * fe_0 + tr_z_xxxz_yyzzzz[i] * pa_x[i];

        tr_z_xxxxz_yzzzzz[i] = 3.0 * tr_z_xxz_yzzzzz[i] * fe_0 + tr_z_xxxz_yzzzzz[i] * pa_x[i];

        tr_z_xxxxz_zzzzzz[i] = 3.0 * tr_z_xxz_zzzzzz[i] * fe_0 + tr_z_xxxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1260-1288 components of targeted buffer : HI

    auto tr_z_xxxyy_xxxxxx = pbuffer.data(idx_dip_hi + 1260);

    auto tr_z_xxxyy_xxxxxy = pbuffer.data(idx_dip_hi + 1261);

    auto tr_z_xxxyy_xxxxxz = pbuffer.data(idx_dip_hi + 1262);

    auto tr_z_xxxyy_xxxxyy = pbuffer.data(idx_dip_hi + 1263);

    auto tr_z_xxxyy_xxxxyz = pbuffer.data(idx_dip_hi + 1264);

    auto tr_z_xxxyy_xxxxzz = pbuffer.data(idx_dip_hi + 1265);

    auto tr_z_xxxyy_xxxyyy = pbuffer.data(idx_dip_hi + 1266);

    auto tr_z_xxxyy_xxxyyz = pbuffer.data(idx_dip_hi + 1267);

    auto tr_z_xxxyy_xxxyzz = pbuffer.data(idx_dip_hi + 1268);

    auto tr_z_xxxyy_xxxzzz = pbuffer.data(idx_dip_hi + 1269);

    auto tr_z_xxxyy_xxyyyy = pbuffer.data(idx_dip_hi + 1270);

    auto tr_z_xxxyy_xxyyyz = pbuffer.data(idx_dip_hi + 1271);

    auto tr_z_xxxyy_xxyyzz = pbuffer.data(idx_dip_hi + 1272);

    auto tr_z_xxxyy_xxyzzz = pbuffer.data(idx_dip_hi + 1273);

    auto tr_z_xxxyy_xxzzzz = pbuffer.data(idx_dip_hi + 1274);

    auto tr_z_xxxyy_xyyyyy = pbuffer.data(idx_dip_hi + 1275);

    auto tr_z_xxxyy_xyyyyz = pbuffer.data(idx_dip_hi + 1276);

    auto tr_z_xxxyy_xyyyzz = pbuffer.data(idx_dip_hi + 1277);

    auto tr_z_xxxyy_xyyzzz = pbuffer.data(idx_dip_hi + 1278);

    auto tr_z_xxxyy_xyzzzz = pbuffer.data(idx_dip_hi + 1279);

    auto tr_z_xxxyy_xzzzzz = pbuffer.data(idx_dip_hi + 1280);

    auto tr_z_xxxyy_yyyyyy = pbuffer.data(idx_dip_hi + 1281);

    auto tr_z_xxxyy_yyyyyz = pbuffer.data(idx_dip_hi + 1282);

    auto tr_z_xxxyy_yyyyzz = pbuffer.data(idx_dip_hi + 1283);

    auto tr_z_xxxyy_yyyzzz = pbuffer.data(idx_dip_hi + 1284);

    auto tr_z_xxxyy_yyzzzz = pbuffer.data(idx_dip_hi + 1285);

    auto tr_z_xxxyy_yzzzzz = pbuffer.data(idx_dip_hi + 1286);

    auto tr_z_xxxyy_zzzzzz = pbuffer.data(idx_dip_hi + 1287);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_z_xxx_xxxxxx,   \
                             tr_z_xxx_xxxxxz,   \
                             tr_z_xxx_xxxxzz,   \
                             tr_z_xxx_xxxzzz,   \
                             tr_z_xxx_xxzzzz,   \
                             tr_z_xxx_xzzzzz,   \
                             tr_z_xxxy_xxxxxx,  \
                             tr_z_xxxy_xxxxxz,  \
                             tr_z_xxxy_xxxxzz,  \
                             tr_z_xxxy_xxxzzz,  \
                             tr_z_xxxy_xxzzzz,  \
                             tr_z_xxxy_xzzzzz,  \
                             tr_z_xxxyy_xxxxxx, \
                             tr_z_xxxyy_xxxxxy, \
                             tr_z_xxxyy_xxxxxz, \
                             tr_z_xxxyy_xxxxyy, \
                             tr_z_xxxyy_xxxxyz, \
                             tr_z_xxxyy_xxxxzz, \
                             tr_z_xxxyy_xxxyyy, \
                             tr_z_xxxyy_xxxyyz, \
                             tr_z_xxxyy_xxxyzz, \
                             tr_z_xxxyy_xxxzzz, \
                             tr_z_xxxyy_xxyyyy, \
                             tr_z_xxxyy_xxyyyz, \
                             tr_z_xxxyy_xxyyzz, \
                             tr_z_xxxyy_xxyzzz, \
                             tr_z_xxxyy_xxzzzz, \
                             tr_z_xxxyy_xyyyyy, \
                             tr_z_xxxyy_xyyyyz, \
                             tr_z_xxxyy_xyyyzz, \
                             tr_z_xxxyy_xyyzzz, \
                             tr_z_xxxyy_xyzzzz, \
                             tr_z_xxxyy_xzzzzz, \
                             tr_z_xxxyy_yyyyyy, \
                             tr_z_xxxyy_yyyyyz, \
                             tr_z_xxxyy_yyyyzz, \
                             tr_z_xxxyy_yyyzzz, \
                             tr_z_xxxyy_yyzzzz, \
                             tr_z_xxxyy_yzzzzz, \
                             tr_z_xxxyy_zzzzzz, \
                             tr_z_xxyy_xxxxxy,  \
                             tr_z_xxyy_xxxxy,   \
                             tr_z_xxyy_xxxxyy,  \
                             tr_z_xxyy_xxxxyz,  \
                             tr_z_xxyy_xxxyy,   \
                             tr_z_xxyy_xxxyyy,  \
                             tr_z_xxyy_xxxyyz,  \
                             tr_z_xxyy_xxxyz,   \
                             tr_z_xxyy_xxxyzz,  \
                             tr_z_xxyy_xxyyy,   \
                             tr_z_xxyy_xxyyyy,  \
                             tr_z_xxyy_xxyyyz,  \
                             tr_z_xxyy_xxyyz,   \
                             tr_z_xxyy_xxyyzz,  \
                             tr_z_xxyy_xxyzz,   \
                             tr_z_xxyy_xxyzzz,  \
                             tr_z_xxyy_xyyyy,   \
                             tr_z_xxyy_xyyyyy,  \
                             tr_z_xxyy_xyyyyz,  \
                             tr_z_xxyy_xyyyz,   \
                             tr_z_xxyy_xyyyzz,  \
                             tr_z_xxyy_xyyzz,   \
                             tr_z_xxyy_xyyzzz,  \
                             tr_z_xxyy_xyzzz,   \
                             tr_z_xxyy_xyzzzz,  \
                             tr_z_xxyy_yyyyy,   \
                             tr_z_xxyy_yyyyyy,  \
                             tr_z_xxyy_yyyyyz,  \
                             tr_z_xxyy_yyyyz,   \
                             tr_z_xxyy_yyyyzz,  \
                             tr_z_xxyy_yyyzz,   \
                             tr_z_xxyy_yyyzzz,  \
                             tr_z_xxyy_yyzzz,   \
                             tr_z_xxyy_yyzzzz,  \
                             tr_z_xxyy_yzzzz,   \
                             tr_z_xxyy_yzzzzz,  \
                             tr_z_xxyy_zzzzzz,  \
                             tr_z_xyy_xxxxxy,   \
                             tr_z_xyy_xxxxyy,   \
                             tr_z_xyy_xxxxyz,   \
                             tr_z_xyy_xxxyyy,   \
                             tr_z_xyy_xxxyyz,   \
                             tr_z_xyy_xxxyzz,   \
                             tr_z_xyy_xxyyyy,   \
                             tr_z_xyy_xxyyyz,   \
                             tr_z_xyy_xxyyzz,   \
                             tr_z_xyy_xxyzzz,   \
                             tr_z_xyy_xyyyyy,   \
                             tr_z_xyy_xyyyyz,   \
                             tr_z_xyy_xyyyzz,   \
                             tr_z_xyy_xyyzzz,   \
                             tr_z_xyy_xyzzzz,   \
                             tr_z_xyy_yyyyyy,   \
                             tr_z_xyy_yyyyyz,   \
                             tr_z_xyy_yyyyzz,   \
                             tr_z_xyy_yyyzzz,   \
                             tr_z_xyy_yyzzzz,   \
                             tr_z_xyy_yzzzzz,   \
                             tr_z_xyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyy_xxxxxx[i] = tr_z_xxx_xxxxxx[i] * fe_0 + tr_z_xxxy_xxxxxx[i] * pa_y[i];

        tr_z_xxxyy_xxxxxy[i] = 2.0 * tr_z_xyy_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxyy_xxxxy[i] * fe_0 + tr_z_xxyy_xxxxxy[i] * pa_x[i];

        tr_z_xxxyy_xxxxxz[i] = tr_z_xxx_xxxxxz[i] * fe_0 + tr_z_xxxy_xxxxxz[i] * pa_y[i];

        tr_z_xxxyy_xxxxyy[i] = 2.0 * tr_z_xyy_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxyy_xxxyy[i] * fe_0 + tr_z_xxyy_xxxxyy[i] * pa_x[i];

        tr_z_xxxyy_xxxxyz[i] = 2.0 * tr_z_xyy_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxyy_xxxyz[i] * fe_0 + tr_z_xxyy_xxxxyz[i] * pa_x[i];

        tr_z_xxxyy_xxxxzz[i] = tr_z_xxx_xxxxzz[i] * fe_0 + tr_z_xxxy_xxxxzz[i] * pa_y[i];

        tr_z_xxxyy_xxxyyy[i] = 2.0 * tr_z_xyy_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxyy_xxyyy[i] * fe_0 + tr_z_xxyy_xxxyyy[i] * pa_x[i];

        tr_z_xxxyy_xxxyyz[i] = 2.0 * tr_z_xyy_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxyy_xxyyz[i] * fe_0 + tr_z_xxyy_xxxyyz[i] * pa_x[i];

        tr_z_xxxyy_xxxyzz[i] = 2.0 * tr_z_xyy_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxyy_xxyzz[i] * fe_0 + tr_z_xxyy_xxxyzz[i] * pa_x[i];

        tr_z_xxxyy_xxxzzz[i] = tr_z_xxx_xxxzzz[i] * fe_0 + tr_z_xxxy_xxxzzz[i] * pa_y[i];

        tr_z_xxxyy_xxyyyy[i] = 2.0 * tr_z_xyy_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxyy_xyyyy[i] * fe_0 + tr_z_xxyy_xxyyyy[i] * pa_x[i];

        tr_z_xxxyy_xxyyyz[i] = 2.0 * tr_z_xyy_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxyy_xyyyz[i] * fe_0 + tr_z_xxyy_xxyyyz[i] * pa_x[i];

        tr_z_xxxyy_xxyyzz[i] = 2.0 * tr_z_xyy_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxyy_xyyzz[i] * fe_0 + tr_z_xxyy_xxyyzz[i] * pa_x[i];

        tr_z_xxxyy_xxyzzz[i] = 2.0 * tr_z_xyy_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxyy_xyzzz[i] * fe_0 + tr_z_xxyy_xxyzzz[i] * pa_x[i];

        tr_z_xxxyy_xxzzzz[i] = tr_z_xxx_xxzzzz[i] * fe_0 + tr_z_xxxy_xxzzzz[i] * pa_y[i];

        tr_z_xxxyy_xyyyyy[i] = 2.0 * tr_z_xyy_xyyyyy[i] * fe_0 + tr_z_xxyy_yyyyy[i] * fe_0 + tr_z_xxyy_xyyyyy[i] * pa_x[i];

        tr_z_xxxyy_xyyyyz[i] = 2.0 * tr_z_xyy_xyyyyz[i] * fe_0 + tr_z_xxyy_yyyyz[i] * fe_0 + tr_z_xxyy_xyyyyz[i] * pa_x[i];

        tr_z_xxxyy_xyyyzz[i] = 2.0 * tr_z_xyy_xyyyzz[i] * fe_0 + tr_z_xxyy_yyyzz[i] * fe_0 + tr_z_xxyy_xyyyzz[i] * pa_x[i];

        tr_z_xxxyy_xyyzzz[i] = 2.0 * tr_z_xyy_xyyzzz[i] * fe_0 + tr_z_xxyy_yyzzz[i] * fe_0 + tr_z_xxyy_xyyzzz[i] * pa_x[i];

        tr_z_xxxyy_xyzzzz[i] = 2.0 * tr_z_xyy_xyzzzz[i] * fe_0 + tr_z_xxyy_yzzzz[i] * fe_0 + tr_z_xxyy_xyzzzz[i] * pa_x[i];

        tr_z_xxxyy_xzzzzz[i] = tr_z_xxx_xzzzzz[i] * fe_0 + tr_z_xxxy_xzzzzz[i] * pa_y[i];

        tr_z_xxxyy_yyyyyy[i] = 2.0 * tr_z_xyy_yyyyyy[i] * fe_0 + tr_z_xxyy_yyyyyy[i] * pa_x[i];

        tr_z_xxxyy_yyyyyz[i] = 2.0 * tr_z_xyy_yyyyyz[i] * fe_0 + tr_z_xxyy_yyyyyz[i] * pa_x[i];

        tr_z_xxxyy_yyyyzz[i] = 2.0 * tr_z_xyy_yyyyzz[i] * fe_0 + tr_z_xxyy_yyyyzz[i] * pa_x[i];

        tr_z_xxxyy_yyyzzz[i] = 2.0 * tr_z_xyy_yyyzzz[i] * fe_0 + tr_z_xxyy_yyyzzz[i] * pa_x[i];

        tr_z_xxxyy_yyzzzz[i] = 2.0 * tr_z_xyy_yyzzzz[i] * fe_0 + tr_z_xxyy_yyzzzz[i] * pa_x[i];

        tr_z_xxxyy_yzzzzz[i] = 2.0 * tr_z_xyy_yzzzzz[i] * fe_0 + tr_z_xxyy_yzzzzz[i] * pa_x[i];

        tr_z_xxxyy_zzzzzz[i] = 2.0 * tr_z_xyy_zzzzzz[i] * fe_0 + tr_z_xxyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1288-1316 components of targeted buffer : HI

    auto tr_z_xxxyz_xxxxxx = pbuffer.data(idx_dip_hi + 1288);

    auto tr_z_xxxyz_xxxxxy = pbuffer.data(idx_dip_hi + 1289);

    auto tr_z_xxxyz_xxxxxz = pbuffer.data(idx_dip_hi + 1290);

    auto tr_z_xxxyz_xxxxyy = pbuffer.data(idx_dip_hi + 1291);

    auto tr_z_xxxyz_xxxxyz = pbuffer.data(idx_dip_hi + 1292);

    auto tr_z_xxxyz_xxxxzz = pbuffer.data(idx_dip_hi + 1293);

    auto tr_z_xxxyz_xxxyyy = pbuffer.data(idx_dip_hi + 1294);

    auto tr_z_xxxyz_xxxyyz = pbuffer.data(idx_dip_hi + 1295);

    auto tr_z_xxxyz_xxxyzz = pbuffer.data(idx_dip_hi + 1296);

    auto tr_z_xxxyz_xxxzzz = pbuffer.data(idx_dip_hi + 1297);

    auto tr_z_xxxyz_xxyyyy = pbuffer.data(idx_dip_hi + 1298);

    auto tr_z_xxxyz_xxyyyz = pbuffer.data(idx_dip_hi + 1299);

    auto tr_z_xxxyz_xxyyzz = pbuffer.data(idx_dip_hi + 1300);

    auto tr_z_xxxyz_xxyzzz = pbuffer.data(idx_dip_hi + 1301);

    auto tr_z_xxxyz_xxzzzz = pbuffer.data(idx_dip_hi + 1302);

    auto tr_z_xxxyz_xyyyyy = pbuffer.data(idx_dip_hi + 1303);

    auto tr_z_xxxyz_xyyyyz = pbuffer.data(idx_dip_hi + 1304);

    auto tr_z_xxxyz_xyyyzz = pbuffer.data(idx_dip_hi + 1305);

    auto tr_z_xxxyz_xyyzzz = pbuffer.data(idx_dip_hi + 1306);

    auto tr_z_xxxyz_xyzzzz = pbuffer.data(idx_dip_hi + 1307);

    auto tr_z_xxxyz_xzzzzz = pbuffer.data(idx_dip_hi + 1308);

    auto tr_z_xxxyz_yyyyyy = pbuffer.data(idx_dip_hi + 1309);

    auto tr_z_xxxyz_yyyyyz = pbuffer.data(idx_dip_hi + 1310);

    auto tr_z_xxxyz_yyyyzz = pbuffer.data(idx_dip_hi + 1311);

    auto tr_z_xxxyz_yyyzzz = pbuffer.data(idx_dip_hi + 1312);

    auto tr_z_xxxyz_yyzzzz = pbuffer.data(idx_dip_hi + 1313);

    auto tr_z_xxxyz_yzzzzz = pbuffer.data(idx_dip_hi + 1314);

    auto tr_z_xxxyz_zzzzzz = pbuffer.data(idx_dip_hi + 1315);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_z_xxxyz_xxxxxx, \
                             tr_z_xxxyz_xxxxxy, \
                             tr_z_xxxyz_xxxxxz, \
                             tr_z_xxxyz_xxxxyy, \
                             tr_z_xxxyz_xxxxyz, \
                             tr_z_xxxyz_xxxxzz, \
                             tr_z_xxxyz_xxxyyy, \
                             tr_z_xxxyz_xxxyyz, \
                             tr_z_xxxyz_xxxyzz, \
                             tr_z_xxxyz_xxxzzz, \
                             tr_z_xxxyz_xxyyyy, \
                             tr_z_xxxyz_xxyyyz, \
                             tr_z_xxxyz_xxyyzz, \
                             tr_z_xxxyz_xxyzzz, \
                             tr_z_xxxyz_xxzzzz, \
                             tr_z_xxxyz_xyyyyy, \
                             tr_z_xxxyz_xyyyyz, \
                             tr_z_xxxyz_xyyyzz, \
                             tr_z_xxxyz_xyyzzz, \
                             tr_z_xxxyz_xyzzzz, \
                             tr_z_xxxyz_xzzzzz, \
                             tr_z_xxxyz_yyyyyy, \
                             tr_z_xxxyz_yyyyyz, \
                             tr_z_xxxyz_yyyyzz, \
                             tr_z_xxxyz_yyyzzz, \
                             tr_z_xxxyz_yyzzzz, \
                             tr_z_xxxyz_yzzzzz, \
                             tr_z_xxxyz_zzzzzz, \
                             tr_z_xxxz_xxxxx,   \
                             tr_z_xxxz_xxxxxx,  \
                             tr_z_xxxz_xxxxxy,  \
                             tr_z_xxxz_xxxxxz,  \
                             tr_z_xxxz_xxxxy,   \
                             tr_z_xxxz_xxxxyy,  \
                             tr_z_xxxz_xxxxyz,  \
                             tr_z_xxxz_xxxxz,   \
                             tr_z_xxxz_xxxxzz,  \
                             tr_z_xxxz_xxxyy,   \
                             tr_z_xxxz_xxxyyy,  \
                             tr_z_xxxz_xxxyyz,  \
                             tr_z_xxxz_xxxyz,   \
                             tr_z_xxxz_xxxyzz,  \
                             tr_z_xxxz_xxxzz,   \
                             tr_z_xxxz_xxxzzz,  \
                             tr_z_xxxz_xxyyy,   \
                             tr_z_xxxz_xxyyyy,  \
                             tr_z_xxxz_xxyyyz,  \
                             tr_z_xxxz_xxyyz,   \
                             tr_z_xxxz_xxyyzz,  \
                             tr_z_xxxz_xxyzz,   \
                             tr_z_xxxz_xxyzzz,  \
                             tr_z_xxxz_xxzzz,   \
                             tr_z_xxxz_xxzzzz,  \
                             tr_z_xxxz_xyyyy,   \
                             tr_z_xxxz_xyyyyy,  \
                             tr_z_xxxz_xyyyyz,  \
                             tr_z_xxxz_xyyyz,   \
                             tr_z_xxxz_xyyyzz,  \
                             tr_z_xxxz_xyyzz,   \
                             tr_z_xxxz_xyyzzz,  \
                             tr_z_xxxz_xyzzz,   \
                             tr_z_xxxz_xyzzzz,  \
                             tr_z_xxxz_xzzzz,   \
                             tr_z_xxxz_xzzzzz,  \
                             tr_z_xxxz_zzzzzz,  \
                             tr_z_xxyz_yyyyyy,  \
                             tr_z_xxyz_yyyyyz,  \
                             tr_z_xxyz_yyyyzz,  \
                             tr_z_xxyz_yyyzzz,  \
                             tr_z_xxyz_yyzzzz,  \
                             tr_z_xxyz_yzzzzz,  \
                             tr_z_xyz_yyyyyy,   \
                             tr_z_xyz_yyyyyz,   \
                             tr_z_xyz_yyyyzz,   \
                             tr_z_xyz_yyyzzz,   \
                             tr_z_xyz_yyzzzz,   \
                             tr_z_xyz_yzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyz_xxxxxx[i] = tr_z_xxxz_xxxxxx[i] * pa_y[i];

        tr_z_xxxyz_xxxxxy[i] = tr_z_xxxz_xxxxx[i] * fe_0 + tr_z_xxxz_xxxxxy[i] * pa_y[i];

        tr_z_xxxyz_xxxxxz[i] = tr_z_xxxz_xxxxxz[i] * pa_y[i];

        tr_z_xxxyz_xxxxyy[i] = 2.0 * tr_z_xxxz_xxxxy[i] * fe_0 + tr_z_xxxz_xxxxyy[i] * pa_y[i];

        tr_z_xxxyz_xxxxyz[i] = tr_z_xxxz_xxxxz[i] * fe_0 + tr_z_xxxz_xxxxyz[i] * pa_y[i];

        tr_z_xxxyz_xxxxzz[i] = tr_z_xxxz_xxxxzz[i] * pa_y[i];

        tr_z_xxxyz_xxxyyy[i] = 3.0 * tr_z_xxxz_xxxyy[i] * fe_0 + tr_z_xxxz_xxxyyy[i] * pa_y[i];

        tr_z_xxxyz_xxxyyz[i] = 2.0 * tr_z_xxxz_xxxyz[i] * fe_0 + tr_z_xxxz_xxxyyz[i] * pa_y[i];

        tr_z_xxxyz_xxxyzz[i] = tr_z_xxxz_xxxzz[i] * fe_0 + tr_z_xxxz_xxxyzz[i] * pa_y[i];

        tr_z_xxxyz_xxxzzz[i] = tr_z_xxxz_xxxzzz[i] * pa_y[i];

        tr_z_xxxyz_xxyyyy[i] = 4.0 * tr_z_xxxz_xxyyy[i] * fe_0 + tr_z_xxxz_xxyyyy[i] * pa_y[i];

        tr_z_xxxyz_xxyyyz[i] = 3.0 * tr_z_xxxz_xxyyz[i] * fe_0 + tr_z_xxxz_xxyyyz[i] * pa_y[i];

        tr_z_xxxyz_xxyyzz[i] = 2.0 * tr_z_xxxz_xxyzz[i] * fe_0 + tr_z_xxxz_xxyyzz[i] * pa_y[i];

        tr_z_xxxyz_xxyzzz[i] = tr_z_xxxz_xxzzz[i] * fe_0 + tr_z_xxxz_xxyzzz[i] * pa_y[i];

        tr_z_xxxyz_xxzzzz[i] = tr_z_xxxz_xxzzzz[i] * pa_y[i];

        tr_z_xxxyz_xyyyyy[i] = 5.0 * tr_z_xxxz_xyyyy[i] * fe_0 + tr_z_xxxz_xyyyyy[i] * pa_y[i];

        tr_z_xxxyz_xyyyyz[i] = 4.0 * tr_z_xxxz_xyyyz[i] * fe_0 + tr_z_xxxz_xyyyyz[i] * pa_y[i];

        tr_z_xxxyz_xyyyzz[i] = 3.0 * tr_z_xxxz_xyyzz[i] * fe_0 + tr_z_xxxz_xyyyzz[i] * pa_y[i];

        tr_z_xxxyz_xyyzzz[i] = 2.0 * tr_z_xxxz_xyzzz[i] * fe_0 + tr_z_xxxz_xyyzzz[i] * pa_y[i];

        tr_z_xxxyz_xyzzzz[i] = tr_z_xxxz_xzzzz[i] * fe_0 + tr_z_xxxz_xyzzzz[i] * pa_y[i];

        tr_z_xxxyz_xzzzzz[i] = tr_z_xxxz_xzzzzz[i] * pa_y[i];

        tr_z_xxxyz_yyyyyy[i] = 2.0 * tr_z_xyz_yyyyyy[i] * fe_0 + tr_z_xxyz_yyyyyy[i] * pa_x[i];

        tr_z_xxxyz_yyyyyz[i] = 2.0 * tr_z_xyz_yyyyyz[i] * fe_0 + tr_z_xxyz_yyyyyz[i] * pa_x[i];

        tr_z_xxxyz_yyyyzz[i] = 2.0 * tr_z_xyz_yyyyzz[i] * fe_0 + tr_z_xxyz_yyyyzz[i] * pa_x[i];

        tr_z_xxxyz_yyyzzz[i] = 2.0 * tr_z_xyz_yyyzzz[i] * fe_0 + tr_z_xxyz_yyyzzz[i] * pa_x[i];

        tr_z_xxxyz_yyzzzz[i] = 2.0 * tr_z_xyz_yyzzzz[i] * fe_0 + tr_z_xxyz_yyzzzz[i] * pa_x[i];

        tr_z_xxxyz_yzzzzz[i] = 2.0 * tr_z_xyz_yzzzzz[i] * fe_0 + tr_z_xxyz_yzzzzz[i] * pa_x[i];

        tr_z_xxxyz_zzzzzz[i] = tr_z_xxxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1316-1344 components of targeted buffer : HI

    auto tr_z_xxxzz_xxxxxx = pbuffer.data(idx_dip_hi + 1316);

    auto tr_z_xxxzz_xxxxxy = pbuffer.data(idx_dip_hi + 1317);

    auto tr_z_xxxzz_xxxxxz = pbuffer.data(idx_dip_hi + 1318);

    auto tr_z_xxxzz_xxxxyy = pbuffer.data(idx_dip_hi + 1319);

    auto tr_z_xxxzz_xxxxyz = pbuffer.data(idx_dip_hi + 1320);

    auto tr_z_xxxzz_xxxxzz = pbuffer.data(idx_dip_hi + 1321);

    auto tr_z_xxxzz_xxxyyy = pbuffer.data(idx_dip_hi + 1322);

    auto tr_z_xxxzz_xxxyyz = pbuffer.data(idx_dip_hi + 1323);

    auto tr_z_xxxzz_xxxyzz = pbuffer.data(idx_dip_hi + 1324);

    auto tr_z_xxxzz_xxxzzz = pbuffer.data(idx_dip_hi + 1325);

    auto tr_z_xxxzz_xxyyyy = pbuffer.data(idx_dip_hi + 1326);

    auto tr_z_xxxzz_xxyyyz = pbuffer.data(idx_dip_hi + 1327);

    auto tr_z_xxxzz_xxyyzz = pbuffer.data(idx_dip_hi + 1328);

    auto tr_z_xxxzz_xxyzzz = pbuffer.data(idx_dip_hi + 1329);

    auto tr_z_xxxzz_xxzzzz = pbuffer.data(idx_dip_hi + 1330);

    auto tr_z_xxxzz_xyyyyy = pbuffer.data(idx_dip_hi + 1331);

    auto tr_z_xxxzz_xyyyyz = pbuffer.data(idx_dip_hi + 1332);

    auto tr_z_xxxzz_xyyyzz = pbuffer.data(idx_dip_hi + 1333);

    auto tr_z_xxxzz_xyyzzz = pbuffer.data(idx_dip_hi + 1334);

    auto tr_z_xxxzz_xyzzzz = pbuffer.data(idx_dip_hi + 1335);

    auto tr_z_xxxzz_xzzzzz = pbuffer.data(idx_dip_hi + 1336);

    auto tr_z_xxxzz_yyyyyy = pbuffer.data(idx_dip_hi + 1337);

    auto tr_z_xxxzz_yyyyyz = pbuffer.data(idx_dip_hi + 1338);

    auto tr_z_xxxzz_yyyyzz = pbuffer.data(idx_dip_hi + 1339);

    auto tr_z_xxxzz_yyyzzz = pbuffer.data(idx_dip_hi + 1340);

    auto tr_z_xxxzz_yyzzzz = pbuffer.data(idx_dip_hi + 1341);

    auto tr_z_xxxzz_yzzzzz = pbuffer.data(idx_dip_hi + 1342);

    auto tr_z_xxxzz_zzzzzz = pbuffer.data(idx_dip_hi + 1343);

#pragma omp simd aligned(pa_x,                  \
                             tr_z_xxxzz_xxxxxx, \
                             tr_z_xxxzz_xxxxxy, \
                             tr_z_xxxzz_xxxxxz, \
                             tr_z_xxxzz_xxxxyy, \
                             tr_z_xxxzz_xxxxyz, \
                             tr_z_xxxzz_xxxxzz, \
                             tr_z_xxxzz_xxxyyy, \
                             tr_z_xxxzz_xxxyyz, \
                             tr_z_xxxzz_xxxyzz, \
                             tr_z_xxxzz_xxxzzz, \
                             tr_z_xxxzz_xxyyyy, \
                             tr_z_xxxzz_xxyyyz, \
                             tr_z_xxxzz_xxyyzz, \
                             tr_z_xxxzz_xxyzzz, \
                             tr_z_xxxzz_xxzzzz, \
                             tr_z_xxxzz_xyyyyy, \
                             tr_z_xxxzz_xyyyyz, \
                             tr_z_xxxzz_xyyyzz, \
                             tr_z_xxxzz_xyyzzz, \
                             tr_z_xxxzz_xyzzzz, \
                             tr_z_xxxzz_xzzzzz, \
                             tr_z_xxxzz_yyyyyy, \
                             tr_z_xxxzz_yyyyyz, \
                             tr_z_xxxzz_yyyyzz, \
                             tr_z_xxxzz_yyyzzz, \
                             tr_z_xxxzz_yyzzzz, \
                             tr_z_xxxzz_yzzzzz, \
                             tr_z_xxxzz_zzzzzz, \
                             tr_z_xxzz_xxxxx,   \
                             tr_z_xxzz_xxxxxx,  \
                             tr_z_xxzz_xxxxxy,  \
                             tr_z_xxzz_xxxxxz,  \
                             tr_z_xxzz_xxxxy,   \
                             tr_z_xxzz_xxxxyy,  \
                             tr_z_xxzz_xxxxyz,  \
                             tr_z_xxzz_xxxxz,   \
                             tr_z_xxzz_xxxxzz,  \
                             tr_z_xxzz_xxxyy,   \
                             tr_z_xxzz_xxxyyy,  \
                             tr_z_xxzz_xxxyyz,  \
                             tr_z_xxzz_xxxyz,   \
                             tr_z_xxzz_xxxyzz,  \
                             tr_z_xxzz_xxxzz,   \
                             tr_z_xxzz_xxxzzz,  \
                             tr_z_xxzz_xxyyy,   \
                             tr_z_xxzz_xxyyyy,  \
                             tr_z_xxzz_xxyyyz,  \
                             tr_z_xxzz_xxyyz,   \
                             tr_z_xxzz_xxyyzz,  \
                             tr_z_xxzz_xxyzz,   \
                             tr_z_xxzz_xxyzzz,  \
                             tr_z_xxzz_xxzzz,   \
                             tr_z_xxzz_xxzzzz,  \
                             tr_z_xxzz_xyyyy,   \
                             tr_z_xxzz_xyyyyy,  \
                             tr_z_xxzz_xyyyyz,  \
                             tr_z_xxzz_xyyyz,   \
                             tr_z_xxzz_xyyyzz,  \
                             tr_z_xxzz_xyyzz,   \
                             tr_z_xxzz_xyyzzz,  \
                             tr_z_xxzz_xyzzz,   \
                             tr_z_xxzz_xyzzzz,  \
                             tr_z_xxzz_xzzzz,   \
                             tr_z_xxzz_xzzzzz,  \
                             tr_z_xxzz_yyyyy,   \
                             tr_z_xxzz_yyyyyy,  \
                             tr_z_xxzz_yyyyyz,  \
                             tr_z_xxzz_yyyyz,   \
                             tr_z_xxzz_yyyyzz,  \
                             tr_z_xxzz_yyyzz,   \
                             tr_z_xxzz_yyyzzz,  \
                             tr_z_xxzz_yyzzz,   \
                             tr_z_xxzz_yyzzzz,  \
                             tr_z_xxzz_yzzzz,   \
                             tr_z_xxzz_yzzzzz,  \
                             tr_z_xxzz_zzzzz,   \
                             tr_z_xxzz_zzzzzz,  \
                             tr_z_xzz_xxxxxx,   \
                             tr_z_xzz_xxxxxy,   \
                             tr_z_xzz_xxxxxz,   \
                             tr_z_xzz_xxxxyy,   \
                             tr_z_xzz_xxxxyz,   \
                             tr_z_xzz_xxxxzz,   \
                             tr_z_xzz_xxxyyy,   \
                             tr_z_xzz_xxxyyz,   \
                             tr_z_xzz_xxxyzz,   \
                             tr_z_xzz_xxxzzz,   \
                             tr_z_xzz_xxyyyy,   \
                             tr_z_xzz_xxyyyz,   \
                             tr_z_xzz_xxyyzz,   \
                             tr_z_xzz_xxyzzz,   \
                             tr_z_xzz_xxzzzz,   \
                             tr_z_xzz_xyyyyy,   \
                             tr_z_xzz_xyyyyz,   \
                             tr_z_xzz_xyyyzz,   \
                             tr_z_xzz_xyyzzz,   \
                             tr_z_xzz_xyzzzz,   \
                             tr_z_xzz_xzzzzz,   \
                             tr_z_xzz_yyyyyy,   \
                             tr_z_xzz_yyyyyz,   \
                             tr_z_xzz_yyyyzz,   \
                             tr_z_xzz_yyyzzz,   \
                             tr_z_xzz_yyzzzz,   \
                             tr_z_xzz_yzzzzz,   \
                             tr_z_xzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzz_xxxxxx[i] = 2.0 * tr_z_xzz_xxxxxx[i] * fe_0 + 6.0 * tr_z_xxzz_xxxxx[i] * fe_0 + tr_z_xxzz_xxxxxx[i] * pa_x[i];

        tr_z_xxxzz_xxxxxy[i] = 2.0 * tr_z_xzz_xxxxxy[i] * fe_0 + 5.0 * tr_z_xxzz_xxxxy[i] * fe_0 + tr_z_xxzz_xxxxxy[i] * pa_x[i];

        tr_z_xxxzz_xxxxxz[i] = 2.0 * tr_z_xzz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xxzz_xxxxz[i] * fe_0 + tr_z_xxzz_xxxxxz[i] * pa_x[i];

        tr_z_xxxzz_xxxxyy[i] = 2.0 * tr_z_xzz_xxxxyy[i] * fe_0 + 4.0 * tr_z_xxzz_xxxyy[i] * fe_0 + tr_z_xxzz_xxxxyy[i] * pa_x[i];

        tr_z_xxxzz_xxxxyz[i] = 2.0 * tr_z_xzz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xxzz_xxxyz[i] * fe_0 + tr_z_xxzz_xxxxyz[i] * pa_x[i];

        tr_z_xxxzz_xxxxzz[i] = 2.0 * tr_z_xzz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xxzz_xxxzz[i] * fe_0 + tr_z_xxzz_xxxxzz[i] * pa_x[i];

        tr_z_xxxzz_xxxyyy[i] = 2.0 * tr_z_xzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_xxzz_xxyyy[i] * fe_0 + tr_z_xxzz_xxxyyy[i] * pa_x[i];

        tr_z_xxxzz_xxxyyz[i] = 2.0 * tr_z_xzz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xxzz_xxyyz[i] * fe_0 + tr_z_xxzz_xxxyyz[i] * pa_x[i];

        tr_z_xxxzz_xxxyzz[i] = 2.0 * tr_z_xzz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xxzz_xxyzz[i] * fe_0 + tr_z_xxzz_xxxyzz[i] * pa_x[i];

        tr_z_xxxzz_xxxzzz[i] = 2.0 * tr_z_xzz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xxzz_xxzzz[i] * fe_0 + tr_z_xxzz_xxxzzz[i] * pa_x[i];

        tr_z_xxxzz_xxyyyy[i] = 2.0 * tr_z_xzz_xxyyyy[i] * fe_0 + 2.0 * tr_z_xxzz_xyyyy[i] * fe_0 + tr_z_xxzz_xxyyyy[i] * pa_x[i];

        tr_z_xxxzz_xxyyyz[i] = 2.0 * tr_z_xzz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xxzz_xyyyz[i] * fe_0 + tr_z_xxzz_xxyyyz[i] * pa_x[i];

        tr_z_xxxzz_xxyyzz[i] = 2.0 * tr_z_xzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xxzz_xyyzz[i] * fe_0 + tr_z_xxzz_xxyyzz[i] * pa_x[i];

        tr_z_xxxzz_xxyzzz[i] = 2.0 * tr_z_xzz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xxzz_xyzzz[i] * fe_0 + tr_z_xxzz_xxyzzz[i] * pa_x[i];

        tr_z_xxxzz_xxzzzz[i] = 2.0 * tr_z_xzz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xxzz_xzzzz[i] * fe_0 + tr_z_xxzz_xxzzzz[i] * pa_x[i];

        tr_z_xxxzz_xyyyyy[i] = 2.0 * tr_z_xzz_xyyyyy[i] * fe_0 + tr_z_xxzz_yyyyy[i] * fe_0 + tr_z_xxzz_xyyyyy[i] * pa_x[i];

        tr_z_xxxzz_xyyyyz[i] = 2.0 * tr_z_xzz_xyyyyz[i] * fe_0 + tr_z_xxzz_yyyyz[i] * fe_0 + tr_z_xxzz_xyyyyz[i] * pa_x[i];

        tr_z_xxxzz_xyyyzz[i] = 2.0 * tr_z_xzz_xyyyzz[i] * fe_0 + tr_z_xxzz_yyyzz[i] * fe_0 + tr_z_xxzz_xyyyzz[i] * pa_x[i];

        tr_z_xxxzz_xyyzzz[i] = 2.0 * tr_z_xzz_xyyzzz[i] * fe_0 + tr_z_xxzz_yyzzz[i] * fe_0 + tr_z_xxzz_xyyzzz[i] * pa_x[i];

        tr_z_xxxzz_xyzzzz[i] = 2.0 * tr_z_xzz_xyzzzz[i] * fe_0 + tr_z_xxzz_yzzzz[i] * fe_0 + tr_z_xxzz_xyzzzz[i] * pa_x[i];

        tr_z_xxxzz_xzzzzz[i] = 2.0 * tr_z_xzz_xzzzzz[i] * fe_0 + tr_z_xxzz_zzzzz[i] * fe_0 + tr_z_xxzz_xzzzzz[i] * pa_x[i];

        tr_z_xxxzz_yyyyyy[i] = 2.0 * tr_z_xzz_yyyyyy[i] * fe_0 + tr_z_xxzz_yyyyyy[i] * pa_x[i];

        tr_z_xxxzz_yyyyyz[i] = 2.0 * tr_z_xzz_yyyyyz[i] * fe_0 + tr_z_xxzz_yyyyyz[i] * pa_x[i];

        tr_z_xxxzz_yyyyzz[i] = 2.0 * tr_z_xzz_yyyyzz[i] * fe_0 + tr_z_xxzz_yyyyzz[i] * pa_x[i];

        tr_z_xxxzz_yyyzzz[i] = 2.0 * tr_z_xzz_yyyzzz[i] * fe_0 + tr_z_xxzz_yyyzzz[i] * pa_x[i];

        tr_z_xxxzz_yyzzzz[i] = 2.0 * tr_z_xzz_yyzzzz[i] * fe_0 + tr_z_xxzz_yyzzzz[i] * pa_x[i];

        tr_z_xxxzz_yzzzzz[i] = 2.0 * tr_z_xzz_yzzzzz[i] * fe_0 + tr_z_xxzz_yzzzzz[i] * pa_x[i];

        tr_z_xxxzz_zzzzzz[i] = 2.0 * tr_z_xzz_zzzzzz[i] * fe_0 + tr_z_xxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1344-1372 components of targeted buffer : HI

    auto tr_z_xxyyy_xxxxxx = pbuffer.data(idx_dip_hi + 1344);

    auto tr_z_xxyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1345);

    auto tr_z_xxyyy_xxxxxz = pbuffer.data(idx_dip_hi + 1346);

    auto tr_z_xxyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1347);

    auto tr_z_xxyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1348);

    auto tr_z_xxyyy_xxxxzz = pbuffer.data(idx_dip_hi + 1349);

    auto tr_z_xxyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1350);

    auto tr_z_xxyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1351);

    auto tr_z_xxyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1352);

    auto tr_z_xxyyy_xxxzzz = pbuffer.data(idx_dip_hi + 1353);

    auto tr_z_xxyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1354);

    auto tr_z_xxyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1355);

    auto tr_z_xxyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1356);

    auto tr_z_xxyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1357);

    auto tr_z_xxyyy_xxzzzz = pbuffer.data(idx_dip_hi + 1358);

    auto tr_z_xxyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1359);

    auto tr_z_xxyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1360);

    auto tr_z_xxyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1361);

    auto tr_z_xxyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1362);

    auto tr_z_xxyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1363);

    auto tr_z_xxyyy_xzzzzz = pbuffer.data(idx_dip_hi + 1364);

    auto tr_z_xxyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1365);

    auto tr_z_xxyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1366);

    auto tr_z_xxyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1367);

    auto tr_z_xxyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1368);

    auto tr_z_xxyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1369);

    auto tr_z_xxyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1370);

    auto tr_z_xxyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1371);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_z_xxy_xxxxxx,   \
                             tr_z_xxy_xxxxxz,   \
                             tr_z_xxy_xxxxzz,   \
                             tr_z_xxy_xxxzzz,   \
                             tr_z_xxy_xxzzzz,   \
                             tr_z_xxy_xzzzzz,   \
                             tr_z_xxyy_xxxxxx,  \
                             tr_z_xxyy_xxxxxz,  \
                             tr_z_xxyy_xxxxzz,  \
                             tr_z_xxyy_xxxzzz,  \
                             tr_z_xxyy_xxzzzz,  \
                             tr_z_xxyy_xzzzzz,  \
                             tr_z_xxyyy_xxxxxx, \
                             tr_z_xxyyy_xxxxxy, \
                             tr_z_xxyyy_xxxxxz, \
                             tr_z_xxyyy_xxxxyy, \
                             tr_z_xxyyy_xxxxyz, \
                             tr_z_xxyyy_xxxxzz, \
                             tr_z_xxyyy_xxxyyy, \
                             tr_z_xxyyy_xxxyyz, \
                             tr_z_xxyyy_xxxyzz, \
                             tr_z_xxyyy_xxxzzz, \
                             tr_z_xxyyy_xxyyyy, \
                             tr_z_xxyyy_xxyyyz, \
                             tr_z_xxyyy_xxyyzz, \
                             tr_z_xxyyy_xxyzzz, \
                             tr_z_xxyyy_xxzzzz, \
                             tr_z_xxyyy_xyyyyy, \
                             tr_z_xxyyy_xyyyyz, \
                             tr_z_xxyyy_xyyyzz, \
                             tr_z_xxyyy_xyyzzz, \
                             tr_z_xxyyy_xyzzzz, \
                             tr_z_xxyyy_xzzzzz, \
                             tr_z_xxyyy_yyyyyy, \
                             tr_z_xxyyy_yyyyyz, \
                             tr_z_xxyyy_yyyyzz, \
                             tr_z_xxyyy_yyyzzz, \
                             tr_z_xxyyy_yyzzzz, \
                             tr_z_xxyyy_yzzzzz, \
                             tr_z_xxyyy_zzzzzz, \
                             tr_z_xyyy_xxxxxy,  \
                             tr_z_xyyy_xxxxy,   \
                             tr_z_xyyy_xxxxyy,  \
                             tr_z_xyyy_xxxxyz,  \
                             tr_z_xyyy_xxxyy,   \
                             tr_z_xyyy_xxxyyy,  \
                             tr_z_xyyy_xxxyyz,  \
                             tr_z_xyyy_xxxyz,   \
                             tr_z_xyyy_xxxyzz,  \
                             tr_z_xyyy_xxyyy,   \
                             tr_z_xyyy_xxyyyy,  \
                             tr_z_xyyy_xxyyyz,  \
                             tr_z_xyyy_xxyyz,   \
                             tr_z_xyyy_xxyyzz,  \
                             tr_z_xyyy_xxyzz,   \
                             tr_z_xyyy_xxyzzz,  \
                             tr_z_xyyy_xyyyy,   \
                             tr_z_xyyy_xyyyyy,  \
                             tr_z_xyyy_xyyyyz,  \
                             tr_z_xyyy_xyyyz,   \
                             tr_z_xyyy_xyyyzz,  \
                             tr_z_xyyy_xyyzz,   \
                             tr_z_xyyy_xyyzzz,  \
                             tr_z_xyyy_xyzzz,   \
                             tr_z_xyyy_xyzzzz,  \
                             tr_z_xyyy_yyyyy,   \
                             tr_z_xyyy_yyyyyy,  \
                             tr_z_xyyy_yyyyyz,  \
                             tr_z_xyyy_yyyyz,   \
                             tr_z_xyyy_yyyyzz,  \
                             tr_z_xyyy_yyyzz,   \
                             tr_z_xyyy_yyyzzz,  \
                             tr_z_xyyy_yyzzz,   \
                             tr_z_xyyy_yyzzzz,  \
                             tr_z_xyyy_yzzzz,   \
                             tr_z_xyyy_yzzzzz,  \
                             tr_z_xyyy_zzzzzz,  \
                             tr_z_yyy_xxxxxy,   \
                             tr_z_yyy_xxxxyy,   \
                             tr_z_yyy_xxxxyz,   \
                             tr_z_yyy_xxxyyy,   \
                             tr_z_yyy_xxxyyz,   \
                             tr_z_yyy_xxxyzz,   \
                             tr_z_yyy_xxyyyy,   \
                             tr_z_yyy_xxyyyz,   \
                             tr_z_yyy_xxyyzz,   \
                             tr_z_yyy_xxyzzz,   \
                             tr_z_yyy_xyyyyy,   \
                             tr_z_yyy_xyyyyz,   \
                             tr_z_yyy_xyyyzz,   \
                             tr_z_yyy_xyyzzz,   \
                             tr_z_yyy_xyzzzz,   \
                             tr_z_yyy_yyyyyy,   \
                             tr_z_yyy_yyyyyz,   \
                             tr_z_yyy_yyyyzz,   \
                             tr_z_yyy_yyyzzz,   \
                             tr_z_yyy_yyzzzz,   \
                             tr_z_yyy_yzzzzz,   \
                             tr_z_yyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyy_xxxxxx[i] = 2.0 * tr_z_xxy_xxxxxx[i] * fe_0 + tr_z_xxyy_xxxxxx[i] * pa_y[i];

        tr_z_xxyyy_xxxxxy[i] = tr_z_yyy_xxxxxy[i] * fe_0 + 5.0 * tr_z_xyyy_xxxxy[i] * fe_0 + tr_z_xyyy_xxxxxy[i] * pa_x[i];

        tr_z_xxyyy_xxxxxz[i] = 2.0 * tr_z_xxy_xxxxxz[i] * fe_0 + tr_z_xxyy_xxxxxz[i] * pa_y[i];

        tr_z_xxyyy_xxxxyy[i] = tr_z_yyy_xxxxyy[i] * fe_0 + 4.0 * tr_z_xyyy_xxxyy[i] * fe_0 + tr_z_xyyy_xxxxyy[i] * pa_x[i];

        tr_z_xxyyy_xxxxyz[i] = tr_z_yyy_xxxxyz[i] * fe_0 + 4.0 * tr_z_xyyy_xxxyz[i] * fe_0 + tr_z_xyyy_xxxxyz[i] * pa_x[i];

        tr_z_xxyyy_xxxxzz[i] = 2.0 * tr_z_xxy_xxxxzz[i] * fe_0 + tr_z_xxyy_xxxxzz[i] * pa_y[i];

        tr_z_xxyyy_xxxyyy[i] = tr_z_yyy_xxxyyy[i] * fe_0 + 3.0 * tr_z_xyyy_xxyyy[i] * fe_0 + tr_z_xyyy_xxxyyy[i] * pa_x[i];

        tr_z_xxyyy_xxxyyz[i] = tr_z_yyy_xxxyyz[i] * fe_0 + 3.0 * tr_z_xyyy_xxyyz[i] * fe_0 + tr_z_xyyy_xxxyyz[i] * pa_x[i];

        tr_z_xxyyy_xxxyzz[i] = tr_z_yyy_xxxyzz[i] * fe_0 + 3.0 * tr_z_xyyy_xxyzz[i] * fe_0 + tr_z_xyyy_xxxyzz[i] * pa_x[i];

        tr_z_xxyyy_xxxzzz[i] = 2.0 * tr_z_xxy_xxxzzz[i] * fe_0 + tr_z_xxyy_xxxzzz[i] * pa_y[i];

        tr_z_xxyyy_xxyyyy[i] = tr_z_yyy_xxyyyy[i] * fe_0 + 2.0 * tr_z_xyyy_xyyyy[i] * fe_0 + tr_z_xyyy_xxyyyy[i] * pa_x[i];

        tr_z_xxyyy_xxyyyz[i] = tr_z_yyy_xxyyyz[i] * fe_0 + 2.0 * tr_z_xyyy_xyyyz[i] * fe_0 + tr_z_xyyy_xxyyyz[i] * pa_x[i];

        tr_z_xxyyy_xxyyzz[i] = tr_z_yyy_xxyyzz[i] * fe_0 + 2.0 * tr_z_xyyy_xyyzz[i] * fe_0 + tr_z_xyyy_xxyyzz[i] * pa_x[i];

        tr_z_xxyyy_xxyzzz[i] = tr_z_yyy_xxyzzz[i] * fe_0 + 2.0 * tr_z_xyyy_xyzzz[i] * fe_0 + tr_z_xyyy_xxyzzz[i] * pa_x[i];

        tr_z_xxyyy_xxzzzz[i] = 2.0 * tr_z_xxy_xxzzzz[i] * fe_0 + tr_z_xxyy_xxzzzz[i] * pa_y[i];

        tr_z_xxyyy_xyyyyy[i] = tr_z_yyy_xyyyyy[i] * fe_0 + tr_z_xyyy_yyyyy[i] * fe_0 + tr_z_xyyy_xyyyyy[i] * pa_x[i];

        tr_z_xxyyy_xyyyyz[i] = tr_z_yyy_xyyyyz[i] * fe_0 + tr_z_xyyy_yyyyz[i] * fe_0 + tr_z_xyyy_xyyyyz[i] * pa_x[i];

        tr_z_xxyyy_xyyyzz[i] = tr_z_yyy_xyyyzz[i] * fe_0 + tr_z_xyyy_yyyzz[i] * fe_0 + tr_z_xyyy_xyyyzz[i] * pa_x[i];

        tr_z_xxyyy_xyyzzz[i] = tr_z_yyy_xyyzzz[i] * fe_0 + tr_z_xyyy_yyzzz[i] * fe_0 + tr_z_xyyy_xyyzzz[i] * pa_x[i];

        tr_z_xxyyy_xyzzzz[i] = tr_z_yyy_xyzzzz[i] * fe_0 + tr_z_xyyy_yzzzz[i] * fe_0 + tr_z_xyyy_xyzzzz[i] * pa_x[i];

        tr_z_xxyyy_xzzzzz[i] = 2.0 * tr_z_xxy_xzzzzz[i] * fe_0 + tr_z_xxyy_xzzzzz[i] * pa_y[i];

        tr_z_xxyyy_yyyyyy[i] = tr_z_yyy_yyyyyy[i] * fe_0 + tr_z_xyyy_yyyyyy[i] * pa_x[i];

        tr_z_xxyyy_yyyyyz[i] = tr_z_yyy_yyyyyz[i] * fe_0 + tr_z_xyyy_yyyyyz[i] * pa_x[i];

        tr_z_xxyyy_yyyyzz[i] = tr_z_yyy_yyyyzz[i] * fe_0 + tr_z_xyyy_yyyyzz[i] * pa_x[i];

        tr_z_xxyyy_yyyzzz[i] = tr_z_yyy_yyyzzz[i] * fe_0 + tr_z_xyyy_yyyzzz[i] * pa_x[i];

        tr_z_xxyyy_yyzzzz[i] = tr_z_yyy_yyzzzz[i] * fe_0 + tr_z_xyyy_yyzzzz[i] * pa_x[i];

        tr_z_xxyyy_yzzzzz[i] = tr_z_yyy_yzzzzz[i] * fe_0 + tr_z_xyyy_yzzzzz[i] * pa_x[i];

        tr_z_xxyyy_zzzzzz[i] = tr_z_yyy_zzzzzz[i] * fe_0 + tr_z_xyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1372-1400 components of targeted buffer : HI

    auto tr_z_xxyyz_xxxxxx = pbuffer.data(idx_dip_hi + 1372);

    auto tr_z_xxyyz_xxxxxy = pbuffer.data(idx_dip_hi + 1373);

    auto tr_z_xxyyz_xxxxxz = pbuffer.data(idx_dip_hi + 1374);

    auto tr_z_xxyyz_xxxxyy = pbuffer.data(idx_dip_hi + 1375);

    auto tr_z_xxyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1376);

    auto tr_z_xxyyz_xxxxzz = pbuffer.data(idx_dip_hi + 1377);

    auto tr_z_xxyyz_xxxyyy = pbuffer.data(idx_dip_hi + 1378);

    auto tr_z_xxyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1379);

    auto tr_z_xxyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1380);

    auto tr_z_xxyyz_xxxzzz = pbuffer.data(idx_dip_hi + 1381);

    auto tr_z_xxyyz_xxyyyy = pbuffer.data(idx_dip_hi + 1382);

    auto tr_z_xxyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1383);

    auto tr_z_xxyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1384);

    auto tr_z_xxyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1385);

    auto tr_z_xxyyz_xxzzzz = pbuffer.data(idx_dip_hi + 1386);

    auto tr_z_xxyyz_xyyyyy = pbuffer.data(idx_dip_hi + 1387);

    auto tr_z_xxyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1388);

    auto tr_z_xxyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1389);

    auto tr_z_xxyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1390);

    auto tr_z_xxyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1391);

    auto tr_z_xxyyz_xzzzzz = pbuffer.data(idx_dip_hi + 1392);

    auto tr_z_xxyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1393);

    auto tr_z_xxyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1394);

    auto tr_z_xxyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1395);

    auto tr_z_xxyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1396);

    auto tr_z_xxyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1397);

    auto tr_z_xxyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1398);

    auto tr_z_xxyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1399);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             tr_z_xxyy_xxxxxy,  \
                             tr_z_xxyy_xxxxyy,  \
                             tr_z_xxyy_xxxyyy,  \
                             tr_z_xxyy_xxyyyy,  \
                             tr_z_xxyy_xyyyyy,  \
                             tr_z_xxyyz_xxxxxx, \
                             tr_z_xxyyz_xxxxxy, \
                             tr_z_xxyyz_xxxxxz, \
                             tr_z_xxyyz_xxxxyy, \
                             tr_z_xxyyz_xxxxyz, \
                             tr_z_xxyyz_xxxxzz, \
                             tr_z_xxyyz_xxxyyy, \
                             tr_z_xxyyz_xxxyyz, \
                             tr_z_xxyyz_xxxyzz, \
                             tr_z_xxyyz_xxxzzz, \
                             tr_z_xxyyz_xxyyyy, \
                             tr_z_xxyyz_xxyyyz, \
                             tr_z_xxyyz_xxyyzz, \
                             tr_z_xxyyz_xxyzzz, \
                             tr_z_xxyyz_xxzzzz, \
                             tr_z_xxyyz_xyyyyy, \
                             tr_z_xxyyz_xyyyyz, \
                             tr_z_xxyyz_xyyyzz, \
                             tr_z_xxyyz_xyyzzz, \
                             tr_z_xxyyz_xyzzzz, \
                             tr_z_xxyyz_xzzzzz, \
                             tr_z_xxyyz_yyyyyy, \
                             tr_z_xxyyz_yyyyyz, \
                             tr_z_xxyyz_yyyyzz, \
                             tr_z_xxyyz_yyyzzz, \
                             tr_z_xxyyz_yyzzzz, \
                             tr_z_xxyyz_yzzzzz, \
                             tr_z_xxyyz_zzzzzz, \
                             tr_z_xxyz_xxxxxx,  \
                             tr_z_xxyz_xxxxxz,  \
                             tr_z_xxyz_xxxxzz,  \
                             tr_z_xxyz_xxxzzz,  \
                             tr_z_xxyz_xxzzzz,  \
                             tr_z_xxyz_xzzzzz,  \
                             tr_z_xxz_xxxxxx,   \
                             tr_z_xxz_xxxxxz,   \
                             tr_z_xxz_xxxxzz,   \
                             tr_z_xxz_xxxzzz,   \
                             tr_z_xxz_xxzzzz,   \
                             tr_z_xxz_xzzzzz,   \
                             tr_z_xyyz_xxxxyz,  \
                             tr_z_xyyz_xxxyyz,  \
                             tr_z_xyyz_xxxyz,   \
                             tr_z_xyyz_xxxyzz,  \
                             tr_z_xyyz_xxyyyz,  \
                             tr_z_xyyz_xxyyz,   \
                             tr_z_xyyz_xxyyzz,  \
                             tr_z_xyyz_xxyzz,   \
                             tr_z_xyyz_xxyzzz,  \
                             tr_z_xyyz_xyyyyz,  \
                             tr_z_xyyz_xyyyz,   \
                             tr_z_xyyz_xyyyzz,  \
                             tr_z_xyyz_xyyzz,   \
                             tr_z_xyyz_xyyzzz,  \
                             tr_z_xyyz_xyzzz,   \
                             tr_z_xyyz_xyzzzz,  \
                             tr_z_xyyz_yyyyyy,  \
                             tr_z_xyyz_yyyyyz,  \
                             tr_z_xyyz_yyyyz,   \
                             tr_z_xyyz_yyyyzz,  \
                             tr_z_xyyz_yyyzz,   \
                             tr_z_xyyz_yyyzzz,  \
                             tr_z_xyyz_yyzzz,   \
                             tr_z_xyyz_yyzzzz,  \
                             tr_z_xyyz_yzzzz,   \
                             tr_z_xyyz_yzzzzz,  \
                             tr_z_xyyz_zzzzzz,  \
                             tr_z_yyz_xxxxyz,   \
                             tr_z_yyz_xxxyyz,   \
                             tr_z_yyz_xxxyzz,   \
                             tr_z_yyz_xxyyyz,   \
                             tr_z_yyz_xxyyzz,   \
                             tr_z_yyz_xxyzzz,   \
                             tr_z_yyz_xyyyyz,   \
                             tr_z_yyz_xyyyzz,   \
                             tr_z_yyz_xyyzzz,   \
                             tr_z_yyz_xyzzzz,   \
                             tr_z_yyz_yyyyyy,   \
                             tr_z_yyz_yyyyyz,   \
                             tr_z_yyz_yyyyzz,   \
                             tr_z_yyz_yyyzzz,   \
                             tr_z_yyz_yyzzzz,   \
                             tr_z_yyz_yzzzzz,   \
                             tr_z_yyz_zzzzzz,   \
                             ts_xxyy_xxxxxy,    \
                             ts_xxyy_xxxxyy,    \
                             ts_xxyy_xxxyyy,    \
                             ts_xxyy_xxyyyy,    \
                             ts_xxyy_xyyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyz_xxxxxx[i] = tr_z_xxz_xxxxxx[i] * fe_0 + tr_z_xxyz_xxxxxx[i] * pa_y[i];

        tr_z_xxyyz_xxxxxy[i] = ts_xxyy_xxxxxy[i] * fe_0 + tr_z_xxyy_xxxxxy[i] * pa_z[i];

        tr_z_xxyyz_xxxxxz[i] = tr_z_xxz_xxxxxz[i] * fe_0 + tr_z_xxyz_xxxxxz[i] * pa_y[i];

        tr_z_xxyyz_xxxxyy[i] = ts_xxyy_xxxxyy[i] * fe_0 + tr_z_xxyy_xxxxyy[i] * pa_z[i];

        tr_z_xxyyz_xxxxyz[i] = tr_z_yyz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xyyz_xxxyz[i] * fe_0 + tr_z_xyyz_xxxxyz[i] * pa_x[i];

        tr_z_xxyyz_xxxxzz[i] = tr_z_xxz_xxxxzz[i] * fe_0 + tr_z_xxyz_xxxxzz[i] * pa_y[i];

        tr_z_xxyyz_xxxyyy[i] = ts_xxyy_xxxyyy[i] * fe_0 + tr_z_xxyy_xxxyyy[i] * pa_z[i];

        tr_z_xxyyz_xxxyyz[i] = tr_z_yyz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xyyz_xxyyz[i] * fe_0 + tr_z_xyyz_xxxyyz[i] * pa_x[i];

        tr_z_xxyyz_xxxyzz[i] = tr_z_yyz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xyyz_xxyzz[i] * fe_0 + tr_z_xyyz_xxxyzz[i] * pa_x[i];

        tr_z_xxyyz_xxxzzz[i] = tr_z_xxz_xxxzzz[i] * fe_0 + tr_z_xxyz_xxxzzz[i] * pa_y[i];

        tr_z_xxyyz_xxyyyy[i] = ts_xxyy_xxyyyy[i] * fe_0 + tr_z_xxyy_xxyyyy[i] * pa_z[i];

        tr_z_xxyyz_xxyyyz[i] = tr_z_yyz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xyyz_xyyyz[i] * fe_0 + tr_z_xyyz_xxyyyz[i] * pa_x[i];

        tr_z_xxyyz_xxyyzz[i] = tr_z_yyz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xyyz_xyyzz[i] * fe_0 + tr_z_xyyz_xxyyzz[i] * pa_x[i];

        tr_z_xxyyz_xxyzzz[i] = tr_z_yyz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xyyz_xyzzz[i] * fe_0 + tr_z_xyyz_xxyzzz[i] * pa_x[i];

        tr_z_xxyyz_xxzzzz[i] = tr_z_xxz_xxzzzz[i] * fe_0 + tr_z_xxyz_xxzzzz[i] * pa_y[i];

        tr_z_xxyyz_xyyyyy[i] = ts_xxyy_xyyyyy[i] * fe_0 + tr_z_xxyy_xyyyyy[i] * pa_z[i];

        tr_z_xxyyz_xyyyyz[i] = tr_z_yyz_xyyyyz[i] * fe_0 + tr_z_xyyz_yyyyz[i] * fe_0 + tr_z_xyyz_xyyyyz[i] * pa_x[i];

        tr_z_xxyyz_xyyyzz[i] = tr_z_yyz_xyyyzz[i] * fe_0 + tr_z_xyyz_yyyzz[i] * fe_0 + tr_z_xyyz_xyyyzz[i] * pa_x[i];

        tr_z_xxyyz_xyyzzz[i] = tr_z_yyz_xyyzzz[i] * fe_0 + tr_z_xyyz_yyzzz[i] * fe_0 + tr_z_xyyz_xyyzzz[i] * pa_x[i];

        tr_z_xxyyz_xyzzzz[i] = tr_z_yyz_xyzzzz[i] * fe_0 + tr_z_xyyz_yzzzz[i] * fe_0 + tr_z_xyyz_xyzzzz[i] * pa_x[i];

        tr_z_xxyyz_xzzzzz[i] = tr_z_xxz_xzzzzz[i] * fe_0 + tr_z_xxyz_xzzzzz[i] * pa_y[i];

        tr_z_xxyyz_yyyyyy[i] = tr_z_yyz_yyyyyy[i] * fe_0 + tr_z_xyyz_yyyyyy[i] * pa_x[i];

        tr_z_xxyyz_yyyyyz[i] = tr_z_yyz_yyyyyz[i] * fe_0 + tr_z_xyyz_yyyyyz[i] * pa_x[i];

        tr_z_xxyyz_yyyyzz[i] = tr_z_yyz_yyyyzz[i] * fe_0 + tr_z_xyyz_yyyyzz[i] * pa_x[i];

        tr_z_xxyyz_yyyzzz[i] = tr_z_yyz_yyyzzz[i] * fe_0 + tr_z_xyyz_yyyzzz[i] * pa_x[i];

        tr_z_xxyyz_yyzzzz[i] = tr_z_yyz_yyzzzz[i] * fe_0 + tr_z_xyyz_yyzzzz[i] * pa_x[i];

        tr_z_xxyyz_yzzzzz[i] = tr_z_yyz_yzzzzz[i] * fe_0 + tr_z_xyyz_yzzzzz[i] * pa_x[i];

        tr_z_xxyyz_zzzzzz[i] = tr_z_yyz_zzzzzz[i] * fe_0 + tr_z_xyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1400-1428 components of targeted buffer : HI

    auto tr_z_xxyzz_xxxxxx = pbuffer.data(idx_dip_hi + 1400);

    auto tr_z_xxyzz_xxxxxy = pbuffer.data(idx_dip_hi + 1401);

    auto tr_z_xxyzz_xxxxxz = pbuffer.data(idx_dip_hi + 1402);

    auto tr_z_xxyzz_xxxxyy = pbuffer.data(idx_dip_hi + 1403);

    auto tr_z_xxyzz_xxxxyz = pbuffer.data(idx_dip_hi + 1404);

    auto tr_z_xxyzz_xxxxzz = pbuffer.data(idx_dip_hi + 1405);

    auto tr_z_xxyzz_xxxyyy = pbuffer.data(idx_dip_hi + 1406);

    auto tr_z_xxyzz_xxxyyz = pbuffer.data(idx_dip_hi + 1407);

    auto tr_z_xxyzz_xxxyzz = pbuffer.data(idx_dip_hi + 1408);

    auto tr_z_xxyzz_xxxzzz = pbuffer.data(idx_dip_hi + 1409);

    auto tr_z_xxyzz_xxyyyy = pbuffer.data(idx_dip_hi + 1410);

    auto tr_z_xxyzz_xxyyyz = pbuffer.data(idx_dip_hi + 1411);

    auto tr_z_xxyzz_xxyyzz = pbuffer.data(idx_dip_hi + 1412);

    auto tr_z_xxyzz_xxyzzz = pbuffer.data(idx_dip_hi + 1413);

    auto tr_z_xxyzz_xxzzzz = pbuffer.data(idx_dip_hi + 1414);

    auto tr_z_xxyzz_xyyyyy = pbuffer.data(idx_dip_hi + 1415);

    auto tr_z_xxyzz_xyyyyz = pbuffer.data(idx_dip_hi + 1416);

    auto tr_z_xxyzz_xyyyzz = pbuffer.data(idx_dip_hi + 1417);

    auto tr_z_xxyzz_xyyzzz = pbuffer.data(idx_dip_hi + 1418);

    auto tr_z_xxyzz_xyzzzz = pbuffer.data(idx_dip_hi + 1419);

    auto tr_z_xxyzz_xzzzzz = pbuffer.data(idx_dip_hi + 1420);

    auto tr_z_xxyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1421);

    auto tr_z_xxyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1422);

    auto tr_z_xxyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1423);

    auto tr_z_xxyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1424);

    auto tr_z_xxyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1425);

    auto tr_z_xxyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1426);

    auto tr_z_xxyzz_zzzzzz = pbuffer.data(idx_dip_hi + 1427);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_z_xxyzz_xxxxxx, \
                             tr_z_xxyzz_xxxxxy, \
                             tr_z_xxyzz_xxxxxz, \
                             tr_z_xxyzz_xxxxyy, \
                             tr_z_xxyzz_xxxxyz, \
                             tr_z_xxyzz_xxxxzz, \
                             tr_z_xxyzz_xxxyyy, \
                             tr_z_xxyzz_xxxyyz, \
                             tr_z_xxyzz_xxxyzz, \
                             tr_z_xxyzz_xxxzzz, \
                             tr_z_xxyzz_xxyyyy, \
                             tr_z_xxyzz_xxyyyz, \
                             tr_z_xxyzz_xxyyzz, \
                             tr_z_xxyzz_xxyzzz, \
                             tr_z_xxyzz_xxzzzz, \
                             tr_z_xxyzz_xyyyyy, \
                             tr_z_xxyzz_xyyyyz, \
                             tr_z_xxyzz_xyyyzz, \
                             tr_z_xxyzz_xyyzzz, \
                             tr_z_xxyzz_xyzzzz, \
                             tr_z_xxyzz_xzzzzz, \
                             tr_z_xxyzz_yyyyyy, \
                             tr_z_xxyzz_yyyyyz, \
                             tr_z_xxyzz_yyyyzz, \
                             tr_z_xxyzz_yyyzzz, \
                             tr_z_xxyzz_yyzzzz, \
                             tr_z_xxyzz_yzzzzz, \
                             tr_z_xxyzz_zzzzzz, \
                             tr_z_xxzz_xxxxx,   \
                             tr_z_xxzz_xxxxxx,  \
                             tr_z_xxzz_xxxxxy,  \
                             tr_z_xxzz_xxxxxz,  \
                             tr_z_xxzz_xxxxy,   \
                             tr_z_xxzz_xxxxyy,  \
                             tr_z_xxzz_xxxxyz,  \
                             tr_z_xxzz_xxxxz,   \
                             tr_z_xxzz_xxxxzz,  \
                             tr_z_xxzz_xxxyy,   \
                             tr_z_xxzz_xxxyyy,  \
                             tr_z_xxzz_xxxyyz,  \
                             tr_z_xxzz_xxxyz,   \
                             tr_z_xxzz_xxxyzz,  \
                             tr_z_xxzz_xxxzz,   \
                             tr_z_xxzz_xxxzzz,  \
                             tr_z_xxzz_xxyyy,   \
                             tr_z_xxzz_xxyyyy,  \
                             tr_z_xxzz_xxyyyz,  \
                             tr_z_xxzz_xxyyz,   \
                             tr_z_xxzz_xxyyzz,  \
                             tr_z_xxzz_xxyzz,   \
                             tr_z_xxzz_xxyzzz,  \
                             tr_z_xxzz_xxzzz,   \
                             tr_z_xxzz_xxzzzz,  \
                             tr_z_xxzz_xyyyy,   \
                             tr_z_xxzz_xyyyyy,  \
                             tr_z_xxzz_xyyyyz,  \
                             tr_z_xxzz_xyyyz,   \
                             tr_z_xxzz_xyyyzz,  \
                             tr_z_xxzz_xyyzz,   \
                             tr_z_xxzz_xyyzzz,  \
                             tr_z_xxzz_xyzzz,   \
                             tr_z_xxzz_xyzzzz,  \
                             tr_z_xxzz_xzzzz,   \
                             tr_z_xxzz_xzzzzz,  \
                             tr_z_xxzz_zzzzzz,  \
                             tr_z_xyzz_yyyyyy,  \
                             tr_z_xyzz_yyyyyz,  \
                             tr_z_xyzz_yyyyzz,  \
                             tr_z_xyzz_yyyzzz,  \
                             tr_z_xyzz_yyzzzz,  \
                             tr_z_xyzz_yzzzzz,  \
                             tr_z_yzz_yyyyyy,   \
                             tr_z_yzz_yyyyyz,   \
                             tr_z_yzz_yyyyzz,   \
                             tr_z_yzz_yyyzzz,   \
                             tr_z_yzz_yyzzzz,   \
                             tr_z_yzz_yzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzz_xxxxxx[i] = tr_z_xxzz_xxxxxx[i] * pa_y[i];

        tr_z_xxyzz_xxxxxy[i] = tr_z_xxzz_xxxxx[i] * fe_0 + tr_z_xxzz_xxxxxy[i] * pa_y[i];

        tr_z_xxyzz_xxxxxz[i] = tr_z_xxzz_xxxxxz[i] * pa_y[i];

        tr_z_xxyzz_xxxxyy[i] = 2.0 * tr_z_xxzz_xxxxy[i] * fe_0 + tr_z_xxzz_xxxxyy[i] * pa_y[i];

        tr_z_xxyzz_xxxxyz[i] = tr_z_xxzz_xxxxz[i] * fe_0 + tr_z_xxzz_xxxxyz[i] * pa_y[i];

        tr_z_xxyzz_xxxxzz[i] = tr_z_xxzz_xxxxzz[i] * pa_y[i];

        tr_z_xxyzz_xxxyyy[i] = 3.0 * tr_z_xxzz_xxxyy[i] * fe_0 + tr_z_xxzz_xxxyyy[i] * pa_y[i];

        tr_z_xxyzz_xxxyyz[i] = 2.0 * tr_z_xxzz_xxxyz[i] * fe_0 + tr_z_xxzz_xxxyyz[i] * pa_y[i];

        tr_z_xxyzz_xxxyzz[i] = tr_z_xxzz_xxxzz[i] * fe_0 + tr_z_xxzz_xxxyzz[i] * pa_y[i];

        tr_z_xxyzz_xxxzzz[i] = tr_z_xxzz_xxxzzz[i] * pa_y[i];

        tr_z_xxyzz_xxyyyy[i] = 4.0 * tr_z_xxzz_xxyyy[i] * fe_0 + tr_z_xxzz_xxyyyy[i] * pa_y[i];

        tr_z_xxyzz_xxyyyz[i] = 3.0 * tr_z_xxzz_xxyyz[i] * fe_0 + tr_z_xxzz_xxyyyz[i] * pa_y[i];

        tr_z_xxyzz_xxyyzz[i] = 2.0 * tr_z_xxzz_xxyzz[i] * fe_0 + tr_z_xxzz_xxyyzz[i] * pa_y[i];

        tr_z_xxyzz_xxyzzz[i] = tr_z_xxzz_xxzzz[i] * fe_0 + tr_z_xxzz_xxyzzz[i] * pa_y[i];

        tr_z_xxyzz_xxzzzz[i] = tr_z_xxzz_xxzzzz[i] * pa_y[i];

        tr_z_xxyzz_xyyyyy[i] = 5.0 * tr_z_xxzz_xyyyy[i] * fe_0 + tr_z_xxzz_xyyyyy[i] * pa_y[i];

        tr_z_xxyzz_xyyyyz[i] = 4.0 * tr_z_xxzz_xyyyz[i] * fe_0 + tr_z_xxzz_xyyyyz[i] * pa_y[i];

        tr_z_xxyzz_xyyyzz[i] = 3.0 * tr_z_xxzz_xyyzz[i] * fe_0 + tr_z_xxzz_xyyyzz[i] * pa_y[i];

        tr_z_xxyzz_xyyzzz[i] = 2.0 * tr_z_xxzz_xyzzz[i] * fe_0 + tr_z_xxzz_xyyzzz[i] * pa_y[i];

        tr_z_xxyzz_xyzzzz[i] = tr_z_xxzz_xzzzz[i] * fe_0 + tr_z_xxzz_xyzzzz[i] * pa_y[i];

        tr_z_xxyzz_xzzzzz[i] = tr_z_xxzz_xzzzzz[i] * pa_y[i];

        tr_z_xxyzz_yyyyyy[i] = tr_z_yzz_yyyyyy[i] * fe_0 + tr_z_xyzz_yyyyyy[i] * pa_x[i];

        tr_z_xxyzz_yyyyyz[i] = tr_z_yzz_yyyyyz[i] * fe_0 + tr_z_xyzz_yyyyyz[i] * pa_x[i];

        tr_z_xxyzz_yyyyzz[i] = tr_z_yzz_yyyyzz[i] * fe_0 + tr_z_xyzz_yyyyzz[i] * pa_x[i];

        tr_z_xxyzz_yyyzzz[i] = tr_z_yzz_yyyzzz[i] * fe_0 + tr_z_xyzz_yyyzzz[i] * pa_x[i];

        tr_z_xxyzz_yyzzzz[i] = tr_z_yzz_yyzzzz[i] * fe_0 + tr_z_xyzz_yyzzzz[i] * pa_x[i];

        tr_z_xxyzz_yzzzzz[i] = tr_z_yzz_yzzzzz[i] * fe_0 + tr_z_xyzz_yzzzzz[i] * pa_x[i];

        tr_z_xxyzz_zzzzzz[i] = tr_z_xxzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1428-1456 components of targeted buffer : HI

    auto tr_z_xxzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1428);

    auto tr_z_xxzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1429);

    auto tr_z_xxzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1430);

    auto tr_z_xxzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1431);

    auto tr_z_xxzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1432);

    auto tr_z_xxzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1433);

    auto tr_z_xxzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1434);

    auto tr_z_xxzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1435);

    auto tr_z_xxzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1436);

    auto tr_z_xxzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1437);

    auto tr_z_xxzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1438);

    auto tr_z_xxzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1439);

    auto tr_z_xxzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1440);

    auto tr_z_xxzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1441);

    auto tr_z_xxzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1442);

    auto tr_z_xxzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1443);

    auto tr_z_xxzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1444);

    auto tr_z_xxzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1445);

    auto tr_z_xxzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1446);

    auto tr_z_xxzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1447);

    auto tr_z_xxzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1448);

    auto tr_z_xxzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1449);

    auto tr_z_xxzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1450);

    auto tr_z_xxzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1451);

    auto tr_z_xxzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1452);

    auto tr_z_xxzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1453);

    auto tr_z_xxzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1454);

    auto tr_z_xxzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1455);

#pragma omp simd aligned(pa_x,                  \
                             tr_z_xxzzz_xxxxxx, \
                             tr_z_xxzzz_xxxxxy, \
                             tr_z_xxzzz_xxxxxz, \
                             tr_z_xxzzz_xxxxyy, \
                             tr_z_xxzzz_xxxxyz, \
                             tr_z_xxzzz_xxxxzz, \
                             tr_z_xxzzz_xxxyyy, \
                             tr_z_xxzzz_xxxyyz, \
                             tr_z_xxzzz_xxxyzz, \
                             tr_z_xxzzz_xxxzzz, \
                             tr_z_xxzzz_xxyyyy, \
                             tr_z_xxzzz_xxyyyz, \
                             tr_z_xxzzz_xxyyzz, \
                             tr_z_xxzzz_xxyzzz, \
                             tr_z_xxzzz_xxzzzz, \
                             tr_z_xxzzz_xyyyyy, \
                             tr_z_xxzzz_xyyyyz, \
                             tr_z_xxzzz_xyyyzz, \
                             tr_z_xxzzz_xyyzzz, \
                             tr_z_xxzzz_xyzzzz, \
                             tr_z_xxzzz_xzzzzz, \
                             tr_z_xxzzz_yyyyyy, \
                             tr_z_xxzzz_yyyyyz, \
                             tr_z_xxzzz_yyyyzz, \
                             tr_z_xxzzz_yyyzzz, \
                             tr_z_xxzzz_yyzzzz, \
                             tr_z_xxzzz_yzzzzz, \
                             tr_z_xxzzz_zzzzzz, \
                             tr_z_xzzz_xxxxx,   \
                             tr_z_xzzz_xxxxxx,  \
                             tr_z_xzzz_xxxxxy,  \
                             tr_z_xzzz_xxxxxz,  \
                             tr_z_xzzz_xxxxy,   \
                             tr_z_xzzz_xxxxyy,  \
                             tr_z_xzzz_xxxxyz,  \
                             tr_z_xzzz_xxxxz,   \
                             tr_z_xzzz_xxxxzz,  \
                             tr_z_xzzz_xxxyy,   \
                             tr_z_xzzz_xxxyyy,  \
                             tr_z_xzzz_xxxyyz,  \
                             tr_z_xzzz_xxxyz,   \
                             tr_z_xzzz_xxxyzz,  \
                             tr_z_xzzz_xxxzz,   \
                             tr_z_xzzz_xxxzzz,  \
                             tr_z_xzzz_xxyyy,   \
                             tr_z_xzzz_xxyyyy,  \
                             tr_z_xzzz_xxyyyz,  \
                             tr_z_xzzz_xxyyz,   \
                             tr_z_xzzz_xxyyzz,  \
                             tr_z_xzzz_xxyzz,   \
                             tr_z_xzzz_xxyzzz,  \
                             tr_z_xzzz_xxzzz,   \
                             tr_z_xzzz_xxzzzz,  \
                             tr_z_xzzz_xyyyy,   \
                             tr_z_xzzz_xyyyyy,  \
                             tr_z_xzzz_xyyyyz,  \
                             tr_z_xzzz_xyyyz,   \
                             tr_z_xzzz_xyyyzz,  \
                             tr_z_xzzz_xyyzz,   \
                             tr_z_xzzz_xyyzzz,  \
                             tr_z_xzzz_xyzzz,   \
                             tr_z_xzzz_xyzzzz,  \
                             tr_z_xzzz_xzzzz,   \
                             tr_z_xzzz_xzzzzz,  \
                             tr_z_xzzz_yyyyy,   \
                             tr_z_xzzz_yyyyyy,  \
                             tr_z_xzzz_yyyyyz,  \
                             tr_z_xzzz_yyyyz,   \
                             tr_z_xzzz_yyyyzz,  \
                             tr_z_xzzz_yyyzz,   \
                             tr_z_xzzz_yyyzzz,  \
                             tr_z_xzzz_yyzzz,   \
                             tr_z_xzzz_yyzzzz,  \
                             tr_z_xzzz_yzzzz,   \
                             tr_z_xzzz_yzzzzz,  \
                             tr_z_xzzz_zzzzz,   \
                             tr_z_xzzz_zzzzzz,  \
                             tr_z_zzz_xxxxxx,   \
                             tr_z_zzz_xxxxxy,   \
                             tr_z_zzz_xxxxxz,   \
                             tr_z_zzz_xxxxyy,   \
                             tr_z_zzz_xxxxyz,   \
                             tr_z_zzz_xxxxzz,   \
                             tr_z_zzz_xxxyyy,   \
                             tr_z_zzz_xxxyyz,   \
                             tr_z_zzz_xxxyzz,   \
                             tr_z_zzz_xxxzzz,   \
                             tr_z_zzz_xxyyyy,   \
                             tr_z_zzz_xxyyyz,   \
                             tr_z_zzz_xxyyzz,   \
                             tr_z_zzz_xxyzzz,   \
                             tr_z_zzz_xxzzzz,   \
                             tr_z_zzz_xyyyyy,   \
                             tr_z_zzz_xyyyyz,   \
                             tr_z_zzz_xyyyzz,   \
                             tr_z_zzz_xyyzzz,   \
                             tr_z_zzz_xyzzzz,   \
                             tr_z_zzz_xzzzzz,   \
                             tr_z_zzz_yyyyyy,   \
                             tr_z_zzz_yyyyyz,   \
                             tr_z_zzz_yyyyzz,   \
                             tr_z_zzz_yyyzzz,   \
                             tr_z_zzz_yyzzzz,   \
                             tr_z_zzz_yzzzzz,   \
                             tr_z_zzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzz_xxxxxx[i] = tr_z_zzz_xxxxxx[i] * fe_0 + 6.0 * tr_z_xzzz_xxxxx[i] * fe_0 + tr_z_xzzz_xxxxxx[i] * pa_x[i];

        tr_z_xxzzz_xxxxxy[i] = tr_z_zzz_xxxxxy[i] * fe_0 + 5.0 * tr_z_xzzz_xxxxy[i] * fe_0 + tr_z_xzzz_xxxxxy[i] * pa_x[i];

        tr_z_xxzzz_xxxxxz[i] = tr_z_zzz_xxxxxz[i] * fe_0 + 5.0 * tr_z_xzzz_xxxxz[i] * fe_0 + tr_z_xzzz_xxxxxz[i] * pa_x[i];

        tr_z_xxzzz_xxxxyy[i] = tr_z_zzz_xxxxyy[i] * fe_0 + 4.0 * tr_z_xzzz_xxxyy[i] * fe_0 + tr_z_xzzz_xxxxyy[i] * pa_x[i];

        tr_z_xxzzz_xxxxyz[i] = tr_z_zzz_xxxxyz[i] * fe_0 + 4.0 * tr_z_xzzz_xxxyz[i] * fe_0 + tr_z_xzzz_xxxxyz[i] * pa_x[i];

        tr_z_xxzzz_xxxxzz[i] = tr_z_zzz_xxxxzz[i] * fe_0 + 4.0 * tr_z_xzzz_xxxzz[i] * fe_0 + tr_z_xzzz_xxxxzz[i] * pa_x[i];

        tr_z_xxzzz_xxxyyy[i] = tr_z_zzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_xzzz_xxyyy[i] * fe_0 + tr_z_xzzz_xxxyyy[i] * pa_x[i];

        tr_z_xxzzz_xxxyyz[i] = tr_z_zzz_xxxyyz[i] * fe_0 + 3.0 * tr_z_xzzz_xxyyz[i] * fe_0 + tr_z_xzzz_xxxyyz[i] * pa_x[i];

        tr_z_xxzzz_xxxyzz[i] = tr_z_zzz_xxxyzz[i] * fe_0 + 3.0 * tr_z_xzzz_xxyzz[i] * fe_0 + tr_z_xzzz_xxxyzz[i] * pa_x[i];

        tr_z_xxzzz_xxxzzz[i] = tr_z_zzz_xxxzzz[i] * fe_0 + 3.0 * tr_z_xzzz_xxzzz[i] * fe_0 + tr_z_xzzz_xxxzzz[i] * pa_x[i];

        tr_z_xxzzz_xxyyyy[i] = tr_z_zzz_xxyyyy[i] * fe_0 + 2.0 * tr_z_xzzz_xyyyy[i] * fe_0 + tr_z_xzzz_xxyyyy[i] * pa_x[i];

        tr_z_xxzzz_xxyyyz[i] = tr_z_zzz_xxyyyz[i] * fe_0 + 2.0 * tr_z_xzzz_xyyyz[i] * fe_0 + tr_z_xzzz_xxyyyz[i] * pa_x[i];

        tr_z_xxzzz_xxyyzz[i] = tr_z_zzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_xzzz_xyyzz[i] * fe_0 + tr_z_xzzz_xxyyzz[i] * pa_x[i];

        tr_z_xxzzz_xxyzzz[i] = tr_z_zzz_xxyzzz[i] * fe_0 + 2.0 * tr_z_xzzz_xyzzz[i] * fe_0 + tr_z_xzzz_xxyzzz[i] * pa_x[i];

        tr_z_xxzzz_xxzzzz[i] = tr_z_zzz_xxzzzz[i] * fe_0 + 2.0 * tr_z_xzzz_xzzzz[i] * fe_0 + tr_z_xzzz_xxzzzz[i] * pa_x[i];

        tr_z_xxzzz_xyyyyy[i] = tr_z_zzz_xyyyyy[i] * fe_0 + tr_z_xzzz_yyyyy[i] * fe_0 + tr_z_xzzz_xyyyyy[i] * pa_x[i];

        tr_z_xxzzz_xyyyyz[i] = tr_z_zzz_xyyyyz[i] * fe_0 + tr_z_xzzz_yyyyz[i] * fe_0 + tr_z_xzzz_xyyyyz[i] * pa_x[i];

        tr_z_xxzzz_xyyyzz[i] = tr_z_zzz_xyyyzz[i] * fe_0 + tr_z_xzzz_yyyzz[i] * fe_0 + tr_z_xzzz_xyyyzz[i] * pa_x[i];

        tr_z_xxzzz_xyyzzz[i] = tr_z_zzz_xyyzzz[i] * fe_0 + tr_z_xzzz_yyzzz[i] * fe_0 + tr_z_xzzz_xyyzzz[i] * pa_x[i];

        tr_z_xxzzz_xyzzzz[i] = tr_z_zzz_xyzzzz[i] * fe_0 + tr_z_xzzz_yzzzz[i] * fe_0 + tr_z_xzzz_xyzzzz[i] * pa_x[i];

        tr_z_xxzzz_xzzzzz[i] = tr_z_zzz_xzzzzz[i] * fe_0 + tr_z_xzzz_zzzzz[i] * fe_0 + tr_z_xzzz_xzzzzz[i] * pa_x[i];

        tr_z_xxzzz_yyyyyy[i] = tr_z_zzz_yyyyyy[i] * fe_0 + tr_z_xzzz_yyyyyy[i] * pa_x[i];

        tr_z_xxzzz_yyyyyz[i] = tr_z_zzz_yyyyyz[i] * fe_0 + tr_z_xzzz_yyyyyz[i] * pa_x[i];

        tr_z_xxzzz_yyyyzz[i] = tr_z_zzz_yyyyzz[i] * fe_0 + tr_z_xzzz_yyyyzz[i] * pa_x[i];

        tr_z_xxzzz_yyyzzz[i] = tr_z_zzz_yyyzzz[i] * fe_0 + tr_z_xzzz_yyyzzz[i] * pa_x[i];

        tr_z_xxzzz_yyzzzz[i] = tr_z_zzz_yyzzzz[i] * fe_0 + tr_z_xzzz_yyzzzz[i] * pa_x[i];

        tr_z_xxzzz_yzzzzz[i] = tr_z_zzz_yzzzzz[i] * fe_0 + tr_z_xzzz_yzzzzz[i] * pa_x[i];

        tr_z_xxzzz_zzzzzz[i] = tr_z_zzz_zzzzzz[i] * fe_0 + tr_z_xzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1456-1484 components of targeted buffer : HI

    auto tr_z_xyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 1456);

    auto tr_z_xyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1457);

    auto tr_z_xyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 1458);

    auto tr_z_xyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1459);

    auto tr_z_xyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1460);

    auto tr_z_xyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 1461);

    auto tr_z_xyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1462);

    auto tr_z_xyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1463);

    auto tr_z_xyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1464);

    auto tr_z_xyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 1465);

    auto tr_z_xyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1466);

    auto tr_z_xyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1467);

    auto tr_z_xyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1468);

    auto tr_z_xyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1469);

    auto tr_z_xyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 1470);

    auto tr_z_xyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1471);

    auto tr_z_xyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1472);

    auto tr_z_xyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1473);

    auto tr_z_xyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1474);

    auto tr_z_xyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1475);

    auto tr_z_xyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 1476);

    auto tr_z_xyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1477);

    auto tr_z_xyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1478);

    auto tr_z_xyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1479);

    auto tr_z_xyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1480);

    auto tr_z_xyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1481);

    auto tr_z_xyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1482);

    auto tr_z_xyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1483);

#pragma omp simd aligned(pa_x,                  \
                             tr_z_xyyyy_xxxxxx, \
                             tr_z_xyyyy_xxxxxy, \
                             tr_z_xyyyy_xxxxxz, \
                             tr_z_xyyyy_xxxxyy, \
                             tr_z_xyyyy_xxxxyz, \
                             tr_z_xyyyy_xxxxzz, \
                             tr_z_xyyyy_xxxyyy, \
                             tr_z_xyyyy_xxxyyz, \
                             tr_z_xyyyy_xxxyzz, \
                             tr_z_xyyyy_xxxzzz, \
                             tr_z_xyyyy_xxyyyy, \
                             tr_z_xyyyy_xxyyyz, \
                             tr_z_xyyyy_xxyyzz, \
                             tr_z_xyyyy_xxyzzz, \
                             tr_z_xyyyy_xxzzzz, \
                             tr_z_xyyyy_xyyyyy, \
                             tr_z_xyyyy_xyyyyz, \
                             tr_z_xyyyy_xyyyzz, \
                             tr_z_xyyyy_xyyzzz, \
                             tr_z_xyyyy_xyzzzz, \
                             tr_z_xyyyy_xzzzzz, \
                             tr_z_xyyyy_yyyyyy, \
                             tr_z_xyyyy_yyyyyz, \
                             tr_z_xyyyy_yyyyzz, \
                             tr_z_xyyyy_yyyzzz, \
                             tr_z_xyyyy_yyzzzz, \
                             tr_z_xyyyy_yzzzzz, \
                             tr_z_xyyyy_zzzzzz, \
                             tr_z_yyyy_xxxxx,   \
                             tr_z_yyyy_xxxxxx,  \
                             tr_z_yyyy_xxxxxy,  \
                             tr_z_yyyy_xxxxxz,  \
                             tr_z_yyyy_xxxxy,   \
                             tr_z_yyyy_xxxxyy,  \
                             tr_z_yyyy_xxxxyz,  \
                             tr_z_yyyy_xxxxz,   \
                             tr_z_yyyy_xxxxzz,  \
                             tr_z_yyyy_xxxyy,   \
                             tr_z_yyyy_xxxyyy,  \
                             tr_z_yyyy_xxxyyz,  \
                             tr_z_yyyy_xxxyz,   \
                             tr_z_yyyy_xxxyzz,  \
                             tr_z_yyyy_xxxzz,   \
                             tr_z_yyyy_xxxzzz,  \
                             tr_z_yyyy_xxyyy,   \
                             tr_z_yyyy_xxyyyy,  \
                             tr_z_yyyy_xxyyyz,  \
                             tr_z_yyyy_xxyyz,   \
                             tr_z_yyyy_xxyyzz,  \
                             tr_z_yyyy_xxyzz,   \
                             tr_z_yyyy_xxyzzz,  \
                             tr_z_yyyy_xxzzz,   \
                             tr_z_yyyy_xxzzzz,  \
                             tr_z_yyyy_xyyyy,   \
                             tr_z_yyyy_xyyyyy,  \
                             tr_z_yyyy_xyyyyz,  \
                             tr_z_yyyy_xyyyz,   \
                             tr_z_yyyy_xyyyzz,  \
                             tr_z_yyyy_xyyzz,   \
                             tr_z_yyyy_xyyzzz,  \
                             tr_z_yyyy_xyzzz,   \
                             tr_z_yyyy_xyzzzz,  \
                             tr_z_yyyy_xzzzz,   \
                             tr_z_yyyy_xzzzzz,  \
                             tr_z_yyyy_yyyyy,   \
                             tr_z_yyyy_yyyyyy,  \
                             tr_z_yyyy_yyyyyz,  \
                             tr_z_yyyy_yyyyz,   \
                             tr_z_yyyy_yyyyzz,  \
                             tr_z_yyyy_yyyzz,   \
                             tr_z_yyyy_yyyzzz,  \
                             tr_z_yyyy_yyzzz,   \
                             tr_z_yyyy_yyzzzz,  \
                             tr_z_yyyy_yzzzz,   \
                             tr_z_yyyy_yzzzzz,  \
                             tr_z_yyyy_zzzzz,   \
                             tr_z_yyyy_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyy_xxxxxx[i] = 6.0 * tr_z_yyyy_xxxxx[i] * fe_0 + tr_z_yyyy_xxxxxx[i] * pa_x[i];

        tr_z_xyyyy_xxxxxy[i] = 5.0 * tr_z_yyyy_xxxxy[i] * fe_0 + tr_z_yyyy_xxxxxy[i] * pa_x[i];

        tr_z_xyyyy_xxxxxz[i] = 5.0 * tr_z_yyyy_xxxxz[i] * fe_0 + tr_z_yyyy_xxxxxz[i] * pa_x[i];

        tr_z_xyyyy_xxxxyy[i] = 4.0 * tr_z_yyyy_xxxyy[i] * fe_0 + tr_z_yyyy_xxxxyy[i] * pa_x[i];

        tr_z_xyyyy_xxxxyz[i] = 4.0 * tr_z_yyyy_xxxyz[i] * fe_0 + tr_z_yyyy_xxxxyz[i] * pa_x[i];

        tr_z_xyyyy_xxxxzz[i] = 4.0 * tr_z_yyyy_xxxzz[i] * fe_0 + tr_z_yyyy_xxxxzz[i] * pa_x[i];

        tr_z_xyyyy_xxxyyy[i] = 3.0 * tr_z_yyyy_xxyyy[i] * fe_0 + tr_z_yyyy_xxxyyy[i] * pa_x[i];

        tr_z_xyyyy_xxxyyz[i] = 3.0 * tr_z_yyyy_xxyyz[i] * fe_0 + tr_z_yyyy_xxxyyz[i] * pa_x[i];

        tr_z_xyyyy_xxxyzz[i] = 3.0 * tr_z_yyyy_xxyzz[i] * fe_0 + tr_z_yyyy_xxxyzz[i] * pa_x[i];

        tr_z_xyyyy_xxxzzz[i] = 3.0 * tr_z_yyyy_xxzzz[i] * fe_0 + tr_z_yyyy_xxxzzz[i] * pa_x[i];

        tr_z_xyyyy_xxyyyy[i] = 2.0 * tr_z_yyyy_xyyyy[i] * fe_0 + tr_z_yyyy_xxyyyy[i] * pa_x[i];

        tr_z_xyyyy_xxyyyz[i] = 2.0 * tr_z_yyyy_xyyyz[i] * fe_0 + tr_z_yyyy_xxyyyz[i] * pa_x[i];

        tr_z_xyyyy_xxyyzz[i] = 2.0 * tr_z_yyyy_xyyzz[i] * fe_0 + tr_z_yyyy_xxyyzz[i] * pa_x[i];

        tr_z_xyyyy_xxyzzz[i] = 2.0 * tr_z_yyyy_xyzzz[i] * fe_0 + tr_z_yyyy_xxyzzz[i] * pa_x[i];

        tr_z_xyyyy_xxzzzz[i] = 2.0 * tr_z_yyyy_xzzzz[i] * fe_0 + tr_z_yyyy_xxzzzz[i] * pa_x[i];

        tr_z_xyyyy_xyyyyy[i] = tr_z_yyyy_yyyyy[i] * fe_0 + tr_z_yyyy_xyyyyy[i] * pa_x[i];

        tr_z_xyyyy_xyyyyz[i] = tr_z_yyyy_yyyyz[i] * fe_0 + tr_z_yyyy_xyyyyz[i] * pa_x[i];

        tr_z_xyyyy_xyyyzz[i] = tr_z_yyyy_yyyzz[i] * fe_0 + tr_z_yyyy_xyyyzz[i] * pa_x[i];

        tr_z_xyyyy_xyyzzz[i] = tr_z_yyyy_yyzzz[i] * fe_0 + tr_z_yyyy_xyyzzz[i] * pa_x[i];

        tr_z_xyyyy_xyzzzz[i] = tr_z_yyyy_yzzzz[i] * fe_0 + tr_z_yyyy_xyzzzz[i] * pa_x[i];

        tr_z_xyyyy_xzzzzz[i] = tr_z_yyyy_zzzzz[i] * fe_0 + tr_z_yyyy_xzzzzz[i] * pa_x[i];

        tr_z_xyyyy_yyyyyy[i] = tr_z_yyyy_yyyyyy[i] * pa_x[i];

        tr_z_xyyyy_yyyyyz[i] = tr_z_yyyy_yyyyyz[i] * pa_x[i];

        tr_z_xyyyy_yyyyzz[i] = tr_z_yyyy_yyyyzz[i] * pa_x[i];

        tr_z_xyyyy_yyyzzz[i] = tr_z_yyyy_yyyzzz[i] * pa_x[i];

        tr_z_xyyyy_yyzzzz[i] = tr_z_yyyy_yyzzzz[i] * pa_x[i];

        tr_z_xyyyy_yzzzzz[i] = tr_z_yyyy_yzzzzz[i] * pa_x[i];

        tr_z_xyyyy_zzzzzz[i] = tr_z_yyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 1484-1512 components of targeted buffer : HI

    auto tr_z_xyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 1484);

    auto tr_z_xyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 1485);

    auto tr_z_xyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 1486);

    auto tr_z_xyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 1487);

    auto tr_z_xyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1488);

    auto tr_z_xyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 1489);

    auto tr_z_xyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 1490);

    auto tr_z_xyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1491);

    auto tr_z_xyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1492);

    auto tr_z_xyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 1493);

    auto tr_z_xyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 1494);

    auto tr_z_xyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1495);

    auto tr_z_xyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1496);

    auto tr_z_xyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1497);

    auto tr_z_xyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 1498);

    auto tr_z_xyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 1499);

    auto tr_z_xyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1500);

    auto tr_z_xyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1501);

    auto tr_z_xyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1502);

    auto tr_z_xyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1503);

    auto tr_z_xyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 1504);

    auto tr_z_xyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1505);

    auto tr_z_xyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1506);

    auto tr_z_xyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1507);

    auto tr_z_xyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1508);

    auto tr_z_xyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1509);

    auto tr_z_xyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1510);

    auto tr_z_xyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1511);

#pragma omp simd aligned(pa_x,                  \
                             tr_z_xyyyz_xxxxxx, \
                             tr_z_xyyyz_xxxxxy, \
                             tr_z_xyyyz_xxxxxz, \
                             tr_z_xyyyz_xxxxyy, \
                             tr_z_xyyyz_xxxxyz, \
                             tr_z_xyyyz_xxxxzz, \
                             tr_z_xyyyz_xxxyyy, \
                             tr_z_xyyyz_xxxyyz, \
                             tr_z_xyyyz_xxxyzz, \
                             tr_z_xyyyz_xxxzzz, \
                             tr_z_xyyyz_xxyyyy, \
                             tr_z_xyyyz_xxyyyz, \
                             tr_z_xyyyz_xxyyzz, \
                             tr_z_xyyyz_xxyzzz, \
                             tr_z_xyyyz_xxzzzz, \
                             tr_z_xyyyz_xyyyyy, \
                             tr_z_xyyyz_xyyyyz, \
                             tr_z_xyyyz_xyyyzz, \
                             tr_z_xyyyz_xyyzzz, \
                             tr_z_xyyyz_xyzzzz, \
                             tr_z_xyyyz_xzzzzz, \
                             tr_z_xyyyz_yyyyyy, \
                             tr_z_xyyyz_yyyyyz, \
                             tr_z_xyyyz_yyyyzz, \
                             tr_z_xyyyz_yyyzzz, \
                             tr_z_xyyyz_yyzzzz, \
                             tr_z_xyyyz_yzzzzz, \
                             tr_z_xyyyz_zzzzzz, \
                             tr_z_yyyz_xxxxx,   \
                             tr_z_yyyz_xxxxxx,  \
                             tr_z_yyyz_xxxxxy,  \
                             tr_z_yyyz_xxxxxz,  \
                             tr_z_yyyz_xxxxy,   \
                             tr_z_yyyz_xxxxyy,  \
                             tr_z_yyyz_xxxxyz,  \
                             tr_z_yyyz_xxxxz,   \
                             tr_z_yyyz_xxxxzz,  \
                             tr_z_yyyz_xxxyy,   \
                             tr_z_yyyz_xxxyyy,  \
                             tr_z_yyyz_xxxyyz,  \
                             tr_z_yyyz_xxxyz,   \
                             tr_z_yyyz_xxxyzz,  \
                             tr_z_yyyz_xxxzz,   \
                             tr_z_yyyz_xxxzzz,  \
                             tr_z_yyyz_xxyyy,   \
                             tr_z_yyyz_xxyyyy,  \
                             tr_z_yyyz_xxyyyz,  \
                             tr_z_yyyz_xxyyz,   \
                             tr_z_yyyz_xxyyzz,  \
                             tr_z_yyyz_xxyzz,   \
                             tr_z_yyyz_xxyzzz,  \
                             tr_z_yyyz_xxzzz,   \
                             tr_z_yyyz_xxzzzz,  \
                             tr_z_yyyz_xyyyy,   \
                             tr_z_yyyz_xyyyyy,  \
                             tr_z_yyyz_xyyyyz,  \
                             tr_z_yyyz_xyyyz,   \
                             tr_z_yyyz_xyyyzz,  \
                             tr_z_yyyz_xyyzz,   \
                             tr_z_yyyz_xyyzzz,  \
                             tr_z_yyyz_xyzzz,   \
                             tr_z_yyyz_xyzzzz,  \
                             tr_z_yyyz_xzzzz,   \
                             tr_z_yyyz_xzzzzz,  \
                             tr_z_yyyz_yyyyy,   \
                             tr_z_yyyz_yyyyyy,  \
                             tr_z_yyyz_yyyyyz,  \
                             tr_z_yyyz_yyyyz,   \
                             tr_z_yyyz_yyyyzz,  \
                             tr_z_yyyz_yyyzz,   \
                             tr_z_yyyz_yyyzzz,  \
                             tr_z_yyyz_yyzzz,   \
                             tr_z_yyyz_yyzzzz,  \
                             tr_z_yyyz_yzzzz,   \
                             tr_z_yyyz_yzzzzz,  \
                             tr_z_yyyz_zzzzz,   \
                             tr_z_yyyz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyz_xxxxxx[i] = 6.0 * tr_z_yyyz_xxxxx[i] * fe_0 + tr_z_yyyz_xxxxxx[i] * pa_x[i];

        tr_z_xyyyz_xxxxxy[i] = 5.0 * tr_z_yyyz_xxxxy[i] * fe_0 + tr_z_yyyz_xxxxxy[i] * pa_x[i];

        tr_z_xyyyz_xxxxxz[i] = 5.0 * tr_z_yyyz_xxxxz[i] * fe_0 + tr_z_yyyz_xxxxxz[i] * pa_x[i];

        tr_z_xyyyz_xxxxyy[i] = 4.0 * tr_z_yyyz_xxxyy[i] * fe_0 + tr_z_yyyz_xxxxyy[i] * pa_x[i];

        tr_z_xyyyz_xxxxyz[i] = 4.0 * tr_z_yyyz_xxxyz[i] * fe_0 + tr_z_yyyz_xxxxyz[i] * pa_x[i];

        tr_z_xyyyz_xxxxzz[i] = 4.0 * tr_z_yyyz_xxxzz[i] * fe_0 + tr_z_yyyz_xxxxzz[i] * pa_x[i];

        tr_z_xyyyz_xxxyyy[i] = 3.0 * tr_z_yyyz_xxyyy[i] * fe_0 + tr_z_yyyz_xxxyyy[i] * pa_x[i];

        tr_z_xyyyz_xxxyyz[i] = 3.0 * tr_z_yyyz_xxyyz[i] * fe_0 + tr_z_yyyz_xxxyyz[i] * pa_x[i];

        tr_z_xyyyz_xxxyzz[i] = 3.0 * tr_z_yyyz_xxyzz[i] * fe_0 + tr_z_yyyz_xxxyzz[i] * pa_x[i];

        tr_z_xyyyz_xxxzzz[i] = 3.0 * tr_z_yyyz_xxzzz[i] * fe_0 + tr_z_yyyz_xxxzzz[i] * pa_x[i];

        tr_z_xyyyz_xxyyyy[i] = 2.0 * tr_z_yyyz_xyyyy[i] * fe_0 + tr_z_yyyz_xxyyyy[i] * pa_x[i];

        tr_z_xyyyz_xxyyyz[i] = 2.0 * tr_z_yyyz_xyyyz[i] * fe_0 + tr_z_yyyz_xxyyyz[i] * pa_x[i];

        tr_z_xyyyz_xxyyzz[i] = 2.0 * tr_z_yyyz_xyyzz[i] * fe_0 + tr_z_yyyz_xxyyzz[i] * pa_x[i];

        tr_z_xyyyz_xxyzzz[i] = 2.0 * tr_z_yyyz_xyzzz[i] * fe_0 + tr_z_yyyz_xxyzzz[i] * pa_x[i];

        tr_z_xyyyz_xxzzzz[i] = 2.0 * tr_z_yyyz_xzzzz[i] * fe_0 + tr_z_yyyz_xxzzzz[i] * pa_x[i];

        tr_z_xyyyz_xyyyyy[i] = tr_z_yyyz_yyyyy[i] * fe_0 + tr_z_yyyz_xyyyyy[i] * pa_x[i];

        tr_z_xyyyz_xyyyyz[i] = tr_z_yyyz_yyyyz[i] * fe_0 + tr_z_yyyz_xyyyyz[i] * pa_x[i];

        tr_z_xyyyz_xyyyzz[i] = tr_z_yyyz_yyyzz[i] * fe_0 + tr_z_yyyz_xyyyzz[i] * pa_x[i];

        tr_z_xyyyz_xyyzzz[i] = tr_z_yyyz_yyzzz[i] * fe_0 + tr_z_yyyz_xyyzzz[i] * pa_x[i];

        tr_z_xyyyz_xyzzzz[i] = tr_z_yyyz_yzzzz[i] * fe_0 + tr_z_yyyz_xyzzzz[i] * pa_x[i];

        tr_z_xyyyz_xzzzzz[i] = tr_z_yyyz_zzzzz[i] * fe_0 + tr_z_yyyz_xzzzzz[i] * pa_x[i];

        tr_z_xyyyz_yyyyyy[i] = tr_z_yyyz_yyyyyy[i] * pa_x[i];

        tr_z_xyyyz_yyyyyz[i] = tr_z_yyyz_yyyyyz[i] * pa_x[i];

        tr_z_xyyyz_yyyyzz[i] = tr_z_yyyz_yyyyzz[i] * pa_x[i];

        tr_z_xyyyz_yyyzzz[i] = tr_z_yyyz_yyyzzz[i] * pa_x[i];

        tr_z_xyyyz_yyzzzz[i] = tr_z_yyyz_yyzzzz[i] * pa_x[i];

        tr_z_xyyyz_yzzzzz[i] = tr_z_yyyz_yzzzzz[i] * pa_x[i];

        tr_z_xyyyz_zzzzzz[i] = tr_z_yyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1512-1540 components of targeted buffer : HI

    auto tr_z_xyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 1512);

    auto tr_z_xyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 1513);

    auto tr_z_xyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 1514);

    auto tr_z_xyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 1515);

    auto tr_z_xyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 1516);

    auto tr_z_xyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 1517);

    auto tr_z_xyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 1518);

    auto tr_z_xyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 1519);

    auto tr_z_xyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 1520);

    auto tr_z_xyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 1521);

    auto tr_z_xyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 1522);

    auto tr_z_xyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 1523);

    auto tr_z_xyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 1524);

    auto tr_z_xyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 1525);

    auto tr_z_xyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 1526);

    auto tr_z_xyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 1527);

    auto tr_z_xyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 1528);

    auto tr_z_xyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 1529);

    auto tr_z_xyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 1530);

    auto tr_z_xyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 1531);

    auto tr_z_xyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 1532);

    auto tr_z_xyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1533);

    auto tr_z_xyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1534);

    auto tr_z_xyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1535);

    auto tr_z_xyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1536);

    auto tr_z_xyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1537);

    auto tr_z_xyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1538);

    auto tr_z_xyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 1539);

#pragma omp simd aligned(pa_x,                  \
                             tr_z_xyyzz_xxxxxx, \
                             tr_z_xyyzz_xxxxxy, \
                             tr_z_xyyzz_xxxxxz, \
                             tr_z_xyyzz_xxxxyy, \
                             tr_z_xyyzz_xxxxyz, \
                             tr_z_xyyzz_xxxxzz, \
                             tr_z_xyyzz_xxxyyy, \
                             tr_z_xyyzz_xxxyyz, \
                             tr_z_xyyzz_xxxyzz, \
                             tr_z_xyyzz_xxxzzz, \
                             tr_z_xyyzz_xxyyyy, \
                             tr_z_xyyzz_xxyyyz, \
                             tr_z_xyyzz_xxyyzz, \
                             tr_z_xyyzz_xxyzzz, \
                             tr_z_xyyzz_xxzzzz, \
                             tr_z_xyyzz_xyyyyy, \
                             tr_z_xyyzz_xyyyyz, \
                             tr_z_xyyzz_xyyyzz, \
                             tr_z_xyyzz_xyyzzz, \
                             tr_z_xyyzz_xyzzzz, \
                             tr_z_xyyzz_xzzzzz, \
                             tr_z_xyyzz_yyyyyy, \
                             tr_z_xyyzz_yyyyyz, \
                             tr_z_xyyzz_yyyyzz, \
                             tr_z_xyyzz_yyyzzz, \
                             tr_z_xyyzz_yyzzzz, \
                             tr_z_xyyzz_yzzzzz, \
                             tr_z_xyyzz_zzzzzz, \
                             tr_z_yyzz_xxxxx,   \
                             tr_z_yyzz_xxxxxx,  \
                             tr_z_yyzz_xxxxxy,  \
                             tr_z_yyzz_xxxxxz,  \
                             tr_z_yyzz_xxxxy,   \
                             tr_z_yyzz_xxxxyy,  \
                             tr_z_yyzz_xxxxyz,  \
                             tr_z_yyzz_xxxxz,   \
                             tr_z_yyzz_xxxxzz,  \
                             tr_z_yyzz_xxxyy,   \
                             tr_z_yyzz_xxxyyy,  \
                             tr_z_yyzz_xxxyyz,  \
                             tr_z_yyzz_xxxyz,   \
                             tr_z_yyzz_xxxyzz,  \
                             tr_z_yyzz_xxxzz,   \
                             tr_z_yyzz_xxxzzz,  \
                             tr_z_yyzz_xxyyy,   \
                             tr_z_yyzz_xxyyyy,  \
                             tr_z_yyzz_xxyyyz,  \
                             tr_z_yyzz_xxyyz,   \
                             tr_z_yyzz_xxyyzz,  \
                             tr_z_yyzz_xxyzz,   \
                             tr_z_yyzz_xxyzzz,  \
                             tr_z_yyzz_xxzzz,   \
                             tr_z_yyzz_xxzzzz,  \
                             tr_z_yyzz_xyyyy,   \
                             tr_z_yyzz_xyyyyy,  \
                             tr_z_yyzz_xyyyyz,  \
                             tr_z_yyzz_xyyyz,   \
                             tr_z_yyzz_xyyyzz,  \
                             tr_z_yyzz_xyyzz,   \
                             tr_z_yyzz_xyyzzz,  \
                             tr_z_yyzz_xyzzz,   \
                             tr_z_yyzz_xyzzzz,  \
                             tr_z_yyzz_xzzzz,   \
                             tr_z_yyzz_xzzzzz,  \
                             tr_z_yyzz_yyyyy,   \
                             tr_z_yyzz_yyyyyy,  \
                             tr_z_yyzz_yyyyyz,  \
                             tr_z_yyzz_yyyyz,   \
                             tr_z_yyzz_yyyyzz,  \
                             tr_z_yyzz_yyyzz,   \
                             tr_z_yyzz_yyyzzz,  \
                             tr_z_yyzz_yyzzz,   \
                             tr_z_yyzz_yyzzzz,  \
                             tr_z_yyzz_yzzzz,   \
                             tr_z_yyzz_yzzzzz,  \
                             tr_z_yyzz_zzzzz,   \
                             tr_z_yyzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzz_xxxxxx[i] = 6.0 * tr_z_yyzz_xxxxx[i] * fe_0 + tr_z_yyzz_xxxxxx[i] * pa_x[i];

        tr_z_xyyzz_xxxxxy[i] = 5.0 * tr_z_yyzz_xxxxy[i] * fe_0 + tr_z_yyzz_xxxxxy[i] * pa_x[i];

        tr_z_xyyzz_xxxxxz[i] = 5.0 * tr_z_yyzz_xxxxz[i] * fe_0 + tr_z_yyzz_xxxxxz[i] * pa_x[i];

        tr_z_xyyzz_xxxxyy[i] = 4.0 * tr_z_yyzz_xxxyy[i] * fe_0 + tr_z_yyzz_xxxxyy[i] * pa_x[i];

        tr_z_xyyzz_xxxxyz[i] = 4.0 * tr_z_yyzz_xxxyz[i] * fe_0 + tr_z_yyzz_xxxxyz[i] * pa_x[i];

        tr_z_xyyzz_xxxxzz[i] = 4.0 * tr_z_yyzz_xxxzz[i] * fe_0 + tr_z_yyzz_xxxxzz[i] * pa_x[i];

        tr_z_xyyzz_xxxyyy[i] = 3.0 * tr_z_yyzz_xxyyy[i] * fe_0 + tr_z_yyzz_xxxyyy[i] * pa_x[i];

        tr_z_xyyzz_xxxyyz[i] = 3.0 * tr_z_yyzz_xxyyz[i] * fe_0 + tr_z_yyzz_xxxyyz[i] * pa_x[i];

        tr_z_xyyzz_xxxyzz[i] = 3.0 * tr_z_yyzz_xxyzz[i] * fe_0 + tr_z_yyzz_xxxyzz[i] * pa_x[i];

        tr_z_xyyzz_xxxzzz[i] = 3.0 * tr_z_yyzz_xxzzz[i] * fe_0 + tr_z_yyzz_xxxzzz[i] * pa_x[i];

        tr_z_xyyzz_xxyyyy[i] = 2.0 * tr_z_yyzz_xyyyy[i] * fe_0 + tr_z_yyzz_xxyyyy[i] * pa_x[i];

        tr_z_xyyzz_xxyyyz[i] = 2.0 * tr_z_yyzz_xyyyz[i] * fe_0 + tr_z_yyzz_xxyyyz[i] * pa_x[i];

        tr_z_xyyzz_xxyyzz[i] = 2.0 * tr_z_yyzz_xyyzz[i] * fe_0 + tr_z_yyzz_xxyyzz[i] * pa_x[i];

        tr_z_xyyzz_xxyzzz[i] = 2.0 * tr_z_yyzz_xyzzz[i] * fe_0 + tr_z_yyzz_xxyzzz[i] * pa_x[i];

        tr_z_xyyzz_xxzzzz[i] = 2.0 * tr_z_yyzz_xzzzz[i] * fe_0 + tr_z_yyzz_xxzzzz[i] * pa_x[i];

        tr_z_xyyzz_xyyyyy[i] = tr_z_yyzz_yyyyy[i] * fe_0 + tr_z_yyzz_xyyyyy[i] * pa_x[i];

        tr_z_xyyzz_xyyyyz[i] = tr_z_yyzz_yyyyz[i] * fe_0 + tr_z_yyzz_xyyyyz[i] * pa_x[i];

        tr_z_xyyzz_xyyyzz[i] = tr_z_yyzz_yyyzz[i] * fe_0 + tr_z_yyzz_xyyyzz[i] * pa_x[i];

        tr_z_xyyzz_xyyzzz[i] = tr_z_yyzz_yyzzz[i] * fe_0 + tr_z_yyzz_xyyzzz[i] * pa_x[i];

        tr_z_xyyzz_xyzzzz[i] = tr_z_yyzz_yzzzz[i] * fe_0 + tr_z_yyzz_xyzzzz[i] * pa_x[i];

        tr_z_xyyzz_xzzzzz[i] = tr_z_yyzz_zzzzz[i] * fe_0 + tr_z_yyzz_xzzzzz[i] * pa_x[i];

        tr_z_xyyzz_yyyyyy[i] = tr_z_yyzz_yyyyyy[i] * pa_x[i];

        tr_z_xyyzz_yyyyyz[i] = tr_z_yyzz_yyyyyz[i] * pa_x[i];

        tr_z_xyyzz_yyyyzz[i] = tr_z_yyzz_yyyyzz[i] * pa_x[i];

        tr_z_xyyzz_yyyzzz[i] = tr_z_yyzz_yyyzzz[i] * pa_x[i];

        tr_z_xyyzz_yyzzzz[i] = tr_z_yyzz_yyzzzz[i] * pa_x[i];

        tr_z_xyyzz_yzzzzz[i] = tr_z_yyzz_yzzzzz[i] * pa_x[i];

        tr_z_xyyzz_zzzzzz[i] = tr_z_yyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1540-1568 components of targeted buffer : HI

    auto tr_z_xyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1540);

    auto tr_z_xyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1541);

    auto tr_z_xyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1542);

    auto tr_z_xyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1543);

    auto tr_z_xyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1544);

    auto tr_z_xyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1545);

    auto tr_z_xyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1546);

    auto tr_z_xyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1547);

    auto tr_z_xyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1548);

    auto tr_z_xyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1549);

    auto tr_z_xyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1550);

    auto tr_z_xyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1551);

    auto tr_z_xyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1552);

    auto tr_z_xyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1553);

    auto tr_z_xyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1554);

    auto tr_z_xyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1555);

    auto tr_z_xyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1556);

    auto tr_z_xyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1557);

    auto tr_z_xyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1558);

    auto tr_z_xyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1559);

    auto tr_z_xyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1560);

    auto tr_z_xyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1561);

    auto tr_z_xyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1562);

    auto tr_z_xyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1563);

    auto tr_z_xyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1564);

    auto tr_z_xyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1565);

    auto tr_z_xyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1566);

    auto tr_z_xyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1567);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             tr_z_xyzzz_xxxxxx, \
                             tr_z_xyzzz_xxxxxy, \
                             tr_z_xyzzz_xxxxxz, \
                             tr_z_xyzzz_xxxxyy, \
                             tr_z_xyzzz_xxxxyz, \
                             tr_z_xyzzz_xxxxzz, \
                             tr_z_xyzzz_xxxyyy, \
                             tr_z_xyzzz_xxxyyz, \
                             tr_z_xyzzz_xxxyzz, \
                             tr_z_xyzzz_xxxzzz, \
                             tr_z_xyzzz_xxyyyy, \
                             tr_z_xyzzz_xxyyyz, \
                             tr_z_xyzzz_xxyyzz, \
                             tr_z_xyzzz_xxyzzz, \
                             tr_z_xyzzz_xxzzzz, \
                             tr_z_xyzzz_xyyyyy, \
                             tr_z_xyzzz_xyyyyz, \
                             tr_z_xyzzz_xyyyzz, \
                             tr_z_xyzzz_xyyzzz, \
                             tr_z_xyzzz_xyzzzz, \
                             tr_z_xyzzz_xzzzzz, \
                             tr_z_xyzzz_yyyyyy, \
                             tr_z_xyzzz_yyyyyz, \
                             tr_z_xyzzz_yyyyzz, \
                             tr_z_xyzzz_yyyzzz, \
                             tr_z_xyzzz_yyzzzz, \
                             tr_z_xyzzz_yzzzzz, \
                             tr_z_xyzzz_zzzzzz, \
                             tr_z_xzzz_xxxxxx,  \
                             tr_z_xzzz_xxxxxz,  \
                             tr_z_xzzz_xxxxzz,  \
                             tr_z_xzzz_xxxzzz,  \
                             tr_z_xzzz_xxzzzz,  \
                             tr_z_xzzz_xzzzzz,  \
                             tr_z_yzzz_xxxxxy,  \
                             tr_z_yzzz_xxxxy,   \
                             tr_z_yzzz_xxxxyy,  \
                             tr_z_yzzz_xxxxyz,  \
                             tr_z_yzzz_xxxyy,   \
                             tr_z_yzzz_xxxyyy,  \
                             tr_z_yzzz_xxxyyz,  \
                             tr_z_yzzz_xxxyz,   \
                             tr_z_yzzz_xxxyzz,  \
                             tr_z_yzzz_xxyyy,   \
                             tr_z_yzzz_xxyyyy,  \
                             tr_z_yzzz_xxyyyz,  \
                             tr_z_yzzz_xxyyz,   \
                             tr_z_yzzz_xxyyzz,  \
                             tr_z_yzzz_xxyzz,   \
                             tr_z_yzzz_xxyzzz,  \
                             tr_z_yzzz_xyyyy,   \
                             tr_z_yzzz_xyyyyy,  \
                             tr_z_yzzz_xyyyyz,  \
                             tr_z_yzzz_xyyyz,   \
                             tr_z_yzzz_xyyyzz,  \
                             tr_z_yzzz_xyyzz,   \
                             tr_z_yzzz_xyyzzz,  \
                             tr_z_yzzz_xyzzz,   \
                             tr_z_yzzz_xyzzzz,  \
                             tr_z_yzzz_yyyyy,   \
                             tr_z_yzzz_yyyyyy,  \
                             tr_z_yzzz_yyyyyz,  \
                             tr_z_yzzz_yyyyz,   \
                             tr_z_yzzz_yyyyzz,  \
                             tr_z_yzzz_yyyzz,   \
                             tr_z_yzzz_yyyzzz,  \
                             tr_z_yzzz_yyzzz,   \
                             tr_z_yzzz_yyzzzz,  \
                             tr_z_yzzz_yzzzz,   \
                             tr_z_yzzz_yzzzzz,  \
                             tr_z_yzzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzz_xxxxxx[i] = tr_z_xzzz_xxxxxx[i] * pa_y[i];

        tr_z_xyzzz_xxxxxy[i] = 5.0 * tr_z_yzzz_xxxxy[i] * fe_0 + tr_z_yzzz_xxxxxy[i] * pa_x[i];

        tr_z_xyzzz_xxxxxz[i] = tr_z_xzzz_xxxxxz[i] * pa_y[i];

        tr_z_xyzzz_xxxxyy[i] = 4.0 * tr_z_yzzz_xxxyy[i] * fe_0 + tr_z_yzzz_xxxxyy[i] * pa_x[i];

        tr_z_xyzzz_xxxxyz[i] = 4.0 * tr_z_yzzz_xxxyz[i] * fe_0 + tr_z_yzzz_xxxxyz[i] * pa_x[i];

        tr_z_xyzzz_xxxxzz[i] = tr_z_xzzz_xxxxzz[i] * pa_y[i];

        tr_z_xyzzz_xxxyyy[i] = 3.0 * tr_z_yzzz_xxyyy[i] * fe_0 + tr_z_yzzz_xxxyyy[i] * pa_x[i];

        tr_z_xyzzz_xxxyyz[i] = 3.0 * tr_z_yzzz_xxyyz[i] * fe_0 + tr_z_yzzz_xxxyyz[i] * pa_x[i];

        tr_z_xyzzz_xxxyzz[i] = 3.0 * tr_z_yzzz_xxyzz[i] * fe_0 + tr_z_yzzz_xxxyzz[i] * pa_x[i];

        tr_z_xyzzz_xxxzzz[i] = tr_z_xzzz_xxxzzz[i] * pa_y[i];

        tr_z_xyzzz_xxyyyy[i] = 2.0 * tr_z_yzzz_xyyyy[i] * fe_0 + tr_z_yzzz_xxyyyy[i] * pa_x[i];

        tr_z_xyzzz_xxyyyz[i] = 2.0 * tr_z_yzzz_xyyyz[i] * fe_0 + tr_z_yzzz_xxyyyz[i] * pa_x[i];

        tr_z_xyzzz_xxyyzz[i] = 2.0 * tr_z_yzzz_xyyzz[i] * fe_0 + tr_z_yzzz_xxyyzz[i] * pa_x[i];

        tr_z_xyzzz_xxyzzz[i] = 2.0 * tr_z_yzzz_xyzzz[i] * fe_0 + tr_z_yzzz_xxyzzz[i] * pa_x[i];

        tr_z_xyzzz_xxzzzz[i] = tr_z_xzzz_xxzzzz[i] * pa_y[i];

        tr_z_xyzzz_xyyyyy[i] = tr_z_yzzz_yyyyy[i] * fe_0 + tr_z_yzzz_xyyyyy[i] * pa_x[i];

        tr_z_xyzzz_xyyyyz[i] = tr_z_yzzz_yyyyz[i] * fe_0 + tr_z_yzzz_xyyyyz[i] * pa_x[i];

        tr_z_xyzzz_xyyyzz[i] = tr_z_yzzz_yyyzz[i] * fe_0 + tr_z_yzzz_xyyyzz[i] * pa_x[i];

        tr_z_xyzzz_xyyzzz[i] = tr_z_yzzz_yyzzz[i] * fe_0 + tr_z_yzzz_xyyzzz[i] * pa_x[i];

        tr_z_xyzzz_xyzzzz[i] = tr_z_yzzz_yzzzz[i] * fe_0 + tr_z_yzzz_xyzzzz[i] * pa_x[i];

        tr_z_xyzzz_xzzzzz[i] = tr_z_xzzz_xzzzzz[i] * pa_y[i];

        tr_z_xyzzz_yyyyyy[i] = tr_z_yzzz_yyyyyy[i] * pa_x[i];

        tr_z_xyzzz_yyyyyz[i] = tr_z_yzzz_yyyyyz[i] * pa_x[i];

        tr_z_xyzzz_yyyyzz[i] = tr_z_yzzz_yyyyzz[i] * pa_x[i];

        tr_z_xyzzz_yyyzzz[i] = tr_z_yzzz_yyyzzz[i] * pa_x[i];

        tr_z_xyzzz_yyzzzz[i] = tr_z_yzzz_yyzzzz[i] * pa_x[i];

        tr_z_xyzzz_yzzzzz[i] = tr_z_yzzz_yzzzzz[i] * pa_x[i];

        tr_z_xyzzz_zzzzzz[i] = tr_z_yzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1568-1596 components of targeted buffer : HI

    auto tr_z_xzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1568);

    auto tr_z_xzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1569);

    auto tr_z_xzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1570);

    auto tr_z_xzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1571);

    auto tr_z_xzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1572);

    auto tr_z_xzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1573);

    auto tr_z_xzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1574);

    auto tr_z_xzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1575);

    auto tr_z_xzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1576);

    auto tr_z_xzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1577);

    auto tr_z_xzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1578);

    auto tr_z_xzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1579);

    auto tr_z_xzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1580);

    auto tr_z_xzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1581);

    auto tr_z_xzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1582);

    auto tr_z_xzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1583);

    auto tr_z_xzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1584);

    auto tr_z_xzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1585);

    auto tr_z_xzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1586);

    auto tr_z_xzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1587);

    auto tr_z_xzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1588);

    auto tr_z_xzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1589);

    auto tr_z_xzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1590);

    auto tr_z_xzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1591);

    auto tr_z_xzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1592);

    auto tr_z_xzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1593);

    auto tr_z_xzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1594);

    auto tr_z_xzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1595);

#pragma omp simd aligned(pa_x,                  \
                             tr_z_xzzzz_xxxxxx, \
                             tr_z_xzzzz_xxxxxy, \
                             tr_z_xzzzz_xxxxxz, \
                             tr_z_xzzzz_xxxxyy, \
                             tr_z_xzzzz_xxxxyz, \
                             tr_z_xzzzz_xxxxzz, \
                             tr_z_xzzzz_xxxyyy, \
                             tr_z_xzzzz_xxxyyz, \
                             tr_z_xzzzz_xxxyzz, \
                             tr_z_xzzzz_xxxzzz, \
                             tr_z_xzzzz_xxyyyy, \
                             tr_z_xzzzz_xxyyyz, \
                             tr_z_xzzzz_xxyyzz, \
                             tr_z_xzzzz_xxyzzz, \
                             tr_z_xzzzz_xxzzzz, \
                             tr_z_xzzzz_xyyyyy, \
                             tr_z_xzzzz_xyyyyz, \
                             tr_z_xzzzz_xyyyzz, \
                             tr_z_xzzzz_xyyzzz, \
                             tr_z_xzzzz_xyzzzz, \
                             tr_z_xzzzz_xzzzzz, \
                             tr_z_xzzzz_yyyyyy, \
                             tr_z_xzzzz_yyyyyz, \
                             tr_z_xzzzz_yyyyzz, \
                             tr_z_xzzzz_yyyzzz, \
                             tr_z_xzzzz_yyzzzz, \
                             tr_z_xzzzz_yzzzzz, \
                             tr_z_xzzzz_zzzzzz, \
                             tr_z_zzzz_xxxxx,   \
                             tr_z_zzzz_xxxxxx,  \
                             tr_z_zzzz_xxxxxy,  \
                             tr_z_zzzz_xxxxxz,  \
                             tr_z_zzzz_xxxxy,   \
                             tr_z_zzzz_xxxxyy,  \
                             tr_z_zzzz_xxxxyz,  \
                             tr_z_zzzz_xxxxz,   \
                             tr_z_zzzz_xxxxzz,  \
                             tr_z_zzzz_xxxyy,   \
                             tr_z_zzzz_xxxyyy,  \
                             tr_z_zzzz_xxxyyz,  \
                             tr_z_zzzz_xxxyz,   \
                             tr_z_zzzz_xxxyzz,  \
                             tr_z_zzzz_xxxzz,   \
                             tr_z_zzzz_xxxzzz,  \
                             tr_z_zzzz_xxyyy,   \
                             tr_z_zzzz_xxyyyy,  \
                             tr_z_zzzz_xxyyyz,  \
                             tr_z_zzzz_xxyyz,   \
                             tr_z_zzzz_xxyyzz,  \
                             tr_z_zzzz_xxyzz,   \
                             tr_z_zzzz_xxyzzz,  \
                             tr_z_zzzz_xxzzz,   \
                             tr_z_zzzz_xxzzzz,  \
                             tr_z_zzzz_xyyyy,   \
                             tr_z_zzzz_xyyyyy,  \
                             tr_z_zzzz_xyyyyz,  \
                             tr_z_zzzz_xyyyz,   \
                             tr_z_zzzz_xyyyzz,  \
                             tr_z_zzzz_xyyzz,   \
                             tr_z_zzzz_xyyzzz,  \
                             tr_z_zzzz_xyzzz,   \
                             tr_z_zzzz_xyzzzz,  \
                             tr_z_zzzz_xzzzz,   \
                             tr_z_zzzz_xzzzzz,  \
                             tr_z_zzzz_yyyyy,   \
                             tr_z_zzzz_yyyyyy,  \
                             tr_z_zzzz_yyyyyz,  \
                             tr_z_zzzz_yyyyz,   \
                             tr_z_zzzz_yyyyzz,  \
                             tr_z_zzzz_yyyzz,   \
                             tr_z_zzzz_yyyzzz,  \
                             tr_z_zzzz_yyzzz,   \
                             tr_z_zzzz_yyzzzz,  \
                             tr_z_zzzz_yzzzz,   \
                             tr_z_zzzz_yzzzzz,  \
                             tr_z_zzzz_zzzzz,   \
                             tr_z_zzzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzz_xxxxxx[i] = 6.0 * tr_z_zzzz_xxxxx[i] * fe_0 + tr_z_zzzz_xxxxxx[i] * pa_x[i];

        tr_z_xzzzz_xxxxxy[i] = 5.0 * tr_z_zzzz_xxxxy[i] * fe_0 + tr_z_zzzz_xxxxxy[i] * pa_x[i];

        tr_z_xzzzz_xxxxxz[i] = 5.0 * tr_z_zzzz_xxxxz[i] * fe_0 + tr_z_zzzz_xxxxxz[i] * pa_x[i];

        tr_z_xzzzz_xxxxyy[i] = 4.0 * tr_z_zzzz_xxxyy[i] * fe_0 + tr_z_zzzz_xxxxyy[i] * pa_x[i];

        tr_z_xzzzz_xxxxyz[i] = 4.0 * tr_z_zzzz_xxxyz[i] * fe_0 + tr_z_zzzz_xxxxyz[i] * pa_x[i];

        tr_z_xzzzz_xxxxzz[i] = 4.0 * tr_z_zzzz_xxxzz[i] * fe_0 + tr_z_zzzz_xxxxzz[i] * pa_x[i];

        tr_z_xzzzz_xxxyyy[i] = 3.0 * tr_z_zzzz_xxyyy[i] * fe_0 + tr_z_zzzz_xxxyyy[i] * pa_x[i];

        tr_z_xzzzz_xxxyyz[i] = 3.0 * tr_z_zzzz_xxyyz[i] * fe_0 + tr_z_zzzz_xxxyyz[i] * pa_x[i];

        tr_z_xzzzz_xxxyzz[i] = 3.0 * tr_z_zzzz_xxyzz[i] * fe_0 + tr_z_zzzz_xxxyzz[i] * pa_x[i];

        tr_z_xzzzz_xxxzzz[i] = 3.0 * tr_z_zzzz_xxzzz[i] * fe_0 + tr_z_zzzz_xxxzzz[i] * pa_x[i];

        tr_z_xzzzz_xxyyyy[i] = 2.0 * tr_z_zzzz_xyyyy[i] * fe_0 + tr_z_zzzz_xxyyyy[i] * pa_x[i];

        tr_z_xzzzz_xxyyyz[i] = 2.0 * tr_z_zzzz_xyyyz[i] * fe_0 + tr_z_zzzz_xxyyyz[i] * pa_x[i];

        tr_z_xzzzz_xxyyzz[i] = 2.0 * tr_z_zzzz_xyyzz[i] * fe_0 + tr_z_zzzz_xxyyzz[i] * pa_x[i];

        tr_z_xzzzz_xxyzzz[i] = 2.0 * tr_z_zzzz_xyzzz[i] * fe_0 + tr_z_zzzz_xxyzzz[i] * pa_x[i];

        tr_z_xzzzz_xxzzzz[i] = 2.0 * tr_z_zzzz_xzzzz[i] * fe_0 + tr_z_zzzz_xxzzzz[i] * pa_x[i];

        tr_z_xzzzz_xyyyyy[i] = tr_z_zzzz_yyyyy[i] * fe_0 + tr_z_zzzz_xyyyyy[i] * pa_x[i];

        tr_z_xzzzz_xyyyyz[i] = tr_z_zzzz_yyyyz[i] * fe_0 + tr_z_zzzz_xyyyyz[i] * pa_x[i];

        tr_z_xzzzz_xyyyzz[i] = tr_z_zzzz_yyyzz[i] * fe_0 + tr_z_zzzz_xyyyzz[i] * pa_x[i];

        tr_z_xzzzz_xyyzzz[i] = tr_z_zzzz_yyzzz[i] * fe_0 + tr_z_zzzz_xyyzzz[i] * pa_x[i];

        tr_z_xzzzz_xyzzzz[i] = tr_z_zzzz_yzzzz[i] * fe_0 + tr_z_zzzz_xyzzzz[i] * pa_x[i];

        tr_z_xzzzz_xzzzzz[i] = tr_z_zzzz_zzzzz[i] * fe_0 + tr_z_zzzz_xzzzzz[i] * pa_x[i];

        tr_z_xzzzz_yyyyyy[i] = tr_z_zzzz_yyyyyy[i] * pa_x[i];

        tr_z_xzzzz_yyyyyz[i] = tr_z_zzzz_yyyyyz[i] * pa_x[i];

        tr_z_xzzzz_yyyyzz[i] = tr_z_zzzz_yyyyzz[i] * pa_x[i];

        tr_z_xzzzz_yyyzzz[i] = tr_z_zzzz_yyyzzz[i] * pa_x[i];

        tr_z_xzzzz_yyzzzz[i] = tr_z_zzzz_yyzzzz[i] * pa_x[i];

        tr_z_xzzzz_yzzzzz[i] = tr_z_zzzz_yzzzzz[i] * pa_x[i];

        tr_z_xzzzz_zzzzzz[i] = tr_z_zzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 1596-1624 components of targeted buffer : HI

    auto tr_z_yyyyy_xxxxxx = pbuffer.data(idx_dip_hi + 1596);

    auto tr_z_yyyyy_xxxxxy = pbuffer.data(idx_dip_hi + 1597);

    auto tr_z_yyyyy_xxxxxz = pbuffer.data(idx_dip_hi + 1598);

    auto tr_z_yyyyy_xxxxyy = pbuffer.data(idx_dip_hi + 1599);

    auto tr_z_yyyyy_xxxxyz = pbuffer.data(idx_dip_hi + 1600);

    auto tr_z_yyyyy_xxxxzz = pbuffer.data(idx_dip_hi + 1601);

    auto tr_z_yyyyy_xxxyyy = pbuffer.data(idx_dip_hi + 1602);

    auto tr_z_yyyyy_xxxyyz = pbuffer.data(idx_dip_hi + 1603);

    auto tr_z_yyyyy_xxxyzz = pbuffer.data(idx_dip_hi + 1604);

    auto tr_z_yyyyy_xxxzzz = pbuffer.data(idx_dip_hi + 1605);

    auto tr_z_yyyyy_xxyyyy = pbuffer.data(idx_dip_hi + 1606);

    auto tr_z_yyyyy_xxyyyz = pbuffer.data(idx_dip_hi + 1607);

    auto tr_z_yyyyy_xxyyzz = pbuffer.data(idx_dip_hi + 1608);

    auto tr_z_yyyyy_xxyzzz = pbuffer.data(idx_dip_hi + 1609);

    auto tr_z_yyyyy_xxzzzz = pbuffer.data(idx_dip_hi + 1610);

    auto tr_z_yyyyy_xyyyyy = pbuffer.data(idx_dip_hi + 1611);

    auto tr_z_yyyyy_xyyyyz = pbuffer.data(idx_dip_hi + 1612);

    auto tr_z_yyyyy_xyyyzz = pbuffer.data(idx_dip_hi + 1613);

    auto tr_z_yyyyy_xyyzzz = pbuffer.data(idx_dip_hi + 1614);

    auto tr_z_yyyyy_xyzzzz = pbuffer.data(idx_dip_hi + 1615);

    auto tr_z_yyyyy_xzzzzz = pbuffer.data(idx_dip_hi + 1616);

    auto tr_z_yyyyy_yyyyyy = pbuffer.data(idx_dip_hi + 1617);

    auto tr_z_yyyyy_yyyyyz = pbuffer.data(idx_dip_hi + 1618);

    auto tr_z_yyyyy_yyyyzz = pbuffer.data(idx_dip_hi + 1619);

    auto tr_z_yyyyy_yyyzzz = pbuffer.data(idx_dip_hi + 1620);

    auto tr_z_yyyyy_yyzzzz = pbuffer.data(idx_dip_hi + 1621);

    auto tr_z_yyyyy_yzzzzz = pbuffer.data(idx_dip_hi + 1622);

    auto tr_z_yyyyy_zzzzzz = pbuffer.data(idx_dip_hi + 1623);

#pragma omp simd aligned(pa_y,                  \
                             tr_z_yyy_xxxxxx,   \
                             tr_z_yyy_xxxxxy,   \
                             tr_z_yyy_xxxxxz,   \
                             tr_z_yyy_xxxxyy,   \
                             tr_z_yyy_xxxxyz,   \
                             tr_z_yyy_xxxxzz,   \
                             tr_z_yyy_xxxyyy,   \
                             tr_z_yyy_xxxyyz,   \
                             tr_z_yyy_xxxyzz,   \
                             tr_z_yyy_xxxzzz,   \
                             tr_z_yyy_xxyyyy,   \
                             tr_z_yyy_xxyyyz,   \
                             tr_z_yyy_xxyyzz,   \
                             tr_z_yyy_xxyzzz,   \
                             tr_z_yyy_xxzzzz,   \
                             tr_z_yyy_xyyyyy,   \
                             tr_z_yyy_xyyyyz,   \
                             tr_z_yyy_xyyyzz,   \
                             tr_z_yyy_xyyzzz,   \
                             tr_z_yyy_xyzzzz,   \
                             tr_z_yyy_xzzzzz,   \
                             tr_z_yyy_yyyyyy,   \
                             tr_z_yyy_yyyyyz,   \
                             tr_z_yyy_yyyyzz,   \
                             tr_z_yyy_yyyzzz,   \
                             tr_z_yyy_yyzzzz,   \
                             tr_z_yyy_yzzzzz,   \
                             tr_z_yyy_zzzzzz,   \
                             tr_z_yyyy_xxxxx,   \
                             tr_z_yyyy_xxxxxx,  \
                             tr_z_yyyy_xxxxxy,  \
                             tr_z_yyyy_xxxxxz,  \
                             tr_z_yyyy_xxxxy,   \
                             tr_z_yyyy_xxxxyy,  \
                             tr_z_yyyy_xxxxyz,  \
                             tr_z_yyyy_xxxxz,   \
                             tr_z_yyyy_xxxxzz,  \
                             tr_z_yyyy_xxxyy,   \
                             tr_z_yyyy_xxxyyy,  \
                             tr_z_yyyy_xxxyyz,  \
                             tr_z_yyyy_xxxyz,   \
                             tr_z_yyyy_xxxyzz,  \
                             tr_z_yyyy_xxxzz,   \
                             tr_z_yyyy_xxxzzz,  \
                             tr_z_yyyy_xxyyy,   \
                             tr_z_yyyy_xxyyyy,  \
                             tr_z_yyyy_xxyyyz,  \
                             tr_z_yyyy_xxyyz,   \
                             tr_z_yyyy_xxyyzz,  \
                             tr_z_yyyy_xxyzz,   \
                             tr_z_yyyy_xxyzzz,  \
                             tr_z_yyyy_xxzzz,   \
                             tr_z_yyyy_xxzzzz,  \
                             tr_z_yyyy_xyyyy,   \
                             tr_z_yyyy_xyyyyy,  \
                             tr_z_yyyy_xyyyyz,  \
                             tr_z_yyyy_xyyyz,   \
                             tr_z_yyyy_xyyyzz,  \
                             tr_z_yyyy_xyyzz,   \
                             tr_z_yyyy_xyyzzz,  \
                             tr_z_yyyy_xyzzz,   \
                             tr_z_yyyy_xyzzzz,  \
                             tr_z_yyyy_xzzzz,   \
                             tr_z_yyyy_xzzzzz,  \
                             tr_z_yyyy_yyyyy,   \
                             tr_z_yyyy_yyyyyy,  \
                             tr_z_yyyy_yyyyyz,  \
                             tr_z_yyyy_yyyyz,   \
                             tr_z_yyyy_yyyyzz,  \
                             tr_z_yyyy_yyyzz,   \
                             tr_z_yyyy_yyyzzz,  \
                             tr_z_yyyy_yyzzz,   \
                             tr_z_yyyy_yyzzzz,  \
                             tr_z_yyyy_yzzzz,   \
                             tr_z_yyyy_yzzzzz,  \
                             tr_z_yyyy_zzzzz,   \
                             tr_z_yyyy_zzzzzz,  \
                             tr_z_yyyyy_xxxxxx, \
                             tr_z_yyyyy_xxxxxy, \
                             tr_z_yyyyy_xxxxxz, \
                             tr_z_yyyyy_xxxxyy, \
                             tr_z_yyyyy_xxxxyz, \
                             tr_z_yyyyy_xxxxzz, \
                             tr_z_yyyyy_xxxyyy, \
                             tr_z_yyyyy_xxxyyz, \
                             tr_z_yyyyy_xxxyzz, \
                             tr_z_yyyyy_xxxzzz, \
                             tr_z_yyyyy_xxyyyy, \
                             tr_z_yyyyy_xxyyyz, \
                             tr_z_yyyyy_xxyyzz, \
                             tr_z_yyyyy_xxyzzz, \
                             tr_z_yyyyy_xxzzzz, \
                             tr_z_yyyyy_xyyyyy, \
                             tr_z_yyyyy_xyyyyz, \
                             tr_z_yyyyy_xyyyzz, \
                             tr_z_yyyyy_xyyzzz, \
                             tr_z_yyyyy_xyzzzz, \
                             tr_z_yyyyy_xzzzzz, \
                             tr_z_yyyyy_yyyyyy, \
                             tr_z_yyyyy_yyyyyz, \
                             tr_z_yyyyy_yyyyzz, \
                             tr_z_yyyyy_yyyzzz, \
                             tr_z_yyyyy_yyzzzz, \
                             tr_z_yyyyy_yzzzzz, \
                             tr_z_yyyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyy_xxxxxx[i] = 4.0 * tr_z_yyy_xxxxxx[i] * fe_0 + tr_z_yyyy_xxxxxx[i] * pa_y[i];

        tr_z_yyyyy_xxxxxy[i] = 4.0 * tr_z_yyy_xxxxxy[i] * fe_0 + tr_z_yyyy_xxxxx[i] * fe_0 + tr_z_yyyy_xxxxxy[i] * pa_y[i];

        tr_z_yyyyy_xxxxxz[i] = 4.0 * tr_z_yyy_xxxxxz[i] * fe_0 + tr_z_yyyy_xxxxxz[i] * pa_y[i];

        tr_z_yyyyy_xxxxyy[i] = 4.0 * tr_z_yyy_xxxxyy[i] * fe_0 + 2.0 * tr_z_yyyy_xxxxy[i] * fe_0 + tr_z_yyyy_xxxxyy[i] * pa_y[i];

        tr_z_yyyyy_xxxxyz[i] = 4.0 * tr_z_yyy_xxxxyz[i] * fe_0 + tr_z_yyyy_xxxxz[i] * fe_0 + tr_z_yyyy_xxxxyz[i] * pa_y[i];

        tr_z_yyyyy_xxxxzz[i] = 4.0 * tr_z_yyy_xxxxzz[i] * fe_0 + tr_z_yyyy_xxxxzz[i] * pa_y[i];

        tr_z_yyyyy_xxxyyy[i] = 4.0 * tr_z_yyy_xxxyyy[i] * fe_0 + 3.0 * tr_z_yyyy_xxxyy[i] * fe_0 + tr_z_yyyy_xxxyyy[i] * pa_y[i];

        tr_z_yyyyy_xxxyyz[i] = 4.0 * tr_z_yyy_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyyy_xxxyz[i] * fe_0 + tr_z_yyyy_xxxyyz[i] * pa_y[i];

        tr_z_yyyyy_xxxyzz[i] = 4.0 * tr_z_yyy_xxxyzz[i] * fe_0 + tr_z_yyyy_xxxzz[i] * fe_0 + tr_z_yyyy_xxxyzz[i] * pa_y[i];

        tr_z_yyyyy_xxxzzz[i] = 4.0 * tr_z_yyy_xxxzzz[i] * fe_0 + tr_z_yyyy_xxxzzz[i] * pa_y[i];

        tr_z_yyyyy_xxyyyy[i] = 4.0 * tr_z_yyy_xxyyyy[i] * fe_0 + 4.0 * tr_z_yyyy_xxyyy[i] * fe_0 + tr_z_yyyy_xxyyyy[i] * pa_y[i];

        tr_z_yyyyy_xxyyyz[i] = 4.0 * tr_z_yyy_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyyy_xxyyz[i] * fe_0 + tr_z_yyyy_xxyyyz[i] * pa_y[i];

        tr_z_yyyyy_xxyyzz[i] = 4.0 * tr_z_yyy_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyyy_xxyzz[i] * fe_0 + tr_z_yyyy_xxyyzz[i] * pa_y[i];

        tr_z_yyyyy_xxyzzz[i] = 4.0 * tr_z_yyy_xxyzzz[i] * fe_0 + tr_z_yyyy_xxzzz[i] * fe_0 + tr_z_yyyy_xxyzzz[i] * pa_y[i];

        tr_z_yyyyy_xxzzzz[i] = 4.0 * tr_z_yyy_xxzzzz[i] * fe_0 + tr_z_yyyy_xxzzzz[i] * pa_y[i];

        tr_z_yyyyy_xyyyyy[i] = 4.0 * tr_z_yyy_xyyyyy[i] * fe_0 + 5.0 * tr_z_yyyy_xyyyy[i] * fe_0 + tr_z_yyyy_xyyyyy[i] * pa_y[i];

        tr_z_yyyyy_xyyyyz[i] = 4.0 * tr_z_yyy_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyyy_xyyyz[i] * fe_0 + tr_z_yyyy_xyyyyz[i] * pa_y[i];

        tr_z_yyyyy_xyyyzz[i] = 4.0 * tr_z_yyy_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyyy_xyyzz[i] * fe_0 + tr_z_yyyy_xyyyzz[i] * pa_y[i];

        tr_z_yyyyy_xyyzzz[i] = 4.0 * tr_z_yyy_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyyy_xyzzz[i] * fe_0 + tr_z_yyyy_xyyzzz[i] * pa_y[i];

        tr_z_yyyyy_xyzzzz[i] = 4.0 * tr_z_yyy_xyzzzz[i] * fe_0 + tr_z_yyyy_xzzzz[i] * fe_0 + tr_z_yyyy_xyzzzz[i] * pa_y[i];

        tr_z_yyyyy_xzzzzz[i] = 4.0 * tr_z_yyy_xzzzzz[i] * fe_0 + tr_z_yyyy_xzzzzz[i] * pa_y[i];

        tr_z_yyyyy_yyyyyy[i] = 4.0 * tr_z_yyy_yyyyyy[i] * fe_0 + 6.0 * tr_z_yyyy_yyyyy[i] * fe_0 + tr_z_yyyy_yyyyyy[i] * pa_y[i];

        tr_z_yyyyy_yyyyyz[i] = 4.0 * tr_z_yyy_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyyy_yyyyz[i] * fe_0 + tr_z_yyyy_yyyyyz[i] * pa_y[i];

        tr_z_yyyyy_yyyyzz[i] = 4.0 * tr_z_yyy_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyyy_yyyzz[i] * fe_0 + tr_z_yyyy_yyyyzz[i] * pa_y[i];

        tr_z_yyyyy_yyyzzz[i] = 4.0 * tr_z_yyy_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyyy_yyzzz[i] * fe_0 + tr_z_yyyy_yyyzzz[i] * pa_y[i];

        tr_z_yyyyy_yyzzzz[i] = 4.0 * tr_z_yyy_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyyy_yzzzz[i] * fe_0 + tr_z_yyyy_yyzzzz[i] * pa_y[i];

        tr_z_yyyyy_yzzzzz[i] = 4.0 * tr_z_yyy_yzzzzz[i] * fe_0 + tr_z_yyyy_zzzzz[i] * fe_0 + tr_z_yyyy_yzzzzz[i] * pa_y[i];

        tr_z_yyyyy_zzzzzz[i] = 4.0 * tr_z_yyy_zzzzzz[i] * fe_0 + tr_z_yyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 1624-1652 components of targeted buffer : HI

    auto tr_z_yyyyz_xxxxxx = pbuffer.data(idx_dip_hi + 1624);

    auto tr_z_yyyyz_xxxxxy = pbuffer.data(idx_dip_hi + 1625);

    auto tr_z_yyyyz_xxxxxz = pbuffer.data(idx_dip_hi + 1626);

    auto tr_z_yyyyz_xxxxyy = pbuffer.data(idx_dip_hi + 1627);

    auto tr_z_yyyyz_xxxxyz = pbuffer.data(idx_dip_hi + 1628);

    auto tr_z_yyyyz_xxxxzz = pbuffer.data(idx_dip_hi + 1629);

    auto tr_z_yyyyz_xxxyyy = pbuffer.data(idx_dip_hi + 1630);

    auto tr_z_yyyyz_xxxyyz = pbuffer.data(idx_dip_hi + 1631);

    auto tr_z_yyyyz_xxxyzz = pbuffer.data(idx_dip_hi + 1632);

    auto tr_z_yyyyz_xxxzzz = pbuffer.data(idx_dip_hi + 1633);

    auto tr_z_yyyyz_xxyyyy = pbuffer.data(idx_dip_hi + 1634);

    auto tr_z_yyyyz_xxyyyz = pbuffer.data(idx_dip_hi + 1635);

    auto tr_z_yyyyz_xxyyzz = pbuffer.data(idx_dip_hi + 1636);

    auto tr_z_yyyyz_xxyzzz = pbuffer.data(idx_dip_hi + 1637);

    auto tr_z_yyyyz_xxzzzz = pbuffer.data(idx_dip_hi + 1638);

    auto tr_z_yyyyz_xyyyyy = pbuffer.data(idx_dip_hi + 1639);

    auto tr_z_yyyyz_xyyyyz = pbuffer.data(idx_dip_hi + 1640);

    auto tr_z_yyyyz_xyyyzz = pbuffer.data(idx_dip_hi + 1641);

    auto tr_z_yyyyz_xyyzzz = pbuffer.data(idx_dip_hi + 1642);

    auto tr_z_yyyyz_xyzzzz = pbuffer.data(idx_dip_hi + 1643);

    auto tr_z_yyyyz_xzzzzz = pbuffer.data(idx_dip_hi + 1644);

    auto tr_z_yyyyz_yyyyyy = pbuffer.data(idx_dip_hi + 1645);

    auto tr_z_yyyyz_yyyyyz = pbuffer.data(idx_dip_hi + 1646);

    auto tr_z_yyyyz_yyyyzz = pbuffer.data(idx_dip_hi + 1647);

    auto tr_z_yyyyz_yyyzzz = pbuffer.data(idx_dip_hi + 1648);

    auto tr_z_yyyyz_yyzzzz = pbuffer.data(idx_dip_hi + 1649);

    auto tr_z_yyyyz_yzzzzz = pbuffer.data(idx_dip_hi + 1650);

    auto tr_z_yyyyz_zzzzzz = pbuffer.data(idx_dip_hi + 1651);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             tr_z_yyyy_xxxxxy,  \
                             tr_z_yyyy_xxxxyy,  \
                             tr_z_yyyy_xxxyyy,  \
                             tr_z_yyyy_xxyyyy,  \
                             tr_z_yyyy_xyyyyy,  \
                             tr_z_yyyy_yyyyyy,  \
                             tr_z_yyyyz_xxxxxx, \
                             tr_z_yyyyz_xxxxxy, \
                             tr_z_yyyyz_xxxxxz, \
                             tr_z_yyyyz_xxxxyy, \
                             tr_z_yyyyz_xxxxyz, \
                             tr_z_yyyyz_xxxxzz, \
                             tr_z_yyyyz_xxxyyy, \
                             tr_z_yyyyz_xxxyyz, \
                             tr_z_yyyyz_xxxyzz, \
                             tr_z_yyyyz_xxxzzz, \
                             tr_z_yyyyz_xxyyyy, \
                             tr_z_yyyyz_xxyyyz, \
                             tr_z_yyyyz_xxyyzz, \
                             tr_z_yyyyz_xxyzzz, \
                             tr_z_yyyyz_xxzzzz, \
                             tr_z_yyyyz_xyyyyy, \
                             tr_z_yyyyz_xyyyyz, \
                             tr_z_yyyyz_xyyyzz, \
                             tr_z_yyyyz_xyyzzz, \
                             tr_z_yyyyz_xyzzzz, \
                             tr_z_yyyyz_xzzzzz, \
                             tr_z_yyyyz_yyyyyy, \
                             tr_z_yyyyz_yyyyyz, \
                             tr_z_yyyyz_yyyyzz, \
                             tr_z_yyyyz_yyyzzz, \
                             tr_z_yyyyz_yyzzzz, \
                             tr_z_yyyyz_yzzzzz, \
                             tr_z_yyyyz_zzzzzz, \
                             tr_z_yyyz_xxxxxx,  \
                             tr_z_yyyz_xxxxxz,  \
                             tr_z_yyyz_xxxxyz,  \
                             tr_z_yyyz_xxxxz,   \
                             tr_z_yyyz_xxxxzz,  \
                             tr_z_yyyz_xxxyyz,  \
                             tr_z_yyyz_xxxyz,   \
                             tr_z_yyyz_xxxyzz,  \
                             tr_z_yyyz_xxxzz,   \
                             tr_z_yyyz_xxxzzz,  \
                             tr_z_yyyz_xxyyyz,  \
                             tr_z_yyyz_xxyyz,   \
                             tr_z_yyyz_xxyyzz,  \
                             tr_z_yyyz_xxyzz,   \
                             tr_z_yyyz_xxyzzz,  \
                             tr_z_yyyz_xxzzz,   \
                             tr_z_yyyz_xxzzzz,  \
                             tr_z_yyyz_xyyyyz,  \
                             tr_z_yyyz_xyyyz,   \
                             tr_z_yyyz_xyyyzz,  \
                             tr_z_yyyz_xyyzz,   \
                             tr_z_yyyz_xyyzzz,  \
                             tr_z_yyyz_xyzzz,   \
                             tr_z_yyyz_xyzzzz,  \
                             tr_z_yyyz_xzzzz,   \
                             tr_z_yyyz_xzzzzz,  \
                             tr_z_yyyz_yyyyyz,  \
                             tr_z_yyyz_yyyyz,   \
                             tr_z_yyyz_yyyyzz,  \
                             tr_z_yyyz_yyyzz,   \
                             tr_z_yyyz_yyyzzz,  \
                             tr_z_yyyz_yyzzz,   \
                             tr_z_yyyz_yyzzzz,  \
                             tr_z_yyyz_yzzzz,   \
                             tr_z_yyyz_yzzzzz,  \
                             tr_z_yyyz_zzzzz,   \
                             tr_z_yyyz_zzzzzz,  \
                             tr_z_yyz_xxxxxx,   \
                             tr_z_yyz_xxxxxz,   \
                             tr_z_yyz_xxxxyz,   \
                             tr_z_yyz_xxxxzz,   \
                             tr_z_yyz_xxxyyz,   \
                             tr_z_yyz_xxxyzz,   \
                             tr_z_yyz_xxxzzz,   \
                             tr_z_yyz_xxyyyz,   \
                             tr_z_yyz_xxyyzz,   \
                             tr_z_yyz_xxyzzz,   \
                             tr_z_yyz_xxzzzz,   \
                             tr_z_yyz_xyyyyz,   \
                             tr_z_yyz_xyyyzz,   \
                             tr_z_yyz_xyyzzz,   \
                             tr_z_yyz_xyzzzz,   \
                             tr_z_yyz_xzzzzz,   \
                             tr_z_yyz_yyyyyz,   \
                             tr_z_yyz_yyyyzz,   \
                             tr_z_yyz_yyyzzz,   \
                             tr_z_yyz_yyzzzz,   \
                             tr_z_yyz_yzzzzz,   \
                             tr_z_yyz_zzzzzz,   \
                             ts_yyyy_xxxxxy,    \
                             ts_yyyy_xxxxyy,    \
                             ts_yyyy_xxxyyy,    \
                             ts_yyyy_xxyyyy,    \
                             ts_yyyy_xyyyyy,    \
                             ts_yyyy_yyyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyz_xxxxxx[i] = 3.0 * tr_z_yyz_xxxxxx[i] * fe_0 + tr_z_yyyz_xxxxxx[i] * pa_y[i];

        tr_z_yyyyz_xxxxxy[i] = ts_yyyy_xxxxxy[i] * fe_0 + tr_z_yyyy_xxxxxy[i] * pa_z[i];

        tr_z_yyyyz_xxxxxz[i] = 3.0 * tr_z_yyz_xxxxxz[i] * fe_0 + tr_z_yyyz_xxxxxz[i] * pa_y[i];

        tr_z_yyyyz_xxxxyy[i] = ts_yyyy_xxxxyy[i] * fe_0 + tr_z_yyyy_xxxxyy[i] * pa_z[i];

        tr_z_yyyyz_xxxxyz[i] = 3.0 * tr_z_yyz_xxxxyz[i] * fe_0 + tr_z_yyyz_xxxxz[i] * fe_0 + tr_z_yyyz_xxxxyz[i] * pa_y[i];

        tr_z_yyyyz_xxxxzz[i] = 3.0 * tr_z_yyz_xxxxzz[i] * fe_0 + tr_z_yyyz_xxxxzz[i] * pa_y[i];

        tr_z_yyyyz_xxxyyy[i] = ts_yyyy_xxxyyy[i] * fe_0 + tr_z_yyyy_xxxyyy[i] * pa_z[i];

        tr_z_yyyyz_xxxyyz[i] = 3.0 * tr_z_yyz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyyz_xxxyz[i] * fe_0 + tr_z_yyyz_xxxyyz[i] * pa_y[i];

        tr_z_yyyyz_xxxyzz[i] = 3.0 * tr_z_yyz_xxxyzz[i] * fe_0 + tr_z_yyyz_xxxzz[i] * fe_0 + tr_z_yyyz_xxxyzz[i] * pa_y[i];

        tr_z_yyyyz_xxxzzz[i] = 3.0 * tr_z_yyz_xxxzzz[i] * fe_0 + tr_z_yyyz_xxxzzz[i] * pa_y[i];

        tr_z_yyyyz_xxyyyy[i] = ts_yyyy_xxyyyy[i] * fe_0 + tr_z_yyyy_xxyyyy[i] * pa_z[i];

        tr_z_yyyyz_xxyyyz[i] = 3.0 * tr_z_yyz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyyz_xxyyz[i] * fe_0 + tr_z_yyyz_xxyyyz[i] * pa_y[i];

        tr_z_yyyyz_xxyyzz[i] = 3.0 * tr_z_yyz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyyz_xxyzz[i] * fe_0 + tr_z_yyyz_xxyyzz[i] * pa_y[i];

        tr_z_yyyyz_xxyzzz[i] = 3.0 * tr_z_yyz_xxyzzz[i] * fe_0 + tr_z_yyyz_xxzzz[i] * fe_0 + tr_z_yyyz_xxyzzz[i] * pa_y[i];

        tr_z_yyyyz_xxzzzz[i] = 3.0 * tr_z_yyz_xxzzzz[i] * fe_0 + tr_z_yyyz_xxzzzz[i] * pa_y[i];

        tr_z_yyyyz_xyyyyy[i] = ts_yyyy_xyyyyy[i] * fe_0 + tr_z_yyyy_xyyyyy[i] * pa_z[i];

        tr_z_yyyyz_xyyyyz[i] = 3.0 * tr_z_yyz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyyz_xyyyz[i] * fe_0 + tr_z_yyyz_xyyyyz[i] * pa_y[i];

        tr_z_yyyyz_xyyyzz[i] = 3.0 * tr_z_yyz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyyz_xyyzz[i] * fe_0 + tr_z_yyyz_xyyyzz[i] * pa_y[i];

        tr_z_yyyyz_xyyzzz[i] = 3.0 * tr_z_yyz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyyz_xyzzz[i] * fe_0 + tr_z_yyyz_xyyzzz[i] * pa_y[i];

        tr_z_yyyyz_xyzzzz[i] = 3.0 * tr_z_yyz_xyzzzz[i] * fe_0 + tr_z_yyyz_xzzzz[i] * fe_0 + tr_z_yyyz_xyzzzz[i] * pa_y[i];

        tr_z_yyyyz_xzzzzz[i] = 3.0 * tr_z_yyz_xzzzzz[i] * fe_0 + tr_z_yyyz_xzzzzz[i] * pa_y[i];

        tr_z_yyyyz_yyyyyy[i] = ts_yyyy_yyyyyy[i] * fe_0 + tr_z_yyyy_yyyyyy[i] * pa_z[i];

        tr_z_yyyyz_yyyyyz[i] = 3.0 * tr_z_yyz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyyz_yyyyz[i] * fe_0 + tr_z_yyyz_yyyyyz[i] * pa_y[i];

        tr_z_yyyyz_yyyyzz[i] = 3.0 * tr_z_yyz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyyz_yyyzz[i] * fe_0 + tr_z_yyyz_yyyyzz[i] * pa_y[i];

        tr_z_yyyyz_yyyzzz[i] = 3.0 * tr_z_yyz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyyz_yyzzz[i] * fe_0 + tr_z_yyyz_yyyzzz[i] * pa_y[i];

        tr_z_yyyyz_yyzzzz[i] = 3.0 * tr_z_yyz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyyz_yzzzz[i] * fe_0 + tr_z_yyyz_yyzzzz[i] * pa_y[i];

        tr_z_yyyyz_yzzzzz[i] = 3.0 * tr_z_yyz_yzzzzz[i] * fe_0 + tr_z_yyyz_zzzzz[i] * fe_0 + tr_z_yyyz_yzzzzz[i] * pa_y[i];

        tr_z_yyyyz_zzzzzz[i] = 3.0 * tr_z_yyz_zzzzzz[i] * fe_0 + tr_z_yyyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1652-1680 components of targeted buffer : HI

    auto tr_z_yyyzz_xxxxxx = pbuffer.data(idx_dip_hi + 1652);

    auto tr_z_yyyzz_xxxxxy = pbuffer.data(idx_dip_hi + 1653);

    auto tr_z_yyyzz_xxxxxz = pbuffer.data(idx_dip_hi + 1654);

    auto tr_z_yyyzz_xxxxyy = pbuffer.data(idx_dip_hi + 1655);

    auto tr_z_yyyzz_xxxxyz = pbuffer.data(idx_dip_hi + 1656);

    auto tr_z_yyyzz_xxxxzz = pbuffer.data(idx_dip_hi + 1657);

    auto tr_z_yyyzz_xxxyyy = pbuffer.data(idx_dip_hi + 1658);

    auto tr_z_yyyzz_xxxyyz = pbuffer.data(idx_dip_hi + 1659);

    auto tr_z_yyyzz_xxxyzz = pbuffer.data(idx_dip_hi + 1660);

    auto tr_z_yyyzz_xxxzzz = pbuffer.data(idx_dip_hi + 1661);

    auto tr_z_yyyzz_xxyyyy = pbuffer.data(idx_dip_hi + 1662);

    auto tr_z_yyyzz_xxyyyz = pbuffer.data(idx_dip_hi + 1663);

    auto tr_z_yyyzz_xxyyzz = pbuffer.data(idx_dip_hi + 1664);

    auto tr_z_yyyzz_xxyzzz = pbuffer.data(idx_dip_hi + 1665);

    auto tr_z_yyyzz_xxzzzz = pbuffer.data(idx_dip_hi + 1666);

    auto tr_z_yyyzz_xyyyyy = pbuffer.data(idx_dip_hi + 1667);

    auto tr_z_yyyzz_xyyyyz = pbuffer.data(idx_dip_hi + 1668);

    auto tr_z_yyyzz_xyyyzz = pbuffer.data(idx_dip_hi + 1669);

    auto tr_z_yyyzz_xyyzzz = pbuffer.data(idx_dip_hi + 1670);

    auto tr_z_yyyzz_xyzzzz = pbuffer.data(idx_dip_hi + 1671);

    auto tr_z_yyyzz_xzzzzz = pbuffer.data(idx_dip_hi + 1672);

    auto tr_z_yyyzz_yyyyyy = pbuffer.data(idx_dip_hi + 1673);

    auto tr_z_yyyzz_yyyyyz = pbuffer.data(idx_dip_hi + 1674);

    auto tr_z_yyyzz_yyyyzz = pbuffer.data(idx_dip_hi + 1675);

    auto tr_z_yyyzz_yyyzzz = pbuffer.data(idx_dip_hi + 1676);

    auto tr_z_yyyzz_yyzzzz = pbuffer.data(idx_dip_hi + 1677);

    auto tr_z_yyyzz_yzzzzz = pbuffer.data(idx_dip_hi + 1678);

    auto tr_z_yyyzz_zzzzzz = pbuffer.data(idx_dip_hi + 1679);

#pragma omp simd aligned(pa_y,                  \
                             tr_z_yyyzz_xxxxxx, \
                             tr_z_yyyzz_xxxxxy, \
                             tr_z_yyyzz_xxxxxz, \
                             tr_z_yyyzz_xxxxyy, \
                             tr_z_yyyzz_xxxxyz, \
                             tr_z_yyyzz_xxxxzz, \
                             tr_z_yyyzz_xxxyyy, \
                             tr_z_yyyzz_xxxyyz, \
                             tr_z_yyyzz_xxxyzz, \
                             tr_z_yyyzz_xxxzzz, \
                             tr_z_yyyzz_xxyyyy, \
                             tr_z_yyyzz_xxyyyz, \
                             tr_z_yyyzz_xxyyzz, \
                             tr_z_yyyzz_xxyzzz, \
                             tr_z_yyyzz_xxzzzz, \
                             tr_z_yyyzz_xyyyyy, \
                             tr_z_yyyzz_xyyyyz, \
                             tr_z_yyyzz_xyyyzz, \
                             tr_z_yyyzz_xyyzzz, \
                             tr_z_yyyzz_xyzzzz, \
                             tr_z_yyyzz_xzzzzz, \
                             tr_z_yyyzz_yyyyyy, \
                             tr_z_yyyzz_yyyyyz, \
                             tr_z_yyyzz_yyyyzz, \
                             tr_z_yyyzz_yyyzzz, \
                             tr_z_yyyzz_yyzzzz, \
                             tr_z_yyyzz_yzzzzz, \
                             tr_z_yyyzz_zzzzzz, \
                             tr_z_yyzz_xxxxx,   \
                             tr_z_yyzz_xxxxxx,  \
                             tr_z_yyzz_xxxxxy,  \
                             tr_z_yyzz_xxxxxz,  \
                             tr_z_yyzz_xxxxy,   \
                             tr_z_yyzz_xxxxyy,  \
                             tr_z_yyzz_xxxxyz,  \
                             tr_z_yyzz_xxxxz,   \
                             tr_z_yyzz_xxxxzz,  \
                             tr_z_yyzz_xxxyy,   \
                             tr_z_yyzz_xxxyyy,  \
                             tr_z_yyzz_xxxyyz,  \
                             tr_z_yyzz_xxxyz,   \
                             tr_z_yyzz_xxxyzz,  \
                             tr_z_yyzz_xxxzz,   \
                             tr_z_yyzz_xxxzzz,  \
                             tr_z_yyzz_xxyyy,   \
                             tr_z_yyzz_xxyyyy,  \
                             tr_z_yyzz_xxyyyz,  \
                             tr_z_yyzz_xxyyz,   \
                             tr_z_yyzz_xxyyzz,  \
                             tr_z_yyzz_xxyzz,   \
                             tr_z_yyzz_xxyzzz,  \
                             tr_z_yyzz_xxzzz,   \
                             tr_z_yyzz_xxzzzz,  \
                             tr_z_yyzz_xyyyy,   \
                             tr_z_yyzz_xyyyyy,  \
                             tr_z_yyzz_xyyyyz,  \
                             tr_z_yyzz_xyyyz,   \
                             tr_z_yyzz_xyyyzz,  \
                             tr_z_yyzz_xyyzz,   \
                             tr_z_yyzz_xyyzzz,  \
                             tr_z_yyzz_xyzzz,   \
                             tr_z_yyzz_xyzzzz,  \
                             tr_z_yyzz_xzzzz,   \
                             tr_z_yyzz_xzzzzz,  \
                             tr_z_yyzz_yyyyy,   \
                             tr_z_yyzz_yyyyyy,  \
                             tr_z_yyzz_yyyyyz,  \
                             tr_z_yyzz_yyyyz,   \
                             tr_z_yyzz_yyyyzz,  \
                             tr_z_yyzz_yyyzz,   \
                             tr_z_yyzz_yyyzzz,  \
                             tr_z_yyzz_yyzzz,   \
                             tr_z_yyzz_yyzzzz,  \
                             tr_z_yyzz_yzzzz,   \
                             tr_z_yyzz_yzzzzz,  \
                             tr_z_yyzz_zzzzz,   \
                             tr_z_yyzz_zzzzzz,  \
                             tr_z_yzz_xxxxxx,   \
                             tr_z_yzz_xxxxxy,   \
                             tr_z_yzz_xxxxxz,   \
                             tr_z_yzz_xxxxyy,   \
                             tr_z_yzz_xxxxyz,   \
                             tr_z_yzz_xxxxzz,   \
                             tr_z_yzz_xxxyyy,   \
                             tr_z_yzz_xxxyyz,   \
                             tr_z_yzz_xxxyzz,   \
                             tr_z_yzz_xxxzzz,   \
                             tr_z_yzz_xxyyyy,   \
                             tr_z_yzz_xxyyyz,   \
                             tr_z_yzz_xxyyzz,   \
                             tr_z_yzz_xxyzzz,   \
                             tr_z_yzz_xxzzzz,   \
                             tr_z_yzz_xyyyyy,   \
                             tr_z_yzz_xyyyyz,   \
                             tr_z_yzz_xyyyzz,   \
                             tr_z_yzz_xyyzzz,   \
                             tr_z_yzz_xyzzzz,   \
                             tr_z_yzz_xzzzzz,   \
                             tr_z_yzz_yyyyyy,   \
                             tr_z_yzz_yyyyyz,   \
                             tr_z_yzz_yyyyzz,   \
                             tr_z_yzz_yyyzzz,   \
                             tr_z_yzz_yyzzzz,   \
                             tr_z_yzz_yzzzzz,   \
                             tr_z_yzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzz_xxxxxx[i] = 2.0 * tr_z_yzz_xxxxxx[i] * fe_0 + tr_z_yyzz_xxxxxx[i] * pa_y[i];

        tr_z_yyyzz_xxxxxy[i] = 2.0 * tr_z_yzz_xxxxxy[i] * fe_0 + tr_z_yyzz_xxxxx[i] * fe_0 + tr_z_yyzz_xxxxxy[i] * pa_y[i];

        tr_z_yyyzz_xxxxxz[i] = 2.0 * tr_z_yzz_xxxxxz[i] * fe_0 + tr_z_yyzz_xxxxxz[i] * pa_y[i];

        tr_z_yyyzz_xxxxyy[i] = 2.0 * tr_z_yzz_xxxxyy[i] * fe_0 + 2.0 * tr_z_yyzz_xxxxy[i] * fe_0 + tr_z_yyzz_xxxxyy[i] * pa_y[i];

        tr_z_yyyzz_xxxxyz[i] = 2.0 * tr_z_yzz_xxxxyz[i] * fe_0 + tr_z_yyzz_xxxxz[i] * fe_0 + tr_z_yyzz_xxxxyz[i] * pa_y[i];

        tr_z_yyyzz_xxxxzz[i] = 2.0 * tr_z_yzz_xxxxzz[i] * fe_0 + tr_z_yyzz_xxxxzz[i] * pa_y[i];

        tr_z_yyyzz_xxxyyy[i] = 2.0 * tr_z_yzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_yyzz_xxxyy[i] * fe_0 + tr_z_yyzz_xxxyyy[i] * pa_y[i];

        tr_z_yyyzz_xxxyyz[i] = 2.0 * tr_z_yzz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yyzz_xxxyz[i] * fe_0 + tr_z_yyzz_xxxyyz[i] * pa_y[i];

        tr_z_yyyzz_xxxyzz[i] = 2.0 * tr_z_yzz_xxxyzz[i] * fe_0 + tr_z_yyzz_xxxzz[i] * fe_0 + tr_z_yyzz_xxxyzz[i] * pa_y[i];

        tr_z_yyyzz_xxxzzz[i] = 2.0 * tr_z_yzz_xxxzzz[i] * fe_0 + tr_z_yyzz_xxxzzz[i] * pa_y[i];

        tr_z_yyyzz_xxyyyy[i] = 2.0 * tr_z_yzz_xxyyyy[i] * fe_0 + 4.0 * tr_z_yyzz_xxyyy[i] * fe_0 + tr_z_yyzz_xxyyyy[i] * pa_y[i];

        tr_z_yyyzz_xxyyyz[i] = 2.0 * tr_z_yzz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yyzz_xxyyz[i] * fe_0 + tr_z_yyzz_xxyyyz[i] * pa_y[i];

        tr_z_yyyzz_xxyyzz[i] = 2.0 * tr_z_yzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yyzz_xxyzz[i] * fe_0 + tr_z_yyzz_xxyyzz[i] * pa_y[i];

        tr_z_yyyzz_xxyzzz[i] = 2.0 * tr_z_yzz_xxyzzz[i] * fe_0 + tr_z_yyzz_xxzzz[i] * fe_0 + tr_z_yyzz_xxyzzz[i] * pa_y[i];

        tr_z_yyyzz_xxzzzz[i] = 2.0 * tr_z_yzz_xxzzzz[i] * fe_0 + tr_z_yyzz_xxzzzz[i] * pa_y[i];

        tr_z_yyyzz_xyyyyy[i] = 2.0 * tr_z_yzz_xyyyyy[i] * fe_0 + 5.0 * tr_z_yyzz_xyyyy[i] * fe_0 + tr_z_yyzz_xyyyyy[i] * pa_y[i];

        tr_z_yyyzz_xyyyyz[i] = 2.0 * tr_z_yzz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yyzz_xyyyz[i] * fe_0 + tr_z_yyzz_xyyyyz[i] * pa_y[i];

        tr_z_yyyzz_xyyyzz[i] = 2.0 * tr_z_yzz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yyzz_xyyzz[i] * fe_0 + tr_z_yyzz_xyyyzz[i] * pa_y[i];

        tr_z_yyyzz_xyyzzz[i] = 2.0 * tr_z_yzz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yyzz_xyzzz[i] * fe_0 + tr_z_yyzz_xyyzzz[i] * pa_y[i];

        tr_z_yyyzz_xyzzzz[i] = 2.0 * tr_z_yzz_xyzzzz[i] * fe_0 + tr_z_yyzz_xzzzz[i] * fe_0 + tr_z_yyzz_xyzzzz[i] * pa_y[i];

        tr_z_yyyzz_xzzzzz[i] = 2.0 * tr_z_yzz_xzzzzz[i] * fe_0 + tr_z_yyzz_xzzzzz[i] * pa_y[i];

        tr_z_yyyzz_yyyyyy[i] = 2.0 * tr_z_yzz_yyyyyy[i] * fe_0 + 6.0 * tr_z_yyzz_yyyyy[i] * fe_0 + tr_z_yyzz_yyyyyy[i] * pa_y[i];

        tr_z_yyyzz_yyyyyz[i] = 2.0 * tr_z_yzz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yyzz_yyyyz[i] * fe_0 + tr_z_yyzz_yyyyyz[i] * pa_y[i];

        tr_z_yyyzz_yyyyzz[i] = 2.0 * tr_z_yzz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yyzz_yyyzz[i] * fe_0 + tr_z_yyzz_yyyyzz[i] * pa_y[i];

        tr_z_yyyzz_yyyzzz[i] = 2.0 * tr_z_yzz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yyzz_yyzzz[i] * fe_0 + tr_z_yyzz_yyyzzz[i] * pa_y[i];

        tr_z_yyyzz_yyzzzz[i] = 2.0 * tr_z_yzz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yyzz_yzzzz[i] * fe_0 + tr_z_yyzz_yyzzzz[i] * pa_y[i];

        tr_z_yyyzz_yzzzzz[i] = 2.0 * tr_z_yzz_yzzzzz[i] * fe_0 + tr_z_yyzz_zzzzz[i] * fe_0 + tr_z_yyzz_yzzzzz[i] * pa_y[i];

        tr_z_yyyzz_zzzzzz[i] = 2.0 * tr_z_yzz_zzzzzz[i] * fe_0 + tr_z_yyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1680-1708 components of targeted buffer : HI

    auto tr_z_yyzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1680);

    auto tr_z_yyzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1681);

    auto tr_z_yyzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1682);

    auto tr_z_yyzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1683);

    auto tr_z_yyzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1684);

    auto tr_z_yyzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1685);

    auto tr_z_yyzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1686);

    auto tr_z_yyzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1687);

    auto tr_z_yyzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1688);

    auto tr_z_yyzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1689);

    auto tr_z_yyzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1690);

    auto tr_z_yyzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1691);

    auto tr_z_yyzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1692);

    auto tr_z_yyzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1693);

    auto tr_z_yyzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1694);

    auto tr_z_yyzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1695);

    auto tr_z_yyzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1696);

    auto tr_z_yyzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1697);

    auto tr_z_yyzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1698);

    auto tr_z_yyzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1699);

    auto tr_z_yyzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1700);

    auto tr_z_yyzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1701);

    auto tr_z_yyzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1702);

    auto tr_z_yyzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1703);

    auto tr_z_yyzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1704);

    auto tr_z_yyzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1705);

    auto tr_z_yyzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1706);

    auto tr_z_yyzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1707);

#pragma omp simd aligned(pa_y,                  \
                             tr_z_yyzzz_xxxxxx, \
                             tr_z_yyzzz_xxxxxy, \
                             tr_z_yyzzz_xxxxxz, \
                             tr_z_yyzzz_xxxxyy, \
                             tr_z_yyzzz_xxxxyz, \
                             tr_z_yyzzz_xxxxzz, \
                             tr_z_yyzzz_xxxyyy, \
                             tr_z_yyzzz_xxxyyz, \
                             tr_z_yyzzz_xxxyzz, \
                             tr_z_yyzzz_xxxzzz, \
                             tr_z_yyzzz_xxyyyy, \
                             tr_z_yyzzz_xxyyyz, \
                             tr_z_yyzzz_xxyyzz, \
                             tr_z_yyzzz_xxyzzz, \
                             tr_z_yyzzz_xxzzzz, \
                             tr_z_yyzzz_xyyyyy, \
                             tr_z_yyzzz_xyyyyz, \
                             tr_z_yyzzz_xyyyzz, \
                             tr_z_yyzzz_xyyzzz, \
                             tr_z_yyzzz_xyzzzz, \
                             tr_z_yyzzz_xzzzzz, \
                             tr_z_yyzzz_yyyyyy, \
                             tr_z_yyzzz_yyyyyz, \
                             tr_z_yyzzz_yyyyzz, \
                             tr_z_yyzzz_yyyzzz, \
                             tr_z_yyzzz_yyzzzz, \
                             tr_z_yyzzz_yzzzzz, \
                             tr_z_yyzzz_zzzzzz, \
                             tr_z_yzzz_xxxxx,   \
                             tr_z_yzzz_xxxxxx,  \
                             tr_z_yzzz_xxxxxy,  \
                             tr_z_yzzz_xxxxxz,  \
                             tr_z_yzzz_xxxxy,   \
                             tr_z_yzzz_xxxxyy,  \
                             tr_z_yzzz_xxxxyz,  \
                             tr_z_yzzz_xxxxz,   \
                             tr_z_yzzz_xxxxzz,  \
                             tr_z_yzzz_xxxyy,   \
                             tr_z_yzzz_xxxyyy,  \
                             tr_z_yzzz_xxxyyz,  \
                             tr_z_yzzz_xxxyz,   \
                             tr_z_yzzz_xxxyzz,  \
                             tr_z_yzzz_xxxzz,   \
                             tr_z_yzzz_xxxzzz,  \
                             tr_z_yzzz_xxyyy,   \
                             tr_z_yzzz_xxyyyy,  \
                             tr_z_yzzz_xxyyyz,  \
                             tr_z_yzzz_xxyyz,   \
                             tr_z_yzzz_xxyyzz,  \
                             tr_z_yzzz_xxyzz,   \
                             tr_z_yzzz_xxyzzz,  \
                             tr_z_yzzz_xxzzz,   \
                             tr_z_yzzz_xxzzzz,  \
                             tr_z_yzzz_xyyyy,   \
                             tr_z_yzzz_xyyyyy,  \
                             tr_z_yzzz_xyyyyz,  \
                             tr_z_yzzz_xyyyz,   \
                             tr_z_yzzz_xyyyzz,  \
                             tr_z_yzzz_xyyzz,   \
                             tr_z_yzzz_xyyzzz,  \
                             tr_z_yzzz_xyzzz,   \
                             tr_z_yzzz_xyzzzz,  \
                             tr_z_yzzz_xzzzz,   \
                             tr_z_yzzz_xzzzzz,  \
                             tr_z_yzzz_yyyyy,   \
                             tr_z_yzzz_yyyyyy,  \
                             tr_z_yzzz_yyyyyz,  \
                             tr_z_yzzz_yyyyz,   \
                             tr_z_yzzz_yyyyzz,  \
                             tr_z_yzzz_yyyzz,   \
                             tr_z_yzzz_yyyzzz,  \
                             tr_z_yzzz_yyzzz,   \
                             tr_z_yzzz_yyzzzz,  \
                             tr_z_yzzz_yzzzz,   \
                             tr_z_yzzz_yzzzzz,  \
                             tr_z_yzzz_zzzzz,   \
                             tr_z_yzzz_zzzzzz,  \
                             tr_z_zzz_xxxxxx,   \
                             tr_z_zzz_xxxxxy,   \
                             tr_z_zzz_xxxxxz,   \
                             tr_z_zzz_xxxxyy,   \
                             tr_z_zzz_xxxxyz,   \
                             tr_z_zzz_xxxxzz,   \
                             tr_z_zzz_xxxyyy,   \
                             tr_z_zzz_xxxyyz,   \
                             tr_z_zzz_xxxyzz,   \
                             tr_z_zzz_xxxzzz,   \
                             tr_z_zzz_xxyyyy,   \
                             tr_z_zzz_xxyyyz,   \
                             tr_z_zzz_xxyyzz,   \
                             tr_z_zzz_xxyzzz,   \
                             tr_z_zzz_xxzzzz,   \
                             tr_z_zzz_xyyyyy,   \
                             tr_z_zzz_xyyyyz,   \
                             tr_z_zzz_xyyyzz,   \
                             tr_z_zzz_xyyzzz,   \
                             tr_z_zzz_xyzzzz,   \
                             tr_z_zzz_xzzzzz,   \
                             tr_z_zzz_yyyyyy,   \
                             tr_z_zzz_yyyyyz,   \
                             tr_z_zzz_yyyyzz,   \
                             tr_z_zzz_yyyzzz,   \
                             tr_z_zzz_yyzzzz,   \
                             tr_z_zzz_yzzzzz,   \
                             tr_z_zzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzz_xxxxxx[i] = tr_z_zzz_xxxxxx[i] * fe_0 + tr_z_yzzz_xxxxxx[i] * pa_y[i];

        tr_z_yyzzz_xxxxxy[i] = tr_z_zzz_xxxxxy[i] * fe_0 + tr_z_yzzz_xxxxx[i] * fe_0 + tr_z_yzzz_xxxxxy[i] * pa_y[i];

        tr_z_yyzzz_xxxxxz[i] = tr_z_zzz_xxxxxz[i] * fe_0 + tr_z_yzzz_xxxxxz[i] * pa_y[i];

        tr_z_yyzzz_xxxxyy[i] = tr_z_zzz_xxxxyy[i] * fe_0 + 2.0 * tr_z_yzzz_xxxxy[i] * fe_0 + tr_z_yzzz_xxxxyy[i] * pa_y[i];

        tr_z_yyzzz_xxxxyz[i] = tr_z_zzz_xxxxyz[i] * fe_0 + tr_z_yzzz_xxxxz[i] * fe_0 + tr_z_yzzz_xxxxyz[i] * pa_y[i];

        tr_z_yyzzz_xxxxzz[i] = tr_z_zzz_xxxxzz[i] * fe_0 + tr_z_yzzz_xxxxzz[i] * pa_y[i];

        tr_z_yyzzz_xxxyyy[i] = tr_z_zzz_xxxyyy[i] * fe_0 + 3.0 * tr_z_yzzz_xxxyy[i] * fe_0 + tr_z_yzzz_xxxyyy[i] * pa_y[i];

        tr_z_yyzzz_xxxyyz[i] = tr_z_zzz_xxxyyz[i] * fe_0 + 2.0 * tr_z_yzzz_xxxyz[i] * fe_0 + tr_z_yzzz_xxxyyz[i] * pa_y[i];

        tr_z_yyzzz_xxxyzz[i] = tr_z_zzz_xxxyzz[i] * fe_0 + tr_z_yzzz_xxxzz[i] * fe_0 + tr_z_yzzz_xxxyzz[i] * pa_y[i];

        tr_z_yyzzz_xxxzzz[i] = tr_z_zzz_xxxzzz[i] * fe_0 + tr_z_yzzz_xxxzzz[i] * pa_y[i];

        tr_z_yyzzz_xxyyyy[i] = tr_z_zzz_xxyyyy[i] * fe_0 + 4.0 * tr_z_yzzz_xxyyy[i] * fe_0 + tr_z_yzzz_xxyyyy[i] * pa_y[i];

        tr_z_yyzzz_xxyyyz[i] = tr_z_zzz_xxyyyz[i] * fe_0 + 3.0 * tr_z_yzzz_xxyyz[i] * fe_0 + tr_z_yzzz_xxyyyz[i] * pa_y[i];

        tr_z_yyzzz_xxyyzz[i] = tr_z_zzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_yzzz_xxyzz[i] * fe_0 + tr_z_yzzz_xxyyzz[i] * pa_y[i];

        tr_z_yyzzz_xxyzzz[i] = tr_z_zzz_xxyzzz[i] * fe_0 + tr_z_yzzz_xxzzz[i] * fe_0 + tr_z_yzzz_xxyzzz[i] * pa_y[i];

        tr_z_yyzzz_xxzzzz[i] = tr_z_zzz_xxzzzz[i] * fe_0 + tr_z_yzzz_xxzzzz[i] * pa_y[i];

        tr_z_yyzzz_xyyyyy[i] = tr_z_zzz_xyyyyy[i] * fe_0 + 5.0 * tr_z_yzzz_xyyyy[i] * fe_0 + tr_z_yzzz_xyyyyy[i] * pa_y[i];

        tr_z_yyzzz_xyyyyz[i] = tr_z_zzz_xyyyyz[i] * fe_0 + 4.0 * tr_z_yzzz_xyyyz[i] * fe_0 + tr_z_yzzz_xyyyyz[i] * pa_y[i];

        tr_z_yyzzz_xyyyzz[i] = tr_z_zzz_xyyyzz[i] * fe_0 + 3.0 * tr_z_yzzz_xyyzz[i] * fe_0 + tr_z_yzzz_xyyyzz[i] * pa_y[i];

        tr_z_yyzzz_xyyzzz[i] = tr_z_zzz_xyyzzz[i] * fe_0 + 2.0 * tr_z_yzzz_xyzzz[i] * fe_0 + tr_z_yzzz_xyyzzz[i] * pa_y[i];

        tr_z_yyzzz_xyzzzz[i] = tr_z_zzz_xyzzzz[i] * fe_0 + tr_z_yzzz_xzzzz[i] * fe_0 + tr_z_yzzz_xyzzzz[i] * pa_y[i];

        tr_z_yyzzz_xzzzzz[i] = tr_z_zzz_xzzzzz[i] * fe_0 + tr_z_yzzz_xzzzzz[i] * pa_y[i];

        tr_z_yyzzz_yyyyyy[i] = tr_z_zzz_yyyyyy[i] * fe_0 + 6.0 * tr_z_yzzz_yyyyy[i] * fe_0 + tr_z_yzzz_yyyyyy[i] * pa_y[i];

        tr_z_yyzzz_yyyyyz[i] = tr_z_zzz_yyyyyz[i] * fe_0 + 5.0 * tr_z_yzzz_yyyyz[i] * fe_0 + tr_z_yzzz_yyyyyz[i] * pa_y[i];

        tr_z_yyzzz_yyyyzz[i] = tr_z_zzz_yyyyzz[i] * fe_0 + 4.0 * tr_z_yzzz_yyyzz[i] * fe_0 + tr_z_yzzz_yyyyzz[i] * pa_y[i];

        tr_z_yyzzz_yyyzzz[i] = tr_z_zzz_yyyzzz[i] * fe_0 + 3.0 * tr_z_yzzz_yyzzz[i] * fe_0 + tr_z_yzzz_yyyzzz[i] * pa_y[i];

        tr_z_yyzzz_yyzzzz[i] = tr_z_zzz_yyzzzz[i] * fe_0 + 2.0 * tr_z_yzzz_yzzzz[i] * fe_0 + tr_z_yzzz_yyzzzz[i] * pa_y[i];

        tr_z_yyzzz_yzzzzz[i] = tr_z_zzz_yzzzzz[i] * fe_0 + tr_z_yzzz_zzzzz[i] * fe_0 + tr_z_yzzz_yzzzzz[i] * pa_y[i];

        tr_z_yyzzz_zzzzzz[i] = tr_z_zzz_zzzzzz[i] * fe_0 + tr_z_yzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1708-1736 components of targeted buffer : HI

    auto tr_z_yzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1708);

    auto tr_z_yzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1709);

    auto tr_z_yzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1710);

    auto tr_z_yzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1711);

    auto tr_z_yzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1712);

    auto tr_z_yzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1713);

    auto tr_z_yzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1714);

    auto tr_z_yzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1715);

    auto tr_z_yzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1716);

    auto tr_z_yzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1717);

    auto tr_z_yzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1718);

    auto tr_z_yzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1719);

    auto tr_z_yzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1720);

    auto tr_z_yzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1721);

    auto tr_z_yzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1722);

    auto tr_z_yzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1723);

    auto tr_z_yzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1724);

    auto tr_z_yzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1725);

    auto tr_z_yzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1726);

    auto tr_z_yzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1727);

    auto tr_z_yzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1728);

    auto tr_z_yzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1729);

    auto tr_z_yzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1730);

    auto tr_z_yzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1731);

    auto tr_z_yzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1732);

    auto tr_z_yzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1733);

    auto tr_z_yzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1734);

    auto tr_z_yzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1735);

#pragma omp simd aligned(pa_y,                  \
                             tr_z_yzzzz_xxxxxx, \
                             tr_z_yzzzz_xxxxxy, \
                             tr_z_yzzzz_xxxxxz, \
                             tr_z_yzzzz_xxxxyy, \
                             tr_z_yzzzz_xxxxyz, \
                             tr_z_yzzzz_xxxxzz, \
                             tr_z_yzzzz_xxxyyy, \
                             tr_z_yzzzz_xxxyyz, \
                             tr_z_yzzzz_xxxyzz, \
                             tr_z_yzzzz_xxxzzz, \
                             tr_z_yzzzz_xxyyyy, \
                             tr_z_yzzzz_xxyyyz, \
                             tr_z_yzzzz_xxyyzz, \
                             tr_z_yzzzz_xxyzzz, \
                             tr_z_yzzzz_xxzzzz, \
                             tr_z_yzzzz_xyyyyy, \
                             tr_z_yzzzz_xyyyyz, \
                             tr_z_yzzzz_xyyyzz, \
                             tr_z_yzzzz_xyyzzz, \
                             tr_z_yzzzz_xyzzzz, \
                             tr_z_yzzzz_xzzzzz, \
                             tr_z_yzzzz_yyyyyy, \
                             tr_z_yzzzz_yyyyyz, \
                             tr_z_yzzzz_yyyyzz, \
                             tr_z_yzzzz_yyyzzz, \
                             tr_z_yzzzz_yyzzzz, \
                             tr_z_yzzzz_yzzzzz, \
                             tr_z_yzzzz_zzzzzz, \
                             tr_z_zzzz_xxxxx,   \
                             tr_z_zzzz_xxxxxx,  \
                             tr_z_zzzz_xxxxxy,  \
                             tr_z_zzzz_xxxxxz,  \
                             tr_z_zzzz_xxxxy,   \
                             tr_z_zzzz_xxxxyy,  \
                             tr_z_zzzz_xxxxyz,  \
                             tr_z_zzzz_xxxxz,   \
                             tr_z_zzzz_xxxxzz,  \
                             tr_z_zzzz_xxxyy,   \
                             tr_z_zzzz_xxxyyy,  \
                             tr_z_zzzz_xxxyyz,  \
                             tr_z_zzzz_xxxyz,   \
                             tr_z_zzzz_xxxyzz,  \
                             tr_z_zzzz_xxxzz,   \
                             tr_z_zzzz_xxxzzz,  \
                             tr_z_zzzz_xxyyy,   \
                             tr_z_zzzz_xxyyyy,  \
                             tr_z_zzzz_xxyyyz,  \
                             tr_z_zzzz_xxyyz,   \
                             tr_z_zzzz_xxyyzz,  \
                             tr_z_zzzz_xxyzz,   \
                             tr_z_zzzz_xxyzzz,  \
                             tr_z_zzzz_xxzzz,   \
                             tr_z_zzzz_xxzzzz,  \
                             tr_z_zzzz_xyyyy,   \
                             tr_z_zzzz_xyyyyy,  \
                             tr_z_zzzz_xyyyyz,  \
                             tr_z_zzzz_xyyyz,   \
                             tr_z_zzzz_xyyyzz,  \
                             tr_z_zzzz_xyyzz,   \
                             tr_z_zzzz_xyyzzz,  \
                             tr_z_zzzz_xyzzz,   \
                             tr_z_zzzz_xyzzzz,  \
                             tr_z_zzzz_xzzzz,   \
                             tr_z_zzzz_xzzzzz,  \
                             tr_z_zzzz_yyyyy,   \
                             tr_z_zzzz_yyyyyy,  \
                             tr_z_zzzz_yyyyyz,  \
                             tr_z_zzzz_yyyyz,   \
                             tr_z_zzzz_yyyyzz,  \
                             tr_z_zzzz_yyyzz,   \
                             tr_z_zzzz_yyyzzz,  \
                             tr_z_zzzz_yyzzz,   \
                             tr_z_zzzz_yyzzzz,  \
                             tr_z_zzzz_yzzzz,   \
                             tr_z_zzzz_yzzzzz,  \
                             tr_z_zzzz_zzzzz,   \
                             tr_z_zzzz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzz_xxxxxx[i] = tr_z_zzzz_xxxxxx[i] * pa_y[i];

        tr_z_yzzzz_xxxxxy[i] = tr_z_zzzz_xxxxx[i] * fe_0 + tr_z_zzzz_xxxxxy[i] * pa_y[i];

        tr_z_yzzzz_xxxxxz[i] = tr_z_zzzz_xxxxxz[i] * pa_y[i];

        tr_z_yzzzz_xxxxyy[i] = 2.0 * tr_z_zzzz_xxxxy[i] * fe_0 + tr_z_zzzz_xxxxyy[i] * pa_y[i];

        tr_z_yzzzz_xxxxyz[i] = tr_z_zzzz_xxxxz[i] * fe_0 + tr_z_zzzz_xxxxyz[i] * pa_y[i];

        tr_z_yzzzz_xxxxzz[i] = tr_z_zzzz_xxxxzz[i] * pa_y[i];

        tr_z_yzzzz_xxxyyy[i] = 3.0 * tr_z_zzzz_xxxyy[i] * fe_0 + tr_z_zzzz_xxxyyy[i] * pa_y[i];

        tr_z_yzzzz_xxxyyz[i] = 2.0 * tr_z_zzzz_xxxyz[i] * fe_0 + tr_z_zzzz_xxxyyz[i] * pa_y[i];

        tr_z_yzzzz_xxxyzz[i] = tr_z_zzzz_xxxzz[i] * fe_0 + tr_z_zzzz_xxxyzz[i] * pa_y[i];

        tr_z_yzzzz_xxxzzz[i] = tr_z_zzzz_xxxzzz[i] * pa_y[i];

        tr_z_yzzzz_xxyyyy[i] = 4.0 * tr_z_zzzz_xxyyy[i] * fe_0 + tr_z_zzzz_xxyyyy[i] * pa_y[i];

        tr_z_yzzzz_xxyyyz[i] = 3.0 * tr_z_zzzz_xxyyz[i] * fe_0 + tr_z_zzzz_xxyyyz[i] * pa_y[i];

        tr_z_yzzzz_xxyyzz[i] = 2.0 * tr_z_zzzz_xxyzz[i] * fe_0 + tr_z_zzzz_xxyyzz[i] * pa_y[i];

        tr_z_yzzzz_xxyzzz[i] = tr_z_zzzz_xxzzz[i] * fe_0 + tr_z_zzzz_xxyzzz[i] * pa_y[i];

        tr_z_yzzzz_xxzzzz[i] = tr_z_zzzz_xxzzzz[i] * pa_y[i];

        tr_z_yzzzz_xyyyyy[i] = 5.0 * tr_z_zzzz_xyyyy[i] * fe_0 + tr_z_zzzz_xyyyyy[i] * pa_y[i];

        tr_z_yzzzz_xyyyyz[i] = 4.0 * tr_z_zzzz_xyyyz[i] * fe_0 + tr_z_zzzz_xyyyyz[i] * pa_y[i];

        tr_z_yzzzz_xyyyzz[i] = 3.0 * tr_z_zzzz_xyyzz[i] * fe_0 + tr_z_zzzz_xyyyzz[i] * pa_y[i];

        tr_z_yzzzz_xyyzzz[i] = 2.0 * tr_z_zzzz_xyzzz[i] * fe_0 + tr_z_zzzz_xyyzzz[i] * pa_y[i];

        tr_z_yzzzz_xyzzzz[i] = tr_z_zzzz_xzzzz[i] * fe_0 + tr_z_zzzz_xyzzzz[i] * pa_y[i];

        tr_z_yzzzz_xzzzzz[i] = tr_z_zzzz_xzzzzz[i] * pa_y[i];

        tr_z_yzzzz_yyyyyy[i] = 6.0 * tr_z_zzzz_yyyyy[i] * fe_0 + tr_z_zzzz_yyyyyy[i] * pa_y[i];

        tr_z_yzzzz_yyyyyz[i] = 5.0 * tr_z_zzzz_yyyyz[i] * fe_0 + tr_z_zzzz_yyyyyz[i] * pa_y[i];

        tr_z_yzzzz_yyyyzz[i] = 4.0 * tr_z_zzzz_yyyzz[i] * fe_0 + tr_z_zzzz_yyyyzz[i] * pa_y[i];

        tr_z_yzzzz_yyyzzz[i] = 3.0 * tr_z_zzzz_yyzzz[i] * fe_0 + tr_z_zzzz_yyyzzz[i] * pa_y[i];

        tr_z_yzzzz_yyzzzz[i] = 2.0 * tr_z_zzzz_yzzzz[i] * fe_0 + tr_z_zzzz_yyzzzz[i] * pa_y[i];

        tr_z_yzzzz_yzzzzz[i] = tr_z_zzzz_zzzzz[i] * fe_0 + tr_z_zzzz_yzzzzz[i] * pa_y[i];

        tr_z_yzzzz_zzzzzz[i] = tr_z_zzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 1736-1764 components of targeted buffer : HI

    auto tr_z_zzzzz_xxxxxx = pbuffer.data(idx_dip_hi + 1736);

    auto tr_z_zzzzz_xxxxxy = pbuffer.data(idx_dip_hi + 1737);

    auto tr_z_zzzzz_xxxxxz = pbuffer.data(idx_dip_hi + 1738);

    auto tr_z_zzzzz_xxxxyy = pbuffer.data(idx_dip_hi + 1739);

    auto tr_z_zzzzz_xxxxyz = pbuffer.data(idx_dip_hi + 1740);

    auto tr_z_zzzzz_xxxxzz = pbuffer.data(idx_dip_hi + 1741);

    auto tr_z_zzzzz_xxxyyy = pbuffer.data(idx_dip_hi + 1742);

    auto tr_z_zzzzz_xxxyyz = pbuffer.data(idx_dip_hi + 1743);

    auto tr_z_zzzzz_xxxyzz = pbuffer.data(idx_dip_hi + 1744);

    auto tr_z_zzzzz_xxxzzz = pbuffer.data(idx_dip_hi + 1745);

    auto tr_z_zzzzz_xxyyyy = pbuffer.data(idx_dip_hi + 1746);

    auto tr_z_zzzzz_xxyyyz = pbuffer.data(idx_dip_hi + 1747);

    auto tr_z_zzzzz_xxyyzz = pbuffer.data(idx_dip_hi + 1748);

    auto tr_z_zzzzz_xxyzzz = pbuffer.data(idx_dip_hi + 1749);

    auto tr_z_zzzzz_xxzzzz = pbuffer.data(idx_dip_hi + 1750);

    auto tr_z_zzzzz_xyyyyy = pbuffer.data(idx_dip_hi + 1751);

    auto tr_z_zzzzz_xyyyyz = pbuffer.data(idx_dip_hi + 1752);

    auto tr_z_zzzzz_xyyyzz = pbuffer.data(idx_dip_hi + 1753);

    auto tr_z_zzzzz_xyyzzz = pbuffer.data(idx_dip_hi + 1754);

    auto tr_z_zzzzz_xyzzzz = pbuffer.data(idx_dip_hi + 1755);

    auto tr_z_zzzzz_xzzzzz = pbuffer.data(idx_dip_hi + 1756);

    auto tr_z_zzzzz_yyyyyy = pbuffer.data(idx_dip_hi + 1757);

    auto tr_z_zzzzz_yyyyyz = pbuffer.data(idx_dip_hi + 1758);

    auto tr_z_zzzzz_yyyyzz = pbuffer.data(idx_dip_hi + 1759);

    auto tr_z_zzzzz_yyyzzz = pbuffer.data(idx_dip_hi + 1760);

    auto tr_z_zzzzz_yyzzzz = pbuffer.data(idx_dip_hi + 1761);

    auto tr_z_zzzzz_yzzzzz = pbuffer.data(idx_dip_hi + 1762);

    auto tr_z_zzzzz_zzzzzz = pbuffer.data(idx_dip_hi + 1763);

#pragma omp simd aligned(pa_z,                  \
                             tr_z_zzz_xxxxxx,   \
                             tr_z_zzz_xxxxxy,   \
                             tr_z_zzz_xxxxxz,   \
                             tr_z_zzz_xxxxyy,   \
                             tr_z_zzz_xxxxyz,   \
                             tr_z_zzz_xxxxzz,   \
                             tr_z_zzz_xxxyyy,   \
                             tr_z_zzz_xxxyyz,   \
                             tr_z_zzz_xxxyzz,   \
                             tr_z_zzz_xxxzzz,   \
                             tr_z_zzz_xxyyyy,   \
                             tr_z_zzz_xxyyyz,   \
                             tr_z_zzz_xxyyzz,   \
                             tr_z_zzz_xxyzzz,   \
                             tr_z_zzz_xxzzzz,   \
                             tr_z_zzz_xyyyyy,   \
                             tr_z_zzz_xyyyyz,   \
                             tr_z_zzz_xyyyzz,   \
                             tr_z_zzz_xyyzzz,   \
                             tr_z_zzz_xyzzzz,   \
                             tr_z_zzz_xzzzzz,   \
                             tr_z_zzz_yyyyyy,   \
                             tr_z_zzz_yyyyyz,   \
                             tr_z_zzz_yyyyzz,   \
                             tr_z_zzz_yyyzzz,   \
                             tr_z_zzz_yyzzzz,   \
                             tr_z_zzz_yzzzzz,   \
                             tr_z_zzz_zzzzzz,   \
                             tr_z_zzzz_xxxxx,   \
                             tr_z_zzzz_xxxxxx,  \
                             tr_z_zzzz_xxxxxy,  \
                             tr_z_zzzz_xxxxxz,  \
                             tr_z_zzzz_xxxxy,   \
                             tr_z_zzzz_xxxxyy,  \
                             tr_z_zzzz_xxxxyz,  \
                             tr_z_zzzz_xxxxz,   \
                             tr_z_zzzz_xxxxzz,  \
                             tr_z_zzzz_xxxyy,   \
                             tr_z_zzzz_xxxyyy,  \
                             tr_z_zzzz_xxxyyz,  \
                             tr_z_zzzz_xxxyz,   \
                             tr_z_zzzz_xxxyzz,  \
                             tr_z_zzzz_xxxzz,   \
                             tr_z_zzzz_xxxzzz,  \
                             tr_z_zzzz_xxyyy,   \
                             tr_z_zzzz_xxyyyy,  \
                             tr_z_zzzz_xxyyyz,  \
                             tr_z_zzzz_xxyyz,   \
                             tr_z_zzzz_xxyyzz,  \
                             tr_z_zzzz_xxyzz,   \
                             tr_z_zzzz_xxyzzz,  \
                             tr_z_zzzz_xxzzz,   \
                             tr_z_zzzz_xxzzzz,  \
                             tr_z_zzzz_xyyyy,   \
                             tr_z_zzzz_xyyyyy,  \
                             tr_z_zzzz_xyyyyz,  \
                             tr_z_zzzz_xyyyz,   \
                             tr_z_zzzz_xyyyzz,  \
                             tr_z_zzzz_xyyzz,   \
                             tr_z_zzzz_xyyzzz,  \
                             tr_z_zzzz_xyzzz,   \
                             tr_z_zzzz_xyzzzz,  \
                             tr_z_zzzz_xzzzz,   \
                             tr_z_zzzz_xzzzzz,  \
                             tr_z_zzzz_yyyyy,   \
                             tr_z_zzzz_yyyyyy,  \
                             tr_z_zzzz_yyyyyz,  \
                             tr_z_zzzz_yyyyz,   \
                             tr_z_zzzz_yyyyzz,  \
                             tr_z_zzzz_yyyzz,   \
                             tr_z_zzzz_yyyzzz,  \
                             tr_z_zzzz_yyzzz,   \
                             tr_z_zzzz_yyzzzz,  \
                             tr_z_zzzz_yzzzz,   \
                             tr_z_zzzz_yzzzzz,  \
                             tr_z_zzzz_zzzzz,   \
                             tr_z_zzzz_zzzzzz,  \
                             tr_z_zzzzz_xxxxxx, \
                             tr_z_zzzzz_xxxxxy, \
                             tr_z_zzzzz_xxxxxz, \
                             tr_z_zzzzz_xxxxyy, \
                             tr_z_zzzzz_xxxxyz, \
                             tr_z_zzzzz_xxxxzz, \
                             tr_z_zzzzz_xxxyyy, \
                             tr_z_zzzzz_xxxyyz, \
                             tr_z_zzzzz_xxxyzz, \
                             tr_z_zzzzz_xxxzzz, \
                             tr_z_zzzzz_xxyyyy, \
                             tr_z_zzzzz_xxyyyz, \
                             tr_z_zzzzz_xxyyzz, \
                             tr_z_zzzzz_xxyzzz, \
                             tr_z_zzzzz_xxzzzz, \
                             tr_z_zzzzz_xyyyyy, \
                             tr_z_zzzzz_xyyyyz, \
                             tr_z_zzzzz_xyyyzz, \
                             tr_z_zzzzz_xyyzzz, \
                             tr_z_zzzzz_xyzzzz, \
                             tr_z_zzzzz_xzzzzz, \
                             tr_z_zzzzz_yyyyyy, \
                             tr_z_zzzzz_yyyyyz, \
                             tr_z_zzzzz_yyyyzz, \
                             tr_z_zzzzz_yyyzzz, \
                             tr_z_zzzzz_yyzzzz, \
                             tr_z_zzzzz_yzzzzz, \
                             tr_z_zzzzz_zzzzzz, \
                             ts_zzzz_xxxxxx,    \
                             ts_zzzz_xxxxxy,    \
                             ts_zzzz_xxxxxz,    \
                             ts_zzzz_xxxxyy,    \
                             ts_zzzz_xxxxyz,    \
                             ts_zzzz_xxxxzz,    \
                             ts_zzzz_xxxyyy,    \
                             ts_zzzz_xxxyyz,    \
                             ts_zzzz_xxxyzz,    \
                             ts_zzzz_xxxzzz,    \
                             ts_zzzz_xxyyyy,    \
                             ts_zzzz_xxyyyz,    \
                             ts_zzzz_xxyyzz,    \
                             ts_zzzz_xxyzzz,    \
                             ts_zzzz_xxzzzz,    \
                             ts_zzzz_xyyyyy,    \
                             ts_zzzz_xyyyyz,    \
                             ts_zzzz_xyyyzz,    \
                             ts_zzzz_xyyzzz,    \
                             ts_zzzz_xyzzzz,    \
                             ts_zzzz_xzzzzz,    \
                             ts_zzzz_yyyyyy,    \
                             ts_zzzz_yyyyyz,    \
                             ts_zzzz_yyyyzz,    \
                             ts_zzzz_yyyzzz,    \
                             ts_zzzz_yyzzzz,    \
                             ts_zzzz_yzzzzz,    \
                             ts_zzzz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzz_xxxxxx[i] = 4.0 * tr_z_zzz_xxxxxx[i] * fe_0 + ts_zzzz_xxxxxx[i] * fe_0 + tr_z_zzzz_xxxxxx[i] * pa_z[i];

        tr_z_zzzzz_xxxxxy[i] = 4.0 * tr_z_zzz_xxxxxy[i] * fe_0 + ts_zzzz_xxxxxy[i] * fe_0 + tr_z_zzzz_xxxxxy[i] * pa_z[i];

        tr_z_zzzzz_xxxxxz[i] = 4.0 * tr_z_zzz_xxxxxz[i] * fe_0 + tr_z_zzzz_xxxxx[i] * fe_0 + ts_zzzz_xxxxxz[i] * fe_0 + tr_z_zzzz_xxxxxz[i] * pa_z[i];

        tr_z_zzzzz_xxxxyy[i] = 4.0 * tr_z_zzz_xxxxyy[i] * fe_0 + ts_zzzz_xxxxyy[i] * fe_0 + tr_z_zzzz_xxxxyy[i] * pa_z[i];

        tr_z_zzzzz_xxxxyz[i] = 4.0 * tr_z_zzz_xxxxyz[i] * fe_0 + tr_z_zzzz_xxxxy[i] * fe_0 + ts_zzzz_xxxxyz[i] * fe_0 + tr_z_zzzz_xxxxyz[i] * pa_z[i];

        tr_z_zzzzz_xxxxzz[i] =
            4.0 * tr_z_zzz_xxxxzz[i] * fe_0 + 2.0 * tr_z_zzzz_xxxxz[i] * fe_0 + ts_zzzz_xxxxzz[i] * fe_0 + tr_z_zzzz_xxxxzz[i] * pa_z[i];

        tr_z_zzzzz_xxxyyy[i] = 4.0 * tr_z_zzz_xxxyyy[i] * fe_0 + ts_zzzz_xxxyyy[i] * fe_0 + tr_z_zzzz_xxxyyy[i] * pa_z[i];

        tr_z_zzzzz_xxxyyz[i] = 4.0 * tr_z_zzz_xxxyyz[i] * fe_0 + tr_z_zzzz_xxxyy[i] * fe_0 + ts_zzzz_xxxyyz[i] * fe_0 + tr_z_zzzz_xxxyyz[i] * pa_z[i];

        tr_z_zzzzz_xxxyzz[i] =
            4.0 * tr_z_zzz_xxxyzz[i] * fe_0 + 2.0 * tr_z_zzzz_xxxyz[i] * fe_0 + ts_zzzz_xxxyzz[i] * fe_0 + tr_z_zzzz_xxxyzz[i] * pa_z[i];

        tr_z_zzzzz_xxxzzz[i] =
            4.0 * tr_z_zzz_xxxzzz[i] * fe_0 + 3.0 * tr_z_zzzz_xxxzz[i] * fe_0 + ts_zzzz_xxxzzz[i] * fe_0 + tr_z_zzzz_xxxzzz[i] * pa_z[i];

        tr_z_zzzzz_xxyyyy[i] = 4.0 * tr_z_zzz_xxyyyy[i] * fe_0 + ts_zzzz_xxyyyy[i] * fe_0 + tr_z_zzzz_xxyyyy[i] * pa_z[i];

        tr_z_zzzzz_xxyyyz[i] = 4.0 * tr_z_zzz_xxyyyz[i] * fe_0 + tr_z_zzzz_xxyyy[i] * fe_0 + ts_zzzz_xxyyyz[i] * fe_0 + tr_z_zzzz_xxyyyz[i] * pa_z[i];

        tr_z_zzzzz_xxyyzz[i] =
            4.0 * tr_z_zzz_xxyyzz[i] * fe_0 + 2.0 * tr_z_zzzz_xxyyz[i] * fe_0 + ts_zzzz_xxyyzz[i] * fe_0 + tr_z_zzzz_xxyyzz[i] * pa_z[i];

        tr_z_zzzzz_xxyzzz[i] =
            4.0 * tr_z_zzz_xxyzzz[i] * fe_0 + 3.0 * tr_z_zzzz_xxyzz[i] * fe_0 + ts_zzzz_xxyzzz[i] * fe_0 + tr_z_zzzz_xxyzzz[i] * pa_z[i];

        tr_z_zzzzz_xxzzzz[i] =
            4.0 * tr_z_zzz_xxzzzz[i] * fe_0 + 4.0 * tr_z_zzzz_xxzzz[i] * fe_0 + ts_zzzz_xxzzzz[i] * fe_0 + tr_z_zzzz_xxzzzz[i] * pa_z[i];

        tr_z_zzzzz_xyyyyy[i] = 4.0 * tr_z_zzz_xyyyyy[i] * fe_0 + ts_zzzz_xyyyyy[i] * fe_0 + tr_z_zzzz_xyyyyy[i] * pa_z[i];

        tr_z_zzzzz_xyyyyz[i] = 4.0 * tr_z_zzz_xyyyyz[i] * fe_0 + tr_z_zzzz_xyyyy[i] * fe_0 + ts_zzzz_xyyyyz[i] * fe_0 + tr_z_zzzz_xyyyyz[i] * pa_z[i];

        tr_z_zzzzz_xyyyzz[i] =
            4.0 * tr_z_zzz_xyyyzz[i] * fe_0 + 2.0 * tr_z_zzzz_xyyyz[i] * fe_0 + ts_zzzz_xyyyzz[i] * fe_0 + tr_z_zzzz_xyyyzz[i] * pa_z[i];

        tr_z_zzzzz_xyyzzz[i] =
            4.0 * tr_z_zzz_xyyzzz[i] * fe_0 + 3.0 * tr_z_zzzz_xyyzz[i] * fe_0 + ts_zzzz_xyyzzz[i] * fe_0 + tr_z_zzzz_xyyzzz[i] * pa_z[i];

        tr_z_zzzzz_xyzzzz[i] =
            4.0 * tr_z_zzz_xyzzzz[i] * fe_0 + 4.0 * tr_z_zzzz_xyzzz[i] * fe_0 + ts_zzzz_xyzzzz[i] * fe_0 + tr_z_zzzz_xyzzzz[i] * pa_z[i];

        tr_z_zzzzz_xzzzzz[i] =
            4.0 * tr_z_zzz_xzzzzz[i] * fe_0 + 5.0 * tr_z_zzzz_xzzzz[i] * fe_0 + ts_zzzz_xzzzzz[i] * fe_0 + tr_z_zzzz_xzzzzz[i] * pa_z[i];

        tr_z_zzzzz_yyyyyy[i] = 4.0 * tr_z_zzz_yyyyyy[i] * fe_0 + ts_zzzz_yyyyyy[i] * fe_0 + tr_z_zzzz_yyyyyy[i] * pa_z[i];

        tr_z_zzzzz_yyyyyz[i] = 4.0 * tr_z_zzz_yyyyyz[i] * fe_0 + tr_z_zzzz_yyyyy[i] * fe_0 + ts_zzzz_yyyyyz[i] * fe_0 + tr_z_zzzz_yyyyyz[i] * pa_z[i];

        tr_z_zzzzz_yyyyzz[i] =
            4.0 * tr_z_zzz_yyyyzz[i] * fe_0 + 2.0 * tr_z_zzzz_yyyyz[i] * fe_0 + ts_zzzz_yyyyzz[i] * fe_0 + tr_z_zzzz_yyyyzz[i] * pa_z[i];

        tr_z_zzzzz_yyyzzz[i] =
            4.0 * tr_z_zzz_yyyzzz[i] * fe_0 + 3.0 * tr_z_zzzz_yyyzz[i] * fe_0 + ts_zzzz_yyyzzz[i] * fe_0 + tr_z_zzzz_yyyzzz[i] * pa_z[i];

        tr_z_zzzzz_yyzzzz[i] =
            4.0 * tr_z_zzz_yyzzzz[i] * fe_0 + 4.0 * tr_z_zzzz_yyzzz[i] * fe_0 + ts_zzzz_yyzzzz[i] * fe_0 + tr_z_zzzz_yyzzzz[i] * pa_z[i];

        tr_z_zzzzz_yzzzzz[i] =
            4.0 * tr_z_zzz_yzzzzz[i] * fe_0 + 5.0 * tr_z_zzzz_yzzzz[i] * fe_0 + ts_zzzz_yzzzzz[i] * fe_0 + tr_z_zzzz_yzzzzz[i] * pa_z[i];

        tr_z_zzzzz_zzzzzz[i] =
            4.0 * tr_z_zzz_zzzzzz[i] * fe_0 + 6.0 * tr_z_zzzz_zzzzz[i] * fe_0 + ts_zzzz_zzzzzz[i] * fe_0 + tr_z_zzzz_zzzzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
