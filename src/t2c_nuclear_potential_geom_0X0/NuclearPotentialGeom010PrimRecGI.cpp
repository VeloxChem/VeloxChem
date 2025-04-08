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

#include "NuclearPotentialGeom010PrimRecGI.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_gi(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_gi,
                                        const size_t              idx_npot_geom_010_0_di,
                                        const size_t              idx_npot_geom_010_1_di,
                                        const size_t              idx_npot_geom_010_0_fh,
                                        const size_t              idx_npot_geom_010_1_fh,
                                        const size_t              idx_npot_1_fi,
                                        const size_t              idx_npot_geom_010_0_fi,
                                        const size_t              idx_npot_geom_010_1_fi,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : DI

    auto ta1_x_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di);

    auto ta1_x_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 1);

    auto ta1_x_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 2);

    auto ta1_x_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 3);

    auto ta1_x_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 4);

    auto ta1_x_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 5);

    auto ta1_x_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 6);

    auto ta1_x_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 7);

    auto ta1_x_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 8);

    auto ta1_x_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 9);

    auto ta1_x_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 10);

    auto ta1_x_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 11);

    auto ta1_x_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 12);

    auto ta1_x_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 13);

    auto ta1_x_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 14);

    auto ta1_x_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 15);

    auto ta1_x_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 16);

    auto ta1_x_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 17);

    auto ta1_x_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 18);

    auto ta1_x_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 19);

    auto ta1_x_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 20);

    auto ta1_x_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 21);

    auto ta1_x_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 22);

    auto ta1_x_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 23);

    auto ta1_x_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 24);

    auto ta1_x_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 25);

    auto ta1_x_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 26);

    auto ta1_x_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 27);

    auto ta1_x_xy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 28);

    auto ta1_x_xy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 30);

    auto ta1_x_xy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 33);

    auto ta1_x_xy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 37);

    auto ta1_x_xy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 42);

    auto ta1_x_xy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 48);

    auto ta1_x_xz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 56);

    auto ta1_x_xz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 57);

    auto ta1_x_xz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 58);

    auto ta1_x_xz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 59);

    auto ta1_x_xz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 61);

    auto ta1_x_xz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 62);

    auto ta1_x_xz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 65);

    auto ta1_x_xz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 66);

    auto ta1_x_xz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 70);

    auto ta1_x_xz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 71);

    auto ta1_x_xz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 76);

    auto ta1_x_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 84);

    auto ta1_x_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 85);

    auto ta1_x_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 86);

    auto ta1_x_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 87);

    auto ta1_x_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 88);

    auto ta1_x_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 89);

    auto ta1_x_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 90);

    auto ta1_x_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 91);

    auto ta1_x_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 92);

    auto ta1_x_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 93);

    auto ta1_x_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 94);

    auto ta1_x_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 95);

    auto ta1_x_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 96);

    auto ta1_x_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 97);

    auto ta1_x_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 98);

    auto ta1_x_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 99);

    auto ta1_x_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 100);

    auto ta1_x_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 101);

    auto ta1_x_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 102);

    auto ta1_x_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 103);

    auto ta1_x_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 104);

    auto ta1_x_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 105);

    auto ta1_x_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 106);

    auto ta1_x_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 107);

    auto ta1_x_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 108);

    auto ta1_x_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 109);

    auto ta1_x_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 110);

    auto ta1_x_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 111);

    auto ta1_x_yz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 114);

    auto ta1_x_yz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 117);

    auto ta1_x_yz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 121);

    auto ta1_x_yz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 126);

    auto ta1_x_yz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 132);

    auto ta1_x_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 139);

    auto ta1_x_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 140);

    auto ta1_x_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 141);

    auto ta1_x_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 142);

    auto ta1_x_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 143);

    auto ta1_x_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 144);

    auto ta1_x_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 145);

    auto ta1_x_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 146);

    auto ta1_x_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 147);

    auto ta1_x_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 148);

    auto ta1_x_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 149);

    auto ta1_x_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 150);

    auto ta1_x_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 151);

    auto ta1_x_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 152);

    auto ta1_x_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 153);

    auto ta1_x_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 154);

    auto ta1_x_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 155);

    auto ta1_x_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 156);

    auto ta1_x_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 157);

    auto ta1_x_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 158);

    auto ta1_x_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 159);

    auto ta1_x_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 160);

    auto ta1_x_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 161);

    auto ta1_x_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 162);

    auto ta1_x_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 163);

    auto ta1_x_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 164);

    auto ta1_x_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 165);

    auto ta1_x_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 166);

    auto ta1_x_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 167);

    auto ta1_y_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 168);

    auto ta1_y_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 169);

    auto ta1_y_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 170);

    auto ta1_y_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 171);

    auto ta1_y_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 172);

    auto ta1_y_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 173);

    auto ta1_y_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 174);

    auto ta1_y_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 175);

    auto ta1_y_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 176);

    auto ta1_y_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 177);

    auto ta1_y_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 178);

    auto ta1_y_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 179);

    auto ta1_y_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 180);

    auto ta1_y_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 181);

    auto ta1_y_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 182);

    auto ta1_y_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 183);

    auto ta1_y_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 184);

    auto ta1_y_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 185);

    auto ta1_y_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 186);

    auto ta1_y_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 187);

    auto ta1_y_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 188);

    auto ta1_y_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 189);

    auto ta1_y_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 190);

    auto ta1_y_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 191);

    auto ta1_y_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 192);

    auto ta1_y_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 193);

    auto ta1_y_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 194);

    auto ta1_y_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 195);

    auto ta1_y_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 217);

    auto ta1_y_xy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 218);

    auto ta1_y_xy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 219);

    auto ta1_y_xy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 220);

    auto ta1_y_xy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 221);

    auto ta1_y_xy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 222);

    auto ta1_y_xz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 246);

    auto ta1_y_xz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 247);

    auto ta1_y_xz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 248);

    auto ta1_y_xz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 249);

    auto ta1_y_xz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 250);

    auto ta1_y_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 251);

    auto ta1_y_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 252);

    auto ta1_y_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 253);

    auto ta1_y_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 254);

    auto ta1_y_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 255);

    auto ta1_y_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 256);

    auto ta1_y_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 257);

    auto ta1_y_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 258);

    auto ta1_y_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 259);

    auto ta1_y_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 260);

    auto ta1_y_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 261);

    auto ta1_y_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 262);

    auto ta1_y_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 263);

    auto ta1_y_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 264);

    auto ta1_y_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 265);

    auto ta1_y_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 266);

    auto ta1_y_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 267);

    auto ta1_y_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 268);

    auto ta1_y_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 269);

    auto ta1_y_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 270);

    auto ta1_y_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 271);

    auto ta1_y_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 272);

    auto ta1_y_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 273);

    auto ta1_y_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 274);

    auto ta1_y_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 275);

    auto ta1_y_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 276);

    auto ta1_y_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 277);

    auto ta1_y_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 278);

    auto ta1_y_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 279);

    auto ta1_y_yz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 281);

    auto ta1_y_yz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 283);

    auto ta1_y_yz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 286);

    auto ta1_y_yz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 290);

    auto ta1_y_yz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 295);

    auto ta1_y_yz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 301);

    auto ta1_y_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 302);

    auto ta1_y_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 303);

    auto ta1_y_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 304);

    auto ta1_y_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 305);

    auto ta1_y_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 306);

    auto ta1_y_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 308);

    auto ta1_y_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 309);

    auto ta1_y_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 310);

    auto ta1_y_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 311);

    auto ta1_y_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 312);

    auto ta1_y_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 313);

    auto ta1_y_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 314);

    auto ta1_y_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 315);

    auto ta1_y_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 316);

    auto ta1_y_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 317);

    auto ta1_y_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 318);

    auto ta1_y_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 319);

    auto ta1_y_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 320);

    auto ta1_y_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 321);

    auto ta1_y_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 322);

    auto ta1_y_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 323);

    auto ta1_y_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 324);

    auto ta1_y_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 325);

    auto ta1_y_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 326);

    auto ta1_y_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 327);

    auto ta1_y_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 328);

    auto ta1_y_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 329);

    auto ta1_y_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 330);

    auto ta1_y_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 331);

    auto ta1_y_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 332);

    auto ta1_y_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 333);

    auto ta1_y_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 334);

    auto ta1_y_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 335);

    auto ta1_z_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 336);

    auto ta1_z_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 337);

    auto ta1_z_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 338);

    auto ta1_z_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 339);

    auto ta1_z_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 340);

    auto ta1_z_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 341);

    auto ta1_z_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 342);

    auto ta1_z_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 343);

    auto ta1_z_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 344);

    auto ta1_z_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 345);

    auto ta1_z_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 346);

    auto ta1_z_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 347);

    auto ta1_z_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 348);

    auto ta1_z_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 349);

    auto ta1_z_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 350);

    auto ta1_z_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 351);

    auto ta1_z_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 352);

    auto ta1_z_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 353);

    auto ta1_z_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 354);

    auto ta1_z_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 355);

    auto ta1_z_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 356);

    auto ta1_z_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 357);

    auto ta1_z_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 358);

    auto ta1_z_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 359);

    auto ta1_z_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 360);

    auto ta1_z_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 361);

    auto ta1_z_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 362);

    auto ta1_z_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 363);

    auto ta1_z_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 385);

    auto ta1_z_xy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 386);

    auto ta1_z_xy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 387);

    auto ta1_z_xy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 388);

    auto ta1_z_xy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 389);

    auto ta1_z_xy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 390);

    auto ta1_z_xz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 414);

    auto ta1_z_xz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 415);

    auto ta1_z_xz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 416);

    auto ta1_z_xz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 417);

    auto ta1_z_xz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 418);

    auto ta1_z_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 419);

    auto ta1_z_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 420);

    auto ta1_z_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 421);

    auto ta1_z_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 422);

    auto ta1_z_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 423);

    auto ta1_z_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 424);

    auto ta1_z_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 425);

    auto ta1_z_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 426);

    auto ta1_z_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 427);

    auto ta1_z_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 428);

    auto ta1_z_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 429);

    auto ta1_z_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 430);

    auto ta1_z_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 431);

    auto ta1_z_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 432);

    auto ta1_z_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 433);

    auto ta1_z_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 434);

    auto ta1_z_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 435);

    auto ta1_z_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 436);

    auto ta1_z_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 437);

    auto ta1_z_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 438);

    auto ta1_z_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 439);

    auto ta1_z_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 440);

    auto ta1_z_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 441);

    auto ta1_z_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 442);

    auto ta1_z_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 443);

    auto ta1_z_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 444);

    auto ta1_z_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 445);

    auto ta1_z_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 446);

    auto ta1_z_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 447);

    auto ta1_z_yz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 450);

    auto ta1_z_yz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 453);

    auto ta1_z_yz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 457);

    auto ta1_z_yz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 462);

    auto ta1_z_yz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 468);

    auto ta1_z_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 470);

    auto ta1_z_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 471);

    auto ta1_z_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 472);

    auto ta1_z_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 473);

    auto ta1_z_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 474);

    auto ta1_z_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 475);

    auto ta1_z_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 476);

    auto ta1_z_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 477);

    auto ta1_z_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 478);

    auto ta1_z_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 479);

    auto ta1_z_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 480);

    auto ta1_z_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 481);

    auto ta1_z_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 482);

    auto ta1_z_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 483);

    auto ta1_z_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 484);

    auto ta1_z_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 485);

    auto ta1_z_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 486);

    auto ta1_z_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 487);

    auto ta1_z_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 488);

    auto ta1_z_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 489);

    auto ta1_z_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 490);

    auto ta1_z_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 491);

    auto ta1_z_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 492);

    auto ta1_z_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 493);

    auto ta1_z_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 494);

    auto ta1_z_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 495);

    auto ta1_z_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 496);

    auto ta1_z_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 497);

    auto ta1_z_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 498);

    auto ta1_z_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 499);

    auto ta1_z_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 500);

    auto ta1_z_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 501);

    auto ta1_z_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 502);

    auto ta1_z_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 503);

    // Set up components of auxiliary buffer : DI

    auto ta1_x_xx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di);

    auto ta1_x_xx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 1);

    auto ta1_x_xx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 2);

    auto ta1_x_xx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 3);

    auto ta1_x_xx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 4);

    auto ta1_x_xx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 5);

    auto ta1_x_xx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 6);

    auto ta1_x_xx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 7);

    auto ta1_x_xx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 8);

    auto ta1_x_xx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 9);

    auto ta1_x_xx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 10);

    auto ta1_x_xx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 11);

    auto ta1_x_xx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 12);

    auto ta1_x_xx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 13);

    auto ta1_x_xx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 14);

    auto ta1_x_xx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 15);

    auto ta1_x_xx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 16);

    auto ta1_x_xx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 17);

    auto ta1_x_xx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 18);

    auto ta1_x_xx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 19);

    auto ta1_x_xx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 20);

    auto ta1_x_xx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 21);

    auto ta1_x_xx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 22);

    auto ta1_x_xx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 23);

    auto ta1_x_xx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 24);

    auto ta1_x_xx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 25);

    auto ta1_x_xx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 26);

    auto ta1_x_xx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 27);

    auto ta1_x_xy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 28);

    auto ta1_x_xy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 30);

    auto ta1_x_xy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 33);

    auto ta1_x_xy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 37);

    auto ta1_x_xy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 42);

    auto ta1_x_xy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 48);

    auto ta1_x_xz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 56);

    auto ta1_x_xz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 57);

    auto ta1_x_xz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 58);

    auto ta1_x_xz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 59);

    auto ta1_x_xz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 61);

    auto ta1_x_xz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 62);

    auto ta1_x_xz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 65);

    auto ta1_x_xz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 66);

    auto ta1_x_xz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 70);

    auto ta1_x_xz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 71);

    auto ta1_x_xz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 76);

    auto ta1_x_yy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 84);

    auto ta1_x_yy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 85);

    auto ta1_x_yy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 86);

    auto ta1_x_yy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 87);

    auto ta1_x_yy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 88);

    auto ta1_x_yy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 89);

    auto ta1_x_yy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 90);

    auto ta1_x_yy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 91);

    auto ta1_x_yy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 92);

    auto ta1_x_yy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 93);

    auto ta1_x_yy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 94);

    auto ta1_x_yy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 95);

    auto ta1_x_yy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 96);

    auto ta1_x_yy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 97);

    auto ta1_x_yy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 98);

    auto ta1_x_yy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 99);

    auto ta1_x_yy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 100);

    auto ta1_x_yy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 101);

    auto ta1_x_yy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 102);

    auto ta1_x_yy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 103);

    auto ta1_x_yy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 104);

    auto ta1_x_yy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 105);

    auto ta1_x_yy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 106);

    auto ta1_x_yy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 107);

    auto ta1_x_yy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 108);

    auto ta1_x_yy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 109);

    auto ta1_x_yy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 110);

    auto ta1_x_yy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 111);

    auto ta1_x_yz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 114);

    auto ta1_x_yz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 117);

    auto ta1_x_yz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 121);

    auto ta1_x_yz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 126);

    auto ta1_x_yz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 132);

    auto ta1_x_yz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 139);

    auto ta1_x_zz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 140);

    auto ta1_x_zz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 141);

    auto ta1_x_zz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 142);

    auto ta1_x_zz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 143);

    auto ta1_x_zz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 144);

    auto ta1_x_zz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 145);

    auto ta1_x_zz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 146);

    auto ta1_x_zz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 147);

    auto ta1_x_zz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 148);

    auto ta1_x_zz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 149);

    auto ta1_x_zz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 150);

    auto ta1_x_zz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 151);

    auto ta1_x_zz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 152);

    auto ta1_x_zz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 153);

    auto ta1_x_zz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 154);

    auto ta1_x_zz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 155);

    auto ta1_x_zz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 156);

    auto ta1_x_zz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 157);

    auto ta1_x_zz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 158);

    auto ta1_x_zz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 159);

    auto ta1_x_zz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 160);

    auto ta1_x_zz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 161);

    auto ta1_x_zz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 162);

    auto ta1_x_zz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 163);

    auto ta1_x_zz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 164);

    auto ta1_x_zz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 165);

    auto ta1_x_zz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 166);

    auto ta1_x_zz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 167);

    auto ta1_y_xx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 168);

    auto ta1_y_xx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 169);

    auto ta1_y_xx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 170);

    auto ta1_y_xx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 171);

    auto ta1_y_xx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 172);

    auto ta1_y_xx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 173);

    auto ta1_y_xx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 174);

    auto ta1_y_xx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 175);

    auto ta1_y_xx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 176);

    auto ta1_y_xx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 177);

    auto ta1_y_xx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 178);

    auto ta1_y_xx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 179);

    auto ta1_y_xx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 180);

    auto ta1_y_xx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 181);

    auto ta1_y_xx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 182);

    auto ta1_y_xx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 183);

    auto ta1_y_xx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 184);

    auto ta1_y_xx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 185);

    auto ta1_y_xx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 186);

    auto ta1_y_xx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 187);

    auto ta1_y_xx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 188);

    auto ta1_y_xx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 189);

    auto ta1_y_xx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 190);

    auto ta1_y_xx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 191);

    auto ta1_y_xx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 192);

    auto ta1_y_xx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 193);

    auto ta1_y_xx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 194);

    auto ta1_y_xx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 195);

    auto ta1_y_xy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 217);

    auto ta1_y_xy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 218);

    auto ta1_y_xy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 219);

    auto ta1_y_xy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 220);

    auto ta1_y_xy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 221);

    auto ta1_y_xy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 222);

    auto ta1_y_xz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 246);

    auto ta1_y_xz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 247);

    auto ta1_y_xz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 248);

    auto ta1_y_xz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 249);

    auto ta1_y_xz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 250);

    auto ta1_y_xz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 251);

    auto ta1_y_yy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 252);

    auto ta1_y_yy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 253);

    auto ta1_y_yy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 254);

    auto ta1_y_yy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 255);

    auto ta1_y_yy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 256);

    auto ta1_y_yy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 257);

    auto ta1_y_yy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 258);

    auto ta1_y_yy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 259);

    auto ta1_y_yy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 260);

    auto ta1_y_yy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 261);

    auto ta1_y_yy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 262);

    auto ta1_y_yy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 263);

    auto ta1_y_yy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 264);

    auto ta1_y_yy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 265);

    auto ta1_y_yy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 266);

    auto ta1_y_yy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 267);

    auto ta1_y_yy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 268);

    auto ta1_y_yy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 269);

    auto ta1_y_yy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 270);

    auto ta1_y_yy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 271);

    auto ta1_y_yy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 272);

    auto ta1_y_yy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 273);

    auto ta1_y_yy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 274);

    auto ta1_y_yy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 275);

    auto ta1_y_yy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 276);

    auto ta1_y_yy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 277);

    auto ta1_y_yy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 278);

    auto ta1_y_yy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 279);

    auto ta1_y_yz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 281);

    auto ta1_y_yz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 283);

    auto ta1_y_yz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 286);

    auto ta1_y_yz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 290);

    auto ta1_y_yz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 295);

    auto ta1_y_yz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 301);

    auto ta1_y_yz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 302);

    auto ta1_y_yz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 303);

    auto ta1_y_yz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 304);

    auto ta1_y_yz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 305);

    auto ta1_y_yz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 306);

    auto ta1_y_zz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 308);

    auto ta1_y_zz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 309);

    auto ta1_y_zz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 310);

    auto ta1_y_zz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 311);

    auto ta1_y_zz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 312);

    auto ta1_y_zz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 313);

    auto ta1_y_zz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 314);

    auto ta1_y_zz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 315);

    auto ta1_y_zz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 316);

    auto ta1_y_zz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 317);

    auto ta1_y_zz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 318);

    auto ta1_y_zz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 319);

    auto ta1_y_zz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 320);

    auto ta1_y_zz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 321);

    auto ta1_y_zz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 322);

    auto ta1_y_zz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 323);

    auto ta1_y_zz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 324);

    auto ta1_y_zz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 325);

    auto ta1_y_zz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 326);

    auto ta1_y_zz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 327);

    auto ta1_y_zz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 328);

    auto ta1_y_zz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 329);

    auto ta1_y_zz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 330);

    auto ta1_y_zz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 331);

    auto ta1_y_zz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 332);

    auto ta1_y_zz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 333);

    auto ta1_y_zz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 334);

    auto ta1_y_zz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 335);

    auto ta1_z_xx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 336);

    auto ta1_z_xx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 337);

    auto ta1_z_xx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 338);

    auto ta1_z_xx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 339);

    auto ta1_z_xx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 340);

    auto ta1_z_xx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 341);

    auto ta1_z_xx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 342);

    auto ta1_z_xx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 343);

    auto ta1_z_xx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 344);

    auto ta1_z_xx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 345);

    auto ta1_z_xx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 346);

    auto ta1_z_xx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 347);

    auto ta1_z_xx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 348);

    auto ta1_z_xx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 349);

    auto ta1_z_xx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 350);

    auto ta1_z_xx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 351);

    auto ta1_z_xx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 352);

    auto ta1_z_xx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 353);

    auto ta1_z_xx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 354);

    auto ta1_z_xx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 355);

    auto ta1_z_xx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 356);

    auto ta1_z_xx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 357);

    auto ta1_z_xx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 358);

    auto ta1_z_xx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 359);

    auto ta1_z_xx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 360);

    auto ta1_z_xx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 361);

    auto ta1_z_xx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 362);

    auto ta1_z_xx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 363);

    auto ta1_z_xy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 385);

    auto ta1_z_xy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 386);

    auto ta1_z_xy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 387);

    auto ta1_z_xy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 388);

    auto ta1_z_xy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 389);

    auto ta1_z_xy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 390);

    auto ta1_z_xz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 414);

    auto ta1_z_xz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 415);

    auto ta1_z_xz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 416);

    auto ta1_z_xz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 417);

    auto ta1_z_xz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 418);

    auto ta1_z_xz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 419);

    auto ta1_z_yy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 420);

    auto ta1_z_yy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 421);

    auto ta1_z_yy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 422);

    auto ta1_z_yy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 423);

    auto ta1_z_yy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 424);

    auto ta1_z_yy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 425);

    auto ta1_z_yy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 426);

    auto ta1_z_yy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 427);

    auto ta1_z_yy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 428);

    auto ta1_z_yy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 429);

    auto ta1_z_yy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 430);

    auto ta1_z_yy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 431);

    auto ta1_z_yy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 432);

    auto ta1_z_yy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 433);

    auto ta1_z_yy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 434);

    auto ta1_z_yy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 435);

    auto ta1_z_yy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 436);

    auto ta1_z_yy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 437);

    auto ta1_z_yy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 438);

    auto ta1_z_yy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 439);

    auto ta1_z_yy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 440);

    auto ta1_z_yy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 441);

    auto ta1_z_yy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 442);

    auto ta1_z_yy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 443);

    auto ta1_z_yy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 444);

    auto ta1_z_yy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 445);

    auto ta1_z_yy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 446);

    auto ta1_z_yy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 447);

    auto ta1_z_yz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 450);

    auto ta1_z_yz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 453);

    auto ta1_z_yz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 457);

    auto ta1_z_yz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 462);

    auto ta1_z_yz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 468);

    auto ta1_z_yz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 470);

    auto ta1_z_yz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 471);

    auto ta1_z_yz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 472);

    auto ta1_z_yz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 473);

    auto ta1_z_yz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 474);

    auto ta1_z_yz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 475);

    auto ta1_z_zz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 476);

    auto ta1_z_zz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 477);

    auto ta1_z_zz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 478);

    auto ta1_z_zz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 479);

    auto ta1_z_zz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 480);

    auto ta1_z_zz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 481);

    auto ta1_z_zz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 482);

    auto ta1_z_zz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 483);

    auto ta1_z_zz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 484);

    auto ta1_z_zz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 485);

    auto ta1_z_zz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 486);

    auto ta1_z_zz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 487);

    auto ta1_z_zz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 488);

    auto ta1_z_zz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 489);

    auto ta1_z_zz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 490);

    auto ta1_z_zz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 491);

    auto ta1_z_zz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 492);

    auto ta1_z_zz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 493);

    auto ta1_z_zz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 494);

    auto ta1_z_zz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 495);

    auto ta1_z_zz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 496);

    auto ta1_z_zz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 497);

    auto ta1_z_zz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 498);

    auto ta1_z_zz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 499);

    auto ta1_z_zz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 500);

    auto ta1_z_zz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 501);

    auto ta1_z_zz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 502);

    auto ta1_z_zz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 503);

    // Set up components of auxiliary buffer : FH

    auto ta1_x_xxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh);

    auto ta1_x_xxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 1);

    auto ta1_x_xxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 2);

    auto ta1_x_xxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 3);

    auto ta1_x_xxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 4);

    auto ta1_x_xxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 5);

    auto ta1_x_xxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 6);

    auto ta1_x_xxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 7);

    auto ta1_x_xxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 8);

    auto ta1_x_xxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 9);

    auto ta1_x_xxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 10);

    auto ta1_x_xxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 11);

    auto ta1_x_xxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 12);

    auto ta1_x_xxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 13);

    auto ta1_x_xxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 14);

    auto ta1_x_xxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 15);

    auto ta1_x_xxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 16);

    auto ta1_x_xxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 17);

    auto ta1_x_xxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 18);

    auto ta1_x_xxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 19);

    auto ta1_x_xxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 20);

    auto ta1_x_xxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 21);

    auto ta1_x_xxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 22);

    auto ta1_x_xxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 23);

    auto ta1_x_xxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 24);

    auto ta1_x_xxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 25);

    auto ta1_x_xxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 26);

    auto ta1_x_xxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 27);

    auto ta1_x_xxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 28);

    auto ta1_x_xxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 29);

    auto ta1_x_xxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 30);

    auto ta1_x_xxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 31);

    auto ta1_x_xxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 32);

    auto ta1_x_xxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 33);

    auto ta1_x_xxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 34);

    auto ta1_x_xxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 35);

    auto ta1_x_xxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 42);

    auto ta1_x_xxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 43);

    auto ta1_x_xxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 44);

    auto ta1_x_xxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 45);

    auto ta1_x_xxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 46);

    auto ta1_x_xxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 47);

    auto ta1_x_xxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 48);

    auto ta1_x_xxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 49);

    auto ta1_x_xxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 50);

    auto ta1_x_xxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 51);

    auto ta1_x_xxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 52);

    auto ta1_x_xxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 53);

    auto ta1_x_xxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 54);

    auto ta1_x_xxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 55);

    auto ta1_x_xxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 56);

    auto ta1_x_xxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 58);

    auto ta1_x_xxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 59);

    auto ta1_x_xxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 60);

    auto ta1_x_xxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 61);

    auto ta1_x_xxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 62);

    auto ta1_x_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 64);

    auto ta1_x_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 66);

    auto ta1_x_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 67);

    auto ta1_x_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 69);

    auto ta1_x_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 70);

    auto ta1_x_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 71);

    auto ta1_x_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 73);

    auto ta1_x_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 74);

    auto ta1_x_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 75);

    auto ta1_x_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 76);

    auto ta1_x_xzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 105);

    auto ta1_x_xzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 106);

    auto ta1_x_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 107);

    auto ta1_x_xzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 108);

    auto ta1_x_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 109);

    auto ta1_x_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 110);

    auto ta1_x_xzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 111);

    auto ta1_x_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 112);

    auto ta1_x_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 113);

    auto ta1_x_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 114);

    auto ta1_x_xzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 115);

    auto ta1_x_xzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 116);

    auto ta1_x_xzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 117);

    auto ta1_x_xzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 118);

    auto ta1_x_xzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 119);

    auto ta1_x_yyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 126);

    auto ta1_x_yyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 127);

    auto ta1_x_yyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 128);

    auto ta1_x_yyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 129);

    auto ta1_x_yyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 130);

    auto ta1_x_yyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 131);

    auto ta1_x_yyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 132);

    auto ta1_x_yyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 133);

    auto ta1_x_yyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 134);

    auto ta1_x_yyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 135);

    auto ta1_x_yyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 136);

    auto ta1_x_yyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 137);

    auto ta1_x_yyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 138);

    auto ta1_x_yyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 139);

    auto ta1_x_yyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 140);

    auto ta1_x_yyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 141);

    auto ta1_x_yyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 142);

    auto ta1_x_yyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 143);

    auto ta1_x_yyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 144);

    auto ta1_x_yyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 145);

    auto ta1_x_yyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 146);

    auto ta1_x_yzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 170);

    auto ta1_x_yzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 172);

    auto ta1_x_yzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 173);

    auto ta1_x_yzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 175);

    auto ta1_x_yzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 176);

    auto ta1_x_yzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 177);

    auto ta1_x_yzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 179);

    auto ta1_x_yzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 180);

    auto ta1_x_yzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 181);

    auto ta1_x_yzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 182);

    auto ta1_x_yzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 184);

    auto ta1_x_yzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 185);

    auto ta1_x_yzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 186);

    auto ta1_x_yzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 187);

    auto ta1_x_yzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 188);

    auto ta1_x_zzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 189);

    auto ta1_x_zzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 190);

    auto ta1_x_zzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 191);

    auto ta1_x_zzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 192);

    auto ta1_x_zzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 193);

    auto ta1_x_zzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 194);

    auto ta1_x_zzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 195);

    auto ta1_x_zzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 196);

    auto ta1_x_zzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 197);

    auto ta1_x_zzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 198);

    auto ta1_x_zzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 199);

    auto ta1_x_zzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 200);

    auto ta1_x_zzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 201);

    auto ta1_x_zzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 202);

    auto ta1_x_zzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 203);

    auto ta1_x_zzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 204);

    auto ta1_x_zzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 205);

    auto ta1_x_zzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 206);

    auto ta1_x_zzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 207);

    auto ta1_x_zzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 208);

    auto ta1_x_zzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 209);

    auto ta1_y_xxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 210);

    auto ta1_y_xxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 211);

    auto ta1_y_xxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 212);

    auto ta1_y_xxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 213);

    auto ta1_y_xxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 214);

    auto ta1_y_xxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 215);

    auto ta1_y_xxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 216);

    auto ta1_y_xxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 217);

    auto ta1_y_xxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 218);

    auto ta1_y_xxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 219);

    auto ta1_y_xxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 220);

    auto ta1_y_xxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 221);

    auto ta1_y_xxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 222);

    auto ta1_y_xxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 223);

    auto ta1_y_xxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 224);

    auto ta1_y_xxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 225);

    auto ta1_y_xxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 226);

    auto ta1_y_xxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 227);

    auto ta1_y_xxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 228);

    auto ta1_y_xxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 229);

    auto ta1_y_xxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 230);

    auto ta1_y_xxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 232);

    auto ta1_y_xxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 234);

    auto ta1_y_xxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 235);

    auto ta1_y_xxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 237);

    auto ta1_y_xxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 238);

    auto ta1_y_xxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 239);

    auto ta1_y_xxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 241);

    auto ta1_y_xxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 242);

    auto ta1_y_xxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 243);

    auto ta1_y_xxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 244);

    auto ta1_y_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 274);

    auto ta1_y_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 276);

    auto ta1_y_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 277);

    auto ta1_y_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 279);

    auto ta1_y_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 280);

    auto ta1_y_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 281);

    auto ta1_y_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 283);

    auto ta1_y_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 284);

    auto ta1_y_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 285);

    auto ta1_y_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 286);

    auto ta1_y_xyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 288);

    auto ta1_y_xyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 289);

    auto ta1_y_xyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 290);

    auto ta1_y_xyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 291);

    auto ta1_y_xyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 292);

    auto ta1_y_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 317);

    auto ta1_y_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 319);

    auto ta1_y_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 320);

    auto ta1_y_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 322);

    auto ta1_y_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 323);

    auto ta1_y_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 324);

    auto ta1_y_xzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 326);

    auto ta1_y_xzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 327);

    auto ta1_y_xzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 328);

    auto ta1_y_xzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 329);

    auto ta1_y_xzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 331);

    auto ta1_y_xzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 332);

    auto ta1_y_xzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 333);

    auto ta1_y_xzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 334);

    auto ta1_y_xzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 335);

    auto ta1_y_yyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 336);

    auto ta1_y_yyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 337);

    auto ta1_y_yyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 338);

    auto ta1_y_yyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 339);

    auto ta1_y_yyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 340);

    auto ta1_y_yyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 341);

    auto ta1_y_yyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 342);

    auto ta1_y_yyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 343);

    auto ta1_y_yyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 344);

    auto ta1_y_yyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 345);

    auto ta1_y_yyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 346);

    auto ta1_y_yyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 347);

    auto ta1_y_yyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 348);

    auto ta1_y_yyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 349);

    auto ta1_y_yyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 350);

    auto ta1_y_yyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 351);

    auto ta1_y_yyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 352);

    auto ta1_y_yyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 353);

    auto ta1_y_yyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 354);

    auto ta1_y_yyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 355);

    auto ta1_y_yyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 356);

    auto ta1_y_yyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 358);

    auto ta1_y_yyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 359);

    auto ta1_y_yyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 360);

    auto ta1_y_yyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 361);

    auto ta1_y_yyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 362);

    auto ta1_y_yyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 363);

    auto ta1_y_yyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 364);

    auto ta1_y_yyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 365);

    auto ta1_y_yyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 366);

    auto ta1_y_yyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 367);

    auto ta1_y_yyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 368);

    auto ta1_y_yyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 369);

    auto ta1_y_yyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 370);

    auto ta1_y_yyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 371);

    auto ta1_y_yyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 372);

    auto ta1_y_yyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 373);

    auto ta1_y_yyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 374);

    auto ta1_y_yyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 375);

    auto ta1_y_yyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 376);

    auto ta1_y_yyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 377);

    auto ta1_y_yzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 379);

    auto ta1_y_yzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 381);

    auto ta1_y_yzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 382);

    auto ta1_y_yzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 384);

    auto ta1_y_yzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 385);

    auto ta1_y_yzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 386);

    auto ta1_y_yzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 388);

    auto ta1_y_yzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 389);

    auto ta1_y_yzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 390);

    auto ta1_y_yzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 391);

    auto ta1_y_yzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 393);

    auto ta1_y_yzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 394);

    auto ta1_y_yzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 395);

    auto ta1_y_yzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 396);

    auto ta1_y_yzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 397);

    auto ta1_y_zzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 399);

    auto ta1_y_zzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 400);

    auto ta1_y_zzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 401);

    auto ta1_y_zzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 402);

    auto ta1_y_zzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 403);

    auto ta1_y_zzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 404);

    auto ta1_y_zzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 405);

    auto ta1_y_zzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 406);

    auto ta1_y_zzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 407);

    auto ta1_y_zzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 408);

    auto ta1_y_zzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 409);

    auto ta1_y_zzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 410);

    auto ta1_y_zzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 411);

    auto ta1_y_zzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 412);

    auto ta1_y_zzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 413);

    auto ta1_y_zzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 414);

    auto ta1_y_zzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 415);

    auto ta1_y_zzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 416);

    auto ta1_y_zzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 417);

    auto ta1_y_zzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 418);

    auto ta1_y_zzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 419);

    auto ta1_z_xxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 420);

    auto ta1_z_xxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 421);

    auto ta1_z_xxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 422);

    auto ta1_z_xxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 423);

    auto ta1_z_xxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 424);

    auto ta1_z_xxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 425);

    auto ta1_z_xxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 426);

    auto ta1_z_xxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 427);

    auto ta1_z_xxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 428);

    auto ta1_z_xxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 429);

    auto ta1_z_xxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 430);

    auto ta1_z_xxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 431);

    auto ta1_z_xxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 432);

    auto ta1_z_xxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 433);

    auto ta1_z_xxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 434);

    auto ta1_z_xxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 435);

    auto ta1_z_xxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 436);

    auto ta1_z_xxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 437);

    auto ta1_z_xxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 438);

    auto ta1_z_xxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 439);

    auto ta1_z_xxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 440);

    auto ta1_z_xxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 464);

    auto ta1_z_xxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 466);

    auto ta1_z_xxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 467);

    auto ta1_z_xxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 469);

    auto ta1_z_xxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 470);

    auto ta1_z_xxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 471);

    auto ta1_z_xxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 473);

    auto ta1_z_xxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 474);

    auto ta1_z_xxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 475);

    auto ta1_z_xxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 476);

    auto ta1_z_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 484);

    auto ta1_z_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 486);

    auto ta1_z_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 487);

    auto ta1_z_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 489);

    auto ta1_z_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 490);

    auto ta1_z_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 491);

    auto ta1_z_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 493);

    auto ta1_z_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 494);

    auto ta1_z_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 495);

    auto ta1_z_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 496);

    auto ta1_z_xyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 498);

    auto ta1_z_xyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 499);

    auto ta1_z_xyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 500);

    auto ta1_z_xyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 501);

    auto ta1_z_xyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 502);

    auto ta1_z_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 527);

    auto ta1_z_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 529);

    auto ta1_z_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 530);

    auto ta1_z_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 532);

    auto ta1_z_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 533);

    auto ta1_z_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 534);

    auto ta1_z_xzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 536);

    auto ta1_z_xzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 537);

    auto ta1_z_xzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 538);

    auto ta1_z_xzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 539);

    auto ta1_z_xzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 541);

    auto ta1_z_xzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 542);

    auto ta1_z_xzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 543);

    auto ta1_z_xzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 544);

    auto ta1_z_xzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 545);

    auto ta1_z_yyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 546);

    auto ta1_z_yyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 547);

    auto ta1_z_yyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 548);

    auto ta1_z_yyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 549);

    auto ta1_z_yyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 550);

    auto ta1_z_yyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 551);

    auto ta1_z_yyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 552);

    auto ta1_z_yyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 553);

    auto ta1_z_yyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 554);

    auto ta1_z_yyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 555);

    auto ta1_z_yyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 556);

    auto ta1_z_yyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 557);

    auto ta1_z_yyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 558);

    auto ta1_z_yyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 559);

    auto ta1_z_yyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 560);

    auto ta1_z_yyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 561);

    auto ta1_z_yyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 562);

    auto ta1_z_yyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 563);

    auto ta1_z_yyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 564);

    auto ta1_z_yyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 565);

    auto ta1_z_yyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 566);

    auto ta1_z_yyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 569);

    auto ta1_z_yyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 571);

    auto ta1_z_yyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 572);

    auto ta1_z_yyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 574);

    auto ta1_z_yyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 575);

    auto ta1_z_yyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 576);

    auto ta1_z_yyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 578);

    auto ta1_z_yyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 579);

    auto ta1_z_yyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 580);

    auto ta1_z_yyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 581);

    auto ta1_z_yyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 583);

    auto ta1_z_yyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 584);

    auto ta1_z_yyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 585);

    auto ta1_z_yyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 586);

    auto ta1_z_yyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 587);

    auto ta1_z_yzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 589);

    auto ta1_z_yzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 590);

    auto ta1_z_yzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 591);

    auto ta1_z_yzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 592);

    auto ta1_z_yzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 593);

    auto ta1_z_yzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 594);

    auto ta1_z_yzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 595);

    auto ta1_z_yzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 596);

    auto ta1_z_yzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 597);

    auto ta1_z_yzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 598);

    auto ta1_z_yzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 599);

    auto ta1_z_yzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 600);

    auto ta1_z_yzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 601);

    auto ta1_z_yzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 602);

    auto ta1_z_yzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 603);

    auto ta1_z_yzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 604);

    auto ta1_z_yzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 605);

    auto ta1_z_yzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 606);

    auto ta1_z_yzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 607);

    auto ta1_z_yzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 608);

    auto ta1_z_zzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 609);

    auto ta1_z_zzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 610);

    auto ta1_z_zzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 611);

    auto ta1_z_zzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 612);

    auto ta1_z_zzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 613);

    auto ta1_z_zzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 614);

    auto ta1_z_zzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 615);

    auto ta1_z_zzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 616);

    auto ta1_z_zzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 617);

    auto ta1_z_zzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 618);

    auto ta1_z_zzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 619);

    auto ta1_z_zzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 620);

    auto ta1_z_zzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 621);

    auto ta1_z_zzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 622);

    auto ta1_z_zzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 623);

    auto ta1_z_zzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 624);

    auto ta1_z_zzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 625);

    auto ta1_z_zzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 626);

    auto ta1_z_zzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 627);

    auto ta1_z_zzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 628);

    auto ta1_z_zzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 629);

    // Set up components of auxiliary buffer : FH

    auto ta1_x_xxx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh);

    auto ta1_x_xxx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 1);

    auto ta1_x_xxx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 2);

    auto ta1_x_xxx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 3);

    auto ta1_x_xxx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 4);

    auto ta1_x_xxx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 5);

    auto ta1_x_xxx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 6);

    auto ta1_x_xxx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 7);

    auto ta1_x_xxx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 8);

    auto ta1_x_xxx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 9);

    auto ta1_x_xxx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 10);

    auto ta1_x_xxx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 11);

    auto ta1_x_xxx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 12);

    auto ta1_x_xxx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 13);

    auto ta1_x_xxx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 14);

    auto ta1_x_xxx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 15);

    auto ta1_x_xxx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 16);

    auto ta1_x_xxx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 17);

    auto ta1_x_xxx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 18);

    auto ta1_x_xxx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 19);

    auto ta1_x_xxx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 20);

    auto ta1_x_xxy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 21);

    auto ta1_x_xxy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 22);

    auto ta1_x_xxy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 23);

    auto ta1_x_xxy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 24);

    auto ta1_x_xxy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 25);

    auto ta1_x_xxy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 26);

    auto ta1_x_xxy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 27);

    auto ta1_x_xxy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 28);

    auto ta1_x_xxy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 29);

    auto ta1_x_xxy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 30);

    auto ta1_x_xxy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 31);

    auto ta1_x_xxy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 32);

    auto ta1_x_xxy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 33);

    auto ta1_x_xxy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 34);

    auto ta1_x_xxy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 35);

    auto ta1_x_xxz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 42);

    auto ta1_x_xxz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 43);

    auto ta1_x_xxz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 44);

    auto ta1_x_xxz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 45);

    auto ta1_x_xxz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 46);

    auto ta1_x_xxz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 47);

    auto ta1_x_xxz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 48);

    auto ta1_x_xxz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 49);

    auto ta1_x_xxz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 50);

    auto ta1_x_xxz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 51);

    auto ta1_x_xxz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 52);

    auto ta1_x_xxz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 53);

    auto ta1_x_xxz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 54);

    auto ta1_x_xxz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 55);

    auto ta1_x_xxz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 56);

    auto ta1_x_xxz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 58);

    auto ta1_x_xxz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 59);

    auto ta1_x_xxz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 60);

    auto ta1_x_xxz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 61);

    auto ta1_x_xxz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 62);

    auto ta1_x_xyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 64);

    auto ta1_x_xyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 66);

    auto ta1_x_xyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 67);

    auto ta1_x_xyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 69);

    auto ta1_x_xyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 70);

    auto ta1_x_xyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 71);

    auto ta1_x_xyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 73);

    auto ta1_x_xyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 74);

    auto ta1_x_xyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 75);

    auto ta1_x_xyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 76);

    auto ta1_x_xzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 105);

    auto ta1_x_xzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 106);

    auto ta1_x_xzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 107);

    auto ta1_x_xzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 108);

    auto ta1_x_xzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 109);

    auto ta1_x_xzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 110);

    auto ta1_x_xzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 111);

    auto ta1_x_xzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 112);

    auto ta1_x_xzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 113);

    auto ta1_x_xzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 114);

    auto ta1_x_xzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 115);

    auto ta1_x_xzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 116);

    auto ta1_x_xzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 117);

    auto ta1_x_xzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 118);

    auto ta1_x_xzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 119);

    auto ta1_x_yyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 126);

    auto ta1_x_yyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 127);

    auto ta1_x_yyy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 128);

    auto ta1_x_yyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 129);

    auto ta1_x_yyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 130);

    auto ta1_x_yyy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 131);

    auto ta1_x_yyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 132);

    auto ta1_x_yyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 133);

    auto ta1_x_yyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 134);

    auto ta1_x_yyy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 135);

    auto ta1_x_yyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 136);

    auto ta1_x_yyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 137);

    auto ta1_x_yyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 138);

    auto ta1_x_yyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 139);

    auto ta1_x_yyy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 140);

    auto ta1_x_yyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 141);

    auto ta1_x_yyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 142);

    auto ta1_x_yyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 143);

    auto ta1_x_yyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 144);

    auto ta1_x_yyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 145);

    auto ta1_x_yyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 146);

    auto ta1_x_yzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 170);

    auto ta1_x_yzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 172);

    auto ta1_x_yzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 173);

    auto ta1_x_yzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 175);

    auto ta1_x_yzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 176);

    auto ta1_x_yzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 177);

    auto ta1_x_yzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 179);

    auto ta1_x_yzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 180);

    auto ta1_x_yzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 181);

    auto ta1_x_yzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 182);

    auto ta1_x_yzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 184);

    auto ta1_x_yzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 185);

    auto ta1_x_yzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 186);

    auto ta1_x_yzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 187);

    auto ta1_x_yzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 188);

    auto ta1_x_zzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 189);

    auto ta1_x_zzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 190);

    auto ta1_x_zzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 191);

    auto ta1_x_zzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 192);

    auto ta1_x_zzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 193);

    auto ta1_x_zzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 194);

    auto ta1_x_zzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 195);

    auto ta1_x_zzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 196);

    auto ta1_x_zzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 197);

    auto ta1_x_zzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 198);

    auto ta1_x_zzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 199);

    auto ta1_x_zzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 200);

    auto ta1_x_zzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 201);

    auto ta1_x_zzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 202);

    auto ta1_x_zzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 203);

    auto ta1_x_zzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 204);

    auto ta1_x_zzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 205);

    auto ta1_x_zzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 206);

    auto ta1_x_zzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 207);

    auto ta1_x_zzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 208);

    auto ta1_x_zzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 209);

    auto ta1_y_xxx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 210);

    auto ta1_y_xxx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 211);

    auto ta1_y_xxx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 212);

    auto ta1_y_xxx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 213);

    auto ta1_y_xxx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 214);

    auto ta1_y_xxx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 215);

    auto ta1_y_xxx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 216);

    auto ta1_y_xxx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 217);

    auto ta1_y_xxx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 218);

    auto ta1_y_xxx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 219);

    auto ta1_y_xxx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 220);

    auto ta1_y_xxx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 221);

    auto ta1_y_xxx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 222);

    auto ta1_y_xxx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 223);

    auto ta1_y_xxx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 224);

    auto ta1_y_xxx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 225);

    auto ta1_y_xxx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 226);

    auto ta1_y_xxx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 227);

    auto ta1_y_xxx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 228);

    auto ta1_y_xxx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 229);

    auto ta1_y_xxx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 230);

    auto ta1_y_xxy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 232);

    auto ta1_y_xxy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 234);

    auto ta1_y_xxy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 235);

    auto ta1_y_xxy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 237);

    auto ta1_y_xxy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 238);

    auto ta1_y_xxy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 239);

    auto ta1_y_xxy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 241);

    auto ta1_y_xxy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 242);

    auto ta1_y_xxy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 243);

    auto ta1_y_xxy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 244);

    auto ta1_y_xyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 274);

    auto ta1_y_xyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 276);

    auto ta1_y_xyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 277);

    auto ta1_y_xyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 279);

    auto ta1_y_xyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 280);

    auto ta1_y_xyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 281);

    auto ta1_y_xyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 283);

    auto ta1_y_xyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 284);

    auto ta1_y_xyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 285);

    auto ta1_y_xyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 286);

    auto ta1_y_xyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 288);

    auto ta1_y_xyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 289);

    auto ta1_y_xyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 290);

    auto ta1_y_xyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 291);

    auto ta1_y_xyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 292);

    auto ta1_y_xzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 317);

    auto ta1_y_xzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 319);

    auto ta1_y_xzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 320);

    auto ta1_y_xzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 322);

    auto ta1_y_xzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 323);

    auto ta1_y_xzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 324);

    auto ta1_y_xzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 326);

    auto ta1_y_xzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 327);

    auto ta1_y_xzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 328);

    auto ta1_y_xzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 329);

    auto ta1_y_xzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 331);

    auto ta1_y_xzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 332);

    auto ta1_y_xzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 333);

    auto ta1_y_xzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 334);

    auto ta1_y_xzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 335);

    auto ta1_y_yyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 336);

    auto ta1_y_yyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 337);

    auto ta1_y_yyy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 338);

    auto ta1_y_yyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 339);

    auto ta1_y_yyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 340);

    auto ta1_y_yyy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 341);

    auto ta1_y_yyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 342);

    auto ta1_y_yyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 343);

    auto ta1_y_yyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 344);

    auto ta1_y_yyy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 345);

    auto ta1_y_yyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 346);

    auto ta1_y_yyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 347);

    auto ta1_y_yyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 348);

    auto ta1_y_yyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 349);

    auto ta1_y_yyy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 350);

    auto ta1_y_yyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 351);

    auto ta1_y_yyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 352);

    auto ta1_y_yyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 353);

    auto ta1_y_yyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 354);

    auto ta1_y_yyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 355);

    auto ta1_y_yyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 356);

    auto ta1_y_yyz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 358);

    auto ta1_y_yyz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 359);

    auto ta1_y_yyz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 360);

    auto ta1_y_yyz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 361);

    auto ta1_y_yyz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 362);

    auto ta1_y_yyz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 363);

    auto ta1_y_yyz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 364);

    auto ta1_y_yyz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 365);

    auto ta1_y_yyz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 366);

    auto ta1_y_yyz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 367);

    auto ta1_y_yyz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 368);

    auto ta1_y_yyz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 369);

    auto ta1_y_yyz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 370);

    auto ta1_y_yyz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 371);

    auto ta1_y_yyz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 372);

    auto ta1_y_yyz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 373);

    auto ta1_y_yyz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 374);

    auto ta1_y_yyz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 375);

    auto ta1_y_yyz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 376);

    auto ta1_y_yyz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 377);

    auto ta1_y_yzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 379);

    auto ta1_y_yzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 381);

    auto ta1_y_yzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 382);

    auto ta1_y_yzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 384);

    auto ta1_y_yzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 385);

    auto ta1_y_yzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 386);

    auto ta1_y_yzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 388);

    auto ta1_y_yzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 389);

    auto ta1_y_yzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 390);

    auto ta1_y_yzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 391);

    auto ta1_y_yzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 393);

    auto ta1_y_yzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 394);

    auto ta1_y_yzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 395);

    auto ta1_y_yzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 396);

    auto ta1_y_yzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 397);

    auto ta1_y_zzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 399);

    auto ta1_y_zzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 400);

    auto ta1_y_zzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 401);

    auto ta1_y_zzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 402);

    auto ta1_y_zzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 403);

    auto ta1_y_zzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 404);

    auto ta1_y_zzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 405);

    auto ta1_y_zzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 406);

    auto ta1_y_zzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 407);

    auto ta1_y_zzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 408);

    auto ta1_y_zzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 409);

    auto ta1_y_zzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 410);

    auto ta1_y_zzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 411);

    auto ta1_y_zzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 412);

    auto ta1_y_zzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 413);

    auto ta1_y_zzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 414);

    auto ta1_y_zzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 415);

    auto ta1_y_zzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 416);

    auto ta1_y_zzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 417);

    auto ta1_y_zzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 418);

    auto ta1_y_zzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 419);

    auto ta1_z_xxx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 420);

    auto ta1_z_xxx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 421);

    auto ta1_z_xxx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 422);

    auto ta1_z_xxx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 423);

    auto ta1_z_xxx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 424);

    auto ta1_z_xxx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 425);

    auto ta1_z_xxx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 426);

    auto ta1_z_xxx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 427);

    auto ta1_z_xxx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 428);

    auto ta1_z_xxx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 429);

    auto ta1_z_xxx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 430);

    auto ta1_z_xxx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 431);

    auto ta1_z_xxx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 432);

    auto ta1_z_xxx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 433);

    auto ta1_z_xxx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 434);

    auto ta1_z_xxx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 435);

    auto ta1_z_xxx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 436);

    auto ta1_z_xxx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 437);

    auto ta1_z_xxx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 438);

    auto ta1_z_xxx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 439);

    auto ta1_z_xxx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 440);

    auto ta1_z_xxz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 464);

    auto ta1_z_xxz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 466);

    auto ta1_z_xxz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 467);

    auto ta1_z_xxz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 469);

    auto ta1_z_xxz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 470);

    auto ta1_z_xxz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 471);

    auto ta1_z_xxz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 473);

    auto ta1_z_xxz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 474);

    auto ta1_z_xxz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 475);

    auto ta1_z_xxz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 476);

    auto ta1_z_xyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 484);

    auto ta1_z_xyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 486);

    auto ta1_z_xyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 487);

    auto ta1_z_xyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 489);

    auto ta1_z_xyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 490);

    auto ta1_z_xyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 491);

    auto ta1_z_xyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 493);

    auto ta1_z_xyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 494);

    auto ta1_z_xyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 495);

    auto ta1_z_xyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 496);

    auto ta1_z_xyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 498);

    auto ta1_z_xyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 499);

    auto ta1_z_xyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 500);

    auto ta1_z_xyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 501);

    auto ta1_z_xyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 502);

    auto ta1_z_xzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 527);

    auto ta1_z_xzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 529);

    auto ta1_z_xzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 530);

    auto ta1_z_xzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 532);

    auto ta1_z_xzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 533);

    auto ta1_z_xzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 534);

    auto ta1_z_xzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 536);

    auto ta1_z_xzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 537);

    auto ta1_z_xzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 538);

    auto ta1_z_xzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 539);

    auto ta1_z_xzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 541);

    auto ta1_z_xzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 542);

    auto ta1_z_xzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 543);

    auto ta1_z_xzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 544);

    auto ta1_z_xzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 545);

    auto ta1_z_yyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 546);

    auto ta1_z_yyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 547);

    auto ta1_z_yyy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 548);

    auto ta1_z_yyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 549);

    auto ta1_z_yyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 550);

    auto ta1_z_yyy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 551);

    auto ta1_z_yyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 552);

    auto ta1_z_yyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 553);

    auto ta1_z_yyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 554);

    auto ta1_z_yyy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 555);

    auto ta1_z_yyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 556);

    auto ta1_z_yyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 557);

    auto ta1_z_yyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 558);

    auto ta1_z_yyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 559);

    auto ta1_z_yyy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 560);

    auto ta1_z_yyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 561);

    auto ta1_z_yyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 562);

    auto ta1_z_yyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 563);

    auto ta1_z_yyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 564);

    auto ta1_z_yyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 565);

    auto ta1_z_yyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 566);

    auto ta1_z_yyz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 569);

    auto ta1_z_yyz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 571);

    auto ta1_z_yyz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 572);

    auto ta1_z_yyz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 574);

    auto ta1_z_yyz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 575);

    auto ta1_z_yyz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 576);

    auto ta1_z_yyz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 578);

    auto ta1_z_yyz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 579);

    auto ta1_z_yyz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 580);

    auto ta1_z_yyz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 581);

    auto ta1_z_yyz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 583);

    auto ta1_z_yyz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 584);

    auto ta1_z_yyz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 585);

    auto ta1_z_yyz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 586);

    auto ta1_z_yyz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 587);

    auto ta1_z_yzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 589);

    auto ta1_z_yzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 590);

    auto ta1_z_yzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 591);

    auto ta1_z_yzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 592);

    auto ta1_z_yzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 593);

    auto ta1_z_yzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 594);

    auto ta1_z_yzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 595);

    auto ta1_z_yzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 596);

    auto ta1_z_yzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 597);

    auto ta1_z_yzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 598);

    auto ta1_z_yzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 599);

    auto ta1_z_yzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 600);

    auto ta1_z_yzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 601);

    auto ta1_z_yzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 602);

    auto ta1_z_yzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 603);

    auto ta1_z_yzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 604);

    auto ta1_z_yzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 605);

    auto ta1_z_yzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 606);

    auto ta1_z_yzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 607);

    auto ta1_z_yzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 608);

    auto ta1_z_zzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 609);

    auto ta1_z_zzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 610);

    auto ta1_z_zzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 611);

    auto ta1_z_zzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 612);

    auto ta1_z_zzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 613);

    auto ta1_z_zzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 614);

    auto ta1_z_zzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 615);

    auto ta1_z_zzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 616);

    auto ta1_z_zzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 617);

    auto ta1_z_zzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 618);

    auto ta1_z_zzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 619);

    auto ta1_z_zzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 620);

    auto ta1_z_zzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 621);

    auto ta1_z_zzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 622);

    auto ta1_z_zzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 623);

    auto ta1_z_zzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 624);

    auto ta1_z_zzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 625);

    auto ta1_z_zzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 626);

    auto ta1_z_zzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 627);

    auto ta1_z_zzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 628);

    auto ta1_z_zzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 629);

    // Set up components of auxiliary buffer : FI

    auto ta_xxx_xxxxxx_1 = pbuffer.data(idx_npot_1_fi);

    auto ta_xxx_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 1);

    auto ta_xxx_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 2);

    auto ta_xxx_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 3);

    auto ta_xxx_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 4);

    auto ta_xxx_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 5);

    auto ta_xxx_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 6);

    auto ta_xxx_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 7);

    auto ta_xxx_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 8);

    auto ta_xxx_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 9);

    auto ta_xxx_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 10);

    auto ta_xxx_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 11);

    auto ta_xxx_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 12);

    auto ta_xxx_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 13);

    auto ta_xxx_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 14);

    auto ta_xxx_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 15);

    auto ta_xxx_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 16);

    auto ta_xxx_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 17);

    auto ta_xxx_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 18);

    auto ta_xxx_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 19);

    auto ta_xxx_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 20);

    auto ta_xxx_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 21);

    auto ta_xxx_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 22);

    auto ta_xxx_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 23);

    auto ta_xxx_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 24);

    auto ta_xxx_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 25);

    auto ta_xxx_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 26);

    auto ta_xxx_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 27);

    auto ta_xxy_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 28);

    auto ta_xxy_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 29);

    auto ta_xxy_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 30);

    auto ta_xxy_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 31);

    auto ta_xxy_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 33);

    auto ta_xxy_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 34);

    auto ta_xxy_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 37);

    auto ta_xxy_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 38);

    auto ta_xxy_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 42);

    auto ta_xxy_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 43);

    auto ta_xxy_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 48);

    auto ta_xxy_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 49);

    auto ta_xxz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 56);

    auto ta_xxz_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 57);

    auto ta_xxz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 58);

    auto ta_xxz_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 59);

    auto ta_xxz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 61);

    auto ta_xxz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 62);

    auto ta_xxz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 65);

    auto ta_xxz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 66);

    auto ta_xxz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 70);

    auto ta_xxz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 71);

    auto ta_xxz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 76);

    auto ta_xxz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 83);

    auto ta_xyy_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 84);

    auto ta_xyy_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 85);

    auto ta_xyy_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 87);

    auto ta_xyy_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 90);

    auto ta_xyy_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 94);

    auto ta_xyy_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 99);

    auto ta_xyy_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 105);

    auto ta_xyy_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 106);

    auto ta_xyy_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 107);

    auto ta_xyy_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 108);

    auto ta_xyy_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 109);

    auto ta_xyy_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 110);

    auto ta_xzz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 140);

    auto ta_xzz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 142);

    auto ta_xzz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 145);

    auto ta_xzz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 149);

    auto ta_xzz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 154);

    auto ta_xzz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 160);

    auto ta_xzz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 162);

    auto ta_xzz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 163);

    auto ta_xzz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 164);

    auto ta_xzz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 165);

    auto ta_xzz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 166);

    auto ta_xzz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 167);

    auto ta_yyy_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 168);

    auto ta_yyy_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 169);

    auto ta_yyy_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 170);

    auto ta_yyy_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 171);

    auto ta_yyy_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 172);

    auto ta_yyy_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 173);

    auto ta_yyy_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 174);

    auto ta_yyy_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 175);

    auto ta_yyy_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 176);

    auto ta_yyy_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 177);

    auto ta_yyy_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 178);

    auto ta_yyy_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 179);

    auto ta_yyy_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 180);

    auto ta_yyy_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 181);

    auto ta_yyy_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 182);

    auto ta_yyy_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 183);

    auto ta_yyy_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 184);

    auto ta_yyy_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 185);

    auto ta_yyy_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 186);

    auto ta_yyy_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 187);

    auto ta_yyy_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 188);

    auto ta_yyy_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 189);

    auto ta_yyy_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 190);

    auto ta_yyy_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 191);

    auto ta_yyy_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 192);

    auto ta_yyy_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 193);

    auto ta_yyy_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 194);

    auto ta_yyy_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 195);

    auto ta_yyz_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 197);

    auto ta_yyz_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 199);

    auto ta_yyz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 202);

    auto ta_yyz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 206);

    auto ta_yyz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 211);

    auto ta_yyz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 217);

    auto ta_yyz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 218);

    auto ta_yyz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 219);

    auto ta_yyz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 220);

    auto ta_yyz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 221);

    auto ta_yyz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 222);

    auto ta_yyz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 223);

    auto ta_yzz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 226);

    auto ta_yzz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 229);

    auto ta_yzz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 233);

    auto ta_yzz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 238);

    auto ta_yzz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 244);

    auto ta_yzz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 245);

    auto ta_yzz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 246);

    auto ta_yzz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 247);

    auto ta_yzz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 248);

    auto ta_yzz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 249);

    auto ta_yzz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 250);

    auto ta_yzz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 251);

    auto ta_zzz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 252);

    auto ta_zzz_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 253);

    auto ta_zzz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 254);

    auto ta_zzz_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 255);

    auto ta_zzz_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 256);

    auto ta_zzz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 257);

    auto ta_zzz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 258);

    auto ta_zzz_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 259);

    auto ta_zzz_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 260);

    auto ta_zzz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 261);

    auto ta_zzz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 262);

    auto ta_zzz_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 263);

    auto ta_zzz_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 264);

    auto ta_zzz_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 265);

    auto ta_zzz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 266);

    auto ta_zzz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 267);

    auto ta_zzz_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 268);

    auto ta_zzz_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 269);

    auto ta_zzz_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 270);

    auto ta_zzz_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 271);

    auto ta_zzz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 272);

    auto ta_zzz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 273);

    auto ta_zzz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 274);

    auto ta_zzz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 275);

    auto ta_zzz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 276);

    auto ta_zzz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 277);

    auto ta_zzz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 278);

    auto ta_zzz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 279);

    // Set up components of auxiliary buffer : FI

    auto ta1_x_xxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi);

    auto ta1_x_xxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 1);

    auto ta1_x_xxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 2);

    auto ta1_x_xxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 3);

    auto ta1_x_xxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 4);

    auto ta1_x_xxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 5);

    auto ta1_x_xxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 6);

    auto ta1_x_xxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 7);

    auto ta1_x_xxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 8);

    auto ta1_x_xxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 9);

    auto ta1_x_xxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 10);

    auto ta1_x_xxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 11);

    auto ta1_x_xxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 12);

    auto ta1_x_xxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 13);

    auto ta1_x_xxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 14);

    auto ta1_x_xxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 15);

    auto ta1_x_xxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 16);

    auto ta1_x_xxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 17);

    auto ta1_x_xxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 18);

    auto ta1_x_xxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 19);

    auto ta1_x_xxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 20);

    auto ta1_x_xxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 21);

    auto ta1_x_xxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 22);

    auto ta1_x_xxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 23);

    auto ta1_x_xxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 24);

    auto ta1_x_xxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 25);

    auto ta1_x_xxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 26);

    auto ta1_x_xxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 27);

    auto ta1_x_xxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 28);

    auto ta1_x_xxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 29);

    auto ta1_x_xxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 30);

    auto ta1_x_xxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 31);

    auto ta1_x_xxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 32);

    auto ta1_x_xxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 33);

    auto ta1_x_xxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 34);

    auto ta1_x_xxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 35);

    auto ta1_x_xxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 36);

    auto ta1_x_xxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 37);

    auto ta1_x_xxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 38);

    auto ta1_x_xxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 39);

    auto ta1_x_xxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 40);

    auto ta1_x_xxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 41);

    auto ta1_x_xxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 42);

    auto ta1_x_xxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 43);

    auto ta1_x_xxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 44);

    auto ta1_x_xxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 45);

    auto ta1_x_xxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 46);

    auto ta1_x_xxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 47);

    auto ta1_x_xxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 48);

    auto ta1_x_xxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 49);

    auto ta1_x_xxy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 55);

    auto ta1_x_xxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 56);

    auto ta1_x_xxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 57);

    auto ta1_x_xxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 58);

    auto ta1_x_xxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 59);

    auto ta1_x_xxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 60);

    auto ta1_x_xxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 61);

    auto ta1_x_xxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 62);

    auto ta1_x_xxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 63);

    auto ta1_x_xxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 64);

    auto ta1_x_xxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 65);

    auto ta1_x_xxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 66);

    auto ta1_x_xxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 67);

    auto ta1_x_xxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 68);

    auto ta1_x_xxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 69);

    auto ta1_x_xxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 70);

    auto ta1_x_xxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 71);

    auto ta1_x_xxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 72);

    auto ta1_x_xxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 73);

    auto ta1_x_xxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 74);

    auto ta1_x_xxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 75);

    auto ta1_x_xxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 76);

    auto ta1_x_xxz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 77);

    auto ta1_x_xxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 78);

    auto ta1_x_xxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 79);

    auto ta1_x_xxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 80);

    auto ta1_x_xxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 81);

    auto ta1_x_xxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 82);

    auto ta1_x_xxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 83);

    auto ta1_x_xyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 84);

    auto ta1_x_xyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 85);

    auto ta1_x_xyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 86);

    auto ta1_x_xyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 87);

    auto ta1_x_xyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 88);

    auto ta1_x_xyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 89);

    auto ta1_x_xyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 90);

    auto ta1_x_xyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 91);

    auto ta1_x_xyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 92);

    auto ta1_x_xyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 93);

    auto ta1_x_xyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 94);

    auto ta1_x_xyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 95);

    auto ta1_x_xyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 96);

    auto ta1_x_xyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 97);

    auto ta1_x_xyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 98);

    auto ta1_x_xyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 99);

    auto ta1_x_xyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 100);

    auto ta1_x_xyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 101);

    auto ta1_x_xyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 102);

    auto ta1_x_xyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 103);

    auto ta1_x_xyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 104);

    auto ta1_x_xyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 105);

    auto ta1_x_xyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 106);

    auto ta1_x_xyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 107);

    auto ta1_x_xyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 108);

    auto ta1_x_xyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 109);

    auto ta1_x_xyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 110);

    auto ta1_x_xyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 114);

    auto ta1_x_xyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 117);

    auto ta1_x_xyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 121);

    auto ta1_x_xyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 126);

    auto ta1_x_xyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 132);

    auto ta1_x_xzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 140);

    auto ta1_x_xzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 141);

    auto ta1_x_xzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 142);

    auto ta1_x_xzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 143);

    auto ta1_x_xzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 144);

    auto ta1_x_xzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 145);

    auto ta1_x_xzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 146);

    auto ta1_x_xzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 147);

    auto ta1_x_xzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 148);

    auto ta1_x_xzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 149);

    auto ta1_x_xzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 150);

    auto ta1_x_xzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 151);

    auto ta1_x_xzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 152);

    auto ta1_x_xzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 153);

    auto ta1_x_xzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 154);

    auto ta1_x_xzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 155);

    auto ta1_x_xzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 156);

    auto ta1_x_xzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 157);

    auto ta1_x_xzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 158);

    auto ta1_x_xzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 159);

    auto ta1_x_xzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 160);

    auto ta1_x_xzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 162);

    auto ta1_x_xzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 163);

    auto ta1_x_xzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 164);

    auto ta1_x_xzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 165);

    auto ta1_x_xzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 166);

    auto ta1_x_xzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 167);

    auto ta1_x_yyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 168);

    auto ta1_x_yyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 169);

    auto ta1_x_yyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 170);

    auto ta1_x_yyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 171);

    auto ta1_x_yyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 172);

    auto ta1_x_yyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 173);

    auto ta1_x_yyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 174);

    auto ta1_x_yyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 175);

    auto ta1_x_yyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 176);

    auto ta1_x_yyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 177);

    auto ta1_x_yyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 178);

    auto ta1_x_yyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 179);

    auto ta1_x_yyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 180);

    auto ta1_x_yyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 181);

    auto ta1_x_yyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 182);

    auto ta1_x_yyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 183);

    auto ta1_x_yyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 184);

    auto ta1_x_yyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 185);

    auto ta1_x_yyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 186);

    auto ta1_x_yyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 187);

    auto ta1_x_yyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 188);

    auto ta1_x_yyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 189);

    auto ta1_x_yyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 190);

    auto ta1_x_yyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 191);

    auto ta1_x_yyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 192);

    auto ta1_x_yyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 193);

    auto ta1_x_yyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 194);

    auto ta1_x_yyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 195);

    auto ta1_x_yyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 197);

    auto ta1_x_yyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 198);

    auto ta1_x_yyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 199);

    auto ta1_x_yyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 201);

    auto ta1_x_yyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 202);

    auto ta1_x_yyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 205);

    auto ta1_x_yyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 206);

    auto ta1_x_yyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 210);

    auto ta1_x_yyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 211);

    auto ta1_x_yyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 216);

    auto ta1_x_yyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 217);

    auto ta1_x_yyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 218);

    auto ta1_x_yyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 219);

    auto ta1_x_yyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 220);

    auto ta1_x_yyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 221);

    auto ta1_x_yyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 222);

    auto ta1_x_yyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 223);

    auto ta1_x_yzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 224);

    auto ta1_x_yzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 226);

    auto ta1_x_yzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 228);

    auto ta1_x_yzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 229);

    auto ta1_x_yzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 231);

    auto ta1_x_yzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 232);

    auto ta1_x_yzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 233);

    auto ta1_x_yzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 235);

    auto ta1_x_yzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 236);

    auto ta1_x_yzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 237);

    auto ta1_x_yzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 238);

    auto ta1_x_yzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 240);

    auto ta1_x_yzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 241);

    auto ta1_x_yzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 242);

    auto ta1_x_yzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 243);

    auto ta1_x_yzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 244);

    auto ta1_x_yzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 245);

    auto ta1_x_yzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 246);

    auto ta1_x_yzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 247);

    auto ta1_x_yzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 248);

    auto ta1_x_yzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 249);

    auto ta1_x_yzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 250);

    auto ta1_x_yzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 251);

    auto ta1_x_zzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 252);

    auto ta1_x_zzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 253);

    auto ta1_x_zzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 254);

    auto ta1_x_zzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 255);

    auto ta1_x_zzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 256);

    auto ta1_x_zzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 257);

    auto ta1_x_zzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 258);

    auto ta1_x_zzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 259);

    auto ta1_x_zzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 260);

    auto ta1_x_zzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 261);

    auto ta1_x_zzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 262);

    auto ta1_x_zzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 263);

    auto ta1_x_zzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 264);

    auto ta1_x_zzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 265);

    auto ta1_x_zzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 266);

    auto ta1_x_zzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 267);

    auto ta1_x_zzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 268);

    auto ta1_x_zzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 269);

    auto ta1_x_zzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 270);

    auto ta1_x_zzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 271);

    auto ta1_x_zzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 272);

    auto ta1_x_zzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 273);

    auto ta1_x_zzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 274);

    auto ta1_x_zzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 275);

    auto ta1_x_zzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 276);

    auto ta1_x_zzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 277);

    auto ta1_x_zzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 278);

    auto ta1_x_zzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 279);

    auto ta1_y_xxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 280);

    auto ta1_y_xxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 281);

    auto ta1_y_xxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 282);

    auto ta1_y_xxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 283);

    auto ta1_y_xxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 284);

    auto ta1_y_xxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 285);

    auto ta1_y_xxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 286);

    auto ta1_y_xxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 287);

    auto ta1_y_xxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 288);

    auto ta1_y_xxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 289);

    auto ta1_y_xxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 290);

    auto ta1_y_xxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 291);

    auto ta1_y_xxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 292);

    auto ta1_y_xxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 293);

    auto ta1_y_xxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 294);

    auto ta1_y_xxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 295);

    auto ta1_y_xxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 296);

    auto ta1_y_xxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 297);

    auto ta1_y_xxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 298);

    auto ta1_y_xxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 299);

    auto ta1_y_xxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 300);

    auto ta1_y_xxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 301);

    auto ta1_y_xxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 302);

    auto ta1_y_xxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 303);

    auto ta1_y_xxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 304);

    auto ta1_y_xxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 305);

    auto ta1_y_xxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 306);

    auto ta1_y_xxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 307);

    auto ta1_y_xxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 308);

    auto ta1_y_xxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 309);

    auto ta1_y_xxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 310);

    auto ta1_y_xxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 311);

    auto ta1_y_xxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 312);

    auto ta1_y_xxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 313);

    auto ta1_y_xxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 314);

    auto ta1_y_xxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 315);

    auto ta1_y_xxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 316);

    auto ta1_y_xxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 317);

    auto ta1_y_xxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 318);

    auto ta1_y_xxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 319);

    auto ta1_y_xxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 320);

    auto ta1_y_xxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 321);

    auto ta1_y_xxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 322);

    auto ta1_y_xxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 323);

    auto ta1_y_xxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 324);

    auto ta1_y_xxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 325);

    auto ta1_y_xxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 326);

    auto ta1_y_xxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 327);

    auto ta1_y_xxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 328);

    auto ta1_y_xxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 329);

    auto ta1_y_xxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 330);

    auto ta1_y_xxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 331);

    auto ta1_y_xxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 332);

    auto ta1_y_xxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 333);

    auto ta1_y_xxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 334);

    auto ta1_y_xxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 336);

    auto ta1_y_xxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 337);

    auto ta1_y_xxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 338);

    auto ta1_y_xxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 339);

    auto ta1_y_xxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 341);

    auto ta1_y_xxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 342);

    auto ta1_y_xxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 345);

    auto ta1_y_xxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 346);

    auto ta1_y_xxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 350);

    auto ta1_y_xxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 351);

    auto ta1_y_xxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 356);

    auto ta1_y_xxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 358);

    auto ta1_y_xxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 359);

    auto ta1_y_xxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 360);

    auto ta1_y_xxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 361);

    auto ta1_y_xxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 362);

    auto ta1_y_xxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 363);

    auto ta1_y_xyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 364);

    auto ta1_y_xyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 365);

    auto ta1_y_xyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 367);

    auto ta1_y_xyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 368);

    auto ta1_y_xyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 370);

    auto ta1_y_xyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 371);

    auto ta1_y_xyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 372);

    auto ta1_y_xyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 374);

    auto ta1_y_xyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 375);

    auto ta1_y_xyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 376);

    auto ta1_y_xyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 377);

    auto ta1_y_xyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 379);

    auto ta1_y_xyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 380);

    auto ta1_y_xyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 381);

    auto ta1_y_xyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 382);

    auto ta1_y_xyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 383);

    auto ta1_y_xyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 385);

    auto ta1_y_xyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 386);

    auto ta1_y_xyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 387);

    auto ta1_y_xyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 388);

    auto ta1_y_xyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 389);

    auto ta1_y_xyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 390);

    auto ta1_y_xyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 391);

    auto ta1_y_xyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 414);

    auto ta1_y_xyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 415);

    auto ta1_y_xyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 416);

    auto ta1_y_xyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 417);

    auto ta1_y_xyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 418);

    auto ta1_y_xzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 420);

    auto ta1_y_xzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 422);

    auto ta1_y_xzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 424);

    auto ta1_y_xzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 425);

    auto ta1_y_xzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 427);

    auto ta1_y_xzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 428);

    auto ta1_y_xzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 429);

    auto ta1_y_xzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 431);

    auto ta1_y_xzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 432);

    auto ta1_y_xzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 433);

    auto ta1_y_xzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 434);

    auto ta1_y_xzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 436);

    auto ta1_y_xzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 437);

    auto ta1_y_xzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 438);

    auto ta1_y_xzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 439);

    auto ta1_y_xzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 440);

    auto ta1_y_xzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 441);

    auto ta1_y_xzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 442);

    auto ta1_y_xzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 443);

    auto ta1_y_xzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 444);

    auto ta1_y_xzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 445);

    auto ta1_y_xzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 446);

    auto ta1_y_xzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 447);

    auto ta1_y_yyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 448);

    auto ta1_y_yyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 449);

    auto ta1_y_yyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 450);

    auto ta1_y_yyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 451);

    auto ta1_y_yyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 452);

    auto ta1_y_yyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 453);

    auto ta1_y_yyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 454);

    auto ta1_y_yyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 455);

    auto ta1_y_yyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 456);

    auto ta1_y_yyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 457);

    auto ta1_y_yyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 458);

    auto ta1_y_yyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 459);

    auto ta1_y_yyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 460);

    auto ta1_y_yyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 461);

    auto ta1_y_yyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 462);

    auto ta1_y_yyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 463);

    auto ta1_y_yyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 464);

    auto ta1_y_yyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 465);

    auto ta1_y_yyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 466);

    auto ta1_y_yyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 467);

    auto ta1_y_yyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 468);

    auto ta1_y_yyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 469);

    auto ta1_y_yyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 470);

    auto ta1_y_yyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 471);

    auto ta1_y_yyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 472);

    auto ta1_y_yyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 473);

    auto ta1_y_yyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 474);

    auto ta1_y_yyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 475);

    auto ta1_y_yyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 476);

    auto ta1_y_yyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 477);

    auto ta1_y_yyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 478);

    auto ta1_y_yyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 479);

    auto ta1_y_yyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 480);

    auto ta1_y_yyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 481);

    auto ta1_y_yyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 482);

    auto ta1_y_yyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 483);

    auto ta1_y_yyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 484);

    auto ta1_y_yyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 485);

    auto ta1_y_yyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 486);

    auto ta1_y_yyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 487);

    auto ta1_y_yyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 488);

    auto ta1_y_yyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 489);

    auto ta1_y_yyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 490);

    auto ta1_y_yyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 491);

    auto ta1_y_yyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 492);

    auto ta1_y_yyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 493);

    auto ta1_y_yyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 494);

    auto ta1_y_yyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 495);

    auto ta1_y_yyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 496);

    auto ta1_y_yyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 497);

    auto ta1_y_yyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 498);

    auto ta1_y_yyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 499);

    auto ta1_y_yyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 500);

    auto ta1_y_yyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 501);

    auto ta1_y_yyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 502);

    auto ta1_y_yyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 503);

    auto ta1_y_yzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 505);

    auto ta1_y_yzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 506);

    auto ta1_y_yzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 507);

    auto ta1_y_yzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 508);

    auto ta1_y_yzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 509);

    auto ta1_y_yzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 510);

    auto ta1_y_yzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 511);

    auto ta1_y_yzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 512);

    auto ta1_y_yzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 513);

    auto ta1_y_yzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 514);

    auto ta1_y_yzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 515);

    auto ta1_y_yzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 516);

    auto ta1_y_yzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 517);

    auto ta1_y_yzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 518);

    auto ta1_y_yzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 519);

    auto ta1_y_yzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 520);

    auto ta1_y_yzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 521);

    auto ta1_y_yzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 522);

    auto ta1_y_yzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 523);

    auto ta1_y_yzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 524);

    auto ta1_y_yzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 525);

    auto ta1_y_yzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 526);

    auto ta1_y_yzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 527);

    auto ta1_y_yzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 528);

    auto ta1_y_yzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 529);

    auto ta1_y_yzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 530);

    auto ta1_y_yzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 531);

    auto ta1_y_zzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 532);

    auto ta1_y_zzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 533);

    auto ta1_y_zzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 534);

    auto ta1_y_zzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 535);

    auto ta1_y_zzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 536);

    auto ta1_y_zzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 537);

    auto ta1_y_zzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 538);

    auto ta1_y_zzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 539);

    auto ta1_y_zzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 540);

    auto ta1_y_zzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 541);

    auto ta1_y_zzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 542);

    auto ta1_y_zzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 543);

    auto ta1_y_zzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 544);

    auto ta1_y_zzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 545);

    auto ta1_y_zzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 546);

    auto ta1_y_zzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 547);

    auto ta1_y_zzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 548);

    auto ta1_y_zzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 549);

    auto ta1_y_zzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 550);

    auto ta1_y_zzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 551);

    auto ta1_y_zzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 552);

    auto ta1_y_zzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 553);

    auto ta1_y_zzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 554);

    auto ta1_y_zzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 555);

    auto ta1_y_zzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 556);

    auto ta1_y_zzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 557);

    auto ta1_y_zzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 558);

    auto ta1_y_zzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 559);

    auto ta1_z_xxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 560);

    auto ta1_z_xxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 561);

    auto ta1_z_xxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 562);

    auto ta1_z_xxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 563);

    auto ta1_z_xxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 564);

    auto ta1_z_xxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 565);

    auto ta1_z_xxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 566);

    auto ta1_z_xxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 567);

    auto ta1_z_xxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 568);

    auto ta1_z_xxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 569);

    auto ta1_z_xxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 570);

    auto ta1_z_xxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 571);

    auto ta1_z_xxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 572);

    auto ta1_z_xxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 573);

    auto ta1_z_xxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 574);

    auto ta1_z_xxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 575);

    auto ta1_z_xxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 576);

    auto ta1_z_xxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 577);

    auto ta1_z_xxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 578);

    auto ta1_z_xxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 579);

    auto ta1_z_xxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 580);

    auto ta1_z_xxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 581);

    auto ta1_z_xxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 582);

    auto ta1_z_xxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 583);

    auto ta1_z_xxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 584);

    auto ta1_z_xxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 585);

    auto ta1_z_xxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 586);

    auto ta1_z_xxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 587);

    auto ta1_z_xxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 588);

    auto ta1_z_xxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 589);

    auto ta1_z_xxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 590);

    auto ta1_z_xxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 591);

    auto ta1_z_xxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 593);

    auto ta1_z_xxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 594);

    auto ta1_z_xxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 597);

    auto ta1_z_xxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 598);

    auto ta1_z_xxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 602);

    auto ta1_z_xxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 603);

    auto ta1_z_xxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 608);

    auto ta1_z_xxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 609);

    auto ta1_z_xxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 610);

    auto ta1_z_xxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 611);

    auto ta1_z_xxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 612);

    auto ta1_z_xxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 613);

    auto ta1_z_xxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 614);

    auto ta1_z_xxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 616);

    auto ta1_z_xxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 617);

    auto ta1_z_xxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 618);

    auto ta1_z_xxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 619);

    auto ta1_z_xxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 620);

    auto ta1_z_xxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 621);

    auto ta1_z_xxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 622);

    auto ta1_z_xxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 623);

    auto ta1_z_xxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 624);

    auto ta1_z_xxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 625);

    auto ta1_z_xxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 626);

    auto ta1_z_xxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 627);

    auto ta1_z_xxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 628);

    auto ta1_z_xxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 629);

    auto ta1_z_xxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 630);

    auto ta1_z_xxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 631);

    auto ta1_z_xxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 632);

    auto ta1_z_xxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 633);

    auto ta1_z_xxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 634);

    auto ta1_z_xxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 635);

    auto ta1_z_xxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 636);

    auto ta1_z_xxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 638);

    auto ta1_z_xxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 639);

    auto ta1_z_xxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 640);

    auto ta1_z_xxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 641);

    auto ta1_z_xxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 642);

    auto ta1_z_xxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 643);

    auto ta1_z_xyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 644);

    auto ta1_z_xyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 645);

    auto ta1_z_xyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 647);

    auto ta1_z_xyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 648);

    auto ta1_z_xyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 650);

    auto ta1_z_xyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 651);

    auto ta1_z_xyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 652);

    auto ta1_z_xyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 654);

    auto ta1_z_xyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 655);

    auto ta1_z_xyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 656);

    auto ta1_z_xyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 657);

    auto ta1_z_xyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 659);

    auto ta1_z_xyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 660);

    auto ta1_z_xyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 661);

    auto ta1_z_xyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 662);

    auto ta1_z_xyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 663);

    auto ta1_z_xyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 665);

    auto ta1_z_xyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 666);

    auto ta1_z_xyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 667);

    auto ta1_z_xyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 668);

    auto ta1_z_xyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 669);

    auto ta1_z_xyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 670);

    auto ta1_z_xyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 671);

    auto ta1_z_xyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 694);

    auto ta1_z_xyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 695);

    auto ta1_z_xyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 696);

    auto ta1_z_xyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 697);

    auto ta1_z_xyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 698);

    auto ta1_z_xzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 700);

    auto ta1_z_xzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 702);

    auto ta1_z_xzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 704);

    auto ta1_z_xzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 705);

    auto ta1_z_xzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 707);

    auto ta1_z_xzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 708);

    auto ta1_z_xzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 709);

    auto ta1_z_xzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 711);

    auto ta1_z_xzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 712);

    auto ta1_z_xzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 713);

    auto ta1_z_xzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 714);

    auto ta1_z_xzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 716);

    auto ta1_z_xzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 717);

    auto ta1_z_xzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 718);

    auto ta1_z_xzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 719);

    auto ta1_z_xzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 720);

    auto ta1_z_xzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 721);

    auto ta1_z_xzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 722);

    auto ta1_z_xzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 723);

    auto ta1_z_xzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 724);

    auto ta1_z_xzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 725);

    auto ta1_z_xzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 726);

    auto ta1_z_xzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 727);

    auto ta1_z_yyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 728);

    auto ta1_z_yyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 729);

    auto ta1_z_yyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 730);

    auto ta1_z_yyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 731);

    auto ta1_z_yyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 732);

    auto ta1_z_yyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 733);

    auto ta1_z_yyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 734);

    auto ta1_z_yyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 735);

    auto ta1_z_yyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 736);

    auto ta1_z_yyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 737);

    auto ta1_z_yyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 738);

    auto ta1_z_yyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 739);

    auto ta1_z_yyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 740);

    auto ta1_z_yyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 741);

    auto ta1_z_yyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 742);

    auto ta1_z_yyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 743);

    auto ta1_z_yyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 744);

    auto ta1_z_yyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 745);

    auto ta1_z_yyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 746);

    auto ta1_z_yyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 747);

    auto ta1_z_yyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 748);

    auto ta1_z_yyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 749);

    auto ta1_z_yyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 750);

    auto ta1_z_yyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 751);

    auto ta1_z_yyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 752);

    auto ta1_z_yyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 753);

    auto ta1_z_yyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 754);

    auto ta1_z_yyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 755);

    auto ta1_z_yyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 757);

    auto ta1_z_yyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 758);

    auto ta1_z_yyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 759);

    auto ta1_z_yyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 760);

    auto ta1_z_yyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 761);

    auto ta1_z_yyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 762);

    auto ta1_z_yyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 763);

    auto ta1_z_yyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 764);

    auto ta1_z_yyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 765);

    auto ta1_z_yyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 766);

    auto ta1_z_yyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 767);

    auto ta1_z_yyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 768);

    auto ta1_z_yyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 769);

    auto ta1_z_yyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 770);

    auto ta1_z_yyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 771);

    auto ta1_z_yyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 772);

    auto ta1_z_yyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 773);

    auto ta1_z_yyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 774);

    auto ta1_z_yyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 775);

    auto ta1_z_yyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 776);

    auto ta1_z_yyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 777);

    auto ta1_z_yyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 778);

    auto ta1_z_yyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 779);

    auto ta1_z_yyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 780);

    auto ta1_z_yyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 781);

    auto ta1_z_yyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 782);

    auto ta1_z_yyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 783);

    auto ta1_z_yzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 784);

    auto ta1_z_yzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 785);

    auto ta1_z_yzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 786);

    auto ta1_z_yzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 787);

    auto ta1_z_yzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 788);

    auto ta1_z_yzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 789);

    auto ta1_z_yzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 790);

    auto ta1_z_yzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 791);

    auto ta1_z_yzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 792);

    auto ta1_z_yzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 793);

    auto ta1_z_yzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 794);

    auto ta1_z_yzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 795);

    auto ta1_z_yzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 796);

    auto ta1_z_yzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 797);

    auto ta1_z_yzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 798);

    auto ta1_z_yzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 799);

    auto ta1_z_yzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 800);

    auto ta1_z_yzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 801);

    auto ta1_z_yzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 802);

    auto ta1_z_yzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 803);

    auto ta1_z_yzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 804);

    auto ta1_z_yzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 805);

    auto ta1_z_yzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 806);

    auto ta1_z_yzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 807);

    auto ta1_z_yzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 808);

    auto ta1_z_yzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 809);

    auto ta1_z_yzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 810);

    auto ta1_z_yzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 811);

    auto ta1_z_zzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 812);

    auto ta1_z_zzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 813);

    auto ta1_z_zzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 814);

    auto ta1_z_zzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 815);

    auto ta1_z_zzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 816);

    auto ta1_z_zzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 817);

    auto ta1_z_zzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 818);

    auto ta1_z_zzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 819);

    auto ta1_z_zzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 820);

    auto ta1_z_zzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 821);

    auto ta1_z_zzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 822);

    auto ta1_z_zzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 823);

    auto ta1_z_zzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 824);

    auto ta1_z_zzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 825);

    auto ta1_z_zzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 826);

    auto ta1_z_zzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 827);

    auto ta1_z_zzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 828);

    auto ta1_z_zzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 829);

    auto ta1_z_zzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 830);

    auto ta1_z_zzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 831);

    auto ta1_z_zzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 832);

    auto ta1_z_zzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 833);

    auto ta1_z_zzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 834);

    auto ta1_z_zzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 835);

    auto ta1_z_zzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 836);

    auto ta1_z_zzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 837);

    auto ta1_z_zzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 838);

    auto ta1_z_zzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 839);

    // Set up components of auxiliary buffer : FI

    auto ta1_x_xxx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi);

    auto ta1_x_xxx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 1);

    auto ta1_x_xxx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 2);

    auto ta1_x_xxx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 3);

    auto ta1_x_xxx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 4);

    auto ta1_x_xxx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 5);

    auto ta1_x_xxx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 6);

    auto ta1_x_xxx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 7);

    auto ta1_x_xxx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 8);

    auto ta1_x_xxx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 9);

    auto ta1_x_xxx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 10);

    auto ta1_x_xxx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 11);

    auto ta1_x_xxx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 12);

    auto ta1_x_xxx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 13);

    auto ta1_x_xxx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 14);

    auto ta1_x_xxx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 15);

    auto ta1_x_xxx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 16);

    auto ta1_x_xxx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 17);

    auto ta1_x_xxx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 18);

    auto ta1_x_xxx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 19);

    auto ta1_x_xxx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 20);

    auto ta1_x_xxx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 21);

    auto ta1_x_xxx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 22);

    auto ta1_x_xxx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 23);

    auto ta1_x_xxx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 24);

    auto ta1_x_xxx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 25);

    auto ta1_x_xxx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 26);

    auto ta1_x_xxx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 27);

    auto ta1_x_xxy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 28);

    auto ta1_x_xxy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 29);

    auto ta1_x_xxy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 30);

    auto ta1_x_xxy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 31);

    auto ta1_x_xxy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 32);

    auto ta1_x_xxy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 33);

    auto ta1_x_xxy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 34);

    auto ta1_x_xxy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 35);

    auto ta1_x_xxy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 36);

    auto ta1_x_xxy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 37);

    auto ta1_x_xxy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 38);

    auto ta1_x_xxy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 39);

    auto ta1_x_xxy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 40);

    auto ta1_x_xxy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 41);

    auto ta1_x_xxy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 42);

    auto ta1_x_xxy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 43);

    auto ta1_x_xxy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 44);

    auto ta1_x_xxy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 45);

    auto ta1_x_xxy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 46);

    auto ta1_x_xxy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 47);

    auto ta1_x_xxy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 48);

    auto ta1_x_xxy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 49);

    auto ta1_x_xxy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 55);

    auto ta1_x_xxz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 56);

    auto ta1_x_xxz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 57);

    auto ta1_x_xxz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 58);

    auto ta1_x_xxz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 59);

    auto ta1_x_xxz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 60);

    auto ta1_x_xxz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 61);

    auto ta1_x_xxz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 62);

    auto ta1_x_xxz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 63);

    auto ta1_x_xxz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 64);

    auto ta1_x_xxz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 65);

    auto ta1_x_xxz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 66);

    auto ta1_x_xxz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 67);

    auto ta1_x_xxz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 68);

    auto ta1_x_xxz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 69);

    auto ta1_x_xxz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 70);

    auto ta1_x_xxz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 71);

    auto ta1_x_xxz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 72);

    auto ta1_x_xxz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 73);

    auto ta1_x_xxz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 74);

    auto ta1_x_xxz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 75);

    auto ta1_x_xxz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 76);

    auto ta1_x_xxz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 77);

    auto ta1_x_xxz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 78);

    auto ta1_x_xxz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 79);

    auto ta1_x_xxz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 80);

    auto ta1_x_xxz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 81);

    auto ta1_x_xxz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 82);

    auto ta1_x_xxz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 83);

    auto ta1_x_xyy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 84);

    auto ta1_x_xyy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 85);

    auto ta1_x_xyy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 86);

    auto ta1_x_xyy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 87);

    auto ta1_x_xyy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 88);

    auto ta1_x_xyy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 89);

    auto ta1_x_xyy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 90);

    auto ta1_x_xyy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 91);

    auto ta1_x_xyy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 92);

    auto ta1_x_xyy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 93);

    auto ta1_x_xyy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 94);

    auto ta1_x_xyy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 95);

    auto ta1_x_xyy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 96);

    auto ta1_x_xyy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 97);

    auto ta1_x_xyy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 98);

    auto ta1_x_xyy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 99);

    auto ta1_x_xyy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 100);

    auto ta1_x_xyy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 101);

    auto ta1_x_xyy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 102);

    auto ta1_x_xyy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 103);

    auto ta1_x_xyy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 104);

    auto ta1_x_xyy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 105);

    auto ta1_x_xyy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 106);

    auto ta1_x_xyy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 107);

    auto ta1_x_xyy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 108);

    auto ta1_x_xyy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 109);

    auto ta1_x_xyy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 110);

    auto ta1_x_xyz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 114);

    auto ta1_x_xyz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 117);

    auto ta1_x_xyz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 121);

    auto ta1_x_xyz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 126);

    auto ta1_x_xyz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 132);

    auto ta1_x_xzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 140);

    auto ta1_x_xzz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 141);

    auto ta1_x_xzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 142);

    auto ta1_x_xzz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 143);

    auto ta1_x_xzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 144);

    auto ta1_x_xzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 145);

    auto ta1_x_xzz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 146);

    auto ta1_x_xzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 147);

    auto ta1_x_xzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 148);

    auto ta1_x_xzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 149);

    auto ta1_x_xzz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 150);

    auto ta1_x_xzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 151);

    auto ta1_x_xzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 152);

    auto ta1_x_xzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 153);

    auto ta1_x_xzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 154);

    auto ta1_x_xzz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 155);

    auto ta1_x_xzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 156);

    auto ta1_x_xzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 157);

    auto ta1_x_xzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 158);

    auto ta1_x_xzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 159);

    auto ta1_x_xzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 160);

    auto ta1_x_xzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 162);

    auto ta1_x_xzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 163);

    auto ta1_x_xzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 164);

    auto ta1_x_xzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 165);

    auto ta1_x_xzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 166);

    auto ta1_x_xzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 167);

    auto ta1_x_yyy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 168);

    auto ta1_x_yyy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 169);

    auto ta1_x_yyy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 170);

    auto ta1_x_yyy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 171);

    auto ta1_x_yyy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 172);

    auto ta1_x_yyy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 173);

    auto ta1_x_yyy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 174);

    auto ta1_x_yyy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 175);

    auto ta1_x_yyy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 176);

    auto ta1_x_yyy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 177);

    auto ta1_x_yyy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 178);

    auto ta1_x_yyy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 179);

    auto ta1_x_yyy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 180);

    auto ta1_x_yyy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 181);

    auto ta1_x_yyy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 182);

    auto ta1_x_yyy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 183);

    auto ta1_x_yyy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 184);

    auto ta1_x_yyy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 185);

    auto ta1_x_yyy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 186);

    auto ta1_x_yyy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 187);

    auto ta1_x_yyy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 188);

    auto ta1_x_yyy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 189);

    auto ta1_x_yyy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 190);

    auto ta1_x_yyy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 191);

    auto ta1_x_yyy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 192);

    auto ta1_x_yyy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 193);

    auto ta1_x_yyy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 194);

    auto ta1_x_yyy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 195);

    auto ta1_x_yyz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 197);

    auto ta1_x_yyz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 198);

    auto ta1_x_yyz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 199);

    auto ta1_x_yyz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 201);

    auto ta1_x_yyz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 202);

    auto ta1_x_yyz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 205);

    auto ta1_x_yyz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 206);

    auto ta1_x_yyz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 210);

    auto ta1_x_yyz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 211);

    auto ta1_x_yyz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 216);

    auto ta1_x_yyz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 217);

    auto ta1_x_yyz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 218);

    auto ta1_x_yyz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 219);

    auto ta1_x_yyz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 220);

    auto ta1_x_yyz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 221);

    auto ta1_x_yyz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 222);

    auto ta1_x_yyz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 223);

    auto ta1_x_yzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 224);

    auto ta1_x_yzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 226);

    auto ta1_x_yzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 228);

    auto ta1_x_yzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 229);

    auto ta1_x_yzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 231);

    auto ta1_x_yzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 232);

    auto ta1_x_yzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 233);

    auto ta1_x_yzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 235);

    auto ta1_x_yzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 236);

    auto ta1_x_yzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 237);

    auto ta1_x_yzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 238);

    auto ta1_x_yzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 240);

    auto ta1_x_yzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 241);

    auto ta1_x_yzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 242);

    auto ta1_x_yzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 243);

    auto ta1_x_yzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 244);

    auto ta1_x_yzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 245);

    auto ta1_x_yzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 246);

    auto ta1_x_yzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 247);

    auto ta1_x_yzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 248);

    auto ta1_x_yzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 249);

    auto ta1_x_yzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 250);

    auto ta1_x_yzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 251);

    auto ta1_x_zzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 252);

    auto ta1_x_zzz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 253);

    auto ta1_x_zzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 254);

    auto ta1_x_zzz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 255);

    auto ta1_x_zzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 256);

    auto ta1_x_zzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 257);

    auto ta1_x_zzz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 258);

    auto ta1_x_zzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 259);

    auto ta1_x_zzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 260);

    auto ta1_x_zzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 261);

    auto ta1_x_zzz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 262);

    auto ta1_x_zzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 263);

    auto ta1_x_zzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 264);

    auto ta1_x_zzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 265);

    auto ta1_x_zzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 266);

    auto ta1_x_zzz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 267);

    auto ta1_x_zzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 268);

    auto ta1_x_zzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 269);

    auto ta1_x_zzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 270);

    auto ta1_x_zzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 271);

    auto ta1_x_zzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 272);

    auto ta1_x_zzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 273);

    auto ta1_x_zzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 274);

    auto ta1_x_zzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 275);

    auto ta1_x_zzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 276);

    auto ta1_x_zzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 277);

    auto ta1_x_zzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 278);

    auto ta1_x_zzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 279);

    auto ta1_y_xxx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 280);

    auto ta1_y_xxx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 281);

    auto ta1_y_xxx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 282);

    auto ta1_y_xxx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 283);

    auto ta1_y_xxx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 284);

    auto ta1_y_xxx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 285);

    auto ta1_y_xxx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 286);

    auto ta1_y_xxx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 287);

    auto ta1_y_xxx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 288);

    auto ta1_y_xxx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 289);

    auto ta1_y_xxx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 290);

    auto ta1_y_xxx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 291);

    auto ta1_y_xxx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 292);

    auto ta1_y_xxx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 293);

    auto ta1_y_xxx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 294);

    auto ta1_y_xxx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 295);

    auto ta1_y_xxx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 296);

    auto ta1_y_xxx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 297);

    auto ta1_y_xxx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 298);

    auto ta1_y_xxx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 299);

    auto ta1_y_xxx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 300);

    auto ta1_y_xxx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 301);

    auto ta1_y_xxx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 302);

    auto ta1_y_xxx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 303);

    auto ta1_y_xxx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 304);

    auto ta1_y_xxx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 305);

    auto ta1_y_xxx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 306);

    auto ta1_y_xxx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 307);

    auto ta1_y_xxy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 308);

    auto ta1_y_xxy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 309);

    auto ta1_y_xxy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 310);

    auto ta1_y_xxy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 311);

    auto ta1_y_xxy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 312);

    auto ta1_y_xxy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 313);

    auto ta1_y_xxy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 314);

    auto ta1_y_xxy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 315);

    auto ta1_y_xxy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 316);

    auto ta1_y_xxy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 317);

    auto ta1_y_xxy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 318);

    auto ta1_y_xxy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 319);

    auto ta1_y_xxy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 320);

    auto ta1_y_xxy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 321);

    auto ta1_y_xxy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 322);

    auto ta1_y_xxy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 323);

    auto ta1_y_xxy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 324);

    auto ta1_y_xxy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 325);

    auto ta1_y_xxy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 326);

    auto ta1_y_xxy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 327);

    auto ta1_y_xxy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 328);

    auto ta1_y_xxy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 329);

    auto ta1_y_xxy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 330);

    auto ta1_y_xxy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 331);

    auto ta1_y_xxy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 332);

    auto ta1_y_xxy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 333);

    auto ta1_y_xxy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 334);

    auto ta1_y_xxz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 336);

    auto ta1_y_xxz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 337);

    auto ta1_y_xxz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 338);

    auto ta1_y_xxz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 339);

    auto ta1_y_xxz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 341);

    auto ta1_y_xxz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 342);

    auto ta1_y_xxz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 345);

    auto ta1_y_xxz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 346);

    auto ta1_y_xxz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 350);

    auto ta1_y_xxz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 351);

    auto ta1_y_xxz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 356);

    auto ta1_y_xxz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 358);

    auto ta1_y_xxz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 359);

    auto ta1_y_xxz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 360);

    auto ta1_y_xxz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 361);

    auto ta1_y_xxz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 362);

    auto ta1_y_xxz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 363);

    auto ta1_y_xyy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 364);

    auto ta1_y_xyy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 365);

    auto ta1_y_xyy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 367);

    auto ta1_y_xyy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 368);

    auto ta1_y_xyy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 370);

    auto ta1_y_xyy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 371);

    auto ta1_y_xyy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 372);

    auto ta1_y_xyy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 374);

    auto ta1_y_xyy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 375);

    auto ta1_y_xyy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 376);

    auto ta1_y_xyy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 377);

    auto ta1_y_xyy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 379);

    auto ta1_y_xyy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 380);

    auto ta1_y_xyy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 381);

    auto ta1_y_xyy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 382);

    auto ta1_y_xyy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 383);

    auto ta1_y_xyy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 385);

    auto ta1_y_xyy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 386);

    auto ta1_y_xyy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 387);

    auto ta1_y_xyy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 388);

    auto ta1_y_xyy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 389);

    auto ta1_y_xyy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 390);

    auto ta1_y_xyy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 391);

    auto ta1_y_xyz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 414);

    auto ta1_y_xyz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 415);

    auto ta1_y_xyz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 416);

    auto ta1_y_xyz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 417);

    auto ta1_y_xyz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 418);

    auto ta1_y_xzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 420);

    auto ta1_y_xzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 422);

    auto ta1_y_xzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 424);

    auto ta1_y_xzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 425);

    auto ta1_y_xzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 427);

    auto ta1_y_xzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 428);

    auto ta1_y_xzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 429);

    auto ta1_y_xzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 431);

    auto ta1_y_xzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 432);

    auto ta1_y_xzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 433);

    auto ta1_y_xzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 434);

    auto ta1_y_xzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 436);

    auto ta1_y_xzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 437);

    auto ta1_y_xzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 438);

    auto ta1_y_xzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 439);

    auto ta1_y_xzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 440);

    auto ta1_y_xzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 441);

    auto ta1_y_xzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 442);

    auto ta1_y_xzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 443);

    auto ta1_y_xzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 444);

    auto ta1_y_xzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 445);

    auto ta1_y_xzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 446);

    auto ta1_y_xzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 447);

    auto ta1_y_yyy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 448);

    auto ta1_y_yyy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 449);

    auto ta1_y_yyy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 450);

    auto ta1_y_yyy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 451);

    auto ta1_y_yyy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 452);

    auto ta1_y_yyy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 453);

    auto ta1_y_yyy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 454);

    auto ta1_y_yyy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 455);

    auto ta1_y_yyy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 456);

    auto ta1_y_yyy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 457);

    auto ta1_y_yyy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 458);

    auto ta1_y_yyy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 459);

    auto ta1_y_yyy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 460);

    auto ta1_y_yyy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 461);

    auto ta1_y_yyy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 462);

    auto ta1_y_yyy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 463);

    auto ta1_y_yyy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 464);

    auto ta1_y_yyy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 465);

    auto ta1_y_yyy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 466);

    auto ta1_y_yyy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 467);

    auto ta1_y_yyy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 468);

    auto ta1_y_yyy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 469);

    auto ta1_y_yyy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 470);

    auto ta1_y_yyy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 471);

    auto ta1_y_yyy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 472);

    auto ta1_y_yyy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 473);

    auto ta1_y_yyy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 474);

    auto ta1_y_yyy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 475);

    auto ta1_y_yyz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 476);

    auto ta1_y_yyz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 477);

    auto ta1_y_yyz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 478);

    auto ta1_y_yyz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 479);

    auto ta1_y_yyz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 480);

    auto ta1_y_yyz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 481);

    auto ta1_y_yyz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 482);

    auto ta1_y_yyz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 483);

    auto ta1_y_yyz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 484);

    auto ta1_y_yyz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 485);

    auto ta1_y_yyz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 486);

    auto ta1_y_yyz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 487);

    auto ta1_y_yyz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 488);

    auto ta1_y_yyz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 489);

    auto ta1_y_yyz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 490);

    auto ta1_y_yyz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 491);

    auto ta1_y_yyz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 492);

    auto ta1_y_yyz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 493);

    auto ta1_y_yyz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 494);

    auto ta1_y_yyz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 495);

    auto ta1_y_yyz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 496);

    auto ta1_y_yyz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 497);

    auto ta1_y_yyz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 498);

    auto ta1_y_yyz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 499);

    auto ta1_y_yyz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 500);

    auto ta1_y_yyz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 501);

    auto ta1_y_yyz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 502);

    auto ta1_y_yyz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 503);

    auto ta1_y_yzz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 505);

    auto ta1_y_yzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 506);

    auto ta1_y_yzz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 507);

    auto ta1_y_yzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 508);

    auto ta1_y_yzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 509);

    auto ta1_y_yzz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 510);

    auto ta1_y_yzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 511);

    auto ta1_y_yzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 512);

    auto ta1_y_yzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 513);

    auto ta1_y_yzz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 514);

    auto ta1_y_yzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 515);

    auto ta1_y_yzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 516);

    auto ta1_y_yzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 517);

    auto ta1_y_yzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 518);

    auto ta1_y_yzz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 519);

    auto ta1_y_yzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 520);

    auto ta1_y_yzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 521);

    auto ta1_y_yzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 522);

    auto ta1_y_yzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 523);

    auto ta1_y_yzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 524);

    auto ta1_y_yzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 525);

    auto ta1_y_yzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 526);

    auto ta1_y_yzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 527);

    auto ta1_y_yzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 528);

    auto ta1_y_yzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 529);

    auto ta1_y_yzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 530);

    auto ta1_y_yzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 531);

    auto ta1_y_zzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 532);

    auto ta1_y_zzz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 533);

    auto ta1_y_zzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 534);

    auto ta1_y_zzz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 535);

    auto ta1_y_zzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 536);

    auto ta1_y_zzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 537);

    auto ta1_y_zzz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 538);

    auto ta1_y_zzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 539);

    auto ta1_y_zzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 540);

    auto ta1_y_zzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 541);

    auto ta1_y_zzz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 542);

    auto ta1_y_zzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 543);

    auto ta1_y_zzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 544);

    auto ta1_y_zzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 545);

    auto ta1_y_zzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 546);

    auto ta1_y_zzz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 547);

    auto ta1_y_zzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 548);

    auto ta1_y_zzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 549);

    auto ta1_y_zzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 550);

    auto ta1_y_zzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 551);

    auto ta1_y_zzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 552);

    auto ta1_y_zzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 553);

    auto ta1_y_zzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 554);

    auto ta1_y_zzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 555);

    auto ta1_y_zzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 556);

    auto ta1_y_zzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 557);

    auto ta1_y_zzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 558);

    auto ta1_y_zzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 559);

    auto ta1_z_xxx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 560);

    auto ta1_z_xxx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 561);

    auto ta1_z_xxx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 562);

    auto ta1_z_xxx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 563);

    auto ta1_z_xxx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 564);

    auto ta1_z_xxx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 565);

    auto ta1_z_xxx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 566);

    auto ta1_z_xxx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 567);

    auto ta1_z_xxx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 568);

    auto ta1_z_xxx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 569);

    auto ta1_z_xxx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 570);

    auto ta1_z_xxx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 571);

    auto ta1_z_xxx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 572);

    auto ta1_z_xxx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 573);

    auto ta1_z_xxx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 574);

    auto ta1_z_xxx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 575);

    auto ta1_z_xxx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 576);

    auto ta1_z_xxx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 577);

    auto ta1_z_xxx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 578);

    auto ta1_z_xxx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 579);

    auto ta1_z_xxx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 580);

    auto ta1_z_xxx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 581);

    auto ta1_z_xxx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 582);

    auto ta1_z_xxx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 583);

    auto ta1_z_xxx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 584);

    auto ta1_z_xxx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 585);

    auto ta1_z_xxx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 586);

    auto ta1_z_xxx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 587);

    auto ta1_z_xxy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 588);

    auto ta1_z_xxy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 589);

    auto ta1_z_xxy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 590);

    auto ta1_z_xxy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 591);

    auto ta1_z_xxy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 593);

    auto ta1_z_xxy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 594);

    auto ta1_z_xxy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 597);

    auto ta1_z_xxy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 598);

    auto ta1_z_xxy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 602);

    auto ta1_z_xxy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 603);

    auto ta1_z_xxy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 608);

    auto ta1_z_xxy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 609);

    auto ta1_z_xxy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 610);

    auto ta1_z_xxy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 611);

    auto ta1_z_xxy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 612);

    auto ta1_z_xxy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 613);

    auto ta1_z_xxy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 614);

    auto ta1_z_xxz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 616);

    auto ta1_z_xxz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 617);

    auto ta1_z_xxz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 618);

    auto ta1_z_xxz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 619);

    auto ta1_z_xxz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 620);

    auto ta1_z_xxz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 621);

    auto ta1_z_xxz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 622);

    auto ta1_z_xxz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 623);

    auto ta1_z_xxz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 624);

    auto ta1_z_xxz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 625);

    auto ta1_z_xxz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 626);

    auto ta1_z_xxz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 627);

    auto ta1_z_xxz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 628);

    auto ta1_z_xxz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 629);

    auto ta1_z_xxz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 630);

    auto ta1_z_xxz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 631);

    auto ta1_z_xxz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 632);

    auto ta1_z_xxz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 633);

    auto ta1_z_xxz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 634);

    auto ta1_z_xxz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 635);

    auto ta1_z_xxz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 636);

    auto ta1_z_xxz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 638);

    auto ta1_z_xxz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 639);

    auto ta1_z_xxz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 640);

    auto ta1_z_xxz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 641);

    auto ta1_z_xxz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 642);

    auto ta1_z_xxz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 643);

    auto ta1_z_xyy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 644);

    auto ta1_z_xyy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 645);

    auto ta1_z_xyy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 647);

    auto ta1_z_xyy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 648);

    auto ta1_z_xyy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 650);

    auto ta1_z_xyy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 651);

    auto ta1_z_xyy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 652);

    auto ta1_z_xyy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 654);

    auto ta1_z_xyy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 655);

    auto ta1_z_xyy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 656);

    auto ta1_z_xyy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 657);

    auto ta1_z_xyy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 659);

    auto ta1_z_xyy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 660);

    auto ta1_z_xyy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 661);

    auto ta1_z_xyy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 662);

    auto ta1_z_xyy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 663);

    auto ta1_z_xyy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 665);

    auto ta1_z_xyy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 666);

    auto ta1_z_xyy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 667);

    auto ta1_z_xyy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 668);

    auto ta1_z_xyy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 669);

    auto ta1_z_xyy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 670);

    auto ta1_z_xyy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 671);

    auto ta1_z_xyz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 694);

    auto ta1_z_xyz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 695);

    auto ta1_z_xyz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 696);

    auto ta1_z_xyz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 697);

    auto ta1_z_xyz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 698);

    auto ta1_z_xzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 700);

    auto ta1_z_xzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 702);

    auto ta1_z_xzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 704);

    auto ta1_z_xzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 705);

    auto ta1_z_xzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 707);

    auto ta1_z_xzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 708);

    auto ta1_z_xzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 709);

    auto ta1_z_xzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 711);

    auto ta1_z_xzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 712);

    auto ta1_z_xzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 713);

    auto ta1_z_xzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 714);

    auto ta1_z_xzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 716);

    auto ta1_z_xzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 717);

    auto ta1_z_xzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 718);

    auto ta1_z_xzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 719);

    auto ta1_z_xzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 720);

    auto ta1_z_xzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 721);

    auto ta1_z_xzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 722);

    auto ta1_z_xzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 723);

    auto ta1_z_xzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 724);

    auto ta1_z_xzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 725);

    auto ta1_z_xzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 726);

    auto ta1_z_xzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 727);

    auto ta1_z_yyy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 728);

    auto ta1_z_yyy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 729);

    auto ta1_z_yyy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 730);

    auto ta1_z_yyy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 731);

    auto ta1_z_yyy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 732);

    auto ta1_z_yyy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 733);

    auto ta1_z_yyy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 734);

    auto ta1_z_yyy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 735);

    auto ta1_z_yyy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 736);

    auto ta1_z_yyy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 737);

    auto ta1_z_yyy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 738);

    auto ta1_z_yyy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 739);

    auto ta1_z_yyy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 740);

    auto ta1_z_yyy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 741);

    auto ta1_z_yyy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 742);

    auto ta1_z_yyy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 743);

    auto ta1_z_yyy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 744);

    auto ta1_z_yyy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 745);

    auto ta1_z_yyy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 746);

    auto ta1_z_yyy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 747);

    auto ta1_z_yyy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 748);

    auto ta1_z_yyy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 749);

    auto ta1_z_yyy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 750);

    auto ta1_z_yyy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 751);

    auto ta1_z_yyy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 752);

    auto ta1_z_yyy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 753);

    auto ta1_z_yyy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 754);

    auto ta1_z_yyy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 755);

    auto ta1_z_yyz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 757);

    auto ta1_z_yyz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 758);

    auto ta1_z_yyz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 759);

    auto ta1_z_yyz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 760);

    auto ta1_z_yyz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 761);

    auto ta1_z_yyz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 762);

    auto ta1_z_yyz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 763);

    auto ta1_z_yyz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 764);

    auto ta1_z_yyz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 765);

    auto ta1_z_yyz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 766);

    auto ta1_z_yyz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 767);

    auto ta1_z_yyz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 768);

    auto ta1_z_yyz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 769);

    auto ta1_z_yyz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 770);

    auto ta1_z_yyz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 771);

    auto ta1_z_yyz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 772);

    auto ta1_z_yyz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 773);

    auto ta1_z_yyz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 774);

    auto ta1_z_yyz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 775);

    auto ta1_z_yyz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 776);

    auto ta1_z_yyz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 777);

    auto ta1_z_yyz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 778);

    auto ta1_z_yyz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 779);

    auto ta1_z_yyz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 780);

    auto ta1_z_yyz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 781);

    auto ta1_z_yyz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 782);

    auto ta1_z_yyz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 783);

    auto ta1_z_yzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 784);

    auto ta1_z_yzz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 785);

    auto ta1_z_yzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 786);

    auto ta1_z_yzz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 787);

    auto ta1_z_yzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 788);

    auto ta1_z_yzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 789);

    auto ta1_z_yzz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 790);

    auto ta1_z_yzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 791);

    auto ta1_z_yzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 792);

    auto ta1_z_yzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 793);

    auto ta1_z_yzz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 794);

    auto ta1_z_yzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 795);

    auto ta1_z_yzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 796);

    auto ta1_z_yzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 797);

    auto ta1_z_yzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 798);

    auto ta1_z_yzz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 799);

    auto ta1_z_yzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 800);

    auto ta1_z_yzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 801);

    auto ta1_z_yzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 802);

    auto ta1_z_yzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 803);

    auto ta1_z_yzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 804);

    auto ta1_z_yzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 805);

    auto ta1_z_yzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 806);

    auto ta1_z_yzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 807);

    auto ta1_z_yzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 808);

    auto ta1_z_yzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 809);

    auto ta1_z_yzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 810);

    auto ta1_z_yzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 811);

    auto ta1_z_zzz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fi + 812);

    auto ta1_z_zzz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 813);

    auto ta1_z_zzz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 814);

    auto ta1_z_zzz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 815);

    auto ta1_z_zzz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 816);

    auto ta1_z_zzz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 817);

    auto ta1_z_zzz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 818);

    auto ta1_z_zzz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 819);

    auto ta1_z_zzz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 820);

    auto ta1_z_zzz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 821);

    auto ta1_z_zzz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 822);

    auto ta1_z_zzz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 823);

    auto ta1_z_zzz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 824);

    auto ta1_z_zzz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 825);

    auto ta1_z_zzz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 826);

    auto ta1_z_zzz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 827);

    auto ta1_z_zzz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 828);

    auto ta1_z_zzz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 829);

    auto ta1_z_zzz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 830);

    auto ta1_z_zzz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 831);

    auto ta1_z_zzz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 832);

    auto ta1_z_zzz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fi + 833);

    auto ta1_z_zzz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 834);

    auto ta1_z_zzz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 835);

    auto ta1_z_zzz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 836);

    auto ta1_z_zzz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 837);

    auto ta1_z_zzz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 838);

    auto ta1_z_zzz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fi + 839);

    // Set up 0-28 components of targeted buffer : GI

    auto ta1_x_xxxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi);

    auto ta1_x_xxxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1);

    auto ta1_x_xxxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 2);

    auto ta1_x_xxxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 3);

    auto ta1_x_xxxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 4);

    auto ta1_x_xxxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 5);

    auto ta1_x_xxxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 6);

    auto ta1_x_xxxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 7);

    auto ta1_x_xxxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 8);

    auto ta1_x_xxxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 9);

    auto ta1_x_xxxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 10);

    auto ta1_x_xxxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 11);

    auto ta1_x_xxxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 12);

    auto ta1_x_xxxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 13);

    auto ta1_x_xxxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 14);

    auto ta1_x_xxxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 15);

    auto ta1_x_xxxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 16);

    auto ta1_x_xxxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 17);

    auto ta1_x_xxxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 18);

    auto ta1_x_xxxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 19);

    auto ta1_x_xxxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 20);

    auto ta1_x_xxxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 21);

    auto ta1_x_xxxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 22);

    auto ta1_x_xxxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 23);

    auto ta1_x_xxxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 24);

    auto ta1_x_xxxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 25);

    auto ta1_x_xxxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 26);

    auto ta1_x_xxxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 27);

#pragma omp simd aligned(pa_x,                    \
                             pc_x,                \
                             ta1_x_xx_xxxxxx_0,   \
                             ta1_x_xx_xxxxxx_1,   \
                             ta1_x_xx_xxxxxy_0,   \
                             ta1_x_xx_xxxxxy_1,   \
                             ta1_x_xx_xxxxxz_0,   \
                             ta1_x_xx_xxxxxz_1,   \
                             ta1_x_xx_xxxxyy_0,   \
                             ta1_x_xx_xxxxyy_1,   \
                             ta1_x_xx_xxxxyz_0,   \
                             ta1_x_xx_xxxxyz_1,   \
                             ta1_x_xx_xxxxzz_0,   \
                             ta1_x_xx_xxxxzz_1,   \
                             ta1_x_xx_xxxyyy_0,   \
                             ta1_x_xx_xxxyyy_1,   \
                             ta1_x_xx_xxxyyz_0,   \
                             ta1_x_xx_xxxyyz_1,   \
                             ta1_x_xx_xxxyzz_0,   \
                             ta1_x_xx_xxxyzz_1,   \
                             ta1_x_xx_xxxzzz_0,   \
                             ta1_x_xx_xxxzzz_1,   \
                             ta1_x_xx_xxyyyy_0,   \
                             ta1_x_xx_xxyyyy_1,   \
                             ta1_x_xx_xxyyyz_0,   \
                             ta1_x_xx_xxyyyz_1,   \
                             ta1_x_xx_xxyyzz_0,   \
                             ta1_x_xx_xxyyzz_1,   \
                             ta1_x_xx_xxyzzz_0,   \
                             ta1_x_xx_xxyzzz_1,   \
                             ta1_x_xx_xxzzzz_0,   \
                             ta1_x_xx_xxzzzz_1,   \
                             ta1_x_xx_xyyyyy_0,   \
                             ta1_x_xx_xyyyyy_1,   \
                             ta1_x_xx_xyyyyz_0,   \
                             ta1_x_xx_xyyyyz_1,   \
                             ta1_x_xx_xyyyzz_0,   \
                             ta1_x_xx_xyyyzz_1,   \
                             ta1_x_xx_xyyzzz_0,   \
                             ta1_x_xx_xyyzzz_1,   \
                             ta1_x_xx_xyzzzz_0,   \
                             ta1_x_xx_xyzzzz_1,   \
                             ta1_x_xx_xzzzzz_0,   \
                             ta1_x_xx_xzzzzz_1,   \
                             ta1_x_xx_yyyyyy_0,   \
                             ta1_x_xx_yyyyyy_1,   \
                             ta1_x_xx_yyyyyz_0,   \
                             ta1_x_xx_yyyyyz_1,   \
                             ta1_x_xx_yyyyzz_0,   \
                             ta1_x_xx_yyyyzz_1,   \
                             ta1_x_xx_yyyzzz_0,   \
                             ta1_x_xx_yyyzzz_1,   \
                             ta1_x_xx_yyzzzz_0,   \
                             ta1_x_xx_yyzzzz_1,   \
                             ta1_x_xx_yzzzzz_0,   \
                             ta1_x_xx_yzzzzz_1,   \
                             ta1_x_xx_zzzzzz_0,   \
                             ta1_x_xx_zzzzzz_1,   \
                             ta1_x_xxx_xxxxx_0,   \
                             ta1_x_xxx_xxxxx_1,   \
                             ta1_x_xxx_xxxxxx_0,  \
                             ta1_x_xxx_xxxxxx_1,  \
                             ta1_x_xxx_xxxxxy_0,  \
                             ta1_x_xxx_xxxxxy_1,  \
                             ta1_x_xxx_xxxxxz_0,  \
                             ta1_x_xxx_xxxxxz_1,  \
                             ta1_x_xxx_xxxxy_0,   \
                             ta1_x_xxx_xxxxy_1,   \
                             ta1_x_xxx_xxxxyy_0,  \
                             ta1_x_xxx_xxxxyy_1,  \
                             ta1_x_xxx_xxxxyz_0,  \
                             ta1_x_xxx_xxxxyz_1,  \
                             ta1_x_xxx_xxxxz_0,   \
                             ta1_x_xxx_xxxxz_1,   \
                             ta1_x_xxx_xxxxzz_0,  \
                             ta1_x_xxx_xxxxzz_1,  \
                             ta1_x_xxx_xxxyy_0,   \
                             ta1_x_xxx_xxxyy_1,   \
                             ta1_x_xxx_xxxyyy_0,  \
                             ta1_x_xxx_xxxyyy_1,  \
                             ta1_x_xxx_xxxyyz_0,  \
                             ta1_x_xxx_xxxyyz_1,  \
                             ta1_x_xxx_xxxyz_0,   \
                             ta1_x_xxx_xxxyz_1,   \
                             ta1_x_xxx_xxxyzz_0,  \
                             ta1_x_xxx_xxxyzz_1,  \
                             ta1_x_xxx_xxxzz_0,   \
                             ta1_x_xxx_xxxzz_1,   \
                             ta1_x_xxx_xxxzzz_0,  \
                             ta1_x_xxx_xxxzzz_1,  \
                             ta1_x_xxx_xxyyy_0,   \
                             ta1_x_xxx_xxyyy_1,   \
                             ta1_x_xxx_xxyyyy_0,  \
                             ta1_x_xxx_xxyyyy_1,  \
                             ta1_x_xxx_xxyyyz_0,  \
                             ta1_x_xxx_xxyyyz_1,  \
                             ta1_x_xxx_xxyyz_0,   \
                             ta1_x_xxx_xxyyz_1,   \
                             ta1_x_xxx_xxyyzz_0,  \
                             ta1_x_xxx_xxyyzz_1,  \
                             ta1_x_xxx_xxyzz_0,   \
                             ta1_x_xxx_xxyzz_1,   \
                             ta1_x_xxx_xxyzzz_0,  \
                             ta1_x_xxx_xxyzzz_1,  \
                             ta1_x_xxx_xxzzz_0,   \
                             ta1_x_xxx_xxzzz_1,   \
                             ta1_x_xxx_xxzzzz_0,  \
                             ta1_x_xxx_xxzzzz_1,  \
                             ta1_x_xxx_xyyyy_0,   \
                             ta1_x_xxx_xyyyy_1,   \
                             ta1_x_xxx_xyyyyy_0,  \
                             ta1_x_xxx_xyyyyy_1,  \
                             ta1_x_xxx_xyyyyz_0,  \
                             ta1_x_xxx_xyyyyz_1,  \
                             ta1_x_xxx_xyyyz_0,   \
                             ta1_x_xxx_xyyyz_1,   \
                             ta1_x_xxx_xyyyzz_0,  \
                             ta1_x_xxx_xyyyzz_1,  \
                             ta1_x_xxx_xyyzz_0,   \
                             ta1_x_xxx_xyyzz_1,   \
                             ta1_x_xxx_xyyzzz_0,  \
                             ta1_x_xxx_xyyzzz_1,  \
                             ta1_x_xxx_xyzzz_0,   \
                             ta1_x_xxx_xyzzz_1,   \
                             ta1_x_xxx_xyzzzz_0,  \
                             ta1_x_xxx_xyzzzz_1,  \
                             ta1_x_xxx_xzzzz_0,   \
                             ta1_x_xxx_xzzzz_1,   \
                             ta1_x_xxx_xzzzzz_0,  \
                             ta1_x_xxx_xzzzzz_1,  \
                             ta1_x_xxx_yyyyy_0,   \
                             ta1_x_xxx_yyyyy_1,   \
                             ta1_x_xxx_yyyyyy_0,  \
                             ta1_x_xxx_yyyyyy_1,  \
                             ta1_x_xxx_yyyyyz_0,  \
                             ta1_x_xxx_yyyyyz_1,  \
                             ta1_x_xxx_yyyyz_0,   \
                             ta1_x_xxx_yyyyz_1,   \
                             ta1_x_xxx_yyyyzz_0,  \
                             ta1_x_xxx_yyyyzz_1,  \
                             ta1_x_xxx_yyyzz_0,   \
                             ta1_x_xxx_yyyzz_1,   \
                             ta1_x_xxx_yyyzzz_0,  \
                             ta1_x_xxx_yyyzzz_1,  \
                             ta1_x_xxx_yyzzz_0,   \
                             ta1_x_xxx_yyzzz_1,   \
                             ta1_x_xxx_yyzzzz_0,  \
                             ta1_x_xxx_yyzzzz_1,  \
                             ta1_x_xxx_yzzzz_0,   \
                             ta1_x_xxx_yzzzz_1,   \
                             ta1_x_xxx_yzzzzz_0,  \
                             ta1_x_xxx_yzzzzz_1,  \
                             ta1_x_xxx_zzzzz_0,   \
                             ta1_x_xxx_zzzzz_1,   \
                             ta1_x_xxx_zzzzzz_0,  \
                             ta1_x_xxx_zzzzzz_1,  \
                             ta1_x_xxxx_xxxxxx_0, \
                             ta1_x_xxxx_xxxxxy_0, \
                             ta1_x_xxxx_xxxxxz_0, \
                             ta1_x_xxxx_xxxxyy_0, \
                             ta1_x_xxxx_xxxxyz_0, \
                             ta1_x_xxxx_xxxxzz_0, \
                             ta1_x_xxxx_xxxyyy_0, \
                             ta1_x_xxxx_xxxyyz_0, \
                             ta1_x_xxxx_xxxyzz_0, \
                             ta1_x_xxxx_xxxzzz_0, \
                             ta1_x_xxxx_xxyyyy_0, \
                             ta1_x_xxxx_xxyyyz_0, \
                             ta1_x_xxxx_xxyyzz_0, \
                             ta1_x_xxxx_xxyzzz_0, \
                             ta1_x_xxxx_xxzzzz_0, \
                             ta1_x_xxxx_xyyyyy_0, \
                             ta1_x_xxxx_xyyyyz_0, \
                             ta1_x_xxxx_xyyyzz_0, \
                             ta1_x_xxxx_xyyzzz_0, \
                             ta1_x_xxxx_xyzzzz_0, \
                             ta1_x_xxxx_xzzzzz_0, \
                             ta1_x_xxxx_yyyyyy_0, \
                             ta1_x_xxxx_yyyyyz_0, \
                             ta1_x_xxxx_yyyyzz_0, \
                             ta1_x_xxxx_yyyzzz_0, \
                             ta1_x_xxxx_yyzzzz_0, \
                             ta1_x_xxxx_yzzzzz_0, \
                             ta1_x_xxxx_zzzzzz_0, \
                             ta_xxx_xxxxxx_1,     \
                             ta_xxx_xxxxxy_1,     \
                             ta_xxx_xxxxxz_1,     \
                             ta_xxx_xxxxyy_1,     \
                             ta_xxx_xxxxyz_1,     \
                             ta_xxx_xxxxzz_1,     \
                             ta_xxx_xxxyyy_1,     \
                             ta_xxx_xxxyyz_1,     \
                             ta_xxx_xxxyzz_1,     \
                             ta_xxx_xxxzzz_1,     \
                             ta_xxx_xxyyyy_1,     \
                             ta_xxx_xxyyyz_1,     \
                             ta_xxx_xxyyzz_1,     \
                             ta_xxx_xxyzzz_1,     \
                             ta_xxx_xxzzzz_1,     \
                             ta_xxx_xyyyyy_1,     \
                             ta_xxx_xyyyyz_1,     \
                             ta_xxx_xyyyzz_1,     \
                             ta_xxx_xyyzzz_1,     \
                             ta_xxx_xyzzzz_1,     \
                             ta_xxx_xzzzzz_1,     \
                             ta_xxx_yyyyyy_1,     \
                             ta_xxx_yyyyyz_1,     \
                             ta_xxx_yyyyzz_1,     \
                             ta_xxx_yyyzzz_1,     \
                             ta_xxx_yyzzzz_1,     \
                             ta_xxx_yzzzzz_1,     \
                             ta_xxx_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxx_xxxxxx_0[i] = 3.0 * ta1_x_xx_xxxxxx_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxxx_1[i] * fe_0 + 6.0 * ta1_x_xxx_xxxxx_0[i] * fe_0 -
                                 6.0 * ta1_x_xxx_xxxxx_1[i] * fe_0 + ta_xxx_xxxxxx_1[i] + ta1_x_xxx_xxxxxx_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxxxx_1[i] * pc_x[i];

        ta1_x_xxxx_xxxxxy_0[i] = 3.0 * ta1_x_xx_xxxxxy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxxy_1[i] * fe_0 + 5.0 * ta1_x_xxx_xxxxy_0[i] * fe_0 -
                                 5.0 * ta1_x_xxx_xxxxy_1[i] * fe_0 + ta_xxx_xxxxxy_1[i] + ta1_x_xxx_xxxxxy_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxxxy_1[i] * pc_x[i];

        ta1_x_xxxx_xxxxxz_0[i] = 3.0 * ta1_x_xx_xxxxxz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxxz_1[i] * fe_0 + 5.0 * ta1_x_xxx_xxxxz_0[i] * fe_0 -
                                 5.0 * ta1_x_xxx_xxxxz_1[i] * fe_0 + ta_xxx_xxxxxz_1[i] + ta1_x_xxx_xxxxxz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxxxz_1[i] * pc_x[i];

        ta1_x_xxxx_xxxxyy_0[i] = 3.0 * ta1_x_xx_xxxxyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxyy_1[i] * fe_0 + 4.0 * ta1_x_xxx_xxxyy_0[i] * fe_0 -
                                 4.0 * ta1_x_xxx_xxxyy_1[i] * fe_0 + ta_xxx_xxxxyy_1[i] + ta1_x_xxx_xxxxyy_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxxyy_1[i] * pc_x[i];

        ta1_x_xxxx_xxxxyz_0[i] = 3.0 * ta1_x_xx_xxxxyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxyz_1[i] * fe_0 + 4.0 * ta1_x_xxx_xxxyz_0[i] * fe_0 -
                                 4.0 * ta1_x_xxx_xxxyz_1[i] * fe_0 + ta_xxx_xxxxyz_1[i] + ta1_x_xxx_xxxxyz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxxyz_1[i] * pc_x[i];

        ta1_x_xxxx_xxxxzz_0[i] = 3.0 * ta1_x_xx_xxxxzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxzz_1[i] * fe_0 + 4.0 * ta1_x_xxx_xxxzz_0[i] * fe_0 -
                                 4.0 * ta1_x_xxx_xxxzz_1[i] * fe_0 + ta_xxx_xxxxzz_1[i] + ta1_x_xxx_xxxxzz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxxzz_1[i] * pc_x[i];

        ta1_x_xxxx_xxxyyy_0[i] = 3.0 * ta1_x_xx_xxxyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxyyy_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxyyy_0[i] * fe_0 -
                                 3.0 * ta1_x_xxx_xxyyy_1[i] * fe_0 + ta_xxx_xxxyyy_1[i] + ta1_x_xxx_xxxyyy_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxyyy_1[i] * pc_x[i];

        ta1_x_xxxx_xxxyyz_0[i] = 3.0 * ta1_x_xx_xxxyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxyyz_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxx_xxyyz_1[i] * fe_0 + ta_xxx_xxxyyz_1[i] + ta1_x_xxx_xxxyyz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxyyz_1[i] * pc_x[i];

        ta1_x_xxxx_xxxyzz_0[i] = 3.0 * ta1_x_xx_xxxyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxyzz_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxx_xxyzz_1[i] * fe_0 + ta_xxx_xxxyzz_1[i] + ta1_x_xxx_xxxyzz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxyzz_1[i] * pc_x[i];

        ta1_x_xxxx_xxxzzz_0[i] = 3.0 * ta1_x_xx_xxxzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxzzz_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxzzz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxx_xxzzz_1[i] * fe_0 + ta_xxx_xxxzzz_1[i] + ta1_x_xxx_xxxzzz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxxzzz_1[i] * pc_x[i];

        ta1_x_xxxx_xxyyyy_0[i] = 3.0 * ta1_x_xx_xxyyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyyyy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyyyy_0[i] * fe_0 -
                                 2.0 * ta1_x_xxx_xyyyy_1[i] * fe_0 + ta_xxx_xxyyyy_1[i] + ta1_x_xxx_xxyyyy_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxyyyy_1[i] * pc_x[i];

        ta1_x_xxxx_xxyyyz_0[i] = 3.0 * ta1_x_xx_xxyyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyyyz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxx_xyyyz_1[i] * fe_0 + ta_xxx_xxyyyz_1[i] + ta1_x_xxx_xxyyyz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxyyyz_1[i] * pc_x[i];

        ta1_x_xxxx_xxyyzz_0[i] = 3.0 * ta1_x_xx_xxyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyyzz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxx_xyyzz_1[i] * fe_0 + ta_xxx_xxyyzz_1[i] + ta1_x_xxx_xxyyzz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxyyzz_1[i] * pc_x[i];

        ta1_x_xxxx_xxyzzz_0[i] = 3.0 * ta1_x_xx_xxyzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyzzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxx_xyzzz_1[i] * fe_0 + ta_xxx_xxyzzz_1[i] + ta1_x_xxx_xxyzzz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxyzzz_1[i] * pc_x[i];

        ta1_x_xxxx_xxzzzz_0[i] = 3.0 * ta1_x_xx_xxzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxzzzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xzzzz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxx_xzzzz_1[i] * fe_0 + ta_xxx_xxzzzz_1[i] + ta1_x_xxx_xxzzzz_0[i] * pa_x[i] -
                                 ta1_x_xxx_xxzzzz_1[i] * pc_x[i];

        ta1_x_xxxx_xyyyyy_0[i] = 3.0 * ta1_x_xx_xyyyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyyyy_1[i] * fe_0 + ta1_x_xxx_yyyyy_0[i] * fe_0 -
                                 ta1_x_xxx_yyyyy_1[i] * fe_0 + ta_xxx_xyyyyy_1[i] + ta1_x_xxx_xyyyyy_0[i] * pa_x[i] - ta1_x_xxx_xyyyyy_1[i] * pc_x[i];

        ta1_x_xxxx_xyyyyz_0[i] = 3.0 * ta1_x_xx_xyyyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyyyz_1[i] * fe_0 + ta1_x_xxx_yyyyz_0[i] * fe_0 -
                                 ta1_x_xxx_yyyyz_1[i] * fe_0 + ta_xxx_xyyyyz_1[i] + ta1_x_xxx_xyyyyz_0[i] * pa_x[i] - ta1_x_xxx_xyyyyz_1[i] * pc_x[i];

        ta1_x_xxxx_xyyyzz_0[i] = 3.0 * ta1_x_xx_xyyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyyzz_1[i] * fe_0 + ta1_x_xxx_yyyzz_0[i] * fe_0 -
                                 ta1_x_xxx_yyyzz_1[i] * fe_0 + ta_xxx_xyyyzz_1[i] + ta1_x_xxx_xyyyzz_0[i] * pa_x[i] - ta1_x_xxx_xyyyzz_1[i] * pc_x[i];

        ta1_x_xxxx_xyyzzz_0[i] = 3.0 * ta1_x_xx_xyyzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyzzz_1[i] * fe_0 + ta1_x_xxx_yyzzz_0[i] * fe_0 -
                                 ta1_x_xxx_yyzzz_1[i] * fe_0 + ta_xxx_xyyzzz_1[i] + ta1_x_xxx_xyyzzz_0[i] * pa_x[i] - ta1_x_xxx_xyyzzz_1[i] * pc_x[i];

        ta1_x_xxxx_xyzzzz_0[i] = 3.0 * ta1_x_xx_xyzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyzzzz_1[i] * fe_0 + ta1_x_xxx_yzzzz_0[i] * fe_0 -
                                 ta1_x_xxx_yzzzz_1[i] * fe_0 + ta_xxx_xyzzzz_1[i] + ta1_x_xxx_xyzzzz_0[i] * pa_x[i] - ta1_x_xxx_xyzzzz_1[i] * pc_x[i];

        ta1_x_xxxx_xzzzzz_0[i] = 3.0 * ta1_x_xx_xzzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xzzzzz_1[i] * fe_0 + ta1_x_xxx_zzzzz_0[i] * fe_0 -
                                 ta1_x_xxx_zzzzz_1[i] * fe_0 + ta_xxx_xzzzzz_1[i] + ta1_x_xxx_xzzzzz_0[i] * pa_x[i] - ta1_x_xxx_xzzzzz_1[i] * pc_x[i];

        ta1_x_xxxx_yyyyyy_0[i] = 3.0 * ta1_x_xx_yyyyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyyyy_1[i] * fe_0 + ta_xxx_yyyyyy_1[i] +
                                 ta1_x_xxx_yyyyyy_0[i] * pa_x[i] - ta1_x_xxx_yyyyyy_1[i] * pc_x[i];

        ta1_x_xxxx_yyyyyz_0[i] = 3.0 * ta1_x_xx_yyyyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyyyz_1[i] * fe_0 + ta_xxx_yyyyyz_1[i] +
                                 ta1_x_xxx_yyyyyz_0[i] * pa_x[i] - ta1_x_xxx_yyyyyz_1[i] * pc_x[i];

        ta1_x_xxxx_yyyyzz_0[i] = 3.0 * ta1_x_xx_yyyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyyzz_1[i] * fe_0 + ta_xxx_yyyyzz_1[i] +
                                 ta1_x_xxx_yyyyzz_0[i] * pa_x[i] - ta1_x_xxx_yyyyzz_1[i] * pc_x[i];

        ta1_x_xxxx_yyyzzz_0[i] = 3.0 * ta1_x_xx_yyyzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyzzz_1[i] * fe_0 + ta_xxx_yyyzzz_1[i] +
                                 ta1_x_xxx_yyyzzz_0[i] * pa_x[i] - ta1_x_xxx_yyyzzz_1[i] * pc_x[i];

        ta1_x_xxxx_yyzzzz_0[i] = 3.0 * ta1_x_xx_yyzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyzzzz_1[i] * fe_0 + ta_xxx_yyzzzz_1[i] +
                                 ta1_x_xxx_yyzzzz_0[i] * pa_x[i] - ta1_x_xxx_yyzzzz_1[i] * pc_x[i];

        ta1_x_xxxx_yzzzzz_0[i] = 3.0 * ta1_x_xx_yzzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yzzzzz_1[i] * fe_0 + ta_xxx_yzzzzz_1[i] +
                                 ta1_x_xxx_yzzzzz_0[i] * pa_x[i] - ta1_x_xxx_yzzzzz_1[i] * pc_x[i];

        ta1_x_xxxx_zzzzzz_0[i] = 3.0 * ta1_x_xx_zzzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_zzzzzz_1[i] * fe_0 + ta_xxx_zzzzzz_1[i] +
                                 ta1_x_xxx_zzzzzz_0[i] * pa_x[i] - ta1_x_xxx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : GI

    auto ta1_x_xxxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 28);

    auto ta1_x_xxxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 29);

    auto ta1_x_xxxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 30);

    auto ta1_x_xxxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 31);

    auto ta1_x_xxxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 32);

    auto ta1_x_xxxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 33);

    auto ta1_x_xxxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 34);

    auto ta1_x_xxxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 35);

    auto ta1_x_xxxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 36);

    auto ta1_x_xxxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 37);

    auto ta1_x_xxxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 38);

    auto ta1_x_xxxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 39);

    auto ta1_x_xxxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 40);

    auto ta1_x_xxxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 41);

    auto ta1_x_xxxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 42);

    auto ta1_x_xxxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 43);

    auto ta1_x_xxxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 44);

    auto ta1_x_xxxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 45);

    auto ta1_x_xxxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 46);

    auto ta1_x_xxxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 47);

    auto ta1_x_xxxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 48);

    auto ta1_x_xxxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 49);

    auto ta1_x_xxxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 50);

    auto ta1_x_xxxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 51);

    auto ta1_x_xxxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 52);

    auto ta1_x_xxxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 53);

    auto ta1_x_xxxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 54);

    auto ta1_x_xxxy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 55);

#pragma omp simd aligned(pa_y,                    \
                             pc_y,                \
                             ta1_x_xxx_xxxxx_0,   \
                             ta1_x_xxx_xxxxx_1,   \
                             ta1_x_xxx_xxxxxx_0,  \
                             ta1_x_xxx_xxxxxx_1,  \
                             ta1_x_xxx_xxxxxy_0,  \
                             ta1_x_xxx_xxxxxy_1,  \
                             ta1_x_xxx_xxxxxz_0,  \
                             ta1_x_xxx_xxxxxz_1,  \
                             ta1_x_xxx_xxxxy_0,   \
                             ta1_x_xxx_xxxxy_1,   \
                             ta1_x_xxx_xxxxyy_0,  \
                             ta1_x_xxx_xxxxyy_1,  \
                             ta1_x_xxx_xxxxyz_0,  \
                             ta1_x_xxx_xxxxyz_1,  \
                             ta1_x_xxx_xxxxz_0,   \
                             ta1_x_xxx_xxxxz_1,   \
                             ta1_x_xxx_xxxxzz_0,  \
                             ta1_x_xxx_xxxxzz_1,  \
                             ta1_x_xxx_xxxyy_0,   \
                             ta1_x_xxx_xxxyy_1,   \
                             ta1_x_xxx_xxxyyy_0,  \
                             ta1_x_xxx_xxxyyy_1,  \
                             ta1_x_xxx_xxxyyz_0,  \
                             ta1_x_xxx_xxxyyz_1,  \
                             ta1_x_xxx_xxxyz_0,   \
                             ta1_x_xxx_xxxyz_1,   \
                             ta1_x_xxx_xxxyzz_0,  \
                             ta1_x_xxx_xxxyzz_1,  \
                             ta1_x_xxx_xxxzz_0,   \
                             ta1_x_xxx_xxxzz_1,   \
                             ta1_x_xxx_xxxzzz_0,  \
                             ta1_x_xxx_xxxzzz_1,  \
                             ta1_x_xxx_xxyyy_0,   \
                             ta1_x_xxx_xxyyy_1,   \
                             ta1_x_xxx_xxyyyy_0,  \
                             ta1_x_xxx_xxyyyy_1,  \
                             ta1_x_xxx_xxyyyz_0,  \
                             ta1_x_xxx_xxyyyz_1,  \
                             ta1_x_xxx_xxyyz_0,   \
                             ta1_x_xxx_xxyyz_1,   \
                             ta1_x_xxx_xxyyzz_0,  \
                             ta1_x_xxx_xxyyzz_1,  \
                             ta1_x_xxx_xxyzz_0,   \
                             ta1_x_xxx_xxyzz_1,   \
                             ta1_x_xxx_xxyzzz_0,  \
                             ta1_x_xxx_xxyzzz_1,  \
                             ta1_x_xxx_xxzzz_0,   \
                             ta1_x_xxx_xxzzz_1,   \
                             ta1_x_xxx_xxzzzz_0,  \
                             ta1_x_xxx_xxzzzz_1,  \
                             ta1_x_xxx_xyyyy_0,   \
                             ta1_x_xxx_xyyyy_1,   \
                             ta1_x_xxx_xyyyyy_0,  \
                             ta1_x_xxx_xyyyyy_1,  \
                             ta1_x_xxx_xyyyyz_0,  \
                             ta1_x_xxx_xyyyyz_1,  \
                             ta1_x_xxx_xyyyz_0,   \
                             ta1_x_xxx_xyyyz_1,   \
                             ta1_x_xxx_xyyyzz_0,  \
                             ta1_x_xxx_xyyyzz_1,  \
                             ta1_x_xxx_xyyzz_0,   \
                             ta1_x_xxx_xyyzz_1,   \
                             ta1_x_xxx_xyyzzz_0,  \
                             ta1_x_xxx_xyyzzz_1,  \
                             ta1_x_xxx_xyzzz_0,   \
                             ta1_x_xxx_xyzzz_1,   \
                             ta1_x_xxx_xyzzzz_0,  \
                             ta1_x_xxx_xyzzzz_1,  \
                             ta1_x_xxx_xzzzz_0,   \
                             ta1_x_xxx_xzzzz_1,   \
                             ta1_x_xxx_xzzzzz_0,  \
                             ta1_x_xxx_xzzzzz_1,  \
                             ta1_x_xxx_yyyyy_0,   \
                             ta1_x_xxx_yyyyy_1,   \
                             ta1_x_xxx_yyyyyy_0,  \
                             ta1_x_xxx_yyyyyy_1,  \
                             ta1_x_xxx_yyyyyz_0,  \
                             ta1_x_xxx_yyyyyz_1,  \
                             ta1_x_xxx_yyyyz_0,   \
                             ta1_x_xxx_yyyyz_1,   \
                             ta1_x_xxx_yyyyzz_0,  \
                             ta1_x_xxx_yyyyzz_1,  \
                             ta1_x_xxx_yyyzz_0,   \
                             ta1_x_xxx_yyyzz_1,   \
                             ta1_x_xxx_yyyzzz_0,  \
                             ta1_x_xxx_yyyzzz_1,  \
                             ta1_x_xxx_yyzzz_0,   \
                             ta1_x_xxx_yyzzz_1,   \
                             ta1_x_xxx_yyzzzz_0,  \
                             ta1_x_xxx_yyzzzz_1,  \
                             ta1_x_xxx_yzzzz_0,   \
                             ta1_x_xxx_yzzzz_1,   \
                             ta1_x_xxx_yzzzzz_0,  \
                             ta1_x_xxx_yzzzzz_1,  \
                             ta1_x_xxx_zzzzz_0,   \
                             ta1_x_xxx_zzzzz_1,   \
                             ta1_x_xxx_zzzzzz_0,  \
                             ta1_x_xxx_zzzzzz_1,  \
                             ta1_x_xxxy_xxxxxx_0, \
                             ta1_x_xxxy_xxxxxy_0, \
                             ta1_x_xxxy_xxxxxz_0, \
                             ta1_x_xxxy_xxxxyy_0, \
                             ta1_x_xxxy_xxxxyz_0, \
                             ta1_x_xxxy_xxxxzz_0, \
                             ta1_x_xxxy_xxxyyy_0, \
                             ta1_x_xxxy_xxxyyz_0, \
                             ta1_x_xxxy_xxxyzz_0, \
                             ta1_x_xxxy_xxxzzz_0, \
                             ta1_x_xxxy_xxyyyy_0, \
                             ta1_x_xxxy_xxyyyz_0, \
                             ta1_x_xxxy_xxyyzz_0, \
                             ta1_x_xxxy_xxyzzz_0, \
                             ta1_x_xxxy_xxzzzz_0, \
                             ta1_x_xxxy_xyyyyy_0, \
                             ta1_x_xxxy_xyyyyz_0, \
                             ta1_x_xxxy_xyyyzz_0, \
                             ta1_x_xxxy_xyyzzz_0, \
                             ta1_x_xxxy_xyzzzz_0, \
                             ta1_x_xxxy_xzzzzz_0, \
                             ta1_x_xxxy_yyyyyy_0, \
                             ta1_x_xxxy_yyyyyz_0, \
                             ta1_x_xxxy_yyyyzz_0, \
                             ta1_x_xxxy_yyyzzz_0, \
                             ta1_x_xxxy_yyzzzz_0, \
                             ta1_x_xxxy_yzzzzz_0, \
                             ta1_x_xxxy_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxy_xxxxxx_0[i] = ta1_x_xxx_xxxxxx_0[i] * pa_y[i] - ta1_x_xxx_xxxxxx_1[i] * pc_y[i];

        ta1_x_xxxy_xxxxxy_0[i] =
            ta1_x_xxx_xxxxx_0[i] * fe_0 - ta1_x_xxx_xxxxx_1[i] * fe_0 + ta1_x_xxx_xxxxxy_0[i] * pa_y[i] - ta1_x_xxx_xxxxxy_1[i] * pc_y[i];

        ta1_x_xxxy_xxxxxz_0[i] = ta1_x_xxx_xxxxxz_0[i] * pa_y[i] - ta1_x_xxx_xxxxxz_1[i] * pc_y[i];

        ta1_x_xxxy_xxxxyy_0[i] =
            2.0 * ta1_x_xxx_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxxxy_1[i] * fe_0 + ta1_x_xxx_xxxxyy_0[i] * pa_y[i] - ta1_x_xxx_xxxxyy_1[i] * pc_y[i];

        ta1_x_xxxy_xxxxyz_0[i] =
            ta1_x_xxx_xxxxz_0[i] * fe_0 - ta1_x_xxx_xxxxz_1[i] * fe_0 + ta1_x_xxx_xxxxyz_0[i] * pa_y[i] - ta1_x_xxx_xxxxyz_1[i] * pc_y[i];

        ta1_x_xxxy_xxxxzz_0[i] = ta1_x_xxx_xxxxzz_0[i] * pa_y[i] - ta1_x_xxx_xxxxzz_1[i] * pc_y[i];

        ta1_x_xxxy_xxxyyy_0[i] =
            3.0 * ta1_x_xxx_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_xxx_xxxyy_1[i] * fe_0 + ta1_x_xxx_xxxyyy_0[i] * pa_y[i] - ta1_x_xxx_xxxyyy_1[i] * pc_y[i];

        ta1_x_xxxy_xxxyyz_0[i] =
            2.0 * ta1_x_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxxyz_1[i] * fe_0 + ta1_x_xxx_xxxyyz_0[i] * pa_y[i] - ta1_x_xxx_xxxyyz_1[i] * pc_y[i];

        ta1_x_xxxy_xxxyzz_0[i] =
            ta1_x_xxx_xxxzz_0[i] * fe_0 - ta1_x_xxx_xxxzz_1[i] * fe_0 + ta1_x_xxx_xxxyzz_0[i] * pa_y[i] - ta1_x_xxx_xxxyzz_1[i] * pc_y[i];

        ta1_x_xxxy_xxxzzz_0[i] = ta1_x_xxx_xxxzzz_0[i] * pa_y[i] - ta1_x_xxx_xxxzzz_1[i] * pc_y[i];

        ta1_x_xxxy_xxyyyy_0[i] =
            4.0 * ta1_x_xxx_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxyyy_1[i] * fe_0 + ta1_x_xxx_xxyyyy_0[i] * pa_y[i] - ta1_x_xxx_xxyyyy_1[i] * pc_y[i];

        ta1_x_xxxy_xxyyyz_0[i] =
            3.0 * ta1_x_xxx_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xxyyz_1[i] * fe_0 + ta1_x_xxx_xxyyyz_0[i] * pa_y[i] - ta1_x_xxx_xxyyyz_1[i] * pc_y[i];

        ta1_x_xxxy_xxyyzz_0[i] =
            2.0 * ta1_x_xxx_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxyzz_1[i] * fe_0 + ta1_x_xxx_xxyyzz_0[i] * pa_y[i] - ta1_x_xxx_xxyyzz_1[i] * pc_y[i];

        ta1_x_xxxy_xxyzzz_0[i] =
            ta1_x_xxx_xxzzz_0[i] * fe_0 - ta1_x_xxx_xxzzz_1[i] * fe_0 + ta1_x_xxx_xxyzzz_0[i] * pa_y[i] - ta1_x_xxx_xxyzzz_1[i] * pc_y[i];

        ta1_x_xxxy_xxzzzz_0[i] = ta1_x_xxx_xxzzzz_0[i] * pa_y[i] - ta1_x_xxx_xxzzzz_1[i] * pc_y[i];

        ta1_x_xxxy_xyyyyy_0[i] =
            5.0 * ta1_x_xxx_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_xxx_xyyyy_1[i] * fe_0 + ta1_x_xxx_xyyyyy_0[i] * pa_y[i] - ta1_x_xxx_xyyyyy_1[i] * pc_y[i];

        ta1_x_xxxy_xyyyyz_0[i] =
            4.0 * ta1_x_xxx_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyyyz_1[i] * fe_0 + ta1_x_xxx_xyyyyz_0[i] * pa_y[i] - ta1_x_xxx_xyyyyz_1[i] * pc_y[i];

        ta1_x_xxxy_xyyyzz_0[i] =
            3.0 * ta1_x_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xyyzz_1[i] * fe_0 + ta1_x_xxx_xyyyzz_0[i] * pa_y[i] - ta1_x_xxx_xyyyzz_1[i] * pc_y[i];

        ta1_x_xxxy_xyyzzz_0[i] =
            2.0 * ta1_x_xxx_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xyzzz_1[i] * fe_0 + ta1_x_xxx_xyyzzz_0[i] * pa_y[i] - ta1_x_xxx_xyyzzz_1[i] * pc_y[i];

        ta1_x_xxxy_xyzzzz_0[i] =
            ta1_x_xxx_xzzzz_0[i] * fe_0 - ta1_x_xxx_xzzzz_1[i] * fe_0 + ta1_x_xxx_xyzzzz_0[i] * pa_y[i] - ta1_x_xxx_xyzzzz_1[i] * pc_y[i];

        ta1_x_xxxy_xzzzzz_0[i] = ta1_x_xxx_xzzzzz_0[i] * pa_y[i] - ta1_x_xxx_xzzzzz_1[i] * pc_y[i];

        ta1_x_xxxy_yyyyyy_0[i] =
            6.0 * ta1_x_xxx_yyyyy_0[i] * fe_0 - 6.0 * ta1_x_xxx_yyyyy_1[i] * fe_0 + ta1_x_xxx_yyyyyy_0[i] * pa_y[i] - ta1_x_xxx_yyyyyy_1[i] * pc_y[i];

        ta1_x_xxxy_yyyyyz_0[i] =
            5.0 * ta1_x_xxx_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_xxx_yyyyz_1[i] * fe_0 + ta1_x_xxx_yyyyyz_0[i] * pa_y[i] - ta1_x_xxx_yyyyyz_1[i] * pc_y[i];

        ta1_x_xxxy_yyyyzz_0[i] =
            4.0 * ta1_x_xxx_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyyzz_1[i] * fe_0 + ta1_x_xxx_yyyyzz_0[i] * pa_y[i] - ta1_x_xxx_yyyyzz_1[i] * pc_y[i];

        ta1_x_xxxy_yyyzzz_0[i] =
            3.0 * ta1_x_xxx_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_yyzzz_1[i] * fe_0 + ta1_x_xxx_yyyzzz_0[i] * pa_y[i] - ta1_x_xxx_yyyzzz_1[i] * pc_y[i];

        ta1_x_xxxy_yyzzzz_0[i] =
            2.0 * ta1_x_xxx_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yzzzz_1[i] * fe_0 + ta1_x_xxx_yyzzzz_0[i] * pa_y[i] - ta1_x_xxx_yyzzzz_1[i] * pc_y[i];

        ta1_x_xxxy_yzzzzz_0[i] =
            ta1_x_xxx_zzzzz_0[i] * fe_0 - ta1_x_xxx_zzzzz_1[i] * fe_0 + ta1_x_xxx_yzzzzz_0[i] * pa_y[i] - ta1_x_xxx_yzzzzz_1[i] * pc_y[i];

        ta1_x_xxxy_zzzzzz_0[i] = ta1_x_xxx_zzzzzz_0[i] * pa_y[i] - ta1_x_xxx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : GI

    auto ta1_x_xxxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 56);

    auto ta1_x_xxxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 57);

    auto ta1_x_xxxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 58);

    auto ta1_x_xxxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 59);

    auto ta1_x_xxxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 60);

    auto ta1_x_xxxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 61);

    auto ta1_x_xxxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 62);

    auto ta1_x_xxxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 63);

    auto ta1_x_xxxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 64);

    auto ta1_x_xxxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 65);

    auto ta1_x_xxxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 66);

    auto ta1_x_xxxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 67);

    auto ta1_x_xxxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 68);

    auto ta1_x_xxxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 69);

    auto ta1_x_xxxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 70);

    auto ta1_x_xxxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 71);

    auto ta1_x_xxxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 72);

    auto ta1_x_xxxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 73);

    auto ta1_x_xxxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 74);

    auto ta1_x_xxxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 75);

    auto ta1_x_xxxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 76);

    auto ta1_x_xxxz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 77);

    auto ta1_x_xxxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 78);

    auto ta1_x_xxxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 79);

    auto ta1_x_xxxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 80);

    auto ta1_x_xxxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 81);

    auto ta1_x_xxxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 82);

    auto ta1_x_xxxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 83);

#pragma omp simd aligned(pa_z,                    \
                             pc_z,                \
                             ta1_x_xxx_xxxxx_0,   \
                             ta1_x_xxx_xxxxx_1,   \
                             ta1_x_xxx_xxxxxx_0,  \
                             ta1_x_xxx_xxxxxx_1,  \
                             ta1_x_xxx_xxxxxy_0,  \
                             ta1_x_xxx_xxxxxy_1,  \
                             ta1_x_xxx_xxxxxz_0,  \
                             ta1_x_xxx_xxxxxz_1,  \
                             ta1_x_xxx_xxxxy_0,   \
                             ta1_x_xxx_xxxxy_1,   \
                             ta1_x_xxx_xxxxyy_0,  \
                             ta1_x_xxx_xxxxyy_1,  \
                             ta1_x_xxx_xxxxyz_0,  \
                             ta1_x_xxx_xxxxyz_1,  \
                             ta1_x_xxx_xxxxz_0,   \
                             ta1_x_xxx_xxxxz_1,   \
                             ta1_x_xxx_xxxxzz_0,  \
                             ta1_x_xxx_xxxxzz_1,  \
                             ta1_x_xxx_xxxyy_0,   \
                             ta1_x_xxx_xxxyy_1,   \
                             ta1_x_xxx_xxxyyy_0,  \
                             ta1_x_xxx_xxxyyy_1,  \
                             ta1_x_xxx_xxxyyz_0,  \
                             ta1_x_xxx_xxxyyz_1,  \
                             ta1_x_xxx_xxxyz_0,   \
                             ta1_x_xxx_xxxyz_1,   \
                             ta1_x_xxx_xxxyzz_0,  \
                             ta1_x_xxx_xxxyzz_1,  \
                             ta1_x_xxx_xxxzz_0,   \
                             ta1_x_xxx_xxxzz_1,   \
                             ta1_x_xxx_xxxzzz_0,  \
                             ta1_x_xxx_xxxzzz_1,  \
                             ta1_x_xxx_xxyyy_0,   \
                             ta1_x_xxx_xxyyy_1,   \
                             ta1_x_xxx_xxyyyy_0,  \
                             ta1_x_xxx_xxyyyy_1,  \
                             ta1_x_xxx_xxyyyz_0,  \
                             ta1_x_xxx_xxyyyz_1,  \
                             ta1_x_xxx_xxyyz_0,   \
                             ta1_x_xxx_xxyyz_1,   \
                             ta1_x_xxx_xxyyzz_0,  \
                             ta1_x_xxx_xxyyzz_1,  \
                             ta1_x_xxx_xxyzz_0,   \
                             ta1_x_xxx_xxyzz_1,   \
                             ta1_x_xxx_xxyzzz_0,  \
                             ta1_x_xxx_xxyzzz_1,  \
                             ta1_x_xxx_xxzzz_0,   \
                             ta1_x_xxx_xxzzz_1,   \
                             ta1_x_xxx_xxzzzz_0,  \
                             ta1_x_xxx_xxzzzz_1,  \
                             ta1_x_xxx_xyyyy_0,   \
                             ta1_x_xxx_xyyyy_1,   \
                             ta1_x_xxx_xyyyyy_0,  \
                             ta1_x_xxx_xyyyyy_1,  \
                             ta1_x_xxx_xyyyyz_0,  \
                             ta1_x_xxx_xyyyyz_1,  \
                             ta1_x_xxx_xyyyz_0,   \
                             ta1_x_xxx_xyyyz_1,   \
                             ta1_x_xxx_xyyyzz_0,  \
                             ta1_x_xxx_xyyyzz_1,  \
                             ta1_x_xxx_xyyzz_0,   \
                             ta1_x_xxx_xyyzz_1,   \
                             ta1_x_xxx_xyyzzz_0,  \
                             ta1_x_xxx_xyyzzz_1,  \
                             ta1_x_xxx_xyzzz_0,   \
                             ta1_x_xxx_xyzzz_1,   \
                             ta1_x_xxx_xyzzzz_0,  \
                             ta1_x_xxx_xyzzzz_1,  \
                             ta1_x_xxx_xzzzz_0,   \
                             ta1_x_xxx_xzzzz_1,   \
                             ta1_x_xxx_xzzzzz_0,  \
                             ta1_x_xxx_xzzzzz_1,  \
                             ta1_x_xxx_yyyyy_0,   \
                             ta1_x_xxx_yyyyy_1,   \
                             ta1_x_xxx_yyyyyy_0,  \
                             ta1_x_xxx_yyyyyy_1,  \
                             ta1_x_xxx_yyyyyz_0,  \
                             ta1_x_xxx_yyyyyz_1,  \
                             ta1_x_xxx_yyyyz_0,   \
                             ta1_x_xxx_yyyyz_1,   \
                             ta1_x_xxx_yyyyzz_0,  \
                             ta1_x_xxx_yyyyzz_1,  \
                             ta1_x_xxx_yyyzz_0,   \
                             ta1_x_xxx_yyyzz_1,   \
                             ta1_x_xxx_yyyzzz_0,  \
                             ta1_x_xxx_yyyzzz_1,  \
                             ta1_x_xxx_yyzzz_0,   \
                             ta1_x_xxx_yyzzz_1,   \
                             ta1_x_xxx_yyzzzz_0,  \
                             ta1_x_xxx_yyzzzz_1,  \
                             ta1_x_xxx_yzzzz_0,   \
                             ta1_x_xxx_yzzzz_1,   \
                             ta1_x_xxx_yzzzzz_0,  \
                             ta1_x_xxx_yzzzzz_1,  \
                             ta1_x_xxx_zzzzz_0,   \
                             ta1_x_xxx_zzzzz_1,   \
                             ta1_x_xxx_zzzzzz_0,  \
                             ta1_x_xxx_zzzzzz_1,  \
                             ta1_x_xxxz_xxxxxx_0, \
                             ta1_x_xxxz_xxxxxy_0, \
                             ta1_x_xxxz_xxxxxz_0, \
                             ta1_x_xxxz_xxxxyy_0, \
                             ta1_x_xxxz_xxxxyz_0, \
                             ta1_x_xxxz_xxxxzz_0, \
                             ta1_x_xxxz_xxxyyy_0, \
                             ta1_x_xxxz_xxxyyz_0, \
                             ta1_x_xxxz_xxxyzz_0, \
                             ta1_x_xxxz_xxxzzz_0, \
                             ta1_x_xxxz_xxyyyy_0, \
                             ta1_x_xxxz_xxyyyz_0, \
                             ta1_x_xxxz_xxyyzz_0, \
                             ta1_x_xxxz_xxyzzz_0, \
                             ta1_x_xxxz_xxzzzz_0, \
                             ta1_x_xxxz_xyyyyy_0, \
                             ta1_x_xxxz_xyyyyz_0, \
                             ta1_x_xxxz_xyyyzz_0, \
                             ta1_x_xxxz_xyyzzz_0, \
                             ta1_x_xxxz_xyzzzz_0, \
                             ta1_x_xxxz_xzzzzz_0, \
                             ta1_x_xxxz_yyyyyy_0, \
                             ta1_x_xxxz_yyyyyz_0, \
                             ta1_x_xxxz_yyyyzz_0, \
                             ta1_x_xxxz_yyyzzz_0, \
                             ta1_x_xxxz_yyzzzz_0, \
                             ta1_x_xxxz_yzzzzz_0, \
                             ta1_x_xxxz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxz_xxxxxx_0[i] = ta1_x_xxx_xxxxxx_0[i] * pa_z[i] - ta1_x_xxx_xxxxxx_1[i] * pc_z[i];

        ta1_x_xxxz_xxxxxy_0[i] = ta1_x_xxx_xxxxxy_0[i] * pa_z[i] - ta1_x_xxx_xxxxxy_1[i] * pc_z[i];

        ta1_x_xxxz_xxxxxz_0[i] =
            ta1_x_xxx_xxxxx_0[i] * fe_0 - ta1_x_xxx_xxxxx_1[i] * fe_0 + ta1_x_xxx_xxxxxz_0[i] * pa_z[i] - ta1_x_xxx_xxxxxz_1[i] * pc_z[i];

        ta1_x_xxxz_xxxxyy_0[i] = ta1_x_xxx_xxxxyy_0[i] * pa_z[i] - ta1_x_xxx_xxxxyy_1[i] * pc_z[i];

        ta1_x_xxxz_xxxxyz_0[i] =
            ta1_x_xxx_xxxxy_0[i] * fe_0 - ta1_x_xxx_xxxxy_1[i] * fe_0 + ta1_x_xxx_xxxxyz_0[i] * pa_z[i] - ta1_x_xxx_xxxxyz_1[i] * pc_z[i];

        ta1_x_xxxz_xxxxzz_0[i] =
            2.0 * ta1_x_xxx_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxxxz_1[i] * fe_0 + ta1_x_xxx_xxxxzz_0[i] * pa_z[i] - ta1_x_xxx_xxxxzz_1[i] * pc_z[i];

        ta1_x_xxxz_xxxyyy_0[i] = ta1_x_xxx_xxxyyy_0[i] * pa_z[i] - ta1_x_xxx_xxxyyy_1[i] * pc_z[i];

        ta1_x_xxxz_xxxyyz_0[i] =
            ta1_x_xxx_xxxyy_0[i] * fe_0 - ta1_x_xxx_xxxyy_1[i] * fe_0 + ta1_x_xxx_xxxyyz_0[i] * pa_z[i] - ta1_x_xxx_xxxyyz_1[i] * pc_z[i];

        ta1_x_xxxz_xxxyzz_0[i] =
            2.0 * ta1_x_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxxyz_1[i] * fe_0 + ta1_x_xxx_xxxyzz_0[i] * pa_z[i] - ta1_x_xxx_xxxyzz_1[i] * pc_z[i];

        ta1_x_xxxz_xxxzzz_0[i] =
            3.0 * ta1_x_xxx_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xxxzz_1[i] * fe_0 + ta1_x_xxx_xxxzzz_0[i] * pa_z[i] - ta1_x_xxx_xxxzzz_1[i] * pc_z[i];

        ta1_x_xxxz_xxyyyy_0[i] = ta1_x_xxx_xxyyyy_0[i] * pa_z[i] - ta1_x_xxx_xxyyyy_1[i] * pc_z[i];

        ta1_x_xxxz_xxyyyz_0[i] =
            ta1_x_xxx_xxyyy_0[i] * fe_0 - ta1_x_xxx_xxyyy_1[i] * fe_0 + ta1_x_xxx_xxyyyz_0[i] * pa_z[i] - ta1_x_xxx_xxyyyz_1[i] * pc_z[i];

        ta1_x_xxxz_xxyyzz_0[i] =
            2.0 * ta1_x_xxx_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxyyz_1[i] * fe_0 + ta1_x_xxx_xxyyzz_0[i] * pa_z[i] - ta1_x_xxx_xxyyzz_1[i] * pc_z[i];

        ta1_x_xxxz_xxyzzz_0[i] =
            3.0 * ta1_x_xxx_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xxyzz_1[i] * fe_0 + ta1_x_xxx_xxyzzz_0[i] * pa_z[i] - ta1_x_xxx_xxyzzz_1[i] * pc_z[i];

        ta1_x_xxxz_xxzzzz_0[i] =
            4.0 * ta1_x_xxx_xxzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xxzzz_1[i] * fe_0 + ta1_x_xxx_xxzzzz_0[i] * pa_z[i] - ta1_x_xxx_xxzzzz_1[i] * pc_z[i];

        ta1_x_xxxz_xyyyyy_0[i] = ta1_x_xxx_xyyyyy_0[i] * pa_z[i] - ta1_x_xxx_xyyyyy_1[i] * pc_z[i];

        ta1_x_xxxz_xyyyyz_0[i] =
            ta1_x_xxx_xyyyy_0[i] * fe_0 - ta1_x_xxx_xyyyy_1[i] * fe_0 + ta1_x_xxx_xyyyyz_0[i] * pa_z[i] - ta1_x_xxx_xyyyyz_1[i] * pc_z[i];

        ta1_x_xxxz_xyyyzz_0[i] =
            2.0 * ta1_x_xxx_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xyyyz_1[i] * fe_0 + ta1_x_xxx_xyyyzz_0[i] * pa_z[i] - ta1_x_xxx_xyyyzz_1[i] * pc_z[i];

        ta1_x_xxxz_xyyzzz_0[i] =
            3.0 * ta1_x_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xyyzz_1[i] * fe_0 + ta1_x_xxx_xyyzzz_0[i] * pa_z[i] - ta1_x_xxx_xyyzzz_1[i] * pc_z[i];

        ta1_x_xxxz_xyzzzz_0[i] =
            4.0 * ta1_x_xxx_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyzzz_1[i] * fe_0 + ta1_x_xxx_xyzzzz_0[i] * pa_z[i] - ta1_x_xxx_xyzzzz_1[i] * pc_z[i];

        ta1_x_xxxz_xzzzzz_0[i] =
            5.0 * ta1_x_xxx_xzzzz_0[i] * fe_0 - 5.0 * ta1_x_xxx_xzzzz_1[i] * fe_0 + ta1_x_xxx_xzzzzz_0[i] * pa_z[i] - ta1_x_xxx_xzzzzz_1[i] * pc_z[i];

        ta1_x_xxxz_yyyyyy_0[i] = ta1_x_xxx_yyyyyy_0[i] * pa_z[i] - ta1_x_xxx_yyyyyy_1[i] * pc_z[i];

        ta1_x_xxxz_yyyyyz_0[i] =
            ta1_x_xxx_yyyyy_0[i] * fe_0 - ta1_x_xxx_yyyyy_1[i] * fe_0 + ta1_x_xxx_yyyyyz_0[i] * pa_z[i] - ta1_x_xxx_yyyyyz_1[i] * pc_z[i];

        ta1_x_xxxz_yyyyzz_0[i] =
            2.0 * ta1_x_xxx_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yyyyz_1[i] * fe_0 + ta1_x_xxx_yyyyzz_0[i] * pa_z[i] - ta1_x_xxx_yyyyzz_1[i] * pc_z[i];

        ta1_x_xxxz_yyyzzz_0[i] =
            3.0 * ta1_x_xxx_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_yyyzz_1[i] * fe_0 + ta1_x_xxx_yyyzzz_0[i] * pa_z[i] - ta1_x_xxx_yyyzzz_1[i] * pc_z[i];

        ta1_x_xxxz_yyzzzz_0[i] =
            4.0 * ta1_x_xxx_yyzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyzzz_1[i] * fe_0 + ta1_x_xxx_yyzzzz_0[i] * pa_z[i] - ta1_x_xxx_yyzzzz_1[i] * pc_z[i];

        ta1_x_xxxz_yzzzzz_0[i] =
            5.0 * ta1_x_xxx_yzzzz_0[i] * fe_0 - 5.0 * ta1_x_xxx_yzzzz_1[i] * fe_0 + ta1_x_xxx_yzzzzz_0[i] * pa_z[i] - ta1_x_xxx_yzzzzz_1[i] * pc_z[i];

        ta1_x_xxxz_zzzzzz_0[i] =
            6.0 * ta1_x_xxx_zzzzz_0[i] * fe_0 - 6.0 * ta1_x_xxx_zzzzz_1[i] * fe_0 + ta1_x_xxx_zzzzzz_0[i] * pa_z[i] - ta1_x_xxx_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 84-112 components of targeted buffer : GI

    auto ta1_x_xxyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 84);

    auto ta1_x_xxyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 85);

    auto ta1_x_xxyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 86);

    auto ta1_x_xxyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 87);

    auto ta1_x_xxyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 88);

    auto ta1_x_xxyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 89);

    auto ta1_x_xxyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 90);

    auto ta1_x_xxyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 91);

    auto ta1_x_xxyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 92);

    auto ta1_x_xxyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 93);

    auto ta1_x_xxyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 94);

    auto ta1_x_xxyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 95);

    auto ta1_x_xxyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 96);

    auto ta1_x_xxyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 97);

    auto ta1_x_xxyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 98);

    auto ta1_x_xxyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 99);

    auto ta1_x_xxyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 100);

    auto ta1_x_xxyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 101);

    auto ta1_x_xxyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 102);

    auto ta1_x_xxyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 103);

    auto ta1_x_xxyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 104);

    auto ta1_x_xxyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 105);

    auto ta1_x_xxyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 106);

    auto ta1_x_xxyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 107);

    auto ta1_x_xxyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 108);

    auto ta1_x_xxyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 109);

    auto ta1_x_xxyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 110);

    auto ta1_x_xxyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 111);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_x_xx_xxxxxx_0,   \
                             ta1_x_xx_xxxxxx_1,   \
                             ta1_x_xx_xxxxxy_0,   \
                             ta1_x_xx_xxxxxy_1,   \
                             ta1_x_xx_xxxxxz_0,   \
                             ta1_x_xx_xxxxxz_1,   \
                             ta1_x_xx_xxxxyy_0,   \
                             ta1_x_xx_xxxxyy_1,   \
                             ta1_x_xx_xxxxyz_0,   \
                             ta1_x_xx_xxxxyz_1,   \
                             ta1_x_xx_xxxxzz_0,   \
                             ta1_x_xx_xxxxzz_1,   \
                             ta1_x_xx_xxxyyy_0,   \
                             ta1_x_xx_xxxyyy_1,   \
                             ta1_x_xx_xxxyyz_0,   \
                             ta1_x_xx_xxxyyz_1,   \
                             ta1_x_xx_xxxyzz_0,   \
                             ta1_x_xx_xxxyzz_1,   \
                             ta1_x_xx_xxxzzz_0,   \
                             ta1_x_xx_xxxzzz_1,   \
                             ta1_x_xx_xxyyyy_0,   \
                             ta1_x_xx_xxyyyy_1,   \
                             ta1_x_xx_xxyyyz_0,   \
                             ta1_x_xx_xxyyyz_1,   \
                             ta1_x_xx_xxyyzz_0,   \
                             ta1_x_xx_xxyyzz_1,   \
                             ta1_x_xx_xxyzzz_0,   \
                             ta1_x_xx_xxyzzz_1,   \
                             ta1_x_xx_xxzzzz_0,   \
                             ta1_x_xx_xxzzzz_1,   \
                             ta1_x_xx_xyyyyy_0,   \
                             ta1_x_xx_xyyyyy_1,   \
                             ta1_x_xx_xyyyyz_0,   \
                             ta1_x_xx_xyyyyz_1,   \
                             ta1_x_xx_xyyyzz_0,   \
                             ta1_x_xx_xyyyzz_1,   \
                             ta1_x_xx_xyyzzz_0,   \
                             ta1_x_xx_xyyzzz_1,   \
                             ta1_x_xx_xyzzzz_0,   \
                             ta1_x_xx_xyzzzz_1,   \
                             ta1_x_xx_xzzzzz_0,   \
                             ta1_x_xx_xzzzzz_1,   \
                             ta1_x_xx_zzzzzz_0,   \
                             ta1_x_xx_zzzzzz_1,   \
                             ta1_x_xxy_xxxxx_0,   \
                             ta1_x_xxy_xxxxx_1,   \
                             ta1_x_xxy_xxxxxx_0,  \
                             ta1_x_xxy_xxxxxx_1,  \
                             ta1_x_xxy_xxxxxy_0,  \
                             ta1_x_xxy_xxxxxy_1,  \
                             ta1_x_xxy_xxxxxz_0,  \
                             ta1_x_xxy_xxxxxz_1,  \
                             ta1_x_xxy_xxxxy_0,   \
                             ta1_x_xxy_xxxxy_1,   \
                             ta1_x_xxy_xxxxyy_0,  \
                             ta1_x_xxy_xxxxyy_1,  \
                             ta1_x_xxy_xxxxyz_0,  \
                             ta1_x_xxy_xxxxyz_1,  \
                             ta1_x_xxy_xxxxz_0,   \
                             ta1_x_xxy_xxxxz_1,   \
                             ta1_x_xxy_xxxxzz_0,  \
                             ta1_x_xxy_xxxxzz_1,  \
                             ta1_x_xxy_xxxyy_0,   \
                             ta1_x_xxy_xxxyy_1,   \
                             ta1_x_xxy_xxxyyy_0,  \
                             ta1_x_xxy_xxxyyy_1,  \
                             ta1_x_xxy_xxxyyz_0,  \
                             ta1_x_xxy_xxxyyz_1,  \
                             ta1_x_xxy_xxxyz_0,   \
                             ta1_x_xxy_xxxyz_1,   \
                             ta1_x_xxy_xxxyzz_0,  \
                             ta1_x_xxy_xxxyzz_1,  \
                             ta1_x_xxy_xxxzz_0,   \
                             ta1_x_xxy_xxxzz_1,   \
                             ta1_x_xxy_xxxzzz_0,  \
                             ta1_x_xxy_xxxzzz_1,  \
                             ta1_x_xxy_xxyyy_0,   \
                             ta1_x_xxy_xxyyy_1,   \
                             ta1_x_xxy_xxyyyy_0,  \
                             ta1_x_xxy_xxyyyy_1,  \
                             ta1_x_xxy_xxyyyz_0,  \
                             ta1_x_xxy_xxyyyz_1,  \
                             ta1_x_xxy_xxyyz_0,   \
                             ta1_x_xxy_xxyyz_1,   \
                             ta1_x_xxy_xxyyzz_0,  \
                             ta1_x_xxy_xxyyzz_1,  \
                             ta1_x_xxy_xxyzz_0,   \
                             ta1_x_xxy_xxyzz_1,   \
                             ta1_x_xxy_xxyzzz_0,  \
                             ta1_x_xxy_xxyzzz_1,  \
                             ta1_x_xxy_xxzzz_0,   \
                             ta1_x_xxy_xxzzz_1,   \
                             ta1_x_xxy_xxzzzz_0,  \
                             ta1_x_xxy_xxzzzz_1,  \
                             ta1_x_xxy_xyyyy_0,   \
                             ta1_x_xxy_xyyyy_1,   \
                             ta1_x_xxy_xyyyyy_0,  \
                             ta1_x_xxy_xyyyyy_1,  \
                             ta1_x_xxy_xyyyyz_0,  \
                             ta1_x_xxy_xyyyyz_1,  \
                             ta1_x_xxy_xyyyz_0,   \
                             ta1_x_xxy_xyyyz_1,   \
                             ta1_x_xxy_xyyyzz_0,  \
                             ta1_x_xxy_xyyyzz_1,  \
                             ta1_x_xxy_xyyzz_0,   \
                             ta1_x_xxy_xyyzz_1,   \
                             ta1_x_xxy_xyyzzz_0,  \
                             ta1_x_xxy_xyyzzz_1,  \
                             ta1_x_xxy_xyzzz_0,   \
                             ta1_x_xxy_xyzzz_1,   \
                             ta1_x_xxy_xyzzzz_0,  \
                             ta1_x_xxy_xyzzzz_1,  \
                             ta1_x_xxy_xzzzz_0,   \
                             ta1_x_xxy_xzzzz_1,   \
                             ta1_x_xxy_xzzzzz_0,  \
                             ta1_x_xxy_xzzzzz_1,  \
                             ta1_x_xxy_zzzzzz_0,  \
                             ta1_x_xxy_zzzzzz_1,  \
                             ta1_x_xxyy_xxxxxx_0, \
                             ta1_x_xxyy_xxxxxy_0, \
                             ta1_x_xxyy_xxxxxz_0, \
                             ta1_x_xxyy_xxxxyy_0, \
                             ta1_x_xxyy_xxxxyz_0, \
                             ta1_x_xxyy_xxxxzz_0, \
                             ta1_x_xxyy_xxxyyy_0, \
                             ta1_x_xxyy_xxxyyz_0, \
                             ta1_x_xxyy_xxxyzz_0, \
                             ta1_x_xxyy_xxxzzz_0, \
                             ta1_x_xxyy_xxyyyy_0, \
                             ta1_x_xxyy_xxyyyz_0, \
                             ta1_x_xxyy_xxyyzz_0, \
                             ta1_x_xxyy_xxyzzz_0, \
                             ta1_x_xxyy_xxzzzz_0, \
                             ta1_x_xxyy_xyyyyy_0, \
                             ta1_x_xxyy_xyyyyz_0, \
                             ta1_x_xxyy_xyyyzz_0, \
                             ta1_x_xxyy_xyyzzz_0, \
                             ta1_x_xxyy_xyzzzz_0, \
                             ta1_x_xxyy_xzzzzz_0, \
                             ta1_x_xxyy_yyyyyy_0, \
                             ta1_x_xxyy_yyyyyz_0, \
                             ta1_x_xxyy_yyyyzz_0, \
                             ta1_x_xxyy_yyyzzz_0, \
                             ta1_x_xxyy_yyzzzz_0, \
                             ta1_x_xxyy_yzzzzz_0, \
                             ta1_x_xxyy_zzzzzz_0, \
                             ta1_x_xyy_yyyyyy_0,  \
                             ta1_x_xyy_yyyyyy_1,  \
                             ta1_x_xyy_yyyyyz_0,  \
                             ta1_x_xyy_yyyyyz_1,  \
                             ta1_x_xyy_yyyyzz_0,  \
                             ta1_x_xyy_yyyyzz_1,  \
                             ta1_x_xyy_yyyzzz_0,  \
                             ta1_x_xyy_yyyzzz_1,  \
                             ta1_x_xyy_yyzzzz_0,  \
                             ta1_x_xyy_yyzzzz_1,  \
                             ta1_x_xyy_yzzzzz_0,  \
                             ta1_x_xyy_yzzzzz_1,  \
                             ta1_x_yy_yyyyyy_0,   \
                             ta1_x_yy_yyyyyy_1,   \
                             ta1_x_yy_yyyyyz_0,   \
                             ta1_x_yy_yyyyyz_1,   \
                             ta1_x_yy_yyyyzz_0,   \
                             ta1_x_yy_yyyyzz_1,   \
                             ta1_x_yy_yyyzzz_0,   \
                             ta1_x_yy_yyyzzz_1,   \
                             ta1_x_yy_yyzzzz_0,   \
                             ta1_x_yy_yyzzzz_1,   \
                             ta1_x_yy_yzzzzz_0,   \
                             ta1_x_yy_yzzzzz_1,   \
                             ta_xyy_yyyyyy_1,     \
                             ta_xyy_yyyyyz_1,     \
                             ta_xyy_yyyyzz_1,     \
                             ta_xyy_yyyzzz_1,     \
                             ta_xyy_yyzzzz_1,     \
                             ta_xyy_yzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyy_xxxxxx_0[i] =
            ta1_x_xx_xxxxxx_0[i] * fe_0 - ta1_x_xx_xxxxxx_1[i] * fe_0 + ta1_x_xxy_xxxxxx_0[i] * pa_y[i] - ta1_x_xxy_xxxxxx_1[i] * pc_y[i];

        ta1_x_xxyy_xxxxxy_0[i] = ta1_x_xx_xxxxxy_0[i] * fe_0 - ta1_x_xx_xxxxxy_1[i] * fe_0 + ta1_x_xxy_xxxxx_0[i] * fe_0 -
                                 ta1_x_xxy_xxxxx_1[i] * fe_0 + ta1_x_xxy_xxxxxy_0[i] * pa_y[i] - ta1_x_xxy_xxxxxy_1[i] * pc_y[i];

        ta1_x_xxyy_xxxxxz_0[i] =
            ta1_x_xx_xxxxxz_0[i] * fe_0 - ta1_x_xx_xxxxxz_1[i] * fe_0 + ta1_x_xxy_xxxxxz_0[i] * pa_y[i] - ta1_x_xxy_xxxxxz_1[i] * pc_y[i];

        ta1_x_xxyy_xxxxyy_0[i] = ta1_x_xx_xxxxyy_0[i] * fe_0 - ta1_x_xx_xxxxyy_1[i] * fe_0 + 2.0 * ta1_x_xxy_xxxxy_0[i] * fe_0 -
                                 2.0 * ta1_x_xxy_xxxxy_1[i] * fe_0 + ta1_x_xxy_xxxxyy_0[i] * pa_y[i] - ta1_x_xxy_xxxxyy_1[i] * pc_y[i];

        ta1_x_xxyy_xxxxyz_0[i] = ta1_x_xx_xxxxyz_0[i] * fe_0 - ta1_x_xx_xxxxyz_1[i] * fe_0 + ta1_x_xxy_xxxxz_0[i] * fe_0 -
                                 ta1_x_xxy_xxxxz_1[i] * fe_0 + ta1_x_xxy_xxxxyz_0[i] * pa_y[i] - ta1_x_xxy_xxxxyz_1[i] * pc_y[i];

        ta1_x_xxyy_xxxxzz_0[i] =
            ta1_x_xx_xxxxzz_0[i] * fe_0 - ta1_x_xx_xxxxzz_1[i] * fe_0 + ta1_x_xxy_xxxxzz_0[i] * pa_y[i] - ta1_x_xxy_xxxxzz_1[i] * pc_y[i];

        ta1_x_xxyy_xxxyyy_0[i] = ta1_x_xx_xxxyyy_0[i] * fe_0 - ta1_x_xx_xxxyyy_1[i] * fe_0 + 3.0 * ta1_x_xxy_xxxyy_0[i] * fe_0 -
                                 3.0 * ta1_x_xxy_xxxyy_1[i] * fe_0 + ta1_x_xxy_xxxyyy_0[i] * pa_y[i] - ta1_x_xxy_xxxyyy_1[i] * pc_y[i];

        ta1_x_xxyy_xxxyyz_0[i] = ta1_x_xx_xxxyyz_0[i] * fe_0 - ta1_x_xx_xxxyyz_1[i] * fe_0 + 2.0 * ta1_x_xxy_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxy_xxxyz_1[i] * fe_0 + ta1_x_xxy_xxxyyz_0[i] * pa_y[i] - ta1_x_xxy_xxxyyz_1[i] * pc_y[i];

        ta1_x_xxyy_xxxyzz_0[i] = ta1_x_xx_xxxyzz_0[i] * fe_0 - ta1_x_xx_xxxyzz_1[i] * fe_0 + ta1_x_xxy_xxxzz_0[i] * fe_0 -
                                 ta1_x_xxy_xxxzz_1[i] * fe_0 + ta1_x_xxy_xxxyzz_0[i] * pa_y[i] - ta1_x_xxy_xxxyzz_1[i] * pc_y[i];

        ta1_x_xxyy_xxxzzz_0[i] =
            ta1_x_xx_xxxzzz_0[i] * fe_0 - ta1_x_xx_xxxzzz_1[i] * fe_0 + ta1_x_xxy_xxxzzz_0[i] * pa_y[i] - ta1_x_xxy_xxxzzz_1[i] * pc_y[i];

        ta1_x_xxyy_xxyyyy_0[i] = ta1_x_xx_xxyyyy_0[i] * fe_0 - ta1_x_xx_xxyyyy_1[i] * fe_0 + 4.0 * ta1_x_xxy_xxyyy_0[i] * fe_0 -
                                 4.0 * ta1_x_xxy_xxyyy_1[i] * fe_0 + ta1_x_xxy_xxyyyy_0[i] * pa_y[i] - ta1_x_xxy_xxyyyy_1[i] * pc_y[i];

        ta1_x_xxyy_xxyyyz_0[i] = ta1_x_xx_xxyyyz_0[i] * fe_0 - ta1_x_xx_xxyyyz_1[i] * fe_0 + 3.0 * ta1_x_xxy_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxy_xxyyz_1[i] * fe_0 + ta1_x_xxy_xxyyyz_0[i] * pa_y[i] - ta1_x_xxy_xxyyyz_1[i] * pc_y[i];

        ta1_x_xxyy_xxyyzz_0[i] = ta1_x_xx_xxyyzz_0[i] * fe_0 - ta1_x_xx_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_xxy_xxyzz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxy_xxyzz_1[i] * fe_0 + ta1_x_xxy_xxyyzz_0[i] * pa_y[i] - ta1_x_xxy_xxyyzz_1[i] * pc_y[i];

        ta1_x_xxyy_xxyzzz_0[i] = ta1_x_xx_xxyzzz_0[i] * fe_0 - ta1_x_xx_xxyzzz_1[i] * fe_0 + ta1_x_xxy_xxzzz_0[i] * fe_0 -
                                 ta1_x_xxy_xxzzz_1[i] * fe_0 + ta1_x_xxy_xxyzzz_0[i] * pa_y[i] - ta1_x_xxy_xxyzzz_1[i] * pc_y[i];

        ta1_x_xxyy_xxzzzz_0[i] =
            ta1_x_xx_xxzzzz_0[i] * fe_0 - ta1_x_xx_xxzzzz_1[i] * fe_0 + ta1_x_xxy_xxzzzz_0[i] * pa_y[i] - ta1_x_xxy_xxzzzz_1[i] * pc_y[i];

        ta1_x_xxyy_xyyyyy_0[i] = ta1_x_xx_xyyyyy_0[i] * fe_0 - ta1_x_xx_xyyyyy_1[i] * fe_0 + 5.0 * ta1_x_xxy_xyyyy_0[i] * fe_0 -
                                 5.0 * ta1_x_xxy_xyyyy_1[i] * fe_0 + ta1_x_xxy_xyyyyy_0[i] * pa_y[i] - ta1_x_xxy_xyyyyy_1[i] * pc_y[i];

        ta1_x_xxyy_xyyyyz_0[i] = ta1_x_xx_xyyyyz_0[i] * fe_0 - ta1_x_xx_xyyyyz_1[i] * fe_0 + 4.0 * ta1_x_xxy_xyyyz_0[i] * fe_0 -
                                 4.0 * ta1_x_xxy_xyyyz_1[i] * fe_0 + ta1_x_xxy_xyyyyz_0[i] * pa_y[i] - ta1_x_xxy_xyyyyz_1[i] * pc_y[i];

        ta1_x_xxyy_xyyyzz_0[i] = ta1_x_xx_xyyyzz_0[i] * fe_0 - ta1_x_xx_xyyyzz_1[i] * fe_0 + 3.0 * ta1_x_xxy_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxy_xyyzz_1[i] * fe_0 + ta1_x_xxy_xyyyzz_0[i] * pa_y[i] - ta1_x_xxy_xyyyzz_1[i] * pc_y[i];

        ta1_x_xxyy_xyyzzz_0[i] = ta1_x_xx_xyyzzz_0[i] * fe_0 - ta1_x_xx_xyyzzz_1[i] * fe_0 + 2.0 * ta1_x_xxy_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxy_xyzzz_1[i] * fe_0 + ta1_x_xxy_xyyzzz_0[i] * pa_y[i] - ta1_x_xxy_xyyzzz_1[i] * pc_y[i];

        ta1_x_xxyy_xyzzzz_0[i] = ta1_x_xx_xyzzzz_0[i] * fe_0 - ta1_x_xx_xyzzzz_1[i] * fe_0 + ta1_x_xxy_xzzzz_0[i] * fe_0 -
                                 ta1_x_xxy_xzzzz_1[i] * fe_0 + ta1_x_xxy_xyzzzz_0[i] * pa_y[i] - ta1_x_xxy_xyzzzz_1[i] * pc_y[i];

        ta1_x_xxyy_xzzzzz_0[i] =
            ta1_x_xx_xzzzzz_0[i] * fe_0 - ta1_x_xx_xzzzzz_1[i] * fe_0 + ta1_x_xxy_xzzzzz_0[i] * pa_y[i] - ta1_x_xxy_xzzzzz_1[i] * pc_y[i];

        ta1_x_xxyy_yyyyyy_0[i] = ta1_x_yy_yyyyyy_0[i] * fe_0 - ta1_x_yy_yyyyyy_1[i] * fe_0 + ta_xyy_yyyyyy_1[i] + ta1_x_xyy_yyyyyy_0[i] * pa_x[i] -
                                 ta1_x_xyy_yyyyyy_1[i] * pc_x[i];

        ta1_x_xxyy_yyyyyz_0[i] = ta1_x_yy_yyyyyz_0[i] * fe_0 - ta1_x_yy_yyyyyz_1[i] * fe_0 + ta_xyy_yyyyyz_1[i] + ta1_x_xyy_yyyyyz_0[i] * pa_x[i] -
                                 ta1_x_xyy_yyyyyz_1[i] * pc_x[i];

        ta1_x_xxyy_yyyyzz_0[i] = ta1_x_yy_yyyyzz_0[i] * fe_0 - ta1_x_yy_yyyyzz_1[i] * fe_0 + ta_xyy_yyyyzz_1[i] + ta1_x_xyy_yyyyzz_0[i] * pa_x[i] -
                                 ta1_x_xyy_yyyyzz_1[i] * pc_x[i];

        ta1_x_xxyy_yyyzzz_0[i] = ta1_x_yy_yyyzzz_0[i] * fe_0 - ta1_x_yy_yyyzzz_1[i] * fe_0 + ta_xyy_yyyzzz_1[i] + ta1_x_xyy_yyyzzz_0[i] * pa_x[i] -
                                 ta1_x_xyy_yyyzzz_1[i] * pc_x[i];

        ta1_x_xxyy_yyzzzz_0[i] = ta1_x_yy_yyzzzz_0[i] * fe_0 - ta1_x_yy_yyzzzz_1[i] * fe_0 + ta_xyy_yyzzzz_1[i] + ta1_x_xyy_yyzzzz_0[i] * pa_x[i] -
                                 ta1_x_xyy_yyzzzz_1[i] * pc_x[i];

        ta1_x_xxyy_yzzzzz_0[i] = ta1_x_yy_yzzzzz_0[i] * fe_0 - ta1_x_yy_yzzzzz_1[i] * fe_0 + ta_xyy_yzzzzz_1[i] + ta1_x_xyy_yzzzzz_0[i] * pa_x[i] -
                                 ta1_x_xyy_yzzzzz_1[i] * pc_x[i];

        ta1_x_xxyy_zzzzzz_0[i] =
            ta1_x_xx_zzzzzz_0[i] * fe_0 - ta1_x_xx_zzzzzz_1[i] * fe_0 + ta1_x_xxy_zzzzzz_0[i] * pa_y[i] - ta1_x_xxy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 112-140 components of targeted buffer : GI

    auto ta1_x_xxyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 112);

    auto ta1_x_xxyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 113);

    auto ta1_x_xxyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 114);

    auto ta1_x_xxyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 115);

    auto ta1_x_xxyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 116);

    auto ta1_x_xxyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 117);

    auto ta1_x_xxyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 118);

    auto ta1_x_xxyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 119);

    auto ta1_x_xxyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 120);

    auto ta1_x_xxyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 121);

    auto ta1_x_xxyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 122);

    auto ta1_x_xxyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 123);

    auto ta1_x_xxyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 124);

    auto ta1_x_xxyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 125);

    auto ta1_x_xxyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 126);

    auto ta1_x_xxyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 127);

    auto ta1_x_xxyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 128);

    auto ta1_x_xxyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 129);

    auto ta1_x_xxyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 130);

    auto ta1_x_xxyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 131);

    auto ta1_x_xxyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 132);

    auto ta1_x_xxyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 133);

    auto ta1_x_xxyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 134);

    auto ta1_x_xxyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 135);

    auto ta1_x_xxyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 136);

    auto ta1_x_xxyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 137);

    auto ta1_x_xxyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 138);

    auto ta1_x_xxyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 139);

#pragma omp simd aligned(pa_y,                    \
                             pa_z,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_x_xxy_xxxxxy_0,  \
                             ta1_x_xxy_xxxxxy_1,  \
                             ta1_x_xxy_xxxxyy_0,  \
                             ta1_x_xxy_xxxxyy_1,  \
                             ta1_x_xxy_xxxyyy_0,  \
                             ta1_x_xxy_xxxyyy_1,  \
                             ta1_x_xxy_xxyyyy_0,  \
                             ta1_x_xxy_xxyyyy_1,  \
                             ta1_x_xxy_xyyyyy_0,  \
                             ta1_x_xxy_xyyyyy_1,  \
                             ta1_x_xxy_yyyyyy_0,  \
                             ta1_x_xxy_yyyyyy_1,  \
                             ta1_x_xxyz_xxxxxx_0, \
                             ta1_x_xxyz_xxxxxy_0, \
                             ta1_x_xxyz_xxxxxz_0, \
                             ta1_x_xxyz_xxxxyy_0, \
                             ta1_x_xxyz_xxxxyz_0, \
                             ta1_x_xxyz_xxxxzz_0, \
                             ta1_x_xxyz_xxxyyy_0, \
                             ta1_x_xxyz_xxxyyz_0, \
                             ta1_x_xxyz_xxxyzz_0, \
                             ta1_x_xxyz_xxxzzz_0, \
                             ta1_x_xxyz_xxyyyy_0, \
                             ta1_x_xxyz_xxyyyz_0, \
                             ta1_x_xxyz_xxyyzz_0, \
                             ta1_x_xxyz_xxyzzz_0, \
                             ta1_x_xxyz_xxzzzz_0, \
                             ta1_x_xxyz_xyyyyy_0, \
                             ta1_x_xxyz_xyyyyz_0, \
                             ta1_x_xxyz_xyyyzz_0, \
                             ta1_x_xxyz_xyyzzz_0, \
                             ta1_x_xxyz_xyzzzz_0, \
                             ta1_x_xxyz_xzzzzz_0, \
                             ta1_x_xxyz_yyyyyy_0, \
                             ta1_x_xxyz_yyyyyz_0, \
                             ta1_x_xxyz_yyyyzz_0, \
                             ta1_x_xxyz_yyyzzz_0, \
                             ta1_x_xxyz_yyzzzz_0, \
                             ta1_x_xxyz_yzzzzz_0, \
                             ta1_x_xxyz_zzzzzz_0, \
                             ta1_x_xxz_xxxxxx_0,  \
                             ta1_x_xxz_xxxxxx_1,  \
                             ta1_x_xxz_xxxxxz_0,  \
                             ta1_x_xxz_xxxxxz_1,  \
                             ta1_x_xxz_xxxxyz_0,  \
                             ta1_x_xxz_xxxxyz_1,  \
                             ta1_x_xxz_xxxxz_0,   \
                             ta1_x_xxz_xxxxz_1,   \
                             ta1_x_xxz_xxxxzz_0,  \
                             ta1_x_xxz_xxxxzz_1,  \
                             ta1_x_xxz_xxxyyz_0,  \
                             ta1_x_xxz_xxxyyz_1,  \
                             ta1_x_xxz_xxxyz_0,   \
                             ta1_x_xxz_xxxyz_1,   \
                             ta1_x_xxz_xxxyzz_0,  \
                             ta1_x_xxz_xxxyzz_1,  \
                             ta1_x_xxz_xxxzz_0,   \
                             ta1_x_xxz_xxxzz_1,   \
                             ta1_x_xxz_xxxzzz_0,  \
                             ta1_x_xxz_xxxzzz_1,  \
                             ta1_x_xxz_xxyyyz_0,  \
                             ta1_x_xxz_xxyyyz_1,  \
                             ta1_x_xxz_xxyyz_0,   \
                             ta1_x_xxz_xxyyz_1,   \
                             ta1_x_xxz_xxyyzz_0,  \
                             ta1_x_xxz_xxyyzz_1,  \
                             ta1_x_xxz_xxyzz_0,   \
                             ta1_x_xxz_xxyzz_1,   \
                             ta1_x_xxz_xxyzzz_0,  \
                             ta1_x_xxz_xxyzzz_1,  \
                             ta1_x_xxz_xxzzz_0,   \
                             ta1_x_xxz_xxzzz_1,   \
                             ta1_x_xxz_xxzzzz_0,  \
                             ta1_x_xxz_xxzzzz_1,  \
                             ta1_x_xxz_xyyyyz_0,  \
                             ta1_x_xxz_xyyyyz_1,  \
                             ta1_x_xxz_xyyyz_0,   \
                             ta1_x_xxz_xyyyz_1,   \
                             ta1_x_xxz_xyyyzz_0,  \
                             ta1_x_xxz_xyyyzz_1,  \
                             ta1_x_xxz_xyyzz_0,   \
                             ta1_x_xxz_xyyzz_1,   \
                             ta1_x_xxz_xyyzzz_0,  \
                             ta1_x_xxz_xyyzzz_1,  \
                             ta1_x_xxz_xyzzz_0,   \
                             ta1_x_xxz_xyzzz_1,   \
                             ta1_x_xxz_xyzzzz_0,  \
                             ta1_x_xxz_xyzzzz_1,  \
                             ta1_x_xxz_xzzzz_0,   \
                             ta1_x_xxz_xzzzz_1,   \
                             ta1_x_xxz_xzzzzz_0,  \
                             ta1_x_xxz_xzzzzz_1,  \
                             ta1_x_xxz_yyyyyz_0,  \
                             ta1_x_xxz_yyyyyz_1,  \
                             ta1_x_xxz_yyyyz_0,   \
                             ta1_x_xxz_yyyyz_1,   \
                             ta1_x_xxz_yyyyzz_0,  \
                             ta1_x_xxz_yyyyzz_1,  \
                             ta1_x_xxz_yyyzz_0,   \
                             ta1_x_xxz_yyyzz_1,   \
                             ta1_x_xxz_yyyzzz_0,  \
                             ta1_x_xxz_yyyzzz_1,  \
                             ta1_x_xxz_yyzzz_0,   \
                             ta1_x_xxz_yyzzz_1,   \
                             ta1_x_xxz_yyzzzz_0,  \
                             ta1_x_xxz_yyzzzz_1,  \
                             ta1_x_xxz_yzzzz_0,   \
                             ta1_x_xxz_yzzzz_1,   \
                             ta1_x_xxz_yzzzzz_0,  \
                             ta1_x_xxz_yzzzzz_1,  \
                             ta1_x_xxz_zzzzz_0,   \
                             ta1_x_xxz_zzzzz_1,   \
                             ta1_x_xxz_zzzzzz_0,  \
                             ta1_x_xxz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyz_xxxxxx_0[i] = ta1_x_xxz_xxxxxx_0[i] * pa_y[i] - ta1_x_xxz_xxxxxx_1[i] * pc_y[i];

        ta1_x_xxyz_xxxxxy_0[i] = ta1_x_xxy_xxxxxy_0[i] * pa_z[i] - ta1_x_xxy_xxxxxy_1[i] * pc_z[i];

        ta1_x_xxyz_xxxxxz_0[i] = ta1_x_xxz_xxxxxz_0[i] * pa_y[i] - ta1_x_xxz_xxxxxz_1[i] * pc_y[i];

        ta1_x_xxyz_xxxxyy_0[i] = ta1_x_xxy_xxxxyy_0[i] * pa_z[i] - ta1_x_xxy_xxxxyy_1[i] * pc_z[i];

        ta1_x_xxyz_xxxxyz_0[i] =
            ta1_x_xxz_xxxxz_0[i] * fe_0 - ta1_x_xxz_xxxxz_1[i] * fe_0 + ta1_x_xxz_xxxxyz_0[i] * pa_y[i] - ta1_x_xxz_xxxxyz_1[i] * pc_y[i];

        ta1_x_xxyz_xxxxzz_0[i] = ta1_x_xxz_xxxxzz_0[i] * pa_y[i] - ta1_x_xxz_xxxxzz_1[i] * pc_y[i];

        ta1_x_xxyz_xxxyyy_0[i] = ta1_x_xxy_xxxyyy_0[i] * pa_z[i] - ta1_x_xxy_xxxyyy_1[i] * pc_z[i];

        ta1_x_xxyz_xxxyyz_0[i] =
            2.0 * ta1_x_xxz_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxxyz_1[i] * fe_0 + ta1_x_xxz_xxxyyz_0[i] * pa_y[i] - ta1_x_xxz_xxxyyz_1[i] * pc_y[i];

        ta1_x_xxyz_xxxyzz_0[i] =
            ta1_x_xxz_xxxzz_0[i] * fe_0 - ta1_x_xxz_xxxzz_1[i] * fe_0 + ta1_x_xxz_xxxyzz_0[i] * pa_y[i] - ta1_x_xxz_xxxyzz_1[i] * pc_y[i];

        ta1_x_xxyz_xxxzzz_0[i] = ta1_x_xxz_xxxzzz_0[i] * pa_y[i] - ta1_x_xxz_xxxzzz_1[i] * pc_y[i];

        ta1_x_xxyz_xxyyyy_0[i] = ta1_x_xxy_xxyyyy_0[i] * pa_z[i] - ta1_x_xxy_xxyyyy_1[i] * pc_z[i];

        ta1_x_xxyz_xxyyyz_0[i] =
            3.0 * ta1_x_xxz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_xxz_xxyyz_1[i] * fe_0 + ta1_x_xxz_xxyyyz_0[i] * pa_y[i] - ta1_x_xxz_xxyyyz_1[i] * pc_y[i];

        ta1_x_xxyz_xxyyzz_0[i] =
            2.0 * ta1_x_xxz_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxyzz_1[i] * fe_0 + ta1_x_xxz_xxyyzz_0[i] * pa_y[i] - ta1_x_xxz_xxyyzz_1[i] * pc_y[i];

        ta1_x_xxyz_xxyzzz_0[i] =
            ta1_x_xxz_xxzzz_0[i] * fe_0 - ta1_x_xxz_xxzzz_1[i] * fe_0 + ta1_x_xxz_xxyzzz_0[i] * pa_y[i] - ta1_x_xxz_xxyzzz_1[i] * pc_y[i];

        ta1_x_xxyz_xxzzzz_0[i] = ta1_x_xxz_xxzzzz_0[i] * pa_y[i] - ta1_x_xxz_xxzzzz_1[i] * pc_y[i];

        ta1_x_xxyz_xyyyyy_0[i] = ta1_x_xxy_xyyyyy_0[i] * pa_z[i] - ta1_x_xxy_xyyyyy_1[i] * pc_z[i];

        ta1_x_xxyz_xyyyyz_0[i] =
            4.0 * ta1_x_xxz_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_xxz_xyyyz_1[i] * fe_0 + ta1_x_xxz_xyyyyz_0[i] * pa_y[i] - ta1_x_xxz_xyyyyz_1[i] * pc_y[i];

        ta1_x_xxyz_xyyyzz_0[i] =
            3.0 * ta1_x_xxz_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xxz_xyyzz_1[i] * fe_0 + ta1_x_xxz_xyyyzz_0[i] * pa_y[i] - ta1_x_xxz_xyyyzz_1[i] * pc_y[i];

        ta1_x_xxyz_xyyzzz_0[i] =
            2.0 * ta1_x_xxz_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyzzz_1[i] * fe_0 + ta1_x_xxz_xyyzzz_0[i] * pa_y[i] - ta1_x_xxz_xyyzzz_1[i] * pc_y[i];

        ta1_x_xxyz_xyzzzz_0[i] =
            ta1_x_xxz_xzzzz_0[i] * fe_0 - ta1_x_xxz_xzzzz_1[i] * fe_0 + ta1_x_xxz_xyzzzz_0[i] * pa_y[i] - ta1_x_xxz_xyzzzz_1[i] * pc_y[i];

        ta1_x_xxyz_xzzzzz_0[i] = ta1_x_xxz_xzzzzz_0[i] * pa_y[i] - ta1_x_xxz_xzzzzz_1[i] * pc_y[i];

        ta1_x_xxyz_yyyyyy_0[i] = ta1_x_xxy_yyyyyy_0[i] * pa_z[i] - ta1_x_xxy_yyyyyy_1[i] * pc_z[i];

        ta1_x_xxyz_yyyyyz_0[i] =
            5.0 * ta1_x_xxz_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_xxz_yyyyz_1[i] * fe_0 + ta1_x_xxz_yyyyyz_0[i] * pa_y[i] - ta1_x_xxz_yyyyyz_1[i] * pc_y[i];

        ta1_x_xxyz_yyyyzz_0[i] =
            4.0 * ta1_x_xxz_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_xxz_yyyzz_1[i] * fe_0 + ta1_x_xxz_yyyyzz_0[i] * pa_y[i] - ta1_x_xxz_yyyyzz_1[i] * pc_y[i];

        ta1_x_xxyz_yyyzzz_0[i] =
            3.0 * ta1_x_xxz_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_xxz_yyzzz_1[i] * fe_0 + ta1_x_xxz_yyyzzz_0[i] * pa_y[i] - ta1_x_xxz_yyyzzz_1[i] * pc_y[i];

        ta1_x_xxyz_yyzzzz_0[i] =
            2.0 * ta1_x_xxz_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_yzzzz_1[i] * fe_0 + ta1_x_xxz_yyzzzz_0[i] * pa_y[i] - ta1_x_xxz_yyzzzz_1[i] * pc_y[i];

        ta1_x_xxyz_yzzzzz_0[i] =
            ta1_x_xxz_zzzzz_0[i] * fe_0 - ta1_x_xxz_zzzzz_1[i] * fe_0 + ta1_x_xxz_yzzzzz_0[i] * pa_y[i] - ta1_x_xxz_yzzzzz_1[i] * pc_y[i];

        ta1_x_xxyz_zzzzzz_0[i] = ta1_x_xxz_zzzzzz_0[i] * pa_y[i] - ta1_x_xxz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : GI

    auto ta1_x_xxzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 140);

    auto ta1_x_xxzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 141);

    auto ta1_x_xxzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 142);

    auto ta1_x_xxzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 143);

    auto ta1_x_xxzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 144);

    auto ta1_x_xxzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 145);

    auto ta1_x_xxzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 146);

    auto ta1_x_xxzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 147);

    auto ta1_x_xxzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 148);

    auto ta1_x_xxzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 149);

    auto ta1_x_xxzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 150);

    auto ta1_x_xxzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 151);

    auto ta1_x_xxzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 152);

    auto ta1_x_xxzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 153);

    auto ta1_x_xxzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 154);

    auto ta1_x_xxzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 155);

    auto ta1_x_xxzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 156);

    auto ta1_x_xxzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 157);

    auto ta1_x_xxzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 158);

    auto ta1_x_xxzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 159);

    auto ta1_x_xxzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 160);

    auto ta1_x_xxzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 161);

    auto ta1_x_xxzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 162);

    auto ta1_x_xxzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 163);

    auto ta1_x_xxzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 164);

    auto ta1_x_xxzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 165);

    auto ta1_x_xxzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 166);

    auto ta1_x_xxzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 167);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_x_xx_xxxxxx_0,   \
                             ta1_x_xx_xxxxxx_1,   \
                             ta1_x_xx_xxxxxy_0,   \
                             ta1_x_xx_xxxxxy_1,   \
                             ta1_x_xx_xxxxxz_0,   \
                             ta1_x_xx_xxxxxz_1,   \
                             ta1_x_xx_xxxxyy_0,   \
                             ta1_x_xx_xxxxyy_1,   \
                             ta1_x_xx_xxxxyz_0,   \
                             ta1_x_xx_xxxxyz_1,   \
                             ta1_x_xx_xxxxzz_0,   \
                             ta1_x_xx_xxxxzz_1,   \
                             ta1_x_xx_xxxyyy_0,   \
                             ta1_x_xx_xxxyyy_1,   \
                             ta1_x_xx_xxxyyz_0,   \
                             ta1_x_xx_xxxyyz_1,   \
                             ta1_x_xx_xxxyzz_0,   \
                             ta1_x_xx_xxxyzz_1,   \
                             ta1_x_xx_xxxzzz_0,   \
                             ta1_x_xx_xxxzzz_1,   \
                             ta1_x_xx_xxyyyy_0,   \
                             ta1_x_xx_xxyyyy_1,   \
                             ta1_x_xx_xxyyyz_0,   \
                             ta1_x_xx_xxyyyz_1,   \
                             ta1_x_xx_xxyyzz_0,   \
                             ta1_x_xx_xxyyzz_1,   \
                             ta1_x_xx_xxyzzz_0,   \
                             ta1_x_xx_xxyzzz_1,   \
                             ta1_x_xx_xxzzzz_0,   \
                             ta1_x_xx_xxzzzz_1,   \
                             ta1_x_xx_xyyyyy_0,   \
                             ta1_x_xx_xyyyyy_1,   \
                             ta1_x_xx_xyyyyz_0,   \
                             ta1_x_xx_xyyyyz_1,   \
                             ta1_x_xx_xyyyzz_0,   \
                             ta1_x_xx_xyyyzz_1,   \
                             ta1_x_xx_xyyzzz_0,   \
                             ta1_x_xx_xyyzzz_1,   \
                             ta1_x_xx_xyzzzz_0,   \
                             ta1_x_xx_xyzzzz_1,   \
                             ta1_x_xx_xzzzzz_0,   \
                             ta1_x_xx_xzzzzz_1,   \
                             ta1_x_xx_yyyyyy_0,   \
                             ta1_x_xx_yyyyyy_1,   \
                             ta1_x_xxz_xxxxx_0,   \
                             ta1_x_xxz_xxxxx_1,   \
                             ta1_x_xxz_xxxxxx_0,  \
                             ta1_x_xxz_xxxxxx_1,  \
                             ta1_x_xxz_xxxxxy_0,  \
                             ta1_x_xxz_xxxxxy_1,  \
                             ta1_x_xxz_xxxxxz_0,  \
                             ta1_x_xxz_xxxxxz_1,  \
                             ta1_x_xxz_xxxxy_0,   \
                             ta1_x_xxz_xxxxy_1,   \
                             ta1_x_xxz_xxxxyy_0,  \
                             ta1_x_xxz_xxxxyy_1,  \
                             ta1_x_xxz_xxxxyz_0,  \
                             ta1_x_xxz_xxxxyz_1,  \
                             ta1_x_xxz_xxxxz_0,   \
                             ta1_x_xxz_xxxxz_1,   \
                             ta1_x_xxz_xxxxzz_0,  \
                             ta1_x_xxz_xxxxzz_1,  \
                             ta1_x_xxz_xxxyy_0,   \
                             ta1_x_xxz_xxxyy_1,   \
                             ta1_x_xxz_xxxyyy_0,  \
                             ta1_x_xxz_xxxyyy_1,  \
                             ta1_x_xxz_xxxyyz_0,  \
                             ta1_x_xxz_xxxyyz_1,  \
                             ta1_x_xxz_xxxyz_0,   \
                             ta1_x_xxz_xxxyz_1,   \
                             ta1_x_xxz_xxxyzz_0,  \
                             ta1_x_xxz_xxxyzz_1,  \
                             ta1_x_xxz_xxxzz_0,   \
                             ta1_x_xxz_xxxzz_1,   \
                             ta1_x_xxz_xxxzzz_0,  \
                             ta1_x_xxz_xxxzzz_1,  \
                             ta1_x_xxz_xxyyy_0,   \
                             ta1_x_xxz_xxyyy_1,   \
                             ta1_x_xxz_xxyyyy_0,  \
                             ta1_x_xxz_xxyyyy_1,  \
                             ta1_x_xxz_xxyyyz_0,  \
                             ta1_x_xxz_xxyyyz_1,  \
                             ta1_x_xxz_xxyyz_0,   \
                             ta1_x_xxz_xxyyz_1,   \
                             ta1_x_xxz_xxyyzz_0,  \
                             ta1_x_xxz_xxyyzz_1,  \
                             ta1_x_xxz_xxyzz_0,   \
                             ta1_x_xxz_xxyzz_1,   \
                             ta1_x_xxz_xxyzzz_0,  \
                             ta1_x_xxz_xxyzzz_1,  \
                             ta1_x_xxz_xxzzz_0,   \
                             ta1_x_xxz_xxzzz_1,   \
                             ta1_x_xxz_xxzzzz_0,  \
                             ta1_x_xxz_xxzzzz_1,  \
                             ta1_x_xxz_xyyyy_0,   \
                             ta1_x_xxz_xyyyy_1,   \
                             ta1_x_xxz_xyyyyy_0,  \
                             ta1_x_xxz_xyyyyy_1,  \
                             ta1_x_xxz_xyyyyz_0,  \
                             ta1_x_xxz_xyyyyz_1,  \
                             ta1_x_xxz_xyyyz_0,   \
                             ta1_x_xxz_xyyyz_1,   \
                             ta1_x_xxz_xyyyzz_0,  \
                             ta1_x_xxz_xyyyzz_1,  \
                             ta1_x_xxz_xyyzz_0,   \
                             ta1_x_xxz_xyyzz_1,   \
                             ta1_x_xxz_xyyzzz_0,  \
                             ta1_x_xxz_xyyzzz_1,  \
                             ta1_x_xxz_xyzzz_0,   \
                             ta1_x_xxz_xyzzz_1,   \
                             ta1_x_xxz_xyzzzz_0,  \
                             ta1_x_xxz_xyzzzz_1,  \
                             ta1_x_xxz_xzzzz_0,   \
                             ta1_x_xxz_xzzzz_1,   \
                             ta1_x_xxz_xzzzzz_0,  \
                             ta1_x_xxz_xzzzzz_1,  \
                             ta1_x_xxz_yyyyyy_0,  \
                             ta1_x_xxz_yyyyyy_1,  \
                             ta1_x_xxzz_xxxxxx_0, \
                             ta1_x_xxzz_xxxxxy_0, \
                             ta1_x_xxzz_xxxxxz_0, \
                             ta1_x_xxzz_xxxxyy_0, \
                             ta1_x_xxzz_xxxxyz_0, \
                             ta1_x_xxzz_xxxxzz_0, \
                             ta1_x_xxzz_xxxyyy_0, \
                             ta1_x_xxzz_xxxyyz_0, \
                             ta1_x_xxzz_xxxyzz_0, \
                             ta1_x_xxzz_xxxzzz_0, \
                             ta1_x_xxzz_xxyyyy_0, \
                             ta1_x_xxzz_xxyyyz_0, \
                             ta1_x_xxzz_xxyyzz_0, \
                             ta1_x_xxzz_xxyzzz_0, \
                             ta1_x_xxzz_xxzzzz_0, \
                             ta1_x_xxzz_xyyyyy_0, \
                             ta1_x_xxzz_xyyyyz_0, \
                             ta1_x_xxzz_xyyyzz_0, \
                             ta1_x_xxzz_xyyzzz_0, \
                             ta1_x_xxzz_xyzzzz_0, \
                             ta1_x_xxzz_xzzzzz_0, \
                             ta1_x_xxzz_yyyyyy_0, \
                             ta1_x_xxzz_yyyyyz_0, \
                             ta1_x_xxzz_yyyyzz_0, \
                             ta1_x_xxzz_yyyzzz_0, \
                             ta1_x_xxzz_yyzzzz_0, \
                             ta1_x_xxzz_yzzzzz_0, \
                             ta1_x_xxzz_zzzzzz_0, \
                             ta1_x_xzz_yyyyyz_0,  \
                             ta1_x_xzz_yyyyyz_1,  \
                             ta1_x_xzz_yyyyzz_0,  \
                             ta1_x_xzz_yyyyzz_1,  \
                             ta1_x_xzz_yyyzzz_0,  \
                             ta1_x_xzz_yyyzzz_1,  \
                             ta1_x_xzz_yyzzzz_0,  \
                             ta1_x_xzz_yyzzzz_1,  \
                             ta1_x_xzz_yzzzzz_0,  \
                             ta1_x_xzz_yzzzzz_1,  \
                             ta1_x_xzz_zzzzzz_0,  \
                             ta1_x_xzz_zzzzzz_1,  \
                             ta1_x_zz_yyyyyz_0,   \
                             ta1_x_zz_yyyyyz_1,   \
                             ta1_x_zz_yyyyzz_0,   \
                             ta1_x_zz_yyyyzz_1,   \
                             ta1_x_zz_yyyzzz_0,   \
                             ta1_x_zz_yyyzzz_1,   \
                             ta1_x_zz_yyzzzz_0,   \
                             ta1_x_zz_yyzzzz_1,   \
                             ta1_x_zz_yzzzzz_0,   \
                             ta1_x_zz_yzzzzz_1,   \
                             ta1_x_zz_zzzzzz_0,   \
                             ta1_x_zz_zzzzzz_1,   \
                             ta_xzz_yyyyyz_1,     \
                             ta_xzz_yyyyzz_1,     \
                             ta_xzz_yyyzzz_1,     \
                             ta_xzz_yyzzzz_1,     \
                             ta_xzz_yzzzzz_1,     \
                             ta_xzz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzz_xxxxxx_0[i] =
            ta1_x_xx_xxxxxx_0[i] * fe_0 - ta1_x_xx_xxxxxx_1[i] * fe_0 + ta1_x_xxz_xxxxxx_0[i] * pa_z[i] - ta1_x_xxz_xxxxxx_1[i] * pc_z[i];

        ta1_x_xxzz_xxxxxy_0[i] =
            ta1_x_xx_xxxxxy_0[i] * fe_0 - ta1_x_xx_xxxxxy_1[i] * fe_0 + ta1_x_xxz_xxxxxy_0[i] * pa_z[i] - ta1_x_xxz_xxxxxy_1[i] * pc_z[i];

        ta1_x_xxzz_xxxxxz_0[i] = ta1_x_xx_xxxxxz_0[i] * fe_0 - ta1_x_xx_xxxxxz_1[i] * fe_0 + ta1_x_xxz_xxxxx_0[i] * fe_0 -
                                 ta1_x_xxz_xxxxx_1[i] * fe_0 + ta1_x_xxz_xxxxxz_0[i] * pa_z[i] - ta1_x_xxz_xxxxxz_1[i] * pc_z[i];

        ta1_x_xxzz_xxxxyy_0[i] =
            ta1_x_xx_xxxxyy_0[i] * fe_0 - ta1_x_xx_xxxxyy_1[i] * fe_0 + ta1_x_xxz_xxxxyy_0[i] * pa_z[i] - ta1_x_xxz_xxxxyy_1[i] * pc_z[i];

        ta1_x_xxzz_xxxxyz_0[i] = ta1_x_xx_xxxxyz_0[i] * fe_0 - ta1_x_xx_xxxxyz_1[i] * fe_0 + ta1_x_xxz_xxxxy_0[i] * fe_0 -
                                 ta1_x_xxz_xxxxy_1[i] * fe_0 + ta1_x_xxz_xxxxyz_0[i] * pa_z[i] - ta1_x_xxz_xxxxyz_1[i] * pc_z[i];

        ta1_x_xxzz_xxxxzz_0[i] = ta1_x_xx_xxxxzz_0[i] * fe_0 - ta1_x_xx_xxxxzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xxxxz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxz_xxxxz_1[i] * fe_0 + ta1_x_xxz_xxxxzz_0[i] * pa_z[i] - ta1_x_xxz_xxxxzz_1[i] * pc_z[i];

        ta1_x_xxzz_xxxyyy_0[i] =
            ta1_x_xx_xxxyyy_0[i] * fe_0 - ta1_x_xx_xxxyyy_1[i] * fe_0 + ta1_x_xxz_xxxyyy_0[i] * pa_z[i] - ta1_x_xxz_xxxyyy_1[i] * pc_z[i];

        ta1_x_xxzz_xxxyyz_0[i] = ta1_x_xx_xxxyyz_0[i] * fe_0 - ta1_x_xx_xxxyyz_1[i] * fe_0 + ta1_x_xxz_xxxyy_0[i] * fe_0 -
                                 ta1_x_xxz_xxxyy_1[i] * fe_0 + ta1_x_xxz_xxxyyz_0[i] * pa_z[i] - ta1_x_xxz_xxxyyz_1[i] * pc_z[i];

        ta1_x_xxzz_xxxyzz_0[i] = ta1_x_xx_xxxyzz_0[i] * fe_0 - ta1_x_xx_xxxyzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxz_xxxyz_1[i] * fe_0 + ta1_x_xxz_xxxyzz_0[i] * pa_z[i] - ta1_x_xxz_xxxyzz_1[i] * pc_z[i];

        ta1_x_xxzz_xxxzzz_0[i] = ta1_x_xx_xxxzzz_0[i] * fe_0 - ta1_x_xx_xxxzzz_1[i] * fe_0 + 3.0 * ta1_x_xxz_xxxzz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxz_xxxzz_1[i] * fe_0 + ta1_x_xxz_xxxzzz_0[i] * pa_z[i] - ta1_x_xxz_xxxzzz_1[i] * pc_z[i];

        ta1_x_xxzz_xxyyyy_0[i] =
            ta1_x_xx_xxyyyy_0[i] * fe_0 - ta1_x_xx_xxyyyy_1[i] * fe_0 + ta1_x_xxz_xxyyyy_0[i] * pa_z[i] - ta1_x_xxz_xxyyyy_1[i] * pc_z[i];

        ta1_x_xxzz_xxyyyz_0[i] = ta1_x_xx_xxyyyz_0[i] * fe_0 - ta1_x_xx_xxyyyz_1[i] * fe_0 + ta1_x_xxz_xxyyy_0[i] * fe_0 -
                                 ta1_x_xxz_xxyyy_1[i] * fe_0 + ta1_x_xxz_xxyyyz_0[i] * pa_z[i] - ta1_x_xxz_xxyyyz_1[i] * pc_z[i];

        ta1_x_xxzz_xxyyzz_0[i] = ta1_x_xx_xxyyzz_0[i] * fe_0 - ta1_x_xx_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xxyyz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxz_xxyyz_1[i] * fe_0 + ta1_x_xxz_xxyyzz_0[i] * pa_z[i] - ta1_x_xxz_xxyyzz_1[i] * pc_z[i];

        ta1_x_xxzz_xxyzzz_0[i] = ta1_x_xx_xxyzzz_0[i] * fe_0 - ta1_x_xx_xxyzzz_1[i] * fe_0 + 3.0 * ta1_x_xxz_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxz_xxyzz_1[i] * fe_0 + ta1_x_xxz_xxyzzz_0[i] * pa_z[i] - ta1_x_xxz_xxyzzz_1[i] * pc_z[i];

        ta1_x_xxzz_xxzzzz_0[i] = ta1_x_xx_xxzzzz_0[i] * fe_0 - ta1_x_xx_xxzzzz_1[i] * fe_0 + 4.0 * ta1_x_xxz_xxzzz_0[i] * fe_0 -
                                 4.0 * ta1_x_xxz_xxzzz_1[i] * fe_0 + ta1_x_xxz_xxzzzz_0[i] * pa_z[i] - ta1_x_xxz_xxzzzz_1[i] * pc_z[i];

        ta1_x_xxzz_xyyyyy_0[i] =
            ta1_x_xx_xyyyyy_0[i] * fe_0 - ta1_x_xx_xyyyyy_1[i] * fe_0 + ta1_x_xxz_xyyyyy_0[i] * pa_z[i] - ta1_x_xxz_xyyyyy_1[i] * pc_z[i];

        ta1_x_xxzz_xyyyyz_0[i] = ta1_x_xx_xyyyyz_0[i] * fe_0 - ta1_x_xx_xyyyyz_1[i] * fe_0 + ta1_x_xxz_xyyyy_0[i] * fe_0 -
                                 ta1_x_xxz_xyyyy_1[i] * fe_0 + ta1_x_xxz_xyyyyz_0[i] * pa_z[i] - ta1_x_xxz_xyyyyz_1[i] * pc_z[i];

        ta1_x_xxzz_xyyyzz_0[i] = ta1_x_xx_xyyyzz_0[i] * fe_0 - ta1_x_xx_xyyyzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_x_xxz_xyyyz_1[i] * fe_0 + ta1_x_xxz_xyyyzz_0[i] * pa_z[i] - ta1_x_xxz_xyyyzz_1[i] * pc_z[i];

        ta1_x_xxzz_xyyzzz_0[i] = ta1_x_xx_xyyzzz_0[i] * fe_0 - ta1_x_xx_xyyzzz_1[i] * fe_0 + 3.0 * ta1_x_xxz_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_xxz_xyyzz_1[i] * fe_0 + ta1_x_xxz_xyyzzz_0[i] * pa_z[i] - ta1_x_xxz_xyyzzz_1[i] * pc_z[i];

        ta1_x_xxzz_xyzzzz_0[i] = ta1_x_xx_xyzzzz_0[i] * fe_0 - ta1_x_xx_xyzzzz_1[i] * fe_0 + 4.0 * ta1_x_xxz_xyzzz_0[i] * fe_0 -
                                 4.0 * ta1_x_xxz_xyzzz_1[i] * fe_0 + ta1_x_xxz_xyzzzz_0[i] * pa_z[i] - ta1_x_xxz_xyzzzz_1[i] * pc_z[i];

        ta1_x_xxzz_xzzzzz_0[i] = ta1_x_xx_xzzzzz_0[i] * fe_0 - ta1_x_xx_xzzzzz_1[i] * fe_0 + 5.0 * ta1_x_xxz_xzzzz_0[i] * fe_0 -
                                 5.0 * ta1_x_xxz_xzzzz_1[i] * fe_0 + ta1_x_xxz_xzzzzz_0[i] * pa_z[i] - ta1_x_xxz_xzzzzz_1[i] * pc_z[i];

        ta1_x_xxzz_yyyyyy_0[i] =
            ta1_x_xx_yyyyyy_0[i] * fe_0 - ta1_x_xx_yyyyyy_1[i] * fe_0 + ta1_x_xxz_yyyyyy_0[i] * pa_z[i] - ta1_x_xxz_yyyyyy_1[i] * pc_z[i];

        ta1_x_xxzz_yyyyyz_0[i] = ta1_x_zz_yyyyyz_0[i] * fe_0 - ta1_x_zz_yyyyyz_1[i] * fe_0 + ta_xzz_yyyyyz_1[i] + ta1_x_xzz_yyyyyz_0[i] * pa_x[i] -
                                 ta1_x_xzz_yyyyyz_1[i] * pc_x[i];

        ta1_x_xxzz_yyyyzz_0[i] = ta1_x_zz_yyyyzz_0[i] * fe_0 - ta1_x_zz_yyyyzz_1[i] * fe_0 + ta_xzz_yyyyzz_1[i] + ta1_x_xzz_yyyyzz_0[i] * pa_x[i] -
                                 ta1_x_xzz_yyyyzz_1[i] * pc_x[i];

        ta1_x_xxzz_yyyzzz_0[i] = ta1_x_zz_yyyzzz_0[i] * fe_0 - ta1_x_zz_yyyzzz_1[i] * fe_0 + ta_xzz_yyyzzz_1[i] + ta1_x_xzz_yyyzzz_0[i] * pa_x[i] -
                                 ta1_x_xzz_yyyzzz_1[i] * pc_x[i];

        ta1_x_xxzz_yyzzzz_0[i] = ta1_x_zz_yyzzzz_0[i] * fe_0 - ta1_x_zz_yyzzzz_1[i] * fe_0 + ta_xzz_yyzzzz_1[i] + ta1_x_xzz_yyzzzz_0[i] * pa_x[i] -
                                 ta1_x_xzz_yyzzzz_1[i] * pc_x[i];

        ta1_x_xxzz_yzzzzz_0[i] = ta1_x_zz_yzzzzz_0[i] * fe_0 - ta1_x_zz_yzzzzz_1[i] * fe_0 + ta_xzz_yzzzzz_1[i] + ta1_x_xzz_yzzzzz_0[i] * pa_x[i] -
                                 ta1_x_xzz_yzzzzz_1[i] * pc_x[i];

        ta1_x_xxzz_zzzzzz_0[i] = ta1_x_zz_zzzzzz_0[i] * fe_0 - ta1_x_zz_zzzzzz_1[i] * fe_0 + ta_xzz_zzzzzz_1[i] + ta1_x_xzz_zzzzzz_0[i] * pa_x[i] -
                                 ta1_x_xzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 168-196 components of targeted buffer : GI

    auto ta1_x_xyyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 168);

    auto ta1_x_xyyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 169);

    auto ta1_x_xyyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 170);

    auto ta1_x_xyyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 171);

    auto ta1_x_xyyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 172);

    auto ta1_x_xyyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 173);

    auto ta1_x_xyyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 174);

    auto ta1_x_xyyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 175);

    auto ta1_x_xyyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 176);

    auto ta1_x_xyyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 177);

    auto ta1_x_xyyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 178);

    auto ta1_x_xyyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 179);

    auto ta1_x_xyyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 180);

    auto ta1_x_xyyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 181);

    auto ta1_x_xyyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 182);

    auto ta1_x_xyyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 183);

    auto ta1_x_xyyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 184);

    auto ta1_x_xyyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 185);

    auto ta1_x_xyyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 186);

    auto ta1_x_xyyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 187);

    auto ta1_x_xyyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 188);

    auto ta1_x_xyyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 189);

    auto ta1_x_xyyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 190);

    auto ta1_x_xyyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 191);

    auto ta1_x_xyyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 192);

    auto ta1_x_xyyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 193);

    auto ta1_x_xyyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 194);

    auto ta1_x_xyyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 195);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_x_xy_xxxxxx_0,   \
                             ta1_x_xy_xxxxxx_1,   \
                             ta1_x_xy_xxxxxz_0,   \
                             ta1_x_xy_xxxxxz_1,   \
                             ta1_x_xy_xxxxzz_0,   \
                             ta1_x_xy_xxxxzz_1,   \
                             ta1_x_xy_xxxzzz_0,   \
                             ta1_x_xy_xxxzzz_1,   \
                             ta1_x_xy_xxzzzz_0,   \
                             ta1_x_xy_xxzzzz_1,   \
                             ta1_x_xy_xzzzzz_0,   \
                             ta1_x_xy_xzzzzz_1,   \
                             ta1_x_xyy_xxxxxx_0,  \
                             ta1_x_xyy_xxxxxx_1,  \
                             ta1_x_xyy_xxxxxz_0,  \
                             ta1_x_xyy_xxxxxz_1,  \
                             ta1_x_xyy_xxxxzz_0,  \
                             ta1_x_xyy_xxxxzz_1,  \
                             ta1_x_xyy_xxxzzz_0,  \
                             ta1_x_xyy_xxxzzz_1,  \
                             ta1_x_xyy_xxzzzz_0,  \
                             ta1_x_xyy_xxzzzz_1,  \
                             ta1_x_xyy_xzzzzz_0,  \
                             ta1_x_xyy_xzzzzz_1,  \
                             ta1_x_xyyy_xxxxxx_0, \
                             ta1_x_xyyy_xxxxxy_0, \
                             ta1_x_xyyy_xxxxxz_0, \
                             ta1_x_xyyy_xxxxyy_0, \
                             ta1_x_xyyy_xxxxyz_0, \
                             ta1_x_xyyy_xxxxzz_0, \
                             ta1_x_xyyy_xxxyyy_0, \
                             ta1_x_xyyy_xxxyyz_0, \
                             ta1_x_xyyy_xxxyzz_0, \
                             ta1_x_xyyy_xxxzzz_0, \
                             ta1_x_xyyy_xxyyyy_0, \
                             ta1_x_xyyy_xxyyyz_0, \
                             ta1_x_xyyy_xxyyzz_0, \
                             ta1_x_xyyy_xxyzzz_0, \
                             ta1_x_xyyy_xxzzzz_0, \
                             ta1_x_xyyy_xyyyyy_0, \
                             ta1_x_xyyy_xyyyyz_0, \
                             ta1_x_xyyy_xyyyzz_0, \
                             ta1_x_xyyy_xyyzzz_0, \
                             ta1_x_xyyy_xyzzzz_0, \
                             ta1_x_xyyy_xzzzzz_0, \
                             ta1_x_xyyy_yyyyyy_0, \
                             ta1_x_xyyy_yyyyyz_0, \
                             ta1_x_xyyy_yyyyzz_0, \
                             ta1_x_xyyy_yyyzzz_0, \
                             ta1_x_xyyy_yyzzzz_0, \
                             ta1_x_xyyy_yzzzzz_0, \
                             ta1_x_xyyy_zzzzzz_0, \
                             ta1_x_yyy_xxxxxy_0,  \
                             ta1_x_yyy_xxxxxy_1,  \
                             ta1_x_yyy_xxxxy_0,   \
                             ta1_x_yyy_xxxxy_1,   \
                             ta1_x_yyy_xxxxyy_0,  \
                             ta1_x_yyy_xxxxyy_1,  \
                             ta1_x_yyy_xxxxyz_0,  \
                             ta1_x_yyy_xxxxyz_1,  \
                             ta1_x_yyy_xxxyy_0,   \
                             ta1_x_yyy_xxxyy_1,   \
                             ta1_x_yyy_xxxyyy_0,  \
                             ta1_x_yyy_xxxyyy_1,  \
                             ta1_x_yyy_xxxyyz_0,  \
                             ta1_x_yyy_xxxyyz_1,  \
                             ta1_x_yyy_xxxyz_0,   \
                             ta1_x_yyy_xxxyz_1,   \
                             ta1_x_yyy_xxxyzz_0,  \
                             ta1_x_yyy_xxxyzz_1,  \
                             ta1_x_yyy_xxyyy_0,   \
                             ta1_x_yyy_xxyyy_1,   \
                             ta1_x_yyy_xxyyyy_0,  \
                             ta1_x_yyy_xxyyyy_1,  \
                             ta1_x_yyy_xxyyyz_0,  \
                             ta1_x_yyy_xxyyyz_1,  \
                             ta1_x_yyy_xxyyz_0,   \
                             ta1_x_yyy_xxyyz_1,   \
                             ta1_x_yyy_xxyyzz_0,  \
                             ta1_x_yyy_xxyyzz_1,  \
                             ta1_x_yyy_xxyzz_0,   \
                             ta1_x_yyy_xxyzz_1,   \
                             ta1_x_yyy_xxyzzz_0,  \
                             ta1_x_yyy_xxyzzz_1,  \
                             ta1_x_yyy_xyyyy_0,   \
                             ta1_x_yyy_xyyyy_1,   \
                             ta1_x_yyy_xyyyyy_0,  \
                             ta1_x_yyy_xyyyyy_1,  \
                             ta1_x_yyy_xyyyyz_0,  \
                             ta1_x_yyy_xyyyyz_1,  \
                             ta1_x_yyy_xyyyz_0,   \
                             ta1_x_yyy_xyyyz_1,   \
                             ta1_x_yyy_xyyyzz_0,  \
                             ta1_x_yyy_xyyyzz_1,  \
                             ta1_x_yyy_xyyzz_0,   \
                             ta1_x_yyy_xyyzz_1,   \
                             ta1_x_yyy_xyyzzz_0,  \
                             ta1_x_yyy_xyyzzz_1,  \
                             ta1_x_yyy_xyzzz_0,   \
                             ta1_x_yyy_xyzzz_1,   \
                             ta1_x_yyy_xyzzzz_0,  \
                             ta1_x_yyy_xyzzzz_1,  \
                             ta1_x_yyy_yyyyy_0,   \
                             ta1_x_yyy_yyyyy_1,   \
                             ta1_x_yyy_yyyyyy_0,  \
                             ta1_x_yyy_yyyyyy_1,  \
                             ta1_x_yyy_yyyyyz_0,  \
                             ta1_x_yyy_yyyyyz_1,  \
                             ta1_x_yyy_yyyyz_0,   \
                             ta1_x_yyy_yyyyz_1,   \
                             ta1_x_yyy_yyyyzz_0,  \
                             ta1_x_yyy_yyyyzz_1,  \
                             ta1_x_yyy_yyyzz_0,   \
                             ta1_x_yyy_yyyzz_1,   \
                             ta1_x_yyy_yyyzzz_0,  \
                             ta1_x_yyy_yyyzzz_1,  \
                             ta1_x_yyy_yyzzz_0,   \
                             ta1_x_yyy_yyzzz_1,   \
                             ta1_x_yyy_yyzzzz_0,  \
                             ta1_x_yyy_yyzzzz_1,  \
                             ta1_x_yyy_yzzzz_0,   \
                             ta1_x_yyy_yzzzz_1,   \
                             ta1_x_yyy_yzzzzz_0,  \
                             ta1_x_yyy_yzzzzz_1,  \
                             ta1_x_yyy_zzzzzz_0,  \
                             ta1_x_yyy_zzzzzz_1,  \
                             ta_yyy_xxxxxy_1,     \
                             ta_yyy_xxxxyy_1,     \
                             ta_yyy_xxxxyz_1,     \
                             ta_yyy_xxxyyy_1,     \
                             ta_yyy_xxxyyz_1,     \
                             ta_yyy_xxxyzz_1,     \
                             ta_yyy_xxyyyy_1,     \
                             ta_yyy_xxyyyz_1,     \
                             ta_yyy_xxyyzz_1,     \
                             ta_yyy_xxyzzz_1,     \
                             ta_yyy_xyyyyy_1,     \
                             ta_yyy_xyyyyz_1,     \
                             ta_yyy_xyyyzz_1,     \
                             ta_yyy_xyyzzz_1,     \
                             ta_yyy_xyzzzz_1,     \
                             ta_yyy_yyyyyy_1,     \
                             ta_yyy_yyyyyz_1,     \
                             ta_yyy_yyyyzz_1,     \
                             ta_yyy_yyyzzz_1,     \
                             ta_yyy_yyzzzz_1,     \
                             ta_yyy_yzzzzz_1,     \
                             ta_yyy_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyy_xxxxxx_0[i] =
            2.0 * ta1_x_xy_xxxxxx_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxxxx_1[i] * fe_0 + ta1_x_xyy_xxxxxx_0[i] * pa_y[i] - ta1_x_xyy_xxxxxx_1[i] * pc_y[i];

        ta1_x_xyyy_xxxxxy_0[i] = 5.0 * ta1_x_yyy_xxxxy_0[i] * fe_0 - 5.0 * ta1_x_yyy_xxxxy_1[i] * fe_0 + ta_yyy_xxxxxy_1[i] +
                                 ta1_x_yyy_xxxxxy_0[i] * pa_x[i] - ta1_x_yyy_xxxxxy_1[i] * pc_x[i];

        ta1_x_xyyy_xxxxxz_0[i] =
            2.0 * ta1_x_xy_xxxxxz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxxxz_1[i] * fe_0 + ta1_x_xyy_xxxxxz_0[i] * pa_y[i] - ta1_x_xyy_xxxxxz_1[i] * pc_y[i];

        ta1_x_xyyy_xxxxyy_0[i] = 4.0 * ta1_x_yyy_xxxyy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxxyy_1[i] * fe_0 + ta_yyy_xxxxyy_1[i] +
                                 ta1_x_yyy_xxxxyy_0[i] * pa_x[i] - ta1_x_yyy_xxxxyy_1[i] * pc_x[i];

        ta1_x_xyyy_xxxxyz_0[i] = 4.0 * ta1_x_yyy_xxxyz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxxyz_1[i] * fe_0 + ta_yyy_xxxxyz_1[i] +
                                 ta1_x_yyy_xxxxyz_0[i] * pa_x[i] - ta1_x_yyy_xxxxyz_1[i] * pc_x[i];

        ta1_x_xyyy_xxxxzz_0[i] =
            2.0 * ta1_x_xy_xxxxzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxxzz_1[i] * fe_0 + ta1_x_xyy_xxxxzz_0[i] * pa_y[i] - ta1_x_xyy_xxxxzz_1[i] * pc_y[i];

        ta1_x_xyyy_xxxyyy_0[i] = 3.0 * ta1_x_yyy_xxyyy_0[i] * fe_0 - 3.0 * ta1_x_yyy_xxyyy_1[i] * fe_0 + ta_yyy_xxxyyy_1[i] +
                                 ta1_x_yyy_xxxyyy_0[i] * pa_x[i] - ta1_x_yyy_xxxyyy_1[i] * pc_x[i];

        ta1_x_xyyy_xxxyyz_0[i] = 3.0 * ta1_x_yyy_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_yyy_xxyyz_1[i] * fe_0 + ta_yyy_xxxyyz_1[i] +
                                 ta1_x_yyy_xxxyyz_0[i] * pa_x[i] - ta1_x_yyy_xxxyyz_1[i] * pc_x[i];

        ta1_x_xyyy_xxxyzz_0[i] = 3.0 * ta1_x_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_yyy_xxyzz_1[i] * fe_0 + ta_yyy_xxxyzz_1[i] +
                                 ta1_x_yyy_xxxyzz_0[i] * pa_x[i] - ta1_x_yyy_xxxyzz_1[i] * pc_x[i];

        ta1_x_xyyy_xxxzzz_0[i] =
            2.0 * ta1_x_xy_xxxzzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxzzz_1[i] * fe_0 + ta1_x_xyy_xxxzzz_0[i] * pa_y[i] - ta1_x_xyy_xxxzzz_1[i] * pc_y[i];

        ta1_x_xyyy_xxyyyy_0[i] = 2.0 * ta1_x_yyy_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyyyy_1[i] * fe_0 + ta_yyy_xxyyyy_1[i] +
                                 ta1_x_yyy_xxyyyy_0[i] * pa_x[i] - ta1_x_yyy_xxyyyy_1[i] * pc_x[i];

        ta1_x_xyyy_xxyyyz_0[i] = 2.0 * ta1_x_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyyyz_1[i] * fe_0 + ta_yyy_xxyyyz_1[i] +
                                 ta1_x_yyy_xxyyyz_0[i] * pa_x[i] - ta1_x_yyy_xxyyyz_1[i] * pc_x[i];

        ta1_x_xyyy_xxyyzz_0[i] = 2.0 * ta1_x_yyy_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyyzz_1[i] * fe_0 + ta_yyy_xxyyzz_1[i] +
                                 ta1_x_yyy_xxyyzz_0[i] * pa_x[i] - ta1_x_yyy_xxyyzz_1[i] * pc_x[i];

        ta1_x_xyyy_xxyzzz_0[i] = 2.0 * ta1_x_yyy_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyzzz_1[i] * fe_0 + ta_yyy_xxyzzz_1[i] +
                                 ta1_x_yyy_xxyzzz_0[i] * pa_x[i] - ta1_x_yyy_xxyzzz_1[i] * pc_x[i];

        ta1_x_xyyy_xxzzzz_0[i] =
            2.0 * ta1_x_xy_xxzzzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxzzzz_1[i] * fe_0 + ta1_x_xyy_xxzzzz_0[i] * pa_y[i] - ta1_x_xyy_xxzzzz_1[i] * pc_y[i];

        ta1_x_xyyy_xyyyyy_0[i] = ta1_x_yyy_yyyyy_0[i] * fe_0 - ta1_x_yyy_yyyyy_1[i] * fe_0 + ta_yyy_xyyyyy_1[i] + ta1_x_yyy_xyyyyy_0[i] * pa_x[i] -
                                 ta1_x_yyy_xyyyyy_1[i] * pc_x[i];

        ta1_x_xyyy_xyyyyz_0[i] = ta1_x_yyy_yyyyz_0[i] * fe_0 - ta1_x_yyy_yyyyz_1[i] * fe_0 + ta_yyy_xyyyyz_1[i] + ta1_x_yyy_xyyyyz_0[i] * pa_x[i] -
                                 ta1_x_yyy_xyyyyz_1[i] * pc_x[i];

        ta1_x_xyyy_xyyyzz_0[i] = ta1_x_yyy_yyyzz_0[i] * fe_0 - ta1_x_yyy_yyyzz_1[i] * fe_0 + ta_yyy_xyyyzz_1[i] + ta1_x_yyy_xyyyzz_0[i] * pa_x[i] -
                                 ta1_x_yyy_xyyyzz_1[i] * pc_x[i];

        ta1_x_xyyy_xyyzzz_0[i] = ta1_x_yyy_yyzzz_0[i] * fe_0 - ta1_x_yyy_yyzzz_1[i] * fe_0 + ta_yyy_xyyzzz_1[i] + ta1_x_yyy_xyyzzz_0[i] * pa_x[i] -
                                 ta1_x_yyy_xyyzzz_1[i] * pc_x[i];

        ta1_x_xyyy_xyzzzz_0[i] = ta1_x_yyy_yzzzz_0[i] * fe_0 - ta1_x_yyy_yzzzz_1[i] * fe_0 + ta_yyy_xyzzzz_1[i] + ta1_x_yyy_xyzzzz_0[i] * pa_x[i] -
                                 ta1_x_yyy_xyzzzz_1[i] * pc_x[i];

        ta1_x_xyyy_xzzzzz_0[i] =
            2.0 * ta1_x_xy_xzzzzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xzzzzz_1[i] * fe_0 + ta1_x_xyy_xzzzzz_0[i] * pa_y[i] - ta1_x_xyy_xzzzzz_1[i] * pc_y[i];

        ta1_x_xyyy_yyyyyy_0[i] = ta_yyy_yyyyyy_1[i] + ta1_x_yyy_yyyyyy_0[i] * pa_x[i] - ta1_x_yyy_yyyyyy_1[i] * pc_x[i];

        ta1_x_xyyy_yyyyyz_0[i] = ta_yyy_yyyyyz_1[i] + ta1_x_yyy_yyyyyz_0[i] * pa_x[i] - ta1_x_yyy_yyyyyz_1[i] * pc_x[i];

        ta1_x_xyyy_yyyyzz_0[i] = ta_yyy_yyyyzz_1[i] + ta1_x_yyy_yyyyzz_0[i] * pa_x[i] - ta1_x_yyy_yyyyzz_1[i] * pc_x[i];

        ta1_x_xyyy_yyyzzz_0[i] = ta_yyy_yyyzzz_1[i] + ta1_x_yyy_yyyzzz_0[i] * pa_x[i] - ta1_x_yyy_yyyzzz_1[i] * pc_x[i];

        ta1_x_xyyy_yyzzzz_0[i] = ta_yyy_yyzzzz_1[i] + ta1_x_yyy_yyzzzz_0[i] * pa_x[i] - ta1_x_yyy_yyzzzz_1[i] * pc_x[i];

        ta1_x_xyyy_yzzzzz_0[i] = ta_yyy_yzzzzz_1[i] + ta1_x_yyy_yzzzzz_0[i] * pa_x[i] - ta1_x_yyy_yzzzzz_1[i] * pc_x[i];

        ta1_x_xyyy_zzzzzz_0[i] = ta_yyy_zzzzzz_1[i] + ta1_x_yyy_zzzzzz_0[i] * pa_x[i] - ta1_x_yyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 196-224 components of targeted buffer : GI

    auto ta1_x_xyyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 196);

    auto ta1_x_xyyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 197);

    auto ta1_x_xyyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 198);

    auto ta1_x_xyyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 199);

    auto ta1_x_xyyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 200);

    auto ta1_x_xyyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 201);

    auto ta1_x_xyyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 202);

    auto ta1_x_xyyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 203);

    auto ta1_x_xyyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 204);

    auto ta1_x_xyyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 205);

    auto ta1_x_xyyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 206);

    auto ta1_x_xyyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 207);

    auto ta1_x_xyyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 208);

    auto ta1_x_xyyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 209);

    auto ta1_x_xyyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 210);

    auto ta1_x_xyyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 211);

    auto ta1_x_xyyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 212);

    auto ta1_x_xyyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 213);

    auto ta1_x_xyyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 214);

    auto ta1_x_xyyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 215);

    auto ta1_x_xyyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 216);

    auto ta1_x_xyyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 217);

    auto ta1_x_xyyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 218);

    auto ta1_x_xyyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 219);

    auto ta1_x_xyyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 220);

    auto ta1_x_xyyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 221);

    auto ta1_x_xyyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 222);

    auto ta1_x_xyyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 223);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pa_z,                \
                             pc_x,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_x_xyy_xxxxxx_0,  \
                             ta1_x_xyy_xxxxxx_1,  \
                             ta1_x_xyy_xxxxxy_0,  \
                             ta1_x_xyy_xxxxxy_1,  \
                             ta1_x_xyy_xxxxy_0,   \
                             ta1_x_xyy_xxxxy_1,   \
                             ta1_x_xyy_xxxxyy_0,  \
                             ta1_x_xyy_xxxxyy_1,  \
                             ta1_x_xyy_xxxxyz_0,  \
                             ta1_x_xyy_xxxxyz_1,  \
                             ta1_x_xyy_xxxyy_0,   \
                             ta1_x_xyy_xxxyy_1,   \
                             ta1_x_xyy_xxxyyy_0,  \
                             ta1_x_xyy_xxxyyy_1,  \
                             ta1_x_xyy_xxxyyz_0,  \
                             ta1_x_xyy_xxxyyz_1,  \
                             ta1_x_xyy_xxxyz_0,   \
                             ta1_x_xyy_xxxyz_1,   \
                             ta1_x_xyy_xxxyzz_0,  \
                             ta1_x_xyy_xxxyzz_1,  \
                             ta1_x_xyy_xxyyy_0,   \
                             ta1_x_xyy_xxyyy_1,   \
                             ta1_x_xyy_xxyyyy_0,  \
                             ta1_x_xyy_xxyyyy_1,  \
                             ta1_x_xyy_xxyyyz_0,  \
                             ta1_x_xyy_xxyyyz_1,  \
                             ta1_x_xyy_xxyyz_0,   \
                             ta1_x_xyy_xxyyz_1,   \
                             ta1_x_xyy_xxyyzz_0,  \
                             ta1_x_xyy_xxyyzz_1,  \
                             ta1_x_xyy_xxyzz_0,   \
                             ta1_x_xyy_xxyzz_1,   \
                             ta1_x_xyy_xxyzzz_0,  \
                             ta1_x_xyy_xxyzzz_1,  \
                             ta1_x_xyy_xyyyy_0,   \
                             ta1_x_xyy_xyyyy_1,   \
                             ta1_x_xyy_xyyyyy_0,  \
                             ta1_x_xyy_xyyyyy_1,  \
                             ta1_x_xyy_xyyyyz_0,  \
                             ta1_x_xyy_xyyyyz_1,  \
                             ta1_x_xyy_xyyyz_0,   \
                             ta1_x_xyy_xyyyz_1,   \
                             ta1_x_xyy_xyyyzz_0,  \
                             ta1_x_xyy_xyyyzz_1,  \
                             ta1_x_xyy_xyyzz_0,   \
                             ta1_x_xyy_xyyzz_1,   \
                             ta1_x_xyy_xyyzzz_0,  \
                             ta1_x_xyy_xyyzzz_1,  \
                             ta1_x_xyy_xyzzz_0,   \
                             ta1_x_xyy_xyzzz_1,   \
                             ta1_x_xyy_xyzzzz_0,  \
                             ta1_x_xyy_xyzzzz_1,  \
                             ta1_x_xyy_yyyyyy_0,  \
                             ta1_x_xyy_yyyyyy_1,  \
                             ta1_x_xyyz_xxxxxx_0, \
                             ta1_x_xyyz_xxxxxy_0, \
                             ta1_x_xyyz_xxxxxz_0, \
                             ta1_x_xyyz_xxxxyy_0, \
                             ta1_x_xyyz_xxxxyz_0, \
                             ta1_x_xyyz_xxxxzz_0, \
                             ta1_x_xyyz_xxxyyy_0, \
                             ta1_x_xyyz_xxxyyz_0, \
                             ta1_x_xyyz_xxxyzz_0, \
                             ta1_x_xyyz_xxxzzz_0, \
                             ta1_x_xyyz_xxyyyy_0, \
                             ta1_x_xyyz_xxyyyz_0, \
                             ta1_x_xyyz_xxyyzz_0, \
                             ta1_x_xyyz_xxyzzz_0, \
                             ta1_x_xyyz_xxzzzz_0, \
                             ta1_x_xyyz_xyyyyy_0, \
                             ta1_x_xyyz_xyyyyz_0, \
                             ta1_x_xyyz_xyyyzz_0, \
                             ta1_x_xyyz_xyyzzz_0, \
                             ta1_x_xyyz_xyzzzz_0, \
                             ta1_x_xyyz_xzzzzz_0, \
                             ta1_x_xyyz_yyyyyy_0, \
                             ta1_x_xyyz_yyyyyz_0, \
                             ta1_x_xyyz_yyyyzz_0, \
                             ta1_x_xyyz_yyyzzz_0, \
                             ta1_x_xyyz_yyzzzz_0, \
                             ta1_x_xyyz_yzzzzz_0, \
                             ta1_x_xyyz_zzzzzz_0, \
                             ta1_x_xyz_xxxxxz_0,  \
                             ta1_x_xyz_xxxxxz_1,  \
                             ta1_x_xyz_xxxxzz_0,  \
                             ta1_x_xyz_xxxxzz_1,  \
                             ta1_x_xyz_xxxzzz_0,  \
                             ta1_x_xyz_xxxzzz_1,  \
                             ta1_x_xyz_xxzzzz_0,  \
                             ta1_x_xyz_xxzzzz_1,  \
                             ta1_x_xyz_xzzzzz_0,  \
                             ta1_x_xyz_xzzzzz_1,  \
                             ta1_x_xz_xxxxxz_0,   \
                             ta1_x_xz_xxxxxz_1,   \
                             ta1_x_xz_xxxxzz_0,   \
                             ta1_x_xz_xxxxzz_1,   \
                             ta1_x_xz_xxxzzz_0,   \
                             ta1_x_xz_xxxzzz_1,   \
                             ta1_x_xz_xxzzzz_0,   \
                             ta1_x_xz_xxzzzz_1,   \
                             ta1_x_xz_xzzzzz_0,   \
                             ta1_x_xz_xzzzzz_1,   \
                             ta1_x_yyz_yyyyyz_0,  \
                             ta1_x_yyz_yyyyyz_1,  \
                             ta1_x_yyz_yyyyzz_0,  \
                             ta1_x_yyz_yyyyzz_1,  \
                             ta1_x_yyz_yyyzzz_0,  \
                             ta1_x_yyz_yyyzzz_1,  \
                             ta1_x_yyz_yyzzzz_0,  \
                             ta1_x_yyz_yyzzzz_1,  \
                             ta1_x_yyz_yzzzzz_0,  \
                             ta1_x_yyz_yzzzzz_1,  \
                             ta1_x_yyz_zzzzzz_0,  \
                             ta1_x_yyz_zzzzzz_1,  \
                             ta_yyz_yyyyyz_1,     \
                             ta_yyz_yyyyzz_1,     \
                             ta_yyz_yyyzzz_1,     \
                             ta_yyz_yyzzzz_1,     \
                             ta_yyz_yzzzzz_1,     \
                             ta_yyz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyz_xxxxxx_0[i] = ta1_x_xyy_xxxxxx_0[i] * pa_z[i] - ta1_x_xyy_xxxxxx_1[i] * pc_z[i];

        ta1_x_xyyz_xxxxxy_0[i] = ta1_x_xyy_xxxxxy_0[i] * pa_z[i] - ta1_x_xyy_xxxxxy_1[i] * pc_z[i];

        ta1_x_xyyz_xxxxxz_0[i] =
            ta1_x_xz_xxxxxz_0[i] * fe_0 - ta1_x_xz_xxxxxz_1[i] * fe_0 + ta1_x_xyz_xxxxxz_0[i] * pa_y[i] - ta1_x_xyz_xxxxxz_1[i] * pc_y[i];

        ta1_x_xyyz_xxxxyy_0[i] = ta1_x_xyy_xxxxyy_0[i] * pa_z[i] - ta1_x_xyy_xxxxyy_1[i] * pc_z[i];

        ta1_x_xyyz_xxxxyz_0[i] =
            ta1_x_xyy_xxxxy_0[i] * fe_0 - ta1_x_xyy_xxxxy_1[i] * fe_0 + ta1_x_xyy_xxxxyz_0[i] * pa_z[i] - ta1_x_xyy_xxxxyz_1[i] * pc_z[i];

        ta1_x_xyyz_xxxxzz_0[i] =
            ta1_x_xz_xxxxzz_0[i] * fe_0 - ta1_x_xz_xxxxzz_1[i] * fe_0 + ta1_x_xyz_xxxxzz_0[i] * pa_y[i] - ta1_x_xyz_xxxxzz_1[i] * pc_y[i];

        ta1_x_xyyz_xxxyyy_0[i] = ta1_x_xyy_xxxyyy_0[i] * pa_z[i] - ta1_x_xyy_xxxyyy_1[i] * pc_z[i];

        ta1_x_xyyz_xxxyyz_0[i] =
            ta1_x_xyy_xxxyy_0[i] * fe_0 - ta1_x_xyy_xxxyy_1[i] * fe_0 + ta1_x_xyy_xxxyyz_0[i] * pa_z[i] - ta1_x_xyy_xxxyyz_1[i] * pc_z[i];

        ta1_x_xyyz_xxxyzz_0[i] =
            2.0 * ta1_x_xyy_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_xxxyz_1[i] * fe_0 + ta1_x_xyy_xxxyzz_0[i] * pa_z[i] - ta1_x_xyy_xxxyzz_1[i] * pc_z[i];

        ta1_x_xyyz_xxxzzz_0[i] =
            ta1_x_xz_xxxzzz_0[i] * fe_0 - ta1_x_xz_xxxzzz_1[i] * fe_0 + ta1_x_xyz_xxxzzz_0[i] * pa_y[i] - ta1_x_xyz_xxxzzz_1[i] * pc_y[i];

        ta1_x_xyyz_xxyyyy_0[i] = ta1_x_xyy_xxyyyy_0[i] * pa_z[i] - ta1_x_xyy_xxyyyy_1[i] * pc_z[i];

        ta1_x_xyyz_xxyyyz_0[i] =
            ta1_x_xyy_xxyyy_0[i] * fe_0 - ta1_x_xyy_xxyyy_1[i] * fe_0 + ta1_x_xyy_xxyyyz_0[i] * pa_z[i] - ta1_x_xyy_xxyyyz_1[i] * pc_z[i];

        ta1_x_xyyz_xxyyzz_0[i] =
            2.0 * ta1_x_xyy_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_xxyyz_1[i] * fe_0 + ta1_x_xyy_xxyyzz_0[i] * pa_z[i] - ta1_x_xyy_xxyyzz_1[i] * pc_z[i];

        ta1_x_xyyz_xxyzzz_0[i] =
            3.0 * ta1_x_xyy_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xxyzz_1[i] * fe_0 + ta1_x_xyy_xxyzzz_0[i] * pa_z[i] - ta1_x_xyy_xxyzzz_1[i] * pc_z[i];

        ta1_x_xyyz_xxzzzz_0[i] =
            ta1_x_xz_xxzzzz_0[i] * fe_0 - ta1_x_xz_xxzzzz_1[i] * fe_0 + ta1_x_xyz_xxzzzz_0[i] * pa_y[i] - ta1_x_xyz_xxzzzz_1[i] * pc_y[i];

        ta1_x_xyyz_xyyyyy_0[i] = ta1_x_xyy_xyyyyy_0[i] * pa_z[i] - ta1_x_xyy_xyyyyy_1[i] * pc_z[i];

        ta1_x_xyyz_xyyyyz_0[i] =
            ta1_x_xyy_xyyyy_0[i] * fe_0 - ta1_x_xyy_xyyyy_1[i] * fe_0 + ta1_x_xyy_xyyyyz_0[i] * pa_z[i] - ta1_x_xyy_xyyyyz_1[i] * pc_z[i];

        ta1_x_xyyz_xyyyzz_0[i] =
            2.0 * ta1_x_xyy_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_xyyyz_1[i] * fe_0 + ta1_x_xyy_xyyyzz_0[i] * pa_z[i] - ta1_x_xyy_xyyyzz_1[i] * pc_z[i];

        ta1_x_xyyz_xyyzzz_0[i] =
            3.0 * ta1_x_xyy_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xyyzz_1[i] * fe_0 + ta1_x_xyy_xyyzzz_0[i] * pa_z[i] - ta1_x_xyy_xyyzzz_1[i] * pc_z[i];

        ta1_x_xyyz_xyzzzz_0[i] =
            4.0 * ta1_x_xyy_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_xyy_xyzzz_1[i] * fe_0 + ta1_x_xyy_xyzzzz_0[i] * pa_z[i] - ta1_x_xyy_xyzzzz_1[i] * pc_z[i];

        ta1_x_xyyz_xzzzzz_0[i] =
            ta1_x_xz_xzzzzz_0[i] * fe_0 - ta1_x_xz_xzzzzz_1[i] * fe_0 + ta1_x_xyz_xzzzzz_0[i] * pa_y[i] - ta1_x_xyz_xzzzzz_1[i] * pc_y[i];

        ta1_x_xyyz_yyyyyy_0[i] = ta1_x_xyy_yyyyyy_0[i] * pa_z[i] - ta1_x_xyy_yyyyyy_1[i] * pc_z[i];

        ta1_x_xyyz_yyyyyz_0[i] = ta_yyz_yyyyyz_1[i] + ta1_x_yyz_yyyyyz_0[i] * pa_x[i] - ta1_x_yyz_yyyyyz_1[i] * pc_x[i];

        ta1_x_xyyz_yyyyzz_0[i] = ta_yyz_yyyyzz_1[i] + ta1_x_yyz_yyyyzz_0[i] * pa_x[i] - ta1_x_yyz_yyyyzz_1[i] * pc_x[i];

        ta1_x_xyyz_yyyzzz_0[i] = ta_yyz_yyyzzz_1[i] + ta1_x_yyz_yyyzzz_0[i] * pa_x[i] - ta1_x_yyz_yyyzzz_1[i] * pc_x[i];

        ta1_x_xyyz_yyzzzz_0[i] = ta_yyz_yyzzzz_1[i] + ta1_x_yyz_yyzzzz_0[i] * pa_x[i] - ta1_x_yyz_yyzzzz_1[i] * pc_x[i];

        ta1_x_xyyz_yzzzzz_0[i] = ta_yyz_yzzzzz_1[i] + ta1_x_yyz_yzzzzz_0[i] * pa_x[i] - ta1_x_yyz_yzzzzz_1[i] * pc_x[i];

        ta1_x_xyyz_zzzzzz_0[i] = ta_yyz_zzzzzz_1[i] + ta1_x_yyz_zzzzzz_0[i] * pa_x[i] - ta1_x_yyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 224-252 components of targeted buffer : GI

    auto ta1_x_xyzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 224);

    auto ta1_x_xyzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 225);

    auto ta1_x_xyzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 226);

    auto ta1_x_xyzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 227);

    auto ta1_x_xyzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 228);

    auto ta1_x_xyzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 229);

    auto ta1_x_xyzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 230);

    auto ta1_x_xyzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 231);

    auto ta1_x_xyzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 232);

    auto ta1_x_xyzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 233);

    auto ta1_x_xyzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 234);

    auto ta1_x_xyzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 235);

    auto ta1_x_xyzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 236);

    auto ta1_x_xyzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 237);

    auto ta1_x_xyzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 238);

    auto ta1_x_xyzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 239);

    auto ta1_x_xyzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 240);

    auto ta1_x_xyzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 241);

    auto ta1_x_xyzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 242);

    auto ta1_x_xyzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 243);

    auto ta1_x_xyzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 244);

    auto ta1_x_xyzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 245);

    auto ta1_x_xyzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 246);

    auto ta1_x_xyzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 247);

    auto ta1_x_xyzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 248);

    auto ta1_x_xyzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 249);

    auto ta1_x_xyzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 250);

    auto ta1_x_xyzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 251);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_x_xyzz_xxxxxx_0, \
                             ta1_x_xyzz_xxxxxy_0, \
                             ta1_x_xyzz_xxxxxz_0, \
                             ta1_x_xyzz_xxxxyy_0, \
                             ta1_x_xyzz_xxxxyz_0, \
                             ta1_x_xyzz_xxxxzz_0, \
                             ta1_x_xyzz_xxxyyy_0, \
                             ta1_x_xyzz_xxxyyz_0, \
                             ta1_x_xyzz_xxxyzz_0, \
                             ta1_x_xyzz_xxxzzz_0, \
                             ta1_x_xyzz_xxyyyy_0, \
                             ta1_x_xyzz_xxyyyz_0, \
                             ta1_x_xyzz_xxyyzz_0, \
                             ta1_x_xyzz_xxyzzz_0, \
                             ta1_x_xyzz_xxzzzz_0, \
                             ta1_x_xyzz_xyyyyy_0, \
                             ta1_x_xyzz_xyyyyz_0, \
                             ta1_x_xyzz_xyyyzz_0, \
                             ta1_x_xyzz_xyyzzz_0, \
                             ta1_x_xyzz_xyzzzz_0, \
                             ta1_x_xyzz_xzzzzz_0, \
                             ta1_x_xyzz_yyyyyy_0, \
                             ta1_x_xyzz_yyyyyz_0, \
                             ta1_x_xyzz_yyyyzz_0, \
                             ta1_x_xyzz_yyyzzz_0, \
                             ta1_x_xyzz_yyzzzz_0, \
                             ta1_x_xyzz_yzzzzz_0, \
                             ta1_x_xyzz_zzzzzz_0, \
                             ta1_x_xzz_xxxxx_0,   \
                             ta1_x_xzz_xxxxx_1,   \
                             ta1_x_xzz_xxxxxx_0,  \
                             ta1_x_xzz_xxxxxx_1,  \
                             ta1_x_xzz_xxxxxy_0,  \
                             ta1_x_xzz_xxxxxy_1,  \
                             ta1_x_xzz_xxxxxz_0,  \
                             ta1_x_xzz_xxxxxz_1,  \
                             ta1_x_xzz_xxxxy_0,   \
                             ta1_x_xzz_xxxxy_1,   \
                             ta1_x_xzz_xxxxyy_0,  \
                             ta1_x_xzz_xxxxyy_1,  \
                             ta1_x_xzz_xxxxyz_0,  \
                             ta1_x_xzz_xxxxyz_1,  \
                             ta1_x_xzz_xxxxz_0,   \
                             ta1_x_xzz_xxxxz_1,   \
                             ta1_x_xzz_xxxxzz_0,  \
                             ta1_x_xzz_xxxxzz_1,  \
                             ta1_x_xzz_xxxyy_0,   \
                             ta1_x_xzz_xxxyy_1,   \
                             ta1_x_xzz_xxxyyy_0,  \
                             ta1_x_xzz_xxxyyy_1,  \
                             ta1_x_xzz_xxxyyz_0,  \
                             ta1_x_xzz_xxxyyz_1,  \
                             ta1_x_xzz_xxxyz_0,   \
                             ta1_x_xzz_xxxyz_1,   \
                             ta1_x_xzz_xxxyzz_0,  \
                             ta1_x_xzz_xxxyzz_1,  \
                             ta1_x_xzz_xxxzz_0,   \
                             ta1_x_xzz_xxxzz_1,   \
                             ta1_x_xzz_xxxzzz_0,  \
                             ta1_x_xzz_xxxzzz_1,  \
                             ta1_x_xzz_xxyyy_0,   \
                             ta1_x_xzz_xxyyy_1,   \
                             ta1_x_xzz_xxyyyy_0,  \
                             ta1_x_xzz_xxyyyy_1,  \
                             ta1_x_xzz_xxyyyz_0,  \
                             ta1_x_xzz_xxyyyz_1,  \
                             ta1_x_xzz_xxyyz_0,   \
                             ta1_x_xzz_xxyyz_1,   \
                             ta1_x_xzz_xxyyzz_0,  \
                             ta1_x_xzz_xxyyzz_1,  \
                             ta1_x_xzz_xxyzz_0,   \
                             ta1_x_xzz_xxyzz_1,   \
                             ta1_x_xzz_xxyzzz_0,  \
                             ta1_x_xzz_xxyzzz_1,  \
                             ta1_x_xzz_xxzzz_0,   \
                             ta1_x_xzz_xxzzz_1,   \
                             ta1_x_xzz_xxzzzz_0,  \
                             ta1_x_xzz_xxzzzz_1,  \
                             ta1_x_xzz_xyyyy_0,   \
                             ta1_x_xzz_xyyyy_1,   \
                             ta1_x_xzz_xyyyyy_0,  \
                             ta1_x_xzz_xyyyyy_1,  \
                             ta1_x_xzz_xyyyyz_0,  \
                             ta1_x_xzz_xyyyyz_1,  \
                             ta1_x_xzz_xyyyz_0,   \
                             ta1_x_xzz_xyyyz_1,   \
                             ta1_x_xzz_xyyyzz_0,  \
                             ta1_x_xzz_xyyyzz_1,  \
                             ta1_x_xzz_xyyzz_0,   \
                             ta1_x_xzz_xyyzz_1,   \
                             ta1_x_xzz_xyyzzz_0,  \
                             ta1_x_xzz_xyyzzz_1,  \
                             ta1_x_xzz_xyzzz_0,   \
                             ta1_x_xzz_xyzzz_1,   \
                             ta1_x_xzz_xyzzzz_0,  \
                             ta1_x_xzz_xyzzzz_1,  \
                             ta1_x_xzz_xzzzz_0,   \
                             ta1_x_xzz_xzzzz_1,   \
                             ta1_x_xzz_xzzzzz_0,  \
                             ta1_x_xzz_xzzzzz_1,  \
                             ta1_x_xzz_zzzzzz_0,  \
                             ta1_x_xzz_zzzzzz_1,  \
                             ta1_x_yzz_yyyyyy_0,  \
                             ta1_x_yzz_yyyyyy_1,  \
                             ta1_x_yzz_yyyyyz_0,  \
                             ta1_x_yzz_yyyyyz_1,  \
                             ta1_x_yzz_yyyyzz_0,  \
                             ta1_x_yzz_yyyyzz_1,  \
                             ta1_x_yzz_yyyzzz_0,  \
                             ta1_x_yzz_yyyzzz_1,  \
                             ta1_x_yzz_yyzzzz_0,  \
                             ta1_x_yzz_yyzzzz_1,  \
                             ta1_x_yzz_yzzzzz_0,  \
                             ta1_x_yzz_yzzzzz_1,  \
                             ta_yzz_yyyyyy_1,     \
                             ta_yzz_yyyyyz_1,     \
                             ta_yzz_yyyyzz_1,     \
                             ta_yzz_yyyzzz_1,     \
                             ta_yzz_yyzzzz_1,     \
                             ta_yzz_yzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzz_xxxxxx_0[i] = ta1_x_xzz_xxxxxx_0[i] * pa_y[i] - ta1_x_xzz_xxxxxx_1[i] * pc_y[i];

        ta1_x_xyzz_xxxxxy_0[i] =
            ta1_x_xzz_xxxxx_0[i] * fe_0 - ta1_x_xzz_xxxxx_1[i] * fe_0 + ta1_x_xzz_xxxxxy_0[i] * pa_y[i] - ta1_x_xzz_xxxxxy_1[i] * pc_y[i];

        ta1_x_xyzz_xxxxxz_0[i] = ta1_x_xzz_xxxxxz_0[i] * pa_y[i] - ta1_x_xzz_xxxxxz_1[i] * pc_y[i];

        ta1_x_xyzz_xxxxyy_0[i] =
            2.0 * ta1_x_xzz_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_xzz_xxxxy_1[i] * fe_0 + ta1_x_xzz_xxxxyy_0[i] * pa_y[i] - ta1_x_xzz_xxxxyy_1[i] * pc_y[i];

        ta1_x_xyzz_xxxxyz_0[i] =
            ta1_x_xzz_xxxxz_0[i] * fe_0 - ta1_x_xzz_xxxxz_1[i] * fe_0 + ta1_x_xzz_xxxxyz_0[i] * pa_y[i] - ta1_x_xzz_xxxxyz_1[i] * pc_y[i];

        ta1_x_xyzz_xxxxzz_0[i] = ta1_x_xzz_xxxxzz_0[i] * pa_y[i] - ta1_x_xzz_xxxxzz_1[i] * pc_y[i];

        ta1_x_xyzz_xxxyyy_0[i] =
            3.0 * ta1_x_xzz_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxxyy_1[i] * fe_0 + ta1_x_xzz_xxxyyy_0[i] * pa_y[i] - ta1_x_xzz_xxxyyy_1[i] * pc_y[i];

        ta1_x_xyzz_xxxyyz_0[i] =
            2.0 * ta1_x_xzz_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xzz_xxxyz_1[i] * fe_0 + ta1_x_xzz_xxxyyz_0[i] * pa_y[i] - ta1_x_xzz_xxxyyz_1[i] * pc_y[i];

        ta1_x_xyzz_xxxyzz_0[i] =
            ta1_x_xzz_xxxzz_0[i] * fe_0 - ta1_x_xzz_xxxzz_1[i] * fe_0 + ta1_x_xzz_xxxyzz_0[i] * pa_y[i] - ta1_x_xzz_xxxyzz_1[i] * pc_y[i];

        ta1_x_xyzz_xxxzzz_0[i] = ta1_x_xzz_xxxzzz_0[i] * pa_y[i] - ta1_x_xzz_xxxzzz_1[i] * pc_y[i];

        ta1_x_xyzz_xxyyyy_0[i] =
            4.0 * ta1_x_xzz_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_xzz_xxyyy_1[i] * fe_0 + ta1_x_xzz_xxyyyy_0[i] * pa_y[i] - ta1_x_xzz_xxyyyy_1[i] * pc_y[i];

        ta1_x_xyzz_xxyyyz_0[i] =
            3.0 * ta1_x_xzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxyyz_1[i] * fe_0 + ta1_x_xzz_xxyyyz_0[i] * pa_y[i] - ta1_x_xzz_xxyyyz_1[i] * pc_y[i];

        ta1_x_xyzz_xxyyzz_0[i] =
            2.0 * ta1_x_xzz_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_xxyzz_1[i] * fe_0 + ta1_x_xzz_xxyyzz_0[i] * pa_y[i] - ta1_x_xzz_xxyyzz_1[i] * pc_y[i];

        ta1_x_xyzz_xxyzzz_0[i] =
            ta1_x_xzz_xxzzz_0[i] * fe_0 - ta1_x_xzz_xxzzz_1[i] * fe_0 + ta1_x_xzz_xxyzzz_0[i] * pa_y[i] - ta1_x_xzz_xxyzzz_1[i] * pc_y[i];

        ta1_x_xyzz_xxzzzz_0[i] = ta1_x_xzz_xxzzzz_0[i] * pa_y[i] - ta1_x_xzz_xxzzzz_1[i] * pc_y[i];

        ta1_x_xyzz_xyyyyy_0[i] =
            5.0 * ta1_x_xzz_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_xzz_xyyyy_1[i] * fe_0 + ta1_x_xzz_xyyyyy_0[i] * pa_y[i] - ta1_x_xzz_xyyyyy_1[i] * pc_y[i];

        ta1_x_xyzz_xyyyyz_0[i] =
            4.0 * ta1_x_xzz_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_xzz_xyyyz_1[i] * fe_0 + ta1_x_xzz_xyyyyz_0[i] * pa_y[i] - ta1_x_xzz_xyyyyz_1[i] * pc_y[i];

        ta1_x_xyzz_xyyyzz_0[i] =
            3.0 * ta1_x_xzz_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xzz_xyyzz_1[i] * fe_0 + ta1_x_xzz_xyyyzz_0[i] * pa_y[i] - ta1_x_xzz_xyyyzz_1[i] * pc_y[i];

        ta1_x_xyzz_xyyzzz_0[i] =
            2.0 * ta1_x_xzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_xyzzz_1[i] * fe_0 + ta1_x_xzz_xyyzzz_0[i] * pa_y[i] - ta1_x_xzz_xyyzzz_1[i] * pc_y[i];

        ta1_x_xyzz_xyzzzz_0[i] =
            ta1_x_xzz_xzzzz_0[i] * fe_0 - ta1_x_xzz_xzzzz_1[i] * fe_0 + ta1_x_xzz_xyzzzz_0[i] * pa_y[i] - ta1_x_xzz_xyzzzz_1[i] * pc_y[i];

        ta1_x_xyzz_xzzzzz_0[i] = ta1_x_xzz_xzzzzz_0[i] * pa_y[i] - ta1_x_xzz_xzzzzz_1[i] * pc_y[i];

        ta1_x_xyzz_yyyyyy_0[i] = ta_yzz_yyyyyy_1[i] + ta1_x_yzz_yyyyyy_0[i] * pa_x[i] - ta1_x_yzz_yyyyyy_1[i] * pc_x[i];

        ta1_x_xyzz_yyyyyz_0[i] = ta_yzz_yyyyyz_1[i] + ta1_x_yzz_yyyyyz_0[i] * pa_x[i] - ta1_x_yzz_yyyyyz_1[i] * pc_x[i];

        ta1_x_xyzz_yyyyzz_0[i] = ta_yzz_yyyyzz_1[i] + ta1_x_yzz_yyyyzz_0[i] * pa_x[i] - ta1_x_yzz_yyyyzz_1[i] * pc_x[i];

        ta1_x_xyzz_yyyzzz_0[i] = ta_yzz_yyyzzz_1[i] + ta1_x_yzz_yyyzzz_0[i] * pa_x[i] - ta1_x_yzz_yyyzzz_1[i] * pc_x[i];

        ta1_x_xyzz_yyzzzz_0[i] = ta_yzz_yyzzzz_1[i] + ta1_x_yzz_yyzzzz_0[i] * pa_x[i] - ta1_x_yzz_yyzzzz_1[i] * pc_x[i];

        ta1_x_xyzz_yzzzzz_0[i] = ta_yzz_yzzzzz_1[i] + ta1_x_yzz_yzzzzz_0[i] * pa_x[i] - ta1_x_yzz_yzzzzz_1[i] * pc_x[i];

        ta1_x_xyzz_zzzzzz_0[i] = ta1_x_xzz_zzzzzz_0[i] * pa_y[i] - ta1_x_xzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 252-280 components of targeted buffer : GI

    auto ta1_x_xzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 252);

    auto ta1_x_xzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 253);

    auto ta1_x_xzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 254);

    auto ta1_x_xzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 255);

    auto ta1_x_xzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 256);

    auto ta1_x_xzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 257);

    auto ta1_x_xzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 258);

    auto ta1_x_xzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 259);

    auto ta1_x_xzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 260);

    auto ta1_x_xzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 261);

    auto ta1_x_xzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 262);

    auto ta1_x_xzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 263);

    auto ta1_x_xzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 264);

    auto ta1_x_xzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 265);

    auto ta1_x_xzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 266);

    auto ta1_x_xzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 267);

    auto ta1_x_xzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 268);

    auto ta1_x_xzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 269);

    auto ta1_x_xzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 270);

    auto ta1_x_xzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 271);

    auto ta1_x_xzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 272);

    auto ta1_x_xzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 273);

    auto ta1_x_xzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 274);

    auto ta1_x_xzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 275);

    auto ta1_x_xzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 276);

    auto ta1_x_xzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 277);

    auto ta1_x_xzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 278);

    auto ta1_x_xzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 279);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_x_xz_xxxxxx_0,   \
                             ta1_x_xz_xxxxxx_1,   \
                             ta1_x_xz_xxxxxy_0,   \
                             ta1_x_xz_xxxxxy_1,   \
                             ta1_x_xz_xxxxyy_0,   \
                             ta1_x_xz_xxxxyy_1,   \
                             ta1_x_xz_xxxyyy_0,   \
                             ta1_x_xz_xxxyyy_1,   \
                             ta1_x_xz_xxyyyy_0,   \
                             ta1_x_xz_xxyyyy_1,   \
                             ta1_x_xz_xyyyyy_0,   \
                             ta1_x_xz_xyyyyy_1,   \
                             ta1_x_xzz_xxxxxx_0,  \
                             ta1_x_xzz_xxxxxx_1,  \
                             ta1_x_xzz_xxxxxy_0,  \
                             ta1_x_xzz_xxxxxy_1,  \
                             ta1_x_xzz_xxxxyy_0,  \
                             ta1_x_xzz_xxxxyy_1,  \
                             ta1_x_xzz_xxxyyy_0,  \
                             ta1_x_xzz_xxxyyy_1,  \
                             ta1_x_xzz_xxyyyy_0,  \
                             ta1_x_xzz_xxyyyy_1,  \
                             ta1_x_xzz_xyyyyy_0,  \
                             ta1_x_xzz_xyyyyy_1,  \
                             ta1_x_xzzz_xxxxxx_0, \
                             ta1_x_xzzz_xxxxxy_0, \
                             ta1_x_xzzz_xxxxxz_0, \
                             ta1_x_xzzz_xxxxyy_0, \
                             ta1_x_xzzz_xxxxyz_0, \
                             ta1_x_xzzz_xxxxzz_0, \
                             ta1_x_xzzz_xxxyyy_0, \
                             ta1_x_xzzz_xxxyyz_0, \
                             ta1_x_xzzz_xxxyzz_0, \
                             ta1_x_xzzz_xxxzzz_0, \
                             ta1_x_xzzz_xxyyyy_0, \
                             ta1_x_xzzz_xxyyyz_0, \
                             ta1_x_xzzz_xxyyzz_0, \
                             ta1_x_xzzz_xxyzzz_0, \
                             ta1_x_xzzz_xxzzzz_0, \
                             ta1_x_xzzz_xyyyyy_0, \
                             ta1_x_xzzz_xyyyyz_0, \
                             ta1_x_xzzz_xyyyzz_0, \
                             ta1_x_xzzz_xyyzzz_0, \
                             ta1_x_xzzz_xyzzzz_0, \
                             ta1_x_xzzz_xzzzzz_0, \
                             ta1_x_xzzz_yyyyyy_0, \
                             ta1_x_xzzz_yyyyyz_0, \
                             ta1_x_xzzz_yyyyzz_0, \
                             ta1_x_xzzz_yyyzzz_0, \
                             ta1_x_xzzz_yyzzzz_0, \
                             ta1_x_xzzz_yzzzzz_0, \
                             ta1_x_xzzz_zzzzzz_0, \
                             ta1_x_zzz_xxxxxz_0,  \
                             ta1_x_zzz_xxxxxz_1,  \
                             ta1_x_zzz_xxxxyz_0,  \
                             ta1_x_zzz_xxxxyz_1,  \
                             ta1_x_zzz_xxxxz_0,   \
                             ta1_x_zzz_xxxxz_1,   \
                             ta1_x_zzz_xxxxzz_0,  \
                             ta1_x_zzz_xxxxzz_1,  \
                             ta1_x_zzz_xxxyyz_0,  \
                             ta1_x_zzz_xxxyyz_1,  \
                             ta1_x_zzz_xxxyz_0,   \
                             ta1_x_zzz_xxxyz_1,   \
                             ta1_x_zzz_xxxyzz_0,  \
                             ta1_x_zzz_xxxyzz_1,  \
                             ta1_x_zzz_xxxzz_0,   \
                             ta1_x_zzz_xxxzz_1,   \
                             ta1_x_zzz_xxxzzz_0,  \
                             ta1_x_zzz_xxxzzz_1,  \
                             ta1_x_zzz_xxyyyz_0,  \
                             ta1_x_zzz_xxyyyz_1,  \
                             ta1_x_zzz_xxyyz_0,   \
                             ta1_x_zzz_xxyyz_1,   \
                             ta1_x_zzz_xxyyzz_0,  \
                             ta1_x_zzz_xxyyzz_1,  \
                             ta1_x_zzz_xxyzz_0,   \
                             ta1_x_zzz_xxyzz_1,   \
                             ta1_x_zzz_xxyzzz_0,  \
                             ta1_x_zzz_xxyzzz_1,  \
                             ta1_x_zzz_xxzzz_0,   \
                             ta1_x_zzz_xxzzz_1,   \
                             ta1_x_zzz_xxzzzz_0,  \
                             ta1_x_zzz_xxzzzz_1,  \
                             ta1_x_zzz_xyyyyz_0,  \
                             ta1_x_zzz_xyyyyz_1,  \
                             ta1_x_zzz_xyyyz_0,   \
                             ta1_x_zzz_xyyyz_1,   \
                             ta1_x_zzz_xyyyzz_0,  \
                             ta1_x_zzz_xyyyzz_1,  \
                             ta1_x_zzz_xyyzz_0,   \
                             ta1_x_zzz_xyyzz_1,   \
                             ta1_x_zzz_xyyzzz_0,  \
                             ta1_x_zzz_xyyzzz_1,  \
                             ta1_x_zzz_xyzzz_0,   \
                             ta1_x_zzz_xyzzz_1,   \
                             ta1_x_zzz_xyzzzz_0,  \
                             ta1_x_zzz_xyzzzz_1,  \
                             ta1_x_zzz_xzzzz_0,   \
                             ta1_x_zzz_xzzzz_1,   \
                             ta1_x_zzz_xzzzzz_0,  \
                             ta1_x_zzz_xzzzzz_1,  \
                             ta1_x_zzz_yyyyyy_0,  \
                             ta1_x_zzz_yyyyyy_1,  \
                             ta1_x_zzz_yyyyyz_0,  \
                             ta1_x_zzz_yyyyyz_1,  \
                             ta1_x_zzz_yyyyz_0,   \
                             ta1_x_zzz_yyyyz_1,   \
                             ta1_x_zzz_yyyyzz_0,  \
                             ta1_x_zzz_yyyyzz_1,  \
                             ta1_x_zzz_yyyzz_0,   \
                             ta1_x_zzz_yyyzz_1,   \
                             ta1_x_zzz_yyyzzz_0,  \
                             ta1_x_zzz_yyyzzz_1,  \
                             ta1_x_zzz_yyzzz_0,   \
                             ta1_x_zzz_yyzzz_1,   \
                             ta1_x_zzz_yyzzzz_0,  \
                             ta1_x_zzz_yyzzzz_1,  \
                             ta1_x_zzz_yzzzz_0,   \
                             ta1_x_zzz_yzzzz_1,   \
                             ta1_x_zzz_yzzzzz_0,  \
                             ta1_x_zzz_yzzzzz_1,  \
                             ta1_x_zzz_zzzzz_0,   \
                             ta1_x_zzz_zzzzz_1,   \
                             ta1_x_zzz_zzzzzz_0,  \
                             ta1_x_zzz_zzzzzz_1,  \
                             ta_zzz_xxxxxz_1,     \
                             ta_zzz_xxxxyz_1,     \
                             ta_zzz_xxxxzz_1,     \
                             ta_zzz_xxxyyz_1,     \
                             ta_zzz_xxxyzz_1,     \
                             ta_zzz_xxxzzz_1,     \
                             ta_zzz_xxyyyz_1,     \
                             ta_zzz_xxyyzz_1,     \
                             ta_zzz_xxyzzz_1,     \
                             ta_zzz_xxzzzz_1,     \
                             ta_zzz_xyyyyz_1,     \
                             ta_zzz_xyyyzz_1,     \
                             ta_zzz_xyyzzz_1,     \
                             ta_zzz_xyzzzz_1,     \
                             ta_zzz_xzzzzz_1,     \
                             ta_zzz_yyyyyy_1,     \
                             ta_zzz_yyyyyz_1,     \
                             ta_zzz_yyyyzz_1,     \
                             ta_zzz_yyyzzz_1,     \
                             ta_zzz_yyzzzz_1,     \
                             ta_zzz_yzzzzz_1,     \
                             ta_zzz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzz_xxxxxx_0[i] =
            2.0 * ta1_x_xz_xxxxxx_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxxxx_1[i] * fe_0 + ta1_x_xzz_xxxxxx_0[i] * pa_z[i] - ta1_x_xzz_xxxxxx_1[i] * pc_z[i];

        ta1_x_xzzz_xxxxxy_0[i] =
            2.0 * ta1_x_xz_xxxxxy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxxxy_1[i] * fe_0 + ta1_x_xzz_xxxxxy_0[i] * pa_z[i] - ta1_x_xzz_xxxxxy_1[i] * pc_z[i];

        ta1_x_xzzz_xxxxxz_0[i] = 5.0 * ta1_x_zzz_xxxxz_0[i] * fe_0 - 5.0 * ta1_x_zzz_xxxxz_1[i] * fe_0 + ta_zzz_xxxxxz_1[i] +
                                 ta1_x_zzz_xxxxxz_0[i] * pa_x[i] - ta1_x_zzz_xxxxxz_1[i] * pc_x[i];

        ta1_x_xzzz_xxxxyy_0[i] =
            2.0 * ta1_x_xz_xxxxyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxxyy_1[i] * fe_0 + ta1_x_xzz_xxxxyy_0[i] * pa_z[i] - ta1_x_xzz_xxxxyy_1[i] * pc_z[i];

        ta1_x_xzzz_xxxxyz_0[i] = 4.0 * ta1_x_zzz_xxxyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxxyz_1[i] * fe_0 + ta_zzz_xxxxyz_1[i] +
                                 ta1_x_zzz_xxxxyz_0[i] * pa_x[i] - ta1_x_zzz_xxxxyz_1[i] * pc_x[i];

        ta1_x_xzzz_xxxxzz_0[i] = 4.0 * ta1_x_zzz_xxxzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxxzz_1[i] * fe_0 + ta_zzz_xxxxzz_1[i] +
                                 ta1_x_zzz_xxxxzz_0[i] * pa_x[i] - ta1_x_zzz_xxxxzz_1[i] * pc_x[i];

        ta1_x_xzzz_xxxyyy_0[i] =
            2.0 * ta1_x_xz_xxxyyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxyyy_1[i] * fe_0 + ta1_x_xzz_xxxyyy_0[i] * pa_z[i] - ta1_x_xzz_xxxyyy_1[i] * pc_z[i];

        ta1_x_xzzz_xxxyyz_0[i] = 3.0 * ta1_x_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxyyz_1[i] * fe_0 + ta_zzz_xxxyyz_1[i] +
                                 ta1_x_zzz_xxxyyz_0[i] * pa_x[i] - ta1_x_zzz_xxxyyz_1[i] * pc_x[i];

        ta1_x_xzzz_xxxyzz_0[i] = 3.0 * ta1_x_zzz_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxyzz_1[i] * fe_0 + ta_zzz_xxxyzz_1[i] +
                                 ta1_x_zzz_xxxyzz_0[i] * pa_x[i] - ta1_x_zzz_xxxyzz_1[i] * pc_x[i];

        ta1_x_xzzz_xxxzzz_0[i] = 3.0 * ta1_x_zzz_xxzzz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxzzz_1[i] * fe_0 + ta_zzz_xxxzzz_1[i] +
                                 ta1_x_zzz_xxxzzz_0[i] * pa_x[i] - ta1_x_zzz_xxxzzz_1[i] * pc_x[i];

        ta1_x_xzzz_xxyyyy_0[i] =
            2.0 * ta1_x_xz_xxyyyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxyyyy_1[i] * fe_0 + ta1_x_xzz_xxyyyy_0[i] * pa_z[i] - ta1_x_xzz_xxyyyy_1[i] * pc_z[i];

        ta1_x_xzzz_xxyyyz_0[i] = 2.0 * ta1_x_zzz_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyyyz_1[i] * fe_0 + ta_zzz_xxyyyz_1[i] +
                                 ta1_x_zzz_xxyyyz_0[i] * pa_x[i] - ta1_x_zzz_xxyyyz_1[i] * pc_x[i];

        ta1_x_xzzz_xxyyzz_0[i] = 2.0 * ta1_x_zzz_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyyzz_1[i] * fe_0 + ta_zzz_xxyyzz_1[i] +
                                 ta1_x_zzz_xxyyzz_0[i] * pa_x[i] - ta1_x_zzz_xxyyzz_1[i] * pc_x[i];

        ta1_x_xzzz_xxyzzz_0[i] = 2.0 * ta1_x_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyzzz_1[i] * fe_0 + ta_zzz_xxyzzz_1[i] +
                                 ta1_x_zzz_xxyzzz_0[i] * pa_x[i] - ta1_x_zzz_xxyzzz_1[i] * pc_x[i];

        ta1_x_xzzz_xxzzzz_0[i] = 2.0 * ta1_x_zzz_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xzzzz_1[i] * fe_0 + ta_zzz_xxzzzz_1[i] +
                                 ta1_x_zzz_xxzzzz_0[i] * pa_x[i] - ta1_x_zzz_xxzzzz_1[i] * pc_x[i];

        ta1_x_xzzz_xyyyyy_0[i] =
            2.0 * ta1_x_xz_xyyyyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xyyyyy_1[i] * fe_0 + ta1_x_xzz_xyyyyy_0[i] * pa_z[i] - ta1_x_xzz_xyyyyy_1[i] * pc_z[i];

        ta1_x_xzzz_xyyyyz_0[i] = ta1_x_zzz_yyyyz_0[i] * fe_0 - ta1_x_zzz_yyyyz_1[i] * fe_0 + ta_zzz_xyyyyz_1[i] + ta1_x_zzz_xyyyyz_0[i] * pa_x[i] -
                                 ta1_x_zzz_xyyyyz_1[i] * pc_x[i];

        ta1_x_xzzz_xyyyzz_0[i] = ta1_x_zzz_yyyzz_0[i] * fe_0 - ta1_x_zzz_yyyzz_1[i] * fe_0 + ta_zzz_xyyyzz_1[i] + ta1_x_zzz_xyyyzz_0[i] * pa_x[i] -
                                 ta1_x_zzz_xyyyzz_1[i] * pc_x[i];

        ta1_x_xzzz_xyyzzz_0[i] = ta1_x_zzz_yyzzz_0[i] * fe_0 - ta1_x_zzz_yyzzz_1[i] * fe_0 + ta_zzz_xyyzzz_1[i] + ta1_x_zzz_xyyzzz_0[i] * pa_x[i] -
                                 ta1_x_zzz_xyyzzz_1[i] * pc_x[i];

        ta1_x_xzzz_xyzzzz_0[i] = ta1_x_zzz_yzzzz_0[i] * fe_0 - ta1_x_zzz_yzzzz_1[i] * fe_0 + ta_zzz_xyzzzz_1[i] + ta1_x_zzz_xyzzzz_0[i] * pa_x[i] -
                                 ta1_x_zzz_xyzzzz_1[i] * pc_x[i];

        ta1_x_xzzz_xzzzzz_0[i] = ta1_x_zzz_zzzzz_0[i] * fe_0 - ta1_x_zzz_zzzzz_1[i] * fe_0 + ta_zzz_xzzzzz_1[i] + ta1_x_zzz_xzzzzz_0[i] * pa_x[i] -
                                 ta1_x_zzz_xzzzzz_1[i] * pc_x[i];

        ta1_x_xzzz_yyyyyy_0[i] = ta_zzz_yyyyyy_1[i] + ta1_x_zzz_yyyyyy_0[i] * pa_x[i] - ta1_x_zzz_yyyyyy_1[i] * pc_x[i];

        ta1_x_xzzz_yyyyyz_0[i] = ta_zzz_yyyyyz_1[i] + ta1_x_zzz_yyyyyz_0[i] * pa_x[i] - ta1_x_zzz_yyyyyz_1[i] * pc_x[i];

        ta1_x_xzzz_yyyyzz_0[i] = ta_zzz_yyyyzz_1[i] + ta1_x_zzz_yyyyzz_0[i] * pa_x[i] - ta1_x_zzz_yyyyzz_1[i] * pc_x[i];

        ta1_x_xzzz_yyyzzz_0[i] = ta_zzz_yyyzzz_1[i] + ta1_x_zzz_yyyzzz_0[i] * pa_x[i] - ta1_x_zzz_yyyzzz_1[i] * pc_x[i];

        ta1_x_xzzz_yyzzzz_0[i] = ta_zzz_yyzzzz_1[i] + ta1_x_zzz_yyzzzz_0[i] * pa_x[i] - ta1_x_zzz_yyzzzz_1[i] * pc_x[i];

        ta1_x_xzzz_yzzzzz_0[i] = ta_zzz_yzzzzz_1[i] + ta1_x_zzz_yzzzzz_0[i] * pa_x[i] - ta1_x_zzz_yzzzzz_1[i] * pc_x[i];

        ta1_x_xzzz_zzzzzz_0[i] = ta_zzz_zzzzzz_1[i] + ta1_x_zzz_zzzzzz_0[i] * pa_x[i] - ta1_x_zzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 280-308 components of targeted buffer : GI

    auto ta1_x_yyyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 280);

    auto ta1_x_yyyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 281);

    auto ta1_x_yyyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 282);

    auto ta1_x_yyyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 283);

    auto ta1_x_yyyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 284);

    auto ta1_x_yyyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 285);

    auto ta1_x_yyyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 286);

    auto ta1_x_yyyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 287);

    auto ta1_x_yyyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 288);

    auto ta1_x_yyyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 289);

    auto ta1_x_yyyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 290);

    auto ta1_x_yyyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 291);

    auto ta1_x_yyyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 292);

    auto ta1_x_yyyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 293);

    auto ta1_x_yyyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 294);

    auto ta1_x_yyyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 295);

    auto ta1_x_yyyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 296);

    auto ta1_x_yyyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 297);

    auto ta1_x_yyyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 298);

    auto ta1_x_yyyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 299);

    auto ta1_x_yyyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 300);

    auto ta1_x_yyyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 301);

    auto ta1_x_yyyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 302);

    auto ta1_x_yyyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 303);

    auto ta1_x_yyyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 304);

    auto ta1_x_yyyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 305);

    auto ta1_x_yyyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 306);

    auto ta1_x_yyyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 307);

#pragma omp simd aligned(pa_y,                    \
                             pc_y,                \
                             ta1_x_yy_xxxxxx_0,   \
                             ta1_x_yy_xxxxxx_1,   \
                             ta1_x_yy_xxxxxy_0,   \
                             ta1_x_yy_xxxxxy_1,   \
                             ta1_x_yy_xxxxxz_0,   \
                             ta1_x_yy_xxxxxz_1,   \
                             ta1_x_yy_xxxxyy_0,   \
                             ta1_x_yy_xxxxyy_1,   \
                             ta1_x_yy_xxxxyz_0,   \
                             ta1_x_yy_xxxxyz_1,   \
                             ta1_x_yy_xxxxzz_0,   \
                             ta1_x_yy_xxxxzz_1,   \
                             ta1_x_yy_xxxyyy_0,   \
                             ta1_x_yy_xxxyyy_1,   \
                             ta1_x_yy_xxxyyz_0,   \
                             ta1_x_yy_xxxyyz_1,   \
                             ta1_x_yy_xxxyzz_0,   \
                             ta1_x_yy_xxxyzz_1,   \
                             ta1_x_yy_xxxzzz_0,   \
                             ta1_x_yy_xxxzzz_1,   \
                             ta1_x_yy_xxyyyy_0,   \
                             ta1_x_yy_xxyyyy_1,   \
                             ta1_x_yy_xxyyyz_0,   \
                             ta1_x_yy_xxyyyz_1,   \
                             ta1_x_yy_xxyyzz_0,   \
                             ta1_x_yy_xxyyzz_1,   \
                             ta1_x_yy_xxyzzz_0,   \
                             ta1_x_yy_xxyzzz_1,   \
                             ta1_x_yy_xxzzzz_0,   \
                             ta1_x_yy_xxzzzz_1,   \
                             ta1_x_yy_xyyyyy_0,   \
                             ta1_x_yy_xyyyyy_1,   \
                             ta1_x_yy_xyyyyz_0,   \
                             ta1_x_yy_xyyyyz_1,   \
                             ta1_x_yy_xyyyzz_0,   \
                             ta1_x_yy_xyyyzz_1,   \
                             ta1_x_yy_xyyzzz_0,   \
                             ta1_x_yy_xyyzzz_1,   \
                             ta1_x_yy_xyzzzz_0,   \
                             ta1_x_yy_xyzzzz_1,   \
                             ta1_x_yy_xzzzzz_0,   \
                             ta1_x_yy_xzzzzz_1,   \
                             ta1_x_yy_yyyyyy_0,   \
                             ta1_x_yy_yyyyyy_1,   \
                             ta1_x_yy_yyyyyz_0,   \
                             ta1_x_yy_yyyyyz_1,   \
                             ta1_x_yy_yyyyzz_0,   \
                             ta1_x_yy_yyyyzz_1,   \
                             ta1_x_yy_yyyzzz_0,   \
                             ta1_x_yy_yyyzzz_1,   \
                             ta1_x_yy_yyzzzz_0,   \
                             ta1_x_yy_yyzzzz_1,   \
                             ta1_x_yy_yzzzzz_0,   \
                             ta1_x_yy_yzzzzz_1,   \
                             ta1_x_yy_zzzzzz_0,   \
                             ta1_x_yy_zzzzzz_1,   \
                             ta1_x_yyy_xxxxx_0,   \
                             ta1_x_yyy_xxxxx_1,   \
                             ta1_x_yyy_xxxxxx_0,  \
                             ta1_x_yyy_xxxxxx_1,  \
                             ta1_x_yyy_xxxxxy_0,  \
                             ta1_x_yyy_xxxxxy_1,  \
                             ta1_x_yyy_xxxxxz_0,  \
                             ta1_x_yyy_xxxxxz_1,  \
                             ta1_x_yyy_xxxxy_0,   \
                             ta1_x_yyy_xxxxy_1,   \
                             ta1_x_yyy_xxxxyy_0,  \
                             ta1_x_yyy_xxxxyy_1,  \
                             ta1_x_yyy_xxxxyz_0,  \
                             ta1_x_yyy_xxxxyz_1,  \
                             ta1_x_yyy_xxxxz_0,   \
                             ta1_x_yyy_xxxxz_1,   \
                             ta1_x_yyy_xxxxzz_0,  \
                             ta1_x_yyy_xxxxzz_1,  \
                             ta1_x_yyy_xxxyy_0,   \
                             ta1_x_yyy_xxxyy_1,   \
                             ta1_x_yyy_xxxyyy_0,  \
                             ta1_x_yyy_xxxyyy_1,  \
                             ta1_x_yyy_xxxyyz_0,  \
                             ta1_x_yyy_xxxyyz_1,  \
                             ta1_x_yyy_xxxyz_0,   \
                             ta1_x_yyy_xxxyz_1,   \
                             ta1_x_yyy_xxxyzz_0,  \
                             ta1_x_yyy_xxxyzz_1,  \
                             ta1_x_yyy_xxxzz_0,   \
                             ta1_x_yyy_xxxzz_1,   \
                             ta1_x_yyy_xxxzzz_0,  \
                             ta1_x_yyy_xxxzzz_1,  \
                             ta1_x_yyy_xxyyy_0,   \
                             ta1_x_yyy_xxyyy_1,   \
                             ta1_x_yyy_xxyyyy_0,  \
                             ta1_x_yyy_xxyyyy_1,  \
                             ta1_x_yyy_xxyyyz_0,  \
                             ta1_x_yyy_xxyyyz_1,  \
                             ta1_x_yyy_xxyyz_0,   \
                             ta1_x_yyy_xxyyz_1,   \
                             ta1_x_yyy_xxyyzz_0,  \
                             ta1_x_yyy_xxyyzz_1,  \
                             ta1_x_yyy_xxyzz_0,   \
                             ta1_x_yyy_xxyzz_1,   \
                             ta1_x_yyy_xxyzzz_0,  \
                             ta1_x_yyy_xxyzzz_1,  \
                             ta1_x_yyy_xxzzz_0,   \
                             ta1_x_yyy_xxzzz_1,   \
                             ta1_x_yyy_xxzzzz_0,  \
                             ta1_x_yyy_xxzzzz_1,  \
                             ta1_x_yyy_xyyyy_0,   \
                             ta1_x_yyy_xyyyy_1,   \
                             ta1_x_yyy_xyyyyy_0,  \
                             ta1_x_yyy_xyyyyy_1,  \
                             ta1_x_yyy_xyyyyz_0,  \
                             ta1_x_yyy_xyyyyz_1,  \
                             ta1_x_yyy_xyyyz_0,   \
                             ta1_x_yyy_xyyyz_1,   \
                             ta1_x_yyy_xyyyzz_0,  \
                             ta1_x_yyy_xyyyzz_1,  \
                             ta1_x_yyy_xyyzz_0,   \
                             ta1_x_yyy_xyyzz_1,   \
                             ta1_x_yyy_xyyzzz_0,  \
                             ta1_x_yyy_xyyzzz_1,  \
                             ta1_x_yyy_xyzzz_0,   \
                             ta1_x_yyy_xyzzz_1,   \
                             ta1_x_yyy_xyzzzz_0,  \
                             ta1_x_yyy_xyzzzz_1,  \
                             ta1_x_yyy_xzzzz_0,   \
                             ta1_x_yyy_xzzzz_1,   \
                             ta1_x_yyy_xzzzzz_0,  \
                             ta1_x_yyy_xzzzzz_1,  \
                             ta1_x_yyy_yyyyy_0,   \
                             ta1_x_yyy_yyyyy_1,   \
                             ta1_x_yyy_yyyyyy_0,  \
                             ta1_x_yyy_yyyyyy_1,  \
                             ta1_x_yyy_yyyyyz_0,  \
                             ta1_x_yyy_yyyyyz_1,  \
                             ta1_x_yyy_yyyyz_0,   \
                             ta1_x_yyy_yyyyz_1,   \
                             ta1_x_yyy_yyyyzz_0,  \
                             ta1_x_yyy_yyyyzz_1,  \
                             ta1_x_yyy_yyyzz_0,   \
                             ta1_x_yyy_yyyzz_1,   \
                             ta1_x_yyy_yyyzzz_0,  \
                             ta1_x_yyy_yyyzzz_1,  \
                             ta1_x_yyy_yyzzz_0,   \
                             ta1_x_yyy_yyzzz_1,   \
                             ta1_x_yyy_yyzzzz_0,  \
                             ta1_x_yyy_yyzzzz_1,  \
                             ta1_x_yyy_yzzzz_0,   \
                             ta1_x_yyy_yzzzz_1,   \
                             ta1_x_yyy_yzzzzz_0,  \
                             ta1_x_yyy_yzzzzz_1,  \
                             ta1_x_yyy_zzzzz_0,   \
                             ta1_x_yyy_zzzzz_1,   \
                             ta1_x_yyy_zzzzzz_0,  \
                             ta1_x_yyy_zzzzzz_1,  \
                             ta1_x_yyyy_xxxxxx_0, \
                             ta1_x_yyyy_xxxxxy_0, \
                             ta1_x_yyyy_xxxxxz_0, \
                             ta1_x_yyyy_xxxxyy_0, \
                             ta1_x_yyyy_xxxxyz_0, \
                             ta1_x_yyyy_xxxxzz_0, \
                             ta1_x_yyyy_xxxyyy_0, \
                             ta1_x_yyyy_xxxyyz_0, \
                             ta1_x_yyyy_xxxyzz_0, \
                             ta1_x_yyyy_xxxzzz_0, \
                             ta1_x_yyyy_xxyyyy_0, \
                             ta1_x_yyyy_xxyyyz_0, \
                             ta1_x_yyyy_xxyyzz_0, \
                             ta1_x_yyyy_xxyzzz_0, \
                             ta1_x_yyyy_xxzzzz_0, \
                             ta1_x_yyyy_xyyyyy_0, \
                             ta1_x_yyyy_xyyyyz_0, \
                             ta1_x_yyyy_xyyyzz_0, \
                             ta1_x_yyyy_xyyzzz_0, \
                             ta1_x_yyyy_xyzzzz_0, \
                             ta1_x_yyyy_xzzzzz_0, \
                             ta1_x_yyyy_yyyyyy_0, \
                             ta1_x_yyyy_yyyyyz_0, \
                             ta1_x_yyyy_yyyyzz_0, \
                             ta1_x_yyyy_yyyzzz_0, \
                             ta1_x_yyyy_yyzzzz_0, \
                             ta1_x_yyyy_yzzzzz_0, \
                             ta1_x_yyyy_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyy_xxxxxx_0[i] =
            3.0 * ta1_x_yy_xxxxxx_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxxx_1[i] * fe_0 + ta1_x_yyy_xxxxxx_0[i] * pa_y[i] - ta1_x_yyy_xxxxxx_1[i] * pc_y[i];

        ta1_x_yyyy_xxxxxy_0[i] = 3.0 * ta1_x_yy_xxxxxy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxxy_1[i] * fe_0 + ta1_x_yyy_xxxxx_0[i] * fe_0 -
                                 ta1_x_yyy_xxxxx_1[i] * fe_0 + ta1_x_yyy_xxxxxy_0[i] * pa_y[i] - ta1_x_yyy_xxxxxy_1[i] * pc_y[i];

        ta1_x_yyyy_xxxxxz_0[i] =
            3.0 * ta1_x_yy_xxxxxz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxxz_1[i] * fe_0 + ta1_x_yyy_xxxxxz_0[i] * pa_y[i] - ta1_x_yyy_xxxxxz_1[i] * pc_y[i];

        ta1_x_yyyy_xxxxyy_0[i] = 3.0 * ta1_x_yy_xxxxyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxyy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxxxy_0[i] * fe_0 -
                                 2.0 * ta1_x_yyy_xxxxy_1[i] * fe_0 + ta1_x_yyy_xxxxyy_0[i] * pa_y[i] - ta1_x_yyy_xxxxyy_1[i] * pc_y[i];

        ta1_x_yyyy_xxxxyz_0[i] = 3.0 * ta1_x_yy_xxxxyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxyz_1[i] * fe_0 + ta1_x_yyy_xxxxz_0[i] * fe_0 -
                                 ta1_x_yyy_xxxxz_1[i] * fe_0 + ta1_x_yyy_xxxxyz_0[i] * pa_y[i] - ta1_x_yyy_xxxxyz_1[i] * pc_y[i];

        ta1_x_yyyy_xxxxzz_0[i] =
            3.0 * ta1_x_yy_xxxxzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxzz_1[i] * fe_0 + ta1_x_yyy_xxxxzz_0[i] * pa_y[i] - ta1_x_yyy_xxxxzz_1[i] * pc_y[i];

        ta1_x_yyyy_xxxyyy_0[i] = 3.0 * ta1_x_yy_xxxyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxyyy_1[i] * fe_0 + 3.0 * ta1_x_yyy_xxxyy_0[i] * fe_0 -
                                 3.0 * ta1_x_yyy_xxxyy_1[i] * fe_0 + ta1_x_yyy_xxxyyy_0[i] * pa_y[i] - ta1_x_yyy_xxxyyy_1[i] * pc_y[i];

        ta1_x_yyyy_xxxyyz_0[i] = 3.0 * ta1_x_yy_xxxyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxyyz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_x_yyy_xxxyz_1[i] * fe_0 + ta1_x_yyy_xxxyyz_0[i] * pa_y[i] - ta1_x_yyy_xxxyyz_1[i] * pc_y[i];

        ta1_x_yyyy_xxxyzz_0[i] = 3.0 * ta1_x_yy_xxxyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxyzz_1[i] * fe_0 + ta1_x_yyy_xxxzz_0[i] * fe_0 -
                                 ta1_x_yyy_xxxzz_1[i] * fe_0 + ta1_x_yyy_xxxyzz_0[i] * pa_y[i] - ta1_x_yyy_xxxyzz_1[i] * pc_y[i];

        ta1_x_yyyy_xxxzzz_0[i] =
            3.0 * ta1_x_yy_xxxzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxzzz_1[i] * fe_0 + ta1_x_yyy_xxxzzz_0[i] * pa_y[i] - ta1_x_yyy_xxxzzz_1[i] * pc_y[i];

        ta1_x_yyyy_xxyyyy_0[i] = 3.0 * ta1_x_yy_xxyyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyyyy_1[i] * fe_0 + 4.0 * ta1_x_yyy_xxyyy_0[i] * fe_0 -
                                 4.0 * ta1_x_yyy_xxyyy_1[i] * fe_0 + ta1_x_yyy_xxyyyy_0[i] * pa_y[i] - ta1_x_yyy_xxyyyy_1[i] * pc_y[i];

        ta1_x_yyyy_xxyyyz_0[i] = 3.0 * ta1_x_yy_xxyyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyyyz_1[i] * fe_0 + 3.0 * ta1_x_yyy_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_x_yyy_xxyyz_1[i] * fe_0 + ta1_x_yyy_xxyyyz_0[i] * pa_y[i] - ta1_x_yyy_xxyyyz_1[i] * pc_y[i];

        ta1_x_yyyy_xxyyzz_0[i] = 3.0 * ta1_x_yy_xxyyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxyzz_0[i] * fe_0 -
                                 2.0 * ta1_x_yyy_xxyzz_1[i] * fe_0 + ta1_x_yyy_xxyyzz_0[i] * pa_y[i] - ta1_x_yyy_xxyyzz_1[i] * pc_y[i];

        ta1_x_yyyy_xxyzzz_0[i] = 3.0 * ta1_x_yy_xxyzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyzzz_1[i] * fe_0 + ta1_x_yyy_xxzzz_0[i] * fe_0 -
                                 ta1_x_yyy_xxzzz_1[i] * fe_0 + ta1_x_yyy_xxyzzz_0[i] * pa_y[i] - ta1_x_yyy_xxyzzz_1[i] * pc_y[i];

        ta1_x_yyyy_xxzzzz_0[i] =
            3.0 * ta1_x_yy_xxzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxzzzz_1[i] * fe_0 + ta1_x_yyy_xxzzzz_0[i] * pa_y[i] - ta1_x_yyy_xxzzzz_1[i] * pc_y[i];

        ta1_x_yyyy_xyyyyy_0[i] = 3.0 * ta1_x_yy_xyyyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyyyy_1[i] * fe_0 + 5.0 * ta1_x_yyy_xyyyy_0[i] * fe_0 -
                                 5.0 * ta1_x_yyy_xyyyy_1[i] * fe_0 + ta1_x_yyy_xyyyyy_0[i] * pa_y[i] - ta1_x_yyy_xyyyyy_1[i] * pc_y[i];

        ta1_x_yyyy_xyyyyz_0[i] = 3.0 * ta1_x_yy_xyyyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyyyz_1[i] * fe_0 + 4.0 * ta1_x_yyy_xyyyz_0[i] * fe_0 -
                                 4.0 * ta1_x_yyy_xyyyz_1[i] * fe_0 + ta1_x_yyy_xyyyyz_0[i] * pa_y[i] - ta1_x_yyy_xyyyyz_1[i] * pc_y[i];

        ta1_x_yyyy_xyyyzz_0[i] = 3.0 * ta1_x_yy_xyyyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyyzz_1[i] * fe_0 + 3.0 * ta1_x_yyy_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_yyy_xyyzz_1[i] * fe_0 + ta1_x_yyy_xyyyzz_0[i] * pa_y[i] - ta1_x_yyy_xyyyzz_1[i] * pc_y[i];

        ta1_x_yyyy_xyyzzz_0[i] = 3.0 * ta1_x_yy_xyyzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyzzz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_x_yyy_xyzzz_1[i] * fe_0 + ta1_x_yyy_xyyzzz_0[i] * pa_y[i] - ta1_x_yyy_xyyzzz_1[i] * pc_y[i];

        ta1_x_yyyy_xyzzzz_0[i] = 3.0 * ta1_x_yy_xyzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyzzzz_1[i] * fe_0 + ta1_x_yyy_xzzzz_0[i] * fe_0 -
                                 ta1_x_yyy_xzzzz_1[i] * fe_0 + ta1_x_yyy_xyzzzz_0[i] * pa_y[i] - ta1_x_yyy_xyzzzz_1[i] * pc_y[i];

        ta1_x_yyyy_xzzzzz_0[i] =
            3.0 * ta1_x_yy_xzzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xzzzzz_1[i] * fe_0 + ta1_x_yyy_xzzzzz_0[i] * pa_y[i] - ta1_x_yyy_xzzzzz_1[i] * pc_y[i];

        ta1_x_yyyy_yyyyyy_0[i] = 3.0 * ta1_x_yy_yyyyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyyyy_1[i] * fe_0 + 6.0 * ta1_x_yyy_yyyyy_0[i] * fe_0 -
                                 6.0 * ta1_x_yyy_yyyyy_1[i] * fe_0 + ta1_x_yyy_yyyyyy_0[i] * pa_y[i] - ta1_x_yyy_yyyyyy_1[i] * pc_y[i];

        ta1_x_yyyy_yyyyyz_0[i] = 3.0 * ta1_x_yy_yyyyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyyyz_1[i] * fe_0 + 5.0 * ta1_x_yyy_yyyyz_0[i] * fe_0 -
                                 5.0 * ta1_x_yyy_yyyyz_1[i] * fe_0 + ta1_x_yyy_yyyyyz_0[i] * pa_y[i] - ta1_x_yyy_yyyyyz_1[i] * pc_y[i];

        ta1_x_yyyy_yyyyzz_0[i] = 3.0 * ta1_x_yy_yyyyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyyzz_1[i] * fe_0 + 4.0 * ta1_x_yyy_yyyzz_0[i] * fe_0 -
                                 4.0 * ta1_x_yyy_yyyzz_1[i] * fe_0 + ta1_x_yyy_yyyyzz_0[i] * pa_y[i] - ta1_x_yyy_yyyyzz_1[i] * pc_y[i];

        ta1_x_yyyy_yyyzzz_0[i] = 3.0 * ta1_x_yy_yyyzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyzzz_1[i] * fe_0 + 3.0 * ta1_x_yyy_yyzzz_0[i] * fe_0 -
                                 3.0 * ta1_x_yyy_yyzzz_1[i] * fe_0 + ta1_x_yyy_yyyzzz_0[i] * pa_y[i] - ta1_x_yyy_yyyzzz_1[i] * pc_y[i];

        ta1_x_yyyy_yyzzzz_0[i] = 3.0 * ta1_x_yy_yyzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyzzzz_1[i] * fe_0 + 2.0 * ta1_x_yyy_yzzzz_0[i] * fe_0 -
                                 2.0 * ta1_x_yyy_yzzzz_1[i] * fe_0 + ta1_x_yyy_yyzzzz_0[i] * pa_y[i] - ta1_x_yyy_yyzzzz_1[i] * pc_y[i];

        ta1_x_yyyy_yzzzzz_0[i] = 3.0 * ta1_x_yy_yzzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yzzzzz_1[i] * fe_0 + ta1_x_yyy_zzzzz_0[i] * fe_0 -
                                 ta1_x_yyy_zzzzz_1[i] * fe_0 + ta1_x_yyy_yzzzzz_0[i] * pa_y[i] - ta1_x_yyy_yzzzzz_1[i] * pc_y[i];

        ta1_x_yyyy_zzzzzz_0[i] =
            3.0 * ta1_x_yy_zzzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_zzzzzz_1[i] * fe_0 + ta1_x_yyy_zzzzzz_0[i] * pa_y[i] - ta1_x_yyy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 308-336 components of targeted buffer : GI

    auto ta1_x_yyyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 308);

    auto ta1_x_yyyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 309);

    auto ta1_x_yyyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 310);

    auto ta1_x_yyyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 311);

    auto ta1_x_yyyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 312);

    auto ta1_x_yyyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 313);

    auto ta1_x_yyyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 314);

    auto ta1_x_yyyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 315);

    auto ta1_x_yyyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 316);

    auto ta1_x_yyyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 317);

    auto ta1_x_yyyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 318);

    auto ta1_x_yyyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 319);

    auto ta1_x_yyyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 320);

    auto ta1_x_yyyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 321);

    auto ta1_x_yyyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 322);

    auto ta1_x_yyyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 323);

    auto ta1_x_yyyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 324);

    auto ta1_x_yyyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 325);

    auto ta1_x_yyyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 326);

    auto ta1_x_yyyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 327);

    auto ta1_x_yyyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 328);

    auto ta1_x_yyyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 329);

    auto ta1_x_yyyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 330);

    auto ta1_x_yyyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 331);

    auto ta1_x_yyyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 332);

    auto ta1_x_yyyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 333);

    auto ta1_x_yyyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 334);

    auto ta1_x_yyyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 335);

#pragma omp simd aligned(pa_y,                    \
                             pa_z,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_x_yyy_xxxxxx_0,  \
                             ta1_x_yyy_xxxxxx_1,  \
                             ta1_x_yyy_xxxxxy_0,  \
                             ta1_x_yyy_xxxxxy_1,  \
                             ta1_x_yyy_xxxxy_0,   \
                             ta1_x_yyy_xxxxy_1,   \
                             ta1_x_yyy_xxxxyy_0,  \
                             ta1_x_yyy_xxxxyy_1,  \
                             ta1_x_yyy_xxxxyz_0,  \
                             ta1_x_yyy_xxxxyz_1,  \
                             ta1_x_yyy_xxxyy_0,   \
                             ta1_x_yyy_xxxyy_1,   \
                             ta1_x_yyy_xxxyyy_0,  \
                             ta1_x_yyy_xxxyyy_1,  \
                             ta1_x_yyy_xxxyyz_0,  \
                             ta1_x_yyy_xxxyyz_1,  \
                             ta1_x_yyy_xxxyz_0,   \
                             ta1_x_yyy_xxxyz_1,   \
                             ta1_x_yyy_xxxyzz_0,  \
                             ta1_x_yyy_xxxyzz_1,  \
                             ta1_x_yyy_xxyyy_0,   \
                             ta1_x_yyy_xxyyy_1,   \
                             ta1_x_yyy_xxyyyy_0,  \
                             ta1_x_yyy_xxyyyy_1,  \
                             ta1_x_yyy_xxyyyz_0,  \
                             ta1_x_yyy_xxyyyz_1,  \
                             ta1_x_yyy_xxyyz_0,   \
                             ta1_x_yyy_xxyyz_1,   \
                             ta1_x_yyy_xxyyzz_0,  \
                             ta1_x_yyy_xxyyzz_1,  \
                             ta1_x_yyy_xxyzz_0,   \
                             ta1_x_yyy_xxyzz_1,   \
                             ta1_x_yyy_xxyzzz_0,  \
                             ta1_x_yyy_xxyzzz_1,  \
                             ta1_x_yyy_xyyyy_0,   \
                             ta1_x_yyy_xyyyy_1,   \
                             ta1_x_yyy_xyyyyy_0,  \
                             ta1_x_yyy_xyyyyy_1,  \
                             ta1_x_yyy_xyyyyz_0,  \
                             ta1_x_yyy_xyyyyz_1,  \
                             ta1_x_yyy_xyyyz_0,   \
                             ta1_x_yyy_xyyyz_1,   \
                             ta1_x_yyy_xyyyzz_0,  \
                             ta1_x_yyy_xyyyzz_1,  \
                             ta1_x_yyy_xyyzz_0,   \
                             ta1_x_yyy_xyyzz_1,   \
                             ta1_x_yyy_xyyzzz_0,  \
                             ta1_x_yyy_xyyzzz_1,  \
                             ta1_x_yyy_xyzzz_0,   \
                             ta1_x_yyy_xyzzz_1,   \
                             ta1_x_yyy_xyzzzz_0,  \
                             ta1_x_yyy_xyzzzz_1,  \
                             ta1_x_yyy_yyyyy_0,   \
                             ta1_x_yyy_yyyyy_1,   \
                             ta1_x_yyy_yyyyyy_0,  \
                             ta1_x_yyy_yyyyyy_1,  \
                             ta1_x_yyy_yyyyyz_0,  \
                             ta1_x_yyy_yyyyyz_1,  \
                             ta1_x_yyy_yyyyz_0,   \
                             ta1_x_yyy_yyyyz_1,   \
                             ta1_x_yyy_yyyyzz_0,  \
                             ta1_x_yyy_yyyyzz_1,  \
                             ta1_x_yyy_yyyzz_0,   \
                             ta1_x_yyy_yyyzz_1,   \
                             ta1_x_yyy_yyyzzz_0,  \
                             ta1_x_yyy_yyyzzz_1,  \
                             ta1_x_yyy_yyzzz_0,   \
                             ta1_x_yyy_yyzzz_1,   \
                             ta1_x_yyy_yyzzzz_0,  \
                             ta1_x_yyy_yyzzzz_1,  \
                             ta1_x_yyy_yzzzz_0,   \
                             ta1_x_yyy_yzzzz_1,   \
                             ta1_x_yyy_yzzzzz_0,  \
                             ta1_x_yyy_yzzzzz_1,  \
                             ta1_x_yyyz_xxxxxx_0, \
                             ta1_x_yyyz_xxxxxy_0, \
                             ta1_x_yyyz_xxxxxz_0, \
                             ta1_x_yyyz_xxxxyy_0, \
                             ta1_x_yyyz_xxxxyz_0, \
                             ta1_x_yyyz_xxxxzz_0, \
                             ta1_x_yyyz_xxxyyy_0, \
                             ta1_x_yyyz_xxxyyz_0, \
                             ta1_x_yyyz_xxxyzz_0, \
                             ta1_x_yyyz_xxxzzz_0, \
                             ta1_x_yyyz_xxyyyy_0, \
                             ta1_x_yyyz_xxyyyz_0, \
                             ta1_x_yyyz_xxyyzz_0, \
                             ta1_x_yyyz_xxyzzz_0, \
                             ta1_x_yyyz_xxzzzz_0, \
                             ta1_x_yyyz_xyyyyy_0, \
                             ta1_x_yyyz_xyyyyz_0, \
                             ta1_x_yyyz_xyyyzz_0, \
                             ta1_x_yyyz_xyyzzz_0, \
                             ta1_x_yyyz_xyzzzz_0, \
                             ta1_x_yyyz_xzzzzz_0, \
                             ta1_x_yyyz_yyyyyy_0, \
                             ta1_x_yyyz_yyyyyz_0, \
                             ta1_x_yyyz_yyyyzz_0, \
                             ta1_x_yyyz_yyyzzz_0, \
                             ta1_x_yyyz_yyzzzz_0, \
                             ta1_x_yyyz_yzzzzz_0, \
                             ta1_x_yyyz_zzzzzz_0, \
                             ta1_x_yyz_xxxxxz_0,  \
                             ta1_x_yyz_xxxxxz_1,  \
                             ta1_x_yyz_xxxxzz_0,  \
                             ta1_x_yyz_xxxxzz_1,  \
                             ta1_x_yyz_xxxzzz_0,  \
                             ta1_x_yyz_xxxzzz_1,  \
                             ta1_x_yyz_xxzzzz_0,  \
                             ta1_x_yyz_xxzzzz_1,  \
                             ta1_x_yyz_xzzzzz_0,  \
                             ta1_x_yyz_xzzzzz_1,  \
                             ta1_x_yyz_zzzzzz_0,  \
                             ta1_x_yyz_zzzzzz_1,  \
                             ta1_x_yz_xxxxxz_0,   \
                             ta1_x_yz_xxxxxz_1,   \
                             ta1_x_yz_xxxxzz_0,   \
                             ta1_x_yz_xxxxzz_1,   \
                             ta1_x_yz_xxxzzz_0,   \
                             ta1_x_yz_xxxzzz_1,   \
                             ta1_x_yz_xxzzzz_0,   \
                             ta1_x_yz_xxzzzz_1,   \
                             ta1_x_yz_xzzzzz_0,   \
                             ta1_x_yz_xzzzzz_1,   \
                             ta1_x_yz_zzzzzz_0,   \
                             ta1_x_yz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyz_xxxxxx_0[i] = ta1_x_yyy_xxxxxx_0[i] * pa_z[i] - ta1_x_yyy_xxxxxx_1[i] * pc_z[i];

        ta1_x_yyyz_xxxxxy_0[i] = ta1_x_yyy_xxxxxy_0[i] * pa_z[i] - ta1_x_yyy_xxxxxy_1[i] * pc_z[i];

        ta1_x_yyyz_xxxxxz_0[i] =
            2.0 * ta1_x_yz_xxxxxz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxxxxz_1[i] * fe_0 + ta1_x_yyz_xxxxxz_0[i] * pa_y[i] - ta1_x_yyz_xxxxxz_1[i] * pc_y[i];

        ta1_x_yyyz_xxxxyy_0[i] = ta1_x_yyy_xxxxyy_0[i] * pa_z[i] - ta1_x_yyy_xxxxyy_1[i] * pc_z[i];

        ta1_x_yyyz_xxxxyz_0[i] =
            ta1_x_yyy_xxxxy_0[i] * fe_0 - ta1_x_yyy_xxxxy_1[i] * fe_0 + ta1_x_yyy_xxxxyz_0[i] * pa_z[i] - ta1_x_yyy_xxxxyz_1[i] * pc_z[i];

        ta1_x_yyyz_xxxxzz_0[i] =
            2.0 * ta1_x_yz_xxxxzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxxxzz_1[i] * fe_0 + ta1_x_yyz_xxxxzz_0[i] * pa_y[i] - ta1_x_yyz_xxxxzz_1[i] * pc_y[i];

        ta1_x_yyyz_xxxyyy_0[i] = ta1_x_yyy_xxxyyy_0[i] * pa_z[i] - ta1_x_yyy_xxxyyy_1[i] * pc_z[i];

        ta1_x_yyyz_xxxyyz_0[i] =
            ta1_x_yyy_xxxyy_0[i] * fe_0 - ta1_x_yyy_xxxyy_1[i] * fe_0 + ta1_x_yyy_xxxyyz_0[i] * pa_z[i] - ta1_x_yyy_xxxyyz_1[i] * pc_z[i];

        ta1_x_yyyz_xxxyzz_0[i] =
            2.0 * ta1_x_yyy_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xxxyz_1[i] * fe_0 + ta1_x_yyy_xxxyzz_0[i] * pa_z[i] - ta1_x_yyy_xxxyzz_1[i] * pc_z[i];

        ta1_x_yyyz_xxxzzz_0[i] =
            2.0 * ta1_x_yz_xxxzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxxzzz_1[i] * fe_0 + ta1_x_yyz_xxxzzz_0[i] * pa_y[i] - ta1_x_yyz_xxxzzz_1[i] * pc_y[i];

        ta1_x_yyyz_xxyyyy_0[i] = ta1_x_yyy_xxyyyy_0[i] * pa_z[i] - ta1_x_yyy_xxyyyy_1[i] * pc_z[i];

        ta1_x_yyyz_xxyyyz_0[i] =
            ta1_x_yyy_xxyyy_0[i] * fe_0 - ta1_x_yyy_xxyyy_1[i] * fe_0 + ta1_x_yyy_xxyyyz_0[i] * pa_z[i] - ta1_x_yyy_xxyyyz_1[i] * pc_z[i];

        ta1_x_yyyz_xxyyzz_0[i] =
            2.0 * ta1_x_yyy_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xxyyz_1[i] * fe_0 + ta1_x_yyy_xxyyzz_0[i] * pa_z[i] - ta1_x_yyy_xxyyzz_1[i] * pc_z[i];

        ta1_x_yyyz_xxyzzz_0[i] =
            3.0 * ta1_x_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_yyy_xxyzz_1[i] * fe_0 + ta1_x_yyy_xxyzzz_0[i] * pa_z[i] - ta1_x_yyy_xxyzzz_1[i] * pc_z[i];

        ta1_x_yyyz_xxzzzz_0[i] =
            2.0 * ta1_x_yz_xxzzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxzzzz_1[i] * fe_0 + ta1_x_yyz_xxzzzz_0[i] * pa_y[i] - ta1_x_yyz_xxzzzz_1[i] * pc_y[i];

        ta1_x_yyyz_xyyyyy_0[i] = ta1_x_yyy_xyyyyy_0[i] * pa_z[i] - ta1_x_yyy_xyyyyy_1[i] * pc_z[i];

        ta1_x_yyyz_xyyyyz_0[i] =
            ta1_x_yyy_xyyyy_0[i] * fe_0 - ta1_x_yyy_xyyyy_1[i] * fe_0 + ta1_x_yyy_xyyyyz_0[i] * pa_z[i] - ta1_x_yyy_xyyyyz_1[i] * pc_z[i];

        ta1_x_yyyz_xyyyzz_0[i] =
            2.0 * ta1_x_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyyyz_1[i] * fe_0 + ta1_x_yyy_xyyyzz_0[i] * pa_z[i] - ta1_x_yyy_xyyyzz_1[i] * pc_z[i];

        ta1_x_yyyz_xyyzzz_0[i] =
            3.0 * ta1_x_yyy_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_yyy_xyyzz_1[i] * fe_0 + ta1_x_yyy_xyyzzz_0[i] * pa_z[i] - ta1_x_yyy_xyyzzz_1[i] * pc_z[i];

        ta1_x_yyyz_xyzzzz_0[i] =
            4.0 * ta1_x_yyy_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xyzzz_1[i] * fe_0 + ta1_x_yyy_xyzzzz_0[i] * pa_z[i] - ta1_x_yyy_xyzzzz_1[i] * pc_z[i];

        ta1_x_yyyz_xzzzzz_0[i] =
            2.0 * ta1_x_yz_xzzzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xzzzzz_1[i] * fe_0 + ta1_x_yyz_xzzzzz_0[i] * pa_y[i] - ta1_x_yyz_xzzzzz_1[i] * pc_y[i];

        ta1_x_yyyz_yyyyyy_0[i] = ta1_x_yyy_yyyyyy_0[i] * pa_z[i] - ta1_x_yyy_yyyyyy_1[i] * pc_z[i];

        ta1_x_yyyz_yyyyyz_0[i] =
            ta1_x_yyy_yyyyy_0[i] * fe_0 - ta1_x_yyy_yyyyy_1[i] * fe_0 + ta1_x_yyy_yyyyyz_0[i] * pa_z[i] - ta1_x_yyy_yyyyyz_1[i] * pc_z[i];

        ta1_x_yyyz_yyyyzz_0[i] =
            2.0 * ta1_x_yyy_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_yyyyz_1[i] * fe_0 + ta1_x_yyy_yyyyzz_0[i] * pa_z[i] - ta1_x_yyy_yyyyzz_1[i] * pc_z[i];

        ta1_x_yyyz_yyyzzz_0[i] =
            3.0 * ta1_x_yyy_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_yyy_yyyzz_1[i] * fe_0 + ta1_x_yyy_yyyzzz_0[i] * pa_z[i] - ta1_x_yyy_yyyzzz_1[i] * pc_z[i];

        ta1_x_yyyz_yyzzzz_0[i] =
            4.0 * ta1_x_yyy_yyzzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yyzzz_1[i] * fe_0 + ta1_x_yyy_yyzzzz_0[i] * pa_z[i] - ta1_x_yyy_yyzzzz_1[i] * pc_z[i];

        ta1_x_yyyz_yzzzzz_0[i] =
            5.0 * ta1_x_yyy_yzzzz_0[i] * fe_0 - 5.0 * ta1_x_yyy_yzzzz_1[i] * fe_0 + ta1_x_yyy_yzzzzz_0[i] * pa_z[i] - ta1_x_yyy_yzzzzz_1[i] * pc_z[i];

        ta1_x_yyyz_zzzzzz_0[i] =
            2.0 * ta1_x_yz_zzzzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_zzzzzz_1[i] * fe_0 + ta1_x_yyz_zzzzzz_0[i] * pa_y[i] - ta1_x_yyz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 336-364 components of targeted buffer : GI

    auto ta1_x_yyzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 336);

    auto ta1_x_yyzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 337);

    auto ta1_x_yyzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 338);

    auto ta1_x_yyzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 339);

    auto ta1_x_yyzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 340);

    auto ta1_x_yyzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 341);

    auto ta1_x_yyzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 342);

    auto ta1_x_yyzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 343);

    auto ta1_x_yyzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 344);

    auto ta1_x_yyzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 345);

    auto ta1_x_yyzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 346);

    auto ta1_x_yyzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 347);

    auto ta1_x_yyzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 348);

    auto ta1_x_yyzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 349);

    auto ta1_x_yyzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 350);

    auto ta1_x_yyzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 351);

    auto ta1_x_yyzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 352);

    auto ta1_x_yyzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 353);

    auto ta1_x_yyzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 354);

    auto ta1_x_yyzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 355);

    auto ta1_x_yyzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 356);

    auto ta1_x_yyzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 357);

    auto ta1_x_yyzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 358);

    auto ta1_x_yyzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 359);

    auto ta1_x_yyzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 360);

    auto ta1_x_yyzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 361);

    auto ta1_x_yyzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 362);

    auto ta1_x_yyzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 363);

#pragma omp simd aligned(pa_y,                    \
                             pa_z,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_x_yy_xxxxxy_0,   \
                             ta1_x_yy_xxxxxy_1,   \
                             ta1_x_yy_xxxxyy_0,   \
                             ta1_x_yy_xxxxyy_1,   \
                             ta1_x_yy_xxxyyy_0,   \
                             ta1_x_yy_xxxyyy_1,   \
                             ta1_x_yy_xxyyyy_0,   \
                             ta1_x_yy_xxyyyy_1,   \
                             ta1_x_yy_xyyyyy_0,   \
                             ta1_x_yy_xyyyyy_1,   \
                             ta1_x_yy_yyyyyy_0,   \
                             ta1_x_yy_yyyyyy_1,   \
                             ta1_x_yyz_xxxxxy_0,  \
                             ta1_x_yyz_xxxxxy_1,  \
                             ta1_x_yyz_xxxxyy_0,  \
                             ta1_x_yyz_xxxxyy_1,  \
                             ta1_x_yyz_xxxyyy_0,  \
                             ta1_x_yyz_xxxyyy_1,  \
                             ta1_x_yyz_xxyyyy_0,  \
                             ta1_x_yyz_xxyyyy_1,  \
                             ta1_x_yyz_xyyyyy_0,  \
                             ta1_x_yyz_xyyyyy_1,  \
                             ta1_x_yyz_yyyyyy_0,  \
                             ta1_x_yyz_yyyyyy_1,  \
                             ta1_x_yyzz_xxxxxx_0, \
                             ta1_x_yyzz_xxxxxy_0, \
                             ta1_x_yyzz_xxxxxz_0, \
                             ta1_x_yyzz_xxxxyy_0, \
                             ta1_x_yyzz_xxxxyz_0, \
                             ta1_x_yyzz_xxxxzz_0, \
                             ta1_x_yyzz_xxxyyy_0, \
                             ta1_x_yyzz_xxxyyz_0, \
                             ta1_x_yyzz_xxxyzz_0, \
                             ta1_x_yyzz_xxxzzz_0, \
                             ta1_x_yyzz_xxyyyy_0, \
                             ta1_x_yyzz_xxyyyz_0, \
                             ta1_x_yyzz_xxyyzz_0, \
                             ta1_x_yyzz_xxyzzz_0, \
                             ta1_x_yyzz_xxzzzz_0, \
                             ta1_x_yyzz_xyyyyy_0, \
                             ta1_x_yyzz_xyyyyz_0, \
                             ta1_x_yyzz_xyyyzz_0, \
                             ta1_x_yyzz_xyyzzz_0, \
                             ta1_x_yyzz_xyzzzz_0, \
                             ta1_x_yyzz_xzzzzz_0, \
                             ta1_x_yyzz_yyyyyy_0, \
                             ta1_x_yyzz_yyyyyz_0, \
                             ta1_x_yyzz_yyyyzz_0, \
                             ta1_x_yyzz_yyyzzz_0, \
                             ta1_x_yyzz_yyzzzz_0, \
                             ta1_x_yyzz_yzzzzz_0, \
                             ta1_x_yyzz_zzzzzz_0, \
                             ta1_x_yzz_xxxxxx_0,  \
                             ta1_x_yzz_xxxxxx_1,  \
                             ta1_x_yzz_xxxxxz_0,  \
                             ta1_x_yzz_xxxxxz_1,  \
                             ta1_x_yzz_xxxxyz_0,  \
                             ta1_x_yzz_xxxxyz_1,  \
                             ta1_x_yzz_xxxxz_0,   \
                             ta1_x_yzz_xxxxz_1,   \
                             ta1_x_yzz_xxxxzz_0,  \
                             ta1_x_yzz_xxxxzz_1,  \
                             ta1_x_yzz_xxxyyz_0,  \
                             ta1_x_yzz_xxxyyz_1,  \
                             ta1_x_yzz_xxxyz_0,   \
                             ta1_x_yzz_xxxyz_1,   \
                             ta1_x_yzz_xxxyzz_0,  \
                             ta1_x_yzz_xxxyzz_1,  \
                             ta1_x_yzz_xxxzz_0,   \
                             ta1_x_yzz_xxxzz_1,   \
                             ta1_x_yzz_xxxzzz_0,  \
                             ta1_x_yzz_xxxzzz_1,  \
                             ta1_x_yzz_xxyyyz_0,  \
                             ta1_x_yzz_xxyyyz_1,  \
                             ta1_x_yzz_xxyyz_0,   \
                             ta1_x_yzz_xxyyz_1,   \
                             ta1_x_yzz_xxyyzz_0,  \
                             ta1_x_yzz_xxyyzz_1,  \
                             ta1_x_yzz_xxyzz_0,   \
                             ta1_x_yzz_xxyzz_1,   \
                             ta1_x_yzz_xxyzzz_0,  \
                             ta1_x_yzz_xxyzzz_1,  \
                             ta1_x_yzz_xxzzz_0,   \
                             ta1_x_yzz_xxzzz_1,   \
                             ta1_x_yzz_xxzzzz_0,  \
                             ta1_x_yzz_xxzzzz_1,  \
                             ta1_x_yzz_xyyyyz_0,  \
                             ta1_x_yzz_xyyyyz_1,  \
                             ta1_x_yzz_xyyyz_0,   \
                             ta1_x_yzz_xyyyz_1,   \
                             ta1_x_yzz_xyyyzz_0,  \
                             ta1_x_yzz_xyyyzz_1,  \
                             ta1_x_yzz_xyyzz_0,   \
                             ta1_x_yzz_xyyzz_1,   \
                             ta1_x_yzz_xyyzzz_0,  \
                             ta1_x_yzz_xyyzzz_1,  \
                             ta1_x_yzz_xyzzz_0,   \
                             ta1_x_yzz_xyzzz_1,   \
                             ta1_x_yzz_xyzzzz_0,  \
                             ta1_x_yzz_xyzzzz_1,  \
                             ta1_x_yzz_xzzzz_0,   \
                             ta1_x_yzz_xzzzz_1,   \
                             ta1_x_yzz_xzzzzz_0,  \
                             ta1_x_yzz_xzzzzz_1,  \
                             ta1_x_yzz_yyyyyz_0,  \
                             ta1_x_yzz_yyyyyz_1,  \
                             ta1_x_yzz_yyyyz_0,   \
                             ta1_x_yzz_yyyyz_1,   \
                             ta1_x_yzz_yyyyzz_0,  \
                             ta1_x_yzz_yyyyzz_1,  \
                             ta1_x_yzz_yyyzz_0,   \
                             ta1_x_yzz_yyyzz_1,   \
                             ta1_x_yzz_yyyzzz_0,  \
                             ta1_x_yzz_yyyzzz_1,  \
                             ta1_x_yzz_yyzzz_0,   \
                             ta1_x_yzz_yyzzz_1,   \
                             ta1_x_yzz_yyzzzz_0,  \
                             ta1_x_yzz_yyzzzz_1,  \
                             ta1_x_yzz_yzzzz_0,   \
                             ta1_x_yzz_yzzzz_1,   \
                             ta1_x_yzz_yzzzzz_0,  \
                             ta1_x_yzz_yzzzzz_1,  \
                             ta1_x_yzz_zzzzz_0,   \
                             ta1_x_yzz_zzzzz_1,   \
                             ta1_x_yzz_zzzzzz_0,  \
                             ta1_x_yzz_zzzzzz_1,  \
                             ta1_x_zz_xxxxxx_0,   \
                             ta1_x_zz_xxxxxx_1,   \
                             ta1_x_zz_xxxxxz_0,   \
                             ta1_x_zz_xxxxxz_1,   \
                             ta1_x_zz_xxxxyz_0,   \
                             ta1_x_zz_xxxxyz_1,   \
                             ta1_x_zz_xxxxzz_0,   \
                             ta1_x_zz_xxxxzz_1,   \
                             ta1_x_zz_xxxyyz_0,   \
                             ta1_x_zz_xxxyyz_1,   \
                             ta1_x_zz_xxxyzz_0,   \
                             ta1_x_zz_xxxyzz_1,   \
                             ta1_x_zz_xxxzzz_0,   \
                             ta1_x_zz_xxxzzz_1,   \
                             ta1_x_zz_xxyyyz_0,   \
                             ta1_x_zz_xxyyyz_1,   \
                             ta1_x_zz_xxyyzz_0,   \
                             ta1_x_zz_xxyyzz_1,   \
                             ta1_x_zz_xxyzzz_0,   \
                             ta1_x_zz_xxyzzz_1,   \
                             ta1_x_zz_xxzzzz_0,   \
                             ta1_x_zz_xxzzzz_1,   \
                             ta1_x_zz_xyyyyz_0,   \
                             ta1_x_zz_xyyyyz_1,   \
                             ta1_x_zz_xyyyzz_0,   \
                             ta1_x_zz_xyyyzz_1,   \
                             ta1_x_zz_xyyzzz_0,   \
                             ta1_x_zz_xyyzzz_1,   \
                             ta1_x_zz_xyzzzz_0,   \
                             ta1_x_zz_xyzzzz_1,   \
                             ta1_x_zz_xzzzzz_0,   \
                             ta1_x_zz_xzzzzz_1,   \
                             ta1_x_zz_yyyyyz_0,   \
                             ta1_x_zz_yyyyyz_1,   \
                             ta1_x_zz_yyyyzz_0,   \
                             ta1_x_zz_yyyyzz_1,   \
                             ta1_x_zz_yyyzzz_0,   \
                             ta1_x_zz_yyyzzz_1,   \
                             ta1_x_zz_yyzzzz_0,   \
                             ta1_x_zz_yyzzzz_1,   \
                             ta1_x_zz_yzzzzz_0,   \
                             ta1_x_zz_yzzzzz_1,   \
                             ta1_x_zz_zzzzzz_0,   \
                             ta1_x_zz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzz_xxxxxx_0[i] =
            ta1_x_zz_xxxxxx_0[i] * fe_0 - ta1_x_zz_xxxxxx_1[i] * fe_0 + ta1_x_yzz_xxxxxx_0[i] * pa_y[i] - ta1_x_yzz_xxxxxx_1[i] * pc_y[i];

        ta1_x_yyzz_xxxxxy_0[i] =
            ta1_x_yy_xxxxxy_0[i] * fe_0 - ta1_x_yy_xxxxxy_1[i] * fe_0 + ta1_x_yyz_xxxxxy_0[i] * pa_z[i] - ta1_x_yyz_xxxxxy_1[i] * pc_z[i];

        ta1_x_yyzz_xxxxxz_0[i] =
            ta1_x_zz_xxxxxz_0[i] * fe_0 - ta1_x_zz_xxxxxz_1[i] * fe_0 + ta1_x_yzz_xxxxxz_0[i] * pa_y[i] - ta1_x_yzz_xxxxxz_1[i] * pc_y[i];

        ta1_x_yyzz_xxxxyy_0[i] =
            ta1_x_yy_xxxxyy_0[i] * fe_0 - ta1_x_yy_xxxxyy_1[i] * fe_0 + ta1_x_yyz_xxxxyy_0[i] * pa_z[i] - ta1_x_yyz_xxxxyy_1[i] * pc_z[i];

        ta1_x_yyzz_xxxxyz_0[i] = ta1_x_zz_xxxxyz_0[i] * fe_0 - ta1_x_zz_xxxxyz_1[i] * fe_0 + ta1_x_yzz_xxxxz_0[i] * fe_0 -
                                 ta1_x_yzz_xxxxz_1[i] * fe_0 + ta1_x_yzz_xxxxyz_0[i] * pa_y[i] - ta1_x_yzz_xxxxyz_1[i] * pc_y[i];

        ta1_x_yyzz_xxxxzz_0[i] =
            ta1_x_zz_xxxxzz_0[i] * fe_0 - ta1_x_zz_xxxxzz_1[i] * fe_0 + ta1_x_yzz_xxxxzz_0[i] * pa_y[i] - ta1_x_yzz_xxxxzz_1[i] * pc_y[i];

        ta1_x_yyzz_xxxyyy_0[i] =
            ta1_x_yy_xxxyyy_0[i] * fe_0 - ta1_x_yy_xxxyyy_1[i] * fe_0 + ta1_x_yyz_xxxyyy_0[i] * pa_z[i] - ta1_x_yyz_xxxyyy_1[i] * pc_z[i];

        ta1_x_yyzz_xxxyyz_0[i] = ta1_x_zz_xxxyyz_0[i] * fe_0 - ta1_x_zz_xxxyyz_1[i] * fe_0 + 2.0 * ta1_x_yzz_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_x_yzz_xxxyz_1[i] * fe_0 + ta1_x_yzz_xxxyyz_0[i] * pa_y[i] - ta1_x_yzz_xxxyyz_1[i] * pc_y[i];

        ta1_x_yyzz_xxxyzz_0[i] = ta1_x_zz_xxxyzz_0[i] * fe_0 - ta1_x_zz_xxxyzz_1[i] * fe_0 + ta1_x_yzz_xxxzz_0[i] * fe_0 -
                                 ta1_x_yzz_xxxzz_1[i] * fe_0 + ta1_x_yzz_xxxyzz_0[i] * pa_y[i] - ta1_x_yzz_xxxyzz_1[i] * pc_y[i];

        ta1_x_yyzz_xxxzzz_0[i] =
            ta1_x_zz_xxxzzz_0[i] * fe_0 - ta1_x_zz_xxxzzz_1[i] * fe_0 + ta1_x_yzz_xxxzzz_0[i] * pa_y[i] - ta1_x_yzz_xxxzzz_1[i] * pc_y[i];

        ta1_x_yyzz_xxyyyy_0[i] =
            ta1_x_yy_xxyyyy_0[i] * fe_0 - ta1_x_yy_xxyyyy_1[i] * fe_0 + ta1_x_yyz_xxyyyy_0[i] * pa_z[i] - ta1_x_yyz_xxyyyy_1[i] * pc_z[i];

        ta1_x_yyzz_xxyyyz_0[i] = ta1_x_zz_xxyyyz_0[i] * fe_0 - ta1_x_zz_xxyyyz_1[i] * fe_0 + 3.0 * ta1_x_yzz_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_x_yzz_xxyyz_1[i] * fe_0 + ta1_x_yzz_xxyyyz_0[i] * pa_y[i] - ta1_x_yzz_xxyyyz_1[i] * pc_y[i];

        ta1_x_yyzz_xxyyzz_0[i] = ta1_x_zz_xxyyzz_0[i] * fe_0 - ta1_x_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_yzz_xxyzz_0[i] * fe_0 -
                                 2.0 * ta1_x_yzz_xxyzz_1[i] * fe_0 + ta1_x_yzz_xxyyzz_0[i] * pa_y[i] - ta1_x_yzz_xxyyzz_1[i] * pc_y[i];

        ta1_x_yyzz_xxyzzz_0[i] = ta1_x_zz_xxyzzz_0[i] * fe_0 - ta1_x_zz_xxyzzz_1[i] * fe_0 + ta1_x_yzz_xxzzz_0[i] * fe_0 -
                                 ta1_x_yzz_xxzzz_1[i] * fe_0 + ta1_x_yzz_xxyzzz_0[i] * pa_y[i] - ta1_x_yzz_xxyzzz_1[i] * pc_y[i];

        ta1_x_yyzz_xxzzzz_0[i] =
            ta1_x_zz_xxzzzz_0[i] * fe_0 - ta1_x_zz_xxzzzz_1[i] * fe_0 + ta1_x_yzz_xxzzzz_0[i] * pa_y[i] - ta1_x_yzz_xxzzzz_1[i] * pc_y[i];

        ta1_x_yyzz_xyyyyy_0[i] =
            ta1_x_yy_xyyyyy_0[i] * fe_0 - ta1_x_yy_xyyyyy_1[i] * fe_0 + ta1_x_yyz_xyyyyy_0[i] * pa_z[i] - ta1_x_yyz_xyyyyy_1[i] * pc_z[i];

        ta1_x_yyzz_xyyyyz_0[i] = ta1_x_zz_xyyyyz_0[i] * fe_0 - ta1_x_zz_xyyyyz_1[i] * fe_0 + 4.0 * ta1_x_yzz_xyyyz_0[i] * fe_0 -
                                 4.0 * ta1_x_yzz_xyyyz_1[i] * fe_0 + ta1_x_yzz_xyyyyz_0[i] * pa_y[i] - ta1_x_yzz_xyyyyz_1[i] * pc_y[i];

        ta1_x_yyzz_xyyyzz_0[i] = ta1_x_zz_xyyyzz_0[i] * fe_0 - ta1_x_zz_xyyyzz_1[i] * fe_0 + 3.0 * ta1_x_yzz_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_yzz_xyyzz_1[i] * fe_0 + ta1_x_yzz_xyyyzz_0[i] * pa_y[i] - ta1_x_yzz_xyyyzz_1[i] * pc_y[i];

        ta1_x_yyzz_xyyzzz_0[i] = ta1_x_zz_xyyzzz_0[i] * fe_0 - ta1_x_zz_xyyzzz_1[i] * fe_0 + 2.0 * ta1_x_yzz_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_x_yzz_xyzzz_1[i] * fe_0 + ta1_x_yzz_xyyzzz_0[i] * pa_y[i] - ta1_x_yzz_xyyzzz_1[i] * pc_y[i];

        ta1_x_yyzz_xyzzzz_0[i] = ta1_x_zz_xyzzzz_0[i] * fe_0 - ta1_x_zz_xyzzzz_1[i] * fe_0 + ta1_x_yzz_xzzzz_0[i] * fe_0 -
                                 ta1_x_yzz_xzzzz_1[i] * fe_0 + ta1_x_yzz_xyzzzz_0[i] * pa_y[i] - ta1_x_yzz_xyzzzz_1[i] * pc_y[i];

        ta1_x_yyzz_xzzzzz_0[i] =
            ta1_x_zz_xzzzzz_0[i] * fe_0 - ta1_x_zz_xzzzzz_1[i] * fe_0 + ta1_x_yzz_xzzzzz_0[i] * pa_y[i] - ta1_x_yzz_xzzzzz_1[i] * pc_y[i];

        ta1_x_yyzz_yyyyyy_0[i] =
            ta1_x_yy_yyyyyy_0[i] * fe_0 - ta1_x_yy_yyyyyy_1[i] * fe_0 + ta1_x_yyz_yyyyyy_0[i] * pa_z[i] - ta1_x_yyz_yyyyyy_1[i] * pc_z[i];

        ta1_x_yyzz_yyyyyz_0[i] = ta1_x_zz_yyyyyz_0[i] * fe_0 - ta1_x_zz_yyyyyz_1[i] * fe_0 + 5.0 * ta1_x_yzz_yyyyz_0[i] * fe_0 -
                                 5.0 * ta1_x_yzz_yyyyz_1[i] * fe_0 + ta1_x_yzz_yyyyyz_0[i] * pa_y[i] - ta1_x_yzz_yyyyyz_1[i] * pc_y[i];

        ta1_x_yyzz_yyyyzz_0[i] = ta1_x_zz_yyyyzz_0[i] * fe_0 - ta1_x_zz_yyyyzz_1[i] * fe_0 + 4.0 * ta1_x_yzz_yyyzz_0[i] * fe_0 -
                                 4.0 * ta1_x_yzz_yyyzz_1[i] * fe_0 + ta1_x_yzz_yyyyzz_0[i] * pa_y[i] - ta1_x_yzz_yyyyzz_1[i] * pc_y[i];

        ta1_x_yyzz_yyyzzz_0[i] = ta1_x_zz_yyyzzz_0[i] * fe_0 - ta1_x_zz_yyyzzz_1[i] * fe_0 + 3.0 * ta1_x_yzz_yyzzz_0[i] * fe_0 -
                                 3.0 * ta1_x_yzz_yyzzz_1[i] * fe_0 + ta1_x_yzz_yyyzzz_0[i] * pa_y[i] - ta1_x_yzz_yyyzzz_1[i] * pc_y[i];

        ta1_x_yyzz_yyzzzz_0[i] = ta1_x_zz_yyzzzz_0[i] * fe_0 - ta1_x_zz_yyzzzz_1[i] * fe_0 + 2.0 * ta1_x_yzz_yzzzz_0[i] * fe_0 -
                                 2.0 * ta1_x_yzz_yzzzz_1[i] * fe_0 + ta1_x_yzz_yyzzzz_0[i] * pa_y[i] - ta1_x_yzz_yyzzzz_1[i] * pc_y[i];

        ta1_x_yyzz_yzzzzz_0[i] = ta1_x_zz_yzzzzz_0[i] * fe_0 - ta1_x_zz_yzzzzz_1[i] * fe_0 + ta1_x_yzz_zzzzz_0[i] * fe_0 -
                                 ta1_x_yzz_zzzzz_1[i] * fe_0 + ta1_x_yzz_yzzzzz_0[i] * pa_y[i] - ta1_x_yzz_yzzzzz_1[i] * pc_y[i];

        ta1_x_yyzz_zzzzzz_0[i] =
            ta1_x_zz_zzzzzz_0[i] * fe_0 - ta1_x_zz_zzzzzz_1[i] * fe_0 + ta1_x_yzz_zzzzzz_0[i] * pa_y[i] - ta1_x_yzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 364-392 components of targeted buffer : GI

    auto ta1_x_yzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 364);

    auto ta1_x_yzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 365);

    auto ta1_x_yzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 366);

    auto ta1_x_yzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 367);

    auto ta1_x_yzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 368);

    auto ta1_x_yzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 369);

    auto ta1_x_yzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 370);

    auto ta1_x_yzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 371);

    auto ta1_x_yzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 372);

    auto ta1_x_yzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 373);

    auto ta1_x_yzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 374);

    auto ta1_x_yzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 375);

    auto ta1_x_yzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 376);

    auto ta1_x_yzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 377);

    auto ta1_x_yzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 378);

    auto ta1_x_yzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 379);

    auto ta1_x_yzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 380);

    auto ta1_x_yzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 381);

    auto ta1_x_yzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 382);

    auto ta1_x_yzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 383);

    auto ta1_x_yzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 384);

    auto ta1_x_yzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 385);

    auto ta1_x_yzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 386);

    auto ta1_x_yzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 387);

    auto ta1_x_yzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 388);

    auto ta1_x_yzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 389);

    auto ta1_x_yzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 390);

    auto ta1_x_yzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 391);

#pragma omp simd aligned(pa_y,                    \
                             pc_y,                \
                             ta1_x_yzzz_xxxxxx_0, \
                             ta1_x_yzzz_xxxxxy_0, \
                             ta1_x_yzzz_xxxxxz_0, \
                             ta1_x_yzzz_xxxxyy_0, \
                             ta1_x_yzzz_xxxxyz_0, \
                             ta1_x_yzzz_xxxxzz_0, \
                             ta1_x_yzzz_xxxyyy_0, \
                             ta1_x_yzzz_xxxyyz_0, \
                             ta1_x_yzzz_xxxyzz_0, \
                             ta1_x_yzzz_xxxzzz_0, \
                             ta1_x_yzzz_xxyyyy_0, \
                             ta1_x_yzzz_xxyyyz_0, \
                             ta1_x_yzzz_xxyyzz_0, \
                             ta1_x_yzzz_xxyzzz_0, \
                             ta1_x_yzzz_xxzzzz_0, \
                             ta1_x_yzzz_xyyyyy_0, \
                             ta1_x_yzzz_xyyyyz_0, \
                             ta1_x_yzzz_xyyyzz_0, \
                             ta1_x_yzzz_xyyzzz_0, \
                             ta1_x_yzzz_xyzzzz_0, \
                             ta1_x_yzzz_xzzzzz_0, \
                             ta1_x_yzzz_yyyyyy_0, \
                             ta1_x_yzzz_yyyyyz_0, \
                             ta1_x_yzzz_yyyyzz_0, \
                             ta1_x_yzzz_yyyzzz_0, \
                             ta1_x_yzzz_yyzzzz_0, \
                             ta1_x_yzzz_yzzzzz_0, \
                             ta1_x_yzzz_zzzzzz_0, \
                             ta1_x_zzz_xxxxx_0,   \
                             ta1_x_zzz_xxxxx_1,   \
                             ta1_x_zzz_xxxxxx_0,  \
                             ta1_x_zzz_xxxxxx_1,  \
                             ta1_x_zzz_xxxxxy_0,  \
                             ta1_x_zzz_xxxxxy_1,  \
                             ta1_x_zzz_xxxxxz_0,  \
                             ta1_x_zzz_xxxxxz_1,  \
                             ta1_x_zzz_xxxxy_0,   \
                             ta1_x_zzz_xxxxy_1,   \
                             ta1_x_zzz_xxxxyy_0,  \
                             ta1_x_zzz_xxxxyy_1,  \
                             ta1_x_zzz_xxxxyz_0,  \
                             ta1_x_zzz_xxxxyz_1,  \
                             ta1_x_zzz_xxxxz_0,   \
                             ta1_x_zzz_xxxxz_1,   \
                             ta1_x_zzz_xxxxzz_0,  \
                             ta1_x_zzz_xxxxzz_1,  \
                             ta1_x_zzz_xxxyy_0,   \
                             ta1_x_zzz_xxxyy_1,   \
                             ta1_x_zzz_xxxyyy_0,  \
                             ta1_x_zzz_xxxyyy_1,  \
                             ta1_x_zzz_xxxyyz_0,  \
                             ta1_x_zzz_xxxyyz_1,  \
                             ta1_x_zzz_xxxyz_0,   \
                             ta1_x_zzz_xxxyz_1,   \
                             ta1_x_zzz_xxxyzz_0,  \
                             ta1_x_zzz_xxxyzz_1,  \
                             ta1_x_zzz_xxxzz_0,   \
                             ta1_x_zzz_xxxzz_1,   \
                             ta1_x_zzz_xxxzzz_0,  \
                             ta1_x_zzz_xxxzzz_1,  \
                             ta1_x_zzz_xxyyy_0,   \
                             ta1_x_zzz_xxyyy_1,   \
                             ta1_x_zzz_xxyyyy_0,  \
                             ta1_x_zzz_xxyyyy_1,  \
                             ta1_x_zzz_xxyyyz_0,  \
                             ta1_x_zzz_xxyyyz_1,  \
                             ta1_x_zzz_xxyyz_0,   \
                             ta1_x_zzz_xxyyz_1,   \
                             ta1_x_zzz_xxyyzz_0,  \
                             ta1_x_zzz_xxyyzz_1,  \
                             ta1_x_zzz_xxyzz_0,   \
                             ta1_x_zzz_xxyzz_1,   \
                             ta1_x_zzz_xxyzzz_0,  \
                             ta1_x_zzz_xxyzzz_1,  \
                             ta1_x_zzz_xxzzz_0,   \
                             ta1_x_zzz_xxzzz_1,   \
                             ta1_x_zzz_xxzzzz_0,  \
                             ta1_x_zzz_xxzzzz_1,  \
                             ta1_x_zzz_xyyyy_0,   \
                             ta1_x_zzz_xyyyy_1,   \
                             ta1_x_zzz_xyyyyy_0,  \
                             ta1_x_zzz_xyyyyy_1,  \
                             ta1_x_zzz_xyyyyz_0,  \
                             ta1_x_zzz_xyyyyz_1,  \
                             ta1_x_zzz_xyyyz_0,   \
                             ta1_x_zzz_xyyyz_1,   \
                             ta1_x_zzz_xyyyzz_0,  \
                             ta1_x_zzz_xyyyzz_1,  \
                             ta1_x_zzz_xyyzz_0,   \
                             ta1_x_zzz_xyyzz_1,   \
                             ta1_x_zzz_xyyzzz_0,  \
                             ta1_x_zzz_xyyzzz_1,  \
                             ta1_x_zzz_xyzzz_0,   \
                             ta1_x_zzz_xyzzz_1,   \
                             ta1_x_zzz_xyzzzz_0,  \
                             ta1_x_zzz_xyzzzz_1,  \
                             ta1_x_zzz_xzzzz_0,   \
                             ta1_x_zzz_xzzzz_1,   \
                             ta1_x_zzz_xzzzzz_0,  \
                             ta1_x_zzz_xzzzzz_1,  \
                             ta1_x_zzz_yyyyy_0,   \
                             ta1_x_zzz_yyyyy_1,   \
                             ta1_x_zzz_yyyyyy_0,  \
                             ta1_x_zzz_yyyyyy_1,  \
                             ta1_x_zzz_yyyyyz_0,  \
                             ta1_x_zzz_yyyyyz_1,  \
                             ta1_x_zzz_yyyyz_0,   \
                             ta1_x_zzz_yyyyz_1,   \
                             ta1_x_zzz_yyyyzz_0,  \
                             ta1_x_zzz_yyyyzz_1,  \
                             ta1_x_zzz_yyyzz_0,   \
                             ta1_x_zzz_yyyzz_1,   \
                             ta1_x_zzz_yyyzzz_0,  \
                             ta1_x_zzz_yyyzzz_1,  \
                             ta1_x_zzz_yyzzz_0,   \
                             ta1_x_zzz_yyzzz_1,   \
                             ta1_x_zzz_yyzzzz_0,  \
                             ta1_x_zzz_yyzzzz_1,  \
                             ta1_x_zzz_yzzzz_0,   \
                             ta1_x_zzz_yzzzz_1,   \
                             ta1_x_zzz_yzzzzz_0,  \
                             ta1_x_zzz_yzzzzz_1,  \
                             ta1_x_zzz_zzzzz_0,   \
                             ta1_x_zzz_zzzzz_1,   \
                             ta1_x_zzz_zzzzzz_0,  \
                             ta1_x_zzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzz_xxxxxx_0[i] = ta1_x_zzz_xxxxxx_0[i] * pa_y[i] - ta1_x_zzz_xxxxxx_1[i] * pc_y[i];

        ta1_x_yzzz_xxxxxy_0[i] =
            ta1_x_zzz_xxxxx_0[i] * fe_0 - ta1_x_zzz_xxxxx_1[i] * fe_0 + ta1_x_zzz_xxxxxy_0[i] * pa_y[i] - ta1_x_zzz_xxxxxy_1[i] * pc_y[i];

        ta1_x_yzzz_xxxxxz_0[i] = ta1_x_zzz_xxxxxz_0[i] * pa_y[i] - ta1_x_zzz_xxxxxz_1[i] * pc_y[i];

        ta1_x_yzzz_xxxxyy_0[i] =
            2.0 * ta1_x_zzz_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_zzz_xxxxy_1[i] * fe_0 + ta1_x_zzz_xxxxyy_0[i] * pa_y[i] - ta1_x_zzz_xxxxyy_1[i] * pc_y[i];

        ta1_x_yzzz_xxxxyz_0[i] =
            ta1_x_zzz_xxxxz_0[i] * fe_0 - ta1_x_zzz_xxxxz_1[i] * fe_0 + ta1_x_zzz_xxxxyz_0[i] * pa_y[i] - ta1_x_zzz_xxxxyz_1[i] * pc_y[i];

        ta1_x_yzzz_xxxxzz_0[i] = ta1_x_zzz_xxxxzz_0[i] * pa_y[i] - ta1_x_zzz_xxxxzz_1[i] * pc_y[i];

        ta1_x_yzzz_xxxyyy_0[i] =
            3.0 * ta1_x_zzz_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxxyy_1[i] * fe_0 + ta1_x_zzz_xxxyyy_0[i] * pa_y[i] - ta1_x_zzz_xxxyyy_1[i] * pc_y[i];

        ta1_x_yzzz_xxxyyz_0[i] =
            2.0 * ta1_x_zzz_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xxxyz_1[i] * fe_0 + ta1_x_zzz_xxxyyz_0[i] * pa_y[i] - ta1_x_zzz_xxxyyz_1[i] * pc_y[i];

        ta1_x_yzzz_xxxyzz_0[i] =
            ta1_x_zzz_xxxzz_0[i] * fe_0 - ta1_x_zzz_xxxzz_1[i] * fe_0 + ta1_x_zzz_xxxyzz_0[i] * pa_y[i] - ta1_x_zzz_xxxyzz_1[i] * pc_y[i];

        ta1_x_yzzz_xxxzzz_0[i] = ta1_x_zzz_xxxzzz_0[i] * pa_y[i] - ta1_x_zzz_xxxzzz_1[i] * pc_y[i];

        ta1_x_yzzz_xxyyyy_0[i] =
            4.0 * ta1_x_zzz_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxyyy_1[i] * fe_0 + ta1_x_zzz_xxyyyy_0[i] * pa_y[i] - ta1_x_zzz_xxyyyy_1[i] * pc_y[i];

        ta1_x_yzzz_xxyyyz_0[i] =
            3.0 * ta1_x_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxyyz_1[i] * fe_0 + ta1_x_zzz_xxyyyz_0[i] * pa_y[i] - ta1_x_zzz_xxyyyz_1[i] * pc_y[i];

        ta1_x_yzzz_xxyyzz_0[i] =
            2.0 * ta1_x_zzz_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xxyzz_1[i] * fe_0 + ta1_x_zzz_xxyyzz_0[i] * pa_y[i] - ta1_x_zzz_xxyyzz_1[i] * pc_y[i];

        ta1_x_yzzz_xxyzzz_0[i] =
            ta1_x_zzz_xxzzz_0[i] * fe_0 - ta1_x_zzz_xxzzz_1[i] * fe_0 + ta1_x_zzz_xxyzzz_0[i] * pa_y[i] - ta1_x_zzz_xxyzzz_1[i] * pc_y[i];

        ta1_x_yzzz_xxzzzz_0[i] = ta1_x_zzz_xxzzzz_0[i] * pa_y[i] - ta1_x_zzz_xxzzzz_1[i] * pc_y[i];

        ta1_x_yzzz_xyyyyy_0[i] =
            5.0 * ta1_x_zzz_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_zzz_xyyyy_1[i] * fe_0 + ta1_x_zzz_xyyyyy_0[i] * pa_y[i] - ta1_x_zzz_xyyyyy_1[i] * pc_y[i];

        ta1_x_yzzz_xyyyyz_0[i] =
            4.0 * ta1_x_zzz_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xyyyz_1[i] * fe_0 + ta1_x_zzz_xyyyyz_0[i] * pa_y[i] - ta1_x_zzz_xyyyyz_1[i] * pc_y[i];

        ta1_x_yzzz_xyyyzz_0[i] =
            3.0 * ta1_x_zzz_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xyyzz_1[i] * fe_0 + ta1_x_zzz_xyyyzz_0[i] * pa_y[i] - ta1_x_zzz_xyyyzz_1[i] * pc_y[i];

        ta1_x_yzzz_xyyzzz_0[i] =
            2.0 * ta1_x_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyzzz_1[i] * fe_0 + ta1_x_zzz_xyyzzz_0[i] * pa_y[i] - ta1_x_zzz_xyyzzz_1[i] * pc_y[i];

        ta1_x_yzzz_xyzzzz_0[i] =
            ta1_x_zzz_xzzzz_0[i] * fe_0 - ta1_x_zzz_xzzzz_1[i] * fe_0 + ta1_x_zzz_xyzzzz_0[i] * pa_y[i] - ta1_x_zzz_xyzzzz_1[i] * pc_y[i];

        ta1_x_yzzz_xzzzzz_0[i] = ta1_x_zzz_xzzzzz_0[i] * pa_y[i] - ta1_x_zzz_xzzzzz_1[i] * pc_y[i];

        ta1_x_yzzz_yyyyyy_0[i] =
            6.0 * ta1_x_zzz_yyyyy_0[i] * fe_0 - 6.0 * ta1_x_zzz_yyyyy_1[i] * fe_0 + ta1_x_zzz_yyyyyy_0[i] * pa_y[i] - ta1_x_zzz_yyyyyy_1[i] * pc_y[i];

        ta1_x_yzzz_yyyyyz_0[i] =
            5.0 * ta1_x_zzz_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_zzz_yyyyz_1[i] * fe_0 + ta1_x_zzz_yyyyyz_0[i] * pa_y[i] - ta1_x_zzz_yyyyyz_1[i] * pc_y[i];

        ta1_x_yzzz_yyyyzz_0[i] =
            4.0 * ta1_x_zzz_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyyzz_1[i] * fe_0 + ta1_x_zzz_yyyyzz_0[i] * pa_y[i] - ta1_x_zzz_yyyyzz_1[i] * pc_y[i];

        ta1_x_yzzz_yyyzzz_0[i] =
            3.0 * ta1_x_zzz_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_zzz_yyzzz_1[i] * fe_0 + ta1_x_zzz_yyyzzz_0[i] * pa_y[i] - ta1_x_zzz_yyyzzz_1[i] * pc_y[i];

        ta1_x_yzzz_yyzzzz_0[i] =
            2.0 * ta1_x_zzz_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_yzzzz_1[i] * fe_0 + ta1_x_zzz_yyzzzz_0[i] * pa_y[i] - ta1_x_zzz_yyzzzz_1[i] * pc_y[i];

        ta1_x_yzzz_yzzzzz_0[i] =
            ta1_x_zzz_zzzzz_0[i] * fe_0 - ta1_x_zzz_zzzzz_1[i] * fe_0 + ta1_x_zzz_yzzzzz_0[i] * pa_y[i] - ta1_x_zzz_yzzzzz_1[i] * pc_y[i];

        ta1_x_yzzz_zzzzzz_0[i] = ta1_x_zzz_zzzzzz_0[i] * pa_y[i] - ta1_x_zzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 392-420 components of targeted buffer : GI

    auto ta1_x_zzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 392);

    auto ta1_x_zzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 393);

    auto ta1_x_zzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 394);

    auto ta1_x_zzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 395);

    auto ta1_x_zzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 396);

    auto ta1_x_zzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 397);

    auto ta1_x_zzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 398);

    auto ta1_x_zzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 399);

    auto ta1_x_zzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 400);

    auto ta1_x_zzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 401);

    auto ta1_x_zzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 402);

    auto ta1_x_zzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 403);

    auto ta1_x_zzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 404);

    auto ta1_x_zzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 405);

    auto ta1_x_zzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 406);

    auto ta1_x_zzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 407);

    auto ta1_x_zzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 408);

    auto ta1_x_zzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 409);

    auto ta1_x_zzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 410);

    auto ta1_x_zzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 411);

    auto ta1_x_zzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 412);

    auto ta1_x_zzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 413);

    auto ta1_x_zzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 414);

    auto ta1_x_zzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 415);

    auto ta1_x_zzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 416);

    auto ta1_x_zzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 417);

    auto ta1_x_zzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 418);

    auto ta1_x_zzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 419);

#pragma omp simd aligned(pa_z,                    \
                             pc_z,                \
                             ta1_x_zz_xxxxxx_0,   \
                             ta1_x_zz_xxxxxx_1,   \
                             ta1_x_zz_xxxxxy_0,   \
                             ta1_x_zz_xxxxxy_1,   \
                             ta1_x_zz_xxxxxz_0,   \
                             ta1_x_zz_xxxxxz_1,   \
                             ta1_x_zz_xxxxyy_0,   \
                             ta1_x_zz_xxxxyy_1,   \
                             ta1_x_zz_xxxxyz_0,   \
                             ta1_x_zz_xxxxyz_1,   \
                             ta1_x_zz_xxxxzz_0,   \
                             ta1_x_zz_xxxxzz_1,   \
                             ta1_x_zz_xxxyyy_0,   \
                             ta1_x_zz_xxxyyy_1,   \
                             ta1_x_zz_xxxyyz_0,   \
                             ta1_x_zz_xxxyyz_1,   \
                             ta1_x_zz_xxxyzz_0,   \
                             ta1_x_zz_xxxyzz_1,   \
                             ta1_x_zz_xxxzzz_0,   \
                             ta1_x_zz_xxxzzz_1,   \
                             ta1_x_zz_xxyyyy_0,   \
                             ta1_x_zz_xxyyyy_1,   \
                             ta1_x_zz_xxyyyz_0,   \
                             ta1_x_zz_xxyyyz_1,   \
                             ta1_x_zz_xxyyzz_0,   \
                             ta1_x_zz_xxyyzz_1,   \
                             ta1_x_zz_xxyzzz_0,   \
                             ta1_x_zz_xxyzzz_1,   \
                             ta1_x_zz_xxzzzz_0,   \
                             ta1_x_zz_xxzzzz_1,   \
                             ta1_x_zz_xyyyyy_0,   \
                             ta1_x_zz_xyyyyy_1,   \
                             ta1_x_zz_xyyyyz_0,   \
                             ta1_x_zz_xyyyyz_1,   \
                             ta1_x_zz_xyyyzz_0,   \
                             ta1_x_zz_xyyyzz_1,   \
                             ta1_x_zz_xyyzzz_0,   \
                             ta1_x_zz_xyyzzz_1,   \
                             ta1_x_zz_xyzzzz_0,   \
                             ta1_x_zz_xyzzzz_1,   \
                             ta1_x_zz_xzzzzz_0,   \
                             ta1_x_zz_xzzzzz_1,   \
                             ta1_x_zz_yyyyyy_0,   \
                             ta1_x_zz_yyyyyy_1,   \
                             ta1_x_zz_yyyyyz_0,   \
                             ta1_x_zz_yyyyyz_1,   \
                             ta1_x_zz_yyyyzz_0,   \
                             ta1_x_zz_yyyyzz_1,   \
                             ta1_x_zz_yyyzzz_0,   \
                             ta1_x_zz_yyyzzz_1,   \
                             ta1_x_zz_yyzzzz_0,   \
                             ta1_x_zz_yyzzzz_1,   \
                             ta1_x_zz_yzzzzz_0,   \
                             ta1_x_zz_yzzzzz_1,   \
                             ta1_x_zz_zzzzzz_0,   \
                             ta1_x_zz_zzzzzz_1,   \
                             ta1_x_zzz_xxxxx_0,   \
                             ta1_x_zzz_xxxxx_1,   \
                             ta1_x_zzz_xxxxxx_0,  \
                             ta1_x_zzz_xxxxxx_1,  \
                             ta1_x_zzz_xxxxxy_0,  \
                             ta1_x_zzz_xxxxxy_1,  \
                             ta1_x_zzz_xxxxxz_0,  \
                             ta1_x_zzz_xxxxxz_1,  \
                             ta1_x_zzz_xxxxy_0,   \
                             ta1_x_zzz_xxxxy_1,   \
                             ta1_x_zzz_xxxxyy_0,  \
                             ta1_x_zzz_xxxxyy_1,  \
                             ta1_x_zzz_xxxxyz_0,  \
                             ta1_x_zzz_xxxxyz_1,  \
                             ta1_x_zzz_xxxxz_0,   \
                             ta1_x_zzz_xxxxz_1,   \
                             ta1_x_zzz_xxxxzz_0,  \
                             ta1_x_zzz_xxxxzz_1,  \
                             ta1_x_zzz_xxxyy_0,   \
                             ta1_x_zzz_xxxyy_1,   \
                             ta1_x_zzz_xxxyyy_0,  \
                             ta1_x_zzz_xxxyyy_1,  \
                             ta1_x_zzz_xxxyyz_0,  \
                             ta1_x_zzz_xxxyyz_1,  \
                             ta1_x_zzz_xxxyz_0,   \
                             ta1_x_zzz_xxxyz_1,   \
                             ta1_x_zzz_xxxyzz_0,  \
                             ta1_x_zzz_xxxyzz_1,  \
                             ta1_x_zzz_xxxzz_0,   \
                             ta1_x_zzz_xxxzz_1,   \
                             ta1_x_zzz_xxxzzz_0,  \
                             ta1_x_zzz_xxxzzz_1,  \
                             ta1_x_zzz_xxyyy_0,   \
                             ta1_x_zzz_xxyyy_1,   \
                             ta1_x_zzz_xxyyyy_0,  \
                             ta1_x_zzz_xxyyyy_1,  \
                             ta1_x_zzz_xxyyyz_0,  \
                             ta1_x_zzz_xxyyyz_1,  \
                             ta1_x_zzz_xxyyz_0,   \
                             ta1_x_zzz_xxyyz_1,   \
                             ta1_x_zzz_xxyyzz_0,  \
                             ta1_x_zzz_xxyyzz_1,  \
                             ta1_x_zzz_xxyzz_0,   \
                             ta1_x_zzz_xxyzz_1,   \
                             ta1_x_zzz_xxyzzz_0,  \
                             ta1_x_zzz_xxyzzz_1,  \
                             ta1_x_zzz_xxzzz_0,   \
                             ta1_x_zzz_xxzzz_1,   \
                             ta1_x_zzz_xxzzzz_0,  \
                             ta1_x_zzz_xxzzzz_1,  \
                             ta1_x_zzz_xyyyy_0,   \
                             ta1_x_zzz_xyyyy_1,   \
                             ta1_x_zzz_xyyyyy_0,  \
                             ta1_x_zzz_xyyyyy_1,  \
                             ta1_x_zzz_xyyyyz_0,  \
                             ta1_x_zzz_xyyyyz_1,  \
                             ta1_x_zzz_xyyyz_0,   \
                             ta1_x_zzz_xyyyz_1,   \
                             ta1_x_zzz_xyyyzz_0,  \
                             ta1_x_zzz_xyyyzz_1,  \
                             ta1_x_zzz_xyyzz_0,   \
                             ta1_x_zzz_xyyzz_1,   \
                             ta1_x_zzz_xyyzzz_0,  \
                             ta1_x_zzz_xyyzzz_1,  \
                             ta1_x_zzz_xyzzz_0,   \
                             ta1_x_zzz_xyzzz_1,   \
                             ta1_x_zzz_xyzzzz_0,  \
                             ta1_x_zzz_xyzzzz_1,  \
                             ta1_x_zzz_xzzzz_0,   \
                             ta1_x_zzz_xzzzz_1,   \
                             ta1_x_zzz_xzzzzz_0,  \
                             ta1_x_zzz_xzzzzz_1,  \
                             ta1_x_zzz_yyyyy_0,   \
                             ta1_x_zzz_yyyyy_1,   \
                             ta1_x_zzz_yyyyyy_0,  \
                             ta1_x_zzz_yyyyyy_1,  \
                             ta1_x_zzz_yyyyyz_0,  \
                             ta1_x_zzz_yyyyyz_1,  \
                             ta1_x_zzz_yyyyz_0,   \
                             ta1_x_zzz_yyyyz_1,   \
                             ta1_x_zzz_yyyyzz_0,  \
                             ta1_x_zzz_yyyyzz_1,  \
                             ta1_x_zzz_yyyzz_0,   \
                             ta1_x_zzz_yyyzz_1,   \
                             ta1_x_zzz_yyyzzz_0,  \
                             ta1_x_zzz_yyyzzz_1,  \
                             ta1_x_zzz_yyzzz_0,   \
                             ta1_x_zzz_yyzzz_1,   \
                             ta1_x_zzz_yyzzzz_0,  \
                             ta1_x_zzz_yyzzzz_1,  \
                             ta1_x_zzz_yzzzz_0,   \
                             ta1_x_zzz_yzzzz_1,   \
                             ta1_x_zzz_yzzzzz_0,  \
                             ta1_x_zzz_yzzzzz_1,  \
                             ta1_x_zzz_zzzzz_0,   \
                             ta1_x_zzz_zzzzz_1,   \
                             ta1_x_zzz_zzzzzz_0,  \
                             ta1_x_zzz_zzzzzz_1,  \
                             ta1_x_zzzz_xxxxxx_0, \
                             ta1_x_zzzz_xxxxxy_0, \
                             ta1_x_zzzz_xxxxxz_0, \
                             ta1_x_zzzz_xxxxyy_0, \
                             ta1_x_zzzz_xxxxyz_0, \
                             ta1_x_zzzz_xxxxzz_0, \
                             ta1_x_zzzz_xxxyyy_0, \
                             ta1_x_zzzz_xxxyyz_0, \
                             ta1_x_zzzz_xxxyzz_0, \
                             ta1_x_zzzz_xxxzzz_0, \
                             ta1_x_zzzz_xxyyyy_0, \
                             ta1_x_zzzz_xxyyyz_0, \
                             ta1_x_zzzz_xxyyzz_0, \
                             ta1_x_zzzz_xxyzzz_0, \
                             ta1_x_zzzz_xxzzzz_0, \
                             ta1_x_zzzz_xyyyyy_0, \
                             ta1_x_zzzz_xyyyyz_0, \
                             ta1_x_zzzz_xyyyzz_0, \
                             ta1_x_zzzz_xyyzzz_0, \
                             ta1_x_zzzz_xyzzzz_0, \
                             ta1_x_zzzz_xzzzzz_0, \
                             ta1_x_zzzz_yyyyyy_0, \
                             ta1_x_zzzz_yyyyyz_0, \
                             ta1_x_zzzz_yyyyzz_0, \
                             ta1_x_zzzz_yyyzzz_0, \
                             ta1_x_zzzz_yyzzzz_0, \
                             ta1_x_zzzz_yzzzzz_0, \
                             ta1_x_zzzz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzz_xxxxxx_0[i] =
            3.0 * ta1_x_zz_xxxxxx_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxxx_1[i] * fe_0 + ta1_x_zzz_xxxxxx_0[i] * pa_z[i] - ta1_x_zzz_xxxxxx_1[i] * pc_z[i];

        ta1_x_zzzz_xxxxxy_0[i] =
            3.0 * ta1_x_zz_xxxxxy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxxy_1[i] * fe_0 + ta1_x_zzz_xxxxxy_0[i] * pa_z[i] - ta1_x_zzz_xxxxxy_1[i] * pc_z[i];

        ta1_x_zzzz_xxxxxz_0[i] = 3.0 * ta1_x_zz_xxxxxz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxxz_1[i] * fe_0 + ta1_x_zzz_xxxxx_0[i] * fe_0 -
                                 ta1_x_zzz_xxxxx_1[i] * fe_0 + ta1_x_zzz_xxxxxz_0[i] * pa_z[i] - ta1_x_zzz_xxxxxz_1[i] * pc_z[i];

        ta1_x_zzzz_xxxxyy_0[i] =
            3.0 * ta1_x_zz_xxxxyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxyy_1[i] * fe_0 + ta1_x_zzz_xxxxyy_0[i] * pa_z[i] - ta1_x_zzz_xxxxyy_1[i] * pc_z[i];

        ta1_x_zzzz_xxxxyz_0[i] = 3.0 * ta1_x_zz_xxxxyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxyz_1[i] * fe_0 + ta1_x_zzz_xxxxy_0[i] * fe_0 -
                                 ta1_x_zzz_xxxxy_1[i] * fe_0 + ta1_x_zzz_xxxxyz_0[i] * pa_z[i] - ta1_x_zzz_xxxxyz_1[i] * pc_z[i];

        ta1_x_zzzz_xxxxzz_0[i] = 3.0 * ta1_x_zz_xxxxzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxxxz_0[i] * fe_0 -
                                 2.0 * ta1_x_zzz_xxxxz_1[i] * fe_0 + ta1_x_zzz_xxxxzz_0[i] * pa_z[i] - ta1_x_zzz_xxxxzz_1[i] * pc_z[i];

        ta1_x_zzzz_xxxyyy_0[i] =
            3.0 * ta1_x_zz_xxxyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxyyy_1[i] * fe_0 + ta1_x_zzz_xxxyyy_0[i] * pa_z[i] - ta1_x_zzz_xxxyyy_1[i] * pc_z[i];

        ta1_x_zzzz_xxxyyz_0[i] = 3.0 * ta1_x_zz_xxxyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxyyz_1[i] * fe_0 + ta1_x_zzz_xxxyy_0[i] * fe_0 -
                                 ta1_x_zzz_xxxyy_1[i] * fe_0 + ta1_x_zzz_xxxyyz_0[i] * pa_z[i] - ta1_x_zzz_xxxyyz_1[i] * pc_z[i];

        ta1_x_zzzz_xxxyzz_0[i] = 3.0 * ta1_x_zz_xxxyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_x_zzz_xxxyz_1[i] * fe_0 + ta1_x_zzz_xxxyzz_0[i] * pa_z[i] - ta1_x_zzz_xxxyzz_1[i] * pc_z[i];

        ta1_x_zzzz_xxxzzz_0[i] = 3.0 * ta1_x_zz_xxxzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_xxxzz_0[i] * fe_0 -
                                 3.0 * ta1_x_zzz_xxxzz_1[i] * fe_0 + ta1_x_zzz_xxxzzz_0[i] * pa_z[i] - ta1_x_zzz_xxxzzz_1[i] * pc_z[i];

        ta1_x_zzzz_xxyyyy_0[i] =
            3.0 * ta1_x_zz_xxyyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyyyy_1[i] * fe_0 + ta1_x_zzz_xxyyyy_0[i] * pa_z[i] - ta1_x_zzz_xxyyyy_1[i] * pc_z[i];

        ta1_x_zzzz_xxyyyz_0[i] = 3.0 * ta1_x_zz_xxyyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyyyz_1[i] * fe_0 + ta1_x_zzz_xxyyy_0[i] * fe_0 -
                                 ta1_x_zzz_xxyyy_1[i] * fe_0 + ta1_x_zzz_xxyyyz_0[i] * pa_z[i] - ta1_x_zzz_xxyyyz_1[i] * pc_z[i];

        ta1_x_zzzz_xxyyzz_0[i] = 3.0 * ta1_x_zz_xxyyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxyyz_0[i] * fe_0 -
                                 2.0 * ta1_x_zzz_xxyyz_1[i] * fe_0 + ta1_x_zzz_xxyyzz_0[i] * pa_z[i] - ta1_x_zzz_xxyyzz_1[i] * pc_z[i];

        ta1_x_zzzz_xxyzzz_0[i] = 3.0 * ta1_x_zz_xxyzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_zzz_xxyzz_1[i] * fe_0 + ta1_x_zzz_xxyzzz_0[i] * pa_z[i] - ta1_x_zzz_xxyzzz_1[i] * pc_z[i];

        ta1_x_zzzz_xxzzzz_0[i] = 3.0 * ta1_x_zz_xxzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxzzzz_1[i] * fe_0 + 4.0 * ta1_x_zzz_xxzzz_0[i] * fe_0 -
                                 4.0 * ta1_x_zzz_xxzzz_1[i] * fe_0 + ta1_x_zzz_xxzzzz_0[i] * pa_z[i] - ta1_x_zzz_xxzzzz_1[i] * pc_z[i];

        ta1_x_zzzz_xyyyyy_0[i] =
            3.0 * ta1_x_zz_xyyyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyyyy_1[i] * fe_0 + ta1_x_zzz_xyyyyy_0[i] * pa_z[i] - ta1_x_zzz_xyyyyy_1[i] * pc_z[i];

        ta1_x_zzzz_xyyyyz_0[i] = 3.0 * ta1_x_zz_xyyyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyyyz_1[i] * fe_0 + ta1_x_zzz_xyyyy_0[i] * fe_0 -
                                 ta1_x_zzz_xyyyy_1[i] * fe_0 + ta1_x_zzz_xyyyyz_0[i] * pa_z[i] - ta1_x_zzz_xyyyyz_1[i] * pc_z[i];

        ta1_x_zzzz_xyyyzz_0[i] = 3.0 * ta1_x_zz_xyyyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_x_zzz_xyyyz_1[i] * fe_0 + ta1_x_zzz_xyyyzz_0[i] * pa_z[i] - ta1_x_zzz_xyyyzz_1[i] * pc_z[i];

        ta1_x_zzzz_xyyzzz_0[i] = 3.0 * ta1_x_zz_xyyzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_zzz_xyyzz_1[i] * fe_0 + ta1_x_zzz_xyyzzz_0[i] * pa_z[i] - ta1_x_zzz_xyyzzz_1[i] * pc_z[i];

        ta1_x_zzzz_xyzzzz_0[i] = 3.0 * ta1_x_zz_xyzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyzzzz_1[i] * fe_0 + 4.0 * ta1_x_zzz_xyzzz_0[i] * fe_0 -
                                 4.0 * ta1_x_zzz_xyzzz_1[i] * fe_0 + ta1_x_zzz_xyzzzz_0[i] * pa_z[i] - ta1_x_zzz_xyzzzz_1[i] * pc_z[i];

        ta1_x_zzzz_xzzzzz_0[i] = 3.0 * ta1_x_zz_xzzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xzzzzz_1[i] * fe_0 + 5.0 * ta1_x_zzz_xzzzz_0[i] * fe_0 -
                                 5.0 * ta1_x_zzz_xzzzz_1[i] * fe_0 + ta1_x_zzz_xzzzzz_0[i] * pa_z[i] - ta1_x_zzz_xzzzzz_1[i] * pc_z[i];

        ta1_x_zzzz_yyyyyy_0[i] =
            3.0 * ta1_x_zz_yyyyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyyyy_1[i] * fe_0 + ta1_x_zzz_yyyyyy_0[i] * pa_z[i] - ta1_x_zzz_yyyyyy_1[i] * pc_z[i];

        ta1_x_zzzz_yyyyyz_0[i] = 3.0 * ta1_x_zz_yyyyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyyyz_1[i] * fe_0 + ta1_x_zzz_yyyyy_0[i] * fe_0 -
                                 ta1_x_zzz_yyyyy_1[i] * fe_0 + ta1_x_zzz_yyyyyz_0[i] * pa_z[i] - ta1_x_zzz_yyyyyz_1[i] * pc_z[i];

        ta1_x_zzzz_yyyyzz_0[i] = 3.0 * ta1_x_zz_yyyyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_yyyyz_0[i] * fe_0 -
                                 2.0 * ta1_x_zzz_yyyyz_1[i] * fe_0 + ta1_x_zzz_yyyyzz_0[i] * pa_z[i] - ta1_x_zzz_yyyyzz_1[i] * pc_z[i];

        ta1_x_zzzz_yyyzzz_0[i] = 3.0 * ta1_x_zz_yyyzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_yyyzz_0[i] * fe_0 -
                                 3.0 * ta1_x_zzz_yyyzz_1[i] * fe_0 + ta1_x_zzz_yyyzzz_0[i] * pa_z[i] - ta1_x_zzz_yyyzzz_1[i] * pc_z[i];

        ta1_x_zzzz_yyzzzz_0[i] = 3.0 * ta1_x_zz_yyzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyzzzz_1[i] * fe_0 + 4.0 * ta1_x_zzz_yyzzz_0[i] * fe_0 -
                                 4.0 * ta1_x_zzz_yyzzz_1[i] * fe_0 + ta1_x_zzz_yyzzzz_0[i] * pa_z[i] - ta1_x_zzz_yyzzzz_1[i] * pc_z[i];

        ta1_x_zzzz_yzzzzz_0[i] = 3.0 * ta1_x_zz_yzzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yzzzzz_1[i] * fe_0 + 5.0 * ta1_x_zzz_yzzzz_0[i] * fe_0 -
                                 5.0 * ta1_x_zzz_yzzzz_1[i] * fe_0 + ta1_x_zzz_yzzzzz_0[i] * pa_z[i] - ta1_x_zzz_yzzzzz_1[i] * pc_z[i];

        ta1_x_zzzz_zzzzzz_0[i] = 3.0 * ta1_x_zz_zzzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_zzzzzz_1[i] * fe_0 + 6.0 * ta1_x_zzz_zzzzz_0[i] * fe_0 -
                                 6.0 * ta1_x_zzz_zzzzz_1[i] * fe_0 + ta1_x_zzz_zzzzzz_0[i] * pa_z[i] - ta1_x_zzz_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 420-448 components of targeted buffer : GI

    auto ta1_y_xxxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 420);

    auto ta1_y_xxxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 421);

    auto ta1_y_xxxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 422);

    auto ta1_y_xxxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 423);

    auto ta1_y_xxxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 424);

    auto ta1_y_xxxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 425);

    auto ta1_y_xxxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 426);

    auto ta1_y_xxxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 427);

    auto ta1_y_xxxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 428);

    auto ta1_y_xxxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 429);

    auto ta1_y_xxxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 430);

    auto ta1_y_xxxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 431);

    auto ta1_y_xxxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 432);

    auto ta1_y_xxxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 433);

    auto ta1_y_xxxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 434);

    auto ta1_y_xxxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 435);

    auto ta1_y_xxxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 436);

    auto ta1_y_xxxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 437);

    auto ta1_y_xxxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 438);

    auto ta1_y_xxxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 439);

    auto ta1_y_xxxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 440);

    auto ta1_y_xxxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 441);

    auto ta1_y_xxxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 442);

    auto ta1_y_xxxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 443);

    auto ta1_y_xxxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 444);

    auto ta1_y_xxxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 445);

    auto ta1_y_xxxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 446);

    auto ta1_y_xxxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 447);

#pragma omp simd aligned(pa_x,                    \
                             pc_x,                \
                             ta1_y_xx_xxxxxx_0,   \
                             ta1_y_xx_xxxxxx_1,   \
                             ta1_y_xx_xxxxxy_0,   \
                             ta1_y_xx_xxxxxy_1,   \
                             ta1_y_xx_xxxxxz_0,   \
                             ta1_y_xx_xxxxxz_1,   \
                             ta1_y_xx_xxxxyy_0,   \
                             ta1_y_xx_xxxxyy_1,   \
                             ta1_y_xx_xxxxyz_0,   \
                             ta1_y_xx_xxxxyz_1,   \
                             ta1_y_xx_xxxxzz_0,   \
                             ta1_y_xx_xxxxzz_1,   \
                             ta1_y_xx_xxxyyy_0,   \
                             ta1_y_xx_xxxyyy_1,   \
                             ta1_y_xx_xxxyyz_0,   \
                             ta1_y_xx_xxxyyz_1,   \
                             ta1_y_xx_xxxyzz_0,   \
                             ta1_y_xx_xxxyzz_1,   \
                             ta1_y_xx_xxxzzz_0,   \
                             ta1_y_xx_xxxzzz_1,   \
                             ta1_y_xx_xxyyyy_0,   \
                             ta1_y_xx_xxyyyy_1,   \
                             ta1_y_xx_xxyyyz_0,   \
                             ta1_y_xx_xxyyyz_1,   \
                             ta1_y_xx_xxyyzz_0,   \
                             ta1_y_xx_xxyyzz_1,   \
                             ta1_y_xx_xxyzzz_0,   \
                             ta1_y_xx_xxyzzz_1,   \
                             ta1_y_xx_xxzzzz_0,   \
                             ta1_y_xx_xxzzzz_1,   \
                             ta1_y_xx_xyyyyy_0,   \
                             ta1_y_xx_xyyyyy_1,   \
                             ta1_y_xx_xyyyyz_0,   \
                             ta1_y_xx_xyyyyz_1,   \
                             ta1_y_xx_xyyyzz_0,   \
                             ta1_y_xx_xyyyzz_1,   \
                             ta1_y_xx_xyyzzz_0,   \
                             ta1_y_xx_xyyzzz_1,   \
                             ta1_y_xx_xyzzzz_0,   \
                             ta1_y_xx_xyzzzz_1,   \
                             ta1_y_xx_xzzzzz_0,   \
                             ta1_y_xx_xzzzzz_1,   \
                             ta1_y_xx_yyyyyy_0,   \
                             ta1_y_xx_yyyyyy_1,   \
                             ta1_y_xx_yyyyyz_0,   \
                             ta1_y_xx_yyyyyz_1,   \
                             ta1_y_xx_yyyyzz_0,   \
                             ta1_y_xx_yyyyzz_1,   \
                             ta1_y_xx_yyyzzz_0,   \
                             ta1_y_xx_yyyzzz_1,   \
                             ta1_y_xx_yyzzzz_0,   \
                             ta1_y_xx_yyzzzz_1,   \
                             ta1_y_xx_yzzzzz_0,   \
                             ta1_y_xx_yzzzzz_1,   \
                             ta1_y_xx_zzzzzz_0,   \
                             ta1_y_xx_zzzzzz_1,   \
                             ta1_y_xxx_xxxxx_0,   \
                             ta1_y_xxx_xxxxx_1,   \
                             ta1_y_xxx_xxxxxx_0,  \
                             ta1_y_xxx_xxxxxx_1,  \
                             ta1_y_xxx_xxxxxy_0,  \
                             ta1_y_xxx_xxxxxy_1,  \
                             ta1_y_xxx_xxxxxz_0,  \
                             ta1_y_xxx_xxxxxz_1,  \
                             ta1_y_xxx_xxxxy_0,   \
                             ta1_y_xxx_xxxxy_1,   \
                             ta1_y_xxx_xxxxyy_0,  \
                             ta1_y_xxx_xxxxyy_1,  \
                             ta1_y_xxx_xxxxyz_0,  \
                             ta1_y_xxx_xxxxyz_1,  \
                             ta1_y_xxx_xxxxz_0,   \
                             ta1_y_xxx_xxxxz_1,   \
                             ta1_y_xxx_xxxxzz_0,  \
                             ta1_y_xxx_xxxxzz_1,  \
                             ta1_y_xxx_xxxyy_0,   \
                             ta1_y_xxx_xxxyy_1,   \
                             ta1_y_xxx_xxxyyy_0,  \
                             ta1_y_xxx_xxxyyy_1,  \
                             ta1_y_xxx_xxxyyz_0,  \
                             ta1_y_xxx_xxxyyz_1,  \
                             ta1_y_xxx_xxxyz_0,   \
                             ta1_y_xxx_xxxyz_1,   \
                             ta1_y_xxx_xxxyzz_0,  \
                             ta1_y_xxx_xxxyzz_1,  \
                             ta1_y_xxx_xxxzz_0,   \
                             ta1_y_xxx_xxxzz_1,   \
                             ta1_y_xxx_xxxzzz_0,  \
                             ta1_y_xxx_xxxzzz_1,  \
                             ta1_y_xxx_xxyyy_0,   \
                             ta1_y_xxx_xxyyy_1,   \
                             ta1_y_xxx_xxyyyy_0,  \
                             ta1_y_xxx_xxyyyy_1,  \
                             ta1_y_xxx_xxyyyz_0,  \
                             ta1_y_xxx_xxyyyz_1,  \
                             ta1_y_xxx_xxyyz_0,   \
                             ta1_y_xxx_xxyyz_1,   \
                             ta1_y_xxx_xxyyzz_0,  \
                             ta1_y_xxx_xxyyzz_1,  \
                             ta1_y_xxx_xxyzz_0,   \
                             ta1_y_xxx_xxyzz_1,   \
                             ta1_y_xxx_xxyzzz_0,  \
                             ta1_y_xxx_xxyzzz_1,  \
                             ta1_y_xxx_xxzzz_0,   \
                             ta1_y_xxx_xxzzz_1,   \
                             ta1_y_xxx_xxzzzz_0,  \
                             ta1_y_xxx_xxzzzz_1,  \
                             ta1_y_xxx_xyyyy_0,   \
                             ta1_y_xxx_xyyyy_1,   \
                             ta1_y_xxx_xyyyyy_0,  \
                             ta1_y_xxx_xyyyyy_1,  \
                             ta1_y_xxx_xyyyyz_0,  \
                             ta1_y_xxx_xyyyyz_1,  \
                             ta1_y_xxx_xyyyz_0,   \
                             ta1_y_xxx_xyyyz_1,   \
                             ta1_y_xxx_xyyyzz_0,  \
                             ta1_y_xxx_xyyyzz_1,  \
                             ta1_y_xxx_xyyzz_0,   \
                             ta1_y_xxx_xyyzz_1,   \
                             ta1_y_xxx_xyyzzz_0,  \
                             ta1_y_xxx_xyyzzz_1,  \
                             ta1_y_xxx_xyzzz_0,   \
                             ta1_y_xxx_xyzzz_1,   \
                             ta1_y_xxx_xyzzzz_0,  \
                             ta1_y_xxx_xyzzzz_1,  \
                             ta1_y_xxx_xzzzz_0,   \
                             ta1_y_xxx_xzzzz_1,   \
                             ta1_y_xxx_xzzzzz_0,  \
                             ta1_y_xxx_xzzzzz_1,  \
                             ta1_y_xxx_yyyyy_0,   \
                             ta1_y_xxx_yyyyy_1,   \
                             ta1_y_xxx_yyyyyy_0,  \
                             ta1_y_xxx_yyyyyy_1,  \
                             ta1_y_xxx_yyyyyz_0,  \
                             ta1_y_xxx_yyyyyz_1,  \
                             ta1_y_xxx_yyyyz_0,   \
                             ta1_y_xxx_yyyyz_1,   \
                             ta1_y_xxx_yyyyzz_0,  \
                             ta1_y_xxx_yyyyzz_1,  \
                             ta1_y_xxx_yyyzz_0,   \
                             ta1_y_xxx_yyyzz_1,   \
                             ta1_y_xxx_yyyzzz_0,  \
                             ta1_y_xxx_yyyzzz_1,  \
                             ta1_y_xxx_yyzzz_0,   \
                             ta1_y_xxx_yyzzz_1,   \
                             ta1_y_xxx_yyzzzz_0,  \
                             ta1_y_xxx_yyzzzz_1,  \
                             ta1_y_xxx_yzzzz_0,   \
                             ta1_y_xxx_yzzzz_1,   \
                             ta1_y_xxx_yzzzzz_0,  \
                             ta1_y_xxx_yzzzzz_1,  \
                             ta1_y_xxx_zzzzz_0,   \
                             ta1_y_xxx_zzzzz_1,   \
                             ta1_y_xxx_zzzzzz_0,  \
                             ta1_y_xxx_zzzzzz_1,  \
                             ta1_y_xxxx_xxxxxx_0, \
                             ta1_y_xxxx_xxxxxy_0, \
                             ta1_y_xxxx_xxxxxz_0, \
                             ta1_y_xxxx_xxxxyy_0, \
                             ta1_y_xxxx_xxxxyz_0, \
                             ta1_y_xxxx_xxxxzz_0, \
                             ta1_y_xxxx_xxxyyy_0, \
                             ta1_y_xxxx_xxxyyz_0, \
                             ta1_y_xxxx_xxxyzz_0, \
                             ta1_y_xxxx_xxxzzz_0, \
                             ta1_y_xxxx_xxyyyy_0, \
                             ta1_y_xxxx_xxyyyz_0, \
                             ta1_y_xxxx_xxyyzz_0, \
                             ta1_y_xxxx_xxyzzz_0, \
                             ta1_y_xxxx_xxzzzz_0, \
                             ta1_y_xxxx_xyyyyy_0, \
                             ta1_y_xxxx_xyyyyz_0, \
                             ta1_y_xxxx_xyyyzz_0, \
                             ta1_y_xxxx_xyyzzz_0, \
                             ta1_y_xxxx_xyzzzz_0, \
                             ta1_y_xxxx_xzzzzz_0, \
                             ta1_y_xxxx_yyyyyy_0, \
                             ta1_y_xxxx_yyyyyz_0, \
                             ta1_y_xxxx_yyyyzz_0, \
                             ta1_y_xxxx_yyyzzz_0, \
                             ta1_y_xxxx_yyzzzz_0, \
                             ta1_y_xxxx_yzzzzz_0, \
                             ta1_y_xxxx_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxx_xxxxxx_0[i] = 3.0 * ta1_y_xx_xxxxxx_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxxx_1[i] * fe_0 + 6.0 * ta1_y_xxx_xxxxx_0[i] * fe_0 -
                                 6.0 * ta1_y_xxx_xxxxx_1[i] * fe_0 + ta1_y_xxx_xxxxxx_0[i] * pa_x[i] - ta1_y_xxx_xxxxxx_1[i] * pc_x[i];

        ta1_y_xxxx_xxxxxy_0[i] = 3.0 * ta1_y_xx_xxxxxy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxxy_1[i] * fe_0 + 5.0 * ta1_y_xxx_xxxxy_0[i] * fe_0 -
                                 5.0 * ta1_y_xxx_xxxxy_1[i] * fe_0 + ta1_y_xxx_xxxxxy_0[i] * pa_x[i] - ta1_y_xxx_xxxxxy_1[i] * pc_x[i];

        ta1_y_xxxx_xxxxxz_0[i] = 3.0 * ta1_y_xx_xxxxxz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxxz_1[i] * fe_0 + 5.0 * ta1_y_xxx_xxxxz_0[i] * fe_0 -
                                 5.0 * ta1_y_xxx_xxxxz_1[i] * fe_0 + ta1_y_xxx_xxxxxz_0[i] * pa_x[i] - ta1_y_xxx_xxxxxz_1[i] * pc_x[i];

        ta1_y_xxxx_xxxxyy_0[i] = 3.0 * ta1_y_xx_xxxxyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxyy_1[i] * fe_0 + 4.0 * ta1_y_xxx_xxxyy_0[i] * fe_0 -
                                 4.0 * ta1_y_xxx_xxxyy_1[i] * fe_0 + ta1_y_xxx_xxxxyy_0[i] * pa_x[i] - ta1_y_xxx_xxxxyy_1[i] * pc_x[i];

        ta1_y_xxxx_xxxxyz_0[i] = 3.0 * ta1_y_xx_xxxxyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxyz_1[i] * fe_0 + 4.0 * ta1_y_xxx_xxxyz_0[i] * fe_0 -
                                 4.0 * ta1_y_xxx_xxxyz_1[i] * fe_0 + ta1_y_xxx_xxxxyz_0[i] * pa_x[i] - ta1_y_xxx_xxxxyz_1[i] * pc_x[i];

        ta1_y_xxxx_xxxxzz_0[i] = 3.0 * ta1_y_xx_xxxxzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxzz_1[i] * fe_0 + 4.0 * ta1_y_xxx_xxxzz_0[i] * fe_0 -
                                 4.0 * ta1_y_xxx_xxxzz_1[i] * fe_0 + ta1_y_xxx_xxxxzz_0[i] * pa_x[i] - ta1_y_xxx_xxxxzz_1[i] * pc_x[i];

        ta1_y_xxxx_xxxyyy_0[i] = 3.0 * ta1_y_xx_xxxyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxyyy_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxyyy_0[i] * fe_0 -
                                 3.0 * ta1_y_xxx_xxyyy_1[i] * fe_0 + ta1_y_xxx_xxxyyy_0[i] * pa_x[i] - ta1_y_xxx_xxxyyy_1[i] * pc_x[i];

        ta1_y_xxxx_xxxyyz_0[i] = 3.0 * ta1_y_xx_xxxyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxyyz_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_y_xxx_xxyyz_1[i] * fe_0 + ta1_y_xxx_xxxyyz_0[i] * pa_x[i] - ta1_y_xxx_xxxyyz_1[i] * pc_x[i];

        ta1_y_xxxx_xxxyzz_0[i] = 3.0 * ta1_y_xx_xxxyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxyzz_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_xxx_xxyzz_1[i] * fe_0 + ta1_y_xxx_xxxyzz_0[i] * pa_x[i] - ta1_y_xxx_xxxyzz_1[i] * pc_x[i];

        ta1_y_xxxx_xxxzzz_0[i] = 3.0 * ta1_y_xx_xxxzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxzzz_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxzzz_0[i] * fe_0 -
                                 3.0 * ta1_y_xxx_xxzzz_1[i] * fe_0 + ta1_y_xxx_xxxzzz_0[i] * pa_x[i] - ta1_y_xxx_xxxzzz_1[i] * pc_x[i];

        ta1_y_xxxx_xxyyyy_0[i] = 3.0 * ta1_y_xx_xxyyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyyyy_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyyyy_0[i] * fe_0 -
                                 2.0 * ta1_y_xxx_xyyyy_1[i] * fe_0 + ta1_y_xxx_xxyyyy_0[i] * pa_x[i] - ta1_y_xxx_xxyyyy_1[i] * pc_x[i];

        ta1_y_xxxx_xxyyyz_0[i] = 3.0 * ta1_y_xx_xxyyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyyyz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_xxx_xyyyz_1[i] * fe_0 + ta1_y_xxx_xxyyyz_0[i] * pa_x[i] - ta1_y_xxx_xxyyyz_1[i] * pc_x[i];

        ta1_y_xxxx_xxyyzz_0[i] = 3.0 * ta1_y_xx_xxyyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyyzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xxx_xyyzz_1[i] * fe_0 + ta1_y_xxx_xxyyzz_0[i] * pa_x[i] - ta1_y_xxx_xxyyzz_1[i] * pc_x[i];

        ta1_y_xxxx_xxyzzz_0[i] = 3.0 * ta1_y_xx_xxyzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyzzz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xxx_xyzzz_1[i] * fe_0 + ta1_y_xxx_xxyzzz_0[i] * pa_x[i] - ta1_y_xxx_xxyzzz_1[i] * pc_x[i];

        ta1_y_xxxx_xxzzzz_0[i] = 3.0 * ta1_y_xx_xxzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxzzzz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xzzzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xxx_xzzzz_1[i] * fe_0 + ta1_y_xxx_xxzzzz_0[i] * pa_x[i] - ta1_y_xxx_xxzzzz_1[i] * pc_x[i];

        ta1_y_xxxx_xyyyyy_0[i] = 3.0 * ta1_y_xx_xyyyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyyyy_1[i] * fe_0 + ta1_y_xxx_yyyyy_0[i] * fe_0 -
                                 ta1_y_xxx_yyyyy_1[i] * fe_0 + ta1_y_xxx_xyyyyy_0[i] * pa_x[i] - ta1_y_xxx_xyyyyy_1[i] * pc_x[i];

        ta1_y_xxxx_xyyyyz_0[i] = 3.0 * ta1_y_xx_xyyyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyyyz_1[i] * fe_0 + ta1_y_xxx_yyyyz_0[i] * fe_0 -
                                 ta1_y_xxx_yyyyz_1[i] * fe_0 + ta1_y_xxx_xyyyyz_0[i] * pa_x[i] - ta1_y_xxx_xyyyyz_1[i] * pc_x[i];

        ta1_y_xxxx_xyyyzz_0[i] = 3.0 * ta1_y_xx_xyyyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyyzz_1[i] * fe_0 + ta1_y_xxx_yyyzz_0[i] * fe_0 -
                                 ta1_y_xxx_yyyzz_1[i] * fe_0 + ta1_y_xxx_xyyyzz_0[i] * pa_x[i] - ta1_y_xxx_xyyyzz_1[i] * pc_x[i];

        ta1_y_xxxx_xyyzzz_0[i] = 3.0 * ta1_y_xx_xyyzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyzzz_1[i] * fe_0 + ta1_y_xxx_yyzzz_0[i] * fe_0 -
                                 ta1_y_xxx_yyzzz_1[i] * fe_0 + ta1_y_xxx_xyyzzz_0[i] * pa_x[i] - ta1_y_xxx_xyyzzz_1[i] * pc_x[i];

        ta1_y_xxxx_xyzzzz_0[i] = 3.0 * ta1_y_xx_xyzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyzzzz_1[i] * fe_0 + ta1_y_xxx_yzzzz_0[i] * fe_0 -
                                 ta1_y_xxx_yzzzz_1[i] * fe_0 + ta1_y_xxx_xyzzzz_0[i] * pa_x[i] - ta1_y_xxx_xyzzzz_1[i] * pc_x[i];

        ta1_y_xxxx_xzzzzz_0[i] = 3.0 * ta1_y_xx_xzzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xzzzzz_1[i] * fe_0 + ta1_y_xxx_zzzzz_0[i] * fe_0 -
                                 ta1_y_xxx_zzzzz_1[i] * fe_0 + ta1_y_xxx_xzzzzz_0[i] * pa_x[i] - ta1_y_xxx_xzzzzz_1[i] * pc_x[i];

        ta1_y_xxxx_yyyyyy_0[i] =
            3.0 * ta1_y_xx_yyyyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyyyy_1[i] * fe_0 + ta1_y_xxx_yyyyyy_0[i] * pa_x[i] - ta1_y_xxx_yyyyyy_1[i] * pc_x[i];

        ta1_y_xxxx_yyyyyz_0[i] =
            3.0 * ta1_y_xx_yyyyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyyyz_1[i] * fe_0 + ta1_y_xxx_yyyyyz_0[i] * pa_x[i] - ta1_y_xxx_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxxx_yyyyzz_0[i] =
            3.0 * ta1_y_xx_yyyyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyyzz_1[i] * fe_0 + ta1_y_xxx_yyyyzz_0[i] * pa_x[i] - ta1_y_xxx_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxxx_yyyzzz_0[i] =
            3.0 * ta1_y_xx_yyyzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyzzz_1[i] * fe_0 + ta1_y_xxx_yyyzzz_0[i] * pa_x[i] - ta1_y_xxx_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxxx_yyzzzz_0[i] =
            3.0 * ta1_y_xx_yyzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyzzzz_1[i] * fe_0 + ta1_y_xxx_yyzzzz_0[i] * pa_x[i] - ta1_y_xxx_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxxx_yzzzzz_0[i] =
            3.0 * ta1_y_xx_yzzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yzzzzz_1[i] * fe_0 + ta1_y_xxx_yzzzzz_0[i] * pa_x[i] - ta1_y_xxx_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxxx_zzzzzz_0[i] =
            3.0 * ta1_y_xx_zzzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_zzzzzz_1[i] * fe_0 + ta1_y_xxx_zzzzzz_0[i] * pa_x[i] - ta1_y_xxx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 448-476 components of targeted buffer : GI

    auto ta1_y_xxxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 448);

    auto ta1_y_xxxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 449);

    auto ta1_y_xxxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 450);

    auto ta1_y_xxxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 451);

    auto ta1_y_xxxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 452);

    auto ta1_y_xxxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 453);

    auto ta1_y_xxxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 454);

    auto ta1_y_xxxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 455);

    auto ta1_y_xxxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 456);

    auto ta1_y_xxxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 457);

    auto ta1_y_xxxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 458);

    auto ta1_y_xxxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 459);

    auto ta1_y_xxxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 460);

    auto ta1_y_xxxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 461);

    auto ta1_y_xxxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 462);

    auto ta1_y_xxxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 463);

    auto ta1_y_xxxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 464);

    auto ta1_y_xxxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 465);

    auto ta1_y_xxxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 466);

    auto ta1_y_xxxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 467);

    auto ta1_y_xxxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 468);

    auto ta1_y_xxxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 469);

    auto ta1_y_xxxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 470);

    auto ta1_y_xxxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 471);

    auto ta1_y_xxxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 472);

    auto ta1_y_xxxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 473);

    auto ta1_y_xxxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 474);

    auto ta1_y_xxxy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 475);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_y_xxx_xxxxx_0,   \
                             ta1_y_xxx_xxxxx_1,   \
                             ta1_y_xxx_xxxxxx_0,  \
                             ta1_y_xxx_xxxxxx_1,  \
                             ta1_y_xxx_xxxxxy_0,  \
                             ta1_y_xxx_xxxxxy_1,  \
                             ta1_y_xxx_xxxxxz_0,  \
                             ta1_y_xxx_xxxxxz_1,  \
                             ta1_y_xxx_xxxxy_0,   \
                             ta1_y_xxx_xxxxy_1,   \
                             ta1_y_xxx_xxxxyy_0,  \
                             ta1_y_xxx_xxxxyy_1,  \
                             ta1_y_xxx_xxxxyz_0,  \
                             ta1_y_xxx_xxxxyz_1,  \
                             ta1_y_xxx_xxxxz_0,   \
                             ta1_y_xxx_xxxxz_1,   \
                             ta1_y_xxx_xxxxzz_0,  \
                             ta1_y_xxx_xxxxzz_1,  \
                             ta1_y_xxx_xxxyy_0,   \
                             ta1_y_xxx_xxxyy_1,   \
                             ta1_y_xxx_xxxyyy_0,  \
                             ta1_y_xxx_xxxyyy_1,  \
                             ta1_y_xxx_xxxyyz_0,  \
                             ta1_y_xxx_xxxyyz_1,  \
                             ta1_y_xxx_xxxyz_0,   \
                             ta1_y_xxx_xxxyz_1,   \
                             ta1_y_xxx_xxxyzz_0,  \
                             ta1_y_xxx_xxxyzz_1,  \
                             ta1_y_xxx_xxxzz_0,   \
                             ta1_y_xxx_xxxzz_1,   \
                             ta1_y_xxx_xxxzzz_0,  \
                             ta1_y_xxx_xxxzzz_1,  \
                             ta1_y_xxx_xxyyy_0,   \
                             ta1_y_xxx_xxyyy_1,   \
                             ta1_y_xxx_xxyyyy_0,  \
                             ta1_y_xxx_xxyyyy_1,  \
                             ta1_y_xxx_xxyyyz_0,  \
                             ta1_y_xxx_xxyyyz_1,  \
                             ta1_y_xxx_xxyyz_0,   \
                             ta1_y_xxx_xxyyz_1,   \
                             ta1_y_xxx_xxyyzz_0,  \
                             ta1_y_xxx_xxyyzz_1,  \
                             ta1_y_xxx_xxyzz_0,   \
                             ta1_y_xxx_xxyzz_1,   \
                             ta1_y_xxx_xxyzzz_0,  \
                             ta1_y_xxx_xxyzzz_1,  \
                             ta1_y_xxx_xxzzz_0,   \
                             ta1_y_xxx_xxzzz_1,   \
                             ta1_y_xxx_xxzzzz_0,  \
                             ta1_y_xxx_xxzzzz_1,  \
                             ta1_y_xxx_xyyyy_0,   \
                             ta1_y_xxx_xyyyy_1,   \
                             ta1_y_xxx_xyyyyy_0,  \
                             ta1_y_xxx_xyyyyy_1,  \
                             ta1_y_xxx_xyyyyz_0,  \
                             ta1_y_xxx_xyyyyz_1,  \
                             ta1_y_xxx_xyyyz_0,   \
                             ta1_y_xxx_xyyyz_1,   \
                             ta1_y_xxx_xyyyzz_0,  \
                             ta1_y_xxx_xyyyzz_1,  \
                             ta1_y_xxx_xyyzz_0,   \
                             ta1_y_xxx_xyyzz_1,   \
                             ta1_y_xxx_xyyzzz_0,  \
                             ta1_y_xxx_xyyzzz_1,  \
                             ta1_y_xxx_xyzzz_0,   \
                             ta1_y_xxx_xyzzz_1,   \
                             ta1_y_xxx_xyzzzz_0,  \
                             ta1_y_xxx_xyzzzz_1,  \
                             ta1_y_xxx_xzzzz_0,   \
                             ta1_y_xxx_xzzzz_1,   \
                             ta1_y_xxx_xzzzzz_0,  \
                             ta1_y_xxx_xzzzzz_1,  \
                             ta1_y_xxx_zzzzzz_0,  \
                             ta1_y_xxx_zzzzzz_1,  \
                             ta1_y_xxxy_xxxxxx_0, \
                             ta1_y_xxxy_xxxxxy_0, \
                             ta1_y_xxxy_xxxxxz_0, \
                             ta1_y_xxxy_xxxxyy_0, \
                             ta1_y_xxxy_xxxxyz_0, \
                             ta1_y_xxxy_xxxxzz_0, \
                             ta1_y_xxxy_xxxyyy_0, \
                             ta1_y_xxxy_xxxyyz_0, \
                             ta1_y_xxxy_xxxyzz_0, \
                             ta1_y_xxxy_xxxzzz_0, \
                             ta1_y_xxxy_xxyyyy_0, \
                             ta1_y_xxxy_xxyyyz_0, \
                             ta1_y_xxxy_xxyyzz_0, \
                             ta1_y_xxxy_xxyzzz_0, \
                             ta1_y_xxxy_xxzzzz_0, \
                             ta1_y_xxxy_xyyyyy_0, \
                             ta1_y_xxxy_xyyyyz_0, \
                             ta1_y_xxxy_xyyyzz_0, \
                             ta1_y_xxxy_xyyzzz_0, \
                             ta1_y_xxxy_xyzzzz_0, \
                             ta1_y_xxxy_xzzzzz_0, \
                             ta1_y_xxxy_yyyyyy_0, \
                             ta1_y_xxxy_yyyyyz_0, \
                             ta1_y_xxxy_yyyyzz_0, \
                             ta1_y_xxxy_yyyzzz_0, \
                             ta1_y_xxxy_yyzzzz_0, \
                             ta1_y_xxxy_yzzzzz_0, \
                             ta1_y_xxxy_zzzzzz_0, \
                             ta1_y_xxy_yyyyyy_0,  \
                             ta1_y_xxy_yyyyyy_1,  \
                             ta1_y_xxy_yyyyyz_0,  \
                             ta1_y_xxy_yyyyyz_1,  \
                             ta1_y_xxy_yyyyzz_0,  \
                             ta1_y_xxy_yyyyzz_1,  \
                             ta1_y_xxy_yyyzzz_0,  \
                             ta1_y_xxy_yyyzzz_1,  \
                             ta1_y_xxy_yyzzzz_0,  \
                             ta1_y_xxy_yyzzzz_1,  \
                             ta1_y_xxy_yzzzzz_0,  \
                             ta1_y_xxy_yzzzzz_1,  \
                             ta1_y_xy_yyyyyy_0,   \
                             ta1_y_xy_yyyyyy_1,   \
                             ta1_y_xy_yyyyyz_0,   \
                             ta1_y_xy_yyyyyz_1,   \
                             ta1_y_xy_yyyyzz_0,   \
                             ta1_y_xy_yyyyzz_1,   \
                             ta1_y_xy_yyyzzz_0,   \
                             ta1_y_xy_yyyzzz_1,   \
                             ta1_y_xy_yyzzzz_0,   \
                             ta1_y_xy_yyzzzz_1,   \
                             ta1_y_xy_yzzzzz_0,   \
                             ta1_y_xy_yzzzzz_1,   \
                             ta_xxx_xxxxxx_1,     \
                             ta_xxx_xxxxxy_1,     \
                             ta_xxx_xxxxxz_1,     \
                             ta_xxx_xxxxyy_1,     \
                             ta_xxx_xxxxyz_1,     \
                             ta_xxx_xxxxzz_1,     \
                             ta_xxx_xxxyyy_1,     \
                             ta_xxx_xxxyyz_1,     \
                             ta_xxx_xxxyzz_1,     \
                             ta_xxx_xxxzzz_1,     \
                             ta_xxx_xxyyyy_1,     \
                             ta_xxx_xxyyyz_1,     \
                             ta_xxx_xxyyzz_1,     \
                             ta_xxx_xxyzzz_1,     \
                             ta_xxx_xxzzzz_1,     \
                             ta_xxx_xyyyyy_1,     \
                             ta_xxx_xyyyyz_1,     \
                             ta_xxx_xyyyzz_1,     \
                             ta_xxx_xyyzzz_1,     \
                             ta_xxx_xyzzzz_1,     \
                             ta_xxx_xzzzzz_1,     \
                             ta_xxx_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxy_xxxxxx_0[i] = ta_xxx_xxxxxx_1[i] + ta1_y_xxx_xxxxxx_0[i] * pa_y[i] - ta1_y_xxx_xxxxxx_1[i] * pc_y[i];

        ta1_y_xxxy_xxxxxy_0[i] = ta1_y_xxx_xxxxx_0[i] * fe_0 - ta1_y_xxx_xxxxx_1[i] * fe_0 + ta_xxx_xxxxxy_1[i] + ta1_y_xxx_xxxxxy_0[i] * pa_y[i] -
                                 ta1_y_xxx_xxxxxy_1[i] * pc_y[i];

        ta1_y_xxxy_xxxxxz_0[i] = ta_xxx_xxxxxz_1[i] + ta1_y_xxx_xxxxxz_0[i] * pa_y[i] - ta1_y_xxx_xxxxxz_1[i] * pc_y[i];

        ta1_y_xxxy_xxxxyy_0[i] = 2.0 * ta1_y_xxx_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxxxy_1[i] * fe_0 + ta_xxx_xxxxyy_1[i] +
                                 ta1_y_xxx_xxxxyy_0[i] * pa_y[i] - ta1_y_xxx_xxxxyy_1[i] * pc_y[i];

        ta1_y_xxxy_xxxxyz_0[i] = ta1_y_xxx_xxxxz_0[i] * fe_0 - ta1_y_xxx_xxxxz_1[i] * fe_0 + ta_xxx_xxxxyz_1[i] + ta1_y_xxx_xxxxyz_0[i] * pa_y[i] -
                                 ta1_y_xxx_xxxxyz_1[i] * pc_y[i];

        ta1_y_xxxy_xxxxzz_0[i] = ta_xxx_xxxxzz_1[i] + ta1_y_xxx_xxxxzz_0[i] * pa_y[i] - ta1_y_xxx_xxxxzz_1[i] * pc_y[i];

        ta1_y_xxxy_xxxyyy_0[i] = 3.0 * ta1_y_xxx_xxxyy_0[i] * fe_0 - 3.0 * ta1_y_xxx_xxxyy_1[i] * fe_0 + ta_xxx_xxxyyy_1[i] +
                                 ta1_y_xxx_xxxyyy_0[i] * pa_y[i] - ta1_y_xxx_xxxyyy_1[i] * pc_y[i];

        ta1_y_xxxy_xxxyyz_0[i] = 2.0 * ta1_y_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxxyz_1[i] * fe_0 + ta_xxx_xxxyyz_1[i] +
                                 ta1_y_xxx_xxxyyz_0[i] * pa_y[i] - ta1_y_xxx_xxxyyz_1[i] * pc_y[i];

        ta1_y_xxxy_xxxyzz_0[i] = ta1_y_xxx_xxxzz_0[i] * fe_0 - ta1_y_xxx_xxxzz_1[i] * fe_0 + ta_xxx_xxxyzz_1[i] + ta1_y_xxx_xxxyzz_0[i] * pa_y[i] -
                                 ta1_y_xxx_xxxyzz_1[i] * pc_y[i];

        ta1_y_xxxy_xxxzzz_0[i] = ta_xxx_xxxzzz_1[i] + ta1_y_xxx_xxxzzz_0[i] * pa_y[i] - ta1_y_xxx_xxxzzz_1[i] * pc_y[i];

        ta1_y_xxxy_xxyyyy_0[i] = 4.0 * ta1_y_xxx_xxyyy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxyyy_1[i] * fe_0 + ta_xxx_xxyyyy_1[i] +
                                 ta1_y_xxx_xxyyyy_0[i] * pa_y[i] - ta1_y_xxx_xxyyyy_1[i] * pc_y[i];

        ta1_y_xxxy_xxyyyz_0[i] = 3.0 * ta1_y_xxx_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xxyyz_1[i] * fe_0 + ta_xxx_xxyyyz_1[i] +
                                 ta1_y_xxx_xxyyyz_0[i] * pa_y[i] - ta1_y_xxx_xxyyyz_1[i] * pc_y[i];

        ta1_y_xxxy_xxyyzz_0[i] = 2.0 * ta1_y_xxx_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxyzz_1[i] * fe_0 + ta_xxx_xxyyzz_1[i] +
                                 ta1_y_xxx_xxyyzz_0[i] * pa_y[i] - ta1_y_xxx_xxyyzz_1[i] * pc_y[i];

        ta1_y_xxxy_xxyzzz_0[i] = ta1_y_xxx_xxzzz_0[i] * fe_0 - ta1_y_xxx_xxzzz_1[i] * fe_0 + ta_xxx_xxyzzz_1[i] + ta1_y_xxx_xxyzzz_0[i] * pa_y[i] -
                                 ta1_y_xxx_xxyzzz_1[i] * pc_y[i];

        ta1_y_xxxy_xxzzzz_0[i] = ta_xxx_xxzzzz_1[i] + ta1_y_xxx_xxzzzz_0[i] * pa_y[i] - ta1_y_xxx_xxzzzz_1[i] * pc_y[i];

        ta1_y_xxxy_xyyyyy_0[i] = 5.0 * ta1_y_xxx_xyyyy_0[i] * fe_0 - 5.0 * ta1_y_xxx_xyyyy_1[i] * fe_0 + ta_xxx_xyyyyy_1[i] +
                                 ta1_y_xxx_xyyyyy_0[i] * pa_y[i] - ta1_y_xxx_xyyyyy_1[i] * pc_y[i];

        ta1_y_xxxy_xyyyyz_0[i] = 4.0 * ta1_y_xxx_xyyyz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyyyz_1[i] * fe_0 + ta_xxx_xyyyyz_1[i] +
                                 ta1_y_xxx_xyyyyz_0[i] * pa_y[i] - ta1_y_xxx_xyyyyz_1[i] * pc_y[i];

        ta1_y_xxxy_xyyyzz_0[i] = 3.0 * ta1_y_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xyyzz_1[i] * fe_0 + ta_xxx_xyyyzz_1[i] +
                                 ta1_y_xxx_xyyyzz_0[i] * pa_y[i] - ta1_y_xxx_xyyyzz_1[i] * pc_y[i];

        ta1_y_xxxy_xyyzzz_0[i] = 2.0 * ta1_y_xxx_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xyzzz_1[i] * fe_0 + ta_xxx_xyyzzz_1[i] +
                                 ta1_y_xxx_xyyzzz_0[i] * pa_y[i] - ta1_y_xxx_xyyzzz_1[i] * pc_y[i];

        ta1_y_xxxy_xyzzzz_0[i] = ta1_y_xxx_xzzzz_0[i] * fe_0 - ta1_y_xxx_xzzzz_1[i] * fe_0 + ta_xxx_xyzzzz_1[i] + ta1_y_xxx_xyzzzz_0[i] * pa_y[i] -
                                 ta1_y_xxx_xyzzzz_1[i] * pc_y[i];

        ta1_y_xxxy_xzzzzz_0[i] = ta_xxx_xzzzzz_1[i] + ta1_y_xxx_xzzzzz_0[i] * pa_y[i] - ta1_y_xxx_xzzzzz_1[i] * pc_y[i];

        ta1_y_xxxy_yyyyyy_0[i] =
            2.0 * ta1_y_xy_yyyyyy_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyyyy_1[i] * fe_0 + ta1_y_xxy_yyyyyy_0[i] * pa_x[i] - ta1_y_xxy_yyyyyy_1[i] * pc_x[i];

        ta1_y_xxxy_yyyyyz_0[i] =
            2.0 * ta1_y_xy_yyyyyz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyyyz_1[i] * fe_0 + ta1_y_xxy_yyyyyz_0[i] * pa_x[i] - ta1_y_xxy_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxxy_yyyyzz_0[i] =
            2.0 * ta1_y_xy_yyyyzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyyzz_1[i] * fe_0 + ta1_y_xxy_yyyyzz_0[i] * pa_x[i] - ta1_y_xxy_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxxy_yyyzzz_0[i] =
            2.0 * ta1_y_xy_yyyzzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyzzz_1[i] * fe_0 + ta1_y_xxy_yyyzzz_0[i] * pa_x[i] - ta1_y_xxy_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxxy_yyzzzz_0[i] =
            2.0 * ta1_y_xy_yyzzzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyzzzz_1[i] * fe_0 + ta1_y_xxy_yyzzzz_0[i] * pa_x[i] - ta1_y_xxy_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxxy_yzzzzz_0[i] =
            2.0 * ta1_y_xy_yzzzzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yzzzzz_1[i] * fe_0 + ta1_y_xxy_yzzzzz_0[i] * pa_x[i] - ta1_y_xxy_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxxy_zzzzzz_0[i] = ta_xxx_zzzzzz_1[i] + ta1_y_xxx_zzzzzz_0[i] * pa_y[i] - ta1_y_xxx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 476-504 components of targeted buffer : GI

    auto ta1_y_xxxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 476);

    auto ta1_y_xxxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 477);

    auto ta1_y_xxxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 478);

    auto ta1_y_xxxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 479);

    auto ta1_y_xxxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 480);

    auto ta1_y_xxxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 481);

    auto ta1_y_xxxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 482);

    auto ta1_y_xxxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 483);

    auto ta1_y_xxxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 484);

    auto ta1_y_xxxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 485);

    auto ta1_y_xxxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 486);

    auto ta1_y_xxxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 487);

    auto ta1_y_xxxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 488);

    auto ta1_y_xxxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 489);

    auto ta1_y_xxxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 490);

    auto ta1_y_xxxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 491);

    auto ta1_y_xxxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 492);

    auto ta1_y_xxxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 493);

    auto ta1_y_xxxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 494);

    auto ta1_y_xxxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 495);

    auto ta1_y_xxxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 496);

    auto ta1_y_xxxz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 497);

    auto ta1_y_xxxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 498);

    auto ta1_y_xxxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 499);

    auto ta1_y_xxxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 500);

    auto ta1_y_xxxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 501);

    auto ta1_y_xxxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 502);

    auto ta1_y_xxxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 503);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_y_xxx_xxxxx_0,   \
                             ta1_y_xxx_xxxxx_1,   \
                             ta1_y_xxx_xxxxxx_0,  \
                             ta1_y_xxx_xxxxxx_1,  \
                             ta1_y_xxx_xxxxxy_0,  \
                             ta1_y_xxx_xxxxxy_1,  \
                             ta1_y_xxx_xxxxxz_0,  \
                             ta1_y_xxx_xxxxxz_1,  \
                             ta1_y_xxx_xxxxy_0,   \
                             ta1_y_xxx_xxxxy_1,   \
                             ta1_y_xxx_xxxxyy_0,  \
                             ta1_y_xxx_xxxxyy_1,  \
                             ta1_y_xxx_xxxxyz_0,  \
                             ta1_y_xxx_xxxxyz_1,  \
                             ta1_y_xxx_xxxxz_0,   \
                             ta1_y_xxx_xxxxz_1,   \
                             ta1_y_xxx_xxxxzz_0,  \
                             ta1_y_xxx_xxxxzz_1,  \
                             ta1_y_xxx_xxxyy_0,   \
                             ta1_y_xxx_xxxyy_1,   \
                             ta1_y_xxx_xxxyyy_0,  \
                             ta1_y_xxx_xxxyyy_1,  \
                             ta1_y_xxx_xxxyyz_0,  \
                             ta1_y_xxx_xxxyyz_1,  \
                             ta1_y_xxx_xxxyz_0,   \
                             ta1_y_xxx_xxxyz_1,   \
                             ta1_y_xxx_xxxyzz_0,  \
                             ta1_y_xxx_xxxyzz_1,  \
                             ta1_y_xxx_xxxzz_0,   \
                             ta1_y_xxx_xxxzz_1,   \
                             ta1_y_xxx_xxxzzz_0,  \
                             ta1_y_xxx_xxxzzz_1,  \
                             ta1_y_xxx_xxyyy_0,   \
                             ta1_y_xxx_xxyyy_1,   \
                             ta1_y_xxx_xxyyyy_0,  \
                             ta1_y_xxx_xxyyyy_1,  \
                             ta1_y_xxx_xxyyyz_0,  \
                             ta1_y_xxx_xxyyyz_1,  \
                             ta1_y_xxx_xxyyz_0,   \
                             ta1_y_xxx_xxyyz_1,   \
                             ta1_y_xxx_xxyyzz_0,  \
                             ta1_y_xxx_xxyyzz_1,  \
                             ta1_y_xxx_xxyzz_0,   \
                             ta1_y_xxx_xxyzz_1,   \
                             ta1_y_xxx_xxyzzz_0,  \
                             ta1_y_xxx_xxyzzz_1,  \
                             ta1_y_xxx_xxzzz_0,   \
                             ta1_y_xxx_xxzzz_1,   \
                             ta1_y_xxx_xxzzzz_0,  \
                             ta1_y_xxx_xxzzzz_1,  \
                             ta1_y_xxx_xyyyy_0,   \
                             ta1_y_xxx_xyyyy_1,   \
                             ta1_y_xxx_xyyyyy_0,  \
                             ta1_y_xxx_xyyyyy_1,  \
                             ta1_y_xxx_xyyyyz_0,  \
                             ta1_y_xxx_xyyyyz_1,  \
                             ta1_y_xxx_xyyyz_0,   \
                             ta1_y_xxx_xyyyz_1,   \
                             ta1_y_xxx_xyyyzz_0,  \
                             ta1_y_xxx_xyyyzz_1,  \
                             ta1_y_xxx_xyyzz_0,   \
                             ta1_y_xxx_xyyzz_1,   \
                             ta1_y_xxx_xyyzzz_0,  \
                             ta1_y_xxx_xyyzzz_1,  \
                             ta1_y_xxx_xyzzz_0,   \
                             ta1_y_xxx_xyzzz_1,   \
                             ta1_y_xxx_xyzzzz_0,  \
                             ta1_y_xxx_xyzzzz_1,  \
                             ta1_y_xxx_xzzzz_0,   \
                             ta1_y_xxx_xzzzz_1,   \
                             ta1_y_xxx_xzzzzz_0,  \
                             ta1_y_xxx_xzzzzz_1,  \
                             ta1_y_xxx_yyyyyy_0,  \
                             ta1_y_xxx_yyyyyy_1,  \
                             ta1_y_xxxz_xxxxxx_0, \
                             ta1_y_xxxz_xxxxxy_0, \
                             ta1_y_xxxz_xxxxxz_0, \
                             ta1_y_xxxz_xxxxyy_0, \
                             ta1_y_xxxz_xxxxyz_0, \
                             ta1_y_xxxz_xxxxzz_0, \
                             ta1_y_xxxz_xxxyyy_0, \
                             ta1_y_xxxz_xxxyyz_0, \
                             ta1_y_xxxz_xxxyzz_0, \
                             ta1_y_xxxz_xxxzzz_0, \
                             ta1_y_xxxz_xxyyyy_0, \
                             ta1_y_xxxz_xxyyyz_0, \
                             ta1_y_xxxz_xxyyzz_0, \
                             ta1_y_xxxz_xxyzzz_0, \
                             ta1_y_xxxz_xxzzzz_0, \
                             ta1_y_xxxz_xyyyyy_0, \
                             ta1_y_xxxz_xyyyyz_0, \
                             ta1_y_xxxz_xyyyzz_0, \
                             ta1_y_xxxz_xyyzzz_0, \
                             ta1_y_xxxz_xyzzzz_0, \
                             ta1_y_xxxz_xzzzzz_0, \
                             ta1_y_xxxz_yyyyyy_0, \
                             ta1_y_xxxz_yyyyyz_0, \
                             ta1_y_xxxz_yyyyzz_0, \
                             ta1_y_xxxz_yyyzzz_0, \
                             ta1_y_xxxz_yyzzzz_0, \
                             ta1_y_xxxz_yzzzzz_0, \
                             ta1_y_xxxz_zzzzzz_0, \
                             ta1_y_xxz_yyyyyz_0,  \
                             ta1_y_xxz_yyyyyz_1,  \
                             ta1_y_xxz_yyyyzz_0,  \
                             ta1_y_xxz_yyyyzz_1,  \
                             ta1_y_xxz_yyyzzz_0,  \
                             ta1_y_xxz_yyyzzz_1,  \
                             ta1_y_xxz_yyzzzz_0,  \
                             ta1_y_xxz_yyzzzz_1,  \
                             ta1_y_xxz_yzzzzz_0,  \
                             ta1_y_xxz_yzzzzz_1,  \
                             ta1_y_xxz_zzzzzz_0,  \
                             ta1_y_xxz_zzzzzz_1,  \
                             ta1_y_xz_yyyyyz_0,   \
                             ta1_y_xz_yyyyyz_1,   \
                             ta1_y_xz_yyyyzz_0,   \
                             ta1_y_xz_yyyyzz_1,   \
                             ta1_y_xz_yyyzzz_0,   \
                             ta1_y_xz_yyyzzz_1,   \
                             ta1_y_xz_yyzzzz_0,   \
                             ta1_y_xz_yyzzzz_1,   \
                             ta1_y_xz_yzzzzz_0,   \
                             ta1_y_xz_yzzzzz_1,   \
                             ta1_y_xz_zzzzzz_0,   \
                             ta1_y_xz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxz_xxxxxx_0[i] = ta1_y_xxx_xxxxxx_0[i] * pa_z[i] - ta1_y_xxx_xxxxxx_1[i] * pc_z[i];

        ta1_y_xxxz_xxxxxy_0[i] = ta1_y_xxx_xxxxxy_0[i] * pa_z[i] - ta1_y_xxx_xxxxxy_1[i] * pc_z[i];

        ta1_y_xxxz_xxxxxz_0[i] =
            ta1_y_xxx_xxxxx_0[i] * fe_0 - ta1_y_xxx_xxxxx_1[i] * fe_0 + ta1_y_xxx_xxxxxz_0[i] * pa_z[i] - ta1_y_xxx_xxxxxz_1[i] * pc_z[i];

        ta1_y_xxxz_xxxxyy_0[i] = ta1_y_xxx_xxxxyy_0[i] * pa_z[i] - ta1_y_xxx_xxxxyy_1[i] * pc_z[i];

        ta1_y_xxxz_xxxxyz_0[i] =
            ta1_y_xxx_xxxxy_0[i] * fe_0 - ta1_y_xxx_xxxxy_1[i] * fe_0 + ta1_y_xxx_xxxxyz_0[i] * pa_z[i] - ta1_y_xxx_xxxxyz_1[i] * pc_z[i];

        ta1_y_xxxz_xxxxzz_0[i] =
            2.0 * ta1_y_xxx_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxxxz_1[i] * fe_0 + ta1_y_xxx_xxxxzz_0[i] * pa_z[i] - ta1_y_xxx_xxxxzz_1[i] * pc_z[i];

        ta1_y_xxxz_xxxyyy_0[i] = ta1_y_xxx_xxxyyy_0[i] * pa_z[i] - ta1_y_xxx_xxxyyy_1[i] * pc_z[i];

        ta1_y_xxxz_xxxyyz_0[i] =
            ta1_y_xxx_xxxyy_0[i] * fe_0 - ta1_y_xxx_xxxyy_1[i] * fe_0 + ta1_y_xxx_xxxyyz_0[i] * pa_z[i] - ta1_y_xxx_xxxyyz_1[i] * pc_z[i];

        ta1_y_xxxz_xxxyzz_0[i] =
            2.0 * ta1_y_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxxyz_1[i] * fe_0 + ta1_y_xxx_xxxyzz_0[i] * pa_z[i] - ta1_y_xxx_xxxyzz_1[i] * pc_z[i];

        ta1_y_xxxz_xxxzzz_0[i] =
            3.0 * ta1_y_xxx_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xxxzz_1[i] * fe_0 + ta1_y_xxx_xxxzzz_0[i] * pa_z[i] - ta1_y_xxx_xxxzzz_1[i] * pc_z[i];

        ta1_y_xxxz_xxyyyy_0[i] = ta1_y_xxx_xxyyyy_0[i] * pa_z[i] - ta1_y_xxx_xxyyyy_1[i] * pc_z[i];

        ta1_y_xxxz_xxyyyz_0[i] =
            ta1_y_xxx_xxyyy_0[i] * fe_0 - ta1_y_xxx_xxyyy_1[i] * fe_0 + ta1_y_xxx_xxyyyz_0[i] * pa_z[i] - ta1_y_xxx_xxyyyz_1[i] * pc_z[i];

        ta1_y_xxxz_xxyyzz_0[i] =
            2.0 * ta1_y_xxx_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxyyz_1[i] * fe_0 + ta1_y_xxx_xxyyzz_0[i] * pa_z[i] - ta1_y_xxx_xxyyzz_1[i] * pc_z[i];

        ta1_y_xxxz_xxyzzz_0[i] =
            3.0 * ta1_y_xxx_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xxyzz_1[i] * fe_0 + ta1_y_xxx_xxyzzz_0[i] * pa_z[i] - ta1_y_xxx_xxyzzz_1[i] * pc_z[i];

        ta1_y_xxxz_xxzzzz_0[i] =
            4.0 * ta1_y_xxx_xxzzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xxzzz_1[i] * fe_0 + ta1_y_xxx_xxzzzz_0[i] * pa_z[i] - ta1_y_xxx_xxzzzz_1[i] * pc_z[i];

        ta1_y_xxxz_xyyyyy_0[i] = ta1_y_xxx_xyyyyy_0[i] * pa_z[i] - ta1_y_xxx_xyyyyy_1[i] * pc_z[i];

        ta1_y_xxxz_xyyyyz_0[i] =
            ta1_y_xxx_xyyyy_0[i] * fe_0 - ta1_y_xxx_xyyyy_1[i] * fe_0 + ta1_y_xxx_xyyyyz_0[i] * pa_z[i] - ta1_y_xxx_xyyyyz_1[i] * pc_z[i];

        ta1_y_xxxz_xyyyzz_0[i] =
            2.0 * ta1_y_xxx_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xyyyz_1[i] * fe_0 + ta1_y_xxx_xyyyzz_0[i] * pa_z[i] - ta1_y_xxx_xyyyzz_1[i] * pc_z[i];

        ta1_y_xxxz_xyyzzz_0[i] =
            3.0 * ta1_y_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xyyzz_1[i] * fe_0 + ta1_y_xxx_xyyzzz_0[i] * pa_z[i] - ta1_y_xxx_xyyzzz_1[i] * pc_z[i];

        ta1_y_xxxz_xyzzzz_0[i] =
            4.0 * ta1_y_xxx_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyzzz_1[i] * fe_0 + ta1_y_xxx_xyzzzz_0[i] * pa_z[i] - ta1_y_xxx_xyzzzz_1[i] * pc_z[i];

        ta1_y_xxxz_xzzzzz_0[i] =
            5.0 * ta1_y_xxx_xzzzz_0[i] * fe_0 - 5.0 * ta1_y_xxx_xzzzz_1[i] * fe_0 + ta1_y_xxx_xzzzzz_0[i] * pa_z[i] - ta1_y_xxx_xzzzzz_1[i] * pc_z[i];

        ta1_y_xxxz_yyyyyy_0[i] = ta1_y_xxx_yyyyyy_0[i] * pa_z[i] - ta1_y_xxx_yyyyyy_1[i] * pc_z[i];

        ta1_y_xxxz_yyyyyz_0[i] =
            2.0 * ta1_y_xz_yyyyyz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyyyyz_1[i] * fe_0 + ta1_y_xxz_yyyyyz_0[i] * pa_x[i] - ta1_y_xxz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxxz_yyyyzz_0[i] =
            2.0 * ta1_y_xz_yyyyzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyyyzz_1[i] * fe_0 + ta1_y_xxz_yyyyzz_0[i] * pa_x[i] - ta1_y_xxz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxxz_yyyzzz_0[i] =
            2.0 * ta1_y_xz_yyyzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyyzzz_1[i] * fe_0 + ta1_y_xxz_yyyzzz_0[i] * pa_x[i] - ta1_y_xxz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxxz_yyzzzz_0[i] =
            2.0 * ta1_y_xz_yyzzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyzzzz_1[i] * fe_0 + ta1_y_xxz_yyzzzz_0[i] * pa_x[i] - ta1_y_xxz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxxz_yzzzzz_0[i] =
            2.0 * ta1_y_xz_yzzzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yzzzzz_1[i] * fe_0 + ta1_y_xxz_yzzzzz_0[i] * pa_x[i] - ta1_y_xxz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxxz_zzzzzz_0[i] =
            2.0 * ta1_y_xz_zzzzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_zzzzzz_1[i] * fe_0 + ta1_y_xxz_zzzzzz_0[i] * pa_x[i] - ta1_y_xxz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 504-532 components of targeted buffer : GI

    auto ta1_y_xxyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 504);

    auto ta1_y_xxyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 505);

    auto ta1_y_xxyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 506);

    auto ta1_y_xxyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 507);

    auto ta1_y_xxyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 508);

    auto ta1_y_xxyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 509);

    auto ta1_y_xxyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 510);

    auto ta1_y_xxyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 511);

    auto ta1_y_xxyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 512);

    auto ta1_y_xxyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 513);

    auto ta1_y_xxyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 514);

    auto ta1_y_xxyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 515);

    auto ta1_y_xxyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 516);

    auto ta1_y_xxyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 517);

    auto ta1_y_xxyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 518);

    auto ta1_y_xxyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 519);

    auto ta1_y_xxyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 520);

    auto ta1_y_xxyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 521);

    auto ta1_y_xxyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 522);

    auto ta1_y_xxyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 523);

    auto ta1_y_xxyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 524);

    auto ta1_y_xxyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 525);

    auto ta1_y_xxyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 526);

    auto ta1_y_xxyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 527);

    auto ta1_y_xxyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 528);

    auto ta1_y_xxyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 529);

    auto ta1_y_xxyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 530);

    auto ta1_y_xxyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 531);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_y_xx_xxxxxx_0,   \
                             ta1_y_xx_xxxxxx_1,   \
                             ta1_y_xx_xxxxxz_0,   \
                             ta1_y_xx_xxxxxz_1,   \
                             ta1_y_xx_xxxxzz_0,   \
                             ta1_y_xx_xxxxzz_1,   \
                             ta1_y_xx_xxxzzz_0,   \
                             ta1_y_xx_xxxzzz_1,   \
                             ta1_y_xx_xxzzzz_0,   \
                             ta1_y_xx_xxzzzz_1,   \
                             ta1_y_xx_xzzzzz_0,   \
                             ta1_y_xx_xzzzzz_1,   \
                             ta1_y_xxy_xxxxxx_0,  \
                             ta1_y_xxy_xxxxxx_1,  \
                             ta1_y_xxy_xxxxxz_0,  \
                             ta1_y_xxy_xxxxxz_1,  \
                             ta1_y_xxy_xxxxzz_0,  \
                             ta1_y_xxy_xxxxzz_1,  \
                             ta1_y_xxy_xxxzzz_0,  \
                             ta1_y_xxy_xxxzzz_1,  \
                             ta1_y_xxy_xxzzzz_0,  \
                             ta1_y_xxy_xxzzzz_1,  \
                             ta1_y_xxy_xzzzzz_0,  \
                             ta1_y_xxy_xzzzzz_1,  \
                             ta1_y_xxyy_xxxxxx_0, \
                             ta1_y_xxyy_xxxxxy_0, \
                             ta1_y_xxyy_xxxxxz_0, \
                             ta1_y_xxyy_xxxxyy_0, \
                             ta1_y_xxyy_xxxxyz_0, \
                             ta1_y_xxyy_xxxxzz_0, \
                             ta1_y_xxyy_xxxyyy_0, \
                             ta1_y_xxyy_xxxyyz_0, \
                             ta1_y_xxyy_xxxyzz_0, \
                             ta1_y_xxyy_xxxzzz_0, \
                             ta1_y_xxyy_xxyyyy_0, \
                             ta1_y_xxyy_xxyyyz_0, \
                             ta1_y_xxyy_xxyyzz_0, \
                             ta1_y_xxyy_xxyzzz_0, \
                             ta1_y_xxyy_xxzzzz_0, \
                             ta1_y_xxyy_xyyyyy_0, \
                             ta1_y_xxyy_xyyyyz_0, \
                             ta1_y_xxyy_xyyyzz_0, \
                             ta1_y_xxyy_xyyzzz_0, \
                             ta1_y_xxyy_xyzzzz_0, \
                             ta1_y_xxyy_xzzzzz_0, \
                             ta1_y_xxyy_yyyyyy_0, \
                             ta1_y_xxyy_yyyyyz_0, \
                             ta1_y_xxyy_yyyyzz_0, \
                             ta1_y_xxyy_yyyzzz_0, \
                             ta1_y_xxyy_yyzzzz_0, \
                             ta1_y_xxyy_yzzzzz_0, \
                             ta1_y_xxyy_zzzzzz_0, \
                             ta1_y_xyy_xxxxxy_0,  \
                             ta1_y_xyy_xxxxxy_1,  \
                             ta1_y_xyy_xxxxy_0,   \
                             ta1_y_xyy_xxxxy_1,   \
                             ta1_y_xyy_xxxxyy_0,  \
                             ta1_y_xyy_xxxxyy_1,  \
                             ta1_y_xyy_xxxxyz_0,  \
                             ta1_y_xyy_xxxxyz_1,  \
                             ta1_y_xyy_xxxyy_0,   \
                             ta1_y_xyy_xxxyy_1,   \
                             ta1_y_xyy_xxxyyy_0,  \
                             ta1_y_xyy_xxxyyy_1,  \
                             ta1_y_xyy_xxxyyz_0,  \
                             ta1_y_xyy_xxxyyz_1,  \
                             ta1_y_xyy_xxxyz_0,   \
                             ta1_y_xyy_xxxyz_1,   \
                             ta1_y_xyy_xxxyzz_0,  \
                             ta1_y_xyy_xxxyzz_1,  \
                             ta1_y_xyy_xxyyy_0,   \
                             ta1_y_xyy_xxyyy_1,   \
                             ta1_y_xyy_xxyyyy_0,  \
                             ta1_y_xyy_xxyyyy_1,  \
                             ta1_y_xyy_xxyyyz_0,  \
                             ta1_y_xyy_xxyyyz_1,  \
                             ta1_y_xyy_xxyyz_0,   \
                             ta1_y_xyy_xxyyz_1,   \
                             ta1_y_xyy_xxyyzz_0,  \
                             ta1_y_xyy_xxyyzz_1,  \
                             ta1_y_xyy_xxyzz_0,   \
                             ta1_y_xyy_xxyzz_1,   \
                             ta1_y_xyy_xxyzzz_0,  \
                             ta1_y_xyy_xxyzzz_1,  \
                             ta1_y_xyy_xyyyy_0,   \
                             ta1_y_xyy_xyyyy_1,   \
                             ta1_y_xyy_xyyyyy_0,  \
                             ta1_y_xyy_xyyyyy_1,  \
                             ta1_y_xyy_xyyyyz_0,  \
                             ta1_y_xyy_xyyyyz_1,  \
                             ta1_y_xyy_xyyyz_0,   \
                             ta1_y_xyy_xyyyz_1,   \
                             ta1_y_xyy_xyyyzz_0,  \
                             ta1_y_xyy_xyyyzz_1,  \
                             ta1_y_xyy_xyyzz_0,   \
                             ta1_y_xyy_xyyzz_1,   \
                             ta1_y_xyy_xyyzzz_0,  \
                             ta1_y_xyy_xyyzzz_1,  \
                             ta1_y_xyy_xyzzz_0,   \
                             ta1_y_xyy_xyzzz_1,   \
                             ta1_y_xyy_xyzzzz_0,  \
                             ta1_y_xyy_xyzzzz_1,  \
                             ta1_y_xyy_yyyyy_0,   \
                             ta1_y_xyy_yyyyy_1,   \
                             ta1_y_xyy_yyyyyy_0,  \
                             ta1_y_xyy_yyyyyy_1,  \
                             ta1_y_xyy_yyyyyz_0,  \
                             ta1_y_xyy_yyyyyz_1,  \
                             ta1_y_xyy_yyyyz_0,   \
                             ta1_y_xyy_yyyyz_1,   \
                             ta1_y_xyy_yyyyzz_0,  \
                             ta1_y_xyy_yyyyzz_1,  \
                             ta1_y_xyy_yyyzz_0,   \
                             ta1_y_xyy_yyyzz_1,   \
                             ta1_y_xyy_yyyzzz_0,  \
                             ta1_y_xyy_yyyzzz_1,  \
                             ta1_y_xyy_yyzzz_0,   \
                             ta1_y_xyy_yyzzz_1,   \
                             ta1_y_xyy_yyzzzz_0,  \
                             ta1_y_xyy_yyzzzz_1,  \
                             ta1_y_xyy_yzzzz_0,   \
                             ta1_y_xyy_yzzzz_1,   \
                             ta1_y_xyy_yzzzzz_0,  \
                             ta1_y_xyy_yzzzzz_1,  \
                             ta1_y_xyy_zzzzzz_0,  \
                             ta1_y_xyy_zzzzzz_1,  \
                             ta1_y_yy_xxxxxy_0,   \
                             ta1_y_yy_xxxxxy_1,   \
                             ta1_y_yy_xxxxyy_0,   \
                             ta1_y_yy_xxxxyy_1,   \
                             ta1_y_yy_xxxxyz_0,   \
                             ta1_y_yy_xxxxyz_1,   \
                             ta1_y_yy_xxxyyy_0,   \
                             ta1_y_yy_xxxyyy_1,   \
                             ta1_y_yy_xxxyyz_0,   \
                             ta1_y_yy_xxxyyz_1,   \
                             ta1_y_yy_xxxyzz_0,   \
                             ta1_y_yy_xxxyzz_1,   \
                             ta1_y_yy_xxyyyy_0,   \
                             ta1_y_yy_xxyyyy_1,   \
                             ta1_y_yy_xxyyyz_0,   \
                             ta1_y_yy_xxyyyz_1,   \
                             ta1_y_yy_xxyyzz_0,   \
                             ta1_y_yy_xxyyzz_1,   \
                             ta1_y_yy_xxyzzz_0,   \
                             ta1_y_yy_xxyzzz_1,   \
                             ta1_y_yy_xyyyyy_0,   \
                             ta1_y_yy_xyyyyy_1,   \
                             ta1_y_yy_xyyyyz_0,   \
                             ta1_y_yy_xyyyyz_1,   \
                             ta1_y_yy_xyyyzz_0,   \
                             ta1_y_yy_xyyyzz_1,   \
                             ta1_y_yy_xyyzzz_0,   \
                             ta1_y_yy_xyyzzz_1,   \
                             ta1_y_yy_xyzzzz_0,   \
                             ta1_y_yy_xyzzzz_1,   \
                             ta1_y_yy_yyyyyy_0,   \
                             ta1_y_yy_yyyyyy_1,   \
                             ta1_y_yy_yyyyyz_0,   \
                             ta1_y_yy_yyyyyz_1,   \
                             ta1_y_yy_yyyyzz_0,   \
                             ta1_y_yy_yyyyzz_1,   \
                             ta1_y_yy_yyyzzz_0,   \
                             ta1_y_yy_yyyzzz_1,   \
                             ta1_y_yy_yyzzzz_0,   \
                             ta1_y_yy_yyzzzz_1,   \
                             ta1_y_yy_yzzzzz_0,   \
                             ta1_y_yy_yzzzzz_1,   \
                             ta1_y_yy_zzzzzz_0,   \
                             ta1_y_yy_zzzzzz_1,   \
                             ta_xxy_xxxxxx_1,     \
                             ta_xxy_xxxxxz_1,     \
                             ta_xxy_xxxxzz_1,     \
                             ta_xxy_xxxzzz_1,     \
                             ta_xxy_xxzzzz_1,     \
                             ta_xxy_xzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyy_xxxxxx_0[i] = ta1_y_xx_xxxxxx_0[i] * fe_0 - ta1_y_xx_xxxxxx_1[i] * fe_0 + ta_xxy_xxxxxx_1[i] + ta1_y_xxy_xxxxxx_0[i] * pa_y[i] -
                                 ta1_y_xxy_xxxxxx_1[i] * pc_y[i];

        ta1_y_xxyy_xxxxxy_0[i] = ta1_y_yy_xxxxxy_0[i] * fe_0 - ta1_y_yy_xxxxxy_1[i] * fe_0 + 5.0 * ta1_y_xyy_xxxxy_0[i] * fe_0 -
                                 5.0 * ta1_y_xyy_xxxxy_1[i] * fe_0 + ta1_y_xyy_xxxxxy_0[i] * pa_x[i] - ta1_y_xyy_xxxxxy_1[i] * pc_x[i];

        ta1_y_xxyy_xxxxxz_0[i] = ta1_y_xx_xxxxxz_0[i] * fe_0 - ta1_y_xx_xxxxxz_1[i] * fe_0 + ta_xxy_xxxxxz_1[i] + ta1_y_xxy_xxxxxz_0[i] * pa_y[i] -
                                 ta1_y_xxy_xxxxxz_1[i] * pc_y[i];

        ta1_y_xxyy_xxxxyy_0[i] = ta1_y_yy_xxxxyy_0[i] * fe_0 - ta1_y_yy_xxxxyy_1[i] * fe_0 + 4.0 * ta1_y_xyy_xxxyy_0[i] * fe_0 -
                                 4.0 * ta1_y_xyy_xxxyy_1[i] * fe_0 + ta1_y_xyy_xxxxyy_0[i] * pa_x[i] - ta1_y_xyy_xxxxyy_1[i] * pc_x[i];

        ta1_y_xxyy_xxxxyz_0[i] = ta1_y_yy_xxxxyz_0[i] * fe_0 - ta1_y_yy_xxxxyz_1[i] * fe_0 + 4.0 * ta1_y_xyy_xxxyz_0[i] * fe_0 -
                                 4.0 * ta1_y_xyy_xxxyz_1[i] * fe_0 + ta1_y_xyy_xxxxyz_0[i] * pa_x[i] - ta1_y_xyy_xxxxyz_1[i] * pc_x[i];

        ta1_y_xxyy_xxxxzz_0[i] = ta1_y_xx_xxxxzz_0[i] * fe_0 - ta1_y_xx_xxxxzz_1[i] * fe_0 + ta_xxy_xxxxzz_1[i] + ta1_y_xxy_xxxxzz_0[i] * pa_y[i] -
                                 ta1_y_xxy_xxxxzz_1[i] * pc_y[i];

        ta1_y_xxyy_xxxyyy_0[i] = ta1_y_yy_xxxyyy_0[i] * fe_0 - ta1_y_yy_xxxyyy_1[i] * fe_0 + 3.0 * ta1_y_xyy_xxyyy_0[i] * fe_0 -
                                 3.0 * ta1_y_xyy_xxyyy_1[i] * fe_0 + ta1_y_xyy_xxxyyy_0[i] * pa_x[i] - ta1_y_xyy_xxxyyy_1[i] * pc_x[i];

        ta1_y_xxyy_xxxyyz_0[i] = ta1_y_yy_xxxyyz_0[i] * fe_0 - ta1_y_yy_xxxyyz_1[i] * fe_0 + 3.0 * ta1_y_xyy_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_y_xyy_xxyyz_1[i] * fe_0 + ta1_y_xyy_xxxyyz_0[i] * pa_x[i] - ta1_y_xyy_xxxyyz_1[i] * pc_x[i];

        ta1_y_xxyy_xxxyzz_0[i] = ta1_y_yy_xxxyzz_0[i] * fe_0 - ta1_y_yy_xxxyzz_1[i] * fe_0 + 3.0 * ta1_y_xyy_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_xyy_xxyzz_1[i] * fe_0 + ta1_y_xyy_xxxyzz_0[i] * pa_x[i] - ta1_y_xyy_xxxyzz_1[i] * pc_x[i];

        ta1_y_xxyy_xxxzzz_0[i] = ta1_y_xx_xxxzzz_0[i] * fe_0 - ta1_y_xx_xxxzzz_1[i] * fe_0 + ta_xxy_xxxzzz_1[i] + ta1_y_xxy_xxxzzz_0[i] * pa_y[i] -
                                 ta1_y_xxy_xxxzzz_1[i] * pc_y[i];

        ta1_y_xxyy_xxyyyy_0[i] = ta1_y_yy_xxyyyy_0[i] * fe_0 - ta1_y_yy_xxyyyy_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyyyy_0[i] * fe_0 -
                                 2.0 * ta1_y_xyy_xyyyy_1[i] * fe_0 + ta1_y_xyy_xxyyyy_0[i] * pa_x[i] - ta1_y_xyy_xxyyyy_1[i] * pc_x[i];

        ta1_y_xxyy_xxyyyz_0[i] = ta1_y_yy_xxyyyz_0[i] * fe_0 - ta1_y_yy_xxyyyz_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_xyy_xyyyz_1[i] * fe_0 + ta1_y_xyy_xxyyyz_0[i] * pa_x[i] - ta1_y_xyy_xxyyyz_1[i] * pc_x[i];

        ta1_y_xxyy_xxyyzz_0[i] = ta1_y_yy_xxyyzz_0[i] * fe_0 - ta1_y_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyyzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xyy_xyyzz_1[i] * fe_0 + ta1_y_xyy_xxyyzz_0[i] * pa_x[i] - ta1_y_xyy_xxyyzz_1[i] * pc_x[i];

        ta1_y_xxyy_xxyzzz_0[i] = ta1_y_yy_xxyzzz_0[i] * fe_0 - ta1_y_yy_xxyzzz_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xyy_xyzzz_1[i] * fe_0 + ta1_y_xyy_xxyzzz_0[i] * pa_x[i] - ta1_y_xyy_xxyzzz_1[i] * pc_x[i];

        ta1_y_xxyy_xxzzzz_0[i] = ta1_y_xx_xxzzzz_0[i] * fe_0 - ta1_y_xx_xxzzzz_1[i] * fe_0 + ta_xxy_xxzzzz_1[i] + ta1_y_xxy_xxzzzz_0[i] * pa_y[i] -
                                 ta1_y_xxy_xxzzzz_1[i] * pc_y[i];

        ta1_y_xxyy_xyyyyy_0[i] = ta1_y_yy_xyyyyy_0[i] * fe_0 - ta1_y_yy_xyyyyy_1[i] * fe_0 + ta1_y_xyy_yyyyy_0[i] * fe_0 -
                                 ta1_y_xyy_yyyyy_1[i] * fe_0 + ta1_y_xyy_xyyyyy_0[i] * pa_x[i] - ta1_y_xyy_xyyyyy_1[i] * pc_x[i];

        ta1_y_xxyy_xyyyyz_0[i] = ta1_y_yy_xyyyyz_0[i] * fe_0 - ta1_y_yy_xyyyyz_1[i] * fe_0 + ta1_y_xyy_yyyyz_0[i] * fe_0 -
                                 ta1_y_xyy_yyyyz_1[i] * fe_0 + ta1_y_xyy_xyyyyz_0[i] * pa_x[i] - ta1_y_xyy_xyyyyz_1[i] * pc_x[i];

        ta1_y_xxyy_xyyyzz_0[i] = ta1_y_yy_xyyyzz_0[i] * fe_0 - ta1_y_yy_xyyyzz_1[i] * fe_0 + ta1_y_xyy_yyyzz_0[i] * fe_0 -
                                 ta1_y_xyy_yyyzz_1[i] * fe_0 + ta1_y_xyy_xyyyzz_0[i] * pa_x[i] - ta1_y_xyy_xyyyzz_1[i] * pc_x[i];

        ta1_y_xxyy_xyyzzz_0[i] = ta1_y_yy_xyyzzz_0[i] * fe_0 - ta1_y_yy_xyyzzz_1[i] * fe_0 + ta1_y_xyy_yyzzz_0[i] * fe_0 -
                                 ta1_y_xyy_yyzzz_1[i] * fe_0 + ta1_y_xyy_xyyzzz_0[i] * pa_x[i] - ta1_y_xyy_xyyzzz_1[i] * pc_x[i];

        ta1_y_xxyy_xyzzzz_0[i] = ta1_y_yy_xyzzzz_0[i] * fe_0 - ta1_y_yy_xyzzzz_1[i] * fe_0 + ta1_y_xyy_yzzzz_0[i] * fe_0 -
                                 ta1_y_xyy_yzzzz_1[i] * fe_0 + ta1_y_xyy_xyzzzz_0[i] * pa_x[i] - ta1_y_xyy_xyzzzz_1[i] * pc_x[i];

        ta1_y_xxyy_xzzzzz_0[i] = ta1_y_xx_xzzzzz_0[i] * fe_0 - ta1_y_xx_xzzzzz_1[i] * fe_0 + ta_xxy_xzzzzz_1[i] + ta1_y_xxy_xzzzzz_0[i] * pa_y[i] -
                                 ta1_y_xxy_xzzzzz_1[i] * pc_y[i];

        ta1_y_xxyy_yyyyyy_0[i] =
            ta1_y_yy_yyyyyy_0[i] * fe_0 - ta1_y_yy_yyyyyy_1[i] * fe_0 + ta1_y_xyy_yyyyyy_0[i] * pa_x[i] - ta1_y_xyy_yyyyyy_1[i] * pc_x[i];

        ta1_y_xxyy_yyyyyz_0[i] =
            ta1_y_yy_yyyyyz_0[i] * fe_0 - ta1_y_yy_yyyyyz_1[i] * fe_0 + ta1_y_xyy_yyyyyz_0[i] * pa_x[i] - ta1_y_xyy_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxyy_yyyyzz_0[i] =
            ta1_y_yy_yyyyzz_0[i] * fe_0 - ta1_y_yy_yyyyzz_1[i] * fe_0 + ta1_y_xyy_yyyyzz_0[i] * pa_x[i] - ta1_y_xyy_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxyy_yyyzzz_0[i] =
            ta1_y_yy_yyyzzz_0[i] * fe_0 - ta1_y_yy_yyyzzz_1[i] * fe_0 + ta1_y_xyy_yyyzzz_0[i] * pa_x[i] - ta1_y_xyy_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxyy_yyzzzz_0[i] =
            ta1_y_yy_yyzzzz_0[i] * fe_0 - ta1_y_yy_yyzzzz_1[i] * fe_0 + ta1_y_xyy_yyzzzz_0[i] * pa_x[i] - ta1_y_xyy_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxyy_yzzzzz_0[i] =
            ta1_y_yy_yzzzzz_0[i] * fe_0 - ta1_y_yy_yzzzzz_1[i] * fe_0 + ta1_y_xyy_yzzzzz_0[i] * pa_x[i] - ta1_y_xyy_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxyy_zzzzzz_0[i] =
            ta1_y_yy_zzzzzz_0[i] * fe_0 - ta1_y_yy_zzzzzz_1[i] * fe_0 + ta1_y_xyy_zzzzzz_0[i] * pa_x[i] - ta1_y_xyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 532-560 components of targeted buffer : GI

    auto ta1_y_xxyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 532);

    auto ta1_y_xxyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 533);

    auto ta1_y_xxyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 534);

    auto ta1_y_xxyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 535);

    auto ta1_y_xxyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 536);

    auto ta1_y_xxyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 537);

    auto ta1_y_xxyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 538);

    auto ta1_y_xxyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 539);

    auto ta1_y_xxyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 540);

    auto ta1_y_xxyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 541);

    auto ta1_y_xxyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 542);

    auto ta1_y_xxyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 543);

    auto ta1_y_xxyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 544);

    auto ta1_y_xxyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 545);

    auto ta1_y_xxyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 546);

    auto ta1_y_xxyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 547);

    auto ta1_y_xxyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 548);

    auto ta1_y_xxyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 549);

    auto ta1_y_xxyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 550);

    auto ta1_y_xxyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 551);

    auto ta1_y_xxyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 552);

    auto ta1_y_xxyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 553);

    auto ta1_y_xxyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 554);

    auto ta1_y_xxyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 555);

    auto ta1_y_xxyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 556);

    auto ta1_y_xxyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 557);

    auto ta1_y_xxyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 558);

    auto ta1_y_xxyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 559);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pa_z,                \
                             pc_x,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_y_xxy_xxxxxx_0,  \
                             ta1_y_xxy_xxxxxx_1,  \
                             ta1_y_xxy_xxxxxy_0,  \
                             ta1_y_xxy_xxxxxy_1,  \
                             ta1_y_xxy_xxxxy_0,   \
                             ta1_y_xxy_xxxxy_1,   \
                             ta1_y_xxy_xxxxyy_0,  \
                             ta1_y_xxy_xxxxyy_1,  \
                             ta1_y_xxy_xxxxyz_0,  \
                             ta1_y_xxy_xxxxyz_1,  \
                             ta1_y_xxy_xxxyy_0,   \
                             ta1_y_xxy_xxxyy_1,   \
                             ta1_y_xxy_xxxyyy_0,  \
                             ta1_y_xxy_xxxyyy_1,  \
                             ta1_y_xxy_xxxyyz_0,  \
                             ta1_y_xxy_xxxyyz_1,  \
                             ta1_y_xxy_xxxyz_0,   \
                             ta1_y_xxy_xxxyz_1,   \
                             ta1_y_xxy_xxxyzz_0,  \
                             ta1_y_xxy_xxxyzz_1,  \
                             ta1_y_xxy_xxyyy_0,   \
                             ta1_y_xxy_xxyyy_1,   \
                             ta1_y_xxy_xxyyyy_0,  \
                             ta1_y_xxy_xxyyyy_1,  \
                             ta1_y_xxy_xxyyyz_0,  \
                             ta1_y_xxy_xxyyyz_1,  \
                             ta1_y_xxy_xxyyz_0,   \
                             ta1_y_xxy_xxyyz_1,   \
                             ta1_y_xxy_xxyyzz_0,  \
                             ta1_y_xxy_xxyyzz_1,  \
                             ta1_y_xxy_xxyzz_0,   \
                             ta1_y_xxy_xxyzz_1,   \
                             ta1_y_xxy_xxyzzz_0,  \
                             ta1_y_xxy_xxyzzz_1,  \
                             ta1_y_xxy_xyyyy_0,   \
                             ta1_y_xxy_xyyyy_1,   \
                             ta1_y_xxy_xyyyyy_0,  \
                             ta1_y_xxy_xyyyyy_1,  \
                             ta1_y_xxy_xyyyyz_0,  \
                             ta1_y_xxy_xyyyyz_1,  \
                             ta1_y_xxy_xyyyz_0,   \
                             ta1_y_xxy_xyyyz_1,   \
                             ta1_y_xxy_xyyyzz_0,  \
                             ta1_y_xxy_xyyyzz_1,  \
                             ta1_y_xxy_xyyzz_0,   \
                             ta1_y_xxy_xyyzz_1,   \
                             ta1_y_xxy_xyyzzz_0,  \
                             ta1_y_xxy_xyyzzz_1,  \
                             ta1_y_xxy_xyzzz_0,   \
                             ta1_y_xxy_xyzzz_1,   \
                             ta1_y_xxy_xyzzzz_0,  \
                             ta1_y_xxy_xyzzzz_1,  \
                             ta1_y_xxy_yyyyyy_0,  \
                             ta1_y_xxy_yyyyyy_1,  \
                             ta1_y_xxyz_xxxxxx_0, \
                             ta1_y_xxyz_xxxxxy_0, \
                             ta1_y_xxyz_xxxxxz_0, \
                             ta1_y_xxyz_xxxxyy_0, \
                             ta1_y_xxyz_xxxxyz_0, \
                             ta1_y_xxyz_xxxxzz_0, \
                             ta1_y_xxyz_xxxyyy_0, \
                             ta1_y_xxyz_xxxyyz_0, \
                             ta1_y_xxyz_xxxyzz_0, \
                             ta1_y_xxyz_xxxzzz_0, \
                             ta1_y_xxyz_xxyyyy_0, \
                             ta1_y_xxyz_xxyyyz_0, \
                             ta1_y_xxyz_xxyyzz_0, \
                             ta1_y_xxyz_xxyzzz_0, \
                             ta1_y_xxyz_xxzzzz_0, \
                             ta1_y_xxyz_xyyyyy_0, \
                             ta1_y_xxyz_xyyyyz_0, \
                             ta1_y_xxyz_xyyyzz_0, \
                             ta1_y_xxyz_xyyzzz_0, \
                             ta1_y_xxyz_xyzzzz_0, \
                             ta1_y_xxyz_xzzzzz_0, \
                             ta1_y_xxyz_yyyyyy_0, \
                             ta1_y_xxyz_yyyyyz_0, \
                             ta1_y_xxyz_yyyyzz_0, \
                             ta1_y_xxyz_yyyzzz_0, \
                             ta1_y_xxyz_yyzzzz_0, \
                             ta1_y_xxyz_yzzzzz_0, \
                             ta1_y_xxyz_zzzzzz_0, \
                             ta1_y_xxz_xxxxxz_0,  \
                             ta1_y_xxz_xxxxxz_1,  \
                             ta1_y_xxz_xxxxzz_0,  \
                             ta1_y_xxz_xxxxzz_1,  \
                             ta1_y_xxz_xxxzzz_0,  \
                             ta1_y_xxz_xxxzzz_1,  \
                             ta1_y_xxz_xxzzzz_0,  \
                             ta1_y_xxz_xxzzzz_1,  \
                             ta1_y_xxz_xzzzzz_0,  \
                             ta1_y_xxz_xzzzzz_1,  \
                             ta1_y_xxz_zzzzzz_0,  \
                             ta1_y_xxz_zzzzzz_1,  \
                             ta1_y_xyz_yyyyyz_0,  \
                             ta1_y_xyz_yyyyyz_1,  \
                             ta1_y_xyz_yyyyzz_0,  \
                             ta1_y_xyz_yyyyzz_1,  \
                             ta1_y_xyz_yyyzzz_0,  \
                             ta1_y_xyz_yyyzzz_1,  \
                             ta1_y_xyz_yyzzzz_0,  \
                             ta1_y_xyz_yyzzzz_1,  \
                             ta1_y_xyz_yzzzzz_0,  \
                             ta1_y_xyz_yzzzzz_1,  \
                             ta1_y_yz_yyyyyz_0,   \
                             ta1_y_yz_yyyyyz_1,   \
                             ta1_y_yz_yyyyzz_0,   \
                             ta1_y_yz_yyyyzz_1,   \
                             ta1_y_yz_yyyzzz_0,   \
                             ta1_y_yz_yyyzzz_1,   \
                             ta1_y_yz_yyzzzz_0,   \
                             ta1_y_yz_yyzzzz_1,   \
                             ta1_y_yz_yzzzzz_0,   \
                             ta1_y_yz_yzzzzz_1,   \
                             ta_xxz_xxxxxz_1,     \
                             ta_xxz_xxxxzz_1,     \
                             ta_xxz_xxxzzz_1,     \
                             ta_xxz_xxzzzz_1,     \
                             ta_xxz_xzzzzz_1,     \
                             ta_xxz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyz_xxxxxx_0[i] = ta1_y_xxy_xxxxxx_0[i] * pa_z[i] - ta1_y_xxy_xxxxxx_1[i] * pc_z[i];

        ta1_y_xxyz_xxxxxy_0[i] = ta1_y_xxy_xxxxxy_0[i] * pa_z[i] - ta1_y_xxy_xxxxxy_1[i] * pc_z[i];

        ta1_y_xxyz_xxxxxz_0[i] = ta_xxz_xxxxxz_1[i] + ta1_y_xxz_xxxxxz_0[i] * pa_y[i] - ta1_y_xxz_xxxxxz_1[i] * pc_y[i];

        ta1_y_xxyz_xxxxyy_0[i] = ta1_y_xxy_xxxxyy_0[i] * pa_z[i] - ta1_y_xxy_xxxxyy_1[i] * pc_z[i];

        ta1_y_xxyz_xxxxyz_0[i] =
            ta1_y_xxy_xxxxy_0[i] * fe_0 - ta1_y_xxy_xxxxy_1[i] * fe_0 + ta1_y_xxy_xxxxyz_0[i] * pa_z[i] - ta1_y_xxy_xxxxyz_1[i] * pc_z[i];

        ta1_y_xxyz_xxxxzz_0[i] = ta_xxz_xxxxzz_1[i] + ta1_y_xxz_xxxxzz_0[i] * pa_y[i] - ta1_y_xxz_xxxxzz_1[i] * pc_y[i];

        ta1_y_xxyz_xxxyyy_0[i] = ta1_y_xxy_xxxyyy_0[i] * pa_z[i] - ta1_y_xxy_xxxyyy_1[i] * pc_z[i];

        ta1_y_xxyz_xxxyyz_0[i] =
            ta1_y_xxy_xxxyy_0[i] * fe_0 - ta1_y_xxy_xxxyy_1[i] * fe_0 + ta1_y_xxy_xxxyyz_0[i] * pa_z[i] - ta1_y_xxy_xxxyyz_1[i] * pc_z[i];

        ta1_y_xxyz_xxxyzz_0[i] =
            2.0 * ta1_y_xxy_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxxyz_1[i] * fe_0 + ta1_y_xxy_xxxyzz_0[i] * pa_z[i] - ta1_y_xxy_xxxyzz_1[i] * pc_z[i];

        ta1_y_xxyz_xxxzzz_0[i] = ta_xxz_xxxzzz_1[i] + ta1_y_xxz_xxxzzz_0[i] * pa_y[i] - ta1_y_xxz_xxxzzz_1[i] * pc_y[i];

        ta1_y_xxyz_xxyyyy_0[i] = ta1_y_xxy_xxyyyy_0[i] * pa_z[i] - ta1_y_xxy_xxyyyy_1[i] * pc_z[i];

        ta1_y_xxyz_xxyyyz_0[i] =
            ta1_y_xxy_xxyyy_0[i] * fe_0 - ta1_y_xxy_xxyyy_1[i] * fe_0 + ta1_y_xxy_xxyyyz_0[i] * pa_z[i] - ta1_y_xxy_xxyyyz_1[i] * pc_z[i];

        ta1_y_xxyz_xxyyzz_0[i] =
            2.0 * ta1_y_xxy_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxyyz_1[i] * fe_0 + ta1_y_xxy_xxyyzz_0[i] * pa_z[i] - ta1_y_xxy_xxyyzz_1[i] * pc_z[i];

        ta1_y_xxyz_xxyzzz_0[i] =
            3.0 * ta1_y_xxy_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_xxy_xxyzz_1[i] * fe_0 + ta1_y_xxy_xxyzzz_0[i] * pa_z[i] - ta1_y_xxy_xxyzzz_1[i] * pc_z[i];

        ta1_y_xxyz_xxzzzz_0[i] = ta_xxz_xxzzzz_1[i] + ta1_y_xxz_xxzzzz_0[i] * pa_y[i] - ta1_y_xxz_xxzzzz_1[i] * pc_y[i];

        ta1_y_xxyz_xyyyyy_0[i] = ta1_y_xxy_xyyyyy_0[i] * pa_z[i] - ta1_y_xxy_xyyyyy_1[i] * pc_z[i];

        ta1_y_xxyz_xyyyyz_0[i] =
            ta1_y_xxy_xyyyy_0[i] * fe_0 - ta1_y_xxy_xyyyy_1[i] * fe_0 + ta1_y_xxy_xyyyyz_0[i] * pa_z[i] - ta1_y_xxy_xyyyyz_1[i] * pc_z[i];

        ta1_y_xxyz_xyyyzz_0[i] =
            2.0 * ta1_y_xxy_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xyyyz_1[i] * fe_0 + ta1_y_xxy_xyyyzz_0[i] * pa_z[i] - ta1_y_xxy_xyyyzz_1[i] * pc_z[i];

        ta1_y_xxyz_xyyzzz_0[i] =
            3.0 * ta1_y_xxy_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_xxy_xyyzz_1[i] * fe_0 + ta1_y_xxy_xyyzzz_0[i] * pa_z[i] - ta1_y_xxy_xyyzzz_1[i] * pc_z[i];

        ta1_y_xxyz_xyzzzz_0[i] =
            4.0 * ta1_y_xxy_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_xxy_xyzzz_1[i] * fe_0 + ta1_y_xxy_xyzzzz_0[i] * pa_z[i] - ta1_y_xxy_xyzzzz_1[i] * pc_z[i];

        ta1_y_xxyz_xzzzzz_0[i] = ta_xxz_xzzzzz_1[i] + ta1_y_xxz_xzzzzz_0[i] * pa_y[i] - ta1_y_xxz_xzzzzz_1[i] * pc_y[i];

        ta1_y_xxyz_yyyyyy_0[i] = ta1_y_xxy_yyyyyy_0[i] * pa_z[i] - ta1_y_xxy_yyyyyy_1[i] * pc_z[i];

        ta1_y_xxyz_yyyyyz_0[i] =
            ta1_y_yz_yyyyyz_0[i] * fe_0 - ta1_y_yz_yyyyyz_1[i] * fe_0 + ta1_y_xyz_yyyyyz_0[i] * pa_x[i] - ta1_y_xyz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxyz_yyyyzz_0[i] =
            ta1_y_yz_yyyyzz_0[i] * fe_0 - ta1_y_yz_yyyyzz_1[i] * fe_0 + ta1_y_xyz_yyyyzz_0[i] * pa_x[i] - ta1_y_xyz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxyz_yyyzzz_0[i] =
            ta1_y_yz_yyyzzz_0[i] * fe_0 - ta1_y_yz_yyyzzz_1[i] * fe_0 + ta1_y_xyz_yyyzzz_0[i] * pa_x[i] - ta1_y_xyz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxyz_yyzzzz_0[i] =
            ta1_y_yz_yyzzzz_0[i] * fe_0 - ta1_y_yz_yyzzzz_1[i] * fe_0 + ta1_y_xyz_yyzzzz_0[i] * pa_x[i] - ta1_y_xyz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxyz_yzzzzz_0[i] =
            ta1_y_yz_yzzzzz_0[i] * fe_0 - ta1_y_yz_yzzzzz_1[i] * fe_0 + ta1_y_xyz_yzzzzz_0[i] * pa_x[i] - ta1_y_xyz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxyz_zzzzzz_0[i] = ta_xxz_zzzzzz_1[i] + ta1_y_xxz_zzzzzz_0[i] * pa_y[i] - ta1_y_xxz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 560-588 components of targeted buffer : GI

    auto ta1_y_xxzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 560);

    auto ta1_y_xxzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 561);

    auto ta1_y_xxzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 562);

    auto ta1_y_xxzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 563);

    auto ta1_y_xxzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 564);

    auto ta1_y_xxzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 565);

    auto ta1_y_xxzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 566);

    auto ta1_y_xxzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 567);

    auto ta1_y_xxzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 568);

    auto ta1_y_xxzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 569);

    auto ta1_y_xxzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 570);

    auto ta1_y_xxzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 571);

    auto ta1_y_xxzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 572);

    auto ta1_y_xxzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 573);

    auto ta1_y_xxzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 574);

    auto ta1_y_xxzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 575);

    auto ta1_y_xxzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 576);

    auto ta1_y_xxzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 577);

    auto ta1_y_xxzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 578);

    auto ta1_y_xxzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 579);

    auto ta1_y_xxzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 580);

    auto ta1_y_xxzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 581);

    auto ta1_y_xxzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 582);

    auto ta1_y_xxzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 583);

    auto ta1_y_xxzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 584);

    auto ta1_y_xxzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 585);

    auto ta1_y_xxzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 586);

    auto ta1_y_xxzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 587);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_y_xx_xxxxxx_0,   \
                             ta1_y_xx_xxxxxx_1,   \
                             ta1_y_xx_xxxxxy_0,   \
                             ta1_y_xx_xxxxxy_1,   \
                             ta1_y_xx_xxxxyy_0,   \
                             ta1_y_xx_xxxxyy_1,   \
                             ta1_y_xx_xxxyyy_0,   \
                             ta1_y_xx_xxxyyy_1,   \
                             ta1_y_xx_xxyyyy_0,   \
                             ta1_y_xx_xxyyyy_1,   \
                             ta1_y_xx_xyyyyy_0,   \
                             ta1_y_xx_xyyyyy_1,   \
                             ta1_y_xxz_xxxxxx_0,  \
                             ta1_y_xxz_xxxxxx_1,  \
                             ta1_y_xxz_xxxxxy_0,  \
                             ta1_y_xxz_xxxxxy_1,  \
                             ta1_y_xxz_xxxxyy_0,  \
                             ta1_y_xxz_xxxxyy_1,  \
                             ta1_y_xxz_xxxyyy_0,  \
                             ta1_y_xxz_xxxyyy_1,  \
                             ta1_y_xxz_xxyyyy_0,  \
                             ta1_y_xxz_xxyyyy_1,  \
                             ta1_y_xxz_xyyyyy_0,  \
                             ta1_y_xxz_xyyyyy_1,  \
                             ta1_y_xxzz_xxxxxx_0, \
                             ta1_y_xxzz_xxxxxy_0, \
                             ta1_y_xxzz_xxxxxz_0, \
                             ta1_y_xxzz_xxxxyy_0, \
                             ta1_y_xxzz_xxxxyz_0, \
                             ta1_y_xxzz_xxxxzz_0, \
                             ta1_y_xxzz_xxxyyy_0, \
                             ta1_y_xxzz_xxxyyz_0, \
                             ta1_y_xxzz_xxxyzz_0, \
                             ta1_y_xxzz_xxxzzz_0, \
                             ta1_y_xxzz_xxyyyy_0, \
                             ta1_y_xxzz_xxyyyz_0, \
                             ta1_y_xxzz_xxyyzz_0, \
                             ta1_y_xxzz_xxyzzz_0, \
                             ta1_y_xxzz_xxzzzz_0, \
                             ta1_y_xxzz_xyyyyy_0, \
                             ta1_y_xxzz_xyyyyz_0, \
                             ta1_y_xxzz_xyyyzz_0, \
                             ta1_y_xxzz_xyyzzz_0, \
                             ta1_y_xxzz_xyzzzz_0, \
                             ta1_y_xxzz_xzzzzz_0, \
                             ta1_y_xxzz_yyyyyy_0, \
                             ta1_y_xxzz_yyyyyz_0, \
                             ta1_y_xxzz_yyyyzz_0, \
                             ta1_y_xxzz_yyyzzz_0, \
                             ta1_y_xxzz_yyzzzz_0, \
                             ta1_y_xxzz_yzzzzz_0, \
                             ta1_y_xxzz_zzzzzz_0, \
                             ta1_y_xzz_xxxxxz_0,  \
                             ta1_y_xzz_xxxxxz_1,  \
                             ta1_y_xzz_xxxxyz_0,  \
                             ta1_y_xzz_xxxxyz_1,  \
                             ta1_y_xzz_xxxxz_0,   \
                             ta1_y_xzz_xxxxz_1,   \
                             ta1_y_xzz_xxxxzz_0,  \
                             ta1_y_xzz_xxxxzz_1,  \
                             ta1_y_xzz_xxxyyz_0,  \
                             ta1_y_xzz_xxxyyz_1,  \
                             ta1_y_xzz_xxxyz_0,   \
                             ta1_y_xzz_xxxyz_1,   \
                             ta1_y_xzz_xxxyzz_0,  \
                             ta1_y_xzz_xxxyzz_1,  \
                             ta1_y_xzz_xxxzz_0,   \
                             ta1_y_xzz_xxxzz_1,   \
                             ta1_y_xzz_xxxzzz_0,  \
                             ta1_y_xzz_xxxzzz_1,  \
                             ta1_y_xzz_xxyyyz_0,  \
                             ta1_y_xzz_xxyyyz_1,  \
                             ta1_y_xzz_xxyyz_0,   \
                             ta1_y_xzz_xxyyz_1,   \
                             ta1_y_xzz_xxyyzz_0,  \
                             ta1_y_xzz_xxyyzz_1,  \
                             ta1_y_xzz_xxyzz_0,   \
                             ta1_y_xzz_xxyzz_1,   \
                             ta1_y_xzz_xxyzzz_0,  \
                             ta1_y_xzz_xxyzzz_1,  \
                             ta1_y_xzz_xxzzz_0,   \
                             ta1_y_xzz_xxzzz_1,   \
                             ta1_y_xzz_xxzzzz_0,  \
                             ta1_y_xzz_xxzzzz_1,  \
                             ta1_y_xzz_xyyyyz_0,  \
                             ta1_y_xzz_xyyyyz_1,  \
                             ta1_y_xzz_xyyyz_0,   \
                             ta1_y_xzz_xyyyz_1,   \
                             ta1_y_xzz_xyyyzz_0,  \
                             ta1_y_xzz_xyyyzz_1,  \
                             ta1_y_xzz_xyyzz_0,   \
                             ta1_y_xzz_xyyzz_1,   \
                             ta1_y_xzz_xyyzzz_0,  \
                             ta1_y_xzz_xyyzzz_1,  \
                             ta1_y_xzz_xyzzz_0,   \
                             ta1_y_xzz_xyzzz_1,   \
                             ta1_y_xzz_xyzzzz_0,  \
                             ta1_y_xzz_xyzzzz_1,  \
                             ta1_y_xzz_xzzzz_0,   \
                             ta1_y_xzz_xzzzz_1,   \
                             ta1_y_xzz_xzzzzz_0,  \
                             ta1_y_xzz_xzzzzz_1,  \
                             ta1_y_xzz_yyyyyy_0,  \
                             ta1_y_xzz_yyyyyy_1,  \
                             ta1_y_xzz_yyyyyz_0,  \
                             ta1_y_xzz_yyyyyz_1,  \
                             ta1_y_xzz_yyyyz_0,   \
                             ta1_y_xzz_yyyyz_1,   \
                             ta1_y_xzz_yyyyzz_0,  \
                             ta1_y_xzz_yyyyzz_1,  \
                             ta1_y_xzz_yyyzz_0,   \
                             ta1_y_xzz_yyyzz_1,   \
                             ta1_y_xzz_yyyzzz_0,  \
                             ta1_y_xzz_yyyzzz_1,  \
                             ta1_y_xzz_yyzzz_0,   \
                             ta1_y_xzz_yyzzz_1,   \
                             ta1_y_xzz_yyzzzz_0,  \
                             ta1_y_xzz_yyzzzz_1,  \
                             ta1_y_xzz_yzzzz_0,   \
                             ta1_y_xzz_yzzzz_1,   \
                             ta1_y_xzz_yzzzzz_0,  \
                             ta1_y_xzz_yzzzzz_1,  \
                             ta1_y_xzz_zzzzz_0,   \
                             ta1_y_xzz_zzzzz_1,   \
                             ta1_y_xzz_zzzzzz_0,  \
                             ta1_y_xzz_zzzzzz_1,  \
                             ta1_y_zz_xxxxxz_0,   \
                             ta1_y_zz_xxxxxz_1,   \
                             ta1_y_zz_xxxxyz_0,   \
                             ta1_y_zz_xxxxyz_1,   \
                             ta1_y_zz_xxxxzz_0,   \
                             ta1_y_zz_xxxxzz_1,   \
                             ta1_y_zz_xxxyyz_0,   \
                             ta1_y_zz_xxxyyz_1,   \
                             ta1_y_zz_xxxyzz_0,   \
                             ta1_y_zz_xxxyzz_1,   \
                             ta1_y_zz_xxxzzz_0,   \
                             ta1_y_zz_xxxzzz_1,   \
                             ta1_y_zz_xxyyyz_0,   \
                             ta1_y_zz_xxyyyz_1,   \
                             ta1_y_zz_xxyyzz_0,   \
                             ta1_y_zz_xxyyzz_1,   \
                             ta1_y_zz_xxyzzz_0,   \
                             ta1_y_zz_xxyzzz_1,   \
                             ta1_y_zz_xxzzzz_0,   \
                             ta1_y_zz_xxzzzz_1,   \
                             ta1_y_zz_xyyyyz_0,   \
                             ta1_y_zz_xyyyyz_1,   \
                             ta1_y_zz_xyyyzz_0,   \
                             ta1_y_zz_xyyyzz_1,   \
                             ta1_y_zz_xyyzzz_0,   \
                             ta1_y_zz_xyyzzz_1,   \
                             ta1_y_zz_xyzzzz_0,   \
                             ta1_y_zz_xyzzzz_1,   \
                             ta1_y_zz_xzzzzz_0,   \
                             ta1_y_zz_xzzzzz_1,   \
                             ta1_y_zz_yyyyyy_0,   \
                             ta1_y_zz_yyyyyy_1,   \
                             ta1_y_zz_yyyyyz_0,   \
                             ta1_y_zz_yyyyyz_1,   \
                             ta1_y_zz_yyyyzz_0,   \
                             ta1_y_zz_yyyyzz_1,   \
                             ta1_y_zz_yyyzzz_0,   \
                             ta1_y_zz_yyyzzz_1,   \
                             ta1_y_zz_yyzzzz_0,   \
                             ta1_y_zz_yyzzzz_1,   \
                             ta1_y_zz_yzzzzz_0,   \
                             ta1_y_zz_yzzzzz_1,   \
                             ta1_y_zz_zzzzzz_0,   \
                             ta1_y_zz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzz_xxxxxx_0[i] =
            ta1_y_xx_xxxxxx_0[i] * fe_0 - ta1_y_xx_xxxxxx_1[i] * fe_0 + ta1_y_xxz_xxxxxx_0[i] * pa_z[i] - ta1_y_xxz_xxxxxx_1[i] * pc_z[i];

        ta1_y_xxzz_xxxxxy_0[i] =
            ta1_y_xx_xxxxxy_0[i] * fe_0 - ta1_y_xx_xxxxxy_1[i] * fe_0 + ta1_y_xxz_xxxxxy_0[i] * pa_z[i] - ta1_y_xxz_xxxxxy_1[i] * pc_z[i];

        ta1_y_xxzz_xxxxxz_0[i] = ta1_y_zz_xxxxxz_0[i] * fe_0 - ta1_y_zz_xxxxxz_1[i] * fe_0 + 5.0 * ta1_y_xzz_xxxxz_0[i] * fe_0 -
                                 5.0 * ta1_y_xzz_xxxxz_1[i] * fe_0 + ta1_y_xzz_xxxxxz_0[i] * pa_x[i] - ta1_y_xzz_xxxxxz_1[i] * pc_x[i];

        ta1_y_xxzz_xxxxyy_0[i] =
            ta1_y_xx_xxxxyy_0[i] * fe_0 - ta1_y_xx_xxxxyy_1[i] * fe_0 + ta1_y_xxz_xxxxyy_0[i] * pa_z[i] - ta1_y_xxz_xxxxyy_1[i] * pc_z[i];

        ta1_y_xxzz_xxxxyz_0[i] = ta1_y_zz_xxxxyz_0[i] * fe_0 - ta1_y_zz_xxxxyz_1[i] * fe_0 + 4.0 * ta1_y_xzz_xxxyz_0[i] * fe_0 -
                                 4.0 * ta1_y_xzz_xxxyz_1[i] * fe_0 + ta1_y_xzz_xxxxyz_0[i] * pa_x[i] - ta1_y_xzz_xxxxyz_1[i] * pc_x[i];

        ta1_y_xxzz_xxxxzz_0[i] = ta1_y_zz_xxxxzz_0[i] * fe_0 - ta1_y_zz_xxxxzz_1[i] * fe_0 + 4.0 * ta1_y_xzz_xxxzz_0[i] * fe_0 -
                                 4.0 * ta1_y_xzz_xxxzz_1[i] * fe_0 + ta1_y_xzz_xxxxzz_0[i] * pa_x[i] - ta1_y_xzz_xxxxzz_1[i] * pc_x[i];

        ta1_y_xxzz_xxxyyy_0[i] =
            ta1_y_xx_xxxyyy_0[i] * fe_0 - ta1_y_xx_xxxyyy_1[i] * fe_0 + ta1_y_xxz_xxxyyy_0[i] * pa_z[i] - ta1_y_xxz_xxxyyy_1[i] * pc_z[i];

        ta1_y_xxzz_xxxyyz_0[i] = ta1_y_zz_xxxyyz_0[i] * fe_0 - ta1_y_zz_xxxyyz_1[i] * fe_0 + 3.0 * ta1_y_xzz_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_y_xzz_xxyyz_1[i] * fe_0 + ta1_y_xzz_xxxyyz_0[i] * pa_x[i] - ta1_y_xzz_xxxyyz_1[i] * pc_x[i];

        ta1_y_xxzz_xxxyzz_0[i] = ta1_y_zz_xxxyzz_0[i] * fe_0 - ta1_y_zz_xxxyzz_1[i] * fe_0 + 3.0 * ta1_y_xzz_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_xzz_xxyzz_1[i] * fe_0 + ta1_y_xzz_xxxyzz_0[i] * pa_x[i] - ta1_y_xzz_xxxyzz_1[i] * pc_x[i];

        ta1_y_xxzz_xxxzzz_0[i] = ta1_y_zz_xxxzzz_0[i] * fe_0 - ta1_y_zz_xxxzzz_1[i] * fe_0 + 3.0 * ta1_y_xzz_xxzzz_0[i] * fe_0 -
                                 3.0 * ta1_y_xzz_xxzzz_1[i] * fe_0 + ta1_y_xzz_xxxzzz_0[i] * pa_x[i] - ta1_y_xzz_xxxzzz_1[i] * pc_x[i];

        ta1_y_xxzz_xxyyyy_0[i] =
            ta1_y_xx_xxyyyy_0[i] * fe_0 - ta1_y_xx_xxyyyy_1[i] * fe_0 + ta1_y_xxz_xxyyyy_0[i] * pa_z[i] - ta1_y_xxz_xxyyyy_1[i] * pc_z[i];

        ta1_y_xxzz_xxyyyz_0[i] = ta1_y_zz_xxyyyz_0[i] * fe_0 - ta1_y_zz_xxyyyz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_xzz_xyyyz_1[i] * fe_0 + ta1_y_xzz_xxyyyz_0[i] * pa_x[i] - ta1_y_xzz_xxyyyz_1[i] * pc_x[i];

        ta1_y_xxzz_xxyyzz_0[i] = ta1_y_zz_xxyyzz_0[i] * fe_0 - ta1_y_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xyyzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xzz_xyyzz_1[i] * fe_0 + ta1_y_xzz_xxyyzz_0[i] * pa_x[i] - ta1_y_xzz_xxyyzz_1[i] * pc_x[i];

        ta1_y_xxzz_xxyzzz_0[i] = ta1_y_zz_xxyzzz_0[i] * fe_0 - ta1_y_zz_xxyzzz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xzz_xyzzz_1[i] * fe_0 + ta1_y_xzz_xxyzzz_0[i] * pa_x[i] - ta1_y_xzz_xxyzzz_1[i] * pc_x[i];

        ta1_y_xxzz_xxzzzz_0[i] = ta1_y_zz_xxzzzz_0[i] * fe_0 - ta1_y_zz_xxzzzz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xzzzz_0[i] * fe_0 -
                                 2.0 * ta1_y_xzz_xzzzz_1[i] * fe_0 + ta1_y_xzz_xxzzzz_0[i] * pa_x[i] - ta1_y_xzz_xxzzzz_1[i] * pc_x[i];

        ta1_y_xxzz_xyyyyy_0[i] =
            ta1_y_xx_xyyyyy_0[i] * fe_0 - ta1_y_xx_xyyyyy_1[i] * fe_0 + ta1_y_xxz_xyyyyy_0[i] * pa_z[i] - ta1_y_xxz_xyyyyy_1[i] * pc_z[i];

        ta1_y_xxzz_xyyyyz_0[i] = ta1_y_zz_xyyyyz_0[i] * fe_0 - ta1_y_zz_xyyyyz_1[i] * fe_0 + ta1_y_xzz_yyyyz_0[i] * fe_0 -
                                 ta1_y_xzz_yyyyz_1[i] * fe_0 + ta1_y_xzz_xyyyyz_0[i] * pa_x[i] - ta1_y_xzz_xyyyyz_1[i] * pc_x[i];

        ta1_y_xxzz_xyyyzz_0[i] = ta1_y_zz_xyyyzz_0[i] * fe_0 - ta1_y_zz_xyyyzz_1[i] * fe_0 + ta1_y_xzz_yyyzz_0[i] * fe_0 -
                                 ta1_y_xzz_yyyzz_1[i] * fe_0 + ta1_y_xzz_xyyyzz_0[i] * pa_x[i] - ta1_y_xzz_xyyyzz_1[i] * pc_x[i];

        ta1_y_xxzz_xyyzzz_0[i] = ta1_y_zz_xyyzzz_0[i] * fe_0 - ta1_y_zz_xyyzzz_1[i] * fe_0 + ta1_y_xzz_yyzzz_0[i] * fe_0 -
                                 ta1_y_xzz_yyzzz_1[i] * fe_0 + ta1_y_xzz_xyyzzz_0[i] * pa_x[i] - ta1_y_xzz_xyyzzz_1[i] * pc_x[i];

        ta1_y_xxzz_xyzzzz_0[i] = ta1_y_zz_xyzzzz_0[i] * fe_0 - ta1_y_zz_xyzzzz_1[i] * fe_0 + ta1_y_xzz_yzzzz_0[i] * fe_0 -
                                 ta1_y_xzz_yzzzz_1[i] * fe_0 + ta1_y_xzz_xyzzzz_0[i] * pa_x[i] - ta1_y_xzz_xyzzzz_1[i] * pc_x[i];

        ta1_y_xxzz_xzzzzz_0[i] = ta1_y_zz_xzzzzz_0[i] * fe_0 - ta1_y_zz_xzzzzz_1[i] * fe_0 + ta1_y_xzz_zzzzz_0[i] * fe_0 -
                                 ta1_y_xzz_zzzzz_1[i] * fe_0 + ta1_y_xzz_xzzzzz_0[i] * pa_x[i] - ta1_y_xzz_xzzzzz_1[i] * pc_x[i];

        ta1_y_xxzz_yyyyyy_0[i] =
            ta1_y_zz_yyyyyy_0[i] * fe_0 - ta1_y_zz_yyyyyy_1[i] * fe_0 + ta1_y_xzz_yyyyyy_0[i] * pa_x[i] - ta1_y_xzz_yyyyyy_1[i] * pc_x[i];

        ta1_y_xxzz_yyyyyz_0[i] =
            ta1_y_zz_yyyyyz_0[i] * fe_0 - ta1_y_zz_yyyyyz_1[i] * fe_0 + ta1_y_xzz_yyyyyz_0[i] * pa_x[i] - ta1_y_xzz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxzz_yyyyzz_0[i] =
            ta1_y_zz_yyyyzz_0[i] * fe_0 - ta1_y_zz_yyyyzz_1[i] * fe_0 + ta1_y_xzz_yyyyzz_0[i] * pa_x[i] - ta1_y_xzz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxzz_yyyzzz_0[i] =
            ta1_y_zz_yyyzzz_0[i] * fe_0 - ta1_y_zz_yyyzzz_1[i] * fe_0 + ta1_y_xzz_yyyzzz_0[i] * pa_x[i] - ta1_y_xzz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxzz_yyzzzz_0[i] =
            ta1_y_zz_yyzzzz_0[i] * fe_0 - ta1_y_zz_yyzzzz_1[i] * fe_0 + ta1_y_xzz_yyzzzz_0[i] * pa_x[i] - ta1_y_xzz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxzz_yzzzzz_0[i] =
            ta1_y_zz_yzzzzz_0[i] * fe_0 - ta1_y_zz_yzzzzz_1[i] * fe_0 + ta1_y_xzz_yzzzzz_0[i] * pa_x[i] - ta1_y_xzz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxzz_zzzzzz_0[i] =
            ta1_y_zz_zzzzzz_0[i] * fe_0 - ta1_y_zz_zzzzzz_1[i] * fe_0 + ta1_y_xzz_zzzzzz_0[i] * pa_x[i] - ta1_y_xzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 588-616 components of targeted buffer : GI

    auto ta1_y_xyyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 588);

    auto ta1_y_xyyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 589);

    auto ta1_y_xyyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 590);

    auto ta1_y_xyyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 591);

    auto ta1_y_xyyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 592);

    auto ta1_y_xyyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 593);

    auto ta1_y_xyyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 594);

    auto ta1_y_xyyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 595);

    auto ta1_y_xyyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 596);

    auto ta1_y_xyyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 597);

    auto ta1_y_xyyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 598);

    auto ta1_y_xyyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 599);

    auto ta1_y_xyyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 600);

    auto ta1_y_xyyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 601);

    auto ta1_y_xyyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 602);

    auto ta1_y_xyyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 603);

    auto ta1_y_xyyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 604);

    auto ta1_y_xyyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 605);

    auto ta1_y_xyyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 606);

    auto ta1_y_xyyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 607);

    auto ta1_y_xyyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 608);

    auto ta1_y_xyyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 609);

    auto ta1_y_xyyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 610);

    auto ta1_y_xyyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 611);

    auto ta1_y_xyyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 612);

    auto ta1_y_xyyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 613);

    auto ta1_y_xyyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 614);

    auto ta1_y_xyyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 615);

#pragma omp simd aligned(pa_x,                    \
                             pc_x,                \
                             ta1_y_xyyy_xxxxxx_0, \
                             ta1_y_xyyy_xxxxxy_0, \
                             ta1_y_xyyy_xxxxxz_0, \
                             ta1_y_xyyy_xxxxyy_0, \
                             ta1_y_xyyy_xxxxyz_0, \
                             ta1_y_xyyy_xxxxzz_0, \
                             ta1_y_xyyy_xxxyyy_0, \
                             ta1_y_xyyy_xxxyyz_0, \
                             ta1_y_xyyy_xxxyzz_0, \
                             ta1_y_xyyy_xxxzzz_0, \
                             ta1_y_xyyy_xxyyyy_0, \
                             ta1_y_xyyy_xxyyyz_0, \
                             ta1_y_xyyy_xxyyzz_0, \
                             ta1_y_xyyy_xxyzzz_0, \
                             ta1_y_xyyy_xxzzzz_0, \
                             ta1_y_xyyy_xyyyyy_0, \
                             ta1_y_xyyy_xyyyyz_0, \
                             ta1_y_xyyy_xyyyzz_0, \
                             ta1_y_xyyy_xyyzzz_0, \
                             ta1_y_xyyy_xyzzzz_0, \
                             ta1_y_xyyy_xzzzzz_0, \
                             ta1_y_xyyy_yyyyyy_0, \
                             ta1_y_xyyy_yyyyyz_0, \
                             ta1_y_xyyy_yyyyzz_0, \
                             ta1_y_xyyy_yyyzzz_0, \
                             ta1_y_xyyy_yyzzzz_0, \
                             ta1_y_xyyy_yzzzzz_0, \
                             ta1_y_xyyy_zzzzzz_0, \
                             ta1_y_yyy_xxxxx_0,   \
                             ta1_y_yyy_xxxxx_1,   \
                             ta1_y_yyy_xxxxxx_0,  \
                             ta1_y_yyy_xxxxxx_1,  \
                             ta1_y_yyy_xxxxxy_0,  \
                             ta1_y_yyy_xxxxxy_1,  \
                             ta1_y_yyy_xxxxxz_0,  \
                             ta1_y_yyy_xxxxxz_1,  \
                             ta1_y_yyy_xxxxy_0,   \
                             ta1_y_yyy_xxxxy_1,   \
                             ta1_y_yyy_xxxxyy_0,  \
                             ta1_y_yyy_xxxxyy_1,  \
                             ta1_y_yyy_xxxxyz_0,  \
                             ta1_y_yyy_xxxxyz_1,  \
                             ta1_y_yyy_xxxxz_0,   \
                             ta1_y_yyy_xxxxz_1,   \
                             ta1_y_yyy_xxxxzz_0,  \
                             ta1_y_yyy_xxxxzz_1,  \
                             ta1_y_yyy_xxxyy_0,   \
                             ta1_y_yyy_xxxyy_1,   \
                             ta1_y_yyy_xxxyyy_0,  \
                             ta1_y_yyy_xxxyyy_1,  \
                             ta1_y_yyy_xxxyyz_0,  \
                             ta1_y_yyy_xxxyyz_1,  \
                             ta1_y_yyy_xxxyz_0,   \
                             ta1_y_yyy_xxxyz_1,   \
                             ta1_y_yyy_xxxyzz_0,  \
                             ta1_y_yyy_xxxyzz_1,  \
                             ta1_y_yyy_xxxzz_0,   \
                             ta1_y_yyy_xxxzz_1,   \
                             ta1_y_yyy_xxxzzz_0,  \
                             ta1_y_yyy_xxxzzz_1,  \
                             ta1_y_yyy_xxyyy_0,   \
                             ta1_y_yyy_xxyyy_1,   \
                             ta1_y_yyy_xxyyyy_0,  \
                             ta1_y_yyy_xxyyyy_1,  \
                             ta1_y_yyy_xxyyyz_0,  \
                             ta1_y_yyy_xxyyyz_1,  \
                             ta1_y_yyy_xxyyz_0,   \
                             ta1_y_yyy_xxyyz_1,   \
                             ta1_y_yyy_xxyyzz_0,  \
                             ta1_y_yyy_xxyyzz_1,  \
                             ta1_y_yyy_xxyzz_0,   \
                             ta1_y_yyy_xxyzz_1,   \
                             ta1_y_yyy_xxyzzz_0,  \
                             ta1_y_yyy_xxyzzz_1,  \
                             ta1_y_yyy_xxzzz_0,   \
                             ta1_y_yyy_xxzzz_1,   \
                             ta1_y_yyy_xxzzzz_0,  \
                             ta1_y_yyy_xxzzzz_1,  \
                             ta1_y_yyy_xyyyy_0,   \
                             ta1_y_yyy_xyyyy_1,   \
                             ta1_y_yyy_xyyyyy_0,  \
                             ta1_y_yyy_xyyyyy_1,  \
                             ta1_y_yyy_xyyyyz_0,  \
                             ta1_y_yyy_xyyyyz_1,  \
                             ta1_y_yyy_xyyyz_0,   \
                             ta1_y_yyy_xyyyz_1,   \
                             ta1_y_yyy_xyyyzz_0,  \
                             ta1_y_yyy_xyyyzz_1,  \
                             ta1_y_yyy_xyyzz_0,   \
                             ta1_y_yyy_xyyzz_1,   \
                             ta1_y_yyy_xyyzzz_0,  \
                             ta1_y_yyy_xyyzzz_1,  \
                             ta1_y_yyy_xyzzz_0,   \
                             ta1_y_yyy_xyzzz_1,   \
                             ta1_y_yyy_xyzzzz_0,  \
                             ta1_y_yyy_xyzzzz_1,  \
                             ta1_y_yyy_xzzzz_0,   \
                             ta1_y_yyy_xzzzz_1,   \
                             ta1_y_yyy_xzzzzz_0,  \
                             ta1_y_yyy_xzzzzz_1,  \
                             ta1_y_yyy_yyyyy_0,   \
                             ta1_y_yyy_yyyyy_1,   \
                             ta1_y_yyy_yyyyyy_0,  \
                             ta1_y_yyy_yyyyyy_1,  \
                             ta1_y_yyy_yyyyyz_0,  \
                             ta1_y_yyy_yyyyyz_1,  \
                             ta1_y_yyy_yyyyz_0,   \
                             ta1_y_yyy_yyyyz_1,   \
                             ta1_y_yyy_yyyyzz_0,  \
                             ta1_y_yyy_yyyyzz_1,  \
                             ta1_y_yyy_yyyzz_0,   \
                             ta1_y_yyy_yyyzz_1,   \
                             ta1_y_yyy_yyyzzz_0,  \
                             ta1_y_yyy_yyyzzz_1,  \
                             ta1_y_yyy_yyzzz_0,   \
                             ta1_y_yyy_yyzzz_1,   \
                             ta1_y_yyy_yyzzzz_0,  \
                             ta1_y_yyy_yyzzzz_1,  \
                             ta1_y_yyy_yzzzz_0,   \
                             ta1_y_yyy_yzzzz_1,   \
                             ta1_y_yyy_yzzzzz_0,  \
                             ta1_y_yyy_yzzzzz_1,  \
                             ta1_y_yyy_zzzzz_0,   \
                             ta1_y_yyy_zzzzz_1,   \
                             ta1_y_yyy_zzzzzz_0,  \
                             ta1_y_yyy_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyy_xxxxxx_0[i] =
            6.0 * ta1_y_yyy_xxxxx_0[i] * fe_0 - 6.0 * ta1_y_yyy_xxxxx_1[i] * fe_0 + ta1_y_yyy_xxxxxx_0[i] * pa_x[i] - ta1_y_yyy_xxxxxx_1[i] * pc_x[i];

        ta1_y_xyyy_xxxxxy_0[i] =
            5.0 * ta1_y_yyy_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_yyy_xxxxy_1[i] * fe_0 + ta1_y_yyy_xxxxxy_0[i] * pa_x[i] - ta1_y_yyy_xxxxxy_1[i] * pc_x[i];

        ta1_y_xyyy_xxxxxz_0[i] =
            5.0 * ta1_y_yyy_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_yyy_xxxxz_1[i] * fe_0 + ta1_y_yyy_xxxxxz_0[i] * pa_x[i] - ta1_y_yyy_xxxxxz_1[i] * pc_x[i];

        ta1_y_xyyy_xxxxyy_0[i] =
            4.0 * ta1_y_yyy_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxyy_1[i] * fe_0 + ta1_y_yyy_xxxxyy_0[i] * pa_x[i] - ta1_y_yyy_xxxxyy_1[i] * pc_x[i];

        ta1_y_xyyy_xxxxyz_0[i] =
            4.0 * ta1_y_yyy_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxyz_1[i] * fe_0 + ta1_y_yyy_xxxxyz_0[i] * pa_x[i] - ta1_y_yyy_xxxxyz_1[i] * pc_x[i];

        ta1_y_xyyy_xxxxzz_0[i] =
            4.0 * ta1_y_yyy_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxzz_1[i] * fe_0 + ta1_y_yyy_xxxxzz_0[i] * pa_x[i] - ta1_y_yyy_xxxxzz_1[i] * pc_x[i];

        ta1_y_xyyy_xxxyyy_0[i] =
            3.0 * ta1_y_yyy_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxyyy_1[i] * fe_0 + ta1_y_yyy_xxxyyy_0[i] * pa_x[i] - ta1_y_yyy_xxxyyy_1[i] * pc_x[i];

        ta1_y_xyyy_xxxyyz_0[i] =
            3.0 * ta1_y_yyy_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxyyz_1[i] * fe_0 + ta1_y_yyy_xxxyyz_0[i] * pa_x[i] - ta1_y_yyy_xxxyyz_1[i] * pc_x[i];

        ta1_y_xyyy_xxxyzz_0[i] =
            3.0 * ta1_y_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxyzz_1[i] * fe_0 + ta1_y_yyy_xxxyzz_0[i] * pa_x[i] - ta1_y_yyy_xxxyzz_1[i] * pc_x[i];

        ta1_y_xyyy_xxxzzz_0[i] =
            3.0 * ta1_y_yyy_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxzzz_1[i] * fe_0 + ta1_y_yyy_xxxzzz_0[i] * pa_x[i] - ta1_y_yyy_xxxzzz_1[i] * pc_x[i];

        ta1_y_xyyy_xxyyyy_0[i] =
            2.0 * ta1_y_yyy_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyyyy_1[i] * fe_0 + ta1_y_yyy_xxyyyy_0[i] * pa_x[i] - ta1_y_yyy_xxyyyy_1[i] * pc_x[i];

        ta1_y_xyyy_xxyyyz_0[i] =
            2.0 * ta1_y_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyyyz_1[i] * fe_0 + ta1_y_yyy_xxyyyz_0[i] * pa_x[i] - ta1_y_yyy_xxyyyz_1[i] * pc_x[i];

        ta1_y_xyyy_xxyyzz_0[i] =
            2.0 * ta1_y_yyy_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyyzz_1[i] * fe_0 + ta1_y_yyy_xxyyzz_0[i] * pa_x[i] - ta1_y_yyy_xxyyzz_1[i] * pc_x[i];

        ta1_y_xyyy_xxyzzz_0[i] =
            2.0 * ta1_y_yyy_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyzzz_1[i] * fe_0 + ta1_y_yyy_xxyzzz_0[i] * pa_x[i] - ta1_y_yyy_xxyzzz_1[i] * pc_x[i];

        ta1_y_xyyy_xxzzzz_0[i] =
            2.0 * ta1_y_yyy_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xzzzz_1[i] * fe_0 + ta1_y_yyy_xxzzzz_0[i] * pa_x[i] - ta1_y_yyy_xxzzzz_1[i] * pc_x[i];

        ta1_y_xyyy_xyyyyy_0[i] =
            ta1_y_yyy_yyyyy_0[i] * fe_0 - ta1_y_yyy_yyyyy_1[i] * fe_0 + ta1_y_yyy_xyyyyy_0[i] * pa_x[i] - ta1_y_yyy_xyyyyy_1[i] * pc_x[i];

        ta1_y_xyyy_xyyyyz_0[i] =
            ta1_y_yyy_yyyyz_0[i] * fe_0 - ta1_y_yyy_yyyyz_1[i] * fe_0 + ta1_y_yyy_xyyyyz_0[i] * pa_x[i] - ta1_y_yyy_xyyyyz_1[i] * pc_x[i];

        ta1_y_xyyy_xyyyzz_0[i] =
            ta1_y_yyy_yyyzz_0[i] * fe_0 - ta1_y_yyy_yyyzz_1[i] * fe_0 + ta1_y_yyy_xyyyzz_0[i] * pa_x[i] - ta1_y_yyy_xyyyzz_1[i] * pc_x[i];

        ta1_y_xyyy_xyyzzz_0[i] =
            ta1_y_yyy_yyzzz_0[i] * fe_0 - ta1_y_yyy_yyzzz_1[i] * fe_0 + ta1_y_yyy_xyyzzz_0[i] * pa_x[i] - ta1_y_yyy_xyyzzz_1[i] * pc_x[i];

        ta1_y_xyyy_xyzzzz_0[i] =
            ta1_y_yyy_yzzzz_0[i] * fe_0 - ta1_y_yyy_yzzzz_1[i] * fe_0 + ta1_y_yyy_xyzzzz_0[i] * pa_x[i] - ta1_y_yyy_xyzzzz_1[i] * pc_x[i];

        ta1_y_xyyy_xzzzzz_0[i] =
            ta1_y_yyy_zzzzz_0[i] * fe_0 - ta1_y_yyy_zzzzz_1[i] * fe_0 + ta1_y_yyy_xzzzzz_0[i] * pa_x[i] - ta1_y_yyy_xzzzzz_1[i] * pc_x[i];

        ta1_y_xyyy_yyyyyy_0[i] = ta1_y_yyy_yyyyyy_0[i] * pa_x[i] - ta1_y_yyy_yyyyyy_1[i] * pc_x[i];

        ta1_y_xyyy_yyyyyz_0[i] = ta1_y_yyy_yyyyyz_0[i] * pa_x[i] - ta1_y_yyy_yyyyyz_1[i] * pc_x[i];

        ta1_y_xyyy_yyyyzz_0[i] = ta1_y_yyy_yyyyzz_0[i] * pa_x[i] - ta1_y_yyy_yyyyzz_1[i] * pc_x[i];

        ta1_y_xyyy_yyyzzz_0[i] = ta1_y_yyy_yyyzzz_0[i] * pa_x[i] - ta1_y_yyy_yyyzzz_1[i] * pc_x[i];

        ta1_y_xyyy_yyzzzz_0[i] = ta1_y_yyy_yyzzzz_0[i] * pa_x[i] - ta1_y_yyy_yyzzzz_1[i] * pc_x[i];

        ta1_y_xyyy_yzzzzz_0[i] = ta1_y_yyy_yzzzzz_0[i] * pa_x[i] - ta1_y_yyy_yzzzzz_1[i] * pc_x[i];

        ta1_y_xyyy_zzzzzz_0[i] = ta1_y_yyy_zzzzzz_0[i] * pa_x[i] - ta1_y_yyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 616-644 components of targeted buffer : GI

    auto ta1_y_xyyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 616);

    auto ta1_y_xyyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 617);

    auto ta1_y_xyyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 618);

    auto ta1_y_xyyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 619);

    auto ta1_y_xyyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 620);

    auto ta1_y_xyyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 621);

    auto ta1_y_xyyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 622);

    auto ta1_y_xyyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 623);

    auto ta1_y_xyyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 624);

    auto ta1_y_xyyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 625);

    auto ta1_y_xyyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 626);

    auto ta1_y_xyyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 627);

    auto ta1_y_xyyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 628);

    auto ta1_y_xyyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 629);

    auto ta1_y_xyyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 630);

    auto ta1_y_xyyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 631);

    auto ta1_y_xyyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 632);

    auto ta1_y_xyyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 633);

    auto ta1_y_xyyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 634);

    auto ta1_y_xyyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 635);

    auto ta1_y_xyyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 636);

    auto ta1_y_xyyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 637);

    auto ta1_y_xyyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 638);

    auto ta1_y_xyyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 639);

    auto ta1_y_xyyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 640);

    auto ta1_y_xyyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 641);

    auto ta1_y_xyyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 642);

    auto ta1_y_xyyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 643);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_y_xyy_xxxxxx_0,  \
                             ta1_y_xyy_xxxxxx_1,  \
                             ta1_y_xyy_xxxxxy_0,  \
                             ta1_y_xyy_xxxxxy_1,  \
                             ta1_y_xyy_xxxxyy_0,  \
                             ta1_y_xyy_xxxxyy_1,  \
                             ta1_y_xyy_xxxyyy_0,  \
                             ta1_y_xyy_xxxyyy_1,  \
                             ta1_y_xyy_xxyyyy_0,  \
                             ta1_y_xyy_xxyyyy_1,  \
                             ta1_y_xyy_xyyyyy_0,  \
                             ta1_y_xyy_xyyyyy_1,  \
                             ta1_y_xyyz_xxxxxx_0, \
                             ta1_y_xyyz_xxxxxy_0, \
                             ta1_y_xyyz_xxxxxz_0, \
                             ta1_y_xyyz_xxxxyy_0, \
                             ta1_y_xyyz_xxxxyz_0, \
                             ta1_y_xyyz_xxxxzz_0, \
                             ta1_y_xyyz_xxxyyy_0, \
                             ta1_y_xyyz_xxxyyz_0, \
                             ta1_y_xyyz_xxxyzz_0, \
                             ta1_y_xyyz_xxxzzz_0, \
                             ta1_y_xyyz_xxyyyy_0, \
                             ta1_y_xyyz_xxyyyz_0, \
                             ta1_y_xyyz_xxyyzz_0, \
                             ta1_y_xyyz_xxyzzz_0, \
                             ta1_y_xyyz_xxzzzz_0, \
                             ta1_y_xyyz_xyyyyy_0, \
                             ta1_y_xyyz_xyyyyz_0, \
                             ta1_y_xyyz_xyyyzz_0, \
                             ta1_y_xyyz_xyyzzz_0, \
                             ta1_y_xyyz_xyzzzz_0, \
                             ta1_y_xyyz_xzzzzz_0, \
                             ta1_y_xyyz_yyyyyy_0, \
                             ta1_y_xyyz_yyyyyz_0, \
                             ta1_y_xyyz_yyyyzz_0, \
                             ta1_y_xyyz_yyyzzz_0, \
                             ta1_y_xyyz_yyzzzz_0, \
                             ta1_y_xyyz_yzzzzz_0, \
                             ta1_y_xyyz_zzzzzz_0, \
                             ta1_y_yyz_xxxxxz_0,  \
                             ta1_y_yyz_xxxxxz_1,  \
                             ta1_y_yyz_xxxxyz_0,  \
                             ta1_y_yyz_xxxxyz_1,  \
                             ta1_y_yyz_xxxxz_0,   \
                             ta1_y_yyz_xxxxz_1,   \
                             ta1_y_yyz_xxxxzz_0,  \
                             ta1_y_yyz_xxxxzz_1,  \
                             ta1_y_yyz_xxxyyz_0,  \
                             ta1_y_yyz_xxxyyz_1,  \
                             ta1_y_yyz_xxxyz_0,   \
                             ta1_y_yyz_xxxyz_1,   \
                             ta1_y_yyz_xxxyzz_0,  \
                             ta1_y_yyz_xxxyzz_1,  \
                             ta1_y_yyz_xxxzz_0,   \
                             ta1_y_yyz_xxxzz_1,   \
                             ta1_y_yyz_xxxzzz_0,  \
                             ta1_y_yyz_xxxzzz_1,  \
                             ta1_y_yyz_xxyyyz_0,  \
                             ta1_y_yyz_xxyyyz_1,  \
                             ta1_y_yyz_xxyyz_0,   \
                             ta1_y_yyz_xxyyz_1,   \
                             ta1_y_yyz_xxyyzz_0,  \
                             ta1_y_yyz_xxyyzz_1,  \
                             ta1_y_yyz_xxyzz_0,   \
                             ta1_y_yyz_xxyzz_1,   \
                             ta1_y_yyz_xxyzzz_0,  \
                             ta1_y_yyz_xxyzzz_1,  \
                             ta1_y_yyz_xxzzz_0,   \
                             ta1_y_yyz_xxzzz_1,   \
                             ta1_y_yyz_xxzzzz_0,  \
                             ta1_y_yyz_xxzzzz_1,  \
                             ta1_y_yyz_xyyyyz_0,  \
                             ta1_y_yyz_xyyyyz_1,  \
                             ta1_y_yyz_xyyyz_0,   \
                             ta1_y_yyz_xyyyz_1,   \
                             ta1_y_yyz_xyyyzz_0,  \
                             ta1_y_yyz_xyyyzz_1,  \
                             ta1_y_yyz_xyyzz_0,   \
                             ta1_y_yyz_xyyzz_1,   \
                             ta1_y_yyz_xyyzzz_0,  \
                             ta1_y_yyz_xyyzzz_1,  \
                             ta1_y_yyz_xyzzz_0,   \
                             ta1_y_yyz_xyzzz_1,   \
                             ta1_y_yyz_xyzzzz_0,  \
                             ta1_y_yyz_xyzzzz_1,  \
                             ta1_y_yyz_xzzzz_0,   \
                             ta1_y_yyz_xzzzz_1,   \
                             ta1_y_yyz_xzzzzz_0,  \
                             ta1_y_yyz_xzzzzz_1,  \
                             ta1_y_yyz_yyyyyy_0,  \
                             ta1_y_yyz_yyyyyy_1,  \
                             ta1_y_yyz_yyyyyz_0,  \
                             ta1_y_yyz_yyyyyz_1,  \
                             ta1_y_yyz_yyyyz_0,   \
                             ta1_y_yyz_yyyyz_1,   \
                             ta1_y_yyz_yyyyzz_0,  \
                             ta1_y_yyz_yyyyzz_1,  \
                             ta1_y_yyz_yyyzz_0,   \
                             ta1_y_yyz_yyyzz_1,   \
                             ta1_y_yyz_yyyzzz_0,  \
                             ta1_y_yyz_yyyzzz_1,  \
                             ta1_y_yyz_yyzzz_0,   \
                             ta1_y_yyz_yyzzz_1,   \
                             ta1_y_yyz_yyzzzz_0,  \
                             ta1_y_yyz_yyzzzz_1,  \
                             ta1_y_yyz_yzzzz_0,   \
                             ta1_y_yyz_yzzzz_1,   \
                             ta1_y_yyz_yzzzzz_0,  \
                             ta1_y_yyz_yzzzzz_1,  \
                             ta1_y_yyz_zzzzz_0,   \
                             ta1_y_yyz_zzzzz_1,   \
                             ta1_y_yyz_zzzzzz_0,  \
                             ta1_y_yyz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyz_xxxxxx_0[i] = ta1_y_xyy_xxxxxx_0[i] * pa_z[i] - ta1_y_xyy_xxxxxx_1[i] * pc_z[i];

        ta1_y_xyyz_xxxxxy_0[i] = ta1_y_xyy_xxxxxy_0[i] * pa_z[i] - ta1_y_xyy_xxxxxy_1[i] * pc_z[i];

        ta1_y_xyyz_xxxxxz_0[i] =
            5.0 * ta1_y_yyz_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_yyz_xxxxz_1[i] * fe_0 + ta1_y_yyz_xxxxxz_0[i] * pa_x[i] - ta1_y_yyz_xxxxxz_1[i] * pc_x[i];

        ta1_y_xyyz_xxxxyy_0[i] = ta1_y_xyy_xxxxyy_0[i] * pa_z[i] - ta1_y_xyy_xxxxyy_1[i] * pc_z[i];

        ta1_y_xyyz_xxxxyz_0[i] =
            4.0 * ta1_y_yyz_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_yyz_xxxyz_1[i] * fe_0 + ta1_y_yyz_xxxxyz_0[i] * pa_x[i] - ta1_y_yyz_xxxxyz_1[i] * pc_x[i];

        ta1_y_xyyz_xxxxzz_0[i] =
            4.0 * ta1_y_yyz_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_yyz_xxxzz_1[i] * fe_0 + ta1_y_yyz_xxxxzz_0[i] * pa_x[i] - ta1_y_yyz_xxxxzz_1[i] * pc_x[i];

        ta1_y_xyyz_xxxyyy_0[i] = ta1_y_xyy_xxxyyy_0[i] * pa_z[i] - ta1_y_xyy_xxxyyy_1[i] * pc_z[i];

        ta1_y_xyyz_xxxyyz_0[i] =
            3.0 * ta1_y_yyz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_yyz_xxyyz_1[i] * fe_0 + ta1_y_yyz_xxxyyz_0[i] * pa_x[i] - ta1_y_yyz_xxxyyz_1[i] * pc_x[i];

        ta1_y_xyyz_xxxyzz_0[i] =
            3.0 * ta1_y_yyz_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yyz_xxyzz_1[i] * fe_0 + ta1_y_yyz_xxxyzz_0[i] * pa_x[i] - ta1_y_yyz_xxxyzz_1[i] * pc_x[i];

        ta1_y_xyyz_xxxzzz_0[i] =
            3.0 * ta1_y_yyz_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_yyz_xxzzz_1[i] * fe_0 + ta1_y_yyz_xxxzzz_0[i] * pa_x[i] - ta1_y_yyz_xxxzzz_1[i] * pc_x[i];

        ta1_y_xyyz_xxyyyy_0[i] = ta1_y_xyy_xxyyyy_0[i] * pa_z[i] - ta1_y_xyy_xxyyyy_1[i] * pc_z[i];

        ta1_y_xyyz_xxyyyz_0[i] =
            2.0 * ta1_y_yyz_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyyyz_1[i] * fe_0 + ta1_y_yyz_xxyyyz_0[i] * pa_x[i] - ta1_y_yyz_xxyyyz_1[i] * pc_x[i];

        ta1_y_xyyz_xxyyzz_0[i] =
            2.0 * ta1_y_yyz_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyyzz_1[i] * fe_0 + ta1_y_yyz_xxyyzz_0[i] * pa_x[i] - ta1_y_yyz_xxyyzz_1[i] * pc_x[i];

        ta1_y_xyyz_xxyzzz_0[i] =
            2.0 * ta1_y_yyz_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyzzz_1[i] * fe_0 + ta1_y_yyz_xxyzzz_0[i] * pa_x[i] - ta1_y_yyz_xxyzzz_1[i] * pc_x[i];

        ta1_y_xyyz_xxzzzz_0[i] =
            2.0 * ta1_y_yyz_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xzzzz_1[i] * fe_0 + ta1_y_yyz_xxzzzz_0[i] * pa_x[i] - ta1_y_yyz_xxzzzz_1[i] * pc_x[i];

        ta1_y_xyyz_xyyyyy_0[i] = ta1_y_xyy_xyyyyy_0[i] * pa_z[i] - ta1_y_xyy_xyyyyy_1[i] * pc_z[i];

        ta1_y_xyyz_xyyyyz_0[i] =
            ta1_y_yyz_yyyyz_0[i] * fe_0 - ta1_y_yyz_yyyyz_1[i] * fe_0 + ta1_y_yyz_xyyyyz_0[i] * pa_x[i] - ta1_y_yyz_xyyyyz_1[i] * pc_x[i];

        ta1_y_xyyz_xyyyzz_0[i] =
            ta1_y_yyz_yyyzz_0[i] * fe_0 - ta1_y_yyz_yyyzz_1[i] * fe_0 + ta1_y_yyz_xyyyzz_0[i] * pa_x[i] - ta1_y_yyz_xyyyzz_1[i] * pc_x[i];

        ta1_y_xyyz_xyyzzz_0[i] =
            ta1_y_yyz_yyzzz_0[i] * fe_0 - ta1_y_yyz_yyzzz_1[i] * fe_0 + ta1_y_yyz_xyyzzz_0[i] * pa_x[i] - ta1_y_yyz_xyyzzz_1[i] * pc_x[i];

        ta1_y_xyyz_xyzzzz_0[i] =
            ta1_y_yyz_yzzzz_0[i] * fe_0 - ta1_y_yyz_yzzzz_1[i] * fe_0 + ta1_y_yyz_xyzzzz_0[i] * pa_x[i] - ta1_y_yyz_xyzzzz_1[i] * pc_x[i];

        ta1_y_xyyz_xzzzzz_0[i] =
            ta1_y_yyz_zzzzz_0[i] * fe_0 - ta1_y_yyz_zzzzz_1[i] * fe_0 + ta1_y_yyz_xzzzzz_0[i] * pa_x[i] - ta1_y_yyz_xzzzzz_1[i] * pc_x[i];

        ta1_y_xyyz_yyyyyy_0[i] = ta1_y_yyz_yyyyyy_0[i] * pa_x[i] - ta1_y_yyz_yyyyyy_1[i] * pc_x[i];

        ta1_y_xyyz_yyyyyz_0[i] = ta1_y_yyz_yyyyyz_0[i] * pa_x[i] - ta1_y_yyz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xyyz_yyyyzz_0[i] = ta1_y_yyz_yyyyzz_0[i] * pa_x[i] - ta1_y_yyz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xyyz_yyyzzz_0[i] = ta1_y_yyz_yyyzzz_0[i] * pa_x[i] - ta1_y_yyz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xyyz_yyzzzz_0[i] = ta1_y_yyz_yyzzzz_0[i] * pa_x[i] - ta1_y_yyz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xyyz_yzzzzz_0[i] = ta1_y_yyz_yzzzzz_0[i] * pa_x[i] - ta1_y_yyz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xyyz_zzzzzz_0[i] = ta1_y_yyz_zzzzzz_0[i] * pa_x[i] - ta1_y_yyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 644-672 components of targeted buffer : GI

    auto ta1_y_xyzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 644);

    auto ta1_y_xyzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 645);

    auto ta1_y_xyzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 646);

    auto ta1_y_xyzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 647);

    auto ta1_y_xyzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 648);

    auto ta1_y_xyzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 649);

    auto ta1_y_xyzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 650);

    auto ta1_y_xyzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 651);

    auto ta1_y_xyzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 652);

    auto ta1_y_xyzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 653);

    auto ta1_y_xyzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 654);

    auto ta1_y_xyzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 655);

    auto ta1_y_xyzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 656);

    auto ta1_y_xyzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 657);

    auto ta1_y_xyzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 658);

    auto ta1_y_xyzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 659);

    auto ta1_y_xyzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 660);

    auto ta1_y_xyzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 661);

    auto ta1_y_xyzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 662);

    auto ta1_y_xyzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 663);

    auto ta1_y_xyzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 664);

    auto ta1_y_xyzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 665);

    auto ta1_y_xyzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 666);

    auto ta1_y_xyzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 667);

    auto ta1_y_xyzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 668);

    auto ta1_y_xyzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 669);

    auto ta1_y_xyzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 670);

    auto ta1_y_xyzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 671);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_y_xyzz_xxxxxx_0, \
                             ta1_y_xyzz_xxxxxy_0, \
                             ta1_y_xyzz_xxxxxz_0, \
                             ta1_y_xyzz_xxxxyy_0, \
                             ta1_y_xyzz_xxxxyz_0, \
                             ta1_y_xyzz_xxxxzz_0, \
                             ta1_y_xyzz_xxxyyy_0, \
                             ta1_y_xyzz_xxxyyz_0, \
                             ta1_y_xyzz_xxxyzz_0, \
                             ta1_y_xyzz_xxxzzz_0, \
                             ta1_y_xyzz_xxyyyy_0, \
                             ta1_y_xyzz_xxyyyz_0, \
                             ta1_y_xyzz_xxyyzz_0, \
                             ta1_y_xyzz_xxyzzz_0, \
                             ta1_y_xyzz_xxzzzz_0, \
                             ta1_y_xyzz_xyyyyy_0, \
                             ta1_y_xyzz_xyyyyz_0, \
                             ta1_y_xyzz_xyyyzz_0, \
                             ta1_y_xyzz_xyyzzz_0, \
                             ta1_y_xyzz_xyzzzz_0, \
                             ta1_y_xyzz_xzzzzz_0, \
                             ta1_y_xyzz_yyyyyy_0, \
                             ta1_y_xyzz_yyyyyz_0, \
                             ta1_y_xyzz_yyyyzz_0, \
                             ta1_y_xyzz_yyyzzz_0, \
                             ta1_y_xyzz_yyzzzz_0, \
                             ta1_y_xyzz_yzzzzz_0, \
                             ta1_y_xyzz_zzzzzz_0, \
                             ta1_y_xzz_xxxxxx_0,  \
                             ta1_y_xzz_xxxxxx_1,  \
                             ta1_y_xzz_xxxxxz_0,  \
                             ta1_y_xzz_xxxxxz_1,  \
                             ta1_y_xzz_xxxxzz_0,  \
                             ta1_y_xzz_xxxxzz_1,  \
                             ta1_y_xzz_xxxzzz_0,  \
                             ta1_y_xzz_xxxzzz_1,  \
                             ta1_y_xzz_xxzzzz_0,  \
                             ta1_y_xzz_xxzzzz_1,  \
                             ta1_y_xzz_xzzzzz_0,  \
                             ta1_y_xzz_xzzzzz_1,  \
                             ta1_y_yzz_xxxxxy_0,  \
                             ta1_y_yzz_xxxxxy_1,  \
                             ta1_y_yzz_xxxxy_0,   \
                             ta1_y_yzz_xxxxy_1,   \
                             ta1_y_yzz_xxxxyy_0,  \
                             ta1_y_yzz_xxxxyy_1,  \
                             ta1_y_yzz_xxxxyz_0,  \
                             ta1_y_yzz_xxxxyz_1,  \
                             ta1_y_yzz_xxxyy_0,   \
                             ta1_y_yzz_xxxyy_1,   \
                             ta1_y_yzz_xxxyyy_0,  \
                             ta1_y_yzz_xxxyyy_1,  \
                             ta1_y_yzz_xxxyyz_0,  \
                             ta1_y_yzz_xxxyyz_1,  \
                             ta1_y_yzz_xxxyz_0,   \
                             ta1_y_yzz_xxxyz_1,   \
                             ta1_y_yzz_xxxyzz_0,  \
                             ta1_y_yzz_xxxyzz_1,  \
                             ta1_y_yzz_xxyyy_0,   \
                             ta1_y_yzz_xxyyy_1,   \
                             ta1_y_yzz_xxyyyy_0,  \
                             ta1_y_yzz_xxyyyy_1,  \
                             ta1_y_yzz_xxyyyz_0,  \
                             ta1_y_yzz_xxyyyz_1,  \
                             ta1_y_yzz_xxyyz_0,   \
                             ta1_y_yzz_xxyyz_1,   \
                             ta1_y_yzz_xxyyzz_0,  \
                             ta1_y_yzz_xxyyzz_1,  \
                             ta1_y_yzz_xxyzz_0,   \
                             ta1_y_yzz_xxyzz_1,   \
                             ta1_y_yzz_xxyzzz_0,  \
                             ta1_y_yzz_xxyzzz_1,  \
                             ta1_y_yzz_xyyyy_0,   \
                             ta1_y_yzz_xyyyy_1,   \
                             ta1_y_yzz_xyyyyy_0,  \
                             ta1_y_yzz_xyyyyy_1,  \
                             ta1_y_yzz_xyyyyz_0,  \
                             ta1_y_yzz_xyyyyz_1,  \
                             ta1_y_yzz_xyyyz_0,   \
                             ta1_y_yzz_xyyyz_1,   \
                             ta1_y_yzz_xyyyzz_0,  \
                             ta1_y_yzz_xyyyzz_1,  \
                             ta1_y_yzz_xyyzz_0,   \
                             ta1_y_yzz_xyyzz_1,   \
                             ta1_y_yzz_xyyzzz_0,  \
                             ta1_y_yzz_xyyzzz_1,  \
                             ta1_y_yzz_xyzzz_0,   \
                             ta1_y_yzz_xyzzz_1,   \
                             ta1_y_yzz_xyzzzz_0,  \
                             ta1_y_yzz_xyzzzz_1,  \
                             ta1_y_yzz_yyyyy_0,   \
                             ta1_y_yzz_yyyyy_1,   \
                             ta1_y_yzz_yyyyyy_0,  \
                             ta1_y_yzz_yyyyyy_1,  \
                             ta1_y_yzz_yyyyyz_0,  \
                             ta1_y_yzz_yyyyyz_1,  \
                             ta1_y_yzz_yyyyz_0,   \
                             ta1_y_yzz_yyyyz_1,   \
                             ta1_y_yzz_yyyyzz_0,  \
                             ta1_y_yzz_yyyyzz_1,  \
                             ta1_y_yzz_yyyzz_0,   \
                             ta1_y_yzz_yyyzz_1,   \
                             ta1_y_yzz_yyyzzz_0,  \
                             ta1_y_yzz_yyyzzz_1,  \
                             ta1_y_yzz_yyzzz_0,   \
                             ta1_y_yzz_yyzzz_1,   \
                             ta1_y_yzz_yyzzzz_0,  \
                             ta1_y_yzz_yyzzzz_1,  \
                             ta1_y_yzz_yzzzz_0,   \
                             ta1_y_yzz_yzzzz_1,   \
                             ta1_y_yzz_yzzzzz_0,  \
                             ta1_y_yzz_yzzzzz_1,  \
                             ta1_y_yzz_zzzzzz_0,  \
                             ta1_y_yzz_zzzzzz_1,  \
                             ta_xzz_xxxxxx_1,     \
                             ta_xzz_xxxxxz_1,     \
                             ta_xzz_xxxxzz_1,     \
                             ta_xzz_xxxzzz_1,     \
                             ta_xzz_xxzzzz_1,     \
                             ta_xzz_xzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzz_xxxxxx_0[i] = ta_xzz_xxxxxx_1[i] + ta1_y_xzz_xxxxxx_0[i] * pa_y[i] - ta1_y_xzz_xxxxxx_1[i] * pc_y[i];

        ta1_y_xyzz_xxxxxy_0[i] =
            5.0 * ta1_y_yzz_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_yzz_xxxxy_1[i] * fe_0 + ta1_y_yzz_xxxxxy_0[i] * pa_x[i] - ta1_y_yzz_xxxxxy_1[i] * pc_x[i];

        ta1_y_xyzz_xxxxxz_0[i] = ta_xzz_xxxxxz_1[i] + ta1_y_xzz_xxxxxz_0[i] * pa_y[i] - ta1_y_xzz_xxxxxz_1[i] * pc_y[i];

        ta1_y_xyzz_xxxxyy_0[i] =
            4.0 * ta1_y_yzz_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_yzz_xxxyy_1[i] * fe_0 + ta1_y_yzz_xxxxyy_0[i] * pa_x[i] - ta1_y_yzz_xxxxyy_1[i] * pc_x[i];

        ta1_y_xyzz_xxxxyz_0[i] =
            4.0 * ta1_y_yzz_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_yzz_xxxyz_1[i] * fe_0 + ta1_y_yzz_xxxxyz_0[i] * pa_x[i] - ta1_y_yzz_xxxxyz_1[i] * pc_x[i];

        ta1_y_xyzz_xxxxzz_0[i] = ta_xzz_xxxxzz_1[i] + ta1_y_xzz_xxxxzz_0[i] * pa_y[i] - ta1_y_xzz_xxxxzz_1[i] * pc_y[i];

        ta1_y_xyzz_xxxyyy_0[i] =
            3.0 * ta1_y_yzz_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxyyy_1[i] * fe_0 + ta1_y_yzz_xxxyyy_0[i] * pa_x[i] - ta1_y_yzz_xxxyyy_1[i] * pc_x[i];

        ta1_y_xyzz_xxxyyz_0[i] =
            3.0 * ta1_y_yzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxyyz_1[i] * fe_0 + ta1_y_yzz_xxxyyz_0[i] * pa_x[i] - ta1_y_yzz_xxxyyz_1[i] * pc_x[i];

        ta1_y_xyzz_xxxyzz_0[i] =
            3.0 * ta1_y_yzz_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxyzz_1[i] * fe_0 + ta1_y_yzz_xxxyzz_0[i] * pa_x[i] - ta1_y_yzz_xxxyzz_1[i] * pc_x[i];

        ta1_y_xyzz_xxxzzz_0[i] = ta_xzz_xxxzzz_1[i] + ta1_y_xzz_xxxzzz_0[i] * pa_y[i] - ta1_y_xzz_xxxzzz_1[i] * pc_y[i];

        ta1_y_xyzz_xxyyyy_0[i] =
            2.0 * ta1_y_yzz_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyyyy_1[i] * fe_0 + ta1_y_yzz_xxyyyy_0[i] * pa_x[i] - ta1_y_yzz_xxyyyy_1[i] * pc_x[i];

        ta1_y_xyzz_xxyyyz_0[i] =
            2.0 * ta1_y_yzz_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyyyz_1[i] * fe_0 + ta1_y_yzz_xxyyyz_0[i] * pa_x[i] - ta1_y_yzz_xxyyyz_1[i] * pc_x[i];

        ta1_y_xyzz_xxyyzz_0[i] =
            2.0 * ta1_y_yzz_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyyzz_1[i] * fe_0 + ta1_y_yzz_xxyyzz_0[i] * pa_x[i] - ta1_y_yzz_xxyyzz_1[i] * pc_x[i];

        ta1_y_xyzz_xxyzzz_0[i] =
            2.0 * ta1_y_yzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyzzz_1[i] * fe_0 + ta1_y_yzz_xxyzzz_0[i] * pa_x[i] - ta1_y_yzz_xxyzzz_1[i] * pc_x[i];

        ta1_y_xyzz_xxzzzz_0[i] = ta_xzz_xxzzzz_1[i] + ta1_y_xzz_xxzzzz_0[i] * pa_y[i] - ta1_y_xzz_xxzzzz_1[i] * pc_y[i];

        ta1_y_xyzz_xyyyyy_0[i] =
            ta1_y_yzz_yyyyy_0[i] * fe_0 - ta1_y_yzz_yyyyy_1[i] * fe_0 + ta1_y_yzz_xyyyyy_0[i] * pa_x[i] - ta1_y_yzz_xyyyyy_1[i] * pc_x[i];

        ta1_y_xyzz_xyyyyz_0[i] =
            ta1_y_yzz_yyyyz_0[i] * fe_0 - ta1_y_yzz_yyyyz_1[i] * fe_0 + ta1_y_yzz_xyyyyz_0[i] * pa_x[i] - ta1_y_yzz_xyyyyz_1[i] * pc_x[i];

        ta1_y_xyzz_xyyyzz_0[i] =
            ta1_y_yzz_yyyzz_0[i] * fe_0 - ta1_y_yzz_yyyzz_1[i] * fe_0 + ta1_y_yzz_xyyyzz_0[i] * pa_x[i] - ta1_y_yzz_xyyyzz_1[i] * pc_x[i];

        ta1_y_xyzz_xyyzzz_0[i] =
            ta1_y_yzz_yyzzz_0[i] * fe_0 - ta1_y_yzz_yyzzz_1[i] * fe_0 + ta1_y_yzz_xyyzzz_0[i] * pa_x[i] - ta1_y_yzz_xyyzzz_1[i] * pc_x[i];

        ta1_y_xyzz_xyzzzz_0[i] =
            ta1_y_yzz_yzzzz_0[i] * fe_0 - ta1_y_yzz_yzzzz_1[i] * fe_0 + ta1_y_yzz_xyzzzz_0[i] * pa_x[i] - ta1_y_yzz_xyzzzz_1[i] * pc_x[i];

        ta1_y_xyzz_xzzzzz_0[i] = ta_xzz_xzzzzz_1[i] + ta1_y_xzz_xzzzzz_0[i] * pa_y[i] - ta1_y_xzz_xzzzzz_1[i] * pc_y[i];

        ta1_y_xyzz_yyyyyy_0[i] = ta1_y_yzz_yyyyyy_0[i] * pa_x[i] - ta1_y_yzz_yyyyyy_1[i] * pc_x[i];

        ta1_y_xyzz_yyyyyz_0[i] = ta1_y_yzz_yyyyyz_0[i] * pa_x[i] - ta1_y_yzz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xyzz_yyyyzz_0[i] = ta1_y_yzz_yyyyzz_0[i] * pa_x[i] - ta1_y_yzz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xyzz_yyyzzz_0[i] = ta1_y_yzz_yyyzzz_0[i] * pa_x[i] - ta1_y_yzz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xyzz_yyzzzz_0[i] = ta1_y_yzz_yyzzzz_0[i] * pa_x[i] - ta1_y_yzz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xyzz_yzzzzz_0[i] = ta1_y_yzz_yzzzzz_0[i] * pa_x[i] - ta1_y_yzz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xyzz_zzzzzz_0[i] = ta1_y_yzz_zzzzzz_0[i] * pa_x[i] - ta1_y_yzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 672-700 components of targeted buffer : GI

    auto ta1_y_xzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 672);

    auto ta1_y_xzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 673);

    auto ta1_y_xzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 674);

    auto ta1_y_xzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 675);

    auto ta1_y_xzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 676);

    auto ta1_y_xzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 677);

    auto ta1_y_xzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 678);

    auto ta1_y_xzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 679);

    auto ta1_y_xzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 680);

    auto ta1_y_xzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 681);

    auto ta1_y_xzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 682);

    auto ta1_y_xzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 683);

    auto ta1_y_xzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 684);

    auto ta1_y_xzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 685);

    auto ta1_y_xzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 686);

    auto ta1_y_xzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 687);

    auto ta1_y_xzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 688);

    auto ta1_y_xzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 689);

    auto ta1_y_xzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 690);

    auto ta1_y_xzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 691);

    auto ta1_y_xzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 692);

    auto ta1_y_xzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 693);

    auto ta1_y_xzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 694);

    auto ta1_y_xzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 695);

    auto ta1_y_xzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 696);

    auto ta1_y_xzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 697);

    auto ta1_y_xzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 698);

    auto ta1_y_xzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 699);

#pragma omp simd aligned(pa_x,                    \
                             pc_x,                \
                             ta1_y_xzzz_xxxxxx_0, \
                             ta1_y_xzzz_xxxxxy_0, \
                             ta1_y_xzzz_xxxxxz_0, \
                             ta1_y_xzzz_xxxxyy_0, \
                             ta1_y_xzzz_xxxxyz_0, \
                             ta1_y_xzzz_xxxxzz_0, \
                             ta1_y_xzzz_xxxyyy_0, \
                             ta1_y_xzzz_xxxyyz_0, \
                             ta1_y_xzzz_xxxyzz_0, \
                             ta1_y_xzzz_xxxzzz_0, \
                             ta1_y_xzzz_xxyyyy_0, \
                             ta1_y_xzzz_xxyyyz_0, \
                             ta1_y_xzzz_xxyyzz_0, \
                             ta1_y_xzzz_xxyzzz_0, \
                             ta1_y_xzzz_xxzzzz_0, \
                             ta1_y_xzzz_xyyyyy_0, \
                             ta1_y_xzzz_xyyyyz_0, \
                             ta1_y_xzzz_xyyyzz_0, \
                             ta1_y_xzzz_xyyzzz_0, \
                             ta1_y_xzzz_xyzzzz_0, \
                             ta1_y_xzzz_xzzzzz_0, \
                             ta1_y_xzzz_yyyyyy_0, \
                             ta1_y_xzzz_yyyyyz_0, \
                             ta1_y_xzzz_yyyyzz_0, \
                             ta1_y_xzzz_yyyzzz_0, \
                             ta1_y_xzzz_yyzzzz_0, \
                             ta1_y_xzzz_yzzzzz_0, \
                             ta1_y_xzzz_zzzzzz_0, \
                             ta1_y_zzz_xxxxx_0,   \
                             ta1_y_zzz_xxxxx_1,   \
                             ta1_y_zzz_xxxxxx_0,  \
                             ta1_y_zzz_xxxxxx_1,  \
                             ta1_y_zzz_xxxxxy_0,  \
                             ta1_y_zzz_xxxxxy_1,  \
                             ta1_y_zzz_xxxxxz_0,  \
                             ta1_y_zzz_xxxxxz_1,  \
                             ta1_y_zzz_xxxxy_0,   \
                             ta1_y_zzz_xxxxy_1,   \
                             ta1_y_zzz_xxxxyy_0,  \
                             ta1_y_zzz_xxxxyy_1,  \
                             ta1_y_zzz_xxxxyz_0,  \
                             ta1_y_zzz_xxxxyz_1,  \
                             ta1_y_zzz_xxxxz_0,   \
                             ta1_y_zzz_xxxxz_1,   \
                             ta1_y_zzz_xxxxzz_0,  \
                             ta1_y_zzz_xxxxzz_1,  \
                             ta1_y_zzz_xxxyy_0,   \
                             ta1_y_zzz_xxxyy_1,   \
                             ta1_y_zzz_xxxyyy_0,  \
                             ta1_y_zzz_xxxyyy_1,  \
                             ta1_y_zzz_xxxyyz_0,  \
                             ta1_y_zzz_xxxyyz_1,  \
                             ta1_y_zzz_xxxyz_0,   \
                             ta1_y_zzz_xxxyz_1,   \
                             ta1_y_zzz_xxxyzz_0,  \
                             ta1_y_zzz_xxxyzz_1,  \
                             ta1_y_zzz_xxxzz_0,   \
                             ta1_y_zzz_xxxzz_1,   \
                             ta1_y_zzz_xxxzzz_0,  \
                             ta1_y_zzz_xxxzzz_1,  \
                             ta1_y_zzz_xxyyy_0,   \
                             ta1_y_zzz_xxyyy_1,   \
                             ta1_y_zzz_xxyyyy_0,  \
                             ta1_y_zzz_xxyyyy_1,  \
                             ta1_y_zzz_xxyyyz_0,  \
                             ta1_y_zzz_xxyyyz_1,  \
                             ta1_y_zzz_xxyyz_0,   \
                             ta1_y_zzz_xxyyz_1,   \
                             ta1_y_zzz_xxyyzz_0,  \
                             ta1_y_zzz_xxyyzz_1,  \
                             ta1_y_zzz_xxyzz_0,   \
                             ta1_y_zzz_xxyzz_1,   \
                             ta1_y_zzz_xxyzzz_0,  \
                             ta1_y_zzz_xxyzzz_1,  \
                             ta1_y_zzz_xxzzz_0,   \
                             ta1_y_zzz_xxzzz_1,   \
                             ta1_y_zzz_xxzzzz_0,  \
                             ta1_y_zzz_xxzzzz_1,  \
                             ta1_y_zzz_xyyyy_0,   \
                             ta1_y_zzz_xyyyy_1,   \
                             ta1_y_zzz_xyyyyy_0,  \
                             ta1_y_zzz_xyyyyy_1,  \
                             ta1_y_zzz_xyyyyz_0,  \
                             ta1_y_zzz_xyyyyz_1,  \
                             ta1_y_zzz_xyyyz_0,   \
                             ta1_y_zzz_xyyyz_1,   \
                             ta1_y_zzz_xyyyzz_0,  \
                             ta1_y_zzz_xyyyzz_1,  \
                             ta1_y_zzz_xyyzz_0,   \
                             ta1_y_zzz_xyyzz_1,   \
                             ta1_y_zzz_xyyzzz_0,  \
                             ta1_y_zzz_xyyzzz_1,  \
                             ta1_y_zzz_xyzzz_0,   \
                             ta1_y_zzz_xyzzz_1,   \
                             ta1_y_zzz_xyzzzz_0,  \
                             ta1_y_zzz_xyzzzz_1,  \
                             ta1_y_zzz_xzzzz_0,   \
                             ta1_y_zzz_xzzzz_1,   \
                             ta1_y_zzz_xzzzzz_0,  \
                             ta1_y_zzz_xzzzzz_1,  \
                             ta1_y_zzz_yyyyy_0,   \
                             ta1_y_zzz_yyyyy_1,   \
                             ta1_y_zzz_yyyyyy_0,  \
                             ta1_y_zzz_yyyyyy_1,  \
                             ta1_y_zzz_yyyyyz_0,  \
                             ta1_y_zzz_yyyyyz_1,  \
                             ta1_y_zzz_yyyyz_0,   \
                             ta1_y_zzz_yyyyz_1,   \
                             ta1_y_zzz_yyyyzz_0,  \
                             ta1_y_zzz_yyyyzz_1,  \
                             ta1_y_zzz_yyyzz_0,   \
                             ta1_y_zzz_yyyzz_1,   \
                             ta1_y_zzz_yyyzzz_0,  \
                             ta1_y_zzz_yyyzzz_1,  \
                             ta1_y_zzz_yyzzz_0,   \
                             ta1_y_zzz_yyzzz_1,   \
                             ta1_y_zzz_yyzzzz_0,  \
                             ta1_y_zzz_yyzzzz_1,  \
                             ta1_y_zzz_yzzzz_0,   \
                             ta1_y_zzz_yzzzz_1,   \
                             ta1_y_zzz_yzzzzz_0,  \
                             ta1_y_zzz_yzzzzz_1,  \
                             ta1_y_zzz_zzzzz_0,   \
                             ta1_y_zzz_zzzzz_1,   \
                             ta1_y_zzz_zzzzzz_0,  \
                             ta1_y_zzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzz_xxxxxx_0[i] =
            6.0 * ta1_y_zzz_xxxxx_0[i] * fe_0 - 6.0 * ta1_y_zzz_xxxxx_1[i] * fe_0 + ta1_y_zzz_xxxxxx_0[i] * pa_x[i] - ta1_y_zzz_xxxxxx_1[i] * pc_x[i];

        ta1_y_xzzz_xxxxxy_0[i] =
            5.0 * ta1_y_zzz_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_zzz_xxxxy_1[i] * fe_0 + ta1_y_zzz_xxxxxy_0[i] * pa_x[i] - ta1_y_zzz_xxxxxy_1[i] * pc_x[i];

        ta1_y_xzzz_xxxxxz_0[i] =
            5.0 * ta1_y_zzz_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_zzz_xxxxz_1[i] * fe_0 + ta1_y_zzz_xxxxxz_0[i] * pa_x[i] - ta1_y_zzz_xxxxxz_1[i] * pc_x[i];

        ta1_y_xzzz_xxxxyy_0[i] =
            4.0 * ta1_y_zzz_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxyy_1[i] * fe_0 + ta1_y_zzz_xxxxyy_0[i] * pa_x[i] - ta1_y_zzz_xxxxyy_1[i] * pc_x[i];

        ta1_y_xzzz_xxxxyz_0[i] =
            4.0 * ta1_y_zzz_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxyz_1[i] * fe_0 + ta1_y_zzz_xxxxyz_0[i] * pa_x[i] - ta1_y_zzz_xxxxyz_1[i] * pc_x[i];

        ta1_y_xzzz_xxxxzz_0[i] =
            4.0 * ta1_y_zzz_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxzz_1[i] * fe_0 + ta1_y_zzz_xxxxzz_0[i] * pa_x[i] - ta1_y_zzz_xxxxzz_1[i] * pc_x[i];

        ta1_y_xzzz_xxxyyy_0[i] =
            3.0 * ta1_y_zzz_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxyyy_1[i] * fe_0 + ta1_y_zzz_xxxyyy_0[i] * pa_x[i] - ta1_y_zzz_xxxyyy_1[i] * pc_x[i];

        ta1_y_xzzz_xxxyyz_0[i] =
            3.0 * ta1_y_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxyyz_1[i] * fe_0 + ta1_y_zzz_xxxyyz_0[i] * pa_x[i] - ta1_y_zzz_xxxyyz_1[i] * pc_x[i];

        ta1_y_xzzz_xxxyzz_0[i] =
            3.0 * ta1_y_zzz_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxyzz_1[i] * fe_0 + ta1_y_zzz_xxxyzz_0[i] * pa_x[i] - ta1_y_zzz_xxxyzz_1[i] * pc_x[i];

        ta1_y_xzzz_xxxzzz_0[i] =
            3.0 * ta1_y_zzz_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxzzz_1[i] * fe_0 + ta1_y_zzz_xxxzzz_0[i] * pa_x[i] - ta1_y_zzz_xxxzzz_1[i] * pc_x[i];

        ta1_y_xzzz_xxyyyy_0[i] =
            2.0 * ta1_y_zzz_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyyyy_1[i] * fe_0 + ta1_y_zzz_xxyyyy_0[i] * pa_x[i] - ta1_y_zzz_xxyyyy_1[i] * pc_x[i];

        ta1_y_xzzz_xxyyyz_0[i] =
            2.0 * ta1_y_zzz_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyyyz_1[i] * fe_0 + ta1_y_zzz_xxyyyz_0[i] * pa_x[i] - ta1_y_zzz_xxyyyz_1[i] * pc_x[i];

        ta1_y_xzzz_xxyyzz_0[i] =
            2.0 * ta1_y_zzz_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyyzz_1[i] * fe_0 + ta1_y_zzz_xxyyzz_0[i] * pa_x[i] - ta1_y_zzz_xxyyzz_1[i] * pc_x[i];

        ta1_y_xzzz_xxyzzz_0[i] =
            2.0 * ta1_y_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyzzz_1[i] * fe_0 + ta1_y_zzz_xxyzzz_0[i] * pa_x[i] - ta1_y_zzz_xxyzzz_1[i] * pc_x[i];

        ta1_y_xzzz_xxzzzz_0[i] =
            2.0 * ta1_y_zzz_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xzzzz_1[i] * fe_0 + ta1_y_zzz_xxzzzz_0[i] * pa_x[i] - ta1_y_zzz_xxzzzz_1[i] * pc_x[i];

        ta1_y_xzzz_xyyyyy_0[i] =
            ta1_y_zzz_yyyyy_0[i] * fe_0 - ta1_y_zzz_yyyyy_1[i] * fe_0 + ta1_y_zzz_xyyyyy_0[i] * pa_x[i] - ta1_y_zzz_xyyyyy_1[i] * pc_x[i];

        ta1_y_xzzz_xyyyyz_0[i] =
            ta1_y_zzz_yyyyz_0[i] * fe_0 - ta1_y_zzz_yyyyz_1[i] * fe_0 + ta1_y_zzz_xyyyyz_0[i] * pa_x[i] - ta1_y_zzz_xyyyyz_1[i] * pc_x[i];

        ta1_y_xzzz_xyyyzz_0[i] =
            ta1_y_zzz_yyyzz_0[i] * fe_0 - ta1_y_zzz_yyyzz_1[i] * fe_0 + ta1_y_zzz_xyyyzz_0[i] * pa_x[i] - ta1_y_zzz_xyyyzz_1[i] * pc_x[i];

        ta1_y_xzzz_xyyzzz_0[i] =
            ta1_y_zzz_yyzzz_0[i] * fe_0 - ta1_y_zzz_yyzzz_1[i] * fe_0 + ta1_y_zzz_xyyzzz_0[i] * pa_x[i] - ta1_y_zzz_xyyzzz_1[i] * pc_x[i];

        ta1_y_xzzz_xyzzzz_0[i] =
            ta1_y_zzz_yzzzz_0[i] * fe_0 - ta1_y_zzz_yzzzz_1[i] * fe_0 + ta1_y_zzz_xyzzzz_0[i] * pa_x[i] - ta1_y_zzz_xyzzzz_1[i] * pc_x[i];

        ta1_y_xzzz_xzzzzz_0[i] =
            ta1_y_zzz_zzzzz_0[i] * fe_0 - ta1_y_zzz_zzzzz_1[i] * fe_0 + ta1_y_zzz_xzzzzz_0[i] * pa_x[i] - ta1_y_zzz_xzzzzz_1[i] * pc_x[i];

        ta1_y_xzzz_yyyyyy_0[i] = ta1_y_zzz_yyyyyy_0[i] * pa_x[i] - ta1_y_zzz_yyyyyy_1[i] * pc_x[i];

        ta1_y_xzzz_yyyyyz_0[i] = ta1_y_zzz_yyyyyz_0[i] * pa_x[i] - ta1_y_zzz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xzzz_yyyyzz_0[i] = ta1_y_zzz_yyyyzz_0[i] * pa_x[i] - ta1_y_zzz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xzzz_yyyzzz_0[i] = ta1_y_zzz_yyyzzz_0[i] * pa_x[i] - ta1_y_zzz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xzzz_yyzzzz_0[i] = ta1_y_zzz_yyzzzz_0[i] * pa_x[i] - ta1_y_zzz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xzzz_yzzzzz_0[i] = ta1_y_zzz_yzzzzz_0[i] * pa_x[i] - ta1_y_zzz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xzzz_zzzzzz_0[i] = ta1_y_zzz_zzzzzz_0[i] * pa_x[i] - ta1_y_zzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 700-728 components of targeted buffer : GI

    auto ta1_y_yyyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 700);

    auto ta1_y_yyyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 701);

    auto ta1_y_yyyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 702);

    auto ta1_y_yyyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 703);

    auto ta1_y_yyyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 704);

    auto ta1_y_yyyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 705);

    auto ta1_y_yyyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 706);

    auto ta1_y_yyyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 707);

    auto ta1_y_yyyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 708);

    auto ta1_y_yyyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 709);

    auto ta1_y_yyyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 710);

    auto ta1_y_yyyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 711);

    auto ta1_y_yyyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 712);

    auto ta1_y_yyyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 713);

    auto ta1_y_yyyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 714);

    auto ta1_y_yyyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 715);

    auto ta1_y_yyyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 716);

    auto ta1_y_yyyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 717);

    auto ta1_y_yyyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 718);

    auto ta1_y_yyyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 719);

    auto ta1_y_yyyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 720);

    auto ta1_y_yyyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 721);

    auto ta1_y_yyyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 722);

    auto ta1_y_yyyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 723);

    auto ta1_y_yyyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 724);

    auto ta1_y_yyyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 725);

    auto ta1_y_yyyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 726);

    auto ta1_y_yyyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 727);

#pragma omp simd aligned(pa_y,                    \
                             pc_y,                \
                             ta1_y_yy_xxxxxx_0,   \
                             ta1_y_yy_xxxxxx_1,   \
                             ta1_y_yy_xxxxxy_0,   \
                             ta1_y_yy_xxxxxy_1,   \
                             ta1_y_yy_xxxxxz_0,   \
                             ta1_y_yy_xxxxxz_1,   \
                             ta1_y_yy_xxxxyy_0,   \
                             ta1_y_yy_xxxxyy_1,   \
                             ta1_y_yy_xxxxyz_0,   \
                             ta1_y_yy_xxxxyz_1,   \
                             ta1_y_yy_xxxxzz_0,   \
                             ta1_y_yy_xxxxzz_1,   \
                             ta1_y_yy_xxxyyy_0,   \
                             ta1_y_yy_xxxyyy_1,   \
                             ta1_y_yy_xxxyyz_0,   \
                             ta1_y_yy_xxxyyz_1,   \
                             ta1_y_yy_xxxyzz_0,   \
                             ta1_y_yy_xxxyzz_1,   \
                             ta1_y_yy_xxxzzz_0,   \
                             ta1_y_yy_xxxzzz_1,   \
                             ta1_y_yy_xxyyyy_0,   \
                             ta1_y_yy_xxyyyy_1,   \
                             ta1_y_yy_xxyyyz_0,   \
                             ta1_y_yy_xxyyyz_1,   \
                             ta1_y_yy_xxyyzz_0,   \
                             ta1_y_yy_xxyyzz_1,   \
                             ta1_y_yy_xxyzzz_0,   \
                             ta1_y_yy_xxyzzz_1,   \
                             ta1_y_yy_xxzzzz_0,   \
                             ta1_y_yy_xxzzzz_1,   \
                             ta1_y_yy_xyyyyy_0,   \
                             ta1_y_yy_xyyyyy_1,   \
                             ta1_y_yy_xyyyyz_0,   \
                             ta1_y_yy_xyyyyz_1,   \
                             ta1_y_yy_xyyyzz_0,   \
                             ta1_y_yy_xyyyzz_1,   \
                             ta1_y_yy_xyyzzz_0,   \
                             ta1_y_yy_xyyzzz_1,   \
                             ta1_y_yy_xyzzzz_0,   \
                             ta1_y_yy_xyzzzz_1,   \
                             ta1_y_yy_xzzzzz_0,   \
                             ta1_y_yy_xzzzzz_1,   \
                             ta1_y_yy_yyyyyy_0,   \
                             ta1_y_yy_yyyyyy_1,   \
                             ta1_y_yy_yyyyyz_0,   \
                             ta1_y_yy_yyyyyz_1,   \
                             ta1_y_yy_yyyyzz_0,   \
                             ta1_y_yy_yyyyzz_1,   \
                             ta1_y_yy_yyyzzz_0,   \
                             ta1_y_yy_yyyzzz_1,   \
                             ta1_y_yy_yyzzzz_0,   \
                             ta1_y_yy_yyzzzz_1,   \
                             ta1_y_yy_yzzzzz_0,   \
                             ta1_y_yy_yzzzzz_1,   \
                             ta1_y_yy_zzzzzz_0,   \
                             ta1_y_yy_zzzzzz_1,   \
                             ta1_y_yyy_xxxxx_0,   \
                             ta1_y_yyy_xxxxx_1,   \
                             ta1_y_yyy_xxxxxx_0,  \
                             ta1_y_yyy_xxxxxx_1,  \
                             ta1_y_yyy_xxxxxy_0,  \
                             ta1_y_yyy_xxxxxy_1,  \
                             ta1_y_yyy_xxxxxz_0,  \
                             ta1_y_yyy_xxxxxz_1,  \
                             ta1_y_yyy_xxxxy_0,   \
                             ta1_y_yyy_xxxxy_1,   \
                             ta1_y_yyy_xxxxyy_0,  \
                             ta1_y_yyy_xxxxyy_1,  \
                             ta1_y_yyy_xxxxyz_0,  \
                             ta1_y_yyy_xxxxyz_1,  \
                             ta1_y_yyy_xxxxz_0,   \
                             ta1_y_yyy_xxxxz_1,   \
                             ta1_y_yyy_xxxxzz_0,  \
                             ta1_y_yyy_xxxxzz_1,  \
                             ta1_y_yyy_xxxyy_0,   \
                             ta1_y_yyy_xxxyy_1,   \
                             ta1_y_yyy_xxxyyy_0,  \
                             ta1_y_yyy_xxxyyy_1,  \
                             ta1_y_yyy_xxxyyz_0,  \
                             ta1_y_yyy_xxxyyz_1,  \
                             ta1_y_yyy_xxxyz_0,   \
                             ta1_y_yyy_xxxyz_1,   \
                             ta1_y_yyy_xxxyzz_0,  \
                             ta1_y_yyy_xxxyzz_1,  \
                             ta1_y_yyy_xxxzz_0,   \
                             ta1_y_yyy_xxxzz_1,   \
                             ta1_y_yyy_xxxzzz_0,  \
                             ta1_y_yyy_xxxzzz_1,  \
                             ta1_y_yyy_xxyyy_0,   \
                             ta1_y_yyy_xxyyy_1,   \
                             ta1_y_yyy_xxyyyy_0,  \
                             ta1_y_yyy_xxyyyy_1,  \
                             ta1_y_yyy_xxyyyz_0,  \
                             ta1_y_yyy_xxyyyz_1,  \
                             ta1_y_yyy_xxyyz_0,   \
                             ta1_y_yyy_xxyyz_1,   \
                             ta1_y_yyy_xxyyzz_0,  \
                             ta1_y_yyy_xxyyzz_1,  \
                             ta1_y_yyy_xxyzz_0,   \
                             ta1_y_yyy_xxyzz_1,   \
                             ta1_y_yyy_xxyzzz_0,  \
                             ta1_y_yyy_xxyzzz_1,  \
                             ta1_y_yyy_xxzzz_0,   \
                             ta1_y_yyy_xxzzz_1,   \
                             ta1_y_yyy_xxzzzz_0,  \
                             ta1_y_yyy_xxzzzz_1,  \
                             ta1_y_yyy_xyyyy_0,   \
                             ta1_y_yyy_xyyyy_1,   \
                             ta1_y_yyy_xyyyyy_0,  \
                             ta1_y_yyy_xyyyyy_1,  \
                             ta1_y_yyy_xyyyyz_0,  \
                             ta1_y_yyy_xyyyyz_1,  \
                             ta1_y_yyy_xyyyz_0,   \
                             ta1_y_yyy_xyyyz_1,   \
                             ta1_y_yyy_xyyyzz_0,  \
                             ta1_y_yyy_xyyyzz_1,  \
                             ta1_y_yyy_xyyzz_0,   \
                             ta1_y_yyy_xyyzz_1,   \
                             ta1_y_yyy_xyyzzz_0,  \
                             ta1_y_yyy_xyyzzz_1,  \
                             ta1_y_yyy_xyzzz_0,   \
                             ta1_y_yyy_xyzzz_1,   \
                             ta1_y_yyy_xyzzzz_0,  \
                             ta1_y_yyy_xyzzzz_1,  \
                             ta1_y_yyy_xzzzz_0,   \
                             ta1_y_yyy_xzzzz_1,   \
                             ta1_y_yyy_xzzzzz_0,  \
                             ta1_y_yyy_xzzzzz_1,  \
                             ta1_y_yyy_yyyyy_0,   \
                             ta1_y_yyy_yyyyy_1,   \
                             ta1_y_yyy_yyyyyy_0,  \
                             ta1_y_yyy_yyyyyy_1,  \
                             ta1_y_yyy_yyyyyz_0,  \
                             ta1_y_yyy_yyyyyz_1,  \
                             ta1_y_yyy_yyyyz_0,   \
                             ta1_y_yyy_yyyyz_1,   \
                             ta1_y_yyy_yyyyzz_0,  \
                             ta1_y_yyy_yyyyzz_1,  \
                             ta1_y_yyy_yyyzz_0,   \
                             ta1_y_yyy_yyyzz_1,   \
                             ta1_y_yyy_yyyzzz_0,  \
                             ta1_y_yyy_yyyzzz_1,  \
                             ta1_y_yyy_yyzzz_0,   \
                             ta1_y_yyy_yyzzz_1,   \
                             ta1_y_yyy_yyzzzz_0,  \
                             ta1_y_yyy_yyzzzz_1,  \
                             ta1_y_yyy_yzzzz_0,   \
                             ta1_y_yyy_yzzzz_1,   \
                             ta1_y_yyy_yzzzzz_0,  \
                             ta1_y_yyy_yzzzzz_1,  \
                             ta1_y_yyy_zzzzz_0,   \
                             ta1_y_yyy_zzzzz_1,   \
                             ta1_y_yyy_zzzzzz_0,  \
                             ta1_y_yyy_zzzzzz_1,  \
                             ta1_y_yyyy_xxxxxx_0, \
                             ta1_y_yyyy_xxxxxy_0, \
                             ta1_y_yyyy_xxxxxz_0, \
                             ta1_y_yyyy_xxxxyy_0, \
                             ta1_y_yyyy_xxxxyz_0, \
                             ta1_y_yyyy_xxxxzz_0, \
                             ta1_y_yyyy_xxxyyy_0, \
                             ta1_y_yyyy_xxxyyz_0, \
                             ta1_y_yyyy_xxxyzz_0, \
                             ta1_y_yyyy_xxxzzz_0, \
                             ta1_y_yyyy_xxyyyy_0, \
                             ta1_y_yyyy_xxyyyz_0, \
                             ta1_y_yyyy_xxyyzz_0, \
                             ta1_y_yyyy_xxyzzz_0, \
                             ta1_y_yyyy_xxzzzz_0, \
                             ta1_y_yyyy_xyyyyy_0, \
                             ta1_y_yyyy_xyyyyz_0, \
                             ta1_y_yyyy_xyyyzz_0, \
                             ta1_y_yyyy_xyyzzz_0, \
                             ta1_y_yyyy_xyzzzz_0, \
                             ta1_y_yyyy_xzzzzz_0, \
                             ta1_y_yyyy_yyyyyy_0, \
                             ta1_y_yyyy_yyyyyz_0, \
                             ta1_y_yyyy_yyyyzz_0, \
                             ta1_y_yyyy_yyyzzz_0, \
                             ta1_y_yyyy_yyzzzz_0, \
                             ta1_y_yyyy_yzzzzz_0, \
                             ta1_y_yyyy_zzzzzz_0, \
                             ta_yyy_xxxxxx_1,     \
                             ta_yyy_xxxxxy_1,     \
                             ta_yyy_xxxxxz_1,     \
                             ta_yyy_xxxxyy_1,     \
                             ta_yyy_xxxxyz_1,     \
                             ta_yyy_xxxxzz_1,     \
                             ta_yyy_xxxyyy_1,     \
                             ta_yyy_xxxyyz_1,     \
                             ta_yyy_xxxyzz_1,     \
                             ta_yyy_xxxzzz_1,     \
                             ta_yyy_xxyyyy_1,     \
                             ta_yyy_xxyyyz_1,     \
                             ta_yyy_xxyyzz_1,     \
                             ta_yyy_xxyzzz_1,     \
                             ta_yyy_xxzzzz_1,     \
                             ta_yyy_xyyyyy_1,     \
                             ta_yyy_xyyyyz_1,     \
                             ta_yyy_xyyyzz_1,     \
                             ta_yyy_xyyzzz_1,     \
                             ta_yyy_xyzzzz_1,     \
                             ta_yyy_xzzzzz_1,     \
                             ta_yyy_yyyyyy_1,     \
                             ta_yyy_yyyyyz_1,     \
                             ta_yyy_yyyyzz_1,     \
                             ta_yyy_yyyzzz_1,     \
                             ta_yyy_yyzzzz_1,     \
                             ta_yyy_yzzzzz_1,     \
                             ta_yyy_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyy_xxxxxx_0[i] = 3.0 * ta1_y_yy_xxxxxx_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxxx_1[i] * fe_0 + ta_yyy_xxxxxx_1[i] +
                                 ta1_y_yyy_xxxxxx_0[i] * pa_y[i] - ta1_y_yyy_xxxxxx_1[i] * pc_y[i];

        ta1_y_yyyy_xxxxxy_0[i] = 3.0 * ta1_y_yy_xxxxxy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxxy_1[i] * fe_0 + ta1_y_yyy_xxxxx_0[i] * fe_0 -
                                 ta1_y_yyy_xxxxx_1[i] * fe_0 + ta_yyy_xxxxxy_1[i] + ta1_y_yyy_xxxxxy_0[i] * pa_y[i] - ta1_y_yyy_xxxxxy_1[i] * pc_y[i];

        ta1_y_yyyy_xxxxxz_0[i] = 3.0 * ta1_y_yy_xxxxxz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxxz_1[i] * fe_0 + ta_yyy_xxxxxz_1[i] +
                                 ta1_y_yyy_xxxxxz_0[i] * pa_y[i] - ta1_y_yyy_xxxxxz_1[i] * pc_y[i];

        ta1_y_yyyy_xxxxyy_0[i] = 3.0 * ta1_y_yy_xxxxyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxyy_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxxxy_0[i] * fe_0 -
                                 2.0 * ta1_y_yyy_xxxxy_1[i] * fe_0 + ta_yyy_xxxxyy_1[i] + ta1_y_yyy_xxxxyy_0[i] * pa_y[i] -
                                 ta1_y_yyy_xxxxyy_1[i] * pc_y[i];

        ta1_y_yyyy_xxxxyz_0[i] = 3.0 * ta1_y_yy_xxxxyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxyz_1[i] * fe_0 + ta1_y_yyy_xxxxz_0[i] * fe_0 -
                                 ta1_y_yyy_xxxxz_1[i] * fe_0 + ta_yyy_xxxxyz_1[i] + ta1_y_yyy_xxxxyz_0[i] * pa_y[i] - ta1_y_yyy_xxxxyz_1[i] * pc_y[i];

        ta1_y_yyyy_xxxxzz_0[i] = 3.0 * ta1_y_yy_xxxxzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxzz_1[i] * fe_0 + ta_yyy_xxxxzz_1[i] +
                                 ta1_y_yyy_xxxxzz_0[i] * pa_y[i] - ta1_y_yyy_xxxxzz_1[i] * pc_y[i];

        ta1_y_yyyy_xxxyyy_0[i] = 3.0 * ta1_y_yy_xxxyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxyyy_1[i] * fe_0 + 3.0 * ta1_y_yyy_xxxyy_0[i] * fe_0 -
                                 3.0 * ta1_y_yyy_xxxyy_1[i] * fe_0 + ta_yyy_xxxyyy_1[i] + ta1_y_yyy_xxxyyy_0[i] * pa_y[i] -
                                 ta1_y_yyy_xxxyyy_1[i] * pc_y[i];

        ta1_y_yyyy_xxxyyz_0[i] = 3.0 * ta1_y_yy_xxxyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxyyz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyy_xxxyz_1[i] * fe_0 + ta_yyy_xxxyyz_1[i] + ta1_y_yyy_xxxyyz_0[i] * pa_y[i] -
                                 ta1_y_yyy_xxxyyz_1[i] * pc_y[i];

        ta1_y_yyyy_xxxyzz_0[i] = 3.0 * ta1_y_yy_xxxyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxyzz_1[i] * fe_0 + ta1_y_yyy_xxxzz_0[i] * fe_0 -
                                 ta1_y_yyy_xxxzz_1[i] * fe_0 + ta_yyy_xxxyzz_1[i] + ta1_y_yyy_xxxyzz_0[i] * pa_y[i] - ta1_y_yyy_xxxyzz_1[i] * pc_y[i];

        ta1_y_yyyy_xxxzzz_0[i] = 3.0 * ta1_y_yy_xxxzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxzzz_1[i] * fe_0 + ta_yyy_xxxzzz_1[i] +
                                 ta1_y_yyy_xxxzzz_0[i] * pa_y[i] - ta1_y_yyy_xxxzzz_1[i] * pc_y[i];

        ta1_y_yyyy_xxyyyy_0[i] = 3.0 * ta1_y_yy_xxyyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyyyy_1[i] * fe_0 + 4.0 * ta1_y_yyy_xxyyy_0[i] * fe_0 -
                                 4.0 * ta1_y_yyy_xxyyy_1[i] * fe_0 + ta_yyy_xxyyyy_1[i] + ta1_y_yyy_xxyyyy_0[i] * pa_y[i] -
                                 ta1_y_yyy_xxyyyy_1[i] * pc_y[i];

        ta1_y_yyyy_xxyyyz_0[i] = 3.0 * ta1_y_yy_xxyyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyyyz_1[i] * fe_0 + 3.0 * ta1_y_yyy_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_y_yyy_xxyyz_1[i] * fe_0 + ta_yyy_xxyyyz_1[i] + ta1_y_yyy_xxyyyz_0[i] * pa_y[i] -
                                 ta1_y_yyy_xxyyyz_1[i] * pc_y[i];

        ta1_y_yyyy_xxyyzz_0[i] = 3.0 * ta1_y_yy_xxyyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxyzz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyy_xxyzz_1[i] * fe_0 + ta_yyy_xxyyzz_1[i] + ta1_y_yyy_xxyyzz_0[i] * pa_y[i] -
                                 ta1_y_yyy_xxyyzz_1[i] * pc_y[i];

        ta1_y_yyyy_xxyzzz_0[i] = 3.0 * ta1_y_yy_xxyzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyzzz_1[i] * fe_0 + ta1_y_yyy_xxzzz_0[i] * fe_0 -
                                 ta1_y_yyy_xxzzz_1[i] * fe_0 + ta_yyy_xxyzzz_1[i] + ta1_y_yyy_xxyzzz_0[i] * pa_y[i] - ta1_y_yyy_xxyzzz_1[i] * pc_y[i];

        ta1_y_yyyy_xxzzzz_0[i] = 3.0 * ta1_y_yy_xxzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxzzzz_1[i] * fe_0 + ta_yyy_xxzzzz_1[i] +
                                 ta1_y_yyy_xxzzzz_0[i] * pa_y[i] - ta1_y_yyy_xxzzzz_1[i] * pc_y[i];

        ta1_y_yyyy_xyyyyy_0[i] = 3.0 * ta1_y_yy_xyyyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyyyy_1[i] * fe_0 + 5.0 * ta1_y_yyy_xyyyy_0[i] * fe_0 -
                                 5.0 * ta1_y_yyy_xyyyy_1[i] * fe_0 + ta_yyy_xyyyyy_1[i] + ta1_y_yyy_xyyyyy_0[i] * pa_y[i] -
                                 ta1_y_yyy_xyyyyy_1[i] * pc_y[i];

        ta1_y_yyyy_xyyyyz_0[i] = 3.0 * ta1_y_yy_xyyyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyyyz_1[i] * fe_0 + 4.0 * ta1_y_yyy_xyyyz_0[i] * fe_0 -
                                 4.0 * ta1_y_yyy_xyyyz_1[i] * fe_0 + ta_yyy_xyyyyz_1[i] + ta1_y_yyy_xyyyyz_0[i] * pa_y[i] -
                                 ta1_y_yyy_xyyyyz_1[i] * pc_y[i];

        ta1_y_yyyy_xyyyzz_0[i] = 3.0 * ta1_y_yy_xyyyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyyzz_1[i] * fe_0 + 3.0 * ta1_y_yyy_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_yyy_xyyzz_1[i] * fe_0 + ta_yyy_xyyyzz_1[i] + ta1_y_yyy_xyyyzz_0[i] * pa_y[i] -
                                 ta1_y_yyy_xyyyzz_1[i] * pc_y[i];

        ta1_y_yyyy_xyyzzz_0[i] = 3.0 * ta1_y_yy_xyyzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyzzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyy_xyzzz_1[i] * fe_0 + ta_yyy_xyyzzz_1[i] + ta1_y_yyy_xyyzzz_0[i] * pa_y[i] -
                                 ta1_y_yyy_xyyzzz_1[i] * pc_y[i];

        ta1_y_yyyy_xyzzzz_0[i] = 3.0 * ta1_y_yy_xyzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyzzzz_1[i] * fe_0 + ta1_y_yyy_xzzzz_0[i] * fe_0 -
                                 ta1_y_yyy_xzzzz_1[i] * fe_0 + ta_yyy_xyzzzz_1[i] + ta1_y_yyy_xyzzzz_0[i] * pa_y[i] - ta1_y_yyy_xyzzzz_1[i] * pc_y[i];

        ta1_y_yyyy_xzzzzz_0[i] = 3.0 * ta1_y_yy_xzzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xzzzzz_1[i] * fe_0 + ta_yyy_xzzzzz_1[i] +
                                 ta1_y_yyy_xzzzzz_0[i] * pa_y[i] - ta1_y_yyy_xzzzzz_1[i] * pc_y[i];

        ta1_y_yyyy_yyyyyy_0[i] = 3.0 * ta1_y_yy_yyyyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyyyy_1[i] * fe_0 + 6.0 * ta1_y_yyy_yyyyy_0[i] * fe_0 -
                                 6.0 * ta1_y_yyy_yyyyy_1[i] * fe_0 + ta_yyy_yyyyyy_1[i] + ta1_y_yyy_yyyyyy_0[i] * pa_y[i] -
                                 ta1_y_yyy_yyyyyy_1[i] * pc_y[i];

        ta1_y_yyyy_yyyyyz_0[i] = 3.0 * ta1_y_yy_yyyyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyyyz_1[i] * fe_0 + 5.0 * ta1_y_yyy_yyyyz_0[i] * fe_0 -
                                 5.0 * ta1_y_yyy_yyyyz_1[i] * fe_0 + ta_yyy_yyyyyz_1[i] + ta1_y_yyy_yyyyyz_0[i] * pa_y[i] -
                                 ta1_y_yyy_yyyyyz_1[i] * pc_y[i];

        ta1_y_yyyy_yyyyzz_0[i] = 3.0 * ta1_y_yy_yyyyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyyzz_1[i] * fe_0 + 4.0 * ta1_y_yyy_yyyzz_0[i] * fe_0 -
                                 4.0 * ta1_y_yyy_yyyzz_1[i] * fe_0 + ta_yyy_yyyyzz_1[i] + ta1_y_yyy_yyyyzz_0[i] * pa_y[i] -
                                 ta1_y_yyy_yyyyzz_1[i] * pc_y[i];

        ta1_y_yyyy_yyyzzz_0[i] = 3.0 * ta1_y_yy_yyyzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyzzz_1[i] * fe_0 + 3.0 * ta1_y_yyy_yyzzz_0[i] * fe_0 -
                                 3.0 * ta1_y_yyy_yyzzz_1[i] * fe_0 + ta_yyy_yyyzzz_1[i] + ta1_y_yyy_yyyzzz_0[i] * pa_y[i] -
                                 ta1_y_yyy_yyyzzz_1[i] * pc_y[i];

        ta1_y_yyyy_yyzzzz_0[i] = 3.0 * ta1_y_yy_yyzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyzzzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yzzzz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyy_yzzzz_1[i] * fe_0 + ta_yyy_yyzzzz_1[i] + ta1_y_yyy_yyzzzz_0[i] * pa_y[i] -
                                 ta1_y_yyy_yyzzzz_1[i] * pc_y[i];

        ta1_y_yyyy_yzzzzz_0[i] = 3.0 * ta1_y_yy_yzzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yzzzzz_1[i] * fe_0 + ta1_y_yyy_zzzzz_0[i] * fe_0 -
                                 ta1_y_yyy_zzzzz_1[i] * fe_0 + ta_yyy_yzzzzz_1[i] + ta1_y_yyy_yzzzzz_0[i] * pa_y[i] - ta1_y_yyy_yzzzzz_1[i] * pc_y[i];

        ta1_y_yyyy_zzzzzz_0[i] = 3.0 * ta1_y_yy_zzzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_zzzzzz_1[i] * fe_0 + ta_yyy_zzzzzz_1[i] +
                                 ta1_y_yyy_zzzzzz_0[i] * pa_y[i] - ta1_y_yyy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 728-756 components of targeted buffer : GI

    auto ta1_y_yyyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 728);

    auto ta1_y_yyyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 729);

    auto ta1_y_yyyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 730);

    auto ta1_y_yyyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 731);

    auto ta1_y_yyyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 732);

    auto ta1_y_yyyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 733);

    auto ta1_y_yyyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 734);

    auto ta1_y_yyyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 735);

    auto ta1_y_yyyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 736);

    auto ta1_y_yyyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 737);

    auto ta1_y_yyyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 738);

    auto ta1_y_yyyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 739);

    auto ta1_y_yyyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 740);

    auto ta1_y_yyyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 741);

    auto ta1_y_yyyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 742);

    auto ta1_y_yyyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 743);

    auto ta1_y_yyyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 744);

    auto ta1_y_yyyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 745);

    auto ta1_y_yyyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 746);

    auto ta1_y_yyyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 747);

    auto ta1_y_yyyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 748);

    auto ta1_y_yyyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 749);

    auto ta1_y_yyyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 750);

    auto ta1_y_yyyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 751);

    auto ta1_y_yyyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 752);

    auto ta1_y_yyyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 753);

    auto ta1_y_yyyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 754);

    auto ta1_y_yyyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 755);

#pragma omp simd aligned(pa_z,                    \
                             pc_z,                \
                             ta1_y_yyy_xxxxx_0,   \
                             ta1_y_yyy_xxxxx_1,   \
                             ta1_y_yyy_xxxxxx_0,  \
                             ta1_y_yyy_xxxxxx_1,  \
                             ta1_y_yyy_xxxxxy_0,  \
                             ta1_y_yyy_xxxxxy_1,  \
                             ta1_y_yyy_xxxxxz_0,  \
                             ta1_y_yyy_xxxxxz_1,  \
                             ta1_y_yyy_xxxxy_0,   \
                             ta1_y_yyy_xxxxy_1,   \
                             ta1_y_yyy_xxxxyy_0,  \
                             ta1_y_yyy_xxxxyy_1,  \
                             ta1_y_yyy_xxxxyz_0,  \
                             ta1_y_yyy_xxxxyz_1,  \
                             ta1_y_yyy_xxxxz_0,   \
                             ta1_y_yyy_xxxxz_1,   \
                             ta1_y_yyy_xxxxzz_0,  \
                             ta1_y_yyy_xxxxzz_1,  \
                             ta1_y_yyy_xxxyy_0,   \
                             ta1_y_yyy_xxxyy_1,   \
                             ta1_y_yyy_xxxyyy_0,  \
                             ta1_y_yyy_xxxyyy_1,  \
                             ta1_y_yyy_xxxyyz_0,  \
                             ta1_y_yyy_xxxyyz_1,  \
                             ta1_y_yyy_xxxyz_0,   \
                             ta1_y_yyy_xxxyz_1,   \
                             ta1_y_yyy_xxxyzz_0,  \
                             ta1_y_yyy_xxxyzz_1,  \
                             ta1_y_yyy_xxxzz_0,   \
                             ta1_y_yyy_xxxzz_1,   \
                             ta1_y_yyy_xxxzzz_0,  \
                             ta1_y_yyy_xxxzzz_1,  \
                             ta1_y_yyy_xxyyy_0,   \
                             ta1_y_yyy_xxyyy_1,   \
                             ta1_y_yyy_xxyyyy_0,  \
                             ta1_y_yyy_xxyyyy_1,  \
                             ta1_y_yyy_xxyyyz_0,  \
                             ta1_y_yyy_xxyyyz_1,  \
                             ta1_y_yyy_xxyyz_0,   \
                             ta1_y_yyy_xxyyz_1,   \
                             ta1_y_yyy_xxyyzz_0,  \
                             ta1_y_yyy_xxyyzz_1,  \
                             ta1_y_yyy_xxyzz_0,   \
                             ta1_y_yyy_xxyzz_1,   \
                             ta1_y_yyy_xxyzzz_0,  \
                             ta1_y_yyy_xxyzzz_1,  \
                             ta1_y_yyy_xxzzz_0,   \
                             ta1_y_yyy_xxzzz_1,   \
                             ta1_y_yyy_xxzzzz_0,  \
                             ta1_y_yyy_xxzzzz_1,  \
                             ta1_y_yyy_xyyyy_0,   \
                             ta1_y_yyy_xyyyy_1,   \
                             ta1_y_yyy_xyyyyy_0,  \
                             ta1_y_yyy_xyyyyy_1,  \
                             ta1_y_yyy_xyyyyz_0,  \
                             ta1_y_yyy_xyyyyz_1,  \
                             ta1_y_yyy_xyyyz_0,   \
                             ta1_y_yyy_xyyyz_1,   \
                             ta1_y_yyy_xyyyzz_0,  \
                             ta1_y_yyy_xyyyzz_1,  \
                             ta1_y_yyy_xyyzz_0,   \
                             ta1_y_yyy_xyyzz_1,   \
                             ta1_y_yyy_xyyzzz_0,  \
                             ta1_y_yyy_xyyzzz_1,  \
                             ta1_y_yyy_xyzzz_0,   \
                             ta1_y_yyy_xyzzz_1,   \
                             ta1_y_yyy_xyzzzz_0,  \
                             ta1_y_yyy_xyzzzz_1,  \
                             ta1_y_yyy_xzzzz_0,   \
                             ta1_y_yyy_xzzzz_1,   \
                             ta1_y_yyy_xzzzzz_0,  \
                             ta1_y_yyy_xzzzzz_1,  \
                             ta1_y_yyy_yyyyy_0,   \
                             ta1_y_yyy_yyyyy_1,   \
                             ta1_y_yyy_yyyyyy_0,  \
                             ta1_y_yyy_yyyyyy_1,  \
                             ta1_y_yyy_yyyyyz_0,  \
                             ta1_y_yyy_yyyyyz_1,  \
                             ta1_y_yyy_yyyyz_0,   \
                             ta1_y_yyy_yyyyz_1,   \
                             ta1_y_yyy_yyyyzz_0,  \
                             ta1_y_yyy_yyyyzz_1,  \
                             ta1_y_yyy_yyyzz_0,   \
                             ta1_y_yyy_yyyzz_1,   \
                             ta1_y_yyy_yyyzzz_0,  \
                             ta1_y_yyy_yyyzzz_1,  \
                             ta1_y_yyy_yyzzz_0,   \
                             ta1_y_yyy_yyzzz_1,   \
                             ta1_y_yyy_yyzzzz_0,  \
                             ta1_y_yyy_yyzzzz_1,  \
                             ta1_y_yyy_yzzzz_0,   \
                             ta1_y_yyy_yzzzz_1,   \
                             ta1_y_yyy_yzzzzz_0,  \
                             ta1_y_yyy_yzzzzz_1,  \
                             ta1_y_yyy_zzzzz_0,   \
                             ta1_y_yyy_zzzzz_1,   \
                             ta1_y_yyy_zzzzzz_0,  \
                             ta1_y_yyy_zzzzzz_1,  \
                             ta1_y_yyyz_xxxxxx_0, \
                             ta1_y_yyyz_xxxxxy_0, \
                             ta1_y_yyyz_xxxxxz_0, \
                             ta1_y_yyyz_xxxxyy_0, \
                             ta1_y_yyyz_xxxxyz_0, \
                             ta1_y_yyyz_xxxxzz_0, \
                             ta1_y_yyyz_xxxyyy_0, \
                             ta1_y_yyyz_xxxyyz_0, \
                             ta1_y_yyyz_xxxyzz_0, \
                             ta1_y_yyyz_xxxzzz_0, \
                             ta1_y_yyyz_xxyyyy_0, \
                             ta1_y_yyyz_xxyyyz_0, \
                             ta1_y_yyyz_xxyyzz_0, \
                             ta1_y_yyyz_xxyzzz_0, \
                             ta1_y_yyyz_xxzzzz_0, \
                             ta1_y_yyyz_xyyyyy_0, \
                             ta1_y_yyyz_xyyyyz_0, \
                             ta1_y_yyyz_xyyyzz_0, \
                             ta1_y_yyyz_xyyzzz_0, \
                             ta1_y_yyyz_xyzzzz_0, \
                             ta1_y_yyyz_xzzzzz_0, \
                             ta1_y_yyyz_yyyyyy_0, \
                             ta1_y_yyyz_yyyyyz_0, \
                             ta1_y_yyyz_yyyyzz_0, \
                             ta1_y_yyyz_yyyzzz_0, \
                             ta1_y_yyyz_yyzzzz_0, \
                             ta1_y_yyyz_yzzzzz_0, \
                             ta1_y_yyyz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyz_xxxxxx_0[i] = ta1_y_yyy_xxxxxx_0[i] * pa_z[i] - ta1_y_yyy_xxxxxx_1[i] * pc_z[i];

        ta1_y_yyyz_xxxxxy_0[i] = ta1_y_yyy_xxxxxy_0[i] * pa_z[i] - ta1_y_yyy_xxxxxy_1[i] * pc_z[i];

        ta1_y_yyyz_xxxxxz_0[i] =
            ta1_y_yyy_xxxxx_0[i] * fe_0 - ta1_y_yyy_xxxxx_1[i] * fe_0 + ta1_y_yyy_xxxxxz_0[i] * pa_z[i] - ta1_y_yyy_xxxxxz_1[i] * pc_z[i];

        ta1_y_yyyz_xxxxyy_0[i] = ta1_y_yyy_xxxxyy_0[i] * pa_z[i] - ta1_y_yyy_xxxxyy_1[i] * pc_z[i];

        ta1_y_yyyz_xxxxyz_0[i] =
            ta1_y_yyy_xxxxy_0[i] * fe_0 - ta1_y_yyy_xxxxy_1[i] * fe_0 + ta1_y_yyy_xxxxyz_0[i] * pa_z[i] - ta1_y_yyy_xxxxyz_1[i] * pc_z[i];

        ta1_y_yyyz_xxxxzz_0[i] =
            2.0 * ta1_y_yyy_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xxxxz_1[i] * fe_0 + ta1_y_yyy_xxxxzz_0[i] * pa_z[i] - ta1_y_yyy_xxxxzz_1[i] * pc_z[i];

        ta1_y_yyyz_xxxyyy_0[i] = ta1_y_yyy_xxxyyy_0[i] * pa_z[i] - ta1_y_yyy_xxxyyy_1[i] * pc_z[i];

        ta1_y_yyyz_xxxyyz_0[i] =
            ta1_y_yyy_xxxyy_0[i] * fe_0 - ta1_y_yyy_xxxyy_1[i] * fe_0 + ta1_y_yyy_xxxyyz_0[i] * pa_z[i] - ta1_y_yyy_xxxyyz_1[i] * pc_z[i];

        ta1_y_yyyz_xxxyzz_0[i] =
            2.0 * ta1_y_yyy_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xxxyz_1[i] * fe_0 + ta1_y_yyy_xxxyzz_0[i] * pa_z[i] - ta1_y_yyy_xxxyzz_1[i] * pc_z[i];

        ta1_y_yyyz_xxxzzz_0[i] =
            3.0 * ta1_y_yyy_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxxzz_1[i] * fe_0 + ta1_y_yyy_xxxzzz_0[i] * pa_z[i] - ta1_y_yyy_xxxzzz_1[i] * pc_z[i];

        ta1_y_yyyz_xxyyyy_0[i] = ta1_y_yyy_xxyyyy_0[i] * pa_z[i] - ta1_y_yyy_xxyyyy_1[i] * pc_z[i];

        ta1_y_yyyz_xxyyyz_0[i] =
            ta1_y_yyy_xxyyy_0[i] * fe_0 - ta1_y_yyy_xxyyy_1[i] * fe_0 + ta1_y_yyy_xxyyyz_0[i] * pa_z[i] - ta1_y_yyy_xxyyyz_1[i] * pc_z[i];

        ta1_y_yyyz_xxyyzz_0[i] =
            2.0 * ta1_y_yyy_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xxyyz_1[i] * fe_0 + ta1_y_yyy_xxyyzz_0[i] * pa_z[i] - ta1_y_yyy_xxyyzz_1[i] * pc_z[i];

        ta1_y_yyyz_xxyzzz_0[i] =
            3.0 * ta1_y_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxyzz_1[i] * fe_0 + ta1_y_yyy_xxyzzz_0[i] * pa_z[i] - ta1_y_yyy_xxyzzz_1[i] * pc_z[i];

        ta1_y_yyyz_xxzzzz_0[i] =
            4.0 * ta1_y_yyy_xxzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxzzz_1[i] * fe_0 + ta1_y_yyy_xxzzzz_0[i] * pa_z[i] - ta1_y_yyy_xxzzzz_1[i] * pc_z[i];

        ta1_y_yyyz_xyyyyy_0[i] = ta1_y_yyy_xyyyyy_0[i] * pa_z[i] - ta1_y_yyy_xyyyyy_1[i] * pc_z[i];

        ta1_y_yyyz_xyyyyz_0[i] =
            ta1_y_yyy_xyyyy_0[i] * fe_0 - ta1_y_yyy_xyyyy_1[i] * fe_0 + ta1_y_yyy_xyyyyz_0[i] * pa_z[i] - ta1_y_yyy_xyyyyz_1[i] * pc_z[i];

        ta1_y_yyyz_xyyyzz_0[i] =
            2.0 * ta1_y_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyyyz_1[i] * fe_0 + ta1_y_yyy_xyyyzz_0[i] * pa_z[i] - ta1_y_yyy_xyyyzz_1[i] * pc_z[i];

        ta1_y_yyyz_xyyzzz_0[i] =
            3.0 * ta1_y_yyy_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xyyzz_1[i] * fe_0 + ta1_y_yyy_xyyzzz_0[i] * pa_z[i] - ta1_y_yyy_xyyzzz_1[i] * pc_z[i];

        ta1_y_yyyz_xyzzzz_0[i] =
            4.0 * ta1_y_yyy_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xyzzz_1[i] * fe_0 + ta1_y_yyy_xyzzzz_0[i] * pa_z[i] - ta1_y_yyy_xyzzzz_1[i] * pc_z[i];

        ta1_y_yyyz_xzzzzz_0[i] =
            5.0 * ta1_y_yyy_xzzzz_0[i] * fe_0 - 5.0 * ta1_y_yyy_xzzzz_1[i] * fe_0 + ta1_y_yyy_xzzzzz_0[i] * pa_z[i] - ta1_y_yyy_xzzzzz_1[i] * pc_z[i];

        ta1_y_yyyz_yyyyyy_0[i] = ta1_y_yyy_yyyyyy_0[i] * pa_z[i] - ta1_y_yyy_yyyyyy_1[i] * pc_z[i];

        ta1_y_yyyz_yyyyyz_0[i] =
            ta1_y_yyy_yyyyy_0[i] * fe_0 - ta1_y_yyy_yyyyy_1[i] * fe_0 + ta1_y_yyy_yyyyyz_0[i] * pa_z[i] - ta1_y_yyy_yyyyyz_1[i] * pc_z[i];

        ta1_y_yyyz_yyyyzz_0[i] =
            2.0 * ta1_y_yyy_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_yyyyz_1[i] * fe_0 + ta1_y_yyy_yyyyzz_0[i] * pa_z[i] - ta1_y_yyy_yyyyzz_1[i] * pc_z[i];

        ta1_y_yyyz_yyyzzz_0[i] =
            3.0 * ta1_y_yyy_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_yyyzz_1[i] * fe_0 + ta1_y_yyy_yyyzzz_0[i] * pa_z[i] - ta1_y_yyy_yyyzzz_1[i] * pc_z[i];

        ta1_y_yyyz_yyzzzz_0[i] =
            4.0 * ta1_y_yyy_yyzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yyzzz_1[i] * fe_0 + ta1_y_yyy_yyzzzz_0[i] * pa_z[i] - ta1_y_yyy_yyzzzz_1[i] * pc_z[i];

        ta1_y_yyyz_yzzzzz_0[i] =
            5.0 * ta1_y_yyy_yzzzz_0[i] * fe_0 - 5.0 * ta1_y_yyy_yzzzz_1[i] * fe_0 + ta1_y_yyy_yzzzzz_0[i] * pa_z[i] - ta1_y_yyy_yzzzzz_1[i] * pc_z[i];

        ta1_y_yyyz_zzzzzz_0[i] =
            6.0 * ta1_y_yyy_zzzzz_0[i] * fe_0 - 6.0 * ta1_y_yyy_zzzzz_1[i] * fe_0 + ta1_y_yyy_zzzzzz_0[i] * pa_z[i] - ta1_y_yyy_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 756-784 components of targeted buffer : GI

    auto ta1_y_yyzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 756);

    auto ta1_y_yyzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 757);

    auto ta1_y_yyzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 758);

    auto ta1_y_yyzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 759);

    auto ta1_y_yyzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 760);

    auto ta1_y_yyzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 761);

    auto ta1_y_yyzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 762);

    auto ta1_y_yyzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 763);

    auto ta1_y_yyzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 764);

    auto ta1_y_yyzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 765);

    auto ta1_y_yyzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 766);

    auto ta1_y_yyzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 767);

    auto ta1_y_yyzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 768);

    auto ta1_y_yyzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 769);

    auto ta1_y_yyzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 770);

    auto ta1_y_yyzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 771);

    auto ta1_y_yyzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 772);

    auto ta1_y_yyzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 773);

    auto ta1_y_yyzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 774);

    auto ta1_y_yyzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 775);

    auto ta1_y_yyzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 776);

    auto ta1_y_yyzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 777);

    auto ta1_y_yyzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 778);

    auto ta1_y_yyzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 779);

    auto ta1_y_yyzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 780);

    auto ta1_y_yyzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 781);

    auto ta1_y_yyzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 782);

    auto ta1_y_yyzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 783);

#pragma omp simd aligned(pa_y,                    \
                             pa_z,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_y_yy_xxxxxx_0,   \
                             ta1_y_yy_xxxxxx_1,   \
                             ta1_y_yy_xxxxxy_0,   \
                             ta1_y_yy_xxxxxy_1,   \
                             ta1_y_yy_xxxxyy_0,   \
                             ta1_y_yy_xxxxyy_1,   \
                             ta1_y_yy_xxxxyz_0,   \
                             ta1_y_yy_xxxxyz_1,   \
                             ta1_y_yy_xxxyyy_0,   \
                             ta1_y_yy_xxxyyy_1,   \
                             ta1_y_yy_xxxyyz_0,   \
                             ta1_y_yy_xxxyyz_1,   \
                             ta1_y_yy_xxxyzz_0,   \
                             ta1_y_yy_xxxyzz_1,   \
                             ta1_y_yy_xxyyyy_0,   \
                             ta1_y_yy_xxyyyy_1,   \
                             ta1_y_yy_xxyyyz_0,   \
                             ta1_y_yy_xxyyyz_1,   \
                             ta1_y_yy_xxyyzz_0,   \
                             ta1_y_yy_xxyyzz_1,   \
                             ta1_y_yy_xxyzzz_0,   \
                             ta1_y_yy_xxyzzz_1,   \
                             ta1_y_yy_xyyyyy_0,   \
                             ta1_y_yy_xyyyyy_1,   \
                             ta1_y_yy_xyyyyz_0,   \
                             ta1_y_yy_xyyyyz_1,   \
                             ta1_y_yy_xyyyzz_0,   \
                             ta1_y_yy_xyyyzz_1,   \
                             ta1_y_yy_xyyzzz_0,   \
                             ta1_y_yy_xyyzzz_1,   \
                             ta1_y_yy_xyzzzz_0,   \
                             ta1_y_yy_xyzzzz_1,   \
                             ta1_y_yy_yyyyyy_0,   \
                             ta1_y_yy_yyyyyy_1,   \
                             ta1_y_yy_yyyyyz_0,   \
                             ta1_y_yy_yyyyyz_1,   \
                             ta1_y_yy_yyyyzz_0,   \
                             ta1_y_yy_yyyyzz_1,   \
                             ta1_y_yy_yyyzzz_0,   \
                             ta1_y_yy_yyyzzz_1,   \
                             ta1_y_yy_yyzzzz_0,   \
                             ta1_y_yy_yyzzzz_1,   \
                             ta1_y_yy_yzzzzz_0,   \
                             ta1_y_yy_yzzzzz_1,   \
                             ta1_y_yyz_xxxxxx_0,  \
                             ta1_y_yyz_xxxxxx_1,  \
                             ta1_y_yyz_xxxxxy_0,  \
                             ta1_y_yyz_xxxxxy_1,  \
                             ta1_y_yyz_xxxxy_0,   \
                             ta1_y_yyz_xxxxy_1,   \
                             ta1_y_yyz_xxxxyy_0,  \
                             ta1_y_yyz_xxxxyy_1,  \
                             ta1_y_yyz_xxxxyz_0,  \
                             ta1_y_yyz_xxxxyz_1,  \
                             ta1_y_yyz_xxxyy_0,   \
                             ta1_y_yyz_xxxyy_1,   \
                             ta1_y_yyz_xxxyyy_0,  \
                             ta1_y_yyz_xxxyyy_1,  \
                             ta1_y_yyz_xxxyyz_0,  \
                             ta1_y_yyz_xxxyyz_1,  \
                             ta1_y_yyz_xxxyz_0,   \
                             ta1_y_yyz_xxxyz_1,   \
                             ta1_y_yyz_xxxyzz_0,  \
                             ta1_y_yyz_xxxyzz_1,  \
                             ta1_y_yyz_xxyyy_0,   \
                             ta1_y_yyz_xxyyy_1,   \
                             ta1_y_yyz_xxyyyy_0,  \
                             ta1_y_yyz_xxyyyy_1,  \
                             ta1_y_yyz_xxyyyz_0,  \
                             ta1_y_yyz_xxyyyz_1,  \
                             ta1_y_yyz_xxyyz_0,   \
                             ta1_y_yyz_xxyyz_1,   \
                             ta1_y_yyz_xxyyzz_0,  \
                             ta1_y_yyz_xxyyzz_1,  \
                             ta1_y_yyz_xxyzz_0,   \
                             ta1_y_yyz_xxyzz_1,   \
                             ta1_y_yyz_xxyzzz_0,  \
                             ta1_y_yyz_xxyzzz_1,  \
                             ta1_y_yyz_xyyyy_0,   \
                             ta1_y_yyz_xyyyy_1,   \
                             ta1_y_yyz_xyyyyy_0,  \
                             ta1_y_yyz_xyyyyy_1,  \
                             ta1_y_yyz_xyyyyz_0,  \
                             ta1_y_yyz_xyyyyz_1,  \
                             ta1_y_yyz_xyyyz_0,   \
                             ta1_y_yyz_xyyyz_1,   \
                             ta1_y_yyz_xyyyzz_0,  \
                             ta1_y_yyz_xyyyzz_1,  \
                             ta1_y_yyz_xyyzz_0,   \
                             ta1_y_yyz_xyyzz_1,   \
                             ta1_y_yyz_xyyzzz_0,  \
                             ta1_y_yyz_xyyzzz_1,  \
                             ta1_y_yyz_xyzzz_0,   \
                             ta1_y_yyz_xyzzz_1,   \
                             ta1_y_yyz_xyzzzz_0,  \
                             ta1_y_yyz_xyzzzz_1,  \
                             ta1_y_yyz_yyyyy_0,   \
                             ta1_y_yyz_yyyyy_1,   \
                             ta1_y_yyz_yyyyyy_0,  \
                             ta1_y_yyz_yyyyyy_1,  \
                             ta1_y_yyz_yyyyyz_0,  \
                             ta1_y_yyz_yyyyyz_1,  \
                             ta1_y_yyz_yyyyz_0,   \
                             ta1_y_yyz_yyyyz_1,   \
                             ta1_y_yyz_yyyyzz_0,  \
                             ta1_y_yyz_yyyyzz_1,  \
                             ta1_y_yyz_yyyzz_0,   \
                             ta1_y_yyz_yyyzz_1,   \
                             ta1_y_yyz_yyyzzz_0,  \
                             ta1_y_yyz_yyyzzz_1,  \
                             ta1_y_yyz_yyzzz_0,   \
                             ta1_y_yyz_yyzzz_1,   \
                             ta1_y_yyz_yyzzzz_0,  \
                             ta1_y_yyz_yyzzzz_1,  \
                             ta1_y_yyz_yzzzz_0,   \
                             ta1_y_yyz_yzzzz_1,   \
                             ta1_y_yyz_yzzzzz_0,  \
                             ta1_y_yyz_yzzzzz_1,  \
                             ta1_y_yyzz_xxxxxx_0, \
                             ta1_y_yyzz_xxxxxy_0, \
                             ta1_y_yyzz_xxxxxz_0, \
                             ta1_y_yyzz_xxxxyy_0, \
                             ta1_y_yyzz_xxxxyz_0, \
                             ta1_y_yyzz_xxxxzz_0, \
                             ta1_y_yyzz_xxxyyy_0, \
                             ta1_y_yyzz_xxxyyz_0, \
                             ta1_y_yyzz_xxxyzz_0, \
                             ta1_y_yyzz_xxxzzz_0, \
                             ta1_y_yyzz_xxyyyy_0, \
                             ta1_y_yyzz_xxyyyz_0, \
                             ta1_y_yyzz_xxyyzz_0, \
                             ta1_y_yyzz_xxyzzz_0, \
                             ta1_y_yyzz_xxzzzz_0, \
                             ta1_y_yyzz_xyyyyy_0, \
                             ta1_y_yyzz_xyyyyz_0, \
                             ta1_y_yyzz_xyyyzz_0, \
                             ta1_y_yyzz_xyyzzz_0, \
                             ta1_y_yyzz_xyzzzz_0, \
                             ta1_y_yyzz_xzzzzz_0, \
                             ta1_y_yyzz_yyyyyy_0, \
                             ta1_y_yyzz_yyyyyz_0, \
                             ta1_y_yyzz_yyyyzz_0, \
                             ta1_y_yyzz_yyyzzz_0, \
                             ta1_y_yyzz_yyzzzz_0, \
                             ta1_y_yyzz_yzzzzz_0, \
                             ta1_y_yyzz_zzzzzz_0, \
                             ta1_y_yzz_xxxxxz_0,  \
                             ta1_y_yzz_xxxxxz_1,  \
                             ta1_y_yzz_xxxxzz_0,  \
                             ta1_y_yzz_xxxxzz_1,  \
                             ta1_y_yzz_xxxzzz_0,  \
                             ta1_y_yzz_xxxzzz_1,  \
                             ta1_y_yzz_xxzzzz_0,  \
                             ta1_y_yzz_xxzzzz_1,  \
                             ta1_y_yzz_xzzzzz_0,  \
                             ta1_y_yzz_xzzzzz_1,  \
                             ta1_y_yzz_zzzzzz_0,  \
                             ta1_y_yzz_zzzzzz_1,  \
                             ta1_y_zz_xxxxxz_0,   \
                             ta1_y_zz_xxxxxz_1,   \
                             ta1_y_zz_xxxxzz_0,   \
                             ta1_y_zz_xxxxzz_1,   \
                             ta1_y_zz_xxxzzz_0,   \
                             ta1_y_zz_xxxzzz_1,   \
                             ta1_y_zz_xxzzzz_0,   \
                             ta1_y_zz_xxzzzz_1,   \
                             ta1_y_zz_xzzzzz_0,   \
                             ta1_y_zz_xzzzzz_1,   \
                             ta1_y_zz_zzzzzz_0,   \
                             ta1_y_zz_zzzzzz_1,   \
                             ta_yzz_xxxxxz_1,     \
                             ta_yzz_xxxxzz_1,     \
                             ta_yzz_xxxzzz_1,     \
                             ta_yzz_xxzzzz_1,     \
                             ta_yzz_xzzzzz_1,     \
                             ta_yzz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzz_xxxxxx_0[i] =
            ta1_y_yy_xxxxxx_0[i] * fe_0 - ta1_y_yy_xxxxxx_1[i] * fe_0 + ta1_y_yyz_xxxxxx_0[i] * pa_z[i] - ta1_y_yyz_xxxxxx_1[i] * pc_z[i];

        ta1_y_yyzz_xxxxxy_0[i] =
            ta1_y_yy_xxxxxy_0[i] * fe_0 - ta1_y_yy_xxxxxy_1[i] * fe_0 + ta1_y_yyz_xxxxxy_0[i] * pa_z[i] - ta1_y_yyz_xxxxxy_1[i] * pc_z[i];

        ta1_y_yyzz_xxxxxz_0[i] = ta1_y_zz_xxxxxz_0[i] * fe_0 - ta1_y_zz_xxxxxz_1[i] * fe_0 + ta_yzz_xxxxxz_1[i] + ta1_y_yzz_xxxxxz_0[i] * pa_y[i] -
                                 ta1_y_yzz_xxxxxz_1[i] * pc_y[i];

        ta1_y_yyzz_xxxxyy_0[i] =
            ta1_y_yy_xxxxyy_0[i] * fe_0 - ta1_y_yy_xxxxyy_1[i] * fe_0 + ta1_y_yyz_xxxxyy_0[i] * pa_z[i] - ta1_y_yyz_xxxxyy_1[i] * pc_z[i];

        ta1_y_yyzz_xxxxyz_0[i] = ta1_y_yy_xxxxyz_0[i] * fe_0 - ta1_y_yy_xxxxyz_1[i] * fe_0 + ta1_y_yyz_xxxxy_0[i] * fe_0 -
                                 ta1_y_yyz_xxxxy_1[i] * fe_0 + ta1_y_yyz_xxxxyz_0[i] * pa_z[i] - ta1_y_yyz_xxxxyz_1[i] * pc_z[i];

        ta1_y_yyzz_xxxxzz_0[i] = ta1_y_zz_xxxxzz_0[i] * fe_0 - ta1_y_zz_xxxxzz_1[i] * fe_0 + ta_yzz_xxxxzz_1[i] + ta1_y_yzz_xxxxzz_0[i] * pa_y[i] -
                                 ta1_y_yzz_xxxxzz_1[i] * pc_y[i];

        ta1_y_yyzz_xxxyyy_0[i] =
            ta1_y_yy_xxxyyy_0[i] * fe_0 - ta1_y_yy_xxxyyy_1[i] * fe_0 + ta1_y_yyz_xxxyyy_0[i] * pa_z[i] - ta1_y_yyz_xxxyyy_1[i] * pc_z[i];

        ta1_y_yyzz_xxxyyz_0[i] = ta1_y_yy_xxxyyz_0[i] * fe_0 - ta1_y_yy_xxxyyz_1[i] * fe_0 + ta1_y_yyz_xxxyy_0[i] * fe_0 -
                                 ta1_y_yyz_xxxyy_1[i] * fe_0 + ta1_y_yyz_xxxyyz_0[i] * pa_z[i] - ta1_y_yyz_xxxyyz_1[i] * pc_z[i];

        ta1_y_yyzz_xxxyzz_0[i] = ta1_y_yy_xxxyzz_0[i] * fe_0 - ta1_y_yy_xxxyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyz_xxxyz_1[i] * fe_0 + ta1_y_yyz_xxxyzz_0[i] * pa_z[i] - ta1_y_yyz_xxxyzz_1[i] * pc_z[i];

        ta1_y_yyzz_xxxzzz_0[i] = ta1_y_zz_xxxzzz_0[i] * fe_0 - ta1_y_zz_xxxzzz_1[i] * fe_0 + ta_yzz_xxxzzz_1[i] + ta1_y_yzz_xxxzzz_0[i] * pa_y[i] -
                                 ta1_y_yzz_xxxzzz_1[i] * pc_y[i];

        ta1_y_yyzz_xxyyyy_0[i] =
            ta1_y_yy_xxyyyy_0[i] * fe_0 - ta1_y_yy_xxyyyy_1[i] * fe_0 + ta1_y_yyz_xxyyyy_0[i] * pa_z[i] - ta1_y_yyz_xxyyyy_1[i] * pc_z[i];

        ta1_y_yyzz_xxyyyz_0[i] = ta1_y_yy_xxyyyz_0[i] * fe_0 - ta1_y_yy_xxyyyz_1[i] * fe_0 + ta1_y_yyz_xxyyy_0[i] * fe_0 -
                                 ta1_y_yyz_xxyyy_1[i] * fe_0 + ta1_y_yyz_xxyyyz_0[i] * pa_z[i] - ta1_y_yyz_xxyyyz_1[i] * pc_z[i];

        ta1_y_yyzz_xxyyzz_0[i] = ta1_y_yy_xxyyzz_0[i] * fe_0 - ta1_y_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_xxyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyz_xxyyz_1[i] * fe_0 + ta1_y_yyz_xxyyzz_0[i] * pa_z[i] - ta1_y_yyz_xxyyzz_1[i] * pc_z[i];

        ta1_y_yyzz_xxyzzz_0[i] = ta1_y_yy_xxyzzz_0[i] * fe_0 - ta1_y_yy_xxyzzz_1[i] * fe_0 + 3.0 * ta1_y_yyz_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_yyz_xxyzz_1[i] * fe_0 + ta1_y_yyz_xxyzzz_0[i] * pa_z[i] - ta1_y_yyz_xxyzzz_1[i] * pc_z[i];

        ta1_y_yyzz_xxzzzz_0[i] = ta1_y_zz_xxzzzz_0[i] * fe_0 - ta1_y_zz_xxzzzz_1[i] * fe_0 + ta_yzz_xxzzzz_1[i] + ta1_y_yzz_xxzzzz_0[i] * pa_y[i] -
                                 ta1_y_yzz_xxzzzz_1[i] * pc_y[i];

        ta1_y_yyzz_xyyyyy_0[i] =
            ta1_y_yy_xyyyyy_0[i] * fe_0 - ta1_y_yy_xyyyyy_1[i] * fe_0 + ta1_y_yyz_xyyyyy_0[i] * pa_z[i] - ta1_y_yyz_xyyyyy_1[i] * pc_z[i];

        ta1_y_yyzz_xyyyyz_0[i] = ta1_y_yy_xyyyyz_0[i] * fe_0 - ta1_y_yy_xyyyyz_1[i] * fe_0 + ta1_y_yyz_xyyyy_0[i] * fe_0 -
                                 ta1_y_yyz_xyyyy_1[i] * fe_0 + ta1_y_yyz_xyyyyz_0[i] * pa_z[i] - ta1_y_yyz_xyyyyz_1[i] * pc_z[i];

        ta1_y_yyzz_xyyyzz_0[i] = ta1_y_yy_xyyyzz_0[i] * fe_0 - ta1_y_yy_xyyyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyz_xyyyz_1[i] * fe_0 + ta1_y_yyz_xyyyzz_0[i] * pa_z[i] - ta1_y_yyz_xyyyzz_1[i] * pc_z[i];

        ta1_y_yyzz_xyyzzz_0[i] = ta1_y_yy_xyyzzz_0[i] * fe_0 - ta1_y_yy_xyyzzz_1[i] * fe_0 + 3.0 * ta1_y_yyz_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_yyz_xyyzz_1[i] * fe_0 + ta1_y_yyz_xyyzzz_0[i] * pa_z[i] - ta1_y_yyz_xyyzzz_1[i] * pc_z[i];

        ta1_y_yyzz_xyzzzz_0[i] = ta1_y_yy_xyzzzz_0[i] * fe_0 - ta1_y_yy_xyzzzz_1[i] * fe_0 + 4.0 * ta1_y_yyz_xyzzz_0[i] * fe_0 -
                                 4.0 * ta1_y_yyz_xyzzz_1[i] * fe_0 + ta1_y_yyz_xyzzzz_0[i] * pa_z[i] - ta1_y_yyz_xyzzzz_1[i] * pc_z[i];

        ta1_y_yyzz_xzzzzz_0[i] = ta1_y_zz_xzzzzz_0[i] * fe_0 - ta1_y_zz_xzzzzz_1[i] * fe_0 + ta_yzz_xzzzzz_1[i] + ta1_y_yzz_xzzzzz_0[i] * pa_y[i] -
                                 ta1_y_yzz_xzzzzz_1[i] * pc_y[i];

        ta1_y_yyzz_yyyyyy_0[i] =
            ta1_y_yy_yyyyyy_0[i] * fe_0 - ta1_y_yy_yyyyyy_1[i] * fe_0 + ta1_y_yyz_yyyyyy_0[i] * pa_z[i] - ta1_y_yyz_yyyyyy_1[i] * pc_z[i];

        ta1_y_yyzz_yyyyyz_0[i] = ta1_y_yy_yyyyyz_0[i] * fe_0 - ta1_y_yy_yyyyyz_1[i] * fe_0 + ta1_y_yyz_yyyyy_0[i] * fe_0 -
                                 ta1_y_yyz_yyyyy_1[i] * fe_0 + ta1_y_yyz_yyyyyz_0[i] * pa_z[i] - ta1_y_yyz_yyyyyz_1[i] * pc_z[i];

        ta1_y_yyzz_yyyyzz_0[i] = ta1_y_yy_yyyyzz_0[i] * fe_0 - ta1_y_yy_yyyyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_yyyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_yyz_yyyyz_1[i] * fe_0 + ta1_y_yyz_yyyyzz_0[i] * pa_z[i] - ta1_y_yyz_yyyyzz_1[i] * pc_z[i];

        ta1_y_yyzz_yyyzzz_0[i] = ta1_y_yy_yyyzzz_0[i] * fe_0 - ta1_y_yy_yyyzzz_1[i] * fe_0 + 3.0 * ta1_y_yyz_yyyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_yyz_yyyzz_1[i] * fe_0 + ta1_y_yyz_yyyzzz_0[i] * pa_z[i] - ta1_y_yyz_yyyzzz_1[i] * pc_z[i];

        ta1_y_yyzz_yyzzzz_0[i] = ta1_y_yy_yyzzzz_0[i] * fe_0 - ta1_y_yy_yyzzzz_1[i] * fe_0 + 4.0 * ta1_y_yyz_yyzzz_0[i] * fe_0 -
                                 4.0 * ta1_y_yyz_yyzzz_1[i] * fe_0 + ta1_y_yyz_yyzzzz_0[i] * pa_z[i] - ta1_y_yyz_yyzzzz_1[i] * pc_z[i];

        ta1_y_yyzz_yzzzzz_0[i] = ta1_y_yy_yzzzzz_0[i] * fe_0 - ta1_y_yy_yzzzzz_1[i] * fe_0 + 5.0 * ta1_y_yyz_yzzzz_0[i] * fe_0 -
                                 5.0 * ta1_y_yyz_yzzzz_1[i] * fe_0 + ta1_y_yyz_yzzzzz_0[i] * pa_z[i] - ta1_y_yyz_yzzzzz_1[i] * pc_z[i];

        ta1_y_yyzz_zzzzzz_0[i] = ta1_y_zz_zzzzzz_0[i] * fe_0 - ta1_y_zz_zzzzzz_1[i] * fe_0 + ta_yzz_zzzzzz_1[i] + ta1_y_yzz_zzzzzz_0[i] * pa_y[i] -
                                 ta1_y_yzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 784-812 components of targeted buffer : GI

    auto ta1_y_yzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 784);

    auto ta1_y_yzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 785);

    auto ta1_y_yzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 786);

    auto ta1_y_yzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 787);

    auto ta1_y_yzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 788);

    auto ta1_y_yzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 789);

    auto ta1_y_yzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 790);

    auto ta1_y_yzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 791);

    auto ta1_y_yzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 792);

    auto ta1_y_yzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 793);

    auto ta1_y_yzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 794);

    auto ta1_y_yzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 795);

    auto ta1_y_yzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 796);

    auto ta1_y_yzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 797);

    auto ta1_y_yzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 798);

    auto ta1_y_yzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 799);

    auto ta1_y_yzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 800);

    auto ta1_y_yzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 801);

    auto ta1_y_yzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 802);

    auto ta1_y_yzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 803);

    auto ta1_y_yzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 804);

    auto ta1_y_yzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 805);

    auto ta1_y_yzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 806);

    auto ta1_y_yzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 807);

    auto ta1_y_yzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 808);

    auto ta1_y_yzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 809);

    auto ta1_y_yzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 810);

    auto ta1_y_yzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 811);

#pragma omp simd aligned(pa_y,                    \
                             pa_z,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_y_yz_xxxxxy_0,   \
                             ta1_y_yz_xxxxxy_1,   \
                             ta1_y_yz_xxxxyy_0,   \
                             ta1_y_yz_xxxxyy_1,   \
                             ta1_y_yz_xxxyyy_0,   \
                             ta1_y_yz_xxxyyy_1,   \
                             ta1_y_yz_xxyyyy_0,   \
                             ta1_y_yz_xxyyyy_1,   \
                             ta1_y_yz_xyyyyy_0,   \
                             ta1_y_yz_xyyyyy_1,   \
                             ta1_y_yz_yyyyyy_0,   \
                             ta1_y_yz_yyyyyy_1,   \
                             ta1_y_yzz_xxxxxy_0,  \
                             ta1_y_yzz_xxxxxy_1,  \
                             ta1_y_yzz_xxxxyy_0,  \
                             ta1_y_yzz_xxxxyy_1,  \
                             ta1_y_yzz_xxxyyy_0,  \
                             ta1_y_yzz_xxxyyy_1,  \
                             ta1_y_yzz_xxyyyy_0,  \
                             ta1_y_yzz_xxyyyy_1,  \
                             ta1_y_yzz_xyyyyy_0,  \
                             ta1_y_yzz_xyyyyy_1,  \
                             ta1_y_yzz_yyyyyy_0,  \
                             ta1_y_yzz_yyyyyy_1,  \
                             ta1_y_yzzz_xxxxxx_0, \
                             ta1_y_yzzz_xxxxxy_0, \
                             ta1_y_yzzz_xxxxxz_0, \
                             ta1_y_yzzz_xxxxyy_0, \
                             ta1_y_yzzz_xxxxyz_0, \
                             ta1_y_yzzz_xxxxzz_0, \
                             ta1_y_yzzz_xxxyyy_0, \
                             ta1_y_yzzz_xxxyyz_0, \
                             ta1_y_yzzz_xxxyzz_0, \
                             ta1_y_yzzz_xxxzzz_0, \
                             ta1_y_yzzz_xxyyyy_0, \
                             ta1_y_yzzz_xxyyyz_0, \
                             ta1_y_yzzz_xxyyzz_0, \
                             ta1_y_yzzz_xxyzzz_0, \
                             ta1_y_yzzz_xxzzzz_0, \
                             ta1_y_yzzz_xyyyyy_0, \
                             ta1_y_yzzz_xyyyyz_0, \
                             ta1_y_yzzz_xyyyzz_0, \
                             ta1_y_yzzz_xyyzzz_0, \
                             ta1_y_yzzz_xyzzzz_0, \
                             ta1_y_yzzz_xzzzzz_0, \
                             ta1_y_yzzz_yyyyyy_0, \
                             ta1_y_yzzz_yyyyyz_0, \
                             ta1_y_yzzz_yyyyzz_0, \
                             ta1_y_yzzz_yyyzzz_0, \
                             ta1_y_yzzz_yyzzzz_0, \
                             ta1_y_yzzz_yzzzzz_0, \
                             ta1_y_yzzz_zzzzzz_0, \
                             ta1_y_zzz_xxxxxx_0,  \
                             ta1_y_zzz_xxxxxx_1,  \
                             ta1_y_zzz_xxxxxz_0,  \
                             ta1_y_zzz_xxxxxz_1,  \
                             ta1_y_zzz_xxxxyz_0,  \
                             ta1_y_zzz_xxxxyz_1,  \
                             ta1_y_zzz_xxxxz_0,   \
                             ta1_y_zzz_xxxxz_1,   \
                             ta1_y_zzz_xxxxzz_0,  \
                             ta1_y_zzz_xxxxzz_1,  \
                             ta1_y_zzz_xxxyyz_0,  \
                             ta1_y_zzz_xxxyyz_1,  \
                             ta1_y_zzz_xxxyz_0,   \
                             ta1_y_zzz_xxxyz_1,   \
                             ta1_y_zzz_xxxyzz_0,  \
                             ta1_y_zzz_xxxyzz_1,  \
                             ta1_y_zzz_xxxzz_0,   \
                             ta1_y_zzz_xxxzz_1,   \
                             ta1_y_zzz_xxxzzz_0,  \
                             ta1_y_zzz_xxxzzz_1,  \
                             ta1_y_zzz_xxyyyz_0,  \
                             ta1_y_zzz_xxyyyz_1,  \
                             ta1_y_zzz_xxyyz_0,   \
                             ta1_y_zzz_xxyyz_1,   \
                             ta1_y_zzz_xxyyzz_0,  \
                             ta1_y_zzz_xxyyzz_1,  \
                             ta1_y_zzz_xxyzz_0,   \
                             ta1_y_zzz_xxyzz_1,   \
                             ta1_y_zzz_xxyzzz_0,  \
                             ta1_y_zzz_xxyzzz_1,  \
                             ta1_y_zzz_xxzzz_0,   \
                             ta1_y_zzz_xxzzz_1,   \
                             ta1_y_zzz_xxzzzz_0,  \
                             ta1_y_zzz_xxzzzz_1,  \
                             ta1_y_zzz_xyyyyz_0,  \
                             ta1_y_zzz_xyyyyz_1,  \
                             ta1_y_zzz_xyyyz_0,   \
                             ta1_y_zzz_xyyyz_1,   \
                             ta1_y_zzz_xyyyzz_0,  \
                             ta1_y_zzz_xyyyzz_1,  \
                             ta1_y_zzz_xyyzz_0,   \
                             ta1_y_zzz_xyyzz_1,   \
                             ta1_y_zzz_xyyzzz_0,  \
                             ta1_y_zzz_xyyzzz_1,  \
                             ta1_y_zzz_xyzzz_0,   \
                             ta1_y_zzz_xyzzz_1,   \
                             ta1_y_zzz_xyzzzz_0,  \
                             ta1_y_zzz_xyzzzz_1,  \
                             ta1_y_zzz_xzzzz_0,   \
                             ta1_y_zzz_xzzzz_1,   \
                             ta1_y_zzz_xzzzzz_0,  \
                             ta1_y_zzz_xzzzzz_1,  \
                             ta1_y_zzz_yyyyyz_0,  \
                             ta1_y_zzz_yyyyyz_1,  \
                             ta1_y_zzz_yyyyz_0,   \
                             ta1_y_zzz_yyyyz_1,   \
                             ta1_y_zzz_yyyyzz_0,  \
                             ta1_y_zzz_yyyyzz_1,  \
                             ta1_y_zzz_yyyzz_0,   \
                             ta1_y_zzz_yyyzz_1,   \
                             ta1_y_zzz_yyyzzz_0,  \
                             ta1_y_zzz_yyyzzz_1,  \
                             ta1_y_zzz_yyzzz_0,   \
                             ta1_y_zzz_yyzzz_1,   \
                             ta1_y_zzz_yyzzzz_0,  \
                             ta1_y_zzz_yyzzzz_1,  \
                             ta1_y_zzz_yzzzz_0,   \
                             ta1_y_zzz_yzzzz_1,   \
                             ta1_y_zzz_yzzzzz_0,  \
                             ta1_y_zzz_yzzzzz_1,  \
                             ta1_y_zzz_zzzzz_0,   \
                             ta1_y_zzz_zzzzz_1,   \
                             ta1_y_zzz_zzzzzz_0,  \
                             ta1_y_zzz_zzzzzz_1,  \
                             ta_zzz_xxxxxx_1,     \
                             ta_zzz_xxxxxz_1,     \
                             ta_zzz_xxxxyz_1,     \
                             ta_zzz_xxxxzz_1,     \
                             ta_zzz_xxxyyz_1,     \
                             ta_zzz_xxxyzz_1,     \
                             ta_zzz_xxxzzz_1,     \
                             ta_zzz_xxyyyz_1,     \
                             ta_zzz_xxyyzz_1,     \
                             ta_zzz_xxyzzz_1,     \
                             ta_zzz_xxzzzz_1,     \
                             ta_zzz_xyyyyz_1,     \
                             ta_zzz_xyyyzz_1,     \
                             ta_zzz_xyyzzz_1,     \
                             ta_zzz_xyzzzz_1,     \
                             ta_zzz_xzzzzz_1,     \
                             ta_zzz_yyyyyz_1,     \
                             ta_zzz_yyyyzz_1,     \
                             ta_zzz_yyyzzz_1,     \
                             ta_zzz_yyzzzz_1,     \
                             ta_zzz_yzzzzz_1,     \
                             ta_zzz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzz_xxxxxx_0[i] = ta_zzz_xxxxxx_1[i] + ta1_y_zzz_xxxxxx_0[i] * pa_y[i] - ta1_y_zzz_xxxxxx_1[i] * pc_y[i];

        ta1_y_yzzz_xxxxxy_0[i] =
            2.0 * ta1_y_yz_xxxxxy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxxxxy_1[i] * fe_0 + ta1_y_yzz_xxxxxy_0[i] * pa_z[i] - ta1_y_yzz_xxxxxy_1[i] * pc_z[i];

        ta1_y_yzzz_xxxxxz_0[i] = ta_zzz_xxxxxz_1[i] + ta1_y_zzz_xxxxxz_0[i] * pa_y[i] - ta1_y_zzz_xxxxxz_1[i] * pc_y[i];

        ta1_y_yzzz_xxxxyy_0[i] =
            2.0 * ta1_y_yz_xxxxyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxxxyy_1[i] * fe_0 + ta1_y_yzz_xxxxyy_0[i] * pa_z[i] - ta1_y_yzz_xxxxyy_1[i] * pc_z[i];

        ta1_y_yzzz_xxxxyz_0[i] = ta1_y_zzz_xxxxz_0[i] * fe_0 - ta1_y_zzz_xxxxz_1[i] * fe_0 + ta_zzz_xxxxyz_1[i] + ta1_y_zzz_xxxxyz_0[i] * pa_y[i] -
                                 ta1_y_zzz_xxxxyz_1[i] * pc_y[i];

        ta1_y_yzzz_xxxxzz_0[i] = ta_zzz_xxxxzz_1[i] + ta1_y_zzz_xxxxzz_0[i] * pa_y[i] - ta1_y_zzz_xxxxzz_1[i] * pc_y[i];

        ta1_y_yzzz_xxxyyy_0[i] =
            2.0 * ta1_y_yz_xxxyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxxyyy_1[i] * fe_0 + ta1_y_yzz_xxxyyy_0[i] * pa_z[i] - ta1_y_yzz_xxxyyy_1[i] * pc_z[i];

        ta1_y_yzzz_xxxyyz_0[i] = 2.0 * ta1_y_zzz_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xxxyz_1[i] * fe_0 + ta_zzz_xxxyyz_1[i] +
                                 ta1_y_zzz_xxxyyz_0[i] * pa_y[i] - ta1_y_zzz_xxxyyz_1[i] * pc_y[i];

        ta1_y_yzzz_xxxyzz_0[i] = ta1_y_zzz_xxxzz_0[i] * fe_0 - ta1_y_zzz_xxxzz_1[i] * fe_0 + ta_zzz_xxxyzz_1[i] + ta1_y_zzz_xxxyzz_0[i] * pa_y[i] -
                                 ta1_y_zzz_xxxyzz_1[i] * pc_y[i];

        ta1_y_yzzz_xxxzzz_0[i] = ta_zzz_xxxzzz_1[i] + ta1_y_zzz_xxxzzz_0[i] * pa_y[i] - ta1_y_zzz_xxxzzz_1[i] * pc_y[i];

        ta1_y_yzzz_xxyyyy_0[i] =
            2.0 * ta1_y_yz_xxyyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxyyyy_1[i] * fe_0 + ta1_y_yzz_xxyyyy_0[i] * pa_z[i] - ta1_y_yzz_xxyyyy_1[i] * pc_z[i];

        ta1_y_yzzz_xxyyyz_0[i] = 3.0 * ta1_y_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxyyz_1[i] * fe_0 + ta_zzz_xxyyyz_1[i] +
                                 ta1_y_zzz_xxyyyz_0[i] * pa_y[i] - ta1_y_zzz_xxyyyz_1[i] * pc_y[i];

        ta1_y_yzzz_xxyyzz_0[i] = 2.0 * ta1_y_zzz_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xxyzz_1[i] * fe_0 + ta_zzz_xxyyzz_1[i] +
                                 ta1_y_zzz_xxyyzz_0[i] * pa_y[i] - ta1_y_zzz_xxyyzz_1[i] * pc_y[i];

        ta1_y_yzzz_xxyzzz_0[i] = ta1_y_zzz_xxzzz_0[i] * fe_0 - ta1_y_zzz_xxzzz_1[i] * fe_0 + ta_zzz_xxyzzz_1[i] + ta1_y_zzz_xxyzzz_0[i] * pa_y[i] -
                                 ta1_y_zzz_xxyzzz_1[i] * pc_y[i];

        ta1_y_yzzz_xxzzzz_0[i] = ta_zzz_xxzzzz_1[i] + ta1_y_zzz_xxzzzz_0[i] * pa_y[i] - ta1_y_zzz_xxzzzz_1[i] * pc_y[i];

        ta1_y_yzzz_xyyyyy_0[i] =
            2.0 * ta1_y_yz_xyyyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xyyyyy_1[i] * fe_0 + ta1_y_yzz_xyyyyy_0[i] * pa_z[i] - ta1_y_yzz_xyyyyy_1[i] * pc_z[i];

        ta1_y_yzzz_xyyyyz_0[i] = 4.0 * ta1_y_zzz_xyyyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xyyyz_1[i] * fe_0 + ta_zzz_xyyyyz_1[i] +
                                 ta1_y_zzz_xyyyyz_0[i] * pa_y[i] - ta1_y_zzz_xyyyyz_1[i] * pc_y[i];

        ta1_y_yzzz_xyyyzz_0[i] = 3.0 * ta1_y_zzz_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xyyzz_1[i] * fe_0 + ta_zzz_xyyyzz_1[i] +
                                 ta1_y_zzz_xyyyzz_0[i] * pa_y[i] - ta1_y_zzz_xyyyzz_1[i] * pc_y[i];

        ta1_y_yzzz_xyyzzz_0[i] = 2.0 * ta1_y_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyzzz_1[i] * fe_0 + ta_zzz_xyyzzz_1[i] +
                                 ta1_y_zzz_xyyzzz_0[i] * pa_y[i] - ta1_y_zzz_xyyzzz_1[i] * pc_y[i];

        ta1_y_yzzz_xyzzzz_0[i] = ta1_y_zzz_xzzzz_0[i] * fe_0 - ta1_y_zzz_xzzzz_1[i] * fe_0 + ta_zzz_xyzzzz_1[i] + ta1_y_zzz_xyzzzz_0[i] * pa_y[i] -
                                 ta1_y_zzz_xyzzzz_1[i] * pc_y[i];

        ta1_y_yzzz_xzzzzz_0[i] = ta_zzz_xzzzzz_1[i] + ta1_y_zzz_xzzzzz_0[i] * pa_y[i] - ta1_y_zzz_xzzzzz_1[i] * pc_y[i];

        ta1_y_yzzz_yyyyyy_0[i] =
            2.0 * ta1_y_yz_yyyyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_yyyyyy_1[i] * fe_0 + ta1_y_yzz_yyyyyy_0[i] * pa_z[i] - ta1_y_yzz_yyyyyy_1[i] * pc_z[i];

        ta1_y_yzzz_yyyyyz_0[i] = 5.0 * ta1_y_zzz_yyyyz_0[i] * fe_0 - 5.0 * ta1_y_zzz_yyyyz_1[i] * fe_0 + ta_zzz_yyyyyz_1[i] +
                                 ta1_y_zzz_yyyyyz_0[i] * pa_y[i] - ta1_y_zzz_yyyyyz_1[i] * pc_y[i];

        ta1_y_yzzz_yyyyzz_0[i] = 4.0 * ta1_y_zzz_yyyzz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yyyzz_1[i] * fe_0 + ta_zzz_yyyyzz_1[i] +
                                 ta1_y_zzz_yyyyzz_0[i] * pa_y[i] - ta1_y_zzz_yyyyzz_1[i] * pc_y[i];

        ta1_y_yzzz_yyyzzz_0[i] = 3.0 * ta1_y_zzz_yyzzz_0[i] * fe_0 - 3.0 * ta1_y_zzz_yyzzz_1[i] * fe_0 + ta_zzz_yyyzzz_1[i] +
                                 ta1_y_zzz_yyyzzz_0[i] * pa_y[i] - ta1_y_zzz_yyyzzz_1[i] * pc_y[i];

        ta1_y_yzzz_yyzzzz_0[i] = 2.0 * ta1_y_zzz_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_yzzzz_1[i] * fe_0 + ta_zzz_yyzzzz_1[i] +
                                 ta1_y_zzz_yyzzzz_0[i] * pa_y[i] - ta1_y_zzz_yyzzzz_1[i] * pc_y[i];

        ta1_y_yzzz_yzzzzz_0[i] = ta1_y_zzz_zzzzz_0[i] * fe_0 - ta1_y_zzz_zzzzz_1[i] * fe_0 + ta_zzz_yzzzzz_1[i] + ta1_y_zzz_yzzzzz_0[i] * pa_y[i] -
                                 ta1_y_zzz_yzzzzz_1[i] * pc_y[i];

        ta1_y_yzzz_zzzzzz_0[i] = ta_zzz_zzzzzz_1[i] + ta1_y_zzz_zzzzzz_0[i] * pa_y[i] - ta1_y_zzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 812-840 components of targeted buffer : GI

    auto ta1_y_zzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 812);

    auto ta1_y_zzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 813);

    auto ta1_y_zzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 814);

    auto ta1_y_zzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 815);

    auto ta1_y_zzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 816);

    auto ta1_y_zzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 817);

    auto ta1_y_zzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 818);

    auto ta1_y_zzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 819);

    auto ta1_y_zzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 820);

    auto ta1_y_zzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 821);

    auto ta1_y_zzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 822);

    auto ta1_y_zzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 823);

    auto ta1_y_zzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 824);

    auto ta1_y_zzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 825);

    auto ta1_y_zzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 826);

    auto ta1_y_zzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 827);

    auto ta1_y_zzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 828);

    auto ta1_y_zzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 829);

    auto ta1_y_zzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 830);

    auto ta1_y_zzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 831);

    auto ta1_y_zzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 832);

    auto ta1_y_zzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 833);

    auto ta1_y_zzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 834);

    auto ta1_y_zzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 835);

    auto ta1_y_zzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 836);

    auto ta1_y_zzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 837);

    auto ta1_y_zzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 838);

    auto ta1_y_zzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 839);

#pragma omp simd aligned(pa_z,                    \
                             pc_z,                \
                             ta1_y_zz_xxxxxx_0,   \
                             ta1_y_zz_xxxxxx_1,   \
                             ta1_y_zz_xxxxxy_0,   \
                             ta1_y_zz_xxxxxy_1,   \
                             ta1_y_zz_xxxxxz_0,   \
                             ta1_y_zz_xxxxxz_1,   \
                             ta1_y_zz_xxxxyy_0,   \
                             ta1_y_zz_xxxxyy_1,   \
                             ta1_y_zz_xxxxyz_0,   \
                             ta1_y_zz_xxxxyz_1,   \
                             ta1_y_zz_xxxxzz_0,   \
                             ta1_y_zz_xxxxzz_1,   \
                             ta1_y_zz_xxxyyy_0,   \
                             ta1_y_zz_xxxyyy_1,   \
                             ta1_y_zz_xxxyyz_0,   \
                             ta1_y_zz_xxxyyz_1,   \
                             ta1_y_zz_xxxyzz_0,   \
                             ta1_y_zz_xxxyzz_1,   \
                             ta1_y_zz_xxxzzz_0,   \
                             ta1_y_zz_xxxzzz_1,   \
                             ta1_y_zz_xxyyyy_0,   \
                             ta1_y_zz_xxyyyy_1,   \
                             ta1_y_zz_xxyyyz_0,   \
                             ta1_y_zz_xxyyyz_1,   \
                             ta1_y_zz_xxyyzz_0,   \
                             ta1_y_zz_xxyyzz_1,   \
                             ta1_y_zz_xxyzzz_0,   \
                             ta1_y_zz_xxyzzz_1,   \
                             ta1_y_zz_xxzzzz_0,   \
                             ta1_y_zz_xxzzzz_1,   \
                             ta1_y_zz_xyyyyy_0,   \
                             ta1_y_zz_xyyyyy_1,   \
                             ta1_y_zz_xyyyyz_0,   \
                             ta1_y_zz_xyyyyz_1,   \
                             ta1_y_zz_xyyyzz_0,   \
                             ta1_y_zz_xyyyzz_1,   \
                             ta1_y_zz_xyyzzz_0,   \
                             ta1_y_zz_xyyzzz_1,   \
                             ta1_y_zz_xyzzzz_0,   \
                             ta1_y_zz_xyzzzz_1,   \
                             ta1_y_zz_xzzzzz_0,   \
                             ta1_y_zz_xzzzzz_1,   \
                             ta1_y_zz_yyyyyy_0,   \
                             ta1_y_zz_yyyyyy_1,   \
                             ta1_y_zz_yyyyyz_0,   \
                             ta1_y_zz_yyyyyz_1,   \
                             ta1_y_zz_yyyyzz_0,   \
                             ta1_y_zz_yyyyzz_1,   \
                             ta1_y_zz_yyyzzz_0,   \
                             ta1_y_zz_yyyzzz_1,   \
                             ta1_y_zz_yyzzzz_0,   \
                             ta1_y_zz_yyzzzz_1,   \
                             ta1_y_zz_yzzzzz_0,   \
                             ta1_y_zz_yzzzzz_1,   \
                             ta1_y_zz_zzzzzz_0,   \
                             ta1_y_zz_zzzzzz_1,   \
                             ta1_y_zzz_xxxxx_0,   \
                             ta1_y_zzz_xxxxx_1,   \
                             ta1_y_zzz_xxxxxx_0,  \
                             ta1_y_zzz_xxxxxx_1,  \
                             ta1_y_zzz_xxxxxy_0,  \
                             ta1_y_zzz_xxxxxy_1,  \
                             ta1_y_zzz_xxxxxz_0,  \
                             ta1_y_zzz_xxxxxz_1,  \
                             ta1_y_zzz_xxxxy_0,   \
                             ta1_y_zzz_xxxxy_1,   \
                             ta1_y_zzz_xxxxyy_0,  \
                             ta1_y_zzz_xxxxyy_1,  \
                             ta1_y_zzz_xxxxyz_0,  \
                             ta1_y_zzz_xxxxyz_1,  \
                             ta1_y_zzz_xxxxz_0,   \
                             ta1_y_zzz_xxxxz_1,   \
                             ta1_y_zzz_xxxxzz_0,  \
                             ta1_y_zzz_xxxxzz_1,  \
                             ta1_y_zzz_xxxyy_0,   \
                             ta1_y_zzz_xxxyy_1,   \
                             ta1_y_zzz_xxxyyy_0,  \
                             ta1_y_zzz_xxxyyy_1,  \
                             ta1_y_zzz_xxxyyz_0,  \
                             ta1_y_zzz_xxxyyz_1,  \
                             ta1_y_zzz_xxxyz_0,   \
                             ta1_y_zzz_xxxyz_1,   \
                             ta1_y_zzz_xxxyzz_0,  \
                             ta1_y_zzz_xxxyzz_1,  \
                             ta1_y_zzz_xxxzz_0,   \
                             ta1_y_zzz_xxxzz_1,   \
                             ta1_y_zzz_xxxzzz_0,  \
                             ta1_y_zzz_xxxzzz_1,  \
                             ta1_y_zzz_xxyyy_0,   \
                             ta1_y_zzz_xxyyy_1,   \
                             ta1_y_zzz_xxyyyy_0,  \
                             ta1_y_zzz_xxyyyy_1,  \
                             ta1_y_zzz_xxyyyz_0,  \
                             ta1_y_zzz_xxyyyz_1,  \
                             ta1_y_zzz_xxyyz_0,   \
                             ta1_y_zzz_xxyyz_1,   \
                             ta1_y_zzz_xxyyzz_0,  \
                             ta1_y_zzz_xxyyzz_1,  \
                             ta1_y_zzz_xxyzz_0,   \
                             ta1_y_zzz_xxyzz_1,   \
                             ta1_y_zzz_xxyzzz_0,  \
                             ta1_y_zzz_xxyzzz_1,  \
                             ta1_y_zzz_xxzzz_0,   \
                             ta1_y_zzz_xxzzz_1,   \
                             ta1_y_zzz_xxzzzz_0,  \
                             ta1_y_zzz_xxzzzz_1,  \
                             ta1_y_zzz_xyyyy_0,   \
                             ta1_y_zzz_xyyyy_1,   \
                             ta1_y_zzz_xyyyyy_0,  \
                             ta1_y_zzz_xyyyyy_1,  \
                             ta1_y_zzz_xyyyyz_0,  \
                             ta1_y_zzz_xyyyyz_1,  \
                             ta1_y_zzz_xyyyz_0,   \
                             ta1_y_zzz_xyyyz_1,   \
                             ta1_y_zzz_xyyyzz_0,  \
                             ta1_y_zzz_xyyyzz_1,  \
                             ta1_y_zzz_xyyzz_0,   \
                             ta1_y_zzz_xyyzz_1,   \
                             ta1_y_zzz_xyyzzz_0,  \
                             ta1_y_zzz_xyyzzz_1,  \
                             ta1_y_zzz_xyzzz_0,   \
                             ta1_y_zzz_xyzzz_1,   \
                             ta1_y_zzz_xyzzzz_0,  \
                             ta1_y_zzz_xyzzzz_1,  \
                             ta1_y_zzz_xzzzz_0,   \
                             ta1_y_zzz_xzzzz_1,   \
                             ta1_y_zzz_xzzzzz_0,  \
                             ta1_y_zzz_xzzzzz_1,  \
                             ta1_y_zzz_yyyyy_0,   \
                             ta1_y_zzz_yyyyy_1,   \
                             ta1_y_zzz_yyyyyy_0,  \
                             ta1_y_zzz_yyyyyy_1,  \
                             ta1_y_zzz_yyyyyz_0,  \
                             ta1_y_zzz_yyyyyz_1,  \
                             ta1_y_zzz_yyyyz_0,   \
                             ta1_y_zzz_yyyyz_1,   \
                             ta1_y_zzz_yyyyzz_0,  \
                             ta1_y_zzz_yyyyzz_1,  \
                             ta1_y_zzz_yyyzz_0,   \
                             ta1_y_zzz_yyyzz_1,   \
                             ta1_y_zzz_yyyzzz_0,  \
                             ta1_y_zzz_yyyzzz_1,  \
                             ta1_y_zzz_yyzzz_0,   \
                             ta1_y_zzz_yyzzz_1,   \
                             ta1_y_zzz_yyzzzz_0,  \
                             ta1_y_zzz_yyzzzz_1,  \
                             ta1_y_zzz_yzzzz_0,   \
                             ta1_y_zzz_yzzzz_1,   \
                             ta1_y_zzz_yzzzzz_0,  \
                             ta1_y_zzz_yzzzzz_1,  \
                             ta1_y_zzz_zzzzz_0,   \
                             ta1_y_zzz_zzzzz_1,   \
                             ta1_y_zzz_zzzzzz_0,  \
                             ta1_y_zzz_zzzzzz_1,  \
                             ta1_y_zzzz_xxxxxx_0, \
                             ta1_y_zzzz_xxxxxy_0, \
                             ta1_y_zzzz_xxxxxz_0, \
                             ta1_y_zzzz_xxxxyy_0, \
                             ta1_y_zzzz_xxxxyz_0, \
                             ta1_y_zzzz_xxxxzz_0, \
                             ta1_y_zzzz_xxxyyy_0, \
                             ta1_y_zzzz_xxxyyz_0, \
                             ta1_y_zzzz_xxxyzz_0, \
                             ta1_y_zzzz_xxxzzz_0, \
                             ta1_y_zzzz_xxyyyy_0, \
                             ta1_y_zzzz_xxyyyz_0, \
                             ta1_y_zzzz_xxyyzz_0, \
                             ta1_y_zzzz_xxyzzz_0, \
                             ta1_y_zzzz_xxzzzz_0, \
                             ta1_y_zzzz_xyyyyy_0, \
                             ta1_y_zzzz_xyyyyz_0, \
                             ta1_y_zzzz_xyyyzz_0, \
                             ta1_y_zzzz_xyyzzz_0, \
                             ta1_y_zzzz_xyzzzz_0, \
                             ta1_y_zzzz_xzzzzz_0, \
                             ta1_y_zzzz_yyyyyy_0, \
                             ta1_y_zzzz_yyyyyz_0, \
                             ta1_y_zzzz_yyyyzz_0, \
                             ta1_y_zzzz_yyyzzz_0, \
                             ta1_y_zzzz_yyzzzz_0, \
                             ta1_y_zzzz_yzzzzz_0, \
                             ta1_y_zzzz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzz_xxxxxx_0[i] =
            3.0 * ta1_y_zz_xxxxxx_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxxx_1[i] * fe_0 + ta1_y_zzz_xxxxxx_0[i] * pa_z[i] - ta1_y_zzz_xxxxxx_1[i] * pc_z[i];

        ta1_y_zzzz_xxxxxy_0[i] =
            3.0 * ta1_y_zz_xxxxxy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxxy_1[i] * fe_0 + ta1_y_zzz_xxxxxy_0[i] * pa_z[i] - ta1_y_zzz_xxxxxy_1[i] * pc_z[i];

        ta1_y_zzzz_xxxxxz_0[i] = 3.0 * ta1_y_zz_xxxxxz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxxz_1[i] * fe_0 + ta1_y_zzz_xxxxx_0[i] * fe_0 -
                                 ta1_y_zzz_xxxxx_1[i] * fe_0 + ta1_y_zzz_xxxxxz_0[i] * pa_z[i] - ta1_y_zzz_xxxxxz_1[i] * pc_z[i];

        ta1_y_zzzz_xxxxyy_0[i] =
            3.0 * ta1_y_zz_xxxxyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxyy_1[i] * fe_0 + ta1_y_zzz_xxxxyy_0[i] * pa_z[i] - ta1_y_zzz_xxxxyy_1[i] * pc_z[i];

        ta1_y_zzzz_xxxxyz_0[i] = 3.0 * ta1_y_zz_xxxxyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxyz_1[i] * fe_0 + ta1_y_zzz_xxxxy_0[i] * fe_0 -
                                 ta1_y_zzz_xxxxy_1[i] * fe_0 + ta1_y_zzz_xxxxyz_0[i] * pa_z[i] - ta1_y_zzz_xxxxyz_1[i] * pc_z[i];

        ta1_y_zzzz_xxxxzz_0[i] = 3.0 * ta1_y_zz_xxxxzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xxxxz_0[i] * fe_0 -
                                 2.0 * ta1_y_zzz_xxxxz_1[i] * fe_0 + ta1_y_zzz_xxxxzz_0[i] * pa_z[i] - ta1_y_zzz_xxxxzz_1[i] * pc_z[i];

        ta1_y_zzzz_xxxyyy_0[i] =
            3.0 * ta1_y_zz_xxxyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxyyy_1[i] * fe_0 + ta1_y_zzz_xxxyyy_0[i] * pa_z[i] - ta1_y_zzz_xxxyyy_1[i] * pc_z[i];

        ta1_y_zzzz_xxxyyz_0[i] = 3.0 * ta1_y_zz_xxxyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxyyz_1[i] * fe_0 + ta1_y_zzz_xxxyy_0[i] * fe_0 -
                                 ta1_y_zzz_xxxyy_1[i] * fe_0 + ta1_y_zzz_xxxyyz_0[i] * pa_z[i] - ta1_y_zzz_xxxyyz_1[i] * pc_z[i];

        ta1_y_zzzz_xxxyzz_0[i] = 3.0 * ta1_y_zz_xxxyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_y_zzz_xxxyz_1[i] * fe_0 + ta1_y_zzz_xxxyzz_0[i] * pa_z[i] - ta1_y_zzz_xxxyzz_1[i] * pc_z[i];

        ta1_y_zzzz_xxxzzz_0[i] = 3.0 * ta1_y_zz_xxxzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_xxxzz_0[i] * fe_0 -
                                 3.0 * ta1_y_zzz_xxxzz_1[i] * fe_0 + ta1_y_zzz_xxxzzz_0[i] * pa_z[i] - ta1_y_zzz_xxxzzz_1[i] * pc_z[i];

        ta1_y_zzzz_xxyyyy_0[i] =
            3.0 * ta1_y_zz_xxyyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyyy_1[i] * fe_0 + ta1_y_zzz_xxyyyy_0[i] * pa_z[i] - ta1_y_zzz_xxyyyy_1[i] * pc_z[i];

        ta1_y_zzzz_xxyyyz_0[i] = 3.0 * ta1_y_zz_xxyyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyyz_1[i] * fe_0 + ta1_y_zzz_xxyyy_0[i] * fe_0 -
                                 ta1_y_zzz_xxyyy_1[i] * fe_0 + ta1_y_zzz_xxyyyz_0[i] * pa_z[i] - ta1_y_zzz_xxyyyz_1[i] * pc_z[i];

        ta1_y_zzzz_xxyyzz_0[i] = 3.0 * ta1_y_zz_xxyyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xxyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_zzz_xxyyz_1[i] * fe_0 + ta1_y_zzz_xxyyzz_0[i] * pa_z[i] - ta1_y_zzz_xxyyzz_1[i] * pc_z[i];

        ta1_y_zzzz_xxyzzz_0[i] = 3.0 * ta1_y_zz_xxyzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_zzz_xxyzz_1[i] * fe_0 + ta1_y_zzz_xxyzzz_0[i] * pa_z[i] - ta1_y_zzz_xxyzzz_1[i] * pc_z[i];

        ta1_y_zzzz_xxzzzz_0[i] = 3.0 * ta1_y_zz_xxzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxzzzz_1[i] * fe_0 + 4.0 * ta1_y_zzz_xxzzz_0[i] * fe_0 -
                                 4.0 * ta1_y_zzz_xxzzz_1[i] * fe_0 + ta1_y_zzz_xxzzzz_0[i] * pa_z[i] - ta1_y_zzz_xxzzzz_1[i] * pc_z[i];

        ta1_y_zzzz_xyyyyy_0[i] =
            3.0 * ta1_y_zz_xyyyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyyyy_1[i] * fe_0 + ta1_y_zzz_xyyyyy_0[i] * pa_z[i] - ta1_y_zzz_xyyyyy_1[i] * pc_z[i];

        ta1_y_zzzz_xyyyyz_0[i] = 3.0 * ta1_y_zz_xyyyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyyyz_1[i] * fe_0 + ta1_y_zzz_xyyyy_0[i] * fe_0 -
                                 ta1_y_zzz_xyyyy_1[i] * fe_0 + ta1_y_zzz_xyyyyz_0[i] * pa_z[i] - ta1_y_zzz_xyyyyz_1[i] * pc_z[i];

        ta1_y_zzzz_xyyyzz_0[i] = 3.0 * ta1_y_zz_xyyyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_zzz_xyyyz_1[i] * fe_0 + ta1_y_zzz_xyyyzz_0[i] * pa_z[i] - ta1_y_zzz_xyyyzz_1[i] * pc_z[i];

        ta1_y_zzzz_xyyzzz_0[i] = 3.0 * ta1_y_zz_xyyzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_zzz_xyyzz_1[i] * fe_0 + ta1_y_zzz_xyyzzz_0[i] * pa_z[i] - ta1_y_zzz_xyyzzz_1[i] * pc_z[i];

        ta1_y_zzzz_xyzzzz_0[i] = 3.0 * ta1_y_zz_xyzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyzzzz_1[i] * fe_0 + 4.0 * ta1_y_zzz_xyzzz_0[i] * fe_0 -
                                 4.0 * ta1_y_zzz_xyzzz_1[i] * fe_0 + ta1_y_zzz_xyzzzz_0[i] * pa_z[i] - ta1_y_zzz_xyzzzz_1[i] * pc_z[i];

        ta1_y_zzzz_xzzzzz_0[i] = 3.0 * ta1_y_zz_xzzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xzzzzz_1[i] * fe_0 + 5.0 * ta1_y_zzz_xzzzz_0[i] * fe_0 -
                                 5.0 * ta1_y_zzz_xzzzz_1[i] * fe_0 + ta1_y_zzz_xzzzzz_0[i] * pa_z[i] - ta1_y_zzz_xzzzzz_1[i] * pc_z[i];

        ta1_y_zzzz_yyyyyy_0[i] =
            3.0 * ta1_y_zz_yyyyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyyyy_1[i] * fe_0 + ta1_y_zzz_yyyyyy_0[i] * pa_z[i] - ta1_y_zzz_yyyyyy_1[i] * pc_z[i];

        ta1_y_zzzz_yyyyyz_0[i] = 3.0 * ta1_y_zz_yyyyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyyyz_1[i] * fe_0 + ta1_y_zzz_yyyyy_0[i] * fe_0 -
                                 ta1_y_zzz_yyyyy_1[i] * fe_0 + ta1_y_zzz_yyyyyz_0[i] * pa_z[i] - ta1_y_zzz_yyyyyz_1[i] * pc_z[i];

        ta1_y_zzzz_yyyyzz_0[i] = 3.0 * ta1_y_zz_yyyyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yyyyz_0[i] * fe_0 -
                                 2.0 * ta1_y_zzz_yyyyz_1[i] * fe_0 + ta1_y_zzz_yyyyzz_0[i] * pa_z[i] - ta1_y_zzz_yyyyzz_1[i] * pc_z[i];

        ta1_y_zzzz_yyyzzz_0[i] = 3.0 * ta1_y_zz_yyyzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_yyyzz_0[i] * fe_0 -
                                 3.0 * ta1_y_zzz_yyyzz_1[i] * fe_0 + ta1_y_zzz_yyyzzz_0[i] * pa_z[i] - ta1_y_zzz_yyyzzz_1[i] * pc_z[i];

        ta1_y_zzzz_yyzzzz_0[i] = 3.0 * ta1_y_zz_yyzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyzzzz_1[i] * fe_0 + 4.0 * ta1_y_zzz_yyzzz_0[i] * fe_0 -
                                 4.0 * ta1_y_zzz_yyzzz_1[i] * fe_0 + ta1_y_zzz_yyzzzz_0[i] * pa_z[i] - ta1_y_zzz_yyzzzz_1[i] * pc_z[i];

        ta1_y_zzzz_yzzzzz_0[i] = 3.0 * ta1_y_zz_yzzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yzzzzz_1[i] * fe_0 + 5.0 * ta1_y_zzz_yzzzz_0[i] * fe_0 -
                                 5.0 * ta1_y_zzz_yzzzz_1[i] * fe_0 + ta1_y_zzz_yzzzzz_0[i] * pa_z[i] - ta1_y_zzz_yzzzzz_1[i] * pc_z[i];

        ta1_y_zzzz_zzzzzz_0[i] = 3.0 * ta1_y_zz_zzzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_zzzzzz_1[i] * fe_0 + 6.0 * ta1_y_zzz_zzzzz_0[i] * fe_0 -
                                 6.0 * ta1_y_zzz_zzzzz_1[i] * fe_0 + ta1_y_zzz_zzzzzz_0[i] * pa_z[i] - ta1_y_zzz_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 840-868 components of targeted buffer : GI

    auto ta1_z_xxxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 840);

    auto ta1_z_xxxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 841);

    auto ta1_z_xxxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 842);

    auto ta1_z_xxxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 843);

    auto ta1_z_xxxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 844);

    auto ta1_z_xxxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 845);

    auto ta1_z_xxxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 846);

    auto ta1_z_xxxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 847);

    auto ta1_z_xxxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 848);

    auto ta1_z_xxxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 849);

    auto ta1_z_xxxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 850);

    auto ta1_z_xxxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 851);

    auto ta1_z_xxxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 852);

    auto ta1_z_xxxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 853);

    auto ta1_z_xxxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 854);

    auto ta1_z_xxxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 855);

    auto ta1_z_xxxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 856);

    auto ta1_z_xxxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 857);

    auto ta1_z_xxxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 858);

    auto ta1_z_xxxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 859);

    auto ta1_z_xxxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 860);

    auto ta1_z_xxxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 861);

    auto ta1_z_xxxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 862);

    auto ta1_z_xxxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 863);

    auto ta1_z_xxxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 864);

    auto ta1_z_xxxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 865);

    auto ta1_z_xxxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 866);

    auto ta1_z_xxxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 867);

#pragma omp simd aligned(pa_x,                    \
                             pc_x,                \
                             ta1_z_xx_xxxxxx_0,   \
                             ta1_z_xx_xxxxxx_1,   \
                             ta1_z_xx_xxxxxy_0,   \
                             ta1_z_xx_xxxxxy_1,   \
                             ta1_z_xx_xxxxxz_0,   \
                             ta1_z_xx_xxxxxz_1,   \
                             ta1_z_xx_xxxxyy_0,   \
                             ta1_z_xx_xxxxyy_1,   \
                             ta1_z_xx_xxxxyz_0,   \
                             ta1_z_xx_xxxxyz_1,   \
                             ta1_z_xx_xxxxzz_0,   \
                             ta1_z_xx_xxxxzz_1,   \
                             ta1_z_xx_xxxyyy_0,   \
                             ta1_z_xx_xxxyyy_1,   \
                             ta1_z_xx_xxxyyz_0,   \
                             ta1_z_xx_xxxyyz_1,   \
                             ta1_z_xx_xxxyzz_0,   \
                             ta1_z_xx_xxxyzz_1,   \
                             ta1_z_xx_xxxzzz_0,   \
                             ta1_z_xx_xxxzzz_1,   \
                             ta1_z_xx_xxyyyy_0,   \
                             ta1_z_xx_xxyyyy_1,   \
                             ta1_z_xx_xxyyyz_0,   \
                             ta1_z_xx_xxyyyz_1,   \
                             ta1_z_xx_xxyyzz_0,   \
                             ta1_z_xx_xxyyzz_1,   \
                             ta1_z_xx_xxyzzz_0,   \
                             ta1_z_xx_xxyzzz_1,   \
                             ta1_z_xx_xxzzzz_0,   \
                             ta1_z_xx_xxzzzz_1,   \
                             ta1_z_xx_xyyyyy_0,   \
                             ta1_z_xx_xyyyyy_1,   \
                             ta1_z_xx_xyyyyz_0,   \
                             ta1_z_xx_xyyyyz_1,   \
                             ta1_z_xx_xyyyzz_0,   \
                             ta1_z_xx_xyyyzz_1,   \
                             ta1_z_xx_xyyzzz_0,   \
                             ta1_z_xx_xyyzzz_1,   \
                             ta1_z_xx_xyzzzz_0,   \
                             ta1_z_xx_xyzzzz_1,   \
                             ta1_z_xx_xzzzzz_0,   \
                             ta1_z_xx_xzzzzz_1,   \
                             ta1_z_xx_yyyyyy_0,   \
                             ta1_z_xx_yyyyyy_1,   \
                             ta1_z_xx_yyyyyz_0,   \
                             ta1_z_xx_yyyyyz_1,   \
                             ta1_z_xx_yyyyzz_0,   \
                             ta1_z_xx_yyyyzz_1,   \
                             ta1_z_xx_yyyzzz_0,   \
                             ta1_z_xx_yyyzzz_1,   \
                             ta1_z_xx_yyzzzz_0,   \
                             ta1_z_xx_yyzzzz_1,   \
                             ta1_z_xx_yzzzzz_0,   \
                             ta1_z_xx_yzzzzz_1,   \
                             ta1_z_xx_zzzzzz_0,   \
                             ta1_z_xx_zzzzzz_1,   \
                             ta1_z_xxx_xxxxx_0,   \
                             ta1_z_xxx_xxxxx_1,   \
                             ta1_z_xxx_xxxxxx_0,  \
                             ta1_z_xxx_xxxxxx_1,  \
                             ta1_z_xxx_xxxxxy_0,  \
                             ta1_z_xxx_xxxxxy_1,  \
                             ta1_z_xxx_xxxxxz_0,  \
                             ta1_z_xxx_xxxxxz_1,  \
                             ta1_z_xxx_xxxxy_0,   \
                             ta1_z_xxx_xxxxy_1,   \
                             ta1_z_xxx_xxxxyy_0,  \
                             ta1_z_xxx_xxxxyy_1,  \
                             ta1_z_xxx_xxxxyz_0,  \
                             ta1_z_xxx_xxxxyz_1,  \
                             ta1_z_xxx_xxxxz_0,   \
                             ta1_z_xxx_xxxxz_1,   \
                             ta1_z_xxx_xxxxzz_0,  \
                             ta1_z_xxx_xxxxzz_1,  \
                             ta1_z_xxx_xxxyy_0,   \
                             ta1_z_xxx_xxxyy_1,   \
                             ta1_z_xxx_xxxyyy_0,  \
                             ta1_z_xxx_xxxyyy_1,  \
                             ta1_z_xxx_xxxyyz_0,  \
                             ta1_z_xxx_xxxyyz_1,  \
                             ta1_z_xxx_xxxyz_0,   \
                             ta1_z_xxx_xxxyz_1,   \
                             ta1_z_xxx_xxxyzz_0,  \
                             ta1_z_xxx_xxxyzz_1,  \
                             ta1_z_xxx_xxxzz_0,   \
                             ta1_z_xxx_xxxzz_1,   \
                             ta1_z_xxx_xxxzzz_0,  \
                             ta1_z_xxx_xxxzzz_1,  \
                             ta1_z_xxx_xxyyy_0,   \
                             ta1_z_xxx_xxyyy_1,   \
                             ta1_z_xxx_xxyyyy_0,  \
                             ta1_z_xxx_xxyyyy_1,  \
                             ta1_z_xxx_xxyyyz_0,  \
                             ta1_z_xxx_xxyyyz_1,  \
                             ta1_z_xxx_xxyyz_0,   \
                             ta1_z_xxx_xxyyz_1,   \
                             ta1_z_xxx_xxyyzz_0,  \
                             ta1_z_xxx_xxyyzz_1,  \
                             ta1_z_xxx_xxyzz_0,   \
                             ta1_z_xxx_xxyzz_1,   \
                             ta1_z_xxx_xxyzzz_0,  \
                             ta1_z_xxx_xxyzzz_1,  \
                             ta1_z_xxx_xxzzz_0,   \
                             ta1_z_xxx_xxzzz_1,   \
                             ta1_z_xxx_xxzzzz_0,  \
                             ta1_z_xxx_xxzzzz_1,  \
                             ta1_z_xxx_xyyyy_0,   \
                             ta1_z_xxx_xyyyy_1,   \
                             ta1_z_xxx_xyyyyy_0,  \
                             ta1_z_xxx_xyyyyy_1,  \
                             ta1_z_xxx_xyyyyz_0,  \
                             ta1_z_xxx_xyyyyz_1,  \
                             ta1_z_xxx_xyyyz_0,   \
                             ta1_z_xxx_xyyyz_1,   \
                             ta1_z_xxx_xyyyzz_0,  \
                             ta1_z_xxx_xyyyzz_1,  \
                             ta1_z_xxx_xyyzz_0,   \
                             ta1_z_xxx_xyyzz_1,   \
                             ta1_z_xxx_xyyzzz_0,  \
                             ta1_z_xxx_xyyzzz_1,  \
                             ta1_z_xxx_xyzzz_0,   \
                             ta1_z_xxx_xyzzz_1,   \
                             ta1_z_xxx_xyzzzz_0,  \
                             ta1_z_xxx_xyzzzz_1,  \
                             ta1_z_xxx_xzzzz_0,   \
                             ta1_z_xxx_xzzzz_1,   \
                             ta1_z_xxx_xzzzzz_0,  \
                             ta1_z_xxx_xzzzzz_1,  \
                             ta1_z_xxx_yyyyy_0,   \
                             ta1_z_xxx_yyyyy_1,   \
                             ta1_z_xxx_yyyyyy_0,  \
                             ta1_z_xxx_yyyyyy_1,  \
                             ta1_z_xxx_yyyyyz_0,  \
                             ta1_z_xxx_yyyyyz_1,  \
                             ta1_z_xxx_yyyyz_0,   \
                             ta1_z_xxx_yyyyz_1,   \
                             ta1_z_xxx_yyyyzz_0,  \
                             ta1_z_xxx_yyyyzz_1,  \
                             ta1_z_xxx_yyyzz_0,   \
                             ta1_z_xxx_yyyzz_1,   \
                             ta1_z_xxx_yyyzzz_0,  \
                             ta1_z_xxx_yyyzzz_1,  \
                             ta1_z_xxx_yyzzz_0,   \
                             ta1_z_xxx_yyzzz_1,   \
                             ta1_z_xxx_yyzzzz_0,  \
                             ta1_z_xxx_yyzzzz_1,  \
                             ta1_z_xxx_yzzzz_0,   \
                             ta1_z_xxx_yzzzz_1,   \
                             ta1_z_xxx_yzzzzz_0,  \
                             ta1_z_xxx_yzzzzz_1,  \
                             ta1_z_xxx_zzzzz_0,   \
                             ta1_z_xxx_zzzzz_1,   \
                             ta1_z_xxx_zzzzzz_0,  \
                             ta1_z_xxx_zzzzzz_1,  \
                             ta1_z_xxxx_xxxxxx_0, \
                             ta1_z_xxxx_xxxxxy_0, \
                             ta1_z_xxxx_xxxxxz_0, \
                             ta1_z_xxxx_xxxxyy_0, \
                             ta1_z_xxxx_xxxxyz_0, \
                             ta1_z_xxxx_xxxxzz_0, \
                             ta1_z_xxxx_xxxyyy_0, \
                             ta1_z_xxxx_xxxyyz_0, \
                             ta1_z_xxxx_xxxyzz_0, \
                             ta1_z_xxxx_xxxzzz_0, \
                             ta1_z_xxxx_xxyyyy_0, \
                             ta1_z_xxxx_xxyyyz_0, \
                             ta1_z_xxxx_xxyyzz_0, \
                             ta1_z_xxxx_xxyzzz_0, \
                             ta1_z_xxxx_xxzzzz_0, \
                             ta1_z_xxxx_xyyyyy_0, \
                             ta1_z_xxxx_xyyyyz_0, \
                             ta1_z_xxxx_xyyyzz_0, \
                             ta1_z_xxxx_xyyzzz_0, \
                             ta1_z_xxxx_xyzzzz_0, \
                             ta1_z_xxxx_xzzzzz_0, \
                             ta1_z_xxxx_yyyyyy_0, \
                             ta1_z_xxxx_yyyyyz_0, \
                             ta1_z_xxxx_yyyyzz_0, \
                             ta1_z_xxxx_yyyzzz_0, \
                             ta1_z_xxxx_yyzzzz_0, \
                             ta1_z_xxxx_yzzzzz_0, \
                             ta1_z_xxxx_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxx_xxxxxx_0[i] = 3.0 * ta1_z_xx_xxxxxx_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxxx_1[i] * fe_0 + 6.0 * ta1_z_xxx_xxxxx_0[i] * fe_0 -
                                 6.0 * ta1_z_xxx_xxxxx_1[i] * fe_0 + ta1_z_xxx_xxxxxx_0[i] * pa_x[i] - ta1_z_xxx_xxxxxx_1[i] * pc_x[i];

        ta1_z_xxxx_xxxxxy_0[i] = 3.0 * ta1_z_xx_xxxxxy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxxy_1[i] * fe_0 + 5.0 * ta1_z_xxx_xxxxy_0[i] * fe_0 -
                                 5.0 * ta1_z_xxx_xxxxy_1[i] * fe_0 + ta1_z_xxx_xxxxxy_0[i] * pa_x[i] - ta1_z_xxx_xxxxxy_1[i] * pc_x[i];

        ta1_z_xxxx_xxxxxz_0[i] = 3.0 * ta1_z_xx_xxxxxz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxxz_1[i] * fe_0 + 5.0 * ta1_z_xxx_xxxxz_0[i] * fe_0 -
                                 5.0 * ta1_z_xxx_xxxxz_1[i] * fe_0 + ta1_z_xxx_xxxxxz_0[i] * pa_x[i] - ta1_z_xxx_xxxxxz_1[i] * pc_x[i];

        ta1_z_xxxx_xxxxyy_0[i] = 3.0 * ta1_z_xx_xxxxyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxyy_1[i] * fe_0 + 4.0 * ta1_z_xxx_xxxyy_0[i] * fe_0 -
                                 4.0 * ta1_z_xxx_xxxyy_1[i] * fe_0 + ta1_z_xxx_xxxxyy_0[i] * pa_x[i] - ta1_z_xxx_xxxxyy_1[i] * pc_x[i];

        ta1_z_xxxx_xxxxyz_0[i] = 3.0 * ta1_z_xx_xxxxyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxyz_1[i] * fe_0 + 4.0 * ta1_z_xxx_xxxyz_0[i] * fe_0 -
                                 4.0 * ta1_z_xxx_xxxyz_1[i] * fe_0 + ta1_z_xxx_xxxxyz_0[i] * pa_x[i] - ta1_z_xxx_xxxxyz_1[i] * pc_x[i];

        ta1_z_xxxx_xxxxzz_0[i] = 3.0 * ta1_z_xx_xxxxzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxzz_1[i] * fe_0 + 4.0 * ta1_z_xxx_xxxzz_0[i] * fe_0 -
                                 4.0 * ta1_z_xxx_xxxzz_1[i] * fe_0 + ta1_z_xxx_xxxxzz_0[i] * pa_x[i] - ta1_z_xxx_xxxxzz_1[i] * pc_x[i];

        ta1_z_xxxx_xxxyyy_0[i] = 3.0 * ta1_z_xx_xxxyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxyyy_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxyyy_0[i] * fe_0 -
                                 3.0 * ta1_z_xxx_xxyyy_1[i] * fe_0 + ta1_z_xxx_xxxyyy_0[i] * pa_x[i] - ta1_z_xxx_xxxyyy_1[i] * pc_x[i];

        ta1_z_xxxx_xxxyyz_0[i] = 3.0 * ta1_z_xx_xxxyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxyyz_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_z_xxx_xxyyz_1[i] * fe_0 + ta1_z_xxx_xxxyyz_0[i] * pa_x[i] - ta1_z_xxx_xxxyyz_1[i] * pc_x[i];

        ta1_z_xxxx_xxxyzz_0[i] = 3.0 * ta1_z_xx_xxxyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxyzz_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_xxx_xxyzz_1[i] * fe_0 + ta1_z_xxx_xxxyzz_0[i] * pa_x[i] - ta1_z_xxx_xxxyzz_1[i] * pc_x[i];

        ta1_z_xxxx_xxxzzz_0[i] = 3.0 * ta1_z_xx_xxxzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxzzz_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxzzz_0[i] * fe_0 -
                                 3.0 * ta1_z_xxx_xxzzz_1[i] * fe_0 + ta1_z_xxx_xxxzzz_0[i] * pa_x[i] - ta1_z_xxx_xxxzzz_1[i] * pc_x[i];

        ta1_z_xxxx_xxyyyy_0[i] = 3.0 * ta1_z_xx_xxyyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyyyy_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyyyy_0[i] * fe_0 -
                                 2.0 * ta1_z_xxx_xyyyy_1[i] * fe_0 + ta1_z_xxx_xxyyyy_0[i] * pa_x[i] - ta1_z_xxx_xxyyyy_1[i] * pc_x[i];

        ta1_z_xxxx_xxyyyz_0[i] = 3.0 * ta1_z_xx_xxyyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyyyz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_z_xxx_xyyyz_1[i] * fe_0 + ta1_z_xxx_xxyyyz_0[i] * pa_x[i] - ta1_z_xxx_xxyyyz_1[i] * pc_x[i];

        ta1_z_xxxx_xxyyzz_0[i] = 3.0 * ta1_z_xx_xxyyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyyzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xxx_xyyzz_1[i] * fe_0 + ta1_z_xxx_xxyyzz_0[i] * pa_x[i] - ta1_z_xxx_xxyyzz_1[i] * pc_x[i];

        ta1_z_xxxx_xxyzzz_0[i] = 3.0 * ta1_z_xx_xxyzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyzzz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xxx_xyzzz_1[i] * fe_0 + ta1_z_xxx_xxyzzz_0[i] * pa_x[i] - ta1_z_xxx_xxyzzz_1[i] * pc_x[i];

        ta1_z_xxxx_xxzzzz_0[i] = 3.0 * ta1_z_xx_xxzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxzzzz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xzzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xxx_xzzzz_1[i] * fe_0 + ta1_z_xxx_xxzzzz_0[i] * pa_x[i] - ta1_z_xxx_xxzzzz_1[i] * pc_x[i];

        ta1_z_xxxx_xyyyyy_0[i] = 3.0 * ta1_z_xx_xyyyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyyyy_1[i] * fe_0 + ta1_z_xxx_yyyyy_0[i] * fe_0 -
                                 ta1_z_xxx_yyyyy_1[i] * fe_0 + ta1_z_xxx_xyyyyy_0[i] * pa_x[i] - ta1_z_xxx_xyyyyy_1[i] * pc_x[i];

        ta1_z_xxxx_xyyyyz_0[i] = 3.0 * ta1_z_xx_xyyyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyyyz_1[i] * fe_0 + ta1_z_xxx_yyyyz_0[i] * fe_0 -
                                 ta1_z_xxx_yyyyz_1[i] * fe_0 + ta1_z_xxx_xyyyyz_0[i] * pa_x[i] - ta1_z_xxx_xyyyyz_1[i] * pc_x[i];

        ta1_z_xxxx_xyyyzz_0[i] = 3.0 * ta1_z_xx_xyyyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyyzz_1[i] * fe_0 + ta1_z_xxx_yyyzz_0[i] * fe_0 -
                                 ta1_z_xxx_yyyzz_1[i] * fe_0 + ta1_z_xxx_xyyyzz_0[i] * pa_x[i] - ta1_z_xxx_xyyyzz_1[i] * pc_x[i];

        ta1_z_xxxx_xyyzzz_0[i] = 3.0 * ta1_z_xx_xyyzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyzzz_1[i] * fe_0 + ta1_z_xxx_yyzzz_0[i] * fe_0 -
                                 ta1_z_xxx_yyzzz_1[i] * fe_0 + ta1_z_xxx_xyyzzz_0[i] * pa_x[i] - ta1_z_xxx_xyyzzz_1[i] * pc_x[i];

        ta1_z_xxxx_xyzzzz_0[i] = 3.0 * ta1_z_xx_xyzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyzzzz_1[i] * fe_0 + ta1_z_xxx_yzzzz_0[i] * fe_0 -
                                 ta1_z_xxx_yzzzz_1[i] * fe_0 + ta1_z_xxx_xyzzzz_0[i] * pa_x[i] - ta1_z_xxx_xyzzzz_1[i] * pc_x[i];

        ta1_z_xxxx_xzzzzz_0[i] = 3.0 * ta1_z_xx_xzzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xzzzzz_1[i] * fe_0 + ta1_z_xxx_zzzzz_0[i] * fe_0 -
                                 ta1_z_xxx_zzzzz_1[i] * fe_0 + ta1_z_xxx_xzzzzz_0[i] * pa_x[i] - ta1_z_xxx_xzzzzz_1[i] * pc_x[i];

        ta1_z_xxxx_yyyyyy_0[i] =
            3.0 * ta1_z_xx_yyyyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyyyy_1[i] * fe_0 + ta1_z_xxx_yyyyyy_0[i] * pa_x[i] - ta1_z_xxx_yyyyyy_1[i] * pc_x[i];

        ta1_z_xxxx_yyyyyz_0[i] =
            3.0 * ta1_z_xx_yyyyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyyyz_1[i] * fe_0 + ta1_z_xxx_yyyyyz_0[i] * pa_x[i] - ta1_z_xxx_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxxx_yyyyzz_0[i] =
            3.0 * ta1_z_xx_yyyyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyyzz_1[i] * fe_0 + ta1_z_xxx_yyyyzz_0[i] * pa_x[i] - ta1_z_xxx_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxxx_yyyzzz_0[i] =
            3.0 * ta1_z_xx_yyyzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyzzz_1[i] * fe_0 + ta1_z_xxx_yyyzzz_0[i] * pa_x[i] - ta1_z_xxx_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxxx_yyzzzz_0[i] =
            3.0 * ta1_z_xx_yyzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyzzzz_1[i] * fe_0 + ta1_z_xxx_yyzzzz_0[i] * pa_x[i] - ta1_z_xxx_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxxx_yzzzzz_0[i] =
            3.0 * ta1_z_xx_yzzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yzzzzz_1[i] * fe_0 + ta1_z_xxx_yzzzzz_0[i] * pa_x[i] - ta1_z_xxx_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxxx_zzzzzz_0[i] =
            3.0 * ta1_z_xx_zzzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_zzzzzz_1[i] * fe_0 + ta1_z_xxx_zzzzzz_0[i] * pa_x[i] - ta1_z_xxx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 868-896 components of targeted buffer : GI

    auto ta1_z_xxxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 868);

    auto ta1_z_xxxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 869);

    auto ta1_z_xxxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 870);

    auto ta1_z_xxxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 871);

    auto ta1_z_xxxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 872);

    auto ta1_z_xxxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 873);

    auto ta1_z_xxxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 874);

    auto ta1_z_xxxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 875);

    auto ta1_z_xxxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 876);

    auto ta1_z_xxxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 877);

    auto ta1_z_xxxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 878);

    auto ta1_z_xxxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 879);

    auto ta1_z_xxxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 880);

    auto ta1_z_xxxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 881);

    auto ta1_z_xxxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 882);

    auto ta1_z_xxxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 883);

    auto ta1_z_xxxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 884);

    auto ta1_z_xxxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 885);

    auto ta1_z_xxxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 886);

    auto ta1_z_xxxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 887);

    auto ta1_z_xxxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 888);

    auto ta1_z_xxxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 889);

    auto ta1_z_xxxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 890);

    auto ta1_z_xxxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 891);

    auto ta1_z_xxxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 892);

    auto ta1_z_xxxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 893);

    auto ta1_z_xxxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 894);

    auto ta1_z_xxxy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 895);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_z_xxx_xxxxx_0,   \
                             ta1_z_xxx_xxxxx_1,   \
                             ta1_z_xxx_xxxxxx_0,  \
                             ta1_z_xxx_xxxxxx_1,  \
                             ta1_z_xxx_xxxxxy_0,  \
                             ta1_z_xxx_xxxxxy_1,  \
                             ta1_z_xxx_xxxxxz_0,  \
                             ta1_z_xxx_xxxxxz_1,  \
                             ta1_z_xxx_xxxxy_0,   \
                             ta1_z_xxx_xxxxy_1,   \
                             ta1_z_xxx_xxxxyy_0,  \
                             ta1_z_xxx_xxxxyy_1,  \
                             ta1_z_xxx_xxxxyz_0,  \
                             ta1_z_xxx_xxxxyz_1,  \
                             ta1_z_xxx_xxxxz_0,   \
                             ta1_z_xxx_xxxxz_1,   \
                             ta1_z_xxx_xxxxzz_0,  \
                             ta1_z_xxx_xxxxzz_1,  \
                             ta1_z_xxx_xxxyy_0,   \
                             ta1_z_xxx_xxxyy_1,   \
                             ta1_z_xxx_xxxyyy_0,  \
                             ta1_z_xxx_xxxyyy_1,  \
                             ta1_z_xxx_xxxyyz_0,  \
                             ta1_z_xxx_xxxyyz_1,  \
                             ta1_z_xxx_xxxyz_0,   \
                             ta1_z_xxx_xxxyz_1,   \
                             ta1_z_xxx_xxxyzz_0,  \
                             ta1_z_xxx_xxxyzz_1,  \
                             ta1_z_xxx_xxxzz_0,   \
                             ta1_z_xxx_xxxzz_1,   \
                             ta1_z_xxx_xxxzzz_0,  \
                             ta1_z_xxx_xxxzzz_1,  \
                             ta1_z_xxx_xxyyy_0,   \
                             ta1_z_xxx_xxyyy_1,   \
                             ta1_z_xxx_xxyyyy_0,  \
                             ta1_z_xxx_xxyyyy_1,  \
                             ta1_z_xxx_xxyyyz_0,  \
                             ta1_z_xxx_xxyyyz_1,  \
                             ta1_z_xxx_xxyyz_0,   \
                             ta1_z_xxx_xxyyz_1,   \
                             ta1_z_xxx_xxyyzz_0,  \
                             ta1_z_xxx_xxyyzz_1,  \
                             ta1_z_xxx_xxyzz_0,   \
                             ta1_z_xxx_xxyzz_1,   \
                             ta1_z_xxx_xxyzzz_0,  \
                             ta1_z_xxx_xxyzzz_1,  \
                             ta1_z_xxx_xxzzz_0,   \
                             ta1_z_xxx_xxzzz_1,   \
                             ta1_z_xxx_xxzzzz_0,  \
                             ta1_z_xxx_xxzzzz_1,  \
                             ta1_z_xxx_xyyyy_0,   \
                             ta1_z_xxx_xyyyy_1,   \
                             ta1_z_xxx_xyyyyy_0,  \
                             ta1_z_xxx_xyyyyy_1,  \
                             ta1_z_xxx_xyyyyz_0,  \
                             ta1_z_xxx_xyyyyz_1,  \
                             ta1_z_xxx_xyyyz_0,   \
                             ta1_z_xxx_xyyyz_1,   \
                             ta1_z_xxx_xyyyzz_0,  \
                             ta1_z_xxx_xyyyzz_1,  \
                             ta1_z_xxx_xyyzz_0,   \
                             ta1_z_xxx_xyyzz_1,   \
                             ta1_z_xxx_xyyzzz_0,  \
                             ta1_z_xxx_xyyzzz_1,  \
                             ta1_z_xxx_xyzzz_0,   \
                             ta1_z_xxx_xyzzz_1,   \
                             ta1_z_xxx_xyzzzz_0,  \
                             ta1_z_xxx_xyzzzz_1,  \
                             ta1_z_xxx_xzzzz_0,   \
                             ta1_z_xxx_xzzzz_1,   \
                             ta1_z_xxx_xzzzzz_0,  \
                             ta1_z_xxx_xzzzzz_1,  \
                             ta1_z_xxx_zzzzzz_0,  \
                             ta1_z_xxx_zzzzzz_1,  \
                             ta1_z_xxxy_xxxxxx_0, \
                             ta1_z_xxxy_xxxxxy_0, \
                             ta1_z_xxxy_xxxxxz_0, \
                             ta1_z_xxxy_xxxxyy_0, \
                             ta1_z_xxxy_xxxxyz_0, \
                             ta1_z_xxxy_xxxxzz_0, \
                             ta1_z_xxxy_xxxyyy_0, \
                             ta1_z_xxxy_xxxyyz_0, \
                             ta1_z_xxxy_xxxyzz_0, \
                             ta1_z_xxxy_xxxzzz_0, \
                             ta1_z_xxxy_xxyyyy_0, \
                             ta1_z_xxxy_xxyyyz_0, \
                             ta1_z_xxxy_xxyyzz_0, \
                             ta1_z_xxxy_xxyzzz_0, \
                             ta1_z_xxxy_xxzzzz_0, \
                             ta1_z_xxxy_xyyyyy_0, \
                             ta1_z_xxxy_xyyyyz_0, \
                             ta1_z_xxxy_xyyyzz_0, \
                             ta1_z_xxxy_xyyzzz_0, \
                             ta1_z_xxxy_xyzzzz_0, \
                             ta1_z_xxxy_xzzzzz_0, \
                             ta1_z_xxxy_yyyyyy_0, \
                             ta1_z_xxxy_yyyyyz_0, \
                             ta1_z_xxxy_yyyyzz_0, \
                             ta1_z_xxxy_yyyzzz_0, \
                             ta1_z_xxxy_yyzzzz_0, \
                             ta1_z_xxxy_yzzzzz_0, \
                             ta1_z_xxxy_zzzzzz_0, \
                             ta1_z_xxy_yyyyyy_0,  \
                             ta1_z_xxy_yyyyyy_1,  \
                             ta1_z_xxy_yyyyyz_0,  \
                             ta1_z_xxy_yyyyyz_1,  \
                             ta1_z_xxy_yyyyzz_0,  \
                             ta1_z_xxy_yyyyzz_1,  \
                             ta1_z_xxy_yyyzzz_0,  \
                             ta1_z_xxy_yyyzzz_1,  \
                             ta1_z_xxy_yyzzzz_0,  \
                             ta1_z_xxy_yyzzzz_1,  \
                             ta1_z_xxy_yzzzzz_0,  \
                             ta1_z_xxy_yzzzzz_1,  \
                             ta1_z_xy_yyyyyy_0,   \
                             ta1_z_xy_yyyyyy_1,   \
                             ta1_z_xy_yyyyyz_0,   \
                             ta1_z_xy_yyyyyz_1,   \
                             ta1_z_xy_yyyyzz_0,   \
                             ta1_z_xy_yyyyzz_1,   \
                             ta1_z_xy_yyyzzz_0,   \
                             ta1_z_xy_yyyzzz_1,   \
                             ta1_z_xy_yyzzzz_0,   \
                             ta1_z_xy_yyzzzz_1,   \
                             ta1_z_xy_yzzzzz_0,   \
                             ta1_z_xy_yzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxy_xxxxxx_0[i] = ta1_z_xxx_xxxxxx_0[i] * pa_y[i] - ta1_z_xxx_xxxxxx_1[i] * pc_y[i];

        ta1_z_xxxy_xxxxxy_0[i] =
            ta1_z_xxx_xxxxx_0[i] * fe_0 - ta1_z_xxx_xxxxx_1[i] * fe_0 + ta1_z_xxx_xxxxxy_0[i] * pa_y[i] - ta1_z_xxx_xxxxxy_1[i] * pc_y[i];

        ta1_z_xxxy_xxxxxz_0[i] = ta1_z_xxx_xxxxxz_0[i] * pa_y[i] - ta1_z_xxx_xxxxxz_1[i] * pc_y[i];

        ta1_z_xxxy_xxxxyy_0[i] =
            2.0 * ta1_z_xxx_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxxxy_1[i] * fe_0 + ta1_z_xxx_xxxxyy_0[i] * pa_y[i] - ta1_z_xxx_xxxxyy_1[i] * pc_y[i];

        ta1_z_xxxy_xxxxyz_0[i] =
            ta1_z_xxx_xxxxz_0[i] * fe_0 - ta1_z_xxx_xxxxz_1[i] * fe_0 + ta1_z_xxx_xxxxyz_0[i] * pa_y[i] - ta1_z_xxx_xxxxyz_1[i] * pc_y[i];

        ta1_z_xxxy_xxxxzz_0[i] = ta1_z_xxx_xxxxzz_0[i] * pa_y[i] - ta1_z_xxx_xxxxzz_1[i] * pc_y[i];

        ta1_z_xxxy_xxxyyy_0[i] =
            3.0 * ta1_z_xxx_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_xxx_xxxyy_1[i] * fe_0 + ta1_z_xxx_xxxyyy_0[i] * pa_y[i] - ta1_z_xxx_xxxyyy_1[i] * pc_y[i];

        ta1_z_xxxy_xxxyyz_0[i] =
            2.0 * ta1_z_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxxyz_1[i] * fe_0 + ta1_z_xxx_xxxyyz_0[i] * pa_y[i] - ta1_z_xxx_xxxyyz_1[i] * pc_y[i];

        ta1_z_xxxy_xxxyzz_0[i] =
            ta1_z_xxx_xxxzz_0[i] * fe_0 - ta1_z_xxx_xxxzz_1[i] * fe_0 + ta1_z_xxx_xxxyzz_0[i] * pa_y[i] - ta1_z_xxx_xxxyzz_1[i] * pc_y[i];

        ta1_z_xxxy_xxxzzz_0[i] = ta1_z_xxx_xxxzzz_0[i] * pa_y[i] - ta1_z_xxx_xxxzzz_1[i] * pc_y[i];

        ta1_z_xxxy_xxyyyy_0[i] =
            4.0 * ta1_z_xxx_xxyyy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxyyy_1[i] * fe_0 + ta1_z_xxx_xxyyyy_0[i] * pa_y[i] - ta1_z_xxx_xxyyyy_1[i] * pc_y[i];

        ta1_z_xxxy_xxyyyz_0[i] =
            3.0 * ta1_z_xxx_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xxyyz_1[i] * fe_0 + ta1_z_xxx_xxyyyz_0[i] * pa_y[i] - ta1_z_xxx_xxyyyz_1[i] * pc_y[i];

        ta1_z_xxxy_xxyyzz_0[i] =
            2.0 * ta1_z_xxx_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxyzz_1[i] * fe_0 + ta1_z_xxx_xxyyzz_0[i] * pa_y[i] - ta1_z_xxx_xxyyzz_1[i] * pc_y[i];

        ta1_z_xxxy_xxyzzz_0[i] =
            ta1_z_xxx_xxzzz_0[i] * fe_0 - ta1_z_xxx_xxzzz_1[i] * fe_0 + ta1_z_xxx_xxyzzz_0[i] * pa_y[i] - ta1_z_xxx_xxyzzz_1[i] * pc_y[i];

        ta1_z_xxxy_xxzzzz_0[i] = ta1_z_xxx_xxzzzz_0[i] * pa_y[i] - ta1_z_xxx_xxzzzz_1[i] * pc_y[i];

        ta1_z_xxxy_xyyyyy_0[i] =
            5.0 * ta1_z_xxx_xyyyy_0[i] * fe_0 - 5.0 * ta1_z_xxx_xyyyy_1[i] * fe_0 + ta1_z_xxx_xyyyyy_0[i] * pa_y[i] - ta1_z_xxx_xyyyyy_1[i] * pc_y[i];

        ta1_z_xxxy_xyyyyz_0[i] =
            4.0 * ta1_z_xxx_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyyyz_1[i] * fe_0 + ta1_z_xxx_xyyyyz_0[i] * pa_y[i] - ta1_z_xxx_xyyyyz_1[i] * pc_y[i];

        ta1_z_xxxy_xyyyzz_0[i] =
            3.0 * ta1_z_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xyyzz_1[i] * fe_0 + ta1_z_xxx_xyyyzz_0[i] * pa_y[i] - ta1_z_xxx_xyyyzz_1[i] * pc_y[i];

        ta1_z_xxxy_xyyzzz_0[i] =
            2.0 * ta1_z_xxx_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xyzzz_1[i] * fe_0 + ta1_z_xxx_xyyzzz_0[i] * pa_y[i] - ta1_z_xxx_xyyzzz_1[i] * pc_y[i];

        ta1_z_xxxy_xyzzzz_0[i] =
            ta1_z_xxx_xzzzz_0[i] * fe_0 - ta1_z_xxx_xzzzz_1[i] * fe_0 + ta1_z_xxx_xyzzzz_0[i] * pa_y[i] - ta1_z_xxx_xyzzzz_1[i] * pc_y[i];

        ta1_z_xxxy_xzzzzz_0[i] = ta1_z_xxx_xzzzzz_0[i] * pa_y[i] - ta1_z_xxx_xzzzzz_1[i] * pc_y[i];

        ta1_z_xxxy_yyyyyy_0[i] =
            2.0 * ta1_z_xy_yyyyyy_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyyyy_1[i] * fe_0 + ta1_z_xxy_yyyyyy_0[i] * pa_x[i] - ta1_z_xxy_yyyyyy_1[i] * pc_x[i];

        ta1_z_xxxy_yyyyyz_0[i] =
            2.0 * ta1_z_xy_yyyyyz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyyyz_1[i] * fe_0 + ta1_z_xxy_yyyyyz_0[i] * pa_x[i] - ta1_z_xxy_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxxy_yyyyzz_0[i] =
            2.0 * ta1_z_xy_yyyyzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyyzz_1[i] * fe_0 + ta1_z_xxy_yyyyzz_0[i] * pa_x[i] - ta1_z_xxy_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxxy_yyyzzz_0[i] =
            2.0 * ta1_z_xy_yyyzzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyzzz_1[i] * fe_0 + ta1_z_xxy_yyyzzz_0[i] * pa_x[i] - ta1_z_xxy_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxxy_yyzzzz_0[i] =
            2.0 * ta1_z_xy_yyzzzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyzzzz_1[i] * fe_0 + ta1_z_xxy_yyzzzz_0[i] * pa_x[i] - ta1_z_xxy_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxxy_yzzzzz_0[i] =
            2.0 * ta1_z_xy_yzzzzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yzzzzz_1[i] * fe_0 + ta1_z_xxy_yzzzzz_0[i] * pa_x[i] - ta1_z_xxy_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxxy_zzzzzz_0[i] = ta1_z_xxx_zzzzzz_0[i] * pa_y[i] - ta1_z_xxx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 896-924 components of targeted buffer : GI

    auto ta1_z_xxxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 896);

    auto ta1_z_xxxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 897);

    auto ta1_z_xxxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 898);

    auto ta1_z_xxxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 899);

    auto ta1_z_xxxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 900);

    auto ta1_z_xxxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 901);

    auto ta1_z_xxxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 902);

    auto ta1_z_xxxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 903);

    auto ta1_z_xxxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 904);

    auto ta1_z_xxxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 905);

    auto ta1_z_xxxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 906);

    auto ta1_z_xxxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 907);

    auto ta1_z_xxxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 908);

    auto ta1_z_xxxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 909);

    auto ta1_z_xxxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 910);

    auto ta1_z_xxxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 911);

    auto ta1_z_xxxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 912);

    auto ta1_z_xxxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 913);

    auto ta1_z_xxxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 914);

    auto ta1_z_xxxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 915);

    auto ta1_z_xxxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 916);

    auto ta1_z_xxxz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 917);

    auto ta1_z_xxxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 918);

    auto ta1_z_xxxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 919);

    auto ta1_z_xxxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 920);

    auto ta1_z_xxxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 921);

    auto ta1_z_xxxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 922);

    auto ta1_z_xxxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 923);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_z_xxx_xxxxx_0,   \
                             ta1_z_xxx_xxxxx_1,   \
                             ta1_z_xxx_xxxxxx_0,  \
                             ta1_z_xxx_xxxxxx_1,  \
                             ta1_z_xxx_xxxxxy_0,  \
                             ta1_z_xxx_xxxxxy_1,  \
                             ta1_z_xxx_xxxxxz_0,  \
                             ta1_z_xxx_xxxxxz_1,  \
                             ta1_z_xxx_xxxxy_0,   \
                             ta1_z_xxx_xxxxy_1,   \
                             ta1_z_xxx_xxxxyy_0,  \
                             ta1_z_xxx_xxxxyy_1,  \
                             ta1_z_xxx_xxxxyz_0,  \
                             ta1_z_xxx_xxxxyz_1,  \
                             ta1_z_xxx_xxxxz_0,   \
                             ta1_z_xxx_xxxxz_1,   \
                             ta1_z_xxx_xxxxzz_0,  \
                             ta1_z_xxx_xxxxzz_1,  \
                             ta1_z_xxx_xxxyy_0,   \
                             ta1_z_xxx_xxxyy_1,   \
                             ta1_z_xxx_xxxyyy_0,  \
                             ta1_z_xxx_xxxyyy_1,  \
                             ta1_z_xxx_xxxyyz_0,  \
                             ta1_z_xxx_xxxyyz_1,  \
                             ta1_z_xxx_xxxyz_0,   \
                             ta1_z_xxx_xxxyz_1,   \
                             ta1_z_xxx_xxxyzz_0,  \
                             ta1_z_xxx_xxxyzz_1,  \
                             ta1_z_xxx_xxxzz_0,   \
                             ta1_z_xxx_xxxzz_1,   \
                             ta1_z_xxx_xxxzzz_0,  \
                             ta1_z_xxx_xxxzzz_1,  \
                             ta1_z_xxx_xxyyy_0,   \
                             ta1_z_xxx_xxyyy_1,   \
                             ta1_z_xxx_xxyyyy_0,  \
                             ta1_z_xxx_xxyyyy_1,  \
                             ta1_z_xxx_xxyyyz_0,  \
                             ta1_z_xxx_xxyyyz_1,  \
                             ta1_z_xxx_xxyyz_0,   \
                             ta1_z_xxx_xxyyz_1,   \
                             ta1_z_xxx_xxyyzz_0,  \
                             ta1_z_xxx_xxyyzz_1,  \
                             ta1_z_xxx_xxyzz_0,   \
                             ta1_z_xxx_xxyzz_1,   \
                             ta1_z_xxx_xxyzzz_0,  \
                             ta1_z_xxx_xxyzzz_1,  \
                             ta1_z_xxx_xxzzz_0,   \
                             ta1_z_xxx_xxzzz_1,   \
                             ta1_z_xxx_xxzzzz_0,  \
                             ta1_z_xxx_xxzzzz_1,  \
                             ta1_z_xxx_xyyyy_0,   \
                             ta1_z_xxx_xyyyy_1,   \
                             ta1_z_xxx_xyyyyy_0,  \
                             ta1_z_xxx_xyyyyy_1,  \
                             ta1_z_xxx_xyyyyz_0,  \
                             ta1_z_xxx_xyyyyz_1,  \
                             ta1_z_xxx_xyyyz_0,   \
                             ta1_z_xxx_xyyyz_1,   \
                             ta1_z_xxx_xyyyzz_0,  \
                             ta1_z_xxx_xyyyzz_1,  \
                             ta1_z_xxx_xyyzz_0,   \
                             ta1_z_xxx_xyyzz_1,   \
                             ta1_z_xxx_xyyzzz_0,  \
                             ta1_z_xxx_xyyzzz_1,  \
                             ta1_z_xxx_xyzzz_0,   \
                             ta1_z_xxx_xyzzz_1,   \
                             ta1_z_xxx_xyzzzz_0,  \
                             ta1_z_xxx_xyzzzz_1,  \
                             ta1_z_xxx_xzzzz_0,   \
                             ta1_z_xxx_xzzzz_1,   \
                             ta1_z_xxx_xzzzzz_0,  \
                             ta1_z_xxx_xzzzzz_1,  \
                             ta1_z_xxx_yyyyyy_0,  \
                             ta1_z_xxx_yyyyyy_1,  \
                             ta1_z_xxxz_xxxxxx_0, \
                             ta1_z_xxxz_xxxxxy_0, \
                             ta1_z_xxxz_xxxxxz_0, \
                             ta1_z_xxxz_xxxxyy_0, \
                             ta1_z_xxxz_xxxxyz_0, \
                             ta1_z_xxxz_xxxxzz_0, \
                             ta1_z_xxxz_xxxyyy_0, \
                             ta1_z_xxxz_xxxyyz_0, \
                             ta1_z_xxxz_xxxyzz_0, \
                             ta1_z_xxxz_xxxzzz_0, \
                             ta1_z_xxxz_xxyyyy_0, \
                             ta1_z_xxxz_xxyyyz_0, \
                             ta1_z_xxxz_xxyyzz_0, \
                             ta1_z_xxxz_xxyzzz_0, \
                             ta1_z_xxxz_xxzzzz_0, \
                             ta1_z_xxxz_xyyyyy_0, \
                             ta1_z_xxxz_xyyyyz_0, \
                             ta1_z_xxxz_xyyyzz_0, \
                             ta1_z_xxxz_xyyzzz_0, \
                             ta1_z_xxxz_xyzzzz_0, \
                             ta1_z_xxxz_xzzzzz_0, \
                             ta1_z_xxxz_yyyyyy_0, \
                             ta1_z_xxxz_yyyyyz_0, \
                             ta1_z_xxxz_yyyyzz_0, \
                             ta1_z_xxxz_yyyzzz_0, \
                             ta1_z_xxxz_yyzzzz_0, \
                             ta1_z_xxxz_yzzzzz_0, \
                             ta1_z_xxxz_zzzzzz_0, \
                             ta1_z_xxz_yyyyyz_0,  \
                             ta1_z_xxz_yyyyyz_1,  \
                             ta1_z_xxz_yyyyzz_0,  \
                             ta1_z_xxz_yyyyzz_1,  \
                             ta1_z_xxz_yyyzzz_0,  \
                             ta1_z_xxz_yyyzzz_1,  \
                             ta1_z_xxz_yyzzzz_0,  \
                             ta1_z_xxz_yyzzzz_1,  \
                             ta1_z_xxz_yzzzzz_0,  \
                             ta1_z_xxz_yzzzzz_1,  \
                             ta1_z_xxz_zzzzzz_0,  \
                             ta1_z_xxz_zzzzzz_1,  \
                             ta1_z_xz_yyyyyz_0,   \
                             ta1_z_xz_yyyyyz_1,   \
                             ta1_z_xz_yyyyzz_0,   \
                             ta1_z_xz_yyyyzz_1,   \
                             ta1_z_xz_yyyzzz_0,   \
                             ta1_z_xz_yyyzzz_1,   \
                             ta1_z_xz_yyzzzz_0,   \
                             ta1_z_xz_yyzzzz_1,   \
                             ta1_z_xz_yzzzzz_0,   \
                             ta1_z_xz_yzzzzz_1,   \
                             ta1_z_xz_zzzzzz_0,   \
                             ta1_z_xz_zzzzzz_1,   \
                             ta_xxx_xxxxxx_1,     \
                             ta_xxx_xxxxxy_1,     \
                             ta_xxx_xxxxxz_1,     \
                             ta_xxx_xxxxyy_1,     \
                             ta_xxx_xxxxyz_1,     \
                             ta_xxx_xxxxzz_1,     \
                             ta_xxx_xxxyyy_1,     \
                             ta_xxx_xxxyyz_1,     \
                             ta_xxx_xxxyzz_1,     \
                             ta_xxx_xxxzzz_1,     \
                             ta_xxx_xxyyyy_1,     \
                             ta_xxx_xxyyyz_1,     \
                             ta_xxx_xxyyzz_1,     \
                             ta_xxx_xxyzzz_1,     \
                             ta_xxx_xxzzzz_1,     \
                             ta_xxx_xyyyyy_1,     \
                             ta_xxx_xyyyyz_1,     \
                             ta_xxx_xyyyzz_1,     \
                             ta_xxx_xyyzzz_1,     \
                             ta_xxx_xyzzzz_1,     \
                             ta_xxx_xzzzzz_1,     \
                             ta_xxx_yyyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxz_xxxxxx_0[i] = ta_xxx_xxxxxx_1[i] + ta1_z_xxx_xxxxxx_0[i] * pa_z[i] - ta1_z_xxx_xxxxxx_1[i] * pc_z[i];

        ta1_z_xxxz_xxxxxy_0[i] = ta_xxx_xxxxxy_1[i] + ta1_z_xxx_xxxxxy_0[i] * pa_z[i] - ta1_z_xxx_xxxxxy_1[i] * pc_z[i];

        ta1_z_xxxz_xxxxxz_0[i] = ta1_z_xxx_xxxxx_0[i] * fe_0 - ta1_z_xxx_xxxxx_1[i] * fe_0 + ta_xxx_xxxxxz_1[i] + ta1_z_xxx_xxxxxz_0[i] * pa_z[i] -
                                 ta1_z_xxx_xxxxxz_1[i] * pc_z[i];

        ta1_z_xxxz_xxxxyy_0[i] = ta_xxx_xxxxyy_1[i] + ta1_z_xxx_xxxxyy_0[i] * pa_z[i] - ta1_z_xxx_xxxxyy_1[i] * pc_z[i];

        ta1_z_xxxz_xxxxyz_0[i] = ta1_z_xxx_xxxxy_0[i] * fe_0 - ta1_z_xxx_xxxxy_1[i] * fe_0 + ta_xxx_xxxxyz_1[i] + ta1_z_xxx_xxxxyz_0[i] * pa_z[i] -
                                 ta1_z_xxx_xxxxyz_1[i] * pc_z[i];

        ta1_z_xxxz_xxxxzz_0[i] = 2.0 * ta1_z_xxx_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxxxz_1[i] * fe_0 + ta_xxx_xxxxzz_1[i] +
                                 ta1_z_xxx_xxxxzz_0[i] * pa_z[i] - ta1_z_xxx_xxxxzz_1[i] * pc_z[i];

        ta1_z_xxxz_xxxyyy_0[i] = ta_xxx_xxxyyy_1[i] + ta1_z_xxx_xxxyyy_0[i] * pa_z[i] - ta1_z_xxx_xxxyyy_1[i] * pc_z[i];

        ta1_z_xxxz_xxxyyz_0[i] = ta1_z_xxx_xxxyy_0[i] * fe_0 - ta1_z_xxx_xxxyy_1[i] * fe_0 + ta_xxx_xxxyyz_1[i] + ta1_z_xxx_xxxyyz_0[i] * pa_z[i] -
                                 ta1_z_xxx_xxxyyz_1[i] * pc_z[i];

        ta1_z_xxxz_xxxyzz_0[i] = 2.0 * ta1_z_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxxyz_1[i] * fe_0 + ta_xxx_xxxyzz_1[i] +
                                 ta1_z_xxx_xxxyzz_0[i] * pa_z[i] - ta1_z_xxx_xxxyzz_1[i] * pc_z[i];

        ta1_z_xxxz_xxxzzz_0[i] = 3.0 * ta1_z_xxx_xxxzz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xxxzz_1[i] * fe_0 + ta_xxx_xxxzzz_1[i] +
                                 ta1_z_xxx_xxxzzz_0[i] * pa_z[i] - ta1_z_xxx_xxxzzz_1[i] * pc_z[i];

        ta1_z_xxxz_xxyyyy_0[i] = ta_xxx_xxyyyy_1[i] + ta1_z_xxx_xxyyyy_0[i] * pa_z[i] - ta1_z_xxx_xxyyyy_1[i] * pc_z[i];

        ta1_z_xxxz_xxyyyz_0[i] = ta1_z_xxx_xxyyy_0[i] * fe_0 - ta1_z_xxx_xxyyy_1[i] * fe_0 + ta_xxx_xxyyyz_1[i] + ta1_z_xxx_xxyyyz_0[i] * pa_z[i] -
                                 ta1_z_xxx_xxyyyz_1[i] * pc_z[i];

        ta1_z_xxxz_xxyyzz_0[i] = 2.0 * ta1_z_xxx_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxyyz_1[i] * fe_0 + ta_xxx_xxyyzz_1[i] +
                                 ta1_z_xxx_xxyyzz_0[i] * pa_z[i] - ta1_z_xxx_xxyyzz_1[i] * pc_z[i];

        ta1_z_xxxz_xxyzzz_0[i] = 3.0 * ta1_z_xxx_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xxyzz_1[i] * fe_0 + ta_xxx_xxyzzz_1[i] +
                                 ta1_z_xxx_xxyzzz_0[i] * pa_z[i] - ta1_z_xxx_xxyzzz_1[i] * pc_z[i];

        ta1_z_xxxz_xxzzzz_0[i] = 4.0 * ta1_z_xxx_xxzzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xxzzz_1[i] * fe_0 + ta_xxx_xxzzzz_1[i] +
                                 ta1_z_xxx_xxzzzz_0[i] * pa_z[i] - ta1_z_xxx_xxzzzz_1[i] * pc_z[i];

        ta1_z_xxxz_xyyyyy_0[i] = ta_xxx_xyyyyy_1[i] + ta1_z_xxx_xyyyyy_0[i] * pa_z[i] - ta1_z_xxx_xyyyyy_1[i] * pc_z[i];

        ta1_z_xxxz_xyyyyz_0[i] = ta1_z_xxx_xyyyy_0[i] * fe_0 - ta1_z_xxx_xyyyy_1[i] * fe_0 + ta_xxx_xyyyyz_1[i] + ta1_z_xxx_xyyyyz_0[i] * pa_z[i] -
                                 ta1_z_xxx_xyyyyz_1[i] * pc_z[i];

        ta1_z_xxxz_xyyyzz_0[i] = 2.0 * ta1_z_xxx_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xyyyz_1[i] * fe_0 + ta_xxx_xyyyzz_1[i] +
                                 ta1_z_xxx_xyyyzz_0[i] * pa_z[i] - ta1_z_xxx_xyyyzz_1[i] * pc_z[i];

        ta1_z_xxxz_xyyzzz_0[i] = 3.0 * ta1_z_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xyyzz_1[i] * fe_0 + ta_xxx_xyyzzz_1[i] +
                                 ta1_z_xxx_xyyzzz_0[i] * pa_z[i] - ta1_z_xxx_xyyzzz_1[i] * pc_z[i];

        ta1_z_xxxz_xyzzzz_0[i] = 4.0 * ta1_z_xxx_xyzzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyzzz_1[i] * fe_0 + ta_xxx_xyzzzz_1[i] +
                                 ta1_z_xxx_xyzzzz_0[i] * pa_z[i] - ta1_z_xxx_xyzzzz_1[i] * pc_z[i];

        ta1_z_xxxz_xzzzzz_0[i] = 5.0 * ta1_z_xxx_xzzzz_0[i] * fe_0 - 5.0 * ta1_z_xxx_xzzzz_1[i] * fe_0 + ta_xxx_xzzzzz_1[i] +
                                 ta1_z_xxx_xzzzzz_0[i] * pa_z[i] - ta1_z_xxx_xzzzzz_1[i] * pc_z[i];

        ta1_z_xxxz_yyyyyy_0[i] = ta_xxx_yyyyyy_1[i] + ta1_z_xxx_yyyyyy_0[i] * pa_z[i] - ta1_z_xxx_yyyyyy_1[i] * pc_z[i];

        ta1_z_xxxz_yyyyyz_0[i] =
            2.0 * ta1_z_xz_yyyyyz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyyyyz_1[i] * fe_0 + ta1_z_xxz_yyyyyz_0[i] * pa_x[i] - ta1_z_xxz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxxz_yyyyzz_0[i] =
            2.0 * ta1_z_xz_yyyyzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyyyzz_1[i] * fe_0 + ta1_z_xxz_yyyyzz_0[i] * pa_x[i] - ta1_z_xxz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxxz_yyyzzz_0[i] =
            2.0 * ta1_z_xz_yyyzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyyzzz_1[i] * fe_0 + ta1_z_xxz_yyyzzz_0[i] * pa_x[i] - ta1_z_xxz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxxz_yyzzzz_0[i] =
            2.0 * ta1_z_xz_yyzzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyzzzz_1[i] * fe_0 + ta1_z_xxz_yyzzzz_0[i] * pa_x[i] - ta1_z_xxz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxxz_yzzzzz_0[i] =
            2.0 * ta1_z_xz_yzzzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yzzzzz_1[i] * fe_0 + ta1_z_xxz_yzzzzz_0[i] * pa_x[i] - ta1_z_xxz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxxz_zzzzzz_0[i] =
            2.0 * ta1_z_xz_zzzzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_zzzzzz_1[i] * fe_0 + ta1_z_xxz_zzzzzz_0[i] * pa_x[i] - ta1_z_xxz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 924-952 components of targeted buffer : GI

    auto ta1_z_xxyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 924);

    auto ta1_z_xxyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 925);

    auto ta1_z_xxyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 926);

    auto ta1_z_xxyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 927);

    auto ta1_z_xxyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 928);

    auto ta1_z_xxyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 929);

    auto ta1_z_xxyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 930);

    auto ta1_z_xxyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 931);

    auto ta1_z_xxyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 932);

    auto ta1_z_xxyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 933);

    auto ta1_z_xxyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 934);

    auto ta1_z_xxyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 935);

    auto ta1_z_xxyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 936);

    auto ta1_z_xxyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 937);

    auto ta1_z_xxyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 938);

    auto ta1_z_xxyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 939);

    auto ta1_z_xxyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 940);

    auto ta1_z_xxyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 941);

    auto ta1_z_xxyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 942);

    auto ta1_z_xxyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 943);

    auto ta1_z_xxyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 944);

    auto ta1_z_xxyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 945);

    auto ta1_z_xxyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 946);

    auto ta1_z_xxyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 947);

    auto ta1_z_xxyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 948);

    auto ta1_z_xxyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 949);

    auto ta1_z_xxyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 950);

    auto ta1_z_xxyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 951);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_z_xx_xxxxxx_0,   \
                             ta1_z_xx_xxxxxx_1,   \
                             ta1_z_xx_xxxxxz_0,   \
                             ta1_z_xx_xxxxxz_1,   \
                             ta1_z_xx_xxxxzz_0,   \
                             ta1_z_xx_xxxxzz_1,   \
                             ta1_z_xx_xxxzzz_0,   \
                             ta1_z_xx_xxxzzz_1,   \
                             ta1_z_xx_xxzzzz_0,   \
                             ta1_z_xx_xxzzzz_1,   \
                             ta1_z_xx_xzzzzz_0,   \
                             ta1_z_xx_xzzzzz_1,   \
                             ta1_z_xxy_xxxxxx_0,  \
                             ta1_z_xxy_xxxxxx_1,  \
                             ta1_z_xxy_xxxxxz_0,  \
                             ta1_z_xxy_xxxxxz_1,  \
                             ta1_z_xxy_xxxxzz_0,  \
                             ta1_z_xxy_xxxxzz_1,  \
                             ta1_z_xxy_xxxzzz_0,  \
                             ta1_z_xxy_xxxzzz_1,  \
                             ta1_z_xxy_xxzzzz_0,  \
                             ta1_z_xxy_xxzzzz_1,  \
                             ta1_z_xxy_xzzzzz_0,  \
                             ta1_z_xxy_xzzzzz_1,  \
                             ta1_z_xxyy_xxxxxx_0, \
                             ta1_z_xxyy_xxxxxy_0, \
                             ta1_z_xxyy_xxxxxz_0, \
                             ta1_z_xxyy_xxxxyy_0, \
                             ta1_z_xxyy_xxxxyz_0, \
                             ta1_z_xxyy_xxxxzz_0, \
                             ta1_z_xxyy_xxxyyy_0, \
                             ta1_z_xxyy_xxxyyz_0, \
                             ta1_z_xxyy_xxxyzz_0, \
                             ta1_z_xxyy_xxxzzz_0, \
                             ta1_z_xxyy_xxyyyy_0, \
                             ta1_z_xxyy_xxyyyz_0, \
                             ta1_z_xxyy_xxyyzz_0, \
                             ta1_z_xxyy_xxyzzz_0, \
                             ta1_z_xxyy_xxzzzz_0, \
                             ta1_z_xxyy_xyyyyy_0, \
                             ta1_z_xxyy_xyyyyz_0, \
                             ta1_z_xxyy_xyyyzz_0, \
                             ta1_z_xxyy_xyyzzz_0, \
                             ta1_z_xxyy_xyzzzz_0, \
                             ta1_z_xxyy_xzzzzz_0, \
                             ta1_z_xxyy_yyyyyy_0, \
                             ta1_z_xxyy_yyyyyz_0, \
                             ta1_z_xxyy_yyyyzz_0, \
                             ta1_z_xxyy_yyyzzz_0, \
                             ta1_z_xxyy_yyzzzz_0, \
                             ta1_z_xxyy_yzzzzz_0, \
                             ta1_z_xxyy_zzzzzz_0, \
                             ta1_z_xyy_xxxxxy_0,  \
                             ta1_z_xyy_xxxxxy_1,  \
                             ta1_z_xyy_xxxxy_0,   \
                             ta1_z_xyy_xxxxy_1,   \
                             ta1_z_xyy_xxxxyy_0,  \
                             ta1_z_xyy_xxxxyy_1,  \
                             ta1_z_xyy_xxxxyz_0,  \
                             ta1_z_xyy_xxxxyz_1,  \
                             ta1_z_xyy_xxxyy_0,   \
                             ta1_z_xyy_xxxyy_1,   \
                             ta1_z_xyy_xxxyyy_0,  \
                             ta1_z_xyy_xxxyyy_1,  \
                             ta1_z_xyy_xxxyyz_0,  \
                             ta1_z_xyy_xxxyyz_1,  \
                             ta1_z_xyy_xxxyz_0,   \
                             ta1_z_xyy_xxxyz_1,   \
                             ta1_z_xyy_xxxyzz_0,  \
                             ta1_z_xyy_xxxyzz_1,  \
                             ta1_z_xyy_xxyyy_0,   \
                             ta1_z_xyy_xxyyy_1,   \
                             ta1_z_xyy_xxyyyy_0,  \
                             ta1_z_xyy_xxyyyy_1,  \
                             ta1_z_xyy_xxyyyz_0,  \
                             ta1_z_xyy_xxyyyz_1,  \
                             ta1_z_xyy_xxyyz_0,   \
                             ta1_z_xyy_xxyyz_1,   \
                             ta1_z_xyy_xxyyzz_0,  \
                             ta1_z_xyy_xxyyzz_1,  \
                             ta1_z_xyy_xxyzz_0,   \
                             ta1_z_xyy_xxyzz_1,   \
                             ta1_z_xyy_xxyzzz_0,  \
                             ta1_z_xyy_xxyzzz_1,  \
                             ta1_z_xyy_xyyyy_0,   \
                             ta1_z_xyy_xyyyy_1,   \
                             ta1_z_xyy_xyyyyy_0,  \
                             ta1_z_xyy_xyyyyy_1,  \
                             ta1_z_xyy_xyyyyz_0,  \
                             ta1_z_xyy_xyyyyz_1,  \
                             ta1_z_xyy_xyyyz_0,   \
                             ta1_z_xyy_xyyyz_1,   \
                             ta1_z_xyy_xyyyzz_0,  \
                             ta1_z_xyy_xyyyzz_1,  \
                             ta1_z_xyy_xyyzz_0,   \
                             ta1_z_xyy_xyyzz_1,   \
                             ta1_z_xyy_xyyzzz_0,  \
                             ta1_z_xyy_xyyzzz_1,  \
                             ta1_z_xyy_xyzzz_0,   \
                             ta1_z_xyy_xyzzz_1,   \
                             ta1_z_xyy_xyzzzz_0,  \
                             ta1_z_xyy_xyzzzz_1,  \
                             ta1_z_xyy_yyyyy_0,   \
                             ta1_z_xyy_yyyyy_1,   \
                             ta1_z_xyy_yyyyyy_0,  \
                             ta1_z_xyy_yyyyyy_1,  \
                             ta1_z_xyy_yyyyyz_0,  \
                             ta1_z_xyy_yyyyyz_1,  \
                             ta1_z_xyy_yyyyz_0,   \
                             ta1_z_xyy_yyyyz_1,   \
                             ta1_z_xyy_yyyyzz_0,  \
                             ta1_z_xyy_yyyyzz_1,  \
                             ta1_z_xyy_yyyzz_0,   \
                             ta1_z_xyy_yyyzz_1,   \
                             ta1_z_xyy_yyyzzz_0,  \
                             ta1_z_xyy_yyyzzz_1,  \
                             ta1_z_xyy_yyzzz_0,   \
                             ta1_z_xyy_yyzzz_1,   \
                             ta1_z_xyy_yyzzzz_0,  \
                             ta1_z_xyy_yyzzzz_1,  \
                             ta1_z_xyy_yzzzz_0,   \
                             ta1_z_xyy_yzzzz_1,   \
                             ta1_z_xyy_yzzzzz_0,  \
                             ta1_z_xyy_yzzzzz_1,  \
                             ta1_z_xyy_zzzzzz_0,  \
                             ta1_z_xyy_zzzzzz_1,  \
                             ta1_z_yy_xxxxxy_0,   \
                             ta1_z_yy_xxxxxy_1,   \
                             ta1_z_yy_xxxxyy_0,   \
                             ta1_z_yy_xxxxyy_1,   \
                             ta1_z_yy_xxxxyz_0,   \
                             ta1_z_yy_xxxxyz_1,   \
                             ta1_z_yy_xxxyyy_0,   \
                             ta1_z_yy_xxxyyy_1,   \
                             ta1_z_yy_xxxyyz_0,   \
                             ta1_z_yy_xxxyyz_1,   \
                             ta1_z_yy_xxxyzz_0,   \
                             ta1_z_yy_xxxyzz_1,   \
                             ta1_z_yy_xxyyyy_0,   \
                             ta1_z_yy_xxyyyy_1,   \
                             ta1_z_yy_xxyyyz_0,   \
                             ta1_z_yy_xxyyyz_1,   \
                             ta1_z_yy_xxyyzz_0,   \
                             ta1_z_yy_xxyyzz_1,   \
                             ta1_z_yy_xxyzzz_0,   \
                             ta1_z_yy_xxyzzz_1,   \
                             ta1_z_yy_xyyyyy_0,   \
                             ta1_z_yy_xyyyyy_1,   \
                             ta1_z_yy_xyyyyz_0,   \
                             ta1_z_yy_xyyyyz_1,   \
                             ta1_z_yy_xyyyzz_0,   \
                             ta1_z_yy_xyyyzz_1,   \
                             ta1_z_yy_xyyzzz_0,   \
                             ta1_z_yy_xyyzzz_1,   \
                             ta1_z_yy_xyzzzz_0,   \
                             ta1_z_yy_xyzzzz_1,   \
                             ta1_z_yy_yyyyyy_0,   \
                             ta1_z_yy_yyyyyy_1,   \
                             ta1_z_yy_yyyyyz_0,   \
                             ta1_z_yy_yyyyyz_1,   \
                             ta1_z_yy_yyyyzz_0,   \
                             ta1_z_yy_yyyyzz_1,   \
                             ta1_z_yy_yyyzzz_0,   \
                             ta1_z_yy_yyyzzz_1,   \
                             ta1_z_yy_yyzzzz_0,   \
                             ta1_z_yy_yyzzzz_1,   \
                             ta1_z_yy_yzzzzz_0,   \
                             ta1_z_yy_yzzzzz_1,   \
                             ta1_z_yy_zzzzzz_0,   \
                             ta1_z_yy_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyy_xxxxxx_0[i] =
            ta1_z_xx_xxxxxx_0[i] * fe_0 - ta1_z_xx_xxxxxx_1[i] * fe_0 + ta1_z_xxy_xxxxxx_0[i] * pa_y[i] - ta1_z_xxy_xxxxxx_1[i] * pc_y[i];

        ta1_z_xxyy_xxxxxy_0[i] = ta1_z_yy_xxxxxy_0[i] * fe_0 - ta1_z_yy_xxxxxy_1[i] * fe_0 + 5.0 * ta1_z_xyy_xxxxy_0[i] * fe_0 -
                                 5.0 * ta1_z_xyy_xxxxy_1[i] * fe_0 + ta1_z_xyy_xxxxxy_0[i] * pa_x[i] - ta1_z_xyy_xxxxxy_1[i] * pc_x[i];

        ta1_z_xxyy_xxxxxz_0[i] =
            ta1_z_xx_xxxxxz_0[i] * fe_0 - ta1_z_xx_xxxxxz_1[i] * fe_0 + ta1_z_xxy_xxxxxz_0[i] * pa_y[i] - ta1_z_xxy_xxxxxz_1[i] * pc_y[i];

        ta1_z_xxyy_xxxxyy_0[i] = ta1_z_yy_xxxxyy_0[i] * fe_0 - ta1_z_yy_xxxxyy_1[i] * fe_0 + 4.0 * ta1_z_xyy_xxxyy_0[i] * fe_0 -
                                 4.0 * ta1_z_xyy_xxxyy_1[i] * fe_0 + ta1_z_xyy_xxxxyy_0[i] * pa_x[i] - ta1_z_xyy_xxxxyy_1[i] * pc_x[i];

        ta1_z_xxyy_xxxxyz_0[i] = ta1_z_yy_xxxxyz_0[i] * fe_0 - ta1_z_yy_xxxxyz_1[i] * fe_0 + 4.0 * ta1_z_xyy_xxxyz_0[i] * fe_0 -
                                 4.0 * ta1_z_xyy_xxxyz_1[i] * fe_0 + ta1_z_xyy_xxxxyz_0[i] * pa_x[i] - ta1_z_xyy_xxxxyz_1[i] * pc_x[i];

        ta1_z_xxyy_xxxxzz_0[i] =
            ta1_z_xx_xxxxzz_0[i] * fe_0 - ta1_z_xx_xxxxzz_1[i] * fe_0 + ta1_z_xxy_xxxxzz_0[i] * pa_y[i] - ta1_z_xxy_xxxxzz_1[i] * pc_y[i];

        ta1_z_xxyy_xxxyyy_0[i] = ta1_z_yy_xxxyyy_0[i] * fe_0 - ta1_z_yy_xxxyyy_1[i] * fe_0 + 3.0 * ta1_z_xyy_xxyyy_0[i] * fe_0 -
                                 3.0 * ta1_z_xyy_xxyyy_1[i] * fe_0 + ta1_z_xyy_xxxyyy_0[i] * pa_x[i] - ta1_z_xyy_xxxyyy_1[i] * pc_x[i];

        ta1_z_xxyy_xxxyyz_0[i] = ta1_z_yy_xxxyyz_0[i] * fe_0 - ta1_z_yy_xxxyyz_1[i] * fe_0 + 3.0 * ta1_z_xyy_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_z_xyy_xxyyz_1[i] * fe_0 + ta1_z_xyy_xxxyyz_0[i] * pa_x[i] - ta1_z_xyy_xxxyyz_1[i] * pc_x[i];

        ta1_z_xxyy_xxxyzz_0[i] = ta1_z_yy_xxxyzz_0[i] * fe_0 - ta1_z_yy_xxxyzz_1[i] * fe_0 + 3.0 * ta1_z_xyy_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_xyy_xxyzz_1[i] * fe_0 + ta1_z_xyy_xxxyzz_0[i] * pa_x[i] - ta1_z_xyy_xxxyzz_1[i] * pc_x[i];

        ta1_z_xxyy_xxxzzz_0[i] =
            ta1_z_xx_xxxzzz_0[i] * fe_0 - ta1_z_xx_xxxzzz_1[i] * fe_0 + ta1_z_xxy_xxxzzz_0[i] * pa_y[i] - ta1_z_xxy_xxxzzz_1[i] * pc_y[i];

        ta1_z_xxyy_xxyyyy_0[i] = ta1_z_yy_xxyyyy_0[i] * fe_0 - ta1_z_yy_xxyyyy_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyyyy_0[i] * fe_0 -
                                 2.0 * ta1_z_xyy_xyyyy_1[i] * fe_0 + ta1_z_xyy_xxyyyy_0[i] * pa_x[i] - ta1_z_xyy_xxyyyy_1[i] * pc_x[i];

        ta1_z_xxyy_xxyyyz_0[i] = ta1_z_yy_xxyyyz_0[i] * fe_0 - ta1_z_yy_xxyyyz_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_z_xyy_xyyyz_1[i] * fe_0 + ta1_z_xyy_xxyyyz_0[i] * pa_x[i] - ta1_z_xyy_xxyyyz_1[i] * pc_x[i];

        ta1_z_xxyy_xxyyzz_0[i] = ta1_z_yy_xxyyzz_0[i] * fe_0 - ta1_z_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyyzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xyy_xyyzz_1[i] * fe_0 + ta1_z_xyy_xxyyzz_0[i] * pa_x[i] - ta1_z_xyy_xxyyzz_1[i] * pc_x[i];

        ta1_z_xxyy_xxyzzz_0[i] = ta1_z_yy_xxyzzz_0[i] * fe_0 - ta1_z_yy_xxyzzz_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xyy_xyzzz_1[i] * fe_0 + ta1_z_xyy_xxyzzz_0[i] * pa_x[i] - ta1_z_xyy_xxyzzz_1[i] * pc_x[i];

        ta1_z_xxyy_xxzzzz_0[i] =
            ta1_z_xx_xxzzzz_0[i] * fe_0 - ta1_z_xx_xxzzzz_1[i] * fe_0 + ta1_z_xxy_xxzzzz_0[i] * pa_y[i] - ta1_z_xxy_xxzzzz_1[i] * pc_y[i];

        ta1_z_xxyy_xyyyyy_0[i] = ta1_z_yy_xyyyyy_0[i] * fe_0 - ta1_z_yy_xyyyyy_1[i] * fe_0 + ta1_z_xyy_yyyyy_0[i] * fe_0 -
                                 ta1_z_xyy_yyyyy_1[i] * fe_0 + ta1_z_xyy_xyyyyy_0[i] * pa_x[i] - ta1_z_xyy_xyyyyy_1[i] * pc_x[i];

        ta1_z_xxyy_xyyyyz_0[i] = ta1_z_yy_xyyyyz_0[i] * fe_0 - ta1_z_yy_xyyyyz_1[i] * fe_0 + ta1_z_xyy_yyyyz_0[i] * fe_0 -
                                 ta1_z_xyy_yyyyz_1[i] * fe_0 + ta1_z_xyy_xyyyyz_0[i] * pa_x[i] - ta1_z_xyy_xyyyyz_1[i] * pc_x[i];

        ta1_z_xxyy_xyyyzz_0[i] = ta1_z_yy_xyyyzz_0[i] * fe_0 - ta1_z_yy_xyyyzz_1[i] * fe_0 + ta1_z_xyy_yyyzz_0[i] * fe_0 -
                                 ta1_z_xyy_yyyzz_1[i] * fe_0 + ta1_z_xyy_xyyyzz_0[i] * pa_x[i] - ta1_z_xyy_xyyyzz_1[i] * pc_x[i];

        ta1_z_xxyy_xyyzzz_0[i] = ta1_z_yy_xyyzzz_0[i] * fe_0 - ta1_z_yy_xyyzzz_1[i] * fe_0 + ta1_z_xyy_yyzzz_0[i] * fe_0 -
                                 ta1_z_xyy_yyzzz_1[i] * fe_0 + ta1_z_xyy_xyyzzz_0[i] * pa_x[i] - ta1_z_xyy_xyyzzz_1[i] * pc_x[i];

        ta1_z_xxyy_xyzzzz_0[i] = ta1_z_yy_xyzzzz_0[i] * fe_0 - ta1_z_yy_xyzzzz_1[i] * fe_0 + ta1_z_xyy_yzzzz_0[i] * fe_0 -
                                 ta1_z_xyy_yzzzz_1[i] * fe_0 + ta1_z_xyy_xyzzzz_0[i] * pa_x[i] - ta1_z_xyy_xyzzzz_1[i] * pc_x[i];

        ta1_z_xxyy_xzzzzz_0[i] =
            ta1_z_xx_xzzzzz_0[i] * fe_0 - ta1_z_xx_xzzzzz_1[i] * fe_0 + ta1_z_xxy_xzzzzz_0[i] * pa_y[i] - ta1_z_xxy_xzzzzz_1[i] * pc_y[i];

        ta1_z_xxyy_yyyyyy_0[i] =
            ta1_z_yy_yyyyyy_0[i] * fe_0 - ta1_z_yy_yyyyyy_1[i] * fe_0 + ta1_z_xyy_yyyyyy_0[i] * pa_x[i] - ta1_z_xyy_yyyyyy_1[i] * pc_x[i];

        ta1_z_xxyy_yyyyyz_0[i] =
            ta1_z_yy_yyyyyz_0[i] * fe_0 - ta1_z_yy_yyyyyz_1[i] * fe_0 + ta1_z_xyy_yyyyyz_0[i] * pa_x[i] - ta1_z_xyy_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxyy_yyyyzz_0[i] =
            ta1_z_yy_yyyyzz_0[i] * fe_0 - ta1_z_yy_yyyyzz_1[i] * fe_0 + ta1_z_xyy_yyyyzz_0[i] * pa_x[i] - ta1_z_xyy_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxyy_yyyzzz_0[i] =
            ta1_z_yy_yyyzzz_0[i] * fe_0 - ta1_z_yy_yyyzzz_1[i] * fe_0 + ta1_z_xyy_yyyzzz_0[i] * pa_x[i] - ta1_z_xyy_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxyy_yyzzzz_0[i] =
            ta1_z_yy_yyzzzz_0[i] * fe_0 - ta1_z_yy_yyzzzz_1[i] * fe_0 + ta1_z_xyy_yyzzzz_0[i] * pa_x[i] - ta1_z_xyy_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxyy_yzzzzz_0[i] =
            ta1_z_yy_yzzzzz_0[i] * fe_0 - ta1_z_yy_yzzzzz_1[i] * fe_0 + ta1_z_xyy_yzzzzz_0[i] * pa_x[i] - ta1_z_xyy_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxyy_zzzzzz_0[i] =
            ta1_z_yy_zzzzzz_0[i] * fe_0 - ta1_z_yy_zzzzzz_1[i] * fe_0 + ta1_z_xyy_zzzzzz_0[i] * pa_x[i] - ta1_z_xyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 952-980 components of targeted buffer : GI

    auto ta1_z_xxyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 952);

    auto ta1_z_xxyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 953);

    auto ta1_z_xxyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 954);

    auto ta1_z_xxyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 955);

    auto ta1_z_xxyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 956);

    auto ta1_z_xxyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 957);

    auto ta1_z_xxyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 958);

    auto ta1_z_xxyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 959);

    auto ta1_z_xxyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 960);

    auto ta1_z_xxyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 961);

    auto ta1_z_xxyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 962);

    auto ta1_z_xxyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 963);

    auto ta1_z_xxyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 964);

    auto ta1_z_xxyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 965);

    auto ta1_z_xxyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 966);

    auto ta1_z_xxyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 967);

    auto ta1_z_xxyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 968);

    auto ta1_z_xxyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 969);

    auto ta1_z_xxyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 970);

    auto ta1_z_xxyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 971);

    auto ta1_z_xxyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 972);

    auto ta1_z_xxyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 973);

    auto ta1_z_xxyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 974);

    auto ta1_z_xxyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 975);

    auto ta1_z_xxyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 976);

    auto ta1_z_xxyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 977);

    auto ta1_z_xxyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 978);

    auto ta1_z_xxyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 979);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pa_z,                \
                             pc_x,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_z_xxy_xxxxxy_0,  \
                             ta1_z_xxy_xxxxxy_1,  \
                             ta1_z_xxy_xxxxyy_0,  \
                             ta1_z_xxy_xxxxyy_1,  \
                             ta1_z_xxy_xxxyyy_0,  \
                             ta1_z_xxy_xxxyyy_1,  \
                             ta1_z_xxy_xxyyyy_0,  \
                             ta1_z_xxy_xxyyyy_1,  \
                             ta1_z_xxy_xyyyyy_0,  \
                             ta1_z_xxy_xyyyyy_1,  \
                             ta1_z_xxy_yyyyyy_0,  \
                             ta1_z_xxy_yyyyyy_1,  \
                             ta1_z_xxyz_xxxxxx_0, \
                             ta1_z_xxyz_xxxxxy_0, \
                             ta1_z_xxyz_xxxxxz_0, \
                             ta1_z_xxyz_xxxxyy_0, \
                             ta1_z_xxyz_xxxxyz_0, \
                             ta1_z_xxyz_xxxxzz_0, \
                             ta1_z_xxyz_xxxyyy_0, \
                             ta1_z_xxyz_xxxyyz_0, \
                             ta1_z_xxyz_xxxyzz_0, \
                             ta1_z_xxyz_xxxzzz_0, \
                             ta1_z_xxyz_xxyyyy_0, \
                             ta1_z_xxyz_xxyyyz_0, \
                             ta1_z_xxyz_xxyyzz_0, \
                             ta1_z_xxyz_xxyzzz_0, \
                             ta1_z_xxyz_xxzzzz_0, \
                             ta1_z_xxyz_xyyyyy_0, \
                             ta1_z_xxyz_xyyyyz_0, \
                             ta1_z_xxyz_xyyyzz_0, \
                             ta1_z_xxyz_xyyzzz_0, \
                             ta1_z_xxyz_xyzzzz_0, \
                             ta1_z_xxyz_xzzzzz_0, \
                             ta1_z_xxyz_yyyyyy_0, \
                             ta1_z_xxyz_yyyyyz_0, \
                             ta1_z_xxyz_yyyyzz_0, \
                             ta1_z_xxyz_yyyzzz_0, \
                             ta1_z_xxyz_yyzzzz_0, \
                             ta1_z_xxyz_yzzzzz_0, \
                             ta1_z_xxyz_zzzzzz_0, \
                             ta1_z_xxz_xxxxxx_0,  \
                             ta1_z_xxz_xxxxxx_1,  \
                             ta1_z_xxz_xxxxxz_0,  \
                             ta1_z_xxz_xxxxxz_1,  \
                             ta1_z_xxz_xxxxyz_0,  \
                             ta1_z_xxz_xxxxyz_1,  \
                             ta1_z_xxz_xxxxz_0,   \
                             ta1_z_xxz_xxxxz_1,   \
                             ta1_z_xxz_xxxxzz_0,  \
                             ta1_z_xxz_xxxxzz_1,  \
                             ta1_z_xxz_xxxyyz_0,  \
                             ta1_z_xxz_xxxyyz_1,  \
                             ta1_z_xxz_xxxyz_0,   \
                             ta1_z_xxz_xxxyz_1,   \
                             ta1_z_xxz_xxxyzz_0,  \
                             ta1_z_xxz_xxxyzz_1,  \
                             ta1_z_xxz_xxxzz_0,   \
                             ta1_z_xxz_xxxzz_1,   \
                             ta1_z_xxz_xxxzzz_0,  \
                             ta1_z_xxz_xxxzzz_1,  \
                             ta1_z_xxz_xxyyyz_0,  \
                             ta1_z_xxz_xxyyyz_1,  \
                             ta1_z_xxz_xxyyz_0,   \
                             ta1_z_xxz_xxyyz_1,   \
                             ta1_z_xxz_xxyyzz_0,  \
                             ta1_z_xxz_xxyyzz_1,  \
                             ta1_z_xxz_xxyzz_0,   \
                             ta1_z_xxz_xxyzz_1,   \
                             ta1_z_xxz_xxyzzz_0,  \
                             ta1_z_xxz_xxyzzz_1,  \
                             ta1_z_xxz_xxzzz_0,   \
                             ta1_z_xxz_xxzzz_1,   \
                             ta1_z_xxz_xxzzzz_0,  \
                             ta1_z_xxz_xxzzzz_1,  \
                             ta1_z_xxz_xyyyyz_0,  \
                             ta1_z_xxz_xyyyyz_1,  \
                             ta1_z_xxz_xyyyz_0,   \
                             ta1_z_xxz_xyyyz_1,   \
                             ta1_z_xxz_xyyyzz_0,  \
                             ta1_z_xxz_xyyyzz_1,  \
                             ta1_z_xxz_xyyzz_0,   \
                             ta1_z_xxz_xyyzz_1,   \
                             ta1_z_xxz_xyyzzz_0,  \
                             ta1_z_xxz_xyyzzz_1,  \
                             ta1_z_xxz_xyzzz_0,   \
                             ta1_z_xxz_xyzzz_1,   \
                             ta1_z_xxz_xyzzzz_0,  \
                             ta1_z_xxz_xyzzzz_1,  \
                             ta1_z_xxz_xzzzz_0,   \
                             ta1_z_xxz_xzzzz_1,   \
                             ta1_z_xxz_xzzzzz_0,  \
                             ta1_z_xxz_xzzzzz_1,  \
                             ta1_z_xxz_zzzzzz_0,  \
                             ta1_z_xxz_zzzzzz_1,  \
                             ta1_z_xyz_yyyyyz_0,  \
                             ta1_z_xyz_yyyyyz_1,  \
                             ta1_z_xyz_yyyyzz_0,  \
                             ta1_z_xyz_yyyyzz_1,  \
                             ta1_z_xyz_yyyzzz_0,  \
                             ta1_z_xyz_yyyzzz_1,  \
                             ta1_z_xyz_yyzzzz_0,  \
                             ta1_z_xyz_yyzzzz_1,  \
                             ta1_z_xyz_yzzzzz_0,  \
                             ta1_z_xyz_yzzzzz_1,  \
                             ta1_z_yz_yyyyyz_0,   \
                             ta1_z_yz_yyyyyz_1,   \
                             ta1_z_yz_yyyyzz_0,   \
                             ta1_z_yz_yyyyzz_1,   \
                             ta1_z_yz_yyyzzz_0,   \
                             ta1_z_yz_yyyzzz_1,   \
                             ta1_z_yz_yyzzzz_0,   \
                             ta1_z_yz_yyzzzz_1,   \
                             ta1_z_yz_yzzzzz_0,   \
                             ta1_z_yz_yzzzzz_1,   \
                             ta_xxy_xxxxxy_1,     \
                             ta_xxy_xxxxyy_1,     \
                             ta_xxy_xxxyyy_1,     \
                             ta_xxy_xxyyyy_1,     \
                             ta_xxy_xyyyyy_1,     \
                             ta_xxy_yyyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyz_xxxxxx_0[i] = ta1_z_xxz_xxxxxx_0[i] * pa_y[i] - ta1_z_xxz_xxxxxx_1[i] * pc_y[i];

        ta1_z_xxyz_xxxxxy_0[i] = ta_xxy_xxxxxy_1[i] + ta1_z_xxy_xxxxxy_0[i] * pa_z[i] - ta1_z_xxy_xxxxxy_1[i] * pc_z[i];

        ta1_z_xxyz_xxxxxz_0[i] = ta1_z_xxz_xxxxxz_0[i] * pa_y[i] - ta1_z_xxz_xxxxxz_1[i] * pc_y[i];

        ta1_z_xxyz_xxxxyy_0[i] = ta_xxy_xxxxyy_1[i] + ta1_z_xxy_xxxxyy_0[i] * pa_z[i] - ta1_z_xxy_xxxxyy_1[i] * pc_z[i];

        ta1_z_xxyz_xxxxyz_0[i] =
            ta1_z_xxz_xxxxz_0[i] * fe_0 - ta1_z_xxz_xxxxz_1[i] * fe_0 + ta1_z_xxz_xxxxyz_0[i] * pa_y[i] - ta1_z_xxz_xxxxyz_1[i] * pc_y[i];

        ta1_z_xxyz_xxxxzz_0[i] = ta1_z_xxz_xxxxzz_0[i] * pa_y[i] - ta1_z_xxz_xxxxzz_1[i] * pc_y[i];

        ta1_z_xxyz_xxxyyy_0[i] = ta_xxy_xxxyyy_1[i] + ta1_z_xxy_xxxyyy_0[i] * pa_z[i] - ta1_z_xxy_xxxyyy_1[i] * pc_z[i];

        ta1_z_xxyz_xxxyyz_0[i] =
            2.0 * ta1_z_xxz_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxxyz_1[i] * fe_0 + ta1_z_xxz_xxxyyz_0[i] * pa_y[i] - ta1_z_xxz_xxxyyz_1[i] * pc_y[i];

        ta1_z_xxyz_xxxyzz_0[i] =
            ta1_z_xxz_xxxzz_0[i] * fe_0 - ta1_z_xxz_xxxzz_1[i] * fe_0 + ta1_z_xxz_xxxyzz_0[i] * pa_y[i] - ta1_z_xxz_xxxyzz_1[i] * pc_y[i];

        ta1_z_xxyz_xxxzzz_0[i] = ta1_z_xxz_xxxzzz_0[i] * pa_y[i] - ta1_z_xxz_xxxzzz_1[i] * pc_y[i];

        ta1_z_xxyz_xxyyyy_0[i] = ta_xxy_xxyyyy_1[i] + ta1_z_xxy_xxyyyy_0[i] * pa_z[i] - ta1_z_xxy_xxyyyy_1[i] * pc_z[i];

        ta1_z_xxyz_xxyyyz_0[i] =
            3.0 * ta1_z_xxz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_xxz_xxyyz_1[i] * fe_0 + ta1_z_xxz_xxyyyz_0[i] * pa_y[i] - ta1_z_xxz_xxyyyz_1[i] * pc_y[i];

        ta1_z_xxyz_xxyyzz_0[i] =
            2.0 * ta1_z_xxz_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxyzz_1[i] * fe_0 + ta1_z_xxz_xxyyzz_0[i] * pa_y[i] - ta1_z_xxz_xxyyzz_1[i] * pc_y[i];

        ta1_z_xxyz_xxyzzz_0[i] =
            ta1_z_xxz_xxzzz_0[i] * fe_0 - ta1_z_xxz_xxzzz_1[i] * fe_0 + ta1_z_xxz_xxyzzz_0[i] * pa_y[i] - ta1_z_xxz_xxyzzz_1[i] * pc_y[i];

        ta1_z_xxyz_xxzzzz_0[i] = ta1_z_xxz_xxzzzz_0[i] * pa_y[i] - ta1_z_xxz_xxzzzz_1[i] * pc_y[i];

        ta1_z_xxyz_xyyyyy_0[i] = ta_xxy_xyyyyy_1[i] + ta1_z_xxy_xyyyyy_0[i] * pa_z[i] - ta1_z_xxy_xyyyyy_1[i] * pc_z[i];

        ta1_z_xxyz_xyyyyz_0[i] =
            4.0 * ta1_z_xxz_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_xxz_xyyyz_1[i] * fe_0 + ta1_z_xxz_xyyyyz_0[i] * pa_y[i] - ta1_z_xxz_xyyyyz_1[i] * pc_y[i];

        ta1_z_xxyz_xyyyzz_0[i] =
            3.0 * ta1_z_xxz_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_xxz_xyyzz_1[i] * fe_0 + ta1_z_xxz_xyyyzz_0[i] * pa_y[i] - ta1_z_xxz_xyyyzz_1[i] * pc_y[i];

        ta1_z_xxyz_xyyzzz_0[i] =
            2.0 * ta1_z_xxz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_xxz_xyzzz_1[i] * fe_0 + ta1_z_xxz_xyyzzz_0[i] * pa_y[i] - ta1_z_xxz_xyyzzz_1[i] * pc_y[i];

        ta1_z_xxyz_xyzzzz_0[i] =
            ta1_z_xxz_xzzzz_0[i] * fe_0 - ta1_z_xxz_xzzzz_1[i] * fe_0 + ta1_z_xxz_xyzzzz_0[i] * pa_y[i] - ta1_z_xxz_xyzzzz_1[i] * pc_y[i];

        ta1_z_xxyz_xzzzzz_0[i] = ta1_z_xxz_xzzzzz_0[i] * pa_y[i] - ta1_z_xxz_xzzzzz_1[i] * pc_y[i];

        ta1_z_xxyz_yyyyyy_0[i] = ta_xxy_yyyyyy_1[i] + ta1_z_xxy_yyyyyy_0[i] * pa_z[i] - ta1_z_xxy_yyyyyy_1[i] * pc_z[i];

        ta1_z_xxyz_yyyyyz_0[i] =
            ta1_z_yz_yyyyyz_0[i] * fe_0 - ta1_z_yz_yyyyyz_1[i] * fe_0 + ta1_z_xyz_yyyyyz_0[i] * pa_x[i] - ta1_z_xyz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxyz_yyyyzz_0[i] =
            ta1_z_yz_yyyyzz_0[i] * fe_0 - ta1_z_yz_yyyyzz_1[i] * fe_0 + ta1_z_xyz_yyyyzz_0[i] * pa_x[i] - ta1_z_xyz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxyz_yyyzzz_0[i] =
            ta1_z_yz_yyyzzz_0[i] * fe_0 - ta1_z_yz_yyyzzz_1[i] * fe_0 + ta1_z_xyz_yyyzzz_0[i] * pa_x[i] - ta1_z_xyz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxyz_yyzzzz_0[i] =
            ta1_z_yz_yyzzzz_0[i] * fe_0 - ta1_z_yz_yyzzzz_1[i] * fe_0 + ta1_z_xyz_yyzzzz_0[i] * pa_x[i] - ta1_z_xyz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxyz_yzzzzz_0[i] =
            ta1_z_yz_yzzzzz_0[i] * fe_0 - ta1_z_yz_yzzzzz_1[i] * fe_0 + ta1_z_xyz_yzzzzz_0[i] * pa_x[i] - ta1_z_xyz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxyz_zzzzzz_0[i] = ta1_z_xxz_zzzzzz_0[i] * pa_y[i] - ta1_z_xxz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 980-1008 components of targeted buffer : GI

    auto ta1_z_xxzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 980);

    auto ta1_z_xxzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 981);

    auto ta1_z_xxzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 982);

    auto ta1_z_xxzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 983);

    auto ta1_z_xxzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 984);

    auto ta1_z_xxzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 985);

    auto ta1_z_xxzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 986);

    auto ta1_z_xxzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 987);

    auto ta1_z_xxzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 988);

    auto ta1_z_xxzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 989);

    auto ta1_z_xxzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 990);

    auto ta1_z_xxzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 991);

    auto ta1_z_xxzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 992);

    auto ta1_z_xxzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 993);

    auto ta1_z_xxzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 994);

    auto ta1_z_xxzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 995);

    auto ta1_z_xxzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 996);

    auto ta1_z_xxzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 997);

    auto ta1_z_xxzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 998);

    auto ta1_z_xxzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 999);

    auto ta1_z_xxzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1000);

    auto ta1_z_xxzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1001);

    auto ta1_z_xxzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1002);

    auto ta1_z_xxzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1003);

    auto ta1_z_xxzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1004);

    auto ta1_z_xxzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1005);

    auto ta1_z_xxzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1006);

    auto ta1_z_xxzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1007);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_z_xx_xxxxxx_0,   \
                             ta1_z_xx_xxxxxx_1,   \
                             ta1_z_xx_xxxxxy_0,   \
                             ta1_z_xx_xxxxxy_1,   \
                             ta1_z_xx_xxxxyy_0,   \
                             ta1_z_xx_xxxxyy_1,   \
                             ta1_z_xx_xxxyyy_0,   \
                             ta1_z_xx_xxxyyy_1,   \
                             ta1_z_xx_xxyyyy_0,   \
                             ta1_z_xx_xxyyyy_1,   \
                             ta1_z_xx_xyyyyy_0,   \
                             ta1_z_xx_xyyyyy_1,   \
                             ta1_z_xxz_xxxxxx_0,  \
                             ta1_z_xxz_xxxxxx_1,  \
                             ta1_z_xxz_xxxxxy_0,  \
                             ta1_z_xxz_xxxxxy_1,  \
                             ta1_z_xxz_xxxxyy_0,  \
                             ta1_z_xxz_xxxxyy_1,  \
                             ta1_z_xxz_xxxyyy_0,  \
                             ta1_z_xxz_xxxyyy_1,  \
                             ta1_z_xxz_xxyyyy_0,  \
                             ta1_z_xxz_xxyyyy_1,  \
                             ta1_z_xxz_xyyyyy_0,  \
                             ta1_z_xxz_xyyyyy_1,  \
                             ta1_z_xxzz_xxxxxx_0, \
                             ta1_z_xxzz_xxxxxy_0, \
                             ta1_z_xxzz_xxxxxz_0, \
                             ta1_z_xxzz_xxxxyy_0, \
                             ta1_z_xxzz_xxxxyz_0, \
                             ta1_z_xxzz_xxxxzz_0, \
                             ta1_z_xxzz_xxxyyy_0, \
                             ta1_z_xxzz_xxxyyz_0, \
                             ta1_z_xxzz_xxxyzz_0, \
                             ta1_z_xxzz_xxxzzz_0, \
                             ta1_z_xxzz_xxyyyy_0, \
                             ta1_z_xxzz_xxyyyz_0, \
                             ta1_z_xxzz_xxyyzz_0, \
                             ta1_z_xxzz_xxyzzz_0, \
                             ta1_z_xxzz_xxzzzz_0, \
                             ta1_z_xxzz_xyyyyy_0, \
                             ta1_z_xxzz_xyyyyz_0, \
                             ta1_z_xxzz_xyyyzz_0, \
                             ta1_z_xxzz_xyyzzz_0, \
                             ta1_z_xxzz_xyzzzz_0, \
                             ta1_z_xxzz_xzzzzz_0, \
                             ta1_z_xxzz_yyyyyy_0, \
                             ta1_z_xxzz_yyyyyz_0, \
                             ta1_z_xxzz_yyyyzz_0, \
                             ta1_z_xxzz_yyyzzz_0, \
                             ta1_z_xxzz_yyzzzz_0, \
                             ta1_z_xxzz_yzzzzz_0, \
                             ta1_z_xxzz_zzzzzz_0, \
                             ta1_z_xzz_xxxxxz_0,  \
                             ta1_z_xzz_xxxxxz_1,  \
                             ta1_z_xzz_xxxxyz_0,  \
                             ta1_z_xzz_xxxxyz_1,  \
                             ta1_z_xzz_xxxxz_0,   \
                             ta1_z_xzz_xxxxz_1,   \
                             ta1_z_xzz_xxxxzz_0,  \
                             ta1_z_xzz_xxxxzz_1,  \
                             ta1_z_xzz_xxxyyz_0,  \
                             ta1_z_xzz_xxxyyz_1,  \
                             ta1_z_xzz_xxxyz_0,   \
                             ta1_z_xzz_xxxyz_1,   \
                             ta1_z_xzz_xxxyzz_0,  \
                             ta1_z_xzz_xxxyzz_1,  \
                             ta1_z_xzz_xxxzz_0,   \
                             ta1_z_xzz_xxxzz_1,   \
                             ta1_z_xzz_xxxzzz_0,  \
                             ta1_z_xzz_xxxzzz_1,  \
                             ta1_z_xzz_xxyyyz_0,  \
                             ta1_z_xzz_xxyyyz_1,  \
                             ta1_z_xzz_xxyyz_0,   \
                             ta1_z_xzz_xxyyz_1,   \
                             ta1_z_xzz_xxyyzz_0,  \
                             ta1_z_xzz_xxyyzz_1,  \
                             ta1_z_xzz_xxyzz_0,   \
                             ta1_z_xzz_xxyzz_1,   \
                             ta1_z_xzz_xxyzzz_0,  \
                             ta1_z_xzz_xxyzzz_1,  \
                             ta1_z_xzz_xxzzz_0,   \
                             ta1_z_xzz_xxzzz_1,   \
                             ta1_z_xzz_xxzzzz_0,  \
                             ta1_z_xzz_xxzzzz_1,  \
                             ta1_z_xzz_xyyyyz_0,  \
                             ta1_z_xzz_xyyyyz_1,  \
                             ta1_z_xzz_xyyyz_0,   \
                             ta1_z_xzz_xyyyz_1,   \
                             ta1_z_xzz_xyyyzz_0,  \
                             ta1_z_xzz_xyyyzz_1,  \
                             ta1_z_xzz_xyyzz_0,   \
                             ta1_z_xzz_xyyzz_1,   \
                             ta1_z_xzz_xyyzzz_0,  \
                             ta1_z_xzz_xyyzzz_1,  \
                             ta1_z_xzz_xyzzz_0,   \
                             ta1_z_xzz_xyzzz_1,   \
                             ta1_z_xzz_xyzzzz_0,  \
                             ta1_z_xzz_xyzzzz_1,  \
                             ta1_z_xzz_xzzzz_0,   \
                             ta1_z_xzz_xzzzz_1,   \
                             ta1_z_xzz_xzzzzz_0,  \
                             ta1_z_xzz_xzzzzz_1,  \
                             ta1_z_xzz_yyyyyy_0,  \
                             ta1_z_xzz_yyyyyy_1,  \
                             ta1_z_xzz_yyyyyz_0,  \
                             ta1_z_xzz_yyyyyz_1,  \
                             ta1_z_xzz_yyyyz_0,   \
                             ta1_z_xzz_yyyyz_1,   \
                             ta1_z_xzz_yyyyzz_0,  \
                             ta1_z_xzz_yyyyzz_1,  \
                             ta1_z_xzz_yyyzz_0,   \
                             ta1_z_xzz_yyyzz_1,   \
                             ta1_z_xzz_yyyzzz_0,  \
                             ta1_z_xzz_yyyzzz_1,  \
                             ta1_z_xzz_yyzzz_0,   \
                             ta1_z_xzz_yyzzz_1,   \
                             ta1_z_xzz_yyzzzz_0,  \
                             ta1_z_xzz_yyzzzz_1,  \
                             ta1_z_xzz_yzzzz_0,   \
                             ta1_z_xzz_yzzzz_1,   \
                             ta1_z_xzz_yzzzzz_0,  \
                             ta1_z_xzz_yzzzzz_1,  \
                             ta1_z_xzz_zzzzz_0,   \
                             ta1_z_xzz_zzzzz_1,   \
                             ta1_z_xzz_zzzzzz_0,  \
                             ta1_z_xzz_zzzzzz_1,  \
                             ta1_z_zz_xxxxxz_0,   \
                             ta1_z_zz_xxxxxz_1,   \
                             ta1_z_zz_xxxxyz_0,   \
                             ta1_z_zz_xxxxyz_1,   \
                             ta1_z_zz_xxxxzz_0,   \
                             ta1_z_zz_xxxxzz_1,   \
                             ta1_z_zz_xxxyyz_0,   \
                             ta1_z_zz_xxxyyz_1,   \
                             ta1_z_zz_xxxyzz_0,   \
                             ta1_z_zz_xxxyzz_1,   \
                             ta1_z_zz_xxxzzz_0,   \
                             ta1_z_zz_xxxzzz_1,   \
                             ta1_z_zz_xxyyyz_0,   \
                             ta1_z_zz_xxyyyz_1,   \
                             ta1_z_zz_xxyyzz_0,   \
                             ta1_z_zz_xxyyzz_1,   \
                             ta1_z_zz_xxyzzz_0,   \
                             ta1_z_zz_xxyzzz_1,   \
                             ta1_z_zz_xxzzzz_0,   \
                             ta1_z_zz_xxzzzz_1,   \
                             ta1_z_zz_xyyyyz_0,   \
                             ta1_z_zz_xyyyyz_1,   \
                             ta1_z_zz_xyyyzz_0,   \
                             ta1_z_zz_xyyyzz_1,   \
                             ta1_z_zz_xyyzzz_0,   \
                             ta1_z_zz_xyyzzz_1,   \
                             ta1_z_zz_xyzzzz_0,   \
                             ta1_z_zz_xyzzzz_1,   \
                             ta1_z_zz_xzzzzz_0,   \
                             ta1_z_zz_xzzzzz_1,   \
                             ta1_z_zz_yyyyyy_0,   \
                             ta1_z_zz_yyyyyy_1,   \
                             ta1_z_zz_yyyyyz_0,   \
                             ta1_z_zz_yyyyyz_1,   \
                             ta1_z_zz_yyyyzz_0,   \
                             ta1_z_zz_yyyyzz_1,   \
                             ta1_z_zz_yyyzzz_0,   \
                             ta1_z_zz_yyyzzz_1,   \
                             ta1_z_zz_yyzzzz_0,   \
                             ta1_z_zz_yyzzzz_1,   \
                             ta1_z_zz_yzzzzz_0,   \
                             ta1_z_zz_yzzzzz_1,   \
                             ta1_z_zz_zzzzzz_0,   \
                             ta1_z_zz_zzzzzz_1,   \
                             ta_xxz_xxxxxx_1,     \
                             ta_xxz_xxxxxy_1,     \
                             ta_xxz_xxxxyy_1,     \
                             ta_xxz_xxxyyy_1,     \
                             ta_xxz_xxyyyy_1,     \
                             ta_xxz_xyyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzz_xxxxxx_0[i] = ta1_z_xx_xxxxxx_0[i] * fe_0 - ta1_z_xx_xxxxxx_1[i] * fe_0 + ta_xxz_xxxxxx_1[i] + ta1_z_xxz_xxxxxx_0[i] * pa_z[i] -
                                 ta1_z_xxz_xxxxxx_1[i] * pc_z[i];

        ta1_z_xxzz_xxxxxy_0[i] = ta1_z_xx_xxxxxy_0[i] * fe_0 - ta1_z_xx_xxxxxy_1[i] * fe_0 + ta_xxz_xxxxxy_1[i] + ta1_z_xxz_xxxxxy_0[i] * pa_z[i] -
                                 ta1_z_xxz_xxxxxy_1[i] * pc_z[i];

        ta1_z_xxzz_xxxxxz_0[i] = ta1_z_zz_xxxxxz_0[i] * fe_0 - ta1_z_zz_xxxxxz_1[i] * fe_0 + 5.0 * ta1_z_xzz_xxxxz_0[i] * fe_0 -
                                 5.0 * ta1_z_xzz_xxxxz_1[i] * fe_0 + ta1_z_xzz_xxxxxz_0[i] * pa_x[i] - ta1_z_xzz_xxxxxz_1[i] * pc_x[i];

        ta1_z_xxzz_xxxxyy_0[i] = ta1_z_xx_xxxxyy_0[i] * fe_0 - ta1_z_xx_xxxxyy_1[i] * fe_0 + ta_xxz_xxxxyy_1[i] + ta1_z_xxz_xxxxyy_0[i] * pa_z[i] -
                                 ta1_z_xxz_xxxxyy_1[i] * pc_z[i];

        ta1_z_xxzz_xxxxyz_0[i] = ta1_z_zz_xxxxyz_0[i] * fe_0 - ta1_z_zz_xxxxyz_1[i] * fe_0 + 4.0 * ta1_z_xzz_xxxyz_0[i] * fe_0 -
                                 4.0 * ta1_z_xzz_xxxyz_1[i] * fe_0 + ta1_z_xzz_xxxxyz_0[i] * pa_x[i] - ta1_z_xzz_xxxxyz_1[i] * pc_x[i];

        ta1_z_xxzz_xxxxzz_0[i] = ta1_z_zz_xxxxzz_0[i] * fe_0 - ta1_z_zz_xxxxzz_1[i] * fe_0 + 4.0 * ta1_z_xzz_xxxzz_0[i] * fe_0 -
                                 4.0 * ta1_z_xzz_xxxzz_1[i] * fe_0 + ta1_z_xzz_xxxxzz_0[i] * pa_x[i] - ta1_z_xzz_xxxxzz_1[i] * pc_x[i];

        ta1_z_xxzz_xxxyyy_0[i] = ta1_z_xx_xxxyyy_0[i] * fe_0 - ta1_z_xx_xxxyyy_1[i] * fe_0 + ta_xxz_xxxyyy_1[i] + ta1_z_xxz_xxxyyy_0[i] * pa_z[i] -
                                 ta1_z_xxz_xxxyyy_1[i] * pc_z[i];

        ta1_z_xxzz_xxxyyz_0[i] = ta1_z_zz_xxxyyz_0[i] * fe_0 - ta1_z_zz_xxxyyz_1[i] * fe_0 + 3.0 * ta1_z_xzz_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_z_xzz_xxyyz_1[i] * fe_0 + ta1_z_xzz_xxxyyz_0[i] * pa_x[i] - ta1_z_xzz_xxxyyz_1[i] * pc_x[i];

        ta1_z_xxzz_xxxyzz_0[i] = ta1_z_zz_xxxyzz_0[i] * fe_0 - ta1_z_zz_xxxyzz_1[i] * fe_0 + 3.0 * ta1_z_xzz_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_xzz_xxyzz_1[i] * fe_0 + ta1_z_xzz_xxxyzz_0[i] * pa_x[i] - ta1_z_xzz_xxxyzz_1[i] * pc_x[i];

        ta1_z_xxzz_xxxzzz_0[i] = ta1_z_zz_xxxzzz_0[i] * fe_0 - ta1_z_zz_xxxzzz_1[i] * fe_0 + 3.0 * ta1_z_xzz_xxzzz_0[i] * fe_0 -
                                 3.0 * ta1_z_xzz_xxzzz_1[i] * fe_0 + ta1_z_xzz_xxxzzz_0[i] * pa_x[i] - ta1_z_xzz_xxxzzz_1[i] * pc_x[i];

        ta1_z_xxzz_xxyyyy_0[i] = ta1_z_xx_xxyyyy_0[i] * fe_0 - ta1_z_xx_xxyyyy_1[i] * fe_0 + ta_xxz_xxyyyy_1[i] + ta1_z_xxz_xxyyyy_0[i] * pa_z[i] -
                                 ta1_z_xxz_xxyyyy_1[i] * pc_z[i];

        ta1_z_xxzz_xxyyyz_0[i] = ta1_z_zz_xxyyyz_0[i] * fe_0 - ta1_z_zz_xxyyyz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_z_xzz_xyyyz_1[i] * fe_0 + ta1_z_xzz_xxyyyz_0[i] * pa_x[i] - ta1_z_xzz_xxyyyz_1[i] * pc_x[i];

        ta1_z_xxzz_xxyyzz_0[i] = ta1_z_zz_xxyyzz_0[i] * fe_0 - ta1_z_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xyyzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xzz_xyyzz_1[i] * fe_0 + ta1_z_xzz_xxyyzz_0[i] * pa_x[i] - ta1_z_xzz_xxyyzz_1[i] * pc_x[i];

        ta1_z_xxzz_xxyzzz_0[i] = ta1_z_zz_xxyzzz_0[i] * fe_0 - ta1_z_zz_xxyzzz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xzz_xyzzz_1[i] * fe_0 + ta1_z_xzz_xxyzzz_0[i] * pa_x[i] - ta1_z_xzz_xxyzzz_1[i] * pc_x[i];

        ta1_z_xxzz_xxzzzz_0[i] = ta1_z_zz_xxzzzz_0[i] * fe_0 - ta1_z_zz_xxzzzz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xzzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_xzz_xzzzz_1[i] * fe_0 + ta1_z_xzz_xxzzzz_0[i] * pa_x[i] - ta1_z_xzz_xxzzzz_1[i] * pc_x[i];

        ta1_z_xxzz_xyyyyy_0[i] = ta1_z_xx_xyyyyy_0[i] * fe_0 - ta1_z_xx_xyyyyy_1[i] * fe_0 + ta_xxz_xyyyyy_1[i] + ta1_z_xxz_xyyyyy_0[i] * pa_z[i] -
                                 ta1_z_xxz_xyyyyy_1[i] * pc_z[i];

        ta1_z_xxzz_xyyyyz_0[i] = ta1_z_zz_xyyyyz_0[i] * fe_0 - ta1_z_zz_xyyyyz_1[i] * fe_0 + ta1_z_xzz_yyyyz_0[i] * fe_0 -
                                 ta1_z_xzz_yyyyz_1[i] * fe_0 + ta1_z_xzz_xyyyyz_0[i] * pa_x[i] - ta1_z_xzz_xyyyyz_1[i] * pc_x[i];

        ta1_z_xxzz_xyyyzz_0[i] = ta1_z_zz_xyyyzz_0[i] * fe_0 - ta1_z_zz_xyyyzz_1[i] * fe_0 + ta1_z_xzz_yyyzz_0[i] * fe_0 -
                                 ta1_z_xzz_yyyzz_1[i] * fe_0 + ta1_z_xzz_xyyyzz_0[i] * pa_x[i] - ta1_z_xzz_xyyyzz_1[i] * pc_x[i];

        ta1_z_xxzz_xyyzzz_0[i] = ta1_z_zz_xyyzzz_0[i] * fe_0 - ta1_z_zz_xyyzzz_1[i] * fe_0 + ta1_z_xzz_yyzzz_0[i] * fe_0 -
                                 ta1_z_xzz_yyzzz_1[i] * fe_0 + ta1_z_xzz_xyyzzz_0[i] * pa_x[i] - ta1_z_xzz_xyyzzz_1[i] * pc_x[i];

        ta1_z_xxzz_xyzzzz_0[i] = ta1_z_zz_xyzzzz_0[i] * fe_0 - ta1_z_zz_xyzzzz_1[i] * fe_0 + ta1_z_xzz_yzzzz_0[i] * fe_0 -
                                 ta1_z_xzz_yzzzz_1[i] * fe_0 + ta1_z_xzz_xyzzzz_0[i] * pa_x[i] - ta1_z_xzz_xyzzzz_1[i] * pc_x[i];

        ta1_z_xxzz_xzzzzz_0[i] = ta1_z_zz_xzzzzz_0[i] * fe_0 - ta1_z_zz_xzzzzz_1[i] * fe_0 + ta1_z_xzz_zzzzz_0[i] * fe_0 -
                                 ta1_z_xzz_zzzzz_1[i] * fe_0 + ta1_z_xzz_xzzzzz_0[i] * pa_x[i] - ta1_z_xzz_xzzzzz_1[i] * pc_x[i];

        ta1_z_xxzz_yyyyyy_0[i] =
            ta1_z_zz_yyyyyy_0[i] * fe_0 - ta1_z_zz_yyyyyy_1[i] * fe_0 + ta1_z_xzz_yyyyyy_0[i] * pa_x[i] - ta1_z_xzz_yyyyyy_1[i] * pc_x[i];

        ta1_z_xxzz_yyyyyz_0[i] =
            ta1_z_zz_yyyyyz_0[i] * fe_0 - ta1_z_zz_yyyyyz_1[i] * fe_0 + ta1_z_xzz_yyyyyz_0[i] * pa_x[i] - ta1_z_xzz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxzz_yyyyzz_0[i] =
            ta1_z_zz_yyyyzz_0[i] * fe_0 - ta1_z_zz_yyyyzz_1[i] * fe_0 + ta1_z_xzz_yyyyzz_0[i] * pa_x[i] - ta1_z_xzz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxzz_yyyzzz_0[i] =
            ta1_z_zz_yyyzzz_0[i] * fe_0 - ta1_z_zz_yyyzzz_1[i] * fe_0 + ta1_z_xzz_yyyzzz_0[i] * pa_x[i] - ta1_z_xzz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxzz_yyzzzz_0[i] =
            ta1_z_zz_yyzzzz_0[i] * fe_0 - ta1_z_zz_yyzzzz_1[i] * fe_0 + ta1_z_xzz_yyzzzz_0[i] * pa_x[i] - ta1_z_xzz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxzz_yzzzzz_0[i] =
            ta1_z_zz_yzzzzz_0[i] * fe_0 - ta1_z_zz_yzzzzz_1[i] * fe_0 + ta1_z_xzz_yzzzzz_0[i] * pa_x[i] - ta1_z_xzz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxzz_zzzzzz_0[i] =
            ta1_z_zz_zzzzzz_0[i] * fe_0 - ta1_z_zz_zzzzzz_1[i] * fe_0 + ta1_z_xzz_zzzzzz_0[i] * pa_x[i] - ta1_z_xzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 1008-1036 components of targeted buffer : GI

    auto ta1_z_xyyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1008);

    auto ta1_z_xyyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1009);

    auto ta1_z_xyyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1010);

    auto ta1_z_xyyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1011);

    auto ta1_z_xyyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1012);

    auto ta1_z_xyyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1013);

    auto ta1_z_xyyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1014);

    auto ta1_z_xyyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1015);

    auto ta1_z_xyyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1016);

    auto ta1_z_xyyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1017);

    auto ta1_z_xyyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1018);

    auto ta1_z_xyyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1019);

    auto ta1_z_xyyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1020);

    auto ta1_z_xyyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1021);

    auto ta1_z_xyyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1022);

    auto ta1_z_xyyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1023);

    auto ta1_z_xyyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1024);

    auto ta1_z_xyyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1025);

    auto ta1_z_xyyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1026);

    auto ta1_z_xyyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1027);

    auto ta1_z_xyyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1028);

    auto ta1_z_xyyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1029);

    auto ta1_z_xyyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1030);

    auto ta1_z_xyyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1031);

    auto ta1_z_xyyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1032);

    auto ta1_z_xyyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1033);

    auto ta1_z_xyyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1034);

    auto ta1_z_xyyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1035);

#pragma omp simd aligned(pa_x,                    \
                             pc_x,                \
                             ta1_z_xyyy_xxxxxx_0, \
                             ta1_z_xyyy_xxxxxy_0, \
                             ta1_z_xyyy_xxxxxz_0, \
                             ta1_z_xyyy_xxxxyy_0, \
                             ta1_z_xyyy_xxxxyz_0, \
                             ta1_z_xyyy_xxxxzz_0, \
                             ta1_z_xyyy_xxxyyy_0, \
                             ta1_z_xyyy_xxxyyz_0, \
                             ta1_z_xyyy_xxxyzz_0, \
                             ta1_z_xyyy_xxxzzz_0, \
                             ta1_z_xyyy_xxyyyy_0, \
                             ta1_z_xyyy_xxyyyz_0, \
                             ta1_z_xyyy_xxyyzz_0, \
                             ta1_z_xyyy_xxyzzz_0, \
                             ta1_z_xyyy_xxzzzz_0, \
                             ta1_z_xyyy_xyyyyy_0, \
                             ta1_z_xyyy_xyyyyz_0, \
                             ta1_z_xyyy_xyyyzz_0, \
                             ta1_z_xyyy_xyyzzz_0, \
                             ta1_z_xyyy_xyzzzz_0, \
                             ta1_z_xyyy_xzzzzz_0, \
                             ta1_z_xyyy_yyyyyy_0, \
                             ta1_z_xyyy_yyyyyz_0, \
                             ta1_z_xyyy_yyyyzz_0, \
                             ta1_z_xyyy_yyyzzz_0, \
                             ta1_z_xyyy_yyzzzz_0, \
                             ta1_z_xyyy_yzzzzz_0, \
                             ta1_z_xyyy_zzzzzz_0, \
                             ta1_z_yyy_xxxxx_0,   \
                             ta1_z_yyy_xxxxx_1,   \
                             ta1_z_yyy_xxxxxx_0,  \
                             ta1_z_yyy_xxxxxx_1,  \
                             ta1_z_yyy_xxxxxy_0,  \
                             ta1_z_yyy_xxxxxy_1,  \
                             ta1_z_yyy_xxxxxz_0,  \
                             ta1_z_yyy_xxxxxz_1,  \
                             ta1_z_yyy_xxxxy_0,   \
                             ta1_z_yyy_xxxxy_1,   \
                             ta1_z_yyy_xxxxyy_0,  \
                             ta1_z_yyy_xxxxyy_1,  \
                             ta1_z_yyy_xxxxyz_0,  \
                             ta1_z_yyy_xxxxyz_1,  \
                             ta1_z_yyy_xxxxz_0,   \
                             ta1_z_yyy_xxxxz_1,   \
                             ta1_z_yyy_xxxxzz_0,  \
                             ta1_z_yyy_xxxxzz_1,  \
                             ta1_z_yyy_xxxyy_0,   \
                             ta1_z_yyy_xxxyy_1,   \
                             ta1_z_yyy_xxxyyy_0,  \
                             ta1_z_yyy_xxxyyy_1,  \
                             ta1_z_yyy_xxxyyz_0,  \
                             ta1_z_yyy_xxxyyz_1,  \
                             ta1_z_yyy_xxxyz_0,   \
                             ta1_z_yyy_xxxyz_1,   \
                             ta1_z_yyy_xxxyzz_0,  \
                             ta1_z_yyy_xxxyzz_1,  \
                             ta1_z_yyy_xxxzz_0,   \
                             ta1_z_yyy_xxxzz_1,   \
                             ta1_z_yyy_xxxzzz_0,  \
                             ta1_z_yyy_xxxzzz_1,  \
                             ta1_z_yyy_xxyyy_0,   \
                             ta1_z_yyy_xxyyy_1,   \
                             ta1_z_yyy_xxyyyy_0,  \
                             ta1_z_yyy_xxyyyy_1,  \
                             ta1_z_yyy_xxyyyz_0,  \
                             ta1_z_yyy_xxyyyz_1,  \
                             ta1_z_yyy_xxyyz_0,   \
                             ta1_z_yyy_xxyyz_1,   \
                             ta1_z_yyy_xxyyzz_0,  \
                             ta1_z_yyy_xxyyzz_1,  \
                             ta1_z_yyy_xxyzz_0,   \
                             ta1_z_yyy_xxyzz_1,   \
                             ta1_z_yyy_xxyzzz_0,  \
                             ta1_z_yyy_xxyzzz_1,  \
                             ta1_z_yyy_xxzzz_0,   \
                             ta1_z_yyy_xxzzz_1,   \
                             ta1_z_yyy_xxzzzz_0,  \
                             ta1_z_yyy_xxzzzz_1,  \
                             ta1_z_yyy_xyyyy_0,   \
                             ta1_z_yyy_xyyyy_1,   \
                             ta1_z_yyy_xyyyyy_0,  \
                             ta1_z_yyy_xyyyyy_1,  \
                             ta1_z_yyy_xyyyyz_0,  \
                             ta1_z_yyy_xyyyyz_1,  \
                             ta1_z_yyy_xyyyz_0,   \
                             ta1_z_yyy_xyyyz_1,   \
                             ta1_z_yyy_xyyyzz_0,  \
                             ta1_z_yyy_xyyyzz_1,  \
                             ta1_z_yyy_xyyzz_0,   \
                             ta1_z_yyy_xyyzz_1,   \
                             ta1_z_yyy_xyyzzz_0,  \
                             ta1_z_yyy_xyyzzz_1,  \
                             ta1_z_yyy_xyzzz_0,   \
                             ta1_z_yyy_xyzzz_1,   \
                             ta1_z_yyy_xyzzzz_0,  \
                             ta1_z_yyy_xyzzzz_1,  \
                             ta1_z_yyy_xzzzz_0,   \
                             ta1_z_yyy_xzzzz_1,   \
                             ta1_z_yyy_xzzzzz_0,  \
                             ta1_z_yyy_xzzzzz_1,  \
                             ta1_z_yyy_yyyyy_0,   \
                             ta1_z_yyy_yyyyy_1,   \
                             ta1_z_yyy_yyyyyy_0,  \
                             ta1_z_yyy_yyyyyy_1,  \
                             ta1_z_yyy_yyyyyz_0,  \
                             ta1_z_yyy_yyyyyz_1,  \
                             ta1_z_yyy_yyyyz_0,   \
                             ta1_z_yyy_yyyyz_1,   \
                             ta1_z_yyy_yyyyzz_0,  \
                             ta1_z_yyy_yyyyzz_1,  \
                             ta1_z_yyy_yyyzz_0,   \
                             ta1_z_yyy_yyyzz_1,   \
                             ta1_z_yyy_yyyzzz_0,  \
                             ta1_z_yyy_yyyzzz_1,  \
                             ta1_z_yyy_yyzzz_0,   \
                             ta1_z_yyy_yyzzz_1,   \
                             ta1_z_yyy_yyzzzz_0,  \
                             ta1_z_yyy_yyzzzz_1,  \
                             ta1_z_yyy_yzzzz_0,   \
                             ta1_z_yyy_yzzzz_1,   \
                             ta1_z_yyy_yzzzzz_0,  \
                             ta1_z_yyy_yzzzzz_1,  \
                             ta1_z_yyy_zzzzz_0,   \
                             ta1_z_yyy_zzzzz_1,   \
                             ta1_z_yyy_zzzzzz_0,  \
                             ta1_z_yyy_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyy_xxxxxx_0[i] =
            6.0 * ta1_z_yyy_xxxxx_0[i] * fe_0 - 6.0 * ta1_z_yyy_xxxxx_1[i] * fe_0 + ta1_z_yyy_xxxxxx_0[i] * pa_x[i] - ta1_z_yyy_xxxxxx_1[i] * pc_x[i];

        ta1_z_xyyy_xxxxxy_0[i] =
            5.0 * ta1_z_yyy_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_yyy_xxxxy_1[i] * fe_0 + ta1_z_yyy_xxxxxy_0[i] * pa_x[i] - ta1_z_yyy_xxxxxy_1[i] * pc_x[i];

        ta1_z_xyyy_xxxxxz_0[i] =
            5.0 * ta1_z_yyy_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_yyy_xxxxz_1[i] * fe_0 + ta1_z_yyy_xxxxxz_0[i] * pa_x[i] - ta1_z_yyy_xxxxxz_1[i] * pc_x[i];

        ta1_z_xyyy_xxxxyy_0[i] =
            4.0 * ta1_z_yyy_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxyy_1[i] * fe_0 + ta1_z_yyy_xxxxyy_0[i] * pa_x[i] - ta1_z_yyy_xxxxyy_1[i] * pc_x[i];

        ta1_z_xyyy_xxxxyz_0[i] =
            4.0 * ta1_z_yyy_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxyz_1[i] * fe_0 + ta1_z_yyy_xxxxyz_0[i] * pa_x[i] - ta1_z_yyy_xxxxyz_1[i] * pc_x[i];

        ta1_z_xyyy_xxxxzz_0[i] =
            4.0 * ta1_z_yyy_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxzz_1[i] * fe_0 + ta1_z_yyy_xxxxzz_0[i] * pa_x[i] - ta1_z_yyy_xxxxzz_1[i] * pc_x[i];

        ta1_z_xyyy_xxxyyy_0[i] =
            3.0 * ta1_z_yyy_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxyyy_1[i] * fe_0 + ta1_z_yyy_xxxyyy_0[i] * pa_x[i] - ta1_z_yyy_xxxyyy_1[i] * pc_x[i];

        ta1_z_xyyy_xxxyyz_0[i] =
            3.0 * ta1_z_yyy_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxyyz_1[i] * fe_0 + ta1_z_yyy_xxxyyz_0[i] * pa_x[i] - ta1_z_yyy_xxxyyz_1[i] * pc_x[i];

        ta1_z_xyyy_xxxyzz_0[i] =
            3.0 * ta1_z_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxyzz_1[i] * fe_0 + ta1_z_yyy_xxxyzz_0[i] * pa_x[i] - ta1_z_yyy_xxxyzz_1[i] * pc_x[i];

        ta1_z_xyyy_xxxzzz_0[i] =
            3.0 * ta1_z_yyy_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxzzz_1[i] * fe_0 + ta1_z_yyy_xxxzzz_0[i] * pa_x[i] - ta1_z_yyy_xxxzzz_1[i] * pc_x[i];

        ta1_z_xyyy_xxyyyy_0[i] =
            2.0 * ta1_z_yyy_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyyyy_1[i] * fe_0 + ta1_z_yyy_xxyyyy_0[i] * pa_x[i] - ta1_z_yyy_xxyyyy_1[i] * pc_x[i];

        ta1_z_xyyy_xxyyyz_0[i] =
            2.0 * ta1_z_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyyyz_1[i] * fe_0 + ta1_z_yyy_xxyyyz_0[i] * pa_x[i] - ta1_z_yyy_xxyyyz_1[i] * pc_x[i];

        ta1_z_xyyy_xxyyzz_0[i] =
            2.0 * ta1_z_yyy_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyyzz_1[i] * fe_0 + ta1_z_yyy_xxyyzz_0[i] * pa_x[i] - ta1_z_yyy_xxyyzz_1[i] * pc_x[i];

        ta1_z_xyyy_xxyzzz_0[i] =
            2.0 * ta1_z_yyy_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyzzz_1[i] * fe_0 + ta1_z_yyy_xxyzzz_0[i] * pa_x[i] - ta1_z_yyy_xxyzzz_1[i] * pc_x[i];

        ta1_z_xyyy_xxzzzz_0[i] =
            2.0 * ta1_z_yyy_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xzzzz_1[i] * fe_0 + ta1_z_yyy_xxzzzz_0[i] * pa_x[i] - ta1_z_yyy_xxzzzz_1[i] * pc_x[i];

        ta1_z_xyyy_xyyyyy_0[i] =
            ta1_z_yyy_yyyyy_0[i] * fe_0 - ta1_z_yyy_yyyyy_1[i] * fe_0 + ta1_z_yyy_xyyyyy_0[i] * pa_x[i] - ta1_z_yyy_xyyyyy_1[i] * pc_x[i];

        ta1_z_xyyy_xyyyyz_0[i] =
            ta1_z_yyy_yyyyz_0[i] * fe_0 - ta1_z_yyy_yyyyz_1[i] * fe_0 + ta1_z_yyy_xyyyyz_0[i] * pa_x[i] - ta1_z_yyy_xyyyyz_1[i] * pc_x[i];

        ta1_z_xyyy_xyyyzz_0[i] =
            ta1_z_yyy_yyyzz_0[i] * fe_0 - ta1_z_yyy_yyyzz_1[i] * fe_0 + ta1_z_yyy_xyyyzz_0[i] * pa_x[i] - ta1_z_yyy_xyyyzz_1[i] * pc_x[i];

        ta1_z_xyyy_xyyzzz_0[i] =
            ta1_z_yyy_yyzzz_0[i] * fe_0 - ta1_z_yyy_yyzzz_1[i] * fe_0 + ta1_z_yyy_xyyzzz_0[i] * pa_x[i] - ta1_z_yyy_xyyzzz_1[i] * pc_x[i];

        ta1_z_xyyy_xyzzzz_0[i] =
            ta1_z_yyy_yzzzz_0[i] * fe_0 - ta1_z_yyy_yzzzz_1[i] * fe_0 + ta1_z_yyy_xyzzzz_0[i] * pa_x[i] - ta1_z_yyy_xyzzzz_1[i] * pc_x[i];

        ta1_z_xyyy_xzzzzz_0[i] =
            ta1_z_yyy_zzzzz_0[i] * fe_0 - ta1_z_yyy_zzzzz_1[i] * fe_0 + ta1_z_yyy_xzzzzz_0[i] * pa_x[i] - ta1_z_yyy_xzzzzz_1[i] * pc_x[i];

        ta1_z_xyyy_yyyyyy_0[i] = ta1_z_yyy_yyyyyy_0[i] * pa_x[i] - ta1_z_yyy_yyyyyy_1[i] * pc_x[i];

        ta1_z_xyyy_yyyyyz_0[i] = ta1_z_yyy_yyyyyz_0[i] * pa_x[i] - ta1_z_yyy_yyyyyz_1[i] * pc_x[i];

        ta1_z_xyyy_yyyyzz_0[i] = ta1_z_yyy_yyyyzz_0[i] * pa_x[i] - ta1_z_yyy_yyyyzz_1[i] * pc_x[i];

        ta1_z_xyyy_yyyzzz_0[i] = ta1_z_yyy_yyyzzz_0[i] * pa_x[i] - ta1_z_yyy_yyyzzz_1[i] * pc_x[i];

        ta1_z_xyyy_yyzzzz_0[i] = ta1_z_yyy_yyzzzz_0[i] * pa_x[i] - ta1_z_yyy_yyzzzz_1[i] * pc_x[i];

        ta1_z_xyyy_yzzzzz_0[i] = ta1_z_yyy_yzzzzz_0[i] * pa_x[i] - ta1_z_yyy_yzzzzz_1[i] * pc_x[i];

        ta1_z_xyyy_zzzzzz_0[i] = ta1_z_yyy_zzzzzz_0[i] * pa_x[i] - ta1_z_yyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 1036-1064 components of targeted buffer : GI

    auto ta1_z_xyyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1036);

    auto ta1_z_xyyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1037);

    auto ta1_z_xyyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1038);

    auto ta1_z_xyyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1039);

    auto ta1_z_xyyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1040);

    auto ta1_z_xyyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1041);

    auto ta1_z_xyyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1042);

    auto ta1_z_xyyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1043);

    auto ta1_z_xyyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1044);

    auto ta1_z_xyyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1045);

    auto ta1_z_xyyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1046);

    auto ta1_z_xyyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1047);

    auto ta1_z_xyyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1048);

    auto ta1_z_xyyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1049);

    auto ta1_z_xyyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1050);

    auto ta1_z_xyyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1051);

    auto ta1_z_xyyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1052);

    auto ta1_z_xyyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1053);

    auto ta1_z_xyyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1054);

    auto ta1_z_xyyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1055);

    auto ta1_z_xyyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1056);

    auto ta1_z_xyyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1057);

    auto ta1_z_xyyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1058);

    auto ta1_z_xyyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1059);

    auto ta1_z_xyyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1060);

    auto ta1_z_xyyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1061);

    auto ta1_z_xyyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1062);

    auto ta1_z_xyyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1063);

#pragma omp simd aligned(pa_x,                    \
                             pa_z,                \
                             pc_x,                \
                             pc_z,                \
                             ta1_z_xyy_xxxxxx_0,  \
                             ta1_z_xyy_xxxxxx_1,  \
                             ta1_z_xyy_xxxxxy_0,  \
                             ta1_z_xyy_xxxxxy_1,  \
                             ta1_z_xyy_xxxxyy_0,  \
                             ta1_z_xyy_xxxxyy_1,  \
                             ta1_z_xyy_xxxyyy_0,  \
                             ta1_z_xyy_xxxyyy_1,  \
                             ta1_z_xyy_xxyyyy_0,  \
                             ta1_z_xyy_xxyyyy_1,  \
                             ta1_z_xyy_xyyyyy_0,  \
                             ta1_z_xyy_xyyyyy_1,  \
                             ta1_z_xyyz_xxxxxx_0, \
                             ta1_z_xyyz_xxxxxy_0, \
                             ta1_z_xyyz_xxxxxz_0, \
                             ta1_z_xyyz_xxxxyy_0, \
                             ta1_z_xyyz_xxxxyz_0, \
                             ta1_z_xyyz_xxxxzz_0, \
                             ta1_z_xyyz_xxxyyy_0, \
                             ta1_z_xyyz_xxxyyz_0, \
                             ta1_z_xyyz_xxxyzz_0, \
                             ta1_z_xyyz_xxxzzz_0, \
                             ta1_z_xyyz_xxyyyy_0, \
                             ta1_z_xyyz_xxyyyz_0, \
                             ta1_z_xyyz_xxyyzz_0, \
                             ta1_z_xyyz_xxyzzz_0, \
                             ta1_z_xyyz_xxzzzz_0, \
                             ta1_z_xyyz_xyyyyy_0, \
                             ta1_z_xyyz_xyyyyz_0, \
                             ta1_z_xyyz_xyyyzz_0, \
                             ta1_z_xyyz_xyyzzz_0, \
                             ta1_z_xyyz_xyzzzz_0, \
                             ta1_z_xyyz_xzzzzz_0, \
                             ta1_z_xyyz_yyyyyy_0, \
                             ta1_z_xyyz_yyyyyz_0, \
                             ta1_z_xyyz_yyyyzz_0, \
                             ta1_z_xyyz_yyyzzz_0, \
                             ta1_z_xyyz_yyzzzz_0, \
                             ta1_z_xyyz_yzzzzz_0, \
                             ta1_z_xyyz_zzzzzz_0, \
                             ta1_z_yyz_xxxxxz_0,  \
                             ta1_z_yyz_xxxxxz_1,  \
                             ta1_z_yyz_xxxxyz_0,  \
                             ta1_z_yyz_xxxxyz_1,  \
                             ta1_z_yyz_xxxxz_0,   \
                             ta1_z_yyz_xxxxz_1,   \
                             ta1_z_yyz_xxxxzz_0,  \
                             ta1_z_yyz_xxxxzz_1,  \
                             ta1_z_yyz_xxxyyz_0,  \
                             ta1_z_yyz_xxxyyz_1,  \
                             ta1_z_yyz_xxxyz_0,   \
                             ta1_z_yyz_xxxyz_1,   \
                             ta1_z_yyz_xxxyzz_0,  \
                             ta1_z_yyz_xxxyzz_1,  \
                             ta1_z_yyz_xxxzz_0,   \
                             ta1_z_yyz_xxxzz_1,   \
                             ta1_z_yyz_xxxzzz_0,  \
                             ta1_z_yyz_xxxzzz_1,  \
                             ta1_z_yyz_xxyyyz_0,  \
                             ta1_z_yyz_xxyyyz_1,  \
                             ta1_z_yyz_xxyyz_0,   \
                             ta1_z_yyz_xxyyz_1,   \
                             ta1_z_yyz_xxyyzz_0,  \
                             ta1_z_yyz_xxyyzz_1,  \
                             ta1_z_yyz_xxyzz_0,   \
                             ta1_z_yyz_xxyzz_1,   \
                             ta1_z_yyz_xxyzzz_0,  \
                             ta1_z_yyz_xxyzzz_1,  \
                             ta1_z_yyz_xxzzz_0,   \
                             ta1_z_yyz_xxzzz_1,   \
                             ta1_z_yyz_xxzzzz_0,  \
                             ta1_z_yyz_xxzzzz_1,  \
                             ta1_z_yyz_xyyyyz_0,  \
                             ta1_z_yyz_xyyyyz_1,  \
                             ta1_z_yyz_xyyyz_0,   \
                             ta1_z_yyz_xyyyz_1,   \
                             ta1_z_yyz_xyyyzz_0,  \
                             ta1_z_yyz_xyyyzz_1,  \
                             ta1_z_yyz_xyyzz_0,   \
                             ta1_z_yyz_xyyzz_1,   \
                             ta1_z_yyz_xyyzzz_0,  \
                             ta1_z_yyz_xyyzzz_1,  \
                             ta1_z_yyz_xyzzz_0,   \
                             ta1_z_yyz_xyzzz_1,   \
                             ta1_z_yyz_xyzzzz_0,  \
                             ta1_z_yyz_xyzzzz_1,  \
                             ta1_z_yyz_xzzzz_0,   \
                             ta1_z_yyz_xzzzz_1,   \
                             ta1_z_yyz_xzzzzz_0,  \
                             ta1_z_yyz_xzzzzz_1,  \
                             ta1_z_yyz_yyyyyy_0,  \
                             ta1_z_yyz_yyyyyy_1,  \
                             ta1_z_yyz_yyyyyz_0,  \
                             ta1_z_yyz_yyyyyz_1,  \
                             ta1_z_yyz_yyyyz_0,   \
                             ta1_z_yyz_yyyyz_1,   \
                             ta1_z_yyz_yyyyzz_0,  \
                             ta1_z_yyz_yyyyzz_1,  \
                             ta1_z_yyz_yyyzz_0,   \
                             ta1_z_yyz_yyyzz_1,   \
                             ta1_z_yyz_yyyzzz_0,  \
                             ta1_z_yyz_yyyzzz_1,  \
                             ta1_z_yyz_yyzzz_0,   \
                             ta1_z_yyz_yyzzz_1,   \
                             ta1_z_yyz_yyzzzz_0,  \
                             ta1_z_yyz_yyzzzz_1,  \
                             ta1_z_yyz_yzzzz_0,   \
                             ta1_z_yyz_yzzzz_1,   \
                             ta1_z_yyz_yzzzzz_0,  \
                             ta1_z_yyz_yzzzzz_1,  \
                             ta1_z_yyz_zzzzz_0,   \
                             ta1_z_yyz_zzzzz_1,   \
                             ta1_z_yyz_zzzzzz_0,  \
                             ta1_z_yyz_zzzzzz_1,  \
                             ta_xyy_xxxxxx_1,     \
                             ta_xyy_xxxxxy_1,     \
                             ta_xyy_xxxxyy_1,     \
                             ta_xyy_xxxyyy_1,     \
                             ta_xyy_xxyyyy_1,     \
                             ta_xyy_xyyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyz_xxxxxx_0[i] = ta_xyy_xxxxxx_1[i] + ta1_z_xyy_xxxxxx_0[i] * pa_z[i] - ta1_z_xyy_xxxxxx_1[i] * pc_z[i];

        ta1_z_xyyz_xxxxxy_0[i] = ta_xyy_xxxxxy_1[i] + ta1_z_xyy_xxxxxy_0[i] * pa_z[i] - ta1_z_xyy_xxxxxy_1[i] * pc_z[i];

        ta1_z_xyyz_xxxxxz_0[i] =
            5.0 * ta1_z_yyz_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_yyz_xxxxz_1[i] * fe_0 + ta1_z_yyz_xxxxxz_0[i] * pa_x[i] - ta1_z_yyz_xxxxxz_1[i] * pc_x[i];

        ta1_z_xyyz_xxxxyy_0[i] = ta_xyy_xxxxyy_1[i] + ta1_z_xyy_xxxxyy_0[i] * pa_z[i] - ta1_z_xyy_xxxxyy_1[i] * pc_z[i];

        ta1_z_xyyz_xxxxyz_0[i] =
            4.0 * ta1_z_yyz_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_yyz_xxxyz_1[i] * fe_0 + ta1_z_yyz_xxxxyz_0[i] * pa_x[i] - ta1_z_yyz_xxxxyz_1[i] * pc_x[i];

        ta1_z_xyyz_xxxxzz_0[i] =
            4.0 * ta1_z_yyz_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_yyz_xxxzz_1[i] * fe_0 + ta1_z_yyz_xxxxzz_0[i] * pa_x[i] - ta1_z_yyz_xxxxzz_1[i] * pc_x[i];

        ta1_z_xyyz_xxxyyy_0[i] = ta_xyy_xxxyyy_1[i] + ta1_z_xyy_xxxyyy_0[i] * pa_z[i] - ta1_z_xyy_xxxyyy_1[i] * pc_z[i];

        ta1_z_xyyz_xxxyyz_0[i] =
            3.0 * ta1_z_yyz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxyyz_1[i] * fe_0 + ta1_z_yyz_xxxyyz_0[i] * pa_x[i] - ta1_z_yyz_xxxyyz_1[i] * pc_x[i];

        ta1_z_xyyz_xxxyzz_0[i] =
            3.0 * ta1_z_yyz_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxyzz_1[i] * fe_0 + ta1_z_yyz_xxxyzz_0[i] * pa_x[i] - ta1_z_yyz_xxxyzz_1[i] * pc_x[i];

        ta1_z_xyyz_xxxzzz_0[i] =
            3.0 * ta1_z_yyz_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxzzz_1[i] * fe_0 + ta1_z_yyz_xxxzzz_0[i] * pa_x[i] - ta1_z_yyz_xxxzzz_1[i] * pc_x[i];

        ta1_z_xyyz_xxyyyy_0[i] = ta_xyy_xxyyyy_1[i] + ta1_z_xyy_xxyyyy_0[i] * pa_z[i] - ta1_z_xyy_xxyyyy_1[i] * pc_z[i];

        ta1_z_xyyz_xxyyyz_0[i] =
            2.0 * ta1_z_yyz_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyyyz_1[i] * fe_0 + ta1_z_yyz_xxyyyz_0[i] * pa_x[i] - ta1_z_yyz_xxyyyz_1[i] * pc_x[i];

        ta1_z_xyyz_xxyyzz_0[i] =
            2.0 * ta1_z_yyz_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyyzz_1[i] * fe_0 + ta1_z_yyz_xxyyzz_0[i] * pa_x[i] - ta1_z_yyz_xxyyzz_1[i] * pc_x[i];

        ta1_z_xyyz_xxyzzz_0[i] =
            2.0 * ta1_z_yyz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyzzz_1[i] * fe_0 + ta1_z_yyz_xxyzzz_0[i] * pa_x[i] - ta1_z_yyz_xxyzzz_1[i] * pc_x[i];

        ta1_z_xyyz_xxzzzz_0[i] =
            2.0 * ta1_z_yyz_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xzzzz_1[i] * fe_0 + ta1_z_yyz_xxzzzz_0[i] * pa_x[i] - ta1_z_yyz_xxzzzz_1[i] * pc_x[i];

        ta1_z_xyyz_xyyyyy_0[i] = ta_xyy_xyyyyy_1[i] + ta1_z_xyy_xyyyyy_0[i] * pa_z[i] - ta1_z_xyy_xyyyyy_1[i] * pc_z[i];

        ta1_z_xyyz_xyyyyz_0[i] =
            ta1_z_yyz_yyyyz_0[i] * fe_0 - ta1_z_yyz_yyyyz_1[i] * fe_0 + ta1_z_yyz_xyyyyz_0[i] * pa_x[i] - ta1_z_yyz_xyyyyz_1[i] * pc_x[i];

        ta1_z_xyyz_xyyyzz_0[i] =
            ta1_z_yyz_yyyzz_0[i] * fe_0 - ta1_z_yyz_yyyzz_1[i] * fe_0 + ta1_z_yyz_xyyyzz_0[i] * pa_x[i] - ta1_z_yyz_xyyyzz_1[i] * pc_x[i];

        ta1_z_xyyz_xyyzzz_0[i] =
            ta1_z_yyz_yyzzz_0[i] * fe_0 - ta1_z_yyz_yyzzz_1[i] * fe_0 + ta1_z_yyz_xyyzzz_0[i] * pa_x[i] - ta1_z_yyz_xyyzzz_1[i] * pc_x[i];

        ta1_z_xyyz_xyzzzz_0[i] =
            ta1_z_yyz_yzzzz_0[i] * fe_0 - ta1_z_yyz_yzzzz_1[i] * fe_0 + ta1_z_yyz_xyzzzz_0[i] * pa_x[i] - ta1_z_yyz_xyzzzz_1[i] * pc_x[i];

        ta1_z_xyyz_xzzzzz_0[i] =
            ta1_z_yyz_zzzzz_0[i] * fe_0 - ta1_z_yyz_zzzzz_1[i] * fe_0 + ta1_z_yyz_xzzzzz_0[i] * pa_x[i] - ta1_z_yyz_xzzzzz_1[i] * pc_x[i];

        ta1_z_xyyz_yyyyyy_0[i] = ta1_z_yyz_yyyyyy_0[i] * pa_x[i] - ta1_z_yyz_yyyyyy_1[i] * pc_x[i];

        ta1_z_xyyz_yyyyyz_0[i] = ta1_z_yyz_yyyyyz_0[i] * pa_x[i] - ta1_z_yyz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xyyz_yyyyzz_0[i] = ta1_z_yyz_yyyyzz_0[i] * pa_x[i] - ta1_z_yyz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xyyz_yyyzzz_0[i] = ta1_z_yyz_yyyzzz_0[i] * pa_x[i] - ta1_z_yyz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xyyz_yyzzzz_0[i] = ta1_z_yyz_yyzzzz_0[i] * pa_x[i] - ta1_z_yyz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xyyz_yzzzzz_0[i] = ta1_z_yyz_yzzzzz_0[i] * pa_x[i] - ta1_z_yyz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xyyz_zzzzzz_0[i] = ta1_z_yyz_zzzzzz_0[i] * pa_x[i] - ta1_z_yyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 1064-1092 components of targeted buffer : GI

    auto ta1_z_xyzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1064);

    auto ta1_z_xyzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1065);

    auto ta1_z_xyzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1066);

    auto ta1_z_xyzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1067);

    auto ta1_z_xyzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1068);

    auto ta1_z_xyzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1069);

    auto ta1_z_xyzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1070);

    auto ta1_z_xyzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1071);

    auto ta1_z_xyzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1072);

    auto ta1_z_xyzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1073);

    auto ta1_z_xyzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1074);

    auto ta1_z_xyzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1075);

    auto ta1_z_xyzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1076);

    auto ta1_z_xyzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1077);

    auto ta1_z_xyzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1078);

    auto ta1_z_xyzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1079);

    auto ta1_z_xyzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1080);

    auto ta1_z_xyzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1081);

    auto ta1_z_xyzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1082);

    auto ta1_z_xyzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1083);

    auto ta1_z_xyzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1084);

    auto ta1_z_xyzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1085);

    auto ta1_z_xyzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1086);

    auto ta1_z_xyzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1087);

    auto ta1_z_xyzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1088);

    auto ta1_z_xyzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1089);

    auto ta1_z_xyzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1090);

    auto ta1_z_xyzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1091);

#pragma omp simd aligned(pa_x,                    \
                             pa_y,                \
                             pc_x,                \
                             pc_y,                \
                             ta1_z_xyzz_xxxxxx_0, \
                             ta1_z_xyzz_xxxxxy_0, \
                             ta1_z_xyzz_xxxxxz_0, \
                             ta1_z_xyzz_xxxxyy_0, \
                             ta1_z_xyzz_xxxxyz_0, \
                             ta1_z_xyzz_xxxxzz_0, \
                             ta1_z_xyzz_xxxyyy_0, \
                             ta1_z_xyzz_xxxyyz_0, \
                             ta1_z_xyzz_xxxyzz_0, \
                             ta1_z_xyzz_xxxzzz_0, \
                             ta1_z_xyzz_xxyyyy_0, \
                             ta1_z_xyzz_xxyyyz_0, \
                             ta1_z_xyzz_xxyyzz_0, \
                             ta1_z_xyzz_xxyzzz_0, \
                             ta1_z_xyzz_xxzzzz_0, \
                             ta1_z_xyzz_xyyyyy_0, \
                             ta1_z_xyzz_xyyyyz_0, \
                             ta1_z_xyzz_xyyyzz_0, \
                             ta1_z_xyzz_xyyzzz_0, \
                             ta1_z_xyzz_xyzzzz_0, \
                             ta1_z_xyzz_xzzzzz_0, \
                             ta1_z_xyzz_yyyyyy_0, \
                             ta1_z_xyzz_yyyyyz_0, \
                             ta1_z_xyzz_yyyyzz_0, \
                             ta1_z_xyzz_yyyzzz_0, \
                             ta1_z_xyzz_yyzzzz_0, \
                             ta1_z_xyzz_yzzzzz_0, \
                             ta1_z_xyzz_zzzzzz_0, \
                             ta1_z_xzz_xxxxxx_0,  \
                             ta1_z_xzz_xxxxxx_1,  \
                             ta1_z_xzz_xxxxxz_0,  \
                             ta1_z_xzz_xxxxxz_1,  \
                             ta1_z_xzz_xxxxzz_0,  \
                             ta1_z_xzz_xxxxzz_1,  \
                             ta1_z_xzz_xxxzzz_0,  \
                             ta1_z_xzz_xxxzzz_1,  \
                             ta1_z_xzz_xxzzzz_0,  \
                             ta1_z_xzz_xxzzzz_1,  \
                             ta1_z_xzz_xzzzzz_0,  \
                             ta1_z_xzz_xzzzzz_1,  \
                             ta1_z_yzz_xxxxxy_0,  \
                             ta1_z_yzz_xxxxxy_1,  \
                             ta1_z_yzz_xxxxy_0,   \
                             ta1_z_yzz_xxxxy_1,   \
                             ta1_z_yzz_xxxxyy_0,  \
                             ta1_z_yzz_xxxxyy_1,  \
                             ta1_z_yzz_xxxxyz_0,  \
                             ta1_z_yzz_xxxxyz_1,  \
                             ta1_z_yzz_xxxyy_0,   \
                             ta1_z_yzz_xxxyy_1,   \
                             ta1_z_yzz_xxxyyy_0,  \
                             ta1_z_yzz_xxxyyy_1,  \
                             ta1_z_yzz_xxxyyz_0,  \
                             ta1_z_yzz_xxxyyz_1,  \
                             ta1_z_yzz_xxxyz_0,   \
                             ta1_z_yzz_xxxyz_1,   \
                             ta1_z_yzz_xxxyzz_0,  \
                             ta1_z_yzz_xxxyzz_1,  \
                             ta1_z_yzz_xxyyy_0,   \
                             ta1_z_yzz_xxyyy_1,   \
                             ta1_z_yzz_xxyyyy_0,  \
                             ta1_z_yzz_xxyyyy_1,  \
                             ta1_z_yzz_xxyyyz_0,  \
                             ta1_z_yzz_xxyyyz_1,  \
                             ta1_z_yzz_xxyyz_0,   \
                             ta1_z_yzz_xxyyz_1,   \
                             ta1_z_yzz_xxyyzz_0,  \
                             ta1_z_yzz_xxyyzz_1,  \
                             ta1_z_yzz_xxyzz_0,   \
                             ta1_z_yzz_xxyzz_1,   \
                             ta1_z_yzz_xxyzzz_0,  \
                             ta1_z_yzz_xxyzzz_1,  \
                             ta1_z_yzz_xyyyy_0,   \
                             ta1_z_yzz_xyyyy_1,   \
                             ta1_z_yzz_xyyyyy_0,  \
                             ta1_z_yzz_xyyyyy_1,  \
                             ta1_z_yzz_xyyyyz_0,  \
                             ta1_z_yzz_xyyyyz_1,  \
                             ta1_z_yzz_xyyyz_0,   \
                             ta1_z_yzz_xyyyz_1,   \
                             ta1_z_yzz_xyyyzz_0,  \
                             ta1_z_yzz_xyyyzz_1,  \
                             ta1_z_yzz_xyyzz_0,   \
                             ta1_z_yzz_xyyzz_1,   \
                             ta1_z_yzz_xyyzzz_0,  \
                             ta1_z_yzz_xyyzzz_1,  \
                             ta1_z_yzz_xyzzz_0,   \
                             ta1_z_yzz_xyzzz_1,   \
                             ta1_z_yzz_xyzzzz_0,  \
                             ta1_z_yzz_xyzzzz_1,  \
                             ta1_z_yzz_yyyyy_0,   \
                             ta1_z_yzz_yyyyy_1,   \
                             ta1_z_yzz_yyyyyy_0,  \
                             ta1_z_yzz_yyyyyy_1,  \
                             ta1_z_yzz_yyyyyz_0,  \
                             ta1_z_yzz_yyyyyz_1,  \
                             ta1_z_yzz_yyyyz_0,   \
                             ta1_z_yzz_yyyyz_1,   \
                             ta1_z_yzz_yyyyzz_0,  \
                             ta1_z_yzz_yyyyzz_1,  \
                             ta1_z_yzz_yyyzz_0,   \
                             ta1_z_yzz_yyyzz_1,   \
                             ta1_z_yzz_yyyzzz_0,  \
                             ta1_z_yzz_yyyzzz_1,  \
                             ta1_z_yzz_yyzzz_0,   \
                             ta1_z_yzz_yyzzz_1,   \
                             ta1_z_yzz_yyzzzz_0,  \
                             ta1_z_yzz_yyzzzz_1,  \
                             ta1_z_yzz_yzzzz_0,   \
                             ta1_z_yzz_yzzzz_1,   \
                             ta1_z_yzz_yzzzzz_0,  \
                             ta1_z_yzz_yzzzzz_1,  \
                             ta1_z_yzz_zzzzzz_0,  \
                             ta1_z_yzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzz_xxxxxx_0[i] = ta1_z_xzz_xxxxxx_0[i] * pa_y[i] - ta1_z_xzz_xxxxxx_1[i] * pc_y[i];

        ta1_z_xyzz_xxxxxy_0[i] =
            5.0 * ta1_z_yzz_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_yzz_xxxxy_1[i] * fe_0 + ta1_z_yzz_xxxxxy_0[i] * pa_x[i] - ta1_z_yzz_xxxxxy_1[i] * pc_x[i];

        ta1_z_xyzz_xxxxxz_0[i] = ta1_z_xzz_xxxxxz_0[i] * pa_y[i] - ta1_z_xzz_xxxxxz_1[i] * pc_y[i];

        ta1_z_xyzz_xxxxyy_0[i] =
            4.0 * ta1_z_yzz_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_yzz_xxxyy_1[i] * fe_0 + ta1_z_yzz_xxxxyy_0[i] * pa_x[i] - ta1_z_yzz_xxxxyy_1[i] * pc_x[i];

        ta1_z_xyzz_xxxxyz_0[i] =
            4.0 * ta1_z_yzz_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_yzz_xxxyz_1[i] * fe_0 + ta1_z_yzz_xxxxyz_0[i] * pa_x[i] - ta1_z_yzz_xxxxyz_1[i] * pc_x[i];

        ta1_z_xyzz_xxxxzz_0[i] = ta1_z_xzz_xxxxzz_0[i] * pa_y[i] - ta1_z_xzz_xxxxzz_1[i] * pc_y[i];

        ta1_z_xyzz_xxxyyy_0[i] =
            3.0 * ta1_z_yzz_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_yzz_xxyyy_1[i] * fe_0 + ta1_z_yzz_xxxyyy_0[i] * pa_x[i] - ta1_z_yzz_xxxyyy_1[i] * pc_x[i];

        ta1_z_xyzz_xxxyyz_0[i] =
            3.0 * ta1_z_yzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_yzz_xxyyz_1[i] * fe_0 + ta1_z_yzz_xxxyyz_0[i] * pa_x[i] - ta1_z_yzz_xxxyyz_1[i] * pc_x[i];

        ta1_z_xyzz_xxxyzz_0[i] =
            3.0 * ta1_z_yzz_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yzz_xxyzz_1[i] * fe_0 + ta1_z_yzz_xxxyzz_0[i] * pa_x[i] - ta1_z_yzz_xxxyzz_1[i] * pc_x[i];

        ta1_z_xyzz_xxxzzz_0[i] = ta1_z_xzz_xxxzzz_0[i] * pa_y[i] - ta1_z_xzz_xxxzzz_1[i] * pc_y[i];

        ta1_z_xyzz_xxyyyy_0[i] =
            2.0 * ta1_z_yzz_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyyyy_1[i] * fe_0 + ta1_z_yzz_xxyyyy_0[i] * pa_x[i] - ta1_z_yzz_xxyyyy_1[i] * pc_x[i];

        ta1_z_xyzz_xxyyyz_0[i] =
            2.0 * ta1_z_yzz_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyyyz_1[i] * fe_0 + ta1_z_yzz_xxyyyz_0[i] * pa_x[i] - ta1_z_yzz_xxyyyz_1[i] * pc_x[i];

        ta1_z_xyzz_xxyyzz_0[i] =
            2.0 * ta1_z_yzz_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyyzz_1[i] * fe_0 + ta1_z_yzz_xxyyzz_0[i] * pa_x[i] - ta1_z_yzz_xxyyzz_1[i] * pc_x[i];

        ta1_z_xyzz_xxyzzz_0[i] =
            2.0 * ta1_z_yzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyzzz_1[i] * fe_0 + ta1_z_yzz_xxyzzz_0[i] * pa_x[i] - ta1_z_yzz_xxyzzz_1[i] * pc_x[i];

        ta1_z_xyzz_xxzzzz_0[i] = ta1_z_xzz_xxzzzz_0[i] * pa_y[i] - ta1_z_xzz_xxzzzz_1[i] * pc_y[i];

        ta1_z_xyzz_xyyyyy_0[i] =
            ta1_z_yzz_yyyyy_0[i] * fe_0 - ta1_z_yzz_yyyyy_1[i] * fe_0 + ta1_z_yzz_xyyyyy_0[i] * pa_x[i] - ta1_z_yzz_xyyyyy_1[i] * pc_x[i];

        ta1_z_xyzz_xyyyyz_0[i] =
            ta1_z_yzz_yyyyz_0[i] * fe_0 - ta1_z_yzz_yyyyz_1[i] * fe_0 + ta1_z_yzz_xyyyyz_0[i] * pa_x[i] - ta1_z_yzz_xyyyyz_1[i] * pc_x[i];

        ta1_z_xyzz_xyyyzz_0[i] =
            ta1_z_yzz_yyyzz_0[i] * fe_0 - ta1_z_yzz_yyyzz_1[i] * fe_0 + ta1_z_yzz_xyyyzz_0[i] * pa_x[i] - ta1_z_yzz_xyyyzz_1[i] * pc_x[i];

        ta1_z_xyzz_xyyzzz_0[i] =
            ta1_z_yzz_yyzzz_0[i] * fe_0 - ta1_z_yzz_yyzzz_1[i] * fe_0 + ta1_z_yzz_xyyzzz_0[i] * pa_x[i] - ta1_z_yzz_xyyzzz_1[i] * pc_x[i];

        ta1_z_xyzz_xyzzzz_0[i] =
            ta1_z_yzz_yzzzz_0[i] * fe_0 - ta1_z_yzz_yzzzz_1[i] * fe_0 + ta1_z_yzz_xyzzzz_0[i] * pa_x[i] - ta1_z_yzz_xyzzzz_1[i] * pc_x[i];

        ta1_z_xyzz_xzzzzz_0[i] = ta1_z_xzz_xzzzzz_0[i] * pa_y[i] - ta1_z_xzz_xzzzzz_1[i] * pc_y[i];

        ta1_z_xyzz_yyyyyy_0[i] = ta1_z_yzz_yyyyyy_0[i] * pa_x[i] - ta1_z_yzz_yyyyyy_1[i] * pc_x[i];

        ta1_z_xyzz_yyyyyz_0[i] = ta1_z_yzz_yyyyyz_0[i] * pa_x[i] - ta1_z_yzz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xyzz_yyyyzz_0[i] = ta1_z_yzz_yyyyzz_0[i] * pa_x[i] - ta1_z_yzz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xyzz_yyyzzz_0[i] = ta1_z_yzz_yyyzzz_0[i] * pa_x[i] - ta1_z_yzz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xyzz_yyzzzz_0[i] = ta1_z_yzz_yyzzzz_0[i] * pa_x[i] - ta1_z_yzz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xyzz_yzzzzz_0[i] = ta1_z_yzz_yzzzzz_0[i] * pa_x[i] - ta1_z_yzz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xyzz_zzzzzz_0[i] = ta1_z_yzz_zzzzzz_0[i] * pa_x[i] - ta1_z_yzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 1092-1120 components of targeted buffer : GI

    auto ta1_z_xzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1092);

    auto ta1_z_xzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1093);

    auto ta1_z_xzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1094);

    auto ta1_z_xzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1095);

    auto ta1_z_xzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1096);

    auto ta1_z_xzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1097);

    auto ta1_z_xzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1098);

    auto ta1_z_xzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1099);

    auto ta1_z_xzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1100);

    auto ta1_z_xzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1101);

    auto ta1_z_xzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1102);

    auto ta1_z_xzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1103);

    auto ta1_z_xzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1104);

    auto ta1_z_xzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1105);

    auto ta1_z_xzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1106);

    auto ta1_z_xzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1107);

    auto ta1_z_xzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1108);

    auto ta1_z_xzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1109);

    auto ta1_z_xzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1110);

    auto ta1_z_xzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1111);

    auto ta1_z_xzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1112);

    auto ta1_z_xzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1113);

    auto ta1_z_xzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1114);

    auto ta1_z_xzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1115);

    auto ta1_z_xzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1116);

    auto ta1_z_xzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1117);

    auto ta1_z_xzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1118);

    auto ta1_z_xzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1119);

#pragma omp simd aligned(pa_x,                    \
                             pc_x,                \
                             ta1_z_xzzz_xxxxxx_0, \
                             ta1_z_xzzz_xxxxxy_0, \
                             ta1_z_xzzz_xxxxxz_0, \
                             ta1_z_xzzz_xxxxyy_0, \
                             ta1_z_xzzz_xxxxyz_0, \
                             ta1_z_xzzz_xxxxzz_0, \
                             ta1_z_xzzz_xxxyyy_0, \
                             ta1_z_xzzz_xxxyyz_0, \
                             ta1_z_xzzz_xxxyzz_0, \
                             ta1_z_xzzz_xxxzzz_0, \
                             ta1_z_xzzz_xxyyyy_0, \
                             ta1_z_xzzz_xxyyyz_0, \
                             ta1_z_xzzz_xxyyzz_0, \
                             ta1_z_xzzz_xxyzzz_0, \
                             ta1_z_xzzz_xxzzzz_0, \
                             ta1_z_xzzz_xyyyyy_0, \
                             ta1_z_xzzz_xyyyyz_0, \
                             ta1_z_xzzz_xyyyzz_0, \
                             ta1_z_xzzz_xyyzzz_0, \
                             ta1_z_xzzz_xyzzzz_0, \
                             ta1_z_xzzz_xzzzzz_0, \
                             ta1_z_xzzz_yyyyyy_0, \
                             ta1_z_xzzz_yyyyyz_0, \
                             ta1_z_xzzz_yyyyzz_0, \
                             ta1_z_xzzz_yyyzzz_0, \
                             ta1_z_xzzz_yyzzzz_0, \
                             ta1_z_xzzz_yzzzzz_0, \
                             ta1_z_xzzz_zzzzzz_0, \
                             ta1_z_zzz_xxxxx_0,   \
                             ta1_z_zzz_xxxxx_1,   \
                             ta1_z_zzz_xxxxxx_0,  \
                             ta1_z_zzz_xxxxxx_1,  \
                             ta1_z_zzz_xxxxxy_0,  \
                             ta1_z_zzz_xxxxxy_1,  \
                             ta1_z_zzz_xxxxxz_0,  \
                             ta1_z_zzz_xxxxxz_1,  \
                             ta1_z_zzz_xxxxy_0,   \
                             ta1_z_zzz_xxxxy_1,   \
                             ta1_z_zzz_xxxxyy_0,  \
                             ta1_z_zzz_xxxxyy_1,  \
                             ta1_z_zzz_xxxxyz_0,  \
                             ta1_z_zzz_xxxxyz_1,  \
                             ta1_z_zzz_xxxxz_0,   \
                             ta1_z_zzz_xxxxz_1,   \
                             ta1_z_zzz_xxxxzz_0,  \
                             ta1_z_zzz_xxxxzz_1,  \
                             ta1_z_zzz_xxxyy_0,   \
                             ta1_z_zzz_xxxyy_1,   \
                             ta1_z_zzz_xxxyyy_0,  \
                             ta1_z_zzz_xxxyyy_1,  \
                             ta1_z_zzz_xxxyyz_0,  \
                             ta1_z_zzz_xxxyyz_1,  \
                             ta1_z_zzz_xxxyz_0,   \
                             ta1_z_zzz_xxxyz_1,   \
                             ta1_z_zzz_xxxyzz_0,  \
                             ta1_z_zzz_xxxyzz_1,  \
                             ta1_z_zzz_xxxzz_0,   \
                             ta1_z_zzz_xxxzz_1,   \
                             ta1_z_zzz_xxxzzz_0,  \
                             ta1_z_zzz_xxxzzz_1,  \
                             ta1_z_zzz_xxyyy_0,   \
                             ta1_z_zzz_xxyyy_1,   \
                             ta1_z_zzz_xxyyyy_0,  \
                             ta1_z_zzz_xxyyyy_1,  \
                             ta1_z_zzz_xxyyyz_0,  \
                             ta1_z_zzz_xxyyyz_1,  \
                             ta1_z_zzz_xxyyz_0,   \
                             ta1_z_zzz_xxyyz_1,   \
                             ta1_z_zzz_xxyyzz_0,  \
                             ta1_z_zzz_xxyyzz_1,  \
                             ta1_z_zzz_xxyzz_0,   \
                             ta1_z_zzz_xxyzz_1,   \
                             ta1_z_zzz_xxyzzz_0,  \
                             ta1_z_zzz_xxyzzz_1,  \
                             ta1_z_zzz_xxzzz_0,   \
                             ta1_z_zzz_xxzzz_1,   \
                             ta1_z_zzz_xxzzzz_0,  \
                             ta1_z_zzz_xxzzzz_1,  \
                             ta1_z_zzz_xyyyy_0,   \
                             ta1_z_zzz_xyyyy_1,   \
                             ta1_z_zzz_xyyyyy_0,  \
                             ta1_z_zzz_xyyyyy_1,  \
                             ta1_z_zzz_xyyyyz_0,  \
                             ta1_z_zzz_xyyyyz_1,  \
                             ta1_z_zzz_xyyyz_0,   \
                             ta1_z_zzz_xyyyz_1,   \
                             ta1_z_zzz_xyyyzz_0,  \
                             ta1_z_zzz_xyyyzz_1,  \
                             ta1_z_zzz_xyyzz_0,   \
                             ta1_z_zzz_xyyzz_1,   \
                             ta1_z_zzz_xyyzzz_0,  \
                             ta1_z_zzz_xyyzzz_1,  \
                             ta1_z_zzz_xyzzz_0,   \
                             ta1_z_zzz_xyzzz_1,   \
                             ta1_z_zzz_xyzzzz_0,  \
                             ta1_z_zzz_xyzzzz_1,  \
                             ta1_z_zzz_xzzzz_0,   \
                             ta1_z_zzz_xzzzz_1,   \
                             ta1_z_zzz_xzzzzz_0,  \
                             ta1_z_zzz_xzzzzz_1,  \
                             ta1_z_zzz_yyyyy_0,   \
                             ta1_z_zzz_yyyyy_1,   \
                             ta1_z_zzz_yyyyyy_0,  \
                             ta1_z_zzz_yyyyyy_1,  \
                             ta1_z_zzz_yyyyyz_0,  \
                             ta1_z_zzz_yyyyyz_1,  \
                             ta1_z_zzz_yyyyz_0,   \
                             ta1_z_zzz_yyyyz_1,   \
                             ta1_z_zzz_yyyyzz_0,  \
                             ta1_z_zzz_yyyyzz_1,  \
                             ta1_z_zzz_yyyzz_0,   \
                             ta1_z_zzz_yyyzz_1,   \
                             ta1_z_zzz_yyyzzz_0,  \
                             ta1_z_zzz_yyyzzz_1,  \
                             ta1_z_zzz_yyzzz_0,   \
                             ta1_z_zzz_yyzzz_1,   \
                             ta1_z_zzz_yyzzzz_0,  \
                             ta1_z_zzz_yyzzzz_1,  \
                             ta1_z_zzz_yzzzz_0,   \
                             ta1_z_zzz_yzzzz_1,   \
                             ta1_z_zzz_yzzzzz_0,  \
                             ta1_z_zzz_yzzzzz_1,  \
                             ta1_z_zzz_zzzzz_0,   \
                             ta1_z_zzz_zzzzz_1,   \
                             ta1_z_zzz_zzzzzz_0,  \
                             ta1_z_zzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzz_xxxxxx_0[i] =
            6.0 * ta1_z_zzz_xxxxx_0[i] * fe_0 - 6.0 * ta1_z_zzz_xxxxx_1[i] * fe_0 + ta1_z_zzz_xxxxxx_0[i] * pa_x[i] - ta1_z_zzz_xxxxxx_1[i] * pc_x[i];

        ta1_z_xzzz_xxxxxy_0[i] =
            5.0 * ta1_z_zzz_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_zzz_xxxxy_1[i] * fe_0 + ta1_z_zzz_xxxxxy_0[i] * pa_x[i] - ta1_z_zzz_xxxxxy_1[i] * pc_x[i];

        ta1_z_xzzz_xxxxxz_0[i] =
            5.0 * ta1_z_zzz_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_zzz_xxxxz_1[i] * fe_0 + ta1_z_zzz_xxxxxz_0[i] * pa_x[i] - ta1_z_zzz_xxxxxz_1[i] * pc_x[i];

        ta1_z_xzzz_xxxxyy_0[i] =
            4.0 * ta1_z_zzz_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxyy_1[i] * fe_0 + ta1_z_zzz_xxxxyy_0[i] * pa_x[i] - ta1_z_zzz_xxxxyy_1[i] * pc_x[i];

        ta1_z_xzzz_xxxxyz_0[i] =
            4.0 * ta1_z_zzz_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxyz_1[i] * fe_0 + ta1_z_zzz_xxxxyz_0[i] * pa_x[i] - ta1_z_zzz_xxxxyz_1[i] * pc_x[i];

        ta1_z_xzzz_xxxxzz_0[i] =
            4.0 * ta1_z_zzz_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxzz_1[i] * fe_0 + ta1_z_zzz_xxxxzz_0[i] * pa_x[i] - ta1_z_zzz_xxxxzz_1[i] * pc_x[i];

        ta1_z_xzzz_xxxyyy_0[i] =
            3.0 * ta1_z_zzz_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxyyy_1[i] * fe_0 + ta1_z_zzz_xxxyyy_0[i] * pa_x[i] - ta1_z_zzz_xxxyyy_1[i] * pc_x[i];

        ta1_z_xzzz_xxxyyz_0[i] =
            3.0 * ta1_z_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxyyz_1[i] * fe_0 + ta1_z_zzz_xxxyyz_0[i] * pa_x[i] - ta1_z_zzz_xxxyyz_1[i] * pc_x[i];

        ta1_z_xzzz_xxxyzz_0[i] =
            3.0 * ta1_z_zzz_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxyzz_1[i] * fe_0 + ta1_z_zzz_xxxyzz_0[i] * pa_x[i] - ta1_z_zzz_xxxyzz_1[i] * pc_x[i];

        ta1_z_xzzz_xxxzzz_0[i] =
            3.0 * ta1_z_zzz_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxzzz_1[i] * fe_0 + ta1_z_zzz_xxxzzz_0[i] * pa_x[i] - ta1_z_zzz_xxxzzz_1[i] * pc_x[i];

        ta1_z_xzzz_xxyyyy_0[i] =
            2.0 * ta1_z_zzz_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyyyy_1[i] * fe_0 + ta1_z_zzz_xxyyyy_0[i] * pa_x[i] - ta1_z_zzz_xxyyyy_1[i] * pc_x[i];

        ta1_z_xzzz_xxyyyz_0[i] =
            2.0 * ta1_z_zzz_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyyyz_1[i] * fe_0 + ta1_z_zzz_xxyyyz_0[i] * pa_x[i] - ta1_z_zzz_xxyyyz_1[i] * pc_x[i];

        ta1_z_xzzz_xxyyzz_0[i] =
            2.0 * ta1_z_zzz_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyyzz_1[i] * fe_0 + ta1_z_zzz_xxyyzz_0[i] * pa_x[i] - ta1_z_zzz_xxyyzz_1[i] * pc_x[i];

        ta1_z_xzzz_xxyzzz_0[i] =
            2.0 * ta1_z_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyzzz_1[i] * fe_0 + ta1_z_zzz_xxyzzz_0[i] * pa_x[i] - ta1_z_zzz_xxyzzz_1[i] * pc_x[i];

        ta1_z_xzzz_xxzzzz_0[i] =
            2.0 * ta1_z_zzz_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xzzzz_1[i] * fe_0 + ta1_z_zzz_xxzzzz_0[i] * pa_x[i] - ta1_z_zzz_xxzzzz_1[i] * pc_x[i];

        ta1_z_xzzz_xyyyyy_0[i] =
            ta1_z_zzz_yyyyy_0[i] * fe_0 - ta1_z_zzz_yyyyy_1[i] * fe_0 + ta1_z_zzz_xyyyyy_0[i] * pa_x[i] - ta1_z_zzz_xyyyyy_1[i] * pc_x[i];

        ta1_z_xzzz_xyyyyz_0[i] =
            ta1_z_zzz_yyyyz_0[i] * fe_0 - ta1_z_zzz_yyyyz_1[i] * fe_0 + ta1_z_zzz_xyyyyz_0[i] * pa_x[i] - ta1_z_zzz_xyyyyz_1[i] * pc_x[i];

        ta1_z_xzzz_xyyyzz_0[i] =
            ta1_z_zzz_yyyzz_0[i] * fe_0 - ta1_z_zzz_yyyzz_1[i] * fe_0 + ta1_z_zzz_xyyyzz_0[i] * pa_x[i] - ta1_z_zzz_xyyyzz_1[i] * pc_x[i];

        ta1_z_xzzz_xyyzzz_0[i] =
            ta1_z_zzz_yyzzz_0[i] * fe_0 - ta1_z_zzz_yyzzz_1[i] * fe_0 + ta1_z_zzz_xyyzzz_0[i] * pa_x[i] - ta1_z_zzz_xyyzzz_1[i] * pc_x[i];

        ta1_z_xzzz_xyzzzz_0[i] =
            ta1_z_zzz_yzzzz_0[i] * fe_0 - ta1_z_zzz_yzzzz_1[i] * fe_0 + ta1_z_zzz_xyzzzz_0[i] * pa_x[i] - ta1_z_zzz_xyzzzz_1[i] * pc_x[i];

        ta1_z_xzzz_xzzzzz_0[i] =
            ta1_z_zzz_zzzzz_0[i] * fe_0 - ta1_z_zzz_zzzzz_1[i] * fe_0 + ta1_z_zzz_xzzzzz_0[i] * pa_x[i] - ta1_z_zzz_xzzzzz_1[i] * pc_x[i];

        ta1_z_xzzz_yyyyyy_0[i] = ta1_z_zzz_yyyyyy_0[i] * pa_x[i] - ta1_z_zzz_yyyyyy_1[i] * pc_x[i];

        ta1_z_xzzz_yyyyyz_0[i] = ta1_z_zzz_yyyyyz_0[i] * pa_x[i] - ta1_z_zzz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xzzz_yyyyzz_0[i] = ta1_z_zzz_yyyyzz_0[i] * pa_x[i] - ta1_z_zzz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xzzz_yyyzzz_0[i] = ta1_z_zzz_yyyzzz_0[i] * pa_x[i] - ta1_z_zzz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xzzz_yyzzzz_0[i] = ta1_z_zzz_yyzzzz_0[i] * pa_x[i] - ta1_z_zzz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xzzz_yzzzzz_0[i] = ta1_z_zzz_yzzzzz_0[i] * pa_x[i] - ta1_z_zzz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xzzz_zzzzzz_0[i] = ta1_z_zzz_zzzzzz_0[i] * pa_x[i] - ta1_z_zzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 1120-1148 components of targeted buffer : GI

    auto ta1_z_yyyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1120);

    auto ta1_z_yyyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1121);

    auto ta1_z_yyyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1122);

    auto ta1_z_yyyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1123);

    auto ta1_z_yyyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1124);

    auto ta1_z_yyyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1125);

    auto ta1_z_yyyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1126);

    auto ta1_z_yyyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1127);

    auto ta1_z_yyyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1128);

    auto ta1_z_yyyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1129);

    auto ta1_z_yyyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1130);

    auto ta1_z_yyyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1131);

    auto ta1_z_yyyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1132);

    auto ta1_z_yyyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1133);

    auto ta1_z_yyyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1134);

    auto ta1_z_yyyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1135);

    auto ta1_z_yyyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1136);

    auto ta1_z_yyyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1137);

    auto ta1_z_yyyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1138);

    auto ta1_z_yyyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1139);

    auto ta1_z_yyyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1140);

    auto ta1_z_yyyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1141);

    auto ta1_z_yyyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1142);

    auto ta1_z_yyyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1143);

    auto ta1_z_yyyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1144);

    auto ta1_z_yyyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1145);

    auto ta1_z_yyyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1146);

    auto ta1_z_yyyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1147);

#pragma omp simd aligned(pa_y,                    \
                             pc_y,                \
                             ta1_z_yy_xxxxxx_0,   \
                             ta1_z_yy_xxxxxx_1,   \
                             ta1_z_yy_xxxxxy_0,   \
                             ta1_z_yy_xxxxxy_1,   \
                             ta1_z_yy_xxxxxz_0,   \
                             ta1_z_yy_xxxxxz_1,   \
                             ta1_z_yy_xxxxyy_0,   \
                             ta1_z_yy_xxxxyy_1,   \
                             ta1_z_yy_xxxxyz_0,   \
                             ta1_z_yy_xxxxyz_1,   \
                             ta1_z_yy_xxxxzz_0,   \
                             ta1_z_yy_xxxxzz_1,   \
                             ta1_z_yy_xxxyyy_0,   \
                             ta1_z_yy_xxxyyy_1,   \
                             ta1_z_yy_xxxyyz_0,   \
                             ta1_z_yy_xxxyyz_1,   \
                             ta1_z_yy_xxxyzz_0,   \
                             ta1_z_yy_xxxyzz_1,   \
                             ta1_z_yy_xxxzzz_0,   \
                             ta1_z_yy_xxxzzz_1,   \
                             ta1_z_yy_xxyyyy_0,   \
                             ta1_z_yy_xxyyyy_1,   \
                             ta1_z_yy_xxyyyz_0,   \
                             ta1_z_yy_xxyyyz_1,   \
                             ta1_z_yy_xxyyzz_0,   \
                             ta1_z_yy_xxyyzz_1,   \
                             ta1_z_yy_xxyzzz_0,   \
                             ta1_z_yy_xxyzzz_1,   \
                             ta1_z_yy_xxzzzz_0,   \
                             ta1_z_yy_xxzzzz_1,   \
                             ta1_z_yy_xyyyyy_0,   \
                             ta1_z_yy_xyyyyy_1,   \
                             ta1_z_yy_xyyyyz_0,   \
                             ta1_z_yy_xyyyyz_1,   \
                             ta1_z_yy_xyyyzz_0,   \
                             ta1_z_yy_xyyyzz_1,   \
                             ta1_z_yy_xyyzzz_0,   \
                             ta1_z_yy_xyyzzz_1,   \
                             ta1_z_yy_xyzzzz_0,   \
                             ta1_z_yy_xyzzzz_1,   \
                             ta1_z_yy_xzzzzz_0,   \
                             ta1_z_yy_xzzzzz_1,   \
                             ta1_z_yy_yyyyyy_0,   \
                             ta1_z_yy_yyyyyy_1,   \
                             ta1_z_yy_yyyyyz_0,   \
                             ta1_z_yy_yyyyyz_1,   \
                             ta1_z_yy_yyyyzz_0,   \
                             ta1_z_yy_yyyyzz_1,   \
                             ta1_z_yy_yyyzzz_0,   \
                             ta1_z_yy_yyyzzz_1,   \
                             ta1_z_yy_yyzzzz_0,   \
                             ta1_z_yy_yyzzzz_1,   \
                             ta1_z_yy_yzzzzz_0,   \
                             ta1_z_yy_yzzzzz_1,   \
                             ta1_z_yy_zzzzzz_0,   \
                             ta1_z_yy_zzzzzz_1,   \
                             ta1_z_yyy_xxxxx_0,   \
                             ta1_z_yyy_xxxxx_1,   \
                             ta1_z_yyy_xxxxxx_0,  \
                             ta1_z_yyy_xxxxxx_1,  \
                             ta1_z_yyy_xxxxxy_0,  \
                             ta1_z_yyy_xxxxxy_1,  \
                             ta1_z_yyy_xxxxxz_0,  \
                             ta1_z_yyy_xxxxxz_1,  \
                             ta1_z_yyy_xxxxy_0,   \
                             ta1_z_yyy_xxxxy_1,   \
                             ta1_z_yyy_xxxxyy_0,  \
                             ta1_z_yyy_xxxxyy_1,  \
                             ta1_z_yyy_xxxxyz_0,  \
                             ta1_z_yyy_xxxxyz_1,  \
                             ta1_z_yyy_xxxxz_0,   \
                             ta1_z_yyy_xxxxz_1,   \
                             ta1_z_yyy_xxxxzz_0,  \
                             ta1_z_yyy_xxxxzz_1,  \
                             ta1_z_yyy_xxxyy_0,   \
                             ta1_z_yyy_xxxyy_1,   \
                             ta1_z_yyy_xxxyyy_0,  \
                             ta1_z_yyy_xxxyyy_1,  \
                             ta1_z_yyy_xxxyyz_0,  \
                             ta1_z_yyy_xxxyyz_1,  \
                             ta1_z_yyy_xxxyz_0,   \
                             ta1_z_yyy_xxxyz_1,   \
                             ta1_z_yyy_xxxyzz_0,  \
                             ta1_z_yyy_xxxyzz_1,  \
                             ta1_z_yyy_xxxzz_0,   \
                             ta1_z_yyy_xxxzz_1,   \
                             ta1_z_yyy_xxxzzz_0,  \
                             ta1_z_yyy_xxxzzz_1,  \
                             ta1_z_yyy_xxyyy_0,   \
                             ta1_z_yyy_xxyyy_1,   \
                             ta1_z_yyy_xxyyyy_0,  \
                             ta1_z_yyy_xxyyyy_1,  \
                             ta1_z_yyy_xxyyyz_0,  \
                             ta1_z_yyy_xxyyyz_1,  \
                             ta1_z_yyy_xxyyz_0,   \
                             ta1_z_yyy_xxyyz_1,   \
                             ta1_z_yyy_xxyyzz_0,  \
                             ta1_z_yyy_xxyyzz_1,  \
                             ta1_z_yyy_xxyzz_0,   \
                             ta1_z_yyy_xxyzz_1,   \
                             ta1_z_yyy_xxyzzz_0,  \
                             ta1_z_yyy_xxyzzz_1,  \
                             ta1_z_yyy_xxzzz_0,   \
                             ta1_z_yyy_xxzzz_1,   \
                             ta1_z_yyy_xxzzzz_0,  \
                             ta1_z_yyy_xxzzzz_1,  \
                             ta1_z_yyy_xyyyy_0,   \
                             ta1_z_yyy_xyyyy_1,   \
                             ta1_z_yyy_xyyyyy_0,  \
                             ta1_z_yyy_xyyyyy_1,  \
                             ta1_z_yyy_xyyyyz_0,  \
                             ta1_z_yyy_xyyyyz_1,  \
                             ta1_z_yyy_xyyyz_0,   \
                             ta1_z_yyy_xyyyz_1,   \
                             ta1_z_yyy_xyyyzz_0,  \
                             ta1_z_yyy_xyyyzz_1,  \
                             ta1_z_yyy_xyyzz_0,   \
                             ta1_z_yyy_xyyzz_1,   \
                             ta1_z_yyy_xyyzzz_0,  \
                             ta1_z_yyy_xyyzzz_1,  \
                             ta1_z_yyy_xyzzz_0,   \
                             ta1_z_yyy_xyzzz_1,   \
                             ta1_z_yyy_xyzzzz_0,  \
                             ta1_z_yyy_xyzzzz_1,  \
                             ta1_z_yyy_xzzzz_0,   \
                             ta1_z_yyy_xzzzz_1,   \
                             ta1_z_yyy_xzzzzz_0,  \
                             ta1_z_yyy_xzzzzz_1,  \
                             ta1_z_yyy_yyyyy_0,   \
                             ta1_z_yyy_yyyyy_1,   \
                             ta1_z_yyy_yyyyyy_0,  \
                             ta1_z_yyy_yyyyyy_1,  \
                             ta1_z_yyy_yyyyyz_0,  \
                             ta1_z_yyy_yyyyyz_1,  \
                             ta1_z_yyy_yyyyz_0,   \
                             ta1_z_yyy_yyyyz_1,   \
                             ta1_z_yyy_yyyyzz_0,  \
                             ta1_z_yyy_yyyyzz_1,  \
                             ta1_z_yyy_yyyzz_0,   \
                             ta1_z_yyy_yyyzz_1,   \
                             ta1_z_yyy_yyyzzz_0,  \
                             ta1_z_yyy_yyyzzz_1,  \
                             ta1_z_yyy_yyzzz_0,   \
                             ta1_z_yyy_yyzzz_1,   \
                             ta1_z_yyy_yyzzzz_0,  \
                             ta1_z_yyy_yyzzzz_1,  \
                             ta1_z_yyy_yzzzz_0,   \
                             ta1_z_yyy_yzzzz_1,   \
                             ta1_z_yyy_yzzzzz_0,  \
                             ta1_z_yyy_yzzzzz_1,  \
                             ta1_z_yyy_zzzzz_0,   \
                             ta1_z_yyy_zzzzz_1,   \
                             ta1_z_yyy_zzzzzz_0,  \
                             ta1_z_yyy_zzzzzz_1,  \
                             ta1_z_yyyy_xxxxxx_0, \
                             ta1_z_yyyy_xxxxxy_0, \
                             ta1_z_yyyy_xxxxxz_0, \
                             ta1_z_yyyy_xxxxyy_0, \
                             ta1_z_yyyy_xxxxyz_0, \
                             ta1_z_yyyy_xxxxzz_0, \
                             ta1_z_yyyy_xxxyyy_0, \
                             ta1_z_yyyy_xxxyyz_0, \
                             ta1_z_yyyy_xxxyzz_0, \
                             ta1_z_yyyy_xxxzzz_0, \
                             ta1_z_yyyy_xxyyyy_0, \
                             ta1_z_yyyy_xxyyyz_0, \
                             ta1_z_yyyy_xxyyzz_0, \
                             ta1_z_yyyy_xxyzzz_0, \
                             ta1_z_yyyy_xxzzzz_0, \
                             ta1_z_yyyy_xyyyyy_0, \
                             ta1_z_yyyy_xyyyyz_0, \
                             ta1_z_yyyy_xyyyzz_0, \
                             ta1_z_yyyy_xyyzzz_0, \
                             ta1_z_yyyy_xyzzzz_0, \
                             ta1_z_yyyy_xzzzzz_0, \
                             ta1_z_yyyy_yyyyyy_0, \
                             ta1_z_yyyy_yyyyyz_0, \
                             ta1_z_yyyy_yyyyzz_0, \
                             ta1_z_yyyy_yyyzzz_0, \
                             ta1_z_yyyy_yyzzzz_0, \
                             ta1_z_yyyy_yzzzzz_0, \
                             ta1_z_yyyy_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyy_xxxxxx_0[i] =
            3.0 * ta1_z_yy_xxxxxx_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxxx_1[i] * fe_0 + ta1_z_yyy_xxxxxx_0[i] * pa_y[i] - ta1_z_yyy_xxxxxx_1[i] * pc_y[i];

        ta1_z_yyyy_xxxxxy_0[i] = 3.0 * ta1_z_yy_xxxxxy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxxy_1[i] * fe_0 + ta1_z_yyy_xxxxx_0[i] * fe_0 -
                                 ta1_z_yyy_xxxxx_1[i] * fe_0 + ta1_z_yyy_xxxxxy_0[i] * pa_y[i] - ta1_z_yyy_xxxxxy_1[i] * pc_y[i];

        ta1_z_yyyy_xxxxxz_0[i] =
            3.0 * ta1_z_yy_xxxxxz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxxz_1[i] * fe_0 + ta1_z_yyy_xxxxxz_0[i] * pa_y[i] - ta1_z_yyy_xxxxxz_1[i] * pc_y[i];

        ta1_z_yyyy_xxxxyy_0[i] = 3.0 * ta1_z_yy_xxxxyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxyy_1[i] * fe_0 + 2.0 * ta1_z_yyy_xxxxy_0[i] * fe_0 -
                                 2.0 * ta1_z_yyy_xxxxy_1[i] * fe_0 + ta1_z_yyy_xxxxyy_0[i] * pa_y[i] - ta1_z_yyy_xxxxyy_1[i] * pc_y[i];

        ta1_z_yyyy_xxxxyz_0[i] = 3.0 * ta1_z_yy_xxxxyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxyz_1[i] * fe_0 + ta1_z_yyy_xxxxz_0[i] * fe_0 -
                                 ta1_z_yyy_xxxxz_1[i] * fe_0 + ta1_z_yyy_xxxxyz_0[i] * pa_y[i] - ta1_z_yyy_xxxxyz_1[i] * pc_y[i];

        ta1_z_yyyy_xxxxzz_0[i] =
            3.0 * ta1_z_yy_xxxxzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxzz_1[i] * fe_0 + ta1_z_yyy_xxxxzz_0[i] * pa_y[i] - ta1_z_yyy_xxxxzz_1[i] * pc_y[i];

        ta1_z_yyyy_xxxyyy_0[i] = 3.0 * ta1_z_yy_xxxyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxyyy_1[i] * fe_0 + 3.0 * ta1_z_yyy_xxxyy_0[i] * fe_0 -
                                 3.0 * ta1_z_yyy_xxxyy_1[i] * fe_0 + ta1_z_yyy_xxxyyy_0[i] * pa_y[i] - ta1_z_yyy_xxxyyy_1[i] * pc_y[i];

        ta1_z_yyyy_xxxyyz_0[i] = 3.0 * ta1_z_yy_xxxyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxyyz_1[i] * fe_0 + 2.0 * ta1_z_yyy_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_z_yyy_xxxyz_1[i] * fe_0 + ta1_z_yyy_xxxyyz_0[i] * pa_y[i] - ta1_z_yyy_xxxyyz_1[i] * pc_y[i];

        ta1_z_yyyy_xxxyzz_0[i] = 3.0 * ta1_z_yy_xxxyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxyzz_1[i] * fe_0 + ta1_z_yyy_xxxzz_0[i] * fe_0 -
                                 ta1_z_yyy_xxxzz_1[i] * fe_0 + ta1_z_yyy_xxxyzz_0[i] * pa_y[i] - ta1_z_yyy_xxxyzz_1[i] * pc_y[i];

        ta1_z_yyyy_xxxzzz_0[i] =
            3.0 * ta1_z_yy_xxxzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxzzz_1[i] * fe_0 + ta1_z_yyy_xxxzzz_0[i] * pa_y[i] - ta1_z_yyy_xxxzzz_1[i] * pc_y[i];

        ta1_z_yyyy_xxyyyy_0[i] = 3.0 * ta1_z_yy_xxyyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyyyy_1[i] * fe_0 + 4.0 * ta1_z_yyy_xxyyy_0[i] * fe_0 -
                                 4.0 * ta1_z_yyy_xxyyy_1[i] * fe_0 + ta1_z_yyy_xxyyyy_0[i] * pa_y[i] - ta1_z_yyy_xxyyyy_1[i] * pc_y[i];

        ta1_z_yyyy_xxyyyz_0[i] = 3.0 * ta1_z_yy_xxyyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyyyz_1[i] * fe_0 + 3.0 * ta1_z_yyy_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_z_yyy_xxyyz_1[i] * fe_0 + ta1_z_yyy_xxyyyz_0[i] * pa_y[i] - ta1_z_yyy_xxyyyz_1[i] * pc_y[i];

        ta1_z_yyyy_xxyyzz_0[i] = 3.0 * ta1_z_yy_xxyyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_yyy_xxyzz_0[i] * fe_0 -
                                 2.0 * ta1_z_yyy_xxyzz_1[i] * fe_0 + ta1_z_yyy_xxyyzz_0[i] * pa_y[i] - ta1_z_yyy_xxyyzz_1[i] * pc_y[i];

        ta1_z_yyyy_xxyzzz_0[i] = 3.0 * ta1_z_yy_xxyzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyzzz_1[i] * fe_0 + ta1_z_yyy_xxzzz_0[i] * fe_0 -
                                 ta1_z_yyy_xxzzz_1[i] * fe_0 + ta1_z_yyy_xxyzzz_0[i] * pa_y[i] - ta1_z_yyy_xxyzzz_1[i] * pc_y[i];

        ta1_z_yyyy_xxzzzz_0[i] =
            3.0 * ta1_z_yy_xxzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxzzzz_1[i] * fe_0 + ta1_z_yyy_xxzzzz_0[i] * pa_y[i] - ta1_z_yyy_xxzzzz_1[i] * pc_y[i];

        ta1_z_yyyy_xyyyyy_0[i] = 3.0 * ta1_z_yy_xyyyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyyyy_1[i] * fe_0 + 5.0 * ta1_z_yyy_xyyyy_0[i] * fe_0 -
                                 5.0 * ta1_z_yyy_xyyyy_1[i] * fe_0 + ta1_z_yyy_xyyyyy_0[i] * pa_y[i] - ta1_z_yyy_xyyyyy_1[i] * pc_y[i];

        ta1_z_yyyy_xyyyyz_0[i] = 3.0 * ta1_z_yy_xyyyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyyyz_1[i] * fe_0 + 4.0 * ta1_z_yyy_xyyyz_0[i] * fe_0 -
                                 4.0 * ta1_z_yyy_xyyyz_1[i] * fe_0 + ta1_z_yyy_xyyyyz_0[i] * pa_y[i] - ta1_z_yyy_xyyyyz_1[i] * pc_y[i];

        ta1_z_yyyy_xyyyzz_0[i] = 3.0 * ta1_z_yy_xyyyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyyzz_1[i] * fe_0 + 3.0 * ta1_z_yyy_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_yyy_xyyzz_1[i] * fe_0 + ta1_z_yyy_xyyyzz_0[i] * pa_y[i] - ta1_z_yyy_xyyyzz_1[i] * pc_y[i];

        ta1_z_yyyy_xyyzzz_0[i] = 3.0 * ta1_z_yy_xyyzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyzzz_1[i] * fe_0 + 2.0 * ta1_z_yyy_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_yyy_xyzzz_1[i] * fe_0 + ta1_z_yyy_xyyzzz_0[i] * pa_y[i] - ta1_z_yyy_xyyzzz_1[i] * pc_y[i];

        ta1_z_yyyy_xyzzzz_0[i] = 3.0 * ta1_z_yy_xyzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyzzzz_1[i] * fe_0 + ta1_z_yyy_xzzzz_0[i] * fe_0 -
                                 ta1_z_yyy_xzzzz_1[i] * fe_0 + ta1_z_yyy_xyzzzz_0[i] * pa_y[i] - ta1_z_yyy_xyzzzz_1[i] * pc_y[i];

        ta1_z_yyyy_xzzzzz_0[i] =
            3.0 * ta1_z_yy_xzzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xzzzzz_1[i] * fe_0 + ta1_z_yyy_xzzzzz_0[i] * pa_y[i] - ta1_z_yyy_xzzzzz_1[i] * pc_y[i];

        ta1_z_yyyy_yyyyyy_0[i] = 3.0 * ta1_z_yy_yyyyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyyyy_1[i] * fe_0 + 6.0 * ta1_z_yyy_yyyyy_0[i] * fe_0 -
                                 6.0 * ta1_z_yyy_yyyyy_1[i] * fe_0 + ta1_z_yyy_yyyyyy_0[i] * pa_y[i] - ta1_z_yyy_yyyyyy_1[i] * pc_y[i];

        ta1_z_yyyy_yyyyyz_0[i] = 3.0 * ta1_z_yy_yyyyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyyyz_1[i] * fe_0 + 5.0 * ta1_z_yyy_yyyyz_0[i] * fe_0 -
                                 5.0 * ta1_z_yyy_yyyyz_1[i] * fe_0 + ta1_z_yyy_yyyyyz_0[i] * pa_y[i] - ta1_z_yyy_yyyyyz_1[i] * pc_y[i];

        ta1_z_yyyy_yyyyzz_0[i] = 3.0 * ta1_z_yy_yyyyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyyzz_1[i] * fe_0 + 4.0 * ta1_z_yyy_yyyzz_0[i] * fe_0 -
                                 4.0 * ta1_z_yyy_yyyzz_1[i] * fe_0 + ta1_z_yyy_yyyyzz_0[i] * pa_y[i] - ta1_z_yyy_yyyyzz_1[i] * pc_y[i];

        ta1_z_yyyy_yyyzzz_0[i] = 3.0 * ta1_z_yy_yyyzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyzzz_1[i] * fe_0 + 3.0 * ta1_z_yyy_yyzzz_0[i] * fe_0 -
                                 3.0 * ta1_z_yyy_yyzzz_1[i] * fe_0 + ta1_z_yyy_yyyzzz_0[i] * pa_y[i] - ta1_z_yyy_yyyzzz_1[i] * pc_y[i];

        ta1_z_yyyy_yyzzzz_0[i] = 3.0 * ta1_z_yy_yyzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyzzzz_1[i] * fe_0 + 2.0 * ta1_z_yyy_yzzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_yyy_yzzzz_1[i] * fe_0 + ta1_z_yyy_yyzzzz_0[i] * pa_y[i] - ta1_z_yyy_yyzzzz_1[i] * pc_y[i];

        ta1_z_yyyy_yzzzzz_0[i] = 3.0 * ta1_z_yy_yzzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yzzzzz_1[i] * fe_0 + ta1_z_yyy_zzzzz_0[i] * fe_0 -
                                 ta1_z_yyy_zzzzz_1[i] * fe_0 + ta1_z_yyy_yzzzzz_0[i] * pa_y[i] - ta1_z_yyy_yzzzzz_1[i] * pc_y[i];

        ta1_z_yyyy_zzzzzz_0[i] =
            3.0 * ta1_z_yy_zzzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_zzzzzz_1[i] * fe_0 + ta1_z_yyy_zzzzzz_0[i] * pa_y[i] - ta1_z_yyy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 1148-1176 components of targeted buffer : GI

    auto ta1_z_yyyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1148);

    auto ta1_z_yyyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1149);

    auto ta1_z_yyyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1150);

    auto ta1_z_yyyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1151);

    auto ta1_z_yyyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1152);

    auto ta1_z_yyyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1153);

    auto ta1_z_yyyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1154);

    auto ta1_z_yyyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1155);

    auto ta1_z_yyyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1156);

    auto ta1_z_yyyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1157);

    auto ta1_z_yyyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1158);

    auto ta1_z_yyyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1159);

    auto ta1_z_yyyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1160);

    auto ta1_z_yyyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1161);

    auto ta1_z_yyyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1162);

    auto ta1_z_yyyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1163);

    auto ta1_z_yyyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1164);

    auto ta1_z_yyyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1165);

    auto ta1_z_yyyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1166);

    auto ta1_z_yyyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1167);

    auto ta1_z_yyyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1168);

    auto ta1_z_yyyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1169);

    auto ta1_z_yyyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1170);

    auto ta1_z_yyyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1171);

    auto ta1_z_yyyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1172);

    auto ta1_z_yyyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1173);

    auto ta1_z_yyyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1174);

    auto ta1_z_yyyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1175);

#pragma omp simd aligned(pa_y,                    \
                             pa_z,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_z_yyy_xxxxxx_0,  \
                             ta1_z_yyy_xxxxxx_1,  \
                             ta1_z_yyy_xxxxxy_0,  \
                             ta1_z_yyy_xxxxxy_1,  \
                             ta1_z_yyy_xxxxy_0,   \
                             ta1_z_yyy_xxxxy_1,   \
                             ta1_z_yyy_xxxxyy_0,  \
                             ta1_z_yyy_xxxxyy_1,  \
                             ta1_z_yyy_xxxxyz_0,  \
                             ta1_z_yyy_xxxxyz_1,  \
                             ta1_z_yyy_xxxyy_0,   \
                             ta1_z_yyy_xxxyy_1,   \
                             ta1_z_yyy_xxxyyy_0,  \
                             ta1_z_yyy_xxxyyy_1,  \
                             ta1_z_yyy_xxxyyz_0,  \
                             ta1_z_yyy_xxxyyz_1,  \
                             ta1_z_yyy_xxxyz_0,   \
                             ta1_z_yyy_xxxyz_1,   \
                             ta1_z_yyy_xxxyzz_0,  \
                             ta1_z_yyy_xxxyzz_1,  \
                             ta1_z_yyy_xxyyy_0,   \
                             ta1_z_yyy_xxyyy_1,   \
                             ta1_z_yyy_xxyyyy_0,  \
                             ta1_z_yyy_xxyyyy_1,  \
                             ta1_z_yyy_xxyyyz_0,  \
                             ta1_z_yyy_xxyyyz_1,  \
                             ta1_z_yyy_xxyyz_0,   \
                             ta1_z_yyy_xxyyz_1,   \
                             ta1_z_yyy_xxyyzz_0,  \
                             ta1_z_yyy_xxyyzz_1,  \
                             ta1_z_yyy_xxyzz_0,   \
                             ta1_z_yyy_xxyzz_1,   \
                             ta1_z_yyy_xxyzzz_0,  \
                             ta1_z_yyy_xxyzzz_1,  \
                             ta1_z_yyy_xyyyy_0,   \
                             ta1_z_yyy_xyyyy_1,   \
                             ta1_z_yyy_xyyyyy_0,  \
                             ta1_z_yyy_xyyyyy_1,  \
                             ta1_z_yyy_xyyyyz_0,  \
                             ta1_z_yyy_xyyyyz_1,  \
                             ta1_z_yyy_xyyyz_0,   \
                             ta1_z_yyy_xyyyz_1,   \
                             ta1_z_yyy_xyyyzz_0,  \
                             ta1_z_yyy_xyyyzz_1,  \
                             ta1_z_yyy_xyyzz_0,   \
                             ta1_z_yyy_xyyzz_1,   \
                             ta1_z_yyy_xyyzzz_0,  \
                             ta1_z_yyy_xyyzzz_1,  \
                             ta1_z_yyy_xyzzz_0,   \
                             ta1_z_yyy_xyzzz_1,   \
                             ta1_z_yyy_xyzzzz_0,  \
                             ta1_z_yyy_xyzzzz_1,  \
                             ta1_z_yyy_yyyyy_0,   \
                             ta1_z_yyy_yyyyy_1,   \
                             ta1_z_yyy_yyyyyy_0,  \
                             ta1_z_yyy_yyyyyy_1,  \
                             ta1_z_yyy_yyyyyz_0,  \
                             ta1_z_yyy_yyyyyz_1,  \
                             ta1_z_yyy_yyyyz_0,   \
                             ta1_z_yyy_yyyyz_1,   \
                             ta1_z_yyy_yyyyzz_0,  \
                             ta1_z_yyy_yyyyzz_1,  \
                             ta1_z_yyy_yyyzz_0,   \
                             ta1_z_yyy_yyyzz_1,   \
                             ta1_z_yyy_yyyzzz_0,  \
                             ta1_z_yyy_yyyzzz_1,  \
                             ta1_z_yyy_yyzzz_0,   \
                             ta1_z_yyy_yyzzz_1,   \
                             ta1_z_yyy_yyzzzz_0,  \
                             ta1_z_yyy_yyzzzz_1,  \
                             ta1_z_yyy_yzzzz_0,   \
                             ta1_z_yyy_yzzzz_1,   \
                             ta1_z_yyy_yzzzzz_0,  \
                             ta1_z_yyy_yzzzzz_1,  \
                             ta1_z_yyyz_xxxxxx_0, \
                             ta1_z_yyyz_xxxxxy_0, \
                             ta1_z_yyyz_xxxxxz_0, \
                             ta1_z_yyyz_xxxxyy_0, \
                             ta1_z_yyyz_xxxxyz_0, \
                             ta1_z_yyyz_xxxxzz_0, \
                             ta1_z_yyyz_xxxyyy_0, \
                             ta1_z_yyyz_xxxyyz_0, \
                             ta1_z_yyyz_xxxyzz_0, \
                             ta1_z_yyyz_xxxzzz_0, \
                             ta1_z_yyyz_xxyyyy_0, \
                             ta1_z_yyyz_xxyyyz_0, \
                             ta1_z_yyyz_xxyyzz_0, \
                             ta1_z_yyyz_xxyzzz_0, \
                             ta1_z_yyyz_xxzzzz_0, \
                             ta1_z_yyyz_xyyyyy_0, \
                             ta1_z_yyyz_xyyyyz_0, \
                             ta1_z_yyyz_xyyyzz_0, \
                             ta1_z_yyyz_xyyzzz_0, \
                             ta1_z_yyyz_xyzzzz_0, \
                             ta1_z_yyyz_xzzzzz_0, \
                             ta1_z_yyyz_yyyyyy_0, \
                             ta1_z_yyyz_yyyyyz_0, \
                             ta1_z_yyyz_yyyyzz_0, \
                             ta1_z_yyyz_yyyzzz_0, \
                             ta1_z_yyyz_yyzzzz_0, \
                             ta1_z_yyyz_yzzzzz_0, \
                             ta1_z_yyyz_zzzzzz_0, \
                             ta1_z_yyz_xxxxxz_0,  \
                             ta1_z_yyz_xxxxxz_1,  \
                             ta1_z_yyz_xxxxzz_0,  \
                             ta1_z_yyz_xxxxzz_1,  \
                             ta1_z_yyz_xxxzzz_0,  \
                             ta1_z_yyz_xxxzzz_1,  \
                             ta1_z_yyz_xxzzzz_0,  \
                             ta1_z_yyz_xxzzzz_1,  \
                             ta1_z_yyz_xzzzzz_0,  \
                             ta1_z_yyz_xzzzzz_1,  \
                             ta1_z_yyz_zzzzzz_0,  \
                             ta1_z_yyz_zzzzzz_1,  \
                             ta1_z_yz_xxxxxz_0,   \
                             ta1_z_yz_xxxxxz_1,   \
                             ta1_z_yz_xxxxzz_0,   \
                             ta1_z_yz_xxxxzz_1,   \
                             ta1_z_yz_xxxzzz_0,   \
                             ta1_z_yz_xxxzzz_1,   \
                             ta1_z_yz_xxzzzz_0,   \
                             ta1_z_yz_xxzzzz_1,   \
                             ta1_z_yz_xzzzzz_0,   \
                             ta1_z_yz_xzzzzz_1,   \
                             ta1_z_yz_zzzzzz_0,   \
                             ta1_z_yz_zzzzzz_1,   \
                             ta_yyy_xxxxxx_1,     \
                             ta_yyy_xxxxxy_1,     \
                             ta_yyy_xxxxyy_1,     \
                             ta_yyy_xxxxyz_1,     \
                             ta_yyy_xxxyyy_1,     \
                             ta_yyy_xxxyyz_1,     \
                             ta_yyy_xxxyzz_1,     \
                             ta_yyy_xxyyyy_1,     \
                             ta_yyy_xxyyyz_1,     \
                             ta_yyy_xxyyzz_1,     \
                             ta_yyy_xxyzzz_1,     \
                             ta_yyy_xyyyyy_1,     \
                             ta_yyy_xyyyyz_1,     \
                             ta_yyy_xyyyzz_1,     \
                             ta_yyy_xyyzzz_1,     \
                             ta_yyy_xyzzzz_1,     \
                             ta_yyy_yyyyyy_1,     \
                             ta_yyy_yyyyyz_1,     \
                             ta_yyy_yyyyzz_1,     \
                             ta_yyy_yyyzzz_1,     \
                             ta_yyy_yyzzzz_1,     \
                             ta_yyy_yzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyz_xxxxxx_0[i] = ta_yyy_xxxxxx_1[i] + ta1_z_yyy_xxxxxx_0[i] * pa_z[i] - ta1_z_yyy_xxxxxx_1[i] * pc_z[i];

        ta1_z_yyyz_xxxxxy_0[i] = ta_yyy_xxxxxy_1[i] + ta1_z_yyy_xxxxxy_0[i] * pa_z[i] - ta1_z_yyy_xxxxxy_1[i] * pc_z[i];

        ta1_z_yyyz_xxxxxz_0[i] =
            2.0 * ta1_z_yz_xxxxxz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxxxxz_1[i] * fe_0 + ta1_z_yyz_xxxxxz_0[i] * pa_y[i] - ta1_z_yyz_xxxxxz_1[i] * pc_y[i];

        ta1_z_yyyz_xxxxyy_0[i] = ta_yyy_xxxxyy_1[i] + ta1_z_yyy_xxxxyy_0[i] * pa_z[i] - ta1_z_yyy_xxxxyy_1[i] * pc_z[i];

        ta1_z_yyyz_xxxxyz_0[i] = ta1_z_yyy_xxxxy_0[i] * fe_0 - ta1_z_yyy_xxxxy_1[i] * fe_0 + ta_yyy_xxxxyz_1[i] + ta1_z_yyy_xxxxyz_0[i] * pa_z[i] -
                                 ta1_z_yyy_xxxxyz_1[i] * pc_z[i];

        ta1_z_yyyz_xxxxzz_0[i] =
            2.0 * ta1_z_yz_xxxxzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxxxzz_1[i] * fe_0 + ta1_z_yyz_xxxxzz_0[i] * pa_y[i] - ta1_z_yyz_xxxxzz_1[i] * pc_y[i];

        ta1_z_yyyz_xxxyyy_0[i] = ta_yyy_xxxyyy_1[i] + ta1_z_yyy_xxxyyy_0[i] * pa_z[i] - ta1_z_yyy_xxxyyy_1[i] * pc_z[i];

        ta1_z_yyyz_xxxyyz_0[i] = ta1_z_yyy_xxxyy_0[i] * fe_0 - ta1_z_yyy_xxxyy_1[i] * fe_0 + ta_yyy_xxxyyz_1[i] + ta1_z_yyy_xxxyyz_0[i] * pa_z[i] -
                                 ta1_z_yyy_xxxyyz_1[i] * pc_z[i];

        ta1_z_yyyz_xxxyzz_0[i] = 2.0 * ta1_z_yyy_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xxxyz_1[i] * fe_0 + ta_yyy_xxxyzz_1[i] +
                                 ta1_z_yyy_xxxyzz_0[i] * pa_z[i] - ta1_z_yyy_xxxyzz_1[i] * pc_z[i];

        ta1_z_yyyz_xxxzzz_0[i] =
            2.0 * ta1_z_yz_xxxzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxxzzz_1[i] * fe_0 + ta1_z_yyz_xxxzzz_0[i] * pa_y[i] - ta1_z_yyz_xxxzzz_1[i] * pc_y[i];

        ta1_z_yyyz_xxyyyy_0[i] = ta_yyy_xxyyyy_1[i] + ta1_z_yyy_xxyyyy_0[i] * pa_z[i] - ta1_z_yyy_xxyyyy_1[i] * pc_z[i];

        ta1_z_yyyz_xxyyyz_0[i] = ta1_z_yyy_xxyyy_0[i] * fe_0 - ta1_z_yyy_xxyyy_1[i] * fe_0 + ta_yyy_xxyyyz_1[i] + ta1_z_yyy_xxyyyz_0[i] * pa_z[i] -
                                 ta1_z_yyy_xxyyyz_1[i] * pc_z[i];

        ta1_z_yyyz_xxyyzz_0[i] = 2.0 * ta1_z_yyy_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xxyyz_1[i] * fe_0 + ta_yyy_xxyyzz_1[i] +
                                 ta1_z_yyy_xxyyzz_0[i] * pa_z[i] - ta1_z_yyy_xxyyzz_1[i] * pc_z[i];

        ta1_z_yyyz_xxyzzz_0[i] = 3.0 * ta1_z_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxyzz_1[i] * fe_0 + ta_yyy_xxyzzz_1[i] +
                                 ta1_z_yyy_xxyzzz_0[i] * pa_z[i] - ta1_z_yyy_xxyzzz_1[i] * pc_z[i];

        ta1_z_yyyz_xxzzzz_0[i] =
            2.0 * ta1_z_yz_xxzzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxzzzz_1[i] * fe_0 + ta1_z_yyz_xxzzzz_0[i] * pa_y[i] - ta1_z_yyz_xxzzzz_1[i] * pc_y[i];

        ta1_z_yyyz_xyyyyy_0[i] = ta_yyy_xyyyyy_1[i] + ta1_z_yyy_xyyyyy_0[i] * pa_z[i] - ta1_z_yyy_xyyyyy_1[i] * pc_z[i];

        ta1_z_yyyz_xyyyyz_0[i] = ta1_z_yyy_xyyyy_0[i] * fe_0 - ta1_z_yyy_xyyyy_1[i] * fe_0 + ta_yyy_xyyyyz_1[i] + ta1_z_yyy_xyyyyz_0[i] * pa_z[i] -
                                 ta1_z_yyy_xyyyyz_1[i] * pc_z[i];

        ta1_z_yyyz_xyyyzz_0[i] = 2.0 * ta1_z_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyyyz_1[i] * fe_0 + ta_yyy_xyyyzz_1[i] +
                                 ta1_z_yyy_xyyyzz_0[i] * pa_z[i] - ta1_z_yyy_xyyyzz_1[i] * pc_z[i];

        ta1_z_yyyz_xyyzzz_0[i] = 3.0 * ta1_z_yyy_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xyyzz_1[i] * fe_0 + ta_yyy_xyyzzz_1[i] +
                                 ta1_z_yyy_xyyzzz_0[i] * pa_z[i] - ta1_z_yyy_xyyzzz_1[i] * pc_z[i];

        ta1_z_yyyz_xyzzzz_0[i] = 4.0 * ta1_z_yyy_xyzzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xyzzz_1[i] * fe_0 + ta_yyy_xyzzzz_1[i] +
                                 ta1_z_yyy_xyzzzz_0[i] * pa_z[i] - ta1_z_yyy_xyzzzz_1[i] * pc_z[i];

        ta1_z_yyyz_xzzzzz_0[i] =
            2.0 * ta1_z_yz_xzzzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xzzzzz_1[i] * fe_0 + ta1_z_yyz_xzzzzz_0[i] * pa_y[i] - ta1_z_yyz_xzzzzz_1[i] * pc_y[i];

        ta1_z_yyyz_yyyyyy_0[i] = ta_yyy_yyyyyy_1[i] + ta1_z_yyy_yyyyyy_0[i] * pa_z[i] - ta1_z_yyy_yyyyyy_1[i] * pc_z[i];

        ta1_z_yyyz_yyyyyz_0[i] = ta1_z_yyy_yyyyy_0[i] * fe_0 - ta1_z_yyy_yyyyy_1[i] * fe_0 + ta_yyy_yyyyyz_1[i] + ta1_z_yyy_yyyyyz_0[i] * pa_z[i] -
                                 ta1_z_yyy_yyyyyz_1[i] * pc_z[i];

        ta1_z_yyyz_yyyyzz_0[i] = 2.0 * ta1_z_yyy_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_yyyyz_1[i] * fe_0 + ta_yyy_yyyyzz_1[i] +
                                 ta1_z_yyy_yyyyzz_0[i] * pa_z[i] - ta1_z_yyy_yyyyzz_1[i] * pc_z[i];

        ta1_z_yyyz_yyyzzz_0[i] = 3.0 * ta1_z_yyy_yyyzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_yyyzz_1[i] * fe_0 + ta_yyy_yyyzzz_1[i] +
                                 ta1_z_yyy_yyyzzz_0[i] * pa_z[i] - ta1_z_yyy_yyyzzz_1[i] * pc_z[i];

        ta1_z_yyyz_yyzzzz_0[i] = 4.0 * ta1_z_yyy_yyzzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yyzzz_1[i] * fe_0 + ta_yyy_yyzzzz_1[i] +
                                 ta1_z_yyy_yyzzzz_0[i] * pa_z[i] - ta1_z_yyy_yyzzzz_1[i] * pc_z[i];

        ta1_z_yyyz_yzzzzz_0[i] = 5.0 * ta1_z_yyy_yzzzz_0[i] * fe_0 - 5.0 * ta1_z_yyy_yzzzz_1[i] * fe_0 + ta_yyy_yzzzzz_1[i] +
                                 ta1_z_yyy_yzzzzz_0[i] * pa_z[i] - ta1_z_yyy_yzzzzz_1[i] * pc_z[i];

        ta1_z_yyyz_zzzzzz_0[i] =
            2.0 * ta1_z_yz_zzzzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_zzzzzz_1[i] * fe_0 + ta1_z_yyz_zzzzzz_0[i] * pa_y[i] - ta1_z_yyz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 1176-1204 components of targeted buffer : GI

    auto ta1_z_yyzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1176);

    auto ta1_z_yyzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1177);

    auto ta1_z_yyzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1178);

    auto ta1_z_yyzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1179);

    auto ta1_z_yyzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1180);

    auto ta1_z_yyzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1181);

    auto ta1_z_yyzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1182);

    auto ta1_z_yyzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1183);

    auto ta1_z_yyzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1184);

    auto ta1_z_yyzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1185);

    auto ta1_z_yyzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1186);

    auto ta1_z_yyzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1187);

    auto ta1_z_yyzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1188);

    auto ta1_z_yyzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1189);

    auto ta1_z_yyzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1190);

    auto ta1_z_yyzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1191);

    auto ta1_z_yyzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1192);

    auto ta1_z_yyzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1193);

    auto ta1_z_yyzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1194);

    auto ta1_z_yyzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1195);

    auto ta1_z_yyzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1196);

    auto ta1_z_yyzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1197);

    auto ta1_z_yyzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1198);

    auto ta1_z_yyzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1199);

    auto ta1_z_yyzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1200);

    auto ta1_z_yyzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1201);

    auto ta1_z_yyzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1202);

    auto ta1_z_yyzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1203);

#pragma omp simd aligned(pa_y,                    \
                             pa_z,                \
                             pc_y,                \
                             pc_z,                \
                             ta1_z_yy_xxxxxy_0,   \
                             ta1_z_yy_xxxxxy_1,   \
                             ta1_z_yy_xxxxyy_0,   \
                             ta1_z_yy_xxxxyy_1,   \
                             ta1_z_yy_xxxyyy_0,   \
                             ta1_z_yy_xxxyyy_1,   \
                             ta1_z_yy_xxyyyy_0,   \
                             ta1_z_yy_xxyyyy_1,   \
                             ta1_z_yy_xyyyyy_0,   \
                             ta1_z_yy_xyyyyy_1,   \
                             ta1_z_yy_yyyyyy_0,   \
                             ta1_z_yy_yyyyyy_1,   \
                             ta1_z_yyz_xxxxxy_0,  \
                             ta1_z_yyz_xxxxxy_1,  \
                             ta1_z_yyz_xxxxyy_0,  \
                             ta1_z_yyz_xxxxyy_1,  \
                             ta1_z_yyz_xxxyyy_0,  \
                             ta1_z_yyz_xxxyyy_1,  \
                             ta1_z_yyz_xxyyyy_0,  \
                             ta1_z_yyz_xxyyyy_1,  \
                             ta1_z_yyz_xyyyyy_0,  \
                             ta1_z_yyz_xyyyyy_1,  \
                             ta1_z_yyz_yyyyyy_0,  \
                             ta1_z_yyz_yyyyyy_1,  \
                             ta1_z_yyzz_xxxxxx_0, \
                             ta1_z_yyzz_xxxxxy_0, \
                             ta1_z_yyzz_xxxxxz_0, \
                             ta1_z_yyzz_xxxxyy_0, \
                             ta1_z_yyzz_xxxxyz_0, \
                             ta1_z_yyzz_xxxxzz_0, \
                             ta1_z_yyzz_xxxyyy_0, \
                             ta1_z_yyzz_xxxyyz_0, \
                             ta1_z_yyzz_xxxyzz_0, \
                             ta1_z_yyzz_xxxzzz_0, \
                             ta1_z_yyzz_xxyyyy_0, \
                             ta1_z_yyzz_xxyyyz_0, \
                             ta1_z_yyzz_xxyyzz_0, \
                             ta1_z_yyzz_xxyzzz_0, \
                             ta1_z_yyzz_xxzzzz_0, \
                             ta1_z_yyzz_xyyyyy_0, \
                             ta1_z_yyzz_xyyyyz_0, \
                             ta1_z_yyzz_xyyyzz_0, \
                             ta1_z_yyzz_xyyzzz_0, \
                             ta1_z_yyzz_xyzzzz_0, \
                             ta1_z_yyzz_xzzzzz_0, \
                             ta1_z_yyzz_yyyyyy_0, \
                             ta1_z_yyzz_yyyyyz_0, \
                             ta1_z_yyzz_yyyyzz_0, \
                             ta1_z_yyzz_yyyzzz_0, \
                             ta1_z_yyzz_yyzzzz_0, \
                             ta1_z_yyzz_yzzzzz_0, \
                             ta1_z_yyzz_zzzzzz_0, \
                             ta1_z_yzz_xxxxxx_0,  \
                             ta1_z_yzz_xxxxxx_1,  \
                             ta1_z_yzz_xxxxxz_0,  \
                             ta1_z_yzz_xxxxxz_1,  \
                             ta1_z_yzz_xxxxyz_0,  \
                             ta1_z_yzz_xxxxyz_1,  \
                             ta1_z_yzz_xxxxz_0,   \
                             ta1_z_yzz_xxxxz_1,   \
                             ta1_z_yzz_xxxxzz_0,  \
                             ta1_z_yzz_xxxxzz_1,  \
                             ta1_z_yzz_xxxyyz_0,  \
                             ta1_z_yzz_xxxyyz_1,  \
                             ta1_z_yzz_xxxyz_0,   \
                             ta1_z_yzz_xxxyz_1,   \
                             ta1_z_yzz_xxxyzz_0,  \
                             ta1_z_yzz_xxxyzz_1,  \
                             ta1_z_yzz_xxxzz_0,   \
                             ta1_z_yzz_xxxzz_1,   \
                             ta1_z_yzz_xxxzzz_0,  \
                             ta1_z_yzz_xxxzzz_1,  \
                             ta1_z_yzz_xxyyyz_0,  \
                             ta1_z_yzz_xxyyyz_1,  \
                             ta1_z_yzz_xxyyz_0,   \
                             ta1_z_yzz_xxyyz_1,   \
                             ta1_z_yzz_xxyyzz_0,  \
                             ta1_z_yzz_xxyyzz_1,  \
                             ta1_z_yzz_xxyzz_0,   \
                             ta1_z_yzz_xxyzz_1,   \
                             ta1_z_yzz_xxyzzz_0,  \
                             ta1_z_yzz_xxyzzz_1,  \
                             ta1_z_yzz_xxzzz_0,   \
                             ta1_z_yzz_xxzzz_1,   \
                             ta1_z_yzz_xxzzzz_0,  \
                             ta1_z_yzz_xxzzzz_1,  \
                             ta1_z_yzz_xyyyyz_0,  \
                             ta1_z_yzz_xyyyyz_1,  \
                             ta1_z_yzz_xyyyz_0,   \
                             ta1_z_yzz_xyyyz_1,   \
                             ta1_z_yzz_xyyyzz_0,  \
                             ta1_z_yzz_xyyyzz_1,  \
                             ta1_z_yzz_xyyzz_0,   \
                             ta1_z_yzz_xyyzz_1,   \
                             ta1_z_yzz_xyyzzz_0,  \
                             ta1_z_yzz_xyyzzz_1,  \
                             ta1_z_yzz_xyzzz_0,   \
                             ta1_z_yzz_xyzzz_1,   \
                             ta1_z_yzz_xyzzzz_0,  \
                             ta1_z_yzz_xyzzzz_1,  \
                             ta1_z_yzz_xzzzz_0,   \
                             ta1_z_yzz_xzzzz_1,   \
                             ta1_z_yzz_xzzzzz_0,  \
                             ta1_z_yzz_xzzzzz_1,  \
                             ta1_z_yzz_yyyyyz_0,  \
                             ta1_z_yzz_yyyyyz_1,  \
                             ta1_z_yzz_yyyyz_0,   \
                             ta1_z_yzz_yyyyz_1,   \
                             ta1_z_yzz_yyyyzz_0,  \
                             ta1_z_yzz_yyyyzz_1,  \
                             ta1_z_yzz_yyyzz_0,   \
                             ta1_z_yzz_yyyzz_1,   \
                             ta1_z_yzz_yyyzzz_0,  \
                             ta1_z_yzz_yyyzzz_1,  \
                             ta1_z_yzz_yyzzz_0,   \
                             ta1_z_yzz_yyzzz_1,   \
                             ta1_z_yzz_yyzzzz_0,  \
                             ta1_z_yzz_yyzzzz_1,  \
                             ta1_z_yzz_yzzzz_0,   \
                             ta1_z_yzz_yzzzz_1,   \
                             ta1_z_yzz_yzzzzz_0,  \
                             ta1_z_yzz_yzzzzz_1,  \
                             ta1_z_yzz_zzzzz_0,   \
                             ta1_z_yzz_zzzzz_1,   \
                             ta1_z_yzz_zzzzzz_0,  \
                             ta1_z_yzz_zzzzzz_1,  \
                             ta1_z_zz_xxxxxx_0,   \
                             ta1_z_zz_xxxxxx_1,   \
                             ta1_z_zz_xxxxxz_0,   \
                             ta1_z_zz_xxxxxz_1,   \
                             ta1_z_zz_xxxxyz_0,   \
                             ta1_z_zz_xxxxyz_1,   \
                             ta1_z_zz_xxxxzz_0,   \
                             ta1_z_zz_xxxxzz_1,   \
                             ta1_z_zz_xxxyyz_0,   \
                             ta1_z_zz_xxxyyz_1,   \
                             ta1_z_zz_xxxyzz_0,   \
                             ta1_z_zz_xxxyzz_1,   \
                             ta1_z_zz_xxxzzz_0,   \
                             ta1_z_zz_xxxzzz_1,   \
                             ta1_z_zz_xxyyyz_0,   \
                             ta1_z_zz_xxyyyz_1,   \
                             ta1_z_zz_xxyyzz_0,   \
                             ta1_z_zz_xxyyzz_1,   \
                             ta1_z_zz_xxyzzz_0,   \
                             ta1_z_zz_xxyzzz_1,   \
                             ta1_z_zz_xxzzzz_0,   \
                             ta1_z_zz_xxzzzz_1,   \
                             ta1_z_zz_xyyyyz_0,   \
                             ta1_z_zz_xyyyyz_1,   \
                             ta1_z_zz_xyyyzz_0,   \
                             ta1_z_zz_xyyyzz_1,   \
                             ta1_z_zz_xyyzzz_0,   \
                             ta1_z_zz_xyyzzz_1,   \
                             ta1_z_zz_xyzzzz_0,   \
                             ta1_z_zz_xyzzzz_1,   \
                             ta1_z_zz_xzzzzz_0,   \
                             ta1_z_zz_xzzzzz_1,   \
                             ta1_z_zz_yyyyyz_0,   \
                             ta1_z_zz_yyyyyz_1,   \
                             ta1_z_zz_yyyyzz_0,   \
                             ta1_z_zz_yyyyzz_1,   \
                             ta1_z_zz_yyyzzz_0,   \
                             ta1_z_zz_yyyzzz_1,   \
                             ta1_z_zz_yyzzzz_0,   \
                             ta1_z_zz_yyzzzz_1,   \
                             ta1_z_zz_yzzzzz_0,   \
                             ta1_z_zz_yzzzzz_1,   \
                             ta1_z_zz_zzzzzz_0,   \
                             ta1_z_zz_zzzzzz_1,   \
                             ta_yyz_xxxxxy_1,     \
                             ta_yyz_xxxxyy_1,     \
                             ta_yyz_xxxyyy_1,     \
                             ta_yyz_xxyyyy_1,     \
                             ta_yyz_xyyyyy_1,     \
                             ta_yyz_yyyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzz_xxxxxx_0[i] =
            ta1_z_zz_xxxxxx_0[i] * fe_0 - ta1_z_zz_xxxxxx_1[i] * fe_0 + ta1_z_yzz_xxxxxx_0[i] * pa_y[i] - ta1_z_yzz_xxxxxx_1[i] * pc_y[i];

        ta1_z_yyzz_xxxxxy_0[i] = ta1_z_yy_xxxxxy_0[i] * fe_0 - ta1_z_yy_xxxxxy_1[i] * fe_0 + ta_yyz_xxxxxy_1[i] + ta1_z_yyz_xxxxxy_0[i] * pa_z[i] -
                                 ta1_z_yyz_xxxxxy_1[i] * pc_z[i];

        ta1_z_yyzz_xxxxxz_0[i] =
            ta1_z_zz_xxxxxz_0[i] * fe_0 - ta1_z_zz_xxxxxz_1[i] * fe_0 + ta1_z_yzz_xxxxxz_0[i] * pa_y[i] - ta1_z_yzz_xxxxxz_1[i] * pc_y[i];

        ta1_z_yyzz_xxxxyy_0[i] = ta1_z_yy_xxxxyy_0[i] * fe_0 - ta1_z_yy_xxxxyy_1[i] * fe_0 + ta_yyz_xxxxyy_1[i] + ta1_z_yyz_xxxxyy_0[i] * pa_z[i] -
                                 ta1_z_yyz_xxxxyy_1[i] * pc_z[i];

        ta1_z_yyzz_xxxxyz_0[i] = ta1_z_zz_xxxxyz_0[i] * fe_0 - ta1_z_zz_xxxxyz_1[i] * fe_0 + ta1_z_yzz_xxxxz_0[i] * fe_0 -
                                 ta1_z_yzz_xxxxz_1[i] * fe_0 + ta1_z_yzz_xxxxyz_0[i] * pa_y[i] - ta1_z_yzz_xxxxyz_1[i] * pc_y[i];

        ta1_z_yyzz_xxxxzz_0[i] =
            ta1_z_zz_xxxxzz_0[i] * fe_0 - ta1_z_zz_xxxxzz_1[i] * fe_0 + ta1_z_yzz_xxxxzz_0[i] * pa_y[i] - ta1_z_yzz_xxxxzz_1[i] * pc_y[i];

        ta1_z_yyzz_xxxyyy_0[i] = ta1_z_yy_xxxyyy_0[i] * fe_0 - ta1_z_yy_xxxyyy_1[i] * fe_0 + ta_yyz_xxxyyy_1[i] + ta1_z_yyz_xxxyyy_0[i] * pa_z[i] -
                                 ta1_z_yyz_xxxyyy_1[i] * pc_z[i];

        ta1_z_yyzz_xxxyyz_0[i] = ta1_z_zz_xxxyyz_0[i] * fe_0 - ta1_z_zz_xxxyyz_1[i] * fe_0 + 2.0 * ta1_z_yzz_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_z_yzz_xxxyz_1[i] * fe_0 + ta1_z_yzz_xxxyyz_0[i] * pa_y[i] - ta1_z_yzz_xxxyyz_1[i] * pc_y[i];

        ta1_z_yyzz_xxxyzz_0[i] = ta1_z_zz_xxxyzz_0[i] * fe_0 - ta1_z_zz_xxxyzz_1[i] * fe_0 + ta1_z_yzz_xxxzz_0[i] * fe_0 -
                                 ta1_z_yzz_xxxzz_1[i] * fe_0 + ta1_z_yzz_xxxyzz_0[i] * pa_y[i] - ta1_z_yzz_xxxyzz_1[i] * pc_y[i];

        ta1_z_yyzz_xxxzzz_0[i] =
            ta1_z_zz_xxxzzz_0[i] * fe_0 - ta1_z_zz_xxxzzz_1[i] * fe_0 + ta1_z_yzz_xxxzzz_0[i] * pa_y[i] - ta1_z_yzz_xxxzzz_1[i] * pc_y[i];

        ta1_z_yyzz_xxyyyy_0[i] = ta1_z_yy_xxyyyy_0[i] * fe_0 - ta1_z_yy_xxyyyy_1[i] * fe_0 + ta_yyz_xxyyyy_1[i] + ta1_z_yyz_xxyyyy_0[i] * pa_z[i] -
                                 ta1_z_yyz_xxyyyy_1[i] * pc_z[i];

        ta1_z_yyzz_xxyyyz_0[i] = ta1_z_zz_xxyyyz_0[i] * fe_0 - ta1_z_zz_xxyyyz_1[i] * fe_0 + 3.0 * ta1_z_yzz_xxyyz_0[i] * fe_0 -
                                 3.0 * ta1_z_yzz_xxyyz_1[i] * fe_0 + ta1_z_yzz_xxyyyz_0[i] * pa_y[i] - ta1_z_yzz_xxyyyz_1[i] * pc_y[i];

        ta1_z_yyzz_xxyyzz_0[i] = ta1_z_zz_xxyyzz_0[i] * fe_0 - ta1_z_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_yzz_xxyzz_0[i] * fe_0 -
                                 2.0 * ta1_z_yzz_xxyzz_1[i] * fe_0 + ta1_z_yzz_xxyyzz_0[i] * pa_y[i] - ta1_z_yzz_xxyyzz_1[i] * pc_y[i];

        ta1_z_yyzz_xxyzzz_0[i] = ta1_z_zz_xxyzzz_0[i] * fe_0 - ta1_z_zz_xxyzzz_1[i] * fe_0 + ta1_z_yzz_xxzzz_0[i] * fe_0 -
                                 ta1_z_yzz_xxzzz_1[i] * fe_0 + ta1_z_yzz_xxyzzz_0[i] * pa_y[i] - ta1_z_yzz_xxyzzz_1[i] * pc_y[i];

        ta1_z_yyzz_xxzzzz_0[i] =
            ta1_z_zz_xxzzzz_0[i] * fe_0 - ta1_z_zz_xxzzzz_1[i] * fe_0 + ta1_z_yzz_xxzzzz_0[i] * pa_y[i] - ta1_z_yzz_xxzzzz_1[i] * pc_y[i];

        ta1_z_yyzz_xyyyyy_0[i] = ta1_z_yy_xyyyyy_0[i] * fe_0 - ta1_z_yy_xyyyyy_1[i] * fe_0 + ta_yyz_xyyyyy_1[i] + ta1_z_yyz_xyyyyy_0[i] * pa_z[i] -
                                 ta1_z_yyz_xyyyyy_1[i] * pc_z[i];

        ta1_z_yyzz_xyyyyz_0[i] = ta1_z_zz_xyyyyz_0[i] * fe_0 - ta1_z_zz_xyyyyz_1[i] * fe_0 + 4.0 * ta1_z_yzz_xyyyz_0[i] * fe_0 -
                                 4.0 * ta1_z_yzz_xyyyz_1[i] * fe_0 + ta1_z_yzz_xyyyyz_0[i] * pa_y[i] - ta1_z_yzz_xyyyyz_1[i] * pc_y[i];

        ta1_z_yyzz_xyyyzz_0[i] = ta1_z_zz_xyyyzz_0[i] * fe_0 - ta1_z_zz_xyyyzz_1[i] * fe_0 + 3.0 * ta1_z_yzz_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_yzz_xyyzz_1[i] * fe_0 + ta1_z_yzz_xyyyzz_0[i] * pa_y[i] - ta1_z_yzz_xyyyzz_1[i] * pc_y[i];

        ta1_z_yyzz_xyyzzz_0[i] = ta1_z_zz_xyyzzz_0[i] * fe_0 - ta1_z_zz_xyyzzz_1[i] * fe_0 + 2.0 * ta1_z_yzz_xyzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_yzz_xyzzz_1[i] * fe_0 + ta1_z_yzz_xyyzzz_0[i] * pa_y[i] - ta1_z_yzz_xyyzzz_1[i] * pc_y[i];

        ta1_z_yyzz_xyzzzz_0[i] = ta1_z_zz_xyzzzz_0[i] * fe_0 - ta1_z_zz_xyzzzz_1[i] * fe_0 + ta1_z_yzz_xzzzz_0[i] * fe_0 -
                                 ta1_z_yzz_xzzzz_1[i] * fe_0 + ta1_z_yzz_xyzzzz_0[i] * pa_y[i] - ta1_z_yzz_xyzzzz_1[i] * pc_y[i];

        ta1_z_yyzz_xzzzzz_0[i] =
            ta1_z_zz_xzzzzz_0[i] * fe_0 - ta1_z_zz_xzzzzz_1[i] * fe_0 + ta1_z_yzz_xzzzzz_0[i] * pa_y[i] - ta1_z_yzz_xzzzzz_1[i] * pc_y[i];

        ta1_z_yyzz_yyyyyy_0[i] = ta1_z_yy_yyyyyy_0[i] * fe_0 - ta1_z_yy_yyyyyy_1[i] * fe_0 + ta_yyz_yyyyyy_1[i] + ta1_z_yyz_yyyyyy_0[i] * pa_z[i] -
                                 ta1_z_yyz_yyyyyy_1[i] * pc_z[i];

        ta1_z_yyzz_yyyyyz_0[i] = ta1_z_zz_yyyyyz_0[i] * fe_0 - ta1_z_zz_yyyyyz_1[i] * fe_0 + 5.0 * ta1_z_yzz_yyyyz_0[i] * fe_0 -
                                 5.0 * ta1_z_yzz_yyyyz_1[i] * fe_0 + ta1_z_yzz_yyyyyz_0[i] * pa_y[i] - ta1_z_yzz_yyyyyz_1[i] * pc_y[i];

        ta1_z_yyzz_yyyyzz_0[i] = ta1_z_zz_yyyyzz_0[i] * fe_0 - ta1_z_zz_yyyyzz_1[i] * fe_0 + 4.0 * ta1_z_yzz_yyyzz_0[i] * fe_0 -
                                 4.0 * ta1_z_yzz_yyyzz_1[i] * fe_0 + ta1_z_yzz_yyyyzz_0[i] * pa_y[i] - ta1_z_yzz_yyyyzz_1[i] * pc_y[i];

        ta1_z_yyzz_yyyzzz_0[i] = ta1_z_zz_yyyzzz_0[i] * fe_0 - ta1_z_zz_yyyzzz_1[i] * fe_0 + 3.0 * ta1_z_yzz_yyzzz_0[i] * fe_0 -
                                 3.0 * ta1_z_yzz_yyzzz_1[i] * fe_0 + ta1_z_yzz_yyyzzz_0[i] * pa_y[i] - ta1_z_yzz_yyyzzz_1[i] * pc_y[i];

        ta1_z_yyzz_yyzzzz_0[i] = ta1_z_zz_yyzzzz_0[i] * fe_0 - ta1_z_zz_yyzzzz_1[i] * fe_0 + 2.0 * ta1_z_yzz_yzzzz_0[i] * fe_0 -
                                 2.0 * ta1_z_yzz_yzzzz_1[i] * fe_0 + ta1_z_yzz_yyzzzz_0[i] * pa_y[i] - ta1_z_yzz_yyzzzz_1[i] * pc_y[i];

        ta1_z_yyzz_yzzzzz_0[i] = ta1_z_zz_yzzzzz_0[i] * fe_0 - ta1_z_zz_yzzzzz_1[i] * fe_0 + ta1_z_yzz_zzzzz_0[i] * fe_0 -
                                 ta1_z_yzz_zzzzz_1[i] * fe_0 + ta1_z_yzz_yzzzzz_0[i] * pa_y[i] - ta1_z_yzz_yzzzzz_1[i] * pc_y[i];

        ta1_z_yyzz_zzzzzz_0[i] =
            ta1_z_zz_zzzzzz_0[i] * fe_0 - ta1_z_zz_zzzzzz_1[i] * fe_0 + ta1_z_yzz_zzzzzz_0[i] * pa_y[i] - ta1_z_yzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 1204-1232 components of targeted buffer : GI

    auto ta1_z_yzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1204);

    auto ta1_z_yzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1205);

    auto ta1_z_yzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1206);

    auto ta1_z_yzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1207);

    auto ta1_z_yzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1208);

    auto ta1_z_yzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1209);

    auto ta1_z_yzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1210);

    auto ta1_z_yzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1211);

    auto ta1_z_yzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1212);

    auto ta1_z_yzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1213);

    auto ta1_z_yzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1214);

    auto ta1_z_yzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1215);

    auto ta1_z_yzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1216);

    auto ta1_z_yzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1217);

    auto ta1_z_yzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1218);

    auto ta1_z_yzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1219);

    auto ta1_z_yzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1220);

    auto ta1_z_yzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1221);

    auto ta1_z_yzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1222);

    auto ta1_z_yzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1223);

    auto ta1_z_yzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1224);

    auto ta1_z_yzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1225);

    auto ta1_z_yzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1226);

    auto ta1_z_yzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1227);

    auto ta1_z_yzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1228);

    auto ta1_z_yzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1229);

    auto ta1_z_yzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1230);

    auto ta1_z_yzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1231);

#pragma omp simd aligned(pa_y,                    \
                             pc_y,                \
                             ta1_z_yzzz_xxxxxx_0, \
                             ta1_z_yzzz_xxxxxy_0, \
                             ta1_z_yzzz_xxxxxz_0, \
                             ta1_z_yzzz_xxxxyy_0, \
                             ta1_z_yzzz_xxxxyz_0, \
                             ta1_z_yzzz_xxxxzz_0, \
                             ta1_z_yzzz_xxxyyy_0, \
                             ta1_z_yzzz_xxxyyz_0, \
                             ta1_z_yzzz_xxxyzz_0, \
                             ta1_z_yzzz_xxxzzz_0, \
                             ta1_z_yzzz_xxyyyy_0, \
                             ta1_z_yzzz_xxyyyz_0, \
                             ta1_z_yzzz_xxyyzz_0, \
                             ta1_z_yzzz_xxyzzz_0, \
                             ta1_z_yzzz_xxzzzz_0, \
                             ta1_z_yzzz_xyyyyy_0, \
                             ta1_z_yzzz_xyyyyz_0, \
                             ta1_z_yzzz_xyyyzz_0, \
                             ta1_z_yzzz_xyyzzz_0, \
                             ta1_z_yzzz_xyzzzz_0, \
                             ta1_z_yzzz_xzzzzz_0, \
                             ta1_z_yzzz_yyyyyy_0, \
                             ta1_z_yzzz_yyyyyz_0, \
                             ta1_z_yzzz_yyyyzz_0, \
                             ta1_z_yzzz_yyyzzz_0, \
                             ta1_z_yzzz_yyzzzz_0, \
                             ta1_z_yzzz_yzzzzz_0, \
                             ta1_z_yzzz_zzzzzz_0, \
                             ta1_z_zzz_xxxxx_0,   \
                             ta1_z_zzz_xxxxx_1,   \
                             ta1_z_zzz_xxxxxx_0,  \
                             ta1_z_zzz_xxxxxx_1,  \
                             ta1_z_zzz_xxxxxy_0,  \
                             ta1_z_zzz_xxxxxy_1,  \
                             ta1_z_zzz_xxxxxz_0,  \
                             ta1_z_zzz_xxxxxz_1,  \
                             ta1_z_zzz_xxxxy_0,   \
                             ta1_z_zzz_xxxxy_1,   \
                             ta1_z_zzz_xxxxyy_0,  \
                             ta1_z_zzz_xxxxyy_1,  \
                             ta1_z_zzz_xxxxyz_0,  \
                             ta1_z_zzz_xxxxyz_1,  \
                             ta1_z_zzz_xxxxz_0,   \
                             ta1_z_zzz_xxxxz_1,   \
                             ta1_z_zzz_xxxxzz_0,  \
                             ta1_z_zzz_xxxxzz_1,  \
                             ta1_z_zzz_xxxyy_0,   \
                             ta1_z_zzz_xxxyy_1,   \
                             ta1_z_zzz_xxxyyy_0,  \
                             ta1_z_zzz_xxxyyy_1,  \
                             ta1_z_zzz_xxxyyz_0,  \
                             ta1_z_zzz_xxxyyz_1,  \
                             ta1_z_zzz_xxxyz_0,   \
                             ta1_z_zzz_xxxyz_1,   \
                             ta1_z_zzz_xxxyzz_0,  \
                             ta1_z_zzz_xxxyzz_1,  \
                             ta1_z_zzz_xxxzz_0,   \
                             ta1_z_zzz_xxxzz_1,   \
                             ta1_z_zzz_xxxzzz_0,  \
                             ta1_z_zzz_xxxzzz_1,  \
                             ta1_z_zzz_xxyyy_0,   \
                             ta1_z_zzz_xxyyy_1,   \
                             ta1_z_zzz_xxyyyy_0,  \
                             ta1_z_zzz_xxyyyy_1,  \
                             ta1_z_zzz_xxyyyz_0,  \
                             ta1_z_zzz_xxyyyz_1,  \
                             ta1_z_zzz_xxyyz_0,   \
                             ta1_z_zzz_xxyyz_1,   \
                             ta1_z_zzz_xxyyzz_0,  \
                             ta1_z_zzz_xxyyzz_1,  \
                             ta1_z_zzz_xxyzz_0,   \
                             ta1_z_zzz_xxyzz_1,   \
                             ta1_z_zzz_xxyzzz_0,  \
                             ta1_z_zzz_xxyzzz_1,  \
                             ta1_z_zzz_xxzzz_0,   \
                             ta1_z_zzz_xxzzz_1,   \
                             ta1_z_zzz_xxzzzz_0,  \
                             ta1_z_zzz_xxzzzz_1,  \
                             ta1_z_zzz_xyyyy_0,   \
                             ta1_z_zzz_xyyyy_1,   \
                             ta1_z_zzz_xyyyyy_0,  \
                             ta1_z_zzz_xyyyyy_1,  \
                             ta1_z_zzz_xyyyyz_0,  \
                             ta1_z_zzz_xyyyyz_1,  \
                             ta1_z_zzz_xyyyz_0,   \
                             ta1_z_zzz_xyyyz_1,   \
                             ta1_z_zzz_xyyyzz_0,  \
                             ta1_z_zzz_xyyyzz_1,  \
                             ta1_z_zzz_xyyzz_0,   \
                             ta1_z_zzz_xyyzz_1,   \
                             ta1_z_zzz_xyyzzz_0,  \
                             ta1_z_zzz_xyyzzz_1,  \
                             ta1_z_zzz_xyzzz_0,   \
                             ta1_z_zzz_xyzzz_1,   \
                             ta1_z_zzz_xyzzzz_0,  \
                             ta1_z_zzz_xyzzzz_1,  \
                             ta1_z_zzz_xzzzz_0,   \
                             ta1_z_zzz_xzzzz_1,   \
                             ta1_z_zzz_xzzzzz_0,  \
                             ta1_z_zzz_xzzzzz_1,  \
                             ta1_z_zzz_yyyyy_0,   \
                             ta1_z_zzz_yyyyy_1,   \
                             ta1_z_zzz_yyyyyy_0,  \
                             ta1_z_zzz_yyyyyy_1,  \
                             ta1_z_zzz_yyyyyz_0,  \
                             ta1_z_zzz_yyyyyz_1,  \
                             ta1_z_zzz_yyyyz_0,   \
                             ta1_z_zzz_yyyyz_1,   \
                             ta1_z_zzz_yyyyzz_0,  \
                             ta1_z_zzz_yyyyzz_1,  \
                             ta1_z_zzz_yyyzz_0,   \
                             ta1_z_zzz_yyyzz_1,   \
                             ta1_z_zzz_yyyzzz_0,  \
                             ta1_z_zzz_yyyzzz_1,  \
                             ta1_z_zzz_yyzzz_0,   \
                             ta1_z_zzz_yyzzz_1,   \
                             ta1_z_zzz_yyzzzz_0,  \
                             ta1_z_zzz_yyzzzz_1,  \
                             ta1_z_zzz_yzzzz_0,   \
                             ta1_z_zzz_yzzzz_1,   \
                             ta1_z_zzz_yzzzzz_0,  \
                             ta1_z_zzz_yzzzzz_1,  \
                             ta1_z_zzz_zzzzz_0,   \
                             ta1_z_zzz_zzzzz_1,   \
                             ta1_z_zzz_zzzzzz_0,  \
                             ta1_z_zzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzz_xxxxxx_0[i] = ta1_z_zzz_xxxxxx_0[i] * pa_y[i] - ta1_z_zzz_xxxxxx_1[i] * pc_y[i];

        ta1_z_yzzz_xxxxxy_0[i] =
            ta1_z_zzz_xxxxx_0[i] * fe_0 - ta1_z_zzz_xxxxx_1[i] * fe_0 + ta1_z_zzz_xxxxxy_0[i] * pa_y[i] - ta1_z_zzz_xxxxxy_1[i] * pc_y[i];

        ta1_z_yzzz_xxxxxz_0[i] = ta1_z_zzz_xxxxxz_0[i] * pa_y[i] - ta1_z_zzz_xxxxxz_1[i] * pc_y[i];

        ta1_z_yzzz_xxxxyy_0[i] =
            2.0 * ta1_z_zzz_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xxxxy_1[i] * fe_0 + ta1_z_zzz_xxxxyy_0[i] * pa_y[i] - ta1_z_zzz_xxxxyy_1[i] * pc_y[i];

        ta1_z_yzzz_xxxxyz_0[i] =
            ta1_z_zzz_xxxxz_0[i] * fe_0 - ta1_z_zzz_xxxxz_1[i] * fe_0 + ta1_z_zzz_xxxxyz_0[i] * pa_y[i] - ta1_z_zzz_xxxxyz_1[i] * pc_y[i];

        ta1_z_yzzz_xxxxzz_0[i] = ta1_z_zzz_xxxxzz_0[i] * pa_y[i] - ta1_z_zzz_xxxxzz_1[i] * pc_y[i];

        ta1_z_yzzz_xxxyyy_0[i] =
            3.0 * ta1_z_zzz_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxxyy_1[i] * fe_0 + ta1_z_zzz_xxxyyy_0[i] * pa_y[i] - ta1_z_zzz_xxxyyy_1[i] * pc_y[i];

        ta1_z_yzzz_xxxyyz_0[i] =
            2.0 * ta1_z_zzz_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xxxyz_1[i] * fe_0 + ta1_z_zzz_xxxyyz_0[i] * pa_y[i] - ta1_z_zzz_xxxyyz_1[i] * pc_y[i];

        ta1_z_yzzz_xxxyzz_0[i] =
            ta1_z_zzz_xxxzz_0[i] * fe_0 - ta1_z_zzz_xxxzz_1[i] * fe_0 + ta1_z_zzz_xxxyzz_0[i] * pa_y[i] - ta1_z_zzz_xxxyzz_1[i] * pc_y[i];

        ta1_z_yzzz_xxxzzz_0[i] = ta1_z_zzz_xxxzzz_0[i] * pa_y[i] - ta1_z_zzz_xxxzzz_1[i] * pc_y[i];

        ta1_z_yzzz_xxyyyy_0[i] =
            4.0 * ta1_z_zzz_xxyyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxyyy_1[i] * fe_0 + ta1_z_zzz_xxyyyy_0[i] * pa_y[i] - ta1_z_zzz_xxyyyy_1[i] * pc_y[i];

        ta1_z_yzzz_xxyyyz_0[i] =
            3.0 * ta1_z_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxyyz_1[i] * fe_0 + ta1_z_zzz_xxyyyz_0[i] * pa_y[i] - ta1_z_zzz_xxyyyz_1[i] * pc_y[i];

        ta1_z_yzzz_xxyyzz_0[i] =
            2.0 * ta1_z_zzz_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xxyzz_1[i] * fe_0 + ta1_z_zzz_xxyyzz_0[i] * pa_y[i] - ta1_z_zzz_xxyyzz_1[i] * pc_y[i];

        ta1_z_yzzz_xxyzzz_0[i] =
            ta1_z_zzz_xxzzz_0[i] * fe_0 - ta1_z_zzz_xxzzz_1[i] * fe_0 + ta1_z_zzz_xxyzzz_0[i] * pa_y[i] - ta1_z_zzz_xxyzzz_1[i] * pc_y[i];

        ta1_z_yzzz_xxzzzz_0[i] = ta1_z_zzz_xxzzzz_0[i] * pa_y[i] - ta1_z_zzz_xxzzzz_1[i] * pc_y[i];

        ta1_z_yzzz_xyyyyy_0[i] =
            5.0 * ta1_z_zzz_xyyyy_0[i] * fe_0 - 5.0 * ta1_z_zzz_xyyyy_1[i] * fe_0 + ta1_z_zzz_xyyyyy_0[i] * pa_y[i] - ta1_z_zzz_xyyyyy_1[i] * pc_y[i];

        ta1_z_yzzz_xyyyyz_0[i] =
            4.0 * ta1_z_zzz_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xyyyz_1[i] * fe_0 + ta1_z_zzz_xyyyyz_0[i] * pa_y[i] - ta1_z_zzz_xyyyyz_1[i] * pc_y[i];

        ta1_z_yzzz_xyyyzz_0[i] =
            3.0 * ta1_z_zzz_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xyyzz_1[i] * fe_0 + ta1_z_zzz_xyyyzz_0[i] * pa_y[i] - ta1_z_zzz_xyyyzz_1[i] * pc_y[i];

        ta1_z_yzzz_xyyzzz_0[i] =
            2.0 * ta1_z_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyzzz_1[i] * fe_0 + ta1_z_zzz_xyyzzz_0[i] * pa_y[i] - ta1_z_zzz_xyyzzz_1[i] * pc_y[i];

        ta1_z_yzzz_xyzzzz_0[i] =
            ta1_z_zzz_xzzzz_0[i] * fe_0 - ta1_z_zzz_xzzzz_1[i] * fe_0 + ta1_z_zzz_xyzzzz_0[i] * pa_y[i] - ta1_z_zzz_xyzzzz_1[i] * pc_y[i];

        ta1_z_yzzz_xzzzzz_0[i] = ta1_z_zzz_xzzzzz_0[i] * pa_y[i] - ta1_z_zzz_xzzzzz_1[i] * pc_y[i];

        ta1_z_yzzz_yyyyyy_0[i] =
            6.0 * ta1_z_zzz_yyyyy_0[i] * fe_0 - 6.0 * ta1_z_zzz_yyyyy_1[i] * fe_0 + ta1_z_zzz_yyyyyy_0[i] * pa_y[i] - ta1_z_zzz_yyyyyy_1[i] * pc_y[i];

        ta1_z_yzzz_yyyyyz_0[i] =
            5.0 * ta1_z_zzz_yyyyz_0[i] * fe_0 - 5.0 * ta1_z_zzz_yyyyz_1[i] * fe_0 + ta1_z_zzz_yyyyyz_0[i] * pa_y[i] - ta1_z_zzz_yyyyyz_1[i] * pc_y[i];

        ta1_z_yzzz_yyyyzz_0[i] =
            4.0 * ta1_z_zzz_yyyzz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyyzz_1[i] * fe_0 + ta1_z_zzz_yyyyzz_0[i] * pa_y[i] - ta1_z_zzz_yyyyzz_1[i] * pc_y[i];

        ta1_z_yzzz_yyyzzz_0[i] =
            3.0 * ta1_z_zzz_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_zzz_yyzzz_1[i] * fe_0 + ta1_z_zzz_yyyzzz_0[i] * pa_y[i] - ta1_z_zzz_yyyzzz_1[i] * pc_y[i];

        ta1_z_yzzz_yyzzzz_0[i] =
            2.0 * ta1_z_zzz_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_yzzzz_1[i] * fe_0 + ta1_z_zzz_yyzzzz_0[i] * pa_y[i] - ta1_z_zzz_yyzzzz_1[i] * pc_y[i];

        ta1_z_yzzz_yzzzzz_0[i] =
            ta1_z_zzz_zzzzz_0[i] * fe_0 - ta1_z_zzz_zzzzz_1[i] * fe_0 + ta1_z_zzz_yzzzzz_0[i] * pa_y[i] - ta1_z_zzz_yzzzzz_1[i] * pc_y[i];

        ta1_z_yzzz_zzzzzz_0[i] = ta1_z_zzz_zzzzzz_0[i] * pa_y[i] - ta1_z_zzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 1232-1260 components of targeted buffer : GI

    auto ta1_z_zzzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1232);

    auto ta1_z_zzzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1233);

    auto ta1_z_zzzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1234);

    auto ta1_z_zzzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1235);

    auto ta1_z_zzzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1236);

    auto ta1_z_zzzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1237);

    auto ta1_z_zzzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1238);

    auto ta1_z_zzzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1239);

    auto ta1_z_zzzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1240);

    auto ta1_z_zzzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1241);

    auto ta1_z_zzzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1242);

    auto ta1_z_zzzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1243);

    auto ta1_z_zzzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1244);

    auto ta1_z_zzzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1245);

    auto ta1_z_zzzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1246);

    auto ta1_z_zzzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1247);

    auto ta1_z_zzzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1248);

    auto ta1_z_zzzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1249);

    auto ta1_z_zzzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1250);

    auto ta1_z_zzzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1251);

    auto ta1_z_zzzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1252);

    auto ta1_z_zzzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1253);

    auto ta1_z_zzzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1254);

    auto ta1_z_zzzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1255);

    auto ta1_z_zzzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1256);

    auto ta1_z_zzzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1257);

    auto ta1_z_zzzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1258);

    auto ta1_z_zzzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gi + 1259);

#pragma omp simd aligned(pa_z,                    \
                             pc_z,                \
                             ta1_z_zz_xxxxxx_0,   \
                             ta1_z_zz_xxxxxx_1,   \
                             ta1_z_zz_xxxxxy_0,   \
                             ta1_z_zz_xxxxxy_1,   \
                             ta1_z_zz_xxxxxz_0,   \
                             ta1_z_zz_xxxxxz_1,   \
                             ta1_z_zz_xxxxyy_0,   \
                             ta1_z_zz_xxxxyy_1,   \
                             ta1_z_zz_xxxxyz_0,   \
                             ta1_z_zz_xxxxyz_1,   \
                             ta1_z_zz_xxxxzz_0,   \
                             ta1_z_zz_xxxxzz_1,   \
                             ta1_z_zz_xxxyyy_0,   \
                             ta1_z_zz_xxxyyy_1,   \
                             ta1_z_zz_xxxyyz_0,   \
                             ta1_z_zz_xxxyyz_1,   \
                             ta1_z_zz_xxxyzz_0,   \
                             ta1_z_zz_xxxyzz_1,   \
                             ta1_z_zz_xxxzzz_0,   \
                             ta1_z_zz_xxxzzz_1,   \
                             ta1_z_zz_xxyyyy_0,   \
                             ta1_z_zz_xxyyyy_1,   \
                             ta1_z_zz_xxyyyz_0,   \
                             ta1_z_zz_xxyyyz_1,   \
                             ta1_z_zz_xxyyzz_0,   \
                             ta1_z_zz_xxyyzz_1,   \
                             ta1_z_zz_xxyzzz_0,   \
                             ta1_z_zz_xxyzzz_1,   \
                             ta1_z_zz_xxzzzz_0,   \
                             ta1_z_zz_xxzzzz_1,   \
                             ta1_z_zz_xyyyyy_0,   \
                             ta1_z_zz_xyyyyy_1,   \
                             ta1_z_zz_xyyyyz_0,   \
                             ta1_z_zz_xyyyyz_1,   \
                             ta1_z_zz_xyyyzz_0,   \
                             ta1_z_zz_xyyyzz_1,   \
                             ta1_z_zz_xyyzzz_0,   \
                             ta1_z_zz_xyyzzz_1,   \
                             ta1_z_zz_xyzzzz_0,   \
                             ta1_z_zz_xyzzzz_1,   \
                             ta1_z_zz_xzzzzz_0,   \
                             ta1_z_zz_xzzzzz_1,   \
                             ta1_z_zz_yyyyyy_0,   \
                             ta1_z_zz_yyyyyy_1,   \
                             ta1_z_zz_yyyyyz_0,   \
                             ta1_z_zz_yyyyyz_1,   \
                             ta1_z_zz_yyyyzz_0,   \
                             ta1_z_zz_yyyyzz_1,   \
                             ta1_z_zz_yyyzzz_0,   \
                             ta1_z_zz_yyyzzz_1,   \
                             ta1_z_zz_yyzzzz_0,   \
                             ta1_z_zz_yyzzzz_1,   \
                             ta1_z_zz_yzzzzz_0,   \
                             ta1_z_zz_yzzzzz_1,   \
                             ta1_z_zz_zzzzzz_0,   \
                             ta1_z_zz_zzzzzz_1,   \
                             ta1_z_zzz_xxxxx_0,   \
                             ta1_z_zzz_xxxxx_1,   \
                             ta1_z_zzz_xxxxxx_0,  \
                             ta1_z_zzz_xxxxxx_1,  \
                             ta1_z_zzz_xxxxxy_0,  \
                             ta1_z_zzz_xxxxxy_1,  \
                             ta1_z_zzz_xxxxxz_0,  \
                             ta1_z_zzz_xxxxxz_1,  \
                             ta1_z_zzz_xxxxy_0,   \
                             ta1_z_zzz_xxxxy_1,   \
                             ta1_z_zzz_xxxxyy_0,  \
                             ta1_z_zzz_xxxxyy_1,  \
                             ta1_z_zzz_xxxxyz_0,  \
                             ta1_z_zzz_xxxxyz_1,  \
                             ta1_z_zzz_xxxxz_0,   \
                             ta1_z_zzz_xxxxz_1,   \
                             ta1_z_zzz_xxxxzz_0,  \
                             ta1_z_zzz_xxxxzz_1,  \
                             ta1_z_zzz_xxxyy_0,   \
                             ta1_z_zzz_xxxyy_1,   \
                             ta1_z_zzz_xxxyyy_0,  \
                             ta1_z_zzz_xxxyyy_1,  \
                             ta1_z_zzz_xxxyyz_0,  \
                             ta1_z_zzz_xxxyyz_1,  \
                             ta1_z_zzz_xxxyz_0,   \
                             ta1_z_zzz_xxxyz_1,   \
                             ta1_z_zzz_xxxyzz_0,  \
                             ta1_z_zzz_xxxyzz_1,  \
                             ta1_z_zzz_xxxzz_0,   \
                             ta1_z_zzz_xxxzz_1,   \
                             ta1_z_zzz_xxxzzz_0,  \
                             ta1_z_zzz_xxxzzz_1,  \
                             ta1_z_zzz_xxyyy_0,   \
                             ta1_z_zzz_xxyyy_1,   \
                             ta1_z_zzz_xxyyyy_0,  \
                             ta1_z_zzz_xxyyyy_1,  \
                             ta1_z_zzz_xxyyyz_0,  \
                             ta1_z_zzz_xxyyyz_1,  \
                             ta1_z_zzz_xxyyz_0,   \
                             ta1_z_zzz_xxyyz_1,   \
                             ta1_z_zzz_xxyyzz_0,  \
                             ta1_z_zzz_xxyyzz_1,  \
                             ta1_z_zzz_xxyzz_0,   \
                             ta1_z_zzz_xxyzz_1,   \
                             ta1_z_zzz_xxyzzz_0,  \
                             ta1_z_zzz_xxyzzz_1,  \
                             ta1_z_zzz_xxzzz_0,   \
                             ta1_z_zzz_xxzzz_1,   \
                             ta1_z_zzz_xxzzzz_0,  \
                             ta1_z_zzz_xxzzzz_1,  \
                             ta1_z_zzz_xyyyy_0,   \
                             ta1_z_zzz_xyyyy_1,   \
                             ta1_z_zzz_xyyyyy_0,  \
                             ta1_z_zzz_xyyyyy_1,  \
                             ta1_z_zzz_xyyyyz_0,  \
                             ta1_z_zzz_xyyyyz_1,  \
                             ta1_z_zzz_xyyyz_0,   \
                             ta1_z_zzz_xyyyz_1,   \
                             ta1_z_zzz_xyyyzz_0,  \
                             ta1_z_zzz_xyyyzz_1,  \
                             ta1_z_zzz_xyyzz_0,   \
                             ta1_z_zzz_xyyzz_1,   \
                             ta1_z_zzz_xyyzzz_0,  \
                             ta1_z_zzz_xyyzzz_1,  \
                             ta1_z_zzz_xyzzz_0,   \
                             ta1_z_zzz_xyzzz_1,   \
                             ta1_z_zzz_xyzzzz_0,  \
                             ta1_z_zzz_xyzzzz_1,  \
                             ta1_z_zzz_xzzzz_0,   \
                             ta1_z_zzz_xzzzz_1,   \
                             ta1_z_zzz_xzzzzz_0,  \
                             ta1_z_zzz_xzzzzz_1,  \
                             ta1_z_zzz_yyyyy_0,   \
                             ta1_z_zzz_yyyyy_1,   \
                             ta1_z_zzz_yyyyyy_0,  \
                             ta1_z_zzz_yyyyyy_1,  \
                             ta1_z_zzz_yyyyyz_0,  \
                             ta1_z_zzz_yyyyyz_1,  \
                             ta1_z_zzz_yyyyz_0,   \
                             ta1_z_zzz_yyyyz_1,   \
                             ta1_z_zzz_yyyyzz_0,  \
                             ta1_z_zzz_yyyyzz_1,  \
                             ta1_z_zzz_yyyzz_0,   \
                             ta1_z_zzz_yyyzz_1,   \
                             ta1_z_zzz_yyyzzz_0,  \
                             ta1_z_zzz_yyyzzz_1,  \
                             ta1_z_zzz_yyzzz_0,   \
                             ta1_z_zzz_yyzzz_1,   \
                             ta1_z_zzz_yyzzzz_0,  \
                             ta1_z_zzz_yyzzzz_1,  \
                             ta1_z_zzz_yzzzz_0,   \
                             ta1_z_zzz_yzzzz_1,   \
                             ta1_z_zzz_yzzzzz_0,  \
                             ta1_z_zzz_yzzzzz_1,  \
                             ta1_z_zzz_zzzzz_0,   \
                             ta1_z_zzz_zzzzz_1,   \
                             ta1_z_zzz_zzzzzz_0,  \
                             ta1_z_zzz_zzzzzz_1,  \
                             ta1_z_zzzz_xxxxxx_0, \
                             ta1_z_zzzz_xxxxxy_0, \
                             ta1_z_zzzz_xxxxxz_0, \
                             ta1_z_zzzz_xxxxyy_0, \
                             ta1_z_zzzz_xxxxyz_0, \
                             ta1_z_zzzz_xxxxzz_0, \
                             ta1_z_zzzz_xxxyyy_0, \
                             ta1_z_zzzz_xxxyyz_0, \
                             ta1_z_zzzz_xxxyzz_0, \
                             ta1_z_zzzz_xxxzzz_0, \
                             ta1_z_zzzz_xxyyyy_0, \
                             ta1_z_zzzz_xxyyyz_0, \
                             ta1_z_zzzz_xxyyzz_0, \
                             ta1_z_zzzz_xxyzzz_0, \
                             ta1_z_zzzz_xxzzzz_0, \
                             ta1_z_zzzz_xyyyyy_0, \
                             ta1_z_zzzz_xyyyyz_0, \
                             ta1_z_zzzz_xyyyzz_0, \
                             ta1_z_zzzz_xyyzzz_0, \
                             ta1_z_zzzz_xyzzzz_0, \
                             ta1_z_zzzz_xzzzzz_0, \
                             ta1_z_zzzz_yyyyyy_0, \
                             ta1_z_zzzz_yyyyyz_0, \
                             ta1_z_zzzz_yyyyzz_0, \
                             ta1_z_zzzz_yyyzzz_0, \
                             ta1_z_zzzz_yyzzzz_0, \
                             ta1_z_zzzz_yzzzzz_0, \
                             ta1_z_zzzz_zzzzzz_0, \
                             ta_zzz_xxxxxx_1,     \
                             ta_zzz_xxxxxy_1,     \
                             ta_zzz_xxxxxz_1,     \
                             ta_zzz_xxxxyy_1,     \
                             ta_zzz_xxxxyz_1,     \
                             ta_zzz_xxxxzz_1,     \
                             ta_zzz_xxxyyy_1,     \
                             ta_zzz_xxxyyz_1,     \
                             ta_zzz_xxxyzz_1,     \
                             ta_zzz_xxxzzz_1,     \
                             ta_zzz_xxyyyy_1,     \
                             ta_zzz_xxyyyz_1,     \
                             ta_zzz_xxyyzz_1,     \
                             ta_zzz_xxyzzz_1,     \
                             ta_zzz_xxzzzz_1,     \
                             ta_zzz_xyyyyy_1,     \
                             ta_zzz_xyyyyz_1,     \
                             ta_zzz_xyyyzz_1,     \
                             ta_zzz_xyyzzz_1,     \
                             ta_zzz_xyzzzz_1,     \
                             ta_zzz_xzzzzz_1,     \
                             ta_zzz_yyyyyy_1,     \
                             ta_zzz_yyyyyz_1,     \
                             ta_zzz_yyyyzz_1,     \
                             ta_zzz_yyyzzz_1,     \
                             ta_zzz_yyzzzz_1,     \
                             ta_zzz_yzzzzz_1,     \
                             ta_zzz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzz_xxxxxx_0[i] = 3.0 * ta1_z_zz_xxxxxx_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxxx_1[i] * fe_0 + ta_zzz_xxxxxx_1[i] +
                                 ta1_z_zzz_xxxxxx_0[i] * pa_z[i] - ta1_z_zzz_xxxxxx_1[i] * pc_z[i];

        ta1_z_zzzz_xxxxxy_0[i] = 3.0 * ta1_z_zz_xxxxxy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxxy_1[i] * fe_0 + ta_zzz_xxxxxy_1[i] +
                                 ta1_z_zzz_xxxxxy_0[i] * pa_z[i] - ta1_z_zzz_xxxxxy_1[i] * pc_z[i];

        ta1_z_zzzz_xxxxxz_0[i] = 3.0 * ta1_z_zz_xxxxxz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxxz_1[i] * fe_0 + ta1_z_zzz_xxxxx_0[i] * fe_0 -
                                 ta1_z_zzz_xxxxx_1[i] * fe_0 + ta_zzz_xxxxxz_1[i] + ta1_z_zzz_xxxxxz_0[i] * pa_z[i] - ta1_z_zzz_xxxxxz_1[i] * pc_z[i];

        ta1_z_zzzz_xxxxyy_0[i] = 3.0 * ta1_z_zz_xxxxyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxyy_1[i] * fe_0 + ta_zzz_xxxxyy_1[i] +
                                 ta1_z_zzz_xxxxyy_0[i] * pa_z[i] - ta1_z_zzz_xxxxyy_1[i] * pc_z[i];

        ta1_z_zzzz_xxxxyz_0[i] = 3.0 * ta1_z_zz_xxxxyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxyz_1[i] * fe_0 + ta1_z_zzz_xxxxy_0[i] * fe_0 -
                                 ta1_z_zzz_xxxxy_1[i] * fe_0 + ta_zzz_xxxxyz_1[i] + ta1_z_zzz_xxxxyz_0[i] * pa_z[i] - ta1_z_zzz_xxxxyz_1[i] * pc_z[i];

        ta1_z_zzzz_xxxxzz_0[i] = 3.0 * ta1_z_zz_xxxxzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxxxz_0[i] * fe_0 -
                                 2.0 * ta1_z_zzz_xxxxz_1[i] * fe_0 + ta_zzz_xxxxzz_1[i] + ta1_z_zzz_xxxxzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xxxxzz_1[i] * pc_z[i];

        ta1_z_zzzz_xxxyyy_0[i] = 3.0 * ta1_z_zz_xxxyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxyyy_1[i] * fe_0 + ta_zzz_xxxyyy_1[i] +
                                 ta1_z_zzz_xxxyyy_0[i] * pa_z[i] - ta1_z_zzz_xxxyyy_1[i] * pc_z[i];

        ta1_z_zzzz_xxxyyz_0[i] = 3.0 * ta1_z_zz_xxxyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxyyz_1[i] * fe_0 + ta1_z_zzz_xxxyy_0[i] * fe_0 -
                                 ta1_z_zzz_xxxyy_1[i] * fe_0 + ta_zzz_xxxyyz_1[i] + ta1_z_zzz_xxxyyz_0[i] * pa_z[i] - ta1_z_zzz_xxxyyz_1[i] * pc_z[i];

        ta1_z_zzzz_xxxyzz_0[i] = 3.0 * ta1_z_zz_xxxyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxxyz_0[i] * fe_0 -
                                 2.0 * ta1_z_zzz_xxxyz_1[i] * fe_0 + ta_zzz_xxxyzz_1[i] + ta1_z_zzz_xxxyzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xxxyzz_1[i] * pc_z[i];

        ta1_z_zzzz_xxxzzz_0[i] = 3.0 * ta1_z_zz_xxxzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_xxxzz_0[i] * fe_0 -
                                 3.0 * ta1_z_zzz_xxxzz_1[i] * fe_0 + ta_zzz_xxxzzz_1[i] + ta1_z_zzz_xxxzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xxxzzz_1[i] * pc_z[i];

        ta1_z_zzzz_xxyyyy_0[i] = 3.0 * ta1_z_zz_xxyyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyyy_1[i] * fe_0 + ta_zzz_xxyyyy_1[i] +
                                 ta1_z_zzz_xxyyyy_0[i] * pa_z[i] - ta1_z_zzz_xxyyyy_1[i] * pc_z[i];

        ta1_z_zzzz_xxyyyz_0[i] = 3.0 * ta1_z_zz_xxyyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyyz_1[i] * fe_0 + ta1_z_zzz_xxyyy_0[i] * fe_0 -
                                 ta1_z_zzz_xxyyy_1[i] * fe_0 + ta_zzz_xxyyyz_1[i] + ta1_z_zzz_xxyyyz_0[i] * pa_z[i] - ta1_z_zzz_xxyyyz_1[i] * pc_z[i];

        ta1_z_zzzz_xxyyzz_0[i] = 3.0 * ta1_z_zz_xxyyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxyyz_0[i] * fe_0 -
                                 2.0 * ta1_z_zzz_xxyyz_1[i] * fe_0 + ta_zzz_xxyyzz_1[i] + ta1_z_zzz_xxyyzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xxyyzz_1[i] * pc_z[i];

        ta1_z_zzzz_xxyzzz_0[i] = 3.0 * ta1_z_zz_xxyzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_xxyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_zzz_xxyzz_1[i] * fe_0 + ta_zzz_xxyzzz_1[i] + ta1_z_zzz_xxyzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xxyzzz_1[i] * pc_z[i];

        ta1_z_zzzz_xxzzzz_0[i] = 3.0 * ta1_z_zz_xxzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxzzzz_1[i] * fe_0 + 4.0 * ta1_z_zzz_xxzzz_0[i] * fe_0 -
                                 4.0 * ta1_z_zzz_xxzzz_1[i] * fe_0 + ta_zzz_xxzzzz_1[i] + ta1_z_zzz_xxzzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xxzzzz_1[i] * pc_z[i];

        ta1_z_zzzz_xyyyyy_0[i] = 3.0 * ta1_z_zz_xyyyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyyyy_1[i] * fe_0 + ta_zzz_xyyyyy_1[i] +
                                 ta1_z_zzz_xyyyyy_0[i] * pa_z[i] - ta1_z_zzz_xyyyyy_1[i] * pc_z[i];

        ta1_z_zzzz_xyyyyz_0[i] = 3.0 * ta1_z_zz_xyyyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyyyz_1[i] * fe_0 + ta1_z_zzz_xyyyy_0[i] * fe_0 -
                                 ta1_z_zzz_xyyyy_1[i] * fe_0 + ta_zzz_xyyyyz_1[i] + ta1_z_zzz_xyyyyz_0[i] * pa_z[i] - ta1_z_zzz_xyyyyz_1[i] * pc_z[i];

        ta1_z_zzzz_xyyyzz_0[i] = 3.0 * ta1_z_zz_xyyyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyyyz_0[i] * fe_0 -
                                 2.0 * ta1_z_zzz_xyyyz_1[i] * fe_0 + ta_zzz_xyyyzz_1[i] + ta1_z_zzz_xyyyzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xyyyzz_1[i] * pc_z[i];

        ta1_z_zzzz_xyyzzz_0[i] = 3.0 * ta1_z_zz_xyyzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_xyyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_zzz_xyyzz_1[i] * fe_0 + ta_zzz_xyyzzz_1[i] + ta1_z_zzz_xyyzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xyyzzz_1[i] * pc_z[i];

        ta1_z_zzzz_xyzzzz_0[i] = 3.0 * ta1_z_zz_xyzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyzzzz_1[i] * fe_0 + 4.0 * ta1_z_zzz_xyzzz_0[i] * fe_0 -
                                 4.0 * ta1_z_zzz_xyzzz_1[i] * fe_0 + ta_zzz_xyzzzz_1[i] + ta1_z_zzz_xyzzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xyzzzz_1[i] * pc_z[i];

        ta1_z_zzzz_xzzzzz_0[i] = 3.0 * ta1_z_zz_xzzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xzzzzz_1[i] * fe_0 + 5.0 * ta1_z_zzz_xzzzz_0[i] * fe_0 -
                                 5.0 * ta1_z_zzz_xzzzz_1[i] * fe_0 + ta_zzz_xzzzzz_1[i] + ta1_z_zzz_xzzzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_xzzzzz_1[i] * pc_z[i];

        ta1_z_zzzz_yyyyyy_0[i] = 3.0 * ta1_z_zz_yyyyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyyyy_1[i] * fe_0 + ta_zzz_yyyyyy_1[i] +
                                 ta1_z_zzz_yyyyyy_0[i] * pa_z[i] - ta1_z_zzz_yyyyyy_1[i] * pc_z[i];

        ta1_z_zzzz_yyyyyz_0[i] = 3.0 * ta1_z_zz_yyyyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyyyz_1[i] * fe_0 + ta1_z_zzz_yyyyy_0[i] * fe_0 -
                                 ta1_z_zzz_yyyyy_1[i] * fe_0 + ta_zzz_yyyyyz_1[i] + ta1_z_zzz_yyyyyz_0[i] * pa_z[i] - ta1_z_zzz_yyyyyz_1[i] * pc_z[i];

        ta1_z_zzzz_yyyyzz_0[i] = 3.0 * ta1_z_zz_yyyyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyyyz_0[i] * fe_0 -
                                 2.0 * ta1_z_zzz_yyyyz_1[i] * fe_0 + ta_zzz_yyyyzz_1[i] + ta1_z_zzz_yyyyzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_yyyyzz_1[i] * pc_z[i];

        ta1_z_zzzz_yyyzzz_0[i] = 3.0 * ta1_z_zz_yyyzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_yyyzz_0[i] * fe_0 -
                                 3.0 * ta1_z_zzz_yyyzz_1[i] * fe_0 + ta_zzz_yyyzzz_1[i] + ta1_z_zzz_yyyzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_yyyzzz_1[i] * pc_z[i];

        ta1_z_zzzz_yyzzzz_0[i] = 3.0 * ta1_z_zz_yyzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyzzzz_1[i] * fe_0 + 4.0 * ta1_z_zzz_yyzzz_0[i] * fe_0 -
                                 4.0 * ta1_z_zzz_yyzzz_1[i] * fe_0 + ta_zzz_yyzzzz_1[i] + ta1_z_zzz_yyzzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_yyzzzz_1[i] * pc_z[i];

        ta1_z_zzzz_yzzzzz_0[i] = 3.0 * ta1_z_zz_yzzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yzzzzz_1[i] * fe_0 + 5.0 * ta1_z_zzz_yzzzz_0[i] * fe_0 -
                                 5.0 * ta1_z_zzz_yzzzz_1[i] * fe_0 + ta_zzz_yzzzzz_1[i] + ta1_z_zzz_yzzzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_yzzzzz_1[i] * pc_z[i];

        ta1_z_zzzz_zzzzzz_0[i] = 3.0 * ta1_z_zz_zzzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_zzzzzz_1[i] * fe_0 + 6.0 * ta1_z_zzz_zzzzz_0[i] * fe_0 -
                                 6.0 * ta1_z_zzz_zzzzz_1[i] * fe_0 + ta_zzz_zzzzzz_1[i] + ta1_z_zzz_zzzzzz_0[i] * pa_z[i] -
                                 ta1_z_zzz_zzzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
