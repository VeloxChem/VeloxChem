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

#include "ElectronRepulsionPrimRecSGSL.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sgsl(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sgsl,
                                  size_t                idx_eri_0_sdsl,
                                  size_t                idx_eri_1_sdsl,
                                  size_t                idx_eri_1_sfsk,
                                  size_t                idx_eri_0_sfsl,
                                  size_t                idx_eri_1_sfsl,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SDSL

    auto g_0_xx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sdsl);

    auto g_0_xx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sdsl + 1);

    auto g_0_xx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sdsl + 2);

    auto g_0_xx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sdsl + 3);

    auto g_0_xx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sdsl + 4);

    auto g_0_xx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sdsl + 5);

    auto g_0_xx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sdsl + 6);

    auto g_0_xx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sdsl + 7);

    auto g_0_xx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sdsl + 8);

    auto g_0_xx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sdsl + 9);

    auto g_0_xx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 10);

    auto g_0_xx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 11);

    auto g_0_xx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 12);

    auto g_0_xx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 13);

    auto g_0_xx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 14);

    auto g_0_xx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 15);

    auto g_0_xx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 16);

    auto g_0_xx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 17);

    auto g_0_xx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 18);

    auto g_0_xx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 19);

    auto g_0_xx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 20);

    auto g_0_xx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 21);

    auto g_0_xx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 22);

    auto g_0_xx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 23);

    auto g_0_xx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 24);

    auto g_0_xx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 25);

    auto g_0_xx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 26);

    auto g_0_xx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 27);

    auto g_0_xx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 28);

    auto g_0_xx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 29);

    auto g_0_xx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 30);

    auto g_0_xx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 31);

    auto g_0_xx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 32);

    auto g_0_xx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 33);

    auto g_0_xx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 34);

    auto g_0_xx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 35);

    auto g_0_xx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 36);

    auto g_0_xx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 37);

    auto g_0_xx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 38);

    auto g_0_xx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 39);

    auto g_0_xx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 40);

    auto g_0_xx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 41);

    auto g_0_xx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 42);

    auto g_0_xx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 43);

    auto g_0_xx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 44);

    auto g_0_yy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sdsl + 135);

    auto g_0_yy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sdsl + 136);

    auto g_0_yy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sdsl + 137);

    auto g_0_yy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sdsl + 138);

    auto g_0_yy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sdsl + 139);

    auto g_0_yy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sdsl + 140);

    auto g_0_yy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sdsl + 141);

    auto g_0_yy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sdsl + 142);

    auto g_0_yy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sdsl + 143);

    auto g_0_yy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sdsl + 144);

    auto g_0_yy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 145);

    auto g_0_yy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 146);

    auto g_0_yy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 147);

    auto g_0_yy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 148);

    auto g_0_yy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 149);

    auto g_0_yy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 150);

    auto g_0_yy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 151);

    auto g_0_yy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 152);

    auto g_0_yy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 153);

    auto g_0_yy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 154);

    auto g_0_yy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 155);

    auto g_0_yy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 156);

    auto g_0_yy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 157);

    auto g_0_yy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 158);

    auto g_0_yy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 159);

    auto g_0_yy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 160);

    auto g_0_yy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 161);

    auto g_0_yy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 162);

    auto g_0_yy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 163);

    auto g_0_yy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 164);

    auto g_0_yy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 165);

    auto g_0_yy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 166);

    auto g_0_yy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 167);

    auto g_0_yy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 168);

    auto g_0_yy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 169);

    auto g_0_yy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 170);

    auto g_0_yy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 171);

    auto g_0_yy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 172);

    auto g_0_yy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 173);

    auto g_0_yy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 174);

    auto g_0_yy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 175);

    auto g_0_yy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 176);

    auto g_0_yy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 177);

    auto g_0_yy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 178);

    auto g_0_yy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 179);

    auto g_0_zz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sdsl + 225);

    auto g_0_zz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sdsl + 226);

    auto g_0_zz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sdsl + 227);

    auto g_0_zz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sdsl + 228);

    auto g_0_zz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sdsl + 229);

    auto g_0_zz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sdsl + 230);

    auto g_0_zz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sdsl + 231);

    auto g_0_zz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sdsl + 232);

    auto g_0_zz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sdsl + 233);

    auto g_0_zz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sdsl + 234);

    auto g_0_zz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 235);

    auto g_0_zz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 236);

    auto g_0_zz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 237);

    auto g_0_zz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 238);

    auto g_0_zz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 239);

    auto g_0_zz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 240);

    auto g_0_zz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 241);

    auto g_0_zz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 242);

    auto g_0_zz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 243);

    auto g_0_zz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 244);

    auto g_0_zz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 245);

    auto g_0_zz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 246);

    auto g_0_zz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 247);

    auto g_0_zz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 248);

    auto g_0_zz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 249);

    auto g_0_zz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 250);

    auto g_0_zz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 251);

    auto g_0_zz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 252);

    auto g_0_zz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 253);

    auto g_0_zz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 254);

    auto g_0_zz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 255);

    auto g_0_zz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 256);

    auto g_0_zz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 257);

    auto g_0_zz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 258);

    auto g_0_zz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 259);

    auto g_0_zz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 260);

    auto g_0_zz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 261);

    auto g_0_zz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 262);

    auto g_0_zz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 263);

    auto g_0_zz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 264);

    auto g_0_zz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 265);

    auto g_0_zz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 266);

    auto g_0_zz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 267);

    auto g_0_zz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 268);

    auto g_0_zz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 269);

    /// Set up components of auxilary buffer : SDSL

    auto g_0_xx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sdsl);

    auto g_0_xx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sdsl + 1);

    auto g_0_xx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sdsl + 2);

    auto g_0_xx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sdsl + 3);

    auto g_0_xx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sdsl + 4);

    auto g_0_xx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sdsl + 5);

    auto g_0_xx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sdsl + 6);

    auto g_0_xx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sdsl + 7);

    auto g_0_xx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sdsl + 8);

    auto g_0_xx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sdsl + 9);

    auto g_0_xx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 10);

    auto g_0_xx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 11);

    auto g_0_xx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 12);

    auto g_0_xx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 13);

    auto g_0_xx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 14);

    auto g_0_xx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 15);

    auto g_0_xx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 16);

    auto g_0_xx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 17);

    auto g_0_xx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 18);

    auto g_0_xx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 19);

    auto g_0_xx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 20);

    auto g_0_xx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 21);

    auto g_0_xx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 22);

    auto g_0_xx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 23);

    auto g_0_xx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 24);

    auto g_0_xx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 25);

    auto g_0_xx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 26);

    auto g_0_xx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 27);

    auto g_0_xx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 28);

    auto g_0_xx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 29);

    auto g_0_xx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 30);

    auto g_0_xx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 31);

    auto g_0_xx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 32);

    auto g_0_xx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 33);

    auto g_0_xx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 34);

    auto g_0_xx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 35);

    auto g_0_xx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 36);

    auto g_0_xx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 37);

    auto g_0_xx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 38);

    auto g_0_xx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 39);

    auto g_0_xx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 40);

    auto g_0_xx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 41);

    auto g_0_xx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 42);

    auto g_0_xx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 43);

    auto g_0_xx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 44);

    auto g_0_yy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sdsl + 135);

    auto g_0_yy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sdsl + 136);

    auto g_0_yy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sdsl + 137);

    auto g_0_yy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sdsl + 138);

    auto g_0_yy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sdsl + 139);

    auto g_0_yy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sdsl + 140);

    auto g_0_yy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sdsl + 141);

    auto g_0_yy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sdsl + 142);

    auto g_0_yy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sdsl + 143);

    auto g_0_yy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sdsl + 144);

    auto g_0_yy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 145);

    auto g_0_yy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 146);

    auto g_0_yy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 147);

    auto g_0_yy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 148);

    auto g_0_yy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 149);

    auto g_0_yy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 150);

    auto g_0_yy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 151);

    auto g_0_yy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 152);

    auto g_0_yy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 153);

    auto g_0_yy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 154);

    auto g_0_yy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 155);

    auto g_0_yy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 156);

    auto g_0_yy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 157);

    auto g_0_yy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 158);

    auto g_0_yy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 159);

    auto g_0_yy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 160);

    auto g_0_yy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 161);

    auto g_0_yy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 162);

    auto g_0_yy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 163);

    auto g_0_yy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 164);

    auto g_0_yy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 165);

    auto g_0_yy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 166);

    auto g_0_yy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 167);

    auto g_0_yy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 168);

    auto g_0_yy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 169);

    auto g_0_yy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 170);

    auto g_0_yy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 171);

    auto g_0_yy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 172);

    auto g_0_yy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 173);

    auto g_0_yy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 174);

    auto g_0_yy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 175);

    auto g_0_yy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 176);

    auto g_0_yy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 177);

    auto g_0_yy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 178);

    auto g_0_yy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 179);

    auto g_0_zz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sdsl + 225);

    auto g_0_zz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sdsl + 226);

    auto g_0_zz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sdsl + 227);

    auto g_0_zz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sdsl + 228);

    auto g_0_zz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sdsl + 229);

    auto g_0_zz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sdsl + 230);

    auto g_0_zz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sdsl + 231);

    auto g_0_zz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sdsl + 232);

    auto g_0_zz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sdsl + 233);

    auto g_0_zz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sdsl + 234);

    auto g_0_zz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 235);

    auto g_0_zz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 236);

    auto g_0_zz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 237);

    auto g_0_zz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 238);

    auto g_0_zz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 239);

    auto g_0_zz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 240);

    auto g_0_zz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 241);

    auto g_0_zz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 242);

    auto g_0_zz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 243);

    auto g_0_zz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 244);

    auto g_0_zz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 245);

    auto g_0_zz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 246);

    auto g_0_zz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 247);

    auto g_0_zz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 248);

    auto g_0_zz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 249);

    auto g_0_zz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 250);

    auto g_0_zz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 251);

    auto g_0_zz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 252);

    auto g_0_zz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 253);

    auto g_0_zz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 254);

    auto g_0_zz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 255);

    auto g_0_zz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 256);

    auto g_0_zz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 257);

    auto g_0_zz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 258);

    auto g_0_zz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 259);

    auto g_0_zz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 260);

    auto g_0_zz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 261);

    auto g_0_zz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 262);

    auto g_0_zz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 263);

    auto g_0_zz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 264);

    auto g_0_zz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 265);

    auto g_0_zz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 266);

    auto g_0_zz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 267);

    auto g_0_zz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 268);

    auto g_0_zz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 269);

    /// Set up components of auxilary buffer : SFSK

    auto g_0_xxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk);

    auto g_0_xxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 1);

    auto g_0_xxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 2);

    auto g_0_xxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 3);

    auto g_0_xxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 4);

    auto g_0_xxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 5);

    auto g_0_xxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 6);

    auto g_0_xxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 7);

    auto g_0_xxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 8);

    auto g_0_xxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 9);

    auto g_0_xxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 10);

    auto g_0_xxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 11);

    auto g_0_xxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 12);

    auto g_0_xxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 13);

    auto g_0_xxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 14);

    auto g_0_xxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 15);

    auto g_0_xxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 16);

    auto g_0_xxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 17);

    auto g_0_xxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 18);

    auto g_0_xxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 19);

    auto g_0_xxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 20);

    auto g_0_xxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 21);

    auto g_0_xxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 22);

    auto g_0_xxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 23);

    auto g_0_xxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 24);

    auto g_0_xxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 25);

    auto g_0_xxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 26);

    auto g_0_xxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 27);

    auto g_0_xxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 28);

    auto g_0_xxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 29);

    auto g_0_xxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 30);

    auto g_0_xxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 31);

    auto g_0_xxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 32);

    auto g_0_xxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 33);

    auto g_0_xxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 34);

    auto g_0_xxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 35);

    auto g_0_xxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 74);

    auto g_0_xxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 76);

    auto g_0_xxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 77);

    auto g_0_xxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 79);

    auto g_0_xxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 80);

    auto g_0_xxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 81);

    auto g_0_xxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 83);

    auto g_0_xxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 84);

    auto g_0_xxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 85);

    auto g_0_xxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 86);

    auto g_0_xxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 88);

    auto g_0_xxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 89);

    auto g_0_xxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 90);

    auto g_0_xxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 91);

    auto g_0_xxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 92);

    auto g_0_xxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 94);

    auto g_0_xxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 95);

    auto g_0_xxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 96);

    auto g_0_xxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 97);

    auto g_0_xxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 98);

    auto g_0_xxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 99);

    auto g_0_xxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 101);

    auto g_0_xxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 102);

    auto g_0_xxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 103);

    auto g_0_xxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 104);

    auto g_0_xxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 105);

    auto g_0_xxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 106);

    auto g_0_xxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 107);

    auto g_0_xyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 109);

    auto g_0_xyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 111);

    auto g_0_xyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 112);

    auto g_0_xyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 114);

    auto g_0_xyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 115);

    auto g_0_xyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 116);

    auto g_0_xyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 118);

    auto g_0_xyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 119);

    auto g_0_xyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 120);

    auto g_0_xyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 121);

    auto g_0_xyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 123);

    auto g_0_xyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 124);

    auto g_0_xyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 125);

    auto g_0_xyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 126);

    auto g_0_xyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 127);

    auto g_0_xyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 129);

    auto g_0_xyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 130);

    auto g_0_xyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 131);

    auto g_0_xyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 132);

    auto g_0_xyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 133);

    auto g_0_xyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 134);

    auto g_0_xyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 136);

    auto g_0_xyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 137);

    auto g_0_xyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 138);

    auto g_0_xyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 139);

    auto g_0_xyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 140);

    auto g_0_xyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 141);

    auto g_0_xyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 142);

    auto g_0_xzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 182);

    auto g_0_xzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 184);

    auto g_0_xzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 185);

    auto g_0_xzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 187);

    auto g_0_xzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 188);

    auto g_0_xzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 189);

    auto g_0_xzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 191);

    auto g_0_xzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 192);

    auto g_0_xzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 193);

    auto g_0_xzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 194);

    auto g_0_xzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 196);

    auto g_0_xzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 197);

    auto g_0_xzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 198);

    auto g_0_xzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 199);

    auto g_0_xzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 200);

    auto g_0_xzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 202);

    auto g_0_xzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 203);

    auto g_0_xzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 204);

    auto g_0_xzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 205);

    auto g_0_xzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 206);

    auto g_0_xzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 207);

    auto g_0_xzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 209);

    auto g_0_xzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 210);

    auto g_0_xzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 211);

    auto g_0_xzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 212);

    auto g_0_xzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 213);

    auto g_0_xzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 214);

    auto g_0_xzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 215);

    auto g_0_yyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk + 216);

    auto g_0_yyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 217);

    auto g_0_yyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 218);

    auto g_0_yyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 219);

    auto g_0_yyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 220);

    auto g_0_yyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 221);

    auto g_0_yyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 222);

    auto g_0_yyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 223);

    auto g_0_yyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 224);

    auto g_0_yyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 225);

    auto g_0_yyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 226);

    auto g_0_yyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 227);

    auto g_0_yyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 228);

    auto g_0_yyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 229);

    auto g_0_yyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 230);

    auto g_0_yyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 231);

    auto g_0_yyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 232);

    auto g_0_yyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 233);

    auto g_0_yyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 234);

    auto g_0_yyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 235);

    auto g_0_yyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 236);

    auto g_0_yyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 237);

    auto g_0_yyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 238);

    auto g_0_yyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 239);

    auto g_0_yyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 240);

    auto g_0_yyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 241);

    auto g_0_yyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 242);

    auto g_0_yyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 243);

    auto g_0_yyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 244);

    auto g_0_yyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 245);

    auto g_0_yyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 246);

    auto g_0_yyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 247);

    auto g_0_yyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 248);

    auto g_0_yyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 249);

    auto g_0_yyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 250);

    auto g_0_yyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 251);

    auto g_0_yyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 254);

    auto g_0_yyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 256);

    auto g_0_yyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 257);

    auto g_0_yyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 259);

    auto g_0_yyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 260);

    auto g_0_yyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 261);

    auto g_0_yyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 263);

    auto g_0_yyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 264);

    auto g_0_yyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 265);

    auto g_0_yyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 266);

    auto g_0_yyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 268);

    auto g_0_yyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 269);

    auto g_0_yyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 270);

    auto g_0_yyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 271);

    auto g_0_yyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 272);

    auto g_0_yyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 274);

    auto g_0_yyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 275);

    auto g_0_yyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 276);

    auto g_0_yyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 277);

    auto g_0_yyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 278);

    auto g_0_yyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 279);

    auto g_0_yyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 281);

    auto g_0_yyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 282);

    auto g_0_yyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 283);

    auto g_0_yyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 284);

    auto g_0_yyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 285);

    auto g_0_yyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 286);

    auto g_0_yyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 287);

    auto g_0_yzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 289);

    auto g_0_yzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 290);

    auto g_0_yzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 291);

    auto g_0_yzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 292);

    auto g_0_yzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 293);

    auto g_0_yzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 294);

    auto g_0_yzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 295);

    auto g_0_yzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 296);

    auto g_0_yzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 297);

    auto g_0_yzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 298);

    auto g_0_yzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 299);

    auto g_0_yzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 300);

    auto g_0_yzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 301);

    auto g_0_yzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 302);

    auto g_0_yzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 303);

    auto g_0_yzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 304);

    auto g_0_yzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 305);

    auto g_0_yzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 306);

    auto g_0_yzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 307);

    auto g_0_yzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 308);

    auto g_0_yzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 309);

    auto g_0_yzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 310);

    auto g_0_yzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 311);

    auto g_0_yzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 312);

    auto g_0_yzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 313);

    auto g_0_yzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 314);

    auto g_0_yzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 315);

    auto g_0_yzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 316);

    auto g_0_yzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 317);

    auto g_0_yzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 318);

    auto g_0_yzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 319);

    auto g_0_yzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 320);

    auto g_0_yzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 321);

    auto g_0_yzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 322);

    auto g_0_yzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 323);

    auto g_0_zzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk + 324);

    auto g_0_zzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 325);

    auto g_0_zzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 326);

    auto g_0_zzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 327);

    auto g_0_zzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 328);

    auto g_0_zzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 329);

    auto g_0_zzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 330);

    auto g_0_zzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 331);

    auto g_0_zzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 332);

    auto g_0_zzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 333);

    auto g_0_zzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 334);

    auto g_0_zzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 335);

    auto g_0_zzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 336);

    auto g_0_zzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 337);

    auto g_0_zzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 338);

    auto g_0_zzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 339);

    auto g_0_zzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 340);

    auto g_0_zzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 341);

    auto g_0_zzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 342);

    auto g_0_zzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 343);

    auto g_0_zzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 344);

    auto g_0_zzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 345);

    auto g_0_zzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 346);

    auto g_0_zzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 347);

    auto g_0_zzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 348);

    auto g_0_zzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 349);

    auto g_0_zzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 350);

    auto g_0_zzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 351);

    auto g_0_zzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 352);

    auto g_0_zzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 353);

    auto g_0_zzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 354);

    auto g_0_zzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 355);

    auto g_0_zzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 356);

    auto g_0_zzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 357);

    auto g_0_zzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 358);

    auto g_0_zzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 359);

    /// Set up components of auxilary buffer : SFSL

    auto g_0_xxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl);

    auto g_0_xxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 1);

    auto g_0_xxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 2);

    auto g_0_xxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 3);

    auto g_0_xxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 4);

    auto g_0_xxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 5);

    auto g_0_xxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 6);

    auto g_0_xxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 7);

    auto g_0_xxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 8);

    auto g_0_xxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 9);

    auto g_0_xxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 10);

    auto g_0_xxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 11);

    auto g_0_xxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 12);

    auto g_0_xxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 13);

    auto g_0_xxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 14);

    auto g_0_xxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 15);

    auto g_0_xxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 16);

    auto g_0_xxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 17);

    auto g_0_xxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 18);

    auto g_0_xxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 19);

    auto g_0_xxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 20);

    auto g_0_xxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 21);

    auto g_0_xxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 22);

    auto g_0_xxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 23);

    auto g_0_xxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 24);

    auto g_0_xxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 25);

    auto g_0_xxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 26);

    auto g_0_xxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 27);

    auto g_0_xxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 28);

    auto g_0_xxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 29);

    auto g_0_xxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 30);

    auto g_0_xxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 31);

    auto g_0_xxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 32);

    auto g_0_xxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 33);

    auto g_0_xxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 34);

    auto g_0_xxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 35);

    auto g_0_xxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 36);

    auto g_0_xxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 37);

    auto g_0_xxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 38);

    auto g_0_xxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 39);

    auto g_0_xxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 40);

    auto g_0_xxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 41);

    auto g_0_xxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 42);

    auto g_0_xxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 43);

    auto g_0_xxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 44);

    auto g_0_xxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 45);

    auto g_0_xxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 46);

    auto g_0_xxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 47);

    auto g_0_xxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 48);

    auto g_0_xxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 50);

    auto g_0_xxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 51);

    auto g_0_xxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 54);

    auto g_0_xxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 55);

    auto g_0_xxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 59);

    auto g_0_xxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 60);

    auto g_0_xxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 65);

    auto g_0_xxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 66);

    auto g_0_xxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 72);

    auto g_0_xxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 73);

    auto g_0_xxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 80);

    auto g_0_xxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 81);

    auto g_0_xxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 90);

    auto g_0_xxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 91);

    auto g_0_xxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 92);

    auto g_0_xxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 93);

    auto g_0_xxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 94);

    auto g_0_xxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 95);

    auto g_0_xxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 96);

    auto g_0_xxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 97);

    auto g_0_xxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 98);

    auto g_0_xxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 99);

    auto g_0_xxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 100);

    auto g_0_xxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 101);

    auto g_0_xxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 102);

    auto g_0_xxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 103);

    auto g_0_xxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 104);

    auto g_0_xxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 105);

    auto g_0_xxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 106);

    auto g_0_xxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 107);

    auto g_0_xxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 108);

    auto g_0_xxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 109);

    auto g_0_xxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 110);

    auto g_0_xxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 111);

    auto g_0_xxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 112);

    auto g_0_xxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 113);

    auto g_0_xxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 114);

    auto g_0_xxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 115);

    auto g_0_xxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 116);

    auto g_0_xxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 117);

    auto g_0_xxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 118);

    auto g_0_xxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 119);

    auto g_0_xxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 120);

    auto g_0_xxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 121);

    auto g_0_xxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 122);

    auto g_0_xxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 123);

    auto g_0_xxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 124);

    auto g_0_xxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 125);

    auto g_0_xxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 127);

    auto g_0_xxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 128);

    auto g_0_xxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 129);

    auto g_0_xxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 130);

    auto g_0_xxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 131);

    auto g_0_xxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 132);

    auto g_0_xxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 133);

    auto g_0_xxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 134);

    auto g_0_xyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 135);

    auto g_0_xyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 136);

    auto g_0_xyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 138);

    auto g_0_xyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 139);

    auto g_0_xyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 141);

    auto g_0_xyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 142);

    auto g_0_xyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 143);

    auto g_0_xyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 145);

    auto g_0_xyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 146);

    auto g_0_xyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 147);

    auto g_0_xyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 148);

    auto g_0_xyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 150);

    auto g_0_xyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 151);

    auto g_0_xyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 152);

    auto g_0_xyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 153);

    auto g_0_xyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 154);

    auto g_0_xyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 156);

    auto g_0_xyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 157);

    auto g_0_xyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 158);

    auto g_0_xyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 159);

    auto g_0_xyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 160);

    auto g_0_xyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 161);

    auto g_0_xyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 163);

    auto g_0_xyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 164);

    auto g_0_xyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 165);

    auto g_0_xyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 166);

    auto g_0_xyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 167);

    auto g_0_xyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 168);

    auto g_0_xyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 169);

    auto g_0_xyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 171);

    auto g_0_xyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 172);

    auto g_0_xyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 173);

    auto g_0_xyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 174);

    auto g_0_xyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 175);

    auto g_0_xyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 176);

    auto g_0_xyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 177);

    auto g_0_xyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 178);

    auto g_0_xyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 179);

    auto g_0_xzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 225);

    auto g_0_xzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 227);

    auto g_0_xzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 229);

    auto g_0_xzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 230);

    auto g_0_xzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 232);

    auto g_0_xzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 233);

    auto g_0_xzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 234);

    auto g_0_xzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 236);

    auto g_0_xzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 237);

    auto g_0_xzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 238);

    auto g_0_xzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 239);

    auto g_0_xzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 241);

    auto g_0_xzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 242);

    auto g_0_xzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 243);

    auto g_0_xzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 244);

    auto g_0_xzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 245);

    auto g_0_xzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 247);

    auto g_0_xzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 248);

    auto g_0_xzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 249);

    auto g_0_xzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 250);

    auto g_0_xzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 251);

    auto g_0_xzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 252);

    auto g_0_xzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 254);

    auto g_0_xzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 255);

    auto g_0_xzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 256);

    auto g_0_xzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 257);

    auto g_0_xzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 258);

    auto g_0_xzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 259);

    auto g_0_xzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 260);

    auto g_0_xzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 261);

    auto g_0_xzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 262);

    auto g_0_xzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 263);

    auto g_0_xzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 264);

    auto g_0_xzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 265);

    auto g_0_xzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 266);

    auto g_0_xzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 267);

    auto g_0_xzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 268);

    auto g_0_xzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 269);

    auto g_0_yyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 270);

    auto g_0_yyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 271);

    auto g_0_yyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 272);

    auto g_0_yyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 273);

    auto g_0_yyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 274);

    auto g_0_yyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 275);

    auto g_0_yyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 276);

    auto g_0_yyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 277);

    auto g_0_yyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 278);

    auto g_0_yyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 279);

    auto g_0_yyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 280);

    auto g_0_yyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 281);

    auto g_0_yyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 282);

    auto g_0_yyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 283);

    auto g_0_yyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 284);

    auto g_0_yyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 285);

    auto g_0_yyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 286);

    auto g_0_yyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 287);

    auto g_0_yyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 288);

    auto g_0_yyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 289);

    auto g_0_yyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 290);

    auto g_0_yyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 291);

    auto g_0_yyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 292);

    auto g_0_yyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 293);

    auto g_0_yyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 294);

    auto g_0_yyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 295);

    auto g_0_yyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 296);

    auto g_0_yyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 297);

    auto g_0_yyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 298);

    auto g_0_yyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 299);

    auto g_0_yyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 300);

    auto g_0_yyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 301);

    auto g_0_yyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 302);

    auto g_0_yyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 303);

    auto g_0_yyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 304);

    auto g_0_yyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 305);

    auto g_0_yyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 306);

    auto g_0_yyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 307);

    auto g_0_yyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 308);

    auto g_0_yyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 309);

    auto g_0_yyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 310);

    auto g_0_yyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 311);

    auto g_0_yyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 312);

    auto g_0_yyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 313);

    auto g_0_yyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 314);

    auto g_0_yyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 316);

    auto g_0_yyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 317);

    auto g_0_yyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 318);

    auto g_0_yyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 319);

    auto g_0_yyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 320);

    auto g_0_yyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 321);

    auto g_0_yyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 322);

    auto g_0_yyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 323);

    auto g_0_yyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 324);

    auto g_0_yyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 325);

    auto g_0_yyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 326);

    auto g_0_yyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 327);

    auto g_0_yyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 328);

    auto g_0_yyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 329);

    auto g_0_yyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 330);

    auto g_0_yyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 331);

    auto g_0_yyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 332);

    auto g_0_yyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 333);

    auto g_0_yyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 334);

    auto g_0_yyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 335);

    auto g_0_yyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 336);

    auto g_0_yyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 337);

    auto g_0_yyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 338);

    auto g_0_yyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 339);

    auto g_0_yyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 340);

    auto g_0_yyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 341);

    auto g_0_yyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 342);

    auto g_0_yyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 343);

    auto g_0_yyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 344);

    auto g_0_yyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 345);

    auto g_0_yyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 346);

    auto g_0_yyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 347);

    auto g_0_yyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 348);

    auto g_0_yyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 349);

    auto g_0_yyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 350);

    auto g_0_yyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 351);

    auto g_0_yyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 352);

    auto g_0_yyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 353);

    auto g_0_yyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 354);

    auto g_0_yyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 355);

    auto g_0_yyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 356);

    auto g_0_yyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 357);

    auto g_0_yyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 358);

    auto g_0_yyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 359);

    auto g_0_yzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 360);

    auto g_0_yzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 361);

    auto g_0_yzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 362);

    auto g_0_yzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 363);

    auto g_0_yzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 364);

    auto g_0_yzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 365);

    auto g_0_yzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 366);

    auto g_0_yzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 367);

    auto g_0_yzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 368);

    auto g_0_yzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 369);

    auto g_0_yzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 370);

    auto g_0_yzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 371);

    auto g_0_yzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 372);

    auto g_0_yzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 373);

    auto g_0_yzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 374);

    auto g_0_yzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 375);

    auto g_0_yzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 376);

    auto g_0_yzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 377);

    auto g_0_yzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 378);

    auto g_0_yzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 379);

    auto g_0_yzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 380);

    auto g_0_yzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 381);

    auto g_0_yzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 382);

    auto g_0_yzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 383);

    auto g_0_yzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 384);

    auto g_0_yzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 385);

    auto g_0_yzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 386);

    auto g_0_yzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 387);

    auto g_0_yzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 388);

    auto g_0_yzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 389);

    auto g_0_yzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 390);

    auto g_0_yzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 391);

    auto g_0_yzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 392);

    auto g_0_yzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 393);

    auto g_0_yzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 394);

    auto g_0_yzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 395);

    auto g_0_yzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 396);

    auto g_0_yzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 397);

    auto g_0_yzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 398);

    auto g_0_yzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 399);

    auto g_0_yzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 400);

    auto g_0_yzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 401);

    auto g_0_yzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 402);

    auto g_0_yzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 403);

    auto g_0_yzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 404);

    auto g_0_zzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 405);

    auto g_0_zzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 406);

    auto g_0_zzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 407);

    auto g_0_zzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 408);

    auto g_0_zzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 409);

    auto g_0_zzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 410);

    auto g_0_zzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 411);

    auto g_0_zzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 412);

    auto g_0_zzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 413);

    auto g_0_zzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 414);

    auto g_0_zzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 415);

    auto g_0_zzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 416);

    auto g_0_zzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 417);

    auto g_0_zzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 418);

    auto g_0_zzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 419);

    auto g_0_zzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 420);

    auto g_0_zzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 421);

    auto g_0_zzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 422);

    auto g_0_zzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 423);

    auto g_0_zzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 424);

    auto g_0_zzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 425);

    auto g_0_zzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 426);

    auto g_0_zzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 427);

    auto g_0_zzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 428);

    auto g_0_zzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 429);

    auto g_0_zzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 430);

    auto g_0_zzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 431);

    auto g_0_zzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 432);

    auto g_0_zzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 433);

    auto g_0_zzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 434);

    auto g_0_zzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 435);

    auto g_0_zzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 436);

    auto g_0_zzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 437);

    auto g_0_zzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 438);

    auto g_0_zzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 439);

    auto g_0_zzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 440);

    auto g_0_zzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 441);

    auto g_0_zzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 442);

    auto g_0_zzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 443);

    auto g_0_zzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 444);

    auto g_0_zzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 445);

    auto g_0_zzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 446);

    auto g_0_zzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 447);

    auto g_0_zzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 448);

    auto g_0_zzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 449);

    /// Set up components of auxilary buffer : SFSL

    auto g_0_xxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl);

    auto g_0_xxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 1);

    auto g_0_xxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 2);

    auto g_0_xxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 3);

    auto g_0_xxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 4);

    auto g_0_xxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 5);

    auto g_0_xxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 6);

    auto g_0_xxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 7);

    auto g_0_xxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 8);

    auto g_0_xxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 9);

    auto g_0_xxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 10);

    auto g_0_xxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 11);

    auto g_0_xxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 12);

    auto g_0_xxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 13);

    auto g_0_xxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 14);

    auto g_0_xxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 15);

    auto g_0_xxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 16);

    auto g_0_xxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 17);

    auto g_0_xxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 18);

    auto g_0_xxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 19);

    auto g_0_xxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 20);

    auto g_0_xxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 21);

    auto g_0_xxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 22);

    auto g_0_xxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 23);

    auto g_0_xxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 24);

    auto g_0_xxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 25);

    auto g_0_xxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 26);

    auto g_0_xxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 27);

    auto g_0_xxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 28);

    auto g_0_xxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 29);

    auto g_0_xxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 30);

    auto g_0_xxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 31);

    auto g_0_xxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 32);

    auto g_0_xxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 33);

    auto g_0_xxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 34);

    auto g_0_xxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 35);

    auto g_0_xxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 36);

    auto g_0_xxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 37);

    auto g_0_xxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 38);

    auto g_0_xxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 39);

    auto g_0_xxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 40);

    auto g_0_xxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 41);

    auto g_0_xxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 42);

    auto g_0_xxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 43);

    auto g_0_xxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 44);

    auto g_0_xxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl + 45);

    auto g_0_xxy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 46);

    auto g_0_xxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 47);

    auto g_0_xxy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 48);

    auto g_0_xxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 50);

    auto g_0_xxy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 51);

    auto g_0_xxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 54);

    auto g_0_xxy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 55);

    auto g_0_xxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 59);

    auto g_0_xxy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 60);

    auto g_0_xxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 65);

    auto g_0_xxy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 66);

    auto g_0_xxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 72);

    auto g_0_xxy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 73);

    auto g_0_xxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 80);

    auto g_0_xxy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 81);

    auto g_0_xxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl + 90);

    auto g_0_xxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 91);

    auto g_0_xxz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 92);

    auto g_0_xxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 93);

    auto g_0_xxz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 94);

    auto g_0_xxz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 95);

    auto g_0_xxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 96);

    auto g_0_xxz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 97);

    auto g_0_xxz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 98);

    auto g_0_xxz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 99);

    auto g_0_xxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 100);

    auto g_0_xxz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 101);

    auto g_0_xxz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 102);

    auto g_0_xxz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 103);

    auto g_0_xxz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 104);

    auto g_0_xxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 105);

    auto g_0_xxz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 106);

    auto g_0_xxz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 107);

    auto g_0_xxz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 108);

    auto g_0_xxz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 109);

    auto g_0_xxz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 110);

    auto g_0_xxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 111);

    auto g_0_xxz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 112);

    auto g_0_xxz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 113);

    auto g_0_xxz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 114);

    auto g_0_xxz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 115);

    auto g_0_xxz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 116);

    auto g_0_xxz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 117);

    auto g_0_xxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 118);

    auto g_0_xxz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 119);

    auto g_0_xxz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 120);

    auto g_0_xxz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 121);

    auto g_0_xxz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 122);

    auto g_0_xxz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 123);

    auto g_0_xxz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 124);

    auto g_0_xxz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 125);

    auto g_0_xxz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 127);

    auto g_0_xxz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 128);

    auto g_0_xxz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 129);

    auto g_0_xxz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 130);

    auto g_0_xxz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 131);

    auto g_0_xxz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 132);

    auto g_0_xxz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 133);

    auto g_0_xxz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 134);

    auto g_0_xyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl + 135);

    auto g_0_xyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 136);

    auto g_0_xyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 138);

    auto g_0_xyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 139);

    auto g_0_xyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 141);

    auto g_0_xyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 142);

    auto g_0_xyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 143);

    auto g_0_xyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 145);

    auto g_0_xyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 146);

    auto g_0_xyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 147);

    auto g_0_xyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 148);

    auto g_0_xyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 150);

    auto g_0_xyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 151);

    auto g_0_xyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 152);

    auto g_0_xyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 153);

    auto g_0_xyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 154);

    auto g_0_xyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 156);

    auto g_0_xyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 157);

    auto g_0_xyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 158);

    auto g_0_xyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 159);

    auto g_0_xyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 160);

    auto g_0_xyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 161);

    auto g_0_xyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 163);

    auto g_0_xyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 164);

    auto g_0_xyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 165);

    auto g_0_xyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 166);

    auto g_0_xyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 167);

    auto g_0_xyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 168);

    auto g_0_xyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 169);

    auto g_0_xyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 171);

    auto g_0_xyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 172);

    auto g_0_xyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 173);

    auto g_0_xyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 174);

    auto g_0_xyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 175);

    auto g_0_xyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 176);

    auto g_0_xyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 177);

    auto g_0_xyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 178);

    auto g_0_xyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 179);

    auto g_0_xzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl + 225);

    auto g_0_xzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 227);

    auto g_0_xzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 229);

    auto g_0_xzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 230);

    auto g_0_xzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 232);

    auto g_0_xzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 233);

    auto g_0_xzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 234);

    auto g_0_xzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 236);

    auto g_0_xzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 237);

    auto g_0_xzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 238);

    auto g_0_xzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 239);

    auto g_0_xzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 241);

    auto g_0_xzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 242);

    auto g_0_xzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 243);

    auto g_0_xzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 244);

    auto g_0_xzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 245);

    auto g_0_xzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 247);

    auto g_0_xzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 248);

    auto g_0_xzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 249);

    auto g_0_xzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 250);

    auto g_0_xzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 251);

    auto g_0_xzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 252);

    auto g_0_xzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 254);

    auto g_0_xzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 255);

    auto g_0_xzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 256);

    auto g_0_xzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 257);

    auto g_0_xzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 258);

    auto g_0_xzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 259);

    auto g_0_xzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 260);

    auto g_0_xzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 261);

    auto g_0_xzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 262);

    auto g_0_xzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 263);

    auto g_0_xzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 264);

    auto g_0_xzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 265);

    auto g_0_xzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 266);

    auto g_0_xzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 267);

    auto g_0_xzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 268);

    auto g_0_xzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 269);

    auto g_0_yyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl + 270);

    auto g_0_yyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 271);

    auto g_0_yyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 272);

    auto g_0_yyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 273);

    auto g_0_yyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 274);

    auto g_0_yyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 275);

    auto g_0_yyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 276);

    auto g_0_yyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 277);

    auto g_0_yyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 278);

    auto g_0_yyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 279);

    auto g_0_yyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 280);

    auto g_0_yyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 281);

    auto g_0_yyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 282);

    auto g_0_yyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 283);

    auto g_0_yyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 284);

    auto g_0_yyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 285);

    auto g_0_yyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 286);

    auto g_0_yyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 287);

    auto g_0_yyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 288);

    auto g_0_yyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 289);

    auto g_0_yyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 290);

    auto g_0_yyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 291);

    auto g_0_yyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 292);

    auto g_0_yyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 293);

    auto g_0_yyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 294);

    auto g_0_yyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 295);

    auto g_0_yyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 296);

    auto g_0_yyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 297);

    auto g_0_yyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 298);

    auto g_0_yyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 299);

    auto g_0_yyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 300);

    auto g_0_yyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 301);

    auto g_0_yyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 302);

    auto g_0_yyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 303);

    auto g_0_yyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 304);

    auto g_0_yyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 305);

    auto g_0_yyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 306);

    auto g_0_yyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 307);

    auto g_0_yyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 308);

    auto g_0_yyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 309);

    auto g_0_yyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 310);

    auto g_0_yyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 311);

    auto g_0_yyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 312);

    auto g_0_yyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 313);

    auto g_0_yyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 314);

    auto g_0_yyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 316);

    auto g_0_yyz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 317);

    auto g_0_yyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 318);

    auto g_0_yyz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 319);

    auto g_0_yyz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 320);

    auto g_0_yyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 321);

    auto g_0_yyz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 322);

    auto g_0_yyz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 323);

    auto g_0_yyz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 324);

    auto g_0_yyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 325);

    auto g_0_yyz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 326);

    auto g_0_yyz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 327);

    auto g_0_yyz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 328);

    auto g_0_yyz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 329);

    auto g_0_yyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 330);

    auto g_0_yyz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 331);

    auto g_0_yyz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 332);

    auto g_0_yyz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 333);

    auto g_0_yyz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 334);

    auto g_0_yyz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 335);

    auto g_0_yyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 336);

    auto g_0_yyz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 337);

    auto g_0_yyz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 338);

    auto g_0_yyz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 339);

    auto g_0_yyz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 340);

    auto g_0_yyz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 341);

    auto g_0_yyz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 342);

    auto g_0_yyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 343);

    auto g_0_yyz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 344);

    auto g_0_yyz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 345);

    auto g_0_yyz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 346);

    auto g_0_yyz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 347);

    auto g_0_yyz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 348);

    auto g_0_yyz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 349);

    auto g_0_yyz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 350);

    auto g_0_yyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 351);

    auto g_0_yyz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 352);

    auto g_0_yyz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 353);

    auto g_0_yyz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 354);

    auto g_0_yyz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 355);

    auto g_0_yyz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 356);

    auto g_0_yyz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 357);

    auto g_0_yyz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 358);

    auto g_0_yyz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 359);

    auto g_0_yzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl + 360);

    auto g_0_yzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 361);

    auto g_0_yzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 362);

    auto g_0_yzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 363);

    auto g_0_yzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 364);

    auto g_0_yzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 365);

    auto g_0_yzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 366);

    auto g_0_yzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 367);

    auto g_0_yzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 368);

    auto g_0_yzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 369);

    auto g_0_yzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 370);

    auto g_0_yzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 371);

    auto g_0_yzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 372);

    auto g_0_yzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 373);

    auto g_0_yzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 374);

    auto g_0_yzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 375);

    auto g_0_yzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 376);

    auto g_0_yzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 377);

    auto g_0_yzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 378);

    auto g_0_yzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 379);

    auto g_0_yzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 380);

    auto g_0_yzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 381);

    auto g_0_yzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 382);

    auto g_0_yzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 383);

    auto g_0_yzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 384);

    auto g_0_yzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 385);

    auto g_0_yzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 386);

    auto g_0_yzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 387);

    auto g_0_yzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 388);

    auto g_0_yzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 389);

    auto g_0_yzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 390);

    auto g_0_yzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 391);

    auto g_0_yzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 392);

    auto g_0_yzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 393);

    auto g_0_yzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 394);

    auto g_0_yzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 395);

    auto g_0_yzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 396);

    auto g_0_yzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 397);

    auto g_0_yzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 398);

    auto g_0_yzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 399);

    auto g_0_yzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 400);

    auto g_0_yzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 401);

    auto g_0_yzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 402);

    auto g_0_yzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 403);

    auto g_0_yzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 404);

    auto g_0_zzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sfsl + 405);

    auto g_0_zzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sfsl + 406);

    auto g_0_zzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sfsl + 407);

    auto g_0_zzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sfsl + 408);

    auto g_0_zzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sfsl + 409);

    auto g_0_zzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sfsl + 410);

    auto g_0_zzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sfsl + 411);

    auto g_0_zzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sfsl + 412);

    auto g_0_zzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sfsl + 413);

    auto g_0_zzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sfsl + 414);

    auto g_0_zzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 415);

    auto g_0_zzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 416);

    auto g_0_zzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 417);

    auto g_0_zzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 418);

    auto g_0_zzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 419);

    auto g_0_zzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 420);

    auto g_0_zzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 421);

    auto g_0_zzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 422);

    auto g_0_zzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 423);

    auto g_0_zzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 424);

    auto g_0_zzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 425);

    auto g_0_zzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 426);

    auto g_0_zzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 427);

    auto g_0_zzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 428);

    auto g_0_zzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 429);

    auto g_0_zzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 430);

    auto g_0_zzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 431);

    auto g_0_zzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 432);

    auto g_0_zzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 433);

    auto g_0_zzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 434);

    auto g_0_zzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 435);

    auto g_0_zzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 436);

    auto g_0_zzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 437);

    auto g_0_zzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 438);

    auto g_0_zzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 439);

    auto g_0_zzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 440);

    auto g_0_zzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sfsl + 441);

    auto g_0_zzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sfsl + 442);

    auto g_0_zzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sfsl + 443);

    auto g_0_zzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sfsl + 444);

    auto g_0_zzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 445);

    auto g_0_zzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 446);

    auto g_0_zzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 447);

    auto g_0_zzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 448);

    auto g_0_zzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sfsl + 449);

    /// Set up 0-45 components of targeted buffer : SGSL

    auto g_0_xxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl);

    auto g_0_xxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 1);

    auto g_0_xxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 2);

    auto g_0_xxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 3);

    auto g_0_xxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 4);

    auto g_0_xxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 5);

    auto g_0_xxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 6);

    auto g_0_xxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 7);

    auto g_0_xxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 8);

    auto g_0_xxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 9);

    auto g_0_xxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 10);

    auto g_0_xxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 11);

    auto g_0_xxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 12);

    auto g_0_xxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 13);

    auto g_0_xxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 14);

    auto g_0_xxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 15);

    auto g_0_xxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 16);

    auto g_0_xxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 17);

    auto g_0_xxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 18);

    auto g_0_xxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 19);

    auto g_0_xxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 20);

    auto g_0_xxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 21);

    auto g_0_xxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 22);

    auto g_0_xxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 23);

    auto g_0_xxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 24);

    auto g_0_xxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 25);

    auto g_0_xxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 26);

    auto g_0_xxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 27);

    auto g_0_xxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 28);

    auto g_0_xxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 29);

    auto g_0_xxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 30);

    auto g_0_xxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 31);

    auto g_0_xxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 32);

    auto g_0_xxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 33);

    auto g_0_xxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 34);

    auto g_0_xxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 35);

    auto g_0_xxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 36);

    auto g_0_xxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 37);

    auto g_0_xxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 38);

    auto g_0_xxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 39);

    auto g_0_xxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 40);

    auto g_0_xxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 41);

    auto g_0_xxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 42);

    auto g_0_xxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 43);

    auto g_0_xxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 44);

#pragma omp simd aligned(g_0_xx_0_xxxxxxxx_0,       \
                             g_0_xx_0_xxxxxxxx_1,   \
                             g_0_xx_0_xxxxxxxy_0,   \
                             g_0_xx_0_xxxxxxxy_1,   \
                             g_0_xx_0_xxxxxxxz_0,   \
                             g_0_xx_0_xxxxxxxz_1,   \
                             g_0_xx_0_xxxxxxyy_0,   \
                             g_0_xx_0_xxxxxxyy_1,   \
                             g_0_xx_0_xxxxxxyz_0,   \
                             g_0_xx_0_xxxxxxyz_1,   \
                             g_0_xx_0_xxxxxxzz_0,   \
                             g_0_xx_0_xxxxxxzz_1,   \
                             g_0_xx_0_xxxxxyyy_0,   \
                             g_0_xx_0_xxxxxyyy_1,   \
                             g_0_xx_0_xxxxxyyz_0,   \
                             g_0_xx_0_xxxxxyyz_1,   \
                             g_0_xx_0_xxxxxyzz_0,   \
                             g_0_xx_0_xxxxxyzz_1,   \
                             g_0_xx_0_xxxxxzzz_0,   \
                             g_0_xx_0_xxxxxzzz_1,   \
                             g_0_xx_0_xxxxyyyy_0,   \
                             g_0_xx_0_xxxxyyyy_1,   \
                             g_0_xx_0_xxxxyyyz_0,   \
                             g_0_xx_0_xxxxyyyz_1,   \
                             g_0_xx_0_xxxxyyzz_0,   \
                             g_0_xx_0_xxxxyyzz_1,   \
                             g_0_xx_0_xxxxyzzz_0,   \
                             g_0_xx_0_xxxxyzzz_1,   \
                             g_0_xx_0_xxxxzzzz_0,   \
                             g_0_xx_0_xxxxzzzz_1,   \
                             g_0_xx_0_xxxyyyyy_0,   \
                             g_0_xx_0_xxxyyyyy_1,   \
                             g_0_xx_0_xxxyyyyz_0,   \
                             g_0_xx_0_xxxyyyyz_1,   \
                             g_0_xx_0_xxxyyyzz_0,   \
                             g_0_xx_0_xxxyyyzz_1,   \
                             g_0_xx_0_xxxyyzzz_0,   \
                             g_0_xx_0_xxxyyzzz_1,   \
                             g_0_xx_0_xxxyzzzz_0,   \
                             g_0_xx_0_xxxyzzzz_1,   \
                             g_0_xx_0_xxxzzzzz_0,   \
                             g_0_xx_0_xxxzzzzz_1,   \
                             g_0_xx_0_xxyyyyyy_0,   \
                             g_0_xx_0_xxyyyyyy_1,   \
                             g_0_xx_0_xxyyyyyz_0,   \
                             g_0_xx_0_xxyyyyyz_1,   \
                             g_0_xx_0_xxyyyyzz_0,   \
                             g_0_xx_0_xxyyyyzz_1,   \
                             g_0_xx_0_xxyyyzzz_0,   \
                             g_0_xx_0_xxyyyzzz_1,   \
                             g_0_xx_0_xxyyzzzz_0,   \
                             g_0_xx_0_xxyyzzzz_1,   \
                             g_0_xx_0_xxyzzzzz_0,   \
                             g_0_xx_0_xxyzzzzz_1,   \
                             g_0_xx_0_xxzzzzzz_0,   \
                             g_0_xx_0_xxzzzzzz_1,   \
                             g_0_xx_0_xyyyyyyy_0,   \
                             g_0_xx_0_xyyyyyyy_1,   \
                             g_0_xx_0_xyyyyyyz_0,   \
                             g_0_xx_0_xyyyyyyz_1,   \
                             g_0_xx_0_xyyyyyzz_0,   \
                             g_0_xx_0_xyyyyyzz_1,   \
                             g_0_xx_0_xyyyyzzz_0,   \
                             g_0_xx_0_xyyyyzzz_1,   \
                             g_0_xx_0_xyyyzzzz_0,   \
                             g_0_xx_0_xyyyzzzz_1,   \
                             g_0_xx_0_xyyzzzzz_0,   \
                             g_0_xx_0_xyyzzzzz_1,   \
                             g_0_xx_0_xyzzzzzz_0,   \
                             g_0_xx_0_xyzzzzzz_1,   \
                             g_0_xx_0_xzzzzzzz_0,   \
                             g_0_xx_0_xzzzzzzz_1,   \
                             g_0_xx_0_yyyyyyyy_0,   \
                             g_0_xx_0_yyyyyyyy_1,   \
                             g_0_xx_0_yyyyyyyz_0,   \
                             g_0_xx_0_yyyyyyyz_1,   \
                             g_0_xx_0_yyyyyyzz_0,   \
                             g_0_xx_0_yyyyyyzz_1,   \
                             g_0_xx_0_yyyyyzzz_0,   \
                             g_0_xx_0_yyyyyzzz_1,   \
                             g_0_xx_0_yyyyzzzz_0,   \
                             g_0_xx_0_yyyyzzzz_1,   \
                             g_0_xx_0_yyyzzzzz_0,   \
                             g_0_xx_0_yyyzzzzz_1,   \
                             g_0_xx_0_yyzzzzzz_0,   \
                             g_0_xx_0_yyzzzzzz_1,   \
                             g_0_xx_0_yzzzzzzz_0,   \
                             g_0_xx_0_yzzzzzzz_1,   \
                             g_0_xx_0_zzzzzzzz_0,   \
                             g_0_xx_0_zzzzzzzz_1,   \
                             g_0_xxx_0_xxxxxxx_1,   \
                             g_0_xxx_0_xxxxxxxx_0,  \
                             g_0_xxx_0_xxxxxxxx_1,  \
                             g_0_xxx_0_xxxxxxxy_0,  \
                             g_0_xxx_0_xxxxxxxy_1,  \
                             g_0_xxx_0_xxxxxxxz_0,  \
                             g_0_xxx_0_xxxxxxxz_1,  \
                             g_0_xxx_0_xxxxxxy_1,   \
                             g_0_xxx_0_xxxxxxyy_0,  \
                             g_0_xxx_0_xxxxxxyy_1,  \
                             g_0_xxx_0_xxxxxxyz_0,  \
                             g_0_xxx_0_xxxxxxyz_1,  \
                             g_0_xxx_0_xxxxxxz_1,   \
                             g_0_xxx_0_xxxxxxzz_0,  \
                             g_0_xxx_0_xxxxxxzz_1,  \
                             g_0_xxx_0_xxxxxyy_1,   \
                             g_0_xxx_0_xxxxxyyy_0,  \
                             g_0_xxx_0_xxxxxyyy_1,  \
                             g_0_xxx_0_xxxxxyyz_0,  \
                             g_0_xxx_0_xxxxxyyz_1,  \
                             g_0_xxx_0_xxxxxyz_1,   \
                             g_0_xxx_0_xxxxxyzz_0,  \
                             g_0_xxx_0_xxxxxyzz_1,  \
                             g_0_xxx_0_xxxxxzz_1,   \
                             g_0_xxx_0_xxxxxzzz_0,  \
                             g_0_xxx_0_xxxxxzzz_1,  \
                             g_0_xxx_0_xxxxyyy_1,   \
                             g_0_xxx_0_xxxxyyyy_0,  \
                             g_0_xxx_0_xxxxyyyy_1,  \
                             g_0_xxx_0_xxxxyyyz_0,  \
                             g_0_xxx_0_xxxxyyyz_1,  \
                             g_0_xxx_0_xxxxyyz_1,   \
                             g_0_xxx_0_xxxxyyzz_0,  \
                             g_0_xxx_0_xxxxyyzz_1,  \
                             g_0_xxx_0_xxxxyzz_1,   \
                             g_0_xxx_0_xxxxyzzz_0,  \
                             g_0_xxx_0_xxxxyzzz_1,  \
                             g_0_xxx_0_xxxxzzz_1,   \
                             g_0_xxx_0_xxxxzzzz_0,  \
                             g_0_xxx_0_xxxxzzzz_1,  \
                             g_0_xxx_0_xxxyyyy_1,   \
                             g_0_xxx_0_xxxyyyyy_0,  \
                             g_0_xxx_0_xxxyyyyy_1,  \
                             g_0_xxx_0_xxxyyyyz_0,  \
                             g_0_xxx_0_xxxyyyyz_1,  \
                             g_0_xxx_0_xxxyyyz_1,   \
                             g_0_xxx_0_xxxyyyzz_0,  \
                             g_0_xxx_0_xxxyyyzz_1,  \
                             g_0_xxx_0_xxxyyzz_1,   \
                             g_0_xxx_0_xxxyyzzz_0,  \
                             g_0_xxx_0_xxxyyzzz_1,  \
                             g_0_xxx_0_xxxyzzz_1,   \
                             g_0_xxx_0_xxxyzzzz_0,  \
                             g_0_xxx_0_xxxyzzzz_1,  \
                             g_0_xxx_0_xxxzzzz_1,   \
                             g_0_xxx_0_xxxzzzzz_0,  \
                             g_0_xxx_0_xxxzzzzz_1,  \
                             g_0_xxx_0_xxyyyyy_1,   \
                             g_0_xxx_0_xxyyyyyy_0,  \
                             g_0_xxx_0_xxyyyyyy_1,  \
                             g_0_xxx_0_xxyyyyyz_0,  \
                             g_0_xxx_0_xxyyyyyz_1,  \
                             g_0_xxx_0_xxyyyyz_1,   \
                             g_0_xxx_0_xxyyyyzz_0,  \
                             g_0_xxx_0_xxyyyyzz_1,  \
                             g_0_xxx_0_xxyyyzz_1,   \
                             g_0_xxx_0_xxyyyzzz_0,  \
                             g_0_xxx_0_xxyyyzzz_1,  \
                             g_0_xxx_0_xxyyzzz_1,   \
                             g_0_xxx_0_xxyyzzzz_0,  \
                             g_0_xxx_0_xxyyzzzz_1,  \
                             g_0_xxx_0_xxyzzzz_1,   \
                             g_0_xxx_0_xxyzzzzz_0,  \
                             g_0_xxx_0_xxyzzzzz_1,  \
                             g_0_xxx_0_xxzzzzz_1,   \
                             g_0_xxx_0_xxzzzzzz_0,  \
                             g_0_xxx_0_xxzzzzzz_1,  \
                             g_0_xxx_0_xyyyyyy_1,   \
                             g_0_xxx_0_xyyyyyyy_0,  \
                             g_0_xxx_0_xyyyyyyy_1,  \
                             g_0_xxx_0_xyyyyyyz_0,  \
                             g_0_xxx_0_xyyyyyyz_1,  \
                             g_0_xxx_0_xyyyyyz_1,   \
                             g_0_xxx_0_xyyyyyzz_0,  \
                             g_0_xxx_0_xyyyyyzz_1,  \
                             g_0_xxx_0_xyyyyzz_1,   \
                             g_0_xxx_0_xyyyyzzz_0,  \
                             g_0_xxx_0_xyyyyzzz_1,  \
                             g_0_xxx_0_xyyyzzz_1,   \
                             g_0_xxx_0_xyyyzzzz_0,  \
                             g_0_xxx_0_xyyyzzzz_1,  \
                             g_0_xxx_0_xyyzzzz_1,   \
                             g_0_xxx_0_xyyzzzzz_0,  \
                             g_0_xxx_0_xyyzzzzz_1,  \
                             g_0_xxx_0_xyzzzzz_1,   \
                             g_0_xxx_0_xyzzzzzz_0,  \
                             g_0_xxx_0_xyzzzzzz_1,  \
                             g_0_xxx_0_xzzzzzz_1,   \
                             g_0_xxx_0_xzzzzzzz_0,  \
                             g_0_xxx_0_xzzzzzzz_1,  \
                             g_0_xxx_0_yyyyyyy_1,   \
                             g_0_xxx_0_yyyyyyyy_0,  \
                             g_0_xxx_0_yyyyyyyy_1,  \
                             g_0_xxx_0_yyyyyyyz_0,  \
                             g_0_xxx_0_yyyyyyyz_1,  \
                             g_0_xxx_0_yyyyyyz_1,   \
                             g_0_xxx_0_yyyyyyzz_0,  \
                             g_0_xxx_0_yyyyyyzz_1,  \
                             g_0_xxx_0_yyyyyzz_1,   \
                             g_0_xxx_0_yyyyyzzz_0,  \
                             g_0_xxx_0_yyyyyzzz_1,  \
                             g_0_xxx_0_yyyyzzz_1,   \
                             g_0_xxx_0_yyyyzzzz_0,  \
                             g_0_xxx_0_yyyyzzzz_1,  \
                             g_0_xxx_0_yyyzzzz_1,   \
                             g_0_xxx_0_yyyzzzzz_0,  \
                             g_0_xxx_0_yyyzzzzz_1,  \
                             g_0_xxx_0_yyzzzzz_1,   \
                             g_0_xxx_0_yyzzzzzz_0,  \
                             g_0_xxx_0_yyzzzzzz_1,  \
                             g_0_xxx_0_yzzzzzz_1,   \
                             g_0_xxx_0_yzzzzzzz_0,  \
                             g_0_xxx_0_yzzzzzzz_1,  \
                             g_0_xxx_0_zzzzzzz_1,   \
                             g_0_xxx_0_zzzzzzzz_0,  \
                             g_0_xxx_0_zzzzzzzz_1,  \
                             g_0_xxxx_0_xxxxxxxx_0, \
                             g_0_xxxx_0_xxxxxxxy_0, \
                             g_0_xxxx_0_xxxxxxxz_0, \
                             g_0_xxxx_0_xxxxxxyy_0, \
                             g_0_xxxx_0_xxxxxxyz_0, \
                             g_0_xxxx_0_xxxxxxzz_0, \
                             g_0_xxxx_0_xxxxxyyy_0, \
                             g_0_xxxx_0_xxxxxyyz_0, \
                             g_0_xxxx_0_xxxxxyzz_0, \
                             g_0_xxxx_0_xxxxxzzz_0, \
                             g_0_xxxx_0_xxxxyyyy_0, \
                             g_0_xxxx_0_xxxxyyyz_0, \
                             g_0_xxxx_0_xxxxyyzz_0, \
                             g_0_xxxx_0_xxxxyzzz_0, \
                             g_0_xxxx_0_xxxxzzzz_0, \
                             g_0_xxxx_0_xxxyyyyy_0, \
                             g_0_xxxx_0_xxxyyyyz_0, \
                             g_0_xxxx_0_xxxyyyzz_0, \
                             g_0_xxxx_0_xxxyyzzz_0, \
                             g_0_xxxx_0_xxxyzzzz_0, \
                             g_0_xxxx_0_xxxzzzzz_0, \
                             g_0_xxxx_0_xxyyyyyy_0, \
                             g_0_xxxx_0_xxyyyyyz_0, \
                             g_0_xxxx_0_xxyyyyzz_0, \
                             g_0_xxxx_0_xxyyyzzz_0, \
                             g_0_xxxx_0_xxyyzzzz_0, \
                             g_0_xxxx_0_xxyzzzzz_0, \
                             g_0_xxxx_0_xxzzzzzz_0, \
                             g_0_xxxx_0_xyyyyyyy_0, \
                             g_0_xxxx_0_xyyyyyyz_0, \
                             g_0_xxxx_0_xyyyyyzz_0, \
                             g_0_xxxx_0_xyyyyzzz_0, \
                             g_0_xxxx_0_xyyyzzzz_0, \
                             g_0_xxxx_0_xyyzzzzz_0, \
                             g_0_xxxx_0_xyzzzzzz_0, \
                             g_0_xxxx_0_xzzzzzzz_0, \
                             g_0_xxxx_0_yyyyyyyy_0, \
                             g_0_xxxx_0_yyyyyyyz_0, \
                             g_0_xxxx_0_yyyyyyzz_0, \
                             g_0_xxxx_0_yyyyyzzz_0, \
                             g_0_xxxx_0_yyyyzzzz_0, \
                             g_0_xxxx_0_yyyzzzzz_0, \
                             g_0_xxxx_0_yyzzzzzz_0, \
                             g_0_xxxx_0_yzzzzzzz_0, \
                             g_0_xxxx_0_zzzzzzzz_0, \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxx_0_xxxxxxxx_0[i] = 3.0 * g_0_xx_0_xxxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxxx_1[i] * fti_ab_0 +
                                   8.0 * g_0_xxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxxx_0[i] * pb_x + g_0_xxx_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxxxy_0[i] = 3.0 * g_0_xx_0_xxxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxxy_1[i] * fti_ab_0 +
                                   7.0 * g_0_xxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxxy_0[i] * pb_x + g_0_xxx_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxxxz_0[i] = 3.0 * g_0_xx_0_xxxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxxz_1[i] * fti_ab_0 +
                                   7.0 * g_0_xxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxxz_0[i] * pb_x + g_0_xxx_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxxyy_0[i] = 3.0 * g_0_xx_0_xxxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxyy_1[i] * fti_ab_0 +
                                   6.0 * g_0_xxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxyy_0[i] * pb_x + g_0_xxx_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxxyz_0[i] = 3.0 * g_0_xx_0_xxxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxyz_1[i] * fti_ab_0 +
                                   6.0 * g_0_xxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxyz_0[i] * pb_x + g_0_xxx_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxxzz_0[i] = 3.0 * g_0_xx_0_xxxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxzz_1[i] * fti_ab_0 +
                                   6.0 * g_0_xxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxzz_0[i] * pb_x + g_0_xxx_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxyyy_0[i] = 3.0 * g_0_xx_0_xxxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxyyy_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyyy_0[i] * pb_x + g_0_xxx_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxyyz_0[i] = 3.0 * g_0_xx_0_xxxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxyyz_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyyz_0[i] * pb_x + g_0_xxx_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxyzz_0[i] = 3.0 * g_0_xx_0_xxxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxyzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyzz_0[i] * pb_x + g_0_xxx_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxzzz_0[i] = 3.0 * g_0_xx_0_xxxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxzzz_0[i] * pb_x + g_0_xxx_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyyyy_0[i] = 3.0 * g_0_xx_0_xxxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyyyy_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyyy_0[i] * pb_x + g_0_xxx_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyyyz_0[i] = 3.0 * g_0_xx_0_xxxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyyyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyyz_0[i] * pb_x + g_0_xxx_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyyzz_0[i] = 3.0 * g_0_xx_0_xxxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyyzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyzz_0[i] * pb_x + g_0_xxx_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyzzz_0[i] = 3.0 * g_0_xx_0_xxxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyzzz_0[i] * pb_x + g_0_xxx_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxzzzz_0[i] = 3.0 * g_0_xx_0_xxxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxzzzz_0[i] * pb_x + g_0_xxx_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyyyy_0[i] = 3.0 * g_0_xx_0_xxxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyyyy_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyyy_0[i] * pb_x + g_0_xxx_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyyyz_0[i] = 3.0 * g_0_xx_0_xxxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyyz_0[i] * pb_x + g_0_xxx_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyyzz_0[i] = 3.0 * g_0_xx_0_xxxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyzz_0[i] * pb_x + g_0_xxx_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyzzz_0[i] = 3.0 * g_0_xx_0_xxxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyzzz_0[i] * pb_x + g_0_xxx_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyzzzz_0[i] = 3.0 * g_0_xx_0_xxxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyzzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzzzz_0[i] * pb_x + g_0_xxx_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxzzzzz_0[i] = 3.0 * g_0_xx_0_xxxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxzzzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxzzzzz_0[i] * pb_x + g_0_xxx_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyyyy_0[i] = 3.0 * g_0_xx_0_xxyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyyyy_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyyy_0[i] * pb_x + g_0_xxx_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyyyz_0[i] = 3.0 * g_0_xx_0_xxyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyyz_0[i] * pb_x + g_0_xxx_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyyzz_0[i] = 3.0 * g_0_xx_0_xxyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyzz_0[i] * pb_x + g_0_xxx_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyzzz_0[i] = 3.0 * g_0_xx_0_xxyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyzzz_0[i] * pb_x + g_0_xxx_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyzzzz_0[i] = 3.0 * g_0_xx_0_xxyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzzzz_0[i] * pb_x + g_0_xxx_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyzzzzz_0[i] = 3.0 * g_0_xx_0_xxyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyzzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzzzz_0[i] * pb_x + g_0_xxx_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxzzzzzz_0[i] = 3.0 * g_0_xx_0_xxzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxzzzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzzzzzz_0[i] * pb_x + g_0_xxx_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyyyy_0[i] = 3.0 * g_0_xx_0_xyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyyyy_1[i] * fti_ab_0 +
                                   g_0_xxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyyy_0[i] * pb_x + g_0_xxx_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyyyz_0[i] = 3.0 * g_0_xx_0_xyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyyyz_1[i] * fti_ab_0 +
                                   g_0_xxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyyz_0[i] * pb_x + g_0_xxx_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyyzz_0[i] = 3.0 * g_0_xx_0_xyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyyzz_1[i] * fti_ab_0 +
                                   g_0_xxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyzz_0[i] * pb_x + g_0_xxx_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyzzz_0[i] = 3.0 * g_0_xx_0_xyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyzzz_1[i] * fti_ab_0 +
                                   g_0_xxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyzzz_0[i] * pb_x + g_0_xxx_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyzzzz_0[i] = 3.0 * g_0_xx_0_xyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyzzzz_1[i] * fti_ab_0 +
                                   g_0_xxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzzzz_0[i] * pb_x + g_0_xxx_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyzzzzz_0[i] = 3.0 * g_0_xx_0_xyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyzzzzz_1[i] * fti_ab_0 +
                                   g_0_xxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzzzz_0[i] * pb_x + g_0_xxx_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyzzzzzz_0[i] = 3.0 * g_0_xx_0_xyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyzzzzzz_1[i] * fti_ab_0 +
                                   g_0_xxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzzzz_0[i] * pb_x + g_0_xxx_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xzzzzzzz_0[i] = 3.0 * g_0_xx_0_xzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xzzzzzzz_1[i] * fti_ab_0 +
                                   g_0_xxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzzzzzz_0[i] * pb_x + g_0_xxx_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyyyy_0[i] = 3.0 * g_0_xx_0_yyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyyyy_0[i] * pb_x +
                                   g_0_xxx_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyyyz_0[i] = 3.0 * g_0_xx_0_yyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyyyz_0[i] * pb_x +
                                   g_0_xxx_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyyzz_0[i] = 3.0 * g_0_xx_0_yyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyyzz_0[i] * pb_x +
                                   g_0_xxx_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyzzz_0[i] = 3.0 * g_0_xx_0_yyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyzzz_0[i] * pb_x +
                                   g_0_xxx_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyzzzz_0[i] = 3.0 * g_0_xx_0_yyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyzzzz_0[i] * pb_x +
                                   g_0_xxx_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyzzzzz_0[i] = 3.0 * g_0_xx_0_yyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyzzzzz_0[i] * pb_x +
                                   g_0_xxx_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyzzzzzz_0[i] = 3.0 * g_0_xx_0_yyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzzzzzz_0[i] * pb_x +
                                   g_0_xxx_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yzzzzzzz_0[i] = 3.0 * g_0_xx_0_yzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzzzzzz_0[i] * pb_x +
                                   g_0_xxx_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_zzzzzzzz_0[i] = 3.0 * g_0_xx_0_zzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzzzzzz_0[i] * pb_x +
                                   g_0_xxx_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 45-90 components of targeted buffer : SGSL

    auto g_0_xxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 45);

    auto g_0_xxxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 46);

    auto g_0_xxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 47);

    auto g_0_xxxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 48);

    auto g_0_xxxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 49);

    auto g_0_xxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 50);

    auto g_0_xxxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 51);

    auto g_0_xxxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 52);

    auto g_0_xxxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 53);

    auto g_0_xxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 54);

    auto g_0_xxxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 55);

    auto g_0_xxxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 56);

    auto g_0_xxxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 57);

    auto g_0_xxxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 58);

    auto g_0_xxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 59);

    auto g_0_xxxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 60);

    auto g_0_xxxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 61);

    auto g_0_xxxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 62);

    auto g_0_xxxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 63);

    auto g_0_xxxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 64);

    auto g_0_xxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 65);

    auto g_0_xxxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 66);

    auto g_0_xxxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 67);

    auto g_0_xxxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 68);

    auto g_0_xxxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 69);

    auto g_0_xxxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 70);

    auto g_0_xxxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 71);

    auto g_0_xxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 72);

    auto g_0_xxxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 73);

    auto g_0_xxxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 74);

    auto g_0_xxxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 75);

    auto g_0_xxxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 76);

    auto g_0_xxxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 77);

    auto g_0_xxxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 78);

    auto g_0_xxxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 79);

    auto g_0_xxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 80);

    auto g_0_xxxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 81);

    auto g_0_xxxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 82);

    auto g_0_xxxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 83);

    auto g_0_xxxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 84);

    auto g_0_xxxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 85);

    auto g_0_xxxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 86);

    auto g_0_xxxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 87);

    auto g_0_xxxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 88);

    auto g_0_xxxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 89);

#pragma omp simd aligned(g_0_xxx_0_xxxxxxx_1,       \
                             g_0_xxx_0_xxxxxxxx_0,  \
                             g_0_xxx_0_xxxxxxxx_1,  \
                             g_0_xxx_0_xxxxxxxy_0,  \
                             g_0_xxx_0_xxxxxxxy_1,  \
                             g_0_xxx_0_xxxxxxxz_0,  \
                             g_0_xxx_0_xxxxxxxz_1,  \
                             g_0_xxx_0_xxxxxxy_1,   \
                             g_0_xxx_0_xxxxxxyy_0,  \
                             g_0_xxx_0_xxxxxxyy_1,  \
                             g_0_xxx_0_xxxxxxyz_0,  \
                             g_0_xxx_0_xxxxxxyz_1,  \
                             g_0_xxx_0_xxxxxxz_1,   \
                             g_0_xxx_0_xxxxxxzz_0,  \
                             g_0_xxx_0_xxxxxxzz_1,  \
                             g_0_xxx_0_xxxxxyy_1,   \
                             g_0_xxx_0_xxxxxyyy_0,  \
                             g_0_xxx_0_xxxxxyyy_1,  \
                             g_0_xxx_0_xxxxxyyz_0,  \
                             g_0_xxx_0_xxxxxyyz_1,  \
                             g_0_xxx_0_xxxxxyz_1,   \
                             g_0_xxx_0_xxxxxyzz_0,  \
                             g_0_xxx_0_xxxxxyzz_1,  \
                             g_0_xxx_0_xxxxxzz_1,   \
                             g_0_xxx_0_xxxxxzzz_0,  \
                             g_0_xxx_0_xxxxxzzz_1,  \
                             g_0_xxx_0_xxxxyyy_1,   \
                             g_0_xxx_0_xxxxyyyy_0,  \
                             g_0_xxx_0_xxxxyyyy_1,  \
                             g_0_xxx_0_xxxxyyyz_0,  \
                             g_0_xxx_0_xxxxyyyz_1,  \
                             g_0_xxx_0_xxxxyyz_1,   \
                             g_0_xxx_0_xxxxyyzz_0,  \
                             g_0_xxx_0_xxxxyyzz_1,  \
                             g_0_xxx_0_xxxxyzz_1,   \
                             g_0_xxx_0_xxxxyzzz_0,  \
                             g_0_xxx_0_xxxxyzzz_1,  \
                             g_0_xxx_0_xxxxzzz_1,   \
                             g_0_xxx_0_xxxxzzzz_0,  \
                             g_0_xxx_0_xxxxzzzz_1,  \
                             g_0_xxx_0_xxxyyyy_1,   \
                             g_0_xxx_0_xxxyyyyy_0,  \
                             g_0_xxx_0_xxxyyyyy_1,  \
                             g_0_xxx_0_xxxyyyyz_0,  \
                             g_0_xxx_0_xxxyyyyz_1,  \
                             g_0_xxx_0_xxxyyyz_1,   \
                             g_0_xxx_0_xxxyyyzz_0,  \
                             g_0_xxx_0_xxxyyyzz_1,  \
                             g_0_xxx_0_xxxyyzz_1,   \
                             g_0_xxx_0_xxxyyzzz_0,  \
                             g_0_xxx_0_xxxyyzzz_1,  \
                             g_0_xxx_0_xxxyzzz_1,   \
                             g_0_xxx_0_xxxyzzzz_0,  \
                             g_0_xxx_0_xxxyzzzz_1,  \
                             g_0_xxx_0_xxxzzzz_1,   \
                             g_0_xxx_0_xxxzzzzz_0,  \
                             g_0_xxx_0_xxxzzzzz_1,  \
                             g_0_xxx_0_xxyyyyy_1,   \
                             g_0_xxx_0_xxyyyyyy_0,  \
                             g_0_xxx_0_xxyyyyyy_1,  \
                             g_0_xxx_0_xxyyyyyz_0,  \
                             g_0_xxx_0_xxyyyyyz_1,  \
                             g_0_xxx_0_xxyyyyz_1,   \
                             g_0_xxx_0_xxyyyyzz_0,  \
                             g_0_xxx_0_xxyyyyzz_1,  \
                             g_0_xxx_0_xxyyyzz_1,   \
                             g_0_xxx_0_xxyyyzzz_0,  \
                             g_0_xxx_0_xxyyyzzz_1,  \
                             g_0_xxx_0_xxyyzzz_1,   \
                             g_0_xxx_0_xxyyzzzz_0,  \
                             g_0_xxx_0_xxyyzzzz_1,  \
                             g_0_xxx_0_xxyzzzz_1,   \
                             g_0_xxx_0_xxyzzzzz_0,  \
                             g_0_xxx_0_xxyzzzzz_1,  \
                             g_0_xxx_0_xxzzzzz_1,   \
                             g_0_xxx_0_xxzzzzzz_0,  \
                             g_0_xxx_0_xxzzzzzz_1,  \
                             g_0_xxx_0_xyyyyyy_1,   \
                             g_0_xxx_0_xyyyyyyy_0,  \
                             g_0_xxx_0_xyyyyyyy_1,  \
                             g_0_xxx_0_xyyyyyyz_0,  \
                             g_0_xxx_0_xyyyyyyz_1,  \
                             g_0_xxx_0_xyyyyyz_1,   \
                             g_0_xxx_0_xyyyyyzz_0,  \
                             g_0_xxx_0_xyyyyyzz_1,  \
                             g_0_xxx_0_xyyyyzz_1,   \
                             g_0_xxx_0_xyyyyzzz_0,  \
                             g_0_xxx_0_xyyyyzzz_1,  \
                             g_0_xxx_0_xyyyzzz_1,   \
                             g_0_xxx_0_xyyyzzzz_0,  \
                             g_0_xxx_0_xyyyzzzz_1,  \
                             g_0_xxx_0_xyyzzzz_1,   \
                             g_0_xxx_0_xyyzzzzz_0,  \
                             g_0_xxx_0_xyyzzzzz_1,  \
                             g_0_xxx_0_xyzzzzz_1,   \
                             g_0_xxx_0_xyzzzzzz_0,  \
                             g_0_xxx_0_xyzzzzzz_1,  \
                             g_0_xxx_0_xzzzzzz_1,   \
                             g_0_xxx_0_xzzzzzzz_0,  \
                             g_0_xxx_0_xzzzzzzz_1,  \
                             g_0_xxx_0_yyyyyyy_1,   \
                             g_0_xxx_0_yyyyyyyy_0,  \
                             g_0_xxx_0_yyyyyyyy_1,  \
                             g_0_xxx_0_yyyyyyyz_0,  \
                             g_0_xxx_0_yyyyyyyz_1,  \
                             g_0_xxx_0_yyyyyyz_1,   \
                             g_0_xxx_0_yyyyyyzz_0,  \
                             g_0_xxx_0_yyyyyyzz_1,  \
                             g_0_xxx_0_yyyyyzz_1,   \
                             g_0_xxx_0_yyyyyzzz_0,  \
                             g_0_xxx_0_yyyyyzzz_1,  \
                             g_0_xxx_0_yyyyzzz_1,   \
                             g_0_xxx_0_yyyyzzzz_0,  \
                             g_0_xxx_0_yyyyzzzz_1,  \
                             g_0_xxx_0_yyyzzzz_1,   \
                             g_0_xxx_0_yyyzzzzz_0,  \
                             g_0_xxx_0_yyyzzzzz_1,  \
                             g_0_xxx_0_yyzzzzz_1,   \
                             g_0_xxx_0_yyzzzzzz_0,  \
                             g_0_xxx_0_yyzzzzzz_1,  \
                             g_0_xxx_0_yzzzzzz_1,   \
                             g_0_xxx_0_yzzzzzzz_0,  \
                             g_0_xxx_0_yzzzzzzz_1,  \
                             g_0_xxx_0_zzzzzzz_1,   \
                             g_0_xxx_0_zzzzzzzz_0,  \
                             g_0_xxx_0_zzzzzzzz_1,  \
                             g_0_xxxy_0_xxxxxxxx_0, \
                             g_0_xxxy_0_xxxxxxxy_0, \
                             g_0_xxxy_0_xxxxxxxz_0, \
                             g_0_xxxy_0_xxxxxxyy_0, \
                             g_0_xxxy_0_xxxxxxyz_0, \
                             g_0_xxxy_0_xxxxxxzz_0, \
                             g_0_xxxy_0_xxxxxyyy_0, \
                             g_0_xxxy_0_xxxxxyyz_0, \
                             g_0_xxxy_0_xxxxxyzz_0, \
                             g_0_xxxy_0_xxxxxzzz_0, \
                             g_0_xxxy_0_xxxxyyyy_0, \
                             g_0_xxxy_0_xxxxyyyz_0, \
                             g_0_xxxy_0_xxxxyyzz_0, \
                             g_0_xxxy_0_xxxxyzzz_0, \
                             g_0_xxxy_0_xxxxzzzz_0, \
                             g_0_xxxy_0_xxxyyyyy_0, \
                             g_0_xxxy_0_xxxyyyyz_0, \
                             g_0_xxxy_0_xxxyyyzz_0, \
                             g_0_xxxy_0_xxxyyzzz_0, \
                             g_0_xxxy_0_xxxyzzzz_0, \
                             g_0_xxxy_0_xxxzzzzz_0, \
                             g_0_xxxy_0_xxyyyyyy_0, \
                             g_0_xxxy_0_xxyyyyyz_0, \
                             g_0_xxxy_0_xxyyyyzz_0, \
                             g_0_xxxy_0_xxyyyzzz_0, \
                             g_0_xxxy_0_xxyyzzzz_0, \
                             g_0_xxxy_0_xxyzzzzz_0, \
                             g_0_xxxy_0_xxzzzzzz_0, \
                             g_0_xxxy_0_xyyyyyyy_0, \
                             g_0_xxxy_0_xyyyyyyz_0, \
                             g_0_xxxy_0_xyyyyyzz_0, \
                             g_0_xxxy_0_xyyyyzzz_0, \
                             g_0_xxxy_0_xyyyzzzz_0, \
                             g_0_xxxy_0_xyyzzzzz_0, \
                             g_0_xxxy_0_xyzzzzzz_0, \
                             g_0_xxxy_0_xzzzzzzz_0, \
                             g_0_xxxy_0_yyyyyyyy_0, \
                             g_0_xxxy_0_yyyyyyyz_0, \
                             g_0_xxxy_0_yyyyyyzz_0, \
                             g_0_xxxy_0_yyyyyzzz_0, \
                             g_0_xxxy_0_yyyyzzzz_0, \
                             g_0_xxxy_0_yyyzzzzz_0, \
                             g_0_xxxy_0_yyzzzzzz_0, \
                             g_0_xxxy_0_yzzzzzzz_0, \
                             g_0_xxxy_0_zzzzzzzz_0, \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxy_0_xxxxxxxx_0[i] = g_0_xxx_0_xxxxxxxx_0[i] * pb_y + g_0_xxx_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxxxy_0[i] = g_0_xxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxxy_0[i] * pb_y + g_0_xxx_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxxxz_0[i] = g_0_xxx_0_xxxxxxxz_0[i] * pb_y + g_0_xxx_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxxyy_0[i] = 2.0 * g_0_xxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxyy_0[i] * pb_y + g_0_xxx_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxxyz_0[i] = g_0_xxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxyz_0[i] * pb_y + g_0_xxx_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxxzz_0[i] = g_0_xxx_0_xxxxxxzz_0[i] * pb_y + g_0_xxx_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxyyy_0[i] = 3.0 * g_0_xxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyyy_0[i] * pb_y + g_0_xxx_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxyyz_0[i] = 2.0 * g_0_xxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyyz_0[i] * pb_y + g_0_xxx_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxyzz_0[i] = g_0_xxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyzz_0[i] * pb_y + g_0_xxx_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxzzz_0[i] = g_0_xxx_0_xxxxxzzz_0[i] * pb_y + g_0_xxx_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyyyy_0[i] = 4.0 * g_0_xxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyyy_0[i] * pb_y + g_0_xxx_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyyyz_0[i] = 3.0 * g_0_xxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyyz_0[i] * pb_y + g_0_xxx_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyyzz_0[i] = 2.0 * g_0_xxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyzz_0[i] * pb_y + g_0_xxx_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyzzz_0[i] = g_0_xxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyzzz_0[i] * pb_y + g_0_xxx_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxzzzz_0[i] = g_0_xxx_0_xxxxzzzz_0[i] * pb_y + g_0_xxx_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyyyy_0[i] = 5.0 * g_0_xxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyyy_0[i] * pb_y + g_0_xxx_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyyyz_0[i] = 4.0 * g_0_xxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyyz_0[i] * pb_y + g_0_xxx_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyyzz_0[i] = 3.0 * g_0_xxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyzz_0[i] * pb_y + g_0_xxx_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyzzz_0[i] = 2.0 * g_0_xxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyzzz_0[i] * pb_y + g_0_xxx_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyzzzz_0[i] = g_0_xxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzzzz_0[i] * pb_y + g_0_xxx_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxzzzzz_0[i] = g_0_xxx_0_xxxzzzzz_0[i] * pb_y + g_0_xxx_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyyyy_0[i] = 6.0 * g_0_xxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyyy_0[i] * pb_y + g_0_xxx_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyyyz_0[i] = 5.0 * g_0_xxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyyz_0[i] * pb_y + g_0_xxx_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyyzz_0[i] = 4.0 * g_0_xxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyzz_0[i] * pb_y + g_0_xxx_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyzzz_0[i] = 3.0 * g_0_xxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyzzz_0[i] * pb_y + g_0_xxx_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyzzzz_0[i] = 2.0 * g_0_xxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzzzz_0[i] * pb_y + g_0_xxx_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyzzzzz_0[i] = g_0_xxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzzzz_0[i] * pb_y + g_0_xxx_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxzzzzzz_0[i] = g_0_xxx_0_xxzzzzzz_0[i] * pb_y + g_0_xxx_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyyyy_0[i] = 7.0 * g_0_xxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyyy_0[i] * pb_y + g_0_xxx_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyyyz_0[i] = 6.0 * g_0_xxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyyz_0[i] * pb_y + g_0_xxx_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyyzz_0[i] = 5.0 * g_0_xxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyzz_0[i] * pb_y + g_0_xxx_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyzzz_0[i] = 4.0 * g_0_xxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyzzz_0[i] * pb_y + g_0_xxx_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyzzzz_0[i] = 3.0 * g_0_xxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzzzz_0[i] * pb_y + g_0_xxx_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyzzzzz_0[i] = 2.0 * g_0_xxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzzzz_0[i] * pb_y + g_0_xxx_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyzzzzzz_0[i] = g_0_xxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzzzz_0[i] * pb_y + g_0_xxx_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xzzzzzzz_0[i] = g_0_xxx_0_xzzzzzzz_0[i] * pb_y + g_0_xxx_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyyyy_0[i] = 8.0 * g_0_xxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyyy_0[i] * pb_y + g_0_xxx_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyyyz_0[i] = 7.0 * g_0_xxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyyz_0[i] * pb_y + g_0_xxx_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyyzz_0[i] = 6.0 * g_0_xxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyzz_0[i] * pb_y + g_0_xxx_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyzzz_0[i] = 5.0 * g_0_xxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyzzz_0[i] * pb_y + g_0_xxx_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyzzzz_0[i] = 4.0 * g_0_xxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyzzzz_0[i] * pb_y + g_0_xxx_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyzzzzz_0[i] = 3.0 * g_0_xxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzzzzz_0[i] * pb_y + g_0_xxx_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyzzzzzz_0[i] = 2.0 * g_0_xxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzzzzz_0[i] * pb_y + g_0_xxx_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yzzzzzzz_0[i] = g_0_xxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzzzzz_0[i] * pb_y + g_0_xxx_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_zzzzzzzz_0[i] = g_0_xxx_0_zzzzzzzz_0[i] * pb_y + g_0_xxx_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 90-135 components of targeted buffer : SGSL

    auto g_0_xxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 90);

    auto g_0_xxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 91);

    auto g_0_xxxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 92);

    auto g_0_xxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 93);

    auto g_0_xxxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 94);

    auto g_0_xxxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 95);

    auto g_0_xxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 96);

    auto g_0_xxxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 97);

    auto g_0_xxxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 98);

    auto g_0_xxxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 99);

    auto g_0_xxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 100);

    auto g_0_xxxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 101);

    auto g_0_xxxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 102);

    auto g_0_xxxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 103);

    auto g_0_xxxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 104);

    auto g_0_xxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 105);

    auto g_0_xxxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 106);

    auto g_0_xxxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 107);

    auto g_0_xxxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 108);

    auto g_0_xxxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 109);

    auto g_0_xxxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 110);

    auto g_0_xxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 111);

    auto g_0_xxxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 112);

    auto g_0_xxxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 113);

    auto g_0_xxxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 114);

    auto g_0_xxxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 115);

    auto g_0_xxxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 116);

    auto g_0_xxxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 117);

    auto g_0_xxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 118);

    auto g_0_xxxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 119);

    auto g_0_xxxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 120);

    auto g_0_xxxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 121);

    auto g_0_xxxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 122);

    auto g_0_xxxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 123);

    auto g_0_xxxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 124);

    auto g_0_xxxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 125);

    auto g_0_xxxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 126);

    auto g_0_xxxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 127);

    auto g_0_xxxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 128);

    auto g_0_xxxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 129);

    auto g_0_xxxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 130);

    auto g_0_xxxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 131);

    auto g_0_xxxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 132);

    auto g_0_xxxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 133);

    auto g_0_xxxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 134);

#pragma omp simd aligned(g_0_xxx_0_xxxxxxx_1,       \
                             g_0_xxx_0_xxxxxxxx_0,  \
                             g_0_xxx_0_xxxxxxxx_1,  \
                             g_0_xxx_0_xxxxxxxy_0,  \
                             g_0_xxx_0_xxxxxxxy_1,  \
                             g_0_xxx_0_xxxxxxxz_0,  \
                             g_0_xxx_0_xxxxxxxz_1,  \
                             g_0_xxx_0_xxxxxxy_1,   \
                             g_0_xxx_0_xxxxxxyy_0,  \
                             g_0_xxx_0_xxxxxxyy_1,  \
                             g_0_xxx_0_xxxxxxyz_0,  \
                             g_0_xxx_0_xxxxxxyz_1,  \
                             g_0_xxx_0_xxxxxxz_1,   \
                             g_0_xxx_0_xxxxxxzz_0,  \
                             g_0_xxx_0_xxxxxxzz_1,  \
                             g_0_xxx_0_xxxxxyy_1,   \
                             g_0_xxx_0_xxxxxyyy_0,  \
                             g_0_xxx_0_xxxxxyyy_1,  \
                             g_0_xxx_0_xxxxxyyz_0,  \
                             g_0_xxx_0_xxxxxyyz_1,  \
                             g_0_xxx_0_xxxxxyz_1,   \
                             g_0_xxx_0_xxxxxyzz_0,  \
                             g_0_xxx_0_xxxxxyzz_1,  \
                             g_0_xxx_0_xxxxxzz_1,   \
                             g_0_xxx_0_xxxxxzzz_0,  \
                             g_0_xxx_0_xxxxxzzz_1,  \
                             g_0_xxx_0_xxxxyyy_1,   \
                             g_0_xxx_0_xxxxyyyy_0,  \
                             g_0_xxx_0_xxxxyyyy_1,  \
                             g_0_xxx_0_xxxxyyyz_0,  \
                             g_0_xxx_0_xxxxyyyz_1,  \
                             g_0_xxx_0_xxxxyyz_1,   \
                             g_0_xxx_0_xxxxyyzz_0,  \
                             g_0_xxx_0_xxxxyyzz_1,  \
                             g_0_xxx_0_xxxxyzz_1,   \
                             g_0_xxx_0_xxxxyzzz_0,  \
                             g_0_xxx_0_xxxxyzzz_1,  \
                             g_0_xxx_0_xxxxzzz_1,   \
                             g_0_xxx_0_xxxxzzzz_0,  \
                             g_0_xxx_0_xxxxzzzz_1,  \
                             g_0_xxx_0_xxxyyyy_1,   \
                             g_0_xxx_0_xxxyyyyy_0,  \
                             g_0_xxx_0_xxxyyyyy_1,  \
                             g_0_xxx_0_xxxyyyyz_0,  \
                             g_0_xxx_0_xxxyyyyz_1,  \
                             g_0_xxx_0_xxxyyyz_1,   \
                             g_0_xxx_0_xxxyyyzz_0,  \
                             g_0_xxx_0_xxxyyyzz_1,  \
                             g_0_xxx_0_xxxyyzz_1,   \
                             g_0_xxx_0_xxxyyzzz_0,  \
                             g_0_xxx_0_xxxyyzzz_1,  \
                             g_0_xxx_0_xxxyzzz_1,   \
                             g_0_xxx_0_xxxyzzzz_0,  \
                             g_0_xxx_0_xxxyzzzz_1,  \
                             g_0_xxx_0_xxxzzzz_1,   \
                             g_0_xxx_0_xxxzzzzz_0,  \
                             g_0_xxx_0_xxxzzzzz_1,  \
                             g_0_xxx_0_xxyyyyy_1,   \
                             g_0_xxx_0_xxyyyyyy_0,  \
                             g_0_xxx_0_xxyyyyyy_1,  \
                             g_0_xxx_0_xxyyyyyz_0,  \
                             g_0_xxx_0_xxyyyyyz_1,  \
                             g_0_xxx_0_xxyyyyz_1,   \
                             g_0_xxx_0_xxyyyyzz_0,  \
                             g_0_xxx_0_xxyyyyzz_1,  \
                             g_0_xxx_0_xxyyyzz_1,   \
                             g_0_xxx_0_xxyyyzzz_0,  \
                             g_0_xxx_0_xxyyyzzz_1,  \
                             g_0_xxx_0_xxyyzzz_1,   \
                             g_0_xxx_0_xxyyzzzz_0,  \
                             g_0_xxx_0_xxyyzzzz_1,  \
                             g_0_xxx_0_xxyzzzz_1,   \
                             g_0_xxx_0_xxyzzzzz_0,  \
                             g_0_xxx_0_xxyzzzzz_1,  \
                             g_0_xxx_0_xxzzzzz_1,   \
                             g_0_xxx_0_xxzzzzzz_0,  \
                             g_0_xxx_0_xxzzzzzz_1,  \
                             g_0_xxx_0_xyyyyyy_1,   \
                             g_0_xxx_0_xyyyyyyy_0,  \
                             g_0_xxx_0_xyyyyyyy_1,  \
                             g_0_xxx_0_xyyyyyyz_0,  \
                             g_0_xxx_0_xyyyyyyz_1,  \
                             g_0_xxx_0_xyyyyyz_1,   \
                             g_0_xxx_0_xyyyyyzz_0,  \
                             g_0_xxx_0_xyyyyyzz_1,  \
                             g_0_xxx_0_xyyyyzz_1,   \
                             g_0_xxx_0_xyyyyzzz_0,  \
                             g_0_xxx_0_xyyyyzzz_1,  \
                             g_0_xxx_0_xyyyzzz_1,   \
                             g_0_xxx_0_xyyyzzzz_0,  \
                             g_0_xxx_0_xyyyzzzz_1,  \
                             g_0_xxx_0_xyyzzzz_1,   \
                             g_0_xxx_0_xyyzzzzz_0,  \
                             g_0_xxx_0_xyyzzzzz_1,  \
                             g_0_xxx_0_xyzzzzz_1,   \
                             g_0_xxx_0_xyzzzzzz_0,  \
                             g_0_xxx_0_xyzzzzzz_1,  \
                             g_0_xxx_0_xzzzzzz_1,   \
                             g_0_xxx_0_xzzzzzzz_0,  \
                             g_0_xxx_0_xzzzzzzz_1,  \
                             g_0_xxx_0_yyyyyyy_1,   \
                             g_0_xxx_0_yyyyyyyy_0,  \
                             g_0_xxx_0_yyyyyyyy_1,  \
                             g_0_xxx_0_yyyyyyyz_0,  \
                             g_0_xxx_0_yyyyyyyz_1,  \
                             g_0_xxx_0_yyyyyyz_1,   \
                             g_0_xxx_0_yyyyyyzz_0,  \
                             g_0_xxx_0_yyyyyyzz_1,  \
                             g_0_xxx_0_yyyyyzz_1,   \
                             g_0_xxx_0_yyyyyzzz_0,  \
                             g_0_xxx_0_yyyyyzzz_1,  \
                             g_0_xxx_0_yyyyzzz_1,   \
                             g_0_xxx_0_yyyyzzzz_0,  \
                             g_0_xxx_0_yyyyzzzz_1,  \
                             g_0_xxx_0_yyyzzzz_1,   \
                             g_0_xxx_0_yyyzzzzz_0,  \
                             g_0_xxx_0_yyyzzzzz_1,  \
                             g_0_xxx_0_yyzzzzz_1,   \
                             g_0_xxx_0_yyzzzzzz_0,  \
                             g_0_xxx_0_yyzzzzzz_1,  \
                             g_0_xxx_0_yzzzzzz_1,   \
                             g_0_xxx_0_yzzzzzzz_0,  \
                             g_0_xxx_0_yzzzzzzz_1,  \
                             g_0_xxx_0_zzzzzzz_1,   \
                             g_0_xxx_0_zzzzzzzz_0,  \
                             g_0_xxx_0_zzzzzzzz_1,  \
                             g_0_xxxz_0_xxxxxxxx_0, \
                             g_0_xxxz_0_xxxxxxxy_0, \
                             g_0_xxxz_0_xxxxxxxz_0, \
                             g_0_xxxz_0_xxxxxxyy_0, \
                             g_0_xxxz_0_xxxxxxyz_0, \
                             g_0_xxxz_0_xxxxxxzz_0, \
                             g_0_xxxz_0_xxxxxyyy_0, \
                             g_0_xxxz_0_xxxxxyyz_0, \
                             g_0_xxxz_0_xxxxxyzz_0, \
                             g_0_xxxz_0_xxxxxzzz_0, \
                             g_0_xxxz_0_xxxxyyyy_0, \
                             g_0_xxxz_0_xxxxyyyz_0, \
                             g_0_xxxz_0_xxxxyyzz_0, \
                             g_0_xxxz_0_xxxxyzzz_0, \
                             g_0_xxxz_0_xxxxzzzz_0, \
                             g_0_xxxz_0_xxxyyyyy_0, \
                             g_0_xxxz_0_xxxyyyyz_0, \
                             g_0_xxxz_0_xxxyyyzz_0, \
                             g_0_xxxz_0_xxxyyzzz_0, \
                             g_0_xxxz_0_xxxyzzzz_0, \
                             g_0_xxxz_0_xxxzzzzz_0, \
                             g_0_xxxz_0_xxyyyyyy_0, \
                             g_0_xxxz_0_xxyyyyyz_0, \
                             g_0_xxxz_0_xxyyyyzz_0, \
                             g_0_xxxz_0_xxyyyzzz_0, \
                             g_0_xxxz_0_xxyyzzzz_0, \
                             g_0_xxxz_0_xxyzzzzz_0, \
                             g_0_xxxz_0_xxzzzzzz_0, \
                             g_0_xxxz_0_xyyyyyyy_0, \
                             g_0_xxxz_0_xyyyyyyz_0, \
                             g_0_xxxz_0_xyyyyyzz_0, \
                             g_0_xxxz_0_xyyyyzzz_0, \
                             g_0_xxxz_0_xyyyzzzz_0, \
                             g_0_xxxz_0_xyyzzzzz_0, \
                             g_0_xxxz_0_xyzzzzzz_0, \
                             g_0_xxxz_0_xzzzzzzz_0, \
                             g_0_xxxz_0_yyyyyyyy_0, \
                             g_0_xxxz_0_yyyyyyyz_0, \
                             g_0_xxxz_0_yyyyyyzz_0, \
                             g_0_xxxz_0_yyyyyzzz_0, \
                             g_0_xxxz_0_yyyyzzzz_0, \
                             g_0_xxxz_0_yyyzzzzz_0, \
                             g_0_xxxz_0_yyzzzzzz_0, \
                             g_0_xxxz_0_yzzzzzzz_0, \
                             g_0_xxxz_0_zzzzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxz_0_xxxxxxxx_0[i] = g_0_xxx_0_xxxxxxxx_0[i] * pb_z + g_0_xxx_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxxxy_0[i] = g_0_xxx_0_xxxxxxxy_0[i] * pb_z + g_0_xxx_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxxxz_0[i] = g_0_xxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxxz_0[i] * pb_z + g_0_xxx_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxxyy_0[i] = g_0_xxx_0_xxxxxxyy_0[i] * pb_z + g_0_xxx_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxxyz_0[i] = g_0_xxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxyz_0[i] * pb_z + g_0_xxx_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxxzz_0[i] = 2.0 * g_0_xxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxzz_0[i] * pb_z + g_0_xxx_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxyyy_0[i] = g_0_xxx_0_xxxxxyyy_0[i] * pb_z + g_0_xxx_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxyyz_0[i] = g_0_xxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyyz_0[i] * pb_z + g_0_xxx_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxyzz_0[i] = 2.0 * g_0_xxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyzz_0[i] * pb_z + g_0_xxx_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxzzz_0[i] = 3.0 * g_0_xxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxzzz_0[i] * pb_z + g_0_xxx_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyyyy_0[i] = g_0_xxx_0_xxxxyyyy_0[i] * pb_z + g_0_xxx_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyyyz_0[i] = g_0_xxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyyz_0[i] * pb_z + g_0_xxx_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyyzz_0[i] = 2.0 * g_0_xxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyzz_0[i] * pb_z + g_0_xxx_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyzzz_0[i] = 3.0 * g_0_xxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyzzz_0[i] * pb_z + g_0_xxx_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxzzzz_0[i] = 4.0 * g_0_xxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxzzzz_0[i] * pb_z + g_0_xxx_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyyyy_0[i] = g_0_xxx_0_xxxyyyyy_0[i] * pb_z + g_0_xxx_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyyyz_0[i] = g_0_xxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyyz_0[i] * pb_z + g_0_xxx_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyyzz_0[i] = 2.0 * g_0_xxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyzz_0[i] * pb_z + g_0_xxx_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyzzz_0[i] = 3.0 * g_0_xxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyzzz_0[i] * pb_z + g_0_xxx_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyzzzz_0[i] = 4.0 * g_0_xxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzzzz_0[i] * pb_z + g_0_xxx_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxzzzzz_0[i] = 5.0 * g_0_xxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxzzzzz_0[i] * pb_z + g_0_xxx_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyyyy_0[i] = g_0_xxx_0_xxyyyyyy_0[i] * pb_z + g_0_xxx_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyyyz_0[i] = g_0_xxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyyz_0[i] * pb_z + g_0_xxx_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyyzz_0[i] = 2.0 * g_0_xxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyzz_0[i] * pb_z + g_0_xxx_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyzzz_0[i] = 3.0 * g_0_xxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyzzz_0[i] * pb_z + g_0_xxx_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyzzzz_0[i] = 4.0 * g_0_xxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzzzz_0[i] * pb_z + g_0_xxx_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyzzzzz_0[i] = 5.0 * g_0_xxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzzzz_0[i] * pb_z + g_0_xxx_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxzzzzzz_0[i] = 6.0 * g_0_xxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzzzzzz_0[i] * pb_z + g_0_xxx_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyyyy_0[i] = g_0_xxx_0_xyyyyyyy_0[i] * pb_z + g_0_xxx_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyyyz_0[i] = g_0_xxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyyz_0[i] * pb_z + g_0_xxx_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyyzz_0[i] = 2.0 * g_0_xxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyzz_0[i] * pb_z + g_0_xxx_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyzzz_0[i] = 3.0 * g_0_xxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyzzz_0[i] * pb_z + g_0_xxx_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyzzzz_0[i] = 4.0 * g_0_xxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzzzz_0[i] * pb_z + g_0_xxx_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyzzzzz_0[i] = 5.0 * g_0_xxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzzzz_0[i] * pb_z + g_0_xxx_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyzzzzzz_0[i] = 6.0 * g_0_xxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzzzz_0[i] * pb_z + g_0_xxx_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xzzzzzzz_0[i] = 7.0 * g_0_xxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzzzzzz_0[i] * pb_z + g_0_xxx_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyyyy_0[i] = g_0_xxx_0_yyyyyyyy_0[i] * pb_z + g_0_xxx_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyyyz_0[i] = g_0_xxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyyz_0[i] * pb_z + g_0_xxx_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyyzz_0[i] = 2.0 * g_0_xxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyzz_0[i] * pb_z + g_0_xxx_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyzzz_0[i] = 3.0 * g_0_xxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyzzz_0[i] * pb_z + g_0_xxx_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyzzzz_0[i] = 4.0 * g_0_xxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyzzzz_0[i] * pb_z + g_0_xxx_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyzzzzz_0[i] = 5.0 * g_0_xxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzzzzz_0[i] * pb_z + g_0_xxx_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyzzzzzz_0[i] = 6.0 * g_0_xxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzzzzz_0[i] * pb_z + g_0_xxx_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yzzzzzzz_0[i] = 7.0 * g_0_xxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzzzzz_0[i] * pb_z + g_0_xxx_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_zzzzzzzz_0[i] = 8.0 * g_0_xxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_zzzzzzzz_0[i] * pb_z + g_0_xxx_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 135-180 components of targeted buffer : SGSL

    auto g_0_xxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 135);

    auto g_0_xxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 136);

    auto g_0_xxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 137);

    auto g_0_xxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 138);

    auto g_0_xxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 139);

    auto g_0_xxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 140);

    auto g_0_xxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 141);

    auto g_0_xxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 142);

    auto g_0_xxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 143);

    auto g_0_xxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 144);

    auto g_0_xxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 145);

    auto g_0_xxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 146);

    auto g_0_xxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 147);

    auto g_0_xxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 148);

    auto g_0_xxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 149);

    auto g_0_xxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 150);

    auto g_0_xxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 151);

    auto g_0_xxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 152);

    auto g_0_xxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 153);

    auto g_0_xxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 154);

    auto g_0_xxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 155);

    auto g_0_xxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 156);

    auto g_0_xxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 157);

    auto g_0_xxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 158);

    auto g_0_xxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 159);

    auto g_0_xxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 160);

    auto g_0_xxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 161);

    auto g_0_xxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 162);

    auto g_0_xxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 163);

    auto g_0_xxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 164);

    auto g_0_xxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 165);

    auto g_0_xxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 166);

    auto g_0_xxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 167);

    auto g_0_xxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 168);

    auto g_0_xxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 169);

    auto g_0_xxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 170);

    auto g_0_xxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 171);

    auto g_0_xxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 172);

    auto g_0_xxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 173);

    auto g_0_xxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 174);

    auto g_0_xxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 175);

    auto g_0_xxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 176);

    auto g_0_xxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 177);

    auto g_0_xxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 178);

    auto g_0_xxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 179);

#pragma omp simd aligned(g_0_xx_0_xxxxxxxx_0,       \
                             g_0_xx_0_xxxxxxxx_1,   \
                             g_0_xx_0_xxxxxxxz_0,   \
                             g_0_xx_0_xxxxxxxz_1,   \
                             g_0_xx_0_xxxxxxzz_0,   \
                             g_0_xx_0_xxxxxxzz_1,   \
                             g_0_xx_0_xxxxxzzz_0,   \
                             g_0_xx_0_xxxxxzzz_1,   \
                             g_0_xx_0_xxxxzzzz_0,   \
                             g_0_xx_0_xxxxzzzz_1,   \
                             g_0_xx_0_xxxzzzzz_0,   \
                             g_0_xx_0_xxxzzzzz_1,   \
                             g_0_xx_0_xxzzzzzz_0,   \
                             g_0_xx_0_xxzzzzzz_1,   \
                             g_0_xx_0_xzzzzzzz_0,   \
                             g_0_xx_0_xzzzzzzz_1,   \
                             g_0_xxy_0_xxxxxxxx_0,  \
                             g_0_xxy_0_xxxxxxxx_1,  \
                             g_0_xxy_0_xxxxxxxz_0,  \
                             g_0_xxy_0_xxxxxxxz_1,  \
                             g_0_xxy_0_xxxxxxzz_0,  \
                             g_0_xxy_0_xxxxxxzz_1,  \
                             g_0_xxy_0_xxxxxzzz_0,  \
                             g_0_xxy_0_xxxxxzzz_1,  \
                             g_0_xxy_0_xxxxzzzz_0,  \
                             g_0_xxy_0_xxxxzzzz_1,  \
                             g_0_xxy_0_xxxzzzzz_0,  \
                             g_0_xxy_0_xxxzzzzz_1,  \
                             g_0_xxy_0_xxzzzzzz_0,  \
                             g_0_xxy_0_xxzzzzzz_1,  \
                             g_0_xxy_0_xzzzzzzz_0,  \
                             g_0_xxy_0_xzzzzzzz_1,  \
                             g_0_xxyy_0_xxxxxxxx_0, \
                             g_0_xxyy_0_xxxxxxxy_0, \
                             g_0_xxyy_0_xxxxxxxz_0, \
                             g_0_xxyy_0_xxxxxxyy_0, \
                             g_0_xxyy_0_xxxxxxyz_0, \
                             g_0_xxyy_0_xxxxxxzz_0, \
                             g_0_xxyy_0_xxxxxyyy_0, \
                             g_0_xxyy_0_xxxxxyyz_0, \
                             g_0_xxyy_0_xxxxxyzz_0, \
                             g_0_xxyy_0_xxxxxzzz_0, \
                             g_0_xxyy_0_xxxxyyyy_0, \
                             g_0_xxyy_0_xxxxyyyz_0, \
                             g_0_xxyy_0_xxxxyyzz_0, \
                             g_0_xxyy_0_xxxxyzzz_0, \
                             g_0_xxyy_0_xxxxzzzz_0, \
                             g_0_xxyy_0_xxxyyyyy_0, \
                             g_0_xxyy_0_xxxyyyyz_0, \
                             g_0_xxyy_0_xxxyyyzz_0, \
                             g_0_xxyy_0_xxxyyzzz_0, \
                             g_0_xxyy_0_xxxyzzzz_0, \
                             g_0_xxyy_0_xxxzzzzz_0, \
                             g_0_xxyy_0_xxyyyyyy_0, \
                             g_0_xxyy_0_xxyyyyyz_0, \
                             g_0_xxyy_0_xxyyyyzz_0, \
                             g_0_xxyy_0_xxyyyzzz_0, \
                             g_0_xxyy_0_xxyyzzzz_0, \
                             g_0_xxyy_0_xxyzzzzz_0, \
                             g_0_xxyy_0_xxzzzzzz_0, \
                             g_0_xxyy_0_xyyyyyyy_0, \
                             g_0_xxyy_0_xyyyyyyz_0, \
                             g_0_xxyy_0_xyyyyyzz_0, \
                             g_0_xxyy_0_xyyyyzzz_0, \
                             g_0_xxyy_0_xyyyzzzz_0, \
                             g_0_xxyy_0_xyyzzzzz_0, \
                             g_0_xxyy_0_xyzzzzzz_0, \
                             g_0_xxyy_0_xzzzzzzz_0, \
                             g_0_xxyy_0_yyyyyyyy_0, \
                             g_0_xxyy_0_yyyyyyyz_0, \
                             g_0_xxyy_0_yyyyyyzz_0, \
                             g_0_xxyy_0_yyyyyzzz_0, \
                             g_0_xxyy_0_yyyyzzzz_0, \
                             g_0_xxyy_0_yyyzzzzz_0, \
                             g_0_xxyy_0_yyzzzzzz_0, \
                             g_0_xxyy_0_yzzzzzzz_0, \
                             g_0_xxyy_0_zzzzzzzz_0, \
                             g_0_xyy_0_xxxxxxxy_0,  \
                             g_0_xyy_0_xxxxxxxy_1,  \
                             g_0_xyy_0_xxxxxxy_1,   \
                             g_0_xyy_0_xxxxxxyy_0,  \
                             g_0_xyy_0_xxxxxxyy_1,  \
                             g_0_xyy_0_xxxxxxyz_0,  \
                             g_0_xyy_0_xxxxxxyz_1,  \
                             g_0_xyy_0_xxxxxyy_1,   \
                             g_0_xyy_0_xxxxxyyy_0,  \
                             g_0_xyy_0_xxxxxyyy_1,  \
                             g_0_xyy_0_xxxxxyyz_0,  \
                             g_0_xyy_0_xxxxxyyz_1,  \
                             g_0_xyy_0_xxxxxyz_1,   \
                             g_0_xyy_0_xxxxxyzz_0,  \
                             g_0_xyy_0_xxxxxyzz_1,  \
                             g_0_xyy_0_xxxxyyy_1,   \
                             g_0_xyy_0_xxxxyyyy_0,  \
                             g_0_xyy_0_xxxxyyyy_1,  \
                             g_0_xyy_0_xxxxyyyz_0,  \
                             g_0_xyy_0_xxxxyyyz_1,  \
                             g_0_xyy_0_xxxxyyz_1,   \
                             g_0_xyy_0_xxxxyyzz_0,  \
                             g_0_xyy_0_xxxxyyzz_1,  \
                             g_0_xyy_0_xxxxyzz_1,   \
                             g_0_xyy_0_xxxxyzzz_0,  \
                             g_0_xyy_0_xxxxyzzz_1,  \
                             g_0_xyy_0_xxxyyyy_1,   \
                             g_0_xyy_0_xxxyyyyy_0,  \
                             g_0_xyy_0_xxxyyyyy_1,  \
                             g_0_xyy_0_xxxyyyyz_0,  \
                             g_0_xyy_0_xxxyyyyz_1,  \
                             g_0_xyy_0_xxxyyyz_1,   \
                             g_0_xyy_0_xxxyyyzz_0,  \
                             g_0_xyy_0_xxxyyyzz_1,  \
                             g_0_xyy_0_xxxyyzz_1,   \
                             g_0_xyy_0_xxxyyzzz_0,  \
                             g_0_xyy_0_xxxyyzzz_1,  \
                             g_0_xyy_0_xxxyzzz_1,   \
                             g_0_xyy_0_xxxyzzzz_0,  \
                             g_0_xyy_0_xxxyzzzz_1,  \
                             g_0_xyy_0_xxyyyyy_1,   \
                             g_0_xyy_0_xxyyyyyy_0,  \
                             g_0_xyy_0_xxyyyyyy_1,  \
                             g_0_xyy_0_xxyyyyyz_0,  \
                             g_0_xyy_0_xxyyyyyz_1,  \
                             g_0_xyy_0_xxyyyyz_1,   \
                             g_0_xyy_0_xxyyyyzz_0,  \
                             g_0_xyy_0_xxyyyyzz_1,  \
                             g_0_xyy_0_xxyyyzz_1,   \
                             g_0_xyy_0_xxyyyzzz_0,  \
                             g_0_xyy_0_xxyyyzzz_1,  \
                             g_0_xyy_0_xxyyzzz_1,   \
                             g_0_xyy_0_xxyyzzzz_0,  \
                             g_0_xyy_0_xxyyzzzz_1,  \
                             g_0_xyy_0_xxyzzzz_1,   \
                             g_0_xyy_0_xxyzzzzz_0,  \
                             g_0_xyy_0_xxyzzzzz_1,  \
                             g_0_xyy_0_xyyyyyy_1,   \
                             g_0_xyy_0_xyyyyyyy_0,  \
                             g_0_xyy_0_xyyyyyyy_1,  \
                             g_0_xyy_0_xyyyyyyz_0,  \
                             g_0_xyy_0_xyyyyyyz_1,  \
                             g_0_xyy_0_xyyyyyz_1,   \
                             g_0_xyy_0_xyyyyyzz_0,  \
                             g_0_xyy_0_xyyyyyzz_1,  \
                             g_0_xyy_0_xyyyyzz_1,   \
                             g_0_xyy_0_xyyyyzzz_0,  \
                             g_0_xyy_0_xyyyyzzz_1,  \
                             g_0_xyy_0_xyyyzzz_1,   \
                             g_0_xyy_0_xyyyzzzz_0,  \
                             g_0_xyy_0_xyyyzzzz_1,  \
                             g_0_xyy_0_xyyzzzz_1,   \
                             g_0_xyy_0_xyyzzzzz_0,  \
                             g_0_xyy_0_xyyzzzzz_1,  \
                             g_0_xyy_0_xyzzzzz_1,   \
                             g_0_xyy_0_xyzzzzzz_0,  \
                             g_0_xyy_0_xyzzzzzz_1,  \
                             g_0_xyy_0_yyyyyyy_1,   \
                             g_0_xyy_0_yyyyyyyy_0,  \
                             g_0_xyy_0_yyyyyyyy_1,  \
                             g_0_xyy_0_yyyyyyyz_0,  \
                             g_0_xyy_0_yyyyyyyz_1,  \
                             g_0_xyy_0_yyyyyyz_1,   \
                             g_0_xyy_0_yyyyyyzz_0,  \
                             g_0_xyy_0_yyyyyyzz_1,  \
                             g_0_xyy_0_yyyyyzz_1,   \
                             g_0_xyy_0_yyyyyzzz_0,  \
                             g_0_xyy_0_yyyyyzzz_1,  \
                             g_0_xyy_0_yyyyzzz_1,   \
                             g_0_xyy_0_yyyyzzzz_0,  \
                             g_0_xyy_0_yyyyzzzz_1,  \
                             g_0_xyy_0_yyyzzzz_1,   \
                             g_0_xyy_0_yyyzzzzz_0,  \
                             g_0_xyy_0_yyyzzzzz_1,  \
                             g_0_xyy_0_yyzzzzz_1,   \
                             g_0_xyy_0_yyzzzzzz_0,  \
                             g_0_xyy_0_yyzzzzzz_1,  \
                             g_0_xyy_0_yzzzzzz_1,   \
                             g_0_xyy_0_yzzzzzzz_0,  \
                             g_0_xyy_0_yzzzzzzz_1,  \
                             g_0_xyy_0_zzzzzzzz_0,  \
                             g_0_xyy_0_zzzzzzzz_1,  \
                             g_0_yy_0_xxxxxxxy_0,   \
                             g_0_yy_0_xxxxxxxy_1,   \
                             g_0_yy_0_xxxxxxyy_0,   \
                             g_0_yy_0_xxxxxxyy_1,   \
                             g_0_yy_0_xxxxxxyz_0,   \
                             g_0_yy_0_xxxxxxyz_1,   \
                             g_0_yy_0_xxxxxyyy_0,   \
                             g_0_yy_0_xxxxxyyy_1,   \
                             g_0_yy_0_xxxxxyyz_0,   \
                             g_0_yy_0_xxxxxyyz_1,   \
                             g_0_yy_0_xxxxxyzz_0,   \
                             g_0_yy_0_xxxxxyzz_1,   \
                             g_0_yy_0_xxxxyyyy_0,   \
                             g_0_yy_0_xxxxyyyy_1,   \
                             g_0_yy_0_xxxxyyyz_0,   \
                             g_0_yy_0_xxxxyyyz_1,   \
                             g_0_yy_0_xxxxyyzz_0,   \
                             g_0_yy_0_xxxxyyzz_1,   \
                             g_0_yy_0_xxxxyzzz_0,   \
                             g_0_yy_0_xxxxyzzz_1,   \
                             g_0_yy_0_xxxyyyyy_0,   \
                             g_0_yy_0_xxxyyyyy_1,   \
                             g_0_yy_0_xxxyyyyz_0,   \
                             g_0_yy_0_xxxyyyyz_1,   \
                             g_0_yy_0_xxxyyyzz_0,   \
                             g_0_yy_0_xxxyyyzz_1,   \
                             g_0_yy_0_xxxyyzzz_0,   \
                             g_0_yy_0_xxxyyzzz_1,   \
                             g_0_yy_0_xxxyzzzz_0,   \
                             g_0_yy_0_xxxyzzzz_1,   \
                             g_0_yy_0_xxyyyyyy_0,   \
                             g_0_yy_0_xxyyyyyy_1,   \
                             g_0_yy_0_xxyyyyyz_0,   \
                             g_0_yy_0_xxyyyyyz_1,   \
                             g_0_yy_0_xxyyyyzz_0,   \
                             g_0_yy_0_xxyyyyzz_1,   \
                             g_0_yy_0_xxyyyzzz_0,   \
                             g_0_yy_0_xxyyyzzz_1,   \
                             g_0_yy_0_xxyyzzzz_0,   \
                             g_0_yy_0_xxyyzzzz_1,   \
                             g_0_yy_0_xxyzzzzz_0,   \
                             g_0_yy_0_xxyzzzzz_1,   \
                             g_0_yy_0_xyyyyyyy_0,   \
                             g_0_yy_0_xyyyyyyy_1,   \
                             g_0_yy_0_xyyyyyyz_0,   \
                             g_0_yy_0_xyyyyyyz_1,   \
                             g_0_yy_0_xyyyyyzz_0,   \
                             g_0_yy_0_xyyyyyzz_1,   \
                             g_0_yy_0_xyyyyzzz_0,   \
                             g_0_yy_0_xyyyyzzz_1,   \
                             g_0_yy_0_xyyyzzzz_0,   \
                             g_0_yy_0_xyyyzzzz_1,   \
                             g_0_yy_0_xyyzzzzz_0,   \
                             g_0_yy_0_xyyzzzzz_1,   \
                             g_0_yy_0_xyzzzzzz_0,   \
                             g_0_yy_0_xyzzzzzz_1,   \
                             g_0_yy_0_yyyyyyyy_0,   \
                             g_0_yy_0_yyyyyyyy_1,   \
                             g_0_yy_0_yyyyyyyz_0,   \
                             g_0_yy_0_yyyyyyyz_1,   \
                             g_0_yy_0_yyyyyyzz_0,   \
                             g_0_yy_0_yyyyyyzz_1,   \
                             g_0_yy_0_yyyyyzzz_0,   \
                             g_0_yy_0_yyyyyzzz_1,   \
                             g_0_yy_0_yyyyzzzz_0,   \
                             g_0_yy_0_yyyyzzzz_1,   \
                             g_0_yy_0_yyyzzzzz_0,   \
                             g_0_yy_0_yyyzzzzz_1,   \
                             g_0_yy_0_yyzzzzzz_0,   \
                             g_0_yy_0_yyzzzzzz_1,   \
                             g_0_yy_0_yzzzzzzz_0,   \
                             g_0_yy_0_yzzzzzzz_1,   \
                             g_0_yy_0_zzzzzzzz_0,   \
                             g_0_yy_0_zzzzzzzz_1,   \
                             wp_x,                  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyy_0_xxxxxxxx_0[i] =
            g_0_xx_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxxxx_0[i] * pb_y + g_0_xxy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxxxxy_0[i] = g_0_yy_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxxxy_1[i] * fti_ab_0 + 7.0 * g_0_xyy_0_xxxxxxy_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxxxxy_0[i] * pb_x + g_0_xyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxxxz_0[i] =
            g_0_xx_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxxxz_0[i] * pb_y + g_0_xxy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxxxyy_0[i] = g_0_yy_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxxyy_1[i] * fti_ab_0 + 6.0 * g_0_xyy_0_xxxxxyy_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxxxyy_0[i] * pb_x + g_0_xyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxxyz_0[i] = g_0_yy_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_yy_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_xyy_0_xxxxxyz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxxxyz_0[i] * pb_x + g_0_xyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxxzz_0[i] =
            g_0_xx_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxxzz_0[i] * pb_y + g_0_xxy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxxyyy_0[i] = g_0_yy_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxyyy_1[i] * fti_ab_0 + 5.0 * g_0_xyy_0_xxxxyyy_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxxyyy_0[i] * pb_x + g_0_xyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxyyz_0[i] = g_0_yy_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_yy_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_xyy_0_xxxxyyz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxxyyz_0[i] * pb_x + g_0_xyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxyzz_0[i] = g_0_yy_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_yy_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_xyy_0_xxxxyzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxxyzz_0[i] * pb_x + g_0_xyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxzzz_0[i] =
            g_0_xx_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_xx_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxzzz_0[i] * pb_y + g_0_xxy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxyyyy_0[i] = g_0_yy_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyyyy_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxyyyy_0[i] * pb_x + g_0_xyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxyyyz_0[i] = g_0_yy_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_yy_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxyyyz_0[i] * pb_x + g_0_xyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxyyzz_0[i] = g_0_yy_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_yy_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxyyzz_0[i] * pb_x + g_0_xyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxyzzz_0[i] = g_0_yy_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_yy_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxxyzzz_0[i] * pb_x + g_0_xyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxzzzz_0[i] =
            g_0_xx_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_xx_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxzzzz_0[i] * pb_y + g_0_xxy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxyyyyy_0[i] = g_0_yy_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyyyy_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxyyyyy_0[i] * pb_x + g_0_xyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyyyyz_0[i] = g_0_yy_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_yy_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxyyyyz_0[i] * pb_x + g_0_xyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyyyzz_0[i] = g_0_yy_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_yy_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxyyyzz_0[i] * pb_x + g_0_xyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyyzzz_0[i] = g_0_yy_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_yy_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxyyzzz_0[i] * pb_x + g_0_xyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyzzzz_0[i] = g_0_yy_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_yy_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxxyzzzz_0[i] * pb_x + g_0_xyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxzzzzz_0[i] =
            g_0_xx_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_xx_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxzzzzz_0[i] * pb_y + g_0_xxy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxyyyyyy_0[i] = g_0_yy_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyyyy_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxyyyyyy_0[i] * pb_x + g_0_xyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyyyyz_0[i] = g_0_yy_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_yy_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxyyyyyz_0[i] * pb_x + g_0_xyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyyyzz_0[i] = g_0_yy_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_yy_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxyyyyzz_0[i] * pb_x + g_0_xyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyyzzz_0[i] = g_0_yy_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_yy_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxyyyzzz_0[i] * pb_x + g_0_xyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyzzzz_0[i] = g_0_yy_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_yy_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxyyzzzz_0[i] * pb_x + g_0_xyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyzzzzz_0[i] = g_0_yy_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_yy_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xxyzzzzz_0[i] * pb_x + g_0_xyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxzzzzzz_0[i] =
            g_0_xx_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_xx_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxzzzzzz_0[i] * pb_y + g_0_xxy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xyyyyyyy_0[i] = g_0_yy_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyyy_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xyyyyyyy_0[i] * pb_x + g_0_xyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyyyyz_0[i] = g_0_yy_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_yy_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xyyyyyyz_0[i] * pb_x + g_0_xyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyyyzz_0[i] = g_0_yy_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_yy_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xyyyyyzz_0[i] * pb_x + g_0_xyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyyzzz_0[i] = g_0_yy_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_yy_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xyyyyzzz_0[i] * pb_x + g_0_xyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyzzzz_0[i] = g_0_yy_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_yy_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xyyyzzzz_0[i] * pb_x + g_0_xyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyzzzzz_0[i] = g_0_yy_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_yy_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xyyzzzzz_0[i] * pb_x + g_0_xyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyzzzzzz_0[i] = g_0_yy_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_yy_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyy_0_xyzzzzzz_0[i] * pb_x + g_0_xyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xzzzzzzz_0[i] =
            g_0_xx_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_xx_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xzzzzzzz_0[i] * pb_y + g_0_xxy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_yyyyyyyy_0[i] =
            g_0_yy_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyyyy_0[i] * pb_x + g_0_xyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyyyyz_0[i] =
            g_0_yy_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_yy_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyyyz_0[i] * pb_x + g_0_xyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyyyzz_0[i] =
            g_0_yy_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_yy_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyyzz_0[i] * pb_x + g_0_xyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyyzzz_0[i] =
            g_0_yy_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_yy_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyzzz_0[i] * pb_x + g_0_xyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyzzzz_0[i] =
            g_0_yy_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_yy_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyzzzz_0[i] * pb_x + g_0_xyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyzzzzz_0[i] =
            g_0_yy_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_yy_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyzzzzz_0[i] * pb_x + g_0_xyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyzzzzzz_0[i] =
            g_0_yy_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_yy_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzzzzzz_0[i] * pb_x + g_0_xyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yzzzzzzz_0[i] =
            g_0_yy_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_yy_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzzzzzz_0[i] * pb_x + g_0_xyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_zzzzzzzz_0[i] =
            g_0_yy_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_yy_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_zzzzzzzz_0[i] * pb_x + g_0_xyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 180-225 components of targeted buffer : SGSL

    auto g_0_xxyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 180);

    auto g_0_xxyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 181);

    auto g_0_xxyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 182);

    auto g_0_xxyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 183);

    auto g_0_xxyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 184);

    auto g_0_xxyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 185);

    auto g_0_xxyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 186);

    auto g_0_xxyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 187);

    auto g_0_xxyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 188);

    auto g_0_xxyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 189);

    auto g_0_xxyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 190);

    auto g_0_xxyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 191);

    auto g_0_xxyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 192);

    auto g_0_xxyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 193);

    auto g_0_xxyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 194);

    auto g_0_xxyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 195);

    auto g_0_xxyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 196);

    auto g_0_xxyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 197);

    auto g_0_xxyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 198);

    auto g_0_xxyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 199);

    auto g_0_xxyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 200);

    auto g_0_xxyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 201);

    auto g_0_xxyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 202);

    auto g_0_xxyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 203);

    auto g_0_xxyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 204);

    auto g_0_xxyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 205);

    auto g_0_xxyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 206);

    auto g_0_xxyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 207);

    auto g_0_xxyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 208);

    auto g_0_xxyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 209);

    auto g_0_xxyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 210);

    auto g_0_xxyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 211);

    auto g_0_xxyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 212);

    auto g_0_xxyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 213);

    auto g_0_xxyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 214);

    auto g_0_xxyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 215);

    auto g_0_xxyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 216);

    auto g_0_xxyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 217);

    auto g_0_xxyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 218);

    auto g_0_xxyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 219);

    auto g_0_xxyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 220);

    auto g_0_xxyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 221);

    auto g_0_xxyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 222);

    auto g_0_xxyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 223);

    auto g_0_xxyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 224);

#pragma omp simd aligned(g_0_xxy_0_xxxxxxxy_0,      \
                             g_0_xxy_0_xxxxxxxy_1,  \
                             g_0_xxy_0_xxxxxxyy_0,  \
                             g_0_xxy_0_xxxxxxyy_1,  \
                             g_0_xxy_0_xxxxxyyy_0,  \
                             g_0_xxy_0_xxxxxyyy_1,  \
                             g_0_xxy_0_xxxxyyyy_0,  \
                             g_0_xxy_0_xxxxyyyy_1,  \
                             g_0_xxy_0_xxxyyyyy_0,  \
                             g_0_xxy_0_xxxyyyyy_1,  \
                             g_0_xxy_0_xxyyyyyy_0,  \
                             g_0_xxy_0_xxyyyyyy_1,  \
                             g_0_xxy_0_xyyyyyyy_0,  \
                             g_0_xxy_0_xyyyyyyy_1,  \
                             g_0_xxy_0_yyyyyyyy_0,  \
                             g_0_xxy_0_yyyyyyyy_1,  \
                             g_0_xxyz_0_xxxxxxxx_0, \
                             g_0_xxyz_0_xxxxxxxy_0, \
                             g_0_xxyz_0_xxxxxxxz_0, \
                             g_0_xxyz_0_xxxxxxyy_0, \
                             g_0_xxyz_0_xxxxxxyz_0, \
                             g_0_xxyz_0_xxxxxxzz_0, \
                             g_0_xxyz_0_xxxxxyyy_0, \
                             g_0_xxyz_0_xxxxxyyz_0, \
                             g_0_xxyz_0_xxxxxyzz_0, \
                             g_0_xxyz_0_xxxxxzzz_0, \
                             g_0_xxyz_0_xxxxyyyy_0, \
                             g_0_xxyz_0_xxxxyyyz_0, \
                             g_0_xxyz_0_xxxxyyzz_0, \
                             g_0_xxyz_0_xxxxyzzz_0, \
                             g_0_xxyz_0_xxxxzzzz_0, \
                             g_0_xxyz_0_xxxyyyyy_0, \
                             g_0_xxyz_0_xxxyyyyz_0, \
                             g_0_xxyz_0_xxxyyyzz_0, \
                             g_0_xxyz_0_xxxyyzzz_0, \
                             g_0_xxyz_0_xxxyzzzz_0, \
                             g_0_xxyz_0_xxxzzzzz_0, \
                             g_0_xxyz_0_xxyyyyyy_0, \
                             g_0_xxyz_0_xxyyyyyz_0, \
                             g_0_xxyz_0_xxyyyyzz_0, \
                             g_0_xxyz_0_xxyyyzzz_0, \
                             g_0_xxyz_0_xxyyzzzz_0, \
                             g_0_xxyz_0_xxyzzzzz_0, \
                             g_0_xxyz_0_xxzzzzzz_0, \
                             g_0_xxyz_0_xyyyyyyy_0, \
                             g_0_xxyz_0_xyyyyyyz_0, \
                             g_0_xxyz_0_xyyyyyzz_0, \
                             g_0_xxyz_0_xyyyyzzz_0, \
                             g_0_xxyz_0_xyyyzzzz_0, \
                             g_0_xxyz_0_xyyzzzzz_0, \
                             g_0_xxyz_0_xyzzzzzz_0, \
                             g_0_xxyz_0_xzzzzzzz_0, \
                             g_0_xxyz_0_yyyyyyyy_0, \
                             g_0_xxyz_0_yyyyyyyz_0, \
                             g_0_xxyz_0_yyyyyyzz_0, \
                             g_0_xxyz_0_yyyyyzzz_0, \
                             g_0_xxyz_0_yyyyzzzz_0, \
                             g_0_xxyz_0_yyyzzzzz_0, \
                             g_0_xxyz_0_yyzzzzzz_0, \
                             g_0_xxyz_0_yzzzzzzz_0, \
                             g_0_xxyz_0_zzzzzzzz_0, \
                             g_0_xxz_0_xxxxxxxx_0,  \
                             g_0_xxz_0_xxxxxxxx_1,  \
                             g_0_xxz_0_xxxxxxxz_0,  \
                             g_0_xxz_0_xxxxxxxz_1,  \
                             g_0_xxz_0_xxxxxxyz_0,  \
                             g_0_xxz_0_xxxxxxyz_1,  \
                             g_0_xxz_0_xxxxxxz_1,   \
                             g_0_xxz_0_xxxxxxzz_0,  \
                             g_0_xxz_0_xxxxxxzz_1,  \
                             g_0_xxz_0_xxxxxyyz_0,  \
                             g_0_xxz_0_xxxxxyyz_1,  \
                             g_0_xxz_0_xxxxxyz_1,   \
                             g_0_xxz_0_xxxxxyzz_0,  \
                             g_0_xxz_0_xxxxxyzz_1,  \
                             g_0_xxz_0_xxxxxzz_1,   \
                             g_0_xxz_0_xxxxxzzz_0,  \
                             g_0_xxz_0_xxxxxzzz_1,  \
                             g_0_xxz_0_xxxxyyyz_0,  \
                             g_0_xxz_0_xxxxyyyz_1,  \
                             g_0_xxz_0_xxxxyyz_1,   \
                             g_0_xxz_0_xxxxyyzz_0,  \
                             g_0_xxz_0_xxxxyyzz_1,  \
                             g_0_xxz_0_xxxxyzz_1,   \
                             g_0_xxz_0_xxxxyzzz_0,  \
                             g_0_xxz_0_xxxxyzzz_1,  \
                             g_0_xxz_0_xxxxzzz_1,   \
                             g_0_xxz_0_xxxxzzzz_0,  \
                             g_0_xxz_0_xxxxzzzz_1,  \
                             g_0_xxz_0_xxxyyyyz_0,  \
                             g_0_xxz_0_xxxyyyyz_1,  \
                             g_0_xxz_0_xxxyyyz_1,   \
                             g_0_xxz_0_xxxyyyzz_0,  \
                             g_0_xxz_0_xxxyyyzz_1,  \
                             g_0_xxz_0_xxxyyzz_1,   \
                             g_0_xxz_0_xxxyyzzz_0,  \
                             g_0_xxz_0_xxxyyzzz_1,  \
                             g_0_xxz_0_xxxyzzz_1,   \
                             g_0_xxz_0_xxxyzzzz_0,  \
                             g_0_xxz_0_xxxyzzzz_1,  \
                             g_0_xxz_0_xxxzzzz_1,   \
                             g_0_xxz_0_xxxzzzzz_0,  \
                             g_0_xxz_0_xxxzzzzz_1,  \
                             g_0_xxz_0_xxyyyyyz_0,  \
                             g_0_xxz_0_xxyyyyyz_1,  \
                             g_0_xxz_0_xxyyyyz_1,   \
                             g_0_xxz_0_xxyyyyzz_0,  \
                             g_0_xxz_0_xxyyyyzz_1,  \
                             g_0_xxz_0_xxyyyzz_1,   \
                             g_0_xxz_0_xxyyyzzz_0,  \
                             g_0_xxz_0_xxyyyzzz_1,  \
                             g_0_xxz_0_xxyyzzz_1,   \
                             g_0_xxz_0_xxyyzzzz_0,  \
                             g_0_xxz_0_xxyyzzzz_1,  \
                             g_0_xxz_0_xxyzzzz_1,   \
                             g_0_xxz_0_xxyzzzzz_0,  \
                             g_0_xxz_0_xxyzzzzz_1,  \
                             g_0_xxz_0_xxzzzzz_1,   \
                             g_0_xxz_0_xxzzzzzz_0,  \
                             g_0_xxz_0_xxzzzzzz_1,  \
                             g_0_xxz_0_xyyyyyyz_0,  \
                             g_0_xxz_0_xyyyyyyz_1,  \
                             g_0_xxz_0_xyyyyyz_1,   \
                             g_0_xxz_0_xyyyyyzz_0,  \
                             g_0_xxz_0_xyyyyyzz_1,  \
                             g_0_xxz_0_xyyyyzz_1,   \
                             g_0_xxz_0_xyyyyzzz_0,  \
                             g_0_xxz_0_xyyyyzzz_1,  \
                             g_0_xxz_0_xyyyzzz_1,   \
                             g_0_xxz_0_xyyyzzzz_0,  \
                             g_0_xxz_0_xyyyzzzz_1,  \
                             g_0_xxz_0_xyyzzzz_1,   \
                             g_0_xxz_0_xyyzzzzz_0,  \
                             g_0_xxz_0_xyyzzzzz_1,  \
                             g_0_xxz_0_xyzzzzz_1,   \
                             g_0_xxz_0_xyzzzzzz_0,  \
                             g_0_xxz_0_xyzzzzzz_1,  \
                             g_0_xxz_0_xzzzzzz_1,   \
                             g_0_xxz_0_xzzzzzzz_0,  \
                             g_0_xxz_0_xzzzzzzz_1,  \
                             g_0_xxz_0_yyyyyyyz_0,  \
                             g_0_xxz_0_yyyyyyyz_1,  \
                             g_0_xxz_0_yyyyyyz_1,   \
                             g_0_xxz_0_yyyyyyzz_0,  \
                             g_0_xxz_0_yyyyyyzz_1,  \
                             g_0_xxz_0_yyyyyzz_1,   \
                             g_0_xxz_0_yyyyyzzz_0,  \
                             g_0_xxz_0_yyyyyzzz_1,  \
                             g_0_xxz_0_yyyyzzz_1,   \
                             g_0_xxz_0_yyyyzzzz_0,  \
                             g_0_xxz_0_yyyyzzzz_1,  \
                             g_0_xxz_0_yyyzzzz_1,   \
                             g_0_xxz_0_yyyzzzzz_0,  \
                             g_0_xxz_0_yyyzzzzz_1,  \
                             g_0_xxz_0_yyzzzzz_1,   \
                             g_0_xxz_0_yyzzzzzz_0,  \
                             g_0_xxz_0_yyzzzzzz_1,  \
                             g_0_xxz_0_yzzzzzz_1,   \
                             g_0_xxz_0_yzzzzzzz_0,  \
                             g_0_xxz_0_yzzzzzzz_1,  \
                             g_0_xxz_0_zzzzzzz_1,   \
                             g_0_xxz_0_zzzzzzzz_0,  \
                             g_0_xxz_0_zzzzzzzz_1,  \
                             wp_y,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyz_0_xxxxxxxx_0[i] = g_0_xxz_0_xxxxxxxx_0[i] * pb_y + g_0_xxz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxxxy_0[i] = g_0_xxy_0_xxxxxxxy_0[i] * pb_z + g_0_xxy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxxxxz_0[i] = g_0_xxz_0_xxxxxxxz_0[i] * pb_y + g_0_xxz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxxyy_0[i] = g_0_xxy_0_xxxxxxyy_0[i] * pb_z + g_0_xxy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxxxyz_0[i] = g_0_xxz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxxxyz_0[i] * pb_y + g_0_xxz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxxzz_0[i] = g_0_xxz_0_xxxxxxzz_0[i] * pb_y + g_0_xxz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxyyy_0[i] = g_0_xxy_0_xxxxxyyy_0[i] * pb_z + g_0_xxy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxxyyz_0[i] = 2.0 * g_0_xxz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxxyyz_0[i] * pb_y + g_0_xxz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxyzz_0[i] = g_0_xxz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxxyzz_0[i] * pb_y + g_0_xxz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxzzz_0[i] = g_0_xxz_0_xxxxxzzz_0[i] * pb_y + g_0_xxz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxyyyy_0[i] = g_0_xxy_0_xxxxyyyy_0[i] * pb_z + g_0_xxy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxyyyz_0[i] = 3.0 * g_0_xxz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxyyyz_0[i] * pb_y + g_0_xxz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxyyzz_0[i] = 2.0 * g_0_xxz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxyyzz_0[i] * pb_y + g_0_xxz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxyzzz_0[i] = g_0_xxz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxyzzz_0[i] * pb_y + g_0_xxz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxzzzz_0[i] = g_0_xxz_0_xxxxzzzz_0[i] * pb_y + g_0_xxz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyyyyy_0[i] = g_0_xxy_0_xxxyyyyy_0[i] * pb_z + g_0_xxy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxyyyyz_0[i] = 4.0 * g_0_xxz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyyyyz_0[i] * pb_y + g_0_xxz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyyyzz_0[i] = 3.0 * g_0_xxz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyyyzz_0[i] * pb_y + g_0_xxz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyyzzz_0[i] = 2.0 * g_0_xxz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyyzzz_0[i] * pb_y + g_0_xxz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyzzzz_0[i] = g_0_xxz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyzzzz_0[i] * pb_y + g_0_xxz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxzzzzz_0[i] = g_0_xxz_0_xxxzzzzz_0[i] * pb_y + g_0_xxz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyyyyy_0[i] = g_0_xxy_0_xxyyyyyy_0[i] * pb_z + g_0_xxy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxyyyyyz_0[i] = 5.0 * g_0_xxz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyyyyz_0[i] * pb_y + g_0_xxz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyyyzz_0[i] = 4.0 * g_0_xxz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyyyzz_0[i] * pb_y + g_0_xxz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyyzzz_0[i] = 3.0 * g_0_xxz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyyzzz_0[i] * pb_y + g_0_xxz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyzzzz_0[i] = 2.0 * g_0_xxz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyzzzz_0[i] * pb_y + g_0_xxz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyzzzzz_0[i] = g_0_xxz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyzzzzz_0[i] * pb_y + g_0_xxz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxzzzzzz_0[i] = g_0_xxz_0_xxzzzzzz_0[i] * pb_y + g_0_xxz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyyyyy_0[i] = g_0_xxy_0_xyyyyyyy_0[i] * pb_z + g_0_xxy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xyyyyyyz_0[i] = 6.0 * g_0_xxz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyyyyz_0[i] * pb_y + g_0_xxz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyyyzz_0[i] = 5.0 * g_0_xxz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyyyzz_0[i] * pb_y + g_0_xxz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyyzzz_0[i] = 4.0 * g_0_xxz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyyzzz_0[i] * pb_y + g_0_xxz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyzzzz_0[i] = 3.0 * g_0_xxz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyzzzz_0[i] * pb_y + g_0_xxz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyzzzzz_0[i] = 2.0 * g_0_xxz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyzzzzz_0[i] * pb_y + g_0_xxz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyzzzzzz_0[i] = g_0_xxz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyzzzzzz_0[i] * pb_y + g_0_xxz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xzzzzzzz_0[i] = g_0_xxz_0_xzzzzzzz_0[i] * pb_y + g_0_xxz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyyyyy_0[i] = g_0_xxy_0_yyyyyyyy_0[i] * pb_z + g_0_xxy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_yyyyyyyz_0[i] = 7.0 * g_0_xxz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyyyyz_0[i] * pb_y + g_0_xxz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyyyzz_0[i] = 6.0 * g_0_xxz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyyyzz_0[i] * pb_y + g_0_xxz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyyzzz_0[i] = 5.0 * g_0_xxz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyyzzz_0[i] * pb_y + g_0_xxz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyzzzz_0[i] = 4.0 * g_0_xxz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyzzzz_0[i] * pb_y + g_0_xxz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyzzzzz_0[i] = 3.0 * g_0_xxz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyzzzzz_0[i] * pb_y + g_0_xxz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyzzzzzz_0[i] = 2.0 * g_0_xxz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyzzzzzz_0[i] * pb_y + g_0_xxz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yzzzzzzz_0[i] = g_0_xxz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yzzzzzzz_0[i] * pb_y + g_0_xxz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_zzzzzzzz_0[i] = g_0_xxz_0_zzzzzzzz_0[i] * pb_y + g_0_xxz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 225-270 components of targeted buffer : SGSL

    auto g_0_xxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 225);

    auto g_0_xxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 226);

    auto g_0_xxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 227);

    auto g_0_xxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 228);

    auto g_0_xxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 229);

    auto g_0_xxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 230);

    auto g_0_xxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 231);

    auto g_0_xxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 232);

    auto g_0_xxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 233);

    auto g_0_xxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 234);

    auto g_0_xxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 235);

    auto g_0_xxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 236);

    auto g_0_xxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 237);

    auto g_0_xxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 238);

    auto g_0_xxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 239);

    auto g_0_xxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 240);

    auto g_0_xxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 241);

    auto g_0_xxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 242);

    auto g_0_xxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 243);

    auto g_0_xxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 244);

    auto g_0_xxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 245);

    auto g_0_xxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 246);

    auto g_0_xxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 247);

    auto g_0_xxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 248);

    auto g_0_xxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 249);

    auto g_0_xxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 250);

    auto g_0_xxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 251);

    auto g_0_xxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 252);

    auto g_0_xxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 253);

    auto g_0_xxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 254);

    auto g_0_xxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 255);

    auto g_0_xxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 256);

    auto g_0_xxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 257);

    auto g_0_xxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 258);

    auto g_0_xxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 259);

    auto g_0_xxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 260);

    auto g_0_xxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 261);

    auto g_0_xxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 262);

    auto g_0_xxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 263);

    auto g_0_xxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 264);

    auto g_0_xxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 265);

    auto g_0_xxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 266);

    auto g_0_xxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 267);

    auto g_0_xxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 268);

    auto g_0_xxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 269);

#pragma omp simd aligned(g_0_xx_0_xxxxxxxx_0,       \
                             g_0_xx_0_xxxxxxxx_1,   \
                             g_0_xx_0_xxxxxxxy_0,   \
                             g_0_xx_0_xxxxxxxy_1,   \
                             g_0_xx_0_xxxxxxyy_0,   \
                             g_0_xx_0_xxxxxxyy_1,   \
                             g_0_xx_0_xxxxxyyy_0,   \
                             g_0_xx_0_xxxxxyyy_1,   \
                             g_0_xx_0_xxxxyyyy_0,   \
                             g_0_xx_0_xxxxyyyy_1,   \
                             g_0_xx_0_xxxyyyyy_0,   \
                             g_0_xx_0_xxxyyyyy_1,   \
                             g_0_xx_0_xxyyyyyy_0,   \
                             g_0_xx_0_xxyyyyyy_1,   \
                             g_0_xx_0_xyyyyyyy_0,   \
                             g_0_xx_0_xyyyyyyy_1,   \
                             g_0_xxz_0_xxxxxxxx_0,  \
                             g_0_xxz_0_xxxxxxxx_1,  \
                             g_0_xxz_0_xxxxxxxy_0,  \
                             g_0_xxz_0_xxxxxxxy_1,  \
                             g_0_xxz_0_xxxxxxyy_0,  \
                             g_0_xxz_0_xxxxxxyy_1,  \
                             g_0_xxz_0_xxxxxyyy_0,  \
                             g_0_xxz_0_xxxxxyyy_1,  \
                             g_0_xxz_0_xxxxyyyy_0,  \
                             g_0_xxz_0_xxxxyyyy_1,  \
                             g_0_xxz_0_xxxyyyyy_0,  \
                             g_0_xxz_0_xxxyyyyy_1,  \
                             g_0_xxz_0_xxyyyyyy_0,  \
                             g_0_xxz_0_xxyyyyyy_1,  \
                             g_0_xxz_0_xyyyyyyy_0,  \
                             g_0_xxz_0_xyyyyyyy_1,  \
                             g_0_xxzz_0_xxxxxxxx_0, \
                             g_0_xxzz_0_xxxxxxxy_0, \
                             g_0_xxzz_0_xxxxxxxz_0, \
                             g_0_xxzz_0_xxxxxxyy_0, \
                             g_0_xxzz_0_xxxxxxyz_0, \
                             g_0_xxzz_0_xxxxxxzz_0, \
                             g_0_xxzz_0_xxxxxyyy_0, \
                             g_0_xxzz_0_xxxxxyyz_0, \
                             g_0_xxzz_0_xxxxxyzz_0, \
                             g_0_xxzz_0_xxxxxzzz_0, \
                             g_0_xxzz_0_xxxxyyyy_0, \
                             g_0_xxzz_0_xxxxyyyz_0, \
                             g_0_xxzz_0_xxxxyyzz_0, \
                             g_0_xxzz_0_xxxxyzzz_0, \
                             g_0_xxzz_0_xxxxzzzz_0, \
                             g_0_xxzz_0_xxxyyyyy_0, \
                             g_0_xxzz_0_xxxyyyyz_0, \
                             g_0_xxzz_0_xxxyyyzz_0, \
                             g_0_xxzz_0_xxxyyzzz_0, \
                             g_0_xxzz_0_xxxyzzzz_0, \
                             g_0_xxzz_0_xxxzzzzz_0, \
                             g_0_xxzz_0_xxyyyyyy_0, \
                             g_0_xxzz_0_xxyyyyyz_0, \
                             g_0_xxzz_0_xxyyyyzz_0, \
                             g_0_xxzz_0_xxyyyzzz_0, \
                             g_0_xxzz_0_xxyyzzzz_0, \
                             g_0_xxzz_0_xxyzzzzz_0, \
                             g_0_xxzz_0_xxzzzzzz_0, \
                             g_0_xxzz_0_xyyyyyyy_0, \
                             g_0_xxzz_0_xyyyyyyz_0, \
                             g_0_xxzz_0_xyyyyyzz_0, \
                             g_0_xxzz_0_xyyyyzzz_0, \
                             g_0_xxzz_0_xyyyzzzz_0, \
                             g_0_xxzz_0_xyyzzzzz_0, \
                             g_0_xxzz_0_xyzzzzzz_0, \
                             g_0_xxzz_0_xzzzzzzz_0, \
                             g_0_xxzz_0_yyyyyyyy_0, \
                             g_0_xxzz_0_yyyyyyyz_0, \
                             g_0_xxzz_0_yyyyyyzz_0, \
                             g_0_xxzz_0_yyyyyzzz_0, \
                             g_0_xxzz_0_yyyyzzzz_0, \
                             g_0_xxzz_0_yyyzzzzz_0, \
                             g_0_xxzz_0_yyzzzzzz_0, \
                             g_0_xxzz_0_yzzzzzzz_0, \
                             g_0_xxzz_0_zzzzzzzz_0, \
                             g_0_xzz_0_xxxxxxxz_0,  \
                             g_0_xzz_0_xxxxxxxz_1,  \
                             g_0_xzz_0_xxxxxxyz_0,  \
                             g_0_xzz_0_xxxxxxyz_1,  \
                             g_0_xzz_0_xxxxxxz_1,   \
                             g_0_xzz_0_xxxxxxzz_0,  \
                             g_0_xzz_0_xxxxxxzz_1,  \
                             g_0_xzz_0_xxxxxyyz_0,  \
                             g_0_xzz_0_xxxxxyyz_1,  \
                             g_0_xzz_0_xxxxxyz_1,   \
                             g_0_xzz_0_xxxxxyzz_0,  \
                             g_0_xzz_0_xxxxxyzz_1,  \
                             g_0_xzz_0_xxxxxzz_1,   \
                             g_0_xzz_0_xxxxxzzz_0,  \
                             g_0_xzz_0_xxxxxzzz_1,  \
                             g_0_xzz_0_xxxxyyyz_0,  \
                             g_0_xzz_0_xxxxyyyz_1,  \
                             g_0_xzz_0_xxxxyyz_1,   \
                             g_0_xzz_0_xxxxyyzz_0,  \
                             g_0_xzz_0_xxxxyyzz_1,  \
                             g_0_xzz_0_xxxxyzz_1,   \
                             g_0_xzz_0_xxxxyzzz_0,  \
                             g_0_xzz_0_xxxxyzzz_1,  \
                             g_0_xzz_0_xxxxzzz_1,   \
                             g_0_xzz_0_xxxxzzzz_0,  \
                             g_0_xzz_0_xxxxzzzz_1,  \
                             g_0_xzz_0_xxxyyyyz_0,  \
                             g_0_xzz_0_xxxyyyyz_1,  \
                             g_0_xzz_0_xxxyyyz_1,   \
                             g_0_xzz_0_xxxyyyzz_0,  \
                             g_0_xzz_0_xxxyyyzz_1,  \
                             g_0_xzz_0_xxxyyzz_1,   \
                             g_0_xzz_0_xxxyyzzz_0,  \
                             g_0_xzz_0_xxxyyzzz_1,  \
                             g_0_xzz_0_xxxyzzz_1,   \
                             g_0_xzz_0_xxxyzzzz_0,  \
                             g_0_xzz_0_xxxyzzzz_1,  \
                             g_0_xzz_0_xxxzzzz_1,   \
                             g_0_xzz_0_xxxzzzzz_0,  \
                             g_0_xzz_0_xxxzzzzz_1,  \
                             g_0_xzz_0_xxyyyyyz_0,  \
                             g_0_xzz_0_xxyyyyyz_1,  \
                             g_0_xzz_0_xxyyyyz_1,   \
                             g_0_xzz_0_xxyyyyzz_0,  \
                             g_0_xzz_0_xxyyyyzz_1,  \
                             g_0_xzz_0_xxyyyzz_1,   \
                             g_0_xzz_0_xxyyyzzz_0,  \
                             g_0_xzz_0_xxyyyzzz_1,  \
                             g_0_xzz_0_xxyyzzz_1,   \
                             g_0_xzz_0_xxyyzzzz_0,  \
                             g_0_xzz_0_xxyyzzzz_1,  \
                             g_0_xzz_0_xxyzzzz_1,   \
                             g_0_xzz_0_xxyzzzzz_0,  \
                             g_0_xzz_0_xxyzzzzz_1,  \
                             g_0_xzz_0_xxzzzzz_1,   \
                             g_0_xzz_0_xxzzzzzz_0,  \
                             g_0_xzz_0_xxzzzzzz_1,  \
                             g_0_xzz_0_xyyyyyyz_0,  \
                             g_0_xzz_0_xyyyyyyz_1,  \
                             g_0_xzz_0_xyyyyyz_1,   \
                             g_0_xzz_0_xyyyyyzz_0,  \
                             g_0_xzz_0_xyyyyyzz_1,  \
                             g_0_xzz_0_xyyyyzz_1,   \
                             g_0_xzz_0_xyyyyzzz_0,  \
                             g_0_xzz_0_xyyyyzzz_1,  \
                             g_0_xzz_0_xyyyzzz_1,   \
                             g_0_xzz_0_xyyyzzzz_0,  \
                             g_0_xzz_0_xyyyzzzz_1,  \
                             g_0_xzz_0_xyyzzzz_1,   \
                             g_0_xzz_0_xyyzzzzz_0,  \
                             g_0_xzz_0_xyyzzzzz_1,  \
                             g_0_xzz_0_xyzzzzz_1,   \
                             g_0_xzz_0_xyzzzzzz_0,  \
                             g_0_xzz_0_xyzzzzzz_1,  \
                             g_0_xzz_0_xzzzzzz_1,   \
                             g_0_xzz_0_xzzzzzzz_0,  \
                             g_0_xzz_0_xzzzzzzz_1,  \
                             g_0_xzz_0_yyyyyyyy_0,  \
                             g_0_xzz_0_yyyyyyyy_1,  \
                             g_0_xzz_0_yyyyyyyz_0,  \
                             g_0_xzz_0_yyyyyyyz_1,  \
                             g_0_xzz_0_yyyyyyz_1,   \
                             g_0_xzz_0_yyyyyyzz_0,  \
                             g_0_xzz_0_yyyyyyzz_1,  \
                             g_0_xzz_0_yyyyyzz_1,   \
                             g_0_xzz_0_yyyyyzzz_0,  \
                             g_0_xzz_0_yyyyyzzz_1,  \
                             g_0_xzz_0_yyyyzzz_1,   \
                             g_0_xzz_0_yyyyzzzz_0,  \
                             g_0_xzz_0_yyyyzzzz_1,  \
                             g_0_xzz_0_yyyzzzz_1,   \
                             g_0_xzz_0_yyyzzzzz_0,  \
                             g_0_xzz_0_yyyzzzzz_1,  \
                             g_0_xzz_0_yyzzzzz_1,   \
                             g_0_xzz_0_yyzzzzzz_0,  \
                             g_0_xzz_0_yyzzzzzz_1,  \
                             g_0_xzz_0_yzzzzzz_1,   \
                             g_0_xzz_0_yzzzzzzz_0,  \
                             g_0_xzz_0_yzzzzzzz_1,  \
                             g_0_xzz_0_zzzzzzz_1,   \
                             g_0_xzz_0_zzzzzzzz_0,  \
                             g_0_xzz_0_zzzzzzzz_1,  \
                             g_0_zz_0_xxxxxxxz_0,   \
                             g_0_zz_0_xxxxxxxz_1,   \
                             g_0_zz_0_xxxxxxyz_0,   \
                             g_0_zz_0_xxxxxxyz_1,   \
                             g_0_zz_0_xxxxxxzz_0,   \
                             g_0_zz_0_xxxxxxzz_1,   \
                             g_0_zz_0_xxxxxyyz_0,   \
                             g_0_zz_0_xxxxxyyz_1,   \
                             g_0_zz_0_xxxxxyzz_0,   \
                             g_0_zz_0_xxxxxyzz_1,   \
                             g_0_zz_0_xxxxxzzz_0,   \
                             g_0_zz_0_xxxxxzzz_1,   \
                             g_0_zz_0_xxxxyyyz_0,   \
                             g_0_zz_0_xxxxyyyz_1,   \
                             g_0_zz_0_xxxxyyzz_0,   \
                             g_0_zz_0_xxxxyyzz_1,   \
                             g_0_zz_0_xxxxyzzz_0,   \
                             g_0_zz_0_xxxxyzzz_1,   \
                             g_0_zz_0_xxxxzzzz_0,   \
                             g_0_zz_0_xxxxzzzz_1,   \
                             g_0_zz_0_xxxyyyyz_0,   \
                             g_0_zz_0_xxxyyyyz_1,   \
                             g_0_zz_0_xxxyyyzz_0,   \
                             g_0_zz_0_xxxyyyzz_1,   \
                             g_0_zz_0_xxxyyzzz_0,   \
                             g_0_zz_0_xxxyyzzz_1,   \
                             g_0_zz_0_xxxyzzzz_0,   \
                             g_0_zz_0_xxxyzzzz_1,   \
                             g_0_zz_0_xxxzzzzz_0,   \
                             g_0_zz_0_xxxzzzzz_1,   \
                             g_0_zz_0_xxyyyyyz_0,   \
                             g_0_zz_0_xxyyyyyz_1,   \
                             g_0_zz_0_xxyyyyzz_0,   \
                             g_0_zz_0_xxyyyyzz_1,   \
                             g_0_zz_0_xxyyyzzz_0,   \
                             g_0_zz_0_xxyyyzzz_1,   \
                             g_0_zz_0_xxyyzzzz_0,   \
                             g_0_zz_0_xxyyzzzz_1,   \
                             g_0_zz_0_xxyzzzzz_0,   \
                             g_0_zz_0_xxyzzzzz_1,   \
                             g_0_zz_0_xxzzzzzz_0,   \
                             g_0_zz_0_xxzzzzzz_1,   \
                             g_0_zz_0_xyyyyyyz_0,   \
                             g_0_zz_0_xyyyyyyz_1,   \
                             g_0_zz_0_xyyyyyzz_0,   \
                             g_0_zz_0_xyyyyyzz_1,   \
                             g_0_zz_0_xyyyyzzz_0,   \
                             g_0_zz_0_xyyyyzzz_1,   \
                             g_0_zz_0_xyyyzzzz_0,   \
                             g_0_zz_0_xyyyzzzz_1,   \
                             g_0_zz_0_xyyzzzzz_0,   \
                             g_0_zz_0_xyyzzzzz_1,   \
                             g_0_zz_0_xyzzzzzz_0,   \
                             g_0_zz_0_xyzzzzzz_1,   \
                             g_0_zz_0_xzzzzzzz_0,   \
                             g_0_zz_0_xzzzzzzz_1,   \
                             g_0_zz_0_yyyyyyyy_0,   \
                             g_0_zz_0_yyyyyyyy_1,   \
                             g_0_zz_0_yyyyyyyz_0,   \
                             g_0_zz_0_yyyyyyyz_1,   \
                             g_0_zz_0_yyyyyyzz_0,   \
                             g_0_zz_0_yyyyyyzz_1,   \
                             g_0_zz_0_yyyyyzzz_0,   \
                             g_0_zz_0_yyyyyzzz_1,   \
                             g_0_zz_0_yyyyzzzz_0,   \
                             g_0_zz_0_yyyyzzzz_1,   \
                             g_0_zz_0_yyyzzzzz_0,   \
                             g_0_zz_0_yyyzzzzz_1,   \
                             g_0_zz_0_yyzzzzzz_0,   \
                             g_0_zz_0_yyzzzzzz_1,   \
                             g_0_zz_0_yzzzzzzz_0,   \
                             g_0_zz_0_yzzzzzzz_1,   \
                             g_0_zz_0_zzzzzzzz_0,   \
                             g_0_zz_0_zzzzzzzz_1,   \
                             wp_x,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzz_0_xxxxxxxx_0[i] =
            g_0_xx_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxxxx_0[i] * pb_z + g_0_xxz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxxxy_0[i] =
            g_0_xx_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxxxy_0[i] * pb_z + g_0_xxz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxxxz_0[i] = g_0_zz_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxxz_1[i] * fti_ab_0 + 7.0 * g_0_xzz_0_xxxxxxz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxxxxz_0[i] * pb_x + g_0_xzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxxxyy_0[i] =
            g_0_xx_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxxyy_0[i] * pb_z + g_0_xxz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxxyz_0[i] = g_0_zz_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_xzz_0_xxxxxyz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxxxyz_0[i] * pb_x + g_0_xzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxxxzz_0[i] = g_0_zz_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxzz_1[i] * fti_ab_0 + 6.0 * g_0_xzz_0_xxxxxzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxxxzz_0[i] * pb_x + g_0_xzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxxyyy_0[i] =
            g_0_xx_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_xx_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxyyy_0[i] * pb_z + g_0_xxz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxyyz_0[i] = g_0_zz_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_xzz_0_xxxxyyz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxxyyz_0[i] * pb_x + g_0_xzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxxyzz_0[i] = g_0_zz_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_xzz_0_xxxxyzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxxyzz_0[i] * pb_x + g_0_xzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxxzzz_0[i] = g_0_zz_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxzzz_1[i] * fti_ab_0 + 5.0 * g_0_xzz_0_xxxxzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxxzzz_0[i] * pb_x + g_0_xzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxyyyy_0[i] =
            g_0_xx_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_xx_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxyyyy_0[i] * pb_z + g_0_xxz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxyyyz_0[i] = g_0_zz_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxyyyz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxyyyz_0[i] * pb_x + g_0_xzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxyyzz_0[i] = g_0_zz_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxyyzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxyyzz_0[i] * pb_x + g_0_xzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxyzzz_0[i] = g_0_zz_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxyzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxyzzz_0[i] * pb_x + g_0_xzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxzzzz_0[i] = g_0_zz_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxxzzzz_0[i] * pb_x + g_0_xzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyyyyy_0[i] =
            g_0_xx_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_xx_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxyyyyy_0[i] * pb_z + g_0_xxz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxyyyyz_0[i] = g_0_zz_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyyyyz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxyyyyz_0[i] * pb_x + g_0_xzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyyyzz_0[i] = g_0_zz_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyyyzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxyyyzz_0[i] * pb_x + g_0_xzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyyzzz_0[i] = g_0_zz_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyyzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxyyzzz_0[i] * pb_x + g_0_xzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyzzzz_0[i] = g_0_zz_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxyzzzz_0[i] * pb_x + g_0_xzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxzzzzz_0[i] = g_0_zz_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxxzzzzz_0[i] * pb_x + g_0_xzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyyyyy_0[i] =
            g_0_xx_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_xx_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxyyyyyy_0[i] * pb_z + g_0_xxz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxyyyyyz_0[i] = g_0_zz_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyyyyz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxyyyyyz_0[i] * pb_x + g_0_xzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyyyzz_0[i] = g_0_zz_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyyyzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxyyyyzz_0[i] * pb_x + g_0_xzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyyzzz_0[i] = g_0_zz_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyyzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxyyyzzz_0[i] * pb_x + g_0_xzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyzzzz_0[i] = g_0_zz_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxyyzzzz_0[i] * pb_x + g_0_xzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyzzzzz_0[i] = g_0_zz_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxyzzzzz_0[i] * pb_x + g_0_xzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxzzzzzz_0[i] = g_0_zz_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xzzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xxzzzzzz_0[i] * pb_x + g_0_xzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyyyyy_0[i] =
            g_0_xx_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_xx_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xyyyyyyy_0[i] * pb_z + g_0_xxz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xyyyyyyz_0[i] = g_0_zz_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyyz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xyyyyyyz_0[i] * pb_x + g_0_xzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyyyzz_0[i] = g_0_zz_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xyyyyyzz_0[i] * pb_x + g_0_xzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyyzzz_0[i] = g_0_zz_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xyyyyzzz_0[i] * pb_x + g_0_xzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyzzzz_0[i] = g_0_zz_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xyyyzzzz_0[i] * pb_x + g_0_xzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyzzzzz_0[i] = g_0_zz_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xyyzzzzz_0[i] * pb_x + g_0_xzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyzzzzzz_0[i] = g_0_zz_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xyzzzzzz_0[i] * pb_x + g_0_xzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xzzzzzzz_0[i] = g_0_zz_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzz_0_xzzzzzzz_0[i] * pb_x + g_0_xzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyyyy_0[i] =
            g_0_zz_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyyyy_0[i] * pb_x + g_0_xzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyyyz_0[i] =
            g_0_zz_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyyyz_0[i] * pb_x + g_0_xzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyyzz_0[i] =
            g_0_zz_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyyzz_0[i] * pb_x + g_0_xzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyzzz_0[i] =
            g_0_zz_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyzzz_0[i] * pb_x + g_0_xzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyzzzz_0[i] =
            g_0_zz_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyzzzz_0[i] * pb_x + g_0_xzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyzzzzz_0[i] =
            g_0_zz_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyzzzzz_0[i] * pb_x + g_0_xzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyzzzzzz_0[i] =
            g_0_zz_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzzzzzz_0[i] * pb_x + g_0_xzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yzzzzzzz_0[i] =
            g_0_zz_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzzzzzz_0[i] * pb_x + g_0_xzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_zzzzzzzz_0[i] =
            g_0_zz_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzzzzzz_0[i] * pb_x + g_0_xzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 270-315 components of targeted buffer : SGSL

    auto g_0_xyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 270);

    auto g_0_xyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 271);

    auto g_0_xyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 272);

    auto g_0_xyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 273);

    auto g_0_xyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 274);

    auto g_0_xyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 275);

    auto g_0_xyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 276);

    auto g_0_xyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 277);

    auto g_0_xyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 278);

    auto g_0_xyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 279);

    auto g_0_xyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 280);

    auto g_0_xyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 281);

    auto g_0_xyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 282);

    auto g_0_xyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 283);

    auto g_0_xyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 284);

    auto g_0_xyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 285);

    auto g_0_xyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 286);

    auto g_0_xyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 287);

    auto g_0_xyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 288);

    auto g_0_xyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 289);

    auto g_0_xyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 290);

    auto g_0_xyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 291);

    auto g_0_xyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 292);

    auto g_0_xyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 293);

    auto g_0_xyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 294);

    auto g_0_xyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 295);

    auto g_0_xyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 296);

    auto g_0_xyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 297);

    auto g_0_xyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 298);

    auto g_0_xyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 299);

    auto g_0_xyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 300);

    auto g_0_xyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 301);

    auto g_0_xyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 302);

    auto g_0_xyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 303);

    auto g_0_xyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 304);

    auto g_0_xyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 305);

    auto g_0_xyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 306);

    auto g_0_xyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 307);

    auto g_0_xyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 308);

    auto g_0_xyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 309);

    auto g_0_xyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 310);

    auto g_0_xyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 311);

    auto g_0_xyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 312);

    auto g_0_xyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 313);

    auto g_0_xyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 314);

#pragma omp simd aligned(g_0_xyyy_0_xxxxxxxx_0,     \
                             g_0_xyyy_0_xxxxxxxy_0, \
                             g_0_xyyy_0_xxxxxxxz_0, \
                             g_0_xyyy_0_xxxxxxyy_0, \
                             g_0_xyyy_0_xxxxxxyz_0, \
                             g_0_xyyy_0_xxxxxxzz_0, \
                             g_0_xyyy_0_xxxxxyyy_0, \
                             g_0_xyyy_0_xxxxxyyz_0, \
                             g_0_xyyy_0_xxxxxyzz_0, \
                             g_0_xyyy_0_xxxxxzzz_0, \
                             g_0_xyyy_0_xxxxyyyy_0, \
                             g_0_xyyy_0_xxxxyyyz_0, \
                             g_0_xyyy_0_xxxxyyzz_0, \
                             g_0_xyyy_0_xxxxyzzz_0, \
                             g_0_xyyy_0_xxxxzzzz_0, \
                             g_0_xyyy_0_xxxyyyyy_0, \
                             g_0_xyyy_0_xxxyyyyz_0, \
                             g_0_xyyy_0_xxxyyyzz_0, \
                             g_0_xyyy_0_xxxyyzzz_0, \
                             g_0_xyyy_0_xxxyzzzz_0, \
                             g_0_xyyy_0_xxxzzzzz_0, \
                             g_0_xyyy_0_xxyyyyyy_0, \
                             g_0_xyyy_0_xxyyyyyz_0, \
                             g_0_xyyy_0_xxyyyyzz_0, \
                             g_0_xyyy_0_xxyyyzzz_0, \
                             g_0_xyyy_0_xxyyzzzz_0, \
                             g_0_xyyy_0_xxyzzzzz_0, \
                             g_0_xyyy_0_xxzzzzzz_0, \
                             g_0_xyyy_0_xyyyyyyy_0, \
                             g_0_xyyy_0_xyyyyyyz_0, \
                             g_0_xyyy_0_xyyyyyzz_0, \
                             g_0_xyyy_0_xyyyyzzz_0, \
                             g_0_xyyy_0_xyyyzzzz_0, \
                             g_0_xyyy_0_xyyzzzzz_0, \
                             g_0_xyyy_0_xyzzzzzz_0, \
                             g_0_xyyy_0_xzzzzzzz_0, \
                             g_0_xyyy_0_yyyyyyyy_0, \
                             g_0_xyyy_0_yyyyyyyz_0, \
                             g_0_xyyy_0_yyyyyyzz_0, \
                             g_0_xyyy_0_yyyyyzzz_0, \
                             g_0_xyyy_0_yyyyzzzz_0, \
                             g_0_xyyy_0_yyyzzzzz_0, \
                             g_0_xyyy_0_yyzzzzzz_0, \
                             g_0_xyyy_0_yzzzzzzz_0, \
                             g_0_xyyy_0_zzzzzzzz_0, \
                             g_0_yyy_0_xxxxxxx_1,   \
                             g_0_yyy_0_xxxxxxxx_0,  \
                             g_0_yyy_0_xxxxxxxx_1,  \
                             g_0_yyy_0_xxxxxxxy_0,  \
                             g_0_yyy_0_xxxxxxxy_1,  \
                             g_0_yyy_0_xxxxxxxz_0,  \
                             g_0_yyy_0_xxxxxxxz_1,  \
                             g_0_yyy_0_xxxxxxy_1,   \
                             g_0_yyy_0_xxxxxxyy_0,  \
                             g_0_yyy_0_xxxxxxyy_1,  \
                             g_0_yyy_0_xxxxxxyz_0,  \
                             g_0_yyy_0_xxxxxxyz_1,  \
                             g_0_yyy_0_xxxxxxz_1,   \
                             g_0_yyy_0_xxxxxxzz_0,  \
                             g_0_yyy_0_xxxxxxzz_1,  \
                             g_0_yyy_0_xxxxxyy_1,   \
                             g_0_yyy_0_xxxxxyyy_0,  \
                             g_0_yyy_0_xxxxxyyy_1,  \
                             g_0_yyy_0_xxxxxyyz_0,  \
                             g_0_yyy_0_xxxxxyyz_1,  \
                             g_0_yyy_0_xxxxxyz_1,   \
                             g_0_yyy_0_xxxxxyzz_0,  \
                             g_0_yyy_0_xxxxxyzz_1,  \
                             g_0_yyy_0_xxxxxzz_1,   \
                             g_0_yyy_0_xxxxxzzz_0,  \
                             g_0_yyy_0_xxxxxzzz_1,  \
                             g_0_yyy_0_xxxxyyy_1,   \
                             g_0_yyy_0_xxxxyyyy_0,  \
                             g_0_yyy_0_xxxxyyyy_1,  \
                             g_0_yyy_0_xxxxyyyz_0,  \
                             g_0_yyy_0_xxxxyyyz_1,  \
                             g_0_yyy_0_xxxxyyz_1,   \
                             g_0_yyy_0_xxxxyyzz_0,  \
                             g_0_yyy_0_xxxxyyzz_1,  \
                             g_0_yyy_0_xxxxyzz_1,   \
                             g_0_yyy_0_xxxxyzzz_0,  \
                             g_0_yyy_0_xxxxyzzz_1,  \
                             g_0_yyy_0_xxxxzzz_1,   \
                             g_0_yyy_0_xxxxzzzz_0,  \
                             g_0_yyy_0_xxxxzzzz_1,  \
                             g_0_yyy_0_xxxyyyy_1,   \
                             g_0_yyy_0_xxxyyyyy_0,  \
                             g_0_yyy_0_xxxyyyyy_1,  \
                             g_0_yyy_0_xxxyyyyz_0,  \
                             g_0_yyy_0_xxxyyyyz_1,  \
                             g_0_yyy_0_xxxyyyz_1,   \
                             g_0_yyy_0_xxxyyyzz_0,  \
                             g_0_yyy_0_xxxyyyzz_1,  \
                             g_0_yyy_0_xxxyyzz_1,   \
                             g_0_yyy_0_xxxyyzzz_0,  \
                             g_0_yyy_0_xxxyyzzz_1,  \
                             g_0_yyy_0_xxxyzzz_1,   \
                             g_0_yyy_0_xxxyzzzz_0,  \
                             g_0_yyy_0_xxxyzzzz_1,  \
                             g_0_yyy_0_xxxzzzz_1,   \
                             g_0_yyy_0_xxxzzzzz_0,  \
                             g_0_yyy_0_xxxzzzzz_1,  \
                             g_0_yyy_0_xxyyyyy_1,   \
                             g_0_yyy_0_xxyyyyyy_0,  \
                             g_0_yyy_0_xxyyyyyy_1,  \
                             g_0_yyy_0_xxyyyyyz_0,  \
                             g_0_yyy_0_xxyyyyyz_1,  \
                             g_0_yyy_0_xxyyyyz_1,   \
                             g_0_yyy_0_xxyyyyzz_0,  \
                             g_0_yyy_0_xxyyyyzz_1,  \
                             g_0_yyy_0_xxyyyzz_1,   \
                             g_0_yyy_0_xxyyyzzz_0,  \
                             g_0_yyy_0_xxyyyzzz_1,  \
                             g_0_yyy_0_xxyyzzz_1,   \
                             g_0_yyy_0_xxyyzzzz_0,  \
                             g_0_yyy_0_xxyyzzzz_1,  \
                             g_0_yyy_0_xxyzzzz_1,   \
                             g_0_yyy_0_xxyzzzzz_0,  \
                             g_0_yyy_0_xxyzzzzz_1,  \
                             g_0_yyy_0_xxzzzzz_1,   \
                             g_0_yyy_0_xxzzzzzz_0,  \
                             g_0_yyy_0_xxzzzzzz_1,  \
                             g_0_yyy_0_xyyyyyy_1,   \
                             g_0_yyy_0_xyyyyyyy_0,  \
                             g_0_yyy_0_xyyyyyyy_1,  \
                             g_0_yyy_0_xyyyyyyz_0,  \
                             g_0_yyy_0_xyyyyyyz_1,  \
                             g_0_yyy_0_xyyyyyz_1,   \
                             g_0_yyy_0_xyyyyyzz_0,  \
                             g_0_yyy_0_xyyyyyzz_1,  \
                             g_0_yyy_0_xyyyyzz_1,   \
                             g_0_yyy_0_xyyyyzzz_0,  \
                             g_0_yyy_0_xyyyyzzz_1,  \
                             g_0_yyy_0_xyyyzzz_1,   \
                             g_0_yyy_0_xyyyzzzz_0,  \
                             g_0_yyy_0_xyyyzzzz_1,  \
                             g_0_yyy_0_xyyzzzz_1,   \
                             g_0_yyy_0_xyyzzzzz_0,  \
                             g_0_yyy_0_xyyzzzzz_1,  \
                             g_0_yyy_0_xyzzzzz_1,   \
                             g_0_yyy_0_xyzzzzzz_0,  \
                             g_0_yyy_0_xyzzzzzz_1,  \
                             g_0_yyy_0_xzzzzzz_1,   \
                             g_0_yyy_0_xzzzzzzz_0,  \
                             g_0_yyy_0_xzzzzzzz_1,  \
                             g_0_yyy_0_yyyyyyy_1,   \
                             g_0_yyy_0_yyyyyyyy_0,  \
                             g_0_yyy_0_yyyyyyyy_1,  \
                             g_0_yyy_0_yyyyyyyz_0,  \
                             g_0_yyy_0_yyyyyyyz_1,  \
                             g_0_yyy_0_yyyyyyz_1,   \
                             g_0_yyy_0_yyyyyyzz_0,  \
                             g_0_yyy_0_yyyyyyzz_1,  \
                             g_0_yyy_0_yyyyyzz_1,   \
                             g_0_yyy_0_yyyyyzzz_0,  \
                             g_0_yyy_0_yyyyyzzz_1,  \
                             g_0_yyy_0_yyyyzzz_1,   \
                             g_0_yyy_0_yyyyzzzz_0,  \
                             g_0_yyy_0_yyyyzzzz_1,  \
                             g_0_yyy_0_yyyzzzz_1,   \
                             g_0_yyy_0_yyyzzzzz_0,  \
                             g_0_yyy_0_yyyzzzzz_1,  \
                             g_0_yyy_0_yyzzzzz_1,   \
                             g_0_yyy_0_yyzzzzzz_0,  \
                             g_0_yyy_0_yyzzzzzz_1,  \
                             g_0_yyy_0_yzzzzzz_1,   \
                             g_0_yyy_0_yzzzzzzz_0,  \
                             g_0_yyy_0_yzzzzzzz_1,  \
                             g_0_yyy_0_zzzzzzz_1,   \
                             g_0_yyy_0_zzzzzzzz_0,  \
                             g_0_yyy_0_zzzzzzzz_1,  \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyy_0_xxxxxxxx_0[i] = 8.0 * g_0_yyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxxx_0[i] * pb_x + g_0_yyy_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxxxy_0[i] = 7.0 * g_0_yyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxxy_0[i] * pb_x + g_0_yyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxxxz_0[i] = 7.0 * g_0_yyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxxz_0[i] * pb_x + g_0_yyy_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxxyy_0[i] = 6.0 * g_0_yyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxyy_0[i] * pb_x + g_0_yyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxxyz_0[i] = 6.0 * g_0_yyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxyz_0[i] * pb_x + g_0_yyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxxzz_0[i] = 6.0 * g_0_yyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxzz_0[i] * pb_x + g_0_yyy_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxyyy_0[i] = 5.0 * g_0_yyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyyy_0[i] * pb_x + g_0_yyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxyyz_0[i] = 5.0 * g_0_yyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyyz_0[i] * pb_x + g_0_yyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxyzz_0[i] = 5.0 * g_0_yyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyzz_0[i] * pb_x + g_0_yyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxzzz_0[i] = 5.0 * g_0_yyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxzzz_0[i] * pb_x + g_0_yyy_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyyyy_0[i] = 4.0 * g_0_yyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyyy_0[i] * pb_x + g_0_yyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyyyz_0[i] = 4.0 * g_0_yyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyyz_0[i] * pb_x + g_0_yyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyyzz_0[i] = 4.0 * g_0_yyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyzz_0[i] * pb_x + g_0_yyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyzzz_0[i] = 4.0 * g_0_yyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyzzz_0[i] * pb_x + g_0_yyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxzzzz_0[i] = 4.0 * g_0_yyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxzzzz_0[i] * pb_x + g_0_yyy_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyyyy_0[i] = 3.0 * g_0_yyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyyy_0[i] * pb_x + g_0_yyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyyyz_0[i] = 3.0 * g_0_yyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyyz_0[i] * pb_x + g_0_yyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyyzz_0[i] = 3.0 * g_0_yyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyzz_0[i] * pb_x + g_0_yyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyzzz_0[i] = 3.0 * g_0_yyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyzzz_0[i] * pb_x + g_0_yyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyzzzz_0[i] = 3.0 * g_0_yyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzzzz_0[i] * pb_x + g_0_yyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxzzzzz_0[i] = 3.0 * g_0_yyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzzzzz_0[i] * pb_x + g_0_yyy_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyyyy_0[i] = 2.0 * g_0_yyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyyy_0[i] * pb_x + g_0_yyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyyyz_0[i] = 2.0 * g_0_yyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyyz_0[i] * pb_x + g_0_yyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyyzz_0[i] = 2.0 * g_0_yyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyzz_0[i] * pb_x + g_0_yyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyzzz_0[i] = 2.0 * g_0_yyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyzzz_0[i] * pb_x + g_0_yyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyzzzz_0[i] = 2.0 * g_0_yyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzzzz_0[i] * pb_x + g_0_yyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyzzzzz_0[i] = 2.0 * g_0_yyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzzzz_0[i] * pb_x + g_0_yyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxzzzzzz_0[i] = 2.0 * g_0_yyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzzzzz_0[i] * pb_x + g_0_yyy_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyyyy_0[i] = g_0_yyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyyy_0[i] * pb_x + g_0_yyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyyyz_0[i] = g_0_yyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyyz_0[i] * pb_x + g_0_yyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyyzz_0[i] = g_0_yyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyzz_0[i] * pb_x + g_0_yyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyzzz_0[i] = g_0_yyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyzzz_0[i] * pb_x + g_0_yyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyzzzz_0[i] = g_0_yyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzzzz_0[i] * pb_x + g_0_yyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyzzzzz_0[i] = g_0_yyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzzzz_0[i] * pb_x + g_0_yyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyzzzzzz_0[i] = g_0_yyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzzzz_0[i] * pb_x + g_0_yyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xzzzzzzz_0[i] = g_0_yyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzzzzz_0[i] * pb_x + g_0_yyy_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyyyy_0[i] = g_0_yyy_0_yyyyyyyy_0[i] * pb_x + g_0_yyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyyyz_0[i] = g_0_yyy_0_yyyyyyyz_0[i] * pb_x + g_0_yyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyyzz_0[i] = g_0_yyy_0_yyyyyyzz_0[i] * pb_x + g_0_yyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyzzz_0[i] = g_0_yyy_0_yyyyyzzz_0[i] * pb_x + g_0_yyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyzzzz_0[i] = g_0_yyy_0_yyyyzzzz_0[i] * pb_x + g_0_yyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyzzzzz_0[i] = g_0_yyy_0_yyyzzzzz_0[i] * pb_x + g_0_yyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyzzzzzz_0[i] = g_0_yyy_0_yyzzzzzz_0[i] * pb_x + g_0_yyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yzzzzzzz_0[i] = g_0_yyy_0_yzzzzzzz_0[i] * pb_x + g_0_yyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_zzzzzzzz_0[i] = g_0_yyy_0_zzzzzzzz_0[i] * pb_x + g_0_yyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 315-360 components of targeted buffer : SGSL

    auto g_0_xyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 315);

    auto g_0_xyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 316);

    auto g_0_xyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 317);

    auto g_0_xyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 318);

    auto g_0_xyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 319);

    auto g_0_xyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 320);

    auto g_0_xyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 321);

    auto g_0_xyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 322);

    auto g_0_xyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 323);

    auto g_0_xyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 324);

    auto g_0_xyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 325);

    auto g_0_xyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 326);

    auto g_0_xyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 327);

    auto g_0_xyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 328);

    auto g_0_xyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 329);

    auto g_0_xyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 330);

    auto g_0_xyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 331);

    auto g_0_xyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 332);

    auto g_0_xyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 333);

    auto g_0_xyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 334);

    auto g_0_xyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 335);

    auto g_0_xyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 336);

    auto g_0_xyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 337);

    auto g_0_xyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 338);

    auto g_0_xyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 339);

    auto g_0_xyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 340);

    auto g_0_xyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 341);

    auto g_0_xyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 342);

    auto g_0_xyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 343);

    auto g_0_xyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 344);

    auto g_0_xyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 345);

    auto g_0_xyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 346);

    auto g_0_xyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 347);

    auto g_0_xyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 348);

    auto g_0_xyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 349);

    auto g_0_xyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 350);

    auto g_0_xyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 351);

    auto g_0_xyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 352);

    auto g_0_xyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 353);

    auto g_0_xyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 354);

    auto g_0_xyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 355);

    auto g_0_xyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 356);

    auto g_0_xyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 357);

    auto g_0_xyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 358);

    auto g_0_xyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 359);

#pragma omp simd aligned(g_0_xyy_0_xxxxxxxx_0,      \
                             g_0_xyy_0_xxxxxxxx_1,  \
                             g_0_xyy_0_xxxxxxxy_0,  \
                             g_0_xyy_0_xxxxxxxy_1,  \
                             g_0_xyy_0_xxxxxxyy_0,  \
                             g_0_xyy_0_xxxxxxyy_1,  \
                             g_0_xyy_0_xxxxxyyy_0,  \
                             g_0_xyy_0_xxxxxyyy_1,  \
                             g_0_xyy_0_xxxxyyyy_0,  \
                             g_0_xyy_0_xxxxyyyy_1,  \
                             g_0_xyy_0_xxxyyyyy_0,  \
                             g_0_xyy_0_xxxyyyyy_1,  \
                             g_0_xyy_0_xxyyyyyy_0,  \
                             g_0_xyy_0_xxyyyyyy_1,  \
                             g_0_xyy_0_xyyyyyyy_0,  \
                             g_0_xyy_0_xyyyyyyy_1,  \
                             g_0_xyyz_0_xxxxxxxx_0, \
                             g_0_xyyz_0_xxxxxxxy_0, \
                             g_0_xyyz_0_xxxxxxxz_0, \
                             g_0_xyyz_0_xxxxxxyy_0, \
                             g_0_xyyz_0_xxxxxxyz_0, \
                             g_0_xyyz_0_xxxxxxzz_0, \
                             g_0_xyyz_0_xxxxxyyy_0, \
                             g_0_xyyz_0_xxxxxyyz_0, \
                             g_0_xyyz_0_xxxxxyzz_0, \
                             g_0_xyyz_0_xxxxxzzz_0, \
                             g_0_xyyz_0_xxxxyyyy_0, \
                             g_0_xyyz_0_xxxxyyyz_0, \
                             g_0_xyyz_0_xxxxyyzz_0, \
                             g_0_xyyz_0_xxxxyzzz_0, \
                             g_0_xyyz_0_xxxxzzzz_0, \
                             g_0_xyyz_0_xxxyyyyy_0, \
                             g_0_xyyz_0_xxxyyyyz_0, \
                             g_0_xyyz_0_xxxyyyzz_0, \
                             g_0_xyyz_0_xxxyyzzz_0, \
                             g_0_xyyz_0_xxxyzzzz_0, \
                             g_0_xyyz_0_xxxzzzzz_0, \
                             g_0_xyyz_0_xxyyyyyy_0, \
                             g_0_xyyz_0_xxyyyyyz_0, \
                             g_0_xyyz_0_xxyyyyzz_0, \
                             g_0_xyyz_0_xxyyyzzz_0, \
                             g_0_xyyz_0_xxyyzzzz_0, \
                             g_0_xyyz_0_xxyzzzzz_0, \
                             g_0_xyyz_0_xxzzzzzz_0, \
                             g_0_xyyz_0_xyyyyyyy_0, \
                             g_0_xyyz_0_xyyyyyyz_0, \
                             g_0_xyyz_0_xyyyyyzz_0, \
                             g_0_xyyz_0_xyyyyzzz_0, \
                             g_0_xyyz_0_xyyyzzzz_0, \
                             g_0_xyyz_0_xyyzzzzz_0, \
                             g_0_xyyz_0_xyzzzzzz_0, \
                             g_0_xyyz_0_xzzzzzzz_0, \
                             g_0_xyyz_0_yyyyyyyy_0, \
                             g_0_xyyz_0_yyyyyyyz_0, \
                             g_0_xyyz_0_yyyyyyzz_0, \
                             g_0_xyyz_0_yyyyyzzz_0, \
                             g_0_xyyz_0_yyyyzzzz_0, \
                             g_0_xyyz_0_yyyzzzzz_0, \
                             g_0_xyyz_0_yyzzzzzz_0, \
                             g_0_xyyz_0_yzzzzzzz_0, \
                             g_0_xyyz_0_zzzzzzzz_0, \
                             g_0_yyz_0_xxxxxxxz_0,  \
                             g_0_yyz_0_xxxxxxxz_1,  \
                             g_0_yyz_0_xxxxxxyz_0,  \
                             g_0_yyz_0_xxxxxxyz_1,  \
                             g_0_yyz_0_xxxxxxz_1,   \
                             g_0_yyz_0_xxxxxxzz_0,  \
                             g_0_yyz_0_xxxxxxzz_1,  \
                             g_0_yyz_0_xxxxxyyz_0,  \
                             g_0_yyz_0_xxxxxyyz_1,  \
                             g_0_yyz_0_xxxxxyz_1,   \
                             g_0_yyz_0_xxxxxyzz_0,  \
                             g_0_yyz_0_xxxxxyzz_1,  \
                             g_0_yyz_0_xxxxxzz_1,   \
                             g_0_yyz_0_xxxxxzzz_0,  \
                             g_0_yyz_0_xxxxxzzz_1,  \
                             g_0_yyz_0_xxxxyyyz_0,  \
                             g_0_yyz_0_xxxxyyyz_1,  \
                             g_0_yyz_0_xxxxyyz_1,   \
                             g_0_yyz_0_xxxxyyzz_0,  \
                             g_0_yyz_0_xxxxyyzz_1,  \
                             g_0_yyz_0_xxxxyzz_1,   \
                             g_0_yyz_0_xxxxyzzz_0,  \
                             g_0_yyz_0_xxxxyzzz_1,  \
                             g_0_yyz_0_xxxxzzz_1,   \
                             g_0_yyz_0_xxxxzzzz_0,  \
                             g_0_yyz_0_xxxxzzzz_1,  \
                             g_0_yyz_0_xxxyyyyz_0,  \
                             g_0_yyz_0_xxxyyyyz_1,  \
                             g_0_yyz_0_xxxyyyz_1,   \
                             g_0_yyz_0_xxxyyyzz_0,  \
                             g_0_yyz_0_xxxyyyzz_1,  \
                             g_0_yyz_0_xxxyyzz_1,   \
                             g_0_yyz_0_xxxyyzzz_0,  \
                             g_0_yyz_0_xxxyyzzz_1,  \
                             g_0_yyz_0_xxxyzzz_1,   \
                             g_0_yyz_0_xxxyzzzz_0,  \
                             g_0_yyz_0_xxxyzzzz_1,  \
                             g_0_yyz_0_xxxzzzz_1,   \
                             g_0_yyz_0_xxxzzzzz_0,  \
                             g_0_yyz_0_xxxzzzzz_1,  \
                             g_0_yyz_0_xxyyyyyz_0,  \
                             g_0_yyz_0_xxyyyyyz_1,  \
                             g_0_yyz_0_xxyyyyz_1,   \
                             g_0_yyz_0_xxyyyyzz_0,  \
                             g_0_yyz_0_xxyyyyzz_1,  \
                             g_0_yyz_0_xxyyyzz_1,   \
                             g_0_yyz_0_xxyyyzzz_0,  \
                             g_0_yyz_0_xxyyyzzz_1,  \
                             g_0_yyz_0_xxyyzzz_1,   \
                             g_0_yyz_0_xxyyzzzz_0,  \
                             g_0_yyz_0_xxyyzzzz_1,  \
                             g_0_yyz_0_xxyzzzz_1,   \
                             g_0_yyz_0_xxyzzzzz_0,  \
                             g_0_yyz_0_xxyzzzzz_1,  \
                             g_0_yyz_0_xxzzzzz_1,   \
                             g_0_yyz_0_xxzzzzzz_0,  \
                             g_0_yyz_0_xxzzzzzz_1,  \
                             g_0_yyz_0_xyyyyyyz_0,  \
                             g_0_yyz_0_xyyyyyyz_1,  \
                             g_0_yyz_0_xyyyyyz_1,   \
                             g_0_yyz_0_xyyyyyzz_0,  \
                             g_0_yyz_0_xyyyyyzz_1,  \
                             g_0_yyz_0_xyyyyzz_1,   \
                             g_0_yyz_0_xyyyyzzz_0,  \
                             g_0_yyz_0_xyyyyzzz_1,  \
                             g_0_yyz_0_xyyyzzz_1,   \
                             g_0_yyz_0_xyyyzzzz_0,  \
                             g_0_yyz_0_xyyyzzzz_1,  \
                             g_0_yyz_0_xyyzzzz_1,   \
                             g_0_yyz_0_xyyzzzzz_0,  \
                             g_0_yyz_0_xyyzzzzz_1,  \
                             g_0_yyz_0_xyzzzzz_1,   \
                             g_0_yyz_0_xyzzzzzz_0,  \
                             g_0_yyz_0_xyzzzzzz_1,  \
                             g_0_yyz_0_xzzzzzz_1,   \
                             g_0_yyz_0_xzzzzzzz_0,  \
                             g_0_yyz_0_xzzzzzzz_1,  \
                             g_0_yyz_0_yyyyyyyy_0,  \
                             g_0_yyz_0_yyyyyyyy_1,  \
                             g_0_yyz_0_yyyyyyyz_0,  \
                             g_0_yyz_0_yyyyyyyz_1,  \
                             g_0_yyz_0_yyyyyyz_1,   \
                             g_0_yyz_0_yyyyyyzz_0,  \
                             g_0_yyz_0_yyyyyyzz_1,  \
                             g_0_yyz_0_yyyyyzz_1,   \
                             g_0_yyz_0_yyyyyzzz_0,  \
                             g_0_yyz_0_yyyyyzzz_1,  \
                             g_0_yyz_0_yyyyzzz_1,   \
                             g_0_yyz_0_yyyyzzzz_0,  \
                             g_0_yyz_0_yyyyzzzz_1,  \
                             g_0_yyz_0_yyyzzzz_1,   \
                             g_0_yyz_0_yyyzzzzz_0,  \
                             g_0_yyz_0_yyyzzzzz_1,  \
                             g_0_yyz_0_yyzzzzz_1,   \
                             g_0_yyz_0_yyzzzzzz_0,  \
                             g_0_yyz_0_yyzzzzzz_1,  \
                             g_0_yyz_0_yzzzzzz_1,   \
                             g_0_yyz_0_yzzzzzzz_0,  \
                             g_0_yyz_0_yzzzzzzz_1,  \
                             g_0_yyz_0_zzzzzzz_1,   \
                             g_0_yyz_0_zzzzzzzz_0,  \
                             g_0_yyz_0_zzzzzzzz_1,  \
                             wp_x,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyz_0_xxxxxxxx_0[i] = g_0_xyy_0_xxxxxxxx_0[i] * pb_z + g_0_xyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxxxy_0[i] = g_0_xyy_0_xxxxxxxy_0[i] * pb_z + g_0_xyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxxxz_0[i] = 7.0 * g_0_yyz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxxxz_0[i] * pb_x + g_0_yyz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxxxyy_0[i] = g_0_xyy_0_xxxxxxyy_0[i] * pb_z + g_0_xyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxxyz_0[i] = 6.0 * g_0_yyz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxxyz_0[i] * pb_x + g_0_yyz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxxxzz_0[i] = 6.0 * g_0_yyz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxxzz_0[i] * pb_x + g_0_yyz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxxyyy_0[i] = g_0_xyy_0_xxxxxyyy_0[i] * pb_z + g_0_xyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxyyz_0[i] = 5.0 * g_0_yyz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxyyz_0[i] * pb_x + g_0_yyz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxxyzz_0[i] = 5.0 * g_0_yyz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxyzz_0[i] * pb_x + g_0_yyz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxxzzz_0[i] = 5.0 * g_0_yyz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxzzz_0[i] * pb_x + g_0_yyz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxyyyy_0[i] = g_0_xyy_0_xxxxyyyy_0[i] * pb_z + g_0_xyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxyyyz_0[i] = 4.0 * g_0_yyz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxyyyz_0[i] * pb_x + g_0_yyz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxyyzz_0[i] = 4.0 * g_0_yyz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxyyzz_0[i] * pb_x + g_0_yyz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxyzzz_0[i] = 4.0 * g_0_yyz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxyzzz_0[i] * pb_x + g_0_yyz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxzzzz_0[i] = 4.0 * g_0_yyz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxzzzz_0[i] * pb_x + g_0_yyz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyyyyy_0[i] = g_0_xyy_0_xxxyyyyy_0[i] * pb_z + g_0_xyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxyyyyz_0[i] = 3.0 * g_0_yyz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyyyyz_0[i] * pb_x + g_0_yyz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyyyzz_0[i] = 3.0 * g_0_yyz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyyyzz_0[i] * pb_x + g_0_yyz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyyzzz_0[i] = 3.0 * g_0_yyz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyyzzz_0[i] * pb_x + g_0_yyz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyzzzz_0[i] = 3.0 * g_0_yyz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyzzzz_0[i] * pb_x + g_0_yyz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxzzzzz_0[i] = 3.0 * g_0_yyz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxzzzzz_0[i] * pb_x + g_0_yyz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyyyyy_0[i] = g_0_xyy_0_xxyyyyyy_0[i] * pb_z + g_0_xyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxyyyyyz_0[i] = 2.0 * g_0_yyz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyyyyz_0[i] * pb_x + g_0_yyz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyyyzz_0[i] = 2.0 * g_0_yyz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyyyzz_0[i] * pb_x + g_0_yyz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyyzzz_0[i] = 2.0 * g_0_yyz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyyzzz_0[i] * pb_x + g_0_yyz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyzzzz_0[i] = 2.0 * g_0_yyz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyzzzz_0[i] * pb_x + g_0_yyz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyzzzzz_0[i] = 2.0 * g_0_yyz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyzzzzz_0[i] * pb_x + g_0_yyz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxzzzzzz_0[i] = 2.0 * g_0_yyz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxzzzzzz_0[i] * pb_x + g_0_yyz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyyyyy_0[i] = g_0_xyy_0_xyyyyyyy_0[i] * pb_z + g_0_xyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xyyyyyyz_0[i] = g_0_yyz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyyyyz_0[i] * pb_x + g_0_yyz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyyyzz_0[i] = g_0_yyz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyyyzz_0[i] * pb_x + g_0_yyz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyyzzz_0[i] = g_0_yyz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyyzzz_0[i] * pb_x + g_0_yyz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyzzzz_0[i] = g_0_yyz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyzzzz_0[i] * pb_x + g_0_yyz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyzzzzz_0[i] = g_0_yyz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyzzzzz_0[i] * pb_x + g_0_yyz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyzzzzzz_0[i] = g_0_yyz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyzzzzzz_0[i] * pb_x + g_0_yyz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xzzzzzzz_0[i] = g_0_yyz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xzzzzzzz_0[i] * pb_x + g_0_yyz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyyyy_0[i] = g_0_yyz_0_yyyyyyyy_0[i] * pb_x + g_0_yyz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyyyz_0[i] = g_0_yyz_0_yyyyyyyz_0[i] * pb_x + g_0_yyz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyyzz_0[i] = g_0_yyz_0_yyyyyyzz_0[i] * pb_x + g_0_yyz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyzzz_0[i] = g_0_yyz_0_yyyyyzzz_0[i] * pb_x + g_0_yyz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyzzzz_0[i] = g_0_yyz_0_yyyyzzzz_0[i] * pb_x + g_0_yyz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyzzzzz_0[i] = g_0_yyz_0_yyyzzzzz_0[i] * pb_x + g_0_yyz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyzzzzzz_0[i] = g_0_yyz_0_yyzzzzzz_0[i] * pb_x + g_0_yyz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yzzzzzzz_0[i] = g_0_yyz_0_yzzzzzzz_0[i] * pb_x + g_0_yyz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_zzzzzzzz_0[i] = g_0_yyz_0_zzzzzzzz_0[i] * pb_x + g_0_yyz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 360-405 components of targeted buffer : SGSL

    auto g_0_xyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 360);

    auto g_0_xyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 361);

    auto g_0_xyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 362);

    auto g_0_xyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 363);

    auto g_0_xyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 364);

    auto g_0_xyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 365);

    auto g_0_xyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 366);

    auto g_0_xyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 367);

    auto g_0_xyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 368);

    auto g_0_xyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 369);

    auto g_0_xyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 370);

    auto g_0_xyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 371);

    auto g_0_xyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 372);

    auto g_0_xyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 373);

    auto g_0_xyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 374);

    auto g_0_xyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 375);

    auto g_0_xyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 376);

    auto g_0_xyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 377);

    auto g_0_xyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 378);

    auto g_0_xyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 379);

    auto g_0_xyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 380);

    auto g_0_xyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 381);

    auto g_0_xyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 382);

    auto g_0_xyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 383);

    auto g_0_xyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 384);

    auto g_0_xyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 385);

    auto g_0_xyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 386);

    auto g_0_xyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 387);

    auto g_0_xyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 388);

    auto g_0_xyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 389);

    auto g_0_xyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 390);

    auto g_0_xyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 391);

    auto g_0_xyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 392);

    auto g_0_xyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 393);

    auto g_0_xyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 394);

    auto g_0_xyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 395);

    auto g_0_xyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 396);

    auto g_0_xyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 397);

    auto g_0_xyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 398);

    auto g_0_xyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 399);

    auto g_0_xyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 400);

    auto g_0_xyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 401);

    auto g_0_xyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 402);

    auto g_0_xyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 403);

    auto g_0_xyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 404);

#pragma omp simd aligned(g_0_xyzz_0_xxxxxxxx_0,     \
                             g_0_xyzz_0_xxxxxxxy_0, \
                             g_0_xyzz_0_xxxxxxxz_0, \
                             g_0_xyzz_0_xxxxxxyy_0, \
                             g_0_xyzz_0_xxxxxxyz_0, \
                             g_0_xyzz_0_xxxxxxzz_0, \
                             g_0_xyzz_0_xxxxxyyy_0, \
                             g_0_xyzz_0_xxxxxyyz_0, \
                             g_0_xyzz_0_xxxxxyzz_0, \
                             g_0_xyzz_0_xxxxxzzz_0, \
                             g_0_xyzz_0_xxxxyyyy_0, \
                             g_0_xyzz_0_xxxxyyyz_0, \
                             g_0_xyzz_0_xxxxyyzz_0, \
                             g_0_xyzz_0_xxxxyzzz_0, \
                             g_0_xyzz_0_xxxxzzzz_0, \
                             g_0_xyzz_0_xxxyyyyy_0, \
                             g_0_xyzz_0_xxxyyyyz_0, \
                             g_0_xyzz_0_xxxyyyzz_0, \
                             g_0_xyzz_0_xxxyyzzz_0, \
                             g_0_xyzz_0_xxxyzzzz_0, \
                             g_0_xyzz_0_xxxzzzzz_0, \
                             g_0_xyzz_0_xxyyyyyy_0, \
                             g_0_xyzz_0_xxyyyyyz_0, \
                             g_0_xyzz_0_xxyyyyzz_0, \
                             g_0_xyzz_0_xxyyyzzz_0, \
                             g_0_xyzz_0_xxyyzzzz_0, \
                             g_0_xyzz_0_xxyzzzzz_0, \
                             g_0_xyzz_0_xxzzzzzz_0, \
                             g_0_xyzz_0_xyyyyyyy_0, \
                             g_0_xyzz_0_xyyyyyyz_0, \
                             g_0_xyzz_0_xyyyyyzz_0, \
                             g_0_xyzz_0_xyyyyzzz_0, \
                             g_0_xyzz_0_xyyyzzzz_0, \
                             g_0_xyzz_0_xyyzzzzz_0, \
                             g_0_xyzz_0_xyzzzzzz_0, \
                             g_0_xyzz_0_xzzzzzzz_0, \
                             g_0_xyzz_0_yyyyyyyy_0, \
                             g_0_xyzz_0_yyyyyyyz_0, \
                             g_0_xyzz_0_yyyyyyzz_0, \
                             g_0_xyzz_0_yyyyyzzz_0, \
                             g_0_xyzz_0_yyyyzzzz_0, \
                             g_0_xyzz_0_yyyzzzzz_0, \
                             g_0_xyzz_0_yyzzzzzz_0, \
                             g_0_xyzz_0_yzzzzzzz_0, \
                             g_0_xyzz_0_zzzzzzzz_0, \
                             g_0_xzz_0_xxxxxxxx_0,  \
                             g_0_xzz_0_xxxxxxxx_1,  \
                             g_0_xzz_0_xxxxxxxz_0,  \
                             g_0_xzz_0_xxxxxxxz_1,  \
                             g_0_xzz_0_xxxxxxzz_0,  \
                             g_0_xzz_0_xxxxxxzz_1,  \
                             g_0_xzz_0_xxxxxzzz_0,  \
                             g_0_xzz_0_xxxxxzzz_1,  \
                             g_0_xzz_0_xxxxzzzz_0,  \
                             g_0_xzz_0_xxxxzzzz_1,  \
                             g_0_xzz_0_xxxzzzzz_0,  \
                             g_0_xzz_0_xxxzzzzz_1,  \
                             g_0_xzz_0_xxzzzzzz_0,  \
                             g_0_xzz_0_xxzzzzzz_1,  \
                             g_0_xzz_0_xzzzzzzz_0,  \
                             g_0_xzz_0_xzzzzzzz_1,  \
                             g_0_yzz_0_xxxxxxxy_0,  \
                             g_0_yzz_0_xxxxxxxy_1,  \
                             g_0_yzz_0_xxxxxxy_1,   \
                             g_0_yzz_0_xxxxxxyy_0,  \
                             g_0_yzz_0_xxxxxxyy_1,  \
                             g_0_yzz_0_xxxxxxyz_0,  \
                             g_0_yzz_0_xxxxxxyz_1,  \
                             g_0_yzz_0_xxxxxyy_1,   \
                             g_0_yzz_0_xxxxxyyy_0,  \
                             g_0_yzz_0_xxxxxyyy_1,  \
                             g_0_yzz_0_xxxxxyyz_0,  \
                             g_0_yzz_0_xxxxxyyz_1,  \
                             g_0_yzz_0_xxxxxyz_1,   \
                             g_0_yzz_0_xxxxxyzz_0,  \
                             g_0_yzz_0_xxxxxyzz_1,  \
                             g_0_yzz_0_xxxxyyy_1,   \
                             g_0_yzz_0_xxxxyyyy_0,  \
                             g_0_yzz_0_xxxxyyyy_1,  \
                             g_0_yzz_0_xxxxyyyz_0,  \
                             g_0_yzz_0_xxxxyyyz_1,  \
                             g_0_yzz_0_xxxxyyz_1,   \
                             g_0_yzz_0_xxxxyyzz_0,  \
                             g_0_yzz_0_xxxxyyzz_1,  \
                             g_0_yzz_0_xxxxyzz_1,   \
                             g_0_yzz_0_xxxxyzzz_0,  \
                             g_0_yzz_0_xxxxyzzz_1,  \
                             g_0_yzz_0_xxxyyyy_1,   \
                             g_0_yzz_0_xxxyyyyy_0,  \
                             g_0_yzz_0_xxxyyyyy_1,  \
                             g_0_yzz_0_xxxyyyyz_0,  \
                             g_0_yzz_0_xxxyyyyz_1,  \
                             g_0_yzz_0_xxxyyyz_1,   \
                             g_0_yzz_0_xxxyyyzz_0,  \
                             g_0_yzz_0_xxxyyyzz_1,  \
                             g_0_yzz_0_xxxyyzz_1,   \
                             g_0_yzz_0_xxxyyzzz_0,  \
                             g_0_yzz_0_xxxyyzzz_1,  \
                             g_0_yzz_0_xxxyzzz_1,   \
                             g_0_yzz_0_xxxyzzzz_0,  \
                             g_0_yzz_0_xxxyzzzz_1,  \
                             g_0_yzz_0_xxyyyyy_1,   \
                             g_0_yzz_0_xxyyyyyy_0,  \
                             g_0_yzz_0_xxyyyyyy_1,  \
                             g_0_yzz_0_xxyyyyyz_0,  \
                             g_0_yzz_0_xxyyyyyz_1,  \
                             g_0_yzz_0_xxyyyyz_1,   \
                             g_0_yzz_0_xxyyyyzz_0,  \
                             g_0_yzz_0_xxyyyyzz_1,  \
                             g_0_yzz_0_xxyyyzz_1,   \
                             g_0_yzz_0_xxyyyzzz_0,  \
                             g_0_yzz_0_xxyyyzzz_1,  \
                             g_0_yzz_0_xxyyzzz_1,   \
                             g_0_yzz_0_xxyyzzzz_0,  \
                             g_0_yzz_0_xxyyzzzz_1,  \
                             g_0_yzz_0_xxyzzzz_1,   \
                             g_0_yzz_0_xxyzzzzz_0,  \
                             g_0_yzz_0_xxyzzzzz_1,  \
                             g_0_yzz_0_xyyyyyy_1,   \
                             g_0_yzz_0_xyyyyyyy_0,  \
                             g_0_yzz_0_xyyyyyyy_1,  \
                             g_0_yzz_0_xyyyyyyz_0,  \
                             g_0_yzz_0_xyyyyyyz_1,  \
                             g_0_yzz_0_xyyyyyz_1,   \
                             g_0_yzz_0_xyyyyyzz_0,  \
                             g_0_yzz_0_xyyyyyzz_1,  \
                             g_0_yzz_0_xyyyyzz_1,   \
                             g_0_yzz_0_xyyyyzzz_0,  \
                             g_0_yzz_0_xyyyyzzz_1,  \
                             g_0_yzz_0_xyyyzzz_1,   \
                             g_0_yzz_0_xyyyzzzz_0,  \
                             g_0_yzz_0_xyyyzzzz_1,  \
                             g_0_yzz_0_xyyzzzz_1,   \
                             g_0_yzz_0_xyyzzzzz_0,  \
                             g_0_yzz_0_xyyzzzzz_1,  \
                             g_0_yzz_0_xyzzzzz_1,   \
                             g_0_yzz_0_xyzzzzzz_0,  \
                             g_0_yzz_0_xyzzzzzz_1,  \
                             g_0_yzz_0_yyyyyyy_1,   \
                             g_0_yzz_0_yyyyyyyy_0,  \
                             g_0_yzz_0_yyyyyyyy_1,  \
                             g_0_yzz_0_yyyyyyyz_0,  \
                             g_0_yzz_0_yyyyyyyz_1,  \
                             g_0_yzz_0_yyyyyyz_1,   \
                             g_0_yzz_0_yyyyyyzz_0,  \
                             g_0_yzz_0_yyyyyyzz_1,  \
                             g_0_yzz_0_yyyyyzz_1,   \
                             g_0_yzz_0_yyyyyzzz_0,  \
                             g_0_yzz_0_yyyyyzzz_1,  \
                             g_0_yzz_0_yyyyzzz_1,   \
                             g_0_yzz_0_yyyyzzzz_0,  \
                             g_0_yzz_0_yyyyzzzz_1,  \
                             g_0_yzz_0_yyyzzzz_1,   \
                             g_0_yzz_0_yyyzzzzz_0,  \
                             g_0_yzz_0_yyyzzzzz_1,  \
                             g_0_yzz_0_yyzzzzz_1,   \
                             g_0_yzz_0_yyzzzzzz_0,  \
                             g_0_yzz_0_yyzzzzzz_1,  \
                             g_0_yzz_0_yzzzzzz_1,   \
                             g_0_yzz_0_yzzzzzzz_0,  \
                             g_0_yzz_0_yzzzzzzz_1,  \
                             g_0_yzz_0_zzzzzzzz_0,  \
                             g_0_yzz_0_zzzzzzzz_1,  \
                             wp_x,                  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzz_0_xxxxxxxx_0[i] = g_0_xzz_0_xxxxxxxx_0[i] * pb_y + g_0_xzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxxxxy_0[i] = 7.0 * g_0_yzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxxxy_0[i] * pb_x + g_0_yzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxxxz_0[i] = g_0_xzz_0_xxxxxxxz_0[i] * pb_y + g_0_xzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxxxyy_0[i] = 6.0 * g_0_yzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxxyy_0[i] * pb_x + g_0_yzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxxyz_0[i] = 6.0 * g_0_yzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxxyz_0[i] * pb_x + g_0_yzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxxzz_0[i] = g_0_xzz_0_xxxxxxzz_0[i] * pb_y + g_0_xzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxxyyy_0[i] = 5.0 * g_0_yzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxyyy_0[i] * pb_x + g_0_yzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxyyz_0[i] = 5.0 * g_0_yzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxyyz_0[i] * pb_x + g_0_yzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxyzz_0[i] = 5.0 * g_0_yzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxyzz_0[i] * pb_x + g_0_yzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxzzz_0[i] = g_0_xzz_0_xxxxxzzz_0[i] * pb_y + g_0_xzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxyyyy_0[i] = 4.0 * g_0_yzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyyyy_0[i] * pb_x + g_0_yzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxyyyz_0[i] = 4.0 * g_0_yzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyyyz_0[i] * pb_x + g_0_yzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxyyzz_0[i] = 4.0 * g_0_yzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyyzz_0[i] * pb_x + g_0_yzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxyzzz_0[i] = 4.0 * g_0_yzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyzzz_0[i] * pb_x + g_0_yzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxzzzz_0[i] = g_0_xzz_0_xxxxzzzz_0[i] * pb_y + g_0_xzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxyyyyy_0[i] = 3.0 * g_0_yzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyyyy_0[i] * pb_x + g_0_yzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyyyyz_0[i] = 3.0 * g_0_yzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyyyz_0[i] * pb_x + g_0_yzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyyyzz_0[i] = 3.0 * g_0_yzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyyzz_0[i] * pb_x + g_0_yzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyyzzz_0[i] = 3.0 * g_0_yzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyzzz_0[i] * pb_x + g_0_yzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyzzzz_0[i] = 3.0 * g_0_yzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyzzzz_0[i] * pb_x + g_0_yzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxzzzzz_0[i] = g_0_xzz_0_xxxzzzzz_0[i] * pb_y + g_0_xzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxyyyyyy_0[i] = 2.0 * g_0_yzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyyyy_0[i] * pb_x + g_0_yzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyyyyz_0[i] = 2.0 * g_0_yzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyyyz_0[i] * pb_x + g_0_yzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyyyzz_0[i] = 2.0 * g_0_yzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyyzz_0[i] * pb_x + g_0_yzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyyzzz_0[i] = 2.0 * g_0_yzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyzzz_0[i] * pb_x + g_0_yzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyzzzz_0[i] = 2.0 * g_0_yzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyzzzz_0[i] * pb_x + g_0_yzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyzzzzz_0[i] = 2.0 * g_0_yzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyzzzzz_0[i] * pb_x + g_0_yzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxzzzzzz_0[i] = g_0_xzz_0_xxzzzzzz_0[i] * pb_y + g_0_xzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xyyyyyyy_0[i] = g_0_yzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyyyy_0[i] * pb_x + g_0_yzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyyyyz_0[i] = g_0_yzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyyyz_0[i] * pb_x + g_0_yzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyyyzz_0[i] = g_0_yzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyyzz_0[i] * pb_x + g_0_yzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyyzzz_0[i] = g_0_yzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyzzz_0[i] * pb_x + g_0_yzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyzzzz_0[i] = g_0_yzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyzzzz_0[i] * pb_x + g_0_yzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyzzzzz_0[i] = g_0_yzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyzzzzz_0[i] * pb_x + g_0_yzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyzzzzzz_0[i] = g_0_yzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzzzzzz_0[i] * pb_x + g_0_yzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xzzzzzzz_0[i] = g_0_xzz_0_xzzzzzzz_0[i] * pb_y + g_0_xzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_yyyyyyyy_0[i] = g_0_yzz_0_yyyyyyyy_0[i] * pb_x + g_0_yzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyyyyz_0[i] = g_0_yzz_0_yyyyyyyz_0[i] * pb_x + g_0_yzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyyyzz_0[i] = g_0_yzz_0_yyyyyyzz_0[i] * pb_x + g_0_yzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyyzzz_0[i] = g_0_yzz_0_yyyyyzzz_0[i] * pb_x + g_0_yzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyzzzz_0[i] = g_0_yzz_0_yyyyzzzz_0[i] * pb_x + g_0_yzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyzzzzz_0[i] = g_0_yzz_0_yyyzzzzz_0[i] * pb_x + g_0_yzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyzzzzzz_0[i] = g_0_yzz_0_yyzzzzzz_0[i] * pb_x + g_0_yzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yzzzzzzz_0[i] = g_0_yzz_0_yzzzzzzz_0[i] * pb_x + g_0_yzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_zzzzzzzz_0[i] = g_0_yzz_0_zzzzzzzz_0[i] * pb_x + g_0_yzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 405-450 components of targeted buffer : SGSL

    auto g_0_xzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 405);

    auto g_0_xzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 406);

    auto g_0_xzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 407);

    auto g_0_xzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 408);

    auto g_0_xzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 409);

    auto g_0_xzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 410);

    auto g_0_xzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 411);

    auto g_0_xzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 412);

    auto g_0_xzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 413);

    auto g_0_xzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 414);

    auto g_0_xzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 415);

    auto g_0_xzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 416);

    auto g_0_xzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 417);

    auto g_0_xzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 418);

    auto g_0_xzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 419);

    auto g_0_xzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 420);

    auto g_0_xzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 421);

    auto g_0_xzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 422);

    auto g_0_xzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 423);

    auto g_0_xzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 424);

    auto g_0_xzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 425);

    auto g_0_xzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 426);

    auto g_0_xzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 427);

    auto g_0_xzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 428);

    auto g_0_xzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 429);

    auto g_0_xzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 430);

    auto g_0_xzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 431);

    auto g_0_xzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 432);

    auto g_0_xzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 433);

    auto g_0_xzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 434);

    auto g_0_xzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 435);

    auto g_0_xzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 436);

    auto g_0_xzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 437);

    auto g_0_xzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 438);

    auto g_0_xzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 439);

    auto g_0_xzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 440);

    auto g_0_xzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 441);

    auto g_0_xzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 442);

    auto g_0_xzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 443);

    auto g_0_xzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 444);

    auto g_0_xzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 445);

    auto g_0_xzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 446);

    auto g_0_xzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 447);

    auto g_0_xzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 448);

    auto g_0_xzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 449);

#pragma omp simd aligned(g_0_xzzz_0_xxxxxxxx_0,     \
                             g_0_xzzz_0_xxxxxxxy_0, \
                             g_0_xzzz_0_xxxxxxxz_0, \
                             g_0_xzzz_0_xxxxxxyy_0, \
                             g_0_xzzz_0_xxxxxxyz_0, \
                             g_0_xzzz_0_xxxxxxzz_0, \
                             g_0_xzzz_0_xxxxxyyy_0, \
                             g_0_xzzz_0_xxxxxyyz_0, \
                             g_0_xzzz_0_xxxxxyzz_0, \
                             g_0_xzzz_0_xxxxxzzz_0, \
                             g_0_xzzz_0_xxxxyyyy_0, \
                             g_0_xzzz_0_xxxxyyyz_0, \
                             g_0_xzzz_0_xxxxyyzz_0, \
                             g_0_xzzz_0_xxxxyzzz_0, \
                             g_0_xzzz_0_xxxxzzzz_0, \
                             g_0_xzzz_0_xxxyyyyy_0, \
                             g_0_xzzz_0_xxxyyyyz_0, \
                             g_0_xzzz_0_xxxyyyzz_0, \
                             g_0_xzzz_0_xxxyyzzz_0, \
                             g_0_xzzz_0_xxxyzzzz_0, \
                             g_0_xzzz_0_xxxzzzzz_0, \
                             g_0_xzzz_0_xxyyyyyy_0, \
                             g_0_xzzz_0_xxyyyyyz_0, \
                             g_0_xzzz_0_xxyyyyzz_0, \
                             g_0_xzzz_0_xxyyyzzz_0, \
                             g_0_xzzz_0_xxyyzzzz_0, \
                             g_0_xzzz_0_xxyzzzzz_0, \
                             g_0_xzzz_0_xxzzzzzz_0, \
                             g_0_xzzz_0_xyyyyyyy_0, \
                             g_0_xzzz_0_xyyyyyyz_0, \
                             g_0_xzzz_0_xyyyyyzz_0, \
                             g_0_xzzz_0_xyyyyzzz_0, \
                             g_0_xzzz_0_xyyyzzzz_0, \
                             g_0_xzzz_0_xyyzzzzz_0, \
                             g_0_xzzz_0_xyzzzzzz_0, \
                             g_0_xzzz_0_xzzzzzzz_0, \
                             g_0_xzzz_0_yyyyyyyy_0, \
                             g_0_xzzz_0_yyyyyyyz_0, \
                             g_0_xzzz_0_yyyyyyzz_0, \
                             g_0_xzzz_0_yyyyyzzz_0, \
                             g_0_xzzz_0_yyyyzzzz_0, \
                             g_0_xzzz_0_yyyzzzzz_0, \
                             g_0_xzzz_0_yyzzzzzz_0, \
                             g_0_xzzz_0_yzzzzzzz_0, \
                             g_0_xzzz_0_zzzzzzzz_0, \
                             g_0_zzz_0_xxxxxxx_1,   \
                             g_0_zzz_0_xxxxxxxx_0,  \
                             g_0_zzz_0_xxxxxxxx_1,  \
                             g_0_zzz_0_xxxxxxxy_0,  \
                             g_0_zzz_0_xxxxxxxy_1,  \
                             g_0_zzz_0_xxxxxxxz_0,  \
                             g_0_zzz_0_xxxxxxxz_1,  \
                             g_0_zzz_0_xxxxxxy_1,   \
                             g_0_zzz_0_xxxxxxyy_0,  \
                             g_0_zzz_0_xxxxxxyy_1,  \
                             g_0_zzz_0_xxxxxxyz_0,  \
                             g_0_zzz_0_xxxxxxyz_1,  \
                             g_0_zzz_0_xxxxxxz_1,   \
                             g_0_zzz_0_xxxxxxzz_0,  \
                             g_0_zzz_0_xxxxxxzz_1,  \
                             g_0_zzz_0_xxxxxyy_1,   \
                             g_0_zzz_0_xxxxxyyy_0,  \
                             g_0_zzz_0_xxxxxyyy_1,  \
                             g_0_zzz_0_xxxxxyyz_0,  \
                             g_0_zzz_0_xxxxxyyz_1,  \
                             g_0_zzz_0_xxxxxyz_1,   \
                             g_0_zzz_0_xxxxxyzz_0,  \
                             g_0_zzz_0_xxxxxyzz_1,  \
                             g_0_zzz_0_xxxxxzz_1,   \
                             g_0_zzz_0_xxxxxzzz_0,  \
                             g_0_zzz_0_xxxxxzzz_1,  \
                             g_0_zzz_0_xxxxyyy_1,   \
                             g_0_zzz_0_xxxxyyyy_0,  \
                             g_0_zzz_0_xxxxyyyy_1,  \
                             g_0_zzz_0_xxxxyyyz_0,  \
                             g_0_zzz_0_xxxxyyyz_1,  \
                             g_0_zzz_0_xxxxyyz_1,   \
                             g_0_zzz_0_xxxxyyzz_0,  \
                             g_0_zzz_0_xxxxyyzz_1,  \
                             g_0_zzz_0_xxxxyzz_1,   \
                             g_0_zzz_0_xxxxyzzz_0,  \
                             g_0_zzz_0_xxxxyzzz_1,  \
                             g_0_zzz_0_xxxxzzz_1,   \
                             g_0_zzz_0_xxxxzzzz_0,  \
                             g_0_zzz_0_xxxxzzzz_1,  \
                             g_0_zzz_0_xxxyyyy_1,   \
                             g_0_zzz_0_xxxyyyyy_0,  \
                             g_0_zzz_0_xxxyyyyy_1,  \
                             g_0_zzz_0_xxxyyyyz_0,  \
                             g_0_zzz_0_xxxyyyyz_1,  \
                             g_0_zzz_0_xxxyyyz_1,   \
                             g_0_zzz_0_xxxyyyzz_0,  \
                             g_0_zzz_0_xxxyyyzz_1,  \
                             g_0_zzz_0_xxxyyzz_1,   \
                             g_0_zzz_0_xxxyyzzz_0,  \
                             g_0_zzz_0_xxxyyzzz_1,  \
                             g_0_zzz_0_xxxyzzz_1,   \
                             g_0_zzz_0_xxxyzzzz_0,  \
                             g_0_zzz_0_xxxyzzzz_1,  \
                             g_0_zzz_0_xxxzzzz_1,   \
                             g_0_zzz_0_xxxzzzzz_0,  \
                             g_0_zzz_0_xxxzzzzz_1,  \
                             g_0_zzz_0_xxyyyyy_1,   \
                             g_0_zzz_0_xxyyyyyy_0,  \
                             g_0_zzz_0_xxyyyyyy_1,  \
                             g_0_zzz_0_xxyyyyyz_0,  \
                             g_0_zzz_0_xxyyyyyz_1,  \
                             g_0_zzz_0_xxyyyyz_1,   \
                             g_0_zzz_0_xxyyyyzz_0,  \
                             g_0_zzz_0_xxyyyyzz_1,  \
                             g_0_zzz_0_xxyyyzz_1,   \
                             g_0_zzz_0_xxyyyzzz_0,  \
                             g_0_zzz_0_xxyyyzzz_1,  \
                             g_0_zzz_0_xxyyzzz_1,   \
                             g_0_zzz_0_xxyyzzzz_0,  \
                             g_0_zzz_0_xxyyzzzz_1,  \
                             g_0_zzz_0_xxyzzzz_1,   \
                             g_0_zzz_0_xxyzzzzz_0,  \
                             g_0_zzz_0_xxyzzzzz_1,  \
                             g_0_zzz_0_xxzzzzz_1,   \
                             g_0_zzz_0_xxzzzzzz_0,  \
                             g_0_zzz_0_xxzzzzzz_1,  \
                             g_0_zzz_0_xyyyyyy_1,   \
                             g_0_zzz_0_xyyyyyyy_0,  \
                             g_0_zzz_0_xyyyyyyy_1,  \
                             g_0_zzz_0_xyyyyyyz_0,  \
                             g_0_zzz_0_xyyyyyyz_1,  \
                             g_0_zzz_0_xyyyyyz_1,   \
                             g_0_zzz_0_xyyyyyzz_0,  \
                             g_0_zzz_0_xyyyyyzz_1,  \
                             g_0_zzz_0_xyyyyzz_1,   \
                             g_0_zzz_0_xyyyyzzz_0,  \
                             g_0_zzz_0_xyyyyzzz_1,  \
                             g_0_zzz_0_xyyyzzz_1,   \
                             g_0_zzz_0_xyyyzzzz_0,  \
                             g_0_zzz_0_xyyyzzzz_1,  \
                             g_0_zzz_0_xyyzzzz_1,   \
                             g_0_zzz_0_xyyzzzzz_0,  \
                             g_0_zzz_0_xyyzzzzz_1,  \
                             g_0_zzz_0_xyzzzzz_1,   \
                             g_0_zzz_0_xyzzzzzz_0,  \
                             g_0_zzz_0_xyzzzzzz_1,  \
                             g_0_zzz_0_xzzzzzz_1,   \
                             g_0_zzz_0_xzzzzzzz_0,  \
                             g_0_zzz_0_xzzzzzzz_1,  \
                             g_0_zzz_0_yyyyyyy_1,   \
                             g_0_zzz_0_yyyyyyyy_0,  \
                             g_0_zzz_0_yyyyyyyy_1,  \
                             g_0_zzz_0_yyyyyyyz_0,  \
                             g_0_zzz_0_yyyyyyyz_1,  \
                             g_0_zzz_0_yyyyyyz_1,   \
                             g_0_zzz_0_yyyyyyzz_0,  \
                             g_0_zzz_0_yyyyyyzz_1,  \
                             g_0_zzz_0_yyyyyzz_1,   \
                             g_0_zzz_0_yyyyyzzz_0,  \
                             g_0_zzz_0_yyyyyzzz_1,  \
                             g_0_zzz_0_yyyyzzz_1,   \
                             g_0_zzz_0_yyyyzzzz_0,  \
                             g_0_zzz_0_yyyyzzzz_1,  \
                             g_0_zzz_0_yyyzzzz_1,   \
                             g_0_zzz_0_yyyzzzzz_0,  \
                             g_0_zzz_0_yyyzzzzz_1,  \
                             g_0_zzz_0_yyzzzzz_1,   \
                             g_0_zzz_0_yyzzzzzz_0,  \
                             g_0_zzz_0_yyzzzzzz_1,  \
                             g_0_zzz_0_yzzzzzz_1,   \
                             g_0_zzz_0_yzzzzzzz_0,  \
                             g_0_zzz_0_yzzzzzzz_1,  \
                             g_0_zzz_0_zzzzzzz_1,   \
                             g_0_zzz_0_zzzzzzzz_0,  \
                             g_0_zzz_0_zzzzzzzz_1,  \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzz_0_xxxxxxxx_0[i] = 8.0 * g_0_zzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxxx_0[i] * pb_x + g_0_zzz_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxxxy_0[i] = 7.0 * g_0_zzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxxy_0[i] * pb_x + g_0_zzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxxxz_0[i] = 7.0 * g_0_zzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxxz_0[i] * pb_x + g_0_zzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxxyy_0[i] = 6.0 * g_0_zzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxyy_0[i] * pb_x + g_0_zzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxxyz_0[i] = 6.0 * g_0_zzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxyz_0[i] * pb_x + g_0_zzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxxzz_0[i] = 6.0 * g_0_zzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxzz_0[i] * pb_x + g_0_zzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxyyy_0[i] = 5.0 * g_0_zzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyyy_0[i] * pb_x + g_0_zzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxyyz_0[i] = 5.0 * g_0_zzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyyz_0[i] * pb_x + g_0_zzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxyzz_0[i] = 5.0 * g_0_zzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyzz_0[i] * pb_x + g_0_zzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxzzz_0[i] = 5.0 * g_0_zzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxzzz_0[i] * pb_x + g_0_zzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyyyy_0[i] = 4.0 * g_0_zzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyyy_0[i] * pb_x + g_0_zzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyyyz_0[i] = 4.0 * g_0_zzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyyz_0[i] * pb_x + g_0_zzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyyzz_0[i] = 4.0 * g_0_zzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyzz_0[i] * pb_x + g_0_zzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyzzz_0[i] = 4.0 * g_0_zzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyzzz_0[i] * pb_x + g_0_zzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxzzzz_0[i] = 4.0 * g_0_zzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxzzzz_0[i] * pb_x + g_0_zzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyyyy_0[i] = 3.0 * g_0_zzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyyy_0[i] * pb_x + g_0_zzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyyyz_0[i] = 3.0 * g_0_zzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyyz_0[i] * pb_x + g_0_zzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyyzz_0[i] = 3.0 * g_0_zzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyzz_0[i] * pb_x + g_0_zzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyzzz_0[i] = 3.0 * g_0_zzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyzzz_0[i] * pb_x + g_0_zzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyzzzz_0[i] = 3.0 * g_0_zzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzzzz_0[i] * pb_x + g_0_zzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxzzzzz_0[i] = 3.0 * g_0_zzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxzzzzz_0[i] * pb_x + g_0_zzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyyyy_0[i] = 2.0 * g_0_zzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyyy_0[i] * pb_x + g_0_zzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyyyz_0[i] = 2.0 * g_0_zzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyyz_0[i] * pb_x + g_0_zzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyyzz_0[i] = 2.0 * g_0_zzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyzz_0[i] * pb_x + g_0_zzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyzzz_0[i] = 2.0 * g_0_zzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyzzz_0[i] * pb_x + g_0_zzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyzzzz_0[i] = 2.0 * g_0_zzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzzzz_0[i] * pb_x + g_0_zzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyzzzzz_0[i] = 2.0 * g_0_zzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzzzz_0[i] * pb_x + g_0_zzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxzzzzzz_0[i] = 2.0 * g_0_zzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzzzzzz_0[i] * pb_x + g_0_zzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyyyy_0[i] = g_0_zzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyyy_0[i] * pb_x + g_0_zzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyyyz_0[i] = g_0_zzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyyz_0[i] * pb_x + g_0_zzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyyzz_0[i] = g_0_zzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyzz_0[i] * pb_x + g_0_zzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyzzz_0[i] = g_0_zzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyzzz_0[i] * pb_x + g_0_zzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyzzzz_0[i] = g_0_zzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzzzz_0[i] * pb_x + g_0_zzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyzzzzz_0[i] = g_0_zzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzzzz_0[i] * pb_x + g_0_zzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyzzzzzz_0[i] = g_0_zzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzzzz_0[i] * pb_x + g_0_zzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xzzzzzzz_0[i] = g_0_zzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzzzzzz_0[i] * pb_x + g_0_zzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyyyy_0[i] = g_0_zzz_0_yyyyyyyy_0[i] * pb_x + g_0_zzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyyyz_0[i] = g_0_zzz_0_yyyyyyyz_0[i] * pb_x + g_0_zzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyyzz_0[i] = g_0_zzz_0_yyyyyyzz_0[i] * pb_x + g_0_zzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyzzz_0[i] = g_0_zzz_0_yyyyyzzz_0[i] * pb_x + g_0_zzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyzzzz_0[i] = g_0_zzz_0_yyyyzzzz_0[i] * pb_x + g_0_zzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyzzzzz_0[i] = g_0_zzz_0_yyyzzzzz_0[i] * pb_x + g_0_zzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyzzzzzz_0[i] = g_0_zzz_0_yyzzzzzz_0[i] * pb_x + g_0_zzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yzzzzzzz_0[i] = g_0_zzz_0_yzzzzzzz_0[i] * pb_x + g_0_zzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_zzzzzzzz_0[i] = g_0_zzz_0_zzzzzzzz_0[i] * pb_x + g_0_zzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 450-495 components of targeted buffer : SGSL

    auto g_0_yyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 450);

    auto g_0_yyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 451);

    auto g_0_yyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 452);

    auto g_0_yyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 453);

    auto g_0_yyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 454);

    auto g_0_yyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 455);

    auto g_0_yyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 456);

    auto g_0_yyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 457);

    auto g_0_yyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 458);

    auto g_0_yyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 459);

    auto g_0_yyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 460);

    auto g_0_yyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 461);

    auto g_0_yyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 462);

    auto g_0_yyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 463);

    auto g_0_yyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 464);

    auto g_0_yyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 465);

    auto g_0_yyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 466);

    auto g_0_yyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 467);

    auto g_0_yyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 468);

    auto g_0_yyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 469);

    auto g_0_yyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 470);

    auto g_0_yyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 471);

    auto g_0_yyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 472);

    auto g_0_yyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 473);

    auto g_0_yyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 474);

    auto g_0_yyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 475);

    auto g_0_yyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 476);

    auto g_0_yyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 477);

    auto g_0_yyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 478);

    auto g_0_yyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 479);

    auto g_0_yyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 480);

    auto g_0_yyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 481);

    auto g_0_yyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 482);

    auto g_0_yyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 483);

    auto g_0_yyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 484);

    auto g_0_yyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 485);

    auto g_0_yyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 486);

    auto g_0_yyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 487);

    auto g_0_yyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 488);

    auto g_0_yyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 489);

    auto g_0_yyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 490);

    auto g_0_yyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 491);

    auto g_0_yyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 492);

    auto g_0_yyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 493);

    auto g_0_yyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 494);

#pragma omp simd aligned(g_0_yy_0_xxxxxxxx_0,       \
                             g_0_yy_0_xxxxxxxx_1,   \
                             g_0_yy_0_xxxxxxxy_0,   \
                             g_0_yy_0_xxxxxxxy_1,   \
                             g_0_yy_0_xxxxxxxz_0,   \
                             g_0_yy_0_xxxxxxxz_1,   \
                             g_0_yy_0_xxxxxxyy_0,   \
                             g_0_yy_0_xxxxxxyy_1,   \
                             g_0_yy_0_xxxxxxyz_0,   \
                             g_0_yy_0_xxxxxxyz_1,   \
                             g_0_yy_0_xxxxxxzz_0,   \
                             g_0_yy_0_xxxxxxzz_1,   \
                             g_0_yy_0_xxxxxyyy_0,   \
                             g_0_yy_0_xxxxxyyy_1,   \
                             g_0_yy_0_xxxxxyyz_0,   \
                             g_0_yy_0_xxxxxyyz_1,   \
                             g_0_yy_0_xxxxxyzz_0,   \
                             g_0_yy_0_xxxxxyzz_1,   \
                             g_0_yy_0_xxxxxzzz_0,   \
                             g_0_yy_0_xxxxxzzz_1,   \
                             g_0_yy_0_xxxxyyyy_0,   \
                             g_0_yy_0_xxxxyyyy_1,   \
                             g_0_yy_0_xxxxyyyz_0,   \
                             g_0_yy_0_xxxxyyyz_1,   \
                             g_0_yy_0_xxxxyyzz_0,   \
                             g_0_yy_0_xxxxyyzz_1,   \
                             g_0_yy_0_xxxxyzzz_0,   \
                             g_0_yy_0_xxxxyzzz_1,   \
                             g_0_yy_0_xxxxzzzz_0,   \
                             g_0_yy_0_xxxxzzzz_1,   \
                             g_0_yy_0_xxxyyyyy_0,   \
                             g_0_yy_0_xxxyyyyy_1,   \
                             g_0_yy_0_xxxyyyyz_0,   \
                             g_0_yy_0_xxxyyyyz_1,   \
                             g_0_yy_0_xxxyyyzz_0,   \
                             g_0_yy_0_xxxyyyzz_1,   \
                             g_0_yy_0_xxxyyzzz_0,   \
                             g_0_yy_0_xxxyyzzz_1,   \
                             g_0_yy_0_xxxyzzzz_0,   \
                             g_0_yy_0_xxxyzzzz_1,   \
                             g_0_yy_0_xxxzzzzz_0,   \
                             g_0_yy_0_xxxzzzzz_1,   \
                             g_0_yy_0_xxyyyyyy_0,   \
                             g_0_yy_0_xxyyyyyy_1,   \
                             g_0_yy_0_xxyyyyyz_0,   \
                             g_0_yy_0_xxyyyyyz_1,   \
                             g_0_yy_0_xxyyyyzz_0,   \
                             g_0_yy_0_xxyyyyzz_1,   \
                             g_0_yy_0_xxyyyzzz_0,   \
                             g_0_yy_0_xxyyyzzz_1,   \
                             g_0_yy_0_xxyyzzzz_0,   \
                             g_0_yy_0_xxyyzzzz_1,   \
                             g_0_yy_0_xxyzzzzz_0,   \
                             g_0_yy_0_xxyzzzzz_1,   \
                             g_0_yy_0_xxzzzzzz_0,   \
                             g_0_yy_0_xxzzzzzz_1,   \
                             g_0_yy_0_xyyyyyyy_0,   \
                             g_0_yy_0_xyyyyyyy_1,   \
                             g_0_yy_0_xyyyyyyz_0,   \
                             g_0_yy_0_xyyyyyyz_1,   \
                             g_0_yy_0_xyyyyyzz_0,   \
                             g_0_yy_0_xyyyyyzz_1,   \
                             g_0_yy_0_xyyyyzzz_0,   \
                             g_0_yy_0_xyyyyzzz_1,   \
                             g_0_yy_0_xyyyzzzz_0,   \
                             g_0_yy_0_xyyyzzzz_1,   \
                             g_0_yy_0_xyyzzzzz_0,   \
                             g_0_yy_0_xyyzzzzz_1,   \
                             g_0_yy_0_xyzzzzzz_0,   \
                             g_0_yy_0_xyzzzzzz_1,   \
                             g_0_yy_0_xzzzzzzz_0,   \
                             g_0_yy_0_xzzzzzzz_1,   \
                             g_0_yy_0_yyyyyyyy_0,   \
                             g_0_yy_0_yyyyyyyy_1,   \
                             g_0_yy_0_yyyyyyyz_0,   \
                             g_0_yy_0_yyyyyyyz_1,   \
                             g_0_yy_0_yyyyyyzz_0,   \
                             g_0_yy_0_yyyyyyzz_1,   \
                             g_0_yy_0_yyyyyzzz_0,   \
                             g_0_yy_0_yyyyyzzz_1,   \
                             g_0_yy_0_yyyyzzzz_0,   \
                             g_0_yy_0_yyyyzzzz_1,   \
                             g_0_yy_0_yyyzzzzz_0,   \
                             g_0_yy_0_yyyzzzzz_1,   \
                             g_0_yy_0_yyzzzzzz_0,   \
                             g_0_yy_0_yyzzzzzz_1,   \
                             g_0_yy_0_yzzzzzzz_0,   \
                             g_0_yy_0_yzzzzzzz_1,   \
                             g_0_yy_0_zzzzzzzz_0,   \
                             g_0_yy_0_zzzzzzzz_1,   \
                             g_0_yyy_0_xxxxxxx_1,   \
                             g_0_yyy_0_xxxxxxxx_0,  \
                             g_0_yyy_0_xxxxxxxx_1,  \
                             g_0_yyy_0_xxxxxxxy_0,  \
                             g_0_yyy_0_xxxxxxxy_1,  \
                             g_0_yyy_0_xxxxxxxz_0,  \
                             g_0_yyy_0_xxxxxxxz_1,  \
                             g_0_yyy_0_xxxxxxy_1,   \
                             g_0_yyy_0_xxxxxxyy_0,  \
                             g_0_yyy_0_xxxxxxyy_1,  \
                             g_0_yyy_0_xxxxxxyz_0,  \
                             g_0_yyy_0_xxxxxxyz_1,  \
                             g_0_yyy_0_xxxxxxz_1,   \
                             g_0_yyy_0_xxxxxxzz_0,  \
                             g_0_yyy_0_xxxxxxzz_1,  \
                             g_0_yyy_0_xxxxxyy_1,   \
                             g_0_yyy_0_xxxxxyyy_0,  \
                             g_0_yyy_0_xxxxxyyy_1,  \
                             g_0_yyy_0_xxxxxyyz_0,  \
                             g_0_yyy_0_xxxxxyyz_1,  \
                             g_0_yyy_0_xxxxxyz_1,   \
                             g_0_yyy_0_xxxxxyzz_0,  \
                             g_0_yyy_0_xxxxxyzz_1,  \
                             g_0_yyy_0_xxxxxzz_1,   \
                             g_0_yyy_0_xxxxxzzz_0,  \
                             g_0_yyy_0_xxxxxzzz_1,  \
                             g_0_yyy_0_xxxxyyy_1,   \
                             g_0_yyy_0_xxxxyyyy_0,  \
                             g_0_yyy_0_xxxxyyyy_1,  \
                             g_0_yyy_0_xxxxyyyz_0,  \
                             g_0_yyy_0_xxxxyyyz_1,  \
                             g_0_yyy_0_xxxxyyz_1,   \
                             g_0_yyy_0_xxxxyyzz_0,  \
                             g_0_yyy_0_xxxxyyzz_1,  \
                             g_0_yyy_0_xxxxyzz_1,   \
                             g_0_yyy_0_xxxxyzzz_0,  \
                             g_0_yyy_0_xxxxyzzz_1,  \
                             g_0_yyy_0_xxxxzzz_1,   \
                             g_0_yyy_0_xxxxzzzz_0,  \
                             g_0_yyy_0_xxxxzzzz_1,  \
                             g_0_yyy_0_xxxyyyy_1,   \
                             g_0_yyy_0_xxxyyyyy_0,  \
                             g_0_yyy_0_xxxyyyyy_1,  \
                             g_0_yyy_0_xxxyyyyz_0,  \
                             g_0_yyy_0_xxxyyyyz_1,  \
                             g_0_yyy_0_xxxyyyz_1,   \
                             g_0_yyy_0_xxxyyyzz_0,  \
                             g_0_yyy_0_xxxyyyzz_1,  \
                             g_0_yyy_0_xxxyyzz_1,   \
                             g_0_yyy_0_xxxyyzzz_0,  \
                             g_0_yyy_0_xxxyyzzz_1,  \
                             g_0_yyy_0_xxxyzzz_1,   \
                             g_0_yyy_0_xxxyzzzz_0,  \
                             g_0_yyy_0_xxxyzzzz_1,  \
                             g_0_yyy_0_xxxzzzz_1,   \
                             g_0_yyy_0_xxxzzzzz_0,  \
                             g_0_yyy_0_xxxzzzzz_1,  \
                             g_0_yyy_0_xxyyyyy_1,   \
                             g_0_yyy_0_xxyyyyyy_0,  \
                             g_0_yyy_0_xxyyyyyy_1,  \
                             g_0_yyy_0_xxyyyyyz_0,  \
                             g_0_yyy_0_xxyyyyyz_1,  \
                             g_0_yyy_0_xxyyyyz_1,   \
                             g_0_yyy_0_xxyyyyzz_0,  \
                             g_0_yyy_0_xxyyyyzz_1,  \
                             g_0_yyy_0_xxyyyzz_1,   \
                             g_0_yyy_0_xxyyyzzz_0,  \
                             g_0_yyy_0_xxyyyzzz_1,  \
                             g_0_yyy_0_xxyyzzz_1,   \
                             g_0_yyy_0_xxyyzzzz_0,  \
                             g_0_yyy_0_xxyyzzzz_1,  \
                             g_0_yyy_0_xxyzzzz_1,   \
                             g_0_yyy_0_xxyzzzzz_0,  \
                             g_0_yyy_0_xxyzzzzz_1,  \
                             g_0_yyy_0_xxzzzzz_1,   \
                             g_0_yyy_0_xxzzzzzz_0,  \
                             g_0_yyy_0_xxzzzzzz_1,  \
                             g_0_yyy_0_xyyyyyy_1,   \
                             g_0_yyy_0_xyyyyyyy_0,  \
                             g_0_yyy_0_xyyyyyyy_1,  \
                             g_0_yyy_0_xyyyyyyz_0,  \
                             g_0_yyy_0_xyyyyyyz_1,  \
                             g_0_yyy_0_xyyyyyz_1,   \
                             g_0_yyy_0_xyyyyyzz_0,  \
                             g_0_yyy_0_xyyyyyzz_1,  \
                             g_0_yyy_0_xyyyyzz_1,   \
                             g_0_yyy_0_xyyyyzzz_0,  \
                             g_0_yyy_0_xyyyyzzz_1,  \
                             g_0_yyy_0_xyyyzzz_1,   \
                             g_0_yyy_0_xyyyzzzz_0,  \
                             g_0_yyy_0_xyyyzzzz_1,  \
                             g_0_yyy_0_xyyzzzz_1,   \
                             g_0_yyy_0_xyyzzzzz_0,  \
                             g_0_yyy_0_xyyzzzzz_1,  \
                             g_0_yyy_0_xyzzzzz_1,   \
                             g_0_yyy_0_xyzzzzzz_0,  \
                             g_0_yyy_0_xyzzzzzz_1,  \
                             g_0_yyy_0_xzzzzzz_1,   \
                             g_0_yyy_0_xzzzzzzz_0,  \
                             g_0_yyy_0_xzzzzzzz_1,  \
                             g_0_yyy_0_yyyyyyy_1,   \
                             g_0_yyy_0_yyyyyyyy_0,  \
                             g_0_yyy_0_yyyyyyyy_1,  \
                             g_0_yyy_0_yyyyyyyz_0,  \
                             g_0_yyy_0_yyyyyyyz_1,  \
                             g_0_yyy_0_yyyyyyz_1,   \
                             g_0_yyy_0_yyyyyyzz_0,  \
                             g_0_yyy_0_yyyyyyzz_1,  \
                             g_0_yyy_0_yyyyyzz_1,   \
                             g_0_yyy_0_yyyyyzzz_0,  \
                             g_0_yyy_0_yyyyyzzz_1,  \
                             g_0_yyy_0_yyyyzzz_1,   \
                             g_0_yyy_0_yyyyzzzz_0,  \
                             g_0_yyy_0_yyyyzzzz_1,  \
                             g_0_yyy_0_yyyzzzz_1,   \
                             g_0_yyy_0_yyyzzzzz_0,  \
                             g_0_yyy_0_yyyzzzzz_1,  \
                             g_0_yyy_0_yyzzzzz_1,   \
                             g_0_yyy_0_yyzzzzzz_0,  \
                             g_0_yyy_0_yyzzzzzz_1,  \
                             g_0_yyy_0_yzzzzzz_1,   \
                             g_0_yyy_0_yzzzzzzz_0,  \
                             g_0_yyy_0_yzzzzzzz_1,  \
                             g_0_yyy_0_zzzzzzz_1,   \
                             g_0_yyy_0_zzzzzzzz_0,  \
                             g_0_yyy_0_zzzzzzzz_1,  \
                             g_0_yyyy_0_xxxxxxxx_0, \
                             g_0_yyyy_0_xxxxxxxy_0, \
                             g_0_yyyy_0_xxxxxxxz_0, \
                             g_0_yyyy_0_xxxxxxyy_0, \
                             g_0_yyyy_0_xxxxxxyz_0, \
                             g_0_yyyy_0_xxxxxxzz_0, \
                             g_0_yyyy_0_xxxxxyyy_0, \
                             g_0_yyyy_0_xxxxxyyz_0, \
                             g_0_yyyy_0_xxxxxyzz_0, \
                             g_0_yyyy_0_xxxxxzzz_0, \
                             g_0_yyyy_0_xxxxyyyy_0, \
                             g_0_yyyy_0_xxxxyyyz_0, \
                             g_0_yyyy_0_xxxxyyzz_0, \
                             g_0_yyyy_0_xxxxyzzz_0, \
                             g_0_yyyy_0_xxxxzzzz_0, \
                             g_0_yyyy_0_xxxyyyyy_0, \
                             g_0_yyyy_0_xxxyyyyz_0, \
                             g_0_yyyy_0_xxxyyyzz_0, \
                             g_0_yyyy_0_xxxyyzzz_0, \
                             g_0_yyyy_0_xxxyzzzz_0, \
                             g_0_yyyy_0_xxxzzzzz_0, \
                             g_0_yyyy_0_xxyyyyyy_0, \
                             g_0_yyyy_0_xxyyyyyz_0, \
                             g_0_yyyy_0_xxyyyyzz_0, \
                             g_0_yyyy_0_xxyyyzzz_0, \
                             g_0_yyyy_0_xxyyzzzz_0, \
                             g_0_yyyy_0_xxyzzzzz_0, \
                             g_0_yyyy_0_xxzzzzzz_0, \
                             g_0_yyyy_0_xyyyyyyy_0, \
                             g_0_yyyy_0_xyyyyyyz_0, \
                             g_0_yyyy_0_xyyyyyzz_0, \
                             g_0_yyyy_0_xyyyyzzz_0, \
                             g_0_yyyy_0_xyyyzzzz_0, \
                             g_0_yyyy_0_xyyzzzzz_0, \
                             g_0_yyyy_0_xyzzzzzz_0, \
                             g_0_yyyy_0_xzzzzzzz_0, \
                             g_0_yyyy_0_yyyyyyyy_0, \
                             g_0_yyyy_0_yyyyyyyz_0, \
                             g_0_yyyy_0_yyyyyyzz_0, \
                             g_0_yyyy_0_yyyyyzzz_0, \
                             g_0_yyyy_0_yyyyzzzz_0, \
                             g_0_yyyy_0_yyyzzzzz_0, \
                             g_0_yyyy_0_yyzzzzzz_0, \
                             g_0_yyyy_0_yzzzzzzz_0, \
                             g_0_yyyy_0_zzzzzzzz_0, \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyy_0_xxxxxxxx_0[i] = 3.0 * g_0_yy_0_xxxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxxxx_0[i] * pb_y +
                                   g_0_yyy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxxxy_0[i] = 3.0 * g_0_yy_0_xxxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxxy_1[i] * fti_ab_0 +
                                   g_0_yyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxxy_0[i] * pb_y + g_0_yyy_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxxxz_0[i] = 3.0 * g_0_yy_0_xxxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxxxz_0[i] * pb_y +
                                   g_0_yyy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxxyy_0[i] = 3.0 * g_0_yy_0_xxxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxyy_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxyy_0[i] * pb_y + g_0_yyy_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxxyz_0[i] = 3.0 * g_0_yy_0_xxxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxyz_1[i] * fti_ab_0 +
                                   g_0_yyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxyz_0[i] * pb_y + g_0_yyy_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxxzz_0[i] = 3.0 * g_0_yy_0_xxxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxxzz_0[i] * pb_y +
                                   g_0_yyy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxyyy_0[i] = 3.0 * g_0_yy_0_xxxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxyyy_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyyy_0[i] * pb_y + g_0_yyy_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxyyz_0[i] = 3.0 * g_0_yy_0_xxxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyyz_0[i] * pb_y + g_0_yyy_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxyzz_0[i] = 3.0 * g_0_yy_0_xxxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxyzz_1[i] * fti_ab_0 +
                                   g_0_yyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyzz_0[i] * pb_y + g_0_yyy_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxzzz_0[i] = 3.0 * g_0_yy_0_xxxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxzzz_0[i] * pb_y +
                                   g_0_yyy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyyyy_0[i] = 3.0 * g_0_yy_0_xxxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyyyy_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyyy_0[i] * pb_y + g_0_yyy_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyyyz_0[i] = 3.0 * g_0_yy_0_xxxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyyz_0[i] * pb_y + g_0_yyy_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyyzz_0[i] = 3.0 * g_0_yy_0_xxxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyzz_0[i] * pb_y + g_0_yyy_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyzzz_0[i] = 3.0 * g_0_yy_0_xxxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyzzz_1[i] * fti_ab_0 +
                                   g_0_yyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyzzz_0[i] * pb_y + g_0_yyy_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxzzzz_0[i] = 3.0 * g_0_yy_0_xxxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxzzzz_0[i] * pb_y +
                                   g_0_yyy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyyyy_0[i] = 3.0 * g_0_yy_0_xxxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyyyy_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyyy_0[i] * pb_y + g_0_yyy_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyyyz_0[i] = 3.0 * g_0_yy_0_xxxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyyyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyyz_0[i] * pb_y + g_0_yyy_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyyzz_0[i] = 3.0 * g_0_yy_0_xxxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyzz_0[i] * pb_y + g_0_yyy_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyzzz_0[i] = 3.0 * g_0_yy_0_xxxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyzzz_0[i] * pb_y + g_0_yyy_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyzzzz_0[i] = 3.0 * g_0_yy_0_xxxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyzzzz_1[i] * fti_ab_0 +
                                   g_0_yyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzzzz_0[i] * pb_y + g_0_yyy_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxzzzzz_0[i] = 3.0 * g_0_yy_0_xxxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxzzzzz_0[i] * pb_y +
                                   g_0_yyy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyyyy_0[i] = 3.0 * g_0_yy_0_xxyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyyyy_1[i] * fti_ab_0 +
                                   6.0 * g_0_yyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyyy_0[i] * pb_y + g_0_yyy_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyyyz_0[i] = 3.0 * g_0_yy_0_xxyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyyyz_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyyz_0[i] * pb_y + g_0_yyy_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyyzz_0[i] = 3.0 * g_0_yy_0_xxyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyyzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyzz_0[i] * pb_y + g_0_yyy_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyzzz_0[i] = 3.0 * g_0_yy_0_xxyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyzzz_0[i] * pb_y + g_0_yyy_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyzzzz_0[i] = 3.0 * g_0_yy_0_xxyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzzzz_0[i] * pb_y + g_0_yyy_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyzzzzz_0[i] = 3.0 * g_0_yy_0_xxyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyzzzzz_1[i] * fti_ab_0 +
                                   g_0_yyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzzzz_0[i] * pb_y + g_0_yyy_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxzzzzzz_0[i] = 3.0 * g_0_yy_0_xxzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzzzzzz_0[i] * pb_y +
                                   g_0_yyy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyyyy_0[i] = 3.0 * g_0_yy_0_xyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyyyy_1[i] * fti_ab_0 +
                                   7.0 * g_0_yyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyyy_0[i] * pb_y + g_0_yyy_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyyyz_0[i] = 3.0 * g_0_yy_0_xyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyyyz_1[i] * fti_ab_0 +
                                   6.0 * g_0_yyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyyz_0[i] * pb_y + g_0_yyy_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyyzz_0[i] = 3.0 * g_0_yy_0_xyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyyzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyzz_0[i] * pb_y + g_0_yyy_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyzzz_0[i] = 3.0 * g_0_yy_0_xyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyzzz_0[i] * pb_y + g_0_yyy_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyzzzz_0[i] = 3.0 * g_0_yy_0_xyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyzzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzzzz_0[i] * pb_y + g_0_yyy_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyzzzzz_0[i] = 3.0 * g_0_yy_0_xyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyzzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzzzz_0[i] * pb_y + g_0_yyy_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyzzzzzz_0[i] = 3.0 * g_0_yy_0_xyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyzzzzzz_1[i] * fti_ab_0 +
                                   g_0_yyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzzzz_0[i] * pb_y + g_0_yyy_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xzzzzzzz_0[i] = 3.0 * g_0_yy_0_xzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzzzzzz_0[i] * pb_y +
                                   g_0_yyy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyyyy_0[i] = 3.0 * g_0_yy_0_yyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyyyy_1[i] * fti_ab_0 +
                                   8.0 * g_0_yyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyyy_0[i] * pb_y + g_0_yyy_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyyyz_0[i] = 3.0 * g_0_yy_0_yyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyyyz_1[i] * fti_ab_0 +
                                   7.0 * g_0_yyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyyz_0[i] * pb_y + g_0_yyy_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyyzz_0[i] = 3.0 * g_0_yy_0_yyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyyzz_1[i] * fti_ab_0 +
                                   6.0 * g_0_yyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyzz_0[i] * pb_y + g_0_yyy_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyzzz_0[i] = 3.0 * g_0_yy_0_yyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyzzz_0[i] * pb_y + g_0_yyy_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyzzzz_0[i] = 3.0 * g_0_yy_0_yyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyzzzz_0[i] * pb_y + g_0_yyy_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyzzzzz_0[i] = 3.0 * g_0_yy_0_yyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyzzzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyzzzzz_0[i] * pb_y + g_0_yyy_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyzzzzzz_0[i] = 3.0 * g_0_yy_0_yyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyzzzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzzzzzz_0[i] * pb_y + g_0_yyy_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yzzzzzzz_0[i] = 3.0 * g_0_yy_0_yzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yzzzzzzz_1[i] * fti_ab_0 +
                                   g_0_yyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzzzzzz_0[i] * pb_y + g_0_yyy_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_zzzzzzzz_0[i] = 3.0 * g_0_yy_0_zzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzzzzzz_0[i] * pb_y +
                                   g_0_yyy_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 495-540 components of targeted buffer : SGSL

    auto g_0_yyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 495);

    auto g_0_yyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 496);

    auto g_0_yyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 497);

    auto g_0_yyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 498);

    auto g_0_yyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 499);

    auto g_0_yyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 500);

    auto g_0_yyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 501);

    auto g_0_yyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 502);

    auto g_0_yyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 503);

    auto g_0_yyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 504);

    auto g_0_yyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 505);

    auto g_0_yyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 506);

    auto g_0_yyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 507);

    auto g_0_yyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 508);

    auto g_0_yyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 509);

    auto g_0_yyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 510);

    auto g_0_yyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 511);

    auto g_0_yyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 512);

    auto g_0_yyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 513);

    auto g_0_yyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 514);

    auto g_0_yyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 515);

    auto g_0_yyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 516);

    auto g_0_yyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 517);

    auto g_0_yyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 518);

    auto g_0_yyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 519);

    auto g_0_yyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 520);

    auto g_0_yyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 521);

    auto g_0_yyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 522);

    auto g_0_yyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 523);

    auto g_0_yyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 524);

    auto g_0_yyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 525);

    auto g_0_yyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 526);

    auto g_0_yyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 527);

    auto g_0_yyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 528);

    auto g_0_yyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 529);

    auto g_0_yyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 530);

    auto g_0_yyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 531);

    auto g_0_yyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 532);

    auto g_0_yyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 533);

    auto g_0_yyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 534);

    auto g_0_yyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 535);

    auto g_0_yyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 536);

    auto g_0_yyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 537);

    auto g_0_yyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 538);

    auto g_0_yyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 539);

#pragma omp simd aligned(g_0_yyy_0_xxxxxxx_1,       \
                             g_0_yyy_0_xxxxxxxx_0,  \
                             g_0_yyy_0_xxxxxxxx_1,  \
                             g_0_yyy_0_xxxxxxxy_0,  \
                             g_0_yyy_0_xxxxxxxy_1,  \
                             g_0_yyy_0_xxxxxxxz_0,  \
                             g_0_yyy_0_xxxxxxxz_1,  \
                             g_0_yyy_0_xxxxxxy_1,   \
                             g_0_yyy_0_xxxxxxyy_0,  \
                             g_0_yyy_0_xxxxxxyy_1,  \
                             g_0_yyy_0_xxxxxxyz_0,  \
                             g_0_yyy_0_xxxxxxyz_1,  \
                             g_0_yyy_0_xxxxxxz_1,   \
                             g_0_yyy_0_xxxxxxzz_0,  \
                             g_0_yyy_0_xxxxxxzz_1,  \
                             g_0_yyy_0_xxxxxyy_1,   \
                             g_0_yyy_0_xxxxxyyy_0,  \
                             g_0_yyy_0_xxxxxyyy_1,  \
                             g_0_yyy_0_xxxxxyyz_0,  \
                             g_0_yyy_0_xxxxxyyz_1,  \
                             g_0_yyy_0_xxxxxyz_1,   \
                             g_0_yyy_0_xxxxxyzz_0,  \
                             g_0_yyy_0_xxxxxyzz_1,  \
                             g_0_yyy_0_xxxxxzz_1,   \
                             g_0_yyy_0_xxxxxzzz_0,  \
                             g_0_yyy_0_xxxxxzzz_1,  \
                             g_0_yyy_0_xxxxyyy_1,   \
                             g_0_yyy_0_xxxxyyyy_0,  \
                             g_0_yyy_0_xxxxyyyy_1,  \
                             g_0_yyy_0_xxxxyyyz_0,  \
                             g_0_yyy_0_xxxxyyyz_1,  \
                             g_0_yyy_0_xxxxyyz_1,   \
                             g_0_yyy_0_xxxxyyzz_0,  \
                             g_0_yyy_0_xxxxyyzz_1,  \
                             g_0_yyy_0_xxxxyzz_1,   \
                             g_0_yyy_0_xxxxyzzz_0,  \
                             g_0_yyy_0_xxxxyzzz_1,  \
                             g_0_yyy_0_xxxxzzz_1,   \
                             g_0_yyy_0_xxxxzzzz_0,  \
                             g_0_yyy_0_xxxxzzzz_1,  \
                             g_0_yyy_0_xxxyyyy_1,   \
                             g_0_yyy_0_xxxyyyyy_0,  \
                             g_0_yyy_0_xxxyyyyy_1,  \
                             g_0_yyy_0_xxxyyyyz_0,  \
                             g_0_yyy_0_xxxyyyyz_1,  \
                             g_0_yyy_0_xxxyyyz_1,   \
                             g_0_yyy_0_xxxyyyzz_0,  \
                             g_0_yyy_0_xxxyyyzz_1,  \
                             g_0_yyy_0_xxxyyzz_1,   \
                             g_0_yyy_0_xxxyyzzz_0,  \
                             g_0_yyy_0_xxxyyzzz_1,  \
                             g_0_yyy_0_xxxyzzz_1,   \
                             g_0_yyy_0_xxxyzzzz_0,  \
                             g_0_yyy_0_xxxyzzzz_1,  \
                             g_0_yyy_0_xxxzzzz_1,   \
                             g_0_yyy_0_xxxzzzzz_0,  \
                             g_0_yyy_0_xxxzzzzz_1,  \
                             g_0_yyy_0_xxyyyyy_1,   \
                             g_0_yyy_0_xxyyyyyy_0,  \
                             g_0_yyy_0_xxyyyyyy_1,  \
                             g_0_yyy_0_xxyyyyyz_0,  \
                             g_0_yyy_0_xxyyyyyz_1,  \
                             g_0_yyy_0_xxyyyyz_1,   \
                             g_0_yyy_0_xxyyyyzz_0,  \
                             g_0_yyy_0_xxyyyyzz_1,  \
                             g_0_yyy_0_xxyyyzz_1,   \
                             g_0_yyy_0_xxyyyzzz_0,  \
                             g_0_yyy_0_xxyyyzzz_1,  \
                             g_0_yyy_0_xxyyzzz_1,   \
                             g_0_yyy_0_xxyyzzzz_0,  \
                             g_0_yyy_0_xxyyzzzz_1,  \
                             g_0_yyy_0_xxyzzzz_1,   \
                             g_0_yyy_0_xxyzzzzz_0,  \
                             g_0_yyy_0_xxyzzzzz_1,  \
                             g_0_yyy_0_xxzzzzz_1,   \
                             g_0_yyy_0_xxzzzzzz_0,  \
                             g_0_yyy_0_xxzzzzzz_1,  \
                             g_0_yyy_0_xyyyyyy_1,   \
                             g_0_yyy_0_xyyyyyyy_0,  \
                             g_0_yyy_0_xyyyyyyy_1,  \
                             g_0_yyy_0_xyyyyyyz_0,  \
                             g_0_yyy_0_xyyyyyyz_1,  \
                             g_0_yyy_0_xyyyyyz_1,   \
                             g_0_yyy_0_xyyyyyzz_0,  \
                             g_0_yyy_0_xyyyyyzz_1,  \
                             g_0_yyy_0_xyyyyzz_1,   \
                             g_0_yyy_0_xyyyyzzz_0,  \
                             g_0_yyy_0_xyyyyzzz_1,  \
                             g_0_yyy_0_xyyyzzz_1,   \
                             g_0_yyy_0_xyyyzzzz_0,  \
                             g_0_yyy_0_xyyyzzzz_1,  \
                             g_0_yyy_0_xyyzzzz_1,   \
                             g_0_yyy_0_xyyzzzzz_0,  \
                             g_0_yyy_0_xyyzzzzz_1,  \
                             g_0_yyy_0_xyzzzzz_1,   \
                             g_0_yyy_0_xyzzzzzz_0,  \
                             g_0_yyy_0_xyzzzzzz_1,  \
                             g_0_yyy_0_xzzzzzz_1,   \
                             g_0_yyy_0_xzzzzzzz_0,  \
                             g_0_yyy_0_xzzzzzzz_1,  \
                             g_0_yyy_0_yyyyyyy_1,   \
                             g_0_yyy_0_yyyyyyyy_0,  \
                             g_0_yyy_0_yyyyyyyy_1,  \
                             g_0_yyy_0_yyyyyyyz_0,  \
                             g_0_yyy_0_yyyyyyyz_1,  \
                             g_0_yyy_0_yyyyyyz_1,   \
                             g_0_yyy_0_yyyyyyzz_0,  \
                             g_0_yyy_0_yyyyyyzz_1,  \
                             g_0_yyy_0_yyyyyzz_1,   \
                             g_0_yyy_0_yyyyyzzz_0,  \
                             g_0_yyy_0_yyyyyzzz_1,  \
                             g_0_yyy_0_yyyyzzz_1,   \
                             g_0_yyy_0_yyyyzzzz_0,  \
                             g_0_yyy_0_yyyyzzzz_1,  \
                             g_0_yyy_0_yyyzzzz_1,   \
                             g_0_yyy_0_yyyzzzzz_0,  \
                             g_0_yyy_0_yyyzzzzz_1,  \
                             g_0_yyy_0_yyzzzzz_1,   \
                             g_0_yyy_0_yyzzzzzz_0,  \
                             g_0_yyy_0_yyzzzzzz_1,  \
                             g_0_yyy_0_yzzzzzz_1,   \
                             g_0_yyy_0_yzzzzzzz_0,  \
                             g_0_yyy_0_yzzzzzzz_1,  \
                             g_0_yyy_0_zzzzzzz_1,   \
                             g_0_yyy_0_zzzzzzzz_0,  \
                             g_0_yyy_0_zzzzzzzz_1,  \
                             g_0_yyyz_0_xxxxxxxx_0, \
                             g_0_yyyz_0_xxxxxxxy_0, \
                             g_0_yyyz_0_xxxxxxxz_0, \
                             g_0_yyyz_0_xxxxxxyy_0, \
                             g_0_yyyz_0_xxxxxxyz_0, \
                             g_0_yyyz_0_xxxxxxzz_0, \
                             g_0_yyyz_0_xxxxxyyy_0, \
                             g_0_yyyz_0_xxxxxyyz_0, \
                             g_0_yyyz_0_xxxxxyzz_0, \
                             g_0_yyyz_0_xxxxxzzz_0, \
                             g_0_yyyz_0_xxxxyyyy_0, \
                             g_0_yyyz_0_xxxxyyyz_0, \
                             g_0_yyyz_0_xxxxyyzz_0, \
                             g_0_yyyz_0_xxxxyzzz_0, \
                             g_0_yyyz_0_xxxxzzzz_0, \
                             g_0_yyyz_0_xxxyyyyy_0, \
                             g_0_yyyz_0_xxxyyyyz_0, \
                             g_0_yyyz_0_xxxyyyzz_0, \
                             g_0_yyyz_0_xxxyyzzz_0, \
                             g_0_yyyz_0_xxxyzzzz_0, \
                             g_0_yyyz_0_xxxzzzzz_0, \
                             g_0_yyyz_0_xxyyyyyy_0, \
                             g_0_yyyz_0_xxyyyyyz_0, \
                             g_0_yyyz_0_xxyyyyzz_0, \
                             g_0_yyyz_0_xxyyyzzz_0, \
                             g_0_yyyz_0_xxyyzzzz_0, \
                             g_0_yyyz_0_xxyzzzzz_0, \
                             g_0_yyyz_0_xxzzzzzz_0, \
                             g_0_yyyz_0_xyyyyyyy_0, \
                             g_0_yyyz_0_xyyyyyyz_0, \
                             g_0_yyyz_0_xyyyyyzz_0, \
                             g_0_yyyz_0_xyyyyzzz_0, \
                             g_0_yyyz_0_xyyyzzzz_0, \
                             g_0_yyyz_0_xyyzzzzz_0, \
                             g_0_yyyz_0_xyzzzzzz_0, \
                             g_0_yyyz_0_xzzzzzzz_0, \
                             g_0_yyyz_0_yyyyyyyy_0, \
                             g_0_yyyz_0_yyyyyyyz_0, \
                             g_0_yyyz_0_yyyyyyzz_0, \
                             g_0_yyyz_0_yyyyyzzz_0, \
                             g_0_yyyz_0_yyyyzzzz_0, \
                             g_0_yyyz_0_yyyzzzzz_0, \
                             g_0_yyyz_0_yyzzzzzz_0, \
                             g_0_yyyz_0_yzzzzzzz_0, \
                             g_0_yyyz_0_zzzzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyz_0_xxxxxxxx_0[i] = g_0_yyy_0_xxxxxxxx_0[i] * pb_z + g_0_yyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxxxy_0[i] = g_0_yyy_0_xxxxxxxy_0[i] * pb_z + g_0_yyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxxxz_0[i] = g_0_yyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxxz_0[i] * pb_z + g_0_yyy_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxxyy_0[i] = g_0_yyy_0_xxxxxxyy_0[i] * pb_z + g_0_yyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxxyz_0[i] = g_0_yyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxyz_0[i] * pb_z + g_0_yyy_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxxzz_0[i] = 2.0 * g_0_yyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxzz_0[i] * pb_z + g_0_yyy_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxyyy_0[i] = g_0_yyy_0_xxxxxyyy_0[i] * pb_z + g_0_yyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxyyz_0[i] = g_0_yyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyyz_0[i] * pb_z + g_0_yyy_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxyzz_0[i] = 2.0 * g_0_yyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyzz_0[i] * pb_z + g_0_yyy_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxzzz_0[i] = 3.0 * g_0_yyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxzzz_0[i] * pb_z + g_0_yyy_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyyyy_0[i] = g_0_yyy_0_xxxxyyyy_0[i] * pb_z + g_0_yyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyyyz_0[i] = g_0_yyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyyz_0[i] * pb_z + g_0_yyy_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyyzz_0[i] = 2.0 * g_0_yyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyzz_0[i] * pb_z + g_0_yyy_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyzzz_0[i] = 3.0 * g_0_yyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyzzz_0[i] * pb_z + g_0_yyy_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxzzzz_0[i] = 4.0 * g_0_yyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxzzzz_0[i] * pb_z + g_0_yyy_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyyyy_0[i] = g_0_yyy_0_xxxyyyyy_0[i] * pb_z + g_0_yyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyyyz_0[i] = g_0_yyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyyz_0[i] * pb_z + g_0_yyy_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyyzz_0[i] = 2.0 * g_0_yyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyzz_0[i] * pb_z + g_0_yyy_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyzzz_0[i] = 3.0 * g_0_yyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyzzz_0[i] * pb_z + g_0_yyy_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyzzzz_0[i] = 4.0 * g_0_yyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzzzz_0[i] * pb_z + g_0_yyy_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxzzzzz_0[i] = 5.0 * g_0_yyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzzzzz_0[i] * pb_z + g_0_yyy_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyyyy_0[i] = g_0_yyy_0_xxyyyyyy_0[i] * pb_z + g_0_yyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyyyz_0[i] = g_0_yyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyyz_0[i] * pb_z + g_0_yyy_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyyzz_0[i] = 2.0 * g_0_yyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyzz_0[i] * pb_z + g_0_yyy_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyzzz_0[i] = 3.0 * g_0_yyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyzzz_0[i] * pb_z + g_0_yyy_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyzzzz_0[i] = 4.0 * g_0_yyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzzzz_0[i] * pb_z + g_0_yyy_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyzzzzz_0[i] = 5.0 * g_0_yyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzzzz_0[i] * pb_z + g_0_yyy_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxzzzzzz_0[i] = 6.0 * g_0_yyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzzzzz_0[i] * pb_z + g_0_yyy_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyyyy_0[i] = g_0_yyy_0_xyyyyyyy_0[i] * pb_z + g_0_yyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyyyz_0[i] = g_0_yyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyyz_0[i] * pb_z + g_0_yyy_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyyzz_0[i] = 2.0 * g_0_yyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyzz_0[i] * pb_z + g_0_yyy_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyzzz_0[i] = 3.0 * g_0_yyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyzzz_0[i] * pb_z + g_0_yyy_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyzzzz_0[i] = 4.0 * g_0_yyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzzzz_0[i] * pb_z + g_0_yyy_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyzzzzz_0[i] = 5.0 * g_0_yyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzzzz_0[i] * pb_z + g_0_yyy_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyzzzzzz_0[i] = 6.0 * g_0_yyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzzzz_0[i] * pb_z + g_0_yyy_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xzzzzzzz_0[i] = 7.0 * g_0_yyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzzzzz_0[i] * pb_z + g_0_yyy_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyyyy_0[i] = g_0_yyy_0_yyyyyyyy_0[i] * pb_z + g_0_yyy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyyyz_0[i] = g_0_yyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyyz_0[i] * pb_z + g_0_yyy_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyyzz_0[i] = 2.0 * g_0_yyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyzz_0[i] * pb_z + g_0_yyy_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyzzz_0[i] = 3.0 * g_0_yyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyzzz_0[i] * pb_z + g_0_yyy_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyzzzz_0[i] = 4.0 * g_0_yyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyzzzz_0[i] * pb_z + g_0_yyy_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyzzzzz_0[i] = 5.0 * g_0_yyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyzzzzz_0[i] * pb_z + g_0_yyy_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyzzzzzz_0[i] = 6.0 * g_0_yyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzzzzzz_0[i] * pb_z + g_0_yyy_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yzzzzzzz_0[i] = 7.0 * g_0_yyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzzzzzz_0[i] * pb_z + g_0_yyy_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_zzzzzzzz_0[i] = 8.0 * g_0_yyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_zzzzzzzz_0[i] * pb_z + g_0_yyy_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 540-585 components of targeted buffer : SGSL

    auto g_0_yyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 540);

    auto g_0_yyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 541);

    auto g_0_yyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 542);

    auto g_0_yyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 543);

    auto g_0_yyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 544);

    auto g_0_yyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 545);

    auto g_0_yyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 546);

    auto g_0_yyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 547);

    auto g_0_yyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 548);

    auto g_0_yyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 549);

    auto g_0_yyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 550);

    auto g_0_yyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 551);

    auto g_0_yyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 552);

    auto g_0_yyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 553);

    auto g_0_yyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 554);

    auto g_0_yyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 555);

    auto g_0_yyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 556);

    auto g_0_yyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 557);

    auto g_0_yyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 558);

    auto g_0_yyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 559);

    auto g_0_yyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 560);

    auto g_0_yyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 561);

    auto g_0_yyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 562);

    auto g_0_yyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 563);

    auto g_0_yyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 564);

    auto g_0_yyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 565);

    auto g_0_yyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 566);

    auto g_0_yyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 567);

    auto g_0_yyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 568);

    auto g_0_yyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 569);

    auto g_0_yyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 570);

    auto g_0_yyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 571);

    auto g_0_yyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 572);

    auto g_0_yyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 573);

    auto g_0_yyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 574);

    auto g_0_yyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 575);

    auto g_0_yyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 576);

    auto g_0_yyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 577);

    auto g_0_yyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 578);

    auto g_0_yyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 579);

    auto g_0_yyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 580);

    auto g_0_yyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 581);

    auto g_0_yyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 582);

    auto g_0_yyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 583);

    auto g_0_yyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 584);

#pragma omp simd aligned(g_0_yy_0_xxxxxxxy_0,       \
                             g_0_yy_0_xxxxxxxy_1,   \
                             g_0_yy_0_xxxxxxyy_0,   \
                             g_0_yy_0_xxxxxxyy_1,   \
                             g_0_yy_0_xxxxxyyy_0,   \
                             g_0_yy_0_xxxxxyyy_1,   \
                             g_0_yy_0_xxxxyyyy_0,   \
                             g_0_yy_0_xxxxyyyy_1,   \
                             g_0_yy_0_xxxyyyyy_0,   \
                             g_0_yy_0_xxxyyyyy_1,   \
                             g_0_yy_0_xxyyyyyy_0,   \
                             g_0_yy_0_xxyyyyyy_1,   \
                             g_0_yy_0_xyyyyyyy_0,   \
                             g_0_yy_0_xyyyyyyy_1,   \
                             g_0_yy_0_yyyyyyyy_0,   \
                             g_0_yy_0_yyyyyyyy_1,   \
                             g_0_yyz_0_xxxxxxxy_0,  \
                             g_0_yyz_0_xxxxxxxy_1,  \
                             g_0_yyz_0_xxxxxxyy_0,  \
                             g_0_yyz_0_xxxxxxyy_1,  \
                             g_0_yyz_0_xxxxxyyy_0,  \
                             g_0_yyz_0_xxxxxyyy_1,  \
                             g_0_yyz_0_xxxxyyyy_0,  \
                             g_0_yyz_0_xxxxyyyy_1,  \
                             g_0_yyz_0_xxxyyyyy_0,  \
                             g_0_yyz_0_xxxyyyyy_1,  \
                             g_0_yyz_0_xxyyyyyy_0,  \
                             g_0_yyz_0_xxyyyyyy_1,  \
                             g_0_yyz_0_xyyyyyyy_0,  \
                             g_0_yyz_0_xyyyyyyy_1,  \
                             g_0_yyz_0_yyyyyyyy_0,  \
                             g_0_yyz_0_yyyyyyyy_1,  \
                             g_0_yyzz_0_xxxxxxxx_0, \
                             g_0_yyzz_0_xxxxxxxy_0, \
                             g_0_yyzz_0_xxxxxxxz_0, \
                             g_0_yyzz_0_xxxxxxyy_0, \
                             g_0_yyzz_0_xxxxxxyz_0, \
                             g_0_yyzz_0_xxxxxxzz_0, \
                             g_0_yyzz_0_xxxxxyyy_0, \
                             g_0_yyzz_0_xxxxxyyz_0, \
                             g_0_yyzz_0_xxxxxyzz_0, \
                             g_0_yyzz_0_xxxxxzzz_0, \
                             g_0_yyzz_0_xxxxyyyy_0, \
                             g_0_yyzz_0_xxxxyyyz_0, \
                             g_0_yyzz_0_xxxxyyzz_0, \
                             g_0_yyzz_0_xxxxyzzz_0, \
                             g_0_yyzz_0_xxxxzzzz_0, \
                             g_0_yyzz_0_xxxyyyyy_0, \
                             g_0_yyzz_0_xxxyyyyz_0, \
                             g_0_yyzz_0_xxxyyyzz_0, \
                             g_0_yyzz_0_xxxyyzzz_0, \
                             g_0_yyzz_0_xxxyzzzz_0, \
                             g_0_yyzz_0_xxxzzzzz_0, \
                             g_0_yyzz_0_xxyyyyyy_0, \
                             g_0_yyzz_0_xxyyyyyz_0, \
                             g_0_yyzz_0_xxyyyyzz_0, \
                             g_0_yyzz_0_xxyyyzzz_0, \
                             g_0_yyzz_0_xxyyzzzz_0, \
                             g_0_yyzz_0_xxyzzzzz_0, \
                             g_0_yyzz_0_xxzzzzzz_0, \
                             g_0_yyzz_0_xyyyyyyy_0, \
                             g_0_yyzz_0_xyyyyyyz_0, \
                             g_0_yyzz_0_xyyyyyzz_0, \
                             g_0_yyzz_0_xyyyyzzz_0, \
                             g_0_yyzz_0_xyyyzzzz_0, \
                             g_0_yyzz_0_xyyzzzzz_0, \
                             g_0_yyzz_0_xyzzzzzz_0, \
                             g_0_yyzz_0_xzzzzzzz_0, \
                             g_0_yyzz_0_yyyyyyyy_0, \
                             g_0_yyzz_0_yyyyyyyz_0, \
                             g_0_yyzz_0_yyyyyyzz_0, \
                             g_0_yyzz_0_yyyyyzzz_0, \
                             g_0_yyzz_0_yyyyzzzz_0, \
                             g_0_yyzz_0_yyyzzzzz_0, \
                             g_0_yyzz_0_yyzzzzzz_0, \
                             g_0_yyzz_0_yzzzzzzz_0, \
                             g_0_yyzz_0_zzzzzzzz_0, \
                             g_0_yzz_0_xxxxxxxx_0,  \
                             g_0_yzz_0_xxxxxxxx_1,  \
                             g_0_yzz_0_xxxxxxxz_0,  \
                             g_0_yzz_0_xxxxxxxz_1,  \
                             g_0_yzz_0_xxxxxxyz_0,  \
                             g_0_yzz_0_xxxxxxyz_1,  \
                             g_0_yzz_0_xxxxxxz_1,   \
                             g_0_yzz_0_xxxxxxzz_0,  \
                             g_0_yzz_0_xxxxxxzz_1,  \
                             g_0_yzz_0_xxxxxyyz_0,  \
                             g_0_yzz_0_xxxxxyyz_1,  \
                             g_0_yzz_0_xxxxxyz_1,   \
                             g_0_yzz_0_xxxxxyzz_0,  \
                             g_0_yzz_0_xxxxxyzz_1,  \
                             g_0_yzz_0_xxxxxzz_1,   \
                             g_0_yzz_0_xxxxxzzz_0,  \
                             g_0_yzz_0_xxxxxzzz_1,  \
                             g_0_yzz_0_xxxxyyyz_0,  \
                             g_0_yzz_0_xxxxyyyz_1,  \
                             g_0_yzz_0_xxxxyyz_1,   \
                             g_0_yzz_0_xxxxyyzz_0,  \
                             g_0_yzz_0_xxxxyyzz_1,  \
                             g_0_yzz_0_xxxxyzz_1,   \
                             g_0_yzz_0_xxxxyzzz_0,  \
                             g_0_yzz_0_xxxxyzzz_1,  \
                             g_0_yzz_0_xxxxzzz_1,   \
                             g_0_yzz_0_xxxxzzzz_0,  \
                             g_0_yzz_0_xxxxzzzz_1,  \
                             g_0_yzz_0_xxxyyyyz_0,  \
                             g_0_yzz_0_xxxyyyyz_1,  \
                             g_0_yzz_0_xxxyyyz_1,   \
                             g_0_yzz_0_xxxyyyzz_0,  \
                             g_0_yzz_0_xxxyyyzz_1,  \
                             g_0_yzz_0_xxxyyzz_1,   \
                             g_0_yzz_0_xxxyyzzz_0,  \
                             g_0_yzz_0_xxxyyzzz_1,  \
                             g_0_yzz_0_xxxyzzz_1,   \
                             g_0_yzz_0_xxxyzzzz_0,  \
                             g_0_yzz_0_xxxyzzzz_1,  \
                             g_0_yzz_0_xxxzzzz_1,   \
                             g_0_yzz_0_xxxzzzzz_0,  \
                             g_0_yzz_0_xxxzzzzz_1,  \
                             g_0_yzz_0_xxyyyyyz_0,  \
                             g_0_yzz_0_xxyyyyyz_1,  \
                             g_0_yzz_0_xxyyyyz_1,   \
                             g_0_yzz_0_xxyyyyzz_0,  \
                             g_0_yzz_0_xxyyyyzz_1,  \
                             g_0_yzz_0_xxyyyzz_1,   \
                             g_0_yzz_0_xxyyyzzz_0,  \
                             g_0_yzz_0_xxyyyzzz_1,  \
                             g_0_yzz_0_xxyyzzz_1,   \
                             g_0_yzz_0_xxyyzzzz_0,  \
                             g_0_yzz_0_xxyyzzzz_1,  \
                             g_0_yzz_0_xxyzzzz_1,   \
                             g_0_yzz_0_xxyzzzzz_0,  \
                             g_0_yzz_0_xxyzzzzz_1,  \
                             g_0_yzz_0_xxzzzzz_1,   \
                             g_0_yzz_0_xxzzzzzz_0,  \
                             g_0_yzz_0_xxzzzzzz_1,  \
                             g_0_yzz_0_xyyyyyyz_0,  \
                             g_0_yzz_0_xyyyyyyz_1,  \
                             g_0_yzz_0_xyyyyyz_1,   \
                             g_0_yzz_0_xyyyyyzz_0,  \
                             g_0_yzz_0_xyyyyyzz_1,  \
                             g_0_yzz_0_xyyyyzz_1,   \
                             g_0_yzz_0_xyyyyzzz_0,  \
                             g_0_yzz_0_xyyyyzzz_1,  \
                             g_0_yzz_0_xyyyzzz_1,   \
                             g_0_yzz_0_xyyyzzzz_0,  \
                             g_0_yzz_0_xyyyzzzz_1,  \
                             g_0_yzz_0_xyyzzzz_1,   \
                             g_0_yzz_0_xyyzzzzz_0,  \
                             g_0_yzz_0_xyyzzzzz_1,  \
                             g_0_yzz_0_xyzzzzz_1,   \
                             g_0_yzz_0_xyzzzzzz_0,  \
                             g_0_yzz_0_xyzzzzzz_1,  \
                             g_0_yzz_0_xzzzzzz_1,   \
                             g_0_yzz_0_xzzzzzzz_0,  \
                             g_0_yzz_0_xzzzzzzz_1,  \
                             g_0_yzz_0_yyyyyyyz_0,  \
                             g_0_yzz_0_yyyyyyyz_1,  \
                             g_0_yzz_0_yyyyyyz_1,   \
                             g_0_yzz_0_yyyyyyzz_0,  \
                             g_0_yzz_0_yyyyyyzz_1,  \
                             g_0_yzz_0_yyyyyzz_1,   \
                             g_0_yzz_0_yyyyyzzz_0,  \
                             g_0_yzz_0_yyyyyzzz_1,  \
                             g_0_yzz_0_yyyyzzz_1,   \
                             g_0_yzz_0_yyyyzzzz_0,  \
                             g_0_yzz_0_yyyyzzzz_1,  \
                             g_0_yzz_0_yyyzzzz_1,   \
                             g_0_yzz_0_yyyzzzzz_0,  \
                             g_0_yzz_0_yyyzzzzz_1,  \
                             g_0_yzz_0_yyzzzzz_1,   \
                             g_0_yzz_0_yyzzzzzz_0,  \
                             g_0_yzz_0_yyzzzzzz_1,  \
                             g_0_yzz_0_yzzzzzz_1,   \
                             g_0_yzz_0_yzzzzzzz_0,  \
                             g_0_yzz_0_yzzzzzzz_1,  \
                             g_0_yzz_0_zzzzzzz_1,   \
                             g_0_yzz_0_zzzzzzzz_0,  \
                             g_0_yzz_0_zzzzzzzz_1,  \
                             g_0_zz_0_xxxxxxxx_0,   \
                             g_0_zz_0_xxxxxxxx_1,   \
                             g_0_zz_0_xxxxxxxz_0,   \
                             g_0_zz_0_xxxxxxxz_1,   \
                             g_0_zz_0_xxxxxxyz_0,   \
                             g_0_zz_0_xxxxxxyz_1,   \
                             g_0_zz_0_xxxxxxzz_0,   \
                             g_0_zz_0_xxxxxxzz_1,   \
                             g_0_zz_0_xxxxxyyz_0,   \
                             g_0_zz_0_xxxxxyyz_1,   \
                             g_0_zz_0_xxxxxyzz_0,   \
                             g_0_zz_0_xxxxxyzz_1,   \
                             g_0_zz_0_xxxxxzzz_0,   \
                             g_0_zz_0_xxxxxzzz_1,   \
                             g_0_zz_0_xxxxyyyz_0,   \
                             g_0_zz_0_xxxxyyyz_1,   \
                             g_0_zz_0_xxxxyyzz_0,   \
                             g_0_zz_0_xxxxyyzz_1,   \
                             g_0_zz_0_xxxxyzzz_0,   \
                             g_0_zz_0_xxxxyzzz_1,   \
                             g_0_zz_0_xxxxzzzz_0,   \
                             g_0_zz_0_xxxxzzzz_1,   \
                             g_0_zz_0_xxxyyyyz_0,   \
                             g_0_zz_0_xxxyyyyz_1,   \
                             g_0_zz_0_xxxyyyzz_0,   \
                             g_0_zz_0_xxxyyyzz_1,   \
                             g_0_zz_0_xxxyyzzz_0,   \
                             g_0_zz_0_xxxyyzzz_1,   \
                             g_0_zz_0_xxxyzzzz_0,   \
                             g_0_zz_0_xxxyzzzz_1,   \
                             g_0_zz_0_xxxzzzzz_0,   \
                             g_0_zz_0_xxxzzzzz_1,   \
                             g_0_zz_0_xxyyyyyz_0,   \
                             g_0_zz_0_xxyyyyyz_1,   \
                             g_0_zz_0_xxyyyyzz_0,   \
                             g_0_zz_0_xxyyyyzz_1,   \
                             g_0_zz_0_xxyyyzzz_0,   \
                             g_0_zz_0_xxyyyzzz_1,   \
                             g_0_zz_0_xxyyzzzz_0,   \
                             g_0_zz_0_xxyyzzzz_1,   \
                             g_0_zz_0_xxyzzzzz_0,   \
                             g_0_zz_0_xxyzzzzz_1,   \
                             g_0_zz_0_xxzzzzzz_0,   \
                             g_0_zz_0_xxzzzzzz_1,   \
                             g_0_zz_0_xyyyyyyz_0,   \
                             g_0_zz_0_xyyyyyyz_1,   \
                             g_0_zz_0_xyyyyyzz_0,   \
                             g_0_zz_0_xyyyyyzz_1,   \
                             g_0_zz_0_xyyyyzzz_0,   \
                             g_0_zz_0_xyyyyzzz_1,   \
                             g_0_zz_0_xyyyzzzz_0,   \
                             g_0_zz_0_xyyyzzzz_1,   \
                             g_0_zz_0_xyyzzzzz_0,   \
                             g_0_zz_0_xyyzzzzz_1,   \
                             g_0_zz_0_xyzzzzzz_0,   \
                             g_0_zz_0_xyzzzzzz_1,   \
                             g_0_zz_0_xzzzzzzz_0,   \
                             g_0_zz_0_xzzzzzzz_1,   \
                             g_0_zz_0_yyyyyyyz_0,   \
                             g_0_zz_0_yyyyyyyz_1,   \
                             g_0_zz_0_yyyyyyzz_0,   \
                             g_0_zz_0_yyyyyyzz_1,   \
                             g_0_zz_0_yyyyyzzz_0,   \
                             g_0_zz_0_yyyyyzzz_1,   \
                             g_0_zz_0_yyyyzzzz_0,   \
                             g_0_zz_0_yyyyzzzz_1,   \
                             g_0_zz_0_yyyzzzzz_0,   \
                             g_0_zz_0_yyyzzzzz_1,   \
                             g_0_zz_0_yyzzzzzz_0,   \
                             g_0_zz_0_yyzzzzzz_1,   \
                             g_0_zz_0_yzzzzzzz_0,   \
                             g_0_zz_0_yzzzzzzz_1,   \
                             g_0_zz_0_zzzzzzzz_0,   \
                             g_0_zz_0_zzzzzzzz_1,   \
                             wp_y,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzz_0_xxxxxxxx_0[i] =
            g_0_zz_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxxxx_0[i] * pb_y + g_0_yzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxxxy_0[i] =
            g_0_yy_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxxxxy_0[i] * pb_z + g_0_yyz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxxxxz_0[i] =
            g_0_zz_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxxxz_0[i] * pb_y + g_0_yzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxxyy_0[i] =
            g_0_yy_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxxxyy_0[i] * pb_z + g_0_yyz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxxxyz_0[i] = g_0_zz_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxxz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxxxxyz_0[i] * pb_y + g_0_yzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxxzz_0[i] =
            g_0_zz_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxxzz_0[i] * pb_y + g_0_yzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxyyy_0[i] =
            g_0_yy_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxxyyy_0[i] * pb_z + g_0_yyz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxxyyz_0[i] = g_0_zz_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxxxxyz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxxxyyz_0[i] * pb_y + g_0_yzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxyzz_0[i] = g_0_zz_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxyzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxxxyzz_0[i] * pb_y + g_0_yzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxzzz_0[i] =
            g_0_zz_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxzzz_0[i] * pb_y + g_0_yzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxyyyy_0[i] =
            g_0_yy_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxyyyy_0[i] * pb_z + g_0_yyz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxyyyz_0[i] = g_0_zz_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xxxxyyz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxxyyyz_0[i] * pb_y + g_0_yzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxyyzz_0[i] = g_0_zz_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxxxyzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxxyyzz_0[i] * pb_y + g_0_yzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxyzzz_0[i] = g_0_zz_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxxyzzz_0[i] * pb_y + g_0_yzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxzzzz_0[i] =
            g_0_zz_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxzzzz_0[i] * pb_y + g_0_yzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyyyyy_0[i] =
            g_0_yy_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxyyyyy_0[i] * pb_z + g_0_yyz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxyyyyz_0[i] = g_0_zz_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_xxxyyyz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxyyyyz_0[i] * pb_y + g_0_yzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyyyzz_0[i] = g_0_zz_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xxxyyzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxyyyzz_0[i] * pb_y + g_0_yzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyyzzz_0[i] = g_0_zz_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxxyzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxyyzzz_0[i] * pb_y + g_0_yzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyzzzz_0[i] = g_0_zz_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxxyzzzz_0[i] * pb_y + g_0_yzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxzzzzz_0[i] =
            g_0_zz_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxzzzzz_0[i] * pb_y + g_0_yzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyyyyy_0[i] =
            g_0_yy_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxyyyyyy_0[i] * pb_z + g_0_yyz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxyyyyyz_0[i] = g_0_zz_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzz_0_xxyyyyz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxyyyyyz_0[i] * pb_y + g_0_yzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyyyzz_0[i] = g_0_zz_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_xxyyyzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxyyyyzz_0[i] * pb_y + g_0_yzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyyzzz_0[i] = g_0_zz_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xxyyzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxyyyzzz_0[i] * pb_y + g_0_yzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyzzzz_0[i] = g_0_zz_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxyzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxyyzzzz_0[i] * pb_y + g_0_yzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyzzzzz_0[i] = g_0_zz_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xxyzzzzz_0[i] * pb_y + g_0_yzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxzzzzzz_0[i] =
            g_0_zz_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzzzzzz_0[i] * pb_y + g_0_yzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyyyyy_0[i] =
            g_0_yy_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xyyyyyyy_0[i] * pb_z + g_0_yyz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xyyyyyyz_0[i] = g_0_zz_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yzz_0_xyyyyyz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xyyyyyyz_0[i] * pb_y + g_0_yzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyyyzz_0[i] = g_0_zz_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yzz_0_xyyyyzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xyyyyyzz_0[i] * pb_y + g_0_yzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyyzzz_0[i] = g_0_zz_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_xyyyzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xyyyyzzz_0[i] * pb_y + g_0_yzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyzzzz_0[i] = g_0_zz_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xyyzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xyyyzzzz_0[i] * pb_y + g_0_yzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyzzzzz_0[i] = g_0_zz_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xyzzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xyyzzzzz_0[i] * pb_y + g_0_yzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyzzzzzz_0[i] = g_0_zz_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_xyzzzzzz_0[i] * pb_y + g_0_yzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xzzzzzzz_0[i] =
            g_0_zz_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzzzzzz_0[i] * pb_y + g_0_yzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyyyyy_0[i] =
            g_0_yy_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_yyyyyyyy_0[i] * pb_z + g_0_yyz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_yyyyyyyz_0[i] = g_0_zz_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyyz_1[i] * fti_ab_0 + 7.0 * g_0_yzz_0_yyyyyyz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_yyyyyyyz_0[i] * pb_y + g_0_yzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyyyzz_0[i] = g_0_zz_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyzz_1[i] * fti_ab_0 + 6.0 * g_0_yzz_0_yyyyyzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_yyyyyyzz_0[i] * pb_y + g_0_yzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyyzzz_0[i] = g_0_zz_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyzzz_1[i] * fti_ab_0 + 5.0 * g_0_yzz_0_yyyyzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_yyyyyzzz_0[i] * pb_y + g_0_yzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyzzzz_0[i] = g_0_zz_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_yyyzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_yyyyzzzz_0[i] * pb_y + g_0_yzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyzzzzz_0[i] = g_0_zz_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_yyzzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_yyyzzzzz_0[i] * pb_y + g_0_yzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyzzzzzz_0[i] = g_0_zz_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_yzzzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_yyzzzzzz_0[i] * pb_y + g_0_yzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yzzzzzzz_0[i] = g_0_zz_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzz_0_yzzzzzzz_0[i] * pb_y + g_0_yzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_zzzzzzzz_0[i] =
            g_0_zz_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzzzzzz_0[i] * pb_y + g_0_yzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 585-630 components of targeted buffer : SGSL

    auto g_0_yzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 585);

    auto g_0_yzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 586);

    auto g_0_yzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 587);

    auto g_0_yzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 588);

    auto g_0_yzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 589);

    auto g_0_yzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 590);

    auto g_0_yzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 591);

    auto g_0_yzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 592);

    auto g_0_yzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 593);

    auto g_0_yzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 594);

    auto g_0_yzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 595);

    auto g_0_yzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 596);

    auto g_0_yzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 597);

    auto g_0_yzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 598);

    auto g_0_yzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 599);

    auto g_0_yzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 600);

    auto g_0_yzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 601);

    auto g_0_yzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 602);

    auto g_0_yzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 603);

    auto g_0_yzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 604);

    auto g_0_yzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 605);

    auto g_0_yzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 606);

    auto g_0_yzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 607);

    auto g_0_yzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 608);

    auto g_0_yzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 609);

    auto g_0_yzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 610);

    auto g_0_yzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 611);

    auto g_0_yzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 612);

    auto g_0_yzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 613);

    auto g_0_yzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 614);

    auto g_0_yzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 615);

    auto g_0_yzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 616);

    auto g_0_yzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 617);

    auto g_0_yzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 618);

    auto g_0_yzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 619);

    auto g_0_yzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 620);

    auto g_0_yzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 621);

    auto g_0_yzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 622);

    auto g_0_yzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 623);

    auto g_0_yzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 624);

    auto g_0_yzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 625);

    auto g_0_yzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 626);

    auto g_0_yzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 627);

    auto g_0_yzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 628);

    auto g_0_yzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 629);

#pragma omp simd aligned(g_0_yzzz_0_xxxxxxxx_0,     \
                             g_0_yzzz_0_xxxxxxxy_0, \
                             g_0_yzzz_0_xxxxxxxz_0, \
                             g_0_yzzz_0_xxxxxxyy_0, \
                             g_0_yzzz_0_xxxxxxyz_0, \
                             g_0_yzzz_0_xxxxxxzz_0, \
                             g_0_yzzz_0_xxxxxyyy_0, \
                             g_0_yzzz_0_xxxxxyyz_0, \
                             g_0_yzzz_0_xxxxxyzz_0, \
                             g_0_yzzz_0_xxxxxzzz_0, \
                             g_0_yzzz_0_xxxxyyyy_0, \
                             g_0_yzzz_0_xxxxyyyz_0, \
                             g_0_yzzz_0_xxxxyyzz_0, \
                             g_0_yzzz_0_xxxxyzzz_0, \
                             g_0_yzzz_0_xxxxzzzz_0, \
                             g_0_yzzz_0_xxxyyyyy_0, \
                             g_0_yzzz_0_xxxyyyyz_0, \
                             g_0_yzzz_0_xxxyyyzz_0, \
                             g_0_yzzz_0_xxxyyzzz_0, \
                             g_0_yzzz_0_xxxyzzzz_0, \
                             g_0_yzzz_0_xxxzzzzz_0, \
                             g_0_yzzz_0_xxyyyyyy_0, \
                             g_0_yzzz_0_xxyyyyyz_0, \
                             g_0_yzzz_0_xxyyyyzz_0, \
                             g_0_yzzz_0_xxyyyzzz_0, \
                             g_0_yzzz_0_xxyyzzzz_0, \
                             g_0_yzzz_0_xxyzzzzz_0, \
                             g_0_yzzz_0_xxzzzzzz_0, \
                             g_0_yzzz_0_xyyyyyyy_0, \
                             g_0_yzzz_0_xyyyyyyz_0, \
                             g_0_yzzz_0_xyyyyyzz_0, \
                             g_0_yzzz_0_xyyyyzzz_0, \
                             g_0_yzzz_0_xyyyzzzz_0, \
                             g_0_yzzz_0_xyyzzzzz_0, \
                             g_0_yzzz_0_xyzzzzzz_0, \
                             g_0_yzzz_0_xzzzzzzz_0, \
                             g_0_yzzz_0_yyyyyyyy_0, \
                             g_0_yzzz_0_yyyyyyyz_0, \
                             g_0_yzzz_0_yyyyyyzz_0, \
                             g_0_yzzz_0_yyyyyzzz_0, \
                             g_0_yzzz_0_yyyyzzzz_0, \
                             g_0_yzzz_0_yyyzzzzz_0, \
                             g_0_yzzz_0_yyzzzzzz_0, \
                             g_0_yzzz_0_yzzzzzzz_0, \
                             g_0_yzzz_0_zzzzzzzz_0, \
                             g_0_zzz_0_xxxxxxx_1,   \
                             g_0_zzz_0_xxxxxxxx_0,  \
                             g_0_zzz_0_xxxxxxxx_1,  \
                             g_0_zzz_0_xxxxxxxy_0,  \
                             g_0_zzz_0_xxxxxxxy_1,  \
                             g_0_zzz_0_xxxxxxxz_0,  \
                             g_0_zzz_0_xxxxxxxz_1,  \
                             g_0_zzz_0_xxxxxxy_1,   \
                             g_0_zzz_0_xxxxxxyy_0,  \
                             g_0_zzz_0_xxxxxxyy_1,  \
                             g_0_zzz_0_xxxxxxyz_0,  \
                             g_0_zzz_0_xxxxxxyz_1,  \
                             g_0_zzz_0_xxxxxxz_1,   \
                             g_0_zzz_0_xxxxxxzz_0,  \
                             g_0_zzz_0_xxxxxxzz_1,  \
                             g_0_zzz_0_xxxxxyy_1,   \
                             g_0_zzz_0_xxxxxyyy_0,  \
                             g_0_zzz_0_xxxxxyyy_1,  \
                             g_0_zzz_0_xxxxxyyz_0,  \
                             g_0_zzz_0_xxxxxyyz_1,  \
                             g_0_zzz_0_xxxxxyz_1,   \
                             g_0_zzz_0_xxxxxyzz_0,  \
                             g_0_zzz_0_xxxxxyzz_1,  \
                             g_0_zzz_0_xxxxxzz_1,   \
                             g_0_zzz_0_xxxxxzzz_0,  \
                             g_0_zzz_0_xxxxxzzz_1,  \
                             g_0_zzz_0_xxxxyyy_1,   \
                             g_0_zzz_0_xxxxyyyy_0,  \
                             g_0_zzz_0_xxxxyyyy_1,  \
                             g_0_zzz_0_xxxxyyyz_0,  \
                             g_0_zzz_0_xxxxyyyz_1,  \
                             g_0_zzz_0_xxxxyyz_1,   \
                             g_0_zzz_0_xxxxyyzz_0,  \
                             g_0_zzz_0_xxxxyyzz_1,  \
                             g_0_zzz_0_xxxxyzz_1,   \
                             g_0_zzz_0_xxxxyzzz_0,  \
                             g_0_zzz_0_xxxxyzzz_1,  \
                             g_0_zzz_0_xxxxzzz_1,   \
                             g_0_zzz_0_xxxxzzzz_0,  \
                             g_0_zzz_0_xxxxzzzz_1,  \
                             g_0_zzz_0_xxxyyyy_1,   \
                             g_0_zzz_0_xxxyyyyy_0,  \
                             g_0_zzz_0_xxxyyyyy_1,  \
                             g_0_zzz_0_xxxyyyyz_0,  \
                             g_0_zzz_0_xxxyyyyz_1,  \
                             g_0_zzz_0_xxxyyyz_1,   \
                             g_0_zzz_0_xxxyyyzz_0,  \
                             g_0_zzz_0_xxxyyyzz_1,  \
                             g_0_zzz_0_xxxyyzz_1,   \
                             g_0_zzz_0_xxxyyzzz_0,  \
                             g_0_zzz_0_xxxyyzzz_1,  \
                             g_0_zzz_0_xxxyzzz_1,   \
                             g_0_zzz_0_xxxyzzzz_0,  \
                             g_0_zzz_0_xxxyzzzz_1,  \
                             g_0_zzz_0_xxxzzzz_1,   \
                             g_0_zzz_0_xxxzzzzz_0,  \
                             g_0_zzz_0_xxxzzzzz_1,  \
                             g_0_zzz_0_xxyyyyy_1,   \
                             g_0_zzz_0_xxyyyyyy_0,  \
                             g_0_zzz_0_xxyyyyyy_1,  \
                             g_0_zzz_0_xxyyyyyz_0,  \
                             g_0_zzz_0_xxyyyyyz_1,  \
                             g_0_zzz_0_xxyyyyz_1,   \
                             g_0_zzz_0_xxyyyyzz_0,  \
                             g_0_zzz_0_xxyyyyzz_1,  \
                             g_0_zzz_0_xxyyyzz_1,   \
                             g_0_zzz_0_xxyyyzzz_0,  \
                             g_0_zzz_0_xxyyyzzz_1,  \
                             g_0_zzz_0_xxyyzzz_1,   \
                             g_0_zzz_0_xxyyzzzz_0,  \
                             g_0_zzz_0_xxyyzzzz_1,  \
                             g_0_zzz_0_xxyzzzz_1,   \
                             g_0_zzz_0_xxyzzzzz_0,  \
                             g_0_zzz_0_xxyzzzzz_1,  \
                             g_0_zzz_0_xxzzzzz_1,   \
                             g_0_zzz_0_xxzzzzzz_0,  \
                             g_0_zzz_0_xxzzzzzz_1,  \
                             g_0_zzz_0_xyyyyyy_1,   \
                             g_0_zzz_0_xyyyyyyy_0,  \
                             g_0_zzz_0_xyyyyyyy_1,  \
                             g_0_zzz_0_xyyyyyyz_0,  \
                             g_0_zzz_0_xyyyyyyz_1,  \
                             g_0_zzz_0_xyyyyyz_1,   \
                             g_0_zzz_0_xyyyyyzz_0,  \
                             g_0_zzz_0_xyyyyyzz_1,  \
                             g_0_zzz_0_xyyyyzz_1,   \
                             g_0_zzz_0_xyyyyzzz_0,  \
                             g_0_zzz_0_xyyyyzzz_1,  \
                             g_0_zzz_0_xyyyzzz_1,   \
                             g_0_zzz_0_xyyyzzzz_0,  \
                             g_0_zzz_0_xyyyzzzz_1,  \
                             g_0_zzz_0_xyyzzzz_1,   \
                             g_0_zzz_0_xyyzzzzz_0,  \
                             g_0_zzz_0_xyyzzzzz_1,  \
                             g_0_zzz_0_xyzzzzz_1,   \
                             g_0_zzz_0_xyzzzzzz_0,  \
                             g_0_zzz_0_xyzzzzzz_1,  \
                             g_0_zzz_0_xzzzzzz_1,   \
                             g_0_zzz_0_xzzzzzzz_0,  \
                             g_0_zzz_0_xzzzzzzz_1,  \
                             g_0_zzz_0_yyyyyyy_1,   \
                             g_0_zzz_0_yyyyyyyy_0,  \
                             g_0_zzz_0_yyyyyyyy_1,  \
                             g_0_zzz_0_yyyyyyyz_0,  \
                             g_0_zzz_0_yyyyyyyz_1,  \
                             g_0_zzz_0_yyyyyyz_1,   \
                             g_0_zzz_0_yyyyyyzz_0,  \
                             g_0_zzz_0_yyyyyyzz_1,  \
                             g_0_zzz_0_yyyyyzz_1,   \
                             g_0_zzz_0_yyyyyzzz_0,  \
                             g_0_zzz_0_yyyyyzzz_1,  \
                             g_0_zzz_0_yyyyzzz_1,   \
                             g_0_zzz_0_yyyyzzzz_0,  \
                             g_0_zzz_0_yyyyzzzz_1,  \
                             g_0_zzz_0_yyyzzzz_1,   \
                             g_0_zzz_0_yyyzzzzz_0,  \
                             g_0_zzz_0_yyyzzzzz_1,  \
                             g_0_zzz_0_yyzzzzz_1,   \
                             g_0_zzz_0_yyzzzzzz_0,  \
                             g_0_zzz_0_yyzzzzzz_1,  \
                             g_0_zzz_0_yzzzzzz_1,   \
                             g_0_zzz_0_yzzzzzzz_0,  \
                             g_0_zzz_0_yzzzzzzz_1,  \
                             g_0_zzz_0_zzzzzzz_1,   \
                             g_0_zzz_0_zzzzzzzz_0,  \
                             g_0_zzz_0_zzzzzzzz_1,  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzz_0_xxxxxxxx_0[i] = g_0_zzz_0_xxxxxxxx_0[i] * pb_y + g_0_zzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxxxy_0[i] = g_0_zzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxxy_0[i] * pb_y + g_0_zzz_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxxxz_0[i] = g_0_zzz_0_xxxxxxxz_0[i] * pb_y + g_0_zzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxxyy_0[i] = 2.0 * g_0_zzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxyy_0[i] * pb_y + g_0_zzz_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxxyz_0[i] = g_0_zzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxyz_0[i] * pb_y + g_0_zzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxxzz_0[i] = g_0_zzz_0_xxxxxxzz_0[i] * pb_y + g_0_zzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxyyy_0[i] = 3.0 * g_0_zzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyyy_0[i] * pb_y + g_0_zzz_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxyyz_0[i] = 2.0 * g_0_zzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyyz_0[i] * pb_y + g_0_zzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxyzz_0[i] = g_0_zzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyzz_0[i] * pb_y + g_0_zzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxzzz_0[i] = g_0_zzz_0_xxxxxzzz_0[i] * pb_y + g_0_zzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyyyy_0[i] = 4.0 * g_0_zzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyyy_0[i] * pb_y + g_0_zzz_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyyyz_0[i] = 3.0 * g_0_zzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyyz_0[i] * pb_y + g_0_zzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyyzz_0[i] = 2.0 * g_0_zzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyzz_0[i] * pb_y + g_0_zzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyzzz_0[i] = g_0_zzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyzzz_0[i] * pb_y + g_0_zzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxzzzz_0[i] = g_0_zzz_0_xxxxzzzz_0[i] * pb_y + g_0_zzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyyyy_0[i] = 5.0 * g_0_zzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyyy_0[i] * pb_y + g_0_zzz_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyyyz_0[i] = 4.0 * g_0_zzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyyz_0[i] * pb_y + g_0_zzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyyzz_0[i] = 3.0 * g_0_zzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyzz_0[i] * pb_y + g_0_zzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyzzz_0[i] = 2.0 * g_0_zzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyzzz_0[i] * pb_y + g_0_zzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyzzzz_0[i] = g_0_zzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzzzz_0[i] * pb_y + g_0_zzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxzzzzz_0[i] = g_0_zzz_0_xxxzzzzz_0[i] * pb_y + g_0_zzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyyyy_0[i] = 6.0 * g_0_zzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyyy_0[i] * pb_y + g_0_zzz_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyyyz_0[i] = 5.0 * g_0_zzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyyz_0[i] * pb_y + g_0_zzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyyzz_0[i] = 4.0 * g_0_zzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyzz_0[i] * pb_y + g_0_zzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyzzz_0[i] = 3.0 * g_0_zzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyzzz_0[i] * pb_y + g_0_zzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyzzzz_0[i] = 2.0 * g_0_zzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzzzz_0[i] * pb_y + g_0_zzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyzzzzz_0[i] = g_0_zzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzzzz_0[i] * pb_y + g_0_zzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxzzzzzz_0[i] = g_0_zzz_0_xxzzzzzz_0[i] * pb_y + g_0_zzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyyyy_0[i] = 7.0 * g_0_zzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyyy_0[i] * pb_y + g_0_zzz_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyyyz_0[i] = 6.0 * g_0_zzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyyz_0[i] * pb_y + g_0_zzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyyzz_0[i] = 5.0 * g_0_zzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyzz_0[i] * pb_y + g_0_zzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyzzz_0[i] = 4.0 * g_0_zzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyzzz_0[i] * pb_y + g_0_zzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyzzzz_0[i] = 3.0 * g_0_zzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzzzz_0[i] * pb_y + g_0_zzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyzzzzz_0[i] = 2.0 * g_0_zzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzzzz_0[i] * pb_y + g_0_zzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyzzzzzz_0[i] = g_0_zzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzzzz_0[i] * pb_y + g_0_zzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xzzzzzzz_0[i] = g_0_zzz_0_xzzzzzzz_0[i] * pb_y + g_0_zzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyyyy_0[i] = 8.0 * g_0_zzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyyy_0[i] * pb_y + g_0_zzz_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyyyz_0[i] = 7.0 * g_0_zzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyyz_0[i] * pb_y + g_0_zzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyyzz_0[i] = 6.0 * g_0_zzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyzz_0[i] * pb_y + g_0_zzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyzzz_0[i] = 5.0 * g_0_zzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyzzz_0[i] * pb_y + g_0_zzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyzzzz_0[i] = 4.0 * g_0_zzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyzzzz_0[i] * pb_y + g_0_zzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyzzzzz_0[i] = 3.0 * g_0_zzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyzzzzz_0[i] * pb_y + g_0_zzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyzzzzzz_0[i] = 2.0 * g_0_zzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzzzzzz_0[i] * pb_y + g_0_zzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yzzzzzzz_0[i] = g_0_zzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzzzzzz_0[i] * pb_y + g_0_zzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_zzzzzzzz_0[i] = g_0_zzz_0_zzzzzzzz_0[i] * pb_y + g_0_zzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 630-675 components of targeted buffer : SGSL

    auto g_0_zzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 630);

    auto g_0_zzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 631);

    auto g_0_zzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 632);

    auto g_0_zzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 633);

    auto g_0_zzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 634);

    auto g_0_zzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 635);

    auto g_0_zzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 636);

    auto g_0_zzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 637);

    auto g_0_zzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 638);

    auto g_0_zzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 639);

    auto g_0_zzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 640);

    auto g_0_zzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 641);

    auto g_0_zzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 642);

    auto g_0_zzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 643);

    auto g_0_zzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 644);

    auto g_0_zzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 645);

    auto g_0_zzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 646);

    auto g_0_zzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 647);

    auto g_0_zzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 648);

    auto g_0_zzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 649);

    auto g_0_zzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 650);

    auto g_0_zzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 651);

    auto g_0_zzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 652);

    auto g_0_zzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 653);

    auto g_0_zzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 654);

    auto g_0_zzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 655);

    auto g_0_zzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 656);

    auto g_0_zzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 657);

    auto g_0_zzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 658);

    auto g_0_zzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 659);

    auto g_0_zzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 660);

    auto g_0_zzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 661);

    auto g_0_zzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 662);

    auto g_0_zzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 663);

    auto g_0_zzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 664);

    auto g_0_zzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 665);

    auto g_0_zzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 666);

    auto g_0_zzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 667);

    auto g_0_zzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 668);

    auto g_0_zzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 669);

    auto g_0_zzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 670);

    auto g_0_zzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 671);

    auto g_0_zzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 672);

    auto g_0_zzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 673);

    auto g_0_zzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 674);

#pragma omp simd aligned(g_0_zz_0_xxxxxxxx_0,       \
                             g_0_zz_0_xxxxxxxx_1,   \
                             g_0_zz_0_xxxxxxxy_0,   \
                             g_0_zz_0_xxxxxxxy_1,   \
                             g_0_zz_0_xxxxxxxz_0,   \
                             g_0_zz_0_xxxxxxxz_1,   \
                             g_0_zz_0_xxxxxxyy_0,   \
                             g_0_zz_0_xxxxxxyy_1,   \
                             g_0_zz_0_xxxxxxyz_0,   \
                             g_0_zz_0_xxxxxxyz_1,   \
                             g_0_zz_0_xxxxxxzz_0,   \
                             g_0_zz_0_xxxxxxzz_1,   \
                             g_0_zz_0_xxxxxyyy_0,   \
                             g_0_zz_0_xxxxxyyy_1,   \
                             g_0_zz_0_xxxxxyyz_0,   \
                             g_0_zz_0_xxxxxyyz_1,   \
                             g_0_zz_0_xxxxxyzz_0,   \
                             g_0_zz_0_xxxxxyzz_1,   \
                             g_0_zz_0_xxxxxzzz_0,   \
                             g_0_zz_0_xxxxxzzz_1,   \
                             g_0_zz_0_xxxxyyyy_0,   \
                             g_0_zz_0_xxxxyyyy_1,   \
                             g_0_zz_0_xxxxyyyz_0,   \
                             g_0_zz_0_xxxxyyyz_1,   \
                             g_0_zz_0_xxxxyyzz_0,   \
                             g_0_zz_0_xxxxyyzz_1,   \
                             g_0_zz_0_xxxxyzzz_0,   \
                             g_0_zz_0_xxxxyzzz_1,   \
                             g_0_zz_0_xxxxzzzz_0,   \
                             g_0_zz_0_xxxxzzzz_1,   \
                             g_0_zz_0_xxxyyyyy_0,   \
                             g_0_zz_0_xxxyyyyy_1,   \
                             g_0_zz_0_xxxyyyyz_0,   \
                             g_0_zz_0_xxxyyyyz_1,   \
                             g_0_zz_0_xxxyyyzz_0,   \
                             g_0_zz_0_xxxyyyzz_1,   \
                             g_0_zz_0_xxxyyzzz_0,   \
                             g_0_zz_0_xxxyyzzz_1,   \
                             g_0_zz_0_xxxyzzzz_0,   \
                             g_0_zz_0_xxxyzzzz_1,   \
                             g_0_zz_0_xxxzzzzz_0,   \
                             g_0_zz_0_xxxzzzzz_1,   \
                             g_0_zz_0_xxyyyyyy_0,   \
                             g_0_zz_0_xxyyyyyy_1,   \
                             g_0_zz_0_xxyyyyyz_0,   \
                             g_0_zz_0_xxyyyyyz_1,   \
                             g_0_zz_0_xxyyyyzz_0,   \
                             g_0_zz_0_xxyyyyzz_1,   \
                             g_0_zz_0_xxyyyzzz_0,   \
                             g_0_zz_0_xxyyyzzz_1,   \
                             g_0_zz_0_xxyyzzzz_0,   \
                             g_0_zz_0_xxyyzzzz_1,   \
                             g_0_zz_0_xxyzzzzz_0,   \
                             g_0_zz_0_xxyzzzzz_1,   \
                             g_0_zz_0_xxzzzzzz_0,   \
                             g_0_zz_0_xxzzzzzz_1,   \
                             g_0_zz_0_xyyyyyyy_0,   \
                             g_0_zz_0_xyyyyyyy_1,   \
                             g_0_zz_0_xyyyyyyz_0,   \
                             g_0_zz_0_xyyyyyyz_1,   \
                             g_0_zz_0_xyyyyyzz_0,   \
                             g_0_zz_0_xyyyyyzz_1,   \
                             g_0_zz_0_xyyyyzzz_0,   \
                             g_0_zz_0_xyyyyzzz_1,   \
                             g_0_zz_0_xyyyzzzz_0,   \
                             g_0_zz_0_xyyyzzzz_1,   \
                             g_0_zz_0_xyyzzzzz_0,   \
                             g_0_zz_0_xyyzzzzz_1,   \
                             g_0_zz_0_xyzzzzzz_0,   \
                             g_0_zz_0_xyzzzzzz_1,   \
                             g_0_zz_0_xzzzzzzz_0,   \
                             g_0_zz_0_xzzzzzzz_1,   \
                             g_0_zz_0_yyyyyyyy_0,   \
                             g_0_zz_0_yyyyyyyy_1,   \
                             g_0_zz_0_yyyyyyyz_0,   \
                             g_0_zz_0_yyyyyyyz_1,   \
                             g_0_zz_0_yyyyyyzz_0,   \
                             g_0_zz_0_yyyyyyzz_1,   \
                             g_0_zz_0_yyyyyzzz_0,   \
                             g_0_zz_0_yyyyyzzz_1,   \
                             g_0_zz_0_yyyyzzzz_0,   \
                             g_0_zz_0_yyyyzzzz_1,   \
                             g_0_zz_0_yyyzzzzz_0,   \
                             g_0_zz_0_yyyzzzzz_1,   \
                             g_0_zz_0_yyzzzzzz_0,   \
                             g_0_zz_0_yyzzzzzz_1,   \
                             g_0_zz_0_yzzzzzzz_0,   \
                             g_0_zz_0_yzzzzzzz_1,   \
                             g_0_zz_0_zzzzzzzz_0,   \
                             g_0_zz_0_zzzzzzzz_1,   \
                             g_0_zzz_0_xxxxxxx_1,   \
                             g_0_zzz_0_xxxxxxxx_0,  \
                             g_0_zzz_0_xxxxxxxx_1,  \
                             g_0_zzz_0_xxxxxxxy_0,  \
                             g_0_zzz_0_xxxxxxxy_1,  \
                             g_0_zzz_0_xxxxxxxz_0,  \
                             g_0_zzz_0_xxxxxxxz_1,  \
                             g_0_zzz_0_xxxxxxy_1,   \
                             g_0_zzz_0_xxxxxxyy_0,  \
                             g_0_zzz_0_xxxxxxyy_1,  \
                             g_0_zzz_0_xxxxxxyz_0,  \
                             g_0_zzz_0_xxxxxxyz_1,  \
                             g_0_zzz_0_xxxxxxz_1,   \
                             g_0_zzz_0_xxxxxxzz_0,  \
                             g_0_zzz_0_xxxxxxzz_1,  \
                             g_0_zzz_0_xxxxxyy_1,   \
                             g_0_zzz_0_xxxxxyyy_0,  \
                             g_0_zzz_0_xxxxxyyy_1,  \
                             g_0_zzz_0_xxxxxyyz_0,  \
                             g_0_zzz_0_xxxxxyyz_1,  \
                             g_0_zzz_0_xxxxxyz_1,   \
                             g_0_zzz_0_xxxxxyzz_0,  \
                             g_0_zzz_0_xxxxxyzz_1,  \
                             g_0_zzz_0_xxxxxzz_1,   \
                             g_0_zzz_0_xxxxxzzz_0,  \
                             g_0_zzz_0_xxxxxzzz_1,  \
                             g_0_zzz_0_xxxxyyy_1,   \
                             g_0_zzz_0_xxxxyyyy_0,  \
                             g_0_zzz_0_xxxxyyyy_1,  \
                             g_0_zzz_0_xxxxyyyz_0,  \
                             g_0_zzz_0_xxxxyyyz_1,  \
                             g_0_zzz_0_xxxxyyz_1,   \
                             g_0_zzz_0_xxxxyyzz_0,  \
                             g_0_zzz_0_xxxxyyzz_1,  \
                             g_0_zzz_0_xxxxyzz_1,   \
                             g_0_zzz_0_xxxxyzzz_0,  \
                             g_0_zzz_0_xxxxyzzz_1,  \
                             g_0_zzz_0_xxxxzzz_1,   \
                             g_0_zzz_0_xxxxzzzz_0,  \
                             g_0_zzz_0_xxxxzzzz_1,  \
                             g_0_zzz_0_xxxyyyy_1,   \
                             g_0_zzz_0_xxxyyyyy_0,  \
                             g_0_zzz_0_xxxyyyyy_1,  \
                             g_0_zzz_0_xxxyyyyz_0,  \
                             g_0_zzz_0_xxxyyyyz_1,  \
                             g_0_zzz_0_xxxyyyz_1,   \
                             g_0_zzz_0_xxxyyyzz_0,  \
                             g_0_zzz_0_xxxyyyzz_1,  \
                             g_0_zzz_0_xxxyyzz_1,   \
                             g_0_zzz_0_xxxyyzzz_0,  \
                             g_0_zzz_0_xxxyyzzz_1,  \
                             g_0_zzz_0_xxxyzzz_1,   \
                             g_0_zzz_0_xxxyzzzz_0,  \
                             g_0_zzz_0_xxxyzzzz_1,  \
                             g_0_zzz_0_xxxzzzz_1,   \
                             g_0_zzz_0_xxxzzzzz_0,  \
                             g_0_zzz_0_xxxzzzzz_1,  \
                             g_0_zzz_0_xxyyyyy_1,   \
                             g_0_zzz_0_xxyyyyyy_0,  \
                             g_0_zzz_0_xxyyyyyy_1,  \
                             g_0_zzz_0_xxyyyyyz_0,  \
                             g_0_zzz_0_xxyyyyyz_1,  \
                             g_0_zzz_0_xxyyyyz_1,   \
                             g_0_zzz_0_xxyyyyzz_0,  \
                             g_0_zzz_0_xxyyyyzz_1,  \
                             g_0_zzz_0_xxyyyzz_1,   \
                             g_0_zzz_0_xxyyyzzz_0,  \
                             g_0_zzz_0_xxyyyzzz_1,  \
                             g_0_zzz_0_xxyyzzz_1,   \
                             g_0_zzz_0_xxyyzzzz_0,  \
                             g_0_zzz_0_xxyyzzzz_1,  \
                             g_0_zzz_0_xxyzzzz_1,   \
                             g_0_zzz_0_xxyzzzzz_0,  \
                             g_0_zzz_0_xxyzzzzz_1,  \
                             g_0_zzz_0_xxzzzzz_1,   \
                             g_0_zzz_0_xxzzzzzz_0,  \
                             g_0_zzz_0_xxzzzzzz_1,  \
                             g_0_zzz_0_xyyyyyy_1,   \
                             g_0_zzz_0_xyyyyyyy_0,  \
                             g_0_zzz_0_xyyyyyyy_1,  \
                             g_0_zzz_0_xyyyyyyz_0,  \
                             g_0_zzz_0_xyyyyyyz_1,  \
                             g_0_zzz_0_xyyyyyz_1,   \
                             g_0_zzz_0_xyyyyyzz_0,  \
                             g_0_zzz_0_xyyyyyzz_1,  \
                             g_0_zzz_0_xyyyyzz_1,   \
                             g_0_zzz_0_xyyyyzzz_0,  \
                             g_0_zzz_0_xyyyyzzz_1,  \
                             g_0_zzz_0_xyyyzzz_1,   \
                             g_0_zzz_0_xyyyzzzz_0,  \
                             g_0_zzz_0_xyyyzzzz_1,  \
                             g_0_zzz_0_xyyzzzz_1,   \
                             g_0_zzz_0_xyyzzzzz_0,  \
                             g_0_zzz_0_xyyzzzzz_1,  \
                             g_0_zzz_0_xyzzzzz_1,   \
                             g_0_zzz_0_xyzzzzzz_0,  \
                             g_0_zzz_0_xyzzzzzz_1,  \
                             g_0_zzz_0_xzzzzzz_1,   \
                             g_0_zzz_0_xzzzzzzz_0,  \
                             g_0_zzz_0_xzzzzzzz_1,  \
                             g_0_zzz_0_yyyyyyy_1,   \
                             g_0_zzz_0_yyyyyyyy_0,  \
                             g_0_zzz_0_yyyyyyyy_1,  \
                             g_0_zzz_0_yyyyyyyz_0,  \
                             g_0_zzz_0_yyyyyyyz_1,  \
                             g_0_zzz_0_yyyyyyz_1,   \
                             g_0_zzz_0_yyyyyyzz_0,  \
                             g_0_zzz_0_yyyyyyzz_1,  \
                             g_0_zzz_0_yyyyyzz_1,   \
                             g_0_zzz_0_yyyyyzzz_0,  \
                             g_0_zzz_0_yyyyyzzz_1,  \
                             g_0_zzz_0_yyyyzzz_1,   \
                             g_0_zzz_0_yyyyzzzz_0,  \
                             g_0_zzz_0_yyyyzzzz_1,  \
                             g_0_zzz_0_yyyzzzz_1,   \
                             g_0_zzz_0_yyyzzzzz_0,  \
                             g_0_zzz_0_yyyzzzzz_1,  \
                             g_0_zzz_0_yyzzzzz_1,   \
                             g_0_zzz_0_yyzzzzzz_0,  \
                             g_0_zzz_0_yyzzzzzz_1,  \
                             g_0_zzz_0_yzzzzzz_1,   \
                             g_0_zzz_0_yzzzzzzz_0,  \
                             g_0_zzz_0_yzzzzzzz_1,  \
                             g_0_zzz_0_zzzzzzz_1,   \
                             g_0_zzz_0_zzzzzzzz_0,  \
                             g_0_zzz_0_zzzzzzzz_1,  \
                             g_0_zzzz_0_xxxxxxxx_0, \
                             g_0_zzzz_0_xxxxxxxy_0, \
                             g_0_zzzz_0_xxxxxxxz_0, \
                             g_0_zzzz_0_xxxxxxyy_0, \
                             g_0_zzzz_0_xxxxxxyz_0, \
                             g_0_zzzz_0_xxxxxxzz_0, \
                             g_0_zzzz_0_xxxxxyyy_0, \
                             g_0_zzzz_0_xxxxxyyz_0, \
                             g_0_zzzz_0_xxxxxyzz_0, \
                             g_0_zzzz_0_xxxxxzzz_0, \
                             g_0_zzzz_0_xxxxyyyy_0, \
                             g_0_zzzz_0_xxxxyyyz_0, \
                             g_0_zzzz_0_xxxxyyzz_0, \
                             g_0_zzzz_0_xxxxyzzz_0, \
                             g_0_zzzz_0_xxxxzzzz_0, \
                             g_0_zzzz_0_xxxyyyyy_0, \
                             g_0_zzzz_0_xxxyyyyz_0, \
                             g_0_zzzz_0_xxxyyyzz_0, \
                             g_0_zzzz_0_xxxyyzzz_0, \
                             g_0_zzzz_0_xxxyzzzz_0, \
                             g_0_zzzz_0_xxxzzzzz_0, \
                             g_0_zzzz_0_xxyyyyyy_0, \
                             g_0_zzzz_0_xxyyyyyz_0, \
                             g_0_zzzz_0_xxyyyyzz_0, \
                             g_0_zzzz_0_xxyyyzzz_0, \
                             g_0_zzzz_0_xxyyzzzz_0, \
                             g_0_zzzz_0_xxyzzzzz_0, \
                             g_0_zzzz_0_xxzzzzzz_0, \
                             g_0_zzzz_0_xyyyyyyy_0, \
                             g_0_zzzz_0_xyyyyyyz_0, \
                             g_0_zzzz_0_xyyyyyzz_0, \
                             g_0_zzzz_0_xyyyyzzz_0, \
                             g_0_zzzz_0_xyyyzzzz_0, \
                             g_0_zzzz_0_xyyzzzzz_0, \
                             g_0_zzzz_0_xyzzzzzz_0, \
                             g_0_zzzz_0_xzzzzzzz_0, \
                             g_0_zzzz_0_yyyyyyyy_0, \
                             g_0_zzzz_0_yyyyyyyz_0, \
                             g_0_zzzz_0_yyyyyyzz_0, \
                             g_0_zzzz_0_yyyyyzzz_0, \
                             g_0_zzzz_0_yyyyzzzz_0, \
                             g_0_zzzz_0_yyyzzzzz_0, \
                             g_0_zzzz_0_yyzzzzzz_0, \
                             g_0_zzzz_0_yzzzzzzz_0, \
                             g_0_zzzz_0_zzzzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzz_0_xxxxxxxx_0[i] = 3.0 * g_0_zz_0_xxxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxxxx_0[i] * pb_z +
                                   g_0_zzz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxxxy_0[i] = 3.0 * g_0_zz_0_xxxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxxxy_0[i] * pb_z +
                                   g_0_zzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxxxz_0[i] = 3.0 * g_0_zz_0_xxxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxxz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxxz_0[i] * pb_z + g_0_zzz_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxxyy_0[i] = 3.0 * g_0_zz_0_xxxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxxyy_0[i] * pb_z +
                                   g_0_zzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxxyz_0[i] = 3.0 * g_0_zz_0_xxxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxyz_0[i] * pb_z + g_0_zzz_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxxzz_0[i] = 3.0 * g_0_zz_0_xxxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxzz_0[i] * pb_z + g_0_zzz_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxyyy_0[i] = 3.0 * g_0_zz_0_xxxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxyyy_0[i] * pb_z +
                                   g_0_zzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxyyz_0[i] = 3.0 * g_0_zz_0_xxxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyyz_0[i] * pb_z + g_0_zzz_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxyzz_0[i] = 3.0 * g_0_zz_0_xxxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyzz_0[i] * pb_z + g_0_zzz_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxzzz_0[i] = 3.0 * g_0_zz_0_xxxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxzzz_0[i] * pb_z + g_0_zzz_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyyyy_0[i] = 3.0 * g_0_zz_0_xxxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxyyyy_0[i] * pb_z +
                                   g_0_zzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyyyz_0[i] = 3.0 * g_0_zz_0_xxxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyyz_0[i] * pb_z + g_0_zzz_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyyzz_0[i] = 3.0 * g_0_zz_0_xxxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyzz_0[i] * pb_z + g_0_zzz_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyzzz_0[i] = 3.0 * g_0_zz_0_xxxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyzzz_0[i] * pb_z + g_0_zzz_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxzzzz_0[i] = 3.0 * g_0_zz_0_xxxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxzzzz_0[i] * pb_z + g_0_zzz_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyyyy_0[i] = 3.0 * g_0_zz_0_xxxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxyyyyy_0[i] * pb_z +
                                   g_0_zzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyyyz_0[i] = 3.0 * g_0_zz_0_xxxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyyz_0[i] * pb_z + g_0_zzz_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyyzz_0[i] = 3.0 * g_0_zz_0_xxxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyzz_0[i] * pb_z + g_0_zzz_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyzzz_0[i] = 3.0 * g_0_zz_0_xxxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyzzz_0[i] * pb_z + g_0_zzz_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyzzzz_0[i] = 3.0 * g_0_zz_0_xxxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzzzz_0[i] * pb_z + g_0_zzz_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxzzzzz_0[i] = 3.0 * g_0_zz_0_xxxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxzzzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_zzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxzzzzz_0[i] * pb_z + g_0_zzz_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyyyy_0[i] = 3.0 * g_0_zz_0_xxyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxyyyyyy_0[i] * pb_z +
                                   g_0_zzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyyyz_0[i] = 3.0 * g_0_zz_0_xxyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyyz_0[i] * pb_z + g_0_zzz_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyyzz_0[i] = 3.0 * g_0_zz_0_xxyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyzz_0[i] * pb_z + g_0_zzz_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyzzz_0[i] = 3.0 * g_0_zz_0_xxyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyzzz_0[i] * pb_z + g_0_zzz_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyzzzz_0[i] = 3.0 * g_0_zz_0_xxyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzzzz_0[i] * pb_z + g_0_zzz_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyzzzzz_0[i] = 3.0 * g_0_zz_0_xxyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_zzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzzzz_0[i] * pb_z + g_0_zzz_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxzzzzzz_0[i] = 3.0 * g_0_zz_0_xxzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxzzzzzz_1[i] * fti_ab_0 +
                                   6.0 * g_0_zzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzzzzzz_0[i] * pb_z + g_0_zzz_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyyyy_0[i] = 3.0 * g_0_zz_0_xyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xyyyyyyy_0[i] * pb_z +
                                   g_0_zzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyyyz_0[i] = 3.0 * g_0_zz_0_xyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyyyz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyyz_0[i] * pb_z + g_0_zzz_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyyzz_0[i] = 3.0 * g_0_zz_0_xyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyzz_0[i] * pb_z + g_0_zzz_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyzzz_0[i] = 3.0 * g_0_zz_0_xyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyzzz_0[i] * pb_z + g_0_zzz_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyzzzz_0[i] = 3.0 * g_0_zz_0_xyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzzzz_0[i] * pb_z + g_0_zzz_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyzzzzz_0[i] = 3.0 * g_0_zz_0_xyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyzzzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_zzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzzzz_0[i] * pb_z + g_0_zzz_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyzzzzzz_0[i] = 3.0 * g_0_zz_0_xyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyzzzzzz_1[i] * fti_ab_0 +
                                   6.0 * g_0_zzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzzzz_0[i] * pb_z + g_0_zzz_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xzzzzzzz_0[i] = 3.0 * g_0_zz_0_xzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xzzzzzzz_1[i] * fti_ab_0 +
                                   7.0 * g_0_zzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzzzzzz_0[i] * pb_z + g_0_zzz_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyyyy_0[i] = 3.0 * g_0_zz_0_yyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_yyyyyyyy_0[i] * pb_z +
                                   g_0_zzz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyyyz_0[i] = 3.0 * g_0_zz_0_yyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyyyz_1[i] * fti_ab_0 +
                                   g_0_zzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyyz_0[i] * pb_z + g_0_zzz_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyyzz_0[i] = 3.0 * g_0_zz_0_yyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyzz_0[i] * pb_z + g_0_zzz_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyzzz_0[i] = 3.0 * g_0_zz_0_yyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyzzz_0[i] * pb_z + g_0_zzz_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyzzzz_0[i] = 3.0 * g_0_zz_0_yyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyzzzz_0[i] * pb_z + g_0_zzz_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyzzzzz_0[i] = 3.0 * g_0_zz_0_yyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyzzzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_zzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyzzzzz_0[i] * pb_z + g_0_zzz_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyzzzzzz_0[i] = 3.0 * g_0_zz_0_yyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyzzzzzz_1[i] * fti_ab_0 +
                                   6.0 * g_0_zzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzzzzzz_0[i] * pb_z + g_0_zzz_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yzzzzzzz_0[i] = 3.0 * g_0_zz_0_yzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yzzzzzzz_1[i] * fti_ab_0 +
                                   7.0 * g_0_zzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzzzzzz_0[i] * pb_z + g_0_zzz_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_zzzzzzzz_0[i] = 3.0 * g_0_zz_0_zzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_zzzzzzzz_1[i] * fti_ab_0 +
                                   8.0 * g_0_zzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_zzzzzzzz_0[i] * pb_z + g_0_zzz_0_zzzzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
