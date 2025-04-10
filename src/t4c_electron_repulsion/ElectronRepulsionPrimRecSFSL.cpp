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

#include "ElectronRepulsionPrimRecSFSL.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sfsl(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sfsl,
                                  size_t                idx_eri_0_spsl,
                                  size_t                idx_eri_1_spsl,
                                  size_t                idx_eri_1_sdsk,
                                  size_t                idx_eri_0_sdsl,
                                  size_t                idx_eri_1_sdsl,
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

    /// Set up components of auxilary buffer : SPSL

    auto g_0_x_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_spsl);

    auto g_0_x_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_spsl + 1);

    auto g_0_x_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_spsl + 2);

    auto g_0_x_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_spsl + 3);

    auto g_0_x_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_spsl + 4);

    auto g_0_x_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_spsl + 5);

    auto g_0_x_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_spsl + 6);

    auto g_0_x_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_spsl + 7);

    auto g_0_x_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_spsl + 8);

    auto g_0_x_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_spsl + 9);

    auto g_0_x_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_spsl + 10);

    auto g_0_x_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_spsl + 11);

    auto g_0_x_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_spsl + 12);

    auto g_0_x_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_spsl + 13);

    auto g_0_x_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_spsl + 14);

    auto g_0_x_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 15);

    auto g_0_x_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 16);

    auto g_0_x_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 17);

    auto g_0_x_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 18);

    auto g_0_x_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 19);

    auto g_0_x_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 20);

    auto g_0_x_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 21);

    auto g_0_x_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 22);

    auto g_0_x_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 23);

    auto g_0_x_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 24);

    auto g_0_x_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 25);

    auto g_0_x_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 26);

    auto g_0_x_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 27);

    auto g_0_x_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 28);

    auto g_0_x_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 29);

    auto g_0_x_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 30);

    auto g_0_x_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 31);

    auto g_0_x_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 32);

    auto g_0_x_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 33);

    auto g_0_x_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 34);

    auto g_0_x_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 35);

    auto g_0_x_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 36);

    auto g_0_x_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 37);

    auto g_0_x_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 38);

    auto g_0_x_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 39);

    auto g_0_x_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 40);

    auto g_0_x_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 41);

    auto g_0_x_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 42);

    auto g_0_x_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 43);

    auto g_0_x_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 44);

    auto g_0_y_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_spsl + 45);

    auto g_0_y_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_spsl + 46);

    auto g_0_y_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_spsl + 47);

    auto g_0_y_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_spsl + 48);

    auto g_0_y_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_spsl + 49);

    auto g_0_y_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_spsl + 50);

    auto g_0_y_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_spsl + 51);

    auto g_0_y_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_spsl + 52);

    auto g_0_y_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_spsl + 53);

    auto g_0_y_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_spsl + 54);

    auto g_0_y_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_spsl + 55);

    auto g_0_y_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_spsl + 56);

    auto g_0_y_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_spsl + 57);

    auto g_0_y_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_spsl + 58);

    auto g_0_y_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_spsl + 59);

    auto g_0_y_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 60);

    auto g_0_y_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 61);

    auto g_0_y_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 62);

    auto g_0_y_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 63);

    auto g_0_y_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 64);

    auto g_0_y_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 65);

    auto g_0_y_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 66);

    auto g_0_y_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 67);

    auto g_0_y_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 68);

    auto g_0_y_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 69);

    auto g_0_y_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 70);

    auto g_0_y_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 71);

    auto g_0_y_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 72);

    auto g_0_y_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 73);

    auto g_0_y_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 74);

    auto g_0_y_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 75);

    auto g_0_y_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 76);

    auto g_0_y_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 77);

    auto g_0_y_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 78);

    auto g_0_y_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 79);

    auto g_0_y_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 80);

    auto g_0_y_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 81);

    auto g_0_y_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 82);

    auto g_0_y_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 83);

    auto g_0_y_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 84);

    auto g_0_y_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 85);

    auto g_0_y_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 86);

    auto g_0_y_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 87);

    auto g_0_y_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 88);

    auto g_0_y_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 89);

    auto g_0_z_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_spsl + 90);

    auto g_0_z_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_spsl + 91);

    auto g_0_z_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_spsl + 92);

    auto g_0_z_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_spsl + 93);

    auto g_0_z_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_spsl + 94);

    auto g_0_z_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_spsl + 95);

    auto g_0_z_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_spsl + 96);

    auto g_0_z_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_spsl + 97);

    auto g_0_z_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_spsl + 98);

    auto g_0_z_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_spsl + 99);

    auto g_0_z_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_spsl + 100);

    auto g_0_z_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_spsl + 101);

    auto g_0_z_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_spsl + 102);

    auto g_0_z_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_spsl + 103);

    auto g_0_z_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_spsl + 104);

    auto g_0_z_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 105);

    auto g_0_z_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 106);

    auto g_0_z_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 107);

    auto g_0_z_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 108);

    auto g_0_z_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 109);

    auto g_0_z_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 110);

    auto g_0_z_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 111);

    auto g_0_z_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 112);

    auto g_0_z_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 113);

    auto g_0_z_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 114);

    auto g_0_z_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 115);

    auto g_0_z_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 116);

    auto g_0_z_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 117);

    auto g_0_z_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 118);

    auto g_0_z_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 119);

    auto g_0_z_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 120);

    auto g_0_z_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 121);

    auto g_0_z_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 122);

    auto g_0_z_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 123);

    auto g_0_z_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 124);

    auto g_0_z_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 125);

    auto g_0_z_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 126);

    auto g_0_z_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 127);

    auto g_0_z_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 128);

    auto g_0_z_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 129);

    auto g_0_z_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 130);

    auto g_0_z_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 131);

    auto g_0_z_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 132);

    auto g_0_z_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 133);

    auto g_0_z_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 134);

    /// Set up components of auxilary buffer : SPSL

    auto g_0_x_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_spsl);

    auto g_0_x_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_spsl + 1);

    auto g_0_x_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_spsl + 2);

    auto g_0_x_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_spsl + 3);

    auto g_0_x_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_spsl + 4);

    auto g_0_x_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_spsl + 5);

    auto g_0_x_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_spsl + 6);

    auto g_0_x_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_spsl + 7);

    auto g_0_x_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_spsl + 8);

    auto g_0_x_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_spsl + 9);

    auto g_0_x_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_spsl + 10);

    auto g_0_x_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_spsl + 11);

    auto g_0_x_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_spsl + 12);

    auto g_0_x_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_spsl + 13);

    auto g_0_x_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_spsl + 14);

    auto g_0_x_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 15);

    auto g_0_x_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 16);

    auto g_0_x_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 17);

    auto g_0_x_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 18);

    auto g_0_x_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 19);

    auto g_0_x_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 20);

    auto g_0_x_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 21);

    auto g_0_x_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 22);

    auto g_0_x_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 23);

    auto g_0_x_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 24);

    auto g_0_x_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 25);

    auto g_0_x_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 26);

    auto g_0_x_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 27);

    auto g_0_x_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 28);

    auto g_0_x_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 29);

    auto g_0_x_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 30);

    auto g_0_x_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 31);

    auto g_0_x_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 32);

    auto g_0_x_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 33);

    auto g_0_x_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 34);

    auto g_0_x_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 35);

    auto g_0_x_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 36);

    auto g_0_x_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 37);

    auto g_0_x_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 38);

    auto g_0_x_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 39);

    auto g_0_x_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 40);

    auto g_0_x_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 41);

    auto g_0_x_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 42);

    auto g_0_x_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 43);

    auto g_0_x_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 44);

    auto g_0_y_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_spsl + 45);

    auto g_0_y_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_spsl + 46);

    auto g_0_y_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_spsl + 47);

    auto g_0_y_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_spsl + 48);

    auto g_0_y_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_spsl + 49);

    auto g_0_y_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_spsl + 50);

    auto g_0_y_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_spsl + 51);

    auto g_0_y_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_spsl + 52);

    auto g_0_y_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_spsl + 53);

    auto g_0_y_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_spsl + 54);

    auto g_0_y_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_spsl + 55);

    auto g_0_y_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_spsl + 56);

    auto g_0_y_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_spsl + 57);

    auto g_0_y_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_spsl + 58);

    auto g_0_y_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_spsl + 59);

    auto g_0_y_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 60);

    auto g_0_y_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 61);

    auto g_0_y_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 62);

    auto g_0_y_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 63);

    auto g_0_y_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 64);

    auto g_0_y_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 65);

    auto g_0_y_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 66);

    auto g_0_y_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 67);

    auto g_0_y_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 68);

    auto g_0_y_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 69);

    auto g_0_y_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 70);

    auto g_0_y_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 71);

    auto g_0_y_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 72);

    auto g_0_y_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 73);

    auto g_0_y_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 74);

    auto g_0_y_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 75);

    auto g_0_y_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 76);

    auto g_0_y_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 77);

    auto g_0_y_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 78);

    auto g_0_y_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 79);

    auto g_0_y_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 80);

    auto g_0_y_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 81);

    auto g_0_y_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 82);

    auto g_0_y_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 83);

    auto g_0_y_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 84);

    auto g_0_y_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 85);

    auto g_0_y_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 86);

    auto g_0_y_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 87);

    auto g_0_y_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 88);

    auto g_0_y_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 89);

    auto g_0_z_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_spsl + 90);

    auto g_0_z_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_spsl + 91);

    auto g_0_z_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_spsl + 92);

    auto g_0_z_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_spsl + 93);

    auto g_0_z_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_spsl + 94);

    auto g_0_z_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_spsl + 95);

    auto g_0_z_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_spsl + 96);

    auto g_0_z_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_spsl + 97);

    auto g_0_z_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_spsl + 98);

    auto g_0_z_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_spsl + 99);

    auto g_0_z_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_spsl + 100);

    auto g_0_z_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_spsl + 101);

    auto g_0_z_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_spsl + 102);

    auto g_0_z_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_spsl + 103);

    auto g_0_z_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_spsl + 104);

    auto g_0_z_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 105);

    auto g_0_z_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 106);

    auto g_0_z_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 107);

    auto g_0_z_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 108);

    auto g_0_z_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 109);

    auto g_0_z_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 110);

    auto g_0_z_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 111);

    auto g_0_z_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 112);

    auto g_0_z_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 113);

    auto g_0_z_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 114);

    auto g_0_z_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 115);

    auto g_0_z_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 116);

    auto g_0_z_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 117);

    auto g_0_z_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 118);

    auto g_0_z_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 119);

    auto g_0_z_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 120);

    auto g_0_z_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 121);

    auto g_0_z_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 122);

    auto g_0_z_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 123);

    auto g_0_z_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 124);

    auto g_0_z_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 125);

    auto g_0_z_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_spsl + 126);

    auto g_0_z_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_spsl + 127);

    auto g_0_z_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_spsl + 128);

    auto g_0_z_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_spsl + 129);

    auto g_0_z_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_spsl + 130);

    auto g_0_z_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 131);

    auto g_0_z_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 132);

    auto g_0_z_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 133);

    auto g_0_z_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_spsl + 134);

    /// Set up components of auxilary buffer : SDSK

    auto g_0_xx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sdsk);

    auto g_0_xx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sdsk + 1);

    auto g_0_xx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sdsk + 2);

    auto g_0_xx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sdsk + 3);

    auto g_0_xx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sdsk + 4);

    auto g_0_xx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sdsk + 5);

    auto g_0_xx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sdsk + 6);

    auto g_0_xx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sdsk + 7);

    auto g_0_xx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sdsk + 8);

    auto g_0_xx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sdsk + 9);

    auto g_0_xx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 10);

    auto g_0_xx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 11);

    auto g_0_xx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 12);

    auto g_0_xx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 13);

    auto g_0_xx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 14);

    auto g_0_xx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 15);

    auto g_0_xx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 16);

    auto g_0_xx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 17);

    auto g_0_xx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 18);

    auto g_0_xx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 19);

    auto g_0_xx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 20);

    auto g_0_xx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 21);

    auto g_0_xx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 22);

    auto g_0_xx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 23);

    auto g_0_xx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 24);

    auto g_0_xx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 25);

    auto g_0_xx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 26);

    auto g_0_xx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 27);

    auto g_0_xx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 28);

    auto g_0_xx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 29);

    auto g_0_xx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 30);

    auto g_0_xx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 31);

    auto g_0_xx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 32);

    auto g_0_xx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 33);

    auto g_0_xx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 34);

    auto g_0_xx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 35);

    auto g_0_yy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sdsk + 108);

    auto g_0_yy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sdsk + 109);

    auto g_0_yy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sdsk + 110);

    auto g_0_yy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sdsk + 111);

    auto g_0_yy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sdsk + 112);

    auto g_0_yy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sdsk + 113);

    auto g_0_yy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sdsk + 114);

    auto g_0_yy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sdsk + 115);

    auto g_0_yy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sdsk + 116);

    auto g_0_yy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sdsk + 117);

    auto g_0_yy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 118);

    auto g_0_yy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 119);

    auto g_0_yy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 120);

    auto g_0_yy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 121);

    auto g_0_yy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 122);

    auto g_0_yy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 123);

    auto g_0_yy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 124);

    auto g_0_yy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 125);

    auto g_0_yy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 126);

    auto g_0_yy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 127);

    auto g_0_yy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 128);

    auto g_0_yy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 129);

    auto g_0_yy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 130);

    auto g_0_yy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 131);

    auto g_0_yy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 132);

    auto g_0_yy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 133);

    auto g_0_yy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 134);

    auto g_0_yy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 135);

    auto g_0_yy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 136);

    auto g_0_yy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 137);

    auto g_0_yy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 138);

    auto g_0_yy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 139);

    auto g_0_yy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 140);

    auto g_0_yy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 141);

    auto g_0_yy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 142);

    auto g_0_yy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 143);

    auto g_0_yz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sdsk + 148);

    auto g_0_yz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sdsk + 151);

    auto g_0_yz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sdsk + 152);

    auto g_0_yz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 155);

    auto g_0_yz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 156);

    auto g_0_yz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 157);

    auto g_0_yz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 160);

    auto g_0_yz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 161);

    auto g_0_yz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 162);

    auto g_0_yz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 163);

    auto g_0_yz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 166);

    auto g_0_yz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 167);

    auto g_0_yz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 168);

    auto g_0_yz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 169);

    auto g_0_yz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 170);

    auto g_0_yz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 173);

    auto g_0_yz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 174);

    auto g_0_yz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 175);

    auto g_0_yz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 176);

    auto g_0_yz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 177);

    auto g_0_yz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 178);

    auto g_0_zz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sdsk + 180);

    auto g_0_zz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sdsk + 181);

    auto g_0_zz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sdsk + 182);

    auto g_0_zz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sdsk + 183);

    auto g_0_zz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sdsk + 184);

    auto g_0_zz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sdsk + 185);

    auto g_0_zz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sdsk + 186);

    auto g_0_zz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sdsk + 187);

    auto g_0_zz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sdsk + 188);

    auto g_0_zz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sdsk + 189);

    auto g_0_zz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 190);

    auto g_0_zz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 191);

    auto g_0_zz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 192);

    auto g_0_zz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 193);

    auto g_0_zz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 194);

    auto g_0_zz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 195);

    auto g_0_zz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 196);

    auto g_0_zz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 197);

    auto g_0_zz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 198);

    auto g_0_zz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 199);

    auto g_0_zz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 200);

    auto g_0_zz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 201);

    auto g_0_zz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 202);

    auto g_0_zz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 203);

    auto g_0_zz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 204);

    auto g_0_zz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 205);

    auto g_0_zz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 206);

    auto g_0_zz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 207);

    auto g_0_zz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 208);

    auto g_0_zz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 209);

    auto g_0_zz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 210);

    auto g_0_zz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 211);

    auto g_0_zz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 212);

    auto g_0_zz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 213);

    auto g_0_zz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 214);

    auto g_0_zz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 215);

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

    auto g_0_xy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sdsl + 46);

    auto g_0_xy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sdsl + 48);

    auto g_0_xy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sdsl + 51);

    auto g_0_xy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 55);

    auto g_0_xy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 60);

    auto g_0_xy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 66);

    auto g_0_xy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 73);

    auto g_0_xz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sdsl + 90);

    auto g_0_xz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sdsl + 92);

    auto g_0_xz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sdsl + 95);

    auto g_0_xz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sdsl + 99);

    auto g_0_xz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 104);

    auto g_0_xz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 110);

    auto g_0_xz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 117);

    auto g_0_xz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 125);

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

    auto g_0_yz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sdsl + 184);

    auto g_0_yz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sdsl + 187);

    auto g_0_yz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sdsl + 188);

    auto g_0_yz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 191);

    auto g_0_yz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 192);

    auto g_0_yz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 193);

    auto g_0_yz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 196);

    auto g_0_yz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 197);

    auto g_0_yz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 198);

    auto g_0_yz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 199);

    auto g_0_yz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 202);

    auto g_0_yz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 203);

    auto g_0_yz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 204);

    auto g_0_yz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 205);

    auto g_0_yz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 206);

    auto g_0_yz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 209);

    auto g_0_yz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 210);

    auto g_0_yz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 211);

    auto g_0_yz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 212);

    auto g_0_yz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 213);

    auto g_0_yz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 214);

    auto g_0_yz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sdsl + 216);

    auto g_0_yz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sdsl + 217);

    auto g_0_yz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sdsl + 218);

    auto g_0_yz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sdsl + 219);

    auto g_0_yz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 220);

    auto g_0_yz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 221);

    auto g_0_yz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 222);

    auto g_0_yz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 223);

    auto g_0_yz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sdsl + 224);

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

    auto g_0_xy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sdsl + 46);

    auto g_0_xy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sdsl + 48);

    auto g_0_xy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sdsl + 51);

    auto g_0_xy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 55);

    auto g_0_xy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 60);

    auto g_0_xy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 66);

    auto g_0_xy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 73);

    auto g_0_xz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sdsl + 90);

    auto g_0_xz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sdsl + 92);

    auto g_0_xz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sdsl + 95);

    auto g_0_xz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sdsl + 99);

    auto g_0_xz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 104);

    auto g_0_xz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 110);

    auto g_0_xz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 117);

    auto g_0_xz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 125);

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

    auto g_0_yz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sdsl + 184);

    auto g_0_yz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sdsl + 187);

    auto g_0_yz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sdsl + 188);

    auto g_0_yz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 191);

    auto g_0_yz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 192);

    auto g_0_yz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 193);

    auto g_0_yz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 196);

    auto g_0_yz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 197);

    auto g_0_yz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 198);

    auto g_0_yz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 199);

    auto g_0_yz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 202);

    auto g_0_yz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 203);

    auto g_0_yz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 204);

    auto g_0_yz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 205);

    auto g_0_yz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 206);

    auto g_0_yz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 209);

    auto g_0_yz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 210);

    auto g_0_yz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 211);

    auto g_0_yz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 212);

    auto g_0_yz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 213);

    auto g_0_yz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 214);

    auto g_0_yz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sdsl + 216);

    auto g_0_yz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sdsl + 217);

    auto g_0_yz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sdsl + 218);

    auto g_0_yz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sdsl + 219);

    auto g_0_yz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 220);

    auto g_0_yz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 221);

    auto g_0_yz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 222);

    auto g_0_yz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 223);

    auto g_0_yz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sdsl + 224);

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

    /// Set up 0-45 components of targeted buffer : SFSL

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

#pragma omp simd aligned(g_0_x_0_xxxxxxxx_0,       \
                             g_0_x_0_xxxxxxxx_1,   \
                             g_0_x_0_xxxxxxxy_0,   \
                             g_0_x_0_xxxxxxxy_1,   \
                             g_0_x_0_xxxxxxxz_0,   \
                             g_0_x_0_xxxxxxxz_1,   \
                             g_0_x_0_xxxxxxyy_0,   \
                             g_0_x_0_xxxxxxyy_1,   \
                             g_0_x_0_xxxxxxyz_0,   \
                             g_0_x_0_xxxxxxyz_1,   \
                             g_0_x_0_xxxxxxzz_0,   \
                             g_0_x_0_xxxxxxzz_1,   \
                             g_0_x_0_xxxxxyyy_0,   \
                             g_0_x_0_xxxxxyyy_1,   \
                             g_0_x_0_xxxxxyyz_0,   \
                             g_0_x_0_xxxxxyyz_1,   \
                             g_0_x_0_xxxxxyzz_0,   \
                             g_0_x_0_xxxxxyzz_1,   \
                             g_0_x_0_xxxxxzzz_0,   \
                             g_0_x_0_xxxxxzzz_1,   \
                             g_0_x_0_xxxxyyyy_0,   \
                             g_0_x_0_xxxxyyyy_1,   \
                             g_0_x_0_xxxxyyyz_0,   \
                             g_0_x_0_xxxxyyyz_1,   \
                             g_0_x_0_xxxxyyzz_0,   \
                             g_0_x_0_xxxxyyzz_1,   \
                             g_0_x_0_xxxxyzzz_0,   \
                             g_0_x_0_xxxxyzzz_1,   \
                             g_0_x_0_xxxxzzzz_0,   \
                             g_0_x_0_xxxxzzzz_1,   \
                             g_0_x_0_xxxyyyyy_0,   \
                             g_0_x_0_xxxyyyyy_1,   \
                             g_0_x_0_xxxyyyyz_0,   \
                             g_0_x_0_xxxyyyyz_1,   \
                             g_0_x_0_xxxyyyzz_0,   \
                             g_0_x_0_xxxyyyzz_1,   \
                             g_0_x_0_xxxyyzzz_0,   \
                             g_0_x_0_xxxyyzzz_1,   \
                             g_0_x_0_xxxyzzzz_0,   \
                             g_0_x_0_xxxyzzzz_1,   \
                             g_0_x_0_xxxzzzzz_0,   \
                             g_0_x_0_xxxzzzzz_1,   \
                             g_0_x_0_xxyyyyyy_0,   \
                             g_0_x_0_xxyyyyyy_1,   \
                             g_0_x_0_xxyyyyyz_0,   \
                             g_0_x_0_xxyyyyyz_1,   \
                             g_0_x_0_xxyyyyzz_0,   \
                             g_0_x_0_xxyyyyzz_1,   \
                             g_0_x_0_xxyyyzzz_0,   \
                             g_0_x_0_xxyyyzzz_1,   \
                             g_0_x_0_xxyyzzzz_0,   \
                             g_0_x_0_xxyyzzzz_1,   \
                             g_0_x_0_xxyzzzzz_0,   \
                             g_0_x_0_xxyzzzzz_1,   \
                             g_0_x_0_xxzzzzzz_0,   \
                             g_0_x_0_xxzzzzzz_1,   \
                             g_0_x_0_xyyyyyyy_0,   \
                             g_0_x_0_xyyyyyyy_1,   \
                             g_0_x_0_xyyyyyyz_0,   \
                             g_0_x_0_xyyyyyyz_1,   \
                             g_0_x_0_xyyyyyzz_0,   \
                             g_0_x_0_xyyyyyzz_1,   \
                             g_0_x_0_xyyyyzzz_0,   \
                             g_0_x_0_xyyyyzzz_1,   \
                             g_0_x_0_xyyyzzzz_0,   \
                             g_0_x_0_xyyyzzzz_1,   \
                             g_0_x_0_xyyzzzzz_0,   \
                             g_0_x_0_xyyzzzzz_1,   \
                             g_0_x_0_xyzzzzzz_0,   \
                             g_0_x_0_xyzzzzzz_1,   \
                             g_0_x_0_xzzzzzzz_0,   \
                             g_0_x_0_xzzzzzzz_1,   \
                             g_0_x_0_yyyyyyyy_0,   \
                             g_0_x_0_yyyyyyyy_1,   \
                             g_0_x_0_yyyyyyyz_0,   \
                             g_0_x_0_yyyyyyyz_1,   \
                             g_0_x_0_yyyyyyzz_0,   \
                             g_0_x_0_yyyyyyzz_1,   \
                             g_0_x_0_yyyyyzzz_0,   \
                             g_0_x_0_yyyyyzzz_1,   \
                             g_0_x_0_yyyyzzzz_0,   \
                             g_0_x_0_yyyyzzzz_1,   \
                             g_0_x_0_yyyzzzzz_0,   \
                             g_0_x_0_yyyzzzzz_1,   \
                             g_0_x_0_yyzzzzzz_0,   \
                             g_0_x_0_yyzzzzzz_1,   \
                             g_0_x_0_yzzzzzzz_0,   \
                             g_0_x_0_yzzzzzzz_1,   \
                             g_0_x_0_zzzzzzzz_0,   \
                             g_0_x_0_zzzzzzzz_1,   \
                             g_0_xx_0_xxxxxxx_1,   \
                             g_0_xx_0_xxxxxxxx_0,  \
                             g_0_xx_0_xxxxxxxx_1,  \
                             g_0_xx_0_xxxxxxxy_0,  \
                             g_0_xx_0_xxxxxxxy_1,  \
                             g_0_xx_0_xxxxxxxz_0,  \
                             g_0_xx_0_xxxxxxxz_1,  \
                             g_0_xx_0_xxxxxxy_1,   \
                             g_0_xx_0_xxxxxxyy_0,  \
                             g_0_xx_0_xxxxxxyy_1,  \
                             g_0_xx_0_xxxxxxyz_0,  \
                             g_0_xx_0_xxxxxxyz_1,  \
                             g_0_xx_0_xxxxxxz_1,   \
                             g_0_xx_0_xxxxxxzz_0,  \
                             g_0_xx_0_xxxxxxzz_1,  \
                             g_0_xx_0_xxxxxyy_1,   \
                             g_0_xx_0_xxxxxyyy_0,  \
                             g_0_xx_0_xxxxxyyy_1,  \
                             g_0_xx_0_xxxxxyyz_0,  \
                             g_0_xx_0_xxxxxyyz_1,  \
                             g_0_xx_0_xxxxxyz_1,   \
                             g_0_xx_0_xxxxxyzz_0,  \
                             g_0_xx_0_xxxxxyzz_1,  \
                             g_0_xx_0_xxxxxzz_1,   \
                             g_0_xx_0_xxxxxzzz_0,  \
                             g_0_xx_0_xxxxxzzz_1,  \
                             g_0_xx_0_xxxxyyy_1,   \
                             g_0_xx_0_xxxxyyyy_0,  \
                             g_0_xx_0_xxxxyyyy_1,  \
                             g_0_xx_0_xxxxyyyz_0,  \
                             g_0_xx_0_xxxxyyyz_1,  \
                             g_0_xx_0_xxxxyyz_1,   \
                             g_0_xx_0_xxxxyyzz_0,  \
                             g_0_xx_0_xxxxyyzz_1,  \
                             g_0_xx_0_xxxxyzz_1,   \
                             g_0_xx_0_xxxxyzzz_0,  \
                             g_0_xx_0_xxxxyzzz_1,  \
                             g_0_xx_0_xxxxzzz_1,   \
                             g_0_xx_0_xxxxzzzz_0,  \
                             g_0_xx_0_xxxxzzzz_1,  \
                             g_0_xx_0_xxxyyyy_1,   \
                             g_0_xx_0_xxxyyyyy_0,  \
                             g_0_xx_0_xxxyyyyy_1,  \
                             g_0_xx_0_xxxyyyyz_0,  \
                             g_0_xx_0_xxxyyyyz_1,  \
                             g_0_xx_0_xxxyyyz_1,   \
                             g_0_xx_0_xxxyyyzz_0,  \
                             g_0_xx_0_xxxyyyzz_1,  \
                             g_0_xx_0_xxxyyzz_1,   \
                             g_0_xx_0_xxxyyzzz_0,  \
                             g_0_xx_0_xxxyyzzz_1,  \
                             g_0_xx_0_xxxyzzz_1,   \
                             g_0_xx_0_xxxyzzzz_0,  \
                             g_0_xx_0_xxxyzzzz_1,  \
                             g_0_xx_0_xxxzzzz_1,   \
                             g_0_xx_0_xxxzzzzz_0,  \
                             g_0_xx_0_xxxzzzzz_1,  \
                             g_0_xx_0_xxyyyyy_1,   \
                             g_0_xx_0_xxyyyyyy_0,  \
                             g_0_xx_0_xxyyyyyy_1,  \
                             g_0_xx_0_xxyyyyyz_0,  \
                             g_0_xx_0_xxyyyyyz_1,  \
                             g_0_xx_0_xxyyyyz_1,   \
                             g_0_xx_0_xxyyyyzz_0,  \
                             g_0_xx_0_xxyyyyzz_1,  \
                             g_0_xx_0_xxyyyzz_1,   \
                             g_0_xx_0_xxyyyzzz_0,  \
                             g_0_xx_0_xxyyyzzz_1,  \
                             g_0_xx_0_xxyyzzz_1,   \
                             g_0_xx_0_xxyyzzzz_0,  \
                             g_0_xx_0_xxyyzzzz_1,  \
                             g_0_xx_0_xxyzzzz_1,   \
                             g_0_xx_0_xxyzzzzz_0,  \
                             g_0_xx_0_xxyzzzzz_1,  \
                             g_0_xx_0_xxzzzzz_1,   \
                             g_0_xx_0_xxzzzzzz_0,  \
                             g_0_xx_0_xxzzzzzz_1,  \
                             g_0_xx_0_xyyyyyy_1,   \
                             g_0_xx_0_xyyyyyyy_0,  \
                             g_0_xx_0_xyyyyyyy_1,  \
                             g_0_xx_0_xyyyyyyz_0,  \
                             g_0_xx_0_xyyyyyyz_1,  \
                             g_0_xx_0_xyyyyyz_1,   \
                             g_0_xx_0_xyyyyyzz_0,  \
                             g_0_xx_0_xyyyyyzz_1,  \
                             g_0_xx_0_xyyyyzz_1,   \
                             g_0_xx_0_xyyyyzzz_0,  \
                             g_0_xx_0_xyyyyzzz_1,  \
                             g_0_xx_0_xyyyzzz_1,   \
                             g_0_xx_0_xyyyzzzz_0,  \
                             g_0_xx_0_xyyyzzzz_1,  \
                             g_0_xx_0_xyyzzzz_1,   \
                             g_0_xx_0_xyyzzzzz_0,  \
                             g_0_xx_0_xyyzzzzz_1,  \
                             g_0_xx_0_xyzzzzz_1,   \
                             g_0_xx_0_xyzzzzzz_0,  \
                             g_0_xx_0_xyzzzzzz_1,  \
                             g_0_xx_0_xzzzzzz_1,   \
                             g_0_xx_0_xzzzzzzz_0,  \
                             g_0_xx_0_xzzzzzzz_1,  \
                             g_0_xx_0_yyyyyyy_1,   \
                             g_0_xx_0_yyyyyyyy_0,  \
                             g_0_xx_0_yyyyyyyy_1,  \
                             g_0_xx_0_yyyyyyyz_0,  \
                             g_0_xx_0_yyyyyyyz_1,  \
                             g_0_xx_0_yyyyyyz_1,   \
                             g_0_xx_0_yyyyyyzz_0,  \
                             g_0_xx_0_yyyyyyzz_1,  \
                             g_0_xx_0_yyyyyzz_1,   \
                             g_0_xx_0_yyyyyzzz_0,  \
                             g_0_xx_0_yyyyyzzz_1,  \
                             g_0_xx_0_yyyyzzz_1,   \
                             g_0_xx_0_yyyyzzzz_0,  \
                             g_0_xx_0_yyyyzzzz_1,  \
                             g_0_xx_0_yyyzzzz_1,   \
                             g_0_xx_0_yyyzzzzz_0,  \
                             g_0_xx_0_yyyzzzzz_1,  \
                             g_0_xx_0_yyzzzzz_1,   \
                             g_0_xx_0_yyzzzzzz_0,  \
                             g_0_xx_0_yyzzzzzz_1,  \
                             g_0_xx_0_yzzzzzz_1,   \
                             g_0_xx_0_yzzzzzzz_0,  \
                             g_0_xx_0_yzzzzzzz_1,  \
                             g_0_xx_0_zzzzzzz_1,   \
                             g_0_xx_0_zzzzzzzz_0,  \
                             g_0_xx_0_zzzzzzzz_1,  \
                             g_0_xxx_0_xxxxxxxx_0, \
                             g_0_xxx_0_xxxxxxxy_0, \
                             g_0_xxx_0_xxxxxxxz_0, \
                             g_0_xxx_0_xxxxxxyy_0, \
                             g_0_xxx_0_xxxxxxyz_0, \
                             g_0_xxx_0_xxxxxxzz_0, \
                             g_0_xxx_0_xxxxxyyy_0, \
                             g_0_xxx_0_xxxxxyyz_0, \
                             g_0_xxx_0_xxxxxyzz_0, \
                             g_0_xxx_0_xxxxxzzz_0, \
                             g_0_xxx_0_xxxxyyyy_0, \
                             g_0_xxx_0_xxxxyyyz_0, \
                             g_0_xxx_0_xxxxyyzz_0, \
                             g_0_xxx_0_xxxxyzzz_0, \
                             g_0_xxx_0_xxxxzzzz_0, \
                             g_0_xxx_0_xxxyyyyy_0, \
                             g_0_xxx_0_xxxyyyyz_0, \
                             g_0_xxx_0_xxxyyyzz_0, \
                             g_0_xxx_0_xxxyyzzz_0, \
                             g_0_xxx_0_xxxyzzzz_0, \
                             g_0_xxx_0_xxxzzzzz_0, \
                             g_0_xxx_0_xxyyyyyy_0, \
                             g_0_xxx_0_xxyyyyyz_0, \
                             g_0_xxx_0_xxyyyyzz_0, \
                             g_0_xxx_0_xxyyyzzz_0, \
                             g_0_xxx_0_xxyyzzzz_0, \
                             g_0_xxx_0_xxyzzzzz_0, \
                             g_0_xxx_0_xxzzzzzz_0, \
                             g_0_xxx_0_xyyyyyyy_0, \
                             g_0_xxx_0_xyyyyyyz_0, \
                             g_0_xxx_0_xyyyyyzz_0, \
                             g_0_xxx_0_xyyyyzzz_0, \
                             g_0_xxx_0_xyyyzzzz_0, \
                             g_0_xxx_0_xyyzzzzz_0, \
                             g_0_xxx_0_xyzzzzzz_0, \
                             g_0_xxx_0_xzzzzzzz_0, \
                             g_0_xxx_0_yyyyyyyy_0, \
                             g_0_xxx_0_yyyyyyyz_0, \
                             g_0_xxx_0_yyyyyyzz_0, \
                             g_0_xxx_0_yyyyyzzz_0, \
                             g_0_xxx_0_yyyyzzzz_0, \
                             g_0_xxx_0_yyyzzzzz_0, \
                             g_0_xxx_0_yyzzzzzz_0, \
                             g_0_xxx_0_yzzzzzzz_0, \
                             g_0_xxx_0_zzzzzzzz_0, \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xxxxxxxx_0[i] = 2.0 * g_0_x_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxxx_1[i] * fti_ab_0 +
                                  8.0 * g_0_xx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxxx_0[i] * pb_x + g_0_xx_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxxxy_0[i] = 2.0 * g_0_x_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxxy_1[i] * fti_ab_0 +
                                  7.0 * g_0_xx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxxy_0[i] * pb_x + g_0_xx_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxxxz_0[i] = 2.0 * g_0_x_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxxz_1[i] * fti_ab_0 +
                                  7.0 * g_0_xx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxxz_0[i] * pb_x + g_0_xx_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxxyy_0[i] = 2.0 * g_0_x_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxyy_1[i] * fti_ab_0 +
                                  6.0 * g_0_xx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxyy_0[i] * pb_x + g_0_xx_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxxyz_0[i] = 2.0 * g_0_x_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxyz_1[i] * fti_ab_0 +
                                  6.0 * g_0_xx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxyz_0[i] * pb_x + g_0_xx_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxxzz_0[i] = 2.0 * g_0_x_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxzz_1[i] * fti_ab_0 +
                                  6.0 * g_0_xx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxzz_0[i] * pb_x + g_0_xx_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxyyy_0[i] = 2.0 * g_0_x_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxyyy_1[i] * fti_ab_0 +
                                  5.0 * g_0_xx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyyy_0[i] * pb_x + g_0_xx_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxyyz_0[i] = 2.0 * g_0_x_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxyyz_1[i] * fti_ab_0 +
                                  5.0 * g_0_xx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyyz_0[i] * pb_x + g_0_xx_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxyzz_0[i] = 2.0 * g_0_x_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxyzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_xx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyzz_0[i] * pb_x + g_0_xx_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxzzz_0[i] = 2.0 * g_0_x_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_xx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxzzz_0[i] * pb_x + g_0_xx_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyyyy_0[i] = 2.0 * g_0_x_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyyyy_1[i] * fti_ab_0 +
                                  4.0 * g_0_xx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyyy_0[i] * pb_x + g_0_xx_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyyyz_0[i] = 2.0 * g_0_x_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyyyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyyz_0[i] * pb_x + g_0_xx_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyyzz_0[i] = 2.0 * g_0_x_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyyzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyzz_0[i] * pb_x + g_0_xx_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyzzz_0[i] = 2.0 * g_0_x_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyzzz_0[i] * pb_x + g_0_xx_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxzzzz_0[i] = 2.0 * g_0_x_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxzzzz_0[i] * pb_x + g_0_xx_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyyyy_0[i] = 2.0 * g_0_x_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyyyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyyy_0[i] * pb_x + g_0_xx_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyyyz_0[i] = 2.0 * g_0_x_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyyz_0[i] * pb_x + g_0_xx_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyyzz_0[i] = 2.0 * g_0_x_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyzz_0[i] * pb_x + g_0_xx_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyzzz_0[i] = 2.0 * g_0_x_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyzzz_0[i] * pb_x + g_0_xx_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyzzzz_0[i] = 2.0 * g_0_x_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyzzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzzzz_0[i] * pb_x + g_0_xx_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxzzzzz_0[i] = 2.0 * g_0_x_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxzzzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxzzzzz_0[i] * pb_x + g_0_xx_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyyyy_0[i] = 2.0 * g_0_x_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyyyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyyy_0[i] * pb_x + g_0_xx_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyyyz_0[i] = 2.0 * g_0_x_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyyz_0[i] * pb_x + g_0_xx_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyyzz_0[i] = 2.0 * g_0_x_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyzz_0[i] * pb_x + g_0_xx_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyzzz_0[i] = 2.0 * g_0_x_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyzzz_0[i] * pb_x + g_0_xx_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyzzzz_0[i] = 2.0 * g_0_x_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzzzz_0[i] * pb_x + g_0_xx_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyzzzzz_0[i] = 2.0 * g_0_x_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyzzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzzzz_0[i] * pb_x + g_0_xx_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxzzzzzz_0[i] = 2.0 * g_0_x_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxzzzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzzzzzz_0[i] * pb_x + g_0_xx_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyyyy_0[i] = 2.0 * g_0_x_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyyyyy_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xyyyyyyy_0[i] * pb_x + g_0_xx_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyyyz_0[i] = 2.0 * g_0_x_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyyz_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xyyyyyyz_0[i] * pb_x + g_0_xx_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyyzz_0[i] = 2.0 * g_0_x_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyzz_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xyyyyyzz_0[i] * pb_x + g_0_xx_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyzzz_0[i] = 2.0 * g_0_x_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyzzz_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xyyyyzzz_0[i] * pb_x + g_0_xx_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyzzzz_0[i] = 2.0 * g_0_x_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyzzzz_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xyyyzzzz_0[i] * pb_x + g_0_xx_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyzzzzz_0[i] = 2.0 * g_0_x_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyzzzzz_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xyyzzzzz_0[i] * pb_x + g_0_xx_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyzzzzzz_0[i] = 2.0 * g_0_x_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzzzzz_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xyzzzzzz_0[i] * pb_x + g_0_xx_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xzzzzzzz_0[i] = 2.0 * g_0_x_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzzzzz_1[i] * fi_abcd_0 +
                                  g_0_xx_0_xzzzzzzz_0[i] * pb_x + g_0_xx_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyyyy_0[i] = 2.0 * g_0_x_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyyyyyy_0[i] * pb_x +
                                  g_0_xx_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyyyz_0[i] = 2.0 * g_0_x_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyyyz_0[i] * pb_x +
                                  g_0_xx_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyyzz_0[i] = 2.0 * g_0_x_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyyzz_0[i] * pb_x +
                                  g_0_xx_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyzzz_0[i] = 2.0 * g_0_x_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyzzz_0[i] * pb_x +
                                  g_0_xx_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyzzzz_0[i] = 2.0 * g_0_x_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyzzzz_0[i] * pb_x +
                                  g_0_xx_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyzzzzz_0[i] = 2.0 * g_0_x_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyzzzzz_0[i] * pb_x +
                                  g_0_xx_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyzzzzzz_0[i] = 2.0 * g_0_x_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyzzzzzz_0[i] * pb_x +
                                  g_0_xx_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yzzzzzzz_0[i] = 2.0 * g_0_x_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzzzzzz_0[i] * pb_x +
                                  g_0_xx_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_zzzzzzzz_0[i] = 2.0 * g_0_x_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzzzzzz_0[i] * pb_x +
                                  g_0_xx_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 45-90 components of targeted buffer : SFSL

    auto g_0_xxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 45);

    auto g_0_xxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 46);

    auto g_0_xxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 47);

    auto g_0_xxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 48);

    auto g_0_xxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 49);

    auto g_0_xxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 50);

    auto g_0_xxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 51);

    auto g_0_xxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 52);

    auto g_0_xxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 53);

    auto g_0_xxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 54);

    auto g_0_xxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 55);

    auto g_0_xxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 56);

    auto g_0_xxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 57);

    auto g_0_xxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 58);

    auto g_0_xxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 59);

    auto g_0_xxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 60);

    auto g_0_xxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 61);

    auto g_0_xxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 62);

    auto g_0_xxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 63);

    auto g_0_xxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 64);

    auto g_0_xxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 65);

    auto g_0_xxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 66);

    auto g_0_xxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 67);

    auto g_0_xxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 68);

    auto g_0_xxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 69);

    auto g_0_xxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 70);

    auto g_0_xxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 71);

    auto g_0_xxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 72);

    auto g_0_xxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 73);

    auto g_0_xxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 74);

    auto g_0_xxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 75);

    auto g_0_xxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 76);

    auto g_0_xxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 77);

    auto g_0_xxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 78);

    auto g_0_xxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 79);

    auto g_0_xxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 80);

    auto g_0_xxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 81);

    auto g_0_xxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 82);

    auto g_0_xxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 83);

    auto g_0_xxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 84);

    auto g_0_xxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 85);

    auto g_0_xxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 86);

    auto g_0_xxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 87);

    auto g_0_xxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 88);

    auto g_0_xxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 89);

#pragma omp simd aligned(g_0_xx_0_xxxxxxx_1,       \
                             g_0_xx_0_xxxxxxxx_0,  \
                             g_0_xx_0_xxxxxxxx_1,  \
                             g_0_xx_0_xxxxxxxy_0,  \
                             g_0_xx_0_xxxxxxxy_1,  \
                             g_0_xx_0_xxxxxxxz_0,  \
                             g_0_xx_0_xxxxxxxz_1,  \
                             g_0_xx_0_xxxxxxy_1,   \
                             g_0_xx_0_xxxxxxyy_0,  \
                             g_0_xx_0_xxxxxxyy_1,  \
                             g_0_xx_0_xxxxxxyz_0,  \
                             g_0_xx_0_xxxxxxyz_1,  \
                             g_0_xx_0_xxxxxxz_1,   \
                             g_0_xx_0_xxxxxxzz_0,  \
                             g_0_xx_0_xxxxxxzz_1,  \
                             g_0_xx_0_xxxxxyy_1,   \
                             g_0_xx_0_xxxxxyyy_0,  \
                             g_0_xx_0_xxxxxyyy_1,  \
                             g_0_xx_0_xxxxxyyz_0,  \
                             g_0_xx_0_xxxxxyyz_1,  \
                             g_0_xx_0_xxxxxyz_1,   \
                             g_0_xx_0_xxxxxyzz_0,  \
                             g_0_xx_0_xxxxxyzz_1,  \
                             g_0_xx_0_xxxxxzz_1,   \
                             g_0_xx_0_xxxxxzzz_0,  \
                             g_0_xx_0_xxxxxzzz_1,  \
                             g_0_xx_0_xxxxyyy_1,   \
                             g_0_xx_0_xxxxyyyy_0,  \
                             g_0_xx_0_xxxxyyyy_1,  \
                             g_0_xx_0_xxxxyyyz_0,  \
                             g_0_xx_0_xxxxyyyz_1,  \
                             g_0_xx_0_xxxxyyz_1,   \
                             g_0_xx_0_xxxxyyzz_0,  \
                             g_0_xx_0_xxxxyyzz_1,  \
                             g_0_xx_0_xxxxyzz_1,   \
                             g_0_xx_0_xxxxyzzz_0,  \
                             g_0_xx_0_xxxxyzzz_1,  \
                             g_0_xx_0_xxxxzzz_1,   \
                             g_0_xx_0_xxxxzzzz_0,  \
                             g_0_xx_0_xxxxzzzz_1,  \
                             g_0_xx_0_xxxyyyy_1,   \
                             g_0_xx_0_xxxyyyyy_0,  \
                             g_0_xx_0_xxxyyyyy_1,  \
                             g_0_xx_0_xxxyyyyz_0,  \
                             g_0_xx_0_xxxyyyyz_1,  \
                             g_0_xx_0_xxxyyyz_1,   \
                             g_0_xx_0_xxxyyyzz_0,  \
                             g_0_xx_0_xxxyyyzz_1,  \
                             g_0_xx_0_xxxyyzz_1,   \
                             g_0_xx_0_xxxyyzzz_0,  \
                             g_0_xx_0_xxxyyzzz_1,  \
                             g_0_xx_0_xxxyzzz_1,   \
                             g_0_xx_0_xxxyzzzz_0,  \
                             g_0_xx_0_xxxyzzzz_1,  \
                             g_0_xx_0_xxxzzzz_1,   \
                             g_0_xx_0_xxxzzzzz_0,  \
                             g_0_xx_0_xxxzzzzz_1,  \
                             g_0_xx_0_xxyyyyy_1,   \
                             g_0_xx_0_xxyyyyyy_0,  \
                             g_0_xx_0_xxyyyyyy_1,  \
                             g_0_xx_0_xxyyyyyz_0,  \
                             g_0_xx_0_xxyyyyyz_1,  \
                             g_0_xx_0_xxyyyyz_1,   \
                             g_0_xx_0_xxyyyyzz_0,  \
                             g_0_xx_0_xxyyyyzz_1,  \
                             g_0_xx_0_xxyyyzz_1,   \
                             g_0_xx_0_xxyyyzzz_0,  \
                             g_0_xx_0_xxyyyzzz_1,  \
                             g_0_xx_0_xxyyzzz_1,   \
                             g_0_xx_0_xxyyzzzz_0,  \
                             g_0_xx_0_xxyyzzzz_1,  \
                             g_0_xx_0_xxyzzzz_1,   \
                             g_0_xx_0_xxyzzzzz_0,  \
                             g_0_xx_0_xxyzzzzz_1,  \
                             g_0_xx_0_xxzzzzz_1,   \
                             g_0_xx_0_xxzzzzzz_0,  \
                             g_0_xx_0_xxzzzzzz_1,  \
                             g_0_xx_0_xyyyyyy_1,   \
                             g_0_xx_0_xyyyyyyy_0,  \
                             g_0_xx_0_xyyyyyyy_1,  \
                             g_0_xx_0_xyyyyyyz_0,  \
                             g_0_xx_0_xyyyyyyz_1,  \
                             g_0_xx_0_xyyyyyz_1,   \
                             g_0_xx_0_xyyyyyzz_0,  \
                             g_0_xx_0_xyyyyyzz_1,  \
                             g_0_xx_0_xyyyyzz_1,   \
                             g_0_xx_0_xyyyyzzz_0,  \
                             g_0_xx_0_xyyyyzzz_1,  \
                             g_0_xx_0_xyyyzzz_1,   \
                             g_0_xx_0_xyyyzzzz_0,  \
                             g_0_xx_0_xyyyzzzz_1,  \
                             g_0_xx_0_xyyzzzz_1,   \
                             g_0_xx_0_xyyzzzzz_0,  \
                             g_0_xx_0_xyyzzzzz_1,  \
                             g_0_xx_0_xyzzzzz_1,   \
                             g_0_xx_0_xyzzzzzz_0,  \
                             g_0_xx_0_xyzzzzzz_1,  \
                             g_0_xx_0_xzzzzzz_1,   \
                             g_0_xx_0_xzzzzzzz_0,  \
                             g_0_xx_0_xzzzzzzz_1,  \
                             g_0_xx_0_yyyyyyy_1,   \
                             g_0_xx_0_yyyyyyyy_0,  \
                             g_0_xx_0_yyyyyyyy_1,  \
                             g_0_xx_0_yyyyyyyz_0,  \
                             g_0_xx_0_yyyyyyyz_1,  \
                             g_0_xx_0_yyyyyyz_1,   \
                             g_0_xx_0_yyyyyyzz_0,  \
                             g_0_xx_0_yyyyyyzz_1,  \
                             g_0_xx_0_yyyyyzz_1,   \
                             g_0_xx_0_yyyyyzzz_0,  \
                             g_0_xx_0_yyyyyzzz_1,  \
                             g_0_xx_0_yyyyzzz_1,   \
                             g_0_xx_0_yyyyzzzz_0,  \
                             g_0_xx_0_yyyyzzzz_1,  \
                             g_0_xx_0_yyyzzzz_1,   \
                             g_0_xx_0_yyyzzzzz_0,  \
                             g_0_xx_0_yyyzzzzz_1,  \
                             g_0_xx_0_yyzzzzz_1,   \
                             g_0_xx_0_yyzzzzzz_0,  \
                             g_0_xx_0_yyzzzzzz_1,  \
                             g_0_xx_0_yzzzzzz_1,   \
                             g_0_xx_0_yzzzzzzz_0,  \
                             g_0_xx_0_yzzzzzzz_1,  \
                             g_0_xx_0_zzzzzzz_1,   \
                             g_0_xx_0_zzzzzzzz_0,  \
                             g_0_xx_0_zzzzzzzz_1,  \
                             g_0_xxy_0_xxxxxxxx_0, \
                             g_0_xxy_0_xxxxxxxy_0, \
                             g_0_xxy_0_xxxxxxxz_0, \
                             g_0_xxy_0_xxxxxxyy_0, \
                             g_0_xxy_0_xxxxxxyz_0, \
                             g_0_xxy_0_xxxxxxzz_0, \
                             g_0_xxy_0_xxxxxyyy_0, \
                             g_0_xxy_0_xxxxxyyz_0, \
                             g_0_xxy_0_xxxxxyzz_0, \
                             g_0_xxy_0_xxxxxzzz_0, \
                             g_0_xxy_0_xxxxyyyy_0, \
                             g_0_xxy_0_xxxxyyyz_0, \
                             g_0_xxy_0_xxxxyyzz_0, \
                             g_0_xxy_0_xxxxyzzz_0, \
                             g_0_xxy_0_xxxxzzzz_0, \
                             g_0_xxy_0_xxxyyyyy_0, \
                             g_0_xxy_0_xxxyyyyz_0, \
                             g_0_xxy_0_xxxyyyzz_0, \
                             g_0_xxy_0_xxxyyzzz_0, \
                             g_0_xxy_0_xxxyzzzz_0, \
                             g_0_xxy_0_xxxzzzzz_0, \
                             g_0_xxy_0_xxyyyyyy_0, \
                             g_0_xxy_0_xxyyyyyz_0, \
                             g_0_xxy_0_xxyyyyzz_0, \
                             g_0_xxy_0_xxyyyzzz_0, \
                             g_0_xxy_0_xxyyzzzz_0, \
                             g_0_xxy_0_xxyzzzzz_0, \
                             g_0_xxy_0_xxzzzzzz_0, \
                             g_0_xxy_0_xyyyyyyy_0, \
                             g_0_xxy_0_xyyyyyyz_0, \
                             g_0_xxy_0_xyyyyyzz_0, \
                             g_0_xxy_0_xyyyyzzz_0, \
                             g_0_xxy_0_xyyyzzzz_0, \
                             g_0_xxy_0_xyyzzzzz_0, \
                             g_0_xxy_0_xyzzzzzz_0, \
                             g_0_xxy_0_xzzzzzzz_0, \
                             g_0_xxy_0_yyyyyyyy_0, \
                             g_0_xxy_0_yyyyyyyz_0, \
                             g_0_xxy_0_yyyyyyzz_0, \
                             g_0_xxy_0_yyyyyzzz_0, \
                             g_0_xxy_0_yyyyzzzz_0, \
                             g_0_xxy_0_yyyzzzzz_0, \
                             g_0_xxy_0_yyzzzzzz_0, \
                             g_0_xxy_0_yzzzzzzz_0, \
                             g_0_xxy_0_zzzzzzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xxxxxxxx_0[i] = g_0_xx_0_xxxxxxxx_0[i] * pb_y + g_0_xx_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxxxy_0[i] = g_0_xx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxxy_0[i] * pb_y + g_0_xx_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxxxz_0[i] = g_0_xx_0_xxxxxxxz_0[i] * pb_y + g_0_xx_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxxyy_0[i] = 2.0 * g_0_xx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxyy_0[i] * pb_y + g_0_xx_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxxyz_0[i] = g_0_xx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxyz_0[i] * pb_y + g_0_xx_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxxzz_0[i] = g_0_xx_0_xxxxxxzz_0[i] * pb_y + g_0_xx_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxyyy_0[i] = 3.0 * g_0_xx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyyy_0[i] * pb_y + g_0_xx_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxyyz_0[i] = 2.0 * g_0_xx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyyz_0[i] * pb_y + g_0_xx_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxyzz_0[i] = g_0_xx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyzz_0[i] * pb_y + g_0_xx_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxzzz_0[i] = g_0_xx_0_xxxxxzzz_0[i] * pb_y + g_0_xx_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyyyy_0[i] = 4.0 * g_0_xx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyyy_0[i] * pb_y + g_0_xx_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyyyz_0[i] = 3.0 * g_0_xx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyyz_0[i] * pb_y + g_0_xx_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyyzz_0[i] = 2.0 * g_0_xx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyzz_0[i] * pb_y + g_0_xx_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyzzz_0[i] = g_0_xx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyzzz_0[i] * pb_y + g_0_xx_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxzzzz_0[i] = g_0_xx_0_xxxxzzzz_0[i] * pb_y + g_0_xx_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyyyy_0[i] = 5.0 * g_0_xx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyyy_0[i] * pb_y + g_0_xx_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyyyz_0[i] = 4.0 * g_0_xx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyyz_0[i] * pb_y + g_0_xx_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyyzz_0[i] = 3.0 * g_0_xx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyzz_0[i] * pb_y + g_0_xx_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyzzz_0[i] = 2.0 * g_0_xx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyzzz_0[i] * pb_y + g_0_xx_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyzzzz_0[i] = g_0_xx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzzzz_0[i] * pb_y + g_0_xx_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxzzzzz_0[i] = g_0_xx_0_xxxzzzzz_0[i] * pb_y + g_0_xx_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyyyy_0[i] = 6.0 * g_0_xx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyyy_0[i] * pb_y + g_0_xx_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyyyz_0[i] = 5.0 * g_0_xx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyyz_0[i] * pb_y + g_0_xx_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyyzz_0[i] = 4.0 * g_0_xx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyzz_0[i] * pb_y + g_0_xx_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyzzz_0[i] = 3.0 * g_0_xx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyzzz_0[i] * pb_y + g_0_xx_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyzzzz_0[i] = 2.0 * g_0_xx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzzzz_0[i] * pb_y + g_0_xx_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyzzzzz_0[i] = g_0_xx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzzzz_0[i] * pb_y + g_0_xx_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxzzzzzz_0[i] = g_0_xx_0_xxzzzzzz_0[i] * pb_y + g_0_xx_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyyyy_0[i] = 7.0 * g_0_xx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyyy_0[i] * pb_y + g_0_xx_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyyyz_0[i] = 6.0 * g_0_xx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyyz_0[i] * pb_y + g_0_xx_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyyzz_0[i] = 5.0 * g_0_xx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyzz_0[i] * pb_y + g_0_xx_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyzzz_0[i] = 4.0 * g_0_xx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyzzz_0[i] * pb_y + g_0_xx_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyzzzz_0[i] = 3.0 * g_0_xx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyzzzz_0[i] * pb_y + g_0_xx_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyzzzzz_0[i] = 2.0 * g_0_xx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzzzzz_0[i] * pb_y + g_0_xx_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyzzzzzz_0[i] = g_0_xx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzzzzz_0[i] * pb_y + g_0_xx_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xzzzzzzz_0[i] = g_0_xx_0_xzzzzzzz_0[i] * pb_y + g_0_xx_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyyyy_0[i] = 8.0 * g_0_xx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyyy_0[i] * pb_y + g_0_xx_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyyyz_0[i] = 7.0 * g_0_xx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyyz_0[i] * pb_y + g_0_xx_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyyzz_0[i] = 6.0 * g_0_xx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyzz_0[i] * pb_y + g_0_xx_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyzzz_0[i] = 5.0 * g_0_xx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyzzz_0[i] * pb_y + g_0_xx_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyzzzz_0[i] = 4.0 * g_0_xx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyzzzz_0[i] * pb_y + g_0_xx_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyzzzzz_0[i] = 3.0 * g_0_xx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzzzzz_0[i] * pb_y + g_0_xx_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyzzzzzz_0[i] = 2.0 * g_0_xx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzzzzz_0[i] * pb_y + g_0_xx_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yzzzzzzz_0[i] = g_0_xx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzzzzz_0[i] * pb_y + g_0_xx_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_zzzzzzzz_0[i] = g_0_xx_0_zzzzzzzz_0[i] * pb_y + g_0_xx_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 90-135 components of targeted buffer : SFSL

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

    auto g_0_xxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 126);

    auto g_0_xxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 127);

    auto g_0_xxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 128);

    auto g_0_xxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 129);

    auto g_0_xxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 130);

    auto g_0_xxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 131);

    auto g_0_xxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 132);

    auto g_0_xxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 133);

    auto g_0_xxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 134);

#pragma omp simd aligned(g_0_xx_0_xxxxxxx_1,       \
                             g_0_xx_0_xxxxxxxx_0,  \
                             g_0_xx_0_xxxxxxxx_1,  \
                             g_0_xx_0_xxxxxxxy_0,  \
                             g_0_xx_0_xxxxxxxy_1,  \
                             g_0_xx_0_xxxxxxxz_0,  \
                             g_0_xx_0_xxxxxxxz_1,  \
                             g_0_xx_0_xxxxxxy_1,   \
                             g_0_xx_0_xxxxxxyy_0,  \
                             g_0_xx_0_xxxxxxyy_1,  \
                             g_0_xx_0_xxxxxxyz_0,  \
                             g_0_xx_0_xxxxxxyz_1,  \
                             g_0_xx_0_xxxxxxz_1,   \
                             g_0_xx_0_xxxxxxzz_0,  \
                             g_0_xx_0_xxxxxxzz_1,  \
                             g_0_xx_0_xxxxxyy_1,   \
                             g_0_xx_0_xxxxxyyy_0,  \
                             g_0_xx_0_xxxxxyyy_1,  \
                             g_0_xx_0_xxxxxyyz_0,  \
                             g_0_xx_0_xxxxxyyz_1,  \
                             g_0_xx_0_xxxxxyz_1,   \
                             g_0_xx_0_xxxxxyzz_0,  \
                             g_0_xx_0_xxxxxyzz_1,  \
                             g_0_xx_0_xxxxxzz_1,   \
                             g_0_xx_0_xxxxxzzz_0,  \
                             g_0_xx_0_xxxxxzzz_1,  \
                             g_0_xx_0_xxxxyyy_1,   \
                             g_0_xx_0_xxxxyyyy_0,  \
                             g_0_xx_0_xxxxyyyy_1,  \
                             g_0_xx_0_xxxxyyyz_0,  \
                             g_0_xx_0_xxxxyyyz_1,  \
                             g_0_xx_0_xxxxyyz_1,   \
                             g_0_xx_0_xxxxyyzz_0,  \
                             g_0_xx_0_xxxxyyzz_1,  \
                             g_0_xx_0_xxxxyzz_1,   \
                             g_0_xx_0_xxxxyzzz_0,  \
                             g_0_xx_0_xxxxyzzz_1,  \
                             g_0_xx_0_xxxxzzz_1,   \
                             g_0_xx_0_xxxxzzzz_0,  \
                             g_0_xx_0_xxxxzzzz_1,  \
                             g_0_xx_0_xxxyyyy_1,   \
                             g_0_xx_0_xxxyyyyy_0,  \
                             g_0_xx_0_xxxyyyyy_1,  \
                             g_0_xx_0_xxxyyyyz_0,  \
                             g_0_xx_0_xxxyyyyz_1,  \
                             g_0_xx_0_xxxyyyz_1,   \
                             g_0_xx_0_xxxyyyzz_0,  \
                             g_0_xx_0_xxxyyyzz_1,  \
                             g_0_xx_0_xxxyyzz_1,   \
                             g_0_xx_0_xxxyyzzz_0,  \
                             g_0_xx_0_xxxyyzzz_1,  \
                             g_0_xx_0_xxxyzzz_1,   \
                             g_0_xx_0_xxxyzzzz_0,  \
                             g_0_xx_0_xxxyzzzz_1,  \
                             g_0_xx_0_xxxzzzz_1,   \
                             g_0_xx_0_xxxzzzzz_0,  \
                             g_0_xx_0_xxxzzzzz_1,  \
                             g_0_xx_0_xxyyyyy_1,   \
                             g_0_xx_0_xxyyyyyy_0,  \
                             g_0_xx_0_xxyyyyyy_1,  \
                             g_0_xx_0_xxyyyyyz_0,  \
                             g_0_xx_0_xxyyyyyz_1,  \
                             g_0_xx_0_xxyyyyz_1,   \
                             g_0_xx_0_xxyyyyzz_0,  \
                             g_0_xx_0_xxyyyyzz_1,  \
                             g_0_xx_0_xxyyyzz_1,   \
                             g_0_xx_0_xxyyyzzz_0,  \
                             g_0_xx_0_xxyyyzzz_1,  \
                             g_0_xx_0_xxyyzzz_1,   \
                             g_0_xx_0_xxyyzzzz_0,  \
                             g_0_xx_0_xxyyzzzz_1,  \
                             g_0_xx_0_xxyzzzz_1,   \
                             g_0_xx_0_xxyzzzzz_0,  \
                             g_0_xx_0_xxyzzzzz_1,  \
                             g_0_xx_0_xxzzzzz_1,   \
                             g_0_xx_0_xxzzzzzz_0,  \
                             g_0_xx_0_xxzzzzzz_1,  \
                             g_0_xx_0_xyyyyyy_1,   \
                             g_0_xx_0_xyyyyyyy_0,  \
                             g_0_xx_0_xyyyyyyy_1,  \
                             g_0_xx_0_xyyyyyyz_0,  \
                             g_0_xx_0_xyyyyyyz_1,  \
                             g_0_xx_0_xyyyyyz_1,   \
                             g_0_xx_0_xyyyyyzz_0,  \
                             g_0_xx_0_xyyyyyzz_1,  \
                             g_0_xx_0_xyyyyzz_1,   \
                             g_0_xx_0_xyyyyzzz_0,  \
                             g_0_xx_0_xyyyyzzz_1,  \
                             g_0_xx_0_xyyyzzz_1,   \
                             g_0_xx_0_xyyyzzzz_0,  \
                             g_0_xx_0_xyyyzzzz_1,  \
                             g_0_xx_0_xyyzzzz_1,   \
                             g_0_xx_0_xyyzzzzz_0,  \
                             g_0_xx_0_xyyzzzzz_1,  \
                             g_0_xx_0_xyzzzzz_1,   \
                             g_0_xx_0_xyzzzzzz_0,  \
                             g_0_xx_0_xyzzzzzz_1,  \
                             g_0_xx_0_xzzzzzz_1,   \
                             g_0_xx_0_xzzzzzzz_0,  \
                             g_0_xx_0_xzzzzzzz_1,  \
                             g_0_xx_0_yyyyyyy_1,   \
                             g_0_xx_0_yyyyyyyy_0,  \
                             g_0_xx_0_yyyyyyyy_1,  \
                             g_0_xx_0_yyyyyyyz_0,  \
                             g_0_xx_0_yyyyyyyz_1,  \
                             g_0_xx_0_yyyyyyz_1,   \
                             g_0_xx_0_yyyyyyzz_0,  \
                             g_0_xx_0_yyyyyyzz_1,  \
                             g_0_xx_0_yyyyyzz_1,   \
                             g_0_xx_0_yyyyyzzz_0,  \
                             g_0_xx_0_yyyyyzzz_1,  \
                             g_0_xx_0_yyyyzzz_1,   \
                             g_0_xx_0_yyyyzzzz_0,  \
                             g_0_xx_0_yyyyzzzz_1,  \
                             g_0_xx_0_yyyzzzz_1,   \
                             g_0_xx_0_yyyzzzzz_0,  \
                             g_0_xx_0_yyyzzzzz_1,  \
                             g_0_xx_0_yyzzzzz_1,   \
                             g_0_xx_0_yyzzzzzz_0,  \
                             g_0_xx_0_yyzzzzzz_1,  \
                             g_0_xx_0_yzzzzzz_1,   \
                             g_0_xx_0_yzzzzzzz_0,  \
                             g_0_xx_0_yzzzzzzz_1,  \
                             g_0_xx_0_zzzzzzz_1,   \
                             g_0_xx_0_zzzzzzzz_0,  \
                             g_0_xx_0_zzzzzzzz_1,  \
                             g_0_xxz_0_xxxxxxxx_0, \
                             g_0_xxz_0_xxxxxxxy_0, \
                             g_0_xxz_0_xxxxxxxz_0, \
                             g_0_xxz_0_xxxxxxyy_0, \
                             g_0_xxz_0_xxxxxxyz_0, \
                             g_0_xxz_0_xxxxxxzz_0, \
                             g_0_xxz_0_xxxxxyyy_0, \
                             g_0_xxz_0_xxxxxyyz_0, \
                             g_0_xxz_0_xxxxxyzz_0, \
                             g_0_xxz_0_xxxxxzzz_0, \
                             g_0_xxz_0_xxxxyyyy_0, \
                             g_0_xxz_0_xxxxyyyz_0, \
                             g_0_xxz_0_xxxxyyzz_0, \
                             g_0_xxz_0_xxxxyzzz_0, \
                             g_0_xxz_0_xxxxzzzz_0, \
                             g_0_xxz_0_xxxyyyyy_0, \
                             g_0_xxz_0_xxxyyyyz_0, \
                             g_0_xxz_0_xxxyyyzz_0, \
                             g_0_xxz_0_xxxyyzzz_0, \
                             g_0_xxz_0_xxxyzzzz_0, \
                             g_0_xxz_0_xxxzzzzz_0, \
                             g_0_xxz_0_xxyyyyyy_0, \
                             g_0_xxz_0_xxyyyyyz_0, \
                             g_0_xxz_0_xxyyyyzz_0, \
                             g_0_xxz_0_xxyyyzzz_0, \
                             g_0_xxz_0_xxyyzzzz_0, \
                             g_0_xxz_0_xxyzzzzz_0, \
                             g_0_xxz_0_xxzzzzzz_0, \
                             g_0_xxz_0_xyyyyyyy_0, \
                             g_0_xxz_0_xyyyyyyz_0, \
                             g_0_xxz_0_xyyyyyzz_0, \
                             g_0_xxz_0_xyyyyzzz_0, \
                             g_0_xxz_0_xyyyzzzz_0, \
                             g_0_xxz_0_xyyzzzzz_0, \
                             g_0_xxz_0_xyzzzzzz_0, \
                             g_0_xxz_0_xzzzzzzz_0, \
                             g_0_xxz_0_yyyyyyyy_0, \
                             g_0_xxz_0_yyyyyyyz_0, \
                             g_0_xxz_0_yyyyyyzz_0, \
                             g_0_xxz_0_yyyyyzzz_0, \
                             g_0_xxz_0_yyyyzzzz_0, \
                             g_0_xxz_0_yyyzzzzz_0, \
                             g_0_xxz_0_yyzzzzzz_0, \
                             g_0_xxz_0_yzzzzzzz_0, \
                             g_0_xxz_0_zzzzzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xxxxxxxx_0[i] = g_0_xx_0_xxxxxxxx_0[i] * pb_z + g_0_xx_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxxxy_0[i] = g_0_xx_0_xxxxxxxy_0[i] * pb_z + g_0_xx_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxxxz_0[i] = g_0_xx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxxz_0[i] * pb_z + g_0_xx_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxxyy_0[i] = g_0_xx_0_xxxxxxyy_0[i] * pb_z + g_0_xx_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxxyz_0[i] = g_0_xx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxyz_0[i] * pb_z + g_0_xx_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxxzz_0[i] = 2.0 * g_0_xx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxzz_0[i] * pb_z + g_0_xx_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxyyy_0[i] = g_0_xx_0_xxxxxyyy_0[i] * pb_z + g_0_xx_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxyyz_0[i] = g_0_xx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyyz_0[i] * pb_z + g_0_xx_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxyzz_0[i] = 2.0 * g_0_xx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyzz_0[i] * pb_z + g_0_xx_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxzzz_0[i] = 3.0 * g_0_xx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxzzz_0[i] * pb_z + g_0_xx_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyyyy_0[i] = g_0_xx_0_xxxxyyyy_0[i] * pb_z + g_0_xx_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyyyz_0[i] = g_0_xx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyyz_0[i] * pb_z + g_0_xx_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyyzz_0[i] = 2.0 * g_0_xx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyzz_0[i] * pb_z + g_0_xx_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyzzz_0[i] = 3.0 * g_0_xx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyzzz_0[i] * pb_z + g_0_xx_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxzzzz_0[i] = 4.0 * g_0_xx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxzzzz_0[i] * pb_z + g_0_xx_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyyyy_0[i] = g_0_xx_0_xxxyyyyy_0[i] * pb_z + g_0_xx_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyyyz_0[i] = g_0_xx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyyz_0[i] * pb_z + g_0_xx_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyyzz_0[i] = 2.0 * g_0_xx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyzz_0[i] * pb_z + g_0_xx_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyzzz_0[i] = 3.0 * g_0_xx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyzzz_0[i] * pb_z + g_0_xx_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyzzzz_0[i] = 4.0 * g_0_xx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzzzz_0[i] * pb_z + g_0_xx_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxzzzzz_0[i] = 5.0 * g_0_xx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxzzzzz_0[i] * pb_z + g_0_xx_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyyyy_0[i] = g_0_xx_0_xxyyyyyy_0[i] * pb_z + g_0_xx_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyyyz_0[i] = g_0_xx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyyz_0[i] * pb_z + g_0_xx_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyyzz_0[i] = 2.0 * g_0_xx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyzz_0[i] * pb_z + g_0_xx_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyzzz_0[i] = 3.0 * g_0_xx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyzzz_0[i] * pb_z + g_0_xx_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyzzzz_0[i] = 4.0 * g_0_xx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzzzz_0[i] * pb_z + g_0_xx_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyzzzzz_0[i] = 5.0 * g_0_xx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzzzz_0[i] * pb_z + g_0_xx_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxzzzzzz_0[i] = 6.0 * g_0_xx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzzzzzz_0[i] * pb_z + g_0_xx_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyyyy_0[i] = g_0_xx_0_xyyyyyyy_0[i] * pb_z + g_0_xx_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyyyz_0[i] = g_0_xx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyyz_0[i] * pb_z + g_0_xx_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyyzz_0[i] = 2.0 * g_0_xx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyzz_0[i] * pb_z + g_0_xx_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyzzz_0[i] = 3.0 * g_0_xx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyzzz_0[i] * pb_z + g_0_xx_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyzzzz_0[i] = 4.0 * g_0_xx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyzzzz_0[i] * pb_z + g_0_xx_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyzzzzz_0[i] = 5.0 * g_0_xx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzzzzz_0[i] * pb_z + g_0_xx_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyzzzzzz_0[i] = 6.0 * g_0_xx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzzzzz_0[i] * pb_z + g_0_xx_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xzzzzzzz_0[i] = 7.0 * g_0_xx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xzzzzzzz_0[i] * pb_z + g_0_xx_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyyyy_0[i] = g_0_xx_0_yyyyyyyy_0[i] * pb_z + g_0_xx_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyyyz_0[i] = g_0_xx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyyz_0[i] * pb_z + g_0_xx_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyyzz_0[i] = 2.0 * g_0_xx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyzz_0[i] * pb_z + g_0_xx_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyzzz_0[i] = 3.0 * g_0_xx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyzzz_0[i] * pb_z + g_0_xx_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyzzzz_0[i] = 4.0 * g_0_xx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyzzzz_0[i] * pb_z + g_0_xx_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyzzzzz_0[i] = 5.0 * g_0_xx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzzzzz_0[i] * pb_z + g_0_xx_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyzzzzzz_0[i] = 6.0 * g_0_xx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzzzzz_0[i] * pb_z + g_0_xx_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yzzzzzzz_0[i] = 7.0 * g_0_xx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzzzzz_0[i] * pb_z + g_0_xx_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_zzzzzzzz_0[i] = 8.0 * g_0_xx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_zzzzzzzz_0[i] * pb_z + g_0_xx_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 135-180 components of targeted buffer : SFSL

    auto g_0_xyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 135);

    auto g_0_xyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 136);

    auto g_0_xyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 137);

    auto g_0_xyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 138);

    auto g_0_xyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 139);

    auto g_0_xyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 140);

    auto g_0_xyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 141);

    auto g_0_xyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 142);

    auto g_0_xyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 143);

    auto g_0_xyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 144);

    auto g_0_xyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 145);

    auto g_0_xyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 146);

    auto g_0_xyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 147);

    auto g_0_xyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 148);

    auto g_0_xyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 149);

    auto g_0_xyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 150);

    auto g_0_xyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 151);

    auto g_0_xyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 152);

    auto g_0_xyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 153);

    auto g_0_xyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 154);

    auto g_0_xyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 155);

    auto g_0_xyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 156);

    auto g_0_xyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 157);

    auto g_0_xyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 158);

    auto g_0_xyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 159);

    auto g_0_xyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 160);

    auto g_0_xyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 161);

    auto g_0_xyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 162);

    auto g_0_xyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 163);

    auto g_0_xyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 164);

    auto g_0_xyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 165);

    auto g_0_xyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 166);

    auto g_0_xyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 167);

    auto g_0_xyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 168);

    auto g_0_xyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 169);

    auto g_0_xyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 170);

    auto g_0_xyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 171);

    auto g_0_xyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 172);

    auto g_0_xyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 173);

    auto g_0_xyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 174);

    auto g_0_xyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 175);

    auto g_0_xyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 176);

    auto g_0_xyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 177);

    auto g_0_xyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 178);

    auto g_0_xyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 179);

#pragma omp simd aligned(g_0_xyy_0_xxxxxxxx_0,     \
                             g_0_xyy_0_xxxxxxxy_0, \
                             g_0_xyy_0_xxxxxxxz_0, \
                             g_0_xyy_0_xxxxxxyy_0, \
                             g_0_xyy_0_xxxxxxyz_0, \
                             g_0_xyy_0_xxxxxxzz_0, \
                             g_0_xyy_0_xxxxxyyy_0, \
                             g_0_xyy_0_xxxxxyyz_0, \
                             g_0_xyy_0_xxxxxyzz_0, \
                             g_0_xyy_0_xxxxxzzz_0, \
                             g_0_xyy_0_xxxxyyyy_0, \
                             g_0_xyy_0_xxxxyyyz_0, \
                             g_0_xyy_0_xxxxyyzz_0, \
                             g_0_xyy_0_xxxxyzzz_0, \
                             g_0_xyy_0_xxxxzzzz_0, \
                             g_0_xyy_0_xxxyyyyy_0, \
                             g_0_xyy_0_xxxyyyyz_0, \
                             g_0_xyy_0_xxxyyyzz_0, \
                             g_0_xyy_0_xxxyyzzz_0, \
                             g_0_xyy_0_xxxyzzzz_0, \
                             g_0_xyy_0_xxxzzzzz_0, \
                             g_0_xyy_0_xxyyyyyy_0, \
                             g_0_xyy_0_xxyyyyyz_0, \
                             g_0_xyy_0_xxyyyyzz_0, \
                             g_0_xyy_0_xxyyyzzz_0, \
                             g_0_xyy_0_xxyyzzzz_0, \
                             g_0_xyy_0_xxyzzzzz_0, \
                             g_0_xyy_0_xxzzzzzz_0, \
                             g_0_xyy_0_xyyyyyyy_0, \
                             g_0_xyy_0_xyyyyyyz_0, \
                             g_0_xyy_0_xyyyyyzz_0, \
                             g_0_xyy_0_xyyyyzzz_0, \
                             g_0_xyy_0_xyyyzzzz_0, \
                             g_0_xyy_0_xyyzzzzz_0, \
                             g_0_xyy_0_xyzzzzzz_0, \
                             g_0_xyy_0_xzzzzzzz_0, \
                             g_0_xyy_0_yyyyyyyy_0, \
                             g_0_xyy_0_yyyyyyyz_0, \
                             g_0_xyy_0_yyyyyyzz_0, \
                             g_0_xyy_0_yyyyyzzz_0, \
                             g_0_xyy_0_yyyyzzzz_0, \
                             g_0_xyy_0_yyyzzzzz_0, \
                             g_0_xyy_0_yyzzzzzz_0, \
                             g_0_xyy_0_yzzzzzzz_0, \
                             g_0_xyy_0_zzzzzzzz_0, \
                             g_0_yy_0_xxxxxxx_1,   \
                             g_0_yy_0_xxxxxxxx_0,  \
                             g_0_yy_0_xxxxxxxx_1,  \
                             g_0_yy_0_xxxxxxxy_0,  \
                             g_0_yy_0_xxxxxxxy_1,  \
                             g_0_yy_0_xxxxxxxz_0,  \
                             g_0_yy_0_xxxxxxxz_1,  \
                             g_0_yy_0_xxxxxxy_1,   \
                             g_0_yy_0_xxxxxxyy_0,  \
                             g_0_yy_0_xxxxxxyy_1,  \
                             g_0_yy_0_xxxxxxyz_0,  \
                             g_0_yy_0_xxxxxxyz_1,  \
                             g_0_yy_0_xxxxxxz_1,   \
                             g_0_yy_0_xxxxxxzz_0,  \
                             g_0_yy_0_xxxxxxzz_1,  \
                             g_0_yy_0_xxxxxyy_1,   \
                             g_0_yy_0_xxxxxyyy_0,  \
                             g_0_yy_0_xxxxxyyy_1,  \
                             g_0_yy_0_xxxxxyyz_0,  \
                             g_0_yy_0_xxxxxyyz_1,  \
                             g_0_yy_0_xxxxxyz_1,   \
                             g_0_yy_0_xxxxxyzz_0,  \
                             g_0_yy_0_xxxxxyzz_1,  \
                             g_0_yy_0_xxxxxzz_1,   \
                             g_0_yy_0_xxxxxzzz_0,  \
                             g_0_yy_0_xxxxxzzz_1,  \
                             g_0_yy_0_xxxxyyy_1,   \
                             g_0_yy_0_xxxxyyyy_0,  \
                             g_0_yy_0_xxxxyyyy_1,  \
                             g_0_yy_0_xxxxyyyz_0,  \
                             g_0_yy_0_xxxxyyyz_1,  \
                             g_0_yy_0_xxxxyyz_1,   \
                             g_0_yy_0_xxxxyyzz_0,  \
                             g_0_yy_0_xxxxyyzz_1,  \
                             g_0_yy_0_xxxxyzz_1,   \
                             g_0_yy_0_xxxxyzzz_0,  \
                             g_0_yy_0_xxxxyzzz_1,  \
                             g_0_yy_0_xxxxzzz_1,   \
                             g_0_yy_0_xxxxzzzz_0,  \
                             g_0_yy_0_xxxxzzzz_1,  \
                             g_0_yy_0_xxxyyyy_1,   \
                             g_0_yy_0_xxxyyyyy_0,  \
                             g_0_yy_0_xxxyyyyy_1,  \
                             g_0_yy_0_xxxyyyyz_0,  \
                             g_0_yy_0_xxxyyyyz_1,  \
                             g_0_yy_0_xxxyyyz_1,   \
                             g_0_yy_0_xxxyyyzz_0,  \
                             g_0_yy_0_xxxyyyzz_1,  \
                             g_0_yy_0_xxxyyzz_1,   \
                             g_0_yy_0_xxxyyzzz_0,  \
                             g_0_yy_0_xxxyyzzz_1,  \
                             g_0_yy_0_xxxyzzz_1,   \
                             g_0_yy_0_xxxyzzzz_0,  \
                             g_0_yy_0_xxxyzzzz_1,  \
                             g_0_yy_0_xxxzzzz_1,   \
                             g_0_yy_0_xxxzzzzz_0,  \
                             g_0_yy_0_xxxzzzzz_1,  \
                             g_0_yy_0_xxyyyyy_1,   \
                             g_0_yy_0_xxyyyyyy_0,  \
                             g_0_yy_0_xxyyyyyy_1,  \
                             g_0_yy_0_xxyyyyyz_0,  \
                             g_0_yy_0_xxyyyyyz_1,  \
                             g_0_yy_0_xxyyyyz_1,   \
                             g_0_yy_0_xxyyyyzz_0,  \
                             g_0_yy_0_xxyyyyzz_1,  \
                             g_0_yy_0_xxyyyzz_1,   \
                             g_0_yy_0_xxyyyzzz_0,  \
                             g_0_yy_0_xxyyyzzz_1,  \
                             g_0_yy_0_xxyyzzz_1,   \
                             g_0_yy_0_xxyyzzzz_0,  \
                             g_0_yy_0_xxyyzzzz_1,  \
                             g_0_yy_0_xxyzzzz_1,   \
                             g_0_yy_0_xxyzzzzz_0,  \
                             g_0_yy_0_xxyzzzzz_1,  \
                             g_0_yy_0_xxzzzzz_1,   \
                             g_0_yy_0_xxzzzzzz_0,  \
                             g_0_yy_0_xxzzzzzz_1,  \
                             g_0_yy_0_xyyyyyy_1,   \
                             g_0_yy_0_xyyyyyyy_0,  \
                             g_0_yy_0_xyyyyyyy_1,  \
                             g_0_yy_0_xyyyyyyz_0,  \
                             g_0_yy_0_xyyyyyyz_1,  \
                             g_0_yy_0_xyyyyyz_1,   \
                             g_0_yy_0_xyyyyyzz_0,  \
                             g_0_yy_0_xyyyyyzz_1,  \
                             g_0_yy_0_xyyyyzz_1,   \
                             g_0_yy_0_xyyyyzzz_0,  \
                             g_0_yy_0_xyyyyzzz_1,  \
                             g_0_yy_0_xyyyzzz_1,   \
                             g_0_yy_0_xyyyzzzz_0,  \
                             g_0_yy_0_xyyyzzzz_1,  \
                             g_0_yy_0_xyyzzzz_1,   \
                             g_0_yy_0_xyyzzzzz_0,  \
                             g_0_yy_0_xyyzzzzz_1,  \
                             g_0_yy_0_xyzzzzz_1,   \
                             g_0_yy_0_xyzzzzzz_0,  \
                             g_0_yy_0_xyzzzzzz_1,  \
                             g_0_yy_0_xzzzzzz_1,   \
                             g_0_yy_0_xzzzzzzz_0,  \
                             g_0_yy_0_xzzzzzzz_1,  \
                             g_0_yy_0_yyyyyyy_1,   \
                             g_0_yy_0_yyyyyyyy_0,  \
                             g_0_yy_0_yyyyyyyy_1,  \
                             g_0_yy_0_yyyyyyyz_0,  \
                             g_0_yy_0_yyyyyyyz_1,  \
                             g_0_yy_0_yyyyyyz_1,   \
                             g_0_yy_0_yyyyyyzz_0,  \
                             g_0_yy_0_yyyyyyzz_1,  \
                             g_0_yy_0_yyyyyzz_1,   \
                             g_0_yy_0_yyyyyzzz_0,  \
                             g_0_yy_0_yyyyyzzz_1,  \
                             g_0_yy_0_yyyyzzz_1,   \
                             g_0_yy_0_yyyyzzzz_0,  \
                             g_0_yy_0_yyyyzzzz_1,  \
                             g_0_yy_0_yyyzzzz_1,   \
                             g_0_yy_0_yyyzzzzz_0,  \
                             g_0_yy_0_yyyzzzzz_1,  \
                             g_0_yy_0_yyzzzzz_1,   \
                             g_0_yy_0_yyzzzzzz_0,  \
                             g_0_yy_0_yyzzzzzz_1,  \
                             g_0_yy_0_yzzzzzz_1,   \
                             g_0_yy_0_yzzzzzzz_0,  \
                             g_0_yy_0_yzzzzzzz_1,  \
                             g_0_yy_0_zzzzzzz_1,   \
                             g_0_yy_0_zzzzzzzz_0,  \
                             g_0_yy_0_zzzzzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xxxxxxxx_0[i] = 8.0 * g_0_yy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxxx_0[i] * pb_x + g_0_yy_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxxxy_0[i] = 7.0 * g_0_yy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxxy_0[i] * pb_x + g_0_yy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxxxz_0[i] = 7.0 * g_0_yy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxxz_0[i] * pb_x + g_0_yy_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxxyy_0[i] = 6.0 * g_0_yy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxyy_0[i] * pb_x + g_0_yy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxxyz_0[i] = 6.0 * g_0_yy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxyz_0[i] * pb_x + g_0_yy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxxzz_0[i] = 6.0 * g_0_yy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxzz_0[i] * pb_x + g_0_yy_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxyyy_0[i] = 5.0 * g_0_yy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyyy_0[i] * pb_x + g_0_yy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxyyz_0[i] = 5.0 * g_0_yy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyyz_0[i] * pb_x + g_0_yy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxyzz_0[i] = 5.0 * g_0_yy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyzz_0[i] * pb_x + g_0_yy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxzzz_0[i] = 5.0 * g_0_yy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxzzz_0[i] * pb_x + g_0_yy_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyyyy_0[i] = 4.0 * g_0_yy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyyy_0[i] * pb_x + g_0_yy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyyyz_0[i] = 4.0 * g_0_yy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyyz_0[i] * pb_x + g_0_yy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyyzz_0[i] = 4.0 * g_0_yy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyzz_0[i] * pb_x + g_0_yy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyzzz_0[i] = 4.0 * g_0_yy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyzzz_0[i] * pb_x + g_0_yy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxzzzz_0[i] = 4.0 * g_0_yy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxzzzz_0[i] * pb_x + g_0_yy_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyyyy_0[i] = 3.0 * g_0_yy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyyy_0[i] * pb_x + g_0_yy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyyyz_0[i] = 3.0 * g_0_yy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyyz_0[i] * pb_x + g_0_yy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyyzz_0[i] = 3.0 * g_0_yy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyzz_0[i] * pb_x + g_0_yy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyzzz_0[i] = 3.0 * g_0_yy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyzzz_0[i] * pb_x + g_0_yy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyzzzz_0[i] = 3.0 * g_0_yy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyzzzz_0[i] * pb_x + g_0_yy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxzzzzz_0[i] = 3.0 * g_0_yy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzzzzz_0[i] * pb_x + g_0_yy_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyyyy_0[i] = 2.0 * g_0_yy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyyy_0[i] * pb_x + g_0_yy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyyyz_0[i] = 2.0 * g_0_yy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyyz_0[i] * pb_x + g_0_yy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyyzz_0[i] = 2.0 * g_0_yy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyzz_0[i] * pb_x + g_0_yy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyzzz_0[i] = 2.0 * g_0_yy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyzzz_0[i] * pb_x + g_0_yy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyzzzz_0[i] = 2.0 * g_0_yy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzzzz_0[i] * pb_x + g_0_yy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyzzzzz_0[i] = 2.0 * g_0_yy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzzzzz_0[i] * pb_x + g_0_yy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxzzzzzz_0[i] = 2.0 * g_0_yy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzzzzz_0[i] * pb_x + g_0_yy_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyyyy_0[i] = g_0_yy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyyy_0[i] * pb_x + g_0_yy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyyyz_0[i] = g_0_yy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyyz_0[i] * pb_x + g_0_yy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyyzz_0[i] = g_0_yy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyzz_0[i] * pb_x + g_0_yy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyzzz_0[i] = g_0_yy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyzzz_0[i] * pb_x + g_0_yy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyzzzz_0[i] = g_0_yy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzzzz_0[i] * pb_x + g_0_yy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyzzzzz_0[i] = g_0_yy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzzzz_0[i] * pb_x + g_0_yy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyzzzzzz_0[i] = g_0_yy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzzzzz_0[i] * pb_x + g_0_yy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xzzzzzzz_0[i] = g_0_yy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzzzzz_0[i] * pb_x + g_0_yy_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyyyy_0[i] = g_0_yy_0_yyyyyyyy_0[i] * pb_x + g_0_yy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyyyz_0[i] = g_0_yy_0_yyyyyyyz_0[i] * pb_x + g_0_yy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyyzz_0[i] = g_0_yy_0_yyyyyyzz_0[i] * pb_x + g_0_yy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyzzz_0[i] = g_0_yy_0_yyyyyzzz_0[i] * pb_x + g_0_yy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyzzzz_0[i] = g_0_yy_0_yyyyzzzz_0[i] * pb_x + g_0_yy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyzzzzz_0[i] = g_0_yy_0_yyyzzzzz_0[i] * pb_x + g_0_yy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyzzzzzz_0[i] = g_0_yy_0_yyzzzzzz_0[i] * pb_x + g_0_yy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yzzzzzzz_0[i] = g_0_yy_0_yzzzzzzz_0[i] * pb_x + g_0_yy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_zzzzzzzz_0[i] = g_0_yy_0_zzzzzzzz_0[i] * pb_x + g_0_yy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 180-225 components of targeted buffer : SFSL

    auto g_0_xyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 180);

    auto g_0_xyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 181);

    auto g_0_xyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 182);

    auto g_0_xyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 183);

    auto g_0_xyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 184);

    auto g_0_xyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 185);

    auto g_0_xyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 186);

    auto g_0_xyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 187);

    auto g_0_xyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 188);

    auto g_0_xyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 189);

    auto g_0_xyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 190);

    auto g_0_xyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 191);

    auto g_0_xyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 192);

    auto g_0_xyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 193);

    auto g_0_xyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 194);

    auto g_0_xyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 195);

    auto g_0_xyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 196);

    auto g_0_xyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 197);

    auto g_0_xyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 198);

    auto g_0_xyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 199);

    auto g_0_xyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 200);

    auto g_0_xyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 201);

    auto g_0_xyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 202);

    auto g_0_xyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 203);

    auto g_0_xyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 204);

    auto g_0_xyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 205);

    auto g_0_xyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 206);

    auto g_0_xyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 207);

    auto g_0_xyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 208);

    auto g_0_xyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 209);

    auto g_0_xyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 210);

    auto g_0_xyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 211);

    auto g_0_xyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 212);

    auto g_0_xyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 213);

    auto g_0_xyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 214);

    auto g_0_xyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 215);

    auto g_0_xyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 216);

    auto g_0_xyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 217);

    auto g_0_xyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 218);

    auto g_0_xyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 219);

    auto g_0_xyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 220);

    auto g_0_xyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 221);

    auto g_0_xyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 222);

    auto g_0_xyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 223);

    auto g_0_xyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 224);

#pragma omp simd aligned(g_0_xy_0_xxxxxxxy_0,      \
                             g_0_xy_0_xxxxxxxy_1,  \
                             g_0_xy_0_xxxxxxyy_0,  \
                             g_0_xy_0_xxxxxxyy_1,  \
                             g_0_xy_0_xxxxxyyy_0,  \
                             g_0_xy_0_xxxxxyyy_1,  \
                             g_0_xy_0_xxxxyyyy_0,  \
                             g_0_xy_0_xxxxyyyy_1,  \
                             g_0_xy_0_xxxyyyyy_0,  \
                             g_0_xy_0_xxxyyyyy_1,  \
                             g_0_xy_0_xxyyyyyy_0,  \
                             g_0_xy_0_xxyyyyyy_1,  \
                             g_0_xy_0_xyyyyyyy_0,  \
                             g_0_xy_0_xyyyyyyy_1,  \
                             g_0_xyz_0_xxxxxxxx_0, \
                             g_0_xyz_0_xxxxxxxy_0, \
                             g_0_xyz_0_xxxxxxxz_0, \
                             g_0_xyz_0_xxxxxxyy_0, \
                             g_0_xyz_0_xxxxxxyz_0, \
                             g_0_xyz_0_xxxxxxzz_0, \
                             g_0_xyz_0_xxxxxyyy_0, \
                             g_0_xyz_0_xxxxxyyz_0, \
                             g_0_xyz_0_xxxxxyzz_0, \
                             g_0_xyz_0_xxxxxzzz_0, \
                             g_0_xyz_0_xxxxyyyy_0, \
                             g_0_xyz_0_xxxxyyyz_0, \
                             g_0_xyz_0_xxxxyyzz_0, \
                             g_0_xyz_0_xxxxyzzz_0, \
                             g_0_xyz_0_xxxxzzzz_0, \
                             g_0_xyz_0_xxxyyyyy_0, \
                             g_0_xyz_0_xxxyyyyz_0, \
                             g_0_xyz_0_xxxyyyzz_0, \
                             g_0_xyz_0_xxxyyzzz_0, \
                             g_0_xyz_0_xxxyzzzz_0, \
                             g_0_xyz_0_xxxzzzzz_0, \
                             g_0_xyz_0_xxyyyyyy_0, \
                             g_0_xyz_0_xxyyyyyz_0, \
                             g_0_xyz_0_xxyyyyzz_0, \
                             g_0_xyz_0_xxyyyzzz_0, \
                             g_0_xyz_0_xxyyzzzz_0, \
                             g_0_xyz_0_xxyzzzzz_0, \
                             g_0_xyz_0_xxzzzzzz_0, \
                             g_0_xyz_0_xyyyyyyy_0, \
                             g_0_xyz_0_xyyyyyyz_0, \
                             g_0_xyz_0_xyyyyyzz_0, \
                             g_0_xyz_0_xyyyyzzz_0, \
                             g_0_xyz_0_xyyyzzzz_0, \
                             g_0_xyz_0_xyyzzzzz_0, \
                             g_0_xyz_0_xyzzzzzz_0, \
                             g_0_xyz_0_xzzzzzzz_0, \
                             g_0_xyz_0_yyyyyyyy_0, \
                             g_0_xyz_0_yyyyyyyz_0, \
                             g_0_xyz_0_yyyyyyzz_0, \
                             g_0_xyz_0_yyyyyzzz_0, \
                             g_0_xyz_0_yyyyzzzz_0, \
                             g_0_xyz_0_yyyzzzzz_0, \
                             g_0_xyz_0_yyzzzzzz_0, \
                             g_0_xyz_0_yzzzzzzz_0, \
                             g_0_xyz_0_zzzzzzzz_0, \
                             g_0_xz_0_xxxxxxxx_0,  \
                             g_0_xz_0_xxxxxxxx_1,  \
                             g_0_xz_0_xxxxxxxz_0,  \
                             g_0_xz_0_xxxxxxxz_1,  \
                             g_0_xz_0_xxxxxxzz_0,  \
                             g_0_xz_0_xxxxxxzz_1,  \
                             g_0_xz_0_xxxxxzzz_0,  \
                             g_0_xz_0_xxxxxzzz_1,  \
                             g_0_xz_0_xxxxzzzz_0,  \
                             g_0_xz_0_xxxxzzzz_1,  \
                             g_0_xz_0_xxxzzzzz_0,  \
                             g_0_xz_0_xxxzzzzz_1,  \
                             g_0_xz_0_xxzzzzzz_0,  \
                             g_0_xz_0_xxzzzzzz_1,  \
                             g_0_xz_0_xzzzzzzz_0,  \
                             g_0_xz_0_xzzzzzzz_1,  \
                             g_0_yz_0_xxxxxxyz_0,  \
                             g_0_yz_0_xxxxxxyz_1,  \
                             g_0_yz_0_xxxxxyyz_0,  \
                             g_0_yz_0_xxxxxyyz_1,  \
                             g_0_yz_0_xxxxxyz_1,   \
                             g_0_yz_0_xxxxxyzz_0,  \
                             g_0_yz_0_xxxxxyzz_1,  \
                             g_0_yz_0_xxxxyyyz_0,  \
                             g_0_yz_0_xxxxyyyz_1,  \
                             g_0_yz_0_xxxxyyz_1,   \
                             g_0_yz_0_xxxxyyzz_0,  \
                             g_0_yz_0_xxxxyyzz_1,  \
                             g_0_yz_0_xxxxyzz_1,   \
                             g_0_yz_0_xxxxyzzz_0,  \
                             g_0_yz_0_xxxxyzzz_1,  \
                             g_0_yz_0_xxxyyyyz_0,  \
                             g_0_yz_0_xxxyyyyz_1,  \
                             g_0_yz_0_xxxyyyz_1,   \
                             g_0_yz_0_xxxyyyzz_0,  \
                             g_0_yz_0_xxxyyyzz_1,  \
                             g_0_yz_0_xxxyyzz_1,   \
                             g_0_yz_0_xxxyyzzz_0,  \
                             g_0_yz_0_xxxyyzzz_1,  \
                             g_0_yz_0_xxxyzzz_1,   \
                             g_0_yz_0_xxxyzzzz_0,  \
                             g_0_yz_0_xxxyzzzz_1,  \
                             g_0_yz_0_xxyyyyyz_0,  \
                             g_0_yz_0_xxyyyyyz_1,  \
                             g_0_yz_0_xxyyyyz_1,   \
                             g_0_yz_0_xxyyyyzz_0,  \
                             g_0_yz_0_xxyyyyzz_1,  \
                             g_0_yz_0_xxyyyzz_1,   \
                             g_0_yz_0_xxyyyzzz_0,  \
                             g_0_yz_0_xxyyyzzz_1,  \
                             g_0_yz_0_xxyyzzz_1,   \
                             g_0_yz_0_xxyyzzzz_0,  \
                             g_0_yz_0_xxyyzzzz_1,  \
                             g_0_yz_0_xxyzzzz_1,   \
                             g_0_yz_0_xxyzzzzz_0,  \
                             g_0_yz_0_xxyzzzzz_1,  \
                             g_0_yz_0_xyyyyyyz_0,  \
                             g_0_yz_0_xyyyyyyz_1,  \
                             g_0_yz_0_xyyyyyz_1,   \
                             g_0_yz_0_xyyyyyzz_0,  \
                             g_0_yz_0_xyyyyyzz_1,  \
                             g_0_yz_0_xyyyyzz_1,   \
                             g_0_yz_0_xyyyyzzz_0,  \
                             g_0_yz_0_xyyyyzzz_1,  \
                             g_0_yz_0_xyyyzzz_1,   \
                             g_0_yz_0_xyyyzzzz_0,  \
                             g_0_yz_0_xyyyzzzz_1,  \
                             g_0_yz_0_xyyzzzz_1,   \
                             g_0_yz_0_xyyzzzzz_0,  \
                             g_0_yz_0_xyyzzzzz_1,  \
                             g_0_yz_0_xyzzzzz_1,   \
                             g_0_yz_0_xyzzzzzz_0,  \
                             g_0_yz_0_xyzzzzzz_1,  \
                             g_0_yz_0_yyyyyyyy_0,  \
                             g_0_yz_0_yyyyyyyy_1,  \
                             g_0_yz_0_yyyyyyyz_0,  \
                             g_0_yz_0_yyyyyyyz_1,  \
                             g_0_yz_0_yyyyyyz_1,   \
                             g_0_yz_0_yyyyyyzz_0,  \
                             g_0_yz_0_yyyyyyzz_1,  \
                             g_0_yz_0_yyyyyzz_1,   \
                             g_0_yz_0_yyyyyzzz_0,  \
                             g_0_yz_0_yyyyyzzz_1,  \
                             g_0_yz_0_yyyyzzz_1,   \
                             g_0_yz_0_yyyyzzzz_0,  \
                             g_0_yz_0_yyyyzzzz_1,  \
                             g_0_yz_0_yyyzzzz_1,   \
                             g_0_yz_0_yyyzzzzz_0,  \
                             g_0_yz_0_yyyzzzzz_1,  \
                             g_0_yz_0_yyzzzzz_1,   \
                             g_0_yz_0_yyzzzzzz_0,  \
                             g_0_yz_0_yyzzzzzz_1,  \
                             g_0_yz_0_yzzzzzz_1,   \
                             g_0_yz_0_yzzzzzzz_0,  \
                             g_0_yz_0_yzzzzzzz_1,  \
                             g_0_yz_0_zzzzzzzz_0,  \
                             g_0_yz_0_zzzzzzzz_1,  \
                             wp_x,                 \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyz_0_xxxxxxxx_0[i] = g_0_xz_0_xxxxxxxx_0[i] * pb_y + g_0_xz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xyz_0_xxxxxxxy_0[i] = g_0_xy_0_xxxxxxxy_0[i] * pb_z + g_0_xy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxxxxz_0[i] = g_0_xz_0_xxxxxxxz_0[i] * pb_y + g_0_xz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xyz_0_xxxxxxyy_0[i] = g_0_xy_0_xxxxxxyy_0[i] * pb_z + g_0_xy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxxxyz_0[i] = 6.0 * g_0_yz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxxxyz_0[i] * pb_x + g_0_yz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxxxzz_0[i] = g_0_xz_0_xxxxxxzz_0[i] * pb_y + g_0_xz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xyz_0_xxxxxyyy_0[i] = g_0_xy_0_xxxxxyyy_0[i] * pb_z + g_0_xy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxxyyz_0[i] = 5.0 * g_0_yz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxxyyz_0[i] * pb_x + g_0_yz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxxyzz_0[i] = 5.0 * g_0_yz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxxyzz_0[i] * pb_x + g_0_yz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxxzzz_0[i] = g_0_xz_0_xxxxxzzz_0[i] * pb_y + g_0_xz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xyz_0_xxxxyyyy_0[i] = g_0_xy_0_xxxxyyyy_0[i] * pb_z + g_0_xy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxyyyz_0[i] = 4.0 * g_0_yz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxyyyz_0[i] * pb_x + g_0_yz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxyyzz_0[i] = 4.0 * g_0_yz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxyyzz_0[i] * pb_x + g_0_yz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxyzzz_0[i] = 4.0 * g_0_yz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxyzzz_0[i] * pb_x + g_0_yz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxzzzz_0[i] = g_0_xz_0_xxxxzzzz_0[i] * pb_y + g_0_xz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xyz_0_xxxyyyyy_0[i] = g_0_xy_0_xxxyyyyy_0[i] * pb_z + g_0_xy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxyyyyz_0[i] = 3.0 * g_0_yz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyyyyz_0[i] * pb_x + g_0_yz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxyyyzz_0[i] = 3.0 * g_0_yz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyyyzz_0[i] * pb_x + g_0_yz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxyyzzz_0[i] = 3.0 * g_0_yz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyyzzz_0[i] * pb_x + g_0_yz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxyzzzz_0[i] = 3.0 * g_0_yz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyzzzz_0[i] * pb_x + g_0_yz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxzzzzz_0[i] = g_0_xz_0_xxxzzzzz_0[i] * pb_y + g_0_xz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xyz_0_xxyyyyyy_0[i] = g_0_xy_0_xxyyyyyy_0[i] * pb_z + g_0_xy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxyyyyyz_0[i] = 2.0 * g_0_yz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyyyyz_0[i] * pb_x + g_0_yz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxyyyyzz_0[i] = 2.0 * g_0_yz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyyyzz_0[i] * pb_x + g_0_yz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxyyyzzz_0[i] = 2.0 * g_0_yz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyyzzz_0[i] * pb_x + g_0_yz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxyyzzzz_0[i] = 2.0 * g_0_yz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyzzzz_0[i] * pb_x + g_0_yz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxyzzzzz_0[i] = 2.0 * g_0_yz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyzzzzz_0[i] * pb_x + g_0_yz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxzzzzzz_0[i] = g_0_xz_0_xxzzzzzz_0[i] * pb_y + g_0_xz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xyz_0_xyyyyyyy_0[i] = g_0_xy_0_xyyyyyyy_0[i] * pb_z + g_0_xy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xyyyyyyz_0[i] = g_0_yz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyyyyz_0[i] * pb_x + g_0_yz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xyyyyyzz_0[i] = g_0_yz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyyyzz_0[i] * pb_x + g_0_yz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xyyyyzzz_0[i] = g_0_yz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyyzzz_0[i] * pb_x + g_0_yz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xyyyzzzz_0[i] = g_0_yz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyzzzz_0[i] * pb_x + g_0_yz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xyyzzzzz_0[i] = g_0_yz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyzzzzz_0[i] * pb_x + g_0_yz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xyzzzzzz_0[i] = g_0_yz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyzzzzzz_0[i] * pb_x + g_0_yz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xzzzzzzz_0[i] = g_0_xz_0_xzzzzzzz_0[i] * pb_y + g_0_xz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xyz_0_yyyyyyyy_0[i] = g_0_yz_0_yyyyyyyy_0[i] * pb_x + g_0_yz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyz_0_yyyyyyyz_0[i] = g_0_yz_0_yyyyyyyz_0[i] * pb_x + g_0_yz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_yyyyyyzz_0[i] = g_0_yz_0_yyyyyyzz_0[i] * pb_x + g_0_yz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_yyyyyzzz_0[i] = g_0_yz_0_yyyyyzzz_0[i] * pb_x + g_0_yz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_yyyyzzzz_0[i] = g_0_yz_0_yyyyzzzz_0[i] * pb_x + g_0_yz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_yyyzzzzz_0[i] = g_0_yz_0_yyyzzzzz_0[i] * pb_x + g_0_yz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_yyzzzzzz_0[i] = g_0_yz_0_yyzzzzzz_0[i] * pb_x + g_0_yz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_yzzzzzzz_0[i] = g_0_yz_0_yzzzzzzz_0[i] * pb_x + g_0_yz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_zzzzzzzz_0[i] = g_0_yz_0_zzzzzzzz_0[i] * pb_x + g_0_yz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 225-270 components of targeted buffer : SFSL

    auto g_0_xzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 225);

    auto g_0_xzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sfsl + 226);

    auto g_0_xzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sfsl + 227);

    auto g_0_xzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sfsl + 228);

    auto g_0_xzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sfsl + 229);

    auto g_0_xzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sfsl + 230);

    auto g_0_xzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sfsl + 231);

    auto g_0_xzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sfsl + 232);

    auto g_0_xzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sfsl + 233);

    auto g_0_xzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sfsl + 234);

    auto g_0_xzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 235);

    auto g_0_xzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 236);

    auto g_0_xzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 237);

    auto g_0_xzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 238);

    auto g_0_xzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 239);

    auto g_0_xzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 240);

    auto g_0_xzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 241);

    auto g_0_xzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 242);

    auto g_0_xzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 243);

    auto g_0_xzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 244);

    auto g_0_xzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 245);

    auto g_0_xzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 246);

    auto g_0_xzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sfsl + 247);

    auto g_0_xzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sfsl + 248);

    auto g_0_xzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sfsl + 249);

    auto g_0_xzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 250);

    auto g_0_xzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 251);

    auto g_0_xzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sfsl + 252);

    auto g_0_xzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sfsl + 253);

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

#pragma omp simd aligned(g_0_xzz_0_xxxxxxxx_0,     \
                             g_0_xzz_0_xxxxxxxy_0, \
                             g_0_xzz_0_xxxxxxxz_0, \
                             g_0_xzz_0_xxxxxxyy_0, \
                             g_0_xzz_0_xxxxxxyz_0, \
                             g_0_xzz_0_xxxxxxzz_0, \
                             g_0_xzz_0_xxxxxyyy_0, \
                             g_0_xzz_0_xxxxxyyz_0, \
                             g_0_xzz_0_xxxxxyzz_0, \
                             g_0_xzz_0_xxxxxzzz_0, \
                             g_0_xzz_0_xxxxyyyy_0, \
                             g_0_xzz_0_xxxxyyyz_0, \
                             g_0_xzz_0_xxxxyyzz_0, \
                             g_0_xzz_0_xxxxyzzz_0, \
                             g_0_xzz_0_xxxxzzzz_0, \
                             g_0_xzz_0_xxxyyyyy_0, \
                             g_0_xzz_0_xxxyyyyz_0, \
                             g_0_xzz_0_xxxyyyzz_0, \
                             g_0_xzz_0_xxxyyzzz_0, \
                             g_0_xzz_0_xxxyzzzz_0, \
                             g_0_xzz_0_xxxzzzzz_0, \
                             g_0_xzz_0_xxyyyyyy_0, \
                             g_0_xzz_0_xxyyyyyz_0, \
                             g_0_xzz_0_xxyyyyzz_0, \
                             g_0_xzz_0_xxyyyzzz_0, \
                             g_0_xzz_0_xxyyzzzz_0, \
                             g_0_xzz_0_xxyzzzzz_0, \
                             g_0_xzz_0_xxzzzzzz_0, \
                             g_0_xzz_0_xyyyyyyy_0, \
                             g_0_xzz_0_xyyyyyyz_0, \
                             g_0_xzz_0_xyyyyyzz_0, \
                             g_0_xzz_0_xyyyyzzz_0, \
                             g_0_xzz_0_xyyyzzzz_0, \
                             g_0_xzz_0_xyyzzzzz_0, \
                             g_0_xzz_0_xyzzzzzz_0, \
                             g_0_xzz_0_xzzzzzzz_0, \
                             g_0_xzz_0_yyyyyyyy_0, \
                             g_0_xzz_0_yyyyyyyz_0, \
                             g_0_xzz_0_yyyyyyzz_0, \
                             g_0_xzz_0_yyyyyzzz_0, \
                             g_0_xzz_0_yyyyzzzz_0, \
                             g_0_xzz_0_yyyzzzzz_0, \
                             g_0_xzz_0_yyzzzzzz_0, \
                             g_0_xzz_0_yzzzzzzz_0, \
                             g_0_xzz_0_zzzzzzzz_0, \
                             g_0_zz_0_xxxxxxx_1,   \
                             g_0_zz_0_xxxxxxxx_0,  \
                             g_0_zz_0_xxxxxxxx_1,  \
                             g_0_zz_0_xxxxxxxy_0,  \
                             g_0_zz_0_xxxxxxxy_1,  \
                             g_0_zz_0_xxxxxxxz_0,  \
                             g_0_zz_0_xxxxxxxz_1,  \
                             g_0_zz_0_xxxxxxy_1,   \
                             g_0_zz_0_xxxxxxyy_0,  \
                             g_0_zz_0_xxxxxxyy_1,  \
                             g_0_zz_0_xxxxxxyz_0,  \
                             g_0_zz_0_xxxxxxyz_1,  \
                             g_0_zz_0_xxxxxxz_1,   \
                             g_0_zz_0_xxxxxxzz_0,  \
                             g_0_zz_0_xxxxxxzz_1,  \
                             g_0_zz_0_xxxxxyy_1,   \
                             g_0_zz_0_xxxxxyyy_0,  \
                             g_0_zz_0_xxxxxyyy_1,  \
                             g_0_zz_0_xxxxxyyz_0,  \
                             g_0_zz_0_xxxxxyyz_1,  \
                             g_0_zz_0_xxxxxyz_1,   \
                             g_0_zz_0_xxxxxyzz_0,  \
                             g_0_zz_0_xxxxxyzz_1,  \
                             g_0_zz_0_xxxxxzz_1,   \
                             g_0_zz_0_xxxxxzzz_0,  \
                             g_0_zz_0_xxxxxzzz_1,  \
                             g_0_zz_0_xxxxyyy_1,   \
                             g_0_zz_0_xxxxyyyy_0,  \
                             g_0_zz_0_xxxxyyyy_1,  \
                             g_0_zz_0_xxxxyyyz_0,  \
                             g_0_zz_0_xxxxyyyz_1,  \
                             g_0_zz_0_xxxxyyz_1,   \
                             g_0_zz_0_xxxxyyzz_0,  \
                             g_0_zz_0_xxxxyyzz_1,  \
                             g_0_zz_0_xxxxyzz_1,   \
                             g_0_zz_0_xxxxyzzz_0,  \
                             g_0_zz_0_xxxxyzzz_1,  \
                             g_0_zz_0_xxxxzzz_1,   \
                             g_0_zz_0_xxxxzzzz_0,  \
                             g_0_zz_0_xxxxzzzz_1,  \
                             g_0_zz_0_xxxyyyy_1,   \
                             g_0_zz_0_xxxyyyyy_0,  \
                             g_0_zz_0_xxxyyyyy_1,  \
                             g_0_zz_0_xxxyyyyz_0,  \
                             g_0_zz_0_xxxyyyyz_1,  \
                             g_0_zz_0_xxxyyyz_1,   \
                             g_0_zz_0_xxxyyyzz_0,  \
                             g_0_zz_0_xxxyyyzz_1,  \
                             g_0_zz_0_xxxyyzz_1,   \
                             g_0_zz_0_xxxyyzzz_0,  \
                             g_0_zz_0_xxxyyzzz_1,  \
                             g_0_zz_0_xxxyzzz_1,   \
                             g_0_zz_0_xxxyzzzz_0,  \
                             g_0_zz_0_xxxyzzzz_1,  \
                             g_0_zz_0_xxxzzzz_1,   \
                             g_0_zz_0_xxxzzzzz_0,  \
                             g_0_zz_0_xxxzzzzz_1,  \
                             g_0_zz_0_xxyyyyy_1,   \
                             g_0_zz_0_xxyyyyyy_0,  \
                             g_0_zz_0_xxyyyyyy_1,  \
                             g_0_zz_0_xxyyyyyz_0,  \
                             g_0_zz_0_xxyyyyyz_1,  \
                             g_0_zz_0_xxyyyyz_1,   \
                             g_0_zz_0_xxyyyyzz_0,  \
                             g_0_zz_0_xxyyyyzz_1,  \
                             g_0_zz_0_xxyyyzz_1,   \
                             g_0_zz_0_xxyyyzzz_0,  \
                             g_0_zz_0_xxyyyzzz_1,  \
                             g_0_zz_0_xxyyzzz_1,   \
                             g_0_zz_0_xxyyzzzz_0,  \
                             g_0_zz_0_xxyyzzzz_1,  \
                             g_0_zz_0_xxyzzzz_1,   \
                             g_0_zz_0_xxyzzzzz_0,  \
                             g_0_zz_0_xxyzzzzz_1,  \
                             g_0_zz_0_xxzzzzz_1,   \
                             g_0_zz_0_xxzzzzzz_0,  \
                             g_0_zz_0_xxzzzzzz_1,  \
                             g_0_zz_0_xyyyyyy_1,   \
                             g_0_zz_0_xyyyyyyy_0,  \
                             g_0_zz_0_xyyyyyyy_1,  \
                             g_0_zz_0_xyyyyyyz_0,  \
                             g_0_zz_0_xyyyyyyz_1,  \
                             g_0_zz_0_xyyyyyz_1,   \
                             g_0_zz_0_xyyyyyzz_0,  \
                             g_0_zz_0_xyyyyyzz_1,  \
                             g_0_zz_0_xyyyyzz_1,   \
                             g_0_zz_0_xyyyyzzz_0,  \
                             g_0_zz_0_xyyyyzzz_1,  \
                             g_0_zz_0_xyyyzzz_1,   \
                             g_0_zz_0_xyyyzzzz_0,  \
                             g_0_zz_0_xyyyzzzz_1,  \
                             g_0_zz_0_xyyzzzz_1,   \
                             g_0_zz_0_xyyzzzzz_0,  \
                             g_0_zz_0_xyyzzzzz_1,  \
                             g_0_zz_0_xyzzzzz_1,   \
                             g_0_zz_0_xyzzzzzz_0,  \
                             g_0_zz_0_xyzzzzzz_1,  \
                             g_0_zz_0_xzzzzzz_1,   \
                             g_0_zz_0_xzzzzzzz_0,  \
                             g_0_zz_0_xzzzzzzz_1,  \
                             g_0_zz_0_yyyyyyy_1,   \
                             g_0_zz_0_yyyyyyyy_0,  \
                             g_0_zz_0_yyyyyyyy_1,  \
                             g_0_zz_0_yyyyyyyz_0,  \
                             g_0_zz_0_yyyyyyyz_1,  \
                             g_0_zz_0_yyyyyyz_1,   \
                             g_0_zz_0_yyyyyyzz_0,  \
                             g_0_zz_0_yyyyyyzz_1,  \
                             g_0_zz_0_yyyyyzz_1,   \
                             g_0_zz_0_yyyyyzzz_0,  \
                             g_0_zz_0_yyyyyzzz_1,  \
                             g_0_zz_0_yyyyzzz_1,   \
                             g_0_zz_0_yyyyzzzz_0,  \
                             g_0_zz_0_yyyyzzzz_1,  \
                             g_0_zz_0_yyyzzzz_1,   \
                             g_0_zz_0_yyyzzzzz_0,  \
                             g_0_zz_0_yyyzzzzz_1,  \
                             g_0_zz_0_yyzzzzz_1,   \
                             g_0_zz_0_yyzzzzzz_0,  \
                             g_0_zz_0_yyzzzzzz_1,  \
                             g_0_zz_0_yzzzzzz_1,   \
                             g_0_zz_0_yzzzzzzz_0,  \
                             g_0_zz_0_yzzzzzzz_1,  \
                             g_0_zz_0_zzzzzzz_1,   \
                             g_0_zz_0_zzzzzzzz_0,  \
                             g_0_zz_0_zzzzzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xxxxxxxx_0[i] = 8.0 * g_0_zz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxxx_0[i] * pb_x + g_0_zz_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxxxy_0[i] = 7.0 * g_0_zz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxxy_0[i] * pb_x + g_0_zz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxxxz_0[i] = 7.0 * g_0_zz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxxz_0[i] * pb_x + g_0_zz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxxyy_0[i] = 6.0 * g_0_zz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxyy_0[i] * pb_x + g_0_zz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxxyz_0[i] = 6.0 * g_0_zz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxyz_0[i] * pb_x + g_0_zz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxxzz_0[i] = 6.0 * g_0_zz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxzz_0[i] * pb_x + g_0_zz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxyyy_0[i] = 5.0 * g_0_zz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyyy_0[i] * pb_x + g_0_zz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxyyz_0[i] = 5.0 * g_0_zz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyyz_0[i] * pb_x + g_0_zz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxyzz_0[i] = 5.0 * g_0_zz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyzz_0[i] * pb_x + g_0_zz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxzzz_0[i] = 5.0 * g_0_zz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxzzz_0[i] * pb_x + g_0_zz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyyyy_0[i] = 4.0 * g_0_zz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyyy_0[i] * pb_x + g_0_zz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyyyz_0[i] = 4.0 * g_0_zz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyyz_0[i] * pb_x + g_0_zz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyyzz_0[i] = 4.0 * g_0_zz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyzz_0[i] * pb_x + g_0_zz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyzzz_0[i] = 4.0 * g_0_zz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyzzz_0[i] * pb_x + g_0_zz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxzzzz_0[i] = 4.0 * g_0_zz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxzzzz_0[i] * pb_x + g_0_zz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyyyy_0[i] = 3.0 * g_0_zz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyyy_0[i] * pb_x + g_0_zz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyyyz_0[i] = 3.0 * g_0_zz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyyz_0[i] * pb_x + g_0_zz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyyzz_0[i] = 3.0 * g_0_zz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyzz_0[i] * pb_x + g_0_zz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyzzz_0[i] = 3.0 * g_0_zz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyzzz_0[i] * pb_x + g_0_zz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyzzzz_0[i] = 3.0 * g_0_zz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzzzz_0[i] * pb_x + g_0_zz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxzzzzz_0[i] = 3.0 * g_0_zz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxzzzzz_0[i] * pb_x + g_0_zz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyyyy_0[i] = 2.0 * g_0_zz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyyy_0[i] * pb_x + g_0_zz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyyyz_0[i] = 2.0 * g_0_zz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyyz_0[i] * pb_x + g_0_zz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyyzz_0[i] = 2.0 * g_0_zz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyzz_0[i] * pb_x + g_0_zz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyzzz_0[i] = 2.0 * g_0_zz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyzzz_0[i] * pb_x + g_0_zz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyzzzz_0[i] = 2.0 * g_0_zz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzzzz_0[i] * pb_x + g_0_zz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyzzzzz_0[i] = 2.0 * g_0_zz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzzzz_0[i] * pb_x + g_0_zz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxzzzzzz_0[i] = 2.0 * g_0_zz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzzzzzz_0[i] * pb_x + g_0_zz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyyyy_0[i] = g_0_zz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyyy_0[i] * pb_x + g_0_zz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyyyz_0[i] = g_0_zz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyyz_0[i] * pb_x + g_0_zz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyyzz_0[i] = g_0_zz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyzz_0[i] * pb_x + g_0_zz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyzzz_0[i] = g_0_zz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyzzz_0[i] * pb_x + g_0_zz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyzzzz_0[i] = g_0_zz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzzzz_0[i] * pb_x + g_0_zz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyzzzzz_0[i] = g_0_zz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzzzz_0[i] * pb_x + g_0_zz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyzzzzzz_0[i] = g_0_zz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzzzz_0[i] * pb_x + g_0_zz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xzzzzzzz_0[i] = g_0_zz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzzzzzz_0[i] * pb_x + g_0_zz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyyyy_0[i] = g_0_zz_0_yyyyyyyy_0[i] * pb_x + g_0_zz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyyyz_0[i] = g_0_zz_0_yyyyyyyz_0[i] * pb_x + g_0_zz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyyzz_0[i] = g_0_zz_0_yyyyyyzz_0[i] * pb_x + g_0_zz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyzzz_0[i] = g_0_zz_0_yyyyyzzz_0[i] * pb_x + g_0_zz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyzzzz_0[i] = g_0_zz_0_yyyyzzzz_0[i] * pb_x + g_0_zz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyzzzzz_0[i] = g_0_zz_0_yyyzzzzz_0[i] * pb_x + g_0_zz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyzzzzzz_0[i] = g_0_zz_0_yyzzzzzz_0[i] * pb_x + g_0_zz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yzzzzzzz_0[i] = g_0_zz_0_yzzzzzzz_0[i] * pb_x + g_0_zz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_zzzzzzzz_0[i] = g_0_zz_0_zzzzzzzz_0[i] * pb_x + g_0_zz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 270-315 components of targeted buffer : SFSL

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

#pragma omp simd aligned(g_0_y_0_xxxxxxxx_0,       \
                             g_0_y_0_xxxxxxxx_1,   \
                             g_0_y_0_xxxxxxxy_0,   \
                             g_0_y_0_xxxxxxxy_1,   \
                             g_0_y_0_xxxxxxxz_0,   \
                             g_0_y_0_xxxxxxxz_1,   \
                             g_0_y_0_xxxxxxyy_0,   \
                             g_0_y_0_xxxxxxyy_1,   \
                             g_0_y_0_xxxxxxyz_0,   \
                             g_0_y_0_xxxxxxyz_1,   \
                             g_0_y_0_xxxxxxzz_0,   \
                             g_0_y_0_xxxxxxzz_1,   \
                             g_0_y_0_xxxxxyyy_0,   \
                             g_0_y_0_xxxxxyyy_1,   \
                             g_0_y_0_xxxxxyyz_0,   \
                             g_0_y_0_xxxxxyyz_1,   \
                             g_0_y_0_xxxxxyzz_0,   \
                             g_0_y_0_xxxxxyzz_1,   \
                             g_0_y_0_xxxxxzzz_0,   \
                             g_0_y_0_xxxxxzzz_1,   \
                             g_0_y_0_xxxxyyyy_0,   \
                             g_0_y_0_xxxxyyyy_1,   \
                             g_0_y_0_xxxxyyyz_0,   \
                             g_0_y_0_xxxxyyyz_1,   \
                             g_0_y_0_xxxxyyzz_0,   \
                             g_0_y_0_xxxxyyzz_1,   \
                             g_0_y_0_xxxxyzzz_0,   \
                             g_0_y_0_xxxxyzzz_1,   \
                             g_0_y_0_xxxxzzzz_0,   \
                             g_0_y_0_xxxxzzzz_1,   \
                             g_0_y_0_xxxyyyyy_0,   \
                             g_0_y_0_xxxyyyyy_1,   \
                             g_0_y_0_xxxyyyyz_0,   \
                             g_0_y_0_xxxyyyyz_1,   \
                             g_0_y_0_xxxyyyzz_0,   \
                             g_0_y_0_xxxyyyzz_1,   \
                             g_0_y_0_xxxyyzzz_0,   \
                             g_0_y_0_xxxyyzzz_1,   \
                             g_0_y_0_xxxyzzzz_0,   \
                             g_0_y_0_xxxyzzzz_1,   \
                             g_0_y_0_xxxzzzzz_0,   \
                             g_0_y_0_xxxzzzzz_1,   \
                             g_0_y_0_xxyyyyyy_0,   \
                             g_0_y_0_xxyyyyyy_1,   \
                             g_0_y_0_xxyyyyyz_0,   \
                             g_0_y_0_xxyyyyyz_1,   \
                             g_0_y_0_xxyyyyzz_0,   \
                             g_0_y_0_xxyyyyzz_1,   \
                             g_0_y_0_xxyyyzzz_0,   \
                             g_0_y_0_xxyyyzzz_1,   \
                             g_0_y_0_xxyyzzzz_0,   \
                             g_0_y_0_xxyyzzzz_1,   \
                             g_0_y_0_xxyzzzzz_0,   \
                             g_0_y_0_xxyzzzzz_1,   \
                             g_0_y_0_xxzzzzzz_0,   \
                             g_0_y_0_xxzzzzzz_1,   \
                             g_0_y_0_xyyyyyyy_0,   \
                             g_0_y_0_xyyyyyyy_1,   \
                             g_0_y_0_xyyyyyyz_0,   \
                             g_0_y_0_xyyyyyyz_1,   \
                             g_0_y_0_xyyyyyzz_0,   \
                             g_0_y_0_xyyyyyzz_1,   \
                             g_0_y_0_xyyyyzzz_0,   \
                             g_0_y_0_xyyyyzzz_1,   \
                             g_0_y_0_xyyyzzzz_0,   \
                             g_0_y_0_xyyyzzzz_1,   \
                             g_0_y_0_xyyzzzzz_0,   \
                             g_0_y_0_xyyzzzzz_1,   \
                             g_0_y_0_xyzzzzzz_0,   \
                             g_0_y_0_xyzzzzzz_1,   \
                             g_0_y_0_xzzzzzzz_0,   \
                             g_0_y_0_xzzzzzzz_1,   \
                             g_0_y_0_yyyyyyyy_0,   \
                             g_0_y_0_yyyyyyyy_1,   \
                             g_0_y_0_yyyyyyyz_0,   \
                             g_0_y_0_yyyyyyyz_1,   \
                             g_0_y_0_yyyyyyzz_0,   \
                             g_0_y_0_yyyyyyzz_1,   \
                             g_0_y_0_yyyyyzzz_0,   \
                             g_0_y_0_yyyyyzzz_1,   \
                             g_0_y_0_yyyyzzzz_0,   \
                             g_0_y_0_yyyyzzzz_1,   \
                             g_0_y_0_yyyzzzzz_0,   \
                             g_0_y_0_yyyzzzzz_1,   \
                             g_0_y_0_yyzzzzzz_0,   \
                             g_0_y_0_yyzzzzzz_1,   \
                             g_0_y_0_yzzzzzzz_0,   \
                             g_0_y_0_yzzzzzzz_1,   \
                             g_0_y_0_zzzzzzzz_0,   \
                             g_0_y_0_zzzzzzzz_1,   \
                             g_0_yy_0_xxxxxxx_1,   \
                             g_0_yy_0_xxxxxxxx_0,  \
                             g_0_yy_0_xxxxxxxx_1,  \
                             g_0_yy_0_xxxxxxxy_0,  \
                             g_0_yy_0_xxxxxxxy_1,  \
                             g_0_yy_0_xxxxxxxz_0,  \
                             g_0_yy_0_xxxxxxxz_1,  \
                             g_0_yy_0_xxxxxxy_1,   \
                             g_0_yy_0_xxxxxxyy_0,  \
                             g_0_yy_0_xxxxxxyy_1,  \
                             g_0_yy_0_xxxxxxyz_0,  \
                             g_0_yy_0_xxxxxxyz_1,  \
                             g_0_yy_0_xxxxxxz_1,   \
                             g_0_yy_0_xxxxxxzz_0,  \
                             g_0_yy_0_xxxxxxzz_1,  \
                             g_0_yy_0_xxxxxyy_1,   \
                             g_0_yy_0_xxxxxyyy_0,  \
                             g_0_yy_0_xxxxxyyy_1,  \
                             g_0_yy_0_xxxxxyyz_0,  \
                             g_0_yy_0_xxxxxyyz_1,  \
                             g_0_yy_0_xxxxxyz_1,   \
                             g_0_yy_0_xxxxxyzz_0,  \
                             g_0_yy_0_xxxxxyzz_1,  \
                             g_0_yy_0_xxxxxzz_1,   \
                             g_0_yy_0_xxxxxzzz_0,  \
                             g_0_yy_0_xxxxxzzz_1,  \
                             g_0_yy_0_xxxxyyy_1,   \
                             g_0_yy_0_xxxxyyyy_0,  \
                             g_0_yy_0_xxxxyyyy_1,  \
                             g_0_yy_0_xxxxyyyz_0,  \
                             g_0_yy_0_xxxxyyyz_1,  \
                             g_0_yy_0_xxxxyyz_1,   \
                             g_0_yy_0_xxxxyyzz_0,  \
                             g_0_yy_0_xxxxyyzz_1,  \
                             g_0_yy_0_xxxxyzz_1,   \
                             g_0_yy_0_xxxxyzzz_0,  \
                             g_0_yy_0_xxxxyzzz_1,  \
                             g_0_yy_0_xxxxzzz_1,   \
                             g_0_yy_0_xxxxzzzz_0,  \
                             g_0_yy_0_xxxxzzzz_1,  \
                             g_0_yy_0_xxxyyyy_1,   \
                             g_0_yy_0_xxxyyyyy_0,  \
                             g_0_yy_0_xxxyyyyy_1,  \
                             g_0_yy_0_xxxyyyyz_0,  \
                             g_0_yy_0_xxxyyyyz_1,  \
                             g_0_yy_0_xxxyyyz_1,   \
                             g_0_yy_0_xxxyyyzz_0,  \
                             g_0_yy_0_xxxyyyzz_1,  \
                             g_0_yy_0_xxxyyzz_1,   \
                             g_0_yy_0_xxxyyzzz_0,  \
                             g_0_yy_0_xxxyyzzz_1,  \
                             g_0_yy_0_xxxyzzz_1,   \
                             g_0_yy_0_xxxyzzzz_0,  \
                             g_0_yy_0_xxxyzzzz_1,  \
                             g_0_yy_0_xxxzzzz_1,   \
                             g_0_yy_0_xxxzzzzz_0,  \
                             g_0_yy_0_xxxzzzzz_1,  \
                             g_0_yy_0_xxyyyyy_1,   \
                             g_0_yy_0_xxyyyyyy_0,  \
                             g_0_yy_0_xxyyyyyy_1,  \
                             g_0_yy_0_xxyyyyyz_0,  \
                             g_0_yy_0_xxyyyyyz_1,  \
                             g_0_yy_0_xxyyyyz_1,   \
                             g_0_yy_0_xxyyyyzz_0,  \
                             g_0_yy_0_xxyyyyzz_1,  \
                             g_0_yy_0_xxyyyzz_1,   \
                             g_0_yy_0_xxyyyzzz_0,  \
                             g_0_yy_0_xxyyyzzz_1,  \
                             g_0_yy_0_xxyyzzz_1,   \
                             g_0_yy_0_xxyyzzzz_0,  \
                             g_0_yy_0_xxyyzzzz_1,  \
                             g_0_yy_0_xxyzzzz_1,   \
                             g_0_yy_0_xxyzzzzz_0,  \
                             g_0_yy_0_xxyzzzzz_1,  \
                             g_0_yy_0_xxzzzzz_1,   \
                             g_0_yy_0_xxzzzzzz_0,  \
                             g_0_yy_0_xxzzzzzz_1,  \
                             g_0_yy_0_xyyyyyy_1,   \
                             g_0_yy_0_xyyyyyyy_0,  \
                             g_0_yy_0_xyyyyyyy_1,  \
                             g_0_yy_0_xyyyyyyz_0,  \
                             g_0_yy_0_xyyyyyyz_1,  \
                             g_0_yy_0_xyyyyyz_1,   \
                             g_0_yy_0_xyyyyyzz_0,  \
                             g_0_yy_0_xyyyyyzz_1,  \
                             g_0_yy_0_xyyyyzz_1,   \
                             g_0_yy_0_xyyyyzzz_0,  \
                             g_0_yy_0_xyyyyzzz_1,  \
                             g_0_yy_0_xyyyzzz_1,   \
                             g_0_yy_0_xyyyzzzz_0,  \
                             g_0_yy_0_xyyyzzzz_1,  \
                             g_0_yy_0_xyyzzzz_1,   \
                             g_0_yy_0_xyyzzzzz_0,  \
                             g_0_yy_0_xyyzzzzz_1,  \
                             g_0_yy_0_xyzzzzz_1,   \
                             g_0_yy_0_xyzzzzzz_0,  \
                             g_0_yy_0_xyzzzzzz_1,  \
                             g_0_yy_0_xzzzzzz_1,   \
                             g_0_yy_0_xzzzzzzz_0,  \
                             g_0_yy_0_xzzzzzzz_1,  \
                             g_0_yy_0_yyyyyyy_1,   \
                             g_0_yy_0_yyyyyyyy_0,  \
                             g_0_yy_0_yyyyyyyy_1,  \
                             g_0_yy_0_yyyyyyyz_0,  \
                             g_0_yy_0_yyyyyyyz_1,  \
                             g_0_yy_0_yyyyyyz_1,   \
                             g_0_yy_0_yyyyyyzz_0,  \
                             g_0_yy_0_yyyyyyzz_1,  \
                             g_0_yy_0_yyyyyzz_1,   \
                             g_0_yy_0_yyyyyzzz_0,  \
                             g_0_yy_0_yyyyyzzz_1,  \
                             g_0_yy_0_yyyyzzz_1,   \
                             g_0_yy_0_yyyyzzzz_0,  \
                             g_0_yy_0_yyyyzzzz_1,  \
                             g_0_yy_0_yyyzzzz_1,   \
                             g_0_yy_0_yyyzzzzz_0,  \
                             g_0_yy_0_yyyzzzzz_1,  \
                             g_0_yy_0_yyzzzzz_1,   \
                             g_0_yy_0_yyzzzzzz_0,  \
                             g_0_yy_0_yyzzzzzz_1,  \
                             g_0_yy_0_yzzzzzz_1,   \
                             g_0_yy_0_yzzzzzzz_0,  \
                             g_0_yy_0_yzzzzzzz_1,  \
                             g_0_yy_0_zzzzzzz_1,   \
                             g_0_yy_0_zzzzzzzz_0,  \
                             g_0_yy_0_zzzzzzzz_1,  \
                             g_0_yyy_0_xxxxxxxx_0, \
                             g_0_yyy_0_xxxxxxxy_0, \
                             g_0_yyy_0_xxxxxxxz_0, \
                             g_0_yyy_0_xxxxxxyy_0, \
                             g_0_yyy_0_xxxxxxyz_0, \
                             g_0_yyy_0_xxxxxxzz_0, \
                             g_0_yyy_0_xxxxxyyy_0, \
                             g_0_yyy_0_xxxxxyyz_0, \
                             g_0_yyy_0_xxxxxyzz_0, \
                             g_0_yyy_0_xxxxxzzz_0, \
                             g_0_yyy_0_xxxxyyyy_0, \
                             g_0_yyy_0_xxxxyyyz_0, \
                             g_0_yyy_0_xxxxyyzz_0, \
                             g_0_yyy_0_xxxxyzzz_0, \
                             g_0_yyy_0_xxxxzzzz_0, \
                             g_0_yyy_0_xxxyyyyy_0, \
                             g_0_yyy_0_xxxyyyyz_0, \
                             g_0_yyy_0_xxxyyyzz_0, \
                             g_0_yyy_0_xxxyyzzz_0, \
                             g_0_yyy_0_xxxyzzzz_0, \
                             g_0_yyy_0_xxxzzzzz_0, \
                             g_0_yyy_0_xxyyyyyy_0, \
                             g_0_yyy_0_xxyyyyyz_0, \
                             g_0_yyy_0_xxyyyyzz_0, \
                             g_0_yyy_0_xxyyyzzz_0, \
                             g_0_yyy_0_xxyyzzzz_0, \
                             g_0_yyy_0_xxyzzzzz_0, \
                             g_0_yyy_0_xxzzzzzz_0, \
                             g_0_yyy_0_xyyyyyyy_0, \
                             g_0_yyy_0_xyyyyyyz_0, \
                             g_0_yyy_0_xyyyyyzz_0, \
                             g_0_yyy_0_xyyyyzzz_0, \
                             g_0_yyy_0_xyyyzzzz_0, \
                             g_0_yyy_0_xyyzzzzz_0, \
                             g_0_yyy_0_xyzzzzzz_0, \
                             g_0_yyy_0_xzzzzzzz_0, \
                             g_0_yyy_0_yyyyyyyy_0, \
                             g_0_yyy_0_yyyyyyyz_0, \
                             g_0_yyy_0_yyyyyyzz_0, \
                             g_0_yyy_0_yyyyyzzz_0, \
                             g_0_yyy_0_yyyyzzzz_0, \
                             g_0_yyy_0_yyyzzzzz_0, \
                             g_0_yyy_0_yyzzzzzz_0, \
                             g_0_yyy_0_yzzzzzzz_0, \
                             g_0_yyy_0_zzzzzzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xxxxxxxx_0[i] = 2.0 * g_0_y_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_yy_0_xxxxxxxx_0[i] * pb_y +
                                  g_0_yy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxxxy_0[i] = 2.0 * g_0_y_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_yy_0_xxxxxxx_1[i] * fi_abcd_0 +
                                  g_0_yy_0_xxxxxxxy_0[i] * pb_y + g_0_yy_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxxxz_0[i] = 2.0 * g_0_y_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxxxz_0[i] * pb_y +
                                  g_0_yy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxxyy_0[i] = 2.0 * g_0_y_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_yy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxyy_0[i] * pb_y + g_0_yy_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxxyz_0[i] = 2.0 * g_0_y_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxxz_1[i] * fi_abcd_0 +
                                  g_0_yy_0_xxxxxxyz_0[i] * pb_y + g_0_yy_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxxzz_0[i] = 2.0 * g_0_y_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxxzz_0[i] * pb_y +
                                  g_0_yy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxyyy_0[i] = 2.0 * g_0_y_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxyyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_yy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyyy_0[i] * pb_y + g_0_yy_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxyyz_0[i] = 2.0 * g_0_y_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyyz_0[i] * pb_y + g_0_yy_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxyzz_0[i] = 2.0 * g_0_y_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxyzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxzz_1[i] * fi_abcd_0 +
                                  g_0_yy_0_xxxxxyzz_0[i] * pb_y + g_0_yy_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxzzz_0[i] = 2.0 * g_0_y_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxzzz_0[i] * pb_y +
                                  g_0_yy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyyyy_0[i] = 2.0 * g_0_y_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyyyy_1[i] * fti_ab_0 +
                                  4.0 * g_0_yy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyyy_0[i] * pb_y + g_0_yy_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyyyz_0[i] = 2.0 * g_0_y_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyyz_0[i] * pb_y + g_0_yy_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyyzz_0[i] = 2.0 * g_0_y_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyzz_0[i] * pb_y + g_0_yy_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyzzz_0[i] = 2.0 * g_0_y_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxzzz_1[i] * fi_abcd_0 +
                                  g_0_yy_0_xxxxyzzz_0[i] * pb_y + g_0_yy_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxzzzz_0[i] = 2.0 * g_0_y_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxzzzz_0[i] * pb_y +
                                  g_0_yy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyyyy_0[i] = 2.0 * g_0_y_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyyyy_1[i] * fti_ab_0 +
                                  5.0 * g_0_yy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyyy_0[i] * pb_y + g_0_yy_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyyyz_0[i] = 2.0 * g_0_y_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyyyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyyz_0[i] * pb_y + g_0_yy_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyyzz_0[i] = 2.0 * g_0_y_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyzz_0[i] * pb_y + g_0_yy_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyzzz_0[i] = 2.0 * g_0_y_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyzzz_0[i] * pb_y + g_0_yy_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyzzzz_0[i] = 2.0 * g_0_y_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxzzzz_1[i] * fi_abcd_0 +
                                  g_0_yy_0_xxxyzzzz_0[i] * pb_y + g_0_yy_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxzzzzz_0[i] = 2.0 * g_0_y_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxzzzzz_0[i] * pb_y +
                                  g_0_yy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyyyy_0[i] = 2.0 * g_0_y_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyyyy_1[i] * fti_ab_0 +
                                  6.0 * g_0_yy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyyy_0[i] * pb_y + g_0_yy_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyyyz_0[i] = 2.0 * g_0_y_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyyyz_1[i] * fti_ab_0 +
                                  5.0 * g_0_yy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyyz_0[i] * pb_y + g_0_yy_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyyzz_0[i] = 2.0 * g_0_y_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyyzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyzz_0[i] * pb_y + g_0_yy_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyzzz_0[i] = 2.0 * g_0_y_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyzzz_0[i] * pb_y + g_0_yy_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyzzzz_0[i] = 2.0 * g_0_y_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzzzz_0[i] * pb_y + g_0_yy_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyzzzzz_0[i] = 2.0 * g_0_y_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxzzzzz_1[i] * fi_abcd_0 +
                                  g_0_yy_0_xxyzzzzz_0[i] * pb_y + g_0_yy_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxzzzzzz_0[i] = 2.0 * g_0_y_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxzzzzzz_0[i] * pb_y +
                                  g_0_yy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyyyy_0[i] = 2.0 * g_0_y_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyyyy_1[i] * fti_ab_0 +
                                  7.0 * g_0_yy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyyy_0[i] * pb_y + g_0_yy_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyyyz_0[i] = 2.0 * g_0_y_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyyyz_1[i] * fti_ab_0 +
                                  6.0 * g_0_yy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyyz_0[i] * pb_y + g_0_yy_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyyzz_0[i] = 2.0 * g_0_y_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyyzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_yy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyzz_0[i] * pb_y + g_0_yy_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyzzz_0[i] = 2.0 * g_0_y_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyzzz_0[i] * pb_y + g_0_yy_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyzzzz_0[i] = 2.0 * g_0_y_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyzzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzzzz_0[i] * pb_y + g_0_yy_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyzzzzz_0[i] = 2.0 * g_0_y_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyzzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzzzz_0[i] * pb_y + g_0_yy_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyzzzzzz_0[i] = 2.0 * g_0_y_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzzzzz_1[i] * fi_abcd_0 +
                                  g_0_yy_0_xyzzzzzz_0[i] * pb_y + g_0_yy_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xzzzzzzz_0[i] = 2.0 * g_0_y_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzzzzzz_0[i] * pb_y +
                                  g_0_yy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyyyy_0[i] = 2.0 * g_0_y_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyyyy_1[i] * fti_ab_0 +
                                  8.0 * g_0_yy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyyy_0[i] * pb_y + g_0_yy_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyyyz_0[i] = 2.0 * g_0_y_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyyyz_1[i] * fti_ab_0 +
                                  7.0 * g_0_yy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyyz_0[i] * pb_y + g_0_yy_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyyzz_0[i] = 2.0 * g_0_y_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyyzz_1[i] * fti_ab_0 +
                                  6.0 * g_0_yy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyzz_0[i] * pb_y + g_0_yy_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyzzz_0[i] = 2.0 * g_0_y_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_yy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyzzz_0[i] * pb_y + g_0_yy_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyzzzz_0[i] = 2.0 * g_0_y_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyzzzz_0[i] * pb_y + g_0_yy_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyzzzzz_0[i] = 2.0 * g_0_y_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyzzzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyzzzzz_0[i] * pb_y + g_0_yy_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyzzzzzz_0[i] = 2.0 * g_0_y_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyzzzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzzzzzz_0[i] * pb_y + g_0_yy_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yzzzzzzz_0[i] = 2.0 * g_0_y_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzzzzz_1[i] * fi_abcd_0 +
                                  g_0_yy_0_yzzzzzzz_0[i] * pb_y + g_0_yy_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_zzzzzzzz_0[i] = 2.0 * g_0_y_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzzzzzz_0[i] * pb_y +
                                  g_0_yy_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 315-360 components of targeted buffer : SFSL

    auto g_0_yyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sfsl + 315);

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

#pragma omp simd aligned(g_0_yy_0_xxxxxxx_1,       \
                             g_0_yy_0_xxxxxxxx_0,  \
                             g_0_yy_0_xxxxxxxx_1,  \
                             g_0_yy_0_xxxxxxxy_0,  \
                             g_0_yy_0_xxxxxxxy_1,  \
                             g_0_yy_0_xxxxxxxz_0,  \
                             g_0_yy_0_xxxxxxxz_1,  \
                             g_0_yy_0_xxxxxxy_1,   \
                             g_0_yy_0_xxxxxxyy_0,  \
                             g_0_yy_0_xxxxxxyy_1,  \
                             g_0_yy_0_xxxxxxyz_0,  \
                             g_0_yy_0_xxxxxxyz_1,  \
                             g_0_yy_0_xxxxxxz_1,   \
                             g_0_yy_0_xxxxxxzz_0,  \
                             g_0_yy_0_xxxxxxzz_1,  \
                             g_0_yy_0_xxxxxyy_1,   \
                             g_0_yy_0_xxxxxyyy_0,  \
                             g_0_yy_0_xxxxxyyy_1,  \
                             g_0_yy_0_xxxxxyyz_0,  \
                             g_0_yy_0_xxxxxyyz_1,  \
                             g_0_yy_0_xxxxxyz_1,   \
                             g_0_yy_0_xxxxxyzz_0,  \
                             g_0_yy_0_xxxxxyzz_1,  \
                             g_0_yy_0_xxxxxzz_1,   \
                             g_0_yy_0_xxxxxzzz_0,  \
                             g_0_yy_0_xxxxxzzz_1,  \
                             g_0_yy_0_xxxxyyy_1,   \
                             g_0_yy_0_xxxxyyyy_0,  \
                             g_0_yy_0_xxxxyyyy_1,  \
                             g_0_yy_0_xxxxyyyz_0,  \
                             g_0_yy_0_xxxxyyyz_1,  \
                             g_0_yy_0_xxxxyyz_1,   \
                             g_0_yy_0_xxxxyyzz_0,  \
                             g_0_yy_0_xxxxyyzz_1,  \
                             g_0_yy_0_xxxxyzz_1,   \
                             g_0_yy_0_xxxxyzzz_0,  \
                             g_0_yy_0_xxxxyzzz_1,  \
                             g_0_yy_0_xxxxzzz_1,   \
                             g_0_yy_0_xxxxzzzz_0,  \
                             g_0_yy_0_xxxxzzzz_1,  \
                             g_0_yy_0_xxxyyyy_1,   \
                             g_0_yy_0_xxxyyyyy_0,  \
                             g_0_yy_0_xxxyyyyy_1,  \
                             g_0_yy_0_xxxyyyyz_0,  \
                             g_0_yy_0_xxxyyyyz_1,  \
                             g_0_yy_0_xxxyyyz_1,   \
                             g_0_yy_0_xxxyyyzz_0,  \
                             g_0_yy_0_xxxyyyzz_1,  \
                             g_0_yy_0_xxxyyzz_1,   \
                             g_0_yy_0_xxxyyzzz_0,  \
                             g_0_yy_0_xxxyyzzz_1,  \
                             g_0_yy_0_xxxyzzz_1,   \
                             g_0_yy_0_xxxyzzzz_0,  \
                             g_0_yy_0_xxxyzzzz_1,  \
                             g_0_yy_0_xxxzzzz_1,   \
                             g_0_yy_0_xxxzzzzz_0,  \
                             g_0_yy_0_xxxzzzzz_1,  \
                             g_0_yy_0_xxyyyyy_1,   \
                             g_0_yy_0_xxyyyyyy_0,  \
                             g_0_yy_0_xxyyyyyy_1,  \
                             g_0_yy_0_xxyyyyyz_0,  \
                             g_0_yy_0_xxyyyyyz_1,  \
                             g_0_yy_0_xxyyyyz_1,   \
                             g_0_yy_0_xxyyyyzz_0,  \
                             g_0_yy_0_xxyyyyzz_1,  \
                             g_0_yy_0_xxyyyzz_1,   \
                             g_0_yy_0_xxyyyzzz_0,  \
                             g_0_yy_0_xxyyyzzz_1,  \
                             g_0_yy_0_xxyyzzz_1,   \
                             g_0_yy_0_xxyyzzzz_0,  \
                             g_0_yy_0_xxyyzzzz_1,  \
                             g_0_yy_0_xxyzzzz_1,   \
                             g_0_yy_0_xxyzzzzz_0,  \
                             g_0_yy_0_xxyzzzzz_1,  \
                             g_0_yy_0_xxzzzzz_1,   \
                             g_0_yy_0_xxzzzzzz_0,  \
                             g_0_yy_0_xxzzzzzz_1,  \
                             g_0_yy_0_xyyyyyy_1,   \
                             g_0_yy_0_xyyyyyyy_0,  \
                             g_0_yy_0_xyyyyyyy_1,  \
                             g_0_yy_0_xyyyyyyz_0,  \
                             g_0_yy_0_xyyyyyyz_1,  \
                             g_0_yy_0_xyyyyyz_1,   \
                             g_0_yy_0_xyyyyyzz_0,  \
                             g_0_yy_0_xyyyyyzz_1,  \
                             g_0_yy_0_xyyyyzz_1,   \
                             g_0_yy_0_xyyyyzzz_0,  \
                             g_0_yy_0_xyyyyzzz_1,  \
                             g_0_yy_0_xyyyzzz_1,   \
                             g_0_yy_0_xyyyzzzz_0,  \
                             g_0_yy_0_xyyyzzzz_1,  \
                             g_0_yy_0_xyyzzzz_1,   \
                             g_0_yy_0_xyyzzzzz_0,  \
                             g_0_yy_0_xyyzzzzz_1,  \
                             g_0_yy_0_xyzzzzz_1,   \
                             g_0_yy_0_xyzzzzzz_0,  \
                             g_0_yy_0_xyzzzzzz_1,  \
                             g_0_yy_0_xzzzzzz_1,   \
                             g_0_yy_0_xzzzzzzz_0,  \
                             g_0_yy_0_xzzzzzzz_1,  \
                             g_0_yy_0_yyyyyyy_1,   \
                             g_0_yy_0_yyyyyyyy_0,  \
                             g_0_yy_0_yyyyyyyy_1,  \
                             g_0_yy_0_yyyyyyyz_0,  \
                             g_0_yy_0_yyyyyyyz_1,  \
                             g_0_yy_0_yyyyyyz_1,   \
                             g_0_yy_0_yyyyyyzz_0,  \
                             g_0_yy_0_yyyyyyzz_1,  \
                             g_0_yy_0_yyyyyzz_1,   \
                             g_0_yy_0_yyyyyzzz_0,  \
                             g_0_yy_0_yyyyyzzz_1,  \
                             g_0_yy_0_yyyyzzz_1,   \
                             g_0_yy_0_yyyyzzzz_0,  \
                             g_0_yy_0_yyyyzzzz_1,  \
                             g_0_yy_0_yyyzzzz_1,   \
                             g_0_yy_0_yyyzzzzz_0,  \
                             g_0_yy_0_yyyzzzzz_1,  \
                             g_0_yy_0_yyzzzzz_1,   \
                             g_0_yy_0_yyzzzzzz_0,  \
                             g_0_yy_0_yyzzzzzz_1,  \
                             g_0_yy_0_yzzzzzz_1,   \
                             g_0_yy_0_yzzzzzzz_0,  \
                             g_0_yy_0_yzzzzzzz_1,  \
                             g_0_yy_0_zzzzzzz_1,   \
                             g_0_yy_0_zzzzzzzz_0,  \
                             g_0_yy_0_zzzzzzzz_1,  \
                             g_0_yyz_0_xxxxxxxx_0, \
                             g_0_yyz_0_xxxxxxxy_0, \
                             g_0_yyz_0_xxxxxxxz_0, \
                             g_0_yyz_0_xxxxxxyy_0, \
                             g_0_yyz_0_xxxxxxyz_0, \
                             g_0_yyz_0_xxxxxxzz_0, \
                             g_0_yyz_0_xxxxxyyy_0, \
                             g_0_yyz_0_xxxxxyyz_0, \
                             g_0_yyz_0_xxxxxyzz_0, \
                             g_0_yyz_0_xxxxxzzz_0, \
                             g_0_yyz_0_xxxxyyyy_0, \
                             g_0_yyz_0_xxxxyyyz_0, \
                             g_0_yyz_0_xxxxyyzz_0, \
                             g_0_yyz_0_xxxxyzzz_0, \
                             g_0_yyz_0_xxxxzzzz_0, \
                             g_0_yyz_0_xxxyyyyy_0, \
                             g_0_yyz_0_xxxyyyyz_0, \
                             g_0_yyz_0_xxxyyyzz_0, \
                             g_0_yyz_0_xxxyyzzz_0, \
                             g_0_yyz_0_xxxyzzzz_0, \
                             g_0_yyz_0_xxxzzzzz_0, \
                             g_0_yyz_0_xxyyyyyy_0, \
                             g_0_yyz_0_xxyyyyyz_0, \
                             g_0_yyz_0_xxyyyyzz_0, \
                             g_0_yyz_0_xxyyyzzz_0, \
                             g_0_yyz_0_xxyyzzzz_0, \
                             g_0_yyz_0_xxyzzzzz_0, \
                             g_0_yyz_0_xxzzzzzz_0, \
                             g_0_yyz_0_xyyyyyyy_0, \
                             g_0_yyz_0_xyyyyyyz_0, \
                             g_0_yyz_0_xyyyyyzz_0, \
                             g_0_yyz_0_xyyyyzzz_0, \
                             g_0_yyz_0_xyyyzzzz_0, \
                             g_0_yyz_0_xyyzzzzz_0, \
                             g_0_yyz_0_xyzzzzzz_0, \
                             g_0_yyz_0_xzzzzzzz_0, \
                             g_0_yyz_0_yyyyyyyy_0, \
                             g_0_yyz_0_yyyyyyyz_0, \
                             g_0_yyz_0_yyyyyyzz_0, \
                             g_0_yyz_0_yyyyyzzz_0, \
                             g_0_yyz_0_yyyyzzzz_0, \
                             g_0_yyz_0_yyyzzzzz_0, \
                             g_0_yyz_0_yyzzzzzz_0, \
                             g_0_yyz_0_yzzzzzzz_0, \
                             g_0_yyz_0_zzzzzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xxxxxxxx_0[i] = g_0_yy_0_xxxxxxxx_0[i] * pb_z + g_0_yy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxxxy_0[i] = g_0_yy_0_xxxxxxxy_0[i] * pb_z + g_0_yy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxxxz_0[i] = g_0_yy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxxz_0[i] * pb_z + g_0_yy_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxxyy_0[i] = g_0_yy_0_xxxxxxyy_0[i] * pb_z + g_0_yy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxxyz_0[i] = g_0_yy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxyz_0[i] * pb_z + g_0_yy_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxxzz_0[i] = 2.0 * g_0_yy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxzz_0[i] * pb_z + g_0_yy_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxyyy_0[i] = g_0_yy_0_xxxxxyyy_0[i] * pb_z + g_0_yy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxyyz_0[i] = g_0_yy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyyz_0[i] * pb_z + g_0_yy_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxyzz_0[i] = 2.0 * g_0_yy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyzz_0[i] * pb_z + g_0_yy_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxzzz_0[i] = 3.0 * g_0_yy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxzzz_0[i] * pb_z + g_0_yy_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyyyy_0[i] = g_0_yy_0_xxxxyyyy_0[i] * pb_z + g_0_yy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyyyz_0[i] = g_0_yy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyyz_0[i] * pb_z + g_0_yy_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyyzz_0[i] = 2.0 * g_0_yy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyzz_0[i] * pb_z + g_0_yy_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyzzz_0[i] = 3.0 * g_0_yy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyzzz_0[i] * pb_z + g_0_yy_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxzzzz_0[i] = 4.0 * g_0_yy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxzzzz_0[i] * pb_z + g_0_yy_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyyyy_0[i] = g_0_yy_0_xxxyyyyy_0[i] * pb_z + g_0_yy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyyyz_0[i] = g_0_yy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyyz_0[i] * pb_z + g_0_yy_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyyzz_0[i] = 2.0 * g_0_yy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyzz_0[i] * pb_z + g_0_yy_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyzzz_0[i] = 3.0 * g_0_yy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyzzz_0[i] * pb_z + g_0_yy_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyzzzz_0[i] = 4.0 * g_0_yy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyzzzz_0[i] * pb_z + g_0_yy_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxzzzzz_0[i] = 5.0 * g_0_yy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzzzzz_0[i] * pb_z + g_0_yy_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyyyy_0[i] = g_0_yy_0_xxyyyyyy_0[i] * pb_z + g_0_yy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyyyz_0[i] = g_0_yy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyyz_0[i] * pb_z + g_0_yy_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyyzz_0[i] = 2.0 * g_0_yy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyzz_0[i] * pb_z + g_0_yy_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyzzz_0[i] = 3.0 * g_0_yy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyzzz_0[i] * pb_z + g_0_yy_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyzzzz_0[i] = 4.0 * g_0_yy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzzzz_0[i] * pb_z + g_0_yy_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyzzzzz_0[i] = 5.0 * g_0_yy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzzzzz_0[i] * pb_z + g_0_yy_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxzzzzzz_0[i] = 6.0 * g_0_yy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzzzzz_0[i] * pb_z + g_0_yy_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyyyy_0[i] = g_0_yy_0_xyyyyyyy_0[i] * pb_z + g_0_yy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyyyz_0[i] = g_0_yy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyyz_0[i] * pb_z + g_0_yy_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyyzz_0[i] = 2.0 * g_0_yy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyzz_0[i] * pb_z + g_0_yy_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyzzz_0[i] = 3.0 * g_0_yy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyzzz_0[i] * pb_z + g_0_yy_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyzzzz_0[i] = 4.0 * g_0_yy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzzzz_0[i] * pb_z + g_0_yy_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyzzzzz_0[i] = 5.0 * g_0_yy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzzzz_0[i] * pb_z + g_0_yy_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyzzzzzz_0[i] = 6.0 * g_0_yy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzzzzz_0[i] * pb_z + g_0_yy_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xzzzzzzz_0[i] = 7.0 * g_0_yy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzzzzz_0[i] * pb_z + g_0_yy_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyyyy_0[i] = g_0_yy_0_yyyyyyyy_0[i] * pb_z + g_0_yy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyyyz_0[i] = g_0_yy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyyz_0[i] * pb_z + g_0_yy_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyyzz_0[i] = 2.0 * g_0_yy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyzz_0[i] * pb_z + g_0_yy_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyzzz_0[i] = 3.0 * g_0_yy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyzzz_0[i] * pb_z + g_0_yy_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyzzzz_0[i] = 4.0 * g_0_yy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyzzzz_0[i] * pb_z + g_0_yy_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyzzzzz_0[i] = 5.0 * g_0_yy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyzzzzz_0[i] * pb_z + g_0_yy_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyzzzzzz_0[i] = 6.0 * g_0_yy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzzzzzz_0[i] * pb_z + g_0_yy_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yzzzzzzz_0[i] = 7.0 * g_0_yy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yzzzzzzz_0[i] * pb_z + g_0_yy_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_zzzzzzzz_0[i] = 8.0 * g_0_yy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_zzzzzzzz_0[i] * pb_z + g_0_yy_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 360-405 components of targeted buffer : SFSL

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

#pragma omp simd aligned(g_0_yzz_0_xxxxxxxx_0,     \
                             g_0_yzz_0_xxxxxxxy_0, \
                             g_0_yzz_0_xxxxxxxz_0, \
                             g_0_yzz_0_xxxxxxyy_0, \
                             g_0_yzz_0_xxxxxxyz_0, \
                             g_0_yzz_0_xxxxxxzz_0, \
                             g_0_yzz_0_xxxxxyyy_0, \
                             g_0_yzz_0_xxxxxyyz_0, \
                             g_0_yzz_0_xxxxxyzz_0, \
                             g_0_yzz_0_xxxxxzzz_0, \
                             g_0_yzz_0_xxxxyyyy_0, \
                             g_0_yzz_0_xxxxyyyz_0, \
                             g_0_yzz_0_xxxxyyzz_0, \
                             g_0_yzz_0_xxxxyzzz_0, \
                             g_0_yzz_0_xxxxzzzz_0, \
                             g_0_yzz_0_xxxyyyyy_0, \
                             g_0_yzz_0_xxxyyyyz_0, \
                             g_0_yzz_0_xxxyyyzz_0, \
                             g_0_yzz_0_xxxyyzzz_0, \
                             g_0_yzz_0_xxxyzzzz_0, \
                             g_0_yzz_0_xxxzzzzz_0, \
                             g_0_yzz_0_xxyyyyyy_0, \
                             g_0_yzz_0_xxyyyyyz_0, \
                             g_0_yzz_0_xxyyyyzz_0, \
                             g_0_yzz_0_xxyyyzzz_0, \
                             g_0_yzz_0_xxyyzzzz_0, \
                             g_0_yzz_0_xxyzzzzz_0, \
                             g_0_yzz_0_xxzzzzzz_0, \
                             g_0_yzz_0_xyyyyyyy_0, \
                             g_0_yzz_0_xyyyyyyz_0, \
                             g_0_yzz_0_xyyyyyzz_0, \
                             g_0_yzz_0_xyyyyzzz_0, \
                             g_0_yzz_0_xyyyzzzz_0, \
                             g_0_yzz_0_xyyzzzzz_0, \
                             g_0_yzz_0_xyzzzzzz_0, \
                             g_0_yzz_0_xzzzzzzz_0, \
                             g_0_yzz_0_yyyyyyyy_0, \
                             g_0_yzz_0_yyyyyyyz_0, \
                             g_0_yzz_0_yyyyyyzz_0, \
                             g_0_yzz_0_yyyyyzzz_0, \
                             g_0_yzz_0_yyyyzzzz_0, \
                             g_0_yzz_0_yyyzzzzz_0, \
                             g_0_yzz_0_yyzzzzzz_0, \
                             g_0_yzz_0_yzzzzzzz_0, \
                             g_0_yzz_0_zzzzzzzz_0, \
                             g_0_zz_0_xxxxxxx_1,   \
                             g_0_zz_0_xxxxxxxx_0,  \
                             g_0_zz_0_xxxxxxxx_1,  \
                             g_0_zz_0_xxxxxxxy_0,  \
                             g_0_zz_0_xxxxxxxy_1,  \
                             g_0_zz_0_xxxxxxxz_0,  \
                             g_0_zz_0_xxxxxxxz_1,  \
                             g_0_zz_0_xxxxxxy_1,   \
                             g_0_zz_0_xxxxxxyy_0,  \
                             g_0_zz_0_xxxxxxyy_1,  \
                             g_0_zz_0_xxxxxxyz_0,  \
                             g_0_zz_0_xxxxxxyz_1,  \
                             g_0_zz_0_xxxxxxz_1,   \
                             g_0_zz_0_xxxxxxzz_0,  \
                             g_0_zz_0_xxxxxxzz_1,  \
                             g_0_zz_0_xxxxxyy_1,   \
                             g_0_zz_0_xxxxxyyy_0,  \
                             g_0_zz_0_xxxxxyyy_1,  \
                             g_0_zz_0_xxxxxyyz_0,  \
                             g_0_zz_0_xxxxxyyz_1,  \
                             g_0_zz_0_xxxxxyz_1,   \
                             g_0_zz_0_xxxxxyzz_0,  \
                             g_0_zz_0_xxxxxyzz_1,  \
                             g_0_zz_0_xxxxxzz_1,   \
                             g_0_zz_0_xxxxxzzz_0,  \
                             g_0_zz_0_xxxxxzzz_1,  \
                             g_0_zz_0_xxxxyyy_1,   \
                             g_0_zz_0_xxxxyyyy_0,  \
                             g_0_zz_0_xxxxyyyy_1,  \
                             g_0_zz_0_xxxxyyyz_0,  \
                             g_0_zz_0_xxxxyyyz_1,  \
                             g_0_zz_0_xxxxyyz_1,   \
                             g_0_zz_0_xxxxyyzz_0,  \
                             g_0_zz_0_xxxxyyzz_1,  \
                             g_0_zz_0_xxxxyzz_1,   \
                             g_0_zz_0_xxxxyzzz_0,  \
                             g_0_zz_0_xxxxyzzz_1,  \
                             g_0_zz_0_xxxxzzz_1,   \
                             g_0_zz_0_xxxxzzzz_0,  \
                             g_0_zz_0_xxxxzzzz_1,  \
                             g_0_zz_0_xxxyyyy_1,   \
                             g_0_zz_0_xxxyyyyy_0,  \
                             g_0_zz_0_xxxyyyyy_1,  \
                             g_0_zz_0_xxxyyyyz_0,  \
                             g_0_zz_0_xxxyyyyz_1,  \
                             g_0_zz_0_xxxyyyz_1,   \
                             g_0_zz_0_xxxyyyzz_0,  \
                             g_0_zz_0_xxxyyyzz_1,  \
                             g_0_zz_0_xxxyyzz_1,   \
                             g_0_zz_0_xxxyyzzz_0,  \
                             g_0_zz_0_xxxyyzzz_1,  \
                             g_0_zz_0_xxxyzzz_1,   \
                             g_0_zz_0_xxxyzzzz_0,  \
                             g_0_zz_0_xxxyzzzz_1,  \
                             g_0_zz_0_xxxzzzz_1,   \
                             g_0_zz_0_xxxzzzzz_0,  \
                             g_0_zz_0_xxxzzzzz_1,  \
                             g_0_zz_0_xxyyyyy_1,   \
                             g_0_zz_0_xxyyyyyy_0,  \
                             g_0_zz_0_xxyyyyyy_1,  \
                             g_0_zz_0_xxyyyyyz_0,  \
                             g_0_zz_0_xxyyyyyz_1,  \
                             g_0_zz_0_xxyyyyz_1,   \
                             g_0_zz_0_xxyyyyzz_0,  \
                             g_0_zz_0_xxyyyyzz_1,  \
                             g_0_zz_0_xxyyyzz_1,   \
                             g_0_zz_0_xxyyyzzz_0,  \
                             g_0_zz_0_xxyyyzzz_1,  \
                             g_0_zz_0_xxyyzzz_1,   \
                             g_0_zz_0_xxyyzzzz_0,  \
                             g_0_zz_0_xxyyzzzz_1,  \
                             g_0_zz_0_xxyzzzz_1,   \
                             g_0_zz_0_xxyzzzzz_0,  \
                             g_0_zz_0_xxyzzzzz_1,  \
                             g_0_zz_0_xxzzzzz_1,   \
                             g_0_zz_0_xxzzzzzz_0,  \
                             g_0_zz_0_xxzzzzzz_1,  \
                             g_0_zz_0_xyyyyyy_1,   \
                             g_0_zz_0_xyyyyyyy_0,  \
                             g_0_zz_0_xyyyyyyy_1,  \
                             g_0_zz_0_xyyyyyyz_0,  \
                             g_0_zz_0_xyyyyyyz_1,  \
                             g_0_zz_0_xyyyyyz_1,   \
                             g_0_zz_0_xyyyyyzz_0,  \
                             g_0_zz_0_xyyyyyzz_1,  \
                             g_0_zz_0_xyyyyzz_1,   \
                             g_0_zz_0_xyyyyzzz_0,  \
                             g_0_zz_0_xyyyyzzz_1,  \
                             g_0_zz_0_xyyyzzz_1,   \
                             g_0_zz_0_xyyyzzzz_0,  \
                             g_0_zz_0_xyyyzzzz_1,  \
                             g_0_zz_0_xyyzzzz_1,   \
                             g_0_zz_0_xyyzzzzz_0,  \
                             g_0_zz_0_xyyzzzzz_1,  \
                             g_0_zz_0_xyzzzzz_1,   \
                             g_0_zz_0_xyzzzzzz_0,  \
                             g_0_zz_0_xyzzzzzz_1,  \
                             g_0_zz_0_xzzzzzz_1,   \
                             g_0_zz_0_xzzzzzzz_0,  \
                             g_0_zz_0_xzzzzzzz_1,  \
                             g_0_zz_0_yyyyyyy_1,   \
                             g_0_zz_0_yyyyyyyy_0,  \
                             g_0_zz_0_yyyyyyyy_1,  \
                             g_0_zz_0_yyyyyyyz_0,  \
                             g_0_zz_0_yyyyyyyz_1,  \
                             g_0_zz_0_yyyyyyz_1,   \
                             g_0_zz_0_yyyyyyzz_0,  \
                             g_0_zz_0_yyyyyyzz_1,  \
                             g_0_zz_0_yyyyyzz_1,   \
                             g_0_zz_0_yyyyyzzz_0,  \
                             g_0_zz_0_yyyyyzzz_1,  \
                             g_0_zz_0_yyyyzzz_1,   \
                             g_0_zz_0_yyyyzzzz_0,  \
                             g_0_zz_0_yyyyzzzz_1,  \
                             g_0_zz_0_yyyzzzz_1,   \
                             g_0_zz_0_yyyzzzzz_0,  \
                             g_0_zz_0_yyyzzzzz_1,  \
                             g_0_zz_0_yyzzzzz_1,   \
                             g_0_zz_0_yyzzzzzz_0,  \
                             g_0_zz_0_yyzzzzzz_1,  \
                             g_0_zz_0_yzzzzzz_1,   \
                             g_0_zz_0_yzzzzzzz_0,  \
                             g_0_zz_0_yzzzzzzz_1,  \
                             g_0_zz_0_zzzzzzz_1,   \
                             g_0_zz_0_zzzzzzzz_0,  \
                             g_0_zz_0_zzzzzzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xxxxxxxx_0[i] = g_0_zz_0_xxxxxxxx_0[i] * pb_y + g_0_zz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxxxy_0[i] = g_0_zz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxxy_0[i] * pb_y + g_0_zz_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxxxz_0[i] = g_0_zz_0_xxxxxxxz_0[i] * pb_y + g_0_zz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxxyy_0[i] = 2.0 * g_0_zz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxyy_0[i] * pb_y + g_0_zz_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxxyz_0[i] = g_0_zz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxyz_0[i] * pb_y + g_0_zz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxxzz_0[i] = g_0_zz_0_xxxxxxzz_0[i] * pb_y + g_0_zz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxyyy_0[i] = 3.0 * g_0_zz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyyy_0[i] * pb_y + g_0_zz_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxyyz_0[i] = 2.0 * g_0_zz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyyz_0[i] * pb_y + g_0_zz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxyzz_0[i] = g_0_zz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyzz_0[i] * pb_y + g_0_zz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxzzz_0[i] = g_0_zz_0_xxxxxzzz_0[i] * pb_y + g_0_zz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyyyy_0[i] = 4.0 * g_0_zz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyyy_0[i] * pb_y + g_0_zz_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyyyz_0[i] = 3.0 * g_0_zz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyyz_0[i] * pb_y + g_0_zz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyyzz_0[i] = 2.0 * g_0_zz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyzz_0[i] * pb_y + g_0_zz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyzzz_0[i] = g_0_zz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyzzz_0[i] * pb_y + g_0_zz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxzzzz_0[i] = g_0_zz_0_xxxxzzzz_0[i] * pb_y + g_0_zz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyyyy_0[i] = 5.0 * g_0_zz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyyy_0[i] * pb_y + g_0_zz_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyyyz_0[i] = 4.0 * g_0_zz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyyz_0[i] * pb_y + g_0_zz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyyzz_0[i] = 3.0 * g_0_zz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyzz_0[i] * pb_y + g_0_zz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyzzz_0[i] = 2.0 * g_0_zz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyzzz_0[i] * pb_y + g_0_zz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyzzzz_0[i] = g_0_zz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzzzz_0[i] * pb_y + g_0_zz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxzzzzz_0[i] = g_0_zz_0_xxxzzzzz_0[i] * pb_y + g_0_zz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyyyy_0[i] = 6.0 * g_0_zz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyyy_0[i] * pb_y + g_0_zz_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyyyz_0[i] = 5.0 * g_0_zz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyyz_0[i] * pb_y + g_0_zz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyyzz_0[i] = 4.0 * g_0_zz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyzz_0[i] * pb_y + g_0_zz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyzzz_0[i] = 3.0 * g_0_zz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyzzz_0[i] * pb_y + g_0_zz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyzzzz_0[i] = 2.0 * g_0_zz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzzzz_0[i] * pb_y + g_0_zz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyzzzzz_0[i] = g_0_zz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzzzz_0[i] * pb_y + g_0_zz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxzzzzzz_0[i] = g_0_zz_0_xxzzzzzz_0[i] * pb_y + g_0_zz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyyyy_0[i] = 7.0 * g_0_zz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyyy_0[i] * pb_y + g_0_zz_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyyyz_0[i] = 6.0 * g_0_zz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyyz_0[i] * pb_y + g_0_zz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyyzz_0[i] = 5.0 * g_0_zz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyzz_0[i] * pb_y + g_0_zz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyzzz_0[i] = 4.0 * g_0_zz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyzzz_0[i] * pb_y + g_0_zz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyzzzz_0[i] = 3.0 * g_0_zz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzzzz_0[i] * pb_y + g_0_zz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyzzzzz_0[i] = 2.0 * g_0_zz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzzzz_0[i] * pb_y + g_0_zz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyzzzzzz_0[i] = g_0_zz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzzzz_0[i] * pb_y + g_0_zz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xzzzzzzz_0[i] = g_0_zz_0_xzzzzzzz_0[i] * pb_y + g_0_zz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyyyy_0[i] = 8.0 * g_0_zz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyyyy_0[i] * pb_y + g_0_zz_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyyyz_0[i] = 7.0 * g_0_zz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyyyz_0[i] * pb_y + g_0_zz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyyzz_0[i] = 6.0 * g_0_zz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyyzz_0[i] * pb_y + g_0_zz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyzzz_0[i] = 5.0 * g_0_zz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyzzz_0[i] * pb_y + g_0_zz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyzzzz_0[i] = 4.0 * g_0_zz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyzzzz_0[i] * pb_y + g_0_zz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyzzzzz_0[i] = 3.0 * g_0_zz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyzzzzz_0[i] * pb_y + g_0_zz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyzzzzzz_0[i] = 2.0 * g_0_zz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzzzzzz_0[i] * pb_y + g_0_zz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yzzzzzzz_0[i] = g_0_zz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzzzzzz_0[i] * pb_y + g_0_zz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_zzzzzzzz_0[i] = g_0_zz_0_zzzzzzzz_0[i] * pb_y + g_0_zz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 405-450 components of targeted buffer : SFSL

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

#pragma omp simd aligned(g_0_z_0_xxxxxxxx_0,       \
                             g_0_z_0_xxxxxxxx_1,   \
                             g_0_z_0_xxxxxxxy_0,   \
                             g_0_z_0_xxxxxxxy_1,   \
                             g_0_z_0_xxxxxxxz_0,   \
                             g_0_z_0_xxxxxxxz_1,   \
                             g_0_z_0_xxxxxxyy_0,   \
                             g_0_z_0_xxxxxxyy_1,   \
                             g_0_z_0_xxxxxxyz_0,   \
                             g_0_z_0_xxxxxxyz_1,   \
                             g_0_z_0_xxxxxxzz_0,   \
                             g_0_z_0_xxxxxxzz_1,   \
                             g_0_z_0_xxxxxyyy_0,   \
                             g_0_z_0_xxxxxyyy_1,   \
                             g_0_z_0_xxxxxyyz_0,   \
                             g_0_z_0_xxxxxyyz_1,   \
                             g_0_z_0_xxxxxyzz_0,   \
                             g_0_z_0_xxxxxyzz_1,   \
                             g_0_z_0_xxxxxzzz_0,   \
                             g_0_z_0_xxxxxzzz_1,   \
                             g_0_z_0_xxxxyyyy_0,   \
                             g_0_z_0_xxxxyyyy_1,   \
                             g_0_z_0_xxxxyyyz_0,   \
                             g_0_z_0_xxxxyyyz_1,   \
                             g_0_z_0_xxxxyyzz_0,   \
                             g_0_z_0_xxxxyyzz_1,   \
                             g_0_z_0_xxxxyzzz_0,   \
                             g_0_z_0_xxxxyzzz_1,   \
                             g_0_z_0_xxxxzzzz_0,   \
                             g_0_z_0_xxxxzzzz_1,   \
                             g_0_z_0_xxxyyyyy_0,   \
                             g_0_z_0_xxxyyyyy_1,   \
                             g_0_z_0_xxxyyyyz_0,   \
                             g_0_z_0_xxxyyyyz_1,   \
                             g_0_z_0_xxxyyyzz_0,   \
                             g_0_z_0_xxxyyyzz_1,   \
                             g_0_z_0_xxxyyzzz_0,   \
                             g_0_z_0_xxxyyzzz_1,   \
                             g_0_z_0_xxxyzzzz_0,   \
                             g_0_z_0_xxxyzzzz_1,   \
                             g_0_z_0_xxxzzzzz_0,   \
                             g_0_z_0_xxxzzzzz_1,   \
                             g_0_z_0_xxyyyyyy_0,   \
                             g_0_z_0_xxyyyyyy_1,   \
                             g_0_z_0_xxyyyyyz_0,   \
                             g_0_z_0_xxyyyyyz_1,   \
                             g_0_z_0_xxyyyyzz_0,   \
                             g_0_z_0_xxyyyyzz_1,   \
                             g_0_z_0_xxyyyzzz_0,   \
                             g_0_z_0_xxyyyzzz_1,   \
                             g_0_z_0_xxyyzzzz_0,   \
                             g_0_z_0_xxyyzzzz_1,   \
                             g_0_z_0_xxyzzzzz_0,   \
                             g_0_z_0_xxyzzzzz_1,   \
                             g_0_z_0_xxzzzzzz_0,   \
                             g_0_z_0_xxzzzzzz_1,   \
                             g_0_z_0_xyyyyyyy_0,   \
                             g_0_z_0_xyyyyyyy_1,   \
                             g_0_z_0_xyyyyyyz_0,   \
                             g_0_z_0_xyyyyyyz_1,   \
                             g_0_z_0_xyyyyyzz_0,   \
                             g_0_z_0_xyyyyyzz_1,   \
                             g_0_z_0_xyyyyzzz_0,   \
                             g_0_z_0_xyyyyzzz_1,   \
                             g_0_z_0_xyyyzzzz_0,   \
                             g_0_z_0_xyyyzzzz_1,   \
                             g_0_z_0_xyyzzzzz_0,   \
                             g_0_z_0_xyyzzzzz_1,   \
                             g_0_z_0_xyzzzzzz_0,   \
                             g_0_z_0_xyzzzzzz_1,   \
                             g_0_z_0_xzzzzzzz_0,   \
                             g_0_z_0_xzzzzzzz_1,   \
                             g_0_z_0_yyyyyyyy_0,   \
                             g_0_z_0_yyyyyyyy_1,   \
                             g_0_z_0_yyyyyyyz_0,   \
                             g_0_z_0_yyyyyyyz_1,   \
                             g_0_z_0_yyyyyyzz_0,   \
                             g_0_z_0_yyyyyyzz_1,   \
                             g_0_z_0_yyyyyzzz_0,   \
                             g_0_z_0_yyyyyzzz_1,   \
                             g_0_z_0_yyyyzzzz_0,   \
                             g_0_z_0_yyyyzzzz_1,   \
                             g_0_z_0_yyyzzzzz_0,   \
                             g_0_z_0_yyyzzzzz_1,   \
                             g_0_z_0_yyzzzzzz_0,   \
                             g_0_z_0_yyzzzzzz_1,   \
                             g_0_z_0_yzzzzzzz_0,   \
                             g_0_z_0_yzzzzzzz_1,   \
                             g_0_z_0_zzzzzzzz_0,   \
                             g_0_z_0_zzzzzzzz_1,   \
                             g_0_zz_0_xxxxxxx_1,   \
                             g_0_zz_0_xxxxxxxx_0,  \
                             g_0_zz_0_xxxxxxxx_1,  \
                             g_0_zz_0_xxxxxxxy_0,  \
                             g_0_zz_0_xxxxxxxy_1,  \
                             g_0_zz_0_xxxxxxxz_0,  \
                             g_0_zz_0_xxxxxxxz_1,  \
                             g_0_zz_0_xxxxxxy_1,   \
                             g_0_zz_0_xxxxxxyy_0,  \
                             g_0_zz_0_xxxxxxyy_1,  \
                             g_0_zz_0_xxxxxxyz_0,  \
                             g_0_zz_0_xxxxxxyz_1,  \
                             g_0_zz_0_xxxxxxz_1,   \
                             g_0_zz_0_xxxxxxzz_0,  \
                             g_0_zz_0_xxxxxxzz_1,  \
                             g_0_zz_0_xxxxxyy_1,   \
                             g_0_zz_0_xxxxxyyy_0,  \
                             g_0_zz_0_xxxxxyyy_1,  \
                             g_0_zz_0_xxxxxyyz_0,  \
                             g_0_zz_0_xxxxxyyz_1,  \
                             g_0_zz_0_xxxxxyz_1,   \
                             g_0_zz_0_xxxxxyzz_0,  \
                             g_0_zz_0_xxxxxyzz_1,  \
                             g_0_zz_0_xxxxxzz_1,   \
                             g_0_zz_0_xxxxxzzz_0,  \
                             g_0_zz_0_xxxxxzzz_1,  \
                             g_0_zz_0_xxxxyyy_1,   \
                             g_0_zz_0_xxxxyyyy_0,  \
                             g_0_zz_0_xxxxyyyy_1,  \
                             g_0_zz_0_xxxxyyyz_0,  \
                             g_0_zz_0_xxxxyyyz_1,  \
                             g_0_zz_0_xxxxyyz_1,   \
                             g_0_zz_0_xxxxyyzz_0,  \
                             g_0_zz_0_xxxxyyzz_1,  \
                             g_0_zz_0_xxxxyzz_1,   \
                             g_0_zz_0_xxxxyzzz_0,  \
                             g_0_zz_0_xxxxyzzz_1,  \
                             g_0_zz_0_xxxxzzz_1,   \
                             g_0_zz_0_xxxxzzzz_0,  \
                             g_0_zz_0_xxxxzzzz_1,  \
                             g_0_zz_0_xxxyyyy_1,   \
                             g_0_zz_0_xxxyyyyy_0,  \
                             g_0_zz_0_xxxyyyyy_1,  \
                             g_0_zz_0_xxxyyyyz_0,  \
                             g_0_zz_0_xxxyyyyz_1,  \
                             g_0_zz_0_xxxyyyz_1,   \
                             g_0_zz_0_xxxyyyzz_0,  \
                             g_0_zz_0_xxxyyyzz_1,  \
                             g_0_zz_0_xxxyyzz_1,   \
                             g_0_zz_0_xxxyyzzz_0,  \
                             g_0_zz_0_xxxyyzzz_1,  \
                             g_0_zz_0_xxxyzzz_1,   \
                             g_0_zz_0_xxxyzzzz_0,  \
                             g_0_zz_0_xxxyzzzz_1,  \
                             g_0_zz_0_xxxzzzz_1,   \
                             g_0_zz_0_xxxzzzzz_0,  \
                             g_0_zz_0_xxxzzzzz_1,  \
                             g_0_zz_0_xxyyyyy_1,   \
                             g_0_zz_0_xxyyyyyy_0,  \
                             g_0_zz_0_xxyyyyyy_1,  \
                             g_0_zz_0_xxyyyyyz_0,  \
                             g_0_zz_0_xxyyyyyz_1,  \
                             g_0_zz_0_xxyyyyz_1,   \
                             g_0_zz_0_xxyyyyzz_0,  \
                             g_0_zz_0_xxyyyyzz_1,  \
                             g_0_zz_0_xxyyyzz_1,   \
                             g_0_zz_0_xxyyyzzz_0,  \
                             g_0_zz_0_xxyyyzzz_1,  \
                             g_0_zz_0_xxyyzzz_1,   \
                             g_0_zz_0_xxyyzzzz_0,  \
                             g_0_zz_0_xxyyzzzz_1,  \
                             g_0_zz_0_xxyzzzz_1,   \
                             g_0_zz_0_xxyzzzzz_0,  \
                             g_0_zz_0_xxyzzzzz_1,  \
                             g_0_zz_0_xxzzzzz_1,   \
                             g_0_zz_0_xxzzzzzz_0,  \
                             g_0_zz_0_xxzzzzzz_1,  \
                             g_0_zz_0_xyyyyyy_1,   \
                             g_0_zz_0_xyyyyyyy_0,  \
                             g_0_zz_0_xyyyyyyy_1,  \
                             g_0_zz_0_xyyyyyyz_0,  \
                             g_0_zz_0_xyyyyyyz_1,  \
                             g_0_zz_0_xyyyyyz_1,   \
                             g_0_zz_0_xyyyyyzz_0,  \
                             g_0_zz_0_xyyyyyzz_1,  \
                             g_0_zz_0_xyyyyzz_1,   \
                             g_0_zz_0_xyyyyzzz_0,  \
                             g_0_zz_0_xyyyyzzz_1,  \
                             g_0_zz_0_xyyyzzz_1,   \
                             g_0_zz_0_xyyyzzzz_0,  \
                             g_0_zz_0_xyyyzzzz_1,  \
                             g_0_zz_0_xyyzzzz_1,   \
                             g_0_zz_0_xyyzzzzz_0,  \
                             g_0_zz_0_xyyzzzzz_1,  \
                             g_0_zz_0_xyzzzzz_1,   \
                             g_0_zz_0_xyzzzzzz_0,  \
                             g_0_zz_0_xyzzzzzz_1,  \
                             g_0_zz_0_xzzzzzz_1,   \
                             g_0_zz_0_xzzzzzzz_0,  \
                             g_0_zz_0_xzzzzzzz_1,  \
                             g_0_zz_0_yyyyyyy_1,   \
                             g_0_zz_0_yyyyyyyy_0,  \
                             g_0_zz_0_yyyyyyyy_1,  \
                             g_0_zz_0_yyyyyyyz_0,  \
                             g_0_zz_0_yyyyyyyz_1,  \
                             g_0_zz_0_yyyyyyz_1,   \
                             g_0_zz_0_yyyyyyzz_0,  \
                             g_0_zz_0_yyyyyyzz_1,  \
                             g_0_zz_0_yyyyyzz_1,   \
                             g_0_zz_0_yyyyyzzz_0,  \
                             g_0_zz_0_yyyyyzzz_1,  \
                             g_0_zz_0_yyyyzzz_1,   \
                             g_0_zz_0_yyyyzzzz_0,  \
                             g_0_zz_0_yyyyzzzz_1,  \
                             g_0_zz_0_yyyzzzz_1,   \
                             g_0_zz_0_yyyzzzzz_0,  \
                             g_0_zz_0_yyyzzzzz_1,  \
                             g_0_zz_0_yyzzzzz_1,   \
                             g_0_zz_0_yyzzzzzz_0,  \
                             g_0_zz_0_yyzzzzzz_1,  \
                             g_0_zz_0_yzzzzzz_1,   \
                             g_0_zz_0_yzzzzzzz_0,  \
                             g_0_zz_0_yzzzzzzz_1,  \
                             g_0_zz_0_zzzzzzz_1,   \
                             g_0_zz_0_zzzzzzzz_0,  \
                             g_0_zz_0_zzzzzzzz_1,  \
                             g_0_zzz_0_xxxxxxxx_0, \
                             g_0_zzz_0_xxxxxxxy_0, \
                             g_0_zzz_0_xxxxxxxz_0, \
                             g_0_zzz_0_xxxxxxyy_0, \
                             g_0_zzz_0_xxxxxxyz_0, \
                             g_0_zzz_0_xxxxxxzz_0, \
                             g_0_zzz_0_xxxxxyyy_0, \
                             g_0_zzz_0_xxxxxyyz_0, \
                             g_0_zzz_0_xxxxxyzz_0, \
                             g_0_zzz_0_xxxxxzzz_0, \
                             g_0_zzz_0_xxxxyyyy_0, \
                             g_0_zzz_0_xxxxyyyz_0, \
                             g_0_zzz_0_xxxxyyzz_0, \
                             g_0_zzz_0_xxxxyzzz_0, \
                             g_0_zzz_0_xxxxzzzz_0, \
                             g_0_zzz_0_xxxyyyyy_0, \
                             g_0_zzz_0_xxxyyyyz_0, \
                             g_0_zzz_0_xxxyyyzz_0, \
                             g_0_zzz_0_xxxyyzzz_0, \
                             g_0_zzz_0_xxxyzzzz_0, \
                             g_0_zzz_0_xxxzzzzz_0, \
                             g_0_zzz_0_xxyyyyyy_0, \
                             g_0_zzz_0_xxyyyyyz_0, \
                             g_0_zzz_0_xxyyyyzz_0, \
                             g_0_zzz_0_xxyyyzzz_0, \
                             g_0_zzz_0_xxyyzzzz_0, \
                             g_0_zzz_0_xxyzzzzz_0, \
                             g_0_zzz_0_xxzzzzzz_0, \
                             g_0_zzz_0_xyyyyyyy_0, \
                             g_0_zzz_0_xyyyyyyz_0, \
                             g_0_zzz_0_xyyyyyzz_0, \
                             g_0_zzz_0_xyyyyzzz_0, \
                             g_0_zzz_0_xyyyzzzz_0, \
                             g_0_zzz_0_xyyzzzzz_0, \
                             g_0_zzz_0_xyzzzzzz_0, \
                             g_0_zzz_0_xzzzzzzz_0, \
                             g_0_zzz_0_yyyyyyyy_0, \
                             g_0_zzz_0_yyyyyyyz_0, \
                             g_0_zzz_0_yyyyyyzz_0, \
                             g_0_zzz_0_yyyyyzzz_0, \
                             g_0_zzz_0_yyyyzzzz_0, \
                             g_0_zzz_0_yyyzzzzz_0, \
                             g_0_zzz_0_yyzzzzzz_0, \
                             g_0_zzz_0_yzzzzzzz_0, \
                             g_0_zzz_0_zzzzzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xxxxxxxx_0[i] = 2.0 * g_0_z_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_zz_0_xxxxxxxx_0[i] * pb_z +
                                  g_0_zz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxxxy_0[i] = 2.0 * g_0_z_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_zz_0_xxxxxxxy_0[i] * pb_z +
                                  g_0_zz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxxxz_0[i] = 2.0 * g_0_z_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_zz_0_xxxxxxx_1[i] * fi_abcd_0 +
                                  g_0_zz_0_xxxxxxxz_0[i] * pb_z + g_0_zz_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxxyy_0[i] = 2.0 * g_0_z_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_zz_0_xxxxxxyy_0[i] * pb_z +
                                  g_0_zz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxxyz_0[i] = 2.0 * g_0_z_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_zz_0_xxxxxxy_1[i] * fi_abcd_0 +
                                  g_0_zz_0_xxxxxxyz_0[i] * pb_z + g_0_zz_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxxzz_0[i] = 2.0 * g_0_z_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxzz_0[i] * pb_z + g_0_zz_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxyyy_0[i] = 2.0 * g_0_z_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_zz_0_xxxxxyyy_0[i] * pb_z +
                                  g_0_zz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxyyz_0[i] = 2.0 * g_0_z_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxyyz_1[i] * fti_ab_0 + g_0_zz_0_xxxxxyy_1[i] * fi_abcd_0 +
                                  g_0_zz_0_xxxxxyyz_0[i] * pb_z + g_0_zz_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxyzz_0[i] = 2.0 * g_0_z_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyzz_0[i] * pb_z + g_0_zz_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxzzz_0[i] = 2.0 * g_0_z_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxzzz_0[i] * pb_z + g_0_zz_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyyyy_0[i] = 2.0 * g_0_z_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_zz_0_xxxxyyyy_0[i] * pb_z +
                                  g_0_zz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyyyz_0[i] = 2.0 * g_0_z_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyyyz_1[i] * fti_ab_0 + g_0_zz_0_xxxxyyy_1[i] * fi_abcd_0 +
                                  g_0_zz_0_xxxxyyyz_0[i] * pb_z + g_0_zz_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyyzz_0[i] = 2.0 * g_0_z_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyzz_0[i] * pb_z + g_0_zz_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyzzz_0[i] = 2.0 * g_0_z_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyzzz_0[i] * pb_z + g_0_zz_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxzzzz_0[i] = 2.0 * g_0_z_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxzzzz_0[i] * pb_z + g_0_zz_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyyyy_0[i] = 2.0 * g_0_z_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_zz_0_xxxyyyyy_0[i] * pb_z +
                                  g_0_zz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyyyz_0[i] = 2.0 * g_0_z_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyyyz_1[i] * fti_ab_0 + g_0_zz_0_xxxyyyy_1[i] * fi_abcd_0 +
                                  g_0_zz_0_xxxyyyyz_0[i] * pb_z + g_0_zz_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyyzz_0[i] = 2.0 * g_0_z_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyzz_0[i] * pb_z + g_0_zz_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyzzz_0[i] = 2.0 * g_0_z_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyzzz_0[i] * pb_z + g_0_zz_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyzzzz_0[i] = 2.0 * g_0_z_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzzzz_0[i] * pb_z + g_0_zz_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxzzzzz_0[i] = 2.0 * g_0_z_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxzzzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_zz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxzzzzz_0[i] * pb_z + g_0_zz_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyyyy_0[i] = 2.0 * g_0_z_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_zz_0_xxyyyyyy_0[i] * pb_z +
                                  g_0_zz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyyyz_0[i] = 2.0 * g_0_z_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyyyz_1[i] * fti_ab_0 + g_0_zz_0_xxyyyyy_1[i] * fi_abcd_0 +
                                  g_0_zz_0_xxyyyyyz_0[i] * pb_z + g_0_zz_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyyzz_0[i] = 2.0 * g_0_z_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyzz_0[i] * pb_z + g_0_zz_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyzzz_0[i] = 2.0 * g_0_z_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyzzz_0[i] * pb_z + g_0_zz_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyzzzz_0[i] = 2.0 * g_0_z_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzzzz_0[i] * pb_z + g_0_zz_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyzzzzz_0[i] = 2.0 * g_0_z_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyzzzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_zz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzzzz_0[i] * pb_z + g_0_zz_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxzzzzzz_0[i] = 2.0 * g_0_z_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxzzzzzz_1[i] * fti_ab_0 +
                                  6.0 * g_0_zz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzzzzzz_0[i] * pb_z + g_0_zz_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyyyy_0[i] = 2.0 * g_0_z_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_zz_0_xyyyyyyy_0[i] * pb_z +
                                  g_0_zz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyyyz_0[i] = 2.0 * g_0_z_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_zz_0_xyyyyyy_1[i] * fi_abcd_0 +
                                  g_0_zz_0_xyyyyyyz_0[i] * pb_z + g_0_zz_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyyzz_0[i] = 2.0 * g_0_z_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyzz_0[i] * pb_z + g_0_zz_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyzzz_0[i] = 2.0 * g_0_z_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyzzz_0[i] * pb_z + g_0_zz_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyzzzz_0[i] = 2.0 * g_0_z_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzzzz_0[i] * pb_z + g_0_zz_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyzzzzz_0[i] = 2.0 * g_0_z_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyzzzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_zz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzzzz_0[i] * pb_z + g_0_zz_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyzzzzzz_0[i] = 2.0 * g_0_z_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyzzzzzz_1[i] * fti_ab_0 +
                                  6.0 * g_0_zz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzzzz_0[i] * pb_z + g_0_zz_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xzzzzzzz_0[i] = 2.0 * g_0_z_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xzzzzzzz_1[i] * fti_ab_0 +
                                  7.0 * g_0_zz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzzzzzz_0[i] * pb_z + g_0_zz_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyyyy_0[i] = 2.0 * g_0_z_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_zz_0_yyyyyyyy_0[i] * pb_z +
                                  g_0_zz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyyyz_0[i] = 2.0 * g_0_z_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_zz_0_yyyyyyy_1[i] * fi_abcd_0 +
                                  g_0_zz_0_yyyyyyyz_0[i] * pb_z + g_0_zz_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyyzz_0[i] = 2.0 * g_0_z_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyyzz_0[i] * pb_z + g_0_zz_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyzzz_0[i] = 2.0 * g_0_z_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyzzz_0[i] * pb_z + g_0_zz_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyzzzz_0[i] = 2.0 * g_0_z_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyzzzz_0[i] * pb_z + g_0_zz_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyzzzzz_0[i] = 2.0 * g_0_z_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyzzzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_zz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyzzzzz_0[i] * pb_z + g_0_zz_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyzzzzzz_0[i] = 2.0 * g_0_z_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyzzzzzz_1[i] * fti_ab_0 +
                                  6.0 * g_0_zz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzzzzzz_0[i] * pb_z + g_0_zz_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yzzzzzzz_0[i] = 2.0 * g_0_z_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yzzzzzzz_1[i] * fti_ab_0 +
                                  7.0 * g_0_zz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzzzzzz_0[i] * pb_z + g_0_zz_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_zzzzzzzz_0[i] = 2.0 * g_0_z_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zzzzzzzz_1[i] * fti_ab_0 +
                                  8.0 * g_0_zz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_zzzzzzzz_0[i] * pb_z + g_0_zz_0_zzzzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
