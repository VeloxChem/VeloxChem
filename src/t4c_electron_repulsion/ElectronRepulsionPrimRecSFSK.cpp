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

#include "ElectronRepulsionPrimRecSFSK.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sfsk(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sfsk,
                                  size_t                idx_eri_0_spsk,
                                  size_t                idx_eri_1_spsk,
                                  size_t                idx_eri_1_sdsi,
                                  size_t                idx_eri_0_sdsk,
                                  size_t                idx_eri_1_sdsk,
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

    /// Set up components of auxilary buffer : SPSK

    auto g_0_x_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_spsk);

    auto g_0_x_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_spsk + 1);

    auto g_0_x_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_spsk + 2);

    auto g_0_x_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_spsk + 3);

    auto g_0_x_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_spsk + 4);

    auto g_0_x_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_spsk + 5);

    auto g_0_x_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_spsk + 6);

    auto g_0_x_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_spsk + 7);

    auto g_0_x_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_spsk + 8);

    auto g_0_x_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_spsk + 9);

    auto g_0_x_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_spsk + 10);

    auto g_0_x_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_spsk + 11);

    auto g_0_x_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_spsk + 12);

    auto g_0_x_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_spsk + 13);

    auto g_0_x_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_spsk + 14);

    auto g_0_x_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 15);

    auto g_0_x_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 16);

    auto g_0_x_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 17);

    auto g_0_x_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 18);

    auto g_0_x_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 19);

    auto g_0_x_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 20);

    auto g_0_x_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 21);

    auto g_0_x_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 22);

    auto g_0_x_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 23);

    auto g_0_x_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 24);

    auto g_0_x_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 25);

    auto g_0_x_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 26);

    auto g_0_x_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 27);

    auto g_0_x_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 28);

    auto g_0_x_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 29);

    auto g_0_x_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 30);

    auto g_0_x_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 31);

    auto g_0_x_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 32);

    auto g_0_x_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 33);

    auto g_0_x_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 34);

    auto g_0_x_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 35);

    auto g_0_y_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_spsk + 36);

    auto g_0_y_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_spsk + 37);

    auto g_0_y_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_spsk + 38);

    auto g_0_y_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_spsk + 39);

    auto g_0_y_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_spsk + 40);

    auto g_0_y_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_spsk + 41);

    auto g_0_y_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_spsk + 42);

    auto g_0_y_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_spsk + 43);

    auto g_0_y_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_spsk + 44);

    auto g_0_y_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_spsk + 45);

    auto g_0_y_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_spsk + 46);

    auto g_0_y_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_spsk + 47);

    auto g_0_y_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_spsk + 48);

    auto g_0_y_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_spsk + 49);

    auto g_0_y_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_spsk + 50);

    auto g_0_y_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 51);

    auto g_0_y_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 52);

    auto g_0_y_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 53);

    auto g_0_y_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 54);

    auto g_0_y_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 55);

    auto g_0_y_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 56);

    auto g_0_y_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 57);

    auto g_0_y_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 58);

    auto g_0_y_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 59);

    auto g_0_y_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 60);

    auto g_0_y_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 61);

    auto g_0_y_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 62);

    auto g_0_y_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 63);

    auto g_0_y_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 64);

    auto g_0_y_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 65);

    auto g_0_y_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 66);

    auto g_0_y_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 67);

    auto g_0_y_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 68);

    auto g_0_y_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 69);

    auto g_0_y_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 70);

    auto g_0_y_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 71);

    auto g_0_z_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_spsk + 72);

    auto g_0_z_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_spsk + 73);

    auto g_0_z_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_spsk + 74);

    auto g_0_z_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_spsk + 75);

    auto g_0_z_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_spsk + 76);

    auto g_0_z_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_spsk + 77);

    auto g_0_z_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_spsk + 78);

    auto g_0_z_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_spsk + 79);

    auto g_0_z_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_spsk + 80);

    auto g_0_z_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_spsk + 81);

    auto g_0_z_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_spsk + 82);

    auto g_0_z_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_spsk + 83);

    auto g_0_z_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_spsk + 84);

    auto g_0_z_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_spsk + 85);

    auto g_0_z_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_spsk + 86);

    auto g_0_z_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 87);

    auto g_0_z_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 88);

    auto g_0_z_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 89);

    auto g_0_z_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 90);

    auto g_0_z_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 91);

    auto g_0_z_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 92);

    auto g_0_z_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 93);

    auto g_0_z_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 94);

    auto g_0_z_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 95);

    auto g_0_z_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 96);

    auto g_0_z_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 97);

    auto g_0_z_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 98);

    auto g_0_z_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 99);

    auto g_0_z_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 100);

    auto g_0_z_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 101);

    auto g_0_z_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 102);

    auto g_0_z_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 103);

    auto g_0_z_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 104);

    auto g_0_z_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 105);

    auto g_0_z_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 106);

    auto g_0_z_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 107);

    /// Set up components of auxilary buffer : SPSK

    auto g_0_x_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_spsk);

    auto g_0_x_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_spsk + 1);

    auto g_0_x_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_spsk + 2);

    auto g_0_x_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_spsk + 3);

    auto g_0_x_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_spsk + 4);

    auto g_0_x_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_spsk + 5);

    auto g_0_x_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_spsk + 6);

    auto g_0_x_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_spsk + 7);

    auto g_0_x_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_spsk + 8);

    auto g_0_x_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_spsk + 9);

    auto g_0_x_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_spsk + 10);

    auto g_0_x_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_spsk + 11);

    auto g_0_x_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_spsk + 12);

    auto g_0_x_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_spsk + 13);

    auto g_0_x_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_spsk + 14);

    auto g_0_x_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 15);

    auto g_0_x_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 16);

    auto g_0_x_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 17);

    auto g_0_x_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 18);

    auto g_0_x_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 19);

    auto g_0_x_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 20);

    auto g_0_x_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 21);

    auto g_0_x_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 22);

    auto g_0_x_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 23);

    auto g_0_x_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 24);

    auto g_0_x_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 25);

    auto g_0_x_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 26);

    auto g_0_x_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 27);

    auto g_0_x_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 28);

    auto g_0_x_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 29);

    auto g_0_x_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 30);

    auto g_0_x_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 31);

    auto g_0_x_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 32);

    auto g_0_x_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 33);

    auto g_0_x_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 34);

    auto g_0_x_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 35);

    auto g_0_y_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_spsk + 36);

    auto g_0_y_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_spsk + 37);

    auto g_0_y_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_spsk + 38);

    auto g_0_y_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_spsk + 39);

    auto g_0_y_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_spsk + 40);

    auto g_0_y_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_spsk + 41);

    auto g_0_y_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_spsk + 42);

    auto g_0_y_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_spsk + 43);

    auto g_0_y_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_spsk + 44);

    auto g_0_y_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_spsk + 45);

    auto g_0_y_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_spsk + 46);

    auto g_0_y_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_spsk + 47);

    auto g_0_y_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_spsk + 48);

    auto g_0_y_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_spsk + 49);

    auto g_0_y_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_spsk + 50);

    auto g_0_y_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 51);

    auto g_0_y_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 52);

    auto g_0_y_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 53);

    auto g_0_y_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 54);

    auto g_0_y_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 55);

    auto g_0_y_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 56);

    auto g_0_y_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 57);

    auto g_0_y_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 58);

    auto g_0_y_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 59);

    auto g_0_y_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 60);

    auto g_0_y_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 61);

    auto g_0_y_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 62);

    auto g_0_y_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 63);

    auto g_0_y_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 64);

    auto g_0_y_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 65);

    auto g_0_y_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 66);

    auto g_0_y_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 67);

    auto g_0_y_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 68);

    auto g_0_y_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 69);

    auto g_0_y_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 70);

    auto g_0_y_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 71);

    auto g_0_z_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_spsk + 72);

    auto g_0_z_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_spsk + 73);

    auto g_0_z_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_spsk + 74);

    auto g_0_z_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_spsk + 75);

    auto g_0_z_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_spsk + 76);

    auto g_0_z_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_spsk + 77);

    auto g_0_z_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_spsk + 78);

    auto g_0_z_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_spsk + 79);

    auto g_0_z_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_spsk + 80);

    auto g_0_z_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_spsk + 81);

    auto g_0_z_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_spsk + 82);

    auto g_0_z_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_spsk + 83);

    auto g_0_z_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_spsk + 84);

    auto g_0_z_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_spsk + 85);

    auto g_0_z_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_spsk + 86);

    auto g_0_z_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 87);

    auto g_0_z_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 88);

    auto g_0_z_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 89);

    auto g_0_z_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 90);

    auto g_0_z_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 91);

    auto g_0_z_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 92);

    auto g_0_z_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 93);

    auto g_0_z_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 94);

    auto g_0_z_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 95);

    auto g_0_z_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 96);

    auto g_0_z_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 97);

    auto g_0_z_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 98);

    auto g_0_z_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 99);

    auto g_0_z_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 100);

    auto g_0_z_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 101);

    auto g_0_z_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 102);

    auto g_0_z_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 103);

    auto g_0_z_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 104);

    auto g_0_z_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 105);

    auto g_0_z_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 106);

    auto g_0_z_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 107);

    /// Set up components of auxilary buffer : SDSI

    auto g_0_xx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sdsi);

    auto g_0_xx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sdsi + 1);

    auto g_0_xx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sdsi + 2);

    auto g_0_xx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sdsi + 3);

    auto g_0_xx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 4);

    auto g_0_xx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sdsi + 5);

    auto g_0_xx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sdsi + 6);

    auto g_0_xx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 7);

    auto g_0_xx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 8);

    auto g_0_xx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sdsi + 9);

    auto g_0_xx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 10);

    auto g_0_xx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 11);

    auto g_0_xx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 12);

    auto g_0_xx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 13);

    auto g_0_xx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 14);

    auto g_0_xx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 15);

    auto g_0_xx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 16);

    auto g_0_xx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 17);

    auto g_0_xx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 18);

    auto g_0_xx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 19);

    auto g_0_xx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 20);

    auto g_0_xx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 21);

    auto g_0_xx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 22);

    auto g_0_xx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 23);

    auto g_0_xx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 24);

    auto g_0_xx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 25);

    auto g_0_xx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 26);

    auto g_0_xx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 27);

    auto g_0_yy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sdsi + 84);

    auto g_0_yy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sdsi + 85);

    auto g_0_yy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sdsi + 86);

    auto g_0_yy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sdsi + 87);

    auto g_0_yy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 88);

    auto g_0_yy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sdsi + 89);

    auto g_0_yy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sdsi + 90);

    auto g_0_yy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 91);

    auto g_0_yy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 92);

    auto g_0_yy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sdsi + 93);

    auto g_0_yy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 94);

    auto g_0_yy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 95);

    auto g_0_yy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 96);

    auto g_0_yy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 97);

    auto g_0_yy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 98);

    auto g_0_yy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 99);

    auto g_0_yy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 100);

    auto g_0_yy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 101);

    auto g_0_yy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 102);

    auto g_0_yy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 103);

    auto g_0_yy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 104);

    auto g_0_yy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 105);

    auto g_0_yy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 106);

    auto g_0_yy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 107);

    auto g_0_yy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 108);

    auto g_0_yy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 109);

    auto g_0_yy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 110);

    auto g_0_yy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 111);

    auto g_0_yz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 116);

    auto g_0_yz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 119);

    auto g_0_yz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 120);

    auto g_0_yz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 123);

    auto g_0_yz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 124);

    auto g_0_yz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 125);

    auto g_0_yz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 128);

    auto g_0_yz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 129);

    auto g_0_yz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 130);

    auto g_0_yz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 131);

    auto g_0_yz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 134);

    auto g_0_yz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 135);

    auto g_0_yz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 136);

    auto g_0_yz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 137);

    auto g_0_yz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 138);

    auto g_0_zz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sdsi + 140);

    auto g_0_zz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sdsi + 141);

    auto g_0_zz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sdsi + 142);

    auto g_0_zz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sdsi + 143);

    auto g_0_zz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 144);

    auto g_0_zz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sdsi + 145);

    auto g_0_zz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sdsi + 146);

    auto g_0_zz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 147);

    auto g_0_zz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 148);

    auto g_0_zz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sdsi + 149);

    auto g_0_zz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 150);

    auto g_0_zz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 151);

    auto g_0_zz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 152);

    auto g_0_zz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 153);

    auto g_0_zz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 154);

    auto g_0_zz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 155);

    auto g_0_zz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 156);

    auto g_0_zz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 157);

    auto g_0_zz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 158);

    auto g_0_zz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 159);

    auto g_0_zz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 160);

    auto g_0_zz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 161);

    auto g_0_zz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 162);

    auto g_0_zz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 163);

    auto g_0_zz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 164);

    auto g_0_zz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 165);

    auto g_0_zz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 166);

    auto g_0_zz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 167);

    /// Set up components of auxilary buffer : SDSK

    auto g_0_xx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk);

    auto g_0_xx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 1);

    auto g_0_xx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 2);

    auto g_0_xx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 3);

    auto g_0_xx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 4);

    auto g_0_xx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 5);

    auto g_0_xx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 6);

    auto g_0_xx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 7);

    auto g_0_xx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 8);

    auto g_0_xx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 9);

    auto g_0_xx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 10);

    auto g_0_xx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 11);

    auto g_0_xx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 12);

    auto g_0_xx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 13);

    auto g_0_xx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 14);

    auto g_0_xx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 15);

    auto g_0_xx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 16);

    auto g_0_xx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 17);

    auto g_0_xx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 18);

    auto g_0_xx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 19);

    auto g_0_xx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 20);

    auto g_0_xx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 21);

    auto g_0_xx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 22);

    auto g_0_xx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 23);

    auto g_0_xx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 24);

    auto g_0_xx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 25);

    auto g_0_xx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 26);

    auto g_0_xx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 27);

    auto g_0_xx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 28);

    auto g_0_xx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 29);

    auto g_0_xx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 30);

    auto g_0_xx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 31);

    auto g_0_xx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 32);

    auto g_0_xx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 33);

    auto g_0_xx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 34);

    auto g_0_xx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 35);

    auto g_0_xy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 37);

    auto g_0_xy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 39);

    auto g_0_xy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 42);

    auto g_0_xy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 46);

    auto g_0_xy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 51);

    auto g_0_xy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 57);

    auto g_0_xz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 72);

    auto g_0_xz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 74);

    auto g_0_xz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 77);

    auto g_0_xz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 81);

    auto g_0_xz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 86);

    auto g_0_xz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 92);

    auto g_0_xz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 99);

    auto g_0_yy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 108);

    auto g_0_yy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 109);

    auto g_0_yy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 110);

    auto g_0_yy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 111);

    auto g_0_yy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 112);

    auto g_0_yy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 113);

    auto g_0_yy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 114);

    auto g_0_yy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 115);

    auto g_0_yy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 116);

    auto g_0_yy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 117);

    auto g_0_yy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 118);

    auto g_0_yy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 119);

    auto g_0_yy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 120);

    auto g_0_yy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 121);

    auto g_0_yy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 122);

    auto g_0_yy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 123);

    auto g_0_yy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 124);

    auto g_0_yy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 125);

    auto g_0_yy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 126);

    auto g_0_yy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 127);

    auto g_0_yy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 128);

    auto g_0_yy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 129);

    auto g_0_yy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 130);

    auto g_0_yy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 131);

    auto g_0_yy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 132);

    auto g_0_yy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 133);

    auto g_0_yy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 134);

    auto g_0_yy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 135);

    auto g_0_yy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 136);

    auto g_0_yy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 137);

    auto g_0_yy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 138);

    auto g_0_yy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 139);

    auto g_0_yy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 140);

    auto g_0_yy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 141);

    auto g_0_yy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 142);

    auto g_0_yy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 143);

    auto g_0_yz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 148);

    auto g_0_yz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 151);

    auto g_0_yz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 152);

    auto g_0_yz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 155);

    auto g_0_yz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 156);

    auto g_0_yz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 157);

    auto g_0_yz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 160);

    auto g_0_yz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 161);

    auto g_0_yz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 162);

    auto g_0_yz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 163);

    auto g_0_yz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 166);

    auto g_0_yz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 167);

    auto g_0_yz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 168);

    auto g_0_yz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 169);

    auto g_0_yz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 170);

    auto g_0_yz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 172);

    auto g_0_yz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 173);

    auto g_0_yz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 174);

    auto g_0_yz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 175);

    auto g_0_yz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 176);

    auto g_0_yz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 177);

    auto g_0_yz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 178);

    auto g_0_yz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 179);

    auto g_0_zz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 180);

    auto g_0_zz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 181);

    auto g_0_zz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 182);

    auto g_0_zz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 183);

    auto g_0_zz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 184);

    auto g_0_zz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 185);

    auto g_0_zz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 186);

    auto g_0_zz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 187);

    auto g_0_zz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 188);

    auto g_0_zz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 189);

    auto g_0_zz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 190);

    auto g_0_zz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 191);

    auto g_0_zz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 192);

    auto g_0_zz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 193);

    auto g_0_zz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 194);

    auto g_0_zz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 195);

    auto g_0_zz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 196);

    auto g_0_zz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 197);

    auto g_0_zz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 198);

    auto g_0_zz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 199);

    auto g_0_zz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 200);

    auto g_0_zz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 201);

    auto g_0_zz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 202);

    auto g_0_zz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 203);

    auto g_0_zz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 204);

    auto g_0_zz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 205);

    auto g_0_zz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 206);

    auto g_0_zz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 207);

    auto g_0_zz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 208);

    auto g_0_zz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 209);

    auto g_0_zz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 210);

    auto g_0_zz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 211);

    auto g_0_zz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 212);

    auto g_0_zz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 213);

    auto g_0_zz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 214);

    auto g_0_zz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 215);

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

    auto g_0_xy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sdsk + 37);

    auto g_0_xy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sdsk + 39);

    auto g_0_xy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sdsk + 42);

    auto g_0_xy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 46);

    auto g_0_xy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 51);

    auto g_0_xy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 57);

    auto g_0_xz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sdsk + 72);

    auto g_0_xz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sdsk + 74);

    auto g_0_xz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sdsk + 77);

    auto g_0_xz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sdsk + 81);

    auto g_0_xz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 86);

    auto g_0_xz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 92);

    auto g_0_xz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 99);

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

    auto g_0_yz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sdsk + 172);

    auto g_0_yz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sdsk + 173);

    auto g_0_yz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sdsk + 174);

    auto g_0_yz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sdsk + 175);

    auto g_0_yz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 176);

    auto g_0_yz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 177);

    auto g_0_yz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 178);

    auto g_0_yz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sdsk + 179);

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

    /// Set up 0-36 components of targeted buffer : SFSK

    auto g_0_xxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk);

    auto g_0_xxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 1);

    auto g_0_xxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 2);

    auto g_0_xxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 3);

    auto g_0_xxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 4);

    auto g_0_xxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 5);

    auto g_0_xxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 6);

    auto g_0_xxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 7);

    auto g_0_xxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 8);

    auto g_0_xxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 9);

    auto g_0_xxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 10);

    auto g_0_xxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 11);

    auto g_0_xxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 12);

    auto g_0_xxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 13);

    auto g_0_xxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 14);

    auto g_0_xxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 15);

    auto g_0_xxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 16);

    auto g_0_xxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 17);

    auto g_0_xxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 18);

    auto g_0_xxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 19);

    auto g_0_xxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 20);

    auto g_0_xxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 21);

    auto g_0_xxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 22);

    auto g_0_xxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 23);

    auto g_0_xxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 24);

    auto g_0_xxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 25);

    auto g_0_xxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 26);

    auto g_0_xxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 27);

    auto g_0_xxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 28);

    auto g_0_xxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 29);

    auto g_0_xxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 30);

    auto g_0_xxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 31);

    auto g_0_xxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 32);

    auto g_0_xxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 33);

    auto g_0_xxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 34);

    auto g_0_xxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 35);

#pragma omp simd aligned(g_0_x_0_xxxxxxx_0,       \
                             g_0_x_0_xxxxxxx_1,   \
                             g_0_x_0_xxxxxxy_0,   \
                             g_0_x_0_xxxxxxy_1,   \
                             g_0_x_0_xxxxxxz_0,   \
                             g_0_x_0_xxxxxxz_1,   \
                             g_0_x_0_xxxxxyy_0,   \
                             g_0_x_0_xxxxxyy_1,   \
                             g_0_x_0_xxxxxyz_0,   \
                             g_0_x_0_xxxxxyz_1,   \
                             g_0_x_0_xxxxxzz_0,   \
                             g_0_x_0_xxxxxzz_1,   \
                             g_0_x_0_xxxxyyy_0,   \
                             g_0_x_0_xxxxyyy_1,   \
                             g_0_x_0_xxxxyyz_0,   \
                             g_0_x_0_xxxxyyz_1,   \
                             g_0_x_0_xxxxyzz_0,   \
                             g_0_x_0_xxxxyzz_1,   \
                             g_0_x_0_xxxxzzz_0,   \
                             g_0_x_0_xxxxzzz_1,   \
                             g_0_x_0_xxxyyyy_0,   \
                             g_0_x_0_xxxyyyy_1,   \
                             g_0_x_0_xxxyyyz_0,   \
                             g_0_x_0_xxxyyyz_1,   \
                             g_0_x_0_xxxyyzz_0,   \
                             g_0_x_0_xxxyyzz_1,   \
                             g_0_x_0_xxxyzzz_0,   \
                             g_0_x_0_xxxyzzz_1,   \
                             g_0_x_0_xxxzzzz_0,   \
                             g_0_x_0_xxxzzzz_1,   \
                             g_0_x_0_xxyyyyy_0,   \
                             g_0_x_0_xxyyyyy_1,   \
                             g_0_x_0_xxyyyyz_0,   \
                             g_0_x_0_xxyyyyz_1,   \
                             g_0_x_0_xxyyyzz_0,   \
                             g_0_x_0_xxyyyzz_1,   \
                             g_0_x_0_xxyyzzz_0,   \
                             g_0_x_0_xxyyzzz_1,   \
                             g_0_x_0_xxyzzzz_0,   \
                             g_0_x_0_xxyzzzz_1,   \
                             g_0_x_0_xxzzzzz_0,   \
                             g_0_x_0_xxzzzzz_1,   \
                             g_0_x_0_xyyyyyy_0,   \
                             g_0_x_0_xyyyyyy_1,   \
                             g_0_x_0_xyyyyyz_0,   \
                             g_0_x_0_xyyyyyz_1,   \
                             g_0_x_0_xyyyyzz_0,   \
                             g_0_x_0_xyyyyzz_1,   \
                             g_0_x_0_xyyyzzz_0,   \
                             g_0_x_0_xyyyzzz_1,   \
                             g_0_x_0_xyyzzzz_0,   \
                             g_0_x_0_xyyzzzz_1,   \
                             g_0_x_0_xyzzzzz_0,   \
                             g_0_x_0_xyzzzzz_1,   \
                             g_0_x_0_xzzzzzz_0,   \
                             g_0_x_0_xzzzzzz_1,   \
                             g_0_x_0_yyyyyyy_0,   \
                             g_0_x_0_yyyyyyy_1,   \
                             g_0_x_0_yyyyyyz_0,   \
                             g_0_x_0_yyyyyyz_1,   \
                             g_0_x_0_yyyyyzz_0,   \
                             g_0_x_0_yyyyyzz_1,   \
                             g_0_x_0_yyyyzzz_0,   \
                             g_0_x_0_yyyyzzz_1,   \
                             g_0_x_0_yyyzzzz_0,   \
                             g_0_x_0_yyyzzzz_1,   \
                             g_0_x_0_yyzzzzz_0,   \
                             g_0_x_0_yyzzzzz_1,   \
                             g_0_x_0_yzzzzzz_0,   \
                             g_0_x_0_yzzzzzz_1,   \
                             g_0_x_0_zzzzzzz_0,   \
                             g_0_x_0_zzzzzzz_1,   \
                             g_0_xx_0_xxxxxx_1,   \
                             g_0_xx_0_xxxxxxx_0,  \
                             g_0_xx_0_xxxxxxx_1,  \
                             g_0_xx_0_xxxxxxy_0,  \
                             g_0_xx_0_xxxxxxy_1,  \
                             g_0_xx_0_xxxxxxz_0,  \
                             g_0_xx_0_xxxxxxz_1,  \
                             g_0_xx_0_xxxxxy_1,   \
                             g_0_xx_0_xxxxxyy_0,  \
                             g_0_xx_0_xxxxxyy_1,  \
                             g_0_xx_0_xxxxxyz_0,  \
                             g_0_xx_0_xxxxxyz_1,  \
                             g_0_xx_0_xxxxxz_1,   \
                             g_0_xx_0_xxxxxzz_0,  \
                             g_0_xx_0_xxxxxzz_1,  \
                             g_0_xx_0_xxxxyy_1,   \
                             g_0_xx_0_xxxxyyy_0,  \
                             g_0_xx_0_xxxxyyy_1,  \
                             g_0_xx_0_xxxxyyz_0,  \
                             g_0_xx_0_xxxxyyz_1,  \
                             g_0_xx_0_xxxxyz_1,   \
                             g_0_xx_0_xxxxyzz_0,  \
                             g_0_xx_0_xxxxyzz_1,  \
                             g_0_xx_0_xxxxzz_1,   \
                             g_0_xx_0_xxxxzzz_0,  \
                             g_0_xx_0_xxxxzzz_1,  \
                             g_0_xx_0_xxxyyy_1,   \
                             g_0_xx_0_xxxyyyy_0,  \
                             g_0_xx_0_xxxyyyy_1,  \
                             g_0_xx_0_xxxyyyz_0,  \
                             g_0_xx_0_xxxyyyz_1,  \
                             g_0_xx_0_xxxyyz_1,   \
                             g_0_xx_0_xxxyyzz_0,  \
                             g_0_xx_0_xxxyyzz_1,  \
                             g_0_xx_0_xxxyzz_1,   \
                             g_0_xx_0_xxxyzzz_0,  \
                             g_0_xx_0_xxxyzzz_1,  \
                             g_0_xx_0_xxxzzz_1,   \
                             g_0_xx_0_xxxzzzz_0,  \
                             g_0_xx_0_xxxzzzz_1,  \
                             g_0_xx_0_xxyyyy_1,   \
                             g_0_xx_0_xxyyyyy_0,  \
                             g_0_xx_0_xxyyyyy_1,  \
                             g_0_xx_0_xxyyyyz_0,  \
                             g_0_xx_0_xxyyyyz_1,  \
                             g_0_xx_0_xxyyyz_1,   \
                             g_0_xx_0_xxyyyzz_0,  \
                             g_0_xx_0_xxyyyzz_1,  \
                             g_0_xx_0_xxyyzz_1,   \
                             g_0_xx_0_xxyyzzz_0,  \
                             g_0_xx_0_xxyyzzz_1,  \
                             g_0_xx_0_xxyzzz_1,   \
                             g_0_xx_0_xxyzzzz_0,  \
                             g_0_xx_0_xxyzzzz_1,  \
                             g_0_xx_0_xxzzzz_1,   \
                             g_0_xx_0_xxzzzzz_0,  \
                             g_0_xx_0_xxzzzzz_1,  \
                             g_0_xx_0_xyyyyy_1,   \
                             g_0_xx_0_xyyyyyy_0,  \
                             g_0_xx_0_xyyyyyy_1,  \
                             g_0_xx_0_xyyyyyz_0,  \
                             g_0_xx_0_xyyyyyz_1,  \
                             g_0_xx_0_xyyyyz_1,   \
                             g_0_xx_0_xyyyyzz_0,  \
                             g_0_xx_0_xyyyyzz_1,  \
                             g_0_xx_0_xyyyzz_1,   \
                             g_0_xx_0_xyyyzzz_0,  \
                             g_0_xx_0_xyyyzzz_1,  \
                             g_0_xx_0_xyyzzz_1,   \
                             g_0_xx_0_xyyzzzz_0,  \
                             g_0_xx_0_xyyzzzz_1,  \
                             g_0_xx_0_xyzzzz_1,   \
                             g_0_xx_0_xyzzzzz_0,  \
                             g_0_xx_0_xyzzzzz_1,  \
                             g_0_xx_0_xzzzzz_1,   \
                             g_0_xx_0_xzzzzzz_0,  \
                             g_0_xx_0_xzzzzzz_1,  \
                             g_0_xx_0_yyyyyy_1,   \
                             g_0_xx_0_yyyyyyy_0,  \
                             g_0_xx_0_yyyyyyy_1,  \
                             g_0_xx_0_yyyyyyz_0,  \
                             g_0_xx_0_yyyyyyz_1,  \
                             g_0_xx_0_yyyyyz_1,   \
                             g_0_xx_0_yyyyyzz_0,  \
                             g_0_xx_0_yyyyyzz_1,  \
                             g_0_xx_0_yyyyzz_1,   \
                             g_0_xx_0_yyyyzzz_0,  \
                             g_0_xx_0_yyyyzzz_1,  \
                             g_0_xx_0_yyyzzz_1,   \
                             g_0_xx_0_yyyzzzz_0,  \
                             g_0_xx_0_yyyzzzz_1,  \
                             g_0_xx_0_yyzzzz_1,   \
                             g_0_xx_0_yyzzzzz_0,  \
                             g_0_xx_0_yyzzzzz_1,  \
                             g_0_xx_0_yzzzzz_1,   \
                             g_0_xx_0_yzzzzzz_0,  \
                             g_0_xx_0_yzzzzzz_1,  \
                             g_0_xx_0_zzzzzz_1,   \
                             g_0_xx_0_zzzzzzz_0,  \
                             g_0_xx_0_zzzzzzz_1,  \
                             g_0_xxx_0_xxxxxxx_0, \
                             g_0_xxx_0_xxxxxxy_0, \
                             g_0_xxx_0_xxxxxxz_0, \
                             g_0_xxx_0_xxxxxyy_0, \
                             g_0_xxx_0_xxxxxyz_0, \
                             g_0_xxx_0_xxxxxzz_0, \
                             g_0_xxx_0_xxxxyyy_0, \
                             g_0_xxx_0_xxxxyyz_0, \
                             g_0_xxx_0_xxxxyzz_0, \
                             g_0_xxx_0_xxxxzzz_0, \
                             g_0_xxx_0_xxxyyyy_0, \
                             g_0_xxx_0_xxxyyyz_0, \
                             g_0_xxx_0_xxxyyzz_0, \
                             g_0_xxx_0_xxxyzzz_0, \
                             g_0_xxx_0_xxxzzzz_0, \
                             g_0_xxx_0_xxyyyyy_0, \
                             g_0_xxx_0_xxyyyyz_0, \
                             g_0_xxx_0_xxyyyzz_0, \
                             g_0_xxx_0_xxyyzzz_0, \
                             g_0_xxx_0_xxyzzzz_0, \
                             g_0_xxx_0_xxzzzzz_0, \
                             g_0_xxx_0_xyyyyyy_0, \
                             g_0_xxx_0_xyyyyyz_0, \
                             g_0_xxx_0_xyyyyzz_0, \
                             g_0_xxx_0_xyyyzzz_0, \
                             g_0_xxx_0_xyyzzzz_0, \
                             g_0_xxx_0_xyzzzzz_0, \
                             g_0_xxx_0_xzzzzzz_0, \
                             g_0_xxx_0_yyyyyyy_0, \
                             g_0_xxx_0_yyyyyyz_0, \
                             g_0_xxx_0_yyyyyzz_0, \
                             g_0_xxx_0_yyyyzzz_0, \
                             g_0_xxx_0_yyyzzzz_0, \
                             g_0_xxx_0_yyzzzzz_0, \
                             g_0_xxx_0_yzzzzzz_0, \
                             g_0_xxx_0_zzzzzzz_0, \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xxxxxxx_0[i] = 2.0 * g_0_x_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxx_1[i] * fti_ab_0 +
                                 7.0 * g_0_xx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxx_0[i] * pb_x + g_0_xx_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxxy_0[i] = 2.0 * g_0_x_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxy_1[i] * fti_ab_0 +
                                 6.0 * g_0_xx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxy_0[i] * pb_x + g_0_xx_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxxz_0[i] = 2.0 * g_0_x_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxxz_1[i] * fti_ab_0 +
                                 6.0 * g_0_xx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxz_0[i] * pb_x + g_0_xx_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxyy_0[i] = 2.0 * g_0_x_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxyy_1[i] * fti_ab_0 +
                                 5.0 * g_0_xx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyy_0[i] * pb_x + g_0_xx_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxyz_0[i] = 2.0 * g_0_x_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxyz_1[i] * fti_ab_0 +
                                 5.0 * g_0_xx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyz_0[i] * pb_x + g_0_xx_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxzz_0[i] = 2.0 * g_0_x_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxzz_1[i] * fti_ab_0 +
                                 5.0 * g_0_xx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxzz_0[i] * pb_x + g_0_xx_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyyy_0[i] = 2.0 * g_0_x_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyyy_1[i] * fti_ab_0 +
                                 4.0 * g_0_xx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyy_0[i] * pb_x + g_0_xx_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyyz_0[i] = 2.0 * g_0_x_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyyz_1[i] * fti_ab_0 +
                                 4.0 * g_0_xx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyz_0[i] * pb_x + g_0_xx_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyzz_0[i] = 2.0 * g_0_x_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_xx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyzz_0[i] * pb_x + g_0_xx_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxzzz_0[i] = 2.0 * g_0_x_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxzzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_xx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxzzz_0[i] * pb_x + g_0_xx_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyyy_0[i] = 2.0 * g_0_x_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyyy_1[i] * fti_ab_0 +
                                 3.0 * g_0_xx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyy_0[i] * pb_x + g_0_xx_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyyz_0[i] = 2.0 * g_0_x_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyyz_1[i] * fti_ab_0 +
                                 3.0 * g_0_xx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyz_0[i] * pb_x + g_0_xx_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyzz_0[i] = 2.0 * g_0_x_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_xx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyzz_0[i] * pb_x + g_0_xx_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyzzz_0[i] = 2.0 * g_0_x_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_xx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzzz_0[i] * pb_x + g_0_xx_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxzzzz_0[i] = 2.0 * g_0_x_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxzzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_xx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxzzzz_0[i] * pb_x + g_0_xx_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyyy_0[i] = 2.0 * g_0_x_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyyy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyy_0[i] * pb_x + g_0_xx_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyyz_0[i] = 2.0 * g_0_x_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyz_0[i] * pb_x + g_0_xx_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyzz_0[i] = 2.0 * g_0_x_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyzz_0[i] * pb_x + g_0_xx_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyzzz_0[i] = 2.0 * g_0_x_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyzzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzzz_0[i] * pb_x + g_0_xx_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyzzzz_0[i] = 2.0 * g_0_x_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyzzzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzzz_0[i] * pb_x + g_0_xx_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxzzzzz_0[i] = 2.0 * g_0_x_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxzzzzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzzzzz_0[i] * pb_x + g_0_xx_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyyy_0[i] = 2.0 * g_0_x_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyyyy_1[i] * fi_abcd_0 +
                                 g_0_xx_0_xyyyyyy_0[i] * pb_x + g_0_xx_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyyz_0[i] = 2.0 * g_0_x_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyz_1[i] * fi_abcd_0 +
                                 g_0_xx_0_xyyyyyz_0[i] * pb_x + g_0_xx_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyzz_0[i] = 2.0 * g_0_x_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyzz_1[i] * fi_abcd_0 +
                                 g_0_xx_0_xyyyyzz_0[i] * pb_x + g_0_xx_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyzzz_0[i] = 2.0 * g_0_x_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyzzz_1[i] * fi_abcd_0 +
                                 g_0_xx_0_xyyyzzz_0[i] * pb_x + g_0_xx_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyzzzz_0[i] = 2.0 * g_0_x_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyzzzz_1[i] * fi_abcd_0 +
                                 g_0_xx_0_xyyzzzz_0[i] * pb_x + g_0_xx_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyzzzzz_0[i] = 2.0 * g_0_x_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzzzz_1[i] * fi_abcd_0 +
                                 g_0_xx_0_xyzzzzz_0[i] * pb_x + g_0_xx_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xzzzzzz_0[i] = 2.0 * g_0_x_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzzzz_1[i] * fi_abcd_0 +
                                 g_0_xx_0_xzzzzzz_0[i] * pb_x + g_0_xx_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyyy_0[i] = 2.0 * g_0_x_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyyyyy_0[i] * pb_x +
                                 g_0_xx_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyyz_0[i] = 2.0 * g_0_x_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyyz_0[i] * pb_x +
                                 g_0_xx_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyzz_0[i] = 2.0 * g_0_x_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyzz_0[i] * pb_x +
                                 g_0_xx_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyzzz_0[i] = 2.0 * g_0_x_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyzzz_0[i] * pb_x +
                                 g_0_xx_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyzzzz_0[i] = 2.0 * g_0_x_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyzzzz_0[i] * pb_x +
                                 g_0_xx_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyzzzzz_0[i] = 2.0 * g_0_x_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyzzzzz_0[i] * pb_x +
                                 g_0_xx_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yzzzzzz_0[i] = 2.0 * g_0_x_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzzzzz_0[i] * pb_x +
                                 g_0_xx_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_zzzzzzz_0[i] = 2.0 * g_0_x_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzzzzz_0[i] * pb_x +
                                 g_0_xx_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : SFSK

    auto g_0_xxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 36);

    auto g_0_xxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 37);

    auto g_0_xxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 38);

    auto g_0_xxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 39);

    auto g_0_xxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 40);

    auto g_0_xxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 41);

    auto g_0_xxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 42);

    auto g_0_xxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 43);

    auto g_0_xxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 44);

    auto g_0_xxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 45);

    auto g_0_xxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 46);

    auto g_0_xxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 47);

    auto g_0_xxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 48);

    auto g_0_xxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 49);

    auto g_0_xxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 50);

    auto g_0_xxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 51);

    auto g_0_xxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 52);

    auto g_0_xxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 53);

    auto g_0_xxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 54);

    auto g_0_xxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 55);

    auto g_0_xxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 56);

    auto g_0_xxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 57);

    auto g_0_xxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 58);

    auto g_0_xxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 59);

    auto g_0_xxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 60);

    auto g_0_xxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 61);

    auto g_0_xxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 62);

    auto g_0_xxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 63);

    auto g_0_xxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 64);

    auto g_0_xxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 65);

    auto g_0_xxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 66);

    auto g_0_xxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 67);

    auto g_0_xxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 68);

    auto g_0_xxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 69);

    auto g_0_xxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 70);

    auto g_0_xxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 71);

#pragma omp simd aligned(g_0_xx_0_xxxxxx_1,       \
                             g_0_xx_0_xxxxxxx_0,  \
                             g_0_xx_0_xxxxxxx_1,  \
                             g_0_xx_0_xxxxxxy_0,  \
                             g_0_xx_0_xxxxxxy_1,  \
                             g_0_xx_0_xxxxxxz_0,  \
                             g_0_xx_0_xxxxxxz_1,  \
                             g_0_xx_0_xxxxxy_1,   \
                             g_0_xx_0_xxxxxyy_0,  \
                             g_0_xx_0_xxxxxyy_1,  \
                             g_0_xx_0_xxxxxyz_0,  \
                             g_0_xx_0_xxxxxyz_1,  \
                             g_0_xx_0_xxxxxz_1,   \
                             g_0_xx_0_xxxxxzz_0,  \
                             g_0_xx_0_xxxxxzz_1,  \
                             g_0_xx_0_xxxxyy_1,   \
                             g_0_xx_0_xxxxyyy_0,  \
                             g_0_xx_0_xxxxyyy_1,  \
                             g_0_xx_0_xxxxyyz_0,  \
                             g_0_xx_0_xxxxyyz_1,  \
                             g_0_xx_0_xxxxyz_1,   \
                             g_0_xx_0_xxxxyzz_0,  \
                             g_0_xx_0_xxxxyzz_1,  \
                             g_0_xx_0_xxxxzz_1,   \
                             g_0_xx_0_xxxxzzz_0,  \
                             g_0_xx_0_xxxxzzz_1,  \
                             g_0_xx_0_xxxyyy_1,   \
                             g_0_xx_0_xxxyyyy_0,  \
                             g_0_xx_0_xxxyyyy_1,  \
                             g_0_xx_0_xxxyyyz_0,  \
                             g_0_xx_0_xxxyyyz_1,  \
                             g_0_xx_0_xxxyyz_1,   \
                             g_0_xx_0_xxxyyzz_0,  \
                             g_0_xx_0_xxxyyzz_1,  \
                             g_0_xx_0_xxxyzz_1,   \
                             g_0_xx_0_xxxyzzz_0,  \
                             g_0_xx_0_xxxyzzz_1,  \
                             g_0_xx_0_xxxzzz_1,   \
                             g_0_xx_0_xxxzzzz_0,  \
                             g_0_xx_0_xxxzzzz_1,  \
                             g_0_xx_0_xxyyyy_1,   \
                             g_0_xx_0_xxyyyyy_0,  \
                             g_0_xx_0_xxyyyyy_1,  \
                             g_0_xx_0_xxyyyyz_0,  \
                             g_0_xx_0_xxyyyyz_1,  \
                             g_0_xx_0_xxyyyz_1,   \
                             g_0_xx_0_xxyyyzz_0,  \
                             g_0_xx_0_xxyyyzz_1,  \
                             g_0_xx_0_xxyyzz_1,   \
                             g_0_xx_0_xxyyzzz_0,  \
                             g_0_xx_0_xxyyzzz_1,  \
                             g_0_xx_0_xxyzzz_1,   \
                             g_0_xx_0_xxyzzzz_0,  \
                             g_0_xx_0_xxyzzzz_1,  \
                             g_0_xx_0_xxzzzz_1,   \
                             g_0_xx_0_xxzzzzz_0,  \
                             g_0_xx_0_xxzzzzz_1,  \
                             g_0_xx_0_xyyyyy_1,   \
                             g_0_xx_0_xyyyyyy_0,  \
                             g_0_xx_0_xyyyyyy_1,  \
                             g_0_xx_0_xyyyyyz_0,  \
                             g_0_xx_0_xyyyyyz_1,  \
                             g_0_xx_0_xyyyyz_1,   \
                             g_0_xx_0_xyyyyzz_0,  \
                             g_0_xx_0_xyyyyzz_1,  \
                             g_0_xx_0_xyyyzz_1,   \
                             g_0_xx_0_xyyyzzz_0,  \
                             g_0_xx_0_xyyyzzz_1,  \
                             g_0_xx_0_xyyzzz_1,   \
                             g_0_xx_0_xyyzzzz_0,  \
                             g_0_xx_0_xyyzzzz_1,  \
                             g_0_xx_0_xyzzzz_1,   \
                             g_0_xx_0_xyzzzzz_0,  \
                             g_0_xx_0_xyzzzzz_1,  \
                             g_0_xx_0_xzzzzz_1,   \
                             g_0_xx_0_xzzzzzz_0,  \
                             g_0_xx_0_xzzzzzz_1,  \
                             g_0_xx_0_yyyyyy_1,   \
                             g_0_xx_0_yyyyyyy_0,  \
                             g_0_xx_0_yyyyyyy_1,  \
                             g_0_xx_0_yyyyyyz_0,  \
                             g_0_xx_0_yyyyyyz_1,  \
                             g_0_xx_0_yyyyyz_1,   \
                             g_0_xx_0_yyyyyzz_0,  \
                             g_0_xx_0_yyyyyzz_1,  \
                             g_0_xx_0_yyyyzz_1,   \
                             g_0_xx_0_yyyyzzz_0,  \
                             g_0_xx_0_yyyyzzz_1,  \
                             g_0_xx_0_yyyzzz_1,   \
                             g_0_xx_0_yyyzzzz_0,  \
                             g_0_xx_0_yyyzzzz_1,  \
                             g_0_xx_0_yyzzzz_1,   \
                             g_0_xx_0_yyzzzzz_0,  \
                             g_0_xx_0_yyzzzzz_1,  \
                             g_0_xx_0_yzzzzz_1,   \
                             g_0_xx_0_yzzzzzz_0,  \
                             g_0_xx_0_yzzzzzz_1,  \
                             g_0_xx_0_zzzzzz_1,   \
                             g_0_xx_0_zzzzzzz_0,  \
                             g_0_xx_0_zzzzzzz_1,  \
                             g_0_xxy_0_xxxxxxx_0, \
                             g_0_xxy_0_xxxxxxy_0, \
                             g_0_xxy_0_xxxxxxz_0, \
                             g_0_xxy_0_xxxxxyy_0, \
                             g_0_xxy_0_xxxxxyz_0, \
                             g_0_xxy_0_xxxxxzz_0, \
                             g_0_xxy_0_xxxxyyy_0, \
                             g_0_xxy_0_xxxxyyz_0, \
                             g_0_xxy_0_xxxxyzz_0, \
                             g_0_xxy_0_xxxxzzz_0, \
                             g_0_xxy_0_xxxyyyy_0, \
                             g_0_xxy_0_xxxyyyz_0, \
                             g_0_xxy_0_xxxyyzz_0, \
                             g_0_xxy_0_xxxyzzz_0, \
                             g_0_xxy_0_xxxzzzz_0, \
                             g_0_xxy_0_xxyyyyy_0, \
                             g_0_xxy_0_xxyyyyz_0, \
                             g_0_xxy_0_xxyyyzz_0, \
                             g_0_xxy_0_xxyyzzz_0, \
                             g_0_xxy_0_xxyzzzz_0, \
                             g_0_xxy_0_xxzzzzz_0, \
                             g_0_xxy_0_xyyyyyy_0, \
                             g_0_xxy_0_xyyyyyz_0, \
                             g_0_xxy_0_xyyyyzz_0, \
                             g_0_xxy_0_xyyyzzz_0, \
                             g_0_xxy_0_xyyzzzz_0, \
                             g_0_xxy_0_xyzzzzz_0, \
                             g_0_xxy_0_xzzzzzz_0, \
                             g_0_xxy_0_yyyyyyy_0, \
                             g_0_xxy_0_yyyyyyz_0, \
                             g_0_xxy_0_yyyyyzz_0, \
                             g_0_xxy_0_yyyyzzz_0, \
                             g_0_xxy_0_yyyzzzz_0, \
                             g_0_xxy_0_yyzzzzz_0, \
                             g_0_xxy_0_yzzzzzz_0, \
                             g_0_xxy_0_zzzzzzz_0, \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xxxxxxx_0[i] = g_0_xx_0_xxxxxxx_0[i] * pb_y + g_0_xx_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxxy_0[i] = g_0_xx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxy_0[i] * pb_y + g_0_xx_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxxz_0[i] = g_0_xx_0_xxxxxxz_0[i] * pb_y + g_0_xx_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxyy_0[i] = 2.0 * g_0_xx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyy_0[i] * pb_y + g_0_xx_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxyz_0[i] = g_0_xx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyz_0[i] * pb_y + g_0_xx_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxzz_0[i] = g_0_xx_0_xxxxxzz_0[i] * pb_y + g_0_xx_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyyy_0[i] = 3.0 * g_0_xx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyy_0[i] * pb_y + g_0_xx_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyyz_0[i] = 2.0 * g_0_xx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyz_0[i] * pb_y + g_0_xx_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyzz_0[i] = g_0_xx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyzz_0[i] * pb_y + g_0_xx_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxzzz_0[i] = g_0_xx_0_xxxxzzz_0[i] * pb_y + g_0_xx_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyyy_0[i] = 4.0 * g_0_xx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyy_0[i] * pb_y + g_0_xx_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyyz_0[i] = 3.0 * g_0_xx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyz_0[i] * pb_y + g_0_xx_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyzz_0[i] = 2.0 * g_0_xx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyzz_0[i] * pb_y + g_0_xx_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyzzz_0[i] = g_0_xx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzzz_0[i] * pb_y + g_0_xx_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxzzzz_0[i] = g_0_xx_0_xxxzzzz_0[i] * pb_y + g_0_xx_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyyy_0[i] = 5.0 * g_0_xx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyy_0[i] * pb_y + g_0_xx_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyyz_0[i] = 4.0 * g_0_xx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyz_0[i] * pb_y + g_0_xx_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyzz_0[i] = 3.0 * g_0_xx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyzz_0[i] * pb_y + g_0_xx_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyzzz_0[i] = 2.0 * g_0_xx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzzz_0[i] * pb_y + g_0_xx_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyzzzz_0[i] = g_0_xx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzzz_0[i] * pb_y + g_0_xx_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxzzzzz_0[i] = g_0_xx_0_xxzzzzz_0[i] * pb_y + g_0_xx_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyyy_0[i] = 6.0 * g_0_xx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyy_0[i] * pb_y + g_0_xx_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyyz_0[i] = 5.0 * g_0_xx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyz_0[i] * pb_y + g_0_xx_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyzz_0[i] = 4.0 * g_0_xx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyzz_0[i] * pb_y + g_0_xx_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyzzz_0[i] = 3.0 * g_0_xx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyzzz_0[i] * pb_y + g_0_xx_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyzzzz_0[i] = 2.0 * g_0_xx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzzzz_0[i] * pb_y + g_0_xx_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyzzzzz_0[i] = g_0_xx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzzzz_0[i] * pb_y + g_0_xx_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xzzzzzz_0[i] = g_0_xx_0_xzzzzzz_0[i] * pb_y + g_0_xx_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyyy_0[i] = 7.0 * g_0_xx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyy_0[i] * pb_y + g_0_xx_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyyz_0[i] = 6.0 * g_0_xx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyz_0[i] * pb_y + g_0_xx_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyzz_0[i] = 5.0 * g_0_xx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyzz_0[i] * pb_y + g_0_xx_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyzzz_0[i] = 4.0 * g_0_xx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyzzz_0[i] * pb_y + g_0_xx_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyzzzz_0[i] = 3.0 * g_0_xx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzzzz_0[i] * pb_y + g_0_xx_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyzzzzz_0[i] = 2.0 * g_0_xx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzzzz_0[i] * pb_y + g_0_xx_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yzzzzzz_0[i] = g_0_xx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzzzz_0[i] * pb_y + g_0_xx_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_zzzzzzz_0[i] = g_0_xx_0_zzzzzzz_0[i] * pb_y + g_0_xx_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 72-108 components of targeted buffer : SFSK

    auto g_0_xxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 72);

    auto g_0_xxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 73);

    auto g_0_xxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 74);

    auto g_0_xxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 75);

    auto g_0_xxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 76);

    auto g_0_xxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 77);

    auto g_0_xxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 78);

    auto g_0_xxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 79);

    auto g_0_xxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 80);

    auto g_0_xxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 81);

    auto g_0_xxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 82);

    auto g_0_xxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 83);

    auto g_0_xxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 84);

    auto g_0_xxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 85);

    auto g_0_xxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 86);

    auto g_0_xxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 87);

    auto g_0_xxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 88);

    auto g_0_xxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 89);

    auto g_0_xxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 90);

    auto g_0_xxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 91);

    auto g_0_xxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 92);

    auto g_0_xxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 93);

    auto g_0_xxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 94);

    auto g_0_xxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 95);

    auto g_0_xxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 96);

    auto g_0_xxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 97);

    auto g_0_xxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 98);

    auto g_0_xxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 99);

    auto g_0_xxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 100);

    auto g_0_xxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 101);

    auto g_0_xxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 102);

    auto g_0_xxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 103);

    auto g_0_xxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 104);

    auto g_0_xxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 105);

    auto g_0_xxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 106);

    auto g_0_xxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 107);

#pragma omp simd aligned(g_0_xx_0_xxxxxx_1,       \
                             g_0_xx_0_xxxxxxx_0,  \
                             g_0_xx_0_xxxxxxx_1,  \
                             g_0_xx_0_xxxxxxy_0,  \
                             g_0_xx_0_xxxxxxy_1,  \
                             g_0_xx_0_xxxxxxz_0,  \
                             g_0_xx_0_xxxxxxz_1,  \
                             g_0_xx_0_xxxxxy_1,   \
                             g_0_xx_0_xxxxxyy_0,  \
                             g_0_xx_0_xxxxxyy_1,  \
                             g_0_xx_0_xxxxxyz_0,  \
                             g_0_xx_0_xxxxxyz_1,  \
                             g_0_xx_0_xxxxxz_1,   \
                             g_0_xx_0_xxxxxzz_0,  \
                             g_0_xx_0_xxxxxzz_1,  \
                             g_0_xx_0_xxxxyy_1,   \
                             g_0_xx_0_xxxxyyy_0,  \
                             g_0_xx_0_xxxxyyy_1,  \
                             g_0_xx_0_xxxxyyz_0,  \
                             g_0_xx_0_xxxxyyz_1,  \
                             g_0_xx_0_xxxxyz_1,   \
                             g_0_xx_0_xxxxyzz_0,  \
                             g_0_xx_0_xxxxyzz_1,  \
                             g_0_xx_0_xxxxzz_1,   \
                             g_0_xx_0_xxxxzzz_0,  \
                             g_0_xx_0_xxxxzzz_1,  \
                             g_0_xx_0_xxxyyy_1,   \
                             g_0_xx_0_xxxyyyy_0,  \
                             g_0_xx_0_xxxyyyy_1,  \
                             g_0_xx_0_xxxyyyz_0,  \
                             g_0_xx_0_xxxyyyz_1,  \
                             g_0_xx_0_xxxyyz_1,   \
                             g_0_xx_0_xxxyyzz_0,  \
                             g_0_xx_0_xxxyyzz_1,  \
                             g_0_xx_0_xxxyzz_1,   \
                             g_0_xx_0_xxxyzzz_0,  \
                             g_0_xx_0_xxxyzzz_1,  \
                             g_0_xx_0_xxxzzz_1,   \
                             g_0_xx_0_xxxzzzz_0,  \
                             g_0_xx_0_xxxzzzz_1,  \
                             g_0_xx_0_xxyyyy_1,   \
                             g_0_xx_0_xxyyyyy_0,  \
                             g_0_xx_0_xxyyyyy_1,  \
                             g_0_xx_0_xxyyyyz_0,  \
                             g_0_xx_0_xxyyyyz_1,  \
                             g_0_xx_0_xxyyyz_1,   \
                             g_0_xx_0_xxyyyzz_0,  \
                             g_0_xx_0_xxyyyzz_1,  \
                             g_0_xx_0_xxyyzz_1,   \
                             g_0_xx_0_xxyyzzz_0,  \
                             g_0_xx_0_xxyyzzz_1,  \
                             g_0_xx_0_xxyzzz_1,   \
                             g_0_xx_0_xxyzzzz_0,  \
                             g_0_xx_0_xxyzzzz_1,  \
                             g_0_xx_0_xxzzzz_1,   \
                             g_0_xx_0_xxzzzzz_0,  \
                             g_0_xx_0_xxzzzzz_1,  \
                             g_0_xx_0_xyyyyy_1,   \
                             g_0_xx_0_xyyyyyy_0,  \
                             g_0_xx_0_xyyyyyy_1,  \
                             g_0_xx_0_xyyyyyz_0,  \
                             g_0_xx_0_xyyyyyz_1,  \
                             g_0_xx_0_xyyyyz_1,   \
                             g_0_xx_0_xyyyyzz_0,  \
                             g_0_xx_0_xyyyyzz_1,  \
                             g_0_xx_0_xyyyzz_1,   \
                             g_0_xx_0_xyyyzzz_0,  \
                             g_0_xx_0_xyyyzzz_1,  \
                             g_0_xx_0_xyyzzz_1,   \
                             g_0_xx_0_xyyzzzz_0,  \
                             g_0_xx_0_xyyzzzz_1,  \
                             g_0_xx_0_xyzzzz_1,   \
                             g_0_xx_0_xyzzzzz_0,  \
                             g_0_xx_0_xyzzzzz_1,  \
                             g_0_xx_0_xzzzzz_1,   \
                             g_0_xx_0_xzzzzzz_0,  \
                             g_0_xx_0_xzzzzzz_1,  \
                             g_0_xx_0_yyyyyy_1,   \
                             g_0_xx_0_yyyyyyy_0,  \
                             g_0_xx_0_yyyyyyy_1,  \
                             g_0_xx_0_yyyyyyz_0,  \
                             g_0_xx_0_yyyyyyz_1,  \
                             g_0_xx_0_yyyyyz_1,   \
                             g_0_xx_0_yyyyyzz_0,  \
                             g_0_xx_0_yyyyyzz_1,  \
                             g_0_xx_0_yyyyzz_1,   \
                             g_0_xx_0_yyyyzzz_0,  \
                             g_0_xx_0_yyyyzzz_1,  \
                             g_0_xx_0_yyyzzz_1,   \
                             g_0_xx_0_yyyzzzz_0,  \
                             g_0_xx_0_yyyzzzz_1,  \
                             g_0_xx_0_yyzzzz_1,   \
                             g_0_xx_0_yyzzzzz_0,  \
                             g_0_xx_0_yyzzzzz_1,  \
                             g_0_xx_0_yzzzzz_1,   \
                             g_0_xx_0_yzzzzzz_0,  \
                             g_0_xx_0_yzzzzzz_1,  \
                             g_0_xx_0_zzzzzz_1,   \
                             g_0_xx_0_zzzzzzz_0,  \
                             g_0_xx_0_zzzzzzz_1,  \
                             g_0_xxz_0_xxxxxxx_0, \
                             g_0_xxz_0_xxxxxxy_0, \
                             g_0_xxz_0_xxxxxxz_0, \
                             g_0_xxz_0_xxxxxyy_0, \
                             g_0_xxz_0_xxxxxyz_0, \
                             g_0_xxz_0_xxxxxzz_0, \
                             g_0_xxz_0_xxxxyyy_0, \
                             g_0_xxz_0_xxxxyyz_0, \
                             g_0_xxz_0_xxxxyzz_0, \
                             g_0_xxz_0_xxxxzzz_0, \
                             g_0_xxz_0_xxxyyyy_0, \
                             g_0_xxz_0_xxxyyyz_0, \
                             g_0_xxz_0_xxxyyzz_0, \
                             g_0_xxz_0_xxxyzzz_0, \
                             g_0_xxz_0_xxxzzzz_0, \
                             g_0_xxz_0_xxyyyyy_0, \
                             g_0_xxz_0_xxyyyyz_0, \
                             g_0_xxz_0_xxyyyzz_0, \
                             g_0_xxz_0_xxyyzzz_0, \
                             g_0_xxz_0_xxyzzzz_0, \
                             g_0_xxz_0_xxzzzzz_0, \
                             g_0_xxz_0_xyyyyyy_0, \
                             g_0_xxz_0_xyyyyyz_0, \
                             g_0_xxz_0_xyyyyzz_0, \
                             g_0_xxz_0_xyyyzzz_0, \
                             g_0_xxz_0_xyyzzzz_0, \
                             g_0_xxz_0_xyzzzzz_0, \
                             g_0_xxz_0_xzzzzzz_0, \
                             g_0_xxz_0_yyyyyyy_0, \
                             g_0_xxz_0_yyyyyyz_0, \
                             g_0_xxz_0_yyyyyzz_0, \
                             g_0_xxz_0_yyyyzzz_0, \
                             g_0_xxz_0_yyyzzzz_0, \
                             g_0_xxz_0_yyzzzzz_0, \
                             g_0_xxz_0_yzzzzzz_0, \
                             g_0_xxz_0_zzzzzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xxxxxxx_0[i] = g_0_xx_0_xxxxxxx_0[i] * pb_z + g_0_xx_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxxy_0[i] = g_0_xx_0_xxxxxxy_0[i] * pb_z + g_0_xx_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxxz_0[i] = g_0_xx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxxz_0[i] * pb_z + g_0_xx_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxyy_0[i] = g_0_xx_0_xxxxxyy_0[i] * pb_z + g_0_xx_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxyz_0[i] = g_0_xx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxyz_0[i] * pb_z + g_0_xx_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxzz_0[i] = 2.0 * g_0_xx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxzz_0[i] * pb_z + g_0_xx_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyyy_0[i] = g_0_xx_0_xxxxyyy_0[i] * pb_z + g_0_xx_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyyz_0[i] = g_0_xx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyyz_0[i] * pb_z + g_0_xx_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyzz_0[i] = 2.0 * g_0_xx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyzz_0[i] * pb_z + g_0_xx_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxzzz_0[i] = 3.0 * g_0_xx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxzzz_0[i] * pb_z + g_0_xx_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyyy_0[i] = g_0_xx_0_xxxyyyy_0[i] * pb_z + g_0_xx_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyyz_0[i] = g_0_xx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyyz_0[i] * pb_z + g_0_xx_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyzz_0[i] = 2.0 * g_0_xx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyzz_0[i] * pb_z + g_0_xx_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyzzz_0[i] = 3.0 * g_0_xx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzzz_0[i] * pb_z + g_0_xx_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxzzzz_0[i] = 4.0 * g_0_xx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxzzzz_0[i] * pb_z + g_0_xx_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyyy_0[i] = g_0_xx_0_xxyyyyy_0[i] * pb_z + g_0_xx_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyyz_0[i] = g_0_xx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyyz_0[i] * pb_z + g_0_xx_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyzz_0[i] = 2.0 * g_0_xx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyzz_0[i] * pb_z + g_0_xx_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyzzz_0[i] = 3.0 * g_0_xx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzzz_0[i] * pb_z + g_0_xx_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyzzzz_0[i] = 4.0 * g_0_xx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzzz_0[i] * pb_z + g_0_xx_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxzzzzz_0[i] = 5.0 * g_0_xx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzzzzz_0[i] * pb_z + g_0_xx_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyyy_0[i] = g_0_xx_0_xyyyyyy_0[i] * pb_z + g_0_xx_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyyz_0[i] = g_0_xx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyyz_0[i] * pb_z + g_0_xx_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyzz_0[i] = 2.0 * g_0_xx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyzz_0[i] * pb_z + g_0_xx_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyzzz_0[i] = 3.0 * g_0_xx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyzzz_0[i] * pb_z + g_0_xx_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyzzzz_0[i] = 4.0 * g_0_xx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzzzz_0[i] * pb_z + g_0_xx_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyzzzzz_0[i] = 5.0 * g_0_xx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzzzz_0[i] * pb_z + g_0_xx_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xzzzzzz_0[i] = 6.0 * g_0_xx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xzzzzzz_0[i] * pb_z + g_0_xx_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyyy_0[i] = g_0_xx_0_yyyyyyy_0[i] * pb_z + g_0_xx_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyyz_0[i] = g_0_xx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyyz_0[i] * pb_z + g_0_xx_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyzz_0[i] = 2.0 * g_0_xx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyzz_0[i] * pb_z + g_0_xx_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyzzz_0[i] = 3.0 * g_0_xx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyzzz_0[i] * pb_z + g_0_xx_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyzzzz_0[i] = 4.0 * g_0_xx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzzzz_0[i] * pb_z + g_0_xx_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyzzzzz_0[i] = 5.0 * g_0_xx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzzzz_0[i] * pb_z + g_0_xx_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yzzzzzz_0[i] = 6.0 * g_0_xx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzzzz_0[i] * pb_z + g_0_xx_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_zzzzzzz_0[i] = 7.0 * g_0_xx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xx_0_zzzzzzz_0[i] * pb_z + g_0_xx_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 108-144 components of targeted buffer : SFSK

    auto g_0_xyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 108);

    auto g_0_xyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 109);

    auto g_0_xyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 110);

    auto g_0_xyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 111);

    auto g_0_xyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 112);

    auto g_0_xyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 113);

    auto g_0_xyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 114);

    auto g_0_xyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 115);

    auto g_0_xyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 116);

    auto g_0_xyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 117);

    auto g_0_xyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 118);

    auto g_0_xyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 119);

    auto g_0_xyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 120);

    auto g_0_xyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 121);

    auto g_0_xyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 122);

    auto g_0_xyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 123);

    auto g_0_xyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 124);

    auto g_0_xyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 125);

    auto g_0_xyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 126);

    auto g_0_xyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 127);

    auto g_0_xyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 128);

    auto g_0_xyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 129);

    auto g_0_xyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 130);

    auto g_0_xyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 131);

    auto g_0_xyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 132);

    auto g_0_xyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 133);

    auto g_0_xyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 134);

    auto g_0_xyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 135);

    auto g_0_xyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 136);

    auto g_0_xyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 137);

    auto g_0_xyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 138);

    auto g_0_xyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 139);

    auto g_0_xyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 140);

    auto g_0_xyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 141);

    auto g_0_xyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 142);

    auto g_0_xyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 143);

#pragma omp simd aligned(g_0_xyy_0_xxxxxxx_0,     \
                             g_0_xyy_0_xxxxxxy_0, \
                             g_0_xyy_0_xxxxxxz_0, \
                             g_0_xyy_0_xxxxxyy_0, \
                             g_0_xyy_0_xxxxxyz_0, \
                             g_0_xyy_0_xxxxxzz_0, \
                             g_0_xyy_0_xxxxyyy_0, \
                             g_0_xyy_0_xxxxyyz_0, \
                             g_0_xyy_0_xxxxyzz_0, \
                             g_0_xyy_0_xxxxzzz_0, \
                             g_0_xyy_0_xxxyyyy_0, \
                             g_0_xyy_0_xxxyyyz_0, \
                             g_0_xyy_0_xxxyyzz_0, \
                             g_0_xyy_0_xxxyzzz_0, \
                             g_0_xyy_0_xxxzzzz_0, \
                             g_0_xyy_0_xxyyyyy_0, \
                             g_0_xyy_0_xxyyyyz_0, \
                             g_0_xyy_0_xxyyyzz_0, \
                             g_0_xyy_0_xxyyzzz_0, \
                             g_0_xyy_0_xxyzzzz_0, \
                             g_0_xyy_0_xxzzzzz_0, \
                             g_0_xyy_0_xyyyyyy_0, \
                             g_0_xyy_0_xyyyyyz_0, \
                             g_0_xyy_0_xyyyyzz_0, \
                             g_0_xyy_0_xyyyzzz_0, \
                             g_0_xyy_0_xyyzzzz_0, \
                             g_0_xyy_0_xyzzzzz_0, \
                             g_0_xyy_0_xzzzzzz_0, \
                             g_0_xyy_0_yyyyyyy_0, \
                             g_0_xyy_0_yyyyyyz_0, \
                             g_0_xyy_0_yyyyyzz_0, \
                             g_0_xyy_0_yyyyzzz_0, \
                             g_0_xyy_0_yyyzzzz_0, \
                             g_0_xyy_0_yyzzzzz_0, \
                             g_0_xyy_0_yzzzzzz_0, \
                             g_0_xyy_0_zzzzzzz_0, \
                             g_0_yy_0_xxxxxx_1,   \
                             g_0_yy_0_xxxxxxx_0,  \
                             g_0_yy_0_xxxxxxx_1,  \
                             g_0_yy_0_xxxxxxy_0,  \
                             g_0_yy_0_xxxxxxy_1,  \
                             g_0_yy_0_xxxxxxz_0,  \
                             g_0_yy_0_xxxxxxz_1,  \
                             g_0_yy_0_xxxxxy_1,   \
                             g_0_yy_0_xxxxxyy_0,  \
                             g_0_yy_0_xxxxxyy_1,  \
                             g_0_yy_0_xxxxxyz_0,  \
                             g_0_yy_0_xxxxxyz_1,  \
                             g_0_yy_0_xxxxxz_1,   \
                             g_0_yy_0_xxxxxzz_0,  \
                             g_0_yy_0_xxxxxzz_1,  \
                             g_0_yy_0_xxxxyy_1,   \
                             g_0_yy_0_xxxxyyy_0,  \
                             g_0_yy_0_xxxxyyy_1,  \
                             g_0_yy_0_xxxxyyz_0,  \
                             g_0_yy_0_xxxxyyz_1,  \
                             g_0_yy_0_xxxxyz_1,   \
                             g_0_yy_0_xxxxyzz_0,  \
                             g_0_yy_0_xxxxyzz_1,  \
                             g_0_yy_0_xxxxzz_1,   \
                             g_0_yy_0_xxxxzzz_0,  \
                             g_0_yy_0_xxxxzzz_1,  \
                             g_0_yy_0_xxxyyy_1,   \
                             g_0_yy_0_xxxyyyy_0,  \
                             g_0_yy_0_xxxyyyy_1,  \
                             g_0_yy_0_xxxyyyz_0,  \
                             g_0_yy_0_xxxyyyz_1,  \
                             g_0_yy_0_xxxyyz_1,   \
                             g_0_yy_0_xxxyyzz_0,  \
                             g_0_yy_0_xxxyyzz_1,  \
                             g_0_yy_0_xxxyzz_1,   \
                             g_0_yy_0_xxxyzzz_0,  \
                             g_0_yy_0_xxxyzzz_1,  \
                             g_0_yy_0_xxxzzz_1,   \
                             g_0_yy_0_xxxzzzz_0,  \
                             g_0_yy_0_xxxzzzz_1,  \
                             g_0_yy_0_xxyyyy_1,   \
                             g_0_yy_0_xxyyyyy_0,  \
                             g_0_yy_0_xxyyyyy_1,  \
                             g_0_yy_0_xxyyyyz_0,  \
                             g_0_yy_0_xxyyyyz_1,  \
                             g_0_yy_0_xxyyyz_1,   \
                             g_0_yy_0_xxyyyzz_0,  \
                             g_0_yy_0_xxyyyzz_1,  \
                             g_0_yy_0_xxyyzz_1,   \
                             g_0_yy_0_xxyyzzz_0,  \
                             g_0_yy_0_xxyyzzz_1,  \
                             g_0_yy_0_xxyzzz_1,   \
                             g_0_yy_0_xxyzzzz_0,  \
                             g_0_yy_0_xxyzzzz_1,  \
                             g_0_yy_0_xxzzzz_1,   \
                             g_0_yy_0_xxzzzzz_0,  \
                             g_0_yy_0_xxzzzzz_1,  \
                             g_0_yy_0_xyyyyy_1,   \
                             g_0_yy_0_xyyyyyy_0,  \
                             g_0_yy_0_xyyyyyy_1,  \
                             g_0_yy_0_xyyyyyz_0,  \
                             g_0_yy_0_xyyyyyz_1,  \
                             g_0_yy_0_xyyyyz_1,   \
                             g_0_yy_0_xyyyyzz_0,  \
                             g_0_yy_0_xyyyyzz_1,  \
                             g_0_yy_0_xyyyzz_1,   \
                             g_0_yy_0_xyyyzzz_0,  \
                             g_0_yy_0_xyyyzzz_1,  \
                             g_0_yy_0_xyyzzz_1,   \
                             g_0_yy_0_xyyzzzz_0,  \
                             g_0_yy_0_xyyzzzz_1,  \
                             g_0_yy_0_xyzzzz_1,   \
                             g_0_yy_0_xyzzzzz_0,  \
                             g_0_yy_0_xyzzzzz_1,  \
                             g_0_yy_0_xzzzzz_1,   \
                             g_0_yy_0_xzzzzzz_0,  \
                             g_0_yy_0_xzzzzzz_1,  \
                             g_0_yy_0_yyyyyy_1,   \
                             g_0_yy_0_yyyyyyy_0,  \
                             g_0_yy_0_yyyyyyy_1,  \
                             g_0_yy_0_yyyyyyz_0,  \
                             g_0_yy_0_yyyyyyz_1,  \
                             g_0_yy_0_yyyyyz_1,   \
                             g_0_yy_0_yyyyyzz_0,  \
                             g_0_yy_0_yyyyyzz_1,  \
                             g_0_yy_0_yyyyzz_1,   \
                             g_0_yy_0_yyyyzzz_0,  \
                             g_0_yy_0_yyyyzzz_1,  \
                             g_0_yy_0_yyyzzz_1,   \
                             g_0_yy_0_yyyzzzz_0,  \
                             g_0_yy_0_yyyzzzz_1,  \
                             g_0_yy_0_yyzzzz_1,   \
                             g_0_yy_0_yyzzzzz_0,  \
                             g_0_yy_0_yyzzzzz_1,  \
                             g_0_yy_0_yzzzzz_1,   \
                             g_0_yy_0_yzzzzzz_0,  \
                             g_0_yy_0_yzzzzzz_1,  \
                             g_0_yy_0_zzzzzz_1,   \
                             g_0_yy_0_zzzzzzz_0,  \
                             g_0_yy_0_zzzzzzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xxxxxxx_0[i] = 7.0 * g_0_yy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxx_0[i] * pb_x + g_0_yy_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxxy_0[i] = 6.0 * g_0_yy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxy_0[i] * pb_x + g_0_yy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxxz_0[i] = 6.0 * g_0_yy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxz_0[i] * pb_x + g_0_yy_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxyy_0[i] = 5.0 * g_0_yy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyy_0[i] * pb_x + g_0_yy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxyz_0[i] = 5.0 * g_0_yy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyz_0[i] * pb_x + g_0_yy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxzz_0[i] = 5.0 * g_0_yy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxzz_0[i] * pb_x + g_0_yy_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyyy_0[i] = 4.0 * g_0_yy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyy_0[i] * pb_x + g_0_yy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyyz_0[i] = 4.0 * g_0_yy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyz_0[i] * pb_x + g_0_yy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyzz_0[i] = 4.0 * g_0_yy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyzz_0[i] * pb_x + g_0_yy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxzzz_0[i] = 4.0 * g_0_yy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxzzz_0[i] * pb_x + g_0_yy_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyyy_0[i] = 3.0 * g_0_yy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyy_0[i] * pb_x + g_0_yy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyyz_0[i] = 3.0 * g_0_yy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyz_0[i] * pb_x + g_0_yy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyzz_0[i] = 3.0 * g_0_yy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyzz_0[i] * pb_x + g_0_yy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyzzz_0[i] = 3.0 * g_0_yy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyzzz_0[i] * pb_x + g_0_yy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxzzzz_0[i] = 3.0 * g_0_yy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzzzz_0[i] * pb_x + g_0_yy_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyyy_0[i] = 2.0 * g_0_yy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyy_0[i] * pb_x + g_0_yy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyyz_0[i] = 2.0 * g_0_yy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyz_0[i] * pb_x + g_0_yy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyzz_0[i] = 2.0 * g_0_yy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyzz_0[i] * pb_x + g_0_yy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyzzz_0[i] = 2.0 * g_0_yy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzzz_0[i] * pb_x + g_0_yy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyzzzz_0[i] = 2.0 * g_0_yy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzzzz_0[i] * pb_x + g_0_yy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxzzzzz_0[i] = 2.0 * g_0_yy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzzzz_0[i] * pb_x + g_0_yy_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyyy_0[i] = g_0_yy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyy_0[i] * pb_x + g_0_yy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyyz_0[i] = g_0_yy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyz_0[i] * pb_x + g_0_yy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyzz_0[i] = g_0_yy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyzz_0[i] * pb_x + g_0_yy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyzzz_0[i] = g_0_yy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzzz_0[i] * pb_x + g_0_yy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyzzzz_0[i] = g_0_yy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzzz_0[i] * pb_x + g_0_yy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyzzzzz_0[i] = g_0_yy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzzzz_0[i] * pb_x + g_0_yy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xzzzzzz_0[i] = g_0_yy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzzzz_0[i] * pb_x + g_0_yy_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyyy_0[i] = g_0_yy_0_yyyyyyy_0[i] * pb_x + g_0_yy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyyz_0[i] = g_0_yy_0_yyyyyyz_0[i] * pb_x + g_0_yy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyzz_0[i] = g_0_yy_0_yyyyyzz_0[i] * pb_x + g_0_yy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyzzz_0[i] = g_0_yy_0_yyyyzzz_0[i] * pb_x + g_0_yy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyzzzz_0[i] = g_0_yy_0_yyyzzzz_0[i] * pb_x + g_0_yy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyzzzzz_0[i] = g_0_yy_0_yyzzzzz_0[i] * pb_x + g_0_yy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yzzzzzz_0[i] = g_0_yy_0_yzzzzzz_0[i] * pb_x + g_0_yy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_zzzzzzz_0[i] = g_0_yy_0_zzzzzzz_0[i] * pb_x + g_0_yy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 144-180 components of targeted buffer : SFSK

    auto g_0_xyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 144);

    auto g_0_xyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 145);

    auto g_0_xyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 146);

    auto g_0_xyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 147);

    auto g_0_xyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 148);

    auto g_0_xyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 149);

    auto g_0_xyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 150);

    auto g_0_xyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 151);

    auto g_0_xyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 152);

    auto g_0_xyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 153);

    auto g_0_xyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 154);

    auto g_0_xyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 155);

    auto g_0_xyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 156);

    auto g_0_xyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 157);

    auto g_0_xyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 158);

    auto g_0_xyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 159);

    auto g_0_xyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 160);

    auto g_0_xyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 161);

    auto g_0_xyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 162);

    auto g_0_xyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 163);

    auto g_0_xyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 164);

    auto g_0_xyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 165);

    auto g_0_xyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 166);

    auto g_0_xyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 167);

    auto g_0_xyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 168);

    auto g_0_xyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 169);

    auto g_0_xyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 170);

    auto g_0_xyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 171);

    auto g_0_xyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 172);

    auto g_0_xyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 173);

    auto g_0_xyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 174);

    auto g_0_xyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 175);

    auto g_0_xyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 176);

    auto g_0_xyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 177);

    auto g_0_xyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 178);

    auto g_0_xyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 179);

#pragma omp simd aligned(g_0_xy_0_xxxxxxy_0,      \
                             g_0_xy_0_xxxxxxy_1,  \
                             g_0_xy_0_xxxxxyy_0,  \
                             g_0_xy_0_xxxxxyy_1,  \
                             g_0_xy_0_xxxxyyy_0,  \
                             g_0_xy_0_xxxxyyy_1,  \
                             g_0_xy_0_xxxyyyy_0,  \
                             g_0_xy_0_xxxyyyy_1,  \
                             g_0_xy_0_xxyyyyy_0,  \
                             g_0_xy_0_xxyyyyy_1,  \
                             g_0_xy_0_xyyyyyy_0,  \
                             g_0_xy_0_xyyyyyy_1,  \
                             g_0_xyz_0_xxxxxxx_0, \
                             g_0_xyz_0_xxxxxxy_0, \
                             g_0_xyz_0_xxxxxxz_0, \
                             g_0_xyz_0_xxxxxyy_0, \
                             g_0_xyz_0_xxxxxyz_0, \
                             g_0_xyz_0_xxxxxzz_0, \
                             g_0_xyz_0_xxxxyyy_0, \
                             g_0_xyz_0_xxxxyyz_0, \
                             g_0_xyz_0_xxxxyzz_0, \
                             g_0_xyz_0_xxxxzzz_0, \
                             g_0_xyz_0_xxxyyyy_0, \
                             g_0_xyz_0_xxxyyyz_0, \
                             g_0_xyz_0_xxxyyzz_0, \
                             g_0_xyz_0_xxxyzzz_0, \
                             g_0_xyz_0_xxxzzzz_0, \
                             g_0_xyz_0_xxyyyyy_0, \
                             g_0_xyz_0_xxyyyyz_0, \
                             g_0_xyz_0_xxyyyzz_0, \
                             g_0_xyz_0_xxyyzzz_0, \
                             g_0_xyz_0_xxyzzzz_0, \
                             g_0_xyz_0_xxzzzzz_0, \
                             g_0_xyz_0_xyyyyyy_0, \
                             g_0_xyz_0_xyyyyyz_0, \
                             g_0_xyz_0_xyyyyzz_0, \
                             g_0_xyz_0_xyyyzzz_0, \
                             g_0_xyz_0_xyyzzzz_0, \
                             g_0_xyz_0_xyzzzzz_0, \
                             g_0_xyz_0_xzzzzzz_0, \
                             g_0_xyz_0_yyyyyyy_0, \
                             g_0_xyz_0_yyyyyyz_0, \
                             g_0_xyz_0_yyyyyzz_0, \
                             g_0_xyz_0_yyyyzzz_0, \
                             g_0_xyz_0_yyyzzzz_0, \
                             g_0_xyz_0_yyzzzzz_0, \
                             g_0_xyz_0_yzzzzzz_0, \
                             g_0_xyz_0_zzzzzzz_0, \
                             g_0_xz_0_xxxxxxx_0,  \
                             g_0_xz_0_xxxxxxx_1,  \
                             g_0_xz_0_xxxxxxz_0,  \
                             g_0_xz_0_xxxxxxz_1,  \
                             g_0_xz_0_xxxxxzz_0,  \
                             g_0_xz_0_xxxxxzz_1,  \
                             g_0_xz_0_xxxxzzz_0,  \
                             g_0_xz_0_xxxxzzz_1,  \
                             g_0_xz_0_xxxzzzz_0,  \
                             g_0_xz_0_xxxzzzz_1,  \
                             g_0_xz_0_xxzzzzz_0,  \
                             g_0_xz_0_xxzzzzz_1,  \
                             g_0_xz_0_xzzzzzz_0,  \
                             g_0_xz_0_xzzzzzz_1,  \
                             g_0_yz_0_xxxxxyz_0,  \
                             g_0_yz_0_xxxxxyz_1,  \
                             g_0_yz_0_xxxxyyz_0,  \
                             g_0_yz_0_xxxxyyz_1,  \
                             g_0_yz_0_xxxxyz_1,   \
                             g_0_yz_0_xxxxyzz_0,  \
                             g_0_yz_0_xxxxyzz_1,  \
                             g_0_yz_0_xxxyyyz_0,  \
                             g_0_yz_0_xxxyyyz_1,  \
                             g_0_yz_0_xxxyyz_1,   \
                             g_0_yz_0_xxxyyzz_0,  \
                             g_0_yz_0_xxxyyzz_1,  \
                             g_0_yz_0_xxxyzz_1,   \
                             g_0_yz_0_xxxyzzz_0,  \
                             g_0_yz_0_xxxyzzz_1,  \
                             g_0_yz_0_xxyyyyz_0,  \
                             g_0_yz_0_xxyyyyz_1,  \
                             g_0_yz_0_xxyyyz_1,   \
                             g_0_yz_0_xxyyyzz_0,  \
                             g_0_yz_0_xxyyyzz_1,  \
                             g_0_yz_0_xxyyzz_1,   \
                             g_0_yz_0_xxyyzzz_0,  \
                             g_0_yz_0_xxyyzzz_1,  \
                             g_0_yz_0_xxyzzz_1,   \
                             g_0_yz_0_xxyzzzz_0,  \
                             g_0_yz_0_xxyzzzz_1,  \
                             g_0_yz_0_xyyyyyz_0,  \
                             g_0_yz_0_xyyyyyz_1,  \
                             g_0_yz_0_xyyyyz_1,   \
                             g_0_yz_0_xyyyyzz_0,  \
                             g_0_yz_0_xyyyyzz_1,  \
                             g_0_yz_0_xyyyzz_1,   \
                             g_0_yz_0_xyyyzzz_0,  \
                             g_0_yz_0_xyyyzzz_1,  \
                             g_0_yz_0_xyyzzz_1,   \
                             g_0_yz_0_xyyzzzz_0,  \
                             g_0_yz_0_xyyzzzz_1,  \
                             g_0_yz_0_xyzzzz_1,   \
                             g_0_yz_0_xyzzzzz_0,  \
                             g_0_yz_0_xyzzzzz_1,  \
                             g_0_yz_0_yyyyyyy_0,  \
                             g_0_yz_0_yyyyyyy_1,  \
                             g_0_yz_0_yyyyyyz_0,  \
                             g_0_yz_0_yyyyyyz_1,  \
                             g_0_yz_0_yyyyyz_1,   \
                             g_0_yz_0_yyyyyzz_0,  \
                             g_0_yz_0_yyyyyzz_1,  \
                             g_0_yz_0_yyyyzz_1,   \
                             g_0_yz_0_yyyyzzz_0,  \
                             g_0_yz_0_yyyyzzz_1,  \
                             g_0_yz_0_yyyzzz_1,   \
                             g_0_yz_0_yyyzzzz_0,  \
                             g_0_yz_0_yyyzzzz_1,  \
                             g_0_yz_0_yyzzzz_1,   \
                             g_0_yz_0_yyzzzzz_0,  \
                             g_0_yz_0_yyzzzzz_1,  \
                             g_0_yz_0_yzzzzz_1,   \
                             g_0_yz_0_yzzzzzz_0,  \
                             g_0_yz_0_yzzzzzz_1,  \
                             g_0_yz_0_zzzzzzz_0,  \
                             g_0_yz_0_zzzzzzz_1,  \
                             wp_x,                \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyz_0_xxxxxxx_0[i] = g_0_xz_0_xxxxxxx_0[i] * pb_y + g_0_xz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xyz_0_xxxxxxy_0[i] = g_0_xy_0_xxxxxxy_0[i] * pb_z + g_0_xy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxxxz_0[i] = g_0_xz_0_xxxxxxz_0[i] * pb_y + g_0_xz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xyz_0_xxxxxyy_0[i] = g_0_xy_0_xxxxxyy_0[i] * pb_z + g_0_xy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxxyz_0[i] = 5.0 * g_0_yz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxxyz_0[i] * pb_x + g_0_yz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxxzz_0[i] = g_0_xz_0_xxxxxzz_0[i] * pb_y + g_0_xz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xyz_0_xxxxyyy_0[i] = g_0_xy_0_xxxxyyy_0[i] * pb_z + g_0_xy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxyyz_0[i] = 4.0 * g_0_yz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxyyz_0[i] * pb_x + g_0_yz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxyzz_0[i] = 4.0 * g_0_yz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxyzz_0[i] * pb_x + g_0_yz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxzzz_0[i] = g_0_xz_0_xxxxzzz_0[i] * pb_y + g_0_xz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xyz_0_xxxyyyy_0[i] = g_0_xy_0_xxxyyyy_0[i] * pb_z + g_0_xy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxyyyz_0[i] = 3.0 * g_0_yz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyyyz_0[i] * pb_x + g_0_yz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxyyzz_0[i] = 3.0 * g_0_yz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyyzz_0[i] * pb_x + g_0_yz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxyzzz_0[i] = 3.0 * g_0_yz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyzzz_0[i] * pb_x + g_0_yz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxzzzz_0[i] = g_0_xz_0_xxxzzzz_0[i] * pb_y + g_0_xz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xyz_0_xxyyyyy_0[i] = g_0_xy_0_xxyyyyy_0[i] * pb_z + g_0_xy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxyyyyz_0[i] = 2.0 * g_0_yz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyyyz_0[i] * pb_x + g_0_yz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxyyyzz_0[i] = 2.0 * g_0_yz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyyzz_0[i] * pb_x + g_0_yz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxyyzzz_0[i] = 2.0 * g_0_yz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyzzz_0[i] * pb_x + g_0_yz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxyzzzz_0[i] = 2.0 * g_0_yz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyzzzz_0[i] * pb_x + g_0_yz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxzzzzz_0[i] = g_0_xz_0_xxzzzzz_0[i] * pb_y + g_0_xz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xyz_0_xyyyyyy_0[i] = g_0_xy_0_xyyyyyy_0[i] * pb_z + g_0_xy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xyyyyyz_0[i] = g_0_yz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyyyz_0[i] * pb_x + g_0_yz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xyyyyzz_0[i] = g_0_yz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyyzz_0[i] * pb_x + g_0_yz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xyyyzzz_0[i] = g_0_yz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyzzz_0[i] * pb_x + g_0_yz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xyyzzzz_0[i] = g_0_yz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyzzzz_0[i] * pb_x + g_0_yz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xyzzzzz_0[i] = g_0_yz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyzzzzz_0[i] * pb_x + g_0_yz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xzzzzzz_0[i] = g_0_xz_0_xzzzzzz_0[i] * pb_y + g_0_xz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xyz_0_yyyyyyy_0[i] = g_0_yz_0_yyyyyyy_0[i] * pb_x + g_0_yz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyz_0_yyyyyyz_0[i] = g_0_yz_0_yyyyyyz_0[i] * pb_x + g_0_yz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_yyyyyzz_0[i] = g_0_yz_0_yyyyyzz_0[i] * pb_x + g_0_yz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_yyyyzzz_0[i] = g_0_yz_0_yyyyzzz_0[i] * pb_x + g_0_yz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_yyyzzzz_0[i] = g_0_yz_0_yyyzzzz_0[i] * pb_x + g_0_yz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_yyzzzzz_0[i] = g_0_yz_0_yyzzzzz_0[i] * pb_x + g_0_yz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_yzzzzzz_0[i] = g_0_yz_0_yzzzzzz_0[i] * pb_x + g_0_yz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_zzzzzzz_0[i] = g_0_yz_0_zzzzzzz_0[i] * pb_x + g_0_yz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 180-216 components of targeted buffer : SFSK

    auto g_0_xzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 180);

    auto g_0_xzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 181);

    auto g_0_xzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 182);

    auto g_0_xzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 183);

    auto g_0_xzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 184);

    auto g_0_xzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 185);

    auto g_0_xzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 186);

    auto g_0_xzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 187);

    auto g_0_xzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 188);

    auto g_0_xzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 189);

    auto g_0_xzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 190);

    auto g_0_xzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 191);

    auto g_0_xzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 192);

    auto g_0_xzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 193);

    auto g_0_xzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 194);

    auto g_0_xzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 195);

    auto g_0_xzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 196);

    auto g_0_xzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 197);

    auto g_0_xzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 198);

    auto g_0_xzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 199);

    auto g_0_xzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 200);

    auto g_0_xzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 201);

    auto g_0_xzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 202);

    auto g_0_xzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 203);

    auto g_0_xzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 204);

    auto g_0_xzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 205);

    auto g_0_xzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 206);

    auto g_0_xzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 207);

    auto g_0_xzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 208);

    auto g_0_xzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 209);

    auto g_0_xzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 210);

    auto g_0_xzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 211);

    auto g_0_xzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 212);

    auto g_0_xzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 213);

    auto g_0_xzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 214);

    auto g_0_xzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 215);

#pragma omp simd aligned(g_0_xzz_0_xxxxxxx_0,     \
                             g_0_xzz_0_xxxxxxy_0, \
                             g_0_xzz_0_xxxxxxz_0, \
                             g_0_xzz_0_xxxxxyy_0, \
                             g_0_xzz_0_xxxxxyz_0, \
                             g_0_xzz_0_xxxxxzz_0, \
                             g_0_xzz_0_xxxxyyy_0, \
                             g_0_xzz_0_xxxxyyz_0, \
                             g_0_xzz_0_xxxxyzz_0, \
                             g_0_xzz_0_xxxxzzz_0, \
                             g_0_xzz_0_xxxyyyy_0, \
                             g_0_xzz_0_xxxyyyz_0, \
                             g_0_xzz_0_xxxyyzz_0, \
                             g_0_xzz_0_xxxyzzz_0, \
                             g_0_xzz_0_xxxzzzz_0, \
                             g_0_xzz_0_xxyyyyy_0, \
                             g_0_xzz_0_xxyyyyz_0, \
                             g_0_xzz_0_xxyyyzz_0, \
                             g_0_xzz_0_xxyyzzz_0, \
                             g_0_xzz_0_xxyzzzz_0, \
                             g_0_xzz_0_xxzzzzz_0, \
                             g_0_xzz_0_xyyyyyy_0, \
                             g_0_xzz_0_xyyyyyz_0, \
                             g_0_xzz_0_xyyyyzz_0, \
                             g_0_xzz_0_xyyyzzz_0, \
                             g_0_xzz_0_xyyzzzz_0, \
                             g_0_xzz_0_xyzzzzz_0, \
                             g_0_xzz_0_xzzzzzz_0, \
                             g_0_xzz_0_yyyyyyy_0, \
                             g_0_xzz_0_yyyyyyz_0, \
                             g_0_xzz_0_yyyyyzz_0, \
                             g_0_xzz_0_yyyyzzz_0, \
                             g_0_xzz_0_yyyzzzz_0, \
                             g_0_xzz_0_yyzzzzz_0, \
                             g_0_xzz_0_yzzzzzz_0, \
                             g_0_xzz_0_zzzzzzz_0, \
                             g_0_zz_0_xxxxxx_1,   \
                             g_0_zz_0_xxxxxxx_0,  \
                             g_0_zz_0_xxxxxxx_1,  \
                             g_0_zz_0_xxxxxxy_0,  \
                             g_0_zz_0_xxxxxxy_1,  \
                             g_0_zz_0_xxxxxxz_0,  \
                             g_0_zz_0_xxxxxxz_1,  \
                             g_0_zz_0_xxxxxy_1,   \
                             g_0_zz_0_xxxxxyy_0,  \
                             g_0_zz_0_xxxxxyy_1,  \
                             g_0_zz_0_xxxxxyz_0,  \
                             g_0_zz_0_xxxxxyz_1,  \
                             g_0_zz_0_xxxxxz_1,   \
                             g_0_zz_0_xxxxxzz_0,  \
                             g_0_zz_0_xxxxxzz_1,  \
                             g_0_zz_0_xxxxyy_1,   \
                             g_0_zz_0_xxxxyyy_0,  \
                             g_0_zz_0_xxxxyyy_1,  \
                             g_0_zz_0_xxxxyyz_0,  \
                             g_0_zz_0_xxxxyyz_1,  \
                             g_0_zz_0_xxxxyz_1,   \
                             g_0_zz_0_xxxxyzz_0,  \
                             g_0_zz_0_xxxxyzz_1,  \
                             g_0_zz_0_xxxxzz_1,   \
                             g_0_zz_0_xxxxzzz_0,  \
                             g_0_zz_0_xxxxzzz_1,  \
                             g_0_zz_0_xxxyyy_1,   \
                             g_0_zz_0_xxxyyyy_0,  \
                             g_0_zz_0_xxxyyyy_1,  \
                             g_0_zz_0_xxxyyyz_0,  \
                             g_0_zz_0_xxxyyyz_1,  \
                             g_0_zz_0_xxxyyz_1,   \
                             g_0_zz_0_xxxyyzz_0,  \
                             g_0_zz_0_xxxyyzz_1,  \
                             g_0_zz_0_xxxyzz_1,   \
                             g_0_zz_0_xxxyzzz_0,  \
                             g_0_zz_0_xxxyzzz_1,  \
                             g_0_zz_0_xxxzzz_1,   \
                             g_0_zz_0_xxxzzzz_0,  \
                             g_0_zz_0_xxxzzzz_1,  \
                             g_0_zz_0_xxyyyy_1,   \
                             g_0_zz_0_xxyyyyy_0,  \
                             g_0_zz_0_xxyyyyy_1,  \
                             g_0_zz_0_xxyyyyz_0,  \
                             g_0_zz_0_xxyyyyz_1,  \
                             g_0_zz_0_xxyyyz_1,   \
                             g_0_zz_0_xxyyyzz_0,  \
                             g_0_zz_0_xxyyyzz_1,  \
                             g_0_zz_0_xxyyzz_1,   \
                             g_0_zz_0_xxyyzzz_0,  \
                             g_0_zz_0_xxyyzzz_1,  \
                             g_0_zz_0_xxyzzz_1,   \
                             g_0_zz_0_xxyzzzz_0,  \
                             g_0_zz_0_xxyzzzz_1,  \
                             g_0_zz_0_xxzzzz_1,   \
                             g_0_zz_0_xxzzzzz_0,  \
                             g_0_zz_0_xxzzzzz_1,  \
                             g_0_zz_0_xyyyyy_1,   \
                             g_0_zz_0_xyyyyyy_0,  \
                             g_0_zz_0_xyyyyyy_1,  \
                             g_0_zz_0_xyyyyyz_0,  \
                             g_0_zz_0_xyyyyyz_1,  \
                             g_0_zz_0_xyyyyz_1,   \
                             g_0_zz_0_xyyyyzz_0,  \
                             g_0_zz_0_xyyyyzz_1,  \
                             g_0_zz_0_xyyyzz_1,   \
                             g_0_zz_0_xyyyzzz_0,  \
                             g_0_zz_0_xyyyzzz_1,  \
                             g_0_zz_0_xyyzzz_1,   \
                             g_0_zz_0_xyyzzzz_0,  \
                             g_0_zz_0_xyyzzzz_1,  \
                             g_0_zz_0_xyzzzz_1,   \
                             g_0_zz_0_xyzzzzz_0,  \
                             g_0_zz_0_xyzzzzz_1,  \
                             g_0_zz_0_xzzzzz_1,   \
                             g_0_zz_0_xzzzzzz_0,  \
                             g_0_zz_0_xzzzzzz_1,  \
                             g_0_zz_0_yyyyyy_1,   \
                             g_0_zz_0_yyyyyyy_0,  \
                             g_0_zz_0_yyyyyyy_1,  \
                             g_0_zz_0_yyyyyyz_0,  \
                             g_0_zz_0_yyyyyyz_1,  \
                             g_0_zz_0_yyyyyz_1,   \
                             g_0_zz_0_yyyyyzz_0,  \
                             g_0_zz_0_yyyyyzz_1,  \
                             g_0_zz_0_yyyyzz_1,   \
                             g_0_zz_0_yyyyzzz_0,  \
                             g_0_zz_0_yyyyzzz_1,  \
                             g_0_zz_0_yyyzzz_1,   \
                             g_0_zz_0_yyyzzzz_0,  \
                             g_0_zz_0_yyyzzzz_1,  \
                             g_0_zz_0_yyzzzz_1,   \
                             g_0_zz_0_yyzzzzz_0,  \
                             g_0_zz_0_yyzzzzz_1,  \
                             g_0_zz_0_yzzzzz_1,   \
                             g_0_zz_0_yzzzzzz_0,  \
                             g_0_zz_0_yzzzzzz_1,  \
                             g_0_zz_0_zzzzzz_1,   \
                             g_0_zz_0_zzzzzzz_0,  \
                             g_0_zz_0_zzzzzzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xxxxxxx_0[i] = 7.0 * g_0_zz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxx_0[i] * pb_x + g_0_zz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxxy_0[i] = 6.0 * g_0_zz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxy_0[i] * pb_x + g_0_zz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxxz_0[i] = 6.0 * g_0_zz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxz_0[i] * pb_x + g_0_zz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxyy_0[i] = 5.0 * g_0_zz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyy_0[i] * pb_x + g_0_zz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxyz_0[i] = 5.0 * g_0_zz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyz_0[i] * pb_x + g_0_zz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxzz_0[i] = 5.0 * g_0_zz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxzz_0[i] * pb_x + g_0_zz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyyy_0[i] = 4.0 * g_0_zz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyy_0[i] * pb_x + g_0_zz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyyz_0[i] = 4.0 * g_0_zz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyz_0[i] * pb_x + g_0_zz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyzz_0[i] = 4.0 * g_0_zz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyzz_0[i] * pb_x + g_0_zz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxzzz_0[i] = 4.0 * g_0_zz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxzzz_0[i] * pb_x + g_0_zz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyyy_0[i] = 3.0 * g_0_zz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyy_0[i] * pb_x + g_0_zz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyyz_0[i] = 3.0 * g_0_zz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyz_0[i] * pb_x + g_0_zz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyzz_0[i] = 3.0 * g_0_zz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyzz_0[i] * pb_x + g_0_zz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyzzz_0[i] = 3.0 * g_0_zz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzzz_0[i] * pb_x + g_0_zz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxzzzz_0[i] = 3.0 * g_0_zz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxzzzz_0[i] * pb_x + g_0_zz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyyy_0[i] = 2.0 * g_0_zz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyy_0[i] * pb_x + g_0_zz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyyz_0[i] = 2.0 * g_0_zz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyz_0[i] * pb_x + g_0_zz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyzz_0[i] = 2.0 * g_0_zz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyzz_0[i] * pb_x + g_0_zz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyzzz_0[i] = 2.0 * g_0_zz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzzz_0[i] * pb_x + g_0_zz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyzzzz_0[i] = 2.0 * g_0_zz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzzz_0[i] * pb_x + g_0_zz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxzzzzz_0[i] = 2.0 * g_0_zz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzzzzz_0[i] * pb_x + g_0_zz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyyy_0[i] = g_0_zz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyy_0[i] * pb_x + g_0_zz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyyz_0[i] = g_0_zz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyz_0[i] * pb_x + g_0_zz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyzz_0[i] = g_0_zz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyzz_0[i] * pb_x + g_0_zz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyzzz_0[i] = g_0_zz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzzz_0[i] * pb_x + g_0_zz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyzzzz_0[i] = g_0_zz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzzz_0[i] * pb_x + g_0_zz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyzzzzz_0[i] = g_0_zz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzzz_0[i] * pb_x + g_0_zz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xzzzzzz_0[i] = g_0_zz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzzzzz_0[i] * pb_x + g_0_zz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyyy_0[i] = g_0_zz_0_yyyyyyy_0[i] * pb_x + g_0_zz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyyz_0[i] = g_0_zz_0_yyyyyyz_0[i] * pb_x + g_0_zz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyzz_0[i] = g_0_zz_0_yyyyyzz_0[i] * pb_x + g_0_zz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyzzz_0[i] = g_0_zz_0_yyyyzzz_0[i] * pb_x + g_0_zz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyzzzz_0[i] = g_0_zz_0_yyyzzzz_0[i] * pb_x + g_0_zz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyzzzzz_0[i] = g_0_zz_0_yyzzzzz_0[i] * pb_x + g_0_zz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yzzzzzz_0[i] = g_0_zz_0_yzzzzzz_0[i] * pb_x + g_0_zz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_zzzzzzz_0[i] = g_0_zz_0_zzzzzzz_0[i] * pb_x + g_0_zz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 216-252 components of targeted buffer : SFSK

    auto g_0_yyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 216);

    auto g_0_yyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 217);

    auto g_0_yyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 218);

    auto g_0_yyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 219);

    auto g_0_yyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 220);

    auto g_0_yyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 221);

    auto g_0_yyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 222);

    auto g_0_yyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 223);

    auto g_0_yyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 224);

    auto g_0_yyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 225);

    auto g_0_yyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 226);

    auto g_0_yyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 227);

    auto g_0_yyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 228);

    auto g_0_yyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 229);

    auto g_0_yyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 230);

    auto g_0_yyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 231);

    auto g_0_yyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 232);

    auto g_0_yyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 233);

    auto g_0_yyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 234);

    auto g_0_yyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 235);

    auto g_0_yyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 236);

    auto g_0_yyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 237);

    auto g_0_yyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 238);

    auto g_0_yyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 239);

    auto g_0_yyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 240);

    auto g_0_yyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 241);

    auto g_0_yyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 242);

    auto g_0_yyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 243);

    auto g_0_yyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 244);

    auto g_0_yyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 245);

    auto g_0_yyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 246);

    auto g_0_yyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 247);

    auto g_0_yyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 248);

    auto g_0_yyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 249);

    auto g_0_yyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 250);

    auto g_0_yyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 251);

#pragma omp simd aligned(g_0_y_0_xxxxxxx_0,       \
                             g_0_y_0_xxxxxxx_1,   \
                             g_0_y_0_xxxxxxy_0,   \
                             g_0_y_0_xxxxxxy_1,   \
                             g_0_y_0_xxxxxxz_0,   \
                             g_0_y_0_xxxxxxz_1,   \
                             g_0_y_0_xxxxxyy_0,   \
                             g_0_y_0_xxxxxyy_1,   \
                             g_0_y_0_xxxxxyz_0,   \
                             g_0_y_0_xxxxxyz_1,   \
                             g_0_y_0_xxxxxzz_0,   \
                             g_0_y_0_xxxxxzz_1,   \
                             g_0_y_0_xxxxyyy_0,   \
                             g_0_y_0_xxxxyyy_1,   \
                             g_0_y_0_xxxxyyz_0,   \
                             g_0_y_0_xxxxyyz_1,   \
                             g_0_y_0_xxxxyzz_0,   \
                             g_0_y_0_xxxxyzz_1,   \
                             g_0_y_0_xxxxzzz_0,   \
                             g_0_y_0_xxxxzzz_1,   \
                             g_0_y_0_xxxyyyy_0,   \
                             g_0_y_0_xxxyyyy_1,   \
                             g_0_y_0_xxxyyyz_0,   \
                             g_0_y_0_xxxyyyz_1,   \
                             g_0_y_0_xxxyyzz_0,   \
                             g_0_y_0_xxxyyzz_1,   \
                             g_0_y_0_xxxyzzz_0,   \
                             g_0_y_0_xxxyzzz_1,   \
                             g_0_y_0_xxxzzzz_0,   \
                             g_0_y_0_xxxzzzz_1,   \
                             g_0_y_0_xxyyyyy_0,   \
                             g_0_y_0_xxyyyyy_1,   \
                             g_0_y_0_xxyyyyz_0,   \
                             g_0_y_0_xxyyyyz_1,   \
                             g_0_y_0_xxyyyzz_0,   \
                             g_0_y_0_xxyyyzz_1,   \
                             g_0_y_0_xxyyzzz_0,   \
                             g_0_y_0_xxyyzzz_1,   \
                             g_0_y_0_xxyzzzz_0,   \
                             g_0_y_0_xxyzzzz_1,   \
                             g_0_y_0_xxzzzzz_0,   \
                             g_0_y_0_xxzzzzz_1,   \
                             g_0_y_0_xyyyyyy_0,   \
                             g_0_y_0_xyyyyyy_1,   \
                             g_0_y_0_xyyyyyz_0,   \
                             g_0_y_0_xyyyyyz_1,   \
                             g_0_y_0_xyyyyzz_0,   \
                             g_0_y_0_xyyyyzz_1,   \
                             g_0_y_0_xyyyzzz_0,   \
                             g_0_y_0_xyyyzzz_1,   \
                             g_0_y_0_xyyzzzz_0,   \
                             g_0_y_0_xyyzzzz_1,   \
                             g_0_y_0_xyzzzzz_0,   \
                             g_0_y_0_xyzzzzz_1,   \
                             g_0_y_0_xzzzzzz_0,   \
                             g_0_y_0_xzzzzzz_1,   \
                             g_0_y_0_yyyyyyy_0,   \
                             g_0_y_0_yyyyyyy_1,   \
                             g_0_y_0_yyyyyyz_0,   \
                             g_0_y_0_yyyyyyz_1,   \
                             g_0_y_0_yyyyyzz_0,   \
                             g_0_y_0_yyyyyzz_1,   \
                             g_0_y_0_yyyyzzz_0,   \
                             g_0_y_0_yyyyzzz_1,   \
                             g_0_y_0_yyyzzzz_0,   \
                             g_0_y_0_yyyzzzz_1,   \
                             g_0_y_0_yyzzzzz_0,   \
                             g_0_y_0_yyzzzzz_1,   \
                             g_0_y_0_yzzzzzz_0,   \
                             g_0_y_0_yzzzzzz_1,   \
                             g_0_y_0_zzzzzzz_0,   \
                             g_0_y_0_zzzzzzz_1,   \
                             g_0_yy_0_xxxxxx_1,   \
                             g_0_yy_0_xxxxxxx_0,  \
                             g_0_yy_0_xxxxxxx_1,  \
                             g_0_yy_0_xxxxxxy_0,  \
                             g_0_yy_0_xxxxxxy_1,  \
                             g_0_yy_0_xxxxxxz_0,  \
                             g_0_yy_0_xxxxxxz_1,  \
                             g_0_yy_0_xxxxxy_1,   \
                             g_0_yy_0_xxxxxyy_0,  \
                             g_0_yy_0_xxxxxyy_1,  \
                             g_0_yy_0_xxxxxyz_0,  \
                             g_0_yy_0_xxxxxyz_1,  \
                             g_0_yy_0_xxxxxz_1,   \
                             g_0_yy_0_xxxxxzz_0,  \
                             g_0_yy_0_xxxxxzz_1,  \
                             g_0_yy_0_xxxxyy_1,   \
                             g_0_yy_0_xxxxyyy_0,  \
                             g_0_yy_0_xxxxyyy_1,  \
                             g_0_yy_0_xxxxyyz_0,  \
                             g_0_yy_0_xxxxyyz_1,  \
                             g_0_yy_0_xxxxyz_1,   \
                             g_0_yy_0_xxxxyzz_0,  \
                             g_0_yy_0_xxxxyzz_1,  \
                             g_0_yy_0_xxxxzz_1,   \
                             g_0_yy_0_xxxxzzz_0,  \
                             g_0_yy_0_xxxxzzz_1,  \
                             g_0_yy_0_xxxyyy_1,   \
                             g_0_yy_0_xxxyyyy_0,  \
                             g_0_yy_0_xxxyyyy_1,  \
                             g_0_yy_0_xxxyyyz_0,  \
                             g_0_yy_0_xxxyyyz_1,  \
                             g_0_yy_0_xxxyyz_1,   \
                             g_0_yy_0_xxxyyzz_0,  \
                             g_0_yy_0_xxxyyzz_1,  \
                             g_0_yy_0_xxxyzz_1,   \
                             g_0_yy_0_xxxyzzz_0,  \
                             g_0_yy_0_xxxyzzz_1,  \
                             g_0_yy_0_xxxzzz_1,   \
                             g_0_yy_0_xxxzzzz_0,  \
                             g_0_yy_0_xxxzzzz_1,  \
                             g_0_yy_0_xxyyyy_1,   \
                             g_0_yy_0_xxyyyyy_0,  \
                             g_0_yy_0_xxyyyyy_1,  \
                             g_0_yy_0_xxyyyyz_0,  \
                             g_0_yy_0_xxyyyyz_1,  \
                             g_0_yy_0_xxyyyz_1,   \
                             g_0_yy_0_xxyyyzz_0,  \
                             g_0_yy_0_xxyyyzz_1,  \
                             g_0_yy_0_xxyyzz_1,   \
                             g_0_yy_0_xxyyzzz_0,  \
                             g_0_yy_0_xxyyzzz_1,  \
                             g_0_yy_0_xxyzzz_1,   \
                             g_0_yy_0_xxyzzzz_0,  \
                             g_0_yy_0_xxyzzzz_1,  \
                             g_0_yy_0_xxzzzz_1,   \
                             g_0_yy_0_xxzzzzz_0,  \
                             g_0_yy_0_xxzzzzz_1,  \
                             g_0_yy_0_xyyyyy_1,   \
                             g_0_yy_0_xyyyyyy_0,  \
                             g_0_yy_0_xyyyyyy_1,  \
                             g_0_yy_0_xyyyyyz_0,  \
                             g_0_yy_0_xyyyyyz_1,  \
                             g_0_yy_0_xyyyyz_1,   \
                             g_0_yy_0_xyyyyzz_0,  \
                             g_0_yy_0_xyyyyzz_1,  \
                             g_0_yy_0_xyyyzz_1,   \
                             g_0_yy_0_xyyyzzz_0,  \
                             g_0_yy_0_xyyyzzz_1,  \
                             g_0_yy_0_xyyzzz_1,   \
                             g_0_yy_0_xyyzzzz_0,  \
                             g_0_yy_0_xyyzzzz_1,  \
                             g_0_yy_0_xyzzzz_1,   \
                             g_0_yy_0_xyzzzzz_0,  \
                             g_0_yy_0_xyzzzzz_1,  \
                             g_0_yy_0_xzzzzz_1,   \
                             g_0_yy_0_xzzzzzz_0,  \
                             g_0_yy_0_xzzzzzz_1,  \
                             g_0_yy_0_yyyyyy_1,   \
                             g_0_yy_0_yyyyyyy_0,  \
                             g_0_yy_0_yyyyyyy_1,  \
                             g_0_yy_0_yyyyyyz_0,  \
                             g_0_yy_0_yyyyyyz_1,  \
                             g_0_yy_0_yyyyyz_1,   \
                             g_0_yy_0_yyyyyzz_0,  \
                             g_0_yy_0_yyyyyzz_1,  \
                             g_0_yy_0_yyyyzz_1,   \
                             g_0_yy_0_yyyyzzz_0,  \
                             g_0_yy_0_yyyyzzz_1,  \
                             g_0_yy_0_yyyzzz_1,   \
                             g_0_yy_0_yyyzzzz_0,  \
                             g_0_yy_0_yyyzzzz_1,  \
                             g_0_yy_0_yyzzzz_1,   \
                             g_0_yy_0_yyzzzzz_0,  \
                             g_0_yy_0_yyzzzzz_1,  \
                             g_0_yy_0_yzzzzz_1,   \
                             g_0_yy_0_yzzzzzz_0,  \
                             g_0_yy_0_yzzzzzz_1,  \
                             g_0_yy_0_zzzzzz_1,   \
                             g_0_yy_0_zzzzzzz_0,  \
                             g_0_yy_0_zzzzzzz_1,  \
                             g_0_yyy_0_xxxxxxx_0, \
                             g_0_yyy_0_xxxxxxy_0, \
                             g_0_yyy_0_xxxxxxz_0, \
                             g_0_yyy_0_xxxxxyy_0, \
                             g_0_yyy_0_xxxxxyz_0, \
                             g_0_yyy_0_xxxxxzz_0, \
                             g_0_yyy_0_xxxxyyy_0, \
                             g_0_yyy_0_xxxxyyz_0, \
                             g_0_yyy_0_xxxxyzz_0, \
                             g_0_yyy_0_xxxxzzz_0, \
                             g_0_yyy_0_xxxyyyy_0, \
                             g_0_yyy_0_xxxyyyz_0, \
                             g_0_yyy_0_xxxyyzz_0, \
                             g_0_yyy_0_xxxyzzz_0, \
                             g_0_yyy_0_xxxzzzz_0, \
                             g_0_yyy_0_xxyyyyy_0, \
                             g_0_yyy_0_xxyyyyz_0, \
                             g_0_yyy_0_xxyyyzz_0, \
                             g_0_yyy_0_xxyyzzz_0, \
                             g_0_yyy_0_xxyzzzz_0, \
                             g_0_yyy_0_xxzzzzz_0, \
                             g_0_yyy_0_xyyyyyy_0, \
                             g_0_yyy_0_xyyyyyz_0, \
                             g_0_yyy_0_xyyyyzz_0, \
                             g_0_yyy_0_xyyyzzz_0, \
                             g_0_yyy_0_xyyzzzz_0, \
                             g_0_yyy_0_xyzzzzz_0, \
                             g_0_yyy_0_xzzzzzz_0, \
                             g_0_yyy_0_yyyyyyy_0, \
                             g_0_yyy_0_yyyyyyz_0, \
                             g_0_yyy_0_yyyyyzz_0, \
                             g_0_yyy_0_yyyyzzz_0, \
                             g_0_yyy_0_yyyzzzz_0, \
                             g_0_yyy_0_yyzzzzz_0, \
                             g_0_yyy_0_yzzzzzz_0, \
                             g_0_yyy_0_zzzzzzz_0, \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xxxxxxx_0[i] = 2.0 * g_0_y_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yy_0_xxxxxxx_0[i] * pb_y +
                                 g_0_yy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxxy_0[i] = 2.0 * g_0_y_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yy_0_xxxxxx_1[i] * fi_abcd_0 +
                                 g_0_yy_0_xxxxxxy_0[i] * pb_y + g_0_yy_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxxz_0[i] = 2.0 * g_0_y_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxxz_0[i] * pb_y +
                                 g_0_yy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxyy_0[i] = 2.0 * g_0_y_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxyy_1[i] * fti_ab_0 +
                                 2.0 * g_0_yy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyy_0[i] * pb_y + g_0_yy_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxyz_0[i] = 2.0 * g_0_y_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxz_1[i] * fi_abcd_0 +
                                 g_0_yy_0_xxxxxyz_0[i] * pb_y + g_0_yy_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxzz_0[i] = 2.0 * g_0_y_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxzz_0[i] * pb_y +
                                 g_0_yy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyyy_0[i] = 2.0 * g_0_y_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyyy_1[i] * fti_ab_0 +
                                 3.0 * g_0_yy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyy_0[i] * pb_y + g_0_yy_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyyz_0[i] = 2.0 * g_0_y_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyz_0[i] * pb_y + g_0_yy_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyzz_0[i] = 2.0 * g_0_y_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxzz_1[i] * fi_abcd_0 +
                                 g_0_yy_0_xxxxyzz_0[i] * pb_y + g_0_yy_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxzzz_0[i] = 2.0 * g_0_y_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxzzz_0[i] * pb_y +
                                 g_0_yy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyyy_0[i] = 2.0 * g_0_y_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyyy_1[i] * fti_ab_0 +
                                 4.0 * g_0_yy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyy_0[i] * pb_y + g_0_yy_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyyz_0[i] = 2.0 * g_0_y_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyyz_1[i] * fti_ab_0 +
                                 3.0 * g_0_yy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyz_0[i] * pb_y + g_0_yy_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyzz_0[i] = 2.0 * g_0_y_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyzz_0[i] * pb_y + g_0_yy_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyzzz_0[i] = 2.0 * g_0_y_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxzzz_1[i] * fi_abcd_0 +
                                 g_0_yy_0_xxxyzzz_0[i] * pb_y + g_0_yy_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxzzzz_0[i] = 2.0 * g_0_y_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxzzzz_0[i] * pb_y +
                                 g_0_yy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyyy_0[i] = 2.0 * g_0_y_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyyy_1[i] * fti_ab_0 +
                                 5.0 * g_0_yy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyy_0[i] * pb_y + g_0_yy_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyyz_0[i] = 2.0 * g_0_y_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyyz_1[i] * fti_ab_0 +
                                 4.0 * g_0_yy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyz_0[i] * pb_y + g_0_yy_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyzz_0[i] = 2.0 * g_0_y_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_yy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyzz_0[i] * pb_y + g_0_yy_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyzzz_0[i] = 2.0 * g_0_y_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyzzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzzz_0[i] * pb_y + g_0_yy_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyzzzz_0[i] = 2.0 * g_0_y_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxzzzz_1[i] * fi_abcd_0 +
                                 g_0_yy_0_xxyzzzz_0[i] * pb_y + g_0_yy_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxzzzzz_0[i] = 2.0 * g_0_y_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxzzzzz_0[i] * pb_y +
                                 g_0_yy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyyy_0[i] = 2.0 * g_0_y_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyyy_1[i] * fti_ab_0 +
                                 6.0 * g_0_yy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyy_0[i] * pb_y + g_0_yy_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyyz_0[i] = 2.0 * g_0_y_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyyz_1[i] * fti_ab_0 +
                                 5.0 * g_0_yy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyz_0[i] * pb_y + g_0_yy_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyzz_0[i] = 2.0 * g_0_y_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_yy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyzz_0[i] * pb_y + g_0_yy_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyzzz_0[i] = 2.0 * g_0_y_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_yy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzzz_0[i] * pb_y + g_0_yy_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyzzzz_0[i] = 2.0 * g_0_y_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyzzzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzzz_0[i] * pb_y + g_0_yy_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyzzzzz_0[i] = 2.0 * g_0_y_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzzzz_1[i] * fi_abcd_0 +
                                 g_0_yy_0_xyzzzzz_0[i] * pb_y + g_0_yy_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xzzzzzz_0[i] = 2.0 * g_0_y_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzzzzz_0[i] * pb_y +
                                 g_0_yy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyyy_0[i] = 2.0 * g_0_y_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyyy_1[i] * fti_ab_0 +
                                 7.0 * g_0_yy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyy_0[i] * pb_y + g_0_yy_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyyz_0[i] = 2.0 * g_0_y_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyyz_1[i] * fti_ab_0 +
                                 6.0 * g_0_yy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyz_0[i] * pb_y + g_0_yy_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyzz_0[i] = 2.0 * g_0_y_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyzz_1[i] * fti_ab_0 +
                                 5.0 * g_0_yy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyzz_0[i] * pb_y + g_0_yy_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyzzz_0[i] = 2.0 * g_0_y_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyzzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_yy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyzzz_0[i] * pb_y + g_0_yy_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyzzzz_0[i] = 2.0 * g_0_y_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyzzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_yy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyzzzz_0[i] * pb_y + g_0_yy_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyzzzzz_0[i] = 2.0 * g_0_y_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyzzzzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzzzzz_0[i] * pb_y + g_0_yy_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yzzzzzz_0[i] = 2.0 * g_0_y_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzzzz_1[i] * fi_abcd_0 +
                                 g_0_yy_0_yzzzzzz_0[i] * pb_y + g_0_yy_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_zzzzzzz_0[i] = 2.0 * g_0_y_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzzzzz_0[i] * pb_y +
                                 g_0_yy_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 252-288 components of targeted buffer : SFSK

    auto g_0_yyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 252);

    auto g_0_yyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 253);

    auto g_0_yyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 254);

    auto g_0_yyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 255);

    auto g_0_yyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 256);

    auto g_0_yyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 257);

    auto g_0_yyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 258);

    auto g_0_yyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 259);

    auto g_0_yyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 260);

    auto g_0_yyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 261);

    auto g_0_yyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 262);

    auto g_0_yyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 263);

    auto g_0_yyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 264);

    auto g_0_yyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 265);

    auto g_0_yyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 266);

    auto g_0_yyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 267);

    auto g_0_yyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 268);

    auto g_0_yyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 269);

    auto g_0_yyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 270);

    auto g_0_yyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 271);

    auto g_0_yyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 272);

    auto g_0_yyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 273);

    auto g_0_yyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 274);

    auto g_0_yyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 275);

    auto g_0_yyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 276);

    auto g_0_yyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 277);

    auto g_0_yyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 278);

    auto g_0_yyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 279);

    auto g_0_yyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 280);

    auto g_0_yyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 281);

    auto g_0_yyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 282);

    auto g_0_yyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 283);

    auto g_0_yyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 284);

    auto g_0_yyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 285);

    auto g_0_yyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 286);

    auto g_0_yyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 287);

#pragma omp simd aligned(g_0_yy_0_xxxxxx_1,       \
                             g_0_yy_0_xxxxxxx_0,  \
                             g_0_yy_0_xxxxxxx_1,  \
                             g_0_yy_0_xxxxxxy_0,  \
                             g_0_yy_0_xxxxxxy_1,  \
                             g_0_yy_0_xxxxxxz_0,  \
                             g_0_yy_0_xxxxxxz_1,  \
                             g_0_yy_0_xxxxxy_1,   \
                             g_0_yy_0_xxxxxyy_0,  \
                             g_0_yy_0_xxxxxyy_1,  \
                             g_0_yy_0_xxxxxyz_0,  \
                             g_0_yy_0_xxxxxyz_1,  \
                             g_0_yy_0_xxxxxz_1,   \
                             g_0_yy_0_xxxxxzz_0,  \
                             g_0_yy_0_xxxxxzz_1,  \
                             g_0_yy_0_xxxxyy_1,   \
                             g_0_yy_0_xxxxyyy_0,  \
                             g_0_yy_0_xxxxyyy_1,  \
                             g_0_yy_0_xxxxyyz_0,  \
                             g_0_yy_0_xxxxyyz_1,  \
                             g_0_yy_0_xxxxyz_1,   \
                             g_0_yy_0_xxxxyzz_0,  \
                             g_0_yy_0_xxxxyzz_1,  \
                             g_0_yy_0_xxxxzz_1,   \
                             g_0_yy_0_xxxxzzz_0,  \
                             g_0_yy_0_xxxxzzz_1,  \
                             g_0_yy_0_xxxyyy_1,   \
                             g_0_yy_0_xxxyyyy_0,  \
                             g_0_yy_0_xxxyyyy_1,  \
                             g_0_yy_0_xxxyyyz_0,  \
                             g_0_yy_0_xxxyyyz_1,  \
                             g_0_yy_0_xxxyyz_1,   \
                             g_0_yy_0_xxxyyzz_0,  \
                             g_0_yy_0_xxxyyzz_1,  \
                             g_0_yy_0_xxxyzz_1,   \
                             g_0_yy_0_xxxyzzz_0,  \
                             g_0_yy_0_xxxyzzz_1,  \
                             g_0_yy_0_xxxzzz_1,   \
                             g_0_yy_0_xxxzzzz_0,  \
                             g_0_yy_0_xxxzzzz_1,  \
                             g_0_yy_0_xxyyyy_1,   \
                             g_0_yy_0_xxyyyyy_0,  \
                             g_0_yy_0_xxyyyyy_1,  \
                             g_0_yy_0_xxyyyyz_0,  \
                             g_0_yy_0_xxyyyyz_1,  \
                             g_0_yy_0_xxyyyz_1,   \
                             g_0_yy_0_xxyyyzz_0,  \
                             g_0_yy_0_xxyyyzz_1,  \
                             g_0_yy_0_xxyyzz_1,   \
                             g_0_yy_0_xxyyzzz_0,  \
                             g_0_yy_0_xxyyzzz_1,  \
                             g_0_yy_0_xxyzzz_1,   \
                             g_0_yy_0_xxyzzzz_0,  \
                             g_0_yy_0_xxyzzzz_1,  \
                             g_0_yy_0_xxzzzz_1,   \
                             g_0_yy_0_xxzzzzz_0,  \
                             g_0_yy_0_xxzzzzz_1,  \
                             g_0_yy_0_xyyyyy_1,   \
                             g_0_yy_0_xyyyyyy_0,  \
                             g_0_yy_0_xyyyyyy_1,  \
                             g_0_yy_0_xyyyyyz_0,  \
                             g_0_yy_0_xyyyyyz_1,  \
                             g_0_yy_0_xyyyyz_1,   \
                             g_0_yy_0_xyyyyzz_0,  \
                             g_0_yy_0_xyyyyzz_1,  \
                             g_0_yy_0_xyyyzz_1,   \
                             g_0_yy_0_xyyyzzz_0,  \
                             g_0_yy_0_xyyyzzz_1,  \
                             g_0_yy_0_xyyzzz_1,   \
                             g_0_yy_0_xyyzzzz_0,  \
                             g_0_yy_0_xyyzzzz_1,  \
                             g_0_yy_0_xyzzzz_1,   \
                             g_0_yy_0_xyzzzzz_0,  \
                             g_0_yy_0_xyzzzzz_1,  \
                             g_0_yy_0_xzzzzz_1,   \
                             g_0_yy_0_xzzzzzz_0,  \
                             g_0_yy_0_xzzzzzz_1,  \
                             g_0_yy_0_yyyyyy_1,   \
                             g_0_yy_0_yyyyyyy_0,  \
                             g_0_yy_0_yyyyyyy_1,  \
                             g_0_yy_0_yyyyyyz_0,  \
                             g_0_yy_0_yyyyyyz_1,  \
                             g_0_yy_0_yyyyyz_1,   \
                             g_0_yy_0_yyyyyzz_0,  \
                             g_0_yy_0_yyyyyzz_1,  \
                             g_0_yy_0_yyyyzz_1,   \
                             g_0_yy_0_yyyyzzz_0,  \
                             g_0_yy_0_yyyyzzz_1,  \
                             g_0_yy_0_yyyzzz_1,   \
                             g_0_yy_0_yyyzzzz_0,  \
                             g_0_yy_0_yyyzzzz_1,  \
                             g_0_yy_0_yyzzzz_1,   \
                             g_0_yy_0_yyzzzzz_0,  \
                             g_0_yy_0_yyzzzzz_1,  \
                             g_0_yy_0_yzzzzz_1,   \
                             g_0_yy_0_yzzzzzz_0,  \
                             g_0_yy_0_yzzzzzz_1,  \
                             g_0_yy_0_zzzzzz_1,   \
                             g_0_yy_0_zzzzzzz_0,  \
                             g_0_yy_0_zzzzzzz_1,  \
                             g_0_yyz_0_xxxxxxx_0, \
                             g_0_yyz_0_xxxxxxy_0, \
                             g_0_yyz_0_xxxxxxz_0, \
                             g_0_yyz_0_xxxxxyy_0, \
                             g_0_yyz_0_xxxxxyz_0, \
                             g_0_yyz_0_xxxxxzz_0, \
                             g_0_yyz_0_xxxxyyy_0, \
                             g_0_yyz_0_xxxxyyz_0, \
                             g_0_yyz_0_xxxxyzz_0, \
                             g_0_yyz_0_xxxxzzz_0, \
                             g_0_yyz_0_xxxyyyy_0, \
                             g_0_yyz_0_xxxyyyz_0, \
                             g_0_yyz_0_xxxyyzz_0, \
                             g_0_yyz_0_xxxyzzz_0, \
                             g_0_yyz_0_xxxzzzz_0, \
                             g_0_yyz_0_xxyyyyy_0, \
                             g_0_yyz_0_xxyyyyz_0, \
                             g_0_yyz_0_xxyyyzz_0, \
                             g_0_yyz_0_xxyyzzz_0, \
                             g_0_yyz_0_xxyzzzz_0, \
                             g_0_yyz_0_xxzzzzz_0, \
                             g_0_yyz_0_xyyyyyy_0, \
                             g_0_yyz_0_xyyyyyz_0, \
                             g_0_yyz_0_xyyyyzz_0, \
                             g_0_yyz_0_xyyyzzz_0, \
                             g_0_yyz_0_xyyzzzz_0, \
                             g_0_yyz_0_xyzzzzz_0, \
                             g_0_yyz_0_xzzzzzz_0, \
                             g_0_yyz_0_yyyyyyy_0, \
                             g_0_yyz_0_yyyyyyz_0, \
                             g_0_yyz_0_yyyyyzz_0, \
                             g_0_yyz_0_yyyyzzz_0, \
                             g_0_yyz_0_yyyzzzz_0, \
                             g_0_yyz_0_yyzzzzz_0, \
                             g_0_yyz_0_yzzzzzz_0, \
                             g_0_yyz_0_zzzzzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xxxxxxx_0[i] = g_0_yy_0_xxxxxxx_0[i] * pb_z + g_0_yy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxxy_0[i] = g_0_yy_0_xxxxxxy_0[i] * pb_z + g_0_yy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxxz_0[i] = g_0_yy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxxz_0[i] * pb_z + g_0_yy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxyy_0[i] = g_0_yy_0_xxxxxyy_0[i] * pb_z + g_0_yy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxyz_0[i] = g_0_yy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxyz_0[i] * pb_z + g_0_yy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxzz_0[i] = 2.0 * g_0_yy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxzz_0[i] * pb_z + g_0_yy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyyy_0[i] = g_0_yy_0_xxxxyyy_0[i] * pb_z + g_0_yy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyyz_0[i] = g_0_yy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyyz_0[i] * pb_z + g_0_yy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyzz_0[i] = 2.0 * g_0_yy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyzz_0[i] * pb_z + g_0_yy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxzzz_0[i] = 3.0 * g_0_yy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxzzz_0[i] * pb_z + g_0_yy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyyy_0[i] = g_0_yy_0_xxxyyyy_0[i] * pb_z + g_0_yy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyyz_0[i] = g_0_yy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyyz_0[i] * pb_z + g_0_yy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyzz_0[i] = 2.0 * g_0_yy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyzz_0[i] * pb_z + g_0_yy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyzzz_0[i] = 3.0 * g_0_yy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyzzz_0[i] * pb_z + g_0_yy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxzzzz_0[i] = 4.0 * g_0_yy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzzzz_0[i] * pb_z + g_0_yy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyyy_0[i] = g_0_yy_0_xxyyyyy_0[i] * pb_z + g_0_yy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyyz_0[i] = g_0_yy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyyz_0[i] * pb_z + g_0_yy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyzz_0[i] = 2.0 * g_0_yy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyzz_0[i] * pb_z + g_0_yy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyzzz_0[i] = 3.0 * g_0_yy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzzz_0[i] * pb_z + g_0_yy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyzzzz_0[i] = 4.0 * g_0_yy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzzzz_0[i] * pb_z + g_0_yy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxzzzzz_0[i] = 5.0 * g_0_yy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzzzz_0[i] * pb_z + g_0_yy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyyy_0[i] = g_0_yy_0_xyyyyyy_0[i] * pb_z + g_0_yy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyyz_0[i] = g_0_yy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyyz_0[i] * pb_z + g_0_yy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyzz_0[i] = 2.0 * g_0_yy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyzz_0[i] * pb_z + g_0_yy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyzzz_0[i] = 3.0 * g_0_yy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzzz_0[i] * pb_z + g_0_yy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyzzzz_0[i] = 4.0 * g_0_yy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzzz_0[i] * pb_z + g_0_yy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyzzzzz_0[i] = 5.0 * g_0_yy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzzzz_0[i] * pb_z + g_0_yy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xzzzzzz_0[i] = 6.0 * g_0_yy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzzzz_0[i] * pb_z + g_0_yy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyyy_0[i] = g_0_yy_0_yyyyyyy_0[i] * pb_z + g_0_yy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyyz_0[i] = g_0_yy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyyz_0[i] * pb_z + g_0_yy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyzz_0[i] = 2.0 * g_0_yy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyzz_0[i] * pb_z + g_0_yy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyzzz_0[i] = 3.0 * g_0_yy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyzzz_0[i] * pb_z + g_0_yy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyzzzz_0[i] = 4.0 * g_0_yy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyzzzz_0[i] * pb_z + g_0_yy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyzzzzz_0[i] = 5.0 * g_0_yy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzzzzz_0[i] * pb_z + g_0_yy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yzzzzzz_0[i] = 6.0 * g_0_yy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yzzzzzz_0[i] * pb_z + g_0_yy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_zzzzzzz_0[i] = 7.0 * g_0_yy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yy_0_zzzzzzz_0[i] * pb_z + g_0_yy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 288-324 components of targeted buffer : SFSK

    auto g_0_yzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 288);

    auto g_0_yzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 289);

    auto g_0_yzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 290);

    auto g_0_yzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 291);

    auto g_0_yzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 292);

    auto g_0_yzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 293);

    auto g_0_yzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 294);

    auto g_0_yzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 295);

    auto g_0_yzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 296);

    auto g_0_yzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 297);

    auto g_0_yzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 298);

    auto g_0_yzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 299);

    auto g_0_yzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 300);

    auto g_0_yzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 301);

    auto g_0_yzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 302);

    auto g_0_yzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 303);

    auto g_0_yzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 304);

    auto g_0_yzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 305);

    auto g_0_yzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 306);

    auto g_0_yzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 307);

    auto g_0_yzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 308);

    auto g_0_yzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 309);

    auto g_0_yzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 310);

    auto g_0_yzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 311);

    auto g_0_yzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 312);

    auto g_0_yzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 313);

    auto g_0_yzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 314);

    auto g_0_yzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 315);

    auto g_0_yzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 316);

    auto g_0_yzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 317);

    auto g_0_yzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 318);

    auto g_0_yzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 319);

    auto g_0_yzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 320);

    auto g_0_yzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 321);

    auto g_0_yzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 322);

    auto g_0_yzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 323);

#pragma omp simd aligned(g_0_yzz_0_xxxxxxx_0,     \
                             g_0_yzz_0_xxxxxxy_0, \
                             g_0_yzz_0_xxxxxxz_0, \
                             g_0_yzz_0_xxxxxyy_0, \
                             g_0_yzz_0_xxxxxyz_0, \
                             g_0_yzz_0_xxxxxzz_0, \
                             g_0_yzz_0_xxxxyyy_0, \
                             g_0_yzz_0_xxxxyyz_0, \
                             g_0_yzz_0_xxxxyzz_0, \
                             g_0_yzz_0_xxxxzzz_0, \
                             g_0_yzz_0_xxxyyyy_0, \
                             g_0_yzz_0_xxxyyyz_0, \
                             g_0_yzz_0_xxxyyzz_0, \
                             g_0_yzz_0_xxxyzzz_0, \
                             g_0_yzz_0_xxxzzzz_0, \
                             g_0_yzz_0_xxyyyyy_0, \
                             g_0_yzz_0_xxyyyyz_0, \
                             g_0_yzz_0_xxyyyzz_0, \
                             g_0_yzz_0_xxyyzzz_0, \
                             g_0_yzz_0_xxyzzzz_0, \
                             g_0_yzz_0_xxzzzzz_0, \
                             g_0_yzz_0_xyyyyyy_0, \
                             g_0_yzz_0_xyyyyyz_0, \
                             g_0_yzz_0_xyyyyzz_0, \
                             g_0_yzz_0_xyyyzzz_0, \
                             g_0_yzz_0_xyyzzzz_0, \
                             g_0_yzz_0_xyzzzzz_0, \
                             g_0_yzz_0_xzzzzzz_0, \
                             g_0_yzz_0_yyyyyyy_0, \
                             g_0_yzz_0_yyyyyyz_0, \
                             g_0_yzz_0_yyyyyzz_0, \
                             g_0_yzz_0_yyyyzzz_0, \
                             g_0_yzz_0_yyyzzzz_0, \
                             g_0_yzz_0_yyzzzzz_0, \
                             g_0_yzz_0_yzzzzzz_0, \
                             g_0_yzz_0_zzzzzzz_0, \
                             g_0_zz_0_xxxxxx_1,   \
                             g_0_zz_0_xxxxxxx_0,  \
                             g_0_zz_0_xxxxxxx_1,  \
                             g_0_zz_0_xxxxxxy_0,  \
                             g_0_zz_0_xxxxxxy_1,  \
                             g_0_zz_0_xxxxxxz_0,  \
                             g_0_zz_0_xxxxxxz_1,  \
                             g_0_zz_0_xxxxxy_1,   \
                             g_0_zz_0_xxxxxyy_0,  \
                             g_0_zz_0_xxxxxyy_1,  \
                             g_0_zz_0_xxxxxyz_0,  \
                             g_0_zz_0_xxxxxyz_1,  \
                             g_0_zz_0_xxxxxz_1,   \
                             g_0_zz_0_xxxxxzz_0,  \
                             g_0_zz_0_xxxxxzz_1,  \
                             g_0_zz_0_xxxxyy_1,   \
                             g_0_zz_0_xxxxyyy_0,  \
                             g_0_zz_0_xxxxyyy_1,  \
                             g_0_zz_0_xxxxyyz_0,  \
                             g_0_zz_0_xxxxyyz_1,  \
                             g_0_zz_0_xxxxyz_1,   \
                             g_0_zz_0_xxxxyzz_0,  \
                             g_0_zz_0_xxxxyzz_1,  \
                             g_0_zz_0_xxxxzz_1,   \
                             g_0_zz_0_xxxxzzz_0,  \
                             g_0_zz_0_xxxxzzz_1,  \
                             g_0_zz_0_xxxyyy_1,   \
                             g_0_zz_0_xxxyyyy_0,  \
                             g_0_zz_0_xxxyyyy_1,  \
                             g_0_zz_0_xxxyyyz_0,  \
                             g_0_zz_0_xxxyyyz_1,  \
                             g_0_zz_0_xxxyyz_1,   \
                             g_0_zz_0_xxxyyzz_0,  \
                             g_0_zz_0_xxxyyzz_1,  \
                             g_0_zz_0_xxxyzz_1,   \
                             g_0_zz_0_xxxyzzz_0,  \
                             g_0_zz_0_xxxyzzz_1,  \
                             g_0_zz_0_xxxzzz_1,   \
                             g_0_zz_0_xxxzzzz_0,  \
                             g_0_zz_0_xxxzzzz_1,  \
                             g_0_zz_0_xxyyyy_1,   \
                             g_0_zz_0_xxyyyyy_0,  \
                             g_0_zz_0_xxyyyyy_1,  \
                             g_0_zz_0_xxyyyyz_0,  \
                             g_0_zz_0_xxyyyyz_1,  \
                             g_0_zz_0_xxyyyz_1,   \
                             g_0_zz_0_xxyyyzz_0,  \
                             g_0_zz_0_xxyyyzz_1,  \
                             g_0_zz_0_xxyyzz_1,   \
                             g_0_zz_0_xxyyzzz_0,  \
                             g_0_zz_0_xxyyzzz_1,  \
                             g_0_zz_0_xxyzzz_1,   \
                             g_0_zz_0_xxyzzzz_0,  \
                             g_0_zz_0_xxyzzzz_1,  \
                             g_0_zz_0_xxzzzz_1,   \
                             g_0_zz_0_xxzzzzz_0,  \
                             g_0_zz_0_xxzzzzz_1,  \
                             g_0_zz_0_xyyyyy_1,   \
                             g_0_zz_0_xyyyyyy_0,  \
                             g_0_zz_0_xyyyyyy_1,  \
                             g_0_zz_0_xyyyyyz_0,  \
                             g_0_zz_0_xyyyyyz_1,  \
                             g_0_zz_0_xyyyyz_1,   \
                             g_0_zz_0_xyyyyzz_0,  \
                             g_0_zz_0_xyyyyzz_1,  \
                             g_0_zz_0_xyyyzz_1,   \
                             g_0_zz_0_xyyyzzz_0,  \
                             g_0_zz_0_xyyyzzz_1,  \
                             g_0_zz_0_xyyzzz_1,   \
                             g_0_zz_0_xyyzzzz_0,  \
                             g_0_zz_0_xyyzzzz_1,  \
                             g_0_zz_0_xyzzzz_1,   \
                             g_0_zz_0_xyzzzzz_0,  \
                             g_0_zz_0_xyzzzzz_1,  \
                             g_0_zz_0_xzzzzz_1,   \
                             g_0_zz_0_xzzzzzz_0,  \
                             g_0_zz_0_xzzzzzz_1,  \
                             g_0_zz_0_yyyyyy_1,   \
                             g_0_zz_0_yyyyyyy_0,  \
                             g_0_zz_0_yyyyyyy_1,  \
                             g_0_zz_0_yyyyyyz_0,  \
                             g_0_zz_0_yyyyyyz_1,  \
                             g_0_zz_0_yyyyyz_1,   \
                             g_0_zz_0_yyyyyzz_0,  \
                             g_0_zz_0_yyyyyzz_1,  \
                             g_0_zz_0_yyyyzz_1,   \
                             g_0_zz_0_yyyyzzz_0,  \
                             g_0_zz_0_yyyyzzz_1,  \
                             g_0_zz_0_yyyzzz_1,   \
                             g_0_zz_0_yyyzzzz_0,  \
                             g_0_zz_0_yyyzzzz_1,  \
                             g_0_zz_0_yyzzzz_1,   \
                             g_0_zz_0_yyzzzzz_0,  \
                             g_0_zz_0_yyzzzzz_1,  \
                             g_0_zz_0_yzzzzz_1,   \
                             g_0_zz_0_yzzzzzz_0,  \
                             g_0_zz_0_yzzzzzz_1,  \
                             g_0_zz_0_zzzzzz_1,   \
                             g_0_zz_0_zzzzzzz_0,  \
                             g_0_zz_0_zzzzzzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xxxxxxx_0[i] = g_0_zz_0_xxxxxxx_0[i] * pb_y + g_0_zz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxxy_0[i] = g_0_zz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxxy_0[i] * pb_y + g_0_zz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxxz_0[i] = g_0_zz_0_xxxxxxz_0[i] * pb_y + g_0_zz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxyy_0[i] = 2.0 * g_0_zz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyy_0[i] * pb_y + g_0_zz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxyz_0[i] = g_0_zz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxyz_0[i] * pb_y + g_0_zz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxzz_0[i] = g_0_zz_0_xxxxxzz_0[i] * pb_y + g_0_zz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyyy_0[i] = 3.0 * g_0_zz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyy_0[i] * pb_y + g_0_zz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyyz_0[i] = 2.0 * g_0_zz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyyz_0[i] * pb_y + g_0_zz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyzz_0[i] = g_0_zz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyzz_0[i] * pb_y + g_0_zz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxzzz_0[i] = g_0_zz_0_xxxxzzz_0[i] * pb_y + g_0_zz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyyy_0[i] = 4.0 * g_0_zz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyy_0[i] * pb_y + g_0_zz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyyz_0[i] = 3.0 * g_0_zz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyyz_0[i] * pb_y + g_0_zz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyzz_0[i] = 2.0 * g_0_zz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyzz_0[i] * pb_y + g_0_zz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyzzz_0[i] = g_0_zz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzzz_0[i] * pb_y + g_0_zz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxzzzz_0[i] = g_0_zz_0_xxxzzzz_0[i] * pb_y + g_0_zz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyyy_0[i] = 5.0 * g_0_zz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyy_0[i] * pb_y + g_0_zz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyyz_0[i] = 4.0 * g_0_zz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyyz_0[i] * pb_y + g_0_zz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyzz_0[i] = 3.0 * g_0_zz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyzz_0[i] * pb_y + g_0_zz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyzzz_0[i] = 2.0 * g_0_zz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzzz_0[i] * pb_y + g_0_zz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyzzzz_0[i] = g_0_zz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzzz_0[i] * pb_y + g_0_zz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxzzzzz_0[i] = g_0_zz_0_xxzzzzz_0[i] * pb_y + g_0_zz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyyy_0[i] = 6.0 * g_0_zz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyy_0[i] * pb_y + g_0_zz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyyz_0[i] = 5.0 * g_0_zz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyyz_0[i] * pb_y + g_0_zz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyzz_0[i] = 4.0 * g_0_zz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyzz_0[i] * pb_y + g_0_zz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyzzz_0[i] = 3.0 * g_0_zz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzzz_0[i] * pb_y + g_0_zz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyzzzz_0[i] = 2.0 * g_0_zz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzzz_0[i] * pb_y + g_0_zz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyzzzzz_0[i] = g_0_zz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzzz_0[i] * pb_y + g_0_zz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xzzzzzz_0[i] = g_0_zz_0_xzzzzzz_0[i] * pb_y + g_0_zz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyyy_0[i] = 7.0 * g_0_zz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyyy_0[i] * pb_y + g_0_zz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyyz_0[i] = 6.0 * g_0_zz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyyz_0[i] * pb_y + g_0_zz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyzz_0[i] = 5.0 * g_0_zz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyzz_0[i] * pb_y + g_0_zz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyzzz_0[i] = 4.0 * g_0_zz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyzzz_0[i] * pb_y + g_0_zz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyzzzz_0[i] = 3.0 * g_0_zz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyzzzz_0[i] * pb_y + g_0_zz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyzzzzz_0[i] = 2.0 * g_0_zz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzzzzz_0[i] * pb_y + g_0_zz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yzzzzzz_0[i] = g_0_zz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzzzzz_0[i] * pb_y + g_0_zz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_zzzzzzz_0[i] = g_0_zz_0_zzzzzzz_0[i] * pb_y + g_0_zz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 324-360 components of targeted buffer : SFSK

    auto g_0_zzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 324);

    auto g_0_zzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 325);

    auto g_0_zzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 326);

    auto g_0_zzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 327);

    auto g_0_zzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 328);

    auto g_0_zzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 329);

    auto g_0_zzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 330);

    auto g_0_zzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 331);

    auto g_0_zzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 332);

    auto g_0_zzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 333);

    auto g_0_zzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 334);

    auto g_0_zzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 335);

    auto g_0_zzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 336);

    auto g_0_zzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 337);

    auto g_0_zzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 338);

    auto g_0_zzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 339);

    auto g_0_zzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 340);

    auto g_0_zzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 341);

    auto g_0_zzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 342);

    auto g_0_zzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 343);

    auto g_0_zzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 344);

    auto g_0_zzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 345);

    auto g_0_zzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 346);

    auto g_0_zzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 347);

    auto g_0_zzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 348);

    auto g_0_zzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 349);

    auto g_0_zzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 350);

    auto g_0_zzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 351);

    auto g_0_zzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 352);

    auto g_0_zzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 353);

    auto g_0_zzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 354);

    auto g_0_zzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 355);

    auto g_0_zzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 356);

    auto g_0_zzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 357);

    auto g_0_zzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 358);

    auto g_0_zzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 359);

#pragma omp simd aligned(g_0_z_0_xxxxxxx_0,       \
                             g_0_z_0_xxxxxxx_1,   \
                             g_0_z_0_xxxxxxy_0,   \
                             g_0_z_0_xxxxxxy_1,   \
                             g_0_z_0_xxxxxxz_0,   \
                             g_0_z_0_xxxxxxz_1,   \
                             g_0_z_0_xxxxxyy_0,   \
                             g_0_z_0_xxxxxyy_1,   \
                             g_0_z_0_xxxxxyz_0,   \
                             g_0_z_0_xxxxxyz_1,   \
                             g_0_z_0_xxxxxzz_0,   \
                             g_0_z_0_xxxxxzz_1,   \
                             g_0_z_0_xxxxyyy_0,   \
                             g_0_z_0_xxxxyyy_1,   \
                             g_0_z_0_xxxxyyz_0,   \
                             g_0_z_0_xxxxyyz_1,   \
                             g_0_z_0_xxxxyzz_0,   \
                             g_0_z_0_xxxxyzz_1,   \
                             g_0_z_0_xxxxzzz_0,   \
                             g_0_z_0_xxxxzzz_1,   \
                             g_0_z_0_xxxyyyy_0,   \
                             g_0_z_0_xxxyyyy_1,   \
                             g_0_z_0_xxxyyyz_0,   \
                             g_0_z_0_xxxyyyz_1,   \
                             g_0_z_0_xxxyyzz_0,   \
                             g_0_z_0_xxxyyzz_1,   \
                             g_0_z_0_xxxyzzz_0,   \
                             g_0_z_0_xxxyzzz_1,   \
                             g_0_z_0_xxxzzzz_0,   \
                             g_0_z_0_xxxzzzz_1,   \
                             g_0_z_0_xxyyyyy_0,   \
                             g_0_z_0_xxyyyyy_1,   \
                             g_0_z_0_xxyyyyz_0,   \
                             g_0_z_0_xxyyyyz_1,   \
                             g_0_z_0_xxyyyzz_0,   \
                             g_0_z_0_xxyyyzz_1,   \
                             g_0_z_0_xxyyzzz_0,   \
                             g_0_z_0_xxyyzzz_1,   \
                             g_0_z_0_xxyzzzz_0,   \
                             g_0_z_0_xxyzzzz_1,   \
                             g_0_z_0_xxzzzzz_0,   \
                             g_0_z_0_xxzzzzz_1,   \
                             g_0_z_0_xyyyyyy_0,   \
                             g_0_z_0_xyyyyyy_1,   \
                             g_0_z_0_xyyyyyz_0,   \
                             g_0_z_0_xyyyyyz_1,   \
                             g_0_z_0_xyyyyzz_0,   \
                             g_0_z_0_xyyyyzz_1,   \
                             g_0_z_0_xyyyzzz_0,   \
                             g_0_z_0_xyyyzzz_1,   \
                             g_0_z_0_xyyzzzz_0,   \
                             g_0_z_0_xyyzzzz_1,   \
                             g_0_z_0_xyzzzzz_0,   \
                             g_0_z_0_xyzzzzz_1,   \
                             g_0_z_0_xzzzzzz_0,   \
                             g_0_z_0_xzzzzzz_1,   \
                             g_0_z_0_yyyyyyy_0,   \
                             g_0_z_0_yyyyyyy_1,   \
                             g_0_z_0_yyyyyyz_0,   \
                             g_0_z_0_yyyyyyz_1,   \
                             g_0_z_0_yyyyyzz_0,   \
                             g_0_z_0_yyyyyzz_1,   \
                             g_0_z_0_yyyyzzz_0,   \
                             g_0_z_0_yyyyzzz_1,   \
                             g_0_z_0_yyyzzzz_0,   \
                             g_0_z_0_yyyzzzz_1,   \
                             g_0_z_0_yyzzzzz_0,   \
                             g_0_z_0_yyzzzzz_1,   \
                             g_0_z_0_yzzzzzz_0,   \
                             g_0_z_0_yzzzzzz_1,   \
                             g_0_z_0_zzzzzzz_0,   \
                             g_0_z_0_zzzzzzz_1,   \
                             g_0_zz_0_xxxxxx_1,   \
                             g_0_zz_0_xxxxxxx_0,  \
                             g_0_zz_0_xxxxxxx_1,  \
                             g_0_zz_0_xxxxxxy_0,  \
                             g_0_zz_0_xxxxxxy_1,  \
                             g_0_zz_0_xxxxxxz_0,  \
                             g_0_zz_0_xxxxxxz_1,  \
                             g_0_zz_0_xxxxxy_1,   \
                             g_0_zz_0_xxxxxyy_0,  \
                             g_0_zz_0_xxxxxyy_1,  \
                             g_0_zz_0_xxxxxyz_0,  \
                             g_0_zz_0_xxxxxyz_1,  \
                             g_0_zz_0_xxxxxz_1,   \
                             g_0_zz_0_xxxxxzz_0,  \
                             g_0_zz_0_xxxxxzz_1,  \
                             g_0_zz_0_xxxxyy_1,   \
                             g_0_zz_0_xxxxyyy_0,  \
                             g_0_zz_0_xxxxyyy_1,  \
                             g_0_zz_0_xxxxyyz_0,  \
                             g_0_zz_0_xxxxyyz_1,  \
                             g_0_zz_0_xxxxyz_1,   \
                             g_0_zz_0_xxxxyzz_0,  \
                             g_0_zz_0_xxxxyzz_1,  \
                             g_0_zz_0_xxxxzz_1,   \
                             g_0_zz_0_xxxxzzz_0,  \
                             g_0_zz_0_xxxxzzz_1,  \
                             g_0_zz_0_xxxyyy_1,   \
                             g_0_zz_0_xxxyyyy_0,  \
                             g_0_zz_0_xxxyyyy_1,  \
                             g_0_zz_0_xxxyyyz_0,  \
                             g_0_zz_0_xxxyyyz_1,  \
                             g_0_zz_0_xxxyyz_1,   \
                             g_0_zz_0_xxxyyzz_0,  \
                             g_0_zz_0_xxxyyzz_1,  \
                             g_0_zz_0_xxxyzz_1,   \
                             g_0_zz_0_xxxyzzz_0,  \
                             g_0_zz_0_xxxyzzz_1,  \
                             g_0_zz_0_xxxzzz_1,   \
                             g_0_zz_0_xxxzzzz_0,  \
                             g_0_zz_0_xxxzzzz_1,  \
                             g_0_zz_0_xxyyyy_1,   \
                             g_0_zz_0_xxyyyyy_0,  \
                             g_0_zz_0_xxyyyyy_1,  \
                             g_0_zz_0_xxyyyyz_0,  \
                             g_0_zz_0_xxyyyyz_1,  \
                             g_0_zz_0_xxyyyz_1,   \
                             g_0_zz_0_xxyyyzz_0,  \
                             g_0_zz_0_xxyyyzz_1,  \
                             g_0_zz_0_xxyyzz_1,   \
                             g_0_zz_0_xxyyzzz_0,  \
                             g_0_zz_0_xxyyzzz_1,  \
                             g_0_zz_0_xxyzzz_1,   \
                             g_0_zz_0_xxyzzzz_0,  \
                             g_0_zz_0_xxyzzzz_1,  \
                             g_0_zz_0_xxzzzz_1,   \
                             g_0_zz_0_xxzzzzz_0,  \
                             g_0_zz_0_xxzzzzz_1,  \
                             g_0_zz_0_xyyyyy_1,   \
                             g_0_zz_0_xyyyyyy_0,  \
                             g_0_zz_0_xyyyyyy_1,  \
                             g_0_zz_0_xyyyyyz_0,  \
                             g_0_zz_0_xyyyyyz_1,  \
                             g_0_zz_0_xyyyyz_1,   \
                             g_0_zz_0_xyyyyzz_0,  \
                             g_0_zz_0_xyyyyzz_1,  \
                             g_0_zz_0_xyyyzz_1,   \
                             g_0_zz_0_xyyyzzz_0,  \
                             g_0_zz_0_xyyyzzz_1,  \
                             g_0_zz_0_xyyzzz_1,   \
                             g_0_zz_0_xyyzzzz_0,  \
                             g_0_zz_0_xyyzzzz_1,  \
                             g_0_zz_0_xyzzzz_1,   \
                             g_0_zz_0_xyzzzzz_0,  \
                             g_0_zz_0_xyzzzzz_1,  \
                             g_0_zz_0_xzzzzz_1,   \
                             g_0_zz_0_xzzzzzz_0,  \
                             g_0_zz_0_xzzzzzz_1,  \
                             g_0_zz_0_yyyyyy_1,   \
                             g_0_zz_0_yyyyyyy_0,  \
                             g_0_zz_0_yyyyyyy_1,  \
                             g_0_zz_0_yyyyyyz_0,  \
                             g_0_zz_0_yyyyyyz_1,  \
                             g_0_zz_0_yyyyyz_1,   \
                             g_0_zz_0_yyyyyzz_0,  \
                             g_0_zz_0_yyyyyzz_1,  \
                             g_0_zz_0_yyyyzz_1,   \
                             g_0_zz_0_yyyyzzz_0,  \
                             g_0_zz_0_yyyyzzz_1,  \
                             g_0_zz_0_yyyzzz_1,   \
                             g_0_zz_0_yyyzzzz_0,  \
                             g_0_zz_0_yyyzzzz_1,  \
                             g_0_zz_0_yyzzzz_1,   \
                             g_0_zz_0_yyzzzzz_0,  \
                             g_0_zz_0_yyzzzzz_1,  \
                             g_0_zz_0_yzzzzz_1,   \
                             g_0_zz_0_yzzzzzz_0,  \
                             g_0_zz_0_yzzzzzz_1,  \
                             g_0_zz_0_zzzzzz_1,   \
                             g_0_zz_0_zzzzzzz_0,  \
                             g_0_zz_0_zzzzzzz_1,  \
                             g_0_zzz_0_xxxxxxx_0, \
                             g_0_zzz_0_xxxxxxy_0, \
                             g_0_zzz_0_xxxxxxz_0, \
                             g_0_zzz_0_xxxxxyy_0, \
                             g_0_zzz_0_xxxxxyz_0, \
                             g_0_zzz_0_xxxxxzz_0, \
                             g_0_zzz_0_xxxxyyy_0, \
                             g_0_zzz_0_xxxxyyz_0, \
                             g_0_zzz_0_xxxxyzz_0, \
                             g_0_zzz_0_xxxxzzz_0, \
                             g_0_zzz_0_xxxyyyy_0, \
                             g_0_zzz_0_xxxyyyz_0, \
                             g_0_zzz_0_xxxyyzz_0, \
                             g_0_zzz_0_xxxyzzz_0, \
                             g_0_zzz_0_xxxzzzz_0, \
                             g_0_zzz_0_xxyyyyy_0, \
                             g_0_zzz_0_xxyyyyz_0, \
                             g_0_zzz_0_xxyyyzz_0, \
                             g_0_zzz_0_xxyyzzz_0, \
                             g_0_zzz_0_xxyzzzz_0, \
                             g_0_zzz_0_xxzzzzz_0, \
                             g_0_zzz_0_xyyyyyy_0, \
                             g_0_zzz_0_xyyyyyz_0, \
                             g_0_zzz_0_xyyyyzz_0, \
                             g_0_zzz_0_xyyyzzz_0, \
                             g_0_zzz_0_xyyzzzz_0, \
                             g_0_zzz_0_xyzzzzz_0, \
                             g_0_zzz_0_xzzzzzz_0, \
                             g_0_zzz_0_yyyyyyy_0, \
                             g_0_zzz_0_yyyyyyz_0, \
                             g_0_zzz_0_yyyyyzz_0, \
                             g_0_zzz_0_yyyyzzz_0, \
                             g_0_zzz_0_yyyzzzz_0, \
                             g_0_zzz_0_yyzzzzz_0, \
                             g_0_zzz_0_yzzzzzz_0, \
                             g_0_zzz_0_zzzzzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xxxxxxx_0[i] = 2.0 * g_0_z_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxx_1[i] * fti_ab_0 + g_0_zz_0_xxxxxxx_0[i] * pb_z +
                                 g_0_zz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxxy_0[i] = 2.0 * g_0_z_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxy_1[i] * fti_ab_0 + g_0_zz_0_xxxxxxy_0[i] * pb_z +
                                 g_0_zz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxxz_0[i] = 2.0 * g_0_z_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxxz_1[i] * fti_ab_0 + g_0_zz_0_xxxxxx_1[i] * fi_abcd_0 +
                                 g_0_zz_0_xxxxxxz_0[i] * pb_z + g_0_zz_0_xxxxxxz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxyy_0[i] = 2.0 * g_0_z_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxyy_1[i] * fti_ab_0 + g_0_zz_0_xxxxxyy_0[i] * pb_z +
                                 g_0_zz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxyz_0[i] = 2.0 * g_0_z_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxyz_1[i] * fti_ab_0 + g_0_zz_0_xxxxxy_1[i] * fi_abcd_0 +
                                 g_0_zz_0_xxxxxyz_0[i] * pb_z + g_0_zz_0_xxxxxyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxzz_0[i] = 2.0 * g_0_z_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxzz_0[i] * pb_z + g_0_zz_0_xxxxxzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyyy_0[i] = 2.0 * g_0_z_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyyy_1[i] * fti_ab_0 + g_0_zz_0_xxxxyyy_0[i] * pb_z +
                                 g_0_zz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyyz_0[i] = 2.0 * g_0_z_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyyz_1[i] * fti_ab_0 + g_0_zz_0_xxxxyy_1[i] * fi_abcd_0 +
                                 g_0_zz_0_xxxxyyz_0[i] * pb_z + g_0_zz_0_xxxxyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyzz_0[i] = 2.0 * g_0_z_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyzz_0[i] * pb_z + g_0_zz_0_xxxxyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxzzz_0[i] = 2.0 * g_0_z_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxzzz_0[i] * pb_z + g_0_zz_0_xxxxzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyyy_0[i] = 2.0 * g_0_z_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyyy_1[i] * fti_ab_0 + g_0_zz_0_xxxyyyy_0[i] * pb_z +
                                 g_0_zz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyyz_0[i] = 2.0 * g_0_z_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyyz_1[i] * fti_ab_0 + g_0_zz_0_xxxyyy_1[i] * fi_abcd_0 +
                                 g_0_zz_0_xxxyyyz_0[i] * pb_z + g_0_zz_0_xxxyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyzz_0[i] = 2.0 * g_0_z_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyzz_0[i] * pb_z + g_0_zz_0_xxxyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyzzz_0[i] = 2.0 * g_0_z_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzzz_0[i] * pb_z + g_0_zz_0_xxxyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxzzzz_0[i] = 2.0 * g_0_z_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxzzzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_zz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxzzzz_0[i] * pb_z + g_0_zz_0_xxxzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyyy_0[i] = 2.0 * g_0_z_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyyy_1[i] * fti_ab_0 + g_0_zz_0_xxyyyyy_0[i] * pb_z +
                                 g_0_zz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyyz_0[i] = 2.0 * g_0_z_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyyz_1[i] * fti_ab_0 + g_0_zz_0_xxyyyy_1[i] * fi_abcd_0 +
                                 g_0_zz_0_xxyyyyz_0[i] * pb_z + g_0_zz_0_xxyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyzz_0[i] = 2.0 * g_0_z_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyzz_0[i] * pb_z + g_0_zz_0_xxyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyzzz_0[i] = 2.0 * g_0_z_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzzz_0[i] * pb_z + g_0_zz_0_xxyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyzzzz_0[i] = 2.0 * g_0_z_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyzzzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_zz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzzz_0[i] * pb_z + g_0_zz_0_xxyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxzzzzz_0[i] = 2.0 * g_0_z_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxzzzzz_1[i] * fti_ab_0 +
                                 5.0 * g_0_zz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzzzzz_0[i] * pb_z + g_0_zz_0_xxzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyyy_0[i] = 2.0 * g_0_z_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyyy_1[i] * fti_ab_0 + g_0_zz_0_xyyyyyy_0[i] * pb_z +
                                 g_0_zz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyyz_0[i] = 2.0 * g_0_z_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyyz_1[i] * fti_ab_0 + g_0_zz_0_xyyyyy_1[i] * fi_abcd_0 +
                                 g_0_zz_0_xyyyyyz_0[i] * pb_z + g_0_zz_0_xyyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyzz_0[i] = 2.0 * g_0_z_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyzz_0[i] * pb_z + g_0_zz_0_xyyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyzzz_0[i] = 2.0 * g_0_z_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzzz_0[i] * pb_z + g_0_zz_0_xyyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyzzzz_0[i] = 2.0 * g_0_z_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyzzzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_zz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzzz_0[i] * pb_z + g_0_zz_0_xyyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyzzzzz_0[i] = 2.0 * g_0_z_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyzzzzz_1[i] * fti_ab_0 +
                                 5.0 * g_0_zz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzzz_0[i] * pb_z + g_0_zz_0_xyzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xzzzzzz_0[i] = 2.0 * g_0_z_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xzzzzzz_1[i] * fti_ab_0 +
                                 6.0 * g_0_zz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzzzzz_0[i] * pb_z + g_0_zz_0_xzzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyyy_0[i] = 2.0 * g_0_z_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyyy_1[i] * fti_ab_0 + g_0_zz_0_yyyyyyy_0[i] * pb_z +
                                 g_0_zz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyyz_0[i] = 2.0 * g_0_z_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyyz_1[i] * fti_ab_0 + g_0_zz_0_yyyyyy_1[i] * fi_abcd_0 +
                                 g_0_zz_0_yyyyyyz_0[i] * pb_z + g_0_zz_0_yyyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyzz_0[i] = 2.0 * g_0_z_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyzz_0[i] * pb_z + g_0_zz_0_yyyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyzzz_0[i] = 2.0 * g_0_z_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyzzz_0[i] * pb_z + g_0_zz_0_yyyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyzzzz_0[i] = 2.0 * g_0_z_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyzzzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_zz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyzzzz_0[i] * pb_z + g_0_zz_0_yyyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyzzzzz_0[i] = 2.0 * g_0_z_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyzzzzz_1[i] * fti_ab_0 +
                                 5.0 * g_0_zz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzzzzz_0[i] * pb_z + g_0_zz_0_yyzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yzzzzzz_0[i] = 2.0 * g_0_z_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yzzzzzz_1[i] * fti_ab_0 +
                                 6.0 * g_0_zz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzzzzz_0[i] * pb_z + g_0_zz_0_yzzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_zzzzzzz_0[i] = 2.0 * g_0_z_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zzzzzzz_1[i] * fti_ab_0 +
                                 7.0 * g_0_zz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zz_0_zzzzzzz_0[i] * pb_z + g_0_zz_0_zzzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
