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

#include "ThreeCenterElectronRepulsionPrimRecFSK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsk,
                                 size_t idx_eri_0_psk,
                                 size_t idx_eri_1_psk,
                                 size_t idx_eri_1_dsi,
                                 size_t idx_eri_1_dsk,
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

    /// Set up components of auxilary buffer : PSK

    auto g_x_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_psk);

    auto g_x_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_psk + 1);

    auto g_x_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_psk + 2);

    auto g_x_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_psk + 3);

    auto g_x_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_psk + 4);

    auto g_x_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_psk + 5);

    auto g_x_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_psk + 6);

    auto g_x_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_psk + 7);

    auto g_x_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_psk + 8);

    auto g_x_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_psk + 9);

    auto g_x_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_psk + 10);

    auto g_x_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_psk + 11);

    auto g_x_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_psk + 12);

    auto g_x_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_psk + 13);

    auto g_x_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_psk + 14);

    auto g_x_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_psk + 15);

    auto g_x_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_psk + 16);

    auto g_x_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_psk + 17);

    auto g_x_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_psk + 18);

    auto g_x_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_psk + 19);

    auto g_x_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_psk + 20);

    auto g_x_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 21);

    auto g_x_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 22);

    auto g_x_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 23);

    auto g_x_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 24);

    auto g_x_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 25);

    auto g_x_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 26);

    auto g_x_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 27);

    auto g_x_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 28);

    auto g_x_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 29);

    auto g_x_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 30);

    auto g_x_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 31);

    auto g_x_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 32);

    auto g_x_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 33);

    auto g_x_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 34);

    auto g_x_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 35);

    auto g_y_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_psk + 36);

    auto g_y_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_psk + 37);

    auto g_y_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_psk + 38);

    auto g_y_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_psk + 39);

    auto g_y_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_psk + 40);

    auto g_y_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_psk + 41);

    auto g_y_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_psk + 42);

    auto g_y_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_psk + 43);

    auto g_y_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_psk + 44);

    auto g_y_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_psk + 45);

    auto g_y_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_psk + 46);

    auto g_y_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_psk + 47);

    auto g_y_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_psk + 48);

    auto g_y_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_psk + 49);

    auto g_y_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_psk + 50);

    auto g_y_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_psk + 51);

    auto g_y_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_psk + 52);

    auto g_y_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_psk + 53);

    auto g_y_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_psk + 54);

    auto g_y_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_psk + 55);

    auto g_y_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_psk + 56);

    auto g_y_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 57);

    auto g_y_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 58);

    auto g_y_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 59);

    auto g_y_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 60);

    auto g_y_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 61);

    auto g_y_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 62);

    auto g_y_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 63);

    auto g_y_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 64);

    auto g_y_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 65);

    auto g_y_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 66);

    auto g_y_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 67);

    auto g_y_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 68);

    auto g_y_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 69);

    auto g_y_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 70);

    auto g_y_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 71);

    auto g_z_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_psk + 72);

    auto g_z_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_psk + 73);

    auto g_z_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_psk + 74);

    auto g_z_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_psk + 75);

    auto g_z_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_psk + 76);

    auto g_z_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_psk + 77);

    auto g_z_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_psk + 78);

    auto g_z_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_psk + 79);

    auto g_z_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_psk + 80);

    auto g_z_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_psk + 81);

    auto g_z_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_psk + 82);

    auto g_z_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_psk + 83);

    auto g_z_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_psk + 84);

    auto g_z_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_psk + 85);

    auto g_z_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_psk + 86);

    auto g_z_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_psk + 87);

    auto g_z_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_psk + 88);

    auto g_z_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_psk + 89);

    auto g_z_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_psk + 90);

    auto g_z_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_psk + 91);

    auto g_z_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_psk + 92);

    auto g_z_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 93);

    auto g_z_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 94);

    auto g_z_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 95);

    auto g_z_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 96);

    auto g_z_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 97);

    auto g_z_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 98);

    auto g_z_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 99);

    auto g_z_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 100);

    auto g_z_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 101);

    auto g_z_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 102);

    auto g_z_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 103);

    auto g_z_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 104);

    auto g_z_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 105);

    auto g_z_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 106);

    auto g_z_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 107);

    /// Set up components of auxilary buffer : PSK

    auto g_x_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_psk);

    auto g_x_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_psk + 1);

    auto g_x_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_psk + 2);

    auto g_x_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_psk + 3);

    auto g_x_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_psk + 4);

    auto g_x_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_psk + 5);

    auto g_x_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_psk + 6);

    auto g_x_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_psk + 7);

    auto g_x_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_psk + 8);

    auto g_x_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_psk + 9);

    auto g_x_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_psk + 10);

    auto g_x_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_psk + 11);

    auto g_x_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_psk + 12);

    auto g_x_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_psk + 13);

    auto g_x_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_psk + 14);

    auto g_x_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_psk + 15);

    auto g_x_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_psk + 16);

    auto g_x_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_psk + 17);

    auto g_x_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_psk + 18);

    auto g_x_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_psk + 19);

    auto g_x_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_psk + 20);

    auto g_x_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_psk + 21);

    auto g_x_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_psk + 22);

    auto g_x_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_psk + 23);

    auto g_x_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_psk + 24);

    auto g_x_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_psk + 25);

    auto g_x_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_psk + 26);

    auto g_x_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 27);

    auto g_x_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_psk + 28);

    auto g_x_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_psk + 29);

    auto g_x_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_psk + 30);

    auto g_x_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_psk + 31);

    auto g_x_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_psk + 32);

    auto g_x_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_psk + 33);

    auto g_x_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 34);

    auto g_x_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 35);

    auto g_y_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_psk + 36);

    auto g_y_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_psk + 37);

    auto g_y_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_psk + 38);

    auto g_y_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_psk + 39);

    auto g_y_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_psk + 40);

    auto g_y_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_psk + 41);

    auto g_y_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_psk + 42);

    auto g_y_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_psk + 43);

    auto g_y_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_psk + 44);

    auto g_y_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_psk + 45);

    auto g_y_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_psk + 46);

    auto g_y_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_psk + 47);

    auto g_y_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_psk + 48);

    auto g_y_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_psk + 49);

    auto g_y_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_psk + 50);

    auto g_y_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_psk + 51);

    auto g_y_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_psk + 52);

    auto g_y_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_psk + 53);

    auto g_y_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_psk + 54);

    auto g_y_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_psk + 55);

    auto g_y_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_psk + 56);

    auto g_y_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_psk + 57);

    auto g_y_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_psk + 58);

    auto g_y_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_psk + 59);

    auto g_y_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_psk + 60);

    auto g_y_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_psk + 61);

    auto g_y_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_psk + 62);

    auto g_y_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 63);

    auto g_y_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_psk + 64);

    auto g_y_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_psk + 65);

    auto g_y_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_psk + 66);

    auto g_y_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_psk + 67);

    auto g_y_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_psk + 68);

    auto g_y_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_psk + 69);

    auto g_y_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 70);

    auto g_y_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 71);

    auto g_z_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_psk + 72);

    auto g_z_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_psk + 73);

    auto g_z_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_psk + 74);

    auto g_z_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_psk + 75);

    auto g_z_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_psk + 76);

    auto g_z_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_psk + 77);

    auto g_z_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_psk + 78);

    auto g_z_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_psk + 79);

    auto g_z_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_psk + 80);

    auto g_z_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_psk + 81);

    auto g_z_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_psk + 82);

    auto g_z_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_psk + 83);

    auto g_z_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_psk + 84);

    auto g_z_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_psk + 85);

    auto g_z_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_psk + 86);

    auto g_z_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_psk + 87);

    auto g_z_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_psk + 88);

    auto g_z_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_psk + 89);

    auto g_z_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_psk + 90);

    auto g_z_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_psk + 91);

    auto g_z_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_psk + 92);

    auto g_z_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_psk + 93);

    auto g_z_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_psk + 94);

    auto g_z_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_psk + 95);

    auto g_z_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_psk + 96);

    auto g_z_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_psk + 97);

    auto g_z_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_psk + 98);

    auto g_z_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 99);

    auto g_z_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_psk + 100);

    auto g_z_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_psk + 101);

    auto g_z_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_psk + 102);

    auto g_z_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_psk + 103);

    auto g_z_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_psk + 104);

    auto g_z_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_psk + 105);

    auto g_z_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 106);

    auto g_z_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_psk + 107);

    /// Set up components of auxilary buffer : DSI

    auto g_xx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_dsi);

    auto g_xx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_dsi + 1);

    auto g_xx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_dsi + 2);

    auto g_xx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_dsi + 3);

    auto g_xx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_dsi + 4);

    auto g_xx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_dsi + 5);

    auto g_xx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_dsi + 6);

    auto g_xx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_dsi + 7);

    auto g_xx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_dsi + 8);

    auto g_xx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_dsi + 9);

    auto g_xx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_dsi + 10);

    auto g_xx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_dsi + 11);

    auto g_xx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_dsi + 12);

    auto g_xx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_dsi + 13);

    auto g_xx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_dsi + 14);

    auto g_xx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 15);

    auto g_xx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 16);

    auto g_xx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 17);

    auto g_xx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 18);

    auto g_xx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 19);

    auto g_xx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 20);

    auto g_xx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 21);

    auto g_xx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 22);

    auto g_xx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 23);

    auto g_xx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 24);

    auto g_xx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 25);

    auto g_xx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 26);

    auto g_xx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 27);

    auto g_yy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_dsi + 84);

    auto g_yy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_dsi + 85);

    auto g_yy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_dsi + 86);

    auto g_yy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_dsi + 87);

    auto g_yy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_dsi + 88);

    auto g_yy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_dsi + 89);

    auto g_yy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_dsi + 90);

    auto g_yy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_dsi + 91);

    auto g_yy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_dsi + 92);

    auto g_yy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_dsi + 93);

    auto g_yy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_dsi + 94);

    auto g_yy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_dsi + 95);

    auto g_yy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_dsi + 96);

    auto g_yy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_dsi + 97);

    auto g_yy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_dsi + 98);

    auto g_yy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 99);

    auto g_yy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 100);

    auto g_yy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 101);

    auto g_yy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 102);

    auto g_yy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 103);

    auto g_yy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 104);

    auto g_yy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 105);

    auto g_yy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 106);

    auto g_yy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 107);

    auto g_yy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 108);

    auto g_yy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 109);

    auto g_yy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 110);

    auto g_yy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 111);

    auto g_yz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_dsi + 116);

    auto g_yz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_dsi + 119);

    auto g_yz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_dsi + 120);

    auto g_yz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_dsi + 123);

    auto g_yz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_dsi + 124);

    auto g_yz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_dsi + 125);

    auto g_yz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 128);

    auto g_yz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 129);

    auto g_yz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 130);

    auto g_yz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 131);

    auto g_yz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 134);

    auto g_yz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 135);

    auto g_yz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 136);

    auto g_yz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 137);

    auto g_yz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 138);

    auto g_zz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_dsi + 140);

    auto g_zz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_dsi + 141);

    auto g_zz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_dsi + 142);

    auto g_zz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_dsi + 143);

    auto g_zz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_dsi + 144);

    auto g_zz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_dsi + 145);

    auto g_zz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_dsi + 146);

    auto g_zz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_dsi + 147);

    auto g_zz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_dsi + 148);

    auto g_zz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_dsi + 149);

    auto g_zz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_dsi + 150);

    auto g_zz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_dsi + 151);

    auto g_zz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_dsi + 152);

    auto g_zz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_dsi + 153);

    auto g_zz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_dsi + 154);

    auto g_zz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 155);

    auto g_zz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 156);

    auto g_zz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 157);

    auto g_zz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 158);

    auto g_zz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 159);

    auto g_zz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 160);

    auto g_zz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 161);

    auto g_zz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 162);

    auto g_zz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 163);

    auto g_zz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 164);

    auto g_zz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 165);

    auto g_zz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 166);

    auto g_zz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 167);

    /// Set up components of auxilary buffer : DSK

    auto g_xx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_dsk);

    auto g_xx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_dsk + 1);

    auto g_xx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_dsk + 2);

    auto g_xx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_dsk + 3);

    auto g_xx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_dsk + 4);

    auto g_xx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_dsk + 5);

    auto g_xx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_dsk + 6);

    auto g_xx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_dsk + 7);

    auto g_xx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_dsk + 8);

    auto g_xx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_dsk + 9);

    auto g_xx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_dsk + 10);

    auto g_xx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_dsk + 11);

    auto g_xx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_dsk + 12);

    auto g_xx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_dsk + 13);

    auto g_xx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_dsk + 14);

    auto g_xx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 15);

    auto g_xx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 16);

    auto g_xx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 17);

    auto g_xx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 18);

    auto g_xx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 19);

    auto g_xx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 20);

    auto g_xx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 21);

    auto g_xx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 22);

    auto g_xx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 23);

    auto g_xx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 24);

    auto g_xx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 25);

    auto g_xx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 26);

    auto g_xx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 27);

    auto g_xx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 28);

    auto g_xx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 29);

    auto g_xx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 30);

    auto g_xx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 31);

    auto g_xx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 32);

    auto g_xx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 33);

    auto g_xx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 34);

    auto g_xx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 35);

    auto g_xy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_dsk + 37);

    auto g_xy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_dsk + 39);

    auto g_xy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_dsk + 42);

    auto g_xy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_dsk + 46);

    auto g_xy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 51);

    auto g_xy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 57);

    auto g_xz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_dsk + 72);

    auto g_xz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_dsk + 74);

    auto g_xz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_dsk + 77);

    auto g_xz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_dsk + 81);

    auto g_xz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_dsk + 86);

    auto g_xz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 92);

    auto g_xz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 99);

    auto g_yy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_dsk + 108);

    auto g_yy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_dsk + 109);

    auto g_yy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_dsk + 110);

    auto g_yy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_dsk + 111);

    auto g_yy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_dsk + 112);

    auto g_yy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_dsk + 113);

    auto g_yy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_dsk + 114);

    auto g_yy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_dsk + 115);

    auto g_yy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_dsk + 116);

    auto g_yy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_dsk + 117);

    auto g_yy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_dsk + 118);

    auto g_yy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_dsk + 119);

    auto g_yy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_dsk + 120);

    auto g_yy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_dsk + 121);

    auto g_yy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_dsk + 122);

    auto g_yy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 123);

    auto g_yy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 124);

    auto g_yy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 125);

    auto g_yy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 126);

    auto g_yy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 127);

    auto g_yy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 128);

    auto g_yy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 129);

    auto g_yy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 130);

    auto g_yy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 131);

    auto g_yy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 132);

    auto g_yy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 133);

    auto g_yy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 134);

    auto g_yy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 135);

    auto g_yy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 136);

    auto g_yy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 137);

    auto g_yy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 138);

    auto g_yy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 139);

    auto g_yy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 140);

    auto g_yy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 141);

    auto g_yy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 142);

    auto g_yy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 143);

    auto g_yz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_dsk + 148);

    auto g_yz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_dsk + 151);

    auto g_yz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_dsk + 152);

    auto g_yz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_dsk + 155);

    auto g_yz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_dsk + 156);

    auto g_yz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_dsk + 157);

    auto g_yz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 160);

    auto g_yz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 161);

    auto g_yz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 162);

    auto g_yz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 163);

    auto g_yz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 166);

    auto g_yz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 167);

    auto g_yz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 168);

    auto g_yz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 169);

    auto g_yz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 170);

    auto g_yz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 172);

    auto g_yz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 173);

    auto g_yz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 174);

    auto g_yz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 175);

    auto g_yz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 176);

    auto g_yz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 177);

    auto g_yz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 178);

    auto g_yz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 179);

    auto g_zz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_dsk + 180);

    auto g_zz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_dsk + 181);

    auto g_zz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_dsk + 182);

    auto g_zz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_dsk + 183);

    auto g_zz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_dsk + 184);

    auto g_zz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_dsk + 185);

    auto g_zz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_dsk + 186);

    auto g_zz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_dsk + 187);

    auto g_zz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_dsk + 188);

    auto g_zz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_dsk + 189);

    auto g_zz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_dsk + 190);

    auto g_zz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_dsk + 191);

    auto g_zz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_dsk + 192);

    auto g_zz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_dsk + 193);

    auto g_zz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_dsk + 194);

    auto g_zz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 195);

    auto g_zz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 196);

    auto g_zz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 197);

    auto g_zz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 198);

    auto g_zz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 199);

    auto g_zz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 200);

    auto g_zz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 201);

    auto g_zz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 202);

    auto g_zz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 203);

    auto g_zz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 204);

    auto g_zz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 205);

    auto g_zz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 206);

    auto g_zz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 207);

    auto g_zz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 208);

    auto g_zz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 209);

    auto g_zz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 210);

    auto g_zz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 211);

    auto g_zz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 212);

    auto g_zz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 213);

    auto g_zz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 214);

    auto g_zz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 215);

    /// Set up 0-36 components of targeted buffer : FSK

    auto g_xxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk);

    auto g_xxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 1);

    auto g_xxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 2);

    auto g_xxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 3);

    auto g_xxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 4);

    auto g_xxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 5);

    auto g_xxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 6);

    auto g_xxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 7);

    auto g_xxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 8);

    auto g_xxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 9);

    auto g_xxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 10);

    auto g_xxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 11);

    auto g_xxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 12);

    auto g_xxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 13);

    auto g_xxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 14);

    auto g_xxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 15);

    auto g_xxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 16);

    auto g_xxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 17);

    auto g_xxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 18);

    auto g_xxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 19);

    auto g_xxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 20);

    auto g_xxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 21);

    auto g_xxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 22);

    auto g_xxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 23);

    auto g_xxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 24);

    auto g_xxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 25);

    auto g_xxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 26);

    auto g_xxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 27);

    auto g_xxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 28);

    auto g_xxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 29);

    auto g_xxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 30);

    auto g_xxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 31);

    auto g_xxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 32);

    auto g_xxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 33);

    auto g_xxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 34);

    auto g_xxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 35);

    #pragma omp simd aligned(g_x_0_xxxxxxx_0, g_x_0_xxxxxxx_1, g_x_0_xxxxxxy_0, g_x_0_xxxxxxy_1, g_x_0_xxxxxxz_0, g_x_0_xxxxxxz_1, g_x_0_xxxxxyy_0, g_x_0_xxxxxyy_1, g_x_0_xxxxxyz_0, g_x_0_xxxxxyz_1, g_x_0_xxxxxzz_0, g_x_0_xxxxxzz_1, g_x_0_xxxxyyy_0, g_x_0_xxxxyyy_1, g_x_0_xxxxyyz_0, g_x_0_xxxxyyz_1, g_x_0_xxxxyzz_0, g_x_0_xxxxyzz_1, g_x_0_xxxxzzz_0, g_x_0_xxxxzzz_1, g_x_0_xxxyyyy_0, g_x_0_xxxyyyy_1, g_x_0_xxxyyyz_0, g_x_0_xxxyyyz_1, g_x_0_xxxyyzz_0, g_x_0_xxxyyzz_1, g_x_0_xxxyzzz_0, g_x_0_xxxyzzz_1, g_x_0_xxxzzzz_0, g_x_0_xxxzzzz_1, g_x_0_xxyyyyy_0, g_x_0_xxyyyyy_1, g_x_0_xxyyyyz_0, g_x_0_xxyyyyz_1, g_x_0_xxyyyzz_0, g_x_0_xxyyyzz_1, g_x_0_xxyyzzz_0, g_x_0_xxyyzzz_1, g_x_0_xxyzzzz_0, g_x_0_xxyzzzz_1, g_x_0_xxzzzzz_0, g_x_0_xxzzzzz_1, g_x_0_xyyyyyy_0, g_x_0_xyyyyyy_1, g_x_0_xyyyyyz_0, g_x_0_xyyyyyz_1, g_x_0_xyyyyzz_0, g_x_0_xyyyyzz_1, g_x_0_xyyyzzz_0, g_x_0_xyyyzzz_1, g_x_0_xyyzzzz_0, g_x_0_xyyzzzz_1, g_x_0_xyzzzzz_0, g_x_0_xyzzzzz_1, g_x_0_xzzzzzz_0, g_x_0_xzzzzzz_1, g_x_0_yyyyyyy_0, g_x_0_yyyyyyy_1, g_x_0_yyyyyyz_0, g_x_0_yyyyyyz_1, g_x_0_yyyyyzz_0, g_x_0_yyyyyzz_1, g_x_0_yyyyzzz_0, g_x_0_yyyyzzz_1, g_x_0_yyyzzzz_0, g_x_0_yyyzzzz_1, g_x_0_yyzzzzz_0, g_x_0_yyzzzzz_1, g_x_0_yzzzzzz_0, g_x_0_yzzzzzz_1, g_x_0_zzzzzzz_0, g_x_0_zzzzzzz_1, g_xx_0_xxxxxx_1, g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxy_1, g_xx_0_xxxxxyy_1, g_xx_0_xxxxxyz_1, g_xx_0_xxxxxz_1, g_xx_0_xxxxxzz_1, g_xx_0_xxxxyy_1, g_xx_0_xxxxyyy_1, g_xx_0_xxxxyyz_1, g_xx_0_xxxxyz_1, g_xx_0_xxxxyzz_1, g_xx_0_xxxxzz_1, g_xx_0_xxxxzzz_1, g_xx_0_xxxyyy_1, g_xx_0_xxxyyyy_1, g_xx_0_xxxyyyz_1, g_xx_0_xxxyyz_1, g_xx_0_xxxyyzz_1, g_xx_0_xxxyzz_1, g_xx_0_xxxyzzz_1, g_xx_0_xxxzzz_1, g_xx_0_xxxzzzz_1, g_xx_0_xxyyyy_1, g_xx_0_xxyyyyy_1, g_xx_0_xxyyyyz_1, g_xx_0_xxyyyz_1, g_xx_0_xxyyyzz_1, g_xx_0_xxyyzz_1, g_xx_0_xxyyzzz_1, g_xx_0_xxyzzz_1, g_xx_0_xxyzzzz_1, g_xx_0_xxzzzz_1, g_xx_0_xxzzzzz_1, g_xx_0_xyyyyy_1, g_xx_0_xyyyyyy_1, g_xx_0_xyyyyyz_1, g_xx_0_xyyyyz_1, g_xx_0_xyyyyzz_1, g_xx_0_xyyyzz_1, g_xx_0_xyyyzzz_1, g_xx_0_xyyzzz_1, g_xx_0_xyyzzzz_1, g_xx_0_xyzzzz_1, g_xx_0_xyzzzzz_1, g_xx_0_xzzzzz_1, g_xx_0_xzzzzzz_1, g_xx_0_yyyyyy_1, g_xx_0_yyyyyyy_1, g_xx_0_yyyyyyz_1, g_xx_0_yyyyyz_1, g_xx_0_yyyyyzz_1, g_xx_0_yyyyzz_1, g_xx_0_yyyyzzz_1, g_xx_0_yyyzzz_1, g_xx_0_yyyzzzz_1, g_xx_0_yyzzzz_1, g_xx_0_yyzzzzz_1, g_xx_0_yzzzzz_1, g_xx_0_yzzzzzz_1, g_xx_0_zzzzzz_1, g_xx_0_zzzzzzz_1, g_xxx_0_xxxxxxx_0, g_xxx_0_xxxxxxy_0, g_xxx_0_xxxxxxz_0, g_xxx_0_xxxxxyy_0, g_xxx_0_xxxxxyz_0, g_xxx_0_xxxxxzz_0, g_xxx_0_xxxxyyy_0, g_xxx_0_xxxxyyz_0, g_xxx_0_xxxxyzz_0, g_xxx_0_xxxxzzz_0, g_xxx_0_xxxyyyy_0, g_xxx_0_xxxyyyz_0, g_xxx_0_xxxyyzz_0, g_xxx_0_xxxyzzz_0, g_xxx_0_xxxzzzz_0, g_xxx_0_xxyyyyy_0, g_xxx_0_xxyyyyz_0, g_xxx_0_xxyyyzz_0, g_xxx_0_xxyyzzz_0, g_xxx_0_xxyzzzz_0, g_xxx_0_xxzzzzz_0, g_xxx_0_xyyyyyy_0, g_xxx_0_xyyyyyz_0, g_xxx_0_xyyyyzz_0, g_xxx_0_xyyyzzz_0, g_xxx_0_xyyzzzz_0, g_xxx_0_xyzzzzz_0, g_xxx_0_xzzzzzz_0, g_xxx_0_yyyyyyy_0, g_xxx_0_yyyyyyz_0, g_xxx_0_yyyyyzz_0, g_xxx_0_yyyyzzz_0, g_xxx_0_yyyzzzz_0, g_xxx_0_yyzzzzz_0, g_xxx_0_yzzzzzz_0, g_xxx_0_zzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xxxxxxx_0[i] = 2.0 * g_x_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxx_1[i] * fz_be_0 + 7.0 * g_xx_0_xxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxx_1[i] * wa_x[i];

        g_xxx_0_xxxxxxy_0[i] = 2.0 * g_x_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxxy_1[i] * wa_x[i];

        g_xxx_0_xxxxxxz_0[i] = 2.0 * g_x_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxxz_1[i] * wa_x[i];

        g_xxx_0_xxxxxyy_0[i] = 2.0 * g_x_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxxyy_1[i] * wa_x[i];

        g_xxx_0_xxxxxyz_0[i] = 2.0 * g_x_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxxyz_1[i] * wa_x[i];

        g_xxx_0_xxxxxzz_0[i] = 2.0 * g_x_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxxzz_1[i] * wa_x[i];

        g_xxx_0_xxxxyyy_0[i] = 2.0 * g_x_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyy_1[i] * wa_x[i];

        g_xxx_0_xxxxyyz_0[i] = 2.0 * g_x_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyz_1[i] * wa_x[i];

        g_xxx_0_xxxxyzz_0[i] = 2.0 * g_x_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzz_1[i] * wa_x[i];

        g_xxx_0_xxxxzzz_0[i] = 2.0 * g_x_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxxzzz_1[i] * wa_x[i];

        g_xxx_0_xxxyyyy_0[i] = 2.0 * g_x_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyy_1[i] * wa_x[i];

        g_xxx_0_xxxyyyz_0[i] = 2.0 * g_x_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyz_1[i] * wa_x[i];

        g_xxx_0_xxxyyzz_0[i] = 2.0 * g_x_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzz_1[i] * wa_x[i];

        g_xxx_0_xxxyzzz_0[i] = 2.0 * g_x_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzz_1[i] * wa_x[i];

        g_xxx_0_xxxzzzz_0[i] = 2.0 * g_x_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxzzzz_1[i] * fi_acd_0 + g_xx_0_xxxzzzz_1[i] * wa_x[i];

        g_xxx_0_xxyyyyy_0[i] = 2.0 * g_x_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyy_1[i] * wa_x[i];

        g_xxx_0_xxyyyyz_0[i] = 2.0 * g_x_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyz_1[i] * wa_x[i];

        g_xxx_0_xxyyyzz_0[i] = 2.0 * g_x_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzz_1[i] * wa_x[i];

        g_xxx_0_xxyyzzz_0[i] = 2.0 * g_x_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzz_1[i] * wa_x[i];

        g_xxx_0_xxyzzzz_0[i] = 2.0 * g_x_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzz_1[i] * wa_x[i];

        g_xxx_0_xxzzzzz_0[i] = 2.0 * g_x_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xzzzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyyyy_0[i] = 2.0 * g_x_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyy_1[i] * wa_x[i];

        g_xxx_0_xyyyyyz_0[i] = 2.0 * g_x_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyz_1[i] * wa_x[i];

        g_xxx_0_xyyyyzz_0[i] = 2.0 * g_x_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyzz_1[i] * fz_be_0 + g_xx_0_yyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzz_1[i] * wa_x[i];

        g_xxx_0_xyyyzzz_0[i] = 2.0 * g_x_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyzzz_1[i] * fz_be_0 + g_xx_0_yyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzz_1[i] * wa_x[i];

        g_xxx_0_xyyzzzz_0[i] = 2.0 * g_x_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyzzzz_1[i] * fz_be_0 + g_xx_0_yyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzz_1[i] * wa_x[i];

        g_xxx_0_xyzzzzz_0[i] = 2.0 * g_x_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyzzzzz_1[i] * fz_be_0 + g_xx_0_yzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzz_1[i] * wa_x[i];

        g_xxx_0_xzzzzzz_0[i] = 2.0 * g_x_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xzzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyyyy_0[i] = 2.0 * g_x_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyyyy_1[i] * wa_x[i];

        g_xxx_0_yyyyyyz_0[i] = 2.0 * g_x_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyyyz_1[i] * wa_x[i];

        g_xxx_0_yyyyyzz_0[i] = 2.0 * g_x_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyzz_1[i] * fz_be_0 + g_xx_0_yyyyyzz_1[i] * wa_x[i];

        g_xxx_0_yyyyzzz_0[i] = 2.0 * g_x_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyzzz_1[i] * fz_be_0 + g_xx_0_yyyyzzz_1[i] * wa_x[i];

        g_xxx_0_yyyzzzz_0[i] = 2.0 * g_x_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyzzzz_1[i] * fz_be_0 + g_xx_0_yyyzzzz_1[i] * wa_x[i];

        g_xxx_0_yyzzzzz_0[i] = 2.0 * g_x_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyzzzzz_1[i] * fz_be_0 + g_xx_0_yyzzzzz_1[i] * wa_x[i];

        g_xxx_0_yzzzzzz_0[i] = 2.0 * g_x_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yzzzzzz_1[i] * fz_be_0 + g_xx_0_yzzzzzz_1[i] * wa_x[i];

        g_xxx_0_zzzzzzz_0[i] = 2.0 * g_x_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_zzzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 36-72 components of targeted buffer : FSK

    auto g_xxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 36);

    auto g_xxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 37);

    auto g_xxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 38);

    auto g_xxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 39);

    auto g_xxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 40);

    auto g_xxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 41);

    auto g_xxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 42);

    auto g_xxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 43);

    auto g_xxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 44);

    auto g_xxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 45);

    auto g_xxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 46);

    auto g_xxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 47);

    auto g_xxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 48);

    auto g_xxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 49);

    auto g_xxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 50);

    auto g_xxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 51);

    auto g_xxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 52);

    auto g_xxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 53);

    auto g_xxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 54);

    auto g_xxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 55);

    auto g_xxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 56);

    auto g_xxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 57);

    auto g_xxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 58);

    auto g_xxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 59);

    auto g_xxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 60);

    auto g_xxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 61);

    auto g_xxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 62);

    auto g_xxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 63);

    auto g_xxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 64);

    auto g_xxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 65);

    auto g_xxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 66);

    auto g_xxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 67);

    auto g_xxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 68);

    auto g_xxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 69);

    auto g_xxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 70);

    auto g_xxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 71);

    #pragma omp simd aligned(g_xx_0_xxxxxx_1, g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxy_1, g_xx_0_xxxxxyy_1, g_xx_0_xxxxxyz_1, g_xx_0_xxxxxz_1, g_xx_0_xxxxxzz_1, g_xx_0_xxxxyy_1, g_xx_0_xxxxyyy_1, g_xx_0_xxxxyyz_1, g_xx_0_xxxxyz_1, g_xx_0_xxxxyzz_1, g_xx_0_xxxxzz_1, g_xx_0_xxxxzzz_1, g_xx_0_xxxyyy_1, g_xx_0_xxxyyyy_1, g_xx_0_xxxyyyz_1, g_xx_0_xxxyyz_1, g_xx_0_xxxyyzz_1, g_xx_0_xxxyzz_1, g_xx_0_xxxyzzz_1, g_xx_0_xxxzzz_1, g_xx_0_xxxzzzz_1, g_xx_0_xxyyyy_1, g_xx_0_xxyyyyy_1, g_xx_0_xxyyyyz_1, g_xx_0_xxyyyz_1, g_xx_0_xxyyyzz_1, g_xx_0_xxyyzz_1, g_xx_0_xxyyzzz_1, g_xx_0_xxyzzz_1, g_xx_0_xxyzzzz_1, g_xx_0_xxzzzz_1, g_xx_0_xxzzzzz_1, g_xx_0_xyyyyy_1, g_xx_0_xyyyyyy_1, g_xx_0_xyyyyyz_1, g_xx_0_xyyyyz_1, g_xx_0_xyyyyzz_1, g_xx_0_xyyyzz_1, g_xx_0_xyyyzzz_1, g_xx_0_xyyzzz_1, g_xx_0_xyyzzzz_1, g_xx_0_xyzzzz_1, g_xx_0_xyzzzzz_1, g_xx_0_xzzzzz_1, g_xx_0_xzzzzzz_1, g_xx_0_yyyyyy_1, g_xx_0_yyyyyyy_1, g_xx_0_yyyyyyz_1, g_xx_0_yyyyyz_1, g_xx_0_yyyyyzz_1, g_xx_0_yyyyzz_1, g_xx_0_yyyyzzz_1, g_xx_0_yyyzzz_1, g_xx_0_yyyzzzz_1, g_xx_0_yyzzzz_1, g_xx_0_yyzzzzz_1, g_xx_0_yzzzzz_1, g_xx_0_yzzzzzz_1, g_xx_0_zzzzzz_1, g_xx_0_zzzzzzz_1, g_xxy_0_xxxxxxx_0, g_xxy_0_xxxxxxy_0, g_xxy_0_xxxxxxz_0, g_xxy_0_xxxxxyy_0, g_xxy_0_xxxxxyz_0, g_xxy_0_xxxxxzz_0, g_xxy_0_xxxxyyy_0, g_xxy_0_xxxxyyz_0, g_xxy_0_xxxxyzz_0, g_xxy_0_xxxxzzz_0, g_xxy_0_xxxyyyy_0, g_xxy_0_xxxyyyz_0, g_xxy_0_xxxyyzz_0, g_xxy_0_xxxyzzz_0, g_xxy_0_xxxzzzz_0, g_xxy_0_xxyyyyy_0, g_xxy_0_xxyyyyz_0, g_xxy_0_xxyyyzz_0, g_xxy_0_xxyyzzz_0, g_xxy_0_xxyzzzz_0, g_xxy_0_xxzzzzz_0, g_xxy_0_xyyyyyy_0, g_xxy_0_xyyyyyz_0, g_xxy_0_xyyyyzz_0, g_xxy_0_xyyyzzz_0, g_xxy_0_xyyzzzz_0, g_xxy_0_xyzzzzz_0, g_xxy_0_xzzzzzz_0, g_xxy_0_yyyyyyy_0, g_xxy_0_yyyyyyz_0, g_xxy_0_yyyyyzz_0, g_xxy_0_yyyyzzz_0, g_xxy_0_yyyzzzz_0, g_xxy_0_yyzzzzz_0, g_xxy_0_yzzzzzz_0, g_xxy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xxxxxxx_0[i] = g_xx_0_xxxxxxx_1[i] * wa_y[i];

        g_xxy_0_xxxxxxy_0[i] = g_xx_0_xxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxy_1[i] * wa_y[i];

        g_xxy_0_xxxxxxz_0[i] = g_xx_0_xxxxxxz_1[i] * wa_y[i];

        g_xxy_0_xxxxxyy_0[i] = 2.0 * g_xx_0_xxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxyy_1[i] * wa_y[i];

        g_xxy_0_xxxxxyz_0[i] = g_xx_0_xxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxyz_1[i] * wa_y[i];

        g_xxy_0_xxxxxzz_0[i] = g_xx_0_xxxxxzz_1[i] * wa_y[i];

        g_xxy_0_xxxxyyy_0[i] = 3.0 * g_xx_0_xxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyy_1[i] * wa_y[i];

        g_xxy_0_xxxxyyz_0[i] = 2.0 * g_xx_0_xxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyz_1[i] * wa_y[i];

        g_xxy_0_xxxxyzz_0[i] = g_xx_0_xxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzz_1[i] * wa_y[i];

        g_xxy_0_xxxxzzz_0[i] = g_xx_0_xxxxzzz_1[i] * wa_y[i];

        g_xxy_0_xxxyyyy_0[i] = 4.0 * g_xx_0_xxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyy_1[i] * wa_y[i];

        g_xxy_0_xxxyyyz_0[i] = 3.0 * g_xx_0_xxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyz_1[i] * wa_y[i];

        g_xxy_0_xxxyyzz_0[i] = 2.0 * g_xx_0_xxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzz_1[i] * wa_y[i];

        g_xxy_0_xxxyzzz_0[i] = g_xx_0_xxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzz_1[i] * wa_y[i];

        g_xxy_0_xxxzzzz_0[i] = g_xx_0_xxxzzzz_1[i] * wa_y[i];

        g_xxy_0_xxyyyyy_0[i] = 5.0 * g_xx_0_xxyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyy_1[i] * wa_y[i];

        g_xxy_0_xxyyyyz_0[i] = 4.0 * g_xx_0_xxyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyz_1[i] * wa_y[i];

        g_xxy_0_xxyyyzz_0[i] = 3.0 * g_xx_0_xxyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzz_1[i] * wa_y[i];

        g_xxy_0_xxyyzzz_0[i] = 2.0 * g_xx_0_xxyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzz_1[i] * wa_y[i];

        g_xxy_0_xxyzzzz_0[i] = g_xx_0_xxzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzz_1[i] * wa_y[i];

        g_xxy_0_xxzzzzz_0[i] = g_xx_0_xxzzzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyyyy_0[i] = 6.0 * g_xx_0_xyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyy_1[i] * wa_y[i];

        g_xxy_0_xyyyyyz_0[i] = 5.0 * g_xx_0_xyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyz_1[i] * wa_y[i];

        g_xxy_0_xyyyyzz_0[i] = 4.0 * g_xx_0_xyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzz_1[i] * wa_y[i];

        g_xxy_0_xyyyzzz_0[i] = 3.0 * g_xx_0_xyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzz_1[i] * wa_y[i];

        g_xxy_0_xyyzzzz_0[i] = 2.0 * g_xx_0_xyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzz_1[i] * wa_y[i];

        g_xxy_0_xyzzzzz_0[i] = g_xx_0_xzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzz_1[i] * wa_y[i];

        g_xxy_0_xzzzzzz_0[i] = g_xx_0_xzzzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyyyy_0[i] = 7.0 * g_xx_0_yyyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyyy_1[i] * wa_y[i];

        g_xxy_0_yyyyyyz_0[i] = 6.0 * g_xx_0_yyyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyyyz_1[i] * wa_y[i];

        g_xxy_0_yyyyyzz_0[i] = 5.0 * g_xx_0_yyyyzz_1[i] * fi_acd_0 + g_xx_0_yyyyyzz_1[i] * wa_y[i];

        g_xxy_0_yyyyzzz_0[i] = 4.0 * g_xx_0_yyyzzz_1[i] * fi_acd_0 + g_xx_0_yyyyzzz_1[i] * wa_y[i];

        g_xxy_0_yyyzzzz_0[i] = 3.0 * g_xx_0_yyzzzz_1[i] * fi_acd_0 + g_xx_0_yyyzzzz_1[i] * wa_y[i];

        g_xxy_0_yyzzzzz_0[i] = 2.0 * g_xx_0_yzzzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzzz_1[i] * wa_y[i];

        g_xxy_0_yzzzzzz_0[i] = g_xx_0_zzzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzzz_1[i] * wa_y[i];

        g_xxy_0_zzzzzzz_0[i] = g_xx_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 72-108 components of targeted buffer : FSK

    auto g_xxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 72);

    auto g_xxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 73);

    auto g_xxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 74);

    auto g_xxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 75);

    auto g_xxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 76);

    auto g_xxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 77);

    auto g_xxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 78);

    auto g_xxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 79);

    auto g_xxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 80);

    auto g_xxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 81);

    auto g_xxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 82);

    auto g_xxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 83);

    auto g_xxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 84);

    auto g_xxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 85);

    auto g_xxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 86);

    auto g_xxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 87);

    auto g_xxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 88);

    auto g_xxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 89);

    auto g_xxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 90);

    auto g_xxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 91);

    auto g_xxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 92);

    auto g_xxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 93);

    auto g_xxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 94);

    auto g_xxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 95);

    auto g_xxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 96);

    auto g_xxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 97);

    auto g_xxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 98);

    auto g_xxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 99);

    auto g_xxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 100);

    auto g_xxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 101);

    auto g_xxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 102);

    auto g_xxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 103);

    auto g_xxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 104);

    auto g_xxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 105);

    auto g_xxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 106);

    auto g_xxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 107);

    #pragma omp simd aligned(g_xx_0_xxxxxx_1, g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxy_1, g_xx_0_xxxxxyy_1, g_xx_0_xxxxxyz_1, g_xx_0_xxxxxz_1, g_xx_0_xxxxxzz_1, g_xx_0_xxxxyy_1, g_xx_0_xxxxyyy_1, g_xx_0_xxxxyyz_1, g_xx_0_xxxxyz_1, g_xx_0_xxxxyzz_1, g_xx_0_xxxxzz_1, g_xx_0_xxxxzzz_1, g_xx_0_xxxyyy_1, g_xx_0_xxxyyyy_1, g_xx_0_xxxyyyz_1, g_xx_0_xxxyyz_1, g_xx_0_xxxyyzz_1, g_xx_0_xxxyzz_1, g_xx_0_xxxyzzz_1, g_xx_0_xxxzzz_1, g_xx_0_xxxzzzz_1, g_xx_0_xxyyyy_1, g_xx_0_xxyyyyy_1, g_xx_0_xxyyyyz_1, g_xx_0_xxyyyz_1, g_xx_0_xxyyyzz_1, g_xx_0_xxyyzz_1, g_xx_0_xxyyzzz_1, g_xx_0_xxyzzz_1, g_xx_0_xxyzzzz_1, g_xx_0_xxzzzz_1, g_xx_0_xxzzzzz_1, g_xx_0_xyyyyy_1, g_xx_0_xyyyyyy_1, g_xx_0_xyyyyyz_1, g_xx_0_xyyyyz_1, g_xx_0_xyyyyzz_1, g_xx_0_xyyyzz_1, g_xx_0_xyyyzzz_1, g_xx_0_xyyzzz_1, g_xx_0_xyyzzzz_1, g_xx_0_xyzzzz_1, g_xx_0_xyzzzzz_1, g_xx_0_xzzzzz_1, g_xx_0_xzzzzzz_1, g_xx_0_yyyyyy_1, g_xx_0_yyyyyyy_1, g_xx_0_yyyyyyz_1, g_xx_0_yyyyyz_1, g_xx_0_yyyyyzz_1, g_xx_0_yyyyzz_1, g_xx_0_yyyyzzz_1, g_xx_0_yyyzzz_1, g_xx_0_yyyzzzz_1, g_xx_0_yyzzzz_1, g_xx_0_yyzzzzz_1, g_xx_0_yzzzzz_1, g_xx_0_yzzzzzz_1, g_xx_0_zzzzzz_1, g_xx_0_zzzzzzz_1, g_xxz_0_xxxxxxx_0, g_xxz_0_xxxxxxy_0, g_xxz_0_xxxxxxz_0, g_xxz_0_xxxxxyy_0, g_xxz_0_xxxxxyz_0, g_xxz_0_xxxxxzz_0, g_xxz_0_xxxxyyy_0, g_xxz_0_xxxxyyz_0, g_xxz_0_xxxxyzz_0, g_xxz_0_xxxxzzz_0, g_xxz_0_xxxyyyy_0, g_xxz_0_xxxyyyz_0, g_xxz_0_xxxyyzz_0, g_xxz_0_xxxyzzz_0, g_xxz_0_xxxzzzz_0, g_xxz_0_xxyyyyy_0, g_xxz_0_xxyyyyz_0, g_xxz_0_xxyyyzz_0, g_xxz_0_xxyyzzz_0, g_xxz_0_xxyzzzz_0, g_xxz_0_xxzzzzz_0, g_xxz_0_xyyyyyy_0, g_xxz_0_xyyyyyz_0, g_xxz_0_xyyyyzz_0, g_xxz_0_xyyyzzz_0, g_xxz_0_xyyzzzz_0, g_xxz_0_xyzzzzz_0, g_xxz_0_xzzzzzz_0, g_xxz_0_yyyyyyy_0, g_xxz_0_yyyyyyz_0, g_xxz_0_yyyyyzz_0, g_xxz_0_yyyyzzz_0, g_xxz_0_yyyzzzz_0, g_xxz_0_yyzzzzz_0, g_xxz_0_yzzzzzz_0, g_xxz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xxxxxxx_0[i] = g_xx_0_xxxxxxx_1[i] * wa_z[i];

        g_xxz_0_xxxxxxy_0[i] = g_xx_0_xxxxxxy_1[i] * wa_z[i];

        g_xxz_0_xxxxxxz_0[i] = g_xx_0_xxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxz_1[i] * wa_z[i];

        g_xxz_0_xxxxxyy_0[i] = g_xx_0_xxxxxyy_1[i] * wa_z[i];

        g_xxz_0_xxxxxyz_0[i] = g_xx_0_xxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxyz_1[i] * wa_z[i];

        g_xxz_0_xxxxxzz_0[i] = 2.0 * g_xx_0_xxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxzz_1[i] * wa_z[i];

        g_xxz_0_xxxxyyy_0[i] = g_xx_0_xxxxyyy_1[i] * wa_z[i];

        g_xxz_0_xxxxyyz_0[i] = g_xx_0_xxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyz_1[i] * wa_z[i];

        g_xxz_0_xxxxyzz_0[i] = 2.0 * g_xx_0_xxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxyzz_1[i] * wa_z[i];

        g_xxz_0_xxxxzzz_0[i] = 3.0 * g_xx_0_xxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxzzz_1[i] * wa_z[i];

        g_xxz_0_xxxyyyy_0[i] = g_xx_0_xxxyyyy_1[i] * wa_z[i];

        g_xxz_0_xxxyyyz_0[i] = g_xx_0_xxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyz_1[i] * wa_z[i];

        g_xxz_0_xxxyyzz_0[i] = 2.0 * g_xx_0_xxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyzz_1[i] * wa_z[i];

        g_xxz_0_xxxyzzz_0[i] = 3.0 * g_xx_0_xxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzz_1[i] * wa_z[i];

        g_xxz_0_xxxzzzz_0[i] = 4.0 * g_xx_0_xxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxzzzz_1[i] * wa_z[i];

        g_xxz_0_xxyyyyy_0[i] = g_xx_0_xxyyyyy_1[i] * wa_z[i];

        g_xxz_0_xxyyyyz_0[i] = g_xx_0_xxyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyz_1[i] * wa_z[i];

        g_xxz_0_xxyyyzz_0[i] = 2.0 * g_xx_0_xxyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyzz_1[i] * wa_z[i];

        g_xxz_0_xxyyzzz_0[i] = 3.0 * g_xx_0_xxyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzz_1[i] * wa_z[i];

        g_xxz_0_xxyzzzz_0[i] = 4.0 * g_xx_0_xxyzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzz_1[i] * wa_z[i];

        g_xxz_0_xxzzzzz_0[i] = 5.0 * g_xx_0_xxzzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyyyy_0[i] = g_xx_0_xyyyyyy_1[i] * wa_z[i];

        g_xxz_0_xyyyyyz_0[i] = g_xx_0_xyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyz_1[i] * wa_z[i];

        g_xxz_0_xyyyyzz_0[i] = 2.0 * g_xx_0_xyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyzz_1[i] * wa_z[i];

        g_xxz_0_xyyyzzz_0[i] = 3.0 * g_xx_0_xyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzz_1[i] * wa_z[i];

        g_xxz_0_xyyzzzz_0[i] = 4.0 * g_xx_0_xyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzz_1[i] * wa_z[i];

        g_xxz_0_xyzzzzz_0[i] = 5.0 * g_xx_0_xyzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzz_1[i] * wa_z[i];

        g_xxz_0_xzzzzzz_0[i] = 6.0 * g_xx_0_xzzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyyyy_0[i] = g_xx_0_yyyyyyy_1[i] * wa_z[i];

        g_xxz_0_yyyyyyz_0[i] = g_xx_0_yyyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyyz_1[i] * wa_z[i];

        g_xxz_0_yyyyyzz_0[i] = 2.0 * g_xx_0_yyyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyyzz_1[i] * wa_z[i];

        g_xxz_0_yyyyzzz_0[i] = 3.0 * g_xx_0_yyyyzz_1[i] * fi_acd_0 + g_xx_0_yyyyzzz_1[i] * wa_z[i];

        g_xxz_0_yyyzzzz_0[i] = 4.0 * g_xx_0_yyyzzz_1[i] * fi_acd_0 + g_xx_0_yyyzzzz_1[i] * wa_z[i];

        g_xxz_0_yyzzzzz_0[i] = 5.0 * g_xx_0_yyzzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzzz_1[i] * wa_z[i];

        g_xxz_0_yzzzzzz_0[i] = 6.0 * g_xx_0_yzzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzzz_1[i] * wa_z[i];

        g_xxz_0_zzzzzzz_0[i] = 7.0 * g_xx_0_zzzzzz_1[i] * fi_acd_0 + g_xx_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 108-144 components of targeted buffer : FSK

    auto g_xyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 108);

    auto g_xyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 109);

    auto g_xyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 110);

    auto g_xyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 111);

    auto g_xyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 112);

    auto g_xyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 113);

    auto g_xyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 114);

    auto g_xyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 115);

    auto g_xyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 116);

    auto g_xyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 117);

    auto g_xyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 118);

    auto g_xyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 119);

    auto g_xyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 120);

    auto g_xyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 121);

    auto g_xyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 122);

    auto g_xyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 123);

    auto g_xyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 124);

    auto g_xyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 125);

    auto g_xyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 126);

    auto g_xyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 127);

    auto g_xyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 128);

    auto g_xyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 129);

    auto g_xyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 130);

    auto g_xyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 131);

    auto g_xyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 132);

    auto g_xyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 133);

    auto g_xyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 134);

    auto g_xyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 135);

    auto g_xyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 136);

    auto g_xyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 137);

    auto g_xyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 138);

    auto g_xyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 139);

    auto g_xyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 140);

    auto g_xyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 141);

    auto g_xyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 142);

    auto g_xyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 143);

    #pragma omp simd aligned(g_xyy_0_xxxxxxx_0, g_xyy_0_xxxxxxy_0, g_xyy_0_xxxxxxz_0, g_xyy_0_xxxxxyy_0, g_xyy_0_xxxxxyz_0, g_xyy_0_xxxxxzz_0, g_xyy_0_xxxxyyy_0, g_xyy_0_xxxxyyz_0, g_xyy_0_xxxxyzz_0, g_xyy_0_xxxxzzz_0, g_xyy_0_xxxyyyy_0, g_xyy_0_xxxyyyz_0, g_xyy_0_xxxyyzz_0, g_xyy_0_xxxyzzz_0, g_xyy_0_xxxzzzz_0, g_xyy_0_xxyyyyy_0, g_xyy_0_xxyyyyz_0, g_xyy_0_xxyyyzz_0, g_xyy_0_xxyyzzz_0, g_xyy_0_xxyzzzz_0, g_xyy_0_xxzzzzz_0, g_xyy_0_xyyyyyy_0, g_xyy_0_xyyyyyz_0, g_xyy_0_xyyyyzz_0, g_xyy_0_xyyyzzz_0, g_xyy_0_xyyzzzz_0, g_xyy_0_xyzzzzz_0, g_xyy_0_xzzzzzz_0, g_xyy_0_yyyyyyy_0, g_xyy_0_yyyyyyz_0, g_xyy_0_yyyyyzz_0, g_xyy_0_yyyyzzz_0, g_xyy_0_yyyzzzz_0, g_xyy_0_yyzzzzz_0, g_xyy_0_yzzzzzz_0, g_xyy_0_zzzzzzz_0, g_yy_0_xxxxxx_1, g_yy_0_xxxxxxx_1, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxxz_1, g_yy_0_xxxxxy_1, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyz_1, g_yy_0_xxxxxz_1, g_yy_0_xxxxxzz_1, g_yy_0_xxxxyy_1, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyz_1, g_yy_0_xxxxyzz_1, g_yy_0_xxxxzz_1, g_yy_0_xxxxzzz_1, g_yy_0_xxxyyy_1, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyz_1, g_yy_0_xxxyyzz_1, g_yy_0_xxxyzz_1, g_yy_0_xxxyzzz_1, g_yy_0_xxxzzz_1, g_yy_0_xxxzzzz_1, g_yy_0_xxyyyy_1, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyz_1, g_yy_0_xxyyyzz_1, g_yy_0_xxyyzz_1, g_yy_0_xxyyzzz_1, g_yy_0_xxyzzz_1, g_yy_0_xxyzzzz_1, g_yy_0_xxzzzz_1, g_yy_0_xxzzzzz_1, g_yy_0_xyyyyy_1, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyz_1, g_yy_0_xyyyyzz_1, g_yy_0_xyyyzz_1, g_yy_0_xyyyzzz_1, g_yy_0_xyyzzz_1, g_yy_0_xyyzzzz_1, g_yy_0_xyzzzz_1, g_yy_0_xyzzzzz_1, g_yy_0_xzzzzz_1, g_yy_0_xzzzzzz_1, g_yy_0_yyyyyy_1, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyz_1, g_yy_0_yyyyyzz_1, g_yy_0_yyyyzz_1, g_yy_0_yyyyzzz_1, g_yy_0_yyyzzz_1, g_yy_0_yyyzzzz_1, g_yy_0_yyzzzz_1, g_yy_0_yyzzzzz_1, g_yy_0_yzzzzz_1, g_yy_0_yzzzzzz_1, g_yy_0_zzzzzz_1, g_yy_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xxxxxxx_0[i] = 7.0 * g_yy_0_xxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxx_1[i] * wa_x[i];

        g_xyy_0_xxxxxxy_0[i] = 6.0 * g_yy_0_xxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxxy_1[i] * wa_x[i];

        g_xyy_0_xxxxxxz_0[i] = 6.0 * g_yy_0_xxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxxz_1[i] * wa_x[i];

        g_xyy_0_xxxxxyy_0[i] = 5.0 * g_yy_0_xxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxxyy_1[i] * wa_x[i];

        g_xyy_0_xxxxxyz_0[i] = 5.0 * g_yy_0_xxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxxyz_1[i] * wa_x[i];

        g_xyy_0_xxxxxzz_0[i] = 5.0 * g_yy_0_xxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxxzz_1[i] * wa_x[i];

        g_xyy_0_xxxxyyy_0[i] = 4.0 * g_yy_0_xxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyy_1[i] * wa_x[i];

        g_xyy_0_xxxxyyz_0[i] = 4.0 * g_yy_0_xxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyz_1[i] * wa_x[i];

        g_xyy_0_xxxxyzz_0[i] = 4.0 * g_yy_0_xxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzz_1[i] * wa_x[i];

        g_xyy_0_xxxxzzz_0[i] = 4.0 * g_yy_0_xxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxxzzz_1[i] * wa_x[i];

        g_xyy_0_xxxyyyy_0[i] = 3.0 * g_yy_0_xxyyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyy_1[i] * wa_x[i];

        g_xyy_0_xxxyyyz_0[i] = 3.0 * g_yy_0_xxyyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyz_1[i] * wa_x[i];

        g_xyy_0_xxxyyzz_0[i] = 3.0 * g_yy_0_xxyyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzz_1[i] * wa_x[i];

        g_xyy_0_xxxyzzz_0[i] = 3.0 * g_yy_0_xxyzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzz_1[i] * wa_x[i];

        g_xyy_0_xxxzzzz_0[i] = 3.0 * g_yy_0_xxzzzz_1[i] * fi_acd_0 + g_yy_0_xxxzzzz_1[i] * wa_x[i];

        g_xyy_0_xxyyyyy_0[i] = 2.0 * g_yy_0_xyyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyy_1[i] * wa_x[i];

        g_xyy_0_xxyyyyz_0[i] = 2.0 * g_yy_0_xyyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyz_1[i] * wa_x[i];

        g_xyy_0_xxyyyzz_0[i] = 2.0 * g_yy_0_xyyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzz_1[i] * wa_x[i];

        g_xyy_0_xxyyzzz_0[i] = 2.0 * g_yy_0_xyyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzz_1[i] * wa_x[i];

        g_xyy_0_xxyzzzz_0[i] = 2.0 * g_yy_0_xyzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzz_1[i] * wa_x[i];

        g_xyy_0_xxzzzzz_0[i] = 2.0 * g_yy_0_xzzzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyyyy_0[i] = g_yy_0_yyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyy_1[i] * wa_x[i];

        g_xyy_0_xyyyyyz_0[i] = g_yy_0_yyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyz_1[i] * wa_x[i];

        g_xyy_0_xyyyyzz_0[i] = g_yy_0_yyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzz_1[i] * wa_x[i];

        g_xyy_0_xyyyzzz_0[i] = g_yy_0_yyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzz_1[i] * wa_x[i];

        g_xyy_0_xyyzzzz_0[i] = g_yy_0_yyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzz_1[i] * wa_x[i];

        g_xyy_0_xyzzzzz_0[i] = g_yy_0_yzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzz_1[i] * wa_x[i];

        g_xyy_0_xzzzzzz_0[i] = g_yy_0_zzzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyyyy_0[i] = g_yy_0_yyyyyyy_1[i] * wa_x[i];

        g_xyy_0_yyyyyyz_0[i] = g_yy_0_yyyyyyz_1[i] * wa_x[i];

        g_xyy_0_yyyyyzz_0[i] = g_yy_0_yyyyyzz_1[i] * wa_x[i];

        g_xyy_0_yyyyzzz_0[i] = g_yy_0_yyyyzzz_1[i] * wa_x[i];

        g_xyy_0_yyyzzzz_0[i] = g_yy_0_yyyzzzz_1[i] * wa_x[i];

        g_xyy_0_yyzzzzz_0[i] = g_yy_0_yyzzzzz_1[i] * wa_x[i];

        g_xyy_0_yzzzzzz_0[i] = g_yy_0_yzzzzzz_1[i] * wa_x[i];

        g_xyy_0_zzzzzzz_0[i] = g_yy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 144-180 components of targeted buffer : FSK

    auto g_xyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 144);

    auto g_xyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 145);

    auto g_xyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 146);

    auto g_xyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 147);

    auto g_xyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 148);

    auto g_xyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 149);

    auto g_xyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 150);

    auto g_xyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 151);

    auto g_xyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 152);

    auto g_xyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 153);

    auto g_xyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 154);

    auto g_xyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 155);

    auto g_xyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 156);

    auto g_xyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 157);

    auto g_xyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 158);

    auto g_xyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 159);

    auto g_xyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 160);

    auto g_xyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 161);

    auto g_xyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 162);

    auto g_xyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 163);

    auto g_xyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 164);

    auto g_xyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 165);

    auto g_xyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 166);

    auto g_xyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 167);

    auto g_xyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 168);

    auto g_xyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 169);

    auto g_xyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 170);

    auto g_xyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 171);

    auto g_xyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 172);

    auto g_xyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 173);

    auto g_xyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 174);

    auto g_xyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 175);

    auto g_xyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 176);

    auto g_xyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 177);

    auto g_xyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 178);

    auto g_xyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 179);

    #pragma omp simd aligned(g_xy_0_xxxxxxy_1, g_xy_0_xxxxxyy_1, g_xy_0_xxxxyyy_1, g_xy_0_xxxyyyy_1, g_xy_0_xxyyyyy_1, g_xy_0_xyyyyyy_1, g_xyz_0_xxxxxxx_0, g_xyz_0_xxxxxxy_0, g_xyz_0_xxxxxxz_0, g_xyz_0_xxxxxyy_0, g_xyz_0_xxxxxyz_0, g_xyz_0_xxxxxzz_0, g_xyz_0_xxxxyyy_0, g_xyz_0_xxxxyyz_0, g_xyz_0_xxxxyzz_0, g_xyz_0_xxxxzzz_0, g_xyz_0_xxxyyyy_0, g_xyz_0_xxxyyyz_0, g_xyz_0_xxxyyzz_0, g_xyz_0_xxxyzzz_0, g_xyz_0_xxxzzzz_0, g_xyz_0_xxyyyyy_0, g_xyz_0_xxyyyyz_0, g_xyz_0_xxyyyzz_0, g_xyz_0_xxyyzzz_0, g_xyz_0_xxyzzzz_0, g_xyz_0_xxzzzzz_0, g_xyz_0_xyyyyyy_0, g_xyz_0_xyyyyyz_0, g_xyz_0_xyyyyzz_0, g_xyz_0_xyyyzzz_0, g_xyz_0_xyyzzzz_0, g_xyz_0_xyzzzzz_0, g_xyz_0_xzzzzzz_0, g_xyz_0_yyyyyyy_0, g_xyz_0_yyyyyyz_0, g_xyz_0_yyyyyzz_0, g_xyz_0_yyyyzzz_0, g_xyz_0_yyyzzzz_0, g_xyz_0_yyzzzzz_0, g_xyz_0_yzzzzzz_0, g_xyz_0_zzzzzzz_0, g_xz_0_xxxxxxx_1, g_xz_0_xxxxxxz_1, g_xz_0_xxxxxzz_1, g_xz_0_xxxxzzz_1, g_xz_0_xxxzzzz_1, g_xz_0_xxzzzzz_1, g_xz_0_xzzzzzz_1, g_yz_0_xxxxxyz_1, g_yz_0_xxxxyyz_1, g_yz_0_xxxxyz_1, g_yz_0_xxxxyzz_1, g_yz_0_xxxyyyz_1, g_yz_0_xxxyyz_1, g_yz_0_xxxyyzz_1, g_yz_0_xxxyzz_1, g_yz_0_xxxyzzz_1, g_yz_0_xxyyyyz_1, g_yz_0_xxyyyz_1, g_yz_0_xxyyyzz_1, g_yz_0_xxyyzz_1, g_yz_0_xxyyzzz_1, g_yz_0_xxyzzz_1, g_yz_0_xxyzzzz_1, g_yz_0_xyyyyyz_1, g_yz_0_xyyyyz_1, g_yz_0_xyyyyzz_1, g_yz_0_xyyyzz_1, g_yz_0_xyyyzzz_1, g_yz_0_xyyzzz_1, g_yz_0_xyyzzzz_1, g_yz_0_xyzzzz_1, g_yz_0_xyzzzzz_1, g_yz_0_yyyyyyy_1, g_yz_0_yyyyyyz_1, g_yz_0_yyyyyz_1, g_yz_0_yyyyyzz_1, g_yz_0_yyyyzz_1, g_yz_0_yyyyzzz_1, g_yz_0_yyyzzz_1, g_yz_0_yyyzzzz_1, g_yz_0_yyzzzz_1, g_yz_0_yyzzzzz_1, g_yz_0_yzzzzz_1, g_yz_0_yzzzzzz_1, g_yz_0_zzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyz_0_xxxxxxx_0[i] = g_xz_0_xxxxxxx_1[i] * wa_y[i];

        g_xyz_0_xxxxxxy_0[i] = g_xy_0_xxxxxxy_1[i] * wa_z[i];

        g_xyz_0_xxxxxxz_0[i] = g_xz_0_xxxxxxz_1[i] * wa_y[i];

        g_xyz_0_xxxxxyy_0[i] = g_xy_0_xxxxxyy_1[i] * wa_z[i];

        g_xyz_0_xxxxxyz_0[i] = 5.0 * g_yz_0_xxxxyz_1[i] * fi_acd_0 + g_yz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyz_0_xxxxxzz_0[i] = g_xz_0_xxxxxzz_1[i] * wa_y[i];

        g_xyz_0_xxxxyyy_0[i] = g_xy_0_xxxxyyy_1[i] * wa_z[i];

        g_xyz_0_xxxxyyz_0[i] = 4.0 * g_yz_0_xxxyyz_1[i] * fi_acd_0 + g_yz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyz_0_xxxxyzz_0[i] = 4.0 * g_yz_0_xxxyzz_1[i] * fi_acd_0 + g_yz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyz_0_xxxxzzz_0[i] = g_xz_0_xxxxzzz_1[i] * wa_y[i];

        g_xyz_0_xxxyyyy_0[i] = g_xy_0_xxxyyyy_1[i] * wa_z[i];

        g_xyz_0_xxxyyyz_0[i] = 3.0 * g_yz_0_xxyyyz_1[i] * fi_acd_0 + g_yz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyz_0_xxxyyzz_0[i] = 3.0 * g_yz_0_xxyyzz_1[i] * fi_acd_0 + g_yz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyz_0_xxxyzzz_0[i] = 3.0 * g_yz_0_xxyzzz_1[i] * fi_acd_0 + g_yz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyz_0_xxxzzzz_0[i] = g_xz_0_xxxzzzz_1[i] * wa_y[i];

        g_xyz_0_xxyyyyy_0[i] = g_xy_0_xxyyyyy_1[i] * wa_z[i];

        g_xyz_0_xxyyyyz_0[i] = 2.0 * g_yz_0_xyyyyz_1[i] * fi_acd_0 + g_yz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyz_0_xxyyyzz_0[i] = 2.0 * g_yz_0_xyyyzz_1[i] * fi_acd_0 + g_yz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyz_0_xxyyzzz_0[i] = 2.0 * g_yz_0_xyyzzz_1[i] * fi_acd_0 + g_yz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyz_0_xxyzzzz_0[i] = 2.0 * g_yz_0_xyzzzz_1[i] * fi_acd_0 + g_yz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyz_0_xxzzzzz_0[i] = g_xz_0_xxzzzzz_1[i] * wa_y[i];

        g_xyz_0_xyyyyyy_0[i] = g_xy_0_xyyyyyy_1[i] * wa_z[i];

        g_xyz_0_xyyyyyz_0[i] = g_yz_0_yyyyyz_1[i] * fi_acd_0 + g_yz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyz_0_xyyyyzz_0[i] = g_yz_0_yyyyzz_1[i] * fi_acd_0 + g_yz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyz_0_xyyyzzz_0[i] = g_yz_0_yyyzzz_1[i] * fi_acd_0 + g_yz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyz_0_xyyzzzz_0[i] = g_yz_0_yyzzzz_1[i] * fi_acd_0 + g_yz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyz_0_xyzzzzz_0[i] = g_yz_0_yzzzzz_1[i] * fi_acd_0 + g_yz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyz_0_xzzzzzz_0[i] = g_xz_0_xzzzzzz_1[i] * wa_y[i];

        g_xyz_0_yyyyyyy_0[i] = g_yz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyz_0_yyyyyyz_0[i] = g_yz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyz_0_yyyyyzz_0[i] = g_yz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyz_0_yyyyzzz_0[i] = g_yz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyz_0_yyyzzzz_0[i] = g_yz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyz_0_yyzzzzz_0[i] = g_yz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyz_0_yzzzzzz_0[i] = g_yz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyz_0_zzzzzzz_0[i] = g_yz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 180-216 components of targeted buffer : FSK

    auto g_xzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 180);

    auto g_xzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 181);

    auto g_xzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 182);

    auto g_xzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 183);

    auto g_xzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 184);

    auto g_xzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 185);

    auto g_xzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 186);

    auto g_xzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 187);

    auto g_xzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 188);

    auto g_xzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 189);

    auto g_xzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 190);

    auto g_xzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 191);

    auto g_xzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 192);

    auto g_xzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 193);

    auto g_xzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 194);

    auto g_xzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 195);

    auto g_xzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 196);

    auto g_xzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 197);

    auto g_xzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 198);

    auto g_xzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 199);

    auto g_xzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 200);

    auto g_xzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 201);

    auto g_xzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 202);

    auto g_xzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 203);

    auto g_xzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 204);

    auto g_xzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 205);

    auto g_xzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 206);

    auto g_xzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 207);

    auto g_xzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 208);

    auto g_xzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 209);

    auto g_xzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 210);

    auto g_xzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 211);

    auto g_xzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 212);

    auto g_xzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 213);

    auto g_xzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 214);

    auto g_xzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 215);

    #pragma omp simd aligned(g_xzz_0_xxxxxxx_0, g_xzz_0_xxxxxxy_0, g_xzz_0_xxxxxxz_0, g_xzz_0_xxxxxyy_0, g_xzz_0_xxxxxyz_0, g_xzz_0_xxxxxzz_0, g_xzz_0_xxxxyyy_0, g_xzz_0_xxxxyyz_0, g_xzz_0_xxxxyzz_0, g_xzz_0_xxxxzzz_0, g_xzz_0_xxxyyyy_0, g_xzz_0_xxxyyyz_0, g_xzz_0_xxxyyzz_0, g_xzz_0_xxxyzzz_0, g_xzz_0_xxxzzzz_0, g_xzz_0_xxyyyyy_0, g_xzz_0_xxyyyyz_0, g_xzz_0_xxyyyzz_0, g_xzz_0_xxyyzzz_0, g_xzz_0_xxyzzzz_0, g_xzz_0_xxzzzzz_0, g_xzz_0_xyyyyyy_0, g_xzz_0_xyyyyyz_0, g_xzz_0_xyyyyzz_0, g_xzz_0_xyyyzzz_0, g_xzz_0_xyyzzzz_0, g_xzz_0_xyzzzzz_0, g_xzz_0_xzzzzzz_0, g_xzz_0_yyyyyyy_0, g_xzz_0_yyyyyyz_0, g_xzz_0_yyyyyzz_0, g_xzz_0_yyyyzzz_0, g_xzz_0_yyyzzzz_0, g_xzz_0_yyzzzzz_0, g_xzz_0_yzzzzzz_0, g_xzz_0_zzzzzzz_0, g_zz_0_xxxxxx_1, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxy_1, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxy_1, g_zz_0_xxxxxyy_1, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxz_1, g_zz_0_xxxxxzz_1, g_zz_0_xxxxyy_1, g_zz_0_xxxxyyy_1, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyz_1, g_zz_0_xxxxyzz_1, g_zz_0_xxxxzz_1, g_zz_0_xxxxzzz_1, g_zz_0_xxxyyy_1, g_zz_0_xxxyyyy_1, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyz_1, g_zz_0_xxxyyzz_1, g_zz_0_xxxyzz_1, g_zz_0_xxxyzzz_1, g_zz_0_xxxzzz_1, g_zz_0_xxxzzzz_1, g_zz_0_xxyyyy_1, g_zz_0_xxyyyyy_1, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyz_1, g_zz_0_xxyyyzz_1, g_zz_0_xxyyzz_1, g_zz_0_xxyyzzz_1, g_zz_0_xxyzzz_1, g_zz_0_xxyzzzz_1, g_zz_0_xxzzzz_1, g_zz_0_xxzzzzz_1, g_zz_0_xyyyyy_1, g_zz_0_xyyyyyy_1, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyz_1, g_zz_0_xyyyyzz_1, g_zz_0_xyyyzz_1, g_zz_0_xyyyzzz_1, g_zz_0_xyyzzz_1, g_zz_0_xyyzzzz_1, g_zz_0_xyzzzz_1, g_zz_0_xyzzzzz_1, g_zz_0_xzzzzz_1, g_zz_0_xzzzzzz_1, g_zz_0_yyyyyy_1, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyz_1, g_zz_0_yyyyyzz_1, g_zz_0_yyyyzz_1, g_zz_0_yyyyzzz_1, g_zz_0_yyyzzz_1, g_zz_0_yyyzzzz_1, g_zz_0_yyzzzz_1, g_zz_0_yyzzzzz_1, g_zz_0_yzzzzz_1, g_zz_0_yzzzzzz_1, g_zz_0_zzzzzz_1, g_zz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xxxxxxx_0[i] = 7.0 * g_zz_0_xxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxx_1[i] * wa_x[i];

        g_xzz_0_xxxxxxy_0[i] = 6.0 * g_zz_0_xxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxxy_1[i] * wa_x[i];

        g_xzz_0_xxxxxxz_0[i] = 6.0 * g_zz_0_xxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxxz_1[i] * wa_x[i];

        g_xzz_0_xxxxxyy_0[i] = 5.0 * g_zz_0_xxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxxyy_1[i] * wa_x[i];

        g_xzz_0_xxxxxyz_0[i] = 5.0 * g_zz_0_xxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxxyz_1[i] * wa_x[i];

        g_xzz_0_xxxxxzz_0[i] = 5.0 * g_zz_0_xxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxxzz_1[i] * wa_x[i];

        g_xzz_0_xxxxyyy_0[i] = 4.0 * g_zz_0_xxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyy_1[i] * wa_x[i];

        g_xzz_0_xxxxyyz_0[i] = 4.0 * g_zz_0_xxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyz_1[i] * wa_x[i];

        g_xzz_0_xxxxyzz_0[i] = 4.0 * g_zz_0_xxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzz_1[i] * wa_x[i];

        g_xzz_0_xxxxzzz_0[i] = 4.0 * g_zz_0_xxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxxzzz_1[i] * wa_x[i];

        g_xzz_0_xxxyyyy_0[i] = 3.0 * g_zz_0_xxyyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyy_1[i] * wa_x[i];

        g_xzz_0_xxxyyyz_0[i] = 3.0 * g_zz_0_xxyyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyz_1[i] * wa_x[i];

        g_xzz_0_xxxyyzz_0[i] = 3.0 * g_zz_0_xxyyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzz_1[i] * wa_x[i];

        g_xzz_0_xxxyzzz_0[i] = 3.0 * g_zz_0_xxyzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzz_1[i] * wa_x[i];

        g_xzz_0_xxxzzzz_0[i] = 3.0 * g_zz_0_xxzzzz_1[i] * fi_acd_0 + g_zz_0_xxxzzzz_1[i] * wa_x[i];

        g_xzz_0_xxyyyyy_0[i] = 2.0 * g_zz_0_xyyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyy_1[i] * wa_x[i];

        g_xzz_0_xxyyyyz_0[i] = 2.0 * g_zz_0_xyyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyz_1[i] * wa_x[i];

        g_xzz_0_xxyyyzz_0[i] = 2.0 * g_zz_0_xyyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzz_1[i] * wa_x[i];

        g_xzz_0_xxyyzzz_0[i] = 2.0 * g_zz_0_xyyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzz_1[i] * wa_x[i];

        g_xzz_0_xxyzzzz_0[i] = 2.0 * g_zz_0_xyzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzz_1[i] * wa_x[i];

        g_xzz_0_xxzzzzz_0[i] = 2.0 * g_zz_0_xzzzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyyyy_0[i] = g_zz_0_yyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyy_1[i] * wa_x[i];

        g_xzz_0_xyyyyyz_0[i] = g_zz_0_yyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyz_1[i] * wa_x[i];

        g_xzz_0_xyyyyzz_0[i] = g_zz_0_yyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzz_1[i] * wa_x[i];

        g_xzz_0_xyyyzzz_0[i] = g_zz_0_yyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzz_1[i] * wa_x[i];

        g_xzz_0_xyyzzzz_0[i] = g_zz_0_yyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzz_1[i] * wa_x[i];

        g_xzz_0_xyzzzzz_0[i] = g_zz_0_yzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzz_1[i] * wa_x[i];

        g_xzz_0_xzzzzzz_0[i] = g_zz_0_zzzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyyyy_0[i] = g_zz_0_yyyyyyy_1[i] * wa_x[i];

        g_xzz_0_yyyyyyz_0[i] = g_zz_0_yyyyyyz_1[i] * wa_x[i];

        g_xzz_0_yyyyyzz_0[i] = g_zz_0_yyyyyzz_1[i] * wa_x[i];

        g_xzz_0_yyyyzzz_0[i] = g_zz_0_yyyyzzz_1[i] * wa_x[i];

        g_xzz_0_yyyzzzz_0[i] = g_zz_0_yyyzzzz_1[i] * wa_x[i];

        g_xzz_0_yyzzzzz_0[i] = g_zz_0_yyzzzzz_1[i] * wa_x[i];

        g_xzz_0_yzzzzzz_0[i] = g_zz_0_yzzzzzz_1[i] * wa_x[i];

        g_xzz_0_zzzzzzz_0[i] = g_zz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 216-252 components of targeted buffer : FSK

    auto g_yyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 216);

    auto g_yyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 217);

    auto g_yyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 218);

    auto g_yyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 219);

    auto g_yyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 220);

    auto g_yyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 221);

    auto g_yyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 222);

    auto g_yyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 223);

    auto g_yyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 224);

    auto g_yyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 225);

    auto g_yyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 226);

    auto g_yyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 227);

    auto g_yyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 228);

    auto g_yyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 229);

    auto g_yyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 230);

    auto g_yyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 231);

    auto g_yyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 232);

    auto g_yyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 233);

    auto g_yyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 234);

    auto g_yyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 235);

    auto g_yyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 236);

    auto g_yyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 237);

    auto g_yyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 238);

    auto g_yyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 239);

    auto g_yyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 240);

    auto g_yyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 241);

    auto g_yyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 242);

    auto g_yyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 243);

    auto g_yyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 244);

    auto g_yyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 245);

    auto g_yyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 246);

    auto g_yyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 247);

    auto g_yyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 248);

    auto g_yyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 249);

    auto g_yyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 250);

    auto g_yyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 251);

    #pragma omp simd aligned(g_y_0_xxxxxxx_0, g_y_0_xxxxxxx_1, g_y_0_xxxxxxy_0, g_y_0_xxxxxxy_1, g_y_0_xxxxxxz_0, g_y_0_xxxxxxz_1, g_y_0_xxxxxyy_0, g_y_0_xxxxxyy_1, g_y_0_xxxxxyz_0, g_y_0_xxxxxyz_1, g_y_0_xxxxxzz_0, g_y_0_xxxxxzz_1, g_y_0_xxxxyyy_0, g_y_0_xxxxyyy_1, g_y_0_xxxxyyz_0, g_y_0_xxxxyyz_1, g_y_0_xxxxyzz_0, g_y_0_xxxxyzz_1, g_y_0_xxxxzzz_0, g_y_0_xxxxzzz_1, g_y_0_xxxyyyy_0, g_y_0_xxxyyyy_1, g_y_0_xxxyyyz_0, g_y_0_xxxyyyz_1, g_y_0_xxxyyzz_0, g_y_0_xxxyyzz_1, g_y_0_xxxyzzz_0, g_y_0_xxxyzzz_1, g_y_0_xxxzzzz_0, g_y_0_xxxzzzz_1, g_y_0_xxyyyyy_0, g_y_0_xxyyyyy_1, g_y_0_xxyyyyz_0, g_y_0_xxyyyyz_1, g_y_0_xxyyyzz_0, g_y_0_xxyyyzz_1, g_y_0_xxyyzzz_0, g_y_0_xxyyzzz_1, g_y_0_xxyzzzz_0, g_y_0_xxyzzzz_1, g_y_0_xxzzzzz_0, g_y_0_xxzzzzz_1, g_y_0_xyyyyyy_0, g_y_0_xyyyyyy_1, g_y_0_xyyyyyz_0, g_y_0_xyyyyyz_1, g_y_0_xyyyyzz_0, g_y_0_xyyyyzz_1, g_y_0_xyyyzzz_0, g_y_0_xyyyzzz_1, g_y_0_xyyzzzz_0, g_y_0_xyyzzzz_1, g_y_0_xyzzzzz_0, g_y_0_xyzzzzz_1, g_y_0_xzzzzzz_0, g_y_0_xzzzzzz_1, g_y_0_yyyyyyy_0, g_y_0_yyyyyyy_1, g_y_0_yyyyyyz_0, g_y_0_yyyyyyz_1, g_y_0_yyyyyzz_0, g_y_0_yyyyyzz_1, g_y_0_yyyyzzz_0, g_y_0_yyyyzzz_1, g_y_0_yyyzzzz_0, g_y_0_yyyzzzz_1, g_y_0_yyzzzzz_0, g_y_0_yyzzzzz_1, g_y_0_yzzzzzz_0, g_y_0_yzzzzzz_1, g_y_0_zzzzzzz_0, g_y_0_zzzzzzz_1, g_yy_0_xxxxxx_1, g_yy_0_xxxxxxx_1, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxxz_1, g_yy_0_xxxxxy_1, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyz_1, g_yy_0_xxxxxz_1, g_yy_0_xxxxxzz_1, g_yy_0_xxxxyy_1, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyz_1, g_yy_0_xxxxyzz_1, g_yy_0_xxxxzz_1, g_yy_0_xxxxzzz_1, g_yy_0_xxxyyy_1, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyz_1, g_yy_0_xxxyyzz_1, g_yy_0_xxxyzz_1, g_yy_0_xxxyzzz_1, g_yy_0_xxxzzz_1, g_yy_0_xxxzzzz_1, g_yy_0_xxyyyy_1, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyz_1, g_yy_0_xxyyyzz_1, g_yy_0_xxyyzz_1, g_yy_0_xxyyzzz_1, g_yy_0_xxyzzz_1, g_yy_0_xxyzzzz_1, g_yy_0_xxzzzz_1, g_yy_0_xxzzzzz_1, g_yy_0_xyyyyy_1, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyz_1, g_yy_0_xyyyyzz_1, g_yy_0_xyyyzz_1, g_yy_0_xyyyzzz_1, g_yy_0_xyyzzz_1, g_yy_0_xyyzzzz_1, g_yy_0_xyzzzz_1, g_yy_0_xyzzzzz_1, g_yy_0_xzzzzz_1, g_yy_0_xzzzzzz_1, g_yy_0_yyyyyy_1, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyz_1, g_yy_0_yyyyyzz_1, g_yy_0_yyyyzz_1, g_yy_0_yyyyzzz_1, g_yy_0_yyyzzz_1, g_yy_0_yyyzzzz_1, g_yy_0_yyzzzz_1, g_yy_0_yyzzzzz_1, g_yy_0_yzzzzz_1, g_yy_0_yzzzzzz_1, g_yy_0_zzzzzz_1, g_yy_0_zzzzzzz_1, g_yyy_0_xxxxxxx_0, g_yyy_0_xxxxxxy_0, g_yyy_0_xxxxxxz_0, g_yyy_0_xxxxxyy_0, g_yyy_0_xxxxxyz_0, g_yyy_0_xxxxxzz_0, g_yyy_0_xxxxyyy_0, g_yyy_0_xxxxyyz_0, g_yyy_0_xxxxyzz_0, g_yyy_0_xxxxzzz_0, g_yyy_0_xxxyyyy_0, g_yyy_0_xxxyyyz_0, g_yyy_0_xxxyyzz_0, g_yyy_0_xxxyzzz_0, g_yyy_0_xxxzzzz_0, g_yyy_0_xxyyyyy_0, g_yyy_0_xxyyyyz_0, g_yyy_0_xxyyyzz_0, g_yyy_0_xxyyzzz_0, g_yyy_0_xxyzzzz_0, g_yyy_0_xxzzzzz_0, g_yyy_0_xyyyyyy_0, g_yyy_0_xyyyyyz_0, g_yyy_0_xyyyyzz_0, g_yyy_0_xyyyzzz_0, g_yyy_0_xyyzzzz_0, g_yyy_0_xyzzzzz_0, g_yyy_0_xzzzzzz_0, g_yyy_0_yyyyyyy_0, g_yyy_0_yyyyyyz_0, g_yyy_0_yyyyyzz_0, g_yyy_0_yyyyzzz_0, g_yyy_0_yyyzzzz_0, g_yyy_0_yyzzzzz_0, g_yyy_0_yzzzzzz_0, g_yyy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xxxxxxx_0[i] = 2.0 * g_y_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxx_1[i] * fz_be_0 + g_yy_0_xxxxxxx_1[i] * wa_y[i];

        g_yyy_0_xxxxxxy_0[i] = 2.0 * g_y_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxy_1[i] * fz_be_0 + g_yy_0_xxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxy_1[i] * wa_y[i];

        g_yyy_0_xxxxxxz_0[i] = 2.0 * g_y_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxz_1[i] * fz_be_0 + g_yy_0_xxxxxxz_1[i] * wa_y[i];

        g_yyy_0_xxxxxyy_0[i] = 2.0 * g_y_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyy_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxyy_1[i] * wa_y[i];

        g_yyy_0_xxxxxyz_0[i] = 2.0 * g_y_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyz_1[i] * fz_be_0 + g_yy_0_xxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxyz_1[i] * wa_y[i];

        g_yyy_0_xxxxxzz_0[i] = 2.0 * g_y_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxzz_1[i] * fz_be_0 + g_yy_0_xxxxxzz_1[i] * wa_y[i];

        g_yyy_0_xxxxyyy_0[i] = 2.0 * g_y_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyy_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyy_1[i] * wa_y[i];

        g_yyy_0_xxxxyyz_0[i] = 2.0 * g_y_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyz_1[i] * wa_y[i];

        g_yyy_0_xxxxyzz_0[i] = 2.0 * g_y_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyzz_1[i] * fz_be_0 + g_yy_0_xxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzz_1[i] * wa_y[i];

        g_yyy_0_xxxxzzz_0[i] = 2.0 * g_y_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxzzz_1[i] * fz_be_0 + g_yy_0_xxxxzzz_1[i] * wa_y[i];

        g_yyy_0_xxxyyyy_0[i] = 2.0 * g_y_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyy_1[i] * fz_be_0 + 4.0 * g_yy_0_xxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyy_1[i] * wa_y[i];

        g_yyy_0_xxxyyyz_0[i] = 2.0 * g_y_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyz_1[i] * wa_y[i];

        g_yyy_0_xxxyyzz_0[i] = 2.0 * g_y_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzz_1[i] * wa_y[i];

        g_yyy_0_xxxyzzz_0[i] = 2.0 * g_y_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyzzz_1[i] * fz_be_0 + g_yy_0_xxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzz_1[i] * wa_y[i];

        g_yyy_0_xxxzzzz_0[i] = 2.0 * g_y_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxzzzz_1[i] * fz_be_0 + g_yy_0_xxxzzzz_1[i] * wa_y[i];

        g_yyy_0_xxyyyyy_0[i] = 2.0 * g_y_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyy_1[i] * fz_be_0 + 5.0 * g_yy_0_xxyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyy_1[i] * wa_y[i];

        g_yyy_0_xxyyyyz_0[i] = 2.0 * g_y_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yy_0_xxyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyz_1[i] * wa_y[i];

        g_yyy_0_xxyyyzz_0[i] = 2.0 * g_y_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzz_1[i] * wa_y[i];

        g_yyy_0_xxyyzzz_0[i] = 2.0 * g_y_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzz_1[i] * wa_y[i];

        g_yyy_0_xxyzzzz_0[i] = 2.0 * g_y_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyzzzz_1[i] * fz_be_0 + g_yy_0_xxzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzz_1[i] * wa_y[i];

        g_yyy_0_xxzzzzz_0[i] = 2.0 * g_y_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxzzzzz_1[i] * fz_be_0 + g_yy_0_xxzzzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyyyy_0[i] = 2.0 * g_y_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyy_1[i] * fz_be_0 + 6.0 * g_yy_0_xyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyy_1[i] * wa_y[i];

        g_yyy_0_xyyyyyz_0[i] = 2.0 * g_y_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yy_0_xyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyz_1[i] * wa_y[i];

        g_yyy_0_xyyyyzz_0[i] = 2.0 * g_y_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yy_0_xyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzz_1[i] * wa_y[i];

        g_yyy_0_xyyyzzz_0[i] = 2.0 * g_y_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzz_1[i] * wa_y[i];

        g_yyy_0_xyyzzzz_0[i] = 2.0 * g_y_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzz_1[i] * wa_y[i];

        g_yyy_0_xyzzzzz_0[i] = 2.0 * g_y_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyzzzzz_1[i] * fz_be_0 + g_yy_0_xzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzz_1[i] * wa_y[i];

        g_yyy_0_xzzzzzz_0[i] = 2.0 * g_y_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xzzzzzz_1[i] * fz_be_0 + g_yy_0_xzzzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyyyy_0[i] = 2.0 * g_y_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyy_1[i] * fz_be_0 + 7.0 * g_yy_0_yyyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyyy_1[i] * wa_y[i];

        g_yyy_0_yyyyyyz_0[i] = 2.0 * g_y_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yy_0_yyyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyyyz_1[i] * wa_y[i];

        g_yyy_0_yyyyyzz_0[i] = 2.0 * g_y_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yy_0_yyyyzz_1[i] * fi_acd_0 + g_yy_0_yyyyyzz_1[i] * wa_y[i];

        g_yyy_0_yyyyzzz_0[i] = 2.0 * g_y_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yy_0_yyyzzz_1[i] * fi_acd_0 + g_yy_0_yyyyzzz_1[i] * wa_y[i];

        g_yyy_0_yyyzzzz_0[i] = 2.0 * g_y_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_yyzzzz_1[i] * fi_acd_0 + g_yy_0_yyyzzzz_1[i] * wa_y[i];

        g_yyy_0_yyzzzzz_0[i] = 2.0 * g_y_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_yzzzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzzz_1[i] * wa_y[i];

        g_yyy_0_yzzzzzz_0[i] = 2.0 * g_y_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yzzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzzz_1[i] * wa_y[i];

        g_yyy_0_zzzzzzz_0[i] = 2.0 * g_y_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_zzzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 252-288 components of targeted buffer : FSK

    auto g_yyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 252);

    auto g_yyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 253);

    auto g_yyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 254);

    auto g_yyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 255);

    auto g_yyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 256);

    auto g_yyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 257);

    auto g_yyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 258);

    auto g_yyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 259);

    auto g_yyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 260);

    auto g_yyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 261);

    auto g_yyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 262);

    auto g_yyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 263);

    auto g_yyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 264);

    auto g_yyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 265);

    auto g_yyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 266);

    auto g_yyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 267);

    auto g_yyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 268);

    auto g_yyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 269);

    auto g_yyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 270);

    auto g_yyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 271);

    auto g_yyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 272);

    auto g_yyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 273);

    auto g_yyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 274);

    auto g_yyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 275);

    auto g_yyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 276);

    auto g_yyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 277);

    auto g_yyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 278);

    auto g_yyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 279);

    auto g_yyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 280);

    auto g_yyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 281);

    auto g_yyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 282);

    auto g_yyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 283);

    auto g_yyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 284);

    auto g_yyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 285);

    auto g_yyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 286);

    auto g_yyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 287);

    #pragma omp simd aligned(g_yy_0_xxxxxx_1, g_yy_0_xxxxxxx_1, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxxz_1, g_yy_0_xxxxxy_1, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyz_1, g_yy_0_xxxxxz_1, g_yy_0_xxxxxzz_1, g_yy_0_xxxxyy_1, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyz_1, g_yy_0_xxxxyzz_1, g_yy_0_xxxxzz_1, g_yy_0_xxxxzzz_1, g_yy_0_xxxyyy_1, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyz_1, g_yy_0_xxxyyzz_1, g_yy_0_xxxyzz_1, g_yy_0_xxxyzzz_1, g_yy_0_xxxzzz_1, g_yy_0_xxxzzzz_1, g_yy_0_xxyyyy_1, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyz_1, g_yy_0_xxyyyzz_1, g_yy_0_xxyyzz_1, g_yy_0_xxyyzzz_1, g_yy_0_xxyzzz_1, g_yy_0_xxyzzzz_1, g_yy_0_xxzzzz_1, g_yy_0_xxzzzzz_1, g_yy_0_xyyyyy_1, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyz_1, g_yy_0_xyyyyzz_1, g_yy_0_xyyyzz_1, g_yy_0_xyyyzzz_1, g_yy_0_xyyzzz_1, g_yy_0_xyyzzzz_1, g_yy_0_xyzzzz_1, g_yy_0_xyzzzzz_1, g_yy_0_xzzzzz_1, g_yy_0_xzzzzzz_1, g_yy_0_yyyyyy_1, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyz_1, g_yy_0_yyyyyzz_1, g_yy_0_yyyyzz_1, g_yy_0_yyyyzzz_1, g_yy_0_yyyzzz_1, g_yy_0_yyyzzzz_1, g_yy_0_yyzzzz_1, g_yy_0_yyzzzzz_1, g_yy_0_yzzzzz_1, g_yy_0_yzzzzzz_1, g_yy_0_zzzzzz_1, g_yy_0_zzzzzzz_1, g_yyz_0_xxxxxxx_0, g_yyz_0_xxxxxxy_0, g_yyz_0_xxxxxxz_0, g_yyz_0_xxxxxyy_0, g_yyz_0_xxxxxyz_0, g_yyz_0_xxxxxzz_0, g_yyz_0_xxxxyyy_0, g_yyz_0_xxxxyyz_0, g_yyz_0_xxxxyzz_0, g_yyz_0_xxxxzzz_0, g_yyz_0_xxxyyyy_0, g_yyz_0_xxxyyyz_0, g_yyz_0_xxxyyzz_0, g_yyz_0_xxxyzzz_0, g_yyz_0_xxxzzzz_0, g_yyz_0_xxyyyyy_0, g_yyz_0_xxyyyyz_0, g_yyz_0_xxyyyzz_0, g_yyz_0_xxyyzzz_0, g_yyz_0_xxyzzzz_0, g_yyz_0_xxzzzzz_0, g_yyz_0_xyyyyyy_0, g_yyz_0_xyyyyyz_0, g_yyz_0_xyyyyzz_0, g_yyz_0_xyyyzzz_0, g_yyz_0_xyyzzzz_0, g_yyz_0_xyzzzzz_0, g_yyz_0_xzzzzzz_0, g_yyz_0_yyyyyyy_0, g_yyz_0_yyyyyyz_0, g_yyz_0_yyyyyzz_0, g_yyz_0_yyyyzzz_0, g_yyz_0_yyyzzzz_0, g_yyz_0_yyzzzzz_0, g_yyz_0_yzzzzzz_0, g_yyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xxxxxxx_0[i] = g_yy_0_xxxxxxx_1[i] * wa_z[i];

        g_yyz_0_xxxxxxy_0[i] = g_yy_0_xxxxxxy_1[i] * wa_z[i];

        g_yyz_0_xxxxxxz_0[i] = g_yy_0_xxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxz_1[i] * wa_z[i];

        g_yyz_0_xxxxxyy_0[i] = g_yy_0_xxxxxyy_1[i] * wa_z[i];

        g_yyz_0_xxxxxyz_0[i] = g_yy_0_xxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxyz_1[i] * wa_z[i];

        g_yyz_0_xxxxxzz_0[i] = 2.0 * g_yy_0_xxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxzz_1[i] * wa_z[i];

        g_yyz_0_xxxxyyy_0[i] = g_yy_0_xxxxyyy_1[i] * wa_z[i];

        g_yyz_0_xxxxyyz_0[i] = g_yy_0_xxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyz_1[i] * wa_z[i];

        g_yyz_0_xxxxyzz_0[i] = 2.0 * g_yy_0_xxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxyzz_1[i] * wa_z[i];

        g_yyz_0_xxxxzzz_0[i] = 3.0 * g_yy_0_xxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxzzz_1[i] * wa_z[i];

        g_yyz_0_xxxyyyy_0[i] = g_yy_0_xxxyyyy_1[i] * wa_z[i];

        g_yyz_0_xxxyyyz_0[i] = g_yy_0_xxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyz_1[i] * wa_z[i];

        g_yyz_0_xxxyyzz_0[i] = 2.0 * g_yy_0_xxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyzz_1[i] * wa_z[i];

        g_yyz_0_xxxyzzz_0[i] = 3.0 * g_yy_0_xxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzz_1[i] * wa_z[i];

        g_yyz_0_xxxzzzz_0[i] = 4.0 * g_yy_0_xxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxzzzz_1[i] * wa_z[i];

        g_yyz_0_xxyyyyy_0[i] = g_yy_0_xxyyyyy_1[i] * wa_z[i];

        g_yyz_0_xxyyyyz_0[i] = g_yy_0_xxyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyz_1[i] * wa_z[i];

        g_yyz_0_xxyyyzz_0[i] = 2.0 * g_yy_0_xxyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyzz_1[i] * wa_z[i];

        g_yyz_0_xxyyzzz_0[i] = 3.0 * g_yy_0_xxyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzz_1[i] * wa_z[i];

        g_yyz_0_xxyzzzz_0[i] = 4.0 * g_yy_0_xxyzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzz_1[i] * wa_z[i];

        g_yyz_0_xxzzzzz_0[i] = 5.0 * g_yy_0_xxzzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyyyy_0[i] = g_yy_0_xyyyyyy_1[i] * wa_z[i];

        g_yyz_0_xyyyyyz_0[i] = g_yy_0_xyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyz_1[i] * wa_z[i];

        g_yyz_0_xyyyyzz_0[i] = 2.0 * g_yy_0_xyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyzz_1[i] * wa_z[i];

        g_yyz_0_xyyyzzz_0[i] = 3.0 * g_yy_0_xyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzz_1[i] * wa_z[i];

        g_yyz_0_xyyzzzz_0[i] = 4.0 * g_yy_0_xyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzz_1[i] * wa_z[i];

        g_yyz_0_xyzzzzz_0[i] = 5.0 * g_yy_0_xyzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzz_1[i] * wa_z[i];

        g_yyz_0_xzzzzzz_0[i] = 6.0 * g_yy_0_xzzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyyyy_0[i] = g_yy_0_yyyyyyy_1[i] * wa_z[i];

        g_yyz_0_yyyyyyz_0[i] = g_yy_0_yyyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyyz_1[i] * wa_z[i];

        g_yyz_0_yyyyyzz_0[i] = 2.0 * g_yy_0_yyyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyyzz_1[i] * wa_z[i];

        g_yyz_0_yyyyzzz_0[i] = 3.0 * g_yy_0_yyyyzz_1[i] * fi_acd_0 + g_yy_0_yyyyzzz_1[i] * wa_z[i];

        g_yyz_0_yyyzzzz_0[i] = 4.0 * g_yy_0_yyyzzz_1[i] * fi_acd_0 + g_yy_0_yyyzzzz_1[i] * wa_z[i];

        g_yyz_0_yyzzzzz_0[i] = 5.0 * g_yy_0_yyzzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzzz_1[i] * wa_z[i];

        g_yyz_0_yzzzzzz_0[i] = 6.0 * g_yy_0_yzzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzzz_1[i] * wa_z[i];

        g_yyz_0_zzzzzzz_0[i] = 7.0 * g_yy_0_zzzzzz_1[i] * fi_acd_0 + g_yy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 288-324 components of targeted buffer : FSK

    auto g_yzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 288);

    auto g_yzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 289);

    auto g_yzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 290);

    auto g_yzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 291);

    auto g_yzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 292);

    auto g_yzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 293);

    auto g_yzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 294);

    auto g_yzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 295);

    auto g_yzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 296);

    auto g_yzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 297);

    auto g_yzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 298);

    auto g_yzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 299);

    auto g_yzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 300);

    auto g_yzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 301);

    auto g_yzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 302);

    auto g_yzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 303);

    auto g_yzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 304);

    auto g_yzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 305);

    auto g_yzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 306);

    auto g_yzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 307);

    auto g_yzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 308);

    auto g_yzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 309);

    auto g_yzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 310);

    auto g_yzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 311);

    auto g_yzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 312);

    auto g_yzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 313);

    auto g_yzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 314);

    auto g_yzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 315);

    auto g_yzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 316);

    auto g_yzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 317);

    auto g_yzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 318);

    auto g_yzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 319);

    auto g_yzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 320);

    auto g_yzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 321);

    auto g_yzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 322);

    auto g_yzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 323);

    #pragma omp simd aligned(g_yzz_0_xxxxxxx_0, g_yzz_0_xxxxxxy_0, g_yzz_0_xxxxxxz_0, g_yzz_0_xxxxxyy_0, g_yzz_0_xxxxxyz_0, g_yzz_0_xxxxxzz_0, g_yzz_0_xxxxyyy_0, g_yzz_0_xxxxyyz_0, g_yzz_0_xxxxyzz_0, g_yzz_0_xxxxzzz_0, g_yzz_0_xxxyyyy_0, g_yzz_0_xxxyyyz_0, g_yzz_0_xxxyyzz_0, g_yzz_0_xxxyzzz_0, g_yzz_0_xxxzzzz_0, g_yzz_0_xxyyyyy_0, g_yzz_0_xxyyyyz_0, g_yzz_0_xxyyyzz_0, g_yzz_0_xxyyzzz_0, g_yzz_0_xxyzzzz_0, g_yzz_0_xxzzzzz_0, g_yzz_0_xyyyyyy_0, g_yzz_0_xyyyyyz_0, g_yzz_0_xyyyyzz_0, g_yzz_0_xyyyzzz_0, g_yzz_0_xyyzzzz_0, g_yzz_0_xyzzzzz_0, g_yzz_0_xzzzzzz_0, g_yzz_0_yyyyyyy_0, g_yzz_0_yyyyyyz_0, g_yzz_0_yyyyyzz_0, g_yzz_0_yyyyzzz_0, g_yzz_0_yyyzzzz_0, g_yzz_0_yyzzzzz_0, g_yzz_0_yzzzzzz_0, g_yzz_0_zzzzzzz_0, g_zz_0_xxxxxx_1, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxy_1, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxy_1, g_zz_0_xxxxxyy_1, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxz_1, g_zz_0_xxxxxzz_1, g_zz_0_xxxxyy_1, g_zz_0_xxxxyyy_1, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyz_1, g_zz_0_xxxxyzz_1, g_zz_0_xxxxzz_1, g_zz_0_xxxxzzz_1, g_zz_0_xxxyyy_1, g_zz_0_xxxyyyy_1, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyz_1, g_zz_0_xxxyyzz_1, g_zz_0_xxxyzz_1, g_zz_0_xxxyzzz_1, g_zz_0_xxxzzz_1, g_zz_0_xxxzzzz_1, g_zz_0_xxyyyy_1, g_zz_0_xxyyyyy_1, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyz_1, g_zz_0_xxyyyzz_1, g_zz_0_xxyyzz_1, g_zz_0_xxyyzzz_1, g_zz_0_xxyzzz_1, g_zz_0_xxyzzzz_1, g_zz_0_xxzzzz_1, g_zz_0_xxzzzzz_1, g_zz_0_xyyyyy_1, g_zz_0_xyyyyyy_1, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyz_1, g_zz_0_xyyyyzz_1, g_zz_0_xyyyzz_1, g_zz_0_xyyyzzz_1, g_zz_0_xyyzzz_1, g_zz_0_xyyzzzz_1, g_zz_0_xyzzzz_1, g_zz_0_xyzzzzz_1, g_zz_0_xzzzzz_1, g_zz_0_xzzzzzz_1, g_zz_0_yyyyyy_1, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyz_1, g_zz_0_yyyyyzz_1, g_zz_0_yyyyzz_1, g_zz_0_yyyyzzz_1, g_zz_0_yyyzzz_1, g_zz_0_yyyzzzz_1, g_zz_0_yyzzzz_1, g_zz_0_yyzzzzz_1, g_zz_0_yzzzzz_1, g_zz_0_yzzzzzz_1, g_zz_0_zzzzzz_1, g_zz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xxxxxxx_0[i] = g_zz_0_xxxxxxx_1[i] * wa_y[i];

        g_yzz_0_xxxxxxy_0[i] = g_zz_0_xxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxy_1[i] * wa_y[i];

        g_yzz_0_xxxxxxz_0[i] = g_zz_0_xxxxxxz_1[i] * wa_y[i];

        g_yzz_0_xxxxxyy_0[i] = 2.0 * g_zz_0_xxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxyy_1[i] * wa_y[i];

        g_yzz_0_xxxxxyz_0[i] = g_zz_0_xxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxyz_1[i] * wa_y[i];

        g_yzz_0_xxxxxzz_0[i] = g_zz_0_xxxxxzz_1[i] * wa_y[i];

        g_yzz_0_xxxxyyy_0[i] = 3.0 * g_zz_0_xxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyy_1[i] * wa_y[i];

        g_yzz_0_xxxxyyz_0[i] = 2.0 * g_zz_0_xxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyz_1[i] * wa_y[i];

        g_yzz_0_xxxxyzz_0[i] = g_zz_0_xxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzz_1[i] * wa_y[i];

        g_yzz_0_xxxxzzz_0[i] = g_zz_0_xxxxzzz_1[i] * wa_y[i];

        g_yzz_0_xxxyyyy_0[i] = 4.0 * g_zz_0_xxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyy_1[i] * wa_y[i];

        g_yzz_0_xxxyyyz_0[i] = 3.0 * g_zz_0_xxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyz_1[i] * wa_y[i];

        g_yzz_0_xxxyyzz_0[i] = 2.0 * g_zz_0_xxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzz_1[i] * wa_y[i];

        g_yzz_0_xxxyzzz_0[i] = g_zz_0_xxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzz_1[i] * wa_y[i];

        g_yzz_0_xxxzzzz_0[i] = g_zz_0_xxxzzzz_1[i] * wa_y[i];

        g_yzz_0_xxyyyyy_0[i] = 5.0 * g_zz_0_xxyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyy_1[i] * wa_y[i];

        g_yzz_0_xxyyyyz_0[i] = 4.0 * g_zz_0_xxyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyz_1[i] * wa_y[i];

        g_yzz_0_xxyyyzz_0[i] = 3.0 * g_zz_0_xxyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzz_1[i] * wa_y[i];

        g_yzz_0_xxyyzzz_0[i] = 2.0 * g_zz_0_xxyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzz_1[i] * wa_y[i];

        g_yzz_0_xxyzzzz_0[i] = g_zz_0_xxzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzz_1[i] * wa_y[i];

        g_yzz_0_xxzzzzz_0[i] = g_zz_0_xxzzzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyyyy_0[i] = 6.0 * g_zz_0_xyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyy_1[i] * wa_y[i];

        g_yzz_0_xyyyyyz_0[i] = 5.0 * g_zz_0_xyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyz_1[i] * wa_y[i];

        g_yzz_0_xyyyyzz_0[i] = 4.0 * g_zz_0_xyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzz_1[i] * wa_y[i];

        g_yzz_0_xyyyzzz_0[i] = 3.0 * g_zz_0_xyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzz_1[i] * wa_y[i];

        g_yzz_0_xyyzzzz_0[i] = 2.0 * g_zz_0_xyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzz_1[i] * wa_y[i];

        g_yzz_0_xyzzzzz_0[i] = g_zz_0_xzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzz_1[i] * wa_y[i];

        g_yzz_0_xzzzzzz_0[i] = g_zz_0_xzzzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyyyy_0[i] = 7.0 * g_zz_0_yyyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyyy_1[i] * wa_y[i];

        g_yzz_0_yyyyyyz_0[i] = 6.0 * g_zz_0_yyyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyyyz_1[i] * wa_y[i];

        g_yzz_0_yyyyyzz_0[i] = 5.0 * g_zz_0_yyyyzz_1[i] * fi_acd_0 + g_zz_0_yyyyyzz_1[i] * wa_y[i];

        g_yzz_0_yyyyzzz_0[i] = 4.0 * g_zz_0_yyyzzz_1[i] * fi_acd_0 + g_zz_0_yyyyzzz_1[i] * wa_y[i];

        g_yzz_0_yyyzzzz_0[i] = 3.0 * g_zz_0_yyzzzz_1[i] * fi_acd_0 + g_zz_0_yyyzzzz_1[i] * wa_y[i];

        g_yzz_0_yyzzzzz_0[i] = 2.0 * g_zz_0_yzzzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzzz_1[i] * wa_y[i];

        g_yzz_0_yzzzzzz_0[i] = g_zz_0_zzzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzzz_1[i] * wa_y[i];

        g_yzz_0_zzzzzzz_0[i] = g_zz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 324-360 components of targeted buffer : FSK

    auto g_zzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 324);

    auto g_zzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 325);

    auto g_zzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 326);

    auto g_zzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 327);

    auto g_zzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 328);

    auto g_zzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 329);

    auto g_zzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 330);

    auto g_zzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 331);

    auto g_zzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 332);

    auto g_zzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 333);

    auto g_zzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 334);

    auto g_zzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 335);

    auto g_zzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 336);

    auto g_zzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 337);

    auto g_zzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 338);

    auto g_zzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 339);

    auto g_zzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 340);

    auto g_zzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 341);

    auto g_zzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 342);

    auto g_zzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 343);

    auto g_zzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 344);

    auto g_zzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 345);

    auto g_zzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 346);

    auto g_zzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 347);

    auto g_zzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 348);

    auto g_zzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 349);

    auto g_zzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 350);

    auto g_zzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 351);

    auto g_zzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 352);

    auto g_zzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 353);

    auto g_zzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 354);

    auto g_zzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 355);

    auto g_zzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 356);

    auto g_zzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 357);

    auto g_zzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 358);

    auto g_zzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 359);

    #pragma omp simd aligned(g_z_0_xxxxxxx_0, g_z_0_xxxxxxx_1, g_z_0_xxxxxxy_0, g_z_0_xxxxxxy_1, g_z_0_xxxxxxz_0, g_z_0_xxxxxxz_1, g_z_0_xxxxxyy_0, g_z_0_xxxxxyy_1, g_z_0_xxxxxyz_0, g_z_0_xxxxxyz_1, g_z_0_xxxxxzz_0, g_z_0_xxxxxzz_1, g_z_0_xxxxyyy_0, g_z_0_xxxxyyy_1, g_z_0_xxxxyyz_0, g_z_0_xxxxyyz_1, g_z_0_xxxxyzz_0, g_z_0_xxxxyzz_1, g_z_0_xxxxzzz_0, g_z_0_xxxxzzz_1, g_z_0_xxxyyyy_0, g_z_0_xxxyyyy_1, g_z_0_xxxyyyz_0, g_z_0_xxxyyyz_1, g_z_0_xxxyyzz_0, g_z_0_xxxyyzz_1, g_z_0_xxxyzzz_0, g_z_0_xxxyzzz_1, g_z_0_xxxzzzz_0, g_z_0_xxxzzzz_1, g_z_0_xxyyyyy_0, g_z_0_xxyyyyy_1, g_z_0_xxyyyyz_0, g_z_0_xxyyyyz_1, g_z_0_xxyyyzz_0, g_z_0_xxyyyzz_1, g_z_0_xxyyzzz_0, g_z_0_xxyyzzz_1, g_z_0_xxyzzzz_0, g_z_0_xxyzzzz_1, g_z_0_xxzzzzz_0, g_z_0_xxzzzzz_1, g_z_0_xyyyyyy_0, g_z_0_xyyyyyy_1, g_z_0_xyyyyyz_0, g_z_0_xyyyyyz_1, g_z_0_xyyyyzz_0, g_z_0_xyyyyzz_1, g_z_0_xyyyzzz_0, g_z_0_xyyyzzz_1, g_z_0_xyyzzzz_0, g_z_0_xyyzzzz_1, g_z_0_xyzzzzz_0, g_z_0_xyzzzzz_1, g_z_0_xzzzzzz_0, g_z_0_xzzzzzz_1, g_z_0_yyyyyyy_0, g_z_0_yyyyyyy_1, g_z_0_yyyyyyz_0, g_z_0_yyyyyyz_1, g_z_0_yyyyyzz_0, g_z_0_yyyyyzz_1, g_z_0_yyyyzzz_0, g_z_0_yyyyzzz_1, g_z_0_yyyzzzz_0, g_z_0_yyyzzzz_1, g_z_0_yyzzzzz_0, g_z_0_yyzzzzz_1, g_z_0_yzzzzzz_0, g_z_0_yzzzzzz_1, g_z_0_zzzzzzz_0, g_z_0_zzzzzzz_1, g_zz_0_xxxxxx_1, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxy_1, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxy_1, g_zz_0_xxxxxyy_1, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxz_1, g_zz_0_xxxxxzz_1, g_zz_0_xxxxyy_1, g_zz_0_xxxxyyy_1, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyz_1, g_zz_0_xxxxyzz_1, g_zz_0_xxxxzz_1, g_zz_0_xxxxzzz_1, g_zz_0_xxxyyy_1, g_zz_0_xxxyyyy_1, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyz_1, g_zz_0_xxxyyzz_1, g_zz_0_xxxyzz_1, g_zz_0_xxxyzzz_1, g_zz_0_xxxzzz_1, g_zz_0_xxxzzzz_1, g_zz_0_xxyyyy_1, g_zz_0_xxyyyyy_1, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyz_1, g_zz_0_xxyyyzz_1, g_zz_0_xxyyzz_1, g_zz_0_xxyyzzz_1, g_zz_0_xxyzzz_1, g_zz_0_xxyzzzz_1, g_zz_0_xxzzzz_1, g_zz_0_xxzzzzz_1, g_zz_0_xyyyyy_1, g_zz_0_xyyyyyy_1, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyz_1, g_zz_0_xyyyyzz_1, g_zz_0_xyyyzz_1, g_zz_0_xyyyzzz_1, g_zz_0_xyyzzz_1, g_zz_0_xyyzzzz_1, g_zz_0_xyzzzz_1, g_zz_0_xyzzzzz_1, g_zz_0_xzzzzz_1, g_zz_0_xzzzzzz_1, g_zz_0_yyyyyy_1, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyz_1, g_zz_0_yyyyyzz_1, g_zz_0_yyyyzz_1, g_zz_0_yyyyzzz_1, g_zz_0_yyyzzz_1, g_zz_0_yyyzzzz_1, g_zz_0_yyzzzz_1, g_zz_0_yyzzzzz_1, g_zz_0_yzzzzz_1, g_zz_0_yzzzzzz_1, g_zz_0_zzzzzz_1, g_zz_0_zzzzzzz_1, g_zzz_0_xxxxxxx_0, g_zzz_0_xxxxxxy_0, g_zzz_0_xxxxxxz_0, g_zzz_0_xxxxxyy_0, g_zzz_0_xxxxxyz_0, g_zzz_0_xxxxxzz_0, g_zzz_0_xxxxyyy_0, g_zzz_0_xxxxyyz_0, g_zzz_0_xxxxyzz_0, g_zzz_0_xxxxzzz_0, g_zzz_0_xxxyyyy_0, g_zzz_0_xxxyyyz_0, g_zzz_0_xxxyyzz_0, g_zzz_0_xxxyzzz_0, g_zzz_0_xxxzzzz_0, g_zzz_0_xxyyyyy_0, g_zzz_0_xxyyyyz_0, g_zzz_0_xxyyyzz_0, g_zzz_0_xxyyzzz_0, g_zzz_0_xxyzzzz_0, g_zzz_0_xxzzzzz_0, g_zzz_0_xyyyyyy_0, g_zzz_0_xyyyyyz_0, g_zzz_0_xyyyyzz_0, g_zzz_0_xyyyzzz_0, g_zzz_0_xyyzzzz_0, g_zzz_0_xyzzzzz_0, g_zzz_0_xzzzzzz_0, g_zzz_0_yyyyyyy_0, g_zzz_0_yyyyyyz_0, g_zzz_0_yyyyyzz_0, g_zzz_0_yyyyzzz_0, g_zzz_0_yyyzzzz_0, g_zzz_0_yyzzzzz_0, g_zzz_0_yzzzzzz_0, g_zzz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xxxxxxx_0[i] = 2.0 * g_z_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxx_1[i] * fz_be_0 + g_zz_0_xxxxxxx_1[i] * wa_z[i];

        g_zzz_0_xxxxxxy_0[i] = 2.0 * g_z_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxy_1[i] * fz_be_0 + g_zz_0_xxxxxxy_1[i] * wa_z[i];

        g_zzz_0_xxxxxxz_0[i] = 2.0 * g_z_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxz_1[i] * fz_be_0 + g_zz_0_xxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxz_1[i] * wa_z[i];

        g_zzz_0_xxxxxyy_0[i] = 2.0 * g_z_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyy_1[i] * fz_be_0 + g_zz_0_xxxxxyy_1[i] * wa_z[i];

        g_zzz_0_xxxxxyz_0[i] = 2.0 * g_z_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyz_1[i] * fz_be_0 + g_zz_0_xxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxyz_1[i] * wa_z[i];

        g_zzz_0_xxxxxzz_0[i] = 2.0 * g_z_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxzz_1[i] * wa_z[i];

        g_zzz_0_xxxxyyy_0[i] = 2.0 * g_z_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyy_1[i] * fz_be_0 + g_zz_0_xxxxyyy_1[i] * wa_z[i];

        g_zzz_0_xxxxyyz_0[i] = 2.0 * g_z_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyz_1[i] * fz_be_0 + g_zz_0_xxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyz_1[i] * wa_z[i];

        g_zzz_0_xxxxyzz_0[i] = 2.0 * g_z_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxyzz_1[i] * wa_z[i];

        g_zzz_0_xxxxzzz_0[i] = 2.0 * g_z_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxzzz_1[i] * wa_z[i];

        g_zzz_0_xxxyyyy_0[i] = 2.0 * g_z_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyy_1[i] * fz_be_0 + g_zz_0_xxxyyyy_1[i] * wa_z[i];

        g_zzz_0_xxxyyyz_0[i] = 2.0 * g_z_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyz_1[i] * fz_be_0 + g_zz_0_xxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyz_1[i] * wa_z[i];

        g_zzz_0_xxxyyzz_0[i] = 2.0 * g_z_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyzz_1[i] * wa_z[i];

        g_zzz_0_xxxyzzz_0[i] = 2.0 * g_z_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzz_1[i] * wa_z[i];

        g_zzz_0_xxxzzzz_0[i] = 2.0 * g_z_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxzzzz_1[i] * wa_z[i];

        g_zzz_0_xxyyyyy_0[i] = 2.0 * g_z_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyy_1[i] * fz_be_0 + g_zz_0_xxyyyyy_1[i] * wa_z[i];

        g_zzz_0_xxyyyyz_0[i] = 2.0 * g_z_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyz_1[i] * fz_be_0 + g_zz_0_xxyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyz_1[i] * wa_z[i];

        g_zzz_0_xxyyyzz_0[i] = 2.0 * g_z_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyzz_1[i] * wa_z[i];

        g_zzz_0_xxyyzzz_0[i] = 2.0 * g_z_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzz_1[i] * wa_z[i];

        g_zzz_0_xxyzzzz_0[i] = 2.0 * g_z_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxyzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzz_1[i] * wa_z[i];

        g_zzz_0_xxzzzzz_0[i] = 2.0 * g_z_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xxzzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyyyy_0[i] = 2.0 * g_z_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyy_1[i] * fz_be_0 + g_zz_0_xyyyyyy_1[i] * wa_z[i];

        g_zzz_0_xyyyyyz_0[i] = 2.0 * g_z_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyz_1[i] * fz_be_0 + g_zz_0_xyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyz_1[i] * wa_z[i];

        g_zzz_0_xyyyyzz_0[i] = 2.0 * g_z_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyzz_1[i] * wa_z[i];

        g_zzz_0_xyyyzzz_0[i] = 2.0 * g_z_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzz_1[i] * wa_z[i];

        g_zzz_0_xyyzzzz_0[i] = 2.0 * g_z_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzz_1[i] * wa_z[i];

        g_zzz_0_xyzzzzz_0[i] = 2.0 * g_z_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xyzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzz_1[i] * wa_z[i];

        g_zzz_0_xzzzzzz_0[i] = 2.0 * g_z_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_xzzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyyyy_0[i] = 2.0 * g_z_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyy_1[i] * fz_be_0 + g_zz_0_yyyyyyy_1[i] * wa_z[i];

        g_zzz_0_yyyyyyz_0[i] = 2.0 * g_z_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyz_1[i] * fz_be_0 + g_zz_0_yyyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyyz_1[i] * wa_z[i];

        g_zzz_0_yyyyyzz_0[i] = 2.0 * g_z_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_yyyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyyzz_1[i] * wa_z[i];

        g_zzz_0_yyyyzzz_0[i] = 2.0 * g_z_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_yyyyzz_1[i] * fi_acd_0 + g_zz_0_yyyyzzz_1[i] * wa_z[i];

        g_zzz_0_yyyzzzz_0[i] = 2.0 * g_z_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_yyyzzz_1[i] * fi_acd_0 + g_zz_0_yyyzzzz_1[i] * wa_z[i];

        g_zzz_0_yyzzzzz_0[i] = 2.0 * g_z_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_yyzzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzzz_1[i] * wa_z[i];

        g_zzz_0_yzzzzzz_0[i] = 2.0 * g_z_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_yzzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzzz_1[i] * wa_z[i];

        g_zzz_0_zzzzzzz_0[i] = 2.0 * g_z_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_zzzzzzz_1[i] * fz_be_0 + 7.0 * g_zz_0_zzzzzz_1[i] * fi_acd_0 + g_zz_0_zzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

