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

#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsl,
                                 size_t idx_eri_0_ssl,
                                 size_t idx_eri_1_ssl,
                                 size_t idx_eri_1_psk,
                                 size_t idx_eri_1_psl,
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

    /// Set up components of auxilary buffer : SSL

    auto g_0_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ssl);

    auto g_0_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ssl + 1);

    auto g_0_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ssl + 2);

    auto g_0_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ssl + 3);

    auto g_0_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ssl + 4);

    auto g_0_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ssl + 5);

    auto g_0_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ssl + 6);

    auto g_0_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ssl + 7);

    auto g_0_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ssl + 8);

    auto g_0_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ssl + 9);

    auto g_0_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ssl + 10);

    auto g_0_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ssl + 11);

    auto g_0_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ssl + 12);

    auto g_0_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ssl + 13);

    auto g_0_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ssl + 14);

    auto g_0_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 15);

    auto g_0_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 16);

    auto g_0_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 17);

    auto g_0_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 18);

    auto g_0_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 19);

    auto g_0_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 20);

    auto g_0_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 21);

    auto g_0_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 22);

    auto g_0_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 23);

    auto g_0_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 24);

    auto g_0_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 25);

    auto g_0_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 26);

    auto g_0_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 27);

    auto g_0_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 28);

    auto g_0_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 29);

    auto g_0_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 30);

    auto g_0_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 31);

    auto g_0_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 32);

    auto g_0_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 33);

    auto g_0_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 34);

    auto g_0_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 35);

    auto g_0_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 36);

    auto g_0_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 37);

    auto g_0_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 38);

    auto g_0_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 39);

    auto g_0_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 40);

    auto g_0_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 41);

    auto g_0_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 42);

    auto g_0_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 43);

    auto g_0_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 44);

    /// Set up components of auxilary buffer : SSL

    auto g_0_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_ssl);

    auto g_0_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_ssl + 1);

    auto g_0_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_ssl + 2);

    auto g_0_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_ssl + 3);

    auto g_0_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_ssl + 4);

    auto g_0_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_ssl + 5);

    auto g_0_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_ssl + 6);

    auto g_0_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_ssl + 7);

    auto g_0_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_ssl + 8);

    auto g_0_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_ssl + 9);

    auto g_0_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_ssl + 10);

    auto g_0_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_ssl + 11);

    auto g_0_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_ssl + 12);

    auto g_0_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_ssl + 13);

    auto g_0_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_ssl + 14);

    auto g_0_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 15);

    auto g_0_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 16);

    auto g_0_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 17);

    auto g_0_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 18);

    auto g_0_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 19);

    auto g_0_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 20);

    auto g_0_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 21);

    auto g_0_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 22);

    auto g_0_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 23);

    auto g_0_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 24);

    auto g_0_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 25);

    auto g_0_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 26);

    auto g_0_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 27);

    auto g_0_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 28);

    auto g_0_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 29);

    auto g_0_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 30);

    auto g_0_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 31);

    auto g_0_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 32);

    auto g_0_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 33);

    auto g_0_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 34);

    auto g_0_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 35);

    auto g_0_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 36);

    auto g_0_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 37);

    auto g_0_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 38);

    auto g_0_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 39);

    auto g_0_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 40);

    auto g_0_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 41);

    auto g_0_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 42);

    auto g_0_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 43);

    auto g_0_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 44);

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

    /// Set up components of auxilary buffer : PSL

    auto g_x_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_psl);

    auto g_x_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_psl + 1);

    auto g_x_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_psl + 2);

    auto g_x_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_psl + 3);

    auto g_x_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_psl + 4);

    auto g_x_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_psl + 5);

    auto g_x_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_psl + 6);

    auto g_x_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_psl + 7);

    auto g_x_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_psl + 8);

    auto g_x_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_psl + 9);

    auto g_x_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_psl + 10);

    auto g_x_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_psl + 11);

    auto g_x_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_psl + 12);

    auto g_x_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_psl + 13);

    auto g_x_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_psl + 14);

    auto g_x_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_psl + 15);

    auto g_x_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_psl + 16);

    auto g_x_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_psl + 17);

    auto g_x_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_psl + 18);

    auto g_x_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_psl + 19);

    auto g_x_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_psl + 20);

    auto g_x_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 21);

    auto g_x_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 22);

    auto g_x_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 23);

    auto g_x_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 24);

    auto g_x_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 25);

    auto g_x_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 26);

    auto g_x_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 27);

    auto g_x_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 28);

    auto g_x_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 29);

    auto g_x_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 30);

    auto g_x_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 31);

    auto g_x_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 32);

    auto g_x_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 33);

    auto g_x_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 34);

    auto g_x_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 35);

    auto g_x_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 36);

    auto g_x_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 37);

    auto g_x_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 38);

    auto g_x_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 39);

    auto g_x_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 40);

    auto g_x_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 41);

    auto g_x_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 42);

    auto g_x_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 43);

    auto g_x_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 44);

    auto g_y_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_psl + 45);

    auto g_y_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_psl + 46);

    auto g_y_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_psl + 47);

    auto g_y_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_psl + 48);

    auto g_y_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_psl + 49);

    auto g_y_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_psl + 50);

    auto g_y_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_psl + 51);

    auto g_y_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_psl + 52);

    auto g_y_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_psl + 53);

    auto g_y_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_psl + 54);

    auto g_y_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_psl + 55);

    auto g_y_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_psl + 56);

    auto g_y_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_psl + 57);

    auto g_y_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_psl + 58);

    auto g_y_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_psl + 59);

    auto g_y_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_psl + 60);

    auto g_y_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_psl + 61);

    auto g_y_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_psl + 62);

    auto g_y_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_psl + 63);

    auto g_y_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_psl + 64);

    auto g_y_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_psl + 65);

    auto g_y_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 66);

    auto g_y_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 67);

    auto g_y_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 68);

    auto g_y_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 69);

    auto g_y_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 70);

    auto g_y_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 71);

    auto g_y_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 72);

    auto g_y_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 73);

    auto g_y_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 74);

    auto g_y_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 75);

    auto g_y_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 76);

    auto g_y_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 77);

    auto g_y_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 78);

    auto g_y_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 79);

    auto g_y_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 80);

    auto g_y_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 81);

    auto g_y_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 82);

    auto g_y_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 83);

    auto g_y_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 84);

    auto g_y_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 85);

    auto g_y_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 86);

    auto g_y_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 87);

    auto g_y_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 88);

    auto g_y_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 89);

    auto g_z_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_psl + 90);

    auto g_z_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_psl + 91);

    auto g_z_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_psl + 92);

    auto g_z_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_psl + 93);

    auto g_z_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_psl + 94);

    auto g_z_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_psl + 95);

    auto g_z_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_psl + 96);

    auto g_z_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_psl + 97);

    auto g_z_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_psl + 98);

    auto g_z_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_psl + 99);

    auto g_z_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_psl + 100);

    auto g_z_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_psl + 101);

    auto g_z_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_psl + 102);

    auto g_z_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_psl + 103);

    auto g_z_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_psl + 104);

    auto g_z_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_psl + 105);

    auto g_z_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_psl + 106);

    auto g_z_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_psl + 107);

    auto g_z_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_psl + 108);

    auto g_z_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_psl + 109);

    auto g_z_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_psl + 110);

    auto g_z_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 111);

    auto g_z_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 112);

    auto g_z_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 113);

    auto g_z_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 114);

    auto g_z_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 115);

    auto g_z_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 116);

    auto g_z_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 117);

    auto g_z_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 118);

    auto g_z_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 119);

    auto g_z_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 120);

    auto g_z_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 121);

    auto g_z_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 122);

    auto g_z_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 123);

    auto g_z_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 124);

    auto g_z_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 125);

    auto g_z_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_psl + 126);

    auto g_z_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_psl + 127);

    auto g_z_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_psl + 128);

    auto g_z_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_psl + 129);

    auto g_z_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_psl + 130);

    auto g_z_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_psl + 131);

    auto g_z_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 132);

    auto g_z_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 133);

    auto g_z_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_psl + 134);

    /// Set up 0-45 components of targeted buffer : DSL

    auto g_xx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl);

    auto g_xx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 1);

    auto g_xx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 2);

    auto g_xx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 3);

    auto g_xx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 4);

    auto g_xx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 5);

    auto g_xx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 6);

    auto g_xx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 7);

    auto g_xx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 8);

    auto g_xx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 9);

    auto g_xx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 10);

    auto g_xx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 11);

    auto g_xx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 12);

    auto g_xx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 13);

    auto g_xx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 14);

    auto g_xx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 15);

    auto g_xx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 16);

    auto g_xx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 17);

    auto g_xx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 18);

    auto g_xx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 19);

    auto g_xx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 20);

    auto g_xx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 21);

    auto g_xx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 22);

    auto g_xx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 23);

    auto g_xx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 24);

    auto g_xx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 25);

    auto g_xx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 26);

    auto g_xx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 27);

    auto g_xx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 28);

    auto g_xx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 29);

    auto g_xx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 30);

    auto g_xx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 31);

    auto g_xx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 32);

    auto g_xx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 33);

    auto g_xx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 34);

    auto g_xx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 35);

    auto g_xx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 36);

    auto g_xx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 37);

    auto g_xx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 38);

    auto g_xx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 39);

    auto g_xx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 40);

    auto g_xx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 41);

    auto g_xx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 42);

    auto g_xx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 43);

    auto g_xx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 44);

    #pragma omp simd aligned(g_0_0_xxxxxxxx_0, g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxy_0, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxz_0, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxyy_0, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyz_0, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxzz_0, g_0_0_xxxxxxzz_1, g_0_0_xxxxxyyy_0, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyz_0, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyzz_0, g_0_0_xxxxxyzz_1, g_0_0_xxxxxzzz_0, g_0_0_xxxxxzzz_1, g_0_0_xxxxyyyy_0, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyz_0, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyzz_0, g_0_0_xxxxyyzz_1, g_0_0_xxxxyzzz_0, g_0_0_xxxxyzzz_1, g_0_0_xxxxzzzz_0, g_0_0_xxxxzzzz_1, g_0_0_xxxyyyyy_0, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyz_0, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyzz_0, g_0_0_xxxyyyzz_1, g_0_0_xxxyyzzz_0, g_0_0_xxxyyzzz_1, g_0_0_xxxyzzzz_0, g_0_0_xxxyzzzz_1, g_0_0_xxxzzzzz_0, g_0_0_xxxzzzzz_1, g_0_0_xxyyyyyy_0, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyz_0, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyzz_0, g_0_0_xxyyyyzz_1, g_0_0_xxyyyzzz_0, g_0_0_xxyyyzzz_1, g_0_0_xxyyzzzz_0, g_0_0_xxyyzzzz_1, g_0_0_xxyzzzzz_0, g_0_0_xxyzzzzz_1, g_0_0_xxzzzzzz_0, g_0_0_xxzzzzzz_1, g_0_0_xyyyyyyy_0, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyz_0, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyzz_0, g_0_0_xyyyyyzz_1, g_0_0_xyyyyzzz_0, g_0_0_xyyyyzzz_1, g_0_0_xyyyzzzz_0, g_0_0_xyyyzzzz_1, g_0_0_xyyzzzzz_0, g_0_0_xyyzzzzz_1, g_0_0_xyzzzzzz_0, g_0_0_xyzzzzzz_1, g_0_0_xzzzzzzz_0, g_0_0_xzzzzzzz_1, g_0_0_yyyyyyyy_0, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyz_0, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyzz_0, g_0_0_yyyyyyzz_1, g_0_0_yyyyyzzz_0, g_0_0_yyyyyzzz_1, g_0_0_yyyyzzzz_0, g_0_0_yyyyzzzz_1, g_0_0_yyyzzzzz_0, g_0_0_yyyzzzzz_1, g_0_0_yyzzzzzz_0, g_0_0_yyzzzzzz_1, g_0_0_yzzzzzzz_0, g_0_0_yzzzzzzz_1, g_0_0_zzzzzzzz_0, g_0_0_zzzzzzzz_1, g_x_0_xxxxxxx_1, g_x_0_xxxxxxxx_1, g_x_0_xxxxxxxy_1, g_x_0_xxxxxxxz_1, g_x_0_xxxxxxy_1, g_x_0_xxxxxxyy_1, g_x_0_xxxxxxyz_1, g_x_0_xxxxxxz_1, g_x_0_xxxxxxzz_1, g_x_0_xxxxxyy_1, g_x_0_xxxxxyyy_1, g_x_0_xxxxxyyz_1, g_x_0_xxxxxyz_1, g_x_0_xxxxxyzz_1, g_x_0_xxxxxzz_1, g_x_0_xxxxxzzz_1, g_x_0_xxxxyyy_1, g_x_0_xxxxyyyy_1, g_x_0_xxxxyyyz_1, g_x_0_xxxxyyz_1, g_x_0_xxxxyyzz_1, g_x_0_xxxxyzz_1, g_x_0_xxxxyzzz_1, g_x_0_xxxxzzz_1, g_x_0_xxxxzzzz_1, g_x_0_xxxyyyy_1, g_x_0_xxxyyyyy_1, g_x_0_xxxyyyyz_1, g_x_0_xxxyyyz_1, g_x_0_xxxyyyzz_1, g_x_0_xxxyyzz_1, g_x_0_xxxyyzzz_1, g_x_0_xxxyzzz_1, g_x_0_xxxyzzzz_1, g_x_0_xxxzzzz_1, g_x_0_xxxzzzzz_1, g_x_0_xxyyyyy_1, g_x_0_xxyyyyyy_1, g_x_0_xxyyyyyz_1, g_x_0_xxyyyyz_1, g_x_0_xxyyyyzz_1, g_x_0_xxyyyzz_1, g_x_0_xxyyyzzz_1, g_x_0_xxyyzzz_1, g_x_0_xxyyzzzz_1, g_x_0_xxyzzzz_1, g_x_0_xxyzzzzz_1, g_x_0_xxzzzzz_1, g_x_0_xxzzzzzz_1, g_x_0_xyyyyyy_1, g_x_0_xyyyyyyy_1, g_x_0_xyyyyyyz_1, g_x_0_xyyyyyz_1, g_x_0_xyyyyyzz_1, g_x_0_xyyyyzz_1, g_x_0_xyyyyzzz_1, g_x_0_xyyyzzz_1, g_x_0_xyyyzzzz_1, g_x_0_xyyzzzz_1, g_x_0_xyyzzzzz_1, g_x_0_xyzzzzz_1, g_x_0_xyzzzzzz_1, g_x_0_xzzzzzz_1, g_x_0_xzzzzzzz_1, g_x_0_yyyyyyy_1, g_x_0_yyyyyyyy_1, g_x_0_yyyyyyyz_1, g_x_0_yyyyyyz_1, g_x_0_yyyyyyzz_1, g_x_0_yyyyyzz_1, g_x_0_yyyyyzzz_1, g_x_0_yyyyzzz_1, g_x_0_yyyyzzzz_1, g_x_0_yyyzzzz_1, g_x_0_yyyzzzzz_1, g_x_0_yyzzzzz_1, g_x_0_yyzzzzzz_1, g_x_0_yzzzzzz_1, g_x_0_yzzzzzzz_1, g_x_0_zzzzzzz_1, g_x_0_zzzzzzzz_1, g_xx_0_xxxxxxxx_0, g_xx_0_xxxxxxxy_0, g_xx_0_xxxxxxxz_0, g_xx_0_xxxxxxyy_0, g_xx_0_xxxxxxyz_0, g_xx_0_xxxxxxzz_0, g_xx_0_xxxxxyyy_0, g_xx_0_xxxxxyyz_0, g_xx_0_xxxxxyzz_0, g_xx_0_xxxxxzzz_0, g_xx_0_xxxxyyyy_0, g_xx_0_xxxxyyyz_0, g_xx_0_xxxxyyzz_0, g_xx_0_xxxxyzzz_0, g_xx_0_xxxxzzzz_0, g_xx_0_xxxyyyyy_0, g_xx_0_xxxyyyyz_0, g_xx_0_xxxyyyzz_0, g_xx_0_xxxyyzzz_0, g_xx_0_xxxyzzzz_0, g_xx_0_xxxzzzzz_0, g_xx_0_xxyyyyyy_0, g_xx_0_xxyyyyyz_0, g_xx_0_xxyyyyzz_0, g_xx_0_xxyyyzzz_0, g_xx_0_xxyyzzzz_0, g_xx_0_xxyzzzzz_0, g_xx_0_xxzzzzzz_0, g_xx_0_xyyyyyyy_0, g_xx_0_xyyyyyyz_0, g_xx_0_xyyyyyzz_0, g_xx_0_xyyyyzzz_0, g_xx_0_xyyyzzzz_0, g_xx_0_xyyzzzzz_0, g_xx_0_xyzzzzzz_0, g_xx_0_xzzzzzzz_0, g_xx_0_yyyyyyyy_0, g_xx_0_yyyyyyyz_0, g_xx_0_yyyyyyzz_0, g_xx_0_yyyyyzzz_0, g_xx_0_yyyyzzzz_0, g_xx_0_yyyzzzzz_0, g_xx_0_yyzzzzzz_0, g_xx_0_yzzzzzzz_0, g_xx_0_zzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xxxxxxxx_0[i] = g_0_0_xxxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxxx_1[i] * fz_be_0 + 8.0 * g_x_0_xxxxxxx_1[i] * fi_acd_0 + g_x_0_xxxxxxxx_1[i] * wa_x[i];

        g_xx_0_xxxxxxxy_0[i] = g_0_0_xxxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_x_0_xxxxxxy_1[i] * fi_acd_0 + g_x_0_xxxxxxxy_1[i] * wa_x[i];

        g_xx_0_xxxxxxxz_0[i] = g_0_0_xxxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_x_0_xxxxxxz_1[i] * fi_acd_0 + g_x_0_xxxxxxxz_1[i] * wa_x[i];

        g_xx_0_xxxxxxyy_0[i] = g_0_0_xxxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxyy_1[i] * fi_acd_0 + g_x_0_xxxxxxyy_1[i] * wa_x[i];

        g_xx_0_xxxxxxyz_0[i] = g_0_0_xxxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxyz_1[i] * fi_acd_0 + g_x_0_xxxxxxyz_1[i] * wa_x[i];

        g_xx_0_xxxxxxzz_0[i] = g_0_0_xxxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxzz_1[i] * fi_acd_0 + g_x_0_xxxxxxzz_1[i] * wa_x[i];

        g_xx_0_xxxxxyyy_0[i] = g_0_0_xxxxxyyy_0[i] * fbe_0 - g_0_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyyy_1[i] * fi_acd_0 + g_x_0_xxxxxyyy_1[i] * wa_x[i];

        g_xx_0_xxxxxyyz_0[i] = g_0_0_xxxxxyyz_0[i] * fbe_0 - g_0_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyyz_1[i] * fi_acd_0 + g_x_0_xxxxxyyz_1[i] * wa_x[i];

        g_xx_0_xxxxxyzz_0[i] = g_0_0_xxxxxyzz_0[i] * fbe_0 - g_0_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyzz_1[i] * fi_acd_0 + g_x_0_xxxxxyzz_1[i] * wa_x[i];

        g_xx_0_xxxxxzzz_0[i] = g_0_0_xxxxxzzz_0[i] * fbe_0 - g_0_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxzzz_1[i] * fi_acd_0 + g_x_0_xxxxxzzz_1[i] * wa_x[i];

        g_xx_0_xxxxyyyy_0[i] = g_0_0_xxxxyyyy_0[i] * fbe_0 - g_0_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyyy_1[i] * fi_acd_0 + g_x_0_xxxxyyyy_1[i] * wa_x[i];

        g_xx_0_xxxxyyyz_0[i] = g_0_0_xxxxyyyz_0[i] * fbe_0 - g_0_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyyz_1[i] * fi_acd_0 + g_x_0_xxxxyyyz_1[i] * wa_x[i];

        g_xx_0_xxxxyyzz_0[i] = g_0_0_xxxxyyzz_0[i] * fbe_0 - g_0_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyzz_1[i] * fi_acd_0 + g_x_0_xxxxyyzz_1[i] * wa_x[i];

        g_xx_0_xxxxyzzz_0[i] = g_0_0_xxxxyzzz_0[i] * fbe_0 - g_0_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyzzz_1[i] * fi_acd_0 + g_x_0_xxxxyzzz_1[i] * wa_x[i];

        g_xx_0_xxxxzzzz_0[i] = g_0_0_xxxxzzzz_0[i] * fbe_0 - g_0_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxzzzz_1[i] * fi_acd_0 + g_x_0_xxxxzzzz_1[i] * wa_x[i];

        g_xx_0_xxxyyyyy_0[i] = g_0_0_xxxyyyyy_0[i] * fbe_0 - g_0_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyyy_1[i] * fi_acd_0 + g_x_0_xxxyyyyy_1[i] * wa_x[i];

        g_xx_0_xxxyyyyz_0[i] = g_0_0_xxxyyyyz_0[i] * fbe_0 - g_0_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyyz_1[i] * fi_acd_0 + g_x_0_xxxyyyyz_1[i] * wa_x[i];

        g_xx_0_xxxyyyzz_0[i] = g_0_0_xxxyyyzz_0[i] * fbe_0 - g_0_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyzz_1[i] * fi_acd_0 + g_x_0_xxxyyyzz_1[i] * wa_x[i];

        g_xx_0_xxxyyzzz_0[i] = g_0_0_xxxyyzzz_0[i] * fbe_0 - g_0_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyzzz_1[i] * fi_acd_0 + g_x_0_xxxyyzzz_1[i] * wa_x[i];

        g_xx_0_xxxyzzzz_0[i] = g_0_0_xxxyzzzz_0[i] * fbe_0 - g_0_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyzzzz_1[i] * fi_acd_0 + g_x_0_xxxyzzzz_1[i] * wa_x[i];

        g_xx_0_xxxzzzzz_0[i] = g_0_0_xxxzzzzz_0[i] * fbe_0 - g_0_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxzzzzz_1[i] * fi_acd_0 + g_x_0_xxxzzzzz_1[i] * wa_x[i];

        g_xx_0_xxyyyyyy_0[i] = g_0_0_xxyyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyyy_1[i] * fi_acd_0 + g_x_0_xxyyyyyy_1[i] * wa_x[i];

        g_xx_0_xxyyyyyz_0[i] = g_0_0_xxyyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyyz_1[i] * fi_acd_0 + g_x_0_xxyyyyyz_1[i] * wa_x[i];

        g_xx_0_xxyyyyzz_0[i] = g_0_0_xxyyyyzz_0[i] * fbe_0 - g_0_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyzz_1[i] * fi_acd_0 + g_x_0_xxyyyyzz_1[i] * wa_x[i];

        g_xx_0_xxyyyzzz_0[i] = g_0_0_xxyyyzzz_0[i] * fbe_0 - g_0_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyzzz_1[i] * fi_acd_0 + g_x_0_xxyyyzzz_1[i] * wa_x[i];

        g_xx_0_xxyyzzzz_0[i] = g_0_0_xxyyzzzz_0[i] * fbe_0 - g_0_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyzzzz_1[i] * fi_acd_0 + g_x_0_xxyyzzzz_1[i] * wa_x[i];

        g_xx_0_xxyzzzzz_0[i] = g_0_0_xxyzzzzz_0[i] * fbe_0 - g_0_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyzzzzz_1[i] * fi_acd_0 + g_x_0_xxyzzzzz_1[i] * wa_x[i];

        g_xx_0_xxzzzzzz_0[i] = g_0_0_xxzzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xzzzzzz_1[i] * fi_acd_0 + g_x_0_xxzzzzzz_1[i] * wa_x[i];

        g_xx_0_xyyyyyyy_0[i] = g_0_0_xyyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyyy_1[i] * fz_be_0 + g_x_0_yyyyyyy_1[i] * fi_acd_0 + g_x_0_xyyyyyyy_1[i] * wa_x[i];

        g_xx_0_xyyyyyyz_0[i] = g_0_0_xyyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyyz_1[i] * fz_be_0 + g_x_0_yyyyyyz_1[i] * fi_acd_0 + g_x_0_xyyyyyyz_1[i] * wa_x[i];

        g_xx_0_xyyyyyzz_0[i] = g_0_0_xyyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyyzz_1[i] * fz_be_0 + g_x_0_yyyyyzz_1[i] * fi_acd_0 + g_x_0_xyyyyyzz_1[i] * wa_x[i];

        g_xx_0_xyyyyzzz_0[i] = g_0_0_xyyyyzzz_0[i] * fbe_0 - g_0_0_xyyyyzzz_1[i] * fz_be_0 + g_x_0_yyyyzzz_1[i] * fi_acd_0 + g_x_0_xyyyyzzz_1[i] * wa_x[i];

        g_xx_0_xyyyzzzz_0[i] = g_0_0_xyyyzzzz_0[i] * fbe_0 - g_0_0_xyyyzzzz_1[i] * fz_be_0 + g_x_0_yyyzzzz_1[i] * fi_acd_0 + g_x_0_xyyyzzzz_1[i] * wa_x[i];

        g_xx_0_xyyzzzzz_0[i] = g_0_0_xyyzzzzz_0[i] * fbe_0 - g_0_0_xyyzzzzz_1[i] * fz_be_0 + g_x_0_yyzzzzz_1[i] * fi_acd_0 + g_x_0_xyyzzzzz_1[i] * wa_x[i];

        g_xx_0_xyzzzzzz_0[i] = g_0_0_xyzzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzzz_1[i] * fz_be_0 + g_x_0_yzzzzzz_1[i] * fi_acd_0 + g_x_0_xyzzzzzz_1[i] * wa_x[i];

        g_xx_0_xzzzzzzz_0[i] = g_0_0_xzzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzzz_1[i] * fz_be_0 + g_x_0_zzzzzzz_1[i] * fi_acd_0 + g_x_0_xzzzzzzz_1[i] * wa_x[i];

        g_xx_0_yyyyyyyy_0[i] = g_0_0_yyyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyyy_1[i] * fz_be_0 + g_x_0_yyyyyyyy_1[i] * wa_x[i];

        g_xx_0_yyyyyyyz_0[i] = g_0_0_yyyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyyz_1[i] * fz_be_0 + g_x_0_yyyyyyyz_1[i] * wa_x[i];

        g_xx_0_yyyyyyzz_0[i] = g_0_0_yyyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyyzz_1[i] * fz_be_0 + g_x_0_yyyyyyzz_1[i] * wa_x[i];

        g_xx_0_yyyyyzzz_0[i] = g_0_0_yyyyyzzz_0[i] * fbe_0 - g_0_0_yyyyyzzz_1[i] * fz_be_0 + g_x_0_yyyyyzzz_1[i] * wa_x[i];

        g_xx_0_yyyyzzzz_0[i] = g_0_0_yyyyzzzz_0[i] * fbe_0 - g_0_0_yyyyzzzz_1[i] * fz_be_0 + g_x_0_yyyyzzzz_1[i] * wa_x[i];

        g_xx_0_yyyzzzzz_0[i] = g_0_0_yyyzzzzz_0[i] * fbe_0 - g_0_0_yyyzzzzz_1[i] * fz_be_0 + g_x_0_yyyzzzzz_1[i] * wa_x[i];

        g_xx_0_yyzzzzzz_0[i] = g_0_0_yyzzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzzz_1[i] * fz_be_0 + g_x_0_yyzzzzzz_1[i] * wa_x[i];

        g_xx_0_yzzzzzzz_0[i] = g_0_0_yzzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzzz_1[i] * fz_be_0 + g_x_0_yzzzzzzz_1[i] * wa_x[i];

        g_xx_0_zzzzzzzz_0[i] = g_0_0_zzzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzzz_1[i] * fz_be_0 + g_x_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 45-90 components of targeted buffer : DSL

    auto g_xy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl + 45);

    auto g_xy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 46);

    auto g_xy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 47);

    auto g_xy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 48);

    auto g_xy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 49);

    auto g_xy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 50);

    auto g_xy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 51);

    auto g_xy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 52);

    auto g_xy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 53);

    auto g_xy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 54);

    auto g_xy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 55);

    auto g_xy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 56);

    auto g_xy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 57);

    auto g_xy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 58);

    auto g_xy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 59);

    auto g_xy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 60);

    auto g_xy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 61);

    auto g_xy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 62);

    auto g_xy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 63);

    auto g_xy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 64);

    auto g_xy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 65);

    auto g_xy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 66);

    auto g_xy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 67);

    auto g_xy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 68);

    auto g_xy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 69);

    auto g_xy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 70);

    auto g_xy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 71);

    auto g_xy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 72);

    auto g_xy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 73);

    auto g_xy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 74);

    auto g_xy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 75);

    auto g_xy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 76);

    auto g_xy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 77);

    auto g_xy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 78);

    auto g_xy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 79);

    auto g_xy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 80);

    auto g_xy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 81);

    auto g_xy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 82);

    auto g_xy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 83);

    auto g_xy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 84);

    auto g_xy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 85);

    auto g_xy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 86);

    auto g_xy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 87);

    auto g_xy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 88);

    auto g_xy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 89);

    #pragma omp simd aligned(g_x_0_xxxxxxxx_1, g_x_0_xxxxxxxz_1, g_x_0_xxxxxxzz_1, g_x_0_xxxxxzzz_1, g_x_0_xxxxzzzz_1, g_x_0_xxxzzzzz_1, g_x_0_xxzzzzzz_1, g_x_0_xzzzzzzz_1, g_xy_0_xxxxxxxx_0, g_xy_0_xxxxxxxy_0, g_xy_0_xxxxxxxz_0, g_xy_0_xxxxxxyy_0, g_xy_0_xxxxxxyz_0, g_xy_0_xxxxxxzz_0, g_xy_0_xxxxxyyy_0, g_xy_0_xxxxxyyz_0, g_xy_0_xxxxxyzz_0, g_xy_0_xxxxxzzz_0, g_xy_0_xxxxyyyy_0, g_xy_0_xxxxyyyz_0, g_xy_0_xxxxyyzz_0, g_xy_0_xxxxyzzz_0, g_xy_0_xxxxzzzz_0, g_xy_0_xxxyyyyy_0, g_xy_0_xxxyyyyz_0, g_xy_0_xxxyyyzz_0, g_xy_0_xxxyyzzz_0, g_xy_0_xxxyzzzz_0, g_xy_0_xxxzzzzz_0, g_xy_0_xxyyyyyy_0, g_xy_0_xxyyyyyz_0, g_xy_0_xxyyyyzz_0, g_xy_0_xxyyyzzz_0, g_xy_0_xxyyzzzz_0, g_xy_0_xxyzzzzz_0, g_xy_0_xxzzzzzz_0, g_xy_0_xyyyyyyy_0, g_xy_0_xyyyyyyz_0, g_xy_0_xyyyyyzz_0, g_xy_0_xyyyyzzz_0, g_xy_0_xyyyzzzz_0, g_xy_0_xyyzzzzz_0, g_xy_0_xyzzzzzz_0, g_xy_0_xzzzzzzz_0, g_xy_0_yyyyyyyy_0, g_xy_0_yyyyyyyz_0, g_xy_0_yyyyyyzz_0, g_xy_0_yyyyyzzz_0, g_xy_0_yyyyzzzz_0, g_xy_0_yyyzzzzz_0, g_xy_0_yyzzzzzz_0, g_xy_0_yzzzzzzz_0, g_xy_0_zzzzzzzz_0, g_y_0_xxxxxxxy_1, g_y_0_xxxxxxy_1, g_y_0_xxxxxxyy_1, g_y_0_xxxxxxyz_1, g_y_0_xxxxxyy_1, g_y_0_xxxxxyyy_1, g_y_0_xxxxxyyz_1, g_y_0_xxxxxyz_1, g_y_0_xxxxxyzz_1, g_y_0_xxxxyyy_1, g_y_0_xxxxyyyy_1, g_y_0_xxxxyyyz_1, g_y_0_xxxxyyz_1, g_y_0_xxxxyyzz_1, g_y_0_xxxxyzz_1, g_y_0_xxxxyzzz_1, g_y_0_xxxyyyy_1, g_y_0_xxxyyyyy_1, g_y_0_xxxyyyyz_1, g_y_0_xxxyyyz_1, g_y_0_xxxyyyzz_1, g_y_0_xxxyyzz_1, g_y_0_xxxyyzzz_1, g_y_0_xxxyzzz_1, g_y_0_xxxyzzzz_1, g_y_0_xxyyyyy_1, g_y_0_xxyyyyyy_1, g_y_0_xxyyyyyz_1, g_y_0_xxyyyyz_1, g_y_0_xxyyyyzz_1, g_y_0_xxyyyzz_1, g_y_0_xxyyyzzz_1, g_y_0_xxyyzzz_1, g_y_0_xxyyzzzz_1, g_y_0_xxyzzzz_1, g_y_0_xxyzzzzz_1, g_y_0_xyyyyyy_1, g_y_0_xyyyyyyy_1, g_y_0_xyyyyyyz_1, g_y_0_xyyyyyz_1, g_y_0_xyyyyyzz_1, g_y_0_xyyyyzz_1, g_y_0_xyyyyzzz_1, g_y_0_xyyyzzz_1, g_y_0_xyyyzzzz_1, g_y_0_xyyzzzz_1, g_y_0_xyyzzzzz_1, g_y_0_xyzzzzz_1, g_y_0_xyzzzzzz_1, g_y_0_yyyyyyy_1, g_y_0_yyyyyyyy_1, g_y_0_yyyyyyyz_1, g_y_0_yyyyyyz_1, g_y_0_yyyyyyzz_1, g_y_0_yyyyyzz_1, g_y_0_yyyyyzzz_1, g_y_0_yyyyzzz_1, g_y_0_yyyyzzzz_1, g_y_0_yyyzzzz_1, g_y_0_yyyzzzzz_1, g_y_0_yyzzzzz_1, g_y_0_yyzzzzzz_1, g_y_0_yzzzzzz_1, g_y_0_yzzzzzzz_1, g_y_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xxxxxxxx_0[i] = g_x_0_xxxxxxxx_1[i] * wa_y[i];

        g_xy_0_xxxxxxxy_0[i] = 7.0 * g_y_0_xxxxxxy_1[i] * fi_acd_0 + g_y_0_xxxxxxxy_1[i] * wa_x[i];

        g_xy_0_xxxxxxxz_0[i] = g_x_0_xxxxxxxz_1[i] * wa_y[i];

        g_xy_0_xxxxxxyy_0[i] = 6.0 * g_y_0_xxxxxyy_1[i] * fi_acd_0 + g_y_0_xxxxxxyy_1[i] * wa_x[i];

        g_xy_0_xxxxxxyz_0[i] = 6.0 * g_y_0_xxxxxyz_1[i] * fi_acd_0 + g_y_0_xxxxxxyz_1[i] * wa_x[i];

        g_xy_0_xxxxxxzz_0[i] = g_x_0_xxxxxxzz_1[i] * wa_y[i];

        g_xy_0_xxxxxyyy_0[i] = 5.0 * g_y_0_xxxxyyy_1[i] * fi_acd_0 + g_y_0_xxxxxyyy_1[i] * wa_x[i];

        g_xy_0_xxxxxyyz_0[i] = 5.0 * g_y_0_xxxxyyz_1[i] * fi_acd_0 + g_y_0_xxxxxyyz_1[i] * wa_x[i];

        g_xy_0_xxxxxyzz_0[i] = 5.0 * g_y_0_xxxxyzz_1[i] * fi_acd_0 + g_y_0_xxxxxyzz_1[i] * wa_x[i];

        g_xy_0_xxxxxzzz_0[i] = g_x_0_xxxxxzzz_1[i] * wa_y[i];

        g_xy_0_xxxxyyyy_0[i] = 4.0 * g_y_0_xxxyyyy_1[i] * fi_acd_0 + g_y_0_xxxxyyyy_1[i] * wa_x[i];

        g_xy_0_xxxxyyyz_0[i] = 4.0 * g_y_0_xxxyyyz_1[i] * fi_acd_0 + g_y_0_xxxxyyyz_1[i] * wa_x[i];

        g_xy_0_xxxxyyzz_0[i] = 4.0 * g_y_0_xxxyyzz_1[i] * fi_acd_0 + g_y_0_xxxxyyzz_1[i] * wa_x[i];

        g_xy_0_xxxxyzzz_0[i] = 4.0 * g_y_0_xxxyzzz_1[i] * fi_acd_0 + g_y_0_xxxxyzzz_1[i] * wa_x[i];

        g_xy_0_xxxxzzzz_0[i] = g_x_0_xxxxzzzz_1[i] * wa_y[i];

        g_xy_0_xxxyyyyy_0[i] = 3.0 * g_y_0_xxyyyyy_1[i] * fi_acd_0 + g_y_0_xxxyyyyy_1[i] * wa_x[i];

        g_xy_0_xxxyyyyz_0[i] = 3.0 * g_y_0_xxyyyyz_1[i] * fi_acd_0 + g_y_0_xxxyyyyz_1[i] * wa_x[i];

        g_xy_0_xxxyyyzz_0[i] = 3.0 * g_y_0_xxyyyzz_1[i] * fi_acd_0 + g_y_0_xxxyyyzz_1[i] * wa_x[i];

        g_xy_0_xxxyyzzz_0[i] = 3.0 * g_y_0_xxyyzzz_1[i] * fi_acd_0 + g_y_0_xxxyyzzz_1[i] * wa_x[i];

        g_xy_0_xxxyzzzz_0[i] = 3.0 * g_y_0_xxyzzzz_1[i] * fi_acd_0 + g_y_0_xxxyzzzz_1[i] * wa_x[i];

        g_xy_0_xxxzzzzz_0[i] = g_x_0_xxxzzzzz_1[i] * wa_y[i];

        g_xy_0_xxyyyyyy_0[i] = 2.0 * g_y_0_xyyyyyy_1[i] * fi_acd_0 + g_y_0_xxyyyyyy_1[i] * wa_x[i];

        g_xy_0_xxyyyyyz_0[i] = 2.0 * g_y_0_xyyyyyz_1[i] * fi_acd_0 + g_y_0_xxyyyyyz_1[i] * wa_x[i];

        g_xy_0_xxyyyyzz_0[i] = 2.0 * g_y_0_xyyyyzz_1[i] * fi_acd_0 + g_y_0_xxyyyyzz_1[i] * wa_x[i];

        g_xy_0_xxyyyzzz_0[i] = 2.0 * g_y_0_xyyyzzz_1[i] * fi_acd_0 + g_y_0_xxyyyzzz_1[i] * wa_x[i];

        g_xy_0_xxyyzzzz_0[i] = 2.0 * g_y_0_xyyzzzz_1[i] * fi_acd_0 + g_y_0_xxyyzzzz_1[i] * wa_x[i];

        g_xy_0_xxyzzzzz_0[i] = 2.0 * g_y_0_xyzzzzz_1[i] * fi_acd_0 + g_y_0_xxyzzzzz_1[i] * wa_x[i];

        g_xy_0_xxzzzzzz_0[i] = g_x_0_xxzzzzzz_1[i] * wa_y[i];

        g_xy_0_xyyyyyyy_0[i] = g_y_0_yyyyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyyyy_1[i] * wa_x[i];

        g_xy_0_xyyyyyyz_0[i] = g_y_0_yyyyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyyyz_1[i] * wa_x[i];

        g_xy_0_xyyyyyzz_0[i] = g_y_0_yyyyyzz_1[i] * fi_acd_0 + g_y_0_xyyyyyzz_1[i] * wa_x[i];

        g_xy_0_xyyyyzzz_0[i] = g_y_0_yyyyzzz_1[i] * fi_acd_0 + g_y_0_xyyyyzzz_1[i] * wa_x[i];

        g_xy_0_xyyyzzzz_0[i] = g_y_0_yyyzzzz_1[i] * fi_acd_0 + g_y_0_xyyyzzzz_1[i] * wa_x[i];

        g_xy_0_xyyzzzzz_0[i] = g_y_0_yyzzzzz_1[i] * fi_acd_0 + g_y_0_xyyzzzzz_1[i] * wa_x[i];

        g_xy_0_xyzzzzzz_0[i] = g_y_0_yzzzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzzzz_1[i] * wa_x[i];

        g_xy_0_xzzzzzzz_0[i] = g_x_0_xzzzzzzz_1[i] * wa_y[i];

        g_xy_0_yyyyyyyy_0[i] = g_y_0_yyyyyyyy_1[i] * wa_x[i];

        g_xy_0_yyyyyyyz_0[i] = g_y_0_yyyyyyyz_1[i] * wa_x[i];

        g_xy_0_yyyyyyzz_0[i] = g_y_0_yyyyyyzz_1[i] * wa_x[i];

        g_xy_0_yyyyyzzz_0[i] = g_y_0_yyyyyzzz_1[i] * wa_x[i];

        g_xy_0_yyyyzzzz_0[i] = g_y_0_yyyyzzzz_1[i] * wa_x[i];

        g_xy_0_yyyzzzzz_0[i] = g_y_0_yyyzzzzz_1[i] * wa_x[i];

        g_xy_0_yyzzzzzz_0[i] = g_y_0_yyzzzzzz_1[i] * wa_x[i];

        g_xy_0_yzzzzzzz_0[i] = g_y_0_yzzzzzzz_1[i] * wa_x[i];

        g_xy_0_zzzzzzzz_0[i] = g_y_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 90-135 components of targeted buffer : DSL

    auto g_xz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl + 90);

    auto g_xz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 91);

    auto g_xz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 92);

    auto g_xz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 93);

    auto g_xz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 94);

    auto g_xz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 95);

    auto g_xz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 96);

    auto g_xz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 97);

    auto g_xz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 98);

    auto g_xz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 99);

    auto g_xz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 100);

    auto g_xz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 101);

    auto g_xz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 102);

    auto g_xz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 103);

    auto g_xz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 104);

    auto g_xz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 105);

    auto g_xz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 106);

    auto g_xz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 107);

    auto g_xz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 108);

    auto g_xz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 109);

    auto g_xz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 110);

    auto g_xz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 111);

    auto g_xz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 112);

    auto g_xz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 113);

    auto g_xz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 114);

    auto g_xz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 115);

    auto g_xz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 116);

    auto g_xz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 117);

    auto g_xz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 118);

    auto g_xz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 119);

    auto g_xz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 120);

    auto g_xz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 121);

    auto g_xz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 122);

    auto g_xz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 123);

    auto g_xz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 124);

    auto g_xz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 125);

    auto g_xz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 126);

    auto g_xz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 127);

    auto g_xz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 128);

    auto g_xz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 129);

    auto g_xz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 130);

    auto g_xz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 131);

    auto g_xz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 132);

    auto g_xz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 133);

    auto g_xz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 134);

    #pragma omp simd aligned(g_x_0_xxxxxxxx_1, g_x_0_xxxxxxxy_1, g_x_0_xxxxxxyy_1, g_x_0_xxxxxyyy_1, g_x_0_xxxxyyyy_1, g_x_0_xxxyyyyy_1, g_x_0_xxyyyyyy_1, g_x_0_xyyyyyyy_1, g_xz_0_xxxxxxxx_0, g_xz_0_xxxxxxxy_0, g_xz_0_xxxxxxxz_0, g_xz_0_xxxxxxyy_0, g_xz_0_xxxxxxyz_0, g_xz_0_xxxxxxzz_0, g_xz_0_xxxxxyyy_0, g_xz_0_xxxxxyyz_0, g_xz_0_xxxxxyzz_0, g_xz_0_xxxxxzzz_0, g_xz_0_xxxxyyyy_0, g_xz_0_xxxxyyyz_0, g_xz_0_xxxxyyzz_0, g_xz_0_xxxxyzzz_0, g_xz_0_xxxxzzzz_0, g_xz_0_xxxyyyyy_0, g_xz_0_xxxyyyyz_0, g_xz_0_xxxyyyzz_0, g_xz_0_xxxyyzzz_0, g_xz_0_xxxyzzzz_0, g_xz_0_xxxzzzzz_0, g_xz_0_xxyyyyyy_0, g_xz_0_xxyyyyyz_0, g_xz_0_xxyyyyzz_0, g_xz_0_xxyyyzzz_0, g_xz_0_xxyyzzzz_0, g_xz_0_xxyzzzzz_0, g_xz_0_xxzzzzzz_0, g_xz_0_xyyyyyyy_0, g_xz_0_xyyyyyyz_0, g_xz_0_xyyyyyzz_0, g_xz_0_xyyyyzzz_0, g_xz_0_xyyyzzzz_0, g_xz_0_xyyzzzzz_0, g_xz_0_xyzzzzzz_0, g_xz_0_xzzzzzzz_0, g_xz_0_yyyyyyyy_0, g_xz_0_yyyyyyyz_0, g_xz_0_yyyyyyzz_0, g_xz_0_yyyyyzzz_0, g_xz_0_yyyyzzzz_0, g_xz_0_yyyzzzzz_0, g_xz_0_yyzzzzzz_0, g_xz_0_yzzzzzzz_0, g_xz_0_zzzzzzzz_0, g_z_0_xxxxxxxz_1, g_z_0_xxxxxxyz_1, g_z_0_xxxxxxz_1, g_z_0_xxxxxxzz_1, g_z_0_xxxxxyyz_1, g_z_0_xxxxxyz_1, g_z_0_xxxxxyzz_1, g_z_0_xxxxxzz_1, g_z_0_xxxxxzzz_1, g_z_0_xxxxyyyz_1, g_z_0_xxxxyyz_1, g_z_0_xxxxyyzz_1, g_z_0_xxxxyzz_1, g_z_0_xxxxyzzz_1, g_z_0_xxxxzzz_1, g_z_0_xxxxzzzz_1, g_z_0_xxxyyyyz_1, g_z_0_xxxyyyz_1, g_z_0_xxxyyyzz_1, g_z_0_xxxyyzz_1, g_z_0_xxxyyzzz_1, g_z_0_xxxyzzz_1, g_z_0_xxxyzzzz_1, g_z_0_xxxzzzz_1, g_z_0_xxxzzzzz_1, g_z_0_xxyyyyyz_1, g_z_0_xxyyyyz_1, g_z_0_xxyyyyzz_1, g_z_0_xxyyyzz_1, g_z_0_xxyyyzzz_1, g_z_0_xxyyzzz_1, g_z_0_xxyyzzzz_1, g_z_0_xxyzzzz_1, g_z_0_xxyzzzzz_1, g_z_0_xxzzzzz_1, g_z_0_xxzzzzzz_1, g_z_0_xyyyyyyz_1, g_z_0_xyyyyyz_1, g_z_0_xyyyyyzz_1, g_z_0_xyyyyzz_1, g_z_0_xyyyyzzz_1, g_z_0_xyyyzzz_1, g_z_0_xyyyzzzz_1, g_z_0_xyyzzzz_1, g_z_0_xyyzzzzz_1, g_z_0_xyzzzzz_1, g_z_0_xyzzzzzz_1, g_z_0_xzzzzzz_1, g_z_0_xzzzzzzz_1, g_z_0_yyyyyyyy_1, g_z_0_yyyyyyyz_1, g_z_0_yyyyyyz_1, g_z_0_yyyyyyzz_1, g_z_0_yyyyyzz_1, g_z_0_yyyyyzzz_1, g_z_0_yyyyzzz_1, g_z_0_yyyyzzzz_1, g_z_0_yyyzzzz_1, g_z_0_yyyzzzzz_1, g_z_0_yyzzzzz_1, g_z_0_yyzzzzzz_1, g_z_0_yzzzzzz_1, g_z_0_yzzzzzzz_1, g_z_0_zzzzzzz_1, g_z_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xxxxxxxx_0[i] = g_x_0_xxxxxxxx_1[i] * wa_z[i];

        g_xz_0_xxxxxxxy_0[i] = g_x_0_xxxxxxxy_1[i] * wa_z[i];

        g_xz_0_xxxxxxxz_0[i] = 7.0 * g_z_0_xxxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxxxz_1[i] * wa_x[i];

        g_xz_0_xxxxxxyy_0[i] = g_x_0_xxxxxxyy_1[i] * wa_z[i];

        g_xz_0_xxxxxxyz_0[i] = 6.0 * g_z_0_xxxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxxxyz_1[i] * wa_x[i];

        g_xz_0_xxxxxxzz_0[i] = 6.0 * g_z_0_xxxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxxxzz_1[i] * wa_x[i];

        g_xz_0_xxxxxyyy_0[i] = g_x_0_xxxxxyyy_1[i] * wa_z[i];

        g_xz_0_xxxxxyyz_0[i] = 5.0 * g_z_0_xxxxyyz_1[i] * fi_acd_0 + g_z_0_xxxxxyyz_1[i] * wa_x[i];

        g_xz_0_xxxxxyzz_0[i] = 5.0 * g_z_0_xxxxyzz_1[i] * fi_acd_0 + g_z_0_xxxxxyzz_1[i] * wa_x[i];

        g_xz_0_xxxxxzzz_0[i] = 5.0 * g_z_0_xxxxzzz_1[i] * fi_acd_0 + g_z_0_xxxxxzzz_1[i] * wa_x[i];

        g_xz_0_xxxxyyyy_0[i] = g_x_0_xxxxyyyy_1[i] * wa_z[i];

        g_xz_0_xxxxyyyz_0[i] = 4.0 * g_z_0_xxxyyyz_1[i] * fi_acd_0 + g_z_0_xxxxyyyz_1[i] * wa_x[i];

        g_xz_0_xxxxyyzz_0[i] = 4.0 * g_z_0_xxxyyzz_1[i] * fi_acd_0 + g_z_0_xxxxyyzz_1[i] * wa_x[i];

        g_xz_0_xxxxyzzz_0[i] = 4.0 * g_z_0_xxxyzzz_1[i] * fi_acd_0 + g_z_0_xxxxyzzz_1[i] * wa_x[i];

        g_xz_0_xxxxzzzz_0[i] = 4.0 * g_z_0_xxxzzzz_1[i] * fi_acd_0 + g_z_0_xxxxzzzz_1[i] * wa_x[i];

        g_xz_0_xxxyyyyy_0[i] = g_x_0_xxxyyyyy_1[i] * wa_z[i];

        g_xz_0_xxxyyyyz_0[i] = 3.0 * g_z_0_xxyyyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyyz_1[i] * wa_x[i];

        g_xz_0_xxxyyyzz_0[i] = 3.0 * g_z_0_xxyyyzz_1[i] * fi_acd_0 + g_z_0_xxxyyyzz_1[i] * wa_x[i];

        g_xz_0_xxxyyzzz_0[i] = 3.0 * g_z_0_xxyyzzz_1[i] * fi_acd_0 + g_z_0_xxxyyzzz_1[i] * wa_x[i];

        g_xz_0_xxxyzzzz_0[i] = 3.0 * g_z_0_xxyzzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzzz_1[i] * wa_x[i];

        g_xz_0_xxxzzzzz_0[i] = 3.0 * g_z_0_xxzzzzz_1[i] * fi_acd_0 + g_z_0_xxxzzzzz_1[i] * wa_x[i];

        g_xz_0_xxyyyyyy_0[i] = g_x_0_xxyyyyyy_1[i] * wa_z[i];

        g_xz_0_xxyyyyyz_0[i] = 2.0 * g_z_0_xyyyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyyz_1[i] * wa_x[i];

        g_xz_0_xxyyyyzz_0[i] = 2.0 * g_z_0_xyyyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyyzz_1[i] * wa_x[i];

        g_xz_0_xxyyyzzz_0[i] = 2.0 * g_z_0_xyyyzzz_1[i] * fi_acd_0 + g_z_0_xxyyyzzz_1[i] * wa_x[i];

        g_xz_0_xxyyzzzz_0[i] = 2.0 * g_z_0_xyyzzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzzz_1[i] * wa_x[i];

        g_xz_0_xxyzzzzz_0[i] = 2.0 * g_z_0_xyzzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzzz_1[i] * wa_x[i];

        g_xz_0_xxzzzzzz_0[i] = 2.0 * g_z_0_xzzzzzz_1[i] * fi_acd_0 + g_z_0_xxzzzzzz_1[i] * wa_x[i];

        g_xz_0_xyyyyyyy_0[i] = g_x_0_xyyyyyyy_1[i] * wa_z[i];

        g_xz_0_xyyyyyyz_0[i] = g_z_0_yyyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyyz_1[i] * wa_x[i];

        g_xz_0_xyyyyyzz_0[i] = g_z_0_yyyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyyzz_1[i] * wa_x[i];

        g_xz_0_xyyyyzzz_0[i] = g_z_0_yyyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyyzzz_1[i] * wa_x[i];

        g_xz_0_xyyyzzzz_0[i] = g_z_0_yyyzzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzzz_1[i] * wa_x[i];

        g_xz_0_xyyzzzzz_0[i] = g_z_0_yyzzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzzz_1[i] * wa_x[i];

        g_xz_0_xyzzzzzz_0[i] = g_z_0_yzzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzzz_1[i] * wa_x[i];

        g_xz_0_xzzzzzzz_0[i] = g_z_0_zzzzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzzzz_1[i] * wa_x[i];

        g_xz_0_yyyyyyyy_0[i] = g_z_0_yyyyyyyy_1[i] * wa_x[i];

        g_xz_0_yyyyyyyz_0[i] = g_z_0_yyyyyyyz_1[i] * wa_x[i];

        g_xz_0_yyyyyyzz_0[i] = g_z_0_yyyyyyzz_1[i] * wa_x[i];

        g_xz_0_yyyyyzzz_0[i] = g_z_0_yyyyyzzz_1[i] * wa_x[i];

        g_xz_0_yyyyzzzz_0[i] = g_z_0_yyyyzzzz_1[i] * wa_x[i];

        g_xz_0_yyyzzzzz_0[i] = g_z_0_yyyzzzzz_1[i] * wa_x[i];

        g_xz_0_yyzzzzzz_0[i] = g_z_0_yyzzzzzz_1[i] * wa_x[i];

        g_xz_0_yzzzzzzz_0[i] = g_z_0_yzzzzzzz_1[i] * wa_x[i];

        g_xz_0_zzzzzzzz_0[i] = g_z_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 135-180 components of targeted buffer : DSL

    auto g_yy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl + 135);

    auto g_yy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 136);

    auto g_yy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 137);

    auto g_yy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 138);

    auto g_yy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 139);

    auto g_yy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 140);

    auto g_yy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 141);

    auto g_yy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 142);

    auto g_yy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 143);

    auto g_yy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 144);

    auto g_yy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 145);

    auto g_yy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 146);

    auto g_yy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 147);

    auto g_yy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 148);

    auto g_yy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 149);

    auto g_yy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 150);

    auto g_yy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 151);

    auto g_yy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 152);

    auto g_yy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 153);

    auto g_yy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 154);

    auto g_yy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 155);

    auto g_yy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 156);

    auto g_yy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 157);

    auto g_yy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 158);

    auto g_yy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 159);

    auto g_yy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 160);

    auto g_yy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 161);

    auto g_yy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 162);

    auto g_yy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 163);

    auto g_yy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 164);

    auto g_yy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 165);

    auto g_yy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 166);

    auto g_yy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 167);

    auto g_yy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 168);

    auto g_yy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 169);

    auto g_yy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 170);

    auto g_yy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 171);

    auto g_yy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 172);

    auto g_yy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 173);

    auto g_yy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 174);

    auto g_yy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 175);

    auto g_yy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 176);

    auto g_yy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 177);

    auto g_yy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 178);

    auto g_yy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 179);

    #pragma omp simd aligned(g_0_0_xxxxxxxx_0, g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxy_0, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxz_0, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxyy_0, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyz_0, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxzz_0, g_0_0_xxxxxxzz_1, g_0_0_xxxxxyyy_0, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyz_0, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyzz_0, g_0_0_xxxxxyzz_1, g_0_0_xxxxxzzz_0, g_0_0_xxxxxzzz_1, g_0_0_xxxxyyyy_0, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyz_0, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyzz_0, g_0_0_xxxxyyzz_1, g_0_0_xxxxyzzz_0, g_0_0_xxxxyzzz_1, g_0_0_xxxxzzzz_0, g_0_0_xxxxzzzz_1, g_0_0_xxxyyyyy_0, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyz_0, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyzz_0, g_0_0_xxxyyyzz_1, g_0_0_xxxyyzzz_0, g_0_0_xxxyyzzz_1, g_0_0_xxxyzzzz_0, g_0_0_xxxyzzzz_1, g_0_0_xxxzzzzz_0, g_0_0_xxxzzzzz_1, g_0_0_xxyyyyyy_0, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyz_0, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyzz_0, g_0_0_xxyyyyzz_1, g_0_0_xxyyyzzz_0, g_0_0_xxyyyzzz_1, g_0_0_xxyyzzzz_0, g_0_0_xxyyzzzz_1, g_0_0_xxyzzzzz_0, g_0_0_xxyzzzzz_1, g_0_0_xxzzzzzz_0, g_0_0_xxzzzzzz_1, g_0_0_xyyyyyyy_0, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyz_0, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyzz_0, g_0_0_xyyyyyzz_1, g_0_0_xyyyyzzz_0, g_0_0_xyyyyzzz_1, g_0_0_xyyyzzzz_0, g_0_0_xyyyzzzz_1, g_0_0_xyyzzzzz_0, g_0_0_xyyzzzzz_1, g_0_0_xyzzzzzz_0, g_0_0_xyzzzzzz_1, g_0_0_xzzzzzzz_0, g_0_0_xzzzzzzz_1, g_0_0_yyyyyyyy_0, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyz_0, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyzz_0, g_0_0_yyyyyyzz_1, g_0_0_yyyyyzzz_0, g_0_0_yyyyyzzz_1, g_0_0_yyyyzzzz_0, g_0_0_yyyyzzzz_1, g_0_0_yyyzzzzz_0, g_0_0_yyyzzzzz_1, g_0_0_yyzzzzzz_0, g_0_0_yyzzzzzz_1, g_0_0_yzzzzzzz_0, g_0_0_yzzzzzzz_1, g_0_0_zzzzzzzz_0, g_0_0_zzzzzzzz_1, g_y_0_xxxxxxx_1, g_y_0_xxxxxxxx_1, g_y_0_xxxxxxxy_1, g_y_0_xxxxxxxz_1, g_y_0_xxxxxxy_1, g_y_0_xxxxxxyy_1, g_y_0_xxxxxxyz_1, g_y_0_xxxxxxz_1, g_y_0_xxxxxxzz_1, g_y_0_xxxxxyy_1, g_y_0_xxxxxyyy_1, g_y_0_xxxxxyyz_1, g_y_0_xxxxxyz_1, g_y_0_xxxxxyzz_1, g_y_0_xxxxxzz_1, g_y_0_xxxxxzzz_1, g_y_0_xxxxyyy_1, g_y_0_xxxxyyyy_1, g_y_0_xxxxyyyz_1, g_y_0_xxxxyyz_1, g_y_0_xxxxyyzz_1, g_y_0_xxxxyzz_1, g_y_0_xxxxyzzz_1, g_y_0_xxxxzzz_1, g_y_0_xxxxzzzz_1, g_y_0_xxxyyyy_1, g_y_0_xxxyyyyy_1, g_y_0_xxxyyyyz_1, g_y_0_xxxyyyz_1, g_y_0_xxxyyyzz_1, g_y_0_xxxyyzz_1, g_y_0_xxxyyzzz_1, g_y_0_xxxyzzz_1, g_y_0_xxxyzzzz_1, g_y_0_xxxzzzz_1, g_y_0_xxxzzzzz_1, g_y_0_xxyyyyy_1, g_y_0_xxyyyyyy_1, g_y_0_xxyyyyyz_1, g_y_0_xxyyyyz_1, g_y_0_xxyyyyzz_1, g_y_0_xxyyyzz_1, g_y_0_xxyyyzzz_1, g_y_0_xxyyzzz_1, g_y_0_xxyyzzzz_1, g_y_0_xxyzzzz_1, g_y_0_xxyzzzzz_1, g_y_0_xxzzzzz_1, g_y_0_xxzzzzzz_1, g_y_0_xyyyyyy_1, g_y_0_xyyyyyyy_1, g_y_0_xyyyyyyz_1, g_y_0_xyyyyyz_1, g_y_0_xyyyyyzz_1, g_y_0_xyyyyzz_1, g_y_0_xyyyyzzz_1, g_y_0_xyyyzzz_1, g_y_0_xyyyzzzz_1, g_y_0_xyyzzzz_1, g_y_0_xyyzzzzz_1, g_y_0_xyzzzzz_1, g_y_0_xyzzzzzz_1, g_y_0_xzzzzzz_1, g_y_0_xzzzzzzz_1, g_y_0_yyyyyyy_1, g_y_0_yyyyyyyy_1, g_y_0_yyyyyyyz_1, g_y_0_yyyyyyz_1, g_y_0_yyyyyyzz_1, g_y_0_yyyyyzz_1, g_y_0_yyyyyzzz_1, g_y_0_yyyyzzz_1, g_y_0_yyyyzzzz_1, g_y_0_yyyzzzz_1, g_y_0_yyyzzzzz_1, g_y_0_yyzzzzz_1, g_y_0_yyzzzzzz_1, g_y_0_yzzzzzz_1, g_y_0_yzzzzzzz_1, g_y_0_zzzzzzz_1, g_y_0_zzzzzzzz_1, g_yy_0_xxxxxxxx_0, g_yy_0_xxxxxxxy_0, g_yy_0_xxxxxxxz_0, g_yy_0_xxxxxxyy_0, g_yy_0_xxxxxxyz_0, g_yy_0_xxxxxxzz_0, g_yy_0_xxxxxyyy_0, g_yy_0_xxxxxyyz_0, g_yy_0_xxxxxyzz_0, g_yy_0_xxxxxzzz_0, g_yy_0_xxxxyyyy_0, g_yy_0_xxxxyyyz_0, g_yy_0_xxxxyyzz_0, g_yy_0_xxxxyzzz_0, g_yy_0_xxxxzzzz_0, g_yy_0_xxxyyyyy_0, g_yy_0_xxxyyyyz_0, g_yy_0_xxxyyyzz_0, g_yy_0_xxxyyzzz_0, g_yy_0_xxxyzzzz_0, g_yy_0_xxxzzzzz_0, g_yy_0_xxyyyyyy_0, g_yy_0_xxyyyyyz_0, g_yy_0_xxyyyyzz_0, g_yy_0_xxyyyzzz_0, g_yy_0_xxyyzzzz_0, g_yy_0_xxyzzzzz_0, g_yy_0_xxzzzzzz_0, g_yy_0_xyyyyyyy_0, g_yy_0_xyyyyyyz_0, g_yy_0_xyyyyyzz_0, g_yy_0_xyyyyzzz_0, g_yy_0_xyyyzzzz_0, g_yy_0_xyyzzzzz_0, g_yy_0_xyzzzzzz_0, g_yy_0_xzzzzzzz_0, g_yy_0_yyyyyyyy_0, g_yy_0_yyyyyyyz_0, g_yy_0_yyyyyyzz_0, g_yy_0_yyyyyzzz_0, g_yy_0_yyyyzzzz_0, g_yy_0_yyyzzzzz_0, g_yy_0_yyzzzzzz_0, g_yy_0_yzzzzzzz_0, g_yy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xxxxxxxx_0[i] = g_0_0_xxxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxxx_1[i] * fz_be_0 + g_y_0_xxxxxxxx_1[i] * wa_y[i];

        g_yy_0_xxxxxxxy_0[i] = g_0_0_xxxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxxy_1[i] * fz_be_0 + g_y_0_xxxxxxx_1[i] * fi_acd_0 + g_y_0_xxxxxxxy_1[i] * wa_y[i];

        g_yy_0_xxxxxxxz_0[i] = g_0_0_xxxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxxz_1[i] * fz_be_0 + g_y_0_xxxxxxxz_1[i] * wa_y[i];

        g_yy_0_xxxxxxyy_0[i] = g_0_0_xxxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxxyy_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxxxy_1[i] * fi_acd_0 + g_y_0_xxxxxxyy_1[i] * wa_y[i];

        g_yy_0_xxxxxxyz_0[i] = g_0_0_xxxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxxyz_1[i] * fz_be_0 + g_y_0_xxxxxxz_1[i] * fi_acd_0 + g_y_0_xxxxxxyz_1[i] * wa_y[i];

        g_yy_0_xxxxxxzz_0[i] = g_0_0_xxxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxxzz_1[i] * fz_be_0 + g_y_0_xxxxxxzz_1[i] * wa_y[i];

        g_yy_0_xxxxxyyy_0[i] = g_0_0_xxxxxyyy_0[i] * fbe_0 - g_0_0_xxxxxyyy_1[i] * fz_be_0 + 3.0 * g_y_0_xxxxxyy_1[i] * fi_acd_0 + g_y_0_xxxxxyyy_1[i] * wa_y[i];

        g_yy_0_xxxxxyyz_0[i] = g_0_0_xxxxxyyz_0[i] * fbe_0 - g_0_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxxyz_1[i] * fi_acd_0 + g_y_0_xxxxxyyz_1[i] * wa_y[i];

        g_yy_0_xxxxxyzz_0[i] = g_0_0_xxxxxyzz_0[i] * fbe_0 - g_0_0_xxxxxyzz_1[i] * fz_be_0 + g_y_0_xxxxxzz_1[i] * fi_acd_0 + g_y_0_xxxxxyzz_1[i] * wa_y[i];

        g_yy_0_xxxxxzzz_0[i] = g_0_0_xxxxxzzz_0[i] * fbe_0 - g_0_0_xxxxxzzz_1[i] * fz_be_0 + g_y_0_xxxxxzzz_1[i] * wa_y[i];

        g_yy_0_xxxxyyyy_0[i] = g_0_0_xxxxyyyy_0[i] * fbe_0 - g_0_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_y_0_xxxxyyy_1[i] * fi_acd_0 + g_y_0_xxxxyyyy_1[i] * wa_y[i];

        g_yy_0_xxxxyyyz_0[i] = g_0_0_xxxxyyyz_0[i] * fbe_0 - g_0_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_y_0_xxxxyyz_1[i] * fi_acd_0 + g_y_0_xxxxyyyz_1[i] * wa_y[i];

        g_yy_0_xxxxyyzz_0[i] = g_0_0_xxxxyyzz_0[i] * fbe_0 - g_0_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxyzz_1[i] * fi_acd_0 + g_y_0_xxxxyyzz_1[i] * wa_y[i];

        g_yy_0_xxxxyzzz_0[i] = g_0_0_xxxxyzzz_0[i] * fbe_0 - g_0_0_xxxxyzzz_1[i] * fz_be_0 + g_y_0_xxxxzzz_1[i] * fi_acd_0 + g_y_0_xxxxyzzz_1[i] * wa_y[i];

        g_yy_0_xxxxzzzz_0[i] = g_0_0_xxxxzzzz_0[i] * fbe_0 - g_0_0_xxxxzzzz_1[i] * fz_be_0 + g_y_0_xxxxzzzz_1[i] * wa_y[i];

        g_yy_0_xxxyyyyy_0[i] = g_0_0_xxxyyyyy_0[i] * fbe_0 - g_0_0_xxxyyyyy_1[i] * fz_be_0 + 5.0 * g_y_0_xxxyyyy_1[i] * fi_acd_0 + g_y_0_xxxyyyyy_1[i] * wa_y[i];

        g_yy_0_xxxyyyyz_0[i] = g_0_0_xxxyyyyz_0[i] * fbe_0 - g_0_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_y_0_xxxyyyz_1[i] * fi_acd_0 + g_y_0_xxxyyyyz_1[i] * wa_y[i];

        g_yy_0_xxxyyyzz_0[i] = g_0_0_xxxyyyzz_0[i] * fbe_0 - g_0_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_y_0_xxxyyzz_1[i] * fi_acd_0 + g_y_0_xxxyyyzz_1[i] * wa_y[i];

        g_yy_0_xxxyyzzz_0[i] = g_0_0_xxxyyzzz_0[i] * fbe_0 - g_0_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxyzzz_1[i] * fi_acd_0 + g_y_0_xxxyyzzz_1[i] * wa_y[i];

        g_yy_0_xxxyzzzz_0[i] = g_0_0_xxxyzzzz_0[i] * fbe_0 - g_0_0_xxxyzzzz_1[i] * fz_be_0 + g_y_0_xxxzzzz_1[i] * fi_acd_0 + g_y_0_xxxyzzzz_1[i] * wa_y[i];

        g_yy_0_xxxzzzzz_0[i] = g_0_0_xxxzzzzz_0[i] * fbe_0 - g_0_0_xxxzzzzz_1[i] * fz_be_0 + g_y_0_xxxzzzzz_1[i] * wa_y[i];

        g_yy_0_xxyyyyyy_0[i] = g_0_0_xxyyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyyy_1[i] * fz_be_0 + 6.0 * g_y_0_xxyyyyy_1[i] * fi_acd_0 + g_y_0_xxyyyyyy_1[i] * wa_y[i];

        g_yy_0_xxyyyyyz_0[i] = g_0_0_xxyyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_y_0_xxyyyyz_1[i] * fi_acd_0 + g_y_0_xxyyyyyz_1[i] * wa_y[i];

        g_yy_0_xxyyyyzz_0[i] = g_0_0_xxyyyyzz_0[i] * fbe_0 - g_0_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_y_0_xxyyyzz_1[i] * fi_acd_0 + g_y_0_xxyyyyzz_1[i] * wa_y[i];

        g_yy_0_xxyyyzzz_0[i] = g_0_0_xxyyyzzz_0[i] * fbe_0 - g_0_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_y_0_xxyyzzz_1[i] * fi_acd_0 + g_y_0_xxyyyzzz_1[i] * wa_y[i];

        g_yy_0_xxyyzzzz_0[i] = g_0_0_xxyyzzzz_0[i] * fbe_0 - g_0_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxyzzzz_1[i] * fi_acd_0 + g_y_0_xxyyzzzz_1[i] * wa_y[i];

        g_yy_0_xxyzzzzz_0[i] = g_0_0_xxyzzzzz_0[i] * fbe_0 - g_0_0_xxyzzzzz_1[i] * fz_be_0 + g_y_0_xxzzzzz_1[i] * fi_acd_0 + g_y_0_xxyzzzzz_1[i] * wa_y[i];

        g_yy_0_xxzzzzzz_0[i] = g_0_0_xxzzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzzz_1[i] * fz_be_0 + g_y_0_xxzzzzzz_1[i] * wa_y[i];

        g_yy_0_xyyyyyyy_0[i] = g_0_0_xyyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyyy_1[i] * fz_be_0 + 7.0 * g_y_0_xyyyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyyyy_1[i] * wa_y[i];

        g_yy_0_xyyyyyyz_0[i] = g_0_0_xyyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_y_0_xyyyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyyyz_1[i] * wa_y[i];

        g_yy_0_xyyyyyzz_0[i] = g_0_0_xyyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_y_0_xyyyyzz_1[i] * fi_acd_0 + g_y_0_xyyyyyzz_1[i] * wa_y[i];

        g_yy_0_xyyyyzzz_0[i] = g_0_0_xyyyyzzz_0[i] * fbe_0 - g_0_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_y_0_xyyyzzz_1[i] * fi_acd_0 + g_y_0_xyyyyzzz_1[i] * wa_y[i];

        g_yy_0_xyyyzzzz_0[i] = g_0_0_xyyyzzzz_0[i] * fbe_0 - g_0_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_y_0_xyyzzzz_1[i] * fi_acd_0 + g_y_0_xyyyzzzz_1[i] * wa_y[i];

        g_yy_0_xyyzzzzz_0[i] = g_0_0_xyyzzzzz_0[i] * fbe_0 - g_0_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xyzzzzz_1[i] * fi_acd_0 + g_y_0_xyyzzzzz_1[i] * wa_y[i];

        g_yy_0_xyzzzzzz_0[i] = g_0_0_xyzzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzzz_1[i] * fz_be_0 + g_y_0_xzzzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzzzz_1[i] * wa_y[i];

        g_yy_0_xzzzzzzz_0[i] = g_0_0_xzzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzzz_1[i] * fz_be_0 + g_y_0_xzzzzzzz_1[i] * wa_y[i];

        g_yy_0_yyyyyyyy_0[i] = g_0_0_yyyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyyy_1[i] * fz_be_0 + 8.0 * g_y_0_yyyyyyy_1[i] * fi_acd_0 + g_y_0_yyyyyyyy_1[i] * wa_y[i];

        g_yy_0_yyyyyyyz_0[i] = g_0_0_yyyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_y_0_yyyyyyz_1[i] * fi_acd_0 + g_y_0_yyyyyyyz_1[i] * wa_y[i];

        g_yy_0_yyyyyyzz_0[i] = g_0_0_yyyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_y_0_yyyyyzz_1[i] * fi_acd_0 + g_y_0_yyyyyyzz_1[i] * wa_y[i];

        g_yy_0_yyyyyzzz_0[i] = g_0_0_yyyyyzzz_0[i] * fbe_0 - g_0_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_y_0_yyyyzzz_1[i] * fi_acd_0 + g_y_0_yyyyyzzz_1[i] * wa_y[i];

        g_yy_0_yyyyzzzz_0[i] = g_0_0_yyyyzzzz_0[i] * fbe_0 - g_0_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_y_0_yyyzzzz_1[i] * fi_acd_0 + g_y_0_yyyyzzzz_1[i] * wa_y[i];

        g_yy_0_yyyzzzzz_0[i] = g_0_0_yyyzzzzz_0[i] * fbe_0 - g_0_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_y_0_yyzzzzz_1[i] * fi_acd_0 + g_y_0_yyyzzzzz_1[i] * wa_y[i];

        g_yy_0_yyzzzzzz_0[i] = g_0_0_yyzzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_yzzzzzz_1[i] * fi_acd_0 + g_y_0_yyzzzzzz_1[i] * wa_y[i];

        g_yy_0_yzzzzzzz_0[i] = g_0_0_yzzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzzz_1[i] * fz_be_0 + g_y_0_zzzzzzz_1[i] * fi_acd_0 + g_y_0_yzzzzzzz_1[i] * wa_y[i];

        g_yy_0_zzzzzzzz_0[i] = g_0_0_zzzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzzz_1[i] * fz_be_0 + g_y_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 180-225 components of targeted buffer : DSL

    auto g_yz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl + 180);

    auto g_yz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 181);

    auto g_yz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 182);

    auto g_yz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 183);

    auto g_yz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 184);

    auto g_yz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 185);

    auto g_yz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 186);

    auto g_yz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 187);

    auto g_yz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 188);

    auto g_yz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 189);

    auto g_yz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 190);

    auto g_yz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 191);

    auto g_yz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 192);

    auto g_yz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 193);

    auto g_yz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 194);

    auto g_yz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 195);

    auto g_yz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 196);

    auto g_yz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 197);

    auto g_yz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 198);

    auto g_yz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 199);

    auto g_yz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 200);

    auto g_yz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 201);

    auto g_yz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 202);

    auto g_yz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 203);

    auto g_yz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 204);

    auto g_yz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 205);

    auto g_yz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 206);

    auto g_yz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 207);

    auto g_yz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 208);

    auto g_yz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 209);

    auto g_yz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 210);

    auto g_yz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 211);

    auto g_yz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 212);

    auto g_yz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 213);

    auto g_yz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 214);

    auto g_yz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 215);

    auto g_yz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 216);

    auto g_yz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 217);

    auto g_yz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 218);

    auto g_yz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 219);

    auto g_yz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 220);

    auto g_yz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 221);

    auto g_yz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 222);

    auto g_yz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 223);

    auto g_yz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 224);

    #pragma omp simd aligned(g_y_0_xxxxxxxy_1, g_y_0_xxxxxxyy_1, g_y_0_xxxxxyyy_1, g_y_0_xxxxyyyy_1, g_y_0_xxxyyyyy_1, g_y_0_xxyyyyyy_1, g_y_0_xyyyyyyy_1, g_y_0_yyyyyyyy_1, g_yz_0_xxxxxxxx_0, g_yz_0_xxxxxxxy_0, g_yz_0_xxxxxxxz_0, g_yz_0_xxxxxxyy_0, g_yz_0_xxxxxxyz_0, g_yz_0_xxxxxxzz_0, g_yz_0_xxxxxyyy_0, g_yz_0_xxxxxyyz_0, g_yz_0_xxxxxyzz_0, g_yz_0_xxxxxzzz_0, g_yz_0_xxxxyyyy_0, g_yz_0_xxxxyyyz_0, g_yz_0_xxxxyyzz_0, g_yz_0_xxxxyzzz_0, g_yz_0_xxxxzzzz_0, g_yz_0_xxxyyyyy_0, g_yz_0_xxxyyyyz_0, g_yz_0_xxxyyyzz_0, g_yz_0_xxxyyzzz_0, g_yz_0_xxxyzzzz_0, g_yz_0_xxxzzzzz_0, g_yz_0_xxyyyyyy_0, g_yz_0_xxyyyyyz_0, g_yz_0_xxyyyyzz_0, g_yz_0_xxyyyzzz_0, g_yz_0_xxyyzzzz_0, g_yz_0_xxyzzzzz_0, g_yz_0_xxzzzzzz_0, g_yz_0_xyyyyyyy_0, g_yz_0_xyyyyyyz_0, g_yz_0_xyyyyyzz_0, g_yz_0_xyyyyzzz_0, g_yz_0_xyyyzzzz_0, g_yz_0_xyyzzzzz_0, g_yz_0_xyzzzzzz_0, g_yz_0_xzzzzzzz_0, g_yz_0_yyyyyyyy_0, g_yz_0_yyyyyyyz_0, g_yz_0_yyyyyyzz_0, g_yz_0_yyyyyzzz_0, g_yz_0_yyyyzzzz_0, g_yz_0_yyyzzzzz_0, g_yz_0_yyzzzzzz_0, g_yz_0_yzzzzzzz_0, g_yz_0_zzzzzzzz_0, g_z_0_xxxxxxxx_1, g_z_0_xxxxxxxz_1, g_z_0_xxxxxxyz_1, g_z_0_xxxxxxz_1, g_z_0_xxxxxxzz_1, g_z_0_xxxxxyyz_1, g_z_0_xxxxxyz_1, g_z_0_xxxxxyzz_1, g_z_0_xxxxxzz_1, g_z_0_xxxxxzzz_1, g_z_0_xxxxyyyz_1, g_z_0_xxxxyyz_1, g_z_0_xxxxyyzz_1, g_z_0_xxxxyzz_1, g_z_0_xxxxyzzz_1, g_z_0_xxxxzzz_1, g_z_0_xxxxzzzz_1, g_z_0_xxxyyyyz_1, g_z_0_xxxyyyz_1, g_z_0_xxxyyyzz_1, g_z_0_xxxyyzz_1, g_z_0_xxxyyzzz_1, g_z_0_xxxyzzz_1, g_z_0_xxxyzzzz_1, g_z_0_xxxzzzz_1, g_z_0_xxxzzzzz_1, g_z_0_xxyyyyyz_1, g_z_0_xxyyyyz_1, g_z_0_xxyyyyzz_1, g_z_0_xxyyyzz_1, g_z_0_xxyyyzzz_1, g_z_0_xxyyzzz_1, g_z_0_xxyyzzzz_1, g_z_0_xxyzzzz_1, g_z_0_xxyzzzzz_1, g_z_0_xxzzzzz_1, g_z_0_xxzzzzzz_1, g_z_0_xyyyyyyz_1, g_z_0_xyyyyyz_1, g_z_0_xyyyyyzz_1, g_z_0_xyyyyzz_1, g_z_0_xyyyyzzz_1, g_z_0_xyyyzzz_1, g_z_0_xyyyzzzz_1, g_z_0_xyyzzzz_1, g_z_0_xyyzzzzz_1, g_z_0_xyzzzzz_1, g_z_0_xyzzzzzz_1, g_z_0_xzzzzzz_1, g_z_0_xzzzzzzz_1, g_z_0_yyyyyyyz_1, g_z_0_yyyyyyz_1, g_z_0_yyyyyyzz_1, g_z_0_yyyyyzz_1, g_z_0_yyyyyzzz_1, g_z_0_yyyyzzz_1, g_z_0_yyyyzzzz_1, g_z_0_yyyzzzz_1, g_z_0_yyyzzzzz_1, g_z_0_yyzzzzz_1, g_z_0_yyzzzzzz_1, g_z_0_yzzzzzz_1, g_z_0_yzzzzzzz_1, g_z_0_zzzzzzz_1, g_z_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xxxxxxxx_0[i] = g_z_0_xxxxxxxx_1[i] * wa_y[i];

        g_yz_0_xxxxxxxy_0[i] = g_y_0_xxxxxxxy_1[i] * wa_z[i];

        g_yz_0_xxxxxxxz_0[i] = g_z_0_xxxxxxxz_1[i] * wa_y[i];

        g_yz_0_xxxxxxyy_0[i] = g_y_0_xxxxxxyy_1[i] * wa_z[i];

        g_yz_0_xxxxxxyz_0[i] = g_z_0_xxxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxxyz_1[i] * wa_y[i];

        g_yz_0_xxxxxxzz_0[i] = g_z_0_xxxxxxzz_1[i] * wa_y[i];

        g_yz_0_xxxxxyyy_0[i] = g_y_0_xxxxxyyy_1[i] * wa_z[i];

        g_yz_0_xxxxxyyz_0[i] = 2.0 * g_z_0_xxxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxxyyz_1[i] * wa_y[i];

        g_yz_0_xxxxxyzz_0[i] = g_z_0_xxxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxxyzz_1[i] * wa_y[i];

        g_yz_0_xxxxxzzz_0[i] = g_z_0_xxxxxzzz_1[i] * wa_y[i];

        g_yz_0_xxxxyyyy_0[i] = g_y_0_xxxxyyyy_1[i] * wa_z[i];

        g_yz_0_xxxxyyyz_0[i] = 3.0 * g_z_0_xxxxyyz_1[i] * fi_acd_0 + g_z_0_xxxxyyyz_1[i] * wa_y[i];

        g_yz_0_xxxxyyzz_0[i] = 2.0 * g_z_0_xxxxyzz_1[i] * fi_acd_0 + g_z_0_xxxxyyzz_1[i] * wa_y[i];

        g_yz_0_xxxxyzzz_0[i] = g_z_0_xxxxzzz_1[i] * fi_acd_0 + g_z_0_xxxxyzzz_1[i] * wa_y[i];

        g_yz_0_xxxxzzzz_0[i] = g_z_0_xxxxzzzz_1[i] * wa_y[i];

        g_yz_0_xxxyyyyy_0[i] = g_y_0_xxxyyyyy_1[i] * wa_z[i];

        g_yz_0_xxxyyyyz_0[i] = 4.0 * g_z_0_xxxyyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyyz_1[i] * wa_y[i];

        g_yz_0_xxxyyyzz_0[i] = 3.0 * g_z_0_xxxyyzz_1[i] * fi_acd_0 + g_z_0_xxxyyyzz_1[i] * wa_y[i];

        g_yz_0_xxxyyzzz_0[i] = 2.0 * g_z_0_xxxyzzz_1[i] * fi_acd_0 + g_z_0_xxxyyzzz_1[i] * wa_y[i];

        g_yz_0_xxxyzzzz_0[i] = g_z_0_xxxzzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzzz_1[i] * wa_y[i];

        g_yz_0_xxxzzzzz_0[i] = g_z_0_xxxzzzzz_1[i] * wa_y[i];

        g_yz_0_xxyyyyyy_0[i] = g_y_0_xxyyyyyy_1[i] * wa_z[i];

        g_yz_0_xxyyyyyz_0[i] = 5.0 * g_z_0_xxyyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyyz_1[i] * wa_y[i];

        g_yz_0_xxyyyyzz_0[i] = 4.0 * g_z_0_xxyyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyyzz_1[i] * wa_y[i];

        g_yz_0_xxyyyzzz_0[i] = 3.0 * g_z_0_xxyyzzz_1[i] * fi_acd_0 + g_z_0_xxyyyzzz_1[i] * wa_y[i];

        g_yz_0_xxyyzzzz_0[i] = 2.0 * g_z_0_xxyzzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzzz_1[i] * wa_y[i];

        g_yz_0_xxyzzzzz_0[i] = g_z_0_xxzzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzzz_1[i] * wa_y[i];

        g_yz_0_xxzzzzzz_0[i] = g_z_0_xxzzzzzz_1[i] * wa_y[i];

        g_yz_0_xyyyyyyy_0[i] = g_y_0_xyyyyyyy_1[i] * wa_z[i];

        g_yz_0_xyyyyyyz_0[i] = 6.0 * g_z_0_xyyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyyz_1[i] * wa_y[i];

        g_yz_0_xyyyyyzz_0[i] = 5.0 * g_z_0_xyyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyyzz_1[i] * wa_y[i];

        g_yz_0_xyyyyzzz_0[i] = 4.0 * g_z_0_xyyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyyzzz_1[i] * wa_y[i];

        g_yz_0_xyyyzzzz_0[i] = 3.0 * g_z_0_xyyzzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzzz_1[i] * wa_y[i];

        g_yz_0_xyyzzzzz_0[i] = 2.0 * g_z_0_xyzzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzzz_1[i] * wa_y[i];

        g_yz_0_xyzzzzzz_0[i] = g_z_0_xzzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzzz_1[i] * wa_y[i];

        g_yz_0_xzzzzzzz_0[i] = g_z_0_xzzzzzzz_1[i] * wa_y[i];

        g_yz_0_yyyyyyyy_0[i] = g_y_0_yyyyyyyy_1[i] * wa_z[i];

        g_yz_0_yyyyyyyz_0[i] = 7.0 * g_z_0_yyyyyyz_1[i] * fi_acd_0 + g_z_0_yyyyyyyz_1[i] * wa_y[i];

        g_yz_0_yyyyyyzz_0[i] = 6.0 * g_z_0_yyyyyzz_1[i] * fi_acd_0 + g_z_0_yyyyyyzz_1[i] * wa_y[i];

        g_yz_0_yyyyyzzz_0[i] = 5.0 * g_z_0_yyyyzzz_1[i] * fi_acd_0 + g_z_0_yyyyyzzz_1[i] * wa_y[i];

        g_yz_0_yyyyzzzz_0[i] = 4.0 * g_z_0_yyyzzzz_1[i] * fi_acd_0 + g_z_0_yyyyzzzz_1[i] * wa_y[i];

        g_yz_0_yyyzzzzz_0[i] = 3.0 * g_z_0_yyzzzzz_1[i] * fi_acd_0 + g_z_0_yyyzzzzz_1[i] * wa_y[i];

        g_yz_0_yyzzzzzz_0[i] = 2.0 * g_z_0_yzzzzzz_1[i] * fi_acd_0 + g_z_0_yyzzzzzz_1[i] * wa_y[i];

        g_yz_0_yzzzzzzz_0[i] = g_z_0_zzzzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzzzz_1[i] * wa_y[i];

        g_yz_0_zzzzzzzz_0[i] = g_z_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 225-270 components of targeted buffer : DSL

    auto g_zz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl + 225);

    auto g_zz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 226);

    auto g_zz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 227);

    auto g_zz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 228);

    auto g_zz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 229);

    auto g_zz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 230);

    auto g_zz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 231);

    auto g_zz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 232);

    auto g_zz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 233);

    auto g_zz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 234);

    auto g_zz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 235);

    auto g_zz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 236);

    auto g_zz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 237);

    auto g_zz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 238);

    auto g_zz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 239);

    auto g_zz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 240);

    auto g_zz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 241);

    auto g_zz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 242);

    auto g_zz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 243);

    auto g_zz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 244);

    auto g_zz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 245);

    auto g_zz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 246);

    auto g_zz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 247);

    auto g_zz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 248);

    auto g_zz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 249);

    auto g_zz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 250);

    auto g_zz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 251);

    auto g_zz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 252);

    auto g_zz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 253);

    auto g_zz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 254);

    auto g_zz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 255);

    auto g_zz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 256);

    auto g_zz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 257);

    auto g_zz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 258);

    auto g_zz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 259);

    auto g_zz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 260);

    auto g_zz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 261);

    auto g_zz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 262);

    auto g_zz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 263);

    auto g_zz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 264);

    auto g_zz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 265);

    auto g_zz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 266);

    auto g_zz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 267);

    auto g_zz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 268);

    auto g_zz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 269);

    #pragma omp simd aligned(g_0_0_xxxxxxxx_0, g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxy_0, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxz_0, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxyy_0, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyz_0, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxzz_0, g_0_0_xxxxxxzz_1, g_0_0_xxxxxyyy_0, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyz_0, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyzz_0, g_0_0_xxxxxyzz_1, g_0_0_xxxxxzzz_0, g_0_0_xxxxxzzz_1, g_0_0_xxxxyyyy_0, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyz_0, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyzz_0, g_0_0_xxxxyyzz_1, g_0_0_xxxxyzzz_0, g_0_0_xxxxyzzz_1, g_0_0_xxxxzzzz_0, g_0_0_xxxxzzzz_1, g_0_0_xxxyyyyy_0, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyz_0, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyzz_0, g_0_0_xxxyyyzz_1, g_0_0_xxxyyzzz_0, g_0_0_xxxyyzzz_1, g_0_0_xxxyzzzz_0, g_0_0_xxxyzzzz_1, g_0_0_xxxzzzzz_0, g_0_0_xxxzzzzz_1, g_0_0_xxyyyyyy_0, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyz_0, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyzz_0, g_0_0_xxyyyyzz_1, g_0_0_xxyyyzzz_0, g_0_0_xxyyyzzz_1, g_0_0_xxyyzzzz_0, g_0_0_xxyyzzzz_1, g_0_0_xxyzzzzz_0, g_0_0_xxyzzzzz_1, g_0_0_xxzzzzzz_0, g_0_0_xxzzzzzz_1, g_0_0_xyyyyyyy_0, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyz_0, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyzz_0, g_0_0_xyyyyyzz_1, g_0_0_xyyyyzzz_0, g_0_0_xyyyyzzz_1, g_0_0_xyyyzzzz_0, g_0_0_xyyyzzzz_1, g_0_0_xyyzzzzz_0, g_0_0_xyyzzzzz_1, g_0_0_xyzzzzzz_0, g_0_0_xyzzzzzz_1, g_0_0_xzzzzzzz_0, g_0_0_xzzzzzzz_1, g_0_0_yyyyyyyy_0, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyz_0, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyzz_0, g_0_0_yyyyyyzz_1, g_0_0_yyyyyzzz_0, g_0_0_yyyyyzzz_1, g_0_0_yyyyzzzz_0, g_0_0_yyyyzzzz_1, g_0_0_yyyzzzzz_0, g_0_0_yyyzzzzz_1, g_0_0_yyzzzzzz_0, g_0_0_yyzzzzzz_1, g_0_0_yzzzzzzz_0, g_0_0_yzzzzzzz_1, g_0_0_zzzzzzzz_0, g_0_0_zzzzzzzz_1, g_z_0_xxxxxxx_1, g_z_0_xxxxxxxx_1, g_z_0_xxxxxxxy_1, g_z_0_xxxxxxxz_1, g_z_0_xxxxxxy_1, g_z_0_xxxxxxyy_1, g_z_0_xxxxxxyz_1, g_z_0_xxxxxxz_1, g_z_0_xxxxxxzz_1, g_z_0_xxxxxyy_1, g_z_0_xxxxxyyy_1, g_z_0_xxxxxyyz_1, g_z_0_xxxxxyz_1, g_z_0_xxxxxyzz_1, g_z_0_xxxxxzz_1, g_z_0_xxxxxzzz_1, g_z_0_xxxxyyy_1, g_z_0_xxxxyyyy_1, g_z_0_xxxxyyyz_1, g_z_0_xxxxyyz_1, g_z_0_xxxxyyzz_1, g_z_0_xxxxyzz_1, g_z_0_xxxxyzzz_1, g_z_0_xxxxzzz_1, g_z_0_xxxxzzzz_1, g_z_0_xxxyyyy_1, g_z_0_xxxyyyyy_1, g_z_0_xxxyyyyz_1, g_z_0_xxxyyyz_1, g_z_0_xxxyyyzz_1, g_z_0_xxxyyzz_1, g_z_0_xxxyyzzz_1, g_z_0_xxxyzzz_1, g_z_0_xxxyzzzz_1, g_z_0_xxxzzzz_1, g_z_0_xxxzzzzz_1, g_z_0_xxyyyyy_1, g_z_0_xxyyyyyy_1, g_z_0_xxyyyyyz_1, g_z_0_xxyyyyz_1, g_z_0_xxyyyyzz_1, g_z_0_xxyyyzz_1, g_z_0_xxyyyzzz_1, g_z_0_xxyyzzz_1, g_z_0_xxyyzzzz_1, g_z_0_xxyzzzz_1, g_z_0_xxyzzzzz_1, g_z_0_xxzzzzz_1, g_z_0_xxzzzzzz_1, g_z_0_xyyyyyy_1, g_z_0_xyyyyyyy_1, g_z_0_xyyyyyyz_1, g_z_0_xyyyyyz_1, g_z_0_xyyyyyzz_1, g_z_0_xyyyyzz_1, g_z_0_xyyyyzzz_1, g_z_0_xyyyzzz_1, g_z_0_xyyyzzzz_1, g_z_0_xyyzzzz_1, g_z_0_xyyzzzzz_1, g_z_0_xyzzzzz_1, g_z_0_xyzzzzzz_1, g_z_0_xzzzzzz_1, g_z_0_xzzzzzzz_1, g_z_0_yyyyyyy_1, g_z_0_yyyyyyyy_1, g_z_0_yyyyyyyz_1, g_z_0_yyyyyyz_1, g_z_0_yyyyyyzz_1, g_z_0_yyyyyzz_1, g_z_0_yyyyyzzz_1, g_z_0_yyyyzzz_1, g_z_0_yyyyzzzz_1, g_z_0_yyyzzzz_1, g_z_0_yyyzzzzz_1, g_z_0_yyzzzzz_1, g_z_0_yyzzzzzz_1, g_z_0_yzzzzzz_1, g_z_0_yzzzzzzz_1, g_z_0_zzzzzzz_1, g_z_0_zzzzzzzz_1, g_zz_0_xxxxxxxx_0, g_zz_0_xxxxxxxy_0, g_zz_0_xxxxxxxz_0, g_zz_0_xxxxxxyy_0, g_zz_0_xxxxxxyz_0, g_zz_0_xxxxxxzz_0, g_zz_0_xxxxxyyy_0, g_zz_0_xxxxxyyz_0, g_zz_0_xxxxxyzz_0, g_zz_0_xxxxxzzz_0, g_zz_0_xxxxyyyy_0, g_zz_0_xxxxyyyz_0, g_zz_0_xxxxyyzz_0, g_zz_0_xxxxyzzz_0, g_zz_0_xxxxzzzz_0, g_zz_0_xxxyyyyy_0, g_zz_0_xxxyyyyz_0, g_zz_0_xxxyyyzz_0, g_zz_0_xxxyyzzz_0, g_zz_0_xxxyzzzz_0, g_zz_0_xxxzzzzz_0, g_zz_0_xxyyyyyy_0, g_zz_0_xxyyyyyz_0, g_zz_0_xxyyyyzz_0, g_zz_0_xxyyyzzz_0, g_zz_0_xxyyzzzz_0, g_zz_0_xxyzzzzz_0, g_zz_0_xxzzzzzz_0, g_zz_0_xyyyyyyy_0, g_zz_0_xyyyyyyz_0, g_zz_0_xyyyyyzz_0, g_zz_0_xyyyyzzz_0, g_zz_0_xyyyzzzz_0, g_zz_0_xyyzzzzz_0, g_zz_0_xyzzzzzz_0, g_zz_0_xzzzzzzz_0, g_zz_0_yyyyyyyy_0, g_zz_0_yyyyyyyz_0, g_zz_0_yyyyyyzz_0, g_zz_0_yyyyyzzz_0, g_zz_0_yyyyzzzz_0, g_zz_0_yyyzzzzz_0, g_zz_0_yyzzzzzz_0, g_zz_0_yzzzzzzz_0, g_zz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xxxxxxxx_0[i] = g_0_0_xxxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxxx_1[i] * fz_be_0 + g_z_0_xxxxxxxx_1[i] * wa_z[i];

        g_zz_0_xxxxxxxy_0[i] = g_0_0_xxxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxxy_1[i] * fz_be_0 + g_z_0_xxxxxxxy_1[i] * wa_z[i];

        g_zz_0_xxxxxxxz_0[i] = g_0_0_xxxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxxz_1[i] * fz_be_0 + g_z_0_xxxxxxx_1[i] * fi_acd_0 + g_z_0_xxxxxxxz_1[i] * wa_z[i];

        g_zz_0_xxxxxxyy_0[i] = g_0_0_xxxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxxyy_1[i] * fz_be_0 + g_z_0_xxxxxxyy_1[i] * wa_z[i];

        g_zz_0_xxxxxxyz_0[i] = g_0_0_xxxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxxyz_1[i] * fz_be_0 + g_z_0_xxxxxxy_1[i] * fi_acd_0 + g_z_0_xxxxxxyz_1[i] * wa_z[i];

        g_zz_0_xxxxxxzz_0[i] = g_0_0_xxxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxxzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxxzz_1[i] * wa_z[i];

        g_zz_0_xxxxxyyy_0[i] = g_0_0_xxxxxyyy_0[i] * fbe_0 - g_0_0_xxxxxyyy_1[i] * fz_be_0 + g_z_0_xxxxxyyy_1[i] * wa_z[i];

        g_zz_0_xxxxxyyz_0[i] = g_0_0_xxxxxyyz_0[i] * fbe_0 - g_0_0_xxxxxyyz_1[i] * fz_be_0 + g_z_0_xxxxxyy_1[i] * fi_acd_0 + g_z_0_xxxxxyyz_1[i] * wa_z[i];

        g_zz_0_xxxxxyzz_0[i] = g_0_0_xxxxxyzz_0[i] * fbe_0 - g_0_0_xxxxxyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxxyzz_1[i] * wa_z[i];

        g_zz_0_xxxxxzzz_0[i] = g_0_0_xxxxxzzz_0[i] * fbe_0 - g_0_0_xxxxxzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxxzzz_1[i] * wa_z[i];

        g_zz_0_xxxxyyyy_0[i] = g_0_0_xxxxyyyy_0[i] * fbe_0 - g_0_0_xxxxyyyy_1[i] * fz_be_0 + g_z_0_xxxxyyyy_1[i] * wa_z[i];

        g_zz_0_xxxxyyyz_0[i] = g_0_0_xxxxyyyz_0[i] * fbe_0 - g_0_0_xxxxyyyz_1[i] * fz_be_0 + g_z_0_xxxxyyy_1[i] * fi_acd_0 + g_z_0_xxxxyyyz_1[i] * wa_z[i];

        g_zz_0_xxxxyyzz_0[i] = g_0_0_xxxxyyzz_0[i] * fbe_0 - g_0_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxyyz_1[i] * fi_acd_0 + g_z_0_xxxxyyzz_1[i] * wa_z[i];

        g_zz_0_xxxxyzzz_0[i] = g_0_0_xxxxyzzz_0[i] * fbe_0 - g_0_0_xxxxyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxxyzz_1[i] * fi_acd_0 + g_z_0_xxxxyzzz_1[i] * wa_z[i];

        g_zz_0_xxxxzzzz_0[i] = g_0_0_xxxxzzzz_0[i] * fbe_0 - g_0_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxxxzzz_1[i] * fi_acd_0 + g_z_0_xxxxzzzz_1[i] * wa_z[i];

        g_zz_0_xxxyyyyy_0[i] = g_0_0_xxxyyyyy_0[i] * fbe_0 - g_0_0_xxxyyyyy_1[i] * fz_be_0 + g_z_0_xxxyyyyy_1[i] * wa_z[i];

        g_zz_0_xxxyyyyz_0[i] = g_0_0_xxxyyyyz_0[i] * fbe_0 - g_0_0_xxxyyyyz_1[i] * fz_be_0 + g_z_0_xxxyyyy_1[i] * fi_acd_0 + g_z_0_xxxyyyyz_1[i] * wa_z[i];

        g_zz_0_xxxyyyzz_0[i] = g_0_0_xxxyyyzz_0[i] * fbe_0 - g_0_0_xxxyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxyyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyzz_1[i] * wa_z[i];

        g_zz_0_xxxyyzzz_0[i] = g_0_0_xxxyyzzz_0[i] * fbe_0 - g_0_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxyyzz_1[i] * fi_acd_0 + g_z_0_xxxyyzzz_1[i] * wa_z[i];

        g_zz_0_xxxyzzzz_0[i] = g_0_0_xxxyzzzz_0[i] * fbe_0 - g_0_0_xxxyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxxyzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzzz_1[i] * wa_z[i];

        g_zz_0_xxxzzzzz_0[i] = g_0_0_xxxzzzzz_0[i] * fbe_0 - g_0_0_xxxzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xxxzzzz_1[i] * fi_acd_0 + g_z_0_xxxzzzzz_1[i] * wa_z[i];

        g_zz_0_xxyyyyyy_0[i] = g_0_0_xxyyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyyy_1[i] * fz_be_0 + g_z_0_xxyyyyyy_1[i] * wa_z[i];

        g_zz_0_xxyyyyyz_0[i] = g_0_0_xxyyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyyz_1[i] * fz_be_0 + g_z_0_xxyyyyy_1[i] * fi_acd_0 + g_z_0_xxyyyyyz_1[i] * wa_z[i];

        g_zz_0_xxyyyyzz_0[i] = g_0_0_xxyyyyzz_0[i] * fbe_0 - g_0_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxyyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyzz_1[i] * wa_z[i];

        g_zz_0_xxyyyzzz_0[i] = g_0_0_xxyyyzzz_0[i] * fbe_0 - g_0_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxyyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyzzz_1[i] * wa_z[i];

        g_zz_0_xxyyzzzz_0[i] = g_0_0_xxyyzzzz_0[i] * fbe_0 - g_0_0_xxyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxyyzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzzz_1[i] * wa_z[i];

        g_zz_0_xxyzzzzz_0[i] = g_0_0_xxyzzzzz_0[i] * fbe_0 - g_0_0_xxyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xxyzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzzz_1[i] * wa_z[i];

        g_zz_0_xxzzzzzz_0[i] = g_0_0_xxzzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_xxzzzzz_1[i] * fi_acd_0 + g_z_0_xxzzzzzz_1[i] * wa_z[i];

        g_zz_0_xyyyyyyy_0[i] = g_0_0_xyyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyyy_1[i] * fz_be_0 + g_z_0_xyyyyyyy_1[i] * wa_z[i];

        g_zz_0_xyyyyyyz_0[i] = g_0_0_xyyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyyz_1[i] * fz_be_0 + g_z_0_xyyyyyy_1[i] * fi_acd_0 + g_z_0_xyyyyyyz_1[i] * wa_z[i];

        g_zz_0_xyyyyyzz_0[i] = g_0_0_xyyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xyyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyzz_1[i] * wa_z[i];

        g_zz_0_xyyyyzzz_0[i] = g_0_0_xyyyyzzz_0[i] * fbe_0 - g_0_0_xyyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xyyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyzzz_1[i] * wa_z[i];

        g_zz_0_xyyyzzzz_0[i] = g_0_0_xyyyzzzz_0[i] * fbe_0 - g_0_0_xyyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xyyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzzz_1[i] * wa_z[i];

        g_zz_0_xyyzzzzz_0[i] = g_0_0_xyyzzzzz_0[i] * fbe_0 - g_0_0_xyyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xyyzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzzz_1[i] * wa_z[i];

        g_zz_0_xyzzzzzz_0[i] = g_0_0_xyzzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_xyzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzzz_1[i] * wa_z[i];

        g_zz_0_xzzzzzzz_0[i] = g_0_0_xzzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzzz_1[i] * fz_be_0 + 7.0 * g_z_0_xzzzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzzzz_1[i] * wa_z[i];

        g_zz_0_yyyyyyyy_0[i] = g_0_0_yyyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyyy_1[i] * fz_be_0 + g_z_0_yyyyyyyy_1[i] * wa_z[i];

        g_zz_0_yyyyyyyz_0[i] = g_0_0_yyyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyyz_1[i] * fz_be_0 + g_z_0_yyyyyyy_1[i] * fi_acd_0 + g_z_0_yyyyyyyz_1[i] * wa_z[i];

        g_zz_0_yyyyyyzz_0[i] = g_0_0_yyyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_yyyyyyz_1[i] * fi_acd_0 + g_z_0_yyyyyyzz_1[i] * wa_z[i];

        g_zz_0_yyyyyzzz_0[i] = g_0_0_yyyyyzzz_0[i] * fbe_0 - g_0_0_yyyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_yyyyyzz_1[i] * fi_acd_0 + g_z_0_yyyyyzzz_1[i] * wa_z[i];

        g_zz_0_yyyyzzzz_0[i] = g_0_0_yyyyzzzz_0[i] * fbe_0 - g_0_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_yyyyzzz_1[i] * fi_acd_0 + g_z_0_yyyyzzzz_1[i] * wa_z[i];

        g_zz_0_yyyzzzzz_0[i] = g_0_0_yyyzzzzz_0[i] * fbe_0 - g_0_0_yyyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_yyyzzzz_1[i] * fi_acd_0 + g_z_0_yyyzzzzz_1[i] * wa_z[i];

        g_zz_0_yyzzzzzz_0[i] = g_0_0_yyzzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_yyzzzzz_1[i] * fi_acd_0 + g_z_0_yyzzzzzz_1[i] * wa_z[i];

        g_zz_0_yzzzzzzz_0[i] = g_0_0_yzzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzzz_1[i] * fz_be_0 + 7.0 * g_z_0_yzzzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzzzz_1[i] * wa_z[i];

        g_zz_0_zzzzzzzz_0[i] = g_0_0_zzzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzzz_1[i] * fz_be_0 + 8.0 * g_z_0_zzzzzzz_1[i] * fi_acd_0 + g_z_0_zzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

