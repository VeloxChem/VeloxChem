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

#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsk,
                                 size_t idx_eri_0_ssk,
                                 size_t idx_eri_1_ssk,
                                 size_t idx_eri_1_psi,
                                 size_t idx_eri_1_psk,
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

    /// Set up components of auxilary buffer : SSK

    auto g_0_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ssk);

    auto g_0_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ssk + 1);

    auto g_0_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ssk + 2);

    auto g_0_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ssk + 3);

    auto g_0_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ssk + 4);

    auto g_0_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ssk + 5);

    auto g_0_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ssk + 6);

    auto g_0_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ssk + 7);

    auto g_0_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ssk + 8);

    auto g_0_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ssk + 9);

    auto g_0_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ssk + 10);

    auto g_0_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ssk + 11);

    auto g_0_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ssk + 12);

    auto g_0_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ssk + 13);

    auto g_0_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ssk + 14);

    auto g_0_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 15);

    auto g_0_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ssk + 16);

    auto g_0_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 17);

    auto g_0_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 18);

    auto g_0_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 19);

    auto g_0_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 20);

    auto g_0_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 21);

    auto g_0_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ssk + 22);

    auto g_0_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 23);

    auto g_0_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 24);

    auto g_0_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 25);

    auto g_0_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 26);

    auto g_0_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 27);

    auto g_0_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 28);

    auto g_0_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ssk + 29);

    auto g_0_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 30);

    auto g_0_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 31);

    auto g_0_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 32);

    auto g_0_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 33);

    auto g_0_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 34);

    auto g_0_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 35);

    /// Set up components of auxilary buffer : SSK

    auto g_0_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_ssk);

    auto g_0_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_ssk + 1);

    auto g_0_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_ssk + 2);

    auto g_0_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_ssk + 3);

    auto g_0_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_ssk + 4);

    auto g_0_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_ssk + 5);

    auto g_0_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_ssk + 6);

    auto g_0_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_ssk + 7);

    auto g_0_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_ssk + 8);

    auto g_0_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_ssk + 9);

    auto g_0_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_ssk + 10);

    auto g_0_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_ssk + 11);

    auto g_0_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_ssk + 12);

    auto g_0_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_ssk + 13);

    auto g_0_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_ssk + 14);

    auto g_0_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 15);

    auto g_0_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_ssk + 16);

    auto g_0_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 17);

    auto g_0_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 18);

    auto g_0_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 19);

    auto g_0_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 20);

    auto g_0_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 21);

    auto g_0_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_ssk + 22);

    auto g_0_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 23);

    auto g_0_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 24);

    auto g_0_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 25);

    auto g_0_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 26);

    auto g_0_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 27);

    auto g_0_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 28);

    auto g_0_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_ssk + 29);

    auto g_0_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 30);

    auto g_0_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 31);

    auto g_0_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 32);

    auto g_0_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 33);

    auto g_0_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 34);

    auto g_0_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 35);

    /// Set up components of auxilary buffer : PSI

    auto g_x_0_xxxxxx_1 = pbuffer.data(idx_eri_1_psi);

    auto g_x_0_xxxxxy_1 = pbuffer.data(idx_eri_1_psi + 1);

    auto g_x_0_xxxxxz_1 = pbuffer.data(idx_eri_1_psi + 2);

    auto g_x_0_xxxxyy_1 = pbuffer.data(idx_eri_1_psi + 3);

    auto g_x_0_xxxxyz_1 = pbuffer.data(idx_eri_1_psi + 4);

    auto g_x_0_xxxxzz_1 = pbuffer.data(idx_eri_1_psi + 5);

    auto g_x_0_xxxyyy_1 = pbuffer.data(idx_eri_1_psi + 6);

    auto g_x_0_xxxyyz_1 = pbuffer.data(idx_eri_1_psi + 7);

    auto g_x_0_xxxyzz_1 = pbuffer.data(idx_eri_1_psi + 8);

    auto g_x_0_xxxzzz_1 = pbuffer.data(idx_eri_1_psi + 9);

    auto g_x_0_xxyyyy_1 = pbuffer.data(idx_eri_1_psi + 10);

    auto g_x_0_xxyyyz_1 = pbuffer.data(idx_eri_1_psi + 11);

    auto g_x_0_xxyyzz_1 = pbuffer.data(idx_eri_1_psi + 12);

    auto g_x_0_xxyzzz_1 = pbuffer.data(idx_eri_1_psi + 13);

    auto g_x_0_xxzzzz_1 = pbuffer.data(idx_eri_1_psi + 14);

    auto g_x_0_xyyyyy_1 = pbuffer.data(idx_eri_1_psi + 15);

    auto g_x_0_xyyyyz_1 = pbuffer.data(idx_eri_1_psi + 16);

    auto g_x_0_xyyyzz_1 = pbuffer.data(idx_eri_1_psi + 17);

    auto g_x_0_xyyzzz_1 = pbuffer.data(idx_eri_1_psi + 18);

    auto g_x_0_xyzzzz_1 = pbuffer.data(idx_eri_1_psi + 19);

    auto g_x_0_xzzzzz_1 = pbuffer.data(idx_eri_1_psi + 20);

    auto g_x_0_yyyyyy_1 = pbuffer.data(idx_eri_1_psi + 21);

    auto g_x_0_yyyyyz_1 = pbuffer.data(idx_eri_1_psi + 22);

    auto g_x_0_yyyyzz_1 = pbuffer.data(idx_eri_1_psi + 23);

    auto g_x_0_yyyzzz_1 = pbuffer.data(idx_eri_1_psi + 24);

    auto g_x_0_yyzzzz_1 = pbuffer.data(idx_eri_1_psi + 25);

    auto g_x_0_yzzzzz_1 = pbuffer.data(idx_eri_1_psi + 26);

    auto g_x_0_zzzzzz_1 = pbuffer.data(idx_eri_1_psi + 27);

    auto g_y_0_xxxxxx_1 = pbuffer.data(idx_eri_1_psi + 28);

    auto g_y_0_xxxxxy_1 = pbuffer.data(idx_eri_1_psi + 29);

    auto g_y_0_xxxxxz_1 = pbuffer.data(idx_eri_1_psi + 30);

    auto g_y_0_xxxxyy_1 = pbuffer.data(idx_eri_1_psi + 31);

    auto g_y_0_xxxxyz_1 = pbuffer.data(idx_eri_1_psi + 32);

    auto g_y_0_xxxxzz_1 = pbuffer.data(idx_eri_1_psi + 33);

    auto g_y_0_xxxyyy_1 = pbuffer.data(idx_eri_1_psi + 34);

    auto g_y_0_xxxyyz_1 = pbuffer.data(idx_eri_1_psi + 35);

    auto g_y_0_xxxyzz_1 = pbuffer.data(idx_eri_1_psi + 36);

    auto g_y_0_xxxzzz_1 = pbuffer.data(idx_eri_1_psi + 37);

    auto g_y_0_xxyyyy_1 = pbuffer.data(idx_eri_1_psi + 38);

    auto g_y_0_xxyyyz_1 = pbuffer.data(idx_eri_1_psi + 39);

    auto g_y_0_xxyyzz_1 = pbuffer.data(idx_eri_1_psi + 40);

    auto g_y_0_xxyzzz_1 = pbuffer.data(idx_eri_1_psi + 41);

    auto g_y_0_xxzzzz_1 = pbuffer.data(idx_eri_1_psi + 42);

    auto g_y_0_xyyyyy_1 = pbuffer.data(idx_eri_1_psi + 43);

    auto g_y_0_xyyyyz_1 = pbuffer.data(idx_eri_1_psi + 44);

    auto g_y_0_xyyyzz_1 = pbuffer.data(idx_eri_1_psi + 45);

    auto g_y_0_xyyzzz_1 = pbuffer.data(idx_eri_1_psi + 46);

    auto g_y_0_xyzzzz_1 = pbuffer.data(idx_eri_1_psi + 47);

    auto g_y_0_xzzzzz_1 = pbuffer.data(idx_eri_1_psi + 48);

    auto g_y_0_yyyyyy_1 = pbuffer.data(idx_eri_1_psi + 49);

    auto g_y_0_yyyyyz_1 = pbuffer.data(idx_eri_1_psi + 50);

    auto g_y_0_yyyyzz_1 = pbuffer.data(idx_eri_1_psi + 51);

    auto g_y_0_yyyzzz_1 = pbuffer.data(idx_eri_1_psi + 52);

    auto g_y_0_yyzzzz_1 = pbuffer.data(idx_eri_1_psi + 53);

    auto g_y_0_yzzzzz_1 = pbuffer.data(idx_eri_1_psi + 54);

    auto g_y_0_zzzzzz_1 = pbuffer.data(idx_eri_1_psi + 55);

    auto g_z_0_xxxxxx_1 = pbuffer.data(idx_eri_1_psi + 56);

    auto g_z_0_xxxxxy_1 = pbuffer.data(idx_eri_1_psi + 57);

    auto g_z_0_xxxxxz_1 = pbuffer.data(idx_eri_1_psi + 58);

    auto g_z_0_xxxxyy_1 = pbuffer.data(idx_eri_1_psi + 59);

    auto g_z_0_xxxxyz_1 = pbuffer.data(idx_eri_1_psi + 60);

    auto g_z_0_xxxxzz_1 = pbuffer.data(idx_eri_1_psi + 61);

    auto g_z_0_xxxyyy_1 = pbuffer.data(idx_eri_1_psi + 62);

    auto g_z_0_xxxyyz_1 = pbuffer.data(idx_eri_1_psi + 63);

    auto g_z_0_xxxyzz_1 = pbuffer.data(idx_eri_1_psi + 64);

    auto g_z_0_xxxzzz_1 = pbuffer.data(idx_eri_1_psi + 65);

    auto g_z_0_xxyyyy_1 = pbuffer.data(idx_eri_1_psi + 66);

    auto g_z_0_xxyyyz_1 = pbuffer.data(idx_eri_1_psi + 67);

    auto g_z_0_xxyyzz_1 = pbuffer.data(idx_eri_1_psi + 68);

    auto g_z_0_xxyzzz_1 = pbuffer.data(idx_eri_1_psi + 69);

    auto g_z_0_xxzzzz_1 = pbuffer.data(idx_eri_1_psi + 70);

    auto g_z_0_xyyyyy_1 = pbuffer.data(idx_eri_1_psi + 71);

    auto g_z_0_xyyyyz_1 = pbuffer.data(idx_eri_1_psi + 72);

    auto g_z_0_xyyyzz_1 = pbuffer.data(idx_eri_1_psi + 73);

    auto g_z_0_xyyzzz_1 = pbuffer.data(idx_eri_1_psi + 74);

    auto g_z_0_xyzzzz_1 = pbuffer.data(idx_eri_1_psi + 75);

    auto g_z_0_xzzzzz_1 = pbuffer.data(idx_eri_1_psi + 76);

    auto g_z_0_yyyyyy_1 = pbuffer.data(idx_eri_1_psi + 77);

    auto g_z_0_yyyyyz_1 = pbuffer.data(idx_eri_1_psi + 78);

    auto g_z_0_yyyyzz_1 = pbuffer.data(idx_eri_1_psi + 79);

    auto g_z_0_yyyzzz_1 = pbuffer.data(idx_eri_1_psi + 80);

    auto g_z_0_yyzzzz_1 = pbuffer.data(idx_eri_1_psi + 81);

    auto g_z_0_yzzzzz_1 = pbuffer.data(idx_eri_1_psi + 82);

    auto g_z_0_zzzzzz_1 = pbuffer.data(idx_eri_1_psi + 83);

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

    /// Set up 0-36 components of targeted buffer : DSK

    auto g_xx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk);

    auto g_xx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 1);

    auto g_xx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 2);

    auto g_xx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 3);

    auto g_xx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 4);

    auto g_xx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 5);

    auto g_xx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 6);

    auto g_xx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 7);

    auto g_xx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 8);

    auto g_xx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 9);

    auto g_xx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 10);

    auto g_xx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 11);

    auto g_xx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 12);

    auto g_xx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 13);

    auto g_xx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 14);

    auto g_xx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 15);

    auto g_xx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 16);

    auto g_xx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 17);

    auto g_xx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 18);

    auto g_xx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 19);

    auto g_xx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 20);

    auto g_xx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 21);

    auto g_xx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 22);

    auto g_xx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 23);

    auto g_xx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 24);

    auto g_xx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 25);

    auto g_xx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 26);

    auto g_xx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 27);

    auto g_xx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 28);

    auto g_xx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 29);

    auto g_xx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 30);

    auto g_xx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 31);

    auto g_xx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 32);

    auto g_xx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 33);

    auto g_xx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 34);

    auto g_xx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 35);

    #pragma omp simd aligned(g_0_0_xxxxxxx_0, g_0_0_xxxxxxx_1, g_0_0_xxxxxxy_0, g_0_0_xxxxxxy_1, g_0_0_xxxxxxz_0, g_0_0_xxxxxxz_1, g_0_0_xxxxxyy_0, g_0_0_xxxxxyy_1, g_0_0_xxxxxyz_0, g_0_0_xxxxxyz_1, g_0_0_xxxxxzz_0, g_0_0_xxxxxzz_1, g_0_0_xxxxyyy_0, g_0_0_xxxxyyy_1, g_0_0_xxxxyyz_0, g_0_0_xxxxyyz_1, g_0_0_xxxxyzz_0, g_0_0_xxxxyzz_1, g_0_0_xxxxzzz_0, g_0_0_xxxxzzz_1, g_0_0_xxxyyyy_0, g_0_0_xxxyyyy_1, g_0_0_xxxyyyz_0, g_0_0_xxxyyyz_1, g_0_0_xxxyyzz_0, g_0_0_xxxyyzz_1, g_0_0_xxxyzzz_0, g_0_0_xxxyzzz_1, g_0_0_xxxzzzz_0, g_0_0_xxxzzzz_1, g_0_0_xxyyyyy_0, g_0_0_xxyyyyy_1, g_0_0_xxyyyyz_0, g_0_0_xxyyyyz_1, g_0_0_xxyyyzz_0, g_0_0_xxyyyzz_1, g_0_0_xxyyzzz_0, g_0_0_xxyyzzz_1, g_0_0_xxyzzzz_0, g_0_0_xxyzzzz_1, g_0_0_xxzzzzz_0, g_0_0_xxzzzzz_1, g_0_0_xyyyyyy_0, g_0_0_xyyyyyy_1, g_0_0_xyyyyyz_0, g_0_0_xyyyyyz_1, g_0_0_xyyyyzz_0, g_0_0_xyyyyzz_1, g_0_0_xyyyzzz_0, g_0_0_xyyyzzz_1, g_0_0_xyyzzzz_0, g_0_0_xyyzzzz_1, g_0_0_xyzzzzz_0, g_0_0_xyzzzzz_1, g_0_0_xzzzzzz_0, g_0_0_xzzzzzz_1, g_0_0_yyyyyyy_0, g_0_0_yyyyyyy_1, g_0_0_yyyyyyz_0, g_0_0_yyyyyyz_1, g_0_0_yyyyyzz_0, g_0_0_yyyyyzz_1, g_0_0_yyyyzzz_0, g_0_0_yyyyzzz_1, g_0_0_yyyzzzz_0, g_0_0_yyyzzzz_1, g_0_0_yyzzzzz_0, g_0_0_yyzzzzz_1, g_0_0_yzzzzzz_0, g_0_0_yzzzzzz_1, g_0_0_zzzzzzz_0, g_0_0_zzzzzzz_1, g_x_0_xxxxxx_1, g_x_0_xxxxxxx_1, g_x_0_xxxxxxy_1, g_x_0_xxxxxxz_1, g_x_0_xxxxxy_1, g_x_0_xxxxxyy_1, g_x_0_xxxxxyz_1, g_x_0_xxxxxz_1, g_x_0_xxxxxzz_1, g_x_0_xxxxyy_1, g_x_0_xxxxyyy_1, g_x_0_xxxxyyz_1, g_x_0_xxxxyz_1, g_x_0_xxxxyzz_1, g_x_0_xxxxzz_1, g_x_0_xxxxzzz_1, g_x_0_xxxyyy_1, g_x_0_xxxyyyy_1, g_x_0_xxxyyyz_1, g_x_0_xxxyyz_1, g_x_0_xxxyyzz_1, g_x_0_xxxyzz_1, g_x_0_xxxyzzz_1, g_x_0_xxxzzz_1, g_x_0_xxxzzzz_1, g_x_0_xxyyyy_1, g_x_0_xxyyyyy_1, g_x_0_xxyyyyz_1, g_x_0_xxyyyz_1, g_x_0_xxyyyzz_1, g_x_0_xxyyzz_1, g_x_0_xxyyzzz_1, g_x_0_xxyzzz_1, g_x_0_xxyzzzz_1, g_x_0_xxzzzz_1, g_x_0_xxzzzzz_1, g_x_0_xyyyyy_1, g_x_0_xyyyyyy_1, g_x_0_xyyyyyz_1, g_x_0_xyyyyz_1, g_x_0_xyyyyzz_1, g_x_0_xyyyzz_1, g_x_0_xyyyzzz_1, g_x_0_xyyzzz_1, g_x_0_xyyzzzz_1, g_x_0_xyzzzz_1, g_x_0_xyzzzzz_1, g_x_0_xzzzzz_1, g_x_0_xzzzzzz_1, g_x_0_yyyyyy_1, g_x_0_yyyyyyy_1, g_x_0_yyyyyyz_1, g_x_0_yyyyyz_1, g_x_0_yyyyyzz_1, g_x_0_yyyyzz_1, g_x_0_yyyyzzz_1, g_x_0_yyyzzz_1, g_x_0_yyyzzzz_1, g_x_0_yyzzzz_1, g_x_0_yyzzzzz_1, g_x_0_yzzzzz_1, g_x_0_yzzzzzz_1, g_x_0_zzzzzz_1, g_x_0_zzzzzzz_1, g_xx_0_xxxxxxx_0, g_xx_0_xxxxxxy_0, g_xx_0_xxxxxxz_0, g_xx_0_xxxxxyy_0, g_xx_0_xxxxxyz_0, g_xx_0_xxxxxzz_0, g_xx_0_xxxxyyy_0, g_xx_0_xxxxyyz_0, g_xx_0_xxxxyzz_0, g_xx_0_xxxxzzz_0, g_xx_0_xxxyyyy_0, g_xx_0_xxxyyyz_0, g_xx_0_xxxyyzz_0, g_xx_0_xxxyzzz_0, g_xx_0_xxxzzzz_0, g_xx_0_xxyyyyy_0, g_xx_0_xxyyyyz_0, g_xx_0_xxyyyzz_0, g_xx_0_xxyyzzz_0, g_xx_0_xxyzzzz_0, g_xx_0_xxzzzzz_0, g_xx_0_xyyyyyy_0, g_xx_0_xyyyyyz_0, g_xx_0_xyyyyzz_0, g_xx_0_xyyyzzz_0, g_xx_0_xyyzzzz_0, g_xx_0_xyzzzzz_0, g_xx_0_xzzzzzz_0, g_xx_0_yyyyyyy_0, g_xx_0_yyyyyyz_0, g_xx_0_yyyyyzz_0, g_xx_0_yyyyzzz_0, g_xx_0_yyyzzzz_0, g_xx_0_yyzzzzz_0, g_xx_0_yzzzzzz_0, g_xx_0_zzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xxxxxxx_0[i] = g_0_0_xxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxx_1[i] * fz_be_0 + 7.0 * g_x_0_xxxxxx_1[i] * fi_acd_0 + g_x_0_xxxxxxx_1[i] * wa_x[i];

        g_xx_0_xxxxxxy_0[i] = g_0_0_xxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxy_1[i] * fi_acd_0 + g_x_0_xxxxxxy_1[i] * wa_x[i];

        g_xx_0_xxxxxxz_0[i] = g_0_0_xxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxz_1[i] * fi_acd_0 + g_x_0_xxxxxxz_1[i] * wa_x[i];

        g_xx_0_xxxxxyy_0[i] = g_0_0_xxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyy_1[i] * fi_acd_0 + g_x_0_xxxxxyy_1[i] * wa_x[i];

        g_xx_0_xxxxxyz_0[i] = g_0_0_xxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyz_1[i] * fi_acd_0 + g_x_0_xxxxxyz_1[i] * wa_x[i];

        g_xx_0_xxxxxzz_0[i] = g_0_0_xxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxzz_1[i] * fi_acd_0 + g_x_0_xxxxxzz_1[i] * wa_x[i];

        g_xx_0_xxxxyyy_0[i] = g_0_0_xxxxyyy_0[i] * fbe_0 - g_0_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyy_1[i] * fi_acd_0 + g_x_0_xxxxyyy_1[i] * wa_x[i];

        g_xx_0_xxxxyyz_0[i] = g_0_0_xxxxyyz_0[i] * fbe_0 - g_0_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyz_1[i] * fi_acd_0 + g_x_0_xxxxyyz_1[i] * wa_x[i];

        g_xx_0_xxxxyzz_0[i] = g_0_0_xxxxyzz_0[i] * fbe_0 - g_0_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyzz_1[i] * fi_acd_0 + g_x_0_xxxxyzz_1[i] * wa_x[i];

        g_xx_0_xxxxzzz_0[i] = g_0_0_xxxxzzz_0[i] * fbe_0 - g_0_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxzzz_1[i] * fi_acd_0 + g_x_0_xxxxzzz_1[i] * wa_x[i];

        g_xx_0_xxxyyyy_0[i] = g_0_0_xxxyyyy_0[i] * fbe_0 - g_0_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyy_1[i] * fi_acd_0 + g_x_0_xxxyyyy_1[i] * wa_x[i];

        g_xx_0_xxxyyyz_0[i] = g_0_0_xxxyyyz_0[i] * fbe_0 - g_0_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyz_1[i] * fi_acd_0 + g_x_0_xxxyyyz_1[i] * wa_x[i];

        g_xx_0_xxxyyzz_0[i] = g_0_0_xxxyyzz_0[i] * fbe_0 - g_0_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyzz_1[i] * fi_acd_0 + g_x_0_xxxyyzz_1[i] * wa_x[i];

        g_xx_0_xxxyzzz_0[i] = g_0_0_xxxyzzz_0[i] * fbe_0 - g_0_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyzzz_1[i] * fi_acd_0 + g_x_0_xxxyzzz_1[i] * wa_x[i];

        g_xx_0_xxxzzzz_0[i] = g_0_0_xxxzzzz_0[i] * fbe_0 - g_0_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxzzzz_1[i] * fi_acd_0 + g_x_0_xxxzzzz_1[i] * wa_x[i];

        g_xx_0_xxyyyyy_0[i] = g_0_0_xxyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyy_1[i] * fi_acd_0 + g_x_0_xxyyyyy_1[i] * wa_x[i];

        g_xx_0_xxyyyyz_0[i] = g_0_0_xxyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyz_1[i] * fi_acd_0 + g_x_0_xxyyyyz_1[i] * wa_x[i];

        g_xx_0_xxyyyzz_0[i] = g_0_0_xxyyyzz_0[i] * fbe_0 - g_0_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyzz_1[i] * fi_acd_0 + g_x_0_xxyyyzz_1[i] * wa_x[i];

        g_xx_0_xxyyzzz_0[i] = g_0_0_xxyyzzz_0[i] * fbe_0 - g_0_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyzzz_1[i] * fi_acd_0 + g_x_0_xxyyzzz_1[i] * wa_x[i];

        g_xx_0_xxyzzzz_0[i] = g_0_0_xxyzzzz_0[i] * fbe_0 - g_0_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyzzzz_1[i] * fi_acd_0 + g_x_0_xxyzzzz_1[i] * wa_x[i];

        g_xx_0_xxzzzzz_0[i] = g_0_0_xxzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xzzzzz_1[i] * fi_acd_0 + g_x_0_xxzzzzz_1[i] * wa_x[i];

        g_xx_0_xyyyyyy_0[i] = g_0_0_xyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyy_1[i] * fz_be_0 + g_x_0_yyyyyy_1[i] * fi_acd_0 + g_x_0_xyyyyyy_1[i] * wa_x[i];

        g_xx_0_xyyyyyz_0[i] = g_0_0_xyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyz_1[i] * fz_be_0 + g_x_0_yyyyyz_1[i] * fi_acd_0 + g_x_0_xyyyyyz_1[i] * wa_x[i];

        g_xx_0_xyyyyzz_0[i] = g_0_0_xyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyzz_1[i] * fz_be_0 + g_x_0_yyyyzz_1[i] * fi_acd_0 + g_x_0_xyyyyzz_1[i] * wa_x[i];

        g_xx_0_xyyyzzz_0[i] = g_0_0_xyyyzzz_0[i] * fbe_0 - g_0_0_xyyyzzz_1[i] * fz_be_0 + g_x_0_yyyzzz_1[i] * fi_acd_0 + g_x_0_xyyyzzz_1[i] * wa_x[i];

        g_xx_0_xyyzzzz_0[i] = g_0_0_xyyzzzz_0[i] * fbe_0 - g_0_0_xyyzzzz_1[i] * fz_be_0 + g_x_0_yyzzzz_1[i] * fi_acd_0 + g_x_0_xyyzzzz_1[i] * wa_x[i];

        g_xx_0_xyzzzzz_0[i] = g_0_0_xyzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzz_1[i] * fz_be_0 + g_x_0_yzzzzz_1[i] * fi_acd_0 + g_x_0_xyzzzzz_1[i] * wa_x[i];

        g_xx_0_xzzzzzz_0[i] = g_0_0_xzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzz_1[i] * fz_be_0 + g_x_0_zzzzzz_1[i] * fi_acd_0 + g_x_0_xzzzzzz_1[i] * wa_x[i];

        g_xx_0_yyyyyyy_0[i] = g_0_0_yyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyy_1[i] * fz_be_0 + g_x_0_yyyyyyy_1[i] * wa_x[i];

        g_xx_0_yyyyyyz_0[i] = g_0_0_yyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyz_1[i] * fz_be_0 + g_x_0_yyyyyyz_1[i] * wa_x[i];

        g_xx_0_yyyyyzz_0[i] = g_0_0_yyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyzz_1[i] * fz_be_0 + g_x_0_yyyyyzz_1[i] * wa_x[i];

        g_xx_0_yyyyzzz_0[i] = g_0_0_yyyyzzz_0[i] * fbe_0 - g_0_0_yyyyzzz_1[i] * fz_be_0 + g_x_0_yyyyzzz_1[i] * wa_x[i];

        g_xx_0_yyyzzzz_0[i] = g_0_0_yyyzzzz_0[i] * fbe_0 - g_0_0_yyyzzzz_1[i] * fz_be_0 + g_x_0_yyyzzzz_1[i] * wa_x[i];

        g_xx_0_yyzzzzz_0[i] = g_0_0_yyzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzz_1[i] * fz_be_0 + g_x_0_yyzzzzz_1[i] * wa_x[i];

        g_xx_0_yzzzzzz_0[i] = g_0_0_yzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzz_1[i] * fz_be_0 + g_x_0_yzzzzzz_1[i] * wa_x[i];

        g_xx_0_zzzzzzz_0[i] = g_0_0_zzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzz_1[i] * fz_be_0 + g_x_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 36-72 components of targeted buffer : DSK

    auto g_xy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk + 36);

    auto g_xy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 37);

    auto g_xy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 38);

    auto g_xy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 39);

    auto g_xy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 40);

    auto g_xy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 41);

    auto g_xy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 42);

    auto g_xy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 43);

    auto g_xy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 44);

    auto g_xy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 45);

    auto g_xy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 46);

    auto g_xy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 47);

    auto g_xy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 48);

    auto g_xy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 49);

    auto g_xy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 50);

    auto g_xy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 51);

    auto g_xy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 52);

    auto g_xy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 53);

    auto g_xy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 54);

    auto g_xy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 55);

    auto g_xy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 56);

    auto g_xy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 57);

    auto g_xy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 58);

    auto g_xy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 59);

    auto g_xy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 60);

    auto g_xy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 61);

    auto g_xy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 62);

    auto g_xy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 63);

    auto g_xy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 64);

    auto g_xy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 65);

    auto g_xy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 66);

    auto g_xy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 67);

    auto g_xy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 68);

    auto g_xy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 69);

    auto g_xy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 70);

    auto g_xy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 71);

    #pragma omp simd aligned(g_x_0_xxxxxxx_1, g_x_0_xxxxxxz_1, g_x_0_xxxxxzz_1, g_x_0_xxxxzzz_1, g_x_0_xxxzzzz_1, g_x_0_xxzzzzz_1, g_x_0_xzzzzzz_1, g_xy_0_xxxxxxx_0, g_xy_0_xxxxxxy_0, g_xy_0_xxxxxxz_0, g_xy_0_xxxxxyy_0, g_xy_0_xxxxxyz_0, g_xy_0_xxxxxzz_0, g_xy_0_xxxxyyy_0, g_xy_0_xxxxyyz_0, g_xy_0_xxxxyzz_0, g_xy_0_xxxxzzz_0, g_xy_0_xxxyyyy_0, g_xy_0_xxxyyyz_0, g_xy_0_xxxyyzz_0, g_xy_0_xxxyzzz_0, g_xy_0_xxxzzzz_0, g_xy_0_xxyyyyy_0, g_xy_0_xxyyyyz_0, g_xy_0_xxyyyzz_0, g_xy_0_xxyyzzz_0, g_xy_0_xxyzzzz_0, g_xy_0_xxzzzzz_0, g_xy_0_xyyyyyy_0, g_xy_0_xyyyyyz_0, g_xy_0_xyyyyzz_0, g_xy_0_xyyyzzz_0, g_xy_0_xyyzzzz_0, g_xy_0_xyzzzzz_0, g_xy_0_xzzzzzz_0, g_xy_0_yyyyyyy_0, g_xy_0_yyyyyyz_0, g_xy_0_yyyyyzz_0, g_xy_0_yyyyzzz_0, g_xy_0_yyyzzzz_0, g_xy_0_yyzzzzz_0, g_xy_0_yzzzzzz_0, g_xy_0_zzzzzzz_0, g_y_0_xxxxxxy_1, g_y_0_xxxxxy_1, g_y_0_xxxxxyy_1, g_y_0_xxxxxyz_1, g_y_0_xxxxyy_1, g_y_0_xxxxyyy_1, g_y_0_xxxxyyz_1, g_y_0_xxxxyz_1, g_y_0_xxxxyzz_1, g_y_0_xxxyyy_1, g_y_0_xxxyyyy_1, g_y_0_xxxyyyz_1, g_y_0_xxxyyz_1, g_y_0_xxxyyzz_1, g_y_0_xxxyzz_1, g_y_0_xxxyzzz_1, g_y_0_xxyyyy_1, g_y_0_xxyyyyy_1, g_y_0_xxyyyyz_1, g_y_0_xxyyyz_1, g_y_0_xxyyyzz_1, g_y_0_xxyyzz_1, g_y_0_xxyyzzz_1, g_y_0_xxyzzz_1, g_y_0_xxyzzzz_1, g_y_0_xyyyyy_1, g_y_0_xyyyyyy_1, g_y_0_xyyyyyz_1, g_y_0_xyyyyz_1, g_y_0_xyyyyzz_1, g_y_0_xyyyzz_1, g_y_0_xyyyzzz_1, g_y_0_xyyzzz_1, g_y_0_xyyzzzz_1, g_y_0_xyzzzz_1, g_y_0_xyzzzzz_1, g_y_0_yyyyyy_1, g_y_0_yyyyyyy_1, g_y_0_yyyyyyz_1, g_y_0_yyyyyz_1, g_y_0_yyyyyzz_1, g_y_0_yyyyzz_1, g_y_0_yyyyzzz_1, g_y_0_yyyzzz_1, g_y_0_yyyzzzz_1, g_y_0_yyzzzz_1, g_y_0_yyzzzzz_1, g_y_0_yzzzzz_1, g_y_0_yzzzzzz_1, g_y_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xxxxxxx_0[i] = g_x_0_xxxxxxx_1[i] * wa_y[i];

        g_xy_0_xxxxxxy_0[i] = 6.0 * g_y_0_xxxxxy_1[i] * fi_acd_0 + g_y_0_xxxxxxy_1[i] * wa_x[i];

        g_xy_0_xxxxxxz_0[i] = g_x_0_xxxxxxz_1[i] * wa_y[i];

        g_xy_0_xxxxxyy_0[i] = 5.0 * g_y_0_xxxxyy_1[i] * fi_acd_0 + g_y_0_xxxxxyy_1[i] * wa_x[i];

        g_xy_0_xxxxxyz_0[i] = 5.0 * g_y_0_xxxxyz_1[i] * fi_acd_0 + g_y_0_xxxxxyz_1[i] * wa_x[i];

        g_xy_0_xxxxxzz_0[i] = g_x_0_xxxxxzz_1[i] * wa_y[i];

        g_xy_0_xxxxyyy_0[i] = 4.0 * g_y_0_xxxyyy_1[i] * fi_acd_0 + g_y_0_xxxxyyy_1[i] * wa_x[i];

        g_xy_0_xxxxyyz_0[i] = 4.0 * g_y_0_xxxyyz_1[i] * fi_acd_0 + g_y_0_xxxxyyz_1[i] * wa_x[i];

        g_xy_0_xxxxyzz_0[i] = 4.0 * g_y_0_xxxyzz_1[i] * fi_acd_0 + g_y_0_xxxxyzz_1[i] * wa_x[i];

        g_xy_0_xxxxzzz_0[i] = g_x_0_xxxxzzz_1[i] * wa_y[i];

        g_xy_0_xxxyyyy_0[i] = 3.0 * g_y_0_xxyyyy_1[i] * fi_acd_0 + g_y_0_xxxyyyy_1[i] * wa_x[i];

        g_xy_0_xxxyyyz_0[i] = 3.0 * g_y_0_xxyyyz_1[i] * fi_acd_0 + g_y_0_xxxyyyz_1[i] * wa_x[i];

        g_xy_0_xxxyyzz_0[i] = 3.0 * g_y_0_xxyyzz_1[i] * fi_acd_0 + g_y_0_xxxyyzz_1[i] * wa_x[i];

        g_xy_0_xxxyzzz_0[i] = 3.0 * g_y_0_xxyzzz_1[i] * fi_acd_0 + g_y_0_xxxyzzz_1[i] * wa_x[i];

        g_xy_0_xxxzzzz_0[i] = g_x_0_xxxzzzz_1[i] * wa_y[i];

        g_xy_0_xxyyyyy_0[i] = 2.0 * g_y_0_xyyyyy_1[i] * fi_acd_0 + g_y_0_xxyyyyy_1[i] * wa_x[i];

        g_xy_0_xxyyyyz_0[i] = 2.0 * g_y_0_xyyyyz_1[i] * fi_acd_0 + g_y_0_xxyyyyz_1[i] * wa_x[i];

        g_xy_0_xxyyyzz_0[i] = 2.0 * g_y_0_xyyyzz_1[i] * fi_acd_0 + g_y_0_xxyyyzz_1[i] * wa_x[i];

        g_xy_0_xxyyzzz_0[i] = 2.0 * g_y_0_xyyzzz_1[i] * fi_acd_0 + g_y_0_xxyyzzz_1[i] * wa_x[i];

        g_xy_0_xxyzzzz_0[i] = 2.0 * g_y_0_xyzzzz_1[i] * fi_acd_0 + g_y_0_xxyzzzz_1[i] * wa_x[i];

        g_xy_0_xxzzzzz_0[i] = g_x_0_xxzzzzz_1[i] * wa_y[i];

        g_xy_0_xyyyyyy_0[i] = g_y_0_yyyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyyy_1[i] * wa_x[i];

        g_xy_0_xyyyyyz_0[i] = g_y_0_yyyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyyz_1[i] * wa_x[i];

        g_xy_0_xyyyyzz_0[i] = g_y_0_yyyyzz_1[i] * fi_acd_0 + g_y_0_xyyyyzz_1[i] * wa_x[i];

        g_xy_0_xyyyzzz_0[i] = g_y_0_yyyzzz_1[i] * fi_acd_0 + g_y_0_xyyyzzz_1[i] * wa_x[i];

        g_xy_0_xyyzzzz_0[i] = g_y_0_yyzzzz_1[i] * fi_acd_0 + g_y_0_xyyzzzz_1[i] * wa_x[i];

        g_xy_0_xyzzzzz_0[i] = g_y_0_yzzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzzz_1[i] * wa_x[i];

        g_xy_0_xzzzzzz_0[i] = g_x_0_xzzzzzz_1[i] * wa_y[i];

        g_xy_0_yyyyyyy_0[i] = g_y_0_yyyyyyy_1[i] * wa_x[i];

        g_xy_0_yyyyyyz_0[i] = g_y_0_yyyyyyz_1[i] * wa_x[i];

        g_xy_0_yyyyyzz_0[i] = g_y_0_yyyyyzz_1[i] * wa_x[i];

        g_xy_0_yyyyzzz_0[i] = g_y_0_yyyyzzz_1[i] * wa_x[i];

        g_xy_0_yyyzzzz_0[i] = g_y_0_yyyzzzz_1[i] * wa_x[i];

        g_xy_0_yyzzzzz_0[i] = g_y_0_yyzzzzz_1[i] * wa_x[i];

        g_xy_0_yzzzzzz_0[i] = g_y_0_yzzzzzz_1[i] * wa_x[i];

        g_xy_0_zzzzzzz_0[i] = g_y_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 72-108 components of targeted buffer : DSK

    auto g_xz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk + 72);

    auto g_xz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 73);

    auto g_xz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 74);

    auto g_xz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 75);

    auto g_xz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 76);

    auto g_xz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 77);

    auto g_xz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 78);

    auto g_xz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 79);

    auto g_xz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 80);

    auto g_xz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 81);

    auto g_xz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 82);

    auto g_xz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 83);

    auto g_xz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 84);

    auto g_xz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 85);

    auto g_xz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 86);

    auto g_xz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 87);

    auto g_xz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 88);

    auto g_xz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 89);

    auto g_xz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 90);

    auto g_xz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 91);

    auto g_xz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 92);

    auto g_xz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 93);

    auto g_xz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 94);

    auto g_xz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 95);

    auto g_xz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 96);

    auto g_xz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 97);

    auto g_xz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 98);

    auto g_xz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 99);

    auto g_xz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 100);

    auto g_xz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 101);

    auto g_xz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 102);

    auto g_xz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 103);

    auto g_xz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 104);

    auto g_xz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 105);

    auto g_xz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 106);

    auto g_xz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 107);

    #pragma omp simd aligned(g_x_0_xxxxxxx_1, g_x_0_xxxxxxy_1, g_x_0_xxxxxyy_1, g_x_0_xxxxyyy_1, g_x_0_xxxyyyy_1, g_x_0_xxyyyyy_1, g_x_0_xyyyyyy_1, g_xz_0_xxxxxxx_0, g_xz_0_xxxxxxy_0, g_xz_0_xxxxxxz_0, g_xz_0_xxxxxyy_0, g_xz_0_xxxxxyz_0, g_xz_0_xxxxxzz_0, g_xz_0_xxxxyyy_0, g_xz_0_xxxxyyz_0, g_xz_0_xxxxyzz_0, g_xz_0_xxxxzzz_0, g_xz_0_xxxyyyy_0, g_xz_0_xxxyyyz_0, g_xz_0_xxxyyzz_0, g_xz_0_xxxyzzz_0, g_xz_0_xxxzzzz_0, g_xz_0_xxyyyyy_0, g_xz_0_xxyyyyz_0, g_xz_0_xxyyyzz_0, g_xz_0_xxyyzzz_0, g_xz_0_xxyzzzz_0, g_xz_0_xxzzzzz_0, g_xz_0_xyyyyyy_0, g_xz_0_xyyyyyz_0, g_xz_0_xyyyyzz_0, g_xz_0_xyyyzzz_0, g_xz_0_xyyzzzz_0, g_xz_0_xyzzzzz_0, g_xz_0_xzzzzzz_0, g_xz_0_yyyyyyy_0, g_xz_0_yyyyyyz_0, g_xz_0_yyyyyzz_0, g_xz_0_yyyyzzz_0, g_xz_0_yyyzzzz_0, g_xz_0_yyzzzzz_0, g_xz_0_yzzzzzz_0, g_xz_0_zzzzzzz_0, g_z_0_xxxxxxz_1, g_z_0_xxxxxyz_1, g_z_0_xxxxxz_1, g_z_0_xxxxxzz_1, g_z_0_xxxxyyz_1, g_z_0_xxxxyz_1, g_z_0_xxxxyzz_1, g_z_0_xxxxzz_1, g_z_0_xxxxzzz_1, g_z_0_xxxyyyz_1, g_z_0_xxxyyz_1, g_z_0_xxxyyzz_1, g_z_0_xxxyzz_1, g_z_0_xxxyzzz_1, g_z_0_xxxzzz_1, g_z_0_xxxzzzz_1, g_z_0_xxyyyyz_1, g_z_0_xxyyyz_1, g_z_0_xxyyyzz_1, g_z_0_xxyyzz_1, g_z_0_xxyyzzz_1, g_z_0_xxyzzz_1, g_z_0_xxyzzzz_1, g_z_0_xxzzzz_1, g_z_0_xxzzzzz_1, g_z_0_xyyyyyz_1, g_z_0_xyyyyz_1, g_z_0_xyyyyzz_1, g_z_0_xyyyzz_1, g_z_0_xyyyzzz_1, g_z_0_xyyzzz_1, g_z_0_xyyzzzz_1, g_z_0_xyzzzz_1, g_z_0_xyzzzzz_1, g_z_0_xzzzzz_1, g_z_0_xzzzzzz_1, g_z_0_yyyyyyy_1, g_z_0_yyyyyyz_1, g_z_0_yyyyyz_1, g_z_0_yyyyyzz_1, g_z_0_yyyyzz_1, g_z_0_yyyyzzz_1, g_z_0_yyyzzz_1, g_z_0_yyyzzzz_1, g_z_0_yyzzzz_1, g_z_0_yyzzzzz_1, g_z_0_yzzzzz_1, g_z_0_yzzzzzz_1, g_z_0_zzzzzz_1, g_z_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xxxxxxx_0[i] = g_x_0_xxxxxxx_1[i] * wa_z[i];

        g_xz_0_xxxxxxy_0[i] = g_x_0_xxxxxxy_1[i] * wa_z[i];

        g_xz_0_xxxxxxz_0[i] = 6.0 * g_z_0_xxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxxz_1[i] * wa_x[i];

        g_xz_0_xxxxxyy_0[i] = g_x_0_xxxxxyy_1[i] * wa_z[i];

        g_xz_0_xxxxxyz_0[i] = 5.0 * g_z_0_xxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxxyz_1[i] * wa_x[i];

        g_xz_0_xxxxxzz_0[i] = 5.0 * g_z_0_xxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxxzz_1[i] * wa_x[i];

        g_xz_0_xxxxyyy_0[i] = g_x_0_xxxxyyy_1[i] * wa_z[i];

        g_xz_0_xxxxyyz_0[i] = 4.0 * g_z_0_xxxyyz_1[i] * fi_acd_0 + g_z_0_xxxxyyz_1[i] * wa_x[i];

        g_xz_0_xxxxyzz_0[i] = 4.0 * g_z_0_xxxyzz_1[i] * fi_acd_0 + g_z_0_xxxxyzz_1[i] * wa_x[i];

        g_xz_0_xxxxzzz_0[i] = 4.0 * g_z_0_xxxzzz_1[i] * fi_acd_0 + g_z_0_xxxxzzz_1[i] * wa_x[i];

        g_xz_0_xxxyyyy_0[i] = g_x_0_xxxyyyy_1[i] * wa_z[i];

        g_xz_0_xxxyyyz_0[i] = 3.0 * g_z_0_xxyyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyz_1[i] * wa_x[i];

        g_xz_0_xxxyyzz_0[i] = 3.0 * g_z_0_xxyyzz_1[i] * fi_acd_0 + g_z_0_xxxyyzz_1[i] * wa_x[i];

        g_xz_0_xxxyzzz_0[i] = 3.0 * g_z_0_xxyzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzz_1[i] * wa_x[i];

        g_xz_0_xxxzzzz_0[i] = 3.0 * g_z_0_xxzzzz_1[i] * fi_acd_0 + g_z_0_xxxzzzz_1[i] * wa_x[i];

        g_xz_0_xxyyyyy_0[i] = g_x_0_xxyyyyy_1[i] * wa_z[i];

        g_xz_0_xxyyyyz_0[i] = 2.0 * g_z_0_xyyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyz_1[i] * wa_x[i];

        g_xz_0_xxyyyzz_0[i] = 2.0 * g_z_0_xyyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyzz_1[i] * wa_x[i];

        g_xz_0_xxyyzzz_0[i] = 2.0 * g_z_0_xyyzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzz_1[i] * wa_x[i];

        g_xz_0_xxyzzzz_0[i] = 2.0 * g_z_0_xyzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzz_1[i] * wa_x[i];

        g_xz_0_xxzzzzz_0[i] = 2.0 * g_z_0_xzzzzz_1[i] * fi_acd_0 + g_z_0_xxzzzzz_1[i] * wa_x[i];

        g_xz_0_xyyyyyy_0[i] = g_x_0_xyyyyyy_1[i] * wa_z[i];

        g_xz_0_xyyyyyz_0[i] = g_z_0_yyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyz_1[i] * wa_x[i];

        g_xz_0_xyyyyzz_0[i] = g_z_0_yyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyzz_1[i] * wa_x[i];

        g_xz_0_xyyyzzz_0[i] = g_z_0_yyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzz_1[i] * wa_x[i];

        g_xz_0_xyyzzzz_0[i] = g_z_0_yyzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzz_1[i] * wa_x[i];

        g_xz_0_xyzzzzz_0[i] = g_z_0_yzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzz_1[i] * wa_x[i];

        g_xz_0_xzzzzzz_0[i] = g_z_0_zzzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzzz_1[i] * wa_x[i];

        g_xz_0_yyyyyyy_0[i] = g_z_0_yyyyyyy_1[i] * wa_x[i];

        g_xz_0_yyyyyyz_0[i] = g_z_0_yyyyyyz_1[i] * wa_x[i];

        g_xz_0_yyyyyzz_0[i] = g_z_0_yyyyyzz_1[i] * wa_x[i];

        g_xz_0_yyyyzzz_0[i] = g_z_0_yyyyzzz_1[i] * wa_x[i];

        g_xz_0_yyyzzzz_0[i] = g_z_0_yyyzzzz_1[i] * wa_x[i];

        g_xz_0_yyzzzzz_0[i] = g_z_0_yyzzzzz_1[i] * wa_x[i];

        g_xz_0_yzzzzzz_0[i] = g_z_0_yzzzzzz_1[i] * wa_x[i];

        g_xz_0_zzzzzzz_0[i] = g_z_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 108-144 components of targeted buffer : DSK

    auto g_yy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk + 108);

    auto g_yy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 109);

    auto g_yy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 110);

    auto g_yy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 111);

    auto g_yy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 112);

    auto g_yy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 113);

    auto g_yy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 114);

    auto g_yy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 115);

    auto g_yy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 116);

    auto g_yy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 117);

    auto g_yy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 118);

    auto g_yy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 119);

    auto g_yy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 120);

    auto g_yy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 121);

    auto g_yy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 122);

    auto g_yy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 123);

    auto g_yy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 124);

    auto g_yy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 125);

    auto g_yy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 126);

    auto g_yy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 127);

    auto g_yy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 128);

    auto g_yy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 129);

    auto g_yy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 130);

    auto g_yy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 131);

    auto g_yy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 132);

    auto g_yy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 133);

    auto g_yy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 134);

    auto g_yy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 135);

    auto g_yy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 136);

    auto g_yy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 137);

    auto g_yy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 138);

    auto g_yy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 139);

    auto g_yy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 140);

    auto g_yy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 141);

    auto g_yy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 142);

    auto g_yy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 143);

    #pragma omp simd aligned(g_0_0_xxxxxxx_0, g_0_0_xxxxxxx_1, g_0_0_xxxxxxy_0, g_0_0_xxxxxxy_1, g_0_0_xxxxxxz_0, g_0_0_xxxxxxz_1, g_0_0_xxxxxyy_0, g_0_0_xxxxxyy_1, g_0_0_xxxxxyz_0, g_0_0_xxxxxyz_1, g_0_0_xxxxxzz_0, g_0_0_xxxxxzz_1, g_0_0_xxxxyyy_0, g_0_0_xxxxyyy_1, g_0_0_xxxxyyz_0, g_0_0_xxxxyyz_1, g_0_0_xxxxyzz_0, g_0_0_xxxxyzz_1, g_0_0_xxxxzzz_0, g_0_0_xxxxzzz_1, g_0_0_xxxyyyy_0, g_0_0_xxxyyyy_1, g_0_0_xxxyyyz_0, g_0_0_xxxyyyz_1, g_0_0_xxxyyzz_0, g_0_0_xxxyyzz_1, g_0_0_xxxyzzz_0, g_0_0_xxxyzzz_1, g_0_0_xxxzzzz_0, g_0_0_xxxzzzz_1, g_0_0_xxyyyyy_0, g_0_0_xxyyyyy_1, g_0_0_xxyyyyz_0, g_0_0_xxyyyyz_1, g_0_0_xxyyyzz_0, g_0_0_xxyyyzz_1, g_0_0_xxyyzzz_0, g_0_0_xxyyzzz_1, g_0_0_xxyzzzz_0, g_0_0_xxyzzzz_1, g_0_0_xxzzzzz_0, g_0_0_xxzzzzz_1, g_0_0_xyyyyyy_0, g_0_0_xyyyyyy_1, g_0_0_xyyyyyz_0, g_0_0_xyyyyyz_1, g_0_0_xyyyyzz_0, g_0_0_xyyyyzz_1, g_0_0_xyyyzzz_0, g_0_0_xyyyzzz_1, g_0_0_xyyzzzz_0, g_0_0_xyyzzzz_1, g_0_0_xyzzzzz_0, g_0_0_xyzzzzz_1, g_0_0_xzzzzzz_0, g_0_0_xzzzzzz_1, g_0_0_yyyyyyy_0, g_0_0_yyyyyyy_1, g_0_0_yyyyyyz_0, g_0_0_yyyyyyz_1, g_0_0_yyyyyzz_0, g_0_0_yyyyyzz_1, g_0_0_yyyyzzz_0, g_0_0_yyyyzzz_1, g_0_0_yyyzzzz_0, g_0_0_yyyzzzz_1, g_0_0_yyzzzzz_0, g_0_0_yyzzzzz_1, g_0_0_yzzzzzz_0, g_0_0_yzzzzzz_1, g_0_0_zzzzzzz_0, g_0_0_zzzzzzz_1, g_y_0_xxxxxx_1, g_y_0_xxxxxxx_1, g_y_0_xxxxxxy_1, g_y_0_xxxxxxz_1, g_y_0_xxxxxy_1, g_y_0_xxxxxyy_1, g_y_0_xxxxxyz_1, g_y_0_xxxxxz_1, g_y_0_xxxxxzz_1, g_y_0_xxxxyy_1, g_y_0_xxxxyyy_1, g_y_0_xxxxyyz_1, g_y_0_xxxxyz_1, g_y_0_xxxxyzz_1, g_y_0_xxxxzz_1, g_y_0_xxxxzzz_1, g_y_0_xxxyyy_1, g_y_0_xxxyyyy_1, g_y_0_xxxyyyz_1, g_y_0_xxxyyz_1, g_y_0_xxxyyzz_1, g_y_0_xxxyzz_1, g_y_0_xxxyzzz_1, g_y_0_xxxzzz_1, g_y_0_xxxzzzz_1, g_y_0_xxyyyy_1, g_y_0_xxyyyyy_1, g_y_0_xxyyyyz_1, g_y_0_xxyyyz_1, g_y_0_xxyyyzz_1, g_y_0_xxyyzz_1, g_y_0_xxyyzzz_1, g_y_0_xxyzzz_1, g_y_0_xxyzzzz_1, g_y_0_xxzzzz_1, g_y_0_xxzzzzz_1, g_y_0_xyyyyy_1, g_y_0_xyyyyyy_1, g_y_0_xyyyyyz_1, g_y_0_xyyyyz_1, g_y_0_xyyyyzz_1, g_y_0_xyyyzz_1, g_y_0_xyyyzzz_1, g_y_0_xyyzzz_1, g_y_0_xyyzzzz_1, g_y_0_xyzzzz_1, g_y_0_xyzzzzz_1, g_y_0_xzzzzz_1, g_y_0_xzzzzzz_1, g_y_0_yyyyyy_1, g_y_0_yyyyyyy_1, g_y_0_yyyyyyz_1, g_y_0_yyyyyz_1, g_y_0_yyyyyzz_1, g_y_0_yyyyzz_1, g_y_0_yyyyzzz_1, g_y_0_yyyzzz_1, g_y_0_yyyzzzz_1, g_y_0_yyzzzz_1, g_y_0_yyzzzzz_1, g_y_0_yzzzzz_1, g_y_0_yzzzzzz_1, g_y_0_zzzzzz_1, g_y_0_zzzzzzz_1, g_yy_0_xxxxxxx_0, g_yy_0_xxxxxxy_0, g_yy_0_xxxxxxz_0, g_yy_0_xxxxxyy_0, g_yy_0_xxxxxyz_0, g_yy_0_xxxxxzz_0, g_yy_0_xxxxyyy_0, g_yy_0_xxxxyyz_0, g_yy_0_xxxxyzz_0, g_yy_0_xxxxzzz_0, g_yy_0_xxxyyyy_0, g_yy_0_xxxyyyz_0, g_yy_0_xxxyyzz_0, g_yy_0_xxxyzzz_0, g_yy_0_xxxzzzz_0, g_yy_0_xxyyyyy_0, g_yy_0_xxyyyyz_0, g_yy_0_xxyyyzz_0, g_yy_0_xxyyzzz_0, g_yy_0_xxyzzzz_0, g_yy_0_xxzzzzz_0, g_yy_0_xyyyyyy_0, g_yy_0_xyyyyyz_0, g_yy_0_xyyyyzz_0, g_yy_0_xyyyzzz_0, g_yy_0_xyyzzzz_0, g_yy_0_xyzzzzz_0, g_yy_0_xzzzzzz_0, g_yy_0_yyyyyyy_0, g_yy_0_yyyyyyz_0, g_yy_0_yyyyyzz_0, g_yy_0_yyyyzzz_0, g_yy_0_yyyzzzz_0, g_yy_0_yyzzzzz_0, g_yy_0_yzzzzzz_0, g_yy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xxxxxxx_0[i] = g_0_0_xxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxx_1[i] * fz_be_0 + g_y_0_xxxxxxx_1[i] * wa_y[i];

        g_yy_0_xxxxxxy_0[i] = g_0_0_xxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxy_1[i] * fz_be_0 + g_y_0_xxxxxx_1[i] * fi_acd_0 + g_y_0_xxxxxxy_1[i] * wa_y[i];

        g_yy_0_xxxxxxz_0[i] = g_0_0_xxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxz_1[i] * fz_be_0 + g_y_0_xxxxxxz_1[i] * wa_y[i];

        g_yy_0_xxxxxyy_0[i] = g_0_0_xxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxyy_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxxy_1[i] * fi_acd_0 + g_y_0_xxxxxyy_1[i] * wa_y[i];

        g_yy_0_xxxxxyz_0[i] = g_0_0_xxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxyz_1[i] * fz_be_0 + g_y_0_xxxxxz_1[i] * fi_acd_0 + g_y_0_xxxxxyz_1[i] * wa_y[i];

        g_yy_0_xxxxxzz_0[i] = g_0_0_xxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxzz_1[i] * fz_be_0 + g_y_0_xxxxxzz_1[i] * wa_y[i];

        g_yy_0_xxxxyyy_0[i] = g_0_0_xxxxyyy_0[i] * fbe_0 - g_0_0_xxxxyyy_1[i] * fz_be_0 + 3.0 * g_y_0_xxxxyy_1[i] * fi_acd_0 + g_y_0_xxxxyyy_1[i] * wa_y[i];

        g_yy_0_xxxxyyz_0[i] = g_0_0_xxxxyyz_0[i] * fbe_0 - g_0_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxyz_1[i] * fi_acd_0 + g_y_0_xxxxyyz_1[i] * wa_y[i];

        g_yy_0_xxxxyzz_0[i] = g_0_0_xxxxyzz_0[i] * fbe_0 - g_0_0_xxxxyzz_1[i] * fz_be_0 + g_y_0_xxxxzz_1[i] * fi_acd_0 + g_y_0_xxxxyzz_1[i] * wa_y[i];

        g_yy_0_xxxxzzz_0[i] = g_0_0_xxxxzzz_0[i] * fbe_0 - g_0_0_xxxxzzz_1[i] * fz_be_0 + g_y_0_xxxxzzz_1[i] * wa_y[i];

        g_yy_0_xxxyyyy_0[i] = g_0_0_xxxyyyy_0[i] * fbe_0 - g_0_0_xxxyyyy_1[i] * fz_be_0 + 4.0 * g_y_0_xxxyyy_1[i] * fi_acd_0 + g_y_0_xxxyyyy_1[i] * wa_y[i];

        g_yy_0_xxxyyyz_0[i] = g_0_0_xxxyyyz_0[i] * fbe_0 - g_0_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_y_0_xxxyyz_1[i] * fi_acd_0 + g_y_0_xxxyyyz_1[i] * wa_y[i];

        g_yy_0_xxxyyzz_0[i] = g_0_0_xxxyyzz_0[i] * fbe_0 - g_0_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxyzz_1[i] * fi_acd_0 + g_y_0_xxxyyzz_1[i] * wa_y[i];

        g_yy_0_xxxyzzz_0[i] = g_0_0_xxxyzzz_0[i] * fbe_0 - g_0_0_xxxyzzz_1[i] * fz_be_0 + g_y_0_xxxzzz_1[i] * fi_acd_0 + g_y_0_xxxyzzz_1[i] * wa_y[i];

        g_yy_0_xxxzzzz_0[i] = g_0_0_xxxzzzz_0[i] * fbe_0 - g_0_0_xxxzzzz_1[i] * fz_be_0 + g_y_0_xxxzzzz_1[i] * wa_y[i];

        g_yy_0_xxyyyyy_0[i] = g_0_0_xxyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyy_1[i] * fz_be_0 + 5.0 * g_y_0_xxyyyy_1[i] * fi_acd_0 + g_y_0_xxyyyyy_1[i] * wa_y[i];

        g_yy_0_xxyyyyz_0[i] = g_0_0_xxyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_y_0_xxyyyz_1[i] * fi_acd_0 + g_y_0_xxyyyyz_1[i] * wa_y[i];

        g_yy_0_xxyyyzz_0[i] = g_0_0_xxyyyzz_0[i] * fbe_0 - g_0_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_y_0_xxyyzz_1[i] * fi_acd_0 + g_y_0_xxyyyzz_1[i] * wa_y[i];

        g_yy_0_xxyyzzz_0[i] = g_0_0_xxyyzzz_0[i] * fbe_0 - g_0_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxyzzz_1[i] * fi_acd_0 + g_y_0_xxyyzzz_1[i] * wa_y[i];

        g_yy_0_xxyzzzz_0[i] = g_0_0_xxyzzzz_0[i] * fbe_0 - g_0_0_xxyzzzz_1[i] * fz_be_0 + g_y_0_xxzzzz_1[i] * fi_acd_0 + g_y_0_xxyzzzz_1[i] * wa_y[i];

        g_yy_0_xxzzzzz_0[i] = g_0_0_xxzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzz_1[i] * fz_be_0 + g_y_0_xxzzzzz_1[i] * wa_y[i];

        g_yy_0_xyyyyyy_0[i] = g_0_0_xyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyy_1[i] * fz_be_0 + 6.0 * g_y_0_xyyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyyy_1[i] * wa_y[i];

        g_yy_0_xyyyyyz_0[i] = g_0_0_xyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_y_0_xyyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyyz_1[i] * wa_y[i];

        g_yy_0_xyyyyzz_0[i] = g_0_0_xyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_y_0_xyyyzz_1[i] * fi_acd_0 + g_y_0_xyyyyzz_1[i] * wa_y[i];

        g_yy_0_xyyyzzz_0[i] = g_0_0_xyyyzzz_0[i] * fbe_0 - g_0_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_y_0_xyyzzz_1[i] * fi_acd_0 + g_y_0_xyyyzzz_1[i] * wa_y[i];

        g_yy_0_xyyzzzz_0[i] = g_0_0_xyyzzzz_0[i] * fbe_0 - g_0_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xyzzzz_1[i] * fi_acd_0 + g_y_0_xyyzzzz_1[i] * wa_y[i];

        g_yy_0_xyzzzzz_0[i] = g_0_0_xyzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzz_1[i] * fz_be_0 + g_y_0_xzzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzzz_1[i] * wa_y[i];

        g_yy_0_xzzzzzz_0[i] = g_0_0_xzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzz_1[i] * fz_be_0 + g_y_0_xzzzzzz_1[i] * wa_y[i];

        g_yy_0_yyyyyyy_0[i] = g_0_0_yyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyy_1[i] * fz_be_0 + 7.0 * g_y_0_yyyyyy_1[i] * fi_acd_0 + g_y_0_yyyyyyy_1[i] * wa_y[i];

        g_yy_0_yyyyyyz_0[i] = g_0_0_yyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_y_0_yyyyyz_1[i] * fi_acd_0 + g_y_0_yyyyyyz_1[i] * wa_y[i];

        g_yy_0_yyyyyzz_0[i] = g_0_0_yyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_y_0_yyyyzz_1[i] * fi_acd_0 + g_y_0_yyyyyzz_1[i] * wa_y[i];

        g_yy_0_yyyyzzz_0[i] = g_0_0_yyyyzzz_0[i] * fbe_0 - g_0_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_y_0_yyyzzz_1[i] * fi_acd_0 + g_y_0_yyyyzzz_1[i] * wa_y[i];

        g_yy_0_yyyzzzz_0[i] = g_0_0_yyyzzzz_0[i] * fbe_0 - g_0_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_y_0_yyzzzz_1[i] * fi_acd_0 + g_y_0_yyyzzzz_1[i] * wa_y[i];

        g_yy_0_yyzzzzz_0[i] = g_0_0_yyzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_yzzzzz_1[i] * fi_acd_0 + g_y_0_yyzzzzz_1[i] * wa_y[i];

        g_yy_0_yzzzzzz_0[i] = g_0_0_yzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzz_1[i] * fz_be_0 + g_y_0_zzzzzz_1[i] * fi_acd_0 + g_y_0_yzzzzzz_1[i] * wa_y[i];

        g_yy_0_zzzzzzz_0[i] = g_0_0_zzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzz_1[i] * fz_be_0 + g_y_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 144-180 components of targeted buffer : DSK

    auto g_yz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk + 144);

    auto g_yz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 145);

    auto g_yz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 146);

    auto g_yz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 147);

    auto g_yz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 148);

    auto g_yz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 149);

    auto g_yz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 150);

    auto g_yz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 151);

    auto g_yz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 152);

    auto g_yz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 153);

    auto g_yz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 154);

    auto g_yz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 155);

    auto g_yz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 156);

    auto g_yz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 157);

    auto g_yz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 158);

    auto g_yz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 159);

    auto g_yz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 160);

    auto g_yz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 161);

    auto g_yz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 162);

    auto g_yz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 163);

    auto g_yz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 164);

    auto g_yz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 165);

    auto g_yz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 166);

    auto g_yz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 167);

    auto g_yz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 168);

    auto g_yz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 169);

    auto g_yz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 170);

    auto g_yz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 171);

    auto g_yz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 172);

    auto g_yz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 173);

    auto g_yz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 174);

    auto g_yz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 175);

    auto g_yz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 176);

    auto g_yz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 177);

    auto g_yz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 178);

    auto g_yz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 179);

    #pragma omp simd aligned(g_y_0_xxxxxxy_1, g_y_0_xxxxxyy_1, g_y_0_xxxxyyy_1, g_y_0_xxxyyyy_1, g_y_0_xxyyyyy_1, g_y_0_xyyyyyy_1, g_y_0_yyyyyyy_1, g_yz_0_xxxxxxx_0, g_yz_0_xxxxxxy_0, g_yz_0_xxxxxxz_0, g_yz_0_xxxxxyy_0, g_yz_0_xxxxxyz_0, g_yz_0_xxxxxzz_0, g_yz_0_xxxxyyy_0, g_yz_0_xxxxyyz_0, g_yz_0_xxxxyzz_0, g_yz_0_xxxxzzz_0, g_yz_0_xxxyyyy_0, g_yz_0_xxxyyyz_0, g_yz_0_xxxyyzz_0, g_yz_0_xxxyzzz_0, g_yz_0_xxxzzzz_0, g_yz_0_xxyyyyy_0, g_yz_0_xxyyyyz_0, g_yz_0_xxyyyzz_0, g_yz_0_xxyyzzz_0, g_yz_0_xxyzzzz_0, g_yz_0_xxzzzzz_0, g_yz_0_xyyyyyy_0, g_yz_0_xyyyyyz_0, g_yz_0_xyyyyzz_0, g_yz_0_xyyyzzz_0, g_yz_0_xyyzzzz_0, g_yz_0_xyzzzzz_0, g_yz_0_xzzzzzz_0, g_yz_0_yyyyyyy_0, g_yz_0_yyyyyyz_0, g_yz_0_yyyyyzz_0, g_yz_0_yyyyzzz_0, g_yz_0_yyyzzzz_0, g_yz_0_yyzzzzz_0, g_yz_0_yzzzzzz_0, g_yz_0_zzzzzzz_0, g_z_0_xxxxxxx_1, g_z_0_xxxxxxz_1, g_z_0_xxxxxyz_1, g_z_0_xxxxxz_1, g_z_0_xxxxxzz_1, g_z_0_xxxxyyz_1, g_z_0_xxxxyz_1, g_z_0_xxxxyzz_1, g_z_0_xxxxzz_1, g_z_0_xxxxzzz_1, g_z_0_xxxyyyz_1, g_z_0_xxxyyz_1, g_z_0_xxxyyzz_1, g_z_0_xxxyzz_1, g_z_0_xxxyzzz_1, g_z_0_xxxzzz_1, g_z_0_xxxzzzz_1, g_z_0_xxyyyyz_1, g_z_0_xxyyyz_1, g_z_0_xxyyyzz_1, g_z_0_xxyyzz_1, g_z_0_xxyyzzz_1, g_z_0_xxyzzz_1, g_z_0_xxyzzzz_1, g_z_0_xxzzzz_1, g_z_0_xxzzzzz_1, g_z_0_xyyyyyz_1, g_z_0_xyyyyz_1, g_z_0_xyyyyzz_1, g_z_0_xyyyzz_1, g_z_0_xyyyzzz_1, g_z_0_xyyzzz_1, g_z_0_xyyzzzz_1, g_z_0_xyzzzz_1, g_z_0_xyzzzzz_1, g_z_0_xzzzzz_1, g_z_0_xzzzzzz_1, g_z_0_yyyyyyz_1, g_z_0_yyyyyz_1, g_z_0_yyyyyzz_1, g_z_0_yyyyzz_1, g_z_0_yyyyzzz_1, g_z_0_yyyzzz_1, g_z_0_yyyzzzz_1, g_z_0_yyzzzz_1, g_z_0_yyzzzzz_1, g_z_0_yzzzzz_1, g_z_0_yzzzzzz_1, g_z_0_zzzzzz_1, g_z_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xxxxxxx_0[i] = g_z_0_xxxxxxx_1[i] * wa_y[i];

        g_yz_0_xxxxxxy_0[i] = g_y_0_xxxxxxy_1[i] * wa_z[i];

        g_yz_0_xxxxxxz_0[i] = g_z_0_xxxxxxz_1[i] * wa_y[i];

        g_yz_0_xxxxxyy_0[i] = g_y_0_xxxxxyy_1[i] * wa_z[i];

        g_yz_0_xxxxxyz_0[i] = g_z_0_xxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxyz_1[i] * wa_y[i];

        g_yz_0_xxxxxzz_0[i] = g_z_0_xxxxxzz_1[i] * wa_y[i];

        g_yz_0_xxxxyyy_0[i] = g_y_0_xxxxyyy_1[i] * wa_z[i];

        g_yz_0_xxxxyyz_0[i] = 2.0 * g_z_0_xxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxyyz_1[i] * wa_y[i];

        g_yz_0_xxxxyzz_0[i] = g_z_0_xxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxyzz_1[i] * wa_y[i];

        g_yz_0_xxxxzzz_0[i] = g_z_0_xxxxzzz_1[i] * wa_y[i];

        g_yz_0_xxxyyyy_0[i] = g_y_0_xxxyyyy_1[i] * wa_z[i];

        g_yz_0_xxxyyyz_0[i] = 3.0 * g_z_0_xxxyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyz_1[i] * wa_y[i];

        g_yz_0_xxxyyzz_0[i] = 2.0 * g_z_0_xxxyzz_1[i] * fi_acd_0 + g_z_0_xxxyyzz_1[i] * wa_y[i];

        g_yz_0_xxxyzzz_0[i] = g_z_0_xxxzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzz_1[i] * wa_y[i];

        g_yz_0_xxxzzzz_0[i] = g_z_0_xxxzzzz_1[i] * wa_y[i];

        g_yz_0_xxyyyyy_0[i] = g_y_0_xxyyyyy_1[i] * wa_z[i];

        g_yz_0_xxyyyyz_0[i] = 4.0 * g_z_0_xxyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyz_1[i] * wa_y[i];

        g_yz_0_xxyyyzz_0[i] = 3.0 * g_z_0_xxyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyzz_1[i] * wa_y[i];

        g_yz_0_xxyyzzz_0[i] = 2.0 * g_z_0_xxyzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzz_1[i] * wa_y[i];

        g_yz_0_xxyzzzz_0[i] = g_z_0_xxzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzz_1[i] * wa_y[i];

        g_yz_0_xxzzzzz_0[i] = g_z_0_xxzzzzz_1[i] * wa_y[i];

        g_yz_0_xyyyyyy_0[i] = g_y_0_xyyyyyy_1[i] * wa_z[i];

        g_yz_0_xyyyyyz_0[i] = 5.0 * g_z_0_xyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyz_1[i] * wa_y[i];

        g_yz_0_xyyyyzz_0[i] = 4.0 * g_z_0_xyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyzz_1[i] * wa_y[i];

        g_yz_0_xyyyzzz_0[i] = 3.0 * g_z_0_xyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzz_1[i] * wa_y[i];

        g_yz_0_xyyzzzz_0[i] = 2.0 * g_z_0_xyzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzz_1[i] * wa_y[i];

        g_yz_0_xyzzzzz_0[i] = g_z_0_xzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzz_1[i] * wa_y[i];

        g_yz_0_xzzzzzz_0[i] = g_z_0_xzzzzzz_1[i] * wa_y[i];

        g_yz_0_yyyyyyy_0[i] = g_y_0_yyyyyyy_1[i] * wa_z[i];

        g_yz_0_yyyyyyz_0[i] = 6.0 * g_z_0_yyyyyz_1[i] * fi_acd_0 + g_z_0_yyyyyyz_1[i] * wa_y[i];

        g_yz_0_yyyyyzz_0[i] = 5.0 * g_z_0_yyyyzz_1[i] * fi_acd_0 + g_z_0_yyyyyzz_1[i] * wa_y[i];

        g_yz_0_yyyyzzz_0[i] = 4.0 * g_z_0_yyyzzz_1[i] * fi_acd_0 + g_z_0_yyyyzzz_1[i] * wa_y[i];

        g_yz_0_yyyzzzz_0[i] = 3.0 * g_z_0_yyzzzz_1[i] * fi_acd_0 + g_z_0_yyyzzzz_1[i] * wa_y[i];

        g_yz_0_yyzzzzz_0[i] = 2.0 * g_z_0_yzzzzz_1[i] * fi_acd_0 + g_z_0_yyzzzzz_1[i] * wa_y[i];

        g_yz_0_yzzzzzz_0[i] = g_z_0_zzzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzzz_1[i] * wa_y[i];

        g_yz_0_zzzzzzz_0[i] = g_z_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 180-216 components of targeted buffer : DSK

    auto g_zz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk + 180);

    auto g_zz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 181);

    auto g_zz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 182);

    auto g_zz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 183);

    auto g_zz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 184);

    auto g_zz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 185);

    auto g_zz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 186);

    auto g_zz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 187);

    auto g_zz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 188);

    auto g_zz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 189);

    auto g_zz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 190);

    auto g_zz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 191);

    auto g_zz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 192);

    auto g_zz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 193);

    auto g_zz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 194);

    auto g_zz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 195);

    auto g_zz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 196);

    auto g_zz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 197);

    auto g_zz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 198);

    auto g_zz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 199);

    auto g_zz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 200);

    auto g_zz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 201);

    auto g_zz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 202);

    auto g_zz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 203);

    auto g_zz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 204);

    auto g_zz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 205);

    auto g_zz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 206);

    auto g_zz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 207);

    auto g_zz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 208);

    auto g_zz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 209);

    auto g_zz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 210);

    auto g_zz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 211);

    auto g_zz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 212);

    auto g_zz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 213);

    auto g_zz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 214);

    auto g_zz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 215);

    #pragma omp simd aligned(g_0_0_xxxxxxx_0, g_0_0_xxxxxxx_1, g_0_0_xxxxxxy_0, g_0_0_xxxxxxy_1, g_0_0_xxxxxxz_0, g_0_0_xxxxxxz_1, g_0_0_xxxxxyy_0, g_0_0_xxxxxyy_1, g_0_0_xxxxxyz_0, g_0_0_xxxxxyz_1, g_0_0_xxxxxzz_0, g_0_0_xxxxxzz_1, g_0_0_xxxxyyy_0, g_0_0_xxxxyyy_1, g_0_0_xxxxyyz_0, g_0_0_xxxxyyz_1, g_0_0_xxxxyzz_0, g_0_0_xxxxyzz_1, g_0_0_xxxxzzz_0, g_0_0_xxxxzzz_1, g_0_0_xxxyyyy_0, g_0_0_xxxyyyy_1, g_0_0_xxxyyyz_0, g_0_0_xxxyyyz_1, g_0_0_xxxyyzz_0, g_0_0_xxxyyzz_1, g_0_0_xxxyzzz_0, g_0_0_xxxyzzz_1, g_0_0_xxxzzzz_0, g_0_0_xxxzzzz_1, g_0_0_xxyyyyy_0, g_0_0_xxyyyyy_1, g_0_0_xxyyyyz_0, g_0_0_xxyyyyz_1, g_0_0_xxyyyzz_0, g_0_0_xxyyyzz_1, g_0_0_xxyyzzz_0, g_0_0_xxyyzzz_1, g_0_0_xxyzzzz_0, g_0_0_xxyzzzz_1, g_0_0_xxzzzzz_0, g_0_0_xxzzzzz_1, g_0_0_xyyyyyy_0, g_0_0_xyyyyyy_1, g_0_0_xyyyyyz_0, g_0_0_xyyyyyz_1, g_0_0_xyyyyzz_0, g_0_0_xyyyyzz_1, g_0_0_xyyyzzz_0, g_0_0_xyyyzzz_1, g_0_0_xyyzzzz_0, g_0_0_xyyzzzz_1, g_0_0_xyzzzzz_0, g_0_0_xyzzzzz_1, g_0_0_xzzzzzz_0, g_0_0_xzzzzzz_1, g_0_0_yyyyyyy_0, g_0_0_yyyyyyy_1, g_0_0_yyyyyyz_0, g_0_0_yyyyyyz_1, g_0_0_yyyyyzz_0, g_0_0_yyyyyzz_1, g_0_0_yyyyzzz_0, g_0_0_yyyyzzz_1, g_0_0_yyyzzzz_0, g_0_0_yyyzzzz_1, g_0_0_yyzzzzz_0, g_0_0_yyzzzzz_1, g_0_0_yzzzzzz_0, g_0_0_yzzzzzz_1, g_0_0_zzzzzzz_0, g_0_0_zzzzzzz_1, g_z_0_xxxxxx_1, g_z_0_xxxxxxx_1, g_z_0_xxxxxxy_1, g_z_0_xxxxxxz_1, g_z_0_xxxxxy_1, g_z_0_xxxxxyy_1, g_z_0_xxxxxyz_1, g_z_0_xxxxxz_1, g_z_0_xxxxxzz_1, g_z_0_xxxxyy_1, g_z_0_xxxxyyy_1, g_z_0_xxxxyyz_1, g_z_0_xxxxyz_1, g_z_0_xxxxyzz_1, g_z_0_xxxxzz_1, g_z_0_xxxxzzz_1, g_z_0_xxxyyy_1, g_z_0_xxxyyyy_1, g_z_0_xxxyyyz_1, g_z_0_xxxyyz_1, g_z_0_xxxyyzz_1, g_z_0_xxxyzz_1, g_z_0_xxxyzzz_1, g_z_0_xxxzzz_1, g_z_0_xxxzzzz_1, g_z_0_xxyyyy_1, g_z_0_xxyyyyy_1, g_z_0_xxyyyyz_1, g_z_0_xxyyyz_1, g_z_0_xxyyyzz_1, g_z_0_xxyyzz_1, g_z_0_xxyyzzz_1, g_z_0_xxyzzz_1, g_z_0_xxyzzzz_1, g_z_0_xxzzzz_1, g_z_0_xxzzzzz_1, g_z_0_xyyyyy_1, g_z_0_xyyyyyy_1, g_z_0_xyyyyyz_1, g_z_0_xyyyyz_1, g_z_0_xyyyyzz_1, g_z_0_xyyyzz_1, g_z_0_xyyyzzz_1, g_z_0_xyyzzz_1, g_z_0_xyyzzzz_1, g_z_0_xyzzzz_1, g_z_0_xyzzzzz_1, g_z_0_xzzzzz_1, g_z_0_xzzzzzz_1, g_z_0_yyyyyy_1, g_z_0_yyyyyyy_1, g_z_0_yyyyyyz_1, g_z_0_yyyyyz_1, g_z_0_yyyyyzz_1, g_z_0_yyyyzz_1, g_z_0_yyyyzzz_1, g_z_0_yyyzzz_1, g_z_0_yyyzzzz_1, g_z_0_yyzzzz_1, g_z_0_yyzzzzz_1, g_z_0_yzzzzz_1, g_z_0_yzzzzzz_1, g_z_0_zzzzzz_1, g_z_0_zzzzzzz_1, g_zz_0_xxxxxxx_0, g_zz_0_xxxxxxy_0, g_zz_0_xxxxxxz_0, g_zz_0_xxxxxyy_0, g_zz_0_xxxxxyz_0, g_zz_0_xxxxxzz_0, g_zz_0_xxxxyyy_0, g_zz_0_xxxxyyz_0, g_zz_0_xxxxyzz_0, g_zz_0_xxxxzzz_0, g_zz_0_xxxyyyy_0, g_zz_0_xxxyyyz_0, g_zz_0_xxxyyzz_0, g_zz_0_xxxyzzz_0, g_zz_0_xxxzzzz_0, g_zz_0_xxyyyyy_0, g_zz_0_xxyyyyz_0, g_zz_0_xxyyyzz_0, g_zz_0_xxyyzzz_0, g_zz_0_xxyzzzz_0, g_zz_0_xxzzzzz_0, g_zz_0_xyyyyyy_0, g_zz_0_xyyyyyz_0, g_zz_0_xyyyyzz_0, g_zz_0_xyyyzzz_0, g_zz_0_xyyzzzz_0, g_zz_0_xyzzzzz_0, g_zz_0_xzzzzzz_0, g_zz_0_yyyyyyy_0, g_zz_0_yyyyyyz_0, g_zz_0_yyyyyzz_0, g_zz_0_yyyyzzz_0, g_zz_0_yyyzzzz_0, g_zz_0_yyzzzzz_0, g_zz_0_yzzzzzz_0, g_zz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xxxxxxx_0[i] = g_0_0_xxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxx_1[i] * fz_be_0 + g_z_0_xxxxxxx_1[i] * wa_z[i];

        g_zz_0_xxxxxxy_0[i] = g_0_0_xxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxy_1[i] * fz_be_0 + g_z_0_xxxxxxy_1[i] * wa_z[i];

        g_zz_0_xxxxxxz_0[i] = g_0_0_xxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxz_1[i] * fz_be_0 + g_z_0_xxxxxx_1[i] * fi_acd_0 + g_z_0_xxxxxxz_1[i] * wa_z[i];

        g_zz_0_xxxxxyy_0[i] = g_0_0_xxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxyy_1[i] * fz_be_0 + g_z_0_xxxxxyy_1[i] * wa_z[i];

        g_zz_0_xxxxxyz_0[i] = g_0_0_xxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxyz_1[i] * fz_be_0 + g_z_0_xxxxxy_1[i] * fi_acd_0 + g_z_0_xxxxxyz_1[i] * wa_z[i];

        g_zz_0_xxxxxzz_0[i] = g_0_0_xxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxzz_1[i] * wa_z[i];

        g_zz_0_xxxxyyy_0[i] = g_0_0_xxxxyyy_0[i] * fbe_0 - g_0_0_xxxxyyy_1[i] * fz_be_0 + g_z_0_xxxxyyy_1[i] * wa_z[i];

        g_zz_0_xxxxyyz_0[i] = g_0_0_xxxxyyz_0[i] * fbe_0 - g_0_0_xxxxyyz_1[i] * fz_be_0 + g_z_0_xxxxyy_1[i] * fi_acd_0 + g_z_0_xxxxyyz_1[i] * wa_z[i];

        g_zz_0_xxxxyzz_0[i] = g_0_0_xxxxyzz_0[i] * fbe_0 - g_0_0_xxxxyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxyzz_1[i] * wa_z[i];

        g_zz_0_xxxxzzz_0[i] = g_0_0_xxxxzzz_0[i] * fbe_0 - g_0_0_xxxxzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxzzz_1[i] * wa_z[i];

        g_zz_0_xxxyyyy_0[i] = g_0_0_xxxyyyy_0[i] * fbe_0 - g_0_0_xxxyyyy_1[i] * fz_be_0 + g_z_0_xxxyyyy_1[i] * wa_z[i];

        g_zz_0_xxxyyyz_0[i] = g_0_0_xxxyyyz_0[i] * fbe_0 - g_0_0_xxxyyyz_1[i] * fz_be_0 + g_z_0_xxxyyy_1[i] * fi_acd_0 + g_z_0_xxxyyyz_1[i] * wa_z[i];

        g_zz_0_xxxyyzz_0[i] = g_0_0_xxxyyzz_0[i] * fbe_0 - g_0_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxyyz_1[i] * fi_acd_0 + g_z_0_xxxyyzz_1[i] * wa_z[i];

        g_zz_0_xxxyzzz_0[i] = g_0_0_xxxyzzz_0[i] * fbe_0 - g_0_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxyzz_1[i] * fi_acd_0 + g_z_0_xxxyzzz_1[i] * wa_z[i];

        g_zz_0_xxxzzzz_0[i] = g_0_0_xxxzzzz_0[i] * fbe_0 - g_0_0_xxxzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxxzzz_1[i] * fi_acd_0 + g_z_0_xxxzzzz_1[i] * wa_z[i];

        g_zz_0_xxyyyyy_0[i] = g_0_0_xxyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyy_1[i] * fz_be_0 + g_z_0_xxyyyyy_1[i] * wa_z[i];

        g_zz_0_xxyyyyz_0[i] = g_0_0_xxyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyz_1[i] * fz_be_0 + g_z_0_xxyyyy_1[i] * fi_acd_0 + g_z_0_xxyyyyz_1[i] * wa_z[i];

        g_zz_0_xxyyyzz_0[i] = g_0_0_xxyyyzz_0[i] * fbe_0 - g_0_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyzz_1[i] * wa_z[i];

        g_zz_0_xxyyzzz_0[i] = g_0_0_xxyyzzz_0[i] * fbe_0 - g_0_0_xxyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxyyzz_1[i] * fi_acd_0 + g_z_0_xxyyzzz_1[i] * wa_z[i];

        g_zz_0_xxyzzzz_0[i] = g_0_0_xxyzzzz_0[i] * fbe_0 - g_0_0_xxyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxyzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzz_1[i] * wa_z[i];

        g_zz_0_xxzzzzz_0[i] = g_0_0_xxzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xxzzzz_1[i] * fi_acd_0 + g_z_0_xxzzzzz_1[i] * wa_z[i];

        g_zz_0_xyyyyyy_0[i] = g_0_0_xyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyy_1[i] * fz_be_0 + g_z_0_xyyyyyy_1[i] * wa_z[i];

        g_zz_0_xyyyyyz_0[i] = g_0_0_xyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyz_1[i] * fz_be_0 + g_z_0_xyyyyy_1[i] * fi_acd_0 + g_z_0_xyyyyyz_1[i] * wa_z[i];

        g_zz_0_xyyyyzz_0[i] = g_0_0_xyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyzz_1[i] * wa_z[i];

        g_zz_0_xyyyzzz_0[i] = g_0_0_xyyyzzz_0[i] * fbe_0 - g_0_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyzzz_1[i] * wa_z[i];

        g_zz_0_xyyzzzz_0[i] = g_0_0_xyyzzzz_0[i] * fbe_0 - g_0_0_xyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xyyzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzz_1[i] * wa_z[i];

        g_zz_0_xyzzzzz_0[i] = g_0_0_xyzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xyzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzz_1[i] * wa_z[i];

        g_zz_0_xzzzzzz_0[i] = g_0_0_xzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_xzzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzzz_1[i] * wa_z[i];

        g_zz_0_yyyyyyy_0[i] = g_0_0_yyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyy_1[i] * fz_be_0 + g_z_0_yyyyyyy_1[i] * wa_z[i];

        g_zz_0_yyyyyyz_0[i] = g_0_0_yyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyz_1[i] * fz_be_0 + g_z_0_yyyyyy_1[i] * fi_acd_0 + g_z_0_yyyyyyz_1[i] * wa_z[i];

        g_zz_0_yyyyyzz_0[i] = g_0_0_yyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_yyyyyz_1[i] * fi_acd_0 + g_z_0_yyyyyzz_1[i] * wa_z[i];

        g_zz_0_yyyyzzz_0[i] = g_0_0_yyyyzzz_0[i] * fbe_0 - g_0_0_yyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_yyyyzz_1[i] * fi_acd_0 + g_z_0_yyyyzzz_1[i] * wa_z[i];

        g_zz_0_yyyzzzz_0[i] = g_0_0_yyyzzzz_0[i] * fbe_0 - g_0_0_yyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_yyyzzz_1[i] * fi_acd_0 + g_z_0_yyyzzzz_1[i] * wa_z[i];

        g_zz_0_yyzzzzz_0[i] = g_0_0_yyzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_yyzzzz_1[i] * fi_acd_0 + g_z_0_yyzzzzz_1[i] * wa_z[i];

        g_zz_0_yzzzzzz_0[i] = g_0_0_yzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_yzzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzzz_1[i] * wa_z[i];

        g_zz_0_zzzzzzz_0[i] = g_0_0_zzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzz_1[i] * fz_be_0 + 7.0 * g_z_0_zzzzzz_1[i] * fi_acd_0 + g_z_0_zzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

