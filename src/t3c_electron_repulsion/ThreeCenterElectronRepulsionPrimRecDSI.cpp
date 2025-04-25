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

#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsi,
                                 size_t idx_eri_0_ssi,
                                 size_t idx_eri_1_ssi,
                                 size_t idx_eri_1_psh,
                                 size_t idx_eri_1_psi,
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

    /// Set up components of auxilary buffer : SSI

    auto g_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ssi);

    auto g_0_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ssi + 1);

    auto g_0_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ssi + 2);

    auto g_0_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ssi + 3);

    auto g_0_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ssi + 4);

    auto g_0_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ssi + 5);

    auto g_0_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ssi + 6);

    auto g_0_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ssi + 7);

    auto g_0_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ssi + 8);

    auto g_0_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ssi + 9);

    auto g_0_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ssi + 10);

    auto g_0_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ssi + 11);

    auto g_0_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ssi + 12);

    auto g_0_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ssi + 13);

    auto g_0_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ssi + 14);

    auto g_0_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ssi + 15);

    auto g_0_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ssi + 16);

    auto g_0_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ssi + 17);

    auto g_0_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ssi + 18);

    auto g_0_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ssi + 19);

    auto g_0_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 20);

    auto g_0_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ssi + 21);

    auto g_0_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ssi + 22);

    auto g_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ssi + 23);

    auto g_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ssi + 24);

    auto g_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ssi + 25);

    auto g_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 26);

    auto g_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 27);

    /// Set up components of auxilary buffer : SSI

    auto g_0_0_xxxxxx_1 = pbuffer.data(idx_eri_1_ssi);

    auto g_0_0_xxxxxy_1 = pbuffer.data(idx_eri_1_ssi + 1);

    auto g_0_0_xxxxxz_1 = pbuffer.data(idx_eri_1_ssi + 2);

    auto g_0_0_xxxxyy_1 = pbuffer.data(idx_eri_1_ssi + 3);

    auto g_0_0_xxxxyz_1 = pbuffer.data(idx_eri_1_ssi + 4);

    auto g_0_0_xxxxzz_1 = pbuffer.data(idx_eri_1_ssi + 5);

    auto g_0_0_xxxyyy_1 = pbuffer.data(idx_eri_1_ssi + 6);

    auto g_0_0_xxxyyz_1 = pbuffer.data(idx_eri_1_ssi + 7);

    auto g_0_0_xxxyzz_1 = pbuffer.data(idx_eri_1_ssi + 8);

    auto g_0_0_xxxzzz_1 = pbuffer.data(idx_eri_1_ssi + 9);

    auto g_0_0_xxyyyy_1 = pbuffer.data(idx_eri_1_ssi + 10);

    auto g_0_0_xxyyyz_1 = pbuffer.data(idx_eri_1_ssi + 11);

    auto g_0_0_xxyyzz_1 = pbuffer.data(idx_eri_1_ssi + 12);

    auto g_0_0_xxyzzz_1 = pbuffer.data(idx_eri_1_ssi + 13);

    auto g_0_0_xxzzzz_1 = pbuffer.data(idx_eri_1_ssi + 14);

    auto g_0_0_xyyyyy_1 = pbuffer.data(idx_eri_1_ssi + 15);

    auto g_0_0_xyyyyz_1 = pbuffer.data(idx_eri_1_ssi + 16);

    auto g_0_0_xyyyzz_1 = pbuffer.data(idx_eri_1_ssi + 17);

    auto g_0_0_xyyzzz_1 = pbuffer.data(idx_eri_1_ssi + 18);

    auto g_0_0_xyzzzz_1 = pbuffer.data(idx_eri_1_ssi + 19);

    auto g_0_0_xzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 20);

    auto g_0_0_yyyyyy_1 = pbuffer.data(idx_eri_1_ssi + 21);

    auto g_0_0_yyyyyz_1 = pbuffer.data(idx_eri_1_ssi + 22);

    auto g_0_0_yyyyzz_1 = pbuffer.data(idx_eri_1_ssi + 23);

    auto g_0_0_yyyzzz_1 = pbuffer.data(idx_eri_1_ssi + 24);

    auto g_0_0_yyzzzz_1 = pbuffer.data(idx_eri_1_ssi + 25);

    auto g_0_0_yzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 26);

    auto g_0_0_zzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 27);

    /// Set up components of auxilary buffer : PSH

    auto g_x_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh);

    auto g_x_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 1);

    auto g_x_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 2);

    auto g_x_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 3);

    auto g_x_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 4);

    auto g_x_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 5);

    auto g_x_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 6);

    auto g_x_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 7);

    auto g_x_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 8);

    auto g_x_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 9);

    auto g_x_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 10);

    auto g_x_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 11);

    auto g_x_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 12);

    auto g_x_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 13);

    auto g_x_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 14);

    auto g_x_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 15);

    auto g_x_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 16);

    auto g_x_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 17);

    auto g_x_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 18);

    auto g_x_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 19);

    auto g_x_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 20);

    auto g_y_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh + 21);

    auto g_y_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 22);

    auto g_y_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 23);

    auto g_y_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 24);

    auto g_y_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 25);

    auto g_y_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 26);

    auto g_y_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 27);

    auto g_y_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 28);

    auto g_y_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 29);

    auto g_y_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 30);

    auto g_y_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 31);

    auto g_y_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 32);

    auto g_y_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 33);

    auto g_y_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 34);

    auto g_y_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 35);

    auto g_y_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 36);

    auto g_y_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 37);

    auto g_y_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 38);

    auto g_y_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 39);

    auto g_y_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 40);

    auto g_y_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 41);

    auto g_z_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh + 42);

    auto g_z_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 43);

    auto g_z_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 44);

    auto g_z_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 45);

    auto g_z_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 46);

    auto g_z_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 47);

    auto g_z_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 48);

    auto g_z_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 49);

    auto g_z_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 50);

    auto g_z_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 51);

    auto g_z_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 52);

    auto g_z_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 53);

    auto g_z_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 54);

    auto g_z_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 55);

    auto g_z_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 56);

    auto g_z_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 57);

    auto g_z_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 58);

    auto g_z_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 59);

    auto g_z_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 60);

    auto g_z_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 61);

    auto g_z_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 62);

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

    /// Set up 0-28 components of targeted buffer : DSI

    auto g_xx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi);

    auto g_xx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 1);

    auto g_xx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 2);

    auto g_xx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 3);

    auto g_xx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 4);

    auto g_xx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 5);

    auto g_xx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 6);

    auto g_xx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 7);

    auto g_xx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 8);

    auto g_xx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 9);

    auto g_xx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 10);

    auto g_xx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 11);

    auto g_xx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 12);

    auto g_xx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 13);

    auto g_xx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 14);

    auto g_xx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 15);

    auto g_xx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 16);

    auto g_xx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 17);

    auto g_xx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 18);

    auto g_xx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 19);

    auto g_xx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 20);

    auto g_xx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 21);

    auto g_xx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 22);

    auto g_xx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 23);

    auto g_xx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 24);

    auto g_xx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 25);

    auto g_xx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 26);

    auto g_xx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 27);

    #pragma omp simd aligned(g_0_0_xxxxxx_0, g_0_0_xxxxxx_1, g_0_0_xxxxxy_0, g_0_0_xxxxxy_1, g_0_0_xxxxxz_0, g_0_0_xxxxxz_1, g_0_0_xxxxyy_0, g_0_0_xxxxyy_1, g_0_0_xxxxyz_0, g_0_0_xxxxyz_1, g_0_0_xxxxzz_0, g_0_0_xxxxzz_1, g_0_0_xxxyyy_0, g_0_0_xxxyyy_1, g_0_0_xxxyyz_0, g_0_0_xxxyyz_1, g_0_0_xxxyzz_0, g_0_0_xxxyzz_1, g_0_0_xxxzzz_0, g_0_0_xxxzzz_1, g_0_0_xxyyyy_0, g_0_0_xxyyyy_1, g_0_0_xxyyyz_0, g_0_0_xxyyyz_1, g_0_0_xxyyzz_0, g_0_0_xxyyzz_1, g_0_0_xxyzzz_0, g_0_0_xxyzzz_1, g_0_0_xxzzzz_0, g_0_0_xxzzzz_1, g_0_0_xyyyyy_0, g_0_0_xyyyyy_1, g_0_0_xyyyyz_0, g_0_0_xyyyyz_1, g_0_0_xyyyzz_0, g_0_0_xyyyzz_1, g_0_0_xyyzzz_0, g_0_0_xyyzzz_1, g_0_0_xyzzzz_0, g_0_0_xyzzzz_1, g_0_0_xzzzzz_0, g_0_0_xzzzzz_1, g_0_0_yyyyyy_0, g_0_0_yyyyyy_1, g_0_0_yyyyyz_0, g_0_0_yyyyyz_1, g_0_0_yyyyzz_0, g_0_0_yyyyzz_1, g_0_0_yyyzzz_0, g_0_0_yyyzzz_1, g_0_0_yyzzzz_0, g_0_0_yyzzzz_1, g_0_0_yzzzzz_0, g_0_0_yzzzzz_1, g_0_0_zzzzzz_0, g_0_0_zzzzzz_1, g_x_0_xxxxx_1, g_x_0_xxxxxx_1, g_x_0_xxxxxy_1, g_x_0_xxxxxz_1, g_x_0_xxxxy_1, g_x_0_xxxxyy_1, g_x_0_xxxxyz_1, g_x_0_xxxxz_1, g_x_0_xxxxzz_1, g_x_0_xxxyy_1, g_x_0_xxxyyy_1, g_x_0_xxxyyz_1, g_x_0_xxxyz_1, g_x_0_xxxyzz_1, g_x_0_xxxzz_1, g_x_0_xxxzzz_1, g_x_0_xxyyy_1, g_x_0_xxyyyy_1, g_x_0_xxyyyz_1, g_x_0_xxyyz_1, g_x_0_xxyyzz_1, g_x_0_xxyzz_1, g_x_0_xxyzzz_1, g_x_0_xxzzz_1, g_x_0_xxzzzz_1, g_x_0_xyyyy_1, g_x_0_xyyyyy_1, g_x_0_xyyyyz_1, g_x_0_xyyyz_1, g_x_0_xyyyzz_1, g_x_0_xyyzz_1, g_x_0_xyyzzz_1, g_x_0_xyzzz_1, g_x_0_xyzzzz_1, g_x_0_xzzzz_1, g_x_0_xzzzzz_1, g_x_0_yyyyy_1, g_x_0_yyyyyy_1, g_x_0_yyyyyz_1, g_x_0_yyyyz_1, g_x_0_yyyyzz_1, g_x_0_yyyzz_1, g_x_0_yyyzzz_1, g_x_0_yyzzz_1, g_x_0_yyzzzz_1, g_x_0_yzzzz_1, g_x_0_yzzzzz_1, g_x_0_zzzzz_1, g_x_0_zzzzzz_1, g_xx_0_xxxxxx_0, g_xx_0_xxxxxy_0, g_xx_0_xxxxxz_0, g_xx_0_xxxxyy_0, g_xx_0_xxxxyz_0, g_xx_0_xxxxzz_0, g_xx_0_xxxyyy_0, g_xx_0_xxxyyz_0, g_xx_0_xxxyzz_0, g_xx_0_xxxzzz_0, g_xx_0_xxyyyy_0, g_xx_0_xxyyyz_0, g_xx_0_xxyyzz_0, g_xx_0_xxyzzz_0, g_xx_0_xxzzzz_0, g_xx_0_xyyyyy_0, g_xx_0_xyyyyz_0, g_xx_0_xyyyzz_0, g_xx_0_xyyzzz_0, g_xx_0_xyzzzz_0, g_xx_0_xzzzzz_0, g_xx_0_yyyyyy_0, g_xx_0_yyyyyz_0, g_xx_0_yyyyzz_0, g_xx_0_yyyzzz_0, g_xx_0_yyzzzz_0, g_xx_0_yzzzzz_0, g_xx_0_zzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xxxxxx_0[i] = g_0_0_xxxxxx_0[i] * fbe_0 - g_0_0_xxxxxx_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxx_1[i] * fi_acd_0 + g_x_0_xxxxxx_1[i] * wa_x[i];

        g_xx_0_xxxxxy_0[i] = g_0_0_xxxxxy_0[i] * fbe_0 - g_0_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxy_1[i] * fi_acd_0 + g_x_0_xxxxxy_1[i] * wa_x[i];

        g_xx_0_xxxxxz_0[i] = g_0_0_xxxxxz_0[i] * fbe_0 - g_0_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxz_1[i] * fi_acd_0 + g_x_0_xxxxxz_1[i] * wa_x[i];

        g_xx_0_xxxxyy_0[i] = g_0_0_xxxxyy_0[i] * fbe_0 - g_0_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyy_1[i] * fi_acd_0 + g_x_0_xxxxyy_1[i] * wa_x[i];

        g_xx_0_xxxxyz_0[i] = g_0_0_xxxxyz_0[i] * fbe_0 - g_0_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyz_1[i] * fi_acd_0 + g_x_0_xxxxyz_1[i] * wa_x[i];

        g_xx_0_xxxxzz_0[i] = g_0_0_xxxxzz_0[i] * fbe_0 - g_0_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxzz_1[i] * fi_acd_0 + g_x_0_xxxxzz_1[i] * wa_x[i];

        g_xx_0_xxxyyy_0[i] = g_0_0_xxxyyy_0[i] * fbe_0 - g_0_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyy_1[i] * fi_acd_0 + g_x_0_xxxyyy_1[i] * wa_x[i];

        g_xx_0_xxxyyz_0[i] = g_0_0_xxxyyz_0[i] * fbe_0 - g_0_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyz_1[i] * fi_acd_0 + g_x_0_xxxyyz_1[i] * wa_x[i];

        g_xx_0_xxxyzz_0[i] = g_0_0_xxxyzz_0[i] * fbe_0 - g_0_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyzz_1[i] * fi_acd_0 + g_x_0_xxxyzz_1[i] * wa_x[i];

        g_xx_0_xxxzzz_0[i] = g_0_0_xxxzzz_0[i] * fbe_0 - g_0_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxzzz_1[i] * fi_acd_0 + g_x_0_xxxzzz_1[i] * wa_x[i];

        g_xx_0_xxyyyy_0[i] = g_0_0_xxyyyy_0[i] * fbe_0 - g_0_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyy_1[i] * fi_acd_0 + g_x_0_xxyyyy_1[i] * wa_x[i];

        g_xx_0_xxyyyz_0[i] = g_0_0_xxyyyz_0[i] * fbe_0 - g_0_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyz_1[i] * fi_acd_0 + g_x_0_xxyyyz_1[i] * wa_x[i];

        g_xx_0_xxyyzz_0[i] = g_0_0_xxyyzz_0[i] * fbe_0 - g_0_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyzz_1[i] * fi_acd_0 + g_x_0_xxyyzz_1[i] * wa_x[i];

        g_xx_0_xxyzzz_0[i] = g_0_0_xxyzzz_0[i] * fbe_0 - g_0_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyzzz_1[i] * fi_acd_0 + g_x_0_xxyzzz_1[i] * wa_x[i];

        g_xx_0_xxzzzz_0[i] = g_0_0_xxzzzz_0[i] * fbe_0 - g_0_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xzzzz_1[i] * fi_acd_0 + g_x_0_xxzzzz_1[i] * wa_x[i];

        g_xx_0_xyyyyy_0[i] = g_0_0_xyyyyy_0[i] * fbe_0 - g_0_0_xyyyyy_1[i] * fz_be_0 + g_x_0_yyyyy_1[i] * fi_acd_0 + g_x_0_xyyyyy_1[i] * wa_x[i];

        g_xx_0_xyyyyz_0[i] = g_0_0_xyyyyz_0[i] * fbe_0 - g_0_0_xyyyyz_1[i] * fz_be_0 + g_x_0_yyyyz_1[i] * fi_acd_0 + g_x_0_xyyyyz_1[i] * wa_x[i];

        g_xx_0_xyyyzz_0[i] = g_0_0_xyyyzz_0[i] * fbe_0 - g_0_0_xyyyzz_1[i] * fz_be_0 + g_x_0_yyyzz_1[i] * fi_acd_0 + g_x_0_xyyyzz_1[i] * wa_x[i];

        g_xx_0_xyyzzz_0[i] = g_0_0_xyyzzz_0[i] * fbe_0 - g_0_0_xyyzzz_1[i] * fz_be_0 + g_x_0_yyzzz_1[i] * fi_acd_0 + g_x_0_xyyzzz_1[i] * wa_x[i];

        g_xx_0_xyzzzz_0[i] = g_0_0_xyzzzz_0[i] * fbe_0 - g_0_0_xyzzzz_1[i] * fz_be_0 + g_x_0_yzzzz_1[i] * fi_acd_0 + g_x_0_xyzzzz_1[i] * wa_x[i];

        g_xx_0_xzzzzz_0[i] = g_0_0_xzzzzz_0[i] * fbe_0 - g_0_0_xzzzzz_1[i] * fz_be_0 + g_x_0_zzzzz_1[i] * fi_acd_0 + g_x_0_xzzzzz_1[i] * wa_x[i];

        g_xx_0_yyyyyy_0[i] = g_0_0_yyyyyy_0[i] * fbe_0 - g_0_0_yyyyyy_1[i] * fz_be_0 + g_x_0_yyyyyy_1[i] * wa_x[i];

        g_xx_0_yyyyyz_0[i] = g_0_0_yyyyyz_0[i] * fbe_0 - g_0_0_yyyyyz_1[i] * fz_be_0 + g_x_0_yyyyyz_1[i] * wa_x[i];

        g_xx_0_yyyyzz_0[i] = g_0_0_yyyyzz_0[i] * fbe_0 - g_0_0_yyyyzz_1[i] * fz_be_0 + g_x_0_yyyyzz_1[i] * wa_x[i];

        g_xx_0_yyyzzz_0[i] = g_0_0_yyyzzz_0[i] * fbe_0 - g_0_0_yyyzzz_1[i] * fz_be_0 + g_x_0_yyyzzz_1[i] * wa_x[i];

        g_xx_0_yyzzzz_0[i] = g_0_0_yyzzzz_0[i] * fbe_0 - g_0_0_yyzzzz_1[i] * fz_be_0 + g_x_0_yyzzzz_1[i] * wa_x[i];

        g_xx_0_yzzzzz_0[i] = g_0_0_yzzzzz_0[i] * fbe_0 - g_0_0_yzzzzz_1[i] * fz_be_0 + g_x_0_yzzzzz_1[i] * wa_x[i];

        g_xx_0_zzzzzz_0[i] = g_0_0_zzzzzz_0[i] * fbe_0 - g_0_0_zzzzzz_1[i] * fz_be_0 + g_x_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 28-56 components of targeted buffer : DSI

    auto g_xy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi + 28);

    auto g_xy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 29);

    auto g_xy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 30);

    auto g_xy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 31);

    auto g_xy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 32);

    auto g_xy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 33);

    auto g_xy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 34);

    auto g_xy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 35);

    auto g_xy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 36);

    auto g_xy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 37);

    auto g_xy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 38);

    auto g_xy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 39);

    auto g_xy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 40);

    auto g_xy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 41);

    auto g_xy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 42);

    auto g_xy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 43);

    auto g_xy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 44);

    auto g_xy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 45);

    auto g_xy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 46);

    auto g_xy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 47);

    auto g_xy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 48);

    auto g_xy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 49);

    auto g_xy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 50);

    auto g_xy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 51);

    auto g_xy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 52);

    auto g_xy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 53);

    auto g_xy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 54);

    auto g_xy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 55);

    #pragma omp simd aligned(g_x_0_xxxxxx_1, g_x_0_xxxxxz_1, g_x_0_xxxxzz_1, g_x_0_xxxzzz_1, g_x_0_xxzzzz_1, g_x_0_xzzzzz_1, g_xy_0_xxxxxx_0, g_xy_0_xxxxxy_0, g_xy_0_xxxxxz_0, g_xy_0_xxxxyy_0, g_xy_0_xxxxyz_0, g_xy_0_xxxxzz_0, g_xy_0_xxxyyy_0, g_xy_0_xxxyyz_0, g_xy_0_xxxyzz_0, g_xy_0_xxxzzz_0, g_xy_0_xxyyyy_0, g_xy_0_xxyyyz_0, g_xy_0_xxyyzz_0, g_xy_0_xxyzzz_0, g_xy_0_xxzzzz_0, g_xy_0_xyyyyy_0, g_xy_0_xyyyyz_0, g_xy_0_xyyyzz_0, g_xy_0_xyyzzz_0, g_xy_0_xyzzzz_0, g_xy_0_xzzzzz_0, g_xy_0_yyyyyy_0, g_xy_0_yyyyyz_0, g_xy_0_yyyyzz_0, g_xy_0_yyyzzz_0, g_xy_0_yyzzzz_0, g_xy_0_yzzzzz_0, g_xy_0_zzzzzz_0, g_y_0_xxxxxy_1, g_y_0_xxxxy_1, g_y_0_xxxxyy_1, g_y_0_xxxxyz_1, g_y_0_xxxyy_1, g_y_0_xxxyyy_1, g_y_0_xxxyyz_1, g_y_0_xxxyz_1, g_y_0_xxxyzz_1, g_y_0_xxyyy_1, g_y_0_xxyyyy_1, g_y_0_xxyyyz_1, g_y_0_xxyyz_1, g_y_0_xxyyzz_1, g_y_0_xxyzz_1, g_y_0_xxyzzz_1, g_y_0_xyyyy_1, g_y_0_xyyyyy_1, g_y_0_xyyyyz_1, g_y_0_xyyyz_1, g_y_0_xyyyzz_1, g_y_0_xyyzz_1, g_y_0_xyyzzz_1, g_y_0_xyzzz_1, g_y_0_xyzzzz_1, g_y_0_yyyyy_1, g_y_0_yyyyyy_1, g_y_0_yyyyyz_1, g_y_0_yyyyz_1, g_y_0_yyyyzz_1, g_y_0_yyyzz_1, g_y_0_yyyzzz_1, g_y_0_yyzzz_1, g_y_0_yyzzzz_1, g_y_0_yzzzz_1, g_y_0_yzzzzz_1, g_y_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xxxxxx_0[i] = g_x_0_xxxxxx_1[i] * wa_y[i];

        g_xy_0_xxxxxy_0[i] = 5.0 * g_y_0_xxxxy_1[i] * fi_acd_0 + g_y_0_xxxxxy_1[i] * wa_x[i];

        g_xy_0_xxxxxz_0[i] = g_x_0_xxxxxz_1[i] * wa_y[i];

        g_xy_0_xxxxyy_0[i] = 4.0 * g_y_0_xxxyy_1[i] * fi_acd_0 + g_y_0_xxxxyy_1[i] * wa_x[i];

        g_xy_0_xxxxyz_0[i] = 4.0 * g_y_0_xxxyz_1[i] * fi_acd_0 + g_y_0_xxxxyz_1[i] * wa_x[i];

        g_xy_0_xxxxzz_0[i] = g_x_0_xxxxzz_1[i] * wa_y[i];

        g_xy_0_xxxyyy_0[i] = 3.0 * g_y_0_xxyyy_1[i] * fi_acd_0 + g_y_0_xxxyyy_1[i] * wa_x[i];

        g_xy_0_xxxyyz_0[i] = 3.0 * g_y_0_xxyyz_1[i] * fi_acd_0 + g_y_0_xxxyyz_1[i] * wa_x[i];

        g_xy_0_xxxyzz_0[i] = 3.0 * g_y_0_xxyzz_1[i] * fi_acd_0 + g_y_0_xxxyzz_1[i] * wa_x[i];

        g_xy_0_xxxzzz_0[i] = g_x_0_xxxzzz_1[i] * wa_y[i];

        g_xy_0_xxyyyy_0[i] = 2.0 * g_y_0_xyyyy_1[i] * fi_acd_0 + g_y_0_xxyyyy_1[i] * wa_x[i];

        g_xy_0_xxyyyz_0[i] = 2.0 * g_y_0_xyyyz_1[i] * fi_acd_0 + g_y_0_xxyyyz_1[i] * wa_x[i];

        g_xy_0_xxyyzz_0[i] = 2.0 * g_y_0_xyyzz_1[i] * fi_acd_0 + g_y_0_xxyyzz_1[i] * wa_x[i];

        g_xy_0_xxyzzz_0[i] = 2.0 * g_y_0_xyzzz_1[i] * fi_acd_0 + g_y_0_xxyzzz_1[i] * wa_x[i];

        g_xy_0_xxzzzz_0[i] = g_x_0_xxzzzz_1[i] * wa_y[i];

        g_xy_0_xyyyyy_0[i] = g_y_0_yyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyy_1[i] * wa_x[i];

        g_xy_0_xyyyyz_0[i] = g_y_0_yyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyz_1[i] * wa_x[i];

        g_xy_0_xyyyzz_0[i] = g_y_0_yyyzz_1[i] * fi_acd_0 + g_y_0_xyyyzz_1[i] * wa_x[i];

        g_xy_0_xyyzzz_0[i] = g_y_0_yyzzz_1[i] * fi_acd_0 + g_y_0_xyyzzz_1[i] * wa_x[i];

        g_xy_0_xyzzzz_0[i] = g_y_0_yzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzz_1[i] * wa_x[i];

        g_xy_0_xzzzzz_0[i] = g_x_0_xzzzzz_1[i] * wa_y[i];

        g_xy_0_yyyyyy_0[i] = g_y_0_yyyyyy_1[i] * wa_x[i];

        g_xy_0_yyyyyz_0[i] = g_y_0_yyyyyz_1[i] * wa_x[i];

        g_xy_0_yyyyzz_0[i] = g_y_0_yyyyzz_1[i] * wa_x[i];

        g_xy_0_yyyzzz_0[i] = g_y_0_yyyzzz_1[i] * wa_x[i];

        g_xy_0_yyzzzz_0[i] = g_y_0_yyzzzz_1[i] * wa_x[i];

        g_xy_0_yzzzzz_0[i] = g_y_0_yzzzzz_1[i] * wa_x[i];

        g_xy_0_zzzzzz_0[i] = g_y_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 56-84 components of targeted buffer : DSI

    auto g_xz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi + 56);

    auto g_xz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 57);

    auto g_xz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 58);

    auto g_xz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 59);

    auto g_xz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 60);

    auto g_xz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 61);

    auto g_xz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 62);

    auto g_xz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 63);

    auto g_xz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 64);

    auto g_xz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 65);

    auto g_xz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 66);

    auto g_xz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 67);

    auto g_xz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 68);

    auto g_xz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 69);

    auto g_xz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 70);

    auto g_xz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 71);

    auto g_xz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 72);

    auto g_xz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 73);

    auto g_xz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 74);

    auto g_xz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 75);

    auto g_xz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 76);

    auto g_xz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 77);

    auto g_xz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 78);

    auto g_xz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 79);

    auto g_xz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 80);

    auto g_xz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 81);

    auto g_xz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 82);

    auto g_xz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 83);

    #pragma omp simd aligned(g_x_0_xxxxxx_1, g_x_0_xxxxxy_1, g_x_0_xxxxyy_1, g_x_0_xxxyyy_1, g_x_0_xxyyyy_1, g_x_0_xyyyyy_1, g_xz_0_xxxxxx_0, g_xz_0_xxxxxy_0, g_xz_0_xxxxxz_0, g_xz_0_xxxxyy_0, g_xz_0_xxxxyz_0, g_xz_0_xxxxzz_0, g_xz_0_xxxyyy_0, g_xz_0_xxxyyz_0, g_xz_0_xxxyzz_0, g_xz_0_xxxzzz_0, g_xz_0_xxyyyy_0, g_xz_0_xxyyyz_0, g_xz_0_xxyyzz_0, g_xz_0_xxyzzz_0, g_xz_0_xxzzzz_0, g_xz_0_xyyyyy_0, g_xz_0_xyyyyz_0, g_xz_0_xyyyzz_0, g_xz_0_xyyzzz_0, g_xz_0_xyzzzz_0, g_xz_0_xzzzzz_0, g_xz_0_yyyyyy_0, g_xz_0_yyyyyz_0, g_xz_0_yyyyzz_0, g_xz_0_yyyzzz_0, g_xz_0_yyzzzz_0, g_xz_0_yzzzzz_0, g_xz_0_zzzzzz_0, g_z_0_xxxxxz_1, g_z_0_xxxxyz_1, g_z_0_xxxxz_1, g_z_0_xxxxzz_1, g_z_0_xxxyyz_1, g_z_0_xxxyz_1, g_z_0_xxxyzz_1, g_z_0_xxxzz_1, g_z_0_xxxzzz_1, g_z_0_xxyyyz_1, g_z_0_xxyyz_1, g_z_0_xxyyzz_1, g_z_0_xxyzz_1, g_z_0_xxyzzz_1, g_z_0_xxzzz_1, g_z_0_xxzzzz_1, g_z_0_xyyyyz_1, g_z_0_xyyyz_1, g_z_0_xyyyzz_1, g_z_0_xyyzz_1, g_z_0_xyyzzz_1, g_z_0_xyzzz_1, g_z_0_xyzzzz_1, g_z_0_xzzzz_1, g_z_0_xzzzzz_1, g_z_0_yyyyyy_1, g_z_0_yyyyyz_1, g_z_0_yyyyz_1, g_z_0_yyyyzz_1, g_z_0_yyyzz_1, g_z_0_yyyzzz_1, g_z_0_yyzzz_1, g_z_0_yyzzzz_1, g_z_0_yzzzz_1, g_z_0_yzzzzz_1, g_z_0_zzzzz_1, g_z_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xxxxxx_0[i] = g_x_0_xxxxxx_1[i] * wa_z[i];

        g_xz_0_xxxxxy_0[i] = g_x_0_xxxxxy_1[i] * wa_z[i];

        g_xz_0_xxxxxz_0[i] = 5.0 * g_z_0_xxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxz_1[i] * wa_x[i];

        g_xz_0_xxxxyy_0[i] = g_x_0_xxxxyy_1[i] * wa_z[i];

        g_xz_0_xxxxyz_0[i] = 4.0 * g_z_0_xxxyz_1[i] * fi_acd_0 + g_z_0_xxxxyz_1[i] * wa_x[i];

        g_xz_0_xxxxzz_0[i] = 4.0 * g_z_0_xxxzz_1[i] * fi_acd_0 + g_z_0_xxxxzz_1[i] * wa_x[i];

        g_xz_0_xxxyyy_0[i] = g_x_0_xxxyyy_1[i] * wa_z[i];

        g_xz_0_xxxyyz_0[i] = 3.0 * g_z_0_xxyyz_1[i] * fi_acd_0 + g_z_0_xxxyyz_1[i] * wa_x[i];

        g_xz_0_xxxyzz_0[i] = 3.0 * g_z_0_xxyzz_1[i] * fi_acd_0 + g_z_0_xxxyzz_1[i] * wa_x[i];

        g_xz_0_xxxzzz_0[i] = 3.0 * g_z_0_xxzzz_1[i] * fi_acd_0 + g_z_0_xxxzzz_1[i] * wa_x[i];

        g_xz_0_xxyyyy_0[i] = g_x_0_xxyyyy_1[i] * wa_z[i];

        g_xz_0_xxyyyz_0[i] = 2.0 * g_z_0_xyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyz_1[i] * wa_x[i];

        g_xz_0_xxyyzz_0[i] = 2.0 * g_z_0_xyyzz_1[i] * fi_acd_0 + g_z_0_xxyyzz_1[i] * wa_x[i];

        g_xz_0_xxyzzz_0[i] = 2.0 * g_z_0_xyzzz_1[i] * fi_acd_0 + g_z_0_xxyzzz_1[i] * wa_x[i];

        g_xz_0_xxzzzz_0[i] = 2.0 * g_z_0_xzzzz_1[i] * fi_acd_0 + g_z_0_xxzzzz_1[i] * wa_x[i];

        g_xz_0_xyyyyy_0[i] = g_x_0_xyyyyy_1[i] * wa_z[i];

        g_xz_0_xyyyyz_0[i] = g_z_0_yyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyz_1[i] * wa_x[i];

        g_xz_0_xyyyzz_0[i] = g_z_0_yyyzz_1[i] * fi_acd_0 + g_z_0_xyyyzz_1[i] * wa_x[i];

        g_xz_0_xyyzzz_0[i] = g_z_0_yyzzz_1[i] * fi_acd_0 + g_z_0_xyyzzz_1[i] * wa_x[i];

        g_xz_0_xyzzzz_0[i] = g_z_0_yzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzz_1[i] * wa_x[i];

        g_xz_0_xzzzzz_0[i] = g_z_0_zzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzz_1[i] * wa_x[i];

        g_xz_0_yyyyyy_0[i] = g_z_0_yyyyyy_1[i] * wa_x[i];

        g_xz_0_yyyyyz_0[i] = g_z_0_yyyyyz_1[i] * wa_x[i];

        g_xz_0_yyyyzz_0[i] = g_z_0_yyyyzz_1[i] * wa_x[i];

        g_xz_0_yyyzzz_0[i] = g_z_0_yyyzzz_1[i] * wa_x[i];

        g_xz_0_yyzzzz_0[i] = g_z_0_yyzzzz_1[i] * wa_x[i];

        g_xz_0_yzzzzz_0[i] = g_z_0_yzzzzz_1[i] * wa_x[i];

        g_xz_0_zzzzzz_0[i] = g_z_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 84-112 components of targeted buffer : DSI

    auto g_yy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi + 84);

    auto g_yy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 85);

    auto g_yy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 86);

    auto g_yy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 87);

    auto g_yy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 88);

    auto g_yy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 89);

    auto g_yy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 90);

    auto g_yy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 91);

    auto g_yy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 92);

    auto g_yy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 93);

    auto g_yy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 94);

    auto g_yy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 95);

    auto g_yy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 96);

    auto g_yy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 97);

    auto g_yy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 98);

    auto g_yy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 99);

    auto g_yy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 100);

    auto g_yy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 101);

    auto g_yy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 102);

    auto g_yy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 103);

    auto g_yy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 104);

    auto g_yy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 105);

    auto g_yy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 106);

    auto g_yy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 107);

    auto g_yy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 108);

    auto g_yy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 109);

    auto g_yy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 110);

    auto g_yy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 111);

    #pragma omp simd aligned(g_0_0_xxxxxx_0, g_0_0_xxxxxx_1, g_0_0_xxxxxy_0, g_0_0_xxxxxy_1, g_0_0_xxxxxz_0, g_0_0_xxxxxz_1, g_0_0_xxxxyy_0, g_0_0_xxxxyy_1, g_0_0_xxxxyz_0, g_0_0_xxxxyz_1, g_0_0_xxxxzz_0, g_0_0_xxxxzz_1, g_0_0_xxxyyy_0, g_0_0_xxxyyy_1, g_0_0_xxxyyz_0, g_0_0_xxxyyz_1, g_0_0_xxxyzz_0, g_0_0_xxxyzz_1, g_0_0_xxxzzz_0, g_0_0_xxxzzz_1, g_0_0_xxyyyy_0, g_0_0_xxyyyy_1, g_0_0_xxyyyz_0, g_0_0_xxyyyz_1, g_0_0_xxyyzz_0, g_0_0_xxyyzz_1, g_0_0_xxyzzz_0, g_0_0_xxyzzz_1, g_0_0_xxzzzz_0, g_0_0_xxzzzz_1, g_0_0_xyyyyy_0, g_0_0_xyyyyy_1, g_0_0_xyyyyz_0, g_0_0_xyyyyz_1, g_0_0_xyyyzz_0, g_0_0_xyyyzz_1, g_0_0_xyyzzz_0, g_0_0_xyyzzz_1, g_0_0_xyzzzz_0, g_0_0_xyzzzz_1, g_0_0_xzzzzz_0, g_0_0_xzzzzz_1, g_0_0_yyyyyy_0, g_0_0_yyyyyy_1, g_0_0_yyyyyz_0, g_0_0_yyyyyz_1, g_0_0_yyyyzz_0, g_0_0_yyyyzz_1, g_0_0_yyyzzz_0, g_0_0_yyyzzz_1, g_0_0_yyzzzz_0, g_0_0_yyzzzz_1, g_0_0_yzzzzz_0, g_0_0_yzzzzz_1, g_0_0_zzzzzz_0, g_0_0_zzzzzz_1, g_y_0_xxxxx_1, g_y_0_xxxxxx_1, g_y_0_xxxxxy_1, g_y_0_xxxxxz_1, g_y_0_xxxxy_1, g_y_0_xxxxyy_1, g_y_0_xxxxyz_1, g_y_0_xxxxz_1, g_y_0_xxxxzz_1, g_y_0_xxxyy_1, g_y_0_xxxyyy_1, g_y_0_xxxyyz_1, g_y_0_xxxyz_1, g_y_0_xxxyzz_1, g_y_0_xxxzz_1, g_y_0_xxxzzz_1, g_y_0_xxyyy_1, g_y_0_xxyyyy_1, g_y_0_xxyyyz_1, g_y_0_xxyyz_1, g_y_0_xxyyzz_1, g_y_0_xxyzz_1, g_y_0_xxyzzz_1, g_y_0_xxzzz_1, g_y_0_xxzzzz_1, g_y_0_xyyyy_1, g_y_0_xyyyyy_1, g_y_0_xyyyyz_1, g_y_0_xyyyz_1, g_y_0_xyyyzz_1, g_y_0_xyyzz_1, g_y_0_xyyzzz_1, g_y_0_xyzzz_1, g_y_0_xyzzzz_1, g_y_0_xzzzz_1, g_y_0_xzzzzz_1, g_y_0_yyyyy_1, g_y_0_yyyyyy_1, g_y_0_yyyyyz_1, g_y_0_yyyyz_1, g_y_0_yyyyzz_1, g_y_0_yyyzz_1, g_y_0_yyyzzz_1, g_y_0_yyzzz_1, g_y_0_yyzzzz_1, g_y_0_yzzzz_1, g_y_0_yzzzzz_1, g_y_0_zzzzz_1, g_y_0_zzzzzz_1, g_yy_0_xxxxxx_0, g_yy_0_xxxxxy_0, g_yy_0_xxxxxz_0, g_yy_0_xxxxyy_0, g_yy_0_xxxxyz_0, g_yy_0_xxxxzz_0, g_yy_0_xxxyyy_0, g_yy_0_xxxyyz_0, g_yy_0_xxxyzz_0, g_yy_0_xxxzzz_0, g_yy_0_xxyyyy_0, g_yy_0_xxyyyz_0, g_yy_0_xxyyzz_0, g_yy_0_xxyzzz_0, g_yy_0_xxzzzz_0, g_yy_0_xyyyyy_0, g_yy_0_xyyyyz_0, g_yy_0_xyyyzz_0, g_yy_0_xyyzzz_0, g_yy_0_xyzzzz_0, g_yy_0_xzzzzz_0, g_yy_0_yyyyyy_0, g_yy_0_yyyyyz_0, g_yy_0_yyyyzz_0, g_yy_0_yyyzzz_0, g_yy_0_yyzzzz_0, g_yy_0_yzzzzz_0, g_yy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xxxxxx_0[i] = g_0_0_xxxxxx_0[i] * fbe_0 - g_0_0_xxxxxx_1[i] * fz_be_0 + g_y_0_xxxxxx_1[i] * wa_y[i];

        g_yy_0_xxxxxy_0[i] = g_0_0_xxxxxy_0[i] * fbe_0 - g_0_0_xxxxxy_1[i] * fz_be_0 + g_y_0_xxxxx_1[i] * fi_acd_0 + g_y_0_xxxxxy_1[i] * wa_y[i];

        g_yy_0_xxxxxz_0[i] = g_0_0_xxxxxz_0[i] * fbe_0 - g_0_0_xxxxxz_1[i] * fz_be_0 + g_y_0_xxxxxz_1[i] * wa_y[i];

        g_yy_0_xxxxyy_0[i] = g_0_0_xxxxyy_0[i] * fbe_0 - g_0_0_xxxxyy_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxy_1[i] * fi_acd_0 + g_y_0_xxxxyy_1[i] * wa_y[i];

        g_yy_0_xxxxyz_0[i] = g_0_0_xxxxyz_0[i] * fbe_0 - g_0_0_xxxxyz_1[i] * fz_be_0 + g_y_0_xxxxz_1[i] * fi_acd_0 + g_y_0_xxxxyz_1[i] * wa_y[i];

        g_yy_0_xxxxzz_0[i] = g_0_0_xxxxzz_0[i] * fbe_0 - g_0_0_xxxxzz_1[i] * fz_be_0 + g_y_0_xxxxzz_1[i] * wa_y[i];

        g_yy_0_xxxyyy_0[i] = g_0_0_xxxyyy_0[i] * fbe_0 - g_0_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_y_0_xxxyy_1[i] * fi_acd_0 + g_y_0_xxxyyy_1[i] * wa_y[i];

        g_yy_0_xxxyyz_0[i] = g_0_0_xxxyyz_0[i] * fbe_0 - g_0_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxyz_1[i] * fi_acd_0 + g_y_0_xxxyyz_1[i] * wa_y[i];

        g_yy_0_xxxyzz_0[i] = g_0_0_xxxyzz_0[i] * fbe_0 - g_0_0_xxxyzz_1[i] * fz_be_0 + g_y_0_xxxzz_1[i] * fi_acd_0 + g_y_0_xxxyzz_1[i] * wa_y[i];

        g_yy_0_xxxzzz_0[i] = g_0_0_xxxzzz_0[i] * fbe_0 - g_0_0_xxxzzz_1[i] * fz_be_0 + g_y_0_xxxzzz_1[i] * wa_y[i];

        g_yy_0_xxyyyy_0[i] = g_0_0_xxyyyy_0[i] * fbe_0 - g_0_0_xxyyyy_1[i] * fz_be_0 + 4.0 * g_y_0_xxyyy_1[i] * fi_acd_0 + g_y_0_xxyyyy_1[i] * wa_y[i];

        g_yy_0_xxyyyz_0[i] = g_0_0_xxyyyz_0[i] * fbe_0 - g_0_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_y_0_xxyyz_1[i] * fi_acd_0 + g_y_0_xxyyyz_1[i] * wa_y[i];

        g_yy_0_xxyyzz_0[i] = g_0_0_xxyyzz_0[i] * fbe_0 - g_0_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxyzz_1[i] * fi_acd_0 + g_y_0_xxyyzz_1[i] * wa_y[i];

        g_yy_0_xxyzzz_0[i] = g_0_0_xxyzzz_0[i] * fbe_0 - g_0_0_xxyzzz_1[i] * fz_be_0 + g_y_0_xxzzz_1[i] * fi_acd_0 + g_y_0_xxyzzz_1[i] * wa_y[i];

        g_yy_0_xxzzzz_0[i] = g_0_0_xxzzzz_0[i] * fbe_0 - g_0_0_xxzzzz_1[i] * fz_be_0 + g_y_0_xxzzzz_1[i] * wa_y[i];

        g_yy_0_xyyyyy_0[i] = g_0_0_xyyyyy_0[i] * fbe_0 - g_0_0_xyyyyy_1[i] * fz_be_0 + 5.0 * g_y_0_xyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyy_1[i] * wa_y[i];

        g_yy_0_xyyyyz_0[i] = g_0_0_xyyyyz_0[i] * fbe_0 - g_0_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_y_0_xyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyz_1[i] * wa_y[i];

        g_yy_0_xyyyzz_0[i] = g_0_0_xyyyzz_0[i] * fbe_0 - g_0_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_y_0_xyyzz_1[i] * fi_acd_0 + g_y_0_xyyyzz_1[i] * wa_y[i];

        g_yy_0_xyyzzz_0[i] = g_0_0_xyyzzz_0[i] * fbe_0 - g_0_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xyzzz_1[i] * fi_acd_0 + g_y_0_xyyzzz_1[i] * wa_y[i];

        g_yy_0_xyzzzz_0[i] = g_0_0_xyzzzz_0[i] * fbe_0 - g_0_0_xyzzzz_1[i] * fz_be_0 + g_y_0_xzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzz_1[i] * wa_y[i];

        g_yy_0_xzzzzz_0[i] = g_0_0_xzzzzz_0[i] * fbe_0 - g_0_0_xzzzzz_1[i] * fz_be_0 + g_y_0_xzzzzz_1[i] * wa_y[i];

        g_yy_0_yyyyyy_0[i] = g_0_0_yyyyyy_0[i] * fbe_0 - g_0_0_yyyyyy_1[i] * fz_be_0 + 6.0 * g_y_0_yyyyy_1[i] * fi_acd_0 + g_y_0_yyyyyy_1[i] * wa_y[i];

        g_yy_0_yyyyyz_0[i] = g_0_0_yyyyyz_0[i] * fbe_0 - g_0_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_y_0_yyyyz_1[i] * fi_acd_0 + g_y_0_yyyyyz_1[i] * wa_y[i];

        g_yy_0_yyyyzz_0[i] = g_0_0_yyyyzz_0[i] * fbe_0 - g_0_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_y_0_yyyzz_1[i] * fi_acd_0 + g_y_0_yyyyzz_1[i] * wa_y[i];

        g_yy_0_yyyzzz_0[i] = g_0_0_yyyzzz_0[i] * fbe_0 - g_0_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_y_0_yyzzz_1[i] * fi_acd_0 + g_y_0_yyyzzz_1[i] * wa_y[i];

        g_yy_0_yyzzzz_0[i] = g_0_0_yyzzzz_0[i] * fbe_0 - g_0_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_yzzzz_1[i] * fi_acd_0 + g_y_0_yyzzzz_1[i] * wa_y[i];

        g_yy_0_yzzzzz_0[i] = g_0_0_yzzzzz_0[i] * fbe_0 - g_0_0_yzzzzz_1[i] * fz_be_0 + g_y_0_zzzzz_1[i] * fi_acd_0 + g_y_0_yzzzzz_1[i] * wa_y[i];

        g_yy_0_zzzzzz_0[i] = g_0_0_zzzzzz_0[i] * fbe_0 - g_0_0_zzzzzz_1[i] * fz_be_0 + g_y_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 112-140 components of targeted buffer : DSI

    auto g_yz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi + 112);

    auto g_yz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 113);

    auto g_yz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 114);

    auto g_yz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 115);

    auto g_yz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 116);

    auto g_yz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 117);

    auto g_yz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 118);

    auto g_yz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 119);

    auto g_yz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 120);

    auto g_yz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 121);

    auto g_yz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 122);

    auto g_yz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 123);

    auto g_yz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 124);

    auto g_yz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 125);

    auto g_yz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 126);

    auto g_yz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 127);

    auto g_yz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 128);

    auto g_yz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 129);

    auto g_yz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 130);

    auto g_yz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 131);

    auto g_yz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 132);

    auto g_yz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 133);

    auto g_yz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 134);

    auto g_yz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 135);

    auto g_yz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 136);

    auto g_yz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 137);

    auto g_yz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 138);

    auto g_yz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 139);

    #pragma omp simd aligned(g_y_0_xxxxxy_1, g_y_0_xxxxyy_1, g_y_0_xxxyyy_1, g_y_0_xxyyyy_1, g_y_0_xyyyyy_1, g_y_0_yyyyyy_1, g_yz_0_xxxxxx_0, g_yz_0_xxxxxy_0, g_yz_0_xxxxxz_0, g_yz_0_xxxxyy_0, g_yz_0_xxxxyz_0, g_yz_0_xxxxzz_0, g_yz_0_xxxyyy_0, g_yz_0_xxxyyz_0, g_yz_0_xxxyzz_0, g_yz_0_xxxzzz_0, g_yz_0_xxyyyy_0, g_yz_0_xxyyyz_0, g_yz_0_xxyyzz_0, g_yz_0_xxyzzz_0, g_yz_0_xxzzzz_0, g_yz_0_xyyyyy_0, g_yz_0_xyyyyz_0, g_yz_0_xyyyzz_0, g_yz_0_xyyzzz_0, g_yz_0_xyzzzz_0, g_yz_0_xzzzzz_0, g_yz_0_yyyyyy_0, g_yz_0_yyyyyz_0, g_yz_0_yyyyzz_0, g_yz_0_yyyzzz_0, g_yz_0_yyzzzz_0, g_yz_0_yzzzzz_0, g_yz_0_zzzzzz_0, g_z_0_xxxxxx_1, g_z_0_xxxxxz_1, g_z_0_xxxxyz_1, g_z_0_xxxxz_1, g_z_0_xxxxzz_1, g_z_0_xxxyyz_1, g_z_0_xxxyz_1, g_z_0_xxxyzz_1, g_z_0_xxxzz_1, g_z_0_xxxzzz_1, g_z_0_xxyyyz_1, g_z_0_xxyyz_1, g_z_0_xxyyzz_1, g_z_0_xxyzz_1, g_z_0_xxyzzz_1, g_z_0_xxzzz_1, g_z_0_xxzzzz_1, g_z_0_xyyyyz_1, g_z_0_xyyyz_1, g_z_0_xyyyzz_1, g_z_0_xyyzz_1, g_z_0_xyyzzz_1, g_z_0_xyzzz_1, g_z_0_xyzzzz_1, g_z_0_xzzzz_1, g_z_0_xzzzzz_1, g_z_0_yyyyyz_1, g_z_0_yyyyz_1, g_z_0_yyyyzz_1, g_z_0_yyyzz_1, g_z_0_yyyzzz_1, g_z_0_yyzzz_1, g_z_0_yyzzzz_1, g_z_0_yzzzz_1, g_z_0_yzzzzz_1, g_z_0_zzzzz_1, g_z_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xxxxxx_0[i] = g_z_0_xxxxxx_1[i] * wa_y[i];

        g_yz_0_xxxxxy_0[i] = g_y_0_xxxxxy_1[i] * wa_z[i];

        g_yz_0_xxxxxz_0[i] = g_z_0_xxxxxz_1[i] * wa_y[i];

        g_yz_0_xxxxyy_0[i] = g_y_0_xxxxyy_1[i] * wa_z[i];

        g_yz_0_xxxxyz_0[i] = g_z_0_xxxxz_1[i] * fi_acd_0 + g_z_0_xxxxyz_1[i] * wa_y[i];

        g_yz_0_xxxxzz_0[i] = g_z_0_xxxxzz_1[i] * wa_y[i];

        g_yz_0_xxxyyy_0[i] = g_y_0_xxxyyy_1[i] * wa_z[i];

        g_yz_0_xxxyyz_0[i] = 2.0 * g_z_0_xxxyz_1[i] * fi_acd_0 + g_z_0_xxxyyz_1[i] * wa_y[i];

        g_yz_0_xxxyzz_0[i] = g_z_0_xxxzz_1[i] * fi_acd_0 + g_z_0_xxxyzz_1[i] * wa_y[i];

        g_yz_0_xxxzzz_0[i] = g_z_0_xxxzzz_1[i] * wa_y[i];

        g_yz_0_xxyyyy_0[i] = g_y_0_xxyyyy_1[i] * wa_z[i];

        g_yz_0_xxyyyz_0[i] = 3.0 * g_z_0_xxyyz_1[i] * fi_acd_0 + g_z_0_xxyyyz_1[i] * wa_y[i];

        g_yz_0_xxyyzz_0[i] = 2.0 * g_z_0_xxyzz_1[i] * fi_acd_0 + g_z_0_xxyyzz_1[i] * wa_y[i];

        g_yz_0_xxyzzz_0[i] = g_z_0_xxzzz_1[i] * fi_acd_0 + g_z_0_xxyzzz_1[i] * wa_y[i];

        g_yz_0_xxzzzz_0[i] = g_z_0_xxzzzz_1[i] * wa_y[i];

        g_yz_0_xyyyyy_0[i] = g_y_0_xyyyyy_1[i] * wa_z[i];

        g_yz_0_xyyyyz_0[i] = 4.0 * g_z_0_xyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyz_1[i] * wa_y[i];

        g_yz_0_xyyyzz_0[i] = 3.0 * g_z_0_xyyzz_1[i] * fi_acd_0 + g_z_0_xyyyzz_1[i] * wa_y[i];

        g_yz_0_xyyzzz_0[i] = 2.0 * g_z_0_xyzzz_1[i] * fi_acd_0 + g_z_0_xyyzzz_1[i] * wa_y[i];

        g_yz_0_xyzzzz_0[i] = g_z_0_xzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzz_1[i] * wa_y[i];

        g_yz_0_xzzzzz_0[i] = g_z_0_xzzzzz_1[i] * wa_y[i];

        g_yz_0_yyyyyy_0[i] = g_y_0_yyyyyy_1[i] * wa_z[i];

        g_yz_0_yyyyyz_0[i] = 5.0 * g_z_0_yyyyz_1[i] * fi_acd_0 + g_z_0_yyyyyz_1[i] * wa_y[i];

        g_yz_0_yyyyzz_0[i] = 4.0 * g_z_0_yyyzz_1[i] * fi_acd_0 + g_z_0_yyyyzz_1[i] * wa_y[i];

        g_yz_0_yyyzzz_0[i] = 3.0 * g_z_0_yyzzz_1[i] * fi_acd_0 + g_z_0_yyyzzz_1[i] * wa_y[i];

        g_yz_0_yyzzzz_0[i] = 2.0 * g_z_0_yzzzz_1[i] * fi_acd_0 + g_z_0_yyzzzz_1[i] * wa_y[i];

        g_yz_0_yzzzzz_0[i] = g_z_0_zzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzz_1[i] * wa_y[i];

        g_yz_0_zzzzzz_0[i] = g_z_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 140-168 components of targeted buffer : DSI

    auto g_zz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi + 140);

    auto g_zz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 141);

    auto g_zz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 142);

    auto g_zz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 143);

    auto g_zz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 144);

    auto g_zz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 145);

    auto g_zz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 146);

    auto g_zz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 147);

    auto g_zz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 148);

    auto g_zz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 149);

    auto g_zz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 150);

    auto g_zz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 151);

    auto g_zz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 152);

    auto g_zz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 153);

    auto g_zz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 154);

    auto g_zz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 155);

    auto g_zz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 156);

    auto g_zz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 157);

    auto g_zz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 158);

    auto g_zz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 159);

    auto g_zz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 160);

    auto g_zz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 161);

    auto g_zz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 162);

    auto g_zz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 163);

    auto g_zz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 164);

    auto g_zz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 165);

    auto g_zz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 166);

    auto g_zz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 167);

    #pragma omp simd aligned(g_0_0_xxxxxx_0, g_0_0_xxxxxx_1, g_0_0_xxxxxy_0, g_0_0_xxxxxy_1, g_0_0_xxxxxz_0, g_0_0_xxxxxz_1, g_0_0_xxxxyy_0, g_0_0_xxxxyy_1, g_0_0_xxxxyz_0, g_0_0_xxxxyz_1, g_0_0_xxxxzz_0, g_0_0_xxxxzz_1, g_0_0_xxxyyy_0, g_0_0_xxxyyy_1, g_0_0_xxxyyz_0, g_0_0_xxxyyz_1, g_0_0_xxxyzz_0, g_0_0_xxxyzz_1, g_0_0_xxxzzz_0, g_0_0_xxxzzz_1, g_0_0_xxyyyy_0, g_0_0_xxyyyy_1, g_0_0_xxyyyz_0, g_0_0_xxyyyz_1, g_0_0_xxyyzz_0, g_0_0_xxyyzz_1, g_0_0_xxyzzz_0, g_0_0_xxyzzz_1, g_0_0_xxzzzz_0, g_0_0_xxzzzz_1, g_0_0_xyyyyy_0, g_0_0_xyyyyy_1, g_0_0_xyyyyz_0, g_0_0_xyyyyz_1, g_0_0_xyyyzz_0, g_0_0_xyyyzz_1, g_0_0_xyyzzz_0, g_0_0_xyyzzz_1, g_0_0_xyzzzz_0, g_0_0_xyzzzz_1, g_0_0_xzzzzz_0, g_0_0_xzzzzz_1, g_0_0_yyyyyy_0, g_0_0_yyyyyy_1, g_0_0_yyyyyz_0, g_0_0_yyyyyz_1, g_0_0_yyyyzz_0, g_0_0_yyyyzz_1, g_0_0_yyyzzz_0, g_0_0_yyyzzz_1, g_0_0_yyzzzz_0, g_0_0_yyzzzz_1, g_0_0_yzzzzz_0, g_0_0_yzzzzz_1, g_0_0_zzzzzz_0, g_0_0_zzzzzz_1, g_z_0_xxxxx_1, g_z_0_xxxxxx_1, g_z_0_xxxxxy_1, g_z_0_xxxxxz_1, g_z_0_xxxxy_1, g_z_0_xxxxyy_1, g_z_0_xxxxyz_1, g_z_0_xxxxz_1, g_z_0_xxxxzz_1, g_z_0_xxxyy_1, g_z_0_xxxyyy_1, g_z_0_xxxyyz_1, g_z_0_xxxyz_1, g_z_0_xxxyzz_1, g_z_0_xxxzz_1, g_z_0_xxxzzz_1, g_z_0_xxyyy_1, g_z_0_xxyyyy_1, g_z_0_xxyyyz_1, g_z_0_xxyyz_1, g_z_0_xxyyzz_1, g_z_0_xxyzz_1, g_z_0_xxyzzz_1, g_z_0_xxzzz_1, g_z_0_xxzzzz_1, g_z_0_xyyyy_1, g_z_0_xyyyyy_1, g_z_0_xyyyyz_1, g_z_0_xyyyz_1, g_z_0_xyyyzz_1, g_z_0_xyyzz_1, g_z_0_xyyzzz_1, g_z_0_xyzzz_1, g_z_0_xyzzzz_1, g_z_0_xzzzz_1, g_z_0_xzzzzz_1, g_z_0_yyyyy_1, g_z_0_yyyyyy_1, g_z_0_yyyyyz_1, g_z_0_yyyyz_1, g_z_0_yyyyzz_1, g_z_0_yyyzz_1, g_z_0_yyyzzz_1, g_z_0_yyzzz_1, g_z_0_yyzzzz_1, g_z_0_yzzzz_1, g_z_0_yzzzzz_1, g_z_0_zzzzz_1, g_z_0_zzzzzz_1, g_zz_0_xxxxxx_0, g_zz_0_xxxxxy_0, g_zz_0_xxxxxz_0, g_zz_0_xxxxyy_0, g_zz_0_xxxxyz_0, g_zz_0_xxxxzz_0, g_zz_0_xxxyyy_0, g_zz_0_xxxyyz_0, g_zz_0_xxxyzz_0, g_zz_0_xxxzzz_0, g_zz_0_xxyyyy_0, g_zz_0_xxyyyz_0, g_zz_0_xxyyzz_0, g_zz_0_xxyzzz_0, g_zz_0_xxzzzz_0, g_zz_0_xyyyyy_0, g_zz_0_xyyyyz_0, g_zz_0_xyyyzz_0, g_zz_0_xyyzzz_0, g_zz_0_xyzzzz_0, g_zz_0_xzzzzz_0, g_zz_0_yyyyyy_0, g_zz_0_yyyyyz_0, g_zz_0_yyyyzz_0, g_zz_0_yyyzzz_0, g_zz_0_yyzzzz_0, g_zz_0_yzzzzz_0, g_zz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xxxxxx_0[i] = g_0_0_xxxxxx_0[i] * fbe_0 - g_0_0_xxxxxx_1[i] * fz_be_0 + g_z_0_xxxxxx_1[i] * wa_z[i];

        g_zz_0_xxxxxy_0[i] = g_0_0_xxxxxy_0[i] * fbe_0 - g_0_0_xxxxxy_1[i] * fz_be_0 + g_z_0_xxxxxy_1[i] * wa_z[i];

        g_zz_0_xxxxxz_0[i] = g_0_0_xxxxxz_0[i] * fbe_0 - g_0_0_xxxxxz_1[i] * fz_be_0 + g_z_0_xxxxx_1[i] * fi_acd_0 + g_z_0_xxxxxz_1[i] * wa_z[i];

        g_zz_0_xxxxyy_0[i] = g_0_0_xxxxyy_0[i] * fbe_0 - g_0_0_xxxxyy_1[i] * fz_be_0 + g_z_0_xxxxyy_1[i] * wa_z[i];

        g_zz_0_xxxxyz_0[i] = g_0_0_xxxxyz_0[i] * fbe_0 - g_0_0_xxxxyz_1[i] * fz_be_0 + g_z_0_xxxxy_1[i] * fi_acd_0 + g_z_0_xxxxyz_1[i] * wa_z[i];

        g_zz_0_xxxxzz_0[i] = g_0_0_xxxxzz_0[i] * fbe_0 - g_0_0_xxxxzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxz_1[i] * fi_acd_0 + g_z_0_xxxxzz_1[i] * wa_z[i];

        g_zz_0_xxxyyy_0[i] = g_0_0_xxxyyy_0[i] * fbe_0 - g_0_0_xxxyyy_1[i] * fz_be_0 + g_z_0_xxxyyy_1[i] * wa_z[i];

        g_zz_0_xxxyyz_0[i] = g_0_0_xxxyyz_0[i] * fbe_0 - g_0_0_xxxyyz_1[i] * fz_be_0 + g_z_0_xxxyy_1[i] * fi_acd_0 + g_z_0_xxxyyz_1[i] * wa_z[i];

        g_zz_0_xxxyzz_0[i] = g_0_0_xxxyzz_0[i] * fbe_0 - g_0_0_xxxyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxyz_1[i] * fi_acd_0 + g_z_0_xxxyzz_1[i] * wa_z[i];

        g_zz_0_xxxzzz_0[i] = g_0_0_xxxzzz_0[i] * fbe_0 - g_0_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxzz_1[i] * fi_acd_0 + g_z_0_xxxzzz_1[i] * wa_z[i];

        g_zz_0_xxyyyy_0[i] = g_0_0_xxyyyy_0[i] * fbe_0 - g_0_0_xxyyyy_1[i] * fz_be_0 + g_z_0_xxyyyy_1[i] * wa_z[i];

        g_zz_0_xxyyyz_0[i] = g_0_0_xxyyyz_0[i] * fbe_0 - g_0_0_xxyyyz_1[i] * fz_be_0 + g_z_0_xxyyy_1[i] * fi_acd_0 + g_z_0_xxyyyz_1[i] * wa_z[i];

        g_zz_0_xxyyzz_0[i] = g_0_0_xxyyzz_0[i] * fbe_0 - g_0_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxyyz_1[i] * fi_acd_0 + g_z_0_xxyyzz_1[i] * wa_z[i];

        g_zz_0_xxyzzz_0[i] = g_0_0_xxyzzz_0[i] * fbe_0 - g_0_0_xxyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxyzz_1[i] * fi_acd_0 + g_z_0_xxyzzz_1[i] * wa_z[i];

        g_zz_0_xxzzzz_0[i] = g_0_0_xxzzzz_0[i] * fbe_0 - g_0_0_xxzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxzzz_1[i] * fi_acd_0 + g_z_0_xxzzzz_1[i] * wa_z[i];

        g_zz_0_xyyyyy_0[i] = g_0_0_xyyyyy_0[i] * fbe_0 - g_0_0_xyyyyy_1[i] * fz_be_0 + g_z_0_xyyyyy_1[i] * wa_z[i];

        g_zz_0_xyyyyz_0[i] = g_0_0_xyyyyz_0[i] * fbe_0 - g_0_0_xyyyyz_1[i] * fz_be_0 + g_z_0_xyyyy_1[i] * fi_acd_0 + g_z_0_xyyyyz_1[i] * wa_z[i];

        g_zz_0_xyyyzz_0[i] = g_0_0_xyyyzz_0[i] * fbe_0 - g_0_0_xyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xyyyz_1[i] * fi_acd_0 + g_z_0_xyyyzz_1[i] * wa_z[i];

        g_zz_0_xyyzzz_0[i] = g_0_0_xyyzzz_0[i] * fbe_0 - g_0_0_xyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xyyzz_1[i] * fi_acd_0 + g_z_0_xyyzzz_1[i] * wa_z[i];

        g_zz_0_xyzzzz_0[i] = g_0_0_xyzzzz_0[i] * fbe_0 - g_0_0_xyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xyzzz_1[i] * fi_acd_0 + g_z_0_xyzzzz_1[i] * wa_z[i];

        g_zz_0_xzzzzz_0[i] = g_0_0_xzzzzz_0[i] * fbe_0 - g_0_0_xzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzz_1[i] * wa_z[i];

        g_zz_0_yyyyyy_0[i] = g_0_0_yyyyyy_0[i] * fbe_0 - g_0_0_yyyyyy_1[i] * fz_be_0 + g_z_0_yyyyyy_1[i] * wa_z[i];

        g_zz_0_yyyyyz_0[i] = g_0_0_yyyyyz_0[i] * fbe_0 - g_0_0_yyyyyz_1[i] * fz_be_0 + g_z_0_yyyyy_1[i] * fi_acd_0 + g_z_0_yyyyyz_1[i] * wa_z[i];

        g_zz_0_yyyyzz_0[i] = g_0_0_yyyyzz_0[i] * fbe_0 - g_0_0_yyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_yyyyz_1[i] * fi_acd_0 + g_z_0_yyyyzz_1[i] * wa_z[i];

        g_zz_0_yyyzzz_0[i] = g_0_0_yyyzzz_0[i] * fbe_0 - g_0_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_yyyzz_1[i] * fi_acd_0 + g_z_0_yyyzzz_1[i] * wa_z[i];

        g_zz_0_yyzzzz_0[i] = g_0_0_yyzzzz_0[i] * fbe_0 - g_0_0_yyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_yyzzz_1[i] * fi_acd_0 + g_z_0_yyzzzz_1[i] * wa_z[i];

        g_zz_0_yzzzzz_0[i] = g_0_0_yzzzzz_0[i] * fbe_0 - g_0_0_yzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_yzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzz_1[i] * wa_z[i];

        g_zz_0_zzzzzz_0[i] = g_0_0_zzzzzz_0[i] * fbe_0 - g_0_0_zzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_zzzzz_1[i] * fi_acd_0 + g_z_0_zzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

