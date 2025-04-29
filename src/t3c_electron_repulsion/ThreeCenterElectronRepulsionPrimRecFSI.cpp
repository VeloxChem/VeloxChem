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

#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsi,
                                 size_t idx_eri_0_psi,
                                 size_t idx_eri_1_psi,
                                 size_t idx_eri_1_dsh,
                                 size_t idx_eri_1_dsi,
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

    /// Set up components of auxilary buffer : PSI

    auto g_x_0_xxxxxx_0 = pbuffer.data(idx_eri_0_psi);

    auto g_x_0_xxxxxy_0 = pbuffer.data(idx_eri_0_psi + 1);

    auto g_x_0_xxxxxz_0 = pbuffer.data(idx_eri_0_psi + 2);

    auto g_x_0_xxxxyy_0 = pbuffer.data(idx_eri_0_psi + 3);

    auto g_x_0_xxxxyz_0 = pbuffer.data(idx_eri_0_psi + 4);

    auto g_x_0_xxxxzz_0 = pbuffer.data(idx_eri_0_psi + 5);

    auto g_x_0_xxxyyy_0 = pbuffer.data(idx_eri_0_psi + 6);

    auto g_x_0_xxxyyz_0 = pbuffer.data(idx_eri_0_psi + 7);

    auto g_x_0_xxxyzz_0 = pbuffer.data(idx_eri_0_psi + 8);

    auto g_x_0_xxxzzz_0 = pbuffer.data(idx_eri_0_psi + 9);

    auto g_x_0_xxyyyy_0 = pbuffer.data(idx_eri_0_psi + 10);

    auto g_x_0_xxyyyz_0 = pbuffer.data(idx_eri_0_psi + 11);

    auto g_x_0_xxyyzz_0 = pbuffer.data(idx_eri_0_psi + 12);

    auto g_x_0_xxyzzz_0 = pbuffer.data(idx_eri_0_psi + 13);

    auto g_x_0_xxzzzz_0 = pbuffer.data(idx_eri_0_psi + 14);

    auto g_x_0_xyyyyy_0 = pbuffer.data(idx_eri_0_psi + 15);

    auto g_x_0_xyyyyz_0 = pbuffer.data(idx_eri_0_psi + 16);

    auto g_x_0_xyyyzz_0 = pbuffer.data(idx_eri_0_psi + 17);

    auto g_x_0_xyyzzz_0 = pbuffer.data(idx_eri_0_psi + 18);

    auto g_x_0_xyzzzz_0 = pbuffer.data(idx_eri_0_psi + 19);

    auto g_x_0_xzzzzz_0 = pbuffer.data(idx_eri_0_psi + 20);

    auto g_x_0_yyyyyy_0 = pbuffer.data(idx_eri_0_psi + 21);

    auto g_x_0_yyyyyz_0 = pbuffer.data(idx_eri_0_psi + 22);

    auto g_x_0_yyyyzz_0 = pbuffer.data(idx_eri_0_psi + 23);

    auto g_x_0_yyyzzz_0 = pbuffer.data(idx_eri_0_psi + 24);

    auto g_x_0_yyzzzz_0 = pbuffer.data(idx_eri_0_psi + 25);

    auto g_x_0_yzzzzz_0 = pbuffer.data(idx_eri_0_psi + 26);

    auto g_x_0_zzzzzz_0 = pbuffer.data(idx_eri_0_psi + 27);

    auto g_y_0_xxxxxx_0 = pbuffer.data(idx_eri_0_psi + 28);

    auto g_y_0_xxxxxy_0 = pbuffer.data(idx_eri_0_psi + 29);

    auto g_y_0_xxxxxz_0 = pbuffer.data(idx_eri_0_psi + 30);

    auto g_y_0_xxxxyy_0 = pbuffer.data(idx_eri_0_psi + 31);

    auto g_y_0_xxxxyz_0 = pbuffer.data(idx_eri_0_psi + 32);

    auto g_y_0_xxxxzz_0 = pbuffer.data(idx_eri_0_psi + 33);

    auto g_y_0_xxxyyy_0 = pbuffer.data(idx_eri_0_psi + 34);

    auto g_y_0_xxxyyz_0 = pbuffer.data(idx_eri_0_psi + 35);

    auto g_y_0_xxxyzz_0 = pbuffer.data(idx_eri_0_psi + 36);

    auto g_y_0_xxxzzz_0 = pbuffer.data(idx_eri_0_psi + 37);

    auto g_y_0_xxyyyy_0 = pbuffer.data(idx_eri_0_psi + 38);

    auto g_y_0_xxyyyz_0 = pbuffer.data(idx_eri_0_psi + 39);

    auto g_y_0_xxyyzz_0 = pbuffer.data(idx_eri_0_psi + 40);

    auto g_y_0_xxyzzz_0 = pbuffer.data(idx_eri_0_psi + 41);

    auto g_y_0_xxzzzz_0 = pbuffer.data(idx_eri_0_psi + 42);

    auto g_y_0_xyyyyy_0 = pbuffer.data(idx_eri_0_psi + 43);

    auto g_y_0_xyyyyz_0 = pbuffer.data(idx_eri_0_psi + 44);

    auto g_y_0_xyyyzz_0 = pbuffer.data(idx_eri_0_psi + 45);

    auto g_y_0_xyyzzz_0 = pbuffer.data(idx_eri_0_psi + 46);

    auto g_y_0_xyzzzz_0 = pbuffer.data(idx_eri_0_psi + 47);

    auto g_y_0_xzzzzz_0 = pbuffer.data(idx_eri_0_psi + 48);

    auto g_y_0_yyyyyy_0 = pbuffer.data(idx_eri_0_psi + 49);

    auto g_y_0_yyyyyz_0 = pbuffer.data(idx_eri_0_psi + 50);

    auto g_y_0_yyyyzz_0 = pbuffer.data(idx_eri_0_psi + 51);

    auto g_y_0_yyyzzz_0 = pbuffer.data(idx_eri_0_psi + 52);

    auto g_y_0_yyzzzz_0 = pbuffer.data(idx_eri_0_psi + 53);

    auto g_y_0_yzzzzz_0 = pbuffer.data(idx_eri_0_psi + 54);

    auto g_y_0_zzzzzz_0 = pbuffer.data(idx_eri_0_psi + 55);

    auto g_z_0_xxxxxx_0 = pbuffer.data(idx_eri_0_psi + 56);

    auto g_z_0_xxxxxy_0 = pbuffer.data(idx_eri_0_psi + 57);

    auto g_z_0_xxxxxz_0 = pbuffer.data(idx_eri_0_psi + 58);

    auto g_z_0_xxxxyy_0 = pbuffer.data(idx_eri_0_psi + 59);

    auto g_z_0_xxxxyz_0 = pbuffer.data(idx_eri_0_psi + 60);

    auto g_z_0_xxxxzz_0 = pbuffer.data(idx_eri_0_psi + 61);

    auto g_z_0_xxxyyy_0 = pbuffer.data(idx_eri_0_psi + 62);

    auto g_z_0_xxxyyz_0 = pbuffer.data(idx_eri_0_psi + 63);

    auto g_z_0_xxxyzz_0 = pbuffer.data(idx_eri_0_psi + 64);

    auto g_z_0_xxxzzz_0 = pbuffer.data(idx_eri_0_psi + 65);

    auto g_z_0_xxyyyy_0 = pbuffer.data(idx_eri_0_psi + 66);

    auto g_z_0_xxyyyz_0 = pbuffer.data(idx_eri_0_psi + 67);

    auto g_z_0_xxyyzz_0 = pbuffer.data(idx_eri_0_psi + 68);

    auto g_z_0_xxyzzz_0 = pbuffer.data(idx_eri_0_psi + 69);

    auto g_z_0_xxzzzz_0 = pbuffer.data(idx_eri_0_psi + 70);

    auto g_z_0_xyyyyy_0 = pbuffer.data(idx_eri_0_psi + 71);

    auto g_z_0_xyyyyz_0 = pbuffer.data(idx_eri_0_psi + 72);

    auto g_z_0_xyyyzz_0 = pbuffer.data(idx_eri_0_psi + 73);

    auto g_z_0_xyyzzz_0 = pbuffer.data(idx_eri_0_psi + 74);

    auto g_z_0_xyzzzz_0 = pbuffer.data(idx_eri_0_psi + 75);

    auto g_z_0_xzzzzz_0 = pbuffer.data(idx_eri_0_psi + 76);

    auto g_z_0_yyyyyy_0 = pbuffer.data(idx_eri_0_psi + 77);

    auto g_z_0_yyyyyz_0 = pbuffer.data(idx_eri_0_psi + 78);

    auto g_z_0_yyyyzz_0 = pbuffer.data(idx_eri_0_psi + 79);

    auto g_z_0_yyyzzz_0 = pbuffer.data(idx_eri_0_psi + 80);

    auto g_z_0_yyzzzz_0 = pbuffer.data(idx_eri_0_psi + 81);

    auto g_z_0_yzzzzz_0 = pbuffer.data(idx_eri_0_psi + 82);

    auto g_z_0_zzzzzz_0 = pbuffer.data(idx_eri_0_psi + 83);

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

    /// Set up components of auxilary buffer : DSH

    auto g_xx_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh);

    auto g_xx_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 1);

    auto g_xx_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 2);

    auto g_xx_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 3);

    auto g_xx_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 4);

    auto g_xx_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 5);

    auto g_xx_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 6);

    auto g_xx_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 7);

    auto g_xx_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 8);

    auto g_xx_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 9);

    auto g_xx_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 10);

    auto g_xx_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 11);

    auto g_xx_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 12);

    auto g_xx_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 13);

    auto g_xx_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 14);

    auto g_xx_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 15);

    auto g_xx_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 16);

    auto g_xx_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 17);

    auto g_xx_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 18);

    auto g_xx_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 19);

    auto g_xx_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 20);

    auto g_yy_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh + 63);

    auto g_yy_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 64);

    auto g_yy_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 65);

    auto g_yy_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 66);

    auto g_yy_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 67);

    auto g_yy_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 68);

    auto g_yy_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 69);

    auto g_yy_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 70);

    auto g_yy_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 71);

    auto g_yy_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 72);

    auto g_yy_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 73);

    auto g_yy_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 74);

    auto g_yy_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 75);

    auto g_yy_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 76);

    auto g_yy_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 77);

    auto g_yy_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 78);

    auto g_yy_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 79);

    auto g_yy_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 80);

    auto g_yy_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 81);

    auto g_yy_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 82);

    auto g_yy_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 83);

    auto g_yz_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 88);

    auto g_yz_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 91);

    auto g_yz_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 92);

    auto g_yz_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 95);

    auto g_yz_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 96);

    auto g_yz_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 97);

    auto g_yz_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 100);

    auto g_yz_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 101);

    auto g_yz_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 102);

    auto g_yz_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 103);

    auto g_zz_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh + 105);

    auto g_zz_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 106);

    auto g_zz_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 107);

    auto g_zz_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 108);

    auto g_zz_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 109);

    auto g_zz_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 110);

    auto g_zz_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 111);

    auto g_zz_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 112);

    auto g_zz_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 113);

    auto g_zz_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 114);

    auto g_zz_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 115);

    auto g_zz_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 116);

    auto g_zz_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 117);

    auto g_zz_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 118);

    auto g_zz_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 119);

    auto g_zz_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 120);

    auto g_zz_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 121);

    auto g_zz_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 122);

    auto g_zz_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 123);

    auto g_zz_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 124);

    auto g_zz_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 125);

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

    auto g_xy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_dsi + 29);

    auto g_xy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_dsi + 31);

    auto g_xy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_dsi + 34);

    auto g_xy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_dsi + 38);

    auto g_xy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 43);

    auto g_xz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_dsi + 56);

    auto g_xz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_dsi + 58);

    auto g_xz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_dsi + 61);

    auto g_xz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_dsi + 65);

    auto g_xz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_dsi + 70);

    auto g_xz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 76);

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

    auto g_yz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 133);

    auto g_yz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 134);

    auto g_yz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 135);

    auto g_yz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 136);

    auto g_yz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 137);

    auto g_yz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 138);

    auto g_yz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 139);

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

    /// Set up 0-28 components of targeted buffer : FSI

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

    #pragma omp simd aligned(g_x_0_xxxxxx_0, g_x_0_xxxxxx_1, g_x_0_xxxxxy_0, g_x_0_xxxxxy_1, g_x_0_xxxxxz_0, g_x_0_xxxxxz_1, g_x_0_xxxxyy_0, g_x_0_xxxxyy_1, g_x_0_xxxxyz_0, g_x_0_xxxxyz_1, g_x_0_xxxxzz_0, g_x_0_xxxxzz_1, g_x_0_xxxyyy_0, g_x_0_xxxyyy_1, g_x_0_xxxyyz_0, g_x_0_xxxyyz_1, g_x_0_xxxyzz_0, g_x_0_xxxyzz_1, g_x_0_xxxzzz_0, g_x_0_xxxzzz_1, g_x_0_xxyyyy_0, g_x_0_xxyyyy_1, g_x_0_xxyyyz_0, g_x_0_xxyyyz_1, g_x_0_xxyyzz_0, g_x_0_xxyyzz_1, g_x_0_xxyzzz_0, g_x_0_xxyzzz_1, g_x_0_xxzzzz_0, g_x_0_xxzzzz_1, g_x_0_xyyyyy_0, g_x_0_xyyyyy_1, g_x_0_xyyyyz_0, g_x_0_xyyyyz_1, g_x_0_xyyyzz_0, g_x_0_xyyyzz_1, g_x_0_xyyzzz_0, g_x_0_xyyzzz_1, g_x_0_xyzzzz_0, g_x_0_xyzzzz_1, g_x_0_xzzzzz_0, g_x_0_xzzzzz_1, g_x_0_yyyyyy_0, g_x_0_yyyyyy_1, g_x_0_yyyyyz_0, g_x_0_yyyyyz_1, g_x_0_yyyyzz_0, g_x_0_yyyyzz_1, g_x_0_yyyzzz_0, g_x_0_yyyzzz_1, g_x_0_yyzzzz_0, g_x_0_yyzzzz_1, g_x_0_yzzzzz_0, g_x_0_yzzzzz_1, g_x_0_zzzzzz_0, g_x_0_zzzzzz_1, g_xx_0_xxxxx_1, g_xx_0_xxxxxx_1, g_xx_0_xxxxxy_1, g_xx_0_xxxxxz_1, g_xx_0_xxxxy_1, g_xx_0_xxxxyy_1, g_xx_0_xxxxyz_1, g_xx_0_xxxxz_1, g_xx_0_xxxxzz_1, g_xx_0_xxxyy_1, g_xx_0_xxxyyy_1, g_xx_0_xxxyyz_1, g_xx_0_xxxyz_1, g_xx_0_xxxyzz_1, g_xx_0_xxxzz_1, g_xx_0_xxxzzz_1, g_xx_0_xxyyy_1, g_xx_0_xxyyyy_1, g_xx_0_xxyyyz_1, g_xx_0_xxyyz_1, g_xx_0_xxyyzz_1, g_xx_0_xxyzz_1, g_xx_0_xxyzzz_1, g_xx_0_xxzzz_1, g_xx_0_xxzzzz_1, g_xx_0_xyyyy_1, g_xx_0_xyyyyy_1, g_xx_0_xyyyyz_1, g_xx_0_xyyyz_1, g_xx_0_xyyyzz_1, g_xx_0_xyyzz_1, g_xx_0_xyyzzz_1, g_xx_0_xyzzz_1, g_xx_0_xyzzzz_1, g_xx_0_xzzzz_1, g_xx_0_xzzzzz_1, g_xx_0_yyyyy_1, g_xx_0_yyyyyy_1, g_xx_0_yyyyyz_1, g_xx_0_yyyyz_1, g_xx_0_yyyyzz_1, g_xx_0_yyyzz_1, g_xx_0_yyyzzz_1, g_xx_0_yyzzz_1, g_xx_0_yyzzzz_1, g_xx_0_yzzzz_1, g_xx_0_yzzzzz_1, g_xx_0_zzzzz_1, g_xx_0_zzzzzz_1, g_xxx_0_xxxxxx_0, g_xxx_0_xxxxxy_0, g_xxx_0_xxxxxz_0, g_xxx_0_xxxxyy_0, g_xxx_0_xxxxyz_0, g_xxx_0_xxxxzz_0, g_xxx_0_xxxyyy_0, g_xxx_0_xxxyyz_0, g_xxx_0_xxxyzz_0, g_xxx_0_xxxzzz_0, g_xxx_0_xxyyyy_0, g_xxx_0_xxyyyz_0, g_xxx_0_xxyyzz_0, g_xxx_0_xxyzzz_0, g_xxx_0_xxzzzz_0, g_xxx_0_xyyyyy_0, g_xxx_0_xyyyyz_0, g_xxx_0_xyyyzz_0, g_xxx_0_xyyzzz_0, g_xxx_0_xyzzzz_0, g_xxx_0_xzzzzz_0, g_xxx_0_yyyyyy_0, g_xxx_0_yyyyyz_0, g_xxx_0_yyyyzz_0, g_xxx_0_yyyzzz_0, g_xxx_0_yyzzzz_0, g_xxx_0_yzzzzz_0, g_xxx_0_zzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xxxxxx_0[i] = 2.0 * g_x_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxx_1[i] * wa_x[i];

        g_xxx_0_xxxxxy_0[i] = 2.0 * g_x_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxy_1[i] * wa_x[i];

        g_xxx_0_xxxxxz_0[i] = 2.0 * g_x_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxz_1[i] * wa_x[i];

        g_xxx_0_xxxxyy_0[i] = 2.0 * g_x_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxyy_1[i] * wa_x[i];

        g_xxx_0_xxxxyz_0[i] = 2.0 * g_x_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxyz_1[i] * wa_x[i];

        g_xxx_0_xxxxzz_0[i] = 2.0 * g_x_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxzz_1[i] * wa_x[i];

        g_xxx_0_xxxyyy_0[i] = 2.0 * g_x_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyy_1[i] * wa_x[i];

        g_xxx_0_xxxyyz_0[i] = 2.0 * g_x_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyz_1[i] * wa_x[i];

        g_xxx_0_xxxyzz_0[i] = 2.0 * g_x_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyzz_1[i] * fi_acd_0 + g_xx_0_xxxyzz_1[i] * wa_x[i];

        g_xxx_0_xxxzzz_0[i] = 2.0 * g_x_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxzzz_1[i] * fi_acd_0 + g_xx_0_xxxzzz_1[i] * wa_x[i];

        g_xxx_0_xxyyyy_0[i] = 2.0 * g_x_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyy_1[i] * wa_x[i];

        g_xxx_0_xxyyyz_0[i] = 2.0 * g_x_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyz_1[i] * wa_x[i];

        g_xxx_0_xxyyzz_0[i] = 2.0 * g_x_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyzz_1[i] * wa_x[i];

        g_xxx_0_xxyzzz_0[i] = 2.0 * g_x_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzz_1[i] * wa_x[i];

        g_xxx_0_xxzzzz_0[i] = 2.0 * g_x_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xzzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyyy_0[i] = 2.0 * g_x_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyy_1[i] * wa_x[i];

        g_xxx_0_xyyyyz_0[i] = 2.0 * g_x_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyz_1[i] * wa_x[i];

        g_xxx_0_xyyyzz_0[i] = 2.0 * g_x_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyzz_1[i] * fz_be_0 + g_xx_0_yyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyzz_1[i] * wa_x[i];

        g_xxx_0_xyyzzz_0[i] = 2.0 * g_x_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyzzz_1[i] * fz_be_0 + g_xx_0_yyzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzz_1[i] * wa_x[i];

        g_xxx_0_xyzzzz_0[i] = 2.0 * g_x_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyzzzz_1[i] * fz_be_0 + g_xx_0_yzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzz_1[i] * wa_x[i];

        g_xxx_0_xzzzzz_0[i] = 2.0 * g_x_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyyy_0[i] = 2.0 * g_x_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyyy_1[i] * wa_x[i];

        g_xxx_0_yyyyyz_0[i] = 2.0 * g_x_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyyz_1[i] * wa_x[i];

        g_xxx_0_yyyyzz_0[i] = 2.0 * g_x_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyzz_1[i] * fz_be_0 + g_xx_0_yyyyzz_1[i] * wa_x[i];

        g_xxx_0_yyyzzz_0[i] = 2.0 * g_x_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyzzz_1[i] * fz_be_0 + g_xx_0_yyyzzz_1[i] * wa_x[i];

        g_xxx_0_yyzzzz_0[i] = 2.0 * g_x_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyzzzz_1[i] * fz_be_0 + g_xx_0_yyzzzz_1[i] * wa_x[i];

        g_xxx_0_yzzzzz_0[i] = 2.0 * g_x_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yzzzzz_1[i] * fz_be_0 + g_xx_0_yzzzzz_1[i] * wa_x[i];

        g_xxx_0_zzzzzz_0[i] = 2.0 * g_x_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_zzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 28-56 components of targeted buffer : FSI

    auto g_xxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 28);

    auto g_xxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 29);

    auto g_xxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 30);

    auto g_xxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 31);

    auto g_xxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 32);

    auto g_xxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 33);

    auto g_xxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 34);

    auto g_xxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 35);

    auto g_xxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 36);

    auto g_xxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 37);

    auto g_xxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 38);

    auto g_xxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 39);

    auto g_xxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 40);

    auto g_xxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 41);

    auto g_xxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 42);

    auto g_xxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 43);

    auto g_xxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 44);

    auto g_xxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 45);

    auto g_xxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 46);

    auto g_xxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 47);

    auto g_xxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 48);

    auto g_xxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 49);

    auto g_xxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 50);

    auto g_xxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 51);

    auto g_xxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 52);

    auto g_xxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 53);

    auto g_xxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 54);

    auto g_xxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 55);

    #pragma omp simd aligned(g_xx_0_xxxxx_1, g_xx_0_xxxxxx_1, g_xx_0_xxxxxy_1, g_xx_0_xxxxxz_1, g_xx_0_xxxxy_1, g_xx_0_xxxxyy_1, g_xx_0_xxxxyz_1, g_xx_0_xxxxz_1, g_xx_0_xxxxzz_1, g_xx_0_xxxyy_1, g_xx_0_xxxyyy_1, g_xx_0_xxxyyz_1, g_xx_0_xxxyz_1, g_xx_0_xxxyzz_1, g_xx_0_xxxzz_1, g_xx_0_xxxzzz_1, g_xx_0_xxyyy_1, g_xx_0_xxyyyy_1, g_xx_0_xxyyyz_1, g_xx_0_xxyyz_1, g_xx_0_xxyyzz_1, g_xx_0_xxyzz_1, g_xx_0_xxyzzz_1, g_xx_0_xxzzz_1, g_xx_0_xxzzzz_1, g_xx_0_xyyyy_1, g_xx_0_xyyyyy_1, g_xx_0_xyyyyz_1, g_xx_0_xyyyz_1, g_xx_0_xyyyzz_1, g_xx_0_xyyzz_1, g_xx_0_xyyzzz_1, g_xx_0_xyzzz_1, g_xx_0_xyzzzz_1, g_xx_0_xzzzz_1, g_xx_0_xzzzzz_1, g_xx_0_yyyyy_1, g_xx_0_yyyyyy_1, g_xx_0_yyyyyz_1, g_xx_0_yyyyz_1, g_xx_0_yyyyzz_1, g_xx_0_yyyzz_1, g_xx_0_yyyzzz_1, g_xx_0_yyzzz_1, g_xx_0_yyzzzz_1, g_xx_0_yzzzz_1, g_xx_0_yzzzzz_1, g_xx_0_zzzzz_1, g_xx_0_zzzzzz_1, g_xxy_0_xxxxxx_0, g_xxy_0_xxxxxy_0, g_xxy_0_xxxxxz_0, g_xxy_0_xxxxyy_0, g_xxy_0_xxxxyz_0, g_xxy_0_xxxxzz_0, g_xxy_0_xxxyyy_0, g_xxy_0_xxxyyz_0, g_xxy_0_xxxyzz_0, g_xxy_0_xxxzzz_0, g_xxy_0_xxyyyy_0, g_xxy_0_xxyyyz_0, g_xxy_0_xxyyzz_0, g_xxy_0_xxyzzz_0, g_xxy_0_xxzzzz_0, g_xxy_0_xyyyyy_0, g_xxy_0_xyyyyz_0, g_xxy_0_xyyyzz_0, g_xxy_0_xyyzzz_0, g_xxy_0_xyzzzz_0, g_xxy_0_xzzzzz_0, g_xxy_0_yyyyyy_0, g_xxy_0_yyyyyz_0, g_xxy_0_yyyyzz_0, g_xxy_0_yyyzzz_0, g_xxy_0_yyzzzz_0, g_xxy_0_yzzzzz_0, g_xxy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xxxxxx_0[i] = g_xx_0_xxxxxx_1[i] * wa_y[i];

        g_xxy_0_xxxxxy_0[i] = g_xx_0_xxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxy_1[i] * wa_y[i];

        g_xxy_0_xxxxxz_0[i] = g_xx_0_xxxxxz_1[i] * wa_y[i];

        g_xxy_0_xxxxyy_0[i] = 2.0 * g_xx_0_xxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxyy_1[i] * wa_y[i];

        g_xxy_0_xxxxyz_0[i] = g_xx_0_xxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxyz_1[i] * wa_y[i];

        g_xxy_0_xxxxzz_0[i] = g_xx_0_xxxxzz_1[i] * wa_y[i];

        g_xxy_0_xxxyyy_0[i] = 3.0 * g_xx_0_xxxyy_1[i] * fi_acd_0 + g_xx_0_xxxyyy_1[i] * wa_y[i];

        g_xxy_0_xxxyyz_0[i] = 2.0 * g_xx_0_xxxyz_1[i] * fi_acd_0 + g_xx_0_xxxyyz_1[i] * wa_y[i];

        g_xxy_0_xxxyzz_0[i] = g_xx_0_xxxzz_1[i] * fi_acd_0 + g_xx_0_xxxyzz_1[i] * wa_y[i];

        g_xxy_0_xxxzzz_0[i] = g_xx_0_xxxzzz_1[i] * wa_y[i];

        g_xxy_0_xxyyyy_0[i] = 4.0 * g_xx_0_xxyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyy_1[i] * wa_y[i];

        g_xxy_0_xxyyyz_0[i] = 3.0 * g_xx_0_xxyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyz_1[i] * wa_y[i];

        g_xxy_0_xxyyzz_0[i] = 2.0 * g_xx_0_xxyzz_1[i] * fi_acd_0 + g_xx_0_xxyyzz_1[i] * wa_y[i];

        g_xxy_0_xxyzzz_0[i] = g_xx_0_xxzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzz_1[i] * wa_y[i];

        g_xxy_0_xxzzzz_0[i] = g_xx_0_xxzzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyyy_0[i] = 5.0 * g_xx_0_xyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyy_1[i] * wa_y[i];

        g_xxy_0_xyyyyz_0[i] = 4.0 * g_xx_0_xyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyz_1[i] * wa_y[i];

        g_xxy_0_xyyyzz_0[i] = 3.0 * g_xx_0_xyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyzz_1[i] * wa_y[i];

        g_xxy_0_xyyzzz_0[i] = 2.0 * g_xx_0_xyzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzz_1[i] * wa_y[i];

        g_xxy_0_xyzzzz_0[i] = g_xx_0_xzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzz_1[i] * wa_y[i];

        g_xxy_0_xzzzzz_0[i] = g_xx_0_xzzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyyy_0[i] = 6.0 * g_xx_0_yyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyy_1[i] * wa_y[i];

        g_xxy_0_yyyyyz_0[i] = 5.0 * g_xx_0_yyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyyz_1[i] * wa_y[i];

        g_xxy_0_yyyyzz_0[i] = 4.0 * g_xx_0_yyyzz_1[i] * fi_acd_0 + g_xx_0_yyyyzz_1[i] * wa_y[i];

        g_xxy_0_yyyzzz_0[i] = 3.0 * g_xx_0_yyzzz_1[i] * fi_acd_0 + g_xx_0_yyyzzz_1[i] * wa_y[i];

        g_xxy_0_yyzzzz_0[i] = 2.0 * g_xx_0_yzzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzz_1[i] * wa_y[i];

        g_xxy_0_yzzzzz_0[i] = g_xx_0_zzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzz_1[i] * wa_y[i];

        g_xxy_0_zzzzzz_0[i] = g_xx_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 56-84 components of targeted buffer : FSI

    auto g_xxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 56);

    auto g_xxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 57);

    auto g_xxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 58);

    auto g_xxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 59);

    auto g_xxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 60);

    auto g_xxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 61);

    auto g_xxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 62);

    auto g_xxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 63);

    auto g_xxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 64);

    auto g_xxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 65);

    auto g_xxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 66);

    auto g_xxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 67);

    auto g_xxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 68);

    auto g_xxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 69);

    auto g_xxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 70);

    auto g_xxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 71);

    auto g_xxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 72);

    auto g_xxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 73);

    auto g_xxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 74);

    auto g_xxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 75);

    auto g_xxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 76);

    auto g_xxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 77);

    auto g_xxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 78);

    auto g_xxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 79);

    auto g_xxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 80);

    auto g_xxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 81);

    auto g_xxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 82);

    auto g_xxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 83);

    #pragma omp simd aligned(g_xx_0_xxxxx_1, g_xx_0_xxxxxx_1, g_xx_0_xxxxxy_1, g_xx_0_xxxxxz_1, g_xx_0_xxxxy_1, g_xx_0_xxxxyy_1, g_xx_0_xxxxyz_1, g_xx_0_xxxxz_1, g_xx_0_xxxxzz_1, g_xx_0_xxxyy_1, g_xx_0_xxxyyy_1, g_xx_0_xxxyyz_1, g_xx_0_xxxyz_1, g_xx_0_xxxyzz_1, g_xx_0_xxxzz_1, g_xx_0_xxxzzz_1, g_xx_0_xxyyy_1, g_xx_0_xxyyyy_1, g_xx_0_xxyyyz_1, g_xx_0_xxyyz_1, g_xx_0_xxyyzz_1, g_xx_0_xxyzz_1, g_xx_0_xxyzzz_1, g_xx_0_xxzzz_1, g_xx_0_xxzzzz_1, g_xx_0_xyyyy_1, g_xx_0_xyyyyy_1, g_xx_0_xyyyyz_1, g_xx_0_xyyyz_1, g_xx_0_xyyyzz_1, g_xx_0_xyyzz_1, g_xx_0_xyyzzz_1, g_xx_0_xyzzz_1, g_xx_0_xyzzzz_1, g_xx_0_xzzzz_1, g_xx_0_xzzzzz_1, g_xx_0_yyyyy_1, g_xx_0_yyyyyy_1, g_xx_0_yyyyyz_1, g_xx_0_yyyyz_1, g_xx_0_yyyyzz_1, g_xx_0_yyyzz_1, g_xx_0_yyyzzz_1, g_xx_0_yyzzz_1, g_xx_0_yyzzzz_1, g_xx_0_yzzzz_1, g_xx_0_yzzzzz_1, g_xx_0_zzzzz_1, g_xx_0_zzzzzz_1, g_xxz_0_xxxxxx_0, g_xxz_0_xxxxxy_0, g_xxz_0_xxxxxz_0, g_xxz_0_xxxxyy_0, g_xxz_0_xxxxyz_0, g_xxz_0_xxxxzz_0, g_xxz_0_xxxyyy_0, g_xxz_0_xxxyyz_0, g_xxz_0_xxxyzz_0, g_xxz_0_xxxzzz_0, g_xxz_0_xxyyyy_0, g_xxz_0_xxyyyz_0, g_xxz_0_xxyyzz_0, g_xxz_0_xxyzzz_0, g_xxz_0_xxzzzz_0, g_xxz_0_xyyyyy_0, g_xxz_0_xyyyyz_0, g_xxz_0_xyyyzz_0, g_xxz_0_xyyzzz_0, g_xxz_0_xyzzzz_0, g_xxz_0_xzzzzz_0, g_xxz_0_yyyyyy_0, g_xxz_0_yyyyyz_0, g_xxz_0_yyyyzz_0, g_xxz_0_yyyzzz_0, g_xxz_0_yyzzzz_0, g_xxz_0_yzzzzz_0, g_xxz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xxxxxx_0[i] = g_xx_0_xxxxxx_1[i] * wa_z[i];

        g_xxz_0_xxxxxy_0[i] = g_xx_0_xxxxxy_1[i] * wa_z[i];

        g_xxz_0_xxxxxz_0[i] = g_xx_0_xxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxz_1[i] * wa_z[i];

        g_xxz_0_xxxxyy_0[i] = g_xx_0_xxxxyy_1[i] * wa_z[i];

        g_xxz_0_xxxxyz_0[i] = g_xx_0_xxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxyz_1[i] * wa_z[i];

        g_xxz_0_xxxxzz_0[i] = 2.0 * g_xx_0_xxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxzz_1[i] * wa_z[i];

        g_xxz_0_xxxyyy_0[i] = g_xx_0_xxxyyy_1[i] * wa_z[i];

        g_xxz_0_xxxyyz_0[i] = g_xx_0_xxxyy_1[i] * fi_acd_0 + g_xx_0_xxxyyz_1[i] * wa_z[i];

        g_xxz_0_xxxyzz_0[i] = 2.0 * g_xx_0_xxxyz_1[i] * fi_acd_0 + g_xx_0_xxxyzz_1[i] * wa_z[i];

        g_xxz_0_xxxzzz_0[i] = 3.0 * g_xx_0_xxxzz_1[i] * fi_acd_0 + g_xx_0_xxxzzz_1[i] * wa_z[i];

        g_xxz_0_xxyyyy_0[i] = g_xx_0_xxyyyy_1[i] * wa_z[i];

        g_xxz_0_xxyyyz_0[i] = g_xx_0_xxyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyz_1[i] * wa_z[i];

        g_xxz_0_xxyyzz_0[i] = 2.0 * g_xx_0_xxyyz_1[i] * fi_acd_0 + g_xx_0_xxyyzz_1[i] * wa_z[i];

        g_xxz_0_xxyzzz_0[i] = 3.0 * g_xx_0_xxyzz_1[i] * fi_acd_0 + g_xx_0_xxyzzz_1[i] * wa_z[i];

        g_xxz_0_xxzzzz_0[i] = 4.0 * g_xx_0_xxzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyyy_0[i] = g_xx_0_xyyyyy_1[i] * wa_z[i];

        g_xxz_0_xyyyyz_0[i] = g_xx_0_xyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyz_1[i] * wa_z[i];

        g_xxz_0_xyyyzz_0[i] = 2.0 * g_xx_0_xyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyzz_1[i] * wa_z[i];

        g_xxz_0_xyyzzz_0[i] = 3.0 * g_xx_0_xyyzz_1[i] * fi_acd_0 + g_xx_0_xyyzzz_1[i] * wa_z[i];

        g_xxz_0_xyzzzz_0[i] = 4.0 * g_xx_0_xyzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzz_1[i] * wa_z[i];

        g_xxz_0_xzzzzz_0[i] = 5.0 * g_xx_0_xzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyyy_0[i] = g_xx_0_yyyyyy_1[i] * wa_z[i];

        g_xxz_0_yyyyyz_0[i] = g_xx_0_yyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyz_1[i] * wa_z[i];

        g_xxz_0_yyyyzz_0[i] = 2.0 * g_xx_0_yyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyzz_1[i] * wa_z[i];

        g_xxz_0_yyyzzz_0[i] = 3.0 * g_xx_0_yyyzz_1[i] * fi_acd_0 + g_xx_0_yyyzzz_1[i] * wa_z[i];

        g_xxz_0_yyzzzz_0[i] = 4.0 * g_xx_0_yyzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzz_1[i] * wa_z[i];

        g_xxz_0_yzzzzz_0[i] = 5.0 * g_xx_0_yzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzz_1[i] * wa_z[i];

        g_xxz_0_zzzzzz_0[i] = 6.0 * g_xx_0_zzzzz_1[i] * fi_acd_0 + g_xx_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 84-112 components of targeted buffer : FSI

    auto g_xyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 84);

    auto g_xyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 85);

    auto g_xyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 86);

    auto g_xyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 87);

    auto g_xyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 88);

    auto g_xyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 89);

    auto g_xyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 90);

    auto g_xyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 91);

    auto g_xyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 92);

    auto g_xyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 93);

    auto g_xyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 94);

    auto g_xyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 95);

    auto g_xyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 96);

    auto g_xyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 97);

    auto g_xyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 98);

    auto g_xyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 99);

    auto g_xyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 100);

    auto g_xyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 101);

    auto g_xyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 102);

    auto g_xyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 103);

    auto g_xyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 104);

    auto g_xyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 105);

    auto g_xyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 106);

    auto g_xyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 107);

    auto g_xyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 108);

    auto g_xyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 109);

    auto g_xyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 110);

    auto g_xyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 111);

    #pragma omp simd aligned(g_xyy_0_xxxxxx_0, g_xyy_0_xxxxxy_0, g_xyy_0_xxxxxz_0, g_xyy_0_xxxxyy_0, g_xyy_0_xxxxyz_0, g_xyy_0_xxxxzz_0, g_xyy_0_xxxyyy_0, g_xyy_0_xxxyyz_0, g_xyy_0_xxxyzz_0, g_xyy_0_xxxzzz_0, g_xyy_0_xxyyyy_0, g_xyy_0_xxyyyz_0, g_xyy_0_xxyyzz_0, g_xyy_0_xxyzzz_0, g_xyy_0_xxzzzz_0, g_xyy_0_xyyyyy_0, g_xyy_0_xyyyyz_0, g_xyy_0_xyyyzz_0, g_xyy_0_xyyzzz_0, g_xyy_0_xyzzzz_0, g_xyy_0_xzzzzz_0, g_xyy_0_yyyyyy_0, g_xyy_0_yyyyyz_0, g_xyy_0_yyyyzz_0, g_xyy_0_yyyzzz_0, g_xyy_0_yyzzzz_0, g_xyy_0_yzzzzz_0, g_xyy_0_zzzzzz_0, g_yy_0_xxxxx_1, g_yy_0_xxxxxx_1, g_yy_0_xxxxxy_1, g_yy_0_xxxxxz_1, g_yy_0_xxxxy_1, g_yy_0_xxxxyy_1, g_yy_0_xxxxyz_1, g_yy_0_xxxxz_1, g_yy_0_xxxxzz_1, g_yy_0_xxxyy_1, g_yy_0_xxxyyy_1, g_yy_0_xxxyyz_1, g_yy_0_xxxyz_1, g_yy_0_xxxyzz_1, g_yy_0_xxxzz_1, g_yy_0_xxxzzz_1, g_yy_0_xxyyy_1, g_yy_0_xxyyyy_1, g_yy_0_xxyyyz_1, g_yy_0_xxyyz_1, g_yy_0_xxyyzz_1, g_yy_0_xxyzz_1, g_yy_0_xxyzzz_1, g_yy_0_xxzzz_1, g_yy_0_xxzzzz_1, g_yy_0_xyyyy_1, g_yy_0_xyyyyy_1, g_yy_0_xyyyyz_1, g_yy_0_xyyyz_1, g_yy_0_xyyyzz_1, g_yy_0_xyyzz_1, g_yy_0_xyyzzz_1, g_yy_0_xyzzz_1, g_yy_0_xyzzzz_1, g_yy_0_xzzzz_1, g_yy_0_xzzzzz_1, g_yy_0_yyyyy_1, g_yy_0_yyyyyy_1, g_yy_0_yyyyyz_1, g_yy_0_yyyyz_1, g_yy_0_yyyyzz_1, g_yy_0_yyyzz_1, g_yy_0_yyyzzz_1, g_yy_0_yyzzz_1, g_yy_0_yyzzzz_1, g_yy_0_yzzzz_1, g_yy_0_yzzzzz_1, g_yy_0_zzzzz_1, g_yy_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xxxxxx_0[i] = 6.0 * g_yy_0_xxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxx_1[i] * wa_x[i];

        g_xyy_0_xxxxxy_0[i] = 5.0 * g_yy_0_xxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxy_1[i] * wa_x[i];

        g_xyy_0_xxxxxz_0[i] = 5.0 * g_yy_0_xxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxz_1[i] * wa_x[i];

        g_xyy_0_xxxxyy_0[i] = 4.0 * g_yy_0_xxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxyy_1[i] * wa_x[i];

        g_xyy_0_xxxxyz_0[i] = 4.0 * g_yy_0_xxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxyz_1[i] * wa_x[i];

        g_xyy_0_xxxxzz_0[i] = 4.0 * g_yy_0_xxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxzz_1[i] * wa_x[i];

        g_xyy_0_xxxyyy_0[i] = 3.0 * g_yy_0_xxyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyy_1[i] * wa_x[i];

        g_xyy_0_xxxyyz_0[i] = 3.0 * g_yy_0_xxyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyz_1[i] * wa_x[i];

        g_xyy_0_xxxyzz_0[i] = 3.0 * g_yy_0_xxyzz_1[i] * fi_acd_0 + g_yy_0_xxxyzz_1[i] * wa_x[i];

        g_xyy_0_xxxzzz_0[i] = 3.0 * g_yy_0_xxzzz_1[i] * fi_acd_0 + g_yy_0_xxxzzz_1[i] * wa_x[i];

        g_xyy_0_xxyyyy_0[i] = 2.0 * g_yy_0_xyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyy_1[i] * wa_x[i];

        g_xyy_0_xxyyyz_0[i] = 2.0 * g_yy_0_xyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyz_1[i] * wa_x[i];

        g_xyy_0_xxyyzz_0[i] = 2.0 * g_yy_0_xyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyzz_1[i] * wa_x[i];

        g_xyy_0_xxyzzz_0[i] = 2.0 * g_yy_0_xyzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzz_1[i] * wa_x[i];

        g_xyy_0_xxzzzz_0[i] = 2.0 * g_yy_0_xzzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyyy_0[i] = g_yy_0_yyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyy_1[i] * wa_x[i];

        g_xyy_0_xyyyyz_0[i] = g_yy_0_yyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyz_1[i] * wa_x[i];

        g_xyy_0_xyyyzz_0[i] = g_yy_0_yyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyzz_1[i] * wa_x[i];

        g_xyy_0_xyyzzz_0[i] = g_yy_0_yyzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzz_1[i] * wa_x[i];

        g_xyy_0_xyzzzz_0[i] = g_yy_0_yzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzz_1[i] * wa_x[i];

        g_xyy_0_xzzzzz_0[i] = g_yy_0_zzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyyy_0[i] = g_yy_0_yyyyyy_1[i] * wa_x[i];

        g_xyy_0_yyyyyz_0[i] = g_yy_0_yyyyyz_1[i] * wa_x[i];

        g_xyy_0_yyyyzz_0[i] = g_yy_0_yyyyzz_1[i] * wa_x[i];

        g_xyy_0_yyyzzz_0[i] = g_yy_0_yyyzzz_1[i] * wa_x[i];

        g_xyy_0_yyzzzz_0[i] = g_yy_0_yyzzzz_1[i] * wa_x[i];

        g_xyy_0_yzzzzz_0[i] = g_yy_0_yzzzzz_1[i] * wa_x[i];

        g_xyy_0_zzzzzz_0[i] = g_yy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 112-140 components of targeted buffer : FSI

    auto g_xyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 112);

    auto g_xyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 113);

    auto g_xyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 114);

    auto g_xyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 115);

    auto g_xyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 116);

    auto g_xyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 117);

    auto g_xyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 118);

    auto g_xyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 119);

    auto g_xyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 120);

    auto g_xyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 121);

    auto g_xyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 122);

    auto g_xyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 123);

    auto g_xyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 124);

    auto g_xyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 125);

    auto g_xyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 126);

    auto g_xyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 127);

    auto g_xyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 128);

    auto g_xyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 129);

    auto g_xyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 130);

    auto g_xyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 131);

    auto g_xyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 132);

    auto g_xyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 133);

    auto g_xyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 134);

    auto g_xyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 135);

    auto g_xyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 136);

    auto g_xyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 137);

    auto g_xyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 138);

    auto g_xyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 139);

    #pragma omp simd aligned(g_xy_0_xxxxxy_1, g_xy_0_xxxxyy_1, g_xy_0_xxxyyy_1, g_xy_0_xxyyyy_1, g_xy_0_xyyyyy_1, g_xyz_0_xxxxxx_0, g_xyz_0_xxxxxy_0, g_xyz_0_xxxxxz_0, g_xyz_0_xxxxyy_0, g_xyz_0_xxxxyz_0, g_xyz_0_xxxxzz_0, g_xyz_0_xxxyyy_0, g_xyz_0_xxxyyz_0, g_xyz_0_xxxyzz_0, g_xyz_0_xxxzzz_0, g_xyz_0_xxyyyy_0, g_xyz_0_xxyyyz_0, g_xyz_0_xxyyzz_0, g_xyz_0_xxyzzz_0, g_xyz_0_xxzzzz_0, g_xyz_0_xyyyyy_0, g_xyz_0_xyyyyz_0, g_xyz_0_xyyyzz_0, g_xyz_0_xyyzzz_0, g_xyz_0_xyzzzz_0, g_xyz_0_xzzzzz_0, g_xyz_0_yyyyyy_0, g_xyz_0_yyyyyz_0, g_xyz_0_yyyyzz_0, g_xyz_0_yyyzzz_0, g_xyz_0_yyzzzz_0, g_xyz_0_yzzzzz_0, g_xyz_0_zzzzzz_0, g_xz_0_xxxxxx_1, g_xz_0_xxxxxz_1, g_xz_0_xxxxzz_1, g_xz_0_xxxzzz_1, g_xz_0_xxzzzz_1, g_xz_0_xzzzzz_1, g_yz_0_xxxxyz_1, g_yz_0_xxxyyz_1, g_yz_0_xxxyz_1, g_yz_0_xxxyzz_1, g_yz_0_xxyyyz_1, g_yz_0_xxyyz_1, g_yz_0_xxyyzz_1, g_yz_0_xxyzz_1, g_yz_0_xxyzzz_1, g_yz_0_xyyyyz_1, g_yz_0_xyyyz_1, g_yz_0_xyyyzz_1, g_yz_0_xyyzz_1, g_yz_0_xyyzzz_1, g_yz_0_xyzzz_1, g_yz_0_xyzzzz_1, g_yz_0_yyyyyy_1, g_yz_0_yyyyyz_1, g_yz_0_yyyyz_1, g_yz_0_yyyyzz_1, g_yz_0_yyyzz_1, g_yz_0_yyyzzz_1, g_yz_0_yyzzz_1, g_yz_0_yyzzzz_1, g_yz_0_yzzzz_1, g_yz_0_yzzzzz_1, g_yz_0_zzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyz_0_xxxxxx_0[i] = g_xz_0_xxxxxx_1[i] * wa_y[i];

        g_xyz_0_xxxxxy_0[i] = g_xy_0_xxxxxy_1[i] * wa_z[i];

        g_xyz_0_xxxxxz_0[i] = g_xz_0_xxxxxz_1[i] * wa_y[i];

        g_xyz_0_xxxxyy_0[i] = g_xy_0_xxxxyy_1[i] * wa_z[i];

        g_xyz_0_xxxxyz_0[i] = 4.0 * g_yz_0_xxxyz_1[i] * fi_acd_0 + g_yz_0_xxxxyz_1[i] * wa_x[i];

        g_xyz_0_xxxxzz_0[i] = g_xz_0_xxxxzz_1[i] * wa_y[i];

        g_xyz_0_xxxyyy_0[i] = g_xy_0_xxxyyy_1[i] * wa_z[i];

        g_xyz_0_xxxyyz_0[i] = 3.0 * g_yz_0_xxyyz_1[i] * fi_acd_0 + g_yz_0_xxxyyz_1[i] * wa_x[i];

        g_xyz_0_xxxyzz_0[i] = 3.0 * g_yz_0_xxyzz_1[i] * fi_acd_0 + g_yz_0_xxxyzz_1[i] * wa_x[i];

        g_xyz_0_xxxzzz_0[i] = g_xz_0_xxxzzz_1[i] * wa_y[i];

        g_xyz_0_xxyyyy_0[i] = g_xy_0_xxyyyy_1[i] * wa_z[i];

        g_xyz_0_xxyyyz_0[i] = 2.0 * g_yz_0_xyyyz_1[i] * fi_acd_0 + g_yz_0_xxyyyz_1[i] * wa_x[i];

        g_xyz_0_xxyyzz_0[i] = 2.0 * g_yz_0_xyyzz_1[i] * fi_acd_0 + g_yz_0_xxyyzz_1[i] * wa_x[i];

        g_xyz_0_xxyzzz_0[i] = 2.0 * g_yz_0_xyzzz_1[i] * fi_acd_0 + g_yz_0_xxyzzz_1[i] * wa_x[i];

        g_xyz_0_xxzzzz_0[i] = g_xz_0_xxzzzz_1[i] * wa_y[i];

        g_xyz_0_xyyyyy_0[i] = g_xy_0_xyyyyy_1[i] * wa_z[i];

        g_xyz_0_xyyyyz_0[i] = g_yz_0_yyyyz_1[i] * fi_acd_0 + g_yz_0_xyyyyz_1[i] * wa_x[i];

        g_xyz_0_xyyyzz_0[i] = g_yz_0_yyyzz_1[i] * fi_acd_0 + g_yz_0_xyyyzz_1[i] * wa_x[i];

        g_xyz_0_xyyzzz_0[i] = g_yz_0_yyzzz_1[i] * fi_acd_0 + g_yz_0_xyyzzz_1[i] * wa_x[i];

        g_xyz_0_xyzzzz_0[i] = g_yz_0_yzzzz_1[i] * fi_acd_0 + g_yz_0_xyzzzz_1[i] * wa_x[i];

        g_xyz_0_xzzzzz_0[i] = g_xz_0_xzzzzz_1[i] * wa_y[i];

        g_xyz_0_yyyyyy_0[i] = g_yz_0_yyyyyy_1[i] * wa_x[i];

        g_xyz_0_yyyyyz_0[i] = g_yz_0_yyyyyz_1[i] * wa_x[i];

        g_xyz_0_yyyyzz_0[i] = g_yz_0_yyyyzz_1[i] * wa_x[i];

        g_xyz_0_yyyzzz_0[i] = g_yz_0_yyyzzz_1[i] * wa_x[i];

        g_xyz_0_yyzzzz_0[i] = g_yz_0_yyzzzz_1[i] * wa_x[i];

        g_xyz_0_yzzzzz_0[i] = g_yz_0_yzzzzz_1[i] * wa_x[i];

        g_xyz_0_zzzzzz_0[i] = g_yz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 140-168 components of targeted buffer : FSI

    auto g_xzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 140);

    auto g_xzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 141);

    auto g_xzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 142);

    auto g_xzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 143);

    auto g_xzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 144);

    auto g_xzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 145);

    auto g_xzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 146);

    auto g_xzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 147);

    auto g_xzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 148);

    auto g_xzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 149);

    auto g_xzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 150);

    auto g_xzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 151);

    auto g_xzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 152);

    auto g_xzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 153);

    auto g_xzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 154);

    auto g_xzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 155);

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

    #pragma omp simd aligned(g_xzz_0_xxxxxx_0, g_xzz_0_xxxxxy_0, g_xzz_0_xxxxxz_0, g_xzz_0_xxxxyy_0, g_xzz_0_xxxxyz_0, g_xzz_0_xxxxzz_0, g_xzz_0_xxxyyy_0, g_xzz_0_xxxyyz_0, g_xzz_0_xxxyzz_0, g_xzz_0_xxxzzz_0, g_xzz_0_xxyyyy_0, g_xzz_0_xxyyyz_0, g_xzz_0_xxyyzz_0, g_xzz_0_xxyzzz_0, g_xzz_0_xxzzzz_0, g_xzz_0_xyyyyy_0, g_xzz_0_xyyyyz_0, g_xzz_0_xyyyzz_0, g_xzz_0_xyyzzz_0, g_xzz_0_xyzzzz_0, g_xzz_0_xzzzzz_0, g_xzz_0_yyyyyy_0, g_xzz_0_yyyyyz_0, g_xzz_0_yyyyzz_0, g_xzz_0_yyyzzz_0, g_xzz_0_yyzzzz_0, g_xzz_0_yzzzzz_0, g_xzz_0_zzzzzz_0, g_zz_0_xxxxx_1, g_zz_0_xxxxxx_1, g_zz_0_xxxxxy_1, g_zz_0_xxxxxz_1, g_zz_0_xxxxy_1, g_zz_0_xxxxyy_1, g_zz_0_xxxxyz_1, g_zz_0_xxxxz_1, g_zz_0_xxxxzz_1, g_zz_0_xxxyy_1, g_zz_0_xxxyyy_1, g_zz_0_xxxyyz_1, g_zz_0_xxxyz_1, g_zz_0_xxxyzz_1, g_zz_0_xxxzz_1, g_zz_0_xxxzzz_1, g_zz_0_xxyyy_1, g_zz_0_xxyyyy_1, g_zz_0_xxyyyz_1, g_zz_0_xxyyz_1, g_zz_0_xxyyzz_1, g_zz_0_xxyzz_1, g_zz_0_xxyzzz_1, g_zz_0_xxzzz_1, g_zz_0_xxzzzz_1, g_zz_0_xyyyy_1, g_zz_0_xyyyyy_1, g_zz_0_xyyyyz_1, g_zz_0_xyyyz_1, g_zz_0_xyyyzz_1, g_zz_0_xyyzz_1, g_zz_0_xyyzzz_1, g_zz_0_xyzzz_1, g_zz_0_xyzzzz_1, g_zz_0_xzzzz_1, g_zz_0_xzzzzz_1, g_zz_0_yyyyy_1, g_zz_0_yyyyyy_1, g_zz_0_yyyyyz_1, g_zz_0_yyyyz_1, g_zz_0_yyyyzz_1, g_zz_0_yyyzz_1, g_zz_0_yyyzzz_1, g_zz_0_yyzzz_1, g_zz_0_yyzzzz_1, g_zz_0_yzzzz_1, g_zz_0_yzzzzz_1, g_zz_0_zzzzz_1, g_zz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xxxxxx_0[i] = 6.0 * g_zz_0_xxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxx_1[i] * wa_x[i];

        g_xzz_0_xxxxxy_0[i] = 5.0 * g_zz_0_xxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxy_1[i] * wa_x[i];

        g_xzz_0_xxxxxz_0[i] = 5.0 * g_zz_0_xxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxz_1[i] * wa_x[i];

        g_xzz_0_xxxxyy_0[i] = 4.0 * g_zz_0_xxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxyy_1[i] * wa_x[i];

        g_xzz_0_xxxxyz_0[i] = 4.0 * g_zz_0_xxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxyz_1[i] * wa_x[i];

        g_xzz_0_xxxxzz_0[i] = 4.0 * g_zz_0_xxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxzz_1[i] * wa_x[i];

        g_xzz_0_xxxyyy_0[i] = 3.0 * g_zz_0_xxyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyy_1[i] * wa_x[i];

        g_xzz_0_xxxyyz_0[i] = 3.0 * g_zz_0_xxyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyz_1[i] * wa_x[i];

        g_xzz_0_xxxyzz_0[i] = 3.0 * g_zz_0_xxyzz_1[i] * fi_acd_0 + g_zz_0_xxxyzz_1[i] * wa_x[i];

        g_xzz_0_xxxzzz_0[i] = 3.0 * g_zz_0_xxzzz_1[i] * fi_acd_0 + g_zz_0_xxxzzz_1[i] * wa_x[i];

        g_xzz_0_xxyyyy_0[i] = 2.0 * g_zz_0_xyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyy_1[i] * wa_x[i];

        g_xzz_0_xxyyyz_0[i] = 2.0 * g_zz_0_xyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyz_1[i] * wa_x[i];

        g_xzz_0_xxyyzz_0[i] = 2.0 * g_zz_0_xyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyzz_1[i] * wa_x[i];

        g_xzz_0_xxyzzz_0[i] = 2.0 * g_zz_0_xyzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzz_1[i] * wa_x[i];

        g_xzz_0_xxzzzz_0[i] = 2.0 * g_zz_0_xzzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyyy_0[i] = g_zz_0_yyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyy_1[i] * wa_x[i];

        g_xzz_0_xyyyyz_0[i] = g_zz_0_yyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyz_1[i] * wa_x[i];

        g_xzz_0_xyyyzz_0[i] = g_zz_0_yyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyzz_1[i] * wa_x[i];

        g_xzz_0_xyyzzz_0[i] = g_zz_0_yyzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzz_1[i] * wa_x[i];

        g_xzz_0_xyzzzz_0[i] = g_zz_0_yzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzz_1[i] * wa_x[i];

        g_xzz_0_xzzzzz_0[i] = g_zz_0_zzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyyy_0[i] = g_zz_0_yyyyyy_1[i] * wa_x[i];

        g_xzz_0_yyyyyz_0[i] = g_zz_0_yyyyyz_1[i] * wa_x[i];

        g_xzz_0_yyyyzz_0[i] = g_zz_0_yyyyzz_1[i] * wa_x[i];

        g_xzz_0_yyyzzz_0[i] = g_zz_0_yyyzzz_1[i] * wa_x[i];

        g_xzz_0_yyzzzz_0[i] = g_zz_0_yyzzzz_1[i] * wa_x[i];

        g_xzz_0_yzzzzz_0[i] = g_zz_0_yzzzzz_1[i] * wa_x[i];

        g_xzz_0_zzzzzz_0[i] = g_zz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 168-196 components of targeted buffer : FSI

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

    #pragma omp simd aligned(g_y_0_xxxxxx_0, g_y_0_xxxxxx_1, g_y_0_xxxxxy_0, g_y_0_xxxxxy_1, g_y_0_xxxxxz_0, g_y_0_xxxxxz_1, g_y_0_xxxxyy_0, g_y_0_xxxxyy_1, g_y_0_xxxxyz_0, g_y_0_xxxxyz_1, g_y_0_xxxxzz_0, g_y_0_xxxxzz_1, g_y_0_xxxyyy_0, g_y_0_xxxyyy_1, g_y_0_xxxyyz_0, g_y_0_xxxyyz_1, g_y_0_xxxyzz_0, g_y_0_xxxyzz_1, g_y_0_xxxzzz_0, g_y_0_xxxzzz_1, g_y_0_xxyyyy_0, g_y_0_xxyyyy_1, g_y_0_xxyyyz_0, g_y_0_xxyyyz_1, g_y_0_xxyyzz_0, g_y_0_xxyyzz_1, g_y_0_xxyzzz_0, g_y_0_xxyzzz_1, g_y_0_xxzzzz_0, g_y_0_xxzzzz_1, g_y_0_xyyyyy_0, g_y_0_xyyyyy_1, g_y_0_xyyyyz_0, g_y_0_xyyyyz_1, g_y_0_xyyyzz_0, g_y_0_xyyyzz_1, g_y_0_xyyzzz_0, g_y_0_xyyzzz_1, g_y_0_xyzzzz_0, g_y_0_xyzzzz_1, g_y_0_xzzzzz_0, g_y_0_xzzzzz_1, g_y_0_yyyyyy_0, g_y_0_yyyyyy_1, g_y_0_yyyyyz_0, g_y_0_yyyyyz_1, g_y_0_yyyyzz_0, g_y_0_yyyyzz_1, g_y_0_yyyzzz_0, g_y_0_yyyzzz_1, g_y_0_yyzzzz_0, g_y_0_yyzzzz_1, g_y_0_yzzzzz_0, g_y_0_yzzzzz_1, g_y_0_zzzzzz_0, g_y_0_zzzzzz_1, g_yy_0_xxxxx_1, g_yy_0_xxxxxx_1, g_yy_0_xxxxxy_1, g_yy_0_xxxxxz_1, g_yy_0_xxxxy_1, g_yy_0_xxxxyy_1, g_yy_0_xxxxyz_1, g_yy_0_xxxxz_1, g_yy_0_xxxxzz_1, g_yy_0_xxxyy_1, g_yy_0_xxxyyy_1, g_yy_0_xxxyyz_1, g_yy_0_xxxyz_1, g_yy_0_xxxyzz_1, g_yy_0_xxxzz_1, g_yy_0_xxxzzz_1, g_yy_0_xxyyy_1, g_yy_0_xxyyyy_1, g_yy_0_xxyyyz_1, g_yy_0_xxyyz_1, g_yy_0_xxyyzz_1, g_yy_0_xxyzz_1, g_yy_0_xxyzzz_1, g_yy_0_xxzzz_1, g_yy_0_xxzzzz_1, g_yy_0_xyyyy_1, g_yy_0_xyyyyy_1, g_yy_0_xyyyyz_1, g_yy_0_xyyyz_1, g_yy_0_xyyyzz_1, g_yy_0_xyyzz_1, g_yy_0_xyyzzz_1, g_yy_0_xyzzz_1, g_yy_0_xyzzzz_1, g_yy_0_xzzzz_1, g_yy_0_xzzzzz_1, g_yy_0_yyyyy_1, g_yy_0_yyyyyy_1, g_yy_0_yyyyyz_1, g_yy_0_yyyyz_1, g_yy_0_yyyyzz_1, g_yy_0_yyyzz_1, g_yy_0_yyyzzz_1, g_yy_0_yyzzz_1, g_yy_0_yyzzzz_1, g_yy_0_yzzzz_1, g_yy_0_yzzzzz_1, g_yy_0_zzzzz_1, g_yy_0_zzzzzz_1, g_yyy_0_xxxxxx_0, g_yyy_0_xxxxxy_0, g_yyy_0_xxxxxz_0, g_yyy_0_xxxxyy_0, g_yyy_0_xxxxyz_0, g_yyy_0_xxxxzz_0, g_yyy_0_xxxyyy_0, g_yyy_0_xxxyyz_0, g_yyy_0_xxxyzz_0, g_yyy_0_xxxzzz_0, g_yyy_0_xxyyyy_0, g_yyy_0_xxyyyz_0, g_yyy_0_xxyyzz_0, g_yyy_0_xxyzzz_0, g_yyy_0_xxzzzz_0, g_yyy_0_xyyyyy_0, g_yyy_0_xyyyyz_0, g_yyy_0_xyyyzz_0, g_yyy_0_xyyzzz_0, g_yyy_0_xyzzzz_0, g_yyy_0_xzzzzz_0, g_yyy_0_yyyyyy_0, g_yyy_0_yyyyyz_0, g_yyy_0_yyyyzz_0, g_yyy_0_yyyzzz_0, g_yyy_0_yyzzzz_0, g_yyy_0_yzzzzz_0, g_yyy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xxxxxx_0[i] = 2.0 * g_y_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxx_1[i] * fz_be_0 + g_yy_0_xxxxxx_1[i] * wa_y[i];

        g_yyy_0_xxxxxy_0[i] = 2.0 * g_y_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxy_1[i] * fz_be_0 + g_yy_0_xxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxy_1[i] * wa_y[i];

        g_yyy_0_xxxxxz_0[i] = 2.0 * g_y_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxz_1[i] * fz_be_0 + g_yy_0_xxxxxz_1[i] * wa_y[i];

        g_yyy_0_xxxxyy_0[i] = 2.0 * g_y_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxyy_1[i] * wa_y[i];

        g_yyy_0_xxxxyz_0[i] = 2.0 * g_y_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyz_1[i] * fz_be_0 + g_yy_0_xxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxyz_1[i] * wa_y[i];

        g_yyy_0_xxxxzz_0[i] = 2.0 * g_y_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxzz_1[i] * fz_be_0 + g_yy_0_xxxxzz_1[i] * wa_y[i];

        g_yyy_0_xxxyyy_0[i] = 2.0 * g_y_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxyy_1[i] * fi_acd_0 + g_yy_0_xxxyyy_1[i] * wa_y[i];

        g_yyy_0_xxxyyz_0[i] = 2.0 * g_y_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxyz_1[i] * fi_acd_0 + g_yy_0_xxxyyz_1[i] * wa_y[i];

        g_yyy_0_xxxyzz_0[i] = 2.0 * g_y_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyzz_1[i] * fz_be_0 + g_yy_0_xxxzz_1[i] * fi_acd_0 + g_yy_0_xxxyzz_1[i] * wa_y[i];

        g_yyy_0_xxxzzz_0[i] = 2.0 * g_y_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxzzz_1[i] * fz_be_0 + g_yy_0_xxxzzz_1[i] * wa_y[i];

        g_yyy_0_xxyyyy_0[i] = 2.0 * g_y_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yy_0_xxyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyy_1[i] * wa_y[i];

        g_yyy_0_xxyyyz_0[i] = 2.0 * g_y_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyz_1[i] * wa_y[i];

        g_yyy_0_xxyyzz_0[i] = 2.0 * g_y_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxyzz_1[i] * fi_acd_0 + g_yy_0_xxyyzz_1[i] * wa_y[i];

        g_yyy_0_xxyzzz_0[i] = 2.0 * g_y_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyzzz_1[i] * fz_be_0 + g_yy_0_xxzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzz_1[i] * wa_y[i];

        g_yyy_0_xxzzzz_0[i] = 2.0 * g_y_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxzzzz_1[i] * fz_be_0 + g_yy_0_xxzzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyyy_0[i] = 2.0 * g_y_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yy_0_xyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyy_1[i] * wa_y[i];

        g_yyy_0_xyyyyz_0[i] = 2.0 * g_y_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yy_0_xyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyz_1[i] * wa_y[i];

        g_yyy_0_xyyyzz_0[i] = 2.0 * g_y_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyzz_1[i] * wa_y[i];

        g_yyy_0_xyyzzz_0[i] = 2.0 * g_y_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xyzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzz_1[i] * wa_y[i];

        g_yyy_0_xyzzzz_0[i] = 2.0 * g_y_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyzzzz_1[i] * fz_be_0 + g_yy_0_xzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzz_1[i] * wa_y[i];

        g_yyy_0_xzzzzz_0[i] = 2.0 * g_y_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xzzzzz_1[i] * fz_be_0 + g_yy_0_xzzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyyy_0[i] = 2.0 * g_y_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yy_0_yyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyy_1[i] * wa_y[i];

        g_yyy_0_yyyyyz_0[i] = 2.0 * g_y_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yy_0_yyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyyz_1[i] * wa_y[i];

        g_yyy_0_yyyyzz_0[i] = 2.0 * g_y_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yy_0_yyyzz_1[i] * fi_acd_0 + g_yy_0_yyyyzz_1[i] * wa_y[i];

        g_yyy_0_yyyzzz_0[i] = 2.0 * g_y_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_yyzzz_1[i] * fi_acd_0 + g_yy_0_yyyzzz_1[i] * wa_y[i];

        g_yyy_0_yyzzzz_0[i] = 2.0 * g_y_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_yzzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzz_1[i] * wa_y[i];

        g_yyy_0_yzzzzz_0[i] = 2.0 * g_y_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzz_1[i] * wa_y[i];

        g_yyy_0_zzzzzz_0[i] = 2.0 * g_y_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_zzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 196-224 components of targeted buffer : FSI

    auto g_yyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 196);

    auto g_yyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 197);

    auto g_yyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 198);

    auto g_yyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 199);

    auto g_yyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 200);

    auto g_yyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 201);

    auto g_yyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 202);

    auto g_yyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 203);

    auto g_yyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 204);

    auto g_yyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 205);

    auto g_yyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 206);

    auto g_yyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 207);

    auto g_yyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 208);

    auto g_yyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 209);

    auto g_yyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 210);

    auto g_yyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 211);

    auto g_yyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 212);

    auto g_yyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 213);

    auto g_yyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 214);

    auto g_yyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 215);

    auto g_yyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 216);

    auto g_yyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 217);

    auto g_yyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 218);

    auto g_yyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 219);

    auto g_yyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 220);

    auto g_yyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 221);

    auto g_yyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 222);

    auto g_yyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 223);

    #pragma omp simd aligned(g_yy_0_xxxxx_1, g_yy_0_xxxxxx_1, g_yy_0_xxxxxy_1, g_yy_0_xxxxxz_1, g_yy_0_xxxxy_1, g_yy_0_xxxxyy_1, g_yy_0_xxxxyz_1, g_yy_0_xxxxz_1, g_yy_0_xxxxzz_1, g_yy_0_xxxyy_1, g_yy_0_xxxyyy_1, g_yy_0_xxxyyz_1, g_yy_0_xxxyz_1, g_yy_0_xxxyzz_1, g_yy_0_xxxzz_1, g_yy_0_xxxzzz_1, g_yy_0_xxyyy_1, g_yy_0_xxyyyy_1, g_yy_0_xxyyyz_1, g_yy_0_xxyyz_1, g_yy_0_xxyyzz_1, g_yy_0_xxyzz_1, g_yy_0_xxyzzz_1, g_yy_0_xxzzz_1, g_yy_0_xxzzzz_1, g_yy_0_xyyyy_1, g_yy_0_xyyyyy_1, g_yy_0_xyyyyz_1, g_yy_0_xyyyz_1, g_yy_0_xyyyzz_1, g_yy_0_xyyzz_1, g_yy_0_xyyzzz_1, g_yy_0_xyzzz_1, g_yy_0_xyzzzz_1, g_yy_0_xzzzz_1, g_yy_0_xzzzzz_1, g_yy_0_yyyyy_1, g_yy_0_yyyyyy_1, g_yy_0_yyyyyz_1, g_yy_0_yyyyz_1, g_yy_0_yyyyzz_1, g_yy_0_yyyzz_1, g_yy_0_yyyzzz_1, g_yy_0_yyzzz_1, g_yy_0_yyzzzz_1, g_yy_0_yzzzz_1, g_yy_0_yzzzzz_1, g_yy_0_zzzzz_1, g_yy_0_zzzzzz_1, g_yyz_0_xxxxxx_0, g_yyz_0_xxxxxy_0, g_yyz_0_xxxxxz_0, g_yyz_0_xxxxyy_0, g_yyz_0_xxxxyz_0, g_yyz_0_xxxxzz_0, g_yyz_0_xxxyyy_0, g_yyz_0_xxxyyz_0, g_yyz_0_xxxyzz_0, g_yyz_0_xxxzzz_0, g_yyz_0_xxyyyy_0, g_yyz_0_xxyyyz_0, g_yyz_0_xxyyzz_0, g_yyz_0_xxyzzz_0, g_yyz_0_xxzzzz_0, g_yyz_0_xyyyyy_0, g_yyz_0_xyyyyz_0, g_yyz_0_xyyyzz_0, g_yyz_0_xyyzzz_0, g_yyz_0_xyzzzz_0, g_yyz_0_xzzzzz_0, g_yyz_0_yyyyyy_0, g_yyz_0_yyyyyz_0, g_yyz_0_yyyyzz_0, g_yyz_0_yyyzzz_0, g_yyz_0_yyzzzz_0, g_yyz_0_yzzzzz_0, g_yyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xxxxxx_0[i] = g_yy_0_xxxxxx_1[i] * wa_z[i];

        g_yyz_0_xxxxxy_0[i] = g_yy_0_xxxxxy_1[i] * wa_z[i];

        g_yyz_0_xxxxxz_0[i] = g_yy_0_xxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxz_1[i] * wa_z[i];

        g_yyz_0_xxxxyy_0[i] = g_yy_0_xxxxyy_1[i] * wa_z[i];

        g_yyz_0_xxxxyz_0[i] = g_yy_0_xxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxyz_1[i] * wa_z[i];

        g_yyz_0_xxxxzz_0[i] = 2.0 * g_yy_0_xxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxzz_1[i] * wa_z[i];

        g_yyz_0_xxxyyy_0[i] = g_yy_0_xxxyyy_1[i] * wa_z[i];

        g_yyz_0_xxxyyz_0[i] = g_yy_0_xxxyy_1[i] * fi_acd_0 + g_yy_0_xxxyyz_1[i] * wa_z[i];

        g_yyz_0_xxxyzz_0[i] = 2.0 * g_yy_0_xxxyz_1[i] * fi_acd_0 + g_yy_0_xxxyzz_1[i] * wa_z[i];

        g_yyz_0_xxxzzz_0[i] = 3.0 * g_yy_0_xxxzz_1[i] * fi_acd_0 + g_yy_0_xxxzzz_1[i] * wa_z[i];

        g_yyz_0_xxyyyy_0[i] = g_yy_0_xxyyyy_1[i] * wa_z[i];

        g_yyz_0_xxyyyz_0[i] = g_yy_0_xxyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyz_1[i] * wa_z[i];

        g_yyz_0_xxyyzz_0[i] = 2.0 * g_yy_0_xxyyz_1[i] * fi_acd_0 + g_yy_0_xxyyzz_1[i] * wa_z[i];

        g_yyz_0_xxyzzz_0[i] = 3.0 * g_yy_0_xxyzz_1[i] * fi_acd_0 + g_yy_0_xxyzzz_1[i] * wa_z[i];

        g_yyz_0_xxzzzz_0[i] = 4.0 * g_yy_0_xxzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyyy_0[i] = g_yy_0_xyyyyy_1[i] * wa_z[i];

        g_yyz_0_xyyyyz_0[i] = g_yy_0_xyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyz_1[i] * wa_z[i];

        g_yyz_0_xyyyzz_0[i] = 2.0 * g_yy_0_xyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyzz_1[i] * wa_z[i];

        g_yyz_0_xyyzzz_0[i] = 3.0 * g_yy_0_xyyzz_1[i] * fi_acd_0 + g_yy_0_xyyzzz_1[i] * wa_z[i];

        g_yyz_0_xyzzzz_0[i] = 4.0 * g_yy_0_xyzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzz_1[i] * wa_z[i];

        g_yyz_0_xzzzzz_0[i] = 5.0 * g_yy_0_xzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyyy_0[i] = g_yy_0_yyyyyy_1[i] * wa_z[i];

        g_yyz_0_yyyyyz_0[i] = g_yy_0_yyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyz_1[i] * wa_z[i];

        g_yyz_0_yyyyzz_0[i] = 2.0 * g_yy_0_yyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyzz_1[i] * wa_z[i];

        g_yyz_0_yyyzzz_0[i] = 3.0 * g_yy_0_yyyzz_1[i] * fi_acd_0 + g_yy_0_yyyzzz_1[i] * wa_z[i];

        g_yyz_0_yyzzzz_0[i] = 4.0 * g_yy_0_yyzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzz_1[i] * wa_z[i];

        g_yyz_0_yzzzzz_0[i] = 5.0 * g_yy_0_yzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzz_1[i] * wa_z[i];

        g_yyz_0_zzzzzz_0[i] = 6.0 * g_yy_0_zzzzz_1[i] * fi_acd_0 + g_yy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 224-252 components of targeted buffer : FSI

    auto g_yzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_fsi + 224);

    auto g_yzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_fsi + 225);

    auto g_yzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_fsi + 226);

    auto g_yzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_fsi + 227);

    auto g_yzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_fsi + 228);

    auto g_yzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_fsi + 229);

    auto g_yzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_fsi + 230);

    auto g_yzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_fsi + 231);

    auto g_yzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_fsi + 232);

    auto g_yzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_fsi + 233);

    auto g_yzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_fsi + 234);

    auto g_yzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_fsi + 235);

    auto g_yzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_fsi + 236);

    auto g_yzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_fsi + 237);

    auto g_yzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_fsi + 238);

    auto g_yzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 239);

    auto g_yzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 240);

    auto g_yzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 241);

    auto g_yzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 242);

    auto g_yzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 243);

    auto g_yzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 244);

    auto g_yzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_fsi + 245);

    auto g_yzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_fsi + 246);

    auto g_yzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_fsi + 247);

    auto g_yzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_fsi + 248);

    auto g_yzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_fsi + 249);

    auto g_yzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 250);

    auto g_yzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_fsi + 251);

    #pragma omp simd aligned(g_yzz_0_xxxxxx_0, g_yzz_0_xxxxxy_0, g_yzz_0_xxxxxz_0, g_yzz_0_xxxxyy_0, g_yzz_0_xxxxyz_0, g_yzz_0_xxxxzz_0, g_yzz_0_xxxyyy_0, g_yzz_0_xxxyyz_0, g_yzz_0_xxxyzz_0, g_yzz_0_xxxzzz_0, g_yzz_0_xxyyyy_0, g_yzz_0_xxyyyz_0, g_yzz_0_xxyyzz_0, g_yzz_0_xxyzzz_0, g_yzz_0_xxzzzz_0, g_yzz_0_xyyyyy_0, g_yzz_0_xyyyyz_0, g_yzz_0_xyyyzz_0, g_yzz_0_xyyzzz_0, g_yzz_0_xyzzzz_0, g_yzz_0_xzzzzz_0, g_yzz_0_yyyyyy_0, g_yzz_0_yyyyyz_0, g_yzz_0_yyyyzz_0, g_yzz_0_yyyzzz_0, g_yzz_0_yyzzzz_0, g_yzz_0_yzzzzz_0, g_yzz_0_zzzzzz_0, g_zz_0_xxxxx_1, g_zz_0_xxxxxx_1, g_zz_0_xxxxxy_1, g_zz_0_xxxxxz_1, g_zz_0_xxxxy_1, g_zz_0_xxxxyy_1, g_zz_0_xxxxyz_1, g_zz_0_xxxxz_1, g_zz_0_xxxxzz_1, g_zz_0_xxxyy_1, g_zz_0_xxxyyy_1, g_zz_0_xxxyyz_1, g_zz_0_xxxyz_1, g_zz_0_xxxyzz_1, g_zz_0_xxxzz_1, g_zz_0_xxxzzz_1, g_zz_0_xxyyy_1, g_zz_0_xxyyyy_1, g_zz_0_xxyyyz_1, g_zz_0_xxyyz_1, g_zz_0_xxyyzz_1, g_zz_0_xxyzz_1, g_zz_0_xxyzzz_1, g_zz_0_xxzzz_1, g_zz_0_xxzzzz_1, g_zz_0_xyyyy_1, g_zz_0_xyyyyy_1, g_zz_0_xyyyyz_1, g_zz_0_xyyyz_1, g_zz_0_xyyyzz_1, g_zz_0_xyyzz_1, g_zz_0_xyyzzz_1, g_zz_0_xyzzz_1, g_zz_0_xyzzzz_1, g_zz_0_xzzzz_1, g_zz_0_xzzzzz_1, g_zz_0_yyyyy_1, g_zz_0_yyyyyy_1, g_zz_0_yyyyyz_1, g_zz_0_yyyyz_1, g_zz_0_yyyyzz_1, g_zz_0_yyyzz_1, g_zz_0_yyyzzz_1, g_zz_0_yyzzz_1, g_zz_0_yyzzzz_1, g_zz_0_yzzzz_1, g_zz_0_yzzzzz_1, g_zz_0_zzzzz_1, g_zz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xxxxxx_0[i] = g_zz_0_xxxxxx_1[i] * wa_y[i];

        g_yzz_0_xxxxxy_0[i] = g_zz_0_xxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxy_1[i] * wa_y[i];

        g_yzz_0_xxxxxz_0[i] = g_zz_0_xxxxxz_1[i] * wa_y[i];

        g_yzz_0_xxxxyy_0[i] = 2.0 * g_zz_0_xxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxyy_1[i] * wa_y[i];

        g_yzz_0_xxxxyz_0[i] = g_zz_0_xxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxyz_1[i] * wa_y[i];

        g_yzz_0_xxxxzz_0[i] = g_zz_0_xxxxzz_1[i] * wa_y[i];

        g_yzz_0_xxxyyy_0[i] = 3.0 * g_zz_0_xxxyy_1[i] * fi_acd_0 + g_zz_0_xxxyyy_1[i] * wa_y[i];

        g_yzz_0_xxxyyz_0[i] = 2.0 * g_zz_0_xxxyz_1[i] * fi_acd_0 + g_zz_0_xxxyyz_1[i] * wa_y[i];

        g_yzz_0_xxxyzz_0[i] = g_zz_0_xxxzz_1[i] * fi_acd_0 + g_zz_0_xxxyzz_1[i] * wa_y[i];

        g_yzz_0_xxxzzz_0[i] = g_zz_0_xxxzzz_1[i] * wa_y[i];

        g_yzz_0_xxyyyy_0[i] = 4.0 * g_zz_0_xxyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyy_1[i] * wa_y[i];

        g_yzz_0_xxyyyz_0[i] = 3.0 * g_zz_0_xxyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyz_1[i] * wa_y[i];

        g_yzz_0_xxyyzz_0[i] = 2.0 * g_zz_0_xxyzz_1[i] * fi_acd_0 + g_zz_0_xxyyzz_1[i] * wa_y[i];

        g_yzz_0_xxyzzz_0[i] = g_zz_0_xxzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzz_1[i] * wa_y[i];

        g_yzz_0_xxzzzz_0[i] = g_zz_0_xxzzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyyy_0[i] = 5.0 * g_zz_0_xyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyy_1[i] * wa_y[i];

        g_yzz_0_xyyyyz_0[i] = 4.0 * g_zz_0_xyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyz_1[i] * wa_y[i];

        g_yzz_0_xyyyzz_0[i] = 3.0 * g_zz_0_xyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyzz_1[i] * wa_y[i];

        g_yzz_0_xyyzzz_0[i] = 2.0 * g_zz_0_xyzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzz_1[i] * wa_y[i];

        g_yzz_0_xyzzzz_0[i] = g_zz_0_xzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzz_1[i] * wa_y[i];

        g_yzz_0_xzzzzz_0[i] = g_zz_0_xzzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyyy_0[i] = 6.0 * g_zz_0_yyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyy_1[i] * wa_y[i];

        g_yzz_0_yyyyyz_0[i] = 5.0 * g_zz_0_yyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyyz_1[i] * wa_y[i];

        g_yzz_0_yyyyzz_0[i] = 4.0 * g_zz_0_yyyzz_1[i] * fi_acd_0 + g_zz_0_yyyyzz_1[i] * wa_y[i];

        g_yzz_0_yyyzzz_0[i] = 3.0 * g_zz_0_yyzzz_1[i] * fi_acd_0 + g_zz_0_yyyzzz_1[i] * wa_y[i];

        g_yzz_0_yyzzzz_0[i] = 2.0 * g_zz_0_yzzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzz_1[i] * wa_y[i];

        g_yzz_0_yzzzzz_0[i] = g_zz_0_zzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzz_1[i] * wa_y[i];

        g_yzz_0_zzzzzz_0[i] = g_zz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 252-280 components of targeted buffer : FSI

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

    #pragma omp simd aligned(g_z_0_xxxxxx_0, g_z_0_xxxxxx_1, g_z_0_xxxxxy_0, g_z_0_xxxxxy_1, g_z_0_xxxxxz_0, g_z_0_xxxxxz_1, g_z_0_xxxxyy_0, g_z_0_xxxxyy_1, g_z_0_xxxxyz_0, g_z_0_xxxxyz_1, g_z_0_xxxxzz_0, g_z_0_xxxxzz_1, g_z_0_xxxyyy_0, g_z_0_xxxyyy_1, g_z_0_xxxyyz_0, g_z_0_xxxyyz_1, g_z_0_xxxyzz_0, g_z_0_xxxyzz_1, g_z_0_xxxzzz_0, g_z_0_xxxzzz_1, g_z_0_xxyyyy_0, g_z_0_xxyyyy_1, g_z_0_xxyyyz_0, g_z_0_xxyyyz_1, g_z_0_xxyyzz_0, g_z_0_xxyyzz_1, g_z_0_xxyzzz_0, g_z_0_xxyzzz_1, g_z_0_xxzzzz_0, g_z_0_xxzzzz_1, g_z_0_xyyyyy_0, g_z_0_xyyyyy_1, g_z_0_xyyyyz_0, g_z_0_xyyyyz_1, g_z_0_xyyyzz_0, g_z_0_xyyyzz_1, g_z_0_xyyzzz_0, g_z_0_xyyzzz_1, g_z_0_xyzzzz_0, g_z_0_xyzzzz_1, g_z_0_xzzzzz_0, g_z_0_xzzzzz_1, g_z_0_yyyyyy_0, g_z_0_yyyyyy_1, g_z_0_yyyyyz_0, g_z_0_yyyyyz_1, g_z_0_yyyyzz_0, g_z_0_yyyyzz_1, g_z_0_yyyzzz_0, g_z_0_yyyzzz_1, g_z_0_yyzzzz_0, g_z_0_yyzzzz_1, g_z_0_yzzzzz_0, g_z_0_yzzzzz_1, g_z_0_zzzzzz_0, g_z_0_zzzzzz_1, g_zz_0_xxxxx_1, g_zz_0_xxxxxx_1, g_zz_0_xxxxxy_1, g_zz_0_xxxxxz_1, g_zz_0_xxxxy_1, g_zz_0_xxxxyy_1, g_zz_0_xxxxyz_1, g_zz_0_xxxxz_1, g_zz_0_xxxxzz_1, g_zz_0_xxxyy_1, g_zz_0_xxxyyy_1, g_zz_0_xxxyyz_1, g_zz_0_xxxyz_1, g_zz_0_xxxyzz_1, g_zz_0_xxxzz_1, g_zz_0_xxxzzz_1, g_zz_0_xxyyy_1, g_zz_0_xxyyyy_1, g_zz_0_xxyyyz_1, g_zz_0_xxyyz_1, g_zz_0_xxyyzz_1, g_zz_0_xxyzz_1, g_zz_0_xxyzzz_1, g_zz_0_xxzzz_1, g_zz_0_xxzzzz_1, g_zz_0_xyyyy_1, g_zz_0_xyyyyy_1, g_zz_0_xyyyyz_1, g_zz_0_xyyyz_1, g_zz_0_xyyyzz_1, g_zz_0_xyyzz_1, g_zz_0_xyyzzz_1, g_zz_0_xyzzz_1, g_zz_0_xyzzzz_1, g_zz_0_xzzzz_1, g_zz_0_xzzzzz_1, g_zz_0_yyyyy_1, g_zz_0_yyyyyy_1, g_zz_0_yyyyyz_1, g_zz_0_yyyyz_1, g_zz_0_yyyyzz_1, g_zz_0_yyyzz_1, g_zz_0_yyyzzz_1, g_zz_0_yyzzz_1, g_zz_0_yyzzzz_1, g_zz_0_yzzzz_1, g_zz_0_yzzzzz_1, g_zz_0_zzzzz_1, g_zz_0_zzzzzz_1, g_zzz_0_xxxxxx_0, g_zzz_0_xxxxxy_0, g_zzz_0_xxxxxz_0, g_zzz_0_xxxxyy_0, g_zzz_0_xxxxyz_0, g_zzz_0_xxxxzz_0, g_zzz_0_xxxyyy_0, g_zzz_0_xxxyyz_0, g_zzz_0_xxxyzz_0, g_zzz_0_xxxzzz_0, g_zzz_0_xxyyyy_0, g_zzz_0_xxyyyz_0, g_zzz_0_xxyyzz_0, g_zzz_0_xxyzzz_0, g_zzz_0_xxzzzz_0, g_zzz_0_xyyyyy_0, g_zzz_0_xyyyyz_0, g_zzz_0_xyyyzz_0, g_zzz_0_xyyzzz_0, g_zzz_0_xyzzzz_0, g_zzz_0_xzzzzz_0, g_zzz_0_yyyyyy_0, g_zzz_0_yyyyyz_0, g_zzz_0_yyyyzz_0, g_zzz_0_yyyzzz_0, g_zzz_0_yyzzzz_0, g_zzz_0_yzzzzz_0, g_zzz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xxxxxx_0[i] = 2.0 * g_z_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxx_1[i] * fz_be_0 + g_zz_0_xxxxxx_1[i] * wa_z[i];

        g_zzz_0_xxxxxy_0[i] = 2.0 * g_z_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxy_1[i] * fz_be_0 + g_zz_0_xxxxxy_1[i] * wa_z[i];

        g_zzz_0_xxxxxz_0[i] = 2.0 * g_z_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxz_1[i] * fz_be_0 + g_zz_0_xxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxz_1[i] * wa_z[i];

        g_zzz_0_xxxxyy_0[i] = 2.0 * g_z_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyy_1[i] * fz_be_0 + g_zz_0_xxxxyy_1[i] * wa_z[i];

        g_zzz_0_xxxxyz_0[i] = 2.0 * g_z_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyz_1[i] * fz_be_0 + g_zz_0_xxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxyz_1[i] * wa_z[i];

        g_zzz_0_xxxxzz_0[i] = 2.0 * g_z_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxzz_1[i] * wa_z[i];

        g_zzz_0_xxxyyy_0[i] = 2.0 * g_z_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyy_1[i] * fz_be_0 + g_zz_0_xxxyyy_1[i] * wa_z[i];

        g_zzz_0_xxxyyz_0[i] = 2.0 * g_z_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyz_1[i] * fz_be_0 + g_zz_0_xxxyy_1[i] * fi_acd_0 + g_zz_0_xxxyyz_1[i] * wa_z[i];

        g_zzz_0_xxxyzz_0[i] = 2.0 * g_z_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxyz_1[i] * fi_acd_0 + g_zz_0_xxxyzz_1[i] * wa_z[i];

        g_zzz_0_xxxzzz_0[i] = 2.0 * g_z_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxzz_1[i] * fi_acd_0 + g_zz_0_xxxzzz_1[i] * wa_z[i];

        g_zzz_0_xxyyyy_0[i] = 2.0 * g_z_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyy_1[i] * fz_be_0 + g_zz_0_xxyyyy_1[i] * wa_z[i];

        g_zzz_0_xxyyyz_0[i] = 2.0 * g_z_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyz_1[i] * fz_be_0 + g_zz_0_xxyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyz_1[i] * wa_z[i];

        g_zzz_0_xxyyzz_0[i] = 2.0 * g_z_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxyyz_1[i] * fi_acd_0 + g_zz_0_xxyyzz_1[i] * wa_z[i];

        g_zzz_0_xxyzzz_0[i] = 2.0 * g_z_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxyzz_1[i] * fi_acd_0 + g_zz_0_xxyzzz_1[i] * wa_z[i];

        g_zzz_0_xxzzzz_0[i] = 2.0 * g_z_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyyy_0[i] = 2.0 * g_z_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyy_1[i] * fz_be_0 + g_zz_0_xyyyyy_1[i] * wa_z[i];

        g_zzz_0_xyyyyz_0[i] = 2.0 * g_z_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyz_1[i] * fz_be_0 + g_zz_0_xyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyz_1[i] * wa_z[i];

        g_zzz_0_xyyyzz_0[i] = 2.0 * g_z_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyzz_1[i] * wa_z[i];

        g_zzz_0_xyyzzz_0[i] = 2.0 * g_z_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xyyzz_1[i] * fi_acd_0 + g_zz_0_xyyzzz_1[i] * wa_z[i];

        g_zzz_0_xyzzzz_0[i] = 2.0 * g_z_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xyzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzz_1[i] * wa_z[i];

        g_zzz_0_xzzzzz_0[i] = 2.0 * g_z_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyyy_0[i] = 2.0 * g_z_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyy_1[i] * fz_be_0 + g_zz_0_yyyyyy_1[i] * wa_z[i];

        g_zzz_0_yyyyyz_0[i] = 2.0 * g_z_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyz_1[i] * fz_be_0 + g_zz_0_yyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyz_1[i] * wa_z[i];

        g_zzz_0_yyyyzz_0[i] = 2.0 * g_z_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_yyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyzz_1[i] * wa_z[i];

        g_zzz_0_yyyzzz_0[i] = 2.0 * g_z_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_yyyzz_1[i] * fi_acd_0 + g_zz_0_yyyzzz_1[i] * wa_z[i];

        g_zzz_0_yyzzzz_0[i] = 2.0 * g_z_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_yyzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzz_1[i] * wa_z[i];

        g_zzz_0_yzzzzz_0[i] = 2.0 * g_z_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_yzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzz_1[i] * wa_z[i];

        g_zzz_0_zzzzzz_0[i] = 2.0 * g_z_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_zzzzz_1[i] * fi_acd_0 + g_zz_0_zzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

