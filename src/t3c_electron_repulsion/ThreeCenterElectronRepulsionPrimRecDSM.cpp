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

#include "ThreeCenterElectronRepulsionPrimRecDSM.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsm,
                                 size_t idx_eri_0_ssm,
                                 size_t idx_eri_1_ssm,
                                 size_t idx_eri_1_psl,
                                 size_t idx_eri_1_psm,
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

    /// Set up components of auxilary buffer : SSM

    auto g_0_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_ssm);

    auto g_0_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_ssm + 1);

    auto g_0_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_ssm + 2);

    auto g_0_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_ssm + 3);

    auto g_0_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_ssm + 4);

    auto g_0_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_ssm + 5);

    auto g_0_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_ssm + 6);

    auto g_0_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_ssm + 7);

    auto g_0_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_ssm + 8);

    auto g_0_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_ssm + 9);

    auto g_0_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_ssm + 10);

    auto g_0_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_ssm + 11);

    auto g_0_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_ssm + 12);

    auto g_0_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_ssm + 13);

    auto g_0_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_ssm + 14);

    auto g_0_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 15);

    auto g_0_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 16);

    auto g_0_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 17);

    auto g_0_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 18);

    auto g_0_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 19);

    auto g_0_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 20);

    auto g_0_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 21);

    auto g_0_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 22);

    auto g_0_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 23);

    auto g_0_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 24);

    auto g_0_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 25);

    auto g_0_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 26);

    auto g_0_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 27);

    auto g_0_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 28);

    auto g_0_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 29);

    auto g_0_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 30);

    auto g_0_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 31);

    auto g_0_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 32);

    auto g_0_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 33);

    auto g_0_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 34);

    auto g_0_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 35);

    auto g_0_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 36);

    auto g_0_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 37);

    auto g_0_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 38);

    auto g_0_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 39);

    auto g_0_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 40);

    auto g_0_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 41);

    auto g_0_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 42);

    auto g_0_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 43);

    auto g_0_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 44);

    auto g_0_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 45);

    auto g_0_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 46);

    auto g_0_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 47);

    auto g_0_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 48);

    auto g_0_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 49);

    auto g_0_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 50);

    auto g_0_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 51);

    auto g_0_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 52);

    auto g_0_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 53);

    auto g_0_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 54);

    /// Set up components of auxilary buffer : SSM

    auto g_0_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_ssm);

    auto g_0_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_ssm + 1);

    auto g_0_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_ssm + 2);

    auto g_0_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_ssm + 3);

    auto g_0_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_ssm + 4);

    auto g_0_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_ssm + 5);

    auto g_0_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_ssm + 6);

    auto g_0_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_ssm + 7);

    auto g_0_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_ssm + 8);

    auto g_0_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_ssm + 9);

    auto g_0_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_ssm + 10);

    auto g_0_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_ssm + 11);

    auto g_0_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_ssm + 12);

    auto g_0_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_ssm + 13);

    auto g_0_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_ssm + 14);

    auto g_0_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 15);

    auto g_0_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 16);

    auto g_0_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 17);

    auto g_0_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 18);

    auto g_0_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 19);

    auto g_0_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 20);

    auto g_0_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 21);

    auto g_0_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 22);

    auto g_0_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 23);

    auto g_0_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 24);

    auto g_0_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 25);

    auto g_0_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 26);

    auto g_0_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 27);

    auto g_0_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 28);

    auto g_0_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 29);

    auto g_0_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 30);

    auto g_0_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 31);

    auto g_0_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 32);

    auto g_0_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 33);

    auto g_0_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 34);

    auto g_0_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 35);

    auto g_0_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 36);

    auto g_0_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 37);

    auto g_0_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 38);

    auto g_0_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 39);

    auto g_0_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 40);

    auto g_0_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 41);

    auto g_0_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 42);

    auto g_0_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 43);

    auto g_0_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 44);

    auto g_0_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 45);

    auto g_0_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 46);

    auto g_0_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 47);

    auto g_0_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 48);

    auto g_0_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 49);

    auto g_0_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 50);

    auto g_0_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 51);

    auto g_0_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 52);

    auto g_0_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 53);

    auto g_0_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 54);

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

    /// Set up components of auxilary buffer : PSM

    auto g_x_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_psm);

    auto g_x_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_psm + 1);

    auto g_x_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_psm + 2);

    auto g_x_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_psm + 3);

    auto g_x_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_psm + 4);

    auto g_x_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_psm + 5);

    auto g_x_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_psm + 6);

    auto g_x_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_psm + 7);

    auto g_x_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_psm + 8);

    auto g_x_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_psm + 9);

    auto g_x_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_psm + 10);

    auto g_x_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_psm + 11);

    auto g_x_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_psm + 12);

    auto g_x_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_psm + 13);

    auto g_x_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_psm + 14);

    auto g_x_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_psm + 15);

    auto g_x_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_psm + 16);

    auto g_x_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_psm + 17);

    auto g_x_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_psm + 18);

    auto g_x_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_psm + 19);

    auto g_x_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_psm + 20);

    auto g_x_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 21);

    auto g_x_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 22);

    auto g_x_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 23);

    auto g_x_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 24);

    auto g_x_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 25);

    auto g_x_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 26);

    auto g_x_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 27);

    auto g_x_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 28);

    auto g_x_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 29);

    auto g_x_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 30);

    auto g_x_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 31);

    auto g_x_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 32);

    auto g_x_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 33);

    auto g_x_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 34);

    auto g_x_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 35);

    auto g_x_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 36);

    auto g_x_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 37);

    auto g_x_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 38);

    auto g_x_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 39);

    auto g_x_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 40);

    auto g_x_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 41);

    auto g_x_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 42);

    auto g_x_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 43);

    auto g_x_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 44);

    auto g_x_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 45);

    auto g_x_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 46);

    auto g_x_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 47);

    auto g_x_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 48);

    auto g_x_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 49);

    auto g_x_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 50);

    auto g_x_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 51);

    auto g_x_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 52);

    auto g_x_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 53);

    auto g_x_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 54);

    auto g_y_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_psm + 55);

    auto g_y_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_psm + 56);

    auto g_y_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_psm + 57);

    auto g_y_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_psm + 58);

    auto g_y_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_psm + 59);

    auto g_y_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_psm + 60);

    auto g_y_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_psm + 61);

    auto g_y_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_psm + 62);

    auto g_y_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_psm + 63);

    auto g_y_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_psm + 64);

    auto g_y_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_psm + 65);

    auto g_y_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_psm + 66);

    auto g_y_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_psm + 67);

    auto g_y_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_psm + 68);

    auto g_y_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_psm + 69);

    auto g_y_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_psm + 70);

    auto g_y_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_psm + 71);

    auto g_y_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_psm + 72);

    auto g_y_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_psm + 73);

    auto g_y_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_psm + 74);

    auto g_y_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_psm + 75);

    auto g_y_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 76);

    auto g_y_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 77);

    auto g_y_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 78);

    auto g_y_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 79);

    auto g_y_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 80);

    auto g_y_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 81);

    auto g_y_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 82);

    auto g_y_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 83);

    auto g_y_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 84);

    auto g_y_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 85);

    auto g_y_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 86);

    auto g_y_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 87);

    auto g_y_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 88);

    auto g_y_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 89);

    auto g_y_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 90);

    auto g_y_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 91);

    auto g_y_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 92);

    auto g_y_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 93);

    auto g_y_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 94);

    auto g_y_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 95);

    auto g_y_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 96);

    auto g_y_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 97);

    auto g_y_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 98);

    auto g_y_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 99);

    auto g_y_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 100);

    auto g_y_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 101);

    auto g_y_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 102);

    auto g_y_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 103);

    auto g_y_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 104);

    auto g_y_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 105);

    auto g_y_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 106);

    auto g_y_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 107);

    auto g_y_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 108);

    auto g_y_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 109);

    auto g_z_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_psm + 110);

    auto g_z_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_psm + 111);

    auto g_z_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_psm + 112);

    auto g_z_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_psm + 113);

    auto g_z_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_psm + 114);

    auto g_z_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_psm + 115);

    auto g_z_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_psm + 116);

    auto g_z_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_psm + 117);

    auto g_z_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_psm + 118);

    auto g_z_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_psm + 119);

    auto g_z_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_psm + 120);

    auto g_z_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_psm + 121);

    auto g_z_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_psm + 122);

    auto g_z_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_psm + 123);

    auto g_z_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_psm + 124);

    auto g_z_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_psm + 125);

    auto g_z_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_psm + 126);

    auto g_z_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_psm + 127);

    auto g_z_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_psm + 128);

    auto g_z_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_psm + 129);

    auto g_z_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_psm + 130);

    auto g_z_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 131);

    auto g_z_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 132);

    auto g_z_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 133);

    auto g_z_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 134);

    auto g_z_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 135);

    auto g_z_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 136);

    auto g_z_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 137);

    auto g_z_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 138);

    auto g_z_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 139);

    auto g_z_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 140);

    auto g_z_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 141);

    auto g_z_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 142);

    auto g_z_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 143);

    auto g_z_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 144);

    auto g_z_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 145);

    auto g_z_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 146);

    auto g_z_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 147);

    auto g_z_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 148);

    auto g_z_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 149);

    auto g_z_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 150);

    auto g_z_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 151);

    auto g_z_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 152);

    auto g_z_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 153);

    auto g_z_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 154);

    auto g_z_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_psm + 155);

    auto g_z_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_psm + 156);

    auto g_z_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_psm + 157);

    auto g_z_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_psm + 158);

    auto g_z_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_psm + 159);

    auto g_z_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_psm + 160);

    auto g_z_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 161);

    auto g_z_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 162);

    auto g_z_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 163);

    auto g_z_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_psm + 164);

    /// Set up 0-55 components of targeted buffer : DSM

    auto g_xx_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm);

    auto g_xx_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 1);

    auto g_xx_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 2);

    auto g_xx_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 3);

    auto g_xx_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 4);

    auto g_xx_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 5);

    auto g_xx_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 6);

    auto g_xx_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 7);

    auto g_xx_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 8);

    auto g_xx_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 9);

    auto g_xx_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 10);

    auto g_xx_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 11);

    auto g_xx_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 12);

    auto g_xx_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 13);

    auto g_xx_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 14);

    auto g_xx_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 15);

    auto g_xx_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 16);

    auto g_xx_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 17);

    auto g_xx_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 18);

    auto g_xx_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 19);

    auto g_xx_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 20);

    auto g_xx_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 21);

    auto g_xx_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 22);

    auto g_xx_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 23);

    auto g_xx_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 24);

    auto g_xx_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 25);

    auto g_xx_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 26);

    auto g_xx_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 27);

    auto g_xx_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 28);

    auto g_xx_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 29);

    auto g_xx_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 30);

    auto g_xx_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 31);

    auto g_xx_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 32);

    auto g_xx_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 33);

    auto g_xx_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 34);

    auto g_xx_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 35);

    auto g_xx_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 36);

    auto g_xx_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 37);

    auto g_xx_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 38);

    auto g_xx_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 39);

    auto g_xx_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 40);

    auto g_xx_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 41);

    auto g_xx_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 42);

    auto g_xx_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 43);

    auto g_xx_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 44);

    auto g_xx_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 45);

    auto g_xx_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 46);

    auto g_xx_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 47);

    auto g_xx_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 48);

    auto g_xx_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 49);

    auto g_xx_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 50);

    auto g_xx_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 51);

    auto g_xx_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 52);

    auto g_xx_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 53);

    auto g_xx_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 54);

    #pragma omp simd aligned(g_0_0_xxxxxxxxx_0, g_0_0_xxxxxxxxx_1, g_0_0_xxxxxxxxy_0, g_0_0_xxxxxxxxy_1, g_0_0_xxxxxxxxz_0, g_0_0_xxxxxxxxz_1, g_0_0_xxxxxxxyy_0, g_0_0_xxxxxxxyy_1, g_0_0_xxxxxxxyz_0, g_0_0_xxxxxxxyz_1, g_0_0_xxxxxxxzz_0, g_0_0_xxxxxxxzz_1, g_0_0_xxxxxxyyy_0, g_0_0_xxxxxxyyy_1, g_0_0_xxxxxxyyz_0, g_0_0_xxxxxxyyz_1, g_0_0_xxxxxxyzz_0, g_0_0_xxxxxxyzz_1, g_0_0_xxxxxxzzz_0, g_0_0_xxxxxxzzz_1, g_0_0_xxxxxyyyy_0, g_0_0_xxxxxyyyy_1, g_0_0_xxxxxyyyz_0, g_0_0_xxxxxyyyz_1, g_0_0_xxxxxyyzz_0, g_0_0_xxxxxyyzz_1, g_0_0_xxxxxyzzz_0, g_0_0_xxxxxyzzz_1, g_0_0_xxxxxzzzz_0, g_0_0_xxxxxzzzz_1, g_0_0_xxxxyyyyy_0, g_0_0_xxxxyyyyy_1, g_0_0_xxxxyyyyz_0, g_0_0_xxxxyyyyz_1, g_0_0_xxxxyyyzz_0, g_0_0_xxxxyyyzz_1, g_0_0_xxxxyyzzz_0, g_0_0_xxxxyyzzz_1, g_0_0_xxxxyzzzz_0, g_0_0_xxxxyzzzz_1, g_0_0_xxxxzzzzz_0, g_0_0_xxxxzzzzz_1, g_0_0_xxxyyyyyy_0, g_0_0_xxxyyyyyy_1, g_0_0_xxxyyyyyz_0, g_0_0_xxxyyyyyz_1, g_0_0_xxxyyyyzz_0, g_0_0_xxxyyyyzz_1, g_0_0_xxxyyyzzz_0, g_0_0_xxxyyyzzz_1, g_0_0_xxxyyzzzz_0, g_0_0_xxxyyzzzz_1, g_0_0_xxxyzzzzz_0, g_0_0_xxxyzzzzz_1, g_0_0_xxxzzzzzz_0, g_0_0_xxxzzzzzz_1, g_0_0_xxyyyyyyy_0, g_0_0_xxyyyyyyy_1, g_0_0_xxyyyyyyz_0, g_0_0_xxyyyyyyz_1, g_0_0_xxyyyyyzz_0, g_0_0_xxyyyyyzz_1, g_0_0_xxyyyyzzz_0, g_0_0_xxyyyyzzz_1, g_0_0_xxyyyzzzz_0, g_0_0_xxyyyzzzz_1, g_0_0_xxyyzzzzz_0, g_0_0_xxyyzzzzz_1, g_0_0_xxyzzzzzz_0, g_0_0_xxyzzzzzz_1, g_0_0_xxzzzzzzz_0, g_0_0_xxzzzzzzz_1, g_0_0_xyyyyyyyy_0, g_0_0_xyyyyyyyy_1, g_0_0_xyyyyyyyz_0, g_0_0_xyyyyyyyz_1, g_0_0_xyyyyyyzz_0, g_0_0_xyyyyyyzz_1, g_0_0_xyyyyyzzz_0, g_0_0_xyyyyyzzz_1, g_0_0_xyyyyzzzz_0, g_0_0_xyyyyzzzz_1, g_0_0_xyyyzzzzz_0, g_0_0_xyyyzzzzz_1, g_0_0_xyyzzzzzz_0, g_0_0_xyyzzzzzz_1, g_0_0_xyzzzzzzz_0, g_0_0_xyzzzzzzz_1, g_0_0_xzzzzzzzz_0, g_0_0_xzzzzzzzz_1, g_0_0_yyyyyyyyy_0, g_0_0_yyyyyyyyy_1, g_0_0_yyyyyyyyz_0, g_0_0_yyyyyyyyz_1, g_0_0_yyyyyyyzz_0, g_0_0_yyyyyyyzz_1, g_0_0_yyyyyyzzz_0, g_0_0_yyyyyyzzz_1, g_0_0_yyyyyzzzz_0, g_0_0_yyyyyzzzz_1, g_0_0_yyyyzzzzz_0, g_0_0_yyyyzzzzz_1, g_0_0_yyyzzzzzz_0, g_0_0_yyyzzzzzz_1, g_0_0_yyzzzzzzz_0, g_0_0_yyzzzzzzz_1, g_0_0_yzzzzzzzz_0, g_0_0_yzzzzzzzz_1, g_0_0_zzzzzzzzz_0, g_0_0_zzzzzzzzz_1, g_x_0_xxxxxxxx_1, g_x_0_xxxxxxxxx_1, g_x_0_xxxxxxxxy_1, g_x_0_xxxxxxxxz_1, g_x_0_xxxxxxxy_1, g_x_0_xxxxxxxyy_1, g_x_0_xxxxxxxyz_1, g_x_0_xxxxxxxz_1, g_x_0_xxxxxxxzz_1, g_x_0_xxxxxxyy_1, g_x_0_xxxxxxyyy_1, g_x_0_xxxxxxyyz_1, g_x_0_xxxxxxyz_1, g_x_0_xxxxxxyzz_1, g_x_0_xxxxxxzz_1, g_x_0_xxxxxxzzz_1, g_x_0_xxxxxyyy_1, g_x_0_xxxxxyyyy_1, g_x_0_xxxxxyyyz_1, g_x_0_xxxxxyyz_1, g_x_0_xxxxxyyzz_1, g_x_0_xxxxxyzz_1, g_x_0_xxxxxyzzz_1, g_x_0_xxxxxzzz_1, g_x_0_xxxxxzzzz_1, g_x_0_xxxxyyyy_1, g_x_0_xxxxyyyyy_1, g_x_0_xxxxyyyyz_1, g_x_0_xxxxyyyz_1, g_x_0_xxxxyyyzz_1, g_x_0_xxxxyyzz_1, g_x_0_xxxxyyzzz_1, g_x_0_xxxxyzzz_1, g_x_0_xxxxyzzzz_1, g_x_0_xxxxzzzz_1, g_x_0_xxxxzzzzz_1, g_x_0_xxxyyyyy_1, g_x_0_xxxyyyyyy_1, g_x_0_xxxyyyyyz_1, g_x_0_xxxyyyyz_1, g_x_0_xxxyyyyzz_1, g_x_0_xxxyyyzz_1, g_x_0_xxxyyyzzz_1, g_x_0_xxxyyzzz_1, g_x_0_xxxyyzzzz_1, g_x_0_xxxyzzzz_1, g_x_0_xxxyzzzzz_1, g_x_0_xxxzzzzz_1, g_x_0_xxxzzzzzz_1, g_x_0_xxyyyyyy_1, g_x_0_xxyyyyyyy_1, g_x_0_xxyyyyyyz_1, g_x_0_xxyyyyyz_1, g_x_0_xxyyyyyzz_1, g_x_0_xxyyyyzz_1, g_x_0_xxyyyyzzz_1, g_x_0_xxyyyzzz_1, g_x_0_xxyyyzzzz_1, g_x_0_xxyyzzzz_1, g_x_0_xxyyzzzzz_1, g_x_0_xxyzzzzz_1, g_x_0_xxyzzzzzz_1, g_x_0_xxzzzzzz_1, g_x_0_xxzzzzzzz_1, g_x_0_xyyyyyyy_1, g_x_0_xyyyyyyyy_1, g_x_0_xyyyyyyyz_1, g_x_0_xyyyyyyz_1, g_x_0_xyyyyyyzz_1, g_x_0_xyyyyyzz_1, g_x_0_xyyyyyzzz_1, g_x_0_xyyyyzzz_1, g_x_0_xyyyyzzzz_1, g_x_0_xyyyzzzz_1, g_x_0_xyyyzzzzz_1, g_x_0_xyyzzzzz_1, g_x_0_xyyzzzzzz_1, g_x_0_xyzzzzzz_1, g_x_0_xyzzzzzzz_1, g_x_0_xzzzzzzz_1, g_x_0_xzzzzzzzz_1, g_x_0_yyyyyyyy_1, g_x_0_yyyyyyyyy_1, g_x_0_yyyyyyyyz_1, g_x_0_yyyyyyyz_1, g_x_0_yyyyyyyzz_1, g_x_0_yyyyyyzz_1, g_x_0_yyyyyyzzz_1, g_x_0_yyyyyzzz_1, g_x_0_yyyyyzzzz_1, g_x_0_yyyyzzzz_1, g_x_0_yyyyzzzzz_1, g_x_0_yyyzzzzz_1, g_x_0_yyyzzzzzz_1, g_x_0_yyzzzzzz_1, g_x_0_yyzzzzzzz_1, g_x_0_yzzzzzzz_1, g_x_0_yzzzzzzzz_1, g_x_0_zzzzzzzz_1, g_x_0_zzzzzzzzz_1, g_xx_0_xxxxxxxxx_0, g_xx_0_xxxxxxxxy_0, g_xx_0_xxxxxxxxz_0, g_xx_0_xxxxxxxyy_0, g_xx_0_xxxxxxxyz_0, g_xx_0_xxxxxxxzz_0, g_xx_0_xxxxxxyyy_0, g_xx_0_xxxxxxyyz_0, g_xx_0_xxxxxxyzz_0, g_xx_0_xxxxxxzzz_0, g_xx_0_xxxxxyyyy_0, g_xx_0_xxxxxyyyz_0, g_xx_0_xxxxxyyzz_0, g_xx_0_xxxxxyzzz_0, g_xx_0_xxxxxzzzz_0, g_xx_0_xxxxyyyyy_0, g_xx_0_xxxxyyyyz_0, g_xx_0_xxxxyyyzz_0, g_xx_0_xxxxyyzzz_0, g_xx_0_xxxxyzzzz_0, g_xx_0_xxxxzzzzz_0, g_xx_0_xxxyyyyyy_0, g_xx_0_xxxyyyyyz_0, g_xx_0_xxxyyyyzz_0, g_xx_0_xxxyyyzzz_0, g_xx_0_xxxyyzzzz_0, g_xx_0_xxxyzzzzz_0, g_xx_0_xxxzzzzzz_0, g_xx_0_xxyyyyyyy_0, g_xx_0_xxyyyyyyz_0, g_xx_0_xxyyyyyzz_0, g_xx_0_xxyyyyzzz_0, g_xx_0_xxyyyzzzz_0, g_xx_0_xxyyzzzzz_0, g_xx_0_xxyzzzzzz_0, g_xx_0_xxzzzzzzz_0, g_xx_0_xyyyyyyyy_0, g_xx_0_xyyyyyyyz_0, g_xx_0_xyyyyyyzz_0, g_xx_0_xyyyyyzzz_0, g_xx_0_xyyyyzzzz_0, g_xx_0_xyyyzzzzz_0, g_xx_0_xyyzzzzzz_0, g_xx_0_xyzzzzzzz_0, g_xx_0_xzzzzzzzz_0, g_xx_0_yyyyyyyyy_0, g_xx_0_yyyyyyyyz_0, g_xx_0_yyyyyyyzz_0, g_xx_0_yyyyyyzzz_0, g_xx_0_yyyyyzzzz_0, g_xx_0_yyyyzzzzz_0, g_xx_0_yyyzzzzzz_0, g_xx_0_yyzzzzzzz_0, g_xx_0_yzzzzzzzz_0, g_xx_0_zzzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xxxxxxxxx_0[i] = g_0_0_xxxxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxxxx_1[i] * fz_be_0 + 9.0 * g_x_0_xxxxxxxx_1[i] * fi_acd_0 + g_x_0_xxxxxxxxx_1[i] * wa_x[i];

        g_xx_0_xxxxxxxxy_0[i] = g_0_0_xxxxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxxxy_1[i] * fz_be_0 + 8.0 * g_x_0_xxxxxxxy_1[i] * fi_acd_0 + g_x_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xx_0_xxxxxxxxz_0[i] = g_0_0_xxxxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxxxz_1[i] * fz_be_0 + 8.0 * g_x_0_xxxxxxxz_1[i] * fi_acd_0 + g_x_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xx_0_xxxxxxxyy_0[i] = g_0_0_xxxxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxxxyy_1[i] * fz_be_0 + 7.0 * g_x_0_xxxxxxyy_1[i] * fi_acd_0 + g_x_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xx_0_xxxxxxxyz_0[i] = g_0_0_xxxxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxxxyz_1[i] * fz_be_0 + 7.0 * g_x_0_xxxxxxyz_1[i] * fi_acd_0 + g_x_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xx_0_xxxxxxxzz_0[i] = g_0_0_xxxxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxxxzz_1[i] * fz_be_0 + 7.0 * g_x_0_xxxxxxzz_1[i] * fi_acd_0 + g_x_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xx_0_xxxxxxyyy_0[i] = g_0_0_xxxxxxyyy_0[i] * fbe_0 - g_0_0_xxxxxxyyy_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxyyy_1[i] * fi_acd_0 + g_x_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xx_0_xxxxxxyyz_0[i] = g_0_0_xxxxxxyyz_0[i] * fbe_0 - g_0_0_xxxxxxyyz_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxyyz_1[i] * fi_acd_0 + g_x_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xx_0_xxxxxxyzz_0[i] = g_0_0_xxxxxxyzz_0[i] * fbe_0 - g_0_0_xxxxxxyzz_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxyzz_1[i] * fi_acd_0 + g_x_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xx_0_xxxxxxzzz_0[i] = g_0_0_xxxxxxzzz_0[i] * fbe_0 - g_0_0_xxxxxxzzz_1[i] * fz_be_0 + 6.0 * g_x_0_xxxxxzzz_1[i] * fi_acd_0 + g_x_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xx_0_xxxxxyyyy_0[i] = g_0_0_xxxxxyyyy_0[i] * fbe_0 - g_0_0_xxxxxyyyy_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyyyy_1[i] * fi_acd_0 + g_x_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xx_0_xxxxxyyyz_0[i] = g_0_0_xxxxxyyyz_0[i] * fbe_0 - g_0_0_xxxxxyyyz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyyyz_1[i] * fi_acd_0 + g_x_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xx_0_xxxxxyyzz_0[i] = g_0_0_xxxxxyyzz_0[i] * fbe_0 - g_0_0_xxxxxyyzz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyyzz_1[i] * fi_acd_0 + g_x_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xx_0_xxxxxyzzz_0[i] = g_0_0_xxxxxyzzz_0[i] * fbe_0 - g_0_0_xxxxxyzzz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxyzzz_1[i] * fi_acd_0 + g_x_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xx_0_xxxxxzzzz_0[i] = g_0_0_xxxxxzzzz_0[i] * fbe_0 - g_0_0_xxxxxzzzz_1[i] * fz_be_0 + 5.0 * g_x_0_xxxxzzzz_1[i] * fi_acd_0 + g_x_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xx_0_xxxxyyyyy_0[i] = g_0_0_xxxxyyyyy_0[i] * fbe_0 - g_0_0_xxxxyyyyy_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyyyy_1[i] * fi_acd_0 + g_x_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xx_0_xxxxyyyyz_0[i] = g_0_0_xxxxyyyyz_0[i] * fbe_0 - g_0_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyyyz_1[i] * fi_acd_0 + g_x_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xx_0_xxxxyyyzz_0[i] = g_0_0_xxxxyyyzz_0[i] * fbe_0 - g_0_0_xxxxyyyzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyyzz_1[i] * fi_acd_0 + g_x_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xx_0_xxxxyyzzz_0[i] = g_0_0_xxxxyyzzz_0[i] * fbe_0 - g_0_0_xxxxyyzzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyyzzz_1[i] * fi_acd_0 + g_x_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xx_0_xxxxyzzzz_0[i] = g_0_0_xxxxyzzzz_0[i] * fbe_0 - g_0_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxyzzzz_1[i] * fi_acd_0 + g_x_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xx_0_xxxxzzzzz_0[i] = g_0_0_xxxxzzzzz_0[i] * fbe_0 - g_0_0_xxxxzzzzz_1[i] * fz_be_0 + 4.0 * g_x_0_xxxzzzzz_1[i] * fi_acd_0 + g_x_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xx_0_xxxyyyyyy_0[i] = g_0_0_xxxyyyyyy_0[i] * fbe_0 - g_0_0_xxxyyyyyy_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyyyy_1[i] * fi_acd_0 + g_x_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xx_0_xxxyyyyyz_0[i] = g_0_0_xxxyyyyyz_0[i] * fbe_0 - g_0_0_xxxyyyyyz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyyyz_1[i] * fi_acd_0 + g_x_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xx_0_xxxyyyyzz_0[i] = g_0_0_xxxyyyyzz_0[i] * fbe_0 - g_0_0_xxxyyyyzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyyzz_1[i] * fi_acd_0 + g_x_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xx_0_xxxyyyzzz_0[i] = g_0_0_xxxyyyzzz_0[i] * fbe_0 - g_0_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyyzzz_1[i] * fi_acd_0 + g_x_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xx_0_xxxyyzzzz_0[i] = g_0_0_xxxyyzzzz_0[i] * fbe_0 - g_0_0_xxxyyzzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyyzzzz_1[i] * fi_acd_0 + g_x_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xx_0_xxxyzzzzz_0[i] = g_0_0_xxxyzzzzz_0[i] * fbe_0 - g_0_0_xxxyzzzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxyzzzzz_1[i] * fi_acd_0 + g_x_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xx_0_xxxzzzzzz_0[i] = g_0_0_xxxzzzzzz_0[i] * fbe_0 - g_0_0_xxxzzzzzz_1[i] * fz_be_0 + 3.0 * g_x_0_xxzzzzzz_1[i] * fi_acd_0 + g_x_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xx_0_xxyyyyyyy_0[i] = g_0_0_xxyyyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyyyy_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyyyy_1[i] * fi_acd_0 + g_x_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xx_0_xxyyyyyyz_0[i] = g_0_0_xxyyyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyyyz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyyyz_1[i] * fi_acd_0 + g_x_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xx_0_xxyyyyyzz_0[i] = g_0_0_xxyyyyyzz_0[i] * fbe_0 - g_0_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyyzz_1[i] * fi_acd_0 + g_x_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xx_0_xxyyyyzzz_0[i] = g_0_0_xxyyyyzzz_0[i] * fbe_0 - g_0_0_xxyyyyzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyyzzz_1[i] * fi_acd_0 + g_x_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xx_0_xxyyyzzzz_0[i] = g_0_0_xxyyyzzzz_0[i] * fbe_0 - g_0_0_xxyyyzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyyzzzz_1[i] * fi_acd_0 + g_x_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xx_0_xxyyzzzzz_0[i] = g_0_0_xxyyzzzzz_0[i] * fbe_0 - g_0_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyyzzzzz_1[i] * fi_acd_0 + g_x_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xx_0_xxyzzzzzz_0[i] = g_0_0_xxyzzzzzz_0[i] * fbe_0 - g_0_0_xxyzzzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xyzzzzzz_1[i] * fi_acd_0 + g_x_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xx_0_xxzzzzzzz_0[i] = g_0_0_xxzzzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzzzz_1[i] * fz_be_0 + 2.0 * g_x_0_xzzzzzzz_1[i] * fi_acd_0 + g_x_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xx_0_xyyyyyyyy_0[i] = g_0_0_xyyyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyyyy_1[i] * fz_be_0 + g_x_0_yyyyyyyy_1[i] * fi_acd_0 + g_x_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xx_0_xyyyyyyyz_0[i] = g_0_0_xyyyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyyyz_1[i] * fz_be_0 + g_x_0_yyyyyyyz_1[i] * fi_acd_0 + g_x_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xx_0_xyyyyyyzz_0[i] = g_0_0_xyyyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyyyzz_1[i] * fz_be_0 + g_x_0_yyyyyyzz_1[i] * fi_acd_0 + g_x_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xx_0_xyyyyyzzz_0[i] = g_0_0_xyyyyyzzz_0[i] * fbe_0 - g_0_0_xyyyyyzzz_1[i] * fz_be_0 + g_x_0_yyyyyzzz_1[i] * fi_acd_0 + g_x_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xx_0_xyyyyzzzz_0[i] = g_0_0_xyyyyzzzz_0[i] * fbe_0 - g_0_0_xyyyyzzzz_1[i] * fz_be_0 + g_x_0_yyyyzzzz_1[i] * fi_acd_0 + g_x_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xx_0_xyyyzzzzz_0[i] = g_0_0_xyyyzzzzz_0[i] * fbe_0 - g_0_0_xyyyzzzzz_1[i] * fz_be_0 + g_x_0_yyyzzzzz_1[i] * fi_acd_0 + g_x_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xx_0_xyyzzzzzz_0[i] = g_0_0_xyyzzzzzz_0[i] * fbe_0 - g_0_0_xyyzzzzzz_1[i] * fz_be_0 + g_x_0_yyzzzzzz_1[i] * fi_acd_0 + g_x_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xx_0_xyzzzzzzz_0[i] = g_0_0_xyzzzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzzzz_1[i] * fz_be_0 + g_x_0_yzzzzzzz_1[i] * fi_acd_0 + g_x_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xx_0_xzzzzzzzz_0[i] = g_0_0_xzzzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzzzz_1[i] * fz_be_0 + g_x_0_zzzzzzzz_1[i] * fi_acd_0 + g_x_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xx_0_yyyyyyyyy_0[i] = g_0_0_yyyyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyyyy_1[i] * fz_be_0 + g_x_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xx_0_yyyyyyyyz_0[i] = g_0_0_yyyyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyyyz_1[i] * fz_be_0 + g_x_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xx_0_yyyyyyyzz_0[i] = g_0_0_yyyyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyyyzz_1[i] * fz_be_0 + g_x_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xx_0_yyyyyyzzz_0[i] = g_0_0_yyyyyyzzz_0[i] * fbe_0 - g_0_0_yyyyyyzzz_1[i] * fz_be_0 + g_x_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xx_0_yyyyyzzzz_0[i] = g_0_0_yyyyyzzzz_0[i] * fbe_0 - g_0_0_yyyyyzzzz_1[i] * fz_be_0 + g_x_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xx_0_yyyyzzzzz_0[i] = g_0_0_yyyyzzzzz_0[i] * fbe_0 - g_0_0_yyyyzzzzz_1[i] * fz_be_0 + g_x_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xx_0_yyyzzzzzz_0[i] = g_0_0_yyyzzzzzz_0[i] * fbe_0 - g_0_0_yyyzzzzzz_1[i] * fz_be_0 + g_x_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xx_0_yyzzzzzzz_0[i] = g_0_0_yyzzzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzzzz_1[i] * fz_be_0 + g_x_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xx_0_yzzzzzzzz_0[i] = g_0_0_yzzzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzzzz_1[i] * fz_be_0 + g_x_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xx_0_zzzzzzzzz_0[i] = g_0_0_zzzzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzzzz_1[i] * fz_be_0 + g_x_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 55-110 components of targeted buffer : DSM

    auto g_xy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm + 55);

    auto g_xy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 56);

    auto g_xy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 57);

    auto g_xy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 58);

    auto g_xy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 59);

    auto g_xy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 60);

    auto g_xy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 61);

    auto g_xy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 62);

    auto g_xy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 63);

    auto g_xy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 64);

    auto g_xy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 65);

    auto g_xy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 66);

    auto g_xy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 67);

    auto g_xy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 68);

    auto g_xy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 69);

    auto g_xy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 70);

    auto g_xy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 71);

    auto g_xy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 72);

    auto g_xy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 73);

    auto g_xy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 74);

    auto g_xy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 75);

    auto g_xy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 76);

    auto g_xy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 77);

    auto g_xy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 78);

    auto g_xy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 79);

    auto g_xy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 80);

    auto g_xy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 81);

    auto g_xy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 82);

    auto g_xy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 83);

    auto g_xy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 84);

    auto g_xy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 85);

    auto g_xy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 86);

    auto g_xy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 87);

    auto g_xy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 88);

    auto g_xy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 89);

    auto g_xy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 90);

    auto g_xy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 91);

    auto g_xy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 92);

    auto g_xy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 93);

    auto g_xy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 94);

    auto g_xy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 95);

    auto g_xy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 96);

    auto g_xy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 97);

    auto g_xy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 98);

    auto g_xy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 99);

    auto g_xy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 100);

    auto g_xy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 101);

    auto g_xy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 102);

    auto g_xy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 103);

    auto g_xy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 104);

    auto g_xy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 105);

    auto g_xy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 106);

    auto g_xy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 107);

    auto g_xy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 108);

    auto g_xy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 109);

    #pragma omp simd aligned(g_x_0_xxxxxxxxx_1, g_x_0_xxxxxxxxz_1, g_x_0_xxxxxxxzz_1, g_x_0_xxxxxxzzz_1, g_x_0_xxxxxzzzz_1, g_x_0_xxxxzzzzz_1, g_x_0_xxxzzzzzz_1, g_x_0_xxzzzzzzz_1, g_x_0_xzzzzzzzz_1, g_xy_0_xxxxxxxxx_0, g_xy_0_xxxxxxxxy_0, g_xy_0_xxxxxxxxz_0, g_xy_0_xxxxxxxyy_0, g_xy_0_xxxxxxxyz_0, g_xy_0_xxxxxxxzz_0, g_xy_0_xxxxxxyyy_0, g_xy_0_xxxxxxyyz_0, g_xy_0_xxxxxxyzz_0, g_xy_0_xxxxxxzzz_0, g_xy_0_xxxxxyyyy_0, g_xy_0_xxxxxyyyz_0, g_xy_0_xxxxxyyzz_0, g_xy_0_xxxxxyzzz_0, g_xy_0_xxxxxzzzz_0, g_xy_0_xxxxyyyyy_0, g_xy_0_xxxxyyyyz_0, g_xy_0_xxxxyyyzz_0, g_xy_0_xxxxyyzzz_0, g_xy_0_xxxxyzzzz_0, g_xy_0_xxxxzzzzz_0, g_xy_0_xxxyyyyyy_0, g_xy_0_xxxyyyyyz_0, g_xy_0_xxxyyyyzz_0, g_xy_0_xxxyyyzzz_0, g_xy_0_xxxyyzzzz_0, g_xy_0_xxxyzzzzz_0, g_xy_0_xxxzzzzzz_0, g_xy_0_xxyyyyyyy_0, g_xy_0_xxyyyyyyz_0, g_xy_0_xxyyyyyzz_0, g_xy_0_xxyyyyzzz_0, g_xy_0_xxyyyzzzz_0, g_xy_0_xxyyzzzzz_0, g_xy_0_xxyzzzzzz_0, g_xy_0_xxzzzzzzz_0, g_xy_0_xyyyyyyyy_0, g_xy_0_xyyyyyyyz_0, g_xy_0_xyyyyyyzz_0, g_xy_0_xyyyyyzzz_0, g_xy_0_xyyyyzzzz_0, g_xy_0_xyyyzzzzz_0, g_xy_0_xyyzzzzzz_0, g_xy_0_xyzzzzzzz_0, g_xy_0_xzzzzzzzz_0, g_xy_0_yyyyyyyyy_0, g_xy_0_yyyyyyyyz_0, g_xy_0_yyyyyyyzz_0, g_xy_0_yyyyyyzzz_0, g_xy_0_yyyyyzzzz_0, g_xy_0_yyyyzzzzz_0, g_xy_0_yyyzzzzzz_0, g_xy_0_yyzzzzzzz_0, g_xy_0_yzzzzzzzz_0, g_xy_0_zzzzzzzzz_0, g_y_0_xxxxxxxxy_1, g_y_0_xxxxxxxy_1, g_y_0_xxxxxxxyy_1, g_y_0_xxxxxxxyz_1, g_y_0_xxxxxxyy_1, g_y_0_xxxxxxyyy_1, g_y_0_xxxxxxyyz_1, g_y_0_xxxxxxyz_1, g_y_0_xxxxxxyzz_1, g_y_0_xxxxxyyy_1, g_y_0_xxxxxyyyy_1, g_y_0_xxxxxyyyz_1, g_y_0_xxxxxyyz_1, g_y_0_xxxxxyyzz_1, g_y_0_xxxxxyzz_1, g_y_0_xxxxxyzzz_1, g_y_0_xxxxyyyy_1, g_y_0_xxxxyyyyy_1, g_y_0_xxxxyyyyz_1, g_y_0_xxxxyyyz_1, g_y_0_xxxxyyyzz_1, g_y_0_xxxxyyzz_1, g_y_0_xxxxyyzzz_1, g_y_0_xxxxyzzz_1, g_y_0_xxxxyzzzz_1, g_y_0_xxxyyyyy_1, g_y_0_xxxyyyyyy_1, g_y_0_xxxyyyyyz_1, g_y_0_xxxyyyyz_1, g_y_0_xxxyyyyzz_1, g_y_0_xxxyyyzz_1, g_y_0_xxxyyyzzz_1, g_y_0_xxxyyzzz_1, g_y_0_xxxyyzzzz_1, g_y_0_xxxyzzzz_1, g_y_0_xxxyzzzzz_1, g_y_0_xxyyyyyy_1, g_y_0_xxyyyyyyy_1, g_y_0_xxyyyyyyz_1, g_y_0_xxyyyyyz_1, g_y_0_xxyyyyyzz_1, g_y_0_xxyyyyzz_1, g_y_0_xxyyyyzzz_1, g_y_0_xxyyyzzz_1, g_y_0_xxyyyzzzz_1, g_y_0_xxyyzzzz_1, g_y_0_xxyyzzzzz_1, g_y_0_xxyzzzzz_1, g_y_0_xxyzzzzzz_1, g_y_0_xyyyyyyy_1, g_y_0_xyyyyyyyy_1, g_y_0_xyyyyyyyz_1, g_y_0_xyyyyyyz_1, g_y_0_xyyyyyyzz_1, g_y_0_xyyyyyzz_1, g_y_0_xyyyyyzzz_1, g_y_0_xyyyyzzz_1, g_y_0_xyyyyzzzz_1, g_y_0_xyyyzzzz_1, g_y_0_xyyyzzzzz_1, g_y_0_xyyzzzzz_1, g_y_0_xyyzzzzzz_1, g_y_0_xyzzzzzz_1, g_y_0_xyzzzzzzz_1, g_y_0_yyyyyyyy_1, g_y_0_yyyyyyyyy_1, g_y_0_yyyyyyyyz_1, g_y_0_yyyyyyyz_1, g_y_0_yyyyyyyzz_1, g_y_0_yyyyyyzz_1, g_y_0_yyyyyyzzz_1, g_y_0_yyyyyzzz_1, g_y_0_yyyyyzzzz_1, g_y_0_yyyyzzzz_1, g_y_0_yyyyzzzzz_1, g_y_0_yyyzzzzz_1, g_y_0_yyyzzzzzz_1, g_y_0_yyzzzzzz_1, g_y_0_yyzzzzzzz_1, g_y_0_yzzzzzzz_1, g_y_0_yzzzzzzzz_1, g_y_0_zzzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xxxxxxxxx_0[i] = g_x_0_xxxxxxxxx_1[i] * wa_y[i];

        g_xy_0_xxxxxxxxy_0[i] = 8.0 * g_y_0_xxxxxxxy_1[i] * fi_acd_0 + g_y_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xy_0_xxxxxxxxz_0[i] = g_x_0_xxxxxxxxz_1[i] * wa_y[i];

        g_xy_0_xxxxxxxyy_0[i] = 7.0 * g_y_0_xxxxxxyy_1[i] * fi_acd_0 + g_y_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xy_0_xxxxxxxyz_0[i] = 7.0 * g_y_0_xxxxxxyz_1[i] * fi_acd_0 + g_y_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xy_0_xxxxxxxzz_0[i] = g_x_0_xxxxxxxzz_1[i] * wa_y[i];

        g_xy_0_xxxxxxyyy_0[i] = 6.0 * g_y_0_xxxxxyyy_1[i] * fi_acd_0 + g_y_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xy_0_xxxxxxyyz_0[i] = 6.0 * g_y_0_xxxxxyyz_1[i] * fi_acd_0 + g_y_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xy_0_xxxxxxyzz_0[i] = 6.0 * g_y_0_xxxxxyzz_1[i] * fi_acd_0 + g_y_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xy_0_xxxxxxzzz_0[i] = g_x_0_xxxxxxzzz_1[i] * wa_y[i];

        g_xy_0_xxxxxyyyy_0[i] = 5.0 * g_y_0_xxxxyyyy_1[i] * fi_acd_0 + g_y_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xy_0_xxxxxyyyz_0[i] = 5.0 * g_y_0_xxxxyyyz_1[i] * fi_acd_0 + g_y_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xy_0_xxxxxyyzz_0[i] = 5.0 * g_y_0_xxxxyyzz_1[i] * fi_acd_0 + g_y_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xy_0_xxxxxyzzz_0[i] = 5.0 * g_y_0_xxxxyzzz_1[i] * fi_acd_0 + g_y_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xy_0_xxxxxzzzz_0[i] = g_x_0_xxxxxzzzz_1[i] * wa_y[i];

        g_xy_0_xxxxyyyyy_0[i] = 4.0 * g_y_0_xxxyyyyy_1[i] * fi_acd_0 + g_y_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xy_0_xxxxyyyyz_0[i] = 4.0 * g_y_0_xxxyyyyz_1[i] * fi_acd_0 + g_y_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xy_0_xxxxyyyzz_0[i] = 4.0 * g_y_0_xxxyyyzz_1[i] * fi_acd_0 + g_y_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xy_0_xxxxyyzzz_0[i] = 4.0 * g_y_0_xxxyyzzz_1[i] * fi_acd_0 + g_y_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xy_0_xxxxyzzzz_0[i] = 4.0 * g_y_0_xxxyzzzz_1[i] * fi_acd_0 + g_y_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xy_0_xxxxzzzzz_0[i] = g_x_0_xxxxzzzzz_1[i] * wa_y[i];

        g_xy_0_xxxyyyyyy_0[i] = 3.0 * g_y_0_xxyyyyyy_1[i] * fi_acd_0 + g_y_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xy_0_xxxyyyyyz_0[i] = 3.0 * g_y_0_xxyyyyyz_1[i] * fi_acd_0 + g_y_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xy_0_xxxyyyyzz_0[i] = 3.0 * g_y_0_xxyyyyzz_1[i] * fi_acd_0 + g_y_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xy_0_xxxyyyzzz_0[i] = 3.0 * g_y_0_xxyyyzzz_1[i] * fi_acd_0 + g_y_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xy_0_xxxyyzzzz_0[i] = 3.0 * g_y_0_xxyyzzzz_1[i] * fi_acd_0 + g_y_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xy_0_xxxyzzzzz_0[i] = 3.0 * g_y_0_xxyzzzzz_1[i] * fi_acd_0 + g_y_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xy_0_xxxzzzzzz_0[i] = g_x_0_xxxzzzzzz_1[i] * wa_y[i];

        g_xy_0_xxyyyyyyy_0[i] = 2.0 * g_y_0_xyyyyyyy_1[i] * fi_acd_0 + g_y_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xy_0_xxyyyyyyz_0[i] = 2.0 * g_y_0_xyyyyyyz_1[i] * fi_acd_0 + g_y_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xy_0_xxyyyyyzz_0[i] = 2.0 * g_y_0_xyyyyyzz_1[i] * fi_acd_0 + g_y_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xy_0_xxyyyyzzz_0[i] = 2.0 * g_y_0_xyyyyzzz_1[i] * fi_acd_0 + g_y_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xy_0_xxyyyzzzz_0[i] = 2.0 * g_y_0_xyyyzzzz_1[i] * fi_acd_0 + g_y_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xy_0_xxyyzzzzz_0[i] = 2.0 * g_y_0_xyyzzzzz_1[i] * fi_acd_0 + g_y_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xy_0_xxyzzzzzz_0[i] = 2.0 * g_y_0_xyzzzzzz_1[i] * fi_acd_0 + g_y_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xy_0_xxzzzzzzz_0[i] = g_x_0_xxzzzzzzz_1[i] * wa_y[i];

        g_xy_0_xyyyyyyyy_0[i] = g_y_0_yyyyyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xy_0_xyyyyyyyz_0[i] = g_y_0_yyyyyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xy_0_xyyyyyyzz_0[i] = g_y_0_yyyyyyzz_1[i] * fi_acd_0 + g_y_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xy_0_xyyyyyzzz_0[i] = g_y_0_yyyyyzzz_1[i] * fi_acd_0 + g_y_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xy_0_xyyyyzzzz_0[i] = g_y_0_yyyyzzzz_1[i] * fi_acd_0 + g_y_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xy_0_xyyyzzzzz_0[i] = g_y_0_yyyzzzzz_1[i] * fi_acd_0 + g_y_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xy_0_xyyzzzzzz_0[i] = g_y_0_yyzzzzzz_1[i] * fi_acd_0 + g_y_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xy_0_xyzzzzzzz_0[i] = g_y_0_yzzzzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xy_0_xzzzzzzzz_0[i] = g_x_0_xzzzzzzzz_1[i] * wa_y[i];

        g_xy_0_yyyyyyyyy_0[i] = g_y_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xy_0_yyyyyyyyz_0[i] = g_y_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xy_0_yyyyyyyzz_0[i] = g_y_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xy_0_yyyyyyzzz_0[i] = g_y_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xy_0_yyyyyzzzz_0[i] = g_y_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xy_0_yyyyzzzzz_0[i] = g_y_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xy_0_yyyzzzzzz_0[i] = g_y_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xy_0_yyzzzzzzz_0[i] = g_y_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xy_0_yzzzzzzzz_0[i] = g_y_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xy_0_zzzzzzzzz_0[i] = g_y_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 110-165 components of targeted buffer : DSM

    auto g_xz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm + 110);

    auto g_xz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 111);

    auto g_xz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 112);

    auto g_xz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 113);

    auto g_xz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 114);

    auto g_xz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 115);

    auto g_xz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 116);

    auto g_xz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 117);

    auto g_xz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 118);

    auto g_xz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 119);

    auto g_xz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 120);

    auto g_xz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 121);

    auto g_xz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 122);

    auto g_xz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 123);

    auto g_xz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 124);

    auto g_xz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 125);

    auto g_xz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 126);

    auto g_xz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 127);

    auto g_xz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 128);

    auto g_xz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 129);

    auto g_xz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 130);

    auto g_xz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 131);

    auto g_xz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 132);

    auto g_xz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 133);

    auto g_xz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 134);

    auto g_xz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 135);

    auto g_xz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 136);

    auto g_xz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 137);

    auto g_xz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 138);

    auto g_xz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 139);

    auto g_xz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 140);

    auto g_xz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 141);

    auto g_xz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 142);

    auto g_xz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 143);

    auto g_xz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 144);

    auto g_xz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 145);

    auto g_xz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 146);

    auto g_xz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 147);

    auto g_xz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 148);

    auto g_xz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 149);

    auto g_xz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 150);

    auto g_xz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 151);

    auto g_xz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 152);

    auto g_xz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 153);

    auto g_xz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 154);

    auto g_xz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 155);

    auto g_xz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 156);

    auto g_xz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 157);

    auto g_xz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 158);

    auto g_xz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 159);

    auto g_xz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 160);

    auto g_xz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 161);

    auto g_xz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 162);

    auto g_xz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 163);

    auto g_xz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 164);

    #pragma omp simd aligned(g_x_0_xxxxxxxxx_1, g_x_0_xxxxxxxxy_1, g_x_0_xxxxxxxyy_1, g_x_0_xxxxxxyyy_1, g_x_0_xxxxxyyyy_1, g_x_0_xxxxyyyyy_1, g_x_0_xxxyyyyyy_1, g_x_0_xxyyyyyyy_1, g_x_0_xyyyyyyyy_1, g_xz_0_xxxxxxxxx_0, g_xz_0_xxxxxxxxy_0, g_xz_0_xxxxxxxxz_0, g_xz_0_xxxxxxxyy_0, g_xz_0_xxxxxxxyz_0, g_xz_0_xxxxxxxzz_0, g_xz_0_xxxxxxyyy_0, g_xz_0_xxxxxxyyz_0, g_xz_0_xxxxxxyzz_0, g_xz_0_xxxxxxzzz_0, g_xz_0_xxxxxyyyy_0, g_xz_0_xxxxxyyyz_0, g_xz_0_xxxxxyyzz_0, g_xz_0_xxxxxyzzz_0, g_xz_0_xxxxxzzzz_0, g_xz_0_xxxxyyyyy_0, g_xz_0_xxxxyyyyz_0, g_xz_0_xxxxyyyzz_0, g_xz_0_xxxxyyzzz_0, g_xz_0_xxxxyzzzz_0, g_xz_0_xxxxzzzzz_0, g_xz_0_xxxyyyyyy_0, g_xz_0_xxxyyyyyz_0, g_xz_0_xxxyyyyzz_0, g_xz_0_xxxyyyzzz_0, g_xz_0_xxxyyzzzz_0, g_xz_0_xxxyzzzzz_0, g_xz_0_xxxzzzzzz_0, g_xz_0_xxyyyyyyy_0, g_xz_0_xxyyyyyyz_0, g_xz_0_xxyyyyyzz_0, g_xz_0_xxyyyyzzz_0, g_xz_0_xxyyyzzzz_0, g_xz_0_xxyyzzzzz_0, g_xz_0_xxyzzzzzz_0, g_xz_0_xxzzzzzzz_0, g_xz_0_xyyyyyyyy_0, g_xz_0_xyyyyyyyz_0, g_xz_0_xyyyyyyzz_0, g_xz_0_xyyyyyzzz_0, g_xz_0_xyyyyzzzz_0, g_xz_0_xyyyzzzzz_0, g_xz_0_xyyzzzzzz_0, g_xz_0_xyzzzzzzz_0, g_xz_0_xzzzzzzzz_0, g_xz_0_yyyyyyyyy_0, g_xz_0_yyyyyyyyz_0, g_xz_0_yyyyyyyzz_0, g_xz_0_yyyyyyzzz_0, g_xz_0_yyyyyzzzz_0, g_xz_0_yyyyzzzzz_0, g_xz_0_yyyzzzzzz_0, g_xz_0_yyzzzzzzz_0, g_xz_0_yzzzzzzzz_0, g_xz_0_zzzzzzzzz_0, g_z_0_xxxxxxxxz_1, g_z_0_xxxxxxxyz_1, g_z_0_xxxxxxxz_1, g_z_0_xxxxxxxzz_1, g_z_0_xxxxxxyyz_1, g_z_0_xxxxxxyz_1, g_z_0_xxxxxxyzz_1, g_z_0_xxxxxxzz_1, g_z_0_xxxxxxzzz_1, g_z_0_xxxxxyyyz_1, g_z_0_xxxxxyyz_1, g_z_0_xxxxxyyzz_1, g_z_0_xxxxxyzz_1, g_z_0_xxxxxyzzz_1, g_z_0_xxxxxzzz_1, g_z_0_xxxxxzzzz_1, g_z_0_xxxxyyyyz_1, g_z_0_xxxxyyyz_1, g_z_0_xxxxyyyzz_1, g_z_0_xxxxyyzz_1, g_z_0_xxxxyyzzz_1, g_z_0_xxxxyzzz_1, g_z_0_xxxxyzzzz_1, g_z_0_xxxxzzzz_1, g_z_0_xxxxzzzzz_1, g_z_0_xxxyyyyyz_1, g_z_0_xxxyyyyz_1, g_z_0_xxxyyyyzz_1, g_z_0_xxxyyyzz_1, g_z_0_xxxyyyzzz_1, g_z_0_xxxyyzzz_1, g_z_0_xxxyyzzzz_1, g_z_0_xxxyzzzz_1, g_z_0_xxxyzzzzz_1, g_z_0_xxxzzzzz_1, g_z_0_xxxzzzzzz_1, g_z_0_xxyyyyyyz_1, g_z_0_xxyyyyyz_1, g_z_0_xxyyyyyzz_1, g_z_0_xxyyyyzz_1, g_z_0_xxyyyyzzz_1, g_z_0_xxyyyzzz_1, g_z_0_xxyyyzzzz_1, g_z_0_xxyyzzzz_1, g_z_0_xxyyzzzzz_1, g_z_0_xxyzzzzz_1, g_z_0_xxyzzzzzz_1, g_z_0_xxzzzzzz_1, g_z_0_xxzzzzzzz_1, g_z_0_xyyyyyyyz_1, g_z_0_xyyyyyyz_1, g_z_0_xyyyyyyzz_1, g_z_0_xyyyyyzz_1, g_z_0_xyyyyyzzz_1, g_z_0_xyyyyzzz_1, g_z_0_xyyyyzzzz_1, g_z_0_xyyyzzzz_1, g_z_0_xyyyzzzzz_1, g_z_0_xyyzzzzz_1, g_z_0_xyyzzzzzz_1, g_z_0_xyzzzzzz_1, g_z_0_xyzzzzzzz_1, g_z_0_xzzzzzzz_1, g_z_0_xzzzzzzzz_1, g_z_0_yyyyyyyyy_1, g_z_0_yyyyyyyyz_1, g_z_0_yyyyyyyz_1, g_z_0_yyyyyyyzz_1, g_z_0_yyyyyyzz_1, g_z_0_yyyyyyzzz_1, g_z_0_yyyyyzzz_1, g_z_0_yyyyyzzzz_1, g_z_0_yyyyzzzz_1, g_z_0_yyyyzzzzz_1, g_z_0_yyyzzzzz_1, g_z_0_yyyzzzzzz_1, g_z_0_yyzzzzzz_1, g_z_0_yyzzzzzzz_1, g_z_0_yzzzzzzz_1, g_z_0_yzzzzzzzz_1, g_z_0_zzzzzzzz_1, g_z_0_zzzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xxxxxxxxx_0[i] = g_x_0_xxxxxxxxx_1[i] * wa_z[i];

        g_xz_0_xxxxxxxxy_0[i] = g_x_0_xxxxxxxxy_1[i] * wa_z[i];

        g_xz_0_xxxxxxxxz_0[i] = 8.0 * g_z_0_xxxxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xz_0_xxxxxxxyy_0[i] = g_x_0_xxxxxxxyy_1[i] * wa_z[i];

        g_xz_0_xxxxxxxyz_0[i] = 7.0 * g_z_0_xxxxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xz_0_xxxxxxxzz_0[i] = 7.0 * g_z_0_xxxxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xz_0_xxxxxxyyy_0[i] = g_x_0_xxxxxxyyy_1[i] * wa_z[i];

        g_xz_0_xxxxxxyyz_0[i] = 6.0 * g_z_0_xxxxxyyz_1[i] * fi_acd_0 + g_z_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xz_0_xxxxxxyzz_0[i] = 6.0 * g_z_0_xxxxxyzz_1[i] * fi_acd_0 + g_z_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xz_0_xxxxxxzzz_0[i] = 6.0 * g_z_0_xxxxxzzz_1[i] * fi_acd_0 + g_z_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xz_0_xxxxxyyyy_0[i] = g_x_0_xxxxxyyyy_1[i] * wa_z[i];

        g_xz_0_xxxxxyyyz_0[i] = 5.0 * g_z_0_xxxxyyyz_1[i] * fi_acd_0 + g_z_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xz_0_xxxxxyyzz_0[i] = 5.0 * g_z_0_xxxxyyzz_1[i] * fi_acd_0 + g_z_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xz_0_xxxxxyzzz_0[i] = 5.0 * g_z_0_xxxxyzzz_1[i] * fi_acd_0 + g_z_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xz_0_xxxxxzzzz_0[i] = 5.0 * g_z_0_xxxxzzzz_1[i] * fi_acd_0 + g_z_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xz_0_xxxxyyyyy_0[i] = g_x_0_xxxxyyyyy_1[i] * wa_z[i];

        g_xz_0_xxxxyyyyz_0[i] = 4.0 * g_z_0_xxxyyyyz_1[i] * fi_acd_0 + g_z_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xz_0_xxxxyyyzz_0[i] = 4.0 * g_z_0_xxxyyyzz_1[i] * fi_acd_0 + g_z_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xz_0_xxxxyyzzz_0[i] = 4.0 * g_z_0_xxxyyzzz_1[i] * fi_acd_0 + g_z_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xz_0_xxxxyzzzz_0[i] = 4.0 * g_z_0_xxxyzzzz_1[i] * fi_acd_0 + g_z_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xz_0_xxxxzzzzz_0[i] = 4.0 * g_z_0_xxxzzzzz_1[i] * fi_acd_0 + g_z_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xz_0_xxxyyyyyy_0[i] = g_x_0_xxxyyyyyy_1[i] * wa_z[i];

        g_xz_0_xxxyyyyyz_0[i] = 3.0 * g_z_0_xxyyyyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xz_0_xxxyyyyzz_0[i] = 3.0 * g_z_0_xxyyyyzz_1[i] * fi_acd_0 + g_z_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xz_0_xxxyyyzzz_0[i] = 3.0 * g_z_0_xxyyyzzz_1[i] * fi_acd_0 + g_z_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xz_0_xxxyyzzzz_0[i] = 3.0 * g_z_0_xxyyzzzz_1[i] * fi_acd_0 + g_z_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xz_0_xxxyzzzzz_0[i] = 3.0 * g_z_0_xxyzzzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xz_0_xxxzzzzzz_0[i] = 3.0 * g_z_0_xxzzzzzz_1[i] * fi_acd_0 + g_z_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xz_0_xxyyyyyyy_0[i] = g_x_0_xxyyyyyyy_1[i] * wa_z[i];

        g_xz_0_xxyyyyyyz_0[i] = 2.0 * g_z_0_xyyyyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xz_0_xxyyyyyzz_0[i] = 2.0 * g_z_0_xyyyyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xz_0_xxyyyyzzz_0[i] = 2.0 * g_z_0_xyyyyzzz_1[i] * fi_acd_0 + g_z_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xz_0_xxyyyzzzz_0[i] = 2.0 * g_z_0_xyyyzzzz_1[i] * fi_acd_0 + g_z_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xz_0_xxyyzzzzz_0[i] = 2.0 * g_z_0_xyyzzzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xz_0_xxyzzzzzz_0[i] = 2.0 * g_z_0_xyzzzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xz_0_xxzzzzzzz_0[i] = 2.0 * g_z_0_xzzzzzzz_1[i] * fi_acd_0 + g_z_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xz_0_xyyyyyyyy_0[i] = g_x_0_xyyyyyyyy_1[i] * wa_z[i];

        g_xz_0_xyyyyyyyz_0[i] = g_z_0_yyyyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xz_0_xyyyyyyzz_0[i] = g_z_0_yyyyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xz_0_xyyyyyzzz_0[i] = g_z_0_yyyyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xz_0_xyyyyzzzz_0[i] = g_z_0_yyyyzzzz_1[i] * fi_acd_0 + g_z_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xz_0_xyyyzzzzz_0[i] = g_z_0_yyyzzzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xz_0_xyyzzzzzz_0[i] = g_z_0_yyzzzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xz_0_xyzzzzzzz_0[i] = g_z_0_yzzzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xz_0_xzzzzzzzz_0[i] = g_z_0_zzzzzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xz_0_yyyyyyyyy_0[i] = g_z_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xz_0_yyyyyyyyz_0[i] = g_z_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xz_0_yyyyyyyzz_0[i] = g_z_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xz_0_yyyyyyzzz_0[i] = g_z_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xz_0_yyyyyzzzz_0[i] = g_z_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xz_0_yyyyzzzzz_0[i] = g_z_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xz_0_yyyzzzzzz_0[i] = g_z_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xz_0_yyzzzzzzz_0[i] = g_z_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xz_0_yzzzzzzzz_0[i] = g_z_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xz_0_zzzzzzzzz_0[i] = g_z_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 165-220 components of targeted buffer : DSM

    auto g_yy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm + 165);

    auto g_yy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 166);

    auto g_yy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 167);

    auto g_yy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 168);

    auto g_yy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 169);

    auto g_yy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 170);

    auto g_yy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 171);

    auto g_yy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 172);

    auto g_yy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 173);

    auto g_yy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 174);

    auto g_yy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 175);

    auto g_yy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 176);

    auto g_yy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 177);

    auto g_yy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 178);

    auto g_yy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 179);

    auto g_yy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 180);

    auto g_yy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 181);

    auto g_yy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 182);

    auto g_yy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 183);

    auto g_yy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 184);

    auto g_yy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 185);

    auto g_yy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 186);

    auto g_yy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 187);

    auto g_yy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 188);

    auto g_yy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 189);

    auto g_yy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 190);

    auto g_yy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 191);

    auto g_yy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 192);

    auto g_yy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 193);

    auto g_yy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 194);

    auto g_yy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 195);

    auto g_yy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 196);

    auto g_yy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 197);

    auto g_yy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 198);

    auto g_yy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 199);

    auto g_yy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 200);

    auto g_yy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 201);

    auto g_yy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 202);

    auto g_yy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 203);

    auto g_yy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 204);

    auto g_yy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 205);

    auto g_yy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 206);

    auto g_yy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 207);

    auto g_yy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 208);

    auto g_yy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 209);

    auto g_yy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 210);

    auto g_yy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 211);

    auto g_yy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 212);

    auto g_yy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 213);

    auto g_yy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 214);

    auto g_yy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 215);

    auto g_yy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 216);

    auto g_yy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 217);

    auto g_yy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 218);

    auto g_yy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 219);

    #pragma omp simd aligned(g_0_0_xxxxxxxxx_0, g_0_0_xxxxxxxxx_1, g_0_0_xxxxxxxxy_0, g_0_0_xxxxxxxxy_1, g_0_0_xxxxxxxxz_0, g_0_0_xxxxxxxxz_1, g_0_0_xxxxxxxyy_0, g_0_0_xxxxxxxyy_1, g_0_0_xxxxxxxyz_0, g_0_0_xxxxxxxyz_1, g_0_0_xxxxxxxzz_0, g_0_0_xxxxxxxzz_1, g_0_0_xxxxxxyyy_0, g_0_0_xxxxxxyyy_1, g_0_0_xxxxxxyyz_0, g_0_0_xxxxxxyyz_1, g_0_0_xxxxxxyzz_0, g_0_0_xxxxxxyzz_1, g_0_0_xxxxxxzzz_0, g_0_0_xxxxxxzzz_1, g_0_0_xxxxxyyyy_0, g_0_0_xxxxxyyyy_1, g_0_0_xxxxxyyyz_0, g_0_0_xxxxxyyyz_1, g_0_0_xxxxxyyzz_0, g_0_0_xxxxxyyzz_1, g_0_0_xxxxxyzzz_0, g_0_0_xxxxxyzzz_1, g_0_0_xxxxxzzzz_0, g_0_0_xxxxxzzzz_1, g_0_0_xxxxyyyyy_0, g_0_0_xxxxyyyyy_1, g_0_0_xxxxyyyyz_0, g_0_0_xxxxyyyyz_1, g_0_0_xxxxyyyzz_0, g_0_0_xxxxyyyzz_1, g_0_0_xxxxyyzzz_0, g_0_0_xxxxyyzzz_1, g_0_0_xxxxyzzzz_0, g_0_0_xxxxyzzzz_1, g_0_0_xxxxzzzzz_0, g_0_0_xxxxzzzzz_1, g_0_0_xxxyyyyyy_0, g_0_0_xxxyyyyyy_1, g_0_0_xxxyyyyyz_0, g_0_0_xxxyyyyyz_1, g_0_0_xxxyyyyzz_0, g_0_0_xxxyyyyzz_1, g_0_0_xxxyyyzzz_0, g_0_0_xxxyyyzzz_1, g_0_0_xxxyyzzzz_0, g_0_0_xxxyyzzzz_1, g_0_0_xxxyzzzzz_0, g_0_0_xxxyzzzzz_1, g_0_0_xxxzzzzzz_0, g_0_0_xxxzzzzzz_1, g_0_0_xxyyyyyyy_0, g_0_0_xxyyyyyyy_1, g_0_0_xxyyyyyyz_0, g_0_0_xxyyyyyyz_1, g_0_0_xxyyyyyzz_0, g_0_0_xxyyyyyzz_1, g_0_0_xxyyyyzzz_0, g_0_0_xxyyyyzzz_1, g_0_0_xxyyyzzzz_0, g_0_0_xxyyyzzzz_1, g_0_0_xxyyzzzzz_0, g_0_0_xxyyzzzzz_1, g_0_0_xxyzzzzzz_0, g_0_0_xxyzzzzzz_1, g_0_0_xxzzzzzzz_0, g_0_0_xxzzzzzzz_1, g_0_0_xyyyyyyyy_0, g_0_0_xyyyyyyyy_1, g_0_0_xyyyyyyyz_0, g_0_0_xyyyyyyyz_1, g_0_0_xyyyyyyzz_0, g_0_0_xyyyyyyzz_1, g_0_0_xyyyyyzzz_0, g_0_0_xyyyyyzzz_1, g_0_0_xyyyyzzzz_0, g_0_0_xyyyyzzzz_1, g_0_0_xyyyzzzzz_0, g_0_0_xyyyzzzzz_1, g_0_0_xyyzzzzzz_0, g_0_0_xyyzzzzzz_1, g_0_0_xyzzzzzzz_0, g_0_0_xyzzzzzzz_1, g_0_0_xzzzzzzzz_0, g_0_0_xzzzzzzzz_1, g_0_0_yyyyyyyyy_0, g_0_0_yyyyyyyyy_1, g_0_0_yyyyyyyyz_0, g_0_0_yyyyyyyyz_1, g_0_0_yyyyyyyzz_0, g_0_0_yyyyyyyzz_1, g_0_0_yyyyyyzzz_0, g_0_0_yyyyyyzzz_1, g_0_0_yyyyyzzzz_0, g_0_0_yyyyyzzzz_1, g_0_0_yyyyzzzzz_0, g_0_0_yyyyzzzzz_1, g_0_0_yyyzzzzzz_0, g_0_0_yyyzzzzzz_1, g_0_0_yyzzzzzzz_0, g_0_0_yyzzzzzzz_1, g_0_0_yzzzzzzzz_0, g_0_0_yzzzzzzzz_1, g_0_0_zzzzzzzzz_0, g_0_0_zzzzzzzzz_1, g_y_0_xxxxxxxx_1, g_y_0_xxxxxxxxx_1, g_y_0_xxxxxxxxy_1, g_y_0_xxxxxxxxz_1, g_y_0_xxxxxxxy_1, g_y_0_xxxxxxxyy_1, g_y_0_xxxxxxxyz_1, g_y_0_xxxxxxxz_1, g_y_0_xxxxxxxzz_1, g_y_0_xxxxxxyy_1, g_y_0_xxxxxxyyy_1, g_y_0_xxxxxxyyz_1, g_y_0_xxxxxxyz_1, g_y_0_xxxxxxyzz_1, g_y_0_xxxxxxzz_1, g_y_0_xxxxxxzzz_1, g_y_0_xxxxxyyy_1, g_y_0_xxxxxyyyy_1, g_y_0_xxxxxyyyz_1, g_y_0_xxxxxyyz_1, g_y_0_xxxxxyyzz_1, g_y_0_xxxxxyzz_1, g_y_0_xxxxxyzzz_1, g_y_0_xxxxxzzz_1, g_y_0_xxxxxzzzz_1, g_y_0_xxxxyyyy_1, g_y_0_xxxxyyyyy_1, g_y_0_xxxxyyyyz_1, g_y_0_xxxxyyyz_1, g_y_0_xxxxyyyzz_1, g_y_0_xxxxyyzz_1, g_y_0_xxxxyyzzz_1, g_y_0_xxxxyzzz_1, g_y_0_xxxxyzzzz_1, g_y_0_xxxxzzzz_1, g_y_0_xxxxzzzzz_1, g_y_0_xxxyyyyy_1, g_y_0_xxxyyyyyy_1, g_y_0_xxxyyyyyz_1, g_y_0_xxxyyyyz_1, g_y_0_xxxyyyyzz_1, g_y_0_xxxyyyzz_1, g_y_0_xxxyyyzzz_1, g_y_0_xxxyyzzz_1, g_y_0_xxxyyzzzz_1, g_y_0_xxxyzzzz_1, g_y_0_xxxyzzzzz_1, g_y_0_xxxzzzzz_1, g_y_0_xxxzzzzzz_1, g_y_0_xxyyyyyy_1, g_y_0_xxyyyyyyy_1, g_y_0_xxyyyyyyz_1, g_y_0_xxyyyyyz_1, g_y_0_xxyyyyyzz_1, g_y_0_xxyyyyzz_1, g_y_0_xxyyyyzzz_1, g_y_0_xxyyyzzz_1, g_y_0_xxyyyzzzz_1, g_y_0_xxyyzzzz_1, g_y_0_xxyyzzzzz_1, g_y_0_xxyzzzzz_1, g_y_0_xxyzzzzzz_1, g_y_0_xxzzzzzz_1, g_y_0_xxzzzzzzz_1, g_y_0_xyyyyyyy_1, g_y_0_xyyyyyyyy_1, g_y_0_xyyyyyyyz_1, g_y_0_xyyyyyyz_1, g_y_0_xyyyyyyzz_1, g_y_0_xyyyyyzz_1, g_y_0_xyyyyyzzz_1, g_y_0_xyyyyzzz_1, g_y_0_xyyyyzzzz_1, g_y_0_xyyyzzzz_1, g_y_0_xyyyzzzzz_1, g_y_0_xyyzzzzz_1, g_y_0_xyyzzzzzz_1, g_y_0_xyzzzzzz_1, g_y_0_xyzzzzzzz_1, g_y_0_xzzzzzzz_1, g_y_0_xzzzzzzzz_1, g_y_0_yyyyyyyy_1, g_y_0_yyyyyyyyy_1, g_y_0_yyyyyyyyz_1, g_y_0_yyyyyyyz_1, g_y_0_yyyyyyyzz_1, g_y_0_yyyyyyzz_1, g_y_0_yyyyyyzzz_1, g_y_0_yyyyyzzz_1, g_y_0_yyyyyzzzz_1, g_y_0_yyyyzzzz_1, g_y_0_yyyyzzzzz_1, g_y_0_yyyzzzzz_1, g_y_0_yyyzzzzzz_1, g_y_0_yyzzzzzz_1, g_y_0_yyzzzzzzz_1, g_y_0_yzzzzzzz_1, g_y_0_yzzzzzzzz_1, g_y_0_zzzzzzzz_1, g_y_0_zzzzzzzzz_1, g_yy_0_xxxxxxxxx_0, g_yy_0_xxxxxxxxy_0, g_yy_0_xxxxxxxxz_0, g_yy_0_xxxxxxxyy_0, g_yy_0_xxxxxxxyz_0, g_yy_0_xxxxxxxzz_0, g_yy_0_xxxxxxyyy_0, g_yy_0_xxxxxxyyz_0, g_yy_0_xxxxxxyzz_0, g_yy_0_xxxxxxzzz_0, g_yy_0_xxxxxyyyy_0, g_yy_0_xxxxxyyyz_0, g_yy_0_xxxxxyyzz_0, g_yy_0_xxxxxyzzz_0, g_yy_0_xxxxxzzzz_0, g_yy_0_xxxxyyyyy_0, g_yy_0_xxxxyyyyz_0, g_yy_0_xxxxyyyzz_0, g_yy_0_xxxxyyzzz_0, g_yy_0_xxxxyzzzz_0, g_yy_0_xxxxzzzzz_0, g_yy_0_xxxyyyyyy_0, g_yy_0_xxxyyyyyz_0, g_yy_0_xxxyyyyzz_0, g_yy_0_xxxyyyzzz_0, g_yy_0_xxxyyzzzz_0, g_yy_0_xxxyzzzzz_0, g_yy_0_xxxzzzzzz_0, g_yy_0_xxyyyyyyy_0, g_yy_0_xxyyyyyyz_0, g_yy_0_xxyyyyyzz_0, g_yy_0_xxyyyyzzz_0, g_yy_0_xxyyyzzzz_0, g_yy_0_xxyyzzzzz_0, g_yy_0_xxyzzzzzz_0, g_yy_0_xxzzzzzzz_0, g_yy_0_xyyyyyyyy_0, g_yy_0_xyyyyyyyz_0, g_yy_0_xyyyyyyzz_0, g_yy_0_xyyyyyzzz_0, g_yy_0_xyyyyzzzz_0, g_yy_0_xyyyzzzzz_0, g_yy_0_xyyzzzzzz_0, g_yy_0_xyzzzzzzz_0, g_yy_0_xzzzzzzzz_0, g_yy_0_yyyyyyyyy_0, g_yy_0_yyyyyyyyz_0, g_yy_0_yyyyyyyzz_0, g_yy_0_yyyyyyzzz_0, g_yy_0_yyyyyzzzz_0, g_yy_0_yyyyzzzzz_0, g_yy_0_yyyzzzzzz_0, g_yy_0_yyzzzzzzz_0, g_yy_0_yzzzzzzzz_0, g_yy_0_zzzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xxxxxxxxx_0[i] = g_0_0_xxxxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxxxx_1[i] * fz_be_0 + g_y_0_xxxxxxxxx_1[i] * wa_y[i];

        g_yy_0_xxxxxxxxy_0[i] = g_0_0_xxxxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxxxy_1[i] * fz_be_0 + g_y_0_xxxxxxxx_1[i] * fi_acd_0 + g_y_0_xxxxxxxxy_1[i] * wa_y[i];

        g_yy_0_xxxxxxxxz_0[i] = g_0_0_xxxxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxxxz_1[i] * fz_be_0 + g_y_0_xxxxxxxxz_1[i] * wa_y[i];

        g_yy_0_xxxxxxxyy_0[i] = g_0_0_xxxxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxxxyy_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxxxxy_1[i] * fi_acd_0 + g_y_0_xxxxxxxyy_1[i] * wa_y[i];

        g_yy_0_xxxxxxxyz_0[i] = g_0_0_xxxxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxxxyz_1[i] * fz_be_0 + g_y_0_xxxxxxxz_1[i] * fi_acd_0 + g_y_0_xxxxxxxyz_1[i] * wa_y[i];

        g_yy_0_xxxxxxxzz_0[i] = g_0_0_xxxxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxxxzz_1[i] * fz_be_0 + g_y_0_xxxxxxxzz_1[i] * wa_y[i];

        g_yy_0_xxxxxxyyy_0[i] = g_0_0_xxxxxxyyy_0[i] * fbe_0 - g_0_0_xxxxxxyyy_1[i] * fz_be_0 + 3.0 * g_y_0_xxxxxxyy_1[i] * fi_acd_0 + g_y_0_xxxxxxyyy_1[i] * wa_y[i];

        g_yy_0_xxxxxxyyz_0[i] = g_0_0_xxxxxxyyz_0[i] * fbe_0 - g_0_0_xxxxxxyyz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxxxyz_1[i] * fi_acd_0 + g_y_0_xxxxxxyyz_1[i] * wa_y[i];

        g_yy_0_xxxxxxyzz_0[i] = g_0_0_xxxxxxyzz_0[i] * fbe_0 - g_0_0_xxxxxxyzz_1[i] * fz_be_0 + g_y_0_xxxxxxzz_1[i] * fi_acd_0 + g_y_0_xxxxxxyzz_1[i] * wa_y[i];

        g_yy_0_xxxxxxzzz_0[i] = g_0_0_xxxxxxzzz_0[i] * fbe_0 - g_0_0_xxxxxxzzz_1[i] * fz_be_0 + g_y_0_xxxxxxzzz_1[i] * wa_y[i];

        g_yy_0_xxxxxyyyy_0[i] = g_0_0_xxxxxyyyy_0[i] * fbe_0 - g_0_0_xxxxxyyyy_1[i] * fz_be_0 + 4.0 * g_y_0_xxxxxyyy_1[i] * fi_acd_0 + g_y_0_xxxxxyyyy_1[i] * wa_y[i];

        g_yy_0_xxxxxyyyz_0[i] = g_0_0_xxxxxyyyz_0[i] * fbe_0 - g_0_0_xxxxxyyyz_1[i] * fz_be_0 + 3.0 * g_y_0_xxxxxyyz_1[i] * fi_acd_0 + g_y_0_xxxxxyyyz_1[i] * wa_y[i];

        g_yy_0_xxxxxyyzz_0[i] = g_0_0_xxxxxyyzz_0[i] * fbe_0 - g_0_0_xxxxxyyzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxxyzz_1[i] * fi_acd_0 + g_y_0_xxxxxyyzz_1[i] * wa_y[i];

        g_yy_0_xxxxxyzzz_0[i] = g_0_0_xxxxxyzzz_0[i] * fbe_0 - g_0_0_xxxxxyzzz_1[i] * fz_be_0 + g_y_0_xxxxxzzz_1[i] * fi_acd_0 + g_y_0_xxxxxyzzz_1[i] * wa_y[i];

        g_yy_0_xxxxxzzzz_0[i] = g_0_0_xxxxxzzzz_0[i] * fbe_0 - g_0_0_xxxxxzzzz_1[i] * fz_be_0 + g_y_0_xxxxxzzzz_1[i] * wa_y[i];

        g_yy_0_xxxxyyyyy_0[i] = g_0_0_xxxxyyyyy_0[i] * fbe_0 - g_0_0_xxxxyyyyy_1[i] * fz_be_0 + 5.0 * g_y_0_xxxxyyyy_1[i] * fi_acd_0 + g_y_0_xxxxyyyyy_1[i] * wa_y[i];

        g_yy_0_xxxxyyyyz_0[i] = g_0_0_xxxxyyyyz_0[i] * fbe_0 - g_0_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_y_0_xxxxyyyz_1[i] * fi_acd_0 + g_y_0_xxxxyyyyz_1[i] * wa_y[i];

        g_yy_0_xxxxyyyzz_0[i] = g_0_0_xxxxyyyzz_0[i] * fbe_0 - g_0_0_xxxxyyyzz_1[i] * fz_be_0 + 3.0 * g_y_0_xxxxyyzz_1[i] * fi_acd_0 + g_y_0_xxxxyyyzz_1[i] * wa_y[i];

        g_yy_0_xxxxyyzzz_0[i] = g_0_0_xxxxyyzzz_0[i] * fbe_0 - g_0_0_xxxxyyzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxxyzzz_1[i] * fi_acd_0 + g_y_0_xxxxyyzzz_1[i] * wa_y[i];

        g_yy_0_xxxxyzzzz_0[i] = g_0_0_xxxxyzzzz_0[i] * fbe_0 - g_0_0_xxxxyzzzz_1[i] * fz_be_0 + g_y_0_xxxxzzzz_1[i] * fi_acd_0 + g_y_0_xxxxyzzzz_1[i] * wa_y[i];

        g_yy_0_xxxxzzzzz_0[i] = g_0_0_xxxxzzzzz_0[i] * fbe_0 - g_0_0_xxxxzzzzz_1[i] * fz_be_0 + g_y_0_xxxxzzzzz_1[i] * wa_y[i];

        g_yy_0_xxxyyyyyy_0[i] = g_0_0_xxxyyyyyy_0[i] * fbe_0 - g_0_0_xxxyyyyyy_1[i] * fz_be_0 + 6.0 * g_y_0_xxxyyyyy_1[i] * fi_acd_0 + g_y_0_xxxyyyyyy_1[i] * wa_y[i];

        g_yy_0_xxxyyyyyz_0[i] = g_0_0_xxxyyyyyz_0[i] * fbe_0 - g_0_0_xxxyyyyyz_1[i] * fz_be_0 + 5.0 * g_y_0_xxxyyyyz_1[i] * fi_acd_0 + g_y_0_xxxyyyyyz_1[i] * wa_y[i];

        g_yy_0_xxxyyyyzz_0[i] = g_0_0_xxxyyyyzz_0[i] * fbe_0 - g_0_0_xxxyyyyzz_1[i] * fz_be_0 + 4.0 * g_y_0_xxxyyyzz_1[i] * fi_acd_0 + g_y_0_xxxyyyyzz_1[i] * wa_y[i];

        g_yy_0_xxxyyyzzz_0[i] = g_0_0_xxxyyyzzz_0[i] * fbe_0 - g_0_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_y_0_xxxyyzzz_1[i] * fi_acd_0 + g_y_0_xxxyyyzzz_1[i] * wa_y[i];

        g_yy_0_xxxyyzzzz_0[i] = g_0_0_xxxyyzzzz_0[i] * fbe_0 - g_0_0_xxxyyzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxxyzzzz_1[i] * fi_acd_0 + g_y_0_xxxyyzzzz_1[i] * wa_y[i];

        g_yy_0_xxxyzzzzz_0[i] = g_0_0_xxxyzzzzz_0[i] * fbe_0 - g_0_0_xxxyzzzzz_1[i] * fz_be_0 + g_y_0_xxxzzzzz_1[i] * fi_acd_0 + g_y_0_xxxyzzzzz_1[i] * wa_y[i];

        g_yy_0_xxxzzzzzz_0[i] = g_0_0_xxxzzzzzz_0[i] * fbe_0 - g_0_0_xxxzzzzzz_1[i] * fz_be_0 + g_y_0_xxxzzzzzz_1[i] * wa_y[i];

        g_yy_0_xxyyyyyyy_0[i] = g_0_0_xxyyyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyyyy_1[i] * fz_be_0 + 7.0 * g_y_0_xxyyyyyy_1[i] * fi_acd_0 + g_y_0_xxyyyyyyy_1[i] * wa_y[i];

        g_yy_0_xxyyyyyyz_0[i] = g_0_0_xxyyyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyyyz_1[i] * fz_be_0 + 6.0 * g_y_0_xxyyyyyz_1[i] * fi_acd_0 + g_y_0_xxyyyyyyz_1[i] * wa_y[i];

        g_yy_0_xxyyyyyzz_0[i] = g_0_0_xxyyyyyzz_0[i] * fbe_0 - g_0_0_xxyyyyyzz_1[i] * fz_be_0 + 5.0 * g_y_0_xxyyyyzz_1[i] * fi_acd_0 + g_y_0_xxyyyyyzz_1[i] * wa_y[i];

        g_yy_0_xxyyyyzzz_0[i] = g_0_0_xxyyyyzzz_0[i] * fbe_0 - g_0_0_xxyyyyzzz_1[i] * fz_be_0 + 4.0 * g_y_0_xxyyyzzz_1[i] * fi_acd_0 + g_y_0_xxyyyyzzz_1[i] * wa_y[i];

        g_yy_0_xxyyyzzzz_0[i] = g_0_0_xxyyyzzzz_0[i] * fbe_0 - g_0_0_xxyyyzzzz_1[i] * fz_be_0 + 3.0 * g_y_0_xxyyzzzz_1[i] * fi_acd_0 + g_y_0_xxyyyzzzz_1[i] * wa_y[i];

        g_yy_0_xxyyzzzzz_0[i] = g_0_0_xxyyzzzzz_0[i] * fbe_0 - g_0_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xxyzzzzz_1[i] * fi_acd_0 + g_y_0_xxyyzzzzz_1[i] * wa_y[i];

        g_yy_0_xxyzzzzzz_0[i] = g_0_0_xxyzzzzzz_0[i] * fbe_0 - g_0_0_xxyzzzzzz_1[i] * fz_be_0 + g_y_0_xxzzzzzz_1[i] * fi_acd_0 + g_y_0_xxyzzzzzz_1[i] * wa_y[i];

        g_yy_0_xxzzzzzzz_0[i] = g_0_0_xxzzzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzzzz_1[i] * fz_be_0 + g_y_0_xxzzzzzzz_1[i] * wa_y[i];

        g_yy_0_xyyyyyyyy_0[i] = g_0_0_xyyyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyyyy_1[i] * fz_be_0 + 8.0 * g_y_0_xyyyyyyy_1[i] * fi_acd_0 + g_y_0_xyyyyyyyy_1[i] * wa_y[i];

        g_yy_0_xyyyyyyyz_0[i] = g_0_0_xyyyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyyyz_1[i] * fz_be_0 + 7.0 * g_y_0_xyyyyyyz_1[i] * fi_acd_0 + g_y_0_xyyyyyyyz_1[i] * wa_y[i];

        g_yy_0_xyyyyyyzz_0[i] = g_0_0_xyyyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyyyzz_1[i] * fz_be_0 + 6.0 * g_y_0_xyyyyyzz_1[i] * fi_acd_0 + g_y_0_xyyyyyyzz_1[i] * wa_y[i];

        g_yy_0_xyyyyyzzz_0[i] = g_0_0_xyyyyyzzz_0[i] * fbe_0 - g_0_0_xyyyyyzzz_1[i] * fz_be_0 + 5.0 * g_y_0_xyyyyzzz_1[i] * fi_acd_0 + g_y_0_xyyyyyzzz_1[i] * wa_y[i];

        g_yy_0_xyyyyzzzz_0[i] = g_0_0_xyyyyzzzz_0[i] * fbe_0 - g_0_0_xyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_y_0_xyyyzzzz_1[i] * fi_acd_0 + g_y_0_xyyyyzzzz_1[i] * wa_y[i];

        g_yy_0_xyyyzzzzz_0[i] = g_0_0_xyyyzzzzz_0[i] * fbe_0 - g_0_0_xyyyzzzzz_1[i] * fz_be_0 + 3.0 * g_y_0_xyyzzzzz_1[i] * fi_acd_0 + g_y_0_xyyyzzzzz_1[i] * wa_y[i];

        g_yy_0_xyyzzzzzz_0[i] = g_0_0_xyyzzzzzz_0[i] * fbe_0 - g_0_0_xyyzzzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_xyzzzzzz_1[i] * fi_acd_0 + g_y_0_xyyzzzzzz_1[i] * wa_y[i];

        g_yy_0_xyzzzzzzz_0[i] = g_0_0_xyzzzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzzzz_1[i] * fz_be_0 + g_y_0_xzzzzzzz_1[i] * fi_acd_0 + g_y_0_xyzzzzzzz_1[i] * wa_y[i];

        g_yy_0_xzzzzzzzz_0[i] = g_0_0_xzzzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzzzz_1[i] * fz_be_0 + g_y_0_xzzzzzzzz_1[i] * wa_y[i];

        g_yy_0_yyyyyyyyy_0[i] = g_0_0_yyyyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyyyy_1[i] * fz_be_0 + 9.0 * g_y_0_yyyyyyyy_1[i] * fi_acd_0 + g_y_0_yyyyyyyyy_1[i] * wa_y[i];

        g_yy_0_yyyyyyyyz_0[i] = g_0_0_yyyyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyyyz_1[i] * fz_be_0 + 8.0 * g_y_0_yyyyyyyz_1[i] * fi_acd_0 + g_y_0_yyyyyyyyz_1[i] * wa_y[i];

        g_yy_0_yyyyyyyzz_0[i] = g_0_0_yyyyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyyyzz_1[i] * fz_be_0 + 7.0 * g_y_0_yyyyyyzz_1[i] * fi_acd_0 + g_y_0_yyyyyyyzz_1[i] * wa_y[i];

        g_yy_0_yyyyyyzzz_0[i] = g_0_0_yyyyyyzzz_0[i] * fbe_0 - g_0_0_yyyyyyzzz_1[i] * fz_be_0 + 6.0 * g_y_0_yyyyyzzz_1[i] * fi_acd_0 + g_y_0_yyyyyyzzz_1[i] * wa_y[i];

        g_yy_0_yyyyyzzzz_0[i] = g_0_0_yyyyyzzzz_0[i] * fbe_0 - g_0_0_yyyyyzzzz_1[i] * fz_be_0 + 5.0 * g_y_0_yyyyzzzz_1[i] * fi_acd_0 + g_y_0_yyyyyzzzz_1[i] * wa_y[i];

        g_yy_0_yyyyzzzzz_0[i] = g_0_0_yyyyzzzzz_0[i] * fbe_0 - g_0_0_yyyyzzzzz_1[i] * fz_be_0 + 4.0 * g_y_0_yyyzzzzz_1[i] * fi_acd_0 + g_y_0_yyyyzzzzz_1[i] * wa_y[i];

        g_yy_0_yyyzzzzzz_0[i] = g_0_0_yyyzzzzzz_0[i] * fbe_0 - g_0_0_yyyzzzzzz_1[i] * fz_be_0 + 3.0 * g_y_0_yyzzzzzz_1[i] * fi_acd_0 + g_y_0_yyyzzzzzz_1[i] * wa_y[i];

        g_yy_0_yyzzzzzzz_0[i] = g_0_0_yyzzzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzzzz_1[i] * fz_be_0 + 2.0 * g_y_0_yzzzzzzz_1[i] * fi_acd_0 + g_y_0_yyzzzzzzz_1[i] * wa_y[i];

        g_yy_0_yzzzzzzzz_0[i] = g_0_0_yzzzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzzzz_1[i] * fz_be_0 + g_y_0_zzzzzzzz_1[i] * fi_acd_0 + g_y_0_yzzzzzzzz_1[i] * wa_y[i];

        g_yy_0_zzzzzzzzz_0[i] = g_0_0_zzzzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzzzz_1[i] * fz_be_0 + g_y_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 220-275 components of targeted buffer : DSM

    auto g_yz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm + 220);

    auto g_yz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 221);

    auto g_yz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 222);

    auto g_yz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 223);

    auto g_yz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 224);

    auto g_yz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 225);

    auto g_yz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 226);

    auto g_yz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 227);

    auto g_yz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 228);

    auto g_yz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 229);

    auto g_yz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 230);

    auto g_yz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 231);

    auto g_yz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 232);

    auto g_yz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 233);

    auto g_yz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 234);

    auto g_yz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 235);

    auto g_yz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 236);

    auto g_yz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 237);

    auto g_yz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 238);

    auto g_yz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 239);

    auto g_yz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 240);

    auto g_yz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 241);

    auto g_yz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 242);

    auto g_yz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 243);

    auto g_yz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 244);

    auto g_yz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 245);

    auto g_yz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 246);

    auto g_yz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 247);

    auto g_yz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 248);

    auto g_yz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 249);

    auto g_yz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 250);

    auto g_yz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 251);

    auto g_yz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 252);

    auto g_yz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 253);

    auto g_yz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 254);

    auto g_yz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 255);

    auto g_yz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 256);

    auto g_yz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 257);

    auto g_yz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 258);

    auto g_yz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 259);

    auto g_yz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 260);

    auto g_yz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 261);

    auto g_yz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 262);

    auto g_yz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 263);

    auto g_yz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 264);

    auto g_yz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 265);

    auto g_yz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 266);

    auto g_yz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 267);

    auto g_yz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 268);

    auto g_yz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 269);

    auto g_yz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 270);

    auto g_yz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 271);

    auto g_yz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 272);

    auto g_yz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 273);

    auto g_yz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 274);

    #pragma omp simd aligned(g_y_0_xxxxxxxxy_1, g_y_0_xxxxxxxyy_1, g_y_0_xxxxxxyyy_1, g_y_0_xxxxxyyyy_1, g_y_0_xxxxyyyyy_1, g_y_0_xxxyyyyyy_1, g_y_0_xxyyyyyyy_1, g_y_0_xyyyyyyyy_1, g_y_0_yyyyyyyyy_1, g_yz_0_xxxxxxxxx_0, g_yz_0_xxxxxxxxy_0, g_yz_0_xxxxxxxxz_0, g_yz_0_xxxxxxxyy_0, g_yz_0_xxxxxxxyz_0, g_yz_0_xxxxxxxzz_0, g_yz_0_xxxxxxyyy_0, g_yz_0_xxxxxxyyz_0, g_yz_0_xxxxxxyzz_0, g_yz_0_xxxxxxzzz_0, g_yz_0_xxxxxyyyy_0, g_yz_0_xxxxxyyyz_0, g_yz_0_xxxxxyyzz_0, g_yz_0_xxxxxyzzz_0, g_yz_0_xxxxxzzzz_0, g_yz_0_xxxxyyyyy_0, g_yz_0_xxxxyyyyz_0, g_yz_0_xxxxyyyzz_0, g_yz_0_xxxxyyzzz_0, g_yz_0_xxxxyzzzz_0, g_yz_0_xxxxzzzzz_0, g_yz_0_xxxyyyyyy_0, g_yz_0_xxxyyyyyz_0, g_yz_0_xxxyyyyzz_0, g_yz_0_xxxyyyzzz_0, g_yz_0_xxxyyzzzz_0, g_yz_0_xxxyzzzzz_0, g_yz_0_xxxzzzzzz_0, g_yz_0_xxyyyyyyy_0, g_yz_0_xxyyyyyyz_0, g_yz_0_xxyyyyyzz_0, g_yz_0_xxyyyyzzz_0, g_yz_0_xxyyyzzzz_0, g_yz_0_xxyyzzzzz_0, g_yz_0_xxyzzzzzz_0, g_yz_0_xxzzzzzzz_0, g_yz_0_xyyyyyyyy_0, g_yz_0_xyyyyyyyz_0, g_yz_0_xyyyyyyzz_0, g_yz_0_xyyyyyzzz_0, g_yz_0_xyyyyzzzz_0, g_yz_0_xyyyzzzzz_0, g_yz_0_xyyzzzzzz_0, g_yz_0_xyzzzzzzz_0, g_yz_0_xzzzzzzzz_0, g_yz_0_yyyyyyyyy_0, g_yz_0_yyyyyyyyz_0, g_yz_0_yyyyyyyzz_0, g_yz_0_yyyyyyzzz_0, g_yz_0_yyyyyzzzz_0, g_yz_0_yyyyzzzzz_0, g_yz_0_yyyzzzzzz_0, g_yz_0_yyzzzzzzz_0, g_yz_0_yzzzzzzzz_0, g_yz_0_zzzzzzzzz_0, g_z_0_xxxxxxxxx_1, g_z_0_xxxxxxxxz_1, g_z_0_xxxxxxxyz_1, g_z_0_xxxxxxxz_1, g_z_0_xxxxxxxzz_1, g_z_0_xxxxxxyyz_1, g_z_0_xxxxxxyz_1, g_z_0_xxxxxxyzz_1, g_z_0_xxxxxxzz_1, g_z_0_xxxxxxzzz_1, g_z_0_xxxxxyyyz_1, g_z_0_xxxxxyyz_1, g_z_0_xxxxxyyzz_1, g_z_0_xxxxxyzz_1, g_z_0_xxxxxyzzz_1, g_z_0_xxxxxzzz_1, g_z_0_xxxxxzzzz_1, g_z_0_xxxxyyyyz_1, g_z_0_xxxxyyyz_1, g_z_0_xxxxyyyzz_1, g_z_0_xxxxyyzz_1, g_z_0_xxxxyyzzz_1, g_z_0_xxxxyzzz_1, g_z_0_xxxxyzzzz_1, g_z_0_xxxxzzzz_1, g_z_0_xxxxzzzzz_1, g_z_0_xxxyyyyyz_1, g_z_0_xxxyyyyz_1, g_z_0_xxxyyyyzz_1, g_z_0_xxxyyyzz_1, g_z_0_xxxyyyzzz_1, g_z_0_xxxyyzzz_1, g_z_0_xxxyyzzzz_1, g_z_0_xxxyzzzz_1, g_z_0_xxxyzzzzz_1, g_z_0_xxxzzzzz_1, g_z_0_xxxzzzzzz_1, g_z_0_xxyyyyyyz_1, g_z_0_xxyyyyyz_1, g_z_0_xxyyyyyzz_1, g_z_0_xxyyyyzz_1, g_z_0_xxyyyyzzz_1, g_z_0_xxyyyzzz_1, g_z_0_xxyyyzzzz_1, g_z_0_xxyyzzzz_1, g_z_0_xxyyzzzzz_1, g_z_0_xxyzzzzz_1, g_z_0_xxyzzzzzz_1, g_z_0_xxzzzzzz_1, g_z_0_xxzzzzzzz_1, g_z_0_xyyyyyyyz_1, g_z_0_xyyyyyyz_1, g_z_0_xyyyyyyzz_1, g_z_0_xyyyyyzz_1, g_z_0_xyyyyyzzz_1, g_z_0_xyyyyzzz_1, g_z_0_xyyyyzzzz_1, g_z_0_xyyyzzzz_1, g_z_0_xyyyzzzzz_1, g_z_0_xyyzzzzz_1, g_z_0_xyyzzzzzz_1, g_z_0_xyzzzzzz_1, g_z_0_xyzzzzzzz_1, g_z_0_xzzzzzzz_1, g_z_0_xzzzzzzzz_1, g_z_0_yyyyyyyyz_1, g_z_0_yyyyyyyz_1, g_z_0_yyyyyyyzz_1, g_z_0_yyyyyyzz_1, g_z_0_yyyyyyzzz_1, g_z_0_yyyyyzzz_1, g_z_0_yyyyyzzzz_1, g_z_0_yyyyzzzz_1, g_z_0_yyyyzzzzz_1, g_z_0_yyyzzzzz_1, g_z_0_yyyzzzzzz_1, g_z_0_yyzzzzzz_1, g_z_0_yyzzzzzzz_1, g_z_0_yzzzzzzz_1, g_z_0_yzzzzzzzz_1, g_z_0_zzzzzzzz_1, g_z_0_zzzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xxxxxxxxx_0[i] = g_z_0_xxxxxxxxx_1[i] * wa_y[i];

        g_yz_0_xxxxxxxxy_0[i] = g_y_0_xxxxxxxxy_1[i] * wa_z[i];

        g_yz_0_xxxxxxxxz_0[i] = g_z_0_xxxxxxxxz_1[i] * wa_y[i];

        g_yz_0_xxxxxxxyy_0[i] = g_y_0_xxxxxxxyy_1[i] * wa_z[i];

        g_yz_0_xxxxxxxyz_0[i] = g_z_0_xxxxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxxxyz_1[i] * wa_y[i];

        g_yz_0_xxxxxxxzz_0[i] = g_z_0_xxxxxxxzz_1[i] * wa_y[i];

        g_yz_0_xxxxxxyyy_0[i] = g_y_0_xxxxxxyyy_1[i] * wa_z[i];

        g_yz_0_xxxxxxyyz_0[i] = 2.0 * g_z_0_xxxxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxxxyyz_1[i] * wa_y[i];

        g_yz_0_xxxxxxyzz_0[i] = g_z_0_xxxxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxxxyzz_1[i] * wa_y[i];

        g_yz_0_xxxxxxzzz_0[i] = g_z_0_xxxxxxzzz_1[i] * wa_y[i];

        g_yz_0_xxxxxyyyy_0[i] = g_y_0_xxxxxyyyy_1[i] * wa_z[i];

        g_yz_0_xxxxxyyyz_0[i] = 3.0 * g_z_0_xxxxxyyz_1[i] * fi_acd_0 + g_z_0_xxxxxyyyz_1[i] * wa_y[i];

        g_yz_0_xxxxxyyzz_0[i] = 2.0 * g_z_0_xxxxxyzz_1[i] * fi_acd_0 + g_z_0_xxxxxyyzz_1[i] * wa_y[i];

        g_yz_0_xxxxxyzzz_0[i] = g_z_0_xxxxxzzz_1[i] * fi_acd_0 + g_z_0_xxxxxyzzz_1[i] * wa_y[i];

        g_yz_0_xxxxxzzzz_0[i] = g_z_0_xxxxxzzzz_1[i] * wa_y[i];

        g_yz_0_xxxxyyyyy_0[i] = g_y_0_xxxxyyyyy_1[i] * wa_z[i];

        g_yz_0_xxxxyyyyz_0[i] = 4.0 * g_z_0_xxxxyyyz_1[i] * fi_acd_0 + g_z_0_xxxxyyyyz_1[i] * wa_y[i];

        g_yz_0_xxxxyyyzz_0[i] = 3.0 * g_z_0_xxxxyyzz_1[i] * fi_acd_0 + g_z_0_xxxxyyyzz_1[i] * wa_y[i];

        g_yz_0_xxxxyyzzz_0[i] = 2.0 * g_z_0_xxxxyzzz_1[i] * fi_acd_0 + g_z_0_xxxxyyzzz_1[i] * wa_y[i];

        g_yz_0_xxxxyzzzz_0[i] = g_z_0_xxxxzzzz_1[i] * fi_acd_0 + g_z_0_xxxxyzzzz_1[i] * wa_y[i];

        g_yz_0_xxxxzzzzz_0[i] = g_z_0_xxxxzzzzz_1[i] * wa_y[i];

        g_yz_0_xxxyyyyyy_0[i] = g_y_0_xxxyyyyyy_1[i] * wa_z[i];

        g_yz_0_xxxyyyyyz_0[i] = 5.0 * g_z_0_xxxyyyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyyyz_1[i] * wa_y[i];

        g_yz_0_xxxyyyyzz_0[i] = 4.0 * g_z_0_xxxyyyzz_1[i] * fi_acd_0 + g_z_0_xxxyyyyzz_1[i] * wa_y[i];

        g_yz_0_xxxyyyzzz_0[i] = 3.0 * g_z_0_xxxyyzzz_1[i] * fi_acd_0 + g_z_0_xxxyyyzzz_1[i] * wa_y[i];

        g_yz_0_xxxyyzzzz_0[i] = 2.0 * g_z_0_xxxyzzzz_1[i] * fi_acd_0 + g_z_0_xxxyyzzzz_1[i] * wa_y[i];

        g_yz_0_xxxyzzzzz_0[i] = g_z_0_xxxzzzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzzzz_1[i] * wa_y[i];

        g_yz_0_xxxzzzzzz_0[i] = g_z_0_xxxzzzzzz_1[i] * wa_y[i];

        g_yz_0_xxyyyyyyy_0[i] = g_y_0_xxyyyyyyy_1[i] * wa_z[i];

        g_yz_0_xxyyyyyyz_0[i] = 6.0 * g_z_0_xxyyyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyyyz_1[i] * wa_y[i];

        g_yz_0_xxyyyyyzz_0[i] = 5.0 * g_z_0_xxyyyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyyyzz_1[i] * wa_y[i];

        g_yz_0_xxyyyyzzz_0[i] = 4.0 * g_z_0_xxyyyzzz_1[i] * fi_acd_0 + g_z_0_xxyyyyzzz_1[i] * wa_y[i];

        g_yz_0_xxyyyzzzz_0[i] = 3.0 * g_z_0_xxyyzzzz_1[i] * fi_acd_0 + g_z_0_xxyyyzzzz_1[i] * wa_y[i];

        g_yz_0_xxyyzzzzz_0[i] = 2.0 * g_z_0_xxyzzzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzzzz_1[i] * wa_y[i];

        g_yz_0_xxyzzzzzz_0[i] = g_z_0_xxzzzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzzzz_1[i] * wa_y[i];

        g_yz_0_xxzzzzzzz_0[i] = g_z_0_xxzzzzzzz_1[i] * wa_y[i];

        g_yz_0_xyyyyyyyy_0[i] = g_y_0_xyyyyyyyy_1[i] * wa_z[i];

        g_yz_0_xyyyyyyyz_0[i] = 7.0 * g_z_0_xyyyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyyyz_1[i] * wa_y[i];

        g_yz_0_xyyyyyyzz_0[i] = 6.0 * g_z_0_xyyyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyyyzz_1[i] * wa_y[i];

        g_yz_0_xyyyyyzzz_0[i] = 5.0 * g_z_0_xyyyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyyyzzz_1[i] * wa_y[i];

        g_yz_0_xyyyyzzzz_0[i] = 4.0 * g_z_0_xyyyzzzz_1[i] * fi_acd_0 + g_z_0_xyyyyzzzz_1[i] * wa_y[i];

        g_yz_0_xyyyzzzzz_0[i] = 3.0 * g_z_0_xyyzzzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzzzz_1[i] * wa_y[i];

        g_yz_0_xyyzzzzzz_0[i] = 2.0 * g_z_0_xyzzzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzzzz_1[i] * wa_y[i];

        g_yz_0_xyzzzzzzz_0[i] = g_z_0_xzzzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzzzz_1[i] * wa_y[i];

        g_yz_0_xzzzzzzzz_0[i] = g_z_0_xzzzzzzzz_1[i] * wa_y[i];

        g_yz_0_yyyyyyyyy_0[i] = g_y_0_yyyyyyyyy_1[i] * wa_z[i];

        g_yz_0_yyyyyyyyz_0[i] = 8.0 * g_z_0_yyyyyyyz_1[i] * fi_acd_0 + g_z_0_yyyyyyyyz_1[i] * wa_y[i];

        g_yz_0_yyyyyyyzz_0[i] = 7.0 * g_z_0_yyyyyyzz_1[i] * fi_acd_0 + g_z_0_yyyyyyyzz_1[i] * wa_y[i];

        g_yz_0_yyyyyyzzz_0[i] = 6.0 * g_z_0_yyyyyzzz_1[i] * fi_acd_0 + g_z_0_yyyyyyzzz_1[i] * wa_y[i];

        g_yz_0_yyyyyzzzz_0[i] = 5.0 * g_z_0_yyyyzzzz_1[i] * fi_acd_0 + g_z_0_yyyyyzzzz_1[i] * wa_y[i];

        g_yz_0_yyyyzzzzz_0[i] = 4.0 * g_z_0_yyyzzzzz_1[i] * fi_acd_0 + g_z_0_yyyyzzzzz_1[i] * wa_y[i];

        g_yz_0_yyyzzzzzz_0[i] = 3.0 * g_z_0_yyzzzzzz_1[i] * fi_acd_0 + g_z_0_yyyzzzzzz_1[i] * wa_y[i];

        g_yz_0_yyzzzzzzz_0[i] = 2.0 * g_z_0_yzzzzzzz_1[i] * fi_acd_0 + g_z_0_yyzzzzzzz_1[i] * wa_y[i];

        g_yz_0_yzzzzzzzz_0[i] = g_z_0_zzzzzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzzzzz_1[i] * wa_y[i];

        g_yz_0_zzzzzzzzz_0[i] = g_z_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 275-330 components of targeted buffer : DSM

    auto g_zz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm + 275);

    auto g_zz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 276);

    auto g_zz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 277);

    auto g_zz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 278);

    auto g_zz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 279);

    auto g_zz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 280);

    auto g_zz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 281);

    auto g_zz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 282);

    auto g_zz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 283);

    auto g_zz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 284);

    auto g_zz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 285);

    auto g_zz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 286);

    auto g_zz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 287);

    auto g_zz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 288);

    auto g_zz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 289);

    auto g_zz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 290);

    auto g_zz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 291);

    auto g_zz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 292);

    auto g_zz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 293);

    auto g_zz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 294);

    auto g_zz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 295);

    auto g_zz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 296);

    auto g_zz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 297);

    auto g_zz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 298);

    auto g_zz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 299);

    auto g_zz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 300);

    auto g_zz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 301);

    auto g_zz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 302);

    auto g_zz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 303);

    auto g_zz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 304);

    auto g_zz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 305);

    auto g_zz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 306);

    auto g_zz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 307);

    auto g_zz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 308);

    auto g_zz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 309);

    auto g_zz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 310);

    auto g_zz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 311);

    auto g_zz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 312);

    auto g_zz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 313);

    auto g_zz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 314);

    auto g_zz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 315);

    auto g_zz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 316);

    auto g_zz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 317);

    auto g_zz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 318);

    auto g_zz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 319);

    auto g_zz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 320);

    auto g_zz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 321);

    auto g_zz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 322);

    auto g_zz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 323);

    auto g_zz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 324);

    auto g_zz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 325);

    auto g_zz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 326);

    auto g_zz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 327);

    auto g_zz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 328);

    auto g_zz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 329);

    #pragma omp simd aligned(g_0_0_xxxxxxxxx_0, g_0_0_xxxxxxxxx_1, g_0_0_xxxxxxxxy_0, g_0_0_xxxxxxxxy_1, g_0_0_xxxxxxxxz_0, g_0_0_xxxxxxxxz_1, g_0_0_xxxxxxxyy_0, g_0_0_xxxxxxxyy_1, g_0_0_xxxxxxxyz_0, g_0_0_xxxxxxxyz_1, g_0_0_xxxxxxxzz_0, g_0_0_xxxxxxxzz_1, g_0_0_xxxxxxyyy_0, g_0_0_xxxxxxyyy_1, g_0_0_xxxxxxyyz_0, g_0_0_xxxxxxyyz_1, g_0_0_xxxxxxyzz_0, g_0_0_xxxxxxyzz_1, g_0_0_xxxxxxzzz_0, g_0_0_xxxxxxzzz_1, g_0_0_xxxxxyyyy_0, g_0_0_xxxxxyyyy_1, g_0_0_xxxxxyyyz_0, g_0_0_xxxxxyyyz_1, g_0_0_xxxxxyyzz_0, g_0_0_xxxxxyyzz_1, g_0_0_xxxxxyzzz_0, g_0_0_xxxxxyzzz_1, g_0_0_xxxxxzzzz_0, g_0_0_xxxxxzzzz_1, g_0_0_xxxxyyyyy_0, g_0_0_xxxxyyyyy_1, g_0_0_xxxxyyyyz_0, g_0_0_xxxxyyyyz_1, g_0_0_xxxxyyyzz_0, g_0_0_xxxxyyyzz_1, g_0_0_xxxxyyzzz_0, g_0_0_xxxxyyzzz_1, g_0_0_xxxxyzzzz_0, g_0_0_xxxxyzzzz_1, g_0_0_xxxxzzzzz_0, g_0_0_xxxxzzzzz_1, g_0_0_xxxyyyyyy_0, g_0_0_xxxyyyyyy_1, g_0_0_xxxyyyyyz_0, g_0_0_xxxyyyyyz_1, g_0_0_xxxyyyyzz_0, g_0_0_xxxyyyyzz_1, g_0_0_xxxyyyzzz_0, g_0_0_xxxyyyzzz_1, g_0_0_xxxyyzzzz_0, g_0_0_xxxyyzzzz_1, g_0_0_xxxyzzzzz_0, g_0_0_xxxyzzzzz_1, g_0_0_xxxzzzzzz_0, g_0_0_xxxzzzzzz_1, g_0_0_xxyyyyyyy_0, g_0_0_xxyyyyyyy_1, g_0_0_xxyyyyyyz_0, g_0_0_xxyyyyyyz_1, g_0_0_xxyyyyyzz_0, g_0_0_xxyyyyyzz_1, g_0_0_xxyyyyzzz_0, g_0_0_xxyyyyzzz_1, g_0_0_xxyyyzzzz_0, g_0_0_xxyyyzzzz_1, g_0_0_xxyyzzzzz_0, g_0_0_xxyyzzzzz_1, g_0_0_xxyzzzzzz_0, g_0_0_xxyzzzzzz_1, g_0_0_xxzzzzzzz_0, g_0_0_xxzzzzzzz_1, g_0_0_xyyyyyyyy_0, g_0_0_xyyyyyyyy_1, g_0_0_xyyyyyyyz_0, g_0_0_xyyyyyyyz_1, g_0_0_xyyyyyyzz_0, g_0_0_xyyyyyyzz_1, g_0_0_xyyyyyzzz_0, g_0_0_xyyyyyzzz_1, g_0_0_xyyyyzzzz_0, g_0_0_xyyyyzzzz_1, g_0_0_xyyyzzzzz_0, g_0_0_xyyyzzzzz_1, g_0_0_xyyzzzzzz_0, g_0_0_xyyzzzzzz_1, g_0_0_xyzzzzzzz_0, g_0_0_xyzzzzzzz_1, g_0_0_xzzzzzzzz_0, g_0_0_xzzzzzzzz_1, g_0_0_yyyyyyyyy_0, g_0_0_yyyyyyyyy_1, g_0_0_yyyyyyyyz_0, g_0_0_yyyyyyyyz_1, g_0_0_yyyyyyyzz_0, g_0_0_yyyyyyyzz_1, g_0_0_yyyyyyzzz_0, g_0_0_yyyyyyzzz_1, g_0_0_yyyyyzzzz_0, g_0_0_yyyyyzzzz_1, g_0_0_yyyyzzzzz_0, g_0_0_yyyyzzzzz_1, g_0_0_yyyzzzzzz_0, g_0_0_yyyzzzzzz_1, g_0_0_yyzzzzzzz_0, g_0_0_yyzzzzzzz_1, g_0_0_yzzzzzzzz_0, g_0_0_yzzzzzzzz_1, g_0_0_zzzzzzzzz_0, g_0_0_zzzzzzzzz_1, g_z_0_xxxxxxxx_1, g_z_0_xxxxxxxxx_1, g_z_0_xxxxxxxxy_1, g_z_0_xxxxxxxxz_1, g_z_0_xxxxxxxy_1, g_z_0_xxxxxxxyy_1, g_z_0_xxxxxxxyz_1, g_z_0_xxxxxxxz_1, g_z_0_xxxxxxxzz_1, g_z_0_xxxxxxyy_1, g_z_0_xxxxxxyyy_1, g_z_0_xxxxxxyyz_1, g_z_0_xxxxxxyz_1, g_z_0_xxxxxxyzz_1, g_z_0_xxxxxxzz_1, g_z_0_xxxxxxzzz_1, g_z_0_xxxxxyyy_1, g_z_0_xxxxxyyyy_1, g_z_0_xxxxxyyyz_1, g_z_0_xxxxxyyz_1, g_z_0_xxxxxyyzz_1, g_z_0_xxxxxyzz_1, g_z_0_xxxxxyzzz_1, g_z_0_xxxxxzzz_1, g_z_0_xxxxxzzzz_1, g_z_0_xxxxyyyy_1, g_z_0_xxxxyyyyy_1, g_z_0_xxxxyyyyz_1, g_z_0_xxxxyyyz_1, g_z_0_xxxxyyyzz_1, g_z_0_xxxxyyzz_1, g_z_0_xxxxyyzzz_1, g_z_0_xxxxyzzz_1, g_z_0_xxxxyzzzz_1, g_z_0_xxxxzzzz_1, g_z_0_xxxxzzzzz_1, g_z_0_xxxyyyyy_1, g_z_0_xxxyyyyyy_1, g_z_0_xxxyyyyyz_1, g_z_0_xxxyyyyz_1, g_z_0_xxxyyyyzz_1, g_z_0_xxxyyyzz_1, g_z_0_xxxyyyzzz_1, g_z_0_xxxyyzzz_1, g_z_0_xxxyyzzzz_1, g_z_0_xxxyzzzz_1, g_z_0_xxxyzzzzz_1, g_z_0_xxxzzzzz_1, g_z_0_xxxzzzzzz_1, g_z_0_xxyyyyyy_1, g_z_0_xxyyyyyyy_1, g_z_0_xxyyyyyyz_1, g_z_0_xxyyyyyz_1, g_z_0_xxyyyyyzz_1, g_z_0_xxyyyyzz_1, g_z_0_xxyyyyzzz_1, g_z_0_xxyyyzzz_1, g_z_0_xxyyyzzzz_1, g_z_0_xxyyzzzz_1, g_z_0_xxyyzzzzz_1, g_z_0_xxyzzzzz_1, g_z_0_xxyzzzzzz_1, g_z_0_xxzzzzzz_1, g_z_0_xxzzzzzzz_1, g_z_0_xyyyyyyy_1, g_z_0_xyyyyyyyy_1, g_z_0_xyyyyyyyz_1, g_z_0_xyyyyyyz_1, g_z_0_xyyyyyyzz_1, g_z_0_xyyyyyzz_1, g_z_0_xyyyyyzzz_1, g_z_0_xyyyyzzz_1, g_z_0_xyyyyzzzz_1, g_z_0_xyyyzzzz_1, g_z_0_xyyyzzzzz_1, g_z_0_xyyzzzzz_1, g_z_0_xyyzzzzzz_1, g_z_0_xyzzzzzz_1, g_z_0_xyzzzzzzz_1, g_z_0_xzzzzzzz_1, g_z_0_xzzzzzzzz_1, g_z_0_yyyyyyyy_1, g_z_0_yyyyyyyyy_1, g_z_0_yyyyyyyyz_1, g_z_0_yyyyyyyz_1, g_z_0_yyyyyyyzz_1, g_z_0_yyyyyyzz_1, g_z_0_yyyyyyzzz_1, g_z_0_yyyyyzzz_1, g_z_0_yyyyyzzzz_1, g_z_0_yyyyzzzz_1, g_z_0_yyyyzzzzz_1, g_z_0_yyyzzzzz_1, g_z_0_yyyzzzzzz_1, g_z_0_yyzzzzzz_1, g_z_0_yyzzzzzzz_1, g_z_0_yzzzzzzz_1, g_z_0_yzzzzzzzz_1, g_z_0_zzzzzzzz_1, g_z_0_zzzzzzzzz_1, g_zz_0_xxxxxxxxx_0, g_zz_0_xxxxxxxxy_0, g_zz_0_xxxxxxxxz_0, g_zz_0_xxxxxxxyy_0, g_zz_0_xxxxxxxyz_0, g_zz_0_xxxxxxxzz_0, g_zz_0_xxxxxxyyy_0, g_zz_0_xxxxxxyyz_0, g_zz_0_xxxxxxyzz_0, g_zz_0_xxxxxxzzz_0, g_zz_0_xxxxxyyyy_0, g_zz_0_xxxxxyyyz_0, g_zz_0_xxxxxyyzz_0, g_zz_0_xxxxxyzzz_0, g_zz_0_xxxxxzzzz_0, g_zz_0_xxxxyyyyy_0, g_zz_0_xxxxyyyyz_0, g_zz_0_xxxxyyyzz_0, g_zz_0_xxxxyyzzz_0, g_zz_0_xxxxyzzzz_0, g_zz_0_xxxxzzzzz_0, g_zz_0_xxxyyyyyy_0, g_zz_0_xxxyyyyyz_0, g_zz_0_xxxyyyyzz_0, g_zz_0_xxxyyyzzz_0, g_zz_0_xxxyyzzzz_0, g_zz_0_xxxyzzzzz_0, g_zz_0_xxxzzzzzz_0, g_zz_0_xxyyyyyyy_0, g_zz_0_xxyyyyyyz_0, g_zz_0_xxyyyyyzz_0, g_zz_0_xxyyyyzzz_0, g_zz_0_xxyyyzzzz_0, g_zz_0_xxyyzzzzz_0, g_zz_0_xxyzzzzzz_0, g_zz_0_xxzzzzzzz_0, g_zz_0_xyyyyyyyy_0, g_zz_0_xyyyyyyyz_0, g_zz_0_xyyyyyyzz_0, g_zz_0_xyyyyyzzz_0, g_zz_0_xyyyyzzzz_0, g_zz_0_xyyyzzzzz_0, g_zz_0_xyyzzzzzz_0, g_zz_0_xyzzzzzzz_0, g_zz_0_xzzzzzzzz_0, g_zz_0_yyyyyyyyy_0, g_zz_0_yyyyyyyyz_0, g_zz_0_yyyyyyyzz_0, g_zz_0_yyyyyyzzz_0, g_zz_0_yyyyyzzzz_0, g_zz_0_yyyyzzzzz_0, g_zz_0_yyyzzzzzz_0, g_zz_0_yyzzzzzzz_0, g_zz_0_yzzzzzzzz_0, g_zz_0_zzzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xxxxxxxxx_0[i] = g_0_0_xxxxxxxxx_0[i] * fbe_0 - g_0_0_xxxxxxxxx_1[i] * fz_be_0 + g_z_0_xxxxxxxxx_1[i] * wa_z[i];

        g_zz_0_xxxxxxxxy_0[i] = g_0_0_xxxxxxxxy_0[i] * fbe_0 - g_0_0_xxxxxxxxy_1[i] * fz_be_0 + g_z_0_xxxxxxxxy_1[i] * wa_z[i];

        g_zz_0_xxxxxxxxz_0[i] = g_0_0_xxxxxxxxz_0[i] * fbe_0 - g_0_0_xxxxxxxxz_1[i] * fz_be_0 + g_z_0_xxxxxxxx_1[i] * fi_acd_0 + g_z_0_xxxxxxxxz_1[i] * wa_z[i];

        g_zz_0_xxxxxxxyy_0[i] = g_0_0_xxxxxxxyy_0[i] * fbe_0 - g_0_0_xxxxxxxyy_1[i] * fz_be_0 + g_z_0_xxxxxxxyy_1[i] * wa_z[i];

        g_zz_0_xxxxxxxyz_0[i] = g_0_0_xxxxxxxyz_0[i] * fbe_0 - g_0_0_xxxxxxxyz_1[i] * fz_be_0 + g_z_0_xxxxxxxy_1[i] * fi_acd_0 + g_z_0_xxxxxxxyz_1[i] * wa_z[i];

        g_zz_0_xxxxxxxzz_0[i] = g_0_0_xxxxxxxzz_0[i] * fbe_0 - g_0_0_xxxxxxxzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxxxxz_1[i] * fi_acd_0 + g_z_0_xxxxxxxzz_1[i] * wa_z[i];

        g_zz_0_xxxxxxyyy_0[i] = g_0_0_xxxxxxyyy_0[i] * fbe_0 - g_0_0_xxxxxxyyy_1[i] * fz_be_0 + g_z_0_xxxxxxyyy_1[i] * wa_z[i];

        g_zz_0_xxxxxxyyz_0[i] = g_0_0_xxxxxxyyz_0[i] * fbe_0 - g_0_0_xxxxxxyyz_1[i] * fz_be_0 + g_z_0_xxxxxxyy_1[i] * fi_acd_0 + g_z_0_xxxxxxyyz_1[i] * wa_z[i];

        g_zz_0_xxxxxxyzz_0[i] = g_0_0_xxxxxxyzz_0[i] * fbe_0 - g_0_0_xxxxxxyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxxxyz_1[i] * fi_acd_0 + g_z_0_xxxxxxyzz_1[i] * wa_z[i];

        g_zz_0_xxxxxxzzz_0[i] = g_0_0_xxxxxxzzz_0[i] * fbe_0 - g_0_0_xxxxxxzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxxxxzz_1[i] * fi_acd_0 + g_z_0_xxxxxxzzz_1[i] * wa_z[i];

        g_zz_0_xxxxxyyyy_0[i] = g_0_0_xxxxxyyyy_0[i] * fbe_0 - g_0_0_xxxxxyyyy_1[i] * fz_be_0 + g_z_0_xxxxxyyyy_1[i] * wa_z[i];

        g_zz_0_xxxxxyyyz_0[i] = g_0_0_xxxxxyyyz_0[i] * fbe_0 - g_0_0_xxxxxyyyz_1[i] * fz_be_0 + g_z_0_xxxxxyyy_1[i] * fi_acd_0 + g_z_0_xxxxxyyyz_1[i] * wa_z[i];

        g_zz_0_xxxxxyyzz_0[i] = g_0_0_xxxxxyyzz_0[i] * fbe_0 - g_0_0_xxxxxyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxxyyz_1[i] * fi_acd_0 + g_z_0_xxxxxyyzz_1[i] * wa_z[i];

        g_zz_0_xxxxxyzzz_0[i] = g_0_0_xxxxxyzzz_0[i] * fbe_0 - g_0_0_xxxxxyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxxxyzz_1[i] * fi_acd_0 + g_z_0_xxxxxyzzz_1[i] * wa_z[i];

        g_zz_0_xxxxxzzzz_0[i] = g_0_0_xxxxxzzzz_0[i] * fbe_0 - g_0_0_xxxxxzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxxxxzzz_1[i] * fi_acd_0 + g_z_0_xxxxxzzzz_1[i] * wa_z[i];

        g_zz_0_xxxxyyyyy_0[i] = g_0_0_xxxxyyyyy_0[i] * fbe_0 - g_0_0_xxxxyyyyy_1[i] * fz_be_0 + g_z_0_xxxxyyyyy_1[i] * wa_z[i];

        g_zz_0_xxxxyyyyz_0[i] = g_0_0_xxxxyyyyz_0[i] * fbe_0 - g_0_0_xxxxyyyyz_1[i] * fz_be_0 + g_z_0_xxxxyyyy_1[i] * fi_acd_0 + g_z_0_xxxxyyyyz_1[i] * wa_z[i];

        g_zz_0_xxxxyyyzz_0[i] = g_0_0_xxxxyyyzz_0[i] * fbe_0 - g_0_0_xxxxyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxxyyyz_1[i] * fi_acd_0 + g_z_0_xxxxyyyzz_1[i] * wa_z[i];

        g_zz_0_xxxxyyzzz_0[i] = g_0_0_xxxxyyzzz_0[i] * fbe_0 - g_0_0_xxxxyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxxyyzz_1[i] * fi_acd_0 + g_z_0_xxxxyyzzz_1[i] * wa_z[i];

        g_zz_0_xxxxyzzzz_0[i] = g_0_0_xxxxyzzzz_0[i] * fbe_0 - g_0_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxxxyzzz_1[i] * fi_acd_0 + g_z_0_xxxxyzzzz_1[i] * wa_z[i];

        g_zz_0_xxxxzzzzz_0[i] = g_0_0_xxxxzzzzz_0[i] * fbe_0 - g_0_0_xxxxzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xxxxzzzz_1[i] * fi_acd_0 + g_z_0_xxxxzzzzz_1[i] * wa_z[i];

        g_zz_0_xxxyyyyyy_0[i] = g_0_0_xxxyyyyyy_0[i] * fbe_0 - g_0_0_xxxyyyyyy_1[i] * fz_be_0 + g_z_0_xxxyyyyyy_1[i] * wa_z[i];

        g_zz_0_xxxyyyyyz_0[i] = g_0_0_xxxyyyyyz_0[i] * fbe_0 - g_0_0_xxxyyyyyz_1[i] * fz_be_0 + g_z_0_xxxyyyyy_1[i] * fi_acd_0 + g_z_0_xxxyyyyyz_1[i] * wa_z[i];

        g_zz_0_xxxyyyyzz_0[i] = g_0_0_xxxyyyyzz_0[i] * fbe_0 - g_0_0_xxxyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxxyyyyz_1[i] * fi_acd_0 + g_z_0_xxxyyyyzz_1[i] * wa_z[i];

        g_zz_0_xxxyyyzzz_0[i] = g_0_0_xxxyyyzzz_0[i] * fbe_0 - g_0_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxxyyyzz_1[i] * fi_acd_0 + g_z_0_xxxyyyzzz_1[i] * wa_z[i];

        g_zz_0_xxxyyzzzz_0[i] = g_0_0_xxxyyzzzz_0[i] * fbe_0 - g_0_0_xxxyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxxyyzzz_1[i] * fi_acd_0 + g_z_0_xxxyyzzzz_1[i] * wa_z[i];

        g_zz_0_xxxyzzzzz_0[i] = g_0_0_xxxyzzzzz_0[i] * fbe_0 - g_0_0_xxxyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xxxyzzzz_1[i] * fi_acd_0 + g_z_0_xxxyzzzzz_1[i] * wa_z[i];

        g_zz_0_xxxzzzzzz_0[i] = g_0_0_xxxzzzzzz_0[i] * fbe_0 - g_0_0_xxxzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_xxxzzzzz_1[i] * fi_acd_0 + g_z_0_xxxzzzzzz_1[i] * wa_z[i];

        g_zz_0_xxyyyyyyy_0[i] = g_0_0_xxyyyyyyy_0[i] * fbe_0 - g_0_0_xxyyyyyyy_1[i] * fz_be_0 + g_z_0_xxyyyyyyy_1[i] * wa_z[i];

        g_zz_0_xxyyyyyyz_0[i] = g_0_0_xxyyyyyyz_0[i] * fbe_0 - g_0_0_xxyyyyyyz_1[i] * fz_be_0 + g_z_0_xxyyyyyy_1[i] * fi_acd_0 + g_z_0_xxyyyyyyz_1[i] * wa_z[i];

        g_zz_0_xxyyyyyzz_0[i] = g_0_0_xxyyyyyzz_0[i] * fbe_0 - g_0_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxyyyyyz_1[i] * fi_acd_0 + g_z_0_xxyyyyyzz_1[i] * wa_z[i];

        g_zz_0_xxyyyyzzz_0[i] = g_0_0_xxyyyyzzz_0[i] * fbe_0 - g_0_0_xxyyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xxyyyyzz_1[i] * fi_acd_0 + g_z_0_xxyyyyzzz_1[i] * wa_z[i];

        g_zz_0_xxyyyzzzz_0[i] = g_0_0_xxyyyzzzz_0[i] * fbe_0 - g_0_0_xxyyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xxyyyzzz_1[i] * fi_acd_0 + g_z_0_xxyyyzzzz_1[i] * wa_z[i];

        g_zz_0_xxyyzzzzz_0[i] = g_0_0_xxyyzzzzz_0[i] * fbe_0 - g_0_0_xxyyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xxyyzzzz_1[i] * fi_acd_0 + g_z_0_xxyyzzzzz_1[i] * wa_z[i];

        g_zz_0_xxyzzzzzz_0[i] = g_0_0_xxyzzzzzz_0[i] * fbe_0 - g_0_0_xxyzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_xxyzzzzz_1[i] * fi_acd_0 + g_z_0_xxyzzzzzz_1[i] * wa_z[i];

        g_zz_0_xxzzzzzzz_0[i] = g_0_0_xxzzzzzzz_0[i] * fbe_0 - g_0_0_xxzzzzzzz_1[i] * fz_be_0 + 7.0 * g_z_0_xxzzzzzz_1[i] * fi_acd_0 + g_z_0_xxzzzzzzz_1[i] * wa_z[i];

        g_zz_0_xyyyyyyyy_0[i] = g_0_0_xyyyyyyyy_0[i] * fbe_0 - g_0_0_xyyyyyyyy_1[i] * fz_be_0 + g_z_0_xyyyyyyyy_1[i] * wa_z[i];

        g_zz_0_xyyyyyyyz_0[i] = g_0_0_xyyyyyyyz_0[i] * fbe_0 - g_0_0_xyyyyyyyz_1[i] * fz_be_0 + g_z_0_xyyyyyyy_1[i] * fi_acd_0 + g_z_0_xyyyyyyyz_1[i] * wa_z[i];

        g_zz_0_xyyyyyyzz_0[i] = g_0_0_xyyyyyyzz_0[i] * fbe_0 - g_0_0_xyyyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xyyyyyyz_1[i] * fi_acd_0 + g_z_0_xyyyyyyzz_1[i] * wa_z[i];

        g_zz_0_xyyyyyzzz_0[i] = g_0_0_xyyyyyzzz_0[i] * fbe_0 - g_0_0_xyyyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xyyyyyzz_1[i] * fi_acd_0 + g_z_0_xyyyyyzzz_1[i] * wa_z[i];

        g_zz_0_xyyyyzzzz_0[i] = g_0_0_xyyyyzzzz_0[i] * fbe_0 - g_0_0_xyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_xyyyyzzz_1[i] * fi_acd_0 + g_z_0_xyyyyzzzz_1[i] * wa_z[i];

        g_zz_0_xyyyzzzzz_0[i] = g_0_0_xyyyzzzzz_0[i] * fbe_0 - g_0_0_xyyyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_xyyyzzzz_1[i] * fi_acd_0 + g_z_0_xyyyzzzzz_1[i] * wa_z[i];

        g_zz_0_xyyzzzzzz_0[i] = g_0_0_xyyzzzzzz_0[i] * fbe_0 - g_0_0_xyyzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_xyyzzzzz_1[i] * fi_acd_0 + g_z_0_xyyzzzzzz_1[i] * wa_z[i];

        g_zz_0_xyzzzzzzz_0[i] = g_0_0_xyzzzzzzz_0[i] * fbe_0 - g_0_0_xyzzzzzzz_1[i] * fz_be_0 + 7.0 * g_z_0_xyzzzzzz_1[i] * fi_acd_0 + g_z_0_xyzzzzzzz_1[i] * wa_z[i];

        g_zz_0_xzzzzzzzz_0[i] = g_0_0_xzzzzzzzz_0[i] * fbe_0 - g_0_0_xzzzzzzzz_1[i] * fz_be_0 + 8.0 * g_z_0_xzzzzzzz_1[i] * fi_acd_0 + g_z_0_xzzzzzzzz_1[i] * wa_z[i];

        g_zz_0_yyyyyyyyy_0[i] = g_0_0_yyyyyyyyy_0[i] * fbe_0 - g_0_0_yyyyyyyyy_1[i] * fz_be_0 + g_z_0_yyyyyyyyy_1[i] * wa_z[i];

        g_zz_0_yyyyyyyyz_0[i] = g_0_0_yyyyyyyyz_0[i] * fbe_0 - g_0_0_yyyyyyyyz_1[i] * fz_be_0 + g_z_0_yyyyyyyy_1[i] * fi_acd_0 + g_z_0_yyyyyyyyz_1[i] * wa_z[i];

        g_zz_0_yyyyyyyzz_0[i] = g_0_0_yyyyyyyzz_0[i] * fbe_0 - g_0_0_yyyyyyyzz_1[i] * fz_be_0 + 2.0 * g_z_0_yyyyyyyz_1[i] * fi_acd_0 + g_z_0_yyyyyyyzz_1[i] * wa_z[i];

        g_zz_0_yyyyyyzzz_0[i] = g_0_0_yyyyyyzzz_0[i] * fbe_0 - g_0_0_yyyyyyzzz_1[i] * fz_be_0 + 3.0 * g_z_0_yyyyyyzz_1[i] * fi_acd_0 + g_z_0_yyyyyyzzz_1[i] * wa_z[i];

        g_zz_0_yyyyyzzzz_0[i] = g_0_0_yyyyyzzzz_0[i] * fbe_0 - g_0_0_yyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_z_0_yyyyyzzz_1[i] * fi_acd_0 + g_z_0_yyyyyzzzz_1[i] * wa_z[i];

        g_zz_0_yyyyzzzzz_0[i] = g_0_0_yyyyzzzzz_0[i] * fbe_0 - g_0_0_yyyyzzzzz_1[i] * fz_be_0 + 5.0 * g_z_0_yyyyzzzz_1[i] * fi_acd_0 + g_z_0_yyyyzzzzz_1[i] * wa_z[i];

        g_zz_0_yyyzzzzzz_0[i] = g_0_0_yyyzzzzzz_0[i] * fbe_0 - g_0_0_yyyzzzzzz_1[i] * fz_be_0 + 6.0 * g_z_0_yyyzzzzz_1[i] * fi_acd_0 + g_z_0_yyyzzzzzz_1[i] * wa_z[i];

        g_zz_0_yyzzzzzzz_0[i] = g_0_0_yyzzzzzzz_0[i] * fbe_0 - g_0_0_yyzzzzzzz_1[i] * fz_be_0 + 7.0 * g_z_0_yyzzzzzz_1[i] * fi_acd_0 + g_z_0_yyzzzzzzz_1[i] * wa_z[i];

        g_zz_0_yzzzzzzzz_0[i] = g_0_0_yzzzzzzzz_0[i] * fbe_0 - g_0_0_yzzzzzzzz_1[i] * fz_be_0 + 8.0 * g_z_0_yzzzzzzz_1[i] * fi_acd_0 + g_z_0_yzzzzzzzz_1[i] * wa_z[i];

        g_zz_0_zzzzzzzzz_0[i] = g_0_0_zzzzzzzzz_0[i] * fbe_0 - g_0_0_zzzzzzzzz_1[i] * fz_be_0 + 9.0 * g_z_0_zzzzzzzz_1[i] * fi_acd_0 + g_z_0_zzzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

