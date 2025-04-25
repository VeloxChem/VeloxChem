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

#include "ThreeCenterElectronRepulsionPrimRecFSL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsl,
                                 size_t idx_eri_0_psl,
                                 size_t idx_eri_1_psl,
                                 size_t idx_eri_1_dsk,
                                 size_t idx_eri_1_dsl,
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

    /// Set up components of auxilary buffer : PSL

    auto g_x_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_psl);

    auto g_x_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_psl + 1);

    auto g_x_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_psl + 2);

    auto g_x_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_psl + 3);

    auto g_x_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_psl + 4);

    auto g_x_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_psl + 5);

    auto g_x_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_psl + 6);

    auto g_x_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_psl + 7);

    auto g_x_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_psl + 8);

    auto g_x_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_psl + 9);

    auto g_x_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_psl + 10);

    auto g_x_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_psl + 11);

    auto g_x_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_psl + 12);

    auto g_x_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_psl + 13);

    auto g_x_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_psl + 14);

    auto g_x_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_psl + 15);

    auto g_x_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_psl + 16);

    auto g_x_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_psl + 17);

    auto g_x_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_psl + 18);

    auto g_x_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_psl + 19);

    auto g_x_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_psl + 20);

    auto g_x_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 21);

    auto g_x_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 22);

    auto g_x_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 23);

    auto g_x_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 24);

    auto g_x_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 25);

    auto g_x_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 26);

    auto g_x_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 27);

    auto g_x_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 28);

    auto g_x_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 29);

    auto g_x_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 30);

    auto g_x_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 31);

    auto g_x_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 32);

    auto g_x_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 33);

    auto g_x_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 34);

    auto g_x_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 35);

    auto g_x_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 36);

    auto g_x_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 37);

    auto g_x_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 38);

    auto g_x_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 39);

    auto g_x_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 40);

    auto g_x_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 41);

    auto g_x_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 42);

    auto g_x_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 43);

    auto g_x_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 44);

    auto g_y_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_psl + 45);

    auto g_y_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_psl + 46);

    auto g_y_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_psl + 47);

    auto g_y_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_psl + 48);

    auto g_y_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_psl + 49);

    auto g_y_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_psl + 50);

    auto g_y_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_psl + 51);

    auto g_y_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_psl + 52);

    auto g_y_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_psl + 53);

    auto g_y_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_psl + 54);

    auto g_y_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_psl + 55);

    auto g_y_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_psl + 56);

    auto g_y_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_psl + 57);

    auto g_y_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_psl + 58);

    auto g_y_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_psl + 59);

    auto g_y_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_psl + 60);

    auto g_y_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_psl + 61);

    auto g_y_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_psl + 62);

    auto g_y_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_psl + 63);

    auto g_y_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_psl + 64);

    auto g_y_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_psl + 65);

    auto g_y_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 66);

    auto g_y_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 67);

    auto g_y_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 68);

    auto g_y_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 69);

    auto g_y_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 70);

    auto g_y_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 71);

    auto g_y_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 72);

    auto g_y_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 73);

    auto g_y_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 74);

    auto g_y_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 75);

    auto g_y_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 76);

    auto g_y_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 77);

    auto g_y_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 78);

    auto g_y_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 79);

    auto g_y_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 80);

    auto g_y_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 81);

    auto g_y_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 82);

    auto g_y_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 83);

    auto g_y_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 84);

    auto g_y_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 85);

    auto g_y_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 86);

    auto g_y_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 87);

    auto g_y_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 88);

    auto g_y_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 89);

    auto g_z_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_psl + 90);

    auto g_z_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_psl + 91);

    auto g_z_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_psl + 92);

    auto g_z_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_psl + 93);

    auto g_z_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_psl + 94);

    auto g_z_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_psl + 95);

    auto g_z_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_psl + 96);

    auto g_z_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_psl + 97);

    auto g_z_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_psl + 98);

    auto g_z_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_psl + 99);

    auto g_z_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_psl + 100);

    auto g_z_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_psl + 101);

    auto g_z_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_psl + 102);

    auto g_z_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_psl + 103);

    auto g_z_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_psl + 104);

    auto g_z_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_psl + 105);

    auto g_z_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_psl + 106);

    auto g_z_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_psl + 107);

    auto g_z_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_psl + 108);

    auto g_z_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_psl + 109);

    auto g_z_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_psl + 110);

    auto g_z_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 111);

    auto g_z_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 112);

    auto g_z_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 113);

    auto g_z_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 114);

    auto g_z_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 115);

    auto g_z_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 116);

    auto g_z_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 117);

    auto g_z_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 118);

    auto g_z_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 119);

    auto g_z_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 120);

    auto g_z_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 121);

    auto g_z_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 122);

    auto g_z_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 123);

    auto g_z_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 124);

    auto g_z_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 125);

    auto g_z_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 126);

    auto g_z_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 127);

    auto g_z_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 128);

    auto g_z_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 129);

    auto g_z_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 130);

    auto g_z_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 131);

    auto g_z_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 132);

    auto g_z_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 133);

    auto g_z_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 134);

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

    auto g_yz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 173);

    auto g_yz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 174);

    auto g_yz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 175);

    auto g_yz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 176);

    auto g_yz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 177);

    auto g_yz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 178);

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

    /// Set up components of auxilary buffer : DSL

    auto g_xx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_dsl);

    auto g_xx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_dsl + 1);

    auto g_xx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_dsl + 2);

    auto g_xx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_dsl + 3);

    auto g_xx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_dsl + 4);

    auto g_xx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_dsl + 5);

    auto g_xx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_dsl + 6);

    auto g_xx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_dsl + 7);

    auto g_xx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_dsl + 8);

    auto g_xx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_dsl + 9);

    auto g_xx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_dsl + 10);

    auto g_xx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_dsl + 11);

    auto g_xx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_dsl + 12);

    auto g_xx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_dsl + 13);

    auto g_xx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_dsl + 14);

    auto g_xx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 15);

    auto g_xx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 16);

    auto g_xx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 17);

    auto g_xx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 18);

    auto g_xx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 19);

    auto g_xx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 20);

    auto g_xx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 21);

    auto g_xx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 22);

    auto g_xx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 23);

    auto g_xx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 24);

    auto g_xx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 25);

    auto g_xx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 26);

    auto g_xx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 27);

    auto g_xx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 28);

    auto g_xx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 29);

    auto g_xx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 30);

    auto g_xx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 31);

    auto g_xx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 32);

    auto g_xx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 33);

    auto g_xx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 34);

    auto g_xx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 35);

    auto g_xx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 36);

    auto g_xx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 37);

    auto g_xx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 38);

    auto g_xx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 39);

    auto g_xx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 40);

    auto g_xx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 41);

    auto g_xx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 42);

    auto g_xx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 43);

    auto g_xx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 44);

    auto g_xy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_dsl + 46);

    auto g_xy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_dsl + 48);

    auto g_xy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_dsl + 51);

    auto g_xy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_dsl + 55);

    auto g_xy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 60);

    auto g_xy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 66);

    auto g_xy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 73);

    auto g_xz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_dsl + 90);

    auto g_xz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_dsl + 92);

    auto g_xz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_dsl + 95);

    auto g_xz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_dsl + 99);

    auto g_xz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_dsl + 104);

    auto g_xz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 110);

    auto g_xz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 117);

    auto g_xz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 125);

    auto g_yy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_dsl + 135);

    auto g_yy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_dsl + 136);

    auto g_yy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_dsl + 137);

    auto g_yy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_dsl + 138);

    auto g_yy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_dsl + 139);

    auto g_yy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_dsl + 140);

    auto g_yy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_dsl + 141);

    auto g_yy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_dsl + 142);

    auto g_yy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_dsl + 143);

    auto g_yy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_dsl + 144);

    auto g_yy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_dsl + 145);

    auto g_yy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_dsl + 146);

    auto g_yy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_dsl + 147);

    auto g_yy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_dsl + 148);

    auto g_yy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_dsl + 149);

    auto g_yy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 150);

    auto g_yy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 151);

    auto g_yy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 152);

    auto g_yy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 153);

    auto g_yy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 154);

    auto g_yy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 155);

    auto g_yy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 156);

    auto g_yy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 157);

    auto g_yy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 158);

    auto g_yy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 159);

    auto g_yy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 160);

    auto g_yy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 161);

    auto g_yy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 162);

    auto g_yy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 163);

    auto g_yy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 164);

    auto g_yy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 165);

    auto g_yy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 166);

    auto g_yy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 167);

    auto g_yy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 168);

    auto g_yy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 169);

    auto g_yy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 170);

    auto g_yy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 171);

    auto g_yy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 172);

    auto g_yy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 173);

    auto g_yy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 174);

    auto g_yy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 175);

    auto g_yy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 176);

    auto g_yy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 177);

    auto g_yy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 178);

    auto g_yy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 179);

    auto g_yz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_dsl + 184);

    auto g_yz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_dsl + 187);

    auto g_yz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_dsl + 188);

    auto g_yz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_dsl + 191);

    auto g_yz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_dsl + 192);

    auto g_yz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_dsl + 193);

    auto g_yz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 196);

    auto g_yz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 197);

    auto g_yz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 198);

    auto g_yz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 199);

    auto g_yz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 202);

    auto g_yz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 203);

    auto g_yz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 204);

    auto g_yz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 205);

    auto g_yz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 206);

    auto g_yz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 209);

    auto g_yz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 210);

    auto g_yz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 211);

    auto g_yz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 212);

    auto g_yz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 213);

    auto g_yz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 214);

    auto g_yz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 216);

    auto g_yz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 217);

    auto g_yz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 218);

    auto g_yz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 219);

    auto g_yz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 220);

    auto g_yz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 221);

    auto g_yz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 222);

    auto g_yz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 223);

    auto g_yz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 224);

    auto g_zz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_dsl + 225);

    auto g_zz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_dsl + 226);

    auto g_zz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_dsl + 227);

    auto g_zz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_dsl + 228);

    auto g_zz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_dsl + 229);

    auto g_zz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_dsl + 230);

    auto g_zz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_dsl + 231);

    auto g_zz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_dsl + 232);

    auto g_zz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_dsl + 233);

    auto g_zz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_dsl + 234);

    auto g_zz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_dsl + 235);

    auto g_zz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_dsl + 236);

    auto g_zz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_dsl + 237);

    auto g_zz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_dsl + 238);

    auto g_zz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_dsl + 239);

    auto g_zz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 240);

    auto g_zz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 241);

    auto g_zz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 242);

    auto g_zz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 243);

    auto g_zz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 244);

    auto g_zz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 245);

    auto g_zz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 246);

    auto g_zz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 247);

    auto g_zz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 248);

    auto g_zz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 249);

    auto g_zz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 250);

    auto g_zz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 251);

    auto g_zz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 252);

    auto g_zz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 253);

    auto g_zz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 254);

    auto g_zz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 255);

    auto g_zz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 256);

    auto g_zz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 257);

    auto g_zz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 258);

    auto g_zz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 259);

    auto g_zz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 260);

    auto g_zz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 261);

    auto g_zz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 262);

    auto g_zz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 263);

    auto g_zz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 264);

    auto g_zz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 265);

    auto g_zz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 266);

    auto g_zz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 267);

    auto g_zz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 268);

    auto g_zz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 269);

    /// Set up 0-45 components of targeted buffer : FSL

    auto g_xxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl);

    auto g_xxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 1);

    auto g_xxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 2);

    auto g_xxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 3);

    auto g_xxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 4);

    auto g_xxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 5);

    auto g_xxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 6);

    auto g_xxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 7);

    auto g_xxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 8);

    auto g_xxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 9);

    auto g_xxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 10);

    auto g_xxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 11);

    auto g_xxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 12);

    auto g_xxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 13);

    auto g_xxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 14);

    auto g_xxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 15);

    auto g_xxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 16);

    auto g_xxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 17);

    auto g_xxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 18);

    auto g_xxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 19);

    auto g_xxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 20);

    auto g_xxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 21);

    auto g_xxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 22);

    auto g_xxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 23);

    auto g_xxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 24);

    auto g_xxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 25);

    auto g_xxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 26);

    auto g_xxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 27);

    auto g_xxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 28);

    auto g_xxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 29);

    auto g_xxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 30);

    auto g_xxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 31);

    auto g_xxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 32);

    auto g_xxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 33);

    auto g_xxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 34);

    auto g_xxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 35);

    auto g_xxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 36);

    auto g_xxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 37);

    auto g_xxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 38);

    auto g_xxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 39);

    auto g_xxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 40);

    auto g_xxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 41);

    auto g_xxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 42);

    auto g_xxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 43);

    auto g_xxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 44);

    #pragma omp simd aligned(g_x_0_xxxxxxxx_0, g_x_0_xxxxxxxx_1, g_x_0_xxxxxxxy_0, g_x_0_xxxxxxxy_1, g_x_0_xxxxxxxz_0, g_x_0_xxxxxxxz_1, g_x_0_xxxxxxyy_0, g_x_0_xxxxxxyy_1, g_x_0_xxxxxxyz_0, g_x_0_xxxxxxyz_1, g_x_0_xxxxxxzz_0, g_x_0_xxxxxxzz_1, g_x_0_xxxxxyyy_0, g_x_0_xxxxxyyy_1, g_x_0_xxxxxyyz_0, g_x_0_xxxxxyyz_1, g_x_0_xxxxxyzz_0, g_x_0_xxxxxyzz_1, g_x_0_xxxxxzzz_0, g_x_0_xxxxxzzz_1, g_x_0_xxxxyyyy_0, g_x_0_xxxxyyyy_1, g_x_0_xxxxyyyz_0, g_x_0_xxxxyyyz_1, g_x_0_xxxxyyzz_0, g_x_0_xxxxyyzz_1, g_x_0_xxxxyzzz_0, g_x_0_xxxxyzzz_1, g_x_0_xxxxzzzz_0, g_x_0_xxxxzzzz_1, g_x_0_xxxyyyyy_0, g_x_0_xxxyyyyy_1, g_x_0_xxxyyyyz_0, g_x_0_xxxyyyyz_1, g_x_0_xxxyyyzz_0, g_x_0_xxxyyyzz_1, g_x_0_xxxyyzzz_0, g_x_0_xxxyyzzz_1, g_x_0_xxxyzzzz_0, g_x_0_xxxyzzzz_1, g_x_0_xxxzzzzz_0, g_x_0_xxxzzzzz_1, g_x_0_xxyyyyyy_0, g_x_0_xxyyyyyy_1, g_x_0_xxyyyyyz_0, g_x_0_xxyyyyyz_1, g_x_0_xxyyyyzz_0, g_x_0_xxyyyyzz_1, g_x_0_xxyyyzzz_0, g_x_0_xxyyyzzz_1, g_x_0_xxyyzzzz_0, g_x_0_xxyyzzzz_1, g_x_0_xxyzzzzz_0, g_x_0_xxyzzzzz_1, g_x_0_xxzzzzzz_0, g_x_0_xxzzzzzz_1, g_x_0_xyyyyyyy_0, g_x_0_xyyyyyyy_1, g_x_0_xyyyyyyz_0, g_x_0_xyyyyyyz_1, g_x_0_xyyyyyzz_0, g_x_0_xyyyyyzz_1, g_x_0_xyyyyzzz_0, g_x_0_xyyyyzzz_1, g_x_0_xyyyzzzz_0, g_x_0_xyyyzzzz_1, g_x_0_xyyzzzzz_0, g_x_0_xyyzzzzz_1, g_x_0_xyzzzzzz_0, g_x_0_xyzzzzzz_1, g_x_0_xzzzzzzz_0, g_x_0_xzzzzzzz_1, g_x_0_yyyyyyyy_0, g_x_0_yyyyyyyy_1, g_x_0_yyyyyyyz_0, g_x_0_yyyyyyyz_1, g_x_0_yyyyyyzz_0, g_x_0_yyyyyyzz_1, g_x_0_yyyyyzzz_0, g_x_0_yyyyyzzz_1, g_x_0_yyyyzzzz_0, g_x_0_yyyyzzzz_1, g_x_0_yyyzzzzz_0, g_x_0_yyyzzzzz_1, g_x_0_yyzzzzzz_0, g_x_0_yyzzzzzz_1, g_x_0_yzzzzzzz_0, g_x_0_yzzzzzzz_1, g_x_0_zzzzzzzz_0, g_x_0_zzzzzzzz_1, g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxxyz_1, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxyy_1, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxxyyz_1, g_xx_0_xxxxxyz_1, g_xx_0_xxxxxyzz_1, g_xx_0_xxxxxzz_1, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxyyy_1, g_xx_0_xxxxyyyy_1, g_xx_0_xxxxyyyz_1, g_xx_0_xxxxyyz_1, g_xx_0_xxxxyyzz_1, g_xx_0_xxxxyzz_1, g_xx_0_xxxxyzzz_1, g_xx_0_xxxxzzz_1, g_xx_0_xxxxzzzz_1, g_xx_0_xxxyyyy_1, g_xx_0_xxxyyyyy_1, g_xx_0_xxxyyyyz_1, g_xx_0_xxxyyyz_1, g_xx_0_xxxyyyzz_1, g_xx_0_xxxyyzz_1, g_xx_0_xxxyyzzz_1, g_xx_0_xxxyzzz_1, g_xx_0_xxxyzzzz_1, g_xx_0_xxxzzzz_1, g_xx_0_xxxzzzzz_1, g_xx_0_xxyyyyy_1, g_xx_0_xxyyyyyy_1, g_xx_0_xxyyyyyz_1, g_xx_0_xxyyyyz_1, g_xx_0_xxyyyyzz_1, g_xx_0_xxyyyzz_1, g_xx_0_xxyyyzzz_1, g_xx_0_xxyyzzz_1, g_xx_0_xxyyzzzz_1, g_xx_0_xxyzzzz_1, g_xx_0_xxyzzzzz_1, g_xx_0_xxzzzzz_1, g_xx_0_xxzzzzzz_1, g_xx_0_xyyyyyy_1, g_xx_0_xyyyyyyy_1, g_xx_0_xyyyyyyz_1, g_xx_0_xyyyyyz_1, g_xx_0_xyyyyyzz_1, g_xx_0_xyyyyzz_1, g_xx_0_xyyyyzzz_1, g_xx_0_xyyyzzz_1, g_xx_0_xyyyzzzz_1, g_xx_0_xyyzzzz_1, g_xx_0_xyyzzzzz_1, g_xx_0_xyzzzzz_1, g_xx_0_xyzzzzzz_1, g_xx_0_xzzzzzz_1, g_xx_0_xzzzzzzz_1, g_xx_0_yyyyyyy_1, g_xx_0_yyyyyyyy_1, g_xx_0_yyyyyyyz_1, g_xx_0_yyyyyyz_1, g_xx_0_yyyyyyzz_1, g_xx_0_yyyyyzz_1, g_xx_0_yyyyyzzz_1, g_xx_0_yyyyzzz_1, g_xx_0_yyyyzzzz_1, g_xx_0_yyyzzzz_1, g_xx_0_yyyzzzzz_1, g_xx_0_yyzzzzz_1, g_xx_0_yyzzzzzz_1, g_xx_0_yzzzzzz_1, g_xx_0_yzzzzzzz_1, g_xx_0_zzzzzzz_1, g_xx_0_zzzzzzzz_1, g_xxx_0_xxxxxxxx_0, g_xxx_0_xxxxxxxy_0, g_xxx_0_xxxxxxxz_0, g_xxx_0_xxxxxxyy_0, g_xxx_0_xxxxxxyz_0, g_xxx_0_xxxxxxzz_0, g_xxx_0_xxxxxyyy_0, g_xxx_0_xxxxxyyz_0, g_xxx_0_xxxxxyzz_0, g_xxx_0_xxxxxzzz_0, g_xxx_0_xxxxyyyy_0, g_xxx_0_xxxxyyyz_0, g_xxx_0_xxxxyyzz_0, g_xxx_0_xxxxyzzz_0, g_xxx_0_xxxxzzzz_0, g_xxx_0_xxxyyyyy_0, g_xxx_0_xxxyyyyz_0, g_xxx_0_xxxyyyzz_0, g_xxx_0_xxxyyzzz_0, g_xxx_0_xxxyzzzz_0, g_xxx_0_xxxzzzzz_0, g_xxx_0_xxyyyyyy_0, g_xxx_0_xxyyyyyz_0, g_xxx_0_xxyyyyzz_0, g_xxx_0_xxyyyzzz_0, g_xxx_0_xxyyzzzz_0, g_xxx_0_xxyzzzzz_0, g_xxx_0_xxzzzzzz_0, g_xxx_0_xyyyyyyy_0, g_xxx_0_xyyyyyyz_0, g_xxx_0_xyyyyyzz_0, g_xxx_0_xyyyyzzz_0, g_xxx_0_xyyyzzzz_0, g_xxx_0_xyyzzzzz_0, g_xxx_0_xyzzzzzz_0, g_xxx_0_xzzzzzzz_0, g_xxx_0_yyyyyyyy_0, g_xxx_0_yyyyyyyz_0, g_xxx_0_yyyyyyzz_0, g_xxx_0_yyyyyzzz_0, g_xxx_0_yyyyzzzz_0, g_xxx_0_yyyzzzzz_0, g_xxx_0_yyzzzzzz_0, g_xxx_0_yzzzzzzz_0, g_xxx_0_zzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xxxxxxxx_0[i] = 2.0 * g_x_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxx_1[i] * fz_be_0 + 8.0 * g_xx_0_xxxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxxx_1[i] * wa_x[i];

        g_xxx_0_xxxxxxxy_0[i] = 2.0 * g_x_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xx_0_xxxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxx_0_xxxxxxxz_0[i] = 2.0 * g_x_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xx_0_xxxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxx_0_xxxxxxyy_0[i] = 2.0 * g_x_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxx_0_xxxxxxyz_0[i] = 2.0 * g_x_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxx_0_xxxxxxzz_0[i] = 2.0 * g_x_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxx_0_xxxxxyyy_0[i] = 2.0 * g_x_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxx_0_xxxxxyyz_0[i] = 2.0 * g_x_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxx_0_xxxxxyzz_0[i] = 2.0 * g_x_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxx_0_xxxxxzzz_0[i] = 2.0 * g_x_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxx_0_xxxxyyyy_0[i] = 2.0 * g_x_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxx_0_xxxxyyyz_0[i] = 2.0 * g_x_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxx_0_xxxxyyzz_0[i] = 2.0 * g_x_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyzz_1[i] * fi_acd_0 + g_xx_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxx_0_xxxxyzzz_0[i] = 2.0 * g_x_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyzzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxx_0_xxxxzzzz_0[i] = 2.0 * g_x_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxzzzz_1[i] * fi_acd_0 + g_xx_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxx_0_xxxyyyyy_0[i] = 2.0 * g_x_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxx_0_xxxyyyyz_0[i] = 2.0 * g_x_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxx_0_xxxyyyzz_0[i] = 2.0 * g_x_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxx_0_xxxyyzzz_0[i] = 2.0 * g_x_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyzzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxx_0_xxxyzzzz_0[i] = 2.0 * g_x_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyzzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxx_0_xxxzzzzz_0[i] = 2.0 * g_x_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxzzzzz_1[i] * fi_acd_0 + g_xx_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxx_0_xxyyyyyy_0[i] = 2.0 * g_x_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxx_0_xxyyyyyz_0[i] = 2.0 * g_x_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxx_0_xxyyyyzz_0[i] = 2.0 * g_x_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxx_0_xxyyyzzz_0[i] = 2.0 * g_x_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxx_0_xxyyzzzz_0[i] = 2.0 * g_x_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyzzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxx_0_xxyzzzzz_0[i] = 2.0 * g_x_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyzzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxx_0_xxzzzzzz_0[i] = 2.0 * g_x_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xzzzzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyyyyy_0[i] = 2.0 * g_x_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxx_0_xyyyyyyz_0[i] = 2.0 * g_x_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxx_0_xyyyyyzz_0[i] = 2.0 * g_x_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyzz_1[i] * fz_be_0 + g_xx_0_yyyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxx_0_xyyyyzzz_0[i] = 2.0 * g_x_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyzzz_1[i] * fz_be_0 + g_xx_0_yyyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyzzzz_0[i] = 2.0 * g_x_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyzzzz_1[i] * fz_be_0 + g_xx_0_yyyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxx_0_xyyzzzzz_0[i] = 2.0 * g_x_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyzzzzz_1[i] * fz_be_0 + g_xx_0_yyzzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxx_0_xyzzzzzz_0[i] = 2.0 * g_x_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyzzzzzz_1[i] * fz_be_0 + g_xx_0_yzzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxx_0_xzzzzzzz_0[i] = 2.0 * g_x_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xzzzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyyyyy_0[i] = 2.0 * g_x_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxx_0_yyyyyyyz_0[i] = 2.0 * g_x_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxx_0_yyyyyyzz_0[i] = 2.0 * g_x_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyzz_1[i] * fz_be_0 + g_xx_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxx_0_yyyyyzzz_0[i] = 2.0 * g_x_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyzzz_1[i] * fz_be_0 + g_xx_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyzzzz_0[i] = 2.0 * g_x_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyzzzz_1[i] * fz_be_0 + g_xx_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyzzzzz_0[i] = 2.0 * g_x_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyzzzzz_1[i] * fz_be_0 + g_xx_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxx_0_yyzzzzzz_0[i] = 2.0 * g_x_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyzzzzzz_1[i] * fz_be_0 + g_xx_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxx_0_yzzzzzzz_0[i] = 2.0 * g_x_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yzzzzzzz_1[i] * fz_be_0 + g_xx_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxx_0_zzzzzzzz_0[i] = 2.0 * g_x_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_zzzzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 45-90 components of targeted buffer : FSL

    auto g_xxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 45);

    auto g_xxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 46);

    auto g_xxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 47);

    auto g_xxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 48);

    auto g_xxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 49);

    auto g_xxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 50);

    auto g_xxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 51);

    auto g_xxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 52);

    auto g_xxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 53);

    auto g_xxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 54);

    auto g_xxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 55);

    auto g_xxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 56);

    auto g_xxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 57);

    auto g_xxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 58);

    auto g_xxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 59);

    auto g_xxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 60);

    auto g_xxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 61);

    auto g_xxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 62);

    auto g_xxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 63);

    auto g_xxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 64);

    auto g_xxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 65);

    auto g_xxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 66);

    auto g_xxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 67);

    auto g_xxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 68);

    auto g_xxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 69);

    auto g_xxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 70);

    auto g_xxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 71);

    auto g_xxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 72);

    auto g_xxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 73);

    auto g_xxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 74);

    auto g_xxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 75);

    auto g_xxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 76);

    auto g_xxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 77);

    auto g_xxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 78);

    auto g_xxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 79);

    auto g_xxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 80);

    auto g_xxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 81);

    auto g_xxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 82);

    auto g_xxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 83);

    auto g_xxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 84);

    auto g_xxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 85);

    auto g_xxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 86);

    auto g_xxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 87);

    auto g_xxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 88);

    auto g_xxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 89);

    #pragma omp simd aligned(g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxxyz_1, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxyy_1, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxxyyz_1, g_xx_0_xxxxxyz_1, g_xx_0_xxxxxyzz_1, g_xx_0_xxxxxzz_1, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxyyy_1, g_xx_0_xxxxyyyy_1, g_xx_0_xxxxyyyz_1, g_xx_0_xxxxyyz_1, g_xx_0_xxxxyyzz_1, g_xx_0_xxxxyzz_1, g_xx_0_xxxxyzzz_1, g_xx_0_xxxxzzz_1, g_xx_0_xxxxzzzz_1, g_xx_0_xxxyyyy_1, g_xx_0_xxxyyyyy_1, g_xx_0_xxxyyyyz_1, g_xx_0_xxxyyyz_1, g_xx_0_xxxyyyzz_1, g_xx_0_xxxyyzz_1, g_xx_0_xxxyyzzz_1, g_xx_0_xxxyzzz_1, g_xx_0_xxxyzzzz_1, g_xx_0_xxxzzzz_1, g_xx_0_xxxzzzzz_1, g_xx_0_xxyyyyy_1, g_xx_0_xxyyyyyy_1, g_xx_0_xxyyyyyz_1, g_xx_0_xxyyyyz_1, g_xx_0_xxyyyyzz_1, g_xx_0_xxyyyzz_1, g_xx_0_xxyyyzzz_1, g_xx_0_xxyyzzz_1, g_xx_0_xxyyzzzz_1, g_xx_0_xxyzzzz_1, g_xx_0_xxyzzzzz_1, g_xx_0_xxzzzzz_1, g_xx_0_xxzzzzzz_1, g_xx_0_xyyyyyy_1, g_xx_0_xyyyyyyy_1, g_xx_0_xyyyyyyz_1, g_xx_0_xyyyyyz_1, g_xx_0_xyyyyyzz_1, g_xx_0_xyyyyzz_1, g_xx_0_xyyyyzzz_1, g_xx_0_xyyyzzz_1, g_xx_0_xyyyzzzz_1, g_xx_0_xyyzzzz_1, g_xx_0_xyyzzzzz_1, g_xx_0_xyzzzzz_1, g_xx_0_xyzzzzzz_1, g_xx_0_xzzzzzz_1, g_xx_0_xzzzzzzz_1, g_xx_0_yyyyyyy_1, g_xx_0_yyyyyyyy_1, g_xx_0_yyyyyyyz_1, g_xx_0_yyyyyyz_1, g_xx_0_yyyyyyzz_1, g_xx_0_yyyyyzz_1, g_xx_0_yyyyyzzz_1, g_xx_0_yyyyzzz_1, g_xx_0_yyyyzzzz_1, g_xx_0_yyyzzzz_1, g_xx_0_yyyzzzzz_1, g_xx_0_yyzzzzz_1, g_xx_0_yyzzzzzz_1, g_xx_0_yzzzzzz_1, g_xx_0_yzzzzzzz_1, g_xx_0_zzzzzzz_1, g_xx_0_zzzzzzzz_1, g_xxy_0_xxxxxxxx_0, g_xxy_0_xxxxxxxy_0, g_xxy_0_xxxxxxxz_0, g_xxy_0_xxxxxxyy_0, g_xxy_0_xxxxxxyz_0, g_xxy_0_xxxxxxzz_0, g_xxy_0_xxxxxyyy_0, g_xxy_0_xxxxxyyz_0, g_xxy_0_xxxxxyzz_0, g_xxy_0_xxxxxzzz_0, g_xxy_0_xxxxyyyy_0, g_xxy_0_xxxxyyyz_0, g_xxy_0_xxxxyyzz_0, g_xxy_0_xxxxyzzz_0, g_xxy_0_xxxxzzzz_0, g_xxy_0_xxxyyyyy_0, g_xxy_0_xxxyyyyz_0, g_xxy_0_xxxyyyzz_0, g_xxy_0_xxxyyzzz_0, g_xxy_0_xxxyzzzz_0, g_xxy_0_xxxzzzzz_0, g_xxy_0_xxyyyyyy_0, g_xxy_0_xxyyyyyz_0, g_xxy_0_xxyyyyzz_0, g_xxy_0_xxyyyzzz_0, g_xxy_0_xxyyzzzz_0, g_xxy_0_xxyzzzzz_0, g_xxy_0_xxzzzzzz_0, g_xxy_0_xyyyyyyy_0, g_xxy_0_xyyyyyyz_0, g_xxy_0_xyyyyyzz_0, g_xxy_0_xyyyyzzz_0, g_xxy_0_xyyyzzzz_0, g_xxy_0_xyyzzzzz_0, g_xxy_0_xyzzzzzz_0, g_xxy_0_xzzzzzzz_0, g_xxy_0_yyyyyyyy_0, g_xxy_0_yyyyyyyz_0, g_xxy_0_yyyyyyzz_0, g_xxy_0_yyyyyzzz_0, g_xxy_0_yyyyzzzz_0, g_xxy_0_yyyzzzzz_0, g_xxy_0_yyzzzzzz_0, g_xxy_0_yzzzzzzz_0, g_xxy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xxxxxxxx_0[i] = g_xx_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxy_0_xxxxxxxy_0[i] = g_xx_0_xxxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxy_0_xxxxxxxz_0[i] = g_xx_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxy_0_xxxxxxyy_0[i] = 2.0 * g_xx_0_xxxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxy_0_xxxxxxyz_0[i] = g_xx_0_xxxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxy_0_xxxxxxzz_0[i] = g_xx_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxy_0_xxxxxyyy_0[i] = 3.0 * g_xx_0_xxxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxy_0_xxxxxyyz_0[i] = 2.0 * g_xx_0_xxxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxy_0_xxxxxyzz_0[i] = g_xx_0_xxxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxy_0_xxxxxzzz_0[i] = g_xx_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxy_0_xxxxyyyy_0[i] = 4.0 * g_xx_0_xxxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxy_0_xxxxyyyz_0[i] = 3.0 * g_xx_0_xxxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxy_0_xxxxyyzz_0[i] = 2.0 * g_xx_0_xxxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxy_0_xxxxyzzz_0[i] = g_xx_0_xxxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxy_0_xxxxzzzz_0[i] = g_xx_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxy_0_xxxyyyyy_0[i] = 5.0 * g_xx_0_xxxyyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxy_0_xxxyyyyz_0[i] = 4.0 * g_xx_0_xxxyyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxy_0_xxxyyyzz_0[i] = 3.0 * g_xx_0_xxxyyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxy_0_xxxyyzzz_0[i] = 2.0 * g_xx_0_xxxyzzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxy_0_xxxyzzzz_0[i] = g_xx_0_xxxzzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxy_0_xxxzzzzz_0[i] = g_xx_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxy_0_xxyyyyyy_0[i] = 6.0 * g_xx_0_xxyyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxy_0_xxyyyyyz_0[i] = 5.0 * g_xx_0_xxyyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxy_0_xxyyyyzz_0[i] = 4.0 * g_xx_0_xxyyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxy_0_xxyyyzzz_0[i] = 3.0 * g_xx_0_xxyyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxy_0_xxyyzzzz_0[i] = 2.0 * g_xx_0_xxyzzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxy_0_xxyzzzzz_0[i] = g_xx_0_xxzzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxy_0_xxzzzzzz_0[i] = g_xx_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyyyyy_0[i] = 7.0 * g_xx_0_xyyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxy_0_xyyyyyyz_0[i] = 6.0 * g_xx_0_xyyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxy_0_xyyyyyzz_0[i] = 5.0 * g_xx_0_xyyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxy_0_xyyyyzzz_0[i] = 4.0 * g_xx_0_xyyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyzzzz_0[i] = 3.0 * g_xx_0_xyyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxy_0_xyyzzzzz_0[i] = 2.0 * g_xx_0_xyzzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxy_0_xyzzzzzz_0[i] = g_xx_0_xzzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxy_0_xzzzzzzz_0[i] = g_xx_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyyyyy_0[i] = 8.0 * g_xx_0_yyyyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxy_0_yyyyyyyz_0[i] = 7.0 * g_xx_0_yyyyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxy_0_yyyyyyzz_0[i] = 6.0 * g_xx_0_yyyyyzz_1[i] * fi_acd_0 + g_xx_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxy_0_yyyyyzzz_0[i] = 5.0 * g_xx_0_yyyyzzz_1[i] * fi_acd_0 + g_xx_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyzzzz_0[i] = 4.0 * g_xx_0_yyyzzzz_1[i] * fi_acd_0 + g_xx_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyzzzzz_0[i] = 3.0 * g_xx_0_yyzzzzz_1[i] * fi_acd_0 + g_xx_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxy_0_yyzzzzzz_0[i] = 2.0 * g_xx_0_yzzzzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxy_0_yzzzzzzz_0[i] = g_xx_0_zzzzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxy_0_zzzzzzzz_0[i] = g_xx_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 90-135 components of targeted buffer : FSL

    auto g_xxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 90);

    auto g_xxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 91);

    auto g_xxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 92);

    auto g_xxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 93);

    auto g_xxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 94);

    auto g_xxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 95);

    auto g_xxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 96);

    auto g_xxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 97);

    auto g_xxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 98);

    auto g_xxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 99);

    auto g_xxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 100);

    auto g_xxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 101);

    auto g_xxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 102);

    auto g_xxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 103);

    auto g_xxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 104);

    auto g_xxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 105);

    auto g_xxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 106);

    auto g_xxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 107);

    auto g_xxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 108);

    auto g_xxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 109);

    auto g_xxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 110);

    auto g_xxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 111);

    auto g_xxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 112);

    auto g_xxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 113);

    auto g_xxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 114);

    auto g_xxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 115);

    auto g_xxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 116);

    auto g_xxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 117);

    auto g_xxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 118);

    auto g_xxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 119);

    auto g_xxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 120);

    auto g_xxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 121);

    auto g_xxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 122);

    auto g_xxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 123);

    auto g_xxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 124);

    auto g_xxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 125);

    auto g_xxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 126);

    auto g_xxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 127);

    auto g_xxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 128);

    auto g_xxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 129);

    auto g_xxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 130);

    auto g_xxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 131);

    auto g_xxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 132);

    auto g_xxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 133);

    auto g_xxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 134);

    #pragma omp simd aligned(g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxxyz_1, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxyy_1, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxxyyz_1, g_xx_0_xxxxxyz_1, g_xx_0_xxxxxyzz_1, g_xx_0_xxxxxzz_1, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxyyy_1, g_xx_0_xxxxyyyy_1, g_xx_0_xxxxyyyz_1, g_xx_0_xxxxyyz_1, g_xx_0_xxxxyyzz_1, g_xx_0_xxxxyzz_1, g_xx_0_xxxxyzzz_1, g_xx_0_xxxxzzz_1, g_xx_0_xxxxzzzz_1, g_xx_0_xxxyyyy_1, g_xx_0_xxxyyyyy_1, g_xx_0_xxxyyyyz_1, g_xx_0_xxxyyyz_1, g_xx_0_xxxyyyzz_1, g_xx_0_xxxyyzz_1, g_xx_0_xxxyyzzz_1, g_xx_0_xxxyzzz_1, g_xx_0_xxxyzzzz_1, g_xx_0_xxxzzzz_1, g_xx_0_xxxzzzzz_1, g_xx_0_xxyyyyy_1, g_xx_0_xxyyyyyy_1, g_xx_0_xxyyyyyz_1, g_xx_0_xxyyyyz_1, g_xx_0_xxyyyyzz_1, g_xx_0_xxyyyzz_1, g_xx_0_xxyyyzzz_1, g_xx_0_xxyyzzz_1, g_xx_0_xxyyzzzz_1, g_xx_0_xxyzzzz_1, g_xx_0_xxyzzzzz_1, g_xx_0_xxzzzzz_1, g_xx_0_xxzzzzzz_1, g_xx_0_xyyyyyy_1, g_xx_0_xyyyyyyy_1, g_xx_0_xyyyyyyz_1, g_xx_0_xyyyyyz_1, g_xx_0_xyyyyyzz_1, g_xx_0_xyyyyzz_1, g_xx_0_xyyyyzzz_1, g_xx_0_xyyyzzz_1, g_xx_0_xyyyzzzz_1, g_xx_0_xyyzzzz_1, g_xx_0_xyyzzzzz_1, g_xx_0_xyzzzzz_1, g_xx_0_xyzzzzzz_1, g_xx_0_xzzzzzz_1, g_xx_0_xzzzzzzz_1, g_xx_0_yyyyyyy_1, g_xx_0_yyyyyyyy_1, g_xx_0_yyyyyyyz_1, g_xx_0_yyyyyyz_1, g_xx_0_yyyyyyzz_1, g_xx_0_yyyyyzz_1, g_xx_0_yyyyyzzz_1, g_xx_0_yyyyzzz_1, g_xx_0_yyyyzzzz_1, g_xx_0_yyyzzzz_1, g_xx_0_yyyzzzzz_1, g_xx_0_yyzzzzz_1, g_xx_0_yyzzzzzz_1, g_xx_0_yzzzzzz_1, g_xx_0_yzzzzzzz_1, g_xx_0_zzzzzzz_1, g_xx_0_zzzzzzzz_1, g_xxz_0_xxxxxxxx_0, g_xxz_0_xxxxxxxy_0, g_xxz_0_xxxxxxxz_0, g_xxz_0_xxxxxxyy_0, g_xxz_0_xxxxxxyz_0, g_xxz_0_xxxxxxzz_0, g_xxz_0_xxxxxyyy_0, g_xxz_0_xxxxxyyz_0, g_xxz_0_xxxxxyzz_0, g_xxz_0_xxxxxzzz_0, g_xxz_0_xxxxyyyy_0, g_xxz_0_xxxxyyyz_0, g_xxz_0_xxxxyyzz_0, g_xxz_0_xxxxyzzz_0, g_xxz_0_xxxxzzzz_0, g_xxz_0_xxxyyyyy_0, g_xxz_0_xxxyyyyz_0, g_xxz_0_xxxyyyzz_0, g_xxz_0_xxxyyzzz_0, g_xxz_0_xxxyzzzz_0, g_xxz_0_xxxzzzzz_0, g_xxz_0_xxyyyyyy_0, g_xxz_0_xxyyyyyz_0, g_xxz_0_xxyyyyzz_0, g_xxz_0_xxyyyzzz_0, g_xxz_0_xxyyzzzz_0, g_xxz_0_xxyzzzzz_0, g_xxz_0_xxzzzzzz_0, g_xxz_0_xyyyyyyy_0, g_xxz_0_xyyyyyyz_0, g_xxz_0_xyyyyyzz_0, g_xxz_0_xyyyyzzz_0, g_xxz_0_xyyyzzzz_0, g_xxz_0_xyyzzzzz_0, g_xxz_0_xyzzzzzz_0, g_xxz_0_xzzzzzzz_0, g_xxz_0_yyyyyyyy_0, g_xxz_0_yyyyyyyz_0, g_xxz_0_yyyyyyzz_0, g_xxz_0_yyyyyzzz_0, g_xxz_0_yyyyzzzz_0, g_xxz_0_yyyzzzzz_0, g_xxz_0_yyzzzzzz_0, g_xxz_0_yzzzzzzz_0, g_xxz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xxxxxxxx_0[i] = g_xx_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxz_0_xxxxxxxy_0[i] = g_xx_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxz_0_xxxxxxxz_0[i] = g_xx_0_xxxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxz_0_xxxxxxyy_0[i] = g_xx_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxz_0_xxxxxxyz_0[i] = g_xx_0_xxxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxz_0_xxxxxxzz_0[i] = 2.0 * g_xx_0_xxxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxz_0_xxxxxyyy_0[i] = g_xx_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxz_0_xxxxxyyz_0[i] = g_xx_0_xxxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxz_0_xxxxxyzz_0[i] = 2.0 * g_xx_0_xxxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxz_0_xxxxxzzz_0[i] = 3.0 * g_xx_0_xxxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxz_0_xxxxyyyy_0[i] = g_xx_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxz_0_xxxxyyyz_0[i] = g_xx_0_xxxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxz_0_xxxxyyzz_0[i] = 2.0 * g_xx_0_xxxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxz_0_xxxxyzzz_0[i] = 3.0 * g_xx_0_xxxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxz_0_xxxxzzzz_0[i] = 4.0 * g_xx_0_xxxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxz_0_xxxyyyyy_0[i] = g_xx_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxz_0_xxxyyyyz_0[i] = g_xx_0_xxxyyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxz_0_xxxyyyzz_0[i] = 2.0 * g_xx_0_xxxyyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxz_0_xxxyyzzz_0[i] = 3.0 * g_xx_0_xxxyyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxz_0_xxxyzzzz_0[i] = 4.0 * g_xx_0_xxxyzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxz_0_xxxzzzzz_0[i] = 5.0 * g_xx_0_xxxzzzz_1[i] * fi_acd_0 + g_xx_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxz_0_xxyyyyyy_0[i] = g_xx_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxz_0_xxyyyyyz_0[i] = g_xx_0_xxyyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxz_0_xxyyyyzz_0[i] = 2.0 * g_xx_0_xxyyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxz_0_xxyyyzzz_0[i] = 3.0 * g_xx_0_xxyyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxz_0_xxyyzzzz_0[i] = 4.0 * g_xx_0_xxyyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxz_0_xxyzzzzz_0[i] = 5.0 * g_xx_0_xxyzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxz_0_xxzzzzzz_0[i] = 6.0 * g_xx_0_xxzzzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyyyyy_0[i] = g_xx_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxz_0_xyyyyyyz_0[i] = g_xx_0_xyyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxz_0_xyyyyyzz_0[i] = 2.0 * g_xx_0_xyyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxz_0_xyyyyzzz_0[i] = 3.0 * g_xx_0_xyyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyzzzz_0[i] = 4.0 * g_xx_0_xyyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxz_0_xyyzzzzz_0[i] = 5.0 * g_xx_0_xyyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxz_0_xyzzzzzz_0[i] = 6.0 * g_xx_0_xyzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxz_0_xzzzzzzz_0[i] = 7.0 * g_xx_0_xzzzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyyyyy_0[i] = g_xx_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxz_0_yyyyyyyz_0[i] = g_xx_0_yyyyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxz_0_yyyyyyzz_0[i] = 2.0 * g_xx_0_yyyyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxz_0_yyyyyzzz_0[i] = 3.0 * g_xx_0_yyyyyzz_1[i] * fi_acd_0 + g_xx_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyzzzz_0[i] = 4.0 * g_xx_0_yyyyzzz_1[i] * fi_acd_0 + g_xx_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyzzzzz_0[i] = 5.0 * g_xx_0_yyyzzzz_1[i] * fi_acd_0 + g_xx_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxz_0_yyzzzzzz_0[i] = 6.0 * g_xx_0_yyzzzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxz_0_yzzzzzzz_0[i] = 7.0 * g_xx_0_yzzzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxz_0_zzzzzzzz_0[i] = 8.0 * g_xx_0_zzzzzzz_1[i] * fi_acd_0 + g_xx_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 135-180 components of targeted buffer : FSL

    auto g_xyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 135);

    auto g_xyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 136);

    auto g_xyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 137);

    auto g_xyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 138);

    auto g_xyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 139);

    auto g_xyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 140);

    auto g_xyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 141);

    auto g_xyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 142);

    auto g_xyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 143);

    auto g_xyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 144);

    auto g_xyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 145);

    auto g_xyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 146);

    auto g_xyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 147);

    auto g_xyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 148);

    auto g_xyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 149);

    auto g_xyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 150);

    auto g_xyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 151);

    auto g_xyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 152);

    auto g_xyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 153);

    auto g_xyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 154);

    auto g_xyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 155);

    auto g_xyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 156);

    auto g_xyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 157);

    auto g_xyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 158);

    auto g_xyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 159);

    auto g_xyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 160);

    auto g_xyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 161);

    auto g_xyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 162);

    auto g_xyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 163);

    auto g_xyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 164);

    auto g_xyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 165);

    auto g_xyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 166);

    auto g_xyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 167);

    auto g_xyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 168);

    auto g_xyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 169);

    auto g_xyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 170);

    auto g_xyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 171);

    auto g_xyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 172);

    auto g_xyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 173);

    auto g_xyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 174);

    auto g_xyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 175);

    auto g_xyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 176);

    auto g_xyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 177);

    auto g_xyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 178);

    auto g_xyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 179);

    #pragma omp simd aligned(g_xyy_0_xxxxxxxx_0, g_xyy_0_xxxxxxxy_0, g_xyy_0_xxxxxxxz_0, g_xyy_0_xxxxxxyy_0, g_xyy_0_xxxxxxyz_0, g_xyy_0_xxxxxxzz_0, g_xyy_0_xxxxxyyy_0, g_xyy_0_xxxxxyyz_0, g_xyy_0_xxxxxyzz_0, g_xyy_0_xxxxxzzz_0, g_xyy_0_xxxxyyyy_0, g_xyy_0_xxxxyyyz_0, g_xyy_0_xxxxyyzz_0, g_xyy_0_xxxxyzzz_0, g_xyy_0_xxxxzzzz_0, g_xyy_0_xxxyyyyy_0, g_xyy_0_xxxyyyyz_0, g_xyy_0_xxxyyyzz_0, g_xyy_0_xxxyyzzz_0, g_xyy_0_xxxyzzzz_0, g_xyy_0_xxxzzzzz_0, g_xyy_0_xxyyyyyy_0, g_xyy_0_xxyyyyyz_0, g_xyy_0_xxyyyyzz_0, g_xyy_0_xxyyyzzz_0, g_xyy_0_xxyyzzzz_0, g_xyy_0_xxyzzzzz_0, g_xyy_0_xxzzzzzz_0, g_xyy_0_xyyyyyyy_0, g_xyy_0_xyyyyyyz_0, g_xyy_0_xyyyyyzz_0, g_xyy_0_xyyyyzzz_0, g_xyy_0_xyyyzzzz_0, g_xyy_0_xyyzzzzz_0, g_xyy_0_xyzzzzzz_0, g_xyy_0_xzzzzzzz_0, g_xyy_0_yyyyyyyy_0, g_xyy_0_yyyyyyyz_0, g_xyy_0_yyyyyyzz_0, g_xyy_0_yyyyyzzz_0, g_xyy_0_yyyyzzzz_0, g_xyy_0_yyyzzzzz_0, g_xyy_0_yyzzzzzz_0, g_xyy_0_yzzzzzzz_0, g_xyy_0_zzzzzzzz_0, g_yy_0_xxxxxxx_1, g_yy_0_xxxxxxxx_1, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxxz_1, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxxz_1, g_yy_0_xxxxxxzz_1, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyz_1, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxxzz_1, g_yy_0_xxxxxzzz_1, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyzz_1, g_yy_0_xxxxyzzz_1, g_yy_0_xxxxzzz_1, g_yy_0_xxxxzzzz_1, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyzz_1, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyzzz_1, g_yy_0_xxxyzzzz_1, g_yy_0_xxxzzzz_1, g_yy_0_xxxzzzzz_1, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyzz_1, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyzzz_1, g_yy_0_xxyyzzzz_1, g_yy_0_xxyzzzz_1, g_yy_0_xxyzzzzz_1, g_yy_0_xxzzzzz_1, g_yy_0_xxzzzzzz_1, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyzz_1, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyzzz_1, g_yy_0_xyyyzzzz_1, g_yy_0_xyyzzzz_1, g_yy_0_xyyzzzzz_1, g_yy_0_xyzzzzz_1, g_yy_0_xyzzzzzz_1, g_yy_0_xzzzzzz_1, g_yy_0_xzzzzzzz_1, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyzz_1, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyzzz_1, g_yy_0_yyyyzzzz_1, g_yy_0_yyyzzzz_1, g_yy_0_yyyzzzzz_1, g_yy_0_yyzzzzz_1, g_yy_0_yyzzzzzz_1, g_yy_0_yzzzzzz_1, g_yy_0_yzzzzzzz_1, g_yy_0_zzzzzzz_1, g_yy_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xxxxxxxx_0[i] = 8.0 * g_yy_0_xxxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyy_0_xxxxxxxy_0[i] = 7.0 * g_yy_0_xxxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyy_0_xxxxxxxz_0[i] = 7.0 * g_yy_0_xxxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyy_0_xxxxxxyy_0[i] = 6.0 * g_yy_0_xxxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyy_0_xxxxxxyz_0[i] = 6.0 * g_yy_0_xxxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyy_0_xxxxxxzz_0[i] = 6.0 * g_yy_0_xxxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyy_0_xxxxxyyy_0[i] = 5.0 * g_yy_0_xxxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyy_0_xxxxxyyz_0[i] = 5.0 * g_yy_0_xxxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyy_0_xxxxxyzz_0[i] = 5.0 * g_yy_0_xxxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyy_0_xxxxxzzz_0[i] = 5.0 * g_yy_0_xxxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyy_0_xxxxyyyy_0[i] = 4.0 * g_yy_0_xxxyyyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyy_0_xxxxyyyz_0[i] = 4.0 * g_yy_0_xxxyyyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyy_0_xxxxyyzz_0[i] = 4.0 * g_yy_0_xxxyyzz_1[i] * fi_acd_0 + g_yy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyy_0_xxxxyzzz_0[i] = 4.0 * g_yy_0_xxxyzzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyy_0_xxxxzzzz_0[i] = 4.0 * g_yy_0_xxxzzzz_1[i] * fi_acd_0 + g_yy_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyy_0_xxxyyyyy_0[i] = 3.0 * g_yy_0_xxyyyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyy_0_xxxyyyyz_0[i] = 3.0 * g_yy_0_xxyyyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyy_0_xxxyyyzz_0[i] = 3.0 * g_yy_0_xxyyyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyy_0_xxxyyzzz_0[i] = 3.0 * g_yy_0_xxyyzzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyy_0_xxxyzzzz_0[i] = 3.0 * g_yy_0_xxyzzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyy_0_xxxzzzzz_0[i] = 3.0 * g_yy_0_xxzzzzz_1[i] * fi_acd_0 + g_yy_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyy_0_xxyyyyyy_0[i] = 2.0 * g_yy_0_xyyyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyy_0_xxyyyyyz_0[i] = 2.0 * g_yy_0_xyyyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyy_0_xxyyyyzz_0[i] = 2.0 * g_yy_0_xyyyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyy_0_xxyyyzzz_0[i] = 2.0 * g_yy_0_xyyyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyy_0_xxyyzzzz_0[i] = 2.0 * g_yy_0_xyyzzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyy_0_xxyzzzzz_0[i] = 2.0 * g_yy_0_xyzzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyy_0_xxzzzzzz_0[i] = 2.0 * g_yy_0_xzzzzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyyyyy_0[i] = g_yy_0_yyyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyy_0_xyyyyyyz_0[i] = g_yy_0_yyyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyy_0_xyyyyyzz_0[i] = g_yy_0_yyyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyy_0_xyyyyzzz_0[i] = g_yy_0_yyyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyzzzz_0[i] = g_yy_0_yyyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyy_0_xyyzzzzz_0[i] = g_yy_0_yyzzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyy_0_xyzzzzzz_0[i] = g_yy_0_yzzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyy_0_xzzzzzzz_0[i] = g_yy_0_zzzzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyyyyy_0[i] = g_yy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyy_0_yyyyyyyz_0[i] = g_yy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyy_0_yyyyyyzz_0[i] = g_yy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyy_0_yyyyyzzz_0[i] = g_yy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyzzzz_0[i] = g_yy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyzzzzz_0[i] = g_yy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyy_0_yyzzzzzz_0[i] = g_yy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyy_0_yzzzzzzz_0[i] = g_yy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyy_0_zzzzzzzz_0[i] = g_yy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 180-225 components of targeted buffer : FSL

    auto g_xyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 180);

    auto g_xyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 181);

    auto g_xyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 182);

    auto g_xyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 183);

    auto g_xyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 184);

    auto g_xyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 185);

    auto g_xyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 186);

    auto g_xyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 187);

    auto g_xyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 188);

    auto g_xyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 189);

    auto g_xyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 190);

    auto g_xyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 191);

    auto g_xyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 192);

    auto g_xyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 193);

    auto g_xyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 194);

    auto g_xyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 195);

    auto g_xyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 196);

    auto g_xyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 197);

    auto g_xyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 198);

    auto g_xyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 199);

    auto g_xyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 200);

    auto g_xyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 201);

    auto g_xyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 202);

    auto g_xyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 203);

    auto g_xyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 204);

    auto g_xyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 205);

    auto g_xyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 206);

    auto g_xyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 207);

    auto g_xyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 208);

    auto g_xyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 209);

    auto g_xyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 210);

    auto g_xyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 211);

    auto g_xyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 212);

    auto g_xyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 213);

    auto g_xyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 214);

    auto g_xyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 215);

    auto g_xyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 216);

    auto g_xyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 217);

    auto g_xyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 218);

    auto g_xyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 219);

    auto g_xyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 220);

    auto g_xyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 221);

    auto g_xyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 222);

    auto g_xyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 223);

    auto g_xyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 224);

    #pragma omp simd aligned(g_xy_0_xxxxxxxy_1, g_xy_0_xxxxxxyy_1, g_xy_0_xxxxxyyy_1, g_xy_0_xxxxyyyy_1, g_xy_0_xxxyyyyy_1, g_xy_0_xxyyyyyy_1, g_xy_0_xyyyyyyy_1, g_xyz_0_xxxxxxxx_0, g_xyz_0_xxxxxxxy_0, g_xyz_0_xxxxxxxz_0, g_xyz_0_xxxxxxyy_0, g_xyz_0_xxxxxxyz_0, g_xyz_0_xxxxxxzz_0, g_xyz_0_xxxxxyyy_0, g_xyz_0_xxxxxyyz_0, g_xyz_0_xxxxxyzz_0, g_xyz_0_xxxxxzzz_0, g_xyz_0_xxxxyyyy_0, g_xyz_0_xxxxyyyz_0, g_xyz_0_xxxxyyzz_0, g_xyz_0_xxxxyzzz_0, g_xyz_0_xxxxzzzz_0, g_xyz_0_xxxyyyyy_0, g_xyz_0_xxxyyyyz_0, g_xyz_0_xxxyyyzz_0, g_xyz_0_xxxyyzzz_0, g_xyz_0_xxxyzzzz_0, g_xyz_0_xxxzzzzz_0, g_xyz_0_xxyyyyyy_0, g_xyz_0_xxyyyyyz_0, g_xyz_0_xxyyyyzz_0, g_xyz_0_xxyyyzzz_0, g_xyz_0_xxyyzzzz_0, g_xyz_0_xxyzzzzz_0, g_xyz_0_xxzzzzzz_0, g_xyz_0_xyyyyyyy_0, g_xyz_0_xyyyyyyz_0, g_xyz_0_xyyyyyzz_0, g_xyz_0_xyyyyzzz_0, g_xyz_0_xyyyzzzz_0, g_xyz_0_xyyzzzzz_0, g_xyz_0_xyzzzzzz_0, g_xyz_0_xzzzzzzz_0, g_xyz_0_yyyyyyyy_0, g_xyz_0_yyyyyyyz_0, g_xyz_0_yyyyyyzz_0, g_xyz_0_yyyyyzzz_0, g_xyz_0_yyyyzzzz_0, g_xyz_0_yyyzzzzz_0, g_xyz_0_yyzzzzzz_0, g_xyz_0_yzzzzzzz_0, g_xyz_0_zzzzzzzz_0, g_xz_0_xxxxxxxx_1, g_xz_0_xxxxxxxz_1, g_xz_0_xxxxxxzz_1, g_xz_0_xxxxxzzz_1, g_xz_0_xxxxzzzz_1, g_xz_0_xxxzzzzz_1, g_xz_0_xxzzzzzz_1, g_xz_0_xzzzzzzz_1, g_yz_0_xxxxxxyz_1, g_yz_0_xxxxxyyz_1, g_yz_0_xxxxxyz_1, g_yz_0_xxxxxyzz_1, g_yz_0_xxxxyyyz_1, g_yz_0_xxxxyyz_1, g_yz_0_xxxxyyzz_1, g_yz_0_xxxxyzz_1, g_yz_0_xxxxyzzz_1, g_yz_0_xxxyyyyz_1, g_yz_0_xxxyyyz_1, g_yz_0_xxxyyyzz_1, g_yz_0_xxxyyzz_1, g_yz_0_xxxyyzzz_1, g_yz_0_xxxyzzz_1, g_yz_0_xxxyzzzz_1, g_yz_0_xxyyyyyz_1, g_yz_0_xxyyyyz_1, g_yz_0_xxyyyyzz_1, g_yz_0_xxyyyzz_1, g_yz_0_xxyyyzzz_1, g_yz_0_xxyyzzz_1, g_yz_0_xxyyzzzz_1, g_yz_0_xxyzzzz_1, g_yz_0_xxyzzzzz_1, g_yz_0_xyyyyyyz_1, g_yz_0_xyyyyyz_1, g_yz_0_xyyyyyzz_1, g_yz_0_xyyyyzz_1, g_yz_0_xyyyyzzz_1, g_yz_0_xyyyzzz_1, g_yz_0_xyyyzzzz_1, g_yz_0_xyyzzzz_1, g_yz_0_xyyzzzzz_1, g_yz_0_xyzzzzz_1, g_yz_0_xyzzzzzz_1, g_yz_0_yyyyyyyy_1, g_yz_0_yyyyyyyz_1, g_yz_0_yyyyyyz_1, g_yz_0_yyyyyyzz_1, g_yz_0_yyyyyzz_1, g_yz_0_yyyyyzzz_1, g_yz_0_yyyyzzz_1, g_yz_0_yyyyzzzz_1, g_yz_0_yyyzzzz_1, g_yz_0_yyyzzzzz_1, g_yz_0_yyzzzzz_1, g_yz_0_yyzzzzzz_1, g_yz_0_yzzzzzz_1, g_yz_0_yzzzzzzz_1, g_yz_0_zzzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyz_0_xxxxxxxx_0[i] = g_xz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xyz_0_xxxxxxxy_0[i] = g_xy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xyz_0_xxxxxxxz_0[i] = g_xz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xyz_0_xxxxxxyy_0[i] = g_xy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xyz_0_xxxxxxyz_0[i] = 6.0 * g_yz_0_xxxxxyz_1[i] * fi_acd_0 + g_yz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyz_0_xxxxxxzz_0[i] = g_xz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xyz_0_xxxxxyyy_0[i] = g_xy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xyz_0_xxxxxyyz_0[i] = 5.0 * g_yz_0_xxxxyyz_1[i] * fi_acd_0 + g_yz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyz_0_xxxxxyzz_0[i] = 5.0 * g_yz_0_xxxxyzz_1[i] * fi_acd_0 + g_yz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyz_0_xxxxxzzz_0[i] = g_xz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xyz_0_xxxxyyyy_0[i] = g_xy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xyz_0_xxxxyyyz_0[i] = 4.0 * g_yz_0_xxxyyyz_1[i] * fi_acd_0 + g_yz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyz_0_xxxxyyzz_0[i] = 4.0 * g_yz_0_xxxyyzz_1[i] * fi_acd_0 + g_yz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyz_0_xxxxyzzz_0[i] = 4.0 * g_yz_0_xxxyzzz_1[i] * fi_acd_0 + g_yz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyz_0_xxxxzzzz_0[i] = g_xz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xyz_0_xxxyyyyy_0[i] = g_xy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xyz_0_xxxyyyyz_0[i] = 3.0 * g_yz_0_xxyyyyz_1[i] * fi_acd_0 + g_yz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyz_0_xxxyyyzz_0[i] = 3.0 * g_yz_0_xxyyyzz_1[i] * fi_acd_0 + g_yz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyz_0_xxxyyzzz_0[i] = 3.0 * g_yz_0_xxyyzzz_1[i] * fi_acd_0 + g_yz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyz_0_xxxyzzzz_0[i] = 3.0 * g_yz_0_xxyzzzz_1[i] * fi_acd_0 + g_yz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyz_0_xxxzzzzz_0[i] = g_xz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xyz_0_xxyyyyyy_0[i] = g_xy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xyz_0_xxyyyyyz_0[i] = 2.0 * g_yz_0_xyyyyyz_1[i] * fi_acd_0 + g_yz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyz_0_xxyyyyzz_0[i] = 2.0 * g_yz_0_xyyyyzz_1[i] * fi_acd_0 + g_yz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyz_0_xxyyyzzz_0[i] = 2.0 * g_yz_0_xyyyzzz_1[i] * fi_acd_0 + g_yz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyz_0_xxyyzzzz_0[i] = 2.0 * g_yz_0_xyyzzzz_1[i] * fi_acd_0 + g_yz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyz_0_xxyzzzzz_0[i] = 2.0 * g_yz_0_xyzzzzz_1[i] * fi_acd_0 + g_yz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyz_0_xxzzzzzz_0[i] = g_xz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xyz_0_xyyyyyyy_0[i] = g_xy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xyz_0_xyyyyyyz_0[i] = g_yz_0_yyyyyyz_1[i] * fi_acd_0 + g_yz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyz_0_xyyyyyzz_0[i] = g_yz_0_yyyyyzz_1[i] * fi_acd_0 + g_yz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyz_0_xyyyyzzz_0[i] = g_yz_0_yyyyzzz_1[i] * fi_acd_0 + g_yz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyz_0_xyyyzzzz_0[i] = g_yz_0_yyyzzzz_1[i] * fi_acd_0 + g_yz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyz_0_xyyzzzzz_0[i] = g_yz_0_yyzzzzz_1[i] * fi_acd_0 + g_yz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyz_0_xyzzzzzz_0[i] = g_yz_0_yzzzzzz_1[i] * fi_acd_0 + g_yz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyz_0_xzzzzzzz_0[i] = g_xz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xyz_0_yyyyyyyy_0[i] = g_yz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyz_0_yyyyyyyz_0[i] = g_yz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyz_0_yyyyyyzz_0[i] = g_yz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyz_0_yyyyyzzz_0[i] = g_yz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyz_0_yyyyzzzz_0[i] = g_yz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyz_0_yyyzzzzz_0[i] = g_yz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyz_0_yyzzzzzz_0[i] = g_yz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyz_0_yzzzzzzz_0[i] = g_yz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyz_0_zzzzzzzz_0[i] = g_yz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 225-270 components of targeted buffer : FSL

    auto g_xzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 225);

    auto g_xzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 226);

    auto g_xzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 227);

    auto g_xzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 228);

    auto g_xzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 229);

    auto g_xzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 230);

    auto g_xzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 231);

    auto g_xzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 232);

    auto g_xzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 233);

    auto g_xzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 234);

    auto g_xzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 235);

    auto g_xzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 236);

    auto g_xzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 237);

    auto g_xzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 238);

    auto g_xzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 239);

    auto g_xzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 240);

    auto g_xzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 241);

    auto g_xzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 242);

    auto g_xzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 243);

    auto g_xzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 244);

    auto g_xzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 245);

    auto g_xzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 246);

    auto g_xzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 247);

    auto g_xzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 248);

    auto g_xzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 249);

    auto g_xzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 250);

    auto g_xzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 251);

    auto g_xzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 252);

    auto g_xzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 253);

    auto g_xzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 254);

    auto g_xzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 255);

    auto g_xzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 256);

    auto g_xzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 257);

    auto g_xzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 258);

    auto g_xzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 259);

    auto g_xzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 260);

    auto g_xzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 261);

    auto g_xzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 262);

    auto g_xzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 263);

    auto g_xzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 264);

    auto g_xzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 265);

    auto g_xzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 266);

    auto g_xzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 267);

    auto g_xzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 268);

    auto g_xzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 269);

    #pragma omp simd aligned(g_xzz_0_xxxxxxxx_0, g_xzz_0_xxxxxxxy_0, g_xzz_0_xxxxxxxz_0, g_xzz_0_xxxxxxyy_0, g_xzz_0_xxxxxxyz_0, g_xzz_0_xxxxxxzz_0, g_xzz_0_xxxxxyyy_0, g_xzz_0_xxxxxyyz_0, g_xzz_0_xxxxxyzz_0, g_xzz_0_xxxxxzzz_0, g_xzz_0_xxxxyyyy_0, g_xzz_0_xxxxyyyz_0, g_xzz_0_xxxxyyzz_0, g_xzz_0_xxxxyzzz_0, g_xzz_0_xxxxzzzz_0, g_xzz_0_xxxyyyyy_0, g_xzz_0_xxxyyyyz_0, g_xzz_0_xxxyyyzz_0, g_xzz_0_xxxyyzzz_0, g_xzz_0_xxxyzzzz_0, g_xzz_0_xxxzzzzz_0, g_xzz_0_xxyyyyyy_0, g_xzz_0_xxyyyyyz_0, g_xzz_0_xxyyyyzz_0, g_xzz_0_xxyyyzzz_0, g_xzz_0_xxyyzzzz_0, g_xzz_0_xxyzzzzz_0, g_xzz_0_xxzzzzzz_0, g_xzz_0_xyyyyyyy_0, g_xzz_0_xyyyyyyz_0, g_xzz_0_xyyyyyzz_0, g_xzz_0_xyyyyzzz_0, g_xzz_0_xyyyzzzz_0, g_xzz_0_xyyzzzzz_0, g_xzz_0_xyzzzzzz_0, g_xzz_0_xzzzzzzz_0, g_xzz_0_yyyyyyyy_0, g_xzz_0_yyyyyyyz_0, g_xzz_0_yyyyyyzz_0, g_xzz_0_yyyyyzzz_0, g_xzz_0_yyyyzzzz_0, g_xzz_0_yyyzzzzz_0, g_xzz_0_yyzzzzzz_0, g_xzz_0_yzzzzzzz_0, g_xzz_0_zzzzzzzz_0, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxy_1, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxy_1, g_zz_0_xxxxxxyy_1, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxyy_1, g_zz_0_xxxxxyyy_1, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxzz_1, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxyyy_1, g_zz_0_xxxxyyyy_1, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyzz_1, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxzzz_1, g_zz_0_xxxxzzzz_1, g_zz_0_xxxyyyy_1, g_zz_0_xxxyyyyy_1, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyzz_1, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyzzz_1, g_zz_0_xxxyzzzz_1, g_zz_0_xxxzzzz_1, g_zz_0_xxxzzzzz_1, g_zz_0_xxyyyyy_1, g_zz_0_xxyyyyyy_1, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyzz_1, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyzzz_1, g_zz_0_xxyyzzzz_1, g_zz_0_xxyzzzz_1, g_zz_0_xxyzzzzz_1, g_zz_0_xxzzzzz_1, g_zz_0_xxzzzzzz_1, g_zz_0_xyyyyyy_1, g_zz_0_xyyyyyyy_1, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyzz_1, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyzzz_1, g_zz_0_xyyyzzzz_1, g_zz_0_xyyzzzz_1, g_zz_0_xyyzzzzz_1, g_zz_0_xyzzzzz_1, g_zz_0_xyzzzzzz_1, g_zz_0_xzzzzzz_1, g_zz_0_xzzzzzzz_1, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyzz_1, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyzzz_1, g_zz_0_yyyyzzzz_1, g_zz_0_yyyzzzz_1, g_zz_0_yyyzzzzz_1, g_zz_0_yyzzzzz_1, g_zz_0_yyzzzzzz_1, g_zz_0_yzzzzzz_1, g_zz_0_yzzzzzzz_1, g_zz_0_zzzzzzz_1, g_zz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xxxxxxxx_0[i] = 8.0 * g_zz_0_xxxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xzz_0_xxxxxxxy_0[i] = 7.0 * g_zz_0_xxxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xzz_0_xxxxxxxz_0[i] = 7.0 * g_zz_0_xxxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xzz_0_xxxxxxyy_0[i] = 6.0 * g_zz_0_xxxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xzz_0_xxxxxxyz_0[i] = 6.0 * g_zz_0_xxxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xzz_0_xxxxxxzz_0[i] = 6.0 * g_zz_0_xxxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xzz_0_xxxxxyyy_0[i] = 5.0 * g_zz_0_xxxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xzz_0_xxxxxyyz_0[i] = 5.0 * g_zz_0_xxxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xzz_0_xxxxxyzz_0[i] = 5.0 * g_zz_0_xxxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xzz_0_xxxxxzzz_0[i] = 5.0 * g_zz_0_xxxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xzz_0_xxxxyyyy_0[i] = 4.0 * g_zz_0_xxxyyyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xzz_0_xxxxyyyz_0[i] = 4.0 * g_zz_0_xxxyyyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xzz_0_xxxxyyzz_0[i] = 4.0 * g_zz_0_xxxyyzz_1[i] * fi_acd_0 + g_zz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xzz_0_xxxxyzzz_0[i] = 4.0 * g_zz_0_xxxyzzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xzz_0_xxxxzzzz_0[i] = 4.0 * g_zz_0_xxxzzzz_1[i] * fi_acd_0 + g_zz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xzz_0_xxxyyyyy_0[i] = 3.0 * g_zz_0_xxyyyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xzz_0_xxxyyyyz_0[i] = 3.0 * g_zz_0_xxyyyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xzz_0_xxxyyyzz_0[i] = 3.0 * g_zz_0_xxyyyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xzz_0_xxxyyzzz_0[i] = 3.0 * g_zz_0_xxyyzzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xzz_0_xxxyzzzz_0[i] = 3.0 * g_zz_0_xxyzzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xzz_0_xxxzzzzz_0[i] = 3.0 * g_zz_0_xxzzzzz_1[i] * fi_acd_0 + g_zz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xzz_0_xxyyyyyy_0[i] = 2.0 * g_zz_0_xyyyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xzz_0_xxyyyyyz_0[i] = 2.0 * g_zz_0_xyyyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xzz_0_xxyyyyzz_0[i] = 2.0 * g_zz_0_xyyyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xzz_0_xxyyyzzz_0[i] = 2.0 * g_zz_0_xyyyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xzz_0_xxyyzzzz_0[i] = 2.0 * g_zz_0_xyyzzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xzz_0_xxyzzzzz_0[i] = 2.0 * g_zz_0_xyzzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xzz_0_xxzzzzzz_0[i] = 2.0 * g_zz_0_xzzzzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyyyyy_0[i] = g_zz_0_yyyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xzz_0_xyyyyyyz_0[i] = g_zz_0_yyyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xzz_0_xyyyyyzz_0[i] = g_zz_0_yyyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xzz_0_xyyyyzzz_0[i] = g_zz_0_yyyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyzzzz_0[i] = g_zz_0_yyyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xzz_0_xyyzzzzz_0[i] = g_zz_0_yyzzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xzz_0_xyzzzzzz_0[i] = g_zz_0_yzzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xzz_0_xzzzzzzz_0[i] = g_zz_0_zzzzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyyyyy_0[i] = g_zz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xzz_0_yyyyyyyz_0[i] = g_zz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xzz_0_yyyyyyzz_0[i] = g_zz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xzz_0_yyyyyzzz_0[i] = g_zz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyzzzz_0[i] = g_zz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyzzzzz_0[i] = g_zz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xzz_0_yyzzzzzz_0[i] = g_zz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xzz_0_yzzzzzzz_0[i] = g_zz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xzz_0_zzzzzzzz_0[i] = g_zz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 270-315 components of targeted buffer : FSL

    auto g_yyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 270);

    auto g_yyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 271);

    auto g_yyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 272);

    auto g_yyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 273);

    auto g_yyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 274);

    auto g_yyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 275);

    auto g_yyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 276);

    auto g_yyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 277);

    auto g_yyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 278);

    auto g_yyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 279);

    auto g_yyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 280);

    auto g_yyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 281);

    auto g_yyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 282);

    auto g_yyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 283);

    auto g_yyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 284);

    auto g_yyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 285);

    auto g_yyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 286);

    auto g_yyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 287);

    auto g_yyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 288);

    auto g_yyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 289);

    auto g_yyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 290);

    auto g_yyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 291);

    auto g_yyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 292);

    auto g_yyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 293);

    auto g_yyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 294);

    auto g_yyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 295);

    auto g_yyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 296);

    auto g_yyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 297);

    auto g_yyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 298);

    auto g_yyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 299);

    auto g_yyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 300);

    auto g_yyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 301);

    auto g_yyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 302);

    auto g_yyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 303);

    auto g_yyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 304);

    auto g_yyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 305);

    auto g_yyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 306);

    auto g_yyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 307);

    auto g_yyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 308);

    auto g_yyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 309);

    auto g_yyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 310);

    auto g_yyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 311);

    auto g_yyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 312);

    auto g_yyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 313);

    auto g_yyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 314);

    #pragma omp simd aligned(g_y_0_xxxxxxxx_0, g_y_0_xxxxxxxx_1, g_y_0_xxxxxxxy_0, g_y_0_xxxxxxxy_1, g_y_0_xxxxxxxz_0, g_y_0_xxxxxxxz_1, g_y_0_xxxxxxyy_0, g_y_0_xxxxxxyy_1, g_y_0_xxxxxxyz_0, g_y_0_xxxxxxyz_1, g_y_0_xxxxxxzz_0, g_y_0_xxxxxxzz_1, g_y_0_xxxxxyyy_0, g_y_0_xxxxxyyy_1, g_y_0_xxxxxyyz_0, g_y_0_xxxxxyyz_1, g_y_0_xxxxxyzz_0, g_y_0_xxxxxyzz_1, g_y_0_xxxxxzzz_0, g_y_0_xxxxxzzz_1, g_y_0_xxxxyyyy_0, g_y_0_xxxxyyyy_1, g_y_0_xxxxyyyz_0, g_y_0_xxxxyyyz_1, g_y_0_xxxxyyzz_0, g_y_0_xxxxyyzz_1, g_y_0_xxxxyzzz_0, g_y_0_xxxxyzzz_1, g_y_0_xxxxzzzz_0, g_y_0_xxxxzzzz_1, g_y_0_xxxyyyyy_0, g_y_0_xxxyyyyy_1, g_y_0_xxxyyyyz_0, g_y_0_xxxyyyyz_1, g_y_0_xxxyyyzz_0, g_y_0_xxxyyyzz_1, g_y_0_xxxyyzzz_0, g_y_0_xxxyyzzz_1, g_y_0_xxxyzzzz_0, g_y_0_xxxyzzzz_1, g_y_0_xxxzzzzz_0, g_y_0_xxxzzzzz_1, g_y_0_xxyyyyyy_0, g_y_0_xxyyyyyy_1, g_y_0_xxyyyyyz_0, g_y_0_xxyyyyyz_1, g_y_0_xxyyyyzz_0, g_y_0_xxyyyyzz_1, g_y_0_xxyyyzzz_0, g_y_0_xxyyyzzz_1, g_y_0_xxyyzzzz_0, g_y_0_xxyyzzzz_1, g_y_0_xxyzzzzz_0, g_y_0_xxyzzzzz_1, g_y_0_xxzzzzzz_0, g_y_0_xxzzzzzz_1, g_y_0_xyyyyyyy_0, g_y_0_xyyyyyyy_1, g_y_0_xyyyyyyz_0, g_y_0_xyyyyyyz_1, g_y_0_xyyyyyzz_0, g_y_0_xyyyyyzz_1, g_y_0_xyyyyzzz_0, g_y_0_xyyyyzzz_1, g_y_0_xyyyzzzz_0, g_y_0_xyyyzzzz_1, g_y_0_xyyzzzzz_0, g_y_0_xyyzzzzz_1, g_y_0_xyzzzzzz_0, g_y_0_xyzzzzzz_1, g_y_0_xzzzzzzz_0, g_y_0_xzzzzzzz_1, g_y_0_yyyyyyyy_0, g_y_0_yyyyyyyy_1, g_y_0_yyyyyyyz_0, g_y_0_yyyyyyyz_1, g_y_0_yyyyyyzz_0, g_y_0_yyyyyyzz_1, g_y_0_yyyyyzzz_0, g_y_0_yyyyyzzz_1, g_y_0_yyyyzzzz_0, g_y_0_yyyyzzzz_1, g_y_0_yyyzzzzz_0, g_y_0_yyyzzzzz_1, g_y_0_yyzzzzzz_0, g_y_0_yyzzzzzz_1, g_y_0_yzzzzzzz_0, g_y_0_yzzzzzzz_1, g_y_0_zzzzzzzz_0, g_y_0_zzzzzzzz_1, g_yy_0_xxxxxxx_1, g_yy_0_xxxxxxxx_1, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxxz_1, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxxz_1, g_yy_0_xxxxxxzz_1, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyz_1, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxxzz_1, g_yy_0_xxxxxzzz_1, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyzz_1, g_yy_0_xxxxyzzz_1, g_yy_0_xxxxzzz_1, g_yy_0_xxxxzzzz_1, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyzz_1, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyzzz_1, g_yy_0_xxxyzzzz_1, g_yy_0_xxxzzzz_1, g_yy_0_xxxzzzzz_1, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyzz_1, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyzzz_1, g_yy_0_xxyyzzzz_1, g_yy_0_xxyzzzz_1, g_yy_0_xxyzzzzz_1, g_yy_0_xxzzzzz_1, g_yy_0_xxzzzzzz_1, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyzz_1, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyzzz_1, g_yy_0_xyyyzzzz_1, g_yy_0_xyyzzzz_1, g_yy_0_xyyzzzzz_1, g_yy_0_xyzzzzz_1, g_yy_0_xyzzzzzz_1, g_yy_0_xzzzzzz_1, g_yy_0_xzzzzzzz_1, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyzz_1, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyzzz_1, g_yy_0_yyyyzzzz_1, g_yy_0_yyyzzzz_1, g_yy_0_yyyzzzzz_1, g_yy_0_yyzzzzz_1, g_yy_0_yyzzzzzz_1, g_yy_0_yzzzzzz_1, g_yy_0_yzzzzzzz_1, g_yy_0_zzzzzzz_1, g_yy_0_zzzzzzzz_1, g_yyy_0_xxxxxxxx_0, g_yyy_0_xxxxxxxy_0, g_yyy_0_xxxxxxxz_0, g_yyy_0_xxxxxxyy_0, g_yyy_0_xxxxxxyz_0, g_yyy_0_xxxxxxzz_0, g_yyy_0_xxxxxyyy_0, g_yyy_0_xxxxxyyz_0, g_yyy_0_xxxxxyzz_0, g_yyy_0_xxxxxzzz_0, g_yyy_0_xxxxyyyy_0, g_yyy_0_xxxxyyyz_0, g_yyy_0_xxxxyyzz_0, g_yyy_0_xxxxyzzz_0, g_yyy_0_xxxxzzzz_0, g_yyy_0_xxxyyyyy_0, g_yyy_0_xxxyyyyz_0, g_yyy_0_xxxyyyzz_0, g_yyy_0_xxxyyzzz_0, g_yyy_0_xxxyzzzz_0, g_yyy_0_xxxzzzzz_0, g_yyy_0_xxyyyyyy_0, g_yyy_0_xxyyyyyz_0, g_yyy_0_xxyyyyzz_0, g_yyy_0_xxyyyzzz_0, g_yyy_0_xxyyzzzz_0, g_yyy_0_xxyzzzzz_0, g_yyy_0_xxzzzzzz_0, g_yyy_0_xyyyyyyy_0, g_yyy_0_xyyyyyyz_0, g_yyy_0_xyyyyyzz_0, g_yyy_0_xyyyyzzz_0, g_yyy_0_xyyyzzzz_0, g_yyy_0_xyyzzzzz_0, g_yyy_0_xyzzzzzz_0, g_yyy_0_xzzzzzzz_0, g_yyy_0_yyyyyyyy_0, g_yyy_0_yyyyyyyz_0, g_yyy_0_yyyyyyzz_0, g_yyy_0_yyyyyzzz_0, g_yyy_0_yyyyzzzz_0, g_yyy_0_yyyzzzzz_0, g_yyy_0_yyzzzzzz_0, g_yyy_0_yzzzzzzz_0, g_yyy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xxxxxxxx_0[i] = 2.0 * g_y_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxx_1[i] * fz_be_0 + g_yy_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyy_0_xxxxxxxy_0[i] = 2.0 * g_y_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxy_1[i] * fz_be_0 + g_yy_0_xxxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxxy_1[i] * wa_y[i];

        g_yyy_0_xxxxxxxz_0[i] = 2.0 * g_y_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxz_1[i] * fz_be_0 + g_yy_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyy_0_xxxxxxyy_0[i] = 2.0 * g_y_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxyy_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxxyy_1[i] * wa_y[i];

        g_yyy_0_xxxxxxyz_0[i] = 2.0 * g_y_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxyz_1[i] * fz_be_0 + g_yy_0_xxxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyy_0_xxxxxxzz_0[i] = 2.0 * g_y_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxzz_1[i] * fz_be_0 + g_yy_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyy_0_xxxxxyyy_0[i] = 2.0 * g_y_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyyy_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxxyyy_1[i] * wa_y[i];

        g_yyy_0_xxxxxyyz_0[i] = 2.0 * g_y_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyy_0_xxxxxyzz_0[i] = 2.0 * g_y_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyzz_1[i] * fz_be_0 + g_yy_0_xxxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyy_0_xxxxxzzz_0[i] = 2.0 * g_y_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxzzz_1[i] * fz_be_0 + g_yy_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyy_0_xxxxyyyy_0[i] = 2.0 * g_y_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_yy_0_xxxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyyy_1[i] * wa_y[i];

        g_yyy_0_xxxxyyyz_0[i] = 2.0 * g_y_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyy_0_xxxxyyzz_0[i] = 2.0 * g_y_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyy_0_xxxxyzzz_0[i] = 2.0 * g_y_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyzzz_1[i] * fz_be_0 + g_yy_0_xxxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyy_0_xxxxzzzz_0[i] = 2.0 * g_y_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxzzzz_1[i] * fz_be_0 + g_yy_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyy_0_xxxyyyyy_0[i] = 2.0 * g_y_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyyy_1[i] * fz_be_0 + 5.0 * g_yy_0_xxxyyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyyy_1[i] * wa_y[i];

        g_yyy_0_xxxyyyyz_0[i] = 2.0 * g_y_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yy_0_xxxyyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyy_0_xxxyyyzz_0[i] = 2.0 * g_y_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxyyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyy_0_xxxyyzzz_0[i] = 2.0 * g_y_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxyzzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyy_0_xxxyzzzz_0[i] = 2.0 * g_y_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyzzzz_1[i] * fz_be_0 + g_yy_0_xxxzzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyy_0_xxxzzzzz_0[i] = 2.0 * g_y_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxzzzzz_1[i] * fz_be_0 + g_yy_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyy_0_xxyyyyyy_0[i] = 2.0 * g_y_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyyy_1[i] * fz_be_0 + 6.0 * g_yy_0_xxyyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyyy_1[i] * wa_y[i];

        g_yyy_0_xxyyyyyz_0[i] = 2.0 * g_y_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yy_0_xxyyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyy_0_xxyyyyzz_0[i] = 2.0 * g_y_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yy_0_xxyyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyy_0_xxyyyzzz_0[i] = 2.0 * g_y_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxyyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyy_0_xxyyzzzz_0[i] = 2.0 * g_y_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxyzzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyy_0_xxyzzzzz_0[i] = 2.0 * g_y_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyzzzzz_1[i] * fz_be_0 + g_yy_0_xxzzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyy_0_xxzzzzzz_0[i] = 2.0 * g_y_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxzzzzzz_1[i] * fz_be_0 + g_yy_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyyyyy_0[i] = 2.0 * g_y_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyyy_1[i] * fz_be_0 + 7.0 * g_yy_0_xyyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyyy_1[i] * wa_y[i];

        g_yyy_0_xyyyyyyz_0[i] = 2.0 * g_y_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yy_0_xyyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyy_0_xyyyyyzz_0[i] = 2.0 * g_y_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yy_0_xyyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyy_0_xyyyyzzz_0[i] = 2.0 * g_y_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yy_0_xyyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyzzzz_0[i] = 2.0 * g_y_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xyyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyy_0_xyyzzzzz_0[i] = 2.0 * g_y_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xyzzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyy_0_xyzzzzzz_0[i] = 2.0 * g_y_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyzzzzzz_1[i] * fz_be_0 + g_yy_0_xzzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyy_0_xzzzzzzz_0[i] = 2.0 * g_y_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xzzzzzzz_1[i] * fz_be_0 + g_yy_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyyyyy_0[i] = 2.0 * g_y_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyyy_1[i] * fz_be_0 + 8.0 * g_yy_0_yyyyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyyyy_1[i] * wa_y[i];

        g_yyy_0_yyyyyyyz_0[i] = 2.0 * g_y_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yy_0_yyyyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyy_0_yyyyyyzz_0[i] = 2.0 * g_y_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yy_0_yyyyyzz_1[i] * fi_acd_0 + g_yy_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyy_0_yyyyyzzz_0[i] = 2.0 * g_y_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yy_0_yyyyzzz_1[i] * fi_acd_0 + g_yy_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyzzzz_0[i] = 2.0 * g_y_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yy_0_yyyzzzz_1[i] * fi_acd_0 + g_yy_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyzzzzz_0[i] = 2.0 * g_y_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_yyzzzzz_1[i] * fi_acd_0 + g_yy_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyy_0_yyzzzzzz_0[i] = 2.0 * g_y_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_yzzzzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyy_0_yzzzzzzz_0[i] = 2.0 * g_y_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yzzzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyy_0_zzzzzzzz_0[i] = 2.0 * g_y_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_zzzzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 315-360 components of targeted buffer : FSL

    auto g_yyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 315);

    auto g_yyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 316);

    auto g_yyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 317);

    auto g_yyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 318);

    auto g_yyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 319);

    auto g_yyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 320);

    auto g_yyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 321);

    auto g_yyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 322);

    auto g_yyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 323);

    auto g_yyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 324);

    auto g_yyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 325);

    auto g_yyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 326);

    auto g_yyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 327);

    auto g_yyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 328);

    auto g_yyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 329);

    auto g_yyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 330);

    auto g_yyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 331);

    auto g_yyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 332);

    auto g_yyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 333);

    auto g_yyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 334);

    auto g_yyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 335);

    auto g_yyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 336);

    auto g_yyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 337);

    auto g_yyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 338);

    auto g_yyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 339);

    auto g_yyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 340);

    auto g_yyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 341);

    auto g_yyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 342);

    auto g_yyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 343);

    auto g_yyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 344);

    auto g_yyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 345);

    auto g_yyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 346);

    auto g_yyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 347);

    auto g_yyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 348);

    auto g_yyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 349);

    auto g_yyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 350);

    auto g_yyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 351);

    auto g_yyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 352);

    auto g_yyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 353);

    auto g_yyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 354);

    auto g_yyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 355);

    auto g_yyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 356);

    auto g_yyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 357);

    auto g_yyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 358);

    auto g_yyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 359);

    #pragma omp simd aligned(g_yy_0_xxxxxxx_1, g_yy_0_xxxxxxxx_1, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxxz_1, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxxz_1, g_yy_0_xxxxxxzz_1, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyz_1, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxxzz_1, g_yy_0_xxxxxzzz_1, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyzz_1, g_yy_0_xxxxyzzz_1, g_yy_0_xxxxzzz_1, g_yy_0_xxxxzzzz_1, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyzz_1, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyzzz_1, g_yy_0_xxxyzzzz_1, g_yy_0_xxxzzzz_1, g_yy_0_xxxzzzzz_1, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyzz_1, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyzzz_1, g_yy_0_xxyyzzzz_1, g_yy_0_xxyzzzz_1, g_yy_0_xxyzzzzz_1, g_yy_0_xxzzzzz_1, g_yy_0_xxzzzzzz_1, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyzz_1, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyzzz_1, g_yy_0_xyyyzzzz_1, g_yy_0_xyyzzzz_1, g_yy_0_xyyzzzzz_1, g_yy_0_xyzzzzz_1, g_yy_0_xyzzzzzz_1, g_yy_0_xzzzzzz_1, g_yy_0_xzzzzzzz_1, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyzz_1, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyzzz_1, g_yy_0_yyyyzzzz_1, g_yy_0_yyyzzzz_1, g_yy_0_yyyzzzzz_1, g_yy_0_yyzzzzz_1, g_yy_0_yyzzzzzz_1, g_yy_0_yzzzzzz_1, g_yy_0_yzzzzzzz_1, g_yy_0_zzzzzzz_1, g_yy_0_zzzzzzzz_1, g_yyz_0_xxxxxxxx_0, g_yyz_0_xxxxxxxy_0, g_yyz_0_xxxxxxxz_0, g_yyz_0_xxxxxxyy_0, g_yyz_0_xxxxxxyz_0, g_yyz_0_xxxxxxzz_0, g_yyz_0_xxxxxyyy_0, g_yyz_0_xxxxxyyz_0, g_yyz_0_xxxxxyzz_0, g_yyz_0_xxxxxzzz_0, g_yyz_0_xxxxyyyy_0, g_yyz_0_xxxxyyyz_0, g_yyz_0_xxxxyyzz_0, g_yyz_0_xxxxyzzz_0, g_yyz_0_xxxxzzzz_0, g_yyz_0_xxxyyyyy_0, g_yyz_0_xxxyyyyz_0, g_yyz_0_xxxyyyzz_0, g_yyz_0_xxxyyzzz_0, g_yyz_0_xxxyzzzz_0, g_yyz_0_xxxzzzzz_0, g_yyz_0_xxyyyyyy_0, g_yyz_0_xxyyyyyz_0, g_yyz_0_xxyyyyzz_0, g_yyz_0_xxyyyzzz_0, g_yyz_0_xxyyzzzz_0, g_yyz_0_xxyzzzzz_0, g_yyz_0_xxzzzzzz_0, g_yyz_0_xyyyyyyy_0, g_yyz_0_xyyyyyyz_0, g_yyz_0_xyyyyyzz_0, g_yyz_0_xyyyyzzz_0, g_yyz_0_xyyyzzzz_0, g_yyz_0_xyyzzzzz_0, g_yyz_0_xyzzzzzz_0, g_yyz_0_xzzzzzzz_0, g_yyz_0_yyyyyyyy_0, g_yyz_0_yyyyyyyz_0, g_yyz_0_yyyyyyzz_0, g_yyz_0_yyyyyzzz_0, g_yyz_0_yyyyzzzz_0, g_yyz_0_yyyzzzzz_0, g_yyz_0_yyzzzzzz_0, g_yyz_0_yzzzzzzz_0, g_yyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xxxxxxxx_0[i] = g_yy_0_xxxxxxxx_1[i] * wa_z[i];

        g_yyz_0_xxxxxxxy_0[i] = g_yy_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyz_0_xxxxxxxz_0[i] = g_yy_0_xxxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxxz_1[i] * wa_z[i];

        g_yyz_0_xxxxxxyy_0[i] = g_yy_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyz_0_xxxxxxyz_0[i] = g_yy_0_xxxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxxyz_1[i] * wa_z[i];

        g_yyz_0_xxxxxxzz_0[i] = 2.0 * g_yy_0_xxxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxxzz_1[i] * wa_z[i];

        g_yyz_0_xxxxxyyy_0[i] = g_yy_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyz_0_xxxxxyyz_0[i] = g_yy_0_xxxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxxyyz_1[i] * wa_z[i];

        g_yyz_0_xxxxxyzz_0[i] = 2.0 * g_yy_0_xxxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxxyzz_1[i] * wa_z[i];

        g_yyz_0_xxxxxzzz_0[i] = 3.0 * g_yy_0_xxxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxxzzz_1[i] * wa_z[i];

        g_yyz_0_xxxxyyyy_0[i] = g_yy_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyz_0_xxxxyyyz_0[i] = g_yy_0_xxxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyyz_1[i] * wa_z[i];

        g_yyz_0_xxxxyyzz_0[i] = 2.0 * g_yy_0_xxxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyzz_1[i] * wa_z[i];

        g_yyz_0_xxxxyzzz_0[i] = 3.0 * g_yy_0_xxxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzzz_1[i] * wa_z[i];

        g_yyz_0_xxxxzzzz_0[i] = 4.0 * g_yy_0_xxxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxxzzzz_1[i] * wa_z[i];

        g_yyz_0_xxxyyyyy_0[i] = g_yy_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyz_0_xxxyyyyz_0[i] = g_yy_0_xxxyyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyyz_1[i] * wa_z[i];

        g_yyz_0_xxxyyyzz_0[i] = 2.0 * g_yy_0_xxxyyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyzz_1[i] * wa_z[i];

        g_yyz_0_xxxyyzzz_0[i] = 3.0 * g_yy_0_xxxyyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzzz_1[i] * wa_z[i];

        g_yyz_0_xxxyzzzz_0[i] = 4.0 * g_yy_0_xxxyzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzzz_1[i] * wa_z[i];

        g_yyz_0_xxxzzzzz_0[i] = 5.0 * g_yy_0_xxxzzzz_1[i] * fi_acd_0 + g_yy_0_xxxzzzzz_1[i] * wa_z[i];

        g_yyz_0_xxyyyyyy_0[i] = g_yy_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyz_0_xxyyyyyz_0[i] = g_yy_0_xxyyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyyz_1[i] * wa_z[i];

        g_yyz_0_xxyyyyzz_0[i] = 2.0 * g_yy_0_xxyyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyzz_1[i] * wa_z[i];

        g_yyz_0_xxyyyzzz_0[i] = 3.0 * g_yy_0_xxyyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzzz_1[i] * wa_z[i];

        g_yyz_0_xxyyzzzz_0[i] = 4.0 * g_yy_0_xxyyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzzz_1[i] * wa_z[i];

        g_yyz_0_xxyzzzzz_0[i] = 5.0 * g_yy_0_xxyzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzzz_1[i] * wa_z[i];

        g_yyz_0_xxzzzzzz_0[i] = 6.0 * g_yy_0_xxzzzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyyyyy_0[i] = g_yy_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyz_0_xyyyyyyz_0[i] = g_yy_0_xyyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyyz_1[i] * wa_z[i];

        g_yyz_0_xyyyyyzz_0[i] = 2.0 * g_yy_0_xyyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyzz_1[i] * wa_z[i];

        g_yyz_0_xyyyyzzz_0[i] = 3.0 * g_yy_0_xyyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyzzzz_0[i] = 4.0 * g_yy_0_xyyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzzz_1[i] * wa_z[i];

        g_yyz_0_xyyzzzzz_0[i] = 5.0 * g_yy_0_xyyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzzz_1[i] * wa_z[i];

        g_yyz_0_xyzzzzzz_0[i] = 6.0 * g_yy_0_xyzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzzz_1[i] * wa_z[i];

        g_yyz_0_xzzzzzzz_0[i] = 7.0 * g_yy_0_xzzzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyyyyy_0[i] = g_yy_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyz_0_yyyyyyyz_0[i] = g_yy_0_yyyyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyyyz_1[i] * wa_z[i];

        g_yyz_0_yyyyyyzz_0[i] = 2.0 * g_yy_0_yyyyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyyyzz_1[i] * wa_z[i];

        g_yyz_0_yyyyyzzz_0[i] = 3.0 * g_yy_0_yyyyyzz_1[i] * fi_acd_0 + g_yy_0_yyyyyzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyzzzz_0[i] = 4.0 * g_yy_0_yyyyzzz_1[i] * fi_acd_0 + g_yy_0_yyyyzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyzzzzz_0[i] = 5.0 * g_yy_0_yyyzzzz_1[i] * fi_acd_0 + g_yy_0_yyyzzzzz_1[i] * wa_z[i];

        g_yyz_0_yyzzzzzz_0[i] = 6.0 * g_yy_0_yyzzzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzzzz_1[i] * wa_z[i];

        g_yyz_0_yzzzzzzz_0[i] = 7.0 * g_yy_0_yzzzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzzzz_1[i] * wa_z[i];

        g_yyz_0_zzzzzzzz_0[i] = 8.0 * g_yy_0_zzzzzzz_1[i] * fi_acd_0 + g_yy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 360-405 components of targeted buffer : FSL

    auto g_yzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 360);

    auto g_yzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 361);

    auto g_yzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 362);

    auto g_yzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 363);

    auto g_yzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 364);

    auto g_yzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 365);

    auto g_yzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 366);

    auto g_yzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 367);

    auto g_yzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 368);

    auto g_yzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 369);

    auto g_yzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 370);

    auto g_yzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 371);

    auto g_yzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 372);

    auto g_yzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 373);

    auto g_yzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 374);

    auto g_yzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 375);

    auto g_yzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 376);

    auto g_yzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 377);

    auto g_yzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 378);

    auto g_yzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 379);

    auto g_yzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 380);

    auto g_yzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 381);

    auto g_yzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 382);

    auto g_yzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 383);

    auto g_yzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 384);

    auto g_yzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 385);

    auto g_yzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 386);

    auto g_yzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 387);

    auto g_yzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 388);

    auto g_yzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 389);

    auto g_yzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 390);

    auto g_yzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 391);

    auto g_yzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 392);

    auto g_yzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 393);

    auto g_yzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 394);

    auto g_yzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 395);

    auto g_yzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 396);

    auto g_yzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 397);

    auto g_yzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 398);

    auto g_yzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 399);

    auto g_yzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 400);

    auto g_yzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 401);

    auto g_yzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 402);

    auto g_yzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 403);

    auto g_yzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 404);

    #pragma omp simd aligned(g_yzz_0_xxxxxxxx_0, g_yzz_0_xxxxxxxy_0, g_yzz_0_xxxxxxxz_0, g_yzz_0_xxxxxxyy_0, g_yzz_0_xxxxxxyz_0, g_yzz_0_xxxxxxzz_0, g_yzz_0_xxxxxyyy_0, g_yzz_0_xxxxxyyz_0, g_yzz_0_xxxxxyzz_0, g_yzz_0_xxxxxzzz_0, g_yzz_0_xxxxyyyy_0, g_yzz_0_xxxxyyyz_0, g_yzz_0_xxxxyyzz_0, g_yzz_0_xxxxyzzz_0, g_yzz_0_xxxxzzzz_0, g_yzz_0_xxxyyyyy_0, g_yzz_0_xxxyyyyz_0, g_yzz_0_xxxyyyzz_0, g_yzz_0_xxxyyzzz_0, g_yzz_0_xxxyzzzz_0, g_yzz_0_xxxzzzzz_0, g_yzz_0_xxyyyyyy_0, g_yzz_0_xxyyyyyz_0, g_yzz_0_xxyyyyzz_0, g_yzz_0_xxyyyzzz_0, g_yzz_0_xxyyzzzz_0, g_yzz_0_xxyzzzzz_0, g_yzz_0_xxzzzzzz_0, g_yzz_0_xyyyyyyy_0, g_yzz_0_xyyyyyyz_0, g_yzz_0_xyyyyyzz_0, g_yzz_0_xyyyyzzz_0, g_yzz_0_xyyyzzzz_0, g_yzz_0_xyyzzzzz_0, g_yzz_0_xyzzzzzz_0, g_yzz_0_xzzzzzzz_0, g_yzz_0_yyyyyyyy_0, g_yzz_0_yyyyyyyz_0, g_yzz_0_yyyyyyzz_0, g_yzz_0_yyyyyzzz_0, g_yzz_0_yyyyzzzz_0, g_yzz_0_yyyzzzzz_0, g_yzz_0_yyzzzzzz_0, g_yzz_0_yzzzzzzz_0, g_yzz_0_zzzzzzzz_0, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxy_1, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxy_1, g_zz_0_xxxxxxyy_1, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxyy_1, g_zz_0_xxxxxyyy_1, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxzz_1, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxyyy_1, g_zz_0_xxxxyyyy_1, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyzz_1, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxzzz_1, g_zz_0_xxxxzzzz_1, g_zz_0_xxxyyyy_1, g_zz_0_xxxyyyyy_1, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyzz_1, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyzzz_1, g_zz_0_xxxyzzzz_1, g_zz_0_xxxzzzz_1, g_zz_0_xxxzzzzz_1, g_zz_0_xxyyyyy_1, g_zz_0_xxyyyyyy_1, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyzz_1, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyzzz_1, g_zz_0_xxyyzzzz_1, g_zz_0_xxyzzzz_1, g_zz_0_xxyzzzzz_1, g_zz_0_xxzzzzz_1, g_zz_0_xxzzzzzz_1, g_zz_0_xyyyyyy_1, g_zz_0_xyyyyyyy_1, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyzz_1, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyzzz_1, g_zz_0_xyyyzzzz_1, g_zz_0_xyyzzzz_1, g_zz_0_xyyzzzzz_1, g_zz_0_xyzzzzz_1, g_zz_0_xyzzzzzz_1, g_zz_0_xzzzzzz_1, g_zz_0_xzzzzzzz_1, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyzz_1, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyzzz_1, g_zz_0_yyyyzzzz_1, g_zz_0_yyyzzzz_1, g_zz_0_yyyzzzzz_1, g_zz_0_yyzzzzz_1, g_zz_0_yyzzzzzz_1, g_zz_0_yzzzzzz_1, g_zz_0_yzzzzzzz_1, g_zz_0_zzzzzzz_1, g_zz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xxxxxxxx_0[i] = g_zz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yzz_0_xxxxxxxy_0[i] = g_zz_0_xxxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxxy_1[i] * wa_y[i];

        g_yzz_0_xxxxxxxz_0[i] = g_zz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yzz_0_xxxxxxyy_0[i] = 2.0 * g_zz_0_xxxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxxyy_1[i] * wa_y[i];

        g_yzz_0_xxxxxxyz_0[i] = g_zz_0_xxxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yzz_0_xxxxxxzz_0[i] = g_zz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yzz_0_xxxxxyyy_0[i] = 3.0 * g_zz_0_xxxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxxyyy_1[i] * wa_y[i];

        g_yzz_0_xxxxxyyz_0[i] = 2.0 * g_zz_0_xxxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yzz_0_xxxxxyzz_0[i] = g_zz_0_xxxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yzz_0_xxxxxzzz_0[i] = g_zz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yzz_0_xxxxyyyy_0[i] = 4.0 * g_zz_0_xxxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyyy_1[i] * wa_y[i];

        g_yzz_0_xxxxyyyz_0[i] = 3.0 * g_zz_0_xxxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yzz_0_xxxxyyzz_0[i] = 2.0 * g_zz_0_xxxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yzz_0_xxxxyzzz_0[i] = g_zz_0_xxxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yzz_0_xxxxzzzz_0[i] = g_zz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yzz_0_xxxyyyyy_0[i] = 5.0 * g_zz_0_xxxyyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyyy_1[i] * wa_y[i];

        g_yzz_0_xxxyyyyz_0[i] = 4.0 * g_zz_0_xxxyyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yzz_0_xxxyyyzz_0[i] = 3.0 * g_zz_0_xxxyyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yzz_0_xxxyyzzz_0[i] = 2.0 * g_zz_0_xxxyzzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yzz_0_xxxyzzzz_0[i] = g_zz_0_xxxzzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yzz_0_xxxzzzzz_0[i] = g_zz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yzz_0_xxyyyyyy_0[i] = 6.0 * g_zz_0_xxyyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyyy_1[i] * wa_y[i];

        g_yzz_0_xxyyyyyz_0[i] = 5.0 * g_zz_0_xxyyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yzz_0_xxyyyyzz_0[i] = 4.0 * g_zz_0_xxyyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yzz_0_xxyyyzzz_0[i] = 3.0 * g_zz_0_xxyyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yzz_0_xxyyzzzz_0[i] = 2.0 * g_zz_0_xxyzzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yzz_0_xxyzzzzz_0[i] = g_zz_0_xxzzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yzz_0_xxzzzzzz_0[i] = g_zz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyyyyy_0[i] = 7.0 * g_zz_0_xyyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyyy_1[i] * wa_y[i];

        g_yzz_0_xyyyyyyz_0[i] = 6.0 * g_zz_0_xyyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yzz_0_xyyyyyzz_0[i] = 5.0 * g_zz_0_xyyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yzz_0_xyyyyzzz_0[i] = 4.0 * g_zz_0_xyyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyzzzz_0[i] = 3.0 * g_zz_0_xyyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yzz_0_xyyzzzzz_0[i] = 2.0 * g_zz_0_xyzzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yzz_0_xyzzzzzz_0[i] = g_zz_0_xzzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yzz_0_xzzzzzzz_0[i] = g_zz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyyyyy_0[i] = 8.0 * g_zz_0_yyyyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyyyy_1[i] * wa_y[i];

        g_yzz_0_yyyyyyyz_0[i] = 7.0 * g_zz_0_yyyyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yzz_0_yyyyyyzz_0[i] = 6.0 * g_zz_0_yyyyyzz_1[i] * fi_acd_0 + g_zz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yzz_0_yyyyyzzz_0[i] = 5.0 * g_zz_0_yyyyzzz_1[i] * fi_acd_0 + g_zz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyzzzz_0[i] = 4.0 * g_zz_0_yyyzzzz_1[i] * fi_acd_0 + g_zz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyzzzzz_0[i] = 3.0 * g_zz_0_yyzzzzz_1[i] * fi_acd_0 + g_zz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yzz_0_yyzzzzzz_0[i] = 2.0 * g_zz_0_yzzzzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yzz_0_yzzzzzzz_0[i] = g_zz_0_zzzzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yzz_0_zzzzzzzz_0[i] = g_zz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 405-450 components of targeted buffer : FSL

    auto g_zzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 405);

    auto g_zzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 406);

    auto g_zzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 407);

    auto g_zzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 408);

    auto g_zzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 409);

    auto g_zzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 410);

    auto g_zzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 411);

    auto g_zzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 412);

    auto g_zzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 413);

    auto g_zzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 414);

    auto g_zzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 415);

    auto g_zzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 416);

    auto g_zzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 417);

    auto g_zzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 418);

    auto g_zzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 419);

    auto g_zzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 420);

    auto g_zzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 421);

    auto g_zzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 422);

    auto g_zzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 423);

    auto g_zzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 424);

    auto g_zzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 425);

    auto g_zzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 426);

    auto g_zzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 427);

    auto g_zzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 428);

    auto g_zzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 429);

    auto g_zzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 430);

    auto g_zzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 431);

    auto g_zzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 432);

    auto g_zzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 433);

    auto g_zzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 434);

    auto g_zzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 435);

    auto g_zzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 436);

    auto g_zzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 437);

    auto g_zzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 438);

    auto g_zzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 439);

    auto g_zzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 440);

    auto g_zzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 441);

    auto g_zzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 442);

    auto g_zzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 443);

    auto g_zzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 444);

    auto g_zzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 445);

    auto g_zzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 446);

    auto g_zzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 447);

    auto g_zzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 448);

    auto g_zzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 449);

    #pragma omp simd aligned(g_z_0_xxxxxxxx_0, g_z_0_xxxxxxxx_1, g_z_0_xxxxxxxy_0, g_z_0_xxxxxxxy_1, g_z_0_xxxxxxxz_0, g_z_0_xxxxxxxz_1, g_z_0_xxxxxxyy_0, g_z_0_xxxxxxyy_1, g_z_0_xxxxxxyz_0, g_z_0_xxxxxxyz_1, g_z_0_xxxxxxzz_0, g_z_0_xxxxxxzz_1, g_z_0_xxxxxyyy_0, g_z_0_xxxxxyyy_1, g_z_0_xxxxxyyz_0, g_z_0_xxxxxyyz_1, g_z_0_xxxxxyzz_0, g_z_0_xxxxxyzz_1, g_z_0_xxxxxzzz_0, g_z_0_xxxxxzzz_1, g_z_0_xxxxyyyy_0, g_z_0_xxxxyyyy_1, g_z_0_xxxxyyyz_0, g_z_0_xxxxyyyz_1, g_z_0_xxxxyyzz_0, g_z_0_xxxxyyzz_1, g_z_0_xxxxyzzz_0, g_z_0_xxxxyzzz_1, g_z_0_xxxxzzzz_0, g_z_0_xxxxzzzz_1, g_z_0_xxxyyyyy_0, g_z_0_xxxyyyyy_1, g_z_0_xxxyyyyz_0, g_z_0_xxxyyyyz_1, g_z_0_xxxyyyzz_0, g_z_0_xxxyyyzz_1, g_z_0_xxxyyzzz_0, g_z_0_xxxyyzzz_1, g_z_0_xxxyzzzz_0, g_z_0_xxxyzzzz_1, g_z_0_xxxzzzzz_0, g_z_0_xxxzzzzz_1, g_z_0_xxyyyyyy_0, g_z_0_xxyyyyyy_1, g_z_0_xxyyyyyz_0, g_z_0_xxyyyyyz_1, g_z_0_xxyyyyzz_0, g_z_0_xxyyyyzz_1, g_z_0_xxyyyzzz_0, g_z_0_xxyyyzzz_1, g_z_0_xxyyzzzz_0, g_z_0_xxyyzzzz_1, g_z_0_xxyzzzzz_0, g_z_0_xxyzzzzz_1, g_z_0_xxzzzzzz_0, g_z_0_xxzzzzzz_1, g_z_0_xyyyyyyy_0, g_z_0_xyyyyyyy_1, g_z_0_xyyyyyyz_0, g_z_0_xyyyyyyz_1, g_z_0_xyyyyyzz_0, g_z_0_xyyyyyzz_1, g_z_0_xyyyyzzz_0, g_z_0_xyyyyzzz_1, g_z_0_xyyyzzzz_0, g_z_0_xyyyzzzz_1, g_z_0_xyyzzzzz_0, g_z_0_xyyzzzzz_1, g_z_0_xyzzzzzz_0, g_z_0_xyzzzzzz_1, g_z_0_xzzzzzzz_0, g_z_0_xzzzzzzz_1, g_z_0_yyyyyyyy_0, g_z_0_yyyyyyyy_1, g_z_0_yyyyyyyz_0, g_z_0_yyyyyyyz_1, g_z_0_yyyyyyzz_0, g_z_0_yyyyyyzz_1, g_z_0_yyyyyzzz_0, g_z_0_yyyyyzzz_1, g_z_0_yyyyzzzz_0, g_z_0_yyyyzzzz_1, g_z_0_yyyzzzzz_0, g_z_0_yyyzzzzz_1, g_z_0_yyzzzzzz_0, g_z_0_yyzzzzzz_1, g_z_0_yzzzzzzz_0, g_z_0_yzzzzzzz_1, g_z_0_zzzzzzzz_0, g_z_0_zzzzzzzz_1, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxy_1, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxy_1, g_zz_0_xxxxxxyy_1, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxyy_1, g_zz_0_xxxxxyyy_1, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxzz_1, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxyyy_1, g_zz_0_xxxxyyyy_1, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyzz_1, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxzzz_1, g_zz_0_xxxxzzzz_1, g_zz_0_xxxyyyy_1, g_zz_0_xxxyyyyy_1, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyzz_1, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyzzz_1, g_zz_0_xxxyzzzz_1, g_zz_0_xxxzzzz_1, g_zz_0_xxxzzzzz_1, g_zz_0_xxyyyyy_1, g_zz_0_xxyyyyyy_1, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyzz_1, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyzzz_1, g_zz_0_xxyyzzzz_1, g_zz_0_xxyzzzz_1, g_zz_0_xxyzzzzz_1, g_zz_0_xxzzzzz_1, g_zz_0_xxzzzzzz_1, g_zz_0_xyyyyyy_1, g_zz_0_xyyyyyyy_1, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyzz_1, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyzzz_1, g_zz_0_xyyyzzzz_1, g_zz_0_xyyzzzz_1, g_zz_0_xyyzzzzz_1, g_zz_0_xyzzzzz_1, g_zz_0_xyzzzzzz_1, g_zz_0_xzzzzzz_1, g_zz_0_xzzzzzzz_1, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyzz_1, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyzzz_1, g_zz_0_yyyyzzzz_1, g_zz_0_yyyzzzz_1, g_zz_0_yyyzzzzz_1, g_zz_0_yyzzzzz_1, g_zz_0_yyzzzzzz_1, g_zz_0_yzzzzzz_1, g_zz_0_yzzzzzzz_1, g_zz_0_zzzzzzz_1, g_zz_0_zzzzzzzz_1, g_zzz_0_xxxxxxxx_0, g_zzz_0_xxxxxxxy_0, g_zzz_0_xxxxxxxz_0, g_zzz_0_xxxxxxyy_0, g_zzz_0_xxxxxxyz_0, g_zzz_0_xxxxxxzz_0, g_zzz_0_xxxxxyyy_0, g_zzz_0_xxxxxyyz_0, g_zzz_0_xxxxxyzz_0, g_zzz_0_xxxxxzzz_0, g_zzz_0_xxxxyyyy_0, g_zzz_0_xxxxyyyz_0, g_zzz_0_xxxxyyzz_0, g_zzz_0_xxxxyzzz_0, g_zzz_0_xxxxzzzz_0, g_zzz_0_xxxyyyyy_0, g_zzz_0_xxxyyyyz_0, g_zzz_0_xxxyyyzz_0, g_zzz_0_xxxyyzzz_0, g_zzz_0_xxxyzzzz_0, g_zzz_0_xxxzzzzz_0, g_zzz_0_xxyyyyyy_0, g_zzz_0_xxyyyyyz_0, g_zzz_0_xxyyyyzz_0, g_zzz_0_xxyyyzzz_0, g_zzz_0_xxyyzzzz_0, g_zzz_0_xxyzzzzz_0, g_zzz_0_xxzzzzzz_0, g_zzz_0_xyyyyyyy_0, g_zzz_0_xyyyyyyz_0, g_zzz_0_xyyyyyzz_0, g_zzz_0_xyyyyzzz_0, g_zzz_0_xyyyzzzz_0, g_zzz_0_xyyzzzzz_0, g_zzz_0_xyzzzzzz_0, g_zzz_0_xzzzzzzz_0, g_zzz_0_yyyyyyyy_0, g_zzz_0_yyyyyyyz_0, g_zzz_0_yyyyyyzz_0, g_zzz_0_yyyyyzzz_0, g_zzz_0_yyyyzzzz_0, g_zzz_0_yyyzzzzz_0, g_zzz_0_yyzzzzzz_0, g_zzz_0_yzzzzzzz_0, g_zzz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xxxxxxxx_0[i] = 2.0 * g_z_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxx_1[i] * fz_be_0 + g_zz_0_xxxxxxxx_1[i] * wa_z[i];

        g_zzz_0_xxxxxxxy_0[i] = 2.0 * g_z_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxy_1[i] * fz_be_0 + g_zz_0_xxxxxxxy_1[i] * wa_z[i];

        g_zzz_0_xxxxxxxz_0[i] = 2.0 * g_z_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxz_1[i] * fz_be_0 + g_zz_0_xxxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxxz_1[i] * wa_z[i];

        g_zzz_0_xxxxxxyy_0[i] = 2.0 * g_z_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxyy_1[i] * fz_be_0 + g_zz_0_xxxxxxyy_1[i] * wa_z[i];

        g_zzz_0_xxxxxxyz_0[i] = 2.0 * g_z_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxyz_1[i] * fz_be_0 + g_zz_0_xxxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxxyz_1[i] * wa_z[i];

        g_zzz_0_xxxxxxzz_0[i] = 2.0 * g_z_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxxzz_1[i] * wa_z[i];

        g_zzz_0_xxxxxyyy_0[i] = 2.0 * g_z_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyyy_1[i] * fz_be_0 + g_zz_0_xxxxxyyy_1[i] * wa_z[i];

        g_zzz_0_xxxxxyyz_0[i] = 2.0 * g_z_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyyz_1[i] * fz_be_0 + g_zz_0_xxxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxxyyz_1[i] * wa_z[i];

        g_zzz_0_xxxxxyzz_0[i] = 2.0 * g_z_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxxyzz_1[i] * wa_z[i];

        g_zzz_0_xxxxxzzz_0[i] = 2.0 * g_z_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxxzzz_1[i] * wa_z[i];

        g_zzz_0_xxxxyyyy_0[i] = 2.0 * g_z_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyyy_1[i] * fz_be_0 + g_zz_0_xxxxyyyy_1[i] * wa_z[i];

        g_zzz_0_xxxxyyyz_0[i] = 2.0 * g_z_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyyz_1[i] * fz_be_0 + g_zz_0_xxxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyyz_1[i] * wa_z[i];

        g_zzz_0_xxxxyyzz_0[i] = 2.0 * g_z_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyzz_1[i] * wa_z[i];

        g_zzz_0_xxxxyzzz_0[i] = 2.0 * g_z_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzzz_1[i] * wa_z[i];

        g_zzz_0_xxxxzzzz_0[i] = 2.0 * g_z_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxxzzzz_1[i] * wa_z[i];

        g_zzz_0_xxxyyyyy_0[i] = 2.0 * g_z_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyyy_1[i] * fz_be_0 + g_zz_0_xxxyyyyy_1[i] * wa_z[i];

        g_zzz_0_xxxyyyyz_0[i] = 2.0 * g_z_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyyz_1[i] * fz_be_0 + g_zz_0_xxxyyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyyz_1[i] * wa_z[i];

        g_zzz_0_xxxyyyzz_0[i] = 2.0 * g_z_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxyyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyzz_1[i] * wa_z[i];

        g_zzz_0_xxxyyzzz_0[i] = 2.0 * g_z_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxyyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzzz_1[i] * wa_z[i];

        g_zzz_0_xxxyzzzz_0[i] = 2.0 * g_z_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxxyzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzzz_1[i] * wa_z[i];

        g_zzz_0_xxxzzzzz_0[i] = 2.0 * g_z_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xxxzzzz_1[i] * fi_acd_0 + g_zz_0_xxxzzzzz_1[i] * wa_z[i];

        g_zzz_0_xxyyyyyy_0[i] = 2.0 * g_z_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyyy_1[i] * fz_be_0 + g_zz_0_xxyyyyyy_1[i] * wa_z[i];

        g_zzz_0_xxyyyyyz_0[i] = 2.0 * g_z_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyyz_1[i] * fz_be_0 + g_zz_0_xxyyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyyz_1[i] * wa_z[i];

        g_zzz_0_xxyyyyzz_0[i] = 2.0 * g_z_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxyyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyzz_1[i] * wa_z[i];

        g_zzz_0_xxyyyzzz_0[i] = 2.0 * g_z_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxyyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzzz_1[i] * wa_z[i];

        g_zzz_0_xxyyzzzz_0[i] = 2.0 * g_z_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxyyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzzz_1[i] * wa_z[i];

        g_zzz_0_xxyzzzzz_0[i] = 2.0 * g_z_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xxyzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzzz_1[i] * wa_z[i];

        g_zzz_0_xxzzzzzz_0[i] = 2.0 * g_z_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_xxzzzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyyyyy_0[i] = 2.0 * g_z_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyyy_1[i] * fz_be_0 + g_zz_0_xyyyyyyy_1[i] * wa_z[i];

        g_zzz_0_xyyyyyyz_0[i] = 2.0 * g_z_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyyz_1[i] * fz_be_0 + g_zz_0_xyyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyyz_1[i] * wa_z[i];

        g_zzz_0_xyyyyyzz_0[i] = 2.0 * g_z_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xyyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyzz_1[i] * wa_z[i];

        g_zzz_0_xyyyyzzz_0[i] = 2.0 * g_z_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xyyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyzzzz_0[i] = 2.0 * g_z_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xyyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzzz_1[i] * wa_z[i];

        g_zzz_0_xyyzzzzz_0[i] = 2.0 * g_z_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xyyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzzz_1[i] * wa_z[i];

        g_zzz_0_xyzzzzzz_0[i] = 2.0 * g_z_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_xyzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzzz_1[i] * wa_z[i];

        g_zzz_0_xzzzzzzz_0[i] = 2.0 * g_z_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zz_0_xzzzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyyyyy_0[i] = 2.0 * g_z_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyyy_1[i] * fz_be_0 + g_zz_0_yyyyyyyy_1[i] * wa_z[i];

        g_zzz_0_yyyyyyyz_0[i] = 2.0 * g_z_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyyz_1[i] * fz_be_0 + g_zz_0_yyyyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyyyz_1[i] * wa_z[i];

        g_zzz_0_yyyyyyzz_0[i] = 2.0 * g_z_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_yyyyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyyyzz_1[i] * wa_z[i];

        g_zzz_0_yyyyyzzz_0[i] = 2.0 * g_z_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_yyyyyzz_1[i] * fi_acd_0 + g_zz_0_yyyyyzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyzzzz_0[i] = 2.0 * g_z_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_yyyyzzz_1[i] * fi_acd_0 + g_zz_0_yyyyzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyzzzzz_0[i] = 2.0 * g_z_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_yyyzzzz_1[i] * fi_acd_0 + g_zz_0_yyyzzzzz_1[i] * wa_z[i];

        g_zzz_0_yyzzzzzz_0[i] = 2.0 * g_z_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_yyzzzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzzzz_1[i] * wa_z[i];

        g_zzz_0_yzzzzzzz_0[i] = 2.0 * g_z_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zz_0_yzzzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzzzz_1[i] * wa_z[i];

        g_zzz_0_zzzzzzzz_0[i] = 2.0 * g_z_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_zzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zz_0_zzzzzzz_1[i] * fi_acd_0 + g_zz_0_zzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

