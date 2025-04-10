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

#include "ThreeCenterElectronRepulsionPrimRecFSM.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsm,
                                 size_t idx_eri_0_psm,
                                 size_t idx_eri_1_psm,
                                 size_t idx_eri_1_dsl,
                                 size_t idx_eri_1_dsm,
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

    /// Set up components of auxilary buffer : PSM

    auto g_x_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_psm);

    auto g_x_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_psm + 1);

    auto g_x_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_psm + 2);

    auto g_x_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_psm + 3);

    auto g_x_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_psm + 4);

    auto g_x_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_psm + 5);

    auto g_x_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_psm + 6);

    auto g_x_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_psm + 7);

    auto g_x_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_psm + 8);

    auto g_x_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_psm + 9);

    auto g_x_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_psm + 10);

    auto g_x_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_psm + 11);

    auto g_x_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_psm + 12);

    auto g_x_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_psm + 13);

    auto g_x_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_psm + 14);

    auto g_x_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_psm + 15);

    auto g_x_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_psm + 16);

    auto g_x_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_psm + 17);

    auto g_x_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_psm + 18);

    auto g_x_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_psm + 19);

    auto g_x_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_psm + 20);

    auto g_x_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 21);

    auto g_x_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 22);

    auto g_x_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 23);

    auto g_x_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 24);

    auto g_x_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 25);

    auto g_x_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 26);

    auto g_x_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 27);

    auto g_x_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 28);

    auto g_x_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 29);

    auto g_x_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 30);

    auto g_x_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 31);

    auto g_x_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 32);

    auto g_x_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 33);

    auto g_x_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 34);

    auto g_x_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 35);

    auto g_x_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 36);

    auto g_x_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 37);

    auto g_x_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 38);

    auto g_x_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 39);

    auto g_x_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 40);

    auto g_x_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 41);

    auto g_x_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 42);

    auto g_x_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 43);

    auto g_x_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 44);

    auto g_x_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 45);

    auto g_x_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 46);

    auto g_x_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 47);

    auto g_x_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 48);

    auto g_x_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 49);

    auto g_x_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 50);

    auto g_x_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 51);

    auto g_x_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 52);

    auto g_x_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 53);

    auto g_x_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 54);

    auto g_y_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_psm + 55);

    auto g_y_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_psm + 56);

    auto g_y_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_psm + 57);

    auto g_y_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_psm + 58);

    auto g_y_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_psm + 59);

    auto g_y_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_psm + 60);

    auto g_y_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_psm + 61);

    auto g_y_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_psm + 62);

    auto g_y_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_psm + 63);

    auto g_y_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_psm + 64);

    auto g_y_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_psm + 65);

    auto g_y_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_psm + 66);

    auto g_y_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_psm + 67);

    auto g_y_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_psm + 68);

    auto g_y_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_psm + 69);

    auto g_y_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_psm + 70);

    auto g_y_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_psm + 71);

    auto g_y_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_psm + 72);

    auto g_y_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_psm + 73);

    auto g_y_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_psm + 74);

    auto g_y_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_psm + 75);

    auto g_y_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 76);

    auto g_y_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 77);

    auto g_y_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 78);

    auto g_y_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 79);

    auto g_y_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 80);

    auto g_y_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 81);

    auto g_y_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 82);

    auto g_y_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 83);

    auto g_y_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 84);

    auto g_y_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 85);

    auto g_y_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 86);

    auto g_y_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 87);

    auto g_y_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 88);

    auto g_y_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 89);

    auto g_y_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 90);

    auto g_y_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 91);

    auto g_y_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 92);

    auto g_y_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 93);

    auto g_y_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 94);

    auto g_y_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 95);

    auto g_y_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 96);

    auto g_y_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 97);

    auto g_y_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 98);

    auto g_y_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 99);

    auto g_y_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 100);

    auto g_y_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 101);

    auto g_y_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 102);

    auto g_y_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 103);

    auto g_y_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 104);

    auto g_y_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 105);

    auto g_y_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 106);

    auto g_y_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 107);

    auto g_y_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 108);

    auto g_y_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 109);

    auto g_z_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_psm + 110);

    auto g_z_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_psm + 111);

    auto g_z_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_psm + 112);

    auto g_z_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_psm + 113);

    auto g_z_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_psm + 114);

    auto g_z_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_psm + 115);

    auto g_z_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_psm + 116);

    auto g_z_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_psm + 117);

    auto g_z_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_psm + 118);

    auto g_z_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_psm + 119);

    auto g_z_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_psm + 120);

    auto g_z_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_psm + 121);

    auto g_z_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_psm + 122);

    auto g_z_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_psm + 123);

    auto g_z_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_psm + 124);

    auto g_z_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_psm + 125);

    auto g_z_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_psm + 126);

    auto g_z_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_psm + 127);

    auto g_z_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_psm + 128);

    auto g_z_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_psm + 129);

    auto g_z_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_psm + 130);

    auto g_z_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 131);

    auto g_z_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 132);

    auto g_z_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 133);

    auto g_z_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 134);

    auto g_z_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 135);

    auto g_z_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 136);

    auto g_z_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 137);

    auto g_z_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 138);

    auto g_z_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 139);

    auto g_z_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 140);

    auto g_z_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 141);

    auto g_z_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 142);

    auto g_z_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 143);

    auto g_z_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 144);

    auto g_z_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 145);

    auto g_z_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 146);

    auto g_z_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 147);

    auto g_z_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 148);

    auto g_z_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 149);

    auto g_z_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 150);

    auto g_z_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 151);

    auto g_z_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 152);

    auto g_z_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 153);

    auto g_z_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 154);

    auto g_z_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 155);

    auto g_z_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 156);

    auto g_z_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 157);

    auto g_z_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 158);

    auto g_z_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 159);

    auto g_z_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 160);

    auto g_z_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 161);

    auto g_z_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 162);

    auto g_z_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 163);

    auto g_z_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 164);

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

    auto g_yz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 217);

    auto g_yz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 218);

    auto g_yz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 219);

    auto g_yz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 220);

    auto g_yz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 221);

    auto g_yz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 222);

    auto g_yz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 223);

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

    /// Set up components of auxilary buffer : DSM

    auto g_xx_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_dsm);

    auto g_xx_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_dsm + 1);

    auto g_xx_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_dsm + 2);

    auto g_xx_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_dsm + 3);

    auto g_xx_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_dsm + 4);

    auto g_xx_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_dsm + 5);

    auto g_xx_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_dsm + 6);

    auto g_xx_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_dsm + 7);

    auto g_xx_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_dsm + 8);

    auto g_xx_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_dsm + 9);

    auto g_xx_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_dsm + 10);

    auto g_xx_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_dsm + 11);

    auto g_xx_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_dsm + 12);

    auto g_xx_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_dsm + 13);

    auto g_xx_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_dsm + 14);

    auto g_xx_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 15);

    auto g_xx_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 16);

    auto g_xx_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 17);

    auto g_xx_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 18);

    auto g_xx_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 19);

    auto g_xx_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 20);

    auto g_xx_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 21);

    auto g_xx_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 22);

    auto g_xx_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 23);

    auto g_xx_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 24);

    auto g_xx_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 25);

    auto g_xx_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 26);

    auto g_xx_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 27);

    auto g_xx_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 28);

    auto g_xx_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 29);

    auto g_xx_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 30);

    auto g_xx_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 31);

    auto g_xx_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 32);

    auto g_xx_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 33);

    auto g_xx_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 34);

    auto g_xx_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 35);

    auto g_xx_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 36);

    auto g_xx_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 37);

    auto g_xx_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 38);

    auto g_xx_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 39);

    auto g_xx_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 40);

    auto g_xx_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 41);

    auto g_xx_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 42);

    auto g_xx_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 43);

    auto g_xx_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 44);

    auto g_xx_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 45);

    auto g_xx_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 46);

    auto g_xx_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 47);

    auto g_xx_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 48);

    auto g_xx_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 49);

    auto g_xx_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 50);

    auto g_xx_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 51);

    auto g_xx_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 52);

    auto g_xx_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 53);

    auto g_xx_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 54);

    auto g_xy_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_dsm + 56);

    auto g_xy_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_dsm + 58);

    auto g_xy_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_dsm + 61);

    auto g_xy_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_dsm + 65);

    auto g_xy_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 70);

    auto g_xy_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 76);

    auto g_xy_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 83);

    auto g_xy_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 91);

    auto g_xz_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_dsm + 110);

    auto g_xz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_dsm + 112);

    auto g_xz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_dsm + 115);

    auto g_xz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_dsm + 119);

    auto g_xz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_dsm + 124);

    auto g_xz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 130);

    auto g_xz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 137);

    auto g_xz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 145);

    auto g_xz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 154);

    auto g_yy_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_dsm + 165);

    auto g_yy_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_dsm + 166);

    auto g_yy_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_dsm + 167);

    auto g_yy_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_dsm + 168);

    auto g_yy_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_dsm + 169);

    auto g_yy_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_dsm + 170);

    auto g_yy_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_dsm + 171);

    auto g_yy_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_dsm + 172);

    auto g_yy_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_dsm + 173);

    auto g_yy_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_dsm + 174);

    auto g_yy_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_dsm + 175);

    auto g_yy_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_dsm + 176);

    auto g_yy_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_dsm + 177);

    auto g_yy_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_dsm + 178);

    auto g_yy_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_dsm + 179);

    auto g_yy_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 180);

    auto g_yy_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 181);

    auto g_yy_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 182);

    auto g_yy_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 183);

    auto g_yy_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 184);

    auto g_yy_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 185);

    auto g_yy_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 186);

    auto g_yy_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 187);

    auto g_yy_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 188);

    auto g_yy_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 189);

    auto g_yy_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 190);

    auto g_yy_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 191);

    auto g_yy_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 192);

    auto g_yy_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 193);

    auto g_yy_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 194);

    auto g_yy_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 195);

    auto g_yy_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 196);

    auto g_yy_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 197);

    auto g_yy_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 198);

    auto g_yy_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 199);

    auto g_yy_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 200);

    auto g_yy_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 201);

    auto g_yy_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 202);

    auto g_yy_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 203);

    auto g_yy_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 204);

    auto g_yy_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 205);

    auto g_yy_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 206);

    auto g_yy_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 207);

    auto g_yy_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 208);

    auto g_yy_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 209);

    auto g_yy_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 210);

    auto g_yy_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 211);

    auto g_yy_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 212);

    auto g_yy_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 213);

    auto g_yy_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 214);

    auto g_yy_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 215);

    auto g_yy_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 216);

    auto g_yy_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 217);

    auto g_yy_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 218);

    auto g_yy_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 219);

    auto g_yz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_dsm + 224);

    auto g_yz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_dsm + 227);

    auto g_yz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_dsm + 228);

    auto g_yz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_dsm + 231);

    auto g_yz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_dsm + 232);

    auto g_yz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_dsm + 233);

    auto g_yz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 236);

    auto g_yz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 237);

    auto g_yz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 238);

    auto g_yz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 239);

    auto g_yz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 242);

    auto g_yz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 243);

    auto g_yz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 244);

    auto g_yz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 245);

    auto g_yz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 246);

    auto g_yz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 249);

    auto g_yz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 250);

    auto g_yz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 251);

    auto g_yz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 252);

    auto g_yz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 253);

    auto g_yz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 254);

    auto g_yz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 257);

    auto g_yz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 258);

    auto g_yz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 259);

    auto g_yz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 260);

    auto g_yz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 261);

    auto g_yz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 262);

    auto g_yz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 263);

    auto g_yz_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 265);

    auto g_yz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 266);

    auto g_yz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 267);

    auto g_yz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 268);

    auto g_yz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 269);

    auto g_yz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 270);

    auto g_yz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 271);

    auto g_yz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 272);

    auto g_yz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 273);

    auto g_yz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 274);

    auto g_zz_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_dsm + 275);

    auto g_zz_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_dsm + 276);

    auto g_zz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_dsm + 277);

    auto g_zz_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_dsm + 278);

    auto g_zz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_dsm + 279);

    auto g_zz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_dsm + 280);

    auto g_zz_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_dsm + 281);

    auto g_zz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_dsm + 282);

    auto g_zz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_dsm + 283);

    auto g_zz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_dsm + 284);

    auto g_zz_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_dsm + 285);

    auto g_zz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_dsm + 286);

    auto g_zz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_dsm + 287);

    auto g_zz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_dsm + 288);

    auto g_zz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_dsm + 289);

    auto g_zz_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 290);

    auto g_zz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 291);

    auto g_zz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 292);

    auto g_zz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 293);

    auto g_zz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 294);

    auto g_zz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 295);

    auto g_zz_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 296);

    auto g_zz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 297);

    auto g_zz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 298);

    auto g_zz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 299);

    auto g_zz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 300);

    auto g_zz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 301);

    auto g_zz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 302);

    auto g_zz_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 303);

    auto g_zz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 304);

    auto g_zz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 305);

    auto g_zz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 306);

    auto g_zz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 307);

    auto g_zz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 308);

    auto g_zz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 309);

    auto g_zz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 310);

    auto g_zz_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 311);

    auto g_zz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 312);

    auto g_zz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 313);

    auto g_zz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 314);

    auto g_zz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 315);

    auto g_zz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 316);

    auto g_zz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 317);

    auto g_zz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 318);

    auto g_zz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 319);

    auto g_zz_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 320);

    auto g_zz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 321);

    auto g_zz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 322);

    auto g_zz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 323);

    auto g_zz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 324);

    auto g_zz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 325);

    auto g_zz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 326);

    auto g_zz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 327);

    auto g_zz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 328);

    auto g_zz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 329);

    /// Set up 0-55 components of targeted buffer : FSM

    auto g_xxx_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm);

    auto g_xxx_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 1);

    auto g_xxx_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 2);

    auto g_xxx_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 3);

    auto g_xxx_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 4);

    auto g_xxx_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 5);

    auto g_xxx_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 6);

    auto g_xxx_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 7);

    auto g_xxx_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 8);

    auto g_xxx_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 9);

    auto g_xxx_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 10);

    auto g_xxx_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 11);

    auto g_xxx_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 12);

    auto g_xxx_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 13);

    auto g_xxx_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 14);

    auto g_xxx_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 15);

    auto g_xxx_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 16);

    auto g_xxx_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 17);

    auto g_xxx_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 18);

    auto g_xxx_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 19);

    auto g_xxx_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 20);

    auto g_xxx_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 21);

    auto g_xxx_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 22);

    auto g_xxx_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 23);

    auto g_xxx_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 24);

    auto g_xxx_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 25);

    auto g_xxx_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 26);

    auto g_xxx_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 27);

    auto g_xxx_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 28);

    auto g_xxx_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 29);

    auto g_xxx_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 30);

    auto g_xxx_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 31);

    auto g_xxx_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 32);

    auto g_xxx_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 33);

    auto g_xxx_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 34);

    auto g_xxx_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 35);

    auto g_xxx_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 36);

    auto g_xxx_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 37);

    auto g_xxx_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 38);

    auto g_xxx_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 39);

    auto g_xxx_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 40);

    auto g_xxx_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 41);

    auto g_xxx_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 42);

    auto g_xxx_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 43);

    auto g_xxx_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 44);

    auto g_xxx_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 45);

    auto g_xxx_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 46);

    auto g_xxx_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 47);

    auto g_xxx_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 48);

    auto g_xxx_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 49);

    auto g_xxx_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 50);

    auto g_xxx_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 51);

    auto g_xxx_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 52);

    auto g_xxx_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 53);

    auto g_xxx_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 54);

    #pragma omp simd aligned(g_x_0_xxxxxxxxx_0, g_x_0_xxxxxxxxx_1, g_x_0_xxxxxxxxy_0, g_x_0_xxxxxxxxy_1, g_x_0_xxxxxxxxz_0, g_x_0_xxxxxxxxz_1, g_x_0_xxxxxxxyy_0, g_x_0_xxxxxxxyy_1, g_x_0_xxxxxxxyz_0, g_x_0_xxxxxxxyz_1, g_x_0_xxxxxxxzz_0, g_x_0_xxxxxxxzz_1, g_x_0_xxxxxxyyy_0, g_x_0_xxxxxxyyy_1, g_x_0_xxxxxxyyz_0, g_x_0_xxxxxxyyz_1, g_x_0_xxxxxxyzz_0, g_x_0_xxxxxxyzz_1, g_x_0_xxxxxxzzz_0, g_x_0_xxxxxxzzz_1, g_x_0_xxxxxyyyy_0, g_x_0_xxxxxyyyy_1, g_x_0_xxxxxyyyz_0, g_x_0_xxxxxyyyz_1, g_x_0_xxxxxyyzz_0, g_x_0_xxxxxyyzz_1, g_x_0_xxxxxyzzz_0, g_x_0_xxxxxyzzz_1, g_x_0_xxxxxzzzz_0, g_x_0_xxxxxzzzz_1, g_x_0_xxxxyyyyy_0, g_x_0_xxxxyyyyy_1, g_x_0_xxxxyyyyz_0, g_x_0_xxxxyyyyz_1, g_x_0_xxxxyyyzz_0, g_x_0_xxxxyyyzz_1, g_x_0_xxxxyyzzz_0, g_x_0_xxxxyyzzz_1, g_x_0_xxxxyzzzz_0, g_x_0_xxxxyzzzz_1, g_x_0_xxxxzzzzz_0, g_x_0_xxxxzzzzz_1, g_x_0_xxxyyyyyy_0, g_x_0_xxxyyyyyy_1, g_x_0_xxxyyyyyz_0, g_x_0_xxxyyyyyz_1, g_x_0_xxxyyyyzz_0, g_x_0_xxxyyyyzz_1, g_x_0_xxxyyyzzz_0, g_x_0_xxxyyyzzz_1, g_x_0_xxxyyzzzz_0, g_x_0_xxxyyzzzz_1, g_x_0_xxxyzzzzz_0, g_x_0_xxxyzzzzz_1, g_x_0_xxxzzzzzz_0, g_x_0_xxxzzzzzz_1, g_x_0_xxyyyyyyy_0, g_x_0_xxyyyyyyy_1, g_x_0_xxyyyyyyz_0, g_x_0_xxyyyyyyz_1, g_x_0_xxyyyyyzz_0, g_x_0_xxyyyyyzz_1, g_x_0_xxyyyyzzz_0, g_x_0_xxyyyyzzz_1, g_x_0_xxyyyzzzz_0, g_x_0_xxyyyzzzz_1, g_x_0_xxyyzzzzz_0, g_x_0_xxyyzzzzz_1, g_x_0_xxyzzzzzz_0, g_x_0_xxyzzzzzz_1, g_x_0_xxzzzzzzz_0, g_x_0_xxzzzzzzz_1, g_x_0_xyyyyyyyy_0, g_x_0_xyyyyyyyy_1, g_x_0_xyyyyyyyz_0, g_x_0_xyyyyyyyz_1, g_x_0_xyyyyyyzz_0, g_x_0_xyyyyyyzz_1, g_x_0_xyyyyyzzz_0, g_x_0_xyyyyyzzz_1, g_x_0_xyyyyzzzz_0, g_x_0_xyyyyzzzz_1, g_x_0_xyyyzzzzz_0, g_x_0_xyyyzzzzz_1, g_x_0_xyyzzzzzz_0, g_x_0_xyyzzzzzz_1, g_x_0_xyzzzzzzz_0, g_x_0_xyzzzzzzz_1, g_x_0_xzzzzzzzz_0, g_x_0_xzzzzzzzz_1, g_x_0_yyyyyyyyy_0, g_x_0_yyyyyyyyy_1, g_x_0_yyyyyyyyz_0, g_x_0_yyyyyyyyz_1, g_x_0_yyyyyyyzz_0, g_x_0_yyyyyyyzz_1, g_x_0_yyyyyyzzz_0, g_x_0_yyyyyyzzz_1, g_x_0_yyyyyzzzz_0, g_x_0_yyyyyzzzz_1, g_x_0_yyyyzzzzz_0, g_x_0_yyyyzzzzz_1, g_x_0_yyyzzzzzz_0, g_x_0_yyyzzzzzz_1, g_x_0_yyzzzzzzz_0, g_x_0_yyzzzzzzz_1, g_x_0_yzzzzzzzz_0, g_x_0_yzzzzzzzz_1, g_x_0_zzzzzzzzz_0, g_x_0_zzzzzzzzz_1, g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxxx_1, g_xx_0_xxxxxxxxy_1, g_xx_0_xxxxxxxxz_1, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxxyy_1, g_xx_0_xxxxxxxyz_1, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxxzz_1, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxxyyy_1, g_xx_0_xxxxxxyyz_1, g_xx_0_xxxxxxyz_1, g_xx_0_xxxxxxyzz_1, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxxzzz_1, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxxyyyy_1, g_xx_0_xxxxxyyyz_1, g_xx_0_xxxxxyyz_1, g_xx_0_xxxxxyyzz_1, g_xx_0_xxxxxyzz_1, g_xx_0_xxxxxyzzz_1, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxxzzzz_1, g_xx_0_xxxxyyyy_1, g_xx_0_xxxxyyyyy_1, g_xx_0_xxxxyyyyz_1, g_xx_0_xxxxyyyz_1, g_xx_0_xxxxyyyzz_1, g_xx_0_xxxxyyzz_1, g_xx_0_xxxxyyzzz_1, g_xx_0_xxxxyzzz_1, g_xx_0_xxxxyzzzz_1, g_xx_0_xxxxzzzz_1, g_xx_0_xxxxzzzzz_1, g_xx_0_xxxyyyyy_1, g_xx_0_xxxyyyyyy_1, g_xx_0_xxxyyyyyz_1, g_xx_0_xxxyyyyz_1, g_xx_0_xxxyyyyzz_1, g_xx_0_xxxyyyzz_1, g_xx_0_xxxyyyzzz_1, g_xx_0_xxxyyzzz_1, g_xx_0_xxxyyzzzz_1, g_xx_0_xxxyzzzz_1, g_xx_0_xxxyzzzzz_1, g_xx_0_xxxzzzzz_1, g_xx_0_xxxzzzzzz_1, g_xx_0_xxyyyyyy_1, g_xx_0_xxyyyyyyy_1, g_xx_0_xxyyyyyyz_1, g_xx_0_xxyyyyyz_1, g_xx_0_xxyyyyyzz_1, g_xx_0_xxyyyyzz_1, g_xx_0_xxyyyyzzz_1, g_xx_0_xxyyyzzz_1, g_xx_0_xxyyyzzzz_1, g_xx_0_xxyyzzzz_1, g_xx_0_xxyyzzzzz_1, g_xx_0_xxyzzzzz_1, g_xx_0_xxyzzzzzz_1, g_xx_0_xxzzzzzz_1, g_xx_0_xxzzzzzzz_1, g_xx_0_xyyyyyyy_1, g_xx_0_xyyyyyyyy_1, g_xx_0_xyyyyyyyz_1, g_xx_0_xyyyyyyz_1, g_xx_0_xyyyyyyzz_1, g_xx_0_xyyyyyzz_1, g_xx_0_xyyyyyzzz_1, g_xx_0_xyyyyzzz_1, g_xx_0_xyyyyzzzz_1, g_xx_0_xyyyzzzz_1, g_xx_0_xyyyzzzzz_1, g_xx_0_xyyzzzzz_1, g_xx_0_xyyzzzzzz_1, g_xx_0_xyzzzzzz_1, g_xx_0_xyzzzzzzz_1, g_xx_0_xzzzzzzz_1, g_xx_0_xzzzzzzzz_1, g_xx_0_yyyyyyyy_1, g_xx_0_yyyyyyyyy_1, g_xx_0_yyyyyyyyz_1, g_xx_0_yyyyyyyz_1, g_xx_0_yyyyyyyzz_1, g_xx_0_yyyyyyzz_1, g_xx_0_yyyyyyzzz_1, g_xx_0_yyyyyzzz_1, g_xx_0_yyyyyzzzz_1, g_xx_0_yyyyzzzz_1, g_xx_0_yyyyzzzzz_1, g_xx_0_yyyzzzzz_1, g_xx_0_yyyzzzzzz_1, g_xx_0_yyzzzzzz_1, g_xx_0_yyzzzzzzz_1, g_xx_0_yzzzzzzz_1, g_xx_0_yzzzzzzzz_1, g_xx_0_zzzzzzzz_1, g_xx_0_zzzzzzzzz_1, g_xxx_0_xxxxxxxxx_0, g_xxx_0_xxxxxxxxy_0, g_xxx_0_xxxxxxxxz_0, g_xxx_0_xxxxxxxyy_0, g_xxx_0_xxxxxxxyz_0, g_xxx_0_xxxxxxxzz_0, g_xxx_0_xxxxxxyyy_0, g_xxx_0_xxxxxxyyz_0, g_xxx_0_xxxxxxyzz_0, g_xxx_0_xxxxxxzzz_0, g_xxx_0_xxxxxyyyy_0, g_xxx_0_xxxxxyyyz_0, g_xxx_0_xxxxxyyzz_0, g_xxx_0_xxxxxyzzz_0, g_xxx_0_xxxxxzzzz_0, g_xxx_0_xxxxyyyyy_0, g_xxx_0_xxxxyyyyz_0, g_xxx_0_xxxxyyyzz_0, g_xxx_0_xxxxyyzzz_0, g_xxx_0_xxxxyzzzz_0, g_xxx_0_xxxxzzzzz_0, g_xxx_0_xxxyyyyyy_0, g_xxx_0_xxxyyyyyz_0, g_xxx_0_xxxyyyyzz_0, g_xxx_0_xxxyyyzzz_0, g_xxx_0_xxxyyzzzz_0, g_xxx_0_xxxyzzzzz_0, g_xxx_0_xxxzzzzzz_0, g_xxx_0_xxyyyyyyy_0, g_xxx_0_xxyyyyyyz_0, g_xxx_0_xxyyyyyzz_0, g_xxx_0_xxyyyyzzz_0, g_xxx_0_xxyyyzzzz_0, g_xxx_0_xxyyzzzzz_0, g_xxx_0_xxyzzzzzz_0, g_xxx_0_xxzzzzzzz_0, g_xxx_0_xyyyyyyyy_0, g_xxx_0_xyyyyyyyz_0, g_xxx_0_xyyyyyyzz_0, g_xxx_0_xyyyyyzzz_0, g_xxx_0_xyyyyzzzz_0, g_xxx_0_xyyyzzzzz_0, g_xxx_0_xyyzzzzzz_0, g_xxx_0_xyzzzzzzz_0, g_xxx_0_xzzzzzzzz_0, g_xxx_0_yyyyyyyyy_0, g_xxx_0_yyyyyyyyz_0, g_xxx_0_yyyyyyyzz_0, g_xxx_0_yyyyyyzzz_0, g_xxx_0_yyyyyzzzz_0, g_xxx_0_yyyyzzzzz_0, g_xxx_0_yyyzzzzzz_0, g_xxx_0_yyzzzzzzz_0, g_xxx_0_yzzzzzzzz_0, g_xxx_0_zzzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xxxxxxxxx_0[i] = 2.0 * g_x_0_xxxxxxxxx_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxxx_1[i] * fz_be_0 + 9.0 * g_xx_0_xxxxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxxxx_1[i] * wa_x[i];

        g_xxx_0_xxxxxxxxy_0[i] = 2.0 * g_x_0_xxxxxxxxy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxxy_1[i] * fz_be_0 + 8.0 * g_xx_0_xxxxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xxx_0_xxxxxxxxz_0[i] = 2.0 * g_x_0_xxxxxxxxz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxxz_1[i] * fz_be_0 + 8.0 * g_xx_0_xxxxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xxx_0_xxxxxxxyy_0[i] = 2.0 * g_x_0_xxxxxxxyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxyy_1[i] * fz_be_0 + 7.0 * g_xx_0_xxxxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xxx_0_xxxxxxxyz_0[i] = 2.0 * g_x_0_xxxxxxxyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxyz_1[i] * fz_be_0 + 7.0 * g_xx_0_xxxxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xxx_0_xxxxxxxzz_0[i] = 2.0 * g_x_0_xxxxxxxzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxxzz_1[i] * fz_be_0 + 7.0 * g_xx_0_xxxxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xxx_0_xxxxxxyyy_0[i] = 2.0 * g_x_0_xxxxxxyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxyyy_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xxx_0_xxxxxxyyz_0[i] = 2.0 * g_x_0_xxxxxxyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxyyz_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xxx_0_xxxxxxyzz_0[i] = 2.0 * g_x_0_xxxxxxyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxyzz_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xxx_0_xxxxxxzzz_0[i] = 2.0 * g_x_0_xxxxxxzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxxzzz_1[i] * fz_be_0 + 6.0 * g_xx_0_xxxxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xxx_0_xxxxxyyyy_0[i] = 2.0 * g_x_0_xxxxxyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyyyy_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyyyy_1[i] * fi_acd_0 + g_xx_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xxx_0_xxxxxyyyz_0[i] = 2.0 * g_x_0_xxxxxyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyyyz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyyyz_1[i] * fi_acd_0 + g_xx_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xxx_0_xxxxxyyzz_0[i] = 2.0 * g_x_0_xxxxxyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyyzz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyyzz_1[i] * fi_acd_0 + g_xx_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xxx_0_xxxxxyzzz_0[i] = 2.0 * g_x_0_xxxxxyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxyzzz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxyzzz_1[i] * fi_acd_0 + g_xx_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xxx_0_xxxxxzzzz_0[i] = 2.0 * g_x_0_xxxxxzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxxzzzz_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxxzzzz_1[i] * fi_acd_0 + g_xx_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xxx_0_xxxxyyyyy_0[i] = 2.0 * g_x_0_xxxxyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyyyy_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyyyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xxx_0_xxxxyyyyz_0[i] = 2.0 * g_x_0_xxxxyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyyyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xxx_0_xxxxyyyzz_0[i] = 2.0 * g_x_0_xxxxyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyyzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyyzz_1[i] * fi_acd_0 + g_xx_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xxx_0_xxxxyyzzz_0[i] = 2.0 * g_x_0_xxxxyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyyzzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyyzzz_1[i] * fi_acd_0 + g_xx_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xxx_0_xxxxyzzzz_0[i] = 2.0 * g_x_0_xxxxyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxyzzzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xxx_0_xxxxzzzzz_0[i] = 2.0 * g_x_0_xxxxzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxzzzzz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxzzzzz_1[i] * fi_acd_0 + g_xx_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xxx_0_xxxyyyyyy_0[i] = 2.0 * g_x_0_xxxyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyyyy_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xxx_0_xxxyyyyyz_0[i] = 2.0 * g_x_0_xxxyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyyyz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xxx_0_xxxyyyyzz_0[i] = 2.0 * g_x_0_xxxyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyyzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xxx_0_xxxyyyzzz_0[i] = 2.0 * g_x_0_xxxyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyyzzz_1[i] * fi_acd_0 + g_xx_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xxx_0_xxxyyzzzz_0[i] = 2.0 * g_x_0_xxxyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyyzzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyyzzzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xxx_0_xxxyzzzzz_0[i] = 2.0 * g_x_0_xxxyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyzzzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyzzzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xxx_0_xxxzzzzzz_0[i] = 2.0 * g_x_0_xxxzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxzzzzzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxzzzzzz_1[i] * fi_acd_0 + g_xx_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xxx_0_xxyyyyyyy_0[i] = 2.0 * g_x_0_xxyyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyyyy_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xxx_0_xxyyyyyyz_0[i] = 2.0 * g_x_0_xxyyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyyyz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xxx_0_xxyyyyyzz_0[i] = 2.0 * g_x_0_xxyyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xxx_0_xxyyyyzzz_0[i] = 2.0 * g_x_0_xxyyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyyzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xxx_0_xxyyyzzzz_0[i] = 2.0 * g_x_0_xxyyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyyzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyyzzzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xxx_0_xxyyzzzzz_0[i] = 2.0 * g_x_0_xxyyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyzzzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xxx_0_xxyzzzzzz_0[i] = 2.0 * g_x_0_xxyzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyzzzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyzzzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xxx_0_xxzzzzzzz_0[i] = 2.0 * g_x_0_xxzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxzzzzzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xzzzzzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyyyyyy_0[i] = 2.0 * g_x_0_xyyyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xxx_0_xyyyyyyyz_0[i] = 2.0 * g_x_0_xyyyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xxx_0_xyyyyyyzz_0[i] = 2.0 * g_x_0_xyyyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyyzz_1[i] * fz_be_0 + g_xx_0_yyyyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xxx_0_xyyyyyzzz_0[i] = 2.0 * g_x_0_xyyyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyyzzz_1[i] * fz_be_0 + g_xx_0_yyyyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyyzzzz_0[i] = 2.0 * g_x_0_xyyyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyyzzzz_1[i] * fz_be_0 + g_xx_0_yyyyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyzzzzz_0[i] = 2.0 * g_x_0_xyyyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyzzzzz_1[i] * fz_be_0 + g_xx_0_yyyzzzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xxx_0_xyyzzzzzz_0[i] = 2.0 * g_x_0_xyyzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyzzzzzz_1[i] * fz_be_0 + g_xx_0_yyzzzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xxx_0_xyzzzzzzz_0[i] = 2.0 * g_x_0_xyzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyzzzzzzz_1[i] * fz_be_0 + g_xx_0_yzzzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xxx_0_xzzzzzzzz_0[i] = 2.0 * g_x_0_xzzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xzzzzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyyyyyy_0[i] = 2.0 * g_x_0_yyyyyyyyy_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyyyy_1[i] * fz_be_0 + g_xx_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xxx_0_yyyyyyyyz_0[i] = 2.0 * g_x_0_yyyyyyyyz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyyyz_1[i] * fz_be_0 + g_xx_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xxx_0_yyyyyyyzz_0[i] = 2.0 * g_x_0_yyyyyyyzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyyzz_1[i] * fz_be_0 + g_xx_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xxx_0_yyyyyyzzz_0[i] = 2.0 * g_x_0_yyyyyyzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyyzzz_1[i] * fz_be_0 + g_xx_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyyzzzz_0[i] = 2.0 * g_x_0_yyyyyzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyyzzzz_1[i] * fz_be_0 + g_xx_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyzzzzz_0[i] = 2.0 * g_x_0_yyyyzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyzzzzz_1[i] * fz_be_0 + g_xx_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyzzzzzz_0[i] = 2.0 * g_x_0_yyyzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyzzzzzz_1[i] * fz_be_0 + g_xx_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xxx_0_yyzzzzzzz_0[i] = 2.0 * g_x_0_yyzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyzzzzzzz_1[i] * fz_be_0 + g_xx_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xxx_0_yzzzzzzzz_0[i] = 2.0 * g_x_0_yzzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yzzzzzzzz_1[i] * fz_be_0 + g_xx_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xxx_0_zzzzzzzzz_0[i] = 2.0 * g_x_0_zzzzzzzzz_0[i] * fbe_0 - 2.0 * g_x_0_zzzzzzzzz_1[i] * fz_be_0 + g_xx_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 55-110 components of targeted buffer : FSM

    auto g_xxy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 55);

    auto g_xxy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 56);

    auto g_xxy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 57);

    auto g_xxy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 58);

    auto g_xxy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 59);

    auto g_xxy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 60);

    auto g_xxy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 61);

    auto g_xxy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 62);

    auto g_xxy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 63);

    auto g_xxy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 64);

    auto g_xxy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 65);

    auto g_xxy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 66);

    auto g_xxy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 67);

    auto g_xxy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 68);

    auto g_xxy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 69);

    auto g_xxy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 70);

    auto g_xxy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 71);

    auto g_xxy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 72);

    auto g_xxy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 73);

    auto g_xxy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 74);

    auto g_xxy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 75);

    auto g_xxy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 76);

    auto g_xxy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 77);

    auto g_xxy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 78);

    auto g_xxy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 79);

    auto g_xxy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 80);

    auto g_xxy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 81);

    auto g_xxy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 82);

    auto g_xxy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 83);

    auto g_xxy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 84);

    auto g_xxy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 85);

    auto g_xxy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 86);

    auto g_xxy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 87);

    auto g_xxy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 88);

    auto g_xxy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 89);

    auto g_xxy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 90);

    auto g_xxy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 91);

    auto g_xxy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 92);

    auto g_xxy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 93);

    auto g_xxy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 94);

    auto g_xxy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 95);

    auto g_xxy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 96);

    auto g_xxy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 97);

    auto g_xxy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 98);

    auto g_xxy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 99);

    auto g_xxy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 100);

    auto g_xxy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 101);

    auto g_xxy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 102);

    auto g_xxy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 103);

    auto g_xxy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 104);

    auto g_xxy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 105);

    auto g_xxy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 106);

    auto g_xxy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 107);

    auto g_xxy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 108);

    auto g_xxy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 109);

    #pragma omp simd aligned(g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxxx_1, g_xx_0_xxxxxxxxy_1, g_xx_0_xxxxxxxxz_1, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxxyy_1, g_xx_0_xxxxxxxyz_1, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxxzz_1, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxxyyy_1, g_xx_0_xxxxxxyyz_1, g_xx_0_xxxxxxyz_1, g_xx_0_xxxxxxyzz_1, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxxzzz_1, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxxyyyy_1, g_xx_0_xxxxxyyyz_1, g_xx_0_xxxxxyyz_1, g_xx_0_xxxxxyyzz_1, g_xx_0_xxxxxyzz_1, g_xx_0_xxxxxyzzz_1, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxxzzzz_1, g_xx_0_xxxxyyyy_1, g_xx_0_xxxxyyyyy_1, g_xx_0_xxxxyyyyz_1, g_xx_0_xxxxyyyz_1, g_xx_0_xxxxyyyzz_1, g_xx_0_xxxxyyzz_1, g_xx_0_xxxxyyzzz_1, g_xx_0_xxxxyzzz_1, g_xx_0_xxxxyzzzz_1, g_xx_0_xxxxzzzz_1, g_xx_0_xxxxzzzzz_1, g_xx_0_xxxyyyyy_1, g_xx_0_xxxyyyyyy_1, g_xx_0_xxxyyyyyz_1, g_xx_0_xxxyyyyz_1, g_xx_0_xxxyyyyzz_1, g_xx_0_xxxyyyzz_1, g_xx_0_xxxyyyzzz_1, g_xx_0_xxxyyzzz_1, g_xx_0_xxxyyzzzz_1, g_xx_0_xxxyzzzz_1, g_xx_0_xxxyzzzzz_1, g_xx_0_xxxzzzzz_1, g_xx_0_xxxzzzzzz_1, g_xx_0_xxyyyyyy_1, g_xx_0_xxyyyyyyy_1, g_xx_0_xxyyyyyyz_1, g_xx_0_xxyyyyyz_1, g_xx_0_xxyyyyyzz_1, g_xx_0_xxyyyyzz_1, g_xx_0_xxyyyyzzz_1, g_xx_0_xxyyyzzz_1, g_xx_0_xxyyyzzzz_1, g_xx_0_xxyyzzzz_1, g_xx_0_xxyyzzzzz_1, g_xx_0_xxyzzzzz_1, g_xx_0_xxyzzzzzz_1, g_xx_0_xxzzzzzz_1, g_xx_0_xxzzzzzzz_1, g_xx_0_xyyyyyyy_1, g_xx_0_xyyyyyyyy_1, g_xx_0_xyyyyyyyz_1, g_xx_0_xyyyyyyz_1, g_xx_0_xyyyyyyzz_1, g_xx_0_xyyyyyzz_1, g_xx_0_xyyyyyzzz_1, g_xx_0_xyyyyzzz_1, g_xx_0_xyyyyzzzz_1, g_xx_0_xyyyzzzz_1, g_xx_0_xyyyzzzzz_1, g_xx_0_xyyzzzzz_1, g_xx_0_xyyzzzzzz_1, g_xx_0_xyzzzzzz_1, g_xx_0_xyzzzzzzz_1, g_xx_0_xzzzzzzz_1, g_xx_0_xzzzzzzzz_1, g_xx_0_yyyyyyyy_1, g_xx_0_yyyyyyyyy_1, g_xx_0_yyyyyyyyz_1, g_xx_0_yyyyyyyz_1, g_xx_0_yyyyyyyzz_1, g_xx_0_yyyyyyzz_1, g_xx_0_yyyyyyzzz_1, g_xx_0_yyyyyzzz_1, g_xx_0_yyyyyzzzz_1, g_xx_0_yyyyzzzz_1, g_xx_0_yyyyzzzzz_1, g_xx_0_yyyzzzzz_1, g_xx_0_yyyzzzzzz_1, g_xx_0_yyzzzzzz_1, g_xx_0_yyzzzzzzz_1, g_xx_0_yzzzzzzz_1, g_xx_0_yzzzzzzzz_1, g_xx_0_zzzzzzzz_1, g_xx_0_zzzzzzzzz_1, g_xxy_0_xxxxxxxxx_0, g_xxy_0_xxxxxxxxy_0, g_xxy_0_xxxxxxxxz_0, g_xxy_0_xxxxxxxyy_0, g_xxy_0_xxxxxxxyz_0, g_xxy_0_xxxxxxxzz_0, g_xxy_0_xxxxxxyyy_0, g_xxy_0_xxxxxxyyz_0, g_xxy_0_xxxxxxyzz_0, g_xxy_0_xxxxxxzzz_0, g_xxy_0_xxxxxyyyy_0, g_xxy_0_xxxxxyyyz_0, g_xxy_0_xxxxxyyzz_0, g_xxy_0_xxxxxyzzz_0, g_xxy_0_xxxxxzzzz_0, g_xxy_0_xxxxyyyyy_0, g_xxy_0_xxxxyyyyz_0, g_xxy_0_xxxxyyyzz_0, g_xxy_0_xxxxyyzzz_0, g_xxy_0_xxxxyzzzz_0, g_xxy_0_xxxxzzzzz_0, g_xxy_0_xxxyyyyyy_0, g_xxy_0_xxxyyyyyz_0, g_xxy_0_xxxyyyyzz_0, g_xxy_0_xxxyyyzzz_0, g_xxy_0_xxxyyzzzz_0, g_xxy_0_xxxyzzzzz_0, g_xxy_0_xxxzzzzzz_0, g_xxy_0_xxyyyyyyy_0, g_xxy_0_xxyyyyyyz_0, g_xxy_0_xxyyyyyzz_0, g_xxy_0_xxyyyyzzz_0, g_xxy_0_xxyyyzzzz_0, g_xxy_0_xxyyzzzzz_0, g_xxy_0_xxyzzzzzz_0, g_xxy_0_xxzzzzzzz_0, g_xxy_0_xyyyyyyyy_0, g_xxy_0_xyyyyyyyz_0, g_xxy_0_xyyyyyyzz_0, g_xxy_0_xyyyyyzzz_0, g_xxy_0_xyyyyzzzz_0, g_xxy_0_xyyyzzzzz_0, g_xxy_0_xyyzzzzzz_0, g_xxy_0_xyzzzzzzz_0, g_xxy_0_xzzzzzzzz_0, g_xxy_0_yyyyyyyyy_0, g_xxy_0_yyyyyyyyz_0, g_xxy_0_yyyyyyyzz_0, g_xxy_0_yyyyyyzzz_0, g_xxy_0_yyyyyzzzz_0, g_xxy_0_yyyyzzzzz_0, g_xxy_0_yyyzzzzzz_0, g_xxy_0_yyzzzzzzz_0, g_xxy_0_yzzzzzzzz_0, g_xxy_0_zzzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xxxxxxxxx_0[i] = g_xx_0_xxxxxxxxx_1[i] * wa_y[i];

        g_xxy_0_xxxxxxxxy_0[i] = g_xx_0_xxxxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxxxy_1[i] * wa_y[i];

        g_xxy_0_xxxxxxxxz_0[i] = g_xx_0_xxxxxxxxz_1[i] * wa_y[i];

        g_xxy_0_xxxxxxxyy_0[i] = 2.0 * g_xx_0_xxxxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxxxyy_1[i] * wa_y[i];

        g_xxy_0_xxxxxxxyz_0[i] = g_xx_0_xxxxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxxxyz_1[i] * wa_y[i];

        g_xxy_0_xxxxxxxzz_0[i] = g_xx_0_xxxxxxxzz_1[i] * wa_y[i];

        g_xxy_0_xxxxxxyyy_0[i] = 3.0 * g_xx_0_xxxxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxxxyyy_1[i] * wa_y[i];

        g_xxy_0_xxxxxxyyz_0[i] = 2.0 * g_xx_0_xxxxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxxxyyz_1[i] * wa_y[i];

        g_xxy_0_xxxxxxyzz_0[i] = g_xx_0_xxxxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxxxyzz_1[i] * wa_y[i];

        g_xxy_0_xxxxxxzzz_0[i] = g_xx_0_xxxxxxzzz_1[i] * wa_y[i];

        g_xxy_0_xxxxxyyyy_0[i] = 4.0 * g_xx_0_xxxxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxxxyyyy_1[i] * wa_y[i];

        g_xxy_0_xxxxxyyyz_0[i] = 3.0 * g_xx_0_xxxxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxxxyyyz_1[i] * wa_y[i];

        g_xxy_0_xxxxxyyzz_0[i] = 2.0 * g_xx_0_xxxxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxxxyyzz_1[i] * wa_y[i];

        g_xxy_0_xxxxxyzzz_0[i] = g_xx_0_xxxxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxxxyzzz_1[i] * wa_y[i];

        g_xxy_0_xxxxxzzzz_0[i] = g_xx_0_xxxxxzzzz_1[i] * wa_y[i];

        g_xxy_0_xxxxyyyyy_0[i] = 5.0 * g_xx_0_xxxxyyyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyyyy_1[i] * wa_y[i];

        g_xxy_0_xxxxyyyyz_0[i] = 4.0 * g_xx_0_xxxxyyyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyyyz_1[i] * wa_y[i];

        g_xxy_0_xxxxyyyzz_0[i] = 3.0 * g_xx_0_xxxxyyzz_1[i] * fi_acd_0 + g_xx_0_xxxxyyyzz_1[i] * wa_y[i];

        g_xxy_0_xxxxyyzzz_0[i] = 2.0 * g_xx_0_xxxxyzzz_1[i] * fi_acd_0 + g_xx_0_xxxxyyzzz_1[i] * wa_y[i];

        g_xxy_0_xxxxyzzzz_0[i] = g_xx_0_xxxxzzzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzzzz_1[i] * wa_y[i];

        g_xxy_0_xxxxzzzzz_0[i] = g_xx_0_xxxxzzzzz_1[i] * wa_y[i];

        g_xxy_0_xxxyyyyyy_0[i] = 6.0 * g_xx_0_xxxyyyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyyyy_1[i] * wa_y[i];

        g_xxy_0_xxxyyyyyz_0[i] = 5.0 * g_xx_0_xxxyyyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyyyz_1[i] * wa_y[i];

        g_xxy_0_xxxyyyyzz_0[i] = 4.0 * g_xx_0_xxxyyyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyyyzz_1[i] * wa_y[i];

        g_xxy_0_xxxyyyzzz_0[i] = 3.0 * g_xx_0_xxxyyzzz_1[i] * fi_acd_0 + g_xx_0_xxxyyyzzz_1[i] * wa_y[i];

        g_xxy_0_xxxyyzzzz_0[i] = 2.0 * g_xx_0_xxxyzzzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzzzz_1[i] * wa_y[i];

        g_xxy_0_xxxyzzzzz_0[i] = g_xx_0_xxxzzzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzzzz_1[i] * wa_y[i];

        g_xxy_0_xxxzzzzzz_0[i] = g_xx_0_xxxzzzzzz_1[i] * wa_y[i];

        g_xxy_0_xxyyyyyyy_0[i] = 7.0 * g_xx_0_xxyyyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyyyy_1[i] * wa_y[i];

        g_xxy_0_xxyyyyyyz_0[i] = 6.0 * g_xx_0_xxyyyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyyyz_1[i] * wa_y[i];

        g_xxy_0_xxyyyyyzz_0[i] = 5.0 * g_xx_0_xxyyyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyyyzz_1[i] * wa_y[i];

        g_xxy_0_xxyyyyzzz_0[i] = 4.0 * g_xx_0_xxyyyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyyyzzz_1[i] * wa_y[i];

        g_xxy_0_xxyyyzzzz_0[i] = 3.0 * g_xx_0_xxyyzzzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzzzz_1[i] * wa_y[i];

        g_xxy_0_xxyyzzzzz_0[i] = 2.0 * g_xx_0_xxyzzzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzzzz_1[i] * wa_y[i];

        g_xxy_0_xxyzzzzzz_0[i] = g_xx_0_xxzzzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzzzz_1[i] * wa_y[i];

        g_xxy_0_xxzzzzzzz_0[i] = g_xx_0_xxzzzzzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyyyyyy_0[i] = 8.0 * g_xx_0_xyyyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyyyy_1[i] * wa_y[i];

        g_xxy_0_xyyyyyyyz_0[i] = 7.0 * g_xx_0_xyyyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyyyz_1[i] * wa_y[i];

        g_xxy_0_xyyyyyyzz_0[i] = 6.0 * g_xx_0_xyyyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyyyzz_1[i] * wa_y[i];

        g_xxy_0_xyyyyyzzz_0[i] = 5.0 * g_xx_0_xyyyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyyyzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyyzzzz_0[i] = 4.0 * g_xx_0_xyyyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyzzzzz_0[i] = 3.0 * g_xx_0_xyyzzzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzzzz_1[i] * wa_y[i];

        g_xxy_0_xyyzzzzzz_0[i] = 2.0 * g_xx_0_xyzzzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzzzz_1[i] * wa_y[i];

        g_xxy_0_xyzzzzzzz_0[i] = g_xx_0_xzzzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzzzz_1[i] * wa_y[i];

        g_xxy_0_xzzzzzzzz_0[i] = g_xx_0_xzzzzzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyyyyyy_0[i] = 9.0 * g_xx_0_yyyyyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyyyyy_1[i] * wa_y[i];

        g_xxy_0_yyyyyyyyz_0[i] = 8.0 * g_xx_0_yyyyyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyyyyyz_1[i] * wa_y[i];

        g_xxy_0_yyyyyyyzz_0[i] = 7.0 * g_xx_0_yyyyyyzz_1[i] * fi_acd_0 + g_xx_0_yyyyyyyzz_1[i] * wa_y[i];

        g_xxy_0_yyyyyyzzz_0[i] = 6.0 * g_xx_0_yyyyyzzz_1[i] * fi_acd_0 + g_xx_0_yyyyyyzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyyzzzz_0[i] = 5.0 * g_xx_0_yyyyzzzz_1[i] * fi_acd_0 + g_xx_0_yyyyyzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyzzzzz_0[i] = 4.0 * g_xx_0_yyyzzzzz_1[i] * fi_acd_0 + g_xx_0_yyyyzzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyzzzzzz_0[i] = 3.0 * g_xx_0_yyzzzzzz_1[i] * fi_acd_0 + g_xx_0_yyyzzzzzz_1[i] * wa_y[i];

        g_xxy_0_yyzzzzzzz_0[i] = 2.0 * g_xx_0_yzzzzzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzzzzz_1[i] * wa_y[i];

        g_xxy_0_yzzzzzzzz_0[i] = g_xx_0_zzzzzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzzzzz_1[i] * wa_y[i];

        g_xxy_0_zzzzzzzzz_0[i] = g_xx_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 110-165 components of targeted buffer : FSM

    auto g_xxz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 110);

    auto g_xxz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 111);

    auto g_xxz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 112);

    auto g_xxz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 113);

    auto g_xxz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 114);

    auto g_xxz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 115);

    auto g_xxz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 116);

    auto g_xxz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 117);

    auto g_xxz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 118);

    auto g_xxz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 119);

    auto g_xxz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 120);

    auto g_xxz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 121);

    auto g_xxz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 122);

    auto g_xxz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 123);

    auto g_xxz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 124);

    auto g_xxz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 125);

    auto g_xxz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 126);

    auto g_xxz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 127);

    auto g_xxz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 128);

    auto g_xxz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 129);

    auto g_xxz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 130);

    auto g_xxz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 131);

    auto g_xxz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 132);

    auto g_xxz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 133);

    auto g_xxz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 134);

    auto g_xxz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 135);

    auto g_xxz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 136);

    auto g_xxz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 137);

    auto g_xxz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 138);

    auto g_xxz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 139);

    auto g_xxz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 140);

    auto g_xxz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 141);

    auto g_xxz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 142);

    auto g_xxz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 143);

    auto g_xxz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 144);

    auto g_xxz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 145);

    auto g_xxz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 146);

    auto g_xxz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 147);

    auto g_xxz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 148);

    auto g_xxz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 149);

    auto g_xxz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 150);

    auto g_xxz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 151);

    auto g_xxz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 152);

    auto g_xxz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 153);

    auto g_xxz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 154);

    auto g_xxz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 155);

    auto g_xxz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 156);

    auto g_xxz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 157);

    auto g_xxz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 158);

    auto g_xxz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 159);

    auto g_xxz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 160);

    auto g_xxz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 161);

    auto g_xxz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 162);

    auto g_xxz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 163);

    auto g_xxz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 164);

    #pragma omp simd aligned(g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxxx_1, g_xx_0_xxxxxxxxy_1, g_xx_0_xxxxxxxxz_1, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxxyy_1, g_xx_0_xxxxxxxyz_1, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxxzz_1, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxxyyy_1, g_xx_0_xxxxxxyyz_1, g_xx_0_xxxxxxyz_1, g_xx_0_xxxxxxyzz_1, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxxzzz_1, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxxyyyy_1, g_xx_0_xxxxxyyyz_1, g_xx_0_xxxxxyyz_1, g_xx_0_xxxxxyyzz_1, g_xx_0_xxxxxyzz_1, g_xx_0_xxxxxyzzz_1, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxxzzzz_1, g_xx_0_xxxxyyyy_1, g_xx_0_xxxxyyyyy_1, g_xx_0_xxxxyyyyz_1, g_xx_0_xxxxyyyz_1, g_xx_0_xxxxyyyzz_1, g_xx_0_xxxxyyzz_1, g_xx_0_xxxxyyzzz_1, g_xx_0_xxxxyzzz_1, g_xx_0_xxxxyzzzz_1, g_xx_0_xxxxzzzz_1, g_xx_0_xxxxzzzzz_1, g_xx_0_xxxyyyyy_1, g_xx_0_xxxyyyyyy_1, g_xx_0_xxxyyyyyz_1, g_xx_0_xxxyyyyz_1, g_xx_0_xxxyyyyzz_1, g_xx_0_xxxyyyzz_1, g_xx_0_xxxyyyzzz_1, g_xx_0_xxxyyzzz_1, g_xx_0_xxxyyzzzz_1, g_xx_0_xxxyzzzz_1, g_xx_0_xxxyzzzzz_1, g_xx_0_xxxzzzzz_1, g_xx_0_xxxzzzzzz_1, g_xx_0_xxyyyyyy_1, g_xx_0_xxyyyyyyy_1, g_xx_0_xxyyyyyyz_1, g_xx_0_xxyyyyyz_1, g_xx_0_xxyyyyyzz_1, g_xx_0_xxyyyyzz_1, g_xx_0_xxyyyyzzz_1, g_xx_0_xxyyyzzz_1, g_xx_0_xxyyyzzzz_1, g_xx_0_xxyyzzzz_1, g_xx_0_xxyyzzzzz_1, g_xx_0_xxyzzzzz_1, g_xx_0_xxyzzzzzz_1, g_xx_0_xxzzzzzz_1, g_xx_0_xxzzzzzzz_1, g_xx_0_xyyyyyyy_1, g_xx_0_xyyyyyyyy_1, g_xx_0_xyyyyyyyz_1, g_xx_0_xyyyyyyz_1, g_xx_0_xyyyyyyzz_1, g_xx_0_xyyyyyzz_1, g_xx_0_xyyyyyzzz_1, g_xx_0_xyyyyzzz_1, g_xx_0_xyyyyzzzz_1, g_xx_0_xyyyzzzz_1, g_xx_0_xyyyzzzzz_1, g_xx_0_xyyzzzzz_1, g_xx_0_xyyzzzzzz_1, g_xx_0_xyzzzzzz_1, g_xx_0_xyzzzzzzz_1, g_xx_0_xzzzzzzz_1, g_xx_0_xzzzzzzzz_1, g_xx_0_yyyyyyyy_1, g_xx_0_yyyyyyyyy_1, g_xx_0_yyyyyyyyz_1, g_xx_0_yyyyyyyz_1, g_xx_0_yyyyyyyzz_1, g_xx_0_yyyyyyzz_1, g_xx_0_yyyyyyzzz_1, g_xx_0_yyyyyzzz_1, g_xx_0_yyyyyzzzz_1, g_xx_0_yyyyzzzz_1, g_xx_0_yyyyzzzzz_1, g_xx_0_yyyzzzzz_1, g_xx_0_yyyzzzzzz_1, g_xx_0_yyzzzzzz_1, g_xx_0_yyzzzzzzz_1, g_xx_0_yzzzzzzz_1, g_xx_0_yzzzzzzzz_1, g_xx_0_zzzzzzzz_1, g_xx_0_zzzzzzzzz_1, g_xxz_0_xxxxxxxxx_0, g_xxz_0_xxxxxxxxy_0, g_xxz_0_xxxxxxxxz_0, g_xxz_0_xxxxxxxyy_0, g_xxz_0_xxxxxxxyz_0, g_xxz_0_xxxxxxxzz_0, g_xxz_0_xxxxxxyyy_0, g_xxz_0_xxxxxxyyz_0, g_xxz_0_xxxxxxyzz_0, g_xxz_0_xxxxxxzzz_0, g_xxz_0_xxxxxyyyy_0, g_xxz_0_xxxxxyyyz_0, g_xxz_0_xxxxxyyzz_0, g_xxz_0_xxxxxyzzz_0, g_xxz_0_xxxxxzzzz_0, g_xxz_0_xxxxyyyyy_0, g_xxz_0_xxxxyyyyz_0, g_xxz_0_xxxxyyyzz_0, g_xxz_0_xxxxyyzzz_0, g_xxz_0_xxxxyzzzz_0, g_xxz_0_xxxxzzzzz_0, g_xxz_0_xxxyyyyyy_0, g_xxz_0_xxxyyyyyz_0, g_xxz_0_xxxyyyyzz_0, g_xxz_0_xxxyyyzzz_0, g_xxz_0_xxxyyzzzz_0, g_xxz_0_xxxyzzzzz_0, g_xxz_0_xxxzzzzzz_0, g_xxz_0_xxyyyyyyy_0, g_xxz_0_xxyyyyyyz_0, g_xxz_0_xxyyyyyzz_0, g_xxz_0_xxyyyyzzz_0, g_xxz_0_xxyyyzzzz_0, g_xxz_0_xxyyzzzzz_0, g_xxz_0_xxyzzzzzz_0, g_xxz_0_xxzzzzzzz_0, g_xxz_0_xyyyyyyyy_0, g_xxz_0_xyyyyyyyz_0, g_xxz_0_xyyyyyyzz_0, g_xxz_0_xyyyyyzzz_0, g_xxz_0_xyyyyzzzz_0, g_xxz_0_xyyyzzzzz_0, g_xxz_0_xyyzzzzzz_0, g_xxz_0_xyzzzzzzz_0, g_xxz_0_xzzzzzzzz_0, g_xxz_0_yyyyyyyyy_0, g_xxz_0_yyyyyyyyz_0, g_xxz_0_yyyyyyyzz_0, g_xxz_0_yyyyyyzzz_0, g_xxz_0_yyyyyzzzz_0, g_xxz_0_yyyyzzzzz_0, g_xxz_0_yyyzzzzzz_0, g_xxz_0_yyzzzzzzz_0, g_xxz_0_yzzzzzzzz_0, g_xxz_0_zzzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xxxxxxxxx_0[i] = g_xx_0_xxxxxxxxx_1[i] * wa_z[i];

        g_xxz_0_xxxxxxxxy_0[i] = g_xx_0_xxxxxxxxy_1[i] * wa_z[i];

        g_xxz_0_xxxxxxxxz_0[i] = g_xx_0_xxxxxxxx_1[i] * fi_acd_0 + g_xx_0_xxxxxxxxz_1[i] * wa_z[i];

        g_xxz_0_xxxxxxxyy_0[i] = g_xx_0_xxxxxxxyy_1[i] * wa_z[i];

        g_xxz_0_xxxxxxxyz_0[i] = g_xx_0_xxxxxxxy_1[i] * fi_acd_0 + g_xx_0_xxxxxxxyz_1[i] * wa_z[i];

        g_xxz_0_xxxxxxxzz_0[i] = 2.0 * g_xx_0_xxxxxxxz_1[i] * fi_acd_0 + g_xx_0_xxxxxxxzz_1[i] * wa_z[i];

        g_xxz_0_xxxxxxyyy_0[i] = g_xx_0_xxxxxxyyy_1[i] * wa_z[i];

        g_xxz_0_xxxxxxyyz_0[i] = g_xx_0_xxxxxxyy_1[i] * fi_acd_0 + g_xx_0_xxxxxxyyz_1[i] * wa_z[i];

        g_xxz_0_xxxxxxyzz_0[i] = 2.0 * g_xx_0_xxxxxxyz_1[i] * fi_acd_0 + g_xx_0_xxxxxxyzz_1[i] * wa_z[i];

        g_xxz_0_xxxxxxzzz_0[i] = 3.0 * g_xx_0_xxxxxxzz_1[i] * fi_acd_0 + g_xx_0_xxxxxxzzz_1[i] * wa_z[i];

        g_xxz_0_xxxxxyyyy_0[i] = g_xx_0_xxxxxyyyy_1[i] * wa_z[i];

        g_xxz_0_xxxxxyyyz_0[i] = g_xx_0_xxxxxyyy_1[i] * fi_acd_0 + g_xx_0_xxxxxyyyz_1[i] * wa_z[i];

        g_xxz_0_xxxxxyyzz_0[i] = 2.0 * g_xx_0_xxxxxyyz_1[i] * fi_acd_0 + g_xx_0_xxxxxyyzz_1[i] * wa_z[i];

        g_xxz_0_xxxxxyzzz_0[i] = 3.0 * g_xx_0_xxxxxyzz_1[i] * fi_acd_0 + g_xx_0_xxxxxyzzz_1[i] * wa_z[i];

        g_xxz_0_xxxxxzzzz_0[i] = 4.0 * g_xx_0_xxxxxzzz_1[i] * fi_acd_0 + g_xx_0_xxxxxzzzz_1[i] * wa_z[i];

        g_xxz_0_xxxxyyyyy_0[i] = g_xx_0_xxxxyyyyy_1[i] * wa_z[i];

        g_xxz_0_xxxxyyyyz_0[i] = g_xx_0_xxxxyyyy_1[i] * fi_acd_0 + g_xx_0_xxxxyyyyz_1[i] * wa_z[i];

        g_xxz_0_xxxxyyyzz_0[i] = 2.0 * g_xx_0_xxxxyyyz_1[i] * fi_acd_0 + g_xx_0_xxxxyyyzz_1[i] * wa_z[i];

        g_xxz_0_xxxxyyzzz_0[i] = 3.0 * g_xx_0_xxxxyyzz_1[i] * fi_acd_0 + g_xx_0_xxxxyyzzz_1[i] * wa_z[i];

        g_xxz_0_xxxxyzzzz_0[i] = 4.0 * g_xx_0_xxxxyzzz_1[i] * fi_acd_0 + g_xx_0_xxxxyzzzz_1[i] * wa_z[i];

        g_xxz_0_xxxxzzzzz_0[i] = 5.0 * g_xx_0_xxxxzzzz_1[i] * fi_acd_0 + g_xx_0_xxxxzzzzz_1[i] * wa_z[i];

        g_xxz_0_xxxyyyyyy_0[i] = g_xx_0_xxxyyyyyy_1[i] * wa_z[i];

        g_xxz_0_xxxyyyyyz_0[i] = g_xx_0_xxxyyyyy_1[i] * fi_acd_0 + g_xx_0_xxxyyyyyz_1[i] * wa_z[i];

        g_xxz_0_xxxyyyyzz_0[i] = 2.0 * g_xx_0_xxxyyyyz_1[i] * fi_acd_0 + g_xx_0_xxxyyyyzz_1[i] * wa_z[i];

        g_xxz_0_xxxyyyzzz_0[i] = 3.0 * g_xx_0_xxxyyyzz_1[i] * fi_acd_0 + g_xx_0_xxxyyyzzz_1[i] * wa_z[i];

        g_xxz_0_xxxyyzzzz_0[i] = 4.0 * g_xx_0_xxxyyzzz_1[i] * fi_acd_0 + g_xx_0_xxxyyzzzz_1[i] * wa_z[i];

        g_xxz_0_xxxyzzzzz_0[i] = 5.0 * g_xx_0_xxxyzzzz_1[i] * fi_acd_0 + g_xx_0_xxxyzzzzz_1[i] * wa_z[i];

        g_xxz_0_xxxzzzzzz_0[i] = 6.0 * g_xx_0_xxxzzzzz_1[i] * fi_acd_0 + g_xx_0_xxxzzzzzz_1[i] * wa_z[i];

        g_xxz_0_xxyyyyyyy_0[i] = g_xx_0_xxyyyyyyy_1[i] * wa_z[i];

        g_xxz_0_xxyyyyyyz_0[i] = g_xx_0_xxyyyyyy_1[i] * fi_acd_0 + g_xx_0_xxyyyyyyz_1[i] * wa_z[i];

        g_xxz_0_xxyyyyyzz_0[i] = 2.0 * g_xx_0_xxyyyyyz_1[i] * fi_acd_0 + g_xx_0_xxyyyyyzz_1[i] * wa_z[i];

        g_xxz_0_xxyyyyzzz_0[i] = 3.0 * g_xx_0_xxyyyyzz_1[i] * fi_acd_0 + g_xx_0_xxyyyyzzz_1[i] * wa_z[i];

        g_xxz_0_xxyyyzzzz_0[i] = 4.0 * g_xx_0_xxyyyzzz_1[i] * fi_acd_0 + g_xx_0_xxyyyzzzz_1[i] * wa_z[i];

        g_xxz_0_xxyyzzzzz_0[i] = 5.0 * g_xx_0_xxyyzzzz_1[i] * fi_acd_0 + g_xx_0_xxyyzzzzz_1[i] * wa_z[i];

        g_xxz_0_xxyzzzzzz_0[i] = 6.0 * g_xx_0_xxyzzzzz_1[i] * fi_acd_0 + g_xx_0_xxyzzzzzz_1[i] * wa_z[i];

        g_xxz_0_xxzzzzzzz_0[i] = 7.0 * g_xx_0_xxzzzzzz_1[i] * fi_acd_0 + g_xx_0_xxzzzzzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyyyyyy_0[i] = g_xx_0_xyyyyyyyy_1[i] * wa_z[i];

        g_xxz_0_xyyyyyyyz_0[i] = g_xx_0_xyyyyyyy_1[i] * fi_acd_0 + g_xx_0_xyyyyyyyz_1[i] * wa_z[i];

        g_xxz_0_xyyyyyyzz_0[i] = 2.0 * g_xx_0_xyyyyyyz_1[i] * fi_acd_0 + g_xx_0_xyyyyyyzz_1[i] * wa_z[i];

        g_xxz_0_xyyyyyzzz_0[i] = 3.0 * g_xx_0_xyyyyyzz_1[i] * fi_acd_0 + g_xx_0_xyyyyyzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyyzzzz_0[i] = 4.0 * g_xx_0_xyyyyzzz_1[i] * fi_acd_0 + g_xx_0_xyyyyzzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyzzzzz_0[i] = 5.0 * g_xx_0_xyyyzzzz_1[i] * fi_acd_0 + g_xx_0_xyyyzzzzz_1[i] * wa_z[i];

        g_xxz_0_xyyzzzzzz_0[i] = 6.0 * g_xx_0_xyyzzzzz_1[i] * fi_acd_0 + g_xx_0_xyyzzzzzz_1[i] * wa_z[i];

        g_xxz_0_xyzzzzzzz_0[i] = 7.0 * g_xx_0_xyzzzzzz_1[i] * fi_acd_0 + g_xx_0_xyzzzzzzz_1[i] * wa_z[i];

        g_xxz_0_xzzzzzzzz_0[i] = 8.0 * g_xx_0_xzzzzzzz_1[i] * fi_acd_0 + g_xx_0_xzzzzzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyyyyyy_0[i] = g_xx_0_yyyyyyyyy_1[i] * wa_z[i];

        g_xxz_0_yyyyyyyyz_0[i] = g_xx_0_yyyyyyyy_1[i] * fi_acd_0 + g_xx_0_yyyyyyyyz_1[i] * wa_z[i];

        g_xxz_0_yyyyyyyzz_0[i] = 2.0 * g_xx_0_yyyyyyyz_1[i] * fi_acd_0 + g_xx_0_yyyyyyyzz_1[i] * wa_z[i];

        g_xxz_0_yyyyyyzzz_0[i] = 3.0 * g_xx_0_yyyyyyzz_1[i] * fi_acd_0 + g_xx_0_yyyyyyzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyyzzzz_0[i] = 4.0 * g_xx_0_yyyyyzzz_1[i] * fi_acd_0 + g_xx_0_yyyyyzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyzzzzz_0[i] = 5.0 * g_xx_0_yyyyzzzz_1[i] * fi_acd_0 + g_xx_0_yyyyzzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyzzzzzz_0[i] = 6.0 * g_xx_0_yyyzzzzz_1[i] * fi_acd_0 + g_xx_0_yyyzzzzzz_1[i] * wa_z[i];

        g_xxz_0_yyzzzzzzz_0[i] = 7.0 * g_xx_0_yyzzzzzz_1[i] * fi_acd_0 + g_xx_0_yyzzzzzzz_1[i] * wa_z[i];

        g_xxz_0_yzzzzzzzz_0[i] = 8.0 * g_xx_0_yzzzzzzz_1[i] * fi_acd_0 + g_xx_0_yzzzzzzzz_1[i] * wa_z[i];

        g_xxz_0_zzzzzzzzz_0[i] = 9.0 * g_xx_0_zzzzzzzz_1[i] * fi_acd_0 + g_xx_0_zzzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 165-220 components of targeted buffer : FSM

    auto g_xyy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 165);

    auto g_xyy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 166);

    auto g_xyy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 167);

    auto g_xyy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 168);

    auto g_xyy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 169);

    auto g_xyy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 170);

    auto g_xyy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 171);

    auto g_xyy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 172);

    auto g_xyy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 173);

    auto g_xyy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 174);

    auto g_xyy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 175);

    auto g_xyy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 176);

    auto g_xyy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 177);

    auto g_xyy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 178);

    auto g_xyy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 179);

    auto g_xyy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 180);

    auto g_xyy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 181);

    auto g_xyy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 182);

    auto g_xyy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 183);

    auto g_xyy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 184);

    auto g_xyy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 185);

    auto g_xyy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 186);

    auto g_xyy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 187);

    auto g_xyy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 188);

    auto g_xyy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 189);

    auto g_xyy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 190);

    auto g_xyy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 191);

    auto g_xyy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 192);

    auto g_xyy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 193);

    auto g_xyy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 194);

    auto g_xyy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 195);

    auto g_xyy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 196);

    auto g_xyy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 197);

    auto g_xyy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 198);

    auto g_xyy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 199);

    auto g_xyy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 200);

    auto g_xyy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 201);

    auto g_xyy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 202);

    auto g_xyy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 203);

    auto g_xyy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 204);

    auto g_xyy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 205);

    auto g_xyy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 206);

    auto g_xyy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 207);

    auto g_xyy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 208);

    auto g_xyy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 209);

    auto g_xyy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 210);

    auto g_xyy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 211);

    auto g_xyy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 212);

    auto g_xyy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 213);

    auto g_xyy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 214);

    auto g_xyy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 215);

    auto g_xyy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 216);

    auto g_xyy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 217);

    auto g_xyy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 218);

    auto g_xyy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 219);

    #pragma omp simd aligned(g_xyy_0_xxxxxxxxx_0, g_xyy_0_xxxxxxxxy_0, g_xyy_0_xxxxxxxxz_0, g_xyy_0_xxxxxxxyy_0, g_xyy_0_xxxxxxxyz_0, g_xyy_0_xxxxxxxzz_0, g_xyy_0_xxxxxxyyy_0, g_xyy_0_xxxxxxyyz_0, g_xyy_0_xxxxxxyzz_0, g_xyy_0_xxxxxxzzz_0, g_xyy_0_xxxxxyyyy_0, g_xyy_0_xxxxxyyyz_0, g_xyy_0_xxxxxyyzz_0, g_xyy_0_xxxxxyzzz_0, g_xyy_0_xxxxxzzzz_0, g_xyy_0_xxxxyyyyy_0, g_xyy_0_xxxxyyyyz_0, g_xyy_0_xxxxyyyzz_0, g_xyy_0_xxxxyyzzz_0, g_xyy_0_xxxxyzzzz_0, g_xyy_0_xxxxzzzzz_0, g_xyy_0_xxxyyyyyy_0, g_xyy_0_xxxyyyyyz_0, g_xyy_0_xxxyyyyzz_0, g_xyy_0_xxxyyyzzz_0, g_xyy_0_xxxyyzzzz_0, g_xyy_0_xxxyzzzzz_0, g_xyy_0_xxxzzzzzz_0, g_xyy_0_xxyyyyyyy_0, g_xyy_0_xxyyyyyyz_0, g_xyy_0_xxyyyyyzz_0, g_xyy_0_xxyyyyzzz_0, g_xyy_0_xxyyyzzzz_0, g_xyy_0_xxyyzzzzz_0, g_xyy_0_xxyzzzzzz_0, g_xyy_0_xxzzzzzzz_0, g_xyy_0_xyyyyyyyy_0, g_xyy_0_xyyyyyyyz_0, g_xyy_0_xyyyyyyzz_0, g_xyy_0_xyyyyyzzz_0, g_xyy_0_xyyyyzzzz_0, g_xyy_0_xyyyzzzzz_0, g_xyy_0_xyyzzzzzz_0, g_xyy_0_xyzzzzzzz_0, g_xyy_0_xzzzzzzzz_0, g_xyy_0_yyyyyyyyy_0, g_xyy_0_yyyyyyyyz_0, g_xyy_0_yyyyyyyzz_0, g_xyy_0_yyyyyyzzz_0, g_xyy_0_yyyyyzzzz_0, g_xyy_0_yyyyzzzzz_0, g_xyy_0_yyyzzzzzz_0, g_xyy_0_yyzzzzzzz_0, g_xyy_0_yzzzzzzzz_0, g_xyy_0_zzzzzzzzz_0, g_yy_0_xxxxxxxx_1, g_yy_0_xxxxxxxxx_1, g_yy_0_xxxxxxxxy_1, g_yy_0_xxxxxxxxz_1, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxxyy_1, g_yy_0_xxxxxxxyz_1, g_yy_0_xxxxxxxz_1, g_yy_0_xxxxxxxzz_1, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyyy_1, g_yy_0_xxxxxxyyz_1, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxxyzz_1, g_yy_0_xxxxxxzz_1, g_yy_0_xxxxxxzzz_1, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyyy_1, g_yy_0_xxxxxyyyz_1, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyyzz_1, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxxyzzz_1, g_yy_0_xxxxxzzz_1, g_yy_0_xxxxxzzzz_1, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyyy_1, g_yy_0_xxxxyyyyz_1, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyyzz_1, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyyzzz_1, g_yy_0_xxxxyzzz_1, g_yy_0_xxxxyzzzz_1, g_yy_0_xxxxzzzz_1, g_yy_0_xxxxzzzzz_1, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyyy_1, g_yy_0_xxxyyyyyz_1, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyyzz_1, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyyzzz_1, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyyzzzz_1, g_yy_0_xxxyzzzz_1, g_yy_0_xxxyzzzzz_1, g_yy_0_xxxzzzzz_1, g_yy_0_xxxzzzzzz_1, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyyy_1, g_yy_0_xxyyyyyyz_1, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyyzz_1, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyyzzz_1, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyyzzzz_1, g_yy_0_xxyyzzzz_1, g_yy_0_xxyyzzzzz_1, g_yy_0_xxyzzzzz_1, g_yy_0_xxyzzzzzz_1, g_yy_0_xxzzzzzz_1, g_yy_0_xxzzzzzzz_1, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyyy_1, g_yy_0_xyyyyyyyz_1, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyyzz_1, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyyzzz_1, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyyzzzz_1, g_yy_0_xyyyzzzz_1, g_yy_0_xyyyzzzzz_1, g_yy_0_xyyzzzzz_1, g_yy_0_xyyzzzzzz_1, g_yy_0_xyzzzzzz_1, g_yy_0_xyzzzzzzz_1, g_yy_0_xzzzzzzz_1, g_yy_0_xzzzzzzzz_1, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyyy_1, g_yy_0_yyyyyyyyz_1, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyyzz_1, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyyzzz_1, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyyzzzz_1, g_yy_0_yyyyzzzz_1, g_yy_0_yyyyzzzzz_1, g_yy_0_yyyzzzzz_1, g_yy_0_yyyzzzzzz_1, g_yy_0_yyzzzzzz_1, g_yy_0_yyzzzzzzz_1, g_yy_0_yzzzzzzz_1, g_yy_0_yzzzzzzzz_1, g_yy_0_zzzzzzzz_1, g_yy_0_zzzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xxxxxxxxx_0[i] = 9.0 * g_yy_0_xxxxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxxxx_1[i] * wa_x[i];

        g_xyy_0_xxxxxxxxy_0[i] = 8.0 * g_yy_0_xxxxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xyy_0_xxxxxxxxz_0[i] = 8.0 * g_yy_0_xxxxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xyy_0_xxxxxxxyy_0[i] = 7.0 * g_yy_0_xxxxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xyy_0_xxxxxxxyz_0[i] = 7.0 * g_yy_0_xxxxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xyy_0_xxxxxxxzz_0[i] = 7.0 * g_yy_0_xxxxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xyy_0_xxxxxxyyy_0[i] = 6.0 * g_yy_0_xxxxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xyy_0_xxxxxxyyz_0[i] = 6.0 * g_yy_0_xxxxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xyy_0_xxxxxxyzz_0[i] = 6.0 * g_yy_0_xxxxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xyy_0_xxxxxxzzz_0[i] = 6.0 * g_yy_0_xxxxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xyy_0_xxxxxyyyy_0[i] = 5.0 * g_yy_0_xxxxyyyy_1[i] * fi_acd_0 + g_yy_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xyy_0_xxxxxyyyz_0[i] = 5.0 * g_yy_0_xxxxyyyz_1[i] * fi_acd_0 + g_yy_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xyy_0_xxxxxyyzz_0[i] = 5.0 * g_yy_0_xxxxyyzz_1[i] * fi_acd_0 + g_yy_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xyy_0_xxxxxyzzz_0[i] = 5.0 * g_yy_0_xxxxyzzz_1[i] * fi_acd_0 + g_yy_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xyy_0_xxxxxzzzz_0[i] = 5.0 * g_yy_0_xxxxzzzz_1[i] * fi_acd_0 + g_yy_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xyy_0_xxxxyyyyy_0[i] = 4.0 * g_yy_0_xxxyyyyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xyy_0_xxxxyyyyz_0[i] = 4.0 * g_yy_0_xxxyyyyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xyy_0_xxxxyyyzz_0[i] = 4.0 * g_yy_0_xxxyyyzz_1[i] * fi_acd_0 + g_yy_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xyy_0_xxxxyyzzz_0[i] = 4.0 * g_yy_0_xxxyyzzz_1[i] * fi_acd_0 + g_yy_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xyy_0_xxxxyzzzz_0[i] = 4.0 * g_yy_0_xxxyzzzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xyy_0_xxxxzzzzz_0[i] = 4.0 * g_yy_0_xxxzzzzz_1[i] * fi_acd_0 + g_yy_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xyy_0_xxxyyyyyy_0[i] = 3.0 * g_yy_0_xxyyyyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xyy_0_xxxyyyyyz_0[i] = 3.0 * g_yy_0_xxyyyyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xyy_0_xxxyyyyzz_0[i] = 3.0 * g_yy_0_xxyyyyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xyy_0_xxxyyyzzz_0[i] = 3.0 * g_yy_0_xxyyyzzz_1[i] * fi_acd_0 + g_yy_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xyy_0_xxxyyzzzz_0[i] = 3.0 * g_yy_0_xxyyzzzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xyy_0_xxxyzzzzz_0[i] = 3.0 * g_yy_0_xxyzzzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xyy_0_xxxzzzzzz_0[i] = 3.0 * g_yy_0_xxzzzzzz_1[i] * fi_acd_0 + g_yy_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xyy_0_xxyyyyyyy_0[i] = 2.0 * g_yy_0_xyyyyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xyy_0_xxyyyyyyz_0[i] = 2.0 * g_yy_0_xyyyyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xyy_0_xxyyyyyzz_0[i] = 2.0 * g_yy_0_xyyyyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xyy_0_xxyyyyzzz_0[i] = 2.0 * g_yy_0_xyyyyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xyy_0_xxyyyzzzz_0[i] = 2.0 * g_yy_0_xyyyzzzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xyy_0_xxyyzzzzz_0[i] = 2.0 * g_yy_0_xyyzzzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xyy_0_xxyzzzzzz_0[i] = 2.0 * g_yy_0_xyzzzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xyy_0_xxzzzzzzz_0[i] = 2.0 * g_yy_0_xzzzzzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyyyyyy_0[i] = g_yy_0_yyyyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xyy_0_xyyyyyyyz_0[i] = g_yy_0_yyyyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xyy_0_xyyyyyyzz_0[i] = g_yy_0_yyyyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xyy_0_xyyyyyzzz_0[i] = g_yy_0_yyyyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyyzzzz_0[i] = g_yy_0_yyyyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyzzzzz_0[i] = g_yy_0_yyyzzzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xyy_0_xyyzzzzzz_0[i] = g_yy_0_yyzzzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xyy_0_xyzzzzzzz_0[i] = g_yy_0_yzzzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xyy_0_xzzzzzzzz_0[i] = g_yy_0_zzzzzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyyyyyy_0[i] = g_yy_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xyy_0_yyyyyyyyz_0[i] = g_yy_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xyy_0_yyyyyyyzz_0[i] = g_yy_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xyy_0_yyyyyyzzz_0[i] = g_yy_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyyzzzz_0[i] = g_yy_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyzzzzz_0[i] = g_yy_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyzzzzzz_0[i] = g_yy_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xyy_0_yyzzzzzzz_0[i] = g_yy_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xyy_0_yzzzzzzzz_0[i] = g_yy_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xyy_0_zzzzzzzzz_0[i] = g_yy_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 220-275 components of targeted buffer : FSM

    auto g_xyz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 220);

    auto g_xyz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 221);

    auto g_xyz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 222);

    auto g_xyz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 223);

    auto g_xyz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 224);

    auto g_xyz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 225);

    auto g_xyz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 226);

    auto g_xyz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 227);

    auto g_xyz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 228);

    auto g_xyz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 229);

    auto g_xyz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 230);

    auto g_xyz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 231);

    auto g_xyz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 232);

    auto g_xyz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 233);

    auto g_xyz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 234);

    auto g_xyz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 235);

    auto g_xyz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 236);

    auto g_xyz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 237);

    auto g_xyz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 238);

    auto g_xyz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 239);

    auto g_xyz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 240);

    auto g_xyz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 241);

    auto g_xyz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 242);

    auto g_xyz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 243);

    auto g_xyz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 244);

    auto g_xyz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 245);

    auto g_xyz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 246);

    auto g_xyz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 247);

    auto g_xyz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 248);

    auto g_xyz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 249);

    auto g_xyz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 250);

    auto g_xyz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 251);

    auto g_xyz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 252);

    auto g_xyz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 253);

    auto g_xyz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 254);

    auto g_xyz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 255);

    auto g_xyz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 256);

    auto g_xyz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 257);

    auto g_xyz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 258);

    auto g_xyz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 259);

    auto g_xyz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 260);

    auto g_xyz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 261);

    auto g_xyz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 262);

    auto g_xyz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 263);

    auto g_xyz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 264);

    auto g_xyz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 265);

    auto g_xyz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 266);

    auto g_xyz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 267);

    auto g_xyz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 268);

    auto g_xyz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 269);

    auto g_xyz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 270);

    auto g_xyz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 271);

    auto g_xyz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 272);

    auto g_xyz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 273);

    auto g_xyz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 274);

    #pragma omp simd aligned(g_xy_0_xxxxxxxxy_1, g_xy_0_xxxxxxxyy_1, g_xy_0_xxxxxxyyy_1, g_xy_0_xxxxxyyyy_1, g_xy_0_xxxxyyyyy_1, g_xy_0_xxxyyyyyy_1, g_xy_0_xxyyyyyyy_1, g_xy_0_xyyyyyyyy_1, g_xyz_0_xxxxxxxxx_0, g_xyz_0_xxxxxxxxy_0, g_xyz_0_xxxxxxxxz_0, g_xyz_0_xxxxxxxyy_0, g_xyz_0_xxxxxxxyz_0, g_xyz_0_xxxxxxxzz_0, g_xyz_0_xxxxxxyyy_0, g_xyz_0_xxxxxxyyz_0, g_xyz_0_xxxxxxyzz_0, g_xyz_0_xxxxxxzzz_0, g_xyz_0_xxxxxyyyy_0, g_xyz_0_xxxxxyyyz_0, g_xyz_0_xxxxxyyzz_0, g_xyz_0_xxxxxyzzz_0, g_xyz_0_xxxxxzzzz_0, g_xyz_0_xxxxyyyyy_0, g_xyz_0_xxxxyyyyz_0, g_xyz_0_xxxxyyyzz_0, g_xyz_0_xxxxyyzzz_0, g_xyz_0_xxxxyzzzz_0, g_xyz_0_xxxxzzzzz_0, g_xyz_0_xxxyyyyyy_0, g_xyz_0_xxxyyyyyz_0, g_xyz_0_xxxyyyyzz_0, g_xyz_0_xxxyyyzzz_0, g_xyz_0_xxxyyzzzz_0, g_xyz_0_xxxyzzzzz_0, g_xyz_0_xxxzzzzzz_0, g_xyz_0_xxyyyyyyy_0, g_xyz_0_xxyyyyyyz_0, g_xyz_0_xxyyyyyzz_0, g_xyz_0_xxyyyyzzz_0, g_xyz_0_xxyyyzzzz_0, g_xyz_0_xxyyzzzzz_0, g_xyz_0_xxyzzzzzz_0, g_xyz_0_xxzzzzzzz_0, g_xyz_0_xyyyyyyyy_0, g_xyz_0_xyyyyyyyz_0, g_xyz_0_xyyyyyyzz_0, g_xyz_0_xyyyyyzzz_0, g_xyz_0_xyyyyzzzz_0, g_xyz_0_xyyyzzzzz_0, g_xyz_0_xyyzzzzzz_0, g_xyz_0_xyzzzzzzz_0, g_xyz_0_xzzzzzzzz_0, g_xyz_0_yyyyyyyyy_0, g_xyz_0_yyyyyyyyz_0, g_xyz_0_yyyyyyyzz_0, g_xyz_0_yyyyyyzzz_0, g_xyz_0_yyyyyzzzz_0, g_xyz_0_yyyyzzzzz_0, g_xyz_0_yyyzzzzzz_0, g_xyz_0_yyzzzzzzz_0, g_xyz_0_yzzzzzzzz_0, g_xyz_0_zzzzzzzzz_0, g_xz_0_xxxxxxxxx_1, g_xz_0_xxxxxxxxz_1, g_xz_0_xxxxxxxzz_1, g_xz_0_xxxxxxzzz_1, g_xz_0_xxxxxzzzz_1, g_xz_0_xxxxzzzzz_1, g_xz_0_xxxzzzzzz_1, g_xz_0_xxzzzzzzz_1, g_xz_0_xzzzzzzzz_1, g_yz_0_xxxxxxxyz_1, g_yz_0_xxxxxxyyz_1, g_yz_0_xxxxxxyz_1, g_yz_0_xxxxxxyzz_1, g_yz_0_xxxxxyyyz_1, g_yz_0_xxxxxyyz_1, g_yz_0_xxxxxyyzz_1, g_yz_0_xxxxxyzz_1, g_yz_0_xxxxxyzzz_1, g_yz_0_xxxxyyyyz_1, g_yz_0_xxxxyyyz_1, g_yz_0_xxxxyyyzz_1, g_yz_0_xxxxyyzz_1, g_yz_0_xxxxyyzzz_1, g_yz_0_xxxxyzzz_1, g_yz_0_xxxxyzzzz_1, g_yz_0_xxxyyyyyz_1, g_yz_0_xxxyyyyz_1, g_yz_0_xxxyyyyzz_1, g_yz_0_xxxyyyzz_1, g_yz_0_xxxyyyzzz_1, g_yz_0_xxxyyzzz_1, g_yz_0_xxxyyzzzz_1, g_yz_0_xxxyzzzz_1, g_yz_0_xxxyzzzzz_1, g_yz_0_xxyyyyyyz_1, g_yz_0_xxyyyyyz_1, g_yz_0_xxyyyyyzz_1, g_yz_0_xxyyyyzz_1, g_yz_0_xxyyyyzzz_1, g_yz_0_xxyyyzzz_1, g_yz_0_xxyyyzzzz_1, g_yz_0_xxyyzzzz_1, g_yz_0_xxyyzzzzz_1, g_yz_0_xxyzzzzz_1, g_yz_0_xxyzzzzzz_1, g_yz_0_xyyyyyyyz_1, g_yz_0_xyyyyyyz_1, g_yz_0_xyyyyyyzz_1, g_yz_0_xyyyyyzz_1, g_yz_0_xyyyyyzzz_1, g_yz_0_xyyyyzzz_1, g_yz_0_xyyyyzzzz_1, g_yz_0_xyyyzzzz_1, g_yz_0_xyyyzzzzz_1, g_yz_0_xyyzzzzz_1, g_yz_0_xyyzzzzzz_1, g_yz_0_xyzzzzzz_1, g_yz_0_xyzzzzzzz_1, g_yz_0_yyyyyyyyy_1, g_yz_0_yyyyyyyyz_1, g_yz_0_yyyyyyyz_1, g_yz_0_yyyyyyyzz_1, g_yz_0_yyyyyyzz_1, g_yz_0_yyyyyyzzz_1, g_yz_0_yyyyyzzz_1, g_yz_0_yyyyyzzzz_1, g_yz_0_yyyyzzzz_1, g_yz_0_yyyyzzzzz_1, g_yz_0_yyyzzzzz_1, g_yz_0_yyyzzzzzz_1, g_yz_0_yyzzzzzz_1, g_yz_0_yyzzzzzzz_1, g_yz_0_yzzzzzzz_1, g_yz_0_yzzzzzzzz_1, g_yz_0_zzzzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyz_0_xxxxxxxxx_0[i] = g_xz_0_xxxxxxxxx_1[i] * wa_y[i];

        g_xyz_0_xxxxxxxxy_0[i] = g_xy_0_xxxxxxxxy_1[i] * wa_z[i];

        g_xyz_0_xxxxxxxxz_0[i] = g_xz_0_xxxxxxxxz_1[i] * wa_y[i];

        g_xyz_0_xxxxxxxyy_0[i] = g_xy_0_xxxxxxxyy_1[i] * wa_z[i];

        g_xyz_0_xxxxxxxyz_0[i] = 7.0 * g_yz_0_xxxxxxyz_1[i] * fi_acd_0 + g_yz_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xyz_0_xxxxxxxzz_0[i] = g_xz_0_xxxxxxxzz_1[i] * wa_y[i];

        g_xyz_0_xxxxxxyyy_0[i] = g_xy_0_xxxxxxyyy_1[i] * wa_z[i];

        g_xyz_0_xxxxxxyyz_0[i] = 6.0 * g_yz_0_xxxxxyyz_1[i] * fi_acd_0 + g_yz_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xyz_0_xxxxxxyzz_0[i] = 6.0 * g_yz_0_xxxxxyzz_1[i] * fi_acd_0 + g_yz_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xyz_0_xxxxxxzzz_0[i] = g_xz_0_xxxxxxzzz_1[i] * wa_y[i];

        g_xyz_0_xxxxxyyyy_0[i] = g_xy_0_xxxxxyyyy_1[i] * wa_z[i];

        g_xyz_0_xxxxxyyyz_0[i] = 5.0 * g_yz_0_xxxxyyyz_1[i] * fi_acd_0 + g_yz_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xyz_0_xxxxxyyzz_0[i] = 5.0 * g_yz_0_xxxxyyzz_1[i] * fi_acd_0 + g_yz_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xyz_0_xxxxxyzzz_0[i] = 5.0 * g_yz_0_xxxxyzzz_1[i] * fi_acd_0 + g_yz_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xyz_0_xxxxxzzzz_0[i] = g_xz_0_xxxxxzzzz_1[i] * wa_y[i];

        g_xyz_0_xxxxyyyyy_0[i] = g_xy_0_xxxxyyyyy_1[i] * wa_z[i];

        g_xyz_0_xxxxyyyyz_0[i] = 4.0 * g_yz_0_xxxyyyyz_1[i] * fi_acd_0 + g_yz_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xyz_0_xxxxyyyzz_0[i] = 4.0 * g_yz_0_xxxyyyzz_1[i] * fi_acd_0 + g_yz_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xyz_0_xxxxyyzzz_0[i] = 4.0 * g_yz_0_xxxyyzzz_1[i] * fi_acd_0 + g_yz_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xyz_0_xxxxyzzzz_0[i] = 4.0 * g_yz_0_xxxyzzzz_1[i] * fi_acd_0 + g_yz_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xyz_0_xxxxzzzzz_0[i] = g_xz_0_xxxxzzzzz_1[i] * wa_y[i];

        g_xyz_0_xxxyyyyyy_0[i] = g_xy_0_xxxyyyyyy_1[i] * wa_z[i];

        g_xyz_0_xxxyyyyyz_0[i] = 3.0 * g_yz_0_xxyyyyyz_1[i] * fi_acd_0 + g_yz_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xyz_0_xxxyyyyzz_0[i] = 3.0 * g_yz_0_xxyyyyzz_1[i] * fi_acd_0 + g_yz_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xyz_0_xxxyyyzzz_0[i] = 3.0 * g_yz_0_xxyyyzzz_1[i] * fi_acd_0 + g_yz_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xyz_0_xxxyyzzzz_0[i] = 3.0 * g_yz_0_xxyyzzzz_1[i] * fi_acd_0 + g_yz_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xyz_0_xxxyzzzzz_0[i] = 3.0 * g_yz_0_xxyzzzzz_1[i] * fi_acd_0 + g_yz_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xyz_0_xxxzzzzzz_0[i] = g_xz_0_xxxzzzzzz_1[i] * wa_y[i];

        g_xyz_0_xxyyyyyyy_0[i] = g_xy_0_xxyyyyyyy_1[i] * wa_z[i];

        g_xyz_0_xxyyyyyyz_0[i] = 2.0 * g_yz_0_xyyyyyyz_1[i] * fi_acd_0 + g_yz_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xyz_0_xxyyyyyzz_0[i] = 2.0 * g_yz_0_xyyyyyzz_1[i] * fi_acd_0 + g_yz_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xyz_0_xxyyyyzzz_0[i] = 2.0 * g_yz_0_xyyyyzzz_1[i] * fi_acd_0 + g_yz_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xyz_0_xxyyyzzzz_0[i] = 2.0 * g_yz_0_xyyyzzzz_1[i] * fi_acd_0 + g_yz_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xyz_0_xxyyzzzzz_0[i] = 2.0 * g_yz_0_xyyzzzzz_1[i] * fi_acd_0 + g_yz_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xyz_0_xxyzzzzzz_0[i] = 2.0 * g_yz_0_xyzzzzzz_1[i] * fi_acd_0 + g_yz_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xyz_0_xxzzzzzzz_0[i] = g_xz_0_xxzzzzzzz_1[i] * wa_y[i];

        g_xyz_0_xyyyyyyyy_0[i] = g_xy_0_xyyyyyyyy_1[i] * wa_z[i];

        g_xyz_0_xyyyyyyyz_0[i] = g_yz_0_yyyyyyyz_1[i] * fi_acd_0 + g_yz_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xyz_0_xyyyyyyzz_0[i] = g_yz_0_yyyyyyzz_1[i] * fi_acd_0 + g_yz_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xyz_0_xyyyyyzzz_0[i] = g_yz_0_yyyyyzzz_1[i] * fi_acd_0 + g_yz_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xyz_0_xyyyyzzzz_0[i] = g_yz_0_yyyyzzzz_1[i] * fi_acd_0 + g_yz_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xyz_0_xyyyzzzzz_0[i] = g_yz_0_yyyzzzzz_1[i] * fi_acd_0 + g_yz_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xyz_0_xyyzzzzzz_0[i] = g_yz_0_yyzzzzzz_1[i] * fi_acd_0 + g_yz_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xyz_0_xyzzzzzzz_0[i] = g_yz_0_yzzzzzzz_1[i] * fi_acd_0 + g_yz_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xyz_0_xzzzzzzzz_0[i] = g_xz_0_xzzzzzzzz_1[i] * wa_y[i];

        g_xyz_0_yyyyyyyyy_0[i] = g_yz_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xyz_0_yyyyyyyyz_0[i] = g_yz_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xyz_0_yyyyyyyzz_0[i] = g_yz_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xyz_0_yyyyyyzzz_0[i] = g_yz_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xyz_0_yyyyyzzzz_0[i] = g_yz_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xyz_0_yyyyzzzzz_0[i] = g_yz_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xyz_0_yyyzzzzzz_0[i] = g_yz_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xyz_0_yyzzzzzzz_0[i] = g_yz_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xyz_0_yzzzzzzzz_0[i] = g_yz_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xyz_0_zzzzzzzzz_0[i] = g_yz_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 275-330 components of targeted buffer : FSM

    auto g_xzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 275);

    auto g_xzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 276);

    auto g_xzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 277);

    auto g_xzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 278);

    auto g_xzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 279);

    auto g_xzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 280);

    auto g_xzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 281);

    auto g_xzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 282);

    auto g_xzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 283);

    auto g_xzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 284);

    auto g_xzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 285);

    auto g_xzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 286);

    auto g_xzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 287);

    auto g_xzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 288);

    auto g_xzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 289);

    auto g_xzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 290);

    auto g_xzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 291);

    auto g_xzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 292);

    auto g_xzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 293);

    auto g_xzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 294);

    auto g_xzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 295);

    auto g_xzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 296);

    auto g_xzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 297);

    auto g_xzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 298);

    auto g_xzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 299);

    auto g_xzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 300);

    auto g_xzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 301);

    auto g_xzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 302);

    auto g_xzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 303);

    auto g_xzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 304);

    auto g_xzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 305);

    auto g_xzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 306);

    auto g_xzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 307);

    auto g_xzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 308);

    auto g_xzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 309);

    auto g_xzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 310);

    auto g_xzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 311);

    auto g_xzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 312);

    auto g_xzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 313);

    auto g_xzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 314);

    auto g_xzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 315);

    auto g_xzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 316);

    auto g_xzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 317);

    auto g_xzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 318);

    auto g_xzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 319);

    auto g_xzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 320);

    auto g_xzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 321);

    auto g_xzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 322);

    auto g_xzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 323);

    auto g_xzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 324);

    auto g_xzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 325);

    auto g_xzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 326);

    auto g_xzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 327);

    auto g_xzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 328);

    auto g_xzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 329);

    #pragma omp simd aligned(g_xzz_0_xxxxxxxxx_0, g_xzz_0_xxxxxxxxy_0, g_xzz_0_xxxxxxxxz_0, g_xzz_0_xxxxxxxyy_0, g_xzz_0_xxxxxxxyz_0, g_xzz_0_xxxxxxxzz_0, g_xzz_0_xxxxxxyyy_0, g_xzz_0_xxxxxxyyz_0, g_xzz_0_xxxxxxyzz_0, g_xzz_0_xxxxxxzzz_0, g_xzz_0_xxxxxyyyy_0, g_xzz_0_xxxxxyyyz_0, g_xzz_0_xxxxxyyzz_0, g_xzz_0_xxxxxyzzz_0, g_xzz_0_xxxxxzzzz_0, g_xzz_0_xxxxyyyyy_0, g_xzz_0_xxxxyyyyz_0, g_xzz_0_xxxxyyyzz_0, g_xzz_0_xxxxyyzzz_0, g_xzz_0_xxxxyzzzz_0, g_xzz_0_xxxxzzzzz_0, g_xzz_0_xxxyyyyyy_0, g_xzz_0_xxxyyyyyz_0, g_xzz_0_xxxyyyyzz_0, g_xzz_0_xxxyyyzzz_0, g_xzz_0_xxxyyzzzz_0, g_xzz_0_xxxyzzzzz_0, g_xzz_0_xxxzzzzzz_0, g_xzz_0_xxyyyyyyy_0, g_xzz_0_xxyyyyyyz_0, g_xzz_0_xxyyyyyzz_0, g_xzz_0_xxyyyyzzz_0, g_xzz_0_xxyyyzzzz_0, g_xzz_0_xxyyzzzzz_0, g_xzz_0_xxyzzzzzz_0, g_xzz_0_xxzzzzzzz_0, g_xzz_0_xyyyyyyyy_0, g_xzz_0_xyyyyyyyz_0, g_xzz_0_xyyyyyyzz_0, g_xzz_0_xyyyyyzzz_0, g_xzz_0_xyyyyzzzz_0, g_xzz_0_xyyyzzzzz_0, g_xzz_0_xyyzzzzzz_0, g_xzz_0_xyzzzzzzz_0, g_xzz_0_xzzzzzzzz_0, g_xzz_0_yyyyyyyyy_0, g_xzz_0_yyyyyyyyz_0, g_xzz_0_yyyyyyyzz_0, g_xzz_0_yyyyyyzzz_0, g_xzz_0_yyyyyzzzz_0, g_xzz_0_yyyyzzzzz_0, g_xzz_0_yyyzzzzzz_0, g_xzz_0_yyzzzzzzz_0, g_xzz_0_yzzzzzzzz_0, g_xzz_0_zzzzzzzzz_0, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxxx_1, g_zz_0_xxxxxxxxy_1, g_zz_0_xxxxxxxxz_1, g_zz_0_xxxxxxxy_1, g_zz_0_xxxxxxxyy_1, g_zz_0_xxxxxxxyz_1, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxxzz_1, g_zz_0_xxxxxxyy_1, g_zz_0_xxxxxxyyy_1, g_zz_0_xxxxxxyyz_1, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxyzz_1, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxxzzz_1, g_zz_0_xxxxxyyy_1, g_zz_0_xxxxxyyyy_1, g_zz_0_xxxxxyyyz_1, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyyzz_1, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxyzzz_1, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxxzzzz_1, g_zz_0_xxxxyyyy_1, g_zz_0_xxxxyyyyy_1, g_zz_0_xxxxyyyyz_1, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyyzz_1, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyyzzz_1, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxyzzzz_1, g_zz_0_xxxxzzzz_1, g_zz_0_xxxxzzzzz_1, g_zz_0_xxxyyyyy_1, g_zz_0_xxxyyyyyy_1, g_zz_0_xxxyyyyyz_1, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyyzz_1, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyyzzz_1, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyyzzzz_1, g_zz_0_xxxyzzzz_1, g_zz_0_xxxyzzzzz_1, g_zz_0_xxxzzzzz_1, g_zz_0_xxxzzzzzz_1, g_zz_0_xxyyyyyy_1, g_zz_0_xxyyyyyyy_1, g_zz_0_xxyyyyyyz_1, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyyzz_1, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyyzzz_1, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyyzzzz_1, g_zz_0_xxyyzzzz_1, g_zz_0_xxyyzzzzz_1, g_zz_0_xxyzzzzz_1, g_zz_0_xxyzzzzzz_1, g_zz_0_xxzzzzzz_1, g_zz_0_xxzzzzzzz_1, g_zz_0_xyyyyyyy_1, g_zz_0_xyyyyyyyy_1, g_zz_0_xyyyyyyyz_1, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyyzz_1, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyyzzz_1, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyyzzzz_1, g_zz_0_xyyyzzzz_1, g_zz_0_xyyyzzzzz_1, g_zz_0_xyyzzzzz_1, g_zz_0_xyyzzzzzz_1, g_zz_0_xyzzzzzz_1, g_zz_0_xyzzzzzzz_1, g_zz_0_xzzzzzzz_1, g_zz_0_xzzzzzzzz_1, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyyy_1, g_zz_0_yyyyyyyyz_1, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyyzz_1, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyyzzz_1, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyyzzzz_1, g_zz_0_yyyyzzzz_1, g_zz_0_yyyyzzzzz_1, g_zz_0_yyyzzzzz_1, g_zz_0_yyyzzzzzz_1, g_zz_0_yyzzzzzz_1, g_zz_0_yyzzzzzzz_1, g_zz_0_yzzzzzzz_1, g_zz_0_yzzzzzzzz_1, g_zz_0_zzzzzzzz_1, g_zz_0_zzzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xxxxxxxxx_0[i] = 9.0 * g_zz_0_xxxxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxxxx_1[i] * wa_x[i];

        g_xzz_0_xxxxxxxxy_0[i] = 8.0 * g_zz_0_xxxxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xzz_0_xxxxxxxxz_0[i] = 8.0 * g_zz_0_xxxxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xzz_0_xxxxxxxyy_0[i] = 7.0 * g_zz_0_xxxxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xzz_0_xxxxxxxyz_0[i] = 7.0 * g_zz_0_xxxxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xzz_0_xxxxxxxzz_0[i] = 7.0 * g_zz_0_xxxxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xzz_0_xxxxxxyyy_0[i] = 6.0 * g_zz_0_xxxxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xzz_0_xxxxxxyyz_0[i] = 6.0 * g_zz_0_xxxxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xzz_0_xxxxxxyzz_0[i] = 6.0 * g_zz_0_xxxxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xzz_0_xxxxxxzzz_0[i] = 6.0 * g_zz_0_xxxxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xzz_0_xxxxxyyyy_0[i] = 5.0 * g_zz_0_xxxxyyyy_1[i] * fi_acd_0 + g_zz_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xzz_0_xxxxxyyyz_0[i] = 5.0 * g_zz_0_xxxxyyyz_1[i] * fi_acd_0 + g_zz_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xzz_0_xxxxxyyzz_0[i] = 5.0 * g_zz_0_xxxxyyzz_1[i] * fi_acd_0 + g_zz_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xzz_0_xxxxxyzzz_0[i] = 5.0 * g_zz_0_xxxxyzzz_1[i] * fi_acd_0 + g_zz_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xzz_0_xxxxxzzzz_0[i] = 5.0 * g_zz_0_xxxxzzzz_1[i] * fi_acd_0 + g_zz_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xzz_0_xxxxyyyyy_0[i] = 4.0 * g_zz_0_xxxyyyyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xzz_0_xxxxyyyyz_0[i] = 4.0 * g_zz_0_xxxyyyyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xzz_0_xxxxyyyzz_0[i] = 4.0 * g_zz_0_xxxyyyzz_1[i] * fi_acd_0 + g_zz_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xzz_0_xxxxyyzzz_0[i] = 4.0 * g_zz_0_xxxyyzzz_1[i] * fi_acd_0 + g_zz_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xzz_0_xxxxyzzzz_0[i] = 4.0 * g_zz_0_xxxyzzzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xzz_0_xxxxzzzzz_0[i] = 4.0 * g_zz_0_xxxzzzzz_1[i] * fi_acd_0 + g_zz_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xzz_0_xxxyyyyyy_0[i] = 3.0 * g_zz_0_xxyyyyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xzz_0_xxxyyyyyz_0[i] = 3.0 * g_zz_0_xxyyyyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xzz_0_xxxyyyyzz_0[i] = 3.0 * g_zz_0_xxyyyyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xzz_0_xxxyyyzzz_0[i] = 3.0 * g_zz_0_xxyyyzzz_1[i] * fi_acd_0 + g_zz_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xzz_0_xxxyyzzzz_0[i] = 3.0 * g_zz_0_xxyyzzzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xzz_0_xxxyzzzzz_0[i] = 3.0 * g_zz_0_xxyzzzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xzz_0_xxxzzzzzz_0[i] = 3.0 * g_zz_0_xxzzzzzz_1[i] * fi_acd_0 + g_zz_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xzz_0_xxyyyyyyy_0[i] = 2.0 * g_zz_0_xyyyyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xzz_0_xxyyyyyyz_0[i] = 2.0 * g_zz_0_xyyyyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xzz_0_xxyyyyyzz_0[i] = 2.0 * g_zz_0_xyyyyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xzz_0_xxyyyyzzz_0[i] = 2.0 * g_zz_0_xyyyyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xzz_0_xxyyyzzzz_0[i] = 2.0 * g_zz_0_xyyyzzzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xzz_0_xxyyzzzzz_0[i] = 2.0 * g_zz_0_xyyzzzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xzz_0_xxyzzzzzz_0[i] = 2.0 * g_zz_0_xyzzzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xzz_0_xxzzzzzzz_0[i] = 2.0 * g_zz_0_xzzzzzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyyyyyy_0[i] = g_zz_0_yyyyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xzz_0_xyyyyyyyz_0[i] = g_zz_0_yyyyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xzz_0_xyyyyyyzz_0[i] = g_zz_0_yyyyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xzz_0_xyyyyyzzz_0[i] = g_zz_0_yyyyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyyzzzz_0[i] = g_zz_0_yyyyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyzzzzz_0[i] = g_zz_0_yyyzzzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xzz_0_xyyzzzzzz_0[i] = g_zz_0_yyzzzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xzz_0_xyzzzzzzz_0[i] = g_zz_0_yzzzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xzz_0_xzzzzzzzz_0[i] = g_zz_0_zzzzzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyyyyyy_0[i] = g_zz_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xzz_0_yyyyyyyyz_0[i] = g_zz_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xzz_0_yyyyyyyzz_0[i] = g_zz_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xzz_0_yyyyyyzzz_0[i] = g_zz_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyyzzzz_0[i] = g_zz_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyzzzzz_0[i] = g_zz_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyzzzzzz_0[i] = g_zz_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xzz_0_yyzzzzzzz_0[i] = g_zz_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xzz_0_yzzzzzzzz_0[i] = g_zz_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xzz_0_zzzzzzzzz_0[i] = g_zz_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 330-385 components of targeted buffer : FSM

    auto g_yyy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 330);

    auto g_yyy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 331);

    auto g_yyy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 332);

    auto g_yyy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 333);

    auto g_yyy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 334);

    auto g_yyy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 335);

    auto g_yyy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 336);

    auto g_yyy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 337);

    auto g_yyy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 338);

    auto g_yyy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 339);

    auto g_yyy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 340);

    auto g_yyy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 341);

    auto g_yyy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 342);

    auto g_yyy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 343);

    auto g_yyy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 344);

    auto g_yyy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 345);

    auto g_yyy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 346);

    auto g_yyy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 347);

    auto g_yyy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 348);

    auto g_yyy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 349);

    auto g_yyy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 350);

    auto g_yyy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 351);

    auto g_yyy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 352);

    auto g_yyy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 353);

    auto g_yyy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 354);

    auto g_yyy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 355);

    auto g_yyy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 356);

    auto g_yyy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 357);

    auto g_yyy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 358);

    auto g_yyy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 359);

    auto g_yyy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 360);

    auto g_yyy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 361);

    auto g_yyy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 362);

    auto g_yyy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 363);

    auto g_yyy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 364);

    auto g_yyy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 365);

    auto g_yyy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 366);

    auto g_yyy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 367);

    auto g_yyy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 368);

    auto g_yyy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 369);

    auto g_yyy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 370);

    auto g_yyy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 371);

    auto g_yyy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 372);

    auto g_yyy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 373);

    auto g_yyy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 374);

    auto g_yyy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 375);

    auto g_yyy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 376);

    auto g_yyy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 377);

    auto g_yyy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 378);

    auto g_yyy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 379);

    auto g_yyy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 380);

    auto g_yyy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 381);

    auto g_yyy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 382);

    auto g_yyy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 383);

    auto g_yyy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 384);

    #pragma omp simd aligned(g_y_0_xxxxxxxxx_0, g_y_0_xxxxxxxxx_1, g_y_0_xxxxxxxxy_0, g_y_0_xxxxxxxxy_1, g_y_0_xxxxxxxxz_0, g_y_0_xxxxxxxxz_1, g_y_0_xxxxxxxyy_0, g_y_0_xxxxxxxyy_1, g_y_0_xxxxxxxyz_0, g_y_0_xxxxxxxyz_1, g_y_0_xxxxxxxzz_0, g_y_0_xxxxxxxzz_1, g_y_0_xxxxxxyyy_0, g_y_0_xxxxxxyyy_1, g_y_0_xxxxxxyyz_0, g_y_0_xxxxxxyyz_1, g_y_0_xxxxxxyzz_0, g_y_0_xxxxxxyzz_1, g_y_0_xxxxxxzzz_0, g_y_0_xxxxxxzzz_1, g_y_0_xxxxxyyyy_0, g_y_0_xxxxxyyyy_1, g_y_0_xxxxxyyyz_0, g_y_0_xxxxxyyyz_1, g_y_0_xxxxxyyzz_0, g_y_0_xxxxxyyzz_1, g_y_0_xxxxxyzzz_0, g_y_0_xxxxxyzzz_1, g_y_0_xxxxxzzzz_0, g_y_0_xxxxxzzzz_1, g_y_0_xxxxyyyyy_0, g_y_0_xxxxyyyyy_1, g_y_0_xxxxyyyyz_0, g_y_0_xxxxyyyyz_1, g_y_0_xxxxyyyzz_0, g_y_0_xxxxyyyzz_1, g_y_0_xxxxyyzzz_0, g_y_0_xxxxyyzzz_1, g_y_0_xxxxyzzzz_0, g_y_0_xxxxyzzzz_1, g_y_0_xxxxzzzzz_0, g_y_0_xxxxzzzzz_1, g_y_0_xxxyyyyyy_0, g_y_0_xxxyyyyyy_1, g_y_0_xxxyyyyyz_0, g_y_0_xxxyyyyyz_1, g_y_0_xxxyyyyzz_0, g_y_0_xxxyyyyzz_1, g_y_0_xxxyyyzzz_0, g_y_0_xxxyyyzzz_1, g_y_0_xxxyyzzzz_0, g_y_0_xxxyyzzzz_1, g_y_0_xxxyzzzzz_0, g_y_0_xxxyzzzzz_1, g_y_0_xxxzzzzzz_0, g_y_0_xxxzzzzzz_1, g_y_0_xxyyyyyyy_0, g_y_0_xxyyyyyyy_1, g_y_0_xxyyyyyyz_0, g_y_0_xxyyyyyyz_1, g_y_0_xxyyyyyzz_0, g_y_0_xxyyyyyzz_1, g_y_0_xxyyyyzzz_0, g_y_0_xxyyyyzzz_1, g_y_0_xxyyyzzzz_0, g_y_0_xxyyyzzzz_1, g_y_0_xxyyzzzzz_0, g_y_0_xxyyzzzzz_1, g_y_0_xxyzzzzzz_0, g_y_0_xxyzzzzzz_1, g_y_0_xxzzzzzzz_0, g_y_0_xxzzzzzzz_1, g_y_0_xyyyyyyyy_0, g_y_0_xyyyyyyyy_1, g_y_0_xyyyyyyyz_0, g_y_0_xyyyyyyyz_1, g_y_0_xyyyyyyzz_0, g_y_0_xyyyyyyzz_1, g_y_0_xyyyyyzzz_0, g_y_0_xyyyyyzzz_1, g_y_0_xyyyyzzzz_0, g_y_0_xyyyyzzzz_1, g_y_0_xyyyzzzzz_0, g_y_0_xyyyzzzzz_1, g_y_0_xyyzzzzzz_0, g_y_0_xyyzzzzzz_1, g_y_0_xyzzzzzzz_0, g_y_0_xyzzzzzzz_1, g_y_0_xzzzzzzzz_0, g_y_0_xzzzzzzzz_1, g_y_0_yyyyyyyyy_0, g_y_0_yyyyyyyyy_1, g_y_0_yyyyyyyyz_0, g_y_0_yyyyyyyyz_1, g_y_0_yyyyyyyzz_0, g_y_0_yyyyyyyzz_1, g_y_0_yyyyyyzzz_0, g_y_0_yyyyyyzzz_1, g_y_0_yyyyyzzzz_0, g_y_0_yyyyyzzzz_1, g_y_0_yyyyzzzzz_0, g_y_0_yyyyzzzzz_1, g_y_0_yyyzzzzzz_0, g_y_0_yyyzzzzzz_1, g_y_0_yyzzzzzzz_0, g_y_0_yyzzzzzzz_1, g_y_0_yzzzzzzzz_0, g_y_0_yzzzzzzzz_1, g_y_0_zzzzzzzzz_0, g_y_0_zzzzzzzzz_1, g_yy_0_xxxxxxxx_1, g_yy_0_xxxxxxxxx_1, g_yy_0_xxxxxxxxy_1, g_yy_0_xxxxxxxxz_1, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxxyy_1, g_yy_0_xxxxxxxyz_1, g_yy_0_xxxxxxxz_1, g_yy_0_xxxxxxxzz_1, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyyy_1, g_yy_0_xxxxxxyyz_1, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxxyzz_1, g_yy_0_xxxxxxzz_1, g_yy_0_xxxxxxzzz_1, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyyy_1, g_yy_0_xxxxxyyyz_1, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyyzz_1, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxxyzzz_1, g_yy_0_xxxxxzzz_1, g_yy_0_xxxxxzzzz_1, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyyy_1, g_yy_0_xxxxyyyyz_1, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyyzz_1, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyyzzz_1, g_yy_0_xxxxyzzz_1, g_yy_0_xxxxyzzzz_1, g_yy_0_xxxxzzzz_1, g_yy_0_xxxxzzzzz_1, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyyy_1, g_yy_0_xxxyyyyyz_1, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyyzz_1, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyyzzz_1, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyyzzzz_1, g_yy_0_xxxyzzzz_1, g_yy_0_xxxyzzzzz_1, g_yy_0_xxxzzzzz_1, g_yy_0_xxxzzzzzz_1, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyyy_1, g_yy_0_xxyyyyyyz_1, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyyzz_1, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyyzzz_1, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyyzzzz_1, g_yy_0_xxyyzzzz_1, g_yy_0_xxyyzzzzz_1, g_yy_0_xxyzzzzz_1, g_yy_0_xxyzzzzzz_1, g_yy_0_xxzzzzzz_1, g_yy_0_xxzzzzzzz_1, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyyy_1, g_yy_0_xyyyyyyyz_1, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyyzz_1, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyyzzz_1, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyyzzzz_1, g_yy_0_xyyyzzzz_1, g_yy_0_xyyyzzzzz_1, g_yy_0_xyyzzzzz_1, g_yy_0_xyyzzzzzz_1, g_yy_0_xyzzzzzz_1, g_yy_0_xyzzzzzzz_1, g_yy_0_xzzzzzzz_1, g_yy_0_xzzzzzzzz_1, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyyy_1, g_yy_0_yyyyyyyyz_1, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyyzz_1, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyyzzz_1, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyyzzzz_1, g_yy_0_yyyyzzzz_1, g_yy_0_yyyyzzzzz_1, g_yy_0_yyyzzzzz_1, g_yy_0_yyyzzzzzz_1, g_yy_0_yyzzzzzz_1, g_yy_0_yyzzzzzzz_1, g_yy_0_yzzzzzzz_1, g_yy_0_yzzzzzzzz_1, g_yy_0_zzzzzzzz_1, g_yy_0_zzzzzzzzz_1, g_yyy_0_xxxxxxxxx_0, g_yyy_0_xxxxxxxxy_0, g_yyy_0_xxxxxxxxz_0, g_yyy_0_xxxxxxxyy_0, g_yyy_0_xxxxxxxyz_0, g_yyy_0_xxxxxxxzz_0, g_yyy_0_xxxxxxyyy_0, g_yyy_0_xxxxxxyyz_0, g_yyy_0_xxxxxxyzz_0, g_yyy_0_xxxxxxzzz_0, g_yyy_0_xxxxxyyyy_0, g_yyy_0_xxxxxyyyz_0, g_yyy_0_xxxxxyyzz_0, g_yyy_0_xxxxxyzzz_0, g_yyy_0_xxxxxzzzz_0, g_yyy_0_xxxxyyyyy_0, g_yyy_0_xxxxyyyyz_0, g_yyy_0_xxxxyyyzz_0, g_yyy_0_xxxxyyzzz_0, g_yyy_0_xxxxyzzzz_0, g_yyy_0_xxxxzzzzz_0, g_yyy_0_xxxyyyyyy_0, g_yyy_0_xxxyyyyyz_0, g_yyy_0_xxxyyyyzz_0, g_yyy_0_xxxyyyzzz_0, g_yyy_0_xxxyyzzzz_0, g_yyy_0_xxxyzzzzz_0, g_yyy_0_xxxzzzzzz_0, g_yyy_0_xxyyyyyyy_0, g_yyy_0_xxyyyyyyz_0, g_yyy_0_xxyyyyyzz_0, g_yyy_0_xxyyyyzzz_0, g_yyy_0_xxyyyzzzz_0, g_yyy_0_xxyyzzzzz_0, g_yyy_0_xxyzzzzzz_0, g_yyy_0_xxzzzzzzz_0, g_yyy_0_xyyyyyyyy_0, g_yyy_0_xyyyyyyyz_0, g_yyy_0_xyyyyyyzz_0, g_yyy_0_xyyyyyzzz_0, g_yyy_0_xyyyyzzzz_0, g_yyy_0_xyyyzzzzz_0, g_yyy_0_xyyzzzzzz_0, g_yyy_0_xyzzzzzzz_0, g_yyy_0_xzzzzzzzz_0, g_yyy_0_yyyyyyyyy_0, g_yyy_0_yyyyyyyyz_0, g_yyy_0_yyyyyyyzz_0, g_yyy_0_yyyyyyzzz_0, g_yyy_0_yyyyyzzzz_0, g_yyy_0_yyyyzzzzz_0, g_yyy_0_yyyzzzzzz_0, g_yyy_0_yyzzzzzzz_0, g_yyy_0_yzzzzzzzz_0, g_yyy_0_zzzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xxxxxxxxx_0[i] = 2.0 * g_y_0_xxxxxxxxx_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxxx_1[i] * fz_be_0 + g_yy_0_xxxxxxxxx_1[i] * wa_y[i];

        g_yyy_0_xxxxxxxxy_0[i] = 2.0 * g_y_0_xxxxxxxxy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxxy_1[i] * fz_be_0 + g_yy_0_xxxxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxxxy_1[i] * wa_y[i];

        g_yyy_0_xxxxxxxxz_0[i] = 2.0 * g_y_0_xxxxxxxxz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxxz_1[i] * fz_be_0 + g_yy_0_xxxxxxxxz_1[i] * wa_y[i];

        g_yyy_0_xxxxxxxyy_0[i] = 2.0 * g_y_0_xxxxxxxyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxyy_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxxxyy_1[i] * wa_y[i];

        g_yyy_0_xxxxxxxyz_0[i] = 2.0 * g_y_0_xxxxxxxyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxyz_1[i] * fz_be_0 + g_yy_0_xxxxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxxxyz_1[i] * wa_y[i];

        g_yyy_0_xxxxxxxzz_0[i] = 2.0 * g_y_0_xxxxxxxzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxxzz_1[i] * fz_be_0 + g_yy_0_xxxxxxxzz_1[i] * wa_y[i];

        g_yyy_0_xxxxxxyyy_0[i] = 2.0 * g_y_0_xxxxxxyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxyyy_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxxxyyy_1[i] * wa_y[i];

        g_yyy_0_xxxxxxyyz_0[i] = 2.0 * g_y_0_xxxxxxyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxxxyyz_1[i] * wa_y[i];

        g_yyy_0_xxxxxxyzz_0[i] = 2.0 * g_y_0_xxxxxxyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxyzz_1[i] * fz_be_0 + g_yy_0_xxxxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxxxyzz_1[i] * wa_y[i];

        g_yyy_0_xxxxxxzzz_0[i] = 2.0 * g_y_0_xxxxxxzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxxzzz_1[i] * fz_be_0 + g_yy_0_xxxxxxzzz_1[i] * wa_y[i];

        g_yyy_0_xxxxxyyyy_0[i] = 2.0 * g_y_0_xxxxxyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyyyy_1[i] * fz_be_0 + 4.0 * g_yy_0_xxxxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxxxyyyy_1[i] * wa_y[i];

        g_yyy_0_xxxxxyyyz_0[i] = 2.0 * g_y_0_xxxxxyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxxxyyyz_1[i] * wa_y[i];

        g_yyy_0_xxxxxyyzz_0[i] = 2.0 * g_y_0_xxxxxyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxxxyyzz_1[i] * wa_y[i];

        g_yyy_0_xxxxxyzzz_0[i] = 2.0 * g_y_0_xxxxxyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxyzzz_1[i] * fz_be_0 + g_yy_0_xxxxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxxxyzzz_1[i] * wa_y[i];

        g_yyy_0_xxxxxzzzz_0[i] = 2.0 * g_y_0_xxxxxzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxxzzzz_1[i] * fz_be_0 + g_yy_0_xxxxxzzzz_1[i] * wa_y[i];

        g_yyy_0_xxxxyyyyy_0[i] = 2.0 * g_y_0_xxxxyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyyyy_1[i] * fz_be_0 + 5.0 * g_yy_0_xxxxyyyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyyyy_1[i] * wa_y[i];

        g_yyy_0_xxxxyyyyz_0[i] = 2.0 * g_y_0_xxxxyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yy_0_xxxxyyyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyyyz_1[i] * wa_y[i];

        g_yyy_0_xxxxyyyzz_0[i] = 2.0 * g_y_0_xxxxyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxxyyzz_1[i] * fi_acd_0 + g_yy_0_xxxxyyyzz_1[i] * wa_y[i];

        g_yyy_0_xxxxyyzzz_0[i] = 2.0 * g_y_0_xxxxyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxxyzzz_1[i] * fi_acd_0 + g_yy_0_xxxxyyzzz_1[i] * wa_y[i];

        g_yyy_0_xxxxyzzzz_0[i] = 2.0 * g_y_0_xxxxyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxyzzzz_1[i] * fz_be_0 + g_yy_0_xxxxzzzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzzzz_1[i] * wa_y[i];

        g_yyy_0_xxxxzzzzz_0[i] = 2.0 * g_y_0_xxxxzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxzzzzz_1[i] * fz_be_0 + g_yy_0_xxxxzzzzz_1[i] * wa_y[i];

        g_yyy_0_xxxyyyyyy_0[i] = 2.0 * g_y_0_xxxyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyyyy_1[i] * fz_be_0 + 6.0 * g_yy_0_xxxyyyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyyyy_1[i] * wa_y[i];

        g_yyy_0_xxxyyyyyz_0[i] = 2.0 * g_y_0_xxxyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yy_0_xxxyyyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyyyz_1[i] * wa_y[i];

        g_yyy_0_xxxyyyyzz_0[i] = 2.0 * g_y_0_xxxyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yy_0_xxxyyyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyyyzz_1[i] * wa_y[i];

        g_yyy_0_xxxyyyzzz_0[i] = 2.0 * g_y_0_xxxyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxxyyzzz_1[i] * fi_acd_0 + g_yy_0_xxxyyyzzz_1[i] * wa_y[i];

        g_yyy_0_xxxyyzzzz_0[i] = 2.0 * g_y_0_xxxyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxyzzzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzzzz_1[i] * wa_y[i];

        g_yyy_0_xxxyzzzzz_0[i] = 2.0 * g_y_0_xxxyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyzzzzz_1[i] * fz_be_0 + g_yy_0_xxxzzzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzzzz_1[i] * wa_y[i];

        g_yyy_0_xxxzzzzzz_0[i] = 2.0 * g_y_0_xxxzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxzzzzzz_1[i] * fz_be_0 + g_yy_0_xxxzzzzzz_1[i] * wa_y[i];

        g_yyy_0_xxyyyyyyy_0[i] = 2.0 * g_y_0_xxyyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyyyy_1[i] * fz_be_0 + 7.0 * g_yy_0_xxyyyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyyyy_1[i] * wa_y[i];

        g_yyy_0_xxyyyyyyz_0[i] = 2.0 * g_y_0_xxyyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yy_0_xxyyyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyyyz_1[i] * wa_y[i];

        g_yyy_0_xxyyyyyzz_0[i] = 2.0 * g_y_0_xxyyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yy_0_xxyyyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyyyzz_1[i] * wa_y[i];

        g_yyy_0_xxyyyyzzz_0[i] = 2.0 * g_y_0_xxyyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yy_0_xxyyyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyyyzzz_1[i] * wa_y[i];

        g_yyy_0_xxyyyzzzz_0[i] = 2.0 * g_y_0_xxyyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xxyyzzzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzzzz_1[i] * wa_y[i];

        g_yyy_0_xxyyzzzzz_0[i] = 2.0 * g_y_0_xxyyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxyzzzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzzzz_1[i] * wa_y[i];

        g_yyy_0_xxyzzzzzz_0[i] = 2.0 * g_y_0_xxyzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyzzzzzz_1[i] * fz_be_0 + g_yy_0_xxzzzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzzzz_1[i] * wa_y[i];

        g_yyy_0_xxzzzzzzz_0[i] = 2.0 * g_y_0_xxzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxzzzzzzz_1[i] * fz_be_0 + g_yy_0_xxzzzzzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyyyyyy_0[i] = 2.0 * g_y_0_xyyyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyyyy_1[i] * fz_be_0 + 8.0 * g_yy_0_xyyyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyyyy_1[i] * wa_y[i];

        g_yyy_0_xyyyyyyyz_0[i] = 2.0 * g_y_0_xyyyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yy_0_xyyyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyyyz_1[i] * wa_y[i];

        g_yyy_0_xyyyyyyzz_0[i] = 2.0 * g_y_0_xyyyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yy_0_xyyyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyyyzz_1[i] * wa_y[i];

        g_yyy_0_xyyyyyzzz_0[i] = 2.0 * g_y_0_xyyyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yy_0_xyyyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyyyzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyyzzzz_0[i] = 2.0 * g_y_0_xyyyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yy_0_xyyyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyzzzzz_0[i] = 2.0 * g_y_0_xyyyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_xyyzzzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzzzz_1[i] * wa_y[i];

        g_yyy_0_xyyzzzzzz_0[i] = 2.0 * g_y_0_xyyzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xyzzzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzzzz_1[i] * wa_y[i];

        g_yyy_0_xyzzzzzzz_0[i] = 2.0 * g_y_0_xyzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyzzzzzzz_1[i] * fz_be_0 + g_yy_0_xzzzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzzzz_1[i] * wa_y[i];

        g_yyy_0_xzzzzzzzz_0[i] = 2.0 * g_y_0_xzzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xzzzzzzzz_1[i] * fz_be_0 + g_yy_0_xzzzzzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyyyyyy_0[i] = 2.0 * g_y_0_yyyyyyyyy_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyyyy_1[i] * fz_be_0 + 9.0 * g_yy_0_yyyyyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyyyyy_1[i] * wa_y[i];

        g_yyy_0_yyyyyyyyz_0[i] = 2.0 * g_y_0_yyyyyyyyz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyyyz_1[i] * fz_be_0 + 8.0 * g_yy_0_yyyyyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyyyyyz_1[i] * wa_y[i];

        g_yyy_0_yyyyyyyzz_0[i] = 2.0 * g_y_0_yyyyyyyzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyyzz_1[i] * fz_be_0 + 7.0 * g_yy_0_yyyyyyzz_1[i] * fi_acd_0 + g_yy_0_yyyyyyyzz_1[i] * wa_y[i];

        g_yyy_0_yyyyyyzzz_0[i] = 2.0 * g_y_0_yyyyyyzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyyzzz_1[i] * fz_be_0 + 6.0 * g_yy_0_yyyyyzzz_1[i] * fi_acd_0 + g_yy_0_yyyyyyzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyyzzzz_0[i] = 2.0 * g_y_0_yyyyyzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyyzzzz_1[i] * fz_be_0 + 5.0 * g_yy_0_yyyyzzzz_1[i] * fi_acd_0 + g_yy_0_yyyyyzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyzzzzz_0[i] = 2.0 * g_y_0_yyyyzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyzzzzz_1[i] * fz_be_0 + 4.0 * g_yy_0_yyyzzzzz_1[i] * fi_acd_0 + g_yy_0_yyyyzzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyzzzzzz_0[i] = 2.0 * g_y_0_yyyzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyzzzzzz_1[i] * fz_be_0 + 3.0 * g_yy_0_yyzzzzzz_1[i] * fi_acd_0 + g_yy_0_yyyzzzzzz_1[i] * wa_y[i];

        g_yyy_0_yyzzzzzzz_0[i] = 2.0 * g_y_0_yyzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyzzzzzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_yzzzzzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzzzzz_1[i] * wa_y[i];

        g_yyy_0_yzzzzzzzz_0[i] = 2.0 * g_y_0_yzzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yzzzzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzzzzz_1[i] * wa_y[i];

        g_yyy_0_zzzzzzzzz_0[i] = 2.0 * g_y_0_zzzzzzzzz_0[i] * fbe_0 - 2.0 * g_y_0_zzzzzzzzz_1[i] * fz_be_0 + g_yy_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 385-440 components of targeted buffer : FSM

    auto g_yyz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 385);

    auto g_yyz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 386);

    auto g_yyz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 387);

    auto g_yyz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 388);

    auto g_yyz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 389);

    auto g_yyz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 390);

    auto g_yyz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 391);

    auto g_yyz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 392);

    auto g_yyz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 393);

    auto g_yyz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 394);

    auto g_yyz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 395);

    auto g_yyz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 396);

    auto g_yyz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 397);

    auto g_yyz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 398);

    auto g_yyz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 399);

    auto g_yyz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 400);

    auto g_yyz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 401);

    auto g_yyz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 402);

    auto g_yyz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 403);

    auto g_yyz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 404);

    auto g_yyz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 405);

    auto g_yyz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 406);

    auto g_yyz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 407);

    auto g_yyz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 408);

    auto g_yyz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 409);

    auto g_yyz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 410);

    auto g_yyz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 411);

    auto g_yyz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 412);

    auto g_yyz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 413);

    auto g_yyz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 414);

    auto g_yyz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 415);

    auto g_yyz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 416);

    auto g_yyz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 417);

    auto g_yyz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 418);

    auto g_yyz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 419);

    auto g_yyz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 420);

    auto g_yyz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 421);

    auto g_yyz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 422);

    auto g_yyz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 423);

    auto g_yyz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 424);

    auto g_yyz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 425);

    auto g_yyz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 426);

    auto g_yyz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 427);

    auto g_yyz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 428);

    auto g_yyz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 429);

    auto g_yyz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 430);

    auto g_yyz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 431);

    auto g_yyz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 432);

    auto g_yyz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 433);

    auto g_yyz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 434);

    auto g_yyz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 435);

    auto g_yyz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 436);

    auto g_yyz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 437);

    auto g_yyz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 438);

    auto g_yyz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 439);

    #pragma omp simd aligned(g_yy_0_xxxxxxxx_1, g_yy_0_xxxxxxxxx_1, g_yy_0_xxxxxxxxy_1, g_yy_0_xxxxxxxxz_1, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxxyy_1, g_yy_0_xxxxxxxyz_1, g_yy_0_xxxxxxxz_1, g_yy_0_xxxxxxxzz_1, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyyy_1, g_yy_0_xxxxxxyyz_1, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxxyzz_1, g_yy_0_xxxxxxzz_1, g_yy_0_xxxxxxzzz_1, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyyy_1, g_yy_0_xxxxxyyyz_1, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyyzz_1, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxxyzzz_1, g_yy_0_xxxxxzzz_1, g_yy_0_xxxxxzzzz_1, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyyy_1, g_yy_0_xxxxyyyyz_1, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyyzz_1, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyyzzz_1, g_yy_0_xxxxyzzz_1, g_yy_0_xxxxyzzzz_1, g_yy_0_xxxxzzzz_1, g_yy_0_xxxxzzzzz_1, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyyy_1, g_yy_0_xxxyyyyyz_1, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyyzz_1, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyyzzz_1, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyyzzzz_1, g_yy_0_xxxyzzzz_1, g_yy_0_xxxyzzzzz_1, g_yy_0_xxxzzzzz_1, g_yy_0_xxxzzzzzz_1, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyyy_1, g_yy_0_xxyyyyyyz_1, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyyzz_1, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyyzzz_1, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyyzzzz_1, g_yy_0_xxyyzzzz_1, g_yy_0_xxyyzzzzz_1, g_yy_0_xxyzzzzz_1, g_yy_0_xxyzzzzzz_1, g_yy_0_xxzzzzzz_1, g_yy_0_xxzzzzzzz_1, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyyy_1, g_yy_0_xyyyyyyyz_1, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyyzz_1, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyyzzz_1, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyyzzzz_1, g_yy_0_xyyyzzzz_1, g_yy_0_xyyyzzzzz_1, g_yy_0_xyyzzzzz_1, g_yy_0_xyyzzzzzz_1, g_yy_0_xyzzzzzz_1, g_yy_0_xyzzzzzzz_1, g_yy_0_xzzzzzzz_1, g_yy_0_xzzzzzzzz_1, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyyy_1, g_yy_0_yyyyyyyyz_1, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyyzz_1, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyyzzz_1, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyyzzzz_1, g_yy_0_yyyyzzzz_1, g_yy_0_yyyyzzzzz_1, g_yy_0_yyyzzzzz_1, g_yy_0_yyyzzzzzz_1, g_yy_0_yyzzzzzz_1, g_yy_0_yyzzzzzzz_1, g_yy_0_yzzzzzzz_1, g_yy_0_yzzzzzzzz_1, g_yy_0_zzzzzzzz_1, g_yy_0_zzzzzzzzz_1, g_yyz_0_xxxxxxxxx_0, g_yyz_0_xxxxxxxxy_0, g_yyz_0_xxxxxxxxz_0, g_yyz_0_xxxxxxxyy_0, g_yyz_0_xxxxxxxyz_0, g_yyz_0_xxxxxxxzz_0, g_yyz_0_xxxxxxyyy_0, g_yyz_0_xxxxxxyyz_0, g_yyz_0_xxxxxxyzz_0, g_yyz_0_xxxxxxzzz_0, g_yyz_0_xxxxxyyyy_0, g_yyz_0_xxxxxyyyz_0, g_yyz_0_xxxxxyyzz_0, g_yyz_0_xxxxxyzzz_0, g_yyz_0_xxxxxzzzz_0, g_yyz_0_xxxxyyyyy_0, g_yyz_0_xxxxyyyyz_0, g_yyz_0_xxxxyyyzz_0, g_yyz_0_xxxxyyzzz_0, g_yyz_0_xxxxyzzzz_0, g_yyz_0_xxxxzzzzz_0, g_yyz_0_xxxyyyyyy_0, g_yyz_0_xxxyyyyyz_0, g_yyz_0_xxxyyyyzz_0, g_yyz_0_xxxyyyzzz_0, g_yyz_0_xxxyyzzzz_0, g_yyz_0_xxxyzzzzz_0, g_yyz_0_xxxzzzzzz_0, g_yyz_0_xxyyyyyyy_0, g_yyz_0_xxyyyyyyz_0, g_yyz_0_xxyyyyyzz_0, g_yyz_0_xxyyyyzzz_0, g_yyz_0_xxyyyzzzz_0, g_yyz_0_xxyyzzzzz_0, g_yyz_0_xxyzzzzzz_0, g_yyz_0_xxzzzzzzz_0, g_yyz_0_xyyyyyyyy_0, g_yyz_0_xyyyyyyyz_0, g_yyz_0_xyyyyyyzz_0, g_yyz_0_xyyyyyzzz_0, g_yyz_0_xyyyyzzzz_0, g_yyz_0_xyyyzzzzz_0, g_yyz_0_xyyzzzzzz_0, g_yyz_0_xyzzzzzzz_0, g_yyz_0_xzzzzzzzz_0, g_yyz_0_yyyyyyyyy_0, g_yyz_0_yyyyyyyyz_0, g_yyz_0_yyyyyyyzz_0, g_yyz_0_yyyyyyzzz_0, g_yyz_0_yyyyyzzzz_0, g_yyz_0_yyyyzzzzz_0, g_yyz_0_yyyzzzzzz_0, g_yyz_0_yyzzzzzzz_0, g_yyz_0_yzzzzzzzz_0, g_yyz_0_zzzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xxxxxxxxx_0[i] = g_yy_0_xxxxxxxxx_1[i] * wa_z[i];

        g_yyz_0_xxxxxxxxy_0[i] = g_yy_0_xxxxxxxxy_1[i] * wa_z[i];

        g_yyz_0_xxxxxxxxz_0[i] = g_yy_0_xxxxxxxx_1[i] * fi_acd_0 + g_yy_0_xxxxxxxxz_1[i] * wa_z[i];

        g_yyz_0_xxxxxxxyy_0[i] = g_yy_0_xxxxxxxyy_1[i] * wa_z[i];

        g_yyz_0_xxxxxxxyz_0[i] = g_yy_0_xxxxxxxy_1[i] * fi_acd_0 + g_yy_0_xxxxxxxyz_1[i] * wa_z[i];

        g_yyz_0_xxxxxxxzz_0[i] = 2.0 * g_yy_0_xxxxxxxz_1[i] * fi_acd_0 + g_yy_0_xxxxxxxzz_1[i] * wa_z[i];

        g_yyz_0_xxxxxxyyy_0[i] = g_yy_0_xxxxxxyyy_1[i] * wa_z[i];

        g_yyz_0_xxxxxxyyz_0[i] = g_yy_0_xxxxxxyy_1[i] * fi_acd_0 + g_yy_0_xxxxxxyyz_1[i] * wa_z[i];

        g_yyz_0_xxxxxxyzz_0[i] = 2.0 * g_yy_0_xxxxxxyz_1[i] * fi_acd_0 + g_yy_0_xxxxxxyzz_1[i] * wa_z[i];

        g_yyz_0_xxxxxxzzz_0[i] = 3.0 * g_yy_0_xxxxxxzz_1[i] * fi_acd_0 + g_yy_0_xxxxxxzzz_1[i] * wa_z[i];

        g_yyz_0_xxxxxyyyy_0[i] = g_yy_0_xxxxxyyyy_1[i] * wa_z[i];

        g_yyz_0_xxxxxyyyz_0[i] = g_yy_0_xxxxxyyy_1[i] * fi_acd_0 + g_yy_0_xxxxxyyyz_1[i] * wa_z[i];

        g_yyz_0_xxxxxyyzz_0[i] = 2.0 * g_yy_0_xxxxxyyz_1[i] * fi_acd_0 + g_yy_0_xxxxxyyzz_1[i] * wa_z[i];

        g_yyz_0_xxxxxyzzz_0[i] = 3.0 * g_yy_0_xxxxxyzz_1[i] * fi_acd_0 + g_yy_0_xxxxxyzzz_1[i] * wa_z[i];

        g_yyz_0_xxxxxzzzz_0[i] = 4.0 * g_yy_0_xxxxxzzz_1[i] * fi_acd_0 + g_yy_0_xxxxxzzzz_1[i] * wa_z[i];

        g_yyz_0_xxxxyyyyy_0[i] = g_yy_0_xxxxyyyyy_1[i] * wa_z[i];

        g_yyz_0_xxxxyyyyz_0[i] = g_yy_0_xxxxyyyy_1[i] * fi_acd_0 + g_yy_0_xxxxyyyyz_1[i] * wa_z[i];

        g_yyz_0_xxxxyyyzz_0[i] = 2.0 * g_yy_0_xxxxyyyz_1[i] * fi_acd_0 + g_yy_0_xxxxyyyzz_1[i] * wa_z[i];

        g_yyz_0_xxxxyyzzz_0[i] = 3.0 * g_yy_0_xxxxyyzz_1[i] * fi_acd_0 + g_yy_0_xxxxyyzzz_1[i] * wa_z[i];

        g_yyz_0_xxxxyzzzz_0[i] = 4.0 * g_yy_0_xxxxyzzz_1[i] * fi_acd_0 + g_yy_0_xxxxyzzzz_1[i] * wa_z[i];

        g_yyz_0_xxxxzzzzz_0[i] = 5.0 * g_yy_0_xxxxzzzz_1[i] * fi_acd_0 + g_yy_0_xxxxzzzzz_1[i] * wa_z[i];

        g_yyz_0_xxxyyyyyy_0[i] = g_yy_0_xxxyyyyyy_1[i] * wa_z[i];

        g_yyz_0_xxxyyyyyz_0[i] = g_yy_0_xxxyyyyy_1[i] * fi_acd_0 + g_yy_0_xxxyyyyyz_1[i] * wa_z[i];

        g_yyz_0_xxxyyyyzz_0[i] = 2.0 * g_yy_0_xxxyyyyz_1[i] * fi_acd_0 + g_yy_0_xxxyyyyzz_1[i] * wa_z[i];

        g_yyz_0_xxxyyyzzz_0[i] = 3.0 * g_yy_0_xxxyyyzz_1[i] * fi_acd_0 + g_yy_0_xxxyyyzzz_1[i] * wa_z[i];

        g_yyz_0_xxxyyzzzz_0[i] = 4.0 * g_yy_0_xxxyyzzz_1[i] * fi_acd_0 + g_yy_0_xxxyyzzzz_1[i] * wa_z[i];

        g_yyz_0_xxxyzzzzz_0[i] = 5.0 * g_yy_0_xxxyzzzz_1[i] * fi_acd_0 + g_yy_0_xxxyzzzzz_1[i] * wa_z[i];

        g_yyz_0_xxxzzzzzz_0[i] = 6.0 * g_yy_0_xxxzzzzz_1[i] * fi_acd_0 + g_yy_0_xxxzzzzzz_1[i] * wa_z[i];

        g_yyz_0_xxyyyyyyy_0[i] = g_yy_0_xxyyyyyyy_1[i] * wa_z[i];

        g_yyz_0_xxyyyyyyz_0[i] = g_yy_0_xxyyyyyy_1[i] * fi_acd_0 + g_yy_0_xxyyyyyyz_1[i] * wa_z[i];

        g_yyz_0_xxyyyyyzz_0[i] = 2.0 * g_yy_0_xxyyyyyz_1[i] * fi_acd_0 + g_yy_0_xxyyyyyzz_1[i] * wa_z[i];

        g_yyz_0_xxyyyyzzz_0[i] = 3.0 * g_yy_0_xxyyyyzz_1[i] * fi_acd_0 + g_yy_0_xxyyyyzzz_1[i] * wa_z[i];

        g_yyz_0_xxyyyzzzz_0[i] = 4.0 * g_yy_0_xxyyyzzz_1[i] * fi_acd_0 + g_yy_0_xxyyyzzzz_1[i] * wa_z[i];

        g_yyz_0_xxyyzzzzz_0[i] = 5.0 * g_yy_0_xxyyzzzz_1[i] * fi_acd_0 + g_yy_0_xxyyzzzzz_1[i] * wa_z[i];

        g_yyz_0_xxyzzzzzz_0[i] = 6.0 * g_yy_0_xxyzzzzz_1[i] * fi_acd_0 + g_yy_0_xxyzzzzzz_1[i] * wa_z[i];

        g_yyz_0_xxzzzzzzz_0[i] = 7.0 * g_yy_0_xxzzzzzz_1[i] * fi_acd_0 + g_yy_0_xxzzzzzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyyyyyy_0[i] = g_yy_0_xyyyyyyyy_1[i] * wa_z[i];

        g_yyz_0_xyyyyyyyz_0[i] = g_yy_0_xyyyyyyy_1[i] * fi_acd_0 + g_yy_0_xyyyyyyyz_1[i] * wa_z[i];

        g_yyz_0_xyyyyyyzz_0[i] = 2.0 * g_yy_0_xyyyyyyz_1[i] * fi_acd_0 + g_yy_0_xyyyyyyzz_1[i] * wa_z[i];

        g_yyz_0_xyyyyyzzz_0[i] = 3.0 * g_yy_0_xyyyyyzz_1[i] * fi_acd_0 + g_yy_0_xyyyyyzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyyzzzz_0[i] = 4.0 * g_yy_0_xyyyyzzz_1[i] * fi_acd_0 + g_yy_0_xyyyyzzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyzzzzz_0[i] = 5.0 * g_yy_0_xyyyzzzz_1[i] * fi_acd_0 + g_yy_0_xyyyzzzzz_1[i] * wa_z[i];

        g_yyz_0_xyyzzzzzz_0[i] = 6.0 * g_yy_0_xyyzzzzz_1[i] * fi_acd_0 + g_yy_0_xyyzzzzzz_1[i] * wa_z[i];

        g_yyz_0_xyzzzzzzz_0[i] = 7.0 * g_yy_0_xyzzzzzz_1[i] * fi_acd_0 + g_yy_0_xyzzzzzzz_1[i] * wa_z[i];

        g_yyz_0_xzzzzzzzz_0[i] = 8.0 * g_yy_0_xzzzzzzz_1[i] * fi_acd_0 + g_yy_0_xzzzzzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyyyyyy_0[i] = g_yy_0_yyyyyyyyy_1[i] * wa_z[i];

        g_yyz_0_yyyyyyyyz_0[i] = g_yy_0_yyyyyyyy_1[i] * fi_acd_0 + g_yy_0_yyyyyyyyz_1[i] * wa_z[i];

        g_yyz_0_yyyyyyyzz_0[i] = 2.0 * g_yy_0_yyyyyyyz_1[i] * fi_acd_0 + g_yy_0_yyyyyyyzz_1[i] * wa_z[i];

        g_yyz_0_yyyyyyzzz_0[i] = 3.0 * g_yy_0_yyyyyyzz_1[i] * fi_acd_0 + g_yy_0_yyyyyyzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyyzzzz_0[i] = 4.0 * g_yy_0_yyyyyzzz_1[i] * fi_acd_0 + g_yy_0_yyyyyzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyzzzzz_0[i] = 5.0 * g_yy_0_yyyyzzzz_1[i] * fi_acd_0 + g_yy_0_yyyyzzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyzzzzzz_0[i] = 6.0 * g_yy_0_yyyzzzzz_1[i] * fi_acd_0 + g_yy_0_yyyzzzzzz_1[i] * wa_z[i];

        g_yyz_0_yyzzzzzzz_0[i] = 7.0 * g_yy_0_yyzzzzzz_1[i] * fi_acd_0 + g_yy_0_yyzzzzzzz_1[i] * wa_z[i];

        g_yyz_0_yzzzzzzzz_0[i] = 8.0 * g_yy_0_yzzzzzzz_1[i] * fi_acd_0 + g_yy_0_yzzzzzzzz_1[i] * wa_z[i];

        g_yyz_0_zzzzzzzzz_0[i] = 9.0 * g_yy_0_zzzzzzzz_1[i] * fi_acd_0 + g_yy_0_zzzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 440-495 components of targeted buffer : FSM

    auto g_yzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 440);

    auto g_yzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 441);

    auto g_yzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 442);

    auto g_yzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 443);

    auto g_yzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 444);

    auto g_yzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 445);

    auto g_yzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 446);

    auto g_yzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 447);

    auto g_yzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 448);

    auto g_yzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 449);

    auto g_yzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 450);

    auto g_yzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 451);

    auto g_yzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 452);

    auto g_yzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 453);

    auto g_yzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 454);

    auto g_yzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 455);

    auto g_yzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 456);

    auto g_yzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 457);

    auto g_yzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 458);

    auto g_yzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 459);

    auto g_yzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 460);

    auto g_yzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 461);

    auto g_yzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 462);

    auto g_yzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 463);

    auto g_yzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 464);

    auto g_yzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 465);

    auto g_yzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 466);

    auto g_yzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 467);

    auto g_yzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 468);

    auto g_yzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 469);

    auto g_yzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 470);

    auto g_yzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 471);

    auto g_yzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 472);

    auto g_yzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 473);

    auto g_yzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 474);

    auto g_yzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 475);

    auto g_yzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 476);

    auto g_yzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 477);

    auto g_yzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 478);

    auto g_yzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 479);

    auto g_yzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 480);

    auto g_yzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 481);

    auto g_yzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 482);

    auto g_yzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 483);

    auto g_yzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 484);

    auto g_yzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 485);

    auto g_yzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 486);

    auto g_yzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 487);

    auto g_yzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 488);

    auto g_yzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 489);

    auto g_yzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 490);

    auto g_yzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 491);

    auto g_yzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 492);

    auto g_yzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 493);

    auto g_yzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 494);

    #pragma omp simd aligned(g_yzz_0_xxxxxxxxx_0, g_yzz_0_xxxxxxxxy_0, g_yzz_0_xxxxxxxxz_0, g_yzz_0_xxxxxxxyy_0, g_yzz_0_xxxxxxxyz_0, g_yzz_0_xxxxxxxzz_0, g_yzz_0_xxxxxxyyy_0, g_yzz_0_xxxxxxyyz_0, g_yzz_0_xxxxxxyzz_0, g_yzz_0_xxxxxxzzz_0, g_yzz_0_xxxxxyyyy_0, g_yzz_0_xxxxxyyyz_0, g_yzz_0_xxxxxyyzz_0, g_yzz_0_xxxxxyzzz_0, g_yzz_0_xxxxxzzzz_0, g_yzz_0_xxxxyyyyy_0, g_yzz_0_xxxxyyyyz_0, g_yzz_0_xxxxyyyzz_0, g_yzz_0_xxxxyyzzz_0, g_yzz_0_xxxxyzzzz_0, g_yzz_0_xxxxzzzzz_0, g_yzz_0_xxxyyyyyy_0, g_yzz_0_xxxyyyyyz_0, g_yzz_0_xxxyyyyzz_0, g_yzz_0_xxxyyyzzz_0, g_yzz_0_xxxyyzzzz_0, g_yzz_0_xxxyzzzzz_0, g_yzz_0_xxxzzzzzz_0, g_yzz_0_xxyyyyyyy_0, g_yzz_0_xxyyyyyyz_0, g_yzz_0_xxyyyyyzz_0, g_yzz_0_xxyyyyzzz_0, g_yzz_0_xxyyyzzzz_0, g_yzz_0_xxyyzzzzz_0, g_yzz_0_xxyzzzzzz_0, g_yzz_0_xxzzzzzzz_0, g_yzz_0_xyyyyyyyy_0, g_yzz_0_xyyyyyyyz_0, g_yzz_0_xyyyyyyzz_0, g_yzz_0_xyyyyyzzz_0, g_yzz_0_xyyyyzzzz_0, g_yzz_0_xyyyzzzzz_0, g_yzz_0_xyyzzzzzz_0, g_yzz_0_xyzzzzzzz_0, g_yzz_0_xzzzzzzzz_0, g_yzz_0_yyyyyyyyy_0, g_yzz_0_yyyyyyyyz_0, g_yzz_0_yyyyyyyzz_0, g_yzz_0_yyyyyyzzz_0, g_yzz_0_yyyyyzzzz_0, g_yzz_0_yyyyzzzzz_0, g_yzz_0_yyyzzzzzz_0, g_yzz_0_yyzzzzzzz_0, g_yzz_0_yzzzzzzzz_0, g_yzz_0_zzzzzzzzz_0, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxxx_1, g_zz_0_xxxxxxxxy_1, g_zz_0_xxxxxxxxz_1, g_zz_0_xxxxxxxy_1, g_zz_0_xxxxxxxyy_1, g_zz_0_xxxxxxxyz_1, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxxzz_1, g_zz_0_xxxxxxyy_1, g_zz_0_xxxxxxyyy_1, g_zz_0_xxxxxxyyz_1, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxyzz_1, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxxzzz_1, g_zz_0_xxxxxyyy_1, g_zz_0_xxxxxyyyy_1, g_zz_0_xxxxxyyyz_1, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyyzz_1, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxyzzz_1, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxxzzzz_1, g_zz_0_xxxxyyyy_1, g_zz_0_xxxxyyyyy_1, g_zz_0_xxxxyyyyz_1, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyyzz_1, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyyzzz_1, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxyzzzz_1, g_zz_0_xxxxzzzz_1, g_zz_0_xxxxzzzzz_1, g_zz_0_xxxyyyyy_1, g_zz_0_xxxyyyyyy_1, g_zz_0_xxxyyyyyz_1, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyyzz_1, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyyzzz_1, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyyzzzz_1, g_zz_0_xxxyzzzz_1, g_zz_0_xxxyzzzzz_1, g_zz_0_xxxzzzzz_1, g_zz_0_xxxzzzzzz_1, g_zz_0_xxyyyyyy_1, g_zz_0_xxyyyyyyy_1, g_zz_0_xxyyyyyyz_1, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyyzz_1, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyyzzz_1, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyyzzzz_1, g_zz_0_xxyyzzzz_1, g_zz_0_xxyyzzzzz_1, g_zz_0_xxyzzzzz_1, g_zz_0_xxyzzzzzz_1, g_zz_0_xxzzzzzz_1, g_zz_0_xxzzzzzzz_1, g_zz_0_xyyyyyyy_1, g_zz_0_xyyyyyyyy_1, g_zz_0_xyyyyyyyz_1, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyyzz_1, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyyzzz_1, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyyzzzz_1, g_zz_0_xyyyzzzz_1, g_zz_0_xyyyzzzzz_1, g_zz_0_xyyzzzzz_1, g_zz_0_xyyzzzzzz_1, g_zz_0_xyzzzzzz_1, g_zz_0_xyzzzzzzz_1, g_zz_0_xzzzzzzz_1, g_zz_0_xzzzzzzzz_1, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyyy_1, g_zz_0_yyyyyyyyz_1, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyyzz_1, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyyzzz_1, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyyzzzz_1, g_zz_0_yyyyzzzz_1, g_zz_0_yyyyzzzzz_1, g_zz_0_yyyzzzzz_1, g_zz_0_yyyzzzzzz_1, g_zz_0_yyzzzzzz_1, g_zz_0_yyzzzzzzz_1, g_zz_0_yzzzzzzz_1, g_zz_0_yzzzzzzzz_1, g_zz_0_zzzzzzzz_1, g_zz_0_zzzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xxxxxxxxx_0[i] = g_zz_0_xxxxxxxxx_1[i] * wa_y[i];

        g_yzz_0_xxxxxxxxy_0[i] = g_zz_0_xxxxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxxxy_1[i] * wa_y[i];

        g_yzz_0_xxxxxxxxz_0[i] = g_zz_0_xxxxxxxxz_1[i] * wa_y[i];

        g_yzz_0_xxxxxxxyy_0[i] = 2.0 * g_zz_0_xxxxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxxxyy_1[i] * wa_y[i];

        g_yzz_0_xxxxxxxyz_0[i] = g_zz_0_xxxxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxxxyz_1[i] * wa_y[i];

        g_yzz_0_xxxxxxxzz_0[i] = g_zz_0_xxxxxxxzz_1[i] * wa_y[i];

        g_yzz_0_xxxxxxyyy_0[i] = 3.0 * g_zz_0_xxxxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxxxyyy_1[i] * wa_y[i];

        g_yzz_0_xxxxxxyyz_0[i] = 2.0 * g_zz_0_xxxxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxxxyyz_1[i] * wa_y[i];

        g_yzz_0_xxxxxxyzz_0[i] = g_zz_0_xxxxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxxxyzz_1[i] * wa_y[i];

        g_yzz_0_xxxxxxzzz_0[i] = g_zz_0_xxxxxxzzz_1[i] * wa_y[i];

        g_yzz_0_xxxxxyyyy_0[i] = 4.0 * g_zz_0_xxxxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxxxyyyy_1[i] * wa_y[i];

        g_yzz_0_xxxxxyyyz_0[i] = 3.0 * g_zz_0_xxxxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxxxyyyz_1[i] * wa_y[i];

        g_yzz_0_xxxxxyyzz_0[i] = 2.0 * g_zz_0_xxxxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxxxyyzz_1[i] * wa_y[i];

        g_yzz_0_xxxxxyzzz_0[i] = g_zz_0_xxxxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxxxyzzz_1[i] * wa_y[i];

        g_yzz_0_xxxxxzzzz_0[i] = g_zz_0_xxxxxzzzz_1[i] * wa_y[i];

        g_yzz_0_xxxxyyyyy_0[i] = 5.0 * g_zz_0_xxxxyyyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyyyy_1[i] * wa_y[i];

        g_yzz_0_xxxxyyyyz_0[i] = 4.0 * g_zz_0_xxxxyyyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyyyz_1[i] * wa_y[i];

        g_yzz_0_xxxxyyyzz_0[i] = 3.0 * g_zz_0_xxxxyyzz_1[i] * fi_acd_0 + g_zz_0_xxxxyyyzz_1[i] * wa_y[i];

        g_yzz_0_xxxxyyzzz_0[i] = 2.0 * g_zz_0_xxxxyzzz_1[i] * fi_acd_0 + g_zz_0_xxxxyyzzz_1[i] * wa_y[i];

        g_yzz_0_xxxxyzzzz_0[i] = g_zz_0_xxxxzzzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzzzz_1[i] * wa_y[i];

        g_yzz_0_xxxxzzzzz_0[i] = g_zz_0_xxxxzzzzz_1[i] * wa_y[i];

        g_yzz_0_xxxyyyyyy_0[i] = 6.0 * g_zz_0_xxxyyyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyyyy_1[i] * wa_y[i];

        g_yzz_0_xxxyyyyyz_0[i] = 5.0 * g_zz_0_xxxyyyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyyyz_1[i] * wa_y[i];

        g_yzz_0_xxxyyyyzz_0[i] = 4.0 * g_zz_0_xxxyyyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyyyzz_1[i] * wa_y[i];

        g_yzz_0_xxxyyyzzz_0[i] = 3.0 * g_zz_0_xxxyyzzz_1[i] * fi_acd_0 + g_zz_0_xxxyyyzzz_1[i] * wa_y[i];

        g_yzz_0_xxxyyzzzz_0[i] = 2.0 * g_zz_0_xxxyzzzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzzzz_1[i] * wa_y[i];

        g_yzz_0_xxxyzzzzz_0[i] = g_zz_0_xxxzzzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzzzz_1[i] * wa_y[i];

        g_yzz_0_xxxzzzzzz_0[i] = g_zz_0_xxxzzzzzz_1[i] * wa_y[i];

        g_yzz_0_xxyyyyyyy_0[i] = 7.0 * g_zz_0_xxyyyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyyyy_1[i] * wa_y[i];

        g_yzz_0_xxyyyyyyz_0[i] = 6.0 * g_zz_0_xxyyyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyyyz_1[i] * wa_y[i];

        g_yzz_0_xxyyyyyzz_0[i] = 5.0 * g_zz_0_xxyyyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyyyzz_1[i] * wa_y[i];

        g_yzz_0_xxyyyyzzz_0[i] = 4.0 * g_zz_0_xxyyyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyyyzzz_1[i] * wa_y[i];

        g_yzz_0_xxyyyzzzz_0[i] = 3.0 * g_zz_0_xxyyzzzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzzzz_1[i] * wa_y[i];

        g_yzz_0_xxyyzzzzz_0[i] = 2.0 * g_zz_0_xxyzzzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzzzz_1[i] * wa_y[i];

        g_yzz_0_xxyzzzzzz_0[i] = g_zz_0_xxzzzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzzzz_1[i] * wa_y[i];

        g_yzz_0_xxzzzzzzz_0[i] = g_zz_0_xxzzzzzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyyyyyy_0[i] = 8.0 * g_zz_0_xyyyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyyyy_1[i] * wa_y[i];

        g_yzz_0_xyyyyyyyz_0[i] = 7.0 * g_zz_0_xyyyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyyyz_1[i] * wa_y[i];

        g_yzz_0_xyyyyyyzz_0[i] = 6.0 * g_zz_0_xyyyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyyyzz_1[i] * wa_y[i];

        g_yzz_0_xyyyyyzzz_0[i] = 5.0 * g_zz_0_xyyyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyyyzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyyzzzz_0[i] = 4.0 * g_zz_0_xyyyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyzzzzz_0[i] = 3.0 * g_zz_0_xyyzzzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzzzz_1[i] * wa_y[i];

        g_yzz_0_xyyzzzzzz_0[i] = 2.0 * g_zz_0_xyzzzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzzzz_1[i] * wa_y[i];

        g_yzz_0_xyzzzzzzz_0[i] = g_zz_0_xzzzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzzzz_1[i] * wa_y[i];

        g_yzz_0_xzzzzzzzz_0[i] = g_zz_0_xzzzzzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyyyyyy_0[i] = 9.0 * g_zz_0_yyyyyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyyyyy_1[i] * wa_y[i];

        g_yzz_0_yyyyyyyyz_0[i] = 8.0 * g_zz_0_yyyyyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyyyyyz_1[i] * wa_y[i];

        g_yzz_0_yyyyyyyzz_0[i] = 7.0 * g_zz_0_yyyyyyzz_1[i] * fi_acd_0 + g_zz_0_yyyyyyyzz_1[i] * wa_y[i];

        g_yzz_0_yyyyyyzzz_0[i] = 6.0 * g_zz_0_yyyyyzzz_1[i] * fi_acd_0 + g_zz_0_yyyyyyzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyyzzzz_0[i] = 5.0 * g_zz_0_yyyyzzzz_1[i] * fi_acd_0 + g_zz_0_yyyyyzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyzzzzz_0[i] = 4.0 * g_zz_0_yyyzzzzz_1[i] * fi_acd_0 + g_zz_0_yyyyzzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyzzzzzz_0[i] = 3.0 * g_zz_0_yyzzzzzz_1[i] * fi_acd_0 + g_zz_0_yyyzzzzzz_1[i] * wa_y[i];

        g_yzz_0_yyzzzzzzz_0[i] = 2.0 * g_zz_0_yzzzzzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzzzzz_1[i] * wa_y[i];

        g_yzz_0_yzzzzzzzz_0[i] = g_zz_0_zzzzzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzzzzz_1[i] * wa_y[i];

        g_yzz_0_zzzzzzzzz_0[i] = g_zz_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 495-550 components of targeted buffer : FSM

    auto g_zzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_fsm + 495);

    auto g_zzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_fsm + 496);

    auto g_zzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_fsm + 497);

    auto g_zzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_fsm + 498);

    auto g_zzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_fsm + 499);

    auto g_zzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_fsm + 500);

    auto g_zzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_fsm + 501);

    auto g_zzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_fsm + 502);

    auto g_zzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_fsm + 503);

    auto g_zzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_fsm + 504);

    auto g_zzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_fsm + 505);

    auto g_zzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_fsm + 506);

    auto g_zzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_fsm + 507);

    auto g_zzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_fsm + 508);

    auto g_zzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_fsm + 509);

    auto g_zzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 510);

    auto g_zzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 511);

    auto g_zzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 512);

    auto g_zzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 513);

    auto g_zzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 514);

    auto g_zzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 515);

    auto g_zzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 516);

    auto g_zzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 517);

    auto g_zzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 518);

    auto g_zzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 519);

    auto g_zzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 520);

    auto g_zzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 521);

    auto g_zzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 522);

    auto g_zzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 523);

    auto g_zzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 524);

    auto g_zzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 525);

    auto g_zzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 526);

    auto g_zzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 527);

    auto g_zzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 528);

    auto g_zzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 529);

    auto g_zzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 530);

    auto g_zzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 531);

    auto g_zzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 532);

    auto g_zzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 533);

    auto g_zzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 534);

    auto g_zzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 535);

    auto g_zzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 536);

    auto g_zzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 537);

    auto g_zzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 538);

    auto g_zzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 539);

    auto g_zzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_fsm + 540);

    auto g_zzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_fsm + 541);

    auto g_zzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_fsm + 542);

    auto g_zzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_fsm + 543);

    auto g_zzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_fsm + 544);

    auto g_zzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 545);

    auto g_zzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 546);

    auto g_zzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 547);

    auto g_zzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 548);

    auto g_zzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_fsm + 549);

    #pragma omp simd aligned(g_z_0_xxxxxxxxx_0, g_z_0_xxxxxxxxx_1, g_z_0_xxxxxxxxy_0, g_z_0_xxxxxxxxy_1, g_z_0_xxxxxxxxz_0, g_z_0_xxxxxxxxz_1, g_z_0_xxxxxxxyy_0, g_z_0_xxxxxxxyy_1, g_z_0_xxxxxxxyz_0, g_z_0_xxxxxxxyz_1, g_z_0_xxxxxxxzz_0, g_z_0_xxxxxxxzz_1, g_z_0_xxxxxxyyy_0, g_z_0_xxxxxxyyy_1, g_z_0_xxxxxxyyz_0, g_z_0_xxxxxxyyz_1, g_z_0_xxxxxxyzz_0, g_z_0_xxxxxxyzz_1, g_z_0_xxxxxxzzz_0, g_z_0_xxxxxxzzz_1, g_z_0_xxxxxyyyy_0, g_z_0_xxxxxyyyy_1, g_z_0_xxxxxyyyz_0, g_z_0_xxxxxyyyz_1, g_z_0_xxxxxyyzz_0, g_z_0_xxxxxyyzz_1, g_z_0_xxxxxyzzz_0, g_z_0_xxxxxyzzz_1, g_z_0_xxxxxzzzz_0, g_z_0_xxxxxzzzz_1, g_z_0_xxxxyyyyy_0, g_z_0_xxxxyyyyy_1, g_z_0_xxxxyyyyz_0, g_z_0_xxxxyyyyz_1, g_z_0_xxxxyyyzz_0, g_z_0_xxxxyyyzz_1, g_z_0_xxxxyyzzz_0, g_z_0_xxxxyyzzz_1, g_z_0_xxxxyzzzz_0, g_z_0_xxxxyzzzz_1, g_z_0_xxxxzzzzz_0, g_z_0_xxxxzzzzz_1, g_z_0_xxxyyyyyy_0, g_z_0_xxxyyyyyy_1, g_z_0_xxxyyyyyz_0, g_z_0_xxxyyyyyz_1, g_z_0_xxxyyyyzz_0, g_z_0_xxxyyyyzz_1, g_z_0_xxxyyyzzz_0, g_z_0_xxxyyyzzz_1, g_z_0_xxxyyzzzz_0, g_z_0_xxxyyzzzz_1, g_z_0_xxxyzzzzz_0, g_z_0_xxxyzzzzz_1, g_z_0_xxxzzzzzz_0, g_z_0_xxxzzzzzz_1, g_z_0_xxyyyyyyy_0, g_z_0_xxyyyyyyy_1, g_z_0_xxyyyyyyz_0, g_z_0_xxyyyyyyz_1, g_z_0_xxyyyyyzz_0, g_z_0_xxyyyyyzz_1, g_z_0_xxyyyyzzz_0, g_z_0_xxyyyyzzz_1, g_z_0_xxyyyzzzz_0, g_z_0_xxyyyzzzz_1, g_z_0_xxyyzzzzz_0, g_z_0_xxyyzzzzz_1, g_z_0_xxyzzzzzz_0, g_z_0_xxyzzzzzz_1, g_z_0_xxzzzzzzz_0, g_z_0_xxzzzzzzz_1, g_z_0_xyyyyyyyy_0, g_z_0_xyyyyyyyy_1, g_z_0_xyyyyyyyz_0, g_z_0_xyyyyyyyz_1, g_z_0_xyyyyyyzz_0, g_z_0_xyyyyyyzz_1, g_z_0_xyyyyyzzz_0, g_z_0_xyyyyyzzz_1, g_z_0_xyyyyzzzz_0, g_z_0_xyyyyzzzz_1, g_z_0_xyyyzzzzz_0, g_z_0_xyyyzzzzz_1, g_z_0_xyyzzzzzz_0, g_z_0_xyyzzzzzz_1, g_z_0_xyzzzzzzz_0, g_z_0_xyzzzzzzz_1, g_z_0_xzzzzzzzz_0, g_z_0_xzzzzzzzz_1, g_z_0_yyyyyyyyy_0, g_z_0_yyyyyyyyy_1, g_z_0_yyyyyyyyz_0, g_z_0_yyyyyyyyz_1, g_z_0_yyyyyyyzz_0, g_z_0_yyyyyyyzz_1, g_z_0_yyyyyyzzz_0, g_z_0_yyyyyyzzz_1, g_z_0_yyyyyzzzz_0, g_z_0_yyyyyzzzz_1, g_z_0_yyyyzzzzz_0, g_z_0_yyyyzzzzz_1, g_z_0_yyyzzzzzz_0, g_z_0_yyyzzzzzz_1, g_z_0_yyzzzzzzz_0, g_z_0_yyzzzzzzz_1, g_z_0_yzzzzzzzz_0, g_z_0_yzzzzzzzz_1, g_z_0_zzzzzzzzz_0, g_z_0_zzzzzzzzz_1, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxxx_1, g_zz_0_xxxxxxxxy_1, g_zz_0_xxxxxxxxz_1, g_zz_0_xxxxxxxy_1, g_zz_0_xxxxxxxyy_1, g_zz_0_xxxxxxxyz_1, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxxzz_1, g_zz_0_xxxxxxyy_1, g_zz_0_xxxxxxyyy_1, g_zz_0_xxxxxxyyz_1, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxyzz_1, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxxzzz_1, g_zz_0_xxxxxyyy_1, g_zz_0_xxxxxyyyy_1, g_zz_0_xxxxxyyyz_1, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyyzz_1, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxyzzz_1, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxxzzzz_1, g_zz_0_xxxxyyyy_1, g_zz_0_xxxxyyyyy_1, g_zz_0_xxxxyyyyz_1, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyyzz_1, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyyzzz_1, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxyzzzz_1, g_zz_0_xxxxzzzz_1, g_zz_0_xxxxzzzzz_1, g_zz_0_xxxyyyyy_1, g_zz_0_xxxyyyyyy_1, g_zz_0_xxxyyyyyz_1, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyyzz_1, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyyzzz_1, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyyzzzz_1, g_zz_0_xxxyzzzz_1, g_zz_0_xxxyzzzzz_1, g_zz_0_xxxzzzzz_1, g_zz_0_xxxzzzzzz_1, g_zz_0_xxyyyyyy_1, g_zz_0_xxyyyyyyy_1, g_zz_0_xxyyyyyyz_1, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyyzz_1, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyyzzz_1, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyyzzzz_1, g_zz_0_xxyyzzzz_1, g_zz_0_xxyyzzzzz_1, g_zz_0_xxyzzzzz_1, g_zz_0_xxyzzzzzz_1, g_zz_0_xxzzzzzz_1, g_zz_0_xxzzzzzzz_1, g_zz_0_xyyyyyyy_1, g_zz_0_xyyyyyyyy_1, g_zz_0_xyyyyyyyz_1, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyyzz_1, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyyzzz_1, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyyzzzz_1, g_zz_0_xyyyzzzz_1, g_zz_0_xyyyzzzzz_1, g_zz_0_xyyzzzzz_1, g_zz_0_xyyzzzzzz_1, g_zz_0_xyzzzzzz_1, g_zz_0_xyzzzzzzz_1, g_zz_0_xzzzzzzz_1, g_zz_0_xzzzzzzzz_1, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyyy_1, g_zz_0_yyyyyyyyz_1, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyyzz_1, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyyzzz_1, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyyzzzz_1, g_zz_0_yyyyzzzz_1, g_zz_0_yyyyzzzzz_1, g_zz_0_yyyzzzzz_1, g_zz_0_yyyzzzzzz_1, g_zz_0_yyzzzzzz_1, g_zz_0_yyzzzzzzz_1, g_zz_0_yzzzzzzz_1, g_zz_0_yzzzzzzzz_1, g_zz_0_zzzzzzzz_1, g_zz_0_zzzzzzzzz_1, g_zzz_0_xxxxxxxxx_0, g_zzz_0_xxxxxxxxy_0, g_zzz_0_xxxxxxxxz_0, g_zzz_0_xxxxxxxyy_0, g_zzz_0_xxxxxxxyz_0, g_zzz_0_xxxxxxxzz_0, g_zzz_0_xxxxxxyyy_0, g_zzz_0_xxxxxxyyz_0, g_zzz_0_xxxxxxyzz_0, g_zzz_0_xxxxxxzzz_0, g_zzz_0_xxxxxyyyy_0, g_zzz_0_xxxxxyyyz_0, g_zzz_0_xxxxxyyzz_0, g_zzz_0_xxxxxyzzz_0, g_zzz_0_xxxxxzzzz_0, g_zzz_0_xxxxyyyyy_0, g_zzz_0_xxxxyyyyz_0, g_zzz_0_xxxxyyyzz_0, g_zzz_0_xxxxyyzzz_0, g_zzz_0_xxxxyzzzz_0, g_zzz_0_xxxxzzzzz_0, g_zzz_0_xxxyyyyyy_0, g_zzz_0_xxxyyyyyz_0, g_zzz_0_xxxyyyyzz_0, g_zzz_0_xxxyyyzzz_0, g_zzz_0_xxxyyzzzz_0, g_zzz_0_xxxyzzzzz_0, g_zzz_0_xxxzzzzzz_0, g_zzz_0_xxyyyyyyy_0, g_zzz_0_xxyyyyyyz_0, g_zzz_0_xxyyyyyzz_0, g_zzz_0_xxyyyyzzz_0, g_zzz_0_xxyyyzzzz_0, g_zzz_0_xxyyzzzzz_0, g_zzz_0_xxyzzzzzz_0, g_zzz_0_xxzzzzzzz_0, g_zzz_0_xyyyyyyyy_0, g_zzz_0_xyyyyyyyz_0, g_zzz_0_xyyyyyyzz_0, g_zzz_0_xyyyyyzzz_0, g_zzz_0_xyyyyzzzz_0, g_zzz_0_xyyyzzzzz_0, g_zzz_0_xyyzzzzzz_0, g_zzz_0_xyzzzzzzz_0, g_zzz_0_xzzzzzzzz_0, g_zzz_0_yyyyyyyyy_0, g_zzz_0_yyyyyyyyz_0, g_zzz_0_yyyyyyyzz_0, g_zzz_0_yyyyyyzzz_0, g_zzz_0_yyyyyzzzz_0, g_zzz_0_yyyyzzzzz_0, g_zzz_0_yyyzzzzzz_0, g_zzz_0_yyzzzzzzz_0, g_zzz_0_yzzzzzzzz_0, g_zzz_0_zzzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xxxxxxxxx_0[i] = 2.0 * g_z_0_xxxxxxxxx_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxxx_1[i] * fz_be_0 + g_zz_0_xxxxxxxxx_1[i] * wa_z[i];

        g_zzz_0_xxxxxxxxy_0[i] = 2.0 * g_z_0_xxxxxxxxy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxxy_1[i] * fz_be_0 + g_zz_0_xxxxxxxxy_1[i] * wa_z[i];

        g_zzz_0_xxxxxxxxz_0[i] = 2.0 * g_z_0_xxxxxxxxz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxxz_1[i] * fz_be_0 + g_zz_0_xxxxxxxx_1[i] * fi_acd_0 + g_zz_0_xxxxxxxxz_1[i] * wa_z[i];

        g_zzz_0_xxxxxxxyy_0[i] = 2.0 * g_z_0_xxxxxxxyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxyy_1[i] * fz_be_0 + g_zz_0_xxxxxxxyy_1[i] * wa_z[i];

        g_zzz_0_xxxxxxxyz_0[i] = 2.0 * g_z_0_xxxxxxxyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxyz_1[i] * fz_be_0 + g_zz_0_xxxxxxxy_1[i] * fi_acd_0 + g_zz_0_xxxxxxxyz_1[i] * wa_z[i];

        g_zzz_0_xxxxxxxzz_0[i] = 2.0 * g_z_0_xxxxxxxzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxxzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxxxxz_1[i] * fi_acd_0 + g_zz_0_xxxxxxxzz_1[i] * wa_z[i];

        g_zzz_0_xxxxxxyyy_0[i] = 2.0 * g_z_0_xxxxxxyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxyyy_1[i] * fz_be_0 + g_zz_0_xxxxxxyyy_1[i] * wa_z[i];

        g_zzz_0_xxxxxxyyz_0[i] = 2.0 * g_z_0_xxxxxxyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxyyz_1[i] * fz_be_0 + g_zz_0_xxxxxxyy_1[i] * fi_acd_0 + g_zz_0_xxxxxxyyz_1[i] * wa_z[i];

        g_zzz_0_xxxxxxyzz_0[i] = 2.0 * g_z_0_xxxxxxyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxxxyz_1[i] * fi_acd_0 + g_zz_0_xxxxxxyzz_1[i] * wa_z[i];

        g_zzz_0_xxxxxxzzz_0[i] = 2.0 * g_z_0_xxxxxxzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxxzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxxxxzz_1[i] * fi_acd_0 + g_zz_0_xxxxxxzzz_1[i] * wa_z[i];

        g_zzz_0_xxxxxyyyy_0[i] = 2.0 * g_z_0_xxxxxyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyyyy_1[i] * fz_be_0 + g_zz_0_xxxxxyyyy_1[i] * wa_z[i];

        g_zzz_0_xxxxxyyyz_0[i] = 2.0 * g_z_0_xxxxxyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyyyz_1[i] * fz_be_0 + g_zz_0_xxxxxyyy_1[i] * fi_acd_0 + g_zz_0_xxxxxyyyz_1[i] * wa_z[i];

        g_zzz_0_xxxxxyyzz_0[i] = 2.0 * g_z_0_xxxxxyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxxyyz_1[i] * fi_acd_0 + g_zz_0_xxxxxyyzz_1[i] * wa_z[i];

        g_zzz_0_xxxxxyzzz_0[i] = 2.0 * g_z_0_xxxxxyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxxxyzz_1[i] * fi_acd_0 + g_zz_0_xxxxxyzzz_1[i] * wa_z[i];

        g_zzz_0_xxxxxzzzz_0[i] = 2.0 * g_z_0_xxxxxzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxxzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxxxxzzz_1[i] * fi_acd_0 + g_zz_0_xxxxxzzzz_1[i] * wa_z[i];

        g_zzz_0_xxxxyyyyy_0[i] = 2.0 * g_z_0_xxxxyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyyyy_1[i] * fz_be_0 + g_zz_0_xxxxyyyyy_1[i] * wa_z[i];

        g_zzz_0_xxxxyyyyz_0[i] = 2.0 * g_z_0_xxxxyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyyyz_1[i] * fz_be_0 + g_zz_0_xxxxyyyy_1[i] * fi_acd_0 + g_zz_0_xxxxyyyyz_1[i] * wa_z[i];

        g_zzz_0_xxxxyyyzz_0[i] = 2.0 * g_z_0_xxxxyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxxyyyz_1[i] * fi_acd_0 + g_zz_0_xxxxyyyzz_1[i] * wa_z[i];

        g_zzz_0_xxxxyyzzz_0[i] = 2.0 * g_z_0_xxxxyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxxyyzz_1[i] * fi_acd_0 + g_zz_0_xxxxyyzzz_1[i] * wa_z[i];

        g_zzz_0_xxxxyzzzz_0[i] = 2.0 * g_z_0_xxxxyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxxxyzzz_1[i] * fi_acd_0 + g_zz_0_xxxxyzzzz_1[i] * wa_z[i];

        g_zzz_0_xxxxzzzzz_0[i] = 2.0 * g_z_0_xxxxzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xxxxzzzz_1[i] * fi_acd_0 + g_zz_0_xxxxzzzzz_1[i] * wa_z[i];

        g_zzz_0_xxxyyyyyy_0[i] = 2.0 * g_z_0_xxxyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyyyy_1[i] * fz_be_0 + g_zz_0_xxxyyyyyy_1[i] * wa_z[i];

        g_zzz_0_xxxyyyyyz_0[i] = 2.0 * g_z_0_xxxyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyyyz_1[i] * fz_be_0 + g_zz_0_xxxyyyyy_1[i] * fi_acd_0 + g_zz_0_xxxyyyyyz_1[i] * wa_z[i];

        g_zzz_0_xxxyyyyzz_0[i] = 2.0 * g_z_0_xxxyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxyyyyz_1[i] * fi_acd_0 + g_zz_0_xxxyyyyzz_1[i] * wa_z[i];

        g_zzz_0_xxxyyyzzz_0[i] = 2.0 * g_z_0_xxxyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxxyyyzz_1[i] * fi_acd_0 + g_zz_0_xxxyyyzzz_1[i] * wa_z[i];

        g_zzz_0_xxxyyzzzz_0[i] = 2.0 * g_z_0_xxxyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxxyyzzz_1[i] * fi_acd_0 + g_zz_0_xxxyyzzzz_1[i] * wa_z[i];

        g_zzz_0_xxxyzzzzz_0[i] = 2.0 * g_z_0_xxxyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xxxyzzzz_1[i] * fi_acd_0 + g_zz_0_xxxyzzzzz_1[i] * wa_z[i];

        g_zzz_0_xxxzzzzzz_0[i] = 2.0 * g_z_0_xxxzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_xxxzzzzz_1[i] * fi_acd_0 + g_zz_0_xxxzzzzzz_1[i] * wa_z[i];

        g_zzz_0_xxyyyyyyy_0[i] = 2.0 * g_z_0_xxyyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyyyy_1[i] * fz_be_0 + g_zz_0_xxyyyyyyy_1[i] * wa_z[i];

        g_zzz_0_xxyyyyyyz_0[i] = 2.0 * g_z_0_xxyyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyyyz_1[i] * fz_be_0 + g_zz_0_xxyyyyyy_1[i] * fi_acd_0 + g_zz_0_xxyyyyyyz_1[i] * wa_z[i];

        g_zzz_0_xxyyyyyzz_0[i] = 2.0 * g_z_0_xxyyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxyyyyyz_1[i] * fi_acd_0 + g_zz_0_xxyyyyyzz_1[i] * wa_z[i];

        g_zzz_0_xxyyyyzzz_0[i] = 2.0 * g_z_0_xxyyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxyyyyzz_1[i] * fi_acd_0 + g_zz_0_xxyyyyzzz_1[i] * wa_z[i];

        g_zzz_0_xxyyyzzzz_0[i] = 2.0 * g_z_0_xxyyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xxyyyzzz_1[i] * fi_acd_0 + g_zz_0_xxyyyzzzz_1[i] * wa_z[i];

        g_zzz_0_xxyyzzzzz_0[i] = 2.0 * g_z_0_xxyyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xxyyzzzz_1[i] * fi_acd_0 + g_zz_0_xxyyzzzzz_1[i] * wa_z[i];

        g_zzz_0_xxyzzzzzz_0[i] = 2.0 * g_z_0_xxyzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_xxyzzzzz_1[i] * fi_acd_0 + g_zz_0_xxyzzzzzz_1[i] * wa_z[i];

        g_zzz_0_xxzzzzzzz_0[i] = 2.0 * g_z_0_xxzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zz_0_xxzzzzzz_1[i] * fi_acd_0 + g_zz_0_xxzzzzzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyyyyyy_0[i] = 2.0 * g_z_0_xyyyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyyyy_1[i] * fz_be_0 + g_zz_0_xyyyyyyyy_1[i] * wa_z[i];

        g_zzz_0_xyyyyyyyz_0[i] = 2.0 * g_z_0_xyyyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyyyz_1[i] * fz_be_0 + g_zz_0_xyyyyyyy_1[i] * fi_acd_0 + g_zz_0_xyyyyyyyz_1[i] * wa_z[i];

        g_zzz_0_xyyyyyyzz_0[i] = 2.0 * g_z_0_xyyyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xyyyyyyz_1[i] * fi_acd_0 + g_zz_0_xyyyyyyzz_1[i] * wa_z[i];

        g_zzz_0_xyyyyyzzz_0[i] = 2.0 * g_z_0_xyyyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xyyyyyzz_1[i] * fi_acd_0 + g_zz_0_xyyyyyzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyyzzzz_0[i] = 2.0 * g_z_0_xyyyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xyyyyzzz_1[i] * fi_acd_0 + g_zz_0_xyyyyzzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyzzzzz_0[i] = 2.0 * g_z_0_xyyyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_xyyyzzzz_1[i] * fi_acd_0 + g_zz_0_xyyyzzzzz_1[i] * wa_z[i];

        g_zzz_0_xyyzzzzzz_0[i] = 2.0 * g_z_0_xyyzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_xyyzzzzz_1[i] * fi_acd_0 + g_zz_0_xyyzzzzzz_1[i] * wa_z[i];

        g_zzz_0_xyzzzzzzz_0[i] = 2.0 * g_z_0_xyzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zz_0_xyzzzzzz_1[i] * fi_acd_0 + g_zz_0_xyzzzzzzz_1[i] * wa_z[i];

        g_zzz_0_xzzzzzzzz_0[i] = 2.0 * g_z_0_xzzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xzzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zz_0_xzzzzzzz_1[i] * fi_acd_0 + g_zz_0_xzzzzzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyyyyyy_0[i] = 2.0 * g_z_0_yyyyyyyyy_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyyyy_1[i] * fz_be_0 + g_zz_0_yyyyyyyyy_1[i] * wa_z[i];

        g_zzz_0_yyyyyyyyz_0[i] = 2.0 * g_z_0_yyyyyyyyz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyyyz_1[i] * fz_be_0 + g_zz_0_yyyyyyyy_1[i] * fi_acd_0 + g_zz_0_yyyyyyyyz_1[i] * wa_z[i];

        g_zzz_0_yyyyyyyzz_0[i] = 2.0 * g_z_0_yyyyyyyzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_yyyyyyyz_1[i] * fi_acd_0 + g_zz_0_yyyyyyyzz_1[i] * wa_z[i];

        g_zzz_0_yyyyyyzzz_0[i] = 2.0 * g_z_0_yyyyyyzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_yyyyyyzz_1[i] * fi_acd_0 + g_zz_0_yyyyyyzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyyzzzz_0[i] = 2.0 * g_z_0_yyyyyzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_yyyyyzzz_1[i] * fi_acd_0 + g_zz_0_yyyyyzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyzzzzz_0[i] = 2.0 * g_z_0_yyyyzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_yyyyzzzz_1[i] * fi_acd_0 + g_zz_0_yyyyzzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyzzzzzz_0[i] = 2.0 * g_z_0_yyyzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_0_yyyzzzzz_1[i] * fi_acd_0 + g_zz_0_yyyzzzzzz_1[i] * wa_z[i];

        g_zzz_0_yyzzzzzzz_0[i] = 2.0 * g_z_0_yyzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zz_0_yyzzzzzz_1[i] * fi_acd_0 + g_zz_0_yyzzzzzzz_1[i] * wa_z[i];

        g_zzz_0_yzzzzzzzz_0[i] = 2.0 * g_z_0_yzzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yzzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zz_0_yzzzzzzz_1[i] * fi_acd_0 + g_zz_0_yzzzzzzzz_1[i] * wa_z[i];

        g_zzz_0_zzzzzzzzz_0[i] = 2.0 * g_z_0_zzzzzzzzz_0[i] * fbe_0 - 2.0 * g_z_0_zzzzzzzzz_1[i] * fz_be_0 + 9.0 * g_zz_0_zzzzzzzz_1[i] * fi_acd_0 + g_zz_0_zzzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

