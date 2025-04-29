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

#include "TwoCenterElectronRepulsionPrimRecIG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ig(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ig,
                                const size_t idx_eri_0_gg,
                                const size_t idx_eri_1_gg,
                                const size_t idx_eri_1_hf,
                                const size_t idx_eri_1_hg,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GG

    auto g_xxxx_xxxx_0 = pbuffer.data(idx_eri_0_gg);

    auto g_xxxx_xxxy_0 = pbuffer.data(idx_eri_0_gg + 1);

    auto g_xxxx_xxxz_0 = pbuffer.data(idx_eri_0_gg + 2);

    auto g_xxxx_xxyy_0 = pbuffer.data(idx_eri_0_gg + 3);

    auto g_xxxx_xxyz_0 = pbuffer.data(idx_eri_0_gg + 4);

    auto g_xxxx_xxzz_0 = pbuffer.data(idx_eri_0_gg + 5);

    auto g_xxxx_xyyy_0 = pbuffer.data(idx_eri_0_gg + 6);

    auto g_xxxx_xyyz_0 = pbuffer.data(idx_eri_0_gg + 7);

    auto g_xxxx_xyzz_0 = pbuffer.data(idx_eri_0_gg + 8);

    auto g_xxxx_xzzz_0 = pbuffer.data(idx_eri_0_gg + 9);

    auto g_xxxx_yyyy_0 = pbuffer.data(idx_eri_0_gg + 10);

    auto g_xxxx_yyyz_0 = pbuffer.data(idx_eri_0_gg + 11);

    auto g_xxxx_yyzz_0 = pbuffer.data(idx_eri_0_gg + 12);

    auto g_xxxx_yzzz_0 = pbuffer.data(idx_eri_0_gg + 13);

    auto g_xxxx_zzzz_0 = pbuffer.data(idx_eri_0_gg + 14);

    auto g_xxxy_xxxx_0 = pbuffer.data(idx_eri_0_gg + 15);

    auto g_xxxy_xxxz_0 = pbuffer.data(idx_eri_0_gg + 17);

    auto g_xxxy_xxzz_0 = pbuffer.data(idx_eri_0_gg + 20);

    auto g_xxxy_xzzz_0 = pbuffer.data(idx_eri_0_gg + 24);

    auto g_xxxz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 30);

    auto g_xxxz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 31);

    auto g_xxxz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 33);

    auto g_xxxz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 36);

    auto g_xxyy_xxxx_0 = pbuffer.data(idx_eri_0_gg + 45);

    auto g_xxyy_xxxy_0 = pbuffer.data(idx_eri_0_gg + 46);

    auto g_xxyy_xxxz_0 = pbuffer.data(idx_eri_0_gg + 47);

    auto g_xxyy_xxyy_0 = pbuffer.data(idx_eri_0_gg + 48);

    auto g_xxyy_xxyz_0 = pbuffer.data(idx_eri_0_gg + 49);

    auto g_xxyy_xxzz_0 = pbuffer.data(idx_eri_0_gg + 50);

    auto g_xxyy_xyyy_0 = pbuffer.data(idx_eri_0_gg + 51);

    auto g_xxyy_xyyz_0 = pbuffer.data(idx_eri_0_gg + 52);

    auto g_xxyy_xyzz_0 = pbuffer.data(idx_eri_0_gg + 53);

    auto g_xxyy_xzzz_0 = pbuffer.data(idx_eri_0_gg + 54);

    auto g_xxyy_yyyy_0 = pbuffer.data(idx_eri_0_gg + 55);

    auto g_xxyy_yyyz_0 = pbuffer.data(idx_eri_0_gg + 56);

    auto g_xxyy_yyzz_0 = pbuffer.data(idx_eri_0_gg + 57);

    auto g_xxyy_yzzz_0 = pbuffer.data(idx_eri_0_gg + 58);

    auto g_xxyy_zzzz_0 = pbuffer.data(idx_eri_0_gg + 59);

    auto g_xxzz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 75);

    auto g_xxzz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 76);

    auto g_xxzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 77);

    auto g_xxzz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 78);

    auto g_xxzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 79);

    auto g_xxzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 80);

    auto g_xxzz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 81);

    auto g_xxzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 82);

    auto g_xxzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 83);

    auto g_xxzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 84);

    auto g_xxzz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 85);

    auto g_xxzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 86);

    auto g_xxzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 87);

    auto g_xxzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 88);

    auto g_xxzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 89);

    auto g_xyyy_xxxy_0 = pbuffer.data(idx_eri_0_gg + 91);

    auto g_xyyy_xxyy_0 = pbuffer.data(idx_eri_0_gg + 93);

    auto g_xyyy_xxyz_0 = pbuffer.data(idx_eri_0_gg + 94);

    auto g_xyyy_xyyy_0 = pbuffer.data(idx_eri_0_gg + 96);

    auto g_xyyy_xyyz_0 = pbuffer.data(idx_eri_0_gg + 97);

    auto g_xyyy_xyzz_0 = pbuffer.data(idx_eri_0_gg + 98);

    auto g_xyyy_yyyy_0 = pbuffer.data(idx_eri_0_gg + 100);

    auto g_xyyy_yyyz_0 = pbuffer.data(idx_eri_0_gg + 101);

    auto g_xyyy_yyzz_0 = pbuffer.data(idx_eri_0_gg + 102);

    auto g_xyyy_yzzz_0 = pbuffer.data(idx_eri_0_gg + 103);

    auto g_xyyy_zzzz_0 = pbuffer.data(idx_eri_0_gg + 104);

    auto g_xzzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 137);

    auto g_xzzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 139);

    auto g_xzzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 140);

    auto g_xzzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 142);

    auto g_xzzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 143);

    auto g_xzzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 144);

    auto g_xzzz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 145);

    auto g_xzzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 146);

    auto g_xzzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 147);

    auto g_xzzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 148);

    auto g_xzzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 149);

    auto g_yyyy_xxxx_0 = pbuffer.data(idx_eri_0_gg + 150);

    auto g_yyyy_xxxy_0 = pbuffer.data(idx_eri_0_gg + 151);

    auto g_yyyy_xxxz_0 = pbuffer.data(idx_eri_0_gg + 152);

    auto g_yyyy_xxyy_0 = pbuffer.data(idx_eri_0_gg + 153);

    auto g_yyyy_xxyz_0 = pbuffer.data(idx_eri_0_gg + 154);

    auto g_yyyy_xxzz_0 = pbuffer.data(idx_eri_0_gg + 155);

    auto g_yyyy_xyyy_0 = pbuffer.data(idx_eri_0_gg + 156);

    auto g_yyyy_xyyz_0 = pbuffer.data(idx_eri_0_gg + 157);

    auto g_yyyy_xyzz_0 = pbuffer.data(idx_eri_0_gg + 158);

    auto g_yyyy_xzzz_0 = pbuffer.data(idx_eri_0_gg + 159);

    auto g_yyyy_yyyy_0 = pbuffer.data(idx_eri_0_gg + 160);

    auto g_yyyy_yyyz_0 = pbuffer.data(idx_eri_0_gg + 161);

    auto g_yyyy_yyzz_0 = pbuffer.data(idx_eri_0_gg + 162);

    auto g_yyyy_yzzz_0 = pbuffer.data(idx_eri_0_gg + 163);

    auto g_yyyy_zzzz_0 = pbuffer.data(idx_eri_0_gg + 164);

    auto g_yyyz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 166);

    auto g_yyyz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 168);

    auto g_yyyz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 171);

    auto g_yyyz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 175);

    auto g_yyzz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 180);

    auto g_yyzz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 181);

    auto g_yyzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 182);

    auto g_yyzz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 183);

    auto g_yyzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 184);

    auto g_yyzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 185);

    auto g_yyzz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 186);

    auto g_yyzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 187);

    auto g_yyzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 188);

    auto g_yyzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 189);

    auto g_yyzz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 190);

    auto g_yyzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 191);

    auto g_yyzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 192);

    auto g_yyzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 193);

    auto g_yyzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 194);

    auto g_yzzz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 195);

    auto g_yzzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 197);

    auto g_yzzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 199);

    auto g_yzzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 200);

    auto g_yzzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 202);

    auto g_yzzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 203);

    auto g_yzzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 204);

    auto g_yzzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 206);

    auto g_yzzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 207);

    auto g_yzzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 208);

    auto g_yzzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 209);

    auto g_zzzz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 210);

    auto g_zzzz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 211);

    auto g_zzzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 212);

    auto g_zzzz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 213);

    auto g_zzzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 214);

    auto g_zzzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 215);

    auto g_zzzz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 216);

    auto g_zzzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 217);

    auto g_zzzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 218);

    auto g_zzzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 219);

    auto g_zzzz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 220);

    auto g_zzzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 221);

    auto g_zzzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 222);

    auto g_zzzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 223);

    auto g_zzzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 224);

    // Set up components of auxiliary buffer : GG

    auto g_xxxx_xxxx_1 = pbuffer.data(idx_eri_1_gg);

    auto g_xxxx_xxxy_1 = pbuffer.data(idx_eri_1_gg + 1);

    auto g_xxxx_xxxz_1 = pbuffer.data(idx_eri_1_gg + 2);

    auto g_xxxx_xxyy_1 = pbuffer.data(idx_eri_1_gg + 3);

    auto g_xxxx_xxyz_1 = pbuffer.data(idx_eri_1_gg + 4);

    auto g_xxxx_xxzz_1 = pbuffer.data(idx_eri_1_gg + 5);

    auto g_xxxx_xyyy_1 = pbuffer.data(idx_eri_1_gg + 6);

    auto g_xxxx_xyyz_1 = pbuffer.data(idx_eri_1_gg + 7);

    auto g_xxxx_xyzz_1 = pbuffer.data(idx_eri_1_gg + 8);

    auto g_xxxx_xzzz_1 = pbuffer.data(idx_eri_1_gg + 9);

    auto g_xxxx_yyyy_1 = pbuffer.data(idx_eri_1_gg + 10);

    auto g_xxxx_yyyz_1 = pbuffer.data(idx_eri_1_gg + 11);

    auto g_xxxx_yyzz_1 = pbuffer.data(idx_eri_1_gg + 12);

    auto g_xxxx_yzzz_1 = pbuffer.data(idx_eri_1_gg + 13);

    auto g_xxxx_zzzz_1 = pbuffer.data(idx_eri_1_gg + 14);

    auto g_xxxy_xxxx_1 = pbuffer.data(idx_eri_1_gg + 15);

    auto g_xxxy_xxxz_1 = pbuffer.data(idx_eri_1_gg + 17);

    auto g_xxxy_xxzz_1 = pbuffer.data(idx_eri_1_gg + 20);

    auto g_xxxy_xzzz_1 = pbuffer.data(idx_eri_1_gg + 24);

    auto g_xxxz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 30);

    auto g_xxxz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 31);

    auto g_xxxz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 33);

    auto g_xxxz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 36);

    auto g_xxyy_xxxx_1 = pbuffer.data(idx_eri_1_gg + 45);

    auto g_xxyy_xxxy_1 = pbuffer.data(idx_eri_1_gg + 46);

    auto g_xxyy_xxxz_1 = pbuffer.data(idx_eri_1_gg + 47);

    auto g_xxyy_xxyy_1 = pbuffer.data(idx_eri_1_gg + 48);

    auto g_xxyy_xxyz_1 = pbuffer.data(idx_eri_1_gg + 49);

    auto g_xxyy_xxzz_1 = pbuffer.data(idx_eri_1_gg + 50);

    auto g_xxyy_xyyy_1 = pbuffer.data(idx_eri_1_gg + 51);

    auto g_xxyy_xyyz_1 = pbuffer.data(idx_eri_1_gg + 52);

    auto g_xxyy_xyzz_1 = pbuffer.data(idx_eri_1_gg + 53);

    auto g_xxyy_xzzz_1 = pbuffer.data(idx_eri_1_gg + 54);

    auto g_xxyy_yyyy_1 = pbuffer.data(idx_eri_1_gg + 55);

    auto g_xxyy_yyyz_1 = pbuffer.data(idx_eri_1_gg + 56);

    auto g_xxyy_yyzz_1 = pbuffer.data(idx_eri_1_gg + 57);

    auto g_xxyy_yzzz_1 = pbuffer.data(idx_eri_1_gg + 58);

    auto g_xxyy_zzzz_1 = pbuffer.data(idx_eri_1_gg + 59);

    auto g_xxzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 75);

    auto g_xxzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 76);

    auto g_xxzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 77);

    auto g_xxzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 78);

    auto g_xxzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 79);

    auto g_xxzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 80);

    auto g_xxzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 81);

    auto g_xxzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 82);

    auto g_xxzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 83);

    auto g_xxzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 84);

    auto g_xxzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 85);

    auto g_xxzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 86);

    auto g_xxzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 87);

    auto g_xxzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 88);

    auto g_xxzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 89);

    auto g_xyyy_xxxy_1 = pbuffer.data(idx_eri_1_gg + 91);

    auto g_xyyy_xxyy_1 = pbuffer.data(idx_eri_1_gg + 93);

    auto g_xyyy_xxyz_1 = pbuffer.data(idx_eri_1_gg + 94);

    auto g_xyyy_xyyy_1 = pbuffer.data(idx_eri_1_gg + 96);

    auto g_xyyy_xyyz_1 = pbuffer.data(idx_eri_1_gg + 97);

    auto g_xyyy_xyzz_1 = pbuffer.data(idx_eri_1_gg + 98);

    auto g_xyyy_yyyy_1 = pbuffer.data(idx_eri_1_gg + 100);

    auto g_xyyy_yyyz_1 = pbuffer.data(idx_eri_1_gg + 101);

    auto g_xyyy_yyzz_1 = pbuffer.data(idx_eri_1_gg + 102);

    auto g_xyyy_yzzz_1 = pbuffer.data(idx_eri_1_gg + 103);

    auto g_xyyy_zzzz_1 = pbuffer.data(idx_eri_1_gg + 104);

    auto g_xzzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 137);

    auto g_xzzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 139);

    auto g_xzzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 140);

    auto g_xzzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 142);

    auto g_xzzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 143);

    auto g_xzzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 144);

    auto g_xzzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 145);

    auto g_xzzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 146);

    auto g_xzzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 147);

    auto g_xzzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 148);

    auto g_xzzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 149);

    auto g_yyyy_xxxx_1 = pbuffer.data(idx_eri_1_gg + 150);

    auto g_yyyy_xxxy_1 = pbuffer.data(idx_eri_1_gg + 151);

    auto g_yyyy_xxxz_1 = pbuffer.data(idx_eri_1_gg + 152);

    auto g_yyyy_xxyy_1 = pbuffer.data(idx_eri_1_gg + 153);

    auto g_yyyy_xxyz_1 = pbuffer.data(idx_eri_1_gg + 154);

    auto g_yyyy_xxzz_1 = pbuffer.data(idx_eri_1_gg + 155);

    auto g_yyyy_xyyy_1 = pbuffer.data(idx_eri_1_gg + 156);

    auto g_yyyy_xyyz_1 = pbuffer.data(idx_eri_1_gg + 157);

    auto g_yyyy_xyzz_1 = pbuffer.data(idx_eri_1_gg + 158);

    auto g_yyyy_xzzz_1 = pbuffer.data(idx_eri_1_gg + 159);

    auto g_yyyy_yyyy_1 = pbuffer.data(idx_eri_1_gg + 160);

    auto g_yyyy_yyyz_1 = pbuffer.data(idx_eri_1_gg + 161);

    auto g_yyyy_yyzz_1 = pbuffer.data(idx_eri_1_gg + 162);

    auto g_yyyy_yzzz_1 = pbuffer.data(idx_eri_1_gg + 163);

    auto g_yyyy_zzzz_1 = pbuffer.data(idx_eri_1_gg + 164);

    auto g_yyyz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 166);

    auto g_yyyz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 168);

    auto g_yyyz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 171);

    auto g_yyyz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 175);

    auto g_yyzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 180);

    auto g_yyzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 181);

    auto g_yyzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 182);

    auto g_yyzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 183);

    auto g_yyzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 184);

    auto g_yyzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 185);

    auto g_yyzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 186);

    auto g_yyzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 187);

    auto g_yyzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 188);

    auto g_yyzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 189);

    auto g_yyzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 190);

    auto g_yyzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 191);

    auto g_yyzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 192);

    auto g_yyzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 193);

    auto g_yyzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 194);

    auto g_yzzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 195);

    auto g_yzzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 197);

    auto g_yzzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 199);

    auto g_yzzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 200);

    auto g_yzzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 202);

    auto g_yzzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 203);

    auto g_yzzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 204);

    auto g_yzzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 206);

    auto g_yzzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 207);

    auto g_yzzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 208);

    auto g_yzzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 209);

    auto g_zzzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 210);

    auto g_zzzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 211);

    auto g_zzzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 212);

    auto g_zzzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 213);

    auto g_zzzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 214);

    auto g_zzzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 215);

    auto g_zzzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 216);

    auto g_zzzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 217);

    auto g_zzzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 218);

    auto g_zzzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 219);

    auto g_zzzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 220);

    auto g_zzzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 221);

    auto g_zzzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 222);

    auto g_zzzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 223);

    auto g_zzzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 224);

    // Set up components of auxiliary buffer : HF

    auto g_xxxxx_xxx_1 = pbuffer.data(idx_eri_1_hf);

    auto g_xxxxx_xxy_1 = pbuffer.data(idx_eri_1_hf + 1);

    auto g_xxxxx_xxz_1 = pbuffer.data(idx_eri_1_hf + 2);

    auto g_xxxxx_xyy_1 = pbuffer.data(idx_eri_1_hf + 3);

    auto g_xxxxx_xyz_1 = pbuffer.data(idx_eri_1_hf + 4);

    auto g_xxxxx_xzz_1 = pbuffer.data(idx_eri_1_hf + 5);

    auto g_xxxxx_yyy_1 = pbuffer.data(idx_eri_1_hf + 6);

    auto g_xxxxx_yyz_1 = pbuffer.data(idx_eri_1_hf + 7);

    auto g_xxxxx_yzz_1 = pbuffer.data(idx_eri_1_hf + 8);

    auto g_xxxxx_zzz_1 = pbuffer.data(idx_eri_1_hf + 9);

    auto g_xxxxz_xxz_1 = pbuffer.data(idx_eri_1_hf + 22);

    auto g_xxxxz_xyz_1 = pbuffer.data(idx_eri_1_hf + 24);

    auto g_xxxxz_xzz_1 = pbuffer.data(idx_eri_1_hf + 25);

    auto g_xxxxz_yyz_1 = pbuffer.data(idx_eri_1_hf + 27);

    auto g_xxxxz_yzz_1 = pbuffer.data(idx_eri_1_hf + 28);

    auto g_xxxxz_zzz_1 = pbuffer.data(idx_eri_1_hf + 29);

    auto g_xxxyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 30);

    auto g_xxxyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 31);

    auto g_xxxyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 32);

    auto g_xxxyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 33);

    auto g_xxxyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 34);

    auto g_xxxyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 35);

    auto g_xxxyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 36);

    auto g_xxxyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 37);

    auto g_xxxyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 38);

    auto g_xxxyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 39);

    auto g_xxxzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 50);

    auto g_xxxzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 51);

    auto g_xxxzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 52);

    auto g_xxxzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 53);

    auto g_xxxzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 54);

    auto g_xxxzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 55);

    auto g_xxxzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 56);

    auto g_xxxzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 57);

    auto g_xxxzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 58);

    auto g_xxxzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 59);

    auto g_xxyyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 60);

    auto g_xxyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 61);

    auto g_xxyyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 62);

    auto g_xxyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 63);

    auto g_xxyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 64);

    auto g_xxyyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 65);

    auto g_xxyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 66);

    auto g_xxyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 67);

    auto g_xxyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 68);

    auto g_xxyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 69);

    auto g_xxzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 90);

    auto g_xxzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 91);

    auto g_xxzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 92);

    auto g_xxzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 93);

    auto g_xxzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 94);

    auto g_xxzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 95);

    auto g_xxzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 96);

    auto g_xxzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 97);

    auto g_xxzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 98);

    auto g_xxzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 99);

    auto g_xyyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 101);

    auto g_xyyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 103);

    auto g_xyyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 104);

    auto g_xyyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 106);

    auto g_xyyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 107);

    auto g_xyyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 108);

    auto g_xyyzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 124);

    auto g_xyyzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 127);

    auto g_xyyzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 128);

    auto g_xzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 142);

    auto g_xzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 144);

    auto g_xzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 145);

    auto g_xzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 147);

    auto g_xzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 148);

    auto g_xzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 149);

    auto g_yyyyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 150);

    auto g_yyyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 151);

    auto g_yyyyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 152);

    auto g_yyyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 153);

    auto g_yyyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 154);

    auto g_yyyyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 155);

    auto g_yyyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 156);

    auto g_yyyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 157);

    auto g_yyyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 158);

    auto g_yyyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 159);

    auto g_yyyyz_xxz_1 = pbuffer.data(idx_eri_1_hf + 162);

    auto g_yyyyz_xyz_1 = pbuffer.data(idx_eri_1_hf + 164);

    auto g_yyyyz_xzz_1 = pbuffer.data(idx_eri_1_hf + 165);

    auto g_yyyyz_yyz_1 = pbuffer.data(idx_eri_1_hf + 167);

    auto g_yyyyz_yzz_1 = pbuffer.data(idx_eri_1_hf + 168);

    auto g_yyyyz_zzz_1 = pbuffer.data(idx_eri_1_hf + 169);

    auto g_yyyzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 170);

    auto g_yyyzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 171);

    auto g_yyyzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 172);

    auto g_yyyzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 173);

    auto g_yyyzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 174);

    auto g_yyyzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 175);

    auto g_yyyzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 176);

    auto g_yyyzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 177);

    auto g_yyyzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 178);

    auto g_yyyzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 179);

    auto g_yyzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 180);

    auto g_yyzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 181);

    auto g_yyzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 182);

    auto g_yyzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 183);

    auto g_yyzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 184);

    auto g_yyzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 185);

    auto g_yyzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 186);

    auto g_yyzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 187);

    auto g_yyzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 188);

    auto g_yyzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 189);

    auto g_yzzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 191);

    auto g_yzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 192);

    auto g_yzzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 193);

    auto g_yzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 194);

    auto g_yzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 195);

    auto g_yzzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 196);

    auto g_yzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 197);

    auto g_yzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 198);

    auto g_yzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 199);

    auto g_zzzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 200);

    auto g_zzzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 201);

    auto g_zzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 202);

    auto g_zzzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 203);

    auto g_zzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 204);

    auto g_zzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 205);

    auto g_zzzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 206);

    auto g_zzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 207);

    auto g_zzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 208);

    auto g_zzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 209);

    // Set up components of auxiliary buffer : HG

    auto g_xxxxx_xxxx_1 = pbuffer.data(idx_eri_1_hg);

    auto g_xxxxx_xxxy_1 = pbuffer.data(idx_eri_1_hg + 1);

    auto g_xxxxx_xxxz_1 = pbuffer.data(idx_eri_1_hg + 2);

    auto g_xxxxx_xxyy_1 = pbuffer.data(idx_eri_1_hg + 3);

    auto g_xxxxx_xxyz_1 = pbuffer.data(idx_eri_1_hg + 4);

    auto g_xxxxx_xxzz_1 = pbuffer.data(idx_eri_1_hg + 5);

    auto g_xxxxx_xyyy_1 = pbuffer.data(idx_eri_1_hg + 6);

    auto g_xxxxx_xyyz_1 = pbuffer.data(idx_eri_1_hg + 7);

    auto g_xxxxx_xyzz_1 = pbuffer.data(idx_eri_1_hg + 8);

    auto g_xxxxx_xzzz_1 = pbuffer.data(idx_eri_1_hg + 9);

    auto g_xxxxx_yyyy_1 = pbuffer.data(idx_eri_1_hg + 10);

    auto g_xxxxx_yyyz_1 = pbuffer.data(idx_eri_1_hg + 11);

    auto g_xxxxx_yyzz_1 = pbuffer.data(idx_eri_1_hg + 12);

    auto g_xxxxx_yzzz_1 = pbuffer.data(idx_eri_1_hg + 13);

    auto g_xxxxx_zzzz_1 = pbuffer.data(idx_eri_1_hg + 14);

    auto g_xxxxy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 15);

    auto g_xxxxy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 16);

    auto g_xxxxy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 17);

    auto g_xxxxy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 18);

    auto g_xxxxy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 20);

    auto g_xxxxy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 21);

    auto g_xxxxy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 24);

    auto g_xxxxy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 25);

    auto g_xxxxz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 30);

    auto g_xxxxz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 31);

    auto g_xxxxz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 32);

    auto g_xxxxz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 33);

    auto g_xxxxz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 34);

    auto g_xxxxz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 35);

    auto g_xxxxz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 36);

    auto g_xxxxz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 37);

    auto g_xxxxz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 38);

    auto g_xxxxz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 39);

    auto g_xxxxz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 41);

    auto g_xxxxz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 42);

    auto g_xxxxz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 43);

    auto g_xxxxz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 44);

    auto g_xxxyy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 45);

    auto g_xxxyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 46);

    auto g_xxxyy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 47);

    auto g_xxxyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 48);

    auto g_xxxyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 49);

    auto g_xxxyy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 50);

    auto g_xxxyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 51);

    auto g_xxxyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 52);

    auto g_xxxyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 53);

    auto g_xxxyy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 54);

    auto g_xxxyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 55);

    auto g_xxxyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 56);

    auto g_xxxyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 57);

    auto g_xxxyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 58);

    auto g_xxxyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 59);

    auto g_xxxzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 75);

    auto g_xxxzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 76);

    auto g_xxxzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 77);

    auto g_xxxzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 78);

    auto g_xxxzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 79);

    auto g_xxxzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 80);

    auto g_xxxzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 81);

    auto g_xxxzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 82);

    auto g_xxxzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 83);

    auto g_xxxzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 84);

    auto g_xxxzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 85);

    auto g_xxxzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 86);

    auto g_xxxzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 87);

    auto g_xxxzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 88);

    auto g_xxxzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 89);

    auto g_xxyyy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 90);

    auto g_xxyyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 91);

    auto g_xxyyy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 92);

    auto g_xxyyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 93);

    auto g_xxyyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 94);

    auto g_xxyyy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 95);

    auto g_xxyyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 96);

    auto g_xxyyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 97);

    auto g_xxyyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 98);

    auto g_xxyyy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 99);

    auto g_xxyyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 100);

    auto g_xxyyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 101);

    auto g_xxyyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 102);

    auto g_xxyyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 103);

    auto g_xxyyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 104);

    auto g_xxyyz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 106);

    auto g_xxyyz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 108);

    auto g_xxyyz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 111);

    auto g_xxyzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 120);

    auto g_xxyzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 122);

    auto g_xxyzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 125);

    auto g_xxyzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 129);

    auto g_xxzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 135);

    auto g_xxzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 136);

    auto g_xxzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 137);

    auto g_xxzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 138);

    auto g_xxzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 139);

    auto g_xxzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 140);

    auto g_xxzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 141);

    auto g_xxzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 142);

    auto g_xxzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 143);

    auto g_xxzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 144);

    auto g_xxzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 145);

    auto g_xxzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 146);

    auto g_xxzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 147);

    auto g_xxzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 148);

    auto g_xxzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 149);

    auto g_xyyyy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 150);

    auto g_xyyyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 151);

    auto g_xyyyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 153);

    auto g_xyyyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 154);

    auto g_xyyyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 156);

    auto g_xyyyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 157);

    auto g_xyyyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 158);

    auto g_xyyyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 160);

    auto g_xyyyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 161);

    auto g_xyyyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 162);

    auto g_xyyyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 163);

    auto g_xyyyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 164);

    auto g_xyyzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 184);

    auto g_xyyzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 187);

    auto g_xyyzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 188);

    auto g_xyyzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 190);

    auto g_xyyzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 191);

    auto g_xyyzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 192);

    auto g_xyyzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 193);

    auto g_xyyzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 194);

    auto g_xzzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 210);

    auto g_xzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 212);

    auto g_xzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 214);

    auto g_xzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 215);

    auto g_xzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 217);

    auto g_xzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 218);

    auto g_xzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 219);

    auto g_xzzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 220);

    auto g_xzzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 221);

    auto g_xzzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 222);

    auto g_xzzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 223);

    auto g_xzzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 224);

    auto g_yyyyy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 225);

    auto g_yyyyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 226);

    auto g_yyyyy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 227);

    auto g_yyyyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 228);

    auto g_yyyyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 229);

    auto g_yyyyy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 230);

    auto g_yyyyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 231);

    auto g_yyyyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 232);

    auto g_yyyyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 233);

    auto g_yyyyy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 234);

    auto g_yyyyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 235);

    auto g_yyyyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 236);

    auto g_yyyyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 237);

    auto g_yyyyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 238);

    auto g_yyyyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 239);

    auto g_yyyyz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 241);

    auto g_yyyyz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 242);

    auto g_yyyyz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 243);

    auto g_yyyyz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 244);

    auto g_yyyyz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 245);

    auto g_yyyyz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 246);

    auto g_yyyyz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 247);

    auto g_yyyyz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 248);

    auto g_yyyyz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 249);

    auto g_yyyyz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 250);

    auto g_yyyyz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 251);

    auto g_yyyyz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 252);

    auto g_yyyyz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 253);

    auto g_yyyyz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 254);

    auto g_yyyzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 255);

    auto g_yyyzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 256);

    auto g_yyyzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 257);

    auto g_yyyzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 258);

    auto g_yyyzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 259);

    auto g_yyyzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 260);

    auto g_yyyzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 261);

    auto g_yyyzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 262);

    auto g_yyyzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 263);

    auto g_yyyzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 264);

    auto g_yyyzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 265);

    auto g_yyyzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 266);

    auto g_yyyzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 267);

    auto g_yyyzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 268);

    auto g_yyyzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 269);

    auto g_yyzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 270);

    auto g_yyzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 271);

    auto g_yyzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 272);

    auto g_yyzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 273);

    auto g_yyzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 274);

    auto g_yyzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 275);

    auto g_yyzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 276);

    auto g_yyzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 277);

    auto g_yyzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 278);

    auto g_yyzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 279);

    auto g_yyzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 280);

    auto g_yyzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 281);

    auto g_yyzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 282);

    auto g_yyzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 283);

    auto g_yyzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 284);

    auto g_yzzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 285);

    auto g_yzzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 286);

    auto g_yzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 287);

    auto g_yzzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 288);

    auto g_yzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 289);

    auto g_yzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 290);

    auto g_yzzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 291);

    auto g_yzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 292);

    auto g_yzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 293);

    auto g_yzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 294);

    auto g_yzzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 295);

    auto g_yzzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 296);

    auto g_yzzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 297);

    auto g_yzzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 298);

    auto g_yzzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 299);

    auto g_zzzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 300);

    auto g_zzzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 301);

    auto g_zzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 302);

    auto g_zzzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 303);

    auto g_zzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 304);

    auto g_zzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 305);

    auto g_zzzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 306);

    auto g_zzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 307);

    auto g_zzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 308);

    auto g_zzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 309);

    auto g_zzzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 310);

    auto g_zzzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 311);

    auto g_zzzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 312);

    auto g_zzzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 313);

    auto g_zzzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 314);

    // Set up 0-15 components of targeted buffer : IG

    auto g_xxxxxx_xxxx_0 = pbuffer.data(idx_eri_0_ig);

    auto g_xxxxxx_xxxy_0 = pbuffer.data(idx_eri_0_ig + 1);

    auto g_xxxxxx_xxxz_0 = pbuffer.data(idx_eri_0_ig + 2);

    auto g_xxxxxx_xxyy_0 = pbuffer.data(idx_eri_0_ig + 3);

    auto g_xxxxxx_xxyz_0 = pbuffer.data(idx_eri_0_ig + 4);

    auto g_xxxxxx_xxzz_0 = pbuffer.data(idx_eri_0_ig + 5);

    auto g_xxxxxx_xyyy_0 = pbuffer.data(idx_eri_0_ig + 6);

    auto g_xxxxxx_xyyz_0 = pbuffer.data(idx_eri_0_ig + 7);

    auto g_xxxxxx_xyzz_0 = pbuffer.data(idx_eri_0_ig + 8);

    auto g_xxxxxx_xzzz_0 = pbuffer.data(idx_eri_0_ig + 9);

    auto g_xxxxxx_yyyy_0 = pbuffer.data(idx_eri_0_ig + 10);

    auto g_xxxxxx_yyyz_0 = pbuffer.data(idx_eri_0_ig + 11);

    auto g_xxxxxx_yyzz_0 = pbuffer.data(idx_eri_0_ig + 12);

    auto g_xxxxxx_yzzz_0 = pbuffer.data(idx_eri_0_ig + 13);

    auto g_xxxxxx_zzzz_0 = pbuffer.data(idx_eri_0_ig + 14);

    #pragma omp simd aligned(g_xxxx_xxxx_0, g_xxxx_xxxx_1, g_xxxx_xxxy_0, g_xxxx_xxxy_1, g_xxxx_xxxz_0, g_xxxx_xxxz_1, g_xxxx_xxyy_0, g_xxxx_xxyy_1, g_xxxx_xxyz_0, g_xxxx_xxyz_1, g_xxxx_xxzz_0, g_xxxx_xxzz_1, g_xxxx_xyyy_0, g_xxxx_xyyy_1, g_xxxx_xyyz_0, g_xxxx_xyyz_1, g_xxxx_xyzz_0, g_xxxx_xyzz_1, g_xxxx_xzzz_0, g_xxxx_xzzz_1, g_xxxx_yyyy_0, g_xxxx_yyyy_1, g_xxxx_yyyz_0, g_xxxx_yyyz_1, g_xxxx_yyzz_0, g_xxxx_yyzz_1, g_xxxx_yzzz_0, g_xxxx_yzzz_1, g_xxxx_zzzz_0, g_xxxx_zzzz_1, g_xxxxx_xxx_1, g_xxxxx_xxxx_1, g_xxxxx_xxxy_1, g_xxxxx_xxxz_1, g_xxxxx_xxy_1, g_xxxxx_xxyy_1, g_xxxxx_xxyz_1, g_xxxxx_xxz_1, g_xxxxx_xxzz_1, g_xxxxx_xyy_1, g_xxxxx_xyyy_1, g_xxxxx_xyyz_1, g_xxxxx_xyz_1, g_xxxxx_xyzz_1, g_xxxxx_xzz_1, g_xxxxx_xzzz_1, g_xxxxx_yyy_1, g_xxxxx_yyyy_1, g_xxxxx_yyyz_1, g_xxxxx_yyz_1, g_xxxxx_yyzz_1, g_xxxxx_yzz_1, g_xxxxx_yzzz_1, g_xxxxx_zzz_1, g_xxxxx_zzzz_1, g_xxxxxx_xxxx_0, g_xxxxxx_xxxy_0, g_xxxxxx_xxxz_0, g_xxxxxx_xxyy_0, g_xxxxxx_xxyz_0, g_xxxxxx_xxzz_0, g_xxxxxx_xyyy_0, g_xxxxxx_xyyz_0, g_xxxxxx_xyzz_0, g_xxxxxx_xzzz_0, g_xxxxxx_yyyy_0, g_xxxxxx_yyyz_0, g_xxxxxx_yyzz_0, g_xxxxxx_yzzz_0, g_xxxxxx_zzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxx_xxxx_0[i] = 5.0 * g_xxxx_xxxx_0[i] * fbe_0 - 5.0 * g_xxxx_xxxx_1[i] * fz_be_0 + 4.0 * g_xxxxx_xxx_1[i] * fe_0 + g_xxxxx_xxxx_1[i] * pa_x[i];

        g_xxxxxx_xxxy_0[i] = 5.0 * g_xxxx_xxxy_0[i] * fbe_0 - 5.0 * g_xxxx_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxy_1[i] * fe_0 + g_xxxxx_xxxy_1[i] * pa_x[i];

        g_xxxxxx_xxxz_0[i] = 5.0 * g_xxxx_xxxz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxz_1[i] * fe_0 + g_xxxxx_xxxz_1[i] * pa_x[i];

        g_xxxxxx_xxyy_0[i] = 5.0 * g_xxxx_xxyy_0[i] * fbe_0 - 5.0 * g_xxxx_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyy_1[i] * fe_0 + g_xxxxx_xxyy_1[i] * pa_x[i];

        g_xxxxxx_xxyz_0[i] = 5.0 * g_xxxx_xxyz_0[i] * fbe_0 - 5.0 * g_xxxx_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyz_1[i] * fe_0 + g_xxxxx_xxyz_1[i] * pa_x[i];

        g_xxxxxx_xxzz_0[i] = 5.0 * g_xxxx_xxzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xzz_1[i] * fe_0 + g_xxxxx_xxzz_1[i] * pa_x[i];

        g_xxxxxx_xyyy_0[i] = 5.0 * g_xxxx_xyyy_0[i] * fbe_0 - 5.0 * g_xxxx_xyyy_1[i] * fz_be_0 + g_xxxxx_yyy_1[i] * fe_0 + g_xxxxx_xyyy_1[i] * pa_x[i];

        g_xxxxxx_xyyz_0[i] = 5.0 * g_xxxx_xyyz_0[i] * fbe_0 - 5.0 * g_xxxx_xyyz_1[i] * fz_be_0 + g_xxxxx_yyz_1[i] * fe_0 + g_xxxxx_xyyz_1[i] * pa_x[i];

        g_xxxxxx_xyzz_0[i] = 5.0 * g_xxxx_xyzz_0[i] * fbe_0 - 5.0 * g_xxxx_xyzz_1[i] * fz_be_0 + g_xxxxx_yzz_1[i] * fe_0 + g_xxxxx_xyzz_1[i] * pa_x[i];

        g_xxxxxx_xzzz_0[i] = 5.0 * g_xxxx_xzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xzzz_1[i] * fz_be_0 + g_xxxxx_zzz_1[i] * fe_0 + g_xxxxx_xzzz_1[i] * pa_x[i];

        g_xxxxxx_yyyy_0[i] = 5.0 * g_xxxx_yyyy_0[i] * fbe_0 - 5.0 * g_xxxx_yyyy_1[i] * fz_be_0 + g_xxxxx_yyyy_1[i] * pa_x[i];

        g_xxxxxx_yyyz_0[i] = 5.0 * g_xxxx_yyyz_0[i] * fbe_0 - 5.0 * g_xxxx_yyyz_1[i] * fz_be_0 + g_xxxxx_yyyz_1[i] * pa_x[i];

        g_xxxxxx_yyzz_0[i] = 5.0 * g_xxxx_yyzz_0[i] * fbe_0 - 5.0 * g_xxxx_yyzz_1[i] * fz_be_0 + g_xxxxx_yyzz_1[i] * pa_x[i];

        g_xxxxxx_yzzz_0[i] = 5.0 * g_xxxx_yzzz_0[i] * fbe_0 - 5.0 * g_xxxx_yzzz_1[i] * fz_be_0 + g_xxxxx_yzzz_1[i] * pa_x[i];

        g_xxxxxx_zzzz_0[i] = 5.0 * g_xxxx_zzzz_0[i] * fbe_0 - 5.0 * g_xxxx_zzzz_1[i] * fz_be_0 + g_xxxxx_zzzz_1[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : IG

    auto g_xxxxxy_xxxx_0 = pbuffer.data(idx_eri_0_ig + 15);

    auto g_xxxxxy_xxxy_0 = pbuffer.data(idx_eri_0_ig + 16);

    auto g_xxxxxy_xxxz_0 = pbuffer.data(idx_eri_0_ig + 17);

    auto g_xxxxxy_xxyy_0 = pbuffer.data(idx_eri_0_ig + 18);

    auto g_xxxxxy_xxyz_0 = pbuffer.data(idx_eri_0_ig + 19);

    auto g_xxxxxy_xxzz_0 = pbuffer.data(idx_eri_0_ig + 20);

    auto g_xxxxxy_xyyy_0 = pbuffer.data(idx_eri_0_ig + 21);

    auto g_xxxxxy_xyyz_0 = pbuffer.data(idx_eri_0_ig + 22);

    auto g_xxxxxy_xyzz_0 = pbuffer.data(idx_eri_0_ig + 23);

    auto g_xxxxxy_xzzz_0 = pbuffer.data(idx_eri_0_ig + 24);

    auto g_xxxxxy_yyyy_0 = pbuffer.data(idx_eri_0_ig + 25);

    auto g_xxxxxy_yyyz_0 = pbuffer.data(idx_eri_0_ig + 26);

    auto g_xxxxxy_yyzz_0 = pbuffer.data(idx_eri_0_ig + 27);

    auto g_xxxxxy_yzzz_0 = pbuffer.data(idx_eri_0_ig + 28);

    auto g_xxxxxy_zzzz_0 = pbuffer.data(idx_eri_0_ig + 29);

    #pragma omp simd aligned(g_xxxxx_xxx_1, g_xxxxx_xxxx_1, g_xxxxx_xxxy_1, g_xxxxx_xxxz_1, g_xxxxx_xxy_1, g_xxxxx_xxyy_1, g_xxxxx_xxyz_1, g_xxxxx_xxz_1, g_xxxxx_xxzz_1, g_xxxxx_xyy_1, g_xxxxx_xyyy_1, g_xxxxx_xyyz_1, g_xxxxx_xyz_1, g_xxxxx_xyzz_1, g_xxxxx_xzz_1, g_xxxxx_xzzz_1, g_xxxxx_yyy_1, g_xxxxx_yyyy_1, g_xxxxx_yyyz_1, g_xxxxx_yyz_1, g_xxxxx_yyzz_1, g_xxxxx_yzz_1, g_xxxxx_yzzz_1, g_xxxxx_zzz_1, g_xxxxx_zzzz_1, g_xxxxxy_xxxx_0, g_xxxxxy_xxxy_0, g_xxxxxy_xxxz_0, g_xxxxxy_xxyy_0, g_xxxxxy_xxyz_0, g_xxxxxy_xxzz_0, g_xxxxxy_xyyy_0, g_xxxxxy_xyyz_0, g_xxxxxy_xyzz_0, g_xxxxxy_xzzz_0, g_xxxxxy_yyyy_0, g_xxxxxy_yyyz_0, g_xxxxxy_yyzz_0, g_xxxxxy_yzzz_0, g_xxxxxy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxy_xxxx_0[i] = g_xxxxx_xxxx_1[i] * pa_y[i];

        g_xxxxxy_xxxy_0[i] = g_xxxxx_xxx_1[i] * fe_0 + g_xxxxx_xxxy_1[i] * pa_y[i];

        g_xxxxxy_xxxz_0[i] = g_xxxxx_xxxz_1[i] * pa_y[i];

        g_xxxxxy_xxyy_0[i] = 2.0 * g_xxxxx_xxy_1[i] * fe_0 + g_xxxxx_xxyy_1[i] * pa_y[i];

        g_xxxxxy_xxyz_0[i] = g_xxxxx_xxz_1[i] * fe_0 + g_xxxxx_xxyz_1[i] * pa_y[i];

        g_xxxxxy_xxzz_0[i] = g_xxxxx_xxzz_1[i] * pa_y[i];

        g_xxxxxy_xyyy_0[i] = 3.0 * g_xxxxx_xyy_1[i] * fe_0 + g_xxxxx_xyyy_1[i] * pa_y[i];

        g_xxxxxy_xyyz_0[i] = 2.0 * g_xxxxx_xyz_1[i] * fe_0 + g_xxxxx_xyyz_1[i] * pa_y[i];

        g_xxxxxy_xyzz_0[i] = g_xxxxx_xzz_1[i] * fe_0 + g_xxxxx_xyzz_1[i] * pa_y[i];

        g_xxxxxy_xzzz_0[i] = g_xxxxx_xzzz_1[i] * pa_y[i];

        g_xxxxxy_yyyy_0[i] = 4.0 * g_xxxxx_yyy_1[i] * fe_0 + g_xxxxx_yyyy_1[i] * pa_y[i];

        g_xxxxxy_yyyz_0[i] = 3.0 * g_xxxxx_yyz_1[i] * fe_0 + g_xxxxx_yyyz_1[i] * pa_y[i];

        g_xxxxxy_yyzz_0[i] = 2.0 * g_xxxxx_yzz_1[i] * fe_0 + g_xxxxx_yyzz_1[i] * pa_y[i];

        g_xxxxxy_yzzz_0[i] = g_xxxxx_zzz_1[i] * fe_0 + g_xxxxx_yzzz_1[i] * pa_y[i];

        g_xxxxxy_zzzz_0[i] = g_xxxxx_zzzz_1[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : IG

    auto g_xxxxxz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 30);

    auto g_xxxxxz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 31);

    auto g_xxxxxz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 32);

    auto g_xxxxxz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 33);

    auto g_xxxxxz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 34);

    auto g_xxxxxz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 35);

    auto g_xxxxxz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 36);

    auto g_xxxxxz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 37);

    auto g_xxxxxz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 38);

    auto g_xxxxxz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 39);

    auto g_xxxxxz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 40);

    auto g_xxxxxz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 41);

    auto g_xxxxxz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 42);

    auto g_xxxxxz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 43);

    auto g_xxxxxz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 44);

    #pragma omp simd aligned(g_xxxxx_xxx_1, g_xxxxx_xxxx_1, g_xxxxx_xxxy_1, g_xxxxx_xxxz_1, g_xxxxx_xxy_1, g_xxxxx_xxyy_1, g_xxxxx_xxyz_1, g_xxxxx_xxz_1, g_xxxxx_xxzz_1, g_xxxxx_xyy_1, g_xxxxx_xyyy_1, g_xxxxx_xyyz_1, g_xxxxx_xyz_1, g_xxxxx_xyzz_1, g_xxxxx_xzz_1, g_xxxxx_xzzz_1, g_xxxxx_yyy_1, g_xxxxx_yyyy_1, g_xxxxx_yyyz_1, g_xxxxx_yyz_1, g_xxxxx_yyzz_1, g_xxxxx_yzz_1, g_xxxxx_yzzz_1, g_xxxxx_zzz_1, g_xxxxx_zzzz_1, g_xxxxxz_xxxx_0, g_xxxxxz_xxxy_0, g_xxxxxz_xxxz_0, g_xxxxxz_xxyy_0, g_xxxxxz_xxyz_0, g_xxxxxz_xxzz_0, g_xxxxxz_xyyy_0, g_xxxxxz_xyyz_0, g_xxxxxz_xyzz_0, g_xxxxxz_xzzz_0, g_xxxxxz_yyyy_0, g_xxxxxz_yyyz_0, g_xxxxxz_yyzz_0, g_xxxxxz_yzzz_0, g_xxxxxz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxz_xxxx_0[i] = g_xxxxx_xxxx_1[i] * pa_z[i];

        g_xxxxxz_xxxy_0[i] = g_xxxxx_xxxy_1[i] * pa_z[i];

        g_xxxxxz_xxxz_0[i] = g_xxxxx_xxx_1[i] * fe_0 + g_xxxxx_xxxz_1[i] * pa_z[i];

        g_xxxxxz_xxyy_0[i] = g_xxxxx_xxyy_1[i] * pa_z[i];

        g_xxxxxz_xxyz_0[i] = g_xxxxx_xxy_1[i] * fe_0 + g_xxxxx_xxyz_1[i] * pa_z[i];

        g_xxxxxz_xxzz_0[i] = 2.0 * g_xxxxx_xxz_1[i] * fe_0 + g_xxxxx_xxzz_1[i] * pa_z[i];

        g_xxxxxz_xyyy_0[i] = g_xxxxx_xyyy_1[i] * pa_z[i];

        g_xxxxxz_xyyz_0[i] = g_xxxxx_xyy_1[i] * fe_0 + g_xxxxx_xyyz_1[i] * pa_z[i];

        g_xxxxxz_xyzz_0[i] = 2.0 * g_xxxxx_xyz_1[i] * fe_0 + g_xxxxx_xyzz_1[i] * pa_z[i];

        g_xxxxxz_xzzz_0[i] = 3.0 * g_xxxxx_xzz_1[i] * fe_0 + g_xxxxx_xzzz_1[i] * pa_z[i];

        g_xxxxxz_yyyy_0[i] = g_xxxxx_yyyy_1[i] * pa_z[i];

        g_xxxxxz_yyyz_0[i] = g_xxxxx_yyy_1[i] * fe_0 + g_xxxxx_yyyz_1[i] * pa_z[i];

        g_xxxxxz_yyzz_0[i] = 2.0 * g_xxxxx_yyz_1[i] * fe_0 + g_xxxxx_yyzz_1[i] * pa_z[i];

        g_xxxxxz_yzzz_0[i] = 3.0 * g_xxxxx_yzz_1[i] * fe_0 + g_xxxxx_yzzz_1[i] * pa_z[i];

        g_xxxxxz_zzzz_0[i] = 4.0 * g_xxxxx_zzz_1[i] * fe_0 + g_xxxxx_zzzz_1[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : IG

    auto g_xxxxyy_xxxx_0 = pbuffer.data(idx_eri_0_ig + 45);

    auto g_xxxxyy_xxxy_0 = pbuffer.data(idx_eri_0_ig + 46);

    auto g_xxxxyy_xxxz_0 = pbuffer.data(idx_eri_0_ig + 47);

    auto g_xxxxyy_xxyy_0 = pbuffer.data(idx_eri_0_ig + 48);

    auto g_xxxxyy_xxyz_0 = pbuffer.data(idx_eri_0_ig + 49);

    auto g_xxxxyy_xxzz_0 = pbuffer.data(idx_eri_0_ig + 50);

    auto g_xxxxyy_xyyy_0 = pbuffer.data(idx_eri_0_ig + 51);

    auto g_xxxxyy_xyyz_0 = pbuffer.data(idx_eri_0_ig + 52);

    auto g_xxxxyy_xyzz_0 = pbuffer.data(idx_eri_0_ig + 53);

    auto g_xxxxyy_xzzz_0 = pbuffer.data(idx_eri_0_ig + 54);

    auto g_xxxxyy_yyyy_0 = pbuffer.data(idx_eri_0_ig + 55);

    auto g_xxxxyy_yyyz_0 = pbuffer.data(idx_eri_0_ig + 56);

    auto g_xxxxyy_yyzz_0 = pbuffer.data(idx_eri_0_ig + 57);

    auto g_xxxxyy_yzzz_0 = pbuffer.data(idx_eri_0_ig + 58);

    auto g_xxxxyy_zzzz_0 = pbuffer.data(idx_eri_0_ig + 59);

    #pragma omp simd aligned(g_xxxx_xxxx_0, g_xxxx_xxxx_1, g_xxxx_xxxz_0, g_xxxx_xxxz_1, g_xxxx_xxzz_0, g_xxxx_xxzz_1, g_xxxx_xzzz_0, g_xxxx_xzzz_1, g_xxxxy_xxxx_1, g_xxxxy_xxxz_1, g_xxxxy_xxzz_1, g_xxxxy_xzzz_1, g_xxxxyy_xxxx_0, g_xxxxyy_xxxy_0, g_xxxxyy_xxxz_0, g_xxxxyy_xxyy_0, g_xxxxyy_xxyz_0, g_xxxxyy_xxzz_0, g_xxxxyy_xyyy_0, g_xxxxyy_xyyz_0, g_xxxxyy_xyzz_0, g_xxxxyy_xzzz_0, g_xxxxyy_yyyy_0, g_xxxxyy_yyyz_0, g_xxxxyy_yyzz_0, g_xxxxyy_yzzz_0, g_xxxxyy_zzzz_0, g_xxxyy_xxxy_1, g_xxxyy_xxy_1, g_xxxyy_xxyy_1, g_xxxyy_xxyz_1, g_xxxyy_xyy_1, g_xxxyy_xyyy_1, g_xxxyy_xyyz_1, g_xxxyy_xyz_1, g_xxxyy_xyzz_1, g_xxxyy_yyy_1, g_xxxyy_yyyy_1, g_xxxyy_yyyz_1, g_xxxyy_yyz_1, g_xxxyy_yyzz_1, g_xxxyy_yzz_1, g_xxxyy_yzzz_1, g_xxxyy_zzzz_1, g_xxyy_xxxy_0, g_xxyy_xxxy_1, g_xxyy_xxyy_0, g_xxyy_xxyy_1, g_xxyy_xxyz_0, g_xxyy_xxyz_1, g_xxyy_xyyy_0, g_xxyy_xyyy_1, g_xxyy_xyyz_0, g_xxyy_xyyz_1, g_xxyy_xyzz_0, g_xxyy_xyzz_1, g_xxyy_yyyy_0, g_xxyy_yyyy_1, g_xxyy_yyyz_0, g_xxyy_yyyz_1, g_xxyy_yyzz_0, g_xxyy_yyzz_1, g_xxyy_yzzz_0, g_xxyy_yzzz_1, g_xxyy_zzzz_0, g_xxyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyy_xxxx_0[i] = g_xxxx_xxxx_0[i] * fbe_0 - g_xxxx_xxxx_1[i] * fz_be_0 + g_xxxxy_xxxx_1[i] * pa_y[i];

        g_xxxxyy_xxxy_0[i] = 3.0 * g_xxyy_xxxy_0[i] * fbe_0 - 3.0 * g_xxyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxyy_xxy_1[i] * fe_0 + g_xxxyy_xxxy_1[i] * pa_x[i];

        g_xxxxyy_xxxz_0[i] = g_xxxx_xxxz_0[i] * fbe_0 - g_xxxx_xxxz_1[i] * fz_be_0 + g_xxxxy_xxxz_1[i] * pa_y[i];

        g_xxxxyy_xxyy_0[i] = 3.0 * g_xxyy_xxyy_0[i] * fbe_0 - 3.0 * g_xxyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyy_1[i] * fe_0 + g_xxxyy_xxyy_1[i] * pa_x[i];

        g_xxxxyy_xxyz_0[i] = 3.0 * g_xxyy_xxyz_0[i] * fbe_0 - 3.0 * g_xxyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyz_1[i] * fe_0 + g_xxxyy_xxyz_1[i] * pa_x[i];

        g_xxxxyy_xxzz_0[i] = g_xxxx_xxzz_0[i] * fbe_0 - g_xxxx_xxzz_1[i] * fz_be_0 + g_xxxxy_xxzz_1[i] * pa_y[i];

        g_xxxxyy_xyyy_0[i] = 3.0 * g_xxyy_xyyy_0[i] * fbe_0 - 3.0 * g_xxyy_xyyy_1[i] * fz_be_0 + g_xxxyy_yyy_1[i] * fe_0 + g_xxxyy_xyyy_1[i] * pa_x[i];

        g_xxxxyy_xyyz_0[i] = 3.0 * g_xxyy_xyyz_0[i] * fbe_0 - 3.0 * g_xxyy_xyyz_1[i] * fz_be_0 + g_xxxyy_yyz_1[i] * fe_0 + g_xxxyy_xyyz_1[i] * pa_x[i];

        g_xxxxyy_xyzz_0[i] = 3.0 * g_xxyy_xyzz_0[i] * fbe_0 - 3.0 * g_xxyy_xyzz_1[i] * fz_be_0 + g_xxxyy_yzz_1[i] * fe_0 + g_xxxyy_xyzz_1[i] * pa_x[i];

        g_xxxxyy_xzzz_0[i] = g_xxxx_xzzz_0[i] * fbe_0 - g_xxxx_xzzz_1[i] * fz_be_0 + g_xxxxy_xzzz_1[i] * pa_y[i];

        g_xxxxyy_yyyy_0[i] = 3.0 * g_xxyy_yyyy_0[i] * fbe_0 - 3.0 * g_xxyy_yyyy_1[i] * fz_be_0 + g_xxxyy_yyyy_1[i] * pa_x[i];

        g_xxxxyy_yyyz_0[i] = 3.0 * g_xxyy_yyyz_0[i] * fbe_0 - 3.0 * g_xxyy_yyyz_1[i] * fz_be_0 + g_xxxyy_yyyz_1[i] * pa_x[i];

        g_xxxxyy_yyzz_0[i] = 3.0 * g_xxyy_yyzz_0[i] * fbe_0 - 3.0 * g_xxyy_yyzz_1[i] * fz_be_0 + g_xxxyy_yyzz_1[i] * pa_x[i];

        g_xxxxyy_yzzz_0[i] = 3.0 * g_xxyy_yzzz_0[i] * fbe_0 - 3.0 * g_xxyy_yzzz_1[i] * fz_be_0 + g_xxxyy_yzzz_1[i] * pa_x[i];

        g_xxxxyy_zzzz_0[i] = 3.0 * g_xxyy_zzzz_0[i] * fbe_0 - 3.0 * g_xxyy_zzzz_1[i] * fz_be_0 + g_xxxyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : IG

    auto g_xxxxyz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 60);

    auto g_xxxxyz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 61);

    auto g_xxxxyz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 62);

    auto g_xxxxyz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 63);

    auto g_xxxxyz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 64);

    auto g_xxxxyz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 65);

    auto g_xxxxyz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 66);

    auto g_xxxxyz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 67);

    auto g_xxxxyz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 68);

    auto g_xxxxyz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 69);

    auto g_xxxxyz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 70);

    auto g_xxxxyz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 71);

    auto g_xxxxyz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 72);

    auto g_xxxxyz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 73);

    auto g_xxxxyz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 74);

    #pragma omp simd aligned(g_xxxxy_xxxy_1, g_xxxxy_xxyy_1, g_xxxxy_xyyy_1, g_xxxxy_yyyy_1, g_xxxxyz_xxxx_0, g_xxxxyz_xxxy_0, g_xxxxyz_xxxz_0, g_xxxxyz_xxyy_0, g_xxxxyz_xxyz_0, g_xxxxyz_xxzz_0, g_xxxxyz_xyyy_0, g_xxxxyz_xyyz_0, g_xxxxyz_xyzz_0, g_xxxxyz_xzzz_0, g_xxxxyz_yyyy_0, g_xxxxyz_yyyz_0, g_xxxxyz_yyzz_0, g_xxxxyz_yzzz_0, g_xxxxyz_zzzz_0, g_xxxxz_xxxx_1, g_xxxxz_xxxz_1, g_xxxxz_xxyz_1, g_xxxxz_xxz_1, g_xxxxz_xxzz_1, g_xxxxz_xyyz_1, g_xxxxz_xyz_1, g_xxxxz_xyzz_1, g_xxxxz_xzz_1, g_xxxxz_xzzz_1, g_xxxxz_yyyz_1, g_xxxxz_yyz_1, g_xxxxz_yyzz_1, g_xxxxz_yzz_1, g_xxxxz_yzzz_1, g_xxxxz_zzz_1, g_xxxxz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyz_xxxx_0[i] = g_xxxxz_xxxx_1[i] * pa_y[i];

        g_xxxxyz_xxxy_0[i] = g_xxxxy_xxxy_1[i] * pa_z[i];

        g_xxxxyz_xxxz_0[i] = g_xxxxz_xxxz_1[i] * pa_y[i];

        g_xxxxyz_xxyy_0[i] = g_xxxxy_xxyy_1[i] * pa_z[i];

        g_xxxxyz_xxyz_0[i] = g_xxxxz_xxz_1[i] * fe_0 + g_xxxxz_xxyz_1[i] * pa_y[i];

        g_xxxxyz_xxzz_0[i] = g_xxxxz_xxzz_1[i] * pa_y[i];

        g_xxxxyz_xyyy_0[i] = g_xxxxy_xyyy_1[i] * pa_z[i];

        g_xxxxyz_xyyz_0[i] = 2.0 * g_xxxxz_xyz_1[i] * fe_0 + g_xxxxz_xyyz_1[i] * pa_y[i];

        g_xxxxyz_xyzz_0[i] = g_xxxxz_xzz_1[i] * fe_0 + g_xxxxz_xyzz_1[i] * pa_y[i];

        g_xxxxyz_xzzz_0[i] = g_xxxxz_xzzz_1[i] * pa_y[i];

        g_xxxxyz_yyyy_0[i] = g_xxxxy_yyyy_1[i] * pa_z[i];

        g_xxxxyz_yyyz_0[i] = 3.0 * g_xxxxz_yyz_1[i] * fe_0 + g_xxxxz_yyyz_1[i] * pa_y[i];

        g_xxxxyz_yyzz_0[i] = 2.0 * g_xxxxz_yzz_1[i] * fe_0 + g_xxxxz_yyzz_1[i] * pa_y[i];

        g_xxxxyz_yzzz_0[i] = g_xxxxz_zzz_1[i] * fe_0 + g_xxxxz_yzzz_1[i] * pa_y[i];

        g_xxxxyz_zzzz_0[i] = g_xxxxz_zzzz_1[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : IG

    auto g_xxxxzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 75);

    auto g_xxxxzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 76);

    auto g_xxxxzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 77);

    auto g_xxxxzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 78);

    auto g_xxxxzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 79);

    auto g_xxxxzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 80);

    auto g_xxxxzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 81);

    auto g_xxxxzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 82);

    auto g_xxxxzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 83);

    auto g_xxxxzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 84);

    auto g_xxxxzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 85);

    auto g_xxxxzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 86);

    auto g_xxxxzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 87);

    auto g_xxxxzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 88);

    auto g_xxxxzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 89);

    #pragma omp simd aligned(g_xxxx_xxxx_0, g_xxxx_xxxx_1, g_xxxx_xxxy_0, g_xxxx_xxxy_1, g_xxxx_xxyy_0, g_xxxx_xxyy_1, g_xxxx_xyyy_0, g_xxxx_xyyy_1, g_xxxxz_xxxx_1, g_xxxxz_xxxy_1, g_xxxxz_xxyy_1, g_xxxxz_xyyy_1, g_xxxxzz_xxxx_0, g_xxxxzz_xxxy_0, g_xxxxzz_xxxz_0, g_xxxxzz_xxyy_0, g_xxxxzz_xxyz_0, g_xxxxzz_xxzz_0, g_xxxxzz_xyyy_0, g_xxxxzz_xyyz_0, g_xxxxzz_xyzz_0, g_xxxxzz_xzzz_0, g_xxxxzz_yyyy_0, g_xxxxzz_yyyz_0, g_xxxxzz_yyzz_0, g_xxxxzz_yzzz_0, g_xxxxzz_zzzz_0, g_xxxzz_xxxz_1, g_xxxzz_xxyz_1, g_xxxzz_xxz_1, g_xxxzz_xxzz_1, g_xxxzz_xyyz_1, g_xxxzz_xyz_1, g_xxxzz_xyzz_1, g_xxxzz_xzz_1, g_xxxzz_xzzz_1, g_xxxzz_yyyy_1, g_xxxzz_yyyz_1, g_xxxzz_yyz_1, g_xxxzz_yyzz_1, g_xxxzz_yzz_1, g_xxxzz_yzzz_1, g_xxxzz_zzz_1, g_xxxzz_zzzz_1, g_xxzz_xxxz_0, g_xxzz_xxxz_1, g_xxzz_xxyz_0, g_xxzz_xxyz_1, g_xxzz_xxzz_0, g_xxzz_xxzz_1, g_xxzz_xyyz_0, g_xxzz_xyyz_1, g_xxzz_xyzz_0, g_xxzz_xyzz_1, g_xxzz_xzzz_0, g_xxzz_xzzz_1, g_xxzz_yyyy_0, g_xxzz_yyyy_1, g_xxzz_yyyz_0, g_xxzz_yyyz_1, g_xxzz_yyzz_0, g_xxzz_yyzz_1, g_xxzz_yzzz_0, g_xxzz_yzzz_1, g_xxzz_zzzz_0, g_xxzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzz_xxxx_0[i] = g_xxxx_xxxx_0[i] * fbe_0 - g_xxxx_xxxx_1[i] * fz_be_0 + g_xxxxz_xxxx_1[i] * pa_z[i];

        g_xxxxzz_xxxy_0[i] = g_xxxx_xxxy_0[i] * fbe_0 - g_xxxx_xxxy_1[i] * fz_be_0 + g_xxxxz_xxxy_1[i] * pa_z[i];

        g_xxxxzz_xxxz_0[i] = 3.0 * g_xxzz_xxxz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxzz_xxz_1[i] * fe_0 + g_xxxzz_xxxz_1[i] * pa_x[i];

        g_xxxxzz_xxyy_0[i] = g_xxxx_xxyy_0[i] * fbe_0 - g_xxxx_xxyy_1[i] * fz_be_0 + g_xxxxz_xxyy_1[i] * pa_z[i];

        g_xxxxzz_xxyz_0[i] = 3.0 * g_xxzz_xxyz_0[i] * fbe_0 - 3.0 * g_xxzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xyz_1[i] * fe_0 + g_xxxzz_xxyz_1[i] * pa_x[i];

        g_xxxxzz_xxzz_0[i] = 3.0 * g_xxzz_xxzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xzz_1[i] * fe_0 + g_xxxzz_xxzz_1[i] * pa_x[i];

        g_xxxxzz_xyyy_0[i] = g_xxxx_xyyy_0[i] * fbe_0 - g_xxxx_xyyy_1[i] * fz_be_0 + g_xxxxz_xyyy_1[i] * pa_z[i];

        g_xxxxzz_xyyz_0[i] = 3.0 * g_xxzz_xyyz_0[i] * fbe_0 - 3.0 * g_xxzz_xyyz_1[i] * fz_be_0 + g_xxxzz_yyz_1[i] * fe_0 + g_xxxzz_xyyz_1[i] * pa_x[i];

        g_xxxxzz_xyzz_0[i] = 3.0 * g_xxzz_xyzz_0[i] * fbe_0 - 3.0 * g_xxzz_xyzz_1[i] * fz_be_0 + g_xxxzz_yzz_1[i] * fe_0 + g_xxxzz_xyzz_1[i] * pa_x[i];

        g_xxxxzz_xzzz_0[i] = 3.0 * g_xxzz_xzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xzzz_1[i] * fz_be_0 + g_xxxzz_zzz_1[i] * fe_0 + g_xxxzz_xzzz_1[i] * pa_x[i];

        g_xxxxzz_yyyy_0[i] = 3.0 * g_xxzz_yyyy_0[i] * fbe_0 - 3.0 * g_xxzz_yyyy_1[i] * fz_be_0 + g_xxxzz_yyyy_1[i] * pa_x[i];

        g_xxxxzz_yyyz_0[i] = 3.0 * g_xxzz_yyyz_0[i] * fbe_0 - 3.0 * g_xxzz_yyyz_1[i] * fz_be_0 + g_xxxzz_yyyz_1[i] * pa_x[i];

        g_xxxxzz_yyzz_0[i] = 3.0 * g_xxzz_yyzz_0[i] * fbe_0 - 3.0 * g_xxzz_yyzz_1[i] * fz_be_0 + g_xxxzz_yyzz_1[i] * pa_x[i];

        g_xxxxzz_yzzz_0[i] = 3.0 * g_xxzz_yzzz_0[i] * fbe_0 - 3.0 * g_xxzz_yzzz_1[i] * fz_be_0 + g_xxxzz_yzzz_1[i] * pa_x[i];

        g_xxxxzz_zzzz_0[i] = 3.0 * g_xxzz_zzzz_0[i] * fbe_0 - 3.0 * g_xxzz_zzzz_1[i] * fz_be_0 + g_xxxzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : IG

    auto g_xxxyyy_xxxx_0 = pbuffer.data(idx_eri_0_ig + 90);

    auto g_xxxyyy_xxxy_0 = pbuffer.data(idx_eri_0_ig + 91);

    auto g_xxxyyy_xxxz_0 = pbuffer.data(idx_eri_0_ig + 92);

    auto g_xxxyyy_xxyy_0 = pbuffer.data(idx_eri_0_ig + 93);

    auto g_xxxyyy_xxyz_0 = pbuffer.data(idx_eri_0_ig + 94);

    auto g_xxxyyy_xxzz_0 = pbuffer.data(idx_eri_0_ig + 95);

    auto g_xxxyyy_xyyy_0 = pbuffer.data(idx_eri_0_ig + 96);

    auto g_xxxyyy_xyyz_0 = pbuffer.data(idx_eri_0_ig + 97);

    auto g_xxxyyy_xyzz_0 = pbuffer.data(idx_eri_0_ig + 98);

    auto g_xxxyyy_xzzz_0 = pbuffer.data(idx_eri_0_ig + 99);

    auto g_xxxyyy_yyyy_0 = pbuffer.data(idx_eri_0_ig + 100);

    auto g_xxxyyy_yyyz_0 = pbuffer.data(idx_eri_0_ig + 101);

    auto g_xxxyyy_yyzz_0 = pbuffer.data(idx_eri_0_ig + 102);

    auto g_xxxyyy_yzzz_0 = pbuffer.data(idx_eri_0_ig + 103);

    auto g_xxxyyy_zzzz_0 = pbuffer.data(idx_eri_0_ig + 104);

    #pragma omp simd aligned(g_xxxy_xxxx_0, g_xxxy_xxxx_1, g_xxxy_xxxz_0, g_xxxy_xxxz_1, g_xxxy_xxzz_0, g_xxxy_xxzz_1, g_xxxy_xzzz_0, g_xxxy_xzzz_1, g_xxxyy_xxxx_1, g_xxxyy_xxxz_1, g_xxxyy_xxzz_1, g_xxxyy_xzzz_1, g_xxxyyy_xxxx_0, g_xxxyyy_xxxy_0, g_xxxyyy_xxxz_0, g_xxxyyy_xxyy_0, g_xxxyyy_xxyz_0, g_xxxyyy_xxzz_0, g_xxxyyy_xyyy_0, g_xxxyyy_xyyz_0, g_xxxyyy_xyzz_0, g_xxxyyy_xzzz_0, g_xxxyyy_yyyy_0, g_xxxyyy_yyyz_0, g_xxxyyy_yyzz_0, g_xxxyyy_yzzz_0, g_xxxyyy_zzzz_0, g_xxyyy_xxxy_1, g_xxyyy_xxy_1, g_xxyyy_xxyy_1, g_xxyyy_xxyz_1, g_xxyyy_xyy_1, g_xxyyy_xyyy_1, g_xxyyy_xyyz_1, g_xxyyy_xyz_1, g_xxyyy_xyzz_1, g_xxyyy_yyy_1, g_xxyyy_yyyy_1, g_xxyyy_yyyz_1, g_xxyyy_yyz_1, g_xxyyy_yyzz_1, g_xxyyy_yzz_1, g_xxyyy_yzzz_1, g_xxyyy_zzzz_1, g_xyyy_xxxy_0, g_xyyy_xxxy_1, g_xyyy_xxyy_0, g_xyyy_xxyy_1, g_xyyy_xxyz_0, g_xyyy_xxyz_1, g_xyyy_xyyy_0, g_xyyy_xyyy_1, g_xyyy_xyyz_0, g_xyyy_xyyz_1, g_xyyy_xyzz_0, g_xyyy_xyzz_1, g_xyyy_yyyy_0, g_xyyy_yyyy_1, g_xyyy_yyyz_0, g_xyyy_yyyz_1, g_xyyy_yyzz_0, g_xyyy_yyzz_1, g_xyyy_yzzz_0, g_xyyy_yzzz_1, g_xyyy_zzzz_0, g_xyyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyy_xxxx_0[i] = 2.0 * g_xxxy_xxxx_0[i] * fbe_0 - 2.0 * g_xxxy_xxxx_1[i] * fz_be_0 + g_xxxyy_xxxx_1[i] * pa_y[i];

        g_xxxyyy_xxxy_0[i] = 2.0 * g_xyyy_xxxy_0[i] * fbe_0 - 2.0 * g_xyyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xxyyy_xxy_1[i] * fe_0 + g_xxyyy_xxxy_1[i] * pa_x[i];

        g_xxxyyy_xxxz_0[i] = 2.0 * g_xxxy_xxxz_0[i] * fbe_0 - 2.0 * g_xxxy_xxxz_1[i] * fz_be_0 + g_xxxyy_xxxz_1[i] * pa_y[i];

        g_xxxyyy_xxyy_0[i] = 2.0 * g_xyyy_xxyy_0[i] * fbe_0 - 2.0 * g_xyyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyy_1[i] * fe_0 + g_xxyyy_xxyy_1[i] * pa_x[i];

        g_xxxyyy_xxyz_0[i] = 2.0 * g_xyyy_xxyz_0[i] * fbe_0 - 2.0 * g_xyyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyz_1[i] * fe_0 + g_xxyyy_xxyz_1[i] * pa_x[i];

        g_xxxyyy_xxzz_0[i] = 2.0 * g_xxxy_xxzz_0[i] * fbe_0 - 2.0 * g_xxxy_xxzz_1[i] * fz_be_0 + g_xxxyy_xxzz_1[i] * pa_y[i];

        g_xxxyyy_xyyy_0[i] = 2.0 * g_xyyy_xyyy_0[i] * fbe_0 - 2.0 * g_xyyy_xyyy_1[i] * fz_be_0 + g_xxyyy_yyy_1[i] * fe_0 + g_xxyyy_xyyy_1[i] * pa_x[i];

        g_xxxyyy_xyyz_0[i] = 2.0 * g_xyyy_xyyz_0[i] * fbe_0 - 2.0 * g_xyyy_xyyz_1[i] * fz_be_0 + g_xxyyy_yyz_1[i] * fe_0 + g_xxyyy_xyyz_1[i] * pa_x[i];

        g_xxxyyy_xyzz_0[i] = 2.0 * g_xyyy_xyzz_0[i] * fbe_0 - 2.0 * g_xyyy_xyzz_1[i] * fz_be_0 + g_xxyyy_yzz_1[i] * fe_0 + g_xxyyy_xyzz_1[i] * pa_x[i];

        g_xxxyyy_xzzz_0[i] = 2.0 * g_xxxy_xzzz_0[i] * fbe_0 - 2.0 * g_xxxy_xzzz_1[i] * fz_be_0 + g_xxxyy_xzzz_1[i] * pa_y[i];

        g_xxxyyy_yyyy_0[i] = 2.0 * g_xyyy_yyyy_0[i] * fbe_0 - 2.0 * g_xyyy_yyyy_1[i] * fz_be_0 + g_xxyyy_yyyy_1[i] * pa_x[i];

        g_xxxyyy_yyyz_0[i] = 2.0 * g_xyyy_yyyz_0[i] * fbe_0 - 2.0 * g_xyyy_yyyz_1[i] * fz_be_0 + g_xxyyy_yyyz_1[i] * pa_x[i];

        g_xxxyyy_yyzz_0[i] = 2.0 * g_xyyy_yyzz_0[i] * fbe_0 - 2.0 * g_xyyy_yyzz_1[i] * fz_be_0 + g_xxyyy_yyzz_1[i] * pa_x[i];

        g_xxxyyy_yzzz_0[i] = 2.0 * g_xyyy_yzzz_0[i] * fbe_0 - 2.0 * g_xyyy_yzzz_1[i] * fz_be_0 + g_xxyyy_yzzz_1[i] * pa_x[i];

        g_xxxyyy_zzzz_0[i] = 2.0 * g_xyyy_zzzz_0[i] * fbe_0 - 2.0 * g_xyyy_zzzz_1[i] * fz_be_0 + g_xxyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : IG

    auto g_xxxyyz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 105);

    auto g_xxxyyz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 106);

    auto g_xxxyyz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 107);

    auto g_xxxyyz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 108);

    auto g_xxxyyz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 109);

    auto g_xxxyyz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 110);

    auto g_xxxyyz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 111);

    auto g_xxxyyz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 112);

    auto g_xxxyyz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 113);

    auto g_xxxyyz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 114);

    auto g_xxxyyz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 115);

    auto g_xxxyyz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 116);

    auto g_xxxyyz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 117);

    auto g_xxxyyz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 118);

    auto g_xxxyyz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 119);

    #pragma omp simd aligned(g_xxxyy_xxx_1, g_xxxyy_xxxx_1, g_xxxyy_xxxy_1, g_xxxyy_xxxz_1, g_xxxyy_xxy_1, g_xxxyy_xxyy_1, g_xxxyy_xxyz_1, g_xxxyy_xxz_1, g_xxxyy_xxzz_1, g_xxxyy_xyy_1, g_xxxyy_xyyy_1, g_xxxyy_xyyz_1, g_xxxyy_xyz_1, g_xxxyy_xyzz_1, g_xxxyy_xzz_1, g_xxxyy_xzzz_1, g_xxxyy_yyy_1, g_xxxyy_yyyy_1, g_xxxyy_yyyz_1, g_xxxyy_yyz_1, g_xxxyy_yyzz_1, g_xxxyy_yzz_1, g_xxxyy_yzzz_1, g_xxxyy_zzz_1, g_xxxyy_zzzz_1, g_xxxyyz_xxxx_0, g_xxxyyz_xxxy_0, g_xxxyyz_xxxz_0, g_xxxyyz_xxyy_0, g_xxxyyz_xxyz_0, g_xxxyyz_xxzz_0, g_xxxyyz_xyyy_0, g_xxxyyz_xyyz_0, g_xxxyyz_xyzz_0, g_xxxyyz_xzzz_0, g_xxxyyz_yyyy_0, g_xxxyyz_yyyz_0, g_xxxyyz_yyzz_0, g_xxxyyz_yzzz_0, g_xxxyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyz_xxxx_0[i] = g_xxxyy_xxxx_1[i] * pa_z[i];

        g_xxxyyz_xxxy_0[i] = g_xxxyy_xxxy_1[i] * pa_z[i];

        g_xxxyyz_xxxz_0[i] = g_xxxyy_xxx_1[i] * fe_0 + g_xxxyy_xxxz_1[i] * pa_z[i];

        g_xxxyyz_xxyy_0[i] = g_xxxyy_xxyy_1[i] * pa_z[i];

        g_xxxyyz_xxyz_0[i] = g_xxxyy_xxy_1[i] * fe_0 + g_xxxyy_xxyz_1[i] * pa_z[i];

        g_xxxyyz_xxzz_0[i] = 2.0 * g_xxxyy_xxz_1[i] * fe_0 + g_xxxyy_xxzz_1[i] * pa_z[i];

        g_xxxyyz_xyyy_0[i] = g_xxxyy_xyyy_1[i] * pa_z[i];

        g_xxxyyz_xyyz_0[i] = g_xxxyy_xyy_1[i] * fe_0 + g_xxxyy_xyyz_1[i] * pa_z[i];

        g_xxxyyz_xyzz_0[i] = 2.0 * g_xxxyy_xyz_1[i] * fe_0 + g_xxxyy_xyzz_1[i] * pa_z[i];

        g_xxxyyz_xzzz_0[i] = 3.0 * g_xxxyy_xzz_1[i] * fe_0 + g_xxxyy_xzzz_1[i] * pa_z[i];

        g_xxxyyz_yyyy_0[i] = g_xxxyy_yyyy_1[i] * pa_z[i];

        g_xxxyyz_yyyz_0[i] = g_xxxyy_yyy_1[i] * fe_0 + g_xxxyy_yyyz_1[i] * pa_z[i];

        g_xxxyyz_yyzz_0[i] = 2.0 * g_xxxyy_yyz_1[i] * fe_0 + g_xxxyy_yyzz_1[i] * pa_z[i];

        g_xxxyyz_yzzz_0[i] = 3.0 * g_xxxyy_yzz_1[i] * fe_0 + g_xxxyy_yzzz_1[i] * pa_z[i];

        g_xxxyyz_zzzz_0[i] = 4.0 * g_xxxyy_zzz_1[i] * fe_0 + g_xxxyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 120-135 components of targeted buffer : IG

    auto g_xxxyzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 120);

    auto g_xxxyzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 121);

    auto g_xxxyzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 122);

    auto g_xxxyzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 123);

    auto g_xxxyzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 124);

    auto g_xxxyzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 125);

    auto g_xxxyzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 126);

    auto g_xxxyzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 127);

    auto g_xxxyzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 128);

    auto g_xxxyzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 129);

    auto g_xxxyzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 130);

    auto g_xxxyzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 131);

    auto g_xxxyzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 132);

    auto g_xxxyzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 133);

    auto g_xxxyzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 134);

    #pragma omp simd aligned(g_xxxyzz_xxxx_0, g_xxxyzz_xxxy_0, g_xxxyzz_xxxz_0, g_xxxyzz_xxyy_0, g_xxxyzz_xxyz_0, g_xxxyzz_xxzz_0, g_xxxyzz_xyyy_0, g_xxxyzz_xyyz_0, g_xxxyzz_xyzz_0, g_xxxyzz_xzzz_0, g_xxxyzz_yyyy_0, g_xxxyzz_yyyz_0, g_xxxyzz_yyzz_0, g_xxxyzz_yzzz_0, g_xxxyzz_zzzz_0, g_xxxzz_xxx_1, g_xxxzz_xxxx_1, g_xxxzz_xxxy_1, g_xxxzz_xxxz_1, g_xxxzz_xxy_1, g_xxxzz_xxyy_1, g_xxxzz_xxyz_1, g_xxxzz_xxz_1, g_xxxzz_xxzz_1, g_xxxzz_xyy_1, g_xxxzz_xyyy_1, g_xxxzz_xyyz_1, g_xxxzz_xyz_1, g_xxxzz_xyzz_1, g_xxxzz_xzz_1, g_xxxzz_xzzz_1, g_xxxzz_yyy_1, g_xxxzz_yyyy_1, g_xxxzz_yyyz_1, g_xxxzz_yyz_1, g_xxxzz_yyzz_1, g_xxxzz_yzz_1, g_xxxzz_yzzz_1, g_xxxzz_zzz_1, g_xxxzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzz_xxxx_0[i] = g_xxxzz_xxxx_1[i] * pa_y[i];

        g_xxxyzz_xxxy_0[i] = g_xxxzz_xxx_1[i] * fe_0 + g_xxxzz_xxxy_1[i] * pa_y[i];

        g_xxxyzz_xxxz_0[i] = g_xxxzz_xxxz_1[i] * pa_y[i];

        g_xxxyzz_xxyy_0[i] = 2.0 * g_xxxzz_xxy_1[i] * fe_0 + g_xxxzz_xxyy_1[i] * pa_y[i];

        g_xxxyzz_xxyz_0[i] = g_xxxzz_xxz_1[i] * fe_0 + g_xxxzz_xxyz_1[i] * pa_y[i];

        g_xxxyzz_xxzz_0[i] = g_xxxzz_xxzz_1[i] * pa_y[i];

        g_xxxyzz_xyyy_0[i] = 3.0 * g_xxxzz_xyy_1[i] * fe_0 + g_xxxzz_xyyy_1[i] * pa_y[i];

        g_xxxyzz_xyyz_0[i] = 2.0 * g_xxxzz_xyz_1[i] * fe_0 + g_xxxzz_xyyz_1[i] * pa_y[i];

        g_xxxyzz_xyzz_0[i] = g_xxxzz_xzz_1[i] * fe_0 + g_xxxzz_xyzz_1[i] * pa_y[i];

        g_xxxyzz_xzzz_0[i] = g_xxxzz_xzzz_1[i] * pa_y[i];

        g_xxxyzz_yyyy_0[i] = 4.0 * g_xxxzz_yyy_1[i] * fe_0 + g_xxxzz_yyyy_1[i] * pa_y[i];

        g_xxxyzz_yyyz_0[i] = 3.0 * g_xxxzz_yyz_1[i] * fe_0 + g_xxxzz_yyyz_1[i] * pa_y[i];

        g_xxxyzz_yyzz_0[i] = 2.0 * g_xxxzz_yzz_1[i] * fe_0 + g_xxxzz_yyzz_1[i] * pa_y[i];

        g_xxxyzz_yzzz_0[i] = g_xxxzz_zzz_1[i] * fe_0 + g_xxxzz_yzzz_1[i] * pa_y[i];

        g_xxxyzz_zzzz_0[i] = g_xxxzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : IG

    auto g_xxxzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 135);

    auto g_xxxzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 136);

    auto g_xxxzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 137);

    auto g_xxxzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 138);

    auto g_xxxzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 139);

    auto g_xxxzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 140);

    auto g_xxxzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 141);

    auto g_xxxzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 142);

    auto g_xxxzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 143);

    auto g_xxxzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 144);

    auto g_xxxzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 145);

    auto g_xxxzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 146);

    auto g_xxxzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 147);

    auto g_xxxzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 148);

    auto g_xxxzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 149);

    #pragma omp simd aligned(g_xxxz_xxxx_0, g_xxxz_xxxx_1, g_xxxz_xxxy_0, g_xxxz_xxxy_1, g_xxxz_xxyy_0, g_xxxz_xxyy_1, g_xxxz_xyyy_0, g_xxxz_xyyy_1, g_xxxzz_xxxx_1, g_xxxzz_xxxy_1, g_xxxzz_xxyy_1, g_xxxzz_xyyy_1, g_xxxzzz_xxxx_0, g_xxxzzz_xxxy_0, g_xxxzzz_xxxz_0, g_xxxzzz_xxyy_0, g_xxxzzz_xxyz_0, g_xxxzzz_xxzz_0, g_xxxzzz_xyyy_0, g_xxxzzz_xyyz_0, g_xxxzzz_xyzz_0, g_xxxzzz_xzzz_0, g_xxxzzz_yyyy_0, g_xxxzzz_yyyz_0, g_xxxzzz_yyzz_0, g_xxxzzz_yzzz_0, g_xxxzzz_zzzz_0, g_xxzzz_xxxz_1, g_xxzzz_xxyz_1, g_xxzzz_xxz_1, g_xxzzz_xxzz_1, g_xxzzz_xyyz_1, g_xxzzz_xyz_1, g_xxzzz_xyzz_1, g_xxzzz_xzz_1, g_xxzzz_xzzz_1, g_xxzzz_yyyy_1, g_xxzzz_yyyz_1, g_xxzzz_yyz_1, g_xxzzz_yyzz_1, g_xxzzz_yzz_1, g_xxzzz_yzzz_1, g_xxzzz_zzz_1, g_xxzzz_zzzz_1, g_xzzz_xxxz_0, g_xzzz_xxxz_1, g_xzzz_xxyz_0, g_xzzz_xxyz_1, g_xzzz_xxzz_0, g_xzzz_xxzz_1, g_xzzz_xyyz_0, g_xzzz_xyyz_1, g_xzzz_xyzz_0, g_xzzz_xyzz_1, g_xzzz_xzzz_0, g_xzzz_xzzz_1, g_xzzz_yyyy_0, g_xzzz_yyyy_1, g_xzzz_yyyz_0, g_xzzz_yyyz_1, g_xzzz_yyzz_0, g_xzzz_yyzz_1, g_xzzz_yzzz_0, g_xzzz_yzzz_1, g_xzzz_zzzz_0, g_xzzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzz_xxxx_0[i] = 2.0 * g_xxxz_xxxx_0[i] * fbe_0 - 2.0 * g_xxxz_xxxx_1[i] * fz_be_0 + g_xxxzz_xxxx_1[i] * pa_z[i];

        g_xxxzzz_xxxy_0[i] = 2.0 * g_xxxz_xxxy_0[i] * fbe_0 - 2.0 * g_xxxz_xxxy_1[i] * fz_be_0 + g_xxxzz_xxxy_1[i] * pa_z[i];

        g_xxxzzz_xxxz_0[i] = 2.0 * g_xzzz_xxxz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xxzzz_xxz_1[i] * fe_0 + g_xxzzz_xxxz_1[i] * pa_x[i];

        g_xxxzzz_xxyy_0[i] = 2.0 * g_xxxz_xxyy_0[i] * fbe_0 - 2.0 * g_xxxz_xxyy_1[i] * fz_be_0 + g_xxxzz_xxyy_1[i] * pa_z[i];

        g_xxxzzz_xxyz_0[i] = 2.0 * g_xzzz_xxyz_0[i] * fbe_0 - 2.0 * g_xzzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xyz_1[i] * fe_0 + g_xxzzz_xxyz_1[i] * pa_x[i];

        g_xxxzzz_xxzz_0[i] = 2.0 * g_xzzz_xxzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xzz_1[i] * fe_0 + g_xxzzz_xxzz_1[i] * pa_x[i];

        g_xxxzzz_xyyy_0[i] = 2.0 * g_xxxz_xyyy_0[i] * fbe_0 - 2.0 * g_xxxz_xyyy_1[i] * fz_be_0 + g_xxxzz_xyyy_1[i] * pa_z[i];

        g_xxxzzz_xyyz_0[i] = 2.0 * g_xzzz_xyyz_0[i] * fbe_0 - 2.0 * g_xzzz_xyyz_1[i] * fz_be_0 + g_xxzzz_yyz_1[i] * fe_0 + g_xxzzz_xyyz_1[i] * pa_x[i];

        g_xxxzzz_xyzz_0[i] = 2.0 * g_xzzz_xyzz_0[i] * fbe_0 - 2.0 * g_xzzz_xyzz_1[i] * fz_be_0 + g_xxzzz_yzz_1[i] * fe_0 + g_xxzzz_xyzz_1[i] * pa_x[i];

        g_xxxzzz_xzzz_0[i] = 2.0 * g_xzzz_xzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xzzz_1[i] * fz_be_0 + g_xxzzz_zzz_1[i] * fe_0 + g_xxzzz_xzzz_1[i] * pa_x[i];

        g_xxxzzz_yyyy_0[i] = 2.0 * g_xzzz_yyyy_0[i] * fbe_0 - 2.0 * g_xzzz_yyyy_1[i] * fz_be_0 + g_xxzzz_yyyy_1[i] * pa_x[i];

        g_xxxzzz_yyyz_0[i] = 2.0 * g_xzzz_yyyz_0[i] * fbe_0 - 2.0 * g_xzzz_yyyz_1[i] * fz_be_0 + g_xxzzz_yyyz_1[i] * pa_x[i];

        g_xxxzzz_yyzz_0[i] = 2.0 * g_xzzz_yyzz_0[i] * fbe_0 - 2.0 * g_xzzz_yyzz_1[i] * fz_be_0 + g_xxzzz_yyzz_1[i] * pa_x[i];

        g_xxxzzz_yzzz_0[i] = 2.0 * g_xzzz_yzzz_0[i] * fbe_0 - 2.0 * g_xzzz_yzzz_1[i] * fz_be_0 + g_xxzzz_yzzz_1[i] * pa_x[i];

        g_xxxzzz_zzzz_0[i] = 2.0 * g_xzzz_zzzz_0[i] * fbe_0 - 2.0 * g_xzzz_zzzz_1[i] * fz_be_0 + g_xxzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : IG

    auto g_xxyyyy_xxxx_0 = pbuffer.data(idx_eri_0_ig + 150);

    auto g_xxyyyy_xxxy_0 = pbuffer.data(idx_eri_0_ig + 151);

    auto g_xxyyyy_xxxz_0 = pbuffer.data(idx_eri_0_ig + 152);

    auto g_xxyyyy_xxyy_0 = pbuffer.data(idx_eri_0_ig + 153);

    auto g_xxyyyy_xxyz_0 = pbuffer.data(idx_eri_0_ig + 154);

    auto g_xxyyyy_xxzz_0 = pbuffer.data(idx_eri_0_ig + 155);

    auto g_xxyyyy_xyyy_0 = pbuffer.data(idx_eri_0_ig + 156);

    auto g_xxyyyy_xyyz_0 = pbuffer.data(idx_eri_0_ig + 157);

    auto g_xxyyyy_xyzz_0 = pbuffer.data(idx_eri_0_ig + 158);

    auto g_xxyyyy_xzzz_0 = pbuffer.data(idx_eri_0_ig + 159);

    auto g_xxyyyy_yyyy_0 = pbuffer.data(idx_eri_0_ig + 160);

    auto g_xxyyyy_yyyz_0 = pbuffer.data(idx_eri_0_ig + 161);

    auto g_xxyyyy_yyzz_0 = pbuffer.data(idx_eri_0_ig + 162);

    auto g_xxyyyy_yzzz_0 = pbuffer.data(idx_eri_0_ig + 163);

    auto g_xxyyyy_zzzz_0 = pbuffer.data(idx_eri_0_ig + 164);

    #pragma omp simd aligned(g_xxyy_xxxx_0, g_xxyy_xxxx_1, g_xxyy_xxxz_0, g_xxyy_xxxz_1, g_xxyy_xxzz_0, g_xxyy_xxzz_1, g_xxyy_xzzz_0, g_xxyy_xzzz_1, g_xxyyy_xxxx_1, g_xxyyy_xxxz_1, g_xxyyy_xxzz_1, g_xxyyy_xzzz_1, g_xxyyyy_xxxx_0, g_xxyyyy_xxxy_0, g_xxyyyy_xxxz_0, g_xxyyyy_xxyy_0, g_xxyyyy_xxyz_0, g_xxyyyy_xxzz_0, g_xxyyyy_xyyy_0, g_xxyyyy_xyyz_0, g_xxyyyy_xyzz_0, g_xxyyyy_xzzz_0, g_xxyyyy_yyyy_0, g_xxyyyy_yyyz_0, g_xxyyyy_yyzz_0, g_xxyyyy_yzzz_0, g_xxyyyy_zzzz_0, g_xyyyy_xxxy_1, g_xyyyy_xxy_1, g_xyyyy_xxyy_1, g_xyyyy_xxyz_1, g_xyyyy_xyy_1, g_xyyyy_xyyy_1, g_xyyyy_xyyz_1, g_xyyyy_xyz_1, g_xyyyy_xyzz_1, g_xyyyy_yyy_1, g_xyyyy_yyyy_1, g_xyyyy_yyyz_1, g_xyyyy_yyz_1, g_xyyyy_yyzz_1, g_xyyyy_yzz_1, g_xyyyy_yzzz_1, g_xyyyy_zzzz_1, g_yyyy_xxxy_0, g_yyyy_xxxy_1, g_yyyy_xxyy_0, g_yyyy_xxyy_1, g_yyyy_xxyz_0, g_yyyy_xxyz_1, g_yyyy_xyyy_0, g_yyyy_xyyy_1, g_yyyy_xyyz_0, g_yyyy_xyyz_1, g_yyyy_xyzz_0, g_yyyy_xyzz_1, g_yyyy_yyyy_0, g_yyyy_yyyy_1, g_yyyy_yyyz_0, g_yyyy_yyyz_1, g_yyyy_yyzz_0, g_yyyy_yyzz_1, g_yyyy_yzzz_0, g_yyyy_yzzz_1, g_yyyy_zzzz_0, g_yyyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyy_xxxx_0[i] = 3.0 * g_xxyy_xxxx_0[i] * fbe_0 - 3.0 * g_xxyy_xxxx_1[i] * fz_be_0 + g_xxyyy_xxxx_1[i] * pa_y[i];

        g_xxyyyy_xxxy_0[i] = g_yyyy_xxxy_0[i] * fbe_0 - g_yyyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xyyyy_xxy_1[i] * fe_0 + g_xyyyy_xxxy_1[i] * pa_x[i];

        g_xxyyyy_xxxz_0[i] = 3.0 * g_xxyy_xxxz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxz_1[i] * fz_be_0 + g_xxyyy_xxxz_1[i] * pa_y[i];

        g_xxyyyy_xxyy_0[i] = g_yyyy_xxyy_0[i] * fbe_0 - g_yyyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyy_1[i] * fe_0 + g_xyyyy_xxyy_1[i] * pa_x[i];

        g_xxyyyy_xxyz_0[i] = g_yyyy_xxyz_0[i] * fbe_0 - g_yyyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyz_1[i] * fe_0 + g_xyyyy_xxyz_1[i] * pa_x[i];

        g_xxyyyy_xxzz_0[i] = 3.0 * g_xxyy_xxzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxzz_1[i] * fz_be_0 + g_xxyyy_xxzz_1[i] * pa_y[i];

        g_xxyyyy_xyyy_0[i] = g_yyyy_xyyy_0[i] * fbe_0 - g_yyyy_xyyy_1[i] * fz_be_0 + g_xyyyy_yyy_1[i] * fe_0 + g_xyyyy_xyyy_1[i] * pa_x[i];

        g_xxyyyy_xyyz_0[i] = g_yyyy_xyyz_0[i] * fbe_0 - g_yyyy_xyyz_1[i] * fz_be_0 + g_xyyyy_yyz_1[i] * fe_0 + g_xyyyy_xyyz_1[i] * pa_x[i];

        g_xxyyyy_xyzz_0[i] = g_yyyy_xyzz_0[i] * fbe_0 - g_yyyy_xyzz_1[i] * fz_be_0 + g_xyyyy_yzz_1[i] * fe_0 + g_xyyyy_xyzz_1[i] * pa_x[i];

        g_xxyyyy_xzzz_0[i] = 3.0 * g_xxyy_xzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xzzz_1[i] * fz_be_0 + g_xxyyy_xzzz_1[i] * pa_y[i];

        g_xxyyyy_yyyy_0[i] = g_yyyy_yyyy_0[i] * fbe_0 - g_yyyy_yyyy_1[i] * fz_be_0 + g_xyyyy_yyyy_1[i] * pa_x[i];

        g_xxyyyy_yyyz_0[i] = g_yyyy_yyyz_0[i] * fbe_0 - g_yyyy_yyyz_1[i] * fz_be_0 + g_xyyyy_yyyz_1[i] * pa_x[i];

        g_xxyyyy_yyzz_0[i] = g_yyyy_yyzz_0[i] * fbe_0 - g_yyyy_yyzz_1[i] * fz_be_0 + g_xyyyy_yyzz_1[i] * pa_x[i];

        g_xxyyyy_yzzz_0[i] = g_yyyy_yzzz_0[i] * fbe_0 - g_yyyy_yzzz_1[i] * fz_be_0 + g_xyyyy_yzzz_1[i] * pa_x[i];

        g_xxyyyy_zzzz_0[i] = g_yyyy_zzzz_0[i] * fbe_0 - g_yyyy_zzzz_1[i] * fz_be_0 + g_xyyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 165-180 components of targeted buffer : IG

    auto g_xxyyyz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 165);

    auto g_xxyyyz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 166);

    auto g_xxyyyz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 167);

    auto g_xxyyyz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 168);

    auto g_xxyyyz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 169);

    auto g_xxyyyz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 170);

    auto g_xxyyyz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 171);

    auto g_xxyyyz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 172);

    auto g_xxyyyz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 173);

    auto g_xxyyyz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 174);

    auto g_xxyyyz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 175);

    auto g_xxyyyz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 176);

    auto g_xxyyyz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 177);

    auto g_xxyyyz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 178);

    auto g_xxyyyz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 179);

    #pragma omp simd aligned(g_xxyyy_xxx_1, g_xxyyy_xxxx_1, g_xxyyy_xxxy_1, g_xxyyy_xxxz_1, g_xxyyy_xxy_1, g_xxyyy_xxyy_1, g_xxyyy_xxyz_1, g_xxyyy_xxz_1, g_xxyyy_xxzz_1, g_xxyyy_xyy_1, g_xxyyy_xyyy_1, g_xxyyy_xyyz_1, g_xxyyy_xyz_1, g_xxyyy_xyzz_1, g_xxyyy_xzz_1, g_xxyyy_xzzz_1, g_xxyyy_yyy_1, g_xxyyy_yyyy_1, g_xxyyy_yyyz_1, g_xxyyy_yyz_1, g_xxyyy_yyzz_1, g_xxyyy_yzz_1, g_xxyyy_yzzz_1, g_xxyyy_zzz_1, g_xxyyy_zzzz_1, g_xxyyyz_xxxx_0, g_xxyyyz_xxxy_0, g_xxyyyz_xxxz_0, g_xxyyyz_xxyy_0, g_xxyyyz_xxyz_0, g_xxyyyz_xxzz_0, g_xxyyyz_xyyy_0, g_xxyyyz_xyyz_0, g_xxyyyz_xyzz_0, g_xxyyyz_xzzz_0, g_xxyyyz_yyyy_0, g_xxyyyz_yyyz_0, g_xxyyyz_yyzz_0, g_xxyyyz_yzzz_0, g_xxyyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyz_xxxx_0[i] = g_xxyyy_xxxx_1[i] * pa_z[i];

        g_xxyyyz_xxxy_0[i] = g_xxyyy_xxxy_1[i] * pa_z[i];

        g_xxyyyz_xxxz_0[i] = g_xxyyy_xxx_1[i] * fe_0 + g_xxyyy_xxxz_1[i] * pa_z[i];

        g_xxyyyz_xxyy_0[i] = g_xxyyy_xxyy_1[i] * pa_z[i];

        g_xxyyyz_xxyz_0[i] = g_xxyyy_xxy_1[i] * fe_0 + g_xxyyy_xxyz_1[i] * pa_z[i];

        g_xxyyyz_xxzz_0[i] = 2.0 * g_xxyyy_xxz_1[i] * fe_0 + g_xxyyy_xxzz_1[i] * pa_z[i];

        g_xxyyyz_xyyy_0[i] = g_xxyyy_xyyy_1[i] * pa_z[i];

        g_xxyyyz_xyyz_0[i] = g_xxyyy_xyy_1[i] * fe_0 + g_xxyyy_xyyz_1[i] * pa_z[i];

        g_xxyyyz_xyzz_0[i] = 2.0 * g_xxyyy_xyz_1[i] * fe_0 + g_xxyyy_xyzz_1[i] * pa_z[i];

        g_xxyyyz_xzzz_0[i] = 3.0 * g_xxyyy_xzz_1[i] * fe_0 + g_xxyyy_xzzz_1[i] * pa_z[i];

        g_xxyyyz_yyyy_0[i] = g_xxyyy_yyyy_1[i] * pa_z[i];

        g_xxyyyz_yyyz_0[i] = g_xxyyy_yyy_1[i] * fe_0 + g_xxyyy_yyyz_1[i] * pa_z[i];

        g_xxyyyz_yyzz_0[i] = 2.0 * g_xxyyy_yyz_1[i] * fe_0 + g_xxyyy_yyzz_1[i] * pa_z[i];

        g_xxyyyz_yzzz_0[i] = 3.0 * g_xxyyy_yzz_1[i] * fe_0 + g_xxyyy_yzzz_1[i] * pa_z[i];

        g_xxyyyz_zzzz_0[i] = 4.0 * g_xxyyy_zzz_1[i] * fe_0 + g_xxyyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 180-195 components of targeted buffer : IG

    auto g_xxyyzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 180);

    auto g_xxyyzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 181);

    auto g_xxyyzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 182);

    auto g_xxyyzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 183);

    auto g_xxyyzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 184);

    auto g_xxyyzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 185);

    auto g_xxyyzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 186);

    auto g_xxyyzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 187);

    auto g_xxyyzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 188);

    auto g_xxyyzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 189);

    auto g_xxyyzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 190);

    auto g_xxyyzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 191);

    auto g_xxyyzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 192);

    auto g_xxyyzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 193);

    auto g_xxyyzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 194);

    #pragma omp simd aligned(g_xxyy_xxxy_0, g_xxyy_xxxy_1, g_xxyy_xxyy_0, g_xxyy_xxyy_1, g_xxyy_xyyy_0, g_xxyy_xyyy_1, g_xxyyz_xxxy_1, g_xxyyz_xxyy_1, g_xxyyz_xyyy_1, g_xxyyzz_xxxx_0, g_xxyyzz_xxxy_0, g_xxyyzz_xxxz_0, g_xxyyzz_xxyy_0, g_xxyyzz_xxyz_0, g_xxyyzz_xxzz_0, g_xxyyzz_xyyy_0, g_xxyyzz_xyyz_0, g_xxyyzz_xyzz_0, g_xxyyzz_xzzz_0, g_xxyyzz_yyyy_0, g_xxyyzz_yyyz_0, g_xxyyzz_yyzz_0, g_xxyyzz_yzzz_0, g_xxyyzz_zzzz_0, g_xxyzz_xxxx_1, g_xxyzz_xxxz_1, g_xxyzz_xxzz_1, g_xxyzz_xzzz_1, g_xxzz_xxxx_0, g_xxzz_xxxx_1, g_xxzz_xxxz_0, g_xxzz_xxxz_1, g_xxzz_xxzz_0, g_xxzz_xxzz_1, g_xxzz_xzzz_0, g_xxzz_xzzz_1, g_xyyzz_xxyz_1, g_xyyzz_xyyz_1, g_xyyzz_xyz_1, g_xyyzz_xyzz_1, g_xyyzz_yyyy_1, g_xyyzz_yyyz_1, g_xyyzz_yyz_1, g_xyyzz_yyzz_1, g_xyyzz_yzz_1, g_xyyzz_yzzz_1, g_xyyzz_zzzz_1, g_yyzz_xxyz_0, g_yyzz_xxyz_1, g_yyzz_xyyz_0, g_yyzz_xyyz_1, g_yyzz_xyzz_0, g_yyzz_xyzz_1, g_yyzz_yyyy_0, g_yyzz_yyyy_1, g_yyzz_yyyz_0, g_yyzz_yyyz_1, g_yyzz_yyzz_0, g_yyzz_yyzz_1, g_yyzz_yzzz_0, g_yyzz_yzzz_1, g_yyzz_zzzz_0, g_yyzz_zzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzz_xxxx_0[i] = g_xxzz_xxxx_0[i] * fbe_0 - g_xxzz_xxxx_1[i] * fz_be_0 + g_xxyzz_xxxx_1[i] * pa_y[i];

        g_xxyyzz_xxxy_0[i] = g_xxyy_xxxy_0[i] * fbe_0 - g_xxyy_xxxy_1[i] * fz_be_0 + g_xxyyz_xxxy_1[i] * pa_z[i];

        g_xxyyzz_xxxz_0[i] = g_xxzz_xxxz_0[i] * fbe_0 - g_xxzz_xxxz_1[i] * fz_be_0 + g_xxyzz_xxxz_1[i] * pa_y[i];

        g_xxyyzz_xxyy_0[i] = g_xxyy_xxyy_0[i] * fbe_0 - g_xxyy_xxyy_1[i] * fz_be_0 + g_xxyyz_xxyy_1[i] * pa_z[i];

        g_xxyyzz_xxyz_0[i] = g_yyzz_xxyz_0[i] * fbe_0 - g_yyzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_xyz_1[i] * fe_0 + g_xyyzz_xxyz_1[i] * pa_x[i];

        g_xxyyzz_xxzz_0[i] = g_xxzz_xxzz_0[i] * fbe_0 - g_xxzz_xxzz_1[i] * fz_be_0 + g_xxyzz_xxzz_1[i] * pa_y[i];

        g_xxyyzz_xyyy_0[i] = g_xxyy_xyyy_0[i] * fbe_0 - g_xxyy_xyyy_1[i] * fz_be_0 + g_xxyyz_xyyy_1[i] * pa_z[i];

        g_xxyyzz_xyyz_0[i] = g_yyzz_xyyz_0[i] * fbe_0 - g_yyzz_xyyz_1[i] * fz_be_0 + g_xyyzz_yyz_1[i] * fe_0 + g_xyyzz_xyyz_1[i] * pa_x[i];

        g_xxyyzz_xyzz_0[i] = g_yyzz_xyzz_0[i] * fbe_0 - g_yyzz_xyzz_1[i] * fz_be_0 + g_xyyzz_yzz_1[i] * fe_0 + g_xyyzz_xyzz_1[i] * pa_x[i];

        g_xxyyzz_xzzz_0[i] = g_xxzz_xzzz_0[i] * fbe_0 - g_xxzz_xzzz_1[i] * fz_be_0 + g_xxyzz_xzzz_1[i] * pa_y[i];

        g_xxyyzz_yyyy_0[i] = g_yyzz_yyyy_0[i] * fbe_0 - g_yyzz_yyyy_1[i] * fz_be_0 + g_xyyzz_yyyy_1[i] * pa_x[i];

        g_xxyyzz_yyyz_0[i] = g_yyzz_yyyz_0[i] * fbe_0 - g_yyzz_yyyz_1[i] * fz_be_0 + g_xyyzz_yyyz_1[i] * pa_x[i];

        g_xxyyzz_yyzz_0[i] = g_yyzz_yyzz_0[i] * fbe_0 - g_yyzz_yyzz_1[i] * fz_be_0 + g_xyyzz_yyzz_1[i] * pa_x[i];

        g_xxyyzz_yzzz_0[i] = g_yyzz_yzzz_0[i] * fbe_0 - g_yyzz_yzzz_1[i] * fz_be_0 + g_xyyzz_yzzz_1[i] * pa_x[i];

        g_xxyyzz_zzzz_0[i] = g_yyzz_zzzz_0[i] * fbe_0 - g_yyzz_zzzz_1[i] * fz_be_0 + g_xyyzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 195-210 components of targeted buffer : IG

    auto g_xxyzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 195);

    auto g_xxyzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 196);

    auto g_xxyzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 197);

    auto g_xxyzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 198);

    auto g_xxyzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 199);

    auto g_xxyzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 200);

    auto g_xxyzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 201);

    auto g_xxyzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 202);

    auto g_xxyzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 203);

    auto g_xxyzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 204);

    auto g_xxyzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 205);

    auto g_xxyzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 206);

    auto g_xxyzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 207);

    auto g_xxyzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 208);

    auto g_xxyzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 209);

    #pragma omp simd aligned(g_xxyzzz_xxxx_0, g_xxyzzz_xxxy_0, g_xxyzzz_xxxz_0, g_xxyzzz_xxyy_0, g_xxyzzz_xxyz_0, g_xxyzzz_xxzz_0, g_xxyzzz_xyyy_0, g_xxyzzz_xyyz_0, g_xxyzzz_xyzz_0, g_xxyzzz_xzzz_0, g_xxyzzz_yyyy_0, g_xxyzzz_yyyz_0, g_xxyzzz_yyzz_0, g_xxyzzz_yzzz_0, g_xxyzzz_zzzz_0, g_xxzzz_xxx_1, g_xxzzz_xxxx_1, g_xxzzz_xxxy_1, g_xxzzz_xxxz_1, g_xxzzz_xxy_1, g_xxzzz_xxyy_1, g_xxzzz_xxyz_1, g_xxzzz_xxz_1, g_xxzzz_xxzz_1, g_xxzzz_xyy_1, g_xxzzz_xyyy_1, g_xxzzz_xyyz_1, g_xxzzz_xyz_1, g_xxzzz_xyzz_1, g_xxzzz_xzz_1, g_xxzzz_xzzz_1, g_xxzzz_yyy_1, g_xxzzz_yyyy_1, g_xxzzz_yyyz_1, g_xxzzz_yyz_1, g_xxzzz_yyzz_1, g_xxzzz_yzz_1, g_xxzzz_yzzz_1, g_xxzzz_zzz_1, g_xxzzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzz_xxxx_0[i] = g_xxzzz_xxxx_1[i] * pa_y[i];

        g_xxyzzz_xxxy_0[i] = g_xxzzz_xxx_1[i] * fe_0 + g_xxzzz_xxxy_1[i] * pa_y[i];

        g_xxyzzz_xxxz_0[i] = g_xxzzz_xxxz_1[i] * pa_y[i];

        g_xxyzzz_xxyy_0[i] = 2.0 * g_xxzzz_xxy_1[i] * fe_0 + g_xxzzz_xxyy_1[i] * pa_y[i];

        g_xxyzzz_xxyz_0[i] = g_xxzzz_xxz_1[i] * fe_0 + g_xxzzz_xxyz_1[i] * pa_y[i];

        g_xxyzzz_xxzz_0[i] = g_xxzzz_xxzz_1[i] * pa_y[i];

        g_xxyzzz_xyyy_0[i] = 3.0 * g_xxzzz_xyy_1[i] * fe_0 + g_xxzzz_xyyy_1[i] * pa_y[i];

        g_xxyzzz_xyyz_0[i] = 2.0 * g_xxzzz_xyz_1[i] * fe_0 + g_xxzzz_xyyz_1[i] * pa_y[i];

        g_xxyzzz_xyzz_0[i] = g_xxzzz_xzz_1[i] * fe_0 + g_xxzzz_xyzz_1[i] * pa_y[i];

        g_xxyzzz_xzzz_0[i] = g_xxzzz_xzzz_1[i] * pa_y[i];

        g_xxyzzz_yyyy_0[i] = 4.0 * g_xxzzz_yyy_1[i] * fe_0 + g_xxzzz_yyyy_1[i] * pa_y[i];

        g_xxyzzz_yyyz_0[i] = 3.0 * g_xxzzz_yyz_1[i] * fe_0 + g_xxzzz_yyyz_1[i] * pa_y[i];

        g_xxyzzz_yyzz_0[i] = 2.0 * g_xxzzz_yzz_1[i] * fe_0 + g_xxzzz_yyzz_1[i] * pa_y[i];

        g_xxyzzz_yzzz_0[i] = g_xxzzz_zzz_1[i] * fe_0 + g_xxzzz_yzzz_1[i] * pa_y[i];

        g_xxyzzz_zzzz_0[i] = g_xxzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 210-225 components of targeted buffer : IG

    auto g_xxzzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 210);

    auto g_xxzzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 211);

    auto g_xxzzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 212);

    auto g_xxzzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 213);

    auto g_xxzzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 214);

    auto g_xxzzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 215);

    auto g_xxzzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 216);

    auto g_xxzzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 217);

    auto g_xxzzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 218);

    auto g_xxzzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 219);

    auto g_xxzzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 220);

    auto g_xxzzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 221);

    auto g_xxzzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 222);

    auto g_xxzzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 223);

    auto g_xxzzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 224);

    #pragma omp simd aligned(g_xxzz_xxxx_0, g_xxzz_xxxx_1, g_xxzz_xxxy_0, g_xxzz_xxxy_1, g_xxzz_xxyy_0, g_xxzz_xxyy_1, g_xxzz_xyyy_0, g_xxzz_xyyy_1, g_xxzzz_xxxx_1, g_xxzzz_xxxy_1, g_xxzzz_xxyy_1, g_xxzzz_xyyy_1, g_xxzzzz_xxxx_0, g_xxzzzz_xxxy_0, g_xxzzzz_xxxz_0, g_xxzzzz_xxyy_0, g_xxzzzz_xxyz_0, g_xxzzzz_xxzz_0, g_xxzzzz_xyyy_0, g_xxzzzz_xyyz_0, g_xxzzzz_xyzz_0, g_xxzzzz_xzzz_0, g_xxzzzz_yyyy_0, g_xxzzzz_yyyz_0, g_xxzzzz_yyzz_0, g_xxzzzz_yzzz_0, g_xxzzzz_zzzz_0, g_xzzzz_xxxz_1, g_xzzzz_xxyz_1, g_xzzzz_xxz_1, g_xzzzz_xxzz_1, g_xzzzz_xyyz_1, g_xzzzz_xyz_1, g_xzzzz_xyzz_1, g_xzzzz_xzz_1, g_xzzzz_xzzz_1, g_xzzzz_yyyy_1, g_xzzzz_yyyz_1, g_xzzzz_yyz_1, g_xzzzz_yyzz_1, g_xzzzz_yzz_1, g_xzzzz_yzzz_1, g_xzzzz_zzz_1, g_xzzzz_zzzz_1, g_zzzz_xxxz_0, g_zzzz_xxxz_1, g_zzzz_xxyz_0, g_zzzz_xxyz_1, g_zzzz_xxzz_0, g_zzzz_xxzz_1, g_zzzz_xyyz_0, g_zzzz_xyyz_1, g_zzzz_xyzz_0, g_zzzz_xyzz_1, g_zzzz_xzzz_0, g_zzzz_xzzz_1, g_zzzz_yyyy_0, g_zzzz_yyyy_1, g_zzzz_yyyz_0, g_zzzz_yyyz_1, g_zzzz_yyzz_0, g_zzzz_yyzz_1, g_zzzz_yzzz_0, g_zzzz_yzzz_1, g_zzzz_zzzz_0, g_zzzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzz_xxxx_0[i] = 3.0 * g_xxzz_xxxx_0[i] * fbe_0 - 3.0 * g_xxzz_xxxx_1[i] * fz_be_0 + g_xxzzz_xxxx_1[i] * pa_z[i];

        g_xxzzzz_xxxy_0[i] = 3.0 * g_xxzz_xxxy_0[i] * fbe_0 - 3.0 * g_xxzz_xxxy_1[i] * fz_be_0 + g_xxzzz_xxxy_1[i] * pa_z[i];

        g_xxzzzz_xxxz_0[i] = g_zzzz_xxxz_0[i] * fbe_0 - g_zzzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xzzzz_xxz_1[i] * fe_0 + g_xzzzz_xxxz_1[i] * pa_x[i];

        g_xxzzzz_xxyy_0[i] = 3.0 * g_xxzz_xxyy_0[i] * fbe_0 - 3.0 * g_xxzz_xxyy_1[i] * fz_be_0 + g_xxzzz_xxyy_1[i] * pa_z[i];

        g_xxzzzz_xxyz_0[i] = g_zzzz_xxyz_0[i] * fbe_0 - g_zzzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xyz_1[i] * fe_0 + g_xzzzz_xxyz_1[i] * pa_x[i];

        g_xxzzzz_xxzz_0[i] = g_zzzz_xxzz_0[i] * fbe_0 - g_zzzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xzz_1[i] * fe_0 + g_xzzzz_xxzz_1[i] * pa_x[i];

        g_xxzzzz_xyyy_0[i] = 3.0 * g_xxzz_xyyy_0[i] * fbe_0 - 3.0 * g_xxzz_xyyy_1[i] * fz_be_0 + g_xxzzz_xyyy_1[i] * pa_z[i];

        g_xxzzzz_xyyz_0[i] = g_zzzz_xyyz_0[i] * fbe_0 - g_zzzz_xyyz_1[i] * fz_be_0 + g_xzzzz_yyz_1[i] * fe_0 + g_xzzzz_xyyz_1[i] * pa_x[i];

        g_xxzzzz_xyzz_0[i] = g_zzzz_xyzz_0[i] * fbe_0 - g_zzzz_xyzz_1[i] * fz_be_0 + g_xzzzz_yzz_1[i] * fe_0 + g_xzzzz_xyzz_1[i] * pa_x[i];

        g_xxzzzz_xzzz_0[i] = g_zzzz_xzzz_0[i] * fbe_0 - g_zzzz_xzzz_1[i] * fz_be_0 + g_xzzzz_zzz_1[i] * fe_0 + g_xzzzz_xzzz_1[i] * pa_x[i];

        g_xxzzzz_yyyy_0[i] = g_zzzz_yyyy_0[i] * fbe_0 - g_zzzz_yyyy_1[i] * fz_be_0 + g_xzzzz_yyyy_1[i] * pa_x[i];

        g_xxzzzz_yyyz_0[i] = g_zzzz_yyyz_0[i] * fbe_0 - g_zzzz_yyyz_1[i] * fz_be_0 + g_xzzzz_yyyz_1[i] * pa_x[i];

        g_xxzzzz_yyzz_0[i] = g_zzzz_yyzz_0[i] * fbe_0 - g_zzzz_yyzz_1[i] * fz_be_0 + g_xzzzz_yyzz_1[i] * pa_x[i];

        g_xxzzzz_yzzz_0[i] = g_zzzz_yzzz_0[i] * fbe_0 - g_zzzz_yzzz_1[i] * fz_be_0 + g_xzzzz_yzzz_1[i] * pa_x[i];

        g_xxzzzz_zzzz_0[i] = g_zzzz_zzzz_0[i] * fbe_0 - g_zzzz_zzzz_1[i] * fz_be_0 + g_xzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : IG

    auto g_xyyyyy_xxxx_0 = pbuffer.data(idx_eri_0_ig + 225);

    auto g_xyyyyy_xxxy_0 = pbuffer.data(idx_eri_0_ig + 226);

    auto g_xyyyyy_xxxz_0 = pbuffer.data(idx_eri_0_ig + 227);

    auto g_xyyyyy_xxyy_0 = pbuffer.data(idx_eri_0_ig + 228);

    auto g_xyyyyy_xxyz_0 = pbuffer.data(idx_eri_0_ig + 229);

    auto g_xyyyyy_xxzz_0 = pbuffer.data(idx_eri_0_ig + 230);

    auto g_xyyyyy_xyyy_0 = pbuffer.data(idx_eri_0_ig + 231);

    auto g_xyyyyy_xyyz_0 = pbuffer.data(idx_eri_0_ig + 232);

    auto g_xyyyyy_xyzz_0 = pbuffer.data(idx_eri_0_ig + 233);

    auto g_xyyyyy_xzzz_0 = pbuffer.data(idx_eri_0_ig + 234);

    auto g_xyyyyy_yyyy_0 = pbuffer.data(idx_eri_0_ig + 235);

    auto g_xyyyyy_yyyz_0 = pbuffer.data(idx_eri_0_ig + 236);

    auto g_xyyyyy_yyzz_0 = pbuffer.data(idx_eri_0_ig + 237);

    auto g_xyyyyy_yzzz_0 = pbuffer.data(idx_eri_0_ig + 238);

    auto g_xyyyyy_zzzz_0 = pbuffer.data(idx_eri_0_ig + 239);

    #pragma omp simd aligned(g_xyyyyy_xxxx_0, g_xyyyyy_xxxy_0, g_xyyyyy_xxxz_0, g_xyyyyy_xxyy_0, g_xyyyyy_xxyz_0, g_xyyyyy_xxzz_0, g_xyyyyy_xyyy_0, g_xyyyyy_xyyz_0, g_xyyyyy_xyzz_0, g_xyyyyy_xzzz_0, g_xyyyyy_yyyy_0, g_xyyyyy_yyyz_0, g_xyyyyy_yyzz_0, g_xyyyyy_yzzz_0, g_xyyyyy_zzzz_0, g_yyyyy_xxx_1, g_yyyyy_xxxx_1, g_yyyyy_xxxy_1, g_yyyyy_xxxz_1, g_yyyyy_xxy_1, g_yyyyy_xxyy_1, g_yyyyy_xxyz_1, g_yyyyy_xxz_1, g_yyyyy_xxzz_1, g_yyyyy_xyy_1, g_yyyyy_xyyy_1, g_yyyyy_xyyz_1, g_yyyyy_xyz_1, g_yyyyy_xyzz_1, g_yyyyy_xzz_1, g_yyyyy_xzzz_1, g_yyyyy_yyy_1, g_yyyyy_yyyy_1, g_yyyyy_yyyz_1, g_yyyyy_yyz_1, g_yyyyy_yyzz_1, g_yyyyy_yzz_1, g_yyyyy_yzzz_1, g_yyyyy_zzz_1, g_yyyyy_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyy_xxxx_0[i] = 4.0 * g_yyyyy_xxx_1[i] * fe_0 + g_yyyyy_xxxx_1[i] * pa_x[i];

        g_xyyyyy_xxxy_0[i] = 3.0 * g_yyyyy_xxy_1[i] * fe_0 + g_yyyyy_xxxy_1[i] * pa_x[i];

        g_xyyyyy_xxxz_0[i] = 3.0 * g_yyyyy_xxz_1[i] * fe_0 + g_yyyyy_xxxz_1[i] * pa_x[i];

        g_xyyyyy_xxyy_0[i] = 2.0 * g_yyyyy_xyy_1[i] * fe_0 + g_yyyyy_xxyy_1[i] * pa_x[i];

        g_xyyyyy_xxyz_0[i] = 2.0 * g_yyyyy_xyz_1[i] * fe_0 + g_yyyyy_xxyz_1[i] * pa_x[i];

        g_xyyyyy_xxzz_0[i] = 2.0 * g_yyyyy_xzz_1[i] * fe_0 + g_yyyyy_xxzz_1[i] * pa_x[i];

        g_xyyyyy_xyyy_0[i] = g_yyyyy_yyy_1[i] * fe_0 + g_yyyyy_xyyy_1[i] * pa_x[i];

        g_xyyyyy_xyyz_0[i] = g_yyyyy_yyz_1[i] * fe_0 + g_yyyyy_xyyz_1[i] * pa_x[i];

        g_xyyyyy_xyzz_0[i] = g_yyyyy_yzz_1[i] * fe_0 + g_yyyyy_xyzz_1[i] * pa_x[i];

        g_xyyyyy_xzzz_0[i] = g_yyyyy_zzz_1[i] * fe_0 + g_yyyyy_xzzz_1[i] * pa_x[i];

        g_xyyyyy_yyyy_0[i] = g_yyyyy_yyyy_1[i] * pa_x[i];

        g_xyyyyy_yyyz_0[i] = g_yyyyy_yyyz_1[i] * pa_x[i];

        g_xyyyyy_yyzz_0[i] = g_yyyyy_yyzz_1[i] * pa_x[i];

        g_xyyyyy_yzzz_0[i] = g_yyyyy_yzzz_1[i] * pa_x[i];

        g_xyyyyy_zzzz_0[i] = g_yyyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 240-255 components of targeted buffer : IG

    auto g_xyyyyz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 240);

    auto g_xyyyyz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 241);

    auto g_xyyyyz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 242);

    auto g_xyyyyz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 243);

    auto g_xyyyyz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 244);

    auto g_xyyyyz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 245);

    auto g_xyyyyz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 246);

    auto g_xyyyyz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 247);

    auto g_xyyyyz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 248);

    auto g_xyyyyz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 249);

    auto g_xyyyyz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 250);

    auto g_xyyyyz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 251);

    auto g_xyyyyz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 252);

    auto g_xyyyyz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 253);

    auto g_xyyyyz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 254);

    #pragma omp simd aligned(g_xyyyy_xxxx_1, g_xyyyy_xxxy_1, g_xyyyy_xxyy_1, g_xyyyy_xyyy_1, g_xyyyyz_xxxx_0, g_xyyyyz_xxxy_0, g_xyyyyz_xxxz_0, g_xyyyyz_xxyy_0, g_xyyyyz_xxyz_0, g_xyyyyz_xxzz_0, g_xyyyyz_xyyy_0, g_xyyyyz_xyyz_0, g_xyyyyz_xyzz_0, g_xyyyyz_xzzz_0, g_xyyyyz_yyyy_0, g_xyyyyz_yyyz_0, g_xyyyyz_yyzz_0, g_xyyyyz_yzzz_0, g_xyyyyz_zzzz_0, g_yyyyz_xxxz_1, g_yyyyz_xxyz_1, g_yyyyz_xxz_1, g_yyyyz_xxzz_1, g_yyyyz_xyyz_1, g_yyyyz_xyz_1, g_yyyyz_xyzz_1, g_yyyyz_xzz_1, g_yyyyz_xzzz_1, g_yyyyz_yyyy_1, g_yyyyz_yyyz_1, g_yyyyz_yyz_1, g_yyyyz_yyzz_1, g_yyyyz_yzz_1, g_yyyyz_yzzz_1, g_yyyyz_zzz_1, g_yyyyz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyz_xxxx_0[i] = g_xyyyy_xxxx_1[i] * pa_z[i];

        g_xyyyyz_xxxy_0[i] = g_xyyyy_xxxy_1[i] * pa_z[i];

        g_xyyyyz_xxxz_0[i] = 3.0 * g_yyyyz_xxz_1[i] * fe_0 + g_yyyyz_xxxz_1[i] * pa_x[i];

        g_xyyyyz_xxyy_0[i] = g_xyyyy_xxyy_1[i] * pa_z[i];

        g_xyyyyz_xxyz_0[i] = 2.0 * g_yyyyz_xyz_1[i] * fe_0 + g_yyyyz_xxyz_1[i] * pa_x[i];

        g_xyyyyz_xxzz_0[i] = 2.0 * g_yyyyz_xzz_1[i] * fe_0 + g_yyyyz_xxzz_1[i] * pa_x[i];

        g_xyyyyz_xyyy_0[i] = g_xyyyy_xyyy_1[i] * pa_z[i];

        g_xyyyyz_xyyz_0[i] = g_yyyyz_yyz_1[i] * fe_0 + g_yyyyz_xyyz_1[i] * pa_x[i];

        g_xyyyyz_xyzz_0[i] = g_yyyyz_yzz_1[i] * fe_0 + g_yyyyz_xyzz_1[i] * pa_x[i];

        g_xyyyyz_xzzz_0[i] = g_yyyyz_zzz_1[i] * fe_0 + g_yyyyz_xzzz_1[i] * pa_x[i];

        g_xyyyyz_yyyy_0[i] = g_yyyyz_yyyy_1[i] * pa_x[i];

        g_xyyyyz_yyyz_0[i] = g_yyyyz_yyyz_1[i] * pa_x[i];

        g_xyyyyz_yyzz_0[i] = g_yyyyz_yyzz_1[i] * pa_x[i];

        g_xyyyyz_yzzz_0[i] = g_yyyyz_yzzz_1[i] * pa_x[i];

        g_xyyyyz_zzzz_0[i] = g_yyyyz_zzzz_1[i] * pa_x[i];
    }

    // Set up 255-270 components of targeted buffer : IG

    auto g_xyyyzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 255);

    auto g_xyyyzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 256);

    auto g_xyyyzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 257);

    auto g_xyyyzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 258);

    auto g_xyyyzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 259);

    auto g_xyyyzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 260);

    auto g_xyyyzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 261);

    auto g_xyyyzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 262);

    auto g_xyyyzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 263);

    auto g_xyyyzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 264);

    auto g_xyyyzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 265);

    auto g_xyyyzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 266);

    auto g_xyyyzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 267);

    auto g_xyyyzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 268);

    auto g_xyyyzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 269);

    #pragma omp simd aligned(g_xyyyzz_xxxx_0, g_xyyyzz_xxxy_0, g_xyyyzz_xxxz_0, g_xyyyzz_xxyy_0, g_xyyyzz_xxyz_0, g_xyyyzz_xxzz_0, g_xyyyzz_xyyy_0, g_xyyyzz_xyyz_0, g_xyyyzz_xyzz_0, g_xyyyzz_xzzz_0, g_xyyyzz_yyyy_0, g_xyyyzz_yyyz_0, g_xyyyzz_yyzz_0, g_xyyyzz_yzzz_0, g_xyyyzz_zzzz_0, g_yyyzz_xxx_1, g_yyyzz_xxxx_1, g_yyyzz_xxxy_1, g_yyyzz_xxxz_1, g_yyyzz_xxy_1, g_yyyzz_xxyy_1, g_yyyzz_xxyz_1, g_yyyzz_xxz_1, g_yyyzz_xxzz_1, g_yyyzz_xyy_1, g_yyyzz_xyyy_1, g_yyyzz_xyyz_1, g_yyyzz_xyz_1, g_yyyzz_xyzz_1, g_yyyzz_xzz_1, g_yyyzz_xzzz_1, g_yyyzz_yyy_1, g_yyyzz_yyyy_1, g_yyyzz_yyyz_1, g_yyyzz_yyz_1, g_yyyzz_yyzz_1, g_yyyzz_yzz_1, g_yyyzz_yzzz_1, g_yyyzz_zzz_1, g_yyyzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzz_xxxx_0[i] = 4.0 * g_yyyzz_xxx_1[i] * fe_0 + g_yyyzz_xxxx_1[i] * pa_x[i];

        g_xyyyzz_xxxy_0[i] = 3.0 * g_yyyzz_xxy_1[i] * fe_0 + g_yyyzz_xxxy_1[i] * pa_x[i];

        g_xyyyzz_xxxz_0[i] = 3.0 * g_yyyzz_xxz_1[i] * fe_0 + g_yyyzz_xxxz_1[i] * pa_x[i];

        g_xyyyzz_xxyy_0[i] = 2.0 * g_yyyzz_xyy_1[i] * fe_0 + g_yyyzz_xxyy_1[i] * pa_x[i];

        g_xyyyzz_xxyz_0[i] = 2.0 * g_yyyzz_xyz_1[i] * fe_0 + g_yyyzz_xxyz_1[i] * pa_x[i];

        g_xyyyzz_xxzz_0[i] = 2.0 * g_yyyzz_xzz_1[i] * fe_0 + g_yyyzz_xxzz_1[i] * pa_x[i];

        g_xyyyzz_xyyy_0[i] = g_yyyzz_yyy_1[i] * fe_0 + g_yyyzz_xyyy_1[i] * pa_x[i];

        g_xyyyzz_xyyz_0[i] = g_yyyzz_yyz_1[i] * fe_0 + g_yyyzz_xyyz_1[i] * pa_x[i];

        g_xyyyzz_xyzz_0[i] = g_yyyzz_yzz_1[i] * fe_0 + g_yyyzz_xyzz_1[i] * pa_x[i];

        g_xyyyzz_xzzz_0[i] = g_yyyzz_zzz_1[i] * fe_0 + g_yyyzz_xzzz_1[i] * pa_x[i];

        g_xyyyzz_yyyy_0[i] = g_yyyzz_yyyy_1[i] * pa_x[i];

        g_xyyyzz_yyyz_0[i] = g_yyyzz_yyyz_1[i] * pa_x[i];

        g_xyyyzz_yyzz_0[i] = g_yyyzz_yyzz_1[i] * pa_x[i];

        g_xyyyzz_yzzz_0[i] = g_yyyzz_yzzz_1[i] * pa_x[i];

        g_xyyyzz_zzzz_0[i] = g_yyyzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 270-285 components of targeted buffer : IG

    auto g_xyyzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 270);

    auto g_xyyzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 271);

    auto g_xyyzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 272);

    auto g_xyyzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 273);

    auto g_xyyzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 274);

    auto g_xyyzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 275);

    auto g_xyyzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 276);

    auto g_xyyzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 277);

    auto g_xyyzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 278);

    auto g_xyyzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 279);

    auto g_xyyzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 280);

    auto g_xyyzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 281);

    auto g_xyyzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 282);

    auto g_xyyzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 283);

    auto g_xyyzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 284);

    #pragma omp simd aligned(g_xyyzzz_xxxx_0, g_xyyzzz_xxxy_0, g_xyyzzz_xxxz_0, g_xyyzzz_xxyy_0, g_xyyzzz_xxyz_0, g_xyyzzz_xxzz_0, g_xyyzzz_xyyy_0, g_xyyzzz_xyyz_0, g_xyyzzz_xyzz_0, g_xyyzzz_xzzz_0, g_xyyzzz_yyyy_0, g_xyyzzz_yyyz_0, g_xyyzzz_yyzz_0, g_xyyzzz_yzzz_0, g_xyyzzz_zzzz_0, g_yyzzz_xxx_1, g_yyzzz_xxxx_1, g_yyzzz_xxxy_1, g_yyzzz_xxxz_1, g_yyzzz_xxy_1, g_yyzzz_xxyy_1, g_yyzzz_xxyz_1, g_yyzzz_xxz_1, g_yyzzz_xxzz_1, g_yyzzz_xyy_1, g_yyzzz_xyyy_1, g_yyzzz_xyyz_1, g_yyzzz_xyz_1, g_yyzzz_xyzz_1, g_yyzzz_xzz_1, g_yyzzz_xzzz_1, g_yyzzz_yyy_1, g_yyzzz_yyyy_1, g_yyzzz_yyyz_1, g_yyzzz_yyz_1, g_yyzzz_yyzz_1, g_yyzzz_yzz_1, g_yyzzz_yzzz_1, g_yyzzz_zzz_1, g_yyzzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzz_xxxx_0[i] = 4.0 * g_yyzzz_xxx_1[i] * fe_0 + g_yyzzz_xxxx_1[i] * pa_x[i];

        g_xyyzzz_xxxy_0[i] = 3.0 * g_yyzzz_xxy_1[i] * fe_0 + g_yyzzz_xxxy_1[i] * pa_x[i];

        g_xyyzzz_xxxz_0[i] = 3.0 * g_yyzzz_xxz_1[i] * fe_0 + g_yyzzz_xxxz_1[i] * pa_x[i];

        g_xyyzzz_xxyy_0[i] = 2.0 * g_yyzzz_xyy_1[i] * fe_0 + g_yyzzz_xxyy_1[i] * pa_x[i];

        g_xyyzzz_xxyz_0[i] = 2.0 * g_yyzzz_xyz_1[i] * fe_0 + g_yyzzz_xxyz_1[i] * pa_x[i];

        g_xyyzzz_xxzz_0[i] = 2.0 * g_yyzzz_xzz_1[i] * fe_0 + g_yyzzz_xxzz_1[i] * pa_x[i];

        g_xyyzzz_xyyy_0[i] = g_yyzzz_yyy_1[i] * fe_0 + g_yyzzz_xyyy_1[i] * pa_x[i];

        g_xyyzzz_xyyz_0[i] = g_yyzzz_yyz_1[i] * fe_0 + g_yyzzz_xyyz_1[i] * pa_x[i];

        g_xyyzzz_xyzz_0[i] = g_yyzzz_yzz_1[i] * fe_0 + g_yyzzz_xyzz_1[i] * pa_x[i];

        g_xyyzzz_xzzz_0[i] = g_yyzzz_zzz_1[i] * fe_0 + g_yyzzz_xzzz_1[i] * pa_x[i];

        g_xyyzzz_yyyy_0[i] = g_yyzzz_yyyy_1[i] * pa_x[i];

        g_xyyzzz_yyyz_0[i] = g_yyzzz_yyyz_1[i] * pa_x[i];

        g_xyyzzz_yyzz_0[i] = g_yyzzz_yyzz_1[i] * pa_x[i];

        g_xyyzzz_yzzz_0[i] = g_yyzzz_yzzz_1[i] * pa_x[i];

        g_xyyzzz_zzzz_0[i] = g_yyzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 285-300 components of targeted buffer : IG

    auto g_xyzzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 285);

    auto g_xyzzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 286);

    auto g_xyzzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 287);

    auto g_xyzzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 288);

    auto g_xyzzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 289);

    auto g_xyzzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 290);

    auto g_xyzzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 291);

    auto g_xyzzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 292);

    auto g_xyzzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 293);

    auto g_xyzzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 294);

    auto g_xyzzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 295);

    auto g_xyzzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 296);

    auto g_xyzzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 297);

    auto g_xyzzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 298);

    auto g_xyzzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 299);

    #pragma omp simd aligned(g_xyzzzz_xxxx_0, g_xyzzzz_xxxy_0, g_xyzzzz_xxxz_0, g_xyzzzz_xxyy_0, g_xyzzzz_xxyz_0, g_xyzzzz_xxzz_0, g_xyzzzz_xyyy_0, g_xyzzzz_xyyz_0, g_xyzzzz_xyzz_0, g_xyzzzz_xzzz_0, g_xyzzzz_yyyy_0, g_xyzzzz_yyyz_0, g_xyzzzz_yyzz_0, g_xyzzzz_yzzz_0, g_xyzzzz_zzzz_0, g_xzzzz_xxxx_1, g_xzzzz_xxxz_1, g_xzzzz_xxzz_1, g_xzzzz_xzzz_1, g_yzzzz_xxxy_1, g_yzzzz_xxy_1, g_yzzzz_xxyy_1, g_yzzzz_xxyz_1, g_yzzzz_xyy_1, g_yzzzz_xyyy_1, g_yzzzz_xyyz_1, g_yzzzz_xyz_1, g_yzzzz_xyzz_1, g_yzzzz_yyy_1, g_yzzzz_yyyy_1, g_yzzzz_yyyz_1, g_yzzzz_yyz_1, g_yzzzz_yyzz_1, g_yzzzz_yzz_1, g_yzzzz_yzzz_1, g_yzzzz_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzz_xxxx_0[i] = g_xzzzz_xxxx_1[i] * pa_y[i];

        g_xyzzzz_xxxy_0[i] = 3.0 * g_yzzzz_xxy_1[i] * fe_0 + g_yzzzz_xxxy_1[i] * pa_x[i];

        g_xyzzzz_xxxz_0[i] = g_xzzzz_xxxz_1[i] * pa_y[i];

        g_xyzzzz_xxyy_0[i] = 2.0 * g_yzzzz_xyy_1[i] * fe_0 + g_yzzzz_xxyy_1[i] * pa_x[i];

        g_xyzzzz_xxyz_0[i] = 2.0 * g_yzzzz_xyz_1[i] * fe_0 + g_yzzzz_xxyz_1[i] * pa_x[i];

        g_xyzzzz_xxzz_0[i] = g_xzzzz_xxzz_1[i] * pa_y[i];

        g_xyzzzz_xyyy_0[i] = g_yzzzz_yyy_1[i] * fe_0 + g_yzzzz_xyyy_1[i] * pa_x[i];

        g_xyzzzz_xyyz_0[i] = g_yzzzz_yyz_1[i] * fe_0 + g_yzzzz_xyyz_1[i] * pa_x[i];

        g_xyzzzz_xyzz_0[i] = g_yzzzz_yzz_1[i] * fe_0 + g_yzzzz_xyzz_1[i] * pa_x[i];

        g_xyzzzz_xzzz_0[i] = g_xzzzz_xzzz_1[i] * pa_y[i];

        g_xyzzzz_yyyy_0[i] = g_yzzzz_yyyy_1[i] * pa_x[i];

        g_xyzzzz_yyyz_0[i] = g_yzzzz_yyyz_1[i] * pa_x[i];

        g_xyzzzz_yyzz_0[i] = g_yzzzz_yyzz_1[i] * pa_x[i];

        g_xyzzzz_yzzz_0[i] = g_yzzzz_yzzz_1[i] * pa_x[i];

        g_xyzzzz_zzzz_0[i] = g_yzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 300-315 components of targeted buffer : IG

    auto g_xzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 300);

    auto g_xzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 301);

    auto g_xzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 302);

    auto g_xzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 303);

    auto g_xzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 304);

    auto g_xzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 305);

    auto g_xzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 306);

    auto g_xzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 307);

    auto g_xzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 308);

    auto g_xzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 309);

    auto g_xzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 310);

    auto g_xzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 311);

    auto g_xzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 312);

    auto g_xzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 313);

    auto g_xzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 314);

    #pragma omp simd aligned(g_xzzzzz_xxxx_0, g_xzzzzz_xxxy_0, g_xzzzzz_xxxz_0, g_xzzzzz_xxyy_0, g_xzzzzz_xxyz_0, g_xzzzzz_xxzz_0, g_xzzzzz_xyyy_0, g_xzzzzz_xyyz_0, g_xzzzzz_xyzz_0, g_xzzzzz_xzzz_0, g_xzzzzz_yyyy_0, g_xzzzzz_yyyz_0, g_xzzzzz_yyzz_0, g_xzzzzz_yzzz_0, g_xzzzzz_zzzz_0, g_zzzzz_xxx_1, g_zzzzz_xxxx_1, g_zzzzz_xxxy_1, g_zzzzz_xxxz_1, g_zzzzz_xxy_1, g_zzzzz_xxyy_1, g_zzzzz_xxyz_1, g_zzzzz_xxz_1, g_zzzzz_xxzz_1, g_zzzzz_xyy_1, g_zzzzz_xyyy_1, g_zzzzz_xyyz_1, g_zzzzz_xyz_1, g_zzzzz_xyzz_1, g_zzzzz_xzz_1, g_zzzzz_xzzz_1, g_zzzzz_yyy_1, g_zzzzz_yyyy_1, g_zzzzz_yyyz_1, g_zzzzz_yyz_1, g_zzzzz_yyzz_1, g_zzzzz_yzz_1, g_zzzzz_yzzz_1, g_zzzzz_zzz_1, g_zzzzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzz_xxxx_0[i] = 4.0 * g_zzzzz_xxx_1[i] * fe_0 + g_zzzzz_xxxx_1[i] * pa_x[i];

        g_xzzzzz_xxxy_0[i] = 3.0 * g_zzzzz_xxy_1[i] * fe_0 + g_zzzzz_xxxy_1[i] * pa_x[i];

        g_xzzzzz_xxxz_0[i] = 3.0 * g_zzzzz_xxz_1[i] * fe_0 + g_zzzzz_xxxz_1[i] * pa_x[i];

        g_xzzzzz_xxyy_0[i] = 2.0 * g_zzzzz_xyy_1[i] * fe_0 + g_zzzzz_xxyy_1[i] * pa_x[i];

        g_xzzzzz_xxyz_0[i] = 2.0 * g_zzzzz_xyz_1[i] * fe_0 + g_zzzzz_xxyz_1[i] * pa_x[i];

        g_xzzzzz_xxzz_0[i] = 2.0 * g_zzzzz_xzz_1[i] * fe_0 + g_zzzzz_xxzz_1[i] * pa_x[i];

        g_xzzzzz_xyyy_0[i] = g_zzzzz_yyy_1[i] * fe_0 + g_zzzzz_xyyy_1[i] * pa_x[i];

        g_xzzzzz_xyyz_0[i] = g_zzzzz_yyz_1[i] * fe_0 + g_zzzzz_xyyz_1[i] * pa_x[i];

        g_xzzzzz_xyzz_0[i] = g_zzzzz_yzz_1[i] * fe_0 + g_zzzzz_xyzz_1[i] * pa_x[i];

        g_xzzzzz_xzzz_0[i] = g_zzzzz_zzz_1[i] * fe_0 + g_zzzzz_xzzz_1[i] * pa_x[i];

        g_xzzzzz_yyyy_0[i] = g_zzzzz_yyyy_1[i] * pa_x[i];

        g_xzzzzz_yyyz_0[i] = g_zzzzz_yyyz_1[i] * pa_x[i];

        g_xzzzzz_yyzz_0[i] = g_zzzzz_yyzz_1[i] * pa_x[i];

        g_xzzzzz_yzzz_0[i] = g_zzzzz_yzzz_1[i] * pa_x[i];

        g_xzzzzz_zzzz_0[i] = g_zzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 315-330 components of targeted buffer : IG

    auto g_yyyyyy_xxxx_0 = pbuffer.data(idx_eri_0_ig + 315);

    auto g_yyyyyy_xxxy_0 = pbuffer.data(idx_eri_0_ig + 316);

    auto g_yyyyyy_xxxz_0 = pbuffer.data(idx_eri_0_ig + 317);

    auto g_yyyyyy_xxyy_0 = pbuffer.data(idx_eri_0_ig + 318);

    auto g_yyyyyy_xxyz_0 = pbuffer.data(idx_eri_0_ig + 319);

    auto g_yyyyyy_xxzz_0 = pbuffer.data(idx_eri_0_ig + 320);

    auto g_yyyyyy_xyyy_0 = pbuffer.data(idx_eri_0_ig + 321);

    auto g_yyyyyy_xyyz_0 = pbuffer.data(idx_eri_0_ig + 322);

    auto g_yyyyyy_xyzz_0 = pbuffer.data(idx_eri_0_ig + 323);

    auto g_yyyyyy_xzzz_0 = pbuffer.data(idx_eri_0_ig + 324);

    auto g_yyyyyy_yyyy_0 = pbuffer.data(idx_eri_0_ig + 325);

    auto g_yyyyyy_yyyz_0 = pbuffer.data(idx_eri_0_ig + 326);

    auto g_yyyyyy_yyzz_0 = pbuffer.data(idx_eri_0_ig + 327);

    auto g_yyyyyy_yzzz_0 = pbuffer.data(idx_eri_0_ig + 328);

    auto g_yyyyyy_zzzz_0 = pbuffer.data(idx_eri_0_ig + 329);

    #pragma omp simd aligned(g_yyyy_xxxx_0, g_yyyy_xxxx_1, g_yyyy_xxxy_0, g_yyyy_xxxy_1, g_yyyy_xxxz_0, g_yyyy_xxxz_1, g_yyyy_xxyy_0, g_yyyy_xxyy_1, g_yyyy_xxyz_0, g_yyyy_xxyz_1, g_yyyy_xxzz_0, g_yyyy_xxzz_1, g_yyyy_xyyy_0, g_yyyy_xyyy_1, g_yyyy_xyyz_0, g_yyyy_xyyz_1, g_yyyy_xyzz_0, g_yyyy_xyzz_1, g_yyyy_xzzz_0, g_yyyy_xzzz_1, g_yyyy_yyyy_0, g_yyyy_yyyy_1, g_yyyy_yyyz_0, g_yyyy_yyyz_1, g_yyyy_yyzz_0, g_yyyy_yyzz_1, g_yyyy_yzzz_0, g_yyyy_yzzz_1, g_yyyy_zzzz_0, g_yyyy_zzzz_1, g_yyyyy_xxx_1, g_yyyyy_xxxx_1, g_yyyyy_xxxy_1, g_yyyyy_xxxz_1, g_yyyyy_xxy_1, g_yyyyy_xxyy_1, g_yyyyy_xxyz_1, g_yyyyy_xxz_1, g_yyyyy_xxzz_1, g_yyyyy_xyy_1, g_yyyyy_xyyy_1, g_yyyyy_xyyz_1, g_yyyyy_xyz_1, g_yyyyy_xyzz_1, g_yyyyy_xzz_1, g_yyyyy_xzzz_1, g_yyyyy_yyy_1, g_yyyyy_yyyy_1, g_yyyyy_yyyz_1, g_yyyyy_yyz_1, g_yyyyy_yyzz_1, g_yyyyy_yzz_1, g_yyyyy_yzzz_1, g_yyyyy_zzz_1, g_yyyyy_zzzz_1, g_yyyyyy_xxxx_0, g_yyyyyy_xxxy_0, g_yyyyyy_xxxz_0, g_yyyyyy_xxyy_0, g_yyyyyy_xxyz_0, g_yyyyyy_xxzz_0, g_yyyyyy_xyyy_0, g_yyyyyy_xyyz_0, g_yyyyyy_xyzz_0, g_yyyyyy_xzzz_0, g_yyyyyy_yyyy_0, g_yyyyyy_yyyz_0, g_yyyyyy_yyzz_0, g_yyyyyy_yzzz_0, g_yyyyyy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyy_xxxx_0[i] = 5.0 * g_yyyy_xxxx_0[i] * fbe_0 - 5.0 * g_yyyy_xxxx_1[i] * fz_be_0 + g_yyyyy_xxxx_1[i] * pa_y[i];

        g_yyyyyy_xxxy_0[i] = 5.0 * g_yyyy_xxxy_0[i] * fbe_0 - 5.0 * g_yyyy_xxxy_1[i] * fz_be_0 + g_yyyyy_xxx_1[i] * fe_0 + g_yyyyy_xxxy_1[i] * pa_y[i];

        g_yyyyyy_xxxz_0[i] = 5.0 * g_yyyy_xxxz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxz_1[i] * fz_be_0 + g_yyyyy_xxxz_1[i] * pa_y[i];

        g_yyyyyy_xxyy_0[i] = 5.0 * g_yyyy_xxyy_0[i] * fbe_0 - 5.0 * g_yyyy_xxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_xxy_1[i] * fe_0 + g_yyyyy_xxyy_1[i] * pa_y[i];

        g_yyyyyy_xxyz_0[i] = 5.0 * g_yyyy_xxyz_0[i] * fbe_0 - 5.0 * g_yyyy_xxyz_1[i] * fz_be_0 + g_yyyyy_xxz_1[i] * fe_0 + g_yyyyy_xxyz_1[i] * pa_y[i];

        g_yyyyyy_xxzz_0[i] = 5.0 * g_yyyy_xxzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxzz_1[i] * fz_be_0 + g_yyyyy_xxzz_1[i] * pa_y[i];

        g_yyyyyy_xyyy_0[i] = 5.0 * g_yyyy_xyyy_0[i] * fbe_0 - 5.0 * g_yyyy_xyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_xyy_1[i] * fe_0 + g_yyyyy_xyyy_1[i] * pa_y[i];

        g_yyyyyy_xyyz_0[i] = 5.0 * g_yyyy_xyyz_0[i] * fbe_0 - 5.0 * g_yyyy_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_xyz_1[i] * fe_0 + g_yyyyy_xyyz_1[i] * pa_y[i];

        g_yyyyyy_xyzz_0[i] = 5.0 * g_yyyy_xyzz_0[i] * fbe_0 - 5.0 * g_yyyy_xyzz_1[i] * fz_be_0 + g_yyyyy_xzz_1[i] * fe_0 + g_yyyyy_xyzz_1[i] * pa_y[i];

        g_yyyyyy_xzzz_0[i] = 5.0 * g_yyyy_xzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xzzz_1[i] * fz_be_0 + g_yyyyy_xzzz_1[i] * pa_y[i];

        g_yyyyyy_yyyy_0[i] = 5.0 * g_yyyy_yyyy_0[i] * fbe_0 - 5.0 * g_yyyy_yyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_yyy_1[i] * fe_0 + g_yyyyy_yyyy_1[i] * pa_y[i];

        g_yyyyyy_yyyz_0[i] = 5.0 * g_yyyy_yyyz_0[i] * fbe_0 - 5.0 * g_yyyy_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_yyz_1[i] * fe_0 + g_yyyyy_yyyz_1[i] * pa_y[i];

        g_yyyyyy_yyzz_0[i] = 5.0 * g_yyyy_yyzz_0[i] * fbe_0 - 5.0 * g_yyyy_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_yzz_1[i] * fe_0 + g_yyyyy_yyzz_1[i] * pa_y[i];

        g_yyyyyy_yzzz_0[i] = 5.0 * g_yyyy_yzzz_0[i] * fbe_0 - 5.0 * g_yyyy_yzzz_1[i] * fz_be_0 + g_yyyyy_zzz_1[i] * fe_0 + g_yyyyy_yzzz_1[i] * pa_y[i];

        g_yyyyyy_zzzz_0[i] = 5.0 * g_yyyy_zzzz_0[i] * fbe_0 - 5.0 * g_yyyy_zzzz_1[i] * fz_be_0 + g_yyyyy_zzzz_1[i] * pa_y[i];
    }

    // Set up 330-345 components of targeted buffer : IG

    auto g_yyyyyz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 330);

    auto g_yyyyyz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 331);

    auto g_yyyyyz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 332);

    auto g_yyyyyz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 333);

    auto g_yyyyyz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 334);

    auto g_yyyyyz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 335);

    auto g_yyyyyz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 336);

    auto g_yyyyyz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 337);

    auto g_yyyyyz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 338);

    auto g_yyyyyz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 339);

    auto g_yyyyyz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 340);

    auto g_yyyyyz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 341);

    auto g_yyyyyz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 342);

    auto g_yyyyyz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 343);

    auto g_yyyyyz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 344);

    #pragma omp simd aligned(g_yyyyy_xxx_1, g_yyyyy_xxxx_1, g_yyyyy_xxxy_1, g_yyyyy_xxxz_1, g_yyyyy_xxy_1, g_yyyyy_xxyy_1, g_yyyyy_xxyz_1, g_yyyyy_xxz_1, g_yyyyy_xxzz_1, g_yyyyy_xyy_1, g_yyyyy_xyyy_1, g_yyyyy_xyyz_1, g_yyyyy_xyz_1, g_yyyyy_xyzz_1, g_yyyyy_xzz_1, g_yyyyy_xzzz_1, g_yyyyy_yyy_1, g_yyyyy_yyyy_1, g_yyyyy_yyyz_1, g_yyyyy_yyz_1, g_yyyyy_yyzz_1, g_yyyyy_yzz_1, g_yyyyy_yzzz_1, g_yyyyy_zzz_1, g_yyyyy_zzzz_1, g_yyyyyz_xxxx_0, g_yyyyyz_xxxy_0, g_yyyyyz_xxxz_0, g_yyyyyz_xxyy_0, g_yyyyyz_xxyz_0, g_yyyyyz_xxzz_0, g_yyyyyz_xyyy_0, g_yyyyyz_xyyz_0, g_yyyyyz_xyzz_0, g_yyyyyz_xzzz_0, g_yyyyyz_yyyy_0, g_yyyyyz_yyyz_0, g_yyyyyz_yyzz_0, g_yyyyyz_yzzz_0, g_yyyyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyz_xxxx_0[i] = g_yyyyy_xxxx_1[i] * pa_z[i];

        g_yyyyyz_xxxy_0[i] = g_yyyyy_xxxy_1[i] * pa_z[i];

        g_yyyyyz_xxxz_0[i] = g_yyyyy_xxx_1[i] * fe_0 + g_yyyyy_xxxz_1[i] * pa_z[i];

        g_yyyyyz_xxyy_0[i] = g_yyyyy_xxyy_1[i] * pa_z[i];

        g_yyyyyz_xxyz_0[i] = g_yyyyy_xxy_1[i] * fe_0 + g_yyyyy_xxyz_1[i] * pa_z[i];

        g_yyyyyz_xxzz_0[i] = 2.0 * g_yyyyy_xxz_1[i] * fe_0 + g_yyyyy_xxzz_1[i] * pa_z[i];

        g_yyyyyz_xyyy_0[i] = g_yyyyy_xyyy_1[i] * pa_z[i];

        g_yyyyyz_xyyz_0[i] = g_yyyyy_xyy_1[i] * fe_0 + g_yyyyy_xyyz_1[i] * pa_z[i];

        g_yyyyyz_xyzz_0[i] = 2.0 * g_yyyyy_xyz_1[i] * fe_0 + g_yyyyy_xyzz_1[i] * pa_z[i];

        g_yyyyyz_xzzz_0[i] = 3.0 * g_yyyyy_xzz_1[i] * fe_0 + g_yyyyy_xzzz_1[i] * pa_z[i];

        g_yyyyyz_yyyy_0[i] = g_yyyyy_yyyy_1[i] * pa_z[i];

        g_yyyyyz_yyyz_0[i] = g_yyyyy_yyy_1[i] * fe_0 + g_yyyyy_yyyz_1[i] * pa_z[i];

        g_yyyyyz_yyzz_0[i] = 2.0 * g_yyyyy_yyz_1[i] * fe_0 + g_yyyyy_yyzz_1[i] * pa_z[i];

        g_yyyyyz_yzzz_0[i] = 3.0 * g_yyyyy_yzz_1[i] * fe_0 + g_yyyyy_yzzz_1[i] * pa_z[i];

        g_yyyyyz_zzzz_0[i] = 4.0 * g_yyyyy_zzz_1[i] * fe_0 + g_yyyyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 345-360 components of targeted buffer : IG

    auto g_yyyyzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 345);

    auto g_yyyyzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 346);

    auto g_yyyyzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 347);

    auto g_yyyyzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 348);

    auto g_yyyyzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 349);

    auto g_yyyyzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 350);

    auto g_yyyyzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 351);

    auto g_yyyyzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 352);

    auto g_yyyyzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 353);

    auto g_yyyyzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 354);

    auto g_yyyyzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 355);

    auto g_yyyyzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 356);

    auto g_yyyyzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 357);

    auto g_yyyyzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 358);

    auto g_yyyyzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 359);

    #pragma omp simd aligned(g_yyyy_xxxy_0, g_yyyy_xxxy_1, g_yyyy_xxyy_0, g_yyyy_xxyy_1, g_yyyy_xyyy_0, g_yyyy_xyyy_1, g_yyyy_yyyy_0, g_yyyy_yyyy_1, g_yyyyz_xxxy_1, g_yyyyz_xxyy_1, g_yyyyz_xyyy_1, g_yyyyz_yyyy_1, g_yyyyzz_xxxx_0, g_yyyyzz_xxxy_0, g_yyyyzz_xxxz_0, g_yyyyzz_xxyy_0, g_yyyyzz_xxyz_0, g_yyyyzz_xxzz_0, g_yyyyzz_xyyy_0, g_yyyyzz_xyyz_0, g_yyyyzz_xyzz_0, g_yyyyzz_xzzz_0, g_yyyyzz_yyyy_0, g_yyyyzz_yyyz_0, g_yyyyzz_yyzz_0, g_yyyyzz_yzzz_0, g_yyyyzz_zzzz_0, g_yyyzz_xxxx_1, g_yyyzz_xxxz_1, g_yyyzz_xxyz_1, g_yyyzz_xxz_1, g_yyyzz_xxzz_1, g_yyyzz_xyyz_1, g_yyyzz_xyz_1, g_yyyzz_xyzz_1, g_yyyzz_xzz_1, g_yyyzz_xzzz_1, g_yyyzz_yyyz_1, g_yyyzz_yyz_1, g_yyyzz_yyzz_1, g_yyyzz_yzz_1, g_yyyzz_yzzz_1, g_yyyzz_zzz_1, g_yyyzz_zzzz_1, g_yyzz_xxxx_0, g_yyzz_xxxx_1, g_yyzz_xxxz_0, g_yyzz_xxxz_1, g_yyzz_xxyz_0, g_yyzz_xxyz_1, g_yyzz_xxzz_0, g_yyzz_xxzz_1, g_yyzz_xyyz_0, g_yyzz_xyyz_1, g_yyzz_xyzz_0, g_yyzz_xyzz_1, g_yyzz_xzzz_0, g_yyzz_xzzz_1, g_yyzz_yyyz_0, g_yyzz_yyyz_1, g_yyzz_yyzz_0, g_yyzz_yyzz_1, g_yyzz_yzzz_0, g_yyzz_yzzz_1, g_yyzz_zzzz_0, g_yyzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzz_xxxx_0[i] = 3.0 * g_yyzz_xxxx_0[i] * fbe_0 - 3.0 * g_yyzz_xxxx_1[i] * fz_be_0 + g_yyyzz_xxxx_1[i] * pa_y[i];

        g_yyyyzz_xxxy_0[i] = g_yyyy_xxxy_0[i] * fbe_0 - g_yyyy_xxxy_1[i] * fz_be_0 + g_yyyyz_xxxy_1[i] * pa_z[i];

        g_yyyyzz_xxxz_0[i] = 3.0 * g_yyzz_xxxz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxz_1[i] * fz_be_0 + g_yyyzz_xxxz_1[i] * pa_y[i];

        g_yyyyzz_xxyy_0[i] = g_yyyy_xxyy_0[i] * fbe_0 - g_yyyy_xxyy_1[i] * fz_be_0 + g_yyyyz_xxyy_1[i] * pa_z[i];

        g_yyyyzz_xxyz_0[i] = 3.0 * g_yyzz_xxyz_0[i] * fbe_0 - 3.0 * g_yyzz_xxyz_1[i] * fz_be_0 + g_yyyzz_xxz_1[i] * fe_0 + g_yyyzz_xxyz_1[i] * pa_y[i];

        g_yyyyzz_xxzz_0[i] = 3.0 * g_yyzz_xxzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxzz_1[i] * fz_be_0 + g_yyyzz_xxzz_1[i] * pa_y[i];

        g_yyyyzz_xyyy_0[i] = g_yyyy_xyyy_0[i] * fbe_0 - g_yyyy_xyyy_1[i] * fz_be_0 + g_yyyyz_xyyy_1[i] * pa_z[i];

        g_yyyyzz_xyyz_0[i] = 3.0 * g_yyzz_xyyz_0[i] * fbe_0 - 3.0 * g_yyzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_xyz_1[i] * fe_0 + g_yyyzz_xyyz_1[i] * pa_y[i];

        g_yyyyzz_xyzz_0[i] = 3.0 * g_yyzz_xyzz_0[i] * fbe_0 - 3.0 * g_yyzz_xyzz_1[i] * fz_be_0 + g_yyyzz_xzz_1[i] * fe_0 + g_yyyzz_xyzz_1[i] * pa_y[i];

        g_yyyyzz_xzzz_0[i] = 3.0 * g_yyzz_xzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xzzz_1[i] * fz_be_0 + g_yyyzz_xzzz_1[i] * pa_y[i];

        g_yyyyzz_yyyy_0[i] = g_yyyy_yyyy_0[i] * fbe_0 - g_yyyy_yyyy_1[i] * fz_be_0 + g_yyyyz_yyyy_1[i] * pa_z[i];

        g_yyyyzz_yyyz_0[i] = 3.0 * g_yyzz_yyyz_0[i] * fbe_0 - 3.0 * g_yyzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_yyz_1[i] * fe_0 + g_yyyzz_yyyz_1[i] * pa_y[i];

        g_yyyyzz_yyzz_0[i] = 3.0 * g_yyzz_yyzz_0[i] * fbe_0 - 3.0 * g_yyzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_yzz_1[i] * fe_0 + g_yyyzz_yyzz_1[i] * pa_y[i];

        g_yyyyzz_yzzz_0[i] = 3.0 * g_yyzz_yzzz_0[i] * fbe_0 - 3.0 * g_yyzz_yzzz_1[i] * fz_be_0 + g_yyyzz_zzz_1[i] * fe_0 + g_yyyzz_yzzz_1[i] * pa_y[i];

        g_yyyyzz_zzzz_0[i] = 3.0 * g_yyzz_zzzz_0[i] * fbe_0 - 3.0 * g_yyzz_zzzz_1[i] * fz_be_0 + g_yyyzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 360-375 components of targeted buffer : IG

    auto g_yyyzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 360);

    auto g_yyyzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 361);

    auto g_yyyzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 362);

    auto g_yyyzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 363);

    auto g_yyyzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 364);

    auto g_yyyzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 365);

    auto g_yyyzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 366);

    auto g_yyyzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 367);

    auto g_yyyzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 368);

    auto g_yyyzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 369);

    auto g_yyyzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 370);

    auto g_yyyzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 371);

    auto g_yyyzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 372);

    auto g_yyyzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 373);

    auto g_yyyzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 374);

    #pragma omp simd aligned(g_yyyz_xxxy_0, g_yyyz_xxxy_1, g_yyyz_xxyy_0, g_yyyz_xxyy_1, g_yyyz_xyyy_0, g_yyyz_xyyy_1, g_yyyz_yyyy_0, g_yyyz_yyyy_1, g_yyyzz_xxxy_1, g_yyyzz_xxyy_1, g_yyyzz_xyyy_1, g_yyyzz_yyyy_1, g_yyyzzz_xxxx_0, g_yyyzzz_xxxy_0, g_yyyzzz_xxxz_0, g_yyyzzz_xxyy_0, g_yyyzzz_xxyz_0, g_yyyzzz_xxzz_0, g_yyyzzz_xyyy_0, g_yyyzzz_xyyz_0, g_yyyzzz_xyzz_0, g_yyyzzz_xzzz_0, g_yyyzzz_yyyy_0, g_yyyzzz_yyyz_0, g_yyyzzz_yyzz_0, g_yyyzzz_yzzz_0, g_yyyzzz_zzzz_0, g_yyzzz_xxxx_1, g_yyzzz_xxxz_1, g_yyzzz_xxyz_1, g_yyzzz_xxz_1, g_yyzzz_xxzz_1, g_yyzzz_xyyz_1, g_yyzzz_xyz_1, g_yyzzz_xyzz_1, g_yyzzz_xzz_1, g_yyzzz_xzzz_1, g_yyzzz_yyyz_1, g_yyzzz_yyz_1, g_yyzzz_yyzz_1, g_yyzzz_yzz_1, g_yyzzz_yzzz_1, g_yyzzz_zzz_1, g_yyzzz_zzzz_1, g_yzzz_xxxx_0, g_yzzz_xxxx_1, g_yzzz_xxxz_0, g_yzzz_xxxz_1, g_yzzz_xxyz_0, g_yzzz_xxyz_1, g_yzzz_xxzz_0, g_yzzz_xxzz_1, g_yzzz_xyyz_0, g_yzzz_xyyz_1, g_yzzz_xyzz_0, g_yzzz_xyzz_1, g_yzzz_xzzz_0, g_yzzz_xzzz_1, g_yzzz_yyyz_0, g_yzzz_yyyz_1, g_yzzz_yyzz_0, g_yzzz_yyzz_1, g_yzzz_yzzz_0, g_yzzz_yzzz_1, g_yzzz_zzzz_0, g_yzzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzz_xxxx_0[i] = 2.0 * g_yzzz_xxxx_0[i] * fbe_0 - 2.0 * g_yzzz_xxxx_1[i] * fz_be_0 + g_yyzzz_xxxx_1[i] * pa_y[i];

        g_yyyzzz_xxxy_0[i] = 2.0 * g_yyyz_xxxy_0[i] * fbe_0 - 2.0 * g_yyyz_xxxy_1[i] * fz_be_0 + g_yyyzz_xxxy_1[i] * pa_z[i];

        g_yyyzzz_xxxz_0[i] = 2.0 * g_yzzz_xxxz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxz_1[i] * fz_be_0 + g_yyzzz_xxxz_1[i] * pa_y[i];

        g_yyyzzz_xxyy_0[i] = 2.0 * g_yyyz_xxyy_0[i] * fbe_0 - 2.0 * g_yyyz_xxyy_1[i] * fz_be_0 + g_yyyzz_xxyy_1[i] * pa_z[i];

        g_yyyzzz_xxyz_0[i] = 2.0 * g_yzzz_xxyz_0[i] * fbe_0 - 2.0 * g_yzzz_xxyz_1[i] * fz_be_0 + g_yyzzz_xxz_1[i] * fe_0 + g_yyzzz_xxyz_1[i] * pa_y[i];

        g_yyyzzz_xxzz_0[i] = 2.0 * g_yzzz_xxzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxzz_1[i] * fz_be_0 + g_yyzzz_xxzz_1[i] * pa_y[i];

        g_yyyzzz_xyyy_0[i] = 2.0 * g_yyyz_xyyy_0[i] * fbe_0 - 2.0 * g_yyyz_xyyy_1[i] * fz_be_0 + g_yyyzz_xyyy_1[i] * pa_z[i];

        g_yyyzzz_xyyz_0[i] = 2.0 * g_yzzz_xyyz_0[i] * fbe_0 - 2.0 * g_yzzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_xyz_1[i] * fe_0 + g_yyzzz_xyyz_1[i] * pa_y[i];

        g_yyyzzz_xyzz_0[i] = 2.0 * g_yzzz_xyzz_0[i] * fbe_0 - 2.0 * g_yzzz_xyzz_1[i] * fz_be_0 + g_yyzzz_xzz_1[i] * fe_0 + g_yyzzz_xyzz_1[i] * pa_y[i];

        g_yyyzzz_xzzz_0[i] = 2.0 * g_yzzz_xzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xzzz_1[i] * fz_be_0 + g_yyzzz_xzzz_1[i] * pa_y[i];

        g_yyyzzz_yyyy_0[i] = 2.0 * g_yyyz_yyyy_0[i] * fbe_0 - 2.0 * g_yyyz_yyyy_1[i] * fz_be_0 + g_yyyzz_yyyy_1[i] * pa_z[i];

        g_yyyzzz_yyyz_0[i] = 2.0 * g_yzzz_yyyz_0[i] * fbe_0 - 2.0 * g_yzzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_yyz_1[i] * fe_0 + g_yyzzz_yyyz_1[i] * pa_y[i];

        g_yyyzzz_yyzz_0[i] = 2.0 * g_yzzz_yyzz_0[i] * fbe_0 - 2.0 * g_yzzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_yzz_1[i] * fe_0 + g_yyzzz_yyzz_1[i] * pa_y[i];

        g_yyyzzz_yzzz_0[i] = 2.0 * g_yzzz_yzzz_0[i] * fbe_0 - 2.0 * g_yzzz_yzzz_1[i] * fz_be_0 + g_yyzzz_zzz_1[i] * fe_0 + g_yyzzz_yzzz_1[i] * pa_y[i];

        g_yyyzzz_zzzz_0[i] = 2.0 * g_yzzz_zzzz_0[i] * fbe_0 - 2.0 * g_yzzz_zzzz_1[i] * fz_be_0 + g_yyzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 375-390 components of targeted buffer : IG

    auto g_yyzzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 375);

    auto g_yyzzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 376);

    auto g_yyzzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 377);

    auto g_yyzzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 378);

    auto g_yyzzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 379);

    auto g_yyzzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 380);

    auto g_yyzzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 381);

    auto g_yyzzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 382);

    auto g_yyzzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 383);

    auto g_yyzzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 384);

    auto g_yyzzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 385);

    auto g_yyzzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 386);

    auto g_yyzzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 387);

    auto g_yyzzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 388);

    auto g_yyzzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 389);

    #pragma omp simd aligned(g_yyzz_xxxy_0, g_yyzz_xxxy_1, g_yyzz_xxyy_0, g_yyzz_xxyy_1, g_yyzz_xyyy_0, g_yyzz_xyyy_1, g_yyzz_yyyy_0, g_yyzz_yyyy_1, g_yyzzz_xxxy_1, g_yyzzz_xxyy_1, g_yyzzz_xyyy_1, g_yyzzz_yyyy_1, g_yyzzzz_xxxx_0, g_yyzzzz_xxxy_0, g_yyzzzz_xxxz_0, g_yyzzzz_xxyy_0, g_yyzzzz_xxyz_0, g_yyzzzz_xxzz_0, g_yyzzzz_xyyy_0, g_yyzzzz_xyyz_0, g_yyzzzz_xyzz_0, g_yyzzzz_xzzz_0, g_yyzzzz_yyyy_0, g_yyzzzz_yyyz_0, g_yyzzzz_yyzz_0, g_yyzzzz_yzzz_0, g_yyzzzz_zzzz_0, g_yzzzz_xxxx_1, g_yzzzz_xxxz_1, g_yzzzz_xxyz_1, g_yzzzz_xxz_1, g_yzzzz_xxzz_1, g_yzzzz_xyyz_1, g_yzzzz_xyz_1, g_yzzzz_xyzz_1, g_yzzzz_xzz_1, g_yzzzz_xzzz_1, g_yzzzz_yyyz_1, g_yzzzz_yyz_1, g_yzzzz_yyzz_1, g_yzzzz_yzz_1, g_yzzzz_yzzz_1, g_yzzzz_zzz_1, g_yzzzz_zzzz_1, g_zzzz_xxxx_0, g_zzzz_xxxx_1, g_zzzz_xxxz_0, g_zzzz_xxxz_1, g_zzzz_xxyz_0, g_zzzz_xxyz_1, g_zzzz_xxzz_0, g_zzzz_xxzz_1, g_zzzz_xyyz_0, g_zzzz_xyyz_1, g_zzzz_xyzz_0, g_zzzz_xyzz_1, g_zzzz_xzzz_0, g_zzzz_xzzz_1, g_zzzz_yyyz_0, g_zzzz_yyyz_1, g_zzzz_yyzz_0, g_zzzz_yyzz_1, g_zzzz_yzzz_0, g_zzzz_yzzz_1, g_zzzz_zzzz_0, g_zzzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzz_xxxx_0[i] = g_zzzz_xxxx_0[i] * fbe_0 - g_zzzz_xxxx_1[i] * fz_be_0 + g_yzzzz_xxxx_1[i] * pa_y[i];

        g_yyzzzz_xxxy_0[i] = 3.0 * g_yyzz_xxxy_0[i] * fbe_0 - 3.0 * g_yyzz_xxxy_1[i] * fz_be_0 + g_yyzzz_xxxy_1[i] * pa_z[i];

        g_yyzzzz_xxxz_0[i] = g_zzzz_xxxz_0[i] * fbe_0 - g_zzzz_xxxz_1[i] * fz_be_0 + g_yzzzz_xxxz_1[i] * pa_y[i];

        g_yyzzzz_xxyy_0[i] = 3.0 * g_yyzz_xxyy_0[i] * fbe_0 - 3.0 * g_yyzz_xxyy_1[i] * fz_be_0 + g_yyzzz_xxyy_1[i] * pa_z[i];

        g_yyzzzz_xxyz_0[i] = g_zzzz_xxyz_0[i] * fbe_0 - g_zzzz_xxyz_1[i] * fz_be_0 + g_yzzzz_xxz_1[i] * fe_0 + g_yzzzz_xxyz_1[i] * pa_y[i];

        g_yyzzzz_xxzz_0[i] = g_zzzz_xxzz_0[i] * fbe_0 - g_zzzz_xxzz_1[i] * fz_be_0 + g_yzzzz_xxzz_1[i] * pa_y[i];

        g_yyzzzz_xyyy_0[i] = 3.0 * g_yyzz_xyyy_0[i] * fbe_0 - 3.0 * g_yyzz_xyyy_1[i] * fz_be_0 + g_yyzzz_xyyy_1[i] * pa_z[i];

        g_yyzzzz_xyyz_0[i] = g_zzzz_xyyz_0[i] * fbe_0 - g_zzzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_xyz_1[i] * fe_0 + g_yzzzz_xyyz_1[i] * pa_y[i];

        g_yyzzzz_xyzz_0[i] = g_zzzz_xyzz_0[i] * fbe_0 - g_zzzz_xyzz_1[i] * fz_be_0 + g_yzzzz_xzz_1[i] * fe_0 + g_yzzzz_xyzz_1[i] * pa_y[i];

        g_yyzzzz_xzzz_0[i] = g_zzzz_xzzz_0[i] * fbe_0 - g_zzzz_xzzz_1[i] * fz_be_0 + g_yzzzz_xzzz_1[i] * pa_y[i];

        g_yyzzzz_yyyy_0[i] = 3.0 * g_yyzz_yyyy_0[i] * fbe_0 - 3.0 * g_yyzz_yyyy_1[i] * fz_be_0 + g_yyzzz_yyyy_1[i] * pa_z[i];

        g_yyzzzz_yyyz_0[i] = g_zzzz_yyyz_0[i] * fbe_0 - g_zzzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_yyz_1[i] * fe_0 + g_yzzzz_yyyz_1[i] * pa_y[i];

        g_yyzzzz_yyzz_0[i] = g_zzzz_yyzz_0[i] * fbe_0 - g_zzzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_yzz_1[i] * fe_0 + g_yzzzz_yyzz_1[i] * pa_y[i];

        g_yyzzzz_yzzz_0[i] = g_zzzz_yzzz_0[i] * fbe_0 - g_zzzz_yzzz_1[i] * fz_be_0 + g_yzzzz_zzz_1[i] * fe_0 + g_yzzzz_yzzz_1[i] * pa_y[i];

        g_yyzzzz_zzzz_0[i] = g_zzzz_zzzz_0[i] * fbe_0 - g_zzzz_zzzz_1[i] * fz_be_0 + g_yzzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 390-405 components of targeted buffer : IG

    auto g_yzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 390);

    auto g_yzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 391);

    auto g_yzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 392);

    auto g_yzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 393);

    auto g_yzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 394);

    auto g_yzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 395);

    auto g_yzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 396);

    auto g_yzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 397);

    auto g_yzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 398);

    auto g_yzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 399);

    auto g_yzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 400);

    auto g_yzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 401);

    auto g_yzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 402);

    auto g_yzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 403);

    auto g_yzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 404);

    #pragma omp simd aligned(g_yzzzzz_xxxx_0, g_yzzzzz_xxxy_0, g_yzzzzz_xxxz_0, g_yzzzzz_xxyy_0, g_yzzzzz_xxyz_0, g_yzzzzz_xxzz_0, g_yzzzzz_xyyy_0, g_yzzzzz_xyyz_0, g_yzzzzz_xyzz_0, g_yzzzzz_xzzz_0, g_yzzzzz_yyyy_0, g_yzzzzz_yyyz_0, g_yzzzzz_yyzz_0, g_yzzzzz_yzzz_0, g_yzzzzz_zzzz_0, g_zzzzz_xxx_1, g_zzzzz_xxxx_1, g_zzzzz_xxxy_1, g_zzzzz_xxxz_1, g_zzzzz_xxy_1, g_zzzzz_xxyy_1, g_zzzzz_xxyz_1, g_zzzzz_xxz_1, g_zzzzz_xxzz_1, g_zzzzz_xyy_1, g_zzzzz_xyyy_1, g_zzzzz_xyyz_1, g_zzzzz_xyz_1, g_zzzzz_xyzz_1, g_zzzzz_xzz_1, g_zzzzz_xzzz_1, g_zzzzz_yyy_1, g_zzzzz_yyyy_1, g_zzzzz_yyyz_1, g_zzzzz_yyz_1, g_zzzzz_yyzz_1, g_zzzzz_yzz_1, g_zzzzz_yzzz_1, g_zzzzz_zzz_1, g_zzzzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzz_xxxx_0[i] = g_zzzzz_xxxx_1[i] * pa_y[i];

        g_yzzzzz_xxxy_0[i] = g_zzzzz_xxx_1[i] * fe_0 + g_zzzzz_xxxy_1[i] * pa_y[i];

        g_yzzzzz_xxxz_0[i] = g_zzzzz_xxxz_1[i] * pa_y[i];

        g_yzzzzz_xxyy_0[i] = 2.0 * g_zzzzz_xxy_1[i] * fe_0 + g_zzzzz_xxyy_1[i] * pa_y[i];

        g_yzzzzz_xxyz_0[i] = g_zzzzz_xxz_1[i] * fe_0 + g_zzzzz_xxyz_1[i] * pa_y[i];

        g_yzzzzz_xxzz_0[i] = g_zzzzz_xxzz_1[i] * pa_y[i];

        g_yzzzzz_xyyy_0[i] = 3.0 * g_zzzzz_xyy_1[i] * fe_0 + g_zzzzz_xyyy_1[i] * pa_y[i];

        g_yzzzzz_xyyz_0[i] = 2.0 * g_zzzzz_xyz_1[i] * fe_0 + g_zzzzz_xyyz_1[i] * pa_y[i];

        g_yzzzzz_xyzz_0[i] = g_zzzzz_xzz_1[i] * fe_0 + g_zzzzz_xyzz_1[i] * pa_y[i];

        g_yzzzzz_xzzz_0[i] = g_zzzzz_xzzz_1[i] * pa_y[i];

        g_yzzzzz_yyyy_0[i] = 4.0 * g_zzzzz_yyy_1[i] * fe_0 + g_zzzzz_yyyy_1[i] * pa_y[i];

        g_yzzzzz_yyyz_0[i] = 3.0 * g_zzzzz_yyz_1[i] * fe_0 + g_zzzzz_yyyz_1[i] * pa_y[i];

        g_yzzzzz_yyzz_0[i] = 2.0 * g_zzzzz_yzz_1[i] * fe_0 + g_zzzzz_yyzz_1[i] * pa_y[i];

        g_yzzzzz_yzzz_0[i] = g_zzzzz_zzz_1[i] * fe_0 + g_zzzzz_yzzz_1[i] * pa_y[i];

        g_yzzzzz_zzzz_0[i] = g_zzzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 405-420 components of targeted buffer : IG

    auto g_zzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_ig + 405);

    auto g_zzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_ig + 406);

    auto g_zzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_ig + 407);

    auto g_zzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_ig + 408);

    auto g_zzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_ig + 409);

    auto g_zzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_ig + 410);

    auto g_zzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_ig + 411);

    auto g_zzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_ig + 412);

    auto g_zzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_ig + 413);

    auto g_zzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_ig + 414);

    auto g_zzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_ig + 415);

    auto g_zzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_ig + 416);

    auto g_zzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_ig + 417);

    auto g_zzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_ig + 418);

    auto g_zzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_ig + 419);

    #pragma omp simd aligned(g_zzzz_xxxx_0, g_zzzz_xxxx_1, g_zzzz_xxxy_0, g_zzzz_xxxy_1, g_zzzz_xxxz_0, g_zzzz_xxxz_1, g_zzzz_xxyy_0, g_zzzz_xxyy_1, g_zzzz_xxyz_0, g_zzzz_xxyz_1, g_zzzz_xxzz_0, g_zzzz_xxzz_1, g_zzzz_xyyy_0, g_zzzz_xyyy_1, g_zzzz_xyyz_0, g_zzzz_xyyz_1, g_zzzz_xyzz_0, g_zzzz_xyzz_1, g_zzzz_xzzz_0, g_zzzz_xzzz_1, g_zzzz_yyyy_0, g_zzzz_yyyy_1, g_zzzz_yyyz_0, g_zzzz_yyyz_1, g_zzzz_yyzz_0, g_zzzz_yyzz_1, g_zzzz_yzzz_0, g_zzzz_yzzz_1, g_zzzz_zzzz_0, g_zzzz_zzzz_1, g_zzzzz_xxx_1, g_zzzzz_xxxx_1, g_zzzzz_xxxy_1, g_zzzzz_xxxz_1, g_zzzzz_xxy_1, g_zzzzz_xxyy_1, g_zzzzz_xxyz_1, g_zzzzz_xxz_1, g_zzzzz_xxzz_1, g_zzzzz_xyy_1, g_zzzzz_xyyy_1, g_zzzzz_xyyz_1, g_zzzzz_xyz_1, g_zzzzz_xyzz_1, g_zzzzz_xzz_1, g_zzzzz_xzzz_1, g_zzzzz_yyy_1, g_zzzzz_yyyy_1, g_zzzzz_yyyz_1, g_zzzzz_yyz_1, g_zzzzz_yyzz_1, g_zzzzz_yzz_1, g_zzzzz_yzzz_1, g_zzzzz_zzz_1, g_zzzzz_zzzz_1, g_zzzzzz_xxxx_0, g_zzzzzz_xxxy_0, g_zzzzzz_xxxz_0, g_zzzzzz_xxyy_0, g_zzzzzz_xxyz_0, g_zzzzzz_xxzz_0, g_zzzzzz_xyyy_0, g_zzzzzz_xyyz_0, g_zzzzzz_xyzz_0, g_zzzzzz_xzzz_0, g_zzzzzz_yyyy_0, g_zzzzzz_yyyz_0, g_zzzzzz_yyzz_0, g_zzzzzz_yzzz_0, g_zzzzzz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzz_xxxx_0[i] = 5.0 * g_zzzz_xxxx_0[i] * fbe_0 - 5.0 * g_zzzz_xxxx_1[i] * fz_be_0 + g_zzzzz_xxxx_1[i] * pa_z[i];

        g_zzzzzz_xxxy_0[i] = 5.0 * g_zzzz_xxxy_0[i] * fbe_0 - 5.0 * g_zzzz_xxxy_1[i] * fz_be_0 + g_zzzzz_xxxy_1[i] * pa_z[i];

        g_zzzzzz_xxxz_0[i] = 5.0 * g_zzzz_xxxz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxz_1[i] * fz_be_0 + g_zzzzz_xxx_1[i] * fe_0 + g_zzzzz_xxxz_1[i] * pa_z[i];

        g_zzzzzz_xxyy_0[i] = 5.0 * g_zzzz_xxyy_0[i] * fbe_0 - 5.0 * g_zzzz_xxyy_1[i] * fz_be_0 + g_zzzzz_xxyy_1[i] * pa_z[i];

        g_zzzzzz_xxyz_0[i] = 5.0 * g_zzzz_xxyz_0[i] * fbe_0 - 5.0 * g_zzzz_xxyz_1[i] * fz_be_0 + g_zzzzz_xxy_1[i] * fe_0 + g_zzzzz_xxyz_1[i] * pa_z[i];

        g_zzzzzz_xxzz_0[i] = 5.0 * g_zzzz_xxzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xxz_1[i] * fe_0 + g_zzzzz_xxzz_1[i] * pa_z[i];

        g_zzzzzz_xyyy_0[i] = 5.0 * g_zzzz_xyyy_0[i] * fbe_0 - 5.0 * g_zzzz_xyyy_1[i] * fz_be_0 + g_zzzzz_xyyy_1[i] * pa_z[i];

        g_zzzzzz_xyyz_0[i] = 5.0 * g_zzzz_xyyz_0[i] * fbe_0 - 5.0 * g_zzzz_xyyz_1[i] * fz_be_0 + g_zzzzz_xyy_1[i] * fe_0 + g_zzzzz_xyyz_1[i] * pa_z[i];

        g_zzzzzz_xyzz_0[i] = 5.0 * g_zzzz_xyzz_0[i] * fbe_0 - 5.0 * g_zzzz_xyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xyz_1[i] * fe_0 + g_zzzzz_xyzz_1[i] * pa_z[i];

        g_zzzzzz_xzzz_0[i] = 5.0 * g_zzzz_xzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_xzz_1[i] * fe_0 + g_zzzzz_xzzz_1[i] * pa_z[i];

        g_zzzzzz_yyyy_0[i] = 5.0 * g_zzzz_yyyy_0[i] * fbe_0 - 5.0 * g_zzzz_yyyy_1[i] * fz_be_0 + g_zzzzz_yyyy_1[i] * pa_z[i];

        g_zzzzzz_yyyz_0[i] = 5.0 * g_zzzz_yyyz_0[i] * fbe_0 - 5.0 * g_zzzz_yyyz_1[i] * fz_be_0 + g_zzzzz_yyy_1[i] * fe_0 + g_zzzzz_yyyz_1[i] * pa_z[i];

        g_zzzzzz_yyzz_0[i] = 5.0 * g_zzzz_yyzz_0[i] * fbe_0 - 5.0 * g_zzzz_yyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_yyz_1[i] * fe_0 + g_zzzzz_yyzz_1[i] * pa_z[i];

        g_zzzzzz_yzzz_0[i] = 5.0 * g_zzzz_yzzz_0[i] * fbe_0 - 5.0 * g_zzzz_yzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_yzz_1[i] * fe_0 + g_zzzzz_yzzz_1[i] * pa_z[i];

        g_zzzzzz_zzzz_0[i] = 5.0 * g_zzzz_zzzz_0[i] * fbe_0 - 5.0 * g_zzzz_zzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_zzz_1[i] * fe_0 + g_zzzzz_zzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

