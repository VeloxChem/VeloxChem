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

#include "ThreeCenterElectronRepulsionPrimRecKSF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksf,
                                 size_t idx_eri_0_hsf,
                                 size_t idx_eri_1_hsf,
                                 size_t idx_eri_1_isd,
                                 size_t idx_eri_1_isf,
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

    /// Set up components of auxilary buffer : HSF

    auto g_xxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_hsf);

    auto g_xxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 1);

    auto g_xxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 2);

    auto g_xxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 3);

    auto g_xxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 4);

    auto g_xxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 5);

    auto g_xxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 6);

    auto g_xxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 7);

    auto g_xxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 8);

    auto g_xxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 9);

    auto g_xxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 10);

    auto g_xxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 12);

    auto g_xxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 15);

    auto g_xxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 20);

    auto g_xxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 21);

    auto g_xxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 23);

    auto g_xxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 30);

    auto g_xxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 31);

    auto g_xxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 32);

    auto g_xxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 33);

    auto g_xxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 34);

    auto g_xxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 35);

    auto g_xxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 36);

    auto g_xxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 37);

    auto g_xxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 38);

    auto g_xxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 39);

    auto g_xxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 50);

    auto g_xxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 51);

    auto g_xxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 52);

    auto g_xxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 53);

    auto g_xxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 54);

    auto g_xxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 55);

    auto g_xxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 56);

    auto g_xxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 57);

    auto g_xxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 58);

    auto g_xxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 59);

    auto g_xxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 60);

    auto g_xxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 61);

    auto g_xxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 62);

    auto g_xxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 63);

    auto g_xxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 64);

    auto g_xxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 65);

    auto g_xxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 66);

    auto g_xxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 67);

    auto g_xxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 68);

    auto g_xxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 69);

    auto g_xxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 71);

    auto g_xxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 73);

    auto g_xxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 80);

    auto g_xxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 82);

    auto g_xxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 85);

    auto g_xxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 90);

    auto g_xxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 91);

    auto g_xxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 92);

    auto g_xxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 93);

    auto g_xxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 94);

    auto g_xxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 95);

    auto g_xxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 96);

    auto g_xxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 97);

    auto g_xxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 98);

    auto g_xxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 99);

    auto g_xyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 101);

    auto g_xyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 103);

    auto g_xyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 104);

    auto g_xyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 106);

    auto g_xyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 107);

    auto g_xyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 108);

    auto g_xyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 109);

    auto g_xyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 124);

    auto g_xyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 126);

    auto g_xyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 127);

    auto g_xyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 128);

    auto g_xyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 129);

    auto g_xzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 142);

    auto g_xzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 144);

    auto g_xzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 145);

    auto g_xzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 146);

    auto g_xzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 147);

    auto g_xzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 148);

    auto g_xzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 149);

    auto g_yyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 150);

    auto g_yyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 151);

    auto g_yyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 152);

    auto g_yyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 153);

    auto g_yyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 154);

    auto g_yyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 155);

    auto g_yyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 156);

    auto g_yyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 157);

    auto g_yyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 158);

    auto g_yyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 159);

    auto g_yyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 161);

    auto g_yyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 163);

    auto g_yyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 166);

    auto g_yyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 170);

    auto g_yyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 171);

    auto g_yyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 172);

    auto g_yyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 173);

    auto g_yyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 174);

    auto g_yyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 175);

    auto g_yyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 176);

    auto g_yyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 177);

    auto g_yyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 178);

    auto g_yyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 179);

    auto g_yyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 180);

    auto g_yyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 181);

    auto g_yyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 182);

    auto g_yyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 183);

    auto g_yyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 184);

    auto g_yyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 185);

    auto g_yyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 186);

    auto g_yyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 187);

    auto g_yyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 188);

    auto g_yyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 189);

    auto g_yzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 190);

    auto g_yzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 192);

    auto g_yzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 194);

    auto g_yzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 195);

    auto g_yzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 197);

    auto g_yzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 198);

    auto g_yzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 199);

    auto g_zzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 200);

    auto g_zzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 201);

    auto g_zzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 202);

    auto g_zzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 203);

    auto g_zzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 204);

    auto g_zzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 205);

    auto g_zzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 206);

    auto g_zzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 207);

    auto g_zzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 208);

    auto g_zzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 209);

    /// Set up components of auxilary buffer : HSF

    auto g_xxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_hsf);

    auto g_xxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 1);

    auto g_xxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 2);

    auto g_xxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 3);

    auto g_xxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 4);

    auto g_xxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 5);

    auto g_xxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 6);

    auto g_xxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 7);

    auto g_xxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 8);

    auto g_xxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 9);

    auto g_xxxxy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 10);

    auto g_xxxxy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 12);

    auto g_xxxxy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 15);

    auto g_xxxxz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 20);

    auto g_xxxxz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 21);

    auto g_xxxxz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 23);

    auto g_xxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 30);

    auto g_xxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 31);

    auto g_xxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 32);

    auto g_xxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 33);

    auto g_xxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 34);

    auto g_xxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 35);

    auto g_xxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 36);

    auto g_xxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 37);

    auto g_xxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 38);

    auto g_xxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 39);

    auto g_xxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 50);

    auto g_xxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 51);

    auto g_xxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 52);

    auto g_xxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 53);

    auto g_xxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 54);

    auto g_xxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 55);

    auto g_xxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 56);

    auto g_xxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 57);

    auto g_xxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 58);

    auto g_xxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 59);

    auto g_xxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 60);

    auto g_xxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 61);

    auto g_xxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 62);

    auto g_xxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 63);

    auto g_xxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 64);

    auto g_xxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 65);

    auto g_xxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 66);

    auto g_xxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 67);

    auto g_xxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 68);

    auto g_xxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 69);

    auto g_xxyyz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 71);

    auto g_xxyyz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 73);

    auto g_xxyzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 80);

    auto g_xxyzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 82);

    auto g_xxyzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 85);

    auto g_xxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 90);

    auto g_xxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 91);

    auto g_xxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 92);

    auto g_xxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 93);

    auto g_xxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 94);

    auto g_xxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 95);

    auto g_xxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 96);

    auto g_xxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 97);

    auto g_xxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 98);

    auto g_xxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 99);

    auto g_xyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 101);

    auto g_xyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 103);

    auto g_xyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 104);

    auto g_xyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 106);

    auto g_xyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 107);

    auto g_xyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 108);

    auto g_xyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 109);

    auto g_xyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 124);

    auto g_xyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 126);

    auto g_xyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 127);

    auto g_xyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 128);

    auto g_xyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 129);

    auto g_xzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 142);

    auto g_xzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 144);

    auto g_xzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 145);

    auto g_xzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 146);

    auto g_xzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 147);

    auto g_xzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 148);

    auto g_xzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 149);

    auto g_yyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 150);

    auto g_yyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 151);

    auto g_yyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 152);

    auto g_yyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 153);

    auto g_yyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 154);

    auto g_yyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 155);

    auto g_yyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 156);

    auto g_yyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 157);

    auto g_yyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 158);

    auto g_yyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 159);

    auto g_yyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 161);

    auto g_yyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 163);

    auto g_yyyyz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 166);

    auto g_yyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 170);

    auto g_yyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 171);

    auto g_yyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 172);

    auto g_yyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 173);

    auto g_yyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 174);

    auto g_yyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 175);

    auto g_yyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 176);

    auto g_yyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 177);

    auto g_yyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 178);

    auto g_yyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 179);

    auto g_yyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 180);

    auto g_yyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 181);

    auto g_yyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 182);

    auto g_yyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 183);

    auto g_yyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 184);

    auto g_yyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 185);

    auto g_yyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 186);

    auto g_yyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 187);

    auto g_yyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 188);

    auto g_yyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 189);

    auto g_yzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 190);

    auto g_yzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 192);

    auto g_yzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 194);

    auto g_yzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 195);

    auto g_yzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 197);

    auto g_yzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 198);

    auto g_yzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 199);

    auto g_zzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 200);

    auto g_zzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 201);

    auto g_zzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 202);

    auto g_zzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 203);

    auto g_zzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 204);

    auto g_zzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 205);

    auto g_zzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 206);

    auto g_zzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 207);

    auto g_zzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 208);

    auto g_zzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 209);

    /// Set up components of auxilary buffer : ISD

    auto g_xxxxxx_0_xx_1 = pbuffer.data(idx_eri_1_isd);

    auto g_xxxxxx_0_xy_1 = pbuffer.data(idx_eri_1_isd + 1);

    auto g_xxxxxx_0_xz_1 = pbuffer.data(idx_eri_1_isd + 2);

    auto g_xxxxxx_0_yy_1 = pbuffer.data(idx_eri_1_isd + 3);

    auto g_xxxxxx_0_yz_1 = pbuffer.data(idx_eri_1_isd + 4);

    auto g_xxxxxx_0_zz_1 = pbuffer.data(idx_eri_1_isd + 5);

    auto g_xxxxxz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 14);

    auto g_xxxxxz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 16);

    auto g_xxxxxz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 17);

    auto g_xxxxyy_0_xx_1 = pbuffer.data(idx_eri_1_isd + 18);

    auto g_xxxxyy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 19);

    auto g_xxxxyy_0_xz_1 = pbuffer.data(idx_eri_1_isd + 20);

    auto g_xxxxyy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 21);

    auto g_xxxxyy_0_yz_1 = pbuffer.data(idx_eri_1_isd + 22);

    auto g_xxxxyy_0_zz_1 = pbuffer.data(idx_eri_1_isd + 23);

    auto g_xxxxzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 30);

    auto g_xxxxzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 31);

    auto g_xxxxzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 32);

    auto g_xxxxzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 33);

    auto g_xxxxzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 34);

    auto g_xxxxzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 35);

    auto g_xxxyyy_0_xx_1 = pbuffer.data(idx_eri_1_isd + 36);

    auto g_xxxyyy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 37);

    auto g_xxxyyy_0_xz_1 = pbuffer.data(idx_eri_1_isd + 38);

    auto g_xxxyyy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 39);

    auto g_xxxyyy_0_yz_1 = pbuffer.data(idx_eri_1_isd + 40);

    auto g_xxxyyy_0_zz_1 = pbuffer.data(idx_eri_1_isd + 41);

    auto g_xxxzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 54);

    auto g_xxxzzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 55);

    auto g_xxxzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 56);

    auto g_xxxzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 57);

    auto g_xxxzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 58);

    auto g_xxxzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 59);

    auto g_xxyyyy_0_xx_1 = pbuffer.data(idx_eri_1_isd + 60);

    auto g_xxyyyy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 61);

    auto g_xxyyyy_0_xz_1 = pbuffer.data(idx_eri_1_isd + 62);

    auto g_xxyyyy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 63);

    auto g_xxyyyy_0_yz_1 = pbuffer.data(idx_eri_1_isd + 64);

    auto g_xxyyyy_0_zz_1 = pbuffer.data(idx_eri_1_isd + 65);

    auto g_xxyyzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 76);

    auto g_xxzzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 84);

    auto g_xxzzzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 85);

    auto g_xxzzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 86);

    auto g_xxzzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 87);

    auto g_xxzzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 88);

    auto g_xxzzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 89);

    auto g_xyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 91);

    auto g_xyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 93);

    auto g_xyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_isd + 94);

    auto g_xyyyzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 106);

    auto g_xyyzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 112);

    auto g_xzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 122);

    auto g_xzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 124);

    auto g_xzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 125);

    auto g_yyyyyy_0_xx_1 = pbuffer.data(idx_eri_1_isd + 126);

    auto g_yyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 127);

    auto g_yyyyyy_0_xz_1 = pbuffer.data(idx_eri_1_isd + 128);

    auto g_yyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 129);

    auto g_yyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_isd + 130);

    auto g_yyyyyy_0_zz_1 = pbuffer.data(idx_eri_1_isd + 131);

    auto g_yyyyyz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 134);

    auto g_yyyyyz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 136);

    auto g_yyyyyz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 137);

    auto g_yyyyzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 138);

    auto g_yyyyzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 139);

    auto g_yyyyzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 140);

    auto g_yyyyzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 141);

    auto g_yyyyzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 142);

    auto g_yyyyzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 143);

    auto g_yyyzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 144);

    auto g_yyyzzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 145);

    auto g_yyyzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 146);

    auto g_yyyzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 147);

    auto g_yyyzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 148);

    auto g_yyyzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 149);

    auto g_yyzzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 150);

    auto g_yyzzzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 151);

    auto g_yyzzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 152);

    auto g_yyzzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 153);

    auto g_yyzzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 154);

    auto g_yyzzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 155);

    auto g_yzzzzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 157);

    auto g_yzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 158);

    auto g_yzzzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 159);

    auto g_yzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 160);

    auto g_yzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 161);

    auto g_zzzzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 162);

    auto g_zzzzzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 163);

    auto g_zzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 164);

    auto g_zzzzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 165);

    auto g_zzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 166);

    auto g_zzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 167);

    /// Set up components of auxilary buffer : ISF

    auto g_xxxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_isf);

    auto g_xxxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 1);

    auto g_xxxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 2);

    auto g_xxxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 3);

    auto g_xxxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 4);

    auto g_xxxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 5);

    auto g_xxxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 6);

    auto g_xxxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 7);

    auto g_xxxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 8);

    auto g_xxxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 9);

    auto g_xxxxxy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 10);

    auto g_xxxxxy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 11);

    auto g_xxxxxy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 12);

    auto g_xxxxxy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 13);

    auto g_xxxxxy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 15);

    auto g_xxxxxy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 16);

    auto g_xxxxxz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 20);

    auto g_xxxxxz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 21);

    auto g_xxxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 22);

    auto g_xxxxxz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 23);

    auto g_xxxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 24);

    auto g_xxxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 25);

    auto g_xxxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 27);

    auto g_xxxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 28);

    auto g_xxxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 29);

    auto g_xxxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 30);

    auto g_xxxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 31);

    auto g_xxxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 32);

    auto g_xxxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 33);

    auto g_xxxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 34);

    auto g_xxxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 35);

    auto g_xxxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 36);

    auto g_xxxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 37);

    auto g_xxxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 38);

    auto g_xxxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 39);

    auto g_xxxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 50);

    auto g_xxxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 51);

    auto g_xxxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 52);

    auto g_xxxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 53);

    auto g_xxxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 54);

    auto g_xxxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 55);

    auto g_xxxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 56);

    auto g_xxxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 57);

    auto g_xxxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 58);

    auto g_xxxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 59);

    auto g_xxxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 60);

    auto g_xxxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 61);

    auto g_xxxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 62);

    auto g_xxxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 63);

    auto g_xxxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 64);

    auto g_xxxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 65);

    auto g_xxxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 66);

    auto g_xxxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 67);

    auto g_xxxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 68);

    auto g_xxxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 69);

    auto g_xxxyyz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 71);

    auto g_xxxyyz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 73);

    auto g_xxxyzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 80);

    auto g_xxxyzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 82);

    auto g_xxxyzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 85);

    auto g_xxxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 90);

    auto g_xxxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 91);

    auto g_xxxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 92);

    auto g_xxxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 93);

    auto g_xxxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 94);

    auto g_xxxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 95);

    auto g_xxxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 96);

    auto g_xxxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 97);

    auto g_xxxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 98);

    auto g_xxxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 99);

    auto g_xxyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 100);

    auto g_xxyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 101);

    auto g_xxyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 102);

    auto g_xxyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 103);

    auto g_xxyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 104);

    auto g_xxyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 105);

    auto g_xxyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 106);

    auto g_xxyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 107);

    auto g_xxyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 108);

    auto g_xxyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 109);

    auto g_xxyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 111);

    auto g_xxyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 113);

    auto g_xxyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 120);

    auto g_xxyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 121);

    auto g_xxyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 122);

    auto g_xxyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 123);

    auto g_xxyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 124);

    auto g_xxyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 125);

    auto g_xxyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 126);

    auto g_xxyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 127);

    auto g_xxyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 128);

    auto g_xxyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 129);

    auto g_xxyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 130);

    auto g_xxyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 132);

    auto g_xxyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 135);

    auto g_xxzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 140);

    auto g_xxzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 141);

    auto g_xxzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 142);

    auto g_xxzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 143);

    auto g_xxzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 144);

    auto g_xxzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 145);

    auto g_xxzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 146);

    auto g_xxzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 147);

    auto g_xxzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 148);

    auto g_xxzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 149);

    auto g_xyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 150);

    auto g_xyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 151);

    auto g_xyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 153);

    auto g_xyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 154);

    auto g_xyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 156);

    auto g_xyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 157);

    auto g_xyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 158);

    auto g_xyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 159);

    auto g_xyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 174);

    auto g_xyyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 176);

    auto g_xyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 177);

    auto g_xyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 178);

    auto g_xyyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 179);

    auto g_xyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 184);

    auto g_xyyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 186);

    auto g_xyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 187);

    auto g_xyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 188);

    auto g_xyyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 189);

    auto g_xzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 200);

    auto g_xzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 202);

    auto g_xzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 204);

    auto g_xzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 205);

    auto g_xzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 206);

    auto g_xzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 207);

    auto g_xzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 208);

    auto g_xzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 209);

    auto g_yyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 210);

    auto g_yyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 211);

    auto g_yyyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 212);

    auto g_yyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 213);

    auto g_yyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 214);

    auto g_yyyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 215);

    auto g_yyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 216);

    auto g_yyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 217);

    auto g_yyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 218);

    auto g_yyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 219);

    auto g_yyyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 221);

    auto g_yyyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 222);

    auto g_yyyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 223);

    auto g_yyyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 224);

    auto g_yyyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 225);

    auto g_yyyyyz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 226);

    auto g_yyyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 227);

    auto g_yyyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 228);

    auto g_yyyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 229);

    auto g_yyyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 230);

    auto g_yyyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 231);

    auto g_yyyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 232);

    auto g_yyyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 233);

    auto g_yyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 234);

    auto g_yyyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 235);

    auto g_yyyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 236);

    auto g_yyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 237);

    auto g_yyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 238);

    auto g_yyyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 239);

    auto g_yyyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 240);

    auto g_yyyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 241);

    auto g_yyyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 242);

    auto g_yyyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 243);

    auto g_yyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 244);

    auto g_yyyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 245);

    auto g_yyyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 246);

    auto g_yyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 247);

    auto g_yyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 248);

    auto g_yyyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 249);

    auto g_yyzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 250);

    auto g_yyzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 251);

    auto g_yyzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 252);

    auto g_yyzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 253);

    auto g_yyzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 254);

    auto g_yyzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 255);

    auto g_yyzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 256);

    auto g_yyzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 257);

    auto g_yyzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 258);

    auto g_yyzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 259);

    auto g_yzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 260);

    auto g_yzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 261);

    auto g_yzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 262);

    auto g_yzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 263);

    auto g_yzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 264);

    auto g_yzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 265);

    auto g_yzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 266);

    auto g_yzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 267);

    auto g_yzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 268);

    auto g_yzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 269);

    auto g_zzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 270);

    auto g_zzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 271);

    auto g_zzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 272);

    auto g_zzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 273);

    auto g_zzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 274);

    auto g_zzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 275);

    auto g_zzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 276);

    auto g_zzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 277);

    auto g_zzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 278);

    auto g_zzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 279);

    /// Set up 0-10 components of targeted buffer : KSF

    auto g_xxxxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_ksf);

    auto g_xxxxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 1);

    auto g_xxxxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 2);

    auto g_xxxxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 3);

    auto g_xxxxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 4);

    auto g_xxxxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 5);

    auto g_xxxxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 6);

    auto g_xxxxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 7);

    auto g_xxxxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 8);

    auto g_xxxxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 9);

    #pragma omp simd aligned(g_xxxxx_0_xxx_0, g_xxxxx_0_xxx_1, g_xxxxx_0_xxy_0, g_xxxxx_0_xxy_1, g_xxxxx_0_xxz_0, g_xxxxx_0_xxz_1, g_xxxxx_0_xyy_0, g_xxxxx_0_xyy_1, g_xxxxx_0_xyz_0, g_xxxxx_0_xyz_1, g_xxxxx_0_xzz_0, g_xxxxx_0_xzz_1, g_xxxxx_0_yyy_0, g_xxxxx_0_yyy_1, g_xxxxx_0_yyz_0, g_xxxxx_0_yyz_1, g_xxxxx_0_yzz_0, g_xxxxx_0_yzz_1, g_xxxxx_0_zzz_0, g_xxxxx_0_zzz_1, g_xxxxxx_0_xx_1, g_xxxxxx_0_xxx_1, g_xxxxxx_0_xxy_1, g_xxxxxx_0_xxz_1, g_xxxxxx_0_xy_1, g_xxxxxx_0_xyy_1, g_xxxxxx_0_xyz_1, g_xxxxxx_0_xz_1, g_xxxxxx_0_xzz_1, g_xxxxxx_0_yy_1, g_xxxxxx_0_yyy_1, g_xxxxxx_0_yyz_1, g_xxxxxx_0_yz_1, g_xxxxxx_0_yzz_1, g_xxxxxx_0_zz_1, g_xxxxxx_0_zzz_1, g_xxxxxxx_0_xxx_0, g_xxxxxxx_0_xxy_0, g_xxxxxxx_0_xxz_0, g_xxxxxxx_0_xyy_0, g_xxxxxxx_0_xyz_0, g_xxxxxxx_0_xzz_0, g_xxxxxxx_0_yyy_0, g_xxxxxxx_0_yyz_0, g_xxxxxxx_0_yzz_0, g_xxxxxxx_0_zzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_xxx_0[i] = 6.0 * g_xxxxx_0_xxx_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxx_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xx_1[i] * fi_acd_0 + g_xxxxxx_0_xxx_1[i] * wa_x[i];

        g_xxxxxxx_0_xxy_0[i] = 6.0 * g_xxxxx_0_xxy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xy_1[i] * fi_acd_0 + g_xxxxxx_0_xxy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxz_0[i] = 6.0 * g_xxxxx_0_xxz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xz_1[i] * fi_acd_0 + g_xxxxxx_0_xxz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyy_0[i] = 6.0 * g_xxxxx_0_xyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyy_1[i] * fz_be_0 + g_xxxxxx_0_yy_1[i] * fi_acd_0 + g_xxxxxx_0_xyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xyz_0[i] = 6.0 * g_xxxxx_0_xyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyz_1[i] * fz_be_0 + g_xxxxxx_0_yz_1[i] * fi_acd_0 + g_xxxxxx_0_xyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xzz_0[i] = 6.0 * g_xxxxx_0_xzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xzz_1[i] * fz_be_0 + g_xxxxxx_0_zz_1[i] * fi_acd_0 + g_xxxxxx_0_xzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyy_0[i] = 6.0 * g_xxxxx_0_yyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyy_1[i] * fz_be_0 + g_xxxxxx_0_yyy_1[i] * wa_x[i];

        g_xxxxxxx_0_yyz_0[i] = 6.0 * g_xxxxx_0_yyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyz_1[i] * fz_be_0 + g_xxxxxx_0_yyz_1[i] * wa_x[i];

        g_xxxxxxx_0_yzz_0[i] = 6.0 * g_xxxxx_0_yzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yzz_1[i] * fz_be_0 + g_xxxxxx_0_yzz_1[i] * wa_x[i];

        g_xxxxxxx_0_zzz_0[i] = 6.0 * g_xxxxx_0_zzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_zzz_1[i] * fz_be_0 + g_xxxxxx_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 10-20 components of targeted buffer : KSF

    auto g_xxxxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 10);

    auto g_xxxxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 11);

    auto g_xxxxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 12);

    auto g_xxxxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 13);

    auto g_xxxxxxy_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 14);

    auto g_xxxxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 15);

    auto g_xxxxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 16);

    auto g_xxxxxxy_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 17);

    auto g_xxxxxxy_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 18);

    auto g_xxxxxxy_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 19);

    #pragma omp simd aligned(g_xxxxxx_0_xx_1, g_xxxxxx_0_xxx_1, g_xxxxxx_0_xxy_1, g_xxxxxx_0_xxz_1, g_xxxxxx_0_xy_1, g_xxxxxx_0_xyy_1, g_xxxxxx_0_xyz_1, g_xxxxxx_0_xz_1, g_xxxxxx_0_xzz_1, g_xxxxxx_0_yy_1, g_xxxxxx_0_yyy_1, g_xxxxxx_0_yyz_1, g_xxxxxx_0_yz_1, g_xxxxxx_0_yzz_1, g_xxxxxx_0_zz_1, g_xxxxxx_0_zzz_1, g_xxxxxxy_0_xxx_0, g_xxxxxxy_0_xxy_0, g_xxxxxxy_0_xxz_0, g_xxxxxxy_0_xyy_0, g_xxxxxxy_0_xyz_0, g_xxxxxxy_0_xzz_0, g_xxxxxxy_0_yyy_0, g_xxxxxxy_0_yyz_0, g_xxxxxxy_0_yzz_0, g_xxxxxxy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_xxx_0[i] = g_xxxxxx_0_xxx_1[i] * wa_y[i];

        g_xxxxxxy_0_xxy_0[i] = g_xxxxxx_0_xx_1[i] * fi_acd_0 + g_xxxxxx_0_xxy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxz_0[i] = g_xxxxxx_0_xxz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyy_0[i] = 2.0 * g_xxxxxx_0_xy_1[i] * fi_acd_0 + g_xxxxxx_0_xyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xyz_0[i] = g_xxxxxx_0_xz_1[i] * fi_acd_0 + g_xxxxxx_0_xyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xzz_0[i] = g_xxxxxx_0_xzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyy_0[i] = 3.0 * g_xxxxxx_0_yy_1[i] * fi_acd_0 + g_xxxxxx_0_yyy_1[i] * wa_y[i];

        g_xxxxxxy_0_yyz_0[i] = 2.0 * g_xxxxxx_0_yz_1[i] * fi_acd_0 + g_xxxxxx_0_yyz_1[i] * wa_y[i];

        g_xxxxxxy_0_yzz_0[i] = g_xxxxxx_0_zz_1[i] * fi_acd_0 + g_xxxxxx_0_yzz_1[i] * wa_y[i];

        g_xxxxxxy_0_zzz_0[i] = g_xxxxxx_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 20-30 components of targeted buffer : KSF

    auto g_xxxxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 20);

    auto g_xxxxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 21);

    auto g_xxxxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 22);

    auto g_xxxxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 23);

    auto g_xxxxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 24);

    auto g_xxxxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 25);

    auto g_xxxxxxz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 26);

    auto g_xxxxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 27);

    auto g_xxxxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 28);

    auto g_xxxxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 29);

    #pragma omp simd aligned(g_xxxxxx_0_xx_1, g_xxxxxx_0_xxx_1, g_xxxxxx_0_xxy_1, g_xxxxxx_0_xxz_1, g_xxxxxx_0_xy_1, g_xxxxxx_0_xyy_1, g_xxxxxx_0_xyz_1, g_xxxxxx_0_xz_1, g_xxxxxx_0_xzz_1, g_xxxxxx_0_yy_1, g_xxxxxx_0_yyy_1, g_xxxxxx_0_yyz_1, g_xxxxxx_0_yz_1, g_xxxxxx_0_yzz_1, g_xxxxxx_0_zz_1, g_xxxxxx_0_zzz_1, g_xxxxxxz_0_xxx_0, g_xxxxxxz_0_xxy_0, g_xxxxxxz_0_xxz_0, g_xxxxxxz_0_xyy_0, g_xxxxxxz_0_xyz_0, g_xxxxxxz_0_xzz_0, g_xxxxxxz_0_yyy_0, g_xxxxxxz_0_yyz_0, g_xxxxxxz_0_yzz_0, g_xxxxxxz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_xxx_0[i] = g_xxxxxx_0_xxx_1[i] * wa_z[i];

        g_xxxxxxz_0_xxy_0[i] = g_xxxxxx_0_xxy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxz_0[i] = g_xxxxxx_0_xx_1[i] * fi_acd_0 + g_xxxxxx_0_xxz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyy_0[i] = g_xxxxxx_0_xyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xyz_0[i] = g_xxxxxx_0_xy_1[i] * fi_acd_0 + g_xxxxxx_0_xyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xzz_0[i] = 2.0 * g_xxxxxx_0_xz_1[i] * fi_acd_0 + g_xxxxxx_0_xzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyy_0[i] = g_xxxxxx_0_yyy_1[i] * wa_z[i];

        g_xxxxxxz_0_yyz_0[i] = g_xxxxxx_0_yy_1[i] * fi_acd_0 + g_xxxxxx_0_yyz_1[i] * wa_z[i];

        g_xxxxxxz_0_yzz_0[i] = 2.0 * g_xxxxxx_0_yz_1[i] * fi_acd_0 + g_xxxxxx_0_yzz_1[i] * wa_z[i];

        g_xxxxxxz_0_zzz_0[i] = 3.0 * g_xxxxxx_0_zz_1[i] * fi_acd_0 + g_xxxxxx_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 30-40 components of targeted buffer : KSF

    auto g_xxxxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 30);

    auto g_xxxxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 31);

    auto g_xxxxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 32);

    auto g_xxxxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 33);

    auto g_xxxxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 34);

    auto g_xxxxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 35);

    auto g_xxxxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 36);

    auto g_xxxxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 37);

    auto g_xxxxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 38);

    auto g_xxxxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 39);

    #pragma omp simd aligned(g_xxxxx_0_xxx_0, g_xxxxx_0_xxx_1, g_xxxxx_0_xxz_0, g_xxxxx_0_xxz_1, g_xxxxx_0_xzz_0, g_xxxxx_0_xzz_1, g_xxxxxy_0_xxx_1, g_xxxxxy_0_xxz_1, g_xxxxxy_0_xzz_1, g_xxxxxyy_0_xxx_0, g_xxxxxyy_0_xxy_0, g_xxxxxyy_0_xxz_0, g_xxxxxyy_0_xyy_0, g_xxxxxyy_0_xyz_0, g_xxxxxyy_0_xzz_0, g_xxxxxyy_0_yyy_0, g_xxxxxyy_0_yyz_0, g_xxxxxyy_0_yzz_0, g_xxxxxyy_0_zzz_0, g_xxxxyy_0_xxy_1, g_xxxxyy_0_xy_1, g_xxxxyy_0_xyy_1, g_xxxxyy_0_xyz_1, g_xxxxyy_0_yy_1, g_xxxxyy_0_yyy_1, g_xxxxyy_0_yyz_1, g_xxxxyy_0_yz_1, g_xxxxyy_0_yzz_1, g_xxxxyy_0_zzz_1, g_xxxyy_0_xxy_0, g_xxxyy_0_xxy_1, g_xxxyy_0_xyy_0, g_xxxyy_0_xyy_1, g_xxxyy_0_xyz_0, g_xxxyy_0_xyz_1, g_xxxyy_0_yyy_0, g_xxxyy_0_yyy_1, g_xxxyy_0_yyz_0, g_xxxyy_0_yyz_1, g_xxxyy_0_yzz_0, g_xxxyy_0_yzz_1, g_xxxyy_0_zzz_0, g_xxxyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxyy_0_xxx_0[i] = g_xxxxx_0_xxx_0[i] * fbe_0 - g_xxxxx_0_xxx_1[i] * fz_be_0 + g_xxxxxy_0_xxx_1[i] * wa_y[i];

        g_xxxxxyy_0_xxy_0[i] = 4.0 * g_xxxyy_0_xxy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xy_1[i] * fi_acd_0 + g_xxxxyy_0_xxy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxz_0[i] = g_xxxxx_0_xxz_0[i] * fbe_0 - g_xxxxx_0_xxz_1[i] * fz_be_0 + g_xxxxxy_0_xxz_1[i] * wa_y[i];

        g_xxxxxyy_0_xyy_0[i] = 4.0 * g_xxxyy_0_xyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyy_1[i] * fz_be_0 + g_xxxxyy_0_yy_1[i] * fi_acd_0 + g_xxxxyy_0_xyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xyz_0[i] = 4.0 * g_xxxyy_0_xyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyz_1[i] * fz_be_0 + g_xxxxyy_0_yz_1[i] * fi_acd_0 + g_xxxxyy_0_xyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xzz_0[i] = g_xxxxx_0_xzz_0[i] * fbe_0 - g_xxxxx_0_xzz_1[i] * fz_be_0 + g_xxxxxy_0_xzz_1[i] * wa_y[i];

        g_xxxxxyy_0_yyy_0[i] = 4.0 * g_xxxyy_0_yyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyy_1[i] * fz_be_0 + g_xxxxyy_0_yyy_1[i] * wa_x[i];

        g_xxxxxyy_0_yyz_0[i] = 4.0 * g_xxxyy_0_yyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyz_1[i] * fz_be_0 + g_xxxxyy_0_yyz_1[i] * wa_x[i];

        g_xxxxxyy_0_yzz_0[i] = 4.0 * g_xxxyy_0_yzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yzz_1[i] * fz_be_0 + g_xxxxyy_0_yzz_1[i] * wa_x[i];

        g_xxxxxyy_0_zzz_0[i] = 4.0 * g_xxxyy_0_zzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_zzz_1[i] * fz_be_0 + g_xxxxyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 40-50 components of targeted buffer : KSF

    auto g_xxxxxyz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 40);

    auto g_xxxxxyz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 41);

    auto g_xxxxxyz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 42);

    auto g_xxxxxyz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 43);

    auto g_xxxxxyz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 44);

    auto g_xxxxxyz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 45);

    auto g_xxxxxyz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 46);

    auto g_xxxxxyz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 47);

    auto g_xxxxxyz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 48);

    auto g_xxxxxyz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 49);

    #pragma omp simd aligned(g_xxxxxy_0_xxy_1, g_xxxxxy_0_xyy_1, g_xxxxxy_0_yyy_1, g_xxxxxyz_0_xxx_0, g_xxxxxyz_0_xxy_0, g_xxxxxyz_0_xxz_0, g_xxxxxyz_0_xyy_0, g_xxxxxyz_0_xyz_0, g_xxxxxyz_0_xzz_0, g_xxxxxyz_0_yyy_0, g_xxxxxyz_0_yyz_0, g_xxxxxyz_0_yzz_0, g_xxxxxyz_0_zzz_0, g_xxxxxz_0_xxx_1, g_xxxxxz_0_xxz_1, g_xxxxxz_0_xyz_1, g_xxxxxz_0_xz_1, g_xxxxxz_0_xzz_1, g_xxxxxz_0_yyz_1, g_xxxxxz_0_yz_1, g_xxxxxz_0_yzz_1, g_xxxxxz_0_zz_1, g_xxxxxz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyz_0_xxx_0[i] = g_xxxxxz_0_xxx_1[i] * wa_y[i];

        g_xxxxxyz_0_xxy_0[i] = g_xxxxxy_0_xxy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxz_0[i] = g_xxxxxz_0_xxz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyy_0[i] = g_xxxxxy_0_xyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xyz_0[i] = g_xxxxxz_0_xz_1[i] * fi_acd_0 + g_xxxxxz_0_xyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xzz_0[i] = g_xxxxxz_0_xzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyy_0[i] = g_xxxxxy_0_yyy_1[i] * wa_z[i];

        g_xxxxxyz_0_yyz_0[i] = 2.0 * g_xxxxxz_0_yz_1[i] * fi_acd_0 + g_xxxxxz_0_yyz_1[i] * wa_y[i];

        g_xxxxxyz_0_yzz_0[i] = g_xxxxxz_0_zz_1[i] * fi_acd_0 + g_xxxxxz_0_yzz_1[i] * wa_y[i];

        g_xxxxxyz_0_zzz_0[i] = g_xxxxxz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 50-60 components of targeted buffer : KSF

    auto g_xxxxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 50);

    auto g_xxxxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 51);

    auto g_xxxxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 52);

    auto g_xxxxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 53);

    auto g_xxxxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 54);

    auto g_xxxxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 55);

    auto g_xxxxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 56);

    auto g_xxxxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 57);

    auto g_xxxxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 58);

    auto g_xxxxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 59);

    #pragma omp simd aligned(g_xxxxx_0_xxx_0, g_xxxxx_0_xxx_1, g_xxxxx_0_xxy_0, g_xxxxx_0_xxy_1, g_xxxxx_0_xyy_0, g_xxxxx_0_xyy_1, g_xxxxxz_0_xxx_1, g_xxxxxz_0_xxy_1, g_xxxxxz_0_xyy_1, g_xxxxxzz_0_xxx_0, g_xxxxxzz_0_xxy_0, g_xxxxxzz_0_xxz_0, g_xxxxxzz_0_xyy_0, g_xxxxxzz_0_xyz_0, g_xxxxxzz_0_xzz_0, g_xxxxxzz_0_yyy_0, g_xxxxxzz_0_yyz_0, g_xxxxxzz_0_yzz_0, g_xxxxxzz_0_zzz_0, g_xxxxzz_0_xxz_1, g_xxxxzz_0_xyz_1, g_xxxxzz_0_xz_1, g_xxxxzz_0_xzz_1, g_xxxxzz_0_yyy_1, g_xxxxzz_0_yyz_1, g_xxxxzz_0_yz_1, g_xxxxzz_0_yzz_1, g_xxxxzz_0_zz_1, g_xxxxzz_0_zzz_1, g_xxxzz_0_xxz_0, g_xxxzz_0_xxz_1, g_xxxzz_0_xyz_0, g_xxxzz_0_xyz_1, g_xxxzz_0_xzz_0, g_xxxzz_0_xzz_1, g_xxxzz_0_yyy_0, g_xxxzz_0_yyy_1, g_xxxzz_0_yyz_0, g_xxxzz_0_yyz_1, g_xxxzz_0_yzz_0, g_xxxzz_0_yzz_1, g_xxxzz_0_zzz_0, g_xxxzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxzz_0_xxx_0[i] = g_xxxxx_0_xxx_0[i] * fbe_0 - g_xxxxx_0_xxx_1[i] * fz_be_0 + g_xxxxxz_0_xxx_1[i] * wa_z[i];

        g_xxxxxzz_0_xxy_0[i] = g_xxxxx_0_xxy_0[i] * fbe_0 - g_xxxxx_0_xxy_1[i] * fz_be_0 + g_xxxxxz_0_xxy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxz_0[i] = 4.0 * g_xxxzz_0_xxz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xz_1[i] * fi_acd_0 + g_xxxxzz_0_xxz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyy_0[i] = g_xxxxx_0_xyy_0[i] * fbe_0 - g_xxxxx_0_xyy_1[i] * fz_be_0 + g_xxxxxz_0_xyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xyz_0[i] = 4.0 * g_xxxzz_0_xyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyz_1[i] * fz_be_0 + g_xxxxzz_0_yz_1[i] * fi_acd_0 + g_xxxxzz_0_xyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xzz_0[i] = 4.0 * g_xxxzz_0_xzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xzz_1[i] * fz_be_0 + g_xxxxzz_0_zz_1[i] * fi_acd_0 + g_xxxxzz_0_xzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyy_0[i] = 4.0 * g_xxxzz_0_yyy_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyy_1[i] * fz_be_0 + g_xxxxzz_0_yyy_1[i] * wa_x[i];

        g_xxxxxzz_0_yyz_0[i] = 4.0 * g_xxxzz_0_yyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyz_1[i] * fz_be_0 + g_xxxxzz_0_yyz_1[i] * wa_x[i];

        g_xxxxxzz_0_yzz_0[i] = 4.0 * g_xxxzz_0_yzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yzz_1[i] * fz_be_0 + g_xxxxzz_0_yzz_1[i] * wa_x[i];

        g_xxxxxzz_0_zzz_0[i] = 4.0 * g_xxxzz_0_zzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_zzz_1[i] * fz_be_0 + g_xxxxzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 60-70 components of targeted buffer : KSF

    auto g_xxxxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 60);

    auto g_xxxxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 61);

    auto g_xxxxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 62);

    auto g_xxxxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 63);

    auto g_xxxxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 64);

    auto g_xxxxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 65);

    auto g_xxxxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 66);

    auto g_xxxxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 67);

    auto g_xxxxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 68);

    auto g_xxxxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 69);

    #pragma omp simd aligned(g_xxxxy_0_xxx_0, g_xxxxy_0_xxx_1, g_xxxxy_0_xxz_0, g_xxxxy_0_xxz_1, g_xxxxy_0_xzz_0, g_xxxxy_0_xzz_1, g_xxxxyy_0_xxx_1, g_xxxxyy_0_xxz_1, g_xxxxyy_0_xzz_1, g_xxxxyyy_0_xxx_0, g_xxxxyyy_0_xxy_0, g_xxxxyyy_0_xxz_0, g_xxxxyyy_0_xyy_0, g_xxxxyyy_0_xyz_0, g_xxxxyyy_0_xzz_0, g_xxxxyyy_0_yyy_0, g_xxxxyyy_0_yyz_0, g_xxxxyyy_0_yzz_0, g_xxxxyyy_0_zzz_0, g_xxxyyy_0_xxy_1, g_xxxyyy_0_xy_1, g_xxxyyy_0_xyy_1, g_xxxyyy_0_xyz_1, g_xxxyyy_0_yy_1, g_xxxyyy_0_yyy_1, g_xxxyyy_0_yyz_1, g_xxxyyy_0_yz_1, g_xxxyyy_0_yzz_1, g_xxxyyy_0_zzz_1, g_xxyyy_0_xxy_0, g_xxyyy_0_xxy_1, g_xxyyy_0_xyy_0, g_xxyyy_0_xyy_1, g_xxyyy_0_xyz_0, g_xxyyy_0_xyz_1, g_xxyyy_0_yyy_0, g_xxyyy_0_yyy_1, g_xxyyy_0_yyz_0, g_xxyyy_0_yyz_1, g_xxyyy_0_yzz_0, g_xxyyy_0_yzz_1, g_xxyyy_0_zzz_0, g_xxyyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyyy_0_xxx_0[i] = 2.0 * g_xxxxy_0_xxx_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxx_1[i] * fz_be_0 + g_xxxxyy_0_xxx_1[i] * wa_y[i];

        g_xxxxyyy_0_xxy_0[i] = 3.0 * g_xxyyy_0_xxy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xy_1[i] * fi_acd_0 + g_xxxyyy_0_xxy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxz_0[i] = 2.0 * g_xxxxy_0_xxz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxz_1[i] * fz_be_0 + g_xxxxyy_0_xxz_1[i] * wa_y[i];

        g_xxxxyyy_0_xyy_0[i] = 3.0 * g_xxyyy_0_xyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyy_1[i] * fz_be_0 + g_xxxyyy_0_yy_1[i] * fi_acd_0 + g_xxxyyy_0_xyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xyz_0[i] = 3.0 * g_xxyyy_0_xyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyz_1[i] * fz_be_0 + g_xxxyyy_0_yz_1[i] * fi_acd_0 + g_xxxyyy_0_xyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xzz_0[i] = 2.0 * g_xxxxy_0_xzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xzz_1[i] * fz_be_0 + g_xxxxyy_0_xzz_1[i] * wa_y[i];

        g_xxxxyyy_0_yyy_0[i] = 3.0 * g_xxyyy_0_yyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyy_1[i] * fz_be_0 + g_xxxyyy_0_yyy_1[i] * wa_x[i];

        g_xxxxyyy_0_yyz_0[i] = 3.0 * g_xxyyy_0_yyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyz_1[i] * fz_be_0 + g_xxxyyy_0_yyz_1[i] * wa_x[i];

        g_xxxxyyy_0_yzz_0[i] = 3.0 * g_xxyyy_0_yzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yzz_1[i] * fz_be_0 + g_xxxyyy_0_yzz_1[i] * wa_x[i];

        g_xxxxyyy_0_zzz_0[i] = 3.0 * g_xxyyy_0_zzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_zzz_1[i] * fz_be_0 + g_xxxyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 70-80 components of targeted buffer : KSF

    auto g_xxxxyyz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 70);

    auto g_xxxxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 71);

    auto g_xxxxyyz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 72);

    auto g_xxxxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 73);

    auto g_xxxxyyz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 74);

    auto g_xxxxyyz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 75);

    auto g_xxxxyyz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 76);

    auto g_xxxxyyz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 77);

    auto g_xxxxyyz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 78);

    auto g_xxxxyyz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 79);

    #pragma omp simd aligned(g_xxxxyy_0_xx_1, g_xxxxyy_0_xxx_1, g_xxxxyy_0_xxy_1, g_xxxxyy_0_xxz_1, g_xxxxyy_0_xy_1, g_xxxxyy_0_xyy_1, g_xxxxyy_0_xyz_1, g_xxxxyy_0_xz_1, g_xxxxyy_0_xzz_1, g_xxxxyy_0_yy_1, g_xxxxyy_0_yyy_1, g_xxxxyy_0_yyz_1, g_xxxxyy_0_yz_1, g_xxxxyy_0_yzz_1, g_xxxxyy_0_zz_1, g_xxxxyy_0_zzz_1, g_xxxxyyz_0_xxx_0, g_xxxxyyz_0_xxy_0, g_xxxxyyz_0_xxz_0, g_xxxxyyz_0_xyy_0, g_xxxxyyz_0_xyz_0, g_xxxxyyz_0_xzz_0, g_xxxxyyz_0_yyy_0, g_xxxxyyz_0_yyz_0, g_xxxxyyz_0_yzz_0, g_xxxxyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_xxx_0[i] = g_xxxxyy_0_xxx_1[i] * wa_z[i];

        g_xxxxyyz_0_xxy_0[i] = g_xxxxyy_0_xxy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxz_0[i] = g_xxxxyy_0_xx_1[i] * fi_acd_0 + g_xxxxyy_0_xxz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyy_0[i] = g_xxxxyy_0_xyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xyz_0[i] = g_xxxxyy_0_xy_1[i] * fi_acd_0 + g_xxxxyy_0_xyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xzz_0[i] = 2.0 * g_xxxxyy_0_xz_1[i] * fi_acd_0 + g_xxxxyy_0_xzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyy_0[i] = g_xxxxyy_0_yyy_1[i] * wa_z[i];

        g_xxxxyyz_0_yyz_0[i] = g_xxxxyy_0_yy_1[i] * fi_acd_0 + g_xxxxyy_0_yyz_1[i] * wa_z[i];

        g_xxxxyyz_0_yzz_0[i] = 2.0 * g_xxxxyy_0_yz_1[i] * fi_acd_0 + g_xxxxyy_0_yzz_1[i] * wa_z[i];

        g_xxxxyyz_0_zzz_0[i] = 3.0 * g_xxxxyy_0_zz_1[i] * fi_acd_0 + g_xxxxyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 80-90 components of targeted buffer : KSF

    auto g_xxxxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 80);

    auto g_xxxxyzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 81);

    auto g_xxxxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 82);

    auto g_xxxxyzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 83);

    auto g_xxxxyzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 84);

    auto g_xxxxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 85);

    auto g_xxxxyzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 86);

    auto g_xxxxyzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 87);

    auto g_xxxxyzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 88);

    auto g_xxxxyzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 89);

    #pragma omp simd aligned(g_xxxxyzz_0_xxx_0, g_xxxxyzz_0_xxy_0, g_xxxxyzz_0_xxz_0, g_xxxxyzz_0_xyy_0, g_xxxxyzz_0_xyz_0, g_xxxxyzz_0_xzz_0, g_xxxxyzz_0_yyy_0, g_xxxxyzz_0_yyz_0, g_xxxxyzz_0_yzz_0, g_xxxxyzz_0_zzz_0, g_xxxxzz_0_xx_1, g_xxxxzz_0_xxx_1, g_xxxxzz_0_xxy_1, g_xxxxzz_0_xxz_1, g_xxxxzz_0_xy_1, g_xxxxzz_0_xyy_1, g_xxxxzz_0_xyz_1, g_xxxxzz_0_xz_1, g_xxxxzz_0_xzz_1, g_xxxxzz_0_yy_1, g_xxxxzz_0_yyy_1, g_xxxxzz_0_yyz_1, g_xxxxzz_0_yz_1, g_xxxxzz_0_yzz_1, g_xxxxzz_0_zz_1, g_xxxxzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_xxx_0[i] = g_xxxxzz_0_xxx_1[i] * wa_y[i];

        g_xxxxyzz_0_xxy_0[i] = g_xxxxzz_0_xx_1[i] * fi_acd_0 + g_xxxxzz_0_xxy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxz_0[i] = g_xxxxzz_0_xxz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyy_0[i] = 2.0 * g_xxxxzz_0_xy_1[i] * fi_acd_0 + g_xxxxzz_0_xyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xyz_0[i] = g_xxxxzz_0_xz_1[i] * fi_acd_0 + g_xxxxzz_0_xyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xzz_0[i] = g_xxxxzz_0_xzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyy_0[i] = 3.0 * g_xxxxzz_0_yy_1[i] * fi_acd_0 + g_xxxxzz_0_yyy_1[i] * wa_y[i];

        g_xxxxyzz_0_yyz_0[i] = 2.0 * g_xxxxzz_0_yz_1[i] * fi_acd_0 + g_xxxxzz_0_yyz_1[i] * wa_y[i];

        g_xxxxyzz_0_yzz_0[i] = g_xxxxzz_0_zz_1[i] * fi_acd_0 + g_xxxxzz_0_yzz_1[i] * wa_y[i];

        g_xxxxyzz_0_zzz_0[i] = g_xxxxzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 90-100 components of targeted buffer : KSF

    auto g_xxxxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 90);

    auto g_xxxxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 91);

    auto g_xxxxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 92);

    auto g_xxxxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 93);

    auto g_xxxxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 94);

    auto g_xxxxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 95);

    auto g_xxxxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 96);

    auto g_xxxxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 97);

    auto g_xxxxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 98);

    auto g_xxxxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 99);

    #pragma omp simd aligned(g_xxxxz_0_xxx_0, g_xxxxz_0_xxx_1, g_xxxxz_0_xxy_0, g_xxxxz_0_xxy_1, g_xxxxz_0_xyy_0, g_xxxxz_0_xyy_1, g_xxxxzz_0_xxx_1, g_xxxxzz_0_xxy_1, g_xxxxzz_0_xyy_1, g_xxxxzzz_0_xxx_0, g_xxxxzzz_0_xxy_0, g_xxxxzzz_0_xxz_0, g_xxxxzzz_0_xyy_0, g_xxxxzzz_0_xyz_0, g_xxxxzzz_0_xzz_0, g_xxxxzzz_0_yyy_0, g_xxxxzzz_0_yyz_0, g_xxxxzzz_0_yzz_0, g_xxxxzzz_0_zzz_0, g_xxxzzz_0_xxz_1, g_xxxzzz_0_xyz_1, g_xxxzzz_0_xz_1, g_xxxzzz_0_xzz_1, g_xxxzzz_0_yyy_1, g_xxxzzz_0_yyz_1, g_xxxzzz_0_yz_1, g_xxxzzz_0_yzz_1, g_xxxzzz_0_zz_1, g_xxxzzz_0_zzz_1, g_xxzzz_0_xxz_0, g_xxzzz_0_xxz_1, g_xxzzz_0_xyz_0, g_xxzzz_0_xyz_1, g_xxzzz_0_xzz_0, g_xxzzz_0_xzz_1, g_xxzzz_0_yyy_0, g_xxzzz_0_yyy_1, g_xxzzz_0_yyz_0, g_xxzzz_0_yyz_1, g_xxzzz_0_yzz_0, g_xxzzz_0_yzz_1, g_xxzzz_0_zzz_0, g_xxzzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzzz_0_xxx_0[i] = 2.0 * g_xxxxz_0_xxx_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxx_1[i] * fz_be_0 + g_xxxxzz_0_xxx_1[i] * wa_z[i];

        g_xxxxzzz_0_xxy_0[i] = 2.0 * g_xxxxz_0_xxy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxy_1[i] * fz_be_0 + g_xxxxzz_0_xxy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxz_0[i] = 3.0 * g_xxzzz_0_xxz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xz_1[i] * fi_acd_0 + g_xxxzzz_0_xxz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyy_0[i] = 2.0 * g_xxxxz_0_xyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xyy_1[i] * fz_be_0 + g_xxxxzz_0_xyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xyz_0[i] = 3.0 * g_xxzzz_0_xyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyz_1[i] * fz_be_0 + g_xxxzzz_0_yz_1[i] * fi_acd_0 + g_xxxzzz_0_xyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xzz_0[i] = 3.0 * g_xxzzz_0_xzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xzz_1[i] * fz_be_0 + g_xxxzzz_0_zz_1[i] * fi_acd_0 + g_xxxzzz_0_xzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyy_0[i] = 3.0 * g_xxzzz_0_yyy_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyy_1[i] * fz_be_0 + g_xxxzzz_0_yyy_1[i] * wa_x[i];

        g_xxxxzzz_0_yyz_0[i] = 3.0 * g_xxzzz_0_yyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyz_1[i] * fz_be_0 + g_xxxzzz_0_yyz_1[i] * wa_x[i];

        g_xxxxzzz_0_yzz_0[i] = 3.0 * g_xxzzz_0_yzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yzz_1[i] * fz_be_0 + g_xxxzzz_0_yzz_1[i] * wa_x[i];

        g_xxxxzzz_0_zzz_0[i] = 3.0 * g_xxzzz_0_zzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_zzz_1[i] * fz_be_0 + g_xxxzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 100-110 components of targeted buffer : KSF

    auto g_xxxyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 100);

    auto g_xxxyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 101);

    auto g_xxxyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 102);

    auto g_xxxyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 103);

    auto g_xxxyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 104);

    auto g_xxxyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 105);

    auto g_xxxyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 106);

    auto g_xxxyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 107);

    auto g_xxxyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 108);

    auto g_xxxyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 109);

    #pragma omp simd aligned(g_xxxyy_0_xxx_0, g_xxxyy_0_xxx_1, g_xxxyy_0_xxz_0, g_xxxyy_0_xxz_1, g_xxxyy_0_xzz_0, g_xxxyy_0_xzz_1, g_xxxyyy_0_xxx_1, g_xxxyyy_0_xxz_1, g_xxxyyy_0_xzz_1, g_xxxyyyy_0_xxx_0, g_xxxyyyy_0_xxy_0, g_xxxyyyy_0_xxz_0, g_xxxyyyy_0_xyy_0, g_xxxyyyy_0_xyz_0, g_xxxyyyy_0_xzz_0, g_xxxyyyy_0_yyy_0, g_xxxyyyy_0_yyz_0, g_xxxyyyy_0_yzz_0, g_xxxyyyy_0_zzz_0, g_xxyyyy_0_xxy_1, g_xxyyyy_0_xy_1, g_xxyyyy_0_xyy_1, g_xxyyyy_0_xyz_1, g_xxyyyy_0_yy_1, g_xxyyyy_0_yyy_1, g_xxyyyy_0_yyz_1, g_xxyyyy_0_yz_1, g_xxyyyy_0_yzz_1, g_xxyyyy_0_zzz_1, g_xyyyy_0_xxy_0, g_xyyyy_0_xxy_1, g_xyyyy_0_xyy_0, g_xyyyy_0_xyy_1, g_xyyyy_0_xyz_0, g_xyyyy_0_xyz_1, g_xyyyy_0_yyy_0, g_xyyyy_0_yyy_1, g_xyyyy_0_yyz_0, g_xyyyy_0_yyz_1, g_xyyyy_0_yzz_0, g_xyyyy_0_yzz_1, g_xyyyy_0_zzz_0, g_xyyyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyyy_0_xxx_0[i] = 3.0 * g_xxxyy_0_xxx_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxx_1[i] * fz_be_0 + g_xxxyyy_0_xxx_1[i] * wa_y[i];

        g_xxxyyyy_0_xxy_0[i] = 2.0 * g_xyyyy_0_xxy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xy_1[i] * fi_acd_0 + g_xxyyyy_0_xxy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxz_0[i] = 3.0 * g_xxxyy_0_xxz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxz_1[i] * fz_be_0 + g_xxxyyy_0_xxz_1[i] * wa_y[i];

        g_xxxyyyy_0_xyy_0[i] = 2.0 * g_xyyyy_0_xyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyy_1[i] * fz_be_0 + g_xxyyyy_0_yy_1[i] * fi_acd_0 + g_xxyyyy_0_xyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xyz_0[i] = 2.0 * g_xyyyy_0_xyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyz_1[i] * fz_be_0 + g_xxyyyy_0_yz_1[i] * fi_acd_0 + g_xxyyyy_0_xyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xzz_0[i] = 3.0 * g_xxxyy_0_xzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xzz_1[i] * fz_be_0 + g_xxxyyy_0_xzz_1[i] * wa_y[i];

        g_xxxyyyy_0_yyy_0[i] = 2.0 * g_xyyyy_0_yyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyy_1[i] * fz_be_0 + g_xxyyyy_0_yyy_1[i] * wa_x[i];

        g_xxxyyyy_0_yyz_0[i] = 2.0 * g_xyyyy_0_yyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyz_1[i] * fz_be_0 + g_xxyyyy_0_yyz_1[i] * wa_x[i];

        g_xxxyyyy_0_yzz_0[i] = 2.0 * g_xyyyy_0_yzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yzz_1[i] * fz_be_0 + g_xxyyyy_0_yzz_1[i] * wa_x[i];

        g_xxxyyyy_0_zzz_0[i] = 2.0 * g_xyyyy_0_zzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_zzz_1[i] * fz_be_0 + g_xxyyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 110-120 components of targeted buffer : KSF

    auto g_xxxyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 110);

    auto g_xxxyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 111);

    auto g_xxxyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 112);

    auto g_xxxyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 113);

    auto g_xxxyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 114);

    auto g_xxxyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 115);

    auto g_xxxyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 116);

    auto g_xxxyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 117);

    auto g_xxxyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 118);

    auto g_xxxyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 119);

    #pragma omp simd aligned(g_xxxyyy_0_xx_1, g_xxxyyy_0_xxx_1, g_xxxyyy_0_xxy_1, g_xxxyyy_0_xxz_1, g_xxxyyy_0_xy_1, g_xxxyyy_0_xyy_1, g_xxxyyy_0_xyz_1, g_xxxyyy_0_xz_1, g_xxxyyy_0_xzz_1, g_xxxyyy_0_yy_1, g_xxxyyy_0_yyy_1, g_xxxyyy_0_yyz_1, g_xxxyyy_0_yz_1, g_xxxyyy_0_yzz_1, g_xxxyyy_0_zz_1, g_xxxyyy_0_zzz_1, g_xxxyyyz_0_xxx_0, g_xxxyyyz_0_xxy_0, g_xxxyyyz_0_xxz_0, g_xxxyyyz_0_xyy_0, g_xxxyyyz_0_xyz_0, g_xxxyyyz_0_xzz_0, g_xxxyyyz_0_yyy_0, g_xxxyyyz_0_yyz_0, g_xxxyyyz_0_yzz_0, g_xxxyyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_xxx_0[i] = g_xxxyyy_0_xxx_1[i] * wa_z[i];

        g_xxxyyyz_0_xxy_0[i] = g_xxxyyy_0_xxy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxz_0[i] = g_xxxyyy_0_xx_1[i] * fi_acd_0 + g_xxxyyy_0_xxz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyy_0[i] = g_xxxyyy_0_xyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xyz_0[i] = g_xxxyyy_0_xy_1[i] * fi_acd_0 + g_xxxyyy_0_xyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xzz_0[i] = 2.0 * g_xxxyyy_0_xz_1[i] * fi_acd_0 + g_xxxyyy_0_xzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyy_0[i] = g_xxxyyy_0_yyy_1[i] * wa_z[i];

        g_xxxyyyz_0_yyz_0[i] = g_xxxyyy_0_yy_1[i] * fi_acd_0 + g_xxxyyy_0_yyz_1[i] * wa_z[i];

        g_xxxyyyz_0_yzz_0[i] = 2.0 * g_xxxyyy_0_yz_1[i] * fi_acd_0 + g_xxxyyy_0_yzz_1[i] * wa_z[i];

        g_xxxyyyz_0_zzz_0[i] = 3.0 * g_xxxyyy_0_zz_1[i] * fi_acd_0 + g_xxxyyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 120-130 components of targeted buffer : KSF

    auto g_xxxyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 120);

    auto g_xxxyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 121);

    auto g_xxxyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 122);

    auto g_xxxyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 123);

    auto g_xxxyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 124);

    auto g_xxxyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 125);

    auto g_xxxyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 126);

    auto g_xxxyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 127);

    auto g_xxxyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 128);

    auto g_xxxyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 129);

    #pragma omp simd aligned(g_xxxyy_0_xxy_0, g_xxxyy_0_xxy_1, g_xxxyy_0_xyy_0, g_xxxyy_0_xyy_1, g_xxxyyz_0_xxy_1, g_xxxyyz_0_xyy_1, g_xxxyyzz_0_xxx_0, g_xxxyyzz_0_xxy_0, g_xxxyyzz_0_xxz_0, g_xxxyyzz_0_xyy_0, g_xxxyyzz_0_xyz_0, g_xxxyyzz_0_xzz_0, g_xxxyyzz_0_yyy_0, g_xxxyyzz_0_yyz_0, g_xxxyyzz_0_yzz_0, g_xxxyyzz_0_zzz_0, g_xxxyzz_0_xxx_1, g_xxxyzz_0_xxz_1, g_xxxyzz_0_xzz_1, g_xxxzz_0_xxx_0, g_xxxzz_0_xxx_1, g_xxxzz_0_xxz_0, g_xxxzz_0_xxz_1, g_xxxzz_0_xzz_0, g_xxxzz_0_xzz_1, g_xxyyzz_0_xyz_1, g_xxyyzz_0_yyy_1, g_xxyyzz_0_yyz_1, g_xxyyzz_0_yz_1, g_xxyyzz_0_yzz_1, g_xxyyzz_0_zzz_1, g_xyyzz_0_xyz_0, g_xyyzz_0_xyz_1, g_xyyzz_0_yyy_0, g_xyyzz_0_yyy_1, g_xyyzz_0_yyz_0, g_xyyzz_0_yyz_1, g_xyyzz_0_yzz_0, g_xyyzz_0_yzz_1, g_xyyzz_0_zzz_0, g_xyyzz_0_zzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyzz_0_xxx_0[i] = g_xxxzz_0_xxx_0[i] * fbe_0 - g_xxxzz_0_xxx_1[i] * fz_be_0 + g_xxxyzz_0_xxx_1[i] * wa_y[i];

        g_xxxyyzz_0_xxy_0[i] = g_xxxyy_0_xxy_0[i] * fbe_0 - g_xxxyy_0_xxy_1[i] * fz_be_0 + g_xxxyyz_0_xxy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxz_0[i] = g_xxxzz_0_xxz_0[i] * fbe_0 - g_xxxzz_0_xxz_1[i] * fz_be_0 + g_xxxyzz_0_xxz_1[i] * wa_y[i];

        g_xxxyyzz_0_xyy_0[i] = g_xxxyy_0_xyy_0[i] * fbe_0 - g_xxxyy_0_xyy_1[i] * fz_be_0 + g_xxxyyz_0_xyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xyz_0[i] = 2.0 * g_xyyzz_0_xyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyz_1[i] * fz_be_0 + g_xxyyzz_0_yz_1[i] * fi_acd_0 + g_xxyyzz_0_xyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xzz_0[i] = g_xxxzz_0_xzz_0[i] * fbe_0 - g_xxxzz_0_xzz_1[i] * fz_be_0 + g_xxxyzz_0_xzz_1[i] * wa_y[i];

        g_xxxyyzz_0_yyy_0[i] = 2.0 * g_xyyzz_0_yyy_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyy_1[i] * fz_be_0 + g_xxyyzz_0_yyy_1[i] * wa_x[i];

        g_xxxyyzz_0_yyz_0[i] = 2.0 * g_xyyzz_0_yyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyz_1[i] * fz_be_0 + g_xxyyzz_0_yyz_1[i] * wa_x[i];

        g_xxxyyzz_0_yzz_0[i] = 2.0 * g_xyyzz_0_yzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yzz_1[i] * fz_be_0 + g_xxyyzz_0_yzz_1[i] * wa_x[i];

        g_xxxyyzz_0_zzz_0[i] = 2.0 * g_xyyzz_0_zzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_zzz_1[i] * fz_be_0 + g_xxyyzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 130-140 components of targeted buffer : KSF

    auto g_xxxyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 130);

    auto g_xxxyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 131);

    auto g_xxxyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 132);

    auto g_xxxyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 133);

    auto g_xxxyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 134);

    auto g_xxxyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 135);

    auto g_xxxyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 136);

    auto g_xxxyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 137);

    auto g_xxxyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 138);

    auto g_xxxyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 139);

    #pragma omp simd aligned(g_xxxyzzz_0_xxx_0, g_xxxyzzz_0_xxy_0, g_xxxyzzz_0_xxz_0, g_xxxyzzz_0_xyy_0, g_xxxyzzz_0_xyz_0, g_xxxyzzz_0_xzz_0, g_xxxyzzz_0_yyy_0, g_xxxyzzz_0_yyz_0, g_xxxyzzz_0_yzz_0, g_xxxyzzz_0_zzz_0, g_xxxzzz_0_xx_1, g_xxxzzz_0_xxx_1, g_xxxzzz_0_xxy_1, g_xxxzzz_0_xxz_1, g_xxxzzz_0_xy_1, g_xxxzzz_0_xyy_1, g_xxxzzz_0_xyz_1, g_xxxzzz_0_xz_1, g_xxxzzz_0_xzz_1, g_xxxzzz_0_yy_1, g_xxxzzz_0_yyy_1, g_xxxzzz_0_yyz_1, g_xxxzzz_0_yz_1, g_xxxzzz_0_yzz_1, g_xxxzzz_0_zz_1, g_xxxzzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_xxx_0[i] = g_xxxzzz_0_xxx_1[i] * wa_y[i];

        g_xxxyzzz_0_xxy_0[i] = g_xxxzzz_0_xx_1[i] * fi_acd_0 + g_xxxzzz_0_xxy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxz_0[i] = g_xxxzzz_0_xxz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyy_0[i] = 2.0 * g_xxxzzz_0_xy_1[i] * fi_acd_0 + g_xxxzzz_0_xyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xyz_0[i] = g_xxxzzz_0_xz_1[i] * fi_acd_0 + g_xxxzzz_0_xyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xzz_0[i] = g_xxxzzz_0_xzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyy_0[i] = 3.0 * g_xxxzzz_0_yy_1[i] * fi_acd_0 + g_xxxzzz_0_yyy_1[i] * wa_y[i];

        g_xxxyzzz_0_yyz_0[i] = 2.0 * g_xxxzzz_0_yz_1[i] * fi_acd_0 + g_xxxzzz_0_yyz_1[i] * wa_y[i];

        g_xxxyzzz_0_yzz_0[i] = g_xxxzzz_0_zz_1[i] * fi_acd_0 + g_xxxzzz_0_yzz_1[i] * wa_y[i];

        g_xxxyzzz_0_zzz_0[i] = g_xxxzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 140-150 components of targeted buffer : KSF

    auto g_xxxzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 140);

    auto g_xxxzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 141);

    auto g_xxxzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 142);

    auto g_xxxzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 143);

    auto g_xxxzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 144);

    auto g_xxxzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 145);

    auto g_xxxzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 146);

    auto g_xxxzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 147);

    auto g_xxxzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 148);

    auto g_xxxzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 149);

    #pragma omp simd aligned(g_xxxzz_0_xxx_0, g_xxxzz_0_xxx_1, g_xxxzz_0_xxy_0, g_xxxzz_0_xxy_1, g_xxxzz_0_xyy_0, g_xxxzz_0_xyy_1, g_xxxzzz_0_xxx_1, g_xxxzzz_0_xxy_1, g_xxxzzz_0_xyy_1, g_xxxzzzz_0_xxx_0, g_xxxzzzz_0_xxy_0, g_xxxzzzz_0_xxz_0, g_xxxzzzz_0_xyy_0, g_xxxzzzz_0_xyz_0, g_xxxzzzz_0_xzz_0, g_xxxzzzz_0_yyy_0, g_xxxzzzz_0_yyz_0, g_xxxzzzz_0_yzz_0, g_xxxzzzz_0_zzz_0, g_xxzzzz_0_xxz_1, g_xxzzzz_0_xyz_1, g_xxzzzz_0_xz_1, g_xxzzzz_0_xzz_1, g_xxzzzz_0_yyy_1, g_xxzzzz_0_yyz_1, g_xxzzzz_0_yz_1, g_xxzzzz_0_yzz_1, g_xxzzzz_0_zz_1, g_xxzzzz_0_zzz_1, g_xzzzz_0_xxz_0, g_xzzzz_0_xxz_1, g_xzzzz_0_xyz_0, g_xzzzz_0_xyz_1, g_xzzzz_0_xzz_0, g_xzzzz_0_xzz_1, g_xzzzz_0_yyy_0, g_xzzzz_0_yyy_1, g_xzzzz_0_yyz_0, g_xzzzz_0_yyz_1, g_xzzzz_0_yzz_0, g_xzzzz_0_yzz_1, g_xzzzz_0_zzz_0, g_xzzzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzzz_0_xxx_0[i] = 3.0 * g_xxxzz_0_xxx_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxx_1[i] * fz_be_0 + g_xxxzzz_0_xxx_1[i] * wa_z[i];

        g_xxxzzzz_0_xxy_0[i] = 3.0 * g_xxxzz_0_xxy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxy_1[i] * fz_be_0 + g_xxxzzz_0_xxy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxz_0[i] = 2.0 * g_xzzzz_0_xxz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xz_1[i] * fi_acd_0 + g_xxzzzz_0_xxz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyy_0[i] = 3.0 * g_xxxzz_0_xyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xyy_1[i] * fz_be_0 + g_xxxzzz_0_xyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xyz_0[i] = 2.0 * g_xzzzz_0_xyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyz_1[i] * fz_be_0 + g_xxzzzz_0_yz_1[i] * fi_acd_0 + g_xxzzzz_0_xyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xzz_0[i] = 2.0 * g_xzzzz_0_xzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xzz_1[i] * fz_be_0 + g_xxzzzz_0_zz_1[i] * fi_acd_0 + g_xxzzzz_0_xzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyy_0[i] = 2.0 * g_xzzzz_0_yyy_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyy_1[i] * fz_be_0 + g_xxzzzz_0_yyy_1[i] * wa_x[i];

        g_xxxzzzz_0_yyz_0[i] = 2.0 * g_xzzzz_0_yyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyz_1[i] * fz_be_0 + g_xxzzzz_0_yyz_1[i] * wa_x[i];

        g_xxxzzzz_0_yzz_0[i] = 2.0 * g_xzzzz_0_yzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yzz_1[i] * fz_be_0 + g_xxzzzz_0_yzz_1[i] * wa_x[i];

        g_xxxzzzz_0_zzz_0[i] = 2.0 * g_xzzzz_0_zzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_zzz_1[i] * fz_be_0 + g_xxzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 150-160 components of targeted buffer : KSF

    auto g_xxyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 150);

    auto g_xxyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 151);

    auto g_xxyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 152);

    auto g_xxyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 153);

    auto g_xxyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 154);

    auto g_xxyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 155);

    auto g_xxyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 156);

    auto g_xxyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 157);

    auto g_xxyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 158);

    auto g_xxyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 159);

    #pragma omp simd aligned(g_xxyyy_0_xxx_0, g_xxyyy_0_xxx_1, g_xxyyy_0_xxz_0, g_xxyyy_0_xxz_1, g_xxyyy_0_xzz_0, g_xxyyy_0_xzz_1, g_xxyyyy_0_xxx_1, g_xxyyyy_0_xxz_1, g_xxyyyy_0_xzz_1, g_xxyyyyy_0_xxx_0, g_xxyyyyy_0_xxy_0, g_xxyyyyy_0_xxz_0, g_xxyyyyy_0_xyy_0, g_xxyyyyy_0_xyz_0, g_xxyyyyy_0_xzz_0, g_xxyyyyy_0_yyy_0, g_xxyyyyy_0_yyz_0, g_xxyyyyy_0_yzz_0, g_xxyyyyy_0_zzz_0, g_xyyyyy_0_xxy_1, g_xyyyyy_0_xy_1, g_xyyyyy_0_xyy_1, g_xyyyyy_0_xyz_1, g_xyyyyy_0_yy_1, g_xyyyyy_0_yyy_1, g_xyyyyy_0_yyz_1, g_xyyyyy_0_yz_1, g_xyyyyy_0_yzz_1, g_xyyyyy_0_zzz_1, g_yyyyy_0_xxy_0, g_yyyyy_0_xxy_1, g_yyyyy_0_xyy_0, g_yyyyy_0_xyy_1, g_yyyyy_0_xyz_0, g_yyyyy_0_xyz_1, g_yyyyy_0_yyy_0, g_yyyyy_0_yyy_1, g_yyyyy_0_yyz_0, g_yyyyy_0_yyz_1, g_yyyyy_0_yzz_0, g_yyyyy_0_yzz_1, g_yyyyy_0_zzz_0, g_yyyyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyyy_0_xxx_0[i] = 4.0 * g_xxyyy_0_xxx_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxx_1[i] * fz_be_0 + g_xxyyyy_0_xxx_1[i] * wa_y[i];

        g_xxyyyyy_0_xxy_0[i] = g_yyyyy_0_xxy_0[i] * fbe_0 - g_yyyyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xy_1[i] * fi_acd_0 + g_xyyyyy_0_xxy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxz_0[i] = 4.0 * g_xxyyy_0_xxz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxz_1[i] * fz_be_0 + g_xxyyyy_0_xxz_1[i] * wa_y[i];

        g_xxyyyyy_0_xyy_0[i] = g_yyyyy_0_xyy_0[i] * fbe_0 - g_yyyyy_0_xyy_1[i] * fz_be_0 + g_xyyyyy_0_yy_1[i] * fi_acd_0 + g_xyyyyy_0_xyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xyz_0[i] = g_yyyyy_0_xyz_0[i] * fbe_0 - g_yyyyy_0_xyz_1[i] * fz_be_0 + g_xyyyyy_0_yz_1[i] * fi_acd_0 + g_xyyyyy_0_xyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xzz_0[i] = 4.0 * g_xxyyy_0_xzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xzz_1[i] * fz_be_0 + g_xxyyyy_0_xzz_1[i] * wa_y[i];

        g_xxyyyyy_0_yyy_0[i] = g_yyyyy_0_yyy_0[i] * fbe_0 - g_yyyyy_0_yyy_1[i] * fz_be_0 + g_xyyyyy_0_yyy_1[i] * wa_x[i];

        g_xxyyyyy_0_yyz_0[i] = g_yyyyy_0_yyz_0[i] * fbe_0 - g_yyyyy_0_yyz_1[i] * fz_be_0 + g_xyyyyy_0_yyz_1[i] * wa_x[i];

        g_xxyyyyy_0_yzz_0[i] = g_yyyyy_0_yzz_0[i] * fbe_0 - g_yyyyy_0_yzz_1[i] * fz_be_0 + g_xyyyyy_0_yzz_1[i] * wa_x[i];

        g_xxyyyyy_0_zzz_0[i] = g_yyyyy_0_zzz_0[i] * fbe_0 - g_yyyyy_0_zzz_1[i] * fz_be_0 + g_xyyyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 160-170 components of targeted buffer : KSF

    auto g_xxyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 160);

    auto g_xxyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 161);

    auto g_xxyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 162);

    auto g_xxyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 163);

    auto g_xxyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 164);

    auto g_xxyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 165);

    auto g_xxyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 166);

    auto g_xxyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 167);

    auto g_xxyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 168);

    auto g_xxyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 169);

    #pragma omp simd aligned(g_xxyyyy_0_xx_1, g_xxyyyy_0_xxx_1, g_xxyyyy_0_xxy_1, g_xxyyyy_0_xxz_1, g_xxyyyy_0_xy_1, g_xxyyyy_0_xyy_1, g_xxyyyy_0_xyz_1, g_xxyyyy_0_xz_1, g_xxyyyy_0_xzz_1, g_xxyyyy_0_yy_1, g_xxyyyy_0_yyy_1, g_xxyyyy_0_yyz_1, g_xxyyyy_0_yz_1, g_xxyyyy_0_yzz_1, g_xxyyyy_0_zz_1, g_xxyyyy_0_zzz_1, g_xxyyyyz_0_xxx_0, g_xxyyyyz_0_xxy_0, g_xxyyyyz_0_xxz_0, g_xxyyyyz_0_xyy_0, g_xxyyyyz_0_xyz_0, g_xxyyyyz_0_xzz_0, g_xxyyyyz_0_yyy_0, g_xxyyyyz_0_yyz_0, g_xxyyyyz_0_yzz_0, g_xxyyyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_xxx_0[i] = g_xxyyyy_0_xxx_1[i] * wa_z[i];

        g_xxyyyyz_0_xxy_0[i] = g_xxyyyy_0_xxy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxz_0[i] = g_xxyyyy_0_xx_1[i] * fi_acd_0 + g_xxyyyy_0_xxz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyy_0[i] = g_xxyyyy_0_xyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xyz_0[i] = g_xxyyyy_0_xy_1[i] * fi_acd_0 + g_xxyyyy_0_xyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xzz_0[i] = 2.0 * g_xxyyyy_0_xz_1[i] * fi_acd_0 + g_xxyyyy_0_xzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyy_0[i] = g_xxyyyy_0_yyy_1[i] * wa_z[i];

        g_xxyyyyz_0_yyz_0[i] = g_xxyyyy_0_yy_1[i] * fi_acd_0 + g_xxyyyy_0_yyz_1[i] * wa_z[i];

        g_xxyyyyz_0_yzz_0[i] = 2.0 * g_xxyyyy_0_yz_1[i] * fi_acd_0 + g_xxyyyy_0_yzz_1[i] * wa_z[i];

        g_xxyyyyz_0_zzz_0[i] = 3.0 * g_xxyyyy_0_zz_1[i] * fi_acd_0 + g_xxyyyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 170-180 components of targeted buffer : KSF

    auto g_xxyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 170);

    auto g_xxyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 171);

    auto g_xxyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 172);

    auto g_xxyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 173);

    auto g_xxyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 174);

    auto g_xxyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 175);

    auto g_xxyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 176);

    auto g_xxyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 177);

    auto g_xxyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 178);

    auto g_xxyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 179);

    #pragma omp simd aligned(g_xxyyy_0_xxy_0, g_xxyyy_0_xxy_1, g_xxyyy_0_xyy_0, g_xxyyy_0_xyy_1, g_xxyyyz_0_xxy_1, g_xxyyyz_0_xyy_1, g_xxyyyzz_0_xxx_0, g_xxyyyzz_0_xxy_0, g_xxyyyzz_0_xxz_0, g_xxyyyzz_0_xyy_0, g_xxyyyzz_0_xyz_0, g_xxyyyzz_0_xzz_0, g_xxyyyzz_0_yyy_0, g_xxyyyzz_0_yyz_0, g_xxyyyzz_0_yzz_0, g_xxyyyzz_0_zzz_0, g_xxyyzz_0_xxx_1, g_xxyyzz_0_xxz_1, g_xxyyzz_0_xzz_1, g_xxyzz_0_xxx_0, g_xxyzz_0_xxx_1, g_xxyzz_0_xxz_0, g_xxyzz_0_xxz_1, g_xxyzz_0_xzz_0, g_xxyzz_0_xzz_1, g_xyyyzz_0_xyz_1, g_xyyyzz_0_yyy_1, g_xyyyzz_0_yyz_1, g_xyyyzz_0_yz_1, g_xyyyzz_0_yzz_1, g_xyyyzz_0_zzz_1, g_yyyzz_0_xyz_0, g_yyyzz_0_xyz_1, g_yyyzz_0_yyy_0, g_yyyzz_0_yyy_1, g_yyyzz_0_yyz_0, g_yyyzz_0_yyz_1, g_yyyzz_0_yzz_0, g_yyyzz_0_yzz_1, g_yyyzz_0_zzz_0, g_yyyzz_0_zzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyzz_0_xxx_0[i] = 2.0 * g_xxyzz_0_xxx_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxx_1[i] * fz_be_0 + g_xxyyzz_0_xxx_1[i] * wa_y[i];

        g_xxyyyzz_0_xxy_0[i] = g_xxyyy_0_xxy_0[i] * fbe_0 - g_xxyyy_0_xxy_1[i] * fz_be_0 + g_xxyyyz_0_xxy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxz_0[i] = 2.0 * g_xxyzz_0_xxz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxz_1[i] * fz_be_0 + g_xxyyzz_0_xxz_1[i] * wa_y[i];

        g_xxyyyzz_0_xyy_0[i] = g_xxyyy_0_xyy_0[i] * fbe_0 - g_xxyyy_0_xyy_1[i] * fz_be_0 + g_xxyyyz_0_xyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xyz_0[i] = g_yyyzz_0_xyz_0[i] * fbe_0 - g_yyyzz_0_xyz_1[i] * fz_be_0 + g_xyyyzz_0_yz_1[i] * fi_acd_0 + g_xyyyzz_0_xyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xzz_0[i] = 2.0 * g_xxyzz_0_xzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xzz_1[i] * fz_be_0 + g_xxyyzz_0_xzz_1[i] * wa_y[i];

        g_xxyyyzz_0_yyy_0[i] = g_yyyzz_0_yyy_0[i] * fbe_0 - g_yyyzz_0_yyy_1[i] * fz_be_0 + g_xyyyzz_0_yyy_1[i] * wa_x[i];

        g_xxyyyzz_0_yyz_0[i] = g_yyyzz_0_yyz_0[i] * fbe_0 - g_yyyzz_0_yyz_1[i] * fz_be_0 + g_xyyyzz_0_yyz_1[i] * wa_x[i];

        g_xxyyyzz_0_yzz_0[i] = g_yyyzz_0_yzz_0[i] * fbe_0 - g_yyyzz_0_yzz_1[i] * fz_be_0 + g_xyyyzz_0_yzz_1[i] * wa_x[i];

        g_xxyyyzz_0_zzz_0[i] = g_yyyzz_0_zzz_0[i] * fbe_0 - g_yyyzz_0_zzz_1[i] * fz_be_0 + g_xyyyzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 180-190 components of targeted buffer : KSF

    auto g_xxyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 180);

    auto g_xxyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 181);

    auto g_xxyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 182);

    auto g_xxyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 183);

    auto g_xxyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 184);

    auto g_xxyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 185);

    auto g_xxyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 186);

    auto g_xxyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 187);

    auto g_xxyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 188);

    auto g_xxyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 189);

    #pragma omp simd aligned(g_xxyyz_0_xxy_0, g_xxyyz_0_xxy_1, g_xxyyz_0_xyy_0, g_xxyyz_0_xyy_1, g_xxyyzz_0_xxy_1, g_xxyyzz_0_xyy_1, g_xxyyzzz_0_xxx_0, g_xxyyzzz_0_xxy_0, g_xxyyzzz_0_xxz_0, g_xxyyzzz_0_xyy_0, g_xxyyzzz_0_xyz_0, g_xxyyzzz_0_xzz_0, g_xxyyzzz_0_yyy_0, g_xxyyzzz_0_yyz_0, g_xxyyzzz_0_yzz_0, g_xxyyzzz_0_zzz_0, g_xxyzzz_0_xxx_1, g_xxyzzz_0_xxz_1, g_xxyzzz_0_xzz_1, g_xxzzz_0_xxx_0, g_xxzzz_0_xxx_1, g_xxzzz_0_xxz_0, g_xxzzz_0_xxz_1, g_xxzzz_0_xzz_0, g_xxzzz_0_xzz_1, g_xyyzzz_0_xyz_1, g_xyyzzz_0_yyy_1, g_xyyzzz_0_yyz_1, g_xyyzzz_0_yz_1, g_xyyzzz_0_yzz_1, g_xyyzzz_0_zzz_1, g_yyzzz_0_xyz_0, g_yyzzz_0_xyz_1, g_yyzzz_0_yyy_0, g_yyzzz_0_yyy_1, g_yyzzz_0_yyz_0, g_yyzzz_0_yyz_1, g_yyzzz_0_yzz_0, g_yyzzz_0_yzz_1, g_yyzzz_0_zzz_0, g_yyzzz_0_zzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzzz_0_xxx_0[i] = g_xxzzz_0_xxx_0[i] * fbe_0 - g_xxzzz_0_xxx_1[i] * fz_be_0 + g_xxyzzz_0_xxx_1[i] * wa_y[i];

        g_xxyyzzz_0_xxy_0[i] = 2.0 * g_xxyyz_0_xxy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxy_1[i] * fz_be_0 + g_xxyyzz_0_xxy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxz_0[i] = g_xxzzz_0_xxz_0[i] * fbe_0 - g_xxzzz_0_xxz_1[i] * fz_be_0 + g_xxyzzz_0_xxz_1[i] * wa_y[i];

        g_xxyyzzz_0_xyy_0[i] = 2.0 * g_xxyyz_0_xyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xyy_1[i] * fz_be_0 + g_xxyyzz_0_xyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xyz_0[i] = g_yyzzz_0_xyz_0[i] * fbe_0 - g_yyzzz_0_xyz_1[i] * fz_be_0 + g_xyyzzz_0_yz_1[i] * fi_acd_0 + g_xyyzzz_0_xyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xzz_0[i] = g_xxzzz_0_xzz_0[i] * fbe_0 - g_xxzzz_0_xzz_1[i] * fz_be_0 + g_xxyzzz_0_xzz_1[i] * wa_y[i];

        g_xxyyzzz_0_yyy_0[i] = g_yyzzz_0_yyy_0[i] * fbe_0 - g_yyzzz_0_yyy_1[i] * fz_be_0 + g_xyyzzz_0_yyy_1[i] * wa_x[i];

        g_xxyyzzz_0_yyz_0[i] = g_yyzzz_0_yyz_0[i] * fbe_0 - g_yyzzz_0_yyz_1[i] * fz_be_0 + g_xyyzzz_0_yyz_1[i] * wa_x[i];

        g_xxyyzzz_0_yzz_0[i] = g_yyzzz_0_yzz_0[i] * fbe_0 - g_yyzzz_0_yzz_1[i] * fz_be_0 + g_xyyzzz_0_yzz_1[i] * wa_x[i];

        g_xxyyzzz_0_zzz_0[i] = g_yyzzz_0_zzz_0[i] * fbe_0 - g_yyzzz_0_zzz_1[i] * fz_be_0 + g_xyyzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 190-200 components of targeted buffer : KSF

    auto g_xxyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 190);

    auto g_xxyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 191);

    auto g_xxyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 192);

    auto g_xxyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 193);

    auto g_xxyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 194);

    auto g_xxyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 195);

    auto g_xxyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 196);

    auto g_xxyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 197);

    auto g_xxyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 198);

    auto g_xxyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 199);

    #pragma omp simd aligned(g_xxyzzzz_0_xxx_0, g_xxyzzzz_0_xxy_0, g_xxyzzzz_0_xxz_0, g_xxyzzzz_0_xyy_0, g_xxyzzzz_0_xyz_0, g_xxyzzzz_0_xzz_0, g_xxyzzzz_0_yyy_0, g_xxyzzzz_0_yyz_0, g_xxyzzzz_0_yzz_0, g_xxyzzzz_0_zzz_0, g_xxzzzz_0_xx_1, g_xxzzzz_0_xxx_1, g_xxzzzz_0_xxy_1, g_xxzzzz_0_xxz_1, g_xxzzzz_0_xy_1, g_xxzzzz_0_xyy_1, g_xxzzzz_0_xyz_1, g_xxzzzz_0_xz_1, g_xxzzzz_0_xzz_1, g_xxzzzz_0_yy_1, g_xxzzzz_0_yyy_1, g_xxzzzz_0_yyz_1, g_xxzzzz_0_yz_1, g_xxzzzz_0_yzz_1, g_xxzzzz_0_zz_1, g_xxzzzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_xxx_0[i] = g_xxzzzz_0_xxx_1[i] * wa_y[i];

        g_xxyzzzz_0_xxy_0[i] = g_xxzzzz_0_xx_1[i] * fi_acd_0 + g_xxzzzz_0_xxy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxz_0[i] = g_xxzzzz_0_xxz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyy_0[i] = 2.0 * g_xxzzzz_0_xy_1[i] * fi_acd_0 + g_xxzzzz_0_xyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xyz_0[i] = g_xxzzzz_0_xz_1[i] * fi_acd_0 + g_xxzzzz_0_xyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xzz_0[i] = g_xxzzzz_0_xzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyy_0[i] = 3.0 * g_xxzzzz_0_yy_1[i] * fi_acd_0 + g_xxzzzz_0_yyy_1[i] * wa_y[i];

        g_xxyzzzz_0_yyz_0[i] = 2.0 * g_xxzzzz_0_yz_1[i] * fi_acd_0 + g_xxzzzz_0_yyz_1[i] * wa_y[i];

        g_xxyzzzz_0_yzz_0[i] = g_xxzzzz_0_zz_1[i] * fi_acd_0 + g_xxzzzz_0_yzz_1[i] * wa_y[i];

        g_xxyzzzz_0_zzz_0[i] = g_xxzzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 200-210 components of targeted buffer : KSF

    auto g_xxzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 200);

    auto g_xxzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 201);

    auto g_xxzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 202);

    auto g_xxzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 203);

    auto g_xxzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 204);

    auto g_xxzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 205);

    auto g_xxzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 206);

    auto g_xxzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 207);

    auto g_xxzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 208);

    auto g_xxzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 209);

    #pragma omp simd aligned(g_xxzzz_0_xxx_0, g_xxzzz_0_xxx_1, g_xxzzz_0_xxy_0, g_xxzzz_0_xxy_1, g_xxzzz_0_xyy_0, g_xxzzz_0_xyy_1, g_xxzzzz_0_xxx_1, g_xxzzzz_0_xxy_1, g_xxzzzz_0_xyy_1, g_xxzzzzz_0_xxx_0, g_xxzzzzz_0_xxy_0, g_xxzzzzz_0_xxz_0, g_xxzzzzz_0_xyy_0, g_xxzzzzz_0_xyz_0, g_xxzzzzz_0_xzz_0, g_xxzzzzz_0_yyy_0, g_xxzzzzz_0_yyz_0, g_xxzzzzz_0_yzz_0, g_xxzzzzz_0_zzz_0, g_xzzzzz_0_xxz_1, g_xzzzzz_0_xyz_1, g_xzzzzz_0_xz_1, g_xzzzzz_0_xzz_1, g_xzzzzz_0_yyy_1, g_xzzzzz_0_yyz_1, g_xzzzzz_0_yz_1, g_xzzzzz_0_yzz_1, g_xzzzzz_0_zz_1, g_xzzzzz_0_zzz_1, g_zzzzz_0_xxz_0, g_zzzzz_0_xxz_1, g_zzzzz_0_xyz_0, g_zzzzz_0_xyz_1, g_zzzzz_0_xzz_0, g_zzzzz_0_xzz_1, g_zzzzz_0_yyy_0, g_zzzzz_0_yyy_1, g_zzzzz_0_yyz_0, g_zzzzz_0_yyz_1, g_zzzzz_0_yzz_0, g_zzzzz_0_yzz_1, g_zzzzz_0_zzz_0, g_zzzzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzzz_0_xxx_0[i] = 4.0 * g_xxzzz_0_xxx_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxx_1[i] * fz_be_0 + g_xxzzzz_0_xxx_1[i] * wa_z[i];

        g_xxzzzzz_0_xxy_0[i] = 4.0 * g_xxzzz_0_xxy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxy_1[i] * fz_be_0 + g_xxzzzz_0_xxy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxz_0[i] = g_zzzzz_0_xxz_0[i] * fbe_0 - g_zzzzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xz_1[i] * fi_acd_0 + g_xzzzzz_0_xxz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyy_0[i] = 4.0 * g_xxzzz_0_xyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xyy_1[i] * fz_be_0 + g_xxzzzz_0_xyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xyz_0[i] = g_zzzzz_0_xyz_0[i] * fbe_0 - g_zzzzz_0_xyz_1[i] * fz_be_0 + g_xzzzzz_0_yz_1[i] * fi_acd_0 + g_xzzzzz_0_xyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xzz_0[i] = g_zzzzz_0_xzz_0[i] * fbe_0 - g_zzzzz_0_xzz_1[i] * fz_be_0 + g_xzzzzz_0_zz_1[i] * fi_acd_0 + g_xzzzzz_0_xzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyy_0[i] = g_zzzzz_0_yyy_0[i] * fbe_0 - g_zzzzz_0_yyy_1[i] * fz_be_0 + g_xzzzzz_0_yyy_1[i] * wa_x[i];

        g_xxzzzzz_0_yyz_0[i] = g_zzzzz_0_yyz_0[i] * fbe_0 - g_zzzzz_0_yyz_1[i] * fz_be_0 + g_xzzzzz_0_yyz_1[i] * wa_x[i];

        g_xxzzzzz_0_yzz_0[i] = g_zzzzz_0_yzz_0[i] * fbe_0 - g_zzzzz_0_yzz_1[i] * fz_be_0 + g_xzzzzz_0_yzz_1[i] * wa_x[i];

        g_xxzzzzz_0_zzz_0[i] = g_zzzzz_0_zzz_0[i] * fbe_0 - g_zzzzz_0_zzz_1[i] * fz_be_0 + g_xzzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 210-220 components of targeted buffer : KSF

    auto g_xyyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 210);

    auto g_xyyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 211);

    auto g_xyyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 212);

    auto g_xyyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 213);

    auto g_xyyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 214);

    auto g_xyyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 215);

    auto g_xyyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 216);

    auto g_xyyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 217);

    auto g_xyyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 218);

    auto g_xyyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 219);

    #pragma omp simd aligned(g_xyyyyyy_0_xxx_0, g_xyyyyyy_0_xxy_0, g_xyyyyyy_0_xxz_0, g_xyyyyyy_0_xyy_0, g_xyyyyyy_0_xyz_0, g_xyyyyyy_0_xzz_0, g_xyyyyyy_0_yyy_0, g_xyyyyyy_0_yyz_0, g_xyyyyyy_0_yzz_0, g_xyyyyyy_0_zzz_0, g_yyyyyy_0_xx_1, g_yyyyyy_0_xxx_1, g_yyyyyy_0_xxy_1, g_yyyyyy_0_xxz_1, g_yyyyyy_0_xy_1, g_yyyyyy_0_xyy_1, g_yyyyyy_0_xyz_1, g_yyyyyy_0_xz_1, g_yyyyyy_0_xzz_1, g_yyyyyy_0_yy_1, g_yyyyyy_0_yyy_1, g_yyyyyy_0_yyz_1, g_yyyyyy_0_yz_1, g_yyyyyy_0_yzz_1, g_yyyyyy_0_zz_1, g_yyyyyy_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_xxx_0[i] = 3.0 * g_yyyyyy_0_xx_1[i] * fi_acd_0 + g_yyyyyy_0_xxx_1[i] * wa_x[i];

        g_xyyyyyy_0_xxy_0[i] = 2.0 * g_yyyyyy_0_xy_1[i] * fi_acd_0 + g_yyyyyy_0_xxy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxz_0[i] = 2.0 * g_yyyyyy_0_xz_1[i] * fi_acd_0 + g_yyyyyy_0_xxz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyy_0[i] = g_yyyyyy_0_yy_1[i] * fi_acd_0 + g_yyyyyy_0_xyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xyz_0[i] = g_yyyyyy_0_yz_1[i] * fi_acd_0 + g_yyyyyy_0_xyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xzz_0[i] = g_yyyyyy_0_zz_1[i] * fi_acd_0 + g_yyyyyy_0_xzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyy_0[i] = g_yyyyyy_0_yyy_1[i] * wa_x[i];

        g_xyyyyyy_0_yyz_0[i] = g_yyyyyy_0_yyz_1[i] * wa_x[i];

        g_xyyyyyy_0_yzz_0[i] = g_yyyyyy_0_yzz_1[i] * wa_x[i];

        g_xyyyyyy_0_zzz_0[i] = g_yyyyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 220-230 components of targeted buffer : KSF

    auto g_xyyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 220);

    auto g_xyyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 221);

    auto g_xyyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 222);

    auto g_xyyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 223);

    auto g_xyyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 224);

    auto g_xyyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 225);

    auto g_xyyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 226);

    auto g_xyyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 227);

    auto g_xyyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 228);

    auto g_xyyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 229);

    #pragma omp simd aligned(g_xyyyyy_0_xxx_1, g_xyyyyy_0_xxy_1, g_xyyyyy_0_xyy_1, g_xyyyyyz_0_xxx_0, g_xyyyyyz_0_xxy_0, g_xyyyyyz_0_xxz_0, g_xyyyyyz_0_xyy_0, g_xyyyyyz_0_xyz_0, g_xyyyyyz_0_xzz_0, g_xyyyyyz_0_yyy_0, g_xyyyyyz_0_yyz_0, g_xyyyyyz_0_yzz_0, g_xyyyyyz_0_zzz_0, g_yyyyyz_0_xxz_1, g_yyyyyz_0_xyz_1, g_yyyyyz_0_xz_1, g_yyyyyz_0_xzz_1, g_yyyyyz_0_yyy_1, g_yyyyyz_0_yyz_1, g_yyyyyz_0_yz_1, g_yyyyyz_0_yzz_1, g_yyyyyz_0_zz_1, g_yyyyyz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyz_0_xxx_0[i] = g_xyyyyy_0_xxx_1[i] * wa_z[i];

        g_xyyyyyz_0_xxy_0[i] = g_xyyyyy_0_xxy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxz_0[i] = 2.0 * g_yyyyyz_0_xz_1[i] * fi_acd_0 + g_yyyyyz_0_xxz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyy_0[i] = g_xyyyyy_0_xyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xyz_0[i] = g_yyyyyz_0_yz_1[i] * fi_acd_0 + g_yyyyyz_0_xyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xzz_0[i] = g_yyyyyz_0_zz_1[i] * fi_acd_0 + g_yyyyyz_0_xzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyy_0[i] = g_yyyyyz_0_yyy_1[i] * wa_x[i];

        g_xyyyyyz_0_yyz_0[i] = g_yyyyyz_0_yyz_1[i] * wa_x[i];

        g_xyyyyyz_0_yzz_0[i] = g_yyyyyz_0_yzz_1[i] * wa_x[i];

        g_xyyyyyz_0_zzz_0[i] = g_yyyyyz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 230-240 components of targeted buffer : KSF

    auto g_xyyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 230);

    auto g_xyyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 231);

    auto g_xyyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 232);

    auto g_xyyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 233);

    auto g_xyyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 234);

    auto g_xyyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 235);

    auto g_xyyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 236);

    auto g_xyyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 237);

    auto g_xyyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 238);

    auto g_xyyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 239);

    #pragma omp simd aligned(g_xyyyyzz_0_xxx_0, g_xyyyyzz_0_xxy_0, g_xyyyyzz_0_xxz_0, g_xyyyyzz_0_xyy_0, g_xyyyyzz_0_xyz_0, g_xyyyyzz_0_xzz_0, g_xyyyyzz_0_yyy_0, g_xyyyyzz_0_yyz_0, g_xyyyyzz_0_yzz_0, g_xyyyyzz_0_zzz_0, g_yyyyzz_0_xx_1, g_yyyyzz_0_xxx_1, g_yyyyzz_0_xxy_1, g_yyyyzz_0_xxz_1, g_yyyyzz_0_xy_1, g_yyyyzz_0_xyy_1, g_yyyyzz_0_xyz_1, g_yyyyzz_0_xz_1, g_yyyyzz_0_xzz_1, g_yyyyzz_0_yy_1, g_yyyyzz_0_yyy_1, g_yyyyzz_0_yyz_1, g_yyyyzz_0_yz_1, g_yyyyzz_0_yzz_1, g_yyyyzz_0_zz_1, g_yyyyzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_xxx_0[i] = 3.0 * g_yyyyzz_0_xx_1[i] * fi_acd_0 + g_yyyyzz_0_xxx_1[i] * wa_x[i];

        g_xyyyyzz_0_xxy_0[i] = 2.0 * g_yyyyzz_0_xy_1[i] * fi_acd_0 + g_yyyyzz_0_xxy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxz_0[i] = 2.0 * g_yyyyzz_0_xz_1[i] * fi_acd_0 + g_yyyyzz_0_xxz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyy_0[i] = g_yyyyzz_0_yy_1[i] * fi_acd_0 + g_yyyyzz_0_xyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xyz_0[i] = g_yyyyzz_0_yz_1[i] * fi_acd_0 + g_yyyyzz_0_xyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xzz_0[i] = g_yyyyzz_0_zz_1[i] * fi_acd_0 + g_yyyyzz_0_xzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyy_0[i] = g_yyyyzz_0_yyy_1[i] * wa_x[i];

        g_xyyyyzz_0_yyz_0[i] = g_yyyyzz_0_yyz_1[i] * wa_x[i];

        g_xyyyyzz_0_yzz_0[i] = g_yyyyzz_0_yzz_1[i] * wa_x[i];

        g_xyyyyzz_0_zzz_0[i] = g_yyyyzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 240-250 components of targeted buffer : KSF

    auto g_xyyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 240);

    auto g_xyyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 241);

    auto g_xyyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 242);

    auto g_xyyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 243);

    auto g_xyyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 244);

    auto g_xyyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 245);

    auto g_xyyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 246);

    auto g_xyyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 247);

    auto g_xyyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 248);

    auto g_xyyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 249);

    #pragma omp simd aligned(g_xyyyzzz_0_xxx_0, g_xyyyzzz_0_xxy_0, g_xyyyzzz_0_xxz_0, g_xyyyzzz_0_xyy_0, g_xyyyzzz_0_xyz_0, g_xyyyzzz_0_xzz_0, g_xyyyzzz_0_yyy_0, g_xyyyzzz_0_yyz_0, g_xyyyzzz_0_yzz_0, g_xyyyzzz_0_zzz_0, g_yyyzzz_0_xx_1, g_yyyzzz_0_xxx_1, g_yyyzzz_0_xxy_1, g_yyyzzz_0_xxz_1, g_yyyzzz_0_xy_1, g_yyyzzz_0_xyy_1, g_yyyzzz_0_xyz_1, g_yyyzzz_0_xz_1, g_yyyzzz_0_xzz_1, g_yyyzzz_0_yy_1, g_yyyzzz_0_yyy_1, g_yyyzzz_0_yyz_1, g_yyyzzz_0_yz_1, g_yyyzzz_0_yzz_1, g_yyyzzz_0_zz_1, g_yyyzzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_xxx_0[i] = 3.0 * g_yyyzzz_0_xx_1[i] * fi_acd_0 + g_yyyzzz_0_xxx_1[i] * wa_x[i];

        g_xyyyzzz_0_xxy_0[i] = 2.0 * g_yyyzzz_0_xy_1[i] * fi_acd_0 + g_yyyzzz_0_xxy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxz_0[i] = 2.0 * g_yyyzzz_0_xz_1[i] * fi_acd_0 + g_yyyzzz_0_xxz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyy_0[i] = g_yyyzzz_0_yy_1[i] * fi_acd_0 + g_yyyzzz_0_xyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xyz_0[i] = g_yyyzzz_0_yz_1[i] * fi_acd_0 + g_yyyzzz_0_xyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xzz_0[i] = g_yyyzzz_0_zz_1[i] * fi_acd_0 + g_yyyzzz_0_xzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyy_0[i] = g_yyyzzz_0_yyy_1[i] * wa_x[i];

        g_xyyyzzz_0_yyz_0[i] = g_yyyzzz_0_yyz_1[i] * wa_x[i];

        g_xyyyzzz_0_yzz_0[i] = g_yyyzzz_0_yzz_1[i] * wa_x[i];

        g_xyyyzzz_0_zzz_0[i] = g_yyyzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 250-260 components of targeted buffer : KSF

    auto g_xyyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 250);

    auto g_xyyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 251);

    auto g_xyyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 252);

    auto g_xyyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 253);

    auto g_xyyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 254);

    auto g_xyyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 255);

    auto g_xyyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 256);

    auto g_xyyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 257);

    auto g_xyyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 258);

    auto g_xyyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 259);

    #pragma omp simd aligned(g_xyyzzzz_0_xxx_0, g_xyyzzzz_0_xxy_0, g_xyyzzzz_0_xxz_0, g_xyyzzzz_0_xyy_0, g_xyyzzzz_0_xyz_0, g_xyyzzzz_0_xzz_0, g_xyyzzzz_0_yyy_0, g_xyyzzzz_0_yyz_0, g_xyyzzzz_0_yzz_0, g_xyyzzzz_0_zzz_0, g_yyzzzz_0_xx_1, g_yyzzzz_0_xxx_1, g_yyzzzz_0_xxy_1, g_yyzzzz_0_xxz_1, g_yyzzzz_0_xy_1, g_yyzzzz_0_xyy_1, g_yyzzzz_0_xyz_1, g_yyzzzz_0_xz_1, g_yyzzzz_0_xzz_1, g_yyzzzz_0_yy_1, g_yyzzzz_0_yyy_1, g_yyzzzz_0_yyz_1, g_yyzzzz_0_yz_1, g_yyzzzz_0_yzz_1, g_yyzzzz_0_zz_1, g_yyzzzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_xxx_0[i] = 3.0 * g_yyzzzz_0_xx_1[i] * fi_acd_0 + g_yyzzzz_0_xxx_1[i] * wa_x[i];

        g_xyyzzzz_0_xxy_0[i] = 2.0 * g_yyzzzz_0_xy_1[i] * fi_acd_0 + g_yyzzzz_0_xxy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxz_0[i] = 2.0 * g_yyzzzz_0_xz_1[i] * fi_acd_0 + g_yyzzzz_0_xxz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyy_0[i] = g_yyzzzz_0_yy_1[i] * fi_acd_0 + g_yyzzzz_0_xyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xyz_0[i] = g_yyzzzz_0_yz_1[i] * fi_acd_0 + g_yyzzzz_0_xyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xzz_0[i] = g_yyzzzz_0_zz_1[i] * fi_acd_0 + g_yyzzzz_0_xzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyy_0[i] = g_yyzzzz_0_yyy_1[i] * wa_x[i];

        g_xyyzzzz_0_yyz_0[i] = g_yyzzzz_0_yyz_1[i] * wa_x[i];

        g_xyyzzzz_0_yzz_0[i] = g_yyzzzz_0_yzz_1[i] * wa_x[i];

        g_xyyzzzz_0_zzz_0[i] = g_yyzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 260-270 components of targeted buffer : KSF

    auto g_xyzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 260);

    auto g_xyzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 261);

    auto g_xyzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 262);

    auto g_xyzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 263);

    auto g_xyzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 264);

    auto g_xyzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 265);

    auto g_xyzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 266);

    auto g_xyzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 267);

    auto g_xyzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 268);

    auto g_xyzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 269);

    #pragma omp simd aligned(g_xyzzzzz_0_xxx_0, g_xyzzzzz_0_xxy_0, g_xyzzzzz_0_xxz_0, g_xyzzzzz_0_xyy_0, g_xyzzzzz_0_xyz_0, g_xyzzzzz_0_xzz_0, g_xyzzzzz_0_yyy_0, g_xyzzzzz_0_yyz_0, g_xyzzzzz_0_yzz_0, g_xyzzzzz_0_zzz_0, g_xzzzzz_0_xxx_1, g_xzzzzz_0_xxz_1, g_xzzzzz_0_xzz_1, g_yzzzzz_0_xxy_1, g_yzzzzz_0_xy_1, g_yzzzzz_0_xyy_1, g_yzzzzz_0_xyz_1, g_yzzzzz_0_yy_1, g_yzzzzz_0_yyy_1, g_yzzzzz_0_yyz_1, g_yzzzzz_0_yz_1, g_yzzzzz_0_yzz_1, g_yzzzzz_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzzz_0_xxx_0[i] = g_xzzzzz_0_xxx_1[i] * wa_y[i];

        g_xyzzzzz_0_xxy_0[i] = 2.0 * g_yzzzzz_0_xy_1[i] * fi_acd_0 + g_yzzzzz_0_xxy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxz_0[i] = g_xzzzzz_0_xxz_1[i] * wa_y[i];

        g_xyzzzzz_0_xyy_0[i] = g_yzzzzz_0_yy_1[i] * fi_acd_0 + g_yzzzzz_0_xyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xyz_0[i] = g_yzzzzz_0_yz_1[i] * fi_acd_0 + g_yzzzzz_0_xyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xzz_0[i] = g_xzzzzz_0_xzz_1[i] * wa_y[i];

        g_xyzzzzz_0_yyy_0[i] = g_yzzzzz_0_yyy_1[i] * wa_x[i];

        g_xyzzzzz_0_yyz_0[i] = g_yzzzzz_0_yyz_1[i] * wa_x[i];

        g_xyzzzzz_0_yzz_0[i] = g_yzzzzz_0_yzz_1[i] * wa_x[i];

        g_xyzzzzz_0_zzz_0[i] = g_yzzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 270-280 components of targeted buffer : KSF

    auto g_xzzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 270);

    auto g_xzzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 271);

    auto g_xzzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 272);

    auto g_xzzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 273);

    auto g_xzzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 274);

    auto g_xzzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 275);

    auto g_xzzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 276);

    auto g_xzzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 277);

    auto g_xzzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 278);

    auto g_xzzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 279);

    #pragma omp simd aligned(g_xzzzzzz_0_xxx_0, g_xzzzzzz_0_xxy_0, g_xzzzzzz_0_xxz_0, g_xzzzzzz_0_xyy_0, g_xzzzzzz_0_xyz_0, g_xzzzzzz_0_xzz_0, g_xzzzzzz_0_yyy_0, g_xzzzzzz_0_yyz_0, g_xzzzzzz_0_yzz_0, g_xzzzzzz_0_zzz_0, g_zzzzzz_0_xx_1, g_zzzzzz_0_xxx_1, g_zzzzzz_0_xxy_1, g_zzzzzz_0_xxz_1, g_zzzzzz_0_xy_1, g_zzzzzz_0_xyy_1, g_zzzzzz_0_xyz_1, g_zzzzzz_0_xz_1, g_zzzzzz_0_xzz_1, g_zzzzzz_0_yy_1, g_zzzzzz_0_yyy_1, g_zzzzzz_0_yyz_1, g_zzzzzz_0_yz_1, g_zzzzzz_0_yzz_1, g_zzzzzz_0_zz_1, g_zzzzzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_xxx_0[i] = 3.0 * g_zzzzzz_0_xx_1[i] * fi_acd_0 + g_zzzzzz_0_xxx_1[i] * wa_x[i];

        g_xzzzzzz_0_xxy_0[i] = 2.0 * g_zzzzzz_0_xy_1[i] * fi_acd_0 + g_zzzzzz_0_xxy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxz_0[i] = 2.0 * g_zzzzzz_0_xz_1[i] * fi_acd_0 + g_zzzzzz_0_xxz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyy_0[i] = g_zzzzzz_0_yy_1[i] * fi_acd_0 + g_zzzzzz_0_xyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xyz_0[i] = g_zzzzzz_0_yz_1[i] * fi_acd_0 + g_zzzzzz_0_xyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xzz_0[i] = g_zzzzzz_0_zz_1[i] * fi_acd_0 + g_zzzzzz_0_xzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyy_0[i] = g_zzzzzz_0_yyy_1[i] * wa_x[i];

        g_xzzzzzz_0_yyz_0[i] = g_zzzzzz_0_yyz_1[i] * wa_x[i];

        g_xzzzzzz_0_yzz_0[i] = g_zzzzzz_0_yzz_1[i] * wa_x[i];

        g_xzzzzzz_0_zzz_0[i] = g_zzzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 280-290 components of targeted buffer : KSF

    auto g_yyyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 280);

    auto g_yyyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 281);

    auto g_yyyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 282);

    auto g_yyyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 283);

    auto g_yyyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 284);

    auto g_yyyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 285);

    auto g_yyyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 286);

    auto g_yyyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 287);

    auto g_yyyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 288);

    auto g_yyyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 289);

    #pragma omp simd aligned(g_yyyyy_0_xxx_0, g_yyyyy_0_xxx_1, g_yyyyy_0_xxy_0, g_yyyyy_0_xxy_1, g_yyyyy_0_xxz_0, g_yyyyy_0_xxz_1, g_yyyyy_0_xyy_0, g_yyyyy_0_xyy_1, g_yyyyy_0_xyz_0, g_yyyyy_0_xyz_1, g_yyyyy_0_xzz_0, g_yyyyy_0_xzz_1, g_yyyyy_0_yyy_0, g_yyyyy_0_yyy_1, g_yyyyy_0_yyz_0, g_yyyyy_0_yyz_1, g_yyyyy_0_yzz_0, g_yyyyy_0_yzz_1, g_yyyyy_0_zzz_0, g_yyyyy_0_zzz_1, g_yyyyyy_0_xx_1, g_yyyyyy_0_xxx_1, g_yyyyyy_0_xxy_1, g_yyyyyy_0_xxz_1, g_yyyyyy_0_xy_1, g_yyyyyy_0_xyy_1, g_yyyyyy_0_xyz_1, g_yyyyyy_0_xz_1, g_yyyyyy_0_xzz_1, g_yyyyyy_0_yy_1, g_yyyyyy_0_yyy_1, g_yyyyyy_0_yyz_1, g_yyyyyy_0_yz_1, g_yyyyyy_0_yzz_1, g_yyyyyy_0_zz_1, g_yyyyyy_0_zzz_1, g_yyyyyyy_0_xxx_0, g_yyyyyyy_0_xxy_0, g_yyyyyyy_0_xxz_0, g_yyyyyyy_0_xyy_0, g_yyyyyyy_0_xyz_0, g_yyyyyyy_0_xzz_0, g_yyyyyyy_0_yyy_0, g_yyyyyyy_0_yyz_0, g_yyyyyyy_0_yzz_0, g_yyyyyyy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_xxx_0[i] = 6.0 * g_yyyyy_0_xxx_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxx_1[i] * fz_be_0 + g_yyyyyy_0_xxx_1[i] * wa_y[i];

        g_yyyyyyy_0_xxy_0[i] = 6.0 * g_yyyyy_0_xxy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxy_1[i] * fz_be_0 + g_yyyyyy_0_xx_1[i] * fi_acd_0 + g_yyyyyy_0_xxy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxz_0[i] = 6.0 * g_yyyyy_0_xxz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxz_1[i] * fz_be_0 + g_yyyyyy_0_xxz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyy_0[i] = 6.0 * g_yyyyy_0_xyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xy_1[i] * fi_acd_0 + g_yyyyyy_0_xyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xyz_0[i] = 6.0 * g_yyyyy_0_xyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyz_1[i] * fz_be_0 + g_yyyyyy_0_xz_1[i] * fi_acd_0 + g_yyyyyy_0_xyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xzz_0[i] = 6.0 * g_yyyyy_0_xzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xzz_1[i] * fz_be_0 + g_yyyyyy_0_xzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyy_0[i] = 6.0 * g_yyyyy_0_yyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_yy_1[i] * fi_acd_0 + g_yyyyyy_0_yyy_1[i] * wa_y[i];

        g_yyyyyyy_0_yyz_0[i] = 6.0 * g_yyyyy_0_yyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_yz_1[i] * fi_acd_0 + g_yyyyyy_0_yyz_1[i] * wa_y[i];

        g_yyyyyyy_0_yzz_0[i] = 6.0 * g_yyyyy_0_yzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yzz_1[i] * fz_be_0 + g_yyyyyy_0_zz_1[i] * fi_acd_0 + g_yyyyyy_0_yzz_1[i] * wa_y[i];

        g_yyyyyyy_0_zzz_0[i] = 6.0 * g_yyyyy_0_zzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_zzz_1[i] * fz_be_0 + g_yyyyyy_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 290-300 components of targeted buffer : KSF

    auto g_yyyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 290);

    auto g_yyyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 291);

    auto g_yyyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 292);

    auto g_yyyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 293);

    auto g_yyyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 294);

    auto g_yyyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 295);

    auto g_yyyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 296);

    auto g_yyyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 297);

    auto g_yyyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 298);

    auto g_yyyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 299);

    #pragma omp simd aligned(g_yyyyyy_0_xx_1, g_yyyyyy_0_xxx_1, g_yyyyyy_0_xxy_1, g_yyyyyy_0_xxz_1, g_yyyyyy_0_xy_1, g_yyyyyy_0_xyy_1, g_yyyyyy_0_xyz_1, g_yyyyyy_0_xz_1, g_yyyyyy_0_xzz_1, g_yyyyyy_0_yy_1, g_yyyyyy_0_yyy_1, g_yyyyyy_0_yyz_1, g_yyyyyy_0_yz_1, g_yyyyyy_0_yzz_1, g_yyyyyy_0_zz_1, g_yyyyyy_0_zzz_1, g_yyyyyyz_0_xxx_0, g_yyyyyyz_0_xxy_0, g_yyyyyyz_0_xxz_0, g_yyyyyyz_0_xyy_0, g_yyyyyyz_0_xyz_0, g_yyyyyyz_0_xzz_0, g_yyyyyyz_0_yyy_0, g_yyyyyyz_0_yyz_0, g_yyyyyyz_0_yzz_0, g_yyyyyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_xxx_0[i] = g_yyyyyy_0_xxx_1[i] * wa_z[i];

        g_yyyyyyz_0_xxy_0[i] = g_yyyyyy_0_xxy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxz_0[i] = g_yyyyyy_0_xx_1[i] * fi_acd_0 + g_yyyyyy_0_xxz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyy_0[i] = g_yyyyyy_0_xyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xyz_0[i] = g_yyyyyy_0_xy_1[i] * fi_acd_0 + g_yyyyyy_0_xyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xzz_0[i] = 2.0 * g_yyyyyy_0_xz_1[i] * fi_acd_0 + g_yyyyyy_0_xzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyy_0[i] = g_yyyyyy_0_yyy_1[i] * wa_z[i];

        g_yyyyyyz_0_yyz_0[i] = g_yyyyyy_0_yy_1[i] * fi_acd_0 + g_yyyyyy_0_yyz_1[i] * wa_z[i];

        g_yyyyyyz_0_yzz_0[i] = 2.0 * g_yyyyyy_0_yz_1[i] * fi_acd_0 + g_yyyyyy_0_yzz_1[i] * wa_z[i];

        g_yyyyyyz_0_zzz_0[i] = 3.0 * g_yyyyyy_0_zz_1[i] * fi_acd_0 + g_yyyyyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 300-310 components of targeted buffer : KSF

    auto g_yyyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 300);

    auto g_yyyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 301);

    auto g_yyyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 302);

    auto g_yyyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 303);

    auto g_yyyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 304);

    auto g_yyyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 305);

    auto g_yyyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 306);

    auto g_yyyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 307);

    auto g_yyyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 308);

    auto g_yyyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 309);

    #pragma omp simd aligned(g_yyyyy_0_xxy_0, g_yyyyy_0_xxy_1, g_yyyyy_0_xyy_0, g_yyyyy_0_xyy_1, g_yyyyy_0_yyy_0, g_yyyyy_0_yyy_1, g_yyyyyz_0_xxy_1, g_yyyyyz_0_xyy_1, g_yyyyyz_0_yyy_1, g_yyyyyzz_0_xxx_0, g_yyyyyzz_0_xxy_0, g_yyyyyzz_0_xxz_0, g_yyyyyzz_0_xyy_0, g_yyyyyzz_0_xyz_0, g_yyyyyzz_0_xzz_0, g_yyyyyzz_0_yyy_0, g_yyyyyzz_0_yyz_0, g_yyyyyzz_0_yzz_0, g_yyyyyzz_0_zzz_0, g_yyyyzz_0_xxx_1, g_yyyyzz_0_xxz_1, g_yyyyzz_0_xyz_1, g_yyyyzz_0_xz_1, g_yyyyzz_0_xzz_1, g_yyyyzz_0_yyz_1, g_yyyyzz_0_yz_1, g_yyyyzz_0_yzz_1, g_yyyyzz_0_zz_1, g_yyyyzz_0_zzz_1, g_yyyzz_0_xxx_0, g_yyyzz_0_xxx_1, g_yyyzz_0_xxz_0, g_yyyzz_0_xxz_1, g_yyyzz_0_xyz_0, g_yyyzz_0_xyz_1, g_yyyzz_0_xzz_0, g_yyyzz_0_xzz_1, g_yyyzz_0_yyz_0, g_yyyzz_0_yyz_1, g_yyyzz_0_yzz_0, g_yyyzz_0_yzz_1, g_yyyzz_0_zzz_0, g_yyyzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyzz_0_xxx_0[i] = 4.0 * g_yyyzz_0_xxx_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxx_1[i] * fz_be_0 + g_yyyyzz_0_xxx_1[i] * wa_y[i];

        g_yyyyyzz_0_xxy_0[i] = g_yyyyy_0_xxy_0[i] * fbe_0 - g_yyyyy_0_xxy_1[i] * fz_be_0 + g_yyyyyz_0_xxy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxz_0[i] = 4.0 * g_yyyzz_0_xxz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxz_1[i] * fz_be_0 + g_yyyyzz_0_xxz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyy_0[i] = g_yyyyy_0_xyy_0[i] * fbe_0 - g_yyyyy_0_xyy_1[i] * fz_be_0 + g_yyyyyz_0_xyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xyz_0[i] = 4.0 * g_yyyzz_0_xyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyz_1[i] * fz_be_0 + g_yyyyzz_0_xz_1[i] * fi_acd_0 + g_yyyyzz_0_xyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xzz_0[i] = 4.0 * g_yyyzz_0_xzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xzz_1[i] * fz_be_0 + g_yyyyzz_0_xzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyy_0[i] = g_yyyyy_0_yyy_0[i] * fbe_0 - g_yyyyy_0_yyy_1[i] * fz_be_0 + g_yyyyyz_0_yyy_1[i] * wa_z[i];

        g_yyyyyzz_0_yyz_0[i] = 4.0 * g_yyyzz_0_yyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_yz_1[i] * fi_acd_0 + g_yyyyzz_0_yyz_1[i] * wa_y[i];

        g_yyyyyzz_0_yzz_0[i] = 4.0 * g_yyyzz_0_yzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yzz_1[i] * fz_be_0 + g_yyyyzz_0_zz_1[i] * fi_acd_0 + g_yyyyzz_0_yzz_1[i] * wa_y[i];

        g_yyyyyzz_0_zzz_0[i] = 4.0 * g_yyyzz_0_zzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_zzz_1[i] * fz_be_0 + g_yyyyzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 310-320 components of targeted buffer : KSF

    auto g_yyyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 310);

    auto g_yyyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 311);

    auto g_yyyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 312);

    auto g_yyyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 313);

    auto g_yyyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 314);

    auto g_yyyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 315);

    auto g_yyyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 316);

    auto g_yyyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 317);

    auto g_yyyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 318);

    auto g_yyyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 319);

    #pragma omp simd aligned(g_yyyyz_0_xxy_0, g_yyyyz_0_xxy_1, g_yyyyz_0_xyy_0, g_yyyyz_0_xyy_1, g_yyyyz_0_yyy_0, g_yyyyz_0_yyy_1, g_yyyyzz_0_xxy_1, g_yyyyzz_0_xyy_1, g_yyyyzz_0_yyy_1, g_yyyyzzz_0_xxx_0, g_yyyyzzz_0_xxy_0, g_yyyyzzz_0_xxz_0, g_yyyyzzz_0_xyy_0, g_yyyyzzz_0_xyz_0, g_yyyyzzz_0_xzz_0, g_yyyyzzz_0_yyy_0, g_yyyyzzz_0_yyz_0, g_yyyyzzz_0_yzz_0, g_yyyyzzz_0_zzz_0, g_yyyzzz_0_xxx_1, g_yyyzzz_0_xxz_1, g_yyyzzz_0_xyz_1, g_yyyzzz_0_xz_1, g_yyyzzz_0_xzz_1, g_yyyzzz_0_yyz_1, g_yyyzzz_0_yz_1, g_yyyzzz_0_yzz_1, g_yyyzzz_0_zz_1, g_yyyzzz_0_zzz_1, g_yyzzz_0_xxx_0, g_yyzzz_0_xxx_1, g_yyzzz_0_xxz_0, g_yyzzz_0_xxz_1, g_yyzzz_0_xyz_0, g_yyzzz_0_xyz_1, g_yyzzz_0_xzz_0, g_yyzzz_0_xzz_1, g_yyzzz_0_yyz_0, g_yyzzz_0_yyz_1, g_yyzzz_0_yzz_0, g_yyzzz_0_yzz_1, g_yyzzz_0_zzz_0, g_yyzzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzzz_0_xxx_0[i] = 3.0 * g_yyzzz_0_xxx_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxx_1[i] * fz_be_0 + g_yyyzzz_0_xxx_1[i] * wa_y[i];

        g_yyyyzzz_0_xxy_0[i] = 2.0 * g_yyyyz_0_xxy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxy_1[i] * fz_be_0 + g_yyyyzz_0_xxy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxz_0[i] = 3.0 * g_yyzzz_0_xxz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxz_1[i] * fz_be_0 + g_yyyzzz_0_xxz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyy_0[i] = 2.0 * g_yyyyz_0_xyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xyy_1[i] * fz_be_0 + g_yyyyzz_0_xyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xyz_0[i] = 3.0 * g_yyzzz_0_xyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyz_1[i] * fz_be_0 + g_yyyzzz_0_xz_1[i] * fi_acd_0 + g_yyyzzz_0_xyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xzz_0[i] = 3.0 * g_yyzzz_0_xzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xzz_1[i] * fz_be_0 + g_yyyzzz_0_xzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyy_0[i] = 2.0 * g_yyyyz_0_yyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_yyy_1[i] * fz_be_0 + g_yyyyzz_0_yyy_1[i] * wa_z[i];

        g_yyyyzzz_0_yyz_0[i] = 3.0 * g_yyzzz_0_yyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_yz_1[i] * fi_acd_0 + g_yyyzzz_0_yyz_1[i] * wa_y[i];

        g_yyyyzzz_0_yzz_0[i] = 3.0 * g_yyzzz_0_yzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yzz_1[i] * fz_be_0 + g_yyyzzz_0_zz_1[i] * fi_acd_0 + g_yyyzzz_0_yzz_1[i] * wa_y[i];

        g_yyyyzzz_0_zzz_0[i] = 3.0 * g_yyzzz_0_zzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_zzz_1[i] * fz_be_0 + g_yyyzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 320-330 components of targeted buffer : KSF

    auto g_yyyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 320);

    auto g_yyyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 321);

    auto g_yyyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 322);

    auto g_yyyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 323);

    auto g_yyyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 324);

    auto g_yyyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 325);

    auto g_yyyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 326);

    auto g_yyyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 327);

    auto g_yyyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 328);

    auto g_yyyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 329);

    #pragma omp simd aligned(g_yyyzz_0_xxy_0, g_yyyzz_0_xxy_1, g_yyyzz_0_xyy_0, g_yyyzz_0_xyy_1, g_yyyzz_0_yyy_0, g_yyyzz_0_yyy_1, g_yyyzzz_0_xxy_1, g_yyyzzz_0_xyy_1, g_yyyzzz_0_yyy_1, g_yyyzzzz_0_xxx_0, g_yyyzzzz_0_xxy_0, g_yyyzzzz_0_xxz_0, g_yyyzzzz_0_xyy_0, g_yyyzzzz_0_xyz_0, g_yyyzzzz_0_xzz_0, g_yyyzzzz_0_yyy_0, g_yyyzzzz_0_yyz_0, g_yyyzzzz_0_yzz_0, g_yyyzzzz_0_zzz_0, g_yyzzzz_0_xxx_1, g_yyzzzz_0_xxz_1, g_yyzzzz_0_xyz_1, g_yyzzzz_0_xz_1, g_yyzzzz_0_xzz_1, g_yyzzzz_0_yyz_1, g_yyzzzz_0_yz_1, g_yyzzzz_0_yzz_1, g_yyzzzz_0_zz_1, g_yyzzzz_0_zzz_1, g_yzzzz_0_xxx_0, g_yzzzz_0_xxx_1, g_yzzzz_0_xxz_0, g_yzzzz_0_xxz_1, g_yzzzz_0_xyz_0, g_yzzzz_0_xyz_1, g_yzzzz_0_xzz_0, g_yzzzz_0_xzz_1, g_yzzzz_0_yyz_0, g_yzzzz_0_yyz_1, g_yzzzz_0_yzz_0, g_yzzzz_0_yzz_1, g_yzzzz_0_zzz_0, g_yzzzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzzz_0_xxx_0[i] = 2.0 * g_yzzzz_0_xxx_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxx_1[i] * fz_be_0 + g_yyzzzz_0_xxx_1[i] * wa_y[i];

        g_yyyzzzz_0_xxy_0[i] = 3.0 * g_yyyzz_0_xxy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxy_1[i] * fz_be_0 + g_yyyzzz_0_xxy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxz_0[i] = 2.0 * g_yzzzz_0_xxz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxz_1[i] * fz_be_0 + g_yyzzzz_0_xxz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyy_0[i] = 3.0 * g_yyyzz_0_xyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xyy_1[i] * fz_be_0 + g_yyyzzz_0_xyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xyz_0[i] = 2.0 * g_yzzzz_0_xyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyz_1[i] * fz_be_0 + g_yyzzzz_0_xz_1[i] * fi_acd_0 + g_yyzzzz_0_xyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xzz_0[i] = 2.0 * g_yzzzz_0_xzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xzz_1[i] * fz_be_0 + g_yyzzzz_0_xzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyy_0[i] = 3.0 * g_yyyzz_0_yyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_yyy_1[i] * fz_be_0 + g_yyyzzz_0_yyy_1[i] * wa_z[i];

        g_yyyzzzz_0_yyz_0[i] = 2.0 * g_yzzzz_0_yyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_yz_1[i] * fi_acd_0 + g_yyzzzz_0_yyz_1[i] * wa_y[i];

        g_yyyzzzz_0_yzz_0[i] = 2.0 * g_yzzzz_0_yzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yzz_1[i] * fz_be_0 + g_yyzzzz_0_zz_1[i] * fi_acd_0 + g_yyzzzz_0_yzz_1[i] * wa_y[i];

        g_yyyzzzz_0_zzz_0[i] = 2.0 * g_yzzzz_0_zzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_zzz_1[i] * fz_be_0 + g_yyzzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 330-340 components of targeted buffer : KSF

    auto g_yyzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 330);

    auto g_yyzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 331);

    auto g_yyzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 332);

    auto g_yyzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 333);

    auto g_yyzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 334);

    auto g_yyzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 335);

    auto g_yyzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 336);

    auto g_yyzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 337);

    auto g_yyzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 338);

    auto g_yyzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 339);

    #pragma omp simd aligned(g_yyzzz_0_xxy_0, g_yyzzz_0_xxy_1, g_yyzzz_0_xyy_0, g_yyzzz_0_xyy_1, g_yyzzz_0_yyy_0, g_yyzzz_0_yyy_1, g_yyzzzz_0_xxy_1, g_yyzzzz_0_xyy_1, g_yyzzzz_0_yyy_1, g_yyzzzzz_0_xxx_0, g_yyzzzzz_0_xxy_0, g_yyzzzzz_0_xxz_0, g_yyzzzzz_0_xyy_0, g_yyzzzzz_0_xyz_0, g_yyzzzzz_0_xzz_0, g_yyzzzzz_0_yyy_0, g_yyzzzzz_0_yyz_0, g_yyzzzzz_0_yzz_0, g_yyzzzzz_0_zzz_0, g_yzzzzz_0_xxx_1, g_yzzzzz_0_xxz_1, g_yzzzzz_0_xyz_1, g_yzzzzz_0_xz_1, g_yzzzzz_0_xzz_1, g_yzzzzz_0_yyz_1, g_yzzzzz_0_yz_1, g_yzzzzz_0_yzz_1, g_yzzzzz_0_zz_1, g_yzzzzz_0_zzz_1, g_zzzzz_0_xxx_0, g_zzzzz_0_xxx_1, g_zzzzz_0_xxz_0, g_zzzzz_0_xxz_1, g_zzzzz_0_xyz_0, g_zzzzz_0_xyz_1, g_zzzzz_0_xzz_0, g_zzzzz_0_xzz_1, g_zzzzz_0_yyz_0, g_zzzzz_0_yyz_1, g_zzzzz_0_yzz_0, g_zzzzz_0_yzz_1, g_zzzzz_0_zzz_0, g_zzzzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzzz_0_xxx_0[i] = g_zzzzz_0_xxx_0[i] * fbe_0 - g_zzzzz_0_xxx_1[i] * fz_be_0 + g_yzzzzz_0_xxx_1[i] * wa_y[i];

        g_yyzzzzz_0_xxy_0[i] = 4.0 * g_yyzzz_0_xxy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxy_1[i] * fz_be_0 + g_yyzzzz_0_xxy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxz_0[i] = g_zzzzz_0_xxz_0[i] * fbe_0 - g_zzzzz_0_xxz_1[i] * fz_be_0 + g_yzzzzz_0_xxz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyy_0[i] = 4.0 * g_yyzzz_0_xyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xyy_1[i] * fz_be_0 + g_yyzzzz_0_xyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xyz_0[i] = g_zzzzz_0_xyz_0[i] * fbe_0 - g_zzzzz_0_xyz_1[i] * fz_be_0 + g_yzzzzz_0_xz_1[i] * fi_acd_0 + g_yzzzzz_0_xyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xzz_0[i] = g_zzzzz_0_xzz_0[i] * fbe_0 - g_zzzzz_0_xzz_1[i] * fz_be_0 + g_yzzzzz_0_xzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyy_0[i] = 4.0 * g_yyzzz_0_yyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_yyy_1[i] * fz_be_0 + g_yyzzzz_0_yyy_1[i] * wa_z[i];

        g_yyzzzzz_0_yyz_0[i] = g_zzzzz_0_yyz_0[i] * fbe_0 - g_zzzzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_yz_1[i] * fi_acd_0 + g_yzzzzz_0_yyz_1[i] * wa_y[i];

        g_yyzzzzz_0_yzz_0[i] = g_zzzzz_0_yzz_0[i] * fbe_0 - g_zzzzz_0_yzz_1[i] * fz_be_0 + g_yzzzzz_0_zz_1[i] * fi_acd_0 + g_yzzzzz_0_yzz_1[i] * wa_y[i];

        g_yyzzzzz_0_zzz_0[i] = g_zzzzz_0_zzz_0[i] * fbe_0 - g_zzzzz_0_zzz_1[i] * fz_be_0 + g_yzzzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 340-350 components of targeted buffer : KSF

    auto g_yzzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 340);

    auto g_yzzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 341);

    auto g_yzzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 342);

    auto g_yzzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 343);

    auto g_yzzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 344);

    auto g_yzzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 345);

    auto g_yzzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 346);

    auto g_yzzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 347);

    auto g_yzzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 348);

    auto g_yzzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 349);

    #pragma omp simd aligned(g_yzzzzzz_0_xxx_0, g_yzzzzzz_0_xxy_0, g_yzzzzzz_0_xxz_0, g_yzzzzzz_0_xyy_0, g_yzzzzzz_0_xyz_0, g_yzzzzzz_0_xzz_0, g_yzzzzzz_0_yyy_0, g_yzzzzzz_0_yyz_0, g_yzzzzzz_0_yzz_0, g_yzzzzzz_0_zzz_0, g_zzzzzz_0_xx_1, g_zzzzzz_0_xxx_1, g_zzzzzz_0_xxy_1, g_zzzzzz_0_xxz_1, g_zzzzzz_0_xy_1, g_zzzzzz_0_xyy_1, g_zzzzzz_0_xyz_1, g_zzzzzz_0_xz_1, g_zzzzzz_0_xzz_1, g_zzzzzz_0_yy_1, g_zzzzzz_0_yyy_1, g_zzzzzz_0_yyz_1, g_zzzzzz_0_yz_1, g_zzzzzz_0_yzz_1, g_zzzzzz_0_zz_1, g_zzzzzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_xxx_0[i] = g_zzzzzz_0_xxx_1[i] * wa_y[i];

        g_yzzzzzz_0_xxy_0[i] = g_zzzzzz_0_xx_1[i] * fi_acd_0 + g_zzzzzz_0_xxy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxz_0[i] = g_zzzzzz_0_xxz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyy_0[i] = 2.0 * g_zzzzzz_0_xy_1[i] * fi_acd_0 + g_zzzzzz_0_xyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xyz_0[i] = g_zzzzzz_0_xz_1[i] * fi_acd_0 + g_zzzzzz_0_xyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xzz_0[i] = g_zzzzzz_0_xzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyy_0[i] = 3.0 * g_zzzzzz_0_yy_1[i] * fi_acd_0 + g_zzzzzz_0_yyy_1[i] * wa_y[i];

        g_yzzzzzz_0_yyz_0[i] = 2.0 * g_zzzzzz_0_yz_1[i] * fi_acd_0 + g_zzzzzz_0_yyz_1[i] * wa_y[i];

        g_yzzzzzz_0_yzz_0[i] = g_zzzzzz_0_zz_1[i] * fi_acd_0 + g_zzzzzz_0_yzz_1[i] * wa_y[i];

        g_yzzzzzz_0_zzz_0[i] = g_zzzzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 350-360 components of targeted buffer : KSF

    auto g_zzzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_ksf + 350);

    auto g_zzzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_ksf + 351);

    auto g_zzzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_ksf + 352);

    auto g_zzzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_ksf + 353);

    auto g_zzzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_ksf + 354);

    auto g_zzzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_ksf + 355);

    auto g_zzzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_ksf + 356);

    auto g_zzzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_ksf + 357);

    auto g_zzzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_ksf + 358);

    auto g_zzzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_ksf + 359);

    #pragma omp simd aligned(g_zzzzz_0_xxx_0, g_zzzzz_0_xxx_1, g_zzzzz_0_xxy_0, g_zzzzz_0_xxy_1, g_zzzzz_0_xxz_0, g_zzzzz_0_xxz_1, g_zzzzz_0_xyy_0, g_zzzzz_0_xyy_1, g_zzzzz_0_xyz_0, g_zzzzz_0_xyz_1, g_zzzzz_0_xzz_0, g_zzzzz_0_xzz_1, g_zzzzz_0_yyy_0, g_zzzzz_0_yyy_1, g_zzzzz_0_yyz_0, g_zzzzz_0_yyz_1, g_zzzzz_0_yzz_0, g_zzzzz_0_yzz_1, g_zzzzz_0_zzz_0, g_zzzzz_0_zzz_1, g_zzzzzz_0_xx_1, g_zzzzzz_0_xxx_1, g_zzzzzz_0_xxy_1, g_zzzzzz_0_xxz_1, g_zzzzzz_0_xy_1, g_zzzzzz_0_xyy_1, g_zzzzzz_0_xyz_1, g_zzzzzz_0_xz_1, g_zzzzzz_0_xzz_1, g_zzzzzz_0_yy_1, g_zzzzzz_0_yyy_1, g_zzzzzz_0_yyz_1, g_zzzzzz_0_yz_1, g_zzzzzz_0_yzz_1, g_zzzzzz_0_zz_1, g_zzzzzz_0_zzz_1, g_zzzzzzz_0_xxx_0, g_zzzzzzz_0_xxy_0, g_zzzzzzz_0_xxz_0, g_zzzzzzz_0_xyy_0, g_zzzzzzz_0_xyz_0, g_zzzzzzz_0_xzz_0, g_zzzzzzz_0_yyy_0, g_zzzzzzz_0_yyz_0, g_zzzzzzz_0_yzz_0, g_zzzzzzz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_xxx_0[i] = 6.0 * g_zzzzz_0_xxx_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxx_1[i] * fz_be_0 + g_zzzzzz_0_xxx_1[i] * wa_z[i];

        g_zzzzzzz_0_xxy_0[i] = 6.0 * g_zzzzz_0_xxy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxy_1[i] * fz_be_0 + g_zzzzzz_0_xxy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxz_0[i] = 6.0 * g_zzzzz_0_xxz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxz_1[i] * fz_be_0 + g_zzzzzz_0_xx_1[i] * fi_acd_0 + g_zzzzzz_0_xxz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyy_0[i] = 6.0 * g_zzzzz_0_xyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyy_1[i] * fz_be_0 + g_zzzzzz_0_xyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xyz_0[i] = 6.0 * g_zzzzz_0_xyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyz_1[i] * fz_be_0 + g_zzzzzz_0_xy_1[i] * fi_acd_0 + g_zzzzzz_0_xyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xzz_0[i] = 6.0 * g_zzzzz_0_xzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xz_1[i] * fi_acd_0 + g_zzzzzz_0_xzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyy_0[i] = 6.0 * g_zzzzz_0_yyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyy_1[i] * fz_be_0 + g_zzzzzz_0_yyy_1[i] * wa_z[i];

        g_zzzzzzz_0_yyz_0[i] = 6.0 * g_zzzzz_0_yyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyz_1[i] * fz_be_0 + g_zzzzzz_0_yy_1[i] * fi_acd_0 + g_zzzzzz_0_yyz_1[i] * wa_z[i];

        g_zzzzzzz_0_yzz_0[i] = 6.0 * g_zzzzz_0_yzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_yz_1[i] * fi_acd_0 + g_zzzzzz_0_yzz_1[i] * wa_z[i];

        g_zzzzzzz_0_zzz_0[i] = 6.0 * g_zzzzz_0_zzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_zzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_zz_1[i] * fi_acd_0 + g_zzzzzz_0_zzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

