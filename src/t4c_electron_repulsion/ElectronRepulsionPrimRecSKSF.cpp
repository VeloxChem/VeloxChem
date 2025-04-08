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

#include "ElectronRepulsionPrimRecSKSF.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sksf(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sksf,
                                  size_t                idx_eri_0_shsf,
                                  size_t                idx_eri_1_shsf,
                                  size_t                idx_eri_1_sisd,
                                  size_t                idx_eri_0_sisf,
                                  size_t                idx_eri_1_sisf,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SHSF

    auto g_0_xxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_shsf);

    auto g_0_xxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 1);

    auto g_0_xxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 2);

    auto g_0_xxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 3);

    auto g_0_xxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 4);

    auto g_0_xxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 5);

    auto g_0_xxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 6);

    auto g_0_xxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 7);

    auto g_0_xxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 8);

    auto g_0_xxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 9);

    auto g_0_xxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 10);

    auto g_0_xxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 12);

    auto g_0_xxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 15);

    auto g_0_xxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 20);

    auto g_0_xxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 21);

    auto g_0_xxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 23);

    auto g_0_xxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 30);

    auto g_0_xxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 31);

    auto g_0_xxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 32);

    auto g_0_xxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 33);

    auto g_0_xxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 34);

    auto g_0_xxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 35);

    auto g_0_xxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 36);

    auto g_0_xxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 37);

    auto g_0_xxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 38);

    auto g_0_xxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 39);

    auto g_0_xxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 50);

    auto g_0_xxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 51);

    auto g_0_xxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 52);

    auto g_0_xxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 53);

    auto g_0_xxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 54);

    auto g_0_xxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 55);

    auto g_0_xxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 56);

    auto g_0_xxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 57);

    auto g_0_xxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 58);

    auto g_0_xxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 59);

    auto g_0_xxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 60);

    auto g_0_xxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 61);

    auto g_0_xxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 62);

    auto g_0_xxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 63);

    auto g_0_xxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 64);

    auto g_0_xxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 65);

    auto g_0_xxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 66);

    auto g_0_xxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 67);

    auto g_0_xxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 68);

    auto g_0_xxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 69);

    auto g_0_xxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 71);

    auto g_0_xxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 73);

    auto g_0_xxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 80);

    auto g_0_xxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 82);

    auto g_0_xxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 85);

    auto g_0_xxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 90);

    auto g_0_xxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 91);

    auto g_0_xxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 92);

    auto g_0_xxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 93);

    auto g_0_xxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 94);

    auto g_0_xxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 95);

    auto g_0_xxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 96);

    auto g_0_xxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 97);

    auto g_0_xxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 98);

    auto g_0_xxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 99);

    auto g_0_xyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 101);

    auto g_0_xyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 103);

    auto g_0_xyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 104);

    auto g_0_xyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 106);

    auto g_0_xyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 107);

    auto g_0_xyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 108);

    auto g_0_xyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 109);

    auto g_0_xyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 124);

    auto g_0_xyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 126);

    auto g_0_xyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 127);

    auto g_0_xyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 128);

    auto g_0_xyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 129);

    auto g_0_xzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 142);

    auto g_0_xzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 144);

    auto g_0_xzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 145);

    auto g_0_xzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 146);

    auto g_0_xzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 147);

    auto g_0_xzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 148);

    auto g_0_xzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 149);

    auto g_0_yyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 150);

    auto g_0_yyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 151);

    auto g_0_yyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 152);

    auto g_0_yyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 153);

    auto g_0_yyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 154);

    auto g_0_yyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 155);

    auto g_0_yyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 156);

    auto g_0_yyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 157);

    auto g_0_yyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 158);

    auto g_0_yyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 159);

    auto g_0_yyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 161);

    auto g_0_yyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 163);

    auto g_0_yyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 166);

    auto g_0_yyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 170);

    auto g_0_yyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 171);

    auto g_0_yyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 172);

    auto g_0_yyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 173);

    auto g_0_yyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 174);

    auto g_0_yyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 175);

    auto g_0_yyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 176);

    auto g_0_yyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 177);

    auto g_0_yyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 178);

    auto g_0_yyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 179);

    auto g_0_yyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 180);

    auto g_0_yyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 181);

    auto g_0_yyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 182);

    auto g_0_yyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 183);

    auto g_0_yyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 184);

    auto g_0_yyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 185);

    auto g_0_yyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 186);

    auto g_0_yyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 187);

    auto g_0_yyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 188);

    auto g_0_yyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 189);

    auto g_0_yzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 190);

    auto g_0_yzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 192);

    auto g_0_yzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 194);

    auto g_0_yzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 195);

    auto g_0_yzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 197);

    auto g_0_yzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 198);

    auto g_0_yzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 199);

    auto g_0_zzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 200);

    auto g_0_zzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 201);

    auto g_0_zzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 202);

    auto g_0_zzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 203);

    auto g_0_zzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 204);

    auto g_0_zzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 205);

    auto g_0_zzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 206);

    auto g_0_zzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 207);

    auto g_0_zzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 208);

    auto g_0_zzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 209);

    /// Set up components of auxilary buffer : SHSF

    auto g_0_xxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_shsf);

    auto g_0_xxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 1);

    auto g_0_xxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 2);

    auto g_0_xxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 3);

    auto g_0_xxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 4);

    auto g_0_xxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 5);

    auto g_0_xxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 6);

    auto g_0_xxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 7);

    auto g_0_xxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 8);

    auto g_0_xxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 9);

    auto g_0_xxxxy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 10);

    auto g_0_xxxxy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 12);

    auto g_0_xxxxy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 15);

    auto g_0_xxxxz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 20);

    auto g_0_xxxxz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 21);

    auto g_0_xxxxz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 23);

    auto g_0_xxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 30);

    auto g_0_xxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 31);

    auto g_0_xxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 32);

    auto g_0_xxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 33);

    auto g_0_xxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 34);

    auto g_0_xxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 35);

    auto g_0_xxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 36);

    auto g_0_xxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 37);

    auto g_0_xxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 38);

    auto g_0_xxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 39);

    auto g_0_xxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 50);

    auto g_0_xxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 51);

    auto g_0_xxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 52);

    auto g_0_xxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 53);

    auto g_0_xxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 54);

    auto g_0_xxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 55);

    auto g_0_xxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 56);

    auto g_0_xxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 57);

    auto g_0_xxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 58);

    auto g_0_xxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 59);

    auto g_0_xxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 60);

    auto g_0_xxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 61);

    auto g_0_xxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 62);

    auto g_0_xxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 63);

    auto g_0_xxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 64);

    auto g_0_xxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 65);

    auto g_0_xxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 66);

    auto g_0_xxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 67);

    auto g_0_xxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 68);

    auto g_0_xxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 69);

    auto g_0_xxyyz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 71);

    auto g_0_xxyyz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 73);

    auto g_0_xxyzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 80);

    auto g_0_xxyzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 82);

    auto g_0_xxyzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 85);

    auto g_0_xxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 90);

    auto g_0_xxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 91);

    auto g_0_xxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 92);

    auto g_0_xxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 93);

    auto g_0_xxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 94);

    auto g_0_xxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 95);

    auto g_0_xxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 96);

    auto g_0_xxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 97);

    auto g_0_xxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 98);

    auto g_0_xxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 99);

    auto g_0_xyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 101);

    auto g_0_xyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 103);

    auto g_0_xyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 104);

    auto g_0_xyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 106);

    auto g_0_xyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 107);

    auto g_0_xyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 108);

    auto g_0_xyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 109);

    auto g_0_xyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 124);

    auto g_0_xyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 126);

    auto g_0_xyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 127);

    auto g_0_xyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 128);

    auto g_0_xyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 129);

    auto g_0_xzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 142);

    auto g_0_xzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 144);

    auto g_0_xzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 145);

    auto g_0_xzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 146);

    auto g_0_xzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 147);

    auto g_0_xzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 148);

    auto g_0_xzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 149);

    auto g_0_yyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 150);

    auto g_0_yyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 151);

    auto g_0_yyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 152);

    auto g_0_yyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 153);

    auto g_0_yyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 154);

    auto g_0_yyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 155);

    auto g_0_yyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 156);

    auto g_0_yyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 157);

    auto g_0_yyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 158);

    auto g_0_yyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 159);

    auto g_0_yyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 161);

    auto g_0_yyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 163);

    auto g_0_yyyyz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 166);

    auto g_0_yyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 170);

    auto g_0_yyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 171);

    auto g_0_yyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 172);

    auto g_0_yyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 173);

    auto g_0_yyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 174);

    auto g_0_yyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 175);

    auto g_0_yyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 176);

    auto g_0_yyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 177);

    auto g_0_yyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 178);

    auto g_0_yyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 179);

    auto g_0_yyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 180);

    auto g_0_yyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 181);

    auto g_0_yyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 182);

    auto g_0_yyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 183);

    auto g_0_yyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 184);

    auto g_0_yyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 185);

    auto g_0_yyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 186);

    auto g_0_yyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 187);

    auto g_0_yyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 188);

    auto g_0_yyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 189);

    auto g_0_yzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 190);

    auto g_0_yzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 192);

    auto g_0_yzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 194);

    auto g_0_yzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 195);

    auto g_0_yzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 197);

    auto g_0_yzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 198);

    auto g_0_yzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 199);

    auto g_0_zzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 200);

    auto g_0_zzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 201);

    auto g_0_zzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 202);

    auto g_0_zzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 203);

    auto g_0_zzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 204);

    auto g_0_zzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 205);

    auto g_0_zzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 206);

    auto g_0_zzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 207);

    auto g_0_zzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 208);

    auto g_0_zzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 209);

    /// Set up components of auxilary buffer : SISD

    auto g_0_xxxxxx_0_xx_1 = pbuffer.data(idx_eri_1_sisd);

    auto g_0_xxxxxx_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 1);

    auto g_0_xxxxxx_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 2);

    auto g_0_xxxxxx_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 3);

    auto g_0_xxxxxx_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 4);

    auto g_0_xxxxxx_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 5);

    auto g_0_xxxxxz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 14);

    auto g_0_xxxxxz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 16);

    auto g_0_xxxxxz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 17);

    auto g_0_xxxxyy_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 18);

    auto g_0_xxxxyy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 19);

    auto g_0_xxxxyy_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 20);

    auto g_0_xxxxyy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 21);

    auto g_0_xxxxyy_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 22);

    auto g_0_xxxxyy_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 23);

    auto g_0_xxxxzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 30);

    auto g_0_xxxxzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 31);

    auto g_0_xxxxzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 32);

    auto g_0_xxxxzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 33);

    auto g_0_xxxxzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 34);

    auto g_0_xxxxzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 35);

    auto g_0_xxxyyy_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 36);

    auto g_0_xxxyyy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 37);

    auto g_0_xxxyyy_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 38);

    auto g_0_xxxyyy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 39);

    auto g_0_xxxyyy_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 40);

    auto g_0_xxxyyy_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 41);

    auto g_0_xxxzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 54);

    auto g_0_xxxzzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 55);

    auto g_0_xxxzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 56);

    auto g_0_xxxzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 57);

    auto g_0_xxxzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 58);

    auto g_0_xxxzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 59);

    auto g_0_xxyyyy_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 60);

    auto g_0_xxyyyy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 61);

    auto g_0_xxyyyy_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 62);

    auto g_0_xxyyyy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 63);

    auto g_0_xxyyyy_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 64);

    auto g_0_xxyyyy_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 65);

    auto g_0_xxyyzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 76);

    auto g_0_xxzzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 84);

    auto g_0_xxzzzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 85);

    auto g_0_xxzzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 86);

    auto g_0_xxzzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 87);

    auto g_0_xxzzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 88);

    auto g_0_xxzzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 89);

    auto g_0_xyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 91);

    auto g_0_xyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 93);

    auto g_0_xyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 94);

    auto g_0_xyyyzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 106);

    auto g_0_xyyzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 112);

    auto g_0_xzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 122);

    auto g_0_xzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 124);

    auto g_0_xzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 125);

    auto g_0_yyyyyy_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 126);

    auto g_0_yyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 127);

    auto g_0_yyyyyy_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 128);

    auto g_0_yyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 129);

    auto g_0_yyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 130);

    auto g_0_yyyyyy_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 131);

    auto g_0_yyyyyz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 134);

    auto g_0_yyyyyz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 136);

    auto g_0_yyyyyz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 137);

    auto g_0_yyyyzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 138);

    auto g_0_yyyyzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 139);

    auto g_0_yyyyzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 140);

    auto g_0_yyyyzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 141);

    auto g_0_yyyyzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 142);

    auto g_0_yyyyzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 143);

    auto g_0_yyyzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 144);

    auto g_0_yyyzzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 145);

    auto g_0_yyyzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 146);

    auto g_0_yyyzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 147);

    auto g_0_yyyzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 148);

    auto g_0_yyyzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 149);

    auto g_0_yyzzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 150);

    auto g_0_yyzzzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 151);

    auto g_0_yyzzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 152);

    auto g_0_yyzzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 153);

    auto g_0_yyzzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 154);

    auto g_0_yyzzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 155);

    auto g_0_yzzzzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 157);

    auto g_0_yzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 158);

    auto g_0_yzzzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 159);

    auto g_0_yzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 160);

    auto g_0_yzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 161);

    auto g_0_zzzzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 162);

    auto g_0_zzzzzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 163);

    auto g_0_zzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 164);

    auto g_0_zzzzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 165);

    auto g_0_zzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 166);

    auto g_0_zzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 167);

    /// Set up components of auxilary buffer : SISF

    auto g_0_xxxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_sisf);

    auto g_0_xxxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 1);

    auto g_0_xxxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 2);

    auto g_0_xxxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 3);

    auto g_0_xxxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 4);

    auto g_0_xxxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 5);

    auto g_0_xxxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 6);

    auto g_0_xxxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 7);

    auto g_0_xxxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 8);

    auto g_0_xxxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 9);

    auto g_0_xxxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 10);

    auto g_0_xxxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 11);

    auto g_0_xxxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 12);

    auto g_0_xxxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 13);

    auto g_0_xxxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 15);

    auto g_0_xxxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 16);

    auto g_0_xxxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 20);

    auto g_0_xxxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 21);

    auto g_0_xxxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 22);

    auto g_0_xxxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 23);

    auto g_0_xxxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 24);

    auto g_0_xxxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 25);

    auto g_0_xxxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 27);

    auto g_0_xxxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 28);

    auto g_0_xxxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 29);

    auto g_0_xxxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 30);

    auto g_0_xxxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 31);

    auto g_0_xxxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 32);

    auto g_0_xxxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 33);

    auto g_0_xxxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 34);

    auto g_0_xxxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 35);

    auto g_0_xxxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 36);

    auto g_0_xxxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 37);

    auto g_0_xxxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 38);

    auto g_0_xxxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 39);

    auto g_0_xxxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 50);

    auto g_0_xxxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 51);

    auto g_0_xxxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 52);

    auto g_0_xxxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 53);

    auto g_0_xxxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 54);

    auto g_0_xxxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 55);

    auto g_0_xxxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 56);

    auto g_0_xxxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 57);

    auto g_0_xxxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 58);

    auto g_0_xxxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 59);

    auto g_0_xxxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 60);

    auto g_0_xxxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 61);

    auto g_0_xxxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 62);

    auto g_0_xxxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 63);

    auto g_0_xxxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 64);

    auto g_0_xxxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 65);

    auto g_0_xxxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 66);

    auto g_0_xxxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 67);

    auto g_0_xxxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 68);

    auto g_0_xxxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 69);

    auto g_0_xxxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 71);

    auto g_0_xxxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 73);

    auto g_0_xxxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 80);

    auto g_0_xxxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 82);

    auto g_0_xxxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 85);

    auto g_0_xxxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 90);

    auto g_0_xxxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 91);

    auto g_0_xxxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 92);

    auto g_0_xxxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 93);

    auto g_0_xxxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 94);

    auto g_0_xxxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 95);

    auto g_0_xxxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 96);

    auto g_0_xxxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 97);

    auto g_0_xxxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 98);

    auto g_0_xxxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 99);

    auto g_0_xxyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 100);

    auto g_0_xxyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 101);

    auto g_0_xxyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 102);

    auto g_0_xxyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 103);

    auto g_0_xxyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 104);

    auto g_0_xxyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 105);

    auto g_0_xxyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 106);

    auto g_0_xxyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 107);

    auto g_0_xxyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 108);

    auto g_0_xxyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 109);

    auto g_0_xxyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 111);

    auto g_0_xxyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 113);

    auto g_0_xxyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 120);

    auto g_0_xxyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 121);

    auto g_0_xxyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 122);

    auto g_0_xxyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 123);

    auto g_0_xxyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 124);

    auto g_0_xxyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 125);

    auto g_0_xxyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 126);

    auto g_0_xxyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 127);

    auto g_0_xxyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 128);

    auto g_0_xxyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 129);

    auto g_0_xxyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 130);

    auto g_0_xxyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 132);

    auto g_0_xxyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 135);

    auto g_0_xxzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 140);

    auto g_0_xxzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 141);

    auto g_0_xxzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 142);

    auto g_0_xxzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 143);

    auto g_0_xxzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 144);

    auto g_0_xxzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 145);

    auto g_0_xxzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 146);

    auto g_0_xxzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 147);

    auto g_0_xxzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 148);

    auto g_0_xxzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 149);

    auto g_0_xyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 150);

    auto g_0_xyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 151);

    auto g_0_xyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 153);

    auto g_0_xyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 154);

    auto g_0_xyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 156);

    auto g_0_xyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 157);

    auto g_0_xyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 158);

    auto g_0_xyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 159);

    auto g_0_xyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 174);

    auto g_0_xyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 176);

    auto g_0_xyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 177);

    auto g_0_xyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 178);

    auto g_0_xyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 179);

    auto g_0_xyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 184);

    auto g_0_xyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 186);

    auto g_0_xyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 187);

    auto g_0_xyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 188);

    auto g_0_xyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 189);

    auto g_0_xzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 200);

    auto g_0_xzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 202);

    auto g_0_xzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 204);

    auto g_0_xzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 205);

    auto g_0_xzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 206);

    auto g_0_xzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 207);

    auto g_0_xzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 208);

    auto g_0_xzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 209);

    auto g_0_yyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 210);

    auto g_0_yyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 211);

    auto g_0_yyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 212);

    auto g_0_yyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 213);

    auto g_0_yyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 214);

    auto g_0_yyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 215);

    auto g_0_yyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 216);

    auto g_0_yyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 217);

    auto g_0_yyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 218);

    auto g_0_yyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 219);

    auto g_0_yyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 221);

    auto g_0_yyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 222);

    auto g_0_yyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 223);

    auto g_0_yyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 224);

    auto g_0_yyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 225);

    auto g_0_yyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 226);

    auto g_0_yyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 227);

    auto g_0_yyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 228);

    auto g_0_yyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 229);

    auto g_0_yyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 230);

    auto g_0_yyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 231);

    auto g_0_yyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 232);

    auto g_0_yyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 233);

    auto g_0_yyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 234);

    auto g_0_yyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 235);

    auto g_0_yyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 236);

    auto g_0_yyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 237);

    auto g_0_yyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 238);

    auto g_0_yyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 239);

    auto g_0_yyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 240);

    auto g_0_yyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 241);

    auto g_0_yyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 242);

    auto g_0_yyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 243);

    auto g_0_yyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 244);

    auto g_0_yyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 245);

    auto g_0_yyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 246);

    auto g_0_yyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 247);

    auto g_0_yyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 248);

    auto g_0_yyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 249);

    auto g_0_yyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 250);

    auto g_0_yyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 251);

    auto g_0_yyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 252);

    auto g_0_yyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 253);

    auto g_0_yyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 254);

    auto g_0_yyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 255);

    auto g_0_yyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 256);

    auto g_0_yyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 257);

    auto g_0_yyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 258);

    auto g_0_yyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 259);

    auto g_0_yzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 260);

    auto g_0_yzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 261);

    auto g_0_yzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 262);

    auto g_0_yzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 263);

    auto g_0_yzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 264);

    auto g_0_yzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 265);

    auto g_0_yzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 266);

    auto g_0_yzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 267);

    auto g_0_yzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 268);

    auto g_0_yzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 269);

    auto g_0_zzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 270);

    auto g_0_zzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 271);

    auto g_0_zzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 272);

    auto g_0_zzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 273);

    auto g_0_zzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 274);

    auto g_0_zzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 275);

    auto g_0_zzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 276);

    auto g_0_zzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 277);

    auto g_0_zzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 278);

    auto g_0_zzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 279);

    /// Set up components of auxilary buffer : SISF

    auto g_0_xxxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_sisf);

    auto g_0_xxxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 1);

    auto g_0_xxxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 2);

    auto g_0_xxxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 3);

    auto g_0_xxxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 4);

    auto g_0_xxxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 5);

    auto g_0_xxxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 6);

    auto g_0_xxxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 7);

    auto g_0_xxxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 8);

    auto g_0_xxxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 9);

    auto g_0_xxxxxy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 10);

    auto g_0_xxxxxy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 11);

    auto g_0_xxxxxy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 12);

    auto g_0_xxxxxy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 13);

    auto g_0_xxxxxy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 15);

    auto g_0_xxxxxy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 16);

    auto g_0_xxxxxz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 20);

    auto g_0_xxxxxz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 21);

    auto g_0_xxxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 22);

    auto g_0_xxxxxz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 23);

    auto g_0_xxxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 24);

    auto g_0_xxxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 25);

    auto g_0_xxxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 27);

    auto g_0_xxxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 28);

    auto g_0_xxxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 29);

    auto g_0_xxxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 30);

    auto g_0_xxxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 31);

    auto g_0_xxxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 32);

    auto g_0_xxxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 33);

    auto g_0_xxxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 34);

    auto g_0_xxxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 35);

    auto g_0_xxxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 36);

    auto g_0_xxxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 37);

    auto g_0_xxxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 38);

    auto g_0_xxxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 39);

    auto g_0_xxxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 50);

    auto g_0_xxxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 51);

    auto g_0_xxxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 52);

    auto g_0_xxxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 53);

    auto g_0_xxxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 54);

    auto g_0_xxxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 55);

    auto g_0_xxxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 56);

    auto g_0_xxxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 57);

    auto g_0_xxxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 58);

    auto g_0_xxxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 59);

    auto g_0_xxxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 60);

    auto g_0_xxxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 61);

    auto g_0_xxxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 62);

    auto g_0_xxxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 63);

    auto g_0_xxxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 64);

    auto g_0_xxxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 65);

    auto g_0_xxxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 66);

    auto g_0_xxxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 67);

    auto g_0_xxxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 68);

    auto g_0_xxxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 69);

    auto g_0_xxxyyz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 71);

    auto g_0_xxxyyz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 73);

    auto g_0_xxxyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 80);

    auto g_0_xxxyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 82);

    auto g_0_xxxyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 85);

    auto g_0_xxxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 90);

    auto g_0_xxxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 91);

    auto g_0_xxxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 92);

    auto g_0_xxxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 93);

    auto g_0_xxxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 94);

    auto g_0_xxxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 95);

    auto g_0_xxxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 96);

    auto g_0_xxxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 97);

    auto g_0_xxxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 98);

    auto g_0_xxxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 99);

    auto g_0_xxyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 100);

    auto g_0_xxyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 101);

    auto g_0_xxyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 102);

    auto g_0_xxyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 103);

    auto g_0_xxyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 104);

    auto g_0_xxyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 105);

    auto g_0_xxyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 106);

    auto g_0_xxyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 107);

    auto g_0_xxyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 108);

    auto g_0_xxyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 109);

    auto g_0_xxyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 111);

    auto g_0_xxyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 113);

    auto g_0_xxyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 120);

    auto g_0_xxyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 121);

    auto g_0_xxyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 122);

    auto g_0_xxyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 123);

    auto g_0_xxyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 124);

    auto g_0_xxyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 125);

    auto g_0_xxyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 126);

    auto g_0_xxyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 127);

    auto g_0_xxyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 128);

    auto g_0_xxyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 129);

    auto g_0_xxyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 130);

    auto g_0_xxyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 132);

    auto g_0_xxyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 135);

    auto g_0_xxzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 140);

    auto g_0_xxzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 141);

    auto g_0_xxzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 142);

    auto g_0_xxzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 143);

    auto g_0_xxzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 144);

    auto g_0_xxzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 145);

    auto g_0_xxzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 146);

    auto g_0_xxzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 147);

    auto g_0_xxzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 148);

    auto g_0_xxzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 149);

    auto g_0_xyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 150);

    auto g_0_xyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 151);

    auto g_0_xyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 153);

    auto g_0_xyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 154);

    auto g_0_xyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 156);

    auto g_0_xyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 157);

    auto g_0_xyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 158);

    auto g_0_xyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 159);

    auto g_0_xyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 174);

    auto g_0_xyyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 176);

    auto g_0_xyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 177);

    auto g_0_xyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 178);

    auto g_0_xyyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 179);

    auto g_0_xyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 184);

    auto g_0_xyyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 186);

    auto g_0_xyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 187);

    auto g_0_xyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 188);

    auto g_0_xyyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 189);

    auto g_0_xzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 200);

    auto g_0_xzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 202);

    auto g_0_xzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 204);

    auto g_0_xzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 205);

    auto g_0_xzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 206);

    auto g_0_xzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 207);

    auto g_0_xzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 208);

    auto g_0_xzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 209);

    auto g_0_yyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 210);

    auto g_0_yyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 211);

    auto g_0_yyyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 212);

    auto g_0_yyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 213);

    auto g_0_yyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 214);

    auto g_0_yyyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 215);

    auto g_0_yyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 216);

    auto g_0_yyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 217);

    auto g_0_yyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 218);

    auto g_0_yyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 219);

    auto g_0_yyyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 221);

    auto g_0_yyyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 222);

    auto g_0_yyyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 223);

    auto g_0_yyyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 224);

    auto g_0_yyyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 225);

    auto g_0_yyyyyz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 226);

    auto g_0_yyyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 227);

    auto g_0_yyyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 228);

    auto g_0_yyyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 229);

    auto g_0_yyyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 230);

    auto g_0_yyyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 231);

    auto g_0_yyyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 232);

    auto g_0_yyyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 233);

    auto g_0_yyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 234);

    auto g_0_yyyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 235);

    auto g_0_yyyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 236);

    auto g_0_yyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 237);

    auto g_0_yyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 238);

    auto g_0_yyyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 239);

    auto g_0_yyyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 240);

    auto g_0_yyyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 241);

    auto g_0_yyyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 242);

    auto g_0_yyyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 243);

    auto g_0_yyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 244);

    auto g_0_yyyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 245);

    auto g_0_yyyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 246);

    auto g_0_yyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 247);

    auto g_0_yyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 248);

    auto g_0_yyyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 249);

    auto g_0_yyzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 250);

    auto g_0_yyzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 251);

    auto g_0_yyzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 252);

    auto g_0_yyzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 253);

    auto g_0_yyzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 254);

    auto g_0_yyzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 255);

    auto g_0_yyzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 256);

    auto g_0_yyzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 257);

    auto g_0_yyzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 258);

    auto g_0_yyzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 259);

    auto g_0_yzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 260);

    auto g_0_yzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 261);

    auto g_0_yzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 262);

    auto g_0_yzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 263);

    auto g_0_yzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 264);

    auto g_0_yzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 265);

    auto g_0_yzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 266);

    auto g_0_yzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 267);

    auto g_0_yzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 268);

    auto g_0_yzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 269);

    auto g_0_zzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sisf + 270);

    auto g_0_zzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sisf + 271);

    auto g_0_zzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sisf + 272);

    auto g_0_zzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sisf + 273);

    auto g_0_zzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sisf + 274);

    auto g_0_zzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sisf + 275);

    auto g_0_zzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sisf + 276);

    auto g_0_zzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sisf + 277);

    auto g_0_zzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sisf + 278);

    auto g_0_zzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sisf + 279);

    /// Set up 0-10 components of targeted buffer : SKSF

    auto g_0_xxxxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_sksf);

    auto g_0_xxxxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 1);

    auto g_0_xxxxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 2);

    auto g_0_xxxxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 3);

    auto g_0_xxxxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 4);

    auto g_0_xxxxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 5);

    auto g_0_xxxxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 6);

    auto g_0_xxxxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 7);

    auto g_0_xxxxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 8);

    auto g_0_xxxxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 9);

#pragma omp simd aligned(g_0_xxxxx_0_xxx_0,       \
                             g_0_xxxxx_0_xxx_1,   \
                             g_0_xxxxx_0_xxy_0,   \
                             g_0_xxxxx_0_xxy_1,   \
                             g_0_xxxxx_0_xxz_0,   \
                             g_0_xxxxx_0_xxz_1,   \
                             g_0_xxxxx_0_xyy_0,   \
                             g_0_xxxxx_0_xyy_1,   \
                             g_0_xxxxx_0_xyz_0,   \
                             g_0_xxxxx_0_xyz_1,   \
                             g_0_xxxxx_0_xzz_0,   \
                             g_0_xxxxx_0_xzz_1,   \
                             g_0_xxxxx_0_yyy_0,   \
                             g_0_xxxxx_0_yyy_1,   \
                             g_0_xxxxx_0_yyz_0,   \
                             g_0_xxxxx_0_yyz_1,   \
                             g_0_xxxxx_0_yzz_0,   \
                             g_0_xxxxx_0_yzz_1,   \
                             g_0_xxxxx_0_zzz_0,   \
                             g_0_xxxxx_0_zzz_1,   \
                             g_0_xxxxxx_0_xx_1,   \
                             g_0_xxxxxx_0_xxx_0,  \
                             g_0_xxxxxx_0_xxx_1,  \
                             g_0_xxxxxx_0_xxy_0,  \
                             g_0_xxxxxx_0_xxy_1,  \
                             g_0_xxxxxx_0_xxz_0,  \
                             g_0_xxxxxx_0_xxz_1,  \
                             g_0_xxxxxx_0_xy_1,   \
                             g_0_xxxxxx_0_xyy_0,  \
                             g_0_xxxxxx_0_xyy_1,  \
                             g_0_xxxxxx_0_xyz_0,  \
                             g_0_xxxxxx_0_xyz_1,  \
                             g_0_xxxxxx_0_xz_1,   \
                             g_0_xxxxxx_0_xzz_0,  \
                             g_0_xxxxxx_0_xzz_1,  \
                             g_0_xxxxxx_0_yy_1,   \
                             g_0_xxxxxx_0_yyy_0,  \
                             g_0_xxxxxx_0_yyy_1,  \
                             g_0_xxxxxx_0_yyz_0,  \
                             g_0_xxxxxx_0_yyz_1,  \
                             g_0_xxxxxx_0_yz_1,   \
                             g_0_xxxxxx_0_yzz_0,  \
                             g_0_xxxxxx_0_yzz_1,  \
                             g_0_xxxxxx_0_zz_1,   \
                             g_0_xxxxxx_0_zzz_0,  \
                             g_0_xxxxxx_0_zzz_1,  \
                             g_0_xxxxxxx_0_xxx_0, \
                             g_0_xxxxxxx_0_xxy_0, \
                             g_0_xxxxxxx_0_xxz_0, \
                             g_0_xxxxxxx_0_xyy_0, \
                             g_0_xxxxxxx_0_xyz_0, \
                             g_0_xxxxxxx_0_xzz_0, \
                             g_0_xxxxxxx_0_yyy_0, \
                             g_0_xxxxxxx_0_yyz_0, \
                             g_0_xxxxxxx_0_yzz_0, \
                             g_0_xxxxxxx_0_zzz_0, \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_xxx_0[i] = 6.0 * g_0_xxxxx_0_xxx_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxx_1[i] * fti_ab_0 +
                                 3.0 * g_0_xxxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxx_0[i] * pb_x + g_0_xxxxxx_0_xxx_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxy_0[i] = 6.0 * g_0_xxxxx_0_xxy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxy_0[i] * pb_x + g_0_xxxxxx_0_xxy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxz_0[i] = 6.0 * g_0_xxxxx_0_xxz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxz_0[i] * pb_x + g_0_xxxxxx_0_xxz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyy_0[i] = 6.0 * g_0_xxxxx_0_xyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yy_1[i] * fi_abcd_0 +
                                 g_0_xxxxxx_0_xyy_0[i] * pb_x + g_0_xxxxxx_0_xyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyz_0[i] = 6.0 * g_0_xxxxx_0_xyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxxxxx_0_xyz_0[i] * pb_x + g_0_xxxxxx_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xzz_0[i] = 6.0 * g_0_xxxxx_0_xzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zz_1[i] * fi_abcd_0 +
                                 g_0_xxxxxx_0_xzz_0[i] * pb_x + g_0_xxxxxx_0_xzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyy_0[i] = 6.0 * g_0_xxxxx_0_yyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyy_0[i] * pb_x +
                                 g_0_xxxxxx_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyz_0[i] = 6.0 * g_0_xxxxx_0_yyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyz_0[i] * pb_x +
                                 g_0_xxxxxx_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yzz_0[i] = 6.0 * g_0_xxxxx_0_yzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yzz_0[i] * pb_x +
                                 g_0_xxxxxx_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_zzz_0[i] = 6.0 * g_0_xxxxx_0_zzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_zzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zzz_0[i] * pb_x +
                                 g_0_xxxxxx_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : SKSF

    auto g_0_xxxxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 10);

    auto g_0_xxxxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 11);

    auto g_0_xxxxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 12);

    auto g_0_xxxxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 13);

    auto g_0_xxxxxxy_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 14);

    auto g_0_xxxxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 15);

    auto g_0_xxxxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 16);

    auto g_0_xxxxxxy_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 17);

    auto g_0_xxxxxxy_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 18);

    auto g_0_xxxxxxy_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 19);

#pragma omp simd aligned(g_0_xxxxxx_0_xx_1,       \
                             g_0_xxxxxx_0_xxx_0,  \
                             g_0_xxxxxx_0_xxx_1,  \
                             g_0_xxxxxx_0_xxy_0,  \
                             g_0_xxxxxx_0_xxy_1,  \
                             g_0_xxxxxx_0_xxz_0,  \
                             g_0_xxxxxx_0_xxz_1,  \
                             g_0_xxxxxx_0_xy_1,   \
                             g_0_xxxxxx_0_xyy_0,  \
                             g_0_xxxxxx_0_xyy_1,  \
                             g_0_xxxxxx_0_xyz_0,  \
                             g_0_xxxxxx_0_xyz_1,  \
                             g_0_xxxxxx_0_xz_1,   \
                             g_0_xxxxxx_0_xzz_0,  \
                             g_0_xxxxxx_0_xzz_1,  \
                             g_0_xxxxxx_0_yy_1,   \
                             g_0_xxxxxx_0_yyy_0,  \
                             g_0_xxxxxx_0_yyy_1,  \
                             g_0_xxxxxx_0_yyz_0,  \
                             g_0_xxxxxx_0_yyz_1,  \
                             g_0_xxxxxx_0_yz_1,   \
                             g_0_xxxxxx_0_yzz_0,  \
                             g_0_xxxxxx_0_yzz_1,  \
                             g_0_xxxxxx_0_zz_1,   \
                             g_0_xxxxxx_0_zzz_0,  \
                             g_0_xxxxxx_0_zzz_1,  \
                             g_0_xxxxxxy_0_xxx_0, \
                             g_0_xxxxxxy_0_xxy_0, \
                             g_0_xxxxxxy_0_xxz_0, \
                             g_0_xxxxxxy_0_xyy_0, \
                             g_0_xxxxxxy_0_xyz_0, \
                             g_0_xxxxxxy_0_xzz_0, \
                             g_0_xxxxxxy_0_yyy_0, \
                             g_0_xxxxxxy_0_yyz_0, \
                             g_0_xxxxxxy_0_yzz_0, \
                             g_0_xxxxxxy_0_zzz_0, \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxy_0_xxx_0[i] = g_0_xxxxxx_0_xxx_0[i] * pb_y + g_0_xxxxxx_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxy_0[i] = g_0_xxxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxy_0[i] * pb_y + g_0_xxxxxx_0_xxy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxz_0[i] = g_0_xxxxxx_0_xxz_0[i] * pb_y + g_0_xxxxxx_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyy_0[i] = 2.0 * g_0_xxxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyy_0[i] * pb_y + g_0_xxxxxx_0_xyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyz_0[i] = g_0_xxxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyz_0[i] * pb_y + g_0_xxxxxx_0_xyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xzz_0[i] = g_0_xxxxxx_0_xzz_0[i] * pb_y + g_0_xxxxxx_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyy_0[i] = 3.0 * g_0_xxxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyy_0[i] * pb_y + g_0_xxxxxx_0_yyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyz_0[i] = 2.0 * g_0_xxxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyz_0[i] * pb_y + g_0_xxxxxx_0_yyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yzz_0[i] = g_0_xxxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzz_0[i] * pb_y + g_0_xxxxxx_0_yzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_zzz_0[i] = g_0_xxxxxx_0_zzz_0[i] * pb_y + g_0_xxxxxx_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : SKSF

    auto g_0_xxxxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 20);

    auto g_0_xxxxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 21);

    auto g_0_xxxxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 22);

    auto g_0_xxxxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 23);

    auto g_0_xxxxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 24);

    auto g_0_xxxxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 25);

    auto g_0_xxxxxxz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 26);

    auto g_0_xxxxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 27);

    auto g_0_xxxxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 28);

    auto g_0_xxxxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 29);

#pragma omp simd aligned(g_0_xxxxxx_0_xx_1,       \
                             g_0_xxxxxx_0_xxx_0,  \
                             g_0_xxxxxx_0_xxx_1,  \
                             g_0_xxxxxx_0_xxy_0,  \
                             g_0_xxxxxx_0_xxy_1,  \
                             g_0_xxxxxx_0_xxz_0,  \
                             g_0_xxxxxx_0_xxz_1,  \
                             g_0_xxxxxx_0_xy_1,   \
                             g_0_xxxxxx_0_xyy_0,  \
                             g_0_xxxxxx_0_xyy_1,  \
                             g_0_xxxxxx_0_xyz_0,  \
                             g_0_xxxxxx_0_xyz_1,  \
                             g_0_xxxxxx_0_xz_1,   \
                             g_0_xxxxxx_0_xzz_0,  \
                             g_0_xxxxxx_0_xzz_1,  \
                             g_0_xxxxxx_0_yy_1,   \
                             g_0_xxxxxx_0_yyy_0,  \
                             g_0_xxxxxx_0_yyy_1,  \
                             g_0_xxxxxx_0_yyz_0,  \
                             g_0_xxxxxx_0_yyz_1,  \
                             g_0_xxxxxx_0_yz_1,   \
                             g_0_xxxxxx_0_yzz_0,  \
                             g_0_xxxxxx_0_yzz_1,  \
                             g_0_xxxxxx_0_zz_1,   \
                             g_0_xxxxxx_0_zzz_0,  \
                             g_0_xxxxxx_0_zzz_1,  \
                             g_0_xxxxxxz_0_xxx_0, \
                             g_0_xxxxxxz_0_xxy_0, \
                             g_0_xxxxxxz_0_xxz_0, \
                             g_0_xxxxxxz_0_xyy_0, \
                             g_0_xxxxxxz_0_xyz_0, \
                             g_0_xxxxxxz_0_xzz_0, \
                             g_0_xxxxxxz_0_yyy_0, \
                             g_0_xxxxxxz_0_yyz_0, \
                             g_0_xxxxxxz_0_yzz_0, \
                             g_0_xxxxxxz_0_zzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxz_0_xxx_0[i] = g_0_xxxxxx_0_xxx_0[i] * pb_z + g_0_xxxxxx_0_xxx_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxy_0[i] = g_0_xxxxxx_0_xxy_0[i] * pb_z + g_0_xxxxxx_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxz_0[i] = g_0_xxxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxz_0[i] * pb_z + g_0_xxxxxx_0_xxz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyy_0[i] = g_0_xxxxxx_0_xyy_0[i] * pb_z + g_0_xxxxxx_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyz_0[i] = g_0_xxxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyz_0[i] * pb_z + g_0_xxxxxx_0_xyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xzz_0[i] = 2.0 * g_0_xxxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzz_0[i] * pb_z + g_0_xxxxxx_0_xzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyy_0[i] = g_0_xxxxxx_0_yyy_0[i] * pb_z + g_0_xxxxxx_0_yyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyz_0[i] = g_0_xxxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyz_0[i] * pb_z + g_0_xxxxxx_0_yyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yzz_0[i] = 2.0 * g_0_xxxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzz_0[i] * pb_z + g_0_xxxxxx_0_yzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_zzz_0[i] = 3.0 * g_0_xxxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_zzz_0[i] * pb_z + g_0_xxxxxx_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 30-40 components of targeted buffer : SKSF

    auto g_0_xxxxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 30);

    auto g_0_xxxxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 31);

    auto g_0_xxxxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 32);

    auto g_0_xxxxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 33);

    auto g_0_xxxxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 34);

    auto g_0_xxxxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 35);

    auto g_0_xxxxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 36);

    auto g_0_xxxxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 37);

    auto g_0_xxxxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 38);

    auto g_0_xxxxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 39);

#pragma omp simd aligned(g_0_xxxxx_0_xxx_0,       \
                             g_0_xxxxx_0_xxx_1,   \
                             g_0_xxxxx_0_xxz_0,   \
                             g_0_xxxxx_0_xxz_1,   \
                             g_0_xxxxx_0_xzz_0,   \
                             g_0_xxxxx_0_xzz_1,   \
                             g_0_xxxxxy_0_xxx_0,  \
                             g_0_xxxxxy_0_xxx_1,  \
                             g_0_xxxxxy_0_xxz_0,  \
                             g_0_xxxxxy_0_xxz_1,  \
                             g_0_xxxxxy_0_xzz_0,  \
                             g_0_xxxxxy_0_xzz_1,  \
                             g_0_xxxxxyy_0_xxx_0, \
                             g_0_xxxxxyy_0_xxy_0, \
                             g_0_xxxxxyy_0_xxz_0, \
                             g_0_xxxxxyy_0_xyy_0, \
                             g_0_xxxxxyy_0_xyz_0, \
                             g_0_xxxxxyy_0_xzz_0, \
                             g_0_xxxxxyy_0_yyy_0, \
                             g_0_xxxxxyy_0_yyz_0, \
                             g_0_xxxxxyy_0_yzz_0, \
                             g_0_xxxxxyy_0_zzz_0, \
                             g_0_xxxxyy_0_xxy_0,  \
                             g_0_xxxxyy_0_xxy_1,  \
                             g_0_xxxxyy_0_xy_1,   \
                             g_0_xxxxyy_0_xyy_0,  \
                             g_0_xxxxyy_0_xyy_1,  \
                             g_0_xxxxyy_0_xyz_0,  \
                             g_0_xxxxyy_0_xyz_1,  \
                             g_0_xxxxyy_0_yy_1,   \
                             g_0_xxxxyy_0_yyy_0,  \
                             g_0_xxxxyy_0_yyy_1,  \
                             g_0_xxxxyy_0_yyz_0,  \
                             g_0_xxxxyy_0_yyz_1,  \
                             g_0_xxxxyy_0_yz_1,   \
                             g_0_xxxxyy_0_yzz_0,  \
                             g_0_xxxxyy_0_yzz_1,  \
                             g_0_xxxxyy_0_zzz_0,  \
                             g_0_xxxxyy_0_zzz_1,  \
                             g_0_xxxyy_0_xxy_0,   \
                             g_0_xxxyy_0_xxy_1,   \
                             g_0_xxxyy_0_xyy_0,   \
                             g_0_xxxyy_0_xyy_1,   \
                             g_0_xxxyy_0_xyz_0,   \
                             g_0_xxxyy_0_xyz_1,   \
                             g_0_xxxyy_0_yyy_0,   \
                             g_0_xxxyy_0_yyy_1,   \
                             g_0_xxxyy_0_yyz_0,   \
                             g_0_xxxyy_0_yyz_1,   \
                             g_0_xxxyy_0_yzz_0,   \
                             g_0_xxxyy_0_yzz_1,   \
                             g_0_xxxyy_0_zzz_0,   \
                             g_0_xxxyy_0_zzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyy_0_xxx_0[i] =
            g_0_xxxxx_0_xxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxx_0[i] * pb_y + g_0_xxxxxy_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxy_0[i] = 4.0 * g_0_xxxyy_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxxyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxy_0[i] * pb_x + g_0_xxxxyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxz_0[i] =
            g_0_xxxxx_0_xxz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxz_0[i] * pb_y + g_0_xxxxxy_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xyy_0[i] = 4.0 * g_0_xxxyy_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yy_1[i] * fi_abcd_0 +
                                 g_0_xxxxyy_0_xyy_0[i] * pb_x + g_0_xxxxyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyz_0[i] = 4.0 * g_0_xxxyy_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxxxyy_0_xyz_0[i] * pb_x + g_0_xxxxyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xzz_0[i] =
            g_0_xxxxx_0_xzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xzz_0[i] * pb_y + g_0_xxxxxy_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_yyy_0[i] = 4.0 * g_0_xxxyy_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyy_0[i] * pb_x +
                                 g_0_xxxxyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyz_0[i] = 4.0 * g_0_xxxyy_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyz_0[i] * pb_x +
                                 g_0_xxxxyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yzz_0[i] = 4.0 * g_0_xxxyy_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yzz_0[i] * pb_x +
                                 g_0_xxxxyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_zzz_0[i] = 4.0 * g_0_xxxyy_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_zzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_zzz_0[i] * pb_x +
                                 g_0_xxxxyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 40-50 components of targeted buffer : SKSF

    auto g_0_xxxxxyz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 40);

    auto g_0_xxxxxyz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 41);

    auto g_0_xxxxxyz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 42);

    auto g_0_xxxxxyz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 43);

    auto g_0_xxxxxyz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 44);

    auto g_0_xxxxxyz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 45);

    auto g_0_xxxxxyz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 46);

    auto g_0_xxxxxyz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 47);

    auto g_0_xxxxxyz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 48);

    auto g_0_xxxxxyz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 49);

#pragma omp simd aligned(g_0_xxxxxy_0_xxy_0,      \
                             g_0_xxxxxy_0_xxy_1,  \
                             g_0_xxxxxy_0_xyy_0,  \
                             g_0_xxxxxy_0_xyy_1,  \
                             g_0_xxxxxy_0_yyy_0,  \
                             g_0_xxxxxy_0_yyy_1,  \
                             g_0_xxxxxyz_0_xxx_0, \
                             g_0_xxxxxyz_0_xxy_0, \
                             g_0_xxxxxyz_0_xxz_0, \
                             g_0_xxxxxyz_0_xyy_0, \
                             g_0_xxxxxyz_0_xyz_0, \
                             g_0_xxxxxyz_0_xzz_0, \
                             g_0_xxxxxyz_0_yyy_0, \
                             g_0_xxxxxyz_0_yyz_0, \
                             g_0_xxxxxyz_0_yzz_0, \
                             g_0_xxxxxyz_0_zzz_0, \
                             g_0_xxxxxz_0_xxx_0,  \
                             g_0_xxxxxz_0_xxx_1,  \
                             g_0_xxxxxz_0_xxz_0,  \
                             g_0_xxxxxz_0_xxz_1,  \
                             g_0_xxxxxz_0_xyz_0,  \
                             g_0_xxxxxz_0_xyz_1,  \
                             g_0_xxxxxz_0_xz_1,   \
                             g_0_xxxxxz_0_xzz_0,  \
                             g_0_xxxxxz_0_xzz_1,  \
                             g_0_xxxxxz_0_yyz_0,  \
                             g_0_xxxxxz_0_yyz_1,  \
                             g_0_xxxxxz_0_yz_1,   \
                             g_0_xxxxxz_0_yzz_0,  \
                             g_0_xxxxxz_0_yzz_1,  \
                             g_0_xxxxxz_0_zz_1,   \
                             g_0_xxxxxz_0_zzz_0,  \
                             g_0_xxxxxz_0_zzz_1,  \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyz_0_xxx_0[i] = g_0_xxxxxz_0_xxx_0[i] * pb_y + g_0_xxxxxz_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxy_0[i] = g_0_xxxxxy_0_xxy_0[i] * pb_z + g_0_xxxxxy_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxz_0[i] = g_0_xxxxxz_0_xxz_0[i] * pb_y + g_0_xxxxxz_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyy_0[i] = g_0_xxxxxy_0_xyy_0[i] * pb_z + g_0_xxxxxy_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xyz_0[i] = g_0_xxxxxz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyz_0[i] * pb_y + g_0_xxxxxz_0_xyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xzz_0[i] = g_0_xxxxxz_0_xzz_0[i] * pb_y + g_0_xxxxxz_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyy_0[i] = g_0_xxxxxy_0_yyy_0[i] * pb_z + g_0_xxxxxy_0_yyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_yyz_0[i] = 2.0 * g_0_xxxxxz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyz_0[i] * pb_y + g_0_xxxxxz_0_yyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yzz_0[i] = g_0_xxxxxz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yzz_0[i] * pb_y + g_0_xxxxxz_0_yzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_zzz_0[i] = g_0_xxxxxz_0_zzz_0[i] * pb_y + g_0_xxxxxz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 50-60 components of targeted buffer : SKSF

    auto g_0_xxxxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 50);

    auto g_0_xxxxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 51);

    auto g_0_xxxxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 52);

    auto g_0_xxxxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 53);

    auto g_0_xxxxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 54);

    auto g_0_xxxxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 55);

    auto g_0_xxxxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 56);

    auto g_0_xxxxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 57);

    auto g_0_xxxxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 58);

    auto g_0_xxxxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 59);

#pragma omp simd aligned(g_0_xxxxx_0_xxx_0,       \
                             g_0_xxxxx_0_xxx_1,   \
                             g_0_xxxxx_0_xxy_0,   \
                             g_0_xxxxx_0_xxy_1,   \
                             g_0_xxxxx_0_xyy_0,   \
                             g_0_xxxxx_0_xyy_1,   \
                             g_0_xxxxxz_0_xxx_0,  \
                             g_0_xxxxxz_0_xxx_1,  \
                             g_0_xxxxxz_0_xxy_0,  \
                             g_0_xxxxxz_0_xxy_1,  \
                             g_0_xxxxxz_0_xyy_0,  \
                             g_0_xxxxxz_0_xyy_1,  \
                             g_0_xxxxxzz_0_xxx_0, \
                             g_0_xxxxxzz_0_xxy_0, \
                             g_0_xxxxxzz_0_xxz_0, \
                             g_0_xxxxxzz_0_xyy_0, \
                             g_0_xxxxxzz_0_xyz_0, \
                             g_0_xxxxxzz_0_xzz_0, \
                             g_0_xxxxxzz_0_yyy_0, \
                             g_0_xxxxxzz_0_yyz_0, \
                             g_0_xxxxxzz_0_yzz_0, \
                             g_0_xxxxxzz_0_zzz_0, \
                             g_0_xxxxzz_0_xxz_0,  \
                             g_0_xxxxzz_0_xxz_1,  \
                             g_0_xxxxzz_0_xyz_0,  \
                             g_0_xxxxzz_0_xyz_1,  \
                             g_0_xxxxzz_0_xz_1,   \
                             g_0_xxxxzz_0_xzz_0,  \
                             g_0_xxxxzz_0_xzz_1,  \
                             g_0_xxxxzz_0_yyy_0,  \
                             g_0_xxxxzz_0_yyy_1,  \
                             g_0_xxxxzz_0_yyz_0,  \
                             g_0_xxxxzz_0_yyz_1,  \
                             g_0_xxxxzz_0_yz_1,   \
                             g_0_xxxxzz_0_yzz_0,  \
                             g_0_xxxxzz_0_yzz_1,  \
                             g_0_xxxxzz_0_zz_1,   \
                             g_0_xxxxzz_0_zzz_0,  \
                             g_0_xxxxzz_0_zzz_1,  \
                             g_0_xxxzz_0_xxz_0,   \
                             g_0_xxxzz_0_xxz_1,   \
                             g_0_xxxzz_0_xyz_0,   \
                             g_0_xxxzz_0_xyz_1,   \
                             g_0_xxxzz_0_xzz_0,   \
                             g_0_xxxzz_0_xzz_1,   \
                             g_0_xxxzz_0_yyy_0,   \
                             g_0_xxxzz_0_yyy_1,   \
                             g_0_xxxzz_0_yyz_0,   \
                             g_0_xxxzz_0_yyz_1,   \
                             g_0_xxxzz_0_yzz_0,   \
                             g_0_xxxzz_0_yzz_1,   \
                             g_0_xxxzz_0_zzz_0,   \
                             g_0_xxxzz_0_zzz_1,   \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzz_0_xxx_0[i] =
            g_0_xxxxx_0_xxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxx_0[i] * pb_z + g_0_xxxxxz_0_xxx_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxy_0[i] =
            g_0_xxxxx_0_xxy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxy_0[i] * pb_z + g_0_xxxxxz_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxz_0[i] = 4.0 * g_0_xxxzz_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxxzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxz_0[i] * pb_x + g_0_xxxxzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyy_0[i] =
            g_0_xxxxx_0_xyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xyy_0[i] * pb_z + g_0_xxxxxz_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xyz_0[i] = 4.0 * g_0_xxxzz_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxxxzz_0_xyz_0[i] * pb_x + g_0_xxxxzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xzz_0[i] = 4.0 * g_0_xxxzz_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_xxxxzz_0_xzz_0[i] * pb_x + g_0_xxxxzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyy_0[i] = 4.0 * g_0_xxxzz_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyy_0[i] * pb_x +
                                 g_0_xxxxzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyz_0[i] = 4.0 * g_0_xxxzz_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyz_0[i] * pb_x +
                                 g_0_xxxxzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yzz_0[i] = 4.0 * g_0_xxxzz_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yzz_0[i] * pb_x +
                                 g_0_xxxxzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_zzz_0[i] = 4.0 * g_0_xxxzz_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_zzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zzz_0[i] * pb_x +
                                 g_0_xxxxzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 60-70 components of targeted buffer : SKSF

    auto g_0_xxxxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 60);

    auto g_0_xxxxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 61);

    auto g_0_xxxxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 62);

    auto g_0_xxxxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 63);

    auto g_0_xxxxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 64);

    auto g_0_xxxxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 65);

    auto g_0_xxxxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 66);

    auto g_0_xxxxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 67);

    auto g_0_xxxxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 68);

    auto g_0_xxxxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 69);

#pragma omp simd aligned(g_0_xxxxy_0_xxx_0,       \
                             g_0_xxxxy_0_xxx_1,   \
                             g_0_xxxxy_0_xxz_0,   \
                             g_0_xxxxy_0_xxz_1,   \
                             g_0_xxxxy_0_xzz_0,   \
                             g_0_xxxxy_0_xzz_1,   \
                             g_0_xxxxyy_0_xxx_0,  \
                             g_0_xxxxyy_0_xxx_1,  \
                             g_0_xxxxyy_0_xxz_0,  \
                             g_0_xxxxyy_0_xxz_1,  \
                             g_0_xxxxyy_0_xzz_0,  \
                             g_0_xxxxyy_0_xzz_1,  \
                             g_0_xxxxyyy_0_xxx_0, \
                             g_0_xxxxyyy_0_xxy_0, \
                             g_0_xxxxyyy_0_xxz_0, \
                             g_0_xxxxyyy_0_xyy_0, \
                             g_0_xxxxyyy_0_xyz_0, \
                             g_0_xxxxyyy_0_xzz_0, \
                             g_0_xxxxyyy_0_yyy_0, \
                             g_0_xxxxyyy_0_yyz_0, \
                             g_0_xxxxyyy_0_yzz_0, \
                             g_0_xxxxyyy_0_zzz_0, \
                             g_0_xxxyyy_0_xxy_0,  \
                             g_0_xxxyyy_0_xxy_1,  \
                             g_0_xxxyyy_0_xy_1,   \
                             g_0_xxxyyy_0_xyy_0,  \
                             g_0_xxxyyy_0_xyy_1,  \
                             g_0_xxxyyy_0_xyz_0,  \
                             g_0_xxxyyy_0_xyz_1,  \
                             g_0_xxxyyy_0_yy_1,   \
                             g_0_xxxyyy_0_yyy_0,  \
                             g_0_xxxyyy_0_yyy_1,  \
                             g_0_xxxyyy_0_yyz_0,  \
                             g_0_xxxyyy_0_yyz_1,  \
                             g_0_xxxyyy_0_yz_1,   \
                             g_0_xxxyyy_0_yzz_0,  \
                             g_0_xxxyyy_0_yzz_1,  \
                             g_0_xxxyyy_0_zzz_0,  \
                             g_0_xxxyyy_0_zzz_1,  \
                             g_0_xxyyy_0_xxy_0,   \
                             g_0_xxyyy_0_xxy_1,   \
                             g_0_xxyyy_0_xyy_0,   \
                             g_0_xxyyy_0_xyy_1,   \
                             g_0_xxyyy_0_xyz_0,   \
                             g_0_xxyyy_0_xyz_1,   \
                             g_0_xxyyy_0_yyy_0,   \
                             g_0_xxyyy_0_yyy_1,   \
                             g_0_xxyyy_0_yyz_0,   \
                             g_0_xxyyy_0_yyz_1,   \
                             g_0_xxyyy_0_yzz_0,   \
                             g_0_xxyyy_0_yzz_1,   \
                             g_0_xxyyy_0_zzz_0,   \
                             g_0_xxyyy_0_zzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyy_0_xxx_0[i] = 2.0 * g_0_xxxxy_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxx_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxx_0[i] * pb_y +
                                 g_0_xxxxyy_0_xxx_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxy_0[i] = 3.0 * g_0_xxyyy_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxy_0[i] * pb_x + g_0_xxxyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxz_0[i] = 2.0 * g_0_xxxxy_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxz_0[i] * pb_y +
                                 g_0_xxxxyy_0_xxz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xyy_0[i] = 3.0 * g_0_xxyyy_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yy_1[i] * fi_abcd_0 +
                                 g_0_xxxyyy_0_xyy_0[i] * pb_x + g_0_xxxyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyz_0[i] = 3.0 * g_0_xxyyy_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxxyyy_0_xyz_0[i] * pb_x + g_0_xxxyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xzz_0[i] = 2.0 * g_0_xxxxy_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xzz_0[i] * pb_y +
                                 g_0_xxxxyy_0_xzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_yyy_0[i] = 3.0 * g_0_xxyyy_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyy_0[i] * pb_x +
                                 g_0_xxxyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyz_0[i] = 3.0 * g_0_xxyyy_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyz_0[i] * pb_x +
                                 g_0_xxxyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yzz_0[i] = 3.0 * g_0_xxyyy_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yzz_0[i] * pb_x +
                                 g_0_xxxyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_zzz_0[i] = 3.0 * g_0_xxyyy_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_zzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_zzz_0[i] * pb_x +
                                 g_0_xxxyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 70-80 components of targeted buffer : SKSF

    auto g_0_xxxxyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 70);

    auto g_0_xxxxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 71);

    auto g_0_xxxxyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 72);

    auto g_0_xxxxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 73);

    auto g_0_xxxxyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 74);

    auto g_0_xxxxyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 75);

    auto g_0_xxxxyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 76);

    auto g_0_xxxxyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 77);

    auto g_0_xxxxyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 78);

    auto g_0_xxxxyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 79);

#pragma omp simd aligned(g_0_xxxxyy_0_xx_1,       \
                             g_0_xxxxyy_0_xxx_0,  \
                             g_0_xxxxyy_0_xxx_1,  \
                             g_0_xxxxyy_0_xxy_0,  \
                             g_0_xxxxyy_0_xxy_1,  \
                             g_0_xxxxyy_0_xxz_0,  \
                             g_0_xxxxyy_0_xxz_1,  \
                             g_0_xxxxyy_0_xy_1,   \
                             g_0_xxxxyy_0_xyy_0,  \
                             g_0_xxxxyy_0_xyy_1,  \
                             g_0_xxxxyy_0_xyz_0,  \
                             g_0_xxxxyy_0_xyz_1,  \
                             g_0_xxxxyy_0_xz_1,   \
                             g_0_xxxxyy_0_xzz_0,  \
                             g_0_xxxxyy_0_xzz_1,  \
                             g_0_xxxxyy_0_yy_1,   \
                             g_0_xxxxyy_0_yyy_0,  \
                             g_0_xxxxyy_0_yyy_1,  \
                             g_0_xxxxyy_0_yyz_0,  \
                             g_0_xxxxyy_0_yyz_1,  \
                             g_0_xxxxyy_0_yz_1,   \
                             g_0_xxxxyy_0_yzz_0,  \
                             g_0_xxxxyy_0_yzz_1,  \
                             g_0_xxxxyy_0_zz_1,   \
                             g_0_xxxxyy_0_zzz_0,  \
                             g_0_xxxxyy_0_zzz_1,  \
                             g_0_xxxxyyz_0_xxx_0, \
                             g_0_xxxxyyz_0_xxy_0, \
                             g_0_xxxxyyz_0_xxz_0, \
                             g_0_xxxxyyz_0_xyy_0, \
                             g_0_xxxxyyz_0_xyz_0, \
                             g_0_xxxxyyz_0_xzz_0, \
                             g_0_xxxxyyz_0_yyy_0, \
                             g_0_xxxxyyz_0_yyz_0, \
                             g_0_xxxxyyz_0_yzz_0, \
                             g_0_xxxxyyz_0_zzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyz_0_xxx_0[i] = g_0_xxxxyy_0_xxx_0[i] * pb_z + g_0_xxxxyy_0_xxx_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxy_0[i] = g_0_xxxxyy_0_xxy_0[i] * pb_z + g_0_xxxxyy_0_xxy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxz_0[i] = g_0_xxxxyy_0_xx_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxz_0[i] * pb_z + g_0_xxxxyy_0_xxz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyy_0[i] = g_0_xxxxyy_0_xyy_0[i] * pb_z + g_0_xxxxyy_0_xyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyz_0[i] = g_0_xxxxyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyz_0[i] * pb_z + g_0_xxxxyy_0_xyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xzz_0[i] = 2.0 * g_0_xxxxyy_0_xz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xzz_0[i] * pb_z + g_0_xxxxyy_0_xzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyy_0[i] = g_0_xxxxyy_0_yyy_0[i] * pb_z + g_0_xxxxyy_0_yyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyz_0[i] = g_0_xxxxyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyz_0[i] * pb_z + g_0_xxxxyy_0_yyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yzz_0[i] = 2.0 * g_0_xxxxyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yzz_0[i] * pb_z + g_0_xxxxyy_0_yzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_zzz_0[i] = 3.0 * g_0_xxxxyy_0_zz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_zzz_0[i] * pb_z + g_0_xxxxyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 80-90 components of targeted buffer : SKSF

    auto g_0_xxxxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 80);

    auto g_0_xxxxyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 81);

    auto g_0_xxxxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 82);

    auto g_0_xxxxyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 83);

    auto g_0_xxxxyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 84);

    auto g_0_xxxxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 85);

    auto g_0_xxxxyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 86);

    auto g_0_xxxxyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 87);

    auto g_0_xxxxyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 88);

    auto g_0_xxxxyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 89);

#pragma omp simd aligned(g_0_xxxxyzz_0_xxx_0,     \
                             g_0_xxxxyzz_0_xxy_0, \
                             g_0_xxxxyzz_0_xxz_0, \
                             g_0_xxxxyzz_0_xyy_0, \
                             g_0_xxxxyzz_0_xyz_0, \
                             g_0_xxxxyzz_0_xzz_0, \
                             g_0_xxxxyzz_0_yyy_0, \
                             g_0_xxxxyzz_0_yyz_0, \
                             g_0_xxxxyzz_0_yzz_0, \
                             g_0_xxxxyzz_0_zzz_0, \
                             g_0_xxxxzz_0_xx_1,   \
                             g_0_xxxxzz_0_xxx_0,  \
                             g_0_xxxxzz_0_xxx_1,  \
                             g_0_xxxxzz_0_xxy_0,  \
                             g_0_xxxxzz_0_xxy_1,  \
                             g_0_xxxxzz_0_xxz_0,  \
                             g_0_xxxxzz_0_xxz_1,  \
                             g_0_xxxxzz_0_xy_1,   \
                             g_0_xxxxzz_0_xyy_0,  \
                             g_0_xxxxzz_0_xyy_1,  \
                             g_0_xxxxzz_0_xyz_0,  \
                             g_0_xxxxzz_0_xyz_1,  \
                             g_0_xxxxzz_0_xz_1,   \
                             g_0_xxxxzz_0_xzz_0,  \
                             g_0_xxxxzz_0_xzz_1,  \
                             g_0_xxxxzz_0_yy_1,   \
                             g_0_xxxxzz_0_yyy_0,  \
                             g_0_xxxxzz_0_yyy_1,  \
                             g_0_xxxxzz_0_yyz_0,  \
                             g_0_xxxxzz_0_yyz_1,  \
                             g_0_xxxxzz_0_yz_1,   \
                             g_0_xxxxzz_0_yzz_0,  \
                             g_0_xxxxzz_0_yzz_1,  \
                             g_0_xxxxzz_0_zz_1,   \
                             g_0_xxxxzz_0_zzz_0,  \
                             g_0_xxxxzz_0_zzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzz_0_xxx_0[i] = g_0_xxxxzz_0_xxx_0[i] * pb_y + g_0_xxxxzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxy_0[i] = g_0_xxxxzz_0_xx_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxy_0[i] * pb_y + g_0_xxxxzz_0_xxy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxz_0[i] = g_0_xxxxzz_0_xxz_0[i] * pb_y + g_0_xxxxzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyy_0[i] = 2.0 * g_0_xxxxzz_0_xy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyy_0[i] * pb_y + g_0_xxxxzz_0_xyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyz_0[i] = g_0_xxxxzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyz_0[i] * pb_y + g_0_xxxxzz_0_xyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xzz_0[i] = g_0_xxxxzz_0_xzz_0[i] * pb_y + g_0_xxxxzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyy_0[i] = 3.0 * g_0_xxxxzz_0_yy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyy_0[i] * pb_y + g_0_xxxxzz_0_yyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyz_0[i] = 2.0 * g_0_xxxxzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyz_0[i] * pb_y + g_0_xxxxzz_0_yyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yzz_0[i] = g_0_xxxxzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yzz_0[i] * pb_y + g_0_xxxxzz_0_yzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_zzz_0[i] = g_0_xxxxzz_0_zzz_0[i] * pb_y + g_0_xxxxzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 90-100 components of targeted buffer : SKSF

    auto g_0_xxxxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 90);

    auto g_0_xxxxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 91);

    auto g_0_xxxxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 92);

    auto g_0_xxxxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 93);

    auto g_0_xxxxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 94);

    auto g_0_xxxxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 95);

    auto g_0_xxxxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 96);

    auto g_0_xxxxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 97);

    auto g_0_xxxxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 98);

    auto g_0_xxxxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 99);

#pragma omp simd aligned(g_0_xxxxz_0_xxx_0,       \
                             g_0_xxxxz_0_xxx_1,   \
                             g_0_xxxxz_0_xxy_0,   \
                             g_0_xxxxz_0_xxy_1,   \
                             g_0_xxxxz_0_xyy_0,   \
                             g_0_xxxxz_0_xyy_1,   \
                             g_0_xxxxzz_0_xxx_0,  \
                             g_0_xxxxzz_0_xxx_1,  \
                             g_0_xxxxzz_0_xxy_0,  \
                             g_0_xxxxzz_0_xxy_1,  \
                             g_0_xxxxzz_0_xyy_0,  \
                             g_0_xxxxzz_0_xyy_1,  \
                             g_0_xxxxzzz_0_xxx_0, \
                             g_0_xxxxzzz_0_xxy_0, \
                             g_0_xxxxzzz_0_xxz_0, \
                             g_0_xxxxzzz_0_xyy_0, \
                             g_0_xxxxzzz_0_xyz_0, \
                             g_0_xxxxzzz_0_xzz_0, \
                             g_0_xxxxzzz_0_yyy_0, \
                             g_0_xxxxzzz_0_yyz_0, \
                             g_0_xxxxzzz_0_yzz_0, \
                             g_0_xxxxzzz_0_zzz_0, \
                             g_0_xxxzzz_0_xxz_0,  \
                             g_0_xxxzzz_0_xxz_1,  \
                             g_0_xxxzzz_0_xyz_0,  \
                             g_0_xxxzzz_0_xyz_1,  \
                             g_0_xxxzzz_0_xz_1,   \
                             g_0_xxxzzz_0_xzz_0,  \
                             g_0_xxxzzz_0_xzz_1,  \
                             g_0_xxxzzz_0_yyy_0,  \
                             g_0_xxxzzz_0_yyy_1,  \
                             g_0_xxxzzz_0_yyz_0,  \
                             g_0_xxxzzz_0_yyz_1,  \
                             g_0_xxxzzz_0_yz_1,   \
                             g_0_xxxzzz_0_yzz_0,  \
                             g_0_xxxzzz_0_yzz_1,  \
                             g_0_xxxzzz_0_zz_1,   \
                             g_0_xxxzzz_0_zzz_0,  \
                             g_0_xxxzzz_0_zzz_1,  \
                             g_0_xxzzz_0_xxz_0,   \
                             g_0_xxzzz_0_xxz_1,   \
                             g_0_xxzzz_0_xyz_0,   \
                             g_0_xxzzz_0_xyz_1,   \
                             g_0_xxzzz_0_xzz_0,   \
                             g_0_xxzzz_0_xzz_1,   \
                             g_0_xxzzz_0_yyy_0,   \
                             g_0_xxzzz_0_yyy_1,   \
                             g_0_xxzzz_0_yyz_0,   \
                             g_0_xxzzz_0_yyz_1,   \
                             g_0_xxzzz_0_yzz_0,   \
                             g_0_xxzzz_0_yzz_1,   \
                             g_0_xxzzz_0_zzz_0,   \
                             g_0_xxzzz_0_zzz_1,   \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzz_0_xxx_0[i] = 2.0 * g_0_xxxxz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxx_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxx_0[i] * pb_z +
                                 g_0_xxxxzz_0_xxx_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxy_0[i] = 2.0 * g_0_xxxxz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxy_0[i] * pb_z +
                                 g_0_xxxxzz_0_xxy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxz_0[i] = 3.0 * g_0_xxzzz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxz_0[i] * pb_x + g_0_xxxzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyy_0[i] = 2.0 * g_0_xxxxz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xyy_0[i] * pb_z +
                                 g_0_xxxxzz_0_xyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xyz_0[i] = 3.0 * g_0_xxzzz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxxzzz_0_xyz_0[i] * pb_x + g_0_xxxzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xzz_0[i] = 3.0 * g_0_xxzzz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_xxxzzz_0_xzz_0[i] * pb_x + g_0_xxxzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyy_0[i] = 3.0 * g_0_xxzzz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyy_0[i] * pb_x +
                                 g_0_xxxzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyz_0[i] = 3.0 * g_0_xxzzz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyz_0[i] * pb_x +
                                 g_0_xxxzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yzz_0[i] = 3.0 * g_0_xxzzz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yzz_0[i] * pb_x +
                                 g_0_xxxzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_zzz_0[i] = 3.0 * g_0_xxzzz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_zzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zzz_0[i] * pb_x +
                                 g_0_xxxzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 100-110 components of targeted buffer : SKSF

    auto g_0_xxxyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 100);

    auto g_0_xxxyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 101);

    auto g_0_xxxyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 102);

    auto g_0_xxxyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 103);

    auto g_0_xxxyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 104);

    auto g_0_xxxyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 105);

    auto g_0_xxxyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 106);

    auto g_0_xxxyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 107);

    auto g_0_xxxyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 108);

    auto g_0_xxxyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 109);

#pragma omp simd aligned(g_0_xxxyy_0_xxx_0,       \
                             g_0_xxxyy_0_xxx_1,   \
                             g_0_xxxyy_0_xxz_0,   \
                             g_0_xxxyy_0_xxz_1,   \
                             g_0_xxxyy_0_xzz_0,   \
                             g_0_xxxyy_0_xzz_1,   \
                             g_0_xxxyyy_0_xxx_0,  \
                             g_0_xxxyyy_0_xxx_1,  \
                             g_0_xxxyyy_0_xxz_0,  \
                             g_0_xxxyyy_0_xxz_1,  \
                             g_0_xxxyyy_0_xzz_0,  \
                             g_0_xxxyyy_0_xzz_1,  \
                             g_0_xxxyyyy_0_xxx_0, \
                             g_0_xxxyyyy_0_xxy_0, \
                             g_0_xxxyyyy_0_xxz_0, \
                             g_0_xxxyyyy_0_xyy_0, \
                             g_0_xxxyyyy_0_xyz_0, \
                             g_0_xxxyyyy_0_xzz_0, \
                             g_0_xxxyyyy_0_yyy_0, \
                             g_0_xxxyyyy_0_yyz_0, \
                             g_0_xxxyyyy_0_yzz_0, \
                             g_0_xxxyyyy_0_zzz_0, \
                             g_0_xxyyyy_0_xxy_0,  \
                             g_0_xxyyyy_0_xxy_1,  \
                             g_0_xxyyyy_0_xy_1,   \
                             g_0_xxyyyy_0_xyy_0,  \
                             g_0_xxyyyy_0_xyy_1,  \
                             g_0_xxyyyy_0_xyz_0,  \
                             g_0_xxyyyy_0_xyz_1,  \
                             g_0_xxyyyy_0_yy_1,   \
                             g_0_xxyyyy_0_yyy_0,  \
                             g_0_xxyyyy_0_yyy_1,  \
                             g_0_xxyyyy_0_yyz_0,  \
                             g_0_xxyyyy_0_yyz_1,  \
                             g_0_xxyyyy_0_yz_1,   \
                             g_0_xxyyyy_0_yzz_0,  \
                             g_0_xxyyyy_0_yzz_1,  \
                             g_0_xxyyyy_0_zzz_0,  \
                             g_0_xxyyyy_0_zzz_1,  \
                             g_0_xyyyy_0_xxy_0,   \
                             g_0_xyyyy_0_xxy_1,   \
                             g_0_xyyyy_0_xyy_0,   \
                             g_0_xyyyy_0_xyy_1,   \
                             g_0_xyyyy_0_xyz_0,   \
                             g_0_xyyyy_0_xyz_1,   \
                             g_0_xyyyy_0_yyy_0,   \
                             g_0_xyyyy_0_yyy_1,   \
                             g_0_xyyyy_0_yyz_0,   \
                             g_0_xyyyy_0_yyz_1,   \
                             g_0_xyyyy_0_yzz_0,   \
                             g_0_xyyyy_0_yzz_1,   \
                             g_0_xyyyy_0_zzz_0,   \
                             g_0_xyyyy_0_zzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyy_0_xxx_0[i] = 3.0 * g_0_xxxyy_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxx_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxx_0[i] * pb_y +
                                 g_0_xxxyyy_0_xxx_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxy_0[i] = 2.0 * g_0_xyyyy_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxy_0[i] * pb_x + g_0_xxyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxz_0[i] = 3.0 * g_0_xxxyy_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxz_0[i] * pb_y +
                                 g_0_xxxyyy_0_xxz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xyy_0[i] = 2.0 * g_0_xyyyy_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yy_1[i] * fi_abcd_0 +
                                 g_0_xxyyyy_0_xyy_0[i] * pb_x + g_0_xxyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyz_0[i] = 2.0 * g_0_xyyyy_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxyyyy_0_xyz_0[i] * pb_x + g_0_xxyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xzz_0[i] = 3.0 * g_0_xxxyy_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xzz_0[i] * pb_y +
                                 g_0_xxxyyy_0_xzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_yyy_0[i] = 2.0 * g_0_xyyyy_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyy_0[i] * pb_x +
                                 g_0_xxyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyz_0[i] = 2.0 * g_0_xyyyy_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyz_0[i] * pb_x +
                                 g_0_xxyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yzz_0[i] = 2.0 * g_0_xyyyy_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yzz_0[i] * pb_x +
                                 g_0_xxyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_zzz_0[i] = 2.0 * g_0_xyyyy_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_zzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_zzz_0[i] * pb_x +
                                 g_0_xxyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 110-120 components of targeted buffer : SKSF

    auto g_0_xxxyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 110);

    auto g_0_xxxyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 111);

    auto g_0_xxxyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 112);

    auto g_0_xxxyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 113);

    auto g_0_xxxyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 114);

    auto g_0_xxxyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 115);

    auto g_0_xxxyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 116);

    auto g_0_xxxyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 117);

    auto g_0_xxxyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 118);

    auto g_0_xxxyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 119);

#pragma omp simd aligned(g_0_xxxyyy_0_xx_1,       \
                             g_0_xxxyyy_0_xxx_0,  \
                             g_0_xxxyyy_0_xxx_1,  \
                             g_0_xxxyyy_0_xxy_0,  \
                             g_0_xxxyyy_0_xxy_1,  \
                             g_0_xxxyyy_0_xxz_0,  \
                             g_0_xxxyyy_0_xxz_1,  \
                             g_0_xxxyyy_0_xy_1,   \
                             g_0_xxxyyy_0_xyy_0,  \
                             g_0_xxxyyy_0_xyy_1,  \
                             g_0_xxxyyy_0_xyz_0,  \
                             g_0_xxxyyy_0_xyz_1,  \
                             g_0_xxxyyy_0_xz_1,   \
                             g_0_xxxyyy_0_xzz_0,  \
                             g_0_xxxyyy_0_xzz_1,  \
                             g_0_xxxyyy_0_yy_1,   \
                             g_0_xxxyyy_0_yyy_0,  \
                             g_0_xxxyyy_0_yyy_1,  \
                             g_0_xxxyyy_0_yyz_0,  \
                             g_0_xxxyyy_0_yyz_1,  \
                             g_0_xxxyyy_0_yz_1,   \
                             g_0_xxxyyy_0_yzz_0,  \
                             g_0_xxxyyy_0_yzz_1,  \
                             g_0_xxxyyy_0_zz_1,   \
                             g_0_xxxyyy_0_zzz_0,  \
                             g_0_xxxyyy_0_zzz_1,  \
                             g_0_xxxyyyz_0_xxx_0, \
                             g_0_xxxyyyz_0_xxy_0, \
                             g_0_xxxyyyz_0_xxz_0, \
                             g_0_xxxyyyz_0_xyy_0, \
                             g_0_xxxyyyz_0_xyz_0, \
                             g_0_xxxyyyz_0_xzz_0, \
                             g_0_xxxyyyz_0_yyy_0, \
                             g_0_xxxyyyz_0_yyz_0, \
                             g_0_xxxyyyz_0_yzz_0, \
                             g_0_xxxyyyz_0_zzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyz_0_xxx_0[i] = g_0_xxxyyy_0_xxx_0[i] * pb_z + g_0_xxxyyy_0_xxx_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxy_0[i] = g_0_xxxyyy_0_xxy_0[i] * pb_z + g_0_xxxyyy_0_xxy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxz_0[i] = g_0_xxxyyy_0_xx_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxz_0[i] * pb_z + g_0_xxxyyy_0_xxz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyy_0[i] = g_0_xxxyyy_0_xyy_0[i] * pb_z + g_0_xxxyyy_0_xyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyz_0[i] = g_0_xxxyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyz_0[i] * pb_z + g_0_xxxyyy_0_xyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xzz_0[i] = 2.0 * g_0_xxxyyy_0_xz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xzz_0[i] * pb_z + g_0_xxxyyy_0_xzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyy_0[i] = g_0_xxxyyy_0_yyy_0[i] * pb_z + g_0_xxxyyy_0_yyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyz_0[i] = g_0_xxxyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyz_0[i] * pb_z + g_0_xxxyyy_0_yyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yzz_0[i] = 2.0 * g_0_xxxyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yzz_0[i] * pb_z + g_0_xxxyyy_0_yzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_zzz_0[i] = 3.0 * g_0_xxxyyy_0_zz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_zzz_0[i] * pb_z + g_0_xxxyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 120-130 components of targeted buffer : SKSF

    auto g_0_xxxyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 120);

    auto g_0_xxxyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 121);

    auto g_0_xxxyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 122);

    auto g_0_xxxyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 123);

    auto g_0_xxxyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 124);

    auto g_0_xxxyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 125);

    auto g_0_xxxyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 126);

    auto g_0_xxxyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 127);

    auto g_0_xxxyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 128);

    auto g_0_xxxyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 129);

#pragma omp simd aligned(g_0_xxxyy_0_xxy_0,       \
                             g_0_xxxyy_0_xxy_1,   \
                             g_0_xxxyy_0_xyy_0,   \
                             g_0_xxxyy_0_xyy_1,   \
                             g_0_xxxyyz_0_xxy_0,  \
                             g_0_xxxyyz_0_xxy_1,  \
                             g_0_xxxyyz_0_xyy_0,  \
                             g_0_xxxyyz_0_xyy_1,  \
                             g_0_xxxyyzz_0_xxx_0, \
                             g_0_xxxyyzz_0_xxy_0, \
                             g_0_xxxyyzz_0_xxz_0, \
                             g_0_xxxyyzz_0_xyy_0, \
                             g_0_xxxyyzz_0_xyz_0, \
                             g_0_xxxyyzz_0_xzz_0, \
                             g_0_xxxyyzz_0_yyy_0, \
                             g_0_xxxyyzz_0_yyz_0, \
                             g_0_xxxyyzz_0_yzz_0, \
                             g_0_xxxyyzz_0_zzz_0, \
                             g_0_xxxyzz_0_xxx_0,  \
                             g_0_xxxyzz_0_xxx_1,  \
                             g_0_xxxyzz_0_xxz_0,  \
                             g_0_xxxyzz_0_xxz_1,  \
                             g_0_xxxyzz_0_xzz_0,  \
                             g_0_xxxyzz_0_xzz_1,  \
                             g_0_xxxzz_0_xxx_0,   \
                             g_0_xxxzz_0_xxx_1,   \
                             g_0_xxxzz_0_xxz_0,   \
                             g_0_xxxzz_0_xxz_1,   \
                             g_0_xxxzz_0_xzz_0,   \
                             g_0_xxxzz_0_xzz_1,   \
                             g_0_xxyyzz_0_xyz_0,  \
                             g_0_xxyyzz_0_xyz_1,  \
                             g_0_xxyyzz_0_yyy_0,  \
                             g_0_xxyyzz_0_yyy_1,  \
                             g_0_xxyyzz_0_yyz_0,  \
                             g_0_xxyyzz_0_yyz_1,  \
                             g_0_xxyyzz_0_yz_1,   \
                             g_0_xxyyzz_0_yzz_0,  \
                             g_0_xxyyzz_0_yzz_1,  \
                             g_0_xxyyzz_0_zzz_0,  \
                             g_0_xxyyzz_0_zzz_1,  \
                             g_0_xyyzz_0_xyz_0,   \
                             g_0_xyyzz_0_xyz_1,   \
                             g_0_xyyzz_0_yyy_0,   \
                             g_0_xyyzz_0_yyy_1,   \
                             g_0_xyyzz_0_yyz_0,   \
                             g_0_xyyzz_0_yyz_1,   \
                             g_0_xyyzz_0_yzz_0,   \
                             g_0_xyyzz_0_yzz_1,   \
                             g_0_xyyzz_0_zzz_0,   \
                             g_0_xyyzz_0_zzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzz_0_xxx_0[i] =
            g_0_xxxzz_0_xxx_0[i] * fi_ab_0 - g_0_xxxzz_0_xxx_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxx_0[i] * pb_y + g_0_xxxyzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxy_0[i] =
            g_0_xxxyy_0_xxy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxy_0[i] * pb_z + g_0_xxxyyz_0_xxy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxz_0[i] =
            g_0_xxxzz_0_xxz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxz_0[i] * pb_y + g_0_xxxyzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xyy_0[i] =
            g_0_xxxyy_0_xyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xyy_0[i] * pb_z + g_0_xxxyyz_0_xyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xyz_0[i] = 2.0 * g_0_xyyzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxyyzz_0_xyz_0[i] * pb_x + g_0_xxyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xzz_0[i] =
            g_0_xxxzz_0_xzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xzz_0[i] * pb_y + g_0_xxxyzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_yyy_0[i] = 2.0 * g_0_xyyzz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyy_0[i] * pb_x +
                                 g_0_xxyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyz_0[i] = 2.0 * g_0_xyyzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyz_0[i] * pb_x +
                                 g_0_xxyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yzz_0[i] = 2.0 * g_0_xyyzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yzz_0[i] * pb_x +
                                 g_0_xxyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_zzz_0[i] = 2.0 * g_0_xyyzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_zzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_zzz_0[i] * pb_x +
                                 g_0_xxyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 130-140 components of targeted buffer : SKSF

    auto g_0_xxxyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 130);

    auto g_0_xxxyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 131);

    auto g_0_xxxyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 132);

    auto g_0_xxxyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 133);

    auto g_0_xxxyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 134);

    auto g_0_xxxyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 135);

    auto g_0_xxxyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 136);

    auto g_0_xxxyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 137);

    auto g_0_xxxyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 138);

    auto g_0_xxxyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 139);

#pragma omp simd aligned(g_0_xxxyzzz_0_xxx_0,     \
                             g_0_xxxyzzz_0_xxy_0, \
                             g_0_xxxyzzz_0_xxz_0, \
                             g_0_xxxyzzz_0_xyy_0, \
                             g_0_xxxyzzz_0_xyz_0, \
                             g_0_xxxyzzz_0_xzz_0, \
                             g_0_xxxyzzz_0_yyy_0, \
                             g_0_xxxyzzz_0_yyz_0, \
                             g_0_xxxyzzz_0_yzz_0, \
                             g_0_xxxyzzz_0_zzz_0, \
                             g_0_xxxzzz_0_xx_1,   \
                             g_0_xxxzzz_0_xxx_0,  \
                             g_0_xxxzzz_0_xxx_1,  \
                             g_0_xxxzzz_0_xxy_0,  \
                             g_0_xxxzzz_0_xxy_1,  \
                             g_0_xxxzzz_0_xxz_0,  \
                             g_0_xxxzzz_0_xxz_1,  \
                             g_0_xxxzzz_0_xy_1,   \
                             g_0_xxxzzz_0_xyy_0,  \
                             g_0_xxxzzz_0_xyy_1,  \
                             g_0_xxxzzz_0_xyz_0,  \
                             g_0_xxxzzz_0_xyz_1,  \
                             g_0_xxxzzz_0_xz_1,   \
                             g_0_xxxzzz_0_xzz_0,  \
                             g_0_xxxzzz_0_xzz_1,  \
                             g_0_xxxzzz_0_yy_1,   \
                             g_0_xxxzzz_0_yyy_0,  \
                             g_0_xxxzzz_0_yyy_1,  \
                             g_0_xxxzzz_0_yyz_0,  \
                             g_0_xxxzzz_0_yyz_1,  \
                             g_0_xxxzzz_0_yz_1,   \
                             g_0_xxxzzz_0_yzz_0,  \
                             g_0_xxxzzz_0_yzz_1,  \
                             g_0_xxxzzz_0_zz_1,   \
                             g_0_xxxzzz_0_zzz_0,  \
                             g_0_xxxzzz_0_zzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzz_0_xxx_0[i] = g_0_xxxzzz_0_xxx_0[i] * pb_y + g_0_xxxzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxy_0[i] = g_0_xxxzzz_0_xx_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxy_0[i] * pb_y + g_0_xxxzzz_0_xxy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxz_0[i] = g_0_xxxzzz_0_xxz_0[i] * pb_y + g_0_xxxzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyy_0[i] = 2.0 * g_0_xxxzzz_0_xy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyy_0[i] * pb_y + g_0_xxxzzz_0_xyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyz_0[i] = g_0_xxxzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyz_0[i] * pb_y + g_0_xxxzzz_0_xyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xzz_0[i] = g_0_xxxzzz_0_xzz_0[i] * pb_y + g_0_xxxzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyy_0[i] = 3.0 * g_0_xxxzzz_0_yy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyy_0[i] * pb_y + g_0_xxxzzz_0_yyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyz_0[i] = 2.0 * g_0_xxxzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyz_0[i] * pb_y + g_0_xxxzzz_0_yyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yzz_0[i] = g_0_xxxzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yzz_0[i] * pb_y + g_0_xxxzzz_0_yzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_zzz_0[i] = g_0_xxxzzz_0_zzz_0[i] * pb_y + g_0_xxxzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 140-150 components of targeted buffer : SKSF

    auto g_0_xxxzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 140);

    auto g_0_xxxzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 141);

    auto g_0_xxxzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 142);

    auto g_0_xxxzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 143);

    auto g_0_xxxzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 144);

    auto g_0_xxxzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 145);

    auto g_0_xxxzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 146);

    auto g_0_xxxzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 147);

    auto g_0_xxxzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 148);

    auto g_0_xxxzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 149);

#pragma omp simd aligned(g_0_xxxzz_0_xxx_0,       \
                             g_0_xxxzz_0_xxx_1,   \
                             g_0_xxxzz_0_xxy_0,   \
                             g_0_xxxzz_0_xxy_1,   \
                             g_0_xxxzz_0_xyy_0,   \
                             g_0_xxxzz_0_xyy_1,   \
                             g_0_xxxzzz_0_xxx_0,  \
                             g_0_xxxzzz_0_xxx_1,  \
                             g_0_xxxzzz_0_xxy_0,  \
                             g_0_xxxzzz_0_xxy_1,  \
                             g_0_xxxzzz_0_xyy_0,  \
                             g_0_xxxzzz_0_xyy_1,  \
                             g_0_xxxzzzz_0_xxx_0, \
                             g_0_xxxzzzz_0_xxy_0, \
                             g_0_xxxzzzz_0_xxz_0, \
                             g_0_xxxzzzz_0_xyy_0, \
                             g_0_xxxzzzz_0_xyz_0, \
                             g_0_xxxzzzz_0_xzz_0, \
                             g_0_xxxzzzz_0_yyy_0, \
                             g_0_xxxzzzz_0_yyz_0, \
                             g_0_xxxzzzz_0_yzz_0, \
                             g_0_xxxzzzz_0_zzz_0, \
                             g_0_xxzzzz_0_xxz_0,  \
                             g_0_xxzzzz_0_xxz_1,  \
                             g_0_xxzzzz_0_xyz_0,  \
                             g_0_xxzzzz_0_xyz_1,  \
                             g_0_xxzzzz_0_xz_1,   \
                             g_0_xxzzzz_0_xzz_0,  \
                             g_0_xxzzzz_0_xzz_1,  \
                             g_0_xxzzzz_0_yyy_0,  \
                             g_0_xxzzzz_0_yyy_1,  \
                             g_0_xxzzzz_0_yyz_0,  \
                             g_0_xxzzzz_0_yyz_1,  \
                             g_0_xxzzzz_0_yz_1,   \
                             g_0_xxzzzz_0_yzz_0,  \
                             g_0_xxzzzz_0_yzz_1,  \
                             g_0_xxzzzz_0_zz_1,   \
                             g_0_xxzzzz_0_zzz_0,  \
                             g_0_xxzzzz_0_zzz_1,  \
                             g_0_xzzzz_0_xxz_0,   \
                             g_0_xzzzz_0_xxz_1,   \
                             g_0_xzzzz_0_xyz_0,   \
                             g_0_xzzzz_0_xyz_1,   \
                             g_0_xzzzz_0_xzz_0,   \
                             g_0_xzzzz_0_xzz_1,   \
                             g_0_xzzzz_0_yyy_0,   \
                             g_0_xzzzz_0_yyy_1,   \
                             g_0_xzzzz_0_yyz_0,   \
                             g_0_xzzzz_0_yyz_1,   \
                             g_0_xzzzz_0_yzz_0,   \
                             g_0_xzzzz_0_yzz_1,   \
                             g_0_xzzzz_0_zzz_0,   \
                             g_0_xzzzz_0_zzz_1,   \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzz_0_xxx_0[i] = 3.0 * g_0_xxxzz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxx_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxx_0[i] * pb_z +
                                 g_0_xxxzzz_0_xxx_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxy_0[i] = 3.0 * g_0_xxxzz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxy_0[i] * pb_z +
                                 g_0_xxxzzz_0_xxy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxz_0[i] = 2.0 * g_0_xzzzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxz_0[i] * pb_x + g_0_xxzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyy_0[i] = 3.0 * g_0_xxxzz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xyy_0[i] * pb_z +
                                 g_0_xxxzzz_0_xyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xyz_0[i] = 2.0 * g_0_xzzzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xxzzzz_0_xyz_0[i] * pb_x + g_0_xxzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xzz_0[i] = 2.0 * g_0_xzzzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_xxzzzz_0_xzz_0[i] * pb_x + g_0_xxzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyy_0[i] = 2.0 * g_0_xzzzz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyy_0[i] * pb_x +
                                 g_0_xxzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyz_0[i] = 2.0 * g_0_xzzzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyz_0[i] * pb_x +
                                 g_0_xxzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yzz_0[i] = 2.0 * g_0_xzzzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yzz_0[i] * pb_x +
                                 g_0_xxzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_zzz_0[i] = 2.0 * g_0_xzzzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_zzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zzz_0[i] * pb_x +
                                 g_0_xxzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 150-160 components of targeted buffer : SKSF

    auto g_0_xxyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 150);

    auto g_0_xxyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 151);

    auto g_0_xxyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 152);

    auto g_0_xxyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 153);

    auto g_0_xxyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 154);

    auto g_0_xxyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 155);

    auto g_0_xxyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 156);

    auto g_0_xxyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 157);

    auto g_0_xxyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 158);

    auto g_0_xxyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 159);

#pragma omp simd aligned(g_0_xxyyy_0_xxx_0,       \
                             g_0_xxyyy_0_xxx_1,   \
                             g_0_xxyyy_0_xxz_0,   \
                             g_0_xxyyy_0_xxz_1,   \
                             g_0_xxyyy_0_xzz_0,   \
                             g_0_xxyyy_0_xzz_1,   \
                             g_0_xxyyyy_0_xxx_0,  \
                             g_0_xxyyyy_0_xxx_1,  \
                             g_0_xxyyyy_0_xxz_0,  \
                             g_0_xxyyyy_0_xxz_1,  \
                             g_0_xxyyyy_0_xzz_0,  \
                             g_0_xxyyyy_0_xzz_1,  \
                             g_0_xxyyyyy_0_xxx_0, \
                             g_0_xxyyyyy_0_xxy_0, \
                             g_0_xxyyyyy_0_xxz_0, \
                             g_0_xxyyyyy_0_xyy_0, \
                             g_0_xxyyyyy_0_xyz_0, \
                             g_0_xxyyyyy_0_xzz_0, \
                             g_0_xxyyyyy_0_yyy_0, \
                             g_0_xxyyyyy_0_yyz_0, \
                             g_0_xxyyyyy_0_yzz_0, \
                             g_0_xxyyyyy_0_zzz_0, \
                             g_0_xyyyyy_0_xxy_0,  \
                             g_0_xyyyyy_0_xxy_1,  \
                             g_0_xyyyyy_0_xy_1,   \
                             g_0_xyyyyy_0_xyy_0,  \
                             g_0_xyyyyy_0_xyy_1,  \
                             g_0_xyyyyy_0_xyz_0,  \
                             g_0_xyyyyy_0_xyz_1,  \
                             g_0_xyyyyy_0_yy_1,   \
                             g_0_xyyyyy_0_yyy_0,  \
                             g_0_xyyyyy_0_yyy_1,  \
                             g_0_xyyyyy_0_yyz_0,  \
                             g_0_xyyyyy_0_yyz_1,  \
                             g_0_xyyyyy_0_yz_1,   \
                             g_0_xyyyyy_0_yzz_0,  \
                             g_0_xyyyyy_0_yzz_1,  \
                             g_0_xyyyyy_0_zzz_0,  \
                             g_0_xyyyyy_0_zzz_1,  \
                             g_0_yyyyy_0_xxy_0,   \
                             g_0_yyyyy_0_xxy_1,   \
                             g_0_yyyyy_0_xyy_0,   \
                             g_0_yyyyy_0_xyy_1,   \
                             g_0_yyyyy_0_xyz_0,   \
                             g_0_yyyyy_0_xyz_1,   \
                             g_0_yyyyy_0_yyy_0,   \
                             g_0_yyyyy_0_yyy_1,   \
                             g_0_yyyyy_0_yyz_0,   \
                             g_0_yyyyy_0_yyz_1,   \
                             g_0_yyyyy_0_yzz_0,   \
                             g_0_yyyyy_0_yzz_1,   \
                             g_0_yyyyy_0_zzz_0,   \
                             g_0_yyyyy_0_zzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyy_0_xxx_0[i] = 4.0 * g_0_xxyyy_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxx_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxx_0[i] * pb_y +
                                 g_0_xxyyyy_0_xxx_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxy_0[i] = g_0_yyyyy_0_xxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xy_1[i] * fi_abcd_0 +
                                 g_0_xyyyyy_0_xxy_0[i] * pb_x + g_0_xyyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxz_0[i] = 4.0 * g_0_xxyyy_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxz_0[i] * pb_y +
                                 g_0_xxyyyy_0_xxz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xyy_0[i] = g_0_yyyyy_0_xyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yy_1[i] * fi_abcd_0 +
                                 g_0_xyyyyy_0_xyy_0[i] * pb_x + g_0_xyyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyz_0[i] = g_0_yyyyy_0_xyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xyyyyy_0_xyz_0[i] * pb_x + g_0_xyyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xzz_0[i] = 4.0 * g_0_xxyyy_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xzz_0[i] * pb_y +
                                 g_0_xxyyyy_0_xzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_yyy_0[i] =
            g_0_yyyyy_0_yyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyy_0[i] * pb_x + g_0_xyyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyz_0[i] =
            g_0_yyyyy_0_yyz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyz_0[i] * pb_x + g_0_xyyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yzz_0[i] =
            g_0_yyyyy_0_yzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzz_0[i] * pb_x + g_0_xyyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_zzz_0[i] =
            g_0_yyyyy_0_zzz_0[i] * fi_ab_0 - g_0_yyyyy_0_zzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_zzz_0[i] * pb_x + g_0_xyyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 160-170 components of targeted buffer : SKSF

    auto g_0_xxyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 160);

    auto g_0_xxyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 161);

    auto g_0_xxyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 162);

    auto g_0_xxyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 163);

    auto g_0_xxyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 164);

    auto g_0_xxyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 165);

    auto g_0_xxyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 166);

    auto g_0_xxyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 167);

    auto g_0_xxyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 168);

    auto g_0_xxyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 169);

#pragma omp simd aligned(g_0_xxyyyy_0_xx_1,       \
                             g_0_xxyyyy_0_xxx_0,  \
                             g_0_xxyyyy_0_xxx_1,  \
                             g_0_xxyyyy_0_xxy_0,  \
                             g_0_xxyyyy_0_xxy_1,  \
                             g_0_xxyyyy_0_xxz_0,  \
                             g_0_xxyyyy_0_xxz_1,  \
                             g_0_xxyyyy_0_xy_1,   \
                             g_0_xxyyyy_0_xyy_0,  \
                             g_0_xxyyyy_0_xyy_1,  \
                             g_0_xxyyyy_0_xyz_0,  \
                             g_0_xxyyyy_0_xyz_1,  \
                             g_0_xxyyyy_0_xz_1,   \
                             g_0_xxyyyy_0_xzz_0,  \
                             g_0_xxyyyy_0_xzz_1,  \
                             g_0_xxyyyy_0_yy_1,   \
                             g_0_xxyyyy_0_yyy_0,  \
                             g_0_xxyyyy_0_yyy_1,  \
                             g_0_xxyyyy_0_yyz_0,  \
                             g_0_xxyyyy_0_yyz_1,  \
                             g_0_xxyyyy_0_yz_1,   \
                             g_0_xxyyyy_0_yzz_0,  \
                             g_0_xxyyyy_0_yzz_1,  \
                             g_0_xxyyyy_0_zz_1,   \
                             g_0_xxyyyy_0_zzz_0,  \
                             g_0_xxyyyy_0_zzz_1,  \
                             g_0_xxyyyyz_0_xxx_0, \
                             g_0_xxyyyyz_0_xxy_0, \
                             g_0_xxyyyyz_0_xxz_0, \
                             g_0_xxyyyyz_0_xyy_0, \
                             g_0_xxyyyyz_0_xyz_0, \
                             g_0_xxyyyyz_0_xzz_0, \
                             g_0_xxyyyyz_0_yyy_0, \
                             g_0_xxyyyyz_0_yyz_0, \
                             g_0_xxyyyyz_0_yzz_0, \
                             g_0_xxyyyyz_0_zzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyz_0_xxx_0[i] = g_0_xxyyyy_0_xxx_0[i] * pb_z + g_0_xxyyyy_0_xxx_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxy_0[i] = g_0_xxyyyy_0_xxy_0[i] * pb_z + g_0_xxyyyy_0_xxy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxz_0[i] = g_0_xxyyyy_0_xx_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxz_0[i] * pb_z + g_0_xxyyyy_0_xxz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyy_0[i] = g_0_xxyyyy_0_xyy_0[i] * pb_z + g_0_xxyyyy_0_xyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyz_0[i] = g_0_xxyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyz_0[i] * pb_z + g_0_xxyyyy_0_xyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xzz_0[i] = 2.0 * g_0_xxyyyy_0_xz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xzz_0[i] * pb_z + g_0_xxyyyy_0_xzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyy_0[i] = g_0_xxyyyy_0_yyy_0[i] * pb_z + g_0_xxyyyy_0_yyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyz_0[i] = g_0_xxyyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyz_0[i] * pb_z + g_0_xxyyyy_0_yyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yzz_0[i] = 2.0 * g_0_xxyyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yzz_0[i] * pb_z + g_0_xxyyyy_0_yzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_zzz_0[i] = 3.0 * g_0_xxyyyy_0_zz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_zzz_0[i] * pb_z + g_0_xxyyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 170-180 components of targeted buffer : SKSF

    auto g_0_xxyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 170);

    auto g_0_xxyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 171);

    auto g_0_xxyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 172);

    auto g_0_xxyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 173);

    auto g_0_xxyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 174);

    auto g_0_xxyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 175);

    auto g_0_xxyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 176);

    auto g_0_xxyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 177);

    auto g_0_xxyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 178);

    auto g_0_xxyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 179);

#pragma omp simd aligned(g_0_xxyyy_0_xxy_0,       \
                             g_0_xxyyy_0_xxy_1,   \
                             g_0_xxyyy_0_xyy_0,   \
                             g_0_xxyyy_0_xyy_1,   \
                             g_0_xxyyyz_0_xxy_0,  \
                             g_0_xxyyyz_0_xxy_1,  \
                             g_0_xxyyyz_0_xyy_0,  \
                             g_0_xxyyyz_0_xyy_1,  \
                             g_0_xxyyyzz_0_xxx_0, \
                             g_0_xxyyyzz_0_xxy_0, \
                             g_0_xxyyyzz_0_xxz_0, \
                             g_0_xxyyyzz_0_xyy_0, \
                             g_0_xxyyyzz_0_xyz_0, \
                             g_0_xxyyyzz_0_xzz_0, \
                             g_0_xxyyyzz_0_yyy_0, \
                             g_0_xxyyyzz_0_yyz_0, \
                             g_0_xxyyyzz_0_yzz_0, \
                             g_0_xxyyyzz_0_zzz_0, \
                             g_0_xxyyzz_0_xxx_0,  \
                             g_0_xxyyzz_0_xxx_1,  \
                             g_0_xxyyzz_0_xxz_0,  \
                             g_0_xxyyzz_0_xxz_1,  \
                             g_0_xxyyzz_0_xzz_0,  \
                             g_0_xxyyzz_0_xzz_1,  \
                             g_0_xxyzz_0_xxx_0,   \
                             g_0_xxyzz_0_xxx_1,   \
                             g_0_xxyzz_0_xxz_0,   \
                             g_0_xxyzz_0_xxz_1,   \
                             g_0_xxyzz_0_xzz_0,   \
                             g_0_xxyzz_0_xzz_1,   \
                             g_0_xyyyzz_0_xyz_0,  \
                             g_0_xyyyzz_0_xyz_1,  \
                             g_0_xyyyzz_0_yyy_0,  \
                             g_0_xyyyzz_0_yyy_1,  \
                             g_0_xyyyzz_0_yyz_0,  \
                             g_0_xyyyzz_0_yyz_1,  \
                             g_0_xyyyzz_0_yz_1,   \
                             g_0_xyyyzz_0_yzz_0,  \
                             g_0_xyyyzz_0_yzz_1,  \
                             g_0_xyyyzz_0_zzz_0,  \
                             g_0_xyyyzz_0_zzz_1,  \
                             g_0_yyyzz_0_xyz_0,   \
                             g_0_yyyzz_0_xyz_1,   \
                             g_0_yyyzz_0_yyy_0,   \
                             g_0_yyyzz_0_yyy_1,   \
                             g_0_yyyzz_0_yyz_0,   \
                             g_0_yyyzz_0_yyz_1,   \
                             g_0_yyyzz_0_yzz_0,   \
                             g_0_yyyzz_0_yzz_1,   \
                             g_0_yyyzz_0_zzz_0,   \
                             g_0_yyyzz_0_zzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzz_0_xxx_0[i] = 2.0 * g_0_xxyzz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxx_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxx_0[i] * pb_y +
                                 g_0_xxyyzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxy_0[i] =
            g_0_xxyyy_0_xxy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxy_0[i] * pb_z + g_0_xxyyyz_0_xxy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxz_0[i] = 2.0 * g_0_xxyzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxz_0[i] * pb_y +
                                 g_0_xxyyzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xyy_0[i] =
            g_0_xxyyy_0_xyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xyy_0[i] * pb_z + g_0_xxyyyz_0_xyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xyz_0[i] = g_0_yyyzz_0_xyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xyyyzz_0_xyz_0[i] * pb_x + g_0_xyyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xzz_0[i] = 2.0 * g_0_xxyzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xzz_0[i] * pb_y +
                                 g_0_xxyyzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_yyy_0[i] =
            g_0_yyyzz_0_yyy_0[i] * fi_ab_0 - g_0_yyyzz_0_yyy_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyy_0[i] * pb_x + g_0_xyyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyz_0[i] =
            g_0_yyyzz_0_yyz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyz_0[i] * pb_x + g_0_xyyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yzz_0[i] =
            g_0_yyyzz_0_yzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzz_0[i] * pb_x + g_0_xyyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_zzz_0[i] =
            g_0_yyyzz_0_zzz_0[i] * fi_ab_0 - g_0_yyyzz_0_zzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_zzz_0[i] * pb_x + g_0_xyyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 180-190 components of targeted buffer : SKSF

    auto g_0_xxyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 180);

    auto g_0_xxyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 181);

    auto g_0_xxyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 182);

    auto g_0_xxyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 183);

    auto g_0_xxyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 184);

    auto g_0_xxyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 185);

    auto g_0_xxyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 186);

    auto g_0_xxyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 187);

    auto g_0_xxyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 188);

    auto g_0_xxyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 189);

#pragma omp simd aligned(g_0_xxyyz_0_xxy_0,       \
                             g_0_xxyyz_0_xxy_1,   \
                             g_0_xxyyz_0_xyy_0,   \
                             g_0_xxyyz_0_xyy_1,   \
                             g_0_xxyyzz_0_xxy_0,  \
                             g_0_xxyyzz_0_xxy_1,  \
                             g_0_xxyyzz_0_xyy_0,  \
                             g_0_xxyyzz_0_xyy_1,  \
                             g_0_xxyyzzz_0_xxx_0, \
                             g_0_xxyyzzz_0_xxy_0, \
                             g_0_xxyyzzz_0_xxz_0, \
                             g_0_xxyyzzz_0_xyy_0, \
                             g_0_xxyyzzz_0_xyz_0, \
                             g_0_xxyyzzz_0_xzz_0, \
                             g_0_xxyyzzz_0_yyy_0, \
                             g_0_xxyyzzz_0_yyz_0, \
                             g_0_xxyyzzz_0_yzz_0, \
                             g_0_xxyyzzz_0_zzz_0, \
                             g_0_xxyzzz_0_xxx_0,  \
                             g_0_xxyzzz_0_xxx_1,  \
                             g_0_xxyzzz_0_xxz_0,  \
                             g_0_xxyzzz_0_xxz_1,  \
                             g_0_xxyzzz_0_xzz_0,  \
                             g_0_xxyzzz_0_xzz_1,  \
                             g_0_xxzzz_0_xxx_0,   \
                             g_0_xxzzz_0_xxx_1,   \
                             g_0_xxzzz_0_xxz_0,   \
                             g_0_xxzzz_0_xxz_1,   \
                             g_0_xxzzz_0_xzz_0,   \
                             g_0_xxzzz_0_xzz_1,   \
                             g_0_xyyzzz_0_xyz_0,  \
                             g_0_xyyzzz_0_xyz_1,  \
                             g_0_xyyzzz_0_yyy_0,  \
                             g_0_xyyzzz_0_yyy_1,  \
                             g_0_xyyzzz_0_yyz_0,  \
                             g_0_xyyzzz_0_yyz_1,  \
                             g_0_xyyzzz_0_yz_1,   \
                             g_0_xyyzzz_0_yzz_0,  \
                             g_0_xyyzzz_0_yzz_1,  \
                             g_0_xyyzzz_0_zzz_0,  \
                             g_0_xyyzzz_0_zzz_1,  \
                             g_0_yyzzz_0_xyz_0,   \
                             g_0_yyzzz_0_xyz_1,   \
                             g_0_yyzzz_0_yyy_0,   \
                             g_0_yyzzz_0_yyy_1,   \
                             g_0_yyzzz_0_yyz_0,   \
                             g_0_yyzzz_0_yyz_1,   \
                             g_0_yyzzz_0_yzz_0,   \
                             g_0_yyzzz_0_yzz_1,   \
                             g_0_yyzzz_0_zzz_0,   \
                             g_0_yyzzz_0_zzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzz_0_xxx_0[i] =
            g_0_xxzzz_0_xxx_0[i] * fi_ab_0 - g_0_xxzzz_0_xxx_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxx_0[i] * pb_y + g_0_xxyzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxy_0[i] = 2.0 * g_0_xxyyz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxy_0[i] * pb_z +
                                 g_0_xxyyzz_0_xxy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxz_0[i] =
            g_0_xxzzz_0_xxz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxz_0[i] * pb_y + g_0_xxyzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xyy_0[i] = 2.0 * g_0_xxyyz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xyy_0[i] * pb_z +
                                 g_0_xxyyzz_0_xyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xyz_0[i] = g_0_yyzzz_0_xyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xyyzzz_0_xyz_0[i] * pb_x + g_0_xyyzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xzz_0[i] =
            g_0_xxzzz_0_xzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xzz_0[i] * pb_y + g_0_xxyzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_yyy_0[i] =
            g_0_yyzzz_0_yyy_0[i] * fi_ab_0 - g_0_yyzzz_0_yyy_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyy_0[i] * pb_x + g_0_xyyzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyz_0[i] =
            g_0_yyzzz_0_yyz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyz_0[i] * pb_x + g_0_xyyzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yzz_0[i] =
            g_0_yyzzz_0_yzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzz_0[i] * pb_x + g_0_xyyzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_zzz_0[i] =
            g_0_yyzzz_0_zzz_0[i] * fi_ab_0 - g_0_yyzzz_0_zzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_zzz_0[i] * pb_x + g_0_xyyzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 190-200 components of targeted buffer : SKSF

    auto g_0_xxyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 190);

    auto g_0_xxyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 191);

    auto g_0_xxyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 192);

    auto g_0_xxyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 193);

    auto g_0_xxyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 194);

    auto g_0_xxyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 195);

    auto g_0_xxyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 196);

    auto g_0_xxyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 197);

    auto g_0_xxyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 198);

    auto g_0_xxyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 199);

#pragma omp simd aligned(g_0_xxyzzzz_0_xxx_0,     \
                             g_0_xxyzzzz_0_xxy_0, \
                             g_0_xxyzzzz_0_xxz_0, \
                             g_0_xxyzzzz_0_xyy_0, \
                             g_0_xxyzzzz_0_xyz_0, \
                             g_0_xxyzzzz_0_xzz_0, \
                             g_0_xxyzzzz_0_yyy_0, \
                             g_0_xxyzzzz_0_yyz_0, \
                             g_0_xxyzzzz_0_yzz_0, \
                             g_0_xxyzzzz_0_zzz_0, \
                             g_0_xxzzzz_0_xx_1,   \
                             g_0_xxzzzz_0_xxx_0,  \
                             g_0_xxzzzz_0_xxx_1,  \
                             g_0_xxzzzz_0_xxy_0,  \
                             g_0_xxzzzz_0_xxy_1,  \
                             g_0_xxzzzz_0_xxz_0,  \
                             g_0_xxzzzz_0_xxz_1,  \
                             g_0_xxzzzz_0_xy_1,   \
                             g_0_xxzzzz_0_xyy_0,  \
                             g_0_xxzzzz_0_xyy_1,  \
                             g_0_xxzzzz_0_xyz_0,  \
                             g_0_xxzzzz_0_xyz_1,  \
                             g_0_xxzzzz_0_xz_1,   \
                             g_0_xxzzzz_0_xzz_0,  \
                             g_0_xxzzzz_0_xzz_1,  \
                             g_0_xxzzzz_0_yy_1,   \
                             g_0_xxzzzz_0_yyy_0,  \
                             g_0_xxzzzz_0_yyy_1,  \
                             g_0_xxzzzz_0_yyz_0,  \
                             g_0_xxzzzz_0_yyz_1,  \
                             g_0_xxzzzz_0_yz_1,   \
                             g_0_xxzzzz_0_yzz_0,  \
                             g_0_xxzzzz_0_yzz_1,  \
                             g_0_xxzzzz_0_zz_1,   \
                             g_0_xxzzzz_0_zzz_0,  \
                             g_0_xxzzzz_0_zzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzz_0_xxx_0[i] = g_0_xxzzzz_0_xxx_0[i] * pb_y + g_0_xxzzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxy_0[i] = g_0_xxzzzz_0_xx_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxy_0[i] * pb_y + g_0_xxzzzz_0_xxy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxz_0[i] = g_0_xxzzzz_0_xxz_0[i] * pb_y + g_0_xxzzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyy_0[i] = 2.0 * g_0_xxzzzz_0_xy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyy_0[i] * pb_y + g_0_xxzzzz_0_xyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyz_0[i] = g_0_xxzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyz_0[i] * pb_y + g_0_xxzzzz_0_xyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xzz_0[i] = g_0_xxzzzz_0_xzz_0[i] * pb_y + g_0_xxzzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyy_0[i] = 3.0 * g_0_xxzzzz_0_yy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyy_0[i] * pb_y + g_0_xxzzzz_0_yyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyz_0[i] = 2.0 * g_0_xxzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyz_0[i] * pb_y + g_0_xxzzzz_0_yyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yzz_0[i] = g_0_xxzzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yzz_0[i] * pb_y + g_0_xxzzzz_0_yzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_zzz_0[i] = g_0_xxzzzz_0_zzz_0[i] * pb_y + g_0_xxzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 200-210 components of targeted buffer : SKSF

    auto g_0_xxzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 200);

    auto g_0_xxzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 201);

    auto g_0_xxzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 202);

    auto g_0_xxzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 203);

    auto g_0_xxzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 204);

    auto g_0_xxzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 205);

    auto g_0_xxzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 206);

    auto g_0_xxzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 207);

    auto g_0_xxzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 208);

    auto g_0_xxzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 209);

#pragma omp simd aligned(g_0_xxzzz_0_xxx_0,       \
                             g_0_xxzzz_0_xxx_1,   \
                             g_0_xxzzz_0_xxy_0,   \
                             g_0_xxzzz_0_xxy_1,   \
                             g_0_xxzzz_0_xyy_0,   \
                             g_0_xxzzz_0_xyy_1,   \
                             g_0_xxzzzz_0_xxx_0,  \
                             g_0_xxzzzz_0_xxx_1,  \
                             g_0_xxzzzz_0_xxy_0,  \
                             g_0_xxzzzz_0_xxy_1,  \
                             g_0_xxzzzz_0_xyy_0,  \
                             g_0_xxzzzz_0_xyy_1,  \
                             g_0_xxzzzzz_0_xxx_0, \
                             g_0_xxzzzzz_0_xxy_0, \
                             g_0_xxzzzzz_0_xxz_0, \
                             g_0_xxzzzzz_0_xyy_0, \
                             g_0_xxzzzzz_0_xyz_0, \
                             g_0_xxzzzzz_0_xzz_0, \
                             g_0_xxzzzzz_0_yyy_0, \
                             g_0_xxzzzzz_0_yyz_0, \
                             g_0_xxzzzzz_0_yzz_0, \
                             g_0_xxzzzzz_0_zzz_0, \
                             g_0_xzzzzz_0_xxz_0,  \
                             g_0_xzzzzz_0_xxz_1,  \
                             g_0_xzzzzz_0_xyz_0,  \
                             g_0_xzzzzz_0_xyz_1,  \
                             g_0_xzzzzz_0_xz_1,   \
                             g_0_xzzzzz_0_xzz_0,  \
                             g_0_xzzzzz_0_xzz_1,  \
                             g_0_xzzzzz_0_yyy_0,  \
                             g_0_xzzzzz_0_yyy_1,  \
                             g_0_xzzzzz_0_yyz_0,  \
                             g_0_xzzzzz_0_yyz_1,  \
                             g_0_xzzzzz_0_yz_1,   \
                             g_0_xzzzzz_0_yzz_0,  \
                             g_0_xzzzzz_0_yzz_1,  \
                             g_0_xzzzzz_0_zz_1,   \
                             g_0_xzzzzz_0_zzz_0,  \
                             g_0_xzzzzz_0_zzz_1,  \
                             g_0_zzzzz_0_xxz_0,   \
                             g_0_zzzzz_0_xxz_1,   \
                             g_0_zzzzz_0_xyz_0,   \
                             g_0_zzzzz_0_xyz_1,   \
                             g_0_zzzzz_0_xzz_0,   \
                             g_0_zzzzz_0_xzz_1,   \
                             g_0_zzzzz_0_yyy_0,   \
                             g_0_zzzzz_0_yyy_1,   \
                             g_0_zzzzz_0_yyz_0,   \
                             g_0_zzzzz_0_yyz_1,   \
                             g_0_zzzzz_0_yzz_0,   \
                             g_0_zzzzz_0_yzz_1,   \
                             g_0_zzzzz_0_zzz_0,   \
                             g_0_zzzzz_0_zzz_1,   \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzz_0_xxx_0[i] = 4.0 * g_0_xxzzz_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxx_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxx_0[i] * pb_z +
                                 g_0_xxzzzz_0_xxx_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxy_0[i] = 4.0 * g_0_xxzzz_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxy_0[i] * pb_z +
                                 g_0_xxzzzz_0_xxy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxz_0[i] = g_0_zzzzz_0_xxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xz_1[i] * fi_abcd_0 +
                                 g_0_xzzzzz_0_xxz_0[i] * pb_x + g_0_xzzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyy_0[i] = 4.0 * g_0_xxzzz_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xyy_0[i] * pb_z +
                                 g_0_xxzzzz_0_xyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xyz_0[i] = g_0_zzzzz_0_xyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_xzzzzz_0_xyz_0[i] * pb_x + g_0_xzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xzz_0[i] = g_0_zzzzz_0_xzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_xzzzzz_0_xzz_0[i] * pb_x + g_0_xzzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyy_0[i] =
            g_0_zzzzz_0_yyy_0[i] * fi_ab_0 - g_0_zzzzz_0_yyy_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyy_0[i] * pb_x + g_0_xzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyz_0[i] =
            g_0_zzzzz_0_yyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyz_0[i] * pb_x + g_0_xzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yzz_0[i] =
            g_0_zzzzz_0_yzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzz_0[i] * pb_x + g_0_xzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_zzz_0[i] =
            g_0_zzzzz_0_zzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzz_0[i] * pb_x + g_0_xzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 210-220 components of targeted buffer : SKSF

    auto g_0_xyyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 210);

    auto g_0_xyyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 211);

    auto g_0_xyyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 212);

    auto g_0_xyyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 213);

    auto g_0_xyyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 214);

    auto g_0_xyyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 215);

    auto g_0_xyyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 216);

    auto g_0_xyyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 217);

    auto g_0_xyyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 218);

    auto g_0_xyyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 219);

#pragma omp simd aligned(g_0_xyyyyyy_0_xxx_0,     \
                             g_0_xyyyyyy_0_xxy_0, \
                             g_0_xyyyyyy_0_xxz_0, \
                             g_0_xyyyyyy_0_xyy_0, \
                             g_0_xyyyyyy_0_xyz_0, \
                             g_0_xyyyyyy_0_xzz_0, \
                             g_0_xyyyyyy_0_yyy_0, \
                             g_0_xyyyyyy_0_yyz_0, \
                             g_0_xyyyyyy_0_yzz_0, \
                             g_0_xyyyyyy_0_zzz_0, \
                             g_0_yyyyyy_0_xx_1,   \
                             g_0_yyyyyy_0_xxx_0,  \
                             g_0_yyyyyy_0_xxx_1,  \
                             g_0_yyyyyy_0_xxy_0,  \
                             g_0_yyyyyy_0_xxy_1,  \
                             g_0_yyyyyy_0_xxz_0,  \
                             g_0_yyyyyy_0_xxz_1,  \
                             g_0_yyyyyy_0_xy_1,   \
                             g_0_yyyyyy_0_xyy_0,  \
                             g_0_yyyyyy_0_xyy_1,  \
                             g_0_yyyyyy_0_xyz_0,  \
                             g_0_yyyyyy_0_xyz_1,  \
                             g_0_yyyyyy_0_xz_1,   \
                             g_0_yyyyyy_0_xzz_0,  \
                             g_0_yyyyyy_0_xzz_1,  \
                             g_0_yyyyyy_0_yy_1,   \
                             g_0_yyyyyy_0_yyy_0,  \
                             g_0_yyyyyy_0_yyy_1,  \
                             g_0_yyyyyy_0_yyz_0,  \
                             g_0_yyyyyy_0_yyz_1,  \
                             g_0_yyyyyy_0_yz_1,   \
                             g_0_yyyyyy_0_yzz_0,  \
                             g_0_yyyyyy_0_yzz_1,  \
                             g_0_yyyyyy_0_zz_1,   \
                             g_0_yyyyyy_0_zzz_0,  \
                             g_0_yyyyyy_0_zzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyy_0_xxx_0[i] = 3.0 * g_0_yyyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxx_0[i] * pb_x + g_0_yyyyyy_0_xxx_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxy_0[i] = 2.0 * g_0_yyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxy_0[i] * pb_x + g_0_yyyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxz_0[i] = 2.0 * g_0_yyyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxz_0[i] * pb_x + g_0_yyyyyy_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyy_0[i] = g_0_yyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyy_0[i] * pb_x + g_0_yyyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyz_0[i] = g_0_yyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyz_0[i] * pb_x + g_0_yyyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xzz_0[i] = g_0_yyyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzz_0[i] * pb_x + g_0_yyyyyy_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyy_0[i] = g_0_yyyyyy_0_yyy_0[i] * pb_x + g_0_yyyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyz_0[i] = g_0_yyyyyy_0_yyz_0[i] * pb_x + g_0_yyyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yzz_0[i] = g_0_yyyyyy_0_yzz_0[i] * pb_x + g_0_yyyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_zzz_0[i] = g_0_yyyyyy_0_zzz_0[i] * pb_x + g_0_yyyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 220-230 components of targeted buffer : SKSF

    auto g_0_xyyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 220);

    auto g_0_xyyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 221);

    auto g_0_xyyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 222);

    auto g_0_xyyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 223);

    auto g_0_xyyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 224);

    auto g_0_xyyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 225);

    auto g_0_xyyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 226);

    auto g_0_xyyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 227);

    auto g_0_xyyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 228);

    auto g_0_xyyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 229);

#pragma omp simd aligned(g_0_xyyyyy_0_xxx_0,      \
                             g_0_xyyyyy_0_xxx_1,  \
                             g_0_xyyyyy_0_xxy_0,  \
                             g_0_xyyyyy_0_xxy_1,  \
                             g_0_xyyyyy_0_xyy_0,  \
                             g_0_xyyyyy_0_xyy_1,  \
                             g_0_xyyyyyz_0_xxx_0, \
                             g_0_xyyyyyz_0_xxy_0, \
                             g_0_xyyyyyz_0_xxz_0, \
                             g_0_xyyyyyz_0_xyy_0, \
                             g_0_xyyyyyz_0_xyz_0, \
                             g_0_xyyyyyz_0_xzz_0, \
                             g_0_xyyyyyz_0_yyy_0, \
                             g_0_xyyyyyz_0_yyz_0, \
                             g_0_xyyyyyz_0_yzz_0, \
                             g_0_xyyyyyz_0_zzz_0, \
                             g_0_yyyyyz_0_xxz_0,  \
                             g_0_yyyyyz_0_xxz_1,  \
                             g_0_yyyyyz_0_xyz_0,  \
                             g_0_yyyyyz_0_xyz_1,  \
                             g_0_yyyyyz_0_xz_1,   \
                             g_0_yyyyyz_0_xzz_0,  \
                             g_0_yyyyyz_0_xzz_1,  \
                             g_0_yyyyyz_0_yyy_0,  \
                             g_0_yyyyyz_0_yyy_1,  \
                             g_0_yyyyyz_0_yyz_0,  \
                             g_0_yyyyyz_0_yyz_1,  \
                             g_0_yyyyyz_0_yz_1,   \
                             g_0_yyyyyz_0_yzz_0,  \
                             g_0_yyyyyz_0_yzz_1,  \
                             g_0_yyyyyz_0_zz_1,   \
                             g_0_yyyyyz_0_zzz_0,  \
                             g_0_yyyyyz_0_zzz_1,  \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyz_0_xxx_0[i] = g_0_xyyyyy_0_xxx_0[i] * pb_z + g_0_xyyyyy_0_xxx_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxy_0[i] = g_0_xyyyyy_0_xxy_0[i] * pb_z + g_0_xyyyyy_0_xxy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxz_0[i] = 2.0 * g_0_yyyyyz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxz_0[i] * pb_x + g_0_yyyyyz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyy_0[i] = g_0_xyyyyy_0_xyy_0[i] * pb_z + g_0_xyyyyy_0_xyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xyz_0[i] = g_0_yyyyyz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyz_0[i] * pb_x + g_0_yyyyyz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xzz_0[i] = g_0_yyyyyz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xzz_0[i] * pb_x + g_0_yyyyyz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyy_0[i] = g_0_yyyyyz_0_yyy_0[i] * pb_x + g_0_yyyyyz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyz_0[i] = g_0_yyyyyz_0_yyz_0[i] * pb_x + g_0_yyyyyz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yzz_0[i] = g_0_yyyyyz_0_yzz_0[i] * pb_x + g_0_yyyyyz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_zzz_0[i] = g_0_yyyyyz_0_zzz_0[i] * pb_x + g_0_yyyyyz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 230-240 components of targeted buffer : SKSF

    auto g_0_xyyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 230);

    auto g_0_xyyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 231);

    auto g_0_xyyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 232);

    auto g_0_xyyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 233);

    auto g_0_xyyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 234);

    auto g_0_xyyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 235);

    auto g_0_xyyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 236);

    auto g_0_xyyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 237);

    auto g_0_xyyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 238);

    auto g_0_xyyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 239);

#pragma omp simd aligned(g_0_xyyyyzz_0_xxx_0,     \
                             g_0_xyyyyzz_0_xxy_0, \
                             g_0_xyyyyzz_0_xxz_0, \
                             g_0_xyyyyzz_0_xyy_0, \
                             g_0_xyyyyzz_0_xyz_0, \
                             g_0_xyyyyzz_0_xzz_0, \
                             g_0_xyyyyzz_0_yyy_0, \
                             g_0_xyyyyzz_0_yyz_0, \
                             g_0_xyyyyzz_0_yzz_0, \
                             g_0_xyyyyzz_0_zzz_0, \
                             g_0_yyyyzz_0_xx_1,   \
                             g_0_yyyyzz_0_xxx_0,  \
                             g_0_yyyyzz_0_xxx_1,  \
                             g_0_yyyyzz_0_xxy_0,  \
                             g_0_yyyyzz_0_xxy_1,  \
                             g_0_yyyyzz_0_xxz_0,  \
                             g_0_yyyyzz_0_xxz_1,  \
                             g_0_yyyyzz_0_xy_1,   \
                             g_0_yyyyzz_0_xyy_0,  \
                             g_0_yyyyzz_0_xyy_1,  \
                             g_0_yyyyzz_0_xyz_0,  \
                             g_0_yyyyzz_0_xyz_1,  \
                             g_0_yyyyzz_0_xz_1,   \
                             g_0_yyyyzz_0_xzz_0,  \
                             g_0_yyyyzz_0_xzz_1,  \
                             g_0_yyyyzz_0_yy_1,   \
                             g_0_yyyyzz_0_yyy_0,  \
                             g_0_yyyyzz_0_yyy_1,  \
                             g_0_yyyyzz_0_yyz_0,  \
                             g_0_yyyyzz_0_yyz_1,  \
                             g_0_yyyyzz_0_yz_1,   \
                             g_0_yyyyzz_0_yzz_0,  \
                             g_0_yyyyzz_0_yzz_1,  \
                             g_0_yyyyzz_0_zz_1,   \
                             g_0_yyyyzz_0_zzz_0,  \
                             g_0_yyyyzz_0_zzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzz_0_xxx_0[i] = 3.0 * g_0_yyyyzz_0_xx_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxx_0[i] * pb_x + g_0_yyyyzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxy_0[i] = 2.0 * g_0_yyyyzz_0_xy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxy_0[i] * pb_x + g_0_yyyyzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxz_0[i] = 2.0 * g_0_yyyyzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxz_0[i] * pb_x + g_0_yyyyzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyy_0[i] = g_0_yyyyzz_0_yy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyy_0[i] * pb_x + g_0_yyyyzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyz_0[i] = g_0_yyyyzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyz_0[i] * pb_x + g_0_yyyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xzz_0[i] = g_0_yyyyzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xzz_0[i] * pb_x + g_0_yyyyzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyy_0[i] = g_0_yyyyzz_0_yyy_0[i] * pb_x + g_0_yyyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyz_0[i] = g_0_yyyyzz_0_yyz_0[i] * pb_x + g_0_yyyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yzz_0[i] = g_0_yyyyzz_0_yzz_0[i] * pb_x + g_0_yyyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_zzz_0[i] = g_0_yyyyzz_0_zzz_0[i] * pb_x + g_0_yyyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 240-250 components of targeted buffer : SKSF

    auto g_0_xyyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 240);

    auto g_0_xyyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 241);

    auto g_0_xyyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 242);

    auto g_0_xyyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 243);

    auto g_0_xyyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 244);

    auto g_0_xyyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 245);

    auto g_0_xyyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 246);

    auto g_0_xyyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 247);

    auto g_0_xyyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 248);

    auto g_0_xyyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 249);

#pragma omp simd aligned(g_0_xyyyzzz_0_xxx_0,     \
                             g_0_xyyyzzz_0_xxy_0, \
                             g_0_xyyyzzz_0_xxz_0, \
                             g_0_xyyyzzz_0_xyy_0, \
                             g_0_xyyyzzz_0_xyz_0, \
                             g_0_xyyyzzz_0_xzz_0, \
                             g_0_xyyyzzz_0_yyy_0, \
                             g_0_xyyyzzz_0_yyz_0, \
                             g_0_xyyyzzz_0_yzz_0, \
                             g_0_xyyyzzz_0_zzz_0, \
                             g_0_yyyzzz_0_xx_1,   \
                             g_0_yyyzzz_0_xxx_0,  \
                             g_0_yyyzzz_0_xxx_1,  \
                             g_0_yyyzzz_0_xxy_0,  \
                             g_0_yyyzzz_0_xxy_1,  \
                             g_0_yyyzzz_0_xxz_0,  \
                             g_0_yyyzzz_0_xxz_1,  \
                             g_0_yyyzzz_0_xy_1,   \
                             g_0_yyyzzz_0_xyy_0,  \
                             g_0_yyyzzz_0_xyy_1,  \
                             g_0_yyyzzz_0_xyz_0,  \
                             g_0_yyyzzz_0_xyz_1,  \
                             g_0_yyyzzz_0_xz_1,   \
                             g_0_yyyzzz_0_xzz_0,  \
                             g_0_yyyzzz_0_xzz_1,  \
                             g_0_yyyzzz_0_yy_1,   \
                             g_0_yyyzzz_0_yyy_0,  \
                             g_0_yyyzzz_0_yyy_1,  \
                             g_0_yyyzzz_0_yyz_0,  \
                             g_0_yyyzzz_0_yyz_1,  \
                             g_0_yyyzzz_0_yz_1,   \
                             g_0_yyyzzz_0_yzz_0,  \
                             g_0_yyyzzz_0_yzz_1,  \
                             g_0_yyyzzz_0_zz_1,   \
                             g_0_yyyzzz_0_zzz_0,  \
                             g_0_yyyzzz_0_zzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzz_0_xxx_0[i] = 3.0 * g_0_yyyzzz_0_xx_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxx_0[i] * pb_x + g_0_yyyzzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxy_0[i] = 2.0 * g_0_yyyzzz_0_xy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxy_0[i] * pb_x + g_0_yyyzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxz_0[i] = 2.0 * g_0_yyyzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxz_0[i] * pb_x + g_0_yyyzzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyy_0[i] = g_0_yyyzzz_0_yy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyy_0[i] * pb_x + g_0_yyyzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyz_0[i] = g_0_yyyzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyz_0[i] * pb_x + g_0_yyyzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xzz_0[i] = g_0_yyyzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xzz_0[i] * pb_x + g_0_yyyzzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyy_0[i] = g_0_yyyzzz_0_yyy_0[i] * pb_x + g_0_yyyzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyz_0[i] = g_0_yyyzzz_0_yyz_0[i] * pb_x + g_0_yyyzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yzz_0[i] = g_0_yyyzzz_0_yzz_0[i] * pb_x + g_0_yyyzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_zzz_0[i] = g_0_yyyzzz_0_zzz_0[i] * pb_x + g_0_yyyzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 250-260 components of targeted buffer : SKSF

    auto g_0_xyyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 250);

    auto g_0_xyyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 251);

    auto g_0_xyyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 252);

    auto g_0_xyyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 253);

    auto g_0_xyyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 254);

    auto g_0_xyyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 255);

    auto g_0_xyyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 256);

    auto g_0_xyyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 257);

    auto g_0_xyyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 258);

    auto g_0_xyyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 259);

#pragma omp simd aligned(g_0_xyyzzzz_0_xxx_0,     \
                             g_0_xyyzzzz_0_xxy_0, \
                             g_0_xyyzzzz_0_xxz_0, \
                             g_0_xyyzzzz_0_xyy_0, \
                             g_0_xyyzzzz_0_xyz_0, \
                             g_0_xyyzzzz_0_xzz_0, \
                             g_0_xyyzzzz_0_yyy_0, \
                             g_0_xyyzzzz_0_yyz_0, \
                             g_0_xyyzzzz_0_yzz_0, \
                             g_0_xyyzzzz_0_zzz_0, \
                             g_0_yyzzzz_0_xx_1,   \
                             g_0_yyzzzz_0_xxx_0,  \
                             g_0_yyzzzz_0_xxx_1,  \
                             g_0_yyzzzz_0_xxy_0,  \
                             g_0_yyzzzz_0_xxy_1,  \
                             g_0_yyzzzz_0_xxz_0,  \
                             g_0_yyzzzz_0_xxz_1,  \
                             g_0_yyzzzz_0_xy_1,   \
                             g_0_yyzzzz_0_xyy_0,  \
                             g_0_yyzzzz_0_xyy_1,  \
                             g_0_yyzzzz_0_xyz_0,  \
                             g_0_yyzzzz_0_xyz_1,  \
                             g_0_yyzzzz_0_xz_1,   \
                             g_0_yyzzzz_0_xzz_0,  \
                             g_0_yyzzzz_0_xzz_1,  \
                             g_0_yyzzzz_0_yy_1,   \
                             g_0_yyzzzz_0_yyy_0,  \
                             g_0_yyzzzz_0_yyy_1,  \
                             g_0_yyzzzz_0_yyz_0,  \
                             g_0_yyzzzz_0_yyz_1,  \
                             g_0_yyzzzz_0_yz_1,   \
                             g_0_yyzzzz_0_yzz_0,  \
                             g_0_yyzzzz_0_yzz_1,  \
                             g_0_yyzzzz_0_zz_1,   \
                             g_0_yyzzzz_0_zzz_0,  \
                             g_0_yyzzzz_0_zzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzz_0_xxx_0[i] = 3.0 * g_0_yyzzzz_0_xx_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxx_0[i] * pb_x + g_0_yyzzzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxy_0[i] = 2.0 * g_0_yyzzzz_0_xy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxy_0[i] * pb_x + g_0_yyzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxz_0[i] = 2.0 * g_0_yyzzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxz_0[i] * pb_x + g_0_yyzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyy_0[i] = g_0_yyzzzz_0_yy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyy_0[i] * pb_x + g_0_yyzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyz_0[i] = g_0_yyzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyz_0[i] * pb_x + g_0_yyzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xzz_0[i] = g_0_yyzzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xzz_0[i] * pb_x + g_0_yyzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyy_0[i] = g_0_yyzzzz_0_yyy_0[i] * pb_x + g_0_yyzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyz_0[i] = g_0_yyzzzz_0_yyz_0[i] * pb_x + g_0_yyzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yzz_0[i] = g_0_yyzzzz_0_yzz_0[i] * pb_x + g_0_yyzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_zzz_0[i] = g_0_yyzzzz_0_zzz_0[i] * pb_x + g_0_yyzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 260-270 components of targeted buffer : SKSF

    auto g_0_xyzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 260);

    auto g_0_xyzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 261);

    auto g_0_xyzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 262);

    auto g_0_xyzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 263);

    auto g_0_xyzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 264);

    auto g_0_xyzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 265);

    auto g_0_xyzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 266);

    auto g_0_xyzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 267);

    auto g_0_xyzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 268);

    auto g_0_xyzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 269);

#pragma omp simd aligned(g_0_xyzzzzz_0_xxx_0,     \
                             g_0_xyzzzzz_0_xxy_0, \
                             g_0_xyzzzzz_0_xxz_0, \
                             g_0_xyzzzzz_0_xyy_0, \
                             g_0_xyzzzzz_0_xyz_0, \
                             g_0_xyzzzzz_0_xzz_0, \
                             g_0_xyzzzzz_0_yyy_0, \
                             g_0_xyzzzzz_0_yyz_0, \
                             g_0_xyzzzzz_0_yzz_0, \
                             g_0_xyzzzzz_0_zzz_0, \
                             g_0_xzzzzz_0_xxx_0,  \
                             g_0_xzzzzz_0_xxx_1,  \
                             g_0_xzzzzz_0_xxz_0,  \
                             g_0_xzzzzz_0_xxz_1,  \
                             g_0_xzzzzz_0_xzz_0,  \
                             g_0_xzzzzz_0_xzz_1,  \
                             g_0_yzzzzz_0_xxy_0,  \
                             g_0_yzzzzz_0_xxy_1,  \
                             g_0_yzzzzz_0_xy_1,   \
                             g_0_yzzzzz_0_xyy_0,  \
                             g_0_yzzzzz_0_xyy_1,  \
                             g_0_yzzzzz_0_xyz_0,  \
                             g_0_yzzzzz_0_xyz_1,  \
                             g_0_yzzzzz_0_yy_1,   \
                             g_0_yzzzzz_0_yyy_0,  \
                             g_0_yzzzzz_0_yyy_1,  \
                             g_0_yzzzzz_0_yyz_0,  \
                             g_0_yzzzzz_0_yyz_1,  \
                             g_0_yzzzzz_0_yz_1,   \
                             g_0_yzzzzz_0_yzz_0,  \
                             g_0_yzzzzz_0_yzz_1,  \
                             g_0_yzzzzz_0_zzz_0,  \
                             g_0_yzzzzz_0_zzz_1,  \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzz_0_xxx_0[i] = g_0_xzzzzz_0_xxx_0[i] * pb_y + g_0_xzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxy_0[i] = 2.0 * g_0_yzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxy_0[i] * pb_x + g_0_yzzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxz_0[i] = g_0_xzzzzz_0_xxz_0[i] * pb_y + g_0_xzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xyy_0[i] = g_0_yzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyy_0[i] * pb_x + g_0_yzzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyz_0[i] = g_0_yzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyz_0[i] * pb_x + g_0_yzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xzz_0[i] = g_0_xzzzzz_0_xzz_0[i] * pb_y + g_0_xzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_yyy_0[i] = g_0_yzzzzz_0_yyy_0[i] * pb_x + g_0_yzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyz_0[i] = g_0_yzzzzz_0_yyz_0[i] * pb_x + g_0_yzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yzz_0[i] = g_0_yzzzzz_0_yzz_0[i] * pb_x + g_0_yzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_zzz_0[i] = g_0_yzzzzz_0_zzz_0[i] * pb_x + g_0_yzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 270-280 components of targeted buffer : SKSF

    auto g_0_xzzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 270);

    auto g_0_xzzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 271);

    auto g_0_xzzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 272);

    auto g_0_xzzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 273);

    auto g_0_xzzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 274);

    auto g_0_xzzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 275);

    auto g_0_xzzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 276);

    auto g_0_xzzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 277);

    auto g_0_xzzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 278);

    auto g_0_xzzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 279);

#pragma omp simd aligned(g_0_xzzzzzz_0_xxx_0,     \
                             g_0_xzzzzzz_0_xxy_0, \
                             g_0_xzzzzzz_0_xxz_0, \
                             g_0_xzzzzzz_0_xyy_0, \
                             g_0_xzzzzzz_0_xyz_0, \
                             g_0_xzzzzzz_0_xzz_0, \
                             g_0_xzzzzzz_0_yyy_0, \
                             g_0_xzzzzzz_0_yyz_0, \
                             g_0_xzzzzzz_0_yzz_0, \
                             g_0_xzzzzzz_0_zzz_0, \
                             g_0_zzzzzz_0_xx_1,   \
                             g_0_zzzzzz_0_xxx_0,  \
                             g_0_zzzzzz_0_xxx_1,  \
                             g_0_zzzzzz_0_xxy_0,  \
                             g_0_zzzzzz_0_xxy_1,  \
                             g_0_zzzzzz_0_xxz_0,  \
                             g_0_zzzzzz_0_xxz_1,  \
                             g_0_zzzzzz_0_xy_1,   \
                             g_0_zzzzzz_0_xyy_0,  \
                             g_0_zzzzzz_0_xyy_1,  \
                             g_0_zzzzzz_0_xyz_0,  \
                             g_0_zzzzzz_0_xyz_1,  \
                             g_0_zzzzzz_0_xz_1,   \
                             g_0_zzzzzz_0_xzz_0,  \
                             g_0_zzzzzz_0_xzz_1,  \
                             g_0_zzzzzz_0_yy_1,   \
                             g_0_zzzzzz_0_yyy_0,  \
                             g_0_zzzzzz_0_yyy_1,  \
                             g_0_zzzzzz_0_yyz_0,  \
                             g_0_zzzzzz_0_yyz_1,  \
                             g_0_zzzzzz_0_yz_1,   \
                             g_0_zzzzzz_0_yzz_0,  \
                             g_0_zzzzzz_0_yzz_1,  \
                             g_0_zzzzzz_0_zz_1,   \
                             g_0_zzzzzz_0_zzz_0,  \
                             g_0_zzzzzz_0_zzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzz_0_xxx_0[i] = 3.0 * g_0_zzzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxx_0[i] * pb_x + g_0_zzzzzz_0_xxx_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxy_0[i] = 2.0 * g_0_zzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxy_0[i] * pb_x + g_0_zzzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxz_0[i] = 2.0 * g_0_zzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxz_0[i] * pb_x + g_0_zzzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyy_0[i] = g_0_zzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyy_0[i] * pb_x + g_0_zzzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyz_0[i] = g_0_zzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyz_0[i] * pb_x + g_0_zzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xzz_0[i] = g_0_zzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzz_0[i] * pb_x + g_0_zzzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyy_0[i] = g_0_zzzzzz_0_yyy_0[i] * pb_x + g_0_zzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyz_0[i] = g_0_zzzzzz_0_yyz_0[i] * pb_x + g_0_zzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yzz_0[i] = g_0_zzzzzz_0_yzz_0[i] * pb_x + g_0_zzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_zzz_0[i] = g_0_zzzzzz_0_zzz_0[i] * pb_x + g_0_zzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 280-290 components of targeted buffer : SKSF

    auto g_0_yyyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 280);

    auto g_0_yyyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 281);

    auto g_0_yyyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 282);

    auto g_0_yyyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 283);

    auto g_0_yyyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 284);

    auto g_0_yyyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 285);

    auto g_0_yyyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 286);

    auto g_0_yyyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 287);

    auto g_0_yyyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 288);

    auto g_0_yyyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 289);

#pragma omp simd aligned(g_0_yyyyy_0_xxx_0,       \
                             g_0_yyyyy_0_xxx_1,   \
                             g_0_yyyyy_0_xxy_0,   \
                             g_0_yyyyy_0_xxy_1,   \
                             g_0_yyyyy_0_xxz_0,   \
                             g_0_yyyyy_0_xxz_1,   \
                             g_0_yyyyy_0_xyy_0,   \
                             g_0_yyyyy_0_xyy_1,   \
                             g_0_yyyyy_0_xyz_0,   \
                             g_0_yyyyy_0_xyz_1,   \
                             g_0_yyyyy_0_xzz_0,   \
                             g_0_yyyyy_0_xzz_1,   \
                             g_0_yyyyy_0_yyy_0,   \
                             g_0_yyyyy_0_yyy_1,   \
                             g_0_yyyyy_0_yyz_0,   \
                             g_0_yyyyy_0_yyz_1,   \
                             g_0_yyyyy_0_yzz_0,   \
                             g_0_yyyyy_0_yzz_1,   \
                             g_0_yyyyy_0_zzz_0,   \
                             g_0_yyyyy_0_zzz_1,   \
                             g_0_yyyyyy_0_xx_1,   \
                             g_0_yyyyyy_0_xxx_0,  \
                             g_0_yyyyyy_0_xxx_1,  \
                             g_0_yyyyyy_0_xxy_0,  \
                             g_0_yyyyyy_0_xxy_1,  \
                             g_0_yyyyyy_0_xxz_0,  \
                             g_0_yyyyyy_0_xxz_1,  \
                             g_0_yyyyyy_0_xy_1,   \
                             g_0_yyyyyy_0_xyy_0,  \
                             g_0_yyyyyy_0_xyy_1,  \
                             g_0_yyyyyy_0_xyz_0,  \
                             g_0_yyyyyy_0_xyz_1,  \
                             g_0_yyyyyy_0_xz_1,   \
                             g_0_yyyyyy_0_xzz_0,  \
                             g_0_yyyyyy_0_xzz_1,  \
                             g_0_yyyyyy_0_yy_1,   \
                             g_0_yyyyyy_0_yyy_0,  \
                             g_0_yyyyyy_0_yyy_1,  \
                             g_0_yyyyyy_0_yyz_0,  \
                             g_0_yyyyyy_0_yyz_1,  \
                             g_0_yyyyyy_0_yz_1,   \
                             g_0_yyyyyy_0_yzz_0,  \
                             g_0_yyyyyy_0_yzz_1,  \
                             g_0_yyyyyy_0_zz_1,   \
                             g_0_yyyyyy_0_zzz_0,  \
                             g_0_yyyyyy_0_zzz_1,  \
                             g_0_yyyyyyy_0_xxx_0, \
                             g_0_yyyyyyy_0_xxy_0, \
                             g_0_yyyyyyy_0_xxz_0, \
                             g_0_yyyyyyy_0_xyy_0, \
                             g_0_yyyyyyy_0_xyz_0, \
                             g_0_yyyyyyy_0_xzz_0, \
                             g_0_yyyyyyy_0_yyy_0, \
                             g_0_yyyyyyy_0_yyz_0, \
                             g_0_yyyyyyy_0_yzz_0, \
                             g_0_yyyyyyy_0_zzz_0, \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyy_0_xxx_0[i] = 6.0 * g_0_yyyyy_0_xxx_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxx_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxx_0[i] * pb_y +
                                 g_0_yyyyyy_0_xxx_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxy_0[i] = 6.0 * g_0_yyyyy_0_xxy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyyyy_0_xx_1[i] * fi_abcd_0 +
                                 g_0_yyyyyy_0_xxy_0[i] * pb_y + g_0_yyyyyy_0_xxy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxz_0[i] = 6.0 * g_0_yyyyy_0_xxz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxz_0[i] * pb_y +
                                 g_0_yyyyyy_0_xxz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyy_0[i] = 6.0 * g_0_yyyyy_0_xyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyy_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyy_0[i] * pb_y + g_0_yyyyyy_0_xyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyz_0[i] = 6.0 * g_0_yyyyy_0_xyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xz_1[i] * fi_abcd_0 +
                                 g_0_yyyyyy_0_xyz_0[i] * pb_y + g_0_yyyyyy_0_xyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xzz_0[i] = 6.0 * g_0_yyyyy_0_xzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xzz_0[i] * pb_y +
                                 g_0_yyyyyy_0_xzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyy_0[i] = 6.0 * g_0_yyyyy_0_yyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyy_1[i] * fti_ab_0 +
                                 3.0 * g_0_yyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyy_0[i] * pb_y + g_0_yyyyyy_0_yyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyz_0[i] = 6.0 * g_0_yyyyy_0_yyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyz_0[i] * pb_y + g_0_yyyyyy_0_yyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yzz_0[i] = 6.0 * g_0_yyyyy_0_yzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zz_1[i] * fi_abcd_0 +
                                 g_0_yyyyyy_0_yzz_0[i] * pb_y + g_0_yyyyyy_0_yzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_zzz_0[i] = 6.0 * g_0_yyyyy_0_zzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_zzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zzz_0[i] * pb_y +
                                 g_0_yyyyyy_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 290-300 components of targeted buffer : SKSF

    auto g_0_yyyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 290);

    auto g_0_yyyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 291);

    auto g_0_yyyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 292);

    auto g_0_yyyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 293);

    auto g_0_yyyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 294);

    auto g_0_yyyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 295);

    auto g_0_yyyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 296);

    auto g_0_yyyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 297);

    auto g_0_yyyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 298);

    auto g_0_yyyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 299);

#pragma omp simd aligned(g_0_yyyyyy_0_xx_1,       \
                             g_0_yyyyyy_0_xxx_0,  \
                             g_0_yyyyyy_0_xxx_1,  \
                             g_0_yyyyyy_0_xxy_0,  \
                             g_0_yyyyyy_0_xxy_1,  \
                             g_0_yyyyyy_0_xxz_0,  \
                             g_0_yyyyyy_0_xxz_1,  \
                             g_0_yyyyyy_0_xy_1,   \
                             g_0_yyyyyy_0_xyy_0,  \
                             g_0_yyyyyy_0_xyy_1,  \
                             g_0_yyyyyy_0_xyz_0,  \
                             g_0_yyyyyy_0_xyz_1,  \
                             g_0_yyyyyy_0_xz_1,   \
                             g_0_yyyyyy_0_xzz_0,  \
                             g_0_yyyyyy_0_xzz_1,  \
                             g_0_yyyyyy_0_yy_1,   \
                             g_0_yyyyyy_0_yyy_0,  \
                             g_0_yyyyyy_0_yyy_1,  \
                             g_0_yyyyyy_0_yyz_0,  \
                             g_0_yyyyyy_0_yyz_1,  \
                             g_0_yyyyyy_0_yz_1,   \
                             g_0_yyyyyy_0_yzz_0,  \
                             g_0_yyyyyy_0_yzz_1,  \
                             g_0_yyyyyy_0_zz_1,   \
                             g_0_yyyyyy_0_zzz_0,  \
                             g_0_yyyyyy_0_zzz_1,  \
                             g_0_yyyyyyz_0_xxx_0, \
                             g_0_yyyyyyz_0_xxy_0, \
                             g_0_yyyyyyz_0_xxz_0, \
                             g_0_yyyyyyz_0_xyy_0, \
                             g_0_yyyyyyz_0_xyz_0, \
                             g_0_yyyyyyz_0_xzz_0, \
                             g_0_yyyyyyz_0_yyy_0, \
                             g_0_yyyyyyz_0_yyz_0, \
                             g_0_yyyyyyz_0_yzz_0, \
                             g_0_yyyyyyz_0_zzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyz_0_xxx_0[i] = g_0_yyyyyy_0_xxx_0[i] * pb_z + g_0_yyyyyy_0_xxx_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxy_0[i] = g_0_yyyyyy_0_xxy_0[i] * pb_z + g_0_yyyyyy_0_xxy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxz_0[i] = g_0_yyyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxz_0[i] * pb_z + g_0_yyyyyy_0_xxz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyy_0[i] = g_0_yyyyyy_0_xyy_0[i] * pb_z + g_0_yyyyyy_0_xyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyz_0[i] = g_0_yyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyz_0[i] * pb_z + g_0_yyyyyy_0_xyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xzz_0[i] = 2.0 * g_0_yyyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzz_0[i] * pb_z + g_0_yyyyyy_0_xzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyy_0[i] = g_0_yyyyyy_0_yyy_0[i] * pb_z + g_0_yyyyyy_0_yyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyz_0[i] = g_0_yyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyz_0[i] * pb_z + g_0_yyyyyy_0_yyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yzz_0[i] = 2.0 * g_0_yyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzz_0[i] * pb_z + g_0_yyyyyy_0_yzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_zzz_0[i] = 3.0 * g_0_yyyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_zzz_0[i] * pb_z + g_0_yyyyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 300-310 components of targeted buffer : SKSF

    auto g_0_yyyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 300);

    auto g_0_yyyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 301);

    auto g_0_yyyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 302);

    auto g_0_yyyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 303);

    auto g_0_yyyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 304);

    auto g_0_yyyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 305);

    auto g_0_yyyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 306);

    auto g_0_yyyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 307);

    auto g_0_yyyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 308);

    auto g_0_yyyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 309);

#pragma omp simd aligned(g_0_yyyyy_0_xxy_0,       \
                             g_0_yyyyy_0_xxy_1,   \
                             g_0_yyyyy_0_xyy_0,   \
                             g_0_yyyyy_0_xyy_1,   \
                             g_0_yyyyy_0_yyy_0,   \
                             g_0_yyyyy_0_yyy_1,   \
                             g_0_yyyyyz_0_xxy_0,  \
                             g_0_yyyyyz_0_xxy_1,  \
                             g_0_yyyyyz_0_xyy_0,  \
                             g_0_yyyyyz_0_xyy_1,  \
                             g_0_yyyyyz_0_yyy_0,  \
                             g_0_yyyyyz_0_yyy_1,  \
                             g_0_yyyyyzz_0_xxx_0, \
                             g_0_yyyyyzz_0_xxy_0, \
                             g_0_yyyyyzz_0_xxz_0, \
                             g_0_yyyyyzz_0_xyy_0, \
                             g_0_yyyyyzz_0_xyz_0, \
                             g_0_yyyyyzz_0_xzz_0, \
                             g_0_yyyyyzz_0_yyy_0, \
                             g_0_yyyyyzz_0_yyz_0, \
                             g_0_yyyyyzz_0_yzz_0, \
                             g_0_yyyyyzz_0_zzz_0, \
                             g_0_yyyyzz_0_xxx_0,  \
                             g_0_yyyyzz_0_xxx_1,  \
                             g_0_yyyyzz_0_xxz_0,  \
                             g_0_yyyyzz_0_xxz_1,  \
                             g_0_yyyyzz_0_xyz_0,  \
                             g_0_yyyyzz_0_xyz_1,  \
                             g_0_yyyyzz_0_xz_1,   \
                             g_0_yyyyzz_0_xzz_0,  \
                             g_0_yyyyzz_0_xzz_1,  \
                             g_0_yyyyzz_0_yyz_0,  \
                             g_0_yyyyzz_0_yyz_1,  \
                             g_0_yyyyzz_0_yz_1,   \
                             g_0_yyyyzz_0_yzz_0,  \
                             g_0_yyyyzz_0_yzz_1,  \
                             g_0_yyyyzz_0_zz_1,   \
                             g_0_yyyyzz_0_zzz_0,  \
                             g_0_yyyyzz_0_zzz_1,  \
                             g_0_yyyzz_0_xxx_0,   \
                             g_0_yyyzz_0_xxx_1,   \
                             g_0_yyyzz_0_xxz_0,   \
                             g_0_yyyzz_0_xxz_1,   \
                             g_0_yyyzz_0_xyz_0,   \
                             g_0_yyyzz_0_xyz_1,   \
                             g_0_yyyzz_0_xzz_0,   \
                             g_0_yyyzz_0_xzz_1,   \
                             g_0_yyyzz_0_yyz_0,   \
                             g_0_yyyzz_0_yyz_1,   \
                             g_0_yyyzz_0_yzz_0,   \
                             g_0_yyyzz_0_yzz_1,   \
                             g_0_yyyzz_0_zzz_0,   \
                             g_0_yyyzz_0_zzz_1,   \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzz_0_xxx_0[i] = 4.0 * g_0_yyyzz_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxx_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxx_0[i] * pb_y +
                                 g_0_yyyyzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxy_0[i] =
            g_0_yyyyy_0_xxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxy_0[i] * pb_z + g_0_yyyyyz_0_xxy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxz_0[i] = 4.0 * g_0_yyyzz_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxz_0[i] * pb_y +
                                 g_0_yyyyzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyy_0[i] =
            g_0_yyyyy_0_xyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xyy_0[i] * pb_z + g_0_yyyyyz_0_xyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xyz_0[i] = 4.0 * g_0_yyyzz_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xz_1[i] * fi_abcd_0 +
                                 g_0_yyyyzz_0_xyz_0[i] * pb_y + g_0_yyyyzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xzz_0[i] = 4.0 * g_0_yyyzz_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xzz_0[i] * pb_y +
                                 g_0_yyyyzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyy_0[i] =
            g_0_yyyyy_0_yyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_yyy_0[i] * pb_z + g_0_yyyyyz_0_yyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_yyz_0[i] = 4.0 * g_0_yyyzz_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyyzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyz_0[i] * pb_y + g_0_yyyyzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yzz_0[i] = 4.0 * g_0_yyyzz_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_yyyyzz_0_yzz_0[i] * pb_y + g_0_yyyyzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_zzz_0[i] = 4.0 * g_0_yyyzz_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_zzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zzz_0[i] * pb_y +
                                 g_0_yyyyzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 310-320 components of targeted buffer : SKSF

    auto g_0_yyyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 310);

    auto g_0_yyyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 311);

    auto g_0_yyyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 312);

    auto g_0_yyyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 313);

    auto g_0_yyyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 314);

    auto g_0_yyyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 315);

    auto g_0_yyyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 316);

    auto g_0_yyyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 317);

    auto g_0_yyyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 318);

    auto g_0_yyyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 319);

#pragma omp simd aligned(g_0_yyyyz_0_xxy_0,       \
                             g_0_yyyyz_0_xxy_1,   \
                             g_0_yyyyz_0_xyy_0,   \
                             g_0_yyyyz_0_xyy_1,   \
                             g_0_yyyyz_0_yyy_0,   \
                             g_0_yyyyz_0_yyy_1,   \
                             g_0_yyyyzz_0_xxy_0,  \
                             g_0_yyyyzz_0_xxy_1,  \
                             g_0_yyyyzz_0_xyy_0,  \
                             g_0_yyyyzz_0_xyy_1,  \
                             g_0_yyyyzz_0_yyy_0,  \
                             g_0_yyyyzz_0_yyy_1,  \
                             g_0_yyyyzzz_0_xxx_0, \
                             g_0_yyyyzzz_0_xxy_0, \
                             g_0_yyyyzzz_0_xxz_0, \
                             g_0_yyyyzzz_0_xyy_0, \
                             g_0_yyyyzzz_0_xyz_0, \
                             g_0_yyyyzzz_0_xzz_0, \
                             g_0_yyyyzzz_0_yyy_0, \
                             g_0_yyyyzzz_0_yyz_0, \
                             g_0_yyyyzzz_0_yzz_0, \
                             g_0_yyyyzzz_0_zzz_0, \
                             g_0_yyyzzz_0_xxx_0,  \
                             g_0_yyyzzz_0_xxx_1,  \
                             g_0_yyyzzz_0_xxz_0,  \
                             g_0_yyyzzz_0_xxz_1,  \
                             g_0_yyyzzz_0_xyz_0,  \
                             g_0_yyyzzz_0_xyz_1,  \
                             g_0_yyyzzz_0_xz_1,   \
                             g_0_yyyzzz_0_xzz_0,  \
                             g_0_yyyzzz_0_xzz_1,  \
                             g_0_yyyzzz_0_yyz_0,  \
                             g_0_yyyzzz_0_yyz_1,  \
                             g_0_yyyzzz_0_yz_1,   \
                             g_0_yyyzzz_0_yzz_0,  \
                             g_0_yyyzzz_0_yzz_1,  \
                             g_0_yyyzzz_0_zz_1,   \
                             g_0_yyyzzz_0_zzz_0,  \
                             g_0_yyyzzz_0_zzz_1,  \
                             g_0_yyzzz_0_xxx_0,   \
                             g_0_yyzzz_0_xxx_1,   \
                             g_0_yyzzz_0_xxz_0,   \
                             g_0_yyzzz_0_xxz_1,   \
                             g_0_yyzzz_0_xyz_0,   \
                             g_0_yyzzz_0_xyz_1,   \
                             g_0_yyzzz_0_xzz_0,   \
                             g_0_yyzzz_0_xzz_1,   \
                             g_0_yyzzz_0_yyz_0,   \
                             g_0_yyzzz_0_yyz_1,   \
                             g_0_yyzzz_0_yzz_0,   \
                             g_0_yyzzz_0_yzz_1,   \
                             g_0_yyzzz_0_zzz_0,   \
                             g_0_yyzzz_0_zzz_1,   \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzz_0_xxx_0[i] = 3.0 * g_0_yyzzz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxx_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxx_0[i] * pb_y +
                                 g_0_yyyzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxy_0[i] = 2.0 * g_0_yyyyz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxy_0[i] * pb_z +
                                 g_0_yyyyzz_0_xxy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxz_0[i] = 3.0 * g_0_yyzzz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxz_0[i] * pb_y +
                                 g_0_yyyzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyy_0[i] = 2.0 * g_0_yyyyz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xyy_0[i] * pb_z +
                                 g_0_yyyyzz_0_xyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xyz_0[i] = 3.0 * g_0_yyzzz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xz_1[i] * fi_abcd_0 +
                                 g_0_yyyzzz_0_xyz_0[i] * pb_y + g_0_yyyzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xzz_0[i] = 3.0 * g_0_yyzzz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xzz_0[i] * pb_y +
                                 g_0_yyyzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyy_0[i] = 2.0 * g_0_yyyyz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_yyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_yyy_0[i] * pb_z +
                                 g_0_yyyyzz_0_yyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_yyz_0[i] = 3.0 * g_0_yyzzz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyz_0[i] * pb_y + g_0_yyyzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yzz_0[i] = 3.0 * g_0_yyzzz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_yyyzzz_0_yzz_0[i] * pb_y + g_0_yyyzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_zzz_0[i] = 3.0 * g_0_yyzzz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_zzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zzz_0[i] * pb_y +
                                 g_0_yyyzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 320-330 components of targeted buffer : SKSF

    auto g_0_yyyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 320);

    auto g_0_yyyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 321);

    auto g_0_yyyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 322);

    auto g_0_yyyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 323);

    auto g_0_yyyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 324);

    auto g_0_yyyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 325);

    auto g_0_yyyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 326);

    auto g_0_yyyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 327);

    auto g_0_yyyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 328);

    auto g_0_yyyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 329);

#pragma omp simd aligned(g_0_yyyzz_0_xxy_0,       \
                             g_0_yyyzz_0_xxy_1,   \
                             g_0_yyyzz_0_xyy_0,   \
                             g_0_yyyzz_0_xyy_1,   \
                             g_0_yyyzz_0_yyy_0,   \
                             g_0_yyyzz_0_yyy_1,   \
                             g_0_yyyzzz_0_xxy_0,  \
                             g_0_yyyzzz_0_xxy_1,  \
                             g_0_yyyzzz_0_xyy_0,  \
                             g_0_yyyzzz_0_xyy_1,  \
                             g_0_yyyzzz_0_yyy_0,  \
                             g_0_yyyzzz_0_yyy_1,  \
                             g_0_yyyzzzz_0_xxx_0, \
                             g_0_yyyzzzz_0_xxy_0, \
                             g_0_yyyzzzz_0_xxz_0, \
                             g_0_yyyzzzz_0_xyy_0, \
                             g_0_yyyzzzz_0_xyz_0, \
                             g_0_yyyzzzz_0_xzz_0, \
                             g_0_yyyzzzz_0_yyy_0, \
                             g_0_yyyzzzz_0_yyz_0, \
                             g_0_yyyzzzz_0_yzz_0, \
                             g_0_yyyzzzz_0_zzz_0, \
                             g_0_yyzzzz_0_xxx_0,  \
                             g_0_yyzzzz_0_xxx_1,  \
                             g_0_yyzzzz_0_xxz_0,  \
                             g_0_yyzzzz_0_xxz_1,  \
                             g_0_yyzzzz_0_xyz_0,  \
                             g_0_yyzzzz_0_xyz_1,  \
                             g_0_yyzzzz_0_xz_1,   \
                             g_0_yyzzzz_0_xzz_0,  \
                             g_0_yyzzzz_0_xzz_1,  \
                             g_0_yyzzzz_0_yyz_0,  \
                             g_0_yyzzzz_0_yyz_1,  \
                             g_0_yyzzzz_0_yz_1,   \
                             g_0_yyzzzz_0_yzz_0,  \
                             g_0_yyzzzz_0_yzz_1,  \
                             g_0_yyzzzz_0_zz_1,   \
                             g_0_yyzzzz_0_zzz_0,  \
                             g_0_yyzzzz_0_zzz_1,  \
                             g_0_yzzzz_0_xxx_0,   \
                             g_0_yzzzz_0_xxx_1,   \
                             g_0_yzzzz_0_xxz_0,   \
                             g_0_yzzzz_0_xxz_1,   \
                             g_0_yzzzz_0_xyz_0,   \
                             g_0_yzzzz_0_xyz_1,   \
                             g_0_yzzzz_0_xzz_0,   \
                             g_0_yzzzz_0_xzz_1,   \
                             g_0_yzzzz_0_yyz_0,   \
                             g_0_yzzzz_0_yyz_1,   \
                             g_0_yzzzz_0_yzz_0,   \
                             g_0_yzzzz_0_yzz_1,   \
                             g_0_yzzzz_0_zzz_0,   \
                             g_0_yzzzz_0_zzz_1,   \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzz_0_xxx_0[i] = 2.0 * g_0_yzzzz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxx_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxx_0[i] * pb_y +
                                 g_0_yyzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxy_0[i] = 3.0 * g_0_yyyzz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxy_0[i] * pb_z +
                                 g_0_yyyzzz_0_xxy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxz_0[i] = 2.0 * g_0_yzzzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxz_0[i] * pb_y +
                                 g_0_yyzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyy_0[i] = 3.0 * g_0_yyyzz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xyy_0[i] * pb_z +
                                 g_0_yyyzzz_0_xyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xyz_0[i] = 2.0 * g_0_yzzzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xz_1[i] * fi_abcd_0 +
                                 g_0_yyzzzz_0_xyz_0[i] * pb_y + g_0_yyzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xzz_0[i] = 2.0 * g_0_yzzzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xzz_0[i] * pb_y +
                                 g_0_yyzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyy_0[i] = 3.0 * g_0_yyyzz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_yyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_yyy_0[i] * pb_z +
                                 g_0_yyyzzz_0_yyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_yyz_0[i] = 2.0 * g_0_yzzzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyz_0[i] * pb_y + g_0_yyzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yzz_0[i] = 2.0 * g_0_yzzzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_yyzzzz_0_yzz_0[i] * pb_y + g_0_yyzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_zzz_0[i] = 2.0 * g_0_yzzzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_zzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zzz_0[i] * pb_y +
                                 g_0_yyzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 330-340 components of targeted buffer : SKSF

    auto g_0_yyzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 330);

    auto g_0_yyzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 331);

    auto g_0_yyzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 332);

    auto g_0_yyzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 333);

    auto g_0_yyzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 334);

    auto g_0_yyzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 335);

    auto g_0_yyzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 336);

    auto g_0_yyzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 337);

    auto g_0_yyzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 338);

    auto g_0_yyzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 339);

#pragma omp simd aligned(g_0_yyzzz_0_xxy_0,       \
                             g_0_yyzzz_0_xxy_1,   \
                             g_0_yyzzz_0_xyy_0,   \
                             g_0_yyzzz_0_xyy_1,   \
                             g_0_yyzzz_0_yyy_0,   \
                             g_0_yyzzz_0_yyy_1,   \
                             g_0_yyzzzz_0_xxy_0,  \
                             g_0_yyzzzz_0_xxy_1,  \
                             g_0_yyzzzz_0_xyy_0,  \
                             g_0_yyzzzz_0_xyy_1,  \
                             g_0_yyzzzz_0_yyy_0,  \
                             g_0_yyzzzz_0_yyy_1,  \
                             g_0_yyzzzzz_0_xxx_0, \
                             g_0_yyzzzzz_0_xxy_0, \
                             g_0_yyzzzzz_0_xxz_0, \
                             g_0_yyzzzzz_0_xyy_0, \
                             g_0_yyzzzzz_0_xyz_0, \
                             g_0_yyzzzzz_0_xzz_0, \
                             g_0_yyzzzzz_0_yyy_0, \
                             g_0_yyzzzzz_0_yyz_0, \
                             g_0_yyzzzzz_0_yzz_0, \
                             g_0_yyzzzzz_0_zzz_0, \
                             g_0_yzzzzz_0_xxx_0,  \
                             g_0_yzzzzz_0_xxx_1,  \
                             g_0_yzzzzz_0_xxz_0,  \
                             g_0_yzzzzz_0_xxz_1,  \
                             g_0_yzzzzz_0_xyz_0,  \
                             g_0_yzzzzz_0_xyz_1,  \
                             g_0_yzzzzz_0_xz_1,   \
                             g_0_yzzzzz_0_xzz_0,  \
                             g_0_yzzzzz_0_xzz_1,  \
                             g_0_yzzzzz_0_yyz_0,  \
                             g_0_yzzzzz_0_yyz_1,  \
                             g_0_yzzzzz_0_yz_1,   \
                             g_0_yzzzzz_0_yzz_0,  \
                             g_0_yzzzzz_0_yzz_1,  \
                             g_0_yzzzzz_0_zz_1,   \
                             g_0_yzzzzz_0_zzz_0,  \
                             g_0_yzzzzz_0_zzz_1,  \
                             g_0_zzzzz_0_xxx_0,   \
                             g_0_zzzzz_0_xxx_1,   \
                             g_0_zzzzz_0_xxz_0,   \
                             g_0_zzzzz_0_xxz_1,   \
                             g_0_zzzzz_0_xyz_0,   \
                             g_0_zzzzz_0_xyz_1,   \
                             g_0_zzzzz_0_xzz_0,   \
                             g_0_zzzzz_0_xzz_1,   \
                             g_0_zzzzz_0_yyz_0,   \
                             g_0_zzzzz_0_yyz_1,   \
                             g_0_zzzzz_0_yzz_0,   \
                             g_0_zzzzz_0_yzz_1,   \
                             g_0_zzzzz_0_zzz_0,   \
                             g_0_zzzzz_0_zzz_1,   \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzz_0_xxx_0[i] =
            g_0_zzzzz_0_xxx_0[i] * fi_ab_0 - g_0_zzzzz_0_xxx_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxx_0[i] * pb_y + g_0_yzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxy_0[i] = 4.0 * g_0_yyzzz_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxy_0[i] * pb_z +
                                 g_0_yyzzzz_0_xxy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxz_0[i] =
            g_0_zzzzz_0_xxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxz_0[i] * pb_y + g_0_yzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyy_0[i] = 4.0 * g_0_yyzzz_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xyy_0[i] * pb_z +
                                 g_0_yyzzzz_0_xyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xyz_0[i] = g_0_zzzzz_0_xyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xz_1[i] * fi_abcd_0 +
                                 g_0_yzzzzz_0_xyz_0[i] * pb_y + g_0_yzzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xzz_0[i] =
            g_0_zzzzz_0_xzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzz_0[i] * pb_y + g_0_yzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyy_0[i] = 4.0 * g_0_yyzzz_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_yyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_yyy_0[i] * pb_z +
                                 g_0_yyzzzz_0_yyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_yyz_0[i] = g_0_zzzzz_0_yyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_yz_1[i] * fi_abcd_0 +
                                 g_0_yzzzzz_0_yyz_0[i] * pb_y + g_0_yzzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yzz_0[i] = g_0_zzzzz_0_yzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zz_1[i] * fi_abcd_0 +
                                 g_0_yzzzzz_0_yzz_0[i] * pb_y + g_0_yzzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_zzz_0[i] =
            g_0_zzzzz_0_zzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzz_0[i] * pb_y + g_0_yzzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 340-350 components of targeted buffer : SKSF

    auto g_0_yzzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 340);

    auto g_0_yzzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 341);

    auto g_0_yzzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 342);

    auto g_0_yzzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 343);

    auto g_0_yzzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 344);

    auto g_0_yzzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 345);

    auto g_0_yzzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 346);

    auto g_0_yzzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 347);

    auto g_0_yzzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 348);

    auto g_0_yzzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 349);

#pragma omp simd aligned(g_0_yzzzzzz_0_xxx_0,     \
                             g_0_yzzzzzz_0_xxy_0, \
                             g_0_yzzzzzz_0_xxz_0, \
                             g_0_yzzzzzz_0_xyy_0, \
                             g_0_yzzzzzz_0_xyz_0, \
                             g_0_yzzzzzz_0_xzz_0, \
                             g_0_yzzzzzz_0_yyy_0, \
                             g_0_yzzzzzz_0_yyz_0, \
                             g_0_yzzzzzz_0_yzz_0, \
                             g_0_yzzzzzz_0_zzz_0, \
                             g_0_zzzzzz_0_xx_1,   \
                             g_0_zzzzzz_0_xxx_0,  \
                             g_0_zzzzzz_0_xxx_1,  \
                             g_0_zzzzzz_0_xxy_0,  \
                             g_0_zzzzzz_0_xxy_1,  \
                             g_0_zzzzzz_0_xxz_0,  \
                             g_0_zzzzzz_0_xxz_1,  \
                             g_0_zzzzzz_0_xy_1,   \
                             g_0_zzzzzz_0_xyy_0,  \
                             g_0_zzzzzz_0_xyy_1,  \
                             g_0_zzzzzz_0_xyz_0,  \
                             g_0_zzzzzz_0_xyz_1,  \
                             g_0_zzzzzz_0_xz_1,   \
                             g_0_zzzzzz_0_xzz_0,  \
                             g_0_zzzzzz_0_xzz_1,  \
                             g_0_zzzzzz_0_yy_1,   \
                             g_0_zzzzzz_0_yyy_0,  \
                             g_0_zzzzzz_0_yyy_1,  \
                             g_0_zzzzzz_0_yyz_0,  \
                             g_0_zzzzzz_0_yyz_1,  \
                             g_0_zzzzzz_0_yz_1,   \
                             g_0_zzzzzz_0_yzz_0,  \
                             g_0_zzzzzz_0_yzz_1,  \
                             g_0_zzzzzz_0_zz_1,   \
                             g_0_zzzzzz_0_zzz_0,  \
                             g_0_zzzzzz_0_zzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzz_0_xxx_0[i] = g_0_zzzzzz_0_xxx_0[i] * pb_y + g_0_zzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxy_0[i] = g_0_zzzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxy_0[i] * pb_y + g_0_zzzzzz_0_xxy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxz_0[i] = g_0_zzzzzz_0_xxz_0[i] * pb_y + g_0_zzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyy_0[i] = 2.0 * g_0_zzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyy_0[i] * pb_y + g_0_zzzzzz_0_xyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyz_0[i] = g_0_zzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyz_0[i] * pb_y + g_0_zzzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xzz_0[i] = g_0_zzzzzz_0_xzz_0[i] * pb_y + g_0_zzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyy_0[i] = 3.0 * g_0_zzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyy_0[i] * pb_y + g_0_zzzzzz_0_yyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyz_0[i] = 2.0 * g_0_zzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyz_0[i] * pb_y + g_0_zzzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yzz_0[i] = g_0_zzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzz_0[i] * pb_y + g_0_zzzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_zzz_0[i] = g_0_zzzzzz_0_zzz_0[i] * pb_y + g_0_zzzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 350-360 components of targeted buffer : SKSF

    auto g_0_zzzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sksf + 350);

    auto g_0_zzzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sksf + 351);

    auto g_0_zzzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sksf + 352);

    auto g_0_zzzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sksf + 353);

    auto g_0_zzzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sksf + 354);

    auto g_0_zzzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sksf + 355);

    auto g_0_zzzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sksf + 356);

    auto g_0_zzzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sksf + 357);

    auto g_0_zzzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sksf + 358);

    auto g_0_zzzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sksf + 359);

#pragma omp simd aligned(g_0_zzzzz_0_xxx_0,       \
                             g_0_zzzzz_0_xxx_1,   \
                             g_0_zzzzz_0_xxy_0,   \
                             g_0_zzzzz_0_xxy_1,   \
                             g_0_zzzzz_0_xxz_0,   \
                             g_0_zzzzz_0_xxz_1,   \
                             g_0_zzzzz_0_xyy_0,   \
                             g_0_zzzzz_0_xyy_1,   \
                             g_0_zzzzz_0_xyz_0,   \
                             g_0_zzzzz_0_xyz_1,   \
                             g_0_zzzzz_0_xzz_0,   \
                             g_0_zzzzz_0_xzz_1,   \
                             g_0_zzzzz_0_yyy_0,   \
                             g_0_zzzzz_0_yyy_1,   \
                             g_0_zzzzz_0_yyz_0,   \
                             g_0_zzzzz_0_yyz_1,   \
                             g_0_zzzzz_0_yzz_0,   \
                             g_0_zzzzz_0_yzz_1,   \
                             g_0_zzzzz_0_zzz_0,   \
                             g_0_zzzzz_0_zzz_1,   \
                             g_0_zzzzzz_0_xx_1,   \
                             g_0_zzzzzz_0_xxx_0,  \
                             g_0_zzzzzz_0_xxx_1,  \
                             g_0_zzzzzz_0_xxy_0,  \
                             g_0_zzzzzz_0_xxy_1,  \
                             g_0_zzzzzz_0_xxz_0,  \
                             g_0_zzzzzz_0_xxz_1,  \
                             g_0_zzzzzz_0_xy_1,   \
                             g_0_zzzzzz_0_xyy_0,  \
                             g_0_zzzzzz_0_xyy_1,  \
                             g_0_zzzzzz_0_xyz_0,  \
                             g_0_zzzzzz_0_xyz_1,  \
                             g_0_zzzzzz_0_xz_1,   \
                             g_0_zzzzzz_0_xzz_0,  \
                             g_0_zzzzzz_0_xzz_1,  \
                             g_0_zzzzzz_0_yy_1,   \
                             g_0_zzzzzz_0_yyy_0,  \
                             g_0_zzzzzz_0_yyy_1,  \
                             g_0_zzzzzz_0_yyz_0,  \
                             g_0_zzzzzz_0_yyz_1,  \
                             g_0_zzzzzz_0_yz_1,   \
                             g_0_zzzzzz_0_yzz_0,  \
                             g_0_zzzzzz_0_yzz_1,  \
                             g_0_zzzzzz_0_zz_1,   \
                             g_0_zzzzzz_0_zzz_0,  \
                             g_0_zzzzzz_0_zzz_1,  \
                             g_0_zzzzzzz_0_xxx_0, \
                             g_0_zzzzzzz_0_xxy_0, \
                             g_0_zzzzzzz_0_xxz_0, \
                             g_0_zzzzzzz_0_xyy_0, \
                             g_0_zzzzzzz_0_xyz_0, \
                             g_0_zzzzzzz_0_xzz_0, \
                             g_0_zzzzzzz_0_yyy_0, \
                             g_0_zzzzzzz_0_yyz_0, \
                             g_0_zzzzzzz_0_yzz_0, \
                             g_0_zzzzzzz_0_zzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzz_0_xxx_0[i] = 6.0 * g_0_zzzzz_0_xxx_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxx_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxx_0[i] * pb_z +
                                 g_0_zzzzzz_0_xxx_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxy_0[i] = 6.0 * g_0_zzzzz_0_xxy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxy_0[i] * pb_z +
                                 g_0_zzzzzz_0_xxy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxz_0[i] = 6.0 * g_0_zzzzz_0_xxz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xx_1[i] * fi_abcd_0 +
                                 g_0_zzzzzz_0_xxz_0[i] * pb_z + g_0_zzzzzz_0_xxz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyy_0[i] = 6.0 * g_0_zzzzz_0_xyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xyy_0[i] * pb_z +
                                 g_0_zzzzzz_0_xyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyz_0[i] = 6.0 * g_0_zzzzz_0_xyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xy_1[i] * fi_abcd_0 +
                                 g_0_zzzzzz_0_xyz_0[i] * pb_z + g_0_zzzzzz_0_xyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xzz_0[i] = 6.0 * g_0_zzzzz_0_xzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzz_0[i] * pb_z + g_0_zzzzzz_0_xzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyy_0[i] = 6.0 * g_0_zzzzz_0_yyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_yyy_0[i] * pb_z +
                                 g_0_zzzzzz_0_yyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyz_0[i] = 6.0 * g_0_zzzzz_0_yyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_yy_1[i] * fi_abcd_0 +
                                 g_0_zzzzzz_0_yyz_0[i] * pb_z + g_0_zzzzzz_0_yyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yzz_0[i] = 6.0 * g_0_zzzzz_0_yzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzz_0[i] * pb_z + g_0_zzzzzz_0_yzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_zzz_0[i] = 6.0 * g_0_zzzzz_0_zzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_zzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_zzz_0[i] * pb_z + g_0_zzzzzz_0_zzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
