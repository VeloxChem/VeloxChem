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

#include "ElectronRepulsionPrimRecSKSD.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sksd(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sksd,
                                  size_t                idx_eri_0_shsd,
                                  size_t                idx_eri_1_shsd,
                                  size_t                idx_eri_1_sisp,
                                  size_t                idx_eri_0_sisd,
                                  size_t                idx_eri_1_sisd,
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

    /// Set up components of auxilary buffer : SHSD

    auto g_0_xxxxx_0_xx_0 = pbuffer.data(idx_eri_0_shsd);

    auto g_0_xxxxx_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 1);

    auto g_0_xxxxx_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 2);

    auto g_0_xxxxx_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 3);

    auto g_0_xxxxx_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 4);

    auto g_0_xxxxx_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 5);

    auto g_0_xxxxy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 6);

    auto g_0_xxxxy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 8);

    auto g_0_xxxxz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 12);

    auto g_0_xxxxz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 13);

    auto g_0_xxxyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 18);

    auto g_0_xxxyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 19);

    auto g_0_xxxyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 20);

    auto g_0_xxxyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 21);

    auto g_0_xxxyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 22);

    auto g_0_xxxyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 23);

    auto g_0_xxxzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 30);

    auto g_0_xxxzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 31);

    auto g_0_xxxzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 32);

    auto g_0_xxxzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 33);

    auto g_0_xxxzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 34);

    auto g_0_xxxzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 35);

    auto g_0_xxyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 36);

    auto g_0_xxyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 37);

    auto g_0_xxyyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 38);

    auto g_0_xxyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 39);

    auto g_0_xxyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 40);

    auto g_0_xxyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 41);

    auto g_0_xxyyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 43);

    auto g_0_xxyzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 48);

    auto g_0_xxyzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 50);

    auto g_0_xxzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 54);

    auto g_0_xxzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 55);

    auto g_0_xxzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 56);

    auto g_0_xxzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 57);

    auto g_0_xxzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 58);

    auto g_0_xxzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 59);

    auto g_0_xyyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 61);

    auto g_0_xyyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 63);

    auto g_0_xyyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 64);

    auto g_0_xyyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 65);

    auto g_0_xyyzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 75);

    auto g_0_xyyzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 76);

    auto g_0_xyyzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 77);

    auto g_0_xzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 86);

    auto g_0_xzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 87);

    auto g_0_xzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 88);

    auto g_0_xzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 89);

    auto g_0_yyyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 90);

    auto g_0_yyyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 91);

    auto g_0_yyyyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 92);

    auto g_0_yyyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 93);

    auto g_0_yyyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 94);

    auto g_0_yyyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 95);

    auto g_0_yyyyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 97);

    auto g_0_yyyyz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 99);

    auto g_0_yyyzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 102);

    auto g_0_yyyzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 103);

    auto g_0_yyyzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 104);

    auto g_0_yyyzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 105);

    auto g_0_yyyzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 106);

    auto g_0_yyyzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 107);

    auto g_0_yyzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 108);

    auto g_0_yyzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 109);

    auto g_0_yyzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 110);

    auto g_0_yyzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 111);

    auto g_0_yyzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 112);

    auto g_0_yyzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 113);

    auto g_0_yzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 114);

    auto g_0_yzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 116);

    auto g_0_yzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 118);

    auto g_0_yzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 119);

    auto g_0_zzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 120);

    auto g_0_zzzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 121);

    auto g_0_zzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 122);

    auto g_0_zzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 123);

    auto g_0_zzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 124);

    auto g_0_zzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 125);

    /// Set up components of auxilary buffer : SHSD

    auto g_0_xxxxx_0_xx_1 = pbuffer.data(idx_eri_1_shsd);

    auto g_0_xxxxx_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 1);

    auto g_0_xxxxx_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 2);

    auto g_0_xxxxx_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 3);

    auto g_0_xxxxx_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 4);

    auto g_0_xxxxx_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 5);

    auto g_0_xxxxy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 6);

    auto g_0_xxxxy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 8);

    auto g_0_xxxxz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 12);

    auto g_0_xxxxz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 13);

    auto g_0_xxxyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 18);

    auto g_0_xxxyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 19);

    auto g_0_xxxyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 20);

    auto g_0_xxxyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 21);

    auto g_0_xxxyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 22);

    auto g_0_xxxyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 23);

    auto g_0_xxxzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 30);

    auto g_0_xxxzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 31);

    auto g_0_xxxzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 32);

    auto g_0_xxxzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 33);

    auto g_0_xxxzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 34);

    auto g_0_xxxzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 35);

    auto g_0_xxyyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 36);

    auto g_0_xxyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 37);

    auto g_0_xxyyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 38);

    auto g_0_xxyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 39);

    auto g_0_xxyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 40);

    auto g_0_xxyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 41);

    auto g_0_xxyyz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 43);

    auto g_0_xxyzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 48);

    auto g_0_xxyzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 50);

    auto g_0_xxzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 54);

    auto g_0_xxzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 55);

    auto g_0_xxzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 56);

    auto g_0_xxzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 57);

    auto g_0_xxzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 58);

    auto g_0_xxzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 59);

    auto g_0_xyyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 61);

    auto g_0_xyyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 63);

    auto g_0_xyyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 64);

    auto g_0_xyyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 65);

    auto g_0_xyyzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 75);

    auto g_0_xyyzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 76);

    auto g_0_xyyzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 77);

    auto g_0_xzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 86);

    auto g_0_xzzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 87);

    auto g_0_xzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 88);

    auto g_0_xzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 89);

    auto g_0_yyyyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 90);

    auto g_0_yyyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 91);

    auto g_0_yyyyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 92);

    auto g_0_yyyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 93);

    auto g_0_yyyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 94);

    auto g_0_yyyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 95);

    auto g_0_yyyyz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 97);

    auto g_0_yyyyz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 99);

    auto g_0_yyyzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 102);

    auto g_0_yyyzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 103);

    auto g_0_yyyzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 104);

    auto g_0_yyyzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 105);

    auto g_0_yyyzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 106);

    auto g_0_yyyzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 107);

    auto g_0_yyzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 108);

    auto g_0_yyzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 109);

    auto g_0_yyzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 110);

    auto g_0_yyzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 111);

    auto g_0_yyzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 112);

    auto g_0_yyzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 113);

    auto g_0_yzzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 114);

    auto g_0_yzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 116);

    auto g_0_yzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 118);

    auto g_0_yzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 119);

    auto g_0_zzzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 120);

    auto g_0_zzzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 121);

    auto g_0_zzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 122);

    auto g_0_zzzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 123);

    auto g_0_zzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 124);

    auto g_0_zzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 125);

    /// Set up components of auxilary buffer : SISP

    auto g_0_xxxxxx_0_x_1 = pbuffer.data(idx_eri_1_sisp);

    auto g_0_xxxxxx_0_y_1 = pbuffer.data(idx_eri_1_sisp + 1);

    auto g_0_xxxxxx_0_z_1 = pbuffer.data(idx_eri_1_sisp + 2);

    auto g_0_xxxxxz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 8);

    auto g_0_xxxxyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 9);

    auto g_0_xxxxyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 10);

    auto g_0_xxxxyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 11);

    auto g_0_xxxxzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 15);

    auto g_0_xxxxzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 16);

    auto g_0_xxxxzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 17);

    auto g_0_xxxyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 18);

    auto g_0_xxxyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 19);

    auto g_0_xxxyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 20);

    auto g_0_xxxzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 27);

    auto g_0_xxxzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 28);

    auto g_0_xxxzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 29);

    auto g_0_xxyyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 30);

    auto g_0_xxyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 31);

    auto g_0_xxyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 32);

    auto g_0_xxzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 42);

    auto g_0_xxzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 43);

    auto g_0_xxzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 44);

    auto g_0_xyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 46);

    auto g_0_xzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 62);

    auto g_0_yyyyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 63);

    auto g_0_yyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 64);

    auto g_0_yyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 65);

    auto g_0_yyyyyz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 68);

    auto g_0_yyyyzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 69);

    auto g_0_yyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 70);

    auto g_0_yyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 71);

    auto g_0_yyyzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 72);

    auto g_0_yyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 73);

    auto g_0_yyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 74);

    auto g_0_yyzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 75);

    auto g_0_yyzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 76);

    auto g_0_yyzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 77);

    auto g_0_yzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 79);

    auto g_0_yzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 80);

    auto g_0_zzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 81);

    auto g_0_zzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 82);

    auto g_0_zzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 83);

    /// Set up components of auxilary buffer : SISD

    auto g_0_xxxxxx_0_xx_0 = pbuffer.data(idx_eri_0_sisd);

    auto g_0_xxxxxx_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 1);

    auto g_0_xxxxxx_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 2);

    auto g_0_xxxxxx_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 3);

    auto g_0_xxxxxx_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 4);

    auto g_0_xxxxxx_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 5);

    auto g_0_xxxxxy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 6);

    auto g_0_xxxxxy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 7);

    auto g_0_xxxxxy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 8);

    auto g_0_xxxxxy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 9);

    auto g_0_xxxxxz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 12);

    auto g_0_xxxxxz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 13);

    auto g_0_xxxxxz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 14);

    auto g_0_xxxxxz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 16);

    auto g_0_xxxxxz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 17);

    auto g_0_xxxxyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 18);

    auto g_0_xxxxyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 19);

    auto g_0_xxxxyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 20);

    auto g_0_xxxxyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 21);

    auto g_0_xxxxyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 22);

    auto g_0_xxxxyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 23);

    auto g_0_xxxxzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 30);

    auto g_0_xxxxzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 31);

    auto g_0_xxxxzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 32);

    auto g_0_xxxxzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 33);

    auto g_0_xxxxzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 34);

    auto g_0_xxxxzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 35);

    auto g_0_xxxyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 36);

    auto g_0_xxxyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 37);

    auto g_0_xxxyyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 38);

    auto g_0_xxxyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 39);

    auto g_0_xxxyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 40);

    auto g_0_xxxyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 41);

    auto g_0_xxxyyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 43);

    auto g_0_xxxyzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 48);

    auto g_0_xxxyzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 50);

    auto g_0_xxxzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 54);

    auto g_0_xxxzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 55);

    auto g_0_xxxzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 56);

    auto g_0_xxxzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 57);

    auto g_0_xxxzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 58);

    auto g_0_xxxzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 59);

    auto g_0_xxyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 60);

    auto g_0_xxyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 61);

    auto g_0_xxyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 62);

    auto g_0_xxyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 63);

    auto g_0_xxyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 64);

    auto g_0_xxyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 65);

    auto g_0_xxyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 67);

    auto g_0_xxyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 72);

    auto g_0_xxyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 73);

    auto g_0_xxyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 74);

    auto g_0_xxyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 75);

    auto g_0_xxyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 76);

    auto g_0_xxyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 77);

    auto g_0_xxyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 78);

    auto g_0_xxyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 80);

    auto g_0_xxzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 84);

    auto g_0_xxzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 85);

    auto g_0_xxzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 86);

    auto g_0_xxzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 87);

    auto g_0_xxzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 88);

    auto g_0_xxzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 89);

    auto g_0_xyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 90);

    auto g_0_xyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 91);

    auto g_0_xyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 93);

    auto g_0_xyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 94);

    auto g_0_xyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 95);

    auto g_0_xyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 105);

    auto g_0_xyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 106);

    auto g_0_xyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 107);

    auto g_0_xyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 111);

    auto g_0_xyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 112);

    auto g_0_xyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 113);

    auto g_0_xzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 120);

    auto g_0_xzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 122);

    auto g_0_xzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 123);

    auto g_0_xzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 124);

    auto g_0_xzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 125);

    auto g_0_yyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 126);

    auto g_0_yyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 127);

    auto g_0_yyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 128);

    auto g_0_yyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 129);

    auto g_0_yyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 130);

    auto g_0_yyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 131);

    auto g_0_yyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 133);

    auto g_0_yyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 134);

    auto g_0_yyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 135);

    auto g_0_yyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 136);

    auto g_0_yyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 137);

    auto g_0_yyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 138);

    auto g_0_yyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 139);

    auto g_0_yyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 140);

    auto g_0_yyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 141);

    auto g_0_yyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 142);

    auto g_0_yyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 143);

    auto g_0_yyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 144);

    auto g_0_yyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 145);

    auto g_0_yyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 146);

    auto g_0_yyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 147);

    auto g_0_yyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 148);

    auto g_0_yyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 149);

    auto g_0_yyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 150);

    auto g_0_yyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 151);

    auto g_0_yyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 152);

    auto g_0_yyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 153);

    auto g_0_yyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 154);

    auto g_0_yyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 155);

    auto g_0_yzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 156);

    auto g_0_yzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 157);

    auto g_0_yzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 158);

    auto g_0_yzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 159);

    auto g_0_yzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 160);

    auto g_0_yzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 161);

    auto g_0_zzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 162);

    auto g_0_zzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 163);

    auto g_0_zzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 164);

    auto g_0_zzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 165);

    auto g_0_zzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 166);

    auto g_0_zzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 167);

    /// Set up components of auxilary buffer : SISD

    auto g_0_xxxxxx_0_xx_1 = pbuffer.data(idx_eri_1_sisd);

    auto g_0_xxxxxx_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 1);

    auto g_0_xxxxxx_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 2);

    auto g_0_xxxxxx_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 3);

    auto g_0_xxxxxx_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 4);

    auto g_0_xxxxxx_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 5);

    auto g_0_xxxxxy_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 6);

    auto g_0_xxxxxy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 7);

    auto g_0_xxxxxy_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 8);

    auto g_0_xxxxxy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 9);

    auto g_0_xxxxxz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 12);

    auto g_0_xxxxxz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 13);

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

    auto g_0_xxxyyz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 43);

    auto g_0_xxxyzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 48);

    auto g_0_xxxyzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 50);

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

    auto g_0_xxyyyz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 67);

    auto g_0_xxyyzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 72);

    auto g_0_xxyyzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 73);

    auto g_0_xxyyzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 74);

    auto g_0_xxyyzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 75);

    auto g_0_xxyyzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 76);

    auto g_0_xxyyzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 77);

    auto g_0_xxyzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 78);

    auto g_0_xxyzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 80);

    auto g_0_xxzzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 84);

    auto g_0_xxzzzz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 85);

    auto g_0_xxzzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 86);

    auto g_0_xxzzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 87);

    auto g_0_xxzzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 88);

    auto g_0_xxzzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 89);

    auto g_0_xyyyyy_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 90);

    auto g_0_xyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 91);

    auto g_0_xyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 93);

    auto g_0_xyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 94);

    auto g_0_xyyyyy_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 95);

    auto g_0_xyyyzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 105);

    auto g_0_xyyyzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 106);

    auto g_0_xyyyzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 107);

    auto g_0_xyyzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 111);

    auto g_0_xyyzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 112);

    auto g_0_xyyzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 113);

    auto g_0_xzzzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 120);

    auto g_0_xzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 122);

    auto g_0_xzzzzz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 123);

    auto g_0_xzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 124);

    auto g_0_xzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 125);

    auto g_0_yyyyyy_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 126);

    auto g_0_yyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 127);

    auto g_0_yyyyyy_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 128);

    auto g_0_yyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 129);

    auto g_0_yyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_sisd + 130);

    auto g_0_yyyyyy_0_zz_1 = pbuffer.data(idx_eri_1_sisd + 131);

    auto g_0_yyyyyz_0_xy_1 = pbuffer.data(idx_eri_1_sisd + 133);

    auto g_0_yyyyyz_0_xz_1 = pbuffer.data(idx_eri_1_sisd + 134);

    auto g_0_yyyyyz_0_yy_1 = pbuffer.data(idx_eri_1_sisd + 135);

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

    auto g_0_yzzzzz_0_xx_1 = pbuffer.data(idx_eri_1_sisd + 156);

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

    /// Set up 0-6 components of targeted buffer : SKSD

    auto g_0_xxxxxxx_0_xx_0 = pbuffer.data(idx_eri_0_sksd);

    auto g_0_xxxxxxx_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 1);

    auto g_0_xxxxxxx_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 2);

    auto g_0_xxxxxxx_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 3);

    auto g_0_xxxxxxx_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 4);

    auto g_0_xxxxxxx_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 5);

#pragma omp simd aligned(g_0_xxxxx_0_xx_0,       \
                             g_0_xxxxx_0_xx_1,   \
                             g_0_xxxxx_0_xy_0,   \
                             g_0_xxxxx_0_xy_1,   \
                             g_0_xxxxx_0_xz_0,   \
                             g_0_xxxxx_0_xz_1,   \
                             g_0_xxxxx_0_yy_0,   \
                             g_0_xxxxx_0_yy_1,   \
                             g_0_xxxxx_0_yz_0,   \
                             g_0_xxxxx_0_yz_1,   \
                             g_0_xxxxx_0_zz_0,   \
                             g_0_xxxxx_0_zz_1,   \
                             g_0_xxxxxx_0_x_1,   \
                             g_0_xxxxxx_0_xx_0,  \
                             g_0_xxxxxx_0_xx_1,  \
                             g_0_xxxxxx_0_xy_0,  \
                             g_0_xxxxxx_0_xy_1,  \
                             g_0_xxxxxx_0_xz_0,  \
                             g_0_xxxxxx_0_xz_1,  \
                             g_0_xxxxxx_0_y_1,   \
                             g_0_xxxxxx_0_yy_0,  \
                             g_0_xxxxxx_0_yy_1,  \
                             g_0_xxxxxx_0_yz_0,  \
                             g_0_xxxxxx_0_yz_1,  \
                             g_0_xxxxxx_0_z_1,   \
                             g_0_xxxxxx_0_zz_0,  \
                             g_0_xxxxxx_0_zz_1,  \
                             g_0_xxxxxxx_0_xx_0, \
                             g_0_xxxxxxx_0_xy_0, \
                             g_0_xxxxxxx_0_xz_0, \
                             g_0_xxxxxxx_0_yy_0, \
                             g_0_xxxxxxx_0_yz_0, \
                             g_0_xxxxxxx_0_zz_0, \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_xx_0[i] = 6.0 * g_0_xxxxx_0_xx_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xx_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_x_1[i] * fi_abcd_0 +
                                g_0_xxxxxx_0_xx_0[i] * pb_x + g_0_xxxxxx_0_xx_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xy_0[i] = 6.0 * g_0_xxxxx_0_xy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xy_1[i] * fti_ab_0 + g_0_xxxxxx_0_y_1[i] * fi_abcd_0 +
                                g_0_xxxxxx_0_xy_0[i] * pb_x + g_0_xxxxxx_0_xy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xz_0[i] = 6.0 * g_0_xxxxx_0_xz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xz_1[i] * fti_ab_0 + g_0_xxxxxx_0_z_1[i] * fi_abcd_0 +
                                g_0_xxxxxx_0_xz_0[i] * pb_x + g_0_xxxxxx_0_xz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yy_0[i] =
            6.0 * g_0_xxxxx_0_yy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yy_0[i] * pb_x + g_0_xxxxxx_0_yy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yz_0[i] =
            6.0 * g_0_xxxxx_0_yz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yz_0[i] * pb_x + g_0_xxxxxx_0_yz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_zz_0[i] =
            6.0 * g_0_xxxxx_0_zz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_zz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zz_0[i] * pb_x + g_0_xxxxxx_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : SKSD

    auto g_0_xxxxxxy_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 6);

    auto g_0_xxxxxxy_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 7);

    auto g_0_xxxxxxy_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 8);

    auto g_0_xxxxxxy_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 9);

    auto g_0_xxxxxxy_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 10);

    auto g_0_xxxxxxy_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 11);

#pragma omp simd aligned(g_0_xxxxxx_0_x_1,       \
                             g_0_xxxxxx_0_xx_0,  \
                             g_0_xxxxxx_0_xx_1,  \
                             g_0_xxxxxx_0_xy_0,  \
                             g_0_xxxxxx_0_xy_1,  \
                             g_0_xxxxxx_0_xz_0,  \
                             g_0_xxxxxx_0_xz_1,  \
                             g_0_xxxxxx_0_y_1,   \
                             g_0_xxxxxx_0_yy_0,  \
                             g_0_xxxxxx_0_yy_1,  \
                             g_0_xxxxxx_0_yz_0,  \
                             g_0_xxxxxx_0_yz_1,  \
                             g_0_xxxxxx_0_z_1,   \
                             g_0_xxxxxx_0_zz_0,  \
                             g_0_xxxxxx_0_zz_1,  \
                             g_0_xxxxxxy_0_xx_0, \
                             g_0_xxxxxxy_0_xy_0, \
                             g_0_xxxxxxy_0_xz_0, \
                             g_0_xxxxxxy_0_yy_0, \
                             g_0_xxxxxxy_0_yz_0, \
                             g_0_xxxxxxy_0_zz_0, \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxy_0_xx_0[i] = g_0_xxxxxx_0_xx_0[i] * pb_y + g_0_xxxxxx_0_xx_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xy_0[i] = g_0_xxxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xy_0[i] * pb_y + g_0_xxxxxx_0_xy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xz_0[i] = g_0_xxxxxx_0_xz_0[i] * pb_y + g_0_xxxxxx_0_xz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yy_0[i] = 2.0 * g_0_xxxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yy_0[i] * pb_y + g_0_xxxxxx_0_yy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yz_0[i] = g_0_xxxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yz_0[i] * pb_y + g_0_xxxxxx_0_yz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_zz_0[i] = g_0_xxxxxx_0_zz_0[i] * pb_y + g_0_xxxxxx_0_zz_1[i] * wp_y[i];
    }

    /// Set up 12-18 components of targeted buffer : SKSD

    auto g_0_xxxxxxz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 12);

    auto g_0_xxxxxxz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 13);

    auto g_0_xxxxxxz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 14);

    auto g_0_xxxxxxz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 15);

    auto g_0_xxxxxxz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 16);

    auto g_0_xxxxxxz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 17);

#pragma omp simd aligned(g_0_xxxxxx_0_x_1,       \
                             g_0_xxxxxx_0_xx_0,  \
                             g_0_xxxxxx_0_xx_1,  \
                             g_0_xxxxxx_0_xy_0,  \
                             g_0_xxxxxx_0_xy_1,  \
                             g_0_xxxxxx_0_xz_0,  \
                             g_0_xxxxxx_0_xz_1,  \
                             g_0_xxxxxx_0_y_1,   \
                             g_0_xxxxxx_0_yy_0,  \
                             g_0_xxxxxx_0_yy_1,  \
                             g_0_xxxxxx_0_yz_0,  \
                             g_0_xxxxxx_0_yz_1,  \
                             g_0_xxxxxx_0_z_1,   \
                             g_0_xxxxxx_0_zz_0,  \
                             g_0_xxxxxx_0_zz_1,  \
                             g_0_xxxxxxz_0_xx_0, \
                             g_0_xxxxxxz_0_xy_0, \
                             g_0_xxxxxxz_0_xz_0, \
                             g_0_xxxxxxz_0_yy_0, \
                             g_0_xxxxxxz_0_yz_0, \
                             g_0_xxxxxxz_0_zz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxz_0_xx_0[i] = g_0_xxxxxx_0_xx_0[i] * pb_z + g_0_xxxxxx_0_xx_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xy_0[i] = g_0_xxxxxx_0_xy_0[i] * pb_z + g_0_xxxxxx_0_xy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xz_0[i] = g_0_xxxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xz_0[i] * pb_z + g_0_xxxxxx_0_xz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yy_0[i] = g_0_xxxxxx_0_yy_0[i] * pb_z + g_0_xxxxxx_0_yy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yz_0[i] = g_0_xxxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yz_0[i] * pb_z + g_0_xxxxxx_0_yz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_zz_0[i] = 2.0 * g_0_xxxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxxx_0_zz_0[i] * pb_z + g_0_xxxxxx_0_zz_1[i] * wp_z[i];
    }

    /// Set up 18-24 components of targeted buffer : SKSD

    auto g_0_xxxxxyy_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 18);

    auto g_0_xxxxxyy_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 19);

    auto g_0_xxxxxyy_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 20);

    auto g_0_xxxxxyy_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 21);

    auto g_0_xxxxxyy_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 22);

    auto g_0_xxxxxyy_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 23);

#pragma omp simd aligned(g_0_xxxxx_0_xx_0,       \
                             g_0_xxxxx_0_xx_1,   \
                             g_0_xxxxx_0_xz_0,   \
                             g_0_xxxxx_0_xz_1,   \
                             g_0_xxxxxy_0_xx_0,  \
                             g_0_xxxxxy_0_xx_1,  \
                             g_0_xxxxxy_0_xz_0,  \
                             g_0_xxxxxy_0_xz_1,  \
                             g_0_xxxxxyy_0_xx_0, \
                             g_0_xxxxxyy_0_xy_0, \
                             g_0_xxxxxyy_0_xz_0, \
                             g_0_xxxxxyy_0_yy_0, \
                             g_0_xxxxxyy_0_yz_0, \
                             g_0_xxxxxyy_0_zz_0, \
                             g_0_xxxxyy_0_xy_0,  \
                             g_0_xxxxyy_0_xy_1,  \
                             g_0_xxxxyy_0_y_1,   \
                             g_0_xxxxyy_0_yy_0,  \
                             g_0_xxxxyy_0_yy_1,  \
                             g_0_xxxxyy_0_yz_0,  \
                             g_0_xxxxyy_0_yz_1,  \
                             g_0_xxxxyy_0_zz_0,  \
                             g_0_xxxxyy_0_zz_1,  \
                             g_0_xxxyy_0_xy_0,   \
                             g_0_xxxyy_0_xy_1,   \
                             g_0_xxxyy_0_yy_0,   \
                             g_0_xxxyy_0_yy_1,   \
                             g_0_xxxyy_0_yz_0,   \
                             g_0_xxxyy_0_yz_1,   \
                             g_0_xxxyy_0_zz_0,   \
                             g_0_xxxyy_0_zz_1,   \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyy_0_xx_0[i] =
            g_0_xxxxx_0_xx_0[i] * fi_ab_0 - g_0_xxxxx_0_xx_1[i] * fti_ab_0 + g_0_xxxxxy_0_xx_0[i] * pb_y + g_0_xxxxxy_0_xx_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xy_0[i] = 4.0 * g_0_xxxyy_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xy_1[i] * fti_ab_0 + g_0_xxxxyy_0_y_1[i] * fi_abcd_0 +
                                g_0_xxxxyy_0_xy_0[i] * pb_x + g_0_xxxxyy_0_xy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xz_0[i] =
            g_0_xxxxx_0_xz_0[i] * fi_ab_0 - g_0_xxxxx_0_xz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xz_0[i] * pb_y + g_0_xxxxxy_0_xz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_yy_0[i] =
            4.0 * g_0_xxxyy_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yy_0[i] * pb_x + g_0_xxxxyy_0_yy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yz_0[i] =
            4.0 * g_0_xxxyy_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yz_0[i] * pb_x + g_0_xxxxyy_0_yz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_zz_0[i] =
            4.0 * g_0_xxxyy_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_zz_1[i] * fti_ab_0 + g_0_xxxxyy_0_zz_0[i] * pb_x + g_0_xxxxyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 24-30 components of targeted buffer : SKSD

    auto g_0_xxxxxyz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 24);

    auto g_0_xxxxxyz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 25);

    auto g_0_xxxxxyz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 26);

    auto g_0_xxxxxyz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 27);

    auto g_0_xxxxxyz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 28);

    auto g_0_xxxxxyz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 29);

#pragma omp simd aligned(g_0_xxxxxy_0_xy_0,      \
                             g_0_xxxxxy_0_xy_1,  \
                             g_0_xxxxxy_0_yy_0,  \
                             g_0_xxxxxy_0_yy_1,  \
                             g_0_xxxxxyz_0_xx_0, \
                             g_0_xxxxxyz_0_xy_0, \
                             g_0_xxxxxyz_0_xz_0, \
                             g_0_xxxxxyz_0_yy_0, \
                             g_0_xxxxxyz_0_yz_0, \
                             g_0_xxxxxyz_0_zz_0, \
                             g_0_xxxxxz_0_xx_0,  \
                             g_0_xxxxxz_0_xx_1,  \
                             g_0_xxxxxz_0_xz_0,  \
                             g_0_xxxxxz_0_xz_1,  \
                             g_0_xxxxxz_0_yz_0,  \
                             g_0_xxxxxz_0_yz_1,  \
                             g_0_xxxxxz_0_z_1,   \
                             g_0_xxxxxz_0_zz_0,  \
                             g_0_xxxxxz_0_zz_1,  \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyz_0_xx_0[i] = g_0_xxxxxz_0_xx_0[i] * pb_y + g_0_xxxxxz_0_xx_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xy_0[i] = g_0_xxxxxy_0_xy_0[i] * pb_z + g_0_xxxxxy_0_xy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xz_0[i] = g_0_xxxxxz_0_xz_0[i] * pb_y + g_0_xxxxxz_0_xz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yy_0[i] = g_0_xxxxxy_0_yy_0[i] * pb_z + g_0_xxxxxy_0_yy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_yz_0[i] = g_0_xxxxxz_0_z_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yz_0[i] * pb_y + g_0_xxxxxz_0_yz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_zz_0[i] = g_0_xxxxxz_0_zz_0[i] * pb_y + g_0_xxxxxz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 30-36 components of targeted buffer : SKSD

    auto g_0_xxxxxzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 30);

    auto g_0_xxxxxzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 31);

    auto g_0_xxxxxzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 32);

    auto g_0_xxxxxzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 33);

    auto g_0_xxxxxzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 34);

    auto g_0_xxxxxzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 35);

#pragma omp simd aligned(g_0_xxxxx_0_xx_0,       \
                             g_0_xxxxx_0_xx_1,   \
                             g_0_xxxxx_0_xy_0,   \
                             g_0_xxxxx_0_xy_1,   \
                             g_0_xxxxxz_0_xx_0,  \
                             g_0_xxxxxz_0_xx_1,  \
                             g_0_xxxxxz_0_xy_0,  \
                             g_0_xxxxxz_0_xy_1,  \
                             g_0_xxxxxzz_0_xx_0, \
                             g_0_xxxxxzz_0_xy_0, \
                             g_0_xxxxxzz_0_xz_0, \
                             g_0_xxxxxzz_0_yy_0, \
                             g_0_xxxxxzz_0_yz_0, \
                             g_0_xxxxxzz_0_zz_0, \
                             g_0_xxxxzz_0_xz_0,  \
                             g_0_xxxxzz_0_xz_1,  \
                             g_0_xxxxzz_0_yy_0,  \
                             g_0_xxxxzz_0_yy_1,  \
                             g_0_xxxxzz_0_yz_0,  \
                             g_0_xxxxzz_0_yz_1,  \
                             g_0_xxxxzz_0_z_1,   \
                             g_0_xxxxzz_0_zz_0,  \
                             g_0_xxxxzz_0_zz_1,  \
                             g_0_xxxzz_0_xz_0,   \
                             g_0_xxxzz_0_xz_1,   \
                             g_0_xxxzz_0_yy_0,   \
                             g_0_xxxzz_0_yy_1,   \
                             g_0_xxxzz_0_yz_0,   \
                             g_0_xxxzz_0_yz_1,   \
                             g_0_xxxzz_0_zz_0,   \
                             g_0_xxxzz_0_zz_1,   \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzz_0_xx_0[i] =
            g_0_xxxxx_0_xx_0[i] * fi_ab_0 - g_0_xxxxx_0_xx_1[i] * fti_ab_0 + g_0_xxxxxz_0_xx_0[i] * pb_z + g_0_xxxxxz_0_xx_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xy_0[i] =
            g_0_xxxxx_0_xy_0[i] * fi_ab_0 - g_0_xxxxx_0_xy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xy_0[i] * pb_z + g_0_xxxxxz_0_xy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xz_0[i] = 4.0 * g_0_xxxzz_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xz_1[i] * fti_ab_0 + g_0_xxxxzz_0_z_1[i] * fi_abcd_0 +
                                g_0_xxxxzz_0_xz_0[i] * pb_x + g_0_xxxxzz_0_xz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yy_0[i] =
            4.0 * g_0_xxxzz_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yy_1[i] * fti_ab_0 + g_0_xxxxzz_0_yy_0[i] * pb_x + g_0_xxxxzz_0_yy_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yz_0[i] =
            4.0 * g_0_xxxzz_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yz_0[i] * pb_x + g_0_xxxxzz_0_yz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_zz_0[i] =
            4.0 * g_0_xxxzz_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_zz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zz_0[i] * pb_x + g_0_xxxxzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 36-42 components of targeted buffer : SKSD

    auto g_0_xxxxyyy_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 36);

    auto g_0_xxxxyyy_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 37);

    auto g_0_xxxxyyy_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 38);

    auto g_0_xxxxyyy_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 39);

    auto g_0_xxxxyyy_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 40);

    auto g_0_xxxxyyy_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 41);

#pragma omp simd aligned(g_0_xxxxy_0_xx_0,       \
                             g_0_xxxxy_0_xx_1,   \
                             g_0_xxxxy_0_xz_0,   \
                             g_0_xxxxy_0_xz_1,   \
                             g_0_xxxxyy_0_xx_0,  \
                             g_0_xxxxyy_0_xx_1,  \
                             g_0_xxxxyy_0_xz_0,  \
                             g_0_xxxxyy_0_xz_1,  \
                             g_0_xxxxyyy_0_xx_0, \
                             g_0_xxxxyyy_0_xy_0, \
                             g_0_xxxxyyy_0_xz_0, \
                             g_0_xxxxyyy_0_yy_0, \
                             g_0_xxxxyyy_0_yz_0, \
                             g_0_xxxxyyy_0_zz_0, \
                             g_0_xxxyyy_0_xy_0,  \
                             g_0_xxxyyy_0_xy_1,  \
                             g_0_xxxyyy_0_y_1,   \
                             g_0_xxxyyy_0_yy_0,  \
                             g_0_xxxyyy_0_yy_1,  \
                             g_0_xxxyyy_0_yz_0,  \
                             g_0_xxxyyy_0_yz_1,  \
                             g_0_xxxyyy_0_zz_0,  \
                             g_0_xxxyyy_0_zz_1,  \
                             g_0_xxyyy_0_xy_0,   \
                             g_0_xxyyy_0_xy_1,   \
                             g_0_xxyyy_0_yy_0,   \
                             g_0_xxyyy_0_yy_1,   \
                             g_0_xxyyy_0_yz_0,   \
                             g_0_xxyyy_0_yz_1,   \
                             g_0_xxyyy_0_zz_0,   \
                             g_0_xxyyy_0_zz_1,   \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyy_0_xx_0[i] =
            2.0 * g_0_xxxxy_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xx_1[i] * fti_ab_0 + g_0_xxxxyy_0_xx_0[i] * pb_y + g_0_xxxxyy_0_xx_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xy_0[i] = 3.0 * g_0_xxyyy_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xy_1[i] * fti_ab_0 + g_0_xxxyyy_0_y_1[i] * fi_abcd_0 +
                                g_0_xxxyyy_0_xy_0[i] * pb_x + g_0_xxxyyy_0_xy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xz_0[i] =
            2.0 * g_0_xxxxy_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xz_0[i] * pb_y + g_0_xxxxyy_0_xz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_yy_0[i] =
            3.0 * g_0_xxyyy_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yy_0[i] * pb_x + g_0_xxxyyy_0_yy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yz_0[i] =
            3.0 * g_0_xxyyy_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yz_0[i] * pb_x + g_0_xxxyyy_0_yz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_zz_0[i] =
            3.0 * g_0_xxyyy_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_zz_1[i] * fti_ab_0 + g_0_xxxyyy_0_zz_0[i] * pb_x + g_0_xxxyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 42-48 components of targeted buffer : SKSD

    auto g_0_xxxxyyz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 42);

    auto g_0_xxxxyyz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 43);

    auto g_0_xxxxyyz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 44);

    auto g_0_xxxxyyz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 45);

    auto g_0_xxxxyyz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 46);

    auto g_0_xxxxyyz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 47);

#pragma omp simd aligned(g_0_xxxxyy_0_x_1,       \
                             g_0_xxxxyy_0_xx_0,  \
                             g_0_xxxxyy_0_xx_1,  \
                             g_0_xxxxyy_0_xy_0,  \
                             g_0_xxxxyy_0_xy_1,  \
                             g_0_xxxxyy_0_xz_0,  \
                             g_0_xxxxyy_0_xz_1,  \
                             g_0_xxxxyy_0_y_1,   \
                             g_0_xxxxyy_0_yy_0,  \
                             g_0_xxxxyy_0_yy_1,  \
                             g_0_xxxxyy_0_yz_0,  \
                             g_0_xxxxyy_0_yz_1,  \
                             g_0_xxxxyy_0_z_1,   \
                             g_0_xxxxyy_0_zz_0,  \
                             g_0_xxxxyy_0_zz_1,  \
                             g_0_xxxxyyz_0_xx_0, \
                             g_0_xxxxyyz_0_xy_0, \
                             g_0_xxxxyyz_0_xz_0, \
                             g_0_xxxxyyz_0_yy_0, \
                             g_0_xxxxyyz_0_yz_0, \
                             g_0_xxxxyyz_0_zz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyz_0_xx_0[i] = g_0_xxxxyy_0_xx_0[i] * pb_z + g_0_xxxxyy_0_xx_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xy_0[i] = g_0_xxxxyy_0_xy_0[i] * pb_z + g_0_xxxxyy_0_xy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xz_0[i] = g_0_xxxxyy_0_x_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xz_0[i] * pb_z + g_0_xxxxyy_0_xz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yy_0[i] = g_0_xxxxyy_0_yy_0[i] * pb_z + g_0_xxxxyy_0_yy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yz_0[i] = g_0_xxxxyy_0_y_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yz_0[i] * pb_z + g_0_xxxxyy_0_yz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_zz_0[i] = 2.0 * g_0_xxxxyy_0_z_1[i] * fi_abcd_0 + g_0_xxxxyy_0_zz_0[i] * pb_z + g_0_xxxxyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 48-54 components of targeted buffer : SKSD

    auto g_0_xxxxyzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 48);

    auto g_0_xxxxyzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 49);

    auto g_0_xxxxyzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 50);

    auto g_0_xxxxyzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 51);

    auto g_0_xxxxyzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 52);

    auto g_0_xxxxyzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 53);

#pragma omp simd aligned(g_0_xxxxyzz_0_xx_0,     \
                             g_0_xxxxyzz_0_xy_0, \
                             g_0_xxxxyzz_0_xz_0, \
                             g_0_xxxxyzz_0_yy_0, \
                             g_0_xxxxyzz_0_yz_0, \
                             g_0_xxxxyzz_0_zz_0, \
                             g_0_xxxxzz_0_x_1,   \
                             g_0_xxxxzz_0_xx_0,  \
                             g_0_xxxxzz_0_xx_1,  \
                             g_0_xxxxzz_0_xy_0,  \
                             g_0_xxxxzz_0_xy_1,  \
                             g_0_xxxxzz_0_xz_0,  \
                             g_0_xxxxzz_0_xz_1,  \
                             g_0_xxxxzz_0_y_1,   \
                             g_0_xxxxzz_0_yy_0,  \
                             g_0_xxxxzz_0_yy_1,  \
                             g_0_xxxxzz_0_yz_0,  \
                             g_0_xxxxzz_0_yz_1,  \
                             g_0_xxxxzz_0_z_1,   \
                             g_0_xxxxzz_0_zz_0,  \
                             g_0_xxxxzz_0_zz_1,  \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzz_0_xx_0[i] = g_0_xxxxzz_0_xx_0[i] * pb_y + g_0_xxxxzz_0_xx_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xy_0[i] = g_0_xxxxzz_0_x_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xy_0[i] * pb_y + g_0_xxxxzz_0_xy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xz_0[i] = g_0_xxxxzz_0_xz_0[i] * pb_y + g_0_xxxxzz_0_xz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yy_0[i] = 2.0 * g_0_xxxxzz_0_y_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yy_0[i] * pb_y + g_0_xxxxzz_0_yy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yz_0[i] = g_0_xxxxzz_0_z_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yz_0[i] * pb_y + g_0_xxxxzz_0_yz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_zz_0[i] = g_0_xxxxzz_0_zz_0[i] * pb_y + g_0_xxxxzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 54-60 components of targeted buffer : SKSD

    auto g_0_xxxxzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 54);

    auto g_0_xxxxzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 55);

    auto g_0_xxxxzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 56);

    auto g_0_xxxxzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 57);

    auto g_0_xxxxzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 58);

    auto g_0_xxxxzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 59);

#pragma omp simd aligned(g_0_xxxxz_0_xx_0,       \
                             g_0_xxxxz_0_xx_1,   \
                             g_0_xxxxz_0_xy_0,   \
                             g_0_xxxxz_0_xy_1,   \
                             g_0_xxxxzz_0_xx_0,  \
                             g_0_xxxxzz_0_xx_1,  \
                             g_0_xxxxzz_0_xy_0,  \
                             g_0_xxxxzz_0_xy_1,  \
                             g_0_xxxxzzz_0_xx_0, \
                             g_0_xxxxzzz_0_xy_0, \
                             g_0_xxxxzzz_0_xz_0, \
                             g_0_xxxxzzz_0_yy_0, \
                             g_0_xxxxzzz_0_yz_0, \
                             g_0_xxxxzzz_0_zz_0, \
                             g_0_xxxzzz_0_xz_0,  \
                             g_0_xxxzzz_0_xz_1,  \
                             g_0_xxxzzz_0_yy_0,  \
                             g_0_xxxzzz_0_yy_1,  \
                             g_0_xxxzzz_0_yz_0,  \
                             g_0_xxxzzz_0_yz_1,  \
                             g_0_xxxzzz_0_z_1,   \
                             g_0_xxxzzz_0_zz_0,  \
                             g_0_xxxzzz_0_zz_1,  \
                             g_0_xxzzz_0_xz_0,   \
                             g_0_xxzzz_0_xz_1,   \
                             g_0_xxzzz_0_yy_0,   \
                             g_0_xxzzz_0_yy_1,   \
                             g_0_xxzzz_0_yz_0,   \
                             g_0_xxzzz_0_yz_1,   \
                             g_0_xxzzz_0_zz_0,   \
                             g_0_xxzzz_0_zz_1,   \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzz_0_xx_0[i] =
            2.0 * g_0_xxxxz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xx_1[i] * fti_ab_0 + g_0_xxxxzz_0_xx_0[i] * pb_z + g_0_xxxxzz_0_xx_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xy_0[i] =
            2.0 * g_0_xxxxz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xy_0[i] * pb_z + g_0_xxxxzz_0_xy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xz_0[i] = 3.0 * g_0_xxzzz_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xz_1[i] * fti_ab_0 + g_0_xxxzzz_0_z_1[i] * fi_abcd_0 +
                                g_0_xxxzzz_0_xz_0[i] * pb_x + g_0_xxxzzz_0_xz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yy_0[i] =
            3.0 * g_0_xxzzz_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yy_1[i] * fti_ab_0 + g_0_xxxzzz_0_yy_0[i] * pb_x + g_0_xxxzzz_0_yy_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yz_0[i] =
            3.0 * g_0_xxzzz_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yz_0[i] * pb_x + g_0_xxxzzz_0_yz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_zz_0[i] =
            3.0 * g_0_xxzzz_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_zz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zz_0[i] * pb_x + g_0_xxxzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 60-66 components of targeted buffer : SKSD

    auto g_0_xxxyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 60);

    auto g_0_xxxyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 61);

    auto g_0_xxxyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 62);

    auto g_0_xxxyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 63);

    auto g_0_xxxyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 64);

    auto g_0_xxxyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 65);

#pragma omp simd aligned(g_0_xxxyy_0_xx_0,       \
                             g_0_xxxyy_0_xx_1,   \
                             g_0_xxxyy_0_xz_0,   \
                             g_0_xxxyy_0_xz_1,   \
                             g_0_xxxyyy_0_xx_0,  \
                             g_0_xxxyyy_0_xx_1,  \
                             g_0_xxxyyy_0_xz_0,  \
                             g_0_xxxyyy_0_xz_1,  \
                             g_0_xxxyyyy_0_xx_0, \
                             g_0_xxxyyyy_0_xy_0, \
                             g_0_xxxyyyy_0_xz_0, \
                             g_0_xxxyyyy_0_yy_0, \
                             g_0_xxxyyyy_0_yz_0, \
                             g_0_xxxyyyy_0_zz_0, \
                             g_0_xxyyyy_0_xy_0,  \
                             g_0_xxyyyy_0_xy_1,  \
                             g_0_xxyyyy_0_y_1,   \
                             g_0_xxyyyy_0_yy_0,  \
                             g_0_xxyyyy_0_yy_1,  \
                             g_0_xxyyyy_0_yz_0,  \
                             g_0_xxyyyy_0_yz_1,  \
                             g_0_xxyyyy_0_zz_0,  \
                             g_0_xxyyyy_0_zz_1,  \
                             g_0_xyyyy_0_xy_0,   \
                             g_0_xyyyy_0_xy_1,   \
                             g_0_xyyyy_0_yy_0,   \
                             g_0_xyyyy_0_yy_1,   \
                             g_0_xyyyy_0_yz_0,   \
                             g_0_xyyyy_0_yz_1,   \
                             g_0_xyyyy_0_zz_0,   \
                             g_0_xyyyy_0_zz_1,   \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyy_0_xx_0[i] =
            3.0 * g_0_xxxyy_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xx_1[i] * fti_ab_0 + g_0_xxxyyy_0_xx_0[i] * pb_y + g_0_xxxyyy_0_xx_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xy_0[i] = 2.0 * g_0_xyyyy_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xy_1[i] * fti_ab_0 + g_0_xxyyyy_0_y_1[i] * fi_abcd_0 +
                                g_0_xxyyyy_0_xy_0[i] * pb_x + g_0_xxyyyy_0_xy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xz_0[i] =
            3.0 * g_0_xxxyy_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xz_0[i] * pb_y + g_0_xxxyyy_0_xz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_yy_0[i] =
            2.0 * g_0_xyyyy_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yy_0[i] * pb_x + g_0_xxyyyy_0_yy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yz_0[i] =
            2.0 * g_0_xyyyy_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yz_0[i] * pb_x + g_0_xxyyyy_0_yz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_zz_0[i] =
            2.0 * g_0_xyyyy_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_zz_1[i] * fti_ab_0 + g_0_xxyyyy_0_zz_0[i] * pb_x + g_0_xxyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 66-72 components of targeted buffer : SKSD

    auto g_0_xxxyyyz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 66);

    auto g_0_xxxyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 67);

    auto g_0_xxxyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 68);

    auto g_0_xxxyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 69);

    auto g_0_xxxyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 70);

    auto g_0_xxxyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 71);

#pragma omp simd aligned(g_0_xxxyyy_0_x_1,       \
                             g_0_xxxyyy_0_xx_0,  \
                             g_0_xxxyyy_0_xx_1,  \
                             g_0_xxxyyy_0_xy_0,  \
                             g_0_xxxyyy_0_xy_1,  \
                             g_0_xxxyyy_0_xz_0,  \
                             g_0_xxxyyy_0_xz_1,  \
                             g_0_xxxyyy_0_y_1,   \
                             g_0_xxxyyy_0_yy_0,  \
                             g_0_xxxyyy_0_yy_1,  \
                             g_0_xxxyyy_0_yz_0,  \
                             g_0_xxxyyy_0_yz_1,  \
                             g_0_xxxyyy_0_z_1,   \
                             g_0_xxxyyy_0_zz_0,  \
                             g_0_xxxyyy_0_zz_1,  \
                             g_0_xxxyyyz_0_xx_0, \
                             g_0_xxxyyyz_0_xy_0, \
                             g_0_xxxyyyz_0_xz_0, \
                             g_0_xxxyyyz_0_yy_0, \
                             g_0_xxxyyyz_0_yz_0, \
                             g_0_xxxyyyz_0_zz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyz_0_xx_0[i] = g_0_xxxyyy_0_xx_0[i] * pb_z + g_0_xxxyyy_0_xx_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xy_0[i] = g_0_xxxyyy_0_xy_0[i] * pb_z + g_0_xxxyyy_0_xy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xz_0[i] = g_0_xxxyyy_0_x_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xz_0[i] * pb_z + g_0_xxxyyy_0_xz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yy_0[i] = g_0_xxxyyy_0_yy_0[i] * pb_z + g_0_xxxyyy_0_yy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yz_0[i] = g_0_xxxyyy_0_y_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yz_0[i] * pb_z + g_0_xxxyyy_0_yz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_zz_0[i] = 2.0 * g_0_xxxyyy_0_z_1[i] * fi_abcd_0 + g_0_xxxyyy_0_zz_0[i] * pb_z + g_0_xxxyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 72-78 components of targeted buffer : SKSD

    auto g_0_xxxyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 72);

    auto g_0_xxxyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 73);

    auto g_0_xxxyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 74);

    auto g_0_xxxyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 75);

    auto g_0_xxxyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 76);

    auto g_0_xxxyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 77);

#pragma omp simd aligned(g_0_xxxyy_0_xy_0,       \
                             g_0_xxxyy_0_xy_1,   \
                             g_0_xxxyyz_0_xy_0,  \
                             g_0_xxxyyz_0_xy_1,  \
                             g_0_xxxyyzz_0_xx_0, \
                             g_0_xxxyyzz_0_xy_0, \
                             g_0_xxxyyzz_0_xz_0, \
                             g_0_xxxyyzz_0_yy_0, \
                             g_0_xxxyyzz_0_yz_0, \
                             g_0_xxxyyzz_0_zz_0, \
                             g_0_xxxyzz_0_xx_0,  \
                             g_0_xxxyzz_0_xx_1,  \
                             g_0_xxxyzz_0_xz_0,  \
                             g_0_xxxyzz_0_xz_1,  \
                             g_0_xxxzz_0_xx_0,   \
                             g_0_xxxzz_0_xx_1,   \
                             g_0_xxxzz_0_xz_0,   \
                             g_0_xxxzz_0_xz_1,   \
                             g_0_xxyyzz_0_yy_0,  \
                             g_0_xxyyzz_0_yy_1,  \
                             g_0_xxyyzz_0_yz_0,  \
                             g_0_xxyyzz_0_yz_1,  \
                             g_0_xxyyzz_0_zz_0,  \
                             g_0_xxyyzz_0_zz_1,  \
                             g_0_xyyzz_0_yy_0,   \
                             g_0_xyyzz_0_yy_1,   \
                             g_0_xyyzz_0_yz_0,   \
                             g_0_xyyzz_0_yz_1,   \
                             g_0_xyyzz_0_zz_0,   \
                             g_0_xyyzz_0_zz_1,   \
                             wp_x,               \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyzz_0_xx_0[i] =
            g_0_xxxzz_0_xx_0[i] * fi_ab_0 - g_0_xxxzz_0_xx_1[i] * fti_ab_0 + g_0_xxxyzz_0_xx_0[i] * pb_y + g_0_xxxyzz_0_xx_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xy_0[i] =
            g_0_xxxyy_0_xy_0[i] * fi_ab_0 - g_0_xxxyy_0_xy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xy_0[i] * pb_z + g_0_xxxyyz_0_xy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xz_0[i] =
            g_0_xxxzz_0_xz_0[i] * fi_ab_0 - g_0_xxxzz_0_xz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xz_0[i] * pb_y + g_0_xxxyzz_0_xz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_yy_0[i] =
            2.0 * g_0_xyyzz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yy_1[i] * fti_ab_0 + g_0_xxyyzz_0_yy_0[i] * pb_x + g_0_xxyyzz_0_yy_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yz_0[i] =
            2.0 * g_0_xyyzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yz_0[i] * pb_x + g_0_xxyyzz_0_yz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_zz_0[i] =
            2.0 * g_0_xyyzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_zz_1[i] * fti_ab_0 + g_0_xxyyzz_0_zz_0[i] * pb_x + g_0_xxyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 78-84 components of targeted buffer : SKSD

    auto g_0_xxxyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 78);

    auto g_0_xxxyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 79);

    auto g_0_xxxyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 80);

    auto g_0_xxxyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 81);

    auto g_0_xxxyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 82);

    auto g_0_xxxyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 83);

#pragma omp simd aligned(g_0_xxxyzzz_0_xx_0,     \
                             g_0_xxxyzzz_0_xy_0, \
                             g_0_xxxyzzz_0_xz_0, \
                             g_0_xxxyzzz_0_yy_0, \
                             g_0_xxxyzzz_0_yz_0, \
                             g_0_xxxyzzz_0_zz_0, \
                             g_0_xxxzzz_0_x_1,   \
                             g_0_xxxzzz_0_xx_0,  \
                             g_0_xxxzzz_0_xx_1,  \
                             g_0_xxxzzz_0_xy_0,  \
                             g_0_xxxzzz_0_xy_1,  \
                             g_0_xxxzzz_0_xz_0,  \
                             g_0_xxxzzz_0_xz_1,  \
                             g_0_xxxzzz_0_y_1,   \
                             g_0_xxxzzz_0_yy_0,  \
                             g_0_xxxzzz_0_yy_1,  \
                             g_0_xxxzzz_0_yz_0,  \
                             g_0_xxxzzz_0_yz_1,  \
                             g_0_xxxzzz_0_z_1,   \
                             g_0_xxxzzz_0_zz_0,  \
                             g_0_xxxzzz_0_zz_1,  \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzz_0_xx_0[i] = g_0_xxxzzz_0_xx_0[i] * pb_y + g_0_xxxzzz_0_xx_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xy_0[i] = g_0_xxxzzz_0_x_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xy_0[i] * pb_y + g_0_xxxzzz_0_xy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xz_0[i] = g_0_xxxzzz_0_xz_0[i] * pb_y + g_0_xxxzzz_0_xz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yy_0[i] = 2.0 * g_0_xxxzzz_0_y_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yy_0[i] * pb_y + g_0_xxxzzz_0_yy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yz_0[i] = g_0_xxxzzz_0_z_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yz_0[i] * pb_y + g_0_xxxzzz_0_yz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_zz_0[i] = g_0_xxxzzz_0_zz_0[i] * pb_y + g_0_xxxzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 84-90 components of targeted buffer : SKSD

    auto g_0_xxxzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 84);

    auto g_0_xxxzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 85);

    auto g_0_xxxzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 86);

    auto g_0_xxxzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 87);

    auto g_0_xxxzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 88);

    auto g_0_xxxzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 89);

#pragma omp simd aligned(g_0_xxxzz_0_xx_0,       \
                             g_0_xxxzz_0_xx_1,   \
                             g_0_xxxzz_0_xy_0,   \
                             g_0_xxxzz_0_xy_1,   \
                             g_0_xxxzzz_0_xx_0,  \
                             g_0_xxxzzz_0_xx_1,  \
                             g_0_xxxzzz_0_xy_0,  \
                             g_0_xxxzzz_0_xy_1,  \
                             g_0_xxxzzzz_0_xx_0, \
                             g_0_xxxzzzz_0_xy_0, \
                             g_0_xxxzzzz_0_xz_0, \
                             g_0_xxxzzzz_0_yy_0, \
                             g_0_xxxzzzz_0_yz_0, \
                             g_0_xxxzzzz_0_zz_0, \
                             g_0_xxzzzz_0_xz_0,  \
                             g_0_xxzzzz_0_xz_1,  \
                             g_0_xxzzzz_0_yy_0,  \
                             g_0_xxzzzz_0_yy_1,  \
                             g_0_xxzzzz_0_yz_0,  \
                             g_0_xxzzzz_0_yz_1,  \
                             g_0_xxzzzz_0_z_1,   \
                             g_0_xxzzzz_0_zz_0,  \
                             g_0_xxzzzz_0_zz_1,  \
                             g_0_xzzzz_0_xz_0,   \
                             g_0_xzzzz_0_xz_1,   \
                             g_0_xzzzz_0_yy_0,   \
                             g_0_xzzzz_0_yy_1,   \
                             g_0_xzzzz_0_yz_0,   \
                             g_0_xzzzz_0_yz_1,   \
                             g_0_xzzzz_0_zz_0,   \
                             g_0_xzzzz_0_zz_1,   \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzz_0_xx_0[i] =
            3.0 * g_0_xxxzz_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xx_1[i] * fti_ab_0 + g_0_xxxzzz_0_xx_0[i] * pb_z + g_0_xxxzzz_0_xx_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xy_0[i] =
            3.0 * g_0_xxxzz_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xy_0[i] * pb_z + g_0_xxxzzz_0_xy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xz_0[i] = 2.0 * g_0_xzzzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xz_1[i] * fti_ab_0 + g_0_xxzzzz_0_z_1[i] * fi_abcd_0 +
                                g_0_xxzzzz_0_xz_0[i] * pb_x + g_0_xxzzzz_0_xz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yy_0[i] =
            2.0 * g_0_xzzzz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yy_1[i] * fti_ab_0 + g_0_xxzzzz_0_yy_0[i] * pb_x + g_0_xxzzzz_0_yy_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yz_0[i] =
            2.0 * g_0_xzzzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yz_0[i] * pb_x + g_0_xxzzzz_0_yz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_zz_0[i] =
            2.0 * g_0_xzzzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_zz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zz_0[i] * pb_x + g_0_xxzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 90-96 components of targeted buffer : SKSD

    auto g_0_xxyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 90);

    auto g_0_xxyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 91);

    auto g_0_xxyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 92);

    auto g_0_xxyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 93);

    auto g_0_xxyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 94);

    auto g_0_xxyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 95);

#pragma omp simd aligned(g_0_xxyyy_0_xx_0,       \
                             g_0_xxyyy_0_xx_1,   \
                             g_0_xxyyy_0_xz_0,   \
                             g_0_xxyyy_0_xz_1,   \
                             g_0_xxyyyy_0_xx_0,  \
                             g_0_xxyyyy_0_xx_1,  \
                             g_0_xxyyyy_0_xz_0,  \
                             g_0_xxyyyy_0_xz_1,  \
                             g_0_xxyyyyy_0_xx_0, \
                             g_0_xxyyyyy_0_xy_0, \
                             g_0_xxyyyyy_0_xz_0, \
                             g_0_xxyyyyy_0_yy_0, \
                             g_0_xxyyyyy_0_yz_0, \
                             g_0_xxyyyyy_0_zz_0, \
                             g_0_xyyyyy_0_xy_0,  \
                             g_0_xyyyyy_0_xy_1,  \
                             g_0_xyyyyy_0_y_1,   \
                             g_0_xyyyyy_0_yy_0,  \
                             g_0_xyyyyy_0_yy_1,  \
                             g_0_xyyyyy_0_yz_0,  \
                             g_0_xyyyyy_0_yz_1,  \
                             g_0_xyyyyy_0_zz_0,  \
                             g_0_xyyyyy_0_zz_1,  \
                             g_0_yyyyy_0_xy_0,   \
                             g_0_yyyyy_0_xy_1,   \
                             g_0_yyyyy_0_yy_0,   \
                             g_0_yyyyy_0_yy_1,   \
                             g_0_yyyyy_0_yz_0,   \
                             g_0_yyyyy_0_yz_1,   \
                             g_0_yyyyy_0_zz_0,   \
                             g_0_yyyyy_0_zz_1,   \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyy_0_xx_0[i] =
            4.0 * g_0_xxyyy_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xx_1[i] * fti_ab_0 + g_0_xxyyyy_0_xx_0[i] * pb_y + g_0_xxyyyy_0_xx_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xy_0[i] = g_0_yyyyy_0_xy_0[i] * fi_ab_0 - g_0_yyyyy_0_xy_1[i] * fti_ab_0 + g_0_xyyyyy_0_y_1[i] * fi_abcd_0 +
                                g_0_xyyyyy_0_xy_0[i] * pb_x + g_0_xyyyyy_0_xy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xz_0[i] =
            4.0 * g_0_xxyyy_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xz_0[i] * pb_y + g_0_xxyyyy_0_xz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_yy_0[i] =
            g_0_yyyyy_0_yy_0[i] * fi_ab_0 - g_0_yyyyy_0_yy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yy_0[i] * pb_x + g_0_xyyyyy_0_yy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yz_0[i] =
            g_0_yyyyy_0_yz_0[i] * fi_ab_0 - g_0_yyyyy_0_yz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yz_0[i] * pb_x + g_0_xyyyyy_0_yz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_zz_0[i] =
            g_0_yyyyy_0_zz_0[i] * fi_ab_0 - g_0_yyyyy_0_zz_1[i] * fti_ab_0 + g_0_xyyyyy_0_zz_0[i] * pb_x + g_0_xyyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 96-102 components of targeted buffer : SKSD

    auto g_0_xxyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 96);

    auto g_0_xxyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 97);

    auto g_0_xxyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 98);

    auto g_0_xxyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 99);

    auto g_0_xxyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 100);

    auto g_0_xxyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 101);

#pragma omp simd aligned(g_0_xxyyyy_0_x_1,       \
                             g_0_xxyyyy_0_xx_0,  \
                             g_0_xxyyyy_0_xx_1,  \
                             g_0_xxyyyy_0_xy_0,  \
                             g_0_xxyyyy_0_xy_1,  \
                             g_0_xxyyyy_0_xz_0,  \
                             g_0_xxyyyy_0_xz_1,  \
                             g_0_xxyyyy_0_y_1,   \
                             g_0_xxyyyy_0_yy_0,  \
                             g_0_xxyyyy_0_yy_1,  \
                             g_0_xxyyyy_0_yz_0,  \
                             g_0_xxyyyy_0_yz_1,  \
                             g_0_xxyyyy_0_z_1,   \
                             g_0_xxyyyy_0_zz_0,  \
                             g_0_xxyyyy_0_zz_1,  \
                             g_0_xxyyyyz_0_xx_0, \
                             g_0_xxyyyyz_0_xy_0, \
                             g_0_xxyyyyz_0_xz_0, \
                             g_0_xxyyyyz_0_yy_0, \
                             g_0_xxyyyyz_0_yz_0, \
                             g_0_xxyyyyz_0_zz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyz_0_xx_0[i] = g_0_xxyyyy_0_xx_0[i] * pb_z + g_0_xxyyyy_0_xx_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xy_0[i] = g_0_xxyyyy_0_xy_0[i] * pb_z + g_0_xxyyyy_0_xy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xz_0[i] = g_0_xxyyyy_0_x_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xz_0[i] * pb_z + g_0_xxyyyy_0_xz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yy_0[i] = g_0_xxyyyy_0_yy_0[i] * pb_z + g_0_xxyyyy_0_yy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yz_0[i] = g_0_xxyyyy_0_y_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yz_0[i] * pb_z + g_0_xxyyyy_0_yz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_zz_0[i] = 2.0 * g_0_xxyyyy_0_z_1[i] * fi_abcd_0 + g_0_xxyyyy_0_zz_0[i] * pb_z + g_0_xxyyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 102-108 components of targeted buffer : SKSD

    auto g_0_xxyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 102);

    auto g_0_xxyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 103);

    auto g_0_xxyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 104);

    auto g_0_xxyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 105);

    auto g_0_xxyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 106);

    auto g_0_xxyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 107);

#pragma omp simd aligned(g_0_xxyyy_0_xy_0,       \
                             g_0_xxyyy_0_xy_1,   \
                             g_0_xxyyyz_0_xy_0,  \
                             g_0_xxyyyz_0_xy_1,  \
                             g_0_xxyyyzz_0_xx_0, \
                             g_0_xxyyyzz_0_xy_0, \
                             g_0_xxyyyzz_0_xz_0, \
                             g_0_xxyyyzz_0_yy_0, \
                             g_0_xxyyyzz_0_yz_0, \
                             g_0_xxyyyzz_0_zz_0, \
                             g_0_xxyyzz_0_xx_0,  \
                             g_0_xxyyzz_0_xx_1,  \
                             g_0_xxyyzz_0_xz_0,  \
                             g_0_xxyyzz_0_xz_1,  \
                             g_0_xxyzz_0_xx_0,   \
                             g_0_xxyzz_0_xx_1,   \
                             g_0_xxyzz_0_xz_0,   \
                             g_0_xxyzz_0_xz_1,   \
                             g_0_xyyyzz_0_yy_0,  \
                             g_0_xyyyzz_0_yy_1,  \
                             g_0_xyyyzz_0_yz_0,  \
                             g_0_xyyyzz_0_yz_1,  \
                             g_0_xyyyzz_0_zz_0,  \
                             g_0_xyyyzz_0_zz_1,  \
                             g_0_yyyzz_0_yy_0,   \
                             g_0_yyyzz_0_yy_1,   \
                             g_0_yyyzz_0_yz_0,   \
                             g_0_yyyzz_0_yz_1,   \
                             g_0_yyyzz_0_zz_0,   \
                             g_0_yyyzz_0_zz_1,   \
                             wp_x,               \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyzz_0_xx_0[i] =
            2.0 * g_0_xxyzz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xx_1[i] * fti_ab_0 + g_0_xxyyzz_0_xx_0[i] * pb_y + g_0_xxyyzz_0_xx_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xy_0[i] =
            g_0_xxyyy_0_xy_0[i] * fi_ab_0 - g_0_xxyyy_0_xy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xy_0[i] * pb_z + g_0_xxyyyz_0_xy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xz_0[i] =
            2.0 * g_0_xxyzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xz_0[i] * pb_y + g_0_xxyyzz_0_xz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_yy_0[i] =
            g_0_yyyzz_0_yy_0[i] * fi_ab_0 - g_0_yyyzz_0_yy_1[i] * fti_ab_0 + g_0_xyyyzz_0_yy_0[i] * pb_x + g_0_xyyyzz_0_yy_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yz_0[i] =
            g_0_yyyzz_0_yz_0[i] * fi_ab_0 - g_0_yyyzz_0_yz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yz_0[i] * pb_x + g_0_xyyyzz_0_yz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_zz_0[i] =
            g_0_yyyzz_0_zz_0[i] * fi_ab_0 - g_0_yyyzz_0_zz_1[i] * fti_ab_0 + g_0_xyyyzz_0_zz_0[i] * pb_x + g_0_xyyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 108-114 components of targeted buffer : SKSD

    auto g_0_xxyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 108);

    auto g_0_xxyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 109);

    auto g_0_xxyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 110);

    auto g_0_xxyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 111);

    auto g_0_xxyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 112);

    auto g_0_xxyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 113);

#pragma omp simd aligned(g_0_xxyyz_0_xy_0,       \
                             g_0_xxyyz_0_xy_1,   \
                             g_0_xxyyzz_0_xy_0,  \
                             g_0_xxyyzz_0_xy_1,  \
                             g_0_xxyyzzz_0_xx_0, \
                             g_0_xxyyzzz_0_xy_0, \
                             g_0_xxyyzzz_0_xz_0, \
                             g_0_xxyyzzz_0_yy_0, \
                             g_0_xxyyzzz_0_yz_0, \
                             g_0_xxyyzzz_0_zz_0, \
                             g_0_xxyzzz_0_xx_0,  \
                             g_0_xxyzzz_0_xx_1,  \
                             g_0_xxyzzz_0_xz_0,  \
                             g_0_xxyzzz_0_xz_1,  \
                             g_0_xxzzz_0_xx_0,   \
                             g_0_xxzzz_0_xx_1,   \
                             g_0_xxzzz_0_xz_0,   \
                             g_0_xxzzz_0_xz_1,   \
                             g_0_xyyzzz_0_yy_0,  \
                             g_0_xyyzzz_0_yy_1,  \
                             g_0_xyyzzz_0_yz_0,  \
                             g_0_xyyzzz_0_yz_1,  \
                             g_0_xyyzzz_0_zz_0,  \
                             g_0_xyyzzz_0_zz_1,  \
                             g_0_yyzzz_0_yy_0,   \
                             g_0_yyzzz_0_yy_1,   \
                             g_0_yyzzz_0_yz_0,   \
                             g_0_yyzzz_0_yz_1,   \
                             g_0_yyzzz_0_zz_0,   \
                             g_0_yyzzz_0_zz_1,   \
                             wp_x,               \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyzzz_0_xx_0[i] =
            g_0_xxzzz_0_xx_0[i] * fi_ab_0 - g_0_xxzzz_0_xx_1[i] * fti_ab_0 + g_0_xxyzzz_0_xx_0[i] * pb_y + g_0_xxyzzz_0_xx_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xy_0[i] =
            2.0 * g_0_xxyyz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xy_0[i] * pb_z + g_0_xxyyzz_0_xy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xz_0[i] =
            g_0_xxzzz_0_xz_0[i] * fi_ab_0 - g_0_xxzzz_0_xz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xz_0[i] * pb_y + g_0_xxyzzz_0_xz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_yy_0[i] =
            g_0_yyzzz_0_yy_0[i] * fi_ab_0 - g_0_yyzzz_0_yy_1[i] * fti_ab_0 + g_0_xyyzzz_0_yy_0[i] * pb_x + g_0_xyyzzz_0_yy_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yz_0[i] =
            g_0_yyzzz_0_yz_0[i] * fi_ab_0 - g_0_yyzzz_0_yz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yz_0[i] * pb_x + g_0_xyyzzz_0_yz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_zz_0[i] =
            g_0_yyzzz_0_zz_0[i] * fi_ab_0 - g_0_yyzzz_0_zz_1[i] * fti_ab_0 + g_0_xyyzzz_0_zz_0[i] * pb_x + g_0_xyyzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 114-120 components of targeted buffer : SKSD

    auto g_0_xxyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 114);

    auto g_0_xxyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 115);

    auto g_0_xxyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 116);

    auto g_0_xxyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 117);

    auto g_0_xxyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 118);

    auto g_0_xxyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 119);

#pragma omp simd aligned(g_0_xxyzzzz_0_xx_0,     \
                             g_0_xxyzzzz_0_xy_0, \
                             g_0_xxyzzzz_0_xz_0, \
                             g_0_xxyzzzz_0_yy_0, \
                             g_0_xxyzzzz_0_yz_0, \
                             g_0_xxyzzzz_0_zz_0, \
                             g_0_xxzzzz_0_x_1,   \
                             g_0_xxzzzz_0_xx_0,  \
                             g_0_xxzzzz_0_xx_1,  \
                             g_0_xxzzzz_0_xy_0,  \
                             g_0_xxzzzz_0_xy_1,  \
                             g_0_xxzzzz_0_xz_0,  \
                             g_0_xxzzzz_0_xz_1,  \
                             g_0_xxzzzz_0_y_1,   \
                             g_0_xxzzzz_0_yy_0,  \
                             g_0_xxzzzz_0_yy_1,  \
                             g_0_xxzzzz_0_yz_0,  \
                             g_0_xxzzzz_0_yz_1,  \
                             g_0_xxzzzz_0_z_1,   \
                             g_0_xxzzzz_0_zz_0,  \
                             g_0_xxzzzz_0_zz_1,  \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzz_0_xx_0[i] = g_0_xxzzzz_0_xx_0[i] * pb_y + g_0_xxzzzz_0_xx_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xy_0[i] = g_0_xxzzzz_0_x_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xy_0[i] * pb_y + g_0_xxzzzz_0_xy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xz_0[i] = g_0_xxzzzz_0_xz_0[i] * pb_y + g_0_xxzzzz_0_xz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yy_0[i] = 2.0 * g_0_xxzzzz_0_y_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yy_0[i] * pb_y + g_0_xxzzzz_0_yy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yz_0[i] = g_0_xxzzzz_0_z_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yz_0[i] * pb_y + g_0_xxzzzz_0_yz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_zz_0[i] = g_0_xxzzzz_0_zz_0[i] * pb_y + g_0_xxzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 120-126 components of targeted buffer : SKSD

    auto g_0_xxzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 120);

    auto g_0_xxzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 121);

    auto g_0_xxzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 122);

    auto g_0_xxzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 123);

    auto g_0_xxzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 124);

    auto g_0_xxzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 125);

#pragma omp simd aligned(g_0_xxzzz_0_xx_0,       \
                             g_0_xxzzz_0_xx_1,   \
                             g_0_xxzzz_0_xy_0,   \
                             g_0_xxzzz_0_xy_1,   \
                             g_0_xxzzzz_0_xx_0,  \
                             g_0_xxzzzz_0_xx_1,  \
                             g_0_xxzzzz_0_xy_0,  \
                             g_0_xxzzzz_0_xy_1,  \
                             g_0_xxzzzzz_0_xx_0, \
                             g_0_xxzzzzz_0_xy_0, \
                             g_0_xxzzzzz_0_xz_0, \
                             g_0_xxzzzzz_0_yy_0, \
                             g_0_xxzzzzz_0_yz_0, \
                             g_0_xxzzzzz_0_zz_0, \
                             g_0_xzzzzz_0_xz_0,  \
                             g_0_xzzzzz_0_xz_1,  \
                             g_0_xzzzzz_0_yy_0,  \
                             g_0_xzzzzz_0_yy_1,  \
                             g_0_xzzzzz_0_yz_0,  \
                             g_0_xzzzzz_0_yz_1,  \
                             g_0_xzzzzz_0_z_1,   \
                             g_0_xzzzzz_0_zz_0,  \
                             g_0_xzzzzz_0_zz_1,  \
                             g_0_zzzzz_0_xz_0,   \
                             g_0_zzzzz_0_xz_1,   \
                             g_0_zzzzz_0_yy_0,   \
                             g_0_zzzzz_0_yy_1,   \
                             g_0_zzzzz_0_yz_0,   \
                             g_0_zzzzz_0_yz_1,   \
                             g_0_zzzzz_0_zz_0,   \
                             g_0_zzzzz_0_zz_1,   \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzz_0_xx_0[i] =
            4.0 * g_0_xxzzz_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xx_1[i] * fti_ab_0 + g_0_xxzzzz_0_xx_0[i] * pb_z + g_0_xxzzzz_0_xx_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xy_0[i] =
            4.0 * g_0_xxzzz_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xy_0[i] * pb_z + g_0_xxzzzz_0_xy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xz_0[i] = g_0_zzzzz_0_xz_0[i] * fi_ab_0 - g_0_zzzzz_0_xz_1[i] * fti_ab_0 + g_0_xzzzzz_0_z_1[i] * fi_abcd_0 +
                                g_0_xzzzzz_0_xz_0[i] * pb_x + g_0_xzzzzz_0_xz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yy_0[i] =
            g_0_zzzzz_0_yy_0[i] * fi_ab_0 - g_0_zzzzz_0_yy_1[i] * fti_ab_0 + g_0_xzzzzz_0_yy_0[i] * pb_x + g_0_xzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yz_0[i] =
            g_0_zzzzz_0_yz_0[i] * fi_ab_0 - g_0_zzzzz_0_yz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yz_0[i] * pb_x + g_0_xzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_zz_0[i] =
            g_0_zzzzz_0_zz_0[i] * fi_ab_0 - g_0_zzzzz_0_zz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zz_0[i] * pb_x + g_0_xzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 126-132 components of targeted buffer : SKSD

    auto g_0_xyyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 126);

    auto g_0_xyyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 127);

    auto g_0_xyyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 128);

    auto g_0_xyyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 129);

    auto g_0_xyyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 130);

    auto g_0_xyyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 131);

#pragma omp simd aligned(g_0_xyyyyyy_0_xx_0,     \
                             g_0_xyyyyyy_0_xy_0, \
                             g_0_xyyyyyy_0_xz_0, \
                             g_0_xyyyyyy_0_yy_0, \
                             g_0_xyyyyyy_0_yz_0, \
                             g_0_xyyyyyy_0_zz_0, \
                             g_0_yyyyyy_0_x_1,   \
                             g_0_yyyyyy_0_xx_0,  \
                             g_0_yyyyyy_0_xx_1,  \
                             g_0_yyyyyy_0_xy_0,  \
                             g_0_yyyyyy_0_xy_1,  \
                             g_0_yyyyyy_0_xz_0,  \
                             g_0_yyyyyy_0_xz_1,  \
                             g_0_yyyyyy_0_y_1,   \
                             g_0_yyyyyy_0_yy_0,  \
                             g_0_yyyyyy_0_yy_1,  \
                             g_0_yyyyyy_0_yz_0,  \
                             g_0_yyyyyy_0_yz_1,  \
                             g_0_yyyyyy_0_z_1,   \
                             g_0_yyyyyy_0_zz_0,  \
                             g_0_yyyyyy_0_zz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyy_0_xx_0[i] = 2.0 * g_0_yyyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xx_0[i] * pb_x + g_0_yyyyyy_0_xx_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xy_0[i] = g_0_yyyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xy_0[i] * pb_x + g_0_yyyyyy_0_xy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xz_0[i] = g_0_yyyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xz_0[i] * pb_x + g_0_yyyyyy_0_xz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yy_0[i] = g_0_yyyyyy_0_yy_0[i] * pb_x + g_0_yyyyyy_0_yy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yz_0[i] = g_0_yyyyyy_0_yz_0[i] * pb_x + g_0_yyyyyy_0_yz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_zz_0[i] = g_0_yyyyyy_0_zz_0[i] * pb_x + g_0_yyyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 132-138 components of targeted buffer : SKSD

    auto g_0_xyyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 132);

    auto g_0_xyyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 133);

    auto g_0_xyyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 134);

    auto g_0_xyyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 135);

    auto g_0_xyyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 136);

    auto g_0_xyyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 137);

#pragma omp simd aligned(g_0_xyyyyy_0_xx_0,      \
                             g_0_xyyyyy_0_xx_1,  \
                             g_0_xyyyyy_0_xy_0,  \
                             g_0_xyyyyy_0_xy_1,  \
                             g_0_xyyyyyz_0_xx_0, \
                             g_0_xyyyyyz_0_xy_0, \
                             g_0_xyyyyyz_0_xz_0, \
                             g_0_xyyyyyz_0_yy_0, \
                             g_0_xyyyyyz_0_yz_0, \
                             g_0_xyyyyyz_0_zz_0, \
                             g_0_yyyyyz_0_xz_0,  \
                             g_0_yyyyyz_0_xz_1,  \
                             g_0_yyyyyz_0_yy_0,  \
                             g_0_yyyyyz_0_yy_1,  \
                             g_0_yyyyyz_0_yz_0,  \
                             g_0_yyyyyz_0_yz_1,  \
                             g_0_yyyyyz_0_z_1,   \
                             g_0_yyyyyz_0_zz_0,  \
                             g_0_yyyyyz_0_zz_1,  \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyz_0_xx_0[i] = g_0_xyyyyy_0_xx_0[i] * pb_z + g_0_xyyyyy_0_xx_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xy_0[i] = g_0_xyyyyy_0_xy_0[i] * pb_z + g_0_xyyyyy_0_xy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xz_0[i] = g_0_yyyyyz_0_z_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xz_0[i] * pb_x + g_0_yyyyyz_0_xz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yy_0[i] = g_0_yyyyyz_0_yy_0[i] * pb_x + g_0_yyyyyz_0_yy_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yz_0[i] = g_0_yyyyyz_0_yz_0[i] * pb_x + g_0_yyyyyz_0_yz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_zz_0[i] = g_0_yyyyyz_0_zz_0[i] * pb_x + g_0_yyyyyz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 138-144 components of targeted buffer : SKSD

    auto g_0_xyyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 138);

    auto g_0_xyyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 139);

    auto g_0_xyyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 140);

    auto g_0_xyyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 141);

    auto g_0_xyyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 142);

    auto g_0_xyyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 143);

#pragma omp simd aligned(g_0_xyyyyzz_0_xx_0,     \
                             g_0_xyyyyzz_0_xy_0, \
                             g_0_xyyyyzz_0_xz_0, \
                             g_0_xyyyyzz_0_yy_0, \
                             g_0_xyyyyzz_0_yz_0, \
                             g_0_xyyyyzz_0_zz_0, \
                             g_0_yyyyzz_0_x_1,   \
                             g_0_yyyyzz_0_xx_0,  \
                             g_0_yyyyzz_0_xx_1,  \
                             g_0_yyyyzz_0_xy_0,  \
                             g_0_yyyyzz_0_xy_1,  \
                             g_0_yyyyzz_0_xz_0,  \
                             g_0_yyyyzz_0_xz_1,  \
                             g_0_yyyyzz_0_y_1,   \
                             g_0_yyyyzz_0_yy_0,  \
                             g_0_yyyyzz_0_yy_1,  \
                             g_0_yyyyzz_0_yz_0,  \
                             g_0_yyyyzz_0_yz_1,  \
                             g_0_yyyyzz_0_z_1,   \
                             g_0_yyyyzz_0_zz_0,  \
                             g_0_yyyyzz_0_zz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzz_0_xx_0[i] = 2.0 * g_0_yyyyzz_0_x_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xx_0[i] * pb_x + g_0_yyyyzz_0_xx_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xy_0[i] = g_0_yyyyzz_0_y_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xy_0[i] * pb_x + g_0_yyyyzz_0_xy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xz_0[i] = g_0_yyyyzz_0_z_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xz_0[i] * pb_x + g_0_yyyyzz_0_xz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yy_0[i] = g_0_yyyyzz_0_yy_0[i] * pb_x + g_0_yyyyzz_0_yy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yz_0[i] = g_0_yyyyzz_0_yz_0[i] * pb_x + g_0_yyyyzz_0_yz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_zz_0[i] = g_0_yyyyzz_0_zz_0[i] * pb_x + g_0_yyyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 144-150 components of targeted buffer : SKSD

    auto g_0_xyyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 144);

    auto g_0_xyyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 145);

    auto g_0_xyyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 146);

    auto g_0_xyyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 147);

    auto g_0_xyyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 148);

    auto g_0_xyyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 149);

#pragma omp simd aligned(g_0_xyyyzzz_0_xx_0,     \
                             g_0_xyyyzzz_0_xy_0, \
                             g_0_xyyyzzz_0_xz_0, \
                             g_0_xyyyzzz_0_yy_0, \
                             g_0_xyyyzzz_0_yz_0, \
                             g_0_xyyyzzz_0_zz_0, \
                             g_0_yyyzzz_0_x_1,   \
                             g_0_yyyzzz_0_xx_0,  \
                             g_0_yyyzzz_0_xx_1,  \
                             g_0_yyyzzz_0_xy_0,  \
                             g_0_yyyzzz_0_xy_1,  \
                             g_0_yyyzzz_0_xz_0,  \
                             g_0_yyyzzz_0_xz_1,  \
                             g_0_yyyzzz_0_y_1,   \
                             g_0_yyyzzz_0_yy_0,  \
                             g_0_yyyzzz_0_yy_1,  \
                             g_0_yyyzzz_0_yz_0,  \
                             g_0_yyyzzz_0_yz_1,  \
                             g_0_yyyzzz_0_z_1,   \
                             g_0_yyyzzz_0_zz_0,  \
                             g_0_yyyzzz_0_zz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzz_0_xx_0[i] = 2.0 * g_0_yyyzzz_0_x_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xx_0[i] * pb_x + g_0_yyyzzz_0_xx_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xy_0[i] = g_0_yyyzzz_0_y_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xy_0[i] * pb_x + g_0_yyyzzz_0_xy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xz_0[i] = g_0_yyyzzz_0_z_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xz_0[i] * pb_x + g_0_yyyzzz_0_xz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yy_0[i] = g_0_yyyzzz_0_yy_0[i] * pb_x + g_0_yyyzzz_0_yy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yz_0[i] = g_0_yyyzzz_0_yz_0[i] * pb_x + g_0_yyyzzz_0_yz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_zz_0[i] = g_0_yyyzzz_0_zz_0[i] * pb_x + g_0_yyyzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 150-156 components of targeted buffer : SKSD

    auto g_0_xyyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 150);

    auto g_0_xyyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 151);

    auto g_0_xyyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 152);

    auto g_0_xyyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 153);

    auto g_0_xyyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 154);

    auto g_0_xyyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 155);

#pragma omp simd aligned(g_0_xyyzzzz_0_xx_0,     \
                             g_0_xyyzzzz_0_xy_0, \
                             g_0_xyyzzzz_0_xz_0, \
                             g_0_xyyzzzz_0_yy_0, \
                             g_0_xyyzzzz_0_yz_0, \
                             g_0_xyyzzzz_0_zz_0, \
                             g_0_yyzzzz_0_x_1,   \
                             g_0_yyzzzz_0_xx_0,  \
                             g_0_yyzzzz_0_xx_1,  \
                             g_0_yyzzzz_0_xy_0,  \
                             g_0_yyzzzz_0_xy_1,  \
                             g_0_yyzzzz_0_xz_0,  \
                             g_0_yyzzzz_0_xz_1,  \
                             g_0_yyzzzz_0_y_1,   \
                             g_0_yyzzzz_0_yy_0,  \
                             g_0_yyzzzz_0_yy_1,  \
                             g_0_yyzzzz_0_yz_0,  \
                             g_0_yyzzzz_0_yz_1,  \
                             g_0_yyzzzz_0_z_1,   \
                             g_0_yyzzzz_0_zz_0,  \
                             g_0_yyzzzz_0_zz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzz_0_xx_0[i] = 2.0 * g_0_yyzzzz_0_x_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xx_0[i] * pb_x + g_0_yyzzzz_0_xx_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xy_0[i] = g_0_yyzzzz_0_y_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xy_0[i] * pb_x + g_0_yyzzzz_0_xy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xz_0[i] = g_0_yyzzzz_0_z_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xz_0[i] * pb_x + g_0_yyzzzz_0_xz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yy_0[i] = g_0_yyzzzz_0_yy_0[i] * pb_x + g_0_yyzzzz_0_yy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yz_0[i] = g_0_yyzzzz_0_yz_0[i] * pb_x + g_0_yyzzzz_0_yz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_zz_0[i] = g_0_yyzzzz_0_zz_0[i] * pb_x + g_0_yyzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 156-162 components of targeted buffer : SKSD

    auto g_0_xyzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 156);

    auto g_0_xyzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 157);

    auto g_0_xyzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 158);

    auto g_0_xyzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 159);

    auto g_0_xyzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 160);

    auto g_0_xyzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 161);

#pragma omp simd aligned(g_0_xyzzzzz_0_xx_0,     \
                             g_0_xyzzzzz_0_xy_0, \
                             g_0_xyzzzzz_0_xz_0, \
                             g_0_xyzzzzz_0_yy_0, \
                             g_0_xyzzzzz_0_yz_0, \
                             g_0_xyzzzzz_0_zz_0, \
                             g_0_xzzzzz_0_xx_0,  \
                             g_0_xzzzzz_0_xx_1,  \
                             g_0_xzzzzz_0_xz_0,  \
                             g_0_xzzzzz_0_xz_1,  \
                             g_0_yzzzzz_0_xy_0,  \
                             g_0_yzzzzz_0_xy_1,  \
                             g_0_yzzzzz_0_y_1,   \
                             g_0_yzzzzz_0_yy_0,  \
                             g_0_yzzzzz_0_yy_1,  \
                             g_0_yzzzzz_0_yz_0,  \
                             g_0_yzzzzz_0_yz_1,  \
                             g_0_yzzzzz_0_zz_0,  \
                             g_0_yzzzzz_0_zz_1,  \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzz_0_xx_0[i] = g_0_xzzzzz_0_xx_0[i] * pb_y + g_0_xzzzzz_0_xx_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xy_0[i] = g_0_yzzzzz_0_y_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xy_0[i] * pb_x + g_0_yzzzzz_0_xy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xz_0[i] = g_0_xzzzzz_0_xz_0[i] * pb_y + g_0_xzzzzz_0_xz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_yy_0[i] = g_0_yzzzzz_0_yy_0[i] * pb_x + g_0_yzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yz_0[i] = g_0_yzzzzz_0_yz_0[i] * pb_x + g_0_yzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_zz_0[i] = g_0_yzzzzz_0_zz_0[i] * pb_x + g_0_yzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 162-168 components of targeted buffer : SKSD

    auto g_0_xzzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 162);

    auto g_0_xzzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 163);

    auto g_0_xzzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 164);

    auto g_0_xzzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 165);

    auto g_0_xzzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 166);

    auto g_0_xzzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 167);

#pragma omp simd aligned(g_0_xzzzzzz_0_xx_0,     \
                             g_0_xzzzzzz_0_xy_0, \
                             g_0_xzzzzzz_0_xz_0, \
                             g_0_xzzzzzz_0_yy_0, \
                             g_0_xzzzzzz_0_yz_0, \
                             g_0_xzzzzzz_0_zz_0, \
                             g_0_zzzzzz_0_x_1,   \
                             g_0_zzzzzz_0_xx_0,  \
                             g_0_zzzzzz_0_xx_1,  \
                             g_0_zzzzzz_0_xy_0,  \
                             g_0_zzzzzz_0_xy_1,  \
                             g_0_zzzzzz_0_xz_0,  \
                             g_0_zzzzzz_0_xz_1,  \
                             g_0_zzzzzz_0_y_1,   \
                             g_0_zzzzzz_0_yy_0,  \
                             g_0_zzzzzz_0_yy_1,  \
                             g_0_zzzzzz_0_yz_0,  \
                             g_0_zzzzzz_0_yz_1,  \
                             g_0_zzzzzz_0_z_1,   \
                             g_0_zzzzzz_0_zz_0,  \
                             g_0_zzzzzz_0_zz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzz_0_xx_0[i] = 2.0 * g_0_zzzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xx_0[i] * pb_x + g_0_zzzzzz_0_xx_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xy_0[i] = g_0_zzzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xy_0[i] * pb_x + g_0_zzzzzz_0_xy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xz_0[i] = g_0_zzzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xz_0[i] * pb_x + g_0_zzzzzz_0_xz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yy_0[i] = g_0_zzzzzz_0_yy_0[i] * pb_x + g_0_zzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yz_0[i] = g_0_zzzzzz_0_yz_0[i] * pb_x + g_0_zzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_zz_0[i] = g_0_zzzzzz_0_zz_0[i] * pb_x + g_0_zzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 168-174 components of targeted buffer : SKSD

    auto g_0_yyyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 168);

    auto g_0_yyyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 169);

    auto g_0_yyyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 170);

    auto g_0_yyyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 171);

    auto g_0_yyyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 172);

    auto g_0_yyyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 173);

#pragma omp simd aligned(g_0_yyyyy_0_xx_0,       \
                             g_0_yyyyy_0_xx_1,   \
                             g_0_yyyyy_0_xy_0,   \
                             g_0_yyyyy_0_xy_1,   \
                             g_0_yyyyy_0_xz_0,   \
                             g_0_yyyyy_0_xz_1,   \
                             g_0_yyyyy_0_yy_0,   \
                             g_0_yyyyy_0_yy_1,   \
                             g_0_yyyyy_0_yz_0,   \
                             g_0_yyyyy_0_yz_1,   \
                             g_0_yyyyy_0_zz_0,   \
                             g_0_yyyyy_0_zz_1,   \
                             g_0_yyyyyy_0_x_1,   \
                             g_0_yyyyyy_0_xx_0,  \
                             g_0_yyyyyy_0_xx_1,  \
                             g_0_yyyyyy_0_xy_0,  \
                             g_0_yyyyyy_0_xy_1,  \
                             g_0_yyyyyy_0_xz_0,  \
                             g_0_yyyyyy_0_xz_1,  \
                             g_0_yyyyyy_0_y_1,   \
                             g_0_yyyyyy_0_yy_0,  \
                             g_0_yyyyyy_0_yy_1,  \
                             g_0_yyyyyy_0_yz_0,  \
                             g_0_yyyyyy_0_yz_1,  \
                             g_0_yyyyyy_0_z_1,   \
                             g_0_yyyyyy_0_zz_0,  \
                             g_0_yyyyyy_0_zz_1,  \
                             g_0_yyyyyyy_0_xx_0, \
                             g_0_yyyyyyy_0_xy_0, \
                             g_0_yyyyyyy_0_xz_0, \
                             g_0_yyyyyyy_0_yy_0, \
                             g_0_yyyyyyy_0_yz_0, \
                             g_0_yyyyyyy_0_zz_0, \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyy_0_xx_0[i] =
            6.0 * g_0_yyyyy_0_xx_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xx_1[i] * fti_ab_0 + g_0_yyyyyy_0_xx_0[i] * pb_y + g_0_yyyyyy_0_xx_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xy_0[i] = 6.0 * g_0_yyyyy_0_xy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xy_1[i] * fti_ab_0 + g_0_yyyyyy_0_x_1[i] * fi_abcd_0 +
                                g_0_yyyyyy_0_xy_0[i] * pb_y + g_0_yyyyyy_0_xy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xz_0[i] =
            6.0 * g_0_yyyyy_0_xz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xz_0[i] * pb_y + g_0_yyyyyy_0_xz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yy_0[i] = 6.0 * g_0_yyyyy_0_yy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_y_1[i] * fi_abcd_0 +
                                g_0_yyyyyy_0_yy_0[i] * pb_y + g_0_yyyyyy_0_yy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yz_0[i] = 6.0 * g_0_yyyyy_0_yz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yz_1[i] * fti_ab_0 + g_0_yyyyyy_0_z_1[i] * fi_abcd_0 +
                                g_0_yyyyyy_0_yz_0[i] * pb_y + g_0_yyyyyy_0_yz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_zz_0[i] =
            6.0 * g_0_yyyyy_0_zz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_zz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zz_0[i] * pb_y + g_0_yyyyyy_0_zz_1[i] * wp_y[i];
    }

    /// Set up 174-180 components of targeted buffer : SKSD

    auto g_0_yyyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 174);

    auto g_0_yyyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 175);

    auto g_0_yyyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 176);

    auto g_0_yyyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 177);

    auto g_0_yyyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 178);

    auto g_0_yyyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 179);

#pragma omp simd aligned(g_0_yyyyyy_0_x_1,       \
                             g_0_yyyyyy_0_xx_0,  \
                             g_0_yyyyyy_0_xx_1,  \
                             g_0_yyyyyy_0_xy_0,  \
                             g_0_yyyyyy_0_xy_1,  \
                             g_0_yyyyyy_0_xz_0,  \
                             g_0_yyyyyy_0_xz_1,  \
                             g_0_yyyyyy_0_y_1,   \
                             g_0_yyyyyy_0_yy_0,  \
                             g_0_yyyyyy_0_yy_1,  \
                             g_0_yyyyyy_0_yz_0,  \
                             g_0_yyyyyy_0_yz_1,  \
                             g_0_yyyyyy_0_z_1,   \
                             g_0_yyyyyy_0_zz_0,  \
                             g_0_yyyyyy_0_zz_1,  \
                             g_0_yyyyyyz_0_xx_0, \
                             g_0_yyyyyyz_0_xy_0, \
                             g_0_yyyyyyz_0_xz_0, \
                             g_0_yyyyyyz_0_yy_0, \
                             g_0_yyyyyyz_0_yz_0, \
                             g_0_yyyyyyz_0_zz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyz_0_xx_0[i] = g_0_yyyyyy_0_xx_0[i] * pb_z + g_0_yyyyyy_0_xx_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xy_0[i] = g_0_yyyyyy_0_xy_0[i] * pb_z + g_0_yyyyyy_0_xy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xz_0[i] = g_0_yyyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xz_0[i] * pb_z + g_0_yyyyyy_0_xz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yy_0[i] = g_0_yyyyyy_0_yy_0[i] * pb_z + g_0_yyyyyy_0_yy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yz_0[i] = g_0_yyyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yz_0[i] * pb_z + g_0_yyyyyy_0_yz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_zz_0[i] = 2.0 * g_0_yyyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyyy_0_zz_0[i] * pb_z + g_0_yyyyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 180-186 components of targeted buffer : SKSD

    auto g_0_yyyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 180);

    auto g_0_yyyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 181);

    auto g_0_yyyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 182);

    auto g_0_yyyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 183);

    auto g_0_yyyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 184);

    auto g_0_yyyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 185);

#pragma omp simd aligned(g_0_yyyyy_0_xy_0,       \
                             g_0_yyyyy_0_xy_1,   \
                             g_0_yyyyy_0_yy_0,   \
                             g_0_yyyyy_0_yy_1,   \
                             g_0_yyyyyz_0_xy_0,  \
                             g_0_yyyyyz_0_xy_1,  \
                             g_0_yyyyyz_0_yy_0,  \
                             g_0_yyyyyz_0_yy_1,  \
                             g_0_yyyyyzz_0_xx_0, \
                             g_0_yyyyyzz_0_xy_0, \
                             g_0_yyyyyzz_0_xz_0, \
                             g_0_yyyyyzz_0_yy_0, \
                             g_0_yyyyyzz_0_yz_0, \
                             g_0_yyyyyzz_0_zz_0, \
                             g_0_yyyyzz_0_xx_0,  \
                             g_0_yyyyzz_0_xx_1,  \
                             g_0_yyyyzz_0_xz_0,  \
                             g_0_yyyyzz_0_xz_1,  \
                             g_0_yyyyzz_0_yz_0,  \
                             g_0_yyyyzz_0_yz_1,  \
                             g_0_yyyyzz_0_z_1,   \
                             g_0_yyyyzz_0_zz_0,  \
                             g_0_yyyyzz_0_zz_1,  \
                             g_0_yyyzz_0_xx_0,   \
                             g_0_yyyzz_0_xx_1,   \
                             g_0_yyyzz_0_xz_0,   \
                             g_0_yyyzz_0_xz_1,   \
                             g_0_yyyzz_0_yz_0,   \
                             g_0_yyyzz_0_yz_1,   \
                             g_0_yyyzz_0_zz_0,   \
                             g_0_yyyzz_0_zz_1,   \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzz_0_xx_0[i] =
            4.0 * g_0_yyyzz_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xx_1[i] * fti_ab_0 + g_0_yyyyzz_0_xx_0[i] * pb_y + g_0_yyyyzz_0_xx_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xy_0[i] =
            g_0_yyyyy_0_xy_0[i] * fi_ab_0 - g_0_yyyyy_0_xy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xy_0[i] * pb_z + g_0_yyyyyz_0_xy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xz_0[i] =
            4.0 * g_0_yyyzz_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xz_0[i] * pb_y + g_0_yyyyzz_0_xz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yy_0[i] =
            g_0_yyyyy_0_yy_0[i] * fi_ab_0 - g_0_yyyyy_0_yy_1[i] * fti_ab_0 + g_0_yyyyyz_0_yy_0[i] * pb_z + g_0_yyyyyz_0_yy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_yz_0[i] = 4.0 * g_0_yyyzz_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yz_1[i] * fti_ab_0 + g_0_yyyyzz_0_z_1[i] * fi_abcd_0 +
                                g_0_yyyyzz_0_yz_0[i] * pb_y + g_0_yyyyzz_0_yz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_zz_0[i] =
            4.0 * g_0_yyyzz_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_zz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zz_0[i] * pb_y + g_0_yyyyzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 186-192 components of targeted buffer : SKSD

    auto g_0_yyyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 186);

    auto g_0_yyyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 187);

    auto g_0_yyyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 188);

    auto g_0_yyyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 189);

    auto g_0_yyyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 190);

    auto g_0_yyyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 191);

#pragma omp simd aligned(g_0_yyyyz_0_xy_0,       \
                             g_0_yyyyz_0_xy_1,   \
                             g_0_yyyyz_0_yy_0,   \
                             g_0_yyyyz_0_yy_1,   \
                             g_0_yyyyzz_0_xy_0,  \
                             g_0_yyyyzz_0_xy_1,  \
                             g_0_yyyyzz_0_yy_0,  \
                             g_0_yyyyzz_0_yy_1,  \
                             g_0_yyyyzzz_0_xx_0, \
                             g_0_yyyyzzz_0_xy_0, \
                             g_0_yyyyzzz_0_xz_0, \
                             g_0_yyyyzzz_0_yy_0, \
                             g_0_yyyyzzz_0_yz_0, \
                             g_0_yyyyzzz_0_zz_0, \
                             g_0_yyyzzz_0_xx_0,  \
                             g_0_yyyzzz_0_xx_1,  \
                             g_0_yyyzzz_0_xz_0,  \
                             g_0_yyyzzz_0_xz_1,  \
                             g_0_yyyzzz_0_yz_0,  \
                             g_0_yyyzzz_0_yz_1,  \
                             g_0_yyyzzz_0_z_1,   \
                             g_0_yyyzzz_0_zz_0,  \
                             g_0_yyyzzz_0_zz_1,  \
                             g_0_yyzzz_0_xx_0,   \
                             g_0_yyzzz_0_xx_1,   \
                             g_0_yyzzz_0_xz_0,   \
                             g_0_yyzzz_0_xz_1,   \
                             g_0_yyzzz_0_yz_0,   \
                             g_0_yyzzz_0_yz_1,   \
                             g_0_yyzzz_0_zz_0,   \
                             g_0_yyzzz_0_zz_1,   \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzz_0_xx_0[i] =
            3.0 * g_0_yyzzz_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xx_1[i] * fti_ab_0 + g_0_yyyzzz_0_xx_0[i] * pb_y + g_0_yyyzzz_0_xx_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xy_0[i] =
            2.0 * g_0_yyyyz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xy_0[i] * pb_z + g_0_yyyyzz_0_xy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xz_0[i] =
            3.0 * g_0_yyzzz_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xz_0[i] * pb_y + g_0_yyyzzz_0_xz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yy_0[i] =
            2.0 * g_0_yyyyz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_yy_1[i] * fti_ab_0 + g_0_yyyyzz_0_yy_0[i] * pb_z + g_0_yyyyzz_0_yy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_yz_0[i] = 3.0 * g_0_yyzzz_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yz_1[i] * fti_ab_0 + g_0_yyyzzz_0_z_1[i] * fi_abcd_0 +
                                g_0_yyyzzz_0_yz_0[i] * pb_y + g_0_yyyzzz_0_yz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_zz_0[i] =
            3.0 * g_0_yyzzz_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_zz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zz_0[i] * pb_y + g_0_yyyzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 192-198 components of targeted buffer : SKSD

    auto g_0_yyyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 192);

    auto g_0_yyyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 193);

    auto g_0_yyyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 194);

    auto g_0_yyyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 195);

    auto g_0_yyyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 196);

    auto g_0_yyyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 197);

#pragma omp simd aligned(g_0_yyyzz_0_xy_0,       \
                             g_0_yyyzz_0_xy_1,   \
                             g_0_yyyzz_0_yy_0,   \
                             g_0_yyyzz_0_yy_1,   \
                             g_0_yyyzzz_0_xy_0,  \
                             g_0_yyyzzz_0_xy_1,  \
                             g_0_yyyzzz_0_yy_0,  \
                             g_0_yyyzzz_0_yy_1,  \
                             g_0_yyyzzzz_0_xx_0, \
                             g_0_yyyzzzz_0_xy_0, \
                             g_0_yyyzzzz_0_xz_0, \
                             g_0_yyyzzzz_0_yy_0, \
                             g_0_yyyzzzz_0_yz_0, \
                             g_0_yyyzzzz_0_zz_0, \
                             g_0_yyzzzz_0_xx_0,  \
                             g_0_yyzzzz_0_xx_1,  \
                             g_0_yyzzzz_0_xz_0,  \
                             g_0_yyzzzz_0_xz_1,  \
                             g_0_yyzzzz_0_yz_0,  \
                             g_0_yyzzzz_0_yz_1,  \
                             g_0_yyzzzz_0_z_1,   \
                             g_0_yyzzzz_0_zz_0,  \
                             g_0_yyzzzz_0_zz_1,  \
                             g_0_yzzzz_0_xx_0,   \
                             g_0_yzzzz_0_xx_1,   \
                             g_0_yzzzz_0_xz_0,   \
                             g_0_yzzzz_0_xz_1,   \
                             g_0_yzzzz_0_yz_0,   \
                             g_0_yzzzz_0_yz_1,   \
                             g_0_yzzzz_0_zz_0,   \
                             g_0_yzzzz_0_zz_1,   \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzz_0_xx_0[i] =
            2.0 * g_0_yzzzz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xx_1[i] * fti_ab_0 + g_0_yyzzzz_0_xx_0[i] * pb_y + g_0_yyzzzz_0_xx_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xy_0[i] =
            3.0 * g_0_yyyzz_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xy_0[i] * pb_z + g_0_yyyzzz_0_xy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xz_0[i] =
            2.0 * g_0_yzzzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xz_0[i] * pb_y + g_0_yyzzzz_0_xz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yy_0[i] =
            3.0 * g_0_yyyzz_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_yy_1[i] * fti_ab_0 + g_0_yyyzzz_0_yy_0[i] * pb_z + g_0_yyyzzz_0_yy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_yz_0[i] = 2.0 * g_0_yzzzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yz_1[i] * fti_ab_0 + g_0_yyzzzz_0_z_1[i] * fi_abcd_0 +
                                g_0_yyzzzz_0_yz_0[i] * pb_y + g_0_yyzzzz_0_yz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_zz_0[i] =
            2.0 * g_0_yzzzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_zz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zz_0[i] * pb_y + g_0_yyzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 198-204 components of targeted buffer : SKSD

    auto g_0_yyzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 198);

    auto g_0_yyzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 199);

    auto g_0_yyzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 200);

    auto g_0_yyzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 201);

    auto g_0_yyzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 202);

    auto g_0_yyzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 203);

#pragma omp simd aligned(g_0_yyzzz_0_xy_0,       \
                             g_0_yyzzz_0_xy_1,   \
                             g_0_yyzzz_0_yy_0,   \
                             g_0_yyzzz_0_yy_1,   \
                             g_0_yyzzzz_0_xy_0,  \
                             g_0_yyzzzz_0_xy_1,  \
                             g_0_yyzzzz_0_yy_0,  \
                             g_0_yyzzzz_0_yy_1,  \
                             g_0_yyzzzzz_0_xx_0, \
                             g_0_yyzzzzz_0_xy_0, \
                             g_0_yyzzzzz_0_xz_0, \
                             g_0_yyzzzzz_0_yy_0, \
                             g_0_yyzzzzz_0_yz_0, \
                             g_0_yyzzzzz_0_zz_0, \
                             g_0_yzzzzz_0_xx_0,  \
                             g_0_yzzzzz_0_xx_1,  \
                             g_0_yzzzzz_0_xz_0,  \
                             g_0_yzzzzz_0_xz_1,  \
                             g_0_yzzzzz_0_yz_0,  \
                             g_0_yzzzzz_0_yz_1,  \
                             g_0_yzzzzz_0_z_1,   \
                             g_0_yzzzzz_0_zz_0,  \
                             g_0_yzzzzz_0_zz_1,  \
                             g_0_zzzzz_0_xx_0,   \
                             g_0_zzzzz_0_xx_1,   \
                             g_0_zzzzz_0_xz_0,   \
                             g_0_zzzzz_0_xz_1,   \
                             g_0_zzzzz_0_yz_0,   \
                             g_0_zzzzz_0_yz_1,   \
                             g_0_zzzzz_0_zz_0,   \
                             g_0_zzzzz_0_zz_1,   \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzz_0_xx_0[i] =
            g_0_zzzzz_0_xx_0[i] * fi_ab_0 - g_0_zzzzz_0_xx_1[i] * fti_ab_0 + g_0_yzzzzz_0_xx_0[i] * pb_y + g_0_yzzzzz_0_xx_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xy_0[i] =
            4.0 * g_0_yyzzz_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xy_0[i] * pb_z + g_0_yyzzzz_0_xy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xz_0[i] =
            g_0_zzzzz_0_xz_0[i] * fi_ab_0 - g_0_zzzzz_0_xz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xz_0[i] * pb_y + g_0_yzzzzz_0_xz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yy_0[i] =
            4.0 * g_0_yyzzz_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_yy_1[i] * fti_ab_0 + g_0_yyzzzz_0_yy_0[i] * pb_z + g_0_yyzzzz_0_yy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_yz_0[i] = g_0_zzzzz_0_yz_0[i] * fi_ab_0 - g_0_zzzzz_0_yz_1[i] * fti_ab_0 + g_0_yzzzzz_0_z_1[i] * fi_abcd_0 +
                                g_0_yzzzzz_0_yz_0[i] * pb_y + g_0_yzzzzz_0_yz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_zz_0[i] =
            g_0_zzzzz_0_zz_0[i] * fi_ab_0 - g_0_zzzzz_0_zz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zz_0[i] * pb_y + g_0_yzzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 204-210 components of targeted buffer : SKSD

    auto g_0_yzzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 204);

    auto g_0_yzzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 205);

    auto g_0_yzzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 206);

    auto g_0_yzzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 207);

    auto g_0_yzzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 208);

    auto g_0_yzzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 209);

#pragma omp simd aligned(g_0_yzzzzzz_0_xx_0,     \
                             g_0_yzzzzzz_0_xy_0, \
                             g_0_yzzzzzz_0_xz_0, \
                             g_0_yzzzzzz_0_yy_0, \
                             g_0_yzzzzzz_0_yz_0, \
                             g_0_yzzzzzz_0_zz_0, \
                             g_0_zzzzzz_0_x_1,   \
                             g_0_zzzzzz_0_xx_0,  \
                             g_0_zzzzzz_0_xx_1,  \
                             g_0_zzzzzz_0_xy_0,  \
                             g_0_zzzzzz_0_xy_1,  \
                             g_0_zzzzzz_0_xz_0,  \
                             g_0_zzzzzz_0_xz_1,  \
                             g_0_zzzzzz_0_y_1,   \
                             g_0_zzzzzz_0_yy_0,  \
                             g_0_zzzzzz_0_yy_1,  \
                             g_0_zzzzzz_0_yz_0,  \
                             g_0_zzzzzz_0_yz_1,  \
                             g_0_zzzzzz_0_z_1,   \
                             g_0_zzzzzz_0_zz_0,  \
                             g_0_zzzzzz_0_zz_1,  \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzz_0_xx_0[i] = g_0_zzzzzz_0_xx_0[i] * pb_y + g_0_zzzzzz_0_xx_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xy_0[i] = g_0_zzzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xy_0[i] * pb_y + g_0_zzzzzz_0_xy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xz_0[i] = g_0_zzzzzz_0_xz_0[i] * pb_y + g_0_zzzzzz_0_xz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yy_0[i] = 2.0 * g_0_zzzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yy_0[i] * pb_y + g_0_zzzzzz_0_yy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yz_0[i] = g_0_zzzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yz_0[i] * pb_y + g_0_zzzzzz_0_yz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_zz_0[i] = g_0_zzzzzz_0_zz_0[i] * pb_y + g_0_zzzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 210-216 components of targeted buffer : SKSD

    auto g_0_zzzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sksd + 210);

    auto g_0_zzzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sksd + 211);

    auto g_0_zzzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sksd + 212);

    auto g_0_zzzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sksd + 213);

    auto g_0_zzzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sksd + 214);

    auto g_0_zzzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sksd + 215);

#pragma omp simd aligned(g_0_zzzzz_0_xx_0,       \
                             g_0_zzzzz_0_xx_1,   \
                             g_0_zzzzz_0_xy_0,   \
                             g_0_zzzzz_0_xy_1,   \
                             g_0_zzzzz_0_xz_0,   \
                             g_0_zzzzz_0_xz_1,   \
                             g_0_zzzzz_0_yy_0,   \
                             g_0_zzzzz_0_yy_1,   \
                             g_0_zzzzz_0_yz_0,   \
                             g_0_zzzzz_0_yz_1,   \
                             g_0_zzzzz_0_zz_0,   \
                             g_0_zzzzz_0_zz_1,   \
                             g_0_zzzzzz_0_x_1,   \
                             g_0_zzzzzz_0_xx_0,  \
                             g_0_zzzzzz_0_xx_1,  \
                             g_0_zzzzzz_0_xy_0,  \
                             g_0_zzzzzz_0_xy_1,  \
                             g_0_zzzzzz_0_xz_0,  \
                             g_0_zzzzzz_0_xz_1,  \
                             g_0_zzzzzz_0_y_1,   \
                             g_0_zzzzzz_0_yy_0,  \
                             g_0_zzzzzz_0_yy_1,  \
                             g_0_zzzzzz_0_yz_0,  \
                             g_0_zzzzzz_0_yz_1,  \
                             g_0_zzzzzz_0_z_1,   \
                             g_0_zzzzzz_0_zz_0,  \
                             g_0_zzzzzz_0_zz_1,  \
                             g_0_zzzzzzz_0_xx_0, \
                             g_0_zzzzzzz_0_xy_0, \
                             g_0_zzzzzzz_0_xz_0, \
                             g_0_zzzzzzz_0_yy_0, \
                             g_0_zzzzzzz_0_yz_0, \
                             g_0_zzzzzzz_0_zz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzz_0_xx_0[i] =
            6.0 * g_0_zzzzz_0_xx_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xx_1[i] * fti_ab_0 + g_0_zzzzzz_0_xx_0[i] * pb_z + g_0_zzzzzz_0_xx_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xy_0[i] =
            6.0 * g_0_zzzzz_0_xy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xy_0[i] * pb_z + g_0_zzzzzz_0_xy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xz_0[i] = 6.0 * g_0_zzzzz_0_xz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xz_1[i] * fti_ab_0 + g_0_zzzzzz_0_x_1[i] * fi_abcd_0 +
                                g_0_zzzzzz_0_xz_0[i] * pb_z + g_0_zzzzzz_0_xz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yy_0[i] =
            6.0 * g_0_zzzzz_0_yy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yy_1[i] * fti_ab_0 + g_0_zzzzzz_0_yy_0[i] * pb_z + g_0_zzzzzz_0_yy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yz_0[i] = 6.0 * g_0_zzzzz_0_yz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yz_1[i] * fti_ab_0 + g_0_zzzzzz_0_y_1[i] * fi_abcd_0 +
                                g_0_zzzzzz_0_yz_0[i] * pb_z + g_0_zzzzzz_0_yz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_zz_0[i] = 6.0 * g_0_zzzzz_0_zz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_zz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_z_1[i] * fi_abcd_0 +
                                g_0_zzzzzz_0_zz_0[i] * pb_z + g_0_zzzzzz_0_zz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
