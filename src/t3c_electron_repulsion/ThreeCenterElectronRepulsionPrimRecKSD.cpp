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

#include "ThreeCenterElectronRepulsionPrimRecKSD.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksd,
                                 size_t idx_eri_0_hsd,
                                 size_t idx_eri_1_hsd,
                                 size_t idx_eri_1_isp,
                                 size_t idx_eri_1_isd,
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

    /// Set up components of auxilary buffer : HSD

    auto g_xxxxx_0_xx_0 = pbuffer.data(idx_eri_0_hsd);

    auto g_xxxxx_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 1);

    auto g_xxxxx_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 2);

    auto g_xxxxx_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 3);

    auto g_xxxxx_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 4);

    auto g_xxxxx_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 5);

    auto g_xxxxy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 6);

    auto g_xxxxy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 8);

    auto g_xxxxz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 12);

    auto g_xxxxz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 13);

    auto g_xxxyy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 18);

    auto g_xxxyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 19);

    auto g_xxxyy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 20);

    auto g_xxxyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 21);

    auto g_xxxyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 22);

    auto g_xxxyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 23);

    auto g_xxxzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 30);

    auto g_xxxzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 31);

    auto g_xxxzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 32);

    auto g_xxxzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 33);

    auto g_xxxzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 34);

    auto g_xxxzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 35);

    auto g_xxyyy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 36);

    auto g_xxyyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 37);

    auto g_xxyyy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 38);

    auto g_xxyyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 39);

    auto g_xxyyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 40);

    auto g_xxyyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 41);

    auto g_xxyyz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 43);

    auto g_xxyzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 48);

    auto g_xxyzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 50);

    auto g_xxzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 54);

    auto g_xxzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 55);

    auto g_xxzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 56);

    auto g_xxzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 57);

    auto g_xxzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 58);

    auto g_xxzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 59);

    auto g_xyyyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 61);

    auto g_xyyyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 63);

    auto g_xyyyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 64);

    auto g_xyyyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 65);

    auto g_xyyzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 75);

    auto g_xyyzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 76);

    auto g_xyyzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 77);

    auto g_xzzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 86);

    auto g_xzzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 87);

    auto g_xzzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 88);

    auto g_xzzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 89);

    auto g_yyyyy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 90);

    auto g_yyyyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 91);

    auto g_yyyyy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 92);

    auto g_yyyyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 93);

    auto g_yyyyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 94);

    auto g_yyyyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 95);

    auto g_yyyyz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 97);

    auto g_yyyyz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 99);

    auto g_yyyzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 102);

    auto g_yyyzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 103);

    auto g_yyyzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 104);

    auto g_yyyzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 105);

    auto g_yyyzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 106);

    auto g_yyyzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 107);

    auto g_yyzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 108);

    auto g_yyzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 109);

    auto g_yyzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 110);

    auto g_yyzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 111);

    auto g_yyzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 112);

    auto g_yyzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 113);

    auto g_yzzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 114);

    auto g_yzzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 116);

    auto g_yzzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 118);

    auto g_yzzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 119);

    auto g_zzzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 120);

    auto g_zzzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 121);

    auto g_zzzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 122);

    auto g_zzzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 123);

    auto g_zzzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 124);

    auto g_zzzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 125);

    /// Set up components of auxilary buffer : HSD

    auto g_xxxxx_0_xx_1 = pbuffer.data(idx_eri_1_hsd);

    auto g_xxxxx_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 1);

    auto g_xxxxx_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 2);

    auto g_xxxxx_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 3);

    auto g_xxxxx_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 4);

    auto g_xxxxx_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 5);

    auto g_xxxxy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 6);

    auto g_xxxxy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 8);

    auto g_xxxxz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 12);

    auto g_xxxxz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 13);

    auto g_xxxyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 18);

    auto g_xxxyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 19);

    auto g_xxxyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 20);

    auto g_xxxyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 21);

    auto g_xxxyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 22);

    auto g_xxxyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 23);

    auto g_xxxzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 30);

    auto g_xxxzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 31);

    auto g_xxxzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 32);

    auto g_xxxzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 33);

    auto g_xxxzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 34);

    auto g_xxxzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 35);

    auto g_xxyyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 36);

    auto g_xxyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 37);

    auto g_xxyyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 38);

    auto g_xxyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 39);

    auto g_xxyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 40);

    auto g_xxyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 41);

    auto g_xxyyz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 43);

    auto g_xxyzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 48);

    auto g_xxyzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 50);

    auto g_xxzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 54);

    auto g_xxzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 55);

    auto g_xxzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 56);

    auto g_xxzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 57);

    auto g_xxzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 58);

    auto g_xxzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 59);

    auto g_xyyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 61);

    auto g_xyyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 63);

    auto g_xyyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 64);

    auto g_xyyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 65);

    auto g_xyyzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 75);

    auto g_xyyzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 76);

    auto g_xyyzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 77);

    auto g_xzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 86);

    auto g_xzzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 87);

    auto g_xzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 88);

    auto g_xzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 89);

    auto g_yyyyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 90);

    auto g_yyyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 91);

    auto g_yyyyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 92);

    auto g_yyyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 93);

    auto g_yyyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 94);

    auto g_yyyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 95);

    auto g_yyyyz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 97);

    auto g_yyyyz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 99);

    auto g_yyyzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 102);

    auto g_yyyzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 103);

    auto g_yyyzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 104);

    auto g_yyyzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 105);

    auto g_yyyzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 106);

    auto g_yyyzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 107);

    auto g_yyzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 108);

    auto g_yyzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 109);

    auto g_yyzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 110);

    auto g_yyzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 111);

    auto g_yyzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 112);

    auto g_yyzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 113);

    auto g_yzzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 114);

    auto g_yzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 116);

    auto g_yzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 118);

    auto g_yzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 119);

    auto g_zzzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 120);

    auto g_zzzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 121);

    auto g_zzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 122);

    auto g_zzzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 123);

    auto g_zzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 124);

    auto g_zzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 125);

    /// Set up components of auxilary buffer : ISP

    auto g_xxxxxx_0_x_1 = pbuffer.data(idx_eri_1_isp);

    auto g_xxxxxx_0_y_1 = pbuffer.data(idx_eri_1_isp + 1);

    auto g_xxxxxx_0_z_1 = pbuffer.data(idx_eri_1_isp + 2);

    auto g_xxxxxz_0_z_1 = pbuffer.data(idx_eri_1_isp + 8);

    auto g_xxxxyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 9);

    auto g_xxxxyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 10);

    auto g_xxxxyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 11);

    auto g_xxxxzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 15);

    auto g_xxxxzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 16);

    auto g_xxxxzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 17);

    auto g_xxxyyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 18);

    auto g_xxxyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 19);

    auto g_xxxyyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 20);

    auto g_xxxzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 27);

    auto g_xxxzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 28);

    auto g_xxxzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 29);

    auto g_xxyyyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 30);

    auto g_xxyyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 31);

    auto g_xxyyyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 32);

    auto g_xxzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 42);

    auto g_xxzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 43);

    auto g_xxzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 44);

    auto g_xyyyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 46);

    auto g_xzzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 62);

    auto g_yyyyyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 63);

    auto g_yyyyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 64);

    auto g_yyyyyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 65);

    auto g_yyyyyz_0_z_1 = pbuffer.data(idx_eri_1_isp + 68);

    auto g_yyyyzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 69);

    auto g_yyyyzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 70);

    auto g_yyyyzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 71);

    auto g_yyyzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 72);

    auto g_yyyzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 73);

    auto g_yyyzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 74);

    auto g_yyzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 75);

    auto g_yyzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 76);

    auto g_yyzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 77);

    auto g_yzzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 79);

    auto g_yzzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 80);

    auto g_zzzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 81);

    auto g_zzzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 82);

    auto g_zzzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 83);

    /// Set up components of auxilary buffer : ISD

    auto g_xxxxxx_0_xx_1 = pbuffer.data(idx_eri_1_isd);

    auto g_xxxxxx_0_xy_1 = pbuffer.data(idx_eri_1_isd + 1);

    auto g_xxxxxx_0_xz_1 = pbuffer.data(idx_eri_1_isd + 2);

    auto g_xxxxxx_0_yy_1 = pbuffer.data(idx_eri_1_isd + 3);

    auto g_xxxxxx_0_yz_1 = pbuffer.data(idx_eri_1_isd + 4);

    auto g_xxxxxx_0_zz_1 = pbuffer.data(idx_eri_1_isd + 5);

    auto g_xxxxxy_0_xx_1 = pbuffer.data(idx_eri_1_isd + 6);

    auto g_xxxxxy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 7);

    auto g_xxxxxy_0_xz_1 = pbuffer.data(idx_eri_1_isd + 8);

    auto g_xxxxxy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 9);

    auto g_xxxxxz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 12);

    auto g_xxxxxz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 13);

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

    auto g_xxxyyz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 43);

    auto g_xxxyzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 48);

    auto g_xxxyzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 50);

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

    auto g_xxyyyz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 67);

    auto g_xxyyzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 72);

    auto g_xxyyzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 73);

    auto g_xxyyzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 74);

    auto g_xxyyzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 75);

    auto g_xxyyzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 76);

    auto g_xxyyzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 77);

    auto g_xxyzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 78);

    auto g_xxyzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 80);

    auto g_xxzzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 84);

    auto g_xxzzzz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 85);

    auto g_xxzzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 86);

    auto g_xxzzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 87);

    auto g_xxzzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 88);

    auto g_xxzzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 89);

    auto g_xyyyyy_0_xx_1 = pbuffer.data(idx_eri_1_isd + 90);

    auto g_xyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 91);

    auto g_xyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 93);

    auto g_xyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_isd + 94);

    auto g_xyyyyy_0_zz_1 = pbuffer.data(idx_eri_1_isd + 95);

    auto g_xyyyzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 105);

    auto g_xyyyzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 106);

    auto g_xyyyzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 107);

    auto g_xyyzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 111);

    auto g_xyyzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 112);

    auto g_xyyzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 113);

    auto g_xzzzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 120);

    auto g_xzzzzz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 122);

    auto g_xzzzzz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 123);

    auto g_xzzzzz_0_yz_1 = pbuffer.data(idx_eri_1_isd + 124);

    auto g_xzzzzz_0_zz_1 = pbuffer.data(idx_eri_1_isd + 125);

    auto g_yyyyyy_0_xx_1 = pbuffer.data(idx_eri_1_isd + 126);

    auto g_yyyyyy_0_xy_1 = pbuffer.data(idx_eri_1_isd + 127);

    auto g_yyyyyy_0_xz_1 = pbuffer.data(idx_eri_1_isd + 128);

    auto g_yyyyyy_0_yy_1 = pbuffer.data(idx_eri_1_isd + 129);

    auto g_yyyyyy_0_yz_1 = pbuffer.data(idx_eri_1_isd + 130);

    auto g_yyyyyy_0_zz_1 = pbuffer.data(idx_eri_1_isd + 131);

    auto g_yyyyyz_0_xy_1 = pbuffer.data(idx_eri_1_isd + 133);

    auto g_yyyyyz_0_xz_1 = pbuffer.data(idx_eri_1_isd + 134);

    auto g_yyyyyz_0_yy_1 = pbuffer.data(idx_eri_1_isd + 135);

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

    auto g_yzzzzz_0_xx_1 = pbuffer.data(idx_eri_1_isd + 156);

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

    /// Set up 0-6 components of targeted buffer : KSD

    auto g_xxxxxxx_0_xx_0 = pbuffer.data(idx_eri_0_ksd);

    auto g_xxxxxxx_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 1);

    auto g_xxxxxxx_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 2);

    auto g_xxxxxxx_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 3);

    auto g_xxxxxxx_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 4);

    auto g_xxxxxxx_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 5);

    #pragma omp simd aligned(g_xxxxx_0_xx_0, g_xxxxx_0_xx_1, g_xxxxx_0_xy_0, g_xxxxx_0_xy_1, g_xxxxx_0_xz_0, g_xxxxx_0_xz_1, g_xxxxx_0_yy_0, g_xxxxx_0_yy_1, g_xxxxx_0_yz_0, g_xxxxx_0_yz_1, g_xxxxx_0_zz_0, g_xxxxx_0_zz_1, g_xxxxxx_0_x_1, g_xxxxxx_0_xx_1, g_xxxxxx_0_xy_1, g_xxxxxx_0_xz_1, g_xxxxxx_0_y_1, g_xxxxxx_0_yy_1, g_xxxxxx_0_yz_1, g_xxxxxx_0_z_1, g_xxxxxx_0_zz_1, g_xxxxxxx_0_xx_0, g_xxxxxxx_0_xy_0, g_xxxxxxx_0_xz_0, g_xxxxxxx_0_yy_0, g_xxxxxxx_0_yz_0, g_xxxxxxx_0_zz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_xx_0[i] = 6.0 * g_xxxxx_0_xx_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xx_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_x_1[i] * fi_acd_0 + g_xxxxxx_0_xx_1[i] * wa_x[i];

        g_xxxxxxx_0_xy_0[i] = 6.0 * g_xxxxx_0_xy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xy_1[i] * fz_be_0 + g_xxxxxx_0_y_1[i] * fi_acd_0 + g_xxxxxx_0_xy_1[i] * wa_x[i];

        g_xxxxxxx_0_xz_0[i] = 6.0 * g_xxxxx_0_xz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xz_1[i] * fz_be_0 + g_xxxxxx_0_z_1[i] * fi_acd_0 + g_xxxxxx_0_xz_1[i] * wa_x[i];

        g_xxxxxxx_0_yy_0[i] = 6.0 * g_xxxxx_0_yy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yy_1[i] * fz_be_0 + g_xxxxxx_0_yy_1[i] * wa_x[i];

        g_xxxxxxx_0_yz_0[i] = 6.0 * g_xxxxx_0_yz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yz_1[i] * fz_be_0 + g_xxxxxx_0_yz_1[i] * wa_x[i];

        g_xxxxxxx_0_zz_0[i] = 6.0 * g_xxxxx_0_zz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_zz_1[i] * fz_be_0 + g_xxxxxx_0_zz_1[i] * wa_x[i];
    }

    /// Set up 6-12 components of targeted buffer : KSD

    auto g_xxxxxxy_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 6);

    auto g_xxxxxxy_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 7);

    auto g_xxxxxxy_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 8);

    auto g_xxxxxxy_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 9);

    auto g_xxxxxxy_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 10);

    auto g_xxxxxxy_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 11);

    #pragma omp simd aligned(g_xxxxxx_0_x_1, g_xxxxxx_0_xx_1, g_xxxxxx_0_xy_1, g_xxxxxx_0_xz_1, g_xxxxxx_0_y_1, g_xxxxxx_0_yy_1, g_xxxxxx_0_yz_1, g_xxxxxx_0_z_1, g_xxxxxx_0_zz_1, g_xxxxxxy_0_xx_0, g_xxxxxxy_0_xy_0, g_xxxxxxy_0_xz_0, g_xxxxxxy_0_yy_0, g_xxxxxxy_0_yz_0, g_xxxxxxy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_xx_0[i] = g_xxxxxx_0_xx_1[i] * wa_y[i];

        g_xxxxxxy_0_xy_0[i] = g_xxxxxx_0_x_1[i] * fi_acd_0 + g_xxxxxx_0_xy_1[i] * wa_y[i];

        g_xxxxxxy_0_xz_0[i] = g_xxxxxx_0_xz_1[i] * wa_y[i];

        g_xxxxxxy_0_yy_0[i] = 2.0 * g_xxxxxx_0_y_1[i] * fi_acd_0 + g_xxxxxx_0_yy_1[i] * wa_y[i];

        g_xxxxxxy_0_yz_0[i] = g_xxxxxx_0_z_1[i] * fi_acd_0 + g_xxxxxx_0_yz_1[i] * wa_y[i];

        g_xxxxxxy_0_zz_0[i] = g_xxxxxx_0_zz_1[i] * wa_y[i];
    }

    /// Set up 12-18 components of targeted buffer : KSD

    auto g_xxxxxxz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 12);

    auto g_xxxxxxz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 13);

    auto g_xxxxxxz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 14);

    auto g_xxxxxxz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 15);

    auto g_xxxxxxz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 16);

    auto g_xxxxxxz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 17);

    #pragma omp simd aligned(g_xxxxxx_0_x_1, g_xxxxxx_0_xx_1, g_xxxxxx_0_xy_1, g_xxxxxx_0_xz_1, g_xxxxxx_0_y_1, g_xxxxxx_0_yy_1, g_xxxxxx_0_yz_1, g_xxxxxx_0_z_1, g_xxxxxx_0_zz_1, g_xxxxxxz_0_xx_0, g_xxxxxxz_0_xy_0, g_xxxxxxz_0_xz_0, g_xxxxxxz_0_yy_0, g_xxxxxxz_0_yz_0, g_xxxxxxz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_xx_0[i] = g_xxxxxx_0_xx_1[i] * wa_z[i];

        g_xxxxxxz_0_xy_0[i] = g_xxxxxx_0_xy_1[i] * wa_z[i];

        g_xxxxxxz_0_xz_0[i] = g_xxxxxx_0_x_1[i] * fi_acd_0 + g_xxxxxx_0_xz_1[i] * wa_z[i];

        g_xxxxxxz_0_yy_0[i] = g_xxxxxx_0_yy_1[i] * wa_z[i];

        g_xxxxxxz_0_yz_0[i] = g_xxxxxx_0_y_1[i] * fi_acd_0 + g_xxxxxx_0_yz_1[i] * wa_z[i];

        g_xxxxxxz_0_zz_0[i] = 2.0 * g_xxxxxx_0_z_1[i] * fi_acd_0 + g_xxxxxx_0_zz_1[i] * wa_z[i];
    }

    /// Set up 18-24 components of targeted buffer : KSD

    auto g_xxxxxyy_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 18);

    auto g_xxxxxyy_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 19);

    auto g_xxxxxyy_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 20);

    auto g_xxxxxyy_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 21);

    auto g_xxxxxyy_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 22);

    auto g_xxxxxyy_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 23);

    #pragma omp simd aligned(g_xxxxx_0_xx_0, g_xxxxx_0_xx_1, g_xxxxx_0_xz_0, g_xxxxx_0_xz_1, g_xxxxxy_0_xx_1, g_xxxxxy_0_xz_1, g_xxxxxyy_0_xx_0, g_xxxxxyy_0_xy_0, g_xxxxxyy_0_xz_0, g_xxxxxyy_0_yy_0, g_xxxxxyy_0_yz_0, g_xxxxxyy_0_zz_0, g_xxxxyy_0_xy_1, g_xxxxyy_0_y_1, g_xxxxyy_0_yy_1, g_xxxxyy_0_yz_1, g_xxxxyy_0_zz_1, g_xxxyy_0_xy_0, g_xxxyy_0_xy_1, g_xxxyy_0_yy_0, g_xxxyy_0_yy_1, g_xxxyy_0_yz_0, g_xxxyy_0_yz_1, g_xxxyy_0_zz_0, g_xxxyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxyy_0_xx_0[i] = g_xxxxx_0_xx_0[i] * fbe_0 - g_xxxxx_0_xx_1[i] * fz_be_0 + g_xxxxxy_0_xx_1[i] * wa_y[i];

        g_xxxxxyy_0_xy_0[i] = 4.0 * g_xxxyy_0_xy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xy_1[i] * fz_be_0 + g_xxxxyy_0_y_1[i] * fi_acd_0 + g_xxxxyy_0_xy_1[i] * wa_x[i];

        g_xxxxxyy_0_xz_0[i] = g_xxxxx_0_xz_0[i] * fbe_0 - g_xxxxx_0_xz_1[i] * fz_be_0 + g_xxxxxy_0_xz_1[i] * wa_y[i];

        g_xxxxxyy_0_yy_0[i] = 4.0 * g_xxxyy_0_yy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yy_1[i] * fz_be_0 + g_xxxxyy_0_yy_1[i] * wa_x[i];

        g_xxxxxyy_0_yz_0[i] = 4.0 * g_xxxyy_0_yz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yz_1[i] * fz_be_0 + g_xxxxyy_0_yz_1[i] * wa_x[i];

        g_xxxxxyy_0_zz_0[i] = 4.0 * g_xxxyy_0_zz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_zz_1[i] * fz_be_0 + g_xxxxyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 24-30 components of targeted buffer : KSD

    auto g_xxxxxyz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 24);

    auto g_xxxxxyz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 25);

    auto g_xxxxxyz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 26);

    auto g_xxxxxyz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 27);

    auto g_xxxxxyz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 28);

    auto g_xxxxxyz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 29);

    #pragma omp simd aligned(g_xxxxxy_0_xy_1, g_xxxxxy_0_yy_1, g_xxxxxyz_0_xx_0, g_xxxxxyz_0_xy_0, g_xxxxxyz_0_xz_0, g_xxxxxyz_0_yy_0, g_xxxxxyz_0_yz_0, g_xxxxxyz_0_zz_0, g_xxxxxz_0_xx_1, g_xxxxxz_0_xz_1, g_xxxxxz_0_yz_1, g_xxxxxz_0_z_1, g_xxxxxz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyz_0_xx_0[i] = g_xxxxxz_0_xx_1[i] * wa_y[i];

        g_xxxxxyz_0_xy_0[i] = g_xxxxxy_0_xy_1[i] * wa_z[i];

        g_xxxxxyz_0_xz_0[i] = g_xxxxxz_0_xz_1[i] * wa_y[i];

        g_xxxxxyz_0_yy_0[i] = g_xxxxxy_0_yy_1[i] * wa_z[i];

        g_xxxxxyz_0_yz_0[i] = g_xxxxxz_0_z_1[i] * fi_acd_0 + g_xxxxxz_0_yz_1[i] * wa_y[i];

        g_xxxxxyz_0_zz_0[i] = g_xxxxxz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 30-36 components of targeted buffer : KSD

    auto g_xxxxxzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 30);

    auto g_xxxxxzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 31);

    auto g_xxxxxzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 32);

    auto g_xxxxxzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 33);

    auto g_xxxxxzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 34);

    auto g_xxxxxzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 35);

    #pragma omp simd aligned(g_xxxxx_0_xx_0, g_xxxxx_0_xx_1, g_xxxxx_0_xy_0, g_xxxxx_0_xy_1, g_xxxxxz_0_xx_1, g_xxxxxz_0_xy_1, g_xxxxxzz_0_xx_0, g_xxxxxzz_0_xy_0, g_xxxxxzz_0_xz_0, g_xxxxxzz_0_yy_0, g_xxxxxzz_0_yz_0, g_xxxxxzz_0_zz_0, g_xxxxzz_0_xz_1, g_xxxxzz_0_yy_1, g_xxxxzz_0_yz_1, g_xxxxzz_0_z_1, g_xxxxzz_0_zz_1, g_xxxzz_0_xz_0, g_xxxzz_0_xz_1, g_xxxzz_0_yy_0, g_xxxzz_0_yy_1, g_xxxzz_0_yz_0, g_xxxzz_0_yz_1, g_xxxzz_0_zz_0, g_xxxzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxzz_0_xx_0[i] = g_xxxxx_0_xx_0[i] * fbe_0 - g_xxxxx_0_xx_1[i] * fz_be_0 + g_xxxxxz_0_xx_1[i] * wa_z[i];

        g_xxxxxzz_0_xy_0[i] = g_xxxxx_0_xy_0[i] * fbe_0 - g_xxxxx_0_xy_1[i] * fz_be_0 + g_xxxxxz_0_xy_1[i] * wa_z[i];

        g_xxxxxzz_0_xz_0[i] = 4.0 * g_xxxzz_0_xz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xz_1[i] * fz_be_0 + g_xxxxzz_0_z_1[i] * fi_acd_0 + g_xxxxzz_0_xz_1[i] * wa_x[i];

        g_xxxxxzz_0_yy_0[i] = 4.0 * g_xxxzz_0_yy_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yy_1[i] * fz_be_0 + g_xxxxzz_0_yy_1[i] * wa_x[i];

        g_xxxxxzz_0_yz_0[i] = 4.0 * g_xxxzz_0_yz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yz_1[i] * fz_be_0 + g_xxxxzz_0_yz_1[i] * wa_x[i];

        g_xxxxxzz_0_zz_0[i] = 4.0 * g_xxxzz_0_zz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_zz_1[i] * fz_be_0 + g_xxxxzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 36-42 components of targeted buffer : KSD

    auto g_xxxxyyy_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 36);

    auto g_xxxxyyy_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 37);

    auto g_xxxxyyy_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 38);

    auto g_xxxxyyy_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 39);

    auto g_xxxxyyy_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 40);

    auto g_xxxxyyy_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 41);

    #pragma omp simd aligned(g_xxxxy_0_xx_0, g_xxxxy_0_xx_1, g_xxxxy_0_xz_0, g_xxxxy_0_xz_1, g_xxxxyy_0_xx_1, g_xxxxyy_0_xz_1, g_xxxxyyy_0_xx_0, g_xxxxyyy_0_xy_0, g_xxxxyyy_0_xz_0, g_xxxxyyy_0_yy_0, g_xxxxyyy_0_yz_0, g_xxxxyyy_0_zz_0, g_xxxyyy_0_xy_1, g_xxxyyy_0_y_1, g_xxxyyy_0_yy_1, g_xxxyyy_0_yz_1, g_xxxyyy_0_zz_1, g_xxyyy_0_xy_0, g_xxyyy_0_xy_1, g_xxyyy_0_yy_0, g_xxyyy_0_yy_1, g_xxyyy_0_yz_0, g_xxyyy_0_yz_1, g_xxyyy_0_zz_0, g_xxyyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyyy_0_xx_0[i] = 2.0 * g_xxxxy_0_xx_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xx_1[i] * fz_be_0 + g_xxxxyy_0_xx_1[i] * wa_y[i];

        g_xxxxyyy_0_xy_0[i] = 3.0 * g_xxyyy_0_xy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xy_1[i] * fz_be_0 + g_xxxyyy_0_y_1[i] * fi_acd_0 + g_xxxyyy_0_xy_1[i] * wa_x[i];

        g_xxxxyyy_0_xz_0[i] = 2.0 * g_xxxxy_0_xz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xz_1[i] * fz_be_0 + g_xxxxyy_0_xz_1[i] * wa_y[i];

        g_xxxxyyy_0_yy_0[i] = 3.0 * g_xxyyy_0_yy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yy_1[i] * fz_be_0 + g_xxxyyy_0_yy_1[i] * wa_x[i];

        g_xxxxyyy_0_yz_0[i] = 3.0 * g_xxyyy_0_yz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yz_1[i] * fz_be_0 + g_xxxyyy_0_yz_1[i] * wa_x[i];

        g_xxxxyyy_0_zz_0[i] = 3.0 * g_xxyyy_0_zz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_zz_1[i] * fz_be_0 + g_xxxyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 42-48 components of targeted buffer : KSD

    auto g_xxxxyyz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 42);

    auto g_xxxxyyz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 43);

    auto g_xxxxyyz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 44);

    auto g_xxxxyyz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 45);

    auto g_xxxxyyz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 46);

    auto g_xxxxyyz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 47);

    #pragma omp simd aligned(g_xxxxyy_0_x_1, g_xxxxyy_0_xx_1, g_xxxxyy_0_xy_1, g_xxxxyy_0_xz_1, g_xxxxyy_0_y_1, g_xxxxyy_0_yy_1, g_xxxxyy_0_yz_1, g_xxxxyy_0_z_1, g_xxxxyy_0_zz_1, g_xxxxyyz_0_xx_0, g_xxxxyyz_0_xy_0, g_xxxxyyz_0_xz_0, g_xxxxyyz_0_yy_0, g_xxxxyyz_0_yz_0, g_xxxxyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_xx_0[i] = g_xxxxyy_0_xx_1[i] * wa_z[i];

        g_xxxxyyz_0_xy_0[i] = g_xxxxyy_0_xy_1[i] * wa_z[i];

        g_xxxxyyz_0_xz_0[i] = g_xxxxyy_0_x_1[i] * fi_acd_0 + g_xxxxyy_0_xz_1[i] * wa_z[i];

        g_xxxxyyz_0_yy_0[i] = g_xxxxyy_0_yy_1[i] * wa_z[i];

        g_xxxxyyz_0_yz_0[i] = g_xxxxyy_0_y_1[i] * fi_acd_0 + g_xxxxyy_0_yz_1[i] * wa_z[i];

        g_xxxxyyz_0_zz_0[i] = 2.0 * g_xxxxyy_0_z_1[i] * fi_acd_0 + g_xxxxyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 48-54 components of targeted buffer : KSD

    auto g_xxxxyzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 48);

    auto g_xxxxyzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 49);

    auto g_xxxxyzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 50);

    auto g_xxxxyzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 51);

    auto g_xxxxyzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 52);

    auto g_xxxxyzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 53);

    #pragma omp simd aligned(g_xxxxyzz_0_xx_0, g_xxxxyzz_0_xy_0, g_xxxxyzz_0_xz_0, g_xxxxyzz_0_yy_0, g_xxxxyzz_0_yz_0, g_xxxxyzz_0_zz_0, g_xxxxzz_0_x_1, g_xxxxzz_0_xx_1, g_xxxxzz_0_xy_1, g_xxxxzz_0_xz_1, g_xxxxzz_0_y_1, g_xxxxzz_0_yy_1, g_xxxxzz_0_yz_1, g_xxxxzz_0_z_1, g_xxxxzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_xx_0[i] = g_xxxxzz_0_xx_1[i] * wa_y[i];

        g_xxxxyzz_0_xy_0[i] = g_xxxxzz_0_x_1[i] * fi_acd_0 + g_xxxxzz_0_xy_1[i] * wa_y[i];

        g_xxxxyzz_0_xz_0[i] = g_xxxxzz_0_xz_1[i] * wa_y[i];

        g_xxxxyzz_0_yy_0[i] = 2.0 * g_xxxxzz_0_y_1[i] * fi_acd_0 + g_xxxxzz_0_yy_1[i] * wa_y[i];

        g_xxxxyzz_0_yz_0[i] = g_xxxxzz_0_z_1[i] * fi_acd_0 + g_xxxxzz_0_yz_1[i] * wa_y[i];

        g_xxxxyzz_0_zz_0[i] = g_xxxxzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 54-60 components of targeted buffer : KSD

    auto g_xxxxzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 54);

    auto g_xxxxzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 55);

    auto g_xxxxzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 56);

    auto g_xxxxzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 57);

    auto g_xxxxzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 58);

    auto g_xxxxzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 59);

    #pragma omp simd aligned(g_xxxxz_0_xx_0, g_xxxxz_0_xx_1, g_xxxxz_0_xy_0, g_xxxxz_0_xy_1, g_xxxxzz_0_xx_1, g_xxxxzz_0_xy_1, g_xxxxzzz_0_xx_0, g_xxxxzzz_0_xy_0, g_xxxxzzz_0_xz_0, g_xxxxzzz_0_yy_0, g_xxxxzzz_0_yz_0, g_xxxxzzz_0_zz_0, g_xxxzzz_0_xz_1, g_xxxzzz_0_yy_1, g_xxxzzz_0_yz_1, g_xxxzzz_0_z_1, g_xxxzzz_0_zz_1, g_xxzzz_0_xz_0, g_xxzzz_0_xz_1, g_xxzzz_0_yy_0, g_xxzzz_0_yy_1, g_xxzzz_0_yz_0, g_xxzzz_0_yz_1, g_xxzzz_0_zz_0, g_xxzzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzzz_0_xx_0[i] = 2.0 * g_xxxxz_0_xx_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xx_1[i] * fz_be_0 + g_xxxxzz_0_xx_1[i] * wa_z[i];

        g_xxxxzzz_0_xy_0[i] = 2.0 * g_xxxxz_0_xy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xy_1[i] * fz_be_0 + g_xxxxzz_0_xy_1[i] * wa_z[i];

        g_xxxxzzz_0_xz_0[i] = 3.0 * g_xxzzz_0_xz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xz_1[i] * fz_be_0 + g_xxxzzz_0_z_1[i] * fi_acd_0 + g_xxxzzz_0_xz_1[i] * wa_x[i];

        g_xxxxzzz_0_yy_0[i] = 3.0 * g_xxzzz_0_yy_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yy_1[i] * fz_be_0 + g_xxxzzz_0_yy_1[i] * wa_x[i];

        g_xxxxzzz_0_yz_0[i] = 3.0 * g_xxzzz_0_yz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yz_1[i] * fz_be_0 + g_xxxzzz_0_yz_1[i] * wa_x[i];

        g_xxxxzzz_0_zz_0[i] = 3.0 * g_xxzzz_0_zz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_zz_1[i] * fz_be_0 + g_xxxzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 60-66 components of targeted buffer : KSD

    auto g_xxxyyyy_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 60);

    auto g_xxxyyyy_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 61);

    auto g_xxxyyyy_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 62);

    auto g_xxxyyyy_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 63);

    auto g_xxxyyyy_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 64);

    auto g_xxxyyyy_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 65);

    #pragma omp simd aligned(g_xxxyy_0_xx_0, g_xxxyy_0_xx_1, g_xxxyy_0_xz_0, g_xxxyy_0_xz_1, g_xxxyyy_0_xx_1, g_xxxyyy_0_xz_1, g_xxxyyyy_0_xx_0, g_xxxyyyy_0_xy_0, g_xxxyyyy_0_xz_0, g_xxxyyyy_0_yy_0, g_xxxyyyy_0_yz_0, g_xxxyyyy_0_zz_0, g_xxyyyy_0_xy_1, g_xxyyyy_0_y_1, g_xxyyyy_0_yy_1, g_xxyyyy_0_yz_1, g_xxyyyy_0_zz_1, g_xyyyy_0_xy_0, g_xyyyy_0_xy_1, g_xyyyy_0_yy_0, g_xyyyy_0_yy_1, g_xyyyy_0_yz_0, g_xyyyy_0_yz_1, g_xyyyy_0_zz_0, g_xyyyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyyy_0_xx_0[i] = 3.0 * g_xxxyy_0_xx_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xx_1[i] * fz_be_0 + g_xxxyyy_0_xx_1[i] * wa_y[i];

        g_xxxyyyy_0_xy_0[i] = 2.0 * g_xyyyy_0_xy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xy_1[i] * fz_be_0 + g_xxyyyy_0_y_1[i] * fi_acd_0 + g_xxyyyy_0_xy_1[i] * wa_x[i];

        g_xxxyyyy_0_xz_0[i] = 3.0 * g_xxxyy_0_xz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xz_1[i] * fz_be_0 + g_xxxyyy_0_xz_1[i] * wa_y[i];

        g_xxxyyyy_0_yy_0[i] = 2.0 * g_xyyyy_0_yy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yy_1[i] * fz_be_0 + g_xxyyyy_0_yy_1[i] * wa_x[i];

        g_xxxyyyy_0_yz_0[i] = 2.0 * g_xyyyy_0_yz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yz_1[i] * fz_be_0 + g_xxyyyy_0_yz_1[i] * wa_x[i];

        g_xxxyyyy_0_zz_0[i] = 2.0 * g_xyyyy_0_zz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_zz_1[i] * fz_be_0 + g_xxyyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 66-72 components of targeted buffer : KSD

    auto g_xxxyyyz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 66);

    auto g_xxxyyyz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 67);

    auto g_xxxyyyz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 68);

    auto g_xxxyyyz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 69);

    auto g_xxxyyyz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 70);

    auto g_xxxyyyz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 71);

    #pragma omp simd aligned(g_xxxyyy_0_x_1, g_xxxyyy_0_xx_1, g_xxxyyy_0_xy_1, g_xxxyyy_0_xz_1, g_xxxyyy_0_y_1, g_xxxyyy_0_yy_1, g_xxxyyy_0_yz_1, g_xxxyyy_0_z_1, g_xxxyyy_0_zz_1, g_xxxyyyz_0_xx_0, g_xxxyyyz_0_xy_0, g_xxxyyyz_0_xz_0, g_xxxyyyz_0_yy_0, g_xxxyyyz_0_yz_0, g_xxxyyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_xx_0[i] = g_xxxyyy_0_xx_1[i] * wa_z[i];

        g_xxxyyyz_0_xy_0[i] = g_xxxyyy_0_xy_1[i] * wa_z[i];

        g_xxxyyyz_0_xz_0[i] = g_xxxyyy_0_x_1[i] * fi_acd_0 + g_xxxyyy_0_xz_1[i] * wa_z[i];

        g_xxxyyyz_0_yy_0[i] = g_xxxyyy_0_yy_1[i] * wa_z[i];

        g_xxxyyyz_0_yz_0[i] = g_xxxyyy_0_y_1[i] * fi_acd_0 + g_xxxyyy_0_yz_1[i] * wa_z[i];

        g_xxxyyyz_0_zz_0[i] = 2.0 * g_xxxyyy_0_z_1[i] * fi_acd_0 + g_xxxyyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 72-78 components of targeted buffer : KSD

    auto g_xxxyyzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 72);

    auto g_xxxyyzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 73);

    auto g_xxxyyzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 74);

    auto g_xxxyyzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 75);

    auto g_xxxyyzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 76);

    auto g_xxxyyzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 77);

    #pragma omp simd aligned(g_xxxyy_0_xy_0, g_xxxyy_0_xy_1, g_xxxyyz_0_xy_1, g_xxxyyzz_0_xx_0, g_xxxyyzz_0_xy_0, g_xxxyyzz_0_xz_0, g_xxxyyzz_0_yy_0, g_xxxyyzz_0_yz_0, g_xxxyyzz_0_zz_0, g_xxxyzz_0_xx_1, g_xxxyzz_0_xz_1, g_xxxzz_0_xx_0, g_xxxzz_0_xx_1, g_xxxzz_0_xz_0, g_xxxzz_0_xz_1, g_xxyyzz_0_yy_1, g_xxyyzz_0_yz_1, g_xxyyzz_0_zz_1, g_xyyzz_0_yy_0, g_xyyzz_0_yy_1, g_xyyzz_0_yz_0, g_xyyzz_0_yz_1, g_xyyzz_0_zz_0, g_xyyzz_0_zz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyzz_0_xx_0[i] = g_xxxzz_0_xx_0[i] * fbe_0 - g_xxxzz_0_xx_1[i] * fz_be_0 + g_xxxyzz_0_xx_1[i] * wa_y[i];

        g_xxxyyzz_0_xy_0[i] = g_xxxyy_0_xy_0[i] * fbe_0 - g_xxxyy_0_xy_1[i] * fz_be_0 + g_xxxyyz_0_xy_1[i] * wa_z[i];

        g_xxxyyzz_0_xz_0[i] = g_xxxzz_0_xz_0[i] * fbe_0 - g_xxxzz_0_xz_1[i] * fz_be_0 + g_xxxyzz_0_xz_1[i] * wa_y[i];

        g_xxxyyzz_0_yy_0[i] = 2.0 * g_xyyzz_0_yy_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yy_1[i] * fz_be_0 + g_xxyyzz_0_yy_1[i] * wa_x[i];

        g_xxxyyzz_0_yz_0[i] = 2.0 * g_xyyzz_0_yz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yz_1[i] * fz_be_0 + g_xxyyzz_0_yz_1[i] * wa_x[i];

        g_xxxyyzz_0_zz_0[i] = 2.0 * g_xyyzz_0_zz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_zz_1[i] * fz_be_0 + g_xxyyzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 78-84 components of targeted buffer : KSD

    auto g_xxxyzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 78);

    auto g_xxxyzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 79);

    auto g_xxxyzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 80);

    auto g_xxxyzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 81);

    auto g_xxxyzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 82);

    auto g_xxxyzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 83);

    #pragma omp simd aligned(g_xxxyzzz_0_xx_0, g_xxxyzzz_0_xy_0, g_xxxyzzz_0_xz_0, g_xxxyzzz_0_yy_0, g_xxxyzzz_0_yz_0, g_xxxyzzz_0_zz_0, g_xxxzzz_0_x_1, g_xxxzzz_0_xx_1, g_xxxzzz_0_xy_1, g_xxxzzz_0_xz_1, g_xxxzzz_0_y_1, g_xxxzzz_0_yy_1, g_xxxzzz_0_yz_1, g_xxxzzz_0_z_1, g_xxxzzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_xx_0[i] = g_xxxzzz_0_xx_1[i] * wa_y[i];

        g_xxxyzzz_0_xy_0[i] = g_xxxzzz_0_x_1[i] * fi_acd_0 + g_xxxzzz_0_xy_1[i] * wa_y[i];

        g_xxxyzzz_0_xz_0[i] = g_xxxzzz_0_xz_1[i] * wa_y[i];

        g_xxxyzzz_0_yy_0[i] = 2.0 * g_xxxzzz_0_y_1[i] * fi_acd_0 + g_xxxzzz_0_yy_1[i] * wa_y[i];

        g_xxxyzzz_0_yz_0[i] = g_xxxzzz_0_z_1[i] * fi_acd_0 + g_xxxzzz_0_yz_1[i] * wa_y[i];

        g_xxxyzzz_0_zz_0[i] = g_xxxzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 84-90 components of targeted buffer : KSD

    auto g_xxxzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 84);

    auto g_xxxzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 85);

    auto g_xxxzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 86);

    auto g_xxxzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 87);

    auto g_xxxzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 88);

    auto g_xxxzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 89);

    #pragma omp simd aligned(g_xxxzz_0_xx_0, g_xxxzz_0_xx_1, g_xxxzz_0_xy_0, g_xxxzz_0_xy_1, g_xxxzzz_0_xx_1, g_xxxzzz_0_xy_1, g_xxxzzzz_0_xx_0, g_xxxzzzz_0_xy_0, g_xxxzzzz_0_xz_0, g_xxxzzzz_0_yy_0, g_xxxzzzz_0_yz_0, g_xxxzzzz_0_zz_0, g_xxzzzz_0_xz_1, g_xxzzzz_0_yy_1, g_xxzzzz_0_yz_1, g_xxzzzz_0_z_1, g_xxzzzz_0_zz_1, g_xzzzz_0_xz_0, g_xzzzz_0_xz_1, g_xzzzz_0_yy_0, g_xzzzz_0_yy_1, g_xzzzz_0_yz_0, g_xzzzz_0_yz_1, g_xzzzz_0_zz_0, g_xzzzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzzz_0_xx_0[i] = 3.0 * g_xxxzz_0_xx_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xx_1[i] * fz_be_0 + g_xxxzzz_0_xx_1[i] * wa_z[i];

        g_xxxzzzz_0_xy_0[i] = 3.0 * g_xxxzz_0_xy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xy_1[i] * fz_be_0 + g_xxxzzz_0_xy_1[i] * wa_z[i];

        g_xxxzzzz_0_xz_0[i] = 2.0 * g_xzzzz_0_xz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xz_1[i] * fz_be_0 + g_xxzzzz_0_z_1[i] * fi_acd_0 + g_xxzzzz_0_xz_1[i] * wa_x[i];

        g_xxxzzzz_0_yy_0[i] = 2.0 * g_xzzzz_0_yy_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yy_1[i] * fz_be_0 + g_xxzzzz_0_yy_1[i] * wa_x[i];

        g_xxxzzzz_0_yz_0[i] = 2.0 * g_xzzzz_0_yz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yz_1[i] * fz_be_0 + g_xxzzzz_0_yz_1[i] * wa_x[i];

        g_xxxzzzz_0_zz_0[i] = 2.0 * g_xzzzz_0_zz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_zz_1[i] * fz_be_0 + g_xxzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 90-96 components of targeted buffer : KSD

    auto g_xxyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 90);

    auto g_xxyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 91);

    auto g_xxyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 92);

    auto g_xxyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 93);

    auto g_xxyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 94);

    auto g_xxyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 95);

    #pragma omp simd aligned(g_xxyyy_0_xx_0, g_xxyyy_0_xx_1, g_xxyyy_0_xz_0, g_xxyyy_0_xz_1, g_xxyyyy_0_xx_1, g_xxyyyy_0_xz_1, g_xxyyyyy_0_xx_0, g_xxyyyyy_0_xy_0, g_xxyyyyy_0_xz_0, g_xxyyyyy_0_yy_0, g_xxyyyyy_0_yz_0, g_xxyyyyy_0_zz_0, g_xyyyyy_0_xy_1, g_xyyyyy_0_y_1, g_xyyyyy_0_yy_1, g_xyyyyy_0_yz_1, g_xyyyyy_0_zz_1, g_yyyyy_0_xy_0, g_yyyyy_0_xy_1, g_yyyyy_0_yy_0, g_yyyyy_0_yy_1, g_yyyyy_0_yz_0, g_yyyyy_0_yz_1, g_yyyyy_0_zz_0, g_yyyyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyyy_0_xx_0[i] = 4.0 * g_xxyyy_0_xx_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xx_1[i] * fz_be_0 + g_xxyyyy_0_xx_1[i] * wa_y[i];

        g_xxyyyyy_0_xy_0[i] = g_yyyyy_0_xy_0[i] * fbe_0 - g_yyyyy_0_xy_1[i] * fz_be_0 + g_xyyyyy_0_y_1[i] * fi_acd_0 + g_xyyyyy_0_xy_1[i] * wa_x[i];

        g_xxyyyyy_0_xz_0[i] = 4.0 * g_xxyyy_0_xz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xz_1[i] * fz_be_0 + g_xxyyyy_0_xz_1[i] * wa_y[i];

        g_xxyyyyy_0_yy_0[i] = g_yyyyy_0_yy_0[i] * fbe_0 - g_yyyyy_0_yy_1[i] * fz_be_0 + g_xyyyyy_0_yy_1[i] * wa_x[i];

        g_xxyyyyy_0_yz_0[i] = g_yyyyy_0_yz_0[i] * fbe_0 - g_yyyyy_0_yz_1[i] * fz_be_0 + g_xyyyyy_0_yz_1[i] * wa_x[i];

        g_xxyyyyy_0_zz_0[i] = g_yyyyy_0_zz_0[i] * fbe_0 - g_yyyyy_0_zz_1[i] * fz_be_0 + g_xyyyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 96-102 components of targeted buffer : KSD

    auto g_xxyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 96);

    auto g_xxyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 97);

    auto g_xxyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 98);

    auto g_xxyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 99);

    auto g_xxyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 100);

    auto g_xxyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 101);

    #pragma omp simd aligned(g_xxyyyy_0_x_1, g_xxyyyy_0_xx_1, g_xxyyyy_0_xy_1, g_xxyyyy_0_xz_1, g_xxyyyy_0_y_1, g_xxyyyy_0_yy_1, g_xxyyyy_0_yz_1, g_xxyyyy_0_z_1, g_xxyyyy_0_zz_1, g_xxyyyyz_0_xx_0, g_xxyyyyz_0_xy_0, g_xxyyyyz_0_xz_0, g_xxyyyyz_0_yy_0, g_xxyyyyz_0_yz_0, g_xxyyyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_xx_0[i] = g_xxyyyy_0_xx_1[i] * wa_z[i];

        g_xxyyyyz_0_xy_0[i] = g_xxyyyy_0_xy_1[i] * wa_z[i];

        g_xxyyyyz_0_xz_0[i] = g_xxyyyy_0_x_1[i] * fi_acd_0 + g_xxyyyy_0_xz_1[i] * wa_z[i];

        g_xxyyyyz_0_yy_0[i] = g_xxyyyy_0_yy_1[i] * wa_z[i];

        g_xxyyyyz_0_yz_0[i] = g_xxyyyy_0_y_1[i] * fi_acd_0 + g_xxyyyy_0_yz_1[i] * wa_z[i];

        g_xxyyyyz_0_zz_0[i] = 2.0 * g_xxyyyy_0_z_1[i] * fi_acd_0 + g_xxyyyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 102-108 components of targeted buffer : KSD

    auto g_xxyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 102);

    auto g_xxyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 103);

    auto g_xxyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 104);

    auto g_xxyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 105);

    auto g_xxyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 106);

    auto g_xxyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 107);

    #pragma omp simd aligned(g_xxyyy_0_xy_0, g_xxyyy_0_xy_1, g_xxyyyz_0_xy_1, g_xxyyyzz_0_xx_0, g_xxyyyzz_0_xy_0, g_xxyyyzz_0_xz_0, g_xxyyyzz_0_yy_0, g_xxyyyzz_0_yz_0, g_xxyyyzz_0_zz_0, g_xxyyzz_0_xx_1, g_xxyyzz_0_xz_1, g_xxyzz_0_xx_0, g_xxyzz_0_xx_1, g_xxyzz_0_xz_0, g_xxyzz_0_xz_1, g_xyyyzz_0_yy_1, g_xyyyzz_0_yz_1, g_xyyyzz_0_zz_1, g_yyyzz_0_yy_0, g_yyyzz_0_yy_1, g_yyyzz_0_yz_0, g_yyyzz_0_yz_1, g_yyyzz_0_zz_0, g_yyyzz_0_zz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyzz_0_xx_0[i] = 2.0 * g_xxyzz_0_xx_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xx_1[i] * fz_be_0 + g_xxyyzz_0_xx_1[i] * wa_y[i];

        g_xxyyyzz_0_xy_0[i] = g_xxyyy_0_xy_0[i] * fbe_0 - g_xxyyy_0_xy_1[i] * fz_be_0 + g_xxyyyz_0_xy_1[i] * wa_z[i];

        g_xxyyyzz_0_xz_0[i] = 2.0 * g_xxyzz_0_xz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xz_1[i] * fz_be_0 + g_xxyyzz_0_xz_1[i] * wa_y[i];

        g_xxyyyzz_0_yy_0[i] = g_yyyzz_0_yy_0[i] * fbe_0 - g_yyyzz_0_yy_1[i] * fz_be_0 + g_xyyyzz_0_yy_1[i] * wa_x[i];

        g_xxyyyzz_0_yz_0[i] = g_yyyzz_0_yz_0[i] * fbe_0 - g_yyyzz_0_yz_1[i] * fz_be_0 + g_xyyyzz_0_yz_1[i] * wa_x[i];

        g_xxyyyzz_0_zz_0[i] = g_yyyzz_0_zz_0[i] * fbe_0 - g_yyyzz_0_zz_1[i] * fz_be_0 + g_xyyyzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 108-114 components of targeted buffer : KSD

    auto g_xxyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 108);

    auto g_xxyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 109);

    auto g_xxyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 110);

    auto g_xxyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 111);

    auto g_xxyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 112);

    auto g_xxyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 113);

    #pragma omp simd aligned(g_xxyyz_0_xy_0, g_xxyyz_0_xy_1, g_xxyyzz_0_xy_1, g_xxyyzzz_0_xx_0, g_xxyyzzz_0_xy_0, g_xxyyzzz_0_xz_0, g_xxyyzzz_0_yy_0, g_xxyyzzz_0_yz_0, g_xxyyzzz_0_zz_0, g_xxyzzz_0_xx_1, g_xxyzzz_0_xz_1, g_xxzzz_0_xx_0, g_xxzzz_0_xx_1, g_xxzzz_0_xz_0, g_xxzzz_0_xz_1, g_xyyzzz_0_yy_1, g_xyyzzz_0_yz_1, g_xyyzzz_0_zz_1, g_yyzzz_0_yy_0, g_yyzzz_0_yy_1, g_yyzzz_0_yz_0, g_yyzzz_0_yz_1, g_yyzzz_0_zz_0, g_yyzzz_0_zz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyzzz_0_xx_0[i] = g_xxzzz_0_xx_0[i] * fbe_0 - g_xxzzz_0_xx_1[i] * fz_be_0 + g_xxyzzz_0_xx_1[i] * wa_y[i];

        g_xxyyzzz_0_xy_0[i] = 2.0 * g_xxyyz_0_xy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xy_1[i] * fz_be_0 + g_xxyyzz_0_xy_1[i] * wa_z[i];

        g_xxyyzzz_0_xz_0[i] = g_xxzzz_0_xz_0[i] * fbe_0 - g_xxzzz_0_xz_1[i] * fz_be_0 + g_xxyzzz_0_xz_1[i] * wa_y[i];

        g_xxyyzzz_0_yy_0[i] = g_yyzzz_0_yy_0[i] * fbe_0 - g_yyzzz_0_yy_1[i] * fz_be_0 + g_xyyzzz_0_yy_1[i] * wa_x[i];

        g_xxyyzzz_0_yz_0[i] = g_yyzzz_0_yz_0[i] * fbe_0 - g_yyzzz_0_yz_1[i] * fz_be_0 + g_xyyzzz_0_yz_1[i] * wa_x[i];

        g_xxyyzzz_0_zz_0[i] = g_yyzzz_0_zz_0[i] * fbe_0 - g_yyzzz_0_zz_1[i] * fz_be_0 + g_xyyzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 114-120 components of targeted buffer : KSD

    auto g_xxyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 114);

    auto g_xxyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 115);

    auto g_xxyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 116);

    auto g_xxyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 117);

    auto g_xxyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 118);

    auto g_xxyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 119);

    #pragma omp simd aligned(g_xxyzzzz_0_xx_0, g_xxyzzzz_0_xy_0, g_xxyzzzz_0_xz_0, g_xxyzzzz_0_yy_0, g_xxyzzzz_0_yz_0, g_xxyzzzz_0_zz_0, g_xxzzzz_0_x_1, g_xxzzzz_0_xx_1, g_xxzzzz_0_xy_1, g_xxzzzz_0_xz_1, g_xxzzzz_0_y_1, g_xxzzzz_0_yy_1, g_xxzzzz_0_yz_1, g_xxzzzz_0_z_1, g_xxzzzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_xx_0[i] = g_xxzzzz_0_xx_1[i] * wa_y[i];

        g_xxyzzzz_0_xy_0[i] = g_xxzzzz_0_x_1[i] * fi_acd_0 + g_xxzzzz_0_xy_1[i] * wa_y[i];

        g_xxyzzzz_0_xz_0[i] = g_xxzzzz_0_xz_1[i] * wa_y[i];

        g_xxyzzzz_0_yy_0[i] = 2.0 * g_xxzzzz_0_y_1[i] * fi_acd_0 + g_xxzzzz_0_yy_1[i] * wa_y[i];

        g_xxyzzzz_0_yz_0[i] = g_xxzzzz_0_z_1[i] * fi_acd_0 + g_xxzzzz_0_yz_1[i] * wa_y[i];

        g_xxyzzzz_0_zz_0[i] = g_xxzzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 120-126 components of targeted buffer : KSD

    auto g_xxzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 120);

    auto g_xxzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 121);

    auto g_xxzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 122);

    auto g_xxzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 123);

    auto g_xxzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 124);

    auto g_xxzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 125);

    #pragma omp simd aligned(g_xxzzz_0_xx_0, g_xxzzz_0_xx_1, g_xxzzz_0_xy_0, g_xxzzz_0_xy_1, g_xxzzzz_0_xx_1, g_xxzzzz_0_xy_1, g_xxzzzzz_0_xx_0, g_xxzzzzz_0_xy_0, g_xxzzzzz_0_xz_0, g_xxzzzzz_0_yy_0, g_xxzzzzz_0_yz_0, g_xxzzzzz_0_zz_0, g_xzzzzz_0_xz_1, g_xzzzzz_0_yy_1, g_xzzzzz_0_yz_1, g_xzzzzz_0_z_1, g_xzzzzz_0_zz_1, g_zzzzz_0_xz_0, g_zzzzz_0_xz_1, g_zzzzz_0_yy_0, g_zzzzz_0_yy_1, g_zzzzz_0_yz_0, g_zzzzz_0_yz_1, g_zzzzz_0_zz_0, g_zzzzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzzz_0_xx_0[i] = 4.0 * g_xxzzz_0_xx_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xx_1[i] * fz_be_0 + g_xxzzzz_0_xx_1[i] * wa_z[i];

        g_xxzzzzz_0_xy_0[i] = 4.0 * g_xxzzz_0_xy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xy_1[i] * fz_be_0 + g_xxzzzz_0_xy_1[i] * wa_z[i];

        g_xxzzzzz_0_xz_0[i] = g_zzzzz_0_xz_0[i] * fbe_0 - g_zzzzz_0_xz_1[i] * fz_be_0 + g_xzzzzz_0_z_1[i] * fi_acd_0 + g_xzzzzz_0_xz_1[i] * wa_x[i];

        g_xxzzzzz_0_yy_0[i] = g_zzzzz_0_yy_0[i] * fbe_0 - g_zzzzz_0_yy_1[i] * fz_be_0 + g_xzzzzz_0_yy_1[i] * wa_x[i];

        g_xxzzzzz_0_yz_0[i] = g_zzzzz_0_yz_0[i] * fbe_0 - g_zzzzz_0_yz_1[i] * fz_be_0 + g_xzzzzz_0_yz_1[i] * wa_x[i];

        g_xxzzzzz_0_zz_0[i] = g_zzzzz_0_zz_0[i] * fbe_0 - g_zzzzz_0_zz_1[i] * fz_be_0 + g_xzzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 126-132 components of targeted buffer : KSD

    auto g_xyyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 126);

    auto g_xyyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 127);

    auto g_xyyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 128);

    auto g_xyyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 129);

    auto g_xyyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 130);

    auto g_xyyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 131);

    #pragma omp simd aligned(g_xyyyyyy_0_xx_0, g_xyyyyyy_0_xy_0, g_xyyyyyy_0_xz_0, g_xyyyyyy_0_yy_0, g_xyyyyyy_0_yz_0, g_xyyyyyy_0_zz_0, g_yyyyyy_0_x_1, g_yyyyyy_0_xx_1, g_yyyyyy_0_xy_1, g_yyyyyy_0_xz_1, g_yyyyyy_0_y_1, g_yyyyyy_0_yy_1, g_yyyyyy_0_yz_1, g_yyyyyy_0_z_1, g_yyyyyy_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_xx_0[i] = 2.0 * g_yyyyyy_0_x_1[i] * fi_acd_0 + g_yyyyyy_0_xx_1[i] * wa_x[i];

        g_xyyyyyy_0_xy_0[i] = g_yyyyyy_0_y_1[i] * fi_acd_0 + g_yyyyyy_0_xy_1[i] * wa_x[i];

        g_xyyyyyy_0_xz_0[i] = g_yyyyyy_0_z_1[i] * fi_acd_0 + g_yyyyyy_0_xz_1[i] * wa_x[i];

        g_xyyyyyy_0_yy_0[i] = g_yyyyyy_0_yy_1[i] * wa_x[i];

        g_xyyyyyy_0_yz_0[i] = g_yyyyyy_0_yz_1[i] * wa_x[i];

        g_xyyyyyy_0_zz_0[i] = g_yyyyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 132-138 components of targeted buffer : KSD

    auto g_xyyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 132);

    auto g_xyyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 133);

    auto g_xyyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 134);

    auto g_xyyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 135);

    auto g_xyyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 136);

    auto g_xyyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 137);

    #pragma omp simd aligned(g_xyyyyy_0_xx_1, g_xyyyyy_0_xy_1, g_xyyyyyz_0_xx_0, g_xyyyyyz_0_xy_0, g_xyyyyyz_0_xz_0, g_xyyyyyz_0_yy_0, g_xyyyyyz_0_yz_0, g_xyyyyyz_0_zz_0, g_yyyyyz_0_xz_1, g_yyyyyz_0_yy_1, g_yyyyyz_0_yz_1, g_yyyyyz_0_z_1, g_yyyyyz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyz_0_xx_0[i] = g_xyyyyy_0_xx_1[i] * wa_z[i];

        g_xyyyyyz_0_xy_0[i] = g_xyyyyy_0_xy_1[i] * wa_z[i];

        g_xyyyyyz_0_xz_0[i] = g_yyyyyz_0_z_1[i] * fi_acd_0 + g_yyyyyz_0_xz_1[i] * wa_x[i];

        g_xyyyyyz_0_yy_0[i] = g_yyyyyz_0_yy_1[i] * wa_x[i];

        g_xyyyyyz_0_yz_0[i] = g_yyyyyz_0_yz_1[i] * wa_x[i];

        g_xyyyyyz_0_zz_0[i] = g_yyyyyz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 138-144 components of targeted buffer : KSD

    auto g_xyyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 138);

    auto g_xyyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 139);

    auto g_xyyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 140);

    auto g_xyyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 141);

    auto g_xyyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 142);

    auto g_xyyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 143);

    #pragma omp simd aligned(g_xyyyyzz_0_xx_0, g_xyyyyzz_0_xy_0, g_xyyyyzz_0_xz_0, g_xyyyyzz_0_yy_0, g_xyyyyzz_0_yz_0, g_xyyyyzz_0_zz_0, g_yyyyzz_0_x_1, g_yyyyzz_0_xx_1, g_yyyyzz_0_xy_1, g_yyyyzz_0_xz_1, g_yyyyzz_0_y_1, g_yyyyzz_0_yy_1, g_yyyyzz_0_yz_1, g_yyyyzz_0_z_1, g_yyyyzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_xx_0[i] = 2.0 * g_yyyyzz_0_x_1[i] * fi_acd_0 + g_yyyyzz_0_xx_1[i] * wa_x[i];

        g_xyyyyzz_0_xy_0[i] = g_yyyyzz_0_y_1[i] * fi_acd_0 + g_yyyyzz_0_xy_1[i] * wa_x[i];

        g_xyyyyzz_0_xz_0[i] = g_yyyyzz_0_z_1[i] * fi_acd_0 + g_yyyyzz_0_xz_1[i] * wa_x[i];

        g_xyyyyzz_0_yy_0[i] = g_yyyyzz_0_yy_1[i] * wa_x[i];

        g_xyyyyzz_0_yz_0[i] = g_yyyyzz_0_yz_1[i] * wa_x[i];

        g_xyyyyzz_0_zz_0[i] = g_yyyyzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 144-150 components of targeted buffer : KSD

    auto g_xyyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 144);

    auto g_xyyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 145);

    auto g_xyyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 146);

    auto g_xyyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 147);

    auto g_xyyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 148);

    auto g_xyyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 149);

    #pragma omp simd aligned(g_xyyyzzz_0_xx_0, g_xyyyzzz_0_xy_0, g_xyyyzzz_0_xz_0, g_xyyyzzz_0_yy_0, g_xyyyzzz_0_yz_0, g_xyyyzzz_0_zz_0, g_yyyzzz_0_x_1, g_yyyzzz_0_xx_1, g_yyyzzz_0_xy_1, g_yyyzzz_0_xz_1, g_yyyzzz_0_y_1, g_yyyzzz_0_yy_1, g_yyyzzz_0_yz_1, g_yyyzzz_0_z_1, g_yyyzzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_xx_0[i] = 2.0 * g_yyyzzz_0_x_1[i] * fi_acd_0 + g_yyyzzz_0_xx_1[i] * wa_x[i];

        g_xyyyzzz_0_xy_0[i] = g_yyyzzz_0_y_1[i] * fi_acd_0 + g_yyyzzz_0_xy_1[i] * wa_x[i];

        g_xyyyzzz_0_xz_0[i] = g_yyyzzz_0_z_1[i] * fi_acd_0 + g_yyyzzz_0_xz_1[i] * wa_x[i];

        g_xyyyzzz_0_yy_0[i] = g_yyyzzz_0_yy_1[i] * wa_x[i];

        g_xyyyzzz_0_yz_0[i] = g_yyyzzz_0_yz_1[i] * wa_x[i];

        g_xyyyzzz_0_zz_0[i] = g_yyyzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 150-156 components of targeted buffer : KSD

    auto g_xyyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 150);

    auto g_xyyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 151);

    auto g_xyyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 152);

    auto g_xyyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 153);

    auto g_xyyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 154);

    auto g_xyyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 155);

    #pragma omp simd aligned(g_xyyzzzz_0_xx_0, g_xyyzzzz_0_xy_0, g_xyyzzzz_0_xz_0, g_xyyzzzz_0_yy_0, g_xyyzzzz_0_yz_0, g_xyyzzzz_0_zz_0, g_yyzzzz_0_x_1, g_yyzzzz_0_xx_1, g_yyzzzz_0_xy_1, g_yyzzzz_0_xz_1, g_yyzzzz_0_y_1, g_yyzzzz_0_yy_1, g_yyzzzz_0_yz_1, g_yyzzzz_0_z_1, g_yyzzzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_xx_0[i] = 2.0 * g_yyzzzz_0_x_1[i] * fi_acd_0 + g_yyzzzz_0_xx_1[i] * wa_x[i];

        g_xyyzzzz_0_xy_0[i] = g_yyzzzz_0_y_1[i] * fi_acd_0 + g_yyzzzz_0_xy_1[i] * wa_x[i];

        g_xyyzzzz_0_xz_0[i] = g_yyzzzz_0_z_1[i] * fi_acd_0 + g_yyzzzz_0_xz_1[i] * wa_x[i];

        g_xyyzzzz_0_yy_0[i] = g_yyzzzz_0_yy_1[i] * wa_x[i];

        g_xyyzzzz_0_yz_0[i] = g_yyzzzz_0_yz_1[i] * wa_x[i];

        g_xyyzzzz_0_zz_0[i] = g_yyzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 156-162 components of targeted buffer : KSD

    auto g_xyzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 156);

    auto g_xyzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 157);

    auto g_xyzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 158);

    auto g_xyzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 159);

    auto g_xyzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 160);

    auto g_xyzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 161);

    #pragma omp simd aligned(g_xyzzzzz_0_xx_0, g_xyzzzzz_0_xy_0, g_xyzzzzz_0_xz_0, g_xyzzzzz_0_yy_0, g_xyzzzzz_0_yz_0, g_xyzzzzz_0_zz_0, g_xzzzzz_0_xx_1, g_xzzzzz_0_xz_1, g_yzzzzz_0_xy_1, g_yzzzzz_0_y_1, g_yzzzzz_0_yy_1, g_yzzzzz_0_yz_1, g_yzzzzz_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzzz_0_xx_0[i] = g_xzzzzz_0_xx_1[i] * wa_y[i];

        g_xyzzzzz_0_xy_0[i] = g_yzzzzz_0_y_1[i] * fi_acd_0 + g_yzzzzz_0_xy_1[i] * wa_x[i];

        g_xyzzzzz_0_xz_0[i] = g_xzzzzz_0_xz_1[i] * wa_y[i];

        g_xyzzzzz_0_yy_0[i] = g_yzzzzz_0_yy_1[i] * wa_x[i];

        g_xyzzzzz_0_yz_0[i] = g_yzzzzz_0_yz_1[i] * wa_x[i];

        g_xyzzzzz_0_zz_0[i] = g_yzzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 162-168 components of targeted buffer : KSD

    auto g_xzzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 162);

    auto g_xzzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 163);

    auto g_xzzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 164);

    auto g_xzzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 165);

    auto g_xzzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 166);

    auto g_xzzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 167);

    #pragma omp simd aligned(g_xzzzzzz_0_xx_0, g_xzzzzzz_0_xy_0, g_xzzzzzz_0_xz_0, g_xzzzzzz_0_yy_0, g_xzzzzzz_0_yz_0, g_xzzzzzz_0_zz_0, g_zzzzzz_0_x_1, g_zzzzzz_0_xx_1, g_zzzzzz_0_xy_1, g_zzzzzz_0_xz_1, g_zzzzzz_0_y_1, g_zzzzzz_0_yy_1, g_zzzzzz_0_yz_1, g_zzzzzz_0_z_1, g_zzzzzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_xx_0[i] = 2.0 * g_zzzzzz_0_x_1[i] * fi_acd_0 + g_zzzzzz_0_xx_1[i] * wa_x[i];

        g_xzzzzzz_0_xy_0[i] = g_zzzzzz_0_y_1[i] * fi_acd_0 + g_zzzzzz_0_xy_1[i] * wa_x[i];

        g_xzzzzzz_0_xz_0[i] = g_zzzzzz_0_z_1[i] * fi_acd_0 + g_zzzzzz_0_xz_1[i] * wa_x[i];

        g_xzzzzzz_0_yy_0[i] = g_zzzzzz_0_yy_1[i] * wa_x[i];

        g_xzzzzzz_0_yz_0[i] = g_zzzzzz_0_yz_1[i] * wa_x[i];

        g_xzzzzzz_0_zz_0[i] = g_zzzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 168-174 components of targeted buffer : KSD

    auto g_yyyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 168);

    auto g_yyyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 169);

    auto g_yyyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 170);

    auto g_yyyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 171);

    auto g_yyyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 172);

    auto g_yyyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 173);

    #pragma omp simd aligned(g_yyyyy_0_xx_0, g_yyyyy_0_xx_1, g_yyyyy_0_xy_0, g_yyyyy_0_xy_1, g_yyyyy_0_xz_0, g_yyyyy_0_xz_1, g_yyyyy_0_yy_0, g_yyyyy_0_yy_1, g_yyyyy_0_yz_0, g_yyyyy_0_yz_1, g_yyyyy_0_zz_0, g_yyyyy_0_zz_1, g_yyyyyy_0_x_1, g_yyyyyy_0_xx_1, g_yyyyyy_0_xy_1, g_yyyyyy_0_xz_1, g_yyyyyy_0_y_1, g_yyyyyy_0_yy_1, g_yyyyyy_0_yz_1, g_yyyyyy_0_z_1, g_yyyyyy_0_zz_1, g_yyyyyyy_0_xx_0, g_yyyyyyy_0_xy_0, g_yyyyyyy_0_xz_0, g_yyyyyyy_0_yy_0, g_yyyyyyy_0_yz_0, g_yyyyyyy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_xx_0[i] = 6.0 * g_yyyyy_0_xx_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xx_1[i] * fz_be_0 + g_yyyyyy_0_xx_1[i] * wa_y[i];

        g_yyyyyyy_0_xy_0[i] = 6.0 * g_yyyyy_0_xy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xy_1[i] * fz_be_0 + g_yyyyyy_0_x_1[i] * fi_acd_0 + g_yyyyyy_0_xy_1[i] * wa_y[i];

        g_yyyyyyy_0_xz_0[i] = 6.0 * g_yyyyy_0_xz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xz_1[i] * fz_be_0 + g_yyyyyy_0_xz_1[i] * wa_y[i];

        g_yyyyyyy_0_yy_0[i] = 6.0 * g_yyyyy_0_yy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_y_1[i] * fi_acd_0 + g_yyyyyy_0_yy_1[i] * wa_y[i];

        g_yyyyyyy_0_yz_0[i] = 6.0 * g_yyyyy_0_yz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yz_1[i] * fz_be_0 + g_yyyyyy_0_z_1[i] * fi_acd_0 + g_yyyyyy_0_yz_1[i] * wa_y[i];

        g_yyyyyyy_0_zz_0[i] = 6.0 * g_yyyyy_0_zz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_zz_1[i] * fz_be_0 + g_yyyyyy_0_zz_1[i] * wa_y[i];
    }

    /// Set up 174-180 components of targeted buffer : KSD

    auto g_yyyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 174);

    auto g_yyyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 175);

    auto g_yyyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 176);

    auto g_yyyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 177);

    auto g_yyyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 178);

    auto g_yyyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 179);

    #pragma omp simd aligned(g_yyyyyy_0_x_1, g_yyyyyy_0_xx_1, g_yyyyyy_0_xy_1, g_yyyyyy_0_xz_1, g_yyyyyy_0_y_1, g_yyyyyy_0_yy_1, g_yyyyyy_0_yz_1, g_yyyyyy_0_z_1, g_yyyyyy_0_zz_1, g_yyyyyyz_0_xx_0, g_yyyyyyz_0_xy_0, g_yyyyyyz_0_xz_0, g_yyyyyyz_0_yy_0, g_yyyyyyz_0_yz_0, g_yyyyyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_xx_0[i] = g_yyyyyy_0_xx_1[i] * wa_z[i];

        g_yyyyyyz_0_xy_0[i] = g_yyyyyy_0_xy_1[i] * wa_z[i];

        g_yyyyyyz_0_xz_0[i] = g_yyyyyy_0_x_1[i] * fi_acd_0 + g_yyyyyy_0_xz_1[i] * wa_z[i];

        g_yyyyyyz_0_yy_0[i] = g_yyyyyy_0_yy_1[i] * wa_z[i];

        g_yyyyyyz_0_yz_0[i] = g_yyyyyy_0_y_1[i] * fi_acd_0 + g_yyyyyy_0_yz_1[i] * wa_z[i];

        g_yyyyyyz_0_zz_0[i] = 2.0 * g_yyyyyy_0_z_1[i] * fi_acd_0 + g_yyyyyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 180-186 components of targeted buffer : KSD

    auto g_yyyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 180);

    auto g_yyyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 181);

    auto g_yyyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 182);

    auto g_yyyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 183);

    auto g_yyyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 184);

    auto g_yyyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 185);

    #pragma omp simd aligned(g_yyyyy_0_xy_0, g_yyyyy_0_xy_1, g_yyyyy_0_yy_0, g_yyyyy_0_yy_1, g_yyyyyz_0_xy_1, g_yyyyyz_0_yy_1, g_yyyyyzz_0_xx_0, g_yyyyyzz_0_xy_0, g_yyyyyzz_0_xz_0, g_yyyyyzz_0_yy_0, g_yyyyyzz_0_yz_0, g_yyyyyzz_0_zz_0, g_yyyyzz_0_xx_1, g_yyyyzz_0_xz_1, g_yyyyzz_0_yz_1, g_yyyyzz_0_z_1, g_yyyyzz_0_zz_1, g_yyyzz_0_xx_0, g_yyyzz_0_xx_1, g_yyyzz_0_xz_0, g_yyyzz_0_xz_1, g_yyyzz_0_yz_0, g_yyyzz_0_yz_1, g_yyyzz_0_zz_0, g_yyyzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyzz_0_xx_0[i] = 4.0 * g_yyyzz_0_xx_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xx_1[i] * fz_be_0 + g_yyyyzz_0_xx_1[i] * wa_y[i];

        g_yyyyyzz_0_xy_0[i] = g_yyyyy_0_xy_0[i] * fbe_0 - g_yyyyy_0_xy_1[i] * fz_be_0 + g_yyyyyz_0_xy_1[i] * wa_z[i];

        g_yyyyyzz_0_xz_0[i] = 4.0 * g_yyyzz_0_xz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xz_1[i] * fz_be_0 + g_yyyyzz_0_xz_1[i] * wa_y[i];

        g_yyyyyzz_0_yy_0[i] = g_yyyyy_0_yy_0[i] * fbe_0 - g_yyyyy_0_yy_1[i] * fz_be_0 + g_yyyyyz_0_yy_1[i] * wa_z[i];

        g_yyyyyzz_0_yz_0[i] = 4.0 * g_yyyzz_0_yz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yz_1[i] * fz_be_0 + g_yyyyzz_0_z_1[i] * fi_acd_0 + g_yyyyzz_0_yz_1[i] * wa_y[i];

        g_yyyyyzz_0_zz_0[i] = 4.0 * g_yyyzz_0_zz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_zz_1[i] * fz_be_0 + g_yyyyzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 186-192 components of targeted buffer : KSD

    auto g_yyyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 186);

    auto g_yyyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 187);

    auto g_yyyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 188);

    auto g_yyyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 189);

    auto g_yyyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 190);

    auto g_yyyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 191);

    #pragma omp simd aligned(g_yyyyz_0_xy_0, g_yyyyz_0_xy_1, g_yyyyz_0_yy_0, g_yyyyz_0_yy_1, g_yyyyzz_0_xy_1, g_yyyyzz_0_yy_1, g_yyyyzzz_0_xx_0, g_yyyyzzz_0_xy_0, g_yyyyzzz_0_xz_0, g_yyyyzzz_0_yy_0, g_yyyyzzz_0_yz_0, g_yyyyzzz_0_zz_0, g_yyyzzz_0_xx_1, g_yyyzzz_0_xz_1, g_yyyzzz_0_yz_1, g_yyyzzz_0_z_1, g_yyyzzz_0_zz_1, g_yyzzz_0_xx_0, g_yyzzz_0_xx_1, g_yyzzz_0_xz_0, g_yyzzz_0_xz_1, g_yyzzz_0_yz_0, g_yyzzz_0_yz_1, g_yyzzz_0_zz_0, g_yyzzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzzz_0_xx_0[i] = 3.0 * g_yyzzz_0_xx_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xx_1[i] * fz_be_0 + g_yyyzzz_0_xx_1[i] * wa_y[i];

        g_yyyyzzz_0_xy_0[i] = 2.0 * g_yyyyz_0_xy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xy_1[i] * fz_be_0 + g_yyyyzz_0_xy_1[i] * wa_z[i];

        g_yyyyzzz_0_xz_0[i] = 3.0 * g_yyzzz_0_xz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xz_1[i] * fz_be_0 + g_yyyzzz_0_xz_1[i] * wa_y[i];

        g_yyyyzzz_0_yy_0[i] = 2.0 * g_yyyyz_0_yy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_yy_1[i] * fz_be_0 + g_yyyyzz_0_yy_1[i] * wa_z[i];

        g_yyyyzzz_0_yz_0[i] = 3.0 * g_yyzzz_0_yz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yz_1[i] * fz_be_0 + g_yyyzzz_0_z_1[i] * fi_acd_0 + g_yyyzzz_0_yz_1[i] * wa_y[i];

        g_yyyyzzz_0_zz_0[i] = 3.0 * g_yyzzz_0_zz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_zz_1[i] * fz_be_0 + g_yyyzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 192-198 components of targeted buffer : KSD

    auto g_yyyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 192);

    auto g_yyyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 193);

    auto g_yyyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 194);

    auto g_yyyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 195);

    auto g_yyyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 196);

    auto g_yyyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 197);

    #pragma omp simd aligned(g_yyyzz_0_xy_0, g_yyyzz_0_xy_1, g_yyyzz_0_yy_0, g_yyyzz_0_yy_1, g_yyyzzz_0_xy_1, g_yyyzzz_0_yy_1, g_yyyzzzz_0_xx_0, g_yyyzzzz_0_xy_0, g_yyyzzzz_0_xz_0, g_yyyzzzz_0_yy_0, g_yyyzzzz_0_yz_0, g_yyyzzzz_0_zz_0, g_yyzzzz_0_xx_1, g_yyzzzz_0_xz_1, g_yyzzzz_0_yz_1, g_yyzzzz_0_z_1, g_yyzzzz_0_zz_1, g_yzzzz_0_xx_0, g_yzzzz_0_xx_1, g_yzzzz_0_xz_0, g_yzzzz_0_xz_1, g_yzzzz_0_yz_0, g_yzzzz_0_yz_1, g_yzzzz_0_zz_0, g_yzzzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzzz_0_xx_0[i] = 2.0 * g_yzzzz_0_xx_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xx_1[i] * fz_be_0 + g_yyzzzz_0_xx_1[i] * wa_y[i];

        g_yyyzzzz_0_xy_0[i] = 3.0 * g_yyyzz_0_xy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xy_1[i] * fz_be_0 + g_yyyzzz_0_xy_1[i] * wa_z[i];

        g_yyyzzzz_0_xz_0[i] = 2.0 * g_yzzzz_0_xz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xz_1[i] * fz_be_0 + g_yyzzzz_0_xz_1[i] * wa_y[i];

        g_yyyzzzz_0_yy_0[i] = 3.0 * g_yyyzz_0_yy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_yy_1[i] * fz_be_0 + g_yyyzzz_0_yy_1[i] * wa_z[i];

        g_yyyzzzz_0_yz_0[i] = 2.0 * g_yzzzz_0_yz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yz_1[i] * fz_be_0 + g_yyzzzz_0_z_1[i] * fi_acd_0 + g_yyzzzz_0_yz_1[i] * wa_y[i];

        g_yyyzzzz_0_zz_0[i] = 2.0 * g_yzzzz_0_zz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_zz_1[i] * fz_be_0 + g_yyzzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 198-204 components of targeted buffer : KSD

    auto g_yyzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 198);

    auto g_yyzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 199);

    auto g_yyzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 200);

    auto g_yyzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 201);

    auto g_yyzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 202);

    auto g_yyzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 203);

    #pragma omp simd aligned(g_yyzzz_0_xy_0, g_yyzzz_0_xy_1, g_yyzzz_0_yy_0, g_yyzzz_0_yy_1, g_yyzzzz_0_xy_1, g_yyzzzz_0_yy_1, g_yyzzzzz_0_xx_0, g_yyzzzzz_0_xy_0, g_yyzzzzz_0_xz_0, g_yyzzzzz_0_yy_0, g_yyzzzzz_0_yz_0, g_yyzzzzz_0_zz_0, g_yzzzzz_0_xx_1, g_yzzzzz_0_xz_1, g_yzzzzz_0_yz_1, g_yzzzzz_0_z_1, g_yzzzzz_0_zz_1, g_zzzzz_0_xx_0, g_zzzzz_0_xx_1, g_zzzzz_0_xz_0, g_zzzzz_0_xz_1, g_zzzzz_0_yz_0, g_zzzzz_0_yz_1, g_zzzzz_0_zz_0, g_zzzzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzzz_0_xx_0[i] = g_zzzzz_0_xx_0[i] * fbe_0 - g_zzzzz_0_xx_1[i] * fz_be_0 + g_yzzzzz_0_xx_1[i] * wa_y[i];

        g_yyzzzzz_0_xy_0[i] = 4.0 * g_yyzzz_0_xy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xy_1[i] * fz_be_0 + g_yyzzzz_0_xy_1[i] * wa_z[i];

        g_yyzzzzz_0_xz_0[i] = g_zzzzz_0_xz_0[i] * fbe_0 - g_zzzzz_0_xz_1[i] * fz_be_0 + g_yzzzzz_0_xz_1[i] * wa_y[i];

        g_yyzzzzz_0_yy_0[i] = 4.0 * g_yyzzz_0_yy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_yy_1[i] * fz_be_0 + g_yyzzzz_0_yy_1[i] * wa_z[i];

        g_yyzzzzz_0_yz_0[i] = g_zzzzz_0_yz_0[i] * fbe_0 - g_zzzzz_0_yz_1[i] * fz_be_0 + g_yzzzzz_0_z_1[i] * fi_acd_0 + g_yzzzzz_0_yz_1[i] * wa_y[i];

        g_yyzzzzz_0_zz_0[i] = g_zzzzz_0_zz_0[i] * fbe_0 - g_zzzzz_0_zz_1[i] * fz_be_0 + g_yzzzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 204-210 components of targeted buffer : KSD

    auto g_yzzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 204);

    auto g_yzzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 205);

    auto g_yzzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 206);

    auto g_yzzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 207);

    auto g_yzzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 208);

    auto g_yzzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 209);

    #pragma omp simd aligned(g_yzzzzzz_0_xx_0, g_yzzzzzz_0_xy_0, g_yzzzzzz_0_xz_0, g_yzzzzzz_0_yy_0, g_yzzzzzz_0_yz_0, g_yzzzzzz_0_zz_0, g_zzzzzz_0_x_1, g_zzzzzz_0_xx_1, g_zzzzzz_0_xy_1, g_zzzzzz_0_xz_1, g_zzzzzz_0_y_1, g_zzzzzz_0_yy_1, g_zzzzzz_0_yz_1, g_zzzzzz_0_z_1, g_zzzzzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_xx_0[i] = g_zzzzzz_0_xx_1[i] * wa_y[i];

        g_yzzzzzz_0_xy_0[i] = g_zzzzzz_0_x_1[i] * fi_acd_0 + g_zzzzzz_0_xy_1[i] * wa_y[i];

        g_yzzzzzz_0_xz_0[i] = g_zzzzzz_0_xz_1[i] * wa_y[i];

        g_yzzzzzz_0_yy_0[i] = 2.0 * g_zzzzzz_0_y_1[i] * fi_acd_0 + g_zzzzzz_0_yy_1[i] * wa_y[i];

        g_yzzzzzz_0_yz_0[i] = g_zzzzzz_0_z_1[i] * fi_acd_0 + g_zzzzzz_0_yz_1[i] * wa_y[i];

        g_yzzzzzz_0_zz_0[i] = g_zzzzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 210-216 components of targeted buffer : KSD

    auto g_zzzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_ksd + 210);

    auto g_zzzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_ksd + 211);

    auto g_zzzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_ksd + 212);

    auto g_zzzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_ksd + 213);

    auto g_zzzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_ksd + 214);

    auto g_zzzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_ksd + 215);

    #pragma omp simd aligned(g_zzzzz_0_xx_0, g_zzzzz_0_xx_1, g_zzzzz_0_xy_0, g_zzzzz_0_xy_1, g_zzzzz_0_xz_0, g_zzzzz_0_xz_1, g_zzzzz_0_yy_0, g_zzzzz_0_yy_1, g_zzzzz_0_yz_0, g_zzzzz_0_yz_1, g_zzzzz_0_zz_0, g_zzzzz_0_zz_1, g_zzzzzz_0_x_1, g_zzzzzz_0_xx_1, g_zzzzzz_0_xy_1, g_zzzzzz_0_xz_1, g_zzzzzz_0_y_1, g_zzzzzz_0_yy_1, g_zzzzzz_0_yz_1, g_zzzzzz_0_z_1, g_zzzzzz_0_zz_1, g_zzzzzzz_0_xx_0, g_zzzzzzz_0_xy_0, g_zzzzzzz_0_xz_0, g_zzzzzzz_0_yy_0, g_zzzzzzz_0_yz_0, g_zzzzzzz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_xx_0[i] = 6.0 * g_zzzzz_0_xx_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xx_1[i] * fz_be_0 + g_zzzzzz_0_xx_1[i] * wa_z[i];

        g_zzzzzzz_0_xy_0[i] = 6.0 * g_zzzzz_0_xy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xy_1[i] * fz_be_0 + g_zzzzzz_0_xy_1[i] * wa_z[i];

        g_zzzzzzz_0_xz_0[i] = 6.0 * g_zzzzz_0_xz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xz_1[i] * fz_be_0 + g_zzzzzz_0_x_1[i] * fi_acd_0 + g_zzzzzz_0_xz_1[i] * wa_z[i];

        g_zzzzzzz_0_yy_0[i] = 6.0 * g_zzzzz_0_yy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yy_1[i] * fz_be_0 + g_zzzzzz_0_yy_1[i] * wa_z[i];

        g_zzzzzzz_0_yz_0[i] = 6.0 * g_zzzzz_0_yz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yz_1[i] * fz_be_0 + g_zzzzzz_0_y_1[i] * fi_acd_0 + g_zzzzzz_0_yz_1[i] * wa_z[i];

        g_zzzzzzz_0_zz_0[i] = 6.0 * g_zzzzz_0_zz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_zz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_z_1[i] * fi_acd_0 + g_zzzzzz_0_zz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

