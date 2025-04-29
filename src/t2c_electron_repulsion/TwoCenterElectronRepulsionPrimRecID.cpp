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

#include "TwoCenterElectronRepulsionPrimRecID.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_id(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_id,
                                const size_t idx_eri_0_gd,
                                const size_t idx_eri_1_gd,
                                const size_t idx_eri_1_hp,
                                const size_t idx_eri_1_hd,
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

    // Set up components of auxiliary buffer : GD

    auto g_xxxx_xx_0 = pbuffer.data(idx_eri_0_gd);

    auto g_xxxx_xy_0 = pbuffer.data(idx_eri_0_gd + 1);

    auto g_xxxx_xz_0 = pbuffer.data(idx_eri_0_gd + 2);

    auto g_xxxx_yy_0 = pbuffer.data(idx_eri_0_gd + 3);

    auto g_xxxx_yz_0 = pbuffer.data(idx_eri_0_gd + 4);

    auto g_xxxx_zz_0 = pbuffer.data(idx_eri_0_gd + 5);

    auto g_xxxy_xx_0 = pbuffer.data(idx_eri_0_gd + 6);

    auto g_xxxy_xz_0 = pbuffer.data(idx_eri_0_gd + 8);

    auto g_xxxz_xx_0 = pbuffer.data(idx_eri_0_gd + 12);

    auto g_xxxz_xy_0 = pbuffer.data(idx_eri_0_gd + 13);

    auto g_xxyy_xx_0 = pbuffer.data(idx_eri_0_gd + 18);

    auto g_xxyy_xy_0 = pbuffer.data(idx_eri_0_gd + 19);

    auto g_xxyy_xz_0 = pbuffer.data(idx_eri_0_gd + 20);

    auto g_xxyy_yy_0 = pbuffer.data(idx_eri_0_gd + 21);

    auto g_xxyy_yz_0 = pbuffer.data(idx_eri_0_gd + 22);

    auto g_xxyy_zz_0 = pbuffer.data(idx_eri_0_gd + 23);

    auto g_xxzz_xx_0 = pbuffer.data(idx_eri_0_gd + 30);

    auto g_xxzz_xy_0 = pbuffer.data(idx_eri_0_gd + 31);

    auto g_xxzz_xz_0 = pbuffer.data(idx_eri_0_gd + 32);

    auto g_xxzz_yy_0 = pbuffer.data(idx_eri_0_gd + 33);

    auto g_xxzz_yz_0 = pbuffer.data(idx_eri_0_gd + 34);

    auto g_xxzz_zz_0 = pbuffer.data(idx_eri_0_gd + 35);

    auto g_xyyy_xy_0 = pbuffer.data(idx_eri_0_gd + 37);

    auto g_xyyy_yy_0 = pbuffer.data(idx_eri_0_gd + 39);

    auto g_xyyy_yz_0 = pbuffer.data(idx_eri_0_gd + 40);

    auto g_xyyy_zz_0 = pbuffer.data(idx_eri_0_gd + 41);

    auto g_xzzz_xz_0 = pbuffer.data(idx_eri_0_gd + 56);

    auto g_xzzz_yy_0 = pbuffer.data(idx_eri_0_gd + 57);

    auto g_xzzz_yz_0 = pbuffer.data(idx_eri_0_gd + 58);

    auto g_xzzz_zz_0 = pbuffer.data(idx_eri_0_gd + 59);

    auto g_yyyy_xx_0 = pbuffer.data(idx_eri_0_gd + 60);

    auto g_yyyy_xy_0 = pbuffer.data(idx_eri_0_gd + 61);

    auto g_yyyy_xz_0 = pbuffer.data(idx_eri_0_gd + 62);

    auto g_yyyy_yy_0 = pbuffer.data(idx_eri_0_gd + 63);

    auto g_yyyy_yz_0 = pbuffer.data(idx_eri_0_gd + 64);

    auto g_yyyy_zz_0 = pbuffer.data(idx_eri_0_gd + 65);

    auto g_yyyz_xy_0 = pbuffer.data(idx_eri_0_gd + 67);

    auto g_yyyz_yy_0 = pbuffer.data(idx_eri_0_gd + 69);

    auto g_yyzz_xx_0 = pbuffer.data(idx_eri_0_gd + 72);

    auto g_yyzz_xy_0 = pbuffer.data(idx_eri_0_gd + 73);

    auto g_yyzz_xz_0 = pbuffer.data(idx_eri_0_gd + 74);

    auto g_yyzz_yy_0 = pbuffer.data(idx_eri_0_gd + 75);

    auto g_yyzz_yz_0 = pbuffer.data(idx_eri_0_gd + 76);

    auto g_yyzz_zz_0 = pbuffer.data(idx_eri_0_gd + 77);

    auto g_yzzz_xx_0 = pbuffer.data(idx_eri_0_gd + 78);

    auto g_yzzz_xz_0 = pbuffer.data(idx_eri_0_gd + 80);

    auto g_yzzz_yz_0 = pbuffer.data(idx_eri_0_gd + 82);

    auto g_yzzz_zz_0 = pbuffer.data(idx_eri_0_gd + 83);

    auto g_zzzz_xx_0 = pbuffer.data(idx_eri_0_gd + 84);

    auto g_zzzz_xy_0 = pbuffer.data(idx_eri_0_gd + 85);

    auto g_zzzz_xz_0 = pbuffer.data(idx_eri_0_gd + 86);

    auto g_zzzz_yy_0 = pbuffer.data(idx_eri_0_gd + 87);

    auto g_zzzz_yz_0 = pbuffer.data(idx_eri_0_gd + 88);

    auto g_zzzz_zz_0 = pbuffer.data(idx_eri_0_gd + 89);

    // Set up components of auxiliary buffer : GD

    auto g_xxxx_xx_1 = pbuffer.data(idx_eri_1_gd);

    auto g_xxxx_xy_1 = pbuffer.data(idx_eri_1_gd + 1);

    auto g_xxxx_xz_1 = pbuffer.data(idx_eri_1_gd + 2);

    auto g_xxxx_yy_1 = pbuffer.data(idx_eri_1_gd + 3);

    auto g_xxxx_yz_1 = pbuffer.data(idx_eri_1_gd + 4);

    auto g_xxxx_zz_1 = pbuffer.data(idx_eri_1_gd + 5);

    auto g_xxxy_xx_1 = pbuffer.data(idx_eri_1_gd + 6);

    auto g_xxxy_xz_1 = pbuffer.data(idx_eri_1_gd + 8);

    auto g_xxxz_xx_1 = pbuffer.data(idx_eri_1_gd + 12);

    auto g_xxxz_xy_1 = pbuffer.data(idx_eri_1_gd + 13);

    auto g_xxyy_xx_1 = pbuffer.data(idx_eri_1_gd + 18);

    auto g_xxyy_xy_1 = pbuffer.data(idx_eri_1_gd + 19);

    auto g_xxyy_xz_1 = pbuffer.data(idx_eri_1_gd + 20);

    auto g_xxyy_yy_1 = pbuffer.data(idx_eri_1_gd + 21);

    auto g_xxyy_yz_1 = pbuffer.data(idx_eri_1_gd + 22);

    auto g_xxyy_zz_1 = pbuffer.data(idx_eri_1_gd + 23);

    auto g_xxzz_xx_1 = pbuffer.data(idx_eri_1_gd + 30);

    auto g_xxzz_xy_1 = pbuffer.data(idx_eri_1_gd + 31);

    auto g_xxzz_xz_1 = pbuffer.data(idx_eri_1_gd + 32);

    auto g_xxzz_yy_1 = pbuffer.data(idx_eri_1_gd + 33);

    auto g_xxzz_yz_1 = pbuffer.data(idx_eri_1_gd + 34);

    auto g_xxzz_zz_1 = pbuffer.data(idx_eri_1_gd + 35);

    auto g_xyyy_xy_1 = pbuffer.data(idx_eri_1_gd + 37);

    auto g_xyyy_yy_1 = pbuffer.data(idx_eri_1_gd + 39);

    auto g_xyyy_yz_1 = pbuffer.data(idx_eri_1_gd + 40);

    auto g_xyyy_zz_1 = pbuffer.data(idx_eri_1_gd + 41);

    auto g_xzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 56);

    auto g_xzzz_yy_1 = pbuffer.data(idx_eri_1_gd + 57);

    auto g_xzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 58);

    auto g_xzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 59);

    auto g_yyyy_xx_1 = pbuffer.data(idx_eri_1_gd + 60);

    auto g_yyyy_xy_1 = pbuffer.data(idx_eri_1_gd + 61);

    auto g_yyyy_xz_1 = pbuffer.data(idx_eri_1_gd + 62);

    auto g_yyyy_yy_1 = pbuffer.data(idx_eri_1_gd + 63);

    auto g_yyyy_yz_1 = pbuffer.data(idx_eri_1_gd + 64);

    auto g_yyyy_zz_1 = pbuffer.data(idx_eri_1_gd + 65);

    auto g_yyyz_xy_1 = pbuffer.data(idx_eri_1_gd + 67);

    auto g_yyyz_yy_1 = pbuffer.data(idx_eri_1_gd + 69);

    auto g_yyzz_xx_1 = pbuffer.data(idx_eri_1_gd + 72);

    auto g_yyzz_xy_1 = pbuffer.data(idx_eri_1_gd + 73);

    auto g_yyzz_xz_1 = pbuffer.data(idx_eri_1_gd + 74);

    auto g_yyzz_yy_1 = pbuffer.data(idx_eri_1_gd + 75);

    auto g_yyzz_yz_1 = pbuffer.data(idx_eri_1_gd + 76);

    auto g_yyzz_zz_1 = pbuffer.data(idx_eri_1_gd + 77);

    auto g_yzzz_xx_1 = pbuffer.data(idx_eri_1_gd + 78);

    auto g_yzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 80);

    auto g_yzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 82);

    auto g_yzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 83);

    auto g_zzzz_xx_1 = pbuffer.data(idx_eri_1_gd + 84);

    auto g_zzzz_xy_1 = pbuffer.data(idx_eri_1_gd + 85);

    auto g_zzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 86);

    auto g_zzzz_yy_1 = pbuffer.data(idx_eri_1_gd + 87);

    auto g_zzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 88);

    auto g_zzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 89);

    // Set up components of auxiliary buffer : HP

    auto g_xxxxx_x_1 = pbuffer.data(idx_eri_1_hp);

    auto g_xxxxx_y_1 = pbuffer.data(idx_eri_1_hp + 1);

    auto g_xxxxx_z_1 = pbuffer.data(idx_eri_1_hp + 2);

    auto g_xxxxz_z_1 = pbuffer.data(idx_eri_1_hp + 8);

    auto g_xxxyy_x_1 = pbuffer.data(idx_eri_1_hp + 9);

    auto g_xxxyy_y_1 = pbuffer.data(idx_eri_1_hp + 10);

    auto g_xxxyy_z_1 = pbuffer.data(idx_eri_1_hp + 11);

    auto g_xxxzz_x_1 = pbuffer.data(idx_eri_1_hp + 15);

    auto g_xxxzz_y_1 = pbuffer.data(idx_eri_1_hp + 16);

    auto g_xxxzz_z_1 = pbuffer.data(idx_eri_1_hp + 17);

    auto g_xxyyy_x_1 = pbuffer.data(idx_eri_1_hp + 18);

    auto g_xxyyy_y_1 = pbuffer.data(idx_eri_1_hp + 19);

    auto g_xxyyy_z_1 = pbuffer.data(idx_eri_1_hp + 20);

    auto g_xxzzz_x_1 = pbuffer.data(idx_eri_1_hp + 27);

    auto g_xxzzz_y_1 = pbuffer.data(idx_eri_1_hp + 28);

    auto g_xxzzz_z_1 = pbuffer.data(idx_eri_1_hp + 29);

    auto g_xyyyy_y_1 = pbuffer.data(idx_eri_1_hp + 31);

    auto g_xzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 44);

    auto g_yyyyy_x_1 = pbuffer.data(idx_eri_1_hp + 45);

    auto g_yyyyy_y_1 = pbuffer.data(idx_eri_1_hp + 46);

    auto g_yyyyy_z_1 = pbuffer.data(idx_eri_1_hp + 47);

    auto g_yyyyz_z_1 = pbuffer.data(idx_eri_1_hp + 50);

    auto g_yyyzz_x_1 = pbuffer.data(idx_eri_1_hp + 51);

    auto g_yyyzz_y_1 = pbuffer.data(idx_eri_1_hp + 52);

    auto g_yyyzz_z_1 = pbuffer.data(idx_eri_1_hp + 53);

    auto g_yyzzz_x_1 = pbuffer.data(idx_eri_1_hp + 54);

    auto g_yyzzz_y_1 = pbuffer.data(idx_eri_1_hp + 55);

    auto g_yyzzz_z_1 = pbuffer.data(idx_eri_1_hp + 56);

    auto g_yzzzz_y_1 = pbuffer.data(idx_eri_1_hp + 58);

    auto g_yzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 59);

    auto g_zzzzz_x_1 = pbuffer.data(idx_eri_1_hp + 60);

    auto g_zzzzz_y_1 = pbuffer.data(idx_eri_1_hp + 61);

    auto g_zzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 62);

    // Set up components of auxiliary buffer : HD

    auto g_xxxxx_xx_1 = pbuffer.data(idx_eri_1_hd);

    auto g_xxxxx_xy_1 = pbuffer.data(idx_eri_1_hd + 1);

    auto g_xxxxx_xz_1 = pbuffer.data(idx_eri_1_hd + 2);

    auto g_xxxxx_yy_1 = pbuffer.data(idx_eri_1_hd + 3);

    auto g_xxxxx_yz_1 = pbuffer.data(idx_eri_1_hd + 4);

    auto g_xxxxx_zz_1 = pbuffer.data(idx_eri_1_hd + 5);

    auto g_xxxxy_xx_1 = pbuffer.data(idx_eri_1_hd + 6);

    auto g_xxxxy_xy_1 = pbuffer.data(idx_eri_1_hd + 7);

    auto g_xxxxy_xz_1 = pbuffer.data(idx_eri_1_hd + 8);

    auto g_xxxxy_yy_1 = pbuffer.data(idx_eri_1_hd + 9);

    auto g_xxxxz_xx_1 = pbuffer.data(idx_eri_1_hd + 12);

    auto g_xxxxz_xy_1 = pbuffer.data(idx_eri_1_hd + 13);

    auto g_xxxxz_xz_1 = pbuffer.data(idx_eri_1_hd + 14);

    auto g_xxxxz_yz_1 = pbuffer.data(idx_eri_1_hd + 16);

    auto g_xxxxz_zz_1 = pbuffer.data(idx_eri_1_hd + 17);

    auto g_xxxyy_xx_1 = pbuffer.data(idx_eri_1_hd + 18);

    auto g_xxxyy_xy_1 = pbuffer.data(idx_eri_1_hd + 19);

    auto g_xxxyy_xz_1 = pbuffer.data(idx_eri_1_hd + 20);

    auto g_xxxyy_yy_1 = pbuffer.data(idx_eri_1_hd + 21);

    auto g_xxxyy_yz_1 = pbuffer.data(idx_eri_1_hd + 22);

    auto g_xxxyy_zz_1 = pbuffer.data(idx_eri_1_hd + 23);

    auto g_xxxzz_xx_1 = pbuffer.data(idx_eri_1_hd + 30);

    auto g_xxxzz_xy_1 = pbuffer.data(idx_eri_1_hd + 31);

    auto g_xxxzz_xz_1 = pbuffer.data(idx_eri_1_hd + 32);

    auto g_xxxzz_yy_1 = pbuffer.data(idx_eri_1_hd + 33);

    auto g_xxxzz_yz_1 = pbuffer.data(idx_eri_1_hd + 34);

    auto g_xxxzz_zz_1 = pbuffer.data(idx_eri_1_hd + 35);

    auto g_xxyyy_xx_1 = pbuffer.data(idx_eri_1_hd + 36);

    auto g_xxyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 37);

    auto g_xxyyy_xz_1 = pbuffer.data(idx_eri_1_hd + 38);

    auto g_xxyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 39);

    auto g_xxyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 40);

    auto g_xxyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 41);

    auto g_xxyyz_xy_1 = pbuffer.data(idx_eri_1_hd + 43);

    auto g_xxyzz_xx_1 = pbuffer.data(idx_eri_1_hd + 48);

    auto g_xxyzz_xz_1 = pbuffer.data(idx_eri_1_hd + 50);

    auto g_xxzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 54);

    auto g_xxzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 55);

    auto g_xxzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 56);

    auto g_xxzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 57);

    auto g_xxzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 58);

    auto g_xxzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 59);

    auto g_xyyyy_xx_1 = pbuffer.data(idx_eri_1_hd + 60);

    auto g_xyyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 61);

    auto g_xyyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 63);

    auto g_xyyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 64);

    auto g_xyyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 65);

    auto g_xyyzz_yy_1 = pbuffer.data(idx_eri_1_hd + 75);

    auto g_xyyzz_yz_1 = pbuffer.data(idx_eri_1_hd + 76);

    auto g_xyyzz_zz_1 = pbuffer.data(idx_eri_1_hd + 77);

    auto g_xzzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 84);

    auto g_xzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 86);

    auto g_xzzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 87);

    auto g_xzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 88);

    auto g_xzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 89);

    auto g_yyyyy_xx_1 = pbuffer.data(idx_eri_1_hd + 90);

    auto g_yyyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 91);

    auto g_yyyyy_xz_1 = pbuffer.data(idx_eri_1_hd + 92);

    auto g_yyyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 93);

    auto g_yyyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 94);

    auto g_yyyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 95);

    auto g_yyyyz_xy_1 = pbuffer.data(idx_eri_1_hd + 97);

    auto g_yyyyz_xz_1 = pbuffer.data(idx_eri_1_hd + 98);

    auto g_yyyyz_yy_1 = pbuffer.data(idx_eri_1_hd + 99);

    auto g_yyyyz_yz_1 = pbuffer.data(idx_eri_1_hd + 100);

    auto g_yyyyz_zz_1 = pbuffer.data(idx_eri_1_hd + 101);

    auto g_yyyzz_xx_1 = pbuffer.data(idx_eri_1_hd + 102);

    auto g_yyyzz_xy_1 = pbuffer.data(idx_eri_1_hd + 103);

    auto g_yyyzz_xz_1 = pbuffer.data(idx_eri_1_hd + 104);

    auto g_yyyzz_yy_1 = pbuffer.data(idx_eri_1_hd + 105);

    auto g_yyyzz_yz_1 = pbuffer.data(idx_eri_1_hd + 106);

    auto g_yyyzz_zz_1 = pbuffer.data(idx_eri_1_hd + 107);

    auto g_yyzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 108);

    auto g_yyzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 109);

    auto g_yyzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 110);

    auto g_yyzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 111);

    auto g_yyzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 112);

    auto g_yyzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 113);

    auto g_yzzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 114);

    auto g_yzzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 115);

    auto g_yzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 116);

    auto g_yzzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 117);

    auto g_yzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 118);

    auto g_yzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 119);

    auto g_zzzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 120);

    auto g_zzzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 121);

    auto g_zzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 122);

    auto g_zzzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 123);

    auto g_zzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 124);

    auto g_zzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 125);

    // Set up 0-6 components of targeted buffer : ID

    auto g_xxxxxx_xx_0 = pbuffer.data(idx_eri_0_id);

    auto g_xxxxxx_xy_0 = pbuffer.data(idx_eri_0_id + 1);

    auto g_xxxxxx_xz_0 = pbuffer.data(idx_eri_0_id + 2);

    auto g_xxxxxx_yy_0 = pbuffer.data(idx_eri_0_id + 3);

    auto g_xxxxxx_yz_0 = pbuffer.data(idx_eri_0_id + 4);

    auto g_xxxxxx_zz_0 = pbuffer.data(idx_eri_0_id + 5);

    #pragma omp simd aligned(g_xxxx_xx_0, g_xxxx_xx_1, g_xxxx_xy_0, g_xxxx_xy_1, g_xxxx_xz_0, g_xxxx_xz_1, g_xxxx_yy_0, g_xxxx_yy_1, g_xxxx_yz_0, g_xxxx_yz_1, g_xxxx_zz_0, g_xxxx_zz_1, g_xxxxx_x_1, g_xxxxx_xx_1, g_xxxxx_xy_1, g_xxxxx_xz_1, g_xxxxx_y_1, g_xxxxx_yy_1, g_xxxxx_yz_1, g_xxxxx_z_1, g_xxxxx_zz_1, g_xxxxxx_xx_0, g_xxxxxx_xy_0, g_xxxxxx_xz_0, g_xxxxxx_yy_0, g_xxxxxx_yz_0, g_xxxxxx_zz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxx_xx_0[i] = 5.0 * g_xxxx_xx_0[i] * fbe_0 - 5.0 * g_xxxx_xx_1[i] * fz_be_0 + 2.0 * g_xxxxx_x_1[i] * fe_0 + g_xxxxx_xx_1[i] * pa_x[i];

        g_xxxxxx_xy_0[i] = 5.0 * g_xxxx_xy_0[i] * fbe_0 - 5.0 * g_xxxx_xy_1[i] * fz_be_0 + g_xxxxx_y_1[i] * fe_0 + g_xxxxx_xy_1[i] * pa_x[i];

        g_xxxxxx_xz_0[i] = 5.0 * g_xxxx_xz_0[i] * fbe_0 - 5.0 * g_xxxx_xz_1[i] * fz_be_0 + g_xxxxx_z_1[i] * fe_0 + g_xxxxx_xz_1[i] * pa_x[i];

        g_xxxxxx_yy_0[i] = 5.0 * g_xxxx_yy_0[i] * fbe_0 - 5.0 * g_xxxx_yy_1[i] * fz_be_0 + g_xxxxx_yy_1[i] * pa_x[i];

        g_xxxxxx_yz_0[i] = 5.0 * g_xxxx_yz_0[i] * fbe_0 - 5.0 * g_xxxx_yz_1[i] * fz_be_0 + g_xxxxx_yz_1[i] * pa_x[i];

        g_xxxxxx_zz_0[i] = 5.0 * g_xxxx_zz_0[i] * fbe_0 - 5.0 * g_xxxx_zz_1[i] * fz_be_0 + g_xxxxx_zz_1[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : ID

    auto g_xxxxxy_xx_0 = pbuffer.data(idx_eri_0_id + 6);

    auto g_xxxxxy_xy_0 = pbuffer.data(idx_eri_0_id + 7);

    auto g_xxxxxy_xz_0 = pbuffer.data(idx_eri_0_id + 8);

    auto g_xxxxxy_yy_0 = pbuffer.data(idx_eri_0_id + 9);

    auto g_xxxxxy_yz_0 = pbuffer.data(idx_eri_0_id + 10);

    auto g_xxxxxy_zz_0 = pbuffer.data(idx_eri_0_id + 11);

    #pragma omp simd aligned(g_xxxxx_x_1, g_xxxxx_xx_1, g_xxxxx_xy_1, g_xxxxx_xz_1, g_xxxxx_y_1, g_xxxxx_yy_1, g_xxxxx_yz_1, g_xxxxx_z_1, g_xxxxx_zz_1, g_xxxxxy_xx_0, g_xxxxxy_xy_0, g_xxxxxy_xz_0, g_xxxxxy_yy_0, g_xxxxxy_yz_0, g_xxxxxy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxy_xx_0[i] = g_xxxxx_xx_1[i] * pa_y[i];

        g_xxxxxy_xy_0[i] = g_xxxxx_x_1[i] * fe_0 + g_xxxxx_xy_1[i] * pa_y[i];

        g_xxxxxy_xz_0[i] = g_xxxxx_xz_1[i] * pa_y[i];

        g_xxxxxy_yy_0[i] = 2.0 * g_xxxxx_y_1[i] * fe_0 + g_xxxxx_yy_1[i] * pa_y[i];

        g_xxxxxy_yz_0[i] = g_xxxxx_z_1[i] * fe_0 + g_xxxxx_yz_1[i] * pa_y[i];

        g_xxxxxy_zz_0[i] = g_xxxxx_zz_1[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : ID

    auto g_xxxxxz_xx_0 = pbuffer.data(idx_eri_0_id + 12);

    auto g_xxxxxz_xy_0 = pbuffer.data(idx_eri_0_id + 13);

    auto g_xxxxxz_xz_0 = pbuffer.data(idx_eri_0_id + 14);

    auto g_xxxxxz_yy_0 = pbuffer.data(idx_eri_0_id + 15);

    auto g_xxxxxz_yz_0 = pbuffer.data(idx_eri_0_id + 16);

    auto g_xxxxxz_zz_0 = pbuffer.data(idx_eri_0_id + 17);

    #pragma omp simd aligned(g_xxxxx_x_1, g_xxxxx_xx_1, g_xxxxx_xy_1, g_xxxxx_xz_1, g_xxxxx_y_1, g_xxxxx_yy_1, g_xxxxx_yz_1, g_xxxxx_z_1, g_xxxxx_zz_1, g_xxxxxz_xx_0, g_xxxxxz_xy_0, g_xxxxxz_xz_0, g_xxxxxz_yy_0, g_xxxxxz_yz_0, g_xxxxxz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxz_xx_0[i] = g_xxxxx_xx_1[i] * pa_z[i];

        g_xxxxxz_xy_0[i] = g_xxxxx_xy_1[i] * pa_z[i];

        g_xxxxxz_xz_0[i] = g_xxxxx_x_1[i] * fe_0 + g_xxxxx_xz_1[i] * pa_z[i];

        g_xxxxxz_yy_0[i] = g_xxxxx_yy_1[i] * pa_z[i];

        g_xxxxxz_yz_0[i] = g_xxxxx_y_1[i] * fe_0 + g_xxxxx_yz_1[i] * pa_z[i];

        g_xxxxxz_zz_0[i] = 2.0 * g_xxxxx_z_1[i] * fe_0 + g_xxxxx_zz_1[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : ID

    auto g_xxxxyy_xx_0 = pbuffer.data(idx_eri_0_id + 18);

    auto g_xxxxyy_xy_0 = pbuffer.data(idx_eri_0_id + 19);

    auto g_xxxxyy_xz_0 = pbuffer.data(idx_eri_0_id + 20);

    auto g_xxxxyy_yy_0 = pbuffer.data(idx_eri_0_id + 21);

    auto g_xxxxyy_yz_0 = pbuffer.data(idx_eri_0_id + 22);

    auto g_xxxxyy_zz_0 = pbuffer.data(idx_eri_0_id + 23);

    #pragma omp simd aligned(g_xxxx_xx_0, g_xxxx_xx_1, g_xxxx_xz_0, g_xxxx_xz_1, g_xxxxy_xx_1, g_xxxxy_xz_1, g_xxxxyy_xx_0, g_xxxxyy_xy_0, g_xxxxyy_xz_0, g_xxxxyy_yy_0, g_xxxxyy_yz_0, g_xxxxyy_zz_0, g_xxxyy_xy_1, g_xxxyy_y_1, g_xxxyy_yy_1, g_xxxyy_yz_1, g_xxxyy_zz_1, g_xxyy_xy_0, g_xxyy_xy_1, g_xxyy_yy_0, g_xxyy_yy_1, g_xxyy_yz_0, g_xxyy_yz_1, g_xxyy_zz_0, g_xxyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyy_xx_0[i] = g_xxxx_xx_0[i] * fbe_0 - g_xxxx_xx_1[i] * fz_be_0 + g_xxxxy_xx_1[i] * pa_y[i];

        g_xxxxyy_xy_0[i] = 3.0 * g_xxyy_xy_0[i] * fbe_0 - 3.0 * g_xxyy_xy_1[i] * fz_be_0 + g_xxxyy_y_1[i] * fe_0 + g_xxxyy_xy_1[i] * pa_x[i];

        g_xxxxyy_xz_0[i] = g_xxxx_xz_0[i] * fbe_0 - g_xxxx_xz_1[i] * fz_be_0 + g_xxxxy_xz_1[i] * pa_y[i];

        g_xxxxyy_yy_0[i] = 3.0 * g_xxyy_yy_0[i] * fbe_0 - 3.0 * g_xxyy_yy_1[i] * fz_be_0 + g_xxxyy_yy_1[i] * pa_x[i];

        g_xxxxyy_yz_0[i] = 3.0 * g_xxyy_yz_0[i] * fbe_0 - 3.0 * g_xxyy_yz_1[i] * fz_be_0 + g_xxxyy_yz_1[i] * pa_x[i];

        g_xxxxyy_zz_0[i] = 3.0 * g_xxyy_zz_0[i] * fbe_0 - 3.0 * g_xxyy_zz_1[i] * fz_be_0 + g_xxxyy_zz_1[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : ID

    auto g_xxxxyz_xx_0 = pbuffer.data(idx_eri_0_id + 24);

    auto g_xxxxyz_xy_0 = pbuffer.data(idx_eri_0_id + 25);

    auto g_xxxxyz_xz_0 = pbuffer.data(idx_eri_0_id + 26);

    auto g_xxxxyz_yy_0 = pbuffer.data(idx_eri_0_id + 27);

    auto g_xxxxyz_yz_0 = pbuffer.data(idx_eri_0_id + 28);

    auto g_xxxxyz_zz_0 = pbuffer.data(idx_eri_0_id + 29);

    #pragma omp simd aligned(g_xxxxy_xy_1, g_xxxxy_yy_1, g_xxxxyz_xx_0, g_xxxxyz_xy_0, g_xxxxyz_xz_0, g_xxxxyz_yy_0, g_xxxxyz_yz_0, g_xxxxyz_zz_0, g_xxxxz_xx_1, g_xxxxz_xz_1, g_xxxxz_yz_1, g_xxxxz_z_1, g_xxxxz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyz_xx_0[i] = g_xxxxz_xx_1[i] * pa_y[i];

        g_xxxxyz_xy_0[i] = g_xxxxy_xy_1[i] * pa_z[i];

        g_xxxxyz_xz_0[i] = g_xxxxz_xz_1[i] * pa_y[i];

        g_xxxxyz_yy_0[i] = g_xxxxy_yy_1[i] * pa_z[i];

        g_xxxxyz_yz_0[i] = g_xxxxz_z_1[i] * fe_0 + g_xxxxz_yz_1[i] * pa_y[i];

        g_xxxxyz_zz_0[i] = g_xxxxz_zz_1[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : ID

    auto g_xxxxzz_xx_0 = pbuffer.data(idx_eri_0_id + 30);

    auto g_xxxxzz_xy_0 = pbuffer.data(idx_eri_0_id + 31);

    auto g_xxxxzz_xz_0 = pbuffer.data(idx_eri_0_id + 32);

    auto g_xxxxzz_yy_0 = pbuffer.data(idx_eri_0_id + 33);

    auto g_xxxxzz_yz_0 = pbuffer.data(idx_eri_0_id + 34);

    auto g_xxxxzz_zz_0 = pbuffer.data(idx_eri_0_id + 35);

    #pragma omp simd aligned(g_xxxx_xx_0, g_xxxx_xx_1, g_xxxx_xy_0, g_xxxx_xy_1, g_xxxxz_xx_1, g_xxxxz_xy_1, g_xxxxzz_xx_0, g_xxxxzz_xy_0, g_xxxxzz_xz_0, g_xxxxzz_yy_0, g_xxxxzz_yz_0, g_xxxxzz_zz_0, g_xxxzz_xz_1, g_xxxzz_yy_1, g_xxxzz_yz_1, g_xxxzz_z_1, g_xxxzz_zz_1, g_xxzz_xz_0, g_xxzz_xz_1, g_xxzz_yy_0, g_xxzz_yy_1, g_xxzz_yz_0, g_xxzz_yz_1, g_xxzz_zz_0, g_xxzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzz_xx_0[i] = g_xxxx_xx_0[i] * fbe_0 - g_xxxx_xx_1[i] * fz_be_0 + g_xxxxz_xx_1[i] * pa_z[i];

        g_xxxxzz_xy_0[i] = g_xxxx_xy_0[i] * fbe_0 - g_xxxx_xy_1[i] * fz_be_0 + g_xxxxz_xy_1[i] * pa_z[i];

        g_xxxxzz_xz_0[i] = 3.0 * g_xxzz_xz_0[i] * fbe_0 - 3.0 * g_xxzz_xz_1[i] * fz_be_0 + g_xxxzz_z_1[i] * fe_0 + g_xxxzz_xz_1[i] * pa_x[i];

        g_xxxxzz_yy_0[i] = 3.0 * g_xxzz_yy_0[i] * fbe_0 - 3.0 * g_xxzz_yy_1[i] * fz_be_0 + g_xxxzz_yy_1[i] * pa_x[i];

        g_xxxxzz_yz_0[i] = 3.0 * g_xxzz_yz_0[i] * fbe_0 - 3.0 * g_xxzz_yz_1[i] * fz_be_0 + g_xxxzz_yz_1[i] * pa_x[i];

        g_xxxxzz_zz_0[i] = 3.0 * g_xxzz_zz_0[i] * fbe_0 - 3.0 * g_xxzz_zz_1[i] * fz_be_0 + g_xxxzz_zz_1[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : ID

    auto g_xxxyyy_xx_0 = pbuffer.data(idx_eri_0_id + 36);

    auto g_xxxyyy_xy_0 = pbuffer.data(idx_eri_0_id + 37);

    auto g_xxxyyy_xz_0 = pbuffer.data(idx_eri_0_id + 38);

    auto g_xxxyyy_yy_0 = pbuffer.data(idx_eri_0_id + 39);

    auto g_xxxyyy_yz_0 = pbuffer.data(idx_eri_0_id + 40);

    auto g_xxxyyy_zz_0 = pbuffer.data(idx_eri_0_id + 41);

    #pragma omp simd aligned(g_xxxy_xx_0, g_xxxy_xx_1, g_xxxy_xz_0, g_xxxy_xz_1, g_xxxyy_xx_1, g_xxxyy_xz_1, g_xxxyyy_xx_0, g_xxxyyy_xy_0, g_xxxyyy_xz_0, g_xxxyyy_yy_0, g_xxxyyy_yz_0, g_xxxyyy_zz_0, g_xxyyy_xy_1, g_xxyyy_y_1, g_xxyyy_yy_1, g_xxyyy_yz_1, g_xxyyy_zz_1, g_xyyy_xy_0, g_xyyy_xy_1, g_xyyy_yy_0, g_xyyy_yy_1, g_xyyy_yz_0, g_xyyy_yz_1, g_xyyy_zz_0, g_xyyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyy_xx_0[i] = 2.0 * g_xxxy_xx_0[i] * fbe_0 - 2.0 * g_xxxy_xx_1[i] * fz_be_0 + g_xxxyy_xx_1[i] * pa_y[i];

        g_xxxyyy_xy_0[i] = 2.0 * g_xyyy_xy_0[i] * fbe_0 - 2.0 * g_xyyy_xy_1[i] * fz_be_0 + g_xxyyy_y_1[i] * fe_0 + g_xxyyy_xy_1[i] * pa_x[i];

        g_xxxyyy_xz_0[i] = 2.0 * g_xxxy_xz_0[i] * fbe_0 - 2.0 * g_xxxy_xz_1[i] * fz_be_0 + g_xxxyy_xz_1[i] * pa_y[i];

        g_xxxyyy_yy_0[i] = 2.0 * g_xyyy_yy_0[i] * fbe_0 - 2.0 * g_xyyy_yy_1[i] * fz_be_0 + g_xxyyy_yy_1[i] * pa_x[i];

        g_xxxyyy_yz_0[i] = 2.0 * g_xyyy_yz_0[i] * fbe_0 - 2.0 * g_xyyy_yz_1[i] * fz_be_0 + g_xxyyy_yz_1[i] * pa_x[i];

        g_xxxyyy_zz_0[i] = 2.0 * g_xyyy_zz_0[i] * fbe_0 - 2.0 * g_xyyy_zz_1[i] * fz_be_0 + g_xxyyy_zz_1[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : ID

    auto g_xxxyyz_xx_0 = pbuffer.data(idx_eri_0_id + 42);

    auto g_xxxyyz_xy_0 = pbuffer.data(idx_eri_0_id + 43);

    auto g_xxxyyz_xz_0 = pbuffer.data(idx_eri_0_id + 44);

    auto g_xxxyyz_yy_0 = pbuffer.data(idx_eri_0_id + 45);

    auto g_xxxyyz_yz_0 = pbuffer.data(idx_eri_0_id + 46);

    auto g_xxxyyz_zz_0 = pbuffer.data(idx_eri_0_id + 47);

    #pragma omp simd aligned(g_xxxyy_x_1, g_xxxyy_xx_1, g_xxxyy_xy_1, g_xxxyy_xz_1, g_xxxyy_y_1, g_xxxyy_yy_1, g_xxxyy_yz_1, g_xxxyy_z_1, g_xxxyy_zz_1, g_xxxyyz_xx_0, g_xxxyyz_xy_0, g_xxxyyz_xz_0, g_xxxyyz_yy_0, g_xxxyyz_yz_0, g_xxxyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyz_xx_0[i] = g_xxxyy_xx_1[i] * pa_z[i];

        g_xxxyyz_xy_0[i] = g_xxxyy_xy_1[i] * pa_z[i];

        g_xxxyyz_xz_0[i] = g_xxxyy_x_1[i] * fe_0 + g_xxxyy_xz_1[i] * pa_z[i];

        g_xxxyyz_yy_0[i] = g_xxxyy_yy_1[i] * pa_z[i];

        g_xxxyyz_yz_0[i] = g_xxxyy_y_1[i] * fe_0 + g_xxxyy_yz_1[i] * pa_z[i];

        g_xxxyyz_zz_0[i] = 2.0 * g_xxxyy_z_1[i] * fe_0 + g_xxxyy_zz_1[i] * pa_z[i];
    }

    // Set up 48-54 components of targeted buffer : ID

    auto g_xxxyzz_xx_0 = pbuffer.data(idx_eri_0_id + 48);

    auto g_xxxyzz_xy_0 = pbuffer.data(idx_eri_0_id + 49);

    auto g_xxxyzz_xz_0 = pbuffer.data(idx_eri_0_id + 50);

    auto g_xxxyzz_yy_0 = pbuffer.data(idx_eri_0_id + 51);

    auto g_xxxyzz_yz_0 = pbuffer.data(idx_eri_0_id + 52);

    auto g_xxxyzz_zz_0 = pbuffer.data(idx_eri_0_id + 53);

    #pragma omp simd aligned(g_xxxyzz_xx_0, g_xxxyzz_xy_0, g_xxxyzz_xz_0, g_xxxyzz_yy_0, g_xxxyzz_yz_0, g_xxxyzz_zz_0, g_xxxzz_x_1, g_xxxzz_xx_1, g_xxxzz_xy_1, g_xxxzz_xz_1, g_xxxzz_y_1, g_xxxzz_yy_1, g_xxxzz_yz_1, g_xxxzz_z_1, g_xxxzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzz_xx_0[i] = g_xxxzz_xx_1[i] * pa_y[i];

        g_xxxyzz_xy_0[i] = g_xxxzz_x_1[i] * fe_0 + g_xxxzz_xy_1[i] * pa_y[i];

        g_xxxyzz_xz_0[i] = g_xxxzz_xz_1[i] * pa_y[i];

        g_xxxyzz_yy_0[i] = 2.0 * g_xxxzz_y_1[i] * fe_0 + g_xxxzz_yy_1[i] * pa_y[i];

        g_xxxyzz_yz_0[i] = g_xxxzz_z_1[i] * fe_0 + g_xxxzz_yz_1[i] * pa_y[i];

        g_xxxyzz_zz_0[i] = g_xxxzz_zz_1[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : ID

    auto g_xxxzzz_xx_0 = pbuffer.data(idx_eri_0_id + 54);

    auto g_xxxzzz_xy_0 = pbuffer.data(idx_eri_0_id + 55);

    auto g_xxxzzz_xz_0 = pbuffer.data(idx_eri_0_id + 56);

    auto g_xxxzzz_yy_0 = pbuffer.data(idx_eri_0_id + 57);

    auto g_xxxzzz_yz_0 = pbuffer.data(idx_eri_0_id + 58);

    auto g_xxxzzz_zz_0 = pbuffer.data(idx_eri_0_id + 59);

    #pragma omp simd aligned(g_xxxz_xx_0, g_xxxz_xx_1, g_xxxz_xy_0, g_xxxz_xy_1, g_xxxzz_xx_1, g_xxxzz_xy_1, g_xxxzzz_xx_0, g_xxxzzz_xy_0, g_xxxzzz_xz_0, g_xxxzzz_yy_0, g_xxxzzz_yz_0, g_xxxzzz_zz_0, g_xxzzz_xz_1, g_xxzzz_yy_1, g_xxzzz_yz_1, g_xxzzz_z_1, g_xxzzz_zz_1, g_xzzz_xz_0, g_xzzz_xz_1, g_xzzz_yy_0, g_xzzz_yy_1, g_xzzz_yz_0, g_xzzz_yz_1, g_xzzz_zz_0, g_xzzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzz_xx_0[i] = 2.0 * g_xxxz_xx_0[i] * fbe_0 - 2.0 * g_xxxz_xx_1[i] * fz_be_0 + g_xxxzz_xx_1[i] * pa_z[i];

        g_xxxzzz_xy_0[i] = 2.0 * g_xxxz_xy_0[i] * fbe_0 - 2.0 * g_xxxz_xy_1[i] * fz_be_0 + g_xxxzz_xy_1[i] * pa_z[i];

        g_xxxzzz_xz_0[i] = 2.0 * g_xzzz_xz_0[i] * fbe_0 - 2.0 * g_xzzz_xz_1[i] * fz_be_0 + g_xxzzz_z_1[i] * fe_0 + g_xxzzz_xz_1[i] * pa_x[i];

        g_xxxzzz_yy_0[i] = 2.0 * g_xzzz_yy_0[i] * fbe_0 - 2.0 * g_xzzz_yy_1[i] * fz_be_0 + g_xxzzz_yy_1[i] * pa_x[i];

        g_xxxzzz_yz_0[i] = 2.0 * g_xzzz_yz_0[i] * fbe_0 - 2.0 * g_xzzz_yz_1[i] * fz_be_0 + g_xxzzz_yz_1[i] * pa_x[i];

        g_xxxzzz_zz_0[i] = 2.0 * g_xzzz_zz_0[i] * fbe_0 - 2.0 * g_xzzz_zz_1[i] * fz_be_0 + g_xxzzz_zz_1[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : ID

    auto g_xxyyyy_xx_0 = pbuffer.data(idx_eri_0_id + 60);

    auto g_xxyyyy_xy_0 = pbuffer.data(idx_eri_0_id + 61);

    auto g_xxyyyy_xz_0 = pbuffer.data(idx_eri_0_id + 62);

    auto g_xxyyyy_yy_0 = pbuffer.data(idx_eri_0_id + 63);

    auto g_xxyyyy_yz_0 = pbuffer.data(idx_eri_0_id + 64);

    auto g_xxyyyy_zz_0 = pbuffer.data(idx_eri_0_id + 65);

    #pragma omp simd aligned(g_xxyy_xx_0, g_xxyy_xx_1, g_xxyy_xz_0, g_xxyy_xz_1, g_xxyyy_xx_1, g_xxyyy_xz_1, g_xxyyyy_xx_0, g_xxyyyy_xy_0, g_xxyyyy_xz_0, g_xxyyyy_yy_0, g_xxyyyy_yz_0, g_xxyyyy_zz_0, g_xyyyy_xy_1, g_xyyyy_y_1, g_xyyyy_yy_1, g_xyyyy_yz_1, g_xyyyy_zz_1, g_yyyy_xy_0, g_yyyy_xy_1, g_yyyy_yy_0, g_yyyy_yy_1, g_yyyy_yz_0, g_yyyy_yz_1, g_yyyy_zz_0, g_yyyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyy_xx_0[i] = 3.0 * g_xxyy_xx_0[i] * fbe_0 - 3.0 * g_xxyy_xx_1[i] * fz_be_0 + g_xxyyy_xx_1[i] * pa_y[i];

        g_xxyyyy_xy_0[i] = g_yyyy_xy_0[i] * fbe_0 - g_yyyy_xy_1[i] * fz_be_0 + g_xyyyy_y_1[i] * fe_0 + g_xyyyy_xy_1[i] * pa_x[i];

        g_xxyyyy_xz_0[i] = 3.0 * g_xxyy_xz_0[i] * fbe_0 - 3.0 * g_xxyy_xz_1[i] * fz_be_0 + g_xxyyy_xz_1[i] * pa_y[i];

        g_xxyyyy_yy_0[i] = g_yyyy_yy_0[i] * fbe_0 - g_yyyy_yy_1[i] * fz_be_0 + g_xyyyy_yy_1[i] * pa_x[i];

        g_xxyyyy_yz_0[i] = g_yyyy_yz_0[i] * fbe_0 - g_yyyy_yz_1[i] * fz_be_0 + g_xyyyy_yz_1[i] * pa_x[i];

        g_xxyyyy_zz_0[i] = g_yyyy_zz_0[i] * fbe_0 - g_yyyy_zz_1[i] * fz_be_0 + g_xyyyy_zz_1[i] * pa_x[i];
    }

    // Set up 66-72 components of targeted buffer : ID

    auto g_xxyyyz_xx_0 = pbuffer.data(idx_eri_0_id + 66);

    auto g_xxyyyz_xy_0 = pbuffer.data(idx_eri_0_id + 67);

    auto g_xxyyyz_xz_0 = pbuffer.data(idx_eri_0_id + 68);

    auto g_xxyyyz_yy_0 = pbuffer.data(idx_eri_0_id + 69);

    auto g_xxyyyz_yz_0 = pbuffer.data(idx_eri_0_id + 70);

    auto g_xxyyyz_zz_0 = pbuffer.data(idx_eri_0_id + 71);

    #pragma omp simd aligned(g_xxyyy_x_1, g_xxyyy_xx_1, g_xxyyy_xy_1, g_xxyyy_xz_1, g_xxyyy_y_1, g_xxyyy_yy_1, g_xxyyy_yz_1, g_xxyyy_z_1, g_xxyyy_zz_1, g_xxyyyz_xx_0, g_xxyyyz_xy_0, g_xxyyyz_xz_0, g_xxyyyz_yy_0, g_xxyyyz_yz_0, g_xxyyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyz_xx_0[i] = g_xxyyy_xx_1[i] * pa_z[i];

        g_xxyyyz_xy_0[i] = g_xxyyy_xy_1[i] * pa_z[i];

        g_xxyyyz_xz_0[i] = g_xxyyy_x_1[i] * fe_0 + g_xxyyy_xz_1[i] * pa_z[i];

        g_xxyyyz_yy_0[i] = g_xxyyy_yy_1[i] * pa_z[i];

        g_xxyyyz_yz_0[i] = g_xxyyy_y_1[i] * fe_0 + g_xxyyy_yz_1[i] * pa_z[i];

        g_xxyyyz_zz_0[i] = 2.0 * g_xxyyy_z_1[i] * fe_0 + g_xxyyy_zz_1[i] * pa_z[i];
    }

    // Set up 72-78 components of targeted buffer : ID

    auto g_xxyyzz_xx_0 = pbuffer.data(idx_eri_0_id + 72);

    auto g_xxyyzz_xy_0 = pbuffer.data(idx_eri_0_id + 73);

    auto g_xxyyzz_xz_0 = pbuffer.data(idx_eri_0_id + 74);

    auto g_xxyyzz_yy_0 = pbuffer.data(idx_eri_0_id + 75);

    auto g_xxyyzz_yz_0 = pbuffer.data(idx_eri_0_id + 76);

    auto g_xxyyzz_zz_0 = pbuffer.data(idx_eri_0_id + 77);

    #pragma omp simd aligned(g_xxyy_xy_0, g_xxyy_xy_1, g_xxyyz_xy_1, g_xxyyzz_xx_0, g_xxyyzz_xy_0, g_xxyyzz_xz_0, g_xxyyzz_yy_0, g_xxyyzz_yz_0, g_xxyyzz_zz_0, g_xxyzz_xx_1, g_xxyzz_xz_1, g_xxzz_xx_0, g_xxzz_xx_1, g_xxzz_xz_0, g_xxzz_xz_1, g_xyyzz_yy_1, g_xyyzz_yz_1, g_xyyzz_zz_1, g_yyzz_yy_0, g_yyzz_yy_1, g_yyzz_yz_0, g_yyzz_yz_1, g_yyzz_zz_0, g_yyzz_zz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyzz_xx_0[i] = g_xxzz_xx_0[i] * fbe_0 - g_xxzz_xx_1[i] * fz_be_0 + g_xxyzz_xx_1[i] * pa_y[i];

        g_xxyyzz_xy_0[i] = g_xxyy_xy_0[i] * fbe_0 - g_xxyy_xy_1[i] * fz_be_0 + g_xxyyz_xy_1[i] * pa_z[i];

        g_xxyyzz_xz_0[i] = g_xxzz_xz_0[i] * fbe_0 - g_xxzz_xz_1[i] * fz_be_0 + g_xxyzz_xz_1[i] * pa_y[i];

        g_xxyyzz_yy_0[i] = g_yyzz_yy_0[i] * fbe_0 - g_yyzz_yy_1[i] * fz_be_0 + g_xyyzz_yy_1[i] * pa_x[i];

        g_xxyyzz_yz_0[i] = g_yyzz_yz_0[i] * fbe_0 - g_yyzz_yz_1[i] * fz_be_0 + g_xyyzz_yz_1[i] * pa_x[i];

        g_xxyyzz_zz_0[i] = g_yyzz_zz_0[i] * fbe_0 - g_yyzz_zz_1[i] * fz_be_0 + g_xyyzz_zz_1[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : ID

    auto g_xxyzzz_xx_0 = pbuffer.data(idx_eri_0_id + 78);

    auto g_xxyzzz_xy_0 = pbuffer.data(idx_eri_0_id + 79);

    auto g_xxyzzz_xz_0 = pbuffer.data(idx_eri_0_id + 80);

    auto g_xxyzzz_yy_0 = pbuffer.data(idx_eri_0_id + 81);

    auto g_xxyzzz_yz_0 = pbuffer.data(idx_eri_0_id + 82);

    auto g_xxyzzz_zz_0 = pbuffer.data(idx_eri_0_id + 83);

    #pragma omp simd aligned(g_xxyzzz_xx_0, g_xxyzzz_xy_0, g_xxyzzz_xz_0, g_xxyzzz_yy_0, g_xxyzzz_yz_0, g_xxyzzz_zz_0, g_xxzzz_x_1, g_xxzzz_xx_1, g_xxzzz_xy_1, g_xxzzz_xz_1, g_xxzzz_y_1, g_xxzzz_yy_1, g_xxzzz_yz_1, g_xxzzz_z_1, g_xxzzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzz_xx_0[i] = g_xxzzz_xx_1[i] * pa_y[i];

        g_xxyzzz_xy_0[i] = g_xxzzz_x_1[i] * fe_0 + g_xxzzz_xy_1[i] * pa_y[i];

        g_xxyzzz_xz_0[i] = g_xxzzz_xz_1[i] * pa_y[i];

        g_xxyzzz_yy_0[i] = 2.0 * g_xxzzz_y_1[i] * fe_0 + g_xxzzz_yy_1[i] * pa_y[i];

        g_xxyzzz_yz_0[i] = g_xxzzz_z_1[i] * fe_0 + g_xxzzz_yz_1[i] * pa_y[i];

        g_xxyzzz_zz_0[i] = g_xxzzz_zz_1[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : ID

    auto g_xxzzzz_xx_0 = pbuffer.data(idx_eri_0_id + 84);

    auto g_xxzzzz_xy_0 = pbuffer.data(idx_eri_0_id + 85);

    auto g_xxzzzz_xz_0 = pbuffer.data(idx_eri_0_id + 86);

    auto g_xxzzzz_yy_0 = pbuffer.data(idx_eri_0_id + 87);

    auto g_xxzzzz_yz_0 = pbuffer.data(idx_eri_0_id + 88);

    auto g_xxzzzz_zz_0 = pbuffer.data(idx_eri_0_id + 89);

    #pragma omp simd aligned(g_xxzz_xx_0, g_xxzz_xx_1, g_xxzz_xy_0, g_xxzz_xy_1, g_xxzzz_xx_1, g_xxzzz_xy_1, g_xxzzzz_xx_0, g_xxzzzz_xy_0, g_xxzzzz_xz_0, g_xxzzzz_yy_0, g_xxzzzz_yz_0, g_xxzzzz_zz_0, g_xzzzz_xz_1, g_xzzzz_yy_1, g_xzzzz_yz_1, g_xzzzz_z_1, g_xzzzz_zz_1, g_zzzz_xz_0, g_zzzz_xz_1, g_zzzz_yy_0, g_zzzz_yy_1, g_zzzz_yz_0, g_zzzz_yz_1, g_zzzz_zz_0, g_zzzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzz_xx_0[i] = 3.0 * g_xxzz_xx_0[i] * fbe_0 - 3.0 * g_xxzz_xx_1[i] * fz_be_0 + g_xxzzz_xx_1[i] * pa_z[i];

        g_xxzzzz_xy_0[i] = 3.0 * g_xxzz_xy_0[i] * fbe_0 - 3.0 * g_xxzz_xy_1[i] * fz_be_0 + g_xxzzz_xy_1[i] * pa_z[i];

        g_xxzzzz_xz_0[i] = g_zzzz_xz_0[i] * fbe_0 - g_zzzz_xz_1[i] * fz_be_0 + g_xzzzz_z_1[i] * fe_0 + g_xzzzz_xz_1[i] * pa_x[i];

        g_xxzzzz_yy_0[i] = g_zzzz_yy_0[i] * fbe_0 - g_zzzz_yy_1[i] * fz_be_0 + g_xzzzz_yy_1[i] * pa_x[i];

        g_xxzzzz_yz_0[i] = g_zzzz_yz_0[i] * fbe_0 - g_zzzz_yz_1[i] * fz_be_0 + g_xzzzz_yz_1[i] * pa_x[i];

        g_xxzzzz_zz_0[i] = g_zzzz_zz_0[i] * fbe_0 - g_zzzz_zz_1[i] * fz_be_0 + g_xzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : ID

    auto g_xyyyyy_xx_0 = pbuffer.data(idx_eri_0_id + 90);

    auto g_xyyyyy_xy_0 = pbuffer.data(idx_eri_0_id + 91);

    auto g_xyyyyy_xz_0 = pbuffer.data(idx_eri_0_id + 92);

    auto g_xyyyyy_yy_0 = pbuffer.data(idx_eri_0_id + 93);

    auto g_xyyyyy_yz_0 = pbuffer.data(idx_eri_0_id + 94);

    auto g_xyyyyy_zz_0 = pbuffer.data(idx_eri_0_id + 95);

    #pragma omp simd aligned(g_xyyyyy_xx_0, g_xyyyyy_xy_0, g_xyyyyy_xz_0, g_xyyyyy_yy_0, g_xyyyyy_yz_0, g_xyyyyy_zz_0, g_yyyyy_x_1, g_yyyyy_xx_1, g_yyyyy_xy_1, g_yyyyy_xz_1, g_yyyyy_y_1, g_yyyyy_yy_1, g_yyyyy_yz_1, g_yyyyy_z_1, g_yyyyy_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyy_xx_0[i] = 2.0 * g_yyyyy_x_1[i] * fe_0 + g_yyyyy_xx_1[i] * pa_x[i];

        g_xyyyyy_xy_0[i] = g_yyyyy_y_1[i] * fe_0 + g_yyyyy_xy_1[i] * pa_x[i];

        g_xyyyyy_xz_0[i] = g_yyyyy_z_1[i] * fe_0 + g_yyyyy_xz_1[i] * pa_x[i];

        g_xyyyyy_yy_0[i] = g_yyyyy_yy_1[i] * pa_x[i];

        g_xyyyyy_yz_0[i] = g_yyyyy_yz_1[i] * pa_x[i];

        g_xyyyyy_zz_0[i] = g_yyyyy_zz_1[i] * pa_x[i];
    }

    // Set up 96-102 components of targeted buffer : ID

    auto g_xyyyyz_xx_0 = pbuffer.data(idx_eri_0_id + 96);

    auto g_xyyyyz_xy_0 = pbuffer.data(idx_eri_0_id + 97);

    auto g_xyyyyz_xz_0 = pbuffer.data(idx_eri_0_id + 98);

    auto g_xyyyyz_yy_0 = pbuffer.data(idx_eri_0_id + 99);

    auto g_xyyyyz_yz_0 = pbuffer.data(idx_eri_0_id + 100);

    auto g_xyyyyz_zz_0 = pbuffer.data(idx_eri_0_id + 101);

    #pragma omp simd aligned(g_xyyyy_xx_1, g_xyyyy_xy_1, g_xyyyyz_xx_0, g_xyyyyz_xy_0, g_xyyyyz_xz_0, g_xyyyyz_yy_0, g_xyyyyz_yz_0, g_xyyyyz_zz_0, g_yyyyz_xz_1, g_yyyyz_yy_1, g_yyyyz_yz_1, g_yyyyz_z_1, g_yyyyz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyz_xx_0[i] = g_xyyyy_xx_1[i] * pa_z[i];

        g_xyyyyz_xy_0[i] = g_xyyyy_xy_1[i] * pa_z[i];

        g_xyyyyz_xz_0[i] = g_yyyyz_z_1[i] * fe_0 + g_yyyyz_xz_1[i] * pa_x[i];

        g_xyyyyz_yy_0[i] = g_yyyyz_yy_1[i] * pa_x[i];

        g_xyyyyz_yz_0[i] = g_yyyyz_yz_1[i] * pa_x[i];

        g_xyyyyz_zz_0[i] = g_yyyyz_zz_1[i] * pa_x[i];
    }

    // Set up 102-108 components of targeted buffer : ID

    auto g_xyyyzz_xx_0 = pbuffer.data(idx_eri_0_id + 102);

    auto g_xyyyzz_xy_0 = pbuffer.data(idx_eri_0_id + 103);

    auto g_xyyyzz_xz_0 = pbuffer.data(idx_eri_0_id + 104);

    auto g_xyyyzz_yy_0 = pbuffer.data(idx_eri_0_id + 105);

    auto g_xyyyzz_yz_0 = pbuffer.data(idx_eri_0_id + 106);

    auto g_xyyyzz_zz_0 = pbuffer.data(idx_eri_0_id + 107);

    #pragma omp simd aligned(g_xyyyzz_xx_0, g_xyyyzz_xy_0, g_xyyyzz_xz_0, g_xyyyzz_yy_0, g_xyyyzz_yz_0, g_xyyyzz_zz_0, g_yyyzz_x_1, g_yyyzz_xx_1, g_yyyzz_xy_1, g_yyyzz_xz_1, g_yyyzz_y_1, g_yyyzz_yy_1, g_yyyzz_yz_1, g_yyyzz_z_1, g_yyyzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzz_xx_0[i] = 2.0 * g_yyyzz_x_1[i] * fe_0 + g_yyyzz_xx_1[i] * pa_x[i];

        g_xyyyzz_xy_0[i] = g_yyyzz_y_1[i] * fe_0 + g_yyyzz_xy_1[i] * pa_x[i];

        g_xyyyzz_xz_0[i] = g_yyyzz_z_1[i] * fe_0 + g_yyyzz_xz_1[i] * pa_x[i];

        g_xyyyzz_yy_0[i] = g_yyyzz_yy_1[i] * pa_x[i];

        g_xyyyzz_yz_0[i] = g_yyyzz_yz_1[i] * pa_x[i];

        g_xyyyzz_zz_0[i] = g_yyyzz_zz_1[i] * pa_x[i];
    }

    // Set up 108-114 components of targeted buffer : ID

    auto g_xyyzzz_xx_0 = pbuffer.data(idx_eri_0_id + 108);

    auto g_xyyzzz_xy_0 = pbuffer.data(idx_eri_0_id + 109);

    auto g_xyyzzz_xz_0 = pbuffer.data(idx_eri_0_id + 110);

    auto g_xyyzzz_yy_0 = pbuffer.data(idx_eri_0_id + 111);

    auto g_xyyzzz_yz_0 = pbuffer.data(idx_eri_0_id + 112);

    auto g_xyyzzz_zz_0 = pbuffer.data(idx_eri_0_id + 113);

    #pragma omp simd aligned(g_xyyzzz_xx_0, g_xyyzzz_xy_0, g_xyyzzz_xz_0, g_xyyzzz_yy_0, g_xyyzzz_yz_0, g_xyyzzz_zz_0, g_yyzzz_x_1, g_yyzzz_xx_1, g_yyzzz_xy_1, g_yyzzz_xz_1, g_yyzzz_y_1, g_yyzzz_yy_1, g_yyzzz_yz_1, g_yyzzz_z_1, g_yyzzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzz_xx_0[i] = 2.0 * g_yyzzz_x_1[i] * fe_0 + g_yyzzz_xx_1[i] * pa_x[i];

        g_xyyzzz_xy_0[i] = g_yyzzz_y_1[i] * fe_0 + g_yyzzz_xy_1[i] * pa_x[i];

        g_xyyzzz_xz_0[i] = g_yyzzz_z_1[i] * fe_0 + g_yyzzz_xz_1[i] * pa_x[i];

        g_xyyzzz_yy_0[i] = g_yyzzz_yy_1[i] * pa_x[i];

        g_xyyzzz_yz_0[i] = g_yyzzz_yz_1[i] * pa_x[i];

        g_xyyzzz_zz_0[i] = g_yyzzz_zz_1[i] * pa_x[i];
    }

    // Set up 114-120 components of targeted buffer : ID

    auto g_xyzzzz_xx_0 = pbuffer.data(idx_eri_0_id + 114);

    auto g_xyzzzz_xy_0 = pbuffer.data(idx_eri_0_id + 115);

    auto g_xyzzzz_xz_0 = pbuffer.data(idx_eri_0_id + 116);

    auto g_xyzzzz_yy_0 = pbuffer.data(idx_eri_0_id + 117);

    auto g_xyzzzz_yz_0 = pbuffer.data(idx_eri_0_id + 118);

    auto g_xyzzzz_zz_0 = pbuffer.data(idx_eri_0_id + 119);

    #pragma omp simd aligned(g_xyzzzz_xx_0, g_xyzzzz_xy_0, g_xyzzzz_xz_0, g_xyzzzz_yy_0, g_xyzzzz_yz_0, g_xyzzzz_zz_0, g_xzzzz_xx_1, g_xzzzz_xz_1, g_yzzzz_xy_1, g_yzzzz_y_1, g_yzzzz_yy_1, g_yzzzz_yz_1, g_yzzzz_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzz_xx_0[i] = g_xzzzz_xx_1[i] * pa_y[i];

        g_xyzzzz_xy_0[i] = g_yzzzz_y_1[i] * fe_0 + g_yzzzz_xy_1[i] * pa_x[i];

        g_xyzzzz_xz_0[i] = g_xzzzz_xz_1[i] * pa_y[i];

        g_xyzzzz_yy_0[i] = g_yzzzz_yy_1[i] * pa_x[i];

        g_xyzzzz_yz_0[i] = g_yzzzz_yz_1[i] * pa_x[i];

        g_xyzzzz_zz_0[i] = g_yzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 120-126 components of targeted buffer : ID

    auto g_xzzzzz_xx_0 = pbuffer.data(idx_eri_0_id + 120);

    auto g_xzzzzz_xy_0 = pbuffer.data(idx_eri_0_id + 121);

    auto g_xzzzzz_xz_0 = pbuffer.data(idx_eri_0_id + 122);

    auto g_xzzzzz_yy_0 = pbuffer.data(idx_eri_0_id + 123);

    auto g_xzzzzz_yz_0 = pbuffer.data(idx_eri_0_id + 124);

    auto g_xzzzzz_zz_0 = pbuffer.data(idx_eri_0_id + 125);

    #pragma omp simd aligned(g_xzzzzz_xx_0, g_xzzzzz_xy_0, g_xzzzzz_xz_0, g_xzzzzz_yy_0, g_xzzzzz_yz_0, g_xzzzzz_zz_0, g_zzzzz_x_1, g_zzzzz_xx_1, g_zzzzz_xy_1, g_zzzzz_xz_1, g_zzzzz_y_1, g_zzzzz_yy_1, g_zzzzz_yz_1, g_zzzzz_z_1, g_zzzzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzz_xx_0[i] = 2.0 * g_zzzzz_x_1[i] * fe_0 + g_zzzzz_xx_1[i] * pa_x[i];

        g_xzzzzz_xy_0[i] = g_zzzzz_y_1[i] * fe_0 + g_zzzzz_xy_1[i] * pa_x[i];

        g_xzzzzz_xz_0[i] = g_zzzzz_z_1[i] * fe_0 + g_zzzzz_xz_1[i] * pa_x[i];

        g_xzzzzz_yy_0[i] = g_zzzzz_yy_1[i] * pa_x[i];

        g_xzzzzz_yz_0[i] = g_zzzzz_yz_1[i] * pa_x[i];

        g_xzzzzz_zz_0[i] = g_zzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 126-132 components of targeted buffer : ID

    auto g_yyyyyy_xx_0 = pbuffer.data(idx_eri_0_id + 126);

    auto g_yyyyyy_xy_0 = pbuffer.data(idx_eri_0_id + 127);

    auto g_yyyyyy_xz_0 = pbuffer.data(idx_eri_0_id + 128);

    auto g_yyyyyy_yy_0 = pbuffer.data(idx_eri_0_id + 129);

    auto g_yyyyyy_yz_0 = pbuffer.data(idx_eri_0_id + 130);

    auto g_yyyyyy_zz_0 = pbuffer.data(idx_eri_0_id + 131);

    #pragma omp simd aligned(g_yyyy_xx_0, g_yyyy_xx_1, g_yyyy_xy_0, g_yyyy_xy_1, g_yyyy_xz_0, g_yyyy_xz_1, g_yyyy_yy_0, g_yyyy_yy_1, g_yyyy_yz_0, g_yyyy_yz_1, g_yyyy_zz_0, g_yyyy_zz_1, g_yyyyy_x_1, g_yyyyy_xx_1, g_yyyyy_xy_1, g_yyyyy_xz_1, g_yyyyy_y_1, g_yyyyy_yy_1, g_yyyyy_yz_1, g_yyyyy_z_1, g_yyyyy_zz_1, g_yyyyyy_xx_0, g_yyyyyy_xy_0, g_yyyyyy_xz_0, g_yyyyyy_yy_0, g_yyyyyy_yz_0, g_yyyyyy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyy_xx_0[i] = 5.0 * g_yyyy_xx_0[i] * fbe_0 - 5.0 * g_yyyy_xx_1[i] * fz_be_0 + g_yyyyy_xx_1[i] * pa_y[i];

        g_yyyyyy_xy_0[i] = 5.0 * g_yyyy_xy_0[i] * fbe_0 - 5.0 * g_yyyy_xy_1[i] * fz_be_0 + g_yyyyy_x_1[i] * fe_0 + g_yyyyy_xy_1[i] * pa_y[i];

        g_yyyyyy_xz_0[i] = 5.0 * g_yyyy_xz_0[i] * fbe_0 - 5.0 * g_yyyy_xz_1[i] * fz_be_0 + g_yyyyy_xz_1[i] * pa_y[i];

        g_yyyyyy_yy_0[i] = 5.0 * g_yyyy_yy_0[i] * fbe_0 - 5.0 * g_yyyy_yy_1[i] * fz_be_0 + 2.0 * g_yyyyy_y_1[i] * fe_0 + g_yyyyy_yy_1[i] * pa_y[i];

        g_yyyyyy_yz_0[i] = 5.0 * g_yyyy_yz_0[i] * fbe_0 - 5.0 * g_yyyy_yz_1[i] * fz_be_0 + g_yyyyy_z_1[i] * fe_0 + g_yyyyy_yz_1[i] * pa_y[i];

        g_yyyyyy_zz_0[i] = 5.0 * g_yyyy_zz_0[i] * fbe_0 - 5.0 * g_yyyy_zz_1[i] * fz_be_0 + g_yyyyy_zz_1[i] * pa_y[i];
    }

    // Set up 132-138 components of targeted buffer : ID

    auto g_yyyyyz_xx_0 = pbuffer.data(idx_eri_0_id + 132);

    auto g_yyyyyz_xy_0 = pbuffer.data(idx_eri_0_id + 133);

    auto g_yyyyyz_xz_0 = pbuffer.data(idx_eri_0_id + 134);

    auto g_yyyyyz_yy_0 = pbuffer.data(idx_eri_0_id + 135);

    auto g_yyyyyz_yz_0 = pbuffer.data(idx_eri_0_id + 136);

    auto g_yyyyyz_zz_0 = pbuffer.data(idx_eri_0_id + 137);

    #pragma omp simd aligned(g_yyyyy_x_1, g_yyyyy_xx_1, g_yyyyy_xy_1, g_yyyyy_xz_1, g_yyyyy_y_1, g_yyyyy_yy_1, g_yyyyy_yz_1, g_yyyyy_z_1, g_yyyyy_zz_1, g_yyyyyz_xx_0, g_yyyyyz_xy_0, g_yyyyyz_xz_0, g_yyyyyz_yy_0, g_yyyyyz_yz_0, g_yyyyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyz_xx_0[i] = g_yyyyy_xx_1[i] * pa_z[i];

        g_yyyyyz_xy_0[i] = g_yyyyy_xy_1[i] * pa_z[i];

        g_yyyyyz_xz_0[i] = g_yyyyy_x_1[i] * fe_0 + g_yyyyy_xz_1[i] * pa_z[i];

        g_yyyyyz_yy_0[i] = g_yyyyy_yy_1[i] * pa_z[i];

        g_yyyyyz_yz_0[i] = g_yyyyy_y_1[i] * fe_0 + g_yyyyy_yz_1[i] * pa_z[i];

        g_yyyyyz_zz_0[i] = 2.0 * g_yyyyy_z_1[i] * fe_0 + g_yyyyy_zz_1[i] * pa_z[i];
    }

    // Set up 138-144 components of targeted buffer : ID

    auto g_yyyyzz_xx_0 = pbuffer.data(idx_eri_0_id + 138);

    auto g_yyyyzz_xy_0 = pbuffer.data(idx_eri_0_id + 139);

    auto g_yyyyzz_xz_0 = pbuffer.data(idx_eri_0_id + 140);

    auto g_yyyyzz_yy_0 = pbuffer.data(idx_eri_0_id + 141);

    auto g_yyyyzz_yz_0 = pbuffer.data(idx_eri_0_id + 142);

    auto g_yyyyzz_zz_0 = pbuffer.data(idx_eri_0_id + 143);

    #pragma omp simd aligned(g_yyyy_xy_0, g_yyyy_xy_1, g_yyyy_yy_0, g_yyyy_yy_1, g_yyyyz_xy_1, g_yyyyz_yy_1, g_yyyyzz_xx_0, g_yyyyzz_xy_0, g_yyyyzz_xz_0, g_yyyyzz_yy_0, g_yyyyzz_yz_0, g_yyyyzz_zz_0, g_yyyzz_xx_1, g_yyyzz_xz_1, g_yyyzz_yz_1, g_yyyzz_z_1, g_yyyzz_zz_1, g_yyzz_xx_0, g_yyzz_xx_1, g_yyzz_xz_0, g_yyzz_xz_1, g_yyzz_yz_0, g_yyzz_yz_1, g_yyzz_zz_0, g_yyzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzz_xx_0[i] = 3.0 * g_yyzz_xx_0[i] * fbe_0 - 3.0 * g_yyzz_xx_1[i] * fz_be_0 + g_yyyzz_xx_1[i] * pa_y[i];

        g_yyyyzz_xy_0[i] = g_yyyy_xy_0[i] * fbe_0 - g_yyyy_xy_1[i] * fz_be_0 + g_yyyyz_xy_1[i] * pa_z[i];

        g_yyyyzz_xz_0[i] = 3.0 * g_yyzz_xz_0[i] * fbe_0 - 3.0 * g_yyzz_xz_1[i] * fz_be_0 + g_yyyzz_xz_1[i] * pa_y[i];

        g_yyyyzz_yy_0[i] = g_yyyy_yy_0[i] * fbe_0 - g_yyyy_yy_1[i] * fz_be_0 + g_yyyyz_yy_1[i] * pa_z[i];

        g_yyyyzz_yz_0[i] = 3.0 * g_yyzz_yz_0[i] * fbe_0 - 3.0 * g_yyzz_yz_1[i] * fz_be_0 + g_yyyzz_z_1[i] * fe_0 + g_yyyzz_yz_1[i] * pa_y[i];

        g_yyyyzz_zz_0[i] = 3.0 * g_yyzz_zz_0[i] * fbe_0 - 3.0 * g_yyzz_zz_1[i] * fz_be_0 + g_yyyzz_zz_1[i] * pa_y[i];
    }

    // Set up 144-150 components of targeted buffer : ID

    auto g_yyyzzz_xx_0 = pbuffer.data(idx_eri_0_id + 144);

    auto g_yyyzzz_xy_0 = pbuffer.data(idx_eri_0_id + 145);

    auto g_yyyzzz_xz_0 = pbuffer.data(idx_eri_0_id + 146);

    auto g_yyyzzz_yy_0 = pbuffer.data(idx_eri_0_id + 147);

    auto g_yyyzzz_yz_0 = pbuffer.data(idx_eri_0_id + 148);

    auto g_yyyzzz_zz_0 = pbuffer.data(idx_eri_0_id + 149);

    #pragma omp simd aligned(g_yyyz_xy_0, g_yyyz_xy_1, g_yyyz_yy_0, g_yyyz_yy_1, g_yyyzz_xy_1, g_yyyzz_yy_1, g_yyyzzz_xx_0, g_yyyzzz_xy_0, g_yyyzzz_xz_0, g_yyyzzz_yy_0, g_yyyzzz_yz_0, g_yyyzzz_zz_0, g_yyzzz_xx_1, g_yyzzz_xz_1, g_yyzzz_yz_1, g_yyzzz_z_1, g_yyzzz_zz_1, g_yzzz_xx_0, g_yzzz_xx_1, g_yzzz_xz_0, g_yzzz_xz_1, g_yzzz_yz_0, g_yzzz_yz_1, g_yzzz_zz_0, g_yzzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzz_xx_0[i] = 2.0 * g_yzzz_xx_0[i] * fbe_0 - 2.0 * g_yzzz_xx_1[i] * fz_be_0 + g_yyzzz_xx_1[i] * pa_y[i];

        g_yyyzzz_xy_0[i] = 2.0 * g_yyyz_xy_0[i] * fbe_0 - 2.0 * g_yyyz_xy_1[i] * fz_be_0 + g_yyyzz_xy_1[i] * pa_z[i];

        g_yyyzzz_xz_0[i] = 2.0 * g_yzzz_xz_0[i] * fbe_0 - 2.0 * g_yzzz_xz_1[i] * fz_be_0 + g_yyzzz_xz_1[i] * pa_y[i];

        g_yyyzzz_yy_0[i] = 2.0 * g_yyyz_yy_0[i] * fbe_0 - 2.0 * g_yyyz_yy_1[i] * fz_be_0 + g_yyyzz_yy_1[i] * pa_z[i];

        g_yyyzzz_yz_0[i] = 2.0 * g_yzzz_yz_0[i] * fbe_0 - 2.0 * g_yzzz_yz_1[i] * fz_be_0 + g_yyzzz_z_1[i] * fe_0 + g_yyzzz_yz_1[i] * pa_y[i];

        g_yyyzzz_zz_0[i] = 2.0 * g_yzzz_zz_0[i] * fbe_0 - 2.0 * g_yzzz_zz_1[i] * fz_be_0 + g_yyzzz_zz_1[i] * pa_y[i];
    }

    // Set up 150-156 components of targeted buffer : ID

    auto g_yyzzzz_xx_0 = pbuffer.data(idx_eri_0_id + 150);

    auto g_yyzzzz_xy_0 = pbuffer.data(idx_eri_0_id + 151);

    auto g_yyzzzz_xz_0 = pbuffer.data(idx_eri_0_id + 152);

    auto g_yyzzzz_yy_0 = pbuffer.data(idx_eri_0_id + 153);

    auto g_yyzzzz_yz_0 = pbuffer.data(idx_eri_0_id + 154);

    auto g_yyzzzz_zz_0 = pbuffer.data(idx_eri_0_id + 155);

    #pragma omp simd aligned(g_yyzz_xy_0, g_yyzz_xy_1, g_yyzz_yy_0, g_yyzz_yy_1, g_yyzzz_xy_1, g_yyzzz_yy_1, g_yyzzzz_xx_0, g_yyzzzz_xy_0, g_yyzzzz_xz_0, g_yyzzzz_yy_0, g_yyzzzz_yz_0, g_yyzzzz_zz_0, g_yzzzz_xx_1, g_yzzzz_xz_1, g_yzzzz_yz_1, g_yzzzz_z_1, g_yzzzz_zz_1, g_zzzz_xx_0, g_zzzz_xx_1, g_zzzz_xz_0, g_zzzz_xz_1, g_zzzz_yz_0, g_zzzz_yz_1, g_zzzz_zz_0, g_zzzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzz_xx_0[i] = g_zzzz_xx_0[i] * fbe_0 - g_zzzz_xx_1[i] * fz_be_0 + g_yzzzz_xx_1[i] * pa_y[i];

        g_yyzzzz_xy_0[i] = 3.0 * g_yyzz_xy_0[i] * fbe_0 - 3.0 * g_yyzz_xy_1[i] * fz_be_0 + g_yyzzz_xy_1[i] * pa_z[i];

        g_yyzzzz_xz_0[i] = g_zzzz_xz_0[i] * fbe_0 - g_zzzz_xz_1[i] * fz_be_0 + g_yzzzz_xz_1[i] * pa_y[i];

        g_yyzzzz_yy_0[i] = 3.0 * g_yyzz_yy_0[i] * fbe_0 - 3.0 * g_yyzz_yy_1[i] * fz_be_0 + g_yyzzz_yy_1[i] * pa_z[i];

        g_yyzzzz_yz_0[i] = g_zzzz_yz_0[i] * fbe_0 - g_zzzz_yz_1[i] * fz_be_0 + g_yzzzz_z_1[i] * fe_0 + g_yzzzz_yz_1[i] * pa_y[i];

        g_yyzzzz_zz_0[i] = g_zzzz_zz_0[i] * fbe_0 - g_zzzz_zz_1[i] * fz_be_0 + g_yzzzz_zz_1[i] * pa_y[i];
    }

    // Set up 156-162 components of targeted buffer : ID

    auto g_yzzzzz_xx_0 = pbuffer.data(idx_eri_0_id + 156);

    auto g_yzzzzz_xy_0 = pbuffer.data(idx_eri_0_id + 157);

    auto g_yzzzzz_xz_0 = pbuffer.data(idx_eri_0_id + 158);

    auto g_yzzzzz_yy_0 = pbuffer.data(idx_eri_0_id + 159);

    auto g_yzzzzz_yz_0 = pbuffer.data(idx_eri_0_id + 160);

    auto g_yzzzzz_zz_0 = pbuffer.data(idx_eri_0_id + 161);

    #pragma omp simd aligned(g_yzzzzz_xx_0, g_yzzzzz_xy_0, g_yzzzzz_xz_0, g_yzzzzz_yy_0, g_yzzzzz_yz_0, g_yzzzzz_zz_0, g_zzzzz_x_1, g_zzzzz_xx_1, g_zzzzz_xy_1, g_zzzzz_xz_1, g_zzzzz_y_1, g_zzzzz_yy_1, g_zzzzz_yz_1, g_zzzzz_z_1, g_zzzzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzz_xx_0[i] = g_zzzzz_xx_1[i] * pa_y[i];

        g_yzzzzz_xy_0[i] = g_zzzzz_x_1[i] * fe_0 + g_zzzzz_xy_1[i] * pa_y[i];

        g_yzzzzz_xz_0[i] = g_zzzzz_xz_1[i] * pa_y[i];

        g_yzzzzz_yy_0[i] = 2.0 * g_zzzzz_y_1[i] * fe_0 + g_zzzzz_yy_1[i] * pa_y[i];

        g_yzzzzz_yz_0[i] = g_zzzzz_z_1[i] * fe_0 + g_zzzzz_yz_1[i] * pa_y[i];

        g_yzzzzz_zz_0[i] = g_zzzzz_zz_1[i] * pa_y[i];
    }

    // Set up 162-168 components of targeted buffer : ID

    auto g_zzzzzz_xx_0 = pbuffer.data(idx_eri_0_id + 162);

    auto g_zzzzzz_xy_0 = pbuffer.data(idx_eri_0_id + 163);

    auto g_zzzzzz_xz_0 = pbuffer.data(idx_eri_0_id + 164);

    auto g_zzzzzz_yy_0 = pbuffer.data(idx_eri_0_id + 165);

    auto g_zzzzzz_yz_0 = pbuffer.data(idx_eri_0_id + 166);

    auto g_zzzzzz_zz_0 = pbuffer.data(idx_eri_0_id + 167);

    #pragma omp simd aligned(g_zzzz_xx_0, g_zzzz_xx_1, g_zzzz_xy_0, g_zzzz_xy_1, g_zzzz_xz_0, g_zzzz_xz_1, g_zzzz_yy_0, g_zzzz_yy_1, g_zzzz_yz_0, g_zzzz_yz_1, g_zzzz_zz_0, g_zzzz_zz_1, g_zzzzz_x_1, g_zzzzz_xx_1, g_zzzzz_xy_1, g_zzzzz_xz_1, g_zzzzz_y_1, g_zzzzz_yy_1, g_zzzzz_yz_1, g_zzzzz_z_1, g_zzzzz_zz_1, g_zzzzzz_xx_0, g_zzzzzz_xy_0, g_zzzzzz_xz_0, g_zzzzzz_yy_0, g_zzzzzz_yz_0, g_zzzzzz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzz_xx_0[i] = 5.0 * g_zzzz_xx_0[i] * fbe_0 - 5.0 * g_zzzz_xx_1[i] * fz_be_0 + g_zzzzz_xx_1[i] * pa_z[i];

        g_zzzzzz_xy_0[i] = 5.0 * g_zzzz_xy_0[i] * fbe_0 - 5.0 * g_zzzz_xy_1[i] * fz_be_0 + g_zzzzz_xy_1[i] * pa_z[i];

        g_zzzzzz_xz_0[i] = 5.0 * g_zzzz_xz_0[i] * fbe_0 - 5.0 * g_zzzz_xz_1[i] * fz_be_0 + g_zzzzz_x_1[i] * fe_0 + g_zzzzz_xz_1[i] * pa_z[i];

        g_zzzzzz_yy_0[i] = 5.0 * g_zzzz_yy_0[i] * fbe_0 - 5.0 * g_zzzz_yy_1[i] * fz_be_0 + g_zzzzz_yy_1[i] * pa_z[i];

        g_zzzzzz_yz_0[i] = 5.0 * g_zzzz_yz_0[i] * fbe_0 - 5.0 * g_zzzz_yz_1[i] * fz_be_0 + g_zzzzz_y_1[i] * fe_0 + g_zzzzz_yz_1[i] * pa_z[i];

        g_zzzzzz_zz_0[i] = 5.0 * g_zzzz_zz_0[i] * fbe_0 - 5.0 * g_zzzz_zz_1[i] * fz_be_0 + 2.0 * g_zzzzz_z_1[i] * fe_0 + g_zzzzz_zz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

