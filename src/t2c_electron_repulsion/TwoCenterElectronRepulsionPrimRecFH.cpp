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

#include "TwoCenterElectronRepulsionPrimRecFH.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_fh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fh,
                                const size_t idx_eri_0_ph,
                                const size_t idx_eri_1_ph,
                                const size_t idx_eri_1_dg,
                                const size_t idx_eri_1_dh,
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

    // Set up components of auxiliary buffer : PH

    auto g_x_xxxxx_0 = pbuffer.data(idx_eri_0_ph);

    auto g_x_xxxxy_0 = pbuffer.data(idx_eri_0_ph + 1);

    auto g_x_xxxxz_0 = pbuffer.data(idx_eri_0_ph + 2);

    auto g_x_xxxyy_0 = pbuffer.data(idx_eri_0_ph + 3);

    auto g_x_xxxyz_0 = pbuffer.data(idx_eri_0_ph + 4);

    auto g_x_xxxzz_0 = pbuffer.data(idx_eri_0_ph + 5);

    auto g_x_xxyyy_0 = pbuffer.data(idx_eri_0_ph + 6);

    auto g_x_xxyyz_0 = pbuffer.data(idx_eri_0_ph + 7);

    auto g_x_xxyzz_0 = pbuffer.data(idx_eri_0_ph + 8);

    auto g_x_xxzzz_0 = pbuffer.data(idx_eri_0_ph + 9);

    auto g_x_xyyyy_0 = pbuffer.data(idx_eri_0_ph + 10);

    auto g_x_xyyyz_0 = pbuffer.data(idx_eri_0_ph + 11);

    auto g_x_xyyzz_0 = pbuffer.data(idx_eri_0_ph + 12);

    auto g_x_xyzzz_0 = pbuffer.data(idx_eri_0_ph + 13);

    auto g_x_xzzzz_0 = pbuffer.data(idx_eri_0_ph + 14);

    auto g_x_yyyyy_0 = pbuffer.data(idx_eri_0_ph + 15);

    auto g_x_yyyyz_0 = pbuffer.data(idx_eri_0_ph + 16);

    auto g_x_yyyzz_0 = pbuffer.data(idx_eri_0_ph + 17);

    auto g_x_yyzzz_0 = pbuffer.data(idx_eri_0_ph + 18);

    auto g_x_yzzzz_0 = pbuffer.data(idx_eri_0_ph + 19);

    auto g_x_zzzzz_0 = pbuffer.data(idx_eri_0_ph + 20);

    auto g_y_xxxxx_0 = pbuffer.data(idx_eri_0_ph + 21);

    auto g_y_xxxxy_0 = pbuffer.data(idx_eri_0_ph + 22);

    auto g_y_xxxxz_0 = pbuffer.data(idx_eri_0_ph + 23);

    auto g_y_xxxyy_0 = pbuffer.data(idx_eri_0_ph + 24);

    auto g_y_xxxyz_0 = pbuffer.data(idx_eri_0_ph + 25);

    auto g_y_xxxzz_0 = pbuffer.data(idx_eri_0_ph + 26);

    auto g_y_xxyyy_0 = pbuffer.data(idx_eri_0_ph + 27);

    auto g_y_xxyyz_0 = pbuffer.data(idx_eri_0_ph + 28);

    auto g_y_xxyzz_0 = pbuffer.data(idx_eri_0_ph + 29);

    auto g_y_xxzzz_0 = pbuffer.data(idx_eri_0_ph + 30);

    auto g_y_xyyyy_0 = pbuffer.data(idx_eri_0_ph + 31);

    auto g_y_xyyyz_0 = pbuffer.data(idx_eri_0_ph + 32);

    auto g_y_xyyzz_0 = pbuffer.data(idx_eri_0_ph + 33);

    auto g_y_xyzzz_0 = pbuffer.data(idx_eri_0_ph + 34);

    auto g_y_xzzzz_0 = pbuffer.data(idx_eri_0_ph + 35);

    auto g_y_yyyyy_0 = pbuffer.data(idx_eri_0_ph + 36);

    auto g_y_yyyyz_0 = pbuffer.data(idx_eri_0_ph + 37);

    auto g_y_yyyzz_0 = pbuffer.data(idx_eri_0_ph + 38);

    auto g_y_yyzzz_0 = pbuffer.data(idx_eri_0_ph + 39);

    auto g_y_yzzzz_0 = pbuffer.data(idx_eri_0_ph + 40);

    auto g_y_zzzzz_0 = pbuffer.data(idx_eri_0_ph + 41);

    auto g_z_xxxxx_0 = pbuffer.data(idx_eri_0_ph + 42);

    auto g_z_xxxxy_0 = pbuffer.data(idx_eri_0_ph + 43);

    auto g_z_xxxxz_0 = pbuffer.data(idx_eri_0_ph + 44);

    auto g_z_xxxyy_0 = pbuffer.data(idx_eri_0_ph + 45);

    auto g_z_xxxyz_0 = pbuffer.data(idx_eri_0_ph + 46);

    auto g_z_xxxzz_0 = pbuffer.data(idx_eri_0_ph + 47);

    auto g_z_xxyyy_0 = pbuffer.data(idx_eri_0_ph + 48);

    auto g_z_xxyyz_0 = pbuffer.data(idx_eri_0_ph + 49);

    auto g_z_xxyzz_0 = pbuffer.data(idx_eri_0_ph + 50);

    auto g_z_xxzzz_0 = pbuffer.data(idx_eri_0_ph + 51);

    auto g_z_xyyyy_0 = pbuffer.data(idx_eri_0_ph + 52);

    auto g_z_xyyyz_0 = pbuffer.data(idx_eri_0_ph + 53);

    auto g_z_xyyzz_0 = pbuffer.data(idx_eri_0_ph + 54);

    auto g_z_xyzzz_0 = pbuffer.data(idx_eri_0_ph + 55);

    auto g_z_xzzzz_0 = pbuffer.data(idx_eri_0_ph + 56);

    auto g_z_yyyyy_0 = pbuffer.data(idx_eri_0_ph + 57);

    auto g_z_yyyyz_0 = pbuffer.data(idx_eri_0_ph + 58);

    auto g_z_yyyzz_0 = pbuffer.data(idx_eri_0_ph + 59);

    auto g_z_yyzzz_0 = pbuffer.data(idx_eri_0_ph + 60);

    auto g_z_yzzzz_0 = pbuffer.data(idx_eri_0_ph + 61);

    auto g_z_zzzzz_0 = pbuffer.data(idx_eri_0_ph + 62);

    // Set up components of auxiliary buffer : PH

    auto g_x_xxxxx_1 = pbuffer.data(idx_eri_1_ph);

    auto g_x_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 1);

    auto g_x_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 2);

    auto g_x_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 3);

    auto g_x_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 4);

    auto g_x_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 5);

    auto g_x_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 6);

    auto g_x_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 7);

    auto g_x_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 8);

    auto g_x_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 9);

    auto g_x_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 10);

    auto g_x_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 11);

    auto g_x_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 12);

    auto g_x_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 13);

    auto g_x_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 14);

    auto g_x_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 15);

    auto g_x_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 16);

    auto g_x_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 17);

    auto g_x_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 18);

    auto g_x_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 19);

    auto g_x_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 20);

    auto g_y_xxxxx_1 = pbuffer.data(idx_eri_1_ph + 21);

    auto g_y_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 22);

    auto g_y_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 23);

    auto g_y_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 24);

    auto g_y_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 25);

    auto g_y_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 26);

    auto g_y_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 27);

    auto g_y_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 28);

    auto g_y_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 29);

    auto g_y_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 30);

    auto g_y_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 31);

    auto g_y_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 32);

    auto g_y_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 33);

    auto g_y_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 34);

    auto g_y_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 35);

    auto g_y_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 36);

    auto g_y_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 37);

    auto g_y_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 38);

    auto g_y_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 39);

    auto g_y_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 40);

    auto g_y_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 41);

    auto g_z_xxxxx_1 = pbuffer.data(idx_eri_1_ph + 42);

    auto g_z_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 43);

    auto g_z_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 44);

    auto g_z_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 45);

    auto g_z_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 46);

    auto g_z_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 47);

    auto g_z_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 48);

    auto g_z_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 49);

    auto g_z_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 50);

    auto g_z_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 51);

    auto g_z_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 52);

    auto g_z_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 53);

    auto g_z_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 54);

    auto g_z_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 55);

    auto g_z_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 56);

    auto g_z_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 57);

    auto g_z_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 58);

    auto g_z_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 59);

    auto g_z_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 60);

    auto g_z_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 61);

    auto g_z_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 62);

    // Set up components of auxiliary buffer : DG

    auto g_xx_xxxx_1 = pbuffer.data(idx_eri_1_dg);

    auto g_xx_xxxy_1 = pbuffer.data(idx_eri_1_dg + 1);

    auto g_xx_xxxz_1 = pbuffer.data(idx_eri_1_dg + 2);

    auto g_xx_xxyy_1 = pbuffer.data(idx_eri_1_dg + 3);

    auto g_xx_xxyz_1 = pbuffer.data(idx_eri_1_dg + 4);

    auto g_xx_xxzz_1 = pbuffer.data(idx_eri_1_dg + 5);

    auto g_xx_xyyy_1 = pbuffer.data(idx_eri_1_dg + 6);

    auto g_xx_xyyz_1 = pbuffer.data(idx_eri_1_dg + 7);

    auto g_xx_xyzz_1 = pbuffer.data(idx_eri_1_dg + 8);

    auto g_xx_xzzz_1 = pbuffer.data(idx_eri_1_dg + 9);

    auto g_xx_yyyy_1 = pbuffer.data(idx_eri_1_dg + 10);

    auto g_xx_yyyz_1 = pbuffer.data(idx_eri_1_dg + 11);

    auto g_xx_yyzz_1 = pbuffer.data(idx_eri_1_dg + 12);

    auto g_xx_yzzz_1 = pbuffer.data(idx_eri_1_dg + 13);

    auto g_xx_zzzz_1 = pbuffer.data(idx_eri_1_dg + 14);

    auto g_yy_xxxx_1 = pbuffer.data(idx_eri_1_dg + 45);

    auto g_yy_xxxy_1 = pbuffer.data(idx_eri_1_dg + 46);

    auto g_yy_xxxz_1 = pbuffer.data(idx_eri_1_dg + 47);

    auto g_yy_xxyy_1 = pbuffer.data(idx_eri_1_dg + 48);

    auto g_yy_xxyz_1 = pbuffer.data(idx_eri_1_dg + 49);

    auto g_yy_xxzz_1 = pbuffer.data(idx_eri_1_dg + 50);

    auto g_yy_xyyy_1 = pbuffer.data(idx_eri_1_dg + 51);

    auto g_yy_xyyz_1 = pbuffer.data(idx_eri_1_dg + 52);

    auto g_yy_xyzz_1 = pbuffer.data(idx_eri_1_dg + 53);

    auto g_yy_xzzz_1 = pbuffer.data(idx_eri_1_dg + 54);

    auto g_yy_yyyy_1 = pbuffer.data(idx_eri_1_dg + 55);

    auto g_yy_yyyz_1 = pbuffer.data(idx_eri_1_dg + 56);

    auto g_yy_yyzz_1 = pbuffer.data(idx_eri_1_dg + 57);

    auto g_yy_yzzz_1 = pbuffer.data(idx_eri_1_dg + 58);

    auto g_yy_zzzz_1 = pbuffer.data(idx_eri_1_dg + 59);

    auto g_yz_xxyz_1 = pbuffer.data(idx_eri_1_dg + 64);

    auto g_yz_xyyz_1 = pbuffer.data(idx_eri_1_dg + 67);

    auto g_yz_xyzz_1 = pbuffer.data(idx_eri_1_dg + 68);

    auto g_yz_yyyz_1 = pbuffer.data(idx_eri_1_dg + 71);

    auto g_yz_yyzz_1 = pbuffer.data(idx_eri_1_dg + 72);

    auto g_yz_yzzz_1 = pbuffer.data(idx_eri_1_dg + 73);

    auto g_zz_xxxx_1 = pbuffer.data(idx_eri_1_dg + 75);

    auto g_zz_xxxy_1 = pbuffer.data(idx_eri_1_dg + 76);

    auto g_zz_xxxz_1 = pbuffer.data(idx_eri_1_dg + 77);

    auto g_zz_xxyy_1 = pbuffer.data(idx_eri_1_dg + 78);

    auto g_zz_xxyz_1 = pbuffer.data(idx_eri_1_dg + 79);

    auto g_zz_xxzz_1 = pbuffer.data(idx_eri_1_dg + 80);

    auto g_zz_xyyy_1 = pbuffer.data(idx_eri_1_dg + 81);

    auto g_zz_xyyz_1 = pbuffer.data(idx_eri_1_dg + 82);

    auto g_zz_xyzz_1 = pbuffer.data(idx_eri_1_dg + 83);

    auto g_zz_xzzz_1 = pbuffer.data(idx_eri_1_dg + 84);

    auto g_zz_yyyy_1 = pbuffer.data(idx_eri_1_dg + 85);

    auto g_zz_yyyz_1 = pbuffer.data(idx_eri_1_dg + 86);

    auto g_zz_yyzz_1 = pbuffer.data(idx_eri_1_dg + 87);

    auto g_zz_yzzz_1 = pbuffer.data(idx_eri_1_dg + 88);

    auto g_zz_zzzz_1 = pbuffer.data(idx_eri_1_dg + 89);

    // Set up components of auxiliary buffer : DH

    auto g_xx_xxxxx_1 = pbuffer.data(idx_eri_1_dh);

    auto g_xx_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 1);

    auto g_xx_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 2);

    auto g_xx_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 3);

    auto g_xx_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 4);

    auto g_xx_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 5);

    auto g_xx_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 6);

    auto g_xx_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 7);

    auto g_xx_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 8);

    auto g_xx_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 9);

    auto g_xx_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 10);

    auto g_xx_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 11);

    auto g_xx_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 12);

    auto g_xx_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 13);

    auto g_xx_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 14);

    auto g_xx_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 15);

    auto g_xx_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 16);

    auto g_xx_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 17);

    auto g_xx_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 18);

    auto g_xx_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 19);

    auto g_xx_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 20);

    auto g_xy_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 22);

    auto g_xy_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 24);

    auto g_xy_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 27);

    auto g_xy_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 31);

    auto g_xz_xxxxx_1 = pbuffer.data(idx_eri_1_dh + 42);

    auto g_xz_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 44);

    auto g_xz_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 47);

    auto g_xz_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 51);

    auto g_xz_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 56);

    auto g_yy_xxxxx_1 = pbuffer.data(idx_eri_1_dh + 63);

    auto g_yy_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 64);

    auto g_yy_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 65);

    auto g_yy_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 66);

    auto g_yy_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 67);

    auto g_yy_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 68);

    auto g_yy_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 69);

    auto g_yy_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 70);

    auto g_yy_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 71);

    auto g_yy_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 72);

    auto g_yy_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 73);

    auto g_yy_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 74);

    auto g_yy_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 75);

    auto g_yy_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 76);

    auto g_yy_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 77);

    auto g_yy_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 78);

    auto g_yy_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 79);

    auto g_yy_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 80);

    auto g_yy_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 81);

    auto g_yy_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 82);

    auto g_yy_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 83);

    auto g_yz_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 88);

    auto g_yz_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 91);

    auto g_yz_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 92);

    auto g_yz_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 95);

    auto g_yz_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 96);

    auto g_yz_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 97);

    auto g_yz_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 99);

    auto g_yz_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 100);

    auto g_yz_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 101);

    auto g_yz_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 102);

    auto g_yz_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 103);

    auto g_yz_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 104);

    auto g_zz_xxxxx_1 = pbuffer.data(idx_eri_1_dh + 105);

    auto g_zz_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 106);

    auto g_zz_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 107);

    auto g_zz_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 108);

    auto g_zz_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 109);

    auto g_zz_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 110);

    auto g_zz_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 111);

    auto g_zz_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 112);

    auto g_zz_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 113);

    auto g_zz_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 114);

    auto g_zz_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 115);

    auto g_zz_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 116);

    auto g_zz_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 117);

    auto g_zz_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 118);

    auto g_zz_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 119);

    auto g_zz_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 120);

    auto g_zz_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 121);

    auto g_zz_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 122);

    auto g_zz_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 123);

    auto g_zz_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 124);

    auto g_zz_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 125);

    // Set up 0-21 components of targeted buffer : FH

    auto g_xxx_xxxxx_0 = pbuffer.data(idx_eri_0_fh);

    auto g_xxx_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 1);

    auto g_xxx_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 2);

    auto g_xxx_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 3);

    auto g_xxx_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 4);

    auto g_xxx_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 5);

    auto g_xxx_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 6);

    auto g_xxx_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 7);

    auto g_xxx_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 8);

    auto g_xxx_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 9);

    auto g_xxx_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 10);

    auto g_xxx_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 11);

    auto g_xxx_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 12);

    auto g_xxx_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 13);

    auto g_xxx_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 14);

    auto g_xxx_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 15);

    auto g_xxx_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 16);

    auto g_xxx_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 17);

    auto g_xxx_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 18);

    auto g_xxx_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 19);

    auto g_xxx_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 20);

    #pragma omp simd aligned(g_x_xxxxx_0, g_x_xxxxx_1, g_x_xxxxy_0, g_x_xxxxy_1, g_x_xxxxz_0, g_x_xxxxz_1, g_x_xxxyy_0, g_x_xxxyy_1, g_x_xxxyz_0, g_x_xxxyz_1, g_x_xxxzz_0, g_x_xxxzz_1, g_x_xxyyy_0, g_x_xxyyy_1, g_x_xxyyz_0, g_x_xxyyz_1, g_x_xxyzz_0, g_x_xxyzz_1, g_x_xxzzz_0, g_x_xxzzz_1, g_x_xyyyy_0, g_x_xyyyy_1, g_x_xyyyz_0, g_x_xyyyz_1, g_x_xyyzz_0, g_x_xyyzz_1, g_x_xyzzz_0, g_x_xyzzz_1, g_x_xzzzz_0, g_x_xzzzz_1, g_x_yyyyy_0, g_x_yyyyy_1, g_x_yyyyz_0, g_x_yyyyz_1, g_x_yyyzz_0, g_x_yyyzz_1, g_x_yyzzz_0, g_x_yyzzz_1, g_x_yzzzz_0, g_x_yzzzz_1, g_x_zzzzz_0, g_x_zzzzz_1, g_xx_xxxx_1, g_xx_xxxxx_1, g_xx_xxxxy_1, g_xx_xxxxz_1, g_xx_xxxy_1, g_xx_xxxyy_1, g_xx_xxxyz_1, g_xx_xxxz_1, g_xx_xxxzz_1, g_xx_xxyy_1, g_xx_xxyyy_1, g_xx_xxyyz_1, g_xx_xxyz_1, g_xx_xxyzz_1, g_xx_xxzz_1, g_xx_xxzzz_1, g_xx_xyyy_1, g_xx_xyyyy_1, g_xx_xyyyz_1, g_xx_xyyz_1, g_xx_xyyzz_1, g_xx_xyzz_1, g_xx_xyzzz_1, g_xx_xzzz_1, g_xx_xzzzz_1, g_xx_yyyy_1, g_xx_yyyyy_1, g_xx_yyyyz_1, g_xx_yyyz_1, g_xx_yyyzz_1, g_xx_yyzz_1, g_xx_yyzzz_1, g_xx_yzzz_1, g_xx_yzzzz_1, g_xx_zzzz_1, g_xx_zzzzz_1, g_xxx_xxxxx_0, g_xxx_xxxxy_0, g_xxx_xxxxz_0, g_xxx_xxxyy_0, g_xxx_xxxyz_0, g_xxx_xxxzz_0, g_xxx_xxyyy_0, g_xxx_xxyyz_0, g_xxx_xxyzz_0, g_xxx_xxzzz_0, g_xxx_xyyyy_0, g_xxx_xyyyz_0, g_xxx_xyyzz_0, g_xxx_xyzzz_0, g_xxx_xzzzz_0, g_xxx_yyyyy_0, g_xxx_yyyyz_0, g_xxx_yyyzz_0, g_xxx_yyzzz_0, g_xxx_yzzzz_0, g_xxx_zzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxx_xxxxx_0[i] = 2.0 * g_x_xxxxx_0[i] * fbe_0 - 2.0 * g_x_xxxxx_1[i] * fz_be_0 + 5.0 * g_xx_xxxx_1[i] * fe_0 + g_xx_xxxxx_1[i] * pa_x[i];

        g_xxx_xxxxy_0[i] = 2.0 * g_x_xxxxy_0[i] * fbe_0 - 2.0 * g_x_xxxxy_1[i] * fz_be_0 + 4.0 * g_xx_xxxy_1[i] * fe_0 + g_xx_xxxxy_1[i] * pa_x[i];

        g_xxx_xxxxz_0[i] = 2.0 * g_x_xxxxz_0[i] * fbe_0 - 2.0 * g_x_xxxxz_1[i] * fz_be_0 + 4.0 * g_xx_xxxz_1[i] * fe_0 + g_xx_xxxxz_1[i] * pa_x[i];

        g_xxx_xxxyy_0[i] = 2.0 * g_x_xxxyy_0[i] * fbe_0 - 2.0 * g_x_xxxyy_1[i] * fz_be_0 + 3.0 * g_xx_xxyy_1[i] * fe_0 + g_xx_xxxyy_1[i] * pa_x[i];

        g_xxx_xxxyz_0[i] = 2.0 * g_x_xxxyz_0[i] * fbe_0 - 2.0 * g_x_xxxyz_1[i] * fz_be_0 + 3.0 * g_xx_xxyz_1[i] * fe_0 + g_xx_xxxyz_1[i] * pa_x[i];

        g_xxx_xxxzz_0[i] = 2.0 * g_x_xxxzz_0[i] * fbe_0 - 2.0 * g_x_xxxzz_1[i] * fz_be_0 + 3.0 * g_xx_xxzz_1[i] * fe_0 + g_xx_xxxzz_1[i] * pa_x[i];

        g_xxx_xxyyy_0[i] = 2.0 * g_x_xxyyy_0[i] * fbe_0 - 2.0 * g_x_xxyyy_1[i] * fz_be_0 + 2.0 * g_xx_xyyy_1[i] * fe_0 + g_xx_xxyyy_1[i] * pa_x[i];

        g_xxx_xxyyz_0[i] = 2.0 * g_x_xxyyz_0[i] * fbe_0 - 2.0 * g_x_xxyyz_1[i] * fz_be_0 + 2.0 * g_xx_xyyz_1[i] * fe_0 + g_xx_xxyyz_1[i] * pa_x[i];

        g_xxx_xxyzz_0[i] = 2.0 * g_x_xxyzz_0[i] * fbe_0 - 2.0 * g_x_xxyzz_1[i] * fz_be_0 + 2.0 * g_xx_xyzz_1[i] * fe_0 + g_xx_xxyzz_1[i] * pa_x[i];

        g_xxx_xxzzz_0[i] = 2.0 * g_x_xxzzz_0[i] * fbe_0 - 2.0 * g_x_xxzzz_1[i] * fz_be_0 + 2.0 * g_xx_xzzz_1[i] * fe_0 + g_xx_xxzzz_1[i] * pa_x[i];

        g_xxx_xyyyy_0[i] = 2.0 * g_x_xyyyy_0[i] * fbe_0 - 2.0 * g_x_xyyyy_1[i] * fz_be_0 + g_xx_yyyy_1[i] * fe_0 + g_xx_xyyyy_1[i] * pa_x[i];

        g_xxx_xyyyz_0[i] = 2.0 * g_x_xyyyz_0[i] * fbe_0 - 2.0 * g_x_xyyyz_1[i] * fz_be_0 + g_xx_yyyz_1[i] * fe_0 + g_xx_xyyyz_1[i] * pa_x[i];

        g_xxx_xyyzz_0[i] = 2.0 * g_x_xyyzz_0[i] * fbe_0 - 2.0 * g_x_xyyzz_1[i] * fz_be_0 + g_xx_yyzz_1[i] * fe_0 + g_xx_xyyzz_1[i] * pa_x[i];

        g_xxx_xyzzz_0[i] = 2.0 * g_x_xyzzz_0[i] * fbe_0 - 2.0 * g_x_xyzzz_1[i] * fz_be_0 + g_xx_yzzz_1[i] * fe_0 + g_xx_xyzzz_1[i] * pa_x[i];

        g_xxx_xzzzz_0[i] = 2.0 * g_x_xzzzz_0[i] * fbe_0 - 2.0 * g_x_xzzzz_1[i] * fz_be_0 + g_xx_zzzz_1[i] * fe_0 + g_xx_xzzzz_1[i] * pa_x[i];

        g_xxx_yyyyy_0[i] = 2.0 * g_x_yyyyy_0[i] * fbe_0 - 2.0 * g_x_yyyyy_1[i] * fz_be_0 + g_xx_yyyyy_1[i] * pa_x[i];

        g_xxx_yyyyz_0[i] = 2.0 * g_x_yyyyz_0[i] * fbe_0 - 2.0 * g_x_yyyyz_1[i] * fz_be_0 + g_xx_yyyyz_1[i] * pa_x[i];

        g_xxx_yyyzz_0[i] = 2.0 * g_x_yyyzz_0[i] * fbe_0 - 2.0 * g_x_yyyzz_1[i] * fz_be_0 + g_xx_yyyzz_1[i] * pa_x[i];

        g_xxx_yyzzz_0[i] = 2.0 * g_x_yyzzz_0[i] * fbe_0 - 2.0 * g_x_yyzzz_1[i] * fz_be_0 + g_xx_yyzzz_1[i] * pa_x[i];

        g_xxx_yzzzz_0[i] = 2.0 * g_x_yzzzz_0[i] * fbe_0 - 2.0 * g_x_yzzzz_1[i] * fz_be_0 + g_xx_yzzzz_1[i] * pa_x[i];

        g_xxx_zzzzz_0[i] = 2.0 * g_x_zzzzz_0[i] * fbe_0 - 2.0 * g_x_zzzzz_1[i] * fz_be_0 + g_xx_zzzzz_1[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : FH

    auto g_xxy_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 21);

    auto g_xxy_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 22);

    auto g_xxy_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 23);

    auto g_xxy_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 24);

    auto g_xxy_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 25);

    auto g_xxy_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 26);

    auto g_xxy_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 27);

    auto g_xxy_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 28);

    auto g_xxy_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 29);

    auto g_xxy_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 30);

    auto g_xxy_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 31);

    auto g_xxy_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 32);

    auto g_xxy_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 33);

    auto g_xxy_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 34);

    auto g_xxy_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 35);

    auto g_xxy_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 36);

    auto g_xxy_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 37);

    auto g_xxy_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 38);

    auto g_xxy_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 39);

    auto g_xxy_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 40);

    auto g_xxy_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 41);

    #pragma omp simd aligned(g_xx_xxxx_1, g_xx_xxxxx_1, g_xx_xxxxy_1, g_xx_xxxxz_1, g_xx_xxxy_1, g_xx_xxxyy_1, g_xx_xxxyz_1, g_xx_xxxz_1, g_xx_xxxzz_1, g_xx_xxyy_1, g_xx_xxyyy_1, g_xx_xxyyz_1, g_xx_xxyz_1, g_xx_xxyzz_1, g_xx_xxzz_1, g_xx_xxzzz_1, g_xx_xyyy_1, g_xx_xyyyy_1, g_xx_xyyyz_1, g_xx_xyyz_1, g_xx_xyyzz_1, g_xx_xyzz_1, g_xx_xyzzz_1, g_xx_xzzz_1, g_xx_xzzzz_1, g_xx_yyyy_1, g_xx_yyyyy_1, g_xx_yyyyz_1, g_xx_yyyz_1, g_xx_yyyzz_1, g_xx_yyzz_1, g_xx_yyzzz_1, g_xx_yzzz_1, g_xx_yzzzz_1, g_xx_zzzz_1, g_xx_zzzzz_1, g_xxy_xxxxx_0, g_xxy_xxxxy_0, g_xxy_xxxxz_0, g_xxy_xxxyy_0, g_xxy_xxxyz_0, g_xxy_xxxzz_0, g_xxy_xxyyy_0, g_xxy_xxyyz_0, g_xxy_xxyzz_0, g_xxy_xxzzz_0, g_xxy_xyyyy_0, g_xxy_xyyyz_0, g_xxy_xyyzz_0, g_xxy_xyzzz_0, g_xxy_xzzzz_0, g_xxy_yyyyy_0, g_xxy_yyyyz_0, g_xxy_yyyzz_0, g_xxy_yyzzz_0, g_xxy_yzzzz_0, g_xxy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxy_xxxxx_0[i] = g_xx_xxxxx_1[i] * pa_y[i];

        g_xxy_xxxxy_0[i] = g_xx_xxxx_1[i] * fe_0 + g_xx_xxxxy_1[i] * pa_y[i];

        g_xxy_xxxxz_0[i] = g_xx_xxxxz_1[i] * pa_y[i];

        g_xxy_xxxyy_0[i] = 2.0 * g_xx_xxxy_1[i] * fe_0 + g_xx_xxxyy_1[i] * pa_y[i];

        g_xxy_xxxyz_0[i] = g_xx_xxxz_1[i] * fe_0 + g_xx_xxxyz_1[i] * pa_y[i];

        g_xxy_xxxzz_0[i] = g_xx_xxxzz_1[i] * pa_y[i];

        g_xxy_xxyyy_0[i] = 3.0 * g_xx_xxyy_1[i] * fe_0 + g_xx_xxyyy_1[i] * pa_y[i];

        g_xxy_xxyyz_0[i] = 2.0 * g_xx_xxyz_1[i] * fe_0 + g_xx_xxyyz_1[i] * pa_y[i];

        g_xxy_xxyzz_0[i] = g_xx_xxzz_1[i] * fe_0 + g_xx_xxyzz_1[i] * pa_y[i];

        g_xxy_xxzzz_0[i] = g_xx_xxzzz_1[i] * pa_y[i];

        g_xxy_xyyyy_0[i] = 4.0 * g_xx_xyyy_1[i] * fe_0 + g_xx_xyyyy_1[i] * pa_y[i];

        g_xxy_xyyyz_0[i] = 3.0 * g_xx_xyyz_1[i] * fe_0 + g_xx_xyyyz_1[i] * pa_y[i];

        g_xxy_xyyzz_0[i] = 2.0 * g_xx_xyzz_1[i] * fe_0 + g_xx_xyyzz_1[i] * pa_y[i];

        g_xxy_xyzzz_0[i] = g_xx_xzzz_1[i] * fe_0 + g_xx_xyzzz_1[i] * pa_y[i];

        g_xxy_xzzzz_0[i] = g_xx_xzzzz_1[i] * pa_y[i];

        g_xxy_yyyyy_0[i] = 5.0 * g_xx_yyyy_1[i] * fe_0 + g_xx_yyyyy_1[i] * pa_y[i];

        g_xxy_yyyyz_0[i] = 4.0 * g_xx_yyyz_1[i] * fe_0 + g_xx_yyyyz_1[i] * pa_y[i];

        g_xxy_yyyzz_0[i] = 3.0 * g_xx_yyzz_1[i] * fe_0 + g_xx_yyyzz_1[i] * pa_y[i];

        g_xxy_yyzzz_0[i] = 2.0 * g_xx_yzzz_1[i] * fe_0 + g_xx_yyzzz_1[i] * pa_y[i];

        g_xxy_yzzzz_0[i] = g_xx_zzzz_1[i] * fe_0 + g_xx_yzzzz_1[i] * pa_y[i];

        g_xxy_zzzzz_0[i] = g_xx_zzzzz_1[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : FH

    auto g_xxz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 42);

    auto g_xxz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 43);

    auto g_xxz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 44);

    auto g_xxz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 45);

    auto g_xxz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 46);

    auto g_xxz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 47);

    auto g_xxz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 48);

    auto g_xxz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 49);

    auto g_xxz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 50);

    auto g_xxz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 51);

    auto g_xxz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 52);

    auto g_xxz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 53);

    auto g_xxz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 54);

    auto g_xxz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 55);

    auto g_xxz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 56);

    auto g_xxz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 57);

    auto g_xxz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 58);

    auto g_xxz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 59);

    auto g_xxz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 60);

    auto g_xxz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 61);

    auto g_xxz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 62);

    #pragma omp simd aligned(g_xx_xxxx_1, g_xx_xxxxx_1, g_xx_xxxxy_1, g_xx_xxxxz_1, g_xx_xxxy_1, g_xx_xxxyy_1, g_xx_xxxyz_1, g_xx_xxxz_1, g_xx_xxxzz_1, g_xx_xxyy_1, g_xx_xxyyy_1, g_xx_xxyyz_1, g_xx_xxyz_1, g_xx_xxyzz_1, g_xx_xxzz_1, g_xx_xxzzz_1, g_xx_xyyy_1, g_xx_xyyyy_1, g_xx_xyyyz_1, g_xx_xyyz_1, g_xx_xyyzz_1, g_xx_xyzz_1, g_xx_xyzzz_1, g_xx_xzzz_1, g_xx_xzzzz_1, g_xx_yyyy_1, g_xx_yyyyy_1, g_xx_yyyyz_1, g_xx_yyyz_1, g_xx_yyyzz_1, g_xx_yyzz_1, g_xx_yyzzz_1, g_xx_yzzz_1, g_xx_yzzzz_1, g_xx_zzzz_1, g_xx_zzzzz_1, g_xxz_xxxxx_0, g_xxz_xxxxy_0, g_xxz_xxxxz_0, g_xxz_xxxyy_0, g_xxz_xxxyz_0, g_xxz_xxxzz_0, g_xxz_xxyyy_0, g_xxz_xxyyz_0, g_xxz_xxyzz_0, g_xxz_xxzzz_0, g_xxz_xyyyy_0, g_xxz_xyyyz_0, g_xxz_xyyzz_0, g_xxz_xyzzz_0, g_xxz_xzzzz_0, g_xxz_yyyyy_0, g_xxz_yyyyz_0, g_xxz_yyyzz_0, g_xxz_yyzzz_0, g_xxz_yzzzz_0, g_xxz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxz_xxxxx_0[i] = g_xx_xxxxx_1[i] * pa_z[i];

        g_xxz_xxxxy_0[i] = g_xx_xxxxy_1[i] * pa_z[i];

        g_xxz_xxxxz_0[i] = g_xx_xxxx_1[i] * fe_0 + g_xx_xxxxz_1[i] * pa_z[i];

        g_xxz_xxxyy_0[i] = g_xx_xxxyy_1[i] * pa_z[i];

        g_xxz_xxxyz_0[i] = g_xx_xxxy_1[i] * fe_0 + g_xx_xxxyz_1[i] * pa_z[i];

        g_xxz_xxxzz_0[i] = 2.0 * g_xx_xxxz_1[i] * fe_0 + g_xx_xxxzz_1[i] * pa_z[i];

        g_xxz_xxyyy_0[i] = g_xx_xxyyy_1[i] * pa_z[i];

        g_xxz_xxyyz_0[i] = g_xx_xxyy_1[i] * fe_0 + g_xx_xxyyz_1[i] * pa_z[i];

        g_xxz_xxyzz_0[i] = 2.0 * g_xx_xxyz_1[i] * fe_0 + g_xx_xxyzz_1[i] * pa_z[i];

        g_xxz_xxzzz_0[i] = 3.0 * g_xx_xxzz_1[i] * fe_0 + g_xx_xxzzz_1[i] * pa_z[i];

        g_xxz_xyyyy_0[i] = g_xx_xyyyy_1[i] * pa_z[i];

        g_xxz_xyyyz_0[i] = g_xx_xyyy_1[i] * fe_0 + g_xx_xyyyz_1[i] * pa_z[i];

        g_xxz_xyyzz_0[i] = 2.0 * g_xx_xyyz_1[i] * fe_0 + g_xx_xyyzz_1[i] * pa_z[i];

        g_xxz_xyzzz_0[i] = 3.0 * g_xx_xyzz_1[i] * fe_0 + g_xx_xyzzz_1[i] * pa_z[i];

        g_xxz_xzzzz_0[i] = 4.0 * g_xx_xzzz_1[i] * fe_0 + g_xx_xzzzz_1[i] * pa_z[i];

        g_xxz_yyyyy_0[i] = g_xx_yyyyy_1[i] * pa_z[i];

        g_xxz_yyyyz_0[i] = g_xx_yyyy_1[i] * fe_0 + g_xx_yyyyz_1[i] * pa_z[i];

        g_xxz_yyyzz_0[i] = 2.0 * g_xx_yyyz_1[i] * fe_0 + g_xx_yyyzz_1[i] * pa_z[i];

        g_xxz_yyzzz_0[i] = 3.0 * g_xx_yyzz_1[i] * fe_0 + g_xx_yyzzz_1[i] * pa_z[i];

        g_xxz_yzzzz_0[i] = 4.0 * g_xx_yzzz_1[i] * fe_0 + g_xx_yzzzz_1[i] * pa_z[i];

        g_xxz_zzzzz_0[i] = 5.0 * g_xx_zzzz_1[i] * fe_0 + g_xx_zzzzz_1[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : FH

    auto g_xyy_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 63);

    auto g_xyy_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 64);

    auto g_xyy_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 65);

    auto g_xyy_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 66);

    auto g_xyy_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 67);

    auto g_xyy_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 68);

    auto g_xyy_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 69);

    auto g_xyy_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 70);

    auto g_xyy_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 71);

    auto g_xyy_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 72);

    auto g_xyy_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 73);

    auto g_xyy_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 74);

    auto g_xyy_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 75);

    auto g_xyy_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 76);

    auto g_xyy_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 77);

    auto g_xyy_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 78);

    auto g_xyy_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 79);

    auto g_xyy_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 80);

    auto g_xyy_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 81);

    auto g_xyy_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 82);

    auto g_xyy_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 83);

    #pragma omp simd aligned(g_xyy_xxxxx_0, g_xyy_xxxxy_0, g_xyy_xxxxz_0, g_xyy_xxxyy_0, g_xyy_xxxyz_0, g_xyy_xxxzz_0, g_xyy_xxyyy_0, g_xyy_xxyyz_0, g_xyy_xxyzz_0, g_xyy_xxzzz_0, g_xyy_xyyyy_0, g_xyy_xyyyz_0, g_xyy_xyyzz_0, g_xyy_xyzzz_0, g_xyy_xzzzz_0, g_xyy_yyyyy_0, g_xyy_yyyyz_0, g_xyy_yyyzz_0, g_xyy_yyzzz_0, g_xyy_yzzzz_0, g_xyy_zzzzz_0, g_yy_xxxx_1, g_yy_xxxxx_1, g_yy_xxxxy_1, g_yy_xxxxz_1, g_yy_xxxy_1, g_yy_xxxyy_1, g_yy_xxxyz_1, g_yy_xxxz_1, g_yy_xxxzz_1, g_yy_xxyy_1, g_yy_xxyyy_1, g_yy_xxyyz_1, g_yy_xxyz_1, g_yy_xxyzz_1, g_yy_xxzz_1, g_yy_xxzzz_1, g_yy_xyyy_1, g_yy_xyyyy_1, g_yy_xyyyz_1, g_yy_xyyz_1, g_yy_xyyzz_1, g_yy_xyzz_1, g_yy_xyzzz_1, g_yy_xzzz_1, g_yy_xzzzz_1, g_yy_yyyy_1, g_yy_yyyyy_1, g_yy_yyyyz_1, g_yy_yyyz_1, g_yy_yyyzz_1, g_yy_yyzz_1, g_yy_yyzzz_1, g_yy_yzzz_1, g_yy_yzzzz_1, g_yy_zzzz_1, g_yy_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyy_xxxxx_0[i] = 5.0 * g_yy_xxxx_1[i] * fe_0 + g_yy_xxxxx_1[i] * pa_x[i];

        g_xyy_xxxxy_0[i] = 4.0 * g_yy_xxxy_1[i] * fe_0 + g_yy_xxxxy_1[i] * pa_x[i];

        g_xyy_xxxxz_0[i] = 4.0 * g_yy_xxxz_1[i] * fe_0 + g_yy_xxxxz_1[i] * pa_x[i];

        g_xyy_xxxyy_0[i] = 3.0 * g_yy_xxyy_1[i] * fe_0 + g_yy_xxxyy_1[i] * pa_x[i];

        g_xyy_xxxyz_0[i] = 3.0 * g_yy_xxyz_1[i] * fe_0 + g_yy_xxxyz_1[i] * pa_x[i];

        g_xyy_xxxzz_0[i] = 3.0 * g_yy_xxzz_1[i] * fe_0 + g_yy_xxxzz_1[i] * pa_x[i];

        g_xyy_xxyyy_0[i] = 2.0 * g_yy_xyyy_1[i] * fe_0 + g_yy_xxyyy_1[i] * pa_x[i];

        g_xyy_xxyyz_0[i] = 2.0 * g_yy_xyyz_1[i] * fe_0 + g_yy_xxyyz_1[i] * pa_x[i];

        g_xyy_xxyzz_0[i] = 2.0 * g_yy_xyzz_1[i] * fe_0 + g_yy_xxyzz_1[i] * pa_x[i];

        g_xyy_xxzzz_0[i] = 2.0 * g_yy_xzzz_1[i] * fe_0 + g_yy_xxzzz_1[i] * pa_x[i];

        g_xyy_xyyyy_0[i] = g_yy_yyyy_1[i] * fe_0 + g_yy_xyyyy_1[i] * pa_x[i];

        g_xyy_xyyyz_0[i] = g_yy_yyyz_1[i] * fe_0 + g_yy_xyyyz_1[i] * pa_x[i];

        g_xyy_xyyzz_0[i] = g_yy_yyzz_1[i] * fe_0 + g_yy_xyyzz_1[i] * pa_x[i];

        g_xyy_xyzzz_0[i] = g_yy_yzzz_1[i] * fe_0 + g_yy_xyzzz_1[i] * pa_x[i];

        g_xyy_xzzzz_0[i] = g_yy_zzzz_1[i] * fe_0 + g_yy_xzzzz_1[i] * pa_x[i];

        g_xyy_yyyyy_0[i] = g_yy_yyyyy_1[i] * pa_x[i];

        g_xyy_yyyyz_0[i] = g_yy_yyyyz_1[i] * pa_x[i];

        g_xyy_yyyzz_0[i] = g_yy_yyyzz_1[i] * pa_x[i];

        g_xyy_yyzzz_0[i] = g_yy_yyzzz_1[i] * pa_x[i];

        g_xyy_yzzzz_0[i] = g_yy_yzzzz_1[i] * pa_x[i];

        g_xyy_zzzzz_0[i] = g_yy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : FH

    auto g_xyz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 84);

    auto g_xyz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 85);

    auto g_xyz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 86);

    auto g_xyz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 87);

    auto g_xyz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 88);

    auto g_xyz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 89);

    auto g_xyz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 90);

    auto g_xyz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 91);

    auto g_xyz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 92);

    auto g_xyz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 93);

    auto g_xyz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 94);

    auto g_xyz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 95);

    auto g_xyz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 96);

    auto g_xyz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 97);

    auto g_xyz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 98);

    auto g_xyz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 99);

    auto g_xyz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 100);

    auto g_xyz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 101);

    auto g_xyz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 102);

    auto g_xyz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 103);

    auto g_xyz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 104);

    #pragma omp simd aligned(g_xy_xxxxy_1, g_xy_xxxyy_1, g_xy_xxyyy_1, g_xy_xyyyy_1, g_xyz_xxxxx_0, g_xyz_xxxxy_0, g_xyz_xxxxz_0, g_xyz_xxxyy_0, g_xyz_xxxyz_0, g_xyz_xxxzz_0, g_xyz_xxyyy_0, g_xyz_xxyyz_0, g_xyz_xxyzz_0, g_xyz_xxzzz_0, g_xyz_xyyyy_0, g_xyz_xyyyz_0, g_xyz_xyyzz_0, g_xyz_xyzzz_0, g_xyz_xzzzz_0, g_xyz_yyyyy_0, g_xyz_yyyyz_0, g_xyz_yyyzz_0, g_xyz_yyzzz_0, g_xyz_yzzzz_0, g_xyz_zzzzz_0, g_xz_xxxxx_1, g_xz_xxxxz_1, g_xz_xxxzz_1, g_xz_xxzzz_1, g_xz_xzzzz_1, g_yz_xxxyz_1, g_yz_xxyyz_1, g_yz_xxyz_1, g_yz_xxyzz_1, g_yz_xyyyz_1, g_yz_xyyz_1, g_yz_xyyzz_1, g_yz_xyzz_1, g_yz_xyzzz_1, g_yz_yyyyy_1, g_yz_yyyyz_1, g_yz_yyyz_1, g_yz_yyyzz_1, g_yz_yyzz_1, g_yz_yyzzz_1, g_yz_yzzz_1, g_yz_yzzzz_1, g_yz_zzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyz_xxxxx_0[i] = g_xz_xxxxx_1[i] * pa_y[i];

        g_xyz_xxxxy_0[i] = g_xy_xxxxy_1[i] * pa_z[i];

        g_xyz_xxxxz_0[i] = g_xz_xxxxz_1[i] * pa_y[i];

        g_xyz_xxxyy_0[i] = g_xy_xxxyy_1[i] * pa_z[i];

        g_xyz_xxxyz_0[i] = 3.0 * g_yz_xxyz_1[i] * fe_0 + g_yz_xxxyz_1[i] * pa_x[i];

        g_xyz_xxxzz_0[i] = g_xz_xxxzz_1[i] * pa_y[i];

        g_xyz_xxyyy_0[i] = g_xy_xxyyy_1[i] * pa_z[i];

        g_xyz_xxyyz_0[i] = 2.0 * g_yz_xyyz_1[i] * fe_0 + g_yz_xxyyz_1[i] * pa_x[i];

        g_xyz_xxyzz_0[i] = 2.0 * g_yz_xyzz_1[i] * fe_0 + g_yz_xxyzz_1[i] * pa_x[i];

        g_xyz_xxzzz_0[i] = g_xz_xxzzz_1[i] * pa_y[i];

        g_xyz_xyyyy_0[i] = g_xy_xyyyy_1[i] * pa_z[i];

        g_xyz_xyyyz_0[i] = g_yz_yyyz_1[i] * fe_0 + g_yz_xyyyz_1[i] * pa_x[i];

        g_xyz_xyyzz_0[i] = g_yz_yyzz_1[i] * fe_0 + g_yz_xyyzz_1[i] * pa_x[i];

        g_xyz_xyzzz_0[i] = g_yz_yzzz_1[i] * fe_0 + g_yz_xyzzz_1[i] * pa_x[i];

        g_xyz_xzzzz_0[i] = g_xz_xzzzz_1[i] * pa_y[i];

        g_xyz_yyyyy_0[i] = g_yz_yyyyy_1[i] * pa_x[i];

        g_xyz_yyyyz_0[i] = g_yz_yyyyz_1[i] * pa_x[i];

        g_xyz_yyyzz_0[i] = g_yz_yyyzz_1[i] * pa_x[i];

        g_xyz_yyzzz_0[i] = g_yz_yyzzz_1[i] * pa_x[i];

        g_xyz_yzzzz_0[i] = g_yz_yzzzz_1[i] * pa_x[i];

        g_xyz_zzzzz_0[i] = g_yz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 105-126 components of targeted buffer : FH

    auto g_xzz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 105);

    auto g_xzz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 106);

    auto g_xzz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 107);

    auto g_xzz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 108);

    auto g_xzz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 109);

    auto g_xzz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 110);

    auto g_xzz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 111);

    auto g_xzz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 112);

    auto g_xzz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 113);

    auto g_xzz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 114);

    auto g_xzz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 115);

    auto g_xzz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 116);

    auto g_xzz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 117);

    auto g_xzz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 118);

    auto g_xzz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 119);

    auto g_xzz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 120);

    auto g_xzz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 121);

    auto g_xzz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 122);

    auto g_xzz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 123);

    auto g_xzz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 124);

    auto g_xzz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 125);

    #pragma omp simd aligned(g_xzz_xxxxx_0, g_xzz_xxxxy_0, g_xzz_xxxxz_0, g_xzz_xxxyy_0, g_xzz_xxxyz_0, g_xzz_xxxzz_0, g_xzz_xxyyy_0, g_xzz_xxyyz_0, g_xzz_xxyzz_0, g_xzz_xxzzz_0, g_xzz_xyyyy_0, g_xzz_xyyyz_0, g_xzz_xyyzz_0, g_xzz_xyzzz_0, g_xzz_xzzzz_0, g_xzz_yyyyy_0, g_xzz_yyyyz_0, g_xzz_yyyzz_0, g_xzz_yyzzz_0, g_xzz_yzzzz_0, g_xzz_zzzzz_0, g_zz_xxxx_1, g_zz_xxxxx_1, g_zz_xxxxy_1, g_zz_xxxxz_1, g_zz_xxxy_1, g_zz_xxxyy_1, g_zz_xxxyz_1, g_zz_xxxz_1, g_zz_xxxzz_1, g_zz_xxyy_1, g_zz_xxyyy_1, g_zz_xxyyz_1, g_zz_xxyz_1, g_zz_xxyzz_1, g_zz_xxzz_1, g_zz_xxzzz_1, g_zz_xyyy_1, g_zz_xyyyy_1, g_zz_xyyyz_1, g_zz_xyyz_1, g_zz_xyyzz_1, g_zz_xyzz_1, g_zz_xyzzz_1, g_zz_xzzz_1, g_zz_xzzzz_1, g_zz_yyyy_1, g_zz_yyyyy_1, g_zz_yyyyz_1, g_zz_yyyz_1, g_zz_yyyzz_1, g_zz_yyzz_1, g_zz_yyzzz_1, g_zz_yzzz_1, g_zz_yzzzz_1, g_zz_zzzz_1, g_zz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzz_xxxxx_0[i] = 5.0 * g_zz_xxxx_1[i] * fe_0 + g_zz_xxxxx_1[i] * pa_x[i];

        g_xzz_xxxxy_0[i] = 4.0 * g_zz_xxxy_1[i] * fe_0 + g_zz_xxxxy_1[i] * pa_x[i];

        g_xzz_xxxxz_0[i] = 4.0 * g_zz_xxxz_1[i] * fe_0 + g_zz_xxxxz_1[i] * pa_x[i];

        g_xzz_xxxyy_0[i] = 3.0 * g_zz_xxyy_1[i] * fe_0 + g_zz_xxxyy_1[i] * pa_x[i];

        g_xzz_xxxyz_0[i] = 3.0 * g_zz_xxyz_1[i] * fe_0 + g_zz_xxxyz_1[i] * pa_x[i];

        g_xzz_xxxzz_0[i] = 3.0 * g_zz_xxzz_1[i] * fe_0 + g_zz_xxxzz_1[i] * pa_x[i];

        g_xzz_xxyyy_0[i] = 2.0 * g_zz_xyyy_1[i] * fe_0 + g_zz_xxyyy_1[i] * pa_x[i];

        g_xzz_xxyyz_0[i] = 2.0 * g_zz_xyyz_1[i] * fe_0 + g_zz_xxyyz_1[i] * pa_x[i];

        g_xzz_xxyzz_0[i] = 2.0 * g_zz_xyzz_1[i] * fe_0 + g_zz_xxyzz_1[i] * pa_x[i];

        g_xzz_xxzzz_0[i] = 2.0 * g_zz_xzzz_1[i] * fe_0 + g_zz_xxzzz_1[i] * pa_x[i];

        g_xzz_xyyyy_0[i] = g_zz_yyyy_1[i] * fe_0 + g_zz_xyyyy_1[i] * pa_x[i];

        g_xzz_xyyyz_0[i] = g_zz_yyyz_1[i] * fe_0 + g_zz_xyyyz_1[i] * pa_x[i];

        g_xzz_xyyzz_0[i] = g_zz_yyzz_1[i] * fe_0 + g_zz_xyyzz_1[i] * pa_x[i];

        g_xzz_xyzzz_0[i] = g_zz_yzzz_1[i] * fe_0 + g_zz_xyzzz_1[i] * pa_x[i];

        g_xzz_xzzzz_0[i] = g_zz_zzzz_1[i] * fe_0 + g_zz_xzzzz_1[i] * pa_x[i];

        g_xzz_yyyyy_0[i] = g_zz_yyyyy_1[i] * pa_x[i];

        g_xzz_yyyyz_0[i] = g_zz_yyyyz_1[i] * pa_x[i];

        g_xzz_yyyzz_0[i] = g_zz_yyyzz_1[i] * pa_x[i];

        g_xzz_yyzzz_0[i] = g_zz_yyzzz_1[i] * pa_x[i];

        g_xzz_yzzzz_0[i] = g_zz_yzzzz_1[i] * pa_x[i];

        g_xzz_zzzzz_0[i] = g_zz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : FH

    auto g_yyy_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 126);

    auto g_yyy_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 127);

    auto g_yyy_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 128);

    auto g_yyy_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 129);

    auto g_yyy_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 130);

    auto g_yyy_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 131);

    auto g_yyy_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 132);

    auto g_yyy_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 133);

    auto g_yyy_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 134);

    auto g_yyy_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 135);

    auto g_yyy_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 136);

    auto g_yyy_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 137);

    auto g_yyy_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 138);

    auto g_yyy_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 139);

    auto g_yyy_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 140);

    auto g_yyy_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 141);

    auto g_yyy_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 142);

    auto g_yyy_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 143);

    auto g_yyy_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 144);

    auto g_yyy_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 145);

    auto g_yyy_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 146);

    #pragma omp simd aligned(g_y_xxxxx_0, g_y_xxxxx_1, g_y_xxxxy_0, g_y_xxxxy_1, g_y_xxxxz_0, g_y_xxxxz_1, g_y_xxxyy_0, g_y_xxxyy_1, g_y_xxxyz_0, g_y_xxxyz_1, g_y_xxxzz_0, g_y_xxxzz_1, g_y_xxyyy_0, g_y_xxyyy_1, g_y_xxyyz_0, g_y_xxyyz_1, g_y_xxyzz_0, g_y_xxyzz_1, g_y_xxzzz_0, g_y_xxzzz_1, g_y_xyyyy_0, g_y_xyyyy_1, g_y_xyyyz_0, g_y_xyyyz_1, g_y_xyyzz_0, g_y_xyyzz_1, g_y_xyzzz_0, g_y_xyzzz_1, g_y_xzzzz_0, g_y_xzzzz_1, g_y_yyyyy_0, g_y_yyyyy_1, g_y_yyyyz_0, g_y_yyyyz_1, g_y_yyyzz_0, g_y_yyyzz_1, g_y_yyzzz_0, g_y_yyzzz_1, g_y_yzzzz_0, g_y_yzzzz_1, g_y_zzzzz_0, g_y_zzzzz_1, g_yy_xxxx_1, g_yy_xxxxx_1, g_yy_xxxxy_1, g_yy_xxxxz_1, g_yy_xxxy_1, g_yy_xxxyy_1, g_yy_xxxyz_1, g_yy_xxxz_1, g_yy_xxxzz_1, g_yy_xxyy_1, g_yy_xxyyy_1, g_yy_xxyyz_1, g_yy_xxyz_1, g_yy_xxyzz_1, g_yy_xxzz_1, g_yy_xxzzz_1, g_yy_xyyy_1, g_yy_xyyyy_1, g_yy_xyyyz_1, g_yy_xyyz_1, g_yy_xyyzz_1, g_yy_xyzz_1, g_yy_xyzzz_1, g_yy_xzzz_1, g_yy_xzzzz_1, g_yy_yyyy_1, g_yy_yyyyy_1, g_yy_yyyyz_1, g_yy_yyyz_1, g_yy_yyyzz_1, g_yy_yyzz_1, g_yy_yyzzz_1, g_yy_yzzz_1, g_yy_yzzzz_1, g_yy_zzzz_1, g_yy_zzzzz_1, g_yyy_xxxxx_0, g_yyy_xxxxy_0, g_yyy_xxxxz_0, g_yyy_xxxyy_0, g_yyy_xxxyz_0, g_yyy_xxxzz_0, g_yyy_xxyyy_0, g_yyy_xxyyz_0, g_yyy_xxyzz_0, g_yyy_xxzzz_0, g_yyy_xyyyy_0, g_yyy_xyyyz_0, g_yyy_xyyzz_0, g_yyy_xyzzz_0, g_yyy_xzzzz_0, g_yyy_yyyyy_0, g_yyy_yyyyz_0, g_yyy_yyyzz_0, g_yyy_yyzzz_0, g_yyy_yzzzz_0, g_yyy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyy_xxxxx_0[i] = 2.0 * g_y_xxxxx_0[i] * fbe_0 - 2.0 * g_y_xxxxx_1[i] * fz_be_0 + g_yy_xxxxx_1[i] * pa_y[i];

        g_yyy_xxxxy_0[i] = 2.0 * g_y_xxxxy_0[i] * fbe_0 - 2.0 * g_y_xxxxy_1[i] * fz_be_0 + g_yy_xxxx_1[i] * fe_0 + g_yy_xxxxy_1[i] * pa_y[i];

        g_yyy_xxxxz_0[i] = 2.0 * g_y_xxxxz_0[i] * fbe_0 - 2.0 * g_y_xxxxz_1[i] * fz_be_0 + g_yy_xxxxz_1[i] * pa_y[i];

        g_yyy_xxxyy_0[i] = 2.0 * g_y_xxxyy_0[i] * fbe_0 - 2.0 * g_y_xxxyy_1[i] * fz_be_0 + 2.0 * g_yy_xxxy_1[i] * fe_0 + g_yy_xxxyy_1[i] * pa_y[i];

        g_yyy_xxxyz_0[i] = 2.0 * g_y_xxxyz_0[i] * fbe_0 - 2.0 * g_y_xxxyz_1[i] * fz_be_0 + g_yy_xxxz_1[i] * fe_0 + g_yy_xxxyz_1[i] * pa_y[i];

        g_yyy_xxxzz_0[i] = 2.0 * g_y_xxxzz_0[i] * fbe_0 - 2.0 * g_y_xxxzz_1[i] * fz_be_0 + g_yy_xxxzz_1[i] * pa_y[i];

        g_yyy_xxyyy_0[i] = 2.0 * g_y_xxyyy_0[i] * fbe_0 - 2.0 * g_y_xxyyy_1[i] * fz_be_0 + 3.0 * g_yy_xxyy_1[i] * fe_0 + g_yy_xxyyy_1[i] * pa_y[i];

        g_yyy_xxyyz_0[i] = 2.0 * g_y_xxyyz_0[i] * fbe_0 - 2.0 * g_y_xxyyz_1[i] * fz_be_0 + 2.0 * g_yy_xxyz_1[i] * fe_0 + g_yy_xxyyz_1[i] * pa_y[i];

        g_yyy_xxyzz_0[i] = 2.0 * g_y_xxyzz_0[i] * fbe_0 - 2.0 * g_y_xxyzz_1[i] * fz_be_0 + g_yy_xxzz_1[i] * fe_0 + g_yy_xxyzz_1[i] * pa_y[i];

        g_yyy_xxzzz_0[i] = 2.0 * g_y_xxzzz_0[i] * fbe_0 - 2.0 * g_y_xxzzz_1[i] * fz_be_0 + g_yy_xxzzz_1[i] * pa_y[i];

        g_yyy_xyyyy_0[i] = 2.0 * g_y_xyyyy_0[i] * fbe_0 - 2.0 * g_y_xyyyy_1[i] * fz_be_0 + 4.0 * g_yy_xyyy_1[i] * fe_0 + g_yy_xyyyy_1[i] * pa_y[i];

        g_yyy_xyyyz_0[i] = 2.0 * g_y_xyyyz_0[i] * fbe_0 - 2.0 * g_y_xyyyz_1[i] * fz_be_0 + 3.0 * g_yy_xyyz_1[i] * fe_0 + g_yy_xyyyz_1[i] * pa_y[i];

        g_yyy_xyyzz_0[i] = 2.0 * g_y_xyyzz_0[i] * fbe_0 - 2.0 * g_y_xyyzz_1[i] * fz_be_0 + 2.0 * g_yy_xyzz_1[i] * fe_0 + g_yy_xyyzz_1[i] * pa_y[i];

        g_yyy_xyzzz_0[i] = 2.0 * g_y_xyzzz_0[i] * fbe_0 - 2.0 * g_y_xyzzz_1[i] * fz_be_0 + g_yy_xzzz_1[i] * fe_0 + g_yy_xyzzz_1[i] * pa_y[i];

        g_yyy_xzzzz_0[i] = 2.0 * g_y_xzzzz_0[i] * fbe_0 - 2.0 * g_y_xzzzz_1[i] * fz_be_0 + g_yy_xzzzz_1[i] * pa_y[i];

        g_yyy_yyyyy_0[i] = 2.0 * g_y_yyyyy_0[i] * fbe_0 - 2.0 * g_y_yyyyy_1[i] * fz_be_0 + 5.0 * g_yy_yyyy_1[i] * fe_0 + g_yy_yyyyy_1[i] * pa_y[i];

        g_yyy_yyyyz_0[i] = 2.0 * g_y_yyyyz_0[i] * fbe_0 - 2.0 * g_y_yyyyz_1[i] * fz_be_0 + 4.0 * g_yy_yyyz_1[i] * fe_0 + g_yy_yyyyz_1[i] * pa_y[i];

        g_yyy_yyyzz_0[i] = 2.0 * g_y_yyyzz_0[i] * fbe_0 - 2.0 * g_y_yyyzz_1[i] * fz_be_0 + 3.0 * g_yy_yyzz_1[i] * fe_0 + g_yy_yyyzz_1[i] * pa_y[i];

        g_yyy_yyzzz_0[i] = 2.0 * g_y_yyzzz_0[i] * fbe_0 - 2.0 * g_y_yyzzz_1[i] * fz_be_0 + 2.0 * g_yy_yzzz_1[i] * fe_0 + g_yy_yyzzz_1[i] * pa_y[i];

        g_yyy_yzzzz_0[i] = 2.0 * g_y_yzzzz_0[i] * fbe_0 - 2.0 * g_y_yzzzz_1[i] * fz_be_0 + g_yy_zzzz_1[i] * fe_0 + g_yy_yzzzz_1[i] * pa_y[i];

        g_yyy_zzzzz_0[i] = 2.0 * g_y_zzzzz_0[i] * fbe_0 - 2.0 * g_y_zzzzz_1[i] * fz_be_0 + g_yy_zzzzz_1[i] * pa_y[i];
    }

    // Set up 147-168 components of targeted buffer : FH

    auto g_yyz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 147);

    auto g_yyz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 148);

    auto g_yyz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 149);

    auto g_yyz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 150);

    auto g_yyz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 151);

    auto g_yyz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 152);

    auto g_yyz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 153);

    auto g_yyz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 154);

    auto g_yyz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 155);

    auto g_yyz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 156);

    auto g_yyz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 157);

    auto g_yyz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 158);

    auto g_yyz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 159);

    auto g_yyz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 160);

    auto g_yyz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 161);

    auto g_yyz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 162);

    auto g_yyz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 163);

    auto g_yyz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 164);

    auto g_yyz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 165);

    auto g_yyz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 166);

    auto g_yyz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 167);

    #pragma omp simd aligned(g_yy_xxxx_1, g_yy_xxxxx_1, g_yy_xxxxy_1, g_yy_xxxxz_1, g_yy_xxxy_1, g_yy_xxxyy_1, g_yy_xxxyz_1, g_yy_xxxz_1, g_yy_xxxzz_1, g_yy_xxyy_1, g_yy_xxyyy_1, g_yy_xxyyz_1, g_yy_xxyz_1, g_yy_xxyzz_1, g_yy_xxzz_1, g_yy_xxzzz_1, g_yy_xyyy_1, g_yy_xyyyy_1, g_yy_xyyyz_1, g_yy_xyyz_1, g_yy_xyyzz_1, g_yy_xyzz_1, g_yy_xyzzz_1, g_yy_xzzz_1, g_yy_xzzzz_1, g_yy_yyyy_1, g_yy_yyyyy_1, g_yy_yyyyz_1, g_yy_yyyz_1, g_yy_yyyzz_1, g_yy_yyzz_1, g_yy_yyzzz_1, g_yy_yzzz_1, g_yy_yzzzz_1, g_yy_zzzz_1, g_yy_zzzzz_1, g_yyz_xxxxx_0, g_yyz_xxxxy_0, g_yyz_xxxxz_0, g_yyz_xxxyy_0, g_yyz_xxxyz_0, g_yyz_xxxzz_0, g_yyz_xxyyy_0, g_yyz_xxyyz_0, g_yyz_xxyzz_0, g_yyz_xxzzz_0, g_yyz_xyyyy_0, g_yyz_xyyyz_0, g_yyz_xyyzz_0, g_yyz_xyzzz_0, g_yyz_xzzzz_0, g_yyz_yyyyy_0, g_yyz_yyyyz_0, g_yyz_yyyzz_0, g_yyz_yyzzz_0, g_yyz_yzzzz_0, g_yyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyz_xxxxx_0[i] = g_yy_xxxxx_1[i] * pa_z[i];

        g_yyz_xxxxy_0[i] = g_yy_xxxxy_1[i] * pa_z[i];

        g_yyz_xxxxz_0[i] = g_yy_xxxx_1[i] * fe_0 + g_yy_xxxxz_1[i] * pa_z[i];

        g_yyz_xxxyy_0[i] = g_yy_xxxyy_1[i] * pa_z[i];

        g_yyz_xxxyz_0[i] = g_yy_xxxy_1[i] * fe_0 + g_yy_xxxyz_1[i] * pa_z[i];

        g_yyz_xxxzz_0[i] = 2.0 * g_yy_xxxz_1[i] * fe_0 + g_yy_xxxzz_1[i] * pa_z[i];

        g_yyz_xxyyy_0[i] = g_yy_xxyyy_1[i] * pa_z[i];

        g_yyz_xxyyz_0[i] = g_yy_xxyy_1[i] * fe_0 + g_yy_xxyyz_1[i] * pa_z[i];

        g_yyz_xxyzz_0[i] = 2.0 * g_yy_xxyz_1[i] * fe_0 + g_yy_xxyzz_1[i] * pa_z[i];

        g_yyz_xxzzz_0[i] = 3.0 * g_yy_xxzz_1[i] * fe_0 + g_yy_xxzzz_1[i] * pa_z[i];

        g_yyz_xyyyy_0[i] = g_yy_xyyyy_1[i] * pa_z[i];

        g_yyz_xyyyz_0[i] = g_yy_xyyy_1[i] * fe_0 + g_yy_xyyyz_1[i] * pa_z[i];

        g_yyz_xyyzz_0[i] = 2.0 * g_yy_xyyz_1[i] * fe_0 + g_yy_xyyzz_1[i] * pa_z[i];

        g_yyz_xyzzz_0[i] = 3.0 * g_yy_xyzz_1[i] * fe_0 + g_yy_xyzzz_1[i] * pa_z[i];

        g_yyz_xzzzz_0[i] = 4.0 * g_yy_xzzz_1[i] * fe_0 + g_yy_xzzzz_1[i] * pa_z[i];

        g_yyz_yyyyy_0[i] = g_yy_yyyyy_1[i] * pa_z[i];

        g_yyz_yyyyz_0[i] = g_yy_yyyy_1[i] * fe_0 + g_yy_yyyyz_1[i] * pa_z[i];

        g_yyz_yyyzz_0[i] = 2.0 * g_yy_yyyz_1[i] * fe_0 + g_yy_yyyzz_1[i] * pa_z[i];

        g_yyz_yyzzz_0[i] = 3.0 * g_yy_yyzz_1[i] * fe_0 + g_yy_yyzzz_1[i] * pa_z[i];

        g_yyz_yzzzz_0[i] = 4.0 * g_yy_yzzz_1[i] * fe_0 + g_yy_yzzzz_1[i] * pa_z[i];

        g_yyz_zzzzz_0[i] = 5.0 * g_yy_zzzz_1[i] * fe_0 + g_yy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 168-189 components of targeted buffer : FH

    auto g_yzz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 168);

    auto g_yzz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 169);

    auto g_yzz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 170);

    auto g_yzz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 171);

    auto g_yzz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 172);

    auto g_yzz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 173);

    auto g_yzz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 174);

    auto g_yzz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 175);

    auto g_yzz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 176);

    auto g_yzz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 177);

    auto g_yzz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 178);

    auto g_yzz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 179);

    auto g_yzz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 180);

    auto g_yzz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 181);

    auto g_yzz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 182);

    auto g_yzz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 183);

    auto g_yzz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 184);

    auto g_yzz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 185);

    auto g_yzz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 186);

    auto g_yzz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 187);

    auto g_yzz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 188);

    #pragma omp simd aligned(g_yzz_xxxxx_0, g_yzz_xxxxy_0, g_yzz_xxxxz_0, g_yzz_xxxyy_0, g_yzz_xxxyz_0, g_yzz_xxxzz_0, g_yzz_xxyyy_0, g_yzz_xxyyz_0, g_yzz_xxyzz_0, g_yzz_xxzzz_0, g_yzz_xyyyy_0, g_yzz_xyyyz_0, g_yzz_xyyzz_0, g_yzz_xyzzz_0, g_yzz_xzzzz_0, g_yzz_yyyyy_0, g_yzz_yyyyz_0, g_yzz_yyyzz_0, g_yzz_yyzzz_0, g_yzz_yzzzz_0, g_yzz_zzzzz_0, g_zz_xxxx_1, g_zz_xxxxx_1, g_zz_xxxxy_1, g_zz_xxxxz_1, g_zz_xxxy_1, g_zz_xxxyy_1, g_zz_xxxyz_1, g_zz_xxxz_1, g_zz_xxxzz_1, g_zz_xxyy_1, g_zz_xxyyy_1, g_zz_xxyyz_1, g_zz_xxyz_1, g_zz_xxyzz_1, g_zz_xxzz_1, g_zz_xxzzz_1, g_zz_xyyy_1, g_zz_xyyyy_1, g_zz_xyyyz_1, g_zz_xyyz_1, g_zz_xyyzz_1, g_zz_xyzz_1, g_zz_xyzzz_1, g_zz_xzzz_1, g_zz_xzzzz_1, g_zz_yyyy_1, g_zz_yyyyy_1, g_zz_yyyyz_1, g_zz_yyyz_1, g_zz_yyyzz_1, g_zz_yyzz_1, g_zz_yyzzz_1, g_zz_yzzz_1, g_zz_yzzzz_1, g_zz_zzzz_1, g_zz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzz_xxxxx_0[i] = g_zz_xxxxx_1[i] * pa_y[i];

        g_yzz_xxxxy_0[i] = g_zz_xxxx_1[i] * fe_0 + g_zz_xxxxy_1[i] * pa_y[i];

        g_yzz_xxxxz_0[i] = g_zz_xxxxz_1[i] * pa_y[i];

        g_yzz_xxxyy_0[i] = 2.0 * g_zz_xxxy_1[i] * fe_0 + g_zz_xxxyy_1[i] * pa_y[i];

        g_yzz_xxxyz_0[i] = g_zz_xxxz_1[i] * fe_0 + g_zz_xxxyz_1[i] * pa_y[i];

        g_yzz_xxxzz_0[i] = g_zz_xxxzz_1[i] * pa_y[i];

        g_yzz_xxyyy_0[i] = 3.0 * g_zz_xxyy_1[i] * fe_0 + g_zz_xxyyy_1[i] * pa_y[i];

        g_yzz_xxyyz_0[i] = 2.0 * g_zz_xxyz_1[i] * fe_0 + g_zz_xxyyz_1[i] * pa_y[i];

        g_yzz_xxyzz_0[i] = g_zz_xxzz_1[i] * fe_0 + g_zz_xxyzz_1[i] * pa_y[i];

        g_yzz_xxzzz_0[i] = g_zz_xxzzz_1[i] * pa_y[i];

        g_yzz_xyyyy_0[i] = 4.0 * g_zz_xyyy_1[i] * fe_0 + g_zz_xyyyy_1[i] * pa_y[i];

        g_yzz_xyyyz_0[i] = 3.0 * g_zz_xyyz_1[i] * fe_0 + g_zz_xyyyz_1[i] * pa_y[i];

        g_yzz_xyyzz_0[i] = 2.0 * g_zz_xyzz_1[i] * fe_0 + g_zz_xyyzz_1[i] * pa_y[i];

        g_yzz_xyzzz_0[i] = g_zz_xzzz_1[i] * fe_0 + g_zz_xyzzz_1[i] * pa_y[i];

        g_yzz_xzzzz_0[i] = g_zz_xzzzz_1[i] * pa_y[i];

        g_yzz_yyyyy_0[i] = 5.0 * g_zz_yyyy_1[i] * fe_0 + g_zz_yyyyy_1[i] * pa_y[i];

        g_yzz_yyyyz_0[i] = 4.0 * g_zz_yyyz_1[i] * fe_0 + g_zz_yyyyz_1[i] * pa_y[i];

        g_yzz_yyyzz_0[i] = 3.0 * g_zz_yyzz_1[i] * fe_0 + g_zz_yyyzz_1[i] * pa_y[i];

        g_yzz_yyzzz_0[i] = 2.0 * g_zz_yzzz_1[i] * fe_0 + g_zz_yyzzz_1[i] * pa_y[i];

        g_yzz_yzzzz_0[i] = g_zz_zzzz_1[i] * fe_0 + g_zz_yzzzz_1[i] * pa_y[i];

        g_yzz_zzzzz_0[i] = g_zz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : FH

    auto g_zzz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 189);

    auto g_zzz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 190);

    auto g_zzz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 191);

    auto g_zzz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 192);

    auto g_zzz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 193);

    auto g_zzz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 194);

    auto g_zzz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 195);

    auto g_zzz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 196);

    auto g_zzz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 197);

    auto g_zzz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 198);

    auto g_zzz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 199);

    auto g_zzz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 200);

    auto g_zzz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 201);

    auto g_zzz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 202);

    auto g_zzz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 203);

    auto g_zzz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 204);

    auto g_zzz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 205);

    auto g_zzz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 206);

    auto g_zzz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 207);

    auto g_zzz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 208);

    auto g_zzz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 209);

    #pragma omp simd aligned(g_z_xxxxx_0, g_z_xxxxx_1, g_z_xxxxy_0, g_z_xxxxy_1, g_z_xxxxz_0, g_z_xxxxz_1, g_z_xxxyy_0, g_z_xxxyy_1, g_z_xxxyz_0, g_z_xxxyz_1, g_z_xxxzz_0, g_z_xxxzz_1, g_z_xxyyy_0, g_z_xxyyy_1, g_z_xxyyz_0, g_z_xxyyz_1, g_z_xxyzz_0, g_z_xxyzz_1, g_z_xxzzz_0, g_z_xxzzz_1, g_z_xyyyy_0, g_z_xyyyy_1, g_z_xyyyz_0, g_z_xyyyz_1, g_z_xyyzz_0, g_z_xyyzz_1, g_z_xyzzz_0, g_z_xyzzz_1, g_z_xzzzz_0, g_z_xzzzz_1, g_z_yyyyy_0, g_z_yyyyy_1, g_z_yyyyz_0, g_z_yyyyz_1, g_z_yyyzz_0, g_z_yyyzz_1, g_z_yyzzz_0, g_z_yyzzz_1, g_z_yzzzz_0, g_z_yzzzz_1, g_z_zzzzz_0, g_z_zzzzz_1, g_zz_xxxx_1, g_zz_xxxxx_1, g_zz_xxxxy_1, g_zz_xxxxz_1, g_zz_xxxy_1, g_zz_xxxyy_1, g_zz_xxxyz_1, g_zz_xxxz_1, g_zz_xxxzz_1, g_zz_xxyy_1, g_zz_xxyyy_1, g_zz_xxyyz_1, g_zz_xxyz_1, g_zz_xxyzz_1, g_zz_xxzz_1, g_zz_xxzzz_1, g_zz_xyyy_1, g_zz_xyyyy_1, g_zz_xyyyz_1, g_zz_xyyz_1, g_zz_xyyzz_1, g_zz_xyzz_1, g_zz_xyzzz_1, g_zz_xzzz_1, g_zz_xzzzz_1, g_zz_yyyy_1, g_zz_yyyyy_1, g_zz_yyyyz_1, g_zz_yyyz_1, g_zz_yyyzz_1, g_zz_yyzz_1, g_zz_yyzzz_1, g_zz_yzzz_1, g_zz_yzzzz_1, g_zz_zzzz_1, g_zz_zzzzz_1, g_zzz_xxxxx_0, g_zzz_xxxxy_0, g_zzz_xxxxz_0, g_zzz_xxxyy_0, g_zzz_xxxyz_0, g_zzz_xxxzz_0, g_zzz_xxyyy_0, g_zzz_xxyyz_0, g_zzz_xxyzz_0, g_zzz_xxzzz_0, g_zzz_xyyyy_0, g_zzz_xyyyz_0, g_zzz_xyyzz_0, g_zzz_xyzzz_0, g_zzz_xzzzz_0, g_zzz_yyyyy_0, g_zzz_yyyyz_0, g_zzz_yyyzz_0, g_zzz_yyzzz_0, g_zzz_yzzzz_0, g_zzz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzz_xxxxx_0[i] = 2.0 * g_z_xxxxx_0[i] * fbe_0 - 2.0 * g_z_xxxxx_1[i] * fz_be_0 + g_zz_xxxxx_1[i] * pa_z[i];

        g_zzz_xxxxy_0[i] = 2.0 * g_z_xxxxy_0[i] * fbe_0 - 2.0 * g_z_xxxxy_1[i] * fz_be_0 + g_zz_xxxxy_1[i] * pa_z[i];

        g_zzz_xxxxz_0[i] = 2.0 * g_z_xxxxz_0[i] * fbe_0 - 2.0 * g_z_xxxxz_1[i] * fz_be_0 + g_zz_xxxx_1[i] * fe_0 + g_zz_xxxxz_1[i] * pa_z[i];

        g_zzz_xxxyy_0[i] = 2.0 * g_z_xxxyy_0[i] * fbe_0 - 2.0 * g_z_xxxyy_1[i] * fz_be_0 + g_zz_xxxyy_1[i] * pa_z[i];

        g_zzz_xxxyz_0[i] = 2.0 * g_z_xxxyz_0[i] * fbe_0 - 2.0 * g_z_xxxyz_1[i] * fz_be_0 + g_zz_xxxy_1[i] * fe_0 + g_zz_xxxyz_1[i] * pa_z[i];

        g_zzz_xxxzz_0[i] = 2.0 * g_z_xxxzz_0[i] * fbe_0 - 2.0 * g_z_xxxzz_1[i] * fz_be_0 + 2.0 * g_zz_xxxz_1[i] * fe_0 + g_zz_xxxzz_1[i] * pa_z[i];

        g_zzz_xxyyy_0[i] = 2.0 * g_z_xxyyy_0[i] * fbe_0 - 2.0 * g_z_xxyyy_1[i] * fz_be_0 + g_zz_xxyyy_1[i] * pa_z[i];

        g_zzz_xxyyz_0[i] = 2.0 * g_z_xxyyz_0[i] * fbe_0 - 2.0 * g_z_xxyyz_1[i] * fz_be_0 + g_zz_xxyy_1[i] * fe_0 + g_zz_xxyyz_1[i] * pa_z[i];

        g_zzz_xxyzz_0[i] = 2.0 * g_z_xxyzz_0[i] * fbe_0 - 2.0 * g_z_xxyzz_1[i] * fz_be_0 + 2.0 * g_zz_xxyz_1[i] * fe_0 + g_zz_xxyzz_1[i] * pa_z[i];

        g_zzz_xxzzz_0[i] = 2.0 * g_z_xxzzz_0[i] * fbe_0 - 2.0 * g_z_xxzzz_1[i] * fz_be_0 + 3.0 * g_zz_xxzz_1[i] * fe_0 + g_zz_xxzzz_1[i] * pa_z[i];

        g_zzz_xyyyy_0[i] = 2.0 * g_z_xyyyy_0[i] * fbe_0 - 2.0 * g_z_xyyyy_1[i] * fz_be_0 + g_zz_xyyyy_1[i] * pa_z[i];

        g_zzz_xyyyz_0[i] = 2.0 * g_z_xyyyz_0[i] * fbe_0 - 2.0 * g_z_xyyyz_1[i] * fz_be_0 + g_zz_xyyy_1[i] * fe_0 + g_zz_xyyyz_1[i] * pa_z[i];

        g_zzz_xyyzz_0[i] = 2.0 * g_z_xyyzz_0[i] * fbe_0 - 2.0 * g_z_xyyzz_1[i] * fz_be_0 + 2.0 * g_zz_xyyz_1[i] * fe_0 + g_zz_xyyzz_1[i] * pa_z[i];

        g_zzz_xyzzz_0[i] = 2.0 * g_z_xyzzz_0[i] * fbe_0 - 2.0 * g_z_xyzzz_1[i] * fz_be_0 + 3.0 * g_zz_xyzz_1[i] * fe_0 + g_zz_xyzzz_1[i] * pa_z[i];

        g_zzz_xzzzz_0[i] = 2.0 * g_z_xzzzz_0[i] * fbe_0 - 2.0 * g_z_xzzzz_1[i] * fz_be_0 + 4.0 * g_zz_xzzz_1[i] * fe_0 + g_zz_xzzzz_1[i] * pa_z[i];

        g_zzz_yyyyy_0[i] = 2.0 * g_z_yyyyy_0[i] * fbe_0 - 2.0 * g_z_yyyyy_1[i] * fz_be_0 + g_zz_yyyyy_1[i] * pa_z[i];

        g_zzz_yyyyz_0[i] = 2.0 * g_z_yyyyz_0[i] * fbe_0 - 2.0 * g_z_yyyyz_1[i] * fz_be_0 + g_zz_yyyy_1[i] * fe_0 + g_zz_yyyyz_1[i] * pa_z[i];

        g_zzz_yyyzz_0[i] = 2.0 * g_z_yyyzz_0[i] * fbe_0 - 2.0 * g_z_yyyzz_1[i] * fz_be_0 + 2.0 * g_zz_yyyz_1[i] * fe_0 + g_zz_yyyzz_1[i] * pa_z[i];

        g_zzz_yyzzz_0[i] = 2.0 * g_z_yyzzz_0[i] * fbe_0 - 2.0 * g_z_yyzzz_1[i] * fz_be_0 + 3.0 * g_zz_yyzz_1[i] * fe_0 + g_zz_yyzzz_1[i] * pa_z[i];

        g_zzz_yzzzz_0[i] = 2.0 * g_z_yzzzz_0[i] * fbe_0 - 2.0 * g_z_yzzzz_1[i] * fz_be_0 + 4.0 * g_zz_yzzz_1[i] * fe_0 + g_zz_yzzzz_1[i] * pa_z[i];

        g_zzz_zzzzz_0[i] = 2.0 * g_z_zzzzz_0[i] * fbe_0 - 2.0 * g_z_zzzzz_1[i] * fz_be_0 + 5.0 * g_zz_zzzz_1[i] * fe_0 + g_zz_zzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

