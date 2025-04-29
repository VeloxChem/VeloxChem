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

#include "TwoCenterElectronRepulsionPrimRecGG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_gg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gg,
                                const size_t idx_eri_0_dg,
                                const size_t idx_eri_1_dg,
                                const size_t idx_eri_1_ff,
                                const size_t idx_eri_1_fg,
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

    // Set up components of auxiliary buffer : DG

    auto g_xx_xxxx_0 = pbuffer.data(idx_eri_0_dg);

    auto g_xx_xxxy_0 = pbuffer.data(idx_eri_0_dg + 1);

    auto g_xx_xxxz_0 = pbuffer.data(idx_eri_0_dg + 2);

    auto g_xx_xxyy_0 = pbuffer.data(idx_eri_0_dg + 3);

    auto g_xx_xxyz_0 = pbuffer.data(idx_eri_0_dg + 4);

    auto g_xx_xxzz_0 = pbuffer.data(idx_eri_0_dg + 5);

    auto g_xx_xyyy_0 = pbuffer.data(idx_eri_0_dg + 6);

    auto g_xx_xyyz_0 = pbuffer.data(idx_eri_0_dg + 7);

    auto g_xx_xyzz_0 = pbuffer.data(idx_eri_0_dg + 8);

    auto g_xx_xzzz_0 = pbuffer.data(idx_eri_0_dg + 9);

    auto g_xx_yyyy_0 = pbuffer.data(idx_eri_0_dg + 10);

    auto g_xx_yyyz_0 = pbuffer.data(idx_eri_0_dg + 11);

    auto g_xx_yyzz_0 = pbuffer.data(idx_eri_0_dg + 12);

    auto g_xx_yzzz_0 = pbuffer.data(idx_eri_0_dg + 13);

    auto g_xx_zzzz_0 = pbuffer.data(idx_eri_0_dg + 14);

    auto g_yy_xxxx_0 = pbuffer.data(idx_eri_0_dg + 45);

    auto g_yy_xxxy_0 = pbuffer.data(idx_eri_0_dg + 46);

    auto g_yy_xxxz_0 = pbuffer.data(idx_eri_0_dg + 47);

    auto g_yy_xxyy_0 = pbuffer.data(idx_eri_0_dg + 48);

    auto g_yy_xxyz_0 = pbuffer.data(idx_eri_0_dg + 49);

    auto g_yy_xxzz_0 = pbuffer.data(idx_eri_0_dg + 50);

    auto g_yy_xyyy_0 = pbuffer.data(idx_eri_0_dg + 51);

    auto g_yy_xyyz_0 = pbuffer.data(idx_eri_0_dg + 52);

    auto g_yy_xyzz_0 = pbuffer.data(idx_eri_0_dg + 53);

    auto g_yy_xzzz_0 = pbuffer.data(idx_eri_0_dg + 54);

    auto g_yy_yyyy_0 = pbuffer.data(idx_eri_0_dg + 55);

    auto g_yy_yyyz_0 = pbuffer.data(idx_eri_0_dg + 56);

    auto g_yy_yyzz_0 = pbuffer.data(idx_eri_0_dg + 57);

    auto g_yy_yzzz_0 = pbuffer.data(idx_eri_0_dg + 58);

    auto g_yy_zzzz_0 = pbuffer.data(idx_eri_0_dg + 59);

    auto g_zz_xxxx_0 = pbuffer.data(idx_eri_0_dg + 75);

    auto g_zz_xxxy_0 = pbuffer.data(idx_eri_0_dg + 76);

    auto g_zz_xxxz_0 = pbuffer.data(idx_eri_0_dg + 77);

    auto g_zz_xxyy_0 = pbuffer.data(idx_eri_0_dg + 78);

    auto g_zz_xxyz_0 = pbuffer.data(idx_eri_0_dg + 79);

    auto g_zz_xxzz_0 = pbuffer.data(idx_eri_0_dg + 80);

    auto g_zz_xyyy_0 = pbuffer.data(idx_eri_0_dg + 81);

    auto g_zz_xyyz_0 = pbuffer.data(idx_eri_0_dg + 82);

    auto g_zz_xyzz_0 = pbuffer.data(idx_eri_0_dg + 83);

    auto g_zz_xzzz_0 = pbuffer.data(idx_eri_0_dg + 84);

    auto g_zz_yyyy_0 = pbuffer.data(idx_eri_0_dg + 85);

    auto g_zz_yyyz_0 = pbuffer.data(idx_eri_0_dg + 86);

    auto g_zz_yyzz_0 = pbuffer.data(idx_eri_0_dg + 87);

    auto g_zz_yzzz_0 = pbuffer.data(idx_eri_0_dg + 88);

    auto g_zz_zzzz_0 = pbuffer.data(idx_eri_0_dg + 89);

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

    // Set up components of auxiliary buffer : FF

    auto g_xxx_xxx_1 = pbuffer.data(idx_eri_1_ff);

    auto g_xxx_xxy_1 = pbuffer.data(idx_eri_1_ff + 1);

    auto g_xxx_xxz_1 = pbuffer.data(idx_eri_1_ff + 2);

    auto g_xxx_xyy_1 = pbuffer.data(idx_eri_1_ff + 3);

    auto g_xxx_xyz_1 = pbuffer.data(idx_eri_1_ff + 4);

    auto g_xxx_xzz_1 = pbuffer.data(idx_eri_1_ff + 5);

    auto g_xxx_yyy_1 = pbuffer.data(idx_eri_1_ff + 6);

    auto g_xxx_yyz_1 = pbuffer.data(idx_eri_1_ff + 7);

    auto g_xxx_yzz_1 = pbuffer.data(idx_eri_1_ff + 8);

    auto g_xxx_zzz_1 = pbuffer.data(idx_eri_1_ff + 9);

    auto g_xxz_xxz_1 = pbuffer.data(idx_eri_1_ff + 22);

    auto g_xxz_xyz_1 = pbuffer.data(idx_eri_1_ff + 24);

    auto g_xxz_xzz_1 = pbuffer.data(idx_eri_1_ff + 25);

    auto g_xxz_yyz_1 = pbuffer.data(idx_eri_1_ff + 27);

    auto g_xxz_yzz_1 = pbuffer.data(idx_eri_1_ff + 28);

    auto g_xxz_zzz_1 = pbuffer.data(idx_eri_1_ff + 29);

    auto g_xyy_xxy_1 = pbuffer.data(idx_eri_1_ff + 31);

    auto g_xyy_xyy_1 = pbuffer.data(idx_eri_1_ff + 33);

    auto g_xyy_xyz_1 = pbuffer.data(idx_eri_1_ff + 34);

    auto g_xyy_yyy_1 = pbuffer.data(idx_eri_1_ff + 36);

    auto g_xyy_yyz_1 = pbuffer.data(idx_eri_1_ff + 37);

    auto g_xyy_yzz_1 = pbuffer.data(idx_eri_1_ff + 38);

    auto g_xzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 52);

    auto g_xzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 54);

    auto g_xzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 55);

    auto g_xzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 57);

    auto g_xzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 58);

    auto g_xzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 59);

    auto g_yyy_xxx_1 = pbuffer.data(idx_eri_1_ff + 60);

    auto g_yyy_xxy_1 = pbuffer.data(idx_eri_1_ff + 61);

    auto g_yyy_xxz_1 = pbuffer.data(idx_eri_1_ff + 62);

    auto g_yyy_xyy_1 = pbuffer.data(idx_eri_1_ff + 63);

    auto g_yyy_xyz_1 = pbuffer.data(idx_eri_1_ff + 64);

    auto g_yyy_xzz_1 = pbuffer.data(idx_eri_1_ff + 65);

    auto g_yyy_yyy_1 = pbuffer.data(idx_eri_1_ff + 66);

    auto g_yyy_yyz_1 = pbuffer.data(idx_eri_1_ff + 67);

    auto g_yyy_yzz_1 = pbuffer.data(idx_eri_1_ff + 68);

    auto g_yyy_zzz_1 = pbuffer.data(idx_eri_1_ff + 69);

    auto g_yyz_xxz_1 = pbuffer.data(idx_eri_1_ff + 72);

    auto g_yyz_xyz_1 = pbuffer.data(idx_eri_1_ff + 74);

    auto g_yyz_xzz_1 = pbuffer.data(idx_eri_1_ff + 75);

    auto g_yyz_yyz_1 = pbuffer.data(idx_eri_1_ff + 77);

    auto g_yyz_yzz_1 = pbuffer.data(idx_eri_1_ff + 78);

    auto g_yyz_zzz_1 = pbuffer.data(idx_eri_1_ff + 79);

    auto g_yzz_xxy_1 = pbuffer.data(idx_eri_1_ff + 81);

    auto g_yzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 82);

    auto g_yzz_xyy_1 = pbuffer.data(idx_eri_1_ff + 83);

    auto g_yzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 84);

    auto g_yzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 85);

    auto g_yzz_yyy_1 = pbuffer.data(idx_eri_1_ff + 86);

    auto g_yzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 87);

    auto g_yzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 88);

    auto g_yzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 89);

    auto g_zzz_xxx_1 = pbuffer.data(idx_eri_1_ff + 90);

    auto g_zzz_xxy_1 = pbuffer.data(idx_eri_1_ff + 91);

    auto g_zzz_xxz_1 = pbuffer.data(idx_eri_1_ff + 92);

    auto g_zzz_xyy_1 = pbuffer.data(idx_eri_1_ff + 93);

    auto g_zzz_xyz_1 = pbuffer.data(idx_eri_1_ff + 94);

    auto g_zzz_xzz_1 = pbuffer.data(idx_eri_1_ff + 95);

    auto g_zzz_yyy_1 = pbuffer.data(idx_eri_1_ff + 96);

    auto g_zzz_yyz_1 = pbuffer.data(idx_eri_1_ff + 97);

    auto g_zzz_yzz_1 = pbuffer.data(idx_eri_1_ff + 98);

    auto g_zzz_zzz_1 = pbuffer.data(idx_eri_1_ff + 99);

    // Set up components of auxiliary buffer : FG

    auto g_xxx_xxxx_1 = pbuffer.data(idx_eri_1_fg);

    auto g_xxx_xxxy_1 = pbuffer.data(idx_eri_1_fg + 1);

    auto g_xxx_xxxz_1 = pbuffer.data(idx_eri_1_fg + 2);

    auto g_xxx_xxyy_1 = pbuffer.data(idx_eri_1_fg + 3);

    auto g_xxx_xxyz_1 = pbuffer.data(idx_eri_1_fg + 4);

    auto g_xxx_xxzz_1 = pbuffer.data(idx_eri_1_fg + 5);

    auto g_xxx_xyyy_1 = pbuffer.data(idx_eri_1_fg + 6);

    auto g_xxx_xyyz_1 = pbuffer.data(idx_eri_1_fg + 7);

    auto g_xxx_xyzz_1 = pbuffer.data(idx_eri_1_fg + 8);

    auto g_xxx_xzzz_1 = pbuffer.data(idx_eri_1_fg + 9);

    auto g_xxx_yyyy_1 = pbuffer.data(idx_eri_1_fg + 10);

    auto g_xxx_yyyz_1 = pbuffer.data(idx_eri_1_fg + 11);

    auto g_xxx_yyzz_1 = pbuffer.data(idx_eri_1_fg + 12);

    auto g_xxx_yzzz_1 = pbuffer.data(idx_eri_1_fg + 13);

    auto g_xxx_zzzz_1 = pbuffer.data(idx_eri_1_fg + 14);

    auto g_xxy_xxxx_1 = pbuffer.data(idx_eri_1_fg + 15);

    auto g_xxy_xxxy_1 = pbuffer.data(idx_eri_1_fg + 16);

    auto g_xxy_xxxz_1 = pbuffer.data(idx_eri_1_fg + 17);

    auto g_xxy_xxyy_1 = pbuffer.data(idx_eri_1_fg + 18);

    auto g_xxy_xxzz_1 = pbuffer.data(idx_eri_1_fg + 20);

    auto g_xxy_xyyy_1 = pbuffer.data(idx_eri_1_fg + 21);

    auto g_xxy_xzzz_1 = pbuffer.data(idx_eri_1_fg + 24);

    auto g_xxy_yyyy_1 = pbuffer.data(idx_eri_1_fg + 25);

    auto g_xxz_xxxx_1 = pbuffer.data(idx_eri_1_fg + 30);

    auto g_xxz_xxxy_1 = pbuffer.data(idx_eri_1_fg + 31);

    auto g_xxz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 32);

    auto g_xxz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 33);

    auto g_xxz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 34);

    auto g_xxz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 35);

    auto g_xxz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 36);

    auto g_xxz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 37);

    auto g_xxz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 38);

    auto g_xxz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 39);

    auto g_xxz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 41);

    auto g_xxz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 42);

    auto g_xxz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 43);

    auto g_xxz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 44);

    auto g_xyy_xxxx_1 = pbuffer.data(idx_eri_1_fg + 45);

    auto g_xyy_xxxy_1 = pbuffer.data(idx_eri_1_fg + 46);

    auto g_xyy_xxyy_1 = pbuffer.data(idx_eri_1_fg + 48);

    auto g_xyy_xxyz_1 = pbuffer.data(idx_eri_1_fg + 49);

    auto g_xyy_xyyy_1 = pbuffer.data(idx_eri_1_fg + 51);

    auto g_xyy_xyyz_1 = pbuffer.data(idx_eri_1_fg + 52);

    auto g_xyy_xyzz_1 = pbuffer.data(idx_eri_1_fg + 53);

    auto g_xyy_yyyy_1 = pbuffer.data(idx_eri_1_fg + 55);

    auto g_xyy_yyyz_1 = pbuffer.data(idx_eri_1_fg + 56);

    auto g_xyy_yyzz_1 = pbuffer.data(idx_eri_1_fg + 57);

    auto g_xyy_yzzz_1 = pbuffer.data(idx_eri_1_fg + 58);

    auto g_xyy_zzzz_1 = pbuffer.data(idx_eri_1_fg + 59);

    auto g_xzz_xxxx_1 = pbuffer.data(idx_eri_1_fg + 75);

    auto g_xzz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 77);

    auto g_xzz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 79);

    auto g_xzz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 80);

    auto g_xzz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 82);

    auto g_xzz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 83);

    auto g_xzz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 84);

    auto g_xzz_yyyy_1 = pbuffer.data(idx_eri_1_fg + 85);

    auto g_xzz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 86);

    auto g_xzz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 87);

    auto g_xzz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 88);

    auto g_xzz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 89);

    auto g_yyy_xxxx_1 = pbuffer.data(idx_eri_1_fg + 90);

    auto g_yyy_xxxy_1 = pbuffer.data(idx_eri_1_fg + 91);

    auto g_yyy_xxxz_1 = pbuffer.data(idx_eri_1_fg + 92);

    auto g_yyy_xxyy_1 = pbuffer.data(idx_eri_1_fg + 93);

    auto g_yyy_xxyz_1 = pbuffer.data(idx_eri_1_fg + 94);

    auto g_yyy_xxzz_1 = pbuffer.data(idx_eri_1_fg + 95);

    auto g_yyy_xyyy_1 = pbuffer.data(idx_eri_1_fg + 96);

    auto g_yyy_xyyz_1 = pbuffer.data(idx_eri_1_fg + 97);

    auto g_yyy_xyzz_1 = pbuffer.data(idx_eri_1_fg + 98);

    auto g_yyy_xzzz_1 = pbuffer.data(idx_eri_1_fg + 99);

    auto g_yyy_yyyy_1 = pbuffer.data(idx_eri_1_fg + 100);

    auto g_yyy_yyyz_1 = pbuffer.data(idx_eri_1_fg + 101);

    auto g_yyy_yyzz_1 = pbuffer.data(idx_eri_1_fg + 102);

    auto g_yyy_yzzz_1 = pbuffer.data(idx_eri_1_fg + 103);

    auto g_yyy_zzzz_1 = pbuffer.data(idx_eri_1_fg + 104);

    auto g_yyz_xxxy_1 = pbuffer.data(idx_eri_1_fg + 106);

    auto g_yyz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 107);

    auto g_yyz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 108);

    auto g_yyz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 109);

    auto g_yyz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 110);

    auto g_yyz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 111);

    auto g_yyz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 112);

    auto g_yyz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 113);

    auto g_yyz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 114);

    auto g_yyz_yyyy_1 = pbuffer.data(idx_eri_1_fg + 115);

    auto g_yyz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 116);

    auto g_yyz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 117);

    auto g_yyz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 118);

    auto g_yyz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 119);

    auto g_yzz_xxxx_1 = pbuffer.data(idx_eri_1_fg + 120);

    auto g_yzz_xxxy_1 = pbuffer.data(idx_eri_1_fg + 121);

    auto g_yzz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 122);

    auto g_yzz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 123);

    auto g_yzz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 124);

    auto g_yzz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 125);

    auto g_yzz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 126);

    auto g_yzz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 127);

    auto g_yzz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 128);

    auto g_yzz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 129);

    auto g_yzz_yyyy_1 = pbuffer.data(idx_eri_1_fg + 130);

    auto g_yzz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 131);

    auto g_yzz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 132);

    auto g_yzz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 133);

    auto g_yzz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 134);

    auto g_zzz_xxxx_1 = pbuffer.data(idx_eri_1_fg + 135);

    auto g_zzz_xxxy_1 = pbuffer.data(idx_eri_1_fg + 136);

    auto g_zzz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 137);

    auto g_zzz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 138);

    auto g_zzz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 139);

    auto g_zzz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 140);

    auto g_zzz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 141);

    auto g_zzz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 142);

    auto g_zzz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 143);

    auto g_zzz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 144);

    auto g_zzz_yyyy_1 = pbuffer.data(idx_eri_1_fg + 145);

    auto g_zzz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 146);

    auto g_zzz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 147);

    auto g_zzz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 148);

    auto g_zzz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 149);

    // Set up 0-15 components of targeted buffer : GG

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

    #pragma omp simd aligned(g_xx_xxxx_0, g_xx_xxxx_1, g_xx_xxxy_0, g_xx_xxxy_1, g_xx_xxxz_0, g_xx_xxxz_1, g_xx_xxyy_0, g_xx_xxyy_1, g_xx_xxyz_0, g_xx_xxyz_1, g_xx_xxzz_0, g_xx_xxzz_1, g_xx_xyyy_0, g_xx_xyyy_1, g_xx_xyyz_0, g_xx_xyyz_1, g_xx_xyzz_0, g_xx_xyzz_1, g_xx_xzzz_0, g_xx_xzzz_1, g_xx_yyyy_0, g_xx_yyyy_1, g_xx_yyyz_0, g_xx_yyyz_1, g_xx_yyzz_0, g_xx_yyzz_1, g_xx_yzzz_0, g_xx_yzzz_1, g_xx_zzzz_0, g_xx_zzzz_1, g_xxx_xxx_1, g_xxx_xxxx_1, g_xxx_xxxy_1, g_xxx_xxxz_1, g_xxx_xxy_1, g_xxx_xxyy_1, g_xxx_xxyz_1, g_xxx_xxz_1, g_xxx_xxzz_1, g_xxx_xyy_1, g_xxx_xyyy_1, g_xxx_xyyz_1, g_xxx_xyz_1, g_xxx_xyzz_1, g_xxx_xzz_1, g_xxx_xzzz_1, g_xxx_yyy_1, g_xxx_yyyy_1, g_xxx_yyyz_1, g_xxx_yyz_1, g_xxx_yyzz_1, g_xxx_yzz_1, g_xxx_yzzz_1, g_xxx_zzz_1, g_xxx_zzzz_1, g_xxxx_xxxx_0, g_xxxx_xxxy_0, g_xxxx_xxxz_0, g_xxxx_xxyy_0, g_xxxx_xxyz_0, g_xxxx_xxzz_0, g_xxxx_xyyy_0, g_xxxx_xyyz_0, g_xxxx_xyzz_0, g_xxxx_xzzz_0, g_xxxx_yyyy_0, g_xxxx_yyyz_0, g_xxxx_yyzz_0, g_xxxx_yzzz_0, g_xxxx_zzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxx_xxxx_0[i] = 3.0 * g_xx_xxxx_0[i] * fbe_0 - 3.0 * g_xx_xxxx_1[i] * fz_be_0 + 4.0 * g_xxx_xxx_1[i] * fe_0 + g_xxx_xxxx_1[i] * pa_x[i];

        g_xxxx_xxxy_0[i] = 3.0 * g_xx_xxxy_0[i] * fbe_0 - 3.0 * g_xx_xxxy_1[i] * fz_be_0 + 3.0 * g_xxx_xxy_1[i] * fe_0 + g_xxx_xxxy_1[i] * pa_x[i];

        g_xxxx_xxxz_0[i] = 3.0 * g_xx_xxxz_0[i] * fbe_0 - 3.0 * g_xx_xxxz_1[i] * fz_be_0 + 3.0 * g_xxx_xxz_1[i] * fe_0 + g_xxx_xxxz_1[i] * pa_x[i];

        g_xxxx_xxyy_0[i] = 3.0 * g_xx_xxyy_0[i] * fbe_0 - 3.0 * g_xx_xxyy_1[i] * fz_be_0 + 2.0 * g_xxx_xyy_1[i] * fe_0 + g_xxx_xxyy_1[i] * pa_x[i];

        g_xxxx_xxyz_0[i] = 3.0 * g_xx_xxyz_0[i] * fbe_0 - 3.0 * g_xx_xxyz_1[i] * fz_be_0 + 2.0 * g_xxx_xyz_1[i] * fe_0 + g_xxx_xxyz_1[i] * pa_x[i];

        g_xxxx_xxzz_0[i] = 3.0 * g_xx_xxzz_0[i] * fbe_0 - 3.0 * g_xx_xxzz_1[i] * fz_be_0 + 2.0 * g_xxx_xzz_1[i] * fe_0 + g_xxx_xxzz_1[i] * pa_x[i];

        g_xxxx_xyyy_0[i] = 3.0 * g_xx_xyyy_0[i] * fbe_0 - 3.0 * g_xx_xyyy_1[i] * fz_be_0 + g_xxx_yyy_1[i] * fe_0 + g_xxx_xyyy_1[i] * pa_x[i];

        g_xxxx_xyyz_0[i] = 3.0 * g_xx_xyyz_0[i] * fbe_0 - 3.0 * g_xx_xyyz_1[i] * fz_be_0 + g_xxx_yyz_1[i] * fe_0 + g_xxx_xyyz_1[i] * pa_x[i];

        g_xxxx_xyzz_0[i] = 3.0 * g_xx_xyzz_0[i] * fbe_0 - 3.0 * g_xx_xyzz_1[i] * fz_be_0 + g_xxx_yzz_1[i] * fe_0 + g_xxx_xyzz_1[i] * pa_x[i];

        g_xxxx_xzzz_0[i] = 3.0 * g_xx_xzzz_0[i] * fbe_0 - 3.0 * g_xx_xzzz_1[i] * fz_be_0 + g_xxx_zzz_1[i] * fe_0 + g_xxx_xzzz_1[i] * pa_x[i];

        g_xxxx_yyyy_0[i] = 3.0 * g_xx_yyyy_0[i] * fbe_0 - 3.0 * g_xx_yyyy_1[i] * fz_be_0 + g_xxx_yyyy_1[i] * pa_x[i];

        g_xxxx_yyyz_0[i] = 3.0 * g_xx_yyyz_0[i] * fbe_0 - 3.0 * g_xx_yyyz_1[i] * fz_be_0 + g_xxx_yyyz_1[i] * pa_x[i];

        g_xxxx_yyzz_0[i] = 3.0 * g_xx_yyzz_0[i] * fbe_0 - 3.0 * g_xx_yyzz_1[i] * fz_be_0 + g_xxx_yyzz_1[i] * pa_x[i];

        g_xxxx_yzzz_0[i] = 3.0 * g_xx_yzzz_0[i] * fbe_0 - 3.0 * g_xx_yzzz_1[i] * fz_be_0 + g_xxx_yzzz_1[i] * pa_x[i];

        g_xxxx_zzzz_0[i] = 3.0 * g_xx_zzzz_0[i] * fbe_0 - 3.0 * g_xx_zzzz_1[i] * fz_be_0 + g_xxx_zzzz_1[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : GG

    auto g_xxxy_xxxx_0 = pbuffer.data(idx_eri_0_gg + 15);

    auto g_xxxy_xxxy_0 = pbuffer.data(idx_eri_0_gg + 16);

    auto g_xxxy_xxxz_0 = pbuffer.data(idx_eri_0_gg + 17);

    auto g_xxxy_xxyy_0 = pbuffer.data(idx_eri_0_gg + 18);

    auto g_xxxy_xxyz_0 = pbuffer.data(idx_eri_0_gg + 19);

    auto g_xxxy_xxzz_0 = pbuffer.data(idx_eri_0_gg + 20);

    auto g_xxxy_xyyy_0 = pbuffer.data(idx_eri_0_gg + 21);

    auto g_xxxy_xyyz_0 = pbuffer.data(idx_eri_0_gg + 22);

    auto g_xxxy_xyzz_0 = pbuffer.data(idx_eri_0_gg + 23);

    auto g_xxxy_xzzz_0 = pbuffer.data(idx_eri_0_gg + 24);

    auto g_xxxy_yyyy_0 = pbuffer.data(idx_eri_0_gg + 25);

    auto g_xxxy_yyyz_0 = pbuffer.data(idx_eri_0_gg + 26);

    auto g_xxxy_yyzz_0 = pbuffer.data(idx_eri_0_gg + 27);

    auto g_xxxy_yzzz_0 = pbuffer.data(idx_eri_0_gg + 28);

    auto g_xxxy_zzzz_0 = pbuffer.data(idx_eri_0_gg + 29);

    #pragma omp simd aligned(g_xxx_xxx_1, g_xxx_xxxx_1, g_xxx_xxxy_1, g_xxx_xxxz_1, g_xxx_xxy_1, g_xxx_xxyy_1, g_xxx_xxyz_1, g_xxx_xxz_1, g_xxx_xxzz_1, g_xxx_xyy_1, g_xxx_xyyy_1, g_xxx_xyyz_1, g_xxx_xyz_1, g_xxx_xyzz_1, g_xxx_xzz_1, g_xxx_xzzz_1, g_xxx_yyy_1, g_xxx_yyyy_1, g_xxx_yyyz_1, g_xxx_yyz_1, g_xxx_yyzz_1, g_xxx_yzz_1, g_xxx_yzzz_1, g_xxx_zzz_1, g_xxx_zzzz_1, g_xxxy_xxxx_0, g_xxxy_xxxy_0, g_xxxy_xxxz_0, g_xxxy_xxyy_0, g_xxxy_xxyz_0, g_xxxy_xxzz_0, g_xxxy_xyyy_0, g_xxxy_xyyz_0, g_xxxy_xyzz_0, g_xxxy_xzzz_0, g_xxxy_yyyy_0, g_xxxy_yyyz_0, g_xxxy_yyzz_0, g_xxxy_yzzz_0, g_xxxy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxy_xxxx_0[i] = g_xxx_xxxx_1[i] * pa_y[i];

        g_xxxy_xxxy_0[i] = g_xxx_xxx_1[i] * fe_0 + g_xxx_xxxy_1[i] * pa_y[i];

        g_xxxy_xxxz_0[i] = g_xxx_xxxz_1[i] * pa_y[i];

        g_xxxy_xxyy_0[i] = 2.0 * g_xxx_xxy_1[i] * fe_0 + g_xxx_xxyy_1[i] * pa_y[i];

        g_xxxy_xxyz_0[i] = g_xxx_xxz_1[i] * fe_0 + g_xxx_xxyz_1[i] * pa_y[i];

        g_xxxy_xxzz_0[i] = g_xxx_xxzz_1[i] * pa_y[i];

        g_xxxy_xyyy_0[i] = 3.0 * g_xxx_xyy_1[i] * fe_0 + g_xxx_xyyy_1[i] * pa_y[i];

        g_xxxy_xyyz_0[i] = 2.0 * g_xxx_xyz_1[i] * fe_0 + g_xxx_xyyz_1[i] * pa_y[i];

        g_xxxy_xyzz_0[i] = g_xxx_xzz_1[i] * fe_0 + g_xxx_xyzz_1[i] * pa_y[i];

        g_xxxy_xzzz_0[i] = g_xxx_xzzz_1[i] * pa_y[i];

        g_xxxy_yyyy_0[i] = 4.0 * g_xxx_yyy_1[i] * fe_0 + g_xxx_yyyy_1[i] * pa_y[i];

        g_xxxy_yyyz_0[i] = 3.0 * g_xxx_yyz_1[i] * fe_0 + g_xxx_yyyz_1[i] * pa_y[i];

        g_xxxy_yyzz_0[i] = 2.0 * g_xxx_yzz_1[i] * fe_0 + g_xxx_yyzz_1[i] * pa_y[i];

        g_xxxy_yzzz_0[i] = g_xxx_zzz_1[i] * fe_0 + g_xxx_yzzz_1[i] * pa_y[i];

        g_xxxy_zzzz_0[i] = g_xxx_zzzz_1[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : GG

    auto g_xxxz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 30);

    auto g_xxxz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 31);

    auto g_xxxz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 32);

    auto g_xxxz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 33);

    auto g_xxxz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 34);

    auto g_xxxz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 35);

    auto g_xxxz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 36);

    auto g_xxxz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 37);

    auto g_xxxz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 38);

    auto g_xxxz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 39);

    auto g_xxxz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 40);

    auto g_xxxz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 41);

    auto g_xxxz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 42);

    auto g_xxxz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 43);

    auto g_xxxz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 44);

    #pragma omp simd aligned(g_xxx_xxx_1, g_xxx_xxxx_1, g_xxx_xxxy_1, g_xxx_xxxz_1, g_xxx_xxy_1, g_xxx_xxyy_1, g_xxx_xxyz_1, g_xxx_xxz_1, g_xxx_xxzz_1, g_xxx_xyy_1, g_xxx_xyyy_1, g_xxx_xyyz_1, g_xxx_xyz_1, g_xxx_xyzz_1, g_xxx_xzz_1, g_xxx_xzzz_1, g_xxx_yyy_1, g_xxx_yyyy_1, g_xxx_yyyz_1, g_xxx_yyz_1, g_xxx_yyzz_1, g_xxx_yzz_1, g_xxx_yzzz_1, g_xxx_zzz_1, g_xxx_zzzz_1, g_xxxz_xxxx_0, g_xxxz_xxxy_0, g_xxxz_xxxz_0, g_xxxz_xxyy_0, g_xxxz_xxyz_0, g_xxxz_xxzz_0, g_xxxz_xyyy_0, g_xxxz_xyyz_0, g_xxxz_xyzz_0, g_xxxz_xzzz_0, g_xxxz_yyyy_0, g_xxxz_yyyz_0, g_xxxz_yyzz_0, g_xxxz_yzzz_0, g_xxxz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxz_xxxx_0[i] = g_xxx_xxxx_1[i] * pa_z[i];

        g_xxxz_xxxy_0[i] = g_xxx_xxxy_1[i] * pa_z[i];

        g_xxxz_xxxz_0[i] = g_xxx_xxx_1[i] * fe_0 + g_xxx_xxxz_1[i] * pa_z[i];

        g_xxxz_xxyy_0[i] = g_xxx_xxyy_1[i] * pa_z[i];

        g_xxxz_xxyz_0[i] = g_xxx_xxy_1[i] * fe_0 + g_xxx_xxyz_1[i] * pa_z[i];

        g_xxxz_xxzz_0[i] = 2.0 * g_xxx_xxz_1[i] * fe_0 + g_xxx_xxzz_1[i] * pa_z[i];

        g_xxxz_xyyy_0[i] = g_xxx_xyyy_1[i] * pa_z[i];

        g_xxxz_xyyz_0[i] = g_xxx_xyy_1[i] * fe_0 + g_xxx_xyyz_1[i] * pa_z[i];

        g_xxxz_xyzz_0[i] = 2.0 * g_xxx_xyz_1[i] * fe_0 + g_xxx_xyzz_1[i] * pa_z[i];

        g_xxxz_xzzz_0[i] = 3.0 * g_xxx_xzz_1[i] * fe_0 + g_xxx_xzzz_1[i] * pa_z[i];

        g_xxxz_yyyy_0[i] = g_xxx_yyyy_1[i] * pa_z[i];

        g_xxxz_yyyz_0[i] = g_xxx_yyy_1[i] * fe_0 + g_xxx_yyyz_1[i] * pa_z[i];

        g_xxxz_yyzz_0[i] = 2.0 * g_xxx_yyz_1[i] * fe_0 + g_xxx_yyzz_1[i] * pa_z[i];

        g_xxxz_yzzz_0[i] = 3.0 * g_xxx_yzz_1[i] * fe_0 + g_xxx_yzzz_1[i] * pa_z[i];

        g_xxxz_zzzz_0[i] = 4.0 * g_xxx_zzz_1[i] * fe_0 + g_xxx_zzzz_1[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : GG

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

    #pragma omp simd aligned(g_xx_xxxx_0, g_xx_xxxx_1, g_xx_xxxz_0, g_xx_xxxz_1, g_xx_xxzz_0, g_xx_xxzz_1, g_xx_xzzz_0, g_xx_xzzz_1, g_xxy_xxxx_1, g_xxy_xxxz_1, g_xxy_xxzz_1, g_xxy_xzzz_1, g_xxyy_xxxx_0, g_xxyy_xxxy_0, g_xxyy_xxxz_0, g_xxyy_xxyy_0, g_xxyy_xxyz_0, g_xxyy_xxzz_0, g_xxyy_xyyy_0, g_xxyy_xyyz_0, g_xxyy_xyzz_0, g_xxyy_xzzz_0, g_xxyy_yyyy_0, g_xxyy_yyyz_0, g_xxyy_yyzz_0, g_xxyy_yzzz_0, g_xxyy_zzzz_0, g_xyy_xxxy_1, g_xyy_xxy_1, g_xyy_xxyy_1, g_xyy_xxyz_1, g_xyy_xyy_1, g_xyy_xyyy_1, g_xyy_xyyz_1, g_xyy_xyz_1, g_xyy_xyzz_1, g_xyy_yyy_1, g_xyy_yyyy_1, g_xyy_yyyz_1, g_xyy_yyz_1, g_xyy_yyzz_1, g_xyy_yzz_1, g_xyy_yzzz_1, g_xyy_zzzz_1, g_yy_xxxy_0, g_yy_xxxy_1, g_yy_xxyy_0, g_yy_xxyy_1, g_yy_xxyz_0, g_yy_xxyz_1, g_yy_xyyy_0, g_yy_xyyy_1, g_yy_xyyz_0, g_yy_xyyz_1, g_yy_xyzz_0, g_yy_xyzz_1, g_yy_yyyy_0, g_yy_yyyy_1, g_yy_yyyz_0, g_yy_yyyz_1, g_yy_yyzz_0, g_yy_yyzz_1, g_yy_yzzz_0, g_yy_yzzz_1, g_yy_zzzz_0, g_yy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyy_xxxx_0[i] = g_xx_xxxx_0[i] * fbe_0 - g_xx_xxxx_1[i] * fz_be_0 + g_xxy_xxxx_1[i] * pa_y[i];

        g_xxyy_xxxy_0[i] = g_yy_xxxy_0[i] * fbe_0 - g_yy_xxxy_1[i] * fz_be_0 + 3.0 * g_xyy_xxy_1[i] * fe_0 + g_xyy_xxxy_1[i] * pa_x[i];

        g_xxyy_xxxz_0[i] = g_xx_xxxz_0[i] * fbe_0 - g_xx_xxxz_1[i] * fz_be_0 + g_xxy_xxxz_1[i] * pa_y[i];

        g_xxyy_xxyy_0[i] = g_yy_xxyy_0[i] * fbe_0 - g_yy_xxyy_1[i] * fz_be_0 + 2.0 * g_xyy_xyy_1[i] * fe_0 + g_xyy_xxyy_1[i] * pa_x[i];

        g_xxyy_xxyz_0[i] = g_yy_xxyz_0[i] * fbe_0 - g_yy_xxyz_1[i] * fz_be_0 + 2.0 * g_xyy_xyz_1[i] * fe_0 + g_xyy_xxyz_1[i] * pa_x[i];

        g_xxyy_xxzz_0[i] = g_xx_xxzz_0[i] * fbe_0 - g_xx_xxzz_1[i] * fz_be_0 + g_xxy_xxzz_1[i] * pa_y[i];

        g_xxyy_xyyy_0[i] = g_yy_xyyy_0[i] * fbe_0 - g_yy_xyyy_1[i] * fz_be_0 + g_xyy_yyy_1[i] * fe_0 + g_xyy_xyyy_1[i] * pa_x[i];

        g_xxyy_xyyz_0[i] = g_yy_xyyz_0[i] * fbe_0 - g_yy_xyyz_1[i] * fz_be_0 + g_xyy_yyz_1[i] * fe_0 + g_xyy_xyyz_1[i] * pa_x[i];

        g_xxyy_xyzz_0[i] = g_yy_xyzz_0[i] * fbe_0 - g_yy_xyzz_1[i] * fz_be_0 + g_xyy_yzz_1[i] * fe_0 + g_xyy_xyzz_1[i] * pa_x[i];

        g_xxyy_xzzz_0[i] = g_xx_xzzz_0[i] * fbe_0 - g_xx_xzzz_1[i] * fz_be_0 + g_xxy_xzzz_1[i] * pa_y[i];

        g_xxyy_yyyy_0[i] = g_yy_yyyy_0[i] * fbe_0 - g_yy_yyyy_1[i] * fz_be_0 + g_xyy_yyyy_1[i] * pa_x[i];

        g_xxyy_yyyz_0[i] = g_yy_yyyz_0[i] * fbe_0 - g_yy_yyyz_1[i] * fz_be_0 + g_xyy_yyyz_1[i] * pa_x[i];

        g_xxyy_yyzz_0[i] = g_yy_yyzz_0[i] * fbe_0 - g_yy_yyzz_1[i] * fz_be_0 + g_xyy_yyzz_1[i] * pa_x[i];

        g_xxyy_yzzz_0[i] = g_yy_yzzz_0[i] * fbe_0 - g_yy_yzzz_1[i] * fz_be_0 + g_xyy_yzzz_1[i] * pa_x[i];

        g_xxyy_zzzz_0[i] = g_yy_zzzz_0[i] * fbe_0 - g_yy_zzzz_1[i] * fz_be_0 + g_xyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : GG

    auto g_xxyz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 60);

    auto g_xxyz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 61);

    auto g_xxyz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 62);

    auto g_xxyz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 63);

    auto g_xxyz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 64);

    auto g_xxyz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 65);

    auto g_xxyz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 66);

    auto g_xxyz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 67);

    auto g_xxyz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 68);

    auto g_xxyz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 69);

    auto g_xxyz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 70);

    auto g_xxyz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 71);

    auto g_xxyz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 72);

    auto g_xxyz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 73);

    auto g_xxyz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 74);

    #pragma omp simd aligned(g_xxy_xxxy_1, g_xxy_xxyy_1, g_xxy_xyyy_1, g_xxy_yyyy_1, g_xxyz_xxxx_0, g_xxyz_xxxy_0, g_xxyz_xxxz_0, g_xxyz_xxyy_0, g_xxyz_xxyz_0, g_xxyz_xxzz_0, g_xxyz_xyyy_0, g_xxyz_xyyz_0, g_xxyz_xyzz_0, g_xxyz_xzzz_0, g_xxyz_yyyy_0, g_xxyz_yyyz_0, g_xxyz_yyzz_0, g_xxyz_yzzz_0, g_xxyz_zzzz_0, g_xxz_xxxx_1, g_xxz_xxxz_1, g_xxz_xxyz_1, g_xxz_xxz_1, g_xxz_xxzz_1, g_xxz_xyyz_1, g_xxz_xyz_1, g_xxz_xyzz_1, g_xxz_xzz_1, g_xxz_xzzz_1, g_xxz_yyyz_1, g_xxz_yyz_1, g_xxz_yyzz_1, g_xxz_yzz_1, g_xxz_yzzz_1, g_xxz_zzz_1, g_xxz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyz_xxxx_0[i] = g_xxz_xxxx_1[i] * pa_y[i];

        g_xxyz_xxxy_0[i] = g_xxy_xxxy_1[i] * pa_z[i];

        g_xxyz_xxxz_0[i] = g_xxz_xxxz_1[i] * pa_y[i];

        g_xxyz_xxyy_0[i] = g_xxy_xxyy_1[i] * pa_z[i];

        g_xxyz_xxyz_0[i] = g_xxz_xxz_1[i] * fe_0 + g_xxz_xxyz_1[i] * pa_y[i];

        g_xxyz_xxzz_0[i] = g_xxz_xxzz_1[i] * pa_y[i];

        g_xxyz_xyyy_0[i] = g_xxy_xyyy_1[i] * pa_z[i];

        g_xxyz_xyyz_0[i] = 2.0 * g_xxz_xyz_1[i] * fe_0 + g_xxz_xyyz_1[i] * pa_y[i];

        g_xxyz_xyzz_0[i] = g_xxz_xzz_1[i] * fe_0 + g_xxz_xyzz_1[i] * pa_y[i];

        g_xxyz_xzzz_0[i] = g_xxz_xzzz_1[i] * pa_y[i];

        g_xxyz_yyyy_0[i] = g_xxy_yyyy_1[i] * pa_z[i];

        g_xxyz_yyyz_0[i] = 3.0 * g_xxz_yyz_1[i] * fe_0 + g_xxz_yyyz_1[i] * pa_y[i];

        g_xxyz_yyzz_0[i] = 2.0 * g_xxz_yzz_1[i] * fe_0 + g_xxz_yyzz_1[i] * pa_y[i];

        g_xxyz_yzzz_0[i] = g_xxz_zzz_1[i] * fe_0 + g_xxz_yzzz_1[i] * pa_y[i];

        g_xxyz_zzzz_0[i] = g_xxz_zzzz_1[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : GG

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

    #pragma omp simd aligned(g_xx_xxxx_0, g_xx_xxxx_1, g_xx_xxxy_0, g_xx_xxxy_1, g_xx_xxyy_0, g_xx_xxyy_1, g_xx_xyyy_0, g_xx_xyyy_1, g_xxz_xxxx_1, g_xxz_xxxy_1, g_xxz_xxyy_1, g_xxz_xyyy_1, g_xxzz_xxxx_0, g_xxzz_xxxy_0, g_xxzz_xxxz_0, g_xxzz_xxyy_0, g_xxzz_xxyz_0, g_xxzz_xxzz_0, g_xxzz_xyyy_0, g_xxzz_xyyz_0, g_xxzz_xyzz_0, g_xxzz_xzzz_0, g_xxzz_yyyy_0, g_xxzz_yyyz_0, g_xxzz_yyzz_0, g_xxzz_yzzz_0, g_xxzz_zzzz_0, g_xzz_xxxz_1, g_xzz_xxyz_1, g_xzz_xxz_1, g_xzz_xxzz_1, g_xzz_xyyz_1, g_xzz_xyz_1, g_xzz_xyzz_1, g_xzz_xzz_1, g_xzz_xzzz_1, g_xzz_yyyy_1, g_xzz_yyyz_1, g_xzz_yyz_1, g_xzz_yyzz_1, g_xzz_yzz_1, g_xzz_yzzz_1, g_xzz_zzz_1, g_xzz_zzzz_1, g_zz_xxxz_0, g_zz_xxxz_1, g_zz_xxyz_0, g_zz_xxyz_1, g_zz_xxzz_0, g_zz_xxzz_1, g_zz_xyyz_0, g_zz_xyyz_1, g_zz_xyzz_0, g_zz_xyzz_1, g_zz_xzzz_0, g_zz_xzzz_1, g_zz_yyyy_0, g_zz_yyyy_1, g_zz_yyyz_0, g_zz_yyyz_1, g_zz_yyzz_0, g_zz_yyzz_1, g_zz_yzzz_0, g_zz_yzzz_1, g_zz_zzzz_0, g_zz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzz_xxxx_0[i] = g_xx_xxxx_0[i] * fbe_0 - g_xx_xxxx_1[i] * fz_be_0 + g_xxz_xxxx_1[i] * pa_z[i];

        g_xxzz_xxxy_0[i] = g_xx_xxxy_0[i] * fbe_0 - g_xx_xxxy_1[i] * fz_be_0 + g_xxz_xxxy_1[i] * pa_z[i];

        g_xxzz_xxxz_0[i] = g_zz_xxxz_0[i] * fbe_0 - g_zz_xxxz_1[i] * fz_be_0 + 3.0 * g_xzz_xxz_1[i] * fe_0 + g_xzz_xxxz_1[i] * pa_x[i];

        g_xxzz_xxyy_0[i] = g_xx_xxyy_0[i] * fbe_0 - g_xx_xxyy_1[i] * fz_be_0 + g_xxz_xxyy_1[i] * pa_z[i];

        g_xxzz_xxyz_0[i] = g_zz_xxyz_0[i] * fbe_0 - g_zz_xxyz_1[i] * fz_be_0 + 2.0 * g_xzz_xyz_1[i] * fe_0 + g_xzz_xxyz_1[i] * pa_x[i];

        g_xxzz_xxzz_0[i] = g_zz_xxzz_0[i] * fbe_0 - g_zz_xxzz_1[i] * fz_be_0 + 2.0 * g_xzz_xzz_1[i] * fe_0 + g_xzz_xxzz_1[i] * pa_x[i];

        g_xxzz_xyyy_0[i] = g_xx_xyyy_0[i] * fbe_0 - g_xx_xyyy_1[i] * fz_be_0 + g_xxz_xyyy_1[i] * pa_z[i];

        g_xxzz_xyyz_0[i] = g_zz_xyyz_0[i] * fbe_0 - g_zz_xyyz_1[i] * fz_be_0 + g_xzz_yyz_1[i] * fe_0 + g_xzz_xyyz_1[i] * pa_x[i];

        g_xxzz_xyzz_0[i] = g_zz_xyzz_0[i] * fbe_0 - g_zz_xyzz_1[i] * fz_be_0 + g_xzz_yzz_1[i] * fe_0 + g_xzz_xyzz_1[i] * pa_x[i];

        g_xxzz_xzzz_0[i] = g_zz_xzzz_0[i] * fbe_0 - g_zz_xzzz_1[i] * fz_be_0 + g_xzz_zzz_1[i] * fe_0 + g_xzz_xzzz_1[i] * pa_x[i];

        g_xxzz_yyyy_0[i] = g_zz_yyyy_0[i] * fbe_0 - g_zz_yyyy_1[i] * fz_be_0 + g_xzz_yyyy_1[i] * pa_x[i];

        g_xxzz_yyyz_0[i] = g_zz_yyyz_0[i] * fbe_0 - g_zz_yyyz_1[i] * fz_be_0 + g_xzz_yyyz_1[i] * pa_x[i];

        g_xxzz_yyzz_0[i] = g_zz_yyzz_0[i] * fbe_0 - g_zz_yyzz_1[i] * fz_be_0 + g_xzz_yyzz_1[i] * pa_x[i];

        g_xxzz_yzzz_0[i] = g_zz_yzzz_0[i] * fbe_0 - g_zz_yzzz_1[i] * fz_be_0 + g_xzz_yzzz_1[i] * pa_x[i];

        g_xxzz_zzzz_0[i] = g_zz_zzzz_0[i] * fbe_0 - g_zz_zzzz_1[i] * fz_be_0 + g_xzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : GG

    auto g_xyyy_xxxx_0 = pbuffer.data(idx_eri_0_gg + 90);

    auto g_xyyy_xxxy_0 = pbuffer.data(idx_eri_0_gg + 91);

    auto g_xyyy_xxxz_0 = pbuffer.data(idx_eri_0_gg + 92);

    auto g_xyyy_xxyy_0 = pbuffer.data(idx_eri_0_gg + 93);

    auto g_xyyy_xxyz_0 = pbuffer.data(idx_eri_0_gg + 94);

    auto g_xyyy_xxzz_0 = pbuffer.data(idx_eri_0_gg + 95);

    auto g_xyyy_xyyy_0 = pbuffer.data(idx_eri_0_gg + 96);

    auto g_xyyy_xyyz_0 = pbuffer.data(idx_eri_0_gg + 97);

    auto g_xyyy_xyzz_0 = pbuffer.data(idx_eri_0_gg + 98);

    auto g_xyyy_xzzz_0 = pbuffer.data(idx_eri_0_gg + 99);

    auto g_xyyy_yyyy_0 = pbuffer.data(idx_eri_0_gg + 100);

    auto g_xyyy_yyyz_0 = pbuffer.data(idx_eri_0_gg + 101);

    auto g_xyyy_yyzz_0 = pbuffer.data(idx_eri_0_gg + 102);

    auto g_xyyy_yzzz_0 = pbuffer.data(idx_eri_0_gg + 103);

    auto g_xyyy_zzzz_0 = pbuffer.data(idx_eri_0_gg + 104);

    #pragma omp simd aligned(g_xyyy_xxxx_0, g_xyyy_xxxy_0, g_xyyy_xxxz_0, g_xyyy_xxyy_0, g_xyyy_xxyz_0, g_xyyy_xxzz_0, g_xyyy_xyyy_0, g_xyyy_xyyz_0, g_xyyy_xyzz_0, g_xyyy_xzzz_0, g_xyyy_yyyy_0, g_xyyy_yyyz_0, g_xyyy_yyzz_0, g_xyyy_yzzz_0, g_xyyy_zzzz_0, g_yyy_xxx_1, g_yyy_xxxx_1, g_yyy_xxxy_1, g_yyy_xxxz_1, g_yyy_xxy_1, g_yyy_xxyy_1, g_yyy_xxyz_1, g_yyy_xxz_1, g_yyy_xxzz_1, g_yyy_xyy_1, g_yyy_xyyy_1, g_yyy_xyyz_1, g_yyy_xyz_1, g_yyy_xyzz_1, g_yyy_xzz_1, g_yyy_xzzz_1, g_yyy_yyy_1, g_yyy_yyyy_1, g_yyy_yyyz_1, g_yyy_yyz_1, g_yyy_yyzz_1, g_yyy_yzz_1, g_yyy_yzzz_1, g_yyy_zzz_1, g_yyy_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyy_xxxx_0[i] = 4.0 * g_yyy_xxx_1[i] * fe_0 + g_yyy_xxxx_1[i] * pa_x[i];

        g_xyyy_xxxy_0[i] = 3.0 * g_yyy_xxy_1[i] * fe_0 + g_yyy_xxxy_1[i] * pa_x[i];

        g_xyyy_xxxz_0[i] = 3.0 * g_yyy_xxz_1[i] * fe_0 + g_yyy_xxxz_1[i] * pa_x[i];

        g_xyyy_xxyy_0[i] = 2.0 * g_yyy_xyy_1[i] * fe_0 + g_yyy_xxyy_1[i] * pa_x[i];

        g_xyyy_xxyz_0[i] = 2.0 * g_yyy_xyz_1[i] * fe_0 + g_yyy_xxyz_1[i] * pa_x[i];

        g_xyyy_xxzz_0[i] = 2.0 * g_yyy_xzz_1[i] * fe_0 + g_yyy_xxzz_1[i] * pa_x[i];

        g_xyyy_xyyy_0[i] = g_yyy_yyy_1[i] * fe_0 + g_yyy_xyyy_1[i] * pa_x[i];

        g_xyyy_xyyz_0[i] = g_yyy_yyz_1[i] * fe_0 + g_yyy_xyyz_1[i] * pa_x[i];

        g_xyyy_xyzz_0[i] = g_yyy_yzz_1[i] * fe_0 + g_yyy_xyzz_1[i] * pa_x[i];

        g_xyyy_xzzz_0[i] = g_yyy_zzz_1[i] * fe_0 + g_yyy_xzzz_1[i] * pa_x[i];

        g_xyyy_yyyy_0[i] = g_yyy_yyyy_1[i] * pa_x[i];

        g_xyyy_yyyz_0[i] = g_yyy_yyyz_1[i] * pa_x[i];

        g_xyyy_yyzz_0[i] = g_yyy_yyzz_1[i] * pa_x[i];

        g_xyyy_yzzz_0[i] = g_yyy_yzzz_1[i] * pa_x[i];

        g_xyyy_zzzz_0[i] = g_yyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : GG

    auto g_xyyz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 105);

    auto g_xyyz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 106);

    auto g_xyyz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 107);

    auto g_xyyz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 108);

    auto g_xyyz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 109);

    auto g_xyyz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 110);

    auto g_xyyz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 111);

    auto g_xyyz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 112);

    auto g_xyyz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 113);

    auto g_xyyz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 114);

    auto g_xyyz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 115);

    auto g_xyyz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 116);

    auto g_xyyz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 117);

    auto g_xyyz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 118);

    auto g_xyyz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 119);

    #pragma omp simd aligned(g_xyy_xxxx_1, g_xyy_xxxy_1, g_xyy_xxyy_1, g_xyy_xyyy_1, g_xyyz_xxxx_0, g_xyyz_xxxy_0, g_xyyz_xxxz_0, g_xyyz_xxyy_0, g_xyyz_xxyz_0, g_xyyz_xxzz_0, g_xyyz_xyyy_0, g_xyyz_xyyz_0, g_xyyz_xyzz_0, g_xyyz_xzzz_0, g_xyyz_yyyy_0, g_xyyz_yyyz_0, g_xyyz_yyzz_0, g_xyyz_yzzz_0, g_xyyz_zzzz_0, g_yyz_xxxz_1, g_yyz_xxyz_1, g_yyz_xxz_1, g_yyz_xxzz_1, g_yyz_xyyz_1, g_yyz_xyz_1, g_yyz_xyzz_1, g_yyz_xzz_1, g_yyz_xzzz_1, g_yyz_yyyy_1, g_yyz_yyyz_1, g_yyz_yyz_1, g_yyz_yyzz_1, g_yyz_yzz_1, g_yyz_yzzz_1, g_yyz_zzz_1, g_yyz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyz_xxxx_0[i] = g_xyy_xxxx_1[i] * pa_z[i];

        g_xyyz_xxxy_0[i] = g_xyy_xxxy_1[i] * pa_z[i];

        g_xyyz_xxxz_0[i] = 3.0 * g_yyz_xxz_1[i] * fe_0 + g_yyz_xxxz_1[i] * pa_x[i];

        g_xyyz_xxyy_0[i] = g_xyy_xxyy_1[i] * pa_z[i];

        g_xyyz_xxyz_0[i] = 2.0 * g_yyz_xyz_1[i] * fe_0 + g_yyz_xxyz_1[i] * pa_x[i];

        g_xyyz_xxzz_0[i] = 2.0 * g_yyz_xzz_1[i] * fe_0 + g_yyz_xxzz_1[i] * pa_x[i];

        g_xyyz_xyyy_0[i] = g_xyy_xyyy_1[i] * pa_z[i];

        g_xyyz_xyyz_0[i] = g_yyz_yyz_1[i] * fe_0 + g_yyz_xyyz_1[i] * pa_x[i];

        g_xyyz_xyzz_0[i] = g_yyz_yzz_1[i] * fe_0 + g_yyz_xyzz_1[i] * pa_x[i];

        g_xyyz_xzzz_0[i] = g_yyz_zzz_1[i] * fe_0 + g_yyz_xzzz_1[i] * pa_x[i];

        g_xyyz_yyyy_0[i] = g_yyz_yyyy_1[i] * pa_x[i];

        g_xyyz_yyyz_0[i] = g_yyz_yyyz_1[i] * pa_x[i];

        g_xyyz_yyzz_0[i] = g_yyz_yyzz_1[i] * pa_x[i];

        g_xyyz_yzzz_0[i] = g_yyz_yzzz_1[i] * pa_x[i];

        g_xyyz_zzzz_0[i] = g_yyz_zzzz_1[i] * pa_x[i];
    }

    // Set up 120-135 components of targeted buffer : GG

    auto g_xyzz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 120);

    auto g_xyzz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 121);

    auto g_xyzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 122);

    auto g_xyzz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 123);

    auto g_xyzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 124);

    auto g_xyzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 125);

    auto g_xyzz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 126);

    auto g_xyzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 127);

    auto g_xyzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 128);

    auto g_xyzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 129);

    auto g_xyzz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 130);

    auto g_xyzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 131);

    auto g_xyzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 132);

    auto g_xyzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 133);

    auto g_xyzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 134);

    #pragma omp simd aligned(g_xyzz_xxxx_0, g_xyzz_xxxy_0, g_xyzz_xxxz_0, g_xyzz_xxyy_0, g_xyzz_xxyz_0, g_xyzz_xxzz_0, g_xyzz_xyyy_0, g_xyzz_xyyz_0, g_xyzz_xyzz_0, g_xyzz_xzzz_0, g_xyzz_yyyy_0, g_xyzz_yyyz_0, g_xyzz_yyzz_0, g_xyzz_yzzz_0, g_xyzz_zzzz_0, g_xzz_xxxx_1, g_xzz_xxxz_1, g_xzz_xxzz_1, g_xzz_xzzz_1, g_yzz_xxxy_1, g_yzz_xxy_1, g_yzz_xxyy_1, g_yzz_xxyz_1, g_yzz_xyy_1, g_yzz_xyyy_1, g_yzz_xyyz_1, g_yzz_xyz_1, g_yzz_xyzz_1, g_yzz_yyy_1, g_yzz_yyyy_1, g_yzz_yyyz_1, g_yzz_yyz_1, g_yzz_yyzz_1, g_yzz_yzz_1, g_yzz_yzzz_1, g_yzz_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzz_xxxx_0[i] = g_xzz_xxxx_1[i] * pa_y[i];

        g_xyzz_xxxy_0[i] = 3.0 * g_yzz_xxy_1[i] * fe_0 + g_yzz_xxxy_1[i] * pa_x[i];

        g_xyzz_xxxz_0[i] = g_xzz_xxxz_1[i] * pa_y[i];

        g_xyzz_xxyy_0[i] = 2.0 * g_yzz_xyy_1[i] * fe_0 + g_yzz_xxyy_1[i] * pa_x[i];

        g_xyzz_xxyz_0[i] = 2.0 * g_yzz_xyz_1[i] * fe_0 + g_yzz_xxyz_1[i] * pa_x[i];

        g_xyzz_xxzz_0[i] = g_xzz_xxzz_1[i] * pa_y[i];

        g_xyzz_xyyy_0[i] = g_yzz_yyy_1[i] * fe_0 + g_yzz_xyyy_1[i] * pa_x[i];

        g_xyzz_xyyz_0[i] = g_yzz_yyz_1[i] * fe_0 + g_yzz_xyyz_1[i] * pa_x[i];

        g_xyzz_xyzz_0[i] = g_yzz_yzz_1[i] * fe_0 + g_yzz_xyzz_1[i] * pa_x[i];

        g_xyzz_xzzz_0[i] = g_xzz_xzzz_1[i] * pa_y[i];

        g_xyzz_yyyy_0[i] = g_yzz_yyyy_1[i] * pa_x[i];

        g_xyzz_yyyz_0[i] = g_yzz_yyyz_1[i] * pa_x[i];

        g_xyzz_yyzz_0[i] = g_yzz_yyzz_1[i] * pa_x[i];

        g_xyzz_yzzz_0[i] = g_yzz_yzzz_1[i] * pa_x[i];

        g_xyzz_zzzz_0[i] = g_yzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 135-150 components of targeted buffer : GG

    auto g_xzzz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 135);

    auto g_xzzz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 136);

    auto g_xzzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 137);

    auto g_xzzz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 138);

    auto g_xzzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 139);

    auto g_xzzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 140);

    auto g_xzzz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 141);

    auto g_xzzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 142);

    auto g_xzzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 143);

    auto g_xzzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 144);

    auto g_xzzz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 145);

    auto g_xzzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 146);

    auto g_xzzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 147);

    auto g_xzzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 148);

    auto g_xzzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 149);

    #pragma omp simd aligned(g_xzzz_xxxx_0, g_xzzz_xxxy_0, g_xzzz_xxxz_0, g_xzzz_xxyy_0, g_xzzz_xxyz_0, g_xzzz_xxzz_0, g_xzzz_xyyy_0, g_xzzz_xyyz_0, g_xzzz_xyzz_0, g_xzzz_xzzz_0, g_xzzz_yyyy_0, g_xzzz_yyyz_0, g_xzzz_yyzz_0, g_xzzz_yzzz_0, g_xzzz_zzzz_0, g_zzz_xxx_1, g_zzz_xxxx_1, g_zzz_xxxy_1, g_zzz_xxxz_1, g_zzz_xxy_1, g_zzz_xxyy_1, g_zzz_xxyz_1, g_zzz_xxz_1, g_zzz_xxzz_1, g_zzz_xyy_1, g_zzz_xyyy_1, g_zzz_xyyz_1, g_zzz_xyz_1, g_zzz_xyzz_1, g_zzz_xzz_1, g_zzz_xzzz_1, g_zzz_yyy_1, g_zzz_yyyy_1, g_zzz_yyyz_1, g_zzz_yyz_1, g_zzz_yyzz_1, g_zzz_yzz_1, g_zzz_yzzz_1, g_zzz_zzz_1, g_zzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzz_xxxx_0[i] = 4.0 * g_zzz_xxx_1[i] * fe_0 + g_zzz_xxxx_1[i] * pa_x[i];

        g_xzzz_xxxy_0[i] = 3.0 * g_zzz_xxy_1[i] * fe_0 + g_zzz_xxxy_1[i] * pa_x[i];

        g_xzzz_xxxz_0[i] = 3.0 * g_zzz_xxz_1[i] * fe_0 + g_zzz_xxxz_1[i] * pa_x[i];

        g_xzzz_xxyy_0[i] = 2.0 * g_zzz_xyy_1[i] * fe_0 + g_zzz_xxyy_1[i] * pa_x[i];

        g_xzzz_xxyz_0[i] = 2.0 * g_zzz_xyz_1[i] * fe_0 + g_zzz_xxyz_1[i] * pa_x[i];

        g_xzzz_xxzz_0[i] = 2.0 * g_zzz_xzz_1[i] * fe_0 + g_zzz_xxzz_1[i] * pa_x[i];

        g_xzzz_xyyy_0[i] = g_zzz_yyy_1[i] * fe_0 + g_zzz_xyyy_1[i] * pa_x[i];

        g_xzzz_xyyz_0[i] = g_zzz_yyz_1[i] * fe_0 + g_zzz_xyyz_1[i] * pa_x[i];

        g_xzzz_xyzz_0[i] = g_zzz_yzz_1[i] * fe_0 + g_zzz_xyzz_1[i] * pa_x[i];

        g_xzzz_xzzz_0[i] = g_zzz_zzz_1[i] * fe_0 + g_zzz_xzzz_1[i] * pa_x[i];

        g_xzzz_yyyy_0[i] = g_zzz_yyyy_1[i] * pa_x[i];

        g_xzzz_yyyz_0[i] = g_zzz_yyyz_1[i] * pa_x[i];

        g_xzzz_yyzz_0[i] = g_zzz_yyzz_1[i] * pa_x[i];

        g_xzzz_yzzz_0[i] = g_zzz_yzzz_1[i] * pa_x[i];

        g_xzzz_zzzz_0[i] = g_zzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : GG

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

    #pragma omp simd aligned(g_yy_xxxx_0, g_yy_xxxx_1, g_yy_xxxy_0, g_yy_xxxy_1, g_yy_xxxz_0, g_yy_xxxz_1, g_yy_xxyy_0, g_yy_xxyy_1, g_yy_xxyz_0, g_yy_xxyz_1, g_yy_xxzz_0, g_yy_xxzz_1, g_yy_xyyy_0, g_yy_xyyy_1, g_yy_xyyz_0, g_yy_xyyz_1, g_yy_xyzz_0, g_yy_xyzz_1, g_yy_xzzz_0, g_yy_xzzz_1, g_yy_yyyy_0, g_yy_yyyy_1, g_yy_yyyz_0, g_yy_yyyz_1, g_yy_yyzz_0, g_yy_yyzz_1, g_yy_yzzz_0, g_yy_yzzz_1, g_yy_zzzz_0, g_yy_zzzz_1, g_yyy_xxx_1, g_yyy_xxxx_1, g_yyy_xxxy_1, g_yyy_xxxz_1, g_yyy_xxy_1, g_yyy_xxyy_1, g_yyy_xxyz_1, g_yyy_xxz_1, g_yyy_xxzz_1, g_yyy_xyy_1, g_yyy_xyyy_1, g_yyy_xyyz_1, g_yyy_xyz_1, g_yyy_xyzz_1, g_yyy_xzz_1, g_yyy_xzzz_1, g_yyy_yyy_1, g_yyy_yyyy_1, g_yyy_yyyz_1, g_yyy_yyz_1, g_yyy_yyzz_1, g_yyy_yzz_1, g_yyy_yzzz_1, g_yyy_zzz_1, g_yyy_zzzz_1, g_yyyy_xxxx_0, g_yyyy_xxxy_0, g_yyyy_xxxz_0, g_yyyy_xxyy_0, g_yyyy_xxyz_0, g_yyyy_xxzz_0, g_yyyy_xyyy_0, g_yyyy_xyyz_0, g_yyyy_xyzz_0, g_yyyy_xzzz_0, g_yyyy_yyyy_0, g_yyyy_yyyz_0, g_yyyy_yyzz_0, g_yyyy_yzzz_0, g_yyyy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyy_xxxx_0[i] = 3.0 * g_yy_xxxx_0[i] * fbe_0 - 3.0 * g_yy_xxxx_1[i] * fz_be_0 + g_yyy_xxxx_1[i] * pa_y[i];

        g_yyyy_xxxy_0[i] = 3.0 * g_yy_xxxy_0[i] * fbe_0 - 3.0 * g_yy_xxxy_1[i] * fz_be_0 + g_yyy_xxx_1[i] * fe_0 + g_yyy_xxxy_1[i] * pa_y[i];

        g_yyyy_xxxz_0[i] = 3.0 * g_yy_xxxz_0[i] * fbe_0 - 3.0 * g_yy_xxxz_1[i] * fz_be_0 + g_yyy_xxxz_1[i] * pa_y[i];

        g_yyyy_xxyy_0[i] = 3.0 * g_yy_xxyy_0[i] * fbe_0 - 3.0 * g_yy_xxyy_1[i] * fz_be_0 + 2.0 * g_yyy_xxy_1[i] * fe_0 + g_yyy_xxyy_1[i] * pa_y[i];

        g_yyyy_xxyz_0[i] = 3.0 * g_yy_xxyz_0[i] * fbe_0 - 3.0 * g_yy_xxyz_1[i] * fz_be_0 + g_yyy_xxz_1[i] * fe_0 + g_yyy_xxyz_1[i] * pa_y[i];

        g_yyyy_xxzz_0[i] = 3.0 * g_yy_xxzz_0[i] * fbe_0 - 3.0 * g_yy_xxzz_1[i] * fz_be_0 + g_yyy_xxzz_1[i] * pa_y[i];

        g_yyyy_xyyy_0[i] = 3.0 * g_yy_xyyy_0[i] * fbe_0 - 3.0 * g_yy_xyyy_1[i] * fz_be_0 + 3.0 * g_yyy_xyy_1[i] * fe_0 + g_yyy_xyyy_1[i] * pa_y[i];

        g_yyyy_xyyz_0[i] = 3.0 * g_yy_xyyz_0[i] * fbe_0 - 3.0 * g_yy_xyyz_1[i] * fz_be_0 + 2.0 * g_yyy_xyz_1[i] * fe_0 + g_yyy_xyyz_1[i] * pa_y[i];

        g_yyyy_xyzz_0[i] = 3.0 * g_yy_xyzz_0[i] * fbe_0 - 3.0 * g_yy_xyzz_1[i] * fz_be_0 + g_yyy_xzz_1[i] * fe_0 + g_yyy_xyzz_1[i] * pa_y[i];

        g_yyyy_xzzz_0[i] = 3.0 * g_yy_xzzz_0[i] * fbe_0 - 3.0 * g_yy_xzzz_1[i] * fz_be_0 + g_yyy_xzzz_1[i] * pa_y[i];

        g_yyyy_yyyy_0[i] = 3.0 * g_yy_yyyy_0[i] * fbe_0 - 3.0 * g_yy_yyyy_1[i] * fz_be_0 + 4.0 * g_yyy_yyy_1[i] * fe_0 + g_yyy_yyyy_1[i] * pa_y[i];

        g_yyyy_yyyz_0[i] = 3.0 * g_yy_yyyz_0[i] * fbe_0 - 3.0 * g_yy_yyyz_1[i] * fz_be_0 + 3.0 * g_yyy_yyz_1[i] * fe_0 + g_yyy_yyyz_1[i] * pa_y[i];

        g_yyyy_yyzz_0[i] = 3.0 * g_yy_yyzz_0[i] * fbe_0 - 3.0 * g_yy_yyzz_1[i] * fz_be_0 + 2.0 * g_yyy_yzz_1[i] * fe_0 + g_yyy_yyzz_1[i] * pa_y[i];

        g_yyyy_yzzz_0[i] = 3.0 * g_yy_yzzz_0[i] * fbe_0 - 3.0 * g_yy_yzzz_1[i] * fz_be_0 + g_yyy_zzz_1[i] * fe_0 + g_yyy_yzzz_1[i] * pa_y[i];

        g_yyyy_zzzz_0[i] = 3.0 * g_yy_zzzz_0[i] * fbe_0 - 3.0 * g_yy_zzzz_1[i] * fz_be_0 + g_yyy_zzzz_1[i] * pa_y[i];
    }

    // Set up 165-180 components of targeted buffer : GG

    auto g_yyyz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 165);

    auto g_yyyz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 166);

    auto g_yyyz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 167);

    auto g_yyyz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 168);

    auto g_yyyz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 169);

    auto g_yyyz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 170);

    auto g_yyyz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 171);

    auto g_yyyz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 172);

    auto g_yyyz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 173);

    auto g_yyyz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 174);

    auto g_yyyz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 175);

    auto g_yyyz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 176);

    auto g_yyyz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 177);

    auto g_yyyz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 178);

    auto g_yyyz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 179);

    #pragma omp simd aligned(g_yyy_xxx_1, g_yyy_xxxx_1, g_yyy_xxxy_1, g_yyy_xxxz_1, g_yyy_xxy_1, g_yyy_xxyy_1, g_yyy_xxyz_1, g_yyy_xxz_1, g_yyy_xxzz_1, g_yyy_xyy_1, g_yyy_xyyy_1, g_yyy_xyyz_1, g_yyy_xyz_1, g_yyy_xyzz_1, g_yyy_xzz_1, g_yyy_xzzz_1, g_yyy_yyy_1, g_yyy_yyyy_1, g_yyy_yyyz_1, g_yyy_yyz_1, g_yyy_yyzz_1, g_yyy_yzz_1, g_yyy_yzzz_1, g_yyy_zzz_1, g_yyy_zzzz_1, g_yyyz_xxxx_0, g_yyyz_xxxy_0, g_yyyz_xxxz_0, g_yyyz_xxyy_0, g_yyyz_xxyz_0, g_yyyz_xxzz_0, g_yyyz_xyyy_0, g_yyyz_xyyz_0, g_yyyz_xyzz_0, g_yyyz_xzzz_0, g_yyyz_yyyy_0, g_yyyz_yyyz_0, g_yyyz_yyzz_0, g_yyyz_yzzz_0, g_yyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyz_xxxx_0[i] = g_yyy_xxxx_1[i] * pa_z[i];

        g_yyyz_xxxy_0[i] = g_yyy_xxxy_1[i] * pa_z[i];

        g_yyyz_xxxz_0[i] = g_yyy_xxx_1[i] * fe_0 + g_yyy_xxxz_1[i] * pa_z[i];

        g_yyyz_xxyy_0[i] = g_yyy_xxyy_1[i] * pa_z[i];

        g_yyyz_xxyz_0[i] = g_yyy_xxy_1[i] * fe_0 + g_yyy_xxyz_1[i] * pa_z[i];

        g_yyyz_xxzz_0[i] = 2.0 * g_yyy_xxz_1[i] * fe_0 + g_yyy_xxzz_1[i] * pa_z[i];

        g_yyyz_xyyy_0[i] = g_yyy_xyyy_1[i] * pa_z[i];

        g_yyyz_xyyz_0[i] = g_yyy_xyy_1[i] * fe_0 + g_yyy_xyyz_1[i] * pa_z[i];

        g_yyyz_xyzz_0[i] = 2.0 * g_yyy_xyz_1[i] * fe_0 + g_yyy_xyzz_1[i] * pa_z[i];

        g_yyyz_xzzz_0[i] = 3.0 * g_yyy_xzz_1[i] * fe_0 + g_yyy_xzzz_1[i] * pa_z[i];

        g_yyyz_yyyy_0[i] = g_yyy_yyyy_1[i] * pa_z[i];

        g_yyyz_yyyz_0[i] = g_yyy_yyy_1[i] * fe_0 + g_yyy_yyyz_1[i] * pa_z[i];

        g_yyyz_yyzz_0[i] = 2.0 * g_yyy_yyz_1[i] * fe_0 + g_yyy_yyzz_1[i] * pa_z[i];

        g_yyyz_yzzz_0[i] = 3.0 * g_yyy_yzz_1[i] * fe_0 + g_yyy_yzzz_1[i] * pa_z[i];

        g_yyyz_zzzz_0[i] = 4.0 * g_yyy_zzz_1[i] * fe_0 + g_yyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 180-195 components of targeted buffer : GG

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

    #pragma omp simd aligned(g_yy_xxxy_0, g_yy_xxxy_1, g_yy_xxyy_0, g_yy_xxyy_1, g_yy_xyyy_0, g_yy_xyyy_1, g_yy_yyyy_0, g_yy_yyyy_1, g_yyz_xxxy_1, g_yyz_xxyy_1, g_yyz_xyyy_1, g_yyz_yyyy_1, g_yyzz_xxxx_0, g_yyzz_xxxy_0, g_yyzz_xxxz_0, g_yyzz_xxyy_0, g_yyzz_xxyz_0, g_yyzz_xxzz_0, g_yyzz_xyyy_0, g_yyzz_xyyz_0, g_yyzz_xyzz_0, g_yyzz_xzzz_0, g_yyzz_yyyy_0, g_yyzz_yyyz_0, g_yyzz_yyzz_0, g_yyzz_yzzz_0, g_yyzz_zzzz_0, g_yzz_xxxx_1, g_yzz_xxxz_1, g_yzz_xxyz_1, g_yzz_xxz_1, g_yzz_xxzz_1, g_yzz_xyyz_1, g_yzz_xyz_1, g_yzz_xyzz_1, g_yzz_xzz_1, g_yzz_xzzz_1, g_yzz_yyyz_1, g_yzz_yyz_1, g_yzz_yyzz_1, g_yzz_yzz_1, g_yzz_yzzz_1, g_yzz_zzz_1, g_yzz_zzzz_1, g_zz_xxxx_0, g_zz_xxxx_1, g_zz_xxxz_0, g_zz_xxxz_1, g_zz_xxyz_0, g_zz_xxyz_1, g_zz_xxzz_0, g_zz_xxzz_1, g_zz_xyyz_0, g_zz_xyyz_1, g_zz_xyzz_0, g_zz_xyzz_1, g_zz_xzzz_0, g_zz_xzzz_1, g_zz_yyyz_0, g_zz_yyyz_1, g_zz_yyzz_0, g_zz_yyzz_1, g_zz_yzzz_0, g_zz_yzzz_1, g_zz_zzzz_0, g_zz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzz_xxxx_0[i] = g_zz_xxxx_0[i] * fbe_0 - g_zz_xxxx_1[i] * fz_be_0 + g_yzz_xxxx_1[i] * pa_y[i];

        g_yyzz_xxxy_0[i] = g_yy_xxxy_0[i] * fbe_0 - g_yy_xxxy_1[i] * fz_be_0 + g_yyz_xxxy_1[i] * pa_z[i];

        g_yyzz_xxxz_0[i] = g_zz_xxxz_0[i] * fbe_0 - g_zz_xxxz_1[i] * fz_be_0 + g_yzz_xxxz_1[i] * pa_y[i];

        g_yyzz_xxyy_0[i] = g_yy_xxyy_0[i] * fbe_0 - g_yy_xxyy_1[i] * fz_be_0 + g_yyz_xxyy_1[i] * pa_z[i];

        g_yyzz_xxyz_0[i] = g_zz_xxyz_0[i] * fbe_0 - g_zz_xxyz_1[i] * fz_be_0 + g_yzz_xxz_1[i] * fe_0 + g_yzz_xxyz_1[i] * pa_y[i];

        g_yyzz_xxzz_0[i] = g_zz_xxzz_0[i] * fbe_0 - g_zz_xxzz_1[i] * fz_be_0 + g_yzz_xxzz_1[i] * pa_y[i];

        g_yyzz_xyyy_0[i] = g_yy_xyyy_0[i] * fbe_0 - g_yy_xyyy_1[i] * fz_be_0 + g_yyz_xyyy_1[i] * pa_z[i];

        g_yyzz_xyyz_0[i] = g_zz_xyyz_0[i] * fbe_0 - g_zz_xyyz_1[i] * fz_be_0 + 2.0 * g_yzz_xyz_1[i] * fe_0 + g_yzz_xyyz_1[i] * pa_y[i];

        g_yyzz_xyzz_0[i] = g_zz_xyzz_0[i] * fbe_0 - g_zz_xyzz_1[i] * fz_be_0 + g_yzz_xzz_1[i] * fe_0 + g_yzz_xyzz_1[i] * pa_y[i];

        g_yyzz_xzzz_0[i] = g_zz_xzzz_0[i] * fbe_0 - g_zz_xzzz_1[i] * fz_be_0 + g_yzz_xzzz_1[i] * pa_y[i];

        g_yyzz_yyyy_0[i] = g_yy_yyyy_0[i] * fbe_0 - g_yy_yyyy_1[i] * fz_be_0 + g_yyz_yyyy_1[i] * pa_z[i];

        g_yyzz_yyyz_0[i] = g_zz_yyyz_0[i] * fbe_0 - g_zz_yyyz_1[i] * fz_be_0 + 3.0 * g_yzz_yyz_1[i] * fe_0 + g_yzz_yyyz_1[i] * pa_y[i];

        g_yyzz_yyzz_0[i] = g_zz_yyzz_0[i] * fbe_0 - g_zz_yyzz_1[i] * fz_be_0 + 2.0 * g_yzz_yzz_1[i] * fe_0 + g_yzz_yyzz_1[i] * pa_y[i];

        g_yyzz_yzzz_0[i] = g_zz_yzzz_0[i] * fbe_0 - g_zz_yzzz_1[i] * fz_be_0 + g_yzz_zzz_1[i] * fe_0 + g_yzz_yzzz_1[i] * pa_y[i];

        g_yyzz_zzzz_0[i] = g_zz_zzzz_0[i] * fbe_0 - g_zz_zzzz_1[i] * fz_be_0 + g_yzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 195-210 components of targeted buffer : GG

    auto g_yzzz_xxxx_0 = pbuffer.data(idx_eri_0_gg + 195);

    auto g_yzzz_xxxy_0 = pbuffer.data(idx_eri_0_gg + 196);

    auto g_yzzz_xxxz_0 = pbuffer.data(idx_eri_0_gg + 197);

    auto g_yzzz_xxyy_0 = pbuffer.data(idx_eri_0_gg + 198);

    auto g_yzzz_xxyz_0 = pbuffer.data(idx_eri_0_gg + 199);

    auto g_yzzz_xxzz_0 = pbuffer.data(idx_eri_0_gg + 200);

    auto g_yzzz_xyyy_0 = pbuffer.data(idx_eri_0_gg + 201);

    auto g_yzzz_xyyz_0 = pbuffer.data(idx_eri_0_gg + 202);

    auto g_yzzz_xyzz_0 = pbuffer.data(idx_eri_0_gg + 203);

    auto g_yzzz_xzzz_0 = pbuffer.data(idx_eri_0_gg + 204);

    auto g_yzzz_yyyy_0 = pbuffer.data(idx_eri_0_gg + 205);

    auto g_yzzz_yyyz_0 = pbuffer.data(idx_eri_0_gg + 206);

    auto g_yzzz_yyzz_0 = pbuffer.data(idx_eri_0_gg + 207);

    auto g_yzzz_yzzz_0 = pbuffer.data(idx_eri_0_gg + 208);

    auto g_yzzz_zzzz_0 = pbuffer.data(idx_eri_0_gg + 209);

    #pragma omp simd aligned(g_yzzz_xxxx_0, g_yzzz_xxxy_0, g_yzzz_xxxz_0, g_yzzz_xxyy_0, g_yzzz_xxyz_0, g_yzzz_xxzz_0, g_yzzz_xyyy_0, g_yzzz_xyyz_0, g_yzzz_xyzz_0, g_yzzz_xzzz_0, g_yzzz_yyyy_0, g_yzzz_yyyz_0, g_yzzz_yyzz_0, g_yzzz_yzzz_0, g_yzzz_zzzz_0, g_zzz_xxx_1, g_zzz_xxxx_1, g_zzz_xxxy_1, g_zzz_xxxz_1, g_zzz_xxy_1, g_zzz_xxyy_1, g_zzz_xxyz_1, g_zzz_xxz_1, g_zzz_xxzz_1, g_zzz_xyy_1, g_zzz_xyyy_1, g_zzz_xyyz_1, g_zzz_xyz_1, g_zzz_xyzz_1, g_zzz_xzz_1, g_zzz_xzzz_1, g_zzz_yyy_1, g_zzz_yyyy_1, g_zzz_yyyz_1, g_zzz_yyz_1, g_zzz_yyzz_1, g_zzz_yzz_1, g_zzz_yzzz_1, g_zzz_zzz_1, g_zzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzz_xxxx_0[i] = g_zzz_xxxx_1[i] * pa_y[i];

        g_yzzz_xxxy_0[i] = g_zzz_xxx_1[i] * fe_0 + g_zzz_xxxy_1[i] * pa_y[i];

        g_yzzz_xxxz_0[i] = g_zzz_xxxz_1[i] * pa_y[i];

        g_yzzz_xxyy_0[i] = 2.0 * g_zzz_xxy_1[i] * fe_0 + g_zzz_xxyy_1[i] * pa_y[i];

        g_yzzz_xxyz_0[i] = g_zzz_xxz_1[i] * fe_0 + g_zzz_xxyz_1[i] * pa_y[i];

        g_yzzz_xxzz_0[i] = g_zzz_xxzz_1[i] * pa_y[i];

        g_yzzz_xyyy_0[i] = 3.0 * g_zzz_xyy_1[i] * fe_0 + g_zzz_xyyy_1[i] * pa_y[i];

        g_yzzz_xyyz_0[i] = 2.0 * g_zzz_xyz_1[i] * fe_0 + g_zzz_xyyz_1[i] * pa_y[i];

        g_yzzz_xyzz_0[i] = g_zzz_xzz_1[i] * fe_0 + g_zzz_xyzz_1[i] * pa_y[i];

        g_yzzz_xzzz_0[i] = g_zzz_xzzz_1[i] * pa_y[i];

        g_yzzz_yyyy_0[i] = 4.0 * g_zzz_yyy_1[i] * fe_0 + g_zzz_yyyy_1[i] * pa_y[i];

        g_yzzz_yyyz_0[i] = 3.0 * g_zzz_yyz_1[i] * fe_0 + g_zzz_yyyz_1[i] * pa_y[i];

        g_yzzz_yyzz_0[i] = 2.0 * g_zzz_yzz_1[i] * fe_0 + g_zzz_yyzz_1[i] * pa_y[i];

        g_yzzz_yzzz_0[i] = g_zzz_zzz_1[i] * fe_0 + g_zzz_yzzz_1[i] * pa_y[i];

        g_yzzz_zzzz_0[i] = g_zzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 210-225 components of targeted buffer : GG

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

    #pragma omp simd aligned(g_zz_xxxx_0, g_zz_xxxx_1, g_zz_xxxy_0, g_zz_xxxy_1, g_zz_xxxz_0, g_zz_xxxz_1, g_zz_xxyy_0, g_zz_xxyy_1, g_zz_xxyz_0, g_zz_xxyz_1, g_zz_xxzz_0, g_zz_xxzz_1, g_zz_xyyy_0, g_zz_xyyy_1, g_zz_xyyz_0, g_zz_xyyz_1, g_zz_xyzz_0, g_zz_xyzz_1, g_zz_xzzz_0, g_zz_xzzz_1, g_zz_yyyy_0, g_zz_yyyy_1, g_zz_yyyz_0, g_zz_yyyz_1, g_zz_yyzz_0, g_zz_yyzz_1, g_zz_yzzz_0, g_zz_yzzz_1, g_zz_zzzz_0, g_zz_zzzz_1, g_zzz_xxx_1, g_zzz_xxxx_1, g_zzz_xxxy_1, g_zzz_xxxz_1, g_zzz_xxy_1, g_zzz_xxyy_1, g_zzz_xxyz_1, g_zzz_xxz_1, g_zzz_xxzz_1, g_zzz_xyy_1, g_zzz_xyyy_1, g_zzz_xyyz_1, g_zzz_xyz_1, g_zzz_xyzz_1, g_zzz_xzz_1, g_zzz_xzzz_1, g_zzz_yyy_1, g_zzz_yyyy_1, g_zzz_yyyz_1, g_zzz_yyz_1, g_zzz_yyzz_1, g_zzz_yzz_1, g_zzz_yzzz_1, g_zzz_zzz_1, g_zzz_zzzz_1, g_zzzz_xxxx_0, g_zzzz_xxxy_0, g_zzzz_xxxz_0, g_zzzz_xxyy_0, g_zzzz_xxyz_0, g_zzzz_xxzz_0, g_zzzz_xyyy_0, g_zzzz_xyyz_0, g_zzzz_xyzz_0, g_zzzz_xzzz_0, g_zzzz_yyyy_0, g_zzzz_yyyz_0, g_zzzz_yyzz_0, g_zzzz_yzzz_0, g_zzzz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzz_xxxx_0[i] = 3.0 * g_zz_xxxx_0[i] * fbe_0 - 3.0 * g_zz_xxxx_1[i] * fz_be_0 + g_zzz_xxxx_1[i] * pa_z[i];

        g_zzzz_xxxy_0[i] = 3.0 * g_zz_xxxy_0[i] * fbe_0 - 3.0 * g_zz_xxxy_1[i] * fz_be_0 + g_zzz_xxxy_1[i] * pa_z[i];

        g_zzzz_xxxz_0[i] = 3.0 * g_zz_xxxz_0[i] * fbe_0 - 3.0 * g_zz_xxxz_1[i] * fz_be_0 + g_zzz_xxx_1[i] * fe_0 + g_zzz_xxxz_1[i] * pa_z[i];

        g_zzzz_xxyy_0[i] = 3.0 * g_zz_xxyy_0[i] * fbe_0 - 3.0 * g_zz_xxyy_1[i] * fz_be_0 + g_zzz_xxyy_1[i] * pa_z[i];

        g_zzzz_xxyz_0[i] = 3.0 * g_zz_xxyz_0[i] * fbe_0 - 3.0 * g_zz_xxyz_1[i] * fz_be_0 + g_zzz_xxy_1[i] * fe_0 + g_zzz_xxyz_1[i] * pa_z[i];

        g_zzzz_xxzz_0[i] = 3.0 * g_zz_xxzz_0[i] * fbe_0 - 3.0 * g_zz_xxzz_1[i] * fz_be_0 + 2.0 * g_zzz_xxz_1[i] * fe_0 + g_zzz_xxzz_1[i] * pa_z[i];

        g_zzzz_xyyy_0[i] = 3.0 * g_zz_xyyy_0[i] * fbe_0 - 3.0 * g_zz_xyyy_1[i] * fz_be_0 + g_zzz_xyyy_1[i] * pa_z[i];

        g_zzzz_xyyz_0[i] = 3.0 * g_zz_xyyz_0[i] * fbe_0 - 3.0 * g_zz_xyyz_1[i] * fz_be_0 + g_zzz_xyy_1[i] * fe_0 + g_zzz_xyyz_1[i] * pa_z[i];

        g_zzzz_xyzz_0[i] = 3.0 * g_zz_xyzz_0[i] * fbe_0 - 3.0 * g_zz_xyzz_1[i] * fz_be_0 + 2.0 * g_zzz_xyz_1[i] * fe_0 + g_zzz_xyzz_1[i] * pa_z[i];

        g_zzzz_xzzz_0[i] = 3.0 * g_zz_xzzz_0[i] * fbe_0 - 3.0 * g_zz_xzzz_1[i] * fz_be_0 + 3.0 * g_zzz_xzz_1[i] * fe_0 + g_zzz_xzzz_1[i] * pa_z[i];

        g_zzzz_yyyy_0[i] = 3.0 * g_zz_yyyy_0[i] * fbe_0 - 3.0 * g_zz_yyyy_1[i] * fz_be_0 + g_zzz_yyyy_1[i] * pa_z[i];

        g_zzzz_yyyz_0[i] = 3.0 * g_zz_yyyz_0[i] * fbe_0 - 3.0 * g_zz_yyyz_1[i] * fz_be_0 + g_zzz_yyy_1[i] * fe_0 + g_zzz_yyyz_1[i] * pa_z[i];

        g_zzzz_yyzz_0[i] = 3.0 * g_zz_yyzz_0[i] * fbe_0 - 3.0 * g_zz_yyzz_1[i] * fz_be_0 + 2.0 * g_zzz_yyz_1[i] * fe_0 + g_zzz_yyzz_1[i] * pa_z[i];

        g_zzzz_yzzz_0[i] = 3.0 * g_zz_yzzz_0[i] * fbe_0 - 3.0 * g_zz_yzzz_1[i] * fz_be_0 + 3.0 * g_zzz_yzz_1[i] * fe_0 + g_zzz_yzzz_1[i] * pa_z[i];

        g_zzzz_zzzz_0[i] = 3.0 * g_zz_zzzz_0[i] * fbe_0 - 3.0 * g_zz_zzzz_1[i] * fz_be_0 + 4.0 * g_zzz_zzz_1[i] * fe_0 + g_zzz_zzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

