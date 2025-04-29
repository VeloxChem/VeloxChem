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

#include "TwoCenterElectronRepulsionPrimRecHG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_hg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_hg,
                                const size_t idx_eri_0_fg,
                                const size_t idx_eri_1_fg,
                                const size_t idx_eri_1_gf,
                                const size_t idx_eri_1_gg,
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

    // Set up components of auxiliary buffer : FG

    auto g_xxx_xxxx_0 = pbuffer.data(idx_eri_0_fg);

    auto g_xxx_xxxy_0 = pbuffer.data(idx_eri_0_fg + 1);

    auto g_xxx_xxxz_0 = pbuffer.data(idx_eri_0_fg + 2);

    auto g_xxx_xxyy_0 = pbuffer.data(idx_eri_0_fg + 3);

    auto g_xxx_xxyz_0 = pbuffer.data(idx_eri_0_fg + 4);

    auto g_xxx_xxzz_0 = pbuffer.data(idx_eri_0_fg + 5);

    auto g_xxx_xyyy_0 = pbuffer.data(idx_eri_0_fg + 6);

    auto g_xxx_xyyz_0 = pbuffer.data(idx_eri_0_fg + 7);

    auto g_xxx_xyzz_0 = pbuffer.data(idx_eri_0_fg + 8);

    auto g_xxx_xzzz_0 = pbuffer.data(idx_eri_0_fg + 9);

    auto g_xxx_yyyy_0 = pbuffer.data(idx_eri_0_fg + 10);

    auto g_xxx_yyyz_0 = pbuffer.data(idx_eri_0_fg + 11);

    auto g_xxx_yyzz_0 = pbuffer.data(idx_eri_0_fg + 12);

    auto g_xxx_yzzz_0 = pbuffer.data(idx_eri_0_fg + 13);

    auto g_xxx_zzzz_0 = pbuffer.data(idx_eri_0_fg + 14);

    auto g_xxy_xxxx_0 = pbuffer.data(idx_eri_0_fg + 15);

    auto g_xxy_xxxz_0 = pbuffer.data(idx_eri_0_fg + 17);

    auto g_xxy_xxzz_0 = pbuffer.data(idx_eri_0_fg + 20);

    auto g_xxy_xzzz_0 = pbuffer.data(idx_eri_0_fg + 24);

    auto g_xxz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 30);

    auto g_xxz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 31);

    auto g_xxz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 33);

    auto g_xxz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 36);

    auto g_xyy_xxxy_0 = pbuffer.data(idx_eri_0_fg + 46);

    auto g_xyy_xxyy_0 = pbuffer.data(idx_eri_0_fg + 48);

    auto g_xyy_xxyz_0 = pbuffer.data(idx_eri_0_fg + 49);

    auto g_xyy_xyyy_0 = pbuffer.data(idx_eri_0_fg + 51);

    auto g_xyy_xyyz_0 = pbuffer.data(idx_eri_0_fg + 52);

    auto g_xyy_xyzz_0 = pbuffer.data(idx_eri_0_fg + 53);

    auto g_xyy_yyyy_0 = pbuffer.data(idx_eri_0_fg + 55);

    auto g_xyy_yyyz_0 = pbuffer.data(idx_eri_0_fg + 56);

    auto g_xyy_yyzz_0 = pbuffer.data(idx_eri_0_fg + 57);

    auto g_xyy_yzzz_0 = pbuffer.data(idx_eri_0_fg + 58);

    auto g_xyy_zzzz_0 = pbuffer.data(idx_eri_0_fg + 59);

    auto g_xzz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 77);

    auto g_xzz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 79);

    auto g_xzz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 80);

    auto g_xzz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 82);

    auto g_xzz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 83);

    auto g_xzz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 84);

    auto g_xzz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 85);

    auto g_xzz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 86);

    auto g_xzz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 87);

    auto g_xzz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 88);

    auto g_xzz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 89);

    auto g_yyy_xxxx_0 = pbuffer.data(idx_eri_0_fg + 90);

    auto g_yyy_xxxy_0 = pbuffer.data(idx_eri_0_fg + 91);

    auto g_yyy_xxxz_0 = pbuffer.data(idx_eri_0_fg + 92);

    auto g_yyy_xxyy_0 = pbuffer.data(idx_eri_0_fg + 93);

    auto g_yyy_xxyz_0 = pbuffer.data(idx_eri_0_fg + 94);

    auto g_yyy_xxzz_0 = pbuffer.data(idx_eri_0_fg + 95);

    auto g_yyy_xyyy_0 = pbuffer.data(idx_eri_0_fg + 96);

    auto g_yyy_xyyz_0 = pbuffer.data(idx_eri_0_fg + 97);

    auto g_yyy_xyzz_0 = pbuffer.data(idx_eri_0_fg + 98);

    auto g_yyy_xzzz_0 = pbuffer.data(idx_eri_0_fg + 99);

    auto g_yyy_yyyy_0 = pbuffer.data(idx_eri_0_fg + 100);

    auto g_yyy_yyyz_0 = pbuffer.data(idx_eri_0_fg + 101);

    auto g_yyy_yyzz_0 = pbuffer.data(idx_eri_0_fg + 102);

    auto g_yyy_yzzz_0 = pbuffer.data(idx_eri_0_fg + 103);

    auto g_yyy_zzzz_0 = pbuffer.data(idx_eri_0_fg + 104);

    auto g_yyz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 106);

    auto g_yyz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 108);

    auto g_yyz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 111);

    auto g_yyz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 115);

    auto g_yzz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 120);

    auto g_yzz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 122);

    auto g_yzz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 124);

    auto g_yzz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 125);

    auto g_yzz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 127);

    auto g_yzz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 128);

    auto g_yzz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 129);

    auto g_yzz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 131);

    auto g_yzz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 132);

    auto g_yzz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 133);

    auto g_yzz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 134);

    auto g_zzz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 135);

    auto g_zzz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 136);

    auto g_zzz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 137);

    auto g_zzz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 138);

    auto g_zzz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 139);

    auto g_zzz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 140);

    auto g_zzz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 141);

    auto g_zzz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 142);

    auto g_zzz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 143);

    auto g_zzz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 144);

    auto g_zzz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 145);

    auto g_zzz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 146);

    auto g_zzz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 147);

    auto g_zzz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 148);

    auto g_zzz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 149);

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

    auto g_xxy_xxxz_1 = pbuffer.data(idx_eri_1_fg + 17);

    auto g_xxy_xxzz_1 = pbuffer.data(idx_eri_1_fg + 20);

    auto g_xxy_xzzz_1 = pbuffer.data(idx_eri_1_fg + 24);

    auto g_xxz_xxxx_1 = pbuffer.data(idx_eri_1_fg + 30);

    auto g_xxz_xxxy_1 = pbuffer.data(idx_eri_1_fg + 31);

    auto g_xxz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 33);

    auto g_xxz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 36);

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

    auto g_yyz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 108);

    auto g_yyz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 111);

    auto g_yyz_yyyy_1 = pbuffer.data(idx_eri_1_fg + 115);

    auto g_yzz_xxxx_1 = pbuffer.data(idx_eri_1_fg + 120);

    auto g_yzz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 122);

    auto g_yzz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 124);

    auto g_yzz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 125);

    auto g_yzz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 127);

    auto g_yzz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 128);

    auto g_yzz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 129);

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

    // Set up components of auxiliary buffer : GF

    auto g_xxxx_xxx_1 = pbuffer.data(idx_eri_1_gf);

    auto g_xxxx_xxy_1 = pbuffer.data(idx_eri_1_gf + 1);

    auto g_xxxx_xxz_1 = pbuffer.data(idx_eri_1_gf + 2);

    auto g_xxxx_xyy_1 = pbuffer.data(idx_eri_1_gf + 3);

    auto g_xxxx_xyz_1 = pbuffer.data(idx_eri_1_gf + 4);

    auto g_xxxx_xzz_1 = pbuffer.data(idx_eri_1_gf + 5);

    auto g_xxxx_yyy_1 = pbuffer.data(idx_eri_1_gf + 6);

    auto g_xxxx_yyz_1 = pbuffer.data(idx_eri_1_gf + 7);

    auto g_xxxx_yzz_1 = pbuffer.data(idx_eri_1_gf + 8);

    auto g_xxxx_zzz_1 = pbuffer.data(idx_eri_1_gf + 9);

    auto g_xxxz_xxz_1 = pbuffer.data(idx_eri_1_gf + 22);

    auto g_xxxz_xyz_1 = pbuffer.data(idx_eri_1_gf + 24);

    auto g_xxxz_xzz_1 = pbuffer.data(idx_eri_1_gf + 25);

    auto g_xxxz_yyz_1 = pbuffer.data(idx_eri_1_gf + 27);

    auto g_xxxz_yzz_1 = pbuffer.data(idx_eri_1_gf + 28);

    auto g_xxxz_zzz_1 = pbuffer.data(idx_eri_1_gf + 29);

    auto g_xxyy_xxx_1 = pbuffer.data(idx_eri_1_gf + 30);

    auto g_xxyy_xxy_1 = pbuffer.data(idx_eri_1_gf + 31);

    auto g_xxyy_xxz_1 = pbuffer.data(idx_eri_1_gf + 32);

    auto g_xxyy_xyy_1 = pbuffer.data(idx_eri_1_gf + 33);

    auto g_xxyy_xyz_1 = pbuffer.data(idx_eri_1_gf + 34);

    auto g_xxyy_xzz_1 = pbuffer.data(idx_eri_1_gf + 35);

    auto g_xxyy_yyy_1 = pbuffer.data(idx_eri_1_gf + 36);

    auto g_xxyy_yyz_1 = pbuffer.data(idx_eri_1_gf + 37);

    auto g_xxyy_yzz_1 = pbuffer.data(idx_eri_1_gf + 38);

    auto g_xxyy_zzz_1 = pbuffer.data(idx_eri_1_gf + 39);

    auto g_xxzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 50);

    auto g_xxzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 51);

    auto g_xxzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 52);

    auto g_xxzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 53);

    auto g_xxzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 54);

    auto g_xxzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 55);

    auto g_xxzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 56);

    auto g_xxzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 57);

    auto g_xxzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 58);

    auto g_xxzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 59);

    auto g_xyyy_xxy_1 = pbuffer.data(idx_eri_1_gf + 61);

    auto g_xyyy_xyy_1 = pbuffer.data(idx_eri_1_gf + 63);

    auto g_xyyy_xyz_1 = pbuffer.data(idx_eri_1_gf + 64);

    auto g_xyyy_yyy_1 = pbuffer.data(idx_eri_1_gf + 66);

    auto g_xyyy_yyz_1 = pbuffer.data(idx_eri_1_gf + 67);

    auto g_xyyy_yzz_1 = pbuffer.data(idx_eri_1_gf + 68);

    auto g_xzzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 92);

    auto g_xzzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 94);

    auto g_xzzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 95);

    auto g_xzzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 97);

    auto g_xzzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 98);

    auto g_xzzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 99);

    auto g_yyyy_xxx_1 = pbuffer.data(idx_eri_1_gf + 100);

    auto g_yyyy_xxy_1 = pbuffer.data(idx_eri_1_gf + 101);

    auto g_yyyy_xxz_1 = pbuffer.data(idx_eri_1_gf + 102);

    auto g_yyyy_xyy_1 = pbuffer.data(idx_eri_1_gf + 103);

    auto g_yyyy_xyz_1 = pbuffer.data(idx_eri_1_gf + 104);

    auto g_yyyy_xzz_1 = pbuffer.data(idx_eri_1_gf + 105);

    auto g_yyyy_yyy_1 = pbuffer.data(idx_eri_1_gf + 106);

    auto g_yyyy_yyz_1 = pbuffer.data(idx_eri_1_gf + 107);

    auto g_yyyy_yzz_1 = pbuffer.data(idx_eri_1_gf + 108);

    auto g_yyyy_zzz_1 = pbuffer.data(idx_eri_1_gf + 109);

    auto g_yyyz_xxz_1 = pbuffer.data(idx_eri_1_gf + 112);

    auto g_yyyz_xyz_1 = pbuffer.data(idx_eri_1_gf + 114);

    auto g_yyyz_xzz_1 = pbuffer.data(idx_eri_1_gf + 115);

    auto g_yyyz_yyz_1 = pbuffer.data(idx_eri_1_gf + 117);

    auto g_yyyz_yzz_1 = pbuffer.data(idx_eri_1_gf + 118);

    auto g_yyyz_zzz_1 = pbuffer.data(idx_eri_1_gf + 119);

    auto g_yyzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 120);

    auto g_yyzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 121);

    auto g_yyzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 122);

    auto g_yyzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 123);

    auto g_yyzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 124);

    auto g_yyzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 125);

    auto g_yyzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 126);

    auto g_yyzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 127);

    auto g_yyzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 128);

    auto g_yyzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 129);

    auto g_yzzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 131);

    auto g_yzzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 132);

    auto g_yzzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 133);

    auto g_yzzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 134);

    auto g_yzzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 135);

    auto g_yzzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 136);

    auto g_yzzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 137);

    auto g_yzzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 138);

    auto g_yzzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 139);

    auto g_zzzz_xxx_1 = pbuffer.data(idx_eri_1_gf + 140);

    auto g_zzzz_xxy_1 = pbuffer.data(idx_eri_1_gf + 141);

    auto g_zzzz_xxz_1 = pbuffer.data(idx_eri_1_gf + 142);

    auto g_zzzz_xyy_1 = pbuffer.data(idx_eri_1_gf + 143);

    auto g_zzzz_xyz_1 = pbuffer.data(idx_eri_1_gf + 144);

    auto g_zzzz_xzz_1 = pbuffer.data(idx_eri_1_gf + 145);

    auto g_zzzz_yyy_1 = pbuffer.data(idx_eri_1_gf + 146);

    auto g_zzzz_yyz_1 = pbuffer.data(idx_eri_1_gf + 147);

    auto g_zzzz_yzz_1 = pbuffer.data(idx_eri_1_gf + 148);

    auto g_zzzz_zzz_1 = pbuffer.data(idx_eri_1_gf + 149);

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

    auto g_xxxy_xxxy_1 = pbuffer.data(idx_eri_1_gg + 16);

    auto g_xxxy_xxxz_1 = pbuffer.data(idx_eri_1_gg + 17);

    auto g_xxxy_xxyy_1 = pbuffer.data(idx_eri_1_gg + 18);

    auto g_xxxy_xxzz_1 = pbuffer.data(idx_eri_1_gg + 20);

    auto g_xxxy_xyyy_1 = pbuffer.data(idx_eri_1_gg + 21);

    auto g_xxxy_xzzz_1 = pbuffer.data(idx_eri_1_gg + 24);

    auto g_xxxy_yyyy_1 = pbuffer.data(idx_eri_1_gg + 25);

    auto g_xxxz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 30);

    auto g_xxxz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 31);

    auto g_xxxz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 32);

    auto g_xxxz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 33);

    auto g_xxxz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 34);

    auto g_xxxz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 35);

    auto g_xxxz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 36);

    auto g_xxxz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 37);

    auto g_xxxz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 38);

    auto g_xxxz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 39);

    auto g_xxxz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 41);

    auto g_xxxz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 42);

    auto g_xxxz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 43);

    auto g_xxxz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 44);

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

    auto g_xyyy_xxxx_1 = pbuffer.data(idx_eri_1_gg + 90);

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

    auto g_xzzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 135);

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

    auto g_yyyz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 167);

    auto g_yyyz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 168);

    auto g_yyyz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 169);

    auto g_yyyz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 170);

    auto g_yyyz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 171);

    auto g_yyyz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 172);

    auto g_yyyz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 173);

    auto g_yyyz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 174);

    auto g_yyyz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 175);

    auto g_yyyz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 176);

    auto g_yyyz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 177);

    auto g_yyyz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 178);

    auto g_yyyz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 179);

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

    auto g_yzzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 196);

    auto g_yzzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 197);

    auto g_yzzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 198);

    auto g_yzzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 199);

    auto g_yzzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 200);

    auto g_yzzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 201);

    auto g_yzzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 202);

    auto g_yzzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 203);

    auto g_yzzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 204);

    auto g_yzzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 205);

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

    // Set up 0-15 components of targeted buffer : HG

    auto g_xxxxx_xxxx_0 = pbuffer.data(idx_eri_0_hg);

    auto g_xxxxx_xxxy_0 = pbuffer.data(idx_eri_0_hg + 1);

    auto g_xxxxx_xxxz_0 = pbuffer.data(idx_eri_0_hg + 2);

    auto g_xxxxx_xxyy_0 = pbuffer.data(idx_eri_0_hg + 3);

    auto g_xxxxx_xxyz_0 = pbuffer.data(idx_eri_0_hg + 4);

    auto g_xxxxx_xxzz_0 = pbuffer.data(idx_eri_0_hg + 5);

    auto g_xxxxx_xyyy_0 = pbuffer.data(idx_eri_0_hg + 6);

    auto g_xxxxx_xyyz_0 = pbuffer.data(idx_eri_0_hg + 7);

    auto g_xxxxx_xyzz_0 = pbuffer.data(idx_eri_0_hg + 8);

    auto g_xxxxx_xzzz_0 = pbuffer.data(idx_eri_0_hg + 9);

    auto g_xxxxx_yyyy_0 = pbuffer.data(idx_eri_0_hg + 10);

    auto g_xxxxx_yyyz_0 = pbuffer.data(idx_eri_0_hg + 11);

    auto g_xxxxx_yyzz_0 = pbuffer.data(idx_eri_0_hg + 12);

    auto g_xxxxx_yzzz_0 = pbuffer.data(idx_eri_0_hg + 13);

    auto g_xxxxx_zzzz_0 = pbuffer.data(idx_eri_0_hg + 14);

    #pragma omp simd aligned(g_xxx_xxxx_0, g_xxx_xxxx_1, g_xxx_xxxy_0, g_xxx_xxxy_1, g_xxx_xxxz_0, g_xxx_xxxz_1, g_xxx_xxyy_0, g_xxx_xxyy_1, g_xxx_xxyz_0, g_xxx_xxyz_1, g_xxx_xxzz_0, g_xxx_xxzz_1, g_xxx_xyyy_0, g_xxx_xyyy_1, g_xxx_xyyz_0, g_xxx_xyyz_1, g_xxx_xyzz_0, g_xxx_xyzz_1, g_xxx_xzzz_0, g_xxx_xzzz_1, g_xxx_yyyy_0, g_xxx_yyyy_1, g_xxx_yyyz_0, g_xxx_yyyz_1, g_xxx_yyzz_0, g_xxx_yyzz_1, g_xxx_yzzz_0, g_xxx_yzzz_1, g_xxx_zzzz_0, g_xxx_zzzz_1, g_xxxx_xxx_1, g_xxxx_xxxx_1, g_xxxx_xxxy_1, g_xxxx_xxxz_1, g_xxxx_xxy_1, g_xxxx_xxyy_1, g_xxxx_xxyz_1, g_xxxx_xxz_1, g_xxxx_xxzz_1, g_xxxx_xyy_1, g_xxxx_xyyy_1, g_xxxx_xyyz_1, g_xxxx_xyz_1, g_xxxx_xyzz_1, g_xxxx_xzz_1, g_xxxx_xzzz_1, g_xxxx_yyy_1, g_xxxx_yyyy_1, g_xxxx_yyyz_1, g_xxxx_yyz_1, g_xxxx_yyzz_1, g_xxxx_yzz_1, g_xxxx_yzzz_1, g_xxxx_zzz_1, g_xxxx_zzzz_1, g_xxxxx_xxxx_0, g_xxxxx_xxxy_0, g_xxxxx_xxxz_0, g_xxxxx_xxyy_0, g_xxxxx_xxyz_0, g_xxxxx_xxzz_0, g_xxxxx_xyyy_0, g_xxxxx_xyyz_0, g_xxxxx_xyzz_0, g_xxxxx_xzzz_0, g_xxxxx_yyyy_0, g_xxxxx_yyyz_0, g_xxxxx_yyzz_0, g_xxxxx_yzzz_0, g_xxxxx_zzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxx_xxxx_0[i] = 4.0 * g_xxx_xxxx_0[i] * fbe_0 - 4.0 * g_xxx_xxxx_1[i] * fz_be_0 + 4.0 * g_xxxx_xxx_1[i] * fe_0 + g_xxxx_xxxx_1[i] * pa_x[i];

        g_xxxxx_xxxy_0[i] = 4.0 * g_xxx_xxxy_0[i] * fbe_0 - 4.0 * g_xxx_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxx_xxy_1[i] * fe_0 + g_xxxx_xxxy_1[i] * pa_x[i];

        g_xxxxx_xxxz_0[i] = 4.0 * g_xxx_xxxz_0[i] * fbe_0 - 4.0 * g_xxx_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxx_xxz_1[i] * fe_0 + g_xxxx_xxxz_1[i] * pa_x[i];

        g_xxxxx_xxyy_0[i] = 4.0 * g_xxx_xxyy_0[i] * fbe_0 - 4.0 * g_xxx_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxx_xyy_1[i] * fe_0 + g_xxxx_xxyy_1[i] * pa_x[i];

        g_xxxxx_xxyz_0[i] = 4.0 * g_xxx_xxyz_0[i] * fbe_0 - 4.0 * g_xxx_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxx_xyz_1[i] * fe_0 + g_xxxx_xxyz_1[i] * pa_x[i];

        g_xxxxx_xxzz_0[i] = 4.0 * g_xxx_xxzz_0[i] * fbe_0 - 4.0 * g_xxx_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxx_xzz_1[i] * fe_0 + g_xxxx_xxzz_1[i] * pa_x[i];

        g_xxxxx_xyyy_0[i] = 4.0 * g_xxx_xyyy_0[i] * fbe_0 - 4.0 * g_xxx_xyyy_1[i] * fz_be_0 + g_xxxx_yyy_1[i] * fe_0 + g_xxxx_xyyy_1[i] * pa_x[i];

        g_xxxxx_xyyz_0[i] = 4.0 * g_xxx_xyyz_0[i] * fbe_0 - 4.0 * g_xxx_xyyz_1[i] * fz_be_0 + g_xxxx_yyz_1[i] * fe_0 + g_xxxx_xyyz_1[i] * pa_x[i];

        g_xxxxx_xyzz_0[i] = 4.0 * g_xxx_xyzz_0[i] * fbe_0 - 4.0 * g_xxx_xyzz_1[i] * fz_be_0 + g_xxxx_yzz_1[i] * fe_0 + g_xxxx_xyzz_1[i] * pa_x[i];

        g_xxxxx_xzzz_0[i] = 4.0 * g_xxx_xzzz_0[i] * fbe_0 - 4.0 * g_xxx_xzzz_1[i] * fz_be_0 + g_xxxx_zzz_1[i] * fe_0 + g_xxxx_xzzz_1[i] * pa_x[i];

        g_xxxxx_yyyy_0[i] = 4.0 * g_xxx_yyyy_0[i] * fbe_0 - 4.0 * g_xxx_yyyy_1[i] * fz_be_0 + g_xxxx_yyyy_1[i] * pa_x[i];

        g_xxxxx_yyyz_0[i] = 4.0 * g_xxx_yyyz_0[i] * fbe_0 - 4.0 * g_xxx_yyyz_1[i] * fz_be_0 + g_xxxx_yyyz_1[i] * pa_x[i];

        g_xxxxx_yyzz_0[i] = 4.0 * g_xxx_yyzz_0[i] * fbe_0 - 4.0 * g_xxx_yyzz_1[i] * fz_be_0 + g_xxxx_yyzz_1[i] * pa_x[i];

        g_xxxxx_yzzz_0[i] = 4.0 * g_xxx_yzzz_0[i] * fbe_0 - 4.0 * g_xxx_yzzz_1[i] * fz_be_0 + g_xxxx_yzzz_1[i] * pa_x[i];

        g_xxxxx_zzzz_0[i] = 4.0 * g_xxx_zzzz_0[i] * fbe_0 - 4.0 * g_xxx_zzzz_1[i] * fz_be_0 + g_xxxx_zzzz_1[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : HG

    auto g_xxxxy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 15);

    auto g_xxxxy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 16);

    auto g_xxxxy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 17);

    auto g_xxxxy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 18);

    auto g_xxxxy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 19);

    auto g_xxxxy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 20);

    auto g_xxxxy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 21);

    auto g_xxxxy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 22);

    auto g_xxxxy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 23);

    auto g_xxxxy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 24);

    auto g_xxxxy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 25);

    auto g_xxxxy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 26);

    auto g_xxxxy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 27);

    auto g_xxxxy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 28);

    auto g_xxxxy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 29);

    #pragma omp simd aligned(g_xxxx_xxx_1, g_xxxx_xxxx_1, g_xxxx_xxxy_1, g_xxxx_xxxz_1, g_xxxx_xxy_1, g_xxxx_xxyy_1, g_xxxx_xxyz_1, g_xxxx_xxz_1, g_xxxx_xxzz_1, g_xxxx_xyy_1, g_xxxx_xyyy_1, g_xxxx_xyyz_1, g_xxxx_xyz_1, g_xxxx_xyzz_1, g_xxxx_xzz_1, g_xxxx_xzzz_1, g_xxxx_yyy_1, g_xxxx_yyyy_1, g_xxxx_yyyz_1, g_xxxx_yyz_1, g_xxxx_yyzz_1, g_xxxx_yzz_1, g_xxxx_yzzz_1, g_xxxx_zzz_1, g_xxxx_zzzz_1, g_xxxxy_xxxx_0, g_xxxxy_xxxy_0, g_xxxxy_xxxz_0, g_xxxxy_xxyy_0, g_xxxxy_xxyz_0, g_xxxxy_xxzz_0, g_xxxxy_xyyy_0, g_xxxxy_xyyz_0, g_xxxxy_xyzz_0, g_xxxxy_xzzz_0, g_xxxxy_yyyy_0, g_xxxxy_yyyz_0, g_xxxxy_yyzz_0, g_xxxxy_yzzz_0, g_xxxxy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxy_xxxx_0[i] = g_xxxx_xxxx_1[i] * pa_y[i];

        g_xxxxy_xxxy_0[i] = g_xxxx_xxx_1[i] * fe_0 + g_xxxx_xxxy_1[i] * pa_y[i];

        g_xxxxy_xxxz_0[i] = g_xxxx_xxxz_1[i] * pa_y[i];

        g_xxxxy_xxyy_0[i] = 2.0 * g_xxxx_xxy_1[i] * fe_0 + g_xxxx_xxyy_1[i] * pa_y[i];

        g_xxxxy_xxyz_0[i] = g_xxxx_xxz_1[i] * fe_0 + g_xxxx_xxyz_1[i] * pa_y[i];

        g_xxxxy_xxzz_0[i] = g_xxxx_xxzz_1[i] * pa_y[i];

        g_xxxxy_xyyy_0[i] = 3.0 * g_xxxx_xyy_1[i] * fe_0 + g_xxxx_xyyy_1[i] * pa_y[i];

        g_xxxxy_xyyz_0[i] = 2.0 * g_xxxx_xyz_1[i] * fe_0 + g_xxxx_xyyz_1[i] * pa_y[i];

        g_xxxxy_xyzz_0[i] = g_xxxx_xzz_1[i] * fe_0 + g_xxxx_xyzz_1[i] * pa_y[i];

        g_xxxxy_xzzz_0[i] = g_xxxx_xzzz_1[i] * pa_y[i];

        g_xxxxy_yyyy_0[i] = 4.0 * g_xxxx_yyy_1[i] * fe_0 + g_xxxx_yyyy_1[i] * pa_y[i];

        g_xxxxy_yyyz_0[i] = 3.0 * g_xxxx_yyz_1[i] * fe_0 + g_xxxx_yyyz_1[i] * pa_y[i];

        g_xxxxy_yyzz_0[i] = 2.0 * g_xxxx_yzz_1[i] * fe_0 + g_xxxx_yyzz_1[i] * pa_y[i];

        g_xxxxy_yzzz_0[i] = g_xxxx_zzz_1[i] * fe_0 + g_xxxx_yzzz_1[i] * pa_y[i];

        g_xxxxy_zzzz_0[i] = g_xxxx_zzzz_1[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : HG

    auto g_xxxxz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 30);

    auto g_xxxxz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 31);

    auto g_xxxxz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 32);

    auto g_xxxxz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 33);

    auto g_xxxxz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 34);

    auto g_xxxxz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 35);

    auto g_xxxxz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 36);

    auto g_xxxxz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 37);

    auto g_xxxxz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 38);

    auto g_xxxxz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 39);

    auto g_xxxxz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 40);

    auto g_xxxxz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 41);

    auto g_xxxxz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 42);

    auto g_xxxxz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 43);

    auto g_xxxxz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 44);

    #pragma omp simd aligned(g_xxxx_xxx_1, g_xxxx_xxxx_1, g_xxxx_xxxy_1, g_xxxx_xxxz_1, g_xxxx_xxy_1, g_xxxx_xxyy_1, g_xxxx_xxyz_1, g_xxxx_xxz_1, g_xxxx_xxzz_1, g_xxxx_xyy_1, g_xxxx_xyyy_1, g_xxxx_xyyz_1, g_xxxx_xyz_1, g_xxxx_xyzz_1, g_xxxx_xzz_1, g_xxxx_xzzz_1, g_xxxx_yyy_1, g_xxxx_yyyy_1, g_xxxx_yyyz_1, g_xxxx_yyz_1, g_xxxx_yyzz_1, g_xxxx_yzz_1, g_xxxx_yzzz_1, g_xxxx_zzz_1, g_xxxx_zzzz_1, g_xxxxz_xxxx_0, g_xxxxz_xxxy_0, g_xxxxz_xxxz_0, g_xxxxz_xxyy_0, g_xxxxz_xxyz_0, g_xxxxz_xxzz_0, g_xxxxz_xyyy_0, g_xxxxz_xyyz_0, g_xxxxz_xyzz_0, g_xxxxz_xzzz_0, g_xxxxz_yyyy_0, g_xxxxz_yyyz_0, g_xxxxz_yyzz_0, g_xxxxz_yzzz_0, g_xxxxz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxz_xxxx_0[i] = g_xxxx_xxxx_1[i] * pa_z[i];

        g_xxxxz_xxxy_0[i] = g_xxxx_xxxy_1[i] * pa_z[i];

        g_xxxxz_xxxz_0[i] = g_xxxx_xxx_1[i] * fe_0 + g_xxxx_xxxz_1[i] * pa_z[i];

        g_xxxxz_xxyy_0[i] = g_xxxx_xxyy_1[i] * pa_z[i];

        g_xxxxz_xxyz_0[i] = g_xxxx_xxy_1[i] * fe_0 + g_xxxx_xxyz_1[i] * pa_z[i];

        g_xxxxz_xxzz_0[i] = 2.0 * g_xxxx_xxz_1[i] * fe_0 + g_xxxx_xxzz_1[i] * pa_z[i];

        g_xxxxz_xyyy_0[i] = g_xxxx_xyyy_1[i] * pa_z[i];

        g_xxxxz_xyyz_0[i] = g_xxxx_xyy_1[i] * fe_0 + g_xxxx_xyyz_1[i] * pa_z[i];

        g_xxxxz_xyzz_0[i] = 2.0 * g_xxxx_xyz_1[i] * fe_0 + g_xxxx_xyzz_1[i] * pa_z[i];

        g_xxxxz_xzzz_0[i] = 3.0 * g_xxxx_xzz_1[i] * fe_0 + g_xxxx_xzzz_1[i] * pa_z[i];

        g_xxxxz_yyyy_0[i] = g_xxxx_yyyy_1[i] * pa_z[i];

        g_xxxxz_yyyz_0[i] = g_xxxx_yyy_1[i] * fe_0 + g_xxxx_yyyz_1[i] * pa_z[i];

        g_xxxxz_yyzz_0[i] = 2.0 * g_xxxx_yyz_1[i] * fe_0 + g_xxxx_yyzz_1[i] * pa_z[i];

        g_xxxxz_yzzz_0[i] = 3.0 * g_xxxx_yzz_1[i] * fe_0 + g_xxxx_yzzz_1[i] * pa_z[i];

        g_xxxxz_zzzz_0[i] = 4.0 * g_xxxx_zzz_1[i] * fe_0 + g_xxxx_zzzz_1[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : HG

    auto g_xxxyy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 45);

    auto g_xxxyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 46);

    auto g_xxxyy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 47);

    auto g_xxxyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 48);

    auto g_xxxyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 49);

    auto g_xxxyy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 50);

    auto g_xxxyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 51);

    auto g_xxxyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 52);

    auto g_xxxyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 53);

    auto g_xxxyy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 54);

    auto g_xxxyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 55);

    auto g_xxxyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 56);

    auto g_xxxyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 57);

    auto g_xxxyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 58);

    auto g_xxxyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 59);

    #pragma omp simd aligned(g_xxx_xxxx_0, g_xxx_xxxx_1, g_xxx_xxxz_0, g_xxx_xxxz_1, g_xxx_xxzz_0, g_xxx_xxzz_1, g_xxx_xzzz_0, g_xxx_xzzz_1, g_xxxy_xxxx_1, g_xxxy_xxxz_1, g_xxxy_xxzz_1, g_xxxy_xzzz_1, g_xxxyy_xxxx_0, g_xxxyy_xxxy_0, g_xxxyy_xxxz_0, g_xxxyy_xxyy_0, g_xxxyy_xxyz_0, g_xxxyy_xxzz_0, g_xxxyy_xyyy_0, g_xxxyy_xyyz_0, g_xxxyy_xyzz_0, g_xxxyy_xzzz_0, g_xxxyy_yyyy_0, g_xxxyy_yyyz_0, g_xxxyy_yyzz_0, g_xxxyy_yzzz_0, g_xxxyy_zzzz_0, g_xxyy_xxxy_1, g_xxyy_xxy_1, g_xxyy_xxyy_1, g_xxyy_xxyz_1, g_xxyy_xyy_1, g_xxyy_xyyy_1, g_xxyy_xyyz_1, g_xxyy_xyz_1, g_xxyy_xyzz_1, g_xxyy_yyy_1, g_xxyy_yyyy_1, g_xxyy_yyyz_1, g_xxyy_yyz_1, g_xxyy_yyzz_1, g_xxyy_yzz_1, g_xxyy_yzzz_1, g_xxyy_zzzz_1, g_xyy_xxxy_0, g_xyy_xxxy_1, g_xyy_xxyy_0, g_xyy_xxyy_1, g_xyy_xxyz_0, g_xyy_xxyz_1, g_xyy_xyyy_0, g_xyy_xyyy_1, g_xyy_xyyz_0, g_xyy_xyyz_1, g_xyy_xyzz_0, g_xyy_xyzz_1, g_xyy_yyyy_0, g_xyy_yyyy_1, g_xyy_yyyz_0, g_xyy_yyyz_1, g_xyy_yyzz_0, g_xyy_yyzz_1, g_xyy_yzzz_0, g_xyy_yzzz_1, g_xyy_zzzz_0, g_xyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyy_xxxx_0[i] = g_xxx_xxxx_0[i] * fbe_0 - g_xxx_xxxx_1[i] * fz_be_0 + g_xxxy_xxxx_1[i] * pa_y[i];

        g_xxxyy_xxxy_0[i] = 2.0 * g_xyy_xxxy_0[i] * fbe_0 - 2.0 * g_xyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xxyy_xxy_1[i] * fe_0 + g_xxyy_xxxy_1[i] * pa_x[i];

        g_xxxyy_xxxz_0[i] = g_xxx_xxxz_0[i] * fbe_0 - g_xxx_xxxz_1[i] * fz_be_0 + g_xxxy_xxxz_1[i] * pa_y[i];

        g_xxxyy_xxyy_0[i] = 2.0 * g_xyy_xxyy_0[i] * fbe_0 - 2.0 * g_xyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xxyy_xyy_1[i] * fe_0 + g_xxyy_xxyy_1[i] * pa_x[i];

        g_xxxyy_xxyz_0[i] = 2.0 * g_xyy_xxyz_0[i] * fbe_0 - 2.0 * g_xyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyy_xyz_1[i] * fe_0 + g_xxyy_xxyz_1[i] * pa_x[i];

        g_xxxyy_xxzz_0[i] = g_xxx_xxzz_0[i] * fbe_0 - g_xxx_xxzz_1[i] * fz_be_0 + g_xxxy_xxzz_1[i] * pa_y[i];

        g_xxxyy_xyyy_0[i] = 2.0 * g_xyy_xyyy_0[i] * fbe_0 - 2.0 * g_xyy_xyyy_1[i] * fz_be_0 + g_xxyy_yyy_1[i] * fe_0 + g_xxyy_xyyy_1[i] * pa_x[i];

        g_xxxyy_xyyz_0[i] = 2.0 * g_xyy_xyyz_0[i] * fbe_0 - 2.0 * g_xyy_xyyz_1[i] * fz_be_0 + g_xxyy_yyz_1[i] * fe_0 + g_xxyy_xyyz_1[i] * pa_x[i];

        g_xxxyy_xyzz_0[i] = 2.0 * g_xyy_xyzz_0[i] * fbe_0 - 2.0 * g_xyy_xyzz_1[i] * fz_be_0 + g_xxyy_yzz_1[i] * fe_0 + g_xxyy_xyzz_1[i] * pa_x[i];

        g_xxxyy_xzzz_0[i] = g_xxx_xzzz_0[i] * fbe_0 - g_xxx_xzzz_1[i] * fz_be_0 + g_xxxy_xzzz_1[i] * pa_y[i];

        g_xxxyy_yyyy_0[i] = 2.0 * g_xyy_yyyy_0[i] * fbe_0 - 2.0 * g_xyy_yyyy_1[i] * fz_be_0 + g_xxyy_yyyy_1[i] * pa_x[i];

        g_xxxyy_yyyz_0[i] = 2.0 * g_xyy_yyyz_0[i] * fbe_0 - 2.0 * g_xyy_yyyz_1[i] * fz_be_0 + g_xxyy_yyyz_1[i] * pa_x[i];

        g_xxxyy_yyzz_0[i] = 2.0 * g_xyy_yyzz_0[i] * fbe_0 - 2.0 * g_xyy_yyzz_1[i] * fz_be_0 + g_xxyy_yyzz_1[i] * pa_x[i];

        g_xxxyy_yzzz_0[i] = 2.0 * g_xyy_yzzz_0[i] * fbe_0 - 2.0 * g_xyy_yzzz_1[i] * fz_be_0 + g_xxyy_yzzz_1[i] * pa_x[i];

        g_xxxyy_zzzz_0[i] = 2.0 * g_xyy_zzzz_0[i] * fbe_0 - 2.0 * g_xyy_zzzz_1[i] * fz_be_0 + g_xxyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : HG

    auto g_xxxyz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 60);

    auto g_xxxyz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 61);

    auto g_xxxyz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 62);

    auto g_xxxyz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 63);

    auto g_xxxyz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 64);

    auto g_xxxyz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 65);

    auto g_xxxyz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 66);

    auto g_xxxyz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 67);

    auto g_xxxyz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 68);

    auto g_xxxyz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 69);

    auto g_xxxyz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 70);

    auto g_xxxyz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 71);

    auto g_xxxyz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 72);

    auto g_xxxyz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 73);

    auto g_xxxyz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 74);

    #pragma omp simd aligned(g_xxxy_xxxy_1, g_xxxy_xxyy_1, g_xxxy_xyyy_1, g_xxxy_yyyy_1, g_xxxyz_xxxx_0, g_xxxyz_xxxy_0, g_xxxyz_xxxz_0, g_xxxyz_xxyy_0, g_xxxyz_xxyz_0, g_xxxyz_xxzz_0, g_xxxyz_xyyy_0, g_xxxyz_xyyz_0, g_xxxyz_xyzz_0, g_xxxyz_xzzz_0, g_xxxyz_yyyy_0, g_xxxyz_yyyz_0, g_xxxyz_yyzz_0, g_xxxyz_yzzz_0, g_xxxyz_zzzz_0, g_xxxz_xxxx_1, g_xxxz_xxxz_1, g_xxxz_xxyz_1, g_xxxz_xxz_1, g_xxxz_xxzz_1, g_xxxz_xyyz_1, g_xxxz_xyz_1, g_xxxz_xyzz_1, g_xxxz_xzz_1, g_xxxz_xzzz_1, g_xxxz_yyyz_1, g_xxxz_yyz_1, g_xxxz_yyzz_1, g_xxxz_yzz_1, g_xxxz_yzzz_1, g_xxxz_zzz_1, g_xxxz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyz_xxxx_0[i] = g_xxxz_xxxx_1[i] * pa_y[i];

        g_xxxyz_xxxy_0[i] = g_xxxy_xxxy_1[i] * pa_z[i];

        g_xxxyz_xxxz_0[i] = g_xxxz_xxxz_1[i] * pa_y[i];

        g_xxxyz_xxyy_0[i] = g_xxxy_xxyy_1[i] * pa_z[i];

        g_xxxyz_xxyz_0[i] = g_xxxz_xxz_1[i] * fe_0 + g_xxxz_xxyz_1[i] * pa_y[i];

        g_xxxyz_xxzz_0[i] = g_xxxz_xxzz_1[i] * pa_y[i];

        g_xxxyz_xyyy_0[i] = g_xxxy_xyyy_1[i] * pa_z[i];

        g_xxxyz_xyyz_0[i] = 2.0 * g_xxxz_xyz_1[i] * fe_0 + g_xxxz_xyyz_1[i] * pa_y[i];

        g_xxxyz_xyzz_0[i] = g_xxxz_xzz_1[i] * fe_0 + g_xxxz_xyzz_1[i] * pa_y[i];

        g_xxxyz_xzzz_0[i] = g_xxxz_xzzz_1[i] * pa_y[i];

        g_xxxyz_yyyy_0[i] = g_xxxy_yyyy_1[i] * pa_z[i];

        g_xxxyz_yyyz_0[i] = 3.0 * g_xxxz_yyz_1[i] * fe_0 + g_xxxz_yyyz_1[i] * pa_y[i];

        g_xxxyz_yyzz_0[i] = 2.0 * g_xxxz_yzz_1[i] * fe_0 + g_xxxz_yyzz_1[i] * pa_y[i];

        g_xxxyz_yzzz_0[i] = g_xxxz_zzz_1[i] * fe_0 + g_xxxz_yzzz_1[i] * pa_y[i];

        g_xxxyz_zzzz_0[i] = g_xxxz_zzzz_1[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : HG

    auto g_xxxzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 75);

    auto g_xxxzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 76);

    auto g_xxxzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 77);

    auto g_xxxzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 78);

    auto g_xxxzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 79);

    auto g_xxxzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 80);

    auto g_xxxzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 81);

    auto g_xxxzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 82);

    auto g_xxxzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 83);

    auto g_xxxzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 84);

    auto g_xxxzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 85);

    auto g_xxxzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 86);

    auto g_xxxzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 87);

    auto g_xxxzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 88);

    auto g_xxxzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 89);

    #pragma omp simd aligned(g_xxx_xxxx_0, g_xxx_xxxx_1, g_xxx_xxxy_0, g_xxx_xxxy_1, g_xxx_xxyy_0, g_xxx_xxyy_1, g_xxx_xyyy_0, g_xxx_xyyy_1, g_xxxz_xxxx_1, g_xxxz_xxxy_1, g_xxxz_xxyy_1, g_xxxz_xyyy_1, g_xxxzz_xxxx_0, g_xxxzz_xxxy_0, g_xxxzz_xxxz_0, g_xxxzz_xxyy_0, g_xxxzz_xxyz_0, g_xxxzz_xxzz_0, g_xxxzz_xyyy_0, g_xxxzz_xyyz_0, g_xxxzz_xyzz_0, g_xxxzz_xzzz_0, g_xxxzz_yyyy_0, g_xxxzz_yyyz_0, g_xxxzz_yyzz_0, g_xxxzz_yzzz_0, g_xxxzz_zzzz_0, g_xxzz_xxxz_1, g_xxzz_xxyz_1, g_xxzz_xxz_1, g_xxzz_xxzz_1, g_xxzz_xyyz_1, g_xxzz_xyz_1, g_xxzz_xyzz_1, g_xxzz_xzz_1, g_xxzz_xzzz_1, g_xxzz_yyyy_1, g_xxzz_yyyz_1, g_xxzz_yyz_1, g_xxzz_yyzz_1, g_xxzz_yzz_1, g_xxzz_yzzz_1, g_xxzz_zzz_1, g_xxzz_zzzz_1, g_xzz_xxxz_0, g_xzz_xxxz_1, g_xzz_xxyz_0, g_xzz_xxyz_1, g_xzz_xxzz_0, g_xzz_xxzz_1, g_xzz_xyyz_0, g_xzz_xyyz_1, g_xzz_xyzz_0, g_xzz_xyzz_1, g_xzz_xzzz_0, g_xzz_xzzz_1, g_xzz_yyyy_0, g_xzz_yyyy_1, g_xzz_yyyz_0, g_xzz_yyyz_1, g_xzz_yyzz_0, g_xzz_yyzz_1, g_xzz_yzzz_0, g_xzz_yzzz_1, g_xzz_zzzz_0, g_xzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzz_xxxx_0[i] = g_xxx_xxxx_0[i] * fbe_0 - g_xxx_xxxx_1[i] * fz_be_0 + g_xxxz_xxxx_1[i] * pa_z[i];

        g_xxxzz_xxxy_0[i] = g_xxx_xxxy_0[i] * fbe_0 - g_xxx_xxxy_1[i] * fz_be_0 + g_xxxz_xxxy_1[i] * pa_z[i];

        g_xxxzz_xxxz_0[i] = 2.0 * g_xzz_xxxz_0[i] * fbe_0 - 2.0 * g_xzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xxzz_xxz_1[i] * fe_0 + g_xxzz_xxxz_1[i] * pa_x[i];

        g_xxxzz_xxyy_0[i] = g_xxx_xxyy_0[i] * fbe_0 - g_xxx_xxyy_1[i] * fz_be_0 + g_xxxz_xxyy_1[i] * pa_z[i];

        g_xxxzz_xxyz_0[i] = 2.0 * g_xzz_xxyz_0[i] * fbe_0 - 2.0 * g_xzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xxzz_xyz_1[i] * fe_0 + g_xxzz_xxyz_1[i] * pa_x[i];

        g_xxxzz_xxzz_0[i] = 2.0 * g_xzz_xxzz_0[i] * fbe_0 - 2.0 * g_xzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xxzz_xzz_1[i] * fe_0 + g_xxzz_xxzz_1[i] * pa_x[i];

        g_xxxzz_xyyy_0[i] = g_xxx_xyyy_0[i] * fbe_0 - g_xxx_xyyy_1[i] * fz_be_0 + g_xxxz_xyyy_1[i] * pa_z[i];

        g_xxxzz_xyyz_0[i] = 2.0 * g_xzz_xyyz_0[i] * fbe_0 - 2.0 * g_xzz_xyyz_1[i] * fz_be_0 + g_xxzz_yyz_1[i] * fe_0 + g_xxzz_xyyz_1[i] * pa_x[i];

        g_xxxzz_xyzz_0[i] = 2.0 * g_xzz_xyzz_0[i] * fbe_0 - 2.0 * g_xzz_xyzz_1[i] * fz_be_0 + g_xxzz_yzz_1[i] * fe_0 + g_xxzz_xyzz_1[i] * pa_x[i];

        g_xxxzz_xzzz_0[i] = 2.0 * g_xzz_xzzz_0[i] * fbe_0 - 2.0 * g_xzz_xzzz_1[i] * fz_be_0 + g_xxzz_zzz_1[i] * fe_0 + g_xxzz_xzzz_1[i] * pa_x[i];

        g_xxxzz_yyyy_0[i] = 2.0 * g_xzz_yyyy_0[i] * fbe_0 - 2.0 * g_xzz_yyyy_1[i] * fz_be_0 + g_xxzz_yyyy_1[i] * pa_x[i];

        g_xxxzz_yyyz_0[i] = 2.0 * g_xzz_yyyz_0[i] * fbe_0 - 2.0 * g_xzz_yyyz_1[i] * fz_be_0 + g_xxzz_yyyz_1[i] * pa_x[i];

        g_xxxzz_yyzz_0[i] = 2.0 * g_xzz_yyzz_0[i] * fbe_0 - 2.0 * g_xzz_yyzz_1[i] * fz_be_0 + g_xxzz_yyzz_1[i] * pa_x[i];

        g_xxxzz_yzzz_0[i] = 2.0 * g_xzz_yzzz_0[i] * fbe_0 - 2.0 * g_xzz_yzzz_1[i] * fz_be_0 + g_xxzz_yzzz_1[i] * pa_x[i];

        g_xxxzz_zzzz_0[i] = 2.0 * g_xzz_zzzz_0[i] * fbe_0 - 2.0 * g_xzz_zzzz_1[i] * fz_be_0 + g_xxzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : HG

    auto g_xxyyy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 90);

    auto g_xxyyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 91);

    auto g_xxyyy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 92);

    auto g_xxyyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 93);

    auto g_xxyyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 94);

    auto g_xxyyy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 95);

    auto g_xxyyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 96);

    auto g_xxyyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 97);

    auto g_xxyyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 98);

    auto g_xxyyy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 99);

    auto g_xxyyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 100);

    auto g_xxyyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 101);

    auto g_xxyyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 102);

    auto g_xxyyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 103);

    auto g_xxyyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 104);

    #pragma omp simd aligned(g_xxy_xxxx_0, g_xxy_xxxx_1, g_xxy_xxxz_0, g_xxy_xxxz_1, g_xxy_xxzz_0, g_xxy_xxzz_1, g_xxy_xzzz_0, g_xxy_xzzz_1, g_xxyy_xxxx_1, g_xxyy_xxxz_1, g_xxyy_xxzz_1, g_xxyy_xzzz_1, g_xxyyy_xxxx_0, g_xxyyy_xxxy_0, g_xxyyy_xxxz_0, g_xxyyy_xxyy_0, g_xxyyy_xxyz_0, g_xxyyy_xxzz_0, g_xxyyy_xyyy_0, g_xxyyy_xyyz_0, g_xxyyy_xyzz_0, g_xxyyy_xzzz_0, g_xxyyy_yyyy_0, g_xxyyy_yyyz_0, g_xxyyy_yyzz_0, g_xxyyy_yzzz_0, g_xxyyy_zzzz_0, g_xyyy_xxxy_1, g_xyyy_xxy_1, g_xyyy_xxyy_1, g_xyyy_xxyz_1, g_xyyy_xyy_1, g_xyyy_xyyy_1, g_xyyy_xyyz_1, g_xyyy_xyz_1, g_xyyy_xyzz_1, g_xyyy_yyy_1, g_xyyy_yyyy_1, g_xyyy_yyyz_1, g_xyyy_yyz_1, g_xyyy_yyzz_1, g_xyyy_yzz_1, g_xyyy_yzzz_1, g_xyyy_zzzz_1, g_yyy_xxxy_0, g_yyy_xxxy_1, g_yyy_xxyy_0, g_yyy_xxyy_1, g_yyy_xxyz_0, g_yyy_xxyz_1, g_yyy_xyyy_0, g_yyy_xyyy_1, g_yyy_xyyz_0, g_yyy_xyyz_1, g_yyy_xyzz_0, g_yyy_xyzz_1, g_yyy_yyyy_0, g_yyy_yyyy_1, g_yyy_yyyz_0, g_yyy_yyyz_1, g_yyy_yyzz_0, g_yyy_yyzz_1, g_yyy_yzzz_0, g_yyy_yzzz_1, g_yyy_zzzz_0, g_yyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyy_xxxx_0[i] = 2.0 * g_xxy_xxxx_0[i] * fbe_0 - 2.0 * g_xxy_xxxx_1[i] * fz_be_0 + g_xxyy_xxxx_1[i] * pa_y[i];

        g_xxyyy_xxxy_0[i] = g_yyy_xxxy_0[i] * fbe_0 - g_yyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xyyy_xxy_1[i] * fe_0 + g_xyyy_xxxy_1[i] * pa_x[i];

        g_xxyyy_xxxz_0[i] = 2.0 * g_xxy_xxxz_0[i] * fbe_0 - 2.0 * g_xxy_xxxz_1[i] * fz_be_0 + g_xxyy_xxxz_1[i] * pa_y[i];

        g_xxyyy_xxyy_0[i] = g_yyy_xxyy_0[i] * fbe_0 - g_yyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xyyy_xyy_1[i] * fe_0 + g_xyyy_xxyy_1[i] * pa_x[i];

        g_xxyyy_xxyz_0[i] = g_yyy_xxyz_0[i] * fbe_0 - g_yyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyy_xyz_1[i] * fe_0 + g_xyyy_xxyz_1[i] * pa_x[i];

        g_xxyyy_xxzz_0[i] = 2.0 * g_xxy_xxzz_0[i] * fbe_0 - 2.0 * g_xxy_xxzz_1[i] * fz_be_0 + g_xxyy_xxzz_1[i] * pa_y[i];

        g_xxyyy_xyyy_0[i] = g_yyy_xyyy_0[i] * fbe_0 - g_yyy_xyyy_1[i] * fz_be_0 + g_xyyy_yyy_1[i] * fe_0 + g_xyyy_xyyy_1[i] * pa_x[i];

        g_xxyyy_xyyz_0[i] = g_yyy_xyyz_0[i] * fbe_0 - g_yyy_xyyz_1[i] * fz_be_0 + g_xyyy_yyz_1[i] * fe_0 + g_xyyy_xyyz_1[i] * pa_x[i];

        g_xxyyy_xyzz_0[i] = g_yyy_xyzz_0[i] * fbe_0 - g_yyy_xyzz_1[i] * fz_be_0 + g_xyyy_yzz_1[i] * fe_0 + g_xyyy_xyzz_1[i] * pa_x[i];

        g_xxyyy_xzzz_0[i] = 2.0 * g_xxy_xzzz_0[i] * fbe_0 - 2.0 * g_xxy_xzzz_1[i] * fz_be_0 + g_xxyy_xzzz_1[i] * pa_y[i];

        g_xxyyy_yyyy_0[i] = g_yyy_yyyy_0[i] * fbe_0 - g_yyy_yyyy_1[i] * fz_be_0 + g_xyyy_yyyy_1[i] * pa_x[i];

        g_xxyyy_yyyz_0[i] = g_yyy_yyyz_0[i] * fbe_0 - g_yyy_yyyz_1[i] * fz_be_0 + g_xyyy_yyyz_1[i] * pa_x[i];

        g_xxyyy_yyzz_0[i] = g_yyy_yyzz_0[i] * fbe_0 - g_yyy_yyzz_1[i] * fz_be_0 + g_xyyy_yyzz_1[i] * pa_x[i];

        g_xxyyy_yzzz_0[i] = g_yyy_yzzz_0[i] * fbe_0 - g_yyy_yzzz_1[i] * fz_be_0 + g_xyyy_yzzz_1[i] * pa_x[i];

        g_xxyyy_zzzz_0[i] = g_yyy_zzzz_0[i] * fbe_0 - g_yyy_zzzz_1[i] * fz_be_0 + g_xyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : HG

    auto g_xxyyz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 105);

    auto g_xxyyz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 106);

    auto g_xxyyz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 107);

    auto g_xxyyz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 108);

    auto g_xxyyz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 109);

    auto g_xxyyz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 110);

    auto g_xxyyz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 111);

    auto g_xxyyz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 112);

    auto g_xxyyz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 113);

    auto g_xxyyz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 114);

    auto g_xxyyz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 115);

    auto g_xxyyz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 116);

    auto g_xxyyz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 117);

    auto g_xxyyz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 118);

    auto g_xxyyz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 119);

    #pragma omp simd aligned(g_xxyy_xxx_1, g_xxyy_xxxx_1, g_xxyy_xxxy_1, g_xxyy_xxxz_1, g_xxyy_xxy_1, g_xxyy_xxyy_1, g_xxyy_xxyz_1, g_xxyy_xxz_1, g_xxyy_xxzz_1, g_xxyy_xyy_1, g_xxyy_xyyy_1, g_xxyy_xyyz_1, g_xxyy_xyz_1, g_xxyy_xyzz_1, g_xxyy_xzz_1, g_xxyy_xzzz_1, g_xxyy_yyy_1, g_xxyy_yyyy_1, g_xxyy_yyyz_1, g_xxyy_yyz_1, g_xxyy_yyzz_1, g_xxyy_yzz_1, g_xxyy_yzzz_1, g_xxyy_zzz_1, g_xxyy_zzzz_1, g_xxyyz_xxxx_0, g_xxyyz_xxxy_0, g_xxyyz_xxxz_0, g_xxyyz_xxyy_0, g_xxyyz_xxyz_0, g_xxyyz_xxzz_0, g_xxyyz_xyyy_0, g_xxyyz_xyyz_0, g_xxyyz_xyzz_0, g_xxyyz_xzzz_0, g_xxyyz_yyyy_0, g_xxyyz_yyyz_0, g_xxyyz_yyzz_0, g_xxyyz_yzzz_0, g_xxyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyz_xxxx_0[i] = g_xxyy_xxxx_1[i] * pa_z[i];

        g_xxyyz_xxxy_0[i] = g_xxyy_xxxy_1[i] * pa_z[i];

        g_xxyyz_xxxz_0[i] = g_xxyy_xxx_1[i] * fe_0 + g_xxyy_xxxz_1[i] * pa_z[i];

        g_xxyyz_xxyy_0[i] = g_xxyy_xxyy_1[i] * pa_z[i];

        g_xxyyz_xxyz_0[i] = g_xxyy_xxy_1[i] * fe_0 + g_xxyy_xxyz_1[i] * pa_z[i];

        g_xxyyz_xxzz_0[i] = 2.0 * g_xxyy_xxz_1[i] * fe_0 + g_xxyy_xxzz_1[i] * pa_z[i];

        g_xxyyz_xyyy_0[i] = g_xxyy_xyyy_1[i] * pa_z[i];

        g_xxyyz_xyyz_0[i] = g_xxyy_xyy_1[i] * fe_0 + g_xxyy_xyyz_1[i] * pa_z[i];

        g_xxyyz_xyzz_0[i] = 2.0 * g_xxyy_xyz_1[i] * fe_0 + g_xxyy_xyzz_1[i] * pa_z[i];

        g_xxyyz_xzzz_0[i] = 3.0 * g_xxyy_xzz_1[i] * fe_0 + g_xxyy_xzzz_1[i] * pa_z[i];

        g_xxyyz_yyyy_0[i] = g_xxyy_yyyy_1[i] * pa_z[i];

        g_xxyyz_yyyz_0[i] = g_xxyy_yyy_1[i] * fe_0 + g_xxyy_yyyz_1[i] * pa_z[i];

        g_xxyyz_yyzz_0[i] = 2.0 * g_xxyy_yyz_1[i] * fe_0 + g_xxyy_yyzz_1[i] * pa_z[i];

        g_xxyyz_yzzz_0[i] = 3.0 * g_xxyy_yzz_1[i] * fe_0 + g_xxyy_yzzz_1[i] * pa_z[i];

        g_xxyyz_zzzz_0[i] = 4.0 * g_xxyy_zzz_1[i] * fe_0 + g_xxyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 120-135 components of targeted buffer : HG

    auto g_xxyzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 120);

    auto g_xxyzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 121);

    auto g_xxyzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 122);

    auto g_xxyzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 123);

    auto g_xxyzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 124);

    auto g_xxyzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 125);

    auto g_xxyzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 126);

    auto g_xxyzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 127);

    auto g_xxyzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 128);

    auto g_xxyzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 129);

    auto g_xxyzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 130);

    auto g_xxyzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 131);

    auto g_xxyzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 132);

    auto g_xxyzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 133);

    auto g_xxyzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 134);

    #pragma omp simd aligned(g_xxyzz_xxxx_0, g_xxyzz_xxxy_0, g_xxyzz_xxxz_0, g_xxyzz_xxyy_0, g_xxyzz_xxyz_0, g_xxyzz_xxzz_0, g_xxyzz_xyyy_0, g_xxyzz_xyyz_0, g_xxyzz_xyzz_0, g_xxyzz_xzzz_0, g_xxyzz_yyyy_0, g_xxyzz_yyyz_0, g_xxyzz_yyzz_0, g_xxyzz_yzzz_0, g_xxyzz_zzzz_0, g_xxzz_xxx_1, g_xxzz_xxxx_1, g_xxzz_xxxy_1, g_xxzz_xxxz_1, g_xxzz_xxy_1, g_xxzz_xxyy_1, g_xxzz_xxyz_1, g_xxzz_xxz_1, g_xxzz_xxzz_1, g_xxzz_xyy_1, g_xxzz_xyyy_1, g_xxzz_xyyz_1, g_xxzz_xyz_1, g_xxzz_xyzz_1, g_xxzz_xzz_1, g_xxzz_xzzz_1, g_xxzz_yyy_1, g_xxzz_yyyy_1, g_xxzz_yyyz_1, g_xxzz_yyz_1, g_xxzz_yyzz_1, g_xxzz_yzz_1, g_xxzz_yzzz_1, g_xxzz_zzz_1, g_xxzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzz_xxxx_0[i] = g_xxzz_xxxx_1[i] * pa_y[i];

        g_xxyzz_xxxy_0[i] = g_xxzz_xxx_1[i] * fe_0 + g_xxzz_xxxy_1[i] * pa_y[i];

        g_xxyzz_xxxz_0[i] = g_xxzz_xxxz_1[i] * pa_y[i];

        g_xxyzz_xxyy_0[i] = 2.0 * g_xxzz_xxy_1[i] * fe_0 + g_xxzz_xxyy_1[i] * pa_y[i];

        g_xxyzz_xxyz_0[i] = g_xxzz_xxz_1[i] * fe_0 + g_xxzz_xxyz_1[i] * pa_y[i];

        g_xxyzz_xxzz_0[i] = g_xxzz_xxzz_1[i] * pa_y[i];

        g_xxyzz_xyyy_0[i] = 3.0 * g_xxzz_xyy_1[i] * fe_0 + g_xxzz_xyyy_1[i] * pa_y[i];

        g_xxyzz_xyyz_0[i] = 2.0 * g_xxzz_xyz_1[i] * fe_0 + g_xxzz_xyyz_1[i] * pa_y[i];

        g_xxyzz_xyzz_0[i] = g_xxzz_xzz_1[i] * fe_0 + g_xxzz_xyzz_1[i] * pa_y[i];

        g_xxyzz_xzzz_0[i] = g_xxzz_xzzz_1[i] * pa_y[i];

        g_xxyzz_yyyy_0[i] = 4.0 * g_xxzz_yyy_1[i] * fe_0 + g_xxzz_yyyy_1[i] * pa_y[i];

        g_xxyzz_yyyz_0[i] = 3.0 * g_xxzz_yyz_1[i] * fe_0 + g_xxzz_yyyz_1[i] * pa_y[i];

        g_xxyzz_yyzz_0[i] = 2.0 * g_xxzz_yzz_1[i] * fe_0 + g_xxzz_yyzz_1[i] * pa_y[i];

        g_xxyzz_yzzz_0[i] = g_xxzz_zzz_1[i] * fe_0 + g_xxzz_yzzz_1[i] * pa_y[i];

        g_xxyzz_zzzz_0[i] = g_xxzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : HG

    auto g_xxzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 135);

    auto g_xxzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 136);

    auto g_xxzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 137);

    auto g_xxzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 138);

    auto g_xxzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 139);

    auto g_xxzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 140);

    auto g_xxzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 141);

    auto g_xxzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 142);

    auto g_xxzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 143);

    auto g_xxzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 144);

    auto g_xxzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 145);

    auto g_xxzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 146);

    auto g_xxzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 147);

    auto g_xxzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 148);

    auto g_xxzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 149);

    #pragma omp simd aligned(g_xxz_xxxx_0, g_xxz_xxxx_1, g_xxz_xxxy_0, g_xxz_xxxy_1, g_xxz_xxyy_0, g_xxz_xxyy_1, g_xxz_xyyy_0, g_xxz_xyyy_1, g_xxzz_xxxx_1, g_xxzz_xxxy_1, g_xxzz_xxyy_1, g_xxzz_xyyy_1, g_xxzzz_xxxx_0, g_xxzzz_xxxy_0, g_xxzzz_xxxz_0, g_xxzzz_xxyy_0, g_xxzzz_xxyz_0, g_xxzzz_xxzz_0, g_xxzzz_xyyy_0, g_xxzzz_xyyz_0, g_xxzzz_xyzz_0, g_xxzzz_xzzz_0, g_xxzzz_yyyy_0, g_xxzzz_yyyz_0, g_xxzzz_yyzz_0, g_xxzzz_yzzz_0, g_xxzzz_zzzz_0, g_xzzz_xxxz_1, g_xzzz_xxyz_1, g_xzzz_xxz_1, g_xzzz_xxzz_1, g_xzzz_xyyz_1, g_xzzz_xyz_1, g_xzzz_xyzz_1, g_xzzz_xzz_1, g_xzzz_xzzz_1, g_xzzz_yyyy_1, g_xzzz_yyyz_1, g_xzzz_yyz_1, g_xzzz_yyzz_1, g_xzzz_yzz_1, g_xzzz_yzzz_1, g_xzzz_zzz_1, g_xzzz_zzzz_1, g_zzz_xxxz_0, g_zzz_xxxz_1, g_zzz_xxyz_0, g_zzz_xxyz_1, g_zzz_xxzz_0, g_zzz_xxzz_1, g_zzz_xyyz_0, g_zzz_xyyz_1, g_zzz_xyzz_0, g_zzz_xyzz_1, g_zzz_xzzz_0, g_zzz_xzzz_1, g_zzz_yyyy_0, g_zzz_yyyy_1, g_zzz_yyyz_0, g_zzz_yyyz_1, g_zzz_yyzz_0, g_zzz_yyzz_1, g_zzz_yzzz_0, g_zzz_yzzz_1, g_zzz_zzzz_0, g_zzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzz_xxxx_0[i] = 2.0 * g_xxz_xxxx_0[i] * fbe_0 - 2.0 * g_xxz_xxxx_1[i] * fz_be_0 + g_xxzz_xxxx_1[i] * pa_z[i];

        g_xxzzz_xxxy_0[i] = 2.0 * g_xxz_xxxy_0[i] * fbe_0 - 2.0 * g_xxz_xxxy_1[i] * fz_be_0 + g_xxzz_xxxy_1[i] * pa_z[i];

        g_xxzzz_xxxz_0[i] = g_zzz_xxxz_0[i] * fbe_0 - g_zzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xzzz_xxz_1[i] * fe_0 + g_xzzz_xxxz_1[i] * pa_x[i];

        g_xxzzz_xxyy_0[i] = 2.0 * g_xxz_xxyy_0[i] * fbe_0 - 2.0 * g_xxz_xxyy_1[i] * fz_be_0 + g_xxzz_xxyy_1[i] * pa_z[i];

        g_xxzzz_xxyz_0[i] = g_zzz_xxyz_0[i] * fbe_0 - g_zzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xzzz_xyz_1[i] * fe_0 + g_xzzz_xxyz_1[i] * pa_x[i];

        g_xxzzz_xxzz_0[i] = g_zzz_xxzz_0[i] * fbe_0 - g_zzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xzzz_xzz_1[i] * fe_0 + g_xzzz_xxzz_1[i] * pa_x[i];

        g_xxzzz_xyyy_0[i] = 2.0 * g_xxz_xyyy_0[i] * fbe_0 - 2.0 * g_xxz_xyyy_1[i] * fz_be_0 + g_xxzz_xyyy_1[i] * pa_z[i];

        g_xxzzz_xyyz_0[i] = g_zzz_xyyz_0[i] * fbe_0 - g_zzz_xyyz_1[i] * fz_be_0 + g_xzzz_yyz_1[i] * fe_0 + g_xzzz_xyyz_1[i] * pa_x[i];

        g_xxzzz_xyzz_0[i] = g_zzz_xyzz_0[i] * fbe_0 - g_zzz_xyzz_1[i] * fz_be_0 + g_xzzz_yzz_1[i] * fe_0 + g_xzzz_xyzz_1[i] * pa_x[i];

        g_xxzzz_xzzz_0[i] = g_zzz_xzzz_0[i] * fbe_0 - g_zzz_xzzz_1[i] * fz_be_0 + g_xzzz_zzz_1[i] * fe_0 + g_xzzz_xzzz_1[i] * pa_x[i];

        g_xxzzz_yyyy_0[i] = g_zzz_yyyy_0[i] * fbe_0 - g_zzz_yyyy_1[i] * fz_be_0 + g_xzzz_yyyy_1[i] * pa_x[i];

        g_xxzzz_yyyz_0[i] = g_zzz_yyyz_0[i] * fbe_0 - g_zzz_yyyz_1[i] * fz_be_0 + g_xzzz_yyyz_1[i] * pa_x[i];

        g_xxzzz_yyzz_0[i] = g_zzz_yyzz_0[i] * fbe_0 - g_zzz_yyzz_1[i] * fz_be_0 + g_xzzz_yyzz_1[i] * pa_x[i];

        g_xxzzz_yzzz_0[i] = g_zzz_yzzz_0[i] * fbe_0 - g_zzz_yzzz_1[i] * fz_be_0 + g_xzzz_yzzz_1[i] * pa_x[i];

        g_xxzzz_zzzz_0[i] = g_zzz_zzzz_0[i] * fbe_0 - g_zzz_zzzz_1[i] * fz_be_0 + g_xzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : HG

    auto g_xyyyy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 150);

    auto g_xyyyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 151);

    auto g_xyyyy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 152);

    auto g_xyyyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 153);

    auto g_xyyyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 154);

    auto g_xyyyy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 155);

    auto g_xyyyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 156);

    auto g_xyyyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 157);

    auto g_xyyyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 158);

    auto g_xyyyy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 159);

    auto g_xyyyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 160);

    auto g_xyyyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 161);

    auto g_xyyyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 162);

    auto g_xyyyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 163);

    auto g_xyyyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 164);

    #pragma omp simd aligned(g_xyyyy_xxxx_0, g_xyyyy_xxxy_0, g_xyyyy_xxxz_0, g_xyyyy_xxyy_0, g_xyyyy_xxyz_0, g_xyyyy_xxzz_0, g_xyyyy_xyyy_0, g_xyyyy_xyyz_0, g_xyyyy_xyzz_0, g_xyyyy_xzzz_0, g_xyyyy_yyyy_0, g_xyyyy_yyyz_0, g_xyyyy_yyzz_0, g_xyyyy_yzzz_0, g_xyyyy_zzzz_0, g_yyyy_xxx_1, g_yyyy_xxxx_1, g_yyyy_xxxy_1, g_yyyy_xxxz_1, g_yyyy_xxy_1, g_yyyy_xxyy_1, g_yyyy_xxyz_1, g_yyyy_xxz_1, g_yyyy_xxzz_1, g_yyyy_xyy_1, g_yyyy_xyyy_1, g_yyyy_xyyz_1, g_yyyy_xyz_1, g_yyyy_xyzz_1, g_yyyy_xzz_1, g_yyyy_xzzz_1, g_yyyy_yyy_1, g_yyyy_yyyy_1, g_yyyy_yyyz_1, g_yyyy_yyz_1, g_yyyy_yyzz_1, g_yyyy_yzz_1, g_yyyy_yzzz_1, g_yyyy_zzz_1, g_yyyy_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyy_xxxx_0[i] = 4.0 * g_yyyy_xxx_1[i] * fe_0 + g_yyyy_xxxx_1[i] * pa_x[i];

        g_xyyyy_xxxy_0[i] = 3.0 * g_yyyy_xxy_1[i] * fe_0 + g_yyyy_xxxy_1[i] * pa_x[i];

        g_xyyyy_xxxz_0[i] = 3.0 * g_yyyy_xxz_1[i] * fe_0 + g_yyyy_xxxz_1[i] * pa_x[i];

        g_xyyyy_xxyy_0[i] = 2.0 * g_yyyy_xyy_1[i] * fe_0 + g_yyyy_xxyy_1[i] * pa_x[i];

        g_xyyyy_xxyz_0[i] = 2.0 * g_yyyy_xyz_1[i] * fe_0 + g_yyyy_xxyz_1[i] * pa_x[i];

        g_xyyyy_xxzz_0[i] = 2.0 * g_yyyy_xzz_1[i] * fe_0 + g_yyyy_xxzz_1[i] * pa_x[i];

        g_xyyyy_xyyy_0[i] = g_yyyy_yyy_1[i] * fe_0 + g_yyyy_xyyy_1[i] * pa_x[i];

        g_xyyyy_xyyz_0[i] = g_yyyy_yyz_1[i] * fe_0 + g_yyyy_xyyz_1[i] * pa_x[i];

        g_xyyyy_xyzz_0[i] = g_yyyy_yzz_1[i] * fe_0 + g_yyyy_xyzz_1[i] * pa_x[i];

        g_xyyyy_xzzz_0[i] = g_yyyy_zzz_1[i] * fe_0 + g_yyyy_xzzz_1[i] * pa_x[i];

        g_xyyyy_yyyy_0[i] = g_yyyy_yyyy_1[i] * pa_x[i];

        g_xyyyy_yyyz_0[i] = g_yyyy_yyyz_1[i] * pa_x[i];

        g_xyyyy_yyzz_0[i] = g_yyyy_yyzz_1[i] * pa_x[i];

        g_xyyyy_yzzz_0[i] = g_yyyy_yzzz_1[i] * pa_x[i];

        g_xyyyy_zzzz_0[i] = g_yyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 165-180 components of targeted buffer : HG

    auto g_xyyyz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 165);

    auto g_xyyyz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 166);

    auto g_xyyyz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 167);

    auto g_xyyyz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 168);

    auto g_xyyyz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 169);

    auto g_xyyyz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 170);

    auto g_xyyyz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 171);

    auto g_xyyyz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 172);

    auto g_xyyyz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 173);

    auto g_xyyyz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 174);

    auto g_xyyyz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 175);

    auto g_xyyyz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 176);

    auto g_xyyyz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 177);

    auto g_xyyyz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 178);

    auto g_xyyyz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 179);

    #pragma omp simd aligned(g_xyyy_xxxx_1, g_xyyy_xxxy_1, g_xyyy_xxyy_1, g_xyyy_xyyy_1, g_xyyyz_xxxx_0, g_xyyyz_xxxy_0, g_xyyyz_xxxz_0, g_xyyyz_xxyy_0, g_xyyyz_xxyz_0, g_xyyyz_xxzz_0, g_xyyyz_xyyy_0, g_xyyyz_xyyz_0, g_xyyyz_xyzz_0, g_xyyyz_xzzz_0, g_xyyyz_yyyy_0, g_xyyyz_yyyz_0, g_xyyyz_yyzz_0, g_xyyyz_yzzz_0, g_xyyyz_zzzz_0, g_yyyz_xxxz_1, g_yyyz_xxyz_1, g_yyyz_xxz_1, g_yyyz_xxzz_1, g_yyyz_xyyz_1, g_yyyz_xyz_1, g_yyyz_xyzz_1, g_yyyz_xzz_1, g_yyyz_xzzz_1, g_yyyz_yyyy_1, g_yyyz_yyyz_1, g_yyyz_yyz_1, g_yyyz_yyzz_1, g_yyyz_yzz_1, g_yyyz_yzzz_1, g_yyyz_zzz_1, g_yyyz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyz_xxxx_0[i] = g_xyyy_xxxx_1[i] * pa_z[i];

        g_xyyyz_xxxy_0[i] = g_xyyy_xxxy_1[i] * pa_z[i];

        g_xyyyz_xxxz_0[i] = 3.0 * g_yyyz_xxz_1[i] * fe_0 + g_yyyz_xxxz_1[i] * pa_x[i];

        g_xyyyz_xxyy_0[i] = g_xyyy_xxyy_1[i] * pa_z[i];

        g_xyyyz_xxyz_0[i] = 2.0 * g_yyyz_xyz_1[i] * fe_0 + g_yyyz_xxyz_1[i] * pa_x[i];

        g_xyyyz_xxzz_0[i] = 2.0 * g_yyyz_xzz_1[i] * fe_0 + g_yyyz_xxzz_1[i] * pa_x[i];

        g_xyyyz_xyyy_0[i] = g_xyyy_xyyy_1[i] * pa_z[i];

        g_xyyyz_xyyz_0[i] = g_yyyz_yyz_1[i] * fe_0 + g_yyyz_xyyz_1[i] * pa_x[i];

        g_xyyyz_xyzz_0[i] = g_yyyz_yzz_1[i] * fe_0 + g_yyyz_xyzz_1[i] * pa_x[i];

        g_xyyyz_xzzz_0[i] = g_yyyz_zzz_1[i] * fe_0 + g_yyyz_xzzz_1[i] * pa_x[i];

        g_xyyyz_yyyy_0[i] = g_yyyz_yyyy_1[i] * pa_x[i];

        g_xyyyz_yyyz_0[i] = g_yyyz_yyyz_1[i] * pa_x[i];

        g_xyyyz_yyzz_0[i] = g_yyyz_yyzz_1[i] * pa_x[i];

        g_xyyyz_yzzz_0[i] = g_yyyz_yzzz_1[i] * pa_x[i];

        g_xyyyz_zzzz_0[i] = g_yyyz_zzzz_1[i] * pa_x[i];
    }

    // Set up 180-195 components of targeted buffer : HG

    auto g_xyyzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 180);

    auto g_xyyzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 181);

    auto g_xyyzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 182);

    auto g_xyyzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 183);

    auto g_xyyzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 184);

    auto g_xyyzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 185);

    auto g_xyyzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 186);

    auto g_xyyzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 187);

    auto g_xyyzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 188);

    auto g_xyyzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 189);

    auto g_xyyzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 190);

    auto g_xyyzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 191);

    auto g_xyyzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 192);

    auto g_xyyzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 193);

    auto g_xyyzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 194);

    #pragma omp simd aligned(g_xyyzz_xxxx_0, g_xyyzz_xxxy_0, g_xyyzz_xxxz_0, g_xyyzz_xxyy_0, g_xyyzz_xxyz_0, g_xyyzz_xxzz_0, g_xyyzz_xyyy_0, g_xyyzz_xyyz_0, g_xyyzz_xyzz_0, g_xyyzz_xzzz_0, g_xyyzz_yyyy_0, g_xyyzz_yyyz_0, g_xyyzz_yyzz_0, g_xyyzz_yzzz_0, g_xyyzz_zzzz_0, g_yyzz_xxx_1, g_yyzz_xxxx_1, g_yyzz_xxxy_1, g_yyzz_xxxz_1, g_yyzz_xxy_1, g_yyzz_xxyy_1, g_yyzz_xxyz_1, g_yyzz_xxz_1, g_yyzz_xxzz_1, g_yyzz_xyy_1, g_yyzz_xyyy_1, g_yyzz_xyyz_1, g_yyzz_xyz_1, g_yyzz_xyzz_1, g_yyzz_xzz_1, g_yyzz_xzzz_1, g_yyzz_yyy_1, g_yyzz_yyyy_1, g_yyzz_yyyz_1, g_yyzz_yyz_1, g_yyzz_yyzz_1, g_yyzz_yzz_1, g_yyzz_yzzz_1, g_yyzz_zzz_1, g_yyzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzz_xxxx_0[i] = 4.0 * g_yyzz_xxx_1[i] * fe_0 + g_yyzz_xxxx_1[i] * pa_x[i];

        g_xyyzz_xxxy_0[i] = 3.0 * g_yyzz_xxy_1[i] * fe_0 + g_yyzz_xxxy_1[i] * pa_x[i];

        g_xyyzz_xxxz_0[i] = 3.0 * g_yyzz_xxz_1[i] * fe_0 + g_yyzz_xxxz_1[i] * pa_x[i];

        g_xyyzz_xxyy_0[i] = 2.0 * g_yyzz_xyy_1[i] * fe_0 + g_yyzz_xxyy_1[i] * pa_x[i];

        g_xyyzz_xxyz_0[i] = 2.0 * g_yyzz_xyz_1[i] * fe_0 + g_yyzz_xxyz_1[i] * pa_x[i];

        g_xyyzz_xxzz_0[i] = 2.0 * g_yyzz_xzz_1[i] * fe_0 + g_yyzz_xxzz_1[i] * pa_x[i];

        g_xyyzz_xyyy_0[i] = g_yyzz_yyy_1[i] * fe_0 + g_yyzz_xyyy_1[i] * pa_x[i];

        g_xyyzz_xyyz_0[i] = g_yyzz_yyz_1[i] * fe_0 + g_yyzz_xyyz_1[i] * pa_x[i];

        g_xyyzz_xyzz_0[i] = g_yyzz_yzz_1[i] * fe_0 + g_yyzz_xyzz_1[i] * pa_x[i];

        g_xyyzz_xzzz_0[i] = g_yyzz_zzz_1[i] * fe_0 + g_yyzz_xzzz_1[i] * pa_x[i];

        g_xyyzz_yyyy_0[i] = g_yyzz_yyyy_1[i] * pa_x[i];

        g_xyyzz_yyyz_0[i] = g_yyzz_yyyz_1[i] * pa_x[i];

        g_xyyzz_yyzz_0[i] = g_yyzz_yyzz_1[i] * pa_x[i];

        g_xyyzz_yzzz_0[i] = g_yyzz_yzzz_1[i] * pa_x[i];

        g_xyyzz_zzzz_0[i] = g_yyzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 195-210 components of targeted buffer : HG

    auto g_xyzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 195);

    auto g_xyzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 196);

    auto g_xyzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 197);

    auto g_xyzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 198);

    auto g_xyzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 199);

    auto g_xyzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 200);

    auto g_xyzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 201);

    auto g_xyzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 202);

    auto g_xyzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 203);

    auto g_xyzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 204);

    auto g_xyzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 205);

    auto g_xyzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 206);

    auto g_xyzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 207);

    auto g_xyzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 208);

    auto g_xyzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 209);

    #pragma omp simd aligned(g_xyzzz_xxxx_0, g_xyzzz_xxxy_0, g_xyzzz_xxxz_0, g_xyzzz_xxyy_0, g_xyzzz_xxyz_0, g_xyzzz_xxzz_0, g_xyzzz_xyyy_0, g_xyzzz_xyyz_0, g_xyzzz_xyzz_0, g_xyzzz_xzzz_0, g_xyzzz_yyyy_0, g_xyzzz_yyyz_0, g_xyzzz_yyzz_0, g_xyzzz_yzzz_0, g_xyzzz_zzzz_0, g_xzzz_xxxx_1, g_xzzz_xxxz_1, g_xzzz_xxzz_1, g_xzzz_xzzz_1, g_yzzz_xxxy_1, g_yzzz_xxy_1, g_yzzz_xxyy_1, g_yzzz_xxyz_1, g_yzzz_xyy_1, g_yzzz_xyyy_1, g_yzzz_xyyz_1, g_yzzz_xyz_1, g_yzzz_xyzz_1, g_yzzz_yyy_1, g_yzzz_yyyy_1, g_yzzz_yyyz_1, g_yzzz_yyz_1, g_yzzz_yyzz_1, g_yzzz_yzz_1, g_yzzz_yzzz_1, g_yzzz_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzz_xxxx_0[i] = g_xzzz_xxxx_1[i] * pa_y[i];

        g_xyzzz_xxxy_0[i] = 3.0 * g_yzzz_xxy_1[i] * fe_0 + g_yzzz_xxxy_1[i] * pa_x[i];

        g_xyzzz_xxxz_0[i] = g_xzzz_xxxz_1[i] * pa_y[i];

        g_xyzzz_xxyy_0[i] = 2.0 * g_yzzz_xyy_1[i] * fe_0 + g_yzzz_xxyy_1[i] * pa_x[i];

        g_xyzzz_xxyz_0[i] = 2.0 * g_yzzz_xyz_1[i] * fe_0 + g_yzzz_xxyz_1[i] * pa_x[i];

        g_xyzzz_xxzz_0[i] = g_xzzz_xxzz_1[i] * pa_y[i];

        g_xyzzz_xyyy_0[i] = g_yzzz_yyy_1[i] * fe_0 + g_yzzz_xyyy_1[i] * pa_x[i];

        g_xyzzz_xyyz_0[i] = g_yzzz_yyz_1[i] * fe_0 + g_yzzz_xyyz_1[i] * pa_x[i];

        g_xyzzz_xyzz_0[i] = g_yzzz_yzz_1[i] * fe_0 + g_yzzz_xyzz_1[i] * pa_x[i];

        g_xyzzz_xzzz_0[i] = g_xzzz_xzzz_1[i] * pa_y[i];

        g_xyzzz_yyyy_0[i] = g_yzzz_yyyy_1[i] * pa_x[i];

        g_xyzzz_yyyz_0[i] = g_yzzz_yyyz_1[i] * pa_x[i];

        g_xyzzz_yyzz_0[i] = g_yzzz_yyzz_1[i] * pa_x[i];

        g_xyzzz_yzzz_0[i] = g_yzzz_yzzz_1[i] * pa_x[i];

        g_xyzzz_zzzz_0[i] = g_yzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 210-225 components of targeted buffer : HG

    auto g_xzzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 210);

    auto g_xzzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 211);

    auto g_xzzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 212);

    auto g_xzzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 213);

    auto g_xzzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 214);

    auto g_xzzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 215);

    auto g_xzzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 216);

    auto g_xzzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 217);

    auto g_xzzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 218);

    auto g_xzzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 219);

    auto g_xzzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 220);

    auto g_xzzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 221);

    auto g_xzzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 222);

    auto g_xzzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 223);

    auto g_xzzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 224);

    #pragma omp simd aligned(g_xzzzz_xxxx_0, g_xzzzz_xxxy_0, g_xzzzz_xxxz_0, g_xzzzz_xxyy_0, g_xzzzz_xxyz_0, g_xzzzz_xxzz_0, g_xzzzz_xyyy_0, g_xzzzz_xyyz_0, g_xzzzz_xyzz_0, g_xzzzz_xzzz_0, g_xzzzz_yyyy_0, g_xzzzz_yyyz_0, g_xzzzz_yyzz_0, g_xzzzz_yzzz_0, g_xzzzz_zzzz_0, g_zzzz_xxx_1, g_zzzz_xxxx_1, g_zzzz_xxxy_1, g_zzzz_xxxz_1, g_zzzz_xxy_1, g_zzzz_xxyy_1, g_zzzz_xxyz_1, g_zzzz_xxz_1, g_zzzz_xxzz_1, g_zzzz_xyy_1, g_zzzz_xyyy_1, g_zzzz_xyyz_1, g_zzzz_xyz_1, g_zzzz_xyzz_1, g_zzzz_xzz_1, g_zzzz_xzzz_1, g_zzzz_yyy_1, g_zzzz_yyyy_1, g_zzzz_yyyz_1, g_zzzz_yyz_1, g_zzzz_yyzz_1, g_zzzz_yzz_1, g_zzzz_yzzz_1, g_zzzz_zzz_1, g_zzzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzz_xxxx_0[i] = 4.0 * g_zzzz_xxx_1[i] * fe_0 + g_zzzz_xxxx_1[i] * pa_x[i];

        g_xzzzz_xxxy_0[i] = 3.0 * g_zzzz_xxy_1[i] * fe_0 + g_zzzz_xxxy_1[i] * pa_x[i];

        g_xzzzz_xxxz_0[i] = 3.0 * g_zzzz_xxz_1[i] * fe_0 + g_zzzz_xxxz_1[i] * pa_x[i];

        g_xzzzz_xxyy_0[i] = 2.0 * g_zzzz_xyy_1[i] * fe_0 + g_zzzz_xxyy_1[i] * pa_x[i];

        g_xzzzz_xxyz_0[i] = 2.0 * g_zzzz_xyz_1[i] * fe_0 + g_zzzz_xxyz_1[i] * pa_x[i];

        g_xzzzz_xxzz_0[i] = 2.0 * g_zzzz_xzz_1[i] * fe_0 + g_zzzz_xxzz_1[i] * pa_x[i];

        g_xzzzz_xyyy_0[i] = g_zzzz_yyy_1[i] * fe_0 + g_zzzz_xyyy_1[i] * pa_x[i];

        g_xzzzz_xyyz_0[i] = g_zzzz_yyz_1[i] * fe_0 + g_zzzz_xyyz_1[i] * pa_x[i];

        g_xzzzz_xyzz_0[i] = g_zzzz_yzz_1[i] * fe_0 + g_zzzz_xyzz_1[i] * pa_x[i];

        g_xzzzz_xzzz_0[i] = g_zzzz_zzz_1[i] * fe_0 + g_zzzz_xzzz_1[i] * pa_x[i];

        g_xzzzz_yyyy_0[i] = g_zzzz_yyyy_1[i] * pa_x[i];

        g_xzzzz_yyyz_0[i] = g_zzzz_yyyz_1[i] * pa_x[i];

        g_xzzzz_yyzz_0[i] = g_zzzz_yyzz_1[i] * pa_x[i];

        g_xzzzz_yzzz_0[i] = g_zzzz_yzzz_1[i] * pa_x[i];

        g_xzzzz_zzzz_0[i] = g_zzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : HG

    auto g_yyyyy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 225);

    auto g_yyyyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 226);

    auto g_yyyyy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 227);

    auto g_yyyyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 228);

    auto g_yyyyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 229);

    auto g_yyyyy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 230);

    auto g_yyyyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 231);

    auto g_yyyyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 232);

    auto g_yyyyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 233);

    auto g_yyyyy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 234);

    auto g_yyyyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 235);

    auto g_yyyyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 236);

    auto g_yyyyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 237);

    auto g_yyyyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 238);

    auto g_yyyyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 239);

    #pragma omp simd aligned(g_yyy_xxxx_0, g_yyy_xxxx_1, g_yyy_xxxy_0, g_yyy_xxxy_1, g_yyy_xxxz_0, g_yyy_xxxz_1, g_yyy_xxyy_0, g_yyy_xxyy_1, g_yyy_xxyz_0, g_yyy_xxyz_1, g_yyy_xxzz_0, g_yyy_xxzz_1, g_yyy_xyyy_0, g_yyy_xyyy_1, g_yyy_xyyz_0, g_yyy_xyyz_1, g_yyy_xyzz_0, g_yyy_xyzz_1, g_yyy_xzzz_0, g_yyy_xzzz_1, g_yyy_yyyy_0, g_yyy_yyyy_1, g_yyy_yyyz_0, g_yyy_yyyz_1, g_yyy_yyzz_0, g_yyy_yyzz_1, g_yyy_yzzz_0, g_yyy_yzzz_1, g_yyy_zzzz_0, g_yyy_zzzz_1, g_yyyy_xxx_1, g_yyyy_xxxx_1, g_yyyy_xxxy_1, g_yyyy_xxxz_1, g_yyyy_xxy_1, g_yyyy_xxyy_1, g_yyyy_xxyz_1, g_yyyy_xxz_1, g_yyyy_xxzz_1, g_yyyy_xyy_1, g_yyyy_xyyy_1, g_yyyy_xyyz_1, g_yyyy_xyz_1, g_yyyy_xyzz_1, g_yyyy_xzz_1, g_yyyy_xzzz_1, g_yyyy_yyy_1, g_yyyy_yyyy_1, g_yyyy_yyyz_1, g_yyyy_yyz_1, g_yyyy_yyzz_1, g_yyyy_yzz_1, g_yyyy_yzzz_1, g_yyyy_zzz_1, g_yyyy_zzzz_1, g_yyyyy_xxxx_0, g_yyyyy_xxxy_0, g_yyyyy_xxxz_0, g_yyyyy_xxyy_0, g_yyyyy_xxyz_0, g_yyyyy_xxzz_0, g_yyyyy_xyyy_0, g_yyyyy_xyyz_0, g_yyyyy_xyzz_0, g_yyyyy_xzzz_0, g_yyyyy_yyyy_0, g_yyyyy_yyyz_0, g_yyyyy_yyzz_0, g_yyyyy_yzzz_0, g_yyyyy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyy_xxxx_0[i] = 4.0 * g_yyy_xxxx_0[i] * fbe_0 - 4.0 * g_yyy_xxxx_1[i] * fz_be_0 + g_yyyy_xxxx_1[i] * pa_y[i];

        g_yyyyy_xxxy_0[i] = 4.0 * g_yyy_xxxy_0[i] * fbe_0 - 4.0 * g_yyy_xxxy_1[i] * fz_be_0 + g_yyyy_xxx_1[i] * fe_0 + g_yyyy_xxxy_1[i] * pa_y[i];

        g_yyyyy_xxxz_0[i] = 4.0 * g_yyy_xxxz_0[i] * fbe_0 - 4.0 * g_yyy_xxxz_1[i] * fz_be_0 + g_yyyy_xxxz_1[i] * pa_y[i];

        g_yyyyy_xxyy_0[i] = 4.0 * g_yyy_xxyy_0[i] * fbe_0 - 4.0 * g_yyy_xxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_xxy_1[i] * fe_0 + g_yyyy_xxyy_1[i] * pa_y[i];

        g_yyyyy_xxyz_0[i] = 4.0 * g_yyy_xxyz_0[i] * fbe_0 - 4.0 * g_yyy_xxyz_1[i] * fz_be_0 + g_yyyy_xxz_1[i] * fe_0 + g_yyyy_xxyz_1[i] * pa_y[i];

        g_yyyyy_xxzz_0[i] = 4.0 * g_yyy_xxzz_0[i] * fbe_0 - 4.0 * g_yyy_xxzz_1[i] * fz_be_0 + g_yyyy_xxzz_1[i] * pa_y[i];

        g_yyyyy_xyyy_0[i] = 4.0 * g_yyy_xyyy_0[i] * fbe_0 - 4.0 * g_yyy_xyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_xyy_1[i] * fe_0 + g_yyyy_xyyy_1[i] * pa_y[i];

        g_yyyyy_xyyz_0[i] = 4.0 * g_yyy_xyyz_0[i] * fbe_0 - 4.0 * g_yyy_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_xyz_1[i] * fe_0 + g_yyyy_xyyz_1[i] * pa_y[i];

        g_yyyyy_xyzz_0[i] = 4.0 * g_yyy_xyzz_0[i] * fbe_0 - 4.0 * g_yyy_xyzz_1[i] * fz_be_0 + g_yyyy_xzz_1[i] * fe_0 + g_yyyy_xyzz_1[i] * pa_y[i];

        g_yyyyy_xzzz_0[i] = 4.0 * g_yyy_xzzz_0[i] * fbe_0 - 4.0 * g_yyy_xzzz_1[i] * fz_be_0 + g_yyyy_xzzz_1[i] * pa_y[i];

        g_yyyyy_yyyy_0[i] = 4.0 * g_yyy_yyyy_0[i] * fbe_0 - 4.0 * g_yyy_yyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_yyy_1[i] * fe_0 + g_yyyy_yyyy_1[i] * pa_y[i];

        g_yyyyy_yyyz_0[i] = 4.0 * g_yyy_yyyz_0[i] * fbe_0 - 4.0 * g_yyy_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_yyz_1[i] * fe_0 + g_yyyy_yyyz_1[i] * pa_y[i];

        g_yyyyy_yyzz_0[i] = 4.0 * g_yyy_yyzz_0[i] * fbe_0 - 4.0 * g_yyy_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_yzz_1[i] * fe_0 + g_yyyy_yyzz_1[i] * pa_y[i];

        g_yyyyy_yzzz_0[i] = 4.0 * g_yyy_yzzz_0[i] * fbe_0 - 4.0 * g_yyy_yzzz_1[i] * fz_be_0 + g_yyyy_zzz_1[i] * fe_0 + g_yyyy_yzzz_1[i] * pa_y[i];

        g_yyyyy_zzzz_0[i] = 4.0 * g_yyy_zzzz_0[i] * fbe_0 - 4.0 * g_yyy_zzzz_1[i] * fz_be_0 + g_yyyy_zzzz_1[i] * pa_y[i];
    }

    // Set up 240-255 components of targeted buffer : HG

    auto g_yyyyz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 240);

    auto g_yyyyz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 241);

    auto g_yyyyz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 242);

    auto g_yyyyz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 243);

    auto g_yyyyz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 244);

    auto g_yyyyz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 245);

    auto g_yyyyz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 246);

    auto g_yyyyz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 247);

    auto g_yyyyz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 248);

    auto g_yyyyz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 249);

    auto g_yyyyz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 250);

    auto g_yyyyz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 251);

    auto g_yyyyz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 252);

    auto g_yyyyz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 253);

    auto g_yyyyz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 254);

    #pragma omp simd aligned(g_yyyy_xxx_1, g_yyyy_xxxx_1, g_yyyy_xxxy_1, g_yyyy_xxxz_1, g_yyyy_xxy_1, g_yyyy_xxyy_1, g_yyyy_xxyz_1, g_yyyy_xxz_1, g_yyyy_xxzz_1, g_yyyy_xyy_1, g_yyyy_xyyy_1, g_yyyy_xyyz_1, g_yyyy_xyz_1, g_yyyy_xyzz_1, g_yyyy_xzz_1, g_yyyy_xzzz_1, g_yyyy_yyy_1, g_yyyy_yyyy_1, g_yyyy_yyyz_1, g_yyyy_yyz_1, g_yyyy_yyzz_1, g_yyyy_yzz_1, g_yyyy_yzzz_1, g_yyyy_zzz_1, g_yyyy_zzzz_1, g_yyyyz_xxxx_0, g_yyyyz_xxxy_0, g_yyyyz_xxxz_0, g_yyyyz_xxyy_0, g_yyyyz_xxyz_0, g_yyyyz_xxzz_0, g_yyyyz_xyyy_0, g_yyyyz_xyyz_0, g_yyyyz_xyzz_0, g_yyyyz_xzzz_0, g_yyyyz_yyyy_0, g_yyyyz_yyyz_0, g_yyyyz_yyzz_0, g_yyyyz_yzzz_0, g_yyyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyz_xxxx_0[i] = g_yyyy_xxxx_1[i] * pa_z[i];

        g_yyyyz_xxxy_0[i] = g_yyyy_xxxy_1[i] * pa_z[i];

        g_yyyyz_xxxz_0[i] = g_yyyy_xxx_1[i] * fe_0 + g_yyyy_xxxz_1[i] * pa_z[i];

        g_yyyyz_xxyy_0[i] = g_yyyy_xxyy_1[i] * pa_z[i];

        g_yyyyz_xxyz_0[i] = g_yyyy_xxy_1[i] * fe_0 + g_yyyy_xxyz_1[i] * pa_z[i];

        g_yyyyz_xxzz_0[i] = 2.0 * g_yyyy_xxz_1[i] * fe_0 + g_yyyy_xxzz_1[i] * pa_z[i];

        g_yyyyz_xyyy_0[i] = g_yyyy_xyyy_1[i] * pa_z[i];

        g_yyyyz_xyyz_0[i] = g_yyyy_xyy_1[i] * fe_0 + g_yyyy_xyyz_1[i] * pa_z[i];

        g_yyyyz_xyzz_0[i] = 2.0 * g_yyyy_xyz_1[i] * fe_0 + g_yyyy_xyzz_1[i] * pa_z[i];

        g_yyyyz_xzzz_0[i] = 3.0 * g_yyyy_xzz_1[i] * fe_0 + g_yyyy_xzzz_1[i] * pa_z[i];

        g_yyyyz_yyyy_0[i] = g_yyyy_yyyy_1[i] * pa_z[i];

        g_yyyyz_yyyz_0[i] = g_yyyy_yyy_1[i] * fe_0 + g_yyyy_yyyz_1[i] * pa_z[i];

        g_yyyyz_yyzz_0[i] = 2.0 * g_yyyy_yyz_1[i] * fe_0 + g_yyyy_yyzz_1[i] * pa_z[i];

        g_yyyyz_yzzz_0[i] = 3.0 * g_yyyy_yzz_1[i] * fe_0 + g_yyyy_yzzz_1[i] * pa_z[i];

        g_yyyyz_zzzz_0[i] = 4.0 * g_yyyy_zzz_1[i] * fe_0 + g_yyyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 255-270 components of targeted buffer : HG

    auto g_yyyzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 255);

    auto g_yyyzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 256);

    auto g_yyyzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 257);

    auto g_yyyzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 258);

    auto g_yyyzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 259);

    auto g_yyyzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 260);

    auto g_yyyzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 261);

    auto g_yyyzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 262);

    auto g_yyyzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 263);

    auto g_yyyzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 264);

    auto g_yyyzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 265);

    auto g_yyyzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 266);

    auto g_yyyzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 267);

    auto g_yyyzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 268);

    auto g_yyyzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 269);

    #pragma omp simd aligned(g_yyy_xxxy_0, g_yyy_xxxy_1, g_yyy_xxyy_0, g_yyy_xxyy_1, g_yyy_xyyy_0, g_yyy_xyyy_1, g_yyy_yyyy_0, g_yyy_yyyy_1, g_yyyz_xxxy_1, g_yyyz_xxyy_1, g_yyyz_xyyy_1, g_yyyz_yyyy_1, g_yyyzz_xxxx_0, g_yyyzz_xxxy_0, g_yyyzz_xxxz_0, g_yyyzz_xxyy_0, g_yyyzz_xxyz_0, g_yyyzz_xxzz_0, g_yyyzz_xyyy_0, g_yyyzz_xyyz_0, g_yyyzz_xyzz_0, g_yyyzz_xzzz_0, g_yyyzz_yyyy_0, g_yyyzz_yyyz_0, g_yyyzz_yyzz_0, g_yyyzz_yzzz_0, g_yyyzz_zzzz_0, g_yyzz_xxxx_1, g_yyzz_xxxz_1, g_yyzz_xxyz_1, g_yyzz_xxz_1, g_yyzz_xxzz_1, g_yyzz_xyyz_1, g_yyzz_xyz_1, g_yyzz_xyzz_1, g_yyzz_xzz_1, g_yyzz_xzzz_1, g_yyzz_yyyz_1, g_yyzz_yyz_1, g_yyzz_yyzz_1, g_yyzz_yzz_1, g_yyzz_yzzz_1, g_yyzz_zzz_1, g_yyzz_zzzz_1, g_yzz_xxxx_0, g_yzz_xxxx_1, g_yzz_xxxz_0, g_yzz_xxxz_1, g_yzz_xxyz_0, g_yzz_xxyz_1, g_yzz_xxzz_0, g_yzz_xxzz_1, g_yzz_xyyz_0, g_yzz_xyyz_1, g_yzz_xyzz_0, g_yzz_xyzz_1, g_yzz_xzzz_0, g_yzz_xzzz_1, g_yzz_yyyz_0, g_yzz_yyyz_1, g_yzz_yyzz_0, g_yzz_yyzz_1, g_yzz_yzzz_0, g_yzz_yzzz_1, g_yzz_zzzz_0, g_yzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzz_xxxx_0[i] = 2.0 * g_yzz_xxxx_0[i] * fbe_0 - 2.0 * g_yzz_xxxx_1[i] * fz_be_0 + g_yyzz_xxxx_1[i] * pa_y[i];

        g_yyyzz_xxxy_0[i] = g_yyy_xxxy_0[i] * fbe_0 - g_yyy_xxxy_1[i] * fz_be_0 + g_yyyz_xxxy_1[i] * pa_z[i];

        g_yyyzz_xxxz_0[i] = 2.0 * g_yzz_xxxz_0[i] * fbe_0 - 2.0 * g_yzz_xxxz_1[i] * fz_be_0 + g_yyzz_xxxz_1[i] * pa_y[i];

        g_yyyzz_xxyy_0[i] = g_yyy_xxyy_0[i] * fbe_0 - g_yyy_xxyy_1[i] * fz_be_0 + g_yyyz_xxyy_1[i] * pa_z[i];

        g_yyyzz_xxyz_0[i] = 2.0 * g_yzz_xxyz_0[i] * fbe_0 - 2.0 * g_yzz_xxyz_1[i] * fz_be_0 + g_yyzz_xxz_1[i] * fe_0 + g_yyzz_xxyz_1[i] * pa_y[i];

        g_yyyzz_xxzz_0[i] = 2.0 * g_yzz_xxzz_0[i] * fbe_0 - 2.0 * g_yzz_xxzz_1[i] * fz_be_0 + g_yyzz_xxzz_1[i] * pa_y[i];

        g_yyyzz_xyyy_0[i] = g_yyy_xyyy_0[i] * fbe_0 - g_yyy_xyyy_1[i] * fz_be_0 + g_yyyz_xyyy_1[i] * pa_z[i];

        g_yyyzz_xyyz_0[i] = 2.0 * g_yzz_xyyz_0[i] * fbe_0 - 2.0 * g_yzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_xyz_1[i] * fe_0 + g_yyzz_xyyz_1[i] * pa_y[i];

        g_yyyzz_xyzz_0[i] = 2.0 * g_yzz_xyzz_0[i] * fbe_0 - 2.0 * g_yzz_xyzz_1[i] * fz_be_0 + g_yyzz_xzz_1[i] * fe_0 + g_yyzz_xyzz_1[i] * pa_y[i];

        g_yyyzz_xzzz_0[i] = 2.0 * g_yzz_xzzz_0[i] * fbe_0 - 2.0 * g_yzz_xzzz_1[i] * fz_be_0 + g_yyzz_xzzz_1[i] * pa_y[i];

        g_yyyzz_yyyy_0[i] = g_yyy_yyyy_0[i] * fbe_0 - g_yyy_yyyy_1[i] * fz_be_0 + g_yyyz_yyyy_1[i] * pa_z[i];

        g_yyyzz_yyyz_0[i] = 2.0 * g_yzz_yyyz_0[i] * fbe_0 - 2.0 * g_yzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_yyz_1[i] * fe_0 + g_yyzz_yyyz_1[i] * pa_y[i];

        g_yyyzz_yyzz_0[i] = 2.0 * g_yzz_yyzz_0[i] * fbe_0 - 2.0 * g_yzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_yzz_1[i] * fe_0 + g_yyzz_yyzz_1[i] * pa_y[i];

        g_yyyzz_yzzz_0[i] = 2.0 * g_yzz_yzzz_0[i] * fbe_0 - 2.0 * g_yzz_yzzz_1[i] * fz_be_0 + g_yyzz_zzz_1[i] * fe_0 + g_yyzz_yzzz_1[i] * pa_y[i];

        g_yyyzz_zzzz_0[i] = 2.0 * g_yzz_zzzz_0[i] * fbe_0 - 2.0 * g_yzz_zzzz_1[i] * fz_be_0 + g_yyzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 270-285 components of targeted buffer : HG

    auto g_yyzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 270);

    auto g_yyzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 271);

    auto g_yyzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 272);

    auto g_yyzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 273);

    auto g_yyzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 274);

    auto g_yyzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 275);

    auto g_yyzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 276);

    auto g_yyzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 277);

    auto g_yyzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 278);

    auto g_yyzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 279);

    auto g_yyzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 280);

    auto g_yyzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 281);

    auto g_yyzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 282);

    auto g_yyzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 283);

    auto g_yyzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 284);

    #pragma omp simd aligned(g_yyz_xxxy_0, g_yyz_xxxy_1, g_yyz_xxyy_0, g_yyz_xxyy_1, g_yyz_xyyy_0, g_yyz_xyyy_1, g_yyz_yyyy_0, g_yyz_yyyy_1, g_yyzz_xxxy_1, g_yyzz_xxyy_1, g_yyzz_xyyy_1, g_yyzz_yyyy_1, g_yyzzz_xxxx_0, g_yyzzz_xxxy_0, g_yyzzz_xxxz_0, g_yyzzz_xxyy_0, g_yyzzz_xxyz_0, g_yyzzz_xxzz_0, g_yyzzz_xyyy_0, g_yyzzz_xyyz_0, g_yyzzz_xyzz_0, g_yyzzz_xzzz_0, g_yyzzz_yyyy_0, g_yyzzz_yyyz_0, g_yyzzz_yyzz_0, g_yyzzz_yzzz_0, g_yyzzz_zzzz_0, g_yzzz_xxxx_1, g_yzzz_xxxz_1, g_yzzz_xxyz_1, g_yzzz_xxz_1, g_yzzz_xxzz_1, g_yzzz_xyyz_1, g_yzzz_xyz_1, g_yzzz_xyzz_1, g_yzzz_xzz_1, g_yzzz_xzzz_1, g_yzzz_yyyz_1, g_yzzz_yyz_1, g_yzzz_yyzz_1, g_yzzz_yzz_1, g_yzzz_yzzz_1, g_yzzz_zzz_1, g_yzzz_zzzz_1, g_zzz_xxxx_0, g_zzz_xxxx_1, g_zzz_xxxz_0, g_zzz_xxxz_1, g_zzz_xxyz_0, g_zzz_xxyz_1, g_zzz_xxzz_0, g_zzz_xxzz_1, g_zzz_xyyz_0, g_zzz_xyyz_1, g_zzz_xyzz_0, g_zzz_xyzz_1, g_zzz_xzzz_0, g_zzz_xzzz_1, g_zzz_yyyz_0, g_zzz_yyyz_1, g_zzz_yyzz_0, g_zzz_yyzz_1, g_zzz_yzzz_0, g_zzz_yzzz_1, g_zzz_zzzz_0, g_zzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzz_xxxx_0[i] = g_zzz_xxxx_0[i] * fbe_0 - g_zzz_xxxx_1[i] * fz_be_0 + g_yzzz_xxxx_1[i] * pa_y[i];

        g_yyzzz_xxxy_0[i] = 2.0 * g_yyz_xxxy_0[i] * fbe_0 - 2.0 * g_yyz_xxxy_1[i] * fz_be_0 + g_yyzz_xxxy_1[i] * pa_z[i];

        g_yyzzz_xxxz_0[i] = g_zzz_xxxz_0[i] * fbe_0 - g_zzz_xxxz_1[i] * fz_be_0 + g_yzzz_xxxz_1[i] * pa_y[i];

        g_yyzzz_xxyy_0[i] = 2.0 * g_yyz_xxyy_0[i] * fbe_0 - 2.0 * g_yyz_xxyy_1[i] * fz_be_0 + g_yyzz_xxyy_1[i] * pa_z[i];

        g_yyzzz_xxyz_0[i] = g_zzz_xxyz_0[i] * fbe_0 - g_zzz_xxyz_1[i] * fz_be_0 + g_yzzz_xxz_1[i] * fe_0 + g_yzzz_xxyz_1[i] * pa_y[i];

        g_yyzzz_xxzz_0[i] = g_zzz_xxzz_0[i] * fbe_0 - g_zzz_xxzz_1[i] * fz_be_0 + g_yzzz_xxzz_1[i] * pa_y[i];

        g_yyzzz_xyyy_0[i] = 2.0 * g_yyz_xyyy_0[i] * fbe_0 - 2.0 * g_yyz_xyyy_1[i] * fz_be_0 + g_yyzz_xyyy_1[i] * pa_z[i];

        g_yyzzz_xyyz_0[i] = g_zzz_xyyz_0[i] * fbe_0 - g_zzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_xyz_1[i] * fe_0 + g_yzzz_xyyz_1[i] * pa_y[i];

        g_yyzzz_xyzz_0[i] = g_zzz_xyzz_0[i] * fbe_0 - g_zzz_xyzz_1[i] * fz_be_0 + g_yzzz_xzz_1[i] * fe_0 + g_yzzz_xyzz_1[i] * pa_y[i];

        g_yyzzz_xzzz_0[i] = g_zzz_xzzz_0[i] * fbe_0 - g_zzz_xzzz_1[i] * fz_be_0 + g_yzzz_xzzz_1[i] * pa_y[i];

        g_yyzzz_yyyy_0[i] = 2.0 * g_yyz_yyyy_0[i] * fbe_0 - 2.0 * g_yyz_yyyy_1[i] * fz_be_0 + g_yyzz_yyyy_1[i] * pa_z[i];

        g_yyzzz_yyyz_0[i] = g_zzz_yyyz_0[i] * fbe_0 - g_zzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_yyz_1[i] * fe_0 + g_yzzz_yyyz_1[i] * pa_y[i];

        g_yyzzz_yyzz_0[i] = g_zzz_yyzz_0[i] * fbe_0 - g_zzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_yzz_1[i] * fe_0 + g_yzzz_yyzz_1[i] * pa_y[i];

        g_yyzzz_yzzz_0[i] = g_zzz_yzzz_0[i] * fbe_0 - g_zzz_yzzz_1[i] * fz_be_0 + g_yzzz_zzz_1[i] * fe_0 + g_yzzz_yzzz_1[i] * pa_y[i];

        g_yyzzz_zzzz_0[i] = g_zzz_zzzz_0[i] * fbe_0 - g_zzz_zzzz_1[i] * fz_be_0 + g_yzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 285-300 components of targeted buffer : HG

    auto g_yzzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 285);

    auto g_yzzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 286);

    auto g_yzzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 287);

    auto g_yzzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 288);

    auto g_yzzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 289);

    auto g_yzzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 290);

    auto g_yzzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 291);

    auto g_yzzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 292);

    auto g_yzzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 293);

    auto g_yzzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 294);

    auto g_yzzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 295);

    auto g_yzzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 296);

    auto g_yzzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 297);

    auto g_yzzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 298);

    auto g_yzzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 299);

    #pragma omp simd aligned(g_yzzzz_xxxx_0, g_yzzzz_xxxy_0, g_yzzzz_xxxz_0, g_yzzzz_xxyy_0, g_yzzzz_xxyz_0, g_yzzzz_xxzz_0, g_yzzzz_xyyy_0, g_yzzzz_xyyz_0, g_yzzzz_xyzz_0, g_yzzzz_xzzz_0, g_yzzzz_yyyy_0, g_yzzzz_yyyz_0, g_yzzzz_yyzz_0, g_yzzzz_yzzz_0, g_yzzzz_zzzz_0, g_zzzz_xxx_1, g_zzzz_xxxx_1, g_zzzz_xxxy_1, g_zzzz_xxxz_1, g_zzzz_xxy_1, g_zzzz_xxyy_1, g_zzzz_xxyz_1, g_zzzz_xxz_1, g_zzzz_xxzz_1, g_zzzz_xyy_1, g_zzzz_xyyy_1, g_zzzz_xyyz_1, g_zzzz_xyz_1, g_zzzz_xyzz_1, g_zzzz_xzz_1, g_zzzz_xzzz_1, g_zzzz_yyy_1, g_zzzz_yyyy_1, g_zzzz_yyyz_1, g_zzzz_yyz_1, g_zzzz_yyzz_1, g_zzzz_yzz_1, g_zzzz_yzzz_1, g_zzzz_zzz_1, g_zzzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzz_xxxx_0[i] = g_zzzz_xxxx_1[i] * pa_y[i];

        g_yzzzz_xxxy_0[i] = g_zzzz_xxx_1[i] * fe_0 + g_zzzz_xxxy_1[i] * pa_y[i];

        g_yzzzz_xxxz_0[i] = g_zzzz_xxxz_1[i] * pa_y[i];

        g_yzzzz_xxyy_0[i] = 2.0 * g_zzzz_xxy_1[i] * fe_0 + g_zzzz_xxyy_1[i] * pa_y[i];

        g_yzzzz_xxyz_0[i] = g_zzzz_xxz_1[i] * fe_0 + g_zzzz_xxyz_1[i] * pa_y[i];

        g_yzzzz_xxzz_0[i] = g_zzzz_xxzz_1[i] * pa_y[i];

        g_yzzzz_xyyy_0[i] = 3.0 * g_zzzz_xyy_1[i] * fe_0 + g_zzzz_xyyy_1[i] * pa_y[i];

        g_yzzzz_xyyz_0[i] = 2.0 * g_zzzz_xyz_1[i] * fe_0 + g_zzzz_xyyz_1[i] * pa_y[i];

        g_yzzzz_xyzz_0[i] = g_zzzz_xzz_1[i] * fe_0 + g_zzzz_xyzz_1[i] * pa_y[i];

        g_yzzzz_xzzz_0[i] = g_zzzz_xzzz_1[i] * pa_y[i];

        g_yzzzz_yyyy_0[i] = 4.0 * g_zzzz_yyy_1[i] * fe_0 + g_zzzz_yyyy_1[i] * pa_y[i];

        g_yzzzz_yyyz_0[i] = 3.0 * g_zzzz_yyz_1[i] * fe_0 + g_zzzz_yyyz_1[i] * pa_y[i];

        g_yzzzz_yyzz_0[i] = 2.0 * g_zzzz_yzz_1[i] * fe_0 + g_zzzz_yyzz_1[i] * pa_y[i];

        g_yzzzz_yzzz_0[i] = g_zzzz_zzz_1[i] * fe_0 + g_zzzz_yzzz_1[i] * pa_y[i];

        g_yzzzz_zzzz_0[i] = g_zzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 300-315 components of targeted buffer : HG

    auto g_zzzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 300);

    auto g_zzzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 301);

    auto g_zzzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 302);

    auto g_zzzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 303);

    auto g_zzzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 304);

    auto g_zzzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 305);

    auto g_zzzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 306);

    auto g_zzzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 307);

    auto g_zzzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 308);

    auto g_zzzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 309);

    auto g_zzzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 310);

    auto g_zzzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 311);

    auto g_zzzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 312);

    auto g_zzzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 313);

    auto g_zzzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 314);

    #pragma omp simd aligned(g_zzz_xxxx_0, g_zzz_xxxx_1, g_zzz_xxxy_0, g_zzz_xxxy_1, g_zzz_xxxz_0, g_zzz_xxxz_1, g_zzz_xxyy_0, g_zzz_xxyy_1, g_zzz_xxyz_0, g_zzz_xxyz_1, g_zzz_xxzz_0, g_zzz_xxzz_1, g_zzz_xyyy_0, g_zzz_xyyy_1, g_zzz_xyyz_0, g_zzz_xyyz_1, g_zzz_xyzz_0, g_zzz_xyzz_1, g_zzz_xzzz_0, g_zzz_xzzz_1, g_zzz_yyyy_0, g_zzz_yyyy_1, g_zzz_yyyz_0, g_zzz_yyyz_1, g_zzz_yyzz_0, g_zzz_yyzz_1, g_zzz_yzzz_0, g_zzz_yzzz_1, g_zzz_zzzz_0, g_zzz_zzzz_1, g_zzzz_xxx_1, g_zzzz_xxxx_1, g_zzzz_xxxy_1, g_zzzz_xxxz_1, g_zzzz_xxy_1, g_zzzz_xxyy_1, g_zzzz_xxyz_1, g_zzzz_xxz_1, g_zzzz_xxzz_1, g_zzzz_xyy_1, g_zzzz_xyyy_1, g_zzzz_xyyz_1, g_zzzz_xyz_1, g_zzzz_xyzz_1, g_zzzz_xzz_1, g_zzzz_xzzz_1, g_zzzz_yyy_1, g_zzzz_yyyy_1, g_zzzz_yyyz_1, g_zzzz_yyz_1, g_zzzz_yyzz_1, g_zzzz_yzz_1, g_zzzz_yzzz_1, g_zzzz_zzz_1, g_zzzz_zzzz_1, g_zzzzz_xxxx_0, g_zzzzz_xxxy_0, g_zzzzz_xxxz_0, g_zzzzz_xxyy_0, g_zzzzz_xxyz_0, g_zzzzz_xxzz_0, g_zzzzz_xyyy_0, g_zzzzz_xyyz_0, g_zzzzz_xyzz_0, g_zzzzz_xzzz_0, g_zzzzz_yyyy_0, g_zzzzz_yyyz_0, g_zzzzz_yyzz_0, g_zzzzz_yzzz_0, g_zzzzz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzz_xxxx_0[i] = 4.0 * g_zzz_xxxx_0[i] * fbe_0 - 4.0 * g_zzz_xxxx_1[i] * fz_be_0 + g_zzzz_xxxx_1[i] * pa_z[i];

        g_zzzzz_xxxy_0[i] = 4.0 * g_zzz_xxxy_0[i] * fbe_0 - 4.0 * g_zzz_xxxy_1[i] * fz_be_0 + g_zzzz_xxxy_1[i] * pa_z[i];

        g_zzzzz_xxxz_0[i] = 4.0 * g_zzz_xxxz_0[i] * fbe_0 - 4.0 * g_zzz_xxxz_1[i] * fz_be_0 + g_zzzz_xxx_1[i] * fe_0 + g_zzzz_xxxz_1[i] * pa_z[i];

        g_zzzzz_xxyy_0[i] = 4.0 * g_zzz_xxyy_0[i] * fbe_0 - 4.0 * g_zzz_xxyy_1[i] * fz_be_0 + g_zzzz_xxyy_1[i] * pa_z[i];

        g_zzzzz_xxyz_0[i] = 4.0 * g_zzz_xxyz_0[i] * fbe_0 - 4.0 * g_zzz_xxyz_1[i] * fz_be_0 + g_zzzz_xxy_1[i] * fe_0 + g_zzzz_xxyz_1[i] * pa_z[i];

        g_zzzzz_xxzz_0[i] = 4.0 * g_zzz_xxzz_0[i] * fbe_0 - 4.0 * g_zzz_xxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xxz_1[i] * fe_0 + g_zzzz_xxzz_1[i] * pa_z[i];

        g_zzzzz_xyyy_0[i] = 4.0 * g_zzz_xyyy_0[i] * fbe_0 - 4.0 * g_zzz_xyyy_1[i] * fz_be_0 + g_zzzz_xyyy_1[i] * pa_z[i];

        g_zzzzz_xyyz_0[i] = 4.0 * g_zzz_xyyz_0[i] * fbe_0 - 4.0 * g_zzz_xyyz_1[i] * fz_be_0 + g_zzzz_xyy_1[i] * fe_0 + g_zzzz_xyyz_1[i] * pa_z[i];

        g_zzzzz_xyzz_0[i] = 4.0 * g_zzz_xyzz_0[i] * fbe_0 - 4.0 * g_zzz_xyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xyz_1[i] * fe_0 + g_zzzz_xyzz_1[i] * pa_z[i];

        g_zzzzz_xzzz_0[i] = 4.0 * g_zzz_xzzz_0[i] * fbe_0 - 4.0 * g_zzz_xzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_xzz_1[i] * fe_0 + g_zzzz_xzzz_1[i] * pa_z[i];

        g_zzzzz_yyyy_0[i] = 4.0 * g_zzz_yyyy_0[i] * fbe_0 - 4.0 * g_zzz_yyyy_1[i] * fz_be_0 + g_zzzz_yyyy_1[i] * pa_z[i];

        g_zzzzz_yyyz_0[i] = 4.0 * g_zzz_yyyz_0[i] * fbe_0 - 4.0 * g_zzz_yyyz_1[i] * fz_be_0 + g_zzzz_yyy_1[i] * fe_0 + g_zzzz_yyyz_1[i] * pa_z[i];

        g_zzzzz_yyzz_0[i] = 4.0 * g_zzz_yyzz_0[i] * fbe_0 - 4.0 * g_zzz_yyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_yyz_1[i] * fe_0 + g_zzzz_yyzz_1[i] * pa_z[i];

        g_zzzzz_yzzz_0[i] = 4.0 * g_zzz_yzzz_0[i] * fbe_0 - 4.0 * g_zzz_yzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_yzz_1[i] * fe_0 + g_zzzz_yzzz_1[i] * pa_z[i];

        g_zzzzz_zzzz_0[i] = 4.0 * g_zzz_zzzz_0[i] * fbe_0 - 4.0 * g_zzz_zzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_zzz_1[i] * fe_0 + g_zzzz_zzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

