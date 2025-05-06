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

#include "KineticEnergyPrimRecGG.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_gg(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_gg,
                            const size_t              idx_ovl_dg,
                            const size_t              idx_kin_dg,
                            const size_t              idx_kin_ff,
                            const size_t              idx_kin_fg,
                            const size_t              idx_ovl_gg,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_ovl_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_ovl_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_ovl_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_ovl_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_ovl_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_ovl_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_ovl_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_ovl_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_ovl_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_ovl_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_ovl_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_ovl_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_ovl_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_ovl_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_ovl_dg + 14);

    auto ts_yy_xxxx = pbuffer.data(idx_ovl_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_ovl_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_ovl_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_ovl_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_ovl_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_ovl_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_ovl_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_ovl_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_ovl_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_ovl_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_ovl_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_ovl_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_ovl_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_ovl_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_ovl_dg + 59);

    auto ts_zz_xxxx = pbuffer.data(idx_ovl_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_ovl_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_ovl_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_ovl_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_ovl_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_ovl_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_ovl_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_ovl_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_ovl_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_ovl_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_ovl_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_ovl_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_ovl_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_ovl_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_ovl_dg + 89);

    // Set up components of auxiliary buffer : DG

    auto tk_xx_xxxx = pbuffer.data(idx_kin_dg);

    auto tk_xx_xxxy = pbuffer.data(idx_kin_dg + 1);

    auto tk_xx_xxxz = pbuffer.data(idx_kin_dg + 2);

    auto tk_xx_xxyy = pbuffer.data(idx_kin_dg + 3);

    auto tk_xx_xxyz = pbuffer.data(idx_kin_dg + 4);

    auto tk_xx_xxzz = pbuffer.data(idx_kin_dg + 5);

    auto tk_xx_xyyy = pbuffer.data(idx_kin_dg + 6);

    auto tk_xx_xyyz = pbuffer.data(idx_kin_dg + 7);

    auto tk_xx_xyzz = pbuffer.data(idx_kin_dg + 8);

    auto tk_xx_xzzz = pbuffer.data(idx_kin_dg + 9);

    auto tk_xx_yyyy = pbuffer.data(idx_kin_dg + 10);

    auto tk_xx_yyyz = pbuffer.data(idx_kin_dg + 11);

    auto tk_xx_yyzz = pbuffer.data(idx_kin_dg + 12);

    auto tk_xx_yzzz = pbuffer.data(idx_kin_dg + 13);

    auto tk_xx_zzzz = pbuffer.data(idx_kin_dg + 14);

    auto tk_yy_xxxx = pbuffer.data(idx_kin_dg + 45);

    auto tk_yy_xxxy = pbuffer.data(idx_kin_dg + 46);

    auto tk_yy_xxxz = pbuffer.data(idx_kin_dg + 47);

    auto tk_yy_xxyy = pbuffer.data(idx_kin_dg + 48);

    auto tk_yy_xxyz = pbuffer.data(idx_kin_dg + 49);

    auto tk_yy_xxzz = pbuffer.data(idx_kin_dg + 50);

    auto tk_yy_xyyy = pbuffer.data(idx_kin_dg + 51);

    auto tk_yy_xyyz = pbuffer.data(idx_kin_dg + 52);

    auto tk_yy_xyzz = pbuffer.data(idx_kin_dg + 53);

    auto tk_yy_xzzz = pbuffer.data(idx_kin_dg + 54);

    auto tk_yy_yyyy = pbuffer.data(idx_kin_dg + 55);

    auto tk_yy_yyyz = pbuffer.data(idx_kin_dg + 56);

    auto tk_yy_yyzz = pbuffer.data(idx_kin_dg + 57);

    auto tk_yy_yzzz = pbuffer.data(idx_kin_dg + 58);

    auto tk_yy_zzzz = pbuffer.data(idx_kin_dg + 59);

    auto tk_zz_xxxx = pbuffer.data(idx_kin_dg + 75);

    auto tk_zz_xxxy = pbuffer.data(idx_kin_dg + 76);

    auto tk_zz_xxxz = pbuffer.data(idx_kin_dg + 77);

    auto tk_zz_xxyy = pbuffer.data(idx_kin_dg + 78);

    auto tk_zz_xxyz = pbuffer.data(idx_kin_dg + 79);

    auto tk_zz_xxzz = pbuffer.data(idx_kin_dg + 80);

    auto tk_zz_xyyy = pbuffer.data(idx_kin_dg + 81);

    auto tk_zz_xyyz = pbuffer.data(idx_kin_dg + 82);

    auto tk_zz_xyzz = pbuffer.data(idx_kin_dg + 83);

    auto tk_zz_xzzz = pbuffer.data(idx_kin_dg + 84);

    auto tk_zz_yyyy = pbuffer.data(idx_kin_dg + 85);

    auto tk_zz_yyyz = pbuffer.data(idx_kin_dg + 86);

    auto tk_zz_yyzz = pbuffer.data(idx_kin_dg + 87);

    auto tk_zz_yzzz = pbuffer.data(idx_kin_dg + 88);

    auto tk_zz_zzzz = pbuffer.data(idx_kin_dg + 89);

    // Set up components of auxiliary buffer : FF

    auto tk_xxx_xxx = pbuffer.data(idx_kin_ff);

    auto tk_xxx_xxy = pbuffer.data(idx_kin_ff + 1);

    auto tk_xxx_xxz = pbuffer.data(idx_kin_ff + 2);

    auto tk_xxx_xyy = pbuffer.data(idx_kin_ff + 3);

    auto tk_xxx_xyz = pbuffer.data(idx_kin_ff + 4);

    auto tk_xxx_xzz = pbuffer.data(idx_kin_ff + 5);

    auto tk_xxx_yyy = pbuffer.data(idx_kin_ff + 6);

    auto tk_xxx_yyz = pbuffer.data(idx_kin_ff + 7);

    auto tk_xxx_yzz = pbuffer.data(idx_kin_ff + 8);

    auto tk_xxx_zzz = pbuffer.data(idx_kin_ff + 9);

    auto tk_xxz_xxz = pbuffer.data(idx_kin_ff + 22);

    auto tk_xxz_xyz = pbuffer.data(idx_kin_ff + 24);

    auto tk_xxz_xzz = pbuffer.data(idx_kin_ff + 25);

    auto tk_xxz_yyz = pbuffer.data(idx_kin_ff + 27);

    auto tk_xxz_yzz = pbuffer.data(idx_kin_ff + 28);

    auto tk_xxz_zzz = pbuffer.data(idx_kin_ff + 29);

    auto tk_xyy_xxy = pbuffer.data(idx_kin_ff + 31);

    auto tk_xyy_xyy = pbuffer.data(idx_kin_ff + 33);

    auto tk_xyy_xyz = pbuffer.data(idx_kin_ff + 34);

    auto tk_xyy_yyy = pbuffer.data(idx_kin_ff + 36);

    auto tk_xyy_yyz = pbuffer.data(idx_kin_ff + 37);

    auto tk_xyy_yzz = pbuffer.data(idx_kin_ff + 38);

    auto tk_xzz_xxz = pbuffer.data(idx_kin_ff + 52);

    auto tk_xzz_xyz = pbuffer.data(idx_kin_ff + 54);

    auto tk_xzz_xzz = pbuffer.data(idx_kin_ff + 55);

    auto tk_xzz_yyz = pbuffer.data(idx_kin_ff + 57);

    auto tk_xzz_yzz = pbuffer.data(idx_kin_ff + 58);

    auto tk_xzz_zzz = pbuffer.data(idx_kin_ff + 59);

    auto tk_yyy_xxx = pbuffer.data(idx_kin_ff + 60);

    auto tk_yyy_xxy = pbuffer.data(idx_kin_ff + 61);

    auto tk_yyy_xxz = pbuffer.data(idx_kin_ff + 62);

    auto tk_yyy_xyy = pbuffer.data(idx_kin_ff + 63);

    auto tk_yyy_xyz = pbuffer.data(idx_kin_ff + 64);

    auto tk_yyy_xzz = pbuffer.data(idx_kin_ff + 65);

    auto tk_yyy_yyy = pbuffer.data(idx_kin_ff + 66);

    auto tk_yyy_yyz = pbuffer.data(idx_kin_ff + 67);

    auto tk_yyy_yzz = pbuffer.data(idx_kin_ff + 68);

    auto tk_yyy_zzz = pbuffer.data(idx_kin_ff + 69);

    auto tk_yyz_xxz = pbuffer.data(idx_kin_ff + 72);

    auto tk_yyz_xyz = pbuffer.data(idx_kin_ff + 74);

    auto tk_yyz_xzz = pbuffer.data(idx_kin_ff + 75);

    auto tk_yyz_yyz = pbuffer.data(idx_kin_ff + 77);

    auto tk_yyz_yzz = pbuffer.data(idx_kin_ff + 78);

    auto tk_yyz_zzz = pbuffer.data(idx_kin_ff + 79);

    auto tk_yzz_xxy = pbuffer.data(idx_kin_ff + 81);

    auto tk_yzz_xxz = pbuffer.data(idx_kin_ff + 82);

    auto tk_yzz_xyy = pbuffer.data(idx_kin_ff + 83);

    auto tk_yzz_xyz = pbuffer.data(idx_kin_ff + 84);

    auto tk_yzz_xzz = pbuffer.data(idx_kin_ff + 85);

    auto tk_yzz_yyy = pbuffer.data(idx_kin_ff + 86);

    auto tk_yzz_yyz = pbuffer.data(idx_kin_ff + 87);

    auto tk_yzz_yzz = pbuffer.data(idx_kin_ff + 88);

    auto tk_yzz_zzz = pbuffer.data(idx_kin_ff + 89);

    auto tk_zzz_xxx = pbuffer.data(idx_kin_ff + 90);

    auto tk_zzz_xxy = pbuffer.data(idx_kin_ff + 91);

    auto tk_zzz_xxz = pbuffer.data(idx_kin_ff + 92);

    auto tk_zzz_xyy = pbuffer.data(idx_kin_ff + 93);

    auto tk_zzz_xyz = pbuffer.data(idx_kin_ff + 94);

    auto tk_zzz_xzz = pbuffer.data(idx_kin_ff + 95);

    auto tk_zzz_yyy = pbuffer.data(idx_kin_ff + 96);

    auto tk_zzz_yyz = pbuffer.data(idx_kin_ff + 97);

    auto tk_zzz_yzz = pbuffer.data(idx_kin_ff + 98);

    auto tk_zzz_zzz = pbuffer.data(idx_kin_ff + 99);

    // Set up components of auxiliary buffer : FG

    auto tk_xxx_xxxx = pbuffer.data(idx_kin_fg);

    auto tk_xxx_xxxy = pbuffer.data(idx_kin_fg + 1);

    auto tk_xxx_xxxz = pbuffer.data(idx_kin_fg + 2);

    auto tk_xxx_xxyy = pbuffer.data(idx_kin_fg + 3);

    auto tk_xxx_xxyz = pbuffer.data(idx_kin_fg + 4);

    auto tk_xxx_xxzz = pbuffer.data(idx_kin_fg + 5);

    auto tk_xxx_xyyy = pbuffer.data(idx_kin_fg + 6);

    auto tk_xxx_xyyz = pbuffer.data(idx_kin_fg + 7);

    auto tk_xxx_xyzz = pbuffer.data(idx_kin_fg + 8);

    auto tk_xxx_xzzz = pbuffer.data(idx_kin_fg + 9);

    auto tk_xxx_yyyy = pbuffer.data(idx_kin_fg + 10);

    auto tk_xxx_yyyz = pbuffer.data(idx_kin_fg + 11);

    auto tk_xxx_yyzz = pbuffer.data(idx_kin_fg + 12);

    auto tk_xxx_yzzz = pbuffer.data(idx_kin_fg + 13);

    auto tk_xxx_zzzz = pbuffer.data(idx_kin_fg + 14);

    auto tk_xxy_xxxx = pbuffer.data(idx_kin_fg + 15);

    auto tk_xxy_xxxy = pbuffer.data(idx_kin_fg + 16);

    auto tk_xxy_xxxz = pbuffer.data(idx_kin_fg + 17);

    auto tk_xxy_xxyy = pbuffer.data(idx_kin_fg + 18);

    auto tk_xxy_xxzz = pbuffer.data(idx_kin_fg + 20);

    auto tk_xxy_xyyy = pbuffer.data(idx_kin_fg + 21);

    auto tk_xxy_xzzz = pbuffer.data(idx_kin_fg + 24);

    auto tk_xxy_yyyy = pbuffer.data(idx_kin_fg + 25);

    auto tk_xxz_xxxx = pbuffer.data(idx_kin_fg + 30);

    auto tk_xxz_xxxy = pbuffer.data(idx_kin_fg + 31);

    auto tk_xxz_xxxz = pbuffer.data(idx_kin_fg + 32);

    auto tk_xxz_xxyy = pbuffer.data(idx_kin_fg + 33);

    auto tk_xxz_xxyz = pbuffer.data(idx_kin_fg + 34);

    auto tk_xxz_xxzz = pbuffer.data(idx_kin_fg + 35);

    auto tk_xxz_xyyy = pbuffer.data(idx_kin_fg + 36);

    auto tk_xxz_xyyz = pbuffer.data(idx_kin_fg + 37);

    auto tk_xxz_xyzz = pbuffer.data(idx_kin_fg + 38);

    auto tk_xxz_xzzz = pbuffer.data(idx_kin_fg + 39);

    auto tk_xxz_yyyz = pbuffer.data(idx_kin_fg + 41);

    auto tk_xxz_yyzz = pbuffer.data(idx_kin_fg + 42);

    auto tk_xxz_yzzz = pbuffer.data(idx_kin_fg + 43);

    auto tk_xxz_zzzz = pbuffer.data(idx_kin_fg + 44);

    auto tk_xyy_xxxx = pbuffer.data(idx_kin_fg + 45);

    auto tk_xyy_xxxy = pbuffer.data(idx_kin_fg + 46);

    auto tk_xyy_xxyy = pbuffer.data(idx_kin_fg + 48);

    auto tk_xyy_xxyz = pbuffer.data(idx_kin_fg + 49);

    auto tk_xyy_xyyy = pbuffer.data(idx_kin_fg + 51);

    auto tk_xyy_xyyz = pbuffer.data(idx_kin_fg + 52);

    auto tk_xyy_xyzz = pbuffer.data(idx_kin_fg + 53);

    auto tk_xyy_yyyy = pbuffer.data(idx_kin_fg + 55);

    auto tk_xyy_yyyz = pbuffer.data(idx_kin_fg + 56);

    auto tk_xyy_yyzz = pbuffer.data(idx_kin_fg + 57);

    auto tk_xyy_yzzz = pbuffer.data(idx_kin_fg + 58);

    auto tk_xyy_zzzz = pbuffer.data(idx_kin_fg + 59);

    auto tk_xzz_xxxx = pbuffer.data(idx_kin_fg + 75);

    auto tk_xzz_xxxz = pbuffer.data(idx_kin_fg + 77);

    auto tk_xzz_xxyz = pbuffer.data(idx_kin_fg + 79);

    auto tk_xzz_xxzz = pbuffer.data(idx_kin_fg + 80);

    auto tk_xzz_xyyz = pbuffer.data(idx_kin_fg + 82);

    auto tk_xzz_xyzz = pbuffer.data(idx_kin_fg + 83);

    auto tk_xzz_xzzz = pbuffer.data(idx_kin_fg + 84);

    auto tk_xzz_yyyy = pbuffer.data(idx_kin_fg + 85);

    auto tk_xzz_yyyz = pbuffer.data(idx_kin_fg + 86);

    auto tk_xzz_yyzz = pbuffer.data(idx_kin_fg + 87);

    auto tk_xzz_yzzz = pbuffer.data(idx_kin_fg + 88);

    auto tk_xzz_zzzz = pbuffer.data(idx_kin_fg + 89);

    auto tk_yyy_xxxx = pbuffer.data(idx_kin_fg + 90);

    auto tk_yyy_xxxy = pbuffer.data(idx_kin_fg + 91);

    auto tk_yyy_xxxz = pbuffer.data(idx_kin_fg + 92);

    auto tk_yyy_xxyy = pbuffer.data(idx_kin_fg + 93);

    auto tk_yyy_xxyz = pbuffer.data(idx_kin_fg + 94);

    auto tk_yyy_xxzz = pbuffer.data(idx_kin_fg + 95);

    auto tk_yyy_xyyy = pbuffer.data(idx_kin_fg + 96);

    auto tk_yyy_xyyz = pbuffer.data(idx_kin_fg + 97);

    auto tk_yyy_xyzz = pbuffer.data(idx_kin_fg + 98);

    auto tk_yyy_xzzz = pbuffer.data(idx_kin_fg + 99);

    auto tk_yyy_yyyy = pbuffer.data(idx_kin_fg + 100);

    auto tk_yyy_yyyz = pbuffer.data(idx_kin_fg + 101);

    auto tk_yyy_yyzz = pbuffer.data(idx_kin_fg + 102);

    auto tk_yyy_yzzz = pbuffer.data(idx_kin_fg + 103);

    auto tk_yyy_zzzz = pbuffer.data(idx_kin_fg + 104);

    auto tk_yyz_xxxy = pbuffer.data(idx_kin_fg + 106);

    auto tk_yyz_xxxz = pbuffer.data(idx_kin_fg + 107);

    auto tk_yyz_xxyy = pbuffer.data(idx_kin_fg + 108);

    auto tk_yyz_xxyz = pbuffer.data(idx_kin_fg + 109);

    auto tk_yyz_xxzz = pbuffer.data(idx_kin_fg + 110);

    auto tk_yyz_xyyy = pbuffer.data(idx_kin_fg + 111);

    auto tk_yyz_xyyz = pbuffer.data(idx_kin_fg + 112);

    auto tk_yyz_xyzz = pbuffer.data(idx_kin_fg + 113);

    auto tk_yyz_xzzz = pbuffer.data(idx_kin_fg + 114);

    auto tk_yyz_yyyy = pbuffer.data(idx_kin_fg + 115);

    auto tk_yyz_yyyz = pbuffer.data(idx_kin_fg + 116);

    auto tk_yyz_yyzz = pbuffer.data(idx_kin_fg + 117);

    auto tk_yyz_yzzz = pbuffer.data(idx_kin_fg + 118);

    auto tk_yyz_zzzz = pbuffer.data(idx_kin_fg + 119);

    auto tk_yzz_xxxx = pbuffer.data(idx_kin_fg + 120);

    auto tk_yzz_xxxy = pbuffer.data(idx_kin_fg + 121);

    auto tk_yzz_xxxz = pbuffer.data(idx_kin_fg + 122);

    auto tk_yzz_xxyy = pbuffer.data(idx_kin_fg + 123);

    auto tk_yzz_xxyz = pbuffer.data(idx_kin_fg + 124);

    auto tk_yzz_xxzz = pbuffer.data(idx_kin_fg + 125);

    auto tk_yzz_xyyy = pbuffer.data(idx_kin_fg + 126);

    auto tk_yzz_xyyz = pbuffer.data(idx_kin_fg + 127);

    auto tk_yzz_xyzz = pbuffer.data(idx_kin_fg + 128);

    auto tk_yzz_xzzz = pbuffer.data(idx_kin_fg + 129);

    auto tk_yzz_yyyy = pbuffer.data(idx_kin_fg + 130);

    auto tk_yzz_yyyz = pbuffer.data(idx_kin_fg + 131);

    auto tk_yzz_yyzz = pbuffer.data(idx_kin_fg + 132);

    auto tk_yzz_yzzz = pbuffer.data(idx_kin_fg + 133);

    auto tk_yzz_zzzz = pbuffer.data(idx_kin_fg + 134);

    auto tk_zzz_xxxx = pbuffer.data(idx_kin_fg + 135);

    auto tk_zzz_xxxy = pbuffer.data(idx_kin_fg + 136);

    auto tk_zzz_xxxz = pbuffer.data(idx_kin_fg + 137);

    auto tk_zzz_xxyy = pbuffer.data(idx_kin_fg + 138);

    auto tk_zzz_xxyz = pbuffer.data(idx_kin_fg + 139);

    auto tk_zzz_xxzz = pbuffer.data(idx_kin_fg + 140);

    auto tk_zzz_xyyy = pbuffer.data(idx_kin_fg + 141);

    auto tk_zzz_xyyz = pbuffer.data(idx_kin_fg + 142);

    auto tk_zzz_xyzz = pbuffer.data(idx_kin_fg + 143);

    auto tk_zzz_xzzz = pbuffer.data(idx_kin_fg + 144);

    auto tk_zzz_yyyy = pbuffer.data(idx_kin_fg + 145);

    auto tk_zzz_yyyz = pbuffer.data(idx_kin_fg + 146);

    auto tk_zzz_yyzz = pbuffer.data(idx_kin_fg + 147);

    auto tk_zzz_yzzz = pbuffer.data(idx_kin_fg + 148);

    auto tk_zzz_zzzz = pbuffer.data(idx_kin_fg + 149);

    // Set up components of auxiliary buffer : GG

    auto ts_xxxx_xxxx = pbuffer.data(idx_ovl_gg);

    auto ts_xxxx_xxxy = pbuffer.data(idx_ovl_gg + 1);

    auto ts_xxxx_xxxz = pbuffer.data(idx_ovl_gg + 2);

    auto ts_xxxx_xxyy = pbuffer.data(idx_ovl_gg + 3);

    auto ts_xxxx_xxyz = pbuffer.data(idx_ovl_gg + 4);

    auto ts_xxxx_xxzz = pbuffer.data(idx_ovl_gg + 5);

    auto ts_xxxx_xyyy = pbuffer.data(idx_ovl_gg + 6);

    auto ts_xxxx_xyyz = pbuffer.data(idx_ovl_gg + 7);

    auto ts_xxxx_xyzz = pbuffer.data(idx_ovl_gg + 8);

    auto ts_xxxx_xzzz = pbuffer.data(idx_ovl_gg + 9);

    auto ts_xxxx_yyyy = pbuffer.data(idx_ovl_gg + 10);

    auto ts_xxxx_yyyz = pbuffer.data(idx_ovl_gg + 11);

    auto ts_xxxx_yyzz = pbuffer.data(idx_ovl_gg + 12);

    auto ts_xxxx_yzzz = pbuffer.data(idx_ovl_gg + 13);

    auto ts_xxxx_zzzz = pbuffer.data(idx_ovl_gg + 14);

    auto ts_xxxy_xxxx = pbuffer.data(idx_ovl_gg + 15);

    auto ts_xxxy_xxxy = pbuffer.data(idx_ovl_gg + 16);

    auto ts_xxxy_xxxz = pbuffer.data(idx_ovl_gg + 17);

    auto ts_xxxy_xxyy = pbuffer.data(idx_ovl_gg + 18);

    auto ts_xxxy_xxyz = pbuffer.data(idx_ovl_gg + 19);

    auto ts_xxxy_xxzz = pbuffer.data(idx_ovl_gg + 20);

    auto ts_xxxy_xyyy = pbuffer.data(idx_ovl_gg + 21);

    auto ts_xxxy_xyyz = pbuffer.data(idx_ovl_gg + 22);

    auto ts_xxxy_xyzz = pbuffer.data(idx_ovl_gg + 23);

    auto ts_xxxy_xzzz = pbuffer.data(idx_ovl_gg + 24);

    auto ts_xxxy_yyyy = pbuffer.data(idx_ovl_gg + 25);

    auto ts_xxxy_yyyz = pbuffer.data(idx_ovl_gg + 26);

    auto ts_xxxy_yyzz = pbuffer.data(idx_ovl_gg + 27);

    auto ts_xxxy_yzzz = pbuffer.data(idx_ovl_gg + 28);

    auto ts_xxxy_zzzz = pbuffer.data(idx_ovl_gg + 29);

    auto ts_xxxz_xxxx = pbuffer.data(idx_ovl_gg + 30);

    auto ts_xxxz_xxxy = pbuffer.data(idx_ovl_gg + 31);

    auto ts_xxxz_xxxz = pbuffer.data(idx_ovl_gg + 32);

    auto ts_xxxz_xxyy = pbuffer.data(idx_ovl_gg + 33);

    auto ts_xxxz_xxyz = pbuffer.data(idx_ovl_gg + 34);

    auto ts_xxxz_xxzz = pbuffer.data(idx_ovl_gg + 35);

    auto ts_xxxz_xyyy = pbuffer.data(idx_ovl_gg + 36);

    auto ts_xxxz_xyyz = pbuffer.data(idx_ovl_gg + 37);

    auto ts_xxxz_xyzz = pbuffer.data(idx_ovl_gg + 38);

    auto ts_xxxz_xzzz = pbuffer.data(idx_ovl_gg + 39);

    auto ts_xxxz_yyyy = pbuffer.data(idx_ovl_gg + 40);

    auto ts_xxxz_yyyz = pbuffer.data(idx_ovl_gg + 41);

    auto ts_xxxz_yyzz = pbuffer.data(idx_ovl_gg + 42);

    auto ts_xxxz_yzzz = pbuffer.data(idx_ovl_gg + 43);

    auto ts_xxxz_zzzz = pbuffer.data(idx_ovl_gg + 44);

    auto ts_xxyy_xxxx = pbuffer.data(idx_ovl_gg + 45);

    auto ts_xxyy_xxxy = pbuffer.data(idx_ovl_gg + 46);

    auto ts_xxyy_xxxz = pbuffer.data(idx_ovl_gg + 47);

    auto ts_xxyy_xxyy = pbuffer.data(idx_ovl_gg + 48);

    auto ts_xxyy_xxyz = pbuffer.data(idx_ovl_gg + 49);

    auto ts_xxyy_xxzz = pbuffer.data(idx_ovl_gg + 50);

    auto ts_xxyy_xyyy = pbuffer.data(idx_ovl_gg + 51);

    auto ts_xxyy_xyyz = pbuffer.data(idx_ovl_gg + 52);

    auto ts_xxyy_xyzz = pbuffer.data(idx_ovl_gg + 53);

    auto ts_xxyy_xzzz = pbuffer.data(idx_ovl_gg + 54);

    auto ts_xxyy_yyyy = pbuffer.data(idx_ovl_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_ovl_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_ovl_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_ovl_gg + 58);

    auto ts_xxyy_zzzz = pbuffer.data(idx_ovl_gg + 59);

    auto ts_xxyz_xxxx = pbuffer.data(idx_ovl_gg + 60);

    auto ts_xxyz_xxxy = pbuffer.data(idx_ovl_gg + 61);

    auto ts_xxyz_xxxz = pbuffer.data(idx_ovl_gg + 62);

    auto ts_xxyz_xxyy = pbuffer.data(idx_ovl_gg + 63);

    auto ts_xxyz_xxyz = pbuffer.data(idx_ovl_gg + 64);

    auto ts_xxyz_xxzz = pbuffer.data(idx_ovl_gg + 65);

    auto ts_xxyz_xyyy = pbuffer.data(idx_ovl_gg + 66);

    auto ts_xxyz_xyyz = pbuffer.data(idx_ovl_gg + 67);

    auto ts_xxyz_xyzz = pbuffer.data(idx_ovl_gg + 68);

    auto ts_xxyz_xzzz = pbuffer.data(idx_ovl_gg + 69);

    auto ts_xxyz_yyyy = pbuffer.data(idx_ovl_gg + 70);

    auto ts_xxyz_yyyz = pbuffer.data(idx_ovl_gg + 71);

    auto ts_xxyz_yyzz = pbuffer.data(idx_ovl_gg + 72);

    auto ts_xxyz_yzzz = pbuffer.data(idx_ovl_gg + 73);

    auto ts_xxyz_zzzz = pbuffer.data(idx_ovl_gg + 74);

    auto ts_xxzz_xxxx = pbuffer.data(idx_ovl_gg + 75);

    auto ts_xxzz_xxxy = pbuffer.data(idx_ovl_gg + 76);

    auto ts_xxzz_xxxz = pbuffer.data(idx_ovl_gg + 77);

    auto ts_xxzz_xxyy = pbuffer.data(idx_ovl_gg + 78);

    auto ts_xxzz_xxyz = pbuffer.data(idx_ovl_gg + 79);

    auto ts_xxzz_xxzz = pbuffer.data(idx_ovl_gg + 80);

    auto ts_xxzz_xyyy = pbuffer.data(idx_ovl_gg + 81);

    auto ts_xxzz_xyyz = pbuffer.data(idx_ovl_gg + 82);

    auto ts_xxzz_xyzz = pbuffer.data(idx_ovl_gg + 83);

    auto ts_xxzz_xzzz = pbuffer.data(idx_ovl_gg + 84);

    auto ts_xxzz_yyyy = pbuffer.data(idx_ovl_gg + 85);

    auto ts_xxzz_yyyz = pbuffer.data(idx_ovl_gg + 86);

    auto ts_xxzz_yyzz = pbuffer.data(idx_ovl_gg + 87);

    auto ts_xxzz_yzzz = pbuffer.data(idx_ovl_gg + 88);

    auto ts_xxzz_zzzz = pbuffer.data(idx_ovl_gg + 89);

    auto ts_xyyy_xxxx = pbuffer.data(idx_ovl_gg + 90);

    auto ts_xyyy_xxxy = pbuffer.data(idx_ovl_gg + 91);

    auto ts_xyyy_xxxz = pbuffer.data(idx_ovl_gg + 92);

    auto ts_xyyy_xxyy = pbuffer.data(idx_ovl_gg + 93);

    auto ts_xyyy_xxyz = pbuffer.data(idx_ovl_gg + 94);

    auto ts_xyyy_xxzz = pbuffer.data(idx_ovl_gg + 95);

    auto ts_xyyy_xyyy = pbuffer.data(idx_ovl_gg + 96);

    auto ts_xyyy_xyyz = pbuffer.data(idx_ovl_gg + 97);

    auto ts_xyyy_xyzz = pbuffer.data(idx_ovl_gg + 98);

    auto ts_xyyy_xzzz = pbuffer.data(idx_ovl_gg + 99);

    auto ts_xyyy_yyyy = pbuffer.data(idx_ovl_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_ovl_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_ovl_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_ovl_gg + 103);

    auto ts_xyyy_zzzz = pbuffer.data(idx_ovl_gg + 104);

    auto ts_xyyz_xxxx = pbuffer.data(idx_ovl_gg + 105);

    auto ts_xyyz_xxxy = pbuffer.data(idx_ovl_gg + 106);

    auto ts_xyyz_xxxz = pbuffer.data(idx_ovl_gg + 107);

    auto ts_xyyz_xxyy = pbuffer.data(idx_ovl_gg + 108);

    auto ts_xyyz_xxyz = pbuffer.data(idx_ovl_gg + 109);

    auto ts_xyyz_xxzz = pbuffer.data(idx_ovl_gg + 110);

    auto ts_xyyz_xyyy = pbuffer.data(idx_ovl_gg + 111);

    auto ts_xyyz_xyyz = pbuffer.data(idx_ovl_gg + 112);

    auto ts_xyyz_xyzz = pbuffer.data(idx_ovl_gg + 113);

    auto ts_xyyz_xzzz = pbuffer.data(idx_ovl_gg + 114);

    auto ts_xyyz_yyyy = pbuffer.data(idx_ovl_gg + 115);

    auto ts_xyyz_yyyz = pbuffer.data(idx_ovl_gg + 116);

    auto ts_xyyz_yyzz = pbuffer.data(idx_ovl_gg + 117);

    auto ts_xyyz_yzzz = pbuffer.data(idx_ovl_gg + 118);

    auto ts_xyyz_zzzz = pbuffer.data(idx_ovl_gg + 119);

    auto ts_xyzz_xxxx = pbuffer.data(idx_ovl_gg + 120);

    auto ts_xyzz_xxxy = pbuffer.data(idx_ovl_gg + 121);

    auto ts_xyzz_xxxz = pbuffer.data(idx_ovl_gg + 122);

    auto ts_xyzz_xxyy = pbuffer.data(idx_ovl_gg + 123);

    auto ts_xyzz_xxyz = pbuffer.data(idx_ovl_gg + 124);

    auto ts_xyzz_xxzz = pbuffer.data(idx_ovl_gg + 125);

    auto ts_xyzz_xyyy = pbuffer.data(idx_ovl_gg + 126);

    auto ts_xyzz_xyyz = pbuffer.data(idx_ovl_gg + 127);

    auto ts_xyzz_xyzz = pbuffer.data(idx_ovl_gg + 128);

    auto ts_xyzz_xzzz = pbuffer.data(idx_ovl_gg + 129);

    auto ts_xyzz_yyyy = pbuffer.data(idx_ovl_gg + 130);

    auto ts_xyzz_yyyz = pbuffer.data(idx_ovl_gg + 131);

    auto ts_xyzz_yyzz = pbuffer.data(idx_ovl_gg + 132);

    auto ts_xyzz_yzzz = pbuffer.data(idx_ovl_gg + 133);

    auto ts_xyzz_zzzz = pbuffer.data(idx_ovl_gg + 134);

    auto ts_xzzz_xxxx = pbuffer.data(idx_ovl_gg + 135);

    auto ts_xzzz_xxxy = pbuffer.data(idx_ovl_gg + 136);

    auto ts_xzzz_xxxz = pbuffer.data(idx_ovl_gg + 137);

    auto ts_xzzz_xxyy = pbuffer.data(idx_ovl_gg + 138);

    auto ts_xzzz_xxyz = pbuffer.data(idx_ovl_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_ovl_gg + 140);

    auto ts_xzzz_xyyy = pbuffer.data(idx_ovl_gg + 141);

    auto ts_xzzz_xyyz = pbuffer.data(idx_ovl_gg + 142);

    auto ts_xzzz_xyzz = pbuffer.data(idx_ovl_gg + 143);

    auto ts_xzzz_xzzz = pbuffer.data(idx_ovl_gg + 144);

    auto ts_xzzz_yyyy = pbuffer.data(idx_ovl_gg + 145);

    auto ts_xzzz_yyyz = pbuffer.data(idx_ovl_gg + 146);

    auto ts_xzzz_yyzz = pbuffer.data(idx_ovl_gg + 147);

    auto ts_xzzz_yzzz = pbuffer.data(idx_ovl_gg + 148);

    auto ts_xzzz_zzzz = pbuffer.data(idx_ovl_gg + 149);

    auto ts_yyyy_xxxx = pbuffer.data(idx_ovl_gg + 150);

    auto ts_yyyy_xxxy = pbuffer.data(idx_ovl_gg + 151);

    auto ts_yyyy_xxxz = pbuffer.data(idx_ovl_gg + 152);

    auto ts_yyyy_xxyy = pbuffer.data(idx_ovl_gg + 153);

    auto ts_yyyy_xxyz = pbuffer.data(idx_ovl_gg + 154);

    auto ts_yyyy_xxzz = pbuffer.data(idx_ovl_gg + 155);

    auto ts_yyyy_xyyy = pbuffer.data(idx_ovl_gg + 156);

    auto ts_yyyy_xyyz = pbuffer.data(idx_ovl_gg + 157);

    auto ts_yyyy_xyzz = pbuffer.data(idx_ovl_gg + 158);

    auto ts_yyyy_xzzz = pbuffer.data(idx_ovl_gg + 159);

    auto ts_yyyy_yyyy = pbuffer.data(idx_ovl_gg + 160);

    auto ts_yyyy_yyyz = pbuffer.data(idx_ovl_gg + 161);

    auto ts_yyyy_yyzz = pbuffer.data(idx_ovl_gg + 162);

    auto ts_yyyy_yzzz = pbuffer.data(idx_ovl_gg + 163);

    auto ts_yyyy_zzzz = pbuffer.data(idx_ovl_gg + 164);

    auto ts_yyyz_xxxx = pbuffer.data(idx_ovl_gg + 165);

    auto ts_yyyz_xxxy = pbuffer.data(idx_ovl_gg + 166);

    auto ts_yyyz_xxxz = pbuffer.data(idx_ovl_gg + 167);

    auto ts_yyyz_xxyy = pbuffer.data(idx_ovl_gg + 168);

    auto ts_yyyz_xxyz = pbuffer.data(idx_ovl_gg + 169);

    auto ts_yyyz_xxzz = pbuffer.data(idx_ovl_gg + 170);

    auto ts_yyyz_xyyy = pbuffer.data(idx_ovl_gg + 171);

    auto ts_yyyz_xyyz = pbuffer.data(idx_ovl_gg + 172);

    auto ts_yyyz_xyzz = pbuffer.data(idx_ovl_gg + 173);

    auto ts_yyyz_xzzz = pbuffer.data(idx_ovl_gg + 174);

    auto ts_yyyz_yyyy = pbuffer.data(idx_ovl_gg + 175);

    auto ts_yyyz_yyyz = pbuffer.data(idx_ovl_gg + 176);

    auto ts_yyyz_yyzz = pbuffer.data(idx_ovl_gg + 177);

    auto ts_yyyz_yzzz = pbuffer.data(idx_ovl_gg + 178);

    auto ts_yyyz_zzzz = pbuffer.data(idx_ovl_gg + 179);

    auto ts_yyzz_xxxx = pbuffer.data(idx_ovl_gg + 180);

    auto ts_yyzz_xxxy = pbuffer.data(idx_ovl_gg + 181);

    auto ts_yyzz_xxxz = pbuffer.data(idx_ovl_gg + 182);

    auto ts_yyzz_xxyy = pbuffer.data(idx_ovl_gg + 183);

    auto ts_yyzz_xxyz = pbuffer.data(idx_ovl_gg + 184);

    auto ts_yyzz_xxzz = pbuffer.data(idx_ovl_gg + 185);

    auto ts_yyzz_xyyy = pbuffer.data(idx_ovl_gg + 186);

    auto ts_yyzz_xyyz = pbuffer.data(idx_ovl_gg + 187);

    auto ts_yyzz_xyzz = pbuffer.data(idx_ovl_gg + 188);

    auto ts_yyzz_xzzz = pbuffer.data(idx_ovl_gg + 189);

    auto ts_yyzz_yyyy = pbuffer.data(idx_ovl_gg + 190);

    auto ts_yyzz_yyyz = pbuffer.data(idx_ovl_gg + 191);

    auto ts_yyzz_yyzz = pbuffer.data(idx_ovl_gg + 192);

    auto ts_yyzz_yzzz = pbuffer.data(idx_ovl_gg + 193);

    auto ts_yyzz_zzzz = pbuffer.data(idx_ovl_gg + 194);

    auto ts_yzzz_xxxx = pbuffer.data(idx_ovl_gg + 195);

    auto ts_yzzz_xxxy = pbuffer.data(idx_ovl_gg + 196);

    auto ts_yzzz_xxxz = pbuffer.data(idx_ovl_gg + 197);

    auto ts_yzzz_xxyy = pbuffer.data(idx_ovl_gg + 198);

    auto ts_yzzz_xxyz = pbuffer.data(idx_ovl_gg + 199);

    auto ts_yzzz_xxzz = pbuffer.data(idx_ovl_gg + 200);

    auto ts_yzzz_xyyy = pbuffer.data(idx_ovl_gg + 201);

    auto ts_yzzz_xyyz = pbuffer.data(idx_ovl_gg + 202);

    auto ts_yzzz_xyzz = pbuffer.data(idx_ovl_gg + 203);

    auto ts_yzzz_xzzz = pbuffer.data(idx_ovl_gg + 204);

    auto ts_yzzz_yyyy = pbuffer.data(idx_ovl_gg + 205);

    auto ts_yzzz_yyyz = pbuffer.data(idx_ovl_gg + 206);

    auto ts_yzzz_yyzz = pbuffer.data(idx_ovl_gg + 207);

    auto ts_yzzz_yzzz = pbuffer.data(idx_ovl_gg + 208);

    auto ts_yzzz_zzzz = pbuffer.data(idx_ovl_gg + 209);

    auto ts_zzzz_xxxx = pbuffer.data(idx_ovl_gg + 210);

    auto ts_zzzz_xxxy = pbuffer.data(idx_ovl_gg + 211);

    auto ts_zzzz_xxxz = pbuffer.data(idx_ovl_gg + 212);

    auto ts_zzzz_xxyy = pbuffer.data(idx_ovl_gg + 213);

    auto ts_zzzz_xxyz = pbuffer.data(idx_ovl_gg + 214);

    auto ts_zzzz_xxzz = pbuffer.data(idx_ovl_gg + 215);

    auto ts_zzzz_xyyy = pbuffer.data(idx_ovl_gg + 216);

    auto ts_zzzz_xyyz = pbuffer.data(idx_ovl_gg + 217);

    auto ts_zzzz_xyzz = pbuffer.data(idx_ovl_gg + 218);

    auto ts_zzzz_xzzz = pbuffer.data(idx_ovl_gg + 219);

    auto ts_zzzz_yyyy = pbuffer.data(idx_ovl_gg + 220);

    auto ts_zzzz_yyyz = pbuffer.data(idx_ovl_gg + 221);

    auto ts_zzzz_yyzz = pbuffer.data(idx_ovl_gg + 222);

    auto ts_zzzz_yzzz = pbuffer.data(idx_ovl_gg + 223);

    auto ts_zzzz_zzzz = pbuffer.data(idx_ovl_gg + 224);

    // Set up 0-15 components of targeted buffer : GG

    auto tk_xxxx_xxxx = pbuffer.data(idx_kin_gg);

    auto tk_xxxx_xxxy = pbuffer.data(idx_kin_gg + 1);

    auto tk_xxxx_xxxz = pbuffer.data(idx_kin_gg + 2);

    auto tk_xxxx_xxyy = pbuffer.data(idx_kin_gg + 3);

    auto tk_xxxx_xxyz = pbuffer.data(idx_kin_gg + 4);

    auto tk_xxxx_xxzz = pbuffer.data(idx_kin_gg + 5);

    auto tk_xxxx_xyyy = pbuffer.data(idx_kin_gg + 6);

    auto tk_xxxx_xyyz = pbuffer.data(idx_kin_gg + 7);

    auto tk_xxxx_xyzz = pbuffer.data(idx_kin_gg + 8);

    auto tk_xxxx_xzzz = pbuffer.data(idx_kin_gg + 9);

    auto tk_xxxx_yyyy = pbuffer.data(idx_kin_gg + 10);

    auto tk_xxxx_yyyz = pbuffer.data(idx_kin_gg + 11);

    auto tk_xxxx_yyzz = pbuffer.data(idx_kin_gg + 12);

    auto tk_xxxx_yzzz = pbuffer.data(idx_kin_gg + 13);

    auto tk_xxxx_zzzz = pbuffer.data(idx_kin_gg + 14);

#pragma omp simd aligned(pa_x,             \
                             tk_xx_xxxx,   \
                             tk_xx_xxxy,   \
                             tk_xx_xxxz,   \
                             tk_xx_xxyy,   \
                             tk_xx_xxyz,   \
                             tk_xx_xxzz,   \
                             tk_xx_xyyy,   \
                             tk_xx_xyyz,   \
                             tk_xx_xyzz,   \
                             tk_xx_xzzz,   \
                             tk_xx_yyyy,   \
                             tk_xx_yyyz,   \
                             tk_xx_yyzz,   \
                             tk_xx_yzzz,   \
                             tk_xx_zzzz,   \
                             tk_xxx_xxx,   \
                             tk_xxx_xxxx,  \
                             tk_xxx_xxxy,  \
                             tk_xxx_xxxz,  \
                             tk_xxx_xxy,   \
                             tk_xxx_xxyy,  \
                             tk_xxx_xxyz,  \
                             tk_xxx_xxz,   \
                             tk_xxx_xxzz,  \
                             tk_xxx_xyy,   \
                             tk_xxx_xyyy,  \
                             tk_xxx_xyyz,  \
                             tk_xxx_xyz,   \
                             tk_xxx_xyzz,  \
                             tk_xxx_xzz,   \
                             tk_xxx_xzzz,  \
                             tk_xxx_yyy,   \
                             tk_xxx_yyyy,  \
                             tk_xxx_yyyz,  \
                             tk_xxx_yyz,   \
                             tk_xxx_yyzz,  \
                             tk_xxx_yzz,   \
                             tk_xxx_yzzz,  \
                             tk_xxx_zzz,   \
                             tk_xxx_zzzz,  \
                             tk_xxxx_xxxx, \
                             tk_xxxx_xxxy, \
                             tk_xxxx_xxxz, \
                             tk_xxxx_xxyy, \
                             tk_xxxx_xxyz, \
                             tk_xxxx_xxzz, \
                             tk_xxxx_xyyy, \
                             tk_xxxx_xyyz, \
                             tk_xxxx_xyzz, \
                             tk_xxxx_xzzz, \
                             tk_xxxx_yyyy, \
                             tk_xxxx_yyyz, \
                             tk_xxxx_yyzz, \
                             tk_xxxx_yzzz, \
                             tk_xxxx_zzzz, \
                             ts_xx_xxxx,   \
                             ts_xx_xxxy,   \
                             ts_xx_xxxz,   \
                             ts_xx_xxyy,   \
                             ts_xx_xxyz,   \
                             ts_xx_xxzz,   \
                             ts_xx_xyyy,   \
                             ts_xx_xyyz,   \
                             ts_xx_xyzz,   \
                             ts_xx_xzzz,   \
                             ts_xx_yyyy,   \
                             ts_xx_yyyz,   \
                             ts_xx_yyzz,   \
                             ts_xx_yzzz,   \
                             ts_xx_zzzz,   \
                             ts_xxxx_xxxx, \
                             ts_xxxx_xxxy, \
                             ts_xxxx_xxxz, \
                             ts_xxxx_xxyy, \
                             ts_xxxx_xxyz, \
                             ts_xxxx_xxzz, \
                             ts_xxxx_xyyy, \
                             ts_xxxx_xyyz, \
                             ts_xxxx_xyzz, \
                             ts_xxxx_xzzz, \
                             ts_xxxx_yyyy, \
                             ts_xxxx_yyyz, \
                             ts_xxxx_yyzz, \
                             ts_xxxx_yzzz, \
                             ts_xxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxx_xxxx[i] = -6.0 * ts_xx_xxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxx[i] * fe_0 + 4.0 * tk_xxx_xxx[i] * fe_0 + tk_xxx_xxxx[i] * pa_x[i] +
                          2.0 * ts_xxxx_xxxx[i] * fz_0;

        tk_xxxx_xxxy[i] = -6.0 * ts_xx_xxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxy[i] * fe_0 + 3.0 * tk_xxx_xxy[i] * fe_0 + tk_xxx_xxxy[i] * pa_x[i] +
                          2.0 * ts_xxxx_xxxy[i] * fz_0;

        tk_xxxx_xxxz[i] = -6.0 * ts_xx_xxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxz[i] * fe_0 + 3.0 * tk_xxx_xxz[i] * fe_0 + tk_xxx_xxxz[i] * pa_x[i] +
                          2.0 * ts_xxxx_xxxz[i] * fz_0;

        tk_xxxx_xxyy[i] = -6.0 * ts_xx_xxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyy[i] * fe_0 + 2.0 * tk_xxx_xyy[i] * fe_0 + tk_xxx_xxyy[i] * pa_x[i] +
                          2.0 * ts_xxxx_xxyy[i] * fz_0;

        tk_xxxx_xxyz[i] = -6.0 * ts_xx_xxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyz[i] * fe_0 + 2.0 * tk_xxx_xyz[i] * fe_0 + tk_xxx_xxyz[i] * pa_x[i] +
                          2.0 * ts_xxxx_xxyz[i] * fz_0;

        tk_xxxx_xxzz[i] = -6.0 * ts_xx_xxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxzz[i] * fe_0 + 2.0 * tk_xxx_xzz[i] * fe_0 + tk_xxx_xxzz[i] * pa_x[i] +
                          2.0 * ts_xxxx_xxzz[i] * fz_0;

        tk_xxxx_xyyy[i] = -6.0 * ts_xx_xyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyy[i] * fe_0 + tk_xxx_yyy[i] * fe_0 + tk_xxx_xyyy[i] * pa_x[i] +
                          2.0 * ts_xxxx_xyyy[i] * fz_0;

        tk_xxxx_xyyz[i] = -6.0 * ts_xx_xyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyz[i] * fe_0 + tk_xxx_yyz[i] * fe_0 + tk_xxx_xyyz[i] * pa_x[i] +
                          2.0 * ts_xxxx_xyyz[i] * fz_0;

        tk_xxxx_xyzz[i] = -6.0 * ts_xx_xyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyzz[i] * fe_0 + tk_xxx_yzz[i] * fe_0 + tk_xxx_xyzz[i] * pa_x[i] +
                          2.0 * ts_xxxx_xyzz[i] * fz_0;

        tk_xxxx_xzzz[i] = -6.0 * ts_xx_xzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xzzz[i] * fe_0 + tk_xxx_zzz[i] * fe_0 + tk_xxx_xzzz[i] * pa_x[i] +
                          2.0 * ts_xxxx_xzzz[i] * fz_0;

        tk_xxxx_yyyy[i] = -6.0 * ts_xx_yyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyy[i] * fe_0 + tk_xxx_yyyy[i] * pa_x[i] + 2.0 * ts_xxxx_yyyy[i] * fz_0;

        tk_xxxx_yyyz[i] = -6.0 * ts_xx_yyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyz[i] * fe_0 + tk_xxx_yyyz[i] * pa_x[i] + 2.0 * ts_xxxx_yyyz[i] * fz_0;

        tk_xxxx_yyzz[i] = -6.0 * ts_xx_yyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyzz[i] * fe_0 + tk_xxx_yyzz[i] * pa_x[i] + 2.0 * ts_xxxx_yyzz[i] * fz_0;

        tk_xxxx_yzzz[i] = -6.0 * ts_xx_yzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yzzz[i] * fe_0 + tk_xxx_yzzz[i] * pa_x[i] + 2.0 * ts_xxxx_yzzz[i] * fz_0;

        tk_xxxx_zzzz[i] = -6.0 * ts_xx_zzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_zzzz[i] * fe_0 + tk_xxx_zzzz[i] * pa_x[i] + 2.0 * ts_xxxx_zzzz[i] * fz_0;
    }

    // Set up 15-30 components of targeted buffer : GG

    auto tk_xxxy_xxxx = pbuffer.data(idx_kin_gg + 15);

    auto tk_xxxy_xxxy = pbuffer.data(idx_kin_gg + 16);

    auto tk_xxxy_xxxz = pbuffer.data(idx_kin_gg + 17);

    auto tk_xxxy_xxyy = pbuffer.data(idx_kin_gg + 18);

    auto tk_xxxy_xxyz = pbuffer.data(idx_kin_gg + 19);

    auto tk_xxxy_xxzz = pbuffer.data(idx_kin_gg + 20);

    auto tk_xxxy_xyyy = pbuffer.data(idx_kin_gg + 21);

    auto tk_xxxy_xyyz = pbuffer.data(idx_kin_gg + 22);

    auto tk_xxxy_xyzz = pbuffer.data(idx_kin_gg + 23);

    auto tk_xxxy_xzzz = pbuffer.data(idx_kin_gg + 24);

    auto tk_xxxy_yyyy = pbuffer.data(idx_kin_gg + 25);

    auto tk_xxxy_yyyz = pbuffer.data(idx_kin_gg + 26);

    auto tk_xxxy_yyzz = pbuffer.data(idx_kin_gg + 27);

    auto tk_xxxy_yzzz = pbuffer.data(idx_kin_gg + 28);

    auto tk_xxxy_zzzz = pbuffer.data(idx_kin_gg + 29);

#pragma omp simd aligned(pa_y,             \
                             tk_xxx_xxx,   \
                             tk_xxx_xxxx,  \
                             tk_xxx_xxxy,  \
                             tk_xxx_xxxz,  \
                             tk_xxx_xxy,   \
                             tk_xxx_xxyy,  \
                             tk_xxx_xxyz,  \
                             tk_xxx_xxz,   \
                             tk_xxx_xxzz,  \
                             tk_xxx_xyy,   \
                             tk_xxx_xyyy,  \
                             tk_xxx_xyyz,  \
                             tk_xxx_xyz,   \
                             tk_xxx_xyzz,  \
                             tk_xxx_xzz,   \
                             tk_xxx_xzzz,  \
                             tk_xxx_yyy,   \
                             tk_xxx_yyyy,  \
                             tk_xxx_yyyz,  \
                             tk_xxx_yyz,   \
                             tk_xxx_yyzz,  \
                             tk_xxx_yzz,   \
                             tk_xxx_yzzz,  \
                             tk_xxx_zzz,   \
                             tk_xxx_zzzz,  \
                             tk_xxxy_xxxx, \
                             tk_xxxy_xxxy, \
                             tk_xxxy_xxxz, \
                             tk_xxxy_xxyy, \
                             tk_xxxy_xxyz, \
                             tk_xxxy_xxzz, \
                             tk_xxxy_xyyy, \
                             tk_xxxy_xyyz, \
                             tk_xxxy_xyzz, \
                             tk_xxxy_xzzz, \
                             tk_xxxy_yyyy, \
                             tk_xxxy_yyyz, \
                             tk_xxxy_yyzz, \
                             tk_xxxy_yzzz, \
                             tk_xxxy_zzzz, \
                             ts_xxxy_xxxx, \
                             ts_xxxy_xxxy, \
                             ts_xxxy_xxxz, \
                             ts_xxxy_xxyy, \
                             ts_xxxy_xxyz, \
                             ts_xxxy_xxzz, \
                             ts_xxxy_xyyy, \
                             ts_xxxy_xyyz, \
                             ts_xxxy_xyzz, \
                             ts_xxxy_xzzz, \
                             ts_xxxy_yyyy, \
                             ts_xxxy_yyyz, \
                             ts_xxxy_yyzz, \
                             ts_xxxy_yzzz, \
                             ts_xxxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxy_xxxx[i] = tk_xxx_xxxx[i] * pa_y[i] + 2.0 * ts_xxxy_xxxx[i] * fz_0;

        tk_xxxy_xxxy[i] = tk_xxx_xxx[i] * fe_0 + tk_xxx_xxxy[i] * pa_y[i] + 2.0 * ts_xxxy_xxxy[i] * fz_0;

        tk_xxxy_xxxz[i] = tk_xxx_xxxz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxz[i] * fz_0;

        tk_xxxy_xxyy[i] = 2.0 * tk_xxx_xxy[i] * fe_0 + tk_xxx_xxyy[i] * pa_y[i] + 2.0 * ts_xxxy_xxyy[i] * fz_0;

        tk_xxxy_xxyz[i] = tk_xxx_xxz[i] * fe_0 + tk_xxx_xxyz[i] * pa_y[i] + 2.0 * ts_xxxy_xxyz[i] * fz_0;

        tk_xxxy_xxzz[i] = tk_xxx_xxzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxzz[i] * fz_0;

        tk_xxxy_xyyy[i] = 3.0 * tk_xxx_xyy[i] * fe_0 + tk_xxx_xyyy[i] * pa_y[i] + 2.0 * ts_xxxy_xyyy[i] * fz_0;

        tk_xxxy_xyyz[i] = 2.0 * tk_xxx_xyz[i] * fe_0 + tk_xxx_xyyz[i] * pa_y[i] + 2.0 * ts_xxxy_xyyz[i] * fz_0;

        tk_xxxy_xyzz[i] = tk_xxx_xzz[i] * fe_0 + tk_xxx_xyzz[i] * pa_y[i] + 2.0 * ts_xxxy_xyzz[i] * fz_0;

        tk_xxxy_xzzz[i] = tk_xxx_xzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xzzz[i] * fz_0;

        tk_xxxy_yyyy[i] = 4.0 * tk_xxx_yyy[i] * fe_0 + tk_xxx_yyyy[i] * pa_y[i] + 2.0 * ts_xxxy_yyyy[i] * fz_0;

        tk_xxxy_yyyz[i] = 3.0 * tk_xxx_yyz[i] * fe_0 + tk_xxx_yyyz[i] * pa_y[i] + 2.0 * ts_xxxy_yyyz[i] * fz_0;

        tk_xxxy_yyzz[i] = 2.0 * tk_xxx_yzz[i] * fe_0 + tk_xxx_yyzz[i] * pa_y[i] + 2.0 * ts_xxxy_yyzz[i] * fz_0;

        tk_xxxy_yzzz[i] = tk_xxx_zzz[i] * fe_0 + tk_xxx_yzzz[i] * pa_y[i] + 2.0 * ts_xxxy_yzzz[i] * fz_0;

        tk_xxxy_zzzz[i] = tk_xxx_zzzz[i] * pa_y[i] + 2.0 * ts_xxxy_zzzz[i] * fz_0;
    }

    // Set up 30-45 components of targeted buffer : GG

    auto tk_xxxz_xxxx = pbuffer.data(idx_kin_gg + 30);

    auto tk_xxxz_xxxy = pbuffer.data(idx_kin_gg + 31);

    auto tk_xxxz_xxxz = pbuffer.data(idx_kin_gg + 32);

    auto tk_xxxz_xxyy = pbuffer.data(idx_kin_gg + 33);

    auto tk_xxxz_xxyz = pbuffer.data(idx_kin_gg + 34);

    auto tk_xxxz_xxzz = pbuffer.data(idx_kin_gg + 35);

    auto tk_xxxz_xyyy = pbuffer.data(idx_kin_gg + 36);

    auto tk_xxxz_xyyz = pbuffer.data(idx_kin_gg + 37);

    auto tk_xxxz_xyzz = pbuffer.data(idx_kin_gg + 38);

    auto tk_xxxz_xzzz = pbuffer.data(idx_kin_gg + 39);

    auto tk_xxxz_yyyy = pbuffer.data(idx_kin_gg + 40);

    auto tk_xxxz_yyyz = pbuffer.data(idx_kin_gg + 41);

    auto tk_xxxz_yyzz = pbuffer.data(idx_kin_gg + 42);

    auto tk_xxxz_yzzz = pbuffer.data(idx_kin_gg + 43);

    auto tk_xxxz_zzzz = pbuffer.data(idx_kin_gg + 44);

#pragma omp simd aligned(pa_z,             \
                             tk_xxx_xxx,   \
                             tk_xxx_xxxx,  \
                             tk_xxx_xxxy,  \
                             tk_xxx_xxxz,  \
                             tk_xxx_xxy,   \
                             tk_xxx_xxyy,  \
                             tk_xxx_xxyz,  \
                             tk_xxx_xxz,   \
                             tk_xxx_xxzz,  \
                             tk_xxx_xyy,   \
                             tk_xxx_xyyy,  \
                             tk_xxx_xyyz,  \
                             tk_xxx_xyz,   \
                             tk_xxx_xyzz,  \
                             tk_xxx_xzz,   \
                             tk_xxx_xzzz,  \
                             tk_xxx_yyy,   \
                             tk_xxx_yyyy,  \
                             tk_xxx_yyyz,  \
                             tk_xxx_yyz,   \
                             tk_xxx_yyzz,  \
                             tk_xxx_yzz,   \
                             tk_xxx_yzzz,  \
                             tk_xxx_zzz,   \
                             tk_xxx_zzzz,  \
                             tk_xxxz_xxxx, \
                             tk_xxxz_xxxy, \
                             tk_xxxz_xxxz, \
                             tk_xxxz_xxyy, \
                             tk_xxxz_xxyz, \
                             tk_xxxz_xxzz, \
                             tk_xxxz_xyyy, \
                             tk_xxxz_xyyz, \
                             tk_xxxz_xyzz, \
                             tk_xxxz_xzzz, \
                             tk_xxxz_yyyy, \
                             tk_xxxz_yyyz, \
                             tk_xxxz_yyzz, \
                             tk_xxxz_yzzz, \
                             tk_xxxz_zzzz, \
                             ts_xxxz_xxxx, \
                             ts_xxxz_xxxy, \
                             ts_xxxz_xxxz, \
                             ts_xxxz_xxyy, \
                             ts_xxxz_xxyz, \
                             ts_xxxz_xxzz, \
                             ts_xxxz_xyyy, \
                             ts_xxxz_xyyz, \
                             ts_xxxz_xyzz, \
                             ts_xxxz_xzzz, \
                             ts_xxxz_yyyy, \
                             ts_xxxz_yyyz, \
                             ts_xxxz_yyzz, \
                             ts_xxxz_yzzz, \
                             ts_xxxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxz_xxxx[i] = tk_xxx_xxxx[i] * pa_z[i] + 2.0 * ts_xxxz_xxxx[i] * fz_0;

        tk_xxxz_xxxy[i] = tk_xxx_xxxy[i] * pa_z[i] + 2.0 * ts_xxxz_xxxy[i] * fz_0;

        tk_xxxz_xxxz[i] = tk_xxx_xxx[i] * fe_0 + tk_xxx_xxxz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxz[i] * fz_0;

        tk_xxxz_xxyy[i] = tk_xxx_xxyy[i] * pa_z[i] + 2.0 * ts_xxxz_xxyy[i] * fz_0;

        tk_xxxz_xxyz[i] = tk_xxx_xxy[i] * fe_0 + tk_xxx_xxyz[i] * pa_z[i] + 2.0 * ts_xxxz_xxyz[i] * fz_0;

        tk_xxxz_xxzz[i] = 2.0 * tk_xxx_xxz[i] * fe_0 + tk_xxx_xxzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxzz[i] * fz_0;

        tk_xxxz_xyyy[i] = tk_xxx_xyyy[i] * pa_z[i] + 2.0 * ts_xxxz_xyyy[i] * fz_0;

        tk_xxxz_xyyz[i] = tk_xxx_xyy[i] * fe_0 + tk_xxx_xyyz[i] * pa_z[i] + 2.0 * ts_xxxz_xyyz[i] * fz_0;

        tk_xxxz_xyzz[i] = 2.0 * tk_xxx_xyz[i] * fe_0 + tk_xxx_xyzz[i] * pa_z[i] + 2.0 * ts_xxxz_xyzz[i] * fz_0;

        tk_xxxz_xzzz[i] = 3.0 * tk_xxx_xzz[i] * fe_0 + tk_xxx_xzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xzzz[i] * fz_0;

        tk_xxxz_yyyy[i] = tk_xxx_yyyy[i] * pa_z[i] + 2.0 * ts_xxxz_yyyy[i] * fz_0;

        tk_xxxz_yyyz[i] = tk_xxx_yyy[i] * fe_0 + tk_xxx_yyyz[i] * pa_z[i] + 2.0 * ts_xxxz_yyyz[i] * fz_0;

        tk_xxxz_yyzz[i] = 2.0 * tk_xxx_yyz[i] * fe_0 + tk_xxx_yyzz[i] * pa_z[i] + 2.0 * ts_xxxz_yyzz[i] * fz_0;

        tk_xxxz_yzzz[i] = 3.0 * tk_xxx_yzz[i] * fe_0 + tk_xxx_yzzz[i] * pa_z[i] + 2.0 * ts_xxxz_yzzz[i] * fz_0;

        tk_xxxz_zzzz[i] = 4.0 * tk_xxx_zzz[i] * fe_0 + tk_xxx_zzzz[i] * pa_z[i] + 2.0 * ts_xxxz_zzzz[i] * fz_0;
    }

    // Set up 45-60 components of targeted buffer : GG

    auto tk_xxyy_xxxx = pbuffer.data(idx_kin_gg + 45);

    auto tk_xxyy_xxxy = pbuffer.data(idx_kin_gg + 46);

    auto tk_xxyy_xxxz = pbuffer.data(idx_kin_gg + 47);

    auto tk_xxyy_xxyy = pbuffer.data(idx_kin_gg + 48);

    auto tk_xxyy_xxyz = pbuffer.data(idx_kin_gg + 49);

    auto tk_xxyy_xxzz = pbuffer.data(idx_kin_gg + 50);

    auto tk_xxyy_xyyy = pbuffer.data(idx_kin_gg + 51);

    auto tk_xxyy_xyyz = pbuffer.data(idx_kin_gg + 52);

    auto tk_xxyy_xyzz = pbuffer.data(idx_kin_gg + 53);

    auto tk_xxyy_xzzz = pbuffer.data(idx_kin_gg + 54);

    auto tk_xxyy_yyyy = pbuffer.data(idx_kin_gg + 55);

    auto tk_xxyy_yyyz = pbuffer.data(idx_kin_gg + 56);

    auto tk_xxyy_yyzz = pbuffer.data(idx_kin_gg + 57);

    auto tk_xxyy_yzzz = pbuffer.data(idx_kin_gg + 58);

    auto tk_xxyy_zzzz = pbuffer.data(idx_kin_gg + 59);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tk_xx_xxxx,   \
                             tk_xx_xxxz,   \
                             tk_xx_xxzz,   \
                             tk_xx_xzzz,   \
                             tk_xxy_xxxx,  \
                             tk_xxy_xxxz,  \
                             tk_xxy_xxzz,  \
                             tk_xxy_xzzz,  \
                             tk_xxyy_xxxx, \
                             tk_xxyy_xxxy, \
                             tk_xxyy_xxxz, \
                             tk_xxyy_xxyy, \
                             tk_xxyy_xxyz, \
                             tk_xxyy_xxzz, \
                             tk_xxyy_xyyy, \
                             tk_xxyy_xyyz, \
                             tk_xxyy_xyzz, \
                             tk_xxyy_xzzz, \
                             tk_xxyy_yyyy, \
                             tk_xxyy_yyyz, \
                             tk_xxyy_yyzz, \
                             tk_xxyy_yzzz, \
                             tk_xxyy_zzzz, \
                             tk_xyy_xxxy,  \
                             tk_xyy_xxy,   \
                             tk_xyy_xxyy,  \
                             tk_xyy_xxyz,  \
                             tk_xyy_xyy,   \
                             tk_xyy_xyyy,  \
                             tk_xyy_xyyz,  \
                             tk_xyy_xyz,   \
                             tk_xyy_xyzz,  \
                             tk_xyy_yyy,   \
                             tk_xyy_yyyy,  \
                             tk_xyy_yyyz,  \
                             tk_xyy_yyz,   \
                             tk_xyy_yyzz,  \
                             tk_xyy_yzz,   \
                             tk_xyy_yzzz,  \
                             tk_xyy_zzzz,  \
                             tk_yy_xxxy,   \
                             tk_yy_xxyy,   \
                             tk_yy_xxyz,   \
                             tk_yy_xyyy,   \
                             tk_yy_xyyz,   \
                             tk_yy_xyzz,   \
                             tk_yy_yyyy,   \
                             tk_yy_yyyz,   \
                             tk_yy_yyzz,   \
                             tk_yy_yzzz,   \
                             tk_yy_zzzz,   \
                             ts_xx_xxxx,   \
                             ts_xx_xxxz,   \
                             ts_xx_xxzz,   \
                             ts_xx_xzzz,   \
                             ts_xxyy_xxxx, \
                             ts_xxyy_xxxy, \
                             ts_xxyy_xxxz, \
                             ts_xxyy_xxyy, \
                             ts_xxyy_xxyz, \
                             ts_xxyy_xxzz, \
                             ts_xxyy_xyyy, \
                             ts_xxyy_xyyz, \
                             ts_xxyy_xyzz, \
                             ts_xxyy_xzzz, \
                             ts_xxyy_yyyy, \
                             ts_xxyy_yyyz, \
                             ts_xxyy_yyzz, \
                             ts_xxyy_yzzz, \
                             ts_xxyy_zzzz, \
                             ts_yy_xxxy,   \
                             ts_yy_xxyy,   \
                             ts_yy_xxyz,   \
                             ts_yy_xyyy,   \
                             ts_yy_xyyz,   \
                             ts_yy_xyzz,   \
                             ts_yy_yyyy,   \
                             ts_yy_yyyz,   \
                             ts_yy_yyzz,   \
                             ts_yy_yzzz,   \
                             ts_yy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyy_xxxx[i] = -2.0 * ts_xx_xxxx[i] * fbe_0 * fz_0 + tk_xx_xxxx[i] * fe_0 + tk_xxy_xxxx[i] * pa_y[i] + 2.0 * ts_xxyy_xxxx[i] * fz_0;

        tk_xxyy_xxxy[i] = -2.0 * ts_yy_xxxy[i] * fbe_0 * fz_0 + tk_yy_xxxy[i] * fe_0 + 3.0 * tk_xyy_xxy[i] * fe_0 + tk_xyy_xxxy[i] * pa_x[i] +
                          2.0 * ts_xxyy_xxxy[i] * fz_0;

        tk_xxyy_xxxz[i] = -2.0 * ts_xx_xxxz[i] * fbe_0 * fz_0 + tk_xx_xxxz[i] * fe_0 + tk_xxy_xxxz[i] * pa_y[i] + 2.0 * ts_xxyy_xxxz[i] * fz_0;

        tk_xxyy_xxyy[i] = -2.0 * ts_yy_xxyy[i] * fbe_0 * fz_0 + tk_yy_xxyy[i] * fe_0 + 2.0 * tk_xyy_xyy[i] * fe_0 + tk_xyy_xxyy[i] * pa_x[i] +
                          2.0 * ts_xxyy_xxyy[i] * fz_0;

        tk_xxyy_xxyz[i] = -2.0 * ts_yy_xxyz[i] * fbe_0 * fz_0 + tk_yy_xxyz[i] * fe_0 + 2.0 * tk_xyy_xyz[i] * fe_0 + tk_xyy_xxyz[i] * pa_x[i] +
                          2.0 * ts_xxyy_xxyz[i] * fz_0;

        tk_xxyy_xxzz[i] = -2.0 * ts_xx_xxzz[i] * fbe_0 * fz_0 + tk_xx_xxzz[i] * fe_0 + tk_xxy_xxzz[i] * pa_y[i] + 2.0 * ts_xxyy_xxzz[i] * fz_0;

        tk_xxyy_xyyy[i] = -2.0 * ts_yy_xyyy[i] * fbe_0 * fz_0 + tk_yy_xyyy[i] * fe_0 + tk_xyy_yyy[i] * fe_0 + tk_xyy_xyyy[i] * pa_x[i] +
                          2.0 * ts_xxyy_xyyy[i] * fz_0;

        tk_xxyy_xyyz[i] = -2.0 * ts_yy_xyyz[i] * fbe_0 * fz_0 + tk_yy_xyyz[i] * fe_0 + tk_xyy_yyz[i] * fe_0 + tk_xyy_xyyz[i] * pa_x[i] +
                          2.0 * ts_xxyy_xyyz[i] * fz_0;

        tk_xxyy_xyzz[i] = -2.0 * ts_yy_xyzz[i] * fbe_0 * fz_0 + tk_yy_xyzz[i] * fe_0 + tk_xyy_yzz[i] * fe_0 + tk_xyy_xyzz[i] * pa_x[i] +
                          2.0 * ts_xxyy_xyzz[i] * fz_0;

        tk_xxyy_xzzz[i] = -2.0 * ts_xx_xzzz[i] * fbe_0 * fz_0 + tk_xx_xzzz[i] * fe_0 + tk_xxy_xzzz[i] * pa_y[i] + 2.0 * ts_xxyy_xzzz[i] * fz_0;

        tk_xxyy_yyyy[i] = -2.0 * ts_yy_yyyy[i] * fbe_0 * fz_0 + tk_yy_yyyy[i] * fe_0 + tk_xyy_yyyy[i] * pa_x[i] + 2.0 * ts_xxyy_yyyy[i] * fz_0;

        tk_xxyy_yyyz[i] = -2.0 * ts_yy_yyyz[i] * fbe_0 * fz_0 + tk_yy_yyyz[i] * fe_0 + tk_xyy_yyyz[i] * pa_x[i] + 2.0 * ts_xxyy_yyyz[i] * fz_0;

        tk_xxyy_yyzz[i] = -2.0 * ts_yy_yyzz[i] * fbe_0 * fz_0 + tk_yy_yyzz[i] * fe_0 + tk_xyy_yyzz[i] * pa_x[i] + 2.0 * ts_xxyy_yyzz[i] * fz_0;

        tk_xxyy_yzzz[i] = -2.0 * ts_yy_yzzz[i] * fbe_0 * fz_0 + tk_yy_yzzz[i] * fe_0 + tk_xyy_yzzz[i] * pa_x[i] + 2.0 * ts_xxyy_yzzz[i] * fz_0;

        tk_xxyy_zzzz[i] = -2.0 * ts_yy_zzzz[i] * fbe_0 * fz_0 + tk_yy_zzzz[i] * fe_0 + tk_xyy_zzzz[i] * pa_x[i] + 2.0 * ts_xxyy_zzzz[i] * fz_0;
    }

    // Set up 60-75 components of targeted buffer : GG

    auto tk_xxyz_xxxx = pbuffer.data(idx_kin_gg + 60);

    auto tk_xxyz_xxxy = pbuffer.data(idx_kin_gg + 61);

    auto tk_xxyz_xxxz = pbuffer.data(idx_kin_gg + 62);

    auto tk_xxyz_xxyy = pbuffer.data(idx_kin_gg + 63);

    auto tk_xxyz_xxyz = pbuffer.data(idx_kin_gg + 64);

    auto tk_xxyz_xxzz = pbuffer.data(idx_kin_gg + 65);

    auto tk_xxyz_xyyy = pbuffer.data(idx_kin_gg + 66);

    auto tk_xxyz_xyyz = pbuffer.data(idx_kin_gg + 67);

    auto tk_xxyz_xyzz = pbuffer.data(idx_kin_gg + 68);

    auto tk_xxyz_xzzz = pbuffer.data(idx_kin_gg + 69);

    auto tk_xxyz_yyyy = pbuffer.data(idx_kin_gg + 70);

    auto tk_xxyz_yyyz = pbuffer.data(idx_kin_gg + 71);

    auto tk_xxyz_yyzz = pbuffer.data(idx_kin_gg + 72);

    auto tk_xxyz_yzzz = pbuffer.data(idx_kin_gg + 73);

    auto tk_xxyz_zzzz = pbuffer.data(idx_kin_gg + 74);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tk_xxy_xxxy,  \
                             tk_xxy_xxyy,  \
                             tk_xxy_xyyy,  \
                             tk_xxy_yyyy,  \
                             tk_xxyz_xxxx, \
                             tk_xxyz_xxxy, \
                             tk_xxyz_xxxz, \
                             tk_xxyz_xxyy, \
                             tk_xxyz_xxyz, \
                             tk_xxyz_xxzz, \
                             tk_xxyz_xyyy, \
                             tk_xxyz_xyyz, \
                             tk_xxyz_xyzz, \
                             tk_xxyz_xzzz, \
                             tk_xxyz_yyyy, \
                             tk_xxyz_yyyz, \
                             tk_xxyz_yyzz, \
                             tk_xxyz_yzzz, \
                             tk_xxyz_zzzz, \
                             tk_xxz_xxxx,  \
                             tk_xxz_xxxz,  \
                             tk_xxz_xxyz,  \
                             tk_xxz_xxz,   \
                             tk_xxz_xxzz,  \
                             tk_xxz_xyyz,  \
                             tk_xxz_xyz,   \
                             tk_xxz_xyzz,  \
                             tk_xxz_xzz,   \
                             tk_xxz_xzzz,  \
                             tk_xxz_yyyz,  \
                             tk_xxz_yyz,   \
                             tk_xxz_yyzz,  \
                             tk_xxz_yzz,   \
                             tk_xxz_yzzz,  \
                             tk_xxz_zzz,   \
                             tk_xxz_zzzz,  \
                             ts_xxyz_xxxx, \
                             ts_xxyz_xxxy, \
                             ts_xxyz_xxxz, \
                             ts_xxyz_xxyy, \
                             ts_xxyz_xxyz, \
                             ts_xxyz_xxzz, \
                             ts_xxyz_xyyy, \
                             ts_xxyz_xyyz, \
                             ts_xxyz_xyzz, \
                             ts_xxyz_xzzz, \
                             ts_xxyz_yyyy, \
                             ts_xxyz_yyyz, \
                             ts_xxyz_yyzz, \
                             ts_xxyz_yzzz, \
                             ts_xxyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyz_xxxx[i] = tk_xxz_xxxx[i] * pa_y[i] + 2.0 * ts_xxyz_xxxx[i] * fz_0;

        tk_xxyz_xxxy[i] = tk_xxy_xxxy[i] * pa_z[i] + 2.0 * ts_xxyz_xxxy[i] * fz_0;

        tk_xxyz_xxxz[i] = tk_xxz_xxxz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxz[i] * fz_0;

        tk_xxyz_xxyy[i] = tk_xxy_xxyy[i] * pa_z[i] + 2.0 * ts_xxyz_xxyy[i] * fz_0;

        tk_xxyz_xxyz[i] = tk_xxz_xxz[i] * fe_0 + tk_xxz_xxyz[i] * pa_y[i] + 2.0 * ts_xxyz_xxyz[i] * fz_0;

        tk_xxyz_xxzz[i] = tk_xxz_xxzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxzz[i] * fz_0;

        tk_xxyz_xyyy[i] = tk_xxy_xyyy[i] * pa_z[i] + 2.0 * ts_xxyz_xyyy[i] * fz_0;

        tk_xxyz_xyyz[i] = 2.0 * tk_xxz_xyz[i] * fe_0 + tk_xxz_xyyz[i] * pa_y[i] + 2.0 * ts_xxyz_xyyz[i] * fz_0;

        tk_xxyz_xyzz[i] = tk_xxz_xzz[i] * fe_0 + tk_xxz_xyzz[i] * pa_y[i] + 2.0 * ts_xxyz_xyzz[i] * fz_0;

        tk_xxyz_xzzz[i] = tk_xxz_xzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xzzz[i] * fz_0;

        tk_xxyz_yyyy[i] = tk_xxy_yyyy[i] * pa_z[i] + 2.0 * ts_xxyz_yyyy[i] * fz_0;

        tk_xxyz_yyyz[i] = 3.0 * tk_xxz_yyz[i] * fe_0 + tk_xxz_yyyz[i] * pa_y[i] + 2.0 * ts_xxyz_yyyz[i] * fz_0;

        tk_xxyz_yyzz[i] = 2.0 * tk_xxz_yzz[i] * fe_0 + tk_xxz_yyzz[i] * pa_y[i] + 2.0 * ts_xxyz_yyzz[i] * fz_0;

        tk_xxyz_yzzz[i] = tk_xxz_zzz[i] * fe_0 + tk_xxz_yzzz[i] * pa_y[i] + 2.0 * ts_xxyz_yzzz[i] * fz_0;

        tk_xxyz_zzzz[i] = tk_xxz_zzzz[i] * pa_y[i] + 2.0 * ts_xxyz_zzzz[i] * fz_0;
    }

    // Set up 75-90 components of targeted buffer : GG

    auto tk_xxzz_xxxx = pbuffer.data(idx_kin_gg + 75);

    auto tk_xxzz_xxxy = pbuffer.data(idx_kin_gg + 76);

    auto tk_xxzz_xxxz = pbuffer.data(idx_kin_gg + 77);

    auto tk_xxzz_xxyy = pbuffer.data(idx_kin_gg + 78);

    auto tk_xxzz_xxyz = pbuffer.data(idx_kin_gg + 79);

    auto tk_xxzz_xxzz = pbuffer.data(idx_kin_gg + 80);

    auto tk_xxzz_xyyy = pbuffer.data(idx_kin_gg + 81);

    auto tk_xxzz_xyyz = pbuffer.data(idx_kin_gg + 82);

    auto tk_xxzz_xyzz = pbuffer.data(idx_kin_gg + 83);

    auto tk_xxzz_xzzz = pbuffer.data(idx_kin_gg + 84);

    auto tk_xxzz_yyyy = pbuffer.data(idx_kin_gg + 85);

    auto tk_xxzz_yyyz = pbuffer.data(idx_kin_gg + 86);

    auto tk_xxzz_yyzz = pbuffer.data(idx_kin_gg + 87);

    auto tk_xxzz_yzzz = pbuffer.data(idx_kin_gg + 88);

    auto tk_xxzz_zzzz = pbuffer.data(idx_kin_gg + 89);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tk_xx_xxxx,   \
                             tk_xx_xxxy,   \
                             tk_xx_xxyy,   \
                             tk_xx_xyyy,   \
                             tk_xxz_xxxx,  \
                             tk_xxz_xxxy,  \
                             tk_xxz_xxyy,  \
                             tk_xxz_xyyy,  \
                             tk_xxzz_xxxx, \
                             tk_xxzz_xxxy, \
                             tk_xxzz_xxxz, \
                             tk_xxzz_xxyy, \
                             tk_xxzz_xxyz, \
                             tk_xxzz_xxzz, \
                             tk_xxzz_xyyy, \
                             tk_xxzz_xyyz, \
                             tk_xxzz_xyzz, \
                             tk_xxzz_xzzz, \
                             tk_xxzz_yyyy, \
                             tk_xxzz_yyyz, \
                             tk_xxzz_yyzz, \
                             tk_xxzz_yzzz, \
                             tk_xxzz_zzzz, \
                             tk_xzz_xxxz,  \
                             tk_xzz_xxyz,  \
                             tk_xzz_xxz,   \
                             tk_xzz_xxzz,  \
                             tk_xzz_xyyz,  \
                             tk_xzz_xyz,   \
                             tk_xzz_xyzz,  \
                             tk_xzz_xzz,   \
                             tk_xzz_xzzz,  \
                             tk_xzz_yyyy,  \
                             tk_xzz_yyyz,  \
                             tk_xzz_yyz,   \
                             tk_xzz_yyzz,  \
                             tk_xzz_yzz,   \
                             tk_xzz_yzzz,  \
                             tk_xzz_zzz,   \
                             tk_xzz_zzzz,  \
                             tk_zz_xxxz,   \
                             tk_zz_xxyz,   \
                             tk_zz_xxzz,   \
                             tk_zz_xyyz,   \
                             tk_zz_xyzz,   \
                             tk_zz_xzzz,   \
                             tk_zz_yyyy,   \
                             tk_zz_yyyz,   \
                             tk_zz_yyzz,   \
                             tk_zz_yzzz,   \
                             tk_zz_zzzz,   \
                             ts_xx_xxxx,   \
                             ts_xx_xxxy,   \
                             ts_xx_xxyy,   \
                             ts_xx_xyyy,   \
                             ts_xxzz_xxxx, \
                             ts_xxzz_xxxy, \
                             ts_xxzz_xxxz, \
                             ts_xxzz_xxyy, \
                             ts_xxzz_xxyz, \
                             ts_xxzz_xxzz, \
                             ts_xxzz_xyyy, \
                             ts_xxzz_xyyz, \
                             ts_xxzz_xyzz, \
                             ts_xxzz_xzzz, \
                             ts_xxzz_yyyy, \
                             ts_xxzz_yyyz, \
                             ts_xxzz_yyzz, \
                             ts_xxzz_yzzz, \
                             ts_xxzz_zzzz, \
                             ts_zz_xxxz,   \
                             ts_zz_xxyz,   \
                             ts_zz_xxzz,   \
                             ts_zz_xyyz,   \
                             ts_zz_xyzz,   \
                             ts_zz_xzzz,   \
                             ts_zz_yyyy,   \
                             ts_zz_yyyz,   \
                             ts_zz_yyzz,   \
                             ts_zz_yzzz,   \
                             ts_zz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzz_xxxx[i] = -2.0 * ts_xx_xxxx[i] * fbe_0 * fz_0 + tk_xx_xxxx[i] * fe_0 + tk_xxz_xxxx[i] * pa_z[i] + 2.0 * ts_xxzz_xxxx[i] * fz_0;

        tk_xxzz_xxxy[i] = -2.0 * ts_xx_xxxy[i] * fbe_0 * fz_0 + tk_xx_xxxy[i] * fe_0 + tk_xxz_xxxy[i] * pa_z[i] + 2.0 * ts_xxzz_xxxy[i] * fz_0;

        tk_xxzz_xxxz[i] = -2.0 * ts_zz_xxxz[i] * fbe_0 * fz_0 + tk_zz_xxxz[i] * fe_0 + 3.0 * tk_xzz_xxz[i] * fe_0 + tk_xzz_xxxz[i] * pa_x[i] +
                          2.0 * ts_xxzz_xxxz[i] * fz_0;

        tk_xxzz_xxyy[i] = -2.0 * ts_xx_xxyy[i] * fbe_0 * fz_0 + tk_xx_xxyy[i] * fe_0 + tk_xxz_xxyy[i] * pa_z[i] + 2.0 * ts_xxzz_xxyy[i] * fz_0;

        tk_xxzz_xxyz[i] = -2.0 * ts_zz_xxyz[i] * fbe_0 * fz_0 + tk_zz_xxyz[i] * fe_0 + 2.0 * tk_xzz_xyz[i] * fe_0 + tk_xzz_xxyz[i] * pa_x[i] +
                          2.0 * ts_xxzz_xxyz[i] * fz_0;

        tk_xxzz_xxzz[i] = -2.0 * ts_zz_xxzz[i] * fbe_0 * fz_0 + tk_zz_xxzz[i] * fe_0 + 2.0 * tk_xzz_xzz[i] * fe_0 + tk_xzz_xxzz[i] * pa_x[i] +
                          2.0 * ts_xxzz_xxzz[i] * fz_0;

        tk_xxzz_xyyy[i] = -2.0 * ts_xx_xyyy[i] * fbe_0 * fz_0 + tk_xx_xyyy[i] * fe_0 + tk_xxz_xyyy[i] * pa_z[i] + 2.0 * ts_xxzz_xyyy[i] * fz_0;

        tk_xxzz_xyyz[i] = -2.0 * ts_zz_xyyz[i] * fbe_0 * fz_0 + tk_zz_xyyz[i] * fe_0 + tk_xzz_yyz[i] * fe_0 + tk_xzz_xyyz[i] * pa_x[i] +
                          2.0 * ts_xxzz_xyyz[i] * fz_0;

        tk_xxzz_xyzz[i] = -2.0 * ts_zz_xyzz[i] * fbe_0 * fz_0 + tk_zz_xyzz[i] * fe_0 + tk_xzz_yzz[i] * fe_0 + tk_xzz_xyzz[i] * pa_x[i] +
                          2.0 * ts_xxzz_xyzz[i] * fz_0;

        tk_xxzz_xzzz[i] = -2.0 * ts_zz_xzzz[i] * fbe_0 * fz_0 + tk_zz_xzzz[i] * fe_0 + tk_xzz_zzz[i] * fe_0 + tk_xzz_xzzz[i] * pa_x[i] +
                          2.0 * ts_xxzz_xzzz[i] * fz_0;

        tk_xxzz_yyyy[i] = -2.0 * ts_zz_yyyy[i] * fbe_0 * fz_0 + tk_zz_yyyy[i] * fe_0 + tk_xzz_yyyy[i] * pa_x[i] + 2.0 * ts_xxzz_yyyy[i] * fz_0;

        tk_xxzz_yyyz[i] = -2.0 * ts_zz_yyyz[i] * fbe_0 * fz_0 + tk_zz_yyyz[i] * fe_0 + tk_xzz_yyyz[i] * pa_x[i] + 2.0 * ts_xxzz_yyyz[i] * fz_0;

        tk_xxzz_yyzz[i] = -2.0 * ts_zz_yyzz[i] * fbe_0 * fz_0 + tk_zz_yyzz[i] * fe_0 + tk_xzz_yyzz[i] * pa_x[i] + 2.0 * ts_xxzz_yyzz[i] * fz_0;

        tk_xxzz_yzzz[i] = -2.0 * ts_zz_yzzz[i] * fbe_0 * fz_0 + tk_zz_yzzz[i] * fe_0 + tk_xzz_yzzz[i] * pa_x[i] + 2.0 * ts_xxzz_yzzz[i] * fz_0;

        tk_xxzz_zzzz[i] = -2.0 * ts_zz_zzzz[i] * fbe_0 * fz_0 + tk_zz_zzzz[i] * fe_0 + tk_xzz_zzzz[i] * pa_x[i] + 2.0 * ts_xxzz_zzzz[i] * fz_0;
    }

    // Set up 90-105 components of targeted buffer : GG

    auto tk_xyyy_xxxx = pbuffer.data(idx_kin_gg + 90);

    auto tk_xyyy_xxxy = pbuffer.data(idx_kin_gg + 91);

    auto tk_xyyy_xxxz = pbuffer.data(idx_kin_gg + 92);

    auto tk_xyyy_xxyy = pbuffer.data(idx_kin_gg + 93);

    auto tk_xyyy_xxyz = pbuffer.data(idx_kin_gg + 94);

    auto tk_xyyy_xxzz = pbuffer.data(idx_kin_gg + 95);

    auto tk_xyyy_xyyy = pbuffer.data(idx_kin_gg + 96);

    auto tk_xyyy_xyyz = pbuffer.data(idx_kin_gg + 97);

    auto tk_xyyy_xyzz = pbuffer.data(idx_kin_gg + 98);

    auto tk_xyyy_xzzz = pbuffer.data(idx_kin_gg + 99);

    auto tk_xyyy_yyyy = pbuffer.data(idx_kin_gg + 100);

    auto tk_xyyy_yyyz = pbuffer.data(idx_kin_gg + 101);

    auto tk_xyyy_yyzz = pbuffer.data(idx_kin_gg + 102);

    auto tk_xyyy_yzzz = pbuffer.data(idx_kin_gg + 103);

    auto tk_xyyy_zzzz = pbuffer.data(idx_kin_gg + 104);

#pragma omp simd aligned(pa_x,             \
                             tk_xyyy_xxxx, \
                             tk_xyyy_xxxy, \
                             tk_xyyy_xxxz, \
                             tk_xyyy_xxyy, \
                             tk_xyyy_xxyz, \
                             tk_xyyy_xxzz, \
                             tk_xyyy_xyyy, \
                             tk_xyyy_xyyz, \
                             tk_xyyy_xyzz, \
                             tk_xyyy_xzzz, \
                             tk_xyyy_yyyy, \
                             tk_xyyy_yyyz, \
                             tk_xyyy_yyzz, \
                             tk_xyyy_yzzz, \
                             tk_xyyy_zzzz, \
                             tk_yyy_xxx,   \
                             tk_yyy_xxxx,  \
                             tk_yyy_xxxy,  \
                             tk_yyy_xxxz,  \
                             tk_yyy_xxy,   \
                             tk_yyy_xxyy,  \
                             tk_yyy_xxyz,  \
                             tk_yyy_xxz,   \
                             tk_yyy_xxzz,  \
                             tk_yyy_xyy,   \
                             tk_yyy_xyyy,  \
                             tk_yyy_xyyz,  \
                             tk_yyy_xyz,   \
                             tk_yyy_xyzz,  \
                             tk_yyy_xzz,   \
                             tk_yyy_xzzz,  \
                             tk_yyy_yyy,   \
                             tk_yyy_yyyy,  \
                             tk_yyy_yyyz,  \
                             tk_yyy_yyz,   \
                             tk_yyy_yyzz,  \
                             tk_yyy_yzz,   \
                             tk_yyy_yzzz,  \
                             tk_yyy_zzz,   \
                             tk_yyy_zzzz,  \
                             ts_xyyy_xxxx, \
                             ts_xyyy_xxxy, \
                             ts_xyyy_xxxz, \
                             ts_xyyy_xxyy, \
                             ts_xyyy_xxyz, \
                             ts_xyyy_xxzz, \
                             ts_xyyy_xyyy, \
                             ts_xyyy_xyyz, \
                             ts_xyyy_xyzz, \
                             ts_xyyy_xzzz, \
                             ts_xyyy_yyyy, \
                             ts_xyyy_yyyz, \
                             ts_xyyy_yyzz, \
                             ts_xyyy_yzzz, \
                             ts_xyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyy_xxxx[i] = 4.0 * tk_yyy_xxx[i] * fe_0 + tk_yyy_xxxx[i] * pa_x[i] + 2.0 * ts_xyyy_xxxx[i] * fz_0;

        tk_xyyy_xxxy[i] = 3.0 * tk_yyy_xxy[i] * fe_0 + tk_yyy_xxxy[i] * pa_x[i] + 2.0 * ts_xyyy_xxxy[i] * fz_0;

        tk_xyyy_xxxz[i] = 3.0 * tk_yyy_xxz[i] * fe_0 + tk_yyy_xxxz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxz[i] * fz_0;

        tk_xyyy_xxyy[i] = 2.0 * tk_yyy_xyy[i] * fe_0 + tk_yyy_xxyy[i] * pa_x[i] + 2.0 * ts_xyyy_xxyy[i] * fz_0;

        tk_xyyy_xxyz[i] = 2.0 * tk_yyy_xyz[i] * fe_0 + tk_yyy_xxyz[i] * pa_x[i] + 2.0 * ts_xyyy_xxyz[i] * fz_0;

        tk_xyyy_xxzz[i] = 2.0 * tk_yyy_xzz[i] * fe_0 + tk_yyy_xxzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxzz[i] * fz_0;

        tk_xyyy_xyyy[i] = tk_yyy_yyy[i] * fe_0 + tk_yyy_xyyy[i] * pa_x[i] + 2.0 * ts_xyyy_xyyy[i] * fz_0;

        tk_xyyy_xyyz[i] = tk_yyy_yyz[i] * fe_0 + tk_yyy_xyyz[i] * pa_x[i] + 2.0 * ts_xyyy_xyyz[i] * fz_0;

        tk_xyyy_xyzz[i] = tk_yyy_yzz[i] * fe_0 + tk_yyy_xyzz[i] * pa_x[i] + 2.0 * ts_xyyy_xyzz[i] * fz_0;

        tk_xyyy_xzzz[i] = tk_yyy_zzz[i] * fe_0 + tk_yyy_xzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xzzz[i] * fz_0;

        tk_xyyy_yyyy[i] = tk_yyy_yyyy[i] * pa_x[i] + 2.0 * ts_xyyy_yyyy[i] * fz_0;

        tk_xyyy_yyyz[i] = tk_yyy_yyyz[i] * pa_x[i] + 2.0 * ts_xyyy_yyyz[i] * fz_0;

        tk_xyyy_yyzz[i] = tk_yyy_yyzz[i] * pa_x[i] + 2.0 * ts_xyyy_yyzz[i] * fz_0;

        tk_xyyy_yzzz[i] = tk_yyy_yzzz[i] * pa_x[i] + 2.0 * ts_xyyy_yzzz[i] * fz_0;

        tk_xyyy_zzzz[i] = tk_yyy_zzzz[i] * pa_x[i] + 2.0 * ts_xyyy_zzzz[i] * fz_0;
    }

    // Set up 105-120 components of targeted buffer : GG

    auto tk_xyyz_xxxx = pbuffer.data(idx_kin_gg + 105);

    auto tk_xyyz_xxxy = pbuffer.data(idx_kin_gg + 106);

    auto tk_xyyz_xxxz = pbuffer.data(idx_kin_gg + 107);

    auto tk_xyyz_xxyy = pbuffer.data(idx_kin_gg + 108);

    auto tk_xyyz_xxyz = pbuffer.data(idx_kin_gg + 109);

    auto tk_xyyz_xxzz = pbuffer.data(idx_kin_gg + 110);

    auto tk_xyyz_xyyy = pbuffer.data(idx_kin_gg + 111);

    auto tk_xyyz_xyyz = pbuffer.data(idx_kin_gg + 112);

    auto tk_xyyz_xyzz = pbuffer.data(idx_kin_gg + 113);

    auto tk_xyyz_xzzz = pbuffer.data(idx_kin_gg + 114);

    auto tk_xyyz_yyyy = pbuffer.data(idx_kin_gg + 115);

    auto tk_xyyz_yyyz = pbuffer.data(idx_kin_gg + 116);

    auto tk_xyyz_yyzz = pbuffer.data(idx_kin_gg + 117);

    auto tk_xyyz_yzzz = pbuffer.data(idx_kin_gg + 118);

    auto tk_xyyz_zzzz = pbuffer.data(idx_kin_gg + 119);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tk_xyy_xxxx,  \
                             tk_xyy_xxxy,  \
                             tk_xyy_xxyy,  \
                             tk_xyy_xyyy,  \
                             tk_xyyz_xxxx, \
                             tk_xyyz_xxxy, \
                             tk_xyyz_xxxz, \
                             tk_xyyz_xxyy, \
                             tk_xyyz_xxyz, \
                             tk_xyyz_xxzz, \
                             tk_xyyz_xyyy, \
                             tk_xyyz_xyyz, \
                             tk_xyyz_xyzz, \
                             tk_xyyz_xzzz, \
                             tk_xyyz_yyyy, \
                             tk_xyyz_yyyz, \
                             tk_xyyz_yyzz, \
                             tk_xyyz_yzzz, \
                             tk_xyyz_zzzz, \
                             tk_yyz_xxxz,  \
                             tk_yyz_xxyz,  \
                             tk_yyz_xxz,   \
                             tk_yyz_xxzz,  \
                             tk_yyz_xyyz,  \
                             tk_yyz_xyz,   \
                             tk_yyz_xyzz,  \
                             tk_yyz_xzz,   \
                             tk_yyz_xzzz,  \
                             tk_yyz_yyyy,  \
                             tk_yyz_yyyz,  \
                             tk_yyz_yyz,   \
                             tk_yyz_yyzz,  \
                             tk_yyz_yzz,   \
                             tk_yyz_yzzz,  \
                             tk_yyz_zzz,   \
                             tk_yyz_zzzz,  \
                             ts_xyyz_xxxx, \
                             ts_xyyz_xxxy, \
                             ts_xyyz_xxxz, \
                             ts_xyyz_xxyy, \
                             ts_xyyz_xxyz, \
                             ts_xyyz_xxzz, \
                             ts_xyyz_xyyy, \
                             ts_xyyz_xyyz, \
                             ts_xyyz_xyzz, \
                             ts_xyyz_xzzz, \
                             ts_xyyz_yyyy, \
                             ts_xyyz_yyyz, \
                             ts_xyyz_yyzz, \
                             ts_xyyz_yzzz, \
                             ts_xyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyz_xxxx[i] = tk_xyy_xxxx[i] * pa_z[i] + 2.0 * ts_xyyz_xxxx[i] * fz_0;

        tk_xyyz_xxxy[i] = tk_xyy_xxxy[i] * pa_z[i] + 2.0 * ts_xyyz_xxxy[i] * fz_0;

        tk_xyyz_xxxz[i] = 3.0 * tk_yyz_xxz[i] * fe_0 + tk_yyz_xxxz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxz[i] * fz_0;

        tk_xyyz_xxyy[i] = tk_xyy_xxyy[i] * pa_z[i] + 2.0 * ts_xyyz_xxyy[i] * fz_0;

        tk_xyyz_xxyz[i] = 2.0 * tk_yyz_xyz[i] * fe_0 + tk_yyz_xxyz[i] * pa_x[i] + 2.0 * ts_xyyz_xxyz[i] * fz_0;

        tk_xyyz_xxzz[i] = 2.0 * tk_yyz_xzz[i] * fe_0 + tk_yyz_xxzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxzz[i] * fz_0;

        tk_xyyz_xyyy[i] = tk_xyy_xyyy[i] * pa_z[i] + 2.0 * ts_xyyz_xyyy[i] * fz_0;

        tk_xyyz_xyyz[i] = tk_yyz_yyz[i] * fe_0 + tk_yyz_xyyz[i] * pa_x[i] + 2.0 * ts_xyyz_xyyz[i] * fz_0;

        tk_xyyz_xyzz[i] = tk_yyz_yzz[i] * fe_0 + tk_yyz_xyzz[i] * pa_x[i] + 2.0 * ts_xyyz_xyzz[i] * fz_0;

        tk_xyyz_xzzz[i] = tk_yyz_zzz[i] * fe_0 + tk_yyz_xzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xzzz[i] * fz_0;

        tk_xyyz_yyyy[i] = tk_yyz_yyyy[i] * pa_x[i] + 2.0 * ts_xyyz_yyyy[i] * fz_0;

        tk_xyyz_yyyz[i] = tk_yyz_yyyz[i] * pa_x[i] + 2.0 * ts_xyyz_yyyz[i] * fz_0;

        tk_xyyz_yyzz[i] = tk_yyz_yyzz[i] * pa_x[i] + 2.0 * ts_xyyz_yyzz[i] * fz_0;

        tk_xyyz_yzzz[i] = tk_yyz_yzzz[i] * pa_x[i] + 2.0 * ts_xyyz_yzzz[i] * fz_0;

        tk_xyyz_zzzz[i] = tk_yyz_zzzz[i] * pa_x[i] + 2.0 * ts_xyyz_zzzz[i] * fz_0;
    }

    // Set up 120-135 components of targeted buffer : GG

    auto tk_xyzz_xxxx = pbuffer.data(idx_kin_gg + 120);

    auto tk_xyzz_xxxy = pbuffer.data(idx_kin_gg + 121);

    auto tk_xyzz_xxxz = pbuffer.data(idx_kin_gg + 122);

    auto tk_xyzz_xxyy = pbuffer.data(idx_kin_gg + 123);

    auto tk_xyzz_xxyz = pbuffer.data(idx_kin_gg + 124);

    auto tk_xyzz_xxzz = pbuffer.data(idx_kin_gg + 125);

    auto tk_xyzz_xyyy = pbuffer.data(idx_kin_gg + 126);

    auto tk_xyzz_xyyz = pbuffer.data(idx_kin_gg + 127);

    auto tk_xyzz_xyzz = pbuffer.data(idx_kin_gg + 128);

    auto tk_xyzz_xzzz = pbuffer.data(idx_kin_gg + 129);

    auto tk_xyzz_yyyy = pbuffer.data(idx_kin_gg + 130);

    auto tk_xyzz_yyyz = pbuffer.data(idx_kin_gg + 131);

    auto tk_xyzz_yyzz = pbuffer.data(idx_kin_gg + 132);

    auto tk_xyzz_yzzz = pbuffer.data(idx_kin_gg + 133);

    auto tk_xyzz_zzzz = pbuffer.data(idx_kin_gg + 134);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tk_xyzz_xxxx, \
                             tk_xyzz_xxxy, \
                             tk_xyzz_xxxz, \
                             tk_xyzz_xxyy, \
                             tk_xyzz_xxyz, \
                             tk_xyzz_xxzz, \
                             tk_xyzz_xyyy, \
                             tk_xyzz_xyyz, \
                             tk_xyzz_xyzz, \
                             tk_xyzz_xzzz, \
                             tk_xyzz_yyyy, \
                             tk_xyzz_yyyz, \
                             tk_xyzz_yyzz, \
                             tk_xyzz_yzzz, \
                             tk_xyzz_zzzz, \
                             tk_xzz_xxxx,  \
                             tk_xzz_xxxz,  \
                             tk_xzz_xxzz,  \
                             tk_xzz_xzzz,  \
                             tk_yzz_xxxy,  \
                             tk_yzz_xxy,   \
                             tk_yzz_xxyy,  \
                             tk_yzz_xxyz,  \
                             tk_yzz_xyy,   \
                             tk_yzz_xyyy,  \
                             tk_yzz_xyyz,  \
                             tk_yzz_xyz,   \
                             tk_yzz_xyzz,  \
                             tk_yzz_yyy,   \
                             tk_yzz_yyyy,  \
                             tk_yzz_yyyz,  \
                             tk_yzz_yyz,   \
                             tk_yzz_yyzz,  \
                             tk_yzz_yzz,   \
                             tk_yzz_yzzz,  \
                             tk_yzz_zzzz,  \
                             ts_xyzz_xxxx, \
                             ts_xyzz_xxxy, \
                             ts_xyzz_xxxz, \
                             ts_xyzz_xxyy, \
                             ts_xyzz_xxyz, \
                             ts_xyzz_xxzz, \
                             ts_xyzz_xyyy, \
                             ts_xyzz_xyyz, \
                             ts_xyzz_xyzz, \
                             ts_xyzz_xzzz, \
                             ts_xyzz_yyyy, \
                             ts_xyzz_yyyz, \
                             ts_xyzz_yyzz, \
                             ts_xyzz_yzzz, \
                             ts_xyzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzz_xxxx[i] = tk_xzz_xxxx[i] * pa_y[i] + 2.0 * ts_xyzz_xxxx[i] * fz_0;

        tk_xyzz_xxxy[i] = 3.0 * tk_yzz_xxy[i] * fe_0 + tk_yzz_xxxy[i] * pa_x[i] + 2.0 * ts_xyzz_xxxy[i] * fz_0;

        tk_xyzz_xxxz[i] = tk_xzz_xxxz[i] * pa_y[i] + 2.0 * ts_xyzz_xxxz[i] * fz_0;

        tk_xyzz_xxyy[i] = 2.0 * tk_yzz_xyy[i] * fe_0 + tk_yzz_xxyy[i] * pa_x[i] + 2.0 * ts_xyzz_xxyy[i] * fz_0;

        tk_xyzz_xxyz[i] = 2.0 * tk_yzz_xyz[i] * fe_0 + tk_yzz_xxyz[i] * pa_x[i] + 2.0 * ts_xyzz_xxyz[i] * fz_0;

        tk_xyzz_xxzz[i] = tk_xzz_xxzz[i] * pa_y[i] + 2.0 * ts_xyzz_xxzz[i] * fz_0;

        tk_xyzz_xyyy[i] = tk_yzz_yyy[i] * fe_0 + tk_yzz_xyyy[i] * pa_x[i] + 2.0 * ts_xyzz_xyyy[i] * fz_0;

        tk_xyzz_xyyz[i] = tk_yzz_yyz[i] * fe_0 + tk_yzz_xyyz[i] * pa_x[i] + 2.0 * ts_xyzz_xyyz[i] * fz_0;

        tk_xyzz_xyzz[i] = tk_yzz_yzz[i] * fe_0 + tk_yzz_xyzz[i] * pa_x[i] + 2.0 * ts_xyzz_xyzz[i] * fz_0;

        tk_xyzz_xzzz[i] = tk_xzz_xzzz[i] * pa_y[i] + 2.0 * ts_xyzz_xzzz[i] * fz_0;

        tk_xyzz_yyyy[i] = tk_yzz_yyyy[i] * pa_x[i] + 2.0 * ts_xyzz_yyyy[i] * fz_0;

        tk_xyzz_yyyz[i] = tk_yzz_yyyz[i] * pa_x[i] + 2.0 * ts_xyzz_yyyz[i] * fz_0;

        tk_xyzz_yyzz[i] = tk_yzz_yyzz[i] * pa_x[i] + 2.0 * ts_xyzz_yyzz[i] * fz_0;

        tk_xyzz_yzzz[i] = tk_yzz_yzzz[i] * pa_x[i] + 2.0 * ts_xyzz_yzzz[i] * fz_0;

        tk_xyzz_zzzz[i] = tk_yzz_zzzz[i] * pa_x[i] + 2.0 * ts_xyzz_zzzz[i] * fz_0;
    }

    // Set up 135-150 components of targeted buffer : GG

    auto tk_xzzz_xxxx = pbuffer.data(idx_kin_gg + 135);

    auto tk_xzzz_xxxy = pbuffer.data(idx_kin_gg + 136);

    auto tk_xzzz_xxxz = pbuffer.data(idx_kin_gg + 137);

    auto tk_xzzz_xxyy = pbuffer.data(idx_kin_gg + 138);

    auto tk_xzzz_xxyz = pbuffer.data(idx_kin_gg + 139);

    auto tk_xzzz_xxzz = pbuffer.data(idx_kin_gg + 140);

    auto tk_xzzz_xyyy = pbuffer.data(idx_kin_gg + 141);

    auto tk_xzzz_xyyz = pbuffer.data(idx_kin_gg + 142);

    auto tk_xzzz_xyzz = pbuffer.data(idx_kin_gg + 143);

    auto tk_xzzz_xzzz = pbuffer.data(idx_kin_gg + 144);

    auto tk_xzzz_yyyy = pbuffer.data(idx_kin_gg + 145);

    auto tk_xzzz_yyyz = pbuffer.data(idx_kin_gg + 146);

    auto tk_xzzz_yyzz = pbuffer.data(idx_kin_gg + 147);

    auto tk_xzzz_yzzz = pbuffer.data(idx_kin_gg + 148);

    auto tk_xzzz_zzzz = pbuffer.data(idx_kin_gg + 149);

#pragma omp simd aligned(pa_x,             \
                             tk_xzzz_xxxx, \
                             tk_xzzz_xxxy, \
                             tk_xzzz_xxxz, \
                             tk_xzzz_xxyy, \
                             tk_xzzz_xxyz, \
                             tk_xzzz_xxzz, \
                             tk_xzzz_xyyy, \
                             tk_xzzz_xyyz, \
                             tk_xzzz_xyzz, \
                             tk_xzzz_xzzz, \
                             tk_xzzz_yyyy, \
                             tk_xzzz_yyyz, \
                             tk_xzzz_yyzz, \
                             tk_xzzz_yzzz, \
                             tk_xzzz_zzzz, \
                             tk_zzz_xxx,   \
                             tk_zzz_xxxx,  \
                             tk_zzz_xxxy,  \
                             tk_zzz_xxxz,  \
                             tk_zzz_xxy,   \
                             tk_zzz_xxyy,  \
                             tk_zzz_xxyz,  \
                             tk_zzz_xxz,   \
                             tk_zzz_xxzz,  \
                             tk_zzz_xyy,   \
                             tk_zzz_xyyy,  \
                             tk_zzz_xyyz,  \
                             tk_zzz_xyz,   \
                             tk_zzz_xyzz,  \
                             tk_zzz_xzz,   \
                             tk_zzz_xzzz,  \
                             tk_zzz_yyy,   \
                             tk_zzz_yyyy,  \
                             tk_zzz_yyyz,  \
                             tk_zzz_yyz,   \
                             tk_zzz_yyzz,  \
                             tk_zzz_yzz,   \
                             tk_zzz_yzzz,  \
                             tk_zzz_zzz,   \
                             tk_zzz_zzzz,  \
                             ts_xzzz_xxxx, \
                             ts_xzzz_xxxy, \
                             ts_xzzz_xxxz, \
                             ts_xzzz_xxyy, \
                             ts_xzzz_xxyz, \
                             ts_xzzz_xxzz, \
                             ts_xzzz_xyyy, \
                             ts_xzzz_xyyz, \
                             ts_xzzz_xyzz, \
                             ts_xzzz_xzzz, \
                             ts_xzzz_yyyy, \
                             ts_xzzz_yyyz, \
                             ts_xzzz_yyzz, \
                             ts_xzzz_yzzz, \
                             ts_xzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzz_xxxx[i] = 4.0 * tk_zzz_xxx[i] * fe_0 + tk_zzz_xxxx[i] * pa_x[i] + 2.0 * ts_xzzz_xxxx[i] * fz_0;

        tk_xzzz_xxxy[i] = 3.0 * tk_zzz_xxy[i] * fe_0 + tk_zzz_xxxy[i] * pa_x[i] + 2.0 * ts_xzzz_xxxy[i] * fz_0;

        tk_xzzz_xxxz[i] = 3.0 * tk_zzz_xxz[i] * fe_0 + tk_zzz_xxxz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxz[i] * fz_0;

        tk_xzzz_xxyy[i] = 2.0 * tk_zzz_xyy[i] * fe_0 + tk_zzz_xxyy[i] * pa_x[i] + 2.0 * ts_xzzz_xxyy[i] * fz_0;

        tk_xzzz_xxyz[i] = 2.0 * tk_zzz_xyz[i] * fe_0 + tk_zzz_xxyz[i] * pa_x[i] + 2.0 * ts_xzzz_xxyz[i] * fz_0;

        tk_xzzz_xxzz[i] = 2.0 * tk_zzz_xzz[i] * fe_0 + tk_zzz_xxzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxzz[i] * fz_0;

        tk_xzzz_xyyy[i] = tk_zzz_yyy[i] * fe_0 + tk_zzz_xyyy[i] * pa_x[i] + 2.0 * ts_xzzz_xyyy[i] * fz_0;

        tk_xzzz_xyyz[i] = tk_zzz_yyz[i] * fe_0 + tk_zzz_xyyz[i] * pa_x[i] + 2.0 * ts_xzzz_xyyz[i] * fz_0;

        tk_xzzz_xyzz[i] = tk_zzz_yzz[i] * fe_0 + tk_zzz_xyzz[i] * pa_x[i] + 2.0 * ts_xzzz_xyzz[i] * fz_0;

        tk_xzzz_xzzz[i] = tk_zzz_zzz[i] * fe_0 + tk_zzz_xzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xzzz[i] * fz_0;

        tk_xzzz_yyyy[i] = tk_zzz_yyyy[i] * pa_x[i] + 2.0 * ts_xzzz_yyyy[i] * fz_0;

        tk_xzzz_yyyz[i] = tk_zzz_yyyz[i] * pa_x[i] + 2.0 * ts_xzzz_yyyz[i] * fz_0;

        tk_xzzz_yyzz[i] = tk_zzz_yyzz[i] * pa_x[i] + 2.0 * ts_xzzz_yyzz[i] * fz_0;

        tk_xzzz_yzzz[i] = tk_zzz_yzzz[i] * pa_x[i] + 2.0 * ts_xzzz_yzzz[i] * fz_0;

        tk_xzzz_zzzz[i] = tk_zzz_zzzz[i] * pa_x[i] + 2.0 * ts_xzzz_zzzz[i] * fz_0;
    }

    // Set up 150-165 components of targeted buffer : GG

    auto tk_yyyy_xxxx = pbuffer.data(idx_kin_gg + 150);

    auto tk_yyyy_xxxy = pbuffer.data(idx_kin_gg + 151);

    auto tk_yyyy_xxxz = pbuffer.data(idx_kin_gg + 152);

    auto tk_yyyy_xxyy = pbuffer.data(idx_kin_gg + 153);

    auto tk_yyyy_xxyz = pbuffer.data(idx_kin_gg + 154);

    auto tk_yyyy_xxzz = pbuffer.data(idx_kin_gg + 155);

    auto tk_yyyy_xyyy = pbuffer.data(idx_kin_gg + 156);

    auto tk_yyyy_xyyz = pbuffer.data(idx_kin_gg + 157);

    auto tk_yyyy_xyzz = pbuffer.data(idx_kin_gg + 158);

    auto tk_yyyy_xzzz = pbuffer.data(idx_kin_gg + 159);

    auto tk_yyyy_yyyy = pbuffer.data(idx_kin_gg + 160);

    auto tk_yyyy_yyyz = pbuffer.data(idx_kin_gg + 161);

    auto tk_yyyy_yyzz = pbuffer.data(idx_kin_gg + 162);

    auto tk_yyyy_yzzz = pbuffer.data(idx_kin_gg + 163);

    auto tk_yyyy_zzzz = pbuffer.data(idx_kin_gg + 164);

#pragma omp simd aligned(pa_y,             \
                             tk_yy_xxxx,   \
                             tk_yy_xxxy,   \
                             tk_yy_xxxz,   \
                             tk_yy_xxyy,   \
                             tk_yy_xxyz,   \
                             tk_yy_xxzz,   \
                             tk_yy_xyyy,   \
                             tk_yy_xyyz,   \
                             tk_yy_xyzz,   \
                             tk_yy_xzzz,   \
                             tk_yy_yyyy,   \
                             tk_yy_yyyz,   \
                             tk_yy_yyzz,   \
                             tk_yy_yzzz,   \
                             tk_yy_zzzz,   \
                             tk_yyy_xxx,   \
                             tk_yyy_xxxx,  \
                             tk_yyy_xxxy,  \
                             tk_yyy_xxxz,  \
                             tk_yyy_xxy,   \
                             tk_yyy_xxyy,  \
                             tk_yyy_xxyz,  \
                             tk_yyy_xxz,   \
                             tk_yyy_xxzz,  \
                             tk_yyy_xyy,   \
                             tk_yyy_xyyy,  \
                             tk_yyy_xyyz,  \
                             tk_yyy_xyz,   \
                             tk_yyy_xyzz,  \
                             tk_yyy_xzz,   \
                             tk_yyy_xzzz,  \
                             tk_yyy_yyy,   \
                             tk_yyy_yyyy,  \
                             tk_yyy_yyyz,  \
                             tk_yyy_yyz,   \
                             tk_yyy_yyzz,  \
                             tk_yyy_yzz,   \
                             tk_yyy_yzzz,  \
                             tk_yyy_zzz,   \
                             tk_yyy_zzzz,  \
                             tk_yyyy_xxxx, \
                             tk_yyyy_xxxy, \
                             tk_yyyy_xxxz, \
                             tk_yyyy_xxyy, \
                             tk_yyyy_xxyz, \
                             tk_yyyy_xxzz, \
                             tk_yyyy_xyyy, \
                             tk_yyyy_xyyz, \
                             tk_yyyy_xyzz, \
                             tk_yyyy_xzzz, \
                             tk_yyyy_yyyy, \
                             tk_yyyy_yyyz, \
                             tk_yyyy_yyzz, \
                             tk_yyyy_yzzz, \
                             tk_yyyy_zzzz, \
                             ts_yy_xxxx,   \
                             ts_yy_xxxy,   \
                             ts_yy_xxxz,   \
                             ts_yy_xxyy,   \
                             ts_yy_xxyz,   \
                             ts_yy_xxzz,   \
                             ts_yy_xyyy,   \
                             ts_yy_xyyz,   \
                             ts_yy_xyzz,   \
                             ts_yy_xzzz,   \
                             ts_yy_yyyy,   \
                             ts_yy_yyyz,   \
                             ts_yy_yyzz,   \
                             ts_yy_yzzz,   \
                             ts_yy_zzzz,   \
                             ts_yyyy_xxxx, \
                             ts_yyyy_xxxy, \
                             ts_yyyy_xxxz, \
                             ts_yyyy_xxyy, \
                             ts_yyyy_xxyz, \
                             ts_yyyy_xxzz, \
                             ts_yyyy_xyyy, \
                             ts_yyyy_xyyz, \
                             ts_yyyy_xyzz, \
                             ts_yyyy_xzzz, \
                             ts_yyyy_yyyy, \
                             ts_yyyy_yyyz, \
                             ts_yyyy_yyzz, \
                             ts_yyyy_yzzz, \
                             ts_yyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyy_xxxx[i] = -6.0 * ts_yy_xxxx[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxx[i] * fe_0 + tk_yyy_xxxx[i] * pa_y[i] + 2.0 * ts_yyyy_xxxx[i] * fz_0;

        tk_yyyy_xxxy[i] = -6.0 * ts_yy_xxxy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxy[i] * fe_0 + tk_yyy_xxx[i] * fe_0 + tk_yyy_xxxy[i] * pa_y[i] +
                          2.0 * ts_yyyy_xxxy[i] * fz_0;

        tk_yyyy_xxxz[i] = -6.0 * ts_yy_xxxz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxz[i] * fe_0 + tk_yyy_xxxz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxz[i] * fz_0;

        tk_yyyy_xxyy[i] = -6.0 * ts_yy_xxyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyy[i] * fe_0 + 2.0 * tk_yyy_xxy[i] * fe_0 + tk_yyy_xxyy[i] * pa_y[i] +
                          2.0 * ts_yyyy_xxyy[i] * fz_0;

        tk_yyyy_xxyz[i] = -6.0 * ts_yy_xxyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyz[i] * fe_0 + tk_yyy_xxz[i] * fe_0 + tk_yyy_xxyz[i] * pa_y[i] +
                          2.0 * ts_yyyy_xxyz[i] * fz_0;

        tk_yyyy_xxzz[i] = -6.0 * ts_yy_xxzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxzz[i] * fe_0 + tk_yyy_xxzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxzz[i] * fz_0;

        tk_yyyy_xyyy[i] = -6.0 * ts_yy_xyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyy[i] * fe_0 + 3.0 * tk_yyy_xyy[i] * fe_0 + tk_yyy_xyyy[i] * pa_y[i] +
                          2.0 * ts_yyyy_xyyy[i] * fz_0;

        tk_yyyy_xyyz[i] = -6.0 * ts_yy_xyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyz[i] * fe_0 + 2.0 * tk_yyy_xyz[i] * fe_0 + tk_yyy_xyyz[i] * pa_y[i] +
                          2.0 * ts_yyyy_xyyz[i] * fz_0;

        tk_yyyy_xyzz[i] = -6.0 * ts_yy_xyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyzz[i] * fe_0 + tk_yyy_xzz[i] * fe_0 + tk_yyy_xyzz[i] * pa_y[i] +
                          2.0 * ts_yyyy_xyzz[i] * fz_0;

        tk_yyyy_xzzz[i] = -6.0 * ts_yy_xzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xzzz[i] * fe_0 + tk_yyy_xzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xzzz[i] * fz_0;

        tk_yyyy_yyyy[i] = -6.0 * ts_yy_yyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyy[i] * fe_0 + 4.0 * tk_yyy_yyy[i] * fe_0 + tk_yyy_yyyy[i] * pa_y[i] +
                          2.0 * ts_yyyy_yyyy[i] * fz_0;

        tk_yyyy_yyyz[i] = -6.0 * ts_yy_yyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyz[i] * fe_0 + 3.0 * tk_yyy_yyz[i] * fe_0 + tk_yyy_yyyz[i] * pa_y[i] +
                          2.0 * ts_yyyy_yyyz[i] * fz_0;

        tk_yyyy_yyzz[i] = -6.0 * ts_yy_yyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyzz[i] * fe_0 + 2.0 * tk_yyy_yzz[i] * fe_0 + tk_yyy_yyzz[i] * pa_y[i] +
                          2.0 * ts_yyyy_yyzz[i] * fz_0;

        tk_yyyy_yzzz[i] = -6.0 * ts_yy_yzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yzzz[i] * fe_0 + tk_yyy_zzz[i] * fe_0 + tk_yyy_yzzz[i] * pa_y[i] +
                          2.0 * ts_yyyy_yzzz[i] * fz_0;

        tk_yyyy_zzzz[i] = -6.0 * ts_yy_zzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_zzzz[i] * fe_0 + tk_yyy_zzzz[i] * pa_y[i] + 2.0 * ts_yyyy_zzzz[i] * fz_0;
    }

    // Set up 165-180 components of targeted buffer : GG

    auto tk_yyyz_xxxx = pbuffer.data(idx_kin_gg + 165);

    auto tk_yyyz_xxxy = pbuffer.data(idx_kin_gg + 166);

    auto tk_yyyz_xxxz = pbuffer.data(idx_kin_gg + 167);

    auto tk_yyyz_xxyy = pbuffer.data(idx_kin_gg + 168);

    auto tk_yyyz_xxyz = pbuffer.data(idx_kin_gg + 169);

    auto tk_yyyz_xxzz = pbuffer.data(idx_kin_gg + 170);

    auto tk_yyyz_xyyy = pbuffer.data(idx_kin_gg + 171);

    auto tk_yyyz_xyyz = pbuffer.data(idx_kin_gg + 172);

    auto tk_yyyz_xyzz = pbuffer.data(idx_kin_gg + 173);

    auto tk_yyyz_xzzz = pbuffer.data(idx_kin_gg + 174);

    auto tk_yyyz_yyyy = pbuffer.data(idx_kin_gg + 175);

    auto tk_yyyz_yyyz = pbuffer.data(idx_kin_gg + 176);

    auto tk_yyyz_yyzz = pbuffer.data(idx_kin_gg + 177);

    auto tk_yyyz_yzzz = pbuffer.data(idx_kin_gg + 178);

    auto tk_yyyz_zzzz = pbuffer.data(idx_kin_gg + 179);

#pragma omp simd aligned(pa_z,             \
                             tk_yyy_xxx,   \
                             tk_yyy_xxxx,  \
                             tk_yyy_xxxy,  \
                             tk_yyy_xxxz,  \
                             tk_yyy_xxy,   \
                             tk_yyy_xxyy,  \
                             tk_yyy_xxyz,  \
                             tk_yyy_xxz,   \
                             tk_yyy_xxzz,  \
                             tk_yyy_xyy,   \
                             tk_yyy_xyyy,  \
                             tk_yyy_xyyz,  \
                             tk_yyy_xyz,   \
                             tk_yyy_xyzz,  \
                             tk_yyy_xzz,   \
                             tk_yyy_xzzz,  \
                             tk_yyy_yyy,   \
                             tk_yyy_yyyy,  \
                             tk_yyy_yyyz,  \
                             tk_yyy_yyz,   \
                             tk_yyy_yyzz,  \
                             tk_yyy_yzz,   \
                             tk_yyy_yzzz,  \
                             tk_yyy_zzz,   \
                             tk_yyy_zzzz,  \
                             tk_yyyz_xxxx, \
                             tk_yyyz_xxxy, \
                             tk_yyyz_xxxz, \
                             tk_yyyz_xxyy, \
                             tk_yyyz_xxyz, \
                             tk_yyyz_xxzz, \
                             tk_yyyz_xyyy, \
                             tk_yyyz_xyyz, \
                             tk_yyyz_xyzz, \
                             tk_yyyz_xzzz, \
                             tk_yyyz_yyyy, \
                             tk_yyyz_yyyz, \
                             tk_yyyz_yyzz, \
                             tk_yyyz_yzzz, \
                             tk_yyyz_zzzz, \
                             ts_yyyz_xxxx, \
                             ts_yyyz_xxxy, \
                             ts_yyyz_xxxz, \
                             ts_yyyz_xxyy, \
                             ts_yyyz_xxyz, \
                             ts_yyyz_xxzz, \
                             ts_yyyz_xyyy, \
                             ts_yyyz_xyyz, \
                             ts_yyyz_xyzz, \
                             ts_yyyz_xzzz, \
                             ts_yyyz_yyyy, \
                             ts_yyyz_yyyz, \
                             ts_yyyz_yyzz, \
                             ts_yyyz_yzzz, \
                             ts_yyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyz_xxxx[i] = tk_yyy_xxxx[i] * pa_z[i] + 2.0 * ts_yyyz_xxxx[i] * fz_0;

        tk_yyyz_xxxy[i] = tk_yyy_xxxy[i] * pa_z[i] + 2.0 * ts_yyyz_xxxy[i] * fz_0;

        tk_yyyz_xxxz[i] = tk_yyy_xxx[i] * fe_0 + tk_yyy_xxxz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxz[i] * fz_0;

        tk_yyyz_xxyy[i] = tk_yyy_xxyy[i] * pa_z[i] + 2.0 * ts_yyyz_xxyy[i] * fz_0;

        tk_yyyz_xxyz[i] = tk_yyy_xxy[i] * fe_0 + tk_yyy_xxyz[i] * pa_z[i] + 2.0 * ts_yyyz_xxyz[i] * fz_0;

        tk_yyyz_xxzz[i] = 2.0 * tk_yyy_xxz[i] * fe_0 + tk_yyy_xxzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxzz[i] * fz_0;

        tk_yyyz_xyyy[i] = tk_yyy_xyyy[i] * pa_z[i] + 2.0 * ts_yyyz_xyyy[i] * fz_0;

        tk_yyyz_xyyz[i] = tk_yyy_xyy[i] * fe_0 + tk_yyy_xyyz[i] * pa_z[i] + 2.0 * ts_yyyz_xyyz[i] * fz_0;

        tk_yyyz_xyzz[i] = 2.0 * tk_yyy_xyz[i] * fe_0 + tk_yyy_xyzz[i] * pa_z[i] + 2.0 * ts_yyyz_xyzz[i] * fz_0;

        tk_yyyz_xzzz[i] = 3.0 * tk_yyy_xzz[i] * fe_0 + tk_yyy_xzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xzzz[i] * fz_0;

        tk_yyyz_yyyy[i] = tk_yyy_yyyy[i] * pa_z[i] + 2.0 * ts_yyyz_yyyy[i] * fz_0;

        tk_yyyz_yyyz[i] = tk_yyy_yyy[i] * fe_0 + tk_yyy_yyyz[i] * pa_z[i] + 2.0 * ts_yyyz_yyyz[i] * fz_0;

        tk_yyyz_yyzz[i] = 2.0 * tk_yyy_yyz[i] * fe_0 + tk_yyy_yyzz[i] * pa_z[i] + 2.0 * ts_yyyz_yyzz[i] * fz_0;

        tk_yyyz_yzzz[i] = 3.0 * tk_yyy_yzz[i] * fe_0 + tk_yyy_yzzz[i] * pa_z[i] + 2.0 * ts_yyyz_yzzz[i] * fz_0;

        tk_yyyz_zzzz[i] = 4.0 * tk_yyy_zzz[i] * fe_0 + tk_yyy_zzzz[i] * pa_z[i] + 2.0 * ts_yyyz_zzzz[i] * fz_0;
    }

    // Set up 180-195 components of targeted buffer : GG

    auto tk_yyzz_xxxx = pbuffer.data(idx_kin_gg + 180);

    auto tk_yyzz_xxxy = pbuffer.data(idx_kin_gg + 181);

    auto tk_yyzz_xxxz = pbuffer.data(idx_kin_gg + 182);

    auto tk_yyzz_xxyy = pbuffer.data(idx_kin_gg + 183);

    auto tk_yyzz_xxyz = pbuffer.data(idx_kin_gg + 184);

    auto tk_yyzz_xxzz = pbuffer.data(idx_kin_gg + 185);

    auto tk_yyzz_xyyy = pbuffer.data(idx_kin_gg + 186);

    auto tk_yyzz_xyyz = pbuffer.data(idx_kin_gg + 187);

    auto tk_yyzz_xyzz = pbuffer.data(idx_kin_gg + 188);

    auto tk_yyzz_xzzz = pbuffer.data(idx_kin_gg + 189);

    auto tk_yyzz_yyyy = pbuffer.data(idx_kin_gg + 190);

    auto tk_yyzz_yyyz = pbuffer.data(idx_kin_gg + 191);

    auto tk_yyzz_yyzz = pbuffer.data(idx_kin_gg + 192);

    auto tk_yyzz_yzzz = pbuffer.data(idx_kin_gg + 193);

    auto tk_yyzz_zzzz = pbuffer.data(idx_kin_gg + 194);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tk_yy_xxxy,   \
                             tk_yy_xxyy,   \
                             tk_yy_xyyy,   \
                             tk_yy_yyyy,   \
                             tk_yyz_xxxy,  \
                             tk_yyz_xxyy,  \
                             tk_yyz_xyyy,  \
                             tk_yyz_yyyy,  \
                             tk_yyzz_xxxx, \
                             tk_yyzz_xxxy, \
                             tk_yyzz_xxxz, \
                             tk_yyzz_xxyy, \
                             tk_yyzz_xxyz, \
                             tk_yyzz_xxzz, \
                             tk_yyzz_xyyy, \
                             tk_yyzz_xyyz, \
                             tk_yyzz_xyzz, \
                             tk_yyzz_xzzz, \
                             tk_yyzz_yyyy, \
                             tk_yyzz_yyyz, \
                             tk_yyzz_yyzz, \
                             tk_yyzz_yzzz, \
                             tk_yyzz_zzzz, \
                             tk_yzz_xxxx,  \
                             tk_yzz_xxxz,  \
                             tk_yzz_xxyz,  \
                             tk_yzz_xxz,   \
                             tk_yzz_xxzz,  \
                             tk_yzz_xyyz,  \
                             tk_yzz_xyz,   \
                             tk_yzz_xyzz,  \
                             tk_yzz_xzz,   \
                             tk_yzz_xzzz,  \
                             tk_yzz_yyyz,  \
                             tk_yzz_yyz,   \
                             tk_yzz_yyzz,  \
                             tk_yzz_yzz,   \
                             tk_yzz_yzzz,  \
                             tk_yzz_zzz,   \
                             tk_yzz_zzzz,  \
                             tk_zz_xxxx,   \
                             tk_zz_xxxz,   \
                             tk_zz_xxyz,   \
                             tk_zz_xxzz,   \
                             tk_zz_xyyz,   \
                             tk_zz_xyzz,   \
                             tk_zz_xzzz,   \
                             tk_zz_yyyz,   \
                             tk_zz_yyzz,   \
                             tk_zz_yzzz,   \
                             tk_zz_zzzz,   \
                             ts_yy_xxxy,   \
                             ts_yy_xxyy,   \
                             ts_yy_xyyy,   \
                             ts_yy_yyyy,   \
                             ts_yyzz_xxxx, \
                             ts_yyzz_xxxy, \
                             ts_yyzz_xxxz, \
                             ts_yyzz_xxyy, \
                             ts_yyzz_xxyz, \
                             ts_yyzz_xxzz, \
                             ts_yyzz_xyyy, \
                             ts_yyzz_xyyz, \
                             ts_yyzz_xyzz, \
                             ts_yyzz_xzzz, \
                             ts_yyzz_yyyy, \
                             ts_yyzz_yyyz, \
                             ts_yyzz_yyzz, \
                             ts_yyzz_yzzz, \
                             ts_yyzz_zzzz, \
                             ts_zz_xxxx,   \
                             ts_zz_xxxz,   \
                             ts_zz_xxyz,   \
                             ts_zz_xxzz,   \
                             ts_zz_xyyz,   \
                             ts_zz_xyzz,   \
                             ts_zz_xzzz,   \
                             ts_zz_yyyz,   \
                             ts_zz_yyzz,   \
                             ts_zz_yzzz,   \
                             ts_zz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzz_xxxx[i] = -2.0 * ts_zz_xxxx[i] * fbe_0 * fz_0 + tk_zz_xxxx[i] * fe_0 + tk_yzz_xxxx[i] * pa_y[i] + 2.0 * ts_yyzz_xxxx[i] * fz_0;

        tk_yyzz_xxxy[i] = -2.0 * ts_yy_xxxy[i] * fbe_0 * fz_0 + tk_yy_xxxy[i] * fe_0 + tk_yyz_xxxy[i] * pa_z[i] + 2.0 * ts_yyzz_xxxy[i] * fz_0;

        tk_yyzz_xxxz[i] = -2.0 * ts_zz_xxxz[i] * fbe_0 * fz_0 + tk_zz_xxxz[i] * fe_0 + tk_yzz_xxxz[i] * pa_y[i] + 2.0 * ts_yyzz_xxxz[i] * fz_0;

        tk_yyzz_xxyy[i] = -2.0 * ts_yy_xxyy[i] * fbe_0 * fz_0 + tk_yy_xxyy[i] * fe_0 + tk_yyz_xxyy[i] * pa_z[i] + 2.0 * ts_yyzz_xxyy[i] * fz_0;

        tk_yyzz_xxyz[i] = -2.0 * ts_zz_xxyz[i] * fbe_0 * fz_0 + tk_zz_xxyz[i] * fe_0 + tk_yzz_xxz[i] * fe_0 + tk_yzz_xxyz[i] * pa_y[i] +
                          2.0 * ts_yyzz_xxyz[i] * fz_0;

        tk_yyzz_xxzz[i] = -2.0 * ts_zz_xxzz[i] * fbe_0 * fz_0 + tk_zz_xxzz[i] * fe_0 + tk_yzz_xxzz[i] * pa_y[i] + 2.0 * ts_yyzz_xxzz[i] * fz_0;

        tk_yyzz_xyyy[i] = -2.0 * ts_yy_xyyy[i] * fbe_0 * fz_0 + tk_yy_xyyy[i] * fe_0 + tk_yyz_xyyy[i] * pa_z[i] + 2.0 * ts_yyzz_xyyy[i] * fz_0;

        tk_yyzz_xyyz[i] = -2.0 * ts_zz_xyyz[i] * fbe_0 * fz_0 + tk_zz_xyyz[i] * fe_0 + 2.0 * tk_yzz_xyz[i] * fe_0 + tk_yzz_xyyz[i] * pa_y[i] +
                          2.0 * ts_yyzz_xyyz[i] * fz_0;

        tk_yyzz_xyzz[i] = -2.0 * ts_zz_xyzz[i] * fbe_0 * fz_0 + tk_zz_xyzz[i] * fe_0 + tk_yzz_xzz[i] * fe_0 + tk_yzz_xyzz[i] * pa_y[i] +
                          2.0 * ts_yyzz_xyzz[i] * fz_0;

        tk_yyzz_xzzz[i] = -2.0 * ts_zz_xzzz[i] * fbe_0 * fz_0 + tk_zz_xzzz[i] * fe_0 + tk_yzz_xzzz[i] * pa_y[i] + 2.0 * ts_yyzz_xzzz[i] * fz_0;

        tk_yyzz_yyyy[i] = -2.0 * ts_yy_yyyy[i] * fbe_0 * fz_0 + tk_yy_yyyy[i] * fe_0 + tk_yyz_yyyy[i] * pa_z[i] + 2.0 * ts_yyzz_yyyy[i] * fz_0;

        tk_yyzz_yyyz[i] = -2.0 * ts_zz_yyyz[i] * fbe_0 * fz_0 + tk_zz_yyyz[i] * fe_0 + 3.0 * tk_yzz_yyz[i] * fe_0 + tk_yzz_yyyz[i] * pa_y[i] +
                          2.0 * ts_yyzz_yyyz[i] * fz_0;

        tk_yyzz_yyzz[i] = -2.0 * ts_zz_yyzz[i] * fbe_0 * fz_0 + tk_zz_yyzz[i] * fe_0 + 2.0 * tk_yzz_yzz[i] * fe_0 + tk_yzz_yyzz[i] * pa_y[i] +
                          2.0 * ts_yyzz_yyzz[i] * fz_0;

        tk_yyzz_yzzz[i] = -2.0 * ts_zz_yzzz[i] * fbe_0 * fz_0 + tk_zz_yzzz[i] * fe_0 + tk_yzz_zzz[i] * fe_0 + tk_yzz_yzzz[i] * pa_y[i] +
                          2.0 * ts_yyzz_yzzz[i] * fz_0;

        tk_yyzz_zzzz[i] = -2.0 * ts_zz_zzzz[i] * fbe_0 * fz_0 + tk_zz_zzzz[i] * fe_0 + tk_yzz_zzzz[i] * pa_y[i] + 2.0 * ts_yyzz_zzzz[i] * fz_0;
    }

    // Set up 195-210 components of targeted buffer : GG

    auto tk_yzzz_xxxx = pbuffer.data(idx_kin_gg + 195);

    auto tk_yzzz_xxxy = pbuffer.data(idx_kin_gg + 196);

    auto tk_yzzz_xxxz = pbuffer.data(idx_kin_gg + 197);

    auto tk_yzzz_xxyy = pbuffer.data(idx_kin_gg + 198);

    auto tk_yzzz_xxyz = pbuffer.data(idx_kin_gg + 199);

    auto tk_yzzz_xxzz = pbuffer.data(idx_kin_gg + 200);

    auto tk_yzzz_xyyy = pbuffer.data(idx_kin_gg + 201);

    auto tk_yzzz_xyyz = pbuffer.data(idx_kin_gg + 202);

    auto tk_yzzz_xyzz = pbuffer.data(idx_kin_gg + 203);

    auto tk_yzzz_xzzz = pbuffer.data(idx_kin_gg + 204);

    auto tk_yzzz_yyyy = pbuffer.data(idx_kin_gg + 205);

    auto tk_yzzz_yyyz = pbuffer.data(idx_kin_gg + 206);

    auto tk_yzzz_yyzz = pbuffer.data(idx_kin_gg + 207);

    auto tk_yzzz_yzzz = pbuffer.data(idx_kin_gg + 208);

    auto tk_yzzz_zzzz = pbuffer.data(idx_kin_gg + 209);

#pragma omp simd aligned(pa_y,             \
                             tk_yzzz_xxxx, \
                             tk_yzzz_xxxy, \
                             tk_yzzz_xxxz, \
                             tk_yzzz_xxyy, \
                             tk_yzzz_xxyz, \
                             tk_yzzz_xxzz, \
                             tk_yzzz_xyyy, \
                             tk_yzzz_xyyz, \
                             tk_yzzz_xyzz, \
                             tk_yzzz_xzzz, \
                             tk_yzzz_yyyy, \
                             tk_yzzz_yyyz, \
                             tk_yzzz_yyzz, \
                             tk_yzzz_yzzz, \
                             tk_yzzz_zzzz, \
                             tk_zzz_xxx,   \
                             tk_zzz_xxxx,  \
                             tk_zzz_xxxy,  \
                             tk_zzz_xxxz,  \
                             tk_zzz_xxy,   \
                             tk_zzz_xxyy,  \
                             tk_zzz_xxyz,  \
                             tk_zzz_xxz,   \
                             tk_zzz_xxzz,  \
                             tk_zzz_xyy,   \
                             tk_zzz_xyyy,  \
                             tk_zzz_xyyz,  \
                             tk_zzz_xyz,   \
                             tk_zzz_xyzz,  \
                             tk_zzz_xzz,   \
                             tk_zzz_xzzz,  \
                             tk_zzz_yyy,   \
                             tk_zzz_yyyy,  \
                             tk_zzz_yyyz,  \
                             tk_zzz_yyz,   \
                             tk_zzz_yyzz,  \
                             tk_zzz_yzz,   \
                             tk_zzz_yzzz,  \
                             tk_zzz_zzz,   \
                             tk_zzz_zzzz,  \
                             ts_yzzz_xxxx, \
                             ts_yzzz_xxxy, \
                             ts_yzzz_xxxz, \
                             ts_yzzz_xxyy, \
                             ts_yzzz_xxyz, \
                             ts_yzzz_xxzz, \
                             ts_yzzz_xyyy, \
                             ts_yzzz_xyyz, \
                             ts_yzzz_xyzz, \
                             ts_yzzz_xzzz, \
                             ts_yzzz_yyyy, \
                             ts_yzzz_yyyz, \
                             ts_yzzz_yyzz, \
                             ts_yzzz_yzzz, \
                             ts_yzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzz_xxxx[i] = tk_zzz_xxxx[i] * pa_y[i] + 2.0 * ts_yzzz_xxxx[i] * fz_0;

        tk_yzzz_xxxy[i] = tk_zzz_xxx[i] * fe_0 + tk_zzz_xxxy[i] * pa_y[i] + 2.0 * ts_yzzz_xxxy[i] * fz_0;

        tk_yzzz_xxxz[i] = tk_zzz_xxxz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxz[i] * fz_0;

        tk_yzzz_xxyy[i] = 2.0 * tk_zzz_xxy[i] * fe_0 + tk_zzz_xxyy[i] * pa_y[i] + 2.0 * ts_yzzz_xxyy[i] * fz_0;

        tk_yzzz_xxyz[i] = tk_zzz_xxz[i] * fe_0 + tk_zzz_xxyz[i] * pa_y[i] + 2.0 * ts_yzzz_xxyz[i] * fz_0;

        tk_yzzz_xxzz[i] = tk_zzz_xxzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxzz[i] * fz_0;

        tk_yzzz_xyyy[i] = 3.0 * tk_zzz_xyy[i] * fe_0 + tk_zzz_xyyy[i] * pa_y[i] + 2.0 * ts_yzzz_xyyy[i] * fz_0;

        tk_yzzz_xyyz[i] = 2.0 * tk_zzz_xyz[i] * fe_0 + tk_zzz_xyyz[i] * pa_y[i] + 2.0 * ts_yzzz_xyyz[i] * fz_0;

        tk_yzzz_xyzz[i] = tk_zzz_xzz[i] * fe_0 + tk_zzz_xyzz[i] * pa_y[i] + 2.0 * ts_yzzz_xyzz[i] * fz_0;

        tk_yzzz_xzzz[i] = tk_zzz_xzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xzzz[i] * fz_0;

        tk_yzzz_yyyy[i] = 4.0 * tk_zzz_yyy[i] * fe_0 + tk_zzz_yyyy[i] * pa_y[i] + 2.0 * ts_yzzz_yyyy[i] * fz_0;

        tk_yzzz_yyyz[i] = 3.0 * tk_zzz_yyz[i] * fe_0 + tk_zzz_yyyz[i] * pa_y[i] + 2.0 * ts_yzzz_yyyz[i] * fz_0;

        tk_yzzz_yyzz[i] = 2.0 * tk_zzz_yzz[i] * fe_0 + tk_zzz_yyzz[i] * pa_y[i] + 2.0 * ts_yzzz_yyzz[i] * fz_0;

        tk_yzzz_yzzz[i] = tk_zzz_zzz[i] * fe_0 + tk_zzz_yzzz[i] * pa_y[i] + 2.0 * ts_yzzz_yzzz[i] * fz_0;

        tk_yzzz_zzzz[i] = tk_zzz_zzzz[i] * pa_y[i] + 2.0 * ts_yzzz_zzzz[i] * fz_0;
    }

    // Set up 210-225 components of targeted buffer : GG

    auto tk_zzzz_xxxx = pbuffer.data(idx_kin_gg + 210);

    auto tk_zzzz_xxxy = pbuffer.data(idx_kin_gg + 211);

    auto tk_zzzz_xxxz = pbuffer.data(idx_kin_gg + 212);

    auto tk_zzzz_xxyy = pbuffer.data(idx_kin_gg + 213);

    auto tk_zzzz_xxyz = pbuffer.data(idx_kin_gg + 214);

    auto tk_zzzz_xxzz = pbuffer.data(idx_kin_gg + 215);

    auto tk_zzzz_xyyy = pbuffer.data(idx_kin_gg + 216);

    auto tk_zzzz_xyyz = pbuffer.data(idx_kin_gg + 217);

    auto tk_zzzz_xyzz = pbuffer.data(idx_kin_gg + 218);

    auto tk_zzzz_xzzz = pbuffer.data(idx_kin_gg + 219);

    auto tk_zzzz_yyyy = pbuffer.data(idx_kin_gg + 220);

    auto tk_zzzz_yyyz = pbuffer.data(idx_kin_gg + 221);

    auto tk_zzzz_yyzz = pbuffer.data(idx_kin_gg + 222);

    auto tk_zzzz_yzzz = pbuffer.data(idx_kin_gg + 223);

    auto tk_zzzz_zzzz = pbuffer.data(idx_kin_gg + 224);

#pragma omp simd aligned(pa_z,             \
                             tk_zz_xxxx,   \
                             tk_zz_xxxy,   \
                             tk_zz_xxxz,   \
                             tk_zz_xxyy,   \
                             tk_zz_xxyz,   \
                             tk_zz_xxzz,   \
                             tk_zz_xyyy,   \
                             tk_zz_xyyz,   \
                             tk_zz_xyzz,   \
                             tk_zz_xzzz,   \
                             tk_zz_yyyy,   \
                             tk_zz_yyyz,   \
                             tk_zz_yyzz,   \
                             tk_zz_yzzz,   \
                             tk_zz_zzzz,   \
                             tk_zzz_xxx,   \
                             tk_zzz_xxxx,  \
                             tk_zzz_xxxy,  \
                             tk_zzz_xxxz,  \
                             tk_zzz_xxy,   \
                             tk_zzz_xxyy,  \
                             tk_zzz_xxyz,  \
                             tk_zzz_xxz,   \
                             tk_zzz_xxzz,  \
                             tk_zzz_xyy,   \
                             tk_zzz_xyyy,  \
                             tk_zzz_xyyz,  \
                             tk_zzz_xyz,   \
                             tk_zzz_xyzz,  \
                             tk_zzz_xzz,   \
                             tk_zzz_xzzz,  \
                             tk_zzz_yyy,   \
                             tk_zzz_yyyy,  \
                             tk_zzz_yyyz,  \
                             tk_zzz_yyz,   \
                             tk_zzz_yyzz,  \
                             tk_zzz_yzz,   \
                             tk_zzz_yzzz,  \
                             tk_zzz_zzz,   \
                             tk_zzz_zzzz,  \
                             tk_zzzz_xxxx, \
                             tk_zzzz_xxxy, \
                             tk_zzzz_xxxz, \
                             tk_zzzz_xxyy, \
                             tk_zzzz_xxyz, \
                             tk_zzzz_xxzz, \
                             tk_zzzz_xyyy, \
                             tk_zzzz_xyyz, \
                             tk_zzzz_xyzz, \
                             tk_zzzz_xzzz, \
                             tk_zzzz_yyyy, \
                             tk_zzzz_yyyz, \
                             tk_zzzz_yyzz, \
                             tk_zzzz_yzzz, \
                             tk_zzzz_zzzz, \
                             ts_zz_xxxx,   \
                             ts_zz_xxxy,   \
                             ts_zz_xxxz,   \
                             ts_zz_xxyy,   \
                             ts_zz_xxyz,   \
                             ts_zz_xxzz,   \
                             ts_zz_xyyy,   \
                             ts_zz_xyyz,   \
                             ts_zz_xyzz,   \
                             ts_zz_xzzz,   \
                             ts_zz_yyyy,   \
                             ts_zz_yyyz,   \
                             ts_zz_yyzz,   \
                             ts_zz_yzzz,   \
                             ts_zz_zzzz,   \
                             ts_zzzz_xxxx, \
                             ts_zzzz_xxxy, \
                             ts_zzzz_xxxz, \
                             ts_zzzz_xxyy, \
                             ts_zzzz_xxyz, \
                             ts_zzzz_xxzz, \
                             ts_zzzz_xyyy, \
                             ts_zzzz_xyyz, \
                             ts_zzzz_xyzz, \
                             ts_zzzz_xzzz, \
                             ts_zzzz_yyyy, \
                             ts_zzzz_yyyz, \
                             ts_zzzz_yyzz, \
                             ts_zzzz_yzzz, \
                             ts_zzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzz_xxxx[i] = -6.0 * ts_zz_xxxx[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxx[i] * fe_0 + tk_zzz_xxxx[i] * pa_z[i] + 2.0 * ts_zzzz_xxxx[i] * fz_0;

        tk_zzzz_xxxy[i] = -6.0 * ts_zz_xxxy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxy[i] * fe_0 + tk_zzz_xxxy[i] * pa_z[i] + 2.0 * ts_zzzz_xxxy[i] * fz_0;

        tk_zzzz_xxxz[i] = -6.0 * ts_zz_xxxz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxz[i] * fe_0 + tk_zzz_xxx[i] * fe_0 + tk_zzz_xxxz[i] * pa_z[i] +
                          2.0 * ts_zzzz_xxxz[i] * fz_0;

        tk_zzzz_xxyy[i] = -6.0 * ts_zz_xxyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyy[i] * fe_0 + tk_zzz_xxyy[i] * pa_z[i] + 2.0 * ts_zzzz_xxyy[i] * fz_0;

        tk_zzzz_xxyz[i] = -6.0 * ts_zz_xxyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyz[i] * fe_0 + tk_zzz_xxy[i] * fe_0 + tk_zzz_xxyz[i] * pa_z[i] +
                          2.0 * ts_zzzz_xxyz[i] * fz_0;

        tk_zzzz_xxzz[i] = -6.0 * ts_zz_xxzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxzz[i] * fe_0 + 2.0 * tk_zzz_xxz[i] * fe_0 + tk_zzz_xxzz[i] * pa_z[i] +
                          2.0 * ts_zzzz_xxzz[i] * fz_0;

        tk_zzzz_xyyy[i] = -6.0 * ts_zz_xyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyy[i] * fe_0 + tk_zzz_xyyy[i] * pa_z[i] + 2.0 * ts_zzzz_xyyy[i] * fz_0;

        tk_zzzz_xyyz[i] = -6.0 * ts_zz_xyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyz[i] * fe_0 + tk_zzz_xyy[i] * fe_0 + tk_zzz_xyyz[i] * pa_z[i] +
                          2.0 * ts_zzzz_xyyz[i] * fz_0;

        tk_zzzz_xyzz[i] = -6.0 * ts_zz_xyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyzz[i] * fe_0 + 2.0 * tk_zzz_xyz[i] * fe_0 + tk_zzz_xyzz[i] * pa_z[i] +
                          2.0 * ts_zzzz_xyzz[i] * fz_0;

        tk_zzzz_xzzz[i] = -6.0 * ts_zz_xzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xzzz[i] * fe_0 + 3.0 * tk_zzz_xzz[i] * fe_0 + tk_zzz_xzzz[i] * pa_z[i] +
                          2.0 * ts_zzzz_xzzz[i] * fz_0;

        tk_zzzz_yyyy[i] = -6.0 * ts_zz_yyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyy[i] * fe_0 + tk_zzz_yyyy[i] * pa_z[i] + 2.0 * ts_zzzz_yyyy[i] * fz_0;

        tk_zzzz_yyyz[i] = -6.0 * ts_zz_yyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyz[i] * fe_0 + tk_zzz_yyy[i] * fe_0 + tk_zzz_yyyz[i] * pa_z[i] +
                          2.0 * ts_zzzz_yyyz[i] * fz_0;

        tk_zzzz_yyzz[i] = -6.0 * ts_zz_yyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyzz[i] * fe_0 + 2.0 * tk_zzz_yyz[i] * fe_0 + tk_zzz_yyzz[i] * pa_z[i] +
                          2.0 * ts_zzzz_yyzz[i] * fz_0;

        tk_zzzz_yzzz[i] = -6.0 * ts_zz_yzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yzzz[i] * fe_0 + 3.0 * tk_zzz_yzz[i] * fe_0 + tk_zzz_yzzz[i] * pa_z[i] +
                          2.0 * ts_zzzz_yzzz[i] * fz_0;

        tk_zzzz_zzzz[i] = -6.0 * ts_zz_zzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_zzzz[i] * fe_0 + 4.0 * tk_zzz_zzz[i] * fe_0 + tk_zzz_zzzz[i] * pa_z[i] +
                          2.0 * ts_zzzz_zzzz[i] * fz_0;
    }
}

}  // namespace kinrec
