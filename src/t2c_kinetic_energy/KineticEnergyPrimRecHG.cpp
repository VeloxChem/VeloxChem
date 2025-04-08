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

#include "KineticEnergyPrimRecHG.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_hg(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_hg,
                            const size_t              idx_ovl_fg,
                            const size_t              idx_kin_fg,
                            const size_t              idx_kin_gf,
                            const size_t              idx_kin_gg,
                            const size_t              idx_ovl_hg,
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

    // Set up components of auxiliary buffer : FG

    auto ts_xxx_xxxx = pbuffer.data(idx_ovl_fg);

    auto ts_xxx_xxxy = pbuffer.data(idx_ovl_fg + 1);

    auto ts_xxx_xxxz = pbuffer.data(idx_ovl_fg + 2);

    auto ts_xxx_xxyy = pbuffer.data(idx_ovl_fg + 3);

    auto ts_xxx_xxyz = pbuffer.data(idx_ovl_fg + 4);

    auto ts_xxx_xxzz = pbuffer.data(idx_ovl_fg + 5);

    auto ts_xxx_xyyy = pbuffer.data(idx_ovl_fg + 6);

    auto ts_xxx_xyyz = pbuffer.data(idx_ovl_fg + 7);

    auto ts_xxx_xyzz = pbuffer.data(idx_ovl_fg + 8);

    auto ts_xxx_xzzz = pbuffer.data(idx_ovl_fg + 9);

    auto ts_xxx_yyyy = pbuffer.data(idx_ovl_fg + 10);

    auto ts_xxx_yyyz = pbuffer.data(idx_ovl_fg + 11);

    auto ts_xxx_yyzz = pbuffer.data(idx_ovl_fg + 12);

    auto ts_xxx_yzzz = pbuffer.data(idx_ovl_fg + 13);

    auto ts_xxx_zzzz = pbuffer.data(idx_ovl_fg + 14);

    auto ts_xxy_xxxx = pbuffer.data(idx_ovl_fg + 15);

    auto ts_xxy_xxxz = pbuffer.data(idx_ovl_fg + 17);

    auto ts_xxy_xxzz = pbuffer.data(idx_ovl_fg + 20);

    auto ts_xxy_xzzz = pbuffer.data(idx_ovl_fg + 24);

    auto ts_xxz_xxxx = pbuffer.data(idx_ovl_fg + 30);

    auto ts_xxz_xxxy = pbuffer.data(idx_ovl_fg + 31);

    auto ts_xxz_xxyy = pbuffer.data(idx_ovl_fg + 33);

    auto ts_xxz_xyyy = pbuffer.data(idx_ovl_fg + 36);

    auto ts_xyy_xxxy = pbuffer.data(idx_ovl_fg + 46);

    auto ts_xyy_xxyy = pbuffer.data(idx_ovl_fg + 48);

    auto ts_xyy_xxyz = pbuffer.data(idx_ovl_fg + 49);

    auto ts_xyy_xyyy = pbuffer.data(idx_ovl_fg + 51);

    auto ts_xyy_xyyz = pbuffer.data(idx_ovl_fg + 52);

    auto ts_xyy_xyzz = pbuffer.data(idx_ovl_fg + 53);

    auto ts_xyy_yyyy = pbuffer.data(idx_ovl_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_ovl_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_ovl_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_ovl_fg + 58);

    auto ts_xyy_zzzz = pbuffer.data(idx_ovl_fg + 59);

    auto ts_xzz_xxxz = pbuffer.data(idx_ovl_fg + 77);

    auto ts_xzz_xxyz = pbuffer.data(idx_ovl_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_ovl_fg + 80);

    auto ts_xzz_xyyz = pbuffer.data(idx_ovl_fg + 82);

    auto ts_xzz_xyzz = pbuffer.data(idx_ovl_fg + 83);

    auto ts_xzz_xzzz = pbuffer.data(idx_ovl_fg + 84);

    auto ts_xzz_yyyy = pbuffer.data(idx_ovl_fg + 85);

    auto ts_xzz_yyyz = pbuffer.data(idx_ovl_fg + 86);

    auto ts_xzz_yyzz = pbuffer.data(idx_ovl_fg + 87);

    auto ts_xzz_yzzz = pbuffer.data(idx_ovl_fg + 88);

    auto ts_xzz_zzzz = pbuffer.data(idx_ovl_fg + 89);

    auto ts_yyy_xxxx = pbuffer.data(idx_ovl_fg + 90);

    auto ts_yyy_xxxy = pbuffer.data(idx_ovl_fg + 91);

    auto ts_yyy_xxxz = pbuffer.data(idx_ovl_fg + 92);

    auto ts_yyy_xxyy = pbuffer.data(idx_ovl_fg + 93);

    auto ts_yyy_xxyz = pbuffer.data(idx_ovl_fg + 94);

    auto ts_yyy_xxzz = pbuffer.data(idx_ovl_fg + 95);

    auto ts_yyy_xyyy = pbuffer.data(idx_ovl_fg + 96);

    auto ts_yyy_xyyz = pbuffer.data(idx_ovl_fg + 97);

    auto ts_yyy_xyzz = pbuffer.data(idx_ovl_fg + 98);

    auto ts_yyy_xzzz = pbuffer.data(idx_ovl_fg + 99);

    auto ts_yyy_yyyy = pbuffer.data(idx_ovl_fg + 100);

    auto ts_yyy_yyyz = pbuffer.data(idx_ovl_fg + 101);

    auto ts_yyy_yyzz = pbuffer.data(idx_ovl_fg + 102);

    auto ts_yyy_yzzz = pbuffer.data(idx_ovl_fg + 103);

    auto ts_yyy_zzzz = pbuffer.data(idx_ovl_fg + 104);

    auto ts_yyz_xxxy = pbuffer.data(idx_ovl_fg + 106);

    auto ts_yyz_xxyy = pbuffer.data(idx_ovl_fg + 108);

    auto ts_yyz_xyyy = pbuffer.data(idx_ovl_fg + 111);

    auto ts_yyz_yyyy = pbuffer.data(idx_ovl_fg + 115);

    auto ts_yzz_xxxx = pbuffer.data(idx_ovl_fg + 120);

    auto ts_yzz_xxxz = pbuffer.data(idx_ovl_fg + 122);

    auto ts_yzz_xxyz = pbuffer.data(idx_ovl_fg + 124);

    auto ts_yzz_xxzz = pbuffer.data(idx_ovl_fg + 125);

    auto ts_yzz_xyyz = pbuffer.data(idx_ovl_fg + 127);

    auto ts_yzz_xyzz = pbuffer.data(idx_ovl_fg + 128);

    auto ts_yzz_xzzz = pbuffer.data(idx_ovl_fg + 129);

    auto ts_yzz_yyyz = pbuffer.data(idx_ovl_fg + 131);

    auto ts_yzz_yyzz = pbuffer.data(idx_ovl_fg + 132);

    auto ts_yzz_yzzz = pbuffer.data(idx_ovl_fg + 133);

    auto ts_yzz_zzzz = pbuffer.data(idx_ovl_fg + 134);

    auto ts_zzz_xxxx = pbuffer.data(idx_ovl_fg + 135);

    auto ts_zzz_xxxy = pbuffer.data(idx_ovl_fg + 136);

    auto ts_zzz_xxxz = pbuffer.data(idx_ovl_fg + 137);

    auto ts_zzz_xxyy = pbuffer.data(idx_ovl_fg + 138);

    auto ts_zzz_xxyz = pbuffer.data(idx_ovl_fg + 139);

    auto ts_zzz_xxzz = pbuffer.data(idx_ovl_fg + 140);

    auto ts_zzz_xyyy = pbuffer.data(idx_ovl_fg + 141);

    auto ts_zzz_xyyz = pbuffer.data(idx_ovl_fg + 142);

    auto ts_zzz_xyzz = pbuffer.data(idx_ovl_fg + 143);

    auto ts_zzz_xzzz = pbuffer.data(idx_ovl_fg + 144);

    auto ts_zzz_yyyy = pbuffer.data(idx_ovl_fg + 145);

    auto ts_zzz_yyyz = pbuffer.data(idx_ovl_fg + 146);

    auto ts_zzz_yyzz = pbuffer.data(idx_ovl_fg + 147);

    auto ts_zzz_yzzz = pbuffer.data(idx_ovl_fg + 148);

    auto ts_zzz_zzzz = pbuffer.data(idx_ovl_fg + 149);

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

    auto tk_xxy_xxxz = pbuffer.data(idx_kin_fg + 17);

    auto tk_xxy_xxzz = pbuffer.data(idx_kin_fg + 20);

    auto tk_xxy_xzzz = pbuffer.data(idx_kin_fg + 24);

    auto tk_xxz_xxxx = pbuffer.data(idx_kin_fg + 30);

    auto tk_xxz_xxxy = pbuffer.data(idx_kin_fg + 31);

    auto tk_xxz_xxyy = pbuffer.data(idx_kin_fg + 33);

    auto tk_xxz_xyyy = pbuffer.data(idx_kin_fg + 36);

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

    auto tk_yyz_xxyy = pbuffer.data(idx_kin_fg + 108);

    auto tk_yyz_xyyy = pbuffer.data(idx_kin_fg + 111);

    auto tk_yyz_yyyy = pbuffer.data(idx_kin_fg + 115);

    auto tk_yzz_xxxx = pbuffer.data(idx_kin_fg + 120);

    auto tk_yzz_xxxz = pbuffer.data(idx_kin_fg + 122);

    auto tk_yzz_xxyz = pbuffer.data(idx_kin_fg + 124);

    auto tk_yzz_xxzz = pbuffer.data(idx_kin_fg + 125);

    auto tk_yzz_xyyz = pbuffer.data(idx_kin_fg + 127);

    auto tk_yzz_xyzz = pbuffer.data(idx_kin_fg + 128);

    auto tk_yzz_xzzz = pbuffer.data(idx_kin_fg + 129);

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

    // Set up components of auxiliary buffer : GF

    auto tk_xxxx_xxx = pbuffer.data(idx_kin_gf);

    auto tk_xxxx_xxy = pbuffer.data(idx_kin_gf + 1);

    auto tk_xxxx_xxz = pbuffer.data(idx_kin_gf + 2);

    auto tk_xxxx_xyy = pbuffer.data(idx_kin_gf + 3);

    auto tk_xxxx_xyz = pbuffer.data(idx_kin_gf + 4);

    auto tk_xxxx_xzz = pbuffer.data(idx_kin_gf + 5);

    auto tk_xxxx_yyy = pbuffer.data(idx_kin_gf + 6);

    auto tk_xxxx_yyz = pbuffer.data(idx_kin_gf + 7);

    auto tk_xxxx_yzz = pbuffer.data(idx_kin_gf + 8);

    auto tk_xxxx_zzz = pbuffer.data(idx_kin_gf + 9);

    auto tk_xxxz_xxz = pbuffer.data(idx_kin_gf + 22);

    auto tk_xxxz_xyz = pbuffer.data(idx_kin_gf + 24);

    auto tk_xxxz_xzz = pbuffer.data(idx_kin_gf + 25);

    auto tk_xxxz_yyz = pbuffer.data(idx_kin_gf + 27);

    auto tk_xxxz_yzz = pbuffer.data(idx_kin_gf + 28);

    auto tk_xxxz_zzz = pbuffer.data(idx_kin_gf + 29);

    auto tk_xxyy_xxx = pbuffer.data(idx_kin_gf + 30);

    auto tk_xxyy_xxy = pbuffer.data(idx_kin_gf + 31);

    auto tk_xxyy_xxz = pbuffer.data(idx_kin_gf + 32);

    auto tk_xxyy_xyy = pbuffer.data(idx_kin_gf + 33);

    auto tk_xxyy_xyz = pbuffer.data(idx_kin_gf + 34);

    auto tk_xxyy_xzz = pbuffer.data(idx_kin_gf + 35);

    auto tk_xxyy_yyy = pbuffer.data(idx_kin_gf + 36);

    auto tk_xxyy_yyz = pbuffer.data(idx_kin_gf + 37);

    auto tk_xxyy_yzz = pbuffer.data(idx_kin_gf + 38);

    auto tk_xxyy_zzz = pbuffer.data(idx_kin_gf + 39);

    auto tk_xxzz_xxx = pbuffer.data(idx_kin_gf + 50);

    auto tk_xxzz_xxy = pbuffer.data(idx_kin_gf + 51);

    auto tk_xxzz_xxz = pbuffer.data(idx_kin_gf + 52);

    auto tk_xxzz_xyy = pbuffer.data(idx_kin_gf + 53);

    auto tk_xxzz_xyz = pbuffer.data(idx_kin_gf + 54);

    auto tk_xxzz_xzz = pbuffer.data(idx_kin_gf + 55);

    auto tk_xxzz_yyy = pbuffer.data(idx_kin_gf + 56);

    auto tk_xxzz_yyz = pbuffer.data(idx_kin_gf + 57);

    auto tk_xxzz_yzz = pbuffer.data(idx_kin_gf + 58);

    auto tk_xxzz_zzz = pbuffer.data(idx_kin_gf + 59);

    auto tk_xyyy_xxy = pbuffer.data(idx_kin_gf + 61);

    auto tk_xyyy_xyy = pbuffer.data(idx_kin_gf + 63);

    auto tk_xyyy_xyz = pbuffer.data(idx_kin_gf + 64);

    auto tk_xyyy_yyy = pbuffer.data(idx_kin_gf + 66);

    auto tk_xyyy_yyz = pbuffer.data(idx_kin_gf + 67);

    auto tk_xyyy_yzz = pbuffer.data(idx_kin_gf + 68);

    auto tk_xzzz_xxz = pbuffer.data(idx_kin_gf + 92);

    auto tk_xzzz_xyz = pbuffer.data(idx_kin_gf + 94);

    auto tk_xzzz_xzz = pbuffer.data(idx_kin_gf + 95);

    auto tk_xzzz_yyz = pbuffer.data(idx_kin_gf + 97);

    auto tk_xzzz_yzz = pbuffer.data(idx_kin_gf + 98);

    auto tk_xzzz_zzz = pbuffer.data(idx_kin_gf + 99);

    auto tk_yyyy_xxx = pbuffer.data(idx_kin_gf + 100);

    auto tk_yyyy_xxy = pbuffer.data(idx_kin_gf + 101);

    auto tk_yyyy_xxz = pbuffer.data(idx_kin_gf + 102);

    auto tk_yyyy_xyy = pbuffer.data(idx_kin_gf + 103);

    auto tk_yyyy_xyz = pbuffer.data(idx_kin_gf + 104);

    auto tk_yyyy_xzz = pbuffer.data(idx_kin_gf + 105);

    auto tk_yyyy_yyy = pbuffer.data(idx_kin_gf + 106);

    auto tk_yyyy_yyz = pbuffer.data(idx_kin_gf + 107);

    auto tk_yyyy_yzz = pbuffer.data(idx_kin_gf + 108);

    auto tk_yyyy_zzz = pbuffer.data(idx_kin_gf + 109);

    auto tk_yyyz_xxz = pbuffer.data(idx_kin_gf + 112);

    auto tk_yyyz_xyz = pbuffer.data(idx_kin_gf + 114);

    auto tk_yyyz_xzz = pbuffer.data(idx_kin_gf + 115);

    auto tk_yyyz_yyz = pbuffer.data(idx_kin_gf + 117);

    auto tk_yyyz_yzz = pbuffer.data(idx_kin_gf + 118);

    auto tk_yyyz_zzz = pbuffer.data(idx_kin_gf + 119);

    auto tk_yyzz_xxx = pbuffer.data(idx_kin_gf + 120);

    auto tk_yyzz_xxy = pbuffer.data(idx_kin_gf + 121);

    auto tk_yyzz_xxz = pbuffer.data(idx_kin_gf + 122);

    auto tk_yyzz_xyy = pbuffer.data(idx_kin_gf + 123);

    auto tk_yyzz_xyz = pbuffer.data(idx_kin_gf + 124);

    auto tk_yyzz_xzz = pbuffer.data(idx_kin_gf + 125);

    auto tk_yyzz_yyy = pbuffer.data(idx_kin_gf + 126);

    auto tk_yyzz_yyz = pbuffer.data(idx_kin_gf + 127);

    auto tk_yyzz_yzz = pbuffer.data(idx_kin_gf + 128);

    auto tk_yyzz_zzz = pbuffer.data(idx_kin_gf + 129);

    auto tk_yzzz_xxy = pbuffer.data(idx_kin_gf + 131);

    auto tk_yzzz_xxz = pbuffer.data(idx_kin_gf + 132);

    auto tk_yzzz_xyy = pbuffer.data(idx_kin_gf + 133);

    auto tk_yzzz_xyz = pbuffer.data(idx_kin_gf + 134);

    auto tk_yzzz_xzz = pbuffer.data(idx_kin_gf + 135);

    auto tk_yzzz_yyy = pbuffer.data(idx_kin_gf + 136);

    auto tk_yzzz_yyz = pbuffer.data(idx_kin_gf + 137);

    auto tk_yzzz_yzz = pbuffer.data(idx_kin_gf + 138);

    auto tk_yzzz_zzz = pbuffer.data(idx_kin_gf + 139);

    auto tk_zzzz_xxx = pbuffer.data(idx_kin_gf + 140);

    auto tk_zzzz_xxy = pbuffer.data(idx_kin_gf + 141);

    auto tk_zzzz_xxz = pbuffer.data(idx_kin_gf + 142);

    auto tk_zzzz_xyy = pbuffer.data(idx_kin_gf + 143);

    auto tk_zzzz_xyz = pbuffer.data(idx_kin_gf + 144);

    auto tk_zzzz_xzz = pbuffer.data(idx_kin_gf + 145);

    auto tk_zzzz_yyy = pbuffer.data(idx_kin_gf + 146);

    auto tk_zzzz_yyz = pbuffer.data(idx_kin_gf + 147);

    auto tk_zzzz_yzz = pbuffer.data(idx_kin_gf + 148);

    auto tk_zzzz_zzz = pbuffer.data(idx_kin_gf + 149);

    // Set up components of auxiliary buffer : GG

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

    auto tk_xxxy_xxxx = pbuffer.data(idx_kin_gg + 15);

    auto tk_xxxy_xxxy = pbuffer.data(idx_kin_gg + 16);

    auto tk_xxxy_xxxz = pbuffer.data(idx_kin_gg + 17);

    auto tk_xxxy_xxyy = pbuffer.data(idx_kin_gg + 18);

    auto tk_xxxy_xxzz = pbuffer.data(idx_kin_gg + 20);

    auto tk_xxxy_xyyy = pbuffer.data(idx_kin_gg + 21);

    auto tk_xxxy_xzzz = pbuffer.data(idx_kin_gg + 24);

    auto tk_xxxy_yyyy = pbuffer.data(idx_kin_gg + 25);

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

    auto tk_xxxz_yyyz = pbuffer.data(idx_kin_gg + 41);

    auto tk_xxxz_yyzz = pbuffer.data(idx_kin_gg + 42);

    auto tk_xxxz_yzzz = pbuffer.data(idx_kin_gg + 43);

    auto tk_xxxz_zzzz = pbuffer.data(idx_kin_gg + 44);

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

    auto tk_xyyy_xxxx = pbuffer.data(idx_kin_gg + 90);

    auto tk_xyyy_xxxy = pbuffer.data(idx_kin_gg + 91);

    auto tk_xyyy_xxyy = pbuffer.data(idx_kin_gg + 93);

    auto tk_xyyy_xxyz = pbuffer.data(idx_kin_gg + 94);

    auto tk_xyyy_xyyy = pbuffer.data(idx_kin_gg + 96);

    auto tk_xyyy_xyyz = pbuffer.data(idx_kin_gg + 97);

    auto tk_xyyy_xyzz = pbuffer.data(idx_kin_gg + 98);

    auto tk_xyyy_yyyy = pbuffer.data(idx_kin_gg + 100);

    auto tk_xyyy_yyyz = pbuffer.data(idx_kin_gg + 101);

    auto tk_xyyy_yyzz = pbuffer.data(idx_kin_gg + 102);

    auto tk_xyyy_yzzz = pbuffer.data(idx_kin_gg + 103);

    auto tk_xyyy_zzzz = pbuffer.data(idx_kin_gg + 104);

    auto tk_xzzz_xxxx = pbuffer.data(idx_kin_gg + 135);

    auto tk_xzzz_xxxz = pbuffer.data(idx_kin_gg + 137);

    auto tk_xzzz_xxyz = pbuffer.data(idx_kin_gg + 139);

    auto tk_xzzz_xxzz = pbuffer.data(idx_kin_gg + 140);

    auto tk_xzzz_xyyz = pbuffer.data(idx_kin_gg + 142);

    auto tk_xzzz_xyzz = pbuffer.data(idx_kin_gg + 143);

    auto tk_xzzz_xzzz = pbuffer.data(idx_kin_gg + 144);

    auto tk_xzzz_yyyy = pbuffer.data(idx_kin_gg + 145);

    auto tk_xzzz_yyyz = pbuffer.data(idx_kin_gg + 146);

    auto tk_xzzz_yyzz = pbuffer.data(idx_kin_gg + 147);

    auto tk_xzzz_yzzz = pbuffer.data(idx_kin_gg + 148);

    auto tk_xzzz_zzzz = pbuffer.data(idx_kin_gg + 149);

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

    // Set up components of auxiliary buffer : HG

    auto ts_xxxxx_xxxx = pbuffer.data(idx_ovl_hg);

    auto ts_xxxxx_xxxy = pbuffer.data(idx_ovl_hg + 1);

    auto ts_xxxxx_xxxz = pbuffer.data(idx_ovl_hg + 2);

    auto ts_xxxxx_xxyy = pbuffer.data(idx_ovl_hg + 3);

    auto ts_xxxxx_xxyz = pbuffer.data(idx_ovl_hg + 4);

    auto ts_xxxxx_xxzz = pbuffer.data(idx_ovl_hg + 5);

    auto ts_xxxxx_xyyy = pbuffer.data(idx_ovl_hg + 6);

    auto ts_xxxxx_xyyz = pbuffer.data(idx_ovl_hg + 7);

    auto ts_xxxxx_xyzz = pbuffer.data(idx_ovl_hg + 8);

    auto ts_xxxxx_xzzz = pbuffer.data(idx_ovl_hg + 9);

    auto ts_xxxxx_yyyy = pbuffer.data(idx_ovl_hg + 10);

    auto ts_xxxxx_yyyz = pbuffer.data(idx_ovl_hg + 11);

    auto ts_xxxxx_yyzz = pbuffer.data(idx_ovl_hg + 12);

    auto ts_xxxxx_yzzz = pbuffer.data(idx_ovl_hg + 13);

    auto ts_xxxxx_zzzz = pbuffer.data(idx_ovl_hg + 14);

    auto ts_xxxxy_xxxx = pbuffer.data(idx_ovl_hg + 15);

    auto ts_xxxxy_xxxy = pbuffer.data(idx_ovl_hg + 16);

    auto ts_xxxxy_xxxz = pbuffer.data(idx_ovl_hg + 17);

    auto ts_xxxxy_xxyy = pbuffer.data(idx_ovl_hg + 18);

    auto ts_xxxxy_xxyz = pbuffer.data(idx_ovl_hg + 19);

    auto ts_xxxxy_xxzz = pbuffer.data(idx_ovl_hg + 20);

    auto ts_xxxxy_xyyy = pbuffer.data(idx_ovl_hg + 21);

    auto ts_xxxxy_xyyz = pbuffer.data(idx_ovl_hg + 22);

    auto ts_xxxxy_xyzz = pbuffer.data(idx_ovl_hg + 23);

    auto ts_xxxxy_xzzz = pbuffer.data(idx_ovl_hg + 24);

    auto ts_xxxxy_yyyy = pbuffer.data(idx_ovl_hg + 25);

    auto ts_xxxxy_yyyz = pbuffer.data(idx_ovl_hg + 26);

    auto ts_xxxxy_yyzz = pbuffer.data(idx_ovl_hg + 27);

    auto ts_xxxxy_yzzz = pbuffer.data(idx_ovl_hg + 28);

    auto ts_xxxxy_zzzz = pbuffer.data(idx_ovl_hg + 29);

    auto ts_xxxxz_xxxx = pbuffer.data(idx_ovl_hg + 30);

    auto ts_xxxxz_xxxy = pbuffer.data(idx_ovl_hg + 31);

    auto ts_xxxxz_xxxz = pbuffer.data(idx_ovl_hg + 32);

    auto ts_xxxxz_xxyy = pbuffer.data(idx_ovl_hg + 33);

    auto ts_xxxxz_xxyz = pbuffer.data(idx_ovl_hg + 34);

    auto ts_xxxxz_xxzz = pbuffer.data(idx_ovl_hg + 35);

    auto ts_xxxxz_xyyy = pbuffer.data(idx_ovl_hg + 36);

    auto ts_xxxxz_xyyz = pbuffer.data(idx_ovl_hg + 37);

    auto ts_xxxxz_xyzz = pbuffer.data(idx_ovl_hg + 38);

    auto ts_xxxxz_xzzz = pbuffer.data(idx_ovl_hg + 39);

    auto ts_xxxxz_yyyy = pbuffer.data(idx_ovl_hg + 40);

    auto ts_xxxxz_yyyz = pbuffer.data(idx_ovl_hg + 41);

    auto ts_xxxxz_yyzz = pbuffer.data(idx_ovl_hg + 42);

    auto ts_xxxxz_yzzz = pbuffer.data(idx_ovl_hg + 43);

    auto ts_xxxxz_zzzz = pbuffer.data(idx_ovl_hg + 44);

    auto ts_xxxyy_xxxx = pbuffer.data(idx_ovl_hg + 45);

    auto ts_xxxyy_xxxy = pbuffer.data(idx_ovl_hg + 46);

    auto ts_xxxyy_xxxz = pbuffer.data(idx_ovl_hg + 47);

    auto ts_xxxyy_xxyy = pbuffer.data(idx_ovl_hg + 48);

    auto ts_xxxyy_xxyz = pbuffer.data(idx_ovl_hg + 49);

    auto ts_xxxyy_xxzz = pbuffer.data(idx_ovl_hg + 50);

    auto ts_xxxyy_xyyy = pbuffer.data(idx_ovl_hg + 51);

    auto ts_xxxyy_xyyz = pbuffer.data(idx_ovl_hg + 52);

    auto ts_xxxyy_xyzz = pbuffer.data(idx_ovl_hg + 53);

    auto ts_xxxyy_xzzz = pbuffer.data(idx_ovl_hg + 54);

    auto ts_xxxyy_yyyy = pbuffer.data(idx_ovl_hg + 55);

    auto ts_xxxyy_yyyz = pbuffer.data(idx_ovl_hg + 56);

    auto ts_xxxyy_yyzz = pbuffer.data(idx_ovl_hg + 57);

    auto ts_xxxyy_yzzz = pbuffer.data(idx_ovl_hg + 58);

    auto ts_xxxyy_zzzz = pbuffer.data(idx_ovl_hg + 59);

    auto ts_xxxyz_xxxx = pbuffer.data(idx_ovl_hg + 60);

    auto ts_xxxyz_xxxy = pbuffer.data(idx_ovl_hg + 61);

    auto ts_xxxyz_xxxz = pbuffer.data(idx_ovl_hg + 62);

    auto ts_xxxyz_xxyy = pbuffer.data(idx_ovl_hg + 63);

    auto ts_xxxyz_xxyz = pbuffer.data(idx_ovl_hg + 64);

    auto ts_xxxyz_xxzz = pbuffer.data(idx_ovl_hg + 65);

    auto ts_xxxyz_xyyy = pbuffer.data(idx_ovl_hg + 66);

    auto ts_xxxyz_xyyz = pbuffer.data(idx_ovl_hg + 67);

    auto ts_xxxyz_xyzz = pbuffer.data(idx_ovl_hg + 68);

    auto ts_xxxyz_xzzz = pbuffer.data(idx_ovl_hg + 69);

    auto ts_xxxyz_yyyy = pbuffer.data(idx_ovl_hg + 70);

    auto ts_xxxyz_yyyz = pbuffer.data(idx_ovl_hg + 71);

    auto ts_xxxyz_yyzz = pbuffer.data(idx_ovl_hg + 72);

    auto ts_xxxyz_yzzz = pbuffer.data(idx_ovl_hg + 73);

    auto ts_xxxyz_zzzz = pbuffer.data(idx_ovl_hg + 74);

    auto ts_xxxzz_xxxx = pbuffer.data(idx_ovl_hg + 75);

    auto ts_xxxzz_xxxy = pbuffer.data(idx_ovl_hg + 76);

    auto ts_xxxzz_xxxz = pbuffer.data(idx_ovl_hg + 77);

    auto ts_xxxzz_xxyy = pbuffer.data(idx_ovl_hg + 78);

    auto ts_xxxzz_xxyz = pbuffer.data(idx_ovl_hg + 79);

    auto ts_xxxzz_xxzz = pbuffer.data(idx_ovl_hg + 80);

    auto ts_xxxzz_xyyy = pbuffer.data(idx_ovl_hg + 81);

    auto ts_xxxzz_xyyz = pbuffer.data(idx_ovl_hg + 82);

    auto ts_xxxzz_xyzz = pbuffer.data(idx_ovl_hg + 83);

    auto ts_xxxzz_xzzz = pbuffer.data(idx_ovl_hg + 84);

    auto ts_xxxzz_yyyy = pbuffer.data(idx_ovl_hg + 85);

    auto ts_xxxzz_yyyz = pbuffer.data(idx_ovl_hg + 86);

    auto ts_xxxzz_yyzz = pbuffer.data(idx_ovl_hg + 87);

    auto ts_xxxzz_yzzz = pbuffer.data(idx_ovl_hg + 88);

    auto ts_xxxzz_zzzz = pbuffer.data(idx_ovl_hg + 89);

    auto ts_xxyyy_xxxx = pbuffer.data(idx_ovl_hg + 90);

    auto ts_xxyyy_xxxy = pbuffer.data(idx_ovl_hg + 91);

    auto ts_xxyyy_xxxz = pbuffer.data(idx_ovl_hg + 92);

    auto ts_xxyyy_xxyy = pbuffer.data(idx_ovl_hg + 93);

    auto ts_xxyyy_xxyz = pbuffer.data(idx_ovl_hg + 94);

    auto ts_xxyyy_xxzz = pbuffer.data(idx_ovl_hg + 95);

    auto ts_xxyyy_xyyy = pbuffer.data(idx_ovl_hg + 96);

    auto ts_xxyyy_xyyz = pbuffer.data(idx_ovl_hg + 97);

    auto ts_xxyyy_xyzz = pbuffer.data(idx_ovl_hg + 98);

    auto ts_xxyyy_xzzz = pbuffer.data(idx_ovl_hg + 99);

    auto ts_xxyyy_yyyy = pbuffer.data(idx_ovl_hg + 100);

    auto ts_xxyyy_yyyz = pbuffer.data(idx_ovl_hg + 101);

    auto ts_xxyyy_yyzz = pbuffer.data(idx_ovl_hg + 102);

    auto ts_xxyyy_yzzz = pbuffer.data(idx_ovl_hg + 103);

    auto ts_xxyyy_zzzz = pbuffer.data(idx_ovl_hg + 104);

    auto ts_xxyyz_xxxx = pbuffer.data(idx_ovl_hg + 105);

    auto ts_xxyyz_xxxy = pbuffer.data(idx_ovl_hg + 106);

    auto ts_xxyyz_xxxz = pbuffer.data(idx_ovl_hg + 107);

    auto ts_xxyyz_xxyy = pbuffer.data(idx_ovl_hg + 108);

    auto ts_xxyyz_xxyz = pbuffer.data(idx_ovl_hg + 109);

    auto ts_xxyyz_xxzz = pbuffer.data(idx_ovl_hg + 110);

    auto ts_xxyyz_xyyy = pbuffer.data(idx_ovl_hg + 111);

    auto ts_xxyyz_xyyz = pbuffer.data(idx_ovl_hg + 112);

    auto ts_xxyyz_xyzz = pbuffer.data(idx_ovl_hg + 113);

    auto ts_xxyyz_xzzz = pbuffer.data(idx_ovl_hg + 114);

    auto ts_xxyyz_yyyy = pbuffer.data(idx_ovl_hg + 115);

    auto ts_xxyyz_yyyz = pbuffer.data(idx_ovl_hg + 116);

    auto ts_xxyyz_yyzz = pbuffer.data(idx_ovl_hg + 117);

    auto ts_xxyyz_yzzz = pbuffer.data(idx_ovl_hg + 118);

    auto ts_xxyyz_zzzz = pbuffer.data(idx_ovl_hg + 119);

    auto ts_xxyzz_xxxx = pbuffer.data(idx_ovl_hg + 120);

    auto ts_xxyzz_xxxy = pbuffer.data(idx_ovl_hg + 121);

    auto ts_xxyzz_xxxz = pbuffer.data(idx_ovl_hg + 122);

    auto ts_xxyzz_xxyy = pbuffer.data(idx_ovl_hg + 123);

    auto ts_xxyzz_xxyz = pbuffer.data(idx_ovl_hg + 124);

    auto ts_xxyzz_xxzz = pbuffer.data(idx_ovl_hg + 125);

    auto ts_xxyzz_xyyy = pbuffer.data(idx_ovl_hg + 126);

    auto ts_xxyzz_xyyz = pbuffer.data(idx_ovl_hg + 127);

    auto ts_xxyzz_xyzz = pbuffer.data(idx_ovl_hg + 128);

    auto ts_xxyzz_xzzz = pbuffer.data(idx_ovl_hg + 129);

    auto ts_xxyzz_yyyy = pbuffer.data(idx_ovl_hg + 130);

    auto ts_xxyzz_yyyz = pbuffer.data(idx_ovl_hg + 131);

    auto ts_xxyzz_yyzz = pbuffer.data(idx_ovl_hg + 132);

    auto ts_xxyzz_yzzz = pbuffer.data(idx_ovl_hg + 133);

    auto ts_xxyzz_zzzz = pbuffer.data(idx_ovl_hg + 134);

    auto ts_xxzzz_xxxx = pbuffer.data(idx_ovl_hg + 135);

    auto ts_xxzzz_xxxy = pbuffer.data(idx_ovl_hg + 136);

    auto ts_xxzzz_xxxz = pbuffer.data(idx_ovl_hg + 137);

    auto ts_xxzzz_xxyy = pbuffer.data(idx_ovl_hg + 138);

    auto ts_xxzzz_xxyz = pbuffer.data(idx_ovl_hg + 139);

    auto ts_xxzzz_xxzz = pbuffer.data(idx_ovl_hg + 140);

    auto ts_xxzzz_xyyy = pbuffer.data(idx_ovl_hg + 141);

    auto ts_xxzzz_xyyz = pbuffer.data(idx_ovl_hg + 142);

    auto ts_xxzzz_xyzz = pbuffer.data(idx_ovl_hg + 143);

    auto ts_xxzzz_xzzz = pbuffer.data(idx_ovl_hg + 144);

    auto ts_xxzzz_yyyy = pbuffer.data(idx_ovl_hg + 145);

    auto ts_xxzzz_yyyz = pbuffer.data(idx_ovl_hg + 146);

    auto ts_xxzzz_yyzz = pbuffer.data(idx_ovl_hg + 147);

    auto ts_xxzzz_yzzz = pbuffer.data(idx_ovl_hg + 148);

    auto ts_xxzzz_zzzz = pbuffer.data(idx_ovl_hg + 149);

    auto ts_xyyyy_xxxx = pbuffer.data(idx_ovl_hg + 150);

    auto ts_xyyyy_xxxy = pbuffer.data(idx_ovl_hg + 151);

    auto ts_xyyyy_xxxz = pbuffer.data(idx_ovl_hg + 152);

    auto ts_xyyyy_xxyy = pbuffer.data(idx_ovl_hg + 153);

    auto ts_xyyyy_xxyz = pbuffer.data(idx_ovl_hg + 154);

    auto ts_xyyyy_xxzz = pbuffer.data(idx_ovl_hg + 155);

    auto ts_xyyyy_xyyy = pbuffer.data(idx_ovl_hg + 156);

    auto ts_xyyyy_xyyz = pbuffer.data(idx_ovl_hg + 157);

    auto ts_xyyyy_xyzz = pbuffer.data(idx_ovl_hg + 158);

    auto ts_xyyyy_xzzz = pbuffer.data(idx_ovl_hg + 159);

    auto ts_xyyyy_yyyy = pbuffer.data(idx_ovl_hg + 160);

    auto ts_xyyyy_yyyz = pbuffer.data(idx_ovl_hg + 161);

    auto ts_xyyyy_yyzz = pbuffer.data(idx_ovl_hg + 162);

    auto ts_xyyyy_yzzz = pbuffer.data(idx_ovl_hg + 163);

    auto ts_xyyyy_zzzz = pbuffer.data(idx_ovl_hg + 164);

    auto ts_xyyyz_xxxx = pbuffer.data(idx_ovl_hg + 165);

    auto ts_xyyyz_xxxy = pbuffer.data(idx_ovl_hg + 166);

    auto ts_xyyyz_xxxz = pbuffer.data(idx_ovl_hg + 167);

    auto ts_xyyyz_xxyy = pbuffer.data(idx_ovl_hg + 168);

    auto ts_xyyyz_xxyz = pbuffer.data(idx_ovl_hg + 169);

    auto ts_xyyyz_xxzz = pbuffer.data(idx_ovl_hg + 170);

    auto ts_xyyyz_xyyy = pbuffer.data(idx_ovl_hg + 171);

    auto ts_xyyyz_xyyz = pbuffer.data(idx_ovl_hg + 172);

    auto ts_xyyyz_xyzz = pbuffer.data(idx_ovl_hg + 173);

    auto ts_xyyyz_xzzz = pbuffer.data(idx_ovl_hg + 174);

    auto ts_xyyyz_yyyy = pbuffer.data(idx_ovl_hg + 175);

    auto ts_xyyyz_yyyz = pbuffer.data(idx_ovl_hg + 176);

    auto ts_xyyyz_yyzz = pbuffer.data(idx_ovl_hg + 177);

    auto ts_xyyyz_yzzz = pbuffer.data(idx_ovl_hg + 178);

    auto ts_xyyyz_zzzz = pbuffer.data(idx_ovl_hg + 179);

    auto ts_xyyzz_xxxx = pbuffer.data(idx_ovl_hg + 180);

    auto ts_xyyzz_xxxy = pbuffer.data(idx_ovl_hg + 181);

    auto ts_xyyzz_xxxz = pbuffer.data(idx_ovl_hg + 182);

    auto ts_xyyzz_xxyy = pbuffer.data(idx_ovl_hg + 183);

    auto ts_xyyzz_xxyz = pbuffer.data(idx_ovl_hg + 184);

    auto ts_xyyzz_xxzz = pbuffer.data(idx_ovl_hg + 185);

    auto ts_xyyzz_xyyy = pbuffer.data(idx_ovl_hg + 186);

    auto ts_xyyzz_xyyz = pbuffer.data(idx_ovl_hg + 187);

    auto ts_xyyzz_xyzz = pbuffer.data(idx_ovl_hg + 188);

    auto ts_xyyzz_xzzz = pbuffer.data(idx_ovl_hg + 189);

    auto ts_xyyzz_yyyy = pbuffer.data(idx_ovl_hg + 190);

    auto ts_xyyzz_yyyz = pbuffer.data(idx_ovl_hg + 191);

    auto ts_xyyzz_yyzz = pbuffer.data(idx_ovl_hg + 192);

    auto ts_xyyzz_yzzz = pbuffer.data(idx_ovl_hg + 193);

    auto ts_xyyzz_zzzz = pbuffer.data(idx_ovl_hg + 194);

    auto ts_xyzzz_xxxx = pbuffer.data(idx_ovl_hg + 195);

    auto ts_xyzzz_xxxy = pbuffer.data(idx_ovl_hg + 196);

    auto ts_xyzzz_xxxz = pbuffer.data(idx_ovl_hg + 197);

    auto ts_xyzzz_xxyy = pbuffer.data(idx_ovl_hg + 198);

    auto ts_xyzzz_xxyz = pbuffer.data(idx_ovl_hg + 199);

    auto ts_xyzzz_xxzz = pbuffer.data(idx_ovl_hg + 200);

    auto ts_xyzzz_xyyy = pbuffer.data(idx_ovl_hg + 201);

    auto ts_xyzzz_xyyz = pbuffer.data(idx_ovl_hg + 202);

    auto ts_xyzzz_xyzz = pbuffer.data(idx_ovl_hg + 203);

    auto ts_xyzzz_xzzz = pbuffer.data(idx_ovl_hg + 204);

    auto ts_xyzzz_yyyy = pbuffer.data(idx_ovl_hg + 205);

    auto ts_xyzzz_yyyz = pbuffer.data(idx_ovl_hg + 206);

    auto ts_xyzzz_yyzz = pbuffer.data(idx_ovl_hg + 207);

    auto ts_xyzzz_yzzz = pbuffer.data(idx_ovl_hg + 208);

    auto ts_xyzzz_zzzz = pbuffer.data(idx_ovl_hg + 209);

    auto ts_xzzzz_xxxx = pbuffer.data(idx_ovl_hg + 210);

    auto ts_xzzzz_xxxy = pbuffer.data(idx_ovl_hg + 211);

    auto ts_xzzzz_xxxz = pbuffer.data(idx_ovl_hg + 212);

    auto ts_xzzzz_xxyy = pbuffer.data(idx_ovl_hg + 213);

    auto ts_xzzzz_xxyz = pbuffer.data(idx_ovl_hg + 214);

    auto ts_xzzzz_xxzz = pbuffer.data(idx_ovl_hg + 215);

    auto ts_xzzzz_xyyy = pbuffer.data(idx_ovl_hg + 216);

    auto ts_xzzzz_xyyz = pbuffer.data(idx_ovl_hg + 217);

    auto ts_xzzzz_xyzz = pbuffer.data(idx_ovl_hg + 218);

    auto ts_xzzzz_xzzz = pbuffer.data(idx_ovl_hg + 219);

    auto ts_xzzzz_yyyy = pbuffer.data(idx_ovl_hg + 220);

    auto ts_xzzzz_yyyz = pbuffer.data(idx_ovl_hg + 221);

    auto ts_xzzzz_yyzz = pbuffer.data(idx_ovl_hg + 222);

    auto ts_xzzzz_yzzz = pbuffer.data(idx_ovl_hg + 223);

    auto ts_xzzzz_zzzz = pbuffer.data(idx_ovl_hg + 224);

    auto ts_yyyyy_xxxx = pbuffer.data(idx_ovl_hg + 225);

    auto ts_yyyyy_xxxy = pbuffer.data(idx_ovl_hg + 226);

    auto ts_yyyyy_xxxz = pbuffer.data(idx_ovl_hg + 227);

    auto ts_yyyyy_xxyy = pbuffer.data(idx_ovl_hg + 228);

    auto ts_yyyyy_xxyz = pbuffer.data(idx_ovl_hg + 229);

    auto ts_yyyyy_xxzz = pbuffer.data(idx_ovl_hg + 230);

    auto ts_yyyyy_xyyy = pbuffer.data(idx_ovl_hg + 231);

    auto ts_yyyyy_xyyz = pbuffer.data(idx_ovl_hg + 232);

    auto ts_yyyyy_xyzz = pbuffer.data(idx_ovl_hg + 233);

    auto ts_yyyyy_xzzz = pbuffer.data(idx_ovl_hg + 234);

    auto ts_yyyyy_yyyy = pbuffer.data(idx_ovl_hg + 235);

    auto ts_yyyyy_yyyz = pbuffer.data(idx_ovl_hg + 236);

    auto ts_yyyyy_yyzz = pbuffer.data(idx_ovl_hg + 237);

    auto ts_yyyyy_yzzz = pbuffer.data(idx_ovl_hg + 238);

    auto ts_yyyyy_zzzz = pbuffer.data(idx_ovl_hg + 239);

    auto ts_yyyyz_xxxx = pbuffer.data(idx_ovl_hg + 240);

    auto ts_yyyyz_xxxy = pbuffer.data(idx_ovl_hg + 241);

    auto ts_yyyyz_xxxz = pbuffer.data(idx_ovl_hg + 242);

    auto ts_yyyyz_xxyy = pbuffer.data(idx_ovl_hg + 243);

    auto ts_yyyyz_xxyz = pbuffer.data(idx_ovl_hg + 244);

    auto ts_yyyyz_xxzz = pbuffer.data(idx_ovl_hg + 245);

    auto ts_yyyyz_xyyy = pbuffer.data(idx_ovl_hg + 246);

    auto ts_yyyyz_xyyz = pbuffer.data(idx_ovl_hg + 247);

    auto ts_yyyyz_xyzz = pbuffer.data(idx_ovl_hg + 248);

    auto ts_yyyyz_xzzz = pbuffer.data(idx_ovl_hg + 249);

    auto ts_yyyyz_yyyy = pbuffer.data(idx_ovl_hg + 250);

    auto ts_yyyyz_yyyz = pbuffer.data(idx_ovl_hg + 251);

    auto ts_yyyyz_yyzz = pbuffer.data(idx_ovl_hg + 252);

    auto ts_yyyyz_yzzz = pbuffer.data(idx_ovl_hg + 253);

    auto ts_yyyyz_zzzz = pbuffer.data(idx_ovl_hg + 254);

    auto ts_yyyzz_xxxx = pbuffer.data(idx_ovl_hg + 255);

    auto ts_yyyzz_xxxy = pbuffer.data(idx_ovl_hg + 256);

    auto ts_yyyzz_xxxz = pbuffer.data(idx_ovl_hg + 257);

    auto ts_yyyzz_xxyy = pbuffer.data(idx_ovl_hg + 258);

    auto ts_yyyzz_xxyz = pbuffer.data(idx_ovl_hg + 259);

    auto ts_yyyzz_xxzz = pbuffer.data(idx_ovl_hg + 260);

    auto ts_yyyzz_xyyy = pbuffer.data(idx_ovl_hg + 261);

    auto ts_yyyzz_xyyz = pbuffer.data(idx_ovl_hg + 262);

    auto ts_yyyzz_xyzz = pbuffer.data(idx_ovl_hg + 263);

    auto ts_yyyzz_xzzz = pbuffer.data(idx_ovl_hg + 264);

    auto ts_yyyzz_yyyy = pbuffer.data(idx_ovl_hg + 265);

    auto ts_yyyzz_yyyz = pbuffer.data(idx_ovl_hg + 266);

    auto ts_yyyzz_yyzz = pbuffer.data(idx_ovl_hg + 267);

    auto ts_yyyzz_yzzz = pbuffer.data(idx_ovl_hg + 268);

    auto ts_yyyzz_zzzz = pbuffer.data(idx_ovl_hg + 269);

    auto ts_yyzzz_xxxx = pbuffer.data(idx_ovl_hg + 270);

    auto ts_yyzzz_xxxy = pbuffer.data(idx_ovl_hg + 271);

    auto ts_yyzzz_xxxz = pbuffer.data(idx_ovl_hg + 272);

    auto ts_yyzzz_xxyy = pbuffer.data(idx_ovl_hg + 273);

    auto ts_yyzzz_xxyz = pbuffer.data(idx_ovl_hg + 274);

    auto ts_yyzzz_xxzz = pbuffer.data(idx_ovl_hg + 275);

    auto ts_yyzzz_xyyy = pbuffer.data(idx_ovl_hg + 276);

    auto ts_yyzzz_xyyz = pbuffer.data(idx_ovl_hg + 277);

    auto ts_yyzzz_xyzz = pbuffer.data(idx_ovl_hg + 278);

    auto ts_yyzzz_xzzz = pbuffer.data(idx_ovl_hg + 279);

    auto ts_yyzzz_yyyy = pbuffer.data(idx_ovl_hg + 280);

    auto ts_yyzzz_yyyz = pbuffer.data(idx_ovl_hg + 281);

    auto ts_yyzzz_yyzz = pbuffer.data(idx_ovl_hg + 282);

    auto ts_yyzzz_yzzz = pbuffer.data(idx_ovl_hg + 283);

    auto ts_yyzzz_zzzz = pbuffer.data(idx_ovl_hg + 284);

    auto ts_yzzzz_xxxx = pbuffer.data(idx_ovl_hg + 285);

    auto ts_yzzzz_xxxy = pbuffer.data(idx_ovl_hg + 286);

    auto ts_yzzzz_xxxz = pbuffer.data(idx_ovl_hg + 287);

    auto ts_yzzzz_xxyy = pbuffer.data(idx_ovl_hg + 288);

    auto ts_yzzzz_xxyz = pbuffer.data(idx_ovl_hg + 289);

    auto ts_yzzzz_xxzz = pbuffer.data(idx_ovl_hg + 290);

    auto ts_yzzzz_xyyy = pbuffer.data(idx_ovl_hg + 291);

    auto ts_yzzzz_xyyz = pbuffer.data(idx_ovl_hg + 292);

    auto ts_yzzzz_xyzz = pbuffer.data(idx_ovl_hg + 293);

    auto ts_yzzzz_xzzz = pbuffer.data(idx_ovl_hg + 294);

    auto ts_yzzzz_yyyy = pbuffer.data(idx_ovl_hg + 295);

    auto ts_yzzzz_yyyz = pbuffer.data(idx_ovl_hg + 296);

    auto ts_yzzzz_yyzz = pbuffer.data(idx_ovl_hg + 297);

    auto ts_yzzzz_yzzz = pbuffer.data(idx_ovl_hg + 298);

    auto ts_yzzzz_zzzz = pbuffer.data(idx_ovl_hg + 299);

    auto ts_zzzzz_xxxx = pbuffer.data(idx_ovl_hg + 300);

    auto ts_zzzzz_xxxy = pbuffer.data(idx_ovl_hg + 301);

    auto ts_zzzzz_xxxz = pbuffer.data(idx_ovl_hg + 302);

    auto ts_zzzzz_xxyy = pbuffer.data(idx_ovl_hg + 303);

    auto ts_zzzzz_xxyz = pbuffer.data(idx_ovl_hg + 304);

    auto ts_zzzzz_xxzz = pbuffer.data(idx_ovl_hg + 305);

    auto ts_zzzzz_xyyy = pbuffer.data(idx_ovl_hg + 306);

    auto ts_zzzzz_xyyz = pbuffer.data(idx_ovl_hg + 307);

    auto ts_zzzzz_xyzz = pbuffer.data(idx_ovl_hg + 308);

    auto ts_zzzzz_xzzz = pbuffer.data(idx_ovl_hg + 309);

    auto ts_zzzzz_yyyy = pbuffer.data(idx_ovl_hg + 310);

    auto ts_zzzzz_yyyz = pbuffer.data(idx_ovl_hg + 311);

    auto ts_zzzzz_yyzz = pbuffer.data(idx_ovl_hg + 312);

    auto ts_zzzzz_yzzz = pbuffer.data(idx_ovl_hg + 313);

    auto ts_zzzzz_zzzz = pbuffer.data(idx_ovl_hg + 314);

    // Set up 0-15 components of targeted buffer : HG

    auto tk_xxxxx_xxxx = pbuffer.data(idx_kin_hg);

    auto tk_xxxxx_xxxy = pbuffer.data(idx_kin_hg + 1);

    auto tk_xxxxx_xxxz = pbuffer.data(idx_kin_hg + 2);

    auto tk_xxxxx_xxyy = pbuffer.data(idx_kin_hg + 3);

    auto tk_xxxxx_xxyz = pbuffer.data(idx_kin_hg + 4);

    auto tk_xxxxx_xxzz = pbuffer.data(idx_kin_hg + 5);

    auto tk_xxxxx_xyyy = pbuffer.data(idx_kin_hg + 6);

    auto tk_xxxxx_xyyz = pbuffer.data(idx_kin_hg + 7);

    auto tk_xxxxx_xyzz = pbuffer.data(idx_kin_hg + 8);

    auto tk_xxxxx_xzzz = pbuffer.data(idx_kin_hg + 9);

    auto tk_xxxxx_yyyy = pbuffer.data(idx_kin_hg + 10);

    auto tk_xxxxx_yyyz = pbuffer.data(idx_kin_hg + 11);

    auto tk_xxxxx_yyzz = pbuffer.data(idx_kin_hg + 12);

    auto tk_xxxxx_yzzz = pbuffer.data(idx_kin_hg + 13);

    auto tk_xxxxx_zzzz = pbuffer.data(idx_kin_hg + 14);

#pragma omp simd aligned(pa_x,              \
                             tk_xxx_xxxx,   \
                             tk_xxx_xxxy,   \
                             tk_xxx_xxxz,   \
                             tk_xxx_xxyy,   \
                             tk_xxx_xxyz,   \
                             tk_xxx_xxzz,   \
                             tk_xxx_xyyy,   \
                             tk_xxx_xyyz,   \
                             tk_xxx_xyzz,   \
                             tk_xxx_xzzz,   \
                             tk_xxx_yyyy,   \
                             tk_xxx_yyyz,   \
                             tk_xxx_yyzz,   \
                             tk_xxx_yzzz,   \
                             tk_xxx_zzzz,   \
                             tk_xxxx_xxx,   \
                             tk_xxxx_xxxx,  \
                             tk_xxxx_xxxy,  \
                             tk_xxxx_xxxz,  \
                             tk_xxxx_xxy,   \
                             tk_xxxx_xxyy,  \
                             tk_xxxx_xxyz,  \
                             tk_xxxx_xxz,   \
                             tk_xxxx_xxzz,  \
                             tk_xxxx_xyy,   \
                             tk_xxxx_xyyy,  \
                             tk_xxxx_xyyz,  \
                             tk_xxxx_xyz,   \
                             tk_xxxx_xyzz,  \
                             tk_xxxx_xzz,   \
                             tk_xxxx_xzzz,  \
                             tk_xxxx_yyy,   \
                             tk_xxxx_yyyy,  \
                             tk_xxxx_yyyz,  \
                             tk_xxxx_yyz,   \
                             tk_xxxx_yyzz,  \
                             tk_xxxx_yzz,   \
                             tk_xxxx_yzzz,  \
                             tk_xxxx_zzz,   \
                             tk_xxxx_zzzz,  \
                             tk_xxxxx_xxxx, \
                             tk_xxxxx_xxxy, \
                             tk_xxxxx_xxxz, \
                             tk_xxxxx_xxyy, \
                             tk_xxxxx_xxyz, \
                             tk_xxxxx_xxzz, \
                             tk_xxxxx_xyyy, \
                             tk_xxxxx_xyyz, \
                             tk_xxxxx_xyzz, \
                             tk_xxxxx_xzzz, \
                             tk_xxxxx_yyyy, \
                             tk_xxxxx_yyyz, \
                             tk_xxxxx_yyzz, \
                             tk_xxxxx_yzzz, \
                             tk_xxxxx_zzzz, \
                             ts_xxx_xxxx,   \
                             ts_xxx_xxxy,   \
                             ts_xxx_xxxz,   \
                             ts_xxx_xxyy,   \
                             ts_xxx_xxyz,   \
                             ts_xxx_xxzz,   \
                             ts_xxx_xyyy,   \
                             ts_xxx_xyyz,   \
                             ts_xxx_xyzz,   \
                             ts_xxx_xzzz,   \
                             ts_xxx_yyyy,   \
                             ts_xxx_yyyz,   \
                             ts_xxx_yyzz,   \
                             ts_xxx_yzzz,   \
                             ts_xxx_zzzz,   \
                             ts_xxxxx_xxxx, \
                             ts_xxxxx_xxxy, \
                             ts_xxxxx_xxxz, \
                             ts_xxxxx_xxyy, \
                             ts_xxxxx_xxyz, \
                             ts_xxxxx_xxzz, \
                             ts_xxxxx_xyyy, \
                             ts_xxxxx_xyyz, \
                             ts_xxxxx_xyzz, \
                             ts_xxxxx_xzzz, \
                             ts_xxxxx_yyyy, \
                             ts_xxxxx_yyyz, \
                             ts_xxxxx_yyzz, \
                             ts_xxxxx_yzzz, \
                             ts_xxxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxx_xxxx[i] = -8.0 * ts_xxx_xxxx[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxx[i] * fe_0 + 4.0 * tk_xxxx_xxx[i] * fe_0 +
                           tk_xxxx_xxxx[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxx[i] * fz_0;

        tk_xxxxx_xxxy[i] = -8.0 * ts_xxx_xxxy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxy[i] * fe_0 + 3.0 * tk_xxxx_xxy[i] * fe_0 +
                           tk_xxxx_xxxy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxy[i] * fz_0;

        tk_xxxxx_xxxz[i] = -8.0 * ts_xxx_xxxz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxz[i] * fe_0 + 3.0 * tk_xxxx_xxz[i] * fe_0 +
                           tk_xxxx_xxxz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxz[i] * fz_0;

        tk_xxxxx_xxyy[i] = -8.0 * ts_xxx_xxyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyy[i] * fe_0 + 2.0 * tk_xxxx_xyy[i] * fe_0 +
                           tk_xxxx_xxyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyy[i] * fz_0;

        tk_xxxxx_xxyz[i] = -8.0 * ts_xxx_xxyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyz[i] * fe_0 + 2.0 * tk_xxxx_xyz[i] * fe_0 +
                           tk_xxxx_xxyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyz[i] * fz_0;

        tk_xxxxx_xxzz[i] = -8.0 * ts_xxx_xxzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxzz[i] * fe_0 + 2.0 * tk_xxxx_xzz[i] * fe_0 +
                           tk_xxxx_xxzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxzz[i] * fz_0;

        tk_xxxxx_xyyy[i] = -8.0 * ts_xxx_xyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyy[i] * fe_0 + tk_xxxx_yyy[i] * fe_0 + tk_xxxx_xyyy[i] * pa_x[i] +
                           2.0 * ts_xxxxx_xyyy[i] * fz_0;

        tk_xxxxx_xyyz[i] = -8.0 * ts_xxx_xyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyz[i] * fe_0 + tk_xxxx_yyz[i] * fe_0 + tk_xxxx_xyyz[i] * pa_x[i] +
                           2.0 * ts_xxxxx_xyyz[i] * fz_0;

        tk_xxxxx_xyzz[i] = -8.0 * ts_xxx_xyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyzz[i] * fe_0 + tk_xxxx_yzz[i] * fe_0 + tk_xxxx_xyzz[i] * pa_x[i] +
                           2.0 * ts_xxxxx_xyzz[i] * fz_0;

        tk_xxxxx_xzzz[i] = -8.0 * ts_xxx_xzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xzzz[i] * fe_0 + tk_xxxx_zzz[i] * fe_0 + tk_xxxx_xzzz[i] * pa_x[i] +
                           2.0 * ts_xxxxx_xzzz[i] * fz_0;

        tk_xxxxx_yyyy[i] =
            -8.0 * ts_xxx_yyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyy[i] * fe_0 + tk_xxxx_yyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyy[i] * fz_0;

        tk_xxxxx_yyyz[i] =
            -8.0 * ts_xxx_yyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyz[i] * fe_0 + tk_xxxx_yyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyz[i] * fz_0;

        tk_xxxxx_yyzz[i] =
            -8.0 * ts_xxx_yyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyzz[i] * fe_0 + tk_xxxx_yyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyzz[i] * fz_0;

        tk_xxxxx_yzzz[i] =
            -8.0 * ts_xxx_yzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yzzz[i] * fe_0 + tk_xxxx_yzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yzzz[i] * fz_0;

        tk_xxxxx_zzzz[i] =
            -8.0 * ts_xxx_zzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_zzzz[i] * fe_0 + tk_xxxx_zzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_zzzz[i] * fz_0;
    }

    // Set up 15-30 components of targeted buffer : HG

    auto tk_xxxxy_xxxx = pbuffer.data(idx_kin_hg + 15);

    auto tk_xxxxy_xxxy = pbuffer.data(idx_kin_hg + 16);

    auto tk_xxxxy_xxxz = pbuffer.data(idx_kin_hg + 17);

    auto tk_xxxxy_xxyy = pbuffer.data(idx_kin_hg + 18);

    auto tk_xxxxy_xxyz = pbuffer.data(idx_kin_hg + 19);

    auto tk_xxxxy_xxzz = pbuffer.data(idx_kin_hg + 20);

    auto tk_xxxxy_xyyy = pbuffer.data(idx_kin_hg + 21);

    auto tk_xxxxy_xyyz = pbuffer.data(idx_kin_hg + 22);

    auto tk_xxxxy_xyzz = pbuffer.data(idx_kin_hg + 23);

    auto tk_xxxxy_xzzz = pbuffer.data(idx_kin_hg + 24);

    auto tk_xxxxy_yyyy = pbuffer.data(idx_kin_hg + 25);

    auto tk_xxxxy_yyyz = pbuffer.data(idx_kin_hg + 26);

    auto tk_xxxxy_yyzz = pbuffer.data(idx_kin_hg + 27);

    auto tk_xxxxy_yzzz = pbuffer.data(idx_kin_hg + 28);

    auto tk_xxxxy_zzzz = pbuffer.data(idx_kin_hg + 29);

#pragma omp simd aligned(pa_y,              \
                             tk_xxxx_xxx,   \
                             tk_xxxx_xxxx,  \
                             tk_xxxx_xxxy,  \
                             tk_xxxx_xxxz,  \
                             tk_xxxx_xxy,   \
                             tk_xxxx_xxyy,  \
                             tk_xxxx_xxyz,  \
                             tk_xxxx_xxz,   \
                             tk_xxxx_xxzz,  \
                             tk_xxxx_xyy,   \
                             tk_xxxx_xyyy,  \
                             tk_xxxx_xyyz,  \
                             tk_xxxx_xyz,   \
                             tk_xxxx_xyzz,  \
                             tk_xxxx_xzz,   \
                             tk_xxxx_xzzz,  \
                             tk_xxxx_yyy,   \
                             tk_xxxx_yyyy,  \
                             tk_xxxx_yyyz,  \
                             tk_xxxx_yyz,   \
                             tk_xxxx_yyzz,  \
                             tk_xxxx_yzz,   \
                             tk_xxxx_yzzz,  \
                             tk_xxxx_zzz,   \
                             tk_xxxx_zzzz,  \
                             tk_xxxxy_xxxx, \
                             tk_xxxxy_xxxy, \
                             tk_xxxxy_xxxz, \
                             tk_xxxxy_xxyy, \
                             tk_xxxxy_xxyz, \
                             tk_xxxxy_xxzz, \
                             tk_xxxxy_xyyy, \
                             tk_xxxxy_xyyz, \
                             tk_xxxxy_xyzz, \
                             tk_xxxxy_xzzz, \
                             tk_xxxxy_yyyy, \
                             tk_xxxxy_yyyz, \
                             tk_xxxxy_yyzz, \
                             tk_xxxxy_yzzz, \
                             tk_xxxxy_zzzz, \
                             ts_xxxxy_xxxx, \
                             ts_xxxxy_xxxy, \
                             ts_xxxxy_xxxz, \
                             ts_xxxxy_xxyy, \
                             ts_xxxxy_xxyz, \
                             ts_xxxxy_xxzz, \
                             ts_xxxxy_xyyy, \
                             ts_xxxxy_xyyz, \
                             ts_xxxxy_xyzz, \
                             ts_xxxxy_xzzz, \
                             ts_xxxxy_yyyy, \
                             ts_xxxxy_yyyz, \
                             ts_xxxxy_yyzz, \
                             ts_xxxxy_yzzz, \
                             ts_xxxxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxy_xxxx[i] = tk_xxxx_xxxx[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxx[i] * fz_0;

        tk_xxxxy_xxxy[i] = tk_xxxx_xxx[i] * fe_0 + tk_xxxx_xxxy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxy[i] * fz_0;

        tk_xxxxy_xxxz[i] = tk_xxxx_xxxz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxz[i] * fz_0;

        tk_xxxxy_xxyy[i] = 2.0 * tk_xxxx_xxy[i] * fe_0 + tk_xxxx_xxyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyy[i] * fz_0;

        tk_xxxxy_xxyz[i] = tk_xxxx_xxz[i] * fe_0 + tk_xxxx_xxyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyz[i] * fz_0;

        tk_xxxxy_xxzz[i] = tk_xxxx_xxzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxzz[i] * fz_0;

        tk_xxxxy_xyyy[i] = 3.0 * tk_xxxx_xyy[i] * fe_0 + tk_xxxx_xyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyy[i] * fz_0;

        tk_xxxxy_xyyz[i] = 2.0 * tk_xxxx_xyz[i] * fe_0 + tk_xxxx_xyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyz[i] * fz_0;

        tk_xxxxy_xyzz[i] = tk_xxxx_xzz[i] * fe_0 + tk_xxxx_xyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyzz[i] * fz_0;

        tk_xxxxy_xzzz[i] = tk_xxxx_xzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xzzz[i] * fz_0;

        tk_xxxxy_yyyy[i] = 4.0 * tk_xxxx_yyy[i] * fe_0 + tk_xxxx_yyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyy[i] * fz_0;

        tk_xxxxy_yyyz[i] = 3.0 * tk_xxxx_yyz[i] * fe_0 + tk_xxxx_yyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyz[i] * fz_0;

        tk_xxxxy_yyzz[i] = 2.0 * tk_xxxx_yzz[i] * fe_0 + tk_xxxx_yyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyzz[i] * fz_0;

        tk_xxxxy_yzzz[i] = tk_xxxx_zzz[i] * fe_0 + tk_xxxx_yzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yzzz[i] * fz_0;

        tk_xxxxy_zzzz[i] = tk_xxxx_zzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_zzzz[i] * fz_0;
    }

    // Set up 30-45 components of targeted buffer : HG

    auto tk_xxxxz_xxxx = pbuffer.data(idx_kin_hg + 30);

    auto tk_xxxxz_xxxy = pbuffer.data(idx_kin_hg + 31);

    auto tk_xxxxz_xxxz = pbuffer.data(idx_kin_hg + 32);

    auto tk_xxxxz_xxyy = pbuffer.data(idx_kin_hg + 33);

    auto tk_xxxxz_xxyz = pbuffer.data(idx_kin_hg + 34);

    auto tk_xxxxz_xxzz = pbuffer.data(idx_kin_hg + 35);

    auto tk_xxxxz_xyyy = pbuffer.data(idx_kin_hg + 36);

    auto tk_xxxxz_xyyz = pbuffer.data(idx_kin_hg + 37);

    auto tk_xxxxz_xyzz = pbuffer.data(idx_kin_hg + 38);

    auto tk_xxxxz_xzzz = pbuffer.data(idx_kin_hg + 39);

    auto tk_xxxxz_yyyy = pbuffer.data(idx_kin_hg + 40);

    auto tk_xxxxz_yyyz = pbuffer.data(idx_kin_hg + 41);

    auto tk_xxxxz_yyzz = pbuffer.data(idx_kin_hg + 42);

    auto tk_xxxxz_yzzz = pbuffer.data(idx_kin_hg + 43);

    auto tk_xxxxz_zzzz = pbuffer.data(idx_kin_hg + 44);

#pragma omp simd aligned(pa_z,              \
                             tk_xxxx_xxx,   \
                             tk_xxxx_xxxx,  \
                             tk_xxxx_xxxy,  \
                             tk_xxxx_xxxz,  \
                             tk_xxxx_xxy,   \
                             tk_xxxx_xxyy,  \
                             tk_xxxx_xxyz,  \
                             tk_xxxx_xxz,   \
                             tk_xxxx_xxzz,  \
                             tk_xxxx_xyy,   \
                             tk_xxxx_xyyy,  \
                             tk_xxxx_xyyz,  \
                             tk_xxxx_xyz,   \
                             tk_xxxx_xyzz,  \
                             tk_xxxx_xzz,   \
                             tk_xxxx_xzzz,  \
                             tk_xxxx_yyy,   \
                             tk_xxxx_yyyy,  \
                             tk_xxxx_yyyz,  \
                             tk_xxxx_yyz,   \
                             tk_xxxx_yyzz,  \
                             tk_xxxx_yzz,   \
                             tk_xxxx_yzzz,  \
                             tk_xxxx_zzz,   \
                             tk_xxxx_zzzz,  \
                             tk_xxxxz_xxxx, \
                             tk_xxxxz_xxxy, \
                             tk_xxxxz_xxxz, \
                             tk_xxxxz_xxyy, \
                             tk_xxxxz_xxyz, \
                             tk_xxxxz_xxzz, \
                             tk_xxxxz_xyyy, \
                             tk_xxxxz_xyyz, \
                             tk_xxxxz_xyzz, \
                             tk_xxxxz_xzzz, \
                             tk_xxxxz_yyyy, \
                             tk_xxxxz_yyyz, \
                             tk_xxxxz_yyzz, \
                             tk_xxxxz_yzzz, \
                             tk_xxxxz_zzzz, \
                             ts_xxxxz_xxxx, \
                             ts_xxxxz_xxxy, \
                             ts_xxxxz_xxxz, \
                             ts_xxxxz_xxyy, \
                             ts_xxxxz_xxyz, \
                             ts_xxxxz_xxzz, \
                             ts_xxxxz_xyyy, \
                             ts_xxxxz_xyyz, \
                             ts_xxxxz_xyzz, \
                             ts_xxxxz_xzzz, \
                             ts_xxxxz_yyyy, \
                             ts_xxxxz_yyyz, \
                             ts_xxxxz_yyzz, \
                             ts_xxxxz_yzzz, \
                             ts_xxxxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxz_xxxx[i] = tk_xxxx_xxxx[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxx[i] * fz_0;

        tk_xxxxz_xxxy[i] = tk_xxxx_xxxy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxy[i] * fz_0;

        tk_xxxxz_xxxz[i] = tk_xxxx_xxx[i] * fe_0 + tk_xxxx_xxxz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxz[i] * fz_0;

        tk_xxxxz_xxyy[i] = tk_xxxx_xxyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyy[i] * fz_0;

        tk_xxxxz_xxyz[i] = tk_xxxx_xxy[i] * fe_0 + tk_xxxx_xxyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyz[i] * fz_0;

        tk_xxxxz_xxzz[i] = 2.0 * tk_xxxx_xxz[i] * fe_0 + tk_xxxx_xxzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxzz[i] * fz_0;

        tk_xxxxz_xyyy[i] = tk_xxxx_xyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyy[i] * fz_0;

        tk_xxxxz_xyyz[i] = tk_xxxx_xyy[i] * fe_0 + tk_xxxx_xyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyz[i] * fz_0;

        tk_xxxxz_xyzz[i] = 2.0 * tk_xxxx_xyz[i] * fe_0 + tk_xxxx_xyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyzz[i] * fz_0;

        tk_xxxxz_xzzz[i] = 3.0 * tk_xxxx_xzz[i] * fe_0 + tk_xxxx_xzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xzzz[i] * fz_0;

        tk_xxxxz_yyyy[i] = tk_xxxx_yyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyy[i] * fz_0;

        tk_xxxxz_yyyz[i] = tk_xxxx_yyy[i] * fe_0 + tk_xxxx_yyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyz[i] * fz_0;

        tk_xxxxz_yyzz[i] = 2.0 * tk_xxxx_yyz[i] * fe_0 + tk_xxxx_yyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyzz[i] * fz_0;

        tk_xxxxz_yzzz[i] = 3.0 * tk_xxxx_yzz[i] * fe_0 + tk_xxxx_yzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yzzz[i] * fz_0;

        tk_xxxxz_zzzz[i] = 4.0 * tk_xxxx_zzz[i] * fe_0 + tk_xxxx_zzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_zzzz[i] * fz_0;
    }

    // Set up 45-60 components of targeted buffer : HG

    auto tk_xxxyy_xxxx = pbuffer.data(idx_kin_hg + 45);

    auto tk_xxxyy_xxxy = pbuffer.data(idx_kin_hg + 46);

    auto tk_xxxyy_xxxz = pbuffer.data(idx_kin_hg + 47);

    auto tk_xxxyy_xxyy = pbuffer.data(idx_kin_hg + 48);

    auto tk_xxxyy_xxyz = pbuffer.data(idx_kin_hg + 49);

    auto tk_xxxyy_xxzz = pbuffer.data(idx_kin_hg + 50);

    auto tk_xxxyy_xyyy = pbuffer.data(idx_kin_hg + 51);

    auto tk_xxxyy_xyyz = pbuffer.data(idx_kin_hg + 52);

    auto tk_xxxyy_xyzz = pbuffer.data(idx_kin_hg + 53);

    auto tk_xxxyy_xzzz = pbuffer.data(idx_kin_hg + 54);

    auto tk_xxxyy_yyyy = pbuffer.data(idx_kin_hg + 55);

    auto tk_xxxyy_yyyz = pbuffer.data(idx_kin_hg + 56);

    auto tk_xxxyy_yyzz = pbuffer.data(idx_kin_hg + 57);

    auto tk_xxxyy_yzzz = pbuffer.data(idx_kin_hg + 58);

    auto tk_xxxyy_zzzz = pbuffer.data(idx_kin_hg + 59);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xxx_xxxx,   \
                             tk_xxx_xxxz,   \
                             tk_xxx_xxzz,   \
                             tk_xxx_xzzz,   \
                             tk_xxxy_xxxx,  \
                             tk_xxxy_xxxz,  \
                             tk_xxxy_xxzz,  \
                             tk_xxxy_xzzz,  \
                             tk_xxxyy_xxxx, \
                             tk_xxxyy_xxxy, \
                             tk_xxxyy_xxxz, \
                             tk_xxxyy_xxyy, \
                             tk_xxxyy_xxyz, \
                             tk_xxxyy_xxzz, \
                             tk_xxxyy_xyyy, \
                             tk_xxxyy_xyyz, \
                             tk_xxxyy_xyzz, \
                             tk_xxxyy_xzzz, \
                             tk_xxxyy_yyyy, \
                             tk_xxxyy_yyyz, \
                             tk_xxxyy_yyzz, \
                             tk_xxxyy_yzzz, \
                             tk_xxxyy_zzzz, \
                             tk_xxyy_xxxy,  \
                             tk_xxyy_xxy,   \
                             tk_xxyy_xxyy,  \
                             tk_xxyy_xxyz,  \
                             tk_xxyy_xyy,   \
                             tk_xxyy_xyyy,  \
                             tk_xxyy_xyyz,  \
                             tk_xxyy_xyz,   \
                             tk_xxyy_xyzz,  \
                             tk_xxyy_yyy,   \
                             tk_xxyy_yyyy,  \
                             tk_xxyy_yyyz,  \
                             tk_xxyy_yyz,   \
                             tk_xxyy_yyzz,  \
                             tk_xxyy_yzz,   \
                             tk_xxyy_yzzz,  \
                             tk_xxyy_zzzz,  \
                             tk_xyy_xxxy,   \
                             tk_xyy_xxyy,   \
                             tk_xyy_xxyz,   \
                             tk_xyy_xyyy,   \
                             tk_xyy_xyyz,   \
                             tk_xyy_xyzz,   \
                             tk_xyy_yyyy,   \
                             tk_xyy_yyyz,   \
                             tk_xyy_yyzz,   \
                             tk_xyy_yzzz,   \
                             tk_xyy_zzzz,   \
                             ts_xxx_xxxx,   \
                             ts_xxx_xxxz,   \
                             ts_xxx_xxzz,   \
                             ts_xxx_xzzz,   \
                             ts_xxxyy_xxxx, \
                             ts_xxxyy_xxxy, \
                             ts_xxxyy_xxxz, \
                             ts_xxxyy_xxyy, \
                             ts_xxxyy_xxyz, \
                             ts_xxxyy_xxzz, \
                             ts_xxxyy_xyyy, \
                             ts_xxxyy_xyyz, \
                             ts_xxxyy_xyzz, \
                             ts_xxxyy_xzzz, \
                             ts_xxxyy_yyyy, \
                             ts_xxxyy_yyyz, \
                             ts_xxxyy_yyzz, \
                             ts_xxxyy_yzzz, \
                             ts_xxxyy_zzzz, \
                             ts_xyy_xxxy,   \
                             ts_xyy_xxyy,   \
                             ts_xyy_xxyz,   \
                             ts_xyy_xyyy,   \
                             ts_xyy_xyyz,   \
                             ts_xyy_xyzz,   \
                             ts_xyy_yyyy,   \
                             ts_xyy_yyyz,   \
                             ts_xyy_yyzz,   \
                             ts_xyy_yzzz,   \
                             ts_xyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyy_xxxx[i] = -2.0 * ts_xxx_xxxx[i] * fbe_0 * fz_0 + tk_xxx_xxxx[i] * fe_0 + tk_xxxy_xxxx[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxx[i] * fz_0;

        tk_xxxyy_xxxy[i] = -4.0 * ts_xyy_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxy[i] * fe_0 + 3.0 * tk_xxyy_xxy[i] * fe_0 +
                           tk_xxyy_xxxy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxy[i] * fz_0;

        tk_xxxyy_xxxz[i] = -2.0 * ts_xxx_xxxz[i] * fbe_0 * fz_0 + tk_xxx_xxxz[i] * fe_0 + tk_xxxy_xxxz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxz[i] * fz_0;

        tk_xxxyy_xxyy[i] = -4.0 * ts_xyy_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyy[i] * fe_0 + 2.0 * tk_xxyy_xyy[i] * fe_0 +
                           tk_xxyy_xxyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyy[i] * fz_0;

        tk_xxxyy_xxyz[i] = -4.0 * ts_xyy_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyz[i] * fe_0 + 2.0 * tk_xxyy_xyz[i] * fe_0 +
                           tk_xxyy_xxyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyz[i] * fz_0;

        tk_xxxyy_xxzz[i] = -2.0 * ts_xxx_xxzz[i] * fbe_0 * fz_0 + tk_xxx_xxzz[i] * fe_0 + tk_xxxy_xxzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxzz[i] * fz_0;

        tk_xxxyy_xyyy[i] = -4.0 * ts_xyy_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyy[i] * fe_0 + tk_xxyy_yyy[i] * fe_0 + tk_xxyy_xyyy[i] * pa_x[i] +
                           2.0 * ts_xxxyy_xyyy[i] * fz_0;

        tk_xxxyy_xyyz[i] = -4.0 * ts_xyy_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyz[i] * fe_0 + tk_xxyy_yyz[i] * fe_0 + tk_xxyy_xyyz[i] * pa_x[i] +
                           2.0 * ts_xxxyy_xyyz[i] * fz_0;

        tk_xxxyy_xyzz[i] = -4.0 * ts_xyy_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyzz[i] * fe_0 + tk_xxyy_yzz[i] * fe_0 + tk_xxyy_xyzz[i] * pa_x[i] +
                           2.0 * ts_xxxyy_xyzz[i] * fz_0;

        tk_xxxyy_xzzz[i] = -2.0 * ts_xxx_xzzz[i] * fbe_0 * fz_0 + tk_xxx_xzzz[i] * fe_0 + tk_xxxy_xzzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xzzz[i] * fz_0;

        tk_xxxyy_yyyy[i] =
            -4.0 * ts_xyy_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyy[i] * fe_0 + tk_xxyy_yyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyy[i] * fz_0;

        tk_xxxyy_yyyz[i] =
            -4.0 * ts_xyy_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyz[i] * fe_0 + tk_xxyy_yyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyz[i] * fz_0;

        tk_xxxyy_yyzz[i] =
            -4.0 * ts_xyy_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyzz[i] * fe_0 + tk_xxyy_yyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyzz[i] * fz_0;

        tk_xxxyy_yzzz[i] =
            -4.0 * ts_xyy_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yzzz[i] * fe_0 + tk_xxyy_yzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yzzz[i] * fz_0;

        tk_xxxyy_zzzz[i] =
            -4.0 * ts_xyy_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_zzzz[i] * fe_0 + tk_xxyy_zzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_zzzz[i] * fz_0;
    }

    // Set up 60-75 components of targeted buffer : HG

    auto tk_xxxyz_xxxx = pbuffer.data(idx_kin_hg + 60);

    auto tk_xxxyz_xxxy = pbuffer.data(idx_kin_hg + 61);

    auto tk_xxxyz_xxxz = pbuffer.data(idx_kin_hg + 62);

    auto tk_xxxyz_xxyy = pbuffer.data(idx_kin_hg + 63);

    auto tk_xxxyz_xxyz = pbuffer.data(idx_kin_hg + 64);

    auto tk_xxxyz_xxzz = pbuffer.data(idx_kin_hg + 65);

    auto tk_xxxyz_xyyy = pbuffer.data(idx_kin_hg + 66);

    auto tk_xxxyz_xyyz = pbuffer.data(idx_kin_hg + 67);

    auto tk_xxxyz_xyzz = pbuffer.data(idx_kin_hg + 68);

    auto tk_xxxyz_xzzz = pbuffer.data(idx_kin_hg + 69);

    auto tk_xxxyz_yyyy = pbuffer.data(idx_kin_hg + 70);

    auto tk_xxxyz_yyyz = pbuffer.data(idx_kin_hg + 71);

    auto tk_xxxyz_yyzz = pbuffer.data(idx_kin_hg + 72);

    auto tk_xxxyz_yzzz = pbuffer.data(idx_kin_hg + 73);

    auto tk_xxxyz_zzzz = pbuffer.data(idx_kin_hg + 74);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_xxxy_xxxy,  \
                             tk_xxxy_xxyy,  \
                             tk_xxxy_xyyy,  \
                             tk_xxxy_yyyy,  \
                             tk_xxxyz_xxxx, \
                             tk_xxxyz_xxxy, \
                             tk_xxxyz_xxxz, \
                             tk_xxxyz_xxyy, \
                             tk_xxxyz_xxyz, \
                             tk_xxxyz_xxzz, \
                             tk_xxxyz_xyyy, \
                             tk_xxxyz_xyyz, \
                             tk_xxxyz_xyzz, \
                             tk_xxxyz_xzzz, \
                             tk_xxxyz_yyyy, \
                             tk_xxxyz_yyyz, \
                             tk_xxxyz_yyzz, \
                             tk_xxxyz_yzzz, \
                             tk_xxxyz_zzzz, \
                             tk_xxxz_xxxx,  \
                             tk_xxxz_xxxz,  \
                             tk_xxxz_xxyz,  \
                             tk_xxxz_xxz,   \
                             tk_xxxz_xxzz,  \
                             tk_xxxz_xyyz,  \
                             tk_xxxz_xyz,   \
                             tk_xxxz_xyzz,  \
                             tk_xxxz_xzz,   \
                             tk_xxxz_xzzz,  \
                             tk_xxxz_yyyz,  \
                             tk_xxxz_yyz,   \
                             tk_xxxz_yyzz,  \
                             tk_xxxz_yzz,   \
                             tk_xxxz_yzzz,  \
                             tk_xxxz_zzz,   \
                             tk_xxxz_zzzz,  \
                             ts_xxxyz_xxxx, \
                             ts_xxxyz_xxxy, \
                             ts_xxxyz_xxxz, \
                             ts_xxxyz_xxyy, \
                             ts_xxxyz_xxyz, \
                             ts_xxxyz_xxzz, \
                             ts_xxxyz_xyyy, \
                             ts_xxxyz_xyyz, \
                             ts_xxxyz_xyzz, \
                             ts_xxxyz_xzzz, \
                             ts_xxxyz_yyyy, \
                             ts_xxxyz_yyyz, \
                             ts_xxxyz_yyzz, \
                             ts_xxxyz_yzzz, \
                             ts_xxxyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyz_xxxx[i] = tk_xxxz_xxxx[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxx[i] * fz_0;

        tk_xxxyz_xxxy[i] = tk_xxxy_xxxy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxxy[i] * fz_0;

        tk_xxxyz_xxxz[i] = tk_xxxz_xxxz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxz[i] * fz_0;

        tk_xxxyz_xxyy[i] = tk_xxxy_xxyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxyy[i] * fz_0;

        tk_xxxyz_xxyz[i] = tk_xxxz_xxz[i] * fe_0 + tk_xxxz_xxyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxyz[i] * fz_0;

        tk_xxxyz_xxzz[i] = tk_xxxz_xxzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxzz[i] * fz_0;

        tk_xxxyz_xyyy[i] = tk_xxxy_xyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xyyy[i] * fz_0;

        tk_xxxyz_xyyz[i] = 2.0 * tk_xxxz_xyz[i] * fe_0 + tk_xxxz_xyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyyz[i] * fz_0;

        tk_xxxyz_xyzz[i] = tk_xxxz_xzz[i] * fe_0 + tk_xxxz_xyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyzz[i] * fz_0;

        tk_xxxyz_xzzz[i] = tk_xxxz_xzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xzzz[i] * fz_0;

        tk_xxxyz_yyyy[i] = tk_xxxy_yyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_yyyy[i] * fz_0;

        tk_xxxyz_yyyz[i] = 3.0 * tk_xxxz_yyz[i] * fe_0 + tk_xxxz_yyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyyz[i] * fz_0;

        tk_xxxyz_yyzz[i] = 2.0 * tk_xxxz_yzz[i] * fe_0 + tk_xxxz_yyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyzz[i] * fz_0;

        tk_xxxyz_yzzz[i] = tk_xxxz_zzz[i] * fe_0 + tk_xxxz_yzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yzzz[i] * fz_0;

        tk_xxxyz_zzzz[i] = tk_xxxz_zzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_zzzz[i] * fz_0;
    }

    // Set up 75-90 components of targeted buffer : HG

    auto tk_xxxzz_xxxx = pbuffer.data(idx_kin_hg + 75);

    auto tk_xxxzz_xxxy = pbuffer.data(idx_kin_hg + 76);

    auto tk_xxxzz_xxxz = pbuffer.data(idx_kin_hg + 77);

    auto tk_xxxzz_xxyy = pbuffer.data(idx_kin_hg + 78);

    auto tk_xxxzz_xxyz = pbuffer.data(idx_kin_hg + 79);

    auto tk_xxxzz_xxzz = pbuffer.data(idx_kin_hg + 80);

    auto tk_xxxzz_xyyy = pbuffer.data(idx_kin_hg + 81);

    auto tk_xxxzz_xyyz = pbuffer.data(idx_kin_hg + 82);

    auto tk_xxxzz_xyzz = pbuffer.data(idx_kin_hg + 83);

    auto tk_xxxzz_xzzz = pbuffer.data(idx_kin_hg + 84);

    auto tk_xxxzz_yyyy = pbuffer.data(idx_kin_hg + 85);

    auto tk_xxxzz_yyyz = pbuffer.data(idx_kin_hg + 86);

    auto tk_xxxzz_yyzz = pbuffer.data(idx_kin_hg + 87);

    auto tk_xxxzz_yzzz = pbuffer.data(idx_kin_hg + 88);

    auto tk_xxxzz_zzzz = pbuffer.data(idx_kin_hg + 89);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xxx_xxxx,   \
                             tk_xxx_xxxy,   \
                             tk_xxx_xxyy,   \
                             tk_xxx_xyyy,   \
                             tk_xxxz_xxxx,  \
                             tk_xxxz_xxxy,  \
                             tk_xxxz_xxyy,  \
                             tk_xxxz_xyyy,  \
                             tk_xxxzz_xxxx, \
                             tk_xxxzz_xxxy, \
                             tk_xxxzz_xxxz, \
                             tk_xxxzz_xxyy, \
                             tk_xxxzz_xxyz, \
                             tk_xxxzz_xxzz, \
                             tk_xxxzz_xyyy, \
                             tk_xxxzz_xyyz, \
                             tk_xxxzz_xyzz, \
                             tk_xxxzz_xzzz, \
                             tk_xxxzz_yyyy, \
                             tk_xxxzz_yyyz, \
                             tk_xxxzz_yyzz, \
                             tk_xxxzz_yzzz, \
                             tk_xxxzz_zzzz, \
                             tk_xxzz_xxxz,  \
                             tk_xxzz_xxyz,  \
                             tk_xxzz_xxz,   \
                             tk_xxzz_xxzz,  \
                             tk_xxzz_xyyz,  \
                             tk_xxzz_xyz,   \
                             tk_xxzz_xyzz,  \
                             tk_xxzz_xzz,   \
                             tk_xxzz_xzzz,  \
                             tk_xxzz_yyyy,  \
                             tk_xxzz_yyyz,  \
                             tk_xxzz_yyz,   \
                             tk_xxzz_yyzz,  \
                             tk_xxzz_yzz,   \
                             tk_xxzz_yzzz,  \
                             tk_xxzz_zzz,   \
                             tk_xxzz_zzzz,  \
                             tk_xzz_xxxz,   \
                             tk_xzz_xxyz,   \
                             tk_xzz_xxzz,   \
                             tk_xzz_xyyz,   \
                             tk_xzz_xyzz,   \
                             tk_xzz_xzzz,   \
                             tk_xzz_yyyy,   \
                             tk_xzz_yyyz,   \
                             tk_xzz_yyzz,   \
                             tk_xzz_yzzz,   \
                             tk_xzz_zzzz,   \
                             ts_xxx_xxxx,   \
                             ts_xxx_xxxy,   \
                             ts_xxx_xxyy,   \
                             ts_xxx_xyyy,   \
                             ts_xxxzz_xxxx, \
                             ts_xxxzz_xxxy, \
                             ts_xxxzz_xxxz, \
                             ts_xxxzz_xxyy, \
                             ts_xxxzz_xxyz, \
                             ts_xxxzz_xxzz, \
                             ts_xxxzz_xyyy, \
                             ts_xxxzz_xyyz, \
                             ts_xxxzz_xyzz, \
                             ts_xxxzz_xzzz, \
                             ts_xxxzz_yyyy, \
                             ts_xxxzz_yyyz, \
                             ts_xxxzz_yyzz, \
                             ts_xxxzz_yzzz, \
                             ts_xxxzz_zzzz, \
                             ts_xzz_xxxz,   \
                             ts_xzz_xxyz,   \
                             ts_xzz_xxzz,   \
                             ts_xzz_xyyz,   \
                             ts_xzz_xyzz,   \
                             ts_xzz_xzzz,   \
                             ts_xzz_yyyy,   \
                             ts_xzz_yyyz,   \
                             ts_xzz_yyzz,   \
                             ts_xzz_yzzz,   \
                             ts_xzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzz_xxxx[i] = -2.0 * ts_xxx_xxxx[i] * fbe_0 * fz_0 + tk_xxx_xxxx[i] * fe_0 + tk_xxxz_xxxx[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxx[i] * fz_0;

        tk_xxxzz_xxxy[i] = -2.0 * ts_xxx_xxxy[i] * fbe_0 * fz_0 + tk_xxx_xxxy[i] * fe_0 + tk_xxxz_xxxy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxy[i] * fz_0;

        tk_xxxzz_xxxz[i] = -4.0 * ts_xzz_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxz[i] * fe_0 + 3.0 * tk_xxzz_xxz[i] * fe_0 +
                           tk_xxzz_xxxz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxz[i] * fz_0;

        tk_xxxzz_xxyy[i] = -2.0 * ts_xxx_xxyy[i] * fbe_0 * fz_0 + tk_xxx_xxyy[i] * fe_0 + tk_xxxz_xxyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxyy[i] * fz_0;

        tk_xxxzz_xxyz[i] = -4.0 * ts_xzz_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxyz[i] * fe_0 + 2.0 * tk_xxzz_xyz[i] * fe_0 +
                           tk_xxzz_xxyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxyz[i] * fz_0;

        tk_xxxzz_xxzz[i] = -4.0 * ts_xzz_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxzz[i] * fe_0 + 2.0 * tk_xxzz_xzz[i] * fe_0 +
                           tk_xxzz_xxzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxzz[i] * fz_0;

        tk_xxxzz_xyyy[i] = -2.0 * ts_xxx_xyyy[i] * fbe_0 * fz_0 + tk_xxx_xyyy[i] * fe_0 + tk_xxxz_xyyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xyyy[i] * fz_0;

        tk_xxxzz_xyyz[i] = -4.0 * ts_xzz_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyyz[i] * fe_0 + tk_xxzz_yyz[i] * fe_0 + tk_xxzz_xyyz[i] * pa_x[i] +
                           2.0 * ts_xxxzz_xyyz[i] * fz_0;

        tk_xxxzz_xyzz[i] = -4.0 * ts_xzz_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyzz[i] * fe_0 + tk_xxzz_yzz[i] * fe_0 + tk_xxzz_xyzz[i] * pa_x[i] +
                           2.0 * ts_xxxzz_xyzz[i] * fz_0;

        tk_xxxzz_xzzz[i] = -4.0 * ts_xzz_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xzzz[i] * fe_0 + tk_xxzz_zzz[i] * fe_0 + tk_xxzz_xzzz[i] * pa_x[i] +
                           2.0 * ts_xxxzz_xzzz[i] * fz_0;

        tk_xxxzz_yyyy[i] =
            -4.0 * ts_xzz_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyy[i] * fe_0 + tk_xxzz_yyyy[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyy[i] * fz_0;

        tk_xxxzz_yyyz[i] =
            -4.0 * ts_xzz_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyz[i] * fe_0 + tk_xxzz_yyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyz[i] * fz_0;

        tk_xxxzz_yyzz[i] =
            -4.0 * ts_xzz_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyzz[i] * fe_0 + tk_xxzz_yyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyzz[i] * fz_0;

        tk_xxxzz_yzzz[i] =
            -4.0 * ts_xzz_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yzzz[i] * fe_0 + tk_xxzz_yzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yzzz[i] * fz_0;

        tk_xxxzz_zzzz[i] =
            -4.0 * ts_xzz_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_zzzz[i] * fe_0 + tk_xxzz_zzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_zzzz[i] * fz_0;
    }

    // Set up 90-105 components of targeted buffer : HG

    auto tk_xxyyy_xxxx = pbuffer.data(idx_kin_hg + 90);

    auto tk_xxyyy_xxxy = pbuffer.data(idx_kin_hg + 91);

    auto tk_xxyyy_xxxz = pbuffer.data(idx_kin_hg + 92);

    auto tk_xxyyy_xxyy = pbuffer.data(idx_kin_hg + 93);

    auto tk_xxyyy_xxyz = pbuffer.data(idx_kin_hg + 94);

    auto tk_xxyyy_xxzz = pbuffer.data(idx_kin_hg + 95);

    auto tk_xxyyy_xyyy = pbuffer.data(idx_kin_hg + 96);

    auto tk_xxyyy_xyyz = pbuffer.data(idx_kin_hg + 97);

    auto tk_xxyyy_xyzz = pbuffer.data(idx_kin_hg + 98);

    auto tk_xxyyy_xzzz = pbuffer.data(idx_kin_hg + 99);

    auto tk_xxyyy_yyyy = pbuffer.data(idx_kin_hg + 100);

    auto tk_xxyyy_yyyz = pbuffer.data(idx_kin_hg + 101);

    auto tk_xxyyy_yyzz = pbuffer.data(idx_kin_hg + 102);

    auto tk_xxyyy_yzzz = pbuffer.data(idx_kin_hg + 103);

    auto tk_xxyyy_zzzz = pbuffer.data(idx_kin_hg + 104);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xxy_xxxx,   \
                             tk_xxy_xxxz,   \
                             tk_xxy_xxzz,   \
                             tk_xxy_xzzz,   \
                             tk_xxyy_xxxx,  \
                             tk_xxyy_xxxz,  \
                             tk_xxyy_xxzz,  \
                             tk_xxyy_xzzz,  \
                             tk_xxyyy_xxxx, \
                             tk_xxyyy_xxxy, \
                             tk_xxyyy_xxxz, \
                             tk_xxyyy_xxyy, \
                             tk_xxyyy_xxyz, \
                             tk_xxyyy_xxzz, \
                             tk_xxyyy_xyyy, \
                             tk_xxyyy_xyyz, \
                             tk_xxyyy_xyzz, \
                             tk_xxyyy_xzzz, \
                             tk_xxyyy_yyyy, \
                             tk_xxyyy_yyyz, \
                             tk_xxyyy_yyzz, \
                             tk_xxyyy_yzzz, \
                             tk_xxyyy_zzzz, \
                             tk_xyyy_xxxy,  \
                             tk_xyyy_xxy,   \
                             tk_xyyy_xxyy,  \
                             tk_xyyy_xxyz,  \
                             tk_xyyy_xyy,   \
                             tk_xyyy_xyyy,  \
                             tk_xyyy_xyyz,  \
                             tk_xyyy_xyz,   \
                             tk_xyyy_xyzz,  \
                             tk_xyyy_yyy,   \
                             tk_xyyy_yyyy,  \
                             tk_xyyy_yyyz,  \
                             tk_xyyy_yyz,   \
                             tk_xyyy_yyzz,  \
                             tk_xyyy_yzz,   \
                             tk_xyyy_yzzz,  \
                             tk_xyyy_zzzz,  \
                             tk_yyy_xxxy,   \
                             tk_yyy_xxyy,   \
                             tk_yyy_xxyz,   \
                             tk_yyy_xyyy,   \
                             tk_yyy_xyyz,   \
                             tk_yyy_xyzz,   \
                             tk_yyy_yyyy,   \
                             tk_yyy_yyyz,   \
                             tk_yyy_yyzz,   \
                             tk_yyy_yzzz,   \
                             tk_yyy_zzzz,   \
                             ts_xxy_xxxx,   \
                             ts_xxy_xxxz,   \
                             ts_xxy_xxzz,   \
                             ts_xxy_xzzz,   \
                             ts_xxyyy_xxxx, \
                             ts_xxyyy_xxxy, \
                             ts_xxyyy_xxxz, \
                             ts_xxyyy_xxyy, \
                             ts_xxyyy_xxyz, \
                             ts_xxyyy_xxzz, \
                             ts_xxyyy_xyyy, \
                             ts_xxyyy_xyyz, \
                             ts_xxyyy_xyzz, \
                             ts_xxyyy_xzzz, \
                             ts_xxyyy_yyyy, \
                             ts_xxyyy_yyyz, \
                             ts_xxyyy_yyzz, \
                             ts_xxyyy_yzzz, \
                             ts_xxyyy_zzzz, \
                             ts_yyy_xxxy,   \
                             ts_yyy_xxyy,   \
                             ts_yyy_xxyz,   \
                             ts_yyy_xyyy,   \
                             ts_yyy_xyyz,   \
                             ts_yyy_xyzz,   \
                             ts_yyy_yyyy,   \
                             ts_yyy_yyyz,   \
                             ts_yyy_yyzz,   \
                             ts_yyy_yzzz,   \
                             ts_yyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyy_xxxx[i] =
            -4.0 * ts_xxy_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxx[i] * fe_0 + tk_xxyy_xxxx[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxx[i] * fz_0;

        tk_xxyyy_xxxy[i] = -2.0 * ts_yyy_xxxy[i] * fbe_0 * fz_0 + tk_yyy_xxxy[i] * fe_0 + 3.0 * tk_xyyy_xxy[i] * fe_0 + tk_xyyy_xxxy[i] * pa_x[i] +
                           2.0 * ts_xxyyy_xxxy[i] * fz_0;

        tk_xxyyy_xxxz[i] =
            -4.0 * ts_xxy_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxz[i] * fe_0 + tk_xxyy_xxxz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxz[i] * fz_0;

        tk_xxyyy_xxyy[i] = -2.0 * ts_yyy_xxyy[i] * fbe_0 * fz_0 + tk_yyy_xxyy[i] * fe_0 + 2.0 * tk_xyyy_xyy[i] * fe_0 + tk_xyyy_xxyy[i] * pa_x[i] +
                           2.0 * ts_xxyyy_xxyy[i] * fz_0;

        tk_xxyyy_xxyz[i] = -2.0 * ts_yyy_xxyz[i] * fbe_0 * fz_0 + tk_yyy_xxyz[i] * fe_0 + 2.0 * tk_xyyy_xyz[i] * fe_0 + tk_xyyy_xxyz[i] * pa_x[i] +
                           2.0 * ts_xxyyy_xxyz[i] * fz_0;

        tk_xxyyy_xxzz[i] =
            -4.0 * ts_xxy_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxzz[i] * fe_0 + tk_xxyy_xxzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxzz[i] * fz_0;

        tk_xxyyy_xyyy[i] = -2.0 * ts_yyy_xyyy[i] * fbe_0 * fz_0 + tk_yyy_xyyy[i] * fe_0 + tk_xyyy_yyy[i] * fe_0 + tk_xyyy_xyyy[i] * pa_x[i] +
                           2.0 * ts_xxyyy_xyyy[i] * fz_0;

        tk_xxyyy_xyyz[i] = -2.0 * ts_yyy_xyyz[i] * fbe_0 * fz_0 + tk_yyy_xyyz[i] * fe_0 + tk_xyyy_yyz[i] * fe_0 + tk_xyyy_xyyz[i] * pa_x[i] +
                           2.0 * ts_xxyyy_xyyz[i] * fz_0;

        tk_xxyyy_xyzz[i] = -2.0 * ts_yyy_xyzz[i] * fbe_0 * fz_0 + tk_yyy_xyzz[i] * fe_0 + tk_xyyy_yzz[i] * fe_0 + tk_xyyy_xyzz[i] * pa_x[i] +
                           2.0 * ts_xxyyy_xyzz[i] * fz_0;

        tk_xxyyy_xzzz[i] =
            -4.0 * ts_xxy_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xzzz[i] * fe_0 + tk_xxyy_xzzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xzzz[i] * fz_0;

        tk_xxyyy_yyyy[i] = -2.0 * ts_yyy_yyyy[i] * fbe_0 * fz_0 + tk_yyy_yyyy[i] * fe_0 + tk_xyyy_yyyy[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyy[i] * fz_0;

        tk_xxyyy_yyyz[i] = -2.0 * ts_yyy_yyyz[i] * fbe_0 * fz_0 + tk_yyy_yyyz[i] * fe_0 + tk_xyyy_yyyz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyz[i] * fz_0;

        tk_xxyyy_yyzz[i] = -2.0 * ts_yyy_yyzz[i] * fbe_0 * fz_0 + tk_yyy_yyzz[i] * fe_0 + tk_xyyy_yyzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyzz[i] * fz_0;

        tk_xxyyy_yzzz[i] = -2.0 * ts_yyy_yzzz[i] * fbe_0 * fz_0 + tk_yyy_yzzz[i] * fe_0 + tk_xyyy_yzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yzzz[i] * fz_0;

        tk_xxyyy_zzzz[i] = -2.0 * ts_yyy_zzzz[i] * fbe_0 * fz_0 + tk_yyy_zzzz[i] * fe_0 + tk_xyyy_zzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_zzzz[i] * fz_0;
    }

    // Set up 105-120 components of targeted buffer : HG

    auto tk_xxyyz_xxxx = pbuffer.data(idx_kin_hg + 105);

    auto tk_xxyyz_xxxy = pbuffer.data(idx_kin_hg + 106);

    auto tk_xxyyz_xxxz = pbuffer.data(idx_kin_hg + 107);

    auto tk_xxyyz_xxyy = pbuffer.data(idx_kin_hg + 108);

    auto tk_xxyyz_xxyz = pbuffer.data(idx_kin_hg + 109);

    auto tk_xxyyz_xxzz = pbuffer.data(idx_kin_hg + 110);

    auto tk_xxyyz_xyyy = pbuffer.data(idx_kin_hg + 111);

    auto tk_xxyyz_xyyz = pbuffer.data(idx_kin_hg + 112);

    auto tk_xxyyz_xyzz = pbuffer.data(idx_kin_hg + 113);

    auto tk_xxyyz_xzzz = pbuffer.data(idx_kin_hg + 114);

    auto tk_xxyyz_yyyy = pbuffer.data(idx_kin_hg + 115);

    auto tk_xxyyz_yyyz = pbuffer.data(idx_kin_hg + 116);

    auto tk_xxyyz_yyzz = pbuffer.data(idx_kin_hg + 117);

    auto tk_xxyyz_yzzz = pbuffer.data(idx_kin_hg + 118);

    auto tk_xxyyz_zzzz = pbuffer.data(idx_kin_hg + 119);

#pragma omp simd aligned(pa_z,              \
                             tk_xxyy_xxx,   \
                             tk_xxyy_xxxx,  \
                             tk_xxyy_xxxy,  \
                             tk_xxyy_xxxz,  \
                             tk_xxyy_xxy,   \
                             tk_xxyy_xxyy,  \
                             tk_xxyy_xxyz,  \
                             tk_xxyy_xxz,   \
                             tk_xxyy_xxzz,  \
                             tk_xxyy_xyy,   \
                             tk_xxyy_xyyy,  \
                             tk_xxyy_xyyz,  \
                             tk_xxyy_xyz,   \
                             tk_xxyy_xyzz,  \
                             tk_xxyy_xzz,   \
                             tk_xxyy_xzzz,  \
                             tk_xxyy_yyy,   \
                             tk_xxyy_yyyy,  \
                             tk_xxyy_yyyz,  \
                             tk_xxyy_yyz,   \
                             tk_xxyy_yyzz,  \
                             tk_xxyy_yzz,   \
                             tk_xxyy_yzzz,  \
                             tk_xxyy_zzz,   \
                             tk_xxyy_zzzz,  \
                             tk_xxyyz_xxxx, \
                             tk_xxyyz_xxxy, \
                             tk_xxyyz_xxxz, \
                             tk_xxyyz_xxyy, \
                             tk_xxyyz_xxyz, \
                             tk_xxyyz_xxzz, \
                             tk_xxyyz_xyyy, \
                             tk_xxyyz_xyyz, \
                             tk_xxyyz_xyzz, \
                             tk_xxyyz_xzzz, \
                             tk_xxyyz_yyyy, \
                             tk_xxyyz_yyyz, \
                             tk_xxyyz_yyzz, \
                             tk_xxyyz_yzzz, \
                             tk_xxyyz_zzzz, \
                             ts_xxyyz_xxxx, \
                             ts_xxyyz_xxxy, \
                             ts_xxyyz_xxxz, \
                             ts_xxyyz_xxyy, \
                             ts_xxyyz_xxyz, \
                             ts_xxyyz_xxzz, \
                             ts_xxyyz_xyyy, \
                             ts_xxyyz_xyyz, \
                             ts_xxyyz_xyzz, \
                             ts_xxyyz_xzzz, \
                             ts_xxyyz_yyyy, \
                             ts_xxyyz_yyyz, \
                             ts_xxyyz_yyzz, \
                             ts_xxyyz_yzzz, \
                             ts_xxyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyz_xxxx[i] = tk_xxyy_xxxx[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxx[i] * fz_0;

        tk_xxyyz_xxxy[i] = tk_xxyy_xxxy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxy[i] * fz_0;

        tk_xxyyz_xxxz[i] = tk_xxyy_xxx[i] * fe_0 + tk_xxyy_xxxz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxz[i] * fz_0;

        tk_xxyyz_xxyy[i] = tk_xxyy_xxyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyy[i] * fz_0;

        tk_xxyyz_xxyz[i] = tk_xxyy_xxy[i] * fe_0 + tk_xxyy_xxyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyz[i] * fz_0;

        tk_xxyyz_xxzz[i] = 2.0 * tk_xxyy_xxz[i] * fe_0 + tk_xxyy_xxzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxzz[i] * fz_0;

        tk_xxyyz_xyyy[i] = tk_xxyy_xyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyy[i] * fz_0;

        tk_xxyyz_xyyz[i] = tk_xxyy_xyy[i] * fe_0 + tk_xxyy_xyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyz[i] * fz_0;

        tk_xxyyz_xyzz[i] = 2.0 * tk_xxyy_xyz[i] * fe_0 + tk_xxyy_xyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyzz[i] * fz_0;

        tk_xxyyz_xzzz[i] = 3.0 * tk_xxyy_xzz[i] * fe_0 + tk_xxyy_xzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xzzz[i] * fz_0;

        tk_xxyyz_yyyy[i] = tk_xxyy_yyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyy[i] * fz_0;

        tk_xxyyz_yyyz[i] = tk_xxyy_yyy[i] * fe_0 + tk_xxyy_yyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyz[i] * fz_0;

        tk_xxyyz_yyzz[i] = 2.0 * tk_xxyy_yyz[i] * fe_0 + tk_xxyy_yyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyzz[i] * fz_0;

        tk_xxyyz_yzzz[i] = 3.0 * tk_xxyy_yzz[i] * fe_0 + tk_xxyy_yzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yzzz[i] * fz_0;

        tk_xxyyz_zzzz[i] = 4.0 * tk_xxyy_zzz[i] * fe_0 + tk_xxyy_zzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_zzzz[i] * fz_0;
    }

    // Set up 120-135 components of targeted buffer : HG

    auto tk_xxyzz_xxxx = pbuffer.data(idx_kin_hg + 120);

    auto tk_xxyzz_xxxy = pbuffer.data(idx_kin_hg + 121);

    auto tk_xxyzz_xxxz = pbuffer.data(idx_kin_hg + 122);

    auto tk_xxyzz_xxyy = pbuffer.data(idx_kin_hg + 123);

    auto tk_xxyzz_xxyz = pbuffer.data(idx_kin_hg + 124);

    auto tk_xxyzz_xxzz = pbuffer.data(idx_kin_hg + 125);

    auto tk_xxyzz_xyyy = pbuffer.data(idx_kin_hg + 126);

    auto tk_xxyzz_xyyz = pbuffer.data(idx_kin_hg + 127);

    auto tk_xxyzz_xyzz = pbuffer.data(idx_kin_hg + 128);

    auto tk_xxyzz_xzzz = pbuffer.data(idx_kin_hg + 129);

    auto tk_xxyzz_yyyy = pbuffer.data(idx_kin_hg + 130);

    auto tk_xxyzz_yyyz = pbuffer.data(idx_kin_hg + 131);

    auto tk_xxyzz_yyzz = pbuffer.data(idx_kin_hg + 132);

    auto tk_xxyzz_yzzz = pbuffer.data(idx_kin_hg + 133);

    auto tk_xxyzz_zzzz = pbuffer.data(idx_kin_hg + 134);

#pragma omp simd aligned(pa_y,              \
                             tk_xxyzz_xxxx, \
                             tk_xxyzz_xxxy, \
                             tk_xxyzz_xxxz, \
                             tk_xxyzz_xxyy, \
                             tk_xxyzz_xxyz, \
                             tk_xxyzz_xxzz, \
                             tk_xxyzz_xyyy, \
                             tk_xxyzz_xyyz, \
                             tk_xxyzz_xyzz, \
                             tk_xxyzz_xzzz, \
                             tk_xxyzz_yyyy, \
                             tk_xxyzz_yyyz, \
                             tk_xxyzz_yyzz, \
                             tk_xxyzz_yzzz, \
                             tk_xxyzz_zzzz, \
                             tk_xxzz_xxx,   \
                             tk_xxzz_xxxx,  \
                             tk_xxzz_xxxy,  \
                             tk_xxzz_xxxz,  \
                             tk_xxzz_xxy,   \
                             tk_xxzz_xxyy,  \
                             tk_xxzz_xxyz,  \
                             tk_xxzz_xxz,   \
                             tk_xxzz_xxzz,  \
                             tk_xxzz_xyy,   \
                             tk_xxzz_xyyy,  \
                             tk_xxzz_xyyz,  \
                             tk_xxzz_xyz,   \
                             tk_xxzz_xyzz,  \
                             tk_xxzz_xzz,   \
                             tk_xxzz_xzzz,  \
                             tk_xxzz_yyy,   \
                             tk_xxzz_yyyy,  \
                             tk_xxzz_yyyz,  \
                             tk_xxzz_yyz,   \
                             tk_xxzz_yyzz,  \
                             tk_xxzz_yzz,   \
                             tk_xxzz_yzzz,  \
                             tk_xxzz_zzz,   \
                             tk_xxzz_zzzz,  \
                             ts_xxyzz_xxxx, \
                             ts_xxyzz_xxxy, \
                             ts_xxyzz_xxxz, \
                             ts_xxyzz_xxyy, \
                             ts_xxyzz_xxyz, \
                             ts_xxyzz_xxzz, \
                             ts_xxyzz_xyyy, \
                             ts_xxyzz_xyyz, \
                             ts_xxyzz_xyzz, \
                             ts_xxyzz_xzzz, \
                             ts_xxyzz_yyyy, \
                             ts_xxyzz_yyyz, \
                             ts_xxyzz_yyzz, \
                             ts_xxyzz_yzzz, \
                             ts_xxyzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzz_xxxx[i] = tk_xxzz_xxxx[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxx[i] * fz_0;

        tk_xxyzz_xxxy[i] = tk_xxzz_xxx[i] * fe_0 + tk_xxzz_xxxy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxy[i] * fz_0;

        tk_xxyzz_xxxz[i] = tk_xxzz_xxxz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxz[i] * fz_0;

        tk_xxyzz_xxyy[i] = 2.0 * tk_xxzz_xxy[i] * fe_0 + tk_xxzz_xxyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyy[i] * fz_0;

        tk_xxyzz_xxyz[i] = tk_xxzz_xxz[i] * fe_0 + tk_xxzz_xxyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyz[i] * fz_0;

        tk_xxyzz_xxzz[i] = tk_xxzz_xxzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxzz[i] * fz_0;

        tk_xxyzz_xyyy[i] = 3.0 * tk_xxzz_xyy[i] * fe_0 + tk_xxzz_xyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyy[i] * fz_0;

        tk_xxyzz_xyyz[i] = 2.0 * tk_xxzz_xyz[i] * fe_0 + tk_xxzz_xyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyz[i] * fz_0;

        tk_xxyzz_xyzz[i] = tk_xxzz_xzz[i] * fe_0 + tk_xxzz_xyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyzz[i] * fz_0;

        tk_xxyzz_xzzz[i] = tk_xxzz_xzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xzzz[i] * fz_0;

        tk_xxyzz_yyyy[i] = 4.0 * tk_xxzz_yyy[i] * fe_0 + tk_xxzz_yyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyy[i] * fz_0;

        tk_xxyzz_yyyz[i] = 3.0 * tk_xxzz_yyz[i] * fe_0 + tk_xxzz_yyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyz[i] * fz_0;

        tk_xxyzz_yyzz[i] = 2.0 * tk_xxzz_yzz[i] * fe_0 + tk_xxzz_yyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyzz[i] * fz_0;

        tk_xxyzz_yzzz[i] = tk_xxzz_zzz[i] * fe_0 + tk_xxzz_yzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yzzz[i] * fz_0;

        tk_xxyzz_zzzz[i] = tk_xxzz_zzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_zzzz[i] * fz_0;
    }

    // Set up 135-150 components of targeted buffer : HG

    auto tk_xxzzz_xxxx = pbuffer.data(idx_kin_hg + 135);

    auto tk_xxzzz_xxxy = pbuffer.data(idx_kin_hg + 136);

    auto tk_xxzzz_xxxz = pbuffer.data(idx_kin_hg + 137);

    auto tk_xxzzz_xxyy = pbuffer.data(idx_kin_hg + 138);

    auto tk_xxzzz_xxyz = pbuffer.data(idx_kin_hg + 139);

    auto tk_xxzzz_xxzz = pbuffer.data(idx_kin_hg + 140);

    auto tk_xxzzz_xyyy = pbuffer.data(idx_kin_hg + 141);

    auto tk_xxzzz_xyyz = pbuffer.data(idx_kin_hg + 142);

    auto tk_xxzzz_xyzz = pbuffer.data(idx_kin_hg + 143);

    auto tk_xxzzz_xzzz = pbuffer.data(idx_kin_hg + 144);

    auto tk_xxzzz_yyyy = pbuffer.data(idx_kin_hg + 145);

    auto tk_xxzzz_yyyz = pbuffer.data(idx_kin_hg + 146);

    auto tk_xxzzz_yyzz = pbuffer.data(idx_kin_hg + 147);

    auto tk_xxzzz_yzzz = pbuffer.data(idx_kin_hg + 148);

    auto tk_xxzzz_zzzz = pbuffer.data(idx_kin_hg + 149);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xxz_xxxx,   \
                             tk_xxz_xxxy,   \
                             tk_xxz_xxyy,   \
                             tk_xxz_xyyy,   \
                             tk_xxzz_xxxx,  \
                             tk_xxzz_xxxy,  \
                             tk_xxzz_xxyy,  \
                             tk_xxzz_xyyy,  \
                             tk_xxzzz_xxxx, \
                             tk_xxzzz_xxxy, \
                             tk_xxzzz_xxxz, \
                             tk_xxzzz_xxyy, \
                             tk_xxzzz_xxyz, \
                             tk_xxzzz_xxzz, \
                             tk_xxzzz_xyyy, \
                             tk_xxzzz_xyyz, \
                             tk_xxzzz_xyzz, \
                             tk_xxzzz_xzzz, \
                             tk_xxzzz_yyyy, \
                             tk_xxzzz_yyyz, \
                             tk_xxzzz_yyzz, \
                             tk_xxzzz_yzzz, \
                             tk_xxzzz_zzzz, \
                             tk_xzzz_xxxz,  \
                             tk_xzzz_xxyz,  \
                             tk_xzzz_xxz,   \
                             tk_xzzz_xxzz,  \
                             tk_xzzz_xyyz,  \
                             tk_xzzz_xyz,   \
                             tk_xzzz_xyzz,  \
                             tk_xzzz_xzz,   \
                             tk_xzzz_xzzz,  \
                             tk_xzzz_yyyy,  \
                             tk_xzzz_yyyz,  \
                             tk_xzzz_yyz,   \
                             tk_xzzz_yyzz,  \
                             tk_xzzz_yzz,   \
                             tk_xzzz_yzzz,  \
                             tk_xzzz_zzz,   \
                             tk_xzzz_zzzz,  \
                             tk_zzz_xxxz,   \
                             tk_zzz_xxyz,   \
                             tk_zzz_xxzz,   \
                             tk_zzz_xyyz,   \
                             tk_zzz_xyzz,   \
                             tk_zzz_xzzz,   \
                             tk_zzz_yyyy,   \
                             tk_zzz_yyyz,   \
                             tk_zzz_yyzz,   \
                             tk_zzz_yzzz,   \
                             tk_zzz_zzzz,   \
                             ts_xxz_xxxx,   \
                             ts_xxz_xxxy,   \
                             ts_xxz_xxyy,   \
                             ts_xxz_xyyy,   \
                             ts_xxzzz_xxxx, \
                             ts_xxzzz_xxxy, \
                             ts_xxzzz_xxxz, \
                             ts_xxzzz_xxyy, \
                             ts_xxzzz_xxyz, \
                             ts_xxzzz_xxzz, \
                             ts_xxzzz_xyyy, \
                             ts_xxzzz_xyyz, \
                             ts_xxzzz_xyzz, \
                             ts_xxzzz_xzzz, \
                             ts_xxzzz_yyyy, \
                             ts_xxzzz_yyyz, \
                             ts_xxzzz_yyzz, \
                             ts_xxzzz_yzzz, \
                             ts_xxzzz_zzzz, \
                             ts_zzz_xxxz,   \
                             ts_zzz_xxyz,   \
                             ts_zzz_xxzz,   \
                             ts_zzz_xyyz,   \
                             ts_zzz_xyzz,   \
                             ts_zzz_xzzz,   \
                             ts_zzz_yyyy,   \
                             ts_zzz_yyyz,   \
                             ts_zzz_yyzz,   \
                             ts_zzz_yzzz,   \
                             ts_zzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzz_xxxx[i] =
            -4.0 * ts_xxz_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxx[i] * fe_0 + tk_xxzz_xxxx[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxx[i] * fz_0;

        tk_xxzzz_xxxy[i] =
            -4.0 * ts_xxz_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxy[i] * fe_0 + tk_xxzz_xxxy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxy[i] * fz_0;

        tk_xxzzz_xxxz[i] = -2.0 * ts_zzz_xxxz[i] * fbe_0 * fz_0 + tk_zzz_xxxz[i] * fe_0 + 3.0 * tk_xzzz_xxz[i] * fe_0 + tk_xzzz_xxxz[i] * pa_x[i] +
                           2.0 * ts_xxzzz_xxxz[i] * fz_0;

        tk_xxzzz_xxyy[i] =
            -4.0 * ts_xxz_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxyy[i] * fe_0 + tk_xxzz_xxyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxyy[i] * fz_0;

        tk_xxzzz_xxyz[i] = -2.0 * ts_zzz_xxyz[i] * fbe_0 * fz_0 + tk_zzz_xxyz[i] * fe_0 + 2.0 * tk_xzzz_xyz[i] * fe_0 + tk_xzzz_xxyz[i] * pa_x[i] +
                           2.0 * ts_xxzzz_xxyz[i] * fz_0;

        tk_xxzzz_xxzz[i] = -2.0 * ts_zzz_xxzz[i] * fbe_0 * fz_0 + tk_zzz_xxzz[i] * fe_0 + 2.0 * tk_xzzz_xzz[i] * fe_0 + tk_xzzz_xxzz[i] * pa_x[i] +
                           2.0 * ts_xxzzz_xxzz[i] * fz_0;

        tk_xxzzz_xyyy[i] =
            -4.0 * ts_xxz_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xyyy[i] * fe_0 + tk_xxzz_xyyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xyyy[i] * fz_0;

        tk_xxzzz_xyyz[i] = -2.0 * ts_zzz_xyyz[i] * fbe_0 * fz_0 + tk_zzz_xyyz[i] * fe_0 + tk_xzzz_yyz[i] * fe_0 + tk_xzzz_xyyz[i] * pa_x[i] +
                           2.0 * ts_xxzzz_xyyz[i] * fz_0;

        tk_xxzzz_xyzz[i] = -2.0 * ts_zzz_xyzz[i] * fbe_0 * fz_0 + tk_zzz_xyzz[i] * fe_0 + tk_xzzz_yzz[i] * fe_0 + tk_xzzz_xyzz[i] * pa_x[i] +
                           2.0 * ts_xxzzz_xyzz[i] * fz_0;

        tk_xxzzz_xzzz[i] = -2.0 * ts_zzz_xzzz[i] * fbe_0 * fz_0 + tk_zzz_xzzz[i] * fe_0 + tk_xzzz_zzz[i] * fe_0 + tk_xzzz_xzzz[i] * pa_x[i] +
                           2.0 * ts_xxzzz_xzzz[i] * fz_0;

        tk_xxzzz_yyyy[i] = -2.0 * ts_zzz_yyyy[i] * fbe_0 * fz_0 + tk_zzz_yyyy[i] * fe_0 + tk_xzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyy[i] * fz_0;

        tk_xxzzz_yyyz[i] = -2.0 * ts_zzz_yyyz[i] * fbe_0 * fz_0 + tk_zzz_yyyz[i] * fe_0 + tk_xzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyz[i] * fz_0;

        tk_xxzzz_yyzz[i] = -2.0 * ts_zzz_yyzz[i] * fbe_0 * fz_0 + tk_zzz_yyzz[i] * fe_0 + tk_xzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyzz[i] * fz_0;

        tk_xxzzz_yzzz[i] = -2.0 * ts_zzz_yzzz[i] * fbe_0 * fz_0 + tk_zzz_yzzz[i] * fe_0 + tk_xzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yzzz[i] * fz_0;

        tk_xxzzz_zzzz[i] = -2.0 * ts_zzz_zzzz[i] * fbe_0 * fz_0 + tk_zzz_zzzz[i] * fe_0 + tk_xzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_zzzz[i] * fz_0;
    }

    // Set up 150-165 components of targeted buffer : HG

    auto tk_xyyyy_xxxx = pbuffer.data(idx_kin_hg + 150);

    auto tk_xyyyy_xxxy = pbuffer.data(idx_kin_hg + 151);

    auto tk_xyyyy_xxxz = pbuffer.data(idx_kin_hg + 152);

    auto tk_xyyyy_xxyy = pbuffer.data(idx_kin_hg + 153);

    auto tk_xyyyy_xxyz = pbuffer.data(idx_kin_hg + 154);

    auto tk_xyyyy_xxzz = pbuffer.data(idx_kin_hg + 155);

    auto tk_xyyyy_xyyy = pbuffer.data(idx_kin_hg + 156);

    auto tk_xyyyy_xyyz = pbuffer.data(idx_kin_hg + 157);

    auto tk_xyyyy_xyzz = pbuffer.data(idx_kin_hg + 158);

    auto tk_xyyyy_xzzz = pbuffer.data(idx_kin_hg + 159);

    auto tk_xyyyy_yyyy = pbuffer.data(idx_kin_hg + 160);

    auto tk_xyyyy_yyyz = pbuffer.data(idx_kin_hg + 161);

    auto tk_xyyyy_yyzz = pbuffer.data(idx_kin_hg + 162);

    auto tk_xyyyy_yzzz = pbuffer.data(idx_kin_hg + 163);

    auto tk_xyyyy_zzzz = pbuffer.data(idx_kin_hg + 164);

#pragma omp simd aligned(pa_x,              \
                             tk_xyyyy_xxxx, \
                             tk_xyyyy_xxxy, \
                             tk_xyyyy_xxxz, \
                             tk_xyyyy_xxyy, \
                             tk_xyyyy_xxyz, \
                             tk_xyyyy_xxzz, \
                             tk_xyyyy_xyyy, \
                             tk_xyyyy_xyyz, \
                             tk_xyyyy_xyzz, \
                             tk_xyyyy_xzzz, \
                             tk_xyyyy_yyyy, \
                             tk_xyyyy_yyyz, \
                             tk_xyyyy_yyzz, \
                             tk_xyyyy_yzzz, \
                             tk_xyyyy_zzzz, \
                             tk_yyyy_xxx,   \
                             tk_yyyy_xxxx,  \
                             tk_yyyy_xxxy,  \
                             tk_yyyy_xxxz,  \
                             tk_yyyy_xxy,   \
                             tk_yyyy_xxyy,  \
                             tk_yyyy_xxyz,  \
                             tk_yyyy_xxz,   \
                             tk_yyyy_xxzz,  \
                             tk_yyyy_xyy,   \
                             tk_yyyy_xyyy,  \
                             tk_yyyy_xyyz,  \
                             tk_yyyy_xyz,   \
                             tk_yyyy_xyzz,  \
                             tk_yyyy_xzz,   \
                             tk_yyyy_xzzz,  \
                             tk_yyyy_yyy,   \
                             tk_yyyy_yyyy,  \
                             tk_yyyy_yyyz,  \
                             tk_yyyy_yyz,   \
                             tk_yyyy_yyzz,  \
                             tk_yyyy_yzz,   \
                             tk_yyyy_yzzz,  \
                             tk_yyyy_zzz,   \
                             tk_yyyy_zzzz,  \
                             ts_xyyyy_xxxx, \
                             ts_xyyyy_xxxy, \
                             ts_xyyyy_xxxz, \
                             ts_xyyyy_xxyy, \
                             ts_xyyyy_xxyz, \
                             ts_xyyyy_xxzz, \
                             ts_xyyyy_xyyy, \
                             ts_xyyyy_xyyz, \
                             ts_xyyyy_xyzz, \
                             ts_xyyyy_xzzz, \
                             ts_xyyyy_yyyy, \
                             ts_xyyyy_yyyz, \
                             ts_xyyyy_yyzz, \
                             ts_xyyyy_yzzz, \
                             ts_xyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyy_xxxx[i] = 4.0 * tk_yyyy_xxx[i] * fe_0 + tk_yyyy_xxxx[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxx[i] * fz_0;

        tk_xyyyy_xxxy[i] = 3.0 * tk_yyyy_xxy[i] * fe_0 + tk_yyyy_xxxy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxy[i] * fz_0;

        tk_xyyyy_xxxz[i] = 3.0 * tk_yyyy_xxz[i] * fe_0 + tk_yyyy_xxxz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxz[i] * fz_0;

        tk_xyyyy_xxyy[i] = 2.0 * tk_yyyy_xyy[i] * fe_0 + tk_yyyy_xxyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyy[i] * fz_0;

        tk_xyyyy_xxyz[i] = 2.0 * tk_yyyy_xyz[i] * fe_0 + tk_yyyy_xxyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyz[i] * fz_0;

        tk_xyyyy_xxzz[i] = 2.0 * tk_yyyy_xzz[i] * fe_0 + tk_yyyy_xxzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxzz[i] * fz_0;

        tk_xyyyy_xyyy[i] = tk_yyyy_yyy[i] * fe_0 + tk_yyyy_xyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyy[i] * fz_0;

        tk_xyyyy_xyyz[i] = tk_yyyy_yyz[i] * fe_0 + tk_yyyy_xyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyz[i] * fz_0;

        tk_xyyyy_xyzz[i] = tk_yyyy_yzz[i] * fe_0 + tk_yyyy_xyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyzz[i] * fz_0;

        tk_xyyyy_xzzz[i] = tk_yyyy_zzz[i] * fe_0 + tk_yyyy_xzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xzzz[i] * fz_0;

        tk_xyyyy_yyyy[i] = tk_yyyy_yyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyy[i] * fz_0;

        tk_xyyyy_yyyz[i] = tk_yyyy_yyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyz[i] * fz_0;

        tk_xyyyy_yyzz[i] = tk_yyyy_yyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyzz[i] * fz_0;

        tk_xyyyy_yzzz[i] = tk_yyyy_yzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yzzz[i] * fz_0;

        tk_xyyyy_zzzz[i] = tk_yyyy_zzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_zzzz[i] * fz_0;
    }

    // Set up 165-180 components of targeted buffer : HG

    auto tk_xyyyz_xxxx = pbuffer.data(idx_kin_hg + 165);

    auto tk_xyyyz_xxxy = pbuffer.data(idx_kin_hg + 166);

    auto tk_xyyyz_xxxz = pbuffer.data(idx_kin_hg + 167);

    auto tk_xyyyz_xxyy = pbuffer.data(idx_kin_hg + 168);

    auto tk_xyyyz_xxyz = pbuffer.data(idx_kin_hg + 169);

    auto tk_xyyyz_xxzz = pbuffer.data(idx_kin_hg + 170);

    auto tk_xyyyz_xyyy = pbuffer.data(idx_kin_hg + 171);

    auto tk_xyyyz_xyyz = pbuffer.data(idx_kin_hg + 172);

    auto tk_xyyyz_xyzz = pbuffer.data(idx_kin_hg + 173);

    auto tk_xyyyz_xzzz = pbuffer.data(idx_kin_hg + 174);

    auto tk_xyyyz_yyyy = pbuffer.data(idx_kin_hg + 175);

    auto tk_xyyyz_yyyz = pbuffer.data(idx_kin_hg + 176);

    auto tk_xyyyz_yyzz = pbuffer.data(idx_kin_hg + 177);

    auto tk_xyyyz_yzzz = pbuffer.data(idx_kin_hg + 178);

    auto tk_xyyyz_zzzz = pbuffer.data(idx_kin_hg + 179);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xyyy_xxxx,  \
                             tk_xyyy_xxxy,  \
                             tk_xyyy_xxyy,  \
                             tk_xyyy_xyyy,  \
                             tk_xyyyz_xxxx, \
                             tk_xyyyz_xxxy, \
                             tk_xyyyz_xxxz, \
                             tk_xyyyz_xxyy, \
                             tk_xyyyz_xxyz, \
                             tk_xyyyz_xxzz, \
                             tk_xyyyz_xyyy, \
                             tk_xyyyz_xyyz, \
                             tk_xyyyz_xyzz, \
                             tk_xyyyz_xzzz, \
                             tk_xyyyz_yyyy, \
                             tk_xyyyz_yyyz, \
                             tk_xyyyz_yyzz, \
                             tk_xyyyz_yzzz, \
                             tk_xyyyz_zzzz, \
                             tk_yyyz_xxxz,  \
                             tk_yyyz_xxyz,  \
                             tk_yyyz_xxz,   \
                             tk_yyyz_xxzz,  \
                             tk_yyyz_xyyz,  \
                             tk_yyyz_xyz,   \
                             tk_yyyz_xyzz,  \
                             tk_yyyz_xzz,   \
                             tk_yyyz_xzzz,  \
                             tk_yyyz_yyyy,  \
                             tk_yyyz_yyyz,  \
                             tk_yyyz_yyz,   \
                             tk_yyyz_yyzz,  \
                             tk_yyyz_yzz,   \
                             tk_yyyz_yzzz,  \
                             tk_yyyz_zzz,   \
                             tk_yyyz_zzzz,  \
                             ts_xyyyz_xxxx, \
                             ts_xyyyz_xxxy, \
                             ts_xyyyz_xxxz, \
                             ts_xyyyz_xxyy, \
                             ts_xyyyz_xxyz, \
                             ts_xyyyz_xxzz, \
                             ts_xyyyz_xyyy, \
                             ts_xyyyz_xyyz, \
                             ts_xyyyz_xyzz, \
                             ts_xyyyz_xzzz, \
                             ts_xyyyz_yyyy, \
                             ts_xyyyz_yyyz, \
                             ts_xyyyz_yyzz, \
                             ts_xyyyz_yzzz, \
                             ts_xyyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyz_xxxx[i] = tk_xyyy_xxxx[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxx[i] * fz_0;

        tk_xyyyz_xxxy[i] = tk_xyyy_xxxy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxy[i] * fz_0;

        tk_xyyyz_xxxz[i] = 3.0 * tk_yyyz_xxz[i] * fe_0 + tk_yyyz_xxxz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxz[i] * fz_0;

        tk_xyyyz_xxyy[i] = tk_xyyy_xxyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxyy[i] * fz_0;

        tk_xyyyz_xxyz[i] = 2.0 * tk_yyyz_xyz[i] * fe_0 + tk_yyyz_xxyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxyz[i] * fz_0;

        tk_xyyyz_xxzz[i] = 2.0 * tk_yyyz_xzz[i] * fe_0 + tk_yyyz_xxzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxzz[i] * fz_0;

        tk_xyyyz_xyyy[i] = tk_xyyy_xyyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xyyy[i] * fz_0;

        tk_xyyyz_xyyz[i] = tk_yyyz_yyz[i] * fe_0 + tk_yyyz_xyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyyz[i] * fz_0;

        tk_xyyyz_xyzz[i] = tk_yyyz_yzz[i] * fe_0 + tk_yyyz_xyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyzz[i] * fz_0;

        tk_xyyyz_xzzz[i] = tk_yyyz_zzz[i] * fe_0 + tk_yyyz_xzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xzzz[i] * fz_0;

        tk_xyyyz_yyyy[i] = tk_yyyz_yyyy[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyy[i] * fz_0;

        tk_xyyyz_yyyz[i] = tk_yyyz_yyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyz[i] * fz_0;

        tk_xyyyz_yyzz[i] = tk_yyyz_yyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyzz[i] * fz_0;

        tk_xyyyz_yzzz[i] = tk_yyyz_yzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yzzz[i] * fz_0;

        tk_xyyyz_zzzz[i] = tk_yyyz_zzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_zzzz[i] * fz_0;
    }

    // Set up 180-195 components of targeted buffer : HG

    auto tk_xyyzz_xxxx = pbuffer.data(idx_kin_hg + 180);

    auto tk_xyyzz_xxxy = pbuffer.data(idx_kin_hg + 181);

    auto tk_xyyzz_xxxz = pbuffer.data(idx_kin_hg + 182);

    auto tk_xyyzz_xxyy = pbuffer.data(idx_kin_hg + 183);

    auto tk_xyyzz_xxyz = pbuffer.data(idx_kin_hg + 184);

    auto tk_xyyzz_xxzz = pbuffer.data(idx_kin_hg + 185);

    auto tk_xyyzz_xyyy = pbuffer.data(idx_kin_hg + 186);

    auto tk_xyyzz_xyyz = pbuffer.data(idx_kin_hg + 187);

    auto tk_xyyzz_xyzz = pbuffer.data(idx_kin_hg + 188);

    auto tk_xyyzz_xzzz = pbuffer.data(idx_kin_hg + 189);

    auto tk_xyyzz_yyyy = pbuffer.data(idx_kin_hg + 190);

    auto tk_xyyzz_yyyz = pbuffer.data(idx_kin_hg + 191);

    auto tk_xyyzz_yyzz = pbuffer.data(idx_kin_hg + 192);

    auto tk_xyyzz_yzzz = pbuffer.data(idx_kin_hg + 193);

    auto tk_xyyzz_zzzz = pbuffer.data(idx_kin_hg + 194);

#pragma omp simd aligned(pa_x,              \
                             tk_xyyzz_xxxx, \
                             tk_xyyzz_xxxy, \
                             tk_xyyzz_xxxz, \
                             tk_xyyzz_xxyy, \
                             tk_xyyzz_xxyz, \
                             tk_xyyzz_xxzz, \
                             tk_xyyzz_xyyy, \
                             tk_xyyzz_xyyz, \
                             tk_xyyzz_xyzz, \
                             tk_xyyzz_xzzz, \
                             tk_xyyzz_yyyy, \
                             tk_xyyzz_yyyz, \
                             tk_xyyzz_yyzz, \
                             tk_xyyzz_yzzz, \
                             tk_xyyzz_zzzz, \
                             tk_yyzz_xxx,   \
                             tk_yyzz_xxxx,  \
                             tk_yyzz_xxxy,  \
                             tk_yyzz_xxxz,  \
                             tk_yyzz_xxy,   \
                             tk_yyzz_xxyy,  \
                             tk_yyzz_xxyz,  \
                             tk_yyzz_xxz,   \
                             tk_yyzz_xxzz,  \
                             tk_yyzz_xyy,   \
                             tk_yyzz_xyyy,  \
                             tk_yyzz_xyyz,  \
                             tk_yyzz_xyz,   \
                             tk_yyzz_xyzz,  \
                             tk_yyzz_xzz,   \
                             tk_yyzz_xzzz,  \
                             tk_yyzz_yyy,   \
                             tk_yyzz_yyyy,  \
                             tk_yyzz_yyyz,  \
                             tk_yyzz_yyz,   \
                             tk_yyzz_yyzz,  \
                             tk_yyzz_yzz,   \
                             tk_yyzz_yzzz,  \
                             tk_yyzz_zzz,   \
                             tk_yyzz_zzzz,  \
                             ts_xyyzz_xxxx, \
                             ts_xyyzz_xxxy, \
                             ts_xyyzz_xxxz, \
                             ts_xyyzz_xxyy, \
                             ts_xyyzz_xxyz, \
                             ts_xyyzz_xxzz, \
                             ts_xyyzz_xyyy, \
                             ts_xyyzz_xyyz, \
                             ts_xyyzz_xyzz, \
                             ts_xyyzz_xzzz, \
                             ts_xyyzz_yyyy, \
                             ts_xyyzz_yyyz, \
                             ts_xyyzz_yyzz, \
                             ts_xyyzz_yzzz, \
                             ts_xyyzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzz_xxxx[i] = 4.0 * tk_yyzz_xxx[i] * fe_0 + tk_yyzz_xxxx[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxx[i] * fz_0;

        tk_xyyzz_xxxy[i] = 3.0 * tk_yyzz_xxy[i] * fe_0 + tk_yyzz_xxxy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxy[i] * fz_0;

        tk_xyyzz_xxxz[i] = 3.0 * tk_yyzz_xxz[i] * fe_0 + tk_yyzz_xxxz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxz[i] * fz_0;

        tk_xyyzz_xxyy[i] = 2.0 * tk_yyzz_xyy[i] * fe_0 + tk_yyzz_xxyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyy[i] * fz_0;

        tk_xyyzz_xxyz[i] = 2.0 * tk_yyzz_xyz[i] * fe_0 + tk_yyzz_xxyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyz[i] * fz_0;

        tk_xyyzz_xxzz[i] = 2.0 * tk_yyzz_xzz[i] * fe_0 + tk_yyzz_xxzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxzz[i] * fz_0;

        tk_xyyzz_xyyy[i] = tk_yyzz_yyy[i] * fe_0 + tk_yyzz_xyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyy[i] * fz_0;

        tk_xyyzz_xyyz[i] = tk_yyzz_yyz[i] * fe_0 + tk_yyzz_xyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyz[i] * fz_0;

        tk_xyyzz_xyzz[i] = tk_yyzz_yzz[i] * fe_0 + tk_yyzz_xyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyzz[i] * fz_0;

        tk_xyyzz_xzzz[i] = tk_yyzz_zzz[i] * fe_0 + tk_yyzz_xzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xzzz[i] * fz_0;

        tk_xyyzz_yyyy[i] = tk_yyzz_yyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyy[i] * fz_0;

        tk_xyyzz_yyyz[i] = tk_yyzz_yyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyz[i] * fz_0;

        tk_xyyzz_yyzz[i] = tk_yyzz_yyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyzz[i] * fz_0;

        tk_xyyzz_yzzz[i] = tk_yyzz_yzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yzzz[i] * fz_0;

        tk_xyyzz_zzzz[i] = tk_yyzz_zzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_zzzz[i] * fz_0;
    }

    // Set up 195-210 components of targeted buffer : HG

    auto tk_xyzzz_xxxx = pbuffer.data(idx_kin_hg + 195);

    auto tk_xyzzz_xxxy = pbuffer.data(idx_kin_hg + 196);

    auto tk_xyzzz_xxxz = pbuffer.data(idx_kin_hg + 197);

    auto tk_xyzzz_xxyy = pbuffer.data(idx_kin_hg + 198);

    auto tk_xyzzz_xxyz = pbuffer.data(idx_kin_hg + 199);

    auto tk_xyzzz_xxzz = pbuffer.data(idx_kin_hg + 200);

    auto tk_xyzzz_xyyy = pbuffer.data(idx_kin_hg + 201);

    auto tk_xyzzz_xyyz = pbuffer.data(idx_kin_hg + 202);

    auto tk_xyzzz_xyzz = pbuffer.data(idx_kin_hg + 203);

    auto tk_xyzzz_xzzz = pbuffer.data(idx_kin_hg + 204);

    auto tk_xyzzz_yyyy = pbuffer.data(idx_kin_hg + 205);

    auto tk_xyzzz_yyyz = pbuffer.data(idx_kin_hg + 206);

    auto tk_xyzzz_yyzz = pbuffer.data(idx_kin_hg + 207);

    auto tk_xyzzz_yzzz = pbuffer.data(idx_kin_hg + 208);

    auto tk_xyzzz_zzzz = pbuffer.data(idx_kin_hg + 209);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xyzzz_xxxx, \
                             tk_xyzzz_xxxy, \
                             tk_xyzzz_xxxz, \
                             tk_xyzzz_xxyy, \
                             tk_xyzzz_xxyz, \
                             tk_xyzzz_xxzz, \
                             tk_xyzzz_xyyy, \
                             tk_xyzzz_xyyz, \
                             tk_xyzzz_xyzz, \
                             tk_xyzzz_xzzz, \
                             tk_xyzzz_yyyy, \
                             tk_xyzzz_yyyz, \
                             tk_xyzzz_yyzz, \
                             tk_xyzzz_yzzz, \
                             tk_xyzzz_zzzz, \
                             tk_xzzz_xxxx,  \
                             tk_xzzz_xxxz,  \
                             tk_xzzz_xxzz,  \
                             tk_xzzz_xzzz,  \
                             tk_yzzz_xxxy,  \
                             tk_yzzz_xxy,   \
                             tk_yzzz_xxyy,  \
                             tk_yzzz_xxyz,  \
                             tk_yzzz_xyy,   \
                             tk_yzzz_xyyy,  \
                             tk_yzzz_xyyz,  \
                             tk_yzzz_xyz,   \
                             tk_yzzz_xyzz,  \
                             tk_yzzz_yyy,   \
                             tk_yzzz_yyyy,  \
                             tk_yzzz_yyyz,  \
                             tk_yzzz_yyz,   \
                             tk_yzzz_yyzz,  \
                             tk_yzzz_yzz,   \
                             tk_yzzz_yzzz,  \
                             tk_yzzz_zzzz,  \
                             ts_xyzzz_xxxx, \
                             ts_xyzzz_xxxy, \
                             ts_xyzzz_xxxz, \
                             ts_xyzzz_xxyy, \
                             ts_xyzzz_xxyz, \
                             ts_xyzzz_xxzz, \
                             ts_xyzzz_xyyy, \
                             ts_xyzzz_xyyz, \
                             ts_xyzzz_xyzz, \
                             ts_xyzzz_xzzz, \
                             ts_xyzzz_yyyy, \
                             ts_xyzzz_yyyz, \
                             ts_xyzzz_yyzz, \
                             ts_xyzzz_yzzz, \
                             ts_xyzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzz_xxxx[i] = tk_xzzz_xxxx[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxx[i] * fz_0;

        tk_xyzzz_xxxy[i] = 3.0 * tk_yzzz_xxy[i] * fe_0 + tk_yzzz_xxxy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxy[i] * fz_0;

        tk_xyzzz_xxxz[i] = tk_xzzz_xxxz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxz[i] * fz_0;

        tk_xyzzz_xxyy[i] = 2.0 * tk_yzzz_xyy[i] * fe_0 + tk_yzzz_xxyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyy[i] * fz_0;

        tk_xyzzz_xxyz[i] = 2.0 * tk_yzzz_xyz[i] * fe_0 + tk_yzzz_xxyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyz[i] * fz_0;

        tk_xyzzz_xxzz[i] = tk_xzzz_xxzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxzz[i] * fz_0;

        tk_xyzzz_xyyy[i] = tk_yzzz_yyy[i] * fe_0 + tk_yzzz_xyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyy[i] * fz_0;

        tk_xyzzz_xyyz[i] = tk_yzzz_yyz[i] * fe_0 + tk_yzzz_xyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyz[i] * fz_0;

        tk_xyzzz_xyzz[i] = tk_yzzz_yzz[i] * fe_0 + tk_yzzz_xyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyzz[i] * fz_0;

        tk_xyzzz_xzzz[i] = tk_xzzz_xzzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xzzz[i] * fz_0;

        tk_xyzzz_yyyy[i] = tk_yzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyy[i] * fz_0;

        tk_xyzzz_yyyz[i] = tk_yzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyz[i] * fz_0;

        tk_xyzzz_yyzz[i] = tk_yzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyzz[i] * fz_0;

        tk_xyzzz_yzzz[i] = tk_yzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yzzz[i] * fz_0;

        tk_xyzzz_zzzz[i] = tk_yzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_zzzz[i] * fz_0;
    }

    // Set up 210-225 components of targeted buffer : HG

    auto tk_xzzzz_xxxx = pbuffer.data(idx_kin_hg + 210);

    auto tk_xzzzz_xxxy = pbuffer.data(idx_kin_hg + 211);

    auto tk_xzzzz_xxxz = pbuffer.data(idx_kin_hg + 212);

    auto tk_xzzzz_xxyy = pbuffer.data(idx_kin_hg + 213);

    auto tk_xzzzz_xxyz = pbuffer.data(idx_kin_hg + 214);

    auto tk_xzzzz_xxzz = pbuffer.data(idx_kin_hg + 215);

    auto tk_xzzzz_xyyy = pbuffer.data(idx_kin_hg + 216);

    auto tk_xzzzz_xyyz = pbuffer.data(idx_kin_hg + 217);

    auto tk_xzzzz_xyzz = pbuffer.data(idx_kin_hg + 218);

    auto tk_xzzzz_xzzz = pbuffer.data(idx_kin_hg + 219);

    auto tk_xzzzz_yyyy = pbuffer.data(idx_kin_hg + 220);

    auto tk_xzzzz_yyyz = pbuffer.data(idx_kin_hg + 221);

    auto tk_xzzzz_yyzz = pbuffer.data(idx_kin_hg + 222);

    auto tk_xzzzz_yzzz = pbuffer.data(idx_kin_hg + 223);

    auto tk_xzzzz_zzzz = pbuffer.data(idx_kin_hg + 224);

#pragma omp simd aligned(pa_x,              \
                             tk_xzzzz_xxxx, \
                             tk_xzzzz_xxxy, \
                             tk_xzzzz_xxxz, \
                             tk_xzzzz_xxyy, \
                             tk_xzzzz_xxyz, \
                             tk_xzzzz_xxzz, \
                             tk_xzzzz_xyyy, \
                             tk_xzzzz_xyyz, \
                             tk_xzzzz_xyzz, \
                             tk_xzzzz_xzzz, \
                             tk_xzzzz_yyyy, \
                             tk_xzzzz_yyyz, \
                             tk_xzzzz_yyzz, \
                             tk_xzzzz_yzzz, \
                             tk_xzzzz_zzzz, \
                             tk_zzzz_xxx,   \
                             tk_zzzz_xxxx,  \
                             tk_zzzz_xxxy,  \
                             tk_zzzz_xxxz,  \
                             tk_zzzz_xxy,   \
                             tk_zzzz_xxyy,  \
                             tk_zzzz_xxyz,  \
                             tk_zzzz_xxz,   \
                             tk_zzzz_xxzz,  \
                             tk_zzzz_xyy,   \
                             tk_zzzz_xyyy,  \
                             tk_zzzz_xyyz,  \
                             tk_zzzz_xyz,   \
                             tk_zzzz_xyzz,  \
                             tk_zzzz_xzz,   \
                             tk_zzzz_xzzz,  \
                             tk_zzzz_yyy,   \
                             tk_zzzz_yyyy,  \
                             tk_zzzz_yyyz,  \
                             tk_zzzz_yyz,   \
                             tk_zzzz_yyzz,  \
                             tk_zzzz_yzz,   \
                             tk_zzzz_yzzz,  \
                             tk_zzzz_zzz,   \
                             tk_zzzz_zzzz,  \
                             ts_xzzzz_xxxx, \
                             ts_xzzzz_xxxy, \
                             ts_xzzzz_xxxz, \
                             ts_xzzzz_xxyy, \
                             ts_xzzzz_xxyz, \
                             ts_xzzzz_xxzz, \
                             ts_xzzzz_xyyy, \
                             ts_xzzzz_xyyz, \
                             ts_xzzzz_xyzz, \
                             ts_xzzzz_xzzz, \
                             ts_xzzzz_yyyy, \
                             ts_xzzzz_yyyz, \
                             ts_xzzzz_yyzz, \
                             ts_xzzzz_yzzz, \
                             ts_xzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzz_xxxx[i] = 4.0 * tk_zzzz_xxx[i] * fe_0 + tk_zzzz_xxxx[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxx[i] * fz_0;

        tk_xzzzz_xxxy[i] = 3.0 * tk_zzzz_xxy[i] * fe_0 + tk_zzzz_xxxy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxy[i] * fz_0;

        tk_xzzzz_xxxz[i] = 3.0 * tk_zzzz_xxz[i] * fe_0 + tk_zzzz_xxxz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxz[i] * fz_0;

        tk_xzzzz_xxyy[i] = 2.0 * tk_zzzz_xyy[i] * fe_0 + tk_zzzz_xxyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyy[i] * fz_0;

        tk_xzzzz_xxyz[i] = 2.0 * tk_zzzz_xyz[i] * fe_0 + tk_zzzz_xxyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyz[i] * fz_0;

        tk_xzzzz_xxzz[i] = 2.0 * tk_zzzz_xzz[i] * fe_0 + tk_zzzz_xxzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxzz[i] * fz_0;

        tk_xzzzz_xyyy[i] = tk_zzzz_yyy[i] * fe_0 + tk_zzzz_xyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyy[i] * fz_0;

        tk_xzzzz_xyyz[i] = tk_zzzz_yyz[i] * fe_0 + tk_zzzz_xyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyz[i] * fz_0;

        tk_xzzzz_xyzz[i] = tk_zzzz_yzz[i] * fe_0 + tk_zzzz_xyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyzz[i] * fz_0;

        tk_xzzzz_xzzz[i] = tk_zzzz_zzz[i] * fe_0 + tk_zzzz_xzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xzzz[i] * fz_0;

        tk_xzzzz_yyyy[i] = tk_zzzz_yyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyy[i] * fz_0;

        tk_xzzzz_yyyz[i] = tk_zzzz_yyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyz[i] * fz_0;

        tk_xzzzz_yyzz[i] = tk_zzzz_yyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyzz[i] * fz_0;

        tk_xzzzz_yzzz[i] = tk_zzzz_yzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yzzz[i] * fz_0;

        tk_xzzzz_zzzz[i] = tk_zzzz_zzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_zzzz[i] * fz_0;
    }

    // Set up 225-240 components of targeted buffer : HG

    auto tk_yyyyy_xxxx = pbuffer.data(idx_kin_hg + 225);

    auto tk_yyyyy_xxxy = pbuffer.data(idx_kin_hg + 226);

    auto tk_yyyyy_xxxz = pbuffer.data(idx_kin_hg + 227);

    auto tk_yyyyy_xxyy = pbuffer.data(idx_kin_hg + 228);

    auto tk_yyyyy_xxyz = pbuffer.data(idx_kin_hg + 229);

    auto tk_yyyyy_xxzz = pbuffer.data(idx_kin_hg + 230);

    auto tk_yyyyy_xyyy = pbuffer.data(idx_kin_hg + 231);

    auto tk_yyyyy_xyyz = pbuffer.data(idx_kin_hg + 232);

    auto tk_yyyyy_xyzz = pbuffer.data(idx_kin_hg + 233);

    auto tk_yyyyy_xzzz = pbuffer.data(idx_kin_hg + 234);

    auto tk_yyyyy_yyyy = pbuffer.data(idx_kin_hg + 235);

    auto tk_yyyyy_yyyz = pbuffer.data(idx_kin_hg + 236);

    auto tk_yyyyy_yyzz = pbuffer.data(idx_kin_hg + 237);

    auto tk_yyyyy_yzzz = pbuffer.data(idx_kin_hg + 238);

    auto tk_yyyyy_zzzz = pbuffer.data(idx_kin_hg + 239);

#pragma omp simd aligned(pa_y,              \
                             tk_yyy_xxxx,   \
                             tk_yyy_xxxy,   \
                             tk_yyy_xxxz,   \
                             tk_yyy_xxyy,   \
                             tk_yyy_xxyz,   \
                             tk_yyy_xxzz,   \
                             tk_yyy_xyyy,   \
                             tk_yyy_xyyz,   \
                             tk_yyy_xyzz,   \
                             tk_yyy_xzzz,   \
                             tk_yyy_yyyy,   \
                             tk_yyy_yyyz,   \
                             tk_yyy_yyzz,   \
                             tk_yyy_yzzz,   \
                             tk_yyy_zzzz,   \
                             tk_yyyy_xxx,   \
                             tk_yyyy_xxxx,  \
                             tk_yyyy_xxxy,  \
                             tk_yyyy_xxxz,  \
                             tk_yyyy_xxy,   \
                             tk_yyyy_xxyy,  \
                             tk_yyyy_xxyz,  \
                             tk_yyyy_xxz,   \
                             tk_yyyy_xxzz,  \
                             tk_yyyy_xyy,   \
                             tk_yyyy_xyyy,  \
                             tk_yyyy_xyyz,  \
                             tk_yyyy_xyz,   \
                             tk_yyyy_xyzz,  \
                             tk_yyyy_xzz,   \
                             tk_yyyy_xzzz,  \
                             tk_yyyy_yyy,   \
                             tk_yyyy_yyyy,  \
                             tk_yyyy_yyyz,  \
                             tk_yyyy_yyz,   \
                             tk_yyyy_yyzz,  \
                             tk_yyyy_yzz,   \
                             tk_yyyy_yzzz,  \
                             tk_yyyy_zzz,   \
                             tk_yyyy_zzzz,  \
                             tk_yyyyy_xxxx, \
                             tk_yyyyy_xxxy, \
                             tk_yyyyy_xxxz, \
                             tk_yyyyy_xxyy, \
                             tk_yyyyy_xxyz, \
                             tk_yyyyy_xxzz, \
                             tk_yyyyy_xyyy, \
                             tk_yyyyy_xyyz, \
                             tk_yyyyy_xyzz, \
                             tk_yyyyy_xzzz, \
                             tk_yyyyy_yyyy, \
                             tk_yyyyy_yyyz, \
                             tk_yyyyy_yyzz, \
                             tk_yyyyy_yzzz, \
                             tk_yyyyy_zzzz, \
                             ts_yyy_xxxx,   \
                             ts_yyy_xxxy,   \
                             ts_yyy_xxxz,   \
                             ts_yyy_xxyy,   \
                             ts_yyy_xxyz,   \
                             ts_yyy_xxzz,   \
                             ts_yyy_xyyy,   \
                             ts_yyy_xyyz,   \
                             ts_yyy_xyzz,   \
                             ts_yyy_xzzz,   \
                             ts_yyy_yyyy,   \
                             ts_yyy_yyyz,   \
                             ts_yyy_yyzz,   \
                             ts_yyy_yzzz,   \
                             ts_yyy_zzzz,   \
                             ts_yyyyy_xxxx, \
                             ts_yyyyy_xxxy, \
                             ts_yyyyy_xxxz, \
                             ts_yyyyy_xxyy, \
                             ts_yyyyy_xxyz, \
                             ts_yyyyy_xxzz, \
                             ts_yyyyy_xyyy, \
                             ts_yyyyy_xyyz, \
                             ts_yyyyy_xyzz, \
                             ts_yyyyy_xzzz, \
                             ts_yyyyy_yyyy, \
                             ts_yyyyy_yyyz, \
                             ts_yyyyy_yyzz, \
                             ts_yyyyy_yzzz, \
                             ts_yyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyy_xxxx[i] =
            -8.0 * ts_yyy_xxxx[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxx[i] * fe_0 + tk_yyyy_xxxx[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxx[i] * fz_0;

        tk_yyyyy_xxxy[i] = -8.0 * ts_yyy_xxxy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxy[i] * fe_0 + tk_yyyy_xxx[i] * fe_0 + tk_yyyy_xxxy[i] * pa_y[i] +
                           2.0 * ts_yyyyy_xxxy[i] * fz_0;

        tk_yyyyy_xxxz[i] =
            -8.0 * ts_yyy_xxxz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxz[i] * fe_0 + tk_yyyy_xxxz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxz[i] * fz_0;

        tk_yyyyy_xxyy[i] = -8.0 * ts_yyy_xxyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyy[i] * fe_0 + 2.0 * tk_yyyy_xxy[i] * fe_0 +
                           tk_yyyy_xxyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyy[i] * fz_0;

        tk_yyyyy_xxyz[i] = -8.0 * ts_yyy_xxyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyz[i] * fe_0 + tk_yyyy_xxz[i] * fe_0 + tk_yyyy_xxyz[i] * pa_y[i] +
                           2.0 * ts_yyyyy_xxyz[i] * fz_0;

        tk_yyyyy_xxzz[i] =
            -8.0 * ts_yyy_xxzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxzz[i] * fe_0 + tk_yyyy_xxzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxzz[i] * fz_0;

        tk_yyyyy_xyyy[i] = -8.0 * ts_yyy_xyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyy[i] * fe_0 + 3.0 * tk_yyyy_xyy[i] * fe_0 +
                           tk_yyyy_xyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyy[i] * fz_0;

        tk_yyyyy_xyyz[i] = -8.0 * ts_yyy_xyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyz[i] * fe_0 + 2.0 * tk_yyyy_xyz[i] * fe_0 +
                           tk_yyyy_xyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyz[i] * fz_0;

        tk_yyyyy_xyzz[i] = -8.0 * ts_yyy_xyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyzz[i] * fe_0 + tk_yyyy_xzz[i] * fe_0 + tk_yyyy_xyzz[i] * pa_y[i] +
                           2.0 * ts_yyyyy_xyzz[i] * fz_0;

        tk_yyyyy_xzzz[i] =
            -8.0 * ts_yyy_xzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xzzz[i] * fe_0 + tk_yyyy_xzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xzzz[i] * fz_0;

        tk_yyyyy_yyyy[i] = -8.0 * ts_yyy_yyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyy[i] * fe_0 + 4.0 * tk_yyyy_yyy[i] * fe_0 +
                           tk_yyyy_yyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyy[i] * fz_0;

        tk_yyyyy_yyyz[i] = -8.0 * ts_yyy_yyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyz[i] * fe_0 + 3.0 * tk_yyyy_yyz[i] * fe_0 +
                           tk_yyyy_yyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyz[i] * fz_0;

        tk_yyyyy_yyzz[i] = -8.0 * ts_yyy_yyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyzz[i] * fe_0 + 2.0 * tk_yyyy_yzz[i] * fe_0 +
                           tk_yyyy_yyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyzz[i] * fz_0;

        tk_yyyyy_yzzz[i] = -8.0 * ts_yyy_yzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yzzz[i] * fe_0 + tk_yyyy_zzz[i] * fe_0 + tk_yyyy_yzzz[i] * pa_y[i] +
                           2.0 * ts_yyyyy_yzzz[i] * fz_0;

        tk_yyyyy_zzzz[i] =
            -8.0 * ts_yyy_zzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_zzzz[i] * fe_0 + tk_yyyy_zzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_zzzz[i] * fz_0;
    }

    // Set up 240-255 components of targeted buffer : HG

    auto tk_yyyyz_xxxx = pbuffer.data(idx_kin_hg + 240);

    auto tk_yyyyz_xxxy = pbuffer.data(idx_kin_hg + 241);

    auto tk_yyyyz_xxxz = pbuffer.data(idx_kin_hg + 242);

    auto tk_yyyyz_xxyy = pbuffer.data(idx_kin_hg + 243);

    auto tk_yyyyz_xxyz = pbuffer.data(idx_kin_hg + 244);

    auto tk_yyyyz_xxzz = pbuffer.data(idx_kin_hg + 245);

    auto tk_yyyyz_xyyy = pbuffer.data(idx_kin_hg + 246);

    auto tk_yyyyz_xyyz = pbuffer.data(idx_kin_hg + 247);

    auto tk_yyyyz_xyzz = pbuffer.data(idx_kin_hg + 248);

    auto tk_yyyyz_xzzz = pbuffer.data(idx_kin_hg + 249);

    auto tk_yyyyz_yyyy = pbuffer.data(idx_kin_hg + 250);

    auto tk_yyyyz_yyyz = pbuffer.data(idx_kin_hg + 251);

    auto tk_yyyyz_yyzz = pbuffer.data(idx_kin_hg + 252);

    auto tk_yyyyz_yzzz = pbuffer.data(idx_kin_hg + 253);

    auto tk_yyyyz_zzzz = pbuffer.data(idx_kin_hg + 254);

#pragma omp simd aligned(pa_z,              \
                             tk_yyyy_xxx,   \
                             tk_yyyy_xxxx,  \
                             tk_yyyy_xxxy,  \
                             tk_yyyy_xxxz,  \
                             tk_yyyy_xxy,   \
                             tk_yyyy_xxyy,  \
                             tk_yyyy_xxyz,  \
                             tk_yyyy_xxz,   \
                             tk_yyyy_xxzz,  \
                             tk_yyyy_xyy,   \
                             tk_yyyy_xyyy,  \
                             tk_yyyy_xyyz,  \
                             tk_yyyy_xyz,   \
                             tk_yyyy_xyzz,  \
                             tk_yyyy_xzz,   \
                             tk_yyyy_xzzz,  \
                             tk_yyyy_yyy,   \
                             tk_yyyy_yyyy,  \
                             tk_yyyy_yyyz,  \
                             tk_yyyy_yyz,   \
                             tk_yyyy_yyzz,  \
                             tk_yyyy_yzz,   \
                             tk_yyyy_yzzz,  \
                             tk_yyyy_zzz,   \
                             tk_yyyy_zzzz,  \
                             tk_yyyyz_xxxx, \
                             tk_yyyyz_xxxy, \
                             tk_yyyyz_xxxz, \
                             tk_yyyyz_xxyy, \
                             tk_yyyyz_xxyz, \
                             tk_yyyyz_xxzz, \
                             tk_yyyyz_xyyy, \
                             tk_yyyyz_xyyz, \
                             tk_yyyyz_xyzz, \
                             tk_yyyyz_xzzz, \
                             tk_yyyyz_yyyy, \
                             tk_yyyyz_yyyz, \
                             tk_yyyyz_yyzz, \
                             tk_yyyyz_yzzz, \
                             tk_yyyyz_zzzz, \
                             ts_yyyyz_xxxx, \
                             ts_yyyyz_xxxy, \
                             ts_yyyyz_xxxz, \
                             ts_yyyyz_xxyy, \
                             ts_yyyyz_xxyz, \
                             ts_yyyyz_xxzz, \
                             ts_yyyyz_xyyy, \
                             ts_yyyyz_xyyz, \
                             ts_yyyyz_xyzz, \
                             ts_yyyyz_xzzz, \
                             ts_yyyyz_yyyy, \
                             ts_yyyyz_yyyz, \
                             ts_yyyyz_yyzz, \
                             ts_yyyyz_yzzz, \
                             ts_yyyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyz_xxxx[i] = tk_yyyy_xxxx[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxx[i] * fz_0;

        tk_yyyyz_xxxy[i] = tk_yyyy_xxxy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxy[i] * fz_0;

        tk_yyyyz_xxxz[i] = tk_yyyy_xxx[i] * fe_0 + tk_yyyy_xxxz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxz[i] * fz_0;

        tk_yyyyz_xxyy[i] = tk_yyyy_xxyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyy[i] * fz_0;

        tk_yyyyz_xxyz[i] = tk_yyyy_xxy[i] * fe_0 + tk_yyyy_xxyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyz[i] * fz_0;

        tk_yyyyz_xxzz[i] = 2.0 * tk_yyyy_xxz[i] * fe_0 + tk_yyyy_xxzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxzz[i] * fz_0;

        tk_yyyyz_xyyy[i] = tk_yyyy_xyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyy[i] * fz_0;

        tk_yyyyz_xyyz[i] = tk_yyyy_xyy[i] * fe_0 + tk_yyyy_xyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyz[i] * fz_0;

        tk_yyyyz_xyzz[i] = 2.0 * tk_yyyy_xyz[i] * fe_0 + tk_yyyy_xyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyzz[i] * fz_0;

        tk_yyyyz_xzzz[i] = 3.0 * tk_yyyy_xzz[i] * fe_0 + tk_yyyy_xzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xzzz[i] * fz_0;

        tk_yyyyz_yyyy[i] = tk_yyyy_yyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyy[i] * fz_0;

        tk_yyyyz_yyyz[i] = tk_yyyy_yyy[i] * fe_0 + tk_yyyy_yyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyz[i] * fz_0;

        tk_yyyyz_yyzz[i] = 2.0 * tk_yyyy_yyz[i] * fe_0 + tk_yyyy_yyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyzz[i] * fz_0;

        tk_yyyyz_yzzz[i] = 3.0 * tk_yyyy_yzz[i] * fe_0 + tk_yyyy_yzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yzzz[i] * fz_0;

        tk_yyyyz_zzzz[i] = 4.0 * tk_yyyy_zzz[i] * fe_0 + tk_yyyy_zzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_zzzz[i] * fz_0;
    }

    // Set up 255-270 components of targeted buffer : HG

    auto tk_yyyzz_xxxx = pbuffer.data(idx_kin_hg + 255);

    auto tk_yyyzz_xxxy = pbuffer.data(idx_kin_hg + 256);

    auto tk_yyyzz_xxxz = pbuffer.data(idx_kin_hg + 257);

    auto tk_yyyzz_xxyy = pbuffer.data(idx_kin_hg + 258);

    auto tk_yyyzz_xxyz = pbuffer.data(idx_kin_hg + 259);

    auto tk_yyyzz_xxzz = pbuffer.data(idx_kin_hg + 260);

    auto tk_yyyzz_xyyy = pbuffer.data(idx_kin_hg + 261);

    auto tk_yyyzz_xyyz = pbuffer.data(idx_kin_hg + 262);

    auto tk_yyyzz_xyzz = pbuffer.data(idx_kin_hg + 263);

    auto tk_yyyzz_xzzz = pbuffer.data(idx_kin_hg + 264);

    auto tk_yyyzz_yyyy = pbuffer.data(idx_kin_hg + 265);

    auto tk_yyyzz_yyyz = pbuffer.data(idx_kin_hg + 266);

    auto tk_yyyzz_yyzz = pbuffer.data(idx_kin_hg + 267);

    auto tk_yyyzz_yzzz = pbuffer.data(idx_kin_hg + 268);

    auto tk_yyyzz_zzzz = pbuffer.data(idx_kin_hg + 269);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_yyy_xxxy,   \
                             tk_yyy_xxyy,   \
                             tk_yyy_xyyy,   \
                             tk_yyy_yyyy,   \
                             tk_yyyz_xxxy,  \
                             tk_yyyz_xxyy,  \
                             tk_yyyz_xyyy,  \
                             tk_yyyz_yyyy,  \
                             tk_yyyzz_xxxx, \
                             tk_yyyzz_xxxy, \
                             tk_yyyzz_xxxz, \
                             tk_yyyzz_xxyy, \
                             tk_yyyzz_xxyz, \
                             tk_yyyzz_xxzz, \
                             tk_yyyzz_xyyy, \
                             tk_yyyzz_xyyz, \
                             tk_yyyzz_xyzz, \
                             tk_yyyzz_xzzz, \
                             tk_yyyzz_yyyy, \
                             tk_yyyzz_yyyz, \
                             tk_yyyzz_yyzz, \
                             tk_yyyzz_yzzz, \
                             tk_yyyzz_zzzz, \
                             tk_yyzz_xxxx,  \
                             tk_yyzz_xxxz,  \
                             tk_yyzz_xxyz,  \
                             tk_yyzz_xxz,   \
                             tk_yyzz_xxzz,  \
                             tk_yyzz_xyyz,  \
                             tk_yyzz_xyz,   \
                             tk_yyzz_xyzz,  \
                             tk_yyzz_xzz,   \
                             tk_yyzz_xzzz,  \
                             tk_yyzz_yyyz,  \
                             tk_yyzz_yyz,   \
                             tk_yyzz_yyzz,  \
                             tk_yyzz_yzz,   \
                             tk_yyzz_yzzz,  \
                             tk_yyzz_zzz,   \
                             tk_yyzz_zzzz,  \
                             tk_yzz_xxxx,   \
                             tk_yzz_xxxz,   \
                             tk_yzz_xxyz,   \
                             tk_yzz_xxzz,   \
                             tk_yzz_xyyz,   \
                             tk_yzz_xyzz,   \
                             tk_yzz_xzzz,   \
                             tk_yzz_yyyz,   \
                             tk_yzz_yyzz,   \
                             tk_yzz_yzzz,   \
                             tk_yzz_zzzz,   \
                             ts_yyy_xxxy,   \
                             ts_yyy_xxyy,   \
                             ts_yyy_xyyy,   \
                             ts_yyy_yyyy,   \
                             ts_yyyzz_xxxx, \
                             ts_yyyzz_xxxy, \
                             ts_yyyzz_xxxz, \
                             ts_yyyzz_xxyy, \
                             ts_yyyzz_xxyz, \
                             ts_yyyzz_xxzz, \
                             ts_yyyzz_xyyy, \
                             ts_yyyzz_xyyz, \
                             ts_yyyzz_xyzz, \
                             ts_yyyzz_xzzz, \
                             ts_yyyzz_yyyy, \
                             ts_yyyzz_yyyz, \
                             ts_yyyzz_yyzz, \
                             ts_yyyzz_yzzz, \
                             ts_yyyzz_zzzz, \
                             ts_yzz_xxxx,   \
                             ts_yzz_xxxz,   \
                             ts_yzz_xxyz,   \
                             ts_yzz_xxzz,   \
                             ts_yzz_xyyz,   \
                             ts_yzz_xyzz,   \
                             ts_yzz_xzzz,   \
                             ts_yzz_yyyz,   \
                             ts_yzz_yyzz,   \
                             ts_yzz_yzzz,   \
                             ts_yzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzz_xxxx[i] =
            -4.0 * ts_yzz_xxxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxx[i] * fe_0 + tk_yyzz_xxxx[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxx[i] * fz_0;

        tk_yyyzz_xxxy[i] = -2.0 * ts_yyy_xxxy[i] * fbe_0 * fz_0 + tk_yyy_xxxy[i] * fe_0 + tk_yyyz_xxxy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxxy[i] * fz_0;

        tk_yyyzz_xxxz[i] =
            -4.0 * ts_yzz_xxxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxz[i] * fe_0 + tk_yyzz_xxxz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxz[i] * fz_0;

        tk_yyyzz_xxyy[i] = -2.0 * ts_yyy_xxyy[i] * fbe_0 * fz_0 + tk_yyy_xxyy[i] * fe_0 + tk_yyyz_xxyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxyy[i] * fz_0;

        tk_yyyzz_xxyz[i] = -4.0 * ts_yzz_xxyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxyz[i] * fe_0 + tk_yyzz_xxz[i] * fe_0 + tk_yyzz_xxyz[i] * pa_y[i] +
                           2.0 * ts_yyyzz_xxyz[i] * fz_0;

        tk_yyyzz_xxzz[i] =
            -4.0 * ts_yzz_xxzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxzz[i] * fe_0 + tk_yyzz_xxzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxzz[i] * fz_0;

        tk_yyyzz_xyyy[i] = -2.0 * ts_yyy_xyyy[i] * fbe_0 * fz_0 + tk_yyy_xyyy[i] * fe_0 + tk_yyyz_xyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xyyy[i] * fz_0;

        tk_yyyzz_xyyz[i] = -4.0 * ts_yzz_xyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyyz[i] * fe_0 + 2.0 * tk_yyzz_xyz[i] * fe_0 +
                           tk_yyzz_xyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyyz[i] * fz_0;

        tk_yyyzz_xyzz[i] = -4.0 * ts_yzz_xyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyzz[i] * fe_0 + tk_yyzz_xzz[i] * fe_0 + tk_yyzz_xyzz[i] * pa_y[i] +
                           2.0 * ts_yyyzz_xyzz[i] * fz_0;

        tk_yyyzz_xzzz[i] =
            -4.0 * ts_yzz_xzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xzzz[i] * fe_0 + tk_yyzz_xzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xzzz[i] * fz_0;

        tk_yyyzz_yyyy[i] = -2.0 * ts_yyy_yyyy[i] * fbe_0 * fz_0 + tk_yyy_yyyy[i] * fe_0 + tk_yyyz_yyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_yyyy[i] * fz_0;

        tk_yyyzz_yyyz[i] = -4.0 * ts_yzz_yyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyyz[i] * fe_0 + 3.0 * tk_yyzz_yyz[i] * fe_0 +
                           tk_yyzz_yyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyyz[i] * fz_0;

        tk_yyyzz_yyzz[i] = -4.0 * ts_yzz_yyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyzz[i] * fe_0 + 2.0 * tk_yyzz_yzz[i] * fe_0 +
                           tk_yyzz_yyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyzz[i] * fz_0;

        tk_yyyzz_yzzz[i] = -4.0 * ts_yzz_yzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yzzz[i] * fe_0 + tk_yyzz_zzz[i] * fe_0 + tk_yyzz_yzzz[i] * pa_y[i] +
                           2.0 * ts_yyyzz_yzzz[i] * fz_0;

        tk_yyyzz_zzzz[i] =
            -4.0 * ts_yzz_zzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_zzzz[i] * fe_0 + tk_yyzz_zzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_zzzz[i] * fz_0;
    }

    // Set up 270-285 components of targeted buffer : HG

    auto tk_yyzzz_xxxx = pbuffer.data(idx_kin_hg + 270);

    auto tk_yyzzz_xxxy = pbuffer.data(idx_kin_hg + 271);

    auto tk_yyzzz_xxxz = pbuffer.data(idx_kin_hg + 272);

    auto tk_yyzzz_xxyy = pbuffer.data(idx_kin_hg + 273);

    auto tk_yyzzz_xxyz = pbuffer.data(idx_kin_hg + 274);

    auto tk_yyzzz_xxzz = pbuffer.data(idx_kin_hg + 275);

    auto tk_yyzzz_xyyy = pbuffer.data(idx_kin_hg + 276);

    auto tk_yyzzz_xyyz = pbuffer.data(idx_kin_hg + 277);

    auto tk_yyzzz_xyzz = pbuffer.data(idx_kin_hg + 278);

    auto tk_yyzzz_xzzz = pbuffer.data(idx_kin_hg + 279);

    auto tk_yyzzz_yyyy = pbuffer.data(idx_kin_hg + 280);

    auto tk_yyzzz_yyyz = pbuffer.data(idx_kin_hg + 281);

    auto tk_yyzzz_yyzz = pbuffer.data(idx_kin_hg + 282);

    auto tk_yyzzz_yzzz = pbuffer.data(idx_kin_hg + 283);

    auto tk_yyzzz_zzzz = pbuffer.data(idx_kin_hg + 284);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_yyz_xxxy,   \
                             tk_yyz_xxyy,   \
                             tk_yyz_xyyy,   \
                             tk_yyz_yyyy,   \
                             tk_yyzz_xxxy,  \
                             tk_yyzz_xxyy,  \
                             tk_yyzz_xyyy,  \
                             tk_yyzz_yyyy,  \
                             tk_yyzzz_xxxx, \
                             tk_yyzzz_xxxy, \
                             tk_yyzzz_xxxz, \
                             tk_yyzzz_xxyy, \
                             tk_yyzzz_xxyz, \
                             tk_yyzzz_xxzz, \
                             tk_yyzzz_xyyy, \
                             tk_yyzzz_xyyz, \
                             tk_yyzzz_xyzz, \
                             tk_yyzzz_xzzz, \
                             tk_yyzzz_yyyy, \
                             tk_yyzzz_yyyz, \
                             tk_yyzzz_yyzz, \
                             tk_yyzzz_yzzz, \
                             tk_yyzzz_zzzz, \
                             tk_yzzz_xxxx,  \
                             tk_yzzz_xxxz,  \
                             tk_yzzz_xxyz,  \
                             tk_yzzz_xxz,   \
                             tk_yzzz_xxzz,  \
                             tk_yzzz_xyyz,  \
                             tk_yzzz_xyz,   \
                             tk_yzzz_xyzz,  \
                             tk_yzzz_xzz,   \
                             tk_yzzz_xzzz,  \
                             tk_yzzz_yyyz,  \
                             tk_yzzz_yyz,   \
                             tk_yzzz_yyzz,  \
                             tk_yzzz_yzz,   \
                             tk_yzzz_yzzz,  \
                             tk_yzzz_zzz,   \
                             tk_yzzz_zzzz,  \
                             tk_zzz_xxxx,   \
                             tk_zzz_xxxz,   \
                             tk_zzz_xxyz,   \
                             tk_zzz_xxzz,   \
                             tk_zzz_xyyz,   \
                             tk_zzz_xyzz,   \
                             tk_zzz_xzzz,   \
                             tk_zzz_yyyz,   \
                             tk_zzz_yyzz,   \
                             tk_zzz_yzzz,   \
                             tk_zzz_zzzz,   \
                             ts_yyz_xxxy,   \
                             ts_yyz_xxyy,   \
                             ts_yyz_xyyy,   \
                             ts_yyz_yyyy,   \
                             ts_yyzzz_xxxx, \
                             ts_yyzzz_xxxy, \
                             ts_yyzzz_xxxz, \
                             ts_yyzzz_xxyy, \
                             ts_yyzzz_xxyz, \
                             ts_yyzzz_xxzz, \
                             ts_yyzzz_xyyy, \
                             ts_yyzzz_xyyz, \
                             ts_yyzzz_xyzz, \
                             ts_yyzzz_xzzz, \
                             ts_yyzzz_yyyy, \
                             ts_yyzzz_yyyz, \
                             ts_yyzzz_yyzz, \
                             ts_yyzzz_yzzz, \
                             ts_yyzzz_zzzz, \
                             ts_zzz_xxxx,   \
                             ts_zzz_xxxz,   \
                             ts_zzz_xxyz,   \
                             ts_zzz_xxzz,   \
                             ts_zzz_xyyz,   \
                             ts_zzz_xyzz,   \
                             ts_zzz_xzzz,   \
                             ts_zzz_yyyz,   \
                             ts_zzz_yyzz,   \
                             ts_zzz_yzzz,   \
                             ts_zzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzz_xxxx[i] = -2.0 * ts_zzz_xxxx[i] * fbe_0 * fz_0 + tk_zzz_xxxx[i] * fe_0 + tk_yzzz_xxxx[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxx[i] * fz_0;

        tk_yyzzz_xxxy[i] =
            -4.0 * ts_yyz_xxxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxxy[i] * fe_0 + tk_yyzz_xxxy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxxy[i] * fz_0;

        tk_yyzzz_xxxz[i] = -2.0 * ts_zzz_xxxz[i] * fbe_0 * fz_0 + tk_zzz_xxxz[i] * fe_0 + tk_yzzz_xxxz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxz[i] * fz_0;

        tk_yyzzz_xxyy[i] =
            -4.0 * ts_yyz_xxyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxyy[i] * fe_0 + tk_yyzz_xxyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxyy[i] * fz_0;

        tk_yyzzz_xxyz[i] = -2.0 * ts_zzz_xxyz[i] * fbe_0 * fz_0 + tk_zzz_xxyz[i] * fe_0 + tk_yzzz_xxz[i] * fe_0 + tk_yzzz_xxyz[i] * pa_y[i] +
                           2.0 * ts_yyzzz_xxyz[i] * fz_0;

        tk_yyzzz_xxzz[i] = -2.0 * ts_zzz_xxzz[i] * fbe_0 * fz_0 + tk_zzz_xxzz[i] * fe_0 + tk_yzzz_xxzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxzz[i] * fz_0;

        tk_yyzzz_xyyy[i] =
            -4.0 * ts_yyz_xyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xyyy[i] * fe_0 + tk_yyzz_xyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xyyy[i] * fz_0;

        tk_yyzzz_xyyz[i] = -2.0 * ts_zzz_xyyz[i] * fbe_0 * fz_0 + tk_zzz_xyyz[i] * fe_0 + 2.0 * tk_yzzz_xyz[i] * fe_0 + tk_yzzz_xyyz[i] * pa_y[i] +
                           2.0 * ts_yyzzz_xyyz[i] * fz_0;

        tk_yyzzz_xyzz[i] = -2.0 * ts_zzz_xyzz[i] * fbe_0 * fz_0 + tk_zzz_xyzz[i] * fe_0 + tk_yzzz_xzz[i] * fe_0 + tk_yzzz_xyzz[i] * pa_y[i] +
                           2.0 * ts_yyzzz_xyzz[i] * fz_0;

        tk_yyzzz_xzzz[i] = -2.0 * ts_zzz_xzzz[i] * fbe_0 * fz_0 + tk_zzz_xzzz[i] * fe_0 + tk_yzzz_xzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xzzz[i] * fz_0;

        tk_yyzzz_yyyy[i] =
            -4.0 * ts_yyz_yyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_yyyy[i] * fe_0 + tk_yyzz_yyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_yyyy[i] * fz_0;

        tk_yyzzz_yyyz[i] = -2.0 * ts_zzz_yyyz[i] * fbe_0 * fz_0 + tk_zzz_yyyz[i] * fe_0 + 3.0 * tk_yzzz_yyz[i] * fe_0 + tk_yzzz_yyyz[i] * pa_y[i] +
                           2.0 * ts_yyzzz_yyyz[i] * fz_0;

        tk_yyzzz_yyzz[i] = -2.0 * ts_zzz_yyzz[i] * fbe_0 * fz_0 + tk_zzz_yyzz[i] * fe_0 + 2.0 * tk_yzzz_yzz[i] * fe_0 + tk_yzzz_yyzz[i] * pa_y[i] +
                           2.0 * ts_yyzzz_yyzz[i] * fz_0;

        tk_yyzzz_yzzz[i] = -2.0 * ts_zzz_yzzz[i] * fbe_0 * fz_0 + tk_zzz_yzzz[i] * fe_0 + tk_yzzz_zzz[i] * fe_0 + tk_yzzz_yzzz[i] * pa_y[i] +
                           2.0 * ts_yyzzz_yzzz[i] * fz_0;

        tk_yyzzz_zzzz[i] = -2.0 * ts_zzz_zzzz[i] * fbe_0 * fz_0 + tk_zzz_zzzz[i] * fe_0 + tk_yzzz_zzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_zzzz[i] * fz_0;
    }

    // Set up 285-300 components of targeted buffer : HG

    auto tk_yzzzz_xxxx = pbuffer.data(idx_kin_hg + 285);

    auto tk_yzzzz_xxxy = pbuffer.data(idx_kin_hg + 286);

    auto tk_yzzzz_xxxz = pbuffer.data(idx_kin_hg + 287);

    auto tk_yzzzz_xxyy = pbuffer.data(idx_kin_hg + 288);

    auto tk_yzzzz_xxyz = pbuffer.data(idx_kin_hg + 289);

    auto tk_yzzzz_xxzz = pbuffer.data(idx_kin_hg + 290);

    auto tk_yzzzz_xyyy = pbuffer.data(idx_kin_hg + 291);

    auto tk_yzzzz_xyyz = pbuffer.data(idx_kin_hg + 292);

    auto tk_yzzzz_xyzz = pbuffer.data(idx_kin_hg + 293);

    auto tk_yzzzz_xzzz = pbuffer.data(idx_kin_hg + 294);

    auto tk_yzzzz_yyyy = pbuffer.data(idx_kin_hg + 295);

    auto tk_yzzzz_yyyz = pbuffer.data(idx_kin_hg + 296);

    auto tk_yzzzz_yyzz = pbuffer.data(idx_kin_hg + 297);

    auto tk_yzzzz_yzzz = pbuffer.data(idx_kin_hg + 298);

    auto tk_yzzzz_zzzz = pbuffer.data(idx_kin_hg + 299);

#pragma omp simd aligned(pa_y,              \
                             tk_yzzzz_xxxx, \
                             tk_yzzzz_xxxy, \
                             tk_yzzzz_xxxz, \
                             tk_yzzzz_xxyy, \
                             tk_yzzzz_xxyz, \
                             tk_yzzzz_xxzz, \
                             tk_yzzzz_xyyy, \
                             tk_yzzzz_xyyz, \
                             tk_yzzzz_xyzz, \
                             tk_yzzzz_xzzz, \
                             tk_yzzzz_yyyy, \
                             tk_yzzzz_yyyz, \
                             tk_yzzzz_yyzz, \
                             tk_yzzzz_yzzz, \
                             tk_yzzzz_zzzz, \
                             tk_zzzz_xxx,   \
                             tk_zzzz_xxxx,  \
                             tk_zzzz_xxxy,  \
                             tk_zzzz_xxxz,  \
                             tk_zzzz_xxy,   \
                             tk_zzzz_xxyy,  \
                             tk_zzzz_xxyz,  \
                             tk_zzzz_xxz,   \
                             tk_zzzz_xxzz,  \
                             tk_zzzz_xyy,   \
                             tk_zzzz_xyyy,  \
                             tk_zzzz_xyyz,  \
                             tk_zzzz_xyz,   \
                             tk_zzzz_xyzz,  \
                             tk_zzzz_xzz,   \
                             tk_zzzz_xzzz,  \
                             tk_zzzz_yyy,   \
                             tk_zzzz_yyyy,  \
                             tk_zzzz_yyyz,  \
                             tk_zzzz_yyz,   \
                             tk_zzzz_yyzz,  \
                             tk_zzzz_yzz,   \
                             tk_zzzz_yzzz,  \
                             tk_zzzz_zzz,   \
                             tk_zzzz_zzzz,  \
                             ts_yzzzz_xxxx, \
                             ts_yzzzz_xxxy, \
                             ts_yzzzz_xxxz, \
                             ts_yzzzz_xxyy, \
                             ts_yzzzz_xxyz, \
                             ts_yzzzz_xxzz, \
                             ts_yzzzz_xyyy, \
                             ts_yzzzz_xyyz, \
                             ts_yzzzz_xyzz, \
                             ts_yzzzz_xzzz, \
                             ts_yzzzz_yyyy, \
                             ts_yzzzz_yyyz, \
                             ts_yzzzz_yyzz, \
                             ts_yzzzz_yzzz, \
                             ts_yzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzz_xxxx[i] = tk_zzzz_xxxx[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxx[i] * fz_0;

        tk_yzzzz_xxxy[i] = tk_zzzz_xxx[i] * fe_0 + tk_zzzz_xxxy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxy[i] * fz_0;

        tk_yzzzz_xxxz[i] = tk_zzzz_xxxz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxz[i] * fz_0;

        tk_yzzzz_xxyy[i] = 2.0 * tk_zzzz_xxy[i] * fe_0 + tk_zzzz_xxyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyy[i] * fz_0;

        tk_yzzzz_xxyz[i] = tk_zzzz_xxz[i] * fe_0 + tk_zzzz_xxyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyz[i] * fz_0;

        tk_yzzzz_xxzz[i] = tk_zzzz_xxzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxzz[i] * fz_0;

        tk_yzzzz_xyyy[i] = 3.0 * tk_zzzz_xyy[i] * fe_0 + tk_zzzz_xyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyy[i] * fz_0;

        tk_yzzzz_xyyz[i] = 2.0 * tk_zzzz_xyz[i] * fe_0 + tk_zzzz_xyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyz[i] * fz_0;

        tk_yzzzz_xyzz[i] = tk_zzzz_xzz[i] * fe_0 + tk_zzzz_xyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyzz[i] * fz_0;

        tk_yzzzz_xzzz[i] = tk_zzzz_xzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xzzz[i] * fz_0;

        tk_yzzzz_yyyy[i] = 4.0 * tk_zzzz_yyy[i] * fe_0 + tk_zzzz_yyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyy[i] * fz_0;

        tk_yzzzz_yyyz[i] = 3.0 * tk_zzzz_yyz[i] * fe_0 + tk_zzzz_yyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyz[i] * fz_0;

        tk_yzzzz_yyzz[i] = 2.0 * tk_zzzz_yzz[i] * fe_0 + tk_zzzz_yyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyzz[i] * fz_0;

        tk_yzzzz_yzzz[i] = tk_zzzz_zzz[i] * fe_0 + tk_zzzz_yzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yzzz[i] * fz_0;

        tk_yzzzz_zzzz[i] = tk_zzzz_zzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_zzzz[i] * fz_0;
    }

    // Set up 300-315 components of targeted buffer : HG

    auto tk_zzzzz_xxxx = pbuffer.data(idx_kin_hg + 300);

    auto tk_zzzzz_xxxy = pbuffer.data(idx_kin_hg + 301);

    auto tk_zzzzz_xxxz = pbuffer.data(idx_kin_hg + 302);

    auto tk_zzzzz_xxyy = pbuffer.data(idx_kin_hg + 303);

    auto tk_zzzzz_xxyz = pbuffer.data(idx_kin_hg + 304);

    auto tk_zzzzz_xxzz = pbuffer.data(idx_kin_hg + 305);

    auto tk_zzzzz_xyyy = pbuffer.data(idx_kin_hg + 306);

    auto tk_zzzzz_xyyz = pbuffer.data(idx_kin_hg + 307);

    auto tk_zzzzz_xyzz = pbuffer.data(idx_kin_hg + 308);

    auto tk_zzzzz_xzzz = pbuffer.data(idx_kin_hg + 309);

    auto tk_zzzzz_yyyy = pbuffer.data(idx_kin_hg + 310);

    auto tk_zzzzz_yyyz = pbuffer.data(idx_kin_hg + 311);

    auto tk_zzzzz_yyzz = pbuffer.data(idx_kin_hg + 312);

    auto tk_zzzzz_yzzz = pbuffer.data(idx_kin_hg + 313);

    auto tk_zzzzz_zzzz = pbuffer.data(idx_kin_hg + 314);

#pragma omp simd aligned(pa_z,              \
                             tk_zzz_xxxx,   \
                             tk_zzz_xxxy,   \
                             tk_zzz_xxxz,   \
                             tk_zzz_xxyy,   \
                             tk_zzz_xxyz,   \
                             tk_zzz_xxzz,   \
                             tk_zzz_xyyy,   \
                             tk_zzz_xyyz,   \
                             tk_zzz_xyzz,   \
                             tk_zzz_xzzz,   \
                             tk_zzz_yyyy,   \
                             tk_zzz_yyyz,   \
                             tk_zzz_yyzz,   \
                             tk_zzz_yzzz,   \
                             tk_zzz_zzzz,   \
                             tk_zzzz_xxx,   \
                             tk_zzzz_xxxx,  \
                             tk_zzzz_xxxy,  \
                             tk_zzzz_xxxz,  \
                             tk_zzzz_xxy,   \
                             tk_zzzz_xxyy,  \
                             tk_zzzz_xxyz,  \
                             tk_zzzz_xxz,   \
                             tk_zzzz_xxzz,  \
                             tk_zzzz_xyy,   \
                             tk_zzzz_xyyy,  \
                             tk_zzzz_xyyz,  \
                             tk_zzzz_xyz,   \
                             tk_zzzz_xyzz,  \
                             tk_zzzz_xzz,   \
                             tk_zzzz_xzzz,  \
                             tk_zzzz_yyy,   \
                             tk_zzzz_yyyy,  \
                             tk_zzzz_yyyz,  \
                             tk_zzzz_yyz,   \
                             tk_zzzz_yyzz,  \
                             tk_zzzz_yzz,   \
                             tk_zzzz_yzzz,  \
                             tk_zzzz_zzz,   \
                             tk_zzzz_zzzz,  \
                             tk_zzzzz_xxxx, \
                             tk_zzzzz_xxxy, \
                             tk_zzzzz_xxxz, \
                             tk_zzzzz_xxyy, \
                             tk_zzzzz_xxyz, \
                             tk_zzzzz_xxzz, \
                             tk_zzzzz_xyyy, \
                             tk_zzzzz_xyyz, \
                             tk_zzzzz_xyzz, \
                             tk_zzzzz_xzzz, \
                             tk_zzzzz_yyyy, \
                             tk_zzzzz_yyyz, \
                             tk_zzzzz_yyzz, \
                             tk_zzzzz_yzzz, \
                             tk_zzzzz_zzzz, \
                             ts_zzz_xxxx,   \
                             ts_zzz_xxxy,   \
                             ts_zzz_xxxz,   \
                             ts_zzz_xxyy,   \
                             ts_zzz_xxyz,   \
                             ts_zzz_xxzz,   \
                             ts_zzz_xyyy,   \
                             ts_zzz_xyyz,   \
                             ts_zzz_xyzz,   \
                             ts_zzz_xzzz,   \
                             ts_zzz_yyyy,   \
                             ts_zzz_yyyz,   \
                             ts_zzz_yyzz,   \
                             ts_zzz_yzzz,   \
                             ts_zzz_zzzz,   \
                             ts_zzzzz_xxxx, \
                             ts_zzzzz_xxxy, \
                             ts_zzzzz_xxxz, \
                             ts_zzzzz_xxyy, \
                             ts_zzzzz_xxyz, \
                             ts_zzzzz_xxzz, \
                             ts_zzzzz_xyyy, \
                             ts_zzzzz_xyyz, \
                             ts_zzzzz_xyzz, \
                             ts_zzzzz_xzzz, \
                             ts_zzzzz_yyyy, \
                             ts_zzzzz_yyyz, \
                             ts_zzzzz_yyzz, \
                             ts_zzzzz_yzzz, \
                             ts_zzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzz_xxxx[i] =
            -8.0 * ts_zzz_xxxx[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxx[i] * fe_0 + tk_zzzz_xxxx[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxx[i] * fz_0;

        tk_zzzzz_xxxy[i] =
            -8.0 * ts_zzz_xxxy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxy[i] * fe_0 + tk_zzzz_xxxy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxy[i] * fz_0;

        tk_zzzzz_xxxz[i] = -8.0 * ts_zzz_xxxz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxz[i] * fe_0 + tk_zzzz_xxx[i] * fe_0 + tk_zzzz_xxxz[i] * pa_z[i] +
                           2.0 * ts_zzzzz_xxxz[i] * fz_0;

        tk_zzzzz_xxyy[i] =
            -8.0 * ts_zzz_xxyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyy[i] * fe_0 + tk_zzzz_xxyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyy[i] * fz_0;

        tk_zzzzz_xxyz[i] = -8.0 * ts_zzz_xxyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyz[i] * fe_0 + tk_zzzz_xxy[i] * fe_0 + tk_zzzz_xxyz[i] * pa_z[i] +
                           2.0 * ts_zzzzz_xxyz[i] * fz_0;

        tk_zzzzz_xxzz[i] = -8.0 * ts_zzz_xxzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxzz[i] * fe_0 + 2.0 * tk_zzzz_xxz[i] * fe_0 +
                           tk_zzzz_xxzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxzz[i] * fz_0;

        tk_zzzzz_xyyy[i] =
            -8.0 * ts_zzz_xyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyy[i] * fe_0 + tk_zzzz_xyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyy[i] * fz_0;

        tk_zzzzz_xyyz[i] = -8.0 * ts_zzz_xyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyz[i] * fe_0 + tk_zzzz_xyy[i] * fe_0 + tk_zzzz_xyyz[i] * pa_z[i] +
                           2.0 * ts_zzzzz_xyyz[i] * fz_0;

        tk_zzzzz_xyzz[i] = -8.0 * ts_zzz_xyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyzz[i] * fe_0 + 2.0 * tk_zzzz_xyz[i] * fe_0 +
                           tk_zzzz_xyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyzz[i] * fz_0;

        tk_zzzzz_xzzz[i] = -8.0 * ts_zzz_xzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xzzz[i] * fe_0 + 3.0 * tk_zzzz_xzz[i] * fe_0 +
                           tk_zzzz_xzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xzzz[i] * fz_0;

        tk_zzzzz_yyyy[i] =
            -8.0 * ts_zzz_yyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyy[i] * fe_0 + tk_zzzz_yyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyy[i] * fz_0;

        tk_zzzzz_yyyz[i] = -8.0 * ts_zzz_yyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyz[i] * fe_0 + tk_zzzz_yyy[i] * fe_0 + tk_zzzz_yyyz[i] * pa_z[i] +
                           2.0 * ts_zzzzz_yyyz[i] * fz_0;

        tk_zzzzz_yyzz[i] = -8.0 * ts_zzz_yyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyzz[i] * fe_0 + 2.0 * tk_zzzz_yyz[i] * fe_0 +
                           tk_zzzz_yyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyzz[i] * fz_0;

        tk_zzzzz_yzzz[i] = -8.0 * ts_zzz_yzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yzzz[i] * fe_0 + 3.0 * tk_zzzz_yzz[i] * fe_0 +
                           tk_zzzz_yzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yzzz[i] * fz_0;

        tk_zzzzz_zzzz[i] = -8.0 * ts_zzz_zzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_zzzz[i] * fe_0 + 4.0 * tk_zzzz_zzz[i] * fe_0 +
                           tk_zzzz_zzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_zzzz[i] * fz_0;
    }
}

}  // namespace kinrec
