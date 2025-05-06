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

#include "KineticEnergyPrimRecGH.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_gh(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_gh,
                            const size_t              idx_ovl_dh,
                            const size_t              idx_kin_dh,
                            const size_t              idx_kin_fg,
                            const size_t              idx_kin_fh,
                            const size_t              idx_ovl_gh,
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

    // Set up components of auxiliary buffer : DH

    auto ts_xx_xxxxx = pbuffer.data(idx_ovl_dh);

    auto ts_xx_xxxxy = pbuffer.data(idx_ovl_dh + 1);

    auto ts_xx_xxxxz = pbuffer.data(idx_ovl_dh + 2);

    auto ts_xx_xxxyy = pbuffer.data(idx_ovl_dh + 3);

    auto ts_xx_xxxyz = pbuffer.data(idx_ovl_dh + 4);

    auto ts_xx_xxxzz = pbuffer.data(idx_ovl_dh + 5);

    auto ts_xx_xxyyy = pbuffer.data(idx_ovl_dh + 6);

    auto ts_xx_xxyyz = pbuffer.data(idx_ovl_dh + 7);

    auto ts_xx_xxyzz = pbuffer.data(idx_ovl_dh + 8);

    auto ts_xx_xxzzz = pbuffer.data(idx_ovl_dh + 9);

    auto ts_xx_xyyyy = pbuffer.data(idx_ovl_dh + 10);

    auto ts_xx_xyyyz = pbuffer.data(idx_ovl_dh + 11);

    auto ts_xx_xyyzz = pbuffer.data(idx_ovl_dh + 12);

    auto ts_xx_xyzzz = pbuffer.data(idx_ovl_dh + 13);

    auto ts_xx_xzzzz = pbuffer.data(idx_ovl_dh + 14);

    auto ts_xx_yyyyy = pbuffer.data(idx_ovl_dh + 15);

    auto ts_xx_yyyyz = pbuffer.data(idx_ovl_dh + 16);

    auto ts_xx_yyyzz = pbuffer.data(idx_ovl_dh + 17);

    auto ts_xx_yyzzz = pbuffer.data(idx_ovl_dh + 18);

    auto ts_xx_yzzzz = pbuffer.data(idx_ovl_dh + 19);

    auto ts_xx_zzzzz = pbuffer.data(idx_ovl_dh + 20);

    auto ts_yy_xxxxx = pbuffer.data(idx_ovl_dh + 63);

    auto ts_yy_xxxxy = pbuffer.data(idx_ovl_dh + 64);

    auto ts_yy_xxxxz = pbuffer.data(idx_ovl_dh + 65);

    auto ts_yy_xxxyy = pbuffer.data(idx_ovl_dh + 66);

    auto ts_yy_xxxyz = pbuffer.data(idx_ovl_dh + 67);

    auto ts_yy_xxxzz = pbuffer.data(idx_ovl_dh + 68);

    auto ts_yy_xxyyy = pbuffer.data(idx_ovl_dh + 69);

    auto ts_yy_xxyyz = pbuffer.data(idx_ovl_dh + 70);

    auto ts_yy_xxyzz = pbuffer.data(idx_ovl_dh + 71);

    auto ts_yy_xxzzz = pbuffer.data(idx_ovl_dh + 72);

    auto ts_yy_xyyyy = pbuffer.data(idx_ovl_dh + 73);

    auto ts_yy_xyyyz = pbuffer.data(idx_ovl_dh + 74);

    auto ts_yy_xyyzz = pbuffer.data(idx_ovl_dh + 75);

    auto ts_yy_xyzzz = pbuffer.data(idx_ovl_dh + 76);

    auto ts_yy_xzzzz = pbuffer.data(idx_ovl_dh + 77);

    auto ts_yy_yyyyy = pbuffer.data(idx_ovl_dh + 78);

    auto ts_yy_yyyyz = pbuffer.data(idx_ovl_dh + 79);

    auto ts_yy_yyyzz = pbuffer.data(idx_ovl_dh + 80);

    auto ts_yy_yyzzz = pbuffer.data(idx_ovl_dh + 81);

    auto ts_yy_yzzzz = pbuffer.data(idx_ovl_dh + 82);

    auto ts_yy_zzzzz = pbuffer.data(idx_ovl_dh + 83);

    auto ts_zz_xxxxx = pbuffer.data(idx_ovl_dh + 105);

    auto ts_zz_xxxxy = pbuffer.data(idx_ovl_dh + 106);

    auto ts_zz_xxxxz = pbuffer.data(idx_ovl_dh + 107);

    auto ts_zz_xxxyy = pbuffer.data(idx_ovl_dh + 108);

    auto ts_zz_xxxyz = pbuffer.data(idx_ovl_dh + 109);

    auto ts_zz_xxxzz = pbuffer.data(idx_ovl_dh + 110);

    auto ts_zz_xxyyy = pbuffer.data(idx_ovl_dh + 111);

    auto ts_zz_xxyyz = pbuffer.data(idx_ovl_dh + 112);

    auto ts_zz_xxyzz = pbuffer.data(idx_ovl_dh + 113);

    auto ts_zz_xxzzz = pbuffer.data(idx_ovl_dh + 114);

    auto ts_zz_xyyyy = pbuffer.data(idx_ovl_dh + 115);

    auto ts_zz_xyyyz = pbuffer.data(idx_ovl_dh + 116);

    auto ts_zz_xyyzz = pbuffer.data(idx_ovl_dh + 117);

    auto ts_zz_xyzzz = pbuffer.data(idx_ovl_dh + 118);

    auto ts_zz_xzzzz = pbuffer.data(idx_ovl_dh + 119);

    auto ts_zz_yyyyy = pbuffer.data(idx_ovl_dh + 120);

    auto ts_zz_yyyyz = pbuffer.data(idx_ovl_dh + 121);

    auto ts_zz_yyyzz = pbuffer.data(idx_ovl_dh + 122);

    auto ts_zz_yyzzz = pbuffer.data(idx_ovl_dh + 123);

    auto ts_zz_yzzzz = pbuffer.data(idx_ovl_dh + 124);

    auto ts_zz_zzzzz = pbuffer.data(idx_ovl_dh + 125);

    // Set up components of auxiliary buffer : DH

    auto tk_xx_xxxxx = pbuffer.data(idx_kin_dh);

    auto tk_xx_xxxxy = pbuffer.data(idx_kin_dh + 1);

    auto tk_xx_xxxxz = pbuffer.data(idx_kin_dh + 2);

    auto tk_xx_xxxyy = pbuffer.data(idx_kin_dh + 3);

    auto tk_xx_xxxyz = pbuffer.data(idx_kin_dh + 4);

    auto tk_xx_xxxzz = pbuffer.data(idx_kin_dh + 5);

    auto tk_xx_xxyyy = pbuffer.data(idx_kin_dh + 6);

    auto tk_xx_xxyyz = pbuffer.data(idx_kin_dh + 7);

    auto tk_xx_xxyzz = pbuffer.data(idx_kin_dh + 8);

    auto tk_xx_xxzzz = pbuffer.data(idx_kin_dh + 9);

    auto tk_xx_xyyyy = pbuffer.data(idx_kin_dh + 10);

    auto tk_xx_xyyyz = pbuffer.data(idx_kin_dh + 11);

    auto tk_xx_xyyzz = pbuffer.data(idx_kin_dh + 12);

    auto tk_xx_xyzzz = pbuffer.data(idx_kin_dh + 13);

    auto tk_xx_xzzzz = pbuffer.data(idx_kin_dh + 14);

    auto tk_xx_yyyyy = pbuffer.data(idx_kin_dh + 15);

    auto tk_xx_yyyyz = pbuffer.data(idx_kin_dh + 16);

    auto tk_xx_yyyzz = pbuffer.data(idx_kin_dh + 17);

    auto tk_xx_yyzzz = pbuffer.data(idx_kin_dh + 18);

    auto tk_xx_yzzzz = pbuffer.data(idx_kin_dh + 19);

    auto tk_xx_zzzzz = pbuffer.data(idx_kin_dh + 20);

    auto tk_yy_xxxxx = pbuffer.data(idx_kin_dh + 63);

    auto tk_yy_xxxxy = pbuffer.data(idx_kin_dh + 64);

    auto tk_yy_xxxxz = pbuffer.data(idx_kin_dh + 65);

    auto tk_yy_xxxyy = pbuffer.data(idx_kin_dh + 66);

    auto tk_yy_xxxyz = pbuffer.data(idx_kin_dh + 67);

    auto tk_yy_xxxzz = pbuffer.data(idx_kin_dh + 68);

    auto tk_yy_xxyyy = pbuffer.data(idx_kin_dh + 69);

    auto tk_yy_xxyyz = pbuffer.data(idx_kin_dh + 70);

    auto tk_yy_xxyzz = pbuffer.data(idx_kin_dh + 71);

    auto tk_yy_xxzzz = pbuffer.data(idx_kin_dh + 72);

    auto tk_yy_xyyyy = pbuffer.data(idx_kin_dh + 73);

    auto tk_yy_xyyyz = pbuffer.data(idx_kin_dh + 74);

    auto tk_yy_xyyzz = pbuffer.data(idx_kin_dh + 75);

    auto tk_yy_xyzzz = pbuffer.data(idx_kin_dh + 76);

    auto tk_yy_xzzzz = pbuffer.data(idx_kin_dh + 77);

    auto tk_yy_yyyyy = pbuffer.data(idx_kin_dh + 78);

    auto tk_yy_yyyyz = pbuffer.data(idx_kin_dh + 79);

    auto tk_yy_yyyzz = pbuffer.data(idx_kin_dh + 80);

    auto tk_yy_yyzzz = pbuffer.data(idx_kin_dh + 81);

    auto tk_yy_yzzzz = pbuffer.data(idx_kin_dh + 82);

    auto tk_yy_zzzzz = pbuffer.data(idx_kin_dh + 83);

    auto tk_zz_xxxxx = pbuffer.data(idx_kin_dh + 105);

    auto tk_zz_xxxxy = pbuffer.data(idx_kin_dh + 106);

    auto tk_zz_xxxxz = pbuffer.data(idx_kin_dh + 107);

    auto tk_zz_xxxyy = pbuffer.data(idx_kin_dh + 108);

    auto tk_zz_xxxyz = pbuffer.data(idx_kin_dh + 109);

    auto tk_zz_xxxzz = pbuffer.data(idx_kin_dh + 110);

    auto tk_zz_xxyyy = pbuffer.data(idx_kin_dh + 111);

    auto tk_zz_xxyyz = pbuffer.data(idx_kin_dh + 112);

    auto tk_zz_xxyzz = pbuffer.data(idx_kin_dh + 113);

    auto tk_zz_xxzzz = pbuffer.data(idx_kin_dh + 114);

    auto tk_zz_xyyyy = pbuffer.data(idx_kin_dh + 115);

    auto tk_zz_xyyyz = pbuffer.data(idx_kin_dh + 116);

    auto tk_zz_xyyzz = pbuffer.data(idx_kin_dh + 117);

    auto tk_zz_xyzzz = pbuffer.data(idx_kin_dh + 118);

    auto tk_zz_xzzzz = pbuffer.data(idx_kin_dh + 119);

    auto tk_zz_yyyyy = pbuffer.data(idx_kin_dh + 120);

    auto tk_zz_yyyyz = pbuffer.data(idx_kin_dh + 121);

    auto tk_zz_yyyzz = pbuffer.data(idx_kin_dh + 122);

    auto tk_zz_yyzzz = pbuffer.data(idx_kin_dh + 123);

    auto tk_zz_yzzzz = pbuffer.data(idx_kin_dh + 124);

    auto tk_zz_zzzzz = pbuffer.data(idx_kin_dh + 125);

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

    auto tk_xxz_xxxz = pbuffer.data(idx_kin_fg + 32);

    auto tk_xxz_xxyz = pbuffer.data(idx_kin_fg + 34);

    auto tk_xxz_xxzz = pbuffer.data(idx_kin_fg + 35);

    auto tk_xxz_xyyz = pbuffer.data(idx_kin_fg + 37);

    auto tk_xxz_xyzz = pbuffer.data(idx_kin_fg + 38);

    auto tk_xxz_xzzz = pbuffer.data(idx_kin_fg + 39);

    auto tk_xxz_yyyz = pbuffer.data(idx_kin_fg + 41);

    auto tk_xxz_yyzz = pbuffer.data(idx_kin_fg + 42);

    auto tk_xxz_yzzz = pbuffer.data(idx_kin_fg + 43);

    auto tk_xxz_zzzz = pbuffer.data(idx_kin_fg + 44);

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

    auto tk_xzz_xxxz = pbuffer.data(idx_kin_fg + 77);

    auto tk_xzz_xxyz = pbuffer.data(idx_kin_fg + 79);

    auto tk_xzz_xxzz = pbuffer.data(idx_kin_fg + 80);

    auto tk_xzz_xyyz = pbuffer.data(idx_kin_fg + 82);

    auto tk_xzz_xyzz = pbuffer.data(idx_kin_fg + 83);

    auto tk_xzz_xzzz = pbuffer.data(idx_kin_fg + 84);

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

    auto tk_yyz_xxxz = pbuffer.data(idx_kin_fg + 107);

    auto tk_yyz_xxyz = pbuffer.data(idx_kin_fg + 109);

    auto tk_yyz_xxzz = pbuffer.data(idx_kin_fg + 110);

    auto tk_yyz_xyyz = pbuffer.data(idx_kin_fg + 112);

    auto tk_yyz_xyzz = pbuffer.data(idx_kin_fg + 113);

    auto tk_yyz_xzzz = pbuffer.data(idx_kin_fg + 114);

    auto tk_yyz_yyyz = pbuffer.data(idx_kin_fg + 116);

    auto tk_yyz_yyzz = pbuffer.data(idx_kin_fg + 117);

    auto tk_yyz_yzzz = pbuffer.data(idx_kin_fg + 118);

    auto tk_yyz_zzzz = pbuffer.data(idx_kin_fg + 119);

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

    // Set up components of auxiliary buffer : FH

    auto tk_xxx_xxxxx = pbuffer.data(idx_kin_fh);

    auto tk_xxx_xxxxy = pbuffer.data(idx_kin_fh + 1);

    auto tk_xxx_xxxxz = pbuffer.data(idx_kin_fh + 2);

    auto tk_xxx_xxxyy = pbuffer.data(idx_kin_fh + 3);

    auto tk_xxx_xxxyz = pbuffer.data(idx_kin_fh + 4);

    auto tk_xxx_xxxzz = pbuffer.data(idx_kin_fh + 5);

    auto tk_xxx_xxyyy = pbuffer.data(idx_kin_fh + 6);

    auto tk_xxx_xxyyz = pbuffer.data(idx_kin_fh + 7);

    auto tk_xxx_xxyzz = pbuffer.data(idx_kin_fh + 8);

    auto tk_xxx_xxzzz = pbuffer.data(idx_kin_fh + 9);

    auto tk_xxx_xyyyy = pbuffer.data(idx_kin_fh + 10);

    auto tk_xxx_xyyyz = pbuffer.data(idx_kin_fh + 11);

    auto tk_xxx_xyyzz = pbuffer.data(idx_kin_fh + 12);

    auto tk_xxx_xyzzz = pbuffer.data(idx_kin_fh + 13);

    auto tk_xxx_xzzzz = pbuffer.data(idx_kin_fh + 14);

    auto tk_xxx_yyyyy = pbuffer.data(idx_kin_fh + 15);

    auto tk_xxx_yyyyz = pbuffer.data(idx_kin_fh + 16);

    auto tk_xxx_yyyzz = pbuffer.data(idx_kin_fh + 17);

    auto tk_xxx_yyzzz = pbuffer.data(idx_kin_fh + 18);

    auto tk_xxx_yzzzz = pbuffer.data(idx_kin_fh + 19);

    auto tk_xxx_zzzzz = pbuffer.data(idx_kin_fh + 20);

    auto tk_xxy_xxxxx = pbuffer.data(idx_kin_fh + 21);

    auto tk_xxy_xxxxy = pbuffer.data(idx_kin_fh + 22);

    auto tk_xxy_xxxxz = pbuffer.data(idx_kin_fh + 23);

    auto tk_xxy_xxxyy = pbuffer.data(idx_kin_fh + 24);

    auto tk_xxy_xxxzz = pbuffer.data(idx_kin_fh + 26);

    auto tk_xxy_xxyyy = pbuffer.data(idx_kin_fh + 27);

    auto tk_xxy_xxzzz = pbuffer.data(idx_kin_fh + 30);

    auto tk_xxy_xyyyy = pbuffer.data(idx_kin_fh + 31);

    auto tk_xxy_xzzzz = pbuffer.data(idx_kin_fh + 35);

    auto tk_xxy_yyyyy = pbuffer.data(idx_kin_fh + 36);

    auto tk_xxz_xxxxx = pbuffer.data(idx_kin_fh + 42);

    auto tk_xxz_xxxxy = pbuffer.data(idx_kin_fh + 43);

    auto tk_xxz_xxxxz = pbuffer.data(idx_kin_fh + 44);

    auto tk_xxz_xxxyy = pbuffer.data(idx_kin_fh + 45);

    auto tk_xxz_xxxyz = pbuffer.data(idx_kin_fh + 46);

    auto tk_xxz_xxxzz = pbuffer.data(idx_kin_fh + 47);

    auto tk_xxz_xxyyy = pbuffer.data(idx_kin_fh + 48);

    auto tk_xxz_xxyyz = pbuffer.data(idx_kin_fh + 49);

    auto tk_xxz_xxyzz = pbuffer.data(idx_kin_fh + 50);

    auto tk_xxz_xxzzz = pbuffer.data(idx_kin_fh + 51);

    auto tk_xxz_xyyyy = pbuffer.data(idx_kin_fh + 52);

    auto tk_xxz_xyyyz = pbuffer.data(idx_kin_fh + 53);

    auto tk_xxz_xyyzz = pbuffer.data(idx_kin_fh + 54);

    auto tk_xxz_xyzzz = pbuffer.data(idx_kin_fh + 55);

    auto tk_xxz_xzzzz = pbuffer.data(idx_kin_fh + 56);

    auto tk_xxz_yyyyz = pbuffer.data(idx_kin_fh + 58);

    auto tk_xxz_yyyzz = pbuffer.data(idx_kin_fh + 59);

    auto tk_xxz_yyzzz = pbuffer.data(idx_kin_fh + 60);

    auto tk_xxz_yzzzz = pbuffer.data(idx_kin_fh + 61);

    auto tk_xxz_zzzzz = pbuffer.data(idx_kin_fh + 62);

    auto tk_xyy_xxxxx = pbuffer.data(idx_kin_fh + 63);

    auto tk_xyy_xxxxy = pbuffer.data(idx_kin_fh + 64);

    auto tk_xyy_xxxyy = pbuffer.data(idx_kin_fh + 66);

    auto tk_xyy_xxxyz = pbuffer.data(idx_kin_fh + 67);

    auto tk_xyy_xxyyy = pbuffer.data(idx_kin_fh + 69);

    auto tk_xyy_xxyyz = pbuffer.data(idx_kin_fh + 70);

    auto tk_xyy_xxyzz = pbuffer.data(idx_kin_fh + 71);

    auto tk_xyy_xyyyy = pbuffer.data(idx_kin_fh + 73);

    auto tk_xyy_xyyyz = pbuffer.data(idx_kin_fh + 74);

    auto tk_xyy_xyyzz = pbuffer.data(idx_kin_fh + 75);

    auto tk_xyy_xyzzz = pbuffer.data(idx_kin_fh + 76);

    auto tk_xyy_yyyyy = pbuffer.data(idx_kin_fh + 78);

    auto tk_xyy_yyyyz = pbuffer.data(idx_kin_fh + 79);

    auto tk_xyy_yyyzz = pbuffer.data(idx_kin_fh + 80);

    auto tk_xyy_yyzzz = pbuffer.data(idx_kin_fh + 81);

    auto tk_xyy_yzzzz = pbuffer.data(idx_kin_fh + 82);

    auto tk_xyy_zzzzz = pbuffer.data(idx_kin_fh + 83);

    auto tk_xzz_xxxxx = pbuffer.data(idx_kin_fh + 105);

    auto tk_xzz_xxxxz = pbuffer.data(idx_kin_fh + 107);

    auto tk_xzz_xxxyz = pbuffer.data(idx_kin_fh + 109);

    auto tk_xzz_xxxzz = pbuffer.data(idx_kin_fh + 110);

    auto tk_xzz_xxyyz = pbuffer.data(idx_kin_fh + 112);

    auto tk_xzz_xxyzz = pbuffer.data(idx_kin_fh + 113);

    auto tk_xzz_xxzzz = pbuffer.data(idx_kin_fh + 114);

    auto tk_xzz_xyyyz = pbuffer.data(idx_kin_fh + 116);

    auto tk_xzz_xyyzz = pbuffer.data(idx_kin_fh + 117);

    auto tk_xzz_xyzzz = pbuffer.data(idx_kin_fh + 118);

    auto tk_xzz_xzzzz = pbuffer.data(idx_kin_fh + 119);

    auto tk_xzz_yyyyy = pbuffer.data(idx_kin_fh + 120);

    auto tk_xzz_yyyyz = pbuffer.data(idx_kin_fh + 121);

    auto tk_xzz_yyyzz = pbuffer.data(idx_kin_fh + 122);

    auto tk_xzz_yyzzz = pbuffer.data(idx_kin_fh + 123);

    auto tk_xzz_yzzzz = pbuffer.data(idx_kin_fh + 124);

    auto tk_xzz_zzzzz = pbuffer.data(idx_kin_fh + 125);

    auto tk_yyy_xxxxx = pbuffer.data(idx_kin_fh + 126);

    auto tk_yyy_xxxxy = pbuffer.data(idx_kin_fh + 127);

    auto tk_yyy_xxxxz = pbuffer.data(idx_kin_fh + 128);

    auto tk_yyy_xxxyy = pbuffer.data(idx_kin_fh + 129);

    auto tk_yyy_xxxyz = pbuffer.data(idx_kin_fh + 130);

    auto tk_yyy_xxxzz = pbuffer.data(idx_kin_fh + 131);

    auto tk_yyy_xxyyy = pbuffer.data(idx_kin_fh + 132);

    auto tk_yyy_xxyyz = pbuffer.data(idx_kin_fh + 133);

    auto tk_yyy_xxyzz = pbuffer.data(idx_kin_fh + 134);

    auto tk_yyy_xxzzz = pbuffer.data(idx_kin_fh + 135);

    auto tk_yyy_xyyyy = pbuffer.data(idx_kin_fh + 136);

    auto tk_yyy_xyyyz = pbuffer.data(idx_kin_fh + 137);

    auto tk_yyy_xyyzz = pbuffer.data(idx_kin_fh + 138);

    auto tk_yyy_xyzzz = pbuffer.data(idx_kin_fh + 139);

    auto tk_yyy_xzzzz = pbuffer.data(idx_kin_fh + 140);

    auto tk_yyy_yyyyy = pbuffer.data(idx_kin_fh + 141);

    auto tk_yyy_yyyyz = pbuffer.data(idx_kin_fh + 142);

    auto tk_yyy_yyyzz = pbuffer.data(idx_kin_fh + 143);

    auto tk_yyy_yyzzz = pbuffer.data(idx_kin_fh + 144);

    auto tk_yyy_yzzzz = pbuffer.data(idx_kin_fh + 145);

    auto tk_yyy_zzzzz = pbuffer.data(idx_kin_fh + 146);

    auto tk_yyz_xxxxy = pbuffer.data(idx_kin_fh + 148);

    auto tk_yyz_xxxxz = pbuffer.data(idx_kin_fh + 149);

    auto tk_yyz_xxxyy = pbuffer.data(idx_kin_fh + 150);

    auto tk_yyz_xxxyz = pbuffer.data(idx_kin_fh + 151);

    auto tk_yyz_xxxzz = pbuffer.data(idx_kin_fh + 152);

    auto tk_yyz_xxyyy = pbuffer.data(idx_kin_fh + 153);

    auto tk_yyz_xxyyz = pbuffer.data(idx_kin_fh + 154);

    auto tk_yyz_xxyzz = pbuffer.data(idx_kin_fh + 155);

    auto tk_yyz_xxzzz = pbuffer.data(idx_kin_fh + 156);

    auto tk_yyz_xyyyy = pbuffer.data(idx_kin_fh + 157);

    auto tk_yyz_xyyyz = pbuffer.data(idx_kin_fh + 158);

    auto tk_yyz_xyyzz = pbuffer.data(idx_kin_fh + 159);

    auto tk_yyz_xyzzz = pbuffer.data(idx_kin_fh + 160);

    auto tk_yyz_xzzzz = pbuffer.data(idx_kin_fh + 161);

    auto tk_yyz_yyyyy = pbuffer.data(idx_kin_fh + 162);

    auto tk_yyz_yyyyz = pbuffer.data(idx_kin_fh + 163);

    auto tk_yyz_yyyzz = pbuffer.data(idx_kin_fh + 164);

    auto tk_yyz_yyzzz = pbuffer.data(idx_kin_fh + 165);

    auto tk_yyz_yzzzz = pbuffer.data(idx_kin_fh + 166);

    auto tk_yyz_zzzzz = pbuffer.data(idx_kin_fh + 167);

    auto tk_yzz_xxxxx = pbuffer.data(idx_kin_fh + 168);

    auto tk_yzz_xxxxy = pbuffer.data(idx_kin_fh + 169);

    auto tk_yzz_xxxxz = pbuffer.data(idx_kin_fh + 170);

    auto tk_yzz_xxxyy = pbuffer.data(idx_kin_fh + 171);

    auto tk_yzz_xxxyz = pbuffer.data(idx_kin_fh + 172);

    auto tk_yzz_xxxzz = pbuffer.data(idx_kin_fh + 173);

    auto tk_yzz_xxyyy = pbuffer.data(idx_kin_fh + 174);

    auto tk_yzz_xxyyz = pbuffer.data(idx_kin_fh + 175);

    auto tk_yzz_xxyzz = pbuffer.data(idx_kin_fh + 176);

    auto tk_yzz_xxzzz = pbuffer.data(idx_kin_fh + 177);

    auto tk_yzz_xyyyy = pbuffer.data(idx_kin_fh + 178);

    auto tk_yzz_xyyyz = pbuffer.data(idx_kin_fh + 179);

    auto tk_yzz_xyyzz = pbuffer.data(idx_kin_fh + 180);

    auto tk_yzz_xyzzz = pbuffer.data(idx_kin_fh + 181);

    auto tk_yzz_xzzzz = pbuffer.data(idx_kin_fh + 182);

    auto tk_yzz_yyyyy = pbuffer.data(idx_kin_fh + 183);

    auto tk_yzz_yyyyz = pbuffer.data(idx_kin_fh + 184);

    auto tk_yzz_yyyzz = pbuffer.data(idx_kin_fh + 185);

    auto tk_yzz_yyzzz = pbuffer.data(idx_kin_fh + 186);

    auto tk_yzz_yzzzz = pbuffer.data(idx_kin_fh + 187);

    auto tk_yzz_zzzzz = pbuffer.data(idx_kin_fh + 188);

    auto tk_zzz_xxxxx = pbuffer.data(idx_kin_fh + 189);

    auto tk_zzz_xxxxy = pbuffer.data(idx_kin_fh + 190);

    auto tk_zzz_xxxxz = pbuffer.data(idx_kin_fh + 191);

    auto tk_zzz_xxxyy = pbuffer.data(idx_kin_fh + 192);

    auto tk_zzz_xxxyz = pbuffer.data(idx_kin_fh + 193);

    auto tk_zzz_xxxzz = pbuffer.data(idx_kin_fh + 194);

    auto tk_zzz_xxyyy = pbuffer.data(idx_kin_fh + 195);

    auto tk_zzz_xxyyz = pbuffer.data(idx_kin_fh + 196);

    auto tk_zzz_xxyzz = pbuffer.data(idx_kin_fh + 197);

    auto tk_zzz_xxzzz = pbuffer.data(idx_kin_fh + 198);

    auto tk_zzz_xyyyy = pbuffer.data(idx_kin_fh + 199);

    auto tk_zzz_xyyyz = pbuffer.data(idx_kin_fh + 200);

    auto tk_zzz_xyyzz = pbuffer.data(idx_kin_fh + 201);

    auto tk_zzz_xyzzz = pbuffer.data(idx_kin_fh + 202);

    auto tk_zzz_xzzzz = pbuffer.data(idx_kin_fh + 203);

    auto tk_zzz_yyyyy = pbuffer.data(idx_kin_fh + 204);

    auto tk_zzz_yyyyz = pbuffer.data(idx_kin_fh + 205);

    auto tk_zzz_yyyzz = pbuffer.data(idx_kin_fh + 206);

    auto tk_zzz_yyzzz = pbuffer.data(idx_kin_fh + 207);

    auto tk_zzz_yzzzz = pbuffer.data(idx_kin_fh + 208);

    auto tk_zzz_zzzzz = pbuffer.data(idx_kin_fh + 209);

    // Set up components of auxiliary buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_ovl_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_ovl_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_ovl_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_ovl_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_ovl_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_ovl_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_ovl_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_ovl_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_ovl_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_ovl_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_ovl_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_ovl_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_ovl_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_ovl_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_ovl_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_ovl_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_ovl_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_ovl_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_ovl_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_ovl_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_ovl_gh + 20);

    auto ts_xxxy_xxxxx = pbuffer.data(idx_ovl_gh + 21);

    auto ts_xxxy_xxxxy = pbuffer.data(idx_ovl_gh + 22);

    auto ts_xxxy_xxxxz = pbuffer.data(idx_ovl_gh + 23);

    auto ts_xxxy_xxxyy = pbuffer.data(idx_ovl_gh + 24);

    auto ts_xxxy_xxxyz = pbuffer.data(idx_ovl_gh + 25);

    auto ts_xxxy_xxxzz = pbuffer.data(idx_ovl_gh + 26);

    auto ts_xxxy_xxyyy = pbuffer.data(idx_ovl_gh + 27);

    auto ts_xxxy_xxyyz = pbuffer.data(idx_ovl_gh + 28);

    auto ts_xxxy_xxyzz = pbuffer.data(idx_ovl_gh + 29);

    auto ts_xxxy_xxzzz = pbuffer.data(idx_ovl_gh + 30);

    auto ts_xxxy_xyyyy = pbuffer.data(idx_ovl_gh + 31);

    auto ts_xxxy_xyyyz = pbuffer.data(idx_ovl_gh + 32);

    auto ts_xxxy_xyyzz = pbuffer.data(idx_ovl_gh + 33);

    auto ts_xxxy_xyzzz = pbuffer.data(idx_ovl_gh + 34);

    auto ts_xxxy_xzzzz = pbuffer.data(idx_ovl_gh + 35);

    auto ts_xxxy_yyyyy = pbuffer.data(idx_ovl_gh + 36);

    auto ts_xxxy_yyyyz = pbuffer.data(idx_ovl_gh + 37);

    auto ts_xxxy_yyyzz = pbuffer.data(idx_ovl_gh + 38);

    auto ts_xxxy_yyzzz = pbuffer.data(idx_ovl_gh + 39);

    auto ts_xxxy_yzzzz = pbuffer.data(idx_ovl_gh + 40);

    auto ts_xxxy_zzzzz = pbuffer.data(idx_ovl_gh + 41);

    auto ts_xxxz_xxxxx = pbuffer.data(idx_ovl_gh + 42);

    auto ts_xxxz_xxxxy = pbuffer.data(idx_ovl_gh + 43);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_ovl_gh + 44);

    auto ts_xxxz_xxxyy = pbuffer.data(idx_ovl_gh + 45);

    auto ts_xxxz_xxxyz = pbuffer.data(idx_ovl_gh + 46);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_ovl_gh + 47);

    auto ts_xxxz_xxyyy = pbuffer.data(idx_ovl_gh + 48);

    auto ts_xxxz_xxyyz = pbuffer.data(idx_ovl_gh + 49);

    auto ts_xxxz_xxyzz = pbuffer.data(idx_ovl_gh + 50);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_ovl_gh + 51);

    auto ts_xxxz_xyyyy = pbuffer.data(idx_ovl_gh + 52);

    auto ts_xxxz_xyyyz = pbuffer.data(idx_ovl_gh + 53);

    auto ts_xxxz_xyyzz = pbuffer.data(idx_ovl_gh + 54);

    auto ts_xxxz_xyzzz = pbuffer.data(idx_ovl_gh + 55);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_ovl_gh + 56);

    auto ts_xxxz_yyyyy = pbuffer.data(idx_ovl_gh + 57);

    auto ts_xxxz_yyyyz = pbuffer.data(idx_ovl_gh + 58);

    auto ts_xxxz_yyyzz = pbuffer.data(idx_ovl_gh + 59);

    auto ts_xxxz_yyzzz = pbuffer.data(idx_ovl_gh + 60);

    auto ts_xxxz_yzzzz = pbuffer.data(idx_ovl_gh + 61);

    auto ts_xxxz_zzzzz = pbuffer.data(idx_ovl_gh + 62);

    auto ts_xxyy_xxxxx = pbuffer.data(idx_ovl_gh + 63);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_ovl_gh + 64);

    auto ts_xxyy_xxxxz = pbuffer.data(idx_ovl_gh + 65);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_ovl_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_ovl_gh + 67);

    auto ts_xxyy_xxxzz = pbuffer.data(idx_ovl_gh + 68);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_ovl_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_ovl_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_ovl_gh + 71);

    auto ts_xxyy_xxzzz = pbuffer.data(idx_ovl_gh + 72);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_ovl_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_ovl_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_ovl_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_ovl_gh + 76);

    auto ts_xxyy_xzzzz = pbuffer.data(idx_ovl_gh + 77);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_ovl_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_ovl_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_ovl_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_ovl_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_ovl_gh + 82);

    auto ts_xxyy_zzzzz = pbuffer.data(idx_ovl_gh + 83);

    auto ts_xxyz_xxxxx = pbuffer.data(idx_ovl_gh + 84);

    auto ts_xxyz_xxxxy = pbuffer.data(idx_ovl_gh + 85);

    auto ts_xxyz_xxxxz = pbuffer.data(idx_ovl_gh + 86);

    auto ts_xxyz_xxxyy = pbuffer.data(idx_ovl_gh + 87);

    auto ts_xxyz_xxxyz = pbuffer.data(idx_ovl_gh + 88);

    auto ts_xxyz_xxxzz = pbuffer.data(idx_ovl_gh + 89);

    auto ts_xxyz_xxyyy = pbuffer.data(idx_ovl_gh + 90);

    auto ts_xxyz_xxyyz = pbuffer.data(idx_ovl_gh + 91);

    auto ts_xxyz_xxyzz = pbuffer.data(idx_ovl_gh + 92);

    auto ts_xxyz_xxzzz = pbuffer.data(idx_ovl_gh + 93);

    auto ts_xxyz_xyyyy = pbuffer.data(idx_ovl_gh + 94);

    auto ts_xxyz_xyyyz = pbuffer.data(idx_ovl_gh + 95);

    auto ts_xxyz_xyyzz = pbuffer.data(idx_ovl_gh + 96);

    auto ts_xxyz_xyzzz = pbuffer.data(idx_ovl_gh + 97);

    auto ts_xxyz_xzzzz = pbuffer.data(idx_ovl_gh + 98);

    auto ts_xxyz_yyyyy = pbuffer.data(idx_ovl_gh + 99);

    auto ts_xxyz_yyyyz = pbuffer.data(idx_ovl_gh + 100);

    auto ts_xxyz_yyyzz = pbuffer.data(idx_ovl_gh + 101);

    auto ts_xxyz_yyzzz = pbuffer.data(idx_ovl_gh + 102);

    auto ts_xxyz_yzzzz = pbuffer.data(idx_ovl_gh + 103);

    auto ts_xxyz_zzzzz = pbuffer.data(idx_ovl_gh + 104);

    auto ts_xxzz_xxxxx = pbuffer.data(idx_ovl_gh + 105);

    auto ts_xxzz_xxxxy = pbuffer.data(idx_ovl_gh + 106);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_ovl_gh + 107);

    auto ts_xxzz_xxxyy = pbuffer.data(idx_ovl_gh + 108);

    auto ts_xxzz_xxxyz = pbuffer.data(idx_ovl_gh + 109);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_ovl_gh + 110);

    auto ts_xxzz_xxyyy = pbuffer.data(idx_ovl_gh + 111);

    auto ts_xxzz_xxyyz = pbuffer.data(idx_ovl_gh + 112);

    auto ts_xxzz_xxyzz = pbuffer.data(idx_ovl_gh + 113);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_ovl_gh + 114);

    auto ts_xxzz_xyyyy = pbuffer.data(idx_ovl_gh + 115);

    auto ts_xxzz_xyyyz = pbuffer.data(idx_ovl_gh + 116);

    auto ts_xxzz_xyyzz = pbuffer.data(idx_ovl_gh + 117);

    auto ts_xxzz_xyzzz = pbuffer.data(idx_ovl_gh + 118);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_ovl_gh + 119);

    auto ts_xxzz_yyyyy = pbuffer.data(idx_ovl_gh + 120);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_ovl_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_ovl_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_ovl_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_ovl_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_ovl_gh + 125);

    auto ts_xyyy_xxxxx = pbuffer.data(idx_ovl_gh + 126);

    auto ts_xyyy_xxxxy = pbuffer.data(idx_ovl_gh + 127);

    auto ts_xyyy_xxxxz = pbuffer.data(idx_ovl_gh + 128);

    auto ts_xyyy_xxxyy = pbuffer.data(idx_ovl_gh + 129);

    auto ts_xyyy_xxxyz = pbuffer.data(idx_ovl_gh + 130);

    auto ts_xyyy_xxxzz = pbuffer.data(idx_ovl_gh + 131);

    auto ts_xyyy_xxyyy = pbuffer.data(idx_ovl_gh + 132);

    auto ts_xyyy_xxyyz = pbuffer.data(idx_ovl_gh + 133);

    auto ts_xyyy_xxyzz = pbuffer.data(idx_ovl_gh + 134);

    auto ts_xyyy_xxzzz = pbuffer.data(idx_ovl_gh + 135);

    auto ts_xyyy_xyyyy = pbuffer.data(idx_ovl_gh + 136);

    auto ts_xyyy_xyyyz = pbuffer.data(idx_ovl_gh + 137);

    auto ts_xyyy_xyyzz = pbuffer.data(idx_ovl_gh + 138);

    auto ts_xyyy_xyzzz = pbuffer.data(idx_ovl_gh + 139);

    auto ts_xyyy_xzzzz = pbuffer.data(idx_ovl_gh + 140);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_ovl_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_ovl_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_ovl_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_ovl_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_ovl_gh + 145);

    auto ts_xyyy_zzzzz = pbuffer.data(idx_ovl_gh + 146);

    auto ts_xyyz_xxxxx = pbuffer.data(idx_ovl_gh + 147);

    auto ts_xyyz_xxxxy = pbuffer.data(idx_ovl_gh + 148);

    auto ts_xyyz_xxxxz = pbuffer.data(idx_ovl_gh + 149);

    auto ts_xyyz_xxxyy = pbuffer.data(idx_ovl_gh + 150);

    auto ts_xyyz_xxxyz = pbuffer.data(idx_ovl_gh + 151);

    auto ts_xyyz_xxxzz = pbuffer.data(idx_ovl_gh + 152);

    auto ts_xyyz_xxyyy = pbuffer.data(idx_ovl_gh + 153);

    auto ts_xyyz_xxyyz = pbuffer.data(idx_ovl_gh + 154);

    auto ts_xyyz_xxyzz = pbuffer.data(idx_ovl_gh + 155);

    auto ts_xyyz_xxzzz = pbuffer.data(idx_ovl_gh + 156);

    auto ts_xyyz_xyyyy = pbuffer.data(idx_ovl_gh + 157);

    auto ts_xyyz_xyyyz = pbuffer.data(idx_ovl_gh + 158);

    auto ts_xyyz_xyyzz = pbuffer.data(idx_ovl_gh + 159);

    auto ts_xyyz_xyzzz = pbuffer.data(idx_ovl_gh + 160);

    auto ts_xyyz_xzzzz = pbuffer.data(idx_ovl_gh + 161);

    auto ts_xyyz_yyyyy = pbuffer.data(idx_ovl_gh + 162);

    auto ts_xyyz_yyyyz = pbuffer.data(idx_ovl_gh + 163);

    auto ts_xyyz_yyyzz = pbuffer.data(idx_ovl_gh + 164);

    auto ts_xyyz_yyzzz = pbuffer.data(idx_ovl_gh + 165);

    auto ts_xyyz_yzzzz = pbuffer.data(idx_ovl_gh + 166);

    auto ts_xyyz_zzzzz = pbuffer.data(idx_ovl_gh + 167);

    auto ts_xyzz_xxxxx = pbuffer.data(idx_ovl_gh + 168);

    auto ts_xyzz_xxxxy = pbuffer.data(idx_ovl_gh + 169);

    auto ts_xyzz_xxxxz = pbuffer.data(idx_ovl_gh + 170);

    auto ts_xyzz_xxxyy = pbuffer.data(idx_ovl_gh + 171);

    auto ts_xyzz_xxxyz = pbuffer.data(idx_ovl_gh + 172);

    auto ts_xyzz_xxxzz = pbuffer.data(idx_ovl_gh + 173);

    auto ts_xyzz_xxyyy = pbuffer.data(idx_ovl_gh + 174);

    auto ts_xyzz_xxyyz = pbuffer.data(idx_ovl_gh + 175);

    auto ts_xyzz_xxyzz = pbuffer.data(idx_ovl_gh + 176);

    auto ts_xyzz_xxzzz = pbuffer.data(idx_ovl_gh + 177);

    auto ts_xyzz_xyyyy = pbuffer.data(idx_ovl_gh + 178);

    auto ts_xyzz_xyyyz = pbuffer.data(idx_ovl_gh + 179);

    auto ts_xyzz_xyyzz = pbuffer.data(idx_ovl_gh + 180);

    auto ts_xyzz_xyzzz = pbuffer.data(idx_ovl_gh + 181);

    auto ts_xyzz_xzzzz = pbuffer.data(idx_ovl_gh + 182);

    auto ts_xyzz_yyyyy = pbuffer.data(idx_ovl_gh + 183);

    auto ts_xyzz_yyyyz = pbuffer.data(idx_ovl_gh + 184);

    auto ts_xyzz_yyyzz = pbuffer.data(idx_ovl_gh + 185);

    auto ts_xyzz_yyzzz = pbuffer.data(idx_ovl_gh + 186);

    auto ts_xyzz_yzzzz = pbuffer.data(idx_ovl_gh + 187);

    auto ts_xyzz_zzzzz = pbuffer.data(idx_ovl_gh + 188);

    auto ts_xzzz_xxxxx = pbuffer.data(idx_ovl_gh + 189);

    auto ts_xzzz_xxxxy = pbuffer.data(idx_ovl_gh + 190);

    auto ts_xzzz_xxxxz = pbuffer.data(idx_ovl_gh + 191);

    auto ts_xzzz_xxxyy = pbuffer.data(idx_ovl_gh + 192);

    auto ts_xzzz_xxxyz = pbuffer.data(idx_ovl_gh + 193);

    auto ts_xzzz_xxxzz = pbuffer.data(idx_ovl_gh + 194);

    auto ts_xzzz_xxyyy = pbuffer.data(idx_ovl_gh + 195);

    auto ts_xzzz_xxyyz = pbuffer.data(idx_ovl_gh + 196);

    auto ts_xzzz_xxyzz = pbuffer.data(idx_ovl_gh + 197);

    auto ts_xzzz_xxzzz = pbuffer.data(idx_ovl_gh + 198);

    auto ts_xzzz_xyyyy = pbuffer.data(idx_ovl_gh + 199);

    auto ts_xzzz_xyyyz = pbuffer.data(idx_ovl_gh + 200);

    auto ts_xzzz_xyyzz = pbuffer.data(idx_ovl_gh + 201);

    auto ts_xzzz_xyzzz = pbuffer.data(idx_ovl_gh + 202);

    auto ts_xzzz_xzzzz = pbuffer.data(idx_ovl_gh + 203);

    auto ts_xzzz_yyyyy = pbuffer.data(idx_ovl_gh + 204);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_ovl_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_ovl_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_ovl_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_ovl_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_ovl_gh + 209);

    auto ts_yyyy_xxxxx = pbuffer.data(idx_ovl_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_ovl_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_ovl_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_ovl_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_ovl_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_ovl_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_ovl_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_ovl_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_ovl_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_ovl_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_ovl_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_ovl_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_ovl_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_ovl_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_ovl_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_ovl_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_ovl_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_ovl_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_ovl_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_ovl_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_ovl_gh + 230);

    auto ts_yyyz_xxxxx = pbuffer.data(idx_ovl_gh + 231);

    auto ts_yyyz_xxxxy = pbuffer.data(idx_ovl_gh + 232);

    auto ts_yyyz_xxxxz = pbuffer.data(idx_ovl_gh + 233);

    auto ts_yyyz_xxxyy = pbuffer.data(idx_ovl_gh + 234);

    auto ts_yyyz_xxxyz = pbuffer.data(idx_ovl_gh + 235);

    auto ts_yyyz_xxxzz = pbuffer.data(idx_ovl_gh + 236);

    auto ts_yyyz_xxyyy = pbuffer.data(idx_ovl_gh + 237);

    auto ts_yyyz_xxyyz = pbuffer.data(idx_ovl_gh + 238);

    auto ts_yyyz_xxyzz = pbuffer.data(idx_ovl_gh + 239);

    auto ts_yyyz_xxzzz = pbuffer.data(idx_ovl_gh + 240);

    auto ts_yyyz_xyyyy = pbuffer.data(idx_ovl_gh + 241);

    auto ts_yyyz_xyyyz = pbuffer.data(idx_ovl_gh + 242);

    auto ts_yyyz_xyyzz = pbuffer.data(idx_ovl_gh + 243);

    auto ts_yyyz_xyzzz = pbuffer.data(idx_ovl_gh + 244);

    auto ts_yyyz_xzzzz = pbuffer.data(idx_ovl_gh + 245);

    auto ts_yyyz_yyyyy = pbuffer.data(idx_ovl_gh + 246);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_ovl_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_ovl_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_ovl_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_ovl_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_ovl_gh + 251);

    auto ts_yyzz_xxxxx = pbuffer.data(idx_ovl_gh + 252);

    auto ts_yyzz_xxxxy = pbuffer.data(idx_ovl_gh + 253);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_ovl_gh + 254);

    auto ts_yyzz_xxxyy = pbuffer.data(idx_ovl_gh + 255);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_ovl_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_ovl_gh + 257);

    auto ts_yyzz_xxyyy = pbuffer.data(idx_ovl_gh + 258);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_ovl_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_ovl_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_ovl_gh + 261);

    auto ts_yyzz_xyyyy = pbuffer.data(idx_ovl_gh + 262);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_ovl_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_ovl_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_ovl_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_ovl_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_ovl_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_ovl_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_ovl_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_ovl_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_ovl_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_ovl_gh + 272);

    auto ts_yzzz_xxxxx = pbuffer.data(idx_ovl_gh + 273);

    auto ts_yzzz_xxxxy = pbuffer.data(idx_ovl_gh + 274);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_ovl_gh + 275);

    auto ts_yzzz_xxxyy = pbuffer.data(idx_ovl_gh + 276);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_ovl_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_ovl_gh + 278);

    auto ts_yzzz_xxyyy = pbuffer.data(idx_ovl_gh + 279);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_ovl_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_ovl_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_ovl_gh + 282);

    auto ts_yzzz_xyyyy = pbuffer.data(idx_ovl_gh + 283);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_ovl_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_ovl_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_ovl_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_ovl_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_ovl_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_ovl_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_ovl_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_ovl_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_ovl_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_ovl_gh + 293);

    auto ts_zzzz_xxxxx = pbuffer.data(idx_ovl_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_ovl_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_ovl_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_ovl_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_ovl_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_ovl_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_ovl_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_ovl_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_ovl_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_ovl_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_ovl_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_ovl_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_ovl_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_ovl_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_ovl_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_ovl_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_ovl_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_ovl_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_ovl_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_ovl_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_ovl_gh + 314);

    // Set up 0-21 components of targeted buffer : GH

    auto tk_xxxx_xxxxx = pbuffer.data(idx_kin_gh);

    auto tk_xxxx_xxxxy = pbuffer.data(idx_kin_gh + 1);

    auto tk_xxxx_xxxxz = pbuffer.data(idx_kin_gh + 2);

    auto tk_xxxx_xxxyy = pbuffer.data(idx_kin_gh + 3);

    auto tk_xxxx_xxxyz = pbuffer.data(idx_kin_gh + 4);

    auto tk_xxxx_xxxzz = pbuffer.data(idx_kin_gh + 5);

    auto tk_xxxx_xxyyy = pbuffer.data(idx_kin_gh + 6);

    auto tk_xxxx_xxyyz = pbuffer.data(idx_kin_gh + 7);

    auto tk_xxxx_xxyzz = pbuffer.data(idx_kin_gh + 8);

    auto tk_xxxx_xxzzz = pbuffer.data(idx_kin_gh + 9);

    auto tk_xxxx_xyyyy = pbuffer.data(idx_kin_gh + 10);

    auto tk_xxxx_xyyyz = pbuffer.data(idx_kin_gh + 11);

    auto tk_xxxx_xyyzz = pbuffer.data(idx_kin_gh + 12);

    auto tk_xxxx_xyzzz = pbuffer.data(idx_kin_gh + 13);

    auto tk_xxxx_xzzzz = pbuffer.data(idx_kin_gh + 14);

    auto tk_xxxx_yyyyy = pbuffer.data(idx_kin_gh + 15);

    auto tk_xxxx_yyyyz = pbuffer.data(idx_kin_gh + 16);

    auto tk_xxxx_yyyzz = pbuffer.data(idx_kin_gh + 17);

    auto tk_xxxx_yyzzz = pbuffer.data(idx_kin_gh + 18);

    auto tk_xxxx_yzzzz = pbuffer.data(idx_kin_gh + 19);

    auto tk_xxxx_zzzzz = pbuffer.data(idx_kin_gh + 20);

#pragma omp simd aligned(pa_x,              \
                             tk_xx_xxxxx,   \
                             tk_xx_xxxxy,   \
                             tk_xx_xxxxz,   \
                             tk_xx_xxxyy,   \
                             tk_xx_xxxyz,   \
                             tk_xx_xxxzz,   \
                             tk_xx_xxyyy,   \
                             tk_xx_xxyyz,   \
                             tk_xx_xxyzz,   \
                             tk_xx_xxzzz,   \
                             tk_xx_xyyyy,   \
                             tk_xx_xyyyz,   \
                             tk_xx_xyyzz,   \
                             tk_xx_xyzzz,   \
                             tk_xx_xzzzz,   \
                             tk_xx_yyyyy,   \
                             tk_xx_yyyyz,   \
                             tk_xx_yyyzz,   \
                             tk_xx_yyzzz,   \
                             tk_xx_yzzzz,   \
                             tk_xx_zzzzz,   \
                             tk_xxx_xxxx,   \
                             tk_xxx_xxxxx,  \
                             tk_xxx_xxxxy,  \
                             tk_xxx_xxxxz,  \
                             tk_xxx_xxxy,   \
                             tk_xxx_xxxyy,  \
                             tk_xxx_xxxyz,  \
                             tk_xxx_xxxz,   \
                             tk_xxx_xxxzz,  \
                             tk_xxx_xxyy,   \
                             tk_xxx_xxyyy,  \
                             tk_xxx_xxyyz,  \
                             tk_xxx_xxyz,   \
                             tk_xxx_xxyzz,  \
                             tk_xxx_xxzz,   \
                             tk_xxx_xxzzz,  \
                             tk_xxx_xyyy,   \
                             tk_xxx_xyyyy,  \
                             tk_xxx_xyyyz,  \
                             tk_xxx_xyyz,   \
                             tk_xxx_xyyzz,  \
                             tk_xxx_xyzz,   \
                             tk_xxx_xyzzz,  \
                             tk_xxx_xzzz,   \
                             tk_xxx_xzzzz,  \
                             tk_xxx_yyyy,   \
                             tk_xxx_yyyyy,  \
                             tk_xxx_yyyyz,  \
                             tk_xxx_yyyz,   \
                             tk_xxx_yyyzz,  \
                             tk_xxx_yyzz,   \
                             tk_xxx_yyzzz,  \
                             tk_xxx_yzzz,   \
                             tk_xxx_yzzzz,  \
                             tk_xxx_zzzz,   \
                             tk_xxx_zzzzz,  \
                             tk_xxxx_xxxxx, \
                             tk_xxxx_xxxxy, \
                             tk_xxxx_xxxxz, \
                             tk_xxxx_xxxyy, \
                             tk_xxxx_xxxyz, \
                             tk_xxxx_xxxzz, \
                             tk_xxxx_xxyyy, \
                             tk_xxxx_xxyyz, \
                             tk_xxxx_xxyzz, \
                             tk_xxxx_xxzzz, \
                             tk_xxxx_xyyyy, \
                             tk_xxxx_xyyyz, \
                             tk_xxxx_xyyzz, \
                             tk_xxxx_xyzzz, \
                             tk_xxxx_xzzzz, \
                             tk_xxxx_yyyyy, \
                             tk_xxxx_yyyyz, \
                             tk_xxxx_yyyzz, \
                             tk_xxxx_yyzzz, \
                             tk_xxxx_yzzzz, \
                             tk_xxxx_zzzzz, \
                             ts_xx_xxxxx,   \
                             ts_xx_xxxxy,   \
                             ts_xx_xxxxz,   \
                             ts_xx_xxxyy,   \
                             ts_xx_xxxyz,   \
                             ts_xx_xxxzz,   \
                             ts_xx_xxyyy,   \
                             ts_xx_xxyyz,   \
                             ts_xx_xxyzz,   \
                             ts_xx_xxzzz,   \
                             ts_xx_xyyyy,   \
                             ts_xx_xyyyz,   \
                             ts_xx_xyyzz,   \
                             ts_xx_xyzzz,   \
                             ts_xx_xzzzz,   \
                             ts_xx_yyyyy,   \
                             ts_xx_yyyyz,   \
                             ts_xx_yyyzz,   \
                             ts_xx_yyzzz,   \
                             ts_xx_yzzzz,   \
                             ts_xx_zzzzz,   \
                             ts_xxxx_xxxxx, \
                             ts_xxxx_xxxxy, \
                             ts_xxxx_xxxxz, \
                             ts_xxxx_xxxyy, \
                             ts_xxxx_xxxyz, \
                             ts_xxxx_xxxzz, \
                             ts_xxxx_xxyyy, \
                             ts_xxxx_xxyyz, \
                             ts_xxxx_xxyzz, \
                             ts_xxxx_xxzzz, \
                             ts_xxxx_xyyyy, \
                             ts_xxxx_xyyyz, \
                             ts_xxxx_xyyzz, \
                             ts_xxxx_xyzzz, \
                             ts_xxxx_xzzzz, \
                             ts_xxxx_yyyyy, \
                             ts_xxxx_yyyyz, \
                             ts_xxxx_yyyzz, \
                             ts_xxxx_yyzzz, \
                             ts_xxxx_yzzzz, \
                             ts_xxxx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxx_xxxxx[i] = -6.0 * ts_xx_xxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxx[i] * fe_0 + 5.0 * tk_xxx_xxxx[i] * fe_0 +
                           tk_xxx_xxxxx[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxx[i] * fz_0;

        tk_xxxx_xxxxy[i] = -6.0 * ts_xx_xxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxy[i] * fe_0 + 4.0 * tk_xxx_xxxy[i] * fe_0 +
                           tk_xxx_xxxxy[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxy[i] * fz_0;

        tk_xxxx_xxxxz[i] = -6.0 * ts_xx_xxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxz[i] * fe_0 + 4.0 * tk_xxx_xxxz[i] * fe_0 +
                           tk_xxx_xxxxz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxz[i] * fz_0;

        tk_xxxx_xxxyy[i] = -6.0 * ts_xx_xxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxyy[i] * fe_0 + 3.0 * tk_xxx_xxyy[i] * fe_0 +
                           tk_xxx_xxxyy[i] * pa_x[i] + 2.0 * ts_xxxx_xxxyy[i] * fz_0;

        tk_xxxx_xxxyz[i] = -6.0 * ts_xx_xxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxyz[i] * fe_0 + 3.0 * tk_xxx_xxyz[i] * fe_0 +
                           tk_xxx_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxyz[i] * fz_0;

        tk_xxxx_xxxzz[i] = -6.0 * ts_xx_xxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxzz[i] * fe_0 + 3.0 * tk_xxx_xxzz[i] * fe_0 +
                           tk_xxx_xxxzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxzz[i] * fz_0;

        tk_xxxx_xxyyy[i] = -6.0 * ts_xx_xxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyyy[i] * fe_0 + 2.0 * tk_xxx_xyyy[i] * fe_0 +
                           tk_xxx_xxyyy[i] * pa_x[i] + 2.0 * ts_xxxx_xxyyy[i] * fz_0;

        tk_xxxx_xxyyz[i] = -6.0 * ts_xx_xxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyyz[i] * fe_0 + 2.0 * tk_xxx_xyyz[i] * fe_0 +
                           tk_xxx_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxx_xxyyz[i] * fz_0;

        tk_xxxx_xxyzz[i] = -6.0 * ts_xx_xxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyzz[i] * fe_0 + 2.0 * tk_xxx_xyzz[i] * fe_0 +
                           tk_xxx_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxyzz[i] * fz_0;

        tk_xxxx_xxzzz[i] = -6.0 * ts_xx_xxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxzzz[i] * fe_0 + 2.0 * tk_xxx_xzzz[i] * fe_0 +
                           tk_xxx_xxzzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxzzz[i] * fz_0;

        tk_xxxx_xyyyy[i] = -6.0 * ts_xx_xyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyyy[i] * fe_0 + tk_xxx_yyyy[i] * fe_0 + tk_xxx_xyyyy[i] * pa_x[i] +
                           2.0 * ts_xxxx_xyyyy[i] * fz_0;

        tk_xxxx_xyyyz[i] = -6.0 * ts_xx_xyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyyz[i] * fe_0 + tk_xxx_yyyz[i] * fe_0 + tk_xxx_xyyyz[i] * pa_x[i] +
                           2.0 * ts_xxxx_xyyyz[i] * fz_0;

        tk_xxxx_xyyzz[i] = -6.0 * ts_xx_xyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyzz[i] * fe_0 + tk_xxx_yyzz[i] * fe_0 + tk_xxx_xyyzz[i] * pa_x[i] +
                           2.0 * ts_xxxx_xyyzz[i] * fz_0;

        tk_xxxx_xyzzz[i] = -6.0 * ts_xx_xyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyzzz[i] * fe_0 + tk_xxx_yzzz[i] * fe_0 + tk_xxx_xyzzz[i] * pa_x[i] +
                           2.0 * ts_xxxx_xyzzz[i] * fz_0;

        tk_xxxx_xzzzz[i] = -6.0 * ts_xx_xzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xzzzz[i] * fe_0 + tk_xxx_zzzz[i] * fe_0 + tk_xxx_xzzzz[i] * pa_x[i] +
                           2.0 * ts_xxxx_xzzzz[i] * fz_0;

        tk_xxxx_yyyyy[i] =
            -6.0 * ts_xx_yyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyyy[i] * fe_0 + tk_xxx_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxx_yyyyy[i] * fz_0;

        tk_xxxx_yyyyz[i] =
            -6.0 * ts_xx_yyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyyz[i] * fe_0 + tk_xxx_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxx_yyyyz[i] * fz_0;

        tk_xxxx_yyyzz[i] =
            -6.0 * ts_xx_yyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyzz[i] * fe_0 + tk_xxx_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxx_yyyzz[i] * fz_0;

        tk_xxxx_yyzzz[i] =
            -6.0 * ts_xx_yyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyzzz[i] * fe_0 + tk_xxx_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxx_yyzzz[i] * fz_0;

        tk_xxxx_yzzzz[i] =
            -6.0 * ts_xx_yzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yzzzz[i] * fe_0 + tk_xxx_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_yzzzz[i] * fz_0;

        tk_xxxx_zzzzz[i] =
            -6.0 * ts_xx_zzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_zzzzz[i] * fe_0 + tk_xxx_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_zzzzz[i] * fz_0;
    }

    // Set up 21-42 components of targeted buffer : GH

    auto tk_xxxy_xxxxx = pbuffer.data(idx_kin_gh + 21);

    auto tk_xxxy_xxxxy = pbuffer.data(idx_kin_gh + 22);

    auto tk_xxxy_xxxxz = pbuffer.data(idx_kin_gh + 23);

    auto tk_xxxy_xxxyy = pbuffer.data(idx_kin_gh + 24);

    auto tk_xxxy_xxxyz = pbuffer.data(idx_kin_gh + 25);

    auto tk_xxxy_xxxzz = pbuffer.data(idx_kin_gh + 26);

    auto tk_xxxy_xxyyy = pbuffer.data(idx_kin_gh + 27);

    auto tk_xxxy_xxyyz = pbuffer.data(idx_kin_gh + 28);

    auto tk_xxxy_xxyzz = pbuffer.data(idx_kin_gh + 29);

    auto tk_xxxy_xxzzz = pbuffer.data(idx_kin_gh + 30);

    auto tk_xxxy_xyyyy = pbuffer.data(idx_kin_gh + 31);

    auto tk_xxxy_xyyyz = pbuffer.data(idx_kin_gh + 32);

    auto tk_xxxy_xyyzz = pbuffer.data(idx_kin_gh + 33);

    auto tk_xxxy_xyzzz = pbuffer.data(idx_kin_gh + 34);

    auto tk_xxxy_xzzzz = pbuffer.data(idx_kin_gh + 35);

    auto tk_xxxy_yyyyy = pbuffer.data(idx_kin_gh + 36);

    auto tk_xxxy_yyyyz = pbuffer.data(idx_kin_gh + 37);

    auto tk_xxxy_yyyzz = pbuffer.data(idx_kin_gh + 38);

    auto tk_xxxy_yyzzz = pbuffer.data(idx_kin_gh + 39);

    auto tk_xxxy_yzzzz = pbuffer.data(idx_kin_gh + 40);

    auto tk_xxxy_zzzzz = pbuffer.data(idx_kin_gh + 41);

#pragma omp simd aligned(pa_y,              \
                             tk_xxx_xxxx,   \
                             tk_xxx_xxxxx,  \
                             tk_xxx_xxxxy,  \
                             tk_xxx_xxxxz,  \
                             tk_xxx_xxxy,   \
                             tk_xxx_xxxyy,  \
                             tk_xxx_xxxyz,  \
                             tk_xxx_xxxz,   \
                             tk_xxx_xxxzz,  \
                             tk_xxx_xxyy,   \
                             tk_xxx_xxyyy,  \
                             tk_xxx_xxyyz,  \
                             tk_xxx_xxyz,   \
                             tk_xxx_xxyzz,  \
                             tk_xxx_xxzz,   \
                             tk_xxx_xxzzz,  \
                             tk_xxx_xyyy,   \
                             tk_xxx_xyyyy,  \
                             tk_xxx_xyyyz,  \
                             tk_xxx_xyyz,   \
                             tk_xxx_xyyzz,  \
                             tk_xxx_xyzz,   \
                             tk_xxx_xyzzz,  \
                             tk_xxx_xzzz,   \
                             tk_xxx_xzzzz,  \
                             tk_xxx_yyyy,   \
                             tk_xxx_yyyyy,  \
                             tk_xxx_yyyyz,  \
                             tk_xxx_yyyz,   \
                             tk_xxx_yyyzz,  \
                             tk_xxx_yyzz,   \
                             tk_xxx_yyzzz,  \
                             tk_xxx_yzzz,   \
                             tk_xxx_yzzzz,  \
                             tk_xxx_zzzz,   \
                             tk_xxx_zzzzz,  \
                             tk_xxxy_xxxxx, \
                             tk_xxxy_xxxxy, \
                             tk_xxxy_xxxxz, \
                             tk_xxxy_xxxyy, \
                             tk_xxxy_xxxyz, \
                             tk_xxxy_xxxzz, \
                             tk_xxxy_xxyyy, \
                             tk_xxxy_xxyyz, \
                             tk_xxxy_xxyzz, \
                             tk_xxxy_xxzzz, \
                             tk_xxxy_xyyyy, \
                             tk_xxxy_xyyyz, \
                             tk_xxxy_xyyzz, \
                             tk_xxxy_xyzzz, \
                             tk_xxxy_xzzzz, \
                             tk_xxxy_yyyyy, \
                             tk_xxxy_yyyyz, \
                             tk_xxxy_yyyzz, \
                             tk_xxxy_yyzzz, \
                             tk_xxxy_yzzzz, \
                             tk_xxxy_zzzzz, \
                             ts_xxxy_xxxxx, \
                             ts_xxxy_xxxxy, \
                             ts_xxxy_xxxxz, \
                             ts_xxxy_xxxyy, \
                             ts_xxxy_xxxyz, \
                             ts_xxxy_xxxzz, \
                             ts_xxxy_xxyyy, \
                             ts_xxxy_xxyyz, \
                             ts_xxxy_xxyzz, \
                             ts_xxxy_xxzzz, \
                             ts_xxxy_xyyyy, \
                             ts_xxxy_xyyyz, \
                             ts_xxxy_xyyzz, \
                             ts_xxxy_xyzzz, \
                             ts_xxxy_xzzzz, \
                             ts_xxxy_yyyyy, \
                             ts_xxxy_yyyyz, \
                             ts_xxxy_yyyzz, \
                             ts_xxxy_yyzzz, \
                             ts_xxxy_yzzzz, \
                             ts_xxxy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxy_xxxxx[i] = tk_xxx_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxx[i] * fz_0;

        tk_xxxy_xxxxy[i] = tk_xxx_xxxx[i] * fe_0 + tk_xxx_xxxxy[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxy[i] * fz_0;

        tk_xxxy_xxxxz[i] = tk_xxx_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxz[i] * fz_0;

        tk_xxxy_xxxyy[i] = 2.0 * tk_xxx_xxxy[i] * fe_0 + tk_xxx_xxxyy[i] * pa_y[i] + 2.0 * ts_xxxy_xxxyy[i] * fz_0;

        tk_xxxy_xxxyz[i] = tk_xxx_xxxz[i] * fe_0 + tk_xxx_xxxyz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxyz[i] * fz_0;

        tk_xxxy_xxxzz[i] = tk_xxx_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxzz[i] * fz_0;

        tk_xxxy_xxyyy[i] = 3.0 * tk_xxx_xxyy[i] * fe_0 + tk_xxx_xxyyy[i] * pa_y[i] + 2.0 * ts_xxxy_xxyyy[i] * fz_0;

        tk_xxxy_xxyyz[i] = 2.0 * tk_xxx_xxyz[i] * fe_0 + tk_xxx_xxyyz[i] * pa_y[i] + 2.0 * ts_xxxy_xxyyz[i] * fz_0;

        tk_xxxy_xxyzz[i] = tk_xxx_xxzz[i] * fe_0 + tk_xxx_xxyzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxyzz[i] * fz_0;

        tk_xxxy_xxzzz[i] = tk_xxx_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxzzz[i] * fz_0;

        tk_xxxy_xyyyy[i] = 4.0 * tk_xxx_xyyy[i] * fe_0 + tk_xxx_xyyyy[i] * pa_y[i] + 2.0 * ts_xxxy_xyyyy[i] * fz_0;

        tk_xxxy_xyyyz[i] = 3.0 * tk_xxx_xyyz[i] * fe_0 + tk_xxx_xyyyz[i] * pa_y[i] + 2.0 * ts_xxxy_xyyyz[i] * fz_0;

        tk_xxxy_xyyzz[i] = 2.0 * tk_xxx_xyzz[i] * fe_0 + tk_xxx_xyyzz[i] * pa_y[i] + 2.0 * ts_xxxy_xyyzz[i] * fz_0;

        tk_xxxy_xyzzz[i] = tk_xxx_xzzz[i] * fe_0 + tk_xxx_xyzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xyzzz[i] * fz_0;

        tk_xxxy_xzzzz[i] = tk_xxx_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xzzzz[i] * fz_0;

        tk_xxxy_yyyyy[i] = 5.0 * tk_xxx_yyyy[i] * fe_0 + tk_xxx_yyyyy[i] * pa_y[i] + 2.0 * ts_xxxy_yyyyy[i] * fz_0;

        tk_xxxy_yyyyz[i] = 4.0 * tk_xxx_yyyz[i] * fe_0 + tk_xxx_yyyyz[i] * pa_y[i] + 2.0 * ts_xxxy_yyyyz[i] * fz_0;

        tk_xxxy_yyyzz[i] = 3.0 * tk_xxx_yyzz[i] * fe_0 + tk_xxx_yyyzz[i] * pa_y[i] + 2.0 * ts_xxxy_yyyzz[i] * fz_0;

        tk_xxxy_yyzzz[i] = 2.0 * tk_xxx_yzzz[i] * fe_0 + tk_xxx_yyzzz[i] * pa_y[i] + 2.0 * ts_xxxy_yyzzz[i] * fz_0;

        tk_xxxy_yzzzz[i] = tk_xxx_zzzz[i] * fe_0 + tk_xxx_yzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_yzzzz[i] * fz_0;

        tk_xxxy_zzzzz[i] = tk_xxx_zzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_zzzzz[i] * fz_0;
    }

    // Set up 42-63 components of targeted buffer : GH

    auto tk_xxxz_xxxxx = pbuffer.data(idx_kin_gh + 42);

    auto tk_xxxz_xxxxy = pbuffer.data(idx_kin_gh + 43);

    auto tk_xxxz_xxxxz = pbuffer.data(idx_kin_gh + 44);

    auto tk_xxxz_xxxyy = pbuffer.data(idx_kin_gh + 45);

    auto tk_xxxz_xxxyz = pbuffer.data(idx_kin_gh + 46);

    auto tk_xxxz_xxxzz = pbuffer.data(idx_kin_gh + 47);

    auto tk_xxxz_xxyyy = pbuffer.data(idx_kin_gh + 48);

    auto tk_xxxz_xxyyz = pbuffer.data(idx_kin_gh + 49);

    auto tk_xxxz_xxyzz = pbuffer.data(idx_kin_gh + 50);

    auto tk_xxxz_xxzzz = pbuffer.data(idx_kin_gh + 51);

    auto tk_xxxz_xyyyy = pbuffer.data(idx_kin_gh + 52);

    auto tk_xxxz_xyyyz = pbuffer.data(idx_kin_gh + 53);

    auto tk_xxxz_xyyzz = pbuffer.data(idx_kin_gh + 54);

    auto tk_xxxz_xyzzz = pbuffer.data(idx_kin_gh + 55);

    auto tk_xxxz_xzzzz = pbuffer.data(idx_kin_gh + 56);

    auto tk_xxxz_yyyyy = pbuffer.data(idx_kin_gh + 57);

    auto tk_xxxz_yyyyz = pbuffer.data(idx_kin_gh + 58);

    auto tk_xxxz_yyyzz = pbuffer.data(idx_kin_gh + 59);

    auto tk_xxxz_yyzzz = pbuffer.data(idx_kin_gh + 60);

    auto tk_xxxz_yzzzz = pbuffer.data(idx_kin_gh + 61);

    auto tk_xxxz_zzzzz = pbuffer.data(idx_kin_gh + 62);

#pragma omp simd aligned(pa_z,              \
                             tk_xxx_xxxx,   \
                             tk_xxx_xxxxx,  \
                             tk_xxx_xxxxy,  \
                             tk_xxx_xxxxz,  \
                             tk_xxx_xxxy,   \
                             tk_xxx_xxxyy,  \
                             tk_xxx_xxxyz,  \
                             tk_xxx_xxxz,   \
                             tk_xxx_xxxzz,  \
                             tk_xxx_xxyy,   \
                             tk_xxx_xxyyy,  \
                             tk_xxx_xxyyz,  \
                             tk_xxx_xxyz,   \
                             tk_xxx_xxyzz,  \
                             tk_xxx_xxzz,   \
                             tk_xxx_xxzzz,  \
                             tk_xxx_xyyy,   \
                             tk_xxx_xyyyy,  \
                             tk_xxx_xyyyz,  \
                             tk_xxx_xyyz,   \
                             tk_xxx_xyyzz,  \
                             tk_xxx_xyzz,   \
                             tk_xxx_xyzzz,  \
                             tk_xxx_xzzz,   \
                             tk_xxx_xzzzz,  \
                             tk_xxx_yyyy,   \
                             tk_xxx_yyyyy,  \
                             tk_xxx_yyyyz,  \
                             tk_xxx_yyyz,   \
                             tk_xxx_yyyzz,  \
                             tk_xxx_yyzz,   \
                             tk_xxx_yyzzz,  \
                             tk_xxx_yzzz,   \
                             tk_xxx_yzzzz,  \
                             tk_xxx_zzzz,   \
                             tk_xxx_zzzzz,  \
                             tk_xxxz_xxxxx, \
                             tk_xxxz_xxxxy, \
                             tk_xxxz_xxxxz, \
                             tk_xxxz_xxxyy, \
                             tk_xxxz_xxxyz, \
                             tk_xxxz_xxxzz, \
                             tk_xxxz_xxyyy, \
                             tk_xxxz_xxyyz, \
                             tk_xxxz_xxyzz, \
                             tk_xxxz_xxzzz, \
                             tk_xxxz_xyyyy, \
                             tk_xxxz_xyyyz, \
                             tk_xxxz_xyyzz, \
                             tk_xxxz_xyzzz, \
                             tk_xxxz_xzzzz, \
                             tk_xxxz_yyyyy, \
                             tk_xxxz_yyyyz, \
                             tk_xxxz_yyyzz, \
                             tk_xxxz_yyzzz, \
                             tk_xxxz_yzzzz, \
                             tk_xxxz_zzzzz, \
                             ts_xxxz_xxxxx, \
                             ts_xxxz_xxxxy, \
                             ts_xxxz_xxxxz, \
                             ts_xxxz_xxxyy, \
                             ts_xxxz_xxxyz, \
                             ts_xxxz_xxxzz, \
                             ts_xxxz_xxyyy, \
                             ts_xxxz_xxyyz, \
                             ts_xxxz_xxyzz, \
                             ts_xxxz_xxzzz, \
                             ts_xxxz_xyyyy, \
                             ts_xxxz_xyyyz, \
                             ts_xxxz_xyyzz, \
                             ts_xxxz_xyzzz, \
                             ts_xxxz_xzzzz, \
                             ts_xxxz_yyyyy, \
                             ts_xxxz_yyyyz, \
                             ts_xxxz_yyyzz, \
                             ts_xxxz_yyzzz, \
                             ts_xxxz_yzzzz, \
                             ts_xxxz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxz_xxxxx[i] = tk_xxx_xxxxx[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxx[i] * fz_0;

        tk_xxxz_xxxxy[i] = tk_xxx_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxy[i] * fz_0;

        tk_xxxz_xxxxz[i] = tk_xxx_xxxx[i] * fe_0 + tk_xxx_xxxxz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxz[i] * fz_0;

        tk_xxxz_xxxyy[i] = tk_xxx_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxz_xxxyy[i] * fz_0;

        tk_xxxz_xxxyz[i] = tk_xxx_xxxy[i] * fe_0 + tk_xxx_xxxyz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxyz[i] * fz_0;

        tk_xxxz_xxxzz[i] = 2.0 * tk_xxx_xxxz[i] * fe_0 + tk_xxx_xxxzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxzz[i] * fz_0;

        tk_xxxz_xxyyy[i] = tk_xxx_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxz_xxyyy[i] * fz_0;

        tk_xxxz_xxyyz[i] = tk_xxx_xxyy[i] * fe_0 + tk_xxx_xxyyz[i] * pa_z[i] + 2.0 * ts_xxxz_xxyyz[i] * fz_0;

        tk_xxxz_xxyzz[i] = 2.0 * tk_xxx_xxyz[i] * fe_0 + tk_xxx_xxyzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxyzz[i] * fz_0;

        tk_xxxz_xxzzz[i] = 3.0 * tk_xxx_xxzz[i] * fe_0 + tk_xxx_xxzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxzzz[i] * fz_0;

        tk_xxxz_xyyyy[i] = tk_xxx_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxz_xyyyy[i] * fz_0;

        tk_xxxz_xyyyz[i] = tk_xxx_xyyy[i] * fe_0 + tk_xxx_xyyyz[i] * pa_z[i] + 2.0 * ts_xxxz_xyyyz[i] * fz_0;

        tk_xxxz_xyyzz[i] = 2.0 * tk_xxx_xyyz[i] * fe_0 + tk_xxx_xyyzz[i] * pa_z[i] + 2.0 * ts_xxxz_xyyzz[i] * fz_0;

        tk_xxxz_xyzzz[i] = 3.0 * tk_xxx_xyzz[i] * fe_0 + tk_xxx_xyzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xyzzz[i] * fz_0;

        tk_xxxz_xzzzz[i] = 4.0 * tk_xxx_xzzz[i] * fe_0 + tk_xxx_xzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xzzzz[i] * fz_0;

        tk_xxxz_yyyyy[i] = tk_xxx_yyyyy[i] * pa_z[i] + 2.0 * ts_xxxz_yyyyy[i] * fz_0;

        tk_xxxz_yyyyz[i] = tk_xxx_yyyy[i] * fe_0 + tk_xxx_yyyyz[i] * pa_z[i] + 2.0 * ts_xxxz_yyyyz[i] * fz_0;

        tk_xxxz_yyyzz[i] = 2.0 * tk_xxx_yyyz[i] * fe_0 + tk_xxx_yyyzz[i] * pa_z[i] + 2.0 * ts_xxxz_yyyzz[i] * fz_0;

        tk_xxxz_yyzzz[i] = 3.0 * tk_xxx_yyzz[i] * fe_0 + tk_xxx_yyzzz[i] * pa_z[i] + 2.0 * ts_xxxz_yyzzz[i] * fz_0;

        tk_xxxz_yzzzz[i] = 4.0 * tk_xxx_yzzz[i] * fe_0 + tk_xxx_yzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_yzzzz[i] * fz_0;

        tk_xxxz_zzzzz[i] = 5.0 * tk_xxx_zzzz[i] * fe_0 + tk_xxx_zzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_zzzzz[i] * fz_0;
    }

    // Set up 63-84 components of targeted buffer : GH

    auto tk_xxyy_xxxxx = pbuffer.data(idx_kin_gh + 63);

    auto tk_xxyy_xxxxy = pbuffer.data(idx_kin_gh + 64);

    auto tk_xxyy_xxxxz = pbuffer.data(idx_kin_gh + 65);

    auto tk_xxyy_xxxyy = pbuffer.data(idx_kin_gh + 66);

    auto tk_xxyy_xxxyz = pbuffer.data(idx_kin_gh + 67);

    auto tk_xxyy_xxxzz = pbuffer.data(idx_kin_gh + 68);

    auto tk_xxyy_xxyyy = pbuffer.data(idx_kin_gh + 69);

    auto tk_xxyy_xxyyz = pbuffer.data(idx_kin_gh + 70);

    auto tk_xxyy_xxyzz = pbuffer.data(idx_kin_gh + 71);

    auto tk_xxyy_xxzzz = pbuffer.data(idx_kin_gh + 72);

    auto tk_xxyy_xyyyy = pbuffer.data(idx_kin_gh + 73);

    auto tk_xxyy_xyyyz = pbuffer.data(idx_kin_gh + 74);

    auto tk_xxyy_xyyzz = pbuffer.data(idx_kin_gh + 75);

    auto tk_xxyy_xyzzz = pbuffer.data(idx_kin_gh + 76);

    auto tk_xxyy_xzzzz = pbuffer.data(idx_kin_gh + 77);

    auto tk_xxyy_yyyyy = pbuffer.data(idx_kin_gh + 78);

    auto tk_xxyy_yyyyz = pbuffer.data(idx_kin_gh + 79);

    auto tk_xxyy_yyyzz = pbuffer.data(idx_kin_gh + 80);

    auto tk_xxyy_yyzzz = pbuffer.data(idx_kin_gh + 81);

    auto tk_xxyy_yzzzz = pbuffer.data(idx_kin_gh + 82);

    auto tk_xxyy_zzzzz = pbuffer.data(idx_kin_gh + 83);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xx_xxxxx,   \
                             tk_xx_xxxxz,   \
                             tk_xx_xxxzz,   \
                             tk_xx_xxzzz,   \
                             tk_xx_xzzzz,   \
                             tk_xxy_xxxxx,  \
                             tk_xxy_xxxxz,  \
                             tk_xxy_xxxzz,  \
                             tk_xxy_xxzzz,  \
                             tk_xxy_xzzzz,  \
                             tk_xxyy_xxxxx, \
                             tk_xxyy_xxxxy, \
                             tk_xxyy_xxxxz, \
                             tk_xxyy_xxxyy, \
                             tk_xxyy_xxxyz, \
                             tk_xxyy_xxxzz, \
                             tk_xxyy_xxyyy, \
                             tk_xxyy_xxyyz, \
                             tk_xxyy_xxyzz, \
                             tk_xxyy_xxzzz, \
                             tk_xxyy_xyyyy, \
                             tk_xxyy_xyyyz, \
                             tk_xxyy_xyyzz, \
                             tk_xxyy_xyzzz, \
                             tk_xxyy_xzzzz, \
                             tk_xxyy_yyyyy, \
                             tk_xxyy_yyyyz, \
                             tk_xxyy_yyyzz, \
                             tk_xxyy_yyzzz, \
                             tk_xxyy_yzzzz, \
                             tk_xxyy_zzzzz, \
                             tk_xyy_xxxxy,  \
                             tk_xyy_xxxy,   \
                             tk_xyy_xxxyy,  \
                             tk_xyy_xxxyz,  \
                             tk_xyy_xxyy,   \
                             tk_xyy_xxyyy,  \
                             tk_xyy_xxyyz,  \
                             tk_xyy_xxyz,   \
                             tk_xyy_xxyzz,  \
                             tk_xyy_xyyy,   \
                             tk_xyy_xyyyy,  \
                             tk_xyy_xyyyz,  \
                             tk_xyy_xyyz,   \
                             tk_xyy_xyyzz,  \
                             tk_xyy_xyzz,   \
                             tk_xyy_xyzzz,  \
                             tk_xyy_yyyy,   \
                             tk_xyy_yyyyy,  \
                             tk_xyy_yyyyz,  \
                             tk_xyy_yyyz,   \
                             tk_xyy_yyyzz,  \
                             tk_xyy_yyzz,   \
                             tk_xyy_yyzzz,  \
                             tk_xyy_yzzz,   \
                             tk_xyy_yzzzz,  \
                             tk_xyy_zzzzz,  \
                             tk_yy_xxxxy,   \
                             tk_yy_xxxyy,   \
                             tk_yy_xxxyz,   \
                             tk_yy_xxyyy,   \
                             tk_yy_xxyyz,   \
                             tk_yy_xxyzz,   \
                             tk_yy_xyyyy,   \
                             tk_yy_xyyyz,   \
                             tk_yy_xyyzz,   \
                             tk_yy_xyzzz,   \
                             tk_yy_yyyyy,   \
                             tk_yy_yyyyz,   \
                             tk_yy_yyyzz,   \
                             tk_yy_yyzzz,   \
                             tk_yy_yzzzz,   \
                             tk_yy_zzzzz,   \
                             ts_xx_xxxxx,   \
                             ts_xx_xxxxz,   \
                             ts_xx_xxxzz,   \
                             ts_xx_xxzzz,   \
                             ts_xx_xzzzz,   \
                             ts_xxyy_xxxxx, \
                             ts_xxyy_xxxxy, \
                             ts_xxyy_xxxxz, \
                             ts_xxyy_xxxyy, \
                             ts_xxyy_xxxyz, \
                             ts_xxyy_xxxzz, \
                             ts_xxyy_xxyyy, \
                             ts_xxyy_xxyyz, \
                             ts_xxyy_xxyzz, \
                             ts_xxyy_xxzzz, \
                             ts_xxyy_xyyyy, \
                             ts_xxyy_xyyyz, \
                             ts_xxyy_xyyzz, \
                             ts_xxyy_xyzzz, \
                             ts_xxyy_xzzzz, \
                             ts_xxyy_yyyyy, \
                             ts_xxyy_yyyyz, \
                             ts_xxyy_yyyzz, \
                             ts_xxyy_yyzzz, \
                             ts_xxyy_yzzzz, \
                             ts_xxyy_zzzzz, \
                             ts_yy_xxxxy,   \
                             ts_yy_xxxyy,   \
                             ts_yy_xxxyz,   \
                             ts_yy_xxyyy,   \
                             ts_yy_xxyyz,   \
                             ts_yy_xxyzz,   \
                             ts_yy_xyyyy,   \
                             ts_yy_xyyyz,   \
                             ts_yy_xyyzz,   \
                             ts_yy_xyzzz,   \
                             ts_yy_yyyyy,   \
                             ts_yy_yyyyz,   \
                             ts_yy_yyyzz,   \
                             ts_yy_yyzzz,   \
                             ts_yy_yzzzz,   \
                             ts_yy_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyy_xxxxx[i] = -2.0 * ts_xx_xxxxx[i] * fbe_0 * fz_0 + tk_xx_xxxxx[i] * fe_0 + tk_xxy_xxxxx[i] * pa_y[i] + 2.0 * ts_xxyy_xxxxx[i] * fz_0;

        tk_xxyy_xxxxy[i] = -2.0 * ts_yy_xxxxy[i] * fbe_0 * fz_0 + tk_yy_xxxxy[i] * fe_0 + 4.0 * tk_xyy_xxxy[i] * fe_0 + tk_xyy_xxxxy[i] * pa_x[i] +
                           2.0 * ts_xxyy_xxxxy[i] * fz_0;

        tk_xxyy_xxxxz[i] = -2.0 * ts_xx_xxxxz[i] * fbe_0 * fz_0 + tk_xx_xxxxz[i] * fe_0 + tk_xxy_xxxxz[i] * pa_y[i] + 2.0 * ts_xxyy_xxxxz[i] * fz_0;

        tk_xxyy_xxxyy[i] = -2.0 * ts_yy_xxxyy[i] * fbe_0 * fz_0 + tk_yy_xxxyy[i] * fe_0 + 3.0 * tk_xyy_xxyy[i] * fe_0 + tk_xyy_xxxyy[i] * pa_x[i] +
                           2.0 * ts_xxyy_xxxyy[i] * fz_0;

        tk_xxyy_xxxyz[i] = -2.0 * ts_yy_xxxyz[i] * fbe_0 * fz_0 + tk_yy_xxxyz[i] * fe_0 + 3.0 * tk_xyy_xxyz[i] * fe_0 + tk_xyy_xxxyz[i] * pa_x[i] +
                           2.0 * ts_xxyy_xxxyz[i] * fz_0;

        tk_xxyy_xxxzz[i] = -2.0 * ts_xx_xxxzz[i] * fbe_0 * fz_0 + tk_xx_xxxzz[i] * fe_0 + tk_xxy_xxxzz[i] * pa_y[i] + 2.0 * ts_xxyy_xxxzz[i] * fz_0;

        tk_xxyy_xxyyy[i] = -2.0 * ts_yy_xxyyy[i] * fbe_0 * fz_0 + tk_yy_xxyyy[i] * fe_0 + 2.0 * tk_xyy_xyyy[i] * fe_0 + tk_xyy_xxyyy[i] * pa_x[i] +
                           2.0 * ts_xxyy_xxyyy[i] * fz_0;

        tk_xxyy_xxyyz[i] = -2.0 * ts_yy_xxyyz[i] * fbe_0 * fz_0 + tk_yy_xxyyz[i] * fe_0 + 2.0 * tk_xyy_xyyz[i] * fe_0 + tk_xyy_xxyyz[i] * pa_x[i] +
                           2.0 * ts_xxyy_xxyyz[i] * fz_0;

        tk_xxyy_xxyzz[i] = -2.0 * ts_yy_xxyzz[i] * fbe_0 * fz_0 + tk_yy_xxyzz[i] * fe_0 + 2.0 * tk_xyy_xyzz[i] * fe_0 + tk_xyy_xxyzz[i] * pa_x[i] +
                           2.0 * ts_xxyy_xxyzz[i] * fz_0;

        tk_xxyy_xxzzz[i] = -2.0 * ts_xx_xxzzz[i] * fbe_0 * fz_0 + tk_xx_xxzzz[i] * fe_0 + tk_xxy_xxzzz[i] * pa_y[i] + 2.0 * ts_xxyy_xxzzz[i] * fz_0;

        tk_xxyy_xyyyy[i] = -2.0 * ts_yy_xyyyy[i] * fbe_0 * fz_0 + tk_yy_xyyyy[i] * fe_0 + tk_xyy_yyyy[i] * fe_0 + tk_xyy_xyyyy[i] * pa_x[i] +
                           2.0 * ts_xxyy_xyyyy[i] * fz_0;

        tk_xxyy_xyyyz[i] = -2.0 * ts_yy_xyyyz[i] * fbe_0 * fz_0 + tk_yy_xyyyz[i] * fe_0 + tk_xyy_yyyz[i] * fe_0 + tk_xyy_xyyyz[i] * pa_x[i] +
                           2.0 * ts_xxyy_xyyyz[i] * fz_0;

        tk_xxyy_xyyzz[i] = -2.0 * ts_yy_xyyzz[i] * fbe_0 * fz_0 + tk_yy_xyyzz[i] * fe_0 + tk_xyy_yyzz[i] * fe_0 + tk_xyy_xyyzz[i] * pa_x[i] +
                           2.0 * ts_xxyy_xyyzz[i] * fz_0;

        tk_xxyy_xyzzz[i] = -2.0 * ts_yy_xyzzz[i] * fbe_0 * fz_0 + tk_yy_xyzzz[i] * fe_0 + tk_xyy_yzzz[i] * fe_0 + tk_xyy_xyzzz[i] * pa_x[i] +
                           2.0 * ts_xxyy_xyzzz[i] * fz_0;

        tk_xxyy_xzzzz[i] = -2.0 * ts_xx_xzzzz[i] * fbe_0 * fz_0 + tk_xx_xzzzz[i] * fe_0 + tk_xxy_xzzzz[i] * pa_y[i] + 2.0 * ts_xxyy_xzzzz[i] * fz_0;

        tk_xxyy_yyyyy[i] = -2.0 * ts_yy_yyyyy[i] * fbe_0 * fz_0 + tk_yy_yyyyy[i] * fe_0 + tk_xyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xxyy_yyyyy[i] * fz_0;

        tk_xxyy_yyyyz[i] = -2.0 * ts_yy_yyyyz[i] * fbe_0 * fz_0 + tk_yy_yyyyz[i] * fe_0 + tk_xyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xxyy_yyyyz[i] * fz_0;

        tk_xxyy_yyyzz[i] = -2.0 * ts_yy_yyyzz[i] * fbe_0 * fz_0 + tk_yy_yyyzz[i] * fe_0 + tk_xyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xxyy_yyyzz[i] * fz_0;

        tk_xxyy_yyzzz[i] = -2.0 * ts_yy_yyzzz[i] * fbe_0 * fz_0 + tk_yy_yyzzz[i] * fe_0 + tk_xyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xxyy_yyzzz[i] * fz_0;

        tk_xxyy_yzzzz[i] = -2.0 * ts_yy_yzzzz[i] * fbe_0 * fz_0 + tk_yy_yzzzz[i] * fe_0 + tk_xyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xxyy_yzzzz[i] * fz_0;

        tk_xxyy_zzzzz[i] = -2.0 * ts_yy_zzzzz[i] * fbe_0 * fz_0 + tk_yy_zzzzz[i] * fe_0 + tk_xyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xxyy_zzzzz[i] * fz_0;
    }

    // Set up 84-105 components of targeted buffer : GH

    auto tk_xxyz_xxxxx = pbuffer.data(idx_kin_gh + 84);

    auto tk_xxyz_xxxxy = pbuffer.data(idx_kin_gh + 85);

    auto tk_xxyz_xxxxz = pbuffer.data(idx_kin_gh + 86);

    auto tk_xxyz_xxxyy = pbuffer.data(idx_kin_gh + 87);

    auto tk_xxyz_xxxyz = pbuffer.data(idx_kin_gh + 88);

    auto tk_xxyz_xxxzz = pbuffer.data(idx_kin_gh + 89);

    auto tk_xxyz_xxyyy = pbuffer.data(idx_kin_gh + 90);

    auto tk_xxyz_xxyyz = pbuffer.data(idx_kin_gh + 91);

    auto tk_xxyz_xxyzz = pbuffer.data(idx_kin_gh + 92);

    auto tk_xxyz_xxzzz = pbuffer.data(idx_kin_gh + 93);

    auto tk_xxyz_xyyyy = pbuffer.data(idx_kin_gh + 94);

    auto tk_xxyz_xyyyz = pbuffer.data(idx_kin_gh + 95);

    auto tk_xxyz_xyyzz = pbuffer.data(idx_kin_gh + 96);

    auto tk_xxyz_xyzzz = pbuffer.data(idx_kin_gh + 97);

    auto tk_xxyz_xzzzz = pbuffer.data(idx_kin_gh + 98);

    auto tk_xxyz_yyyyy = pbuffer.data(idx_kin_gh + 99);

    auto tk_xxyz_yyyyz = pbuffer.data(idx_kin_gh + 100);

    auto tk_xxyz_yyyzz = pbuffer.data(idx_kin_gh + 101);

    auto tk_xxyz_yyzzz = pbuffer.data(idx_kin_gh + 102);

    auto tk_xxyz_yzzzz = pbuffer.data(idx_kin_gh + 103);

    auto tk_xxyz_zzzzz = pbuffer.data(idx_kin_gh + 104);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_xxy_xxxxy,  \
                             tk_xxy_xxxyy,  \
                             tk_xxy_xxyyy,  \
                             tk_xxy_xyyyy,  \
                             tk_xxy_yyyyy,  \
                             tk_xxyz_xxxxx, \
                             tk_xxyz_xxxxy, \
                             tk_xxyz_xxxxz, \
                             tk_xxyz_xxxyy, \
                             tk_xxyz_xxxyz, \
                             tk_xxyz_xxxzz, \
                             tk_xxyz_xxyyy, \
                             tk_xxyz_xxyyz, \
                             tk_xxyz_xxyzz, \
                             tk_xxyz_xxzzz, \
                             tk_xxyz_xyyyy, \
                             tk_xxyz_xyyyz, \
                             tk_xxyz_xyyzz, \
                             tk_xxyz_xyzzz, \
                             tk_xxyz_xzzzz, \
                             tk_xxyz_yyyyy, \
                             tk_xxyz_yyyyz, \
                             tk_xxyz_yyyzz, \
                             tk_xxyz_yyzzz, \
                             tk_xxyz_yzzzz, \
                             tk_xxyz_zzzzz, \
                             tk_xxz_xxxxx,  \
                             tk_xxz_xxxxz,  \
                             tk_xxz_xxxyz,  \
                             tk_xxz_xxxz,   \
                             tk_xxz_xxxzz,  \
                             tk_xxz_xxyyz,  \
                             tk_xxz_xxyz,   \
                             tk_xxz_xxyzz,  \
                             tk_xxz_xxzz,   \
                             tk_xxz_xxzzz,  \
                             tk_xxz_xyyyz,  \
                             tk_xxz_xyyz,   \
                             tk_xxz_xyyzz,  \
                             tk_xxz_xyzz,   \
                             tk_xxz_xyzzz,  \
                             tk_xxz_xzzz,   \
                             tk_xxz_xzzzz,  \
                             tk_xxz_yyyyz,  \
                             tk_xxz_yyyz,   \
                             tk_xxz_yyyzz,  \
                             tk_xxz_yyzz,   \
                             tk_xxz_yyzzz,  \
                             tk_xxz_yzzz,   \
                             tk_xxz_yzzzz,  \
                             tk_xxz_zzzz,   \
                             tk_xxz_zzzzz,  \
                             ts_xxyz_xxxxx, \
                             ts_xxyz_xxxxy, \
                             ts_xxyz_xxxxz, \
                             ts_xxyz_xxxyy, \
                             ts_xxyz_xxxyz, \
                             ts_xxyz_xxxzz, \
                             ts_xxyz_xxyyy, \
                             ts_xxyz_xxyyz, \
                             ts_xxyz_xxyzz, \
                             ts_xxyz_xxzzz, \
                             ts_xxyz_xyyyy, \
                             ts_xxyz_xyyyz, \
                             ts_xxyz_xyyzz, \
                             ts_xxyz_xyzzz, \
                             ts_xxyz_xzzzz, \
                             ts_xxyz_yyyyy, \
                             ts_xxyz_yyyyz, \
                             ts_xxyz_yyyzz, \
                             ts_xxyz_yyzzz, \
                             ts_xxyz_yzzzz, \
                             ts_xxyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyz_xxxxx[i] = tk_xxz_xxxxx[i] * pa_y[i] + 2.0 * ts_xxyz_xxxxx[i] * fz_0;

        tk_xxyz_xxxxy[i] = tk_xxy_xxxxy[i] * pa_z[i] + 2.0 * ts_xxyz_xxxxy[i] * fz_0;

        tk_xxyz_xxxxz[i] = tk_xxz_xxxxz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxxz[i] * fz_0;

        tk_xxyz_xxxyy[i] = tk_xxy_xxxyy[i] * pa_z[i] + 2.0 * ts_xxyz_xxxyy[i] * fz_0;

        tk_xxyz_xxxyz[i] = tk_xxz_xxxz[i] * fe_0 + tk_xxz_xxxyz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxyz[i] * fz_0;

        tk_xxyz_xxxzz[i] = tk_xxz_xxxzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxzz[i] * fz_0;

        tk_xxyz_xxyyy[i] = tk_xxy_xxyyy[i] * pa_z[i] + 2.0 * ts_xxyz_xxyyy[i] * fz_0;

        tk_xxyz_xxyyz[i] = 2.0 * tk_xxz_xxyz[i] * fe_0 + tk_xxz_xxyyz[i] * pa_y[i] + 2.0 * ts_xxyz_xxyyz[i] * fz_0;

        tk_xxyz_xxyzz[i] = tk_xxz_xxzz[i] * fe_0 + tk_xxz_xxyzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxyzz[i] * fz_0;

        tk_xxyz_xxzzz[i] = tk_xxz_xxzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxzzz[i] * fz_0;

        tk_xxyz_xyyyy[i] = tk_xxy_xyyyy[i] * pa_z[i] + 2.0 * ts_xxyz_xyyyy[i] * fz_0;

        tk_xxyz_xyyyz[i] = 3.0 * tk_xxz_xyyz[i] * fe_0 + tk_xxz_xyyyz[i] * pa_y[i] + 2.0 * ts_xxyz_xyyyz[i] * fz_0;

        tk_xxyz_xyyzz[i] = 2.0 * tk_xxz_xyzz[i] * fe_0 + tk_xxz_xyyzz[i] * pa_y[i] + 2.0 * ts_xxyz_xyyzz[i] * fz_0;

        tk_xxyz_xyzzz[i] = tk_xxz_xzzz[i] * fe_0 + tk_xxz_xyzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xyzzz[i] * fz_0;

        tk_xxyz_xzzzz[i] = tk_xxz_xzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xzzzz[i] * fz_0;

        tk_xxyz_yyyyy[i] = tk_xxy_yyyyy[i] * pa_z[i] + 2.0 * ts_xxyz_yyyyy[i] * fz_0;

        tk_xxyz_yyyyz[i] = 4.0 * tk_xxz_yyyz[i] * fe_0 + tk_xxz_yyyyz[i] * pa_y[i] + 2.0 * ts_xxyz_yyyyz[i] * fz_0;

        tk_xxyz_yyyzz[i] = 3.0 * tk_xxz_yyzz[i] * fe_0 + tk_xxz_yyyzz[i] * pa_y[i] + 2.0 * ts_xxyz_yyyzz[i] * fz_0;

        tk_xxyz_yyzzz[i] = 2.0 * tk_xxz_yzzz[i] * fe_0 + tk_xxz_yyzzz[i] * pa_y[i] + 2.0 * ts_xxyz_yyzzz[i] * fz_0;

        tk_xxyz_yzzzz[i] = tk_xxz_zzzz[i] * fe_0 + tk_xxz_yzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_yzzzz[i] * fz_0;

        tk_xxyz_zzzzz[i] = tk_xxz_zzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_zzzzz[i] * fz_0;
    }

    // Set up 105-126 components of targeted buffer : GH

    auto tk_xxzz_xxxxx = pbuffer.data(idx_kin_gh + 105);

    auto tk_xxzz_xxxxy = pbuffer.data(idx_kin_gh + 106);

    auto tk_xxzz_xxxxz = pbuffer.data(idx_kin_gh + 107);

    auto tk_xxzz_xxxyy = pbuffer.data(idx_kin_gh + 108);

    auto tk_xxzz_xxxyz = pbuffer.data(idx_kin_gh + 109);

    auto tk_xxzz_xxxzz = pbuffer.data(idx_kin_gh + 110);

    auto tk_xxzz_xxyyy = pbuffer.data(idx_kin_gh + 111);

    auto tk_xxzz_xxyyz = pbuffer.data(idx_kin_gh + 112);

    auto tk_xxzz_xxyzz = pbuffer.data(idx_kin_gh + 113);

    auto tk_xxzz_xxzzz = pbuffer.data(idx_kin_gh + 114);

    auto tk_xxzz_xyyyy = pbuffer.data(idx_kin_gh + 115);

    auto tk_xxzz_xyyyz = pbuffer.data(idx_kin_gh + 116);

    auto tk_xxzz_xyyzz = pbuffer.data(idx_kin_gh + 117);

    auto tk_xxzz_xyzzz = pbuffer.data(idx_kin_gh + 118);

    auto tk_xxzz_xzzzz = pbuffer.data(idx_kin_gh + 119);

    auto tk_xxzz_yyyyy = pbuffer.data(idx_kin_gh + 120);

    auto tk_xxzz_yyyyz = pbuffer.data(idx_kin_gh + 121);

    auto tk_xxzz_yyyzz = pbuffer.data(idx_kin_gh + 122);

    auto tk_xxzz_yyzzz = pbuffer.data(idx_kin_gh + 123);

    auto tk_xxzz_yzzzz = pbuffer.data(idx_kin_gh + 124);

    auto tk_xxzz_zzzzz = pbuffer.data(idx_kin_gh + 125);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xx_xxxxx,   \
                             tk_xx_xxxxy,   \
                             tk_xx_xxxyy,   \
                             tk_xx_xxyyy,   \
                             tk_xx_xyyyy,   \
                             tk_xxz_xxxxx,  \
                             tk_xxz_xxxxy,  \
                             tk_xxz_xxxyy,  \
                             tk_xxz_xxyyy,  \
                             tk_xxz_xyyyy,  \
                             tk_xxzz_xxxxx, \
                             tk_xxzz_xxxxy, \
                             tk_xxzz_xxxxz, \
                             tk_xxzz_xxxyy, \
                             tk_xxzz_xxxyz, \
                             tk_xxzz_xxxzz, \
                             tk_xxzz_xxyyy, \
                             tk_xxzz_xxyyz, \
                             tk_xxzz_xxyzz, \
                             tk_xxzz_xxzzz, \
                             tk_xxzz_xyyyy, \
                             tk_xxzz_xyyyz, \
                             tk_xxzz_xyyzz, \
                             tk_xxzz_xyzzz, \
                             tk_xxzz_xzzzz, \
                             tk_xxzz_yyyyy, \
                             tk_xxzz_yyyyz, \
                             tk_xxzz_yyyzz, \
                             tk_xxzz_yyzzz, \
                             tk_xxzz_yzzzz, \
                             tk_xxzz_zzzzz, \
                             tk_xzz_xxxxz,  \
                             tk_xzz_xxxyz,  \
                             tk_xzz_xxxz,   \
                             tk_xzz_xxxzz,  \
                             tk_xzz_xxyyz,  \
                             tk_xzz_xxyz,   \
                             tk_xzz_xxyzz,  \
                             tk_xzz_xxzz,   \
                             tk_xzz_xxzzz,  \
                             tk_xzz_xyyyz,  \
                             tk_xzz_xyyz,   \
                             tk_xzz_xyyzz,  \
                             tk_xzz_xyzz,   \
                             tk_xzz_xyzzz,  \
                             tk_xzz_xzzz,   \
                             tk_xzz_xzzzz,  \
                             tk_xzz_yyyyy,  \
                             tk_xzz_yyyyz,  \
                             tk_xzz_yyyz,   \
                             tk_xzz_yyyzz,  \
                             tk_xzz_yyzz,   \
                             tk_xzz_yyzzz,  \
                             tk_xzz_yzzz,   \
                             tk_xzz_yzzzz,  \
                             tk_xzz_zzzz,   \
                             tk_xzz_zzzzz,  \
                             tk_zz_xxxxz,   \
                             tk_zz_xxxyz,   \
                             tk_zz_xxxzz,   \
                             tk_zz_xxyyz,   \
                             tk_zz_xxyzz,   \
                             tk_zz_xxzzz,   \
                             tk_zz_xyyyz,   \
                             tk_zz_xyyzz,   \
                             tk_zz_xyzzz,   \
                             tk_zz_xzzzz,   \
                             tk_zz_yyyyy,   \
                             tk_zz_yyyyz,   \
                             tk_zz_yyyzz,   \
                             tk_zz_yyzzz,   \
                             tk_zz_yzzzz,   \
                             tk_zz_zzzzz,   \
                             ts_xx_xxxxx,   \
                             ts_xx_xxxxy,   \
                             ts_xx_xxxyy,   \
                             ts_xx_xxyyy,   \
                             ts_xx_xyyyy,   \
                             ts_xxzz_xxxxx, \
                             ts_xxzz_xxxxy, \
                             ts_xxzz_xxxxz, \
                             ts_xxzz_xxxyy, \
                             ts_xxzz_xxxyz, \
                             ts_xxzz_xxxzz, \
                             ts_xxzz_xxyyy, \
                             ts_xxzz_xxyyz, \
                             ts_xxzz_xxyzz, \
                             ts_xxzz_xxzzz, \
                             ts_xxzz_xyyyy, \
                             ts_xxzz_xyyyz, \
                             ts_xxzz_xyyzz, \
                             ts_xxzz_xyzzz, \
                             ts_xxzz_xzzzz, \
                             ts_xxzz_yyyyy, \
                             ts_xxzz_yyyyz, \
                             ts_xxzz_yyyzz, \
                             ts_xxzz_yyzzz, \
                             ts_xxzz_yzzzz, \
                             ts_xxzz_zzzzz, \
                             ts_zz_xxxxz,   \
                             ts_zz_xxxyz,   \
                             ts_zz_xxxzz,   \
                             ts_zz_xxyyz,   \
                             ts_zz_xxyzz,   \
                             ts_zz_xxzzz,   \
                             ts_zz_xyyyz,   \
                             ts_zz_xyyzz,   \
                             ts_zz_xyzzz,   \
                             ts_zz_xzzzz,   \
                             ts_zz_yyyyy,   \
                             ts_zz_yyyyz,   \
                             ts_zz_yyyzz,   \
                             ts_zz_yyzzz,   \
                             ts_zz_yzzzz,   \
                             ts_zz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzz_xxxxx[i] = -2.0 * ts_xx_xxxxx[i] * fbe_0 * fz_0 + tk_xx_xxxxx[i] * fe_0 + tk_xxz_xxxxx[i] * pa_z[i] + 2.0 * ts_xxzz_xxxxx[i] * fz_0;

        tk_xxzz_xxxxy[i] = -2.0 * ts_xx_xxxxy[i] * fbe_0 * fz_0 + tk_xx_xxxxy[i] * fe_0 + tk_xxz_xxxxy[i] * pa_z[i] + 2.0 * ts_xxzz_xxxxy[i] * fz_0;

        tk_xxzz_xxxxz[i] = -2.0 * ts_zz_xxxxz[i] * fbe_0 * fz_0 + tk_zz_xxxxz[i] * fe_0 + 4.0 * tk_xzz_xxxz[i] * fe_0 + tk_xzz_xxxxz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xxxxz[i] * fz_0;

        tk_xxzz_xxxyy[i] = -2.0 * ts_xx_xxxyy[i] * fbe_0 * fz_0 + tk_xx_xxxyy[i] * fe_0 + tk_xxz_xxxyy[i] * pa_z[i] + 2.0 * ts_xxzz_xxxyy[i] * fz_0;

        tk_xxzz_xxxyz[i] = -2.0 * ts_zz_xxxyz[i] * fbe_0 * fz_0 + tk_zz_xxxyz[i] * fe_0 + 3.0 * tk_xzz_xxyz[i] * fe_0 + tk_xzz_xxxyz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xxxyz[i] * fz_0;

        tk_xxzz_xxxzz[i] = -2.0 * ts_zz_xxxzz[i] * fbe_0 * fz_0 + tk_zz_xxxzz[i] * fe_0 + 3.0 * tk_xzz_xxzz[i] * fe_0 + tk_xzz_xxxzz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xxxzz[i] * fz_0;

        tk_xxzz_xxyyy[i] = -2.0 * ts_xx_xxyyy[i] * fbe_0 * fz_0 + tk_xx_xxyyy[i] * fe_0 + tk_xxz_xxyyy[i] * pa_z[i] + 2.0 * ts_xxzz_xxyyy[i] * fz_0;

        tk_xxzz_xxyyz[i] = -2.0 * ts_zz_xxyyz[i] * fbe_0 * fz_0 + tk_zz_xxyyz[i] * fe_0 + 2.0 * tk_xzz_xyyz[i] * fe_0 + tk_xzz_xxyyz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xxyyz[i] * fz_0;

        tk_xxzz_xxyzz[i] = -2.0 * ts_zz_xxyzz[i] * fbe_0 * fz_0 + tk_zz_xxyzz[i] * fe_0 + 2.0 * tk_xzz_xyzz[i] * fe_0 + tk_xzz_xxyzz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xxyzz[i] * fz_0;

        tk_xxzz_xxzzz[i] = -2.0 * ts_zz_xxzzz[i] * fbe_0 * fz_0 + tk_zz_xxzzz[i] * fe_0 + 2.0 * tk_xzz_xzzz[i] * fe_0 + tk_xzz_xxzzz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xxzzz[i] * fz_0;

        tk_xxzz_xyyyy[i] = -2.0 * ts_xx_xyyyy[i] * fbe_0 * fz_0 + tk_xx_xyyyy[i] * fe_0 + tk_xxz_xyyyy[i] * pa_z[i] + 2.0 * ts_xxzz_xyyyy[i] * fz_0;

        tk_xxzz_xyyyz[i] = -2.0 * ts_zz_xyyyz[i] * fbe_0 * fz_0 + tk_zz_xyyyz[i] * fe_0 + tk_xzz_yyyz[i] * fe_0 + tk_xzz_xyyyz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xyyyz[i] * fz_0;

        tk_xxzz_xyyzz[i] = -2.0 * ts_zz_xyyzz[i] * fbe_0 * fz_0 + tk_zz_xyyzz[i] * fe_0 + tk_xzz_yyzz[i] * fe_0 + tk_xzz_xyyzz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xyyzz[i] * fz_0;

        tk_xxzz_xyzzz[i] = -2.0 * ts_zz_xyzzz[i] * fbe_0 * fz_0 + tk_zz_xyzzz[i] * fe_0 + tk_xzz_yzzz[i] * fe_0 + tk_xzz_xyzzz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xyzzz[i] * fz_0;

        tk_xxzz_xzzzz[i] = -2.0 * ts_zz_xzzzz[i] * fbe_0 * fz_0 + tk_zz_xzzzz[i] * fe_0 + tk_xzz_zzzz[i] * fe_0 + tk_xzz_xzzzz[i] * pa_x[i] +
                           2.0 * ts_xxzz_xzzzz[i] * fz_0;

        tk_xxzz_yyyyy[i] = -2.0 * ts_zz_yyyyy[i] * fbe_0 * fz_0 + tk_zz_yyyyy[i] * fe_0 + tk_xzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xxzz_yyyyy[i] * fz_0;

        tk_xxzz_yyyyz[i] = -2.0 * ts_zz_yyyyz[i] * fbe_0 * fz_0 + tk_zz_yyyyz[i] * fe_0 + tk_xzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xxzz_yyyyz[i] * fz_0;

        tk_xxzz_yyyzz[i] = -2.0 * ts_zz_yyyzz[i] * fbe_0 * fz_0 + tk_zz_yyyzz[i] * fe_0 + tk_xzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xxzz_yyyzz[i] * fz_0;

        tk_xxzz_yyzzz[i] = -2.0 * ts_zz_yyzzz[i] * fbe_0 * fz_0 + tk_zz_yyzzz[i] * fe_0 + tk_xzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xxzz_yyzzz[i] * fz_0;

        tk_xxzz_yzzzz[i] = -2.0 * ts_zz_yzzzz[i] * fbe_0 * fz_0 + tk_zz_yzzzz[i] * fe_0 + tk_xzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xxzz_yzzzz[i] * fz_0;

        tk_xxzz_zzzzz[i] = -2.0 * ts_zz_zzzzz[i] * fbe_0 * fz_0 + tk_zz_zzzzz[i] * fe_0 + tk_xzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xxzz_zzzzz[i] * fz_0;
    }

    // Set up 126-147 components of targeted buffer : GH

    auto tk_xyyy_xxxxx = pbuffer.data(idx_kin_gh + 126);

    auto tk_xyyy_xxxxy = pbuffer.data(idx_kin_gh + 127);

    auto tk_xyyy_xxxxz = pbuffer.data(idx_kin_gh + 128);

    auto tk_xyyy_xxxyy = pbuffer.data(idx_kin_gh + 129);

    auto tk_xyyy_xxxyz = pbuffer.data(idx_kin_gh + 130);

    auto tk_xyyy_xxxzz = pbuffer.data(idx_kin_gh + 131);

    auto tk_xyyy_xxyyy = pbuffer.data(idx_kin_gh + 132);

    auto tk_xyyy_xxyyz = pbuffer.data(idx_kin_gh + 133);

    auto tk_xyyy_xxyzz = pbuffer.data(idx_kin_gh + 134);

    auto tk_xyyy_xxzzz = pbuffer.data(idx_kin_gh + 135);

    auto tk_xyyy_xyyyy = pbuffer.data(idx_kin_gh + 136);

    auto tk_xyyy_xyyyz = pbuffer.data(idx_kin_gh + 137);

    auto tk_xyyy_xyyzz = pbuffer.data(idx_kin_gh + 138);

    auto tk_xyyy_xyzzz = pbuffer.data(idx_kin_gh + 139);

    auto tk_xyyy_xzzzz = pbuffer.data(idx_kin_gh + 140);

    auto tk_xyyy_yyyyy = pbuffer.data(idx_kin_gh + 141);

    auto tk_xyyy_yyyyz = pbuffer.data(idx_kin_gh + 142);

    auto tk_xyyy_yyyzz = pbuffer.data(idx_kin_gh + 143);

    auto tk_xyyy_yyzzz = pbuffer.data(idx_kin_gh + 144);

    auto tk_xyyy_yzzzz = pbuffer.data(idx_kin_gh + 145);

    auto tk_xyyy_zzzzz = pbuffer.data(idx_kin_gh + 146);

#pragma omp simd aligned(pa_x,              \
                             tk_xyyy_xxxxx, \
                             tk_xyyy_xxxxy, \
                             tk_xyyy_xxxxz, \
                             tk_xyyy_xxxyy, \
                             tk_xyyy_xxxyz, \
                             tk_xyyy_xxxzz, \
                             tk_xyyy_xxyyy, \
                             tk_xyyy_xxyyz, \
                             tk_xyyy_xxyzz, \
                             tk_xyyy_xxzzz, \
                             tk_xyyy_xyyyy, \
                             tk_xyyy_xyyyz, \
                             tk_xyyy_xyyzz, \
                             tk_xyyy_xyzzz, \
                             tk_xyyy_xzzzz, \
                             tk_xyyy_yyyyy, \
                             tk_xyyy_yyyyz, \
                             tk_xyyy_yyyzz, \
                             tk_xyyy_yyzzz, \
                             tk_xyyy_yzzzz, \
                             tk_xyyy_zzzzz, \
                             tk_yyy_xxxx,   \
                             tk_yyy_xxxxx,  \
                             tk_yyy_xxxxy,  \
                             tk_yyy_xxxxz,  \
                             tk_yyy_xxxy,   \
                             tk_yyy_xxxyy,  \
                             tk_yyy_xxxyz,  \
                             tk_yyy_xxxz,   \
                             tk_yyy_xxxzz,  \
                             tk_yyy_xxyy,   \
                             tk_yyy_xxyyy,  \
                             tk_yyy_xxyyz,  \
                             tk_yyy_xxyz,   \
                             tk_yyy_xxyzz,  \
                             tk_yyy_xxzz,   \
                             tk_yyy_xxzzz,  \
                             tk_yyy_xyyy,   \
                             tk_yyy_xyyyy,  \
                             tk_yyy_xyyyz,  \
                             tk_yyy_xyyz,   \
                             tk_yyy_xyyzz,  \
                             tk_yyy_xyzz,   \
                             tk_yyy_xyzzz,  \
                             tk_yyy_xzzz,   \
                             tk_yyy_xzzzz,  \
                             tk_yyy_yyyy,   \
                             tk_yyy_yyyyy,  \
                             tk_yyy_yyyyz,  \
                             tk_yyy_yyyz,   \
                             tk_yyy_yyyzz,  \
                             tk_yyy_yyzz,   \
                             tk_yyy_yyzzz,  \
                             tk_yyy_yzzz,   \
                             tk_yyy_yzzzz,  \
                             tk_yyy_zzzz,   \
                             tk_yyy_zzzzz,  \
                             ts_xyyy_xxxxx, \
                             ts_xyyy_xxxxy, \
                             ts_xyyy_xxxxz, \
                             ts_xyyy_xxxyy, \
                             ts_xyyy_xxxyz, \
                             ts_xyyy_xxxzz, \
                             ts_xyyy_xxyyy, \
                             ts_xyyy_xxyyz, \
                             ts_xyyy_xxyzz, \
                             ts_xyyy_xxzzz, \
                             ts_xyyy_xyyyy, \
                             ts_xyyy_xyyyz, \
                             ts_xyyy_xyyzz, \
                             ts_xyyy_xyzzz, \
                             ts_xyyy_xzzzz, \
                             ts_xyyy_yyyyy, \
                             ts_xyyy_yyyyz, \
                             ts_xyyy_yyyzz, \
                             ts_xyyy_yyzzz, \
                             ts_xyyy_yzzzz, \
                             ts_xyyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyy_xxxxx[i] = 5.0 * tk_yyy_xxxx[i] * fe_0 + tk_yyy_xxxxx[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxx[i] * fz_0;

        tk_xyyy_xxxxy[i] = 4.0 * tk_yyy_xxxy[i] * fe_0 + tk_yyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxy[i] * fz_0;

        tk_xyyy_xxxxz[i] = 4.0 * tk_yyy_xxxz[i] * fe_0 + tk_yyy_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxz[i] * fz_0;

        tk_xyyy_xxxyy[i] = 3.0 * tk_yyy_xxyy[i] * fe_0 + tk_yyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xyyy_xxxyy[i] * fz_0;

        tk_xyyy_xxxyz[i] = 3.0 * tk_yyy_xxyz[i] * fe_0 + tk_yyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxyz[i] * fz_0;

        tk_xyyy_xxxzz[i] = 3.0 * tk_yyy_xxzz[i] * fe_0 + tk_yyy_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxzz[i] * fz_0;

        tk_xyyy_xxyyy[i] = 2.0 * tk_yyy_xyyy[i] * fe_0 + tk_yyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xyyy_xxyyy[i] * fz_0;

        tk_xyyy_xxyyz[i] = 2.0 * tk_yyy_xyyz[i] * fe_0 + tk_yyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyy_xxyyz[i] * fz_0;

        tk_xyyy_xxyzz[i] = 2.0 * tk_yyy_xyzz[i] * fe_0 + tk_yyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxyzz[i] * fz_0;

        tk_xyyy_xxzzz[i] = 2.0 * tk_yyy_xzzz[i] * fe_0 + tk_yyy_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxzzz[i] * fz_0;

        tk_xyyy_xyyyy[i] = tk_yyy_yyyy[i] * fe_0 + tk_yyy_xyyyy[i] * pa_x[i] + 2.0 * ts_xyyy_xyyyy[i] * fz_0;

        tk_xyyy_xyyyz[i] = tk_yyy_yyyz[i] * fe_0 + tk_yyy_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyy_xyyyz[i] * fz_0;

        tk_xyyy_xyyzz[i] = tk_yyy_yyzz[i] * fe_0 + tk_yyy_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyy_xyyzz[i] * fz_0;

        tk_xyyy_xyzzz[i] = tk_yyy_yzzz[i] * fe_0 + tk_yyy_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xyzzz[i] * fz_0;

        tk_xyyy_xzzzz[i] = tk_yyy_zzzz[i] * fe_0 + tk_yyy_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xzzzz[i] * fz_0;

        tk_xyyy_yyyyy[i] = tk_yyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyy_yyyyy[i] * fz_0;

        tk_xyyy_yyyyz[i] = tk_yyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyy_yyyyz[i] * fz_0;

        tk_xyyy_yyyzz[i] = tk_yyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyy_yyyzz[i] * fz_0;

        tk_xyyy_yyzzz[i] = tk_yyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyy_yyzzz[i] * fz_0;

        tk_xyyy_yzzzz[i] = tk_yyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_yzzzz[i] * fz_0;

        tk_xyyy_zzzzz[i] = tk_yyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_zzzzz[i] * fz_0;
    }

    // Set up 147-168 components of targeted buffer : GH

    auto tk_xyyz_xxxxx = pbuffer.data(idx_kin_gh + 147);

    auto tk_xyyz_xxxxy = pbuffer.data(idx_kin_gh + 148);

    auto tk_xyyz_xxxxz = pbuffer.data(idx_kin_gh + 149);

    auto tk_xyyz_xxxyy = pbuffer.data(idx_kin_gh + 150);

    auto tk_xyyz_xxxyz = pbuffer.data(idx_kin_gh + 151);

    auto tk_xyyz_xxxzz = pbuffer.data(idx_kin_gh + 152);

    auto tk_xyyz_xxyyy = pbuffer.data(idx_kin_gh + 153);

    auto tk_xyyz_xxyyz = pbuffer.data(idx_kin_gh + 154);

    auto tk_xyyz_xxyzz = pbuffer.data(idx_kin_gh + 155);

    auto tk_xyyz_xxzzz = pbuffer.data(idx_kin_gh + 156);

    auto tk_xyyz_xyyyy = pbuffer.data(idx_kin_gh + 157);

    auto tk_xyyz_xyyyz = pbuffer.data(idx_kin_gh + 158);

    auto tk_xyyz_xyyzz = pbuffer.data(idx_kin_gh + 159);

    auto tk_xyyz_xyzzz = pbuffer.data(idx_kin_gh + 160);

    auto tk_xyyz_xzzzz = pbuffer.data(idx_kin_gh + 161);

    auto tk_xyyz_yyyyy = pbuffer.data(idx_kin_gh + 162);

    auto tk_xyyz_yyyyz = pbuffer.data(idx_kin_gh + 163);

    auto tk_xyyz_yyyzz = pbuffer.data(idx_kin_gh + 164);

    auto tk_xyyz_yyzzz = pbuffer.data(idx_kin_gh + 165);

    auto tk_xyyz_yzzzz = pbuffer.data(idx_kin_gh + 166);

    auto tk_xyyz_zzzzz = pbuffer.data(idx_kin_gh + 167);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tk_xyy_xxxxx,  \
                             tk_xyy_xxxxy,  \
                             tk_xyy_xxxyy,  \
                             tk_xyy_xxyyy,  \
                             tk_xyy_xyyyy,  \
                             tk_xyyz_xxxxx, \
                             tk_xyyz_xxxxy, \
                             tk_xyyz_xxxxz, \
                             tk_xyyz_xxxyy, \
                             tk_xyyz_xxxyz, \
                             tk_xyyz_xxxzz, \
                             tk_xyyz_xxyyy, \
                             tk_xyyz_xxyyz, \
                             tk_xyyz_xxyzz, \
                             tk_xyyz_xxzzz, \
                             tk_xyyz_xyyyy, \
                             tk_xyyz_xyyyz, \
                             tk_xyyz_xyyzz, \
                             tk_xyyz_xyzzz, \
                             tk_xyyz_xzzzz, \
                             tk_xyyz_yyyyy, \
                             tk_xyyz_yyyyz, \
                             tk_xyyz_yyyzz, \
                             tk_xyyz_yyzzz, \
                             tk_xyyz_yzzzz, \
                             tk_xyyz_zzzzz, \
                             tk_yyz_xxxxz,  \
                             tk_yyz_xxxyz,  \
                             tk_yyz_xxxz,   \
                             tk_yyz_xxxzz,  \
                             tk_yyz_xxyyz,  \
                             tk_yyz_xxyz,   \
                             tk_yyz_xxyzz,  \
                             tk_yyz_xxzz,   \
                             tk_yyz_xxzzz,  \
                             tk_yyz_xyyyz,  \
                             tk_yyz_xyyz,   \
                             tk_yyz_xyyzz,  \
                             tk_yyz_xyzz,   \
                             tk_yyz_xyzzz,  \
                             tk_yyz_xzzz,   \
                             tk_yyz_xzzzz,  \
                             tk_yyz_yyyyy,  \
                             tk_yyz_yyyyz,  \
                             tk_yyz_yyyz,   \
                             tk_yyz_yyyzz,  \
                             tk_yyz_yyzz,   \
                             tk_yyz_yyzzz,  \
                             tk_yyz_yzzz,   \
                             tk_yyz_yzzzz,  \
                             tk_yyz_zzzz,   \
                             tk_yyz_zzzzz,  \
                             ts_xyyz_xxxxx, \
                             ts_xyyz_xxxxy, \
                             ts_xyyz_xxxxz, \
                             ts_xyyz_xxxyy, \
                             ts_xyyz_xxxyz, \
                             ts_xyyz_xxxzz, \
                             ts_xyyz_xxyyy, \
                             ts_xyyz_xxyyz, \
                             ts_xyyz_xxyzz, \
                             ts_xyyz_xxzzz, \
                             ts_xyyz_xyyyy, \
                             ts_xyyz_xyyyz, \
                             ts_xyyz_xyyzz, \
                             ts_xyyz_xyzzz, \
                             ts_xyyz_xzzzz, \
                             ts_xyyz_yyyyy, \
                             ts_xyyz_yyyyz, \
                             ts_xyyz_yyyzz, \
                             ts_xyyz_yyzzz, \
                             ts_xyyz_yzzzz, \
                             ts_xyyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyz_xxxxx[i] = tk_xyy_xxxxx[i] * pa_z[i] + 2.0 * ts_xyyz_xxxxx[i] * fz_0;

        tk_xyyz_xxxxy[i] = tk_xyy_xxxxy[i] * pa_z[i] + 2.0 * ts_xyyz_xxxxy[i] * fz_0;

        tk_xyyz_xxxxz[i] = 4.0 * tk_yyz_xxxz[i] * fe_0 + tk_yyz_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxxz[i] * fz_0;

        tk_xyyz_xxxyy[i] = tk_xyy_xxxyy[i] * pa_z[i] + 2.0 * ts_xyyz_xxxyy[i] * fz_0;

        tk_xyyz_xxxyz[i] = 3.0 * tk_yyz_xxyz[i] * fe_0 + tk_yyz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxyz[i] * fz_0;

        tk_xyyz_xxxzz[i] = 3.0 * tk_yyz_xxzz[i] * fe_0 + tk_yyz_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxzz[i] * fz_0;

        tk_xyyz_xxyyy[i] = tk_xyy_xxyyy[i] * pa_z[i] + 2.0 * ts_xyyz_xxyyy[i] * fz_0;

        tk_xyyz_xxyyz[i] = 2.0 * tk_yyz_xyyz[i] * fe_0 + tk_yyz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyz_xxyyz[i] * fz_0;

        tk_xyyz_xxyzz[i] = 2.0 * tk_yyz_xyzz[i] * fe_0 + tk_yyz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxyzz[i] * fz_0;

        tk_xyyz_xxzzz[i] = 2.0 * tk_yyz_xzzz[i] * fe_0 + tk_yyz_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxzzz[i] * fz_0;

        tk_xyyz_xyyyy[i] = tk_xyy_xyyyy[i] * pa_z[i] + 2.0 * ts_xyyz_xyyyy[i] * fz_0;

        tk_xyyz_xyyyz[i] = tk_yyz_yyyz[i] * fe_0 + tk_yyz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyz_xyyyz[i] * fz_0;

        tk_xyyz_xyyzz[i] = tk_yyz_yyzz[i] * fe_0 + tk_yyz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyz_xyyzz[i] * fz_0;

        tk_xyyz_xyzzz[i] = tk_yyz_yzzz[i] * fe_0 + tk_yyz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xyzzz[i] * fz_0;

        tk_xyyz_xzzzz[i] = tk_yyz_zzzz[i] * fe_0 + tk_yyz_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xzzzz[i] * fz_0;

        tk_xyyz_yyyyy[i] = tk_yyz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyz_yyyyy[i] * fz_0;

        tk_xyyz_yyyyz[i] = tk_yyz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyz_yyyyz[i] * fz_0;

        tk_xyyz_yyyzz[i] = tk_yyz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyz_yyyzz[i] * fz_0;

        tk_xyyz_yyzzz[i] = tk_yyz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyz_yyzzz[i] * fz_0;

        tk_xyyz_yzzzz[i] = tk_yyz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_yzzzz[i] * fz_0;

        tk_xyyz_zzzzz[i] = tk_yyz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_zzzzz[i] * fz_0;
    }

    // Set up 168-189 components of targeted buffer : GH

    auto tk_xyzz_xxxxx = pbuffer.data(idx_kin_gh + 168);

    auto tk_xyzz_xxxxy = pbuffer.data(idx_kin_gh + 169);

    auto tk_xyzz_xxxxz = pbuffer.data(idx_kin_gh + 170);

    auto tk_xyzz_xxxyy = pbuffer.data(idx_kin_gh + 171);

    auto tk_xyzz_xxxyz = pbuffer.data(idx_kin_gh + 172);

    auto tk_xyzz_xxxzz = pbuffer.data(idx_kin_gh + 173);

    auto tk_xyzz_xxyyy = pbuffer.data(idx_kin_gh + 174);

    auto tk_xyzz_xxyyz = pbuffer.data(idx_kin_gh + 175);

    auto tk_xyzz_xxyzz = pbuffer.data(idx_kin_gh + 176);

    auto tk_xyzz_xxzzz = pbuffer.data(idx_kin_gh + 177);

    auto tk_xyzz_xyyyy = pbuffer.data(idx_kin_gh + 178);

    auto tk_xyzz_xyyyz = pbuffer.data(idx_kin_gh + 179);

    auto tk_xyzz_xyyzz = pbuffer.data(idx_kin_gh + 180);

    auto tk_xyzz_xyzzz = pbuffer.data(idx_kin_gh + 181);

    auto tk_xyzz_xzzzz = pbuffer.data(idx_kin_gh + 182);

    auto tk_xyzz_yyyyy = pbuffer.data(idx_kin_gh + 183);

    auto tk_xyzz_yyyyz = pbuffer.data(idx_kin_gh + 184);

    auto tk_xyzz_yyyzz = pbuffer.data(idx_kin_gh + 185);

    auto tk_xyzz_yyzzz = pbuffer.data(idx_kin_gh + 186);

    auto tk_xyzz_yzzzz = pbuffer.data(idx_kin_gh + 187);

    auto tk_xyzz_zzzzz = pbuffer.data(idx_kin_gh + 188);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tk_xyzz_xxxxx, \
                             tk_xyzz_xxxxy, \
                             tk_xyzz_xxxxz, \
                             tk_xyzz_xxxyy, \
                             tk_xyzz_xxxyz, \
                             tk_xyzz_xxxzz, \
                             tk_xyzz_xxyyy, \
                             tk_xyzz_xxyyz, \
                             tk_xyzz_xxyzz, \
                             tk_xyzz_xxzzz, \
                             tk_xyzz_xyyyy, \
                             tk_xyzz_xyyyz, \
                             tk_xyzz_xyyzz, \
                             tk_xyzz_xyzzz, \
                             tk_xyzz_xzzzz, \
                             tk_xyzz_yyyyy, \
                             tk_xyzz_yyyyz, \
                             tk_xyzz_yyyzz, \
                             tk_xyzz_yyzzz, \
                             tk_xyzz_yzzzz, \
                             tk_xyzz_zzzzz, \
                             tk_xzz_xxxxx,  \
                             tk_xzz_xxxxz,  \
                             tk_xzz_xxxzz,  \
                             tk_xzz_xxzzz,  \
                             tk_xzz_xzzzz,  \
                             tk_yzz_xxxxy,  \
                             tk_yzz_xxxy,   \
                             tk_yzz_xxxyy,  \
                             tk_yzz_xxxyz,  \
                             tk_yzz_xxyy,   \
                             tk_yzz_xxyyy,  \
                             tk_yzz_xxyyz,  \
                             tk_yzz_xxyz,   \
                             tk_yzz_xxyzz,  \
                             tk_yzz_xyyy,   \
                             tk_yzz_xyyyy,  \
                             tk_yzz_xyyyz,  \
                             tk_yzz_xyyz,   \
                             tk_yzz_xyyzz,  \
                             tk_yzz_xyzz,   \
                             tk_yzz_xyzzz,  \
                             tk_yzz_yyyy,   \
                             tk_yzz_yyyyy,  \
                             tk_yzz_yyyyz,  \
                             tk_yzz_yyyz,   \
                             tk_yzz_yyyzz,  \
                             tk_yzz_yyzz,   \
                             tk_yzz_yyzzz,  \
                             tk_yzz_yzzz,   \
                             tk_yzz_yzzzz,  \
                             tk_yzz_zzzzz,  \
                             ts_xyzz_xxxxx, \
                             ts_xyzz_xxxxy, \
                             ts_xyzz_xxxxz, \
                             ts_xyzz_xxxyy, \
                             ts_xyzz_xxxyz, \
                             ts_xyzz_xxxzz, \
                             ts_xyzz_xxyyy, \
                             ts_xyzz_xxyyz, \
                             ts_xyzz_xxyzz, \
                             ts_xyzz_xxzzz, \
                             ts_xyzz_xyyyy, \
                             ts_xyzz_xyyyz, \
                             ts_xyzz_xyyzz, \
                             ts_xyzz_xyzzz, \
                             ts_xyzz_xzzzz, \
                             ts_xyzz_yyyyy, \
                             ts_xyzz_yyyyz, \
                             ts_xyzz_yyyzz, \
                             ts_xyzz_yyzzz, \
                             ts_xyzz_yzzzz, \
                             ts_xyzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzz_xxxxx[i] = tk_xzz_xxxxx[i] * pa_y[i] + 2.0 * ts_xyzz_xxxxx[i] * fz_0;

        tk_xyzz_xxxxy[i] = 4.0 * tk_yzz_xxxy[i] * fe_0 + tk_yzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xyzz_xxxxy[i] * fz_0;

        tk_xyzz_xxxxz[i] = tk_xzz_xxxxz[i] * pa_y[i] + 2.0 * ts_xyzz_xxxxz[i] * fz_0;

        tk_xyzz_xxxyy[i] = 3.0 * tk_yzz_xxyy[i] * fe_0 + tk_yzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xyzz_xxxyy[i] * fz_0;

        tk_xyzz_xxxyz[i] = 3.0 * tk_yzz_xxyz[i] * fe_0 + tk_yzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyzz_xxxyz[i] * fz_0;

        tk_xyzz_xxxzz[i] = tk_xzz_xxxzz[i] * pa_y[i] + 2.0 * ts_xyzz_xxxzz[i] * fz_0;

        tk_xyzz_xxyyy[i] = 2.0 * tk_yzz_xyyy[i] * fe_0 + tk_yzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xyzz_xxyyy[i] * fz_0;

        tk_xyzz_xxyyz[i] = 2.0 * tk_yzz_xyyz[i] * fe_0 + tk_yzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyzz_xxyyz[i] * fz_0;

        tk_xyzz_xxyzz[i] = 2.0 * tk_yzz_xyzz[i] * fe_0 + tk_yzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyzz_xxyzz[i] * fz_0;

        tk_xyzz_xxzzz[i] = tk_xzz_xxzzz[i] * pa_y[i] + 2.0 * ts_xyzz_xxzzz[i] * fz_0;

        tk_xyzz_xyyyy[i] = tk_yzz_yyyy[i] * fe_0 + tk_yzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xyzz_xyyyy[i] * fz_0;

        tk_xyzz_xyyyz[i] = tk_yzz_yyyz[i] * fe_0 + tk_yzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyzz_xyyyz[i] * fz_0;

        tk_xyzz_xyyzz[i] = tk_yzz_yyzz[i] * fe_0 + tk_yzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyzz_xyyzz[i] * fz_0;

        tk_xyzz_xyzzz[i] = tk_yzz_yzzz[i] * fe_0 + tk_yzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyzz_xyzzz[i] * fz_0;

        tk_xyzz_xzzzz[i] = tk_xzz_xzzzz[i] * pa_y[i] + 2.0 * ts_xyzz_xzzzz[i] * fz_0;

        tk_xyzz_yyyyy[i] = tk_yzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyzz_yyyyy[i] * fz_0;

        tk_xyzz_yyyyz[i] = tk_yzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyzz_yyyyz[i] * fz_0;

        tk_xyzz_yyyzz[i] = tk_yzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyzz_yyyzz[i] * fz_0;

        tk_xyzz_yyzzz[i] = tk_yzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyzz_yyzzz[i] * fz_0;

        tk_xyzz_yzzzz[i] = tk_yzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyzz_yzzzz[i] * fz_0;

        tk_xyzz_zzzzz[i] = tk_yzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyzz_zzzzz[i] * fz_0;
    }

    // Set up 189-210 components of targeted buffer : GH

    auto tk_xzzz_xxxxx = pbuffer.data(idx_kin_gh + 189);

    auto tk_xzzz_xxxxy = pbuffer.data(idx_kin_gh + 190);

    auto tk_xzzz_xxxxz = pbuffer.data(idx_kin_gh + 191);

    auto tk_xzzz_xxxyy = pbuffer.data(idx_kin_gh + 192);

    auto tk_xzzz_xxxyz = pbuffer.data(idx_kin_gh + 193);

    auto tk_xzzz_xxxzz = pbuffer.data(idx_kin_gh + 194);

    auto tk_xzzz_xxyyy = pbuffer.data(idx_kin_gh + 195);

    auto tk_xzzz_xxyyz = pbuffer.data(idx_kin_gh + 196);

    auto tk_xzzz_xxyzz = pbuffer.data(idx_kin_gh + 197);

    auto tk_xzzz_xxzzz = pbuffer.data(idx_kin_gh + 198);

    auto tk_xzzz_xyyyy = pbuffer.data(idx_kin_gh + 199);

    auto tk_xzzz_xyyyz = pbuffer.data(idx_kin_gh + 200);

    auto tk_xzzz_xyyzz = pbuffer.data(idx_kin_gh + 201);

    auto tk_xzzz_xyzzz = pbuffer.data(idx_kin_gh + 202);

    auto tk_xzzz_xzzzz = pbuffer.data(idx_kin_gh + 203);

    auto tk_xzzz_yyyyy = pbuffer.data(idx_kin_gh + 204);

    auto tk_xzzz_yyyyz = pbuffer.data(idx_kin_gh + 205);

    auto tk_xzzz_yyyzz = pbuffer.data(idx_kin_gh + 206);

    auto tk_xzzz_yyzzz = pbuffer.data(idx_kin_gh + 207);

    auto tk_xzzz_yzzzz = pbuffer.data(idx_kin_gh + 208);

    auto tk_xzzz_zzzzz = pbuffer.data(idx_kin_gh + 209);

#pragma omp simd aligned(pa_x,              \
                             tk_xzzz_xxxxx, \
                             tk_xzzz_xxxxy, \
                             tk_xzzz_xxxxz, \
                             tk_xzzz_xxxyy, \
                             tk_xzzz_xxxyz, \
                             tk_xzzz_xxxzz, \
                             tk_xzzz_xxyyy, \
                             tk_xzzz_xxyyz, \
                             tk_xzzz_xxyzz, \
                             tk_xzzz_xxzzz, \
                             tk_xzzz_xyyyy, \
                             tk_xzzz_xyyyz, \
                             tk_xzzz_xyyzz, \
                             tk_xzzz_xyzzz, \
                             tk_xzzz_xzzzz, \
                             tk_xzzz_yyyyy, \
                             tk_xzzz_yyyyz, \
                             tk_xzzz_yyyzz, \
                             tk_xzzz_yyzzz, \
                             tk_xzzz_yzzzz, \
                             tk_xzzz_zzzzz, \
                             tk_zzz_xxxx,   \
                             tk_zzz_xxxxx,  \
                             tk_zzz_xxxxy,  \
                             tk_zzz_xxxxz,  \
                             tk_zzz_xxxy,   \
                             tk_zzz_xxxyy,  \
                             tk_zzz_xxxyz,  \
                             tk_zzz_xxxz,   \
                             tk_zzz_xxxzz,  \
                             tk_zzz_xxyy,   \
                             tk_zzz_xxyyy,  \
                             tk_zzz_xxyyz,  \
                             tk_zzz_xxyz,   \
                             tk_zzz_xxyzz,  \
                             tk_zzz_xxzz,   \
                             tk_zzz_xxzzz,  \
                             tk_zzz_xyyy,   \
                             tk_zzz_xyyyy,  \
                             tk_zzz_xyyyz,  \
                             tk_zzz_xyyz,   \
                             tk_zzz_xyyzz,  \
                             tk_zzz_xyzz,   \
                             tk_zzz_xyzzz,  \
                             tk_zzz_xzzz,   \
                             tk_zzz_xzzzz,  \
                             tk_zzz_yyyy,   \
                             tk_zzz_yyyyy,  \
                             tk_zzz_yyyyz,  \
                             tk_zzz_yyyz,   \
                             tk_zzz_yyyzz,  \
                             tk_zzz_yyzz,   \
                             tk_zzz_yyzzz,  \
                             tk_zzz_yzzz,   \
                             tk_zzz_yzzzz,  \
                             tk_zzz_zzzz,   \
                             tk_zzz_zzzzz,  \
                             ts_xzzz_xxxxx, \
                             ts_xzzz_xxxxy, \
                             ts_xzzz_xxxxz, \
                             ts_xzzz_xxxyy, \
                             ts_xzzz_xxxyz, \
                             ts_xzzz_xxxzz, \
                             ts_xzzz_xxyyy, \
                             ts_xzzz_xxyyz, \
                             ts_xzzz_xxyzz, \
                             ts_xzzz_xxzzz, \
                             ts_xzzz_xyyyy, \
                             ts_xzzz_xyyyz, \
                             ts_xzzz_xyyzz, \
                             ts_xzzz_xyzzz, \
                             ts_xzzz_xzzzz, \
                             ts_xzzz_yyyyy, \
                             ts_xzzz_yyyyz, \
                             ts_xzzz_yyyzz, \
                             ts_xzzz_yyzzz, \
                             ts_xzzz_yzzzz, \
                             ts_xzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzz_xxxxx[i] = 5.0 * tk_zzz_xxxx[i] * fe_0 + tk_zzz_xxxxx[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxx[i] * fz_0;

        tk_xzzz_xxxxy[i] = 4.0 * tk_zzz_xxxy[i] * fe_0 + tk_zzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxy[i] * fz_0;

        tk_xzzz_xxxxz[i] = 4.0 * tk_zzz_xxxz[i] * fe_0 + tk_zzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxz[i] * fz_0;

        tk_xzzz_xxxyy[i] = 3.0 * tk_zzz_xxyy[i] * fe_0 + tk_zzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xzzz_xxxyy[i] * fz_0;

        tk_xzzz_xxxyz[i] = 3.0 * tk_zzz_xxyz[i] * fe_0 + tk_zzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxyz[i] * fz_0;

        tk_xzzz_xxxzz[i] = 3.0 * tk_zzz_xxzz[i] * fe_0 + tk_zzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxzz[i] * fz_0;

        tk_xzzz_xxyyy[i] = 2.0 * tk_zzz_xyyy[i] * fe_0 + tk_zzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xzzz_xxyyy[i] * fz_0;

        tk_xzzz_xxyyz[i] = 2.0 * tk_zzz_xyyz[i] * fe_0 + tk_zzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xzzz_xxyyz[i] * fz_0;

        tk_xzzz_xxyzz[i] = 2.0 * tk_zzz_xyzz[i] * fe_0 + tk_zzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxyzz[i] * fz_0;

        tk_xzzz_xxzzz[i] = 2.0 * tk_zzz_xzzz[i] * fe_0 + tk_zzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxzzz[i] * fz_0;

        tk_xzzz_xyyyy[i] = tk_zzz_yyyy[i] * fe_0 + tk_zzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xzzz_xyyyy[i] * fz_0;

        tk_xzzz_xyyyz[i] = tk_zzz_yyyz[i] * fe_0 + tk_zzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xzzz_xyyyz[i] * fz_0;

        tk_xzzz_xyyzz[i] = tk_zzz_yyzz[i] * fe_0 + tk_zzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xzzz_xyyzz[i] * fz_0;

        tk_xzzz_xyzzz[i] = tk_zzz_yzzz[i] * fe_0 + tk_zzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xyzzz[i] * fz_0;

        tk_xzzz_xzzzz[i] = tk_zzz_zzzz[i] * fe_0 + tk_zzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xzzzz[i] * fz_0;

        tk_xzzz_yyyyy[i] = tk_zzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xzzz_yyyyy[i] * fz_0;

        tk_xzzz_yyyyz[i] = tk_zzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xzzz_yyyyz[i] * fz_0;

        tk_xzzz_yyyzz[i] = tk_zzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xzzz_yyyzz[i] * fz_0;

        tk_xzzz_yyzzz[i] = tk_zzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xzzz_yyzzz[i] * fz_0;

        tk_xzzz_yzzzz[i] = tk_zzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_yzzzz[i] * fz_0;

        tk_xzzz_zzzzz[i] = tk_zzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_zzzzz[i] * fz_0;
    }

    // Set up 210-231 components of targeted buffer : GH

    auto tk_yyyy_xxxxx = pbuffer.data(idx_kin_gh + 210);

    auto tk_yyyy_xxxxy = pbuffer.data(idx_kin_gh + 211);

    auto tk_yyyy_xxxxz = pbuffer.data(idx_kin_gh + 212);

    auto tk_yyyy_xxxyy = pbuffer.data(idx_kin_gh + 213);

    auto tk_yyyy_xxxyz = pbuffer.data(idx_kin_gh + 214);

    auto tk_yyyy_xxxzz = pbuffer.data(idx_kin_gh + 215);

    auto tk_yyyy_xxyyy = pbuffer.data(idx_kin_gh + 216);

    auto tk_yyyy_xxyyz = pbuffer.data(idx_kin_gh + 217);

    auto tk_yyyy_xxyzz = pbuffer.data(idx_kin_gh + 218);

    auto tk_yyyy_xxzzz = pbuffer.data(idx_kin_gh + 219);

    auto tk_yyyy_xyyyy = pbuffer.data(idx_kin_gh + 220);

    auto tk_yyyy_xyyyz = pbuffer.data(idx_kin_gh + 221);

    auto tk_yyyy_xyyzz = pbuffer.data(idx_kin_gh + 222);

    auto tk_yyyy_xyzzz = pbuffer.data(idx_kin_gh + 223);

    auto tk_yyyy_xzzzz = pbuffer.data(idx_kin_gh + 224);

    auto tk_yyyy_yyyyy = pbuffer.data(idx_kin_gh + 225);

    auto tk_yyyy_yyyyz = pbuffer.data(idx_kin_gh + 226);

    auto tk_yyyy_yyyzz = pbuffer.data(idx_kin_gh + 227);

    auto tk_yyyy_yyzzz = pbuffer.data(idx_kin_gh + 228);

    auto tk_yyyy_yzzzz = pbuffer.data(idx_kin_gh + 229);

    auto tk_yyyy_zzzzz = pbuffer.data(idx_kin_gh + 230);

#pragma omp simd aligned(pa_y,              \
                             tk_yy_xxxxx,   \
                             tk_yy_xxxxy,   \
                             tk_yy_xxxxz,   \
                             tk_yy_xxxyy,   \
                             tk_yy_xxxyz,   \
                             tk_yy_xxxzz,   \
                             tk_yy_xxyyy,   \
                             tk_yy_xxyyz,   \
                             tk_yy_xxyzz,   \
                             tk_yy_xxzzz,   \
                             tk_yy_xyyyy,   \
                             tk_yy_xyyyz,   \
                             tk_yy_xyyzz,   \
                             tk_yy_xyzzz,   \
                             tk_yy_xzzzz,   \
                             tk_yy_yyyyy,   \
                             tk_yy_yyyyz,   \
                             tk_yy_yyyzz,   \
                             tk_yy_yyzzz,   \
                             tk_yy_yzzzz,   \
                             tk_yy_zzzzz,   \
                             tk_yyy_xxxx,   \
                             tk_yyy_xxxxx,  \
                             tk_yyy_xxxxy,  \
                             tk_yyy_xxxxz,  \
                             tk_yyy_xxxy,   \
                             tk_yyy_xxxyy,  \
                             tk_yyy_xxxyz,  \
                             tk_yyy_xxxz,   \
                             tk_yyy_xxxzz,  \
                             tk_yyy_xxyy,   \
                             tk_yyy_xxyyy,  \
                             tk_yyy_xxyyz,  \
                             tk_yyy_xxyz,   \
                             tk_yyy_xxyzz,  \
                             tk_yyy_xxzz,   \
                             tk_yyy_xxzzz,  \
                             tk_yyy_xyyy,   \
                             tk_yyy_xyyyy,  \
                             tk_yyy_xyyyz,  \
                             tk_yyy_xyyz,   \
                             tk_yyy_xyyzz,  \
                             tk_yyy_xyzz,   \
                             tk_yyy_xyzzz,  \
                             tk_yyy_xzzz,   \
                             tk_yyy_xzzzz,  \
                             tk_yyy_yyyy,   \
                             tk_yyy_yyyyy,  \
                             tk_yyy_yyyyz,  \
                             tk_yyy_yyyz,   \
                             tk_yyy_yyyzz,  \
                             tk_yyy_yyzz,   \
                             tk_yyy_yyzzz,  \
                             tk_yyy_yzzz,   \
                             tk_yyy_yzzzz,  \
                             tk_yyy_zzzz,   \
                             tk_yyy_zzzzz,  \
                             tk_yyyy_xxxxx, \
                             tk_yyyy_xxxxy, \
                             tk_yyyy_xxxxz, \
                             tk_yyyy_xxxyy, \
                             tk_yyyy_xxxyz, \
                             tk_yyyy_xxxzz, \
                             tk_yyyy_xxyyy, \
                             tk_yyyy_xxyyz, \
                             tk_yyyy_xxyzz, \
                             tk_yyyy_xxzzz, \
                             tk_yyyy_xyyyy, \
                             tk_yyyy_xyyyz, \
                             tk_yyyy_xyyzz, \
                             tk_yyyy_xyzzz, \
                             tk_yyyy_xzzzz, \
                             tk_yyyy_yyyyy, \
                             tk_yyyy_yyyyz, \
                             tk_yyyy_yyyzz, \
                             tk_yyyy_yyzzz, \
                             tk_yyyy_yzzzz, \
                             tk_yyyy_zzzzz, \
                             ts_yy_xxxxx,   \
                             ts_yy_xxxxy,   \
                             ts_yy_xxxxz,   \
                             ts_yy_xxxyy,   \
                             ts_yy_xxxyz,   \
                             ts_yy_xxxzz,   \
                             ts_yy_xxyyy,   \
                             ts_yy_xxyyz,   \
                             ts_yy_xxyzz,   \
                             ts_yy_xxzzz,   \
                             ts_yy_xyyyy,   \
                             ts_yy_xyyyz,   \
                             ts_yy_xyyzz,   \
                             ts_yy_xyzzz,   \
                             ts_yy_xzzzz,   \
                             ts_yy_yyyyy,   \
                             ts_yy_yyyyz,   \
                             ts_yy_yyyzz,   \
                             ts_yy_yyzzz,   \
                             ts_yy_yzzzz,   \
                             ts_yy_zzzzz,   \
                             ts_yyyy_xxxxx, \
                             ts_yyyy_xxxxy, \
                             ts_yyyy_xxxxz, \
                             ts_yyyy_xxxyy, \
                             ts_yyyy_xxxyz, \
                             ts_yyyy_xxxzz, \
                             ts_yyyy_xxyyy, \
                             ts_yyyy_xxyyz, \
                             ts_yyyy_xxyzz, \
                             ts_yyyy_xxzzz, \
                             ts_yyyy_xyyyy, \
                             ts_yyyy_xyyyz, \
                             ts_yyyy_xyyzz, \
                             ts_yyyy_xyzzz, \
                             ts_yyyy_xzzzz, \
                             ts_yyyy_yyyyy, \
                             ts_yyyy_yyyyz, \
                             ts_yyyy_yyyzz, \
                             ts_yyyy_yyzzz, \
                             ts_yyyy_yzzzz, \
                             ts_yyyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyy_xxxxx[i] =
            -6.0 * ts_yy_xxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxx[i] * fe_0 + tk_yyy_xxxxx[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxx[i] * fz_0;

        tk_yyyy_xxxxy[i] = -6.0 * ts_yy_xxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxy[i] * fe_0 + tk_yyy_xxxx[i] * fe_0 + tk_yyy_xxxxy[i] * pa_y[i] +
                           2.0 * ts_yyyy_xxxxy[i] * fz_0;

        tk_yyyy_xxxxz[i] =
            -6.0 * ts_yy_xxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxz[i] * fe_0 + tk_yyy_xxxxz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxz[i] * fz_0;

        tk_yyyy_xxxyy[i] = -6.0 * ts_yy_xxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxyy[i] * fe_0 + 2.0 * tk_yyy_xxxy[i] * fe_0 +
                           tk_yyy_xxxyy[i] * pa_y[i] + 2.0 * ts_yyyy_xxxyy[i] * fz_0;

        tk_yyyy_xxxyz[i] = -6.0 * ts_yy_xxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxyz[i] * fe_0 + tk_yyy_xxxz[i] * fe_0 + tk_yyy_xxxyz[i] * pa_y[i] +
                           2.0 * ts_yyyy_xxxyz[i] * fz_0;

        tk_yyyy_xxxzz[i] =
            -6.0 * ts_yy_xxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxzz[i] * fe_0 + tk_yyy_xxxzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxzz[i] * fz_0;

        tk_yyyy_xxyyy[i] = -6.0 * ts_yy_xxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyyy[i] * fe_0 + 3.0 * tk_yyy_xxyy[i] * fe_0 +
                           tk_yyy_xxyyy[i] * pa_y[i] + 2.0 * ts_yyyy_xxyyy[i] * fz_0;

        tk_yyyy_xxyyz[i] = -6.0 * ts_yy_xxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyyz[i] * fe_0 + 2.0 * tk_yyy_xxyz[i] * fe_0 +
                           tk_yyy_xxyyz[i] * pa_y[i] + 2.0 * ts_yyyy_xxyyz[i] * fz_0;

        tk_yyyy_xxyzz[i] = -6.0 * ts_yy_xxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyzz[i] * fe_0 + tk_yyy_xxzz[i] * fe_0 + tk_yyy_xxyzz[i] * pa_y[i] +
                           2.0 * ts_yyyy_xxyzz[i] * fz_0;

        tk_yyyy_xxzzz[i] =
            -6.0 * ts_yy_xxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxzzz[i] * fe_0 + tk_yyy_xxzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxzzz[i] * fz_0;

        tk_yyyy_xyyyy[i] = -6.0 * ts_yy_xyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyyy[i] * fe_0 + 4.0 * tk_yyy_xyyy[i] * fe_0 +
                           tk_yyy_xyyyy[i] * pa_y[i] + 2.0 * ts_yyyy_xyyyy[i] * fz_0;

        tk_yyyy_xyyyz[i] = -6.0 * ts_yy_xyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyyz[i] * fe_0 + 3.0 * tk_yyy_xyyz[i] * fe_0 +
                           tk_yyy_xyyyz[i] * pa_y[i] + 2.0 * ts_yyyy_xyyyz[i] * fz_0;

        tk_yyyy_xyyzz[i] = -6.0 * ts_yy_xyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyzz[i] * fe_0 + 2.0 * tk_yyy_xyzz[i] * fe_0 +
                           tk_yyy_xyyzz[i] * pa_y[i] + 2.0 * ts_yyyy_xyyzz[i] * fz_0;

        tk_yyyy_xyzzz[i] = -6.0 * ts_yy_xyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyzzz[i] * fe_0 + tk_yyy_xzzz[i] * fe_0 + tk_yyy_xyzzz[i] * pa_y[i] +
                           2.0 * ts_yyyy_xyzzz[i] * fz_0;

        tk_yyyy_xzzzz[i] =
            -6.0 * ts_yy_xzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xzzzz[i] * fe_0 + tk_yyy_xzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xzzzz[i] * fz_0;

        tk_yyyy_yyyyy[i] = -6.0 * ts_yy_yyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyyy[i] * fe_0 + 5.0 * tk_yyy_yyyy[i] * fe_0 +
                           tk_yyy_yyyyy[i] * pa_y[i] + 2.0 * ts_yyyy_yyyyy[i] * fz_0;

        tk_yyyy_yyyyz[i] = -6.0 * ts_yy_yyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyyz[i] * fe_0 + 4.0 * tk_yyy_yyyz[i] * fe_0 +
                           tk_yyy_yyyyz[i] * pa_y[i] + 2.0 * ts_yyyy_yyyyz[i] * fz_0;

        tk_yyyy_yyyzz[i] = -6.0 * ts_yy_yyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyzz[i] * fe_0 + 3.0 * tk_yyy_yyzz[i] * fe_0 +
                           tk_yyy_yyyzz[i] * pa_y[i] + 2.0 * ts_yyyy_yyyzz[i] * fz_0;

        tk_yyyy_yyzzz[i] = -6.0 * ts_yy_yyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyzzz[i] * fe_0 + 2.0 * tk_yyy_yzzz[i] * fe_0 +
                           tk_yyy_yyzzz[i] * pa_y[i] + 2.0 * ts_yyyy_yyzzz[i] * fz_0;

        tk_yyyy_yzzzz[i] = -6.0 * ts_yy_yzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yzzzz[i] * fe_0 + tk_yyy_zzzz[i] * fe_0 + tk_yyy_yzzzz[i] * pa_y[i] +
                           2.0 * ts_yyyy_yzzzz[i] * fz_0;

        tk_yyyy_zzzzz[i] =
            -6.0 * ts_yy_zzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_zzzzz[i] * fe_0 + tk_yyy_zzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_zzzzz[i] * fz_0;
    }

    // Set up 231-252 components of targeted buffer : GH

    auto tk_yyyz_xxxxx = pbuffer.data(idx_kin_gh + 231);

    auto tk_yyyz_xxxxy = pbuffer.data(idx_kin_gh + 232);

    auto tk_yyyz_xxxxz = pbuffer.data(idx_kin_gh + 233);

    auto tk_yyyz_xxxyy = pbuffer.data(idx_kin_gh + 234);

    auto tk_yyyz_xxxyz = pbuffer.data(idx_kin_gh + 235);

    auto tk_yyyz_xxxzz = pbuffer.data(idx_kin_gh + 236);

    auto tk_yyyz_xxyyy = pbuffer.data(idx_kin_gh + 237);

    auto tk_yyyz_xxyyz = pbuffer.data(idx_kin_gh + 238);

    auto tk_yyyz_xxyzz = pbuffer.data(idx_kin_gh + 239);

    auto tk_yyyz_xxzzz = pbuffer.data(idx_kin_gh + 240);

    auto tk_yyyz_xyyyy = pbuffer.data(idx_kin_gh + 241);

    auto tk_yyyz_xyyyz = pbuffer.data(idx_kin_gh + 242);

    auto tk_yyyz_xyyzz = pbuffer.data(idx_kin_gh + 243);

    auto tk_yyyz_xyzzz = pbuffer.data(idx_kin_gh + 244);

    auto tk_yyyz_xzzzz = pbuffer.data(idx_kin_gh + 245);

    auto tk_yyyz_yyyyy = pbuffer.data(idx_kin_gh + 246);

    auto tk_yyyz_yyyyz = pbuffer.data(idx_kin_gh + 247);

    auto tk_yyyz_yyyzz = pbuffer.data(idx_kin_gh + 248);

    auto tk_yyyz_yyzzz = pbuffer.data(idx_kin_gh + 249);

    auto tk_yyyz_yzzzz = pbuffer.data(idx_kin_gh + 250);

    auto tk_yyyz_zzzzz = pbuffer.data(idx_kin_gh + 251);

#pragma omp simd aligned(pa_z,              \
                             tk_yyy_xxxx,   \
                             tk_yyy_xxxxx,  \
                             tk_yyy_xxxxy,  \
                             tk_yyy_xxxxz,  \
                             tk_yyy_xxxy,   \
                             tk_yyy_xxxyy,  \
                             tk_yyy_xxxyz,  \
                             tk_yyy_xxxz,   \
                             tk_yyy_xxxzz,  \
                             tk_yyy_xxyy,   \
                             tk_yyy_xxyyy,  \
                             tk_yyy_xxyyz,  \
                             tk_yyy_xxyz,   \
                             tk_yyy_xxyzz,  \
                             tk_yyy_xxzz,   \
                             tk_yyy_xxzzz,  \
                             tk_yyy_xyyy,   \
                             tk_yyy_xyyyy,  \
                             tk_yyy_xyyyz,  \
                             tk_yyy_xyyz,   \
                             tk_yyy_xyyzz,  \
                             tk_yyy_xyzz,   \
                             tk_yyy_xyzzz,  \
                             tk_yyy_xzzz,   \
                             tk_yyy_xzzzz,  \
                             tk_yyy_yyyy,   \
                             tk_yyy_yyyyy,  \
                             tk_yyy_yyyyz,  \
                             tk_yyy_yyyz,   \
                             tk_yyy_yyyzz,  \
                             tk_yyy_yyzz,   \
                             tk_yyy_yyzzz,  \
                             tk_yyy_yzzz,   \
                             tk_yyy_yzzzz,  \
                             tk_yyy_zzzz,   \
                             tk_yyy_zzzzz,  \
                             tk_yyyz_xxxxx, \
                             tk_yyyz_xxxxy, \
                             tk_yyyz_xxxxz, \
                             tk_yyyz_xxxyy, \
                             tk_yyyz_xxxyz, \
                             tk_yyyz_xxxzz, \
                             tk_yyyz_xxyyy, \
                             tk_yyyz_xxyyz, \
                             tk_yyyz_xxyzz, \
                             tk_yyyz_xxzzz, \
                             tk_yyyz_xyyyy, \
                             tk_yyyz_xyyyz, \
                             tk_yyyz_xyyzz, \
                             tk_yyyz_xyzzz, \
                             tk_yyyz_xzzzz, \
                             tk_yyyz_yyyyy, \
                             tk_yyyz_yyyyz, \
                             tk_yyyz_yyyzz, \
                             tk_yyyz_yyzzz, \
                             tk_yyyz_yzzzz, \
                             tk_yyyz_zzzzz, \
                             ts_yyyz_xxxxx, \
                             ts_yyyz_xxxxy, \
                             ts_yyyz_xxxxz, \
                             ts_yyyz_xxxyy, \
                             ts_yyyz_xxxyz, \
                             ts_yyyz_xxxzz, \
                             ts_yyyz_xxyyy, \
                             ts_yyyz_xxyyz, \
                             ts_yyyz_xxyzz, \
                             ts_yyyz_xxzzz, \
                             ts_yyyz_xyyyy, \
                             ts_yyyz_xyyyz, \
                             ts_yyyz_xyyzz, \
                             ts_yyyz_xyzzz, \
                             ts_yyyz_xzzzz, \
                             ts_yyyz_yyyyy, \
                             ts_yyyz_yyyyz, \
                             ts_yyyz_yyyzz, \
                             ts_yyyz_yyzzz, \
                             ts_yyyz_yzzzz, \
                             ts_yyyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyz_xxxxx[i] = tk_yyy_xxxxx[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxx[i] * fz_0;

        tk_yyyz_xxxxy[i] = tk_yyy_xxxxy[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxy[i] * fz_0;

        tk_yyyz_xxxxz[i] = tk_yyy_xxxx[i] * fe_0 + tk_yyy_xxxxz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxz[i] * fz_0;

        tk_yyyz_xxxyy[i] = tk_yyy_xxxyy[i] * pa_z[i] + 2.0 * ts_yyyz_xxxyy[i] * fz_0;

        tk_yyyz_xxxyz[i] = tk_yyy_xxxy[i] * fe_0 + tk_yyy_xxxyz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxyz[i] * fz_0;

        tk_yyyz_xxxzz[i] = 2.0 * tk_yyy_xxxz[i] * fe_0 + tk_yyy_xxxzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxzz[i] * fz_0;

        tk_yyyz_xxyyy[i] = tk_yyy_xxyyy[i] * pa_z[i] + 2.0 * ts_yyyz_xxyyy[i] * fz_0;

        tk_yyyz_xxyyz[i] = tk_yyy_xxyy[i] * fe_0 + tk_yyy_xxyyz[i] * pa_z[i] + 2.0 * ts_yyyz_xxyyz[i] * fz_0;

        tk_yyyz_xxyzz[i] = 2.0 * tk_yyy_xxyz[i] * fe_0 + tk_yyy_xxyzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxyzz[i] * fz_0;

        tk_yyyz_xxzzz[i] = 3.0 * tk_yyy_xxzz[i] * fe_0 + tk_yyy_xxzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxzzz[i] * fz_0;

        tk_yyyz_xyyyy[i] = tk_yyy_xyyyy[i] * pa_z[i] + 2.0 * ts_yyyz_xyyyy[i] * fz_0;

        tk_yyyz_xyyyz[i] = tk_yyy_xyyy[i] * fe_0 + tk_yyy_xyyyz[i] * pa_z[i] + 2.0 * ts_yyyz_xyyyz[i] * fz_0;

        tk_yyyz_xyyzz[i] = 2.0 * tk_yyy_xyyz[i] * fe_0 + tk_yyy_xyyzz[i] * pa_z[i] + 2.0 * ts_yyyz_xyyzz[i] * fz_0;

        tk_yyyz_xyzzz[i] = 3.0 * tk_yyy_xyzz[i] * fe_0 + tk_yyy_xyzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xyzzz[i] * fz_0;

        tk_yyyz_xzzzz[i] = 4.0 * tk_yyy_xzzz[i] * fe_0 + tk_yyy_xzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xzzzz[i] * fz_0;

        tk_yyyz_yyyyy[i] = tk_yyy_yyyyy[i] * pa_z[i] + 2.0 * ts_yyyz_yyyyy[i] * fz_0;

        tk_yyyz_yyyyz[i] = tk_yyy_yyyy[i] * fe_0 + tk_yyy_yyyyz[i] * pa_z[i] + 2.0 * ts_yyyz_yyyyz[i] * fz_0;

        tk_yyyz_yyyzz[i] = 2.0 * tk_yyy_yyyz[i] * fe_0 + tk_yyy_yyyzz[i] * pa_z[i] + 2.0 * ts_yyyz_yyyzz[i] * fz_0;

        tk_yyyz_yyzzz[i] = 3.0 * tk_yyy_yyzz[i] * fe_0 + tk_yyy_yyzzz[i] * pa_z[i] + 2.0 * ts_yyyz_yyzzz[i] * fz_0;

        tk_yyyz_yzzzz[i] = 4.0 * tk_yyy_yzzz[i] * fe_0 + tk_yyy_yzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_yzzzz[i] * fz_0;

        tk_yyyz_zzzzz[i] = 5.0 * tk_yyy_zzzz[i] * fe_0 + tk_yyy_zzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_zzzzz[i] * fz_0;
    }

    // Set up 252-273 components of targeted buffer : GH

    auto tk_yyzz_xxxxx = pbuffer.data(idx_kin_gh + 252);

    auto tk_yyzz_xxxxy = pbuffer.data(idx_kin_gh + 253);

    auto tk_yyzz_xxxxz = pbuffer.data(idx_kin_gh + 254);

    auto tk_yyzz_xxxyy = pbuffer.data(idx_kin_gh + 255);

    auto tk_yyzz_xxxyz = pbuffer.data(idx_kin_gh + 256);

    auto tk_yyzz_xxxzz = pbuffer.data(idx_kin_gh + 257);

    auto tk_yyzz_xxyyy = pbuffer.data(idx_kin_gh + 258);

    auto tk_yyzz_xxyyz = pbuffer.data(idx_kin_gh + 259);

    auto tk_yyzz_xxyzz = pbuffer.data(idx_kin_gh + 260);

    auto tk_yyzz_xxzzz = pbuffer.data(idx_kin_gh + 261);

    auto tk_yyzz_xyyyy = pbuffer.data(idx_kin_gh + 262);

    auto tk_yyzz_xyyyz = pbuffer.data(idx_kin_gh + 263);

    auto tk_yyzz_xyyzz = pbuffer.data(idx_kin_gh + 264);

    auto tk_yyzz_xyzzz = pbuffer.data(idx_kin_gh + 265);

    auto tk_yyzz_xzzzz = pbuffer.data(idx_kin_gh + 266);

    auto tk_yyzz_yyyyy = pbuffer.data(idx_kin_gh + 267);

    auto tk_yyzz_yyyyz = pbuffer.data(idx_kin_gh + 268);

    auto tk_yyzz_yyyzz = pbuffer.data(idx_kin_gh + 269);

    auto tk_yyzz_yyzzz = pbuffer.data(idx_kin_gh + 270);

    auto tk_yyzz_yzzzz = pbuffer.data(idx_kin_gh + 271);

    auto tk_yyzz_zzzzz = pbuffer.data(idx_kin_gh + 272);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tk_yy_xxxxy,   \
                             tk_yy_xxxyy,   \
                             tk_yy_xxyyy,   \
                             tk_yy_xyyyy,   \
                             tk_yy_yyyyy,   \
                             tk_yyz_xxxxy,  \
                             tk_yyz_xxxyy,  \
                             tk_yyz_xxyyy,  \
                             tk_yyz_xyyyy,  \
                             tk_yyz_yyyyy,  \
                             tk_yyzz_xxxxx, \
                             tk_yyzz_xxxxy, \
                             tk_yyzz_xxxxz, \
                             tk_yyzz_xxxyy, \
                             tk_yyzz_xxxyz, \
                             tk_yyzz_xxxzz, \
                             tk_yyzz_xxyyy, \
                             tk_yyzz_xxyyz, \
                             tk_yyzz_xxyzz, \
                             tk_yyzz_xxzzz, \
                             tk_yyzz_xyyyy, \
                             tk_yyzz_xyyyz, \
                             tk_yyzz_xyyzz, \
                             tk_yyzz_xyzzz, \
                             tk_yyzz_xzzzz, \
                             tk_yyzz_yyyyy, \
                             tk_yyzz_yyyyz, \
                             tk_yyzz_yyyzz, \
                             tk_yyzz_yyzzz, \
                             tk_yyzz_yzzzz, \
                             tk_yyzz_zzzzz, \
                             tk_yzz_xxxxx,  \
                             tk_yzz_xxxxz,  \
                             tk_yzz_xxxyz,  \
                             tk_yzz_xxxz,   \
                             tk_yzz_xxxzz,  \
                             tk_yzz_xxyyz,  \
                             tk_yzz_xxyz,   \
                             tk_yzz_xxyzz,  \
                             tk_yzz_xxzz,   \
                             tk_yzz_xxzzz,  \
                             tk_yzz_xyyyz,  \
                             tk_yzz_xyyz,   \
                             tk_yzz_xyyzz,  \
                             tk_yzz_xyzz,   \
                             tk_yzz_xyzzz,  \
                             tk_yzz_xzzz,   \
                             tk_yzz_xzzzz,  \
                             tk_yzz_yyyyz,  \
                             tk_yzz_yyyz,   \
                             tk_yzz_yyyzz,  \
                             tk_yzz_yyzz,   \
                             tk_yzz_yyzzz,  \
                             tk_yzz_yzzz,   \
                             tk_yzz_yzzzz,  \
                             tk_yzz_zzzz,   \
                             tk_yzz_zzzzz,  \
                             tk_zz_xxxxx,   \
                             tk_zz_xxxxz,   \
                             tk_zz_xxxyz,   \
                             tk_zz_xxxzz,   \
                             tk_zz_xxyyz,   \
                             tk_zz_xxyzz,   \
                             tk_zz_xxzzz,   \
                             tk_zz_xyyyz,   \
                             tk_zz_xyyzz,   \
                             tk_zz_xyzzz,   \
                             tk_zz_xzzzz,   \
                             tk_zz_yyyyz,   \
                             tk_zz_yyyzz,   \
                             tk_zz_yyzzz,   \
                             tk_zz_yzzzz,   \
                             tk_zz_zzzzz,   \
                             ts_yy_xxxxy,   \
                             ts_yy_xxxyy,   \
                             ts_yy_xxyyy,   \
                             ts_yy_xyyyy,   \
                             ts_yy_yyyyy,   \
                             ts_yyzz_xxxxx, \
                             ts_yyzz_xxxxy, \
                             ts_yyzz_xxxxz, \
                             ts_yyzz_xxxyy, \
                             ts_yyzz_xxxyz, \
                             ts_yyzz_xxxzz, \
                             ts_yyzz_xxyyy, \
                             ts_yyzz_xxyyz, \
                             ts_yyzz_xxyzz, \
                             ts_yyzz_xxzzz, \
                             ts_yyzz_xyyyy, \
                             ts_yyzz_xyyyz, \
                             ts_yyzz_xyyzz, \
                             ts_yyzz_xyzzz, \
                             ts_yyzz_xzzzz, \
                             ts_yyzz_yyyyy, \
                             ts_yyzz_yyyyz, \
                             ts_yyzz_yyyzz, \
                             ts_yyzz_yyzzz, \
                             ts_yyzz_yzzzz, \
                             ts_yyzz_zzzzz, \
                             ts_zz_xxxxx,   \
                             ts_zz_xxxxz,   \
                             ts_zz_xxxyz,   \
                             ts_zz_xxxzz,   \
                             ts_zz_xxyyz,   \
                             ts_zz_xxyzz,   \
                             ts_zz_xxzzz,   \
                             ts_zz_xyyyz,   \
                             ts_zz_xyyzz,   \
                             ts_zz_xyzzz,   \
                             ts_zz_xzzzz,   \
                             ts_zz_yyyyz,   \
                             ts_zz_yyyzz,   \
                             ts_zz_yyzzz,   \
                             ts_zz_yzzzz,   \
                             ts_zz_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzz_xxxxx[i] = -2.0 * ts_zz_xxxxx[i] * fbe_0 * fz_0 + tk_zz_xxxxx[i] * fe_0 + tk_yzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yyzz_xxxxx[i] * fz_0;

        tk_yyzz_xxxxy[i] = -2.0 * ts_yy_xxxxy[i] * fbe_0 * fz_0 + tk_yy_xxxxy[i] * fe_0 + tk_yyz_xxxxy[i] * pa_z[i] + 2.0 * ts_yyzz_xxxxy[i] * fz_0;

        tk_yyzz_xxxxz[i] = -2.0 * ts_zz_xxxxz[i] * fbe_0 * fz_0 + tk_zz_xxxxz[i] * fe_0 + tk_yzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yyzz_xxxxz[i] * fz_0;

        tk_yyzz_xxxyy[i] = -2.0 * ts_yy_xxxyy[i] * fbe_0 * fz_0 + tk_yy_xxxyy[i] * fe_0 + tk_yyz_xxxyy[i] * pa_z[i] + 2.0 * ts_yyzz_xxxyy[i] * fz_0;

        tk_yyzz_xxxyz[i] = -2.0 * ts_zz_xxxyz[i] * fbe_0 * fz_0 + tk_zz_xxxyz[i] * fe_0 + tk_yzz_xxxz[i] * fe_0 + tk_yzz_xxxyz[i] * pa_y[i] +
                           2.0 * ts_yyzz_xxxyz[i] * fz_0;

        tk_yyzz_xxxzz[i] = -2.0 * ts_zz_xxxzz[i] * fbe_0 * fz_0 + tk_zz_xxxzz[i] * fe_0 + tk_yzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yyzz_xxxzz[i] * fz_0;

        tk_yyzz_xxyyy[i] = -2.0 * ts_yy_xxyyy[i] * fbe_0 * fz_0 + tk_yy_xxyyy[i] * fe_0 + tk_yyz_xxyyy[i] * pa_z[i] + 2.0 * ts_yyzz_xxyyy[i] * fz_0;

        tk_yyzz_xxyyz[i] = -2.0 * ts_zz_xxyyz[i] * fbe_0 * fz_0 + tk_zz_xxyyz[i] * fe_0 + 2.0 * tk_yzz_xxyz[i] * fe_0 + tk_yzz_xxyyz[i] * pa_y[i] +
                           2.0 * ts_yyzz_xxyyz[i] * fz_0;

        tk_yyzz_xxyzz[i] = -2.0 * ts_zz_xxyzz[i] * fbe_0 * fz_0 + tk_zz_xxyzz[i] * fe_0 + tk_yzz_xxzz[i] * fe_0 + tk_yzz_xxyzz[i] * pa_y[i] +
                           2.0 * ts_yyzz_xxyzz[i] * fz_0;

        tk_yyzz_xxzzz[i] = -2.0 * ts_zz_xxzzz[i] * fbe_0 * fz_0 + tk_zz_xxzzz[i] * fe_0 + tk_yzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yyzz_xxzzz[i] * fz_0;

        tk_yyzz_xyyyy[i] = -2.0 * ts_yy_xyyyy[i] * fbe_0 * fz_0 + tk_yy_xyyyy[i] * fe_0 + tk_yyz_xyyyy[i] * pa_z[i] + 2.0 * ts_yyzz_xyyyy[i] * fz_0;

        tk_yyzz_xyyyz[i] = -2.0 * ts_zz_xyyyz[i] * fbe_0 * fz_0 + tk_zz_xyyyz[i] * fe_0 + 3.0 * tk_yzz_xyyz[i] * fe_0 + tk_yzz_xyyyz[i] * pa_y[i] +
                           2.0 * ts_yyzz_xyyyz[i] * fz_0;

        tk_yyzz_xyyzz[i] = -2.0 * ts_zz_xyyzz[i] * fbe_0 * fz_0 + tk_zz_xyyzz[i] * fe_0 + 2.0 * tk_yzz_xyzz[i] * fe_0 + tk_yzz_xyyzz[i] * pa_y[i] +
                           2.0 * ts_yyzz_xyyzz[i] * fz_0;

        tk_yyzz_xyzzz[i] = -2.0 * ts_zz_xyzzz[i] * fbe_0 * fz_0 + tk_zz_xyzzz[i] * fe_0 + tk_yzz_xzzz[i] * fe_0 + tk_yzz_xyzzz[i] * pa_y[i] +
                           2.0 * ts_yyzz_xyzzz[i] * fz_0;

        tk_yyzz_xzzzz[i] = -2.0 * ts_zz_xzzzz[i] * fbe_0 * fz_0 + tk_zz_xzzzz[i] * fe_0 + tk_yzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yyzz_xzzzz[i] * fz_0;

        tk_yyzz_yyyyy[i] = -2.0 * ts_yy_yyyyy[i] * fbe_0 * fz_0 + tk_yy_yyyyy[i] * fe_0 + tk_yyz_yyyyy[i] * pa_z[i] + 2.0 * ts_yyzz_yyyyy[i] * fz_0;

        tk_yyzz_yyyyz[i] = -2.0 * ts_zz_yyyyz[i] * fbe_0 * fz_0 + tk_zz_yyyyz[i] * fe_0 + 4.0 * tk_yzz_yyyz[i] * fe_0 + tk_yzz_yyyyz[i] * pa_y[i] +
                           2.0 * ts_yyzz_yyyyz[i] * fz_0;

        tk_yyzz_yyyzz[i] = -2.0 * ts_zz_yyyzz[i] * fbe_0 * fz_0 + tk_zz_yyyzz[i] * fe_0 + 3.0 * tk_yzz_yyzz[i] * fe_0 + tk_yzz_yyyzz[i] * pa_y[i] +
                           2.0 * ts_yyzz_yyyzz[i] * fz_0;

        tk_yyzz_yyzzz[i] = -2.0 * ts_zz_yyzzz[i] * fbe_0 * fz_0 + tk_zz_yyzzz[i] * fe_0 + 2.0 * tk_yzz_yzzz[i] * fe_0 + tk_yzz_yyzzz[i] * pa_y[i] +
                           2.0 * ts_yyzz_yyzzz[i] * fz_0;

        tk_yyzz_yzzzz[i] = -2.0 * ts_zz_yzzzz[i] * fbe_0 * fz_0 + tk_zz_yzzzz[i] * fe_0 + tk_yzz_zzzz[i] * fe_0 + tk_yzz_yzzzz[i] * pa_y[i] +
                           2.0 * ts_yyzz_yzzzz[i] * fz_0;

        tk_yyzz_zzzzz[i] = -2.0 * ts_zz_zzzzz[i] * fbe_0 * fz_0 + tk_zz_zzzzz[i] * fe_0 + tk_yzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yyzz_zzzzz[i] * fz_0;
    }

    // Set up 273-294 components of targeted buffer : GH

    auto tk_yzzz_xxxxx = pbuffer.data(idx_kin_gh + 273);

    auto tk_yzzz_xxxxy = pbuffer.data(idx_kin_gh + 274);

    auto tk_yzzz_xxxxz = pbuffer.data(idx_kin_gh + 275);

    auto tk_yzzz_xxxyy = pbuffer.data(idx_kin_gh + 276);

    auto tk_yzzz_xxxyz = pbuffer.data(idx_kin_gh + 277);

    auto tk_yzzz_xxxzz = pbuffer.data(idx_kin_gh + 278);

    auto tk_yzzz_xxyyy = pbuffer.data(idx_kin_gh + 279);

    auto tk_yzzz_xxyyz = pbuffer.data(idx_kin_gh + 280);

    auto tk_yzzz_xxyzz = pbuffer.data(idx_kin_gh + 281);

    auto tk_yzzz_xxzzz = pbuffer.data(idx_kin_gh + 282);

    auto tk_yzzz_xyyyy = pbuffer.data(idx_kin_gh + 283);

    auto tk_yzzz_xyyyz = pbuffer.data(idx_kin_gh + 284);

    auto tk_yzzz_xyyzz = pbuffer.data(idx_kin_gh + 285);

    auto tk_yzzz_xyzzz = pbuffer.data(idx_kin_gh + 286);

    auto tk_yzzz_xzzzz = pbuffer.data(idx_kin_gh + 287);

    auto tk_yzzz_yyyyy = pbuffer.data(idx_kin_gh + 288);

    auto tk_yzzz_yyyyz = pbuffer.data(idx_kin_gh + 289);

    auto tk_yzzz_yyyzz = pbuffer.data(idx_kin_gh + 290);

    auto tk_yzzz_yyzzz = pbuffer.data(idx_kin_gh + 291);

    auto tk_yzzz_yzzzz = pbuffer.data(idx_kin_gh + 292);

    auto tk_yzzz_zzzzz = pbuffer.data(idx_kin_gh + 293);

#pragma omp simd aligned(pa_y,              \
                             tk_yzzz_xxxxx, \
                             tk_yzzz_xxxxy, \
                             tk_yzzz_xxxxz, \
                             tk_yzzz_xxxyy, \
                             tk_yzzz_xxxyz, \
                             tk_yzzz_xxxzz, \
                             tk_yzzz_xxyyy, \
                             tk_yzzz_xxyyz, \
                             tk_yzzz_xxyzz, \
                             tk_yzzz_xxzzz, \
                             tk_yzzz_xyyyy, \
                             tk_yzzz_xyyyz, \
                             tk_yzzz_xyyzz, \
                             tk_yzzz_xyzzz, \
                             tk_yzzz_xzzzz, \
                             tk_yzzz_yyyyy, \
                             tk_yzzz_yyyyz, \
                             tk_yzzz_yyyzz, \
                             tk_yzzz_yyzzz, \
                             tk_yzzz_yzzzz, \
                             tk_yzzz_zzzzz, \
                             tk_zzz_xxxx,   \
                             tk_zzz_xxxxx,  \
                             tk_zzz_xxxxy,  \
                             tk_zzz_xxxxz,  \
                             tk_zzz_xxxy,   \
                             tk_zzz_xxxyy,  \
                             tk_zzz_xxxyz,  \
                             tk_zzz_xxxz,   \
                             tk_zzz_xxxzz,  \
                             tk_zzz_xxyy,   \
                             tk_zzz_xxyyy,  \
                             tk_zzz_xxyyz,  \
                             tk_zzz_xxyz,   \
                             tk_zzz_xxyzz,  \
                             tk_zzz_xxzz,   \
                             tk_zzz_xxzzz,  \
                             tk_zzz_xyyy,   \
                             tk_zzz_xyyyy,  \
                             tk_zzz_xyyyz,  \
                             tk_zzz_xyyz,   \
                             tk_zzz_xyyzz,  \
                             tk_zzz_xyzz,   \
                             tk_zzz_xyzzz,  \
                             tk_zzz_xzzz,   \
                             tk_zzz_xzzzz,  \
                             tk_zzz_yyyy,   \
                             tk_zzz_yyyyy,  \
                             tk_zzz_yyyyz,  \
                             tk_zzz_yyyz,   \
                             tk_zzz_yyyzz,  \
                             tk_zzz_yyzz,   \
                             tk_zzz_yyzzz,  \
                             tk_zzz_yzzz,   \
                             tk_zzz_yzzzz,  \
                             tk_zzz_zzzz,   \
                             tk_zzz_zzzzz,  \
                             ts_yzzz_xxxxx, \
                             ts_yzzz_xxxxy, \
                             ts_yzzz_xxxxz, \
                             ts_yzzz_xxxyy, \
                             ts_yzzz_xxxyz, \
                             ts_yzzz_xxxzz, \
                             ts_yzzz_xxyyy, \
                             ts_yzzz_xxyyz, \
                             ts_yzzz_xxyzz, \
                             ts_yzzz_xxzzz, \
                             ts_yzzz_xyyyy, \
                             ts_yzzz_xyyyz, \
                             ts_yzzz_xyyzz, \
                             ts_yzzz_xyzzz, \
                             ts_yzzz_xzzzz, \
                             ts_yzzz_yyyyy, \
                             ts_yzzz_yyyyz, \
                             ts_yzzz_yyyzz, \
                             ts_yzzz_yyzzz, \
                             ts_yzzz_yzzzz, \
                             ts_yzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzz_xxxxx[i] = tk_zzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxx[i] * fz_0;

        tk_yzzz_xxxxy[i] = tk_zzz_xxxx[i] * fe_0 + tk_zzz_xxxxy[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxy[i] * fz_0;

        tk_yzzz_xxxxz[i] = tk_zzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxz[i] * fz_0;

        tk_yzzz_xxxyy[i] = 2.0 * tk_zzz_xxxy[i] * fe_0 + tk_zzz_xxxyy[i] * pa_y[i] + 2.0 * ts_yzzz_xxxyy[i] * fz_0;

        tk_yzzz_xxxyz[i] = tk_zzz_xxxz[i] * fe_0 + tk_zzz_xxxyz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxyz[i] * fz_0;

        tk_yzzz_xxxzz[i] = tk_zzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxzz[i] * fz_0;

        tk_yzzz_xxyyy[i] = 3.0 * tk_zzz_xxyy[i] * fe_0 + tk_zzz_xxyyy[i] * pa_y[i] + 2.0 * ts_yzzz_xxyyy[i] * fz_0;

        tk_yzzz_xxyyz[i] = 2.0 * tk_zzz_xxyz[i] * fe_0 + tk_zzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yzzz_xxyyz[i] * fz_0;

        tk_yzzz_xxyzz[i] = tk_zzz_xxzz[i] * fe_0 + tk_zzz_xxyzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxyzz[i] * fz_0;

        tk_yzzz_xxzzz[i] = tk_zzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxzzz[i] * fz_0;

        tk_yzzz_xyyyy[i] = 4.0 * tk_zzz_xyyy[i] * fe_0 + tk_zzz_xyyyy[i] * pa_y[i] + 2.0 * ts_yzzz_xyyyy[i] * fz_0;

        tk_yzzz_xyyyz[i] = 3.0 * tk_zzz_xyyz[i] * fe_0 + tk_zzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yzzz_xyyyz[i] * fz_0;

        tk_yzzz_xyyzz[i] = 2.0 * tk_zzz_xyzz[i] * fe_0 + tk_zzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yzzz_xyyzz[i] * fz_0;

        tk_yzzz_xyzzz[i] = tk_zzz_xzzz[i] * fe_0 + tk_zzz_xyzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xyzzz[i] * fz_0;

        tk_yzzz_xzzzz[i] = tk_zzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xzzzz[i] * fz_0;

        tk_yzzz_yyyyy[i] = 5.0 * tk_zzz_yyyy[i] * fe_0 + tk_zzz_yyyyy[i] * pa_y[i] + 2.0 * ts_yzzz_yyyyy[i] * fz_0;

        tk_yzzz_yyyyz[i] = 4.0 * tk_zzz_yyyz[i] * fe_0 + tk_zzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yzzz_yyyyz[i] * fz_0;

        tk_yzzz_yyyzz[i] = 3.0 * tk_zzz_yyzz[i] * fe_0 + tk_zzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yzzz_yyyzz[i] * fz_0;

        tk_yzzz_yyzzz[i] = 2.0 * tk_zzz_yzzz[i] * fe_0 + tk_zzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yzzz_yyzzz[i] * fz_0;

        tk_yzzz_yzzzz[i] = tk_zzz_zzzz[i] * fe_0 + tk_zzz_yzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_yzzzz[i] * fz_0;

        tk_yzzz_zzzzz[i] = tk_zzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_zzzzz[i] * fz_0;
    }

    // Set up 294-315 components of targeted buffer : GH

    auto tk_zzzz_xxxxx = pbuffer.data(idx_kin_gh + 294);

    auto tk_zzzz_xxxxy = pbuffer.data(idx_kin_gh + 295);

    auto tk_zzzz_xxxxz = pbuffer.data(idx_kin_gh + 296);

    auto tk_zzzz_xxxyy = pbuffer.data(idx_kin_gh + 297);

    auto tk_zzzz_xxxyz = pbuffer.data(idx_kin_gh + 298);

    auto tk_zzzz_xxxzz = pbuffer.data(idx_kin_gh + 299);

    auto tk_zzzz_xxyyy = pbuffer.data(idx_kin_gh + 300);

    auto tk_zzzz_xxyyz = pbuffer.data(idx_kin_gh + 301);

    auto tk_zzzz_xxyzz = pbuffer.data(idx_kin_gh + 302);

    auto tk_zzzz_xxzzz = pbuffer.data(idx_kin_gh + 303);

    auto tk_zzzz_xyyyy = pbuffer.data(idx_kin_gh + 304);

    auto tk_zzzz_xyyyz = pbuffer.data(idx_kin_gh + 305);

    auto tk_zzzz_xyyzz = pbuffer.data(idx_kin_gh + 306);

    auto tk_zzzz_xyzzz = pbuffer.data(idx_kin_gh + 307);

    auto tk_zzzz_xzzzz = pbuffer.data(idx_kin_gh + 308);

    auto tk_zzzz_yyyyy = pbuffer.data(idx_kin_gh + 309);

    auto tk_zzzz_yyyyz = pbuffer.data(idx_kin_gh + 310);

    auto tk_zzzz_yyyzz = pbuffer.data(idx_kin_gh + 311);

    auto tk_zzzz_yyzzz = pbuffer.data(idx_kin_gh + 312);

    auto tk_zzzz_yzzzz = pbuffer.data(idx_kin_gh + 313);

    auto tk_zzzz_zzzzz = pbuffer.data(idx_kin_gh + 314);

#pragma omp simd aligned(pa_z,              \
                             tk_zz_xxxxx,   \
                             tk_zz_xxxxy,   \
                             tk_zz_xxxxz,   \
                             tk_zz_xxxyy,   \
                             tk_zz_xxxyz,   \
                             tk_zz_xxxzz,   \
                             tk_zz_xxyyy,   \
                             tk_zz_xxyyz,   \
                             tk_zz_xxyzz,   \
                             tk_zz_xxzzz,   \
                             tk_zz_xyyyy,   \
                             tk_zz_xyyyz,   \
                             tk_zz_xyyzz,   \
                             tk_zz_xyzzz,   \
                             tk_zz_xzzzz,   \
                             tk_zz_yyyyy,   \
                             tk_zz_yyyyz,   \
                             tk_zz_yyyzz,   \
                             tk_zz_yyzzz,   \
                             tk_zz_yzzzz,   \
                             tk_zz_zzzzz,   \
                             tk_zzz_xxxx,   \
                             tk_zzz_xxxxx,  \
                             tk_zzz_xxxxy,  \
                             tk_zzz_xxxxz,  \
                             tk_zzz_xxxy,   \
                             tk_zzz_xxxyy,  \
                             tk_zzz_xxxyz,  \
                             tk_zzz_xxxz,   \
                             tk_zzz_xxxzz,  \
                             tk_zzz_xxyy,   \
                             tk_zzz_xxyyy,  \
                             tk_zzz_xxyyz,  \
                             tk_zzz_xxyz,   \
                             tk_zzz_xxyzz,  \
                             tk_zzz_xxzz,   \
                             tk_zzz_xxzzz,  \
                             tk_zzz_xyyy,   \
                             tk_zzz_xyyyy,  \
                             tk_zzz_xyyyz,  \
                             tk_zzz_xyyz,   \
                             tk_zzz_xyyzz,  \
                             tk_zzz_xyzz,   \
                             tk_zzz_xyzzz,  \
                             tk_zzz_xzzz,   \
                             tk_zzz_xzzzz,  \
                             tk_zzz_yyyy,   \
                             tk_zzz_yyyyy,  \
                             tk_zzz_yyyyz,  \
                             tk_zzz_yyyz,   \
                             tk_zzz_yyyzz,  \
                             tk_zzz_yyzz,   \
                             tk_zzz_yyzzz,  \
                             tk_zzz_yzzz,   \
                             tk_zzz_yzzzz,  \
                             tk_zzz_zzzz,   \
                             tk_zzz_zzzzz,  \
                             tk_zzzz_xxxxx, \
                             tk_zzzz_xxxxy, \
                             tk_zzzz_xxxxz, \
                             tk_zzzz_xxxyy, \
                             tk_zzzz_xxxyz, \
                             tk_zzzz_xxxzz, \
                             tk_zzzz_xxyyy, \
                             tk_zzzz_xxyyz, \
                             tk_zzzz_xxyzz, \
                             tk_zzzz_xxzzz, \
                             tk_zzzz_xyyyy, \
                             tk_zzzz_xyyyz, \
                             tk_zzzz_xyyzz, \
                             tk_zzzz_xyzzz, \
                             tk_zzzz_xzzzz, \
                             tk_zzzz_yyyyy, \
                             tk_zzzz_yyyyz, \
                             tk_zzzz_yyyzz, \
                             tk_zzzz_yyzzz, \
                             tk_zzzz_yzzzz, \
                             tk_zzzz_zzzzz, \
                             ts_zz_xxxxx,   \
                             ts_zz_xxxxy,   \
                             ts_zz_xxxxz,   \
                             ts_zz_xxxyy,   \
                             ts_zz_xxxyz,   \
                             ts_zz_xxxzz,   \
                             ts_zz_xxyyy,   \
                             ts_zz_xxyyz,   \
                             ts_zz_xxyzz,   \
                             ts_zz_xxzzz,   \
                             ts_zz_xyyyy,   \
                             ts_zz_xyyyz,   \
                             ts_zz_xyyzz,   \
                             ts_zz_xyzzz,   \
                             ts_zz_xzzzz,   \
                             ts_zz_yyyyy,   \
                             ts_zz_yyyyz,   \
                             ts_zz_yyyzz,   \
                             ts_zz_yyzzz,   \
                             ts_zz_yzzzz,   \
                             ts_zz_zzzzz,   \
                             ts_zzzz_xxxxx, \
                             ts_zzzz_xxxxy, \
                             ts_zzzz_xxxxz, \
                             ts_zzzz_xxxyy, \
                             ts_zzzz_xxxyz, \
                             ts_zzzz_xxxzz, \
                             ts_zzzz_xxyyy, \
                             ts_zzzz_xxyyz, \
                             ts_zzzz_xxyzz, \
                             ts_zzzz_xxzzz, \
                             ts_zzzz_xyyyy, \
                             ts_zzzz_xyyyz, \
                             ts_zzzz_xyyzz, \
                             ts_zzzz_xyzzz, \
                             ts_zzzz_xzzzz, \
                             ts_zzzz_yyyyy, \
                             ts_zzzz_yyyyz, \
                             ts_zzzz_yyyzz, \
                             ts_zzzz_yyzzz, \
                             ts_zzzz_yzzzz, \
                             ts_zzzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzz_xxxxx[i] =
            -6.0 * ts_zz_xxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxx[i] * fe_0 + tk_zzz_xxxxx[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxx[i] * fz_0;

        tk_zzzz_xxxxy[i] =
            -6.0 * ts_zz_xxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxy[i] * fe_0 + tk_zzz_xxxxy[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxy[i] * fz_0;

        tk_zzzz_xxxxz[i] = -6.0 * ts_zz_xxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxz[i] * fe_0 + tk_zzz_xxxx[i] * fe_0 + tk_zzz_xxxxz[i] * pa_z[i] +
                           2.0 * ts_zzzz_xxxxz[i] * fz_0;

        tk_zzzz_xxxyy[i] =
            -6.0 * ts_zz_xxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxyy[i] * fe_0 + tk_zzz_xxxyy[i] * pa_z[i] + 2.0 * ts_zzzz_xxxyy[i] * fz_0;

        tk_zzzz_xxxyz[i] = -6.0 * ts_zz_xxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxyz[i] * fe_0 + tk_zzz_xxxy[i] * fe_0 + tk_zzz_xxxyz[i] * pa_z[i] +
                           2.0 * ts_zzzz_xxxyz[i] * fz_0;

        tk_zzzz_xxxzz[i] = -6.0 * ts_zz_xxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxzz[i] * fe_0 + 2.0 * tk_zzz_xxxz[i] * fe_0 +
                           tk_zzz_xxxzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxxzz[i] * fz_0;

        tk_zzzz_xxyyy[i] =
            -6.0 * ts_zz_xxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyyy[i] * fe_0 + tk_zzz_xxyyy[i] * pa_z[i] + 2.0 * ts_zzzz_xxyyy[i] * fz_0;

        tk_zzzz_xxyyz[i] = -6.0 * ts_zz_xxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyyz[i] * fe_0 + tk_zzz_xxyy[i] * fe_0 + tk_zzz_xxyyz[i] * pa_z[i] +
                           2.0 * ts_zzzz_xxyyz[i] * fz_0;

        tk_zzzz_xxyzz[i] = -6.0 * ts_zz_xxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyzz[i] * fe_0 + 2.0 * tk_zzz_xxyz[i] * fe_0 +
                           tk_zzz_xxyzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxyzz[i] * fz_0;

        tk_zzzz_xxzzz[i] = -6.0 * ts_zz_xxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxzzz[i] * fe_0 + 3.0 * tk_zzz_xxzz[i] * fe_0 +
                           tk_zzz_xxzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxzzz[i] * fz_0;

        tk_zzzz_xyyyy[i] =
            -6.0 * ts_zz_xyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyyy[i] * fe_0 + tk_zzz_xyyyy[i] * pa_z[i] + 2.0 * ts_zzzz_xyyyy[i] * fz_0;

        tk_zzzz_xyyyz[i] = -6.0 * ts_zz_xyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyyz[i] * fe_0 + tk_zzz_xyyy[i] * fe_0 + tk_zzz_xyyyz[i] * pa_z[i] +
                           2.0 * ts_zzzz_xyyyz[i] * fz_0;

        tk_zzzz_xyyzz[i] = -6.0 * ts_zz_xyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyzz[i] * fe_0 + 2.0 * tk_zzz_xyyz[i] * fe_0 +
                           tk_zzz_xyyzz[i] * pa_z[i] + 2.0 * ts_zzzz_xyyzz[i] * fz_0;

        tk_zzzz_xyzzz[i] = -6.0 * ts_zz_xyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyzzz[i] * fe_0 + 3.0 * tk_zzz_xyzz[i] * fe_0 +
                           tk_zzz_xyzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xyzzz[i] * fz_0;

        tk_zzzz_xzzzz[i] = -6.0 * ts_zz_xzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xzzzz[i] * fe_0 + 4.0 * tk_zzz_xzzz[i] * fe_0 +
                           tk_zzz_xzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xzzzz[i] * fz_0;

        tk_zzzz_yyyyy[i] =
            -6.0 * ts_zz_yyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyyy[i] * fe_0 + tk_zzz_yyyyy[i] * pa_z[i] + 2.0 * ts_zzzz_yyyyy[i] * fz_0;

        tk_zzzz_yyyyz[i] = -6.0 * ts_zz_yyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyyz[i] * fe_0 + tk_zzz_yyyy[i] * fe_0 + tk_zzz_yyyyz[i] * pa_z[i] +
                           2.0 * ts_zzzz_yyyyz[i] * fz_0;

        tk_zzzz_yyyzz[i] = -6.0 * ts_zz_yyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyzz[i] * fe_0 + 2.0 * tk_zzz_yyyz[i] * fe_0 +
                           tk_zzz_yyyzz[i] * pa_z[i] + 2.0 * ts_zzzz_yyyzz[i] * fz_0;

        tk_zzzz_yyzzz[i] = -6.0 * ts_zz_yyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyzzz[i] * fe_0 + 3.0 * tk_zzz_yyzz[i] * fe_0 +
                           tk_zzz_yyzzz[i] * pa_z[i] + 2.0 * ts_zzzz_yyzzz[i] * fz_0;

        tk_zzzz_yzzzz[i] = -6.0 * ts_zz_yzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yzzzz[i] * fe_0 + 4.0 * tk_zzz_yzzz[i] * fe_0 +
                           tk_zzz_yzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_yzzzz[i] * fz_0;

        tk_zzzz_zzzzz[i] = -6.0 * ts_zz_zzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_zzzzz[i] * fe_0 + 5.0 * tk_zzz_zzzz[i] * fe_0 +
                           tk_zzz_zzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_zzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
