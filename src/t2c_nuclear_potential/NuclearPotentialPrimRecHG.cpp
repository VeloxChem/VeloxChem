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

#include "NuclearPotentialPrimRecHG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_hg(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_hg,
                               const size_t              idx_npot_0_fg,
                               const size_t              idx_npot_1_fg,
                               const size_t              idx_npot_0_gf,
                               const size_t              idx_npot_1_gf,
                               const size_t              idx_npot_0_gg,
                               const size_t              idx_npot_1_gg,
                               const CSimdArray<double>& factors,
                               const size_t              idx_rpa,
                               const size_t              idx_rpc,
                               const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : FG

    auto ta_xxx_xxxx_0 = pbuffer.data(idx_npot_0_fg);

    auto ta_xxx_xxxy_0 = pbuffer.data(idx_npot_0_fg + 1);

    auto ta_xxx_xxxz_0 = pbuffer.data(idx_npot_0_fg + 2);

    auto ta_xxx_xxyy_0 = pbuffer.data(idx_npot_0_fg + 3);

    auto ta_xxx_xxyz_0 = pbuffer.data(idx_npot_0_fg + 4);

    auto ta_xxx_xxzz_0 = pbuffer.data(idx_npot_0_fg + 5);

    auto ta_xxx_xyyy_0 = pbuffer.data(idx_npot_0_fg + 6);

    auto ta_xxx_xyyz_0 = pbuffer.data(idx_npot_0_fg + 7);

    auto ta_xxx_xyzz_0 = pbuffer.data(idx_npot_0_fg + 8);

    auto ta_xxx_xzzz_0 = pbuffer.data(idx_npot_0_fg + 9);

    auto ta_xxx_yyyy_0 = pbuffer.data(idx_npot_0_fg + 10);

    auto ta_xxx_yyyz_0 = pbuffer.data(idx_npot_0_fg + 11);

    auto ta_xxx_yyzz_0 = pbuffer.data(idx_npot_0_fg + 12);

    auto ta_xxx_yzzz_0 = pbuffer.data(idx_npot_0_fg + 13);

    auto ta_xxx_zzzz_0 = pbuffer.data(idx_npot_0_fg + 14);

    auto ta_xxy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 15);

    auto ta_xxy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 17);

    auto ta_xxy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 20);

    auto ta_xxy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 24);

    auto ta_xxy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 25);

    auto ta_xxy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 26);

    auto ta_xxy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 27);

    auto ta_xxy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 28);

    auto ta_xxz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 30);

    auto ta_xxz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 31);

    auto ta_xxz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 32);

    auto ta_xxz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 33);

    auto ta_xxz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 35);

    auto ta_xxz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 36);

    auto ta_xxz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 39);

    auto ta_xxz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 41);

    auto ta_xxz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 42);

    auto ta_xxz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 43);

    auto ta_xxz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 44);

    auto ta_xyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 46);

    auto ta_xyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 48);

    auto ta_xyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 49);

    auto ta_xyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 51);

    auto ta_xyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 52);

    auto ta_xyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 53);

    auto ta_xyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 55);

    auto ta_xyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 56);

    auto ta_xyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 57);

    auto ta_xyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 58);

    auto ta_xyy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 59);

    auto ta_xyz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 71);

    auto ta_xyz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 72);

    auto ta_xyz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 73);

    auto ta_xzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 77);

    auto ta_xzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 79);

    auto ta_xzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 80);

    auto ta_xzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 82);

    auto ta_xzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 83);

    auto ta_xzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 84);

    auto ta_xzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 85);

    auto ta_xzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 86);

    auto ta_xzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 87);

    auto ta_xzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 88);

    auto ta_xzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 89);

    auto ta_yyy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 90);

    auto ta_yyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 91);

    auto ta_yyy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 92);

    auto ta_yyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 93);

    auto ta_yyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 94);

    auto ta_yyy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 95);

    auto ta_yyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 96);

    auto ta_yyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 97);

    auto ta_yyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 98);

    auto ta_yyy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 99);

    auto ta_yyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 100);

    auto ta_yyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 101);

    auto ta_yyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 102);

    auto ta_yyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 103);

    auto ta_yyy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 104);

    auto ta_yyz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 106);

    auto ta_yyz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 107);

    auto ta_yyz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 108);

    auto ta_yyz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 110);

    auto ta_yyz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 111);

    auto ta_yyz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 114);

    auto ta_yyz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 115);

    auto ta_yyz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 116);

    auto ta_yyz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 117);

    auto ta_yyz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 118);

    auto ta_yyz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 119);

    auto ta_yzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 120);

    auto ta_yzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 122);

    auto ta_yzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 124);

    auto ta_yzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 125);

    auto ta_yzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 127);

    auto ta_yzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 128);

    auto ta_yzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 129);

    auto ta_yzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 130);

    auto ta_yzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 131);

    auto ta_yzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 132);

    auto ta_yzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 133);

    auto ta_yzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 134);

    auto ta_zzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 135);

    auto ta_zzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 136);

    auto ta_zzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 137);

    auto ta_zzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 138);

    auto ta_zzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 139);

    auto ta_zzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 140);

    auto ta_zzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 141);

    auto ta_zzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 142);

    auto ta_zzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 143);

    auto ta_zzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 144);

    auto ta_zzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 145);

    auto ta_zzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 146);

    auto ta_zzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 147);

    auto ta_zzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 148);

    auto ta_zzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 149);

    // Set up components of auxiliary buffer : FG

    auto ta_xxx_xxxx_1 = pbuffer.data(idx_npot_1_fg);

    auto ta_xxx_xxxy_1 = pbuffer.data(idx_npot_1_fg + 1);

    auto ta_xxx_xxxz_1 = pbuffer.data(idx_npot_1_fg + 2);

    auto ta_xxx_xxyy_1 = pbuffer.data(idx_npot_1_fg + 3);

    auto ta_xxx_xxyz_1 = pbuffer.data(idx_npot_1_fg + 4);

    auto ta_xxx_xxzz_1 = pbuffer.data(idx_npot_1_fg + 5);

    auto ta_xxx_xyyy_1 = pbuffer.data(idx_npot_1_fg + 6);

    auto ta_xxx_xyyz_1 = pbuffer.data(idx_npot_1_fg + 7);

    auto ta_xxx_xyzz_1 = pbuffer.data(idx_npot_1_fg + 8);

    auto ta_xxx_xzzz_1 = pbuffer.data(idx_npot_1_fg + 9);

    auto ta_xxx_yyyy_1 = pbuffer.data(idx_npot_1_fg + 10);

    auto ta_xxx_yyyz_1 = pbuffer.data(idx_npot_1_fg + 11);

    auto ta_xxx_yyzz_1 = pbuffer.data(idx_npot_1_fg + 12);

    auto ta_xxx_yzzz_1 = pbuffer.data(idx_npot_1_fg + 13);

    auto ta_xxx_zzzz_1 = pbuffer.data(idx_npot_1_fg + 14);

    auto ta_xxy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 15);

    auto ta_xxy_xxxz_1 = pbuffer.data(idx_npot_1_fg + 17);

    auto ta_xxy_xxzz_1 = pbuffer.data(idx_npot_1_fg + 20);

    auto ta_xxy_xzzz_1 = pbuffer.data(idx_npot_1_fg + 24);

    auto ta_xxy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 25);

    auto ta_xxy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 26);

    auto ta_xxy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 27);

    auto ta_xxy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 28);

    auto ta_xxz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 30);

    auto ta_xxz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 31);

    auto ta_xxz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 32);

    auto ta_xxz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 33);

    auto ta_xxz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 35);

    auto ta_xxz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 36);

    auto ta_xxz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 39);

    auto ta_xxz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 41);

    auto ta_xxz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 42);

    auto ta_xxz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 43);

    auto ta_xxz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 44);

    auto ta_xyy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 46);

    auto ta_xyy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 48);

    auto ta_xyy_xxyz_1 = pbuffer.data(idx_npot_1_fg + 49);

    auto ta_xyy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 51);

    auto ta_xyy_xyyz_1 = pbuffer.data(idx_npot_1_fg + 52);

    auto ta_xyy_xyzz_1 = pbuffer.data(idx_npot_1_fg + 53);

    auto ta_xyy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 55);

    auto ta_xyy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 56);

    auto ta_xyy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 57);

    auto ta_xyy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 58);

    auto ta_xyy_zzzz_1 = pbuffer.data(idx_npot_1_fg + 59);

    auto ta_xyz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 71);

    auto ta_xyz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 72);

    auto ta_xyz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 73);

    auto ta_xzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 77);

    auto ta_xzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 79);

    auto ta_xzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 80);

    auto ta_xzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 82);

    auto ta_xzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 83);

    auto ta_xzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 84);

    auto ta_xzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 85);

    auto ta_xzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 86);

    auto ta_xzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 87);

    auto ta_xzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 88);

    auto ta_xzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 89);

    auto ta_yyy_xxxx_1 = pbuffer.data(idx_npot_1_fg + 90);

    auto ta_yyy_xxxy_1 = pbuffer.data(idx_npot_1_fg + 91);

    auto ta_yyy_xxxz_1 = pbuffer.data(idx_npot_1_fg + 92);

    auto ta_yyy_xxyy_1 = pbuffer.data(idx_npot_1_fg + 93);

    auto ta_yyy_xxyz_1 = pbuffer.data(idx_npot_1_fg + 94);

    auto ta_yyy_xxzz_1 = pbuffer.data(idx_npot_1_fg + 95);

    auto ta_yyy_xyyy_1 = pbuffer.data(idx_npot_1_fg + 96);

    auto ta_yyy_xyyz_1 = pbuffer.data(idx_npot_1_fg + 97);

    auto ta_yyy_xyzz_1 = pbuffer.data(idx_npot_1_fg + 98);

    auto ta_yyy_xzzz_1 = pbuffer.data(idx_npot_1_fg + 99);

    auto ta_yyy_yyyy_1 = pbuffer.data(idx_npot_1_fg + 100);

    auto ta_yyy_yyyz_1 = pbuffer.data(idx_npot_1_fg + 101);

    auto ta_yyy_yyzz_1 = pbuffer.data(idx_npot_1_fg + 102);

    auto ta_yyy_yzzz_1 = pbuffer.data(idx_npot_1_fg + 103);

    auto ta_yyy_zzzz_1 = pbuffer.data(idx_npot_1_fg + 104);

    auto ta_yyz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 106);

    auto ta_yyz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 107);

    auto ta_yyz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 108);

    auto ta_yyz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 110);

    auto ta_yyz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 111);

    auto ta_yyz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 114);

    auto ta_yyz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 115);

    auto ta_yyz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 116);

    auto ta_yyz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 117);

    auto ta_yyz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 118);

    auto ta_yyz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 119);

    auto ta_yzz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 120);

    auto ta_yzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 122);

    auto ta_yzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 124);

    auto ta_yzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 125);

    auto ta_yzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 127);

    auto ta_yzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 128);

    auto ta_yzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 129);

    auto ta_yzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 130);

    auto ta_yzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 131);

    auto ta_yzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 132);

    auto ta_yzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 133);

    auto ta_yzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 134);

    auto ta_zzz_xxxx_1 = pbuffer.data(idx_npot_1_fg + 135);

    auto ta_zzz_xxxy_1 = pbuffer.data(idx_npot_1_fg + 136);

    auto ta_zzz_xxxz_1 = pbuffer.data(idx_npot_1_fg + 137);

    auto ta_zzz_xxyy_1 = pbuffer.data(idx_npot_1_fg + 138);

    auto ta_zzz_xxyz_1 = pbuffer.data(idx_npot_1_fg + 139);

    auto ta_zzz_xxzz_1 = pbuffer.data(idx_npot_1_fg + 140);

    auto ta_zzz_xyyy_1 = pbuffer.data(idx_npot_1_fg + 141);

    auto ta_zzz_xyyz_1 = pbuffer.data(idx_npot_1_fg + 142);

    auto ta_zzz_xyzz_1 = pbuffer.data(idx_npot_1_fg + 143);

    auto ta_zzz_xzzz_1 = pbuffer.data(idx_npot_1_fg + 144);

    auto ta_zzz_yyyy_1 = pbuffer.data(idx_npot_1_fg + 145);

    auto ta_zzz_yyyz_1 = pbuffer.data(idx_npot_1_fg + 146);

    auto ta_zzz_yyzz_1 = pbuffer.data(idx_npot_1_fg + 147);

    auto ta_zzz_yzzz_1 = pbuffer.data(idx_npot_1_fg + 148);

    auto ta_zzz_zzzz_1 = pbuffer.data(idx_npot_1_fg + 149);

    // Set up components of auxiliary buffer : GF

    auto ta_xxxx_xxx_0 = pbuffer.data(idx_npot_0_gf);

    auto ta_xxxx_xxy_0 = pbuffer.data(idx_npot_0_gf + 1);

    auto ta_xxxx_xxz_0 = pbuffer.data(idx_npot_0_gf + 2);

    auto ta_xxxx_xyy_0 = pbuffer.data(idx_npot_0_gf + 3);

    auto ta_xxxx_xyz_0 = pbuffer.data(idx_npot_0_gf + 4);

    auto ta_xxxx_xzz_0 = pbuffer.data(idx_npot_0_gf + 5);

    auto ta_xxxx_yyy_0 = pbuffer.data(idx_npot_0_gf + 6);

    auto ta_xxxx_yyz_0 = pbuffer.data(idx_npot_0_gf + 7);

    auto ta_xxxx_yzz_0 = pbuffer.data(idx_npot_0_gf + 8);

    auto ta_xxxx_zzz_0 = pbuffer.data(idx_npot_0_gf + 9);

    auto ta_xxxz_xxz_0 = pbuffer.data(idx_npot_0_gf + 22);

    auto ta_xxxz_xyz_0 = pbuffer.data(idx_npot_0_gf + 24);

    auto ta_xxxz_xzz_0 = pbuffer.data(idx_npot_0_gf + 25);

    auto ta_xxyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 31);

    auto ta_xxyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 33);

    auto ta_xxyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 34);

    auto ta_xxyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 36);

    auto ta_xxyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 37);

    auto ta_xxyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 38);

    auto ta_xxzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 50);

    auto ta_xxzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 51);

    auto ta_xxzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 52);

    auto ta_xxzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 53);

    auto ta_xxzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 54);

    auto ta_xxzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 55);

    auto ta_xxzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 57);

    auto ta_xxzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 58);

    auto ta_xxzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 59);

    auto ta_xyyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 61);

    auto ta_xyyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 63);

    auto ta_xyyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 64);

    auto ta_xyyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 66);

    auto ta_xyyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 67);

    auto ta_xyyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 68);

    auto ta_xzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 92);

    auto ta_xzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 94);

    auto ta_xzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 95);

    auto ta_xzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 97);

    auto ta_xzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 98);

    auto ta_xzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 99);

    auto ta_yyyy_xxx_0 = pbuffer.data(idx_npot_0_gf + 100);

    auto ta_yyyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 101);

    auto ta_yyyy_xxz_0 = pbuffer.data(idx_npot_0_gf + 102);

    auto ta_yyyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 103);

    auto ta_yyyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 104);

    auto ta_yyyy_xzz_0 = pbuffer.data(idx_npot_0_gf + 105);

    auto ta_yyyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 106);

    auto ta_yyyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 107);

    auto ta_yyyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 108);

    auto ta_yyyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 109);

    auto ta_yyyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 112);

    auto ta_yyyz_xyz_0 = pbuffer.data(idx_npot_0_gf + 114);

    auto ta_yyyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 115);

    auto ta_yyyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 117);

    auto ta_yyyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 118);

    auto ta_yyyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 119);

    auto ta_yyzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 120);

    auto ta_yyzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 121);

    auto ta_yyzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 122);

    auto ta_yyzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 123);

    auto ta_yyzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 124);

    auto ta_yyzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 125);

    auto ta_yyzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 126);

    auto ta_yyzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 127);

    auto ta_yyzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 128);

    auto ta_yyzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 129);

    auto ta_yzzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 131);

    auto ta_yzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 132);

    auto ta_yzzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 133);

    auto ta_yzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 134);

    auto ta_yzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 135);

    auto ta_yzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 136);

    auto ta_yzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 137);

    auto ta_yzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 138);

    auto ta_yzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 139);

    auto ta_zzzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 140);

    auto ta_zzzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 141);

    auto ta_zzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 142);

    auto ta_zzzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 143);

    auto ta_zzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 144);

    auto ta_zzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 145);

    auto ta_zzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 146);

    auto ta_zzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 147);

    auto ta_zzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 148);

    auto ta_zzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 149);

    // Set up components of auxiliary buffer : GF

    auto ta_xxxx_xxx_1 = pbuffer.data(idx_npot_1_gf);

    auto ta_xxxx_xxy_1 = pbuffer.data(idx_npot_1_gf + 1);

    auto ta_xxxx_xxz_1 = pbuffer.data(idx_npot_1_gf + 2);

    auto ta_xxxx_xyy_1 = pbuffer.data(idx_npot_1_gf + 3);

    auto ta_xxxx_xyz_1 = pbuffer.data(idx_npot_1_gf + 4);

    auto ta_xxxx_xzz_1 = pbuffer.data(idx_npot_1_gf + 5);

    auto ta_xxxx_yyy_1 = pbuffer.data(idx_npot_1_gf + 6);

    auto ta_xxxx_yyz_1 = pbuffer.data(idx_npot_1_gf + 7);

    auto ta_xxxx_yzz_1 = pbuffer.data(idx_npot_1_gf + 8);

    auto ta_xxxx_zzz_1 = pbuffer.data(idx_npot_1_gf + 9);

    auto ta_xxxz_xxz_1 = pbuffer.data(idx_npot_1_gf + 22);

    auto ta_xxxz_xyz_1 = pbuffer.data(idx_npot_1_gf + 24);

    auto ta_xxxz_xzz_1 = pbuffer.data(idx_npot_1_gf + 25);

    auto ta_xxyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 31);

    auto ta_xxyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 33);

    auto ta_xxyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 34);

    auto ta_xxyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 36);

    auto ta_xxyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 37);

    auto ta_xxyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 38);

    auto ta_xxzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 50);

    auto ta_xxzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 51);

    auto ta_xxzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 52);

    auto ta_xxzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 53);

    auto ta_xxzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 54);

    auto ta_xxzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 55);

    auto ta_xxzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 57);

    auto ta_xxzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 58);

    auto ta_xxzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 59);

    auto ta_xyyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 61);

    auto ta_xyyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 63);

    auto ta_xyyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 64);

    auto ta_xyyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 66);

    auto ta_xyyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 67);

    auto ta_xyyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 68);

    auto ta_xzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 92);

    auto ta_xzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 94);

    auto ta_xzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 95);

    auto ta_xzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 97);

    auto ta_xzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 98);

    auto ta_xzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 99);

    auto ta_yyyy_xxx_1 = pbuffer.data(idx_npot_1_gf + 100);

    auto ta_yyyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 101);

    auto ta_yyyy_xxz_1 = pbuffer.data(idx_npot_1_gf + 102);

    auto ta_yyyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 103);

    auto ta_yyyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 104);

    auto ta_yyyy_xzz_1 = pbuffer.data(idx_npot_1_gf + 105);

    auto ta_yyyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 106);

    auto ta_yyyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 107);

    auto ta_yyyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 108);

    auto ta_yyyy_zzz_1 = pbuffer.data(idx_npot_1_gf + 109);

    auto ta_yyyz_xxz_1 = pbuffer.data(idx_npot_1_gf + 112);

    auto ta_yyyz_xyz_1 = pbuffer.data(idx_npot_1_gf + 114);

    auto ta_yyyz_xzz_1 = pbuffer.data(idx_npot_1_gf + 115);

    auto ta_yyyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 117);

    auto ta_yyyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 118);

    auto ta_yyyz_zzz_1 = pbuffer.data(idx_npot_1_gf + 119);

    auto ta_yyzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 120);

    auto ta_yyzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 121);

    auto ta_yyzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 122);

    auto ta_yyzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 123);

    auto ta_yyzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 124);

    auto ta_yyzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 125);

    auto ta_yyzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 126);

    auto ta_yyzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 127);

    auto ta_yyzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 128);

    auto ta_yyzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 129);

    auto ta_yzzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 131);

    auto ta_yzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 132);

    auto ta_yzzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 133);

    auto ta_yzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 134);

    auto ta_yzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 135);

    auto ta_yzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 136);

    auto ta_yzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 137);

    auto ta_yzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 138);

    auto ta_yzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 139);

    auto ta_zzzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 140);

    auto ta_zzzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 141);

    auto ta_zzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 142);

    auto ta_zzzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 143);

    auto ta_zzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 144);

    auto ta_zzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 145);

    auto ta_zzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 146);

    auto ta_zzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 147);

    auto ta_zzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 148);

    auto ta_zzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 149);

    // Set up components of auxiliary buffer : GG

    auto ta_xxxx_xxxx_0 = pbuffer.data(idx_npot_0_gg);

    auto ta_xxxx_xxxy_0 = pbuffer.data(idx_npot_0_gg + 1);

    auto ta_xxxx_xxxz_0 = pbuffer.data(idx_npot_0_gg + 2);

    auto ta_xxxx_xxyy_0 = pbuffer.data(idx_npot_0_gg + 3);

    auto ta_xxxx_xxyz_0 = pbuffer.data(idx_npot_0_gg + 4);

    auto ta_xxxx_xxzz_0 = pbuffer.data(idx_npot_0_gg + 5);

    auto ta_xxxx_xyyy_0 = pbuffer.data(idx_npot_0_gg + 6);

    auto ta_xxxx_xyyz_0 = pbuffer.data(idx_npot_0_gg + 7);

    auto ta_xxxx_xyzz_0 = pbuffer.data(idx_npot_0_gg + 8);

    auto ta_xxxx_xzzz_0 = pbuffer.data(idx_npot_0_gg + 9);

    auto ta_xxxx_yyyy_0 = pbuffer.data(idx_npot_0_gg + 10);

    auto ta_xxxx_yyyz_0 = pbuffer.data(idx_npot_0_gg + 11);

    auto ta_xxxx_yyzz_0 = pbuffer.data(idx_npot_0_gg + 12);

    auto ta_xxxx_yzzz_0 = pbuffer.data(idx_npot_0_gg + 13);

    auto ta_xxxx_zzzz_0 = pbuffer.data(idx_npot_0_gg + 14);

    auto ta_xxxy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 15);

    auto ta_xxxy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 16);

    auto ta_xxxy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 17);

    auto ta_xxxy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 18);

    auto ta_xxxy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 20);

    auto ta_xxxy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 21);

    auto ta_xxxy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 24);

    auto ta_xxxy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 25);

    auto ta_xxxy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 26);

    auto ta_xxxy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 27);

    auto ta_xxxy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 28);

    auto ta_xxxz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 30);

    auto ta_xxxz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 31);

    auto ta_xxxz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 32);

    auto ta_xxxz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 33);

    auto ta_xxxz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 34);

    auto ta_xxxz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 35);

    auto ta_xxxz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 36);

    auto ta_xxxz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 37);

    auto ta_xxxz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 38);

    auto ta_xxxz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 39);

    auto ta_xxxz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 41);

    auto ta_xxxz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 42);

    auto ta_xxxz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 43);

    auto ta_xxxz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 44);

    auto ta_xxyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 45);

    auto ta_xxyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 46);

    auto ta_xxyy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 47);

    auto ta_xxyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 48);

    auto ta_xxyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 49);

    auto ta_xxyy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 50);

    auto ta_xxyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 51);

    auto ta_xxyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 52);

    auto ta_xxyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 53);

    auto ta_xxyy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 54);

    auto ta_xxyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 55);

    auto ta_xxyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 56);

    auto ta_xxyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 57);

    auto ta_xxyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 58);

    auto ta_xxyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 59);

    auto ta_xxyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 62);

    auto ta_xxyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 65);

    auto ta_xxyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 69);

    auto ta_xxyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 71);

    auto ta_xxyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 72);

    auto ta_xxyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 73);

    auto ta_xxzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 75);

    auto ta_xxzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 76);

    auto ta_xxzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 77);

    auto ta_xxzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 78);

    auto ta_xxzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 79);

    auto ta_xxzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 80);

    auto ta_xxzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 81);

    auto ta_xxzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 82);

    auto ta_xxzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 83);

    auto ta_xxzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 84);

    auto ta_xxzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 85);

    auto ta_xxzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 86);

    auto ta_xxzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 87);

    auto ta_xxzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 88);

    auto ta_xxzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 89);

    auto ta_xyyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 90);

    auto ta_xyyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 91);

    auto ta_xyyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 93);

    auto ta_xyyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 94);

    auto ta_xyyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 96);

    auto ta_xyyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 97);

    auto ta_xyyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 98);

    auto ta_xyyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 100);

    auto ta_xyyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 101);

    auto ta_xyyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 102);

    auto ta_xyyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 103);

    auto ta_xyyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 104);

    auto ta_xyyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 116);

    auto ta_xyyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 117);

    auto ta_xyyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 118);

    auto ta_xyyz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 119);

    auto ta_xyzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 130);

    auto ta_xyzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 131);

    auto ta_xyzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 132);

    auto ta_xyzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 133);

    auto ta_xzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 135);

    auto ta_xzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 137);

    auto ta_xzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 139);

    auto ta_xzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 140);

    auto ta_xzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 142);

    auto ta_xzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 143);

    auto ta_xzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 144);

    auto ta_xzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 145);

    auto ta_xzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 146);

    auto ta_xzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 147);

    auto ta_xzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 148);

    auto ta_xzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 149);

    auto ta_yyyy_xxxx_0 = pbuffer.data(idx_npot_0_gg + 150);

    auto ta_yyyy_xxxy_0 = pbuffer.data(idx_npot_0_gg + 151);

    auto ta_yyyy_xxxz_0 = pbuffer.data(idx_npot_0_gg + 152);

    auto ta_yyyy_xxyy_0 = pbuffer.data(idx_npot_0_gg + 153);

    auto ta_yyyy_xxyz_0 = pbuffer.data(idx_npot_0_gg + 154);

    auto ta_yyyy_xxzz_0 = pbuffer.data(idx_npot_0_gg + 155);

    auto ta_yyyy_xyyy_0 = pbuffer.data(idx_npot_0_gg + 156);

    auto ta_yyyy_xyyz_0 = pbuffer.data(idx_npot_0_gg + 157);

    auto ta_yyyy_xyzz_0 = pbuffer.data(idx_npot_0_gg + 158);

    auto ta_yyyy_xzzz_0 = pbuffer.data(idx_npot_0_gg + 159);

    auto ta_yyyy_yyyy_0 = pbuffer.data(idx_npot_0_gg + 160);

    auto ta_yyyy_yyyz_0 = pbuffer.data(idx_npot_0_gg + 161);

    auto ta_yyyy_yyzz_0 = pbuffer.data(idx_npot_0_gg + 162);

    auto ta_yyyy_yzzz_0 = pbuffer.data(idx_npot_0_gg + 163);

    auto ta_yyyy_zzzz_0 = pbuffer.data(idx_npot_0_gg + 164);

    auto ta_yyyz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 166);

    auto ta_yyyz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 167);

    auto ta_yyyz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 168);

    auto ta_yyyz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 169);

    auto ta_yyyz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 170);

    auto ta_yyyz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 171);

    auto ta_yyyz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 172);

    auto ta_yyyz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 173);

    auto ta_yyyz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 174);

    auto ta_yyyz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 175);

    auto ta_yyyz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 176);

    auto ta_yyyz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 177);

    auto ta_yyyz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 178);

    auto ta_yyyz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 179);

    auto ta_yyzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 180);

    auto ta_yyzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 181);

    auto ta_yyzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 182);

    auto ta_yyzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 183);

    auto ta_yyzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 184);

    auto ta_yyzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 185);

    auto ta_yyzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 186);

    auto ta_yyzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 187);

    auto ta_yyzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 188);

    auto ta_yyzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 189);

    auto ta_yyzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 190);

    auto ta_yyzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 191);

    auto ta_yyzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 192);

    auto ta_yyzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 193);

    auto ta_yyzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 194);

    auto ta_yzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 195);

    auto ta_yzzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 196);

    auto ta_yzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 197);

    auto ta_yzzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 198);

    auto ta_yzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 199);

    auto ta_yzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 200);

    auto ta_yzzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 201);

    auto ta_yzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 202);

    auto ta_yzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 203);

    auto ta_yzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 204);

    auto ta_yzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 205);

    auto ta_yzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 206);

    auto ta_yzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 207);

    auto ta_yzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 208);

    auto ta_yzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 209);

    auto ta_zzzz_xxxx_0 = pbuffer.data(idx_npot_0_gg + 210);

    auto ta_zzzz_xxxy_0 = pbuffer.data(idx_npot_0_gg + 211);

    auto ta_zzzz_xxxz_0 = pbuffer.data(idx_npot_0_gg + 212);

    auto ta_zzzz_xxyy_0 = pbuffer.data(idx_npot_0_gg + 213);

    auto ta_zzzz_xxyz_0 = pbuffer.data(idx_npot_0_gg + 214);

    auto ta_zzzz_xxzz_0 = pbuffer.data(idx_npot_0_gg + 215);

    auto ta_zzzz_xyyy_0 = pbuffer.data(idx_npot_0_gg + 216);

    auto ta_zzzz_xyyz_0 = pbuffer.data(idx_npot_0_gg + 217);

    auto ta_zzzz_xyzz_0 = pbuffer.data(idx_npot_0_gg + 218);

    auto ta_zzzz_xzzz_0 = pbuffer.data(idx_npot_0_gg + 219);

    auto ta_zzzz_yyyy_0 = pbuffer.data(idx_npot_0_gg + 220);

    auto ta_zzzz_yyyz_0 = pbuffer.data(idx_npot_0_gg + 221);

    auto ta_zzzz_yyzz_0 = pbuffer.data(idx_npot_0_gg + 222);

    auto ta_zzzz_yzzz_0 = pbuffer.data(idx_npot_0_gg + 223);

    auto ta_zzzz_zzzz_0 = pbuffer.data(idx_npot_0_gg + 224);

    // Set up components of auxiliary buffer : GG

    auto ta_xxxx_xxxx_1 = pbuffer.data(idx_npot_1_gg);

    auto ta_xxxx_xxxy_1 = pbuffer.data(idx_npot_1_gg + 1);

    auto ta_xxxx_xxxz_1 = pbuffer.data(idx_npot_1_gg + 2);

    auto ta_xxxx_xxyy_1 = pbuffer.data(idx_npot_1_gg + 3);

    auto ta_xxxx_xxyz_1 = pbuffer.data(idx_npot_1_gg + 4);

    auto ta_xxxx_xxzz_1 = pbuffer.data(idx_npot_1_gg + 5);

    auto ta_xxxx_xyyy_1 = pbuffer.data(idx_npot_1_gg + 6);

    auto ta_xxxx_xyyz_1 = pbuffer.data(idx_npot_1_gg + 7);

    auto ta_xxxx_xyzz_1 = pbuffer.data(idx_npot_1_gg + 8);

    auto ta_xxxx_xzzz_1 = pbuffer.data(idx_npot_1_gg + 9);

    auto ta_xxxx_yyyy_1 = pbuffer.data(idx_npot_1_gg + 10);

    auto ta_xxxx_yyyz_1 = pbuffer.data(idx_npot_1_gg + 11);

    auto ta_xxxx_yyzz_1 = pbuffer.data(idx_npot_1_gg + 12);

    auto ta_xxxx_yzzz_1 = pbuffer.data(idx_npot_1_gg + 13);

    auto ta_xxxx_zzzz_1 = pbuffer.data(idx_npot_1_gg + 14);

    auto ta_xxxy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 15);

    auto ta_xxxy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 16);

    auto ta_xxxy_xxxz_1 = pbuffer.data(idx_npot_1_gg + 17);

    auto ta_xxxy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 18);

    auto ta_xxxy_xxzz_1 = pbuffer.data(idx_npot_1_gg + 20);

    auto ta_xxxy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 21);

    auto ta_xxxy_xzzz_1 = pbuffer.data(idx_npot_1_gg + 24);

    auto ta_xxxy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 25);

    auto ta_xxxy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 26);

    auto ta_xxxy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 27);

    auto ta_xxxy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 28);

    auto ta_xxxz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 30);

    auto ta_xxxz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 31);

    auto ta_xxxz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 32);

    auto ta_xxxz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 33);

    auto ta_xxxz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 34);

    auto ta_xxxz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 35);

    auto ta_xxxz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 36);

    auto ta_xxxz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 37);

    auto ta_xxxz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 38);

    auto ta_xxxz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 39);

    auto ta_xxxz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 41);

    auto ta_xxxz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 42);

    auto ta_xxxz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 43);

    auto ta_xxxz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 44);

    auto ta_xxyy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 45);

    auto ta_xxyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 46);

    auto ta_xxyy_xxxz_1 = pbuffer.data(idx_npot_1_gg + 47);

    auto ta_xxyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 48);

    auto ta_xxyy_xxyz_1 = pbuffer.data(idx_npot_1_gg + 49);

    auto ta_xxyy_xxzz_1 = pbuffer.data(idx_npot_1_gg + 50);

    auto ta_xxyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 51);

    auto ta_xxyy_xyyz_1 = pbuffer.data(idx_npot_1_gg + 52);

    auto ta_xxyy_xyzz_1 = pbuffer.data(idx_npot_1_gg + 53);

    auto ta_xxyy_xzzz_1 = pbuffer.data(idx_npot_1_gg + 54);

    auto ta_xxyy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 55);

    auto ta_xxyy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 56);

    auto ta_xxyy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 57);

    auto ta_xxyy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 58);

    auto ta_xxyy_zzzz_1 = pbuffer.data(idx_npot_1_gg + 59);

    auto ta_xxyz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 62);

    auto ta_xxyz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 65);

    auto ta_xxyz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 69);

    auto ta_xxyz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 71);

    auto ta_xxyz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 72);

    auto ta_xxyz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 73);

    auto ta_xxzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 75);

    auto ta_xxzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 76);

    auto ta_xxzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 77);

    auto ta_xxzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 78);

    auto ta_xxzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 79);

    auto ta_xxzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 80);

    auto ta_xxzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 81);

    auto ta_xxzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 82);

    auto ta_xxzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 83);

    auto ta_xxzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 84);

    auto ta_xxzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 85);

    auto ta_xxzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 86);

    auto ta_xxzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 87);

    auto ta_xxzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 88);

    auto ta_xxzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 89);

    auto ta_xyyy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 90);

    auto ta_xyyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 91);

    auto ta_xyyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 93);

    auto ta_xyyy_xxyz_1 = pbuffer.data(idx_npot_1_gg + 94);

    auto ta_xyyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 96);

    auto ta_xyyy_xyyz_1 = pbuffer.data(idx_npot_1_gg + 97);

    auto ta_xyyy_xyzz_1 = pbuffer.data(idx_npot_1_gg + 98);

    auto ta_xyyy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 100);

    auto ta_xyyy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 101);

    auto ta_xyyy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 102);

    auto ta_xyyy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 103);

    auto ta_xyyy_zzzz_1 = pbuffer.data(idx_npot_1_gg + 104);

    auto ta_xyyz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 116);

    auto ta_xyyz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 117);

    auto ta_xyyz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 118);

    auto ta_xyyz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 119);

    auto ta_xyzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 130);

    auto ta_xyzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 131);

    auto ta_xyzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 132);

    auto ta_xyzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 133);

    auto ta_xzzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 135);

    auto ta_xzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 137);

    auto ta_xzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 139);

    auto ta_xzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 140);

    auto ta_xzzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 142);

    auto ta_xzzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 143);

    auto ta_xzzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 144);

    auto ta_xzzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 145);

    auto ta_xzzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 146);

    auto ta_xzzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 147);

    auto ta_xzzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 148);

    auto ta_xzzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 149);

    auto ta_yyyy_xxxx_1 = pbuffer.data(idx_npot_1_gg + 150);

    auto ta_yyyy_xxxy_1 = pbuffer.data(idx_npot_1_gg + 151);

    auto ta_yyyy_xxxz_1 = pbuffer.data(idx_npot_1_gg + 152);

    auto ta_yyyy_xxyy_1 = pbuffer.data(idx_npot_1_gg + 153);

    auto ta_yyyy_xxyz_1 = pbuffer.data(idx_npot_1_gg + 154);

    auto ta_yyyy_xxzz_1 = pbuffer.data(idx_npot_1_gg + 155);

    auto ta_yyyy_xyyy_1 = pbuffer.data(idx_npot_1_gg + 156);

    auto ta_yyyy_xyyz_1 = pbuffer.data(idx_npot_1_gg + 157);

    auto ta_yyyy_xyzz_1 = pbuffer.data(idx_npot_1_gg + 158);

    auto ta_yyyy_xzzz_1 = pbuffer.data(idx_npot_1_gg + 159);

    auto ta_yyyy_yyyy_1 = pbuffer.data(idx_npot_1_gg + 160);

    auto ta_yyyy_yyyz_1 = pbuffer.data(idx_npot_1_gg + 161);

    auto ta_yyyy_yyzz_1 = pbuffer.data(idx_npot_1_gg + 162);

    auto ta_yyyy_yzzz_1 = pbuffer.data(idx_npot_1_gg + 163);

    auto ta_yyyy_zzzz_1 = pbuffer.data(idx_npot_1_gg + 164);

    auto ta_yyyz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 166);

    auto ta_yyyz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 167);

    auto ta_yyyz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 168);

    auto ta_yyyz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 169);

    auto ta_yyyz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 170);

    auto ta_yyyz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 171);

    auto ta_yyyz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 172);

    auto ta_yyyz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 173);

    auto ta_yyyz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 174);

    auto ta_yyyz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 175);

    auto ta_yyyz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 176);

    auto ta_yyyz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 177);

    auto ta_yyyz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 178);

    auto ta_yyyz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 179);

    auto ta_yyzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 180);

    auto ta_yyzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 181);

    auto ta_yyzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 182);

    auto ta_yyzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 183);

    auto ta_yyzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 184);

    auto ta_yyzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 185);

    auto ta_yyzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 186);

    auto ta_yyzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 187);

    auto ta_yyzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 188);

    auto ta_yyzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 189);

    auto ta_yyzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 190);

    auto ta_yyzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 191);

    auto ta_yyzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 192);

    auto ta_yyzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 193);

    auto ta_yyzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 194);

    auto ta_yzzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 195);

    auto ta_yzzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 196);

    auto ta_yzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 197);

    auto ta_yzzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 198);

    auto ta_yzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 199);

    auto ta_yzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 200);

    auto ta_yzzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 201);

    auto ta_yzzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 202);

    auto ta_yzzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 203);

    auto ta_yzzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 204);

    auto ta_yzzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 205);

    auto ta_yzzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 206);

    auto ta_yzzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 207);

    auto ta_yzzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 208);

    auto ta_yzzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 209);

    auto ta_zzzz_xxxx_1 = pbuffer.data(idx_npot_1_gg + 210);

    auto ta_zzzz_xxxy_1 = pbuffer.data(idx_npot_1_gg + 211);

    auto ta_zzzz_xxxz_1 = pbuffer.data(idx_npot_1_gg + 212);

    auto ta_zzzz_xxyy_1 = pbuffer.data(idx_npot_1_gg + 213);

    auto ta_zzzz_xxyz_1 = pbuffer.data(idx_npot_1_gg + 214);

    auto ta_zzzz_xxzz_1 = pbuffer.data(idx_npot_1_gg + 215);

    auto ta_zzzz_xyyy_1 = pbuffer.data(idx_npot_1_gg + 216);

    auto ta_zzzz_xyyz_1 = pbuffer.data(idx_npot_1_gg + 217);

    auto ta_zzzz_xyzz_1 = pbuffer.data(idx_npot_1_gg + 218);

    auto ta_zzzz_xzzz_1 = pbuffer.data(idx_npot_1_gg + 219);

    auto ta_zzzz_yyyy_1 = pbuffer.data(idx_npot_1_gg + 220);

    auto ta_zzzz_yyyz_1 = pbuffer.data(idx_npot_1_gg + 221);

    auto ta_zzzz_yyzz_1 = pbuffer.data(idx_npot_1_gg + 222);

    auto ta_zzzz_yzzz_1 = pbuffer.data(idx_npot_1_gg + 223);

    auto ta_zzzz_zzzz_1 = pbuffer.data(idx_npot_1_gg + 224);

    // Set up 0-15 components of targeted buffer : HG

    auto ta_xxxxx_xxxx_0 = pbuffer.data(idx_npot_0_hg);

    auto ta_xxxxx_xxxy_0 = pbuffer.data(idx_npot_0_hg + 1);

    auto ta_xxxxx_xxxz_0 = pbuffer.data(idx_npot_0_hg + 2);

    auto ta_xxxxx_xxyy_0 = pbuffer.data(idx_npot_0_hg + 3);

    auto ta_xxxxx_xxyz_0 = pbuffer.data(idx_npot_0_hg + 4);

    auto ta_xxxxx_xxzz_0 = pbuffer.data(idx_npot_0_hg + 5);

    auto ta_xxxxx_xyyy_0 = pbuffer.data(idx_npot_0_hg + 6);

    auto ta_xxxxx_xyyz_0 = pbuffer.data(idx_npot_0_hg + 7);

    auto ta_xxxxx_xyzz_0 = pbuffer.data(idx_npot_0_hg + 8);

    auto ta_xxxxx_xzzz_0 = pbuffer.data(idx_npot_0_hg + 9);

    auto ta_xxxxx_yyyy_0 = pbuffer.data(idx_npot_0_hg + 10);

    auto ta_xxxxx_yyyz_0 = pbuffer.data(idx_npot_0_hg + 11);

    auto ta_xxxxx_yyzz_0 = pbuffer.data(idx_npot_0_hg + 12);

    auto ta_xxxxx_yzzz_0 = pbuffer.data(idx_npot_0_hg + 13);

    auto ta_xxxxx_zzzz_0 = pbuffer.data(idx_npot_0_hg + 14);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xxx_xxxx_0,   \
                             ta_xxx_xxxx_1,   \
                             ta_xxx_xxxy_0,   \
                             ta_xxx_xxxy_1,   \
                             ta_xxx_xxxz_0,   \
                             ta_xxx_xxxz_1,   \
                             ta_xxx_xxyy_0,   \
                             ta_xxx_xxyy_1,   \
                             ta_xxx_xxyz_0,   \
                             ta_xxx_xxyz_1,   \
                             ta_xxx_xxzz_0,   \
                             ta_xxx_xxzz_1,   \
                             ta_xxx_xyyy_0,   \
                             ta_xxx_xyyy_1,   \
                             ta_xxx_xyyz_0,   \
                             ta_xxx_xyyz_1,   \
                             ta_xxx_xyzz_0,   \
                             ta_xxx_xyzz_1,   \
                             ta_xxx_xzzz_0,   \
                             ta_xxx_xzzz_1,   \
                             ta_xxx_yyyy_0,   \
                             ta_xxx_yyyy_1,   \
                             ta_xxx_yyyz_0,   \
                             ta_xxx_yyyz_1,   \
                             ta_xxx_yyzz_0,   \
                             ta_xxx_yyzz_1,   \
                             ta_xxx_yzzz_0,   \
                             ta_xxx_yzzz_1,   \
                             ta_xxx_zzzz_0,   \
                             ta_xxx_zzzz_1,   \
                             ta_xxxx_xxx_0,   \
                             ta_xxxx_xxx_1,   \
                             ta_xxxx_xxxx_0,  \
                             ta_xxxx_xxxx_1,  \
                             ta_xxxx_xxxy_0,  \
                             ta_xxxx_xxxy_1,  \
                             ta_xxxx_xxxz_0,  \
                             ta_xxxx_xxxz_1,  \
                             ta_xxxx_xxy_0,   \
                             ta_xxxx_xxy_1,   \
                             ta_xxxx_xxyy_0,  \
                             ta_xxxx_xxyy_1,  \
                             ta_xxxx_xxyz_0,  \
                             ta_xxxx_xxyz_1,  \
                             ta_xxxx_xxz_0,   \
                             ta_xxxx_xxz_1,   \
                             ta_xxxx_xxzz_0,  \
                             ta_xxxx_xxzz_1,  \
                             ta_xxxx_xyy_0,   \
                             ta_xxxx_xyy_1,   \
                             ta_xxxx_xyyy_0,  \
                             ta_xxxx_xyyy_1,  \
                             ta_xxxx_xyyz_0,  \
                             ta_xxxx_xyyz_1,  \
                             ta_xxxx_xyz_0,   \
                             ta_xxxx_xyz_1,   \
                             ta_xxxx_xyzz_0,  \
                             ta_xxxx_xyzz_1,  \
                             ta_xxxx_xzz_0,   \
                             ta_xxxx_xzz_1,   \
                             ta_xxxx_xzzz_0,  \
                             ta_xxxx_xzzz_1,  \
                             ta_xxxx_yyy_0,   \
                             ta_xxxx_yyy_1,   \
                             ta_xxxx_yyyy_0,  \
                             ta_xxxx_yyyy_1,  \
                             ta_xxxx_yyyz_0,  \
                             ta_xxxx_yyyz_1,  \
                             ta_xxxx_yyz_0,   \
                             ta_xxxx_yyz_1,   \
                             ta_xxxx_yyzz_0,  \
                             ta_xxxx_yyzz_1,  \
                             ta_xxxx_yzz_0,   \
                             ta_xxxx_yzz_1,   \
                             ta_xxxx_yzzz_0,  \
                             ta_xxxx_yzzz_1,  \
                             ta_xxxx_zzz_0,   \
                             ta_xxxx_zzz_1,   \
                             ta_xxxx_zzzz_0,  \
                             ta_xxxx_zzzz_1,  \
                             ta_xxxxx_xxxx_0, \
                             ta_xxxxx_xxxy_0, \
                             ta_xxxxx_xxxz_0, \
                             ta_xxxxx_xxyy_0, \
                             ta_xxxxx_xxyz_0, \
                             ta_xxxxx_xxzz_0, \
                             ta_xxxxx_xyyy_0, \
                             ta_xxxxx_xyyz_0, \
                             ta_xxxxx_xyzz_0, \
                             ta_xxxxx_xzzz_0, \
                             ta_xxxxx_yyyy_0, \
                             ta_xxxxx_yyyz_0, \
                             ta_xxxxx_yyzz_0, \
                             ta_xxxxx_yzzz_0, \
                             ta_xxxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxx_xxxx_0[i] = 4.0 * ta_xxx_xxxx_0[i] * fe_0 - 4.0 * ta_xxx_xxxx_1[i] * fe_0 + 4.0 * ta_xxxx_xxx_0[i] * fe_0 -
                             4.0 * ta_xxxx_xxx_1[i] * fe_0 + ta_xxxx_xxxx_0[i] * pa_x[i] - ta_xxxx_xxxx_1[i] * pc_x[i];

        ta_xxxxx_xxxy_0[i] = 4.0 * ta_xxx_xxxy_0[i] * fe_0 - 4.0 * ta_xxx_xxxy_1[i] * fe_0 + 3.0 * ta_xxxx_xxy_0[i] * fe_0 -
                             3.0 * ta_xxxx_xxy_1[i] * fe_0 + ta_xxxx_xxxy_0[i] * pa_x[i] - ta_xxxx_xxxy_1[i] * pc_x[i];

        ta_xxxxx_xxxz_0[i] = 4.0 * ta_xxx_xxxz_0[i] * fe_0 - 4.0 * ta_xxx_xxxz_1[i] * fe_0 + 3.0 * ta_xxxx_xxz_0[i] * fe_0 -
                             3.0 * ta_xxxx_xxz_1[i] * fe_0 + ta_xxxx_xxxz_0[i] * pa_x[i] - ta_xxxx_xxxz_1[i] * pc_x[i];

        ta_xxxxx_xxyy_0[i] = 4.0 * ta_xxx_xxyy_0[i] * fe_0 - 4.0 * ta_xxx_xxyy_1[i] * fe_0 + 2.0 * ta_xxxx_xyy_0[i] * fe_0 -
                             2.0 * ta_xxxx_xyy_1[i] * fe_0 + ta_xxxx_xxyy_0[i] * pa_x[i] - ta_xxxx_xxyy_1[i] * pc_x[i];

        ta_xxxxx_xxyz_0[i] = 4.0 * ta_xxx_xxyz_0[i] * fe_0 - 4.0 * ta_xxx_xxyz_1[i] * fe_0 + 2.0 * ta_xxxx_xyz_0[i] * fe_0 -
                             2.0 * ta_xxxx_xyz_1[i] * fe_0 + ta_xxxx_xxyz_0[i] * pa_x[i] - ta_xxxx_xxyz_1[i] * pc_x[i];

        ta_xxxxx_xxzz_0[i] = 4.0 * ta_xxx_xxzz_0[i] * fe_0 - 4.0 * ta_xxx_xxzz_1[i] * fe_0 + 2.0 * ta_xxxx_xzz_0[i] * fe_0 -
                             2.0 * ta_xxxx_xzz_1[i] * fe_0 + ta_xxxx_xxzz_0[i] * pa_x[i] - ta_xxxx_xxzz_1[i] * pc_x[i];

        ta_xxxxx_xyyy_0[i] = 4.0 * ta_xxx_xyyy_0[i] * fe_0 - 4.0 * ta_xxx_xyyy_1[i] * fe_0 + ta_xxxx_yyy_0[i] * fe_0 - ta_xxxx_yyy_1[i] * fe_0 +
                             ta_xxxx_xyyy_0[i] * pa_x[i] - ta_xxxx_xyyy_1[i] * pc_x[i];

        ta_xxxxx_xyyz_0[i] = 4.0 * ta_xxx_xyyz_0[i] * fe_0 - 4.0 * ta_xxx_xyyz_1[i] * fe_0 + ta_xxxx_yyz_0[i] * fe_0 - ta_xxxx_yyz_1[i] * fe_0 +
                             ta_xxxx_xyyz_0[i] * pa_x[i] - ta_xxxx_xyyz_1[i] * pc_x[i];

        ta_xxxxx_xyzz_0[i] = 4.0 * ta_xxx_xyzz_0[i] * fe_0 - 4.0 * ta_xxx_xyzz_1[i] * fe_0 + ta_xxxx_yzz_0[i] * fe_0 - ta_xxxx_yzz_1[i] * fe_0 +
                             ta_xxxx_xyzz_0[i] * pa_x[i] - ta_xxxx_xyzz_1[i] * pc_x[i];

        ta_xxxxx_xzzz_0[i] = 4.0 * ta_xxx_xzzz_0[i] * fe_0 - 4.0 * ta_xxx_xzzz_1[i] * fe_0 + ta_xxxx_zzz_0[i] * fe_0 - ta_xxxx_zzz_1[i] * fe_0 +
                             ta_xxxx_xzzz_0[i] * pa_x[i] - ta_xxxx_xzzz_1[i] * pc_x[i];

        ta_xxxxx_yyyy_0[i] =
            4.0 * ta_xxx_yyyy_0[i] * fe_0 - 4.0 * ta_xxx_yyyy_1[i] * fe_0 + ta_xxxx_yyyy_0[i] * pa_x[i] - ta_xxxx_yyyy_1[i] * pc_x[i];

        ta_xxxxx_yyyz_0[i] =
            4.0 * ta_xxx_yyyz_0[i] * fe_0 - 4.0 * ta_xxx_yyyz_1[i] * fe_0 + ta_xxxx_yyyz_0[i] * pa_x[i] - ta_xxxx_yyyz_1[i] * pc_x[i];

        ta_xxxxx_yyzz_0[i] =
            4.0 * ta_xxx_yyzz_0[i] * fe_0 - 4.0 * ta_xxx_yyzz_1[i] * fe_0 + ta_xxxx_yyzz_0[i] * pa_x[i] - ta_xxxx_yyzz_1[i] * pc_x[i];

        ta_xxxxx_yzzz_0[i] =
            4.0 * ta_xxx_yzzz_0[i] * fe_0 - 4.0 * ta_xxx_yzzz_1[i] * fe_0 + ta_xxxx_yzzz_0[i] * pa_x[i] - ta_xxxx_yzzz_1[i] * pc_x[i];

        ta_xxxxx_zzzz_0[i] =
            4.0 * ta_xxx_zzzz_0[i] * fe_0 - 4.0 * ta_xxx_zzzz_1[i] * fe_0 + ta_xxxx_zzzz_0[i] * pa_x[i] - ta_xxxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : HG

    auto ta_xxxxy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 15);

    auto ta_xxxxy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 16);

    auto ta_xxxxy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 17);

    auto ta_xxxxy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 18);

    auto ta_xxxxy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 19);

    auto ta_xxxxy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 20);

    auto ta_xxxxy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 21);

    auto ta_xxxxy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 22);

    auto ta_xxxxy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 23);

    auto ta_xxxxy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 24);

    auto ta_xxxxy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 25);

    auto ta_xxxxy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 26);

    auto ta_xxxxy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 27);

    auto ta_xxxxy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 28);

    auto ta_xxxxy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 29);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxxx_xxx_0,   \
                             ta_xxxx_xxx_1,   \
                             ta_xxxx_xxxx_0,  \
                             ta_xxxx_xxxx_1,  \
                             ta_xxxx_xxxy_0,  \
                             ta_xxxx_xxxy_1,  \
                             ta_xxxx_xxxz_0,  \
                             ta_xxxx_xxxz_1,  \
                             ta_xxxx_xxy_0,   \
                             ta_xxxx_xxy_1,   \
                             ta_xxxx_xxyy_0,  \
                             ta_xxxx_xxyy_1,  \
                             ta_xxxx_xxyz_0,  \
                             ta_xxxx_xxyz_1,  \
                             ta_xxxx_xxz_0,   \
                             ta_xxxx_xxz_1,   \
                             ta_xxxx_xxzz_0,  \
                             ta_xxxx_xxzz_1,  \
                             ta_xxxx_xyy_0,   \
                             ta_xxxx_xyy_1,   \
                             ta_xxxx_xyyy_0,  \
                             ta_xxxx_xyyy_1,  \
                             ta_xxxx_xyyz_0,  \
                             ta_xxxx_xyyz_1,  \
                             ta_xxxx_xyz_0,   \
                             ta_xxxx_xyz_1,   \
                             ta_xxxx_xyzz_0,  \
                             ta_xxxx_xyzz_1,  \
                             ta_xxxx_xzz_0,   \
                             ta_xxxx_xzz_1,   \
                             ta_xxxx_xzzz_0,  \
                             ta_xxxx_xzzz_1,  \
                             ta_xxxx_zzzz_0,  \
                             ta_xxxx_zzzz_1,  \
                             ta_xxxxy_xxxx_0, \
                             ta_xxxxy_xxxy_0, \
                             ta_xxxxy_xxxz_0, \
                             ta_xxxxy_xxyy_0, \
                             ta_xxxxy_xxyz_0, \
                             ta_xxxxy_xxzz_0, \
                             ta_xxxxy_xyyy_0, \
                             ta_xxxxy_xyyz_0, \
                             ta_xxxxy_xyzz_0, \
                             ta_xxxxy_xzzz_0, \
                             ta_xxxxy_yyyy_0, \
                             ta_xxxxy_yyyz_0, \
                             ta_xxxxy_yyzz_0, \
                             ta_xxxxy_yzzz_0, \
                             ta_xxxxy_zzzz_0, \
                             ta_xxxy_yyyy_0,  \
                             ta_xxxy_yyyy_1,  \
                             ta_xxxy_yyyz_0,  \
                             ta_xxxy_yyyz_1,  \
                             ta_xxxy_yyzz_0,  \
                             ta_xxxy_yyzz_1,  \
                             ta_xxxy_yzzz_0,  \
                             ta_xxxy_yzzz_1,  \
                             ta_xxy_yyyy_0,   \
                             ta_xxy_yyyy_1,   \
                             ta_xxy_yyyz_0,   \
                             ta_xxy_yyyz_1,   \
                             ta_xxy_yyzz_0,   \
                             ta_xxy_yyzz_1,   \
                             ta_xxy_yzzz_0,   \
                             ta_xxy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxy_xxxx_0[i] = ta_xxxx_xxxx_0[i] * pa_y[i] - ta_xxxx_xxxx_1[i] * pc_y[i];

        ta_xxxxy_xxxy_0[i] = ta_xxxx_xxx_0[i] * fe_0 - ta_xxxx_xxx_1[i] * fe_0 + ta_xxxx_xxxy_0[i] * pa_y[i] - ta_xxxx_xxxy_1[i] * pc_y[i];

        ta_xxxxy_xxxz_0[i] = ta_xxxx_xxxz_0[i] * pa_y[i] - ta_xxxx_xxxz_1[i] * pc_y[i];

        ta_xxxxy_xxyy_0[i] =
            2.0 * ta_xxxx_xxy_0[i] * fe_0 - 2.0 * ta_xxxx_xxy_1[i] * fe_0 + ta_xxxx_xxyy_0[i] * pa_y[i] - ta_xxxx_xxyy_1[i] * pc_y[i];

        ta_xxxxy_xxyz_0[i] = ta_xxxx_xxz_0[i] * fe_0 - ta_xxxx_xxz_1[i] * fe_0 + ta_xxxx_xxyz_0[i] * pa_y[i] - ta_xxxx_xxyz_1[i] * pc_y[i];

        ta_xxxxy_xxzz_0[i] = ta_xxxx_xxzz_0[i] * pa_y[i] - ta_xxxx_xxzz_1[i] * pc_y[i];

        ta_xxxxy_xyyy_0[i] =
            3.0 * ta_xxxx_xyy_0[i] * fe_0 - 3.0 * ta_xxxx_xyy_1[i] * fe_0 + ta_xxxx_xyyy_0[i] * pa_y[i] - ta_xxxx_xyyy_1[i] * pc_y[i];

        ta_xxxxy_xyyz_0[i] =
            2.0 * ta_xxxx_xyz_0[i] * fe_0 - 2.0 * ta_xxxx_xyz_1[i] * fe_0 + ta_xxxx_xyyz_0[i] * pa_y[i] - ta_xxxx_xyyz_1[i] * pc_y[i];

        ta_xxxxy_xyzz_0[i] = ta_xxxx_xzz_0[i] * fe_0 - ta_xxxx_xzz_1[i] * fe_0 + ta_xxxx_xyzz_0[i] * pa_y[i] - ta_xxxx_xyzz_1[i] * pc_y[i];

        ta_xxxxy_xzzz_0[i] = ta_xxxx_xzzz_0[i] * pa_y[i] - ta_xxxx_xzzz_1[i] * pc_y[i];

        ta_xxxxy_yyyy_0[i] =
            3.0 * ta_xxy_yyyy_0[i] * fe_0 - 3.0 * ta_xxy_yyyy_1[i] * fe_0 + ta_xxxy_yyyy_0[i] * pa_x[i] - ta_xxxy_yyyy_1[i] * pc_x[i];

        ta_xxxxy_yyyz_0[i] =
            3.0 * ta_xxy_yyyz_0[i] * fe_0 - 3.0 * ta_xxy_yyyz_1[i] * fe_0 + ta_xxxy_yyyz_0[i] * pa_x[i] - ta_xxxy_yyyz_1[i] * pc_x[i];

        ta_xxxxy_yyzz_0[i] =
            3.0 * ta_xxy_yyzz_0[i] * fe_0 - 3.0 * ta_xxy_yyzz_1[i] * fe_0 + ta_xxxy_yyzz_0[i] * pa_x[i] - ta_xxxy_yyzz_1[i] * pc_x[i];

        ta_xxxxy_yzzz_0[i] =
            3.0 * ta_xxy_yzzz_0[i] * fe_0 - 3.0 * ta_xxy_yzzz_1[i] * fe_0 + ta_xxxy_yzzz_0[i] * pa_x[i] - ta_xxxy_yzzz_1[i] * pc_x[i];

        ta_xxxxy_zzzz_0[i] = ta_xxxx_zzzz_0[i] * pa_y[i] - ta_xxxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : HG

    auto ta_xxxxz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 30);

    auto ta_xxxxz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 31);

    auto ta_xxxxz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 32);

    auto ta_xxxxz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 33);

    auto ta_xxxxz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 34);

    auto ta_xxxxz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 35);

    auto ta_xxxxz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 36);

    auto ta_xxxxz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 37);

    auto ta_xxxxz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 38);

    auto ta_xxxxz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 39);

    auto ta_xxxxz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 40);

    auto ta_xxxxz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 41);

    auto ta_xxxxz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 42);

    auto ta_xxxxz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 43);

    auto ta_xxxxz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 44);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xxxx_xxx_0,   \
                             ta_xxxx_xxx_1,   \
                             ta_xxxx_xxxx_0,  \
                             ta_xxxx_xxxx_1,  \
                             ta_xxxx_xxxy_0,  \
                             ta_xxxx_xxxy_1,  \
                             ta_xxxx_xxxz_0,  \
                             ta_xxxx_xxxz_1,  \
                             ta_xxxx_xxy_0,   \
                             ta_xxxx_xxy_1,   \
                             ta_xxxx_xxyy_0,  \
                             ta_xxxx_xxyy_1,  \
                             ta_xxxx_xxyz_0,  \
                             ta_xxxx_xxyz_1,  \
                             ta_xxxx_xxz_0,   \
                             ta_xxxx_xxz_1,   \
                             ta_xxxx_xxzz_0,  \
                             ta_xxxx_xxzz_1,  \
                             ta_xxxx_xyy_0,   \
                             ta_xxxx_xyy_1,   \
                             ta_xxxx_xyyy_0,  \
                             ta_xxxx_xyyy_1,  \
                             ta_xxxx_xyyz_0,  \
                             ta_xxxx_xyyz_1,  \
                             ta_xxxx_xyz_0,   \
                             ta_xxxx_xyz_1,   \
                             ta_xxxx_xyzz_0,  \
                             ta_xxxx_xyzz_1,  \
                             ta_xxxx_xzz_0,   \
                             ta_xxxx_xzz_1,   \
                             ta_xxxx_xzzz_0,  \
                             ta_xxxx_xzzz_1,  \
                             ta_xxxx_yyyy_0,  \
                             ta_xxxx_yyyy_1,  \
                             ta_xxxxz_xxxx_0, \
                             ta_xxxxz_xxxy_0, \
                             ta_xxxxz_xxxz_0, \
                             ta_xxxxz_xxyy_0, \
                             ta_xxxxz_xxyz_0, \
                             ta_xxxxz_xxzz_0, \
                             ta_xxxxz_xyyy_0, \
                             ta_xxxxz_xyyz_0, \
                             ta_xxxxz_xyzz_0, \
                             ta_xxxxz_xzzz_0, \
                             ta_xxxxz_yyyy_0, \
                             ta_xxxxz_yyyz_0, \
                             ta_xxxxz_yyzz_0, \
                             ta_xxxxz_yzzz_0, \
                             ta_xxxxz_zzzz_0, \
                             ta_xxxz_yyyz_0,  \
                             ta_xxxz_yyyz_1,  \
                             ta_xxxz_yyzz_0,  \
                             ta_xxxz_yyzz_1,  \
                             ta_xxxz_yzzz_0,  \
                             ta_xxxz_yzzz_1,  \
                             ta_xxxz_zzzz_0,  \
                             ta_xxxz_zzzz_1,  \
                             ta_xxz_yyyz_0,   \
                             ta_xxz_yyyz_1,   \
                             ta_xxz_yyzz_0,   \
                             ta_xxz_yyzz_1,   \
                             ta_xxz_yzzz_0,   \
                             ta_xxz_yzzz_1,   \
                             ta_xxz_zzzz_0,   \
                             ta_xxz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxz_xxxx_0[i] = ta_xxxx_xxxx_0[i] * pa_z[i] - ta_xxxx_xxxx_1[i] * pc_z[i];

        ta_xxxxz_xxxy_0[i] = ta_xxxx_xxxy_0[i] * pa_z[i] - ta_xxxx_xxxy_1[i] * pc_z[i];

        ta_xxxxz_xxxz_0[i] = ta_xxxx_xxx_0[i] * fe_0 - ta_xxxx_xxx_1[i] * fe_0 + ta_xxxx_xxxz_0[i] * pa_z[i] - ta_xxxx_xxxz_1[i] * pc_z[i];

        ta_xxxxz_xxyy_0[i] = ta_xxxx_xxyy_0[i] * pa_z[i] - ta_xxxx_xxyy_1[i] * pc_z[i];

        ta_xxxxz_xxyz_0[i] = ta_xxxx_xxy_0[i] * fe_0 - ta_xxxx_xxy_1[i] * fe_0 + ta_xxxx_xxyz_0[i] * pa_z[i] - ta_xxxx_xxyz_1[i] * pc_z[i];

        ta_xxxxz_xxzz_0[i] =
            2.0 * ta_xxxx_xxz_0[i] * fe_0 - 2.0 * ta_xxxx_xxz_1[i] * fe_0 + ta_xxxx_xxzz_0[i] * pa_z[i] - ta_xxxx_xxzz_1[i] * pc_z[i];

        ta_xxxxz_xyyy_0[i] = ta_xxxx_xyyy_0[i] * pa_z[i] - ta_xxxx_xyyy_1[i] * pc_z[i];

        ta_xxxxz_xyyz_0[i] = ta_xxxx_xyy_0[i] * fe_0 - ta_xxxx_xyy_1[i] * fe_0 + ta_xxxx_xyyz_0[i] * pa_z[i] - ta_xxxx_xyyz_1[i] * pc_z[i];

        ta_xxxxz_xyzz_0[i] =
            2.0 * ta_xxxx_xyz_0[i] * fe_0 - 2.0 * ta_xxxx_xyz_1[i] * fe_0 + ta_xxxx_xyzz_0[i] * pa_z[i] - ta_xxxx_xyzz_1[i] * pc_z[i];

        ta_xxxxz_xzzz_0[i] =
            3.0 * ta_xxxx_xzz_0[i] * fe_0 - 3.0 * ta_xxxx_xzz_1[i] * fe_0 + ta_xxxx_xzzz_0[i] * pa_z[i] - ta_xxxx_xzzz_1[i] * pc_z[i];

        ta_xxxxz_yyyy_0[i] = ta_xxxx_yyyy_0[i] * pa_z[i] - ta_xxxx_yyyy_1[i] * pc_z[i];

        ta_xxxxz_yyyz_0[i] =
            3.0 * ta_xxz_yyyz_0[i] * fe_0 - 3.0 * ta_xxz_yyyz_1[i] * fe_0 + ta_xxxz_yyyz_0[i] * pa_x[i] - ta_xxxz_yyyz_1[i] * pc_x[i];

        ta_xxxxz_yyzz_0[i] =
            3.0 * ta_xxz_yyzz_0[i] * fe_0 - 3.0 * ta_xxz_yyzz_1[i] * fe_0 + ta_xxxz_yyzz_0[i] * pa_x[i] - ta_xxxz_yyzz_1[i] * pc_x[i];

        ta_xxxxz_yzzz_0[i] =
            3.0 * ta_xxz_yzzz_0[i] * fe_0 - 3.0 * ta_xxz_yzzz_1[i] * fe_0 + ta_xxxz_yzzz_0[i] * pa_x[i] - ta_xxxz_yzzz_1[i] * pc_x[i];

        ta_xxxxz_zzzz_0[i] =
            3.0 * ta_xxz_zzzz_0[i] * fe_0 - 3.0 * ta_xxz_zzzz_1[i] * fe_0 + ta_xxxz_zzzz_0[i] * pa_x[i] - ta_xxxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 45-60 components of targeted buffer : HG

    auto ta_xxxyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 45);

    auto ta_xxxyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 46);

    auto ta_xxxyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 47);

    auto ta_xxxyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 48);

    auto ta_xxxyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 49);

    auto ta_xxxyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 50);

    auto ta_xxxyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 51);

    auto ta_xxxyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 52);

    auto ta_xxxyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 53);

    auto ta_xxxyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 54);

    auto ta_xxxyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 55);

    auto ta_xxxyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 56);

    auto ta_xxxyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 57);

    auto ta_xxxyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 58);

    auto ta_xxxyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 59);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxx_xxxx_0,   \
                             ta_xxx_xxxx_1,   \
                             ta_xxx_xxxz_0,   \
                             ta_xxx_xxxz_1,   \
                             ta_xxx_xxzz_0,   \
                             ta_xxx_xxzz_1,   \
                             ta_xxx_xzzz_0,   \
                             ta_xxx_xzzz_1,   \
                             ta_xxxy_xxxx_0,  \
                             ta_xxxy_xxxx_1,  \
                             ta_xxxy_xxxz_0,  \
                             ta_xxxy_xxxz_1,  \
                             ta_xxxy_xxzz_0,  \
                             ta_xxxy_xxzz_1,  \
                             ta_xxxy_xzzz_0,  \
                             ta_xxxy_xzzz_1,  \
                             ta_xxxyy_xxxx_0, \
                             ta_xxxyy_xxxy_0, \
                             ta_xxxyy_xxxz_0, \
                             ta_xxxyy_xxyy_0, \
                             ta_xxxyy_xxyz_0, \
                             ta_xxxyy_xxzz_0, \
                             ta_xxxyy_xyyy_0, \
                             ta_xxxyy_xyyz_0, \
                             ta_xxxyy_xyzz_0, \
                             ta_xxxyy_xzzz_0, \
                             ta_xxxyy_yyyy_0, \
                             ta_xxxyy_yyyz_0, \
                             ta_xxxyy_yyzz_0, \
                             ta_xxxyy_yzzz_0, \
                             ta_xxxyy_zzzz_0, \
                             ta_xxyy_xxxy_0,  \
                             ta_xxyy_xxxy_1,  \
                             ta_xxyy_xxy_0,   \
                             ta_xxyy_xxy_1,   \
                             ta_xxyy_xxyy_0,  \
                             ta_xxyy_xxyy_1,  \
                             ta_xxyy_xxyz_0,  \
                             ta_xxyy_xxyz_1,  \
                             ta_xxyy_xyy_0,   \
                             ta_xxyy_xyy_1,   \
                             ta_xxyy_xyyy_0,  \
                             ta_xxyy_xyyy_1,  \
                             ta_xxyy_xyyz_0,  \
                             ta_xxyy_xyyz_1,  \
                             ta_xxyy_xyz_0,   \
                             ta_xxyy_xyz_1,   \
                             ta_xxyy_xyzz_0,  \
                             ta_xxyy_xyzz_1,  \
                             ta_xxyy_yyy_0,   \
                             ta_xxyy_yyy_1,   \
                             ta_xxyy_yyyy_0,  \
                             ta_xxyy_yyyy_1,  \
                             ta_xxyy_yyyz_0,  \
                             ta_xxyy_yyyz_1,  \
                             ta_xxyy_yyz_0,   \
                             ta_xxyy_yyz_1,   \
                             ta_xxyy_yyzz_0,  \
                             ta_xxyy_yyzz_1,  \
                             ta_xxyy_yzz_0,   \
                             ta_xxyy_yzz_1,   \
                             ta_xxyy_yzzz_0,  \
                             ta_xxyy_yzzz_1,  \
                             ta_xxyy_zzzz_0,  \
                             ta_xxyy_zzzz_1,  \
                             ta_xyy_xxxy_0,   \
                             ta_xyy_xxxy_1,   \
                             ta_xyy_xxyy_0,   \
                             ta_xyy_xxyy_1,   \
                             ta_xyy_xxyz_0,   \
                             ta_xyy_xxyz_1,   \
                             ta_xyy_xyyy_0,   \
                             ta_xyy_xyyy_1,   \
                             ta_xyy_xyyz_0,   \
                             ta_xyy_xyyz_1,   \
                             ta_xyy_xyzz_0,   \
                             ta_xyy_xyzz_1,   \
                             ta_xyy_yyyy_0,   \
                             ta_xyy_yyyy_1,   \
                             ta_xyy_yyyz_0,   \
                             ta_xyy_yyyz_1,   \
                             ta_xyy_yyzz_0,   \
                             ta_xyy_yyzz_1,   \
                             ta_xyy_yzzz_0,   \
                             ta_xyy_yzzz_1,   \
                             ta_xyy_zzzz_0,   \
                             ta_xyy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyy_xxxx_0[i] = ta_xxx_xxxx_0[i] * fe_0 - ta_xxx_xxxx_1[i] * fe_0 + ta_xxxy_xxxx_0[i] * pa_y[i] - ta_xxxy_xxxx_1[i] * pc_y[i];

        ta_xxxyy_xxxy_0[i] = 2.0 * ta_xyy_xxxy_0[i] * fe_0 - 2.0 * ta_xyy_xxxy_1[i] * fe_0 + 3.0 * ta_xxyy_xxy_0[i] * fe_0 -
                             3.0 * ta_xxyy_xxy_1[i] * fe_0 + ta_xxyy_xxxy_0[i] * pa_x[i] - ta_xxyy_xxxy_1[i] * pc_x[i];

        ta_xxxyy_xxxz_0[i] = ta_xxx_xxxz_0[i] * fe_0 - ta_xxx_xxxz_1[i] * fe_0 + ta_xxxy_xxxz_0[i] * pa_y[i] - ta_xxxy_xxxz_1[i] * pc_y[i];

        ta_xxxyy_xxyy_0[i] = 2.0 * ta_xyy_xxyy_0[i] * fe_0 - 2.0 * ta_xyy_xxyy_1[i] * fe_0 + 2.0 * ta_xxyy_xyy_0[i] * fe_0 -
                             2.0 * ta_xxyy_xyy_1[i] * fe_0 + ta_xxyy_xxyy_0[i] * pa_x[i] - ta_xxyy_xxyy_1[i] * pc_x[i];

        ta_xxxyy_xxyz_0[i] = 2.0 * ta_xyy_xxyz_0[i] * fe_0 - 2.0 * ta_xyy_xxyz_1[i] * fe_0 + 2.0 * ta_xxyy_xyz_0[i] * fe_0 -
                             2.0 * ta_xxyy_xyz_1[i] * fe_0 + ta_xxyy_xxyz_0[i] * pa_x[i] - ta_xxyy_xxyz_1[i] * pc_x[i];

        ta_xxxyy_xxzz_0[i] = ta_xxx_xxzz_0[i] * fe_0 - ta_xxx_xxzz_1[i] * fe_0 + ta_xxxy_xxzz_0[i] * pa_y[i] - ta_xxxy_xxzz_1[i] * pc_y[i];

        ta_xxxyy_xyyy_0[i] = 2.0 * ta_xyy_xyyy_0[i] * fe_0 - 2.0 * ta_xyy_xyyy_1[i] * fe_0 + ta_xxyy_yyy_0[i] * fe_0 - ta_xxyy_yyy_1[i] * fe_0 +
                             ta_xxyy_xyyy_0[i] * pa_x[i] - ta_xxyy_xyyy_1[i] * pc_x[i];

        ta_xxxyy_xyyz_0[i] = 2.0 * ta_xyy_xyyz_0[i] * fe_0 - 2.0 * ta_xyy_xyyz_1[i] * fe_0 + ta_xxyy_yyz_0[i] * fe_0 - ta_xxyy_yyz_1[i] * fe_0 +
                             ta_xxyy_xyyz_0[i] * pa_x[i] - ta_xxyy_xyyz_1[i] * pc_x[i];

        ta_xxxyy_xyzz_0[i] = 2.0 * ta_xyy_xyzz_0[i] * fe_0 - 2.0 * ta_xyy_xyzz_1[i] * fe_0 + ta_xxyy_yzz_0[i] * fe_0 - ta_xxyy_yzz_1[i] * fe_0 +
                             ta_xxyy_xyzz_0[i] * pa_x[i] - ta_xxyy_xyzz_1[i] * pc_x[i];

        ta_xxxyy_xzzz_0[i] = ta_xxx_xzzz_0[i] * fe_0 - ta_xxx_xzzz_1[i] * fe_0 + ta_xxxy_xzzz_0[i] * pa_y[i] - ta_xxxy_xzzz_1[i] * pc_y[i];

        ta_xxxyy_yyyy_0[i] =
            2.0 * ta_xyy_yyyy_0[i] * fe_0 - 2.0 * ta_xyy_yyyy_1[i] * fe_0 + ta_xxyy_yyyy_0[i] * pa_x[i] - ta_xxyy_yyyy_1[i] * pc_x[i];

        ta_xxxyy_yyyz_0[i] =
            2.0 * ta_xyy_yyyz_0[i] * fe_0 - 2.0 * ta_xyy_yyyz_1[i] * fe_0 + ta_xxyy_yyyz_0[i] * pa_x[i] - ta_xxyy_yyyz_1[i] * pc_x[i];

        ta_xxxyy_yyzz_0[i] =
            2.0 * ta_xyy_yyzz_0[i] * fe_0 - 2.0 * ta_xyy_yyzz_1[i] * fe_0 + ta_xxyy_yyzz_0[i] * pa_x[i] - ta_xxyy_yyzz_1[i] * pc_x[i];

        ta_xxxyy_yzzz_0[i] =
            2.0 * ta_xyy_yzzz_0[i] * fe_0 - 2.0 * ta_xyy_yzzz_1[i] * fe_0 + ta_xxyy_yzzz_0[i] * pa_x[i] - ta_xxyy_yzzz_1[i] * pc_x[i];

        ta_xxxyy_zzzz_0[i] =
            2.0 * ta_xyy_zzzz_0[i] * fe_0 - 2.0 * ta_xyy_zzzz_1[i] * fe_0 + ta_xxyy_zzzz_0[i] * pa_x[i] - ta_xxyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : HG

    auto ta_xxxyz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 60);

    auto ta_xxxyz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 61);

    auto ta_xxxyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 62);

    auto ta_xxxyz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 63);

    auto ta_xxxyz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 64);

    auto ta_xxxyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 65);

    auto ta_xxxyz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 66);

    auto ta_xxxyz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 67);

    auto ta_xxxyz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 68);

    auto ta_xxxyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 69);

    auto ta_xxxyz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 70);

    auto ta_xxxyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 71);

    auto ta_xxxyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 72);

    auto ta_xxxyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 73);

    auto ta_xxxyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 74);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta_xxxy_xxxy_0,  \
                             ta_xxxy_xxxy_1,  \
                             ta_xxxy_xxyy_0,  \
                             ta_xxxy_xxyy_1,  \
                             ta_xxxy_xyyy_0,  \
                             ta_xxxy_xyyy_1,  \
                             ta_xxxy_yyyy_0,  \
                             ta_xxxy_yyyy_1,  \
                             ta_xxxyz_xxxx_0, \
                             ta_xxxyz_xxxy_0, \
                             ta_xxxyz_xxxz_0, \
                             ta_xxxyz_xxyy_0, \
                             ta_xxxyz_xxyz_0, \
                             ta_xxxyz_xxzz_0, \
                             ta_xxxyz_xyyy_0, \
                             ta_xxxyz_xyyz_0, \
                             ta_xxxyz_xyzz_0, \
                             ta_xxxyz_xzzz_0, \
                             ta_xxxyz_yyyy_0, \
                             ta_xxxyz_yyyz_0, \
                             ta_xxxyz_yyzz_0, \
                             ta_xxxyz_yzzz_0, \
                             ta_xxxyz_zzzz_0, \
                             ta_xxxz_xxxx_0,  \
                             ta_xxxz_xxxx_1,  \
                             ta_xxxz_xxxz_0,  \
                             ta_xxxz_xxxz_1,  \
                             ta_xxxz_xxyz_0,  \
                             ta_xxxz_xxyz_1,  \
                             ta_xxxz_xxz_0,   \
                             ta_xxxz_xxz_1,   \
                             ta_xxxz_xxzz_0,  \
                             ta_xxxz_xxzz_1,  \
                             ta_xxxz_xyyz_0,  \
                             ta_xxxz_xyyz_1,  \
                             ta_xxxz_xyz_0,   \
                             ta_xxxz_xyz_1,   \
                             ta_xxxz_xyzz_0,  \
                             ta_xxxz_xyzz_1,  \
                             ta_xxxz_xzz_0,   \
                             ta_xxxz_xzz_1,   \
                             ta_xxxz_xzzz_0,  \
                             ta_xxxz_xzzz_1,  \
                             ta_xxxz_zzzz_0,  \
                             ta_xxxz_zzzz_1,  \
                             ta_xxyz_yyyz_0,  \
                             ta_xxyz_yyyz_1,  \
                             ta_xxyz_yyzz_0,  \
                             ta_xxyz_yyzz_1,  \
                             ta_xxyz_yzzz_0,  \
                             ta_xxyz_yzzz_1,  \
                             ta_xyz_yyyz_0,   \
                             ta_xyz_yyyz_1,   \
                             ta_xyz_yyzz_0,   \
                             ta_xyz_yyzz_1,   \
                             ta_xyz_yzzz_0,   \
                             ta_xyz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyz_xxxx_0[i] = ta_xxxz_xxxx_0[i] * pa_y[i] - ta_xxxz_xxxx_1[i] * pc_y[i];

        ta_xxxyz_xxxy_0[i] = ta_xxxy_xxxy_0[i] * pa_z[i] - ta_xxxy_xxxy_1[i] * pc_z[i];

        ta_xxxyz_xxxz_0[i] = ta_xxxz_xxxz_0[i] * pa_y[i] - ta_xxxz_xxxz_1[i] * pc_y[i];

        ta_xxxyz_xxyy_0[i] = ta_xxxy_xxyy_0[i] * pa_z[i] - ta_xxxy_xxyy_1[i] * pc_z[i];

        ta_xxxyz_xxyz_0[i] = ta_xxxz_xxz_0[i] * fe_0 - ta_xxxz_xxz_1[i] * fe_0 + ta_xxxz_xxyz_0[i] * pa_y[i] - ta_xxxz_xxyz_1[i] * pc_y[i];

        ta_xxxyz_xxzz_0[i] = ta_xxxz_xxzz_0[i] * pa_y[i] - ta_xxxz_xxzz_1[i] * pc_y[i];

        ta_xxxyz_xyyy_0[i] = ta_xxxy_xyyy_0[i] * pa_z[i] - ta_xxxy_xyyy_1[i] * pc_z[i];

        ta_xxxyz_xyyz_0[i] =
            2.0 * ta_xxxz_xyz_0[i] * fe_0 - 2.0 * ta_xxxz_xyz_1[i] * fe_0 + ta_xxxz_xyyz_0[i] * pa_y[i] - ta_xxxz_xyyz_1[i] * pc_y[i];

        ta_xxxyz_xyzz_0[i] = ta_xxxz_xzz_0[i] * fe_0 - ta_xxxz_xzz_1[i] * fe_0 + ta_xxxz_xyzz_0[i] * pa_y[i] - ta_xxxz_xyzz_1[i] * pc_y[i];

        ta_xxxyz_xzzz_0[i] = ta_xxxz_xzzz_0[i] * pa_y[i] - ta_xxxz_xzzz_1[i] * pc_y[i];

        ta_xxxyz_yyyy_0[i] = ta_xxxy_yyyy_0[i] * pa_z[i] - ta_xxxy_yyyy_1[i] * pc_z[i];

        ta_xxxyz_yyyz_0[i] =
            2.0 * ta_xyz_yyyz_0[i] * fe_0 - 2.0 * ta_xyz_yyyz_1[i] * fe_0 + ta_xxyz_yyyz_0[i] * pa_x[i] - ta_xxyz_yyyz_1[i] * pc_x[i];

        ta_xxxyz_yyzz_0[i] =
            2.0 * ta_xyz_yyzz_0[i] * fe_0 - 2.0 * ta_xyz_yyzz_1[i] * fe_0 + ta_xxyz_yyzz_0[i] * pa_x[i] - ta_xxyz_yyzz_1[i] * pc_x[i];

        ta_xxxyz_yzzz_0[i] =
            2.0 * ta_xyz_yzzz_0[i] * fe_0 - 2.0 * ta_xyz_yzzz_1[i] * fe_0 + ta_xxyz_yzzz_0[i] * pa_x[i] - ta_xxyz_yzzz_1[i] * pc_x[i];

        ta_xxxyz_zzzz_0[i] = ta_xxxz_zzzz_0[i] * pa_y[i] - ta_xxxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : HG

    auto ta_xxxzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 75);

    auto ta_xxxzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 76);

    auto ta_xxxzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 77);

    auto ta_xxxzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 78);

    auto ta_xxxzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 79);

    auto ta_xxxzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 80);

    auto ta_xxxzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 81);

    auto ta_xxxzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 82);

    auto ta_xxxzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 83);

    auto ta_xxxzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 84);

    auto ta_xxxzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 85);

    auto ta_xxxzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 86);

    auto ta_xxxzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 87);

    auto ta_xxxzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 88);

    auto ta_xxxzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 89);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xxx_xxxx_0,   \
                             ta_xxx_xxxx_1,   \
                             ta_xxx_xxxy_0,   \
                             ta_xxx_xxxy_1,   \
                             ta_xxx_xxyy_0,   \
                             ta_xxx_xxyy_1,   \
                             ta_xxx_xyyy_0,   \
                             ta_xxx_xyyy_1,   \
                             ta_xxxz_xxxx_0,  \
                             ta_xxxz_xxxx_1,  \
                             ta_xxxz_xxxy_0,  \
                             ta_xxxz_xxxy_1,  \
                             ta_xxxz_xxyy_0,  \
                             ta_xxxz_xxyy_1,  \
                             ta_xxxz_xyyy_0,  \
                             ta_xxxz_xyyy_1,  \
                             ta_xxxzz_xxxx_0, \
                             ta_xxxzz_xxxy_0, \
                             ta_xxxzz_xxxz_0, \
                             ta_xxxzz_xxyy_0, \
                             ta_xxxzz_xxyz_0, \
                             ta_xxxzz_xxzz_0, \
                             ta_xxxzz_xyyy_0, \
                             ta_xxxzz_xyyz_0, \
                             ta_xxxzz_xyzz_0, \
                             ta_xxxzz_xzzz_0, \
                             ta_xxxzz_yyyy_0, \
                             ta_xxxzz_yyyz_0, \
                             ta_xxxzz_yyzz_0, \
                             ta_xxxzz_yzzz_0, \
                             ta_xxxzz_zzzz_0, \
                             ta_xxzz_xxxz_0,  \
                             ta_xxzz_xxxz_1,  \
                             ta_xxzz_xxyz_0,  \
                             ta_xxzz_xxyz_1,  \
                             ta_xxzz_xxz_0,   \
                             ta_xxzz_xxz_1,   \
                             ta_xxzz_xxzz_0,  \
                             ta_xxzz_xxzz_1,  \
                             ta_xxzz_xyyz_0,  \
                             ta_xxzz_xyyz_1,  \
                             ta_xxzz_xyz_0,   \
                             ta_xxzz_xyz_1,   \
                             ta_xxzz_xyzz_0,  \
                             ta_xxzz_xyzz_1,  \
                             ta_xxzz_xzz_0,   \
                             ta_xxzz_xzz_1,   \
                             ta_xxzz_xzzz_0,  \
                             ta_xxzz_xzzz_1,  \
                             ta_xxzz_yyyy_0,  \
                             ta_xxzz_yyyy_1,  \
                             ta_xxzz_yyyz_0,  \
                             ta_xxzz_yyyz_1,  \
                             ta_xxzz_yyz_0,   \
                             ta_xxzz_yyz_1,   \
                             ta_xxzz_yyzz_0,  \
                             ta_xxzz_yyzz_1,  \
                             ta_xxzz_yzz_0,   \
                             ta_xxzz_yzz_1,   \
                             ta_xxzz_yzzz_0,  \
                             ta_xxzz_yzzz_1,  \
                             ta_xxzz_zzz_0,   \
                             ta_xxzz_zzz_1,   \
                             ta_xxzz_zzzz_0,  \
                             ta_xxzz_zzzz_1,  \
                             ta_xzz_xxxz_0,   \
                             ta_xzz_xxxz_1,   \
                             ta_xzz_xxyz_0,   \
                             ta_xzz_xxyz_1,   \
                             ta_xzz_xxzz_0,   \
                             ta_xzz_xxzz_1,   \
                             ta_xzz_xyyz_0,   \
                             ta_xzz_xyyz_1,   \
                             ta_xzz_xyzz_0,   \
                             ta_xzz_xyzz_1,   \
                             ta_xzz_xzzz_0,   \
                             ta_xzz_xzzz_1,   \
                             ta_xzz_yyyy_0,   \
                             ta_xzz_yyyy_1,   \
                             ta_xzz_yyyz_0,   \
                             ta_xzz_yyyz_1,   \
                             ta_xzz_yyzz_0,   \
                             ta_xzz_yyzz_1,   \
                             ta_xzz_yzzz_0,   \
                             ta_xzz_yzzz_1,   \
                             ta_xzz_zzzz_0,   \
                             ta_xzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzz_xxxx_0[i] = ta_xxx_xxxx_0[i] * fe_0 - ta_xxx_xxxx_1[i] * fe_0 + ta_xxxz_xxxx_0[i] * pa_z[i] - ta_xxxz_xxxx_1[i] * pc_z[i];

        ta_xxxzz_xxxy_0[i] = ta_xxx_xxxy_0[i] * fe_0 - ta_xxx_xxxy_1[i] * fe_0 + ta_xxxz_xxxy_0[i] * pa_z[i] - ta_xxxz_xxxy_1[i] * pc_z[i];

        ta_xxxzz_xxxz_0[i] = 2.0 * ta_xzz_xxxz_0[i] * fe_0 - 2.0 * ta_xzz_xxxz_1[i] * fe_0 + 3.0 * ta_xxzz_xxz_0[i] * fe_0 -
                             3.0 * ta_xxzz_xxz_1[i] * fe_0 + ta_xxzz_xxxz_0[i] * pa_x[i] - ta_xxzz_xxxz_1[i] * pc_x[i];

        ta_xxxzz_xxyy_0[i] = ta_xxx_xxyy_0[i] * fe_0 - ta_xxx_xxyy_1[i] * fe_0 + ta_xxxz_xxyy_0[i] * pa_z[i] - ta_xxxz_xxyy_1[i] * pc_z[i];

        ta_xxxzz_xxyz_0[i] = 2.0 * ta_xzz_xxyz_0[i] * fe_0 - 2.0 * ta_xzz_xxyz_1[i] * fe_0 + 2.0 * ta_xxzz_xyz_0[i] * fe_0 -
                             2.0 * ta_xxzz_xyz_1[i] * fe_0 + ta_xxzz_xxyz_0[i] * pa_x[i] - ta_xxzz_xxyz_1[i] * pc_x[i];

        ta_xxxzz_xxzz_0[i] = 2.0 * ta_xzz_xxzz_0[i] * fe_0 - 2.0 * ta_xzz_xxzz_1[i] * fe_0 + 2.0 * ta_xxzz_xzz_0[i] * fe_0 -
                             2.0 * ta_xxzz_xzz_1[i] * fe_0 + ta_xxzz_xxzz_0[i] * pa_x[i] - ta_xxzz_xxzz_1[i] * pc_x[i];

        ta_xxxzz_xyyy_0[i] = ta_xxx_xyyy_0[i] * fe_0 - ta_xxx_xyyy_1[i] * fe_0 + ta_xxxz_xyyy_0[i] * pa_z[i] - ta_xxxz_xyyy_1[i] * pc_z[i];

        ta_xxxzz_xyyz_0[i] = 2.0 * ta_xzz_xyyz_0[i] * fe_0 - 2.0 * ta_xzz_xyyz_1[i] * fe_0 + ta_xxzz_yyz_0[i] * fe_0 - ta_xxzz_yyz_1[i] * fe_0 +
                             ta_xxzz_xyyz_0[i] * pa_x[i] - ta_xxzz_xyyz_1[i] * pc_x[i];

        ta_xxxzz_xyzz_0[i] = 2.0 * ta_xzz_xyzz_0[i] * fe_0 - 2.0 * ta_xzz_xyzz_1[i] * fe_0 + ta_xxzz_yzz_0[i] * fe_0 - ta_xxzz_yzz_1[i] * fe_0 +
                             ta_xxzz_xyzz_0[i] * pa_x[i] - ta_xxzz_xyzz_1[i] * pc_x[i];

        ta_xxxzz_xzzz_0[i] = 2.0 * ta_xzz_xzzz_0[i] * fe_0 - 2.0 * ta_xzz_xzzz_1[i] * fe_0 + ta_xxzz_zzz_0[i] * fe_0 - ta_xxzz_zzz_1[i] * fe_0 +
                             ta_xxzz_xzzz_0[i] * pa_x[i] - ta_xxzz_xzzz_1[i] * pc_x[i];

        ta_xxxzz_yyyy_0[i] =
            2.0 * ta_xzz_yyyy_0[i] * fe_0 - 2.0 * ta_xzz_yyyy_1[i] * fe_0 + ta_xxzz_yyyy_0[i] * pa_x[i] - ta_xxzz_yyyy_1[i] * pc_x[i];

        ta_xxxzz_yyyz_0[i] =
            2.0 * ta_xzz_yyyz_0[i] * fe_0 - 2.0 * ta_xzz_yyyz_1[i] * fe_0 + ta_xxzz_yyyz_0[i] * pa_x[i] - ta_xxzz_yyyz_1[i] * pc_x[i];

        ta_xxxzz_yyzz_0[i] =
            2.0 * ta_xzz_yyzz_0[i] * fe_0 - 2.0 * ta_xzz_yyzz_1[i] * fe_0 + ta_xxzz_yyzz_0[i] * pa_x[i] - ta_xxzz_yyzz_1[i] * pc_x[i];

        ta_xxxzz_yzzz_0[i] =
            2.0 * ta_xzz_yzzz_0[i] * fe_0 - 2.0 * ta_xzz_yzzz_1[i] * fe_0 + ta_xxzz_yzzz_0[i] * pa_x[i] - ta_xxzz_yzzz_1[i] * pc_x[i];

        ta_xxxzz_zzzz_0[i] =
            2.0 * ta_xzz_zzzz_0[i] * fe_0 - 2.0 * ta_xzz_zzzz_1[i] * fe_0 + ta_xxzz_zzzz_0[i] * pa_x[i] - ta_xxzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : HG

    auto ta_xxyyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 90);

    auto ta_xxyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 91);

    auto ta_xxyyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 92);

    auto ta_xxyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 93);

    auto ta_xxyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 94);

    auto ta_xxyyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 95);

    auto ta_xxyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 96);

    auto ta_xxyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 97);

    auto ta_xxyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 98);

    auto ta_xxyyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 99);

    auto ta_xxyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 100);

    auto ta_xxyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 101);

    auto ta_xxyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 102);

    auto ta_xxyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 103);

    auto ta_xxyyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 104);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxy_xxxx_0,   \
                             ta_xxy_xxxx_1,   \
                             ta_xxy_xxxz_0,   \
                             ta_xxy_xxxz_1,   \
                             ta_xxy_xxzz_0,   \
                             ta_xxy_xxzz_1,   \
                             ta_xxy_xzzz_0,   \
                             ta_xxy_xzzz_1,   \
                             ta_xxyy_xxxx_0,  \
                             ta_xxyy_xxxx_1,  \
                             ta_xxyy_xxxz_0,  \
                             ta_xxyy_xxxz_1,  \
                             ta_xxyy_xxzz_0,  \
                             ta_xxyy_xxzz_1,  \
                             ta_xxyy_xzzz_0,  \
                             ta_xxyy_xzzz_1,  \
                             ta_xxyyy_xxxx_0, \
                             ta_xxyyy_xxxy_0, \
                             ta_xxyyy_xxxz_0, \
                             ta_xxyyy_xxyy_0, \
                             ta_xxyyy_xxyz_0, \
                             ta_xxyyy_xxzz_0, \
                             ta_xxyyy_xyyy_0, \
                             ta_xxyyy_xyyz_0, \
                             ta_xxyyy_xyzz_0, \
                             ta_xxyyy_xzzz_0, \
                             ta_xxyyy_yyyy_0, \
                             ta_xxyyy_yyyz_0, \
                             ta_xxyyy_yyzz_0, \
                             ta_xxyyy_yzzz_0, \
                             ta_xxyyy_zzzz_0, \
                             ta_xyyy_xxxy_0,  \
                             ta_xyyy_xxxy_1,  \
                             ta_xyyy_xxy_0,   \
                             ta_xyyy_xxy_1,   \
                             ta_xyyy_xxyy_0,  \
                             ta_xyyy_xxyy_1,  \
                             ta_xyyy_xxyz_0,  \
                             ta_xyyy_xxyz_1,  \
                             ta_xyyy_xyy_0,   \
                             ta_xyyy_xyy_1,   \
                             ta_xyyy_xyyy_0,  \
                             ta_xyyy_xyyy_1,  \
                             ta_xyyy_xyyz_0,  \
                             ta_xyyy_xyyz_1,  \
                             ta_xyyy_xyz_0,   \
                             ta_xyyy_xyz_1,   \
                             ta_xyyy_xyzz_0,  \
                             ta_xyyy_xyzz_1,  \
                             ta_xyyy_yyy_0,   \
                             ta_xyyy_yyy_1,   \
                             ta_xyyy_yyyy_0,  \
                             ta_xyyy_yyyy_1,  \
                             ta_xyyy_yyyz_0,  \
                             ta_xyyy_yyyz_1,  \
                             ta_xyyy_yyz_0,   \
                             ta_xyyy_yyz_1,   \
                             ta_xyyy_yyzz_0,  \
                             ta_xyyy_yyzz_1,  \
                             ta_xyyy_yzz_0,   \
                             ta_xyyy_yzz_1,   \
                             ta_xyyy_yzzz_0,  \
                             ta_xyyy_yzzz_1,  \
                             ta_xyyy_zzzz_0,  \
                             ta_xyyy_zzzz_1,  \
                             ta_yyy_xxxy_0,   \
                             ta_yyy_xxxy_1,   \
                             ta_yyy_xxyy_0,   \
                             ta_yyy_xxyy_1,   \
                             ta_yyy_xxyz_0,   \
                             ta_yyy_xxyz_1,   \
                             ta_yyy_xyyy_0,   \
                             ta_yyy_xyyy_1,   \
                             ta_yyy_xyyz_0,   \
                             ta_yyy_xyyz_1,   \
                             ta_yyy_xyzz_0,   \
                             ta_yyy_xyzz_1,   \
                             ta_yyy_yyyy_0,   \
                             ta_yyy_yyyy_1,   \
                             ta_yyy_yyyz_0,   \
                             ta_yyy_yyyz_1,   \
                             ta_yyy_yyzz_0,   \
                             ta_yyy_yyzz_1,   \
                             ta_yyy_yzzz_0,   \
                             ta_yyy_yzzz_1,   \
                             ta_yyy_zzzz_0,   \
                             ta_yyy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyy_xxxx_0[i] =
            2.0 * ta_xxy_xxxx_0[i] * fe_0 - 2.0 * ta_xxy_xxxx_1[i] * fe_0 + ta_xxyy_xxxx_0[i] * pa_y[i] - ta_xxyy_xxxx_1[i] * pc_y[i];

        ta_xxyyy_xxxy_0[i] = ta_yyy_xxxy_0[i] * fe_0 - ta_yyy_xxxy_1[i] * fe_0 + 3.0 * ta_xyyy_xxy_0[i] * fe_0 - 3.0 * ta_xyyy_xxy_1[i] * fe_0 +
                             ta_xyyy_xxxy_0[i] * pa_x[i] - ta_xyyy_xxxy_1[i] * pc_x[i];

        ta_xxyyy_xxxz_0[i] =
            2.0 * ta_xxy_xxxz_0[i] * fe_0 - 2.0 * ta_xxy_xxxz_1[i] * fe_0 + ta_xxyy_xxxz_0[i] * pa_y[i] - ta_xxyy_xxxz_1[i] * pc_y[i];

        ta_xxyyy_xxyy_0[i] = ta_yyy_xxyy_0[i] * fe_0 - ta_yyy_xxyy_1[i] * fe_0 + 2.0 * ta_xyyy_xyy_0[i] * fe_0 - 2.0 * ta_xyyy_xyy_1[i] * fe_0 +
                             ta_xyyy_xxyy_0[i] * pa_x[i] - ta_xyyy_xxyy_1[i] * pc_x[i];

        ta_xxyyy_xxyz_0[i] = ta_yyy_xxyz_0[i] * fe_0 - ta_yyy_xxyz_1[i] * fe_0 + 2.0 * ta_xyyy_xyz_0[i] * fe_0 - 2.0 * ta_xyyy_xyz_1[i] * fe_0 +
                             ta_xyyy_xxyz_0[i] * pa_x[i] - ta_xyyy_xxyz_1[i] * pc_x[i];

        ta_xxyyy_xxzz_0[i] =
            2.0 * ta_xxy_xxzz_0[i] * fe_0 - 2.0 * ta_xxy_xxzz_1[i] * fe_0 + ta_xxyy_xxzz_0[i] * pa_y[i] - ta_xxyy_xxzz_1[i] * pc_y[i];

        ta_xxyyy_xyyy_0[i] = ta_yyy_xyyy_0[i] * fe_0 - ta_yyy_xyyy_1[i] * fe_0 + ta_xyyy_yyy_0[i] * fe_0 - ta_xyyy_yyy_1[i] * fe_0 +
                             ta_xyyy_xyyy_0[i] * pa_x[i] - ta_xyyy_xyyy_1[i] * pc_x[i];

        ta_xxyyy_xyyz_0[i] = ta_yyy_xyyz_0[i] * fe_0 - ta_yyy_xyyz_1[i] * fe_0 + ta_xyyy_yyz_0[i] * fe_0 - ta_xyyy_yyz_1[i] * fe_0 +
                             ta_xyyy_xyyz_0[i] * pa_x[i] - ta_xyyy_xyyz_1[i] * pc_x[i];

        ta_xxyyy_xyzz_0[i] = ta_yyy_xyzz_0[i] * fe_0 - ta_yyy_xyzz_1[i] * fe_0 + ta_xyyy_yzz_0[i] * fe_0 - ta_xyyy_yzz_1[i] * fe_0 +
                             ta_xyyy_xyzz_0[i] * pa_x[i] - ta_xyyy_xyzz_1[i] * pc_x[i];

        ta_xxyyy_xzzz_0[i] =
            2.0 * ta_xxy_xzzz_0[i] * fe_0 - 2.0 * ta_xxy_xzzz_1[i] * fe_0 + ta_xxyy_xzzz_0[i] * pa_y[i] - ta_xxyy_xzzz_1[i] * pc_y[i];

        ta_xxyyy_yyyy_0[i] = ta_yyy_yyyy_0[i] * fe_0 - ta_yyy_yyyy_1[i] * fe_0 + ta_xyyy_yyyy_0[i] * pa_x[i] - ta_xyyy_yyyy_1[i] * pc_x[i];

        ta_xxyyy_yyyz_0[i] = ta_yyy_yyyz_0[i] * fe_0 - ta_yyy_yyyz_1[i] * fe_0 + ta_xyyy_yyyz_0[i] * pa_x[i] - ta_xyyy_yyyz_1[i] * pc_x[i];

        ta_xxyyy_yyzz_0[i] = ta_yyy_yyzz_0[i] * fe_0 - ta_yyy_yyzz_1[i] * fe_0 + ta_xyyy_yyzz_0[i] * pa_x[i] - ta_xyyy_yyzz_1[i] * pc_x[i];

        ta_xxyyy_yzzz_0[i] = ta_yyy_yzzz_0[i] * fe_0 - ta_yyy_yzzz_1[i] * fe_0 + ta_xyyy_yzzz_0[i] * pa_x[i] - ta_xyyy_yzzz_1[i] * pc_x[i];

        ta_xxyyy_zzzz_0[i] = ta_yyy_zzzz_0[i] * fe_0 - ta_yyy_zzzz_1[i] * fe_0 + ta_xyyy_zzzz_0[i] * pa_x[i] - ta_xyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : HG

    auto ta_xxyyz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 105);

    auto ta_xxyyz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 106);

    auto ta_xxyyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 107);

    auto ta_xxyyz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 108);

    auto ta_xxyyz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 109);

    auto ta_xxyyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 110);

    auto ta_xxyyz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 111);

    auto ta_xxyyz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 112);

    auto ta_xxyyz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 113);

    auto ta_xxyyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 114);

    auto ta_xxyyz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 115);

    auto ta_xxyyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 116);

    auto ta_xxyyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 117);

    auto ta_xxyyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 118);

    auto ta_xxyyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 119);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta_xxyy_xxxx_0,  \
                             ta_xxyy_xxxx_1,  \
                             ta_xxyy_xxxy_0,  \
                             ta_xxyy_xxxy_1,  \
                             ta_xxyy_xxy_0,   \
                             ta_xxyy_xxy_1,   \
                             ta_xxyy_xxyy_0,  \
                             ta_xxyy_xxyy_1,  \
                             ta_xxyy_xxyz_0,  \
                             ta_xxyy_xxyz_1,  \
                             ta_xxyy_xyy_0,   \
                             ta_xxyy_xyy_1,   \
                             ta_xxyy_xyyy_0,  \
                             ta_xxyy_xyyy_1,  \
                             ta_xxyy_xyyz_0,  \
                             ta_xxyy_xyyz_1,  \
                             ta_xxyy_xyz_0,   \
                             ta_xxyy_xyz_1,   \
                             ta_xxyy_xyzz_0,  \
                             ta_xxyy_xyzz_1,  \
                             ta_xxyy_yyyy_0,  \
                             ta_xxyy_yyyy_1,  \
                             ta_xxyyz_xxxx_0, \
                             ta_xxyyz_xxxy_0, \
                             ta_xxyyz_xxxz_0, \
                             ta_xxyyz_xxyy_0, \
                             ta_xxyyz_xxyz_0, \
                             ta_xxyyz_xxzz_0, \
                             ta_xxyyz_xyyy_0, \
                             ta_xxyyz_xyyz_0, \
                             ta_xxyyz_xyzz_0, \
                             ta_xxyyz_xzzz_0, \
                             ta_xxyyz_yyyy_0, \
                             ta_xxyyz_yyyz_0, \
                             ta_xxyyz_yyzz_0, \
                             ta_xxyyz_yzzz_0, \
                             ta_xxyyz_zzzz_0, \
                             ta_xxyz_xxxz_0,  \
                             ta_xxyz_xxxz_1,  \
                             ta_xxyz_xxzz_0,  \
                             ta_xxyz_xxzz_1,  \
                             ta_xxyz_xzzz_0,  \
                             ta_xxyz_xzzz_1,  \
                             ta_xxz_xxxz_0,   \
                             ta_xxz_xxxz_1,   \
                             ta_xxz_xxzz_0,   \
                             ta_xxz_xxzz_1,   \
                             ta_xxz_xzzz_0,   \
                             ta_xxz_xzzz_1,   \
                             ta_xyyz_yyyz_0,  \
                             ta_xyyz_yyyz_1,  \
                             ta_xyyz_yyzz_0,  \
                             ta_xyyz_yyzz_1,  \
                             ta_xyyz_yzzz_0,  \
                             ta_xyyz_yzzz_1,  \
                             ta_xyyz_zzzz_0,  \
                             ta_xyyz_zzzz_1,  \
                             ta_yyz_yyyz_0,   \
                             ta_yyz_yyyz_1,   \
                             ta_yyz_yyzz_0,   \
                             ta_yyz_yyzz_1,   \
                             ta_yyz_yzzz_0,   \
                             ta_yyz_yzzz_1,   \
                             ta_yyz_zzzz_0,   \
                             ta_yyz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyz_xxxx_0[i] = ta_xxyy_xxxx_0[i] * pa_z[i] - ta_xxyy_xxxx_1[i] * pc_z[i];

        ta_xxyyz_xxxy_0[i] = ta_xxyy_xxxy_0[i] * pa_z[i] - ta_xxyy_xxxy_1[i] * pc_z[i];

        ta_xxyyz_xxxz_0[i] = ta_xxz_xxxz_0[i] * fe_0 - ta_xxz_xxxz_1[i] * fe_0 + ta_xxyz_xxxz_0[i] * pa_y[i] - ta_xxyz_xxxz_1[i] * pc_y[i];

        ta_xxyyz_xxyy_0[i] = ta_xxyy_xxyy_0[i] * pa_z[i] - ta_xxyy_xxyy_1[i] * pc_z[i];

        ta_xxyyz_xxyz_0[i] = ta_xxyy_xxy_0[i] * fe_0 - ta_xxyy_xxy_1[i] * fe_0 + ta_xxyy_xxyz_0[i] * pa_z[i] - ta_xxyy_xxyz_1[i] * pc_z[i];

        ta_xxyyz_xxzz_0[i] = ta_xxz_xxzz_0[i] * fe_0 - ta_xxz_xxzz_1[i] * fe_0 + ta_xxyz_xxzz_0[i] * pa_y[i] - ta_xxyz_xxzz_1[i] * pc_y[i];

        ta_xxyyz_xyyy_0[i] = ta_xxyy_xyyy_0[i] * pa_z[i] - ta_xxyy_xyyy_1[i] * pc_z[i];

        ta_xxyyz_xyyz_0[i] = ta_xxyy_xyy_0[i] * fe_0 - ta_xxyy_xyy_1[i] * fe_0 + ta_xxyy_xyyz_0[i] * pa_z[i] - ta_xxyy_xyyz_1[i] * pc_z[i];

        ta_xxyyz_xyzz_0[i] =
            2.0 * ta_xxyy_xyz_0[i] * fe_0 - 2.0 * ta_xxyy_xyz_1[i] * fe_0 + ta_xxyy_xyzz_0[i] * pa_z[i] - ta_xxyy_xyzz_1[i] * pc_z[i];

        ta_xxyyz_xzzz_0[i] = ta_xxz_xzzz_0[i] * fe_0 - ta_xxz_xzzz_1[i] * fe_0 + ta_xxyz_xzzz_0[i] * pa_y[i] - ta_xxyz_xzzz_1[i] * pc_y[i];

        ta_xxyyz_yyyy_0[i] = ta_xxyy_yyyy_0[i] * pa_z[i] - ta_xxyy_yyyy_1[i] * pc_z[i];

        ta_xxyyz_yyyz_0[i] = ta_yyz_yyyz_0[i] * fe_0 - ta_yyz_yyyz_1[i] * fe_0 + ta_xyyz_yyyz_0[i] * pa_x[i] - ta_xyyz_yyyz_1[i] * pc_x[i];

        ta_xxyyz_yyzz_0[i] = ta_yyz_yyzz_0[i] * fe_0 - ta_yyz_yyzz_1[i] * fe_0 + ta_xyyz_yyzz_0[i] * pa_x[i] - ta_xyyz_yyzz_1[i] * pc_x[i];

        ta_xxyyz_yzzz_0[i] = ta_yyz_yzzz_0[i] * fe_0 - ta_yyz_yzzz_1[i] * fe_0 + ta_xyyz_yzzz_0[i] * pa_x[i] - ta_xyyz_yzzz_1[i] * pc_x[i];

        ta_xxyyz_zzzz_0[i] = ta_yyz_zzzz_0[i] * fe_0 - ta_yyz_zzzz_1[i] * fe_0 + ta_xyyz_zzzz_0[i] * pa_x[i] - ta_xyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 120-135 components of targeted buffer : HG

    auto ta_xxyzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 120);

    auto ta_xxyzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 121);

    auto ta_xxyzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 122);

    auto ta_xxyzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 123);

    auto ta_xxyzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 124);

    auto ta_xxyzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 125);

    auto ta_xxyzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 126);

    auto ta_xxyzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 127);

    auto ta_xxyzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 128);

    auto ta_xxyzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 129);

    auto ta_xxyzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 130);

    auto ta_xxyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 131);

    auto ta_xxyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 132);

    auto ta_xxyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 133);

    auto ta_xxyzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 134);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xxyzz_xxxx_0, \
                             ta_xxyzz_xxxy_0, \
                             ta_xxyzz_xxxz_0, \
                             ta_xxyzz_xxyy_0, \
                             ta_xxyzz_xxyz_0, \
                             ta_xxyzz_xxzz_0, \
                             ta_xxyzz_xyyy_0, \
                             ta_xxyzz_xyyz_0, \
                             ta_xxyzz_xyzz_0, \
                             ta_xxyzz_xzzz_0, \
                             ta_xxyzz_yyyy_0, \
                             ta_xxyzz_yyyz_0, \
                             ta_xxyzz_yyzz_0, \
                             ta_xxyzz_yzzz_0, \
                             ta_xxyzz_zzzz_0, \
                             ta_xxzz_xxx_0,   \
                             ta_xxzz_xxx_1,   \
                             ta_xxzz_xxxx_0,  \
                             ta_xxzz_xxxx_1,  \
                             ta_xxzz_xxxy_0,  \
                             ta_xxzz_xxxy_1,  \
                             ta_xxzz_xxxz_0,  \
                             ta_xxzz_xxxz_1,  \
                             ta_xxzz_xxy_0,   \
                             ta_xxzz_xxy_1,   \
                             ta_xxzz_xxyy_0,  \
                             ta_xxzz_xxyy_1,  \
                             ta_xxzz_xxyz_0,  \
                             ta_xxzz_xxyz_1,  \
                             ta_xxzz_xxz_0,   \
                             ta_xxzz_xxz_1,   \
                             ta_xxzz_xxzz_0,  \
                             ta_xxzz_xxzz_1,  \
                             ta_xxzz_xyy_0,   \
                             ta_xxzz_xyy_1,   \
                             ta_xxzz_xyyy_0,  \
                             ta_xxzz_xyyy_1,  \
                             ta_xxzz_xyyz_0,  \
                             ta_xxzz_xyyz_1,  \
                             ta_xxzz_xyz_0,   \
                             ta_xxzz_xyz_1,   \
                             ta_xxzz_xyzz_0,  \
                             ta_xxzz_xyzz_1,  \
                             ta_xxzz_xzz_0,   \
                             ta_xxzz_xzz_1,   \
                             ta_xxzz_xzzz_0,  \
                             ta_xxzz_xzzz_1,  \
                             ta_xxzz_zzzz_0,  \
                             ta_xxzz_zzzz_1,  \
                             ta_xyzz_yyyy_0,  \
                             ta_xyzz_yyyy_1,  \
                             ta_xyzz_yyyz_0,  \
                             ta_xyzz_yyyz_1,  \
                             ta_xyzz_yyzz_0,  \
                             ta_xyzz_yyzz_1,  \
                             ta_xyzz_yzzz_0,  \
                             ta_xyzz_yzzz_1,  \
                             ta_yzz_yyyy_0,   \
                             ta_yzz_yyyy_1,   \
                             ta_yzz_yyyz_0,   \
                             ta_yzz_yyyz_1,   \
                             ta_yzz_yyzz_0,   \
                             ta_yzz_yyzz_1,   \
                             ta_yzz_yzzz_0,   \
                             ta_yzz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzz_xxxx_0[i] = ta_xxzz_xxxx_0[i] * pa_y[i] - ta_xxzz_xxxx_1[i] * pc_y[i];

        ta_xxyzz_xxxy_0[i] = ta_xxzz_xxx_0[i] * fe_0 - ta_xxzz_xxx_1[i] * fe_0 + ta_xxzz_xxxy_0[i] * pa_y[i] - ta_xxzz_xxxy_1[i] * pc_y[i];

        ta_xxyzz_xxxz_0[i] = ta_xxzz_xxxz_0[i] * pa_y[i] - ta_xxzz_xxxz_1[i] * pc_y[i];

        ta_xxyzz_xxyy_0[i] =
            2.0 * ta_xxzz_xxy_0[i] * fe_0 - 2.0 * ta_xxzz_xxy_1[i] * fe_0 + ta_xxzz_xxyy_0[i] * pa_y[i] - ta_xxzz_xxyy_1[i] * pc_y[i];

        ta_xxyzz_xxyz_0[i] = ta_xxzz_xxz_0[i] * fe_0 - ta_xxzz_xxz_1[i] * fe_0 + ta_xxzz_xxyz_0[i] * pa_y[i] - ta_xxzz_xxyz_1[i] * pc_y[i];

        ta_xxyzz_xxzz_0[i] = ta_xxzz_xxzz_0[i] * pa_y[i] - ta_xxzz_xxzz_1[i] * pc_y[i];

        ta_xxyzz_xyyy_0[i] =
            3.0 * ta_xxzz_xyy_0[i] * fe_0 - 3.0 * ta_xxzz_xyy_1[i] * fe_0 + ta_xxzz_xyyy_0[i] * pa_y[i] - ta_xxzz_xyyy_1[i] * pc_y[i];

        ta_xxyzz_xyyz_0[i] =
            2.0 * ta_xxzz_xyz_0[i] * fe_0 - 2.0 * ta_xxzz_xyz_1[i] * fe_0 + ta_xxzz_xyyz_0[i] * pa_y[i] - ta_xxzz_xyyz_1[i] * pc_y[i];

        ta_xxyzz_xyzz_0[i] = ta_xxzz_xzz_0[i] * fe_0 - ta_xxzz_xzz_1[i] * fe_0 + ta_xxzz_xyzz_0[i] * pa_y[i] - ta_xxzz_xyzz_1[i] * pc_y[i];

        ta_xxyzz_xzzz_0[i] = ta_xxzz_xzzz_0[i] * pa_y[i] - ta_xxzz_xzzz_1[i] * pc_y[i];

        ta_xxyzz_yyyy_0[i] = ta_yzz_yyyy_0[i] * fe_0 - ta_yzz_yyyy_1[i] * fe_0 + ta_xyzz_yyyy_0[i] * pa_x[i] - ta_xyzz_yyyy_1[i] * pc_x[i];

        ta_xxyzz_yyyz_0[i] = ta_yzz_yyyz_0[i] * fe_0 - ta_yzz_yyyz_1[i] * fe_0 + ta_xyzz_yyyz_0[i] * pa_x[i] - ta_xyzz_yyyz_1[i] * pc_x[i];

        ta_xxyzz_yyzz_0[i] = ta_yzz_yyzz_0[i] * fe_0 - ta_yzz_yyzz_1[i] * fe_0 + ta_xyzz_yyzz_0[i] * pa_x[i] - ta_xyzz_yyzz_1[i] * pc_x[i];

        ta_xxyzz_yzzz_0[i] = ta_yzz_yzzz_0[i] * fe_0 - ta_yzz_yzzz_1[i] * fe_0 + ta_xyzz_yzzz_0[i] * pa_x[i] - ta_xyzz_yzzz_1[i] * pc_x[i];

        ta_xxyzz_zzzz_0[i] = ta_xxzz_zzzz_0[i] * pa_y[i] - ta_xxzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : HG

    auto ta_xxzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 135);

    auto ta_xxzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 136);

    auto ta_xxzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 137);

    auto ta_xxzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 138);

    auto ta_xxzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 139);

    auto ta_xxzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 140);

    auto ta_xxzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 141);

    auto ta_xxzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 142);

    auto ta_xxzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 143);

    auto ta_xxzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 144);

    auto ta_xxzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 145);

    auto ta_xxzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 146);

    auto ta_xxzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 147);

    auto ta_xxzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 148);

    auto ta_xxzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 149);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xxz_xxxx_0,   \
                             ta_xxz_xxxx_1,   \
                             ta_xxz_xxxy_0,   \
                             ta_xxz_xxxy_1,   \
                             ta_xxz_xxyy_0,   \
                             ta_xxz_xxyy_1,   \
                             ta_xxz_xyyy_0,   \
                             ta_xxz_xyyy_1,   \
                             ta_xxzz_xxxx_0,  \
                             ta_xxzz_xxxx_1,  \
                             ta_xxzz_xxxy_0,  \
                             ta_xxzz_xxxy_1,  \
                             ta_xxzz_xxyy_0,  \
                             ta_xxzz_xxyy_1,  \
                             ta_xxzz_xyyy_0,  \
                             ta_xxzz_xyyy_1,  \
                             ta_xxzzz_xxxx_0, \
                             ta_xxzzz_xxxy_0, \
                             ta_xxzzz_xxxz_0, \
                             ta_xxzzz_xxyy_0, \
                             ta_xxzzz_xxyz_0, \
                             ta_xxzzz_xxzz_0, \
                             ta_xxzzz_xyyy_0, \
                             ta_xxzzz_xyyz_0, \
                             ta_xxzzz_xyzz_0, \
                             ta_xxzzz_xzzz_0, \
                             ta_xxzzz_yyyy_0, \
                             ta_xxzzz_yyyz_0, \
                             ta_xxzzz_yyzz_0, \
                             ta_xxzzz_yzzz_0, \
                             ta_xxzzz_zzzz_0, \
                             ta_xzzz_xxxz_0,  \
                             ta_xzzz_xxxz_1,  \
                             ta_xzzz_xxyz_0,  \
                             ta_xzzz_xxyz_1,  \
                             ta_xzzz_xxz_0,   \
                             ta_xzzz_xxz_1,   \
                             ta_xzzz_xxzz_0,  \
                             ta_xzzz_xxzz_1,  \
                             ta_xzzz_xyyz_0,  \
                             ta_xzzz_xyyz_1,  \
                             ta_xzzz_xyz_0,   \
                             ta_xzzz_xyz_1,   \
                             ta_xzzz_xyzz_0,  \
                             ta_xzzz_xyzz_1,  \
                             ta_xzzz_xzz_0,   \
                             ta_xzzz_xzz_1,   \
                             ta_xzzz_xzzz_0,  \
                             ta_xzzz_xzzz_1,  \
                             ta_xzzz_yyyy_0,  \
                             ta_xzzz_yyyy_1,  \
                             ta_xzzz_yyyz_0,  \
                             ta_xzzz_yyyz_1,  \
                             ta_xzzz_yyz_0,   \
                             ta_xzzz_yyz_1,   \
                             ta_xzzz_yyzz_0,  \
                             ta_xzzz_yyzz_1,  \
                             ta_xzzz_yzz_0,   \
                             ta_xzzz_yzz_1,   \
                             ta_xzzz_yzzz_0,  \
                             ta_xzzz_yzzz_1,  \
                             ta_xzzz_zzz_0,   \
                             ta_xzzz_zzz_1,   \
                             ta_xzzz_zzzz_0,  \
                             ta_xzzz_zzzz_1,  \
                             ta_zzz_xxxz_0,   \
                             ta_zzz_xxxz_1,   \
                             ta_zzz_xxyz_0,   \
                             ta_zzz_xxyz_1,   \
                             ta_zzz_xxzz_0,   \
                             ta_zzz_xxzz_1,   \
                             ta_zzz_xyyz_0,   \
                             ta_zzz_xyyz_1,   \
                             ta_zzz_xyzz_0,   \
                             ta_zzz_xyzz_1,   \
                             ta_zzz_xzzz_0,   \
                             ta_zzz_xzzz_1,   \
                             ta_zzz_yyyy_0,   \
                             ta_zzz_yyyy_1,   \
                             ta_zzz_yyyz_0,   \
                             ta_zzz_yyyz_1,   \
                             ta_zzz_yyzz_0,   \
                             ta_zzz_yyzz_1,   \
                             ta_zzz_yzzz_0,   \
                             ta_zzz_yzzz_1,   \
                             ta_zzz_zzzz_0,   \
                             ta_zzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzz_xxxx_0[i] =
            2.0 * ta_xxz_xxxx_0[i] * fe_0 - 2.0 * ta_xxz_xxxx_1[i] * fe_0 + ta_xxzz_xxxx_0[i] * pa_z[i] - ta_xxzz_xxxx_1[i] * pc_z[i];

        ta_xxzzz_xxxy_0[i] =
            2.0 * ta_xxz_xxxy_0[i] * fe_0 - 2.0 * ta_xxz_xxxy_1[i] * fe_0 + ta_xxzz_xxxy_0[i] * pa_z[i] - ta_xxzz_xxxy_1[i] * pc_z[i];

        ta_xxzzz_xxxz_0[i] = ta_zzz_xxxz_0[i] * fe_0 - ta_zzz_xxxz_1[i] * fe_0 + 3.0 * ta_xzzz_xxz_0[i] * fe_0 - 3.0 * ta_xzzz_xxz_1[i] * fe_0 +
                             ta_xzzz_xxxz_0[i] * pa_x[i] - ta_xzzz_xxxz_1[i] * pc_x[i];

        ta_xxzzz_xxyy_0[i] =
            2.0 * ta_xxz_xxyy_0[i] * fe_0 - 2.0 * ta_xxz_xxyy_1[i] * fe_0 + ta_xxzz_xxyy_0[i] * pa_z[i] - ta_xxzz_xxyy_1[i] * pc_z[i];

        ta_xxzzz_xxyz_0[i] = ta_zzz_xxyz_0[i] * fe_0 - ta_zzz_xxyz_1[i] * fe_0 + 2.0 * ta_xzzz_xyz_0[i] * fe_0 - 2.0 * ta_xzzz_xyz_1[i] * fe_0 +
                             ta_xzzz_xxyz_0[i] * pa_x[i] - ta_xzzz_xxyz_1[i] * pc_x[i];

        ta_xxzzz_xxzz_0[i] = ta_zzz_xxzz_0[i] * fe_0 - ta_zzz_xxzz_1[i] * fe_0 + 2.0 * ta_xzzz_xzz_0[i] * fe_0 - 2.0 * ta_xzzz_xzz_1[i] * fe_0 +
                             ta_xzzz_xxzz_0[i] * pa_x[i] - ta_xzzz_xxzz_1[i] * pc_x[i];

        ta_xxzzz_xyyy_0[i] =
            2.0 * ta_xxz_xyyy_0[i] * fe_0 - 2.0 * ta_xxz_xyyy_1[i] * fe_0 + ta_xxzz_xyyy_0[i] * pa_z[i] - ta_xxzz_xyyy_1[i] * pc_z[i];

        ta_xxzzz_xyyz_0[i] = ta_zzz_xyyz_0[i] * fe_0 - ta_zzz_xyyz_1[i] * fe_0 + ta_xzzz_yyz_0[i] * fe_0 - ta_xzzz_yyz_1[i] * fe_0 +
                             ta_xzzz_xyyz_0[i] * pa_x[i] - ta_xzzz_xyyz_1[i] * pc_x[i];

        ta_xxzzz_xyzz_0[i] = ta_zzz_xyzz_0[i] * fe_0 - ta_zzz_xyzz_1[i] * fe_0 + ta_xzzz_yzz_0[i] * fe_0 - ta_xzzz_yzz_1[i] * fe_0 +
                             ta_xzzz_xyzz_0[i] * pa_x[i] - ta_xzzz_xyzz_1[i] * pc_x[i];

        ta_xxzzz_xzzz_0[i] = ta_zzz_xzzz_0[i] * fe_0 - ta_zzz_xzzz_1[i] * fe_0 + ta_xzzz_zzz_0[i] * fe_0 - ta_xzzz_zzz_1[i] * fe_0 +
                             ta_xzzz_xzzz_0[i] * pa_x[i] - ta_xzzz_xzzz_1[i] * pc_x[i];

        ta_xxzzz_yyyy_0[i] = ta_zzz_yyyy_0[i] * fe_0 - ta_zzz_yyyy_1[i] * fe_0 + ta_xzzz_yyyy_0[i] * pa_x[i] - ta_xzzz_yyyy_1[i] * pc_x[i];

        ta_xxzzz_yyyz_0[i] = ta_zzz_yyyz_0[i] * fe_0 - ta_zzz_yyyz_1[i] * fe_0 + ta_xzzz_yyyz_0[i] * pa_x[i] - ta_xzzz_yyyz_1[i] * pc_x[i];

        ta_xxzzz_yyzz_0[i] = ta_zzz_yyzz_0[i] * fe_0 - ta_zzz_yyzz_1[i] * fe_0 + ta_xzzz_yyzz_0[i] * pa_x[i] - ta_xzzz_yyzz_1[i] * pc_x[i];

        ta_xxzzz_yzzz_0[i] = ta_zzz_yzzz_0[i] * fe_0 - ta_zzz_yzzz_1[i] * fe_0 + ta_xzzz_yzzz_0[i] * pa_x[i] - ta_xzzz_yzzz_1[i] * pc_x[i];

        ta_xxzzz_zzzz_0[i] = ta_zzz_zzzz_0[i] * fe_0 - ta_zzz_zzzz_1[i] * fe_0 + ta_xzzz_zzzz_0[i] * pa_x[i] - ta_xzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : HG

    auto ta_xyyyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 150);

    auto ta_xyyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 151);

    auto ta_xyyyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 152);

    auto ta_xyyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 153);

    auto ta_xyyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 154);

    auto ta_xyyyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 155);

    auto ta_xyyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 156);

    auto ta_xyyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 157);

    auto ta_xyyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 158);

    auto ta_xyyyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 159);

    auto ta_xyyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 160);

    auto ta_xyyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 161);

    auto ta_xyyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 162);

    auto ta_xyyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 163);

    auto ta_xyyyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 164);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xyyyy_xxxx_0, \
                             ta_xyyyy_xxxy_0, \
                             ta_xyyyy_xxxz_0, \
                             ta_xyyyy_xxyy_0, \
                             ta_xyyyy_xxyz_0, \
                             ta_xyyyy_xxzz_0, \
                             ta_xyyyy_xyyy_0, \
                             ta_xyyyy_xyyz_0, \
                             ta_xyyyy_xyzz_0, \
                             ta_xyyyy_xzzz_0, \
                             ta_xyyyy_yyyy_0, \
                             ta_xyyyy_yyyz_0, \
                             ta_xyyyy_yyzz_0, \
                             ta_xyyyy_yzzz_0, \
                             ta_xyyyy_zzzz_0, \
                             ta_yyyy_xxx_0,   \
                             ta_yyyy_xxx_1,   \
                             ta_yyyy_xxxx_0,  \
                             ta_yyyy_xxxx_1,  \
                             ta_yyyy_xxxy_0,  \
                             ta_yyyy_xxxy_1,  \
                             ta_yyyy_xxxz_0,  \
                             ta_yyyy_xxxz_1,  \
                             ta_yyyy_xxy_0,   \
                             ta_yyyy_xxy_1,   \
                             ta_yyyy_xxyy_0,  \
                             ta_yyyy_xxyy_1,  \
                             ta_yyyy_xxyz_0,  \
                             ta_yyyy_xxyz_1,  \
                             ta_yyyy_xxz_0,   \
                             ta_yyyy_xxz_1,   \
                             ta_yyyy_xxzz_0,  \
                             ta_yyyy_xxzz_1,  \
                             ta_yyyy_xyy_0,   \
                             ta_yyyy_xyy_1,   \
                             ta_yyyy_xyyy_0,  \
                             ta_yyyy_xyyy_1,  \
                             ta_yyyy_xyyz_0,  \
                             ta_yyyy_xyyz_1,  \
                             ta_yyyy_xyz_0,   \
                             ta_yyyy_xyz_1,   \
                             ta_yyyy_xyzz_0,  \
                             ta_yyyy_xyzz_1,  \
                             ta_yyyy_xzz_0,   \
                             ta_yyyy_xzz_1,   \
                             ta_yyyy_xzzz_0,  \
                             ta_yyyy_xzzz_1,  \
                             ta_yyyy_yyy_0,   \
                             ta_yyyy_yyy_1,   \
                             ta_yyyy_yyyy_0,  \
                             ta_yyyy_yyyy_1,  \
                             ta_yyyy_yyyz_0,  \
                             ta_yyyy_yyyz_1,  \
                             ta_yyyy_yyz_0,   \
                             ta_yyyy_yyz_1,   \
                             ta_yyyy_yyzz_0,  \
                             ta_yyyy_yyzz_1,  \
                             ta_yyyy_yzz_0,   \
                             ta_yyyy_yzz_1,   \
                             ta_yyyy_yzzz_0,  \
                             ta_yyyy_yzzz_1,  \
                             ta_yyyy_zzz_0,   \
                             ta_yyyy_zzz_1,   \
                             ta_yyyy_zzzz_0,  \
                             ta_yyyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyy_xxxx_0[i] =
            4.0 * ta_yyyy_xxx_0[i] * fe_0 - 4.0 * ta_yyyy_xxx_1[i] * fe_0 + ta_yyyy_xxxx_0[i] * pa_x[i] - ta_yyyy_xxxx_1[i] * pc_x[i];

        ta_xyyyy_xxxy_0[i] =
            3.0 * ta_yyyy_xxy_0[i] * fe_0 - 3.0 * ta_yyyy_xxy_1[i] * fe_0 + ta_yyyy_xxxy_0[i] * pa_x[i] - ta_yyyy_xxxy_1[i] * pc_x[i];

        ta_xyyyy_xxxz_0[i] =
            3.0 * ta_yyyy_xxz_0[i] * fe_0 - 3.0 * ta_yyyy_xxz_1[i] * fe_0 + ta_yyyy_xxxz_0[i] * pa_x[i] - ta_yyyy_xxxz_1[i] * pc_x[i];

        ta_xyyyy_xxyy_0[i] =
            2.0 * ta_yyyy_xyy_0[i] * fe_0 - 2.0 * ta_yyyy_xyy_1[i] * fe_0 + ta_yyyy_xxyy_0[i] * pa_x[i] - ta_yyyy_xxyy_1[i] * pc_x[i];

        ta_xyyyy_xxyz_0[i] =
            2.0 * ta_yyyy_xyz_0[i] * fe_0 - 2.0 * ta_yyyy_xyz_1[i] * fe_0 + ta_yyyy_xxyz_0[i] * pa_x[i] - ta_yyyy_xxyz_1[i] * pc_x[i];

        ta_xyyyy_xxzz_0[i] =
            2.0 * ta_yyyy_xzz_0[i] * fe_0 - 2.0 * ta_yyyy_xzz_1[i] * fe_0 + ta_yyyy_xxzz_0[i] * pa_x[i] - ta_yyyy_xxzz_1[i] * pc_x[i];

        ta_xyyyy_xyyy_0[i] = ta_yyyy_yyy_0[i] * fe_0 - ta_yyyy_yyy_1[i] * fe_0 + ta_yyyy_xyyy_0[i] * pa_x[i] - ta_yyyy_xyyy_1[i] * pc_x[i];

        ta_xyyyy_xyyz_0[i] = ta_yyyy_yyz_0[i] * fe_0 - ta_yyyy_yyz_1[i] * fe_0 + ta_yyyy_xyyz_0[i] * pa_x[i] - ta_yyyy_xyyz_1[i] * pc_x[i];

        ta_xyyyy_xyzz_0[i] = ta_yyyy_yzz_0[i] * fe_0 - ta_yyyy_yzz_1[i] * fe_0 + ta_yyyy_xyzz_0[i] * pa_x[i] - ta_yyyy_xyzz_1[i] * pc_x[i];

        ta_xyyyy_xzzz_0[i] = ta_yyyy_zzz_0[i] * fe_0 - ta_yyyy_zzz_1[i] * fe_0 + ta_yyyy_xzzz_0[i] * pa_x[i] - ta_yyyy_xzzz_1[i] * pc_x[i];

        ta_xyyyy_yyyy_0[i] = ta_yyyy_yyyy_0[i] * pa_x[i] - ta_yyyy_yyyy_1[i] * pc_x[i];

        ta_xyyyy_yyyz_0[i] = ta_yyyy_yyyz_0[i] * pa_x[i] - ta_yyyy_yyyz_1[i] * pc_x[i];

        ta_xyyyy_yyzz_0[i] = ta_yyyy_yyzz_0[i] * pa_x[i] - ta_yyyy_yyzz_1[i] * pc_x[i];

        ta_xyyyy_yzzz_0[i] = ta_yyyy_yzzz_0[i] * pa_x[i] - ta_yyyy_yzzz_1[i] * pc_x[i];

        ta_xyyyy_zzzz_0[i] = ta_yyyy_zzzz_0[i] * pa_x[i] - ta_yyyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 165-180 components of targeted buffer : HG

    auto ta_xyyyz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 165);

    auto ta_xyyyz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 166);

    auto ta_xyyyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 167);

    auto ta_xyyyz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 168);

    auto ta_xyyyz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 169);

    auto ta_xyyyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 170);

    auto ta_xyyyz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 171);

    auto ta_xyyyz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 172);

    auto ta_xyyyz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 173);

    auto ta_xyyyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 174);

    auto ta_xyyyz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 175);

    auto ta_xyyyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 176);

    auto ta_xyyyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 177);

    auto ta_xyyyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 178);

    auto ta_xyyyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 179);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta_xyyy_xxxx_0,  \
                             ta_xyyy_xxxx_1,  \
                             ta_xyyy_xxxy_0,  \
                             ta_xyyy_xxxy_1,  \
                             ta_xyyy_xxyy_0,  \
                             ta_xyyy_xxyy_1,  \
                             ta_xyyy_xyyy_0,  \
                             ta_xyyy_xyyy_1,  \
                             ta_xyyyz_xxxx_0, \
                             ta_xyyyz_xxxy_0, \
                             ta_xyyyz_xxxz_0, \
                             ta_xyyyz_xxyy_0, \
                             ta_xyyyz_xxyz_0, \
                             ta_xyyyz_xxzz_0, \
                             ta_xyyyz_xyyy_0, \
                             ta_xyyyz_xyyz_0, \
                             ta_xyyyz_xyzz_0, \
                             ta_xyyyz_xzzz_0, \
                             ta_xyyyz_yyyy_0, \
                             ta_xyyyz_yyyz_0, \
                             ta_xyyyz_yyzz_0, \
                             ta_xyyyz_yzzz_0, \
                             ta_xyyyz_zzzz_0, \
                             ta_yyyz_xxxz_0,  \
                             ta_yyyz_xxxz_1,  \
                             ta_yyyz_xxyz_0,  \
                             ta_yyyz_xxyz_1,  \
                             ta_yyyz_xxz_0,   \
                             ta_yyyz_xxz_1,   \
                             ta_yyyz_xxzz_0,  \
                             ta_yyyz_xxzz_1,  \
                             ta_yyyz_xyyz_0,  \
                             ta_yyyz_xyyz_1,  \
                             ta_yyyz_xyz_0,   \
                             ta_yyyz_xyz_1,   \
                             ta_yyyz_xyzz_0,  \
                             ta_yyyz_xyzz_1,  \
                             ta_yyyz_xzz_0,   \
                             ta_yyyz_xzz_1,   \
                             ta_yyyz_xzzz_0,  \
                             ta_yyyz_xzzz_1,  \
                             ta_yyyz_yyyy_0,  \
                             ta_yyyz_yyyy_1,  \
                             ta_yyyz_yyyz_0,  \
                             ta_yyyz_yyyz_1,  \
                             ta_yyyz_yyz_0,   \
                             ta_yyyz_yyz_1,   \
                             ta_yyyz_yyzz_0,  \
                             ta_yyyz_yyzz_1,  \
                             ta_yyyz_yzz_0,   \
                             ta_yyyz_yzz_1,   \
                             ta_yyyz_yzzz_0,  \
                             ta_yyyz_yzzz_1,  \
                             ta_yyyz_zzz_0,   \
                             ta_yyyz_zzz_1,   \
                             ta_yyyz_zzzz_0,  \
                             ta_yyyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyz_xxxx_0[i] = ta_xyyy_xxxx_0[i] * pa_z[i] - ta_xyyy_xxxx_1[i] * pc_z[i];

        ta_xyyyz_xxxy_0[i] = ta_xyyy_xxxy_0[i] * pa_z[i] - ta_xyyy_xxxy_1[i] * pc_z[i];

        ta_xyyyz_xxxz_0[i] =
            3.0 * ta_yyyz_xxz_0[i] * fe_0 - 3.0 * ta_yyyz_xxz_1[i] * fe_0 + ta_yyyz_xxxz_0[i] * pa_x[i] - ta_yyyz_xxxz_1[i] * pc_x[i];

        ta_xyyyz_xxyy_0[i] = ta_xyyy_xxyy_0[i] * pa_z[i] - ta_xyyy_xxyy_1[i] * pc_z[i];

        ta_xyyyz_xxyz_0[i] =
            2.0 * ta_yyyz_xyz_0[i] * fe_0 - 2.0 * ta_yyyz_xyz_1[i] * fe_0 + ta_yyyz_xxyz_0[i] * pa_x[i] - ta_yyyz_xxyz_1[i] * pc_x[i];

        ta_xyyyz_xxzz_0[i] =
            2.0 * ta_yyyz_xzz_0[i] * fe_0 - 2.0 * ta_yyyz_xzz_1[i] * fe_0 + ta_yyyz_xxzz_0[i] * pa_x[i] - ta_yyyz_xxzz_1[i] * pc_x[i];

        ta_xyyyz_xyyy_0[i] = ta_xyyy_xyyy_0[i] * pa_z[i] - ta_xyyy_xyyy_1[i] * pc_z[i];

        ta_xyyyz_xyyz_0[i] = ta_yyyz_yyz_0[i] * fe_0 - ta_yyyz_yyz_1[i] * fe_0 + ta_yyyz_xyyz_0[i] * pa_x[i] - ta_yyyz_xyyz_1[i] * pc_x[i];

        ta_xyyyz_xyzz_0[i] = ta_yyyz_yzz_0[i] * fe_0 - ta_yyyz_yzz_1[i] * fe_0 + ta_yyyz_xyzz_0[i] * pa_x[i] - ta_yyyz_xyzz_1[i] * pc_x[i];

        ta_xyyyz_xzzz_0[i] = ta_yyyz_zzz_0[i] * fe_0 - ta_yyyz_zzz_1[i] * fe_0 + ta_yyyz_xzzz_0[i] * pa_x[i] - ta_yyyz_xzzz_1[i] * pc_x[i];

        ta_xyyyz_yyyy_0[i] = ta_yyyz_yyyy_0[i] * pa_x[i] - ta_yyyz_yyyy_1[i] * pc_x[i];

        ta_xyyyz_yyyz_0[i] = ta_yyyz_yyyz_0[i] * pa_x[i] - ta_yyyz_yyyz_1[i] * pc_x[i];

        ta_xyyyz_yyzz_0[i] = ta_yyyz_yyzz_0[i] * pa_x[i] - ta_yyyz_yyzz_1[i] * pc_x[i];

        ta_xyyyz_yzzz_0[i] = ta_yyyz_yzzz_0[i] * pa_x[i] - ta_yyyz_yzzz_1[i] * pc_x[i];

        ta_xyyyz_zzzz_0[i] = ta_yyyz_zzzz_0[i] * pa_x[i] - ta_yyyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 180-195 components of targeted buffer : HG

    auto ta_xyyzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 180);

    auto ta_xyyzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 181);

    auto ta_xyyzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 182);

    auto ta_xyyzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 183);

    auto ta_xyyzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 184);

    auto ta_xyyzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 185);

    auto ta_xyyzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 186);

    auto ta_xyyzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 187);

    auto ta_xyyzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 188);

    auto ta_xyyzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 189);

    auto ta_xyyzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 190);

    auto ta_xyyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 191);

    auto ta_xyyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 192);

    auto ta_xyyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 193);

    auto ta_xyyzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 194);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xyyzz_xxxx_0, \
                             ta_xyyzz_xxxy_0, \
                             ta_xyyzz_xxxz_0, \
                             ta_xyyzz_xxyy_0, \
                             ta_xyyzz_xxyz_0, \
                             ta_xyyzz_xxzz_0, \
                             ta_xyyzz_xyyy_0, \
                             ta_xyyzz_xyyz_0, \
                             ta_xyyzz_xyzz_0, \
                             ta_xyyzz_xzzz_0, \
                             ta_xyyzz_yyyy_0, \
                             ta_xyyzz_yyyz_0, \
                             ta_xyyzz_yyzz_0, \
                             ta_xyyzz_yzzz_0, \
                             ta_xyyzz_zzzz_0, \
                             ta_yyzz_xxx_0,   \
                             ta_yyzz_xxx_1,   \
                             ta_yyzz_xxxx_0,  \
                             ta_yyzz_xxxx_1,  \
                             ta_yyzz_xxxy_0,  \
                             ta_yyzz_xxxy_1,  \
                             ta_yyzz_xxxz_0,  \
                             ta_yyzz_xxxz_1,  \
                             ta_yyzz_xxy_0,   \
                             ta_yyzz_xxy_1,   \
                             ta_yyzz_xxyy_0,  \
                             ta_yyzz_xxyy_1,  \
                             ta_yyzz_xxyz_0,  \
                             ta_yyzz_xxyz_1,  \
                             ta_yyzz_xxz_0,   \
                             ta_yyzz_xxz_1,   \
                             ta_yyzz_xxzz_0,  \
                             ta_yyzz_xxzz_1,  \
                             ta_yyzz_xyy_0,   \
                             ta_yyzz_xyy_1,   \
                             ta_yyzz_xyyy_0,  \
                             ta_yyzz_xyyy_1,  \
                             ta_yyzz_xyyz_0,  \
                             ta_yyzz_xyyz_1,  \
                             ta_yyzz_xyz_0,   \
                             ta_yyzz_xyz_1,   \
                             ta_yyzz_xyzz_0,  \
                             ta_yyzz_xyzz_1,  \
                             ta_yyzz_xzz_0,   \
                             ta_yyzz_xzz_1,   \
                             ta_yyzz_xzzz_0,  \
                             ta_yyzz_xzzz_1,  \
                             ta_yyzz_yyy_0,   \
                             ta_yyzz_yyy_1,   \
                             ta_yyzz_yyyy_0,  \
                             ta_yyzz_yyyy_1,  \
                             ta_yyzz_yyyz_0,  \
                             ta_yyzz_yyyz_1,  \
                             ta_yyzz_yyz_0,   \
                             ta_yyzz_yyz_1,   \
                             ta_yyzz_yyzz_0,  \
                             ta_yyzz_yyzz_1,  \
                             ta_yyzz_yzz_0,   \
                             ta_yyzz_yzz_1,   \
                             ta_yyzz_yzzz_0,  \
                             ta_yyzz_yzzz_1,  \
                             ta_yyzz_zzz_0,   \
                             ta_yyzz_zzz_1,   \
                             ta_yyzz_zzzz_0,  \
                             ta_yyzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzz_xxxx_0[i] =
            4.0 * ta_yyzz_xxx_0[i] * fe_0 - 4.0 * ta_yyzz_xxx_1[i] * fe_0 + ta_yyzz_xxxx_0[i] * pa_x[i] - ta_yyzz_xxxx_1[i] * pc_x[i];

        ta_xyyzz_xxxy_0[i] =
            3.0 * ta_yyzz_xxy_0[i] * fe_0 - 3.0 * ta_yyzz_xxy_1[i] * fe_0 + ta_yyzz_xxxy_0[i] * pa_x[i] - ta_yyzz_xxxy_1[i] * pc_x[i];

        ta_xyyzz_xxxz_0[i] =
            3.0 * ta_yyzz_xxz_0[i] * fe_0 - 3.0 * ta_yyzz_xxz_1[i] * fe_0 + ta_yyzz_xxxz_0[i] * pa_x[i] - ta_yyzz_xxxz_1[i] * pc_x[i];

        ta_xyyzz_xxyy_0[i] =
            2.0 * ta_yyzz_xyy_0[i] * fe_0 - 2.0 * ta_yyzz_xyy_1[i] * fe_0 + ta_yyzz_xxyy_0[i] * pa_x[i] - ta_yyzz_xxyy_1[i] * pc_x[i];

        ta_xyyzz_xxyz_0[i] =
            2.0 * ta_yyzz_xyz_0[i] * fe_0 - 2.0 * ta_yyzz_xyz_1[i] * fe_0 + ta_yyzz_xxyz_0[i] * pa_x[i] - ta_yyzz_xxyz_1[i] * pc_x[i];

        ta_xyyzz_xxzz_0[i] =
            2.0 * ta_yyzz_xzz_0[i] * fe_0 - 2.0 * ta_yyzz_xzz_1[i] * fe_0 + ta_yyzz_xxzz_0[i] * pa_x[i] - ta_yyzz_xxzz_1[i] * pc_x[i];

        ta_xyyzz_xyyy_0[i] = ta_yyzz_yyy_0[i] * fe_0 - ta_yyzz_yyy_1[i] * fe_0 + ta_yyzz_xyyy_0[i] * pa_x[i] - ta_yyzz_xyyy_1[i] * pc_x[i];

        ta_xyyzz_xyyz_0[i] = ta_yyzz_yyz_0[i] * fe_0 - ta_yyzz_yyz_1[i] * fe_0 + ta_yyzz_xyyz_0[i] * pa_x[i] - ta_yyzz_xyyz_1[i] * pc_x[i];

        ta_xyyzz_xyzz_0[i] = ta_yyzz_yzz_0[i] * fe_0 - ta_yyzz_yzz_1[i] * fe_0 + ta_yyzz_xyzz_0[i] * pa_x[i] - ta_yyzz_xyzz_1[i] * pc_x[i];

        ta_xyyzz_xzzz_0[i] = ta_yyzz_zzz_0[i] * fe_0 - ta_yyzz_zzz_1[i] * fe_0 + ta_yyzz_xzzz_0[i] * pa_x[i] - ta_yyzz_xzzz_1[i] * pc_x[i];

        ta_xyyzz_yyyy_0[i] = ta_yyzz_yyyy_0[i] * pa_x[i] - ta_yyzz_yyyy_1[i] * pc_x[i];

        ta_xyyzz_yyyz_0[i] = ta_yyzz_yyyz_0[i] * pa_x[i] - ta_yyzz_yyyz_1[i] * pc_x[i];

        ta_xyyzz_yyzz_0[i] = ta_yyzz_yyzz_0[i] * pa_x[i] - ta_yyzz_yyzz_1[i] * pc_x[i];

        ta_xyyzz_yzzz_0[i] = ta_yyzz_yzzz_0[i] * pa_x[i] - ta_yyzz_yzzz_1[i] * pc_x[i];

        ta_xyyzz_zzzz_0[i] = ta_yyzz_zzzz_0[i] * pa_x[i] - ta_yyzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 195-210 components of targeted buffer : HG

    auto ta_xyzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 195);

    auto ta_xyzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 196);

    auto ta_xyzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 197);

    auto ta_xyzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 198);

    auto ta_xyzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 199);

    auto ta_xyzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 200);

    auto ta_xyzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 201);

    auto ta_xyzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 202);

    auto ta_xyzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 203);

    auto ta_xyzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 204);

    auto ta_xyzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 205);

    auto ta_xyzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 206);

    auto ta_xyzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 207);

    auto ta_xyzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 208);

    auto ta_xyzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 209);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta_xyzzz_xxxx_0, \
                             ta_xyzzz_xxxy_0, \
                             ta_xyzzz_xxxz_0, \
                             ta_xyzzz_xxyy_0, \
                             ta_xyzzz_xxyz_0, \
                             ta_xyzzz_xxzz_0, \
                             ta_xyzzz_xyyy_0, \
                             ta_xyzzz_xyyz_0, \
                             ta_xyzzz_xyzz_0, \
                             ta_xyzzz_xzzz_0, \
                             ta_xyzzz_yyyy_0, \
                             ta_xyzzz_yyyz_0, \
                             ta_xyzzz_yyzz_0, \
                             ta_xyzzz_yzzz_0, \
                             ta_xyzzz_zzzz_0, \
                             ta_xzzz_xxxx_0,  \
                             ta_xzzz_xxxx_1,  \
                             ta_xzzz_xxxz_0,  \
                             ta_xzzz_xxxz_1,  \
                             ta_xzzz_xxzz_0,  \
                             ta_xzzz_xxzz_1,  \
                             ta_xzzz_xzzz_0,  \
                             ta_xzzz_xzzz_1,  \
                             ta_yzzz_xxxy_0,  \
                             ta_yzzz_xxxy_1,  \
                             ta_yzzz_xxy_0,   \
                             ta_yzzz_xxy_1,   \
                             ta_yzzz_xxyy_0,  \
                             ta_yzzz_xxyy_1,  \
                             ta_yzzz_xxyz_0,  \
                             ta_yzzz_xxyz_1,  \
                             ta_yzzz_xyy_0,   \
                             ta_yzzz_xyy_1,   \
                             ta_yzzz_xyyy_0,  \
                             ta_yzzz_xyyy_1,  \
                             ta_yzzz_xyyz_0,  \
                             ta_yzzz_xyyz_1,  \
                             ta_yzzz_xyz_0,   \
                             ta_yzzz_xyz_1,   \
                             ta_yzzz_xyzz_0,  \
                             ta_yzzz_xyzz_1,  \
                             ta_yzzz_yyy_0,   \
                             ta_yzzz_yyy_1,   \
                             ta_yzzz_yyyy_0,  \
                             ta_yzzz_yyyy_1,  \
                             ta_yzzz_yyyz_0,  \
                             ta_yzzz_yyyz_1,  \
                             ta_yzzz_yyz_0,   \
                             ta_yzzz_yyz_1,   \
                             ta_yzzz_yyzz_0,  \
                             ta_yzzz_yyzz_1,  \
                             ta_yzzz_yzz_0,   \
                             ta_yzzz_yzz_1,   \
                             ta_yzzz_yzzz_0,  \
                             ta_yzzz_yzzz_1,  \
                             ta_yzzz_zzzz_0,  \
                             ta_yzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzz_xxxx_0[i] = ta_xzzz_xxxx_0[i] * pa_y[i] - ta_xzzz_xxxx_1[i] * pc_y[i];

        ta_xyzzz_xxxy_0[i] =
            3.0 * ta_yzzz_xxy_0[i] * fe_0 - 3.0 * ta_yzzz_xxy_1[i] * fe_0 + ta_yzzz_xxxy_0[i] * pa_x[i] - ta_yzzz_xxxy_1[i] * pc_x[i];

        ta_xyzzz_xxxz_0[i] = ta_xzzz_xxxz_0[i] * pa_y[i] - ta_xzzz_xxxz_1[i] * pc_y[i];

        ta_xyzzz_xxyy_0[i] =
            2.0 * ta_yzzz_xyy_0[i] * fe_0 - 2.0 * ta_yzzz_xyy_1[i] * fe_0 + ta_yzzz_xxyy_0[i] * pa_x[i] - ta_yzzz_xxyy_1[i] * pc_x[i];

        ta_xyzzz_xxyz_0[i] =
            2.0 * ta_yzzz_xyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyz_1[i] * fe_0 + ta_yzzz_xxyz_0[i] * pa_x[i] - ta_yzzz_xxyz_1[i] * pc_x[i];

        ta_xyzzz_xxzz_0[i] = ta_xzzz_xxzz_0[i] * pa_y[i] - ta_xzzz_xxzz_1[i] * pc_y[i];

        ta_xyzzz_xyyy_0[i] = ta_yzzz_yyy_0[i] * fe_0 - ta_yzzz_yyy_1[i] * fe_0 + ta_yzzz_xyyy_0[i] * pa_x[i] - ta_yzzz_xyyy_1[i] * pc_x[i];

        ta_xyzzz_xyyz_0[i] = ta_yzzz_yyz_0[i] * fe_0 - ta_yzzz_yyz_1[i] * fe_0 + ta_yzzz_xyyz_0[i] * pa_x[i] - ta_yzzz_xyyz_1[i] * pc_x[i];

        ta_xyzzz_xyzz_0[i] = ta_yzzz_yzz_0[i] * fe_0 - ta_yzzz_yzz_1[i] * fe_0 + ta_yzzz_xyzz_0[i] * pa_x[i] - ta_yzzz_xyzz_1[i] * pc_x[i];

        ta_xyzzz_xzzz_0[i] = ta_xzzz_xzzz_0[i] * pa_y[i] - ta_xzzz_xzzz_1[i] * pc_y[i];

        ta_xyzzz_yyyy_0[i] = ta_yzzz_yyyy_0[i] * pa_x[i] - ta_yzzz_yyyy_1[i] * pc_x[i];

        ta_xyzzz_yyyz_0[i] = ta_yzzz_yyyz_0[i] * pa_x[i] - ta_yzzz_yyyz_1[i] * pc_x[i];

        ta_xyzzz_yyzz_0[i] = ta_yzzz_yyzz_0[i] * pa_x[i] - ta_yzzz_yyzz_1[i] * pc_x[i];

        ta_xyzzz_yzzz_0[i] = ta_yzzz_yzzz_0[i] * pa_x[i] - ta_yzzz_yzzz_1[i] * pc_x[i];

        ta_xyzzz_zzzz_0[i] = ta_yzzz_zzzz_0[i] * pa_x[i] - ta_yzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 210-225 components of targeted buffer : HG

    auto ta_xzzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 210);

    auto ta_xzzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 211);

    auto ta_xzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 212);

    auto ta_xzzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 213);

    auto ta_xzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 214);

    auto ta_xzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 215);

    auto ta_xzzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 216);

    auto ta_xzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 217);

    auto ta_xzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 218);

    auto ta_xzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 219);

    auto ta_xzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 220);

    auto ta_xzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 221);

    auto ta_xzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 222);

    auto ta_xzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 223);

    auto ta_xzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 224);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta_xzzzz_xxxx_0, \
                             ta_xzzzz_xxxy_0, \
                             ta_xzzzz_xxxz_0, \
                             ta_xzzzz_xxyy_0, \
                             ta_xzzzz_xxyz_0, \
                             ta_xzzzz_xxzz_0, \
                             ta_xzzzz_xyyy_0, \
                             ta_xzzzz_xyyz_0, \
                             ta_xzzzz_xyzz_0, \
                             ta_xzzzz_xzzz_0, \
                             ta_xzzzz_yyyy_0, \
                             ta_xzzzz_yyyz_0, \
                             ta_xzzzz_yyzz_0, \
                             ta_xzzzz_yzzz_0, \
                             ta_xzzzz_zzzz_0, \
                             ta_zzzz_xxx_0,   \
                             ta_zzzz_xxx_1,   \
                             ta_zzzz_xxxx_0,  \
                             ta_zzzz_xxxx_1,  \
                             ta_zzzz_xxxy_0,  \
                             ta_zzzz_xxxy_1,  \
                             ta_zzzz_xxxz_0,  \
                             ta_zzzz_xxxz_1,  \
                             ta_zzzz_xxy_0,   \
                             ta_zzzz_xxy_1,   \
                             ta_zzzz_xxyy_0,  \
                             ta_zzzz_xxyy_1,  \
                             ta_zzzz_xxyz_0,  \
                             ta_zzzz_xxyz_1,  \
                             ta_zzzz_xxz_0,   \
                             ta_zzzz_xxz_1,   \
                             ta_zzzz_xxzz_0,  \
                             ta_zzzz_xxzz_1,  \
                             ta_zzzz_xyy_0,   \
                             ta_zzzz_xyy_1,   \
                             ta_zzzz_xyyy_0,  \
                             ta_zzzz_xyyy_1,  \
                             ta_zzzz_xyyz_0,  \
                             ta_zzzz_xyyz_1,  \
                             ta_zzzz_xyz_0,   \
                             ta_zzzz_xyz_1,   \
                             ta_zzzz_xyzz_0,  \
                             ta_zzzz_xyzz_1,  \
                             ta_zzzz_xzz_0,   \
                             ta_zzzz_xzz_1,   \
                             ta_zzzz_xzzz_0,  \
                             ta_zzzz_xzzz_1,  \
                             ta_zzzz_yyy_0,   \
                             ta_zzzz_yyy_1,   \
                             ta_zzzz_yyyy_0,  \
                             ta_zzzz_yyyy_1,  \
                             ta_zzzz_yyyz_0,  \
                             ta_zzzz_yyyz_1,  \
                             ta_zzzz_yyz_0,   \
                             ta_zzzz_yyz_1,   \
                             ta_zzzz_yyzz_0,  \
                             ta_zzzz_yyzz_1,  \
                             ta_zzzz_yzz_0,   \
                             ta_zzzz_yzz_1,   \
                             ta_zzzz_yzzz_0,  \
                             ta_zzzz_yzzz_1,  \
                             ta_zzzz_zzz_0,   \
                             ta_zzzz_zzz_1,   \
                             ta_zzzz_zzzz_0,  \
                             ta_zzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzz_xxxx_0[i] =
            4.0 * ta_zzzz_xxx_0[i] * fe_0 - 4.0 * ta_zzzz_xxx_1[i] * fe_0 + ta_zzzz_xxxx_0[i] * pa_x[i] - ta_zzzz_xxxx_1[i] * pc_x[i];

        ta_xzzzz_xxxy_0[i] =
            3.0 * ta_zzzz_xxy_0[i] * fe_0 - 3.0 * ta_zzzz_xxy_1[i] * fe_0 + ta_zzzz_xxxy_0[i] * pa_x[i] - ta_zzzz_xxxy_1[i] * pc_x[i];

        ta_xzzzz_xxxz_0[i] =
            3.0 * ta_zzzz_xxz_0[i] * fe_0 - 3.0 * ta_zzzz_xxz_1[i] * fe_0 + ta_zzzz_xxxz_0[i] * pa_x[i] - ta_zzzz_xxxz_1[i] * pc_x[i];

        ta_xzzzz_xxyy_0[i] =
            2.0 * ta_zzzz_xyy_0[i] * fe_0 - 2.0 * ta_zzzz_xyy_1[i] * fe_0 + ta_zzzz_xxyy_0[i] * pa_x[i] - ta_zzzz_xxyy_1[i] * pc_x[i];

        ta_xzzzz_xxyz_0[i] =
            2.0 * ta_zzzz_xyz_0[i] * fe_0 - 2.0 * ta_zzzz_xyz_1[i] * fe_0 + ta_zzzz_xxyz_0[i] * pa_x[i] - ta_zzzz_xxyz_1[i] * pc_x[i];

        ta_xzzzz_xxzz_0[i] =
            2.0 * ta_zzzz_xzz_0[i] * fe_0 - 2.0 * ta_zzzz_xzz_1[i] * fe_0 + ta_zzzz_xxzz_0[i] * pa_x[i] - ta_zzzz_xxzz_1[i] * pc_x[i];

        ta_xzzzz_xyyy_0[i] = ta_zzzz_yyy_0[i] * fe_0 - ta_zzzz_yyy_1[i] * fe_0 + ta_zzzz_xyyy_0[i] * pa_x[i] - ta_zzzz_xyyy_1[i] * pc_x[i];

        ta_xzzzz_xyyz_0[i] = ta_zzzz_yyz_0[i] * fe_0 - ta_zzzz_yyz_1[i] * fe_0 + ta_zzzz_xyyz_0[i] * pa_x[i] - ta_zzzz_xyyz_1[i] * pc_x[i];

        ta_xzzzz_xyzz_0[i] = ta_zzzz_yzz_0[i] * fe_0 - ta_zzzz_yzz_1[i] * fe_0 + ta_zzzz_xyzz_0[i] * pa_x[i] - ta_zzzz_xyzz_1[i] * pc_x[i];

        ta_xzzzz_xzzz_0[i] = ta_zzzz_zzz_0[i] * fe_0 - ta_zzzz_zzz_1[i] * fe_0 + ta_zzzz_xzzz_0[i] * pa_x[i] - ta_zzzz_xzzz_1[i] * pc_x[i];

        ta_xzzzz_yyyy_0[i] = ta_zzzz_yyyy_0[i] * pa_x[i] - ta_zzzz_yyyy_1[i] * pc_x[i];

        ta_xzzzz_yyyz_0[i] = ta_zzzz_yyyz_0[i] * pa_x[i] - ta_zzzz_yyyz_1[i] * pc_x[i];

        ta_xzzzz_yyzz_0[i] = ta_zzzz_yyzz_0[i] * pa_x[i] - ta_zzzz_yyzz_1[i] * pc_x[i];

        ta_xzzzz_yzzz_0[i] = ta_zzzz_yzzz_0[i] * pa_x[i] - ta_zzzz_yzzz_1[i] * pc_x[i];

        ta_xzzzz_zzzz_0[i] = ta_zzzz_zzzz_0[i] * pa_x[i] - ta_zzzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : HG

    auto ta_yyyyy_xxxx_0 = pbuffer.data(idx_npot_0_hg + 225);

    auto ta_yyyyy_xxxy_0 = pbuffer.data(idx_npot_0_hg + 226);

    auto ta_yyyyy_xxxz_0 = pbuffer.data(idx_npot_0_hg + 227);

    auto ta_yyyyy_xxyy_0 = pbuffer.data(idx_npot_0_hg + 228);

    auto ta_yyyyy_xxyz_0 = pbuffer.data(idx_npot_0_hg + 229);

    auto ta_yyyyy_xxzz_0 = pbuffer.data(idx_npot_0_hg + 230);

    auto ta_yyyyy_xyyy_0 = pbuffer.data(idx_npot_0_hg + 231);

    auto ta_yyyyy_xyyz_0 = pbuffer.data(idx_npot_0_hg + 232);

    auto ta_yyyyy_xyzz_0 = pbuffer.data(idx_npot_0_hg + 233);

    auto ta_yyyyy_xzzz_0 = pbuffer.data(idx_npot_0_hg + 234);

    auto ta_yyyyy_yyyy_0 = pbuffer.data(idx_npot_0_hg + 235);

    auto ta_yyyyy_yyyz_0 = pbuffer.data(idx_npot_0_hg + 236);

    auto ta_yyyyy_yyzz_0 = pbuffer.data(idx_npot_0_hg + 237);

    auto ta_yyyyy_yzzz_0 = pbuffer.data(idx_npot_0_hg + 238);

    auto ta_yyyyy_zzzz_0 = pbuffer.data(idx_npot_0_hg + 239);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta_yyy_xxxx_0,   \
                             ta_yyy_xxxx_1,   \
                             ta_yyy_xxxy_0,   \
                             ta_yyy_xxxy_1,   \
                             ta_yyy_xxxz_0,   \
                             ta_yyy_xxxz_1,   \
                             ta_yyy_xxyy_0,   \
                             ta_yyy_xxyy_1,   \
                             ta_yyy_xxyz_0,   \
                             ta_yyy_xxyz_1,   \
                             ta_yyy_xxzz_0,   \
                             ta_yyy_xxzz_1,   \
                             ta_yyy_xyyy_0,   \
                             ta_yyy_xyyy_1,   \
                             ta_yyy_xyyz_0,   \
                             ta_yyy_xyyz_1,   \
                             ta_yyy_xyzz_0,   \
                             ta_yyy_xyzz_1,   \
                             ta_yyy_xzzz_0,   \
                             ta_yyy_xzzz_1,   \
                             ta_yyy_yyyy_0,   \
                             ta_yyy_yyyy_1,   \
                             ta_yyy_yyyz_0,   \
                             ta_yyy_yyyz_1,   \
                             ta_yyy_yyzz_0,   \
                             ta_yyy_yyzz_1,   \
                             ta_yyy_yzzz_0,   \
                             ta_yyy_yzzz_1,   \
                             ta_yyy_zzzz_0,   \
                             ta_yyy_zzzz_1,   \
                             ta_yyyy_xxx_0,   \
                             ta_yyyy_xxx_1,   \
                             ta_yyyy_xxxx_0,  \
                             ta_yyyy_xxxx_1,  \
                             ta_yyyy_xxxy_0,  \
                             ta_yyyy_xxxy_1,  \
                             ta_yyyy_xxxz_0,  \
                             ta_yyyy_xxxz_1,  \
                             ta_yyyy_xxy_0,   \
                             ta_yyyy_xxy_1,   \
                             ta_yyyy_xxyy_0,  \
                             ta_yyyy_xxyy_1,  \
                             ta_yyyy_xxyz_0,  \
                             ta_yyyy_xxyz_1,  \
                             ta_yyyy_xxz_0,   \
                             ta_yyyy_xxz_1,   \
                             ta_yyyy_xxzz_0,  \
                             ta_yyyy_xxzz_1,  \
                             ta_yyyy_xyy_0,   \
                             ta_yyyy_xyy_1,   \
                             ta_yyyy_xyyy_0,  \
                             ta_yyyy_xyyy_1,  \
                             ta_yyyy_xyyz_0,  \
                             ta_yyyy_xyyz_1,  \
                             ta_yyyy_xyz_0,   \
                             ta_yyyy_xyz_1,   \
                             ta_yyyy_xyzz_0,  \
                             ta_yyyy_xyzz_1,  \
                             ta_yyyy_xzz_0,   \
                             ta_yyyy_xzz_1,   \
                             ta_yyyy_xzzz_0,  \
                             ta_yyyy_xzzz_1,  \
                             ta_yyyy_yyy_0,   \
                             ta_yyyy_yyy_1,   \
                             ta_yyyy_yyyy_0,  \
                             ta_yyyy_yyyy_1,  \
                             ta_yyyy_yyyz_0,  \
                             ta_yyyy_yyyz_1,  \
                             ta_yyyy_yyz_0,   \
                             ta_yyyy_yyz_1,   \
                             ta_yyyy_yyzz_0,  \
                             ta_yyyy_yyzz_1,  \
                             ta_yyyy_yzz_0,   \
                             ta_yyyy_yzz_1,   \
                             ta_yyyy_yzzz_0,  \
                             ta_yyyy_yzzz_1,  \
                             ta_yyyy_zzz_0,   \
                             ta_yyyy_zzz_1,   \
                             ta_yyyy_zzzz_0,  \
                             ta_yyyy_zzzz_1,  \
                             ta_yyyyy_xxxx_0, \
                             ta_yyyyy_xxxy_0, \
                             ta_yyyyy_xxxz_0, \
                             ta_yyyyy_xxyy_0, \
                             ta_yyyyy_xxyz_0, \
                             ta_yyyyy_xxzz_0, \
                             ta_yyyyy_xyyy_0, \
                             ta_yyyyy_xyyz_0, \
                             ta_yyyyy_xyzz_0, \
                             ta_yyyyy_xzzz_0, \
                             ta_yyyyy_yyyy_0, \
                             ta_yyyyy_yyyz_0, \
                             ta_yyyyy_yyzz_0, \
                             ta_yyyyy_yzzz_0, \
                             ta_yyyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyy_xxxx_0[i] =
            4.0 * ta_yyy_xxxx_0[i] * fe_0 - 4.0 * ta_yyy_xxxx_1[i] * fe_0 + ta_yyyy_xxxx_0[i] * pa_y[i] - ta_yyyy_xxxx_1[i] * pc_y[i];

        ta_yyyyy_xxxy_0[i] = 4.0 * ta_yyy_xxxy_0[i] * fe_0 - 4.0 * ta_yyy_xxxy_1[i] * fe_0 + ta_yyyy_xxx_0[i] * fe_0 - ta_yyyy_xxx_1[i] * fe_0 +
                             ta_yyyy_xxxy_0[i] * pa_y[i] - ta_yyyy_xxxy_1[i] * pc_y[i];

        ta_yyyyy_xxxz_0[i] =
            4.0 * ta_yyy_xxxz_0[i] * fe_0 - 4.0 * ta_yyy_xxxz_1[i] * fe_0 + ta_yyyy_xxxz_0[i] * pa_y[i] - ta_yyyy_xxxz_1[i] * pc_y[i];

        ta_yyyyy_xxyy_0[i] = 4.0 * ta_yyy_xxyy_0[i] * fe_0 - 4.0 * ta_yyy_xxyy_1[i] * fe_0 + 2.0 * ta_yyyy_xxy_0[i] * fe_0 -
                             2.0 * ta_yyyy_xxy_1[i] * fe_0 + ta_yyyy_xxyy_0[i] * pa_y[i] - ta_yyyy_xxyy_1[i] * pc_y[i];

        ta_yyyyy_xxyz_0[i] = 4.0 * ta_yyy_xxyz_0[i] * fe_0 - 4.0 * ta_yyy_xxyz_1[i] * fe_0 + ta_yyyy_xxz_0[i] * fe_0 - ta_yyyy_xxz_1[i] * fe_0 +
                             ta_yyyy_xxyz_0[i] * pa_y[i] - ta_yyyy_xxyz_1[i] * pc_y[i];

        ta_yyyyy_xxzz_0[i] =
            4.0 * ta_yyy_xxzz_0[i] * fe_0 - 4.0 * ta_yyy_xxzz_1[i] * fe_0 + ta_yyyy_xxzz_0[i] * pa_y[i] - ta_yyyy_xxzz_1[i] * pc_y[i];

        ta_yyyyy_xyyy_0[i] = 4.0 * ta_yyy_xyyy_0[i] * fe_0 - 4.0 * ta_yyy_xyyy_1[i] * fe_0 + 3.0 * ta_yyyy_xyy_0[i] * fe_0 -
                             3.0 * ta_yyyy_xyy_1[i] * fe_0 + ta_yyyy_xyyy_0[i] * pa_y[i] - ta_yyyy_xyyy_1[i] * pc_y[i];

        ta_yyyyy_xyyz_0[i] = 4.0 * ta_yyy_xyyz_0[i] * fe_0 - 4.0 * ta_yyy_xyyz_1[i] * fe_0 + 2.0 * ta_yyyy_xyz_0[i] * fe_0 -
                             2.0 * ta_yyyy_xyz_1[i] * fe_0 + ta_yyyy_xyyz_0[i] * pa_y[i] - ta_yyyy_xyyz_1[i] * pc_y[i];

        ta_yyyyy_xyzz_0[i] = 4.0 * ta_yyy_xyzz_0[i] * fe_0 - 4.0 * ta_yyy_xyzz_1[i] * fe_0 + ta_yyyy_xzz_0[i] * fe_0 - ta_yyyy_xzz_1[i] * fe_0 +
                             ta_yyyy_xyzz_0[i] * pa_y[i] - ta_yyyy_xyzz_1[i] * pc_y[i];

        ta_yyyyy_xzzz_0[i] =
            4.0 * ta_yyy_xzzz_0[i] * fe_0 - 4.0 * ta_yyy_xzzz_1[i] * fe_0 + ta_yyyy_xzzz_0[i] * pa_y[i] - ta_yyyy_xzzz_1[i] * pc_y[i];

        ta_yyyyy_yyyy_0[i] = 4.0 * ta_yyy_yyyy_0[i] * fe_0 - 4.0 * ta_yyy_yyyy_1[i] * fe_0 + 4.0 * ta_yyyy_yyy_0[i] * fe_0 -
                             4.0 * ta_yyyy_yyy_1[i] * fe_0 + ta_yyyy_yyyy_0[i] * pa_y[i] - ta_yyyy_yyyy_1[i] * pc_y[i];

        ta_yyyyy_yyyz_0[i] = 4.0 * ta_yyy_yyyz_0[i] * fe_0 - 4.0 * ta_yyy_yyyz_1[i] * fe_0 + 3.0 * ta_yyyy_yyz_0[i] * fe_0 -
                             3.0 * ta_yyyy_yyz_1[i] * fe_0 + ta_yyyy_yyyz_0[i] * pa_y[i] - ta_yyyy_yyyz_1[i] * pc_y[i];

        ta_yyyyy_yyzz_0[i] = 4.0 * ta_yyy_yyzz_0[i] * fe_0 - 4.0 * ta_yyy_yyzz_1[i] * fe_0 + 2.0 * ta_yyyy_yzz_0[i] * fe_0 -
                             2.0 * ta_yyyy_yzz_1[i] * fe_0 + ta_yyyy_yyzz_0[i] * pa_y[i] - ta_yyyy_yyzz_1[i] * pc_y[i];

        ta_yyyyy_yzzz_0[i] = 4.0 * ta_yyy_yzzz_0[i] * fe_0 - 4.0 * ta_yyy_yzzz_1[i] * fe_0 + ta_yyyy_zzz_0[i] * fe_0 - ta_yyyy_zzz_1[i] * fe_0 +
                             ta_yyyy_yzzz_0[i] * pa_y[i] - ta_yyyy_yzzz_1[i] * pc_y[i];

        ta_yyyyy_zzzz_0[i] =
            4.0 * ta_yyy_zzzz_0[i] * fe_0 - 4.0 * ta_yyy_zzzz_1[i] * fe_0 + ta_yyyy_zzzz_0[i] * pa_y[i] - ta_yyyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 240-255 components of targeted buffer : HG

    auto ta_yyyyz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 240);

    auto ta_yyyyz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 241);

    auto ta_yyyyz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 242);

    auto ta_yyyyz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 243);

    auto ta_yyyyz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 244);

    auto ta_yyyyz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 245);

    auto ta_yyyyz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 246);

    auto ta_yyyyz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 247);

    auto ta_yyyyz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 248);

    auto ta_yyyyz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 249);

    auto ta_yyyyz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 250);

    auto ta_yyyyz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 251);

    auto ta_yyyyz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 252);

    auto ta_yyyyz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 253);

    auto ta_yyyyz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 254);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta_yyyy_xxxx_0,  \
                             ta_yyyy_xxxx_1,  \
                             ta_yyyy_xxxy_0,  \
                             ta_yyyy_xxxy_1,  \
                             ta_yyyy_xxy_0,   \
                             ta_yyyy_xxy_1,   \
                             ta_yyyy_xxyy_0,  \
                             ta_yyyy_xxyy_1,  \
                             ta_yyyy_xxyz_0,  \
                             ta_yyyy_xxyz_1,  \
                             ta_yyyy_xyy_0,   \
                             ta_yyyy_xyy_1,   \
                             ta_yyyy_xyyy_0,  \
                             ta_yyyy_xyyy_1,  \
                             ta_yyyy_xyyz_0,  \
                             ta_yyyy_xyyz_1,  \
                             ta_yyyy_xyz_0,   \
                             ta_yyyy_xyz_1,   \
                             ta_yyyy_xyzz_0,  \
                             ta_yyyy_xyzz_1,  \
                             ta_yyyy_yyy_0,   \
                             ta_yyyy_yyy_1,   \
                             ta_yyyy_yyyy_0,  \
                             ta_yyyy_yyyy_1,  \
                             ta_yyyy_yyyz_0,  \
                             ta_yyyy_yyyz_1,  \
                             ta_yyyy_yyz_0,   \
                             ta_yyyy_yyz_1,   \
                             ta_yyyy_yyzz_0,  \
                             ta_yyyy_yyzz_1,  \
                             ta_yyyy_yzz_0,   \
                             ta_yyyy_yzz_1,   \
                             ta_yyyy_yzzz_0,  \
                             ta_yyyy_yzzz_1,  \
                             ta_yyyyz_xxxx_0, \
                             ta_yyyyz_xxxy_0, \
                             ta_yyyyz_xxxz_0, \
                             ta_yyyyz_xxyy_0, \
                             ta_yyyyz_xxyz_0, \
                             ta_yyyyz_xxzz_0, \
                             ta_yyyyz_xyyy_0, \
                             ta_yyyyz_xyyz_0, \
                             ta_yyyyz_xyzz_0, \
                             ta_yyyyz_xzzz_0, \
                             ta_yyyyz_yyyy_0, \
                             ta_yyyyz_yyyz_0, \
                             ta_yyyyz_yyzz_0, \
                             ta_yyyyz_yzzz_0, \
                             ta_yyyyz_zzzz_0, \
                             ta_yyyz_xxxz_0,  \
                             ta_yyyz_xxxz_1,  \
                             ta_yyyz_xxzz_0,  \
                             ta_yyyz_xxzz_1,  \
                             ta_yyyz_xzzz_0,  \
                             ta_yyyz_xzzz_1,  \
                             ta_yyyz_zzzz_0,  \
                             ta_yyyz_zzzz_1,  \
                             ta_yyz_xxxz_0,   \
                             ta_yyz_xxxz_1,   \
                             ta_yyz_xxzz_0,   \
                             ta_yyz_xxzz_1,   \
                             ta_yyz_xzzz_0,   \
                             ta_yyz_xzzz_1,   \
                             ta_yyz_zzzz_0,   \
                             ta_yyz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyz_xxxx_0[i] = ta_yyyy_xxxx_0[i] * pa_z[i] - ta_yyyy_xxxx_1[i] * pc_z[i];

        ta_yyyyz_xxxy_0[i] = ta_yyyy_xxxy_0[i] * pa_z[i] - ta_yyyy_xxxy_1[i] * pc_z[i];

        ta_yyyyz_xxxz_0[i] =
            3.0 * ta_yyz_xxxz_0[i] * fe_0 - 3.0 * ta_yyz_xxxz_1[i] * fe_0 + ta_yyyz_xxxz_0[i] * pa_y[i] - ta_yyyz_xxxz_1[i] * pc_y[i];

        ta_yyyyz_xxyy_0[i] = ta_yyyy_xxyy_0[i] * pa_z[i] - ta_yyyy_xxyy_1[i] * pc_z[i];

        ta_yyyyz_xxyz_0[i] = ta_yyyy_xxy_0[i] * fe_0 - ta_yyyy_xxy_1[i] * fe_0 + ta_yyyy_xxyz_0[i] * pa_z[i] - ta_yyyy_xxyz_1[i] * pc_z[i];

        ta_yyyyz_xxzz_0[i] =
            3.0 * ta_yyz_xxzz_0[i] * fe_0 - 3.0 * ta_yyz_xxzz_1[i] * fe_0 + ta_yyyz_xxzz_0[i] * pa_y[i] - ta_yyyz_xxzz_1[i] * pc_y[i];

        ta_yyyyz_xyyy_0[i] = ta_yyyy_xyyy_0[i] * pa_z[i] - ta_yyyy_xyyy_1[i] * pc_z[i];

        ta_yyyyz_xyyz_0[i] = ta_yyyy_xyy_0[i] * fe_0 - ta_yyyy_xyy_1[i] * fe_0 + ta_yyyy_xyyz_0[i] * pa_z[i] - ta_yyyy_xyyz_1[i] * pc_z[i];

        ta_yyyyz_xyzz_0[i] =
            2.0 * ta_yyyy_xyz_0[i] * fe_0 - 2.0 * ta_yyyy_xyz_1[i] * fe_0 + ta_yyyy_xyzz_0[i] * pa_z[i] - ta_yyyy_xyzz_1[i] * pc_z[i];

        ta_yyyyz_xzzz_0[i] =
            3.0 * ta_yyz_xzzz_0[i] * fe_0 - 3.0 * ta_yyz_xzzz_1[i] * fe_0 + ta_yyyz_xzzz_0[i] * pa_y[i] - ta_yyyz_xzzz_1[i] * pc_y[i];

        ta_yyyyz_yyyy_0[i] = ta_yyyy_yyyy_0[i] * pa_z[i] - ta_yyyy_yyyy_1[i] * pc_z[i];

        ta_yyyyz_yyyz_0[i] = ta_yyyy_yyy_0[i] * fe_0 - ta_yyyy_yyy_1[i] * fe_0 + ta_yyyy_yyyz_0[i] * pa_z[i] - ta_yyyy_yyyz_1[i] * pc_z[i];

        ta_yyyyz_yyzz_0[i] =
            2.0 * ta_yyyy_yyz_0[i] * fe_0 - 2.0 * ta_yyyy_yyz_1[i] * fe_0 + ta_yyyy_yyzz_0[i] * pa_z[i] - ta_yyyy_yyzz_1[i] * pc_z[i];

        ta_yyyyz_yzzz_0[i] =
            3.0 * ta_yyyy_yzz_0[i] * fe_0 - 3.0 * ta_yyyy_yzz_1[i] * fe_0 + ta_yyyy_yzzz_0[i] * pa_z[i] - ta_yyyy_yzzz_1[i] * pc_z[i];

        ta_yyyyz_zzzz_0[i] =
            3.0 * ta_yyz_zzzz_0[i] * fe_0 - 3.0 * ta_yyz_zzzz_1[i] * fe_0 + ta_yyyz_zzzz_0[i] * pa_y[i] - ta_yyyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : HG

    auto ta_yyyzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 255);

    auto ta_yyyzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 256);

    auto ta_yyyzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 257);

    auto ta_yyyzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 258);

    auto ta_yyyzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 259);

    auto ta_yyyzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 260);

    auto ta_yyyzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 261);

    auto ta_yyyzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 262);

    auto ta_yyyzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 263);

    auto ta_yyyzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 264);

    auto ta_yyyzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 265);

    auto ta_yyyzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 266);

    auto ta_yyyzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 267);

    auto ta_yyyzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 268);

    auto ta_yyyzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 269);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta_yyy_xxxy_0,   \
                             ta_yyy_xxxy_1,   \
                             ta_yyy_xxyy_0,   \
                             ta_yyy_xxyy_1,   \
                             ta_yyy_xyyy_0,   \
                             ta_yyy_xyyy_1,   \
                             ta_yyy_yyyy_0,   \
                             ta_yyy_yyyy_1,   \
                             ta_yyyz_xxxy_0,  \
                             ta_yyyz_xxxy_1,  \
                             ta_yyyz_xxyy_0,  \
                             ta_yyyz_xxyy_1,  \
                             ta_yyyz_xyyy_0,  \
                             ta_yyyz_xyyy_1,  \
                             ta_yyyz_yyyy_0,  \
                             ta_yyyz_yyyy_1,  \
                             ta_yyyzz_xxxx_0, \
                             ta_yyyzz_xxxy_0, \
                             ta_yyyzz_xxxz_0, \
                             ta_yyyzz_xxyy_0, \
                             ta_yyyzz_xxyz_0, \
                             ta_yyyzz_xxzz_0, \
                             ta_yyyzz_xyyy_0, \
                             ta_yyyzz_xyyz_0, \
                             ta_yyyzz_xyzz_0, \
                             ta_yyyzz_xzzz_0, \
                             ta_yyyzz_yyyy_0, \
                             ta_yyyzz_yyyz_0, \
                             ta_yyyzz_yyzz_0, \
                             ta_yyyzz_yzzz_0, \
                             ta_yyyzz_zzzz_0, \
                             ta_yyzz_xxxx_0,  \
                             ta_yyzz_xxxx_1,  \
                             ta_yyzz_xxxz_0,  \
                             ta_yyzz_xxxz_1,  \
                             ta_yyzz_xxyz_0,  \
                             ta_yyzz_xxyz_1,  \
                             ta_yyzz_xxz_0,   \
                             ta_yyzz_xxz_1,   \
                             ta_yyzz_xxzz_0,  \
                             ta_yyzz_xxzz_1,  \
                             ta_yyzz_xyyz_0,  \
                             ta_yyzz_xyyz_1,  \
                             ta_yyzz_xyz_0,   \
                             ta_yyzz_xyz_1,   \
                             ta_yyzz_xyzz_0,  \
                             ta_yyzz_xyzz_1,  \
                             ta_yyzz_xzz_0,   \
                             ta_yyzz_xzz_1,   \
                             ta_yyzz_xzzz_0,  \
                             ta_yyzz_xzzz_1,  \
                             ta_yyzz_yyyz_0,  \
                             ta_yyzz_yyyz_1,  \
                             ta_yyzz_yyz_0,   \
                             ta_yyzz_yyz_1,   \
                             ta_yyzz_yyzz_0,  \
                             ta_yyzz_yyzz_1,  \
                             ta_yyzz_yzz_0,   \
                             ta_yyzz_yzz_1,   \
                             ta_yyzz_yzzz_0,  \
                             ta_yyzz_yzzz_1,  \
                             ta_yyzz_zzz_0,   \
                             ta_yyzz_zzz_1,   \
                             ta_yyzz_zzzz_0,  \
                             ta_yyzz_zzzz_1,  \
                             ta_yzz_xxxx_0,   \
                             ta_yzz_xxxx_1,   \
                             ta_yzz_xxxz_0,   \
                             ta_yzz_xxxz_1,   \
                             ta_yzz_xxyz_0,   \
                             ta_yzz_xxyz_1,   \
                             ta_yzz_xxzz_0,   \
                             ta_yzz_xxzz_1,   \
                             ta_yzz_xyyz_0,   \
                             ta_yzz_xyyz_1,   \
                             ta_yzz_xyzz_0,   \
                             ta_yzz_xyzz_1,   \
                             ta_yzz_xzzz_0,   \
                             ta_yzz_xzzz_1,   \
                             ta_yzz_yyyz_0,   \
                             ta_yzz_yyyz_1,   \
                             ta_yzz_yyzz_0,   \
                             ta_yzz_yyzz_1,   \
                             ta_yzz_yzzz_0,   \
                             ta_yzz_yzzz_1,   \
                             ta_yzz_zzzz_0,   \
                             ta_yzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzz_xxxx_0[i] =
            2.0 * ta_yzz_xxxx_0[i] * fe_0 - 2.0 * ta_yzz_xxxx_1[i] * fe_0 + ta_yyzz_xxxx_0[i] * pa_y[i] - ta_yyzz_xxxx_1[i] * pc_y[i];

        ta_yyyzz_xxxy_0[i] = ta_yyy_xxxy_0[i] * fe_0 - ta_yyy_xxxy_1[i] * fe_0 + ta_yyyz_xxxy_0[i] * pa_z[i] - ta_yyyz_xxxy_1[i] * pc_z[i];

        ta_yyyzz_xxxz_0[i] =
            2.0 * ta_yzz_xxxz_0[i] * fe_0 - 2.0 * ta_yzz_xxxz_1[i] * fe_0 + ta_yyzz_xxxz_0[i] * pa_y[i] - ta_yyzz_xxxz_1[i] * pc_y[i];

        ta_yyyzz_xxyy_0[i] = ta_yyy_xxyy_0[i] * fe_0 - ta_yyy_xxyy_1[i] * fe_0 + ta_yyyz_xxyy_0[i] * pa_z[i] - ta_yyyz_xxyy_1[i] * pc_z[i];

        ta_yyyzz_xxyz_0[i] = 2.0 * ta_yzz_xxyz_0[i] * fe_0 - 2.0 * ta_yzz_xxyz_1[i] * fe_0 + ta_yyzz_xxz_0[i] * fe_0 - ta_yyzz_xxz_1[i] * fe_0 +
                             ta_yyzz_xxyz_0[i] * pa_y[i] - ta_yyzz_xxyz_1[i] * pc_y[i];

        ta_yyyzz_xxzz_0[i] =
            2.0 * ta_yzz_xxzz_0[i] * fe_0 - 2.0 * ta_yzz_xxzz_1[i] * fe_0 + ta_yyzz_xxzz_0[i] * pa_y[i] - ta_yyzz_xxzz_1[i] * pc_y[i];

        ta_yyyzz_xyyy_0[i] = ta_yyy_xyyy_0[i] * fe_0 - ta_yyy_xyyy_1[i] * fe_0 + ta_yyyz_xyyy_0[i] * pa_z[i] - ta_yyyz_xyyy_1[i] * pc_z[i];

        ta_yyyzz_xyyz_0[i] = 2.0 * ta_yzz_xyyz_0[i] * fe_0 - 2.0 * ta_yzz_xyyz_1[i] * fe_0 + 2.0 * ta_yyzz_xyz_0[i] * fe_0 -
                             2.0 * ta_yyzz_xyz_1[i] * fe_0 + ta_yyzz_xyyz_0[i] * pa_y[i] - ta_yyzz_xyyz_1[i] * pc_y[i];

        ta_yyyzz_xyzz_0[i] = 2.0 * ta_yzz_xyzz_0[i] * fe_0 - 2.0 * ta_yzz_xyzz_1[i] * fe_0 + ta_yyzz_xzz_0[i] * fe_0 - ta_yyzz_xzz_1[i] * fe_0 +
                             ta_yyzz_xyzz_0[i] * pa_y[i] - ta_yyzz_xyzz_1[i] * pc_y[i];

        ta_yyyzz_xzzz_0[i] =
            2.0 * ta_yzz_xzzz_0[i] * fe_0 - 2.0 * ta_yzz_xzzz_1[i] * fe_0 + ta_yyzz_xzzz_0[i] * pa_y[i] - ta_yyzz_xzzz_1[i] * pc_y[i];

        ta_yyyzz_yyyy_0[i] = ta_yyy_yyyy_0[i] * fe_0 - ta_yyy_yyyy_1[i] * fe_0 + ta_yyyz_yyyy_0[i] * pa_z[i] - ta_yyyz_yyyy_1[i] * pc_z[i];

        ta_yyyzz_yyyz_0[i] = 2.0 * ta_yzz_yyyz_0[i] * fe_0 - 2.0 * ta_yzz_yyyz_1[i] * fe_0 + 3.0 * ta_yyzz_yyz_0[i] * fe_0 -
                             3.0 * ta_yyzz_yyz_1[i] * fe_0 + ta_yyzz_yyyz_0[i] * pa_y[i] - ta_yyzz_yyyz_1[i] * pc_y[i];

        ta_yyyzz_yyzz_0[i] = 2.0 * ta_yzz_yyzz_0[i] * fe_0 - 2.0 * ta_yzz_yyzz_1[i] * fe_0 + 2.0 * ta_yyzz_yzz_0[i] * fe_0 -
                             2.0 * ta_yyzz_yzz_1[i] * fe_0 + ta_yyzz_yyzz_0[i] * pa_y[i] - ta_yyzz_yyzz_1[i] * pc_y[i];

        ta_yyyzz_yzzz_0[i] = 2.0 * ta_yzz_yzzz_0[i] * fe_0 - 2.0 * ta_yzz_yzzz_1[i] * fe_0 + ta_yyzz_zzz_0[i] * fe_0 - ta_yyzz_zzz_1[i] * fe_0 +
                             ta_yyzz_yzzz_0[i] * pa_y[i] - ta_yyzz_yzzz_1[i] * pc_y[i];

        ta_yyyzz_zzzz_0[i] =
            2.0 * ta_yzz_zzzz_0[i] * fe_0 - 2.0 * ta_yzz_zzzz_1[i] * fe_0 + ta_yyzz_zzzz_0[i] * pa_y[i] - ta_yyzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 270-285 components of targeted buffer : HG

    auto ta_yyzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 270);

    auto ta_yyzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 271);

    auto ta_yyzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 272);

    auto ta_yyzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 273);

    auto ta_yyzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 274);

    auto ta_yyzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 275);

    auto ta_yyzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 276);

    auto ta_yyzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 277);

    auto ta_yyzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 278);

    auto ta_yyzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 279);

    auto ta_yyzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 280);

    auto ta_yyzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 281);

    auto ta_yyzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 282);

    auto ta_yyzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 283);

    auto ta_yyzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 284);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta_yyz_xxxy_0,   \
                             ta_yyz_xxxy_1,   \
                             ta_yyz_xxyy_0,   \
                             ta_yyz_xxyy_1,   \
                             ta_yyz_xyyy_0,   \
                             ta_yyz_xyyy_1,   \
                             ta_yyz_yyyy_0,   \
                             ta_yyz_yyyy_1,   \
                             ta_yyzz_xxxy_0,  \
                             ta_yyzz_xxxy_1,  \
                             ta_yyzz_xxyy_0,  \
                             ta_yyzz_xxyy_1,  \
                             ta_yyzz_xyyy_0,  \
                             ta_yyzz_xyyy_1,  \
                             ta_yyzz_yyyy_0,  \
                             ta_yyzz_yyyy_1,  \
                             ta_yyzzz_xxxx_0, \
                             ta_yyzzz_xxxy_0, \
                             ta_yyzzz_xxxz_0, \
                             ta_yyzzz_xxyy_0, \
                             ta_yyzzz_xxyz_0, \
                             ta_yyzzz_xxzz_0, \
                             ta_yyzzz_xyyy_0, \
                             ta_yyzzz_xyyz_0, \
                             ta_yyzzz_xyzz_0, \
                             ta_yyzzz_xzzz_0, \
                             ta_yyzzz_yyyy_0, \
                             ta_yyzzz_yyyz_0, \
                             ta_yyzzz_yyzz_0, \
                             ta_yyzzz_yzzz_0, \
                             ta_yyzzz_zzzz_0, \
                             ta_yzzz_xxxx_0,  \
                             ta_yzzz_xxxx_1,  \
                             ta_yzzz_xxxz_0,  \
                             ta_yzzz_xxxz_1,  \
                             ta_yzzz_xxyz_0,  \
                             ta_yzzz_xxyz_1,  \
                             ta_yzzz_xxz_0,   \
                             ta_yzzz_xxz_1,   \
                             ta_yzzz_xxzz_0,  \
                             ta_yzzz_xxzz_1,  \
                             ta_yzzz_xyyz_0,  \
                             ta_yzzz_xyyz_1,  \
                             ta_yzzz_xyz_0,   \
                             ta_yzzz_xyz_1,   \
                             ta_yzzz_xyzz_0,  \
                             ta_yzzz_xyzz_1,  \
                             ta_yzzz_xzz_0,   \
                             ta_yzzz_xzz_1,   \
                             ta_yzzz_xzzz_0,  \
                             ta_yzzz_xzzz_1,  \
                             ta_yzzz_yyyz_0,  \
                             ta_yzzz_yyyz_1,  \
                             ta_yzzz_yyz_0,   \
                             ta_yzzz_yyz_1,   \
                             ta_yzzz_yyzz_0,  \
                             ta_yzzz_yyzz_1,  \
                             ta_yzzz_yzz_0,   \
                             ta_yzzz_yzz_1,   \
                             ta_yzzz_yzzz_0,  \
                             ta_yzzz_yzzz_1,  \
                             ta_yzzz_zzz_0,   \
                             ta_yzzz_zzz_1,   \
                             ta_yzzz_zzzz_0,  \
                             ta_yzzz_zzzz_1,  \
                             ta_zzz_xxxx_0,   \
                             ta_zzz_xxxx_1,   \
                             ta_zzz_xxxz_0,   \
                             ta_zzz_xxxz_1,   \
                             ta_zzz_xxyz_0,   \
                             ta_zzz_xxyz_1,   \
                             ta_zzz_xxzz_0,   \
                             ta_zzz_xxzz_1,   \
                             ta_zzz_xyyz_0,   \
                             ta_zzz_xyyz_1,   \
                             ta_zzz_xyzz_0,   \
                             ta_zzz_xyzz_1,   \
                             ta_zzz_xzzz_0,   \
                             ta_zzz_xzzz_1,   \
                             ta_zzz_yyyz_0,   \
                             ta_zzz_yyyz_1,   \
                             ta_zzz_yyzz_0,   \
                             ta_zzz_yyzz_1,   \
                             ta_zzz_yzzz_0,   \
                             ta_zzz_yzzz_1,   \
                             ta_zzz_zzzz_0,   \
                             ta_zzz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzz_xxxx_0[i] = ta_zzz_xxxx_0[i] * fe_0 - ta_zzz_xxxx_1[i] * fe_0 + ta_yzzz_xxxx_0[i] * pa_y[i] - ta_yzzz_xxxx_1[i] * pc_y[i];

        ta_yyzzz_xxxy_0[i] =
            2.0 * ta_yyz_xxxy_0[i] * fe_0 - 2.0 * ta_yyz_xxxy_1[i] * fe_0 + ta_yyzz_xxxy_0[i] * pa_z[i] - ta_yyzz_xxxy_1[i] * pc_z[i];

        ta_yyzzz_xxxz_0[i] = ta_zzz_xxxz_0[i] * fe_0 - ta_zzz_xxxz_1[i] * fe_0 + ta_yzzz_xxxz_0[i] * pa_y[i] - ta_yzzz_xxxz_1[i] * pc_y[i];

        ta_yyzzz_xxyy_0[i] =
            2.0 * ta_yyz_xxyy_0[i] * fe_0 - 2.0 * ta_yyz_xxyy_1[i] * fe_0 + ta_yyzz_xxyy_0[i] * pa_z[i] - ta_yyzz_xxyy_1[i] * pc_z[i];

        ta_yyzzz_xxyz_0[i] = ta_zzz_xxyz_0[i] * fe_0 - ta_zzz_xxyz_1[i] * fe_0 + ta_yzzz_xxz_0[i] * fe_0 - ta_yzzz_xxz_1[i] * fe_0 +
                             ta_yzzz_xxyz_0[i] * pa_y[i] - ta_yzzz_xxyz_1[i] * pc_y[i];

        ta_yyzzz_xxzz_0[i] = ta_zzz_xxzz_0[i] * fe_0 - ta_zzz_xxzz_1[i] * fe_0 + ta_yzzz_xxzz_0[i] * pa_y[i] - ta_yzzz_xxzz_1[i] * pc_y[i];

        ta_yyzzz_xyyy_0[i] =
            2.0 * ta_yyz_xyyy_0[i] * fe_0 - 2.0 * ta_yyz_xyyy_1[i] * fe_0 + ta_yyzz_xyyy_0[i] * pa_z[i] - ta_yyzz_xyyy_1[i] * pc_z[i];

        ta_yyzzz_xyyz_0[i] = ta_zzz_xyyz_0[i] * fe_0 - ta_zzz_xyyz_1[i] * fe_0 + 2.0 * ta_yzzz_xyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyz_1[i] * fe_0 +
                             ta_yzzz_xyyz_0[i] * pa_y[i] - ta_yzzz_xyyz_1[i] * pc_y[i];

        ta_yyzzz_xyzz_0[i] = ta_zzz_xyzz_0[i] * fe_0 - ta_zzz_xyzz_1[i] * fe_0 + ta_yzzz_xzz_0[i] * fe_0 - ta_yzzz_xzz_1[i] * fe_0 +
                             ta_yzzz_xyzz_0[i] * pa_y[i] - ta_yzzz_xyzz_1[i] * pc_y[i];

        ta_yyzzz_xzzz_0[i] = ta_zzz_xzzz_0[i] * fe_0 - ta_zzz_xzzz_1[i] * fe_0 + ta_yzzz_xzzz_0[i] * pa_y[i] - ta_yzzz_xzzz_1[i] * pc_y[i];

        ta_yyzzz_yyyy_0[i] =
            2.0 * ta_yyz_yyyy_0[i] * fe_0 - 2.0 * ta_yyz_yyyy_1[i] * fe_0 + ta_yyzz_yyyy_0[i] * pa_z[i] - ta_yyzz_yyyy_1[i] * pc_z[i];

        ta_yyzzz_yyyz_0[i] = ta_zzz_yyyz_0[i] * fe_0 - ta_zzz_yyyz_1[i] * fe_0 + 3.0 * ta_yzzz_yyz_0[i] * fe_0 - 3.0 * ta_yzzz_yyz_1[i] * fe_0 +
                             ta_yzzz_yyyz_0[i] * pa_y[i] - ta_yzzz_yyyz_1[i] * pc_y[i];

        ta_yyzzz_yyzz_0[i] = ta_zzz_yyzz_0[i] * fe_0 - ta_zzz_yyzz_1[i] * fe_0 + 2.0 * ta_yzzz_yzz_0[i] * fe_0 - 2.0 * ta_yzzz_yzz_1[i] * fe_0 +
                             ta_yzzz_yyzz_0[i] * pa_y[i] - ta_yzzz_yyzz_1[i] * pc_y[i];

        ta_yyzzz_yzzz_0[i] = ta_zzz_yzzz_0[i] * fe_0 - ta_zzz_yzzz_1[i] * fe_0 + ta_yzzz_zzz_0[i] * fe_0 - ta_yzzz_zzz_1[i] * fe_0 +
                             ta_yzzz_yzzz_0[i] * pa_y[i] - ta_yzzz_yzzz_1[i] * pc_y[i];

        ta_yyzzz_zzzz_0[i] = ta_zzz_zzzz_0[i] * fe_0 - ta_zzz_zzzz_1[i] * fe_0 + ta_yzzz_zzzz_0[i] * pa_y[i] - ta_yzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 285-300 components of targeted buffer : HG

    auto ta_yzzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 285);

    auto ta_yzzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 286);

    auto ta_yzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 287);

    auto ta_yzzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 288);

    auto ta_yzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 289);

    auto ta_yzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 290);

    auto ta_yzzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 291);

    auto ta_yzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 292);

    auto ta_yzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 293);

    auto ta_yzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 294);

    auto ta_yzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 295);

    auto ta_yzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 296);

    auto ta_yzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 297);

    auto ta_yzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 298);

    auto ta_yzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 299);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta_yzzzz_xxxx_0, \
                             ta_yzzzz_xxxy_0, \
                             ta_yzzzz_xxxz_0, \
                             ta_yzzzz_xxyy_0, \
                             ta_yzzzz_xxyz_0, \
                             ta_yzzzz_xxzz_0, \
                             ta_yzzzz_xyyy_0, \
                             ta_yzzzz_xyyz_0, \
                             ta_yzzzz_xyzz_0, \
                             ta_yzzzz_xzzz_0, \
                             ta_yzzzz_yyyy_0, \
                             ta_yzzzz_yyyz_0, \
                             ta_yzzzz_yyzz_0, \
                             ta_yzzzz_yzzz_0, \
                             ta_yzzzz_zzzz_0, \
                             ta_zzzz_xxx_0,   \
                             ta_zzzz_xxx_1,   \
                             ta_zzzz_xxxx_0,  \
                             ta_zzzz_xxxx_1,  \
                             ta_zzzz_xxxy_0,  \
                             ta_zzzz_xxxy_1,  \
                             ta_zzzz_xxxz_0,  \
                             ta_zzzz_xxxz_1,  \
                             ta_zzzz_xxy_0,   \
                             ta_zzzz_xxy_1,   \
                             ta_zzzz_xxyy_0,  \
                             ta_zzzz_xxyy_1,  \
                             ta_zzzz_xxyz_0,  \
                             ta_zzzz_xxyz_1,  \
                             ta_zzzz_xxz_0,   \
                             ta_zzzz_xxz_1,   \
                             ta_zzzz_xxzz_0,  \
                             ta_zzzz_xxzz_1,  \
                             ta_zzzz_xyy_0,   \
                             ta_zzzz_xyy_1,   \
                             ta_zzzz_xyyy_0,  \
                             ta_zzzz_xyyy_1,  \
                             ta_zzzz_xyyz_0,  \
                             ta_zzzz_xyyz_1,  \
                             ta_zzzz_xyz_0,   \
                             ta_zzzz_xyz_1,   \
                             ta_zzzz_xyzz_0,  \
                             ta_zzzz_xyzz_1,  \
                             ta_zzzz_xzz_0,   \
                             ta_zzzz_xzz_1,   \
                             ta_zzzz_xzzz_0,  \
                             ta_zzzz_xzzz_1,  \
                             ta_zzzz_yyy_0,   \
                             ta_zzzz_yyy_1,   \
                             ta_zzzz_yyyy_0,  \
                             ta_zzzz_yyyy_1,  \
                             ta_zzzz_yyyz_0,  \
                             ta_zzzz_yyyz_1,  \
                             ta_zzzz_yyz_0,   \
                             ta_zzzz_yyz_1,   \
                             ta_zzzz_yyzz_0,  \
                             ta_zzzz_yyzz_1,  \
                             ta_zzzz_yzz_0,   \
                             ta_zzzz_yzz_1,   \
                             ta_zzzz_yzzz_0,  \
                             ta_zzzz_yzzz_1,  \
                             ta_zzzz_zzz_0,   \
                             ta_zzzz_zzz_1,   \
                             ta_zzzz_zzzz_0,  \
                             ta_zzzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzz_xxxx_0[i] = ta_zzzz_xxxx_0[i] * pa_y[i] - ta_zzzz_xxxx_1[i] * pc_y[i];

        ta_yzzzz_xxxy_0[i] = ta_zzzz_xxx_0[i] * fe_0 - ta_zzzz_xxx_1[i] * fe_0 + ta_zzzz_xxxy_0[i] * pa_y[i] - ta_zzzz_xxxy_1[i] * pc_y[i];

        ta_yzzzz_xxxz_0[i] = ta_zzzz_xxxz_0[i] * pa_y[i] - ta_zzzz_xxxz_1[i] * pc_y[i];

        ta_yzzzz_xxyy_0[i] =
            2.0 * ta_zzzz_xxy_0[i] * fe_0 - 2.0 * ta_zzzz_xxy_1[i] * fe_0 + ta_zzzz_xxyy_0[i] * pa_y[i] - ta_zzzz_xxyy_1[i] * pc_y[i];

        ta_yzzzz_xxyz_0[i] = ta_zzzz_xxz_0[i] * fe_0 - ta_zzzz_xxz_1[i] * fe_0 + ta_zzzz_xxyz_0[i] * pa_y[i] - ta_zzzz_xxyz_1[i] * pc_y[i];

        ta_yzzzz_xxzz_0[i] = ta_zzzz_xxzz_0[i] * pa_y[i] - ta_zzzz_xxzz_1[i] * pc_y[i];

        ta_yzzzz_xyyy_0[i] =
            3.0 * ta_zzzz_xyy_0[i] * fe_0 - 3.0 * ta_zzzz_xyy_1[i] * fe_0 + ta_zzzz_xyyy_0[i] * pa_y[i] - ta_zzzz_xyyy_1[i] * pc_y[i];

        ta_yzzzz_xyyz_0[i] =
            2.0 * ta_zzzz_xyz_0[i] * fe_0 - 2.0 * ta_zzzz_xyz_1[i] * fe_0 + ta_zzzz_xyyz_0[i] * pa_y[i] - ta_zzzz_xyyz_1[i] * pc_y[i];

        ta_yzzzz_xyzz_0[i] = ta_zzzz_xzz_0[i] * fe_0 - ta_zzzz_xzz_1[i] * fe_0 + ta_zzzz_xyzz_0[i] * pa_y[i] - ta_zzzz_xyzz_1[i] * pc_y[i];

        ta_yzzzz_xzzz_0[i] = ta_zzzz_xzzz_0[i] * pa_y[i] - ta_zzzz_xzzz_1[i] * pc_y[i];

        ta_yzzzz_yyyy_0[i] =
            4.0 * ta_zzzz_yyy_0[i] * fe_0 - 4.0 * ta_zzzz_yyy_1[i] * fe_0 + ta_zzzz_yyyy_0[i] * pa_y[i] - ta_zzzz_yyyy_1[i] * pc_y[i];

        ta_yzzzz_yyyz_0[i] =
            3.0 * ta_zzzz_yyz_0[i] * fe_0 - 3.0 * ta_zzzz_yyz_1[i] * fe_0 + ta_zzzz_yyyz_0[i] * pa_y[i] - ta_zzzz_yyyz_1[i] * pc_y[i];

        ta_yzzzz_yyzz_0[i] =
            2.0 * ta_zzzz_yzz_0[i] * fe_0 - 2.0 * ta_zzzz_yzz_1[i] * fe_0 + ta_zzzz_yyzz_0[i] * pa_y[i] - ta_zzzz_yyzz_1[i] * pc_y[i];

        ta_yzzzz_yzzz_0[i] = ta_zzzz_zzz_0[i] * fe_0 - ta_zzzz_zzz_1[i] * fe_0 + ta_zzzz_yzzz_0[i] * pa_y[i] - ta_zzzz_yzzz_1[i] * pc_y[i];

        ta_yzzzz_zzzz_0[i] = ta_zzzz_zzzz_0[i] * pa_y[i] - ta_zzzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 300-315 components of targeted buffer : HG

    auto ta_zzzzz_xxxx_0 = pbuffer.data(idx_npot_0_hg + 300);

    auto ta_zzzzz_xxxy_0 = pbuffer.data(idx_npot_0_hg + 301);

    auto ta_zzzzz_xxxz_0 = pbuffer.data(idx_npot_0_hg + 302);

    auto ta_zzzzz_xxyy_0 = pbuffer.data(idx_npot_0_hg + 303);

    auto ta_zzzzz_xxyz_0 = pbuffer.data(idx_npot_0_hg + 304);

    auto ta_zzzzz_xxzz_0 = pbuffer.data(idx_npot_0_hg + 305);

    auto ta_zzzzz_xyyy_0 = pbuffer.data(idx_npot_0_hg + 306);

    auto ta_zzzzz_xyyz_0 = pbuffer.data(idx_npot_0_hg + 307);

    auto ta_zzzzz_xyzz_0 = pbuffer.data(idx_npot_0_hg + 308);

    auto ta_zzzzz_xzzz_0 = pbuffer.data(idx_npot_0_hg + 309);

    auto ta_zzzzz_yyyy_0 = pbuffer.data(idx_npot_0_hg + 310);

    auto ta_zzzzz_yyyz_0 = pbuffer.data(idx_npot_0_hg + 311);

    auto ta_zzzzz_yyzz_0 = pbuffer.data(idx_npot_0_hg + 312);

    auto ta_zzzzz_yzzz_0 = pbuffer.data(idx_npot_0_hg + 313);

    auto ta_zzzzz_zzzz_0 = pbuffer.data(idx_npot_0_hg + 314);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta_zzz_xxxx_0,   \
                             ta_zzz_xxxx_1,   \
                             ta_zzz_xxxy_0,   \
                             ta_zzz_xxxy_1,   \
                             ta_zzz_xxxz_0,   \
                             ta_zzz_xxxz_1,   \
                             ta_zzz_xxyy_0,   \
                             ta_zzz_xxyy_1,   \
                             ta_zzz_xxyz_0,   \
                             ta_zzz_xxyz_1,   \
                             ta_zzz_xxzz_0,   \
                             ta_zzz_xxzz_1,   \
                             ta_zzz_xyyy_0,   \
                             ta_zzz_xyyy_1,   \
                             ta_zzz_xyyz_0,   \
                             ta_zzz_xyyz_1,   \
                             ta_zzz_xyzz_0,   \
                             ta_zzz_xyzz_1,   \
                             ta_zzz_xzzz_0,   \
                             ta_zzz_xzzz_1,   \
                             ta_zzz_yyyy_0,   \
                             ta_zzz_yyyy_1,   \
                             ta_zzz_yyyz_0,   \
                             ta_zzz_yyyz_1,   \
                             ta_zzz_yyzz_0,   \
                             ta_zzz_yyzz_1,   \
                             ta_zzz_yzzz_0,   \
                             ta_zzz_yzzz_1,   \
                             ta_zzz_zzzz_0,   \
                             ta_zzz_zzzz_1,   \
                             ta_zzzz_xxx_0,   \
                             ta_zzzz_xxx_1,   \
                             ta_zzzz_xxxx_0,  \
                             ta_zzzz_xxxx_1,  \
                             ta_zzzz_xxxy_0,  \
                             ta_zzzz_xxxy_1,  \
                             ta_zzzz_xxxz_0,  \
                             ta_zzzz_xxxz_1,  \
                             ta_zzzz_xxy_0,   \
                             ta_zzzz_xxy_1,   \
                             ta_zzzz_xxyy_0,  \
                             ta_zzzz_xxyy_1,  \
                             ta_zzzz_xxyz_0,  \
                             ta_zzzz_xxyz_1,  \
                             ta_zzzz_xxz_0,   \
                             ta_zzzz_xxz_1,   \
                             ta_zzzz_xxzz_0,  \
                             ta_zzzz_xxzz_1,  \
                             ta_zzzz_xyy_0,   \
                             ta_zzzz_xyy_1,   \
                             ta_zzzz_xyyy_0,  \
                             ta_zzzz_xyyy_1,  \
                             ta_zzzz_xyyz_0,  \
                             ta_zzzz_xyyz_1,  \
                             ta_zzzz_xyz_0,   \
                             ta_zzzz_xyz_1,   \
                             ta_zzzz_xyzz_0,  \
                             ta_zzzz_xyzz_1,  \
                             ta_zzzz_xzz_0,   \
                             ta_zzzz_xzz_1,   \
                             ta_zzzz_xzzz_0,  \
                             ta_zzzz_xzzz_1,  \
                             ta_zzzz_yyy_0,   \
                             ta_zzzz_yyy_1,   \
                             ta_zzzz_yyyy_0,  \
                             ta_zzzz_yyyy_1,  \
                             ta_zzzz_yyyz_0,  \
                             ta_zzzz_yyyz_1,  \
                             ta_zzzz_yyz_0,   \
                             ta_zzzz_yyz_1,   \
                             ta_zzzz_yyzz_0,  \
                             ta_zzzz_yyzz_1,  \
                             ta_zzzz_yzz_0,   \
                             ta_zzzz_yzz_1,   \
                             ta_zzzz_yzzz_0,  \
                             ta_zzzz_yzzz_1,  \
                             ta_zzzz_zzz_0,   \
                             ta_zzzz_zzz_1,   \
                             ta_zzzz_zzzz_0,  \
                             ta_zzzz_zzzz_1,  \
                             ta_zzzzz_xxxx_0, \
                             ta_zzzzz_xxxy_0, \
                             ta_zzzzz_xxxz_0, \
                             ta_zzzzz_xxyy_0, \
                             ta_zzzzz_xxyz_0, \
                             ta_zzzzz_xxzz_0, \
                             ta_zzzzz_xyyy_0, \
                             ta_zzzzz_xyyz_0, \
                             ta_zzzzz_xyzz_0, \
                             ta_zzzzz_xzzz_0, \
                             ta_zzzzz_yyyy_0, \
                             ta_zzzzz_yyyz_0, \
                             ta_zzzzz_yyzz_0, \
                             ta_zzzzz_yzzz_0, \
                             ta_zzzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzz_xxxx_0[i] =
            4.0 * ta_zzz_xxxx_0[i] * fe_0 - 4.0 * ta_zzz_xxxx_1[i] * fe_0 + ta_zzzz_xxxx_0[i] * pa_z[i] - ta_zzzz_xxxx_1[i] * pc_z[i];

        ta_zzzzz_xxxy_0[i] =
            4.0 * ta_zzz_xxxy_0[i] * fe_0 - 4.0 * ta_zzz_xxxy_1[i] * fe_0 + ta_zzzz_xxxy_0[i] * pa_z[i] - ta_zzzz_xxxy_1[i] * pc_z[i];

        ta_zzzzz_xxxz_0[i] = 4.0 * ta_zzz_xxxz_0[i] * fe_0 - 4.0 * ta_zzz_xxxz_1[i] * fe_0 + ta_zzzz_xxx_0[i] * fe_0 - ta_zzzz_xxx_1[i] * fe_0 +
                             ta_zzzz_xxxz_0[i] * pa_z[i] - ta_zzzz_xxxz_1[i] * pc_z[i];

        ta_zzzzz_xxyy_0[i] =
            4.0 * ta_zzz_xxyy_0[i] * fe_0 - 4.0 * ta_zzz_xxyy_1[i] * fe_0 + ta_zzzz_xxyy_0[i] * pa_z[i] - ta_zzzz_xxyy_1[i] * pc_z[i];

        ta_zzzzz_xxyz_0[i] = 4.0 * ta_zzz_xxyz_0[i] * fe_0 - 4.0 * ta_zzz_xxyz_1[i] * fe_0 + ta_zzzz_xxy_0[i] * fe_0 - ta_zzzz_xxy_1[i] * fe_0 +
                             ta_zzzz_xxyz_0[i] * pa_z[i] - ta_zzzz_xxyz_1[i] * pc_z[i];

        ta_zzzzz_xxzz_0[i] = 4.0 * ta_zzz_xxzz_0[i] * fe_0 - 4.0 * ta_zzz_xxzz_1[i] * fe_0 + 2.0 * ta_zzzz_xxz_0[i] * fe_0 -
                             2.0 * ta_zzzz_xxz_1[i] * fe_0 + ta_zzzz_xxzz_0[i] * pa_z[i] - ta_zzzz_xxzz_1[i] * pc_z[i];

        ta_zzzzz_xyyy_0[i] =
            4.0 * ta_zzz_xyyy_0[i] * fe_0 - 4.0 * ta_zzz_xyyy_1[i] * fe_0 + ta_zzzz_xyyy_0[i] * pa_z[i] - ta_zzzz_xyyy_1[i] * pc_z[i];

        ta_zzzzz_xyyz_0[i] = 4.0 * ta_zzz_xyyz_0[i] * fe_0 - 4.0 * ta_zzz_xyyz_1[i] * fe_0 + ta_zzzz_xyy_0[i] * fe_0 - ta_zzzz_xyy_1[i] * fe_0 +
                             ta_zzzz_xyyz_0[i] * pa_z[i] - ta_zzzz_xyyz_1[i] * pc_z[i];

        ta_zzzzz_xyzz_0[i] = 4.0 * ta_zzz_xyzz_0[i] * fe_0 - 4.0 * ta_zzz_xyzz_1[i] * fe_0 + 2.0 * ta_zzzz_xyz_0[i] * fe_0 -
                             2.0 * ta_zzzz_xyz_1[i] * fe_0 + ta_zzzz_xyzz_0[i] * pa_z[i] - ta_zzzz_xyzz_1[i] * pc_z[i];

        ta_zzzzz_xzzz_0[i] = 4.0 * ta_zzz_xzzz_0[i] * fe_0 - 4.0 * ta_zzz_xzzz_1[i] * fe_0 + 3.0 * ta_zzzz_xzz_0[i] * fe_0 -
                             3.0 * ta_zzzz_xzz_1[i] * fe_0 + ta_zzzz_xzzz_0[i] * pa_z[i] - ta_zzzz_xzzz_1[i] * pc_z[i];

        ta_zzzzz_yyyy_0[i] =
            4.0 * ta_zzz_yyyy_0[i] * fe_0 - 4.0 * ta_zzz_yyyy_1[i] * fe_0 + ta_zzzz_yyyy_0[i] * pa_z[i] - ta_zzzz_yyyy_1[i] * pc_z[i];

        ta_zzzzz_yyyz_0[i] = 4.0 * ta_zzz_yyyz_0[i] * fe_0 - 4.0 * ta_zzz_yyyz_1[i] * fe_0 + ta_zzzz_yyy_0[i] * fe_0 - ta_zzzz_yyy_1[i] * fe_0 +
                             ta_zzzz_yyyz_0[i] * pa_z[i] - ta_zzzz_yyyz_1[i] * pc_z[i];

        ta_zzzzz_yyzz_0[i] = 4.0 * ta_zzz_yyzz_0[i] * fe_0 - 4.0 * ta_zzz_yyzz_1[i] * fe_0 + 2.0 * ta_zzzz_yyz_0[i] * fe_0 -
                             2.0 * ta_zzzz_yyz_1[i] * fe_0 + ta_zzzz_yyzz_0[i] * pa_z[i] - ta_zzzz_yyzz_1[i] * pc_z[i];

        ta_zzzzz_yzzz_0[i] = 4.0 * ta_zzz_yzzz_0[i] * fe_0 - 4.0 * ta_zzz_yzzz_1[i] * fe_0 + 3.0 * ta_zzzz_yzz_0[i] * fe_0 -
                             3.0 * ta_zzzz_yzz_1[i] * fe_0 + ta_zzzz_yzzz_0[i] * pa_z[i] - ta_zzzz_yzzz_1[i] * pc_z[i];

        ta_zzzzz_zzzz_0[i] = 4.0 * ta_zzz_zzzz_0[i] * fe_0 - 4.0 * ta_zzz_zzzz_1[i] * fe_0 + 4.0 * ta_zzzz_zzz_0[i] * fe_0 -
                             4.0 * ta_zzzz_zzz_1[i] * fe_0 + ta_zzzz_zzzz_0[i] * pa_z[i] - ta_zzzz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
