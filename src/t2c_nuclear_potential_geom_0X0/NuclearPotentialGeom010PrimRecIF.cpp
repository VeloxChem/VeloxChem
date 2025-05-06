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

#include "NuclearPotentialGeom010PrimRecIF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_if(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_if,
                                        const size_t              idx_npot_geom_010_0_gf,
                                        const size_t              idx_npot_geom_010_1_gf,
                                        const size_t              idx_npot_geom_010_0_hd,
                                        const size_t              idx_npot_geom_010_1_hd,
                                        const size_t              idx_npot_1_hf,
                                        const size_t              idx_npot_geom_010_0_hf,
                                        const size_t              idx_npot_geom_010_1_hf,
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

    // Set up components of auxiliary buffer : GF

    auto ta1_x_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf);

    auto ta1_x_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 1);

    auto ta1_x_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 2);

    auto ta1_x_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 3);

    auto ta1_x_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 4);

    auto ta1_x_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 5);

    auto ta1_x_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 6);

    auto ta1_x_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 7);

    auto ta1_x_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 8);

    auto ta1_x_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 9);

    auto ta1_x_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 10);

    auto ta1_x_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 11);

    auto ta1_x_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 12);

    auto ta1_x_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 13);

    auto ta1_x_xxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 14);

    auto ta1_x_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 15);

    auto ta1_x_xxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 19);

    auto ta1_x_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 20);

    auto ta1_x_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 21);

    auto ta1_x_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 22);

    auto ta1_x_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 23);

    auto ta1_x_xxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 24);

    auto ta1_x_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 25);

    auto ta1_x_xxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 26);

    auto ta1_x_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 29);

    auto ta1_x_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 30);

    auto ta1_x_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 31);

    auto ta1_x_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 32);

    auto ta1_x_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 33);

    auto ta1_x_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 34);

    auto ta1_x_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 35);

    auto ta1_x_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 36);

    auto ta1_x_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 37);

    auto ta1_x_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 38);

    auto ta1_x_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 39);

    auto ta1_x_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 42);

    auto ta1_x_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 45);

    auto ta1_x_xxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 49);

    auto ta1_x_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 50);

    auto ta1_x_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 51);

    auto ta1_x_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 52);

    auto ta1_x_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 53);

    auto ta1_x_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 54);

    auto ta1_x_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 55);

    auto ta1_x_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 56);

    auto ta1_x_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 57);

    auto ta1_x_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 58);

    auto ta1_x_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 59);

    auto ta1_x_xyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 60);

    auto ta1_x_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 61);

    auto ta1_x_xyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 62);

    auto ta1_x_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 63);

    auto ta1_x_xyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 65);

    auto ta1_x_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 66);

    auto ta1_x_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 67);

    auto ta1_x_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 68);

    auto ta1_x_xyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 71);

    auto ta1_x_xyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 72);

    auto ta1_x_xyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 73);

    auto ta1_x_xyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 75);

    auto ta1_x_xyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 80);

    auto ta1_x_xyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 82);

    auto ta1_x_xyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 85);

    auto ta1_x_xzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 90);

    auto ta1_x_xzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 91);

    auto ta1_x_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 92);

    auto ta1_x_xzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 93);

    auto ta1_x_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 95);

    auto ta1_x_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 97);

    auto ta1_x_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 98);

    auto ta1_x_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 99);

    auto ta1_x_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 100);

    auto ta1_x_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 101);

    auto ta1_x_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 102);

    auto ta1_x_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 103);

    auto ta1_x_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 104);

    auto ta1_x_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 105);

    auto ta1_x_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 106);

    auto ta1_x_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 107);

    auto ta1_x_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 108);

    auto ta1_x_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 109);

    auto ta1_x_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 111);

    auto ta1_x_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 112);

    auto ta1_x_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 113);

    auto ta1_x_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 115);

    auto ta1_x_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 116);

    auto ta1_x_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 119);

    auto ta1_x_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 120);

    auto ta1_x_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 121);

    auto ta1_x_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 122);

    auto ta1_x_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 123);

    auto ta1_x_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 124);

    auto ta1_x_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 125);

    auto ta1_x_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 126);

    auto ta1_x_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 127);

    auto ta1_x_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 128);

    auto ta1_x_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 129);

    auto ta1_x_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 130);

    auto ta1_x_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 132);

    auto ta1_x_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 134);

    auto ta1_x_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 135);

    auto ta1_x_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 137);

    auto ta1_x_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 138);

    auto ta1_x_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 139);

    auto ta1_x_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 140);

    auto ta1_x_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 141);

    auto ta1_x_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 142);

    auto ta1_x_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 143);

    auto ta1_x_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 144);

    auto ta1_x_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 145);

    auto ta1_x_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 146);

    auto ta1_x_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 147);

    auto ta1_x_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 148);

    auto ta1_x_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 149);

    auto ta1_y_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 150);

    auto ta1_y_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 151);

    auto ta1_y_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 152);

    auto ta1_y_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 153);

    auto ta1_y_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 154);

    auto ta1_y_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 155);

    auto ta1_y_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 156);

    auto ta1_y_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 157);

    auto ta1_y_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 158);

    auto ta1_y_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 159);

    auto ta1_y_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 160);

    auto ta1_y_xxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 161);

    auto ta1_y_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 162);

    auto ta1_y_xxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 163);

    auto ta1_y_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 165);

    auto ta1_y_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 166);

    auto ta1_y_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 167);

    auto ta1_y_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 168);

    auto ta1_y_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 170);

    auto ta1_y_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 171);

    auto ta1_y_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 173);

    auto ta1_y_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 177);

    auto ta1_y_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 178);

    auto ta1_y_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 179);

    auto ta1_y_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 180);

    auto ta1_y_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 181);

    auto ta1_y_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 182);

    auto ta1_y_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 183);

    auto ta1_y_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 184);

    auto ta1_y_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 185);

    auto ta1_y_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 186);

    auto ta1_y_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 187);

    auto ta1_y_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 188);

    auto ta1_y_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 189);

    auto ta1_y_xxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 191);

    auto ta1_y_xxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 193);

    auto ta1_y_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 197);

    auto ta1_y_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 198);

    auto ta1_y_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 200);

    auto ta1_y_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 201);

    auto ta1_y_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 202);

    auto ta1_y_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 203);

    auto ta1_y_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 204);

    auto ta1_y_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 205);

    auto ta1_y_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 206);

    auto ta1_y_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 207);

    auto ta1_y_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 208);

    auto ta1_y_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 209);

    auto ta1_y_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 211);

    auto ta1_y_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 213);

    auto ta1_y_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 214);

    auto ta1_y_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 216);

    auto ta1_y_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 217);

    auto ta1_y_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 218);

    auto ta1_y_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 219);

    auto ta1_y_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 227);

    auto ta1_y_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 228);

    auto ta1_y_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 229);

    auto ta1_y_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 236);

    auto ta1_y_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 237);

    auto ta1_y_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 238);

    auto ta1_y_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 242);

    auto ta1_y_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 244);

    auto ta1_y_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 245);

    auto ta1_y_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 246);

    auto ta1_y_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 247);

    auto ta1_y_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 248);

    auto ta1_y_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 249);

    auto ta1_y_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 250);

    auto ta1_y_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 251);

    auto ta1_y_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 252);

    auto ta1_y_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 253);

    auto ta1_y_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 254);

    auto ta1_y_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 255);

    auto ta1_y_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 256);

    auto ta1_y_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 257);

    auto ta1_y_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 258);

    auto ta1_y_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 259);

    auto ta1_y_yyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 260);

    auto ta1_y_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 261);

    auto ta1_y_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 263);

    auto ta1_y_yyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 264);

    auto ta1_y_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 266);

    auto ta1_y_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 267);

    auto ta1_y_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 268);

    auto ta1_y_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 269);

    auto ta1_y_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 270);

    auto ta1_y_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 271);

    auto ta1_y_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 272);

    auto ta1_y_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 273);

    auto ta1_y_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 274);

    auto ta1_y_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 275);

    auto ta1_y_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 276);

    auto ta1_y_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 277);

    auto ta1_y_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 278);

    auto ta1_y_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 279);

    auto ta1_y_yzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 281);

    auto ta1_y_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 282);

    auto ta1_y_yzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 283);

    auto ta1_y_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 285);

    auto ta1_y_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 286);

    auto ta1_y_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 287);

    auto ta1_y_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 288);

    auto ta1_y_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 289);

    auto ta1_y_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 290);

    auto ta1_y_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 291);

    auto ta1_y_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 292);

    auto ta1_y_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 293);

    auto ta1_y_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 294);

    auto ta1_y_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 295);

    auto ta1_y_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 296);

    auto ta1_y_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 297);

    auto ta1_y_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 298);

    auto ta1_y_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 299);

    auto ta1_z_xxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 300);

    auto ta1_z_xxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 301);

    auto ta1_z_xxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 302);

    auto ta1_z_xxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 303);

    auto ta1_z_xxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 304);

    auto ta1_z_xxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 305);

    auto ta1_z_xxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 306);

    auto ta1_z_xxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 307);

    auto ta1_z_xxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 308);

    auto ta1_z_xxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 309);

    auto ta1_z_xxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 310);

    auto ta1_z_xxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 312);

    auto ta1_z_xxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 315);

    auto ta1_z_xxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 316);

    auto ta1_z_xxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 317);

    auto ta1_z_xxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 318);

    auto ta1_z_xxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 320);

    auto ta1_z_xxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 321);

    auto ta1_z_xxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 322);

    auto ta1_z_xxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 323);

    auto ta1_z_xxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 325);

    auto ta1_z_xxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 327);

    auto ta1_z_xxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 328);

    auto ta1_z_xxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 329);

    auto ta1_z_xxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 330);

    auto ta1_z_xxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 331);

    auto ta1_z_xxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 332);

    auto ta1_z_xxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 333);

    auto ta1_z_xxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 334);

    auto ta1_z_xxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 335);

    auto ta1_z_xxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 336);

    auto ta1_z_xxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 337);

    auto ta1_z_xxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 338);

    auto ta1_z_xxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 339);

    auto ta1_z_xxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 342);

    auto ta1_z_xxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 345);

    auto ta1_z_xxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 347);

    auto ta1_z_xxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 348);

    auto ta1_z_xxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 350);

    auto ta1_z_xxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 351);

    auto ta1_z_xxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 352);

    auto ta1_z_xxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 353);

    auto ta1_z_xxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 354);

    auto ta1_z_xxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 355);

    auto ta1_z_xxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 356);

    auto ta1_z_xxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 357);

    auto ta1_z_xxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 358);

    auto ta1_z_xxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 359);

    auto ta1_z_xyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 361);

    auto ta1_z_xyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 363);

    auto ta1_z_xyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 364);

    auto ta1_z_xyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 366);

    auto ta1_z_xyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 367);

    auto ta1_z_xyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 368);

    auto ta1_z_xyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 369);

    auto ta1_z_xyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 377);

    auto ta1_z_xyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 378);

    auto ta1_z_xyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 379);

    auto ta1_z_xyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 386);

    auto ta1_z_xyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 387);

    auto ta1_z_xyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 388);

    auto ta1_z_xzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 392);

    auto ta1_z_xzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 394);

    auto ta1_z_xzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 395);

    auto ta1_z_xzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 396);

    auto ta1_z_xzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 397);

    auto ta1_z_xzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 398);

    auto ta1_z_xzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 399);

    auto ta1_z_yyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 400);

    auto ta1_z_yyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 401);

    auto ta1_z_yyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 402);

    auto ta1_z_yyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 403);

    auto ta1_z_yyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 404);

    auto ta1_z_yyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 405);

    auto ta1_z_yyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 406);

    auto ta1_z_yyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 407);

    auto ta1_z_yyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 408);

    auto ta1_z_yyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 409);

    auto ta1_z_yyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 411);

    auto ta1_z_yyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 412);

    auto ta1_z_yyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 413);

    auto ta1_z_yyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 415);

    auto ta1_z_yyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 416);

    auto ta1_z_yyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 417);

    auto ta1_z_yyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 418);

    auto ta1_z_yyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 419);

    auto ta1_z_yyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 420);

    auto ta1_z_yyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 421);

    auto ta1_z_yyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 422);

    auto ta1_z_yyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 423);

    auto ta1_z_yyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 424);

    auto ta1_z_yyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 425);

    auto ta1_z_yyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 426);

    auto ta1_z_yyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 427);

    auto ta1_z_yyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 428);

    auto ta1_z_yyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 429);

    auto ta1_z_yzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 430);

    auto ta1_z_yzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 432);

    auto ta1_z_yzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 434);

    auto ta1_z_yzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 435);

    auto ta1_z_yzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 436);

    auto ta1_z_yzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 437);

    auto ta1_z_yzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 438);

    auto ta1_z_yzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 439);

    auto ta1_z_zzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_gf + 440);

    auto ta1_z_zzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 441);

    auto ta1_z_zzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 442);

    auto ta1_z_zzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 443);

    auto ta1_z_zzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 444);

    auto ta1_z_zzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 445);

    auto ta1_z_zzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_gf + 446);

    auto ta1_z_zzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 447);

    auto ta1_z_zzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 448);

    auto ta1_z_zzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_gf + 449);

    // Set up components of auxiliary buffer : GF

    auto ta1_x_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf);

    auto ta1_x_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 1);

    auto ta1_x_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 2);

    auto ta1_x_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 3);

    auto ta1_x_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 4);

    auto ta1_x_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 5);

    auto ta1_x_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 6);

    auto ta1_x_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 7);

    auto ta1_x_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 8);

    auto ta1_x_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 9);

    auto ta1_x_xxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 10);

    auto ta1_x_xxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 11);

    auto ta1_x_xxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 12);

    auto ta1_x_xxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 13);

    auto ta1_x_xxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 14);

    auto ta1_x_xxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 15);

    auto ta1_x_xxxy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 19);

    auto ta1_x_xxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 20);

    auto ta1_x_xxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 21);

    auto ta1_x_xxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 22);

    auto ta1_x_xxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 23);

    auto ta1_x_xxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 24);

    auto ta1_x_xxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 25);

    auto ta1_x_xxxz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 26);

    auto ta1_x_xxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 29);

    auto ta1_x_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 30);

    auto ta1_x_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 31);

    auto ta1_x_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 32);

    auto ta1_x_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 33);

    auto ta1_x_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 34);

    auto ta1_x_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 35);

    auto ta1_x_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 36);

    auto ta1_x_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 37);

    auto ta1_x_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 38);

    auto ta1_x_xxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 39);

    auto ta1_x_xxyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 42);

    auto ta1_x_xxyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 45);

    auto ta1_x_xxyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 49);

    auto ta1_x_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 50);

    auto ta1_x_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 51);

    auto ta1_x_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 52);

    auto ta1_x_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 53);

    auto ta1_x_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 54);

    auto ta1_x_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 55);

    auto ta1_x_xxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 56);

    auto ta1_x_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 57);

    auto ta1_x_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 58);

    auto ta1_x_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 59);

    auto ta1_x_xyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 60);

    auto ta1_x_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 61);

    auto ta1_x_xyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 62);

    auto ta1_x_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 63);

    auto ta1_x_xyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 65);

    auto ta1_x_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 66);

    auto ta1_x_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 67);

    auto ta1_x_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 68);

    auto ta1_x_xyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 71);

    auto ta1_x_xyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 72);

    auto ta1_x_xyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 73);

    auto ta1_x_xyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 75);

    auto ta1_x_xyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 80);

    auto ta1_x_xyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 82);

    auto ta1_x_xyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 85);

    auto ta1_x_xzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 90);

    auto ta1_x_xzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 91);

    auto ta1_x_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 92);

    auto ta1_x_xzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 93);

    auto ta1_x_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 95);

    auto ta1_x_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 97);

    auto ta1_x_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 98);

    auto ta1_x_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 99);

    auto ta1_x_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 100);

    auto ta1_x_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 101);

    auto ta1_x_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 102);

    auto ta1_x_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 103);

    auto ta1_x_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 104);

    auto ta1_x_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 105);

    auto ta1_x_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 106);

    auto ta1_x_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 107);

    auto ta1_x_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 108);

    auto ta1_x_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 109);

    auto ta1_x_yyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 111);

    auto ta1_x_yyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 112);

    auto ta1_x_yyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 113);

    auto ta1_x_yyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 115);

    auto ta1_x_yyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 116);

    auto ta1_x_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 119);

    auto ta1_x_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 120);

    auto ta1_x_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 121);

    auto ta1_x_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 122);

    auto ta1_x_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 123);

    auto ta1_x_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 124);

    auto ta1_x_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 125);

    auto ta1_x_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 126);

    auto ta1_x_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 127);

    auto ta1_x_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 128);

    auto ta1_x_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 129);

    auto ta1_x_yzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 130);

    auto ta1_x_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 132);

    auto ta1_x_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 134);

    auto ta1_x_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 135);

    auto ta1_x_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 137);

    auto ta1_x_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 138);

    auto ta1_x_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 139);

    auto ta1_x_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 140);

    auto ta1_x_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 141);

    auto ta1_x_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 142);

    auto ta1_x_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 143);

    auto ta1_x_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 144);

    auto ta1_x_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 145);

    auto ta1_x_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 146);

    auto ta1_x_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 147);

    auto ta1_x_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 148);

    auto ta1_x_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 149);

    auto ta1_y_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 150);

    auto ta1_y_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 151);

    auto ta1_y_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 152);

    auto ta1_y_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 153);

    auto ta1_y_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 154);

    auto ta1_y_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 155);

    auto ta1_y_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 156);

    auto ta1_y_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 157);

    auto ta1_y_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 158);

    auto ta1_y_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 159);

    auto ta1_y_xxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 160);

    auto ta1_y_xxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 161);

    auto ta1_y_xxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 162);

    auto ta1_y_xxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 163);

    auto ta1_y_xxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 165);

    auto ta1_y_xxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 166);

    auto ta1_y_xxxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 167);

    auto ta1_y_xxxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 168);

    auto ta1_y_xxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 170);

    auto ta1_y_xxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 171);

    auto ta1_y_xxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 173);

    auto ta1_y_xxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 177);

    auto ta1_y_xxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 178);

    auto ta1_y_xxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 179);

    auto ta1_y_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 180);

    auto ta1_y_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 181);

    auto ta1_y_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 182);

    auto ta1_y_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 183);

    auto ta1_y_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 184);

    auto ta1_y_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 185);

    auto ta1_y_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 186);

    auto ta1_y_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 187);

    auto ta1_y_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 188);

    auto ta1_y_xxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 189);

    auto ta1_y_xxyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 191);

    auto ta1_y_xxyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 193);

    auto ta1_y_xxyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 197);

    auto ta1_y_xxyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 198);

    auto ta1_y_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 200);

    auto ta1_y_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 201);

    auto ta1_y_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 202);

    auto ta1_y_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 203);

    auto ta1_y_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 204);

    auto ta1_y_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 205);

    auto ta1_y_xxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 206);

    auto ta1_y_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 207);

    auto ta1_y_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 208);

    auto ta1_y_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 209);

    auto ta1_y_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 211);

    auto ta1_y_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 213);

    auto ta1_y_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 214);

    auto ta1_y_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 216);

    auto ta1_y_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 217);

    auto ta1_y_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 218);

    auto ta1_y_xyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 219);

    auto ta1_y_xyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 227);

    auto ta1_y_xyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 228);

    auto ta1_y_xyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 229);

    auto ta1_y_xyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 236);

    auto ta1_y_xyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 237);

    auto ta1_y_xyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 238);

    auto ta1_y_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 242);

    auto ta1_y_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 244);

    auto ta1_y_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 245);

    auto ta1_y_xzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 246);

    auto ta1_y_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 247);

    auto ta1_y_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 248);

    auto ta1_y_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 249);

    auto ta1_y_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 250);

    auto ta1_y_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 251);

    auto ta1_y_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 252);

    auto ta1_y_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 253);

    auto ta1_y_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 254);

    auto ta1_y_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 255);

    auto ta1_y_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 256);

    auto ta1_y_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 257);

    auto ta1_y_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 258);

    auto ta1_y_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 259);

    auto ta1_y_yyyz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 260);

    auto ta1_y_yyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 261);

    auto ta1_y_yyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 263);

    auto ta1_y_yyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 264);

    auto ta1_y_yyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 266);

    auto ta1_y_yyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 267);

    auto ta1_y_yyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 268);

    auto ta1_y_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 269);

    auto ta1_y_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 270);

    auto ta1_y_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 271);

    auto ta1_y_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 272);

    auto ta1_y_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 273);

    auto ta1_y_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 274);

    auto ta1_y_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 275);

    auto ta1_y_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 276);

    auto ta1_y_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 277);

    auto ta1_y_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 278);

    auto ta1_y_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 279);

    auto ta1_y_yzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 281);

    auto ta1_y_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 282);

    auto ta1_y_yzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 283);

    auto ta1_y_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 285);

    auto ta1_y_yzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 286);

    auto ta1_y_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 287);

    auto ta1_y_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 288);

    auto ta1_y_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 289);

    auto ta1_y_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 290);

    auto ta1_y_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 291);

    auto ta1_y_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 292);

    auto ta1_y_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 293);

    auto ta1_y_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 294);

    auto ta1_y_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 295);

    auto ta1_y_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 296);

    auto ta1_y_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 297);

    auto ta1_y_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 298);

    auto ta1_y_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 299);

    auto ta1_z_xxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 300);

    auto ta1_z_xxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 301);

    auto ta1_z_xxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 302);

    auto ta1_z_xxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 303);

    auto ta1_z_xxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 304);

    auto ta1_z_xxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 305);

    auto ta1_z_xxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 306);

    auto ta1_z_xxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 307);

    auto ta1_z_xxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 308);

    auto ta1_z_xxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 309);

    auto ta1_z_xxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 310);

    auto ta1_z_xxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 312);

    auto ta1_z_xxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 315);

    auto ta1_z_xxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 316);

    auto ta1_z_xxxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 317);

    auto ta1_z_xxxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 318);

    auto ta1_z_xxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 320);

    auto ta1_z_xxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 321);

    auto ta1_z_xxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 322);

    auto ta1_z_xxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 323);

    auto ta1_z_xxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 325);

    auto ta1_z_xxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 327);

    auto ta1_z_xxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 328);

    auto ta1_z_xxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 329);

    auto ta1_z_xxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 330);

    auto ta1_z_xxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 331);

    auto ta1_z_xxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 332);

    auto ta1_z_xxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 333);

    auto ta1_z_xxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 334);

    auto ta1_z_xxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 335);

    auto ta1_z_xxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 336);

    auto ta1_z_xxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 337);

    auto ta1_z_xxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 338);

    auto ta1_z_xxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 339);

    auto ta1_z_xxyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 342);

    auto ta1_z_xxyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 345);

    auto ta1_z_xxyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 347);

    auto ta1_z_xxyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 348);

    auto ta1_z_xxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 350);

    auto ta1_z_xxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 351);

    auto ta1_z_xxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 352);

    auto ta1_z_xxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 353);

    auto ta1_z_xxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 354);

    auto ta1_z_xxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 355);

    auto ta1_z_xxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 356);

    auto ta1_z_xxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 357);

    auto ta1_z_xxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 358);

    auto ta1_z_xxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 359);

    auto ta1_z_xyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 361);

    auto ta1_z_xyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 363);

    auto ta1_z_xyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 364);

    auto ta1_z_xyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 366);

    auto ta1_z_xyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 367);

    auto ta1_z_xyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 368);

    auto ta1_z_xyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 369);

    auto ta1_z_xyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 377);

    auto ta1_z_xyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 378);

    auto ta1_z_xyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 379);

    auto ta1_z_xyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 386);

    auto ta1_z_xyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 387);

    auto ta1_z_xyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 388);

    auto ta1_z_xzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 392);

    auto ta1_z_xzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 394);

    auto ta1_z_xzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 395);

    auto ta1_z_xzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 396);

    auto ta1_z_xzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 397);

    auto ta1_z_xzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 398);

    auto ta1_z_xzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 399);

    auto ta1_z_yyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 400);

    auto ta1_z_yyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 401);

    auto ta1_z_yyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 402);

    auto ta1_z_yyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 403);

    auto ta1_z_yyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 404);

    auto ta1_z_yyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 405);

    auto ta1_z_yyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 406);

    auto ta1_z_yyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 407);

    auto ta1_z_yyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 408);

    auto ta1_z_yyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 409);

    auto ta1_z_yyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 411);

    auto ta1_z_yyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 412);

    auto ta1_z_yyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 413);

    auto ta1_z_yyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 415);

    auto ta1_z_yyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 416);

    auto ta1_z_yyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 417);

    auto ta1_z_yyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 418);

    auto ta1_z_yyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 419);

    auto ta1_z_yyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 420);

    auto ta1_z_yyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 421);

    auto ta1_z_yyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 422);

    auto ta1_z_yyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 423);

    auto ta1_z_yyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 424);

    auto ta1_z_yyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 425);

    auto ta1_z_yyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 426);

    auto ta1_z_yyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 427);

    auto ta1_z_yyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 428);

    auto ta1_z_yyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 429);

    auto ta1_z_yzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 430);

    auto ta1_z_yzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 432);

    auto ta1_z_yzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 434);

    auto ta1_z_yzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 435);

    auto ta1_z_yzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 436);

    auto ta1_z_yzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 437);

    auto ta1_z_yzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 438);

    auto ta1_z_yzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 439);

    auto ta1_z_zzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_gf + 440);

    auto ta1_z_zzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 441);

    auto ta1_z_zzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 442);

    auto ta1_z_zzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 443);

    auto ta1_z_zzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 444);

    auto ta1_z_zzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 445);

    auto ta1_z_zzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_gf + 446);

    auto ta1_z_zzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 447);

    auto ta1_z_zzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 448);

    auto ta1_z_zzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_gf + 449);

    // Set up components of auxiliary buffer : HD

    auto ta1_x_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd);

    auto ta1_x_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 1);

    auto ta1_x_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 2);

    auto ta1_x_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 3);

    auto ta1_x_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 4);

    auto ta1_x_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 5);

    auto ta1_x_xxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 6);

    auto ta1_x_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 7);

    auto ta1_x_xxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 8);

    auto ta1_x_xxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 12);

    auto ta1_x_xxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 13);

    auto ta1_x_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 14);

    auto ta1_x_xxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 16);

    auto ta1_x_xxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 17);

    auto ta1_x_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 18);

    auto ta1_x_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 19);

    auto ta1_x_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 20);

    auto ta1_x_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 21);

    auto ta1_x_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 22);

    auto ta1_x_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 30);

    auto ta1_x_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 31);

    auto ta1_x_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 32);

    auto ta1_x_xxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 33);

    auto ta1_x_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 34);

    auto ta1_x_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 35);

    auto ta1_x_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 36);

    auto ta1_x_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 37);

    auto ta1_x_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 38);

    auto ta1_x_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 39);

    auto ta1_x_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 40);

    auto ta1_x_xxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 50);

    auto ta1_x_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 54);

    auto ta1_x_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 55);

    auto ta1_x_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 56);

    auto ta1_x_xxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 57);

    auto ta1_x_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 58);

    auto ta1_x_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 59);

    auto ta1_x_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 61);

    auto ta1_x_xzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 84);

    auto ta1_x_xzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 85);

    auto ta1_x_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 86);

    auto ta1_x_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 90);

    auto ta1_x_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 91);

    auto ta1_x_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 92);

    auto ta1_x_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 93);

    auto ta1_x_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 94);

    auto ta1_x_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 95);

    auto ta1_x_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 104);

    auto ta1_x_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 106);

    auto ta1_x_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 107);

    auto ta1_x_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 110);

    auto ta1_x_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 112);

    auto ta1_x_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 113);

    auto ta1_x_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 116);

    auto ta1_x_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 118);

    auto ta1_x_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 119);

    auto ta1_x_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 120);

    auto ta1_x_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 121);

    auto ta1_x_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 122);

    auto ta1_x_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 123);

    auto ta1_x_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 124);

    auto ta1_x_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 125);

    auto ta1_y_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 126);

    auto ta1_y_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 127);

    auto ta1_y_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 128);

    auto ta1_y_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 129);

    auto ta1_y_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 130);

    auto ta1_y_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 131);

    auto ta1_y_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 133);

    auto ta1_y_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 144);

    auto ta1_y_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 145);

    auto ta1_y_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 146);

    auto ta1_y_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 147);

    auto ta1_y_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 148);

    auto ta1_y_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 158);

    auto ta1_y_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 160);

    auto ta1_y_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 161);

    auto ta1_y_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 162);

    auto ta1_y_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 163);

    auto ta1_y_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 164);

    auto ta1_y_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 165);

    auto ta1_y_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 166);

    auto ta1_y_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 182);

    auto ta1_y_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 184);

    auto ta1_y_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 185);

    auto ta1_y_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 187);

    auto ta1_y_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 189);

    auto ta1_y_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 190);

    auto ta1_y_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 202);

    auto ta1_y_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 212);

    auto ta1_y_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 214);

    auto ta1_y_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 215);

    auto ta1_y_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 216);

    auto ta1_y_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 217);

    auto ta1_y_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 218);

    auto ta1_y_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 219);

    auto ta1_y_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 220);

    auto ta1_y_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 221);

    auto ta1_y_yyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 223);

    auto ta1_y_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 224);

    auto ta1_y_yyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 225);

    auto ta1_y_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 226);

    auto ta1_y_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 227);

    auto ta1_y_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 228);

    auto ta1_y_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 229);

    auto ta1_y_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 230);

    auto ta1_y_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 231);

    auto ta1_y_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 232);

    auto ta1_y_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 233);

    auto ta1_y_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 234);

    auto ta1_y_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 235);

    auto ta1_y_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 236);

    auto ta1_y_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 237);

    auto ta1_y_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 238);

    auto ta1_y_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 239);

    auto ta1_y_yzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 241);

    auto ta1_y_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 243);

    auto ta1_y_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 244);

    auto ta1_y_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 246);

    auto ta1_y_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 247);

    auto ta1_y_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 248);

    auto ta1_y_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 249);

    auto ta1_y_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 250);

    auto ta1_y_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 251);

    auto ta1_z_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 252);

    auto ta1_z_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 253);

    auto ta1_z_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 254);

    auto ta1_z_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 255);

    auto ta1_z_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 256);

    auto ta1_z_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 257);

    auto ta1_z_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 266);

    auto ta1_z_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 271);

    auto ta1_z_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 273);

    auto ta1_z_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 274);

    auto ta1_z_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 282);

    auto ta1_z_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 283);

    auto ta1_z_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 284);

    auto ta1_z_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 286);

    auto ta1_z_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 287);

    auto ta1_z_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 289);

    auto ta1_z_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 291);

    auto ta1_z_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 292);

    auto ta1_z_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 306);

    auto ta1_z_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 307);

    auto ta1_z_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 308);

    auto ta1_z_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 310);

    auto ta1_z_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 311);

    auto ta1_z_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 313);

    auto ta1_z_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 315);

    auto ta1_z_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 316);

    auto ta1_z_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 328);

    auto ta1_z_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 338);

    auto ta1_z_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 340);

    auto ta1_z_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 341);

    auto ta1_z_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 342);

    auto ta1_z_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 343);

    auto ta1_z_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 344);

    auto ta1_z_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 345);

    auto ta1_z_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 346);

    auto ta1_z_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 347);

    auto ta1_z_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 350);

    auto ta1_z_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 352);

    auto ta1_z_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 353);

    auto ta1_z_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 354);

    auto ta1_z_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 355);

    auto ta1_z_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 356);

    auto ta1_z_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 357);

    auto ta1_z_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 358);

    auto ta1_z_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 359);

    auto ta1_z_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 360);

    auto ta1_z_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 361);

    auto ta1_z_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 362);

    auto ta1_z_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 363);

    auto ta1_z_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 364);

    auto ta1_z_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 365);

    auto ta1_z_yzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 367);

    auto ta1_z_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 368);

    auto ta1_z_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 369);

    auto ta1_z_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 370);

    auto ta1_z_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 371);

    auto ta1_z_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 372);

    auto ta1_z_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 373);

    auto ta1_z_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 374);

    auto ta1_z_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 375);

    auto ta1_z_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 376);

    auto ta1_z_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 377);

    // Set up components of auxiliary buffer : HD

    auto ta1_x_xxxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd);

    auto ta1_x_xxxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 1);

    auto ta1_x_xxxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 2);

    auto ta1_x_xxxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 3);

    auto ta1_x_xxxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 4);

    auto ta1_x_xxxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 5);

    auto ta1_x_xxxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 6);

    auto ta1_x_xxxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 7);

    auto ta1_x_xxxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 8);

    auto ta1_x_xxxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 12);

    auto ta1_x_xxxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 13);

    auto ta1_x_xxxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 14);

    auto ta1_x_xxxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 16);

    auto ta1_x_xxxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 17);

    auto ta1_x_xxxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 18);

    auto ta1_x_xxxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 19);

    auto ta1_x_xxxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 20);

    auto ta1_x_xxxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 21);

    auto ta1_x_xxxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 22);

    auto ta1_x_xxxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 30);

    auto ta1_x_xxxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 31);

    auto ta1_x_xxxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 32);

    auto ta1_x_xxxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 33);

    auto ta1_x_xxxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 34);

    auto ta1_x_xxxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 35);

    auto ta1_x_xxyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 36);

    auto ta1_x_xxyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 37);

    auto ta1_x_xxyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 38);

    auto ta1_x_xxyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 39);

    auto ta1_x_xxyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 40);

    auto ta1_x_xxyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 50);

    auto ta1_x_xxzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 54);

    auto ta1_x_xxzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 55);

    auto ta1_x_xxzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 56);

    auto ta1_x_xxzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 57);

    auto ta1_x_xxzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 58);

    auto ta1_x_xxzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 59);

    auto ta1_x_xyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 61);

    auto ta1_x_xzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 84);

    auto ta1_x_xzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 85);

    auto ta1_x_xzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 86);

    auto ta1_x_yyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 90);

    auto ta1_x_yyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 91);

    auto ta1_x_yyyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 92);

    auto ta1_x_yyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 93);

    auto ta1_x_yyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 94);

    auto ta1_x_yyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 95);

    auto ta1_x_yyyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 104);

    auto ta1_x_yyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 106);

    auto ta1_x_yyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 107);

    auto ta1_x_yyzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 110);

    auto ta1_x_yyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 112);

    auto ta1_x_yyzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 113);

    auto ta1_x_yzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 116);

    auto ta1_x_yzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 118);

    auto ta1_x_yzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 119);

    auto ta1_x_zzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 120);

    auto ta1_x_zzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 121);

    auto ta1_x_zzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 122);

    auto ta1_x_zzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 123);

    auto ta1_x_zzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 124);

    auto ta1_x_zzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 125);

    auto ta1_y_xxxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 126);

    auto ta1_y_xxxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 127);

    auto ta1_y_xxxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 128);

    auto ta1_y_xxxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 129);

    auto ta1_y_xxxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 130);

    auto ta1_y_xxxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 131);

    auto ta1_y_xxxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 133);

    auto ta1_y_xxxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 144);

    auto ta1_y_xxxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 145);

    auto ta1_y_xxxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 146);

    auto ta1_y_xxxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 147);

    auto ta1_y_xxxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 148);

    auto ta1_y_xxxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 158);

    auto ta1_y_xxxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 160);

    auto ta1_y_xxxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 161);

    auto ta1_y_xxyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 162);

    auto ta1_y_xxyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 163);

    auto ta1_y_xxyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 164);

    auto ta1_y_xxyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 165);

    auto ta1_y_xxyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 166);

    auto ta1_y_xxzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 182);

    auto ta1_y_xxzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 184);

    auto ta1_y_xxzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 185);

    auto ta1_y_xyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 187);

    auto ta1_y_xyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 189);

    auto ta1_y_xyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 190);

    auto ta1_y_xyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 202);

    auto ta1_y_xzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 212);

    auto ta1_y_xzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 214);

    auto ta1_y_xzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 215);

    auto ta1_y_yyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 216);

    auto ta1_y_yyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 217);

    auto ta1_y_yyyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 218);

    auto ta1_y_yyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 219);

    auto ta1_y_yyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 220);

    auto ta1_y_yyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 221);

    auto ta1_y_yyyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 223);

    auto ta1_y_yyyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 224);

    auto ta1_y_yyyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 225);

    auto ta1_y_yyyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 226);

    auto ta1_y_yyyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 227);

    auto ta1_y_yyyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 228);

    auto ta1_y_yyyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 229);

    auto ta1_y_yyyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 230);

    auto ta1_y_yyyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 231);

    auto ta1_y_yyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 232);

    auto ta1_y_yyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 233);

    auto ta1_y_yyzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 234);

    auto ta1_y_yyzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 235);

    auto ta1_y_yyzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 236);

    auto ta1_y_yyzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 237);

    auto ta1_y_yyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 238);

    auto ta1_y_yyzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 239);

    auto ta1_y_yzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 241);

    auto ta1_y_yzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 243);

    auto ta1_y_yzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 244);

    auto ta1_y_zzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 246);

    auto ta1_y_zzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 247);

    auto ta1_y_zzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 248);

    auto ta1_y_zzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 249);

    auto ta1_y_zzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 250);

    auto ta1_y_zzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 251);

    auto ta1_z_xxxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 252);

    auto ta1_z_xxxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 253);

    auto ta1_z_xxxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 254);

    auto ta1_z_xxxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 255);

    auto ta1_z_xxxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 256);

    auto ta1_z_xxxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 257);

    auto ta1_z_xxxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 266);

    auto ta1_z_xxxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 271);

    auto ta1_z_xxxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 273);

    auto ta1_z_xxxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 274);

    auto ta1_z_xxxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 282);

    auto ta1_z_xxxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 283);

    auto ta1_z_xxxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 284);

    auto ta1_z_xxxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 286);

    auto ta1_z_xxxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 287);

    auto ta1_z_xxyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 289);

    auto ta1_z_xxyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 291);

    auto ta1_z_xxyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 292);

    auto ta1_z_xxzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 306);

    auto ta1_z_xxzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 307);

    auto ta1_z_xxzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 308);

    auto ta1_z_xxzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 310);

    auto ta1_z_xxzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 311);

    auto ta1_z_xyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 313);

    auto ta1_z_xyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 315);

    auto ta1_z_xyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 316);

    auto ta1_z_xyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 328);

    auto ta1_z_xzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 338);

    auto ta1_z_xzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 340);

    auto ta1_z_xzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 341);

    auto ta1_z_yyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 342);

    auto ta1_z_yyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 343);

    auto ta1_z_yyyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 344);

    auto ta1_z_yyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 345);

    auto ta1_z_yyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 346);

    auto ta1_z_yyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 347);

    auto ta1_z_yyyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 350);

    auto ta1_z_yyyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 352);

    auto ta1_z_yyyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 353);

    auto ta1_z_yyyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 354);

    auto ta1_z_yyyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 355);

    auto ta1_z_yyyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 356);

    auto ta1_z_yyyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 357);

    auto ta1_z_yyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 358);

    auto ta1_z_yyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 359);

    auto ta1_z_yyzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 360);

    auto ta1_z_yyzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 361);

    auto ta1_z_yyzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 362);

    auto ta1_z_yyzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 363);

    auto ta1_z_yyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 364);

    auto ta1_z_yyzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 365);

    auto ta1_z_yzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 367);

    auto ta1_z_yzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 368);

    auto ta1_z_yzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 369);

    auto ta1_z_yzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 370);

    auto ta1_z_yzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 371);

    auto ta1_z_zzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 372);

    auto ta1_z_zzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 373);

    auto ta1_z_zzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 374);

    auto ta1_z_zzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 375);

    auto ta1_z_zzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 376);

    auto ta1_z_zzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 377);

    // Set up components of auxiliary buffer : HF

    auto ta_xxxxx_xxx_1 = pbuffer.data(idx_npot_1_hf);

    auto ta_xxxxx_xxy_1 = pbuffer.data(idx_npot_1_hf + 1);

    auto ta_xxxxx_xxz_1 = pbuffer.data(idx_npot_1_hf + 2);

    auto ta_xxxxx_xyy_1 = pbuffer.data(idx_npot_1_hf + 3);

    auto ta_xxxxx_xyz_1 = pbuffer.data(idx_npot_1_hf + 4);

    auto ta_xxxxx_xzz_1 = pbuffer.data(idx_npot_1_hf + 5);

    auto ta_xxxxx_yyy_1 = pbuffer.data(idx_npot_1_hf + 6);

    auto ta_xxxxx_yyz_1 = pbuffer.data(idx_npot_1_hf + 7);

    auto ta_xxxxx_yzz_1 = pbuffer.data(idx_npot_1_hf + 8);

    auto ta_xxxxx_zzz_1 = pbuffer.data(idx_npot_1_hf + 9);

    auto ta_xxxxy_xxx_1 = pbuffer.data(idx_npot_1_hf + 10);

    auto ta_xxxxy_xxy_1 = pbuffer.data(idx_npot_1_hf + 11);

    auto ta_xxxxy_xxz_1 = pbuffer.data(idx_npot_1_hf + 12);

    auto ta_xxxxy_xyy_1 = pbuffer.data(idx_npot_1_hf + 13);

    auto ta_xxxxy_xzz_1 = pbuffer.data(idx_npot_1_hf + 15);

    auto ta_xxxxy_yyy_1 = pbuffer.data(idx_npot_1_hf + 16);

    auto ta_xxxxz_xxx_1 = pbuffer.data(idx_npot_1_hf + 20);

    auto ta_xxxxz_xxy_1 = pbuffer.data(idx_npot_1_hf + 21);

    auto ta_xxxxz_xxz_1 = pbuffer.data(idx_npot_1_hf + 22);

    auto ta_xxxxz_xyy_1 = pbuffer.data(idx_npot_1_hf + 23);

    auto ta_xxxxz_xzz_1 = pbuffer.data(idx_npot_1_hf + 25);

    auto ta_xxxxz_zzz_1 = pbuffer.data(idx_npot_1_hf + 29);

    auto ta_xxxyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 30);

    auto ta_xxxyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 31);

    auto ta_xxxyy_xxz_1 = pbuffer.data(idx_npot_1_hf + 32);

    auto ta_xxxyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 33);

    auto ta_xxxyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 34);

    auto ta_xxxyy_xzz_1 = pbuffer.data(idx_npot_1_hf + 35);

    auto ta_xxxyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 36);

    auto ta_xxxyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 37);

    auto ta_xxxyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 38);

    auto ta_xxxzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 50);

    auto ta_xxxzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 51);

    auto ta_xxxzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 52);

    auto ta_xxxzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 53);

    auto ta_xxxzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 54);

    auto ta_xxxzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 55);

    auto ta_xxxzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 57);

    auto ta_xxxzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 58);

    auto ta_xxxzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 59);

    auto ta_xxyyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 60);

    auto ta_xxyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 61);

    auto ta_xxyyy_xxz_1 = pbuffer.data(idx_npot_1_hf + 62);

    auto ta_xxyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 63);

    auto ta_xxyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 64);

    auto ta_xxyyy_xzz_1 = pbuffer.data(idx_npot_1_hf + 65);

    auto ta_xxyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 66);

    auto ta_xxyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 67);

    auto ta_xxyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 68);

    auto ta_xxyyz_xxy_1 = pbuffer.data(idx_npot_1_hf + 71);

    auto ta_xxyyz_xyy_1 = pbuffer.data(idx_npot_1_hf + 73);

    auto ta_xxyzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 82);

    auto ta_xxyzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 85);

    auto ta_xxzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 90);

    auto ta_xxzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 91);

    auto ta_xxzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 92);

    auto ta_xxzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 93);

    auto ta_xxzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 94);

    auto ta_xxzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 95);

    auto ta_xxzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 97);

    auto ta_xxzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 98);

    auto ta_xxzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 99);

    auto ta_xyyyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 100);

    auto ta_xyyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 101);

    auto ta_xyyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 103);

    auto ta_xyyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 106);

    auto ta_xyyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 107);

    auto ta_xyyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 108);

    auto ta_xyyzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 127);

    auto ta_xyyzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 128);

    auto ta_xzzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 140);

    auto ta_xzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 142);

    auto ta_xzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 145);

    auto ta_xzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 147);

    auto ta_xzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 148);

    auto ta_xzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 149);

    auto ta_yyyyy_xxx_1 = pbuffer.data(idx_npot_1_hf + 150);

    auto ta_yyyyy_xxy_1 = pbuffer.data(idx_npot_1_hf + 151);

    auto ta_yyyyy_xxz_1 = pbuffer.data(idx_npot_1_hf + 152);

    auto ta_yyyyy_xyy_1 = pbuffer.data(idx_npot_1_hf + 153);

    auto ta_yyyyy_xyz_1 = pbuffer.data(idx_npot_1_hf + 154);

    auto ta_yyyyy_xzz_1 = pbuffer.data(idx_npot_1_hf + 155);

    auto ta_yyyyy_yyy_1 = pbuffer.data(idx_npot_1_hf + 156);

    auto ta_yyyyy_yyz_1 = pbuffer.data(idx_npot_1_hf + 157);

    auto ta_yyyyy_yzz_1 = pbuffer.data(idx_npot_1_hf + 158);

    auto ta_yyyyy_zzz_1 = pbuffer.data(idx_npot_1_hf + 159);

    auto ta_yyyyz_xxy_1 = pbuffer.data(idx_npot_1_hf + 161);

    auto ta_yyyyz_xyy_1 = pbuffer.data(idx_npot_1_hf + 163);

    auto ta_yyyyz_yyy_1 = pbuffer.data(idx_npot_1_hf + 166);

    auto ta_yyyyz_yyz_1 = pbuffer.data(idx_npot_1_hf + 167);

    auto ta_yyyyz_yzz_1 = pbuffer.data(idx_npot_1_hf + 168);

    auto ta_yyyyz_zzz_1 = pbuffer.data(idx_npot_1_hf + 169);

    auto ta_yyyzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 171);

    auto ta_yyyzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 172);

    auto ta_yyyzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 173);

    auto ta_yyyzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 174);

    auto ta_yyyzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 175);

    auto ta_yyyzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 176);

    auto ta_yyyzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 177);

    auto ta_yyyzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 178);

    auto ta_yyyzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 179);

    auto ta_yyzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 181);

    auto ta_yyzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 182);

    auto ta_yyzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 183);

    auto ta_yyzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 184);

    auto ta_yyzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 185);

    auto ta_yyzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 186);

    auto ta_yyzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 187);

    auto ta_yyzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 188);

    auto ta_yyzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 189);

    auto ta_yzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 192);

    auto ta_yzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 195);

    auto ta_yzzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 196);

    auto ta_yzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 197);

    auto ta_yzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 198);

    auto ta_yzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 199);

    auto ta_zzzzz_xxx_1 = pbuffer.data(idx_npot_1_hf + 200);

    auto ta_zzzzz_xxy_1 = pbuffer.data(idx_npot_1_hf + 201);

    auto ta_zzzzz_xxz_1 = pbuffer.data(idx_npot_1_hf + 202);

    auto ta_zzzzz_xyy_1 = pbuffer.data(idx_npot_1_hf + 203);

    auto ta_zzzzz_xyz_1 = pbuffer.data(idx_npot_1_hf + 204);

    auto ta_zzzzz_xzz_1 = pbuffer.data(idx_npot_1_hf + 205);

    auto ta_zzzzz_yyy_1 = pbuffer.data(idx_npot_1_hf + 206);

    auto ta_zzzzz_yyz_1 = pbuffer.data(idx_npot_1_hf + 207);

    auto ta_zzzzz_yzz_1 = pbuffer.data(idx_npot_1_hf + 208);

    auto ta_zzzzz_zzz_1 = pbuffer.data(idx_npot_1_hf + 209);

    // Set up components of auxiliary buffer : HF

    auto ta1_x_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf);

    auto ta1_x_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 1);

    auto ta1_x_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 2);

    auto ta1_x_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 3);

    auto ta1_x_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 4);

    auto ta1_x_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 5);

    auto ta1_x_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 6);

    auto ta1_x_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 7);

    auto ta1_x_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 8);

    auto ta1_x_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 9);

    auto ta1_x_xxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 10);

    auto ta1_x_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 11);

    auto ta1_x_xxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 12);

    auto ta1_x_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 13);

    auto ta1_x_xxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 14);

    auto ta1_x_xxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 15);

    auto ta1_x_xxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 16);

    auto ta1_x_xxxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 19);

    auto ta1_x_xxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 20);

    auto ta1_x_xxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 21);

    auto ta1_x_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 22);

    auto ta1_x_xxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 23);

    auto ta1_x_xxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 24);

    auto ta1_x_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 25);

    auto ta1_x_xxxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 26);

    auto ta1_x_xxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 27);

    auto ta1_x_xxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 28);

    auto ta1_x_xxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 29);

    auto ta1_x_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 30);

    auto ta1_x_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 31);

    auto ta1_x_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 32);

    auto ta1_x_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 33);

    auto ta1_x_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 34);

    auto ta1_x_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 35);

    auto ta1_x_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 36);

    auto ta1_x_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 37);

    auto ta1_x_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 38);

    auto ta1_x_xxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 39);

    auto ta1_x_xxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 42);

    auto ta1_x_xxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 45);

    auto ta1_x_xxxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 49);

    auto ta1_x_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 50);

    auto ta1_x_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 51);

    auto ta1_x_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 52);

    auto ta1_x_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 53);

    auto ta1_x_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 54);

    auto ta1_x_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 55);

    auto ta1_x_xxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 56);

    auto ta1_x_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 57);

    auto ta1_x_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 58);

    auto ta1_x_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 59);

    auto ta1_x_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 60);

    auto ta1_x_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 61);

    auto ta1_x_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 62);

    auto ta1_x_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 63);

    auto ta1_x_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 64);

    auto ta1_x_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 65);

    auto ta1_x_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 66);

    auto ta1_x_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 67);

    auto ta1_x_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 68);

    auto ta1_x_xxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 69);

    auto ta1_x_xxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 71);

    auto ta1_x_xxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 72);

    auto ta1_x_xxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 73);

    auto ta1_x_xxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 75);

    auto ta1_x_xxyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 76);

    auto ta1_x_xxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 79);

    auto ta1_x_xxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 80);

    auto ta1_x_xxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 82);

    auto ta1_x_xxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 84);

    auto ta1_x_xxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 85);

    auto ta1_x_xxyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 89);

    auto ta1_x_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 90);

    auto ta1_x_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 91);

    auto ta1_x_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 92);

    auto ta1_x_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 93);

    auto ta1_x_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 94);

    auto ta1_x_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 95);

    auto ta1_x_xxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 96);

    auto ta1_x_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 97);

    auto ta1_x_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 98);

    auto ta1_x_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 99);

    auto ta1_x_xyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 100);

    auto ta1_x_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 101);

    auto ta1_x_xyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 102);

    auto ta1_x_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 103);

    auto ta1_x_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 104);

    auto ta1_x_xyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 105);

    auto ta1_x_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 106);

    auto ta1_x_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 107);

    auto ta1_x_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 108);

    auto ta1_x_xyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 111);

    auto ta1_x_xyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 112);

    auto ta1_x_xyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 113);

    auto ta1_x_xyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 115);

    auto ta1_x_xyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 120);

    auto ta1_x_xyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 121);

    auto ta1_x_xyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 122);

    auto ta1_x_xyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 123);

    auto ta1_x_xyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 125);

    auto ta1_x_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 127);

    auto ta1_x_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 128);

    auto ta1_x_xyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 130);

    auto ta1_x_xyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 132);

    auto ta1_x_xyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 135);

    auto ta1_x_xzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 140);

    auto ta1_x_xzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 141);

    auto ta1_x_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 142);

    auto ta1_x_xzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 143);

    auto ta1_x_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 144);

    auto ta1_x_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 145);

    auto ta1_x_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 147);

    auto ta1_x_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 148);

    auto ta1_x_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 149);

    auto ta1_x_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 150);

    auto ta1_x_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 151);

    auto ta1_x_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 152);

    auto ta1_x_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 153);

    auto ta1_x_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 154);

    auto ta1_x_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 155);

    auto ta1_x_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 156);

    auto ta1_x_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 157);

    auto ta1_x_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 158);

    auto ta1_x_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 159);

    auto ta1_x_yyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 161);

    auto ta1_x_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 162);

    auto ta1_x_yyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 163);

    auto ta1_x_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 165);

    auto ta1_x_yyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 166);

    auto ta1_x_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 167);

    auto ta1_x_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 168);

    auto ta1_x_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 169);

    auto ta1_x_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 170);

    auto ta1_x_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 171);

    auto ta1_x_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 172);

    auto ta1_x_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 173);

    auto ta1_x_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 174);

    auto ta1_x_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 175);

    auto ta1_x_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 176);

    auto ta1_x_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 177);

    auto ta1_x_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 178);

    auto ta1_x_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 179);

    auto ta1_x_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 180);

    auto ta1_x_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 181);

    auto ta1_x_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 182);

    auto ta1_x_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 183);

    auto ta1_x_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 184);

    auto ta1_x_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 185);

    auto ta1_x_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 186);

    auto ta1_x_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 187);

    auto ta1_x_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 188);

    auto ta1_x_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 189);

    auto ta1_x_yzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 190);

    auto ta1_x_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 192);

    auto ta1_x_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 194);

    auto ta1_x_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 195);

    auto ta1_x_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 196);

    auto ta1_x_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 197);

    auto ta1_x_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 198);

    auto ta1_x_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 199);

    auto ta1_x_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 200);

    auto ta1_x_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 201);

    auto ta1_x_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 202);

    auto ta1_x_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 203);

    auto ta1_x_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 204);

    auto ta1_x_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 205);

    auto ta1_x_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 206);

    auto ta1_x_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 207);

    auto ta1_x_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 208);

    auto ta1_x_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 209);

    auto ta1_y_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 210);

    auto ta1_y_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 211);

    auto ta1_y_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 212);

    auto ta1_y_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 213);

    auto ta1_y_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 214);

    auto ta1_y_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 215);

    auto ta1_y_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 216);

    auto ta1_y_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 217);

    auto ta1_y_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 218);

    auto ta1_y_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 219);

    auto ta1_y_xxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 220);

    auto ta1_y_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 221);

    auto ta1_y_xxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 222);

    auto ta1_y_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 223);

    auto ta1_y_xxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 224);

    auto ta1_y_xxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 225);

    auto ta1_y_xxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 226);

    auto ta1_y_xxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 227);

    auto ta1_y_xxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 228);

    auto ta1_y_xxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 230);

    auto ta1_y_xxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 231);

    auto ta1_y_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 232);

    auto ta1_y_xxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 233);

    auto ta1_y_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 235);

    auto ta1_y_xxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 237);

    auto ta1_y_xxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 238);

    auto ta1_y_xxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 239);

    auto ta1_y_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 240);

    auto ta1_y_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 241);

    auto ta1_y_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 242);

    auto ta1_y_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 243);

    auto ta1_y_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 244);

    auto ta1_y_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 245);

    auto ta1_y_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 246);

    auto ta1_y_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 247);

    auto ta1_y_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 248);

    auto ta1_y_xxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 249);

    auto ta1_y_xxxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 251);

    auto ta1_y_xxxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 253);

    auto ta1_y_xxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 257);

    auto ta1_y_xxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 258);

    auto ta1_y_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 260);

    auto ta1_y_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 261);

    auto ta1_y_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 262);

    auto ta1_y_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 263);

    auto ta1_y_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 264);

    auto ta1_y_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 265);

    auto ta1_y_xxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 266);

    auto ta1_y_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 267);

    auto ta1_y_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 268);

    auto ta1_y_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 269);

    auto ta1_y_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 270);

    auto ta1_y_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 271);

    auto ta1_y_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 272);

    auto ta1_y_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 273);

    auto ta1_y_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 274);

    auto ta1_y_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 275);

    auto ta1_y_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 276);

    auto ta1_y_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 277);

    auto ta1_y_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 278);

    auto ta1_y_xxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 279);

    auto ta1_y_xxyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 280);

    auto ta1_y_xxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 281);

    auto ta1_y_xxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 283);

    auto ta1_y_xxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 287);

    auto ta1_y_xxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 288);

    auto ta1_y_xxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 289);

    auto ta1_y_xxyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 291);

    auto ta1_y_xxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 292);

    auto ta1_y_xxyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 293);

    auto ta1_y_xxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 295);

    auto ta1_y_xxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 296);

    auto ta1_y_xxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 297);

    auto ta1_y_xxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 298);

    auto ta1_y_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 300);

    auto ta1_y_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 301);

    auto ta1_y_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 302);

    auto ta1_y_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 303);

    auto ta1_y_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 304);

    auto ta1_y_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 305);

    auto ta1_y_xxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 306);

    auto ta1_y_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 307);

    auto ta1_y_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 308);

    auto ta1_y_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 309);

    auto ta1_y_xyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 310);

    auto ta1_y_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 311);

    auto ta1_y_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 313);

    auto ta1_y_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 314);

    auto ta1_y_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 316);

    auto ta1_y_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 317);

    auto ta1_y_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 318);

    auto ta1_y_xyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 319);

    auto ta1_y_xyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 327);

    auto ta1_y_xyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 328);

    auto ta1_y_xyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 329);

    auto ta1_y_xyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 334);

    auto ta1_y_xyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 336);

    auto ta1_y_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 337);

    auto ta1_y_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 338);

    auto ta1_y_xyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 339);

    auto ta1_y_xyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 346);

    auto ta1_y_xyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 347);

    auto ta1_y_xyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 348);

    auto ta1_y_xzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 350);

    auto ta1_y_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 352);

    auto ta1_y_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 354);

    auto ta1_y_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 355);

    auto ta1_y_xzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 356);

    auto ta1_y_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 357);

    auto ta1_y_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 358);

    auto ta1_y_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 359);

    auto ta1_y_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 360);

    auto ta1_y_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 361);

    auto ta1_y_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 362);

    auto ta1_y_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 363);

    auto ta1_y_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 364);

    auto ta1_y_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 365);

    auto ta1_y_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 366);

    auto ta1_y_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 367);

    auto ta1_y_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 368);

    auto ta1_y_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 369);

    auto ta1_y_yyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 370);

    auto ta1_y_yyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 371);

    auto ta1_y_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 372);

    auto ta1_y_yyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 373);

    auto ta1_y_yyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 374);

    auto ta1_y_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 375);

    auto ta1_y_yyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 376);

    auto ta1_y_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 377);

    auto ta1_y_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 378);

    auto ta1_y_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 379);

    auto ta1_y_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 380);

    auto ta1_y_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 381);

    auto ta1_y_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 382);

    auto ta1_y_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 383);

    auto ta1_y_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 384);

    auto ta1_y_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 385);

    auto ta1_y_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 386);

    auto ta1_y_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 387);

    auto ta1_y_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 388);

    auto ta1_y_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 389);

    auto ta1_y_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 390);

    auto ta1_y_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 391);

    auto ta1_y_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 392);

    auto ta1_y_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 393);

    auto ta1_y_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 394);

    auto ta1_y_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 395);

    auto ta1_y_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 396);

    auto ta1_y_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 397);

    auto ta1_y_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 398);

    auto ta1_y_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 399);

    auto ta1_y_yzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 401);

    auto ta1_y_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 402);

    auto ta1_y_yzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 403);

    auto ta1_y_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 404);

    auto ta1_y_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 405);

    auto ta1_y_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 406);

    auto ta1_y_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 407);

    auto ta1_y_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 408);

    auto ta1_y_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 409);

    auto ta1_y_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 410);

    auto ta1_y_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 411);

    auto ta1_y_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 412);

    auto ta1_y_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 413);

    auto ta1_y_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 414);

    auto ta1_y_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 415);

    auto ta1_y_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 416);

    auto ta1_y_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 417);

    auto ta1_y_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 418);

    auto ta1_y_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 419);

    auto ta1_z_xxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 420);

    auto ta1_z_xxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 421);

    auto ta1_z_xxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 422);

    auto ta1_z_xxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 423);

    auto ta1_z_xxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 424);

    auto ta1_z_xxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 425);

    auto ta1_z_xxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 426);

    auto ta1_z_xxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 427);

    auto ta1_z_xxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 428);

    auto ta1_z_xxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 429);

    auto ta1_z_xxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 430);

    auto ta1_z_xxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 431);

    auto ta1_z_xxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 432);

    auto ta1_z_xxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 433);

    auto ta1_z_xxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 435);

    auto ta1_z_xxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 436);

    auto ta1_z_xxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 437);

    auto ta1_z_xxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 438);

    auto ta1_z_xxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 440);

    auto ta1_z_xxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 441);

    auto ta1_z_xxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 442);

    auto ta1_z_xxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 443);

    auto ta1_z_xxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 444);

    auto ta1_z_xxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 445);

    auto ta1_z_xxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 447);

    auto ta1_z_xxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 448);

    auto ta1_z_xxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 449);

    auto ta1_z_xxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 450);

    auto ta1_z_xxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 451);

    auto ta1_z_xxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 452);

    auto ta1_z_xxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 453);

    auto ta1_z_xxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 454);

    auto ta1_z_xxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 455);

    auto ta1_z_xxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 456);

    auto ta1_z_xxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 457);

    auto ta1_z_xxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 458);

    auto ta1_z_xxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 459);

    auto ta1_z_xxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 462);

    auto ta1_z_xxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 465);

    auto ta1_z_xxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 467);

    auto ta1_z_xxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 468);

    auto ta1_z_xxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 470);

    auto ta1_z_xxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 471);

    auto ta1_z_xxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 472);

    auto ta1_z_xxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 473);

    auto ta1_z_xxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 474);

    auto ta1_z_xxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 475);

    auto ta1_z_xxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 476);

    auto ta1_z_xxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 477);

    auto ta1_z_xxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 478);

    auto ta1_z_xxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 479);

    auto ta1_z_xxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 480);

    auto ta1_z_xxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 481);

    auto ta1_z_xxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 482);

    auto ta1_z_xxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 483);

    auto ta1_z_xxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 484);

    auto ta1_z_xxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 485);

    auto ta1_z_xxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 486);

    auto ta1_z_xxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 487);

    auto ta1_z_xxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 488);

    auto ta1_z_xxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 489);

    auto ta1_z_xxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 491);

    auto ta1_z_xxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 492);

    auto ta1_z_xxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 493);

    auto ta1_z_xxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 495);

    auto ta1_z_xxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 497);

    auto ta1_z_xxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 498);

    auto ta1_z_xxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 499);

    auto ta1_z_xxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 500);

    auto ta1_z_xxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 502);

    auto ta1_z_xxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 505);

    auto ta1_z_xxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 506);

    auto ta1_z_xxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 507);

    auto ta1_z_xxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 508);

    auto ta1_z_xxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 510);

    auto ta1_z_xxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 511);

    auto ta1_z_xxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 512);

    auto ta1_z_xxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 513);

    auto ta1_z_xxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 514);

    auto ta1_z_xxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 515);

    auto ta1_z_xxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 516);

    auto ta1_z_xxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 517);

    auto ta1_z_xxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 518);

    auto ta1_z_xxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 519);

    auto ta1_z_xyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 520);

    auto ta1_z_xyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 521);

    auto ta1_z_xyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 523);

    auto ta1_z_xyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 524);

    auto ta1_z_xyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 526);

    auto ta1_z_xyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 527);

    auto ta1_z_xyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 528);

    auto ta1_z_xyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 529);

    auto ta1_z_xyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 537);

    auto ta1_z_xyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 538);

    auto ta1_z_xyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 539);

    auto ta1_z_xyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 544);

    auto ta1_z_xyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 546);

    auto ta1_z_xyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 547);

    auto ta1_z_xyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 548);

    auto ta1_z_xyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 549);

    auto ta1_z_xyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 556);

    auto ta1_z_xyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 557);

    auto ta1_z_xyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 558);

    auto ta1_z_xzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 560);

    auto ta1_z_xzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 562);

    auto ta1_z_xzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 564);

    auto ta1_z_xzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 565);

    auto ta1_z_xzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 566);

    auto ta1_z_xzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 567);

    auto ta1_z_xzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 568);

    auto ta1_z_xzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 569);

    auto ta1_z_yyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 570);

    auto ta1_z_yyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 571);

    auto ta1_z_yyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 572);

    auto ta1_z_yyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 573);

    auto ta1_z_yyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 574);

    auto ta1_z_yyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 575);

    auto ta1_z_yyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 576);

    auto ta1_z_yyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 577);

    auto ta1_z_yyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 578);

    auto ta1_z_yyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 579);

    auto ta1_z_yyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 581);

    auto ta1_z_yyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 582);

    auto ta1_z_yyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 583);

    auto ta1_z_yyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 584);

    auto ta1_z_yyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 585);

    auto ta1_z_yyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 586);

    auto ta1_z_yyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 587);

    auto ta1_z_yyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 588);

    auto ta1_z_yyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 589);

    auto ta1_z_yyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 590);

    auto ta1_z_yyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 591);

    auto ta1_z_yyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 592);

    auto ta1_z_yyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 593);

    auto ta1_z_yyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 594);

    auto ta1_z_yyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 595);

    auto ta1_z_yyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 596);

    auto ta1_z_yyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 597);

    auto ta1_z_yyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 598);

    auto ta1_z_yyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 599);

    auto ta1_z_yyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 600);

    auto ta1_z_yyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 601);

    auto ta1_z_yyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 602);

    auto ta1_z_yyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 603);

    auto ta1_z_yyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 604);

    auto ta1_z_yyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 605);

    auto ta1_z_yyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 606);

    auto ta1_z_yyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 607);

    auto ta1_z_yyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 608);

    auto ta1_z_yyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 609);

    auto ta1_z_yzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 610);

    auto ta1_z_yzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 611);

    auto ta1_z_yzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 612);

    auto ta1_z_yzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 613);

    auto ta1_z_yzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 614);

    auto ta1_z_yzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 615);

    auto ta1_z_yzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 616);

    auto ta1_z_yzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 617);

    auto ta1_z_yzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 618);

    auto ta1_z_yzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 619);

    auto ta1_z_zzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_hf + 620);

    auto ta1_z_zzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 621);

    auto ta1_z_zzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 622);

    auto ta1_z_zzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 623);

    auto ta1_z_zzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 624);

    auto ta1_z_zzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 625);

    auto ta1_z_zzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_hf + 626);

    auto ta1_z_zzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 627);

    auto ta1_z_zzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 628);

    auto ta1_z_zzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_hf + 629);

    // Set up components of auxiliary buffer : HF

    auto ta1_x_xxxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf);

    auto ta1_x_xxxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 1);

    auto ta1_x_xxxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 2);

    auto ta1_x_xxxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 3);

    auto ta1_x_xxxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 4);

    auto ta1_x_xxxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 5);

    auto ta1_x_xxxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 6);

    auto ta1_x_xxxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 7);

    auto ta1_x_xxxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 8);

    auto ta1_x_xxxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 9);

    auto ta1_x_xxxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 10);

    auto ta1_x_xxxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 11);

    auto ta1_x_xxxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 12);

    auto ta1_x_xxxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 13);

    auto ta1_x_xxxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 14);

    auto ta1_x_xxxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 15);

    auto ta1_x_xxxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 16);

    auto ta1_x_xxxxy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 19);

    auto ta1_x_xxxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 20);

    auto ta1_x_xxxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 21);

    auto ta1_x_xxxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 22);

    auto ta1_x_xxxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 23);

    auto ta1_x_xxxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 24);

    auto ta1_x_xxxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 25);

    auto ta1_x_xxxxz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 26);

    auto ta1_x_xxxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 27);

    auto ta1_x_xxxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 28);

    auto ta1_x_xxxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 29);

    auto ta1_x_xxxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 30);

    auto ta1_x_xxxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 31);

    auto ta1_x_xxxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 32);

    auto ta1_x_xxxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 33);

    auto ta1_x_xxxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 34);

    auto ta1_x_xxxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 35);

    auto ta1_x_xxxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 36);

    auto ta1_x_xxxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 37);

    auto ta1_x_xxxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 38);

    auto ta1_x_xxxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 39);

    auto ta1_x_xxxyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 42);

    auto ta1_x_xxxyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 45);

    auto ta1_x_xxxyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 49);

    auto ta1_x_xxxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 50);

    auto ta1_x_xxxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 51);

    auto ta1_x_xxxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 52);

    auto ta1_x_xxxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 53);

    auto ta1_x_xxxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 54);

    auto ta1_x_xxxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 55);

    auto ta1_x_xxxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 56);

    auto ta1_x_xxxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 57);

    auto ta1_x_xxxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 58);

    auto ta1_x_xxxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 59);

    auto ta1_x_xxyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 60);

    auto ta1_x_xxyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 61);

    auto ta1_x_xxyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 62);

    auto ta1_x_xxyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 63);

    auto ta1_x_xxyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 64);

    auto ta1_x_xxyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 65);

    auto ta1_x_xxyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 66);

    auto ta1_x_xxyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 67);

    auto ta1_x_xxyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 68);

    auto ta1_x_xxyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 69);

    auto ta1_x_xxyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 71);

    auto ta1_x_xxyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 72);

    auto ta1_x_xxyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 73);

    auto ta1_x_xxyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 75);

    auto ta1_x_xxyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 76);

    auto ta1_x_xxyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 79);

    auto ta1_x_xxyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 80);

    auto ta1_x_xxyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 82);

    auto ta1_x_xxyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 84);

    auto ta1_x_xxyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 85);

    auto ta1_x_xxyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 89);

    auto ta1_x_xxzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 90);

    auto ta1_x_xxzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 91);

    auto ta1_x_xxzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 92);

    auto ta1_x_xxzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 93);

    auto ta1_x_xxzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 94);

    auto ta1_x_xxzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 95);

    auto ta1_x_xxzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 96);

    auto ta1_x_xxzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 97);

    auto ta1_x_xxzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 98);

    auto ta1_x_xxzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 99);

    auto ta1_x_xyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 100);

    auto ta1_x_xyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 101);

    auto ta1_x_xyyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 102);

    auto ta1_x_xyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 103);

    auto ta1_x_xyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 104);

    auto ta1_x_xyyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 105);

    auto ta1_x_xyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 106);

    auto ta1_x_xyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 107);

    auto ta1_x_xyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 108);

    auto ta1_x_xyyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 111);

    auto ta1_x_xyyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 112);

    auto ta1_x_xyyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 113);

    auto ta1_x_xyyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 115);

    auto ta1_x_xyyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 120);

    auto ta1_x_xyyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 121);

    auto ta1_x_xyyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 122);

    auto ta1_x_xyyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 123);

    auto ta1_x_xyyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 125);

    auto ta1_x_xyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 127);

    auto ta1_x_xyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 128);

    auto ta1_x_xyzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 130);

    auto ta1_x_xyzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 132);

    auto ta1_x_xyzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 135);

    auto ta1_x_xzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 140);

    auto ta1_x_xzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 141);

    auto ta1_x_xzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 142);

    auto ta1_x_xzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 143);

    auto ta1_x_xzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 144);

    auto ta1_x_xzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 145);

    auto ta1_x_xzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 147);

    auto ta1_x_xzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 148);

    auto ta1_x_xzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 149);

    auto ta1_x_yyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 150);

    auto ta1_x_yyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 151);

    auto ta1_x_yyyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 152);

    auto ta1_x_yyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 153);

    auto ta1_x_yyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 154);

    auto ta1_x_yyyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 155);

    auto ta1_x_yyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 156);

    auto ta1_x_yyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 157);

    auto ta1_x_yyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 158);

    auto ta1_x_yyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 159);

    auto ta1_x_yyyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 161);

    auto ta1_x_yyyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 162);

    auto ta1_x_yyyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 163);

    auto ta1_x_yyyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 165);

    auto ta1_x_yyyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 166);

    auto ta1_x_yyyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 167);

    auto ta1_x_yyyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 168);

    auto ta1_x_yyyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 169);

    auto ta1_x_yyyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 170);

    auto ta1_x_yyyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 171);

    auto ta1_x_yyyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 172);

    auto ta1_x_yyyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 173);

    auto ta1_x_yyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 174);

    auto ta1_x_yyyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 175);

    auto ta1_x_yyyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 176);

    auto ta1_x_yyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 177);

    auto ta1_x_yyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 178);

    auto ta1_x_yyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 179);

    auto ta1_x_yyzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 180);

    auto ta1_x_yyzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 181);

    auto ta1_x_yyzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 182);

    auto ta1_x_yyzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 183);

    auto ta1_x_yyzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 184);

    auto ta1_x_yyzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 185);

    auto ta1_x_yyzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 186);

    auto ta1_x_yyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 187);

    auto ta1_x_yyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 188);

    auto ta1_x_yyzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 189);

    auto ta1_x_yzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 190);

    auto ta1_x_yzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 192);

    auto ta1_x_yzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 194);

    auto ta1_x_yzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 195);

    auto ta1_x_yzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 196);

    auto ta1_x_yzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 197);

    auto ta1_x_yzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 198);

    auto ta1_x_yzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 199);

    auto ta1_x_zzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 200);

    auto ta1_x_zzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 201);

    auto ta1_x_zzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 202);

    auto ta1_x_zzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 203);

    auto ta1_x_zzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 204);

    auto ta1_x_zzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 205);

    auto ta1_x_zzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 206);

    auto ta1_x_zzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 207);

    auto ta1_x_zzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 208);

    auto ta1_x_zzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 209);

    auto ta1_y_xxxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 210);

    auto ta1_y_xxxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 211);

    auto ta1_y_xxxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 212);

    auto ta1_y_xxxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 213);

    auto ta1_y_xxxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 214);

    auto ta1_y_xxxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 215);

    auto ta1_y_xxxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 216);

    auto ta1_y_xxxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 217);

    auto ta1_y_xxxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 218);

    auto ta1_y_xxxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 219);

    auto ta1_y_xxxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 220);

    auto ta1_y_xxxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 221);

    auto ta1_y_xxxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 222);

    auto ta1_y_xxxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 223);

    auto ta1_y_xxxxy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 224);

    auto ta1_y_xxxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 225);

    auto ta1_y_xxxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 226);

    auto ta1_y_xxxxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 227);

    auto ta1_y_xxxxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 228);

    auto ta1_y_xxxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 230);

    auto ta1_y_xxxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 231);

    auto ta1_y_xxxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 232);

    auto ta1_y_xxxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 233);

    auto ta1_y_xxxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 235);

    auto ta1_y_xxxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 237);

    auto ta1_y_xxxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 238);

    auto ta1_y_xxxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 239);

    auto ta1_y_xxxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 240);

    auto ta1_y_xxxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 241);

    auto ta1_y_xxxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 242);

    auto ta1_y_xxxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 243);

    auto ta1_y_xxxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 244);

    auto ta1_y_xxxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 245);

    auto ta1_y_xxxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 246);

    auto ta1_y_xxxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 247);

    auto ta1_y_xxxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 248);

    auto ta1_y_xxxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 249);

    auto ta1_y_xxxyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 251);

    auto ta1_y_xxxyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 253);

    auto ta1_y_xxxyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 257);

    auto ta1_y_xxxyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 258);

    auto ta1_y_xxxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 260);

    auto ta1_y_xxxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 261);

    auto ta1_y_xxxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 262);

    auto ta1_y_xxxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 263);

    auto ta1_y_xxxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 264);

    auto ta1_y_xxxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 265);

    auto ta1_y_xxxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 266);

    auto ta1_y_xxxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 267);

    auto ta1_y_xxxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 268);

    auto ta1_y_xxxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 269);

    auto ta1_y_xxyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 270);

    auto ta1_y_xxyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 271);

    auto ta1_y_xxyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 272);

    auto ta1_y_xxyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 273);

    auto ta1_y_xxyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 274);

    auto ta1_y_xxyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 275);

    auto ta1_y_xxyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 276);

    auto ta1_y_xxyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 277);

    auto ta1_y_xxyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 278);

    auto ta1_y_xxyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 279);

    auto ta1_y_xxyyz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 280);

    auto ta1_y_xxyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 281);

    auto ta1_y_xxyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 283);

    auto ta1_y_xxyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 287);

    auto ta1_y_xxyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 288);

    auto ta1_y_xxyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 289);

    auto ta1_y_xxyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 291);

    auto ta1_y_xxyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 292);

    auto ta1_y_xxyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 293);

    auto ta1_y_xxyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 295);

    auto ta1_y_xxyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 296);

    auto ta1_y_xxyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 297);

    auto ta1_y_xxyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 298);

    auto ta1_y_xxzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 300);

    auto ta1_y_xxzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 301);

    auto ta1_y_xxzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 302);

    auto ta1_y_xxzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 303);

    auto ta1_y_xxzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 304);

    auto ta1_y_xxzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 305);

    auto ta1_y_xxzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 306);

    auto ta1_y_xxzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 307);

    auto ta1_y_xxzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 308);

    auto ta1_y_xxzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 309);

    auto ta1_y_xyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 310);

    auto ta1_y_xyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 311);

    auto ta1_y_xyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 313);

    auto ta1_y_xyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 314);

    auto ta1_y_xyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 316);

    auto ta1_y_xyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 317);

    auto ta1_y_xyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 318);

    auto ta1_y_xyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 319);

    auto ta1_y_xyyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 327);

    auto ta1_y_xyyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 328);

    auto ta1_y_xyyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 329);

    auto ta1_y_xyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 334);

    auto ta1_y_xyyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 336);

    auto ta1_y_xyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 337);

    auto ta1_y_xyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 338);

    auto ta1_y_xyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 339);

    auto ta1_y_xyzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 346);

    auto ta1_y_xyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 347);

    auto ta1_y_xyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 348);

    auto ta1_y_xzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 350);

    auto ta1_y_xzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 352);

    auto ta1_y_xzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 354);

    auto ta1_y_xzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 355);

    auto ta1_y_xzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 356);

    auto ta1_y_xzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 357);

    auto ta1_y_xzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 358);

    auto ta1_y_xzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 359);

    auto ta1_y_yyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 360);

    auto ta1_y_yyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 361);

    auto ta1_y_yyyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 362);

    auto ta1_y_yyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 363);

    auto ta1_y_yyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 364);

    auto ta1_y_yyyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 365);

    auto ta1_y_yyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 366);

    auto ta1_y_yyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 367);

    auto ta1_y_yyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 368);

    auto ta1_y_yyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 369);

    auto ta1_y_yyyyz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 370);

    auto ta1_y_yyyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 371);

    auto ta1_y_yyyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 372);

    auto ta1_y_yyyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 373);

    auto ta1_y_yyyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 374);

    auto ta1_y_yyyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 375);

    auto ta1_y_yyyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 376);

    auto ta1_y_yyyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 377);

    auto ta1_y_yyyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 378);

    auto ta1_y_yyyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 379);

    auto ta1_y_yyyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 380);

    auto ta1_y_yyyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 381);

    auto ta1_y_yyyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 382);

    auto ta1_y_yyyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 383);

    auto ta1_y_yyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 384);

    auto ta1_y_yyyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 385);

    auto ta1_y_yyyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 386);

    auto ta1_y_yyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 387);

    auto ta1_y_yyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 388);

    auto ta1_y_yyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 389);

    auto ta1_y_yyzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 390);

    auto ta1_y_yyzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 391);

    auto ta1_y_yyzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 392);

    auto ta1_y_yyzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 393);

    auto ta1_y_yyzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 394);

    auto ta1_y_yyzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 395);

    auto ta1_y_yyzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 396);

    auto ta1_y_yyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 397);

    auto ta1_y_yyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 398);

    auto ta1_y_yyzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 399);

    auto ta1_y_yzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 401);

    auto ta1_y_yzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 402);

    auto ta1_y_yzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 403);

    auto ta1_y_yzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 404);

    auto ta1_y_yzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 405);

    auto ta1_y_yzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 406);

    auto ta1_y_yzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 407);

    auto ta1_y_yzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 408);

    auto ta1_y_yzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 409);

    auto ta1_y_zzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 410);

    auto ta1_y_zzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 411);

    auto ta1_y_zzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 412);

    auto ta1_y_zzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 413);

    auto ta1_y_zzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 414);

    auto ta1_y_zzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 415);

    auto ta1_y_zzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 416);

    auto ta1_y_zzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 417);

    auto ta1_y_zzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 418);

    auto ta1_y_zzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 419);

    auto ta1_z_xxxxx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 420);

    auto ta1_z_xxxxx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 421);

    auto ta1_z_xxxxx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 422);

    auto ta1_z_xxxxx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 423);

    auto ta1_z_xxxxx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 424);

    auto ta1_z_xxxxx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 425);

    auto ta1_z_xxxxx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 426);

    auto ta1_z_xxxxx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 427);

    auto ta1_z_xxxxx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 428);

    auto ta1_z_xxxxx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 429);

    auto ta1_z_xxxxy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 430);

    auto ta1_z_xxxxy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 431);

    auto ta1_z_xxxxy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 432);

    auto ta1_z_xxxxy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 433);

    auto ta1_z_xxxxy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 435);

    auto ta1_z_xxxxy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 436);

    auto ta1_z_xxxxy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 437);

    auto ta1_z_xxxxy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 438);

    auto ta1_z_xxxxz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 440);

    auto ta1_z_xxxxz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 441);

    auto ta1_z_xxxxz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 442);

    auto ta1_z_xxxxz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 443);

    auto ta1_z_xxxxz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 444);

    auto ta1_z_xxxxz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 445);

    auto ta1_z_xxxxz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 447);

    auto ta1_z_xxxxz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 448);

    auto ta1_z_xxxxz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 449);

    auto ta1_z_xxxyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 450);

    auto ta1_z_xxxyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 451);

    auto ta1_z_xxxyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 452);

    auto ta1_z_xxxyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 453);

    auto ta1_z_xxxyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 454);

    auto ta1_z_xxxyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 455);

    auto ta1_z_xxxyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 456);

    auto ta1_z_xxxyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 457);

    auto ta1_z_xxxyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 458);

    auto ta1_z_xxxyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 459);

    auto ta1_z_xxxyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 462);

    auto ta1_z_xxxyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 465);

    auto ta1_z_xxxyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 467);

    auto ta1_z_xxxyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 468);

    auto ta1_z_xxxzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 470);

    auto ta1_z_xxxzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 471);

    auto ta1_z_xxxzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 472);

    auto ta1_z_xxxzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 473);

    auto ta1_z_xxxzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 474);

    auto ta1_z_xxxzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 475);

    auto ta1_z_xxxzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 476);

    auto ta1_z_xxxzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 477);

    auto ta1_z_xxxzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 478);

    auto ta1_z_xxxzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 479);

    auto ta1_z_xxyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 480);

    auto ta1_z_xxyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 481);

    auto ta1_z_xxyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 482);

    auto ta1_z_xxyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 483);

    auto ta1_z_xxyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 484);

    auto ta1_z_xxyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 485);

    auto ta1_z_xxyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 486);

    auto ta1_z_xxyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 487);

    auto ta1_z_xxyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 488);

    auto ta1_z_xxyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 489);

    auto ta1_z_xxyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 491);

    auto ta1_z_xxyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 492);

    auto ta1_z_xxyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 493);

    auto ta1_z_xxyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 495);

    auto ta1_z_xxyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 497);

    auto ta1_z_xxyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 498);

    auto ta1_z_xxyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 499);

    auto ta1_z_xxyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 500);

    auto ta1_z_xxyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 502);

    auto ta1_z_xxyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 505);

    auto ta1_z_xxyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 506);

    auto ta1_z_xxyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 507);

    auto ta1_z_xxyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 508);

    auto ta1_z_xxzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 510);

    auto ta1_z_xxzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 511);

    auto ta1_z_xxzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 512);

    auto ta1_z_xxzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 513);

    auto ta1_z_xxzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 514);

    auto ta1_z_xxzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 515);

    auto ta1_z_xxzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 516);

    auto ta1_z_xxzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 517);

    auto ta1_z_xxzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 518);

    auto ta1_z_xxzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 519);

    auto ta1_z_xyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 520);

    auto ta1_z_xyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 521);

    auto ta1_z_xyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 523);

    auto ta1_z_xyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 524);

    auto ta1_z_xyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 526);

    auto ta1_z_xyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 527);

    auto ta1_z_xyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 528);

    auto ta1_z_xyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 529);

    auto ta1_z_xyyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 537);

    auto ta1_z_xyyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 538);

    auto ta1_z_xyyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 539);

    auto ta1_z_xyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 544);

    auto ta1_z_xyyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 546);

    auto ta1_z_xyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 547);

    auto ta1_z_xyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 548);

    auto ta1_z_xyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 549);

    auto ta1_z_xyzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 556);

    auto ta1_z_xyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 557);

    auto ta1_z_xyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 558);

    auto ta1_z_xzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 560);

    auto ta1_z_xzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 562);

    auto ta1_z_xzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 564);

    auto ta1_z_xzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 565);

    auto ta1_z_xzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 566);

    auto ta1_z_xzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 567);

    auto ta1_z_xzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 568);

    auto ta1_z_xzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 569);

    auto ta1_z_yyyyy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 570);

    auto ta1_z_yyyyy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 571);

    auto ta1_z_yyyyy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 572);

    auto ta1_z_yyyyy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 573);

    auto ta1_z_yyyyy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 574);

    auto ta1_z_yyyyy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 575);

    auto ta1_z_yyyyy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 576);

    auto ta1_z_yyyyy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 577);

    auto ta1_z_yyyyy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 578);

    auto ta1_z_yyyyy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 579);

    auto ta1_z_yyyyz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 581);

    auto ta1_z_yyyyz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 582);

    auto ta1_z_yyyyz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 583);

    auto ta1_z_yyyyz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 584);

    auto ta1_z_yyyyz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 585);

    auto ta1_z_yyyyz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 586);

    auto ta1_z_yyyyz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 587);

    auto ta1_z_yyyyz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 588);

    auto ta1_z_yyyyz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 589);

    auto ta1_z_yyyzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 590);

    auto ta1_z_yyyzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 591);

    auto ta1_z_yyyzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 592);

    auto ta1_z_yyyzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 593);

    auto ta1_z_yyyzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 594);

    auto ta1_z_yyyzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 595);

    auto ta1_z_yyyzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 596);

    auto ta1_z_yyyzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 597);

    auto ta1_z_yyyzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 598);

    auto ta1_z_yyyzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 599);

    auto ta1_z_yyzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 600);

    auto ta1_z_yyzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 601);

    auto ta1_z_yyzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 602);

    auto ta1_z_yyzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 603);

    auto ta1_z_yyzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 604);

    auto ta1_z_yyzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 605);

    auto ta1_z_yyzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 606);

    auto ta1_z_yyzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 607);

    auto ta1_z_yyzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 608);

    auto ta1_z_yyzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 609);

    auto ta1_z_yzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 610);

    auto ta1_z_yzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 611);

    auto ta1_z_yzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 612);

    auto ta1_z_yzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 613);

    auto ta1_z_yzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 614);

    auto ta1_z_yzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 615);

    auto ta1_z_yzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 616);

    auto ta1_z_yzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 617);

    auto ta1_z_yzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 618);

    auto ta1_z_yzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 619);

    auto ta1_z_zzzzz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_hf + 620);

    auto ta1_z_zzzzz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 621);

    auto ta1_z_zzzzz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 622);

    auto ta1_z_zzzzz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 623);

    auto ta1_z_zzzzz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 624);

    auto ta1_z_zzzzz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 625);

    auto ta1_z_zzzzz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_hf + 626);

    auto ta1_z_zzzzz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 627);

    auto ta1_z_zzzzz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 628);

    auto ta1_z_zzzzz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_hf + 629);

    // Set up 0-10 components of targeted buffer : IF

    auto ta1_x_xxxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if);

    auto ta1_x_xxxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 1);

    auto ta1_x_xxxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 2);

    auto ta1_x_xxxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 3);

    auto ta1_x_xxxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 4);

    auto ta1_x_xxxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 5);

    auto ta1_x_xxxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 6);

    auto ta1_x_xxxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 7);

    auto ta1_x_xxxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 8);

    auto ta1_x_xxxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 9);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_x_xxxx_xxx_0,   \
                             ta1_x_xxxx_xxx_1,   \
                             ta1_x_xxxx_xxy_0,   \
                             ta1_x_xxxx_xxy_1,   \
                             ta1_x_xxxx_xxz_0,   \
                             ta1_x_xxxx_xxz_1,   \
                             ta1_x_xxxx_xyy_0,   \
                             ta1_x_xxxx_xyy_1,   \
                             ta1_x_xxxx_xyz_0,   \
                             ta1_x_xxxx_xyz_1,   \
                             ta1_x_xxxx_xzz_0,   \
                             ta1_x_xxxx_xzz_1,   \
                             ta1_x_xxxx_yyy_0,   \
                             ta1_x_xxxx_yyy_1,   \
                             ta1_x_xxxx_yyz_0,   \
                             ta1_x_xxxx_yyz_1,   \
                             ta1_x_xxxx_yzz_0,   \
                             ta1_x_xxxx_yzz_1,   \
                             ta1_x_xxxx_zzz_0,   \
                             ta1_x_xxxx_zzz_1,   \
                             ta1_x_xxxxx_xx_0,   \
                             ta1_x_xxxxx_xx_1,   \
                             ta1_x_xxxxx_xxx_0,  \
                             ta1_x_xxxxx_xxx_1,  \
                             ta1_x_xxxxx_xxy_0,  \
                             ta1_x_xxxxx_xxy_1,  \
                             ta1_x_xxxxx_xxz_0,  \
                             ta1_x_xxxxx_xxz_1,  \
                             ta1_x_xxxxx_xy_0,   \
                             ta1_x_xxxxx_xy_1,   \
                             ta1_x_xxxxx_xyy_0,  \
                             ta1_x_xxxxx_xyy_1,  \
                             ta1_x_xxxxx_xyz_0,  \
                             ta1_x_xxxxx_xyz_1,  \
                             ta1_x_xxxxx_xz_0,   \
                             ta1_x_xxxxx_xz_1,   \
                             ta1_x_xxxxx_xzz_0,  \
                             ta1_x_xxxxx_xzz_1,  \
                             ta1_x_xxxxx_yy_0,   \
                             ta1_x_xxxxx_yy_1,   \
                             ta1_x_xxxxx_yyy_0,  \
                             ta1_x_xxxxx_yyy_1,  \
                             ta1_x_xxxxx_yyz_0,  \
                             ta1_x_xxxxx_yyz_1,  \
                             ta1_x_xxxxx_yz_0,   \
                             ta1_x_xxxxx_yz_1,   \
                             ta1_x_xxxxx_yzz_0,  \
                             ta1_x_xxxxx_yzz_1,  \
                             ta1_x_xxxxx_zz_0,   \
                             ta1_x_xxxxx_zz_1,   \
                             ta1_x_xxxxx_zzz_0,  \
                             ta1_x_xxxxx_zzz_1,  \
                             ta1_x_xxxxxx_xxx_0, \
                             ta1_x_xxxxxx_xxy_0, \
                             ta1_x_xxxxxx_xxz_0, \
                             ta1_x_xxxxxx_xyy_0, \
                             ta1_x_xxxxxx_xyz_0, \
                             ta1_x_xxxxxx_xzz_0, \
                             ta1_x_xxxxxx_yyy_0, \
                             ta1_x_xxxxxx_yyz_0, \
                             ta1_x_xxxxxx_yzz_0, \
                             ta1_x_xxxxxx_zzz_0, \
                             ta_xxxxx_xxx_1,     \
                             ta_xxxxx_xxy_1,     \
                             ta_xxxxx_xxz_1,     \
                             ta_xxxxx_xyy_1,     \
                             ta_xxxxx_xyz_1,     \
                             ta_xxxxx_xzz_1,     \
                             ta_xxxxx_yyy_1,     \
                             ta_xxxxx_yyz_1,     \
                             ta_xxxxx_yzz_1,     \
                             ta_xxxxx_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxx_xxx_0[i] = 5.0 * ta1_x_xxxx_xxx_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxx_1[i] * fe_0 + 3.0 * ta1_x_xxxxx_xx_0[i] * fe_0 -
                                3.0 * ta1_x_xxxxx_xx_1[i] * fe_0 + ta_xxxxx_xxx_1[i] + ta1_x_xxxxx_xxx_0[i] * pa_x[i] -
                                ta1_x_xxxxx_xxx_1[i] * pc_x[i];

        ta1_x_xxxxxx_xxy_0[i] = 5.0 * ta1_x_xxxx_xxy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxy_1[i] * fe_0 + 2.0 * ta1_x_xxxxx_xy_0[i] * fe_0 -
                                2.0 * ta1_x_xxxxx_xy_1[i] * fe_0 + ta_xxxxx_xxy_1[i] + ta1_x_xxxxx_xxy_0[i] * pa_x[i] -
                                ta1_x_xxxxx_xxy_1[i] * pc_x[i];

        ta1_x_xxxxxx_xxz_0[i] = 5.0 * ta1_x_xxxx_xxz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xxz_1[i] * fe_0 + 2.0 * ta1_x_xxxxx_xz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxxx_xz_1[i] * fe_0 + ta_xxxxx_xxz_1[i] + ta1_x_xxxxx_xxz_0[i] * pa_x[i] -
                                ta1_x_xxxxx_xxz_1[i] * pc_x[i];

        ta1_x_xxxxxx_xyy_0[i] = 5.0 * ta1_x_xxxx_xyy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xyy_1[i] * fe_0 + ta1_x_xxxxx_yy_0[i] * fe_0 -
                                ta1_x_xxxxx_yy_1[i] * fe_0 + ta_xxxxx_xyy_1[i] + ta1_x_xxxxx_xyy_0[i] * pa_x[i] - ta1_x_xxxxx_xyy_1[i] * pc_x[i];

        ta1_x_xxxxxx_xyz_0[i] = 5.0 * ta1_x_xxxx_xyz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xyz_1[i] * fe_0 + ta1_x_xxxxx_yz_0[i] * fe_0 -
                                ta1_x_xxxxx_yz_1[i] * fe_0 + ta_xxxxx_xyz_1[i] + ta1_x_xxxxx_xyz_0[i] * pa_x[i] - ta1_x_xxxxx_xyz_1[i] * pc_x[i];

        ta1_x_xxxxxx_xzz_0[i] = 5.0 * ta1_x_xxxx_xzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xzz_1[i] * fe_0 + ta1_x_xxxxx_zz_0[i] * fe_0 -
                                ta1_x_xxxxx_zz_1[i] * fe_0 + ta_xxxxx_xzz_1[i] + ta1_x_xxxxx_xzz_0[i] * pa_x[i] - ta1_x_xxxxx_xzz_1[i] * pc_x[i];

        ta1_x_xxxxxx_yyy_0[i] = 5.0 * ta1_x_xxxx_yyy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yyy_1[i] * fe_0 + ta_xxxxx_yyy_1[i] +
                                ta1_x_xxxxx_yyy_0[i] * pa_x[i] - ta1_x_xxxxx_yyy_1[i] * pc_x[i];

        ta1_x_xxxxxx_yyz_0[i] = 5.0 * ta1_x_xxxx_yyz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yyz_1[i] * fe_0 + ta_xxxxx_yyz_1[i] +
                                ta1_x_xxxxx_yyz_0[i] * pa_x[i] - ta1_x_xxxxx_yyz_1[i] * pc_x[i];

        ta1_x_xxxxxx_yzz_0[i] = 5.0 * ta1_x_xxxx_yzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yzz_1[i] * fe_0 + ta_xxxxx_yzz_1[i] +
                                ta1_x_xxxxx_yzz_0[i] * pa_x[i] - ta1_x_xxxxx_yzz_1[i] * pc_x[i];

        ta1_x_xxxxxx_zzz_0[i] = 5.0 * ta1_x_xxxx_zzz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_zzz_1[i] * fe_0 + ta_xxxxx_zzz_1[i] +
                                ta1_x_xxxxx_zzz_0[i] * pa_x[i] - ta1_x_xxxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : IF

    auto ta1_x_xxxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 10);

    auto ta1_x_xxxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 11);

    auto ta1_x_xxxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 12);

    auto ta1_x_xxxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 13);

    auto ta1_x_xxxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 14);

    auto ta1_x_xxxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 15);

    auto ta1_x_xxxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 16);

    auto ta1_x_xxxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 17);

    auto ta1_x_xxxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 18);

    auto ta1_x_xxxxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 19);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_xxxxx_xx_0,   \
                             ta1_x_xxxxx_xx_1,   \
                             ta1_x_xxxxx_xxx_0,  \
                             ta1_x_xxxxx_xxx_1,  \
                             ta1_x_xxxxx_xxy_0,  \
                             ta1_x_xxxxx_xxy_1,  \
                             ta1_x_xxxxx_xxz_0,  \
                             ta1_x_xxxxx_xxz_1,  \
                             ta1_x_xxxxx_xy_0,   \
                             ta1_x_xxxxx_xy_1,   \
                             ta1_x_xxxxx_xyy_0,  \
                             ta1_x_xxxxx_xyy_1,  \
                             ta1_x_xxxxx_xyz_0,  \
                             ta1_x_xxxxx_xyz_1,  \
                             ta1_x_xxxxx_xz_0,   \
                             ta1_x_xxxxx_xz_1,   \
                             ta1_x_xxxxx_xzz_0,  \
                             ta1_x_xxxxx_xzz_1,  \
                             ta1_x_xxxxx_yy_0,   \
                             ta1_x_xxxxx_yy_1,   \
                             ta1_x_xxxxx_yyy_0,  \
                             ta1_x_xxxxx_yyy_1,  \
                             ta1_x_xxxxx_yyz_0,  \
                             ta1_x_xxxxx_yyz_1,  \
                             ta1_x_xxxxx_yz_0,   \
                             ta1_x_xxxxx_yz_1,   \
                             ta1_x_xxxxx_yzz_0,  \
                             ta1_x_xxxxx_yzz_1,  \
                             ta1_x_xxxxx_zz_0,   \
                             ta1_x_xxxxx_zz_1,   \
                             ta1_x_xxxxx_zzz_0,  \
                             ta1_x_xxxxx_zzz_1,  \
                             ta1_x_xxxxxy_xxx_0, \
                             ta1_x_xxxxxy_xxy_0, \
                             ta1_x_xxxxxy_xxz_0, \
                             ta1_x_xxxxxy_xyy_0, \
                             ta1_x_xxxxxy_xyz_0, \
                             ta1_x_xxxxxy_xzz_0, \
                             ta1_x_xxxxxy_yyy_0, \
                             ta1_x_xxxxxy_yyz_0, \
                             ta1_x_xxxxxy_yzz_0, \
                             ta1_x_xxxxxy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxy_xxx_0[i] = ta1_x_xxxxx_xxx_0[i] * pa_y[i] - ta1_x_xxxxx_xxx_1[i] * pc_y[i];

        ta1_x_xxxxxy_xxy_0[i] =
            ta1_x_xxxxx_xx_0[i] * fe_0 - ta1_x_xxxxx_xx_1[i] * fe_0 + ta1_x_xxxxx_xxy_0[i] * pa_y[i] - ta1_x_xxxxx_xxy_1[i] * pc_y[i];

        ta1_x_xxxxxy_xxz_0[i] = ta1_x_xxxxx_xxz_0[i] * pa_y[i] - ta1_x_xxxxx_xxz_1[i] * pc_y[i];

        ta1_x_xxxxxy_xyy_0[i] =
            2.0 * ta1_x_xxxxx_xy_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xy_1[i] * fe_0 + ta1_x_xxxxx_xyy_0[i] * pa_y[i] - ta1_x_xxxxx_xyy_1[i] * pc_y[i];

        ta1_x_xxxxxy_xyz_0[i] =
            ta1_x_xxxxx_xz_0[i] * fe_0 - ta1_x_xxxxx_xz_1[i] * fe_0 + ta1_x_xxxxx_xyz_0[i] * pa_y[i] - ta1_x_xxxxx_xyz_1[i] * pc_y[i];

        ta1_x_xxxxxy_xzz_0[i] = ta1_x_xxxxx_xzz_0[i] * pa_y[i] - ta1_x_xxxxx_xzz_1[i] * pc_y[i];

        ta1_x_xxxxxy_yyy_0[i] =
            3.0 * ta1_x_xxxxx_yy_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_yy_1[i] * fe_0 + ta1_x_xxxxx_yyy_0[i] * pa_y[i] - ta1_x_xxxxx_yyy_1[i] * pc_y[i];

        ta1_x_xxxxxy_yyz_0[i] =
            2.0 * ta1_x_xxxxx_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_yz_1[i] * fe_0 + ta1_x_xxxxx_yyz_0[i] * pa_y[i] - ta1_x_xxxxx_yyz_1[i] * pc_y[i];

        ta1_x_xxxxxy_yzz_0[i] =
            ta1_x_xxxxx_zz_0[i] * fe_0 - ta1_x_xxxxx_zz_1[i] * fe_0 + ta1_x_xxxxx_yzz_0[i] * pa_y[i] - ta1_x_xxxxx_yzz_1[i] * pc_y[i];

        ta1_x_xxxxxy_zzz_0[i] = ta1_x_xxxxx_zzz_0[i] * pa_y[i] - ta1_x_xxxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : IF

    auto ta1_x_xxxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 20);

    auto ta1_x_xxxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 21);

    auto ta1_x_xxxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 22);

    auto ta1_x_xxxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 23);

    auto ta1_x_xxxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 24);

    auto ta1_x_xxxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 25);

    auto ta1_x_xxxxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 26);

    auto ta1_x_xxxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 27);

    auto ta1_x_xxxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 28);

    auto ta1_x_xxxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 29);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_xxxxx_xx_0,   \
                             ta1_x_xxxxx_xx_1,   \
                             ta1_x_xxxxx_xxx_0,  \
                             ta1_x_xxxxx_xxx_1,  \
                             ta1_x_xxxxx_xxy_0,  \
                             ta1_x_xxxxx_xxy_1,  \
                             ta1_x_xxxxx_xxz_0,  \
                             ta1_x_xxxxx_xxz_1,  \
                             ta1_x_xxxxx_xy_0,   \
                             ta1_x_xxxxx_xy_1,   \
                             ta1_x_xxxxx_xyy_0,  \
                             ta1_x_xxxxx_xyy_1,  \
                             ta1_x_xxxxx_xyz_0,  \
                             ta1_x_xxxxx_xyz_1,  \
                             ta1_x_xxxxx_xz_0,   \
                             ta1_x_xxxxx_xz_1,   \
                             ta1_x_xxxxx_xzz_0,  \
                             ta1_x_xxxxx_xzz_1,  \
                             ta1_x_xxxxx_yy_0,   \
                             ta1_x_xxxxx_yy_1,   \
                             ta1_x_xxxxx_yyy_0,  \
                             ta1_x_xxxxx_yyy_1,  \
                             ta1_x_xxxxx_yyz_0,  \
                             ta1_x_xxxxx_yyz_1,  \
                             ta1_x_xxxxx_yz_0,   \
                             ta1_x_xxxxx_yz_1,   \
                             ta1_x_xxxxx_yzz_0,  \
                             ta1_x_xxxxx_yzz_1,  \
                             ta1_x_xxxxx_zz_0,   \
                             ta1_x_xxxxx_zz_1,   \
                             ta1_x_xxxxx_zzz_0,  \
                             ta1_x_xxxxx_zzz_1,  \
                             ta1_x_xxxxxz_xxx_0, \
                             ta1_x_xxxxxz_xxy_0, \
                             ta1_x_xxxxxz_xxz_0, \
                             ta1_x_xxxxxz_xyy_0, \
                             ta1_x_xxxxxz_xyz_0, \
                             ta1_x_xxxxxz_xzz_0, \
                             ta1_x_xxxxxz_yyy_0, \
                             ta1_x_xxxxxz_yyz_0, \
                             ta1_x_xxxxxz_yzz_0, \
                             ta1_x_xxxxxz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxz_xxx_0[i] = ta1_x_xxxxx_xxx_0[i] * pa_z[i] - ta1_x_xxxxx_xxx_1[i] * pc_z[i];

        ta1_x_xxxxxz_xxy_0[i] = ta1_x_xxxxx_xxy_0[i] * pa_z[i] - ta1_x_xxxxx_xxy_1[i] * pc_z[i];

        ta1_x_xxxxxz_xxz_0[i] =
            ta1_x_xxxxx_xx_0[i] * fe_0 - ta1_x_xxxxx_xx_1[i] * fe_0 + ta1_x_xxxxx_xxz_0[i] * pa_z[i] - ta1_x_xxxxx_xxz_1[i] * pc_z[i];

        ta1_x_xxxxxz_xyy_0[i] = ta1_x_xxxxx_xyy_0[i] * pa_z[i] - ta1_x_xxxxx_xyy_1[i] * pc_z[i];

        ta1_x_xxxxxz_xyz_0[i] =
            ta1_x_xxxxx_xy_0[i] * fe_0 - ta1_x_xxxxx_xy_1[i] * fe_0 + ta1_x_xxxxx_xyz_0[i] * pa_z[i] - ta1_x_xxxxx_xyz_1[i] * pc_z[i];

        ta1_x_xxxxxz_xzz_0[i] =
            2.0 * ta1_x_xxxxx_xz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_xz_1[i] * fe_0 + ta1_x_xxxxx_xzz_0[i] * pa_z[i] - ta1_x_xxxxx_xzz_1[i] * pc_z[i];

        ta1_x_xxxxxz_yyy_0[i] = ta1_x_xxxxx_yyy_0[i] * pa_z[i] - ta1_x_xxxxx_yyy_1[i] * pc_z[i];

        ta1_x_xxxxxz_yyz_0[i] =
            ta1_x_xxxxx_yy_0[i] * fe_0 - ta1_x_xxxxx_yy_1[i] * fe_0 + ta1_x_xxxxx_yyz_0[i] * pa_z[i] - ta1_x_xxxxx_yyz_1[i] * pc_z[i];

        ta1_x_xxxxxz_yzz_0[i] =
            2.0 * ta1_x_xxxxx_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_yz_1[i] * fe_0 + ta1_x_xxxxx_yzz_0[i] * pa_z[i] - ta1_x_xxxxx_yzz_1[i] * pc_z[i];

        ta1_x_xxxxxz_zzz_0[i] =
            3.0 * ta1_x_xxxxx_zz_0[i] * fe_0 - 3.0 * ta1_x_xxxxx_zz_1[i] * fe_0 + ta1_x_xxxxx_zzz_0[i] * pa_z[i] - ta1_x_xxxxx_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : IF

    auto ta1_x_xxxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 30);

    auto ta1_x_xxxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 31);

    auto ta1_x_xxxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 32);

    auto ta1_x_xxxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 33);

    auto ta1_x_xxxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 34);

    auto ta1_x_xxxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 35);

    auto ta1_x_xxxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 36);

    auto ta1_x_xxxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 37);

    auto ta1_x_xxxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 38);

    auto ta1_x_xxxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 39);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xxxx_xxx_0,   \
                             ta1_x_xxxx_xxx_1,   \
                             ta1_x_xxxx_xxy_0,   \
                             ta1_x_xxxx_xxy_1,   \
                             ta1_x_xxxx_xxz_0,   \
                             ta1_x_xxxx_xxz_1,   \
                             ta1_x_xxxx_xyy_0,   \
                             ta1_x_xxxx_xyy_1,   \
                             ta1_x_xxxx_xyz_0,   \
                             ta1_x_xxxx_xyz_1,   \
                             ta1_x_xxxx_xzz_0,   \
                             ta1_x_xxxx_xzz_1,   \
                             ta1_x_xxxx_zzz_0,   \
                             ta1_x_xxxx_zzz_1,   \
                             ta1_x_xxxxy_xx_0,   \
                             ta1_x_xxxxy_xx_1,   \
                             ta1_x_xxxxy_xxx_0,  \
                             ta1_x_xxxxy_xxx_1,  \
                             ta1_x_xxxxy_xxy_0,  \
                             ta1_x_xxxxy_xxy_1,  \
                             ta1_x_xxxxy_xxz_0,  \
                             ta1_x_xxxxy_xxz_1,  \
                             ta1_x_xxxxy_xy_0,   \
                             ta1_x_xxxxy_xy_1,   \
                             ta1_x_xxxxy_xyy_0,  \
                             ta1_x_xxxxy_xyy_1,  \
                             ta1_x_xxxxy_xyz_0,  \
                             ta1_x_xxxxy_xyz_1,  \
                             ta1_x_xxxxy_xz_0,   \
                             ta1_x_xxxxy_xz_1,   \
                             ta1_x_xxxxy_xzz_0,  \
                             ta1_x_xxxxy_xzz_1,  \
                             ta1_x_xxxxy_zzz_0,  \
                             ta1_x_xxxxy_zzz_1,  \
                             ta1_x_xxxxyy_xxx_0, \
                             ta1_x_xxxxyy_xxy_0, \
                             ta1_x_xxxxyy_xxz_0, \
                             ta1_x_xxxxyy_xyy_0, \
                             ta1_x_xxxxyy_xyz_0, \
                             ta1_x_xxxxyy_xzz_0, \
                             ta1_x_xxxxyy_yyy_0, \
                             ta1_x_xxxxyy_yyz_0, \
                             ta1_x_xxxxyy_yzz_0, \
                             ta1_x_xxxxyy_zzz_0, \
                             ta1_x_xxxyy_yyy_0,  \
                             ta1_x_xxxyy_yyy_1,  \
                             ta1_x_xxxyy_yyz_0,  \
                             ta1_x_xxxyy_yyz_1,  \
                             ta1_x_xxxyy_yzz_0,  \
                             ta1_x_xxxyy_yzz_1,  \
                             ta1_x_xxyy_yyy_0,   \
                             ta1_x_xxyy_yyy_1,   \
                             ta1_x_xxyy_yyz_0,   \
                             ta1_x_xxyy_yyz_1,   \
                             ta1_x_xxyy_yzz_0,   \
                             ta1_x_xxyy_yzz_1,   \
                             ta_xxxyy_yyy_1,     \
                             ta_xxxyy_yyz_1,     \
                             ta_xxxyy_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxyy_xxx_0[i] =
            ta1_x_xxxx_xxx_0[i] * fe_0 - ta1_x_xxxx_xxx_1[i] * fe_0 + ta1_x_xxxxy_xxx_0[i] * pa_y[i] - ta1_x_xxxxy_xxx_1[i] * pc_y[i];

        ta1_x_xxxxyy_xxy_0[i] = ta1_x_xxxx_xxy_0[i] * fe_0 - ta1_x_xxxx_xxy_1[i] * fe_0 + ta1_x_xxxxy_xx_0[i] * fe_0 - ta1_x_xxxxy_xx_1[i] * fe_0 +
                                ta1_x_xxxxy_xxy_0[i] * pa_y[i] - ta1_x_xxxxy_xxy_1[i] * pc_y[i];

        ta1_x_xxxxyy_xxz_0[i] =
            ta1_x_xxxx_xxz_0[i] * fe_0 - ta1_x_xxxx_xxz_1[i] * fe_0 + ta1_x_xxxxy_xxz_0[i] * pa_y[i] - ta1_x_xxxxy_xxz_1[i] * pc_y[i];

        ta1_x_xxxxyy_xyy_0[i] = ta1_x_xxxx_xyy_0[i] * fe_0 - ta1_x_xxxx_xyy_1[i] * fe_0 + 2.0 * ta1_x_xxxxy_xy_0[i] * fe_0 -
                                2.0 * ta1_x_xxxxy_xy_1[i] * fe_0 + ta1_x_xxxxy_xyy_0[i] * pa_y[i] - ta1_x_xxxxy_xyy_1[i] * pc_y[i];

        ta1_x_xxxxyy_xyz_0[i] = ta1_x_xxxx_xyz_0[i] * fe_0 - ta1_x_xxxx_xyz_1[i] * fe_0 + ta1_x_xxxxy_xz_0[i] * fe_0 - ta1_x_xxxxy_xz_1[i] * fe_0 +
                                ta1_x_xxxxy_xyz_0[i] * pa_y[i] - ta1_x_xxxxy_xyz_1[i] * pc_y[i];

        ta1_x_xxxxyy_xzz_0[i] =
            ta1_x_xxxx_xzz_0[i] * fe_0 - ta1_x_xxxx_xzz_1[i] * fe_0 + ta1_x_xxxxy_xzz_0[i] * pa_y[i] - ta1_x_xxxxy_xzz_1[i] * pc_y[i];

        ta1_x_xxxxyy_yyy_0[i] = 3.0 * ta1_x_xxyy_yyy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yyy_1[i] * fe_0 + ta_xxxyy_yyy_1[i] +
                                ta1_x_xxxyy_yyy_0[i] * pa_x[i] - ta1_x_xxxyy_yyy_1[i] * pc_x[i];

        ta1_x_xxxxyy_yyz_0[i] = 3.0 * ta1_x_xxyy_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yyz_1[i] * fe_0 + ta_xxxyy_yyz_1[i] +
                                ta1_x_xxxyy_yyz_0[i] * pa_x[i] - ta1_x_xxxyy_yyz_1[i] * pc_x[i];

        ta1_x_xxxxyy_yzz_0[i] = 3.0 * ta1_x_xxyy_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yzz_1[i] * fe_0 + ta_xxxyy_yzz_1[i] +
                                ta1_x_xxxyy_yzz_0[i] * pa_x[i] - ta1_x_xxxyy_yzz_1[i] * pc_x[i];

        ta1_x_xxxxyy_zzz_0[i] =
            ta1_x_xxxx_zzz_0[i] * fe_0 - ta1_x_xxxx_zzz_1[i] * fe_0 + ta1_x_xxxxy_zzz_0[i] * pa_y[i] - ta1_x_xxxxy_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : IF

    auto ta1_x_xxxxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 40);

    auto ta1_x_xxxxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 41);

    auto ta1_x_xxxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 42);

    auto ta1_x_xxxxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 43);

    auto ta1_x_xxxxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 44);

    auto ta1_x_xxxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 45);

    auto ta1_x_xxxxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 46);

    auto ta1_x_xxxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 47);

    auto ta1_x_xxxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 48);

    auto ta1_x_xxxxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 49);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxxxy_xxy_0,  \
                             ta1_x_xxxxy_xxy_1,  \
                             ta1_x_xxxxy_xyy_0,  \
                             ta1_x_xxxxy_xyy_1,  \
                             ta1_x_xxxxy_yyy_0,  \
                             ta1_x_xxxxy_yyy_1,  \
                             ta1_x_xxxxyz_xxx_0, \
                             ta1_x_xxxxyz_xxy_0, \
                             ta1_x_xxxxyz_xxz_0, \
                             ta1_x_xxxxyz_xyy_0, \
                             ta1_x_xxxxyz_xyz_0, \
                             ta1_x_xxxxyz_xzz_0, \
                             ta1_x_xxxxyz_yyy_0, \
                             ta1_x_xxxxyz_yyz_0, \
                             ta1_x_xxxxyz_yzz_0, \
                             ta1_x_xxxxyz_zzz_0, \
                             ta1_x_xxxxz_xxx_0,  \
                             ta1_x_xxxxz_xxx_1,  \
                             ta1_x_xxxxz_xxz_0,  \
                             ta1_x_xxxxz_xxz_1,  \
                             ta1_x_xxxxz_xyz_0,  \
                             ta1_x_xxxxz_xyz_1,  \
                             ta1_x_xxxxz_xz_0,   \
                             ta1_x_xxxxz_xz_1,   \
                             ta1_x_xxxxz_xzz_0,  \
                             ta1_x_xxxxz_xzz_1,  \
                             ta1_x_xxxxz_yyz_0,  \
                             ta1_x_xxxxz_yyz_1,  \
                             ta1_x_xxxxz_yz_0,   \
                             ta1_x_xxxxz_yz_1,   \
                             ta1_x_xxxxz_yzz_0,  \
                             ta1_x_xxxxz_yzz_1,  \
                             ta1_x_xxxxz_zz_0,   \
                             ta1_x_xxxxz_zz_1,   \
                             ta1_x_xxxxz_zzz_0,  \
                             ta1_x_xxxxz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxyz_xxx_0[i] = ta1_x_xxxxz_xxx_0[i] * pa_y[i] - ta1_x_xxxxz_xxx_1[i] * pc_y[i];

        ta1_x_xxxxyz_xxy_0[i] = ta1_x_xxxxy_xxy_0[i] * pa_z[i] - ta1_x_xxxxy_xxy_1[i] * pc_z[i];

        ta1_x_xxxxyz_xxz_0[i] = ta1_x_xxxxz_xxz_0[i] * pa_y[i] - ta1_x_xxxxz_xxz_1[i] * pc_y[i];

        ta1_x_xxxxyz_xyy_0[i] = ta1_x_xxxxy_xyy_0[i] * pa_z[i] - ta1_x_xxxxy_xyy_1[i] * pc_z[i];

        ta1_x_xxxxyz_xyz_0[i] =
            ta1_x_xxxxz_xz_0[i] * fe_0 - ta1_x_xxxxz_xz_1[i] * fe_0 + ta1_x_xxxxz_xyz_0[i] * pa_y[i] - ta1_x_xxxxz_xyz_1[i] * pc_y[i];

        ta1_x_xxxxyz_xzz_0[i] = ta1_x_xxxxz_xzz_0[i] * pa_y[i] - ta1_x_xxxxz_xzz_1[i] * pc_y[i];

        ta1_x_xxxxyz_yyy_0[i] = ta1_x_xxxxy_yyy_0[i] * pa_z[i] - ta1_x_xxxxy_yyy_1[i] * pc_z[i];

        ta1_x_xxxxyz_yyz_0[i] =
            2.0 * ta1_x_xxxxz_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxxz_yz_1[i] * fe_0 + ta1_x_xxxxz_yyz_0[i] * pa_y[i] - ta1_x_xxxxz_yyz_1[i] * pc_y[i];

        ta1_x_xxxxyz_yzz_0[i] =
            ta1_x_xxxxz_zz_0[i] * fe_0 - ta1_x_xxxxz_zz_1[i] * fe_0 + ta1_x_xxxxz_yzz_0[i] * pa_y[i] - ta1_x_xxxxz_yzz_1[i] * pc_y[i];

        ta1_x_xxxxyz_zzz_0[i] = ta1_x_xxxxz_zzz_0[i] * pa_y[i] - ta1_x_xxxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : IF

    auto ta1_x_xxxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 50);

    auto ta1_x_xxxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 51);

    auto ta1_x_xxxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 52);

    auto ta1_x_xxxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 53);

    auto ta1_x_xxxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 54);

    auto ta1_x_xxxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 55);

    auto ta1_x_xxxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 56);

    auto ta1_x_xxxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 57);

    auto ta1_x_xxxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 58);

    auto ta1_x_xxxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 59);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xxxx_xxx_0,   \
                             ta1_x_xxxx_xxx_1,   \
                             ta1_x_xxxx_xxy_0,   \
                             ta1_x_xxxx_xxy_1,   \
                             ta1_x_xxxx_xxz_0,   \
                             ta1_x_xxxx_xxz_1,   \
                             ta1_x_xxxx_xyy_0,   \
                             ta1_x_xxxx_xyy_1,   \
                             ta1_x_xxxx_xyz_0,   \
                             ta1_x_xxxx_xyz_1,   \
                             ta1_x_xxxx_xzz_0,   \
                             ta1_x_xxxx_xzz_1,   \
                             ta1_x_xxxx_yyy_0,   \
                             ta1_x_xxxx_yyy_1,   \
                             ta1_x_xxxxz_xx_0,   \
                             ta1_x_xxxxz_xx_1,   \
                             ta1_x_xxxxz_xxx_0,  \
                             ta1_x_xxxxz_xxx_1,  \
                             ta1_x_xxxxz_xxy_0,  \
                             ta1_x_xxxxz_xxy_1,  \
                             ta1_x_xxxxz_xxz_0,  \
                             ta1_x_xxxxz_xxz_1,  \
                             ta1_x_xxxxz_xy_0,   \
                             ta1_x_xxxxz_xy_1,   \
                             ta1_x_xxxxz_xyy_0,  \
                             ta1_x_xxxxz_xyy_1,  \
                             ta1_x_xxxxz_xyz_0,  \
                             ta1_x_xxxxz_xyz_1,  \
                             ta1_x_xxxxz_xz_0,   \
                             ta1_x_xxxxz_xz_1,   \
                             ta1_x_xxxxz_xzz_0,  \
                             ta1_x_xxxxz_xzz_1,  \
                             ta1_x_xxxxz_yyy_0,  \
                             ta1_x_xxxxz_yyy_1,  \
                             ta1_x_xxxxzz_xxx_0, \
                             ta1_x_xxxxzz_xxy_0, \
                             ta1_x_xxxxzz_xxz_0, \
                             ta1_x_xxxxzz_xyy_0, \
                             ta1_x_xxxxzz_xyz_0, \
                             ta1_x_xxxxzz_xzz_0, \
                             ta1_x_xxxxzz_yyy_0, \
                             ta1_x_xxxxzz_yyz_0, \
                             ta1_x_xxxxzz_yzz_0, \
                             ta1_x_xxxxzz_zzz_0, \
                             ta1_x_xxxzz_yyz_0,  \
                             ta1_x_xxxzz_yyz_1,  \
                             ta1_x_xxxzz_yzz_0,  \
                             ta1_x_xxxzz_yzz_1,  \
                             ta1_x_xxxzz_zzz_0,  \
                             ta1_x_xxxzz_zzz_1,  \
                             ta1_x_xxzz_yyz_0,   \
                             ta1_x_xxzz_yyz_1,   \
                             ta1_x_xxzz_yzz_0,   \
                             ta1_x_xxzz_yzz_1,   \
                             ta1_x_xxzz_zzz_0,   \
                             ta1_x_xxzz_zzz_1,   \
                             ta_xxxzz_yyz_1,     \
                             ta_xxxzz_yzz_1,     \
                             ta_xxxzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxzz_xxx_0[i] =
            ta1_x_xxxx_xxx_0[i] * fe_0 - ta1_x_xxxx_xxx_1[i] * fe_0 + ta1_x_xxxxz_xxx_0[i] * pa_z[i] - ta1_x_xxxxz_xxx_1[i] * pc_z[i];

        ta1_x_xxxxzz_xxy_0[i] =
            ta1_x_xxxx_xxy_0[i] * fe_0 - ta1_x_xxxx_xxy_1[i] * fe_0 + ta1_x_xxxxz_xxy_0[i] * pa_z[i] - ta1_x_xxxxz_xxy_1[i] * pc_z[i];

        ta1_x_xxxxzz_xxz_0[i] = ta1_x_xxxx_xxz_0[i] * fe_0 - ta1_x_xxxx_xxz_1[i] * fe_0 + ta1_x_xxxxz_xx_0[i] * fe_0 - ta1_x_xxxxz_xx_1[i] * fe_0 +
                                ta1_x_xxxxz_xxz_0[i] * pa_z[i] - ta1_x_xxxxz_xxz_1[i] * pc_z[i];

        ta1_x_xxxxzz_xyy_0[i] =
            ta1_x_xxxx_xyy_0[i] * fe_0 - ta1_x_xxxx_xyy_1[i] * fe_0 + ta1_x_xxxxz_xyy_0[i] * pa_z[i] - ta1_x_xxxxz_xyy_1[i] * pc_z[i];

        ta1_x_xxxxzz_xyz_0[i] = ta1_x_xxxx_xyz_0[i] * fe_0 - ta1_x_xxxx_xyz_1[i] * fe_0 + ta1_x_xxxxz_xy_0[i] * fe_0 - ta1_x_xxxxz_xy_1[i] * fe_0 +
                                ta1_x_xxxxz_xyz_0[i] * pa_z[i] - ta1_x_xxxxz_xyz_1[i] * pc_z[i];

        ta1_x_xxxxzz_xzz_0[i] = ta1_x_xxxx_xzz_0[i] * fe_0 - ta1_x_xxxx_xzz_1[i] * fe_0 + 2.0 * ta1_x_xxxxz_xz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxxz_xz_1[i] * fe_0 + ta1_x_xxxxz_xzz_0[i] * pa_z[i] - ta1_x_xxxxz_xzz_1[i] * pc_z[i];

        ta1_x_xxxxzz_yyy_0[i] =
            ta1_x_xxxx_yyy_0[i] * fe_0 - ta1_x_xxxx_yyy_1[i] * fe_0 + ta1_x_xxxxz_yyy_0[i] * pa_z[i] - ta1_x_xxxxz_yyy_1[i] * pc_z[i];

        ta1_x_xxxxzz_yyz_0[i] = 3.0 * ta1_x_xxzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yyz_1[i] * fe_0 + ta_xxxzz_yyz_1[i] +
                                ta1_x_xxxzz_yyz_0[i] * pa_x[i] - ta1_x_xxxzz_yyz_1[i] * pc_x[i];

        ta1_x_xxxxzz_yzz_0[i] = 3.0 * ta1_x_xxzz_yzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yzz_1[i] * fe_0 + ta_xxxzz_yzz_1[i] +
                                ta1_x_xxxzz_yzz_0[i] * pa_x[i] - ta1_x_xxxzz_yzz_1[i] * pc_x[i];

        ta1_x_xxxxzz_zzz_0[i] = 3.0 * ta1_x_xxzz_zzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_zzz_1[i] * fe_0 + ta_xxxzz_zzz_1[i] +
                                ta1_x_xxxzz_zzz_0[i] * pa_x[i] - ta1_x_xxxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : IF

    auto ta1_x_xxxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 60);

    auto ta1_x_xxxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 61);

    auto ta1_x_xxxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 62);

    auto ta1_x_xxxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 63);

    auto ta1_x_xxxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 64);

    auto ta1_x_xxxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 65);

    auto ta1_x_xxxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 66);

    auto ta1_x_xxxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 67);

    auto ta1_x_xxxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 68);

    auto ta1_x_xxxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 69);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xxxy_xxx_0,   \
                             ta1_x_xxxy_xxx_1,   \
                             ta1_x_xxxy_xxy_0,   \
                             ta1_x_xxxy_xxy_1,   \
                             ta1_x_xxxy_xxz_0,   \
                             ta1_x_xxxy_xxz_1,   \
                             ta1_x_xxxy_xyy_0,   \
                             ta1_x_xxxy_xyy_1,   \
                             ta1_x_xxxy_xyz_0,   \
                             ta1_x_xxxy_xyz_1,   \
                             ta1_x_xxxy_xzz_0,   \
                             ta1_x_xxxy_xzz_1,   \
                             ta1_x_xxxy_zzz_0,   \
                             ta1_x_xxxy_zzz_1,   \
                             ta1_x_xxxyy_xx_0,   \
                             ta1_x_xxxyy_xx_1,   \
                             ta1_x_xxxyy_xxx_0,  \
                             ta1_x_xxxyy_xxx_1,  \
                             ta1_x_xxxyy_xxy_0,  \
                             ta1_x_xxxyy_xxy_1,  \
                             ta1_x_xxxyy_xxz_0,  \
                             ta1_x_xxxyy_xxz_1,  \
                             ta1_x_xxxyy_xy_0,   \
                             ta1_x_xxxyy_xy_1,   \
                             ta1_x_xxxyy_xyy_0,  \
                             ta1_x_xxxyy_xyy_1,  \
                             ta1_x_xxxyy_xyz_0,  \
                             ta1_x_xxxyy_xyz_1,  \
                             ta1_x_xxxyy_xz_0,   \
                             ta1_x_xxxyy_xz_1,   \
                             ta1_x_xxxyy_xzz_0,  \
                             ta1_x_xxxyy_xzz_1,  \
                             ta1_x_xxxyy_zzz_0,  \
                             ta1_x_xxxyy_zzz_1,  \
                             ta1_x_xxxyyy_xxx_0, \
                             ta1_x_xxxyyy_xxy_0, \
                             ta1_x_xxxyyy_xxz_0, \
                             ta1_x_xxxyyy_xyy_0, \
                             ta1_x_xxxyyy_xyz_0, \
                             ta1_x_xxxyyy_xzz_0, \
                             ta1_x_xxxyyy_yyy_0, \
                             ta1_x_xxxyyy_yyz_0, \
                             ta1_x_xxxyyy_yzz_0, \
                             ta1_x_xxxyyy_zzz_0, \
                             ta1_x_xxyyy_yyy_0,  \
                             ta1_x_xxyyy_yyy_1,  \
                             ta1_x_xxyyy_yyz_0,  \
                             ta1_x_xxyyy_yyz_1,  \
                             ta1_x_xxyyy_yzz_0,  \
                             ta1_x_xxyyy_yzz_1,  \
                             ta1_x_xyyy_yyy_0,   \
                             ta1_x_xyyy_yyy_1,   \
                             ta1_x_xyyy_yyz_0,   \
                             ta1_x_xyyy_yyz_1,   \
                             ta1_x_xyyy_yzz_0,   \
                             ta1_x_xyyy_yzz_1,   \
                             ta_xxyyy_yyy_1,     \
                             ta_xxyyy_yyz_1,     \
                             ta_xxyyy_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyy_xxx_0[i] =
            2.0 * ta1_x_xxxy_xxx_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxx_1[i] * fe_0 + ta1_x_xxxyy_xxx_0[i] * pa_y[i] - ta1_x_xxxyy_xxx_1[i] * pc_y[i];

        ta1_x_xxxyyy_xxy_0[i] = 2.0 * ta1_x_xxxy_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxy_1[i] * fe_0 + ta1_x_xxxyy_xx_0[i] * fe_0 -
                                ta1_x_xxxyy_xx_1[i] * fe_0 + ta1_x_xxxyy_xxy_0[i] * pa_y[i] - ta1_x_xxxyy_xxy_1[i] * pc_y[i];

        ta1_x_xxxyyy_xxz_0[i] =
            2.0 * ta1_x_xxxy_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xxz_1[i] * fe_0 + ta1_x_xxxyy_xxz_0[i] * pa_y[i] - ta1_x_xxxyy_xxz_1[i] * pc_y[i];

        ta1_x_xxxyyy_xyy_0[i] = 2.0 * ta1_x_xxxy_xyy_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xyy_1[i] * fe_0 + 2.0 * ta1_x_xxxyy_xy_0[i] * fe_0 -
                                2.0 * ta1_x_xxxyy_xy_1[i] * fe_0 + ta1_x_xxxyy_xyy_0[i] * pa_y[i] - ta1_x_xxxyy_xyy_1[i] * pc_y[i];

        ta1_x_xxxyyy_xyz_0[i] = 2.0 * ta1_x_xxxy_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xyz_1[i] * fe_0 + ta1_x_xxxyy_xz_0[i] * fe_0 -
                                ta1_x_xxxyy_xz_1[i] * fe_0 + ta1_x_xxxyy_xyz_0[i] * pa_y[i] - ta1_x_xxxyy_xyz_1[i] * pc_y[i];

        ta1_x_xxxyyy_xzz_0[i] =
            2.0 * ta1_x_xxxy_xzz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xzz_1[i] * fe_0 + ta1_x_xxxyy_xzz_0[i] * pa_y[i] - ta1_x_xxxyy_xzz_1[i] * pc_y[i];

        ta1_x_xxxyyy_yyy_0[i] = 2.0 * ta1_x_xyyy_yyy_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yyy_1[i] * fe_0 + ta_xxyyy_yyy_1[i] +
                                ta1_x_xxyyy_yyy_0[i] * pa_x[i] - ta1_x_xxyyy_yyy_1[i] * pc_x[i];

        ta1_x_xxxyyy_yyz_0[i] = 2.0 * ta1_x_xyyy_yyz_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yyz_1[i] * fe_0 + ta_xxyyy_yyz_1[i] +
                                ta1_x_xxyyy_yyz_0[i] * pa_x[i] - ta1_x_xxyyy_yyz_1[i] * pc_x[i];

        ta1_x_xxxyyy_yzz_0[i] = 2.0 * ta1_x_xyyy_yzz_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yzz_1[i] * fe_0 + ta_xxyyy_yzz_1[i] +
                                ta1_x_xxyyy_yzz_0[i] * pa_x[i] - ta1_x_xxyyy_yzz_1[i] * pc_x[i];

        ta1_x_xxxyyy_zzz_0[i] =
            2.0 * ta1_x_xxxy_zzz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_zzz_1[i] * fe_0 + ta1_x_xxxyy_zzz_0[i] * pa_y[i] - ta1_x_xxxyy_zzz_1[i] * pc_y[i];
    }

    // Set up 70-80 components of targeted buffer : IF

    auto ta1_x_xxxyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 70);

    auto ta1_x_xxxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 71);

    auto ta1_x_xxxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 72);

    auto ta1_x_xxxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 73);

    auto ta1_x_xxxyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 74);

    auto ta1_x_xxxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 75);

    auto ta1_x_xxxyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 76);

    auto ta1_x_xxxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 77);

    auto ta1_x_xxxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 78);

    auto ta1_x_xxxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 79);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxxyy_xxx_0,  \
                             ta1_x_xxxyy_xxx_1,  \
                             ta1_x_xxxyy_xxy_0,  \
                             ta1_x_xxxyy_xxy_1,  \
                             ta1_x_xxxyy_xy_0,   \
                             ta1_x_xxxyy_xy_1,   \
                             ta1_x_xxxyy_xyy_0,  \
                             ta1_x_xxxyy_xyy_1,  \
                             ta1_x_xxxyy_xyz_0,  \
                             ta1_x_xxxyy_xyz_1,  \
                             ta1_x_xxxyy_yy_0,   \
                             ta1_x_xxxyy_yy_1,   \
                             ta1_x_xxxyy_yyy_0,  \
                             ta1_x_xxxyy_yyy_1,  \
                             ta1_x_xxxyy_yyz_0,  \
                             ta1_x_xxxyy_yyz_1,  \
                             ta1_x_xxxyy_yz_0,   \
                             ta1_x_xxxyy_yz_1,   \
                             ta1_x_xxxyy_yzz_0,  \
                             ta1_x_xxxyy_yzz_1,  \
                             ta1_x_xxxyyz_xxx_0, \
                             ta1_x_xxxyyz_xxy_0, \
                             ta1_x_xxxyyz_xxz_0, \
                             ta1_x_xxxyyz_xyy_0, \
                             ta1_x_xxxyyz_xyz_0, \
                             ta1_x_xxxyyz_xzz_0, \
                             ta1_x_xxxyyz_yyy_0, \
                             ta1_x_xxxyyz_yyz_0, \
                             ta1_x_xxxyyz_yzz_0, \
                             ta1_x_xxxyyz_zzz_0, \
                             ta1_x_xxxyz_xxz_0,  \
                             ta1_x_xxxyz_xxz_1,  \
                             ta1_x_xxxyz_xzz_0,  \
                             ta1_x_xxxyz_xzz_1,  \
                             ta1_x_xxxyz_zzz_0,  \
                             ta1_x_xxxyz_zzz_1,  \
                             ta1_x_xxxz_xxz_0,   \
                             ta1_x_xxxz_xxz_1,   \
                             ta1_x_xxxz_xzz_0,   \
                             ta1_x_xxxz_xzz_1,   \
                             ta1_x_xxxz_zzz_0,   \
                             ta1_x_xxxz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyz_xxx_0[i] = ta1_x_xxxyy_xxx_0[i] * pa_z[i] - ta1_x_xxxyy_xxx_1[i] * pc_z[i];

        ta1_x_xxxyyz_xxy_0[i] = ta1_x_xxxyy_xxy_0[i] * pa_z[i] - ta1_x_xxxyy_xxy_1[i] * pc_z[i];

        ta1_x_xxxyyz_xxz_0[i] =
            ta1_x_xxxz_xxz_0[i] * fe_0 - ta1_x_xxxz_xxz_1[i] * fe_0 + ta1_x_xxxyz_xxz_0[i] * pa_y[i] - ta1_x_xxxyz_xxz_1[i] * pc_y[i];

        ta1_x_xxxyyz_xyy_0[i] = ta1_x_xxxyy_xyy_0[i] * pa_z[i] - ta1_x_xxxyy_xyy_1[i] * pc_z[i];

        ta1_x_xxxyyz_xyz_0[i] =
            ta1_x_xxxyy_xy_0[i] * fe_0 - ta1_x_xxxyy_xy_1[i] * fe_0 + ta1_x_xxxyy_xyz_0[i] * pa_z[i] - ta1_x_xxxyy_xyz_1[i] * pc_z[i];

        ta1_x_xxxyyz_xzz_0[i] =
            ta1_x_xxxz_xzz_0[i] * fe_0 - ta1_x_xxxz_xzz_1[i] * fe_0 + ta1_x_xxxyz_xzz_0[i] * pa_y[i] - ta1_x_xxxyz_xzz_1[i] * pc_y[i];

        ta1_x_xxxyyz_yyy_0[i] = ta1_x_xxxyy_yyy_0[i] * pa_z[i] - ta1_x_xxxyy_yyy_1[i] * pc_z[i];

        ta1_x_xxxyyz_yyz_0[i] =
            ta1_x_xxxyy_yy_0[i] * fe_0 - ta1_x_xxxyy_yy_1[i] * fe_0 + ta1_x_xxxyy_yyz_0[i] * pa_z[i] - ta1_x_xxxyy_yyz_1[i] * pc_z[i];

        ta1_x_xxxyyz_yzz_0[i] =
            2.0 * ta1_x_xxxyy_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxyy_yz_1[i] * fe_0 + ta1_x_xxxyy_yzz_0[i] * pa_z[i] - ta1_x_xxxyy_yzz_1[i] * pc_z[i];

        ta1_x_xxxyyz_zzz_0[i] =
            ta1_x_xxxz_zzz_0[i] * fe_0 - ta1_x_xxxz_zzz_1[i] * fe_0 + ta1_x_xxxyz_zzz_0[i] * pa_y[i] - ta1_x_xxxyz_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : IF

    auto ta1_x_xxxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 80);

    auto ta1_x_xxxyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 81);

    auto ta1_x_xxxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 82);

    auto ta1_x_xxxyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 83);

    auto ta1_x_xxxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 84);

    auto ta1_x_xxxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 85);

    auto ta1_x_xxxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 86);

    auto ta1_x_xxxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 87);

    auto ta1_x_xxxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 88);

    auto ta1_x_xxxyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 89);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_xxxyzz_xxx_0, \
                             ta1_x_xxxyzz_xxy_0, \
                             ta1_x_xxxyzz_xxz_0, \
                             ta1_x_xxxyzz_xyy_0, \
                             ta1_x_xxxyzz_xyz_0, \
                             ta1_x_xxxyzz_xzz_0, \
                             ta1_x_xxxyzz_yyy_0, \
                             ta1_x_xxxyzz_yyz_0, \
                             ta1_x_xxxyzz_yzz_0, \
                             ta1_x_xxxyzz_zzz_0, \
                             ta1_x_xxxzz_xx_0,   \
                             ta1_x_xxxzz_xx_1,   \
                             ta1_x_xxxzz_xxx_0,  \
                             ta1_x_xxxzz_xxx_1,  \
                             ta1_x_xxxzz_xxy_0,  \
                             ta1_x_xxxzz_xxy_1,  \
                             ta1_x_xxxzz_xxz_0,  \
                             ta1_x_xxxzz_xxz_1,  \
                             ta1_x_xxxzz_xy_0,   \
                             ta1_x_xxxzz_xy_1,   \
                             ta1_x_xxxzz_xyy_0,  \
                             ta1_x_xxxzz_xyy_1,  \
                             ta1_x_xxxzz_xyz_0,  \
                             ta1_x_xxxzz_xyz_1,  \
                             ta1_x_xxxzz_xz_0,   \
                             ta1_x_xxxzz_xz_1,   \
                             ta1_x_xxxzz_xzz_0,  \
                             ta1_x_xxxzz_xzz_1,  \
                             ta1_x_xxxzz_yy_0,   \
                             ta1_x_xxxzz_yy_1,   \
                             ta1_x_xxxzz_yyy_0,  \
                             ta1_x_xxxzz_yyy_1,  \
                             ta1_x_xxxzz_yyz_0,  \
                             ta1_x_xxxzz_yyz_1,  \
                             ta1_x_xxxzz_yz_0,   \
                             ta1_x_xxxzz_yz_1,   \
                             ta1_x_xxxzz_yzz_0,  \
                             ta1_x_xxxzz_yzz_1,  \
                             ta1_x_xxxzz_zz_0,   \
                             ta1_x_xxxzz_zz_1,   \
                             ta1_x_xxxzz_zzz_0,  \
                             ta1_x_xxxzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyzz_xxx_0[i] = ta1_x_xxxzz_xxx_0[i] * pa_y[i] - ta1_x_xxxzz_xxx_1[i] * pc_y[i];

        ta1_x_xxxyzz_xxy_0[i] =
            ta1_x_xxxzz_xx_0[i] * fe_0 - ta1_x_xxxzz_xx_1[i] * fe_0 + ta1_x_xxxzz_xxy_0[i] * pa_y[i] - ta1_x_xxxzz_xxy_1[i] * pc_y[i];

        ta1_x_xxxyzz_xxz_0[i] = ta1_x_xxxzz_xxz_0[i] * pa_y[i] - ta1_x_xxxzz_xxz_1[i] * pc_y[i];

        ta1_x_xxxyzz_xyy_0[i] =
            2.0 * ta1_x_xxxzz_xy_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_xy_1[i] * fe_0 + ta1_x_xxxzz_xyy_0[i] * pa_y[i] - ta1_x_xxxzz_xyy_1[i] * pc_y[i];

        ta1_x_xxxyzz_xyz_0[i] =
            ta1_x_xxxzz_xz_0[i] * fe_0 - ta1_x_xxxzz_xz_1[i] * fe_0 + ta1_x_xxxzz_xyz_0[i] * pa_y[i] - ta1_x_xxxzz_xyz_1[i] * pc_y[i];

        ta1_x_xxxyzz_xzz_0[i] = ta1_x_xxxzz_xzz_0[i] * pa_y[i] - ta1_x_xxxzz_xzz_1[i] * pc_y[i];

        ta1_x_xxxyzz_yyy_0[i] =
            3.0 * ta1_x_xxxzz_yy_0[i] * fe_0 - 3.0 * ta1_x_xxxzz_yy_1[i] * fe_0 + ta1_x_xxxzz_yyy_0[i] * pa_y[i] - ta1_x_xxxzz_yyy_1[i] * pc_y[i];

        ta1_x_xxxyzz_yyz_0[i] =
            2.0 * ta1_x_xxxzz_yz_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_yz_1[i] * fe_0 + ta1_x_xxxzz_yyz_0[i] * pa_y[i] - ta1_x_xxxzz_yyz_1[i] * pc_y[i];

        ta1_x_xxxyzz_yzz_0[i] =
            ta1_x_xxxzz_zz_0[i] * fe_0 - ta1_x_xxxzz_zz_1[i] * fe_0 + ta1_x_xxxzz_yzz_0[i] * pa_y[i] - ta1_x_xxxzz_yzz_1[i] * pc_y[i];

        ta1_x_xxxyzz_zzz_0[i] = ta1_x_xxxzz_zzz_0[i] * pa_y[i] - ta1_x_xxxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : IF

    auto ta1_x_xxxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 90);

    auto ta1_x_xxxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 91);

    auto ta1_x_xxxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 92);

    auto ta1_x_xxxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 93);

    auto ta1_x_xxxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 94);

    auto ta1_x_xxxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 95);

    auto ta1_x_xxxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 96);

    auto ta1_x_xxxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 97);

    auto ta1_x_xxxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 98);

    auto ta1_x_xxxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 99);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xxxz_xxx_0,   \
                             ta1_x_xxxz_xxx_1,   \
                             ta1_x_xxxz_xxy_0,   \
                             ta1_x_xxxz_xxy_1,   \
                             ta1_x_xxxz_xxz_0,   \
                             ta1_x_xxxz_xxz_1,   \
                             ta1_x_xxxz_xyy_0,   \
                             ta1_x_xxxz_xyy_1,   \
                             ta1_x_xxxz_xyz_0,   \
                             ta1_x_xxxz_xyz_1,   \
                             ta1_x_xxxz_xzz_0,   \
                             ta1_x_xxxz_xzz_1,   \
                             ta1_x_xxxz_yyy_0,   \
                             ta1_x_xxxz_yyy_1,   \
                             ta1_x_xxxzz_xx_0,   \
                             ta1_x_xxxzz_xx_1,   \
                             ta1_x_xxxzz_xxx_0,  \
                             ta1_x_xxxzz_xxx_1,  \
                             ta1_x_xxxzz_xxy_0,  \
                             ta1_x_xxxzz_xxy_1,  \
                             ta1_x_xxxzz_xxz_0,  \
                             ta1_x_xxxzz_xxz_1,  \
                             ta1_x_xxxzz_xy_0,   \
                             ta1_x_xxxzz_xy_1,   \
                             ta1_x_xxxzz_xyy_0,  \
                             ta1_x_xxxzz_xyy_1,  \
                             ta1_x_xxxzz_xyz_0,  \
                             ta1_x_xxxzz_xyz_1,  \
                             ta1_x_xxxzz_xz_0,   \
                             ta1_x_xxxzz_xz_1,   \
                             ta1_x_xxxzz_xzz_0,  \
                             ta1_x_xxxzz_xzz_1,  \
                             ta1_x_xxxzz_yyy_0,  \
                             ta1_x_xxxzz_yyy_1,  \
                             ta1_x_xxxzzz_xxx_0, \
                             ta1_x_xxxzzz_xxy_0, \
                             ta1_x_xxxzzz_xxz_0, \
                             ta1_x_xxxzzz_xyy_0, \
                             ta1_x_xxxzzz_xyz_0, \
                             ta1_x_xxxzzz_xzz_0, \
                             ta1_x_xxxzzz_yyy_0, \
                             ta1_x_xxxzzz_yyz_0, \
                             ta1_x_xxxzzz_yzz_0, \
                             ta1_x_xxxzzz_zzz_0, \
                             ta1_x_xxzzz_yyz_0,  \
                             ta1_x_xxzzz_yyz_1,  \
                             ta1_x_xxzzz_yzz_0,  \
                             ta1_x_xxzzz_yzz_1,  \
                             ta1_x_xxzzz_zzz_0,  \
                             ta1_x_xxzzz_zzz_1,  \
                             ta1_x_xzzz_yyz_0,   \
                             ta1_x_xzzz_yyz_1,   \
                             ta1_x_xzzz_yzz_0,   \
                             ta1_x_xzzz_yzz_1,   \
                             ta1_x_xzzz_zzz_0,   \
                             ta1_x_xzzz_zzz_1,   \
                             ta_xxzzz_yyz_1,     \
                             ta_xxzzz_yzz_1,     \
                             ta_xxzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzzz_xxx_0[i] =
            2.0 * ta1_x_xxxz_xxx_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxx_1[i] * fe_0 + ta1_x_xxxzz_xxx_0[i] * pa_z[i] - ta1_x_xxxzz_xxx_1[i] * pc_z[i];

        ta1_x_xxxzzz_xxy_0[i] =
            2.0 * ta1_x_xxxz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxy_1[i] * fe_0 + ta1_x_xxxzz_xxy_0[i] * pa_z[i] - ta1_x_xxxzz_xxy_1[i] * pc_z[i];

        ta1_x_xxxzzz_xxz_0[i] = 2.0 * ta1_x_xxxz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xxz_1[i] * fe_0 + ta1_x_xxxzz_xx_0[i] * fe_0 -
                                ta1_x_xxxzz_xx_1[i] * fe_0 + ta1_x_xxxzz_xxz_0[i] * pa_z[i] - ta1_x_xxxzz_xxz_1[i] * pc_z[i];

        ta1_x_xxxzzz_xyy_0[i] =
            2.0 * ta1_x_xxxz_xyy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xyy_1[i] * fe_0 + ta1_x_xxxzz_xyy_0[i] * pa_z[i] - ta1_x_xxxzz_xyy_1[i] * pc_z[i];

        ta1_x_xxxzzz_xyz_0[i] = 2.0 * ta1_x_xxxz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xyz_1[i] * fe_0 + ta1_x_xxxzz_xy_0[i] * fe_0 -
                                ta1_x_xxxzz_xy_1[i] * fe_0 + ta1_x_xxxzz_xyz_0[i] * pa_z[i] - ta1_x_xxxzz_xyz_1[i] * pc_z[i];

        ta1_x_xxxzzz_xzz_0[i] = 2.0 * ta1_x_xxxz_xzz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xzz_1[i] * fe_0 + 2.0 * ta1_x_xxxzz_xz_0[i] * fe_0 -
                                2.0 * ta1_x_xxxzz_xz_1[i] * fe_0 + ta1_x_xxxzz_xzz_0[i] * pa_z[i] - ta1_x_xxxzz_xzz_1[i] * pc_z[i];

        ta1_x_xxxzzz_yyy_0[i] =
            2.0 * ta1_x_xxxz_yyy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_yyy_1[i] * fe_0 + ta1_x_xxxzz_yyy_0[i] * pa_z[i] - ta1_x_xxxzz_yyy_1[i] * pc_z[i];

        ta1_x_xxxzzz_yyz_0[i] = 2.0 * ta1_x_xzzz_yyz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_yyz_1[i] * fe_0 + ta_xxzzz_yyz_1[i] +
                                ta1_x_xxzzz_yyz_0[i] * pa_x[i] - ta1_x_xxzzz_yyz_1[i] * pc_x[i];

        ta1_x_xxxzzz_yzz_0[i] = 2.0 * ta1_x_xzzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_yzz_1[i] * fe_0 + ta_xxzzz_yzz_1[i] +
                                ta1_x_xxzzz_yzz_0[i] * pa_x[i] - ta1_x_xxzzz_yzz_1[i] * pc_x[i];

        ta1_x_xxxzzz_zzz_0[i] = 2.0 * ta1_x_xzzz_zzz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_zzz_1[i] * fe_0 + ta_xxzzz_zzz_1[i] +
                                ta1_x_xxzzz_zzz_0[i] * pa_x[i] - ta1_x_xxzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : IF

    auto ta1_x_xxyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 100);

    auto ta1_x_xxyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 101);

    auto ta1_x_xxyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 102);

    auto ta1_x_xxyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 103);

    auto ta1_x_xxyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 104);

    auto ta1_x_xxyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 105);

    auto ta1_x_xxyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 106);

    auto ta1_x_xxyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 107);

    auto ta1_x_xxyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 108);

    auto ta1_x_xxyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 109);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xxyy_xxx_0,   \
                             ta1_x_xxyy_xxx_1,   \
                             ta1_x_xxyy_xxy_0,   \
                             ta1_x_xxyy_xxy_1,   \
                             ta1_x_xxyy_xxz_0,   \
                             ta1_x_xxyy_xxz_1,   \
                             ta1_x_xxyy_xyy_0,   \
                             ta1_x_xxyy_xyy_1,   \
                             ta1_x_xxyy_xyz_0,   \
                             ta1_x_xxyy_xyz_1,   \
                             ta1_x_xxyy_xzz_0,   \
                             ta1_x_xxyy_xzz_1,   \
                             ta1_x_xxyy_zzz_0,   \
                             ta1_x_xxyy_zzz_1,   \
                             ta1_x_xxyyy_xx_0,   \
                             ta1_x_xxyyy_xx_1,   \
                             ta1_x_xxyyy_xxx_0,  \
                             ta1_x_xxyyy_xxx_1,  \
                             ta1_x_xxyyy_xxy_0,  \
                             ta1_x_xxyyy_xxy_1,  \
                             ta1_x_xxyyy_xxz_0,  \
                             ta1_x_xxyyy_xxz_1,  \
                             ta1_x_xxyyy_xy_0,   \
                             ta1_x_xxyyy_xy_1,   \
                             ta1_x_xxyyy_xyy_0,  \
                             ta1_x_xxyyy_xyy_1,  \
                             ta1_x_xxyyy_xyz_0,  \
                             ta1_x_xxyyy_xyz_1,  \
                             ta1_x_xxyyy_xz_0,   \
                             ta1_x_xxyyy_xz_1,   \
                             ta1_x_xxyyy_xzz_0,  \
                             ta1_x_xxyyy_xzz_1,  \
                             ta1_x_xxyyy_zzz_0,  \
                             ta1_x_xxyyy_zzz_1,  \
                             ta1_x_xxyyyy_xxx_0, \
                             ta1_x_xxyyyy_xxy_0, \
                             ta1_x_xxyyyy_xxz_0, \
                             ta1_x_xxyyyy_xyy_0, \
                             ta1_x_xxyyyy_xyz_0, \
                             ta1_x_xxyyyy_xzz_0, \
                             ta1_x_xxyyyy_yyy_0, \
                             ta1_x_xxyyyy_yyz_0, \
                             ta1_x_xxyyyy_yzz_0, \
                             ta1_x_xxyyyy_zzz_0, \
                             ta1_x_xyyyy_yyy_0,  \
                             ta1_x_xyyyy_yyy_1,  \
                             ta1_x_xyyyy_yyz_0,  \
                             ta1_x_xyyyy_yyz_1,  \
                             ta1_x_xyyyy_yzz_0,  \
                             ta1_x_xyyyy_yzz_1,  \
                             ta1_x_yyyy_yyy_0,   \
                             ta1_x_yyyy_yyy_1,   \
                             ta1_x_yyyy_yyz_0,   \
                             ta1_x_yyyy_yyz_1,   \
                             ta1_x_yyyy_yzz_0,   \
                             ta1_x_yyyy_yzz_1,   \
                             ta_xyyyy_yyy_1,     \
                             ta_xyyyy_yyz_1,     \
                             ta_xyyyy_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyy_xxx_0[i] =
            3.0 * ta1_x_xxyy_xxx_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxx_1[i] * fe_0 + ta1_x_xxyyy_xxx_0[i] * pa_y[i] - ta1_x_xxyyy_xxx_1[i] * pc_y[i];

        ta1_x_xxyyyy_xxy_0[i] = 3.0 * ta1_x_xxyy_xxy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxy_1[i] * fe_0 + ta1_x_xxyyy_xx_0[i] * fe_0 -
                                ta1_x_xxyyy_xx_1[i] * fe_0 + ta1_x_xxyyy_xxy_0[i] * pa_y[i] - ta1_x_xxyyy_xxy_1[i] * pc_y[i];

        ta1_x_xxyyyy_xxz_0[i] =
            3.0 * ta1_x_xxyy_xxz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xxz_1[i] * fe_0 + ta1_x_xxyyy_xxz_0[i] * pa_y[i] - ta1_x_xxyyy_xxz_1[i] * pc_y[i];

        ta1_x_xxyyyy_xyy_0[i] = 3.0 * ta1_x_xxyy_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xyy_1[i] * fe_0 + 2.0 * ta1_x_xxyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_x_xxyyy_xy_1[i] * fe_0 + ta1_x_xxyyy_xyy_0[i] * pa_y[i] - ta1_x_xxyyy_xyy_1[i] * pc_y[i];

        ta1_x_xxyyyy_xyz_0[i] = 3.0 * ta1_x_xxyy_xyz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xyz_1[i] * fe_0 + ta1_x_xxyyy_xz_0[i] * fe_0 -
                                ta1_x_xxyyy_xz_1[i] * fe_0 + ta1_x_xxyyy_xyz_0[i] * pa_y[i] - ta1_x_xxyyy_xyz_1[i] * pc_y[i];

        ta1_x_xxyyyy_xzz_0[i] =
            3.0 * ta1_x_xxyy_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xzz_1[i] * fe_0 + ta1_x_xxyyy_xzz_0[i] * pa_y[i] - ta1_x_xxyyy_xzz_1[i] * pc_y[i];

        ta1_x_xxyyyy_yyy_0[i] = ta1_x_yyyy_yyy_0[i] * fe_0 - ta1_x_yyyy_yyy_1[i] * fe_0 + ta_xyyyy_yyy_1[i] + ta1_x_xyyyy_yyy_0[i] * pa_x[i] -
                                ta1_x_xyyyy_yyy_1[i] * pc_x[i];

        ta1_x_xxyyyy_yyz_0[i] = ta1_x_yyyy_yyz_0[i] * fe_0 - ta1_x_yyyy_yyz_1[i] * fe_0 + ta_xyyyy_yyz_1[i] + ta1_x_xyyyy_yyz_0[i] * pa_x[i] -
                                ta1_x_xyyyy_yyz_1[i] * pc_x[i];

        ta1_x_xxyyyy_yzz_0[i] = ta1_x_yyyy_yzz_0[i] * fe_0 - ta1_x_yyyy_yzz_1[i] * fe_0 + ta_xyyyy_yzz_1[i] + ta1_x_xyyyy_yzz_0[i] * pa_x[i] -
                                ta1_x_xyyyy_yzz_1[i] * pc_x[i];

        ta1_x_xxyyyy_zzz_0[i] =
            3.0 * ta1_x_xxyy_zzz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_zzz_1[i] * fe_0 + ta1_x_xxyyy_zzz_0[i] * pa_y[i] - ta1_x_xxyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 110-120 components of targeted buffer : IF

    auto ta1_x_xxyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 110);

    auto ta1_x_xxyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 111);

    auto ta1_x_xxyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 112);

    auto ta1_x_xxyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 113);

    auto ta1_x_xxyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 114);

    auto ta1_x_xxyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 115);

    auto ta1_x_xxyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 116);

    auto ta1_x_xxyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 117);

    auto ta1_x_xxyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 118);

    auto ta1_x_xxyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 119);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxyyy_xxx_0,  \
                             ta1_x_xxyyy_xxx_1,  \
                             ta1_x_xxyyy_xxy_0,  \
                             ta1_x_xxyyy_xxy_1,  \
                             ta1_x_xxyyy_xy_0,   \
                             ta1_x_xxyyy_xy_1,   \
                             ta1_x_xxyyy_xyy_0,  \
                             ta1_x_xxyyy_xyy_1,  \
                             ta1_x_xxyyy_xyz_0,  \
                             ta1_x_xxyyy_xyz_1,  \
                             ta1_x_xxyyy_yy_0,   \
                             ta1_x_xxyyy_yy_1,   \
                             ta1_x_xxyyy_yyy_0,  \
                             ta1_x_xxyyy_yyy_1,  \
                             ta1_x_xxyyy_yyz_0,  \
                             ta1_x_xxyyy_yyz_1,  \
                             ta1_x_xxyyy_yz_0,   \
                             ta1_x_xxyyy_yz_1,   \
                             ta1_x_xxyyy_yzz_0,  \
                             ta1_x_xxyyy_yzz_1,  \
                             ta1_x_xxyyyz_xxx_0, \
                             ta1_x_xxyyyz_xxy_0, \
                             ta1_x_xxyyyz_xxz_0, \
                             ta1_x_xxyyyz_xyy_0, \
                             ta1_x_xxyyyz_xyz_0, \
                             ta1_x_xxyyyz_xzz_0, \
                             ta1_x_xxyyyz_yyy_0, \
                             ta1_x_xxyyyz_yyz_0, \
                             ta1_x_xxyyyz_yzz_0, \
                             ta1_x_xxyyyz_zzz_0, \
                             ta1_x_xxyyz_xxz_0,  \
                             ta1_x_xxyyz_xxz_1,  \
                             ta1_x_xxyyz_xzz_0,  \
                             ta1_x_xxyyz_xzz_1,  \
                             ta1_x_xxyyz_zzz_0,  \
                             ta1_x_xxyyz_zzz_1,  \
                             ta1_x_xxyz_xxz_0,   \
                             ta1_x_xxyz_xxz_1,   \
                             ta1_x_xxyz_xzz_0,   \
                             ta1_x_xxyz_xzz_1,   \
                             ta1_x_xxyz_zzz_0,   \
                             ta1_x_xxyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyz_xxx_0[i] = ta1_x_xxyyy_xxx_0[i] * pa_z[i] - ta1_x_xxyyy_xxx_1[i] * pc_z[i];

        ta1_x_xxyyyz_xxy_0[i] = ta1_x_xxyyy_xxy_0[i] * pa_z[i] - ta1_x_xxyyy_xxy_1[i] * pc_z[i];

        ta1_x_xxyyyz_xxz_0[i] =
            2.0 * ta1_x_xxyz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_xxz_1[i] * fe_0 + ta1_x_xxyyz_xxz_0[i] * pa_y[i] - ta1_x_xxyyz_xxz_1[i] * pc_y[i];

        ta1_x_xxyyyz_xyy_0[i] = ta1_x_xxyyy_xyy_0[i] * pa_z[i] - ta1_x_xxyyy_xyy_1[i] * pc_z[i];

        ta1_x_xxyyyz_xyz_0[i] =
            ta1_x_xxyyy_xy_0[i] * fe_0 - ta1_x_xxyyy_xy_1[i] * fe_0 + ta1_x_xxyyy_xyz_0[i] * pa_z[i] - ta1_x_xxyyy_xyz_1[i] * pc_z[i];

        ta1_x_xxyyyz_xzz_0[i] =
            2.0 * ta1_x_xxyz_xzz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_xzz_1[i] * fe_0 + ta1_x_xxyyz_xzz_0[i] * pa_y[i] - ta1_x_xxyyz_xzz_1[i] * pc_y[i];

        ta1_x_xxyyyz_yyy_0[i] = ta1_x_xxyyy_yyy_0[i] * pa_z[i] - ta1_x_xxyyy_yyy_1[i] * pc_z[i];

        ta1_x_xxyyyz_yyz_0[i] =
            ta1_x_xxyyy_yy_0[i] * fe_0 - ta1_x_xxyyy_yy_1[i] * fe_0 + ta1_x_xxyyy_yyz_0[i] * pa_z[i] - ta1_x_xxyyy_yyz_1[i] * pc_z[i];

        ta1_x_xxyyyz_yzz_0[i] =
            2.0 * ta1_x_xxyyy_yz_0[i] * fe_0 - 2.0 * ta1_x_xxyyy_yz_1[i] * fe_0 + ta1_x_xxyyy_yzz_0[i] * pa_z[i] - ta1_x_xxyyy_yzz_1[i] * pc_z[i];

        ta1_x_xxyyyz_zzz_0[i] =
            2.0 * ta1_x_xxyz_zzz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_zzz_1[i] * fe_0 + ta1_x_xxyyz_zzz_0[i] * pa_y[i] - ta1_x_xxyyz_zzz_1[i] * pc_y[i];
    }

    // Set up 120-130 components of targeted buffer : IF

    auto ta1_x_xxyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 120);

    auto ta1_x_xxyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 121);

    auto ta1_x_xxyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 122);

    auto ta1_x_xxyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 123);

    auto ta1_x_xxyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 124);

    auto ta1_x_xxyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 125);

    auto ta1_x_xxyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 126);

    auto ta1_x_xxyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 127);

    auto ta1_x_xxyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 128);

    auto ta1_x_xxyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 129);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxyy_xxy_0,   \
                             ta1_x_xxyy_xxy_1,   \
                             ta1_x_xxyy_xyy_0,   \
                             ta1_x_xxyy_xyy_1,   \
                             ta1_x_xxyy_yyy_0,   \
                             ta1_x_xxyy_yyy_1,   \
                             ta1_x_xxyyz_xxy_0,  \
                             ta1_x_xxyyz_xxy_1,  \
                             ta1_x_xxyyz_xyy_0,  \
                             ta1_x_xxyyz_xyy_1,  \
                             ta1_x_xxyyz_yyy_0,  \
                             ta1_x_xxyyz_yyy_1,  \
                             ta1_x_xxyyzz_xxx_0, \
                             ta1_x_xxyyzz_xxy_0, \
                             ta1_x_xxyyzz_xxz_0, \
                             ta1_x_xxyyzz_xyy_0, \
                             ta1_x_xxyyzz_xyz_0, \
                             ta1_x_xxyyzz_xzz_0, \
                             ta1_x_xxyyzz_yyy_0, \
                             ta1_x_xxyyzz_yyz_0, \
                             ta1_x_xxyyzz_yzz_0, \
                             ta1_x_xxyyzz_zzz_0, \
                             ta1_x_xxyzz_xxx_0,  \
                             ta1_x_xxyzz_xxx_1,  \
                             ta1_x_xxyzz_xxz_0,  \
                             ta1_x_xxyzz_xxz_1,  \
                             ta1_x_xxyzz_xyz_0,  \
                             ta1_x_xxyzz_xyz_1,  \
                             ta1_x_xxyzz_xz_0,   \
                             ta1_x_xxyzz_xz_1,   \
                             ta1_x_xxyzz_xzz_0,  \
                             ta1_x_xxyzz_xzz_1,  \
                             ta1_x_xxyzz_zzz_0,  \
                             ta1_x_xxyzz_zzz_1,  \
                             ta1_x_xxzz_xxx_0,   \
                             ta1_x_xxzz_xxx_1,   \
                             ta1_x_xxzz_xxz_0,   \
                             ta1_x_xxzz_xxz_1,   \
                             ta1_x_xxzz_xyz_0,   \
                             ta1_x_xxzz_xyz_1,   \
                             ta1_x_xxzz_xzz_0,   \
                             ta1_x_xxzz_xzz_1,   \
                             ta1_x_xxzz_zzz_0,   \
                             ta1_x_xxzz_zzz_1,   \
                             ta1_x_xyyzz_yyz_0,  \
                             ta1_x_xyyzz_yyz_1,  \
                             ta1_x_xyyzz_yzz_0,  \
                             ta1_x_xyyzz_yzz_1,  \
                             ta1_x_yyzz_yyz_0,   \
                             ta1_x_yyzz_yyz_1,   \
                             ta1_x_yyzz_yzz_0,   \
                             ta1_x_yyzz_yzz_1,   \
                             ta_xyyzz_yyz_1,     \
                             ta_xyyzz_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyzz_xxx_0[i] =
            ta1_x_xxzz_xxx_0[i] * fe_0 - ta1_x_xxzz_xxx_1[i] * fe_0 + ta1_x_xxyzz_xxx_0[i] * pa_y[i] - ta1_x_xxyzz_xxx_1[i] * pc_y[i];

        ta1_x_xxyyzz_xxy_0[i] =
            ta1_x_xxyy_xxy_0[i] * fe_0 - ta1_x_xxyy_xxy_1[i] * fe_0 + ta1_x_xxyyz_xxy_0[i] * pa_z[i] - ta1_x_xxyyz_xxy_1[i] * pc_z[i];

        ta1_x_xxyyzz_xxz_0[i] =
            ta1_x_xxzz_xxz_0[i] * fe_0 - ta1_x_xxzz_xxz_1[i] * fe_0 + ta1_x_xxyzz_xxz_0[i] * pa_y[i] - ta1_x_xxyzz_xxz_1[i] * pc_y[i];

        ta1_x_xxyyzz_xyy_0[i] =
            ta1_x_xxyy_xyy_0[i] * fe_0 - ta1_x_xxyy_xyy_1[i] * fe_0 + ta1_x_xxyyz_xyy_0[i] * pa_z[i] - ta1_x_xxyyz_xyy_1[i] * pc_z[i];

        ta1_x_xxyyzz_xyz_0[i] = ta1_x_xxzz_xyz_0[i] * fe_0 - ta1_x_xxzz_xyz_1[i] * fe_0 + ta1_x_xxyzz_xz_0[i] * fe_0 - ta1_x_xxyzz_xz_1[i] * fe_0 +
                                ta1_x_xxyzz_xyz_0[i] * pa_y[i] - ta1_x_xxyzz_xyz_1[i] * pc_y[i];

        ta1_x_xxyyzz_xzz_0[i] =
            ta1_x_xxzz_xzz_0[i] * fe_0 - ta1_x_xxzz_xzz_1[i] * fe_0 + ta1_x_xxyzz_xzz_0[i] * pa_y[i] - ta1_x_xxyzz_xzz_1[i] * pc_y[i];

        ta1_x_xxyyzz_yyy_0[i] =
            ta1_x_xxyy_yyy_0[i] * fe_0 - ta1_x_xxyy_yyy_1[i] * fe_0 + ta1_x_xxyyz_yyy_0[i] * pa_z[i] - ta1_x_xxyyz_yyy_1[i] * pc_z[i];

        ta1_x_xxyyzz_yyz_0[i] = ta1_x_yyzz_yyz_0[i] * fe_0 - ta1_x_yyzz_yyz_1[i] * fe_0 + ta_xyyzz_yyz_1[i] + ta1_x_xyyzz_yyz_0[i] * pa_x[i] -
                                ta1_x_xyyzz_yyz_1[i] * pc_x[i];

        ta1_x_xxyyzz_yzz_0[i] = ta1_x_yyzz_yzz_0[i] * fe_0 - ta1_x_yyzz_yzz_1[i] * fe_0 + ta_xyyzz_yzz_1[i] + ta1_x_xyyzz_yzz_0[i] * pa_x[i] -
                                ta1_x_xyyzz_yzz_1[i] * pc_x[i];

        ta1_x_xxyyzz_zzz_0[i] =
            ta1_x_xxzz_zzz_0[i] * fe_0 - ta1_x_xxzz_zzz_1[i] * fe_0 + ta1_x_xxyzz_zzz_0[i] * pa_y[i] - ta1_x_xxyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 130-140 components of targeted buffer : IF

    auto ta1_x_xxyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 130);

    auto ta1_x_xxyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 131);

    auto ta1_x_xxyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 132);

    auto ta1_x_xxyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 133);

    auto ta1_x_xxyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 134);

    auto ta1_x_xxyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 135);

    auto ta1_x_xxyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 136);

    auto ta1_x_xxyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 137);

    auto ta1_x_xxyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 138);

    auto ta1_x_xxyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 139);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_xxyzzz_xxx_0, \
                             ta1_x_xxyzzz_xxy_0, \
                             ta1_x_xxyzzz_xxz_0, \
                             ta1_x_xxyzzz_xyy_0, \
                             ta1_x_xxyzzz_xyz_0, \
                             ta1_x_xxyzzz_xzz_0, \
                             ta1_x_xxyzzz_yyy_0, \
                             ta1_x_xxyzzz_yyz_0, \
                             ta1_x_xxyzzz_yzz_0, \
                             ta1_x_xxyzzz_zzz_0, \
                             ta1_x_xxzzz_xx_0,   \
                             ta1_x_xxzzz_xx_1,   \
                             ta1_x_xxzzz_xxx_0,  \
                             ta1_x_xxzzz_xxx_1,  \
                             ta1_x_xxzzz_xxy_0,  \
                             ta1_x_xxzzz_xxy_1,  \
                             ta1_x_xxzzz_xxz_0,  \
                             ta1_x_xxzzz_xxz_1,  \
                             ta1_x_xxzzz_xy_0,   \
                             ta1_x_xxzzz_xy_1,   \
                             ta1_x_xxzzz_xyy_0,  \
                             ta1_x_xxzzz_xyy_1,  \
                             ta1_x_xxzzz_xyz_0,  \
                             ta1_x_xxzzz_xyz_1,  \
                             ta1_x_xxzzz_xz_0,   \
                             ta1_x_xxzzz_xz_1,   \
                             ta1_x_xxzzz_xzz_0,  \
                             ta1_x_xxzzz_xzz_1,  \
                             ta1_x_xxzzz_yy_0,   \
                             ta1_x_xxzzz_yy_1,   \
                             ta1_x_xxzzz_yyy_0,  \
                             ta1_x_xxzzz_yyy_1,  \
                             ta1_x_xxzzz_yyz_0,  \
                             ta1_x_xxzzz_yyz_1,  \
                             ta1_x_xxzzz_yz_0,   \
                             ta1_x_xxzzz_yz_1,   \
                             ta1_x_xxzzz_yzz_0,  \
                             ta1_x_xxzzz_yzz_1,  \
                             ta1_x_xxzzz_zz_0,   \
                             ta1_x_xxzzz_zz_1,   \
                             ta1_x_xxzzz_zzz_0,  \
                             ta1_x_xxzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzzz_xxx_0[i] = ta1_x_xxzzz_xxx_0[i] * pa_y[i] - ta1_x_xxzzz_xxx_1[i] * pc_y[i];

        ta1_x_xxyzzz_xxy_0[i] =
            ta1_x_xxzzz_xx_0[i] * fe_0 - ta1_x_xxzzz_xx_1[i] * fe_0 + ta1_x_xxzzz_xxy_0[i] * pa_y[i] - ta1_x_xxzzz_xxy_1[i] * pc_y[i];

        ta1_x_xxyzzz_xxz_0[i] = ta1_x_xxzzz_xxz_0[i] * pa_y[i] - ta1_x_xxzzz_xxz_1[i] * pc_y[i];

        ta1_x_xxyzzz_xyy_0[i] =
            2.0 * ta1_x_xxzzz_xy_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_xy_1[i] * fe_0 + ta1_x_xxzzz_xyy_0[i] * pa_y[i] - ta1_x_xxzzz_xyy_1[i] * pc_y[i];

        ta1_x_xxyzzz_xyz_0[i] =
            ta1_x_xxzzz_xz_0[i] * fe_0 - ta1_x_xxzzz_xz_1[i] * fe_0 + ta1_x_xxzzz_xyz_0[i] * pa_y[i] - ta1_x_xxzzz_xyz_1[i] * pc_y[i];

        ta1_x_xxyzzz_xzz_0[i] = ta1_x_xxzzz_xzz_0[i] * pa_y[i] - ta1_x_xxzzz_xzz_1[i] * pc_y[i];

        ta1_x_xxyzzz_yyy_0[i] =
            3.0 * ta1_x_xxzzz_yy_0[i] * fe_0 - 3.0 * ta1_x_xxzzz_yy_1[i] * fe_0 + ta1_x_xxzzz_yyy_0[i] * pa_y[i] - ta1_x_xxzzz_yyy_1[i] * pc_y[i];

        ta1_x_xxyzzz_yyz_0[i] =
            2.0 * ta1_x_xxzzz_yz_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_yz_1[i] * fe_0 + ta1_x_xxzzz_yyz_0[i] * pa_y[i] - ta1_x_xxzzz_yyz_1[i] * pc_y[i];

        ta1_x_xxyzzz_yzz_0[i] =
            ta1_x_xxzzz_zz_0[i] * fe_0 - ta1_x_xxzzz_zz_1[i] * fe_0 + ta1_x_xxzzz_yzz_0[i] * pa_y[i] - ta1_x_xxzzz_yzz_1[i] * pc_y[i];

        ta1_x_xxyzzz_zzz_0[i] = ta1_x_xxzzz_zzz_0[i] * pa_y[i] - ta1_x_xxzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : IF

    auto ta1_x_xxzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 140);

    auto ta1_x_xxzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 141);

    auto ta1_x_xxzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 142);

    auto ta1_x_xxzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 143);

    auto ta1_x_xxzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 144);

    auto ta1_x_xxzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 145);

    auto ta1_x_xxzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 146);

    auto ta1_x_xxzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 147);

    auto ta1_x_xxzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 148);

    auto ta1_x_xxzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 149);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xxzz_xxx_0,   \
                             ta1_x_xxzz_xxx_1,   \
                             ta1_x_xxzz_xxy_0,   \
                             ta1_x_xxzz_xxy_1,   \
                             ta1_x_xxzz_xxz_0,   \
                             ta1_x_xxzz_xxz_1,   \
                             ta1_x_xxzz_xyy_0,   \
                             ta1_x_xxzz_xyy_1,   \
                             ta1_x_xxzz_xyz_0,   \
                             ta1_x_xxzz_xyz_1,   \
                             ta1_x_xxzz_xzz_0,   \
                             ta1_x_xxzz_xzz_1,   \
                             ta1_x_xxzz_yyy_0,   \
                             ta1_x_xxzz_yyy_1,   \
                             ta1_x_xxzzz_xx_0,   \
                             ta1_x_xxzzz_xx_1,   \
                             ta1_x_xxzzz_xxx_0,  \
                             ta1_x_xxzzz_xxx_1,  \
                             ta1_x_xxzzz_xxy_0,  \
                             ta1_x_xxzzz_xxy_1,  \
                             ta1_x_xxzzz_xxz_0,  \
                             ta1_x_xxzzz_xxz_1,  \
                             ta1_x_xxzzz_xy_0,   \
                             ta1_x_xxzzz_xy_1,   \
                             ta1_x_xxzzz_xyy_0,  \
                             ta1_x_xxzzz_xyy_1,  \
                             ta1_x_xxzzz_xyz_0,  \
                             ta1_x_xxzzz_xyz_1,  \
                             ta1_x_xxzzz_xz_0,   \
                             ta1_x_xxzzz_xz_1,   \
                             ta1_x_xxzzz_xzz_0,  \
                             ta1_x_xxzzz_xzz_1,  \
                             ta1_x_xxzzz_yyy_0,  \
                             ta1_x_xxzzz_yyy_1,  \
                             ta1_x_xxzzzz_xxx_0, \
                             ta1_x_xxzzzz_xxy_0, \
                             ta1_x_xxzzzz_xxz_0, \
                             ta1_x_xxzzzz_xyy_0, \
                             ta1_x_xxzzzz_xyz_0, \
                             ta1_x_xxzzzz_xzz_0, \
                             ta1_x_xxzzzz_yyy_0, \
                             ta1_x_xxzzzz_yyz_0, \
                             ta1_x_xxzzzz_yzz_0, \
                             ta1_x_xxzzzz_zzz_0, \
                             ta1_x_xzzzz_yyz_0,  \
                             ta1_x_xzzzz_yyz_1,  \
                             ta1_x_xzzzz_yzz_0,  \
                             ta1_x_xzzzz_yzz_1,  \
                             ta1_x_xzzzz_zzz_0,  \
                             ta1_x_xzzzz_zzz_1,  \
                             ta1_x_zzzz_yyz_0,   \
                             ta1_x_zzzz_yyz_1,   \
                             ta1_x_zzzz_yzz_0,   \
                             ta1_x_zzzz_yzz_1,   \
                             ta1_x_zzzz_zzz_0,   \
                             ta1_x_zzzz_zzz_1,   \
                             ta_xzzzz_yyz_1,     \
                             ta_xzzzz_yzz_1,     \
                             ta_xzzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzzz_xxx_0[i] =
            3.0 * ta1_x_xxzz_xxx_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxx_1[i] * fe_0 + ta1_x_xxzzz_xxx_0[i] * pa_z[i] - ta1_x_xxzzz_xxx_1[i] * pc_z[i];

        ta1_x_xxzzzz_xxy_0[i] =
            3.0 * ta1_x_xxzz_xxy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxy_1[i] * fe_0 + ta1_x_xxzzz_xxy_0[i] * pa_z[i] - ta1_x_xxzzz_xxy_1[i] * pc_z[i];

        ta1_x_xxzzzz_xxz_0[i] = 3.0 * ta1_x_xxzz_xxz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xxz_1[i] * fe_0 + ta1_x_xxzzz_xx_0[i] * fe_0 -
                                ta1_x_xxzzz_xx_1[i] * fe_0 + ta1_x_xxzzz_xxz_0[i] * pa_z[i] - ta1_x_xxzzz_xxz_1[i] * pc_z[i];

        ta1_x_xxzzzz_xyy_0[i] =
            3.0 * ta1_x_xxzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xyy_1[i] * fe_0 + ta1_x_xxzzz_xyy_0[i] * pa_z[i] - ta1_x_xxzzz_xyy_1[i] * pc_z[i];

        ta1_x_xxzzzz_xyz_0[i] = 3.0 * ta1_x_xxzz_xyz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xyz_1[i] * fe_0 + ta1_x_xxzzz_xy_0[i] * fe_0 -
                                ta1_x_xxzzz_xy_1[i] * fe_0 + ta1_x_xxzzz_xyz_0[i] * pa_z[i] - ta1_x_xxzzz_xyz_1[i] * pc_z[i];

        ta1_x_xxzzzz_xzz_0[i] = 3.0 * ta1_x_xxzz_xzz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xzz_1[i] * fe_0 + 2.0 * ta1_x_xxzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_x_xxzzz_xz_1[i] * fe_0 + ta1_x_xxzzz_xzz_0[i] * pa_z[i] - ta1_x_xxzzz_xzz_1[i] * pc_z[i];

        ta1_x_xxzzzz_yyy_0[i] =
            3.0 * ta1_x_xxzz_yyy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yyy_1[i] * fe_0 + ta1_x_xxzzz_yyy_0[i] * pa_z[i] - ta1_x_xxzzz_yyy_1[i] * pc_z[i];

        ta1_x_xxzzzz_yyz_0[i] = ta1_x_zzzz_yyz_0[i] * fe_0 - ta1_x_zzzz_yyz_1[i] * fe_0 + ta_xzzzz_yyz_1[i] + ta1_x_xzzzz_yyz_0[i] * pa_x[i] -
                                ta1_x_xzzzz_yyz_1[i] * pc_x[i];

        ta1_x_xxzzzz_yzz_0[i] = ta1_x_zzzz_yzz_0[i] * fe_0 - ta1_x_zzzz_yzz_1[i] * fe_0 + ta_xzzzz_yzz_1[i] + ta1_x_xzzzz_yzz_0[i] * pa_x[i] -
                                ta1_x_xzzzz_yzz_1[i] * pc_x[i];

        ta1_x_xxzzzz_zzz_0[i] = ta1_x_zzzz_zzz_0[i] * fe_0 - ta1_x_zzzz_zzz_1[i] * fe_0 + ta_xzzzz_zzz_1[i] + ta1_x_xzzzz_zzz_0[i] * pa_x[i] -
                                ta1_x_xzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : IF

    auto ta1_x_xyyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 150);

    auto ta1_x_xyyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 151);

    auto ta1_x_xyyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 152);

    auto ta1_x_xyyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 153);

    auto ta1_x_xyyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 154);

    auto ta1_x_xyyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 155);

    auto ta1_x_xyyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 156);

    auto ta1_x_xyyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 157);

    auto ta1_x_xyyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 158);

    auto ta1_x_xyyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 159);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xyyy_xxx_0,   \
                             ta1_x_xyyy_xxx_1,   \
                             ta1_x_xyyy_xxz_0,   \
                             ta1_x_xyyy_xxz_1,   \
                             ta1_x_xyyy_xzz_0,   \
                             ta1_x_xyyy_xzz_1,   \
                             ta1_x_xyyyy_xxx_0,  \
                             ta1_x_xyyyy_xxx_1,  \
                             ta1_x_xyyyy_xxz_0,  \
                             ta1_x_xyyyy_xxz_1,  \
                             ta1_x_xyyyy_xzz_0,  \
                             ta1_x_xyyyy_xzz_1,  \
                             ta1_x_xyyyyy_xxx_0, \
                             ta1_x_xyyyyy_xxy_0, \
                             ta1_x_xyyyyy_xxz_0, \
                             ta1_x_xyyyyy_xyy_0, \
                             ta1_x_xyyyyy_xyz_0, \
                             ta1_x_xyyyyy_xzz_0, \
                             ta1_x_xyyyyy_yyy_0, \
                             ta1_x_xyyyyy_yyz_0, \
                             ta1_x_xyyyyy_yzz_0, \
                             ta1_x_xyyyyy_zzz_0, \
                             ta1_x_yyyyy_xxy_0,  \
                             ta1_x_yyyyy_xxy_1,  \
                             ta1_x_yyyyy_xy_0,   \
                             ta1_x_yyyyy_xy_1,   \
                             ta1_x_yyyyy_xyy_0,  \
                             ta1_x_yyyyy_xyy_1,  \
                             ta1_x_yyyyy_xyz_0,  \
                             ta1_x_yyyyy_xyz_1,  \
                             ta1_x_yyyyy_yy_0,   \
                             ta1_x_yyyyy_yy_1,   \
                             ta1_x_yyyyy_yyy_0,  \
                             ta1_x_yyyyy_yyy_1,  \
                             ta1_x_yyyyy_yyz_0,  \
                             ta1_x_yyyyy_yyz_1,  \
                             ta1_x_yyyyy_yz_0,   \
                             ta1_x_yyyyy_yz_1,   \
                             ta1_x_yyyyy_yzz_0,  \
                             ta1_x_yyyyy_yzz_1,  \
                             ta1_x_yyyyy_zzz_0,  \
                             ta1_x_yyyyy_zzz_1,  \
                             ta_yyyyy_xxy_1,     \
                             ta_yyyyy_xyy_1,     \
                             ta_yyyyy_xyz_1,     \
                             ta_yyyyy_yyy_1,     \
                             ta_yyyyy_yyz_1,     \
                             ta_yyyyy_yzz_1,     \
                             ta_yyyyy_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyyy_xxx_0[i] =
            4.0 * ta1_x_xyyy_xxx_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xxx_1[i] * fe_0 + ta1_x_xyyyy_xxx_0[i] * pa_y[i] - ta1_x_xyyyy_xxx_1[i] * pc_y[i];

        ta1_x_xyyyyy_xxy_0[i] = 2.0 * ta1_x_yyyyy_xy_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_xy_1[i] * fe_0 + ta_yyyyy_xxy_1[i] +
                                ta1_x_yyyyy_xxy_0[i] * pa_x[i] - ta1_x_yyyyy_xxy_1[i] * pc_x[i];

        ta1_x_xyyyyy_xxz_0[i] =
            4.0 * ta1_x_xyyy_xxz_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xxz_1[i] * fe_0 + ta1_x_xyyyy_xxz_0[i] * pa_y[i] - ta1_x_xyyyy_xxz_1[i] * pc_y[i];

        ta1_x_xyyyyy_xyy_0[i] = ta1_x_yyyyy_yy_0[i] * fe_0 - ta1_x_yyyyy_yy_1[i] * fe_0 + ta_yyyyy_xyy_1[i] + ta1_x_yyyyy_xyy_0[i] * pa_x[i] -
                                ta1_x_yyyyy_xyy_1[i] * pc_x[i];

        ta1_x_xyyyyy_xyz_0[i] = ta1_x_yyyyy_yz_0[i] * fe_0 - ta1_x_yyyyy_yz_1[i] * fe_0 + ta_yyyyy_xyz_1[i] + ta1_x_yyyyy_xyz_0[i] * pa_x[i] -
                                ta1_x_yyyyy_xyz_1[i] * pc_x[i];

        ta1_x_xyyyyy_xzz_0[i] =
            4.0 * ta1_x_xyyy_xzz_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xzz_1[i] * fe_0 + ta1_x_xyyyy_xzz_0[i] * pa_y[i] - ta1_x_xyyyy_xzz_1[i] * pc_y[i];

        ta1_x_xyyyyy_yyy_0[i] = ta_yyyyy_yyy_1[i] + ta1_x_yyyyy_yyy_0[i] * pa_x[i] - ta1_x_yyyyy_yyy_1[i] * pc_x[i];

        ta1_x_xyyyyy_yyz_0[i] = ta_yyyyy_yyz_1[i] + ta1_x_yyyyy_yyz_0[i] * pa_x[i] - ta1_x_yyyyy_yyz_1[i] * pc_x[i];

        ta1_x_xyyyyy_yzz_0[i] = ta_yyyyy_yzz_1[i] + ta1_x_yyyyy_yzz_0[i] * pa_x[i] - ta1_x_yyyyy_yzz_1[i] * pc_x[i];

        ta1_x_xyyyyy_zzz_0[i] = ta_yyyyy_zzz_1[i] + ta1_x_yyyyy_zzz_0[i] * pa_x[i] - ta1_x_yyyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 160-170 components of targeted buffer : IF

    auto ta1_x_xyyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 160);

    auto ta1_x_xyyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 161);

    auto ta1_x_xyyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 162);

    auto ta1_x_xyyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 163);

    auto ta1_x_xyyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 164);

    auto ta1_x_xyyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 165);

    auto ta1_x_xyyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 166);

    auto ta1_x_xyyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 167);

    auto ta1_x_xyyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 168);

    auto ta1_x_xyyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 169);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xyyyy_xxx_0,  \
                             ta1_x_xyyyy_xxx_1,  \
                             ta1_x_xyyyy_xxy_0,  \
                             ta1_x_xyyyy_xxy_1,  \
                             ta1_x_xyyyy_xy_0,   \
                             ta1_x_xyyyy_xy_1,   \
                             ta1_x_xyyyy_xyy_0,  \
                             ta1_x_xyyyy_xyy_1,  \
                             ta1_x_xyyyy_xyz_0,  \
                             ta1_x_xyyyy_xyz_1,  \
                             ta1_x_xyyyy_yyy_0,  \
                             ta1_x_xyyyy_yyy_1,  \
                             ta1_x_xyyyyz_xxx_0, \
                             ta1_x_xyyyyz_xxy_0, \
                             ta1_x_xyyyyz_xxz_0, \
                             ta1_x_xyyyyz_xyy_0, \
                             ta1_x_xyyyyz_xyz_0, \
                             ta1_x_xyyyyz_xzz_0, \
                             ta1_x_xyyyyz_yyy_0, \
                             ta1_x_xyyyyz_yyz_0, \
                             ta1_x_xyyyyz_yzz_0, \
                             ta1_x_xyyyyz_zzz_0, \
                             ta1_x_xyyyz_xxz_0,  \
                             ta1_x_xyyyz_xxz_1,  \
                             ta1_x_xyyyz_xzz_0,  \
                             ta1_x_xyyyz_xzz_1,  \
                             ta1_x_xyyz_xxz_0,   \
                             ta1_x_xyyz_xxz_1,   \
                             ta1_x_xyyz_xzz_0,   \
                             ta1_x_xyyz_xzz_1,   \
                             ta1_x_yyyyz_yyz_0,  \
                             ta1_x_yyyyz_yyz_1,  \
                             ta1_x_yyyyz_yzz_0,  \
                             ta1_x_yyyyz_yzz_1,  \
                             ta1_x_yyyyz_zzz_0,  \
                             ta1_x_yyyyz_zzz_1,  \
                             ta_yyyyz_yyz_1,     \
                             ta_yyyyz_yzz_1,     \
                             ta_yyyyz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyyz_xxx_0[i] = ta1_x_xyyyy_xxx_0[i] * pa_z[i] - ta1_x_xyyyy_xxx_1[i] * pc_z[i];

        ta1_x_xyyyyz_xxy_0[i] = ta1_x_xyyyy_xxy_0[i] * pa_z[i] - ta1_x_xyyyy_xxy_1[i] * pc_z[i];

        ta1_x_xyyyyz_xxz_0[i] =
            3.0 * ta1_x_xyyz_xxz_0[i] * fe_0 - 3.0 * ta1_x_xyyz_xxz_1[i] * fe_0 + ta1_x_xyyyz_xxz_0[i] * pa_y[i] - ta1_x_xyyyz_xxz_1[i] * pc_y[i];

        ta1_x_xyyyyz_xyy_0[i] = ta1_x_xyyyy_xyy_0[i] * pa_z[i] - ta1_x_xyyyy_xyy_1[i] * pc_z[i];

        ta1_x_xyyyyz_xyz_0[i] =
            ta1_x_xyyyy_xy_0[i] * fe_0 - ta1_x_xyyyy_xy_1[i] * fe_0 + ta1_x_xyyyy_xyz_0[i] * pa_z[i] - ta1_x_xyyyy_xyz_1[i] * pc_z[i];

        ta1_x_xyyyyz_xzz_0[i] =
            3.0 * ta1_x_xyyz_xzz_0[i] * fe_0 - 3.0 * ta1_x_xyyz_xzz_1[i] * fe_0 + ta1_x_xyyyz_xzz_0[i] * pa_y[i] - ta1_x_xyyyz_xzz_1[i] * pc_y[i];

        ta1_x_xyyyyz_yyy_0[i] = ta1_x_xyyyy_yyy_0[i] * pa_z[i] - ta1_x_xyyyy_yyy_1[i] * pc_z[i];

        ta1_x_xyyyyz_yyz_0[i] = ta_yyyyz_yyz_1[i] + ta1_x_yyyyz_yyz_0[i] * pa_x[i] - ta1_x_yyyyz_yyz_1[i] * pc_x[i];

        ta1_x_xyyyyz_yzz_0[i] = ta_yyyyz_yzz_1[i] + ta1_x_yyyyz_yzz_0[i] * pa_x[i] - ta1_x_yyyyz_yzz_1[i] * pc_x[i];

        ta1_x_xyyyyz_zzz_0[i] = ta_yyyyz_zzz_1[i] + ta1_x_yyyyz_zzz_0[i] * pa_x[i] - ta1_x_yyyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 170-180 components of targeted buffer : IF

    auto ta1_x_xyyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 170);

    auto ta1_x_xyyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 171);

    auto ta1_x_xyyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 172);

    auto ta1_x_xyyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 173);

    auto ta1_x_xyyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 174);

    auto ta1_x_xyyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 175);

    auto ta1_x_xyyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 176);

    auto ta1_x_xyyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 177);

    auto ta1_x_xyyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 178);

    auto ta1_x_xyyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 179);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xyyy_xxy_0,   \
                             ta1_x_xyyy_xxy_1,   \
                             ta1_x_xyyy_xyy_0,   \
                             ta1_x_xyyy_xyy_1,   \
                             ta1_x_xyyyz_xxy_0,  \
                             ta1_x_xyyyz_xxy_1,  \
                             ta1_x_xyyyz_xyy_0,  \
                             ta1_x_xyyyz_xyy_1,  \
                             ta1_x_xyyyzz_xxx_0, \
                             ta1_x_xyyyzz_xxy_0, \
                             ta1_x_xyyyzz_xxz_0, \
                             ta1_x_xyyyzz_xyy_0, \
                             ta1_x_xyyyzz_xyz_0, \
                             ta1_x_xyyyzz_xzz_0, \
                             ta1_x_xyyyzz_yyy_0, \
                             ta1_x_xyyyzz_yyz_0, \
                             ta1_x_xyyyzz_yzz_0, \
                             ta1_x_xyyyzz_zzz_0, \
                             ta1_x_xyyzz_xxx_0,  \
                             ta1_x_xyyzz_xxx_1,  \
                             ta1_x_xyyzz_xxz_0,  \
                             ta1_x_xyyzz_xxz_1,  \
                             ta1_x_xyyzz_xzz_0,  \
                             ta1_x_xyyzz_xzz_1,  \
                             ta1_x_xyzz_xxx_0,   \
                             ta1_x_xyzz_xxx_1,   \
                             ta1_x_xyzz_xxz_0,   \
                             ta1_x_xyzz_xxz_1,   \
                             ta1_x_xyzz_xzz_0,   \
                             ta1_x_xyzz_xzz_1,   \
                             ta1_x_yyyzz_xyz_0,  \
                             ta1_x_yyyzz_xyz_1,  \
                             ta1_x_yyyzz_yyy_0,  \
                             ta1_x_yyyzz_yyy_1,  \
                             ta1_x_yyyzz_yyz_0,  \
                             ta1_x_yyyzz_yyz_1,  \
                             ta1_x_yyyzz_yz_0,   \
                             ta1_x_yyyzz_yz_1,   \
                             ta1_x_yyyzz_yzz_0,  \
                             ta1_x_yyyzz_yzz_1,  \
                             ta1_x_yyyzz_zzz_0,  \
                             ta1_x_yyyzz_zzz_1,  \
                             ta_yyyzz_xyz_1,     \
                             ta_yyyzz_yyy_1,     \
                             ta_yyyzz_yyz_1,     \
                             ta_yyyzz_yzz_1,     \
                             ta_yyyzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyzz_xxx_0[i] =
            2.0 * ta1_x_xyzz_xxx_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xxx_1[i] * fe_0 + ta1_x_xyyzz_xxx_0[i] * pa_y[i] - ta1_x_xyyzz_xxx_1[i] * pc_y[i];

        ta1_x_xyyyzz_xxy_0[i] =
            ta1_x_xyyy_xxy_0[i] * fe_0 - ta1_x_xyyy_xxy_1[i] * fe_0 + ta1_x_xyyyz_xxy_0[i] * pa_z[i] - ta1_x_xyyyz_xxy_1[i] * pc_z[i];

        ta1_x_xyyyzz_xxz_0[i] =
            2.0 * ta1_x_xyzz_xxz_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xxz_1[i] * fe_0 + ta1_x_xyyzz_xxz_0[i] * pa_y[i] - ta1_x_xyyzz_xxz_1[i] * pc_y[i];

        ta1_x_xyyyzz_xyy_0[i] =
            ta1_x_xyyy_xyy_0[i] * fe_0 - ta1_x_xyyy_xyy_1[i] * fe_0 + ta1_x_xyyyz_xyy_0[i] * pa_z[i] - ta1_x_xyyyz_xyy_1[i] * pc_z[i];

        ta1_x_xyyyzz_xyz_0[i] = ta1_x_yyyzz_yz_0[i] * fe_0 - ta1_x_yyyzz_yz_1[i] * fe_0 + ta_yyyzz_xyz_1[i] + ta1_x_yyyzz_xyz_0[i] * pa_x[i] -
                                ta1_x_yyyzz_xyz_1[i] * pc_x[i];

        ta1_x_xyyyzz_xzz_0[i] =
            2.0 * ta1_x_xyzz_xzz_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xzz_1[i] * fe_0 + ta1_x_xyyzz_xzz_0[i] * pa_y[i] - ta1_x_xyyzz_xzz_1[i] * pc_y[i];

        ta1_x_xyyyzz_yyy_0[i] = ta_yyyzz_yyy_1[i] + ta1_x_yyyzz_yyy_0[i] * pa_x[i] - ta1_x_yyyzz_yyy_1[i] * pc_x[i];

        ta1_x_xyyyzz_yyz_0[i] = ta_yyyzz_yyz_1[i] + ta1_x_yyyzz_yyz_0[i] * pa_x[i] - ta1_x_yyyzz_yyz_1[i] * pc_x[i];

        ta1_x_xyyyzz_yzz_0[i] = ta_yyyzz_yzz_1[i] + ta1_x_yyyzz_yzz_0[i] * pa_x[i] - ta1_x_yyyzz_yzz_1[i] * pc_x[i];

        ta1_x_xyyyzz_zzz_0[i] = ta_yyyzz_zzz_1[i] + ta1_x_yyyzz_zzz_0[i] * pa_x[i] - ta1_x_yyyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 180-190 components of targeted buffer : IF

    auto ta1_x_xyyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 180);

    auto ta1_x_xyyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 181);

    auto ta1_x_xyyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 182);

    auto ta1_x_xyyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 183);

    auto ta1_x_xyyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 184);

    auto ta1_x_xyyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 185);

    auto ta1_x_xyyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 186);

    auto ta1_x_xyyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 187);

    auto ta1_x_xyyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 188);

    auto ta1_x_xyyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 189);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xyyz_xxy_0,   \
                             ta1_x_xyyz_xxy_1,   \
                             ta1_x_xyyz_xyy_0,   \
                             ta1_x_xyyz_xyy_1,   \
                             ta1_x_xyyzz_xxy_0,  \
                             ta1_x_xyyzz_xxy_1,  \
                             ta1_x_xyyzz_xyy_0,  \
                             ta1_x_xyyzz_xyy_1,  \
                             ta1_x_xyyzzz_xxx_0, \
                             ta1_x_xyyzzz_xxy_0, \
                             ta1_x_xyyzzz_xxz_0, \
                             ta1_x_xyyzzz_xyy_0, \
                             ta1_x_xyyzzz_xyz_0, \
                             ta1_x_xyyzzz_xzz_0, \
                             ta1_x_xyyzzz_yyy_0, \
                             ta1_x_xyyzzz_yyz_0, \
                             ta1_x_xyyzzz_yzz_0, \
                             ta1_x_xyyzzz_zzz_0, \
                             ta1_x_xyzzz_xxx_0,  \
                             ta1_x_xyzzz_xxx_1,  \
                             ta1_x_xyzzz_xxz_0,  \
                             ta1_x_xyzzz_xxz_1,  \
                             ta1_x_xyzzz_xzz_0,  \
                             ta1_x_xyzzz_xzz_1,  \
                             ta1_x_xzzz_xxx_0,   \
                             ta1_x_xzzz_xxx_1,   \
                             ta1_x_xzzz_xxz_0,   \
                             ta1_x_xzzz_xxz_1,   \
                             ta1_x_xzzz_xzz_0,   \
                             ta1_x_xzzz_xzz_1,   \
                             ta1_x_yyzzz_xyz_0,  \
                             ta1_x_yyzzz_xyz_1,  \
                             ta1_x_yyzzz_yyy_0,  \
                             ta1_x_yyzzz_yyy_1,  \
                             ta1_x_yyzzz_yyz_0,  \
                             ta1_x_yyzzz_yyz_1,  \
                             ta1_x_yyzzz_yz_0,   \
                             ta1_x_yyzzz_yz_1,   \
                             ta1_x_yyzzz_yzz_0,  \
                             ta1_x_yyzzz_yzz_1,  \
                             ta1_x_yyzzz_zzz_0,  \
                             ta1_x_yyzzz_zzz_1,  \
                             ta_yyzzz_xyz_1,     \
                             ta_yyzzz_yyy_1,     \
                             ta_yyzzz_yyz_1,     \
                             ta_yyzzz_yzz_1,     \
                             ta_yyzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzzz_xxx_0[i] =
            ta1_x_xzzz_xxx_0[i] * fe_0 - ta1_x_xzzz_xxx_1[i] * fe_0 + ta1_x_xyzzz_xxx_0[i] * pa_y[i] - ta1_x_xyzzz_xxx_1[i] * pc_y[i];

        ta1_x_xyyzzz_xxy_0[i] =
            2.0 * ta1_x_xyyz_xxy_0[i] * fe_0 - 2.0 * ta1_x_xyyz_xxy_1[i] * fe_0 + ta1_x_xyyzz_xxy_0[i] * pa_z[i] - ta1_x_xyyzz_xxy_1[i] * pc_z[i];

        ta1_x_xyyzzz_xxz_0[i] =
            ta1_x_xzzz_xxz_0[i] * fe_0 - ta1_x_xzzz_xxz_1[i] * fe_0 + ta1_x_xyzzz_xxz_0[i] * pa_y[i] - ta1_x_xyzzz_xxz_1[i] * pc_y[i];

        ta1_x_xyyzzz_xyy_0[i] =
            2.0 * ta1_x_xyyz_xyy_0[i] * fe_0 - 2.0 * ta1_x_xyyz_xyy_1[i] * fe_0 + ta1_x_xyyzz_xyy_0[i] * pa_z[i] - ta1_x_xyyzz_xyy_1[i] * pc_z[i];

        ta1_x_xyyzzz_xyz_0[i] = ta1_x_yyzzz_yz_0[i] * fe_0 - ta1_x_yyzzz_yz_1[i] * fe_0 + ta_yyzzz_xyz_1[i] + ta1_x_yyzzz_xyz_0[i] * pa_x[i] -
                                ta1_x_yyzzz_xyz_1[i] * pc_x[i];

        ta1_x_xyyzzz_xzz_0[i] =
            ta1_x_xzzz_xzz_0[i] * fe_0 - ta1_x_xzzz_xzz_1[i] * fe_0 + ta1_x_xyzzz_xzz_0[i] * pa_y[i] - ta1_x_xyzzz_xzz_1[i] * pc_y[i];

        ta1_x_xyyzzz_yyy_0[i] = ta_yyzzz_yyy_1[i] + ta1_x_yyzzz_yyy_0[i] * pa_x[i] - ta1_x_yyzzz_yyy_1[i] * pc_x[i];

        ta1_x_xyyzzz_yyz_0[i] = ta_yyzzz_yyz_1[i] + ta1_x_yyzzz_yyz_0[i] * pa_x[i] - ta1_x_yyzzz_yyz_1[i] * pc_x[i];

        ta1_x_xyyzzz_yzz_0[i] = ta_yyzzz_yzz_1[i] + ta1_x_yyzzz_yzz_0[i] * pa_x[i] - ta1_x_yyzzz_yzz_1[i] * pc_x[i];

        ta1_x_xyyzzz_zzz_0[i] = ta_yyzzz_zzz_1[i] + ta1_x_yyzzz_zzz_0[i] * pa_x[i] - ta1_x_yyzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 190-200 components of targeted buffer : IF

    auto ta1_x_xyzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 190);

    auto ta1_x_xyzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 191);

    auto ta1_x_xyzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 192);

    auto ta1_x_xyzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 193);

    auto ta1_x_xyzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 194);

    auto ta1_x_xyzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 195);

    auto ta1_x_xyzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 196);

    auto ta1_x_xyzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 197);

    auto ta1_x_xyzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 198);

    auto ta1_x_xyzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 199);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xyzzzz_xxx_0, \
                             ta1_x_xyzzzz_xxy_0, \
                             ta1_x_xyzzzz_xxz_0, \
                             ta1_x_xyzzzz_xyy_0, \
                             ta1_x_xyzzzz_xyz_0, \
                             ta1_x_xyzzzz_xzz_0, \
                             ta1_x_xyzzzz_yyy_0, \
                             ta1_x_xyzzzz_yyz_0, \
                             ta1_x_xyzzzz_yzz_0, \
                             ta1_x_xyzzzz_zzz_0, \
                             ta1_x_xzzzz_xx_0,   \
                             ta1_x_xzzzz_xx_1,   \
                             ta1_x_xzzzz_xxx_0,  \
                             ta1_x_xzzzz_xxx_1,  \
                             ta1_x_xzzzz_xxy_0,  \
                             ta1_x_xzzzz_xxy_1,  \
                             ta1_x_xzzzz_xxz_0,  \
                             ta1_x_xzzzz_xxz_1,  \
                             ta1_x_xzzzz_xy_0,   \
                             ta1_x_xzzzz_xy_1,   \
                             ta1_x_xzzzz_xyy_0,  \
                             ta1_x_xzzzz_xyy_1,  \
                             ta1_x_xzzzz_xyz_0,  \
                             ta1_x_xzzzz_xyz_1,  \
                             ta1_x_xzzzz_xz_0,   \
                             ta1_x_xzzzz_xz_1,   \
                             ta1_x_xzzzz_xzz_0,  \
                             ta1_x_xzzzz_xzz_1,  \
                             ta1_x_xzzzz_zzz_0,  \
                             ta1_x_xzzzz_zzz_1,  \
                             ta1_x_yzzzz_yyy_0,  \
                             ta1_x_yzzzz_yyy_1,  \
                             ta1_x_yzzzz_yyz_0,  \
                             ta1_x_yzzzz_yyz_1,  \
                             ta1_x_yzzzz_yzz_0,  \
                             ta1_x_yzzzz_yzz_1,  \
                             ta_yzzzz_yyy_1,     \
                             ta_yzzzz_yyz_1,     \
                             ta_yzzzz_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzzzz_xxx_0[i] = ta1_x_xzzzz_xxx_0[i] * pa_y[i] - ta1_x_xzzzz_xxx_1[i] * pc_y[i];

        ta1_x_xyzzzz_xxy_0[i] =
            ta1_x_xzzzz_xx_0[i] * fe_0 - ta1_x_xzzzz_xx_1[i] * fe_0 + ta1_x_xzzzz_xxy_0[i] * pa_y[i] - ta1_x_xzzzz_xxy_1[i] * pc_y[i];

        ta1_x_xyzzzz_xxz_0[i] = ta1_x_xzzzz_xxz_0[i] * pa_y[i] - ta1_x_xzzzz_xxz_1[i] * pc_y[i];

        ta1_x_xyzzzz_xyy_0[i] =
            2.0 * ta1_x_xzzzz_xy_0[i] * fe_0 - 2.0 * ta1_x_xzzzz_xy_1[i] * fe_0 + ta1_x_xzzzz_xyy_0[i] * pa_y[i] - ta1_x_xzzzz_xyy_1[i] * pc_y[i];

        ta1_x_xyzzzz_xyz_0[i] =
            ta1_x_xzzzz_xz_0[i] * fe_0 - ta1_x_xzzzz_xz_1[i] * fe_0 + ta1_x_xzzzz_xyz_0[i] * pa_y[i] - ta1_x_xzzzz_xyz_1[i] * pc_y[i];

        ta1_x_xyzzzz_xzz_0[i] = ta1_x_xzzzz_xzz_0[i] * pa_y[i] - ta1_x_xzzzz_xzz_1[i] * pc_y[i];

        ta1_x_xyzzzz_yyy_0[i] = ta_yzzzz_yyy_1[i] + ta1_x_yzzzz_yyy_0[i] * pa_x[i] - ta1_x_yzzzz_yyy_1[i] * pc_x[i];

        ta1_x_xyzzzz_yyz_0[i] = ta_yzzzz_yyz_1[i] + ta1_x_yzzzz_yyz_0[i] * pa_x[i] - ta1_x_yzzzz_yyz_1[i] * pc_x[i];

        ta1_x_xyzzzz_yzz_0[i] = ta_yzzzz_yzz_1[i] + ta1_x_yzzzz_yzz_0[i] * pa_x[i] - ta1_x_yzzzz_yzz_1[i] * pc_x[i];

        ta1_x_xyzzzz_zzz_0[i] = ta1_x_xzzzz_zzz_0[i] * pa_y[i] - ta1_x_xzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 200-210 components of targeted buffer : IF

    auto ta1_x_xzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 200);

    auto ta1_x_xzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 201);

    auto ta1_x_xzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 202);

    auto ta1_x_xzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 203);

    auto ta1_x_xzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 204);

    auto ta1_x_xzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 205);

    auto ta1_x_xzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 206);

    auto ta1_x_xzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 207);

    auto ta1_x_xzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 208);

    auto ta1_x_xzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 209);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xzzz_xxx_0,   \
                             ta1_x_xzzz_xxx_1,   \
                             ta1_x_xzzz_xxy_0,   \
                             ta1_x_xzzz_xxy_1,   \
                             ta1_x_xzzz_xyy_0,   \
                             ta1_x_xzzz_xyy_1,   \
                             ta1_x_xzzzz_xxx_0,  \
                             ta1_x_xzzzz_xxx_1,  \
                             ta1_x_xzzzz_xxy_0,  \
                             ta1_x_xzzzz_xxy_1,  \
                             ta1_x_xzzzz_xyy_0,  \
                             ta1_x_xzzzz_xyy_1,  \
                             ta1_x_xzzzzz_xxx_0, \
                             ta1_x_xzzzzz_xxy_0, \
                             ta1_x_xzzzzz_xxz_0, \
                             ta1_x_xzzzzz_xyy_0, \
                             ta1_x_xzzzzz_xyz_0, \
                             ta1_x_xzzzzz_xzz_0, \
                             ta1_x_xzzzzz_yyy_0, \
                             ta1_x_xzzzzz_yyz_0, \
                             ta1_x_xzzzzz_yzz_0, \
                             ta1_x_xzzzzz_zzz_0, \
                             ta1_x_zzzzz_xxz_0,  \
                             ta1_x_zzzzz_xxz_1,  \
                             ta1_x_zzzzz_xyz_0,  \
                             ta1_x_zzzzz_xyz_1,  \
                             ta1_x_zzzzz_xz_0,   \
                             ta1_x_zzzzz_xz_1,   \
                             ta1_x_zzzzz_xzz_0,  \
                             ta1_x_zzzzz_xzz_1,  \
                             ta1_x_zzzzz_yyy_0,  \
                             ta1_x_zzzzz_yyy_1,  \
                             ta1_x_zzzzz_yyz_0,  \
                             ta1_x_zzzzz_yyz_1,  \
                             ta1_x_zzzzz_yz_0,   \
                             ta1_x_zzzzz_yz_1,   \
                             ta1_x_zzzzz_yzz_0,  \
                             ta1_x_zzzzz_yzz_1,  \
                             ta1_x_zzzzz_zz_0,   \
                             ta1_x_zzzzz_zz_1,   \
                             ta1_x_zzzzz_zzz_0,  \
                             ta1_x_zzzzz_zzz_1,  \
                             ta_zzzzz_xxz_1,     \
                             ta_zzzzz_xyz_1,     \
                             ta_zzzzz_xzz_1,     \
                             ta_zzzzz_yyy_1,     \
                             ta_zzzzz_yyz_1,     \
                             ta_zzzzz_yzz_1,     \
                             ta_zzzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzzz_xxx_0[i] =
            4.0 * ta1_x_xzzz_xxx_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xxx_1[i] * fe_0 + ta1_x_xzzzz_xxx_0[i] * pa_z[i] - ta1_x_xzzzz_xxx_1[i] * pc_z[i];

        ta1_x_xzzzzz_xxy_0[i] =
            4.0 * ta1_x_xzzz_xxy_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xxy_1[i] * fe_0 + ta1_x_xzzzz_xxy_0[i] * pa_z[i] - ta1_x_xzzzz_xxy_1[i] * pc_z[i];

        ta1_x_xzzzzz_xxz_0[i] = 2.0 * ta1_x_zzzzz_xz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xz_1[i] * fe_0 + ta_zzzzz_xxz_1[i] +
                                ta1_x_zzzzz_xxz_0[i] * pa_x[i] - ta1_x_zzzzz_xxz_1[i] * pc_x[i];

        ta1_x_xzzzzz_xyy_0[i] =
            4.0 * ta1_x_xzzz_xyy_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xyy_1[i] * fe_0 + ta1_x_xzzzz_xyy_0[i] * pa_z[i] - ta1_x_xzzzz_xyy_1[i] * pc_z[i];

        ta1_x_xzzzzz_xyz_0[i] = ta1_x_zzzzz_yz_0[i] * fe_0 - ta1_x_zzzzz_yz_1[i] * fe_0 + ta_zzzzz_xyz_1[i] + ta1_x_zzzzz_xyz_0[i] * pa_x[i] -
                                ta1_x_zzzzz_xyz_1[i] * pc_x[i];

        ta1_x_xzzzzz_xzz_0[i] = ta1_x_zzzzz_zz_0[i] * fe_0 - ta1_x_zzzzz_zz_1[i] * fe_0 + ta_zzzzz_xzz_1[i] + ta1_x_zzzzz_xzz_0[i] * pa_x[i] -
                                ta1_x_zzzzz_xzz_1[i] * pc_x[i];

        ta1_x_xzzzzz_yyy_0[i] = ta_zzzzz_yyy_1[i] + ta1_x_zzzzz_yyy_0[i] * pa_x[i] - ta1_x_zzzzz_yyy_1[i] * pc_x[i];

        ta1_x_xzzzzz_yyz_0[i] = ta_zzzzz_yyz_1[i] + ta1_x_zzzzz_yyz_0[i] * pa_x[i] - ta1_x_zzzzz_yyz_1[i] * pc_x[i];

        ta1_x_xzzzzz_yzz_0[i] = ta_zzzzz_yzz_1[i] + ta1_x_zzzzz_yzz_0[i] * pa_x[i] - ta1_x_zzzzz_yzz_1[i] * pc_x[i];

        ta1_x_xzzzzz_zzz_0[i] = ta_zzzzz_zzz_1[i] + ta1_x_zzzzz_zzz_0[i] * pa_x[i] - ta1_x_zzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 210-220 components of targeted buffer : IF

    auto ta1_x_yyyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 210);

    auto ta1_x_yyyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 211);

    auto ta1_x_yyyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 212);

    auto ta1_x_yyyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 213);

    auto ta1_x_yyyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 214);

    auto ta1_x_yyyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 215);

    auto ta1_x_yyyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 216);

    auto ta1_x_yyyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 217);

    auto ta1_x_yyyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 218);

    auto ta1_x_yyyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 219);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yyyy_xxx_0,   \
                             ta1_x_yyyy_xxx_1,   \
                             ta1_x_yyyy_xxy_0,   \
                             ta1_x_yyyy_xxy_1,   \
                             ta1_x_yyyy_xxz_0,   \
                             ta1_x_yyyy_xxz_1,   \
                             ta1_x_yyyy_xyy_0,   \
                             ta1_x_yyyy_xyy_1,   \
                             ta1_x_yyyy_xyz_0,   \
                             ta1_x_yyyy_xyz_1,   \
                             ta1_x_yyyy_xzz_0,   \
                             ta1_x_yyyy_xzz_1,   \
                             ta1_x_yyyy_yyy_0,   \
                             ta1_x_yyyy_yyy_1,   \
                             ta1_x_yyyy_yyz_0,   \
                             ta1_x_yyyy_yyz_1,   \
                             ta1_x_yyyy_yzz_0,   \
                             ta1_x_yyyy_yzz_1,   \
                             ta1_x_yyyy_zzz_0,   \
                             ta1_x_yyyy_zzz_1,   \
                             ta1_x_yyyyy_xx_0,   \
                             ta1_x_yyyyy_xx_1,   \
                             ta1_x_yyyyy_xxx_0,  \
                             ta1_x_yyyyy_xxx_1,  \
                             ta1_x_yyyyy_xxy_0,  \
                             ta1_x_yyyyy_xxy_1,  \
                             ta1_x_yyyyy_xxz_0,  \
                             ta1_x_yyyyy_xxz_1,  \
                             ta1_x_yyyyy_xy_0,   \
                             ta1_x_yyyyy_xy_1,   \
                             ta1_x_yyyyy_xyy_0,  \
                             ta1_x_yyyyy_xyy_1,  \
                             ta1_x_yyyyy_xyz_0,  \
                             ta1_x_yyyyy_xyz_1,  \
                             ta1_x_yyyyy_xz_0,   \
                             ta1_x_yyyyy_xz_1,   \
                             ta1_x_yyyyy_xzz_0,  \
                             ta1_x_yyyyy_xzz_1,  \
                             ta1_x_yyyyy_yy_0,   \
                             ta1_x_yyyyy_yy_1,   \
                             ta1_x_yyyyy_yyy_0,  \
                             ta1_x_yyyyy_yyy_1,  \
                             ta1_x_yyyyy_yyz_0,  \
                             ta1_x_yyyyy_yyz_1,  \
                             ta1_x_yyyyy_yz_0,   \
                             ta1_x_yyyyy_yz_1,   \
                             ta1_x_yyyyy_yzz_0,  \
                             ta1_x_yyyyy_yzz_1,  \
                             ta1_x_yyyyy_zz_0,   \
                             ta1_x_yyyyy_zz_1,   \
                             ta1_x_yyyyy_zzz_0,  \
                             ta1_x_yyyyy_zzz_1,  \
                             ta1_x_yyyyyy_xxx_0, \
                             ta1_x_yyyyyy_xxy_0, \
                             ta1_x_yyyyyy_xxz_0, \
                             ta1_x_yyyyyy_xyy_0, \
                             ta1_x_yyyyyy_xyz_0, \
                             ta1_x_yyyyyy_xzz_0, \
                             ta1_x_yyyyyy_yyy_0, \
                             ta1_x_yyyyyy_yyz_0, \
                             ta1_x_yyyyyy_yzz_0, \
                             ta1_x_yyyyyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyy_xxx_0[i] =
            5.0 * ta1_x_yyyy_xxx_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxx_1[i] * fe_0 + ta1_x_yyyyy_xxx_0[i] * pa_y[i] - ta1_x_yyyyy_xxx_1[i] * pc_y[i];

        ta1_x_yyyyyy_xxy_0[i] = 5.0 * ta1_x_yyyy_xxy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxy_1[i] * fe_0 + ta1_x_yyyyy_xx_0[i] * fe_0 -
                                ta1_x_yyyyy_xx_1[i] * fe_0 + ta1_x_yyyyy_xxy_0[i] * pa_y[i] - ta1_x_yyyyy_xxy_1[i] * pc_y[i];

        ta1_x_yyyyyy_xxz_0[i] =
            5.0 * ta1_x_yyyy_xxz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xxz_1[i] * fe_0 + ta1_x_yyyyy_xxz_0[i] * pa_y[i] - ta1_x_yyyyy_xxz_1[i] * pc_y[i];

        ta1_x_yyyyyy_xyy_0[i] = 5.0 * ta1_x_yyyy_xyy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xyy_1[i] * fe_0 + 2.0 * ta1_x_yyyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_x_yyyyy_xy_1[i] * fe_0 + ta1_x_yyyyy_xyy_0[i] * pa_y[i] - ta1_x_yyyyy_xyy_1[i] * pc_y[i];

        ta1_x_yyyyyy_xyz_0[i] = 5.0 * ta1_x_yyyy_xyz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xyz_1[i] * fe_0 + ta1_x_yyyyy_xz_0[i] * fe_0 -
                                ta1_x_yyyyy_xz_1[i] * fe_0 + ta1_x_yyyyy_xyz_0[i] * pa_y[i] - ta1_x_yyyyy_xyz_1[i] * pc_y[i];

        ta1_x_yyyyyy_xzz_0[i] =
            5.0 * ta1_x_yyyy_xzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xzz_1[i] * fe_0 + ta1_x_yyyyy_xzz_0[i] * pa_y[i] - ta1_x_yyyyy_xzz_1[i] * pc_y[i];

        ta1_x_yyyyyy_yyy_0[i] = 5.0 * ta1_x_yyyy_yyy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yyy_1[i] * fe_0 + 3.0 * ta1_x_yyyyy_yy_0[i] * fe_0 -
                                3.0 * ta1_x_yyyyy_yy_1[i] * fe_0 + ta1_x_yyyyy_yyy_0[i] * pa_y[i] - ta1_x_yyyyy_yyy_1[i] * pc_y[i];

        ta1_x_yyyyyy_yyz_0[i] = 5.0 * ta1_x_yyyy_yyz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yyz_1[i] * fe_0 + 2.0 * ta1_x_yyyyy_yz_0[i] * fe_0 -
                                2.0 * ta1_x_yyyyy_yz_1[i] * fe_0 + ta1_x_yyyyy_yyz_0[i] * pa_y[i] - ta1_x_yyyyy_yyz_1[i] * pc_y[i];

        ta1_x_yyyyyy_yzz_0[i] = 5.0 * ta1_x_yyyy_yzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yzz_1[i] * fe_0 + ta1_x_yyyyy_zz_0[i] * fe_0 -
                                ta1_x_yyyyy_zz_1[i] * fe_0 + ta1_x_yyyyy_yzz_0[i] * pa_y[i] - ta1_x_yyyyy_yzz_1[i] * pc_y[i];

        ta1_x_yyyyyy_zzz_0[i] =
            5.0 * ta1_x_yyyy_zzz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_zzz_1[i] * fe_0 + ta1_x_yyyyy_zzz_0[i] * pa_y[i] - ta1_x_yyyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 220-230 components of targeted buffer : IF

    auto ta1_x_yyyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 220);

    auto ta1_x_yyyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 221);

    auto ta1_x_yyyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 222);

    auto ta1_x_yyyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 223);

    auto ta1_x_yyyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 224);

    auto ta1_x_yyyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 225);

    auto ta1_x_yyyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 226);

    auto ta1_x_yyyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 227);

    auto ta1_x_yyyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 228);

    auto ta1_x_yyyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 229);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyyyy_xxx_0,  \
                             ta1_x_yyyyy_xxx_1,  \
                             ta1_x_yyyyy_xxy_0,  \
                             ta1_x_yyyyy_xxy_1,  \
                             ta1_x_yyyyy_xy_0,   \
                             ta1_x_yyyyy_xy_1,   \
                             ta1_x_yyyyy_xyy_0,  \
                             ta1_x_yyyyy_xyy_1,  \
                             ta1_x_yyyyy_xyz_0,  \
                             ta1_x_yyyyy_xyz_1,  \
                             ta1_x_yyyyy_yy_0,   \
                             ta1_x_yyyyy_yy_1,   \
                             ta1_x_yyyyy_yyy_0,  \
                             ta1_x_yyyyy_yyy_1,  \
                             ta1_x_yyyyy_yyz_0,  \
                             ta1_x_yyyyy_yyz_1,  \
                             ta1_x_yyyyy_yz_0,   \
                             ta1_x_yyyyy_yz_1,   \
                             ta1_x_yyyyy_yzz_0,  \
                             ta1_x_yyyyy_yzz_1,  \
                             ta1_x_yyyyyz_xxx_0, \
                             ta1_x_yyyyyz_xxy_0, \
                             ta1_x_yyyyyz_xxz_0, \
                             ta1_x_yyyyyz_xyy_0, \
                             ta1_x_yyyyyz_xyz_0, \
                             ta1_x_yyyyyz_xzz_0, \
                             ta1_x_yyyyyz_yyy_0, \
                             ta1_x_yyyyyz_yyz_0, \
                             ta1_x_yyyyyz_yzz_0, \
                             ta1_x_yyyyyz_zzz_0, \
                             ta1_x_yyyyz_xxz_0,  \
                             ta1_x_yyyyz_xxz_1,  \
                             ta1_x_yyyyz_xzz_0,  \
                             ta1_x_yyyyz_xzz_1,  \
                             ta1_x_yyyyz_zzz_0,  \
                             ta1_x_yyyyz_zzz_1,  \
                             ta1_x_yyyz_xxz_0,   \
                             ta1_x_yyyz_xxz_1,   \
                             ta1_x_yyyz_xzz_0,   \
                             ta1_x_yyyz_xzz_1,   \
                             ta1_x_yyyz_zzz_0,   \
                             ta1_x_yyyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyz_xxx_0[i] = ta1_x_yyyyy_xxx_0[i] * pa_z[i] - ta1_x_yyyyy_xxx_1[i] * pc_z[i];

        ta1_x_yyyyyz_xxy_0[i] = ta1_x_yyyyy_xxy_0[i] * pa_z[i] - ta1_x_yyyyy_xxy_1[i] * pc_z[i];

        ta1_x_yyyyyz_xxz_0[i] =
            4.0 * ta1_x_yyyz_xxz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_xxz_1[i] * fe_0 + ta1_x_yyyyz_xxz_0[i] * pa_y[i] - ta1_x_yyyyz_xxz_1[i] * pc_y[i];

        ta1_x_yyyyyz_xyy_0[i] = ta1_x_yyyyy_xyy_0[i] * pa_z[i] - ta1_x_yyyyy_xyy_1[i] * pc_z[i];

        ta1_x_yyyyyz_xyz_0[i] =
            ta1_x_yyyyy_xy_0[i] * fe_0 - ta1_x_yyyyy_xy_1[i] * fe_0 + ta1_x_yyyyy_xyz_0[i] * pa_z[i] - ta1_x_yyyyy_xyz_1[i] * pc_z[i];

        ta1_x_yyyyyz_xzz_0[i] =
            4.0 * ta1_x_yyyz_xzz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_xzz_1[i] * fe_0 + ta1_x_yyyyz_xzz_0[i] * pa_y[i] - ta1_x_yyyyz_xzz_1[i] * pc_y[i];

        ta1_x_yyyyyz_yyy_0[i] = ta1_x_yyyyy_yyy_0[i] * pa_z[i] - ta1_x_yyyyy_yyy_1[i] * pc_z[i];

        ta1_x_yyyyyz_yyz_0[i] =
            ta1_x_yyyyy_yy_0[i] * fe_0 - ta1_x_yyyyy_yy_1[i] * fe_0 + ta1_x_yyyyy_yyz_0[i] * pa_z[i] - ta1_x_yyyyy_yyz_1[i] * pc_z[i];

        ta1_x_yyyyyz_yzz_0[i] =
            2.0 * ta1_x_yyyyy_yz_0[i] * fe_0 - 2.0 * ta1_x_yyyyy_yz_1[i] * fe_0 + ta1_x_yyyyy_yzz_0[i] * pa_z[i] - ta1_x_yyyyy_yzz_1[i] * pc_z[i];

        ta1_x_yyyyyz_zzz_0[i] =
            4.0 * ta1_x_yyyz_zzz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_zzz_1[i] * fe_0 + ta1_x_yyyyz_zzz_0[i] * pa_y[i] - ta1_x_yyyyz_zzz_1[i] * pc_y[i];
    }

    // Set up 230-240 components of targeted buffer : IF

    auto ta1_x_yyyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 230);

    auto ta1_x_yyyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 231);

    auto ta1_x_yyyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 232);

    auto ta1_x_yyyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 233);

    auto ta1_x_yyyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 234);

    auto ta1_x_yyyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 235);

    auto ta1_x_yyyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 236);

    auto ta1_x_yyyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 237);

    auto ta1_x_yyyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 238);

    auto ta1_x_yyyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 239);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyyy_xxy_0,   \
                             ta1_x_yyyy_xxy_1,   \
                             ta1_x_yyyy_xyy_0,   \
                             ta1_x_yyyy_xyy_1,   \
                             ta1_x_yyyy_yyy_0,   \
                             ta1_x_yyyy_yyy_1,   \
                             ta1_x_yyyyz_xxy_0,  \
                             ta1_x_yyyyz_xxy_1,  \
                             ta1_x_yyyyz_xyy_0,  \
                             ta1_x_yyyyz_xyy_1,  \
                             ta1_x_yyyyz_yyy_0,  \
                             ta1_x_yyyyz_yyy_1,  \
                             ta1_x_yyyyzz_xxx_0, \
                             ta1_x_yyyyzz_xxy_0, \
                             ta1_x_yyyyzz_xxz_0, \
                             ta1_x_yyyyzz_xyy_0, \
                             ta1_x_yyyyzz_xyz_0, \
                             ta1_x_yyyyzz_xzz_0, \
                             ta1_x_yyyyzz_yyy_0, \
                             ta1_x_yyyyzz_yyz_0, \
                             ta1_x_yyyyzz_yzz_0, \
                             ta1_x_yyyyzz_zzz_0, \
                             ta1_x_yyyzz_xxx_0,  \
                             ta1_x_yyyzz_xxx_1,  \
                             ta1_x_yyyzz_xxz_0,  \
                             ta1_x_yyyzz_xxz_1,  \
                             ta1_x_yyyzz_xyz_0,  \
                             ta1_x_yyyzz_xyz_1,  \
                             ta1_x_yyyzz_xz_0,   \
                             ta1_x_yyyzz_xz_1,   \
                             ta1_x_yyyzz_xzz_0,  \
                             ta1_x_yyyzz_xzz_1,  \
                             ta1_x_yyyzz_yyz_0,  \
                             ta1_x_yyyzz_yyz_1,  \
                             ta1_x_yyyzz_yz_0,   \
                             ta1_x_yyyzz_yz_1,   \
                             ta1_x_yyyzz_yzz_0,  \
                             ta1_x_yyyzz_yzz_1,  \
                             ta1_x_yyyzz_zz_0,   \
                             ta1_x_yyyzz_zz_1,   \
                             ta1_x_yyyzz_zzz_0,  \
                             ta1_x_yyyzz_zzz_1,  \
                             ta1_x_yyzz_xxx_0,   \
                             ta1_x_yyzz_xxx_1,   \
                             ta1_x_yyzz_xxz_0,   \
                             ta1_x_yyzz_xxz_1,   \
                             ta1_x_yyzz_xyz_0,   \
                             ta1_x_yyzz_xyz_1,   \
                             ta1_x_yyzz_xzz_0,   \
                             ta1_x_yyzz_xzz_1,   \
                             ta1_x_yyzz_yyz_0,   \
                             ta1_x_yyzz_yyz_1,   \
                             ta1_x_yyzz_yzz_0,   \
                             ta1_x_yyzz_yzz_1,   \
                             ta1_x_yyzz_zzz_0,   \
                             ta1_x_yyzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyzz_xxx_0[i] =
            3.0 * ta1_x_yyzz_xxx_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxx_1[i] * fe_0 + ta1_x_yyyzz_xxx_0[i] * pa_y[i] - ta1_x_yyyzz_xxx_1[i] * pc_y[i];

        ta1_x_yyyyzz_xxy_0[i] =
            ta1_x_yyyy_xxy_0[i] * fe_0 - ta1_x_yyyy_xxy_1[i] * fe_0 + ta1_x_yyyyz_xxy_0[i] * pa_z[i] - ta1_x_yyyyz_xxy_1[i] * pc_z[i];

        ta1_x_yyyyzz_xxz_0[i] =
            3.0 * ta1_x_yyzz_xxz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxz_1[i] * fe_0 + ta1_x_yyyzz_xxz_0[i] * pa_y[i] - ta1_x_yyyzz_xxz_1[i] * pc_y[i];

        ta1_x_yyyyzz_xyy_0[i] =
            ta1_x_yyyy_xyy_0[i] * fe_0 - ta1_x_yyyy_xyy_1[i] * fe_0 + ta1_x_yyyyz_xyy_0[i] * pa_z[i] - ta1_x_yyyyz_xyy_1[i] * pc_z[i];

        ta1_x_yyyyzz_xyz_0[i] = 3.0 * ta1_x_yyzz_xyz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xyz_1[i] * fe_0 + ta1_x_yyyzz_xz_0[i] * fe_0 -
                                ta1_x_yyyzz_xz_1[i] * fe_0 + ta1_x_yyyzz_xyz_0[i] * pa_y[i] - ta1_x_yyyzz_xyz_1[i] * pc_y[i];

        ta1_x_yyyyzz_xzz_0[i] =
            3.0 * ta1_x_yyzz_xzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xzz_1[i] * fe_0 + ta1_x_yyyzz_xzz_0[i] * pa_y[i] - ta1_x_yyyzz_xzz_1[i] * pc_y[i];

        ta1_x_yyyyzz_yyy_0[i] =
            ta1_x_yyyy_yyy_0[i] * fe_0 - ta1_x_yyyy_yyy_1[i] * fe_0 + ta1_x_yyyyz_yyy_0[i] * pa_z[i] - ta1_x_yyyyz_yyy_1[i] * pc_z[i];

        ta1_x_yyyyzz_yyz_0[i] = 3.0 * ta1_x_yyzz_yyz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yyz_1[i] * fe_0 + 2.0 * ta1_x_yyyzz_yz_0[i] * fe_0 -
                                2.0 * ta1_x_yyyzz_yz_1[i] * fe_0 + ta1_x_yyyzz_yyz_0[i] * pa_y[i] - ta1_x_yyyzz_yyz_1[i] * pc_y[i];

        ta1_x_yyyyzz_yzz_0[i] = 3.0 * ta1_x_yyzz_yzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yzz_1[i] * fe_0 + ta1_x_yyyzz_zz_0[i] * fe_0 -
                                ta1_x_yyyzz_zz_1[i] * fe_0 + ta1_x_yyyzz_yzz_0[i] * pa_y[i] - ta1_x_yyyzz_yzz_1[i] * pc_y[i];

        ta1_x_yyyyzz_zzz_0[i] =
            3.0 * ta1_x_yyzz_zzz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_zzz_1[i] * fe_0 + ta1_x_yyyzz_zzz_0[i] * pa_y[i] - ta1_x_yyyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 240-250 components of targeted buffer : IF

    auto ta1_x_yyyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 240);

    auto ta1_x_yyyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 241);

    auto ta1_x_yyyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 242);

    auto ta1_x_yyyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 243);

    auto ta1_x_yyyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 244);

    auto ta1_x_yyyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 245);

    auto ta1_x_yyyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 246);

    auto ta1_x_yyyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 247);

    auto ta1_x_yyyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 248);

    auto ta1_x_yyyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 249);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyyz_xxy_0,   \
                             ta1_x_yyyz_xxy_1,   \
                             ta1_x_yyyz_xyy_0,   \
                             ta1_x_yyyz_xyy_1,   \
                             ta1_x_yyyz_yyy_0,   \
                             ta1_x_yyyz_yyy_1,   \
                             ta1_x_yyyzz_xxy_0,  \
                             ta1_x_yyyzz_xxy_1,  \
                             ta1_x_yyyzz_xyy_0,  \
                             ta1_x_yyyzz_xyy_1,  \
                             ta1_x_yyyzz_yyy_0,  \
                             ta1_x_yyyzz_yyy_1,  \
                             ta1_x_yyyzzz_xxx_0, \
                             ta1_x_yyyzzz_xxy_0, \
                             ta1_x_yyyzzz_xxz_0, \
                             ta1_x_yyyzzz_xyy_0, \
                             ta1_x_yyyzzz_xyz_0, \
                             ta1_x_yyyzzz_xzz_0, \
                             ta1_x_yyyzzz_yyy_0, \
                             ta1_x_yyyzzz_yyz_0, \
                             ta1_x_yyyzzz_yzz_0, \
                             ta1_x_yyyzzz_zzz_0, \
                             ta1_x_yyzzz_xxx_0,  \
                             ta1_x_yyzzz_xxx_1,  \
                             ta1_x_yyzzz_xxz_0,  \
                             ta1_x_yyzzz_xxz_1,  \
                             ta1_x_yyzzz_xyz_0,  \
                             ta1_x_yyzzz_xyz_1,  \
                             ta1_x_yyzzz_xz_0,   \
                             ta1_x_yyzzz_xz_1,   \
                             ta1_x_yyzzz_xzz_0,  \
                             ta1_x_yyzzz_xzz_1,  \
                             ta1_x_yyzzz_yyz_0,  \
                             ta1_x_yyzzz_yyz_1,  \
                             ta1_x_yyzzz_yz_0,   \
                             ta1_x_yyzzz_yz_1,   \
                             ta1_x_yyzzz_yzz_0,  \
                             ta1_x_yyzzz_yzz_1,  \
                             ta1_x_yyzzz_zz_0,   \
                             ta1_x_yyzzz_zz_1,   \
                             ta1_x_yyzzz_zzz_0,  \
                             ta1_x_yyzzz_zzz_1,  \
                             ta1_x_yzzz_xxx_0,   \
                             ta1_x_yzzz_xxx_1,   \
                             ta1_x_yzzz_xxz_0,   \
                             ta1_x_yzzz_xxz_1,   \
                             ta1_x_yzzz_xyz_0,   \
                             ta1_x_yzzz_xyz_1,   \
                             ta1_x_yzzz_xzz_0,   \
                             ta1_x_yzzz_xzz_1,   \
                             ta1_x_yzzz_yyz_0,   \
                             ta1_x_yzzz_yyz_1,   \
                             ta1_x_yzzz_yzz_0,   \
                             ta1_x_yzzz_yzz_1,   \
                             ta1_x_yzzz_zzz_0,   \
                             ta1_x_yzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzzz_xxx_0[i] =
            2.0 * ta1_x_yzzz_xxx_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xxx_1[i] * fe_0 + ta1_x_yyzzz_xxx_0[i] * pa_y[i] - ta1_x_yyzzz_xxx_1[i] * pc_y[i];

        ta1_x_yyyzzz_xxy_0[i] =
            2.0 * ta1_x_yyyz_xxy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_xxy_1[i] * fe_0 + ta1_x_yyyzz_xxy_0[i] * pa_z[i] - ta1_x_yyyzz_xxy_1[i] * pc_z[i];

        ta1_x_yyyzzz_xxz_0[i] =
            2.0 * ta1_x_yzzz_xxz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xxz_1[i] * fe_0 + ta1_x_yyzzz_xxz_0[i] * pa_y[i] - ta1_x_yyzzz_xxz_1[i] * pc_y[i];

        ta1_x_yyyzzz_xyy_0[i] =
            2.0 * ta1_x_yyyz_xyy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_xyy_1[i] * fe_0 + ta1_x_yyyzz_xyy_0[i] * pa_z[i] - ta1_x_yyyzz_xyy_1[i] * pc_z[i];

        ta1_x_yyyzzz_xyz_0[i] = 2.0 * ta1_x_yzzz_xyz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xyz_1[i] * fe_0 + ta1_x_yyzzz_xz_0[i] * fe_0 -
                                ta1_x_yyzzz_xz_1[i] * fe_0 + ta1_x_yyzzz_xyz_0[i] * pa_y[i] - ta1_x_yyzzz_xyz_1[i] * pc_y[i];

        ta1_x_yyyzzz_xzz_0[i] =
            2.0 * ta1_x_yzzz_xzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xzz_1[i] * fe_0 + ta1_x_yyzzz_xzz_0[i] * pa_y[i] - ta1_x_yyzzz_xzz_1[i] * pc_y[i];

        ta1_x_yyyzzz_yyy_0[i] =
            2.0 * ta1_x_yyyz_yyy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_yyy_1[i] * fe_0 + ta1_x_yyyzz_yyy_0[i] * pa_z[i] - ta1_x_yyyzz_yyy_1[i] * pc_z[i];

        ta1_x_yyyzzz_yyz_0[i] = 2.0 * ta1_x_yzzz_yyz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_yyz_1[i] * fe_0 + 2.0 * ta1_x_yyzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_x_yyzzz_yz_1[i] * fe_0 + ta1_x_yyzzz_yyz_0[i] * pa_y[i] - ta1_x_yyzzz_yyz_1[i] * pc_y[i];

        ta1_x_yyyzzz_yzz_0[i] = 2.0 * ta1_x_yzzz_yzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_yzz_1[i] * fe_0 + ta1_x_yyzzz_zz_0[i] * fe_0 -
                                ta1_x_yyzzz_zz_1[i] * fe_0 + ta1_x_yyzzz_yzz_0[i] * pa_y[i] - ta1_x_yyzzz_yzz_1[i] * pc_y[i];

        ta1_x_yyyzzz_zzz_0[i] =
            2.0 * ta1_x_yzzz_zzz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_zzz_1[i] * fe_0 + ta1_x_yyzzz_zzz_0[i] * pa_y[i] - ta1_x_yyzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 250-260 components of targeted buffer : IF

    auto ta1_x_yyzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 250);

    auto ta1_x_yyzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 251);

    auto ta1_x_yyzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 252);

    auto ta1_x_yyzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 253);

    auto ta1_x_yyzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 254);

    auto ta1_x_yyzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 255);

    auto ta1_x_yyzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 256);

    auto ta1_x_yyzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 257);

    auto ta1_x_yyzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 258);

    auto ta1_x_yyzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 259);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyzz_xxy_0,   \
                             ta1_x_yyzz_xxy_1,   \
                             ta1_x_yyzz_xyy_0,   \
                             ta1_x_yyzz_xyy_1,   \
                             ta1_x_yyzz_yyy_0,   \
                             ta1_x_yyzz_yyy_1,   \
                             ta1_x_yyzzz_xxy_0,  \
                             ta1_x_yyzzz_xxy_1,  \
                             ta1_x_yyzzz_xyy_0,  \
                             ta1_x_yyzzz_xyy_1,  \
                             ta1_x_yyzzz_yyy_0,  \
                             ta1_x_yyzzz_yyy_1,  \
                             ta1_x_yyzzzz_xxx_0, \
                             ta1_x_yyzzzz_xxy_0, \
                             ta1_x_yyzzzz_xxz_0, \
                             ta1_x_yyzzzz_xyy_0, \
                             ta1_x_yyzzzz_xyz_0, \
                             ta1_x_yyzzzz_xzz_0, \
                             ta1_x_yyzzzz_yyy_0, \
                             ta1_x_yyzzzz_yyz_0, \
                             ta1_x_yyzzzz_yzz_0, \
                             ta1_x_yyzzzz_zzz_0, \
                             ta1_x_yzzzz_xxx_0,  \
                             ta1_x_yzzzz_xxx_1,  \
                             ta1_x_yzzzz_xxz_0,  \
                             ta1_x_yzzzz_xxz_1,  \
                             ta1_x_yzzzz_xyz_0,  \
                             ta1_x_yzzzz_xyz_1,  \
                             ta1_x_yzzzz_xz_0,   \
                             ta1_x_yzzzz_xz_1,   \
                             ta1_x_yzzzz_xzz_0,  \
                             ta1_x_yzzzz_xzz_1,  \
                             ta1_x_yzzzz_yyz_0,  \
                             ta1_x_yzzzz_yyz_1,  \
                             ta1_x_yzzzz_yz_0,   \
                             ta1_x_yzzzz_yz_1,   \
                             ta1_x_yzzzz_yzz_0,  \
                             ta1_x_yzzzz_yzz_1,  \
                             ta1_x_yzzzz_zz_0,   \
                             ta1_x_yzzzz_zz_1,   \
                             ta1_x_yzzzz_zzz_0,  \
                             ta1_x_yzzzz_zzz_1,  \
                             ta1_x_zzzz_xxx_0,   \
                             ta1_x_zzzz_xxx_1,   \
                             ta1_x_zzzz_xxz_0,   \
                             ta1_x_zzzz_xxz_1,   \
                             ta1_x_zzzz_xyz_0,   \
                             ta1_x_zzzz_xyz_1,   \
                             ta1_x_zzzz_xzz_0,   \
                             ta1_x_zzzz_xzz_1,   \
                             ta1_x_zzzz_yyz_0,   \
                             ta1_x_zzzz_yyz_1,   \
                             ta1_x_zzzz_yzz_0,   \
                             ta1_x_zzzz_yzz_1,   \
                             ta1_x_zzzz_zzz_0,   \
                             ta1_x_zzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzzz_xxx_0[i] =
            ta1_x_zzzz_xxx_0[i] * fe_0 - ta1_x_zzzz_xxx_1[i] * fe_0 + ta1_x_yzzzz_xxx_0[i] * pa_y[i] - ta1_x_yzzzz_xxx_1[i] * pc_y[i];

        ta1_x_yyzzzz_xxy_0[i] =
            3.0 * ta1_x_yyzz_xxy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xxy_1[i] * fe_0 + ta1_x_yyzzz_xxy_0[i] * pa_z[i] - ta1_x_yyzzz_xxy_1[i] * pc_z[i];

        ta1_x_yyzzzz_xxz_0[i] =
            ta1_x_zzzz_xxz_0[i] * fe_0 - ta1_x_zzzz_xxz_1[i] * fe_0 + ta1_x_yzzzz_xxz_0[i] * pa_y[i] - ta1_x_yzzzz_xxz_1[i] * pc_y[i];

        ta1_x_yyzzzz_xyy_0[i] =
            3.0 * ta1_x_yyzz_xyy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xyy_1[i] * fe_0 + ta1_x_yyzzz_xyy_0[i] * pa_z[i] - ta1_x_yyzzz_xyy_1[i] * pc_z[i];

        ta1_x_yyzzzz_xyz_0[i] = ta1_x_zzzz_xyz_0[i] * fe_0 - ta1_x_zzzz_xyz_1[i] * fe_0 + ta1_x_yzzzz_xz_0[i] * fe_0 - ta1_x_yzzzz_xz_1[i] * fe_0 +
                                ta1_x_yzzzz_xyz_0[i] * pa_y[i] - ta1_x_yzzzz_xyz_1[i] * pc_y[i];

        ta1_x_yyzzzz_xzz_0[i] =
            ta1_x_zzzz_xzz_0[i] * fe_0 - ta1_x_zzzz_xzz_1[i] * fe_0 + ta1_x_yzzzz_xzz_0[i] * pa_y[i] - ta1_x_yzzzz_xzz_1[i] * pc_y[i];

        ta1_x_yyzzzz_yyy_0[i] =
            3.0 * ta1_x_yyzz_yyy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yyy_1[i] * fe_0 + ta1_x_yyzzz_yyy_0[i] * pa_z[i] - ta1_x_yyzzz_yyy_1[i] * pc_z[i];

        ta1_x_yyzzzz_yyz_0[i] = ta1_x_zzzz_yyz_0[i] * fe_0 - ta1_x_zzzz_yyz_1[i] * fe_0 + 2.0 * ta1_x_yzzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_x_yzzzz_yz_1[i] * fe_0 + ta1_x_yzzzz_yyz_0[i] * pa_y[i] - ta1_x_yzzzz_yyz_1[i] * pc_y[i];

        ta1_x_yyzzzz_yzz_0[i] = ta1_x_zzzz_yzz_0[i] * fe_0 - ta1_x_zzzz_yzz_1[i] * fe_0 + ta1_x_yzzzz_zz_0[i] * fe_0 - ta1_x_yzzzz_zz_1[i] * fe_0 +
                                ta1_x_yzzzz_yzz_0[i] * pa_y[i] - ta1_x_yzzzz_yzz_1[i] * pc_y[i];

        ta1_x_yyzzzz_zzz_0[i] =
            ta1_x_zzzz_zzz_0[i] * fe_0 - ta1_x_zzzz_zzz_1[i] * fe_0 + ta1_x_yzzzz_zzz_0[i] * pa_y[i] - ta1_x_yzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 260-270 components of targeted buffer : IF

    auto ta1_x_yzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 260);

    auto ta1_x_yzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 261);

    auto ta1_x_yzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 262);

    auto ta1_x_yzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 263);

    auto ta1_x_yzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 264);

    auto ta1_x_yzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 265);

    auto ta1_x_yzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 266);

    auto ta1_x_yzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 267);

    auto ta1_x_yzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 268);

    auto ta1_x_yzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 269);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yzzzzz_xxx_0, \
                             ta1_x_yzzzzz_xxy_0, \
                             ta1_x_yzzzzz_xxz_0, \
                             ta1_x_yzzzzz_xyy_0, \
                             ta1_x_yzzzzz_xyz_0, \
                             ta1_x_yzzzzz_xzz_0, \
                             ta1_x_yzzzzz_yyy_0, \
                             ta1_x_yzzzzz_yyz_0, \
                             ta1_x_yzzzzz_yzz_0, \
                             ta1_x_yzzzzz_zzz_0, \
                             ta1_x_zzzzz_xx_0,   \
                             ta1_x_zzzzz_xx_1,   \
                             ta1_x_zzzzz_xxx_0,  \
                             ta1_x_zzzzz_xxx_1,  \
                             ta1_x_zzzzz_xxy_0,  \
                             ta1_x_zzzzz_xxy_1,  \
                             ta1_x_zzzzz_xxz_0,  \
                             ta1_x_zzzzz_xxz_1,  \
                             ta1_x_zzzzz_xy_0,   \
                             ta1_x_zzzzz_xy_1,   \
                             ta1_x_zzzzz_xyy_0,  \
                             ta1_x_zzzzz_xyy_1,  \
                             ta1_x_zzzzz_xyz_0,  \
                             ta1_x_zzzzz_xyz_1,  \
                             ta1_x_zzzzz_xz_0,   \
                             ta1_x_zzzzz_xz_1,   \
                             ta1_x_zzzzz_xzz_0,  \
                             ta1_x_zzzzz_xzz_1,  \
                             ta1_x_zzzzz_yy_0,   \
                             ta1_x_zzzzz_yy_1,   \
                             ta1_x_zzzzz_yyy_0,  \
                             ta1_x_zzzzz_yyy_1,  \
                             ta1_x_zzzzz_yyz_0,  \
                             ta1_x_zzzzz_yyz_1,  \
                             ta1_x_zzzzz_yz_0,   \
                             ta1_x_zzzzz_yz_1,   \
                             ta1_x_zzzzz_yzz_0,  \
                             ta1_x_zzzzz_yzz_1,  \
                             ta1_x_zzzzz_zz_0,   \
                             ta1_x_zzzzz_zz_1,   \
                             ta1_x_zzzzz_zzz_0,  \
                             ta1_x_zzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzzz_xxx_0[i] = ta1_x_zzzzz_xxx_0[i] * pa_y[i] - ta1_x_zzzzz_xxx_1[i] * pc_y[i];

        ta1_x_yzzzzz_xxy_0[i] =
            ta1_x_zzzzz_xx_0[i] * fe_0 - ta1_x_zzzzz_xx_1[i] * fe_0 + ta1_x_zzzzz_xxy_0[i] * pa_y[i] - ta1_x_zzzzz_xxy_1[i] * pc_y[i];

        ta1_x_yzzzzz_xxz_0[i] = ta1_x_zzzzz_xxz_0[i] * pa_y[i] - ta1_x_zzzzz_xxz_1[i] * pc_y[i];

        ta1_x_yzzzzz_xyy_0[i] =
            2.0 * ta1_x_zzzzz_xy_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_xy_1[i] * fe_0 + ta1_x_zzzzz_xyy_0[i] * pa_y[i] - ta1_x_zzzzz_xyy_1[i] * pc_y[i];

        ta1_x_yzzzzz_xyz_0[i] =
            ta1_x_zzzzz_xz_0[i] * fe_0 - ta1_x_zzzzz_xz_1[i] * fe_0 + ta1_x_zzzzz_xyz_0[i] * pa_y[i] - ta1_x_zzzzz_xyz_1[i] * pc_y[i];

        ta1_x_yzzzzz_xzz_0[i] = ta1_x_zzzzz_xzz_0[i] * pa_y[i] - ta1_x_zzzzz_xzz_1[i] * pc_y[i];

        ta1_x_yzzzzz_yyy_0[i] =
            3.0 * ta1_x_zzzzz_yy_0[i] * fe_0 - 3.0 * ta1_x_zzzzz_yy_1[i] * fe_0 + ta1_x_zzzzz_yyy_0[i] * pa_y[i] - ta1_x_zzzzz_yyy_1[i] * pc_y[i];

        ta1_x_yzzzzz_yyz_0[i] =
            2.0 * ta1_x_zzzzz_yz_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_yz_1[i] * fe_0 + ta1_x_zzzzz_yyz_0[i] * pa_y[i] - ta1_x_zzzzz_yyz_1[i] * pc_y[i];

        ta1_x_yzzzzz_yzz_0[i] =
            ta1_x_zzzzz_zz_0[i] * fe_0 - ta1_x_zzzzz_zz_1[i] * fe_0 + ta1_x_zzzzz_yzz_0[i] * pa_y[i] - ta1_x_zzzzz_yzz_1[i] * pc_y[i];

        ta1_x_yzzzzz_zzz_0[i] = ta1_x_zzzzz_zzz_0[i] * pa_y[i] - ta1_x_zzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 270-280 components of targeted buffer : IF

    auto ta1_x_zzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 270);

    auto ta1_x_zzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 271);

    auto ta1_x_zzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 272);

    auto ta1_x_zzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 273);

    auto ta1_x_zzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 274);

    auto ta1_x_zzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 275);

    auto ta1_x_zzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 276);

    auto ta1_x_zzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 277);

    auto ta1_x_zzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 278);

    auto ta1_x_zzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 279);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_zzzz_xxx_0,   \
                             ta1_x_zzzz_xxx_1,   \
                             ta1_x_zzzz_xxy_0,   \
                             ta1_x_zzzz_xxy_1,   \
                             ta1_x_zzzz_xxz_0,   \
                             ta1_x_zzzz_xxz_1,   \
                             ta1_x_zzzz_xyy_0,   \
                             ta1_x_zzzz_xyy_1,   \
                             ta1_x_zzzz_xyz_0,   \
                             ta1_x_zzzz_xyz_1,   \
                             ta1_x_zzzz_xzz_0,   \
                             ta1_x_zzzz_xzz_1,   \
                             ta1_x_zzzz_yyy_0,   \
                             ta1_x_zzzz_yyy_1,   \
                             ta1_x_zzzz_yyz_0,   \
                             ta1_x_zzzz_yyz_1,   \
                             ta1_x_zzzz_yzz_0,   \
                             ta1_x_zzzz_yzz_1,   \
                             ta1_x_zzzz_zzz_0,   \
                             ta1_x_zzzz_zzz_1,   \
                             ta1_x_zzzzz_xx_0,   \
                             ta1_x_zzzzz_xx_1,   \
                             ta1_x_zzzzz_xxx_0,  \
                             ta1_x_zzzzz_xxx_1,  \
                             ta1_x_zzzzz_xxy_0,  \
                             ta1_x_zzzzz_xxy_1,  \
                             ta1_x_zzzzz_xxz_0,  \
                             ta1_x_zzzzz_xxz_1,  \
                             ta1_x_zzzzz_xy_0,   \
                             ta1_x_zzzzz_xy_1,   \
                             ta1_x_zzzzz_xyy_0,  \
                             ta1_x_zzzzz_xyy_1,  \
                             ta1_x_zzzzz_xyz_0,  \
                             ta1_x_zzzzz_xyz_1,  \
                             ta1_x_zzzzz_xz_0,   \
                             ta1_x_zzzzz_xz_1,   \
                             ta1_x_zzzzz_xzz_0,  \
                             ta1_x_zzzzz_xzz_1,  \
                             ta1_x_zzzzz_yy_0,   \
                             ta1_x_zzzzz_yy_1,   \
                             ta1_x_zzzzz_yyy_0,  \
                             ta1_x_zzzzz_yyy_1,  \
                             ta1_x_zzzzz_yyz_0,  \
                             ta1_x_zzzzz_yyz_1,  \
                             ta1_x_zzzzz_yz_0,   \
                             ta1_x_zzzzz_yz_1,   \
                             ta1_x_zzzzz_yzz_0,  \
                             ta1_x_zzzzz_yzz_1,  \
                             ta1_x_zzzzz_zz_0,   \
                             ta1_x_zzzzz_zz_1,   \
                             ta1_x_zzzzz_zzz_0,  \
                             ta1_x_zzzzz_zzz_1,  \
                             ta1_x_zzzzzz_xxx_0, \
                             ta1_x_zzzzzz_xxy_0, \
                             ta1_x_zzzzzz_xxz_0, \
                             ta1_x_zzzzzz_xyy_0, \
                             ta1_x_zzzzzz_xyz_0, \
                             ta1_x_zzzzzz_xzz_0, \
                             ta1_x_zzzzzz_yyy_0, \
                             ta1_x_zzzzzz_yyz_0, \
                             ta1_x_zzzzzz_yzz_0, \
                             ta1_x_zzzzzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzzz_xxx_0[i] =
            5.0 * ta1_x_zzzz_xxx_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxx_1[i] * fe_0 + ta1_x_zzzzz_xxx_0[i] * pa_z[i] - ta1_x_zzzzz_xxx_1[i] * pc_z[i];

        ta1_x_zzzzzz_xxy_0[i] =
            5.0 * ta1_x_zzzz_xxy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxy_1[i] * fe_0 + ta1_x_zzzzz_xxy_0[i] * pa_z[i] - ta1_x_zzzzz_xxy_1[i] * pc_z[i];

        ta1_x_zzzzzz_xxz_0[i] = 5.0 * ta1_x_zzzz_xxz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xxz_1[i] * fe_0 + ta1_x_zzzzz_xx_0[i] * fe_0 -
                                ta1_x_zzzzz_xx_1[i] * fe_0 + ta1_x_zzzzz_xxz_0[i] * pa_z[i] - ta1_x_zzzzz_xxz_1[i] * pc_z[i];

        ta1_x_zzzzzz_xyy_0[i] =
            5.0 * ta1_x_zzzz_xyy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xyy_1[i] * fe_0 + ta1_x_zzzzz_xyy_0[i] * pa_z[i] - ta1_x_zzzzz_xyy_1[i] * pc_z[i];

        ta1_x_zzzzzz_xyz_0[i] = 5.0 * ta1_x_zzzz_xyz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xyz_1[i] * fe_0 + ta1_x_zzzzz_xy_0[i] * fe_0 -
                                ta1_x_zzzzz_xy_1[i] * fe_0 + ta1_x_zzzzz_xyz_0[i] * pa_z[i] - ta1_x_zzzzz_xyz_1[i] * pc_z[i];

        ta1_x_zzzzzz_xzz_0[i] = 5.0 * ta1_x_zzzz_xzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xzz_1[i] * fe_0 + 2.0 * ta1_x_zzzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_x_zzzzz_xz_1[i] * fe_0 + ta1_x_zzzzz_xzz_0[i] * pa_z[i] - ta1_x_zzzzz_xzz_1[i] * pc_z[i];

        ta1_x_zzzzzz_yyy_0[i] =
            5.0 * ta1_x_zzzz_yyy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yyy_1[i] * fe_0 + ta1_x_zzzzz_yyy_0[i] * pa_z[i] - ta1_x_zzzzz_yyy_1[i] * pc_z[i];

        ta1_x_zzzzzz_yyz_0[i] = 5.0 * ta1_x_zzzz_yyz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yyz_1[i] * fe_0 + ta1_x_zzzzz_yy_0[i] * fe_0 -
                                ta1_x_zzzzz_yy_1[i] * fe_0 + ta1_x_zzzzz_yyz_0[i] * pa_z[i] - ta1_x_zzzzz_yyz_1[i] * pc_z[i];

        ta1_x_zzzzzz_yzz_0[i] = 5.0 * ta1_x_zzzz_yzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yzz_1[i] * fe_0 + 2.0 * ta1_x_zzzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_x_zzzzz_yz_1[i] * fe_0 + ta1_x_zzzzz_yzz_0[i] * pa_z[i] - ta1_x_zzzzz_yzz_1[i] * pc_z[i];

        ta1_x_zzzzzz_zzz_0[i] = 5.0 * ta1_x_zzzz_zzz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_zzz_1[i] * fe_0 + 3.0 * ta1_x_zzzzz_zz_0[i] * fe_0 -
                                3.0 * ta1_x_zzzzz_zz_1[i] * fe_0 + ta1_x_zzzzz_zzz_0[i] * pa_z[i] - ta1_x_zzzzz_zzz_1[i] * pc_z[i];
    }

    // Set up 280-290 components of targeted buffer : IF

    auto ta1_y_xxxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 280);

    auto ta1_y_xxxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 281);

    auto ta1_y_xxxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 282);

    auto ta1_y_xxxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 283);

    auto ta1_y_xxxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 284);

    auto ta1_y_xxxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 285);

    auto ta1_y_xxxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 286);

    auto ta1_y_xxxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 287);

    auto ta1_y_xxxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 288);

    auto ta1_y_xxxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 289);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xxxx_xxx_0,   \
                             ta1_y_xxxx_xxx_1,   \
                             ta1_y_xxxx_xxy_0,   \
                             ta1_y_xxxx_xxy_1,   \
                             ta1_y_xxxx_xxz_0,   \
                             ta1_y_xxxx_xxz_1,   \
                             ta1_y_xxxx_xyy_0,   \
                             ta1_y_xxxx_xyy_1,   \
                             ta1_y_xxxx_xyz_0,   \
                             ta1_y_xxxx_xyz_1,   \
                             ta1_y_xxxx_xzz_0,   \
                             ta1_y_xxxx_xzz_1,   \
                             ta1_y_xxxx_yyy_0,   \
                             ta1_y_xxxx_yyy_1,   \
                             ta1_y_xxxx_yyz_0,   \
                             ta1_y_xxxx_yyz_1,   \
                             ta1_y_xxxx_yzz_0,   \
                             ta1_y_xxxx_yzz_1,   \
                             ta1_y_xxxx_zzz_0,   \
                             ta1_y_xxxx_zzz_1,   \
                             ta1_y_xxxxx_xx_0,   \
                             ta1_y_xxxxx_xx_1,   \
                             ta1_y_xxxxx_xxx_0,  \
                             ta1_y_xxxxx_xxx_1,  \
                             ta1_y_xxxxx_xxy_0,  \
                             ta1_y_xxxxx_xxy_1,  \
                             ta1_y_xxxxx_xxz_0,  \
                             ta1_y_xxxxx_xxz_1,  \
                             ta1_y_xxxxx_xy_0,   \
                             ta1_y_xxxxx_xy_1,   \
                             ta1_y_xxxxx_xyy_0,  \
                             ta1_y_xxxxx_xyy_1,  \
                             ta1_y_xxxxx_xyz_0,  \
                             ta1_y_xxxxx_xyz_1,  \
                             ta1_y_xxxxx_xz_0,   \
                             ta1_y_xxxxx_xz_1,   \
                             ta1_y_xxxxx_xzz_0,  \
                             ta1_y_xxxxx_xzz_1,  \
                             ta1_y_xxxxx_yy_0,   \
                             ta1_y_xxxxx_yy_1,   \
                             ta1_y_xxxxx_yyy_0,  \
                             ta1_y_xxxxx_yyy_1,  \
                             ta1_y_xxxxx_yyz_0,  \
                             ta1_y_xxxxx_yyz_1,  \
                             ta1_y_xxxxx_yz_0,   \
                             ta1_y_xxxxx_yz_1,   \
                             ta1_y_xxxxx_yzz_0,  \
                             ta1_y_xxxxx_yzz_1,  \
                             ta1_y_xxxxx_zz_0,   \
                             ta1_y_xxxxx_zz_1,   \
                             ta1_y_xxxxx_zzz_0,  \
                             ta1_y_xxxxx_zzz_1,  \
                             ta1_y_xxxxxx_xxx_0, \
                             ta1_y_xxxxxx_xxy_0, \
                             ta1_y_xxxxxx_xxz_0, \
                             ta1_y_xxxxxx_xyy_0, \
                             ta1_y_xxxxxx_xyz_0, \
                             ta1_y_xxxxxx_xzz_0, \
                             ta1_y_xxxxxx_yyy_0, \
                             ta1_y_xxxxxx_yyz_0, \
                             ta1_y_xxxxxx_yzz_0, \
                             ta1_y_xxxxxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxx_xxx_0[i] = 5.0 * ta1_y_xxxx_xxx_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxx_1[i] * fe_0 + 3.0 * ta1_y_xxxxx_xx_0[i] * fe_0 -
                                3.0 * ta1_y_xxxxx_xx_1[i] * fe_0 + ta1_y_xxxxx_xxx_0[i] * pa_x[i] - ta1_y_xxxxx_xxx_1[i] * pc_x[i];

        ta1_y_xxxxxx_xxy_0[i] = 5.0 * ta1_y_xxxx_xxy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxy_1[i] * fe_0 + 2.0 * ta1_y_xxxxx_xy_0[i] * fe_0 -
                                2.0 * ta1_y_xxxxx_xy_1[i] * fe_0 + ta1_y_xxxxx_xxy_0[i] * pa_x[i] - ta1_y_xxxxx_xxy_1[i] * pc_x[i];

        ta1_y_xxxxxx_xxz_0[i] = 5.0 * ta1_y_xxxx_xxz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxxxx_xz_0[i] * fe_0 -
                                2.0 * ta1_y_xxxxx_xz_1[i] * fe_0 + ta1_y_xxxxx_xxz_0[i] * pa_x[i] - ta1_y_xxxxx_xxz_1[i] * pc_x[i];

        ta1_y_xxxxxx_xyy_0[i] = 5.0 * ta1_y_xxxx_xyy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xyy_1[i] * fe_0 + ta1_y_xxxxx_yy_0[i] * fe_0 -
                                ta1_y_xxxxx_yy_1[i] * fe_0 + ta1_y_xxxxx_xyy_0[i] * pa_x[i] - ta1_y_xxxxx_xyy_1[i] * pc_x[i];

        ta1_y_xxxxxx_xyz_0[i] = 5.0 * ta1_y_xxxx_xyz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xyz_1[i] * fe_0 + ta1_y_xxxxx_yz_0[i] * fe_0 -
                                ta1_y_xxxxx_yz_1[i] * fe_0 + ta1_y_xxxxx_xyz_0[i] * pa_x[i] - ta1_y_xxxxx_xyz_1[i] * pc_x[i];

        ta1_y_xxxxxx_xzz_0[i] = 5.0 * ta1_y_xxxx_xzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xzz_1[i] * fe_0 + ta1_y_xxxxx_zz_0[i] * fe_0 -
                                ta1_y_xxxxx_zz_1[i] * fe_0 + ta1_y_xxxxx_xzz_0[i] * pa_x[i] - ta1_y_xxxxx_xzz_1[i] * pc_x[i];

        ta1_y_xxxxxx_yyy_0[i] =
            5.0 * ta1_y_xxxx_yyy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yyy_1[i] * fe_0 + ta1_y_xxxxx_yyy_0[i] * pa_x[i] - ta1_y_xxxxx_yyy_1[i] * pc_x[i];

        ta1_y_xxxxxx_yyz_0[i] =
            5.0 * ta1_y_xxxx_yyz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yyz_1[i] * fe_0 + ta1_y_xxxxx_yyz_0[i] * pa_x[i] - ta1_y_xxxxx_yyz_1[i] * pc_x[i];

        ta1_y_xxxxxx_yzz_0[i] =
            5.0 * ta1_y_xxxx_yzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yzz_1[i] * fe_0 + ta1_y_xxxxx_yzz_0[i] * pa_x[i] - ta1_y_xxxxx_yzz_1[i] * pc_x[i];

        ta1_y_xxxxxx_zzz_0[i] =
            5.0 * ta1_y_xxxx_zzz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_zzz_1[i] * fe_0 + ta1_y_xxxxx_zzz_0[i] * pa_x[i] - ta1_y_xxxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 290-300 components of targeted buffer : IF

    auto ta1_y_xxxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 290);

    auto ta1_y_xxxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 291);

    auto ta1_y_xxxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 292);

    auto ta1_y_xxxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 293);

    auto ta1_y_xxxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 294);

    auto ta1_y_xxxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 295);

    auto ta1_y_xxxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 296);

    auto ta1_y_xxxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 297);

    auto ta1_y_xxxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 298);

    auto ta1_y_xxxxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 299);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxxxx_xx_0,   \
                             ta1_y_xxxxx_xx_1,   \
                             ta1_y_xxxxx_xxx_0,  \
                             ta1_y_xxxxx_xxx_1,  \
                             ta1_y_xxxxx_xxy_0,  \
                             ta1_y_xxxxx_xxy_1,  \
                             ta1_y_xxxxx_xxz_0,  \
                             ta1_y_xxxxx_xxz_1,  \
                             ta1_y_xxxxx_xy_0,   \
                             ta1_y_xxxxx_xy_1,   \
                             ta1_y_xxxxx_xyy_0,  \
                             ta1_y_xxxxx_xyy_1,  \
                             ta1_y_xxxxx_xyz_0,  \
                             ta1_y_xxxxx_xyz_1,  \
                             ta1_y_xxxxx_xz_0,   \
                             ta1_y_xxxxx_xz_1,   \
                             ta1_y_xxxxx_xzz_0,  \
                             ta1_y_xxxxx_xzz_1,  \
                             ta1_y_xxxxx_zzz_0,  \
                             ta1_y_xxxxx_zzz_1,  \
                             ta1_y_xxxxxy_xxx_0, \
                             ta1_y_xxxxxy_xxy_0, \
                             ta1_y_xxxxxy_xxz_0, \
                             ta1_y_xxxxxy_xyy_0, \
                             ta1_y_xxxxxy_xyz_0, \
                             ta1_y_xxxxxy_xzz_0, \
                             ta1_y_xxxxxy_yyy_0, \
                             ta1_y_xxxxxy_yyz_0, \
                             ta1_y_xxxxxy_yzz_0, \
                             ta1_y_xxxxxy_zzz_0, \
                             ta1_y_xxxxy_yyy_0,  \
                             ta1_y_xxxxy_yyy_1,  \
                             ta1_y_xxxxy_yyz_0,  \
                             ta1_y_xxxxy_yyz_1,  \
                             ta1_y_xxxxy_yzz_0,  \
                             ta1_y_xxxxy_yzz_1,  \
                             ta1_y_xxxy_yyy_0,   \
                             ta1_y_xxxy_yyy_1,   \
                             ta1_y_xxxy_yyz_0,   \
                             ta1_y_xxxy_yyz_1,   \
                             ta1_y_xxxy_yzz_0,   \
                             ta1_y_xxxy_yzz_1,   \
                             ta_xxxxx_xxx_1,     \
                             ta_xxxxx_xxy_1,     \
                             ta_xxxxx_xxz_1,     \
                             ta_xxxxx_xyy_1,     \
                             ta_xxxxx_xyz_1,     \
                             ta_xxxxx_xzz_1,     \
                             ta_xxxxx_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxy_xxx_0[i] = ta_xxxxx_xxx_1[i] + ta1_y_xxxxx_xxx_0[i] * pa_y[i] - ta1_y_xxxxx_xxx_1[i] * pc_y[i];

        ta1_y_xxxxxy_xxy_0[i] = ta1_y_xxxxx_xx_0[i] * fe_0 - ta1_y_xxxxx_xx_1[i] * fe_0 + ta_xxxxx_xxy_1[i] + ta1_y_xxxxx_xxy_0[i] * pa_y[i] -
                                ta1_y_xxxxx_xxy_1[i] * pc_y[i];

        ta1_y_xxxxxy_xxz_0[i] = ta_xxxxx_xxz_1[i] + ta1_y_xxxxx_xxz_0[i] * pa_y[i] - ta1_y_xxxxx_xxz_1[i] * pc_y[i];

        ta1_y_xxxxxy_xyy_0[i] = 2.0 * ta1_y_xxxxx_xy_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xy_1[i] * fe_0 + ta_xxxxx_xyy_1[i] +
                                ta1_y_xxxxx_xyy_0[i] * pa_y[i] - ta1_y_xxxxx_xyy_1[i] * pc_y[i];

        ta1_y_xxxxxy_xyz_0[i] = ta1_y_xxxxx_xz_0[i] * fe_0 - ta1_y_xxxxx_xz_1[i] * fe_0 + ta_xxxxx_xyz_1[i] + ta1_y_xxxxx_xyz_0[i] * pa_y[i] -
                                ta1_y_xxxxx_xyz_1[i] * pc_y[i];

        ta1_y_xxxxxy_xzz_0[i] = ta_xxxxx_xzz_1[i] + ta1_y_xxxxx_xzz_0[i] * pa_y[i] - ta1_y_xxxxx_xzz_1[i] * pc_y[i];

        ta1_y_xxxxxy_yyy_0[i] =
            4.0 * ta1_y_xxxy_yyy_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yyy_1[i] * fe_0 + ta1_y_xxxxy_yyy_0[i] * pa_x[i] - ta1_y_xxxxy_yyy_1[i] * pc_x[i];

        ta1_y_xxxxxy_yyz_0[i] =
            4.0 * ta1_y_xxxy_yyz_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yyz_1[i] * fe_0 + ta1_y_xxxxy_yyz_0[i] * pa_x[i] - ta1_y_xxxxy_yyz_1[i] * pc_x[i];

        ta1_y_xxxxxy_yzz_0[i] =
            4.0 * ta1_y_xxxy_yzz_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yzz_1[i] * fe_0 + ta1_y_xxxxy_yzz_0[i] * pa_x[i] - ta1_y_xxxxy_yzz_1[i] * pc_x[i];

        ta1_y_xxxxxy_zzz_0[i] = ta_xxxxx_zzz_1[i] + ta1_y_xxxxx_zzz_0[i] * pa_y[i] - ta1_y_xxxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 300-310 components of targeted buffer : IF

    auto ta1_y_xxxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 300);

    auto ta1_y_xxxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 301);

    auto ta1_y_xxxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 302);

    auto ta1_y_xxxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 303);

    auto ta1_y_xxxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 304);

    auto ta1_y_xxxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 305);

    auto ta1_y_xxxxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 306);

    auto ta1_y_xxxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 307);

    auto ta1_y_xxxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 308);

    auto ta1_y_xxxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 309);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxxxx_xx_0,   \
                             ta1_y_xxxxx_xx_1,   \
                             ta1_y_xxxxx_xxx_0,  \
                             ta1_y_xxxxx_xxx_1,  \
                             ta1_y_xxxxx_xxy_0,  \
                             ta1_y_xxxxx_xxy_1,  \
                             ta1_y_xxxxx_xxz_0,  \
                             ta1_y_xxxxx_xxz_1,  \
                             ta1_y_xxxxx_xy_0,   \
                             ta1_y_xxxxx_xy_1,   \
                             ta1_y_xxxxx_xyy_0,  \
                             ta1_y_xxxxx_xyy_1,  \
                             ta1_y_xxxxx_xyz_0,  \
                             ta1_y_xxxxx_xyz_1,  \
                             ta1_y_xxxxx_xz_0,   \
                             ta1_y_xxxxx_xz_1,   \
                             ta1_y_xxxxx_xzz_0,  \
                             ta1_y_xxxxx_xzz_1,  \
                             ta1_y_xxxxx_yyy_0,  \
                             ta1_y_xxxxx_yyy_1,  \
                             ta1_y_xxxxxz_xxx_0, \
                             ta1_y_xxxxxz_xxy_0, \
                             ta1_y_xxxxxz_xxz_0, \
                             ta1_y_xxxxxz_xyy_0, \
                             ta1_y_xxxxxz_xyz_0, \
                             ta1_y_xxxxxz_xzz_0, \
                             ta1_y_xxxxxz_yyy_0, \
                             ta1_y_xxxxxz_yyz_0, \
                             ta1_y_xxxxxz_yzz_0, \
                             ta1_y_xxxxxz_zzz_0, \
                             ta1_y_xxxxz_yyz_0,  \
                             ta1_y_xxxxz_yyz_1,  \
                             ta1_y_xxxxz_yzz_0,  \
                             ta1_y_xxxxz_yzz_1,  \
                             ta1_y_xxxxz_zzz_0,  \
                             ta1_y_xxxxz_zzz_1,  \
                             ta1_y_xxxz_yyz_0,   \
                             ta1_y_xxxz_yyz_1,   \
                             ta1_y_xxxz_yzz_0,   \
                             ta1_y_xxxz_yzz_1,   \
                             ta1_y_xxxz_zzz_0,   \
                             ta1_y_xxxz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxz_xxx_0[i] = ta1_y_xxxxx_xxx_0[i] * pa_z[i] - ta1_y_xxxxx_xxx_1[i] * pc_z[i];

        ta1_y_xxxxxz_xxy_0[i] = ta1_y_xxxxx_xxy_0[i] * pa_z[i] - ta1_y_xxxxx_xxy_1[i] * pc_z[i];

        ta1_y_xxxxxz_xxz_0[i] =
            ta1_y_xxxxx_xx_0[i] * fe_0 - ta1_y_xxxxx_xx_1[i] * fe_0 + ta1_y_xxxxx_xxz_0[i] * pa_z[i] - ta1_y_xxxxx_xxz_1[i] * pc_z[i];

        ta1_y_xxxxxz_xyy_0[i] = ta1_y_xxxxx_xyy_0[i] * pa_z[i] - ta1_y_xxxxx_xyy_1[i] * pc_z[i];

        ta1_y_xxxxxz_xyz_0[i] =
            ta1_y_xxxxx_xy_0[i] * fe_0 - ta1_y_xxxxx_xy_1[i] * fe_0 + ta1_y_xxxxx_xyz_0[i] * pa_z[i] - ta1_y_xxxxx_xyz_1[i] * pc_z[i];

        ta1_y_xxxxxz_xzz_0[i] =
            2.0 * ta1_y_xxxxx_xz_0[i] * fe_0 - 2.0 * ta1_y_xxxxx_xz_1[i] * fe_0 + ta1_y_xxxxx_xzz_0[i] * pa_z[i] - ta1_y_xxxxx_xzz_1[i] * pc_z[i];

        ta1_y_xxxxxz_yyy_0[i] = ta1_y_xxxxx_yyy_0[i] * pa_z[i] - ta1_y_xxxxx_yyy_1[i] * pc_z[i];

        ta1_y_xxxxxz_yyz_0[i] =
            4.0 * ta1_y_xxxz_yyz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_yyz_1[i] * fe_0 + ta1_y_xxxxz_yyz_0[i] * pa_x[i] - ta1_y_xxxxz_yyz_1[i] * pc_x[i];

        ta1_y_xxxxxz_yzz_0[i] =
            4.0 * ta1_y_xxxz_yzz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_yzz_1[i] * fe_0 + ta1_y_xxxxz_yzz_0[i] * pa_x[i] - ta1_y_xxxxz_yzz_1[i] * pc_x[i];

        ta1_y_xxxxxz_zzz_0[i] =
            4.0 * ta1_y_xxxz_zzz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_zzz_1[i] * fe_0 + ta1_y_xxxxz_zzz_0[i] * pa_x[i] - ta1_y_xxxxz_zzz_1[i] * pc_x[i];
    }

    // Set up 310-320 components of targeted buffer : IF

    auto ta1_y_xxxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 310);

    auto ta1_y_xxxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 311);

    auto ta1_y_xxxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 312);

    auto ta1_y_xxxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 313);

    auto ta1_y_xxxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 314);

    auto ta1_y_xxxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 315);

    auto ta1_y_xxxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 316);

    auto ta1_y_xxxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 317);

    auto ta1_y_xxxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 318);

    auto ta1_y_xxxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 319);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxxx_xxx_0,   \
                             ta1_y_xxxx_xxx_1,   \
                             ta1_y_xxxx_xxz_0,   \
                             ta1_y_xxxx_xxz_1,   \
                             ta1_y_xxxx_xzz_0,   \
                             ta1_y_xxxx_xzz_1,   \
                             ta1_y_xxxxy_xxx_0,  \
                             ta1_y_xxxxy_xxx_1,  \
                             ta1_y_xxxxy_xxz_0,  \
                             ta1_y_xxxxy_xxz_1,  \
                             ta1_y_xxxxy_xzz_0,  \
                             ta1_y_xxxxy_xzz_1,  \
                             ta1_y_xxxxyy_xxx_0, \
                             ta1_y_xxxxyy_xxy_0, \
                             ta1_y_xxxxyy_xxz_0, \
                             ta1_y_xxxxyy_xyy_0, \
                             ta1_y_xxxxyy_xyz_0, \
                             ta1_y_xxxxyy_xzz_0, \
                             ta1_y_xxxxyy_yyy_0, \
                             ta1_y_xxxxyy_yyz_0, \
                             ta1_y_xxxxyy_yzz_0, \
                             ta1_y_xxxxyy_zzz_0, \
                             ta1_y_xxxyy_xxy_0,  \
                             ta1_y_xxxyy_xxy_1,  \
                             ta1_y_xxxyy_xy_0,   \
                             ta1_y_xxxyy_xy_1,   \
                             ta1_y_xxxyy_xyy_0,  \
                             ta1_y_xxxyy_xyy_1,  \
                             ta1_y_xxxyy_xyz_0,  \
                             ta1_y_xxxyy_xyz_1,  \
                             ta1_y_xxxyy_yy_0,   \
                             ta1_y_xxxyy_yy_1,   \
                             ta1_y_xxxyy_yyy_0,  \
                             ta1_y_xxxyy_yyy_1,  \
                             ta1_y_xxxyy_yyz_0,  \
                             ta1_y_xxxyy_yyz_1,  \
                             ta1_y_xxxyy_yz_0,   \
                             ta1_y_xxxyy_yz_1,   \
                             ta1_y_xxxyy_yzz_0,  \
                             ta1_y_xxxyy_yzz_1,  \
                             ta1_y_xxxyy_zzz_0,  \
                             ta1_y_xxxyy_zzz_1,  \
                             ta1_y_xxyy_xxy_0,   \
                             ta1_y_xxyy_xxy_1,   \
                             ta1_y_xxyy_xyy_0,   \
                             ta1_y_xxyy_xyy_1,   \
                             ta1_y_xxyy_xyz_0,   \
                             ta1_y_xxyy_xyz_1,   \
                             ta1_y_xxyy_yyy_0,   \
                             ta1_y_xxyy_yyy_1,   \
                             ta1_y_xxyy_yyz_0,   \
                             ta1_y_xxyy_yyz_1,   \
                             ta1_y_xxyy_yzz_0,   \
                             ta1_y_xxyy_yzz_1,   \
                             ta1_y_xxyy_zzz_0,   \
                             ta1_y_xxyy_zzz_1,   \
                             ta_xxxxy_xxx_1,     \
                             ta_xxxxy_xxz_1,     \
                             ta_xxxxy_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxyy_xxx_0[i] = ta1_y_xxxx_xxx_0[i] * fe_0 - ta1_y_xxxx_xxx_1[i] * fe_0 + ta_xxxxy_xxx_1[i] + ta1_y_xxxxy_xxx_0[i] * pa_y[i] -
                                ta1_y_xxxxy_xxx_1[i] * pc_y[i];

        ta1_y_xxxxyy_xxy_0[i] = 3.0 * ta1_y_xxyy_xxy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxy_1[i] * fe_0 + 2.0 * ta1_y_xxxyy_xy_0[i] * fe_0 -
                                2.0 * ta1_y_xxxyy_xy_1[i] * fe_0 + ta1_y_xxxyy_xxy_0[i] * pa_x[i] - ta1_y_xxxyy_xxy_1[i] * pc_x[i];

        ta1_y_xxxxyy_xxz_0[i] = ta1_y_xxxx_xxz_0[i] * fe_0 - ta1_y_xxxx_xxz_1[i] * fe_0 + ta_xxxxy_xxz_1[i] + ta1_y_xxxxy_xxz_0[i] * pa_y[i] -
                                ta1_y_xxxxy_xxz_1[i] * pc_y[i];

        ta1_y_xxxxyy_xyy_0[i] = 3.0 * ta1_y_xxyy_xyy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xyy_1[i] * fe_0 + ta1_y_xxxyy_yy_0[i] * fe_0 -
                                ta1_y_xxxyy_yy_1[i] * fe_0 + ta1_y_xxxyy_xyy_0[i] * pa_x[i] - ta1_y_xxxyy_xyy_1[i] * pc_x[i];

        ta1_y_xxxxyy_xyz_0[i] = 3.0 * ta1_y_xxyy_xyz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xyz_1[i] * fe_0 + ta1_y_xxxyy_yz_0[i] * fe_0 -
                                ta1_y_xxxyy_yz_1[i] * fe_0 + ta1_y_xxxyy_xyz_0[i] * pa_x[i] - ta1_y_xxxyy_xyz_1[i] * pc_x[i];

        ta1_y_xxxxyy_xzz_0[i] = ta1_y_xxxx_xzz_0[i] * fe_0 - ta1_y_xxxx_xzz_1[i] * fe_0 + ta_xxxxy_xzz_1[i] + ta1_y_xxxxy_xzz_0[i] * pa_y[i] -
                                ta1_y_xxxxy_xzz_1[i] * pc_y[i];

        ta1_y_xxxxyy_yyy_0[i] =
            3.0 * ta1_y_xxyy_yyy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yyy_1[i] * fe_0 + ta1_y_xxxyy_yyy_0[i] * pa_x[i] - ta1_y_xxxyy_yyy_1[i] * pc_x[i];

        ta1_y_xxxxyy_yyz_0[i] =
            3.0 * ta1_y_xxyy_yyz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yyz_1[i] * fe_0 + ta1_y_xxxyy_yyz_0[i] * pa_x[i] - ta1_y_xxxyy_yyz_1[i] * pc_x[i];

        ta1_y_xxxxyy_yzz_0[i] =
            3.0 * ta1_y_xxyy_yzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yzz_1[i] * fe_0 + ta1_y_xxxyy_yzz_0[i] * pa_x[i] - ta1_y_xxxyy_yzz_1[i] * pc_x[i];

        ta1_y_xxxxyy_zzz_0[i] =
            3.0 * ta1_y_xxyy_zzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_zzz_1[i] * fe_0 + ta1_y_xxxyy_zzz_0[i] * pa_x[i] - ta1_y_xxxyy_zzz_1[i] * pc_x[i];
    }

    // Set up 320-330 components of targeted buffer : IF

    auto ta1_y_xxxxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 320);

    auto ta1_y_xxxxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 321);

    auto ta1_y_xxxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 322);

    auto ta1_y_xxxxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 323);

    auto ta1_y_xxxxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 324);

    auto ta1_y_xxxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 325);

    auto ta1_y_xxxxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 326);

    auto ta1_y_xxxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 327);

    auto ta1_y_xxxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 328);

    auto ta1_y_xxxxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 329);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxxxy_xxx_0,  \
                             ta1_y_xxxxy_xxx_1,  \
                             ta1_y_xxxxy_xxy_0,  \
                             ta1_y_xxxxy_xxy_1,  \
                             ta1_y_xxxxy_xy_0,   \
                             ta1_y_xxxxy_xy_1,   \
                             ta1_y_xxxxy_xyy_0,  \
                             ta1_y_xxxxy_xyy_1,  \
                             ta1_y_xxxxy_xyz_0,  \
                             ta1_y_xxxxy_xyz_1,  \
                             ta1_y_xxxxy_yyy_0,  \
                             ta1_y_xxxxy_yyy_1,  \
                             ta1_y_xxxxyz_xxx_0, \
                             ta1_y_xxxxyz_xxy_0, \
                             ta1_y_xxxxyz_xxz_0, \
                             ta1_y_xxxxyz_xyy_0, \
                             ta1_y_xxxxyz_xyz_0, \
                             ta1_y_xxxxyz_xzz_0, \
                             ta1_y_xxxxyz_yyy_0, \
                             ta1_y_xxxxyz_yyz_0, \
                             ta1_y_xxxxyz_yzz_0, \
                             ta1_y_xxxxyz_zzz_0, \
                             ta1_y_xxxxz_xxz_0,  \
                             ta1_y_xxxxz_xxz_1,  \
                             ta1_y_xxxxz_xzz_0,  \
                             ta1_y_xxxxz_xzz_1,  \
                             ta1_y_xxxxz_zzz_0,  \
                             ta1_y_xxxxz_zzz_1,  \
                             ta1_y_xxxyz_yyz_0,  \
                             ta1_y_xxxyz_yyz_1,  \
                             ta1_y_xxxyz_yzz_0,  \
                             ta1_y_xxxyz_yzz_1,  \
                             ta1_y_xxyz_yyz_0,   \
                             ta1_y_xxyz_yyz_1,   \
                             ta1_y_xxyz_yzz_0,   \
                             ta1_y_xxyz_yzz_1,   \
                             ta_xxxxz_xxz_1,     \
                             ta_xxxxz_xzz_1,     \
                             ta_xxxxz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxyz_xxx_0[i] = ta1_y_xxxxy_xxx_0[i] * pa_z[i] - ta1_y_xxxxy_xxx_1[i] * pc_z[i];

        ta1_y_xxxxyz_xxy_0[i] = ta1_y_xxxxy_xxy_0[i] * pa_z[i] - ta1_y_xxxxy_xxy_1[i] * pc_z[i];

        ta1_y_xxxxyz_xxz_0[i] = ta_xxxxz_xxz_1[i] + ta1_y_xxxxz_xxz_0[i] * pa_y[i] - ta1_y_xxxxz_xxz_1[i] * pc_y[i];

        ta1_y_xxxxyz_xyy_0[i] = ta1_y_xxxxy_xyy_0[i] * pa_z[i] - ta1_y_xxxxy_xyy_1[i] * pc_z[i];

        ta1_y_xxxxyz_xyz_0[i] =
            ta1_y_xxxxy_xy_0[i] * fe_0 - ta1_y_xxxxy_xy_1[i] * fe_0 + ta1_y_xxxxy_xyz_0[i] * pa_z[i] - ta1_y_xxxxy_xyz_1[i] * pc_z[i];

        ta1_y_xxxxyz_xzz_0[i] = ta_xxxxz_xzz_1[i] + ta1_y_xxxxz_xzz_0[i] * pa_y[i] - ta1_y_xxxxz_xzz_1[i] * pc_y[i];

        ta1_y_xxxxyz_yyy_0[i] = ta1_y_xxxxy_yyy_0[i] * pa_z[i] - ta1_y_xxxxy_yyy_1[i] * pc_z[i];

        ta1_y_xxxxyz_yyz_0[i] =
            3.0 * ta1_y_xxyz_yyz_0[i] * fe_0 - 3.0 * ta1_y_xxyz_yyz_1[i] * fe_0 + ta1_y_xxxyz_yyz_0[i] * pa_x[i] - ta1_y_xxxyz_yyz_1[i] * pc_x[i];

        ta1_y_xxxxyz_yzz_0[i] =
            3.0 * ta1_y_xxyz_yzz_0[i] * fe_0 - 3.0 * ta1_y_xxyz_yzz_1[i] * fe_0 + ta1_y_xxxyz_yzz_0[i] * pa_x[i] - ta1_y_xxxyz_yzz_1[i] * pc_x[i];

        ta1_y_xxxxyz_zzz_0[i] = ta_xxxxz_zzz_1[i] + ta1_y_xxxxz_zzz_0[i] * pa_y[i] - ta1_y_xxxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 330-340 components of targeted buffer : IF

    auto ta1_y_xxxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 330);

    auto ta1_y_xxxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 331);

    auto ta1_y_xxxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 332);

    auto ta1_y_xxxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 333);

    auto ta1_y_xxxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 334);

    auto ta1_y_xxxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 335);

    auto ta1_y_xxxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 336);

    auto ta1_y_xxxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 337);

    auto ta1_y_xxxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 338);

    auto ta1_y_xxxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 339);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxxx_xxx_0,   \
                             ta1_y_xxxx_xxx_1,   \
                             ta1_y_xxxx_xxy_0,   \
                             ta1_y_xxxx_xxy_1,   \
                             ta1_y_xxxx_xyy_0,   \
                             ta1_y_xxxx_xyy_1,   \
                             ta1_y_xxxxz_xxx_0,  \
                             ta1_y_xxxxz_xxx_1,  \
                             ta1_y_xxxxz_xxy_0,  \
                             ta1_y_xxxxz_xxy_1,  \
                             ta1_y_xxxxz_xyy_0,  \
                             ta1_y_xxxxz_xyy_1,  \
                             ta1_y_xxxxzz_xxx_0, \
                             ta1_y_xxxxzz_xxy_0, \
                             ta1_y_xxxxzz_xxz_0, \
                             ta1_y_xxxxzz_xyy_0, \
                             ta1_y_xxxxzz_xyz_0, \
                             ta1_y_xxxxzz_xzz_0, \
                             ta1_y_xxxxzz_yyy_0, \
                             ta1_y_xxxxzz_yyz_0, \
                             ta1_y_xxxxzz_yzz_0, \
                             ta1_y_xxxxzz_zzz_0, \
                             ta1_y_xxxzz_xxz_0,  \
                             ta1_y_xxxzz_xxz_1,  \
                             ta1_y_xxxzz_xyz_0,  \
                             ta1_y_xxxzz_xyz_1,  \
                             ta1_y_xxxzz_xz_0,   \
                             ta1_y_xxxzz_xz_1,   \
                             ta1_y_xxxzz_xzz_0,  \
                             ta1_y_xxxzz_xzz_1,  \
                             ta1_y_xxxzz_yyy_0,  \
                             ta1_y_xxxzz_yyy_1,  \
                             ta1_y_xxxzz_yyz_0,  \
                             ta1_y_xxxzz_yyz_1,  \
                             ta1_y_xxxzz_yz_0,   \
                             ta1_y_xxxzz_yz_1,   \
                             ta1_y_xxxzz_yzz_0,  \
                             ta1_y_xxxzz_yzz_1,  \
                             ta1_y_xxxzz_zz_0,   \
                             ta1_y_xxxzz_zz_1,   \
                             ta1_y_xxxzz_zzz_0,  \
                             ta1_y_xxxzz_zzz_1,  \
                             ta1_y_xxzz_xxz_0,   \
                             ta1_y_xxzz_xxz_1,   \
                             ta1_y_xxzz_xyz_0,   \
                             ta1_y_xxzz_xyz_1,   \
                             ta1_y_xxzz_xzz_0,   \
                             ta1_y_xxzz_xzz_1,   \
                             ta1_y_xxzz_yyy_0,   \
                             ta1_y_xxzz_yyy_1,   \
                             ta1_y_xxzz_yyz_0,   \
                             ta1_y_xxzz_yyz_1,   \
                             ta1_y_xxzz_yzz_0,   \
                             ta1_y_xxzz_yzz_1,   \
                             ta1_y_xxzz_zzz_0,   \
                             ta1_y_xxzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxzz_xxx_0[i] =
            ta1_y_xxxx_xxx_0[i] * fe_0 - ta1_y_xxxx_xxx_1[i] * fe_0 + ta1_y_xxxxz_xxx_0[i] * pa_z[i] - ta1_y_xxxxz_xxx_1[i] * pc_z[i];

        ta1_y_xxxxzz_xxy_0[i] =
            ta1_y_xxxx_xxy_0[i] * fe_0 - ta1_y_xxxx_xxy_1[i] * fe_0 + ta1_y_xxxxz_xxy_0[i] * pa_z[i] - ta1_y_xxxxz_xxy_1[i] * pc_z[i];

        ta1_y_xxxxzz_xxz_0[i] = 3.0 * ta1_y_xxzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxxzz_xz_0[i] * fe_0 -
                                2.0 * ta1_y_xxxzz_xz_1[i] * fe_0 + ta1_y_xxxzz_xxz_0[i] * pa_x[i] - ta1_y_xxxzz_xxz_1[i] * pc_x[i];

        ta1_y_xxxxzz_xyy_0[i] =
            ta1_y_xxxx_xyy_0[i] * fe_0 - ta1_y_xxxx_xyy_1[i] * fe_0 + ta1_y_xxxxz_xyy_0[i] * pa_z[i] - ta1_y_xxxxz_xyy_1[i] * pc_z[i];

        ta1_y_xxxxzz_xyz_0[i] = 3.0 * ta1_y_xxzz_xyz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xyz_1[i] * fe_0 + ta1_y_xxxzz_yz_0[i] * fe_0 -
                                ta1_y_xxxzz_yz_1[i] * fe_0 + ta1_y_xxxzz_xyz_0[i] * pa_x[i] - ta1_y_xxxzz_xyz_1[i] * pc_x[i];

        ta1_y_xxxxzz_xzz_0[i] = 3.0 * ta1_y_xxzz_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xzz_1[i] * fe_0 + ta1_y_xxxzz_zz_0[i] * fe_0 -
                                ta1_y_xxxzz_zz_1[i] * fe_0 + ta1_y_xxxzz_xzz_0[i] * pa_x[i] - ta1_y_xxxzz_xzz_1[i] * pc_x[i];

        ta1_y_xxxxzz_yyy_0[i] =
            3.0 * ta1_y_xxzz_yyy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yyy_1[i] * fe_0 + ta1_y_xxxzz_yyy_0[i] * pa_x[i] - ta1_y_xxxzz_yyy_1[i] * pc_x[i];

        ta1_y_xxxxzz_yyz_0[i] =
            3.0 * ta1_y_xxzz_yyz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yyz_1[i] * fe_0 + ta1_y_xxxzz_yyz_0[i] * pa_x[i] - ta1_y_xxxzz_yyz_1[i] * pc_x[i];

        ta1_y_xxxxzz_yzz_0[i] =
            3.0 * ta1_y_xxzz_yzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yzz_1[i] * fe_0 + ta1_y_xxxzz_yzz_0[i] * pa_x[i] - ta1_y_xxxzz_yzz_1[i] * pc_x[i];

        ta1_y_xxxxzz_zzz_0[i] =
            3.0 * ta1_y_xxzz_zzz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_zzz_1[i] * fe_0 + ta1_y_xxxzz_zzz_0[i] * pa_x[i] - ta1_y_xxxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 340-350 components of targeted buffer : IF

    auto ta1_y_xxxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 340);

    auto ta1_y_xxxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 341);

    auto ta1_y_xxxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 342);

    auto ta1_y_xxxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 343);

    auto ta1_y_xxxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 344);

    auto ta1_y_xxxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 345);

    auto ta1_y_xxxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 346);

    auto ta1_y_xxxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 347);

    auto ta1_y_xxxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 348);

    auto ta1_y_xxxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 349);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxxy_xxx_0,   \
                             ta1_y_xxxy_xxx_1,   \
                             ta1_y_xxxy_xxz_0,   \
                             ta1_y_xxxy_xxz_1,   \
                             ta1_y_xxxy_xzz_0,   \
                             ta1_y_xxxy_xzz_1,   \
                             ta1_y_xxxyy_xxx_0,  \
                             ta1_y_xxxyy_xxx_1,  \
                             ta1_y_xxxyy_xxz_0,  \
                             ta1_y_xxxyy_xxz_1,  \
                             ta1_y_xxxyy_xzz_0,  \
                             ta1_y_xxxyy_xzz_1,  \
                             ta1_y_xxxyyy_xxx_0, \
                             ta1_y_xxxyyy_xxy_0, \
                             ta1_y_xxxyyy_xxz_0, \
                             ta1_y_xxxyyy_xyy_0, \
                             ta1_y_xxxyyy_xyz_0, \
                             ta1_y_xxxyyy_xzz_0, \
                             ta1_y_xxxyyy_yyy_0, \
                             ta1_y_xxxyyy_yyz_0, \
                             ta1_y_xxxyyy_yzz_0, \
                             ta1_y_xxxyyy_zzz_0, \
                             ta1_y_xxyyy_xxy_0,  \
                             ta1_y_xxyyy_xxy_1,  \
                             ta1_y_xxyyy_xy_0,   \
                             ta1_y_xxyyy_xy_1,   \
                             ta1_y_xxyyy_xyy_0,  \
                             ta1_y_xxyyy_xyy_1,  \
                             ta1_y_xxyyy_xyz_0,  \
                             ta1_y_xxyyy_xyz_1,  \
                             ta1_y_xxyyy_yy_0,   \
                             ta1_y_xxyyy_yy_1,   \
                             ta1_y_xxyyy_yyy_0,  \
                             ta1_y_xxyyy_yyy_1,  \
                             ta1_y_xxyyy_yyz_0,  \
                             ta1_y_xxyyy_yyz_1,  \
                             ta1_y_xxyyy_yz_0,   \
                             ta1_y_xxyyy_yz_1,   \
                             ta1_y_xxyyy_yzz_0,  \
                             ta1_y_xxyyy_yzz_1,  \
                             ta1_y_xxyyy_zzz_0,  \
                             ta1_y_xxyyy_zzz_1,  \
                             ta1_y_xyyy_xxy_0,   \
                             ta1_y_xyyy_xxy_1,   \
                             ta1_y_xyyy_xyy_0,   \
                             ta1_y_xyyy_xyy_1,   \
                             ta1_y_xyyy_xyz_0,   \
                             ta1_y_xyyy_xyz_1,   \
                             ta1_y_xyyy_yyy_0,   \
                             ta1_y_xyyy_yyy_1,   \
                             ta1_y_xyyy_yyz_0,   \
                             ta1_y_xyyy_yyz_1,   \
                             ta1_y_xyyy_yzz_0,   \
                             ta1_y_xyyy_yzz_1,   \
                             ta1_y_xyyy_zzz_0,   \
                             ta1_y_xyyy_zzz_1,   \
                             ta_xxxyy_xxx_1,     \
                             ta_xxxyy_xxz_1,     \
                             ta_xxxyy_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyy_xxx_0[i] = 2.0 * ta1_y_xxxy_xxx_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xxx_1[i] * fe_0 + ta_xxxyy_xxx_1[i] +
                                ta1_y_xxxyy_xxx_0[i] * pa_y[i] - ta1_y_xxxyy_xxx_1[i] * pc_y[i];

        ta1_y_xxxyyy_xxy_0[i] = 2.0 * ta1_y_xyyy_xxy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xxy_1[i] * fe_0 + 2.0 * ta1_y_xxyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_y_xxyyy_xy_1[i] * fe_0 + ta1_y_xxyyy_xxy_0[i] * pa_x[i] - ta1_y_xxyyy_xxy_1[i] * pc_x[i];

        ta1_y_xxxyyy_xxz_0[i] = 2.0 * ta1_y_xxxy_xxz_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xxz_1[i] * fe_0 + ta_xxxyy_xxz_1[i] +
                                ta1_y_xxxyy_xxz_0[i] * pa_y[i] - ta1_y_xxxyy_xxz_1[i] * pc_y[i];

        ta1_y_xxxyyy_xyy_0[i] = 2.0 * ta1_y_xyyy_xyy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xyy_1[i] * fe_0 + ta1_y_xxyyy_yy_0[i] * fe_0 -
                                ta1_y_xxyyy_yy_1[i] * fe_0 + ta1_y_xxyyy_xyy_0[i] * pa_x[i] - ta1_y_xxyyy_xyy_1[i] * pc_x[i];

        ta1_y_xxxyyy_xyz_0[i] = 2.0 * ta1_y_xyyy_xyz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xyz_1[i] * fe_0 + ta1_y_xxyyy_yz_0[i] * fe_0 -
                                ta1_y_xxyyy_yz_1[i] * fe_0 + ta1_y_xxyyy_xyz_0[i] * pa_x[i] - ta1_y_xxyyy_xyz_1[i] * pc_x[i];

        ta1_y_xxxyyy_xzz_0[i] = 2.0 * ta1_y_xxxy_xzz_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xzz_1[i] * fe_0 + ta_xxxyy_xzz_1[i] +
                                ta1_y_xxxyy_xzz_0[i] * pa_y[i] - ta1_y_xxxyy_xzz_1[i] * pc_y[i];

        ta1_y_xxxyyy_yyy_0[i] =
            2.0 * ta1_y_xyyy_yyy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yyy_1[i] * fe_0 + ta1_y_xxyyy_yyy_0[i] * pa_x[i] - ta1_y_xxyyy_yyy_1[i] * pc_x[i];

        ta1_y_xxxyyy_yyz_0[i] =
            2.0 * ta1_y_xyyy_yyz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yyz_1[i] * fe_0 + ta1_y_xxyyy_yyz_0[i] * pa_x[i] - ta1_y_xxyyy_yyz_1[i] * pc_x[i];

        ta1_y_xxxyyy_yzz_0[i] =
            2.0 * ta1_y_xyyy_yzz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yzz_1[i] * fe_0 + ta1_y_xxyyy_yzz_0[i] * pa_x[i] - ta1_y_xxyyy_yzz_1[i] * pc_x[i];

        ta1_y_xxxyyy_zzz_0[i] =
            2.0 * ta1_y_xyyy_zzz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_zzz_1[i] * fe_0 + ta1_y_xxyyy_zzz_0[i] * pa_x[i] - ta1_y_xxyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 350-360 components of targeted buffer : IF

    auto ta1_y_xxxyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 350);

    auto ta1_y_xxxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 351);

    auto ta1_y_xxxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 352);

    auto ta1_y_xxxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 353);

    auto ta1_y_xxxyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 354);

    auto ta1_y_xxxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 355);

    auto ta1_y_xxxyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 356);

    auto ta1_y_xxxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 357);

    auto ta1_y_xxxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 358);

    auto ta1_y_xxxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 359);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxxyy_xx_0,   \
                             ta1_y_xxxyy_xx_1,   \
                             ta1_y_xxxyy_xxx_0,  \
                             ta1_y_xxxyy_xxx_1,  \
                             ta1_y_xxxyy_xxy_0,  \
                             ta1_y_xxxyy_xxy_1,  \
                             ta1_y_xxxyy_xxz_0,  \
                             ta1_y_xxxyy_xxz_1,  \
                             ta1_y_xxxyy_xy_0,   \
                             ta1_y_xxxyy_xy_1,   \
                             ta1_y_xxxyy_xyy_0,  \
                             ta1_y_xxxyy_xyy_1,  \
                             ta1_y_xxxyy_xyz_0,  \
                             ta1_y_xxxyy_xyz_1,  \
                             ta1_y_xxxyy_xz_0,   \
                             ta1_y_xxxyy_xz_1,   \
                             ta1_y_xxxyy_xzz_0,  \
                             ta1_y_xxxyy_xzz_1,  \
                             ta1_y_xxxyy_yyy_0,  \
                             ta1_y_xxxyy_yyy_1,  \
                             ta1_y_xxxyyz_xxx_0, \
                             ta1_y_xxxyyz_xxy_0, \
                             ta1_y_xxxyyz_xxz_0, \
                             ta1_y_xxxyyz_xyy_0, \
                             ta1_y_xxxyyz_xyz_0, \
                             ta1_y_xxxyyz_xzz_0, \
                             ta1_y_xxxyyz_yyy_0, \
                             ta1_y_xxxyyz_yyz_0, \
                             ta1_y_xxxyyz_yzz_0, \
                             ta1_y_xxxyyz_zzz_0, \
                             ta1_y_xxyyz_yyz_0,  \
                             ta1_y_xxyyz_yyz_1,  \
                             ta1_y_xxyyz_yzz_0,  \
                             ta1_y_xxyyz_yzz_1,  \
                             ta1_y_xxyyz_zzz_0,  \
                             ta1_y_xxyyz_zzz_1,  \
                             ta1_y_xyyz_yyz_0,   \
                             ta1_y_xyyz_yyz_1,   \
                             ta1_y_xyyz_yzz_0,   \
                             ta1_y_xyyz_yzz_1,   \
                             ta1_y_xyyz_zzz_0,   \
                             ta1_y_xyyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyz_xxx_0[i] = ta1_y_xxxyy_xxx_0[i] * pa_z[i] - ta1_y_xxxyy_xxx_1[i] * pc_z[i];

        ta1_y_xxxyyz_xxy_0[i] = ta1_y_xxxyy_xxy_0[i] * pa_z[i] - ta1_y_xxxyy_xxy_1[i] * pc_z[i];

        ta1_y_xxxyyz_xxz_0[i] =
            ta1_y_xxxyy_xx_0[i] * fe_0 - ta1_y_xxxyy_xx_1[i] * fe_0 + ta1_y_xxxyy_xxz_0[i] * pa_z[i] - ta1_y_xxxyy_xxz_1[i] * pc_z[i];

        ta1_y_xxxyyz_xyy_0[i] = ta1_y_xxxyy_xyy_0[i] * pa_z[i] - ta1_y_xxxyy_xyy_1[i] * pc_z[i];

        ta1_y_xxxyyz_xyz_0[i] =
            ta1_y_xxxyy_xy_0[i] * fe_0 - ta1_y_xxxyy_xy_1[i] * fe_0 + ta1_y_xxxyy_xyz_0[i] * pa_z[i] - ta1_y_xxxyy_xyz_1[i] * pc_z[i];

        ta1_y_xxxyyz_xzz_0[i] =
            2.0 * ta1_y_xxxyy_xz_0[i] * fe_0 - 2.0 * ta1_y_xxxyy_xz_1[i] * fe_0 + ta1_y_xxxyy_xzz_0[i] * pa_z[i] - ta1_y_xxxyy_xzz_1[i] * pc_z[i];

        ta1_y_xxxyyz_yyy_0[i] = ta1_y_xxxyy_yyy_0[i] * pa_z[i] - ta1_y_xxxyy_yyy_1[i] * pc_z[i];

        ta1_y_xxxyyz_yyz_0[i] =
            2.0 * ta1_y_xyyz_yyz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_yyz_1[i] * fe_0 + ta1_y_xxyyz_yyz_0[i] * pa_x[i] - ta1_y_xxyyz_yyz_1[i] * pc_x[i];

        ta1_y_xxxyyz_yzz_0[i] =
            2.0 * ta1_y_xyyz_yzz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_yzz_1[i] * fe_0 + ta1_y_xxyyz_yzz_0[i] * pa_x[i] - ta1_y_xxyyz_yzz_1[i] * pc_x[i];

        ta1_y_xxxyyz_zzz_0[i] =
            2.0 * ta1_y_xyyz_zzz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_zzz_1[i] * fe_0 + ta1_y_xxyyz_zzz_0[i] * pa_x[i] - ta1_y_xxyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 360-370 components of targeted buffer : IF

    auto ta1_y_xxxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 360);

    auto ta1_y_xxxyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 361);

    auto ta1_y_xxxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 362);

    auto ta1_y_xxxyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 363);

    auto ta1_y_xxxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 364);

    auto ta1_y_xxxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 365);

    auto ta1_y_xxxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 366);

    auto ta1_y_xxxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 367);

    auto ta1_y_xxxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 368);

    auto ta1_y_xxxyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 369);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxxy_xxy_0,   \
                             ta1_y_xxxy_xxy_1,   \
                             ta1_y_xxxy_xyy_0,   \
                             ta1_y_xxxy_xyy_1,   \
                             ta1_y_xxxyz_xxy_0,  \
                             ta1_y_xxxyz_xxy_1,  \
                             ta1_y_xxxyz_xyy_0,  \
                             ta1_y_xxxyz_xyy_1,  \
                             ta1_y_xxxyzz_xxx_0, \
                             ta1_y_xxxyzz_xxy_0, \
                             ta1_y_xxxyzz_xxz_0, \
                             ta1_y_xxxyzz_xyy_0, \
                             ta1_y_xxxyzz_xyz_0, \
                             ta1_y_xxxyzz_xzz_0, \
                             ta1_y_xxxyzz_yyy_0, \
                             ta1_y_xxxyzz_yyz_0, \
                             ta1_y_xxxyzz_yzz_0, \
                             ta1_y_xxxyzz_zzz_0, \
                             ta1_y_xxxzz_xxx_0,  \
                             ta1_y_xxxzz_xxx_1,  \
                             ta1_y_xxxzz_xxz_0,  \
                             ta1_y_xxxzz_xxz_1,  \
                             ta1_y_xxxzz_xyz_0,  \
                             ta1_y_xxxzz_xyz_1,  \
                             ta1_y_xxxzz_xz_0,   \
                             ta1_y_xxxzz_xz_1,   \
                             ta1_y_xxxzz_xzz_0,  \
                             ta1_y_xxxzz_xzz_1,  \
                             ta1_y_xxxzz_zzz_0,  \
                             ta1_y_xxxzz_zzz_1,  \
                             ta1_y_xxyzz_yyy_0,  \
                             ta1_y_xxyzz_yyy_1,  \
                             ta1_y_xxyzz_yyz_0,  \
                             ta1_y_xxyzz_yyz_1,  \
                             ta1_y_xxyzz_yzz_0,  \
                             ta1_y_xxyzz_yzz_1,  \
                             ta1_y_xyzz_yyy_0,   \
                             ta1_y_xyzz_yyy_1,   \
                             ta1_y_xyzz_yyz_0,   \
                             ta1_y_xyzz_yyz_1,   \
                             ta1_y_xyzz_yzz_0,   \
                             ta1_y_xyzz_yzz_1,   \
                             ta_xxxzz_xxx_1,     \
                             ta_xxxzz_xxz_1,     \
                             ta_xxxzz_xyz_1,     \
                             ta_xxxzz_xzz_1,     \
                             ta_xxxzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyzz_xxx_0[i] = ta_xxxzz_xxx_1[i] + ta1_y_xxxzz_xxx_0[i] * pa_y[i] - ta1_y_xxxzz_xxx_1[i] * pc_y[i];

        ta1_y_xxxyzz_xxy_0[i] =
            ta1_y_xxxy_xxy_0[i] * fe_0 - ta1_y_xxxy_xxy_1[i] * fe_0 + ta1_y_xxxyz_xxy_0[i] * pa_z[i] - ta1_y_xxxyz_xxy_1[i] * pc_z[i];

        ta1_y_xxxyzz_xxz_0[i] = ta_xxxzz_xxz_1[i] + ta1_y_xxxzz_xxz_0[i] * pa_y[i] - ta1_y_xxxzz_xxz_1[i] * pc_y[i];

        ta1_y_xxxyzz_xyy_0[i] =
            ta1_y_xxxy_xyy_0[i] * fe_0 - ta1_y_xxxy_xyy_1[i] * fe_0 + ta1_y_xxxyz_xyy_0[i] * pa_z[i] - ta1_y_xxxyz_xyy_1[i] * pc_z[i];

        ta1_y_xxxyzz_xyz_0[i] = ta1_y_xxxzz_xz_0[i] * fe_0 - ta1_y_xxxzz_xz_1[i] * fe_0 + ta_xxxzz_xyz_1[i] + ta1_y_xxxzz_xyz_0[i] * pa_y[i] -
                                ta1_y_xxxzz_xyz_1[i] * pc_y[i];

        ta1_y_xxxyzz_xzz_0[i] = ta_xxxzz_xzz_1[i] + ta1_y_xxxzz_xzz_0[i] * pa_y[i] - ta1_y_xxxzz_xzz_1[i] * pc_y[i];

        ta1_y_xxxyzz_yyy_0[i] =
            2.0 * ta1_y_xyzz_yyy_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yyy_1[i] * fe_0 + ta1_y_xxyzz_yyy_0[i] * pa_x[i] - ta1_y_xxyzz_yyy_1[i] * pc_x[i];

        ta1_y_xxxyzz_yyz_0[i] =
            2.0 * ta1_y_xyzz_yyz_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yyz_1[i] * fe_0 + ta1_y_xxyzz_yyz_0[i] * pa_x[i] - ta1_y_xxyzz_yyz_1[i] * pc_x[i];

        ta1_y_xxxyzz_yzz_0[i] =
            2.0 * ta1_y_xyzz_yzz_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yzz_1[i] * fe_0 + ta1_y_xxyzz_yzz_0[i] * pa_x[i] - ta1_y_xxyzz_yzz_1[i] * pc_x[i];

        ta1_y_xxxyzz_zzz_0[i] = ta_xxxzz_zzz_1[i] + ta1_y_xxxzz_zzz_0[i] * pa_y[i] - ta1_y_xxxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 370-380 components of targeted buffer : IF

    auto ta1_y_xxxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 370);

    auto ta1_y_xxxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 371);

    auto ta1_y_xxxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 372);

    auto ta1_y_xxxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 373);

    auto ta1_y_xxxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 374);

    auto ta1_y_xxxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 375);

    auto ta1_y_xxxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 376);

    auto ta1_y_xxxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 377);

    auto ta1_y_xxxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 378);

    auto ta1_y_xxxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 379);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxxz_xxx_0,   \
                             ta1_y_xxxz_xxx_1,   \
                             ta1_y_xxxz_xxy_0,   \
                             ta1_y_xxxz_xxy_1,   \
                             ta1_y_xxxz_xyy_0,   \
                             ta1_y_xxxz_xyy_1,   \
                             ta1_y_xxxzz_xxx_0,  \
                             ta1_y_xxxzz_xxx_1,  \
                             ta1_y_xxxzz_xxy_0,  \
                             ta1_y_xxxzz_xxy_1,  \
                             ta1_y_xxxzz_xyy_0,  \
                             ta1_y_xxxzz_xyy_1,  \
                             ta1_y_xxxzzz_xxx_0, \
                             ta1_y_xxxzzz_xxy_0, \
                             ta1_y_xxxzzz_xxz_0, \
                             ta1_y_xxxzzz_xyy_0, \
                             ta1_y_xxxzzz_xyz_0, \
                             ta1_y_xxxzzz_xzz_0, \
                             ta1_y_xxxzzz_yyy_0, \
                             ta1_y_xxxzzz_yyz_0, \
                             ta1_y_xxxzzz_yzz_0, \
                             ta1_y_xxxzzz_zzz_0, \
                             ta1_y_xxzzz_xxz_0,  \
                             ta1_y_xxzzz_xxz_1,  \
                             ta1_y_xxzzz_xyz_0,  \
                             ta1_y_xxzzz_xyz_1,  \
                             ta1_y_xxzzz_xz_0,   \
                             ta1_y_xxzzz_xz_1,   \
                             ta1_y_xxzzz_xzz_0,  \
                             ta1_y_xxzzz_xzz_1,  \
                             ta1_y_xxzzz_yyy_0,  \
                             ta1_y_xxzzz_yyy_1,  \
                             ta1_y_xxzzz_yyz_0,  \
                             ta1_y_xxzzz_yyz_1,  \
                             ta1_y_xxzzz_yz_0,   \
                             ta1_y_xxzzz_yz_1,   \
                             ta1_y_xxzzz_yzz_0,  \
                             ta1_y_xxzzz_yzz_1,  \
                             ta1_y_xxzzz_zz_0,   \
                             ta1_y_xxzzz_zz_1,   \
                             ta1_y_xxzzz_zzz_0,  \
                             ta1_y_xxzzz_zzz_1,  \
                             ta1_y_xzzz_xxz_0,   \
                             ta1_y_xzzz_xxz_1,   \
                             ta1_y_xzzz_xyz_0,   \
                             ta1_y_xzzz_xyz_1,   \
                             ta1_y_xzzz_xzz_0,   \
                             ta1_y_xzzz_xzz_1,   \
                             ta1_y_xzzz_yyy_0,   \
                             ta1_y_xzzz_yyy_1,   \
                             ta1_y_xzzz_yyz_0,   \
                             ta1_y_xzzz_yyz_1,   \
                             ta1_y_xzzz_yzz_0,   \
                             ta1_y_xzzz_yzz_1,   \
                             ta1_y_xzzz_zzz_0,   \
                             ta1_y_xzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzzz_xxx_0[i] =
            2.0 * ta1_y_xxxz_xxx_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xxx_1[i] * fe_0 + ta1_y_xxxzz_xxx_0[i] * pa_z[i] - ta1_y_xxxzz_xxx_1[i] * pc_z[i];

        ta1_y_xxxzzz_xxy_0[i] =
            2.0 * ta1_y_xxxz_xxy_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xxy_1[i] * fe_0 + ta1_y_xxxzz_xxy_0[i] * pa_z[i] - ta1_y_xxxzz_xxy_1[i] * pc_z[i];

        ta1_y_xxxzzz_xxz_0[i] = 2.0 * ta1_y_xzzz_xxz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_y_xxzzz_xz_1[i] * fe_0 + ta1_y_xxzzz_xxz_0[i] * pa_x[i] - ta1_y_xxzzz_xxz_1[i] * pc_x[i];

        ta1_y_xxxzzz_xyy_0[i] =
            2.0 * ta1_y_xxxz_xyy_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xyy_1[i] * fe_0 + ta1_y_xxxzz_xyy_0[i] * pa_z[i] - ta1_y_xxxzz_xyy_1[i] * pc_z[i];

        ta1_y_xxxzzz_xyz_0[i] = 2.0 * ta1_y_xzzz_xyz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xyz_1[i] * fe_0 + ta1_y_xxzzz_yz_0[i] * fe_0 -
                                ta1_y_xxzzz_yz_1[i] * fe_0 + ta1_y_xxzzz_xyz_0[i] * pa_x[i] - ta1_y_xxzzz_xyz_1[i] * pc_x[i];

        ta1_y_xxxzzz_xzz_0[i] = 2.0 * ta1_y_xzzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xzz_1[i] * fe_0 + ta1_y_xxzzz_zz_0[i] * fe_0 -
                                ta1_y_xxzzz_zz_1[i] * fe_0 + ta1_y_xxzzz_xzz_0[i] * pa_x[i] - ta1_y_xxzzz_xzz_1[i] * pc_x[i];

        ta1_y_xxxzzz_yyy_0[i] =
            2.0 * ta1_y_xzzz_yyy_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yyy_1[i] * fe_0 + ta1_y_xxzzz_yyy_0[i] * pa_x[i] - ta1_y_xxzzz_yyy_1[i] * pc_x[i];

        ta1_y_xxxzzz_yyz_0[i] =
            2.0 * ta1_y_xzzz_yyz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yyz_1[i] * fe_0 + ta1_y_xxzzz_yyz_0[i] * pa_x[i] - ta1_y_xxzzz_yyz_1[i] * pc_x[i];

        ta1_y_xxxzzz_yzz_0[i] =
            2.0 * ta1_y_xzzz_yzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yzz_1[i] * fe_0 + ta1_y_xxzzz_yzz_0[i] * pa_x[i] - ta1_y_xxzzz_yzz_1[i] * pc_x[i];

        ta1_y_xxxzzz_zzz_0[i] =
            2.0 * ta1_y_xzzz_zzz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_zzz_1[i] * fe_0 + ta1_y_xxzzz_zzz_0[i] * pa_x[i] - ta1_y_xxzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 380-390 components of targeted buffer : IF

    auto ta1_y_xxyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 380);

    auto ta1_y_xxyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 381);

    auto ta1_y_xxyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 382);

    auto ta1_y_xxyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 383);

    auto ta1_y_xxyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 384);

    auto ta1_y_xxyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 385);

    auto ta1_y_xxyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 386);

    auto ta1_y_xxyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 387);

    auto ta1_y_xxyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 388);

    auto ta1_y_xxyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 389);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxyy_xxx_0,   \
                             ta1_y_xxyy_xxx_1,   \
                             ta1_y_xxyy_xxz_0,   \
                             ta1_y_xxyy_xxz_1,   \
                             ta1_y_xxyy_xzz_0,   \
                             ta1_y_xxyy_xzz_1,   \
                             ta1_y_xxyyy_xxx_0,  \
                             ta1_y_xxyyy_xxx_1,  \
                             ta1_y_xxyyy_xxz_0,  \
                             ta1_y_xxyyy_xxz_1,  \
                             ta1_y_xxyyy_xzz_0,  \
                             ta1_y_xxyyy_xzz_1,  \
                             ta1_y_xxyyyy_xxx_0, \
                             ta1_y_xxyyyy_xxy_0, \
                             ta1_y_xxyyyy_xxz_0, \
                             ta1_y_xxyyyy_xyy_0, \
                             ta1_y_xxyyyy_xyz_0, \
                             ta1_y_xxyyyy_xzz_0, \
                             ta1_y_xxyyyy_yyy_0, \
                             ta1_y_xxyyyy_yyz_0, \
                             ta1_y_xxyyyy_yzz_0, \
                             ta1_y_xxyyyy_zzz_0, \
                             ta1_y_xyyyy_xxy_0,  \
                             ta1_y_xyyyy_xxy_1,  \
                             ta1_y_xyyyy_xy_0,   \
                             ta1_y_xyyyy_xy_1,   \
                             ta1_y_xyyyy_xyy_0,  \
                             ta1_y_xyyyy_xyy_1,  \
                             ta1_y_xyyyy_xyz_0,  \
                             ta1_y_xyyyy_xyz_1,  \
                             ta1_y_xyyyy_yy_0,   \
                             ta1_y_xyyyy_yy_1,   \
                             ta1_y_xyyyy_yyy_0,  \
                             ta1_y_xyyyy_yyy_1,  \
                             ta1_y_xyyyy_yyz_0,  \
                             ta1_y_xyyyy_yyz_1,  \
                             ta1_y_xyyyy_yz_0,   \
                             ta1_y_xyyyy_yz_1,   \
                             ta1_y_xyyyy_yzz_0,  \
                             ta1_y_xyyyy_yzz_1,  \
                             ta1_y_xyyyy_zzz_0,  \
                             ta1_y_xyyyy_zzz_1,  \
                             ta1_y_yyyy_xxy_0,   \
                             ta1_y_yyyy_xxy_1,   \
                             ta1_y_yyyy_xyy_0,   \
                             ta1_y_yyyy_xyy_1,   \
                             ta1_y_yyyy_xyz_0,   \
                             ta1_y_yyyy_xyz_1,   \
                             ta1_y_yyyy_yyy_0,   \
                             ta1_y_yyyy_yyy_1,   \
                             ta1_y_yyyy_yyz_0,   \
                             ta1_y_yyyy_yyz_1,   \
                             ta1_y_yyyy_yzz_0,   \
                             ta1_y_yyyy_yzz_1,   \
                             ta1_y_yyyy_zzz_0,   \
                             ta1_y_yyyy_zzz_1,   \
                             ta_xxyyy_xxx_1,     \
                             ta_xxyyy_xxz_1,     \
                             ta_xxyyy_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyy_xxx_0[i] = 3.0 * ta1_y_xxyy_xxx_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxx_1[i] * fe_0 + ta_xxyyy_xxx_1[i] +
                                ta1_y_xxyyy_xxx_0[i] * pa_y[i] - ta1_y_xxyyy_xxx_1[i] * pc_y[i];

        ta1_y_xxyyyy_xxy_0[i] = ta1_y_yyyy_xxy_0[i] * fe_0 - ta1_y_yyyy_xxy_1[i] * fe_0 + 2.0 * ta1_y_xyyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_y_xyyyy_xy_1[i] * fe_0 + ta1_y_xyyyy_xxy_0[i] * pa_x[i] - ta1_y_xyyyy_xxy_1[i] * pc_x[i];

        ta1_y_xxyyyy_xxz_0[i] = 3.0 * ta1_y_xxyy_xxz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xxz_1[i] * fe_0 + ta_xxyyy_xxz_1[i] +
                                ta1_y_xxyyy_xxz_0[i] * pa_y[i] - ta1_y_xxyyy_xxz_1[i] * pc_y[i];

        ta1_y_xxyyyy_xyy_0[i] = ta1_y_yyyy_xyy_0[i] * fe_0 - ta1_y_yyyy_xyy_1[i] * fe_0 + ta1_y_xyyyy_yy_0[i] * fe_0 - ta1_y_xyyyy_yy_1[i] * fe_0 +
                                ta1_y_xyyyy_xyy_0[i] * pa_x[i] - ta1_y_xyyyy_xyy_1[i] * pc_x[i];

        ta1_y_xxyyyy_xyz_0[i] = ta1_y_yyyy_xyz_0[i] * fe_0 - ta1_y_yyyy_xyz_1[i] * fe_0 + ta1_y_xyyyy_yz_0[i] * fe_0 - ta1_y_xyyyy_yz_1[i] * fe_0 +
                                ta1_y_xyyyy_xyz_0[i] * pa_x[i] - ta1_y_xyyyy_xyz_1[i] * pc_x[i];

        ta1_y_xxyyyy_xzz_0[i] = 3.0 * ta1_y_xxyy_xzz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xzz_1[i] * fe_0 + ta_xxyyy_xzz_1[i] +
                                ta1_y_xxyyy_xzz_0[i] * pa_y[i] - ta1_y_xxyyy_xzz_1[i] * pc_y[i];

        ta1_y_xxyyyy_yyy_0[i] =
            ta1_y_yyyy_yyy_0[i] * fe_0 - ta1_y_yyyy_yyy_1[i] * fe_0 + ta1_y_xyyyy_yyy_0[i] * pa_x[i] - ta1_y_xyyyy_yyy_1[i] * pc_x[i];

        ta1_y_xxyyyy_yyz_0[i] =
            ta1_y_yyyy_yyz_0[i] * fe_0 - ta1_y_yyyy_yyz_1[i] * fe_0 + ta1_y_xyyyy_yyz_0[i] * pa_x[i] - ta1_y_xyyyy_yyz_1[i] * pc_x[i];

        ta1_y_xxyyyy_yzz_0[i] =
            ta1_y_yyyy_yzz_0[i] * fe_0 - ta1_y_yyyy_yzz_1[i] * fe_0 + ta1_y_xyyyy_yzz_0[i] * pa_x[i] - ta1_y_xyyyy_yzz_1[i] * pc_x[i];

        ta1_y_xxyyyy_zzz_0[i] =
            ta1_y_yyyy_zzz_0[i] * fe_0 - ta1_y_yyyy_zzz_1[i] * fe_0 + ta1_y_xyyyy_zzz_0[i] * pa_x[i] - ta1_y_xyyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 390-400 components of targeted buffer : IF

    auto ta1_y_xxyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 390);

    auto ta1_y_xxyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 391);

    auto ta1_y_xxyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 392);

    auto ta1_y_xxyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 393);

    auto ta1_y_xxyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 394);

    auto ta1_y_xxyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 395);

    auto ta1_y_xxyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 396);

    auto ta1_y_xxyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 397);

    auto ta1_y_xxyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 398);

    auto ta1_y_xxyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 399);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxyyy_xx_0,   \
                             ta1_y_xxyyy_xx_1,   \
                             ta1_y_xxyyy_xxx_0,  \
                             ta1_y_xxyyy_xxx_1,  \
                             ta1_y_xxyyy_xxy_0,  \
                             ta1_y_xxyyy_xxy_1,  \
                             ta1_y_xxyyy_xxz_0,  \
                             ta1_y_xxyyy_xxz_1,  \
                             ta1_y_xxyyy_xy_0,   \
                             ta1_y_xxyyy_xy_1,   \
                             ta1_y_xxyyy_xyy_0,  \
                             ta1_y_xxyyy_xyy_1,  \
                             ta1_y_xxyyy_xyz_0,  \
                             ta1_y_xxyyy_xyz_1,  \
                             ta1_y_xxyyy_xz_0,   \
                             ta1_y_xxyyy_xz_1,   \
                             ta1_y_xxyyy_xzz_0,  \
                             ta1_y_xxyyy_xzz_1,  \
                             ta1_y_xxyyy_yyy_0,  \
                             ta1_y_xxyyy_yyy_1,  \
                             ta1_y_xxyyyz_xxx_0, \
                             ta1_y_xxyyyz_xxy_0, \
                             ta1_y_xxyyyz_xxz_0, \
                             ta1_y_xxyyyz_xyy_0, \
                             ta1_y_xxyyyz_xyz_0, \
                             ta1_y_xxyyyz_xzz_0, \
                             ta1_y_xxyyyz_yyy_0, \
                             ta1_y_xxyyyz_yyz_0, \
                             ta1_y_xxyyyz_yzz_0, \
                             ta1_y_xxyyyz_zzz_0, \
                             ta1_y_xyyyz_yyz_0,  \
                             ta1_y_xyyyz_yyz_1,  \
                             ta1_y_xyyyz_yzz_0,  \
                             ta1_y_xyyyz_yzz_1,  \
                             ta1_y_xyyyz_zzz_0,  \
                             ta1_y_xyyyz_zzz_1,  \
                             ta1_y_yyyz_yyz_0,   \
                             ta1_y_yyyz_yyz_1,   \
                             ta1_y_yyyz_yzz_0,   \
                             ta1_y_yyyz_yzz_1,   \
                             ta1_y_yyyz_zzz_0,   \
                             ta1_y_yyyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyz_xxx_0[i] = ta1_y_xxyyy_xxx_0[i] * pa_z[i] - ta1_y_xxyyy_xxx_1[i] * pc_z[i];

        ta1_y_xxyyyz_xxy_0[i] = ta1_y_xxyyy_xxy_0[i] * pa_z[i] - ta1_y_xxyyy_xxy_1[i] * pc_z[i];

        ta1_y_xxyyyz_xxz_0[i] =
            ta1_y_xxyyy_xx_0[i] * fe_0 - ta1_y_xxyyy_xx_1[i] * fe_0 + ta1_y_xxyyy_xxz_0[i] * pa_z[i] - ta1_y_xxyyy_xxz_1[i] * pc_z[i];

        ta1_y_xxyyyz_xyy_0[i] = ta1_y_xxyyy_xyy_0[i] * pa_z[i] - ta1_y_xxyyy_xyy_1[i] * pc_z[i];

        ta1_y_xxyyyz_xyz_0[i] =
            ta1_y_xxyyy_xy_0[i] * fe_0 - ta1_y_xxyyy_xy_1[i] * fe_0 + ta1_y_xxyyy_xyz_0[i] * pa_z[i] - ta1_y_xxyyy_xyz_1[i] * pc_z[i];

        ta1_y_xxyyyz_xzz_0[i] =
            2.0 * ta1_y_xxyyy_xz_0[i] * fe_0 - 2.0 * ta1_y_xxyyy_xz_1[i] * fe_0 + ta1_y_xxyyy_xzz_0[i] * pa_z[i] - ta1_y_xxyyy_xzz_1[i] * pc_z[i];

        ta1_y_xxyyyz_yyy_0[i] = ta1_y_xxyyy_yyy_0[i] * pa_z[i] - ta1_y_xxyyy_yyy_1[i] * pc_z[i];

        ta1_y_xxyyyz_yyz_0[i] =
            ta1_y_yyyz_yyz_0[i] * fe_0 - ta1_y_yyyz_yyz_1[i] * fe_0 + ta1_y_xyyyz_yyz_0[i] * pa_x[i] - ta1_y_xyyyz_yyz_1[i] * pc_x[i];

        ta1_y_xxyyyz_yzz_0[i] =
            ta1_y_yyyz_yzz_0[i] * fe_0 - ta1_y_yyyz_yzz_1[i] * fe_0 + ta1_y_xyyyz_yzz_0[i] * pa_x[i] - ta1_y_xyyyz_yzz_1[i] * pc_x[i];

        ta1_y_xxyyyz_zzz_0[i] =
            ta1_y_yyyz_zzz_0[i] * fe_0 - ta1_y_yyyz_zzz_1[i] * fe_0 + ta1_y_xyyyz_zzz_0[i] * pa_x[i] - ta1_y_xyyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 400-410 components of targeted buffer : IF

    auto ta1_y_xxyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 400);

    auto ta1_y_xxyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 401);

    auto ta1_y_xxyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 402);

    auto ta1_y_xxyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 403);

    auto ta1_y_xxyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 404);

    auto ta1_y_xxyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 405);

    auto ta1_y_xxyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 406);

    auto ta1_y_xxyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 407);

    auto ta1_y_xxyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 408);

    auto ta1_y_xxyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 409);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxyy_xxx_0,   \
                             ta1_y_xxyy_xxx_1,   \
                             ta1_y_xxyy_xxy_0,   \
                             ta1_y_xxyy_xxy_1,   \
                             ta1_y_xxyy_xyy_0,   \
                             ta1_y_xxyy_xyy_1,   \
                             ta1_y_xxyyz_xxx_0,  \
                             ta1_y_xxyyz_xxx_1,  \
                             ta1_y_xxyyz_xxy_0,  \
                             ta1_y_xxyyz_xxy_1,  \
                             ta1_y_xxyyz_xyy_0,  \
                             ta1_y_xxyyz_xyy_1,  \
                             ta1_y_xxyyzz_xxx_0, \
                             ta1_y_xxyyzz_xxy_0, \
                             ta1_y_xxyyzz_xxz_0, \
                             ta1_y_xxyyzz_xyy_0, \
                             ta1_y_xxyyzz_xyz_0, \
                             ta1_y_xxyyzz_xzz_0, \
                             ta1_y_xxyyzz_yyy_0, \
                             ta1_y_xxyyzz_yyz_0, \
                             ta1_y_xxyyzz_yzz_0, \
                             ta1_y_xxyyzz_zzz_0, \
                             ta1_y_xxyzz_xxz_0,  \
                             ta1_y_xxyzz_xxz_1,  \
                             ta1_y_xxyzz_xzz_0,  \
                             ta1_y_xxyzz_xzz_1,  \
                             ta1_y_xxzz_xxz_0,   \
                             ta1_y_xxzz_xxz_1,   \
                             ta1_y_xxzz_xzz_0,   \
                             ta1_y_xxzz_xzz_1,   \
                             ta1_y_xyyzz_xyz_0,  \
                             ta1_y_xyyzz_xyz_1,  \
                             ta1_y_xyyzz_yyy_0,  \
                             ta1_y_xyyzz_yyy_1,  \
                             ta1_y_xyyzz_yyz_0,  \
                             ta1_y_xyyzz_yyz_1,  \
                             ta1_y_xyyzz_yz_0,   \
                             ta1_y_xyyzz_yz_1,   \
                             ta1_y_xyyzz_yzz_0,  \
                             ta1_y_xyyzz_yzz_1,  \
                             ta1_y_xyyzz_zzz_0,  \
                             ta1_y_xyyzz_zzz_1,  \
                             ta1_y_yyzz_xyz_0,   \
                             ta1_y_yyzz_xyz_1,   \
                             ta1_y_yyzz_yyy_0,   \
                             ta1_y_yyzz_yyy_1,   \
                             ta1_y_yyzz_yyz_0,   \
                             ta1_y_yyzz_yyz_1,   \
                             ta1_y_yyzz_yzz_0,   \
                             ta1_y_yyzz_yzz_1,   \
                             ta1_y_yyzz_zzz_0,   \
                             ta1_y_yyzz_zzz_1,   \
                             ta_xxyzz_xxz_1,     \
                             ta_xxyzz_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyzz_xxx_0[i] =
            ta1_y_xxyy_xxx_0[i] * fe_0 - ta1_y_xxyy_xxx_1[i] * fe_0 + ta1_y_xxyyz_xxx_0[i] * pa_z[i] - ta1_y_xxyyz_xxx_1[i] * pc_z[i];

        ta1_y_xxyyzz_xxy_0[i] =
            ta1_y_xxyy_xxy_0[i] * fe_0 - ta1_y_xxyy_xxy_1[i] * fe_0 + ta1_y_xxyyz_xxy_0[i] * pa_z[i] - ta1_y_xxyyz_xxy_1[i] * pc_z[i];

        ta1_y_xxyyzz_xxz_0[i] = ta1_y_xxzz_xxz_0[i] * fe_0 - ta1_y_xxzz_xxz_1[i] * fe_0 + ta_xxyzz_xxz_1[i] + ta1_y_xxyzz_xxz_0[i] * pa_y[i] -
                                ta1_y_xxyzz_xxz_1[i] * pc_y[i];

        ta1_y_xxyyzz_xyy_0[i] =
            ta1_y_xxyy_xyy_0[i] * fe_0 - ta1_y_xxyy_xyy_1[i] * fe_0 + ta1_y_xxyyz_xyy_0[i] * pa_z[i] - ta1_y_xxyyz_xyy_1[i] * pc_z[i];

        ta1_y_xxyyzz_xyz_0[i] = ta1_y_yyzz_xyz_0[i] * fe_0 - ta1_y_yyzz_xyz_1[i] * fe_0 + ta1_y_xyyzz_yz_0[i] * fe_0 - ta1_y_xyyzz_yz_1[i] * fe_0 +
                                ta1_y_xyyzz_xyz_0[i] * pa_x[i] - ta1_y_xyyzz_xyz_1[i] * pc_x[i];

        ta1_y_xxyyzz_xzz_0[i] = ta1_y_xxzz_xzz_0[i] * fe_0 - ta1_y_xxzz_xzz_1[i] * fe_0 + ta_xxyzz_xzz_1[i] + ta1_y_xxyzz_xzz_0[i] * pa_y[i] -
                                ta1_y_xxyzz_xzz_1[i] * pc_y[i];

        ta1_y_xxyyzz_yyy_0[i] =
            ta1_y_yyzz_yyy_0[i] * fe_0 - ta1_y_yyzz_yyy_1[i] * fe_0 + ta1_y_xyyzz_yyy_0[i] * pa_x[i] - ta1_y_xyyzz_yyy_1[i] * pc_x[i];

        ta1_y_xxyyzz_yyz_0[i] =
            ta1_y_yyzz_yyz_0[i] * fe_0 - ta1_y_yyzz_yyz_1[i] * fe_0 + ta1_y_xyyzz_yyz_0[i] * pa_x[i] - ta1_y_xyyzz_yyz_1[i] * pc_x[i];

        ta1_y_xxyyzz_yzz_0[i] =
            ta1_y_yyzz_yzz_0[i] * fe_0 - ta1_y_yyzz_yzz_1[i] * fe_0 + ta1_y_xyyzz_yzz_0[i] * pa_x[i] - ta1_y_xyyzz_yzz_1[i] * pc_x[i];

        ta1_y_xxyyzz_zzz_0[i] =
            ta1_y_yyzz_zzz_0[i] * fe_0 - ta1_y_yyzz_zzz_1[i] * fe_0 + ta1_y_xyyzz_zzz_0[i] * pa_x[i] - ta1_y_xyyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 410-420 components of targeted buffer : IF

    auto ta1_y_xxyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 410);

    auto ta1_y_xxyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 411);

    auto ta1_y_xxyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 412);

    auto ta1_y_xxyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 413);

    auto ta1_y_xxyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 414);

    auto ta1_y_xxyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 415);

    auto ta1_y_xxyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 416);

    auto ta1_y_xxyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 417);

    auto ta1_y_xxyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 418);

    auto ta1_y_xxyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 419);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxyz_xxy_0,   \
                             ta1_y_xxyz_xxy_1,   \
                             ta1_y_xxyz_xyy_0,   \
                             ta1_y_xxyz_xyy_1,   \
                             ta1_y_xxyzz_xxy_0,  \
                             ta1_y_xxyzz_xxy_1,  \
                             ta1_y_xxyzz_xyy_0,  \
                             ta1_y_xxyzz_xyy_1,  \
                             ta1_y_xxyzzz_xxx_0, \
                             ta1_y_xxyzzz_xxy_0, \
                             ta1_y_xxyzzz_xxz_0, \
                             ta1_y_xxyzzz_xyy_0, \
                             ta1_y_xxyzzz_xyz_0, \
                             ta1_y_xxyzzz_xzz_0, \
                             ta1_y_xxyzzz_yyy_0, \
                             ta1_y_xxyzzz_yyz_0, \
                             ta1_y_xxyzzz_yzz_0, \
                             ta1_y_xxyzzz_zzz_0, \
                             ta1_y_xxzzz_xxx_0,  \
                             ta1_y_xxzzz_xxx_1,  \
                             ta1_y_xxzzz_xxz_0,  \
                             ta1_y_xxzzz_xxz_1,  \
                             ta1_y_xxzzz_xyz_0,  \
                             ta1_y_xxzzz_xyz_1,  \
                             ta1_y_xxzzz_xz_0,   \
                             ta1_y_xxzzz_xz_1,   \
                             ta1_y_xxzzz_xzz_0,  \
                             ta1_y_xxzzz_xzz_1,  \
                             ta1_y_xxzzz_zzz_0,  \
                             ta1_y_xxzzz_zzz_1,  \
                             ta1_y_xyzzz_yyy_0,  \
                             ta1_y_xyzzz_yyy_1,  \
                             ta1_y_xyzzz_yyz_0,  \
                             ta1_y_xyzzz_yyz_1,  \
                             ta1_y_xyzzz_yzz_0,  \
                             ta1_y_xyzzz_yzz_1,  \
                             ta1_y_yzzz_yyy_0,   \
                             ta1_y_yzzz_yyy_1,   \
                             ta1_y_yzzz_yyz_0,   \
                             ta1_y_yzzz_yyz_1,   \
                             ta1_y_yzzz_yzz_0,   \
                             ta1_y_yzzz_yzz_1,   \
                             ta_xxzzz_xxx_1,     \
                             ta_xxzzz_xxz_1,     \
                             ta_xxzzz_xyz_1,     \
                             ta_xxzzz_xzz_1,     \
                             ta_xxzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzzz_xxx_0[i] = ta_xxzzz_xxx_1[i] + ta1_y_xxzzz_xxx_0[i] * pa_y[i] - ta1_y_xxzzz_xxx_1[i] * pc_y[i];

        ta1_y_xxyzzz_xxy_0[i] =
            2.0 * ta1_y_xxyz_xxy_0[i] * fe_0 - 2.0 * ta1_y_xxyz_xxy_1[i] * fe_0 + ta1_y_xxyzz_xxy_0[i] * pa_z[i] - ta1_y_xxyzz_xxy_1[i] * pc_z[i];

        ta1_y_xxyzzz_xxz_0[i] = ta_xxzzz_xxz_1[i] + ta1_y_xxzzz_xxz_0[i] * pa_y[i] - ta1_y_xxzzz_xxz_1[i] * pc_y[i];

        ta1_y_xxyzzz_xyy_0[i] =
            2.0 * ta1_y_xxyz_xyy_0[i] * fe_0 - 2.0 * ta1_y_xxyz_xyy_1[i] * fe_0 + ta1_y_xxyzz_xyy_0[i] * pa_z[i] - ta1_y_xxyzz_xyy_1[i] * pc_z[i];

        ta1_y_xxyzzz_xyz_0[i] = ta1_y_xxzzz_xz_0[i] * fe_0 - ta1_y_xxzzz_xz_1[i] * fe_0 + ta_xxzzz_xyz_1[i] + ta1_y_xxzzz_xyz_0[i] * pa_y[i] -
                                ta1_y_xxzzz_xyz_1[i] * pc_y[i];

        ta1_y_xxyzzz_xzz_0[i] = ta_xxzzz_xzz_1[i] + ta1_y_xxzzz_xzz_0[i] * pa_y[i] - ta1_y_xxzzz_xzz_1[i] * pc_y[i];

        ta1_y_xxyzzz_yyy_0[i] =
            ta1_y_yzzz_yyy_0[i] * fe_0 - ta1_y_yzzz_yyy_1[i] * fe_0 + ta1_y_xyzzz_yyy_0[i] * pa_x[i] - ta1_y_xyzzz_yyy_1[i] * pc_x[i];

        ta1_y_xxyzzz_yyz_0[i] =
            ta1_y_yzzz_yyz_0[i] * fe_0 - ta1_y_yzzz_yyz_1[i] * fe_0 + ta1_y_xyzzz_yyz_0[i] * pa_x[i] - ta1_y_xyzzz_yyz_1[i] * pc_x[i];

        ta1_y_xxyzzz_yzz_0[i] =
            ta1_y_yzzz_yzz_0[i] * fe_0 - ta1_y_yzzz_yzz_1[i] * fe_0 + ta1_y_xyzzz_yzz_0[i] * pa_x[i] - ta1_y_xyzzz_yzz_1[i] * pc_x[i];

        ta1_y_xxyzzz_zzz_0[i] = ta_xxzzz_zzz_1[i] + ta1_y_xxzzz_zzz_0[i] * pa_y[i] - ta1_y_xxzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 420-430 components of targeted buffer : IF

    auto ta1_y_xxzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 420);

    auto ta1_y_xxzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 421);

    auto ta1_y_xxzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 422);

    auto ta1_y_xxzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 423);

    auto ta1_y_xxzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 424);

    auto ta1_y_xxzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 425);

    auto ta1_y_xxzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 426);

    auto ta1_y_xxzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 427);

    auto ta1_y_xxzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 428);

    auto ta1_y_xxzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 429);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxzz_xxx_0,   \
                             ta1_y_xxzz_xxx_1,   \
                             ta1_y_xxzz_xxy_0,   \
                             ta1_y_xxzz_xxy_1,   \
                             ta1_y_xxzz_xyy_0,   \
                             ta1_y_xxzz_xyy_1,   \
                             ta1_y_xxzzz_xxx_0,  \
                             ta1_y_xxzzz_xxx_1,  \
                             ta1_y_xxzzz_xxy_0,  \
                             ta1_y_xxzzz_xxy_1,  \
                             ta1_y_xxzzz_xyy_0,  \
                             ta1_y_xxzzz_xyy_1,  \
                             ta1_y_xxzzzz_xxx_0, \
                             ta1_y_xxzzzz_xxy_0, \
                             ta1_y_xxzzzz_xxz_0, \
                             ta1_y_xxzzzz_xyy_0, \
                             ta1_y_xxzzzz_xyz_0, \
                             ta1_y_xxzzzz_xzz_0, \
                             ta1_y_xxzzzz_yyy_0, \
                             ta1_y_xxzzzz_yyz_0, \
                             ta1_y_xxzzzz_yzz_0, \
                             ta1_y_xxzzzz_zzz_0, \
                             ta1_y_xzzzz_xxz_0,  \
                             ta1_y_xzzzz_xxz_1,  \
                             ta1_y_xzzzz_xyz_0,  \
                             ta1_y_xzzzz_xyz_1,  \
                             ta1_y_xzzzz_xz_0,   \
                             ta1_y_xzzzz_xz_1,   \
                             ta1_y_xzzzz_xzz_0,  \
                             ta1_y_xzzzz_xzz_1,  \
                             ta1_y_xzzzz_yyy_0,  \
                             ta1_y_xzzzz_yyy_1,  \
                             ta1_y_xzzzz_yyz_0,  \
                             ta1_y_xzzzz_yyz_1,  \
                             ta1_y_xzzzz_yz_0,   \
                             ta1_y_xzzzz_yz_1,   \
                             ta1_y_xzzzz_yzz_0,  \
                             ta1_y_xzzzz_yzz_1,  \
                             ta1_y_xzzzz_zz_0,   \
                             ta1_y_xzzzz_zz_1,   \
                             ta1_y_xzzzz_zzz_0,  \
                             ta1_y_xzzzz_zzz_1,  \
                             ta1_y_zzzz_xxz_0,   \
                             ta1_y_zzzz_xxz_1,   \
                             ta1_y_zzzz_xyz_0,   \
                             ta1_y_zzzz_xyz_1,   \
                             ta1_y_zzzz_xzz_0,   \
                             ta1_y_zzzz_xzz_1,   \
                             ta1_y_zzzz_yyy_0,   \
                             ta1_y_zzzz_yyy_1,   \
                             ta1_y_zzzz_yyz_0,   \
                             ta1_y_zzzz_yyz_1,   \
                             ta1_y_zzzz_yzz_0,   \
                             ta1_y_zzzz_yzz_1,   \
                             ta1_y_zzzz_zzz_0,   \
                             ta1_y_zzzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzzz_xxx_0[i] =
            3.0 * ta1_y_xxzz_xxx_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxx_1[i] * fe_0 + ta1_y_xxzzz_xxx_0[i] * pa_z[i] - ta1_y_xxzzz_xxx_1[i] * pc_z[i];

        ta1_y_xxzzzz_xxy_0[i] =
            3.0 * ta1_y_xxzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xxy_1[i] * fe_0 + ta1_y_xxzzz_xxy_0[i] * pa_z[i] - ta1_y_xxzzz_xxy_1[i] * pc_z[i];

        ta1_y_xxzzzz_xxz_0[i] = ta1_y_zzzz_xxz_0[i] * fe_0 - ta1_y_zzzz_xxz_1[i] * fe_0 + 2.0 * ta1_y_xzzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_y_xzzzz_xz_1[i] * fe_0 + ta1_y_xzzzz_xxz_0[i] * pa_x[i] - ta1_y_xzzzz_xxz_1[i] * pc_x[i];

        ta1_y_xxzzzz_xyy_0[i] =
            3.0 * ta1_y_xxzz_xyy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xyy_1[i] * fe_0 + ta1_y_xxzzz_xyy_0[i] * pa_z[i] - ta1_y_xxzzz_xyy_1[i] * pc_z[i];

        ta1_y_xxzzzz_xyz_0[i] = ta1_y_zzzz_xyz_0[i] * fe_0 - ta1_y_zzzz_xyz_1[i] * fe_0 + ta1_y_xzzzz_yz_0[i] * fe_0 - ta1_y_xzzzz_yz_1[i] * fe_0 +
                                ta1_y_xzzzz_xyz_0[i] * pa_x[i] - ta1_y_xzzzz_xyz_1[i] * pc_x[i];

        ta1_y_xxzzzz_xzz_0[i] = ta1_y_zzzz_xzz_0[i] * fe_0 - ta1_y_zzzz_xzz_1[i] * fe_0 + ta1_y_xzzzz_zz_0[i] * fe_0 - ta1_y_xzzzz_zz_1[i] * fe_0 +
                                ta1_y_xzzzz_xzz_0[i] * pa_x[i] - ta1_y_xzzzz_xzz_1[i] * pc_x[i];

        ta1_y_xxzzzz_yyy_0[i] =
            ta1_y_zzzz_yyy_0[i] * fe_0 - ta1_y_zzzz_yyy_1[i] * fe_0 + ta1_y_xzzzz_yyy_0[i] * pa_x[i] - ta1_y_xzzzz_yyy_1[i] * pc_x[i];

        ta1_y_xxzzzz_yyz_0[i] =
            ta1_y_zzzz_yyz_0[i] * fe_0 - ta1_y_zzzz_yyz_1[i] * fe_0 + ta1_y_xzzzz_yyz_0[i] * pa_x[i] - ta1_y_xzzzz_yyz_1[i] * pc_x[i];

        ta1_y_xxzzzz_yzz_0[i] =
            ta1_y_zzzz_yzz_0[i] * fe_0 - ta1_y_zzzz_yzz_1[i] * fe_0 + ta1_y_xzzzz_yzz_0[i] * pa_x[i] - ta1_y_xzzzz_yzz_1[i] * pc_x[i];

        ta1_y_xxzzzz_zzz_0[i] =
            ta1_y_zzzz_zzz_0[i] * fe_0 - ta1_y_zzzz_zzz_1[i] * fe_0 + ta1_y_xzzzz_zzz_0[i] * pa_x[i] - ta1_y_xzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 430-440 components of targeted buffer : IF

    auto ta1_y_xyyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 430);

    auto ta1_y_xyyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 431);

    auto ta1_y_xyyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 432);

    auto ta1_y_xyyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 433);

    auto ta1_y_xyyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 434);

    auto ta1_y_xyyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 435);

    auto ta1_y_xyyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 436);

    auto ta1_y_xyyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 437);

    auto ta1_y_xyyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 438);

    auto ta1_y_xyyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 439);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xyyyyy_xxx_0, \
                             ta1_y_xyyyyy_xxy_0, \
                             ta1_y_xyyyyy_xxz_0, \
                             ta1_y_xyyyyy_xyy_0, \
                             ta1_y_xyyyyy_xyz_0, \
                             ta1_y_xyyyyy_xzz_0, \
                             ta1_y_xyyyyy_yyy_0, \
                             ta1_y_xyyyyy_yyz_0, \
                             ta1_y_xyyyyy_yzz_0, \
                             ta1_y_xyyyyy_zzz_0, \
                             ta1_y_yyyyy_xx_0,   \
                             ta1_y_yyyyy_xx_1,   \
                             ta1_y_yyyyy_xxx_0,  \
                             ta1_y_yyyyy_xxx_1,  \
                             ta1_y_yyyyy_xxy_0,  \
                             ta1_y_yyyyy_xxy_1,  \
                             ta1_y_yyyyy_xxz_0,  \
                             ta1_y_yyyyy_xxz_1,  \
                             ta1_y_yyyyy_xy_0,   \
                             ta1_y_yyyyy_xy_1,   \
                             ta1_y_yyyyy_xyy_0,  \
                             ta1_y_yyyyy_xyy_1,  \
                             ta1_y_yyyyy_xyz_0,  \
                             ta1_y_yyyyy_xyz_1,  \
                             ta1_y_yyyyy_xz_0,   \
                             ta1_y_yyyyy_xz_1,   \
                             ta1_y_yyyyy_xzz_0,  \
                             ta1_y_yyyyy_xzz_1,  \
                             ta1_y_yyyyy_yy_0,   \
                             ta1_y_yyyyy_yy_1,   \
                             ta1_y_yyyyy_yyy_0,  \
                             ta1_y_yyyyy_yyy_1,  \
                             ta1_y_yyyyy_yyz_0,  \
                             ta1_y_yyyyy_yyz_1,  \
                             ta1_y_yyyyy_yz_0,   \
                             ta1_y_yyyyy_yz_1,   \
                             ta1_y_yyyyy_yzz_0,  \
                             ta1_y_yyyyy_yzz_1,  \
                             ta1_y_yyyyy_zz_0,   \
                             ta1_y_yyyyy_zz_1,   \
                             ta1_y_yyyyy_zzz_0,  \
                             ta1_y_yyyyy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyyy_xxx_0[i] =
            3.0 * ta1_y_yyyyy_xx_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_xx_1[i] * fe_0 + ta1_y_yyyyy_xxx_0[i] * pa_x[i] - ta1_y_yyyyy_xxx_1[i] * pc_x[i];

        ta1_y_xyyyyy_xxy_0[i] =
            2.0 * ta1_y_yyyyy_xy_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xy_1[i] * fe_0 + ta1_y_yyyyy_xxy_0[i] * pa_x[i] - ta1_y_yyyyy_xxy_1[i] * pc_x[i];

        ta1_y_xyyyyy_xxz_0[i] =
            2.0 * ta1_y_yyyyy_xz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xz_1[i] * fe_0 + ta1_y_yyyyy_xxz_0[i] * pa_x[i] - ta1_y_yyyyy_xxz_1[i] * pc_x[i];

        ta1_y_xyyyyy_xyy_0[i] =
            ta1_y_yyyyy_yy_0[i] * fe_0 - ta1_y_yyyyy_yy_1[i] * fe_0 + ta1_y_yyyyy_xyy_0[i] * pa_x[i] - ta1_y_yyyyy_xyy_1[i] * pc_x[i];

        ta1_y_xyyyyy_xyz_0[i] =
            ta1_y_yyyyy_yz_0[i] * fe_0 - ta1_y_yyyyy_yz_1[i] * fe_0 + ta1_y_yyyyy_xyz_0[i] * pa_x[i] - ta1_y_yyyyy_xyz_1[i] * pc_x[i];

        ta1_y_xyyyyy_xzz_0[i] =
            ta1_y_yyyyy_zz_0[i] * fe_0 - ta1_y_yyyyy_zz_1[i] * fe_0 + ta1_y_yyyyy_xzz_0[i] * pa_x[i] - ta1_y_yyyyy_xzz_1[i] * pc_x[i];

        ta1_y_xyyyyy_yyy_0[i] = ta1_y_yyyyy_yyy_0[i] * pa_x[i] - ta1_y_yyyyy_yyy_1[i] * pc_x[i];

        ta1_y_xyyyyy_yyz_0[i] = ta1_y_yyyyy_yyz_0[i] * pa_x[i] - ta1_y_yyyyy_yyz_1[i] * pc_x[i];

        ta1_y_xyyyyy_yzz_0[i] = ta1_y_yyyyy_yzz_0[i] * pa_x[i] - ta1_y_yyyyy_yzz_1[i] * pc_x[i];

        ta1_y_xyyyyy_zzz_0[i] = ta1_y_yyyyy_zzz_0[i] * pa_x[i] - ta1_y_yyyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 440-450 components of targeted buffer : IF

    auto ta1_y_xyyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 440);

    auto ta1_y_xyyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 441);

    auto ta1_y_xyyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 442);

    auto ta1_y_xyyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 443);

    auto ta1_y_xyyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 444);

    auto ta1_y_xyyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 445);

    auto ta1_y_xyyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 446);

    auto ta1_y_xyyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 447);

    auto ta1_y_xyyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 448);

    auto ta1_y_xyyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 449);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xyyyy_xxx_0,  \
                             ta1_y_xyyyy_xxx_1,  \
                             ta1_y_xyyyy_xxy_0,  \
                             ta1_y_xyyyy_xxy_1,  \
                             ta1_y_xyyyy_xyy_0,  \
                             ta1_y_xyyyy_xyy_1,  \
                             ta1_y_xyyyyz_xxx_0, \
                             ta1_y_xyyyyz_xxy_0, \
                             ta1_y_xyyyyz_xxz_0, \
                             ta1_y_xyyyyz_xyy_0, \
                             ta1_y_xyyyyz_xyz_0, \
                             ta1_y_xyyyyz_xzz_0, \
                             ta1_y_xyyyyz_yyy_0, \
                             ta1_y_xyyyyz_yyz_0, \
                             ta1_y_xyyyyz_yzz_0, \
                             ta1_y_xyyyyz_zzz_0, \
                             ta1_y_yyyyz_xxz_0,  \
                             ta1_y_yyyyz_xxz_1,  \
                             ta1_y_yyyyz_xyz_0,  \
                             ta1_y_yyyyz_xyz_1,  \
                             ta1_y_yyyyz_xz_0,   \
                             ta1_y_yyyyz_xz_1,   \
                             ta1_y_yyyyz_xzz_0,  \
                             ta1_y_yyyyz_xzz_1,  \
                             ta1_y_yyyyz_yyy_0,  \
                             ta1_y_yyyyz_yyy_1,  \
                             ta1_y_yyyyz_yyz_0,  \
                             ta1_y_yyyyz_yyz_1,  \
                             ta1_y_yyyyz_yz_0,   \
                             ta1_y_yyyyz_yz_1,   \
                             ta1_y_yyyyz_yzz_0,  \
                             ta1_y_yyyyz_yzz_1,  \
                             ta1_y_yyyyz_zz_0,   \
                             ta1_y_yyyyz_zz_1,   \
                             ta1_y_yyyyz_zzz_0,  \
                             ta1_y_yyyyz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyyz_xxx_0[i] = ta1_y_xyyyy_xxx_0[i] * pa_z[i] - ta1_y_xyyyy_xxx_1[i] * pc_z[i];

        ta1_y_xyyyyz_xxy_0[i] = ta1_y_xyyyy_xxy_0[i] * pa_z[i] - ta1_y_xyyyy_xxy_1[i] * pc_z[i];

        ta1_y_xyyyyz_xxz_0[i] =
            2.0 * ta1_y_yyyyz_xz_0[i] * fe_0 - 2.0 * ta1_y_yyyyz_xz_1[i] * fe_0 + ta1_y_yyyyz_xxz_0[i] * pa_x[i] - ta1_y_yyyyz_xxz_1[i] * pc_x[i];

        ta1_y_xyyyyz_xyy_0[i] = ta1_y_xyyyy_xyy_0[i] * pa_z[i] - ta1_y_xyyyy_xyy_1[i] * pc_z[i];

        ta1_y_xyyyyz_xyz_0[i] =
            ta1_y_yyyyz_yz_0[i] * fe_0 - ta1_y_yyyyz_yz_1[i] * fe_0 + ta1_y_yyyyz_xyz_0[i] * pa_x[i] - ta1_y_yyyyz_xyz_1[i] * pc_x[i];

        ta1_y_xyyyyz_xzz_0[i] =
            ta1_y_yyyyz_zz_0[i] * fe_0 - ta1_y_yyyyz_zz_1[i] * fe_0 + ta1_y_yyyyz_xzz_0[i] * pa_x[i] - ta1_y_yyyyz_xzz_1[i] * pc_x[i];

        ta1_y_xyyyyz_yyy_0[i] = ta1_y_yyyyz_yyy_0[i] * pa_x[i] - ta1_y_yyyyz_yyy_1[i] * pc_x[i];

        ta1_y_xyyyyz_yyz_0[i] = ta1_y_yyyyz_yyz_0[i] * pa_x[i] - ta1_y_yyyyz_yyz_1[i] * pc_x[i];

        ta1_y_xyyyyz_yzz_0[i] = ta1_y_yyyyz_yzz_0[i] * pa_x[i] - ta1_y_yyyyz_yzz_1[i] * pc_x[i];

        ta1_y_xyyyyz_zzz_0[i] = ta1_y_yyyyz_zzz_0[i] * pa_x[i] - ta1_y_yyyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 450-460 components of targeted buffer : IF

    auto ta1_y_xyyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 450);

    auto ta1_y_xyyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 451);

    auto ta1_y_xyyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 452);

    auto ta1_y_xyyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 453);

    auto ta1_y_xyyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 454);

    auto ta1_y_xyyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 455);

    auto ta1_y_xyyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 456);

    auto ta1_y_xyyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 457);

    auto ta1_y_xyyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 458);

    auto ta1_y_xyyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 459);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xyyyzz_xxx_0, \
                             ta1_y_xyyyzz_xxy_0, \
                             ta1_y_xyyyzz_xxz_0, \
                             ta1_y_xyyyzz_xyy_0, \
                             ta1_y_xyyyzz_xyz_0, \
                             ta1_y_xyyyzz_xzz_0, \
                             ta1_y_xyyyzz_yyy_0, \
                             ta1_y_xyyyzz_yyz_0, \
                             ta1_y_xyyyzz_yzz_0, \
                             ta1_y_xyyyzz_zzz_0, \
                             ta1_y_yyyzz_xx_0,   \
                             ta1_y_yyyzz_xx_1,   \
                             ta1_y_yyyzz_xxx_0,  \
                             ta1_y_yyyzz_xxx_1,  \
                             ta1_y_yyyzz_xxy_0,  \
                             ta1_y_yyyzz_xxy_1,  \
                             ta1_y_yyyzz_xxz_0,  \
                             ta1_y_yyyzz_xxz_1,  \
                             ta1_y_yyyzz_xy_0,   \
                             ta1_y_yyyzz_xy_1,   \
                             ta1_y_yyyzz_xyy_0,  \
                             ta1_y_yyyzz_xyy_1,  \
                             ta1_y_yyyzz_xyz_0,  \
                             ta1_y_yyyzz_xyz_1,  \
                             ta1_y_yyyzz_xz_0,   \
                             ta1_y_yyyzz_xz_1,   \
                             ta1_y_yyyzz_xzz_0,  \
                             ta1_y_yyyzz_xzz_1,  \
                             ta1_y_yyyzz_yy_0,   \
                             ta1_y_yyyzz_yy_1,   \
                             ta1_y_yyyzz_yyy_0,  \
                             ta1_y_yyyzz_yyy_1,  \
                             ta1_y_yyyzz_yyz_0,  \
                             ta1_y_yyyzz_yyz_1,  \
                             ta1_y_yyyzz_yz_0,   \
                             ta1_y_yyyzz_yz_1,   \
                             ta1_y_yyyzz_yzz_0,  \
                             ta1_y_yyyzz_yzz_1,  \
                             ta1_y_yyyzz_zz_0,   \
                             ta1_y_yyyzz_zz_1,   \
                             ta1_y_yyyzz_zzz_0,  \
                             ta1_y_yyyzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyzz_xxx_0[i] =
            3.0 * ta1_y_yyyzz_xx_0[i] * fe_0 - 3.0 * ta1_y_yyyzz_xx_1[i] * fe_0 + ta1_y_yyyzz_xxx_0[i] * pa_x[i] - ta1_y_yyyzz_xxx_1[i] * pc_x[i];

        ta1_y_xyyyzz_xxy_0[i] =
            2.0 * ta1_y_yyyzz_xy_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_xy_1[i] * fe_0 + ta1_y_yyyzz_xxy_0[i] * pa_x[i] - ta1_y_yyyzz_xxy_1[i] * pc_x[i];

        ta1_y_xyyyzz_xxz_0[i] =
            2.0 * ta1_y_yyyzz_xz_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_xz_1[i] * fe_0 + ta1_y_yyyzz_xxz_0[i] * pa_x[i] - ta1_y_yyyzz_xxz_1[i] * pc_x[i];

        ta1_y_xyyyzz_xyy_0[i] =
            ta1_y_yyyzz_yy_0[i] * fe_0 - ta1_y_yyyzz_yy_1[i] * fe_0 + ta1_y_yyyzz_xyy_0[i] * pa_x[i] - ta1_y_yyyzz_xyy_1[i] * pc_x[i];

        ta1_y_xyyyzz_xyz_0[i] =
            ta1_y_yyyzz_yz_0[i] * fe_0 - ta1_y_yyyzz_yz_1[i] * fe_0 + ta1_y_yyyzz_xyz_0[i] * pa_x[i] - ta1_y_yyyzz_xyz_1[i] * pc_x[i];

        ta1_y_xyyyzz_xzz_0[i] =
            ta1_y_yyyzz_zz_0[i] * fe_0 - ta1_y_yyyzz_zz_1[i] * fe_0 + ta1_y_yyyzz_xzz_0[i] * pa_x[i] - ta1_y_yyyzz_xzz_1[i] * pc_x[i];

        ta1_y_xyyyzz_yyy_0[i] = ta1_y_yyyzz_yyy_0[i] * pa_x[i] - ta1_y_yyyzz_yyy_1[i] * pc_x[i];

        ta1_y_xyyyzz_yyz_0[i] = ta1_y_yyyzz_yyz_0[i] * pa_x[i] - ta1_y_yyyzz_yyz_1[i] * pc_x[i];

        ta1_y_xyyyzz_yzz_0[i] = ta1_y_yyyzz_yzz_0[i] * pa_x[i] - ta1_y_yyyzz_yzz_1[i] * pc_x[i];

        ta1_y_xyyyzz_zzz_0[i] = ta1_y_yyyzz_zzz_0[i] * pa_x[i] - ta1_y_yyyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 460-470 components of targeted buffer : IF

    auto ta1_y_xyyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 460);

    auto ta1_y_xyyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 461);

    auto ta1_y_xyyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 462);

    auto ta1_y_xyyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 463);

    auto ta1_y_xyyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 464);

    auto ta1_y_xyyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 465);

    auto ta1_y_xyyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 466);

    auto ta1_y_xyyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 467);

    auto ta1_y_xyyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 468);

    auto ta1_y_xyyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 469);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xyyzzz_xxx_0, \
                             ta1_y_xyyzzz_xxy_0, \
                             ta1_y_xyyzzz_xxz_0, \
                             ta1_y_xyyzzz_xyy_0, \
                             ta1_y_xyyzzz_xyz_0, \
                             ta1_y_xyyzzz_xzz_0, \
                             ta1_y_xyyzzz_yyy_0, \
                             ta1_y_xyyzzz_yyz_0, \
                             ta1_y_xyyzzz_yzz_0, \
                             ta1_y_xyyzzz_zzz_0, \
                             ta1_y_yyzzz_xx_0,   \
                             ta1_y_yyzzz_xx_1,   \
                             ta1_y_yyzzz_xxx_0,  \
                             ta1_y_yyzzz_xxx_1,  \
                             ta1_y_yyzzz_xxy_0,  \
                             ta1_y_yyzzz_xxy_1,  \
                             ta1_y_yyzzz_xxz_0,  \
                             ta1_y_yyzzz_xxz_1,  \
                             ta1_y_yyzzz_xy_0,   \
                             ta1_y_yyzzz_xy_1,   \
                             ta1_y_yyzzz_xyy_0,  \
                             ta1_y_yyzzz_xyy_1,  \
                             ta1_y_yyzzz_xyz_0,  \
                             ta1_y_yyzzz_xyz_1,  \
                             ta1_y_yyzzz_xz_0,   \
                             ta1_y_yyzzz_xz_1,   \
                             ta1_y_yyzzz_xzz_0,  \
                             ta1_y_yyzzz_xzz_1,  \
                             ta1_y_yyzzz_yy_0,   \
                             ta1_y_yyzzz_yy_1,   \
                             ta1_y_yyzzz_yyy_0,  \
                             ta1_y_yyzzz_yyy_1,  \
                             ta1_y_yyzzz_yyz_0,  \
                             ta1_y_yyzzz_yyz_1,  \
                             ta1_y_yyzzz_yz_0,   \
                             ta1_y_yyzzz_yz_1,   \
                             ta1_y_yyzzz_yzz_0,  \
                             ta1_y_yyzzz_yzz_1,  \
                             ta1_y_yyzzz_zz_0,   \
                             ta1_y_yyzzz_zz_1,   \
                             ta1_y_yyzzz_zzz_0,  \
                             ta1_y_yyzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzzz_xxx_0[i] =
            3.0 * ta1_y_yyzzz_xx_0[i] * fe_0 - 3.0 * ta1_y_yyzzz_xx_1[i] * fe_0 + ta1_y_yyzzz_xxx_0[i] * pa_x[i] - ta1_y_yyzzz_xxx_1[i] * pc_x[i];

        ta1_y_xyyzzz_xxy_0[i] =
            2.0 * ta1_y_yyzzz_xy_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_xy_1[i] * fe_0 + ta1_y_yyzzz_xxy_0[i] * pa_x[i] - ta1_y_yyzzz_xxy_1[i] * pc_x[i];

        ta1_y_xyyzzz_xxz_0[i] =
            2.0 * ta1_y_yyzzz_xz_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_xz_1[i] * fe_0 + ta1_y_yyzzz_xxz_0[i] * pa_x[i] - ta1_y_yyzzz_xxz_1[i] * pc_x[i];

        ta1_y_xyyzzz_xyy_0[i] =
            ta1_y_yyzzz_yy_0[i] * fe_0 - ta1_y_yyzzz_yy_1[i] * fe_0 + ta1_y_yyzzz_xyy_0[i] * pa_x[i] - ta1_y_yyzzz_xyy_1[i] * pc_x[i];

        ta1_y_xyyzzz_xyz_0[i] =
            ta1_y_yyzzz_yz_0[i] * fe_0 - ta1_y_yyzzz_yz_1[i] * fe_0 + ta1_y_yyzzz_xyz_0[i] * pa_x[i] - ta1_y_yyzzz_xyz_1[i] * pc_x[i];

        ta1_y_xyyzzz_xzz_0[i] =
            ta1_y_yyzzz_zz_0[i] * fe_0 - ta1_y_yyzzz_zz_1[i] * fe_0 + ta1_y_yyzzz_xzz_0[i] * pa_x[i] - ta1_y_yyzzz_xzz_1[i] * pc_x[i];

        ta1_y_xyyzzz_yyy_0[i] = ta1_y_yyzzz_yyy_0[i] * pa_x[i] - ta1_y_yyzzz_yyy_1[i] * pc_x[i];

        ta1_y_xyyzzz_yyz_0[i] = ta1_y_yyzzz_yyz_0[i] * pa_x[i] - ta1_y_yyzzz_yyz_1[i] * pc_x[i];

        ta1_y_xyyzzz_yzz_0[i] = ta1_y_yyzzz_yzz_0[i] * pa_x[i] - ta1_y_yyzzz_yzz_1[i] * pc_x[i];

        ta1_y_xyyzzz_zzz_0[i] = ta1_y_yyzzz_zzz_0[i] * pa_x[i] - ta1_y_yyzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 470-480 components of targeted buffer : IF

    auto ta1_y_xyzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 470);

    auto ta1_y_xyzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 471);

    auto ta1_y_xyzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 472);

    auto ta1_y_xyzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 473);

    auto ta1_y_xyzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 474);

    auto ta1_y_xyzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 475);

    auto ta1_y_xyzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 476);

    auto ta1_y_xyzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 477);

    auto ta1_y_xyzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 478);

    auto ta1_y_xyzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 479);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xyzzzz_xxx_0, \
                             ta1_y_xyzzzz_xxy_0, \
                             ta1_y_xyzzzz_xxz_0, \
                             ta1_y_xyzzzz_xyy_0, \
                             ta1_y_xyzzzz_xyz_0, \
                             ta1_y_xyzzzz_xzz_0, \
                             ta1_y_xyzzzz_yyy_0, \
                             ta1_y_xyzzzz_yyz_0, \
                             ta1_y_xyzzzz_yzz_0, \
                             ta1_y_xyzzzz_zzz_0, \
                             ta1_y_xzzzz_xxx_0,  \
                             ta1_y_xzzzz_xxx_1,  \
                             ta1_y_xzzzz_xxz_0,  \
                             ta1_y_xzzzz_xxz_1,  \
                             ta1_y_xzzzz_xzz_0,  \
                             ta1_y_xzzzz_xzz_1,  \
                             ta1_y_yzzzz_xxy_0,  \
                             ta1_y_yzzzz_xxy_1,  \
                             ta1_y_yzzzz_xy_0,   \
                             ta1_y_yzzzz_xy_1,   \
                             ta1_y_yzzzz_xyy_0,  \
                             ta1_y_yzzzz_xyy_1,  \
                             ta1_y_yzzzz_xyz_0,  \
                             ta1_y_yzzzz_xyz_1,  \
                             ta1_y_yzzzz_yy_0,   \
                             ta1_y_yzzzz_yy_1,   \
                             ta1_y_yzzzz_yyy_0,  \
                             ta1_y_yzzzz_yyy_1,  \
                             ta1_y_yzzzz_yyz_0,  \
                             ta1_y_yzzzz_yyz_1,  \
                             ta1_y_yzzzz_yz_0,   \
                             ta1_y_yzzzz_yz_1,   \
                             ta1_y_yzzzz_yzz_0,  \
                             ta1_y_yzzzz_yzz_1,  \
                             ta1_y_yzzzz_zzz_0,  \
                             ta1_y_yzzzz_zzz_1,  \
                             ta_xzzzz_xxx_1,     \
                             ta_xzzzz_xxz_1,     \
                             ta_xzzzz_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzzzz_xxx_0[i] = ta_xzzzz_xxx_1[i] + ta1_y_xzzzz_xxx_0[i] * pa_y[i] - ta1_y_xzzzz_xxx_1[i] * pc_y[i];

        ta1_y_xyzzzz_xxy_0[i] =
            2.0 * ta1_y_yzzzz_xy_0[i] * fe_0 - 2.0 * ta1_y_yzzzz_xy_1[i] * fe_0 + ta1_y_yzzzz_xxy_0[i] * pa_x[i] - ta1_y_yzzzz_xxy_1[i] * pc_x[i];

        ta1_y_xyzzzz_xxz_0[i] = ta_xzzzz_xxz_1[i] + ta1_y_xzzzz_xxz_0[i] * pa_y[i] - ta1_y_xzzzz_xxz_1[i] * pc_y[i];

        ta1_y_xyzzzz_xyy_0[i] =
            ta1_y_yzzzz_yy_0[i] * fe_0 - ta1_y_yzzzz_yy_1[i] * fe_0 + ta1_y_yzzzz_xyy_0[i] * pa_x[i] - ta1_y_yzzzz_xyy_1[i] * pc_x[i];

        ta1_y_xyzzzz_xyz_0[i] =
            ta1_y_yzzzz_yz_0[i] * fe_0 - ta1_y_yzzzz_yz_1[i] * fe_0 + ta1_y_yzzzz_xyz_0[i] * pa_x[i] - ta1_y_yzzzz_xyz_1[i] * pc_x[i];

        ta1_y_xyzzzz_xzz_0[i] = ta_xzzzz_xzz_1[i] + ta1_y_xzzzz_xzz_0[i] * pa_y[i] - ta1_y_xzzzz_xzz_1[i] * pc_y[i];

        ta1_y_xyzzzz_yyy_0[i] = ta1_y_yzzzz_yyy_0[i] * pa_x[i] - ta1_y_yzzzz_yyy_1[i] * pc_x[i];

        ta1_y_xyzzzz_yyz_0[i] = ta1_y_yzzzz_yyz_0[i] * pa_x[i] - ta1_y_yzzzz_yyz_1[i] * pc_x[i];

        ta1_y_xyzzzz_yzz_0[i] = ta1_y_yzzzz_yzz_0[i] * pa_x[i] - ta1_y_yzzzz_yzz_1[i] * pc_x[i];

        ta1_y_xyzzzz_zzz_0[i] = ta1_y_yzzzz_zzz_0[i] * pa_x[i] - ta1_y_yzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 480-490 components of targeted buffer : IF

    auto ta1_y_xzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 480);

    auto ta1_y_xzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 481);

    auto ta1_y_xzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 482);

    auto ta1_y_xzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 483);

    auto ta1_y_xzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 484);

    auto ta1_y_xzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 485);

    auto ta1_y_xzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 486);

    auto ta1_y_xzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 487);

    auto ta1_y_xzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 488);

    auto ta1_y_xzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 489);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xzzzzz_xxx_0, \
                             ta1_y_xzzzzz_xxy_0, \
                             ta1_y_xzzzzz_xxz_0, \
                             ta1_y_xzzzzz_xyy_0, \
                             ta1_y_xzzzzz_xyz_0, \
                             ta1_y_xzzzzz_xzz_0, \
                             ta1_y_xzzzzz_yyy_0, \
                             ta1_y_xzzzzz_yyz_0, \
                             ta1_y_xzzzzz_yzz_0, \
                             ta1_y_xzzzzz_zzz_0, \
                             ta1_y_zzzzz_xx_0,   \
                             ta1_y_zzzzz_xx_1,   \
                             ta1_y_zzzzz_xxx_0,  \
                             ta1_y_zzzzz_xxx_1,  \
                             ta1_y_zzzzz_xxy_0,  \
                             ta1_y_zzzzz_xxy_1,  \
                             ta1_y_zzzzz_xxz_0,  \
                             ta1_y_zzzzz_xxz_1,  \
                             ta1_y_zzzzz_xy_0,   \
                             ta1_y_zzzzz_xy_1,   \
                             ta1_y_zzzzz_xyy_0,  \
                             ta1_y_zzzzz_xyy_1,  \
                             ta1_y_zzzzz_xyz_0,  \
                             ta1_y_zzzzz_xyz_1,  \
                             ta1_y_zzzzz_xz_0,   \
                             ta1_y_zzzzz_xz_1,   \
                             ta1_y_zzzzz_xzz_0,  \
                             ta1_y_zzzzz_xzz_1,  \
                             ta1_y_zzzzz_yy_0,   \
                             ta1_y_zzzzz_yy_1,   \
                             ta1_y_zzzzz_yyy_0,  \
                             ta1_y_zzzzz_yyy_1,  \
                             ta1_y_zzzzz_yyz_0,  \
                             ta1_y_zzzzz_yyz_1,  \
                             ta1_y_zzzzz_yz_0,   \
                             ta1_y_zzzzz_yz_1,   \
                             ta1_y_zzzzz_yzz_0,  \
                             ta1_y_zzzzz_yzz_1,  \
                             ta1_y_zzzzz_zz_0,   \
                             ta1_y_zzzzz_zz_1,   \
                             ta1_y_zzzzz_zzz_0,  \
                             ta1_y_zzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzzz_xxx_0[i] =
            3.0 * ta1_y_zzzzz_xx_0[i] * fe_0 - 3.0 * ta1_y_zzzzz_xx_1[i] * fe_0 + ta1_y_zzzzz_xxx_0[i] * pa_x[i] - ta1_y_zzzzz_xxx_1[i] * pc_x[i];

        ta1_y_xzzzzz_xxy_0[i] =
            2.0 * ta1_y_zzzzz_xy_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xy_1[i] * fe_0 + ta1_y_zzzzz_xxy_0[i] * pa_x[i] - ta1_y_zzzzz_xxy_1[i] * pc_x[i];

        ta1_y_xzzzzz_xxz_0[i] =
            2.0 * ta1_y_zzzzz_xz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_xz_1[i] * fe_0 + ta1_y_zzzzz_xxz_0[i] * pa_x[i] - ta1_y_zzzzz_xxz_1[i] * pc_x[i];

        ta1_y_xzzzzz_xyy_0[i] =
            ta1_y_zzzzz_yy_0[i] * fe_0 - ta1_y_zzzzz_yy_1[i] * fe_0 + ta1_y_zzzzz_xyy_0[i] * pa_x[i] - ta1_y_zzzzz_xyy_1[i] * pc_x[i];

        ta1_y_xzzzzz_xyz_0[i] =
            ta1_y_zzzzz_yz_0[i] * fe_0 - ta1_y_zzzzz_yz_1[i] * fe_0 + ta1_y_zzzzz_xyz_0[i] * pa_x[i] - ta1_y_zzzzz_xyz_1[i] * pc_x[i];

        ta1_y_xzzzzz_xzz_0[i] =
            ta1_y_zzzzz_zz_0[i] * fe_0 - ta1_y_zzzzz_zz_1[i] * fe_0 + ta1_y_zzzzz_xzz_0[i] * pa_x[i] - ta1_y_zzzzz_xzz_1[i] * pc_x[i];

        ta1_y_xzzzzz_yyy_0[i] = ta1_y_zzzzz_yyy_0[i] * pa_x[i] - ta1_y_zzzzz_yyy_1[i] * pc_x[i];

        ta1_y_xzzzzz_yyz_0[i] = ta1_y_zzzzz_yyz_0[i] * pa_x[i] - ta1_y_zzzzz_yyz_1[i] * pc_x[i];

        ta1_y_xzzzzz_yzz_0[i] = ta1_y_zzzzz_yzz_0[i] * pa_x[i] - ta1_y_zzzzz_yzz_1[i] * pc_x[i];

        ta1_y_xzzzzz_zzz_0[i] = ta1_y_zzzzz_zzz_0[i] * pa_x[i] - ta1_y_zzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 490-500 components of targeted buffer : IF

    auto ta1_y_yyyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 490);

    auto ta1_y_yyyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 491);

    auto ta1_y_yyyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 492);

    auto ta1_y_yyyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 493);

    auto ta1_y_yyyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 494);

    auto ta1_y_yyyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 495);

    auto ta1_y_yyyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 496);

    auto ta1_y_yyyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 497);

    auto ta1_y_yyyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 498);

    auto ta1_y_yyyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 499);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_y_yyyy_xxx_0,   \
                             ta1_y_yyyy_xxx_1,   \
                             ta1_y_yyyy_xxy_0,   \
                             ta1_y_yyyy_xxy_1,   \
                             ta1_y_yyyy_xxz_0,   \
                             ta1_y_yyyy_xxz_1,   \
                             ta1_y_yyyy_xyy_0,   \
                             ta1_y_yyyy_xyy_1,   \
                             ta1_y_yyyy_xyz_0,   \
                             ta1_y_yyyy_xyz_1,   \
                             ta1_y_yyyy_xzz_0,   \
                             ta1_y_yyyy_xzz_1,   \
                             ta1_y_yyyy_yyy_0,   \
                             ta1_y_yyyy_yyy_1,   \
                             ta1_y_yyyy_yyz_0,   \
                             ta1_y_yyyy_yyz_1,   \
                             ta1_y_yyyy_yzz_0,   \
                             ta1_y_yyyy_yzz_1,   \
                             ta1_y_yyyy_zzz_0,   \
                             ta1_y_yyyy_zzz_1,   \
                             ta1_y_yyyyy_xx_0,   \
                             ta1_y_yyyyy_xx_1,   \
                             ta1_y_yyyyy_xxx_0,  \
                             ta1_y_yyyyy_xxx_1,  \
                             ta1_y_yyyyy_xxy_0,  \
                             ta1_y_yyyyy_xxy_1,  \
                             ta1_y_yyyyy_xxz_0,  \
                             ta1_y_yyyyy_xxz_1,  \
                             ta1_y_yyyyy_xy_0,   \
                             ta1_y_yyyyy_xy_1,   \
                             ta1_y_yyyyy_xyy_0,  \
                             ta1_y_yyyyy_xyy_1,  \
                             ta1_y_yyyyy_xyz_0,  \
                             ta1_y_yyyyy_xyz_1,  \
                             ta1_y_yyyyy_xz_0,   \
                             ta1_y_yyyyy_xz_1,   \
                             ta1_y_yyyyy_xzz_0,  \
                             ta1_y_yyyyy_xzz_1,  \
                             ta1_y_yyyyy_yy_0,   \
                             ta1_y_yyyyy_yy_1,   \
                             ta1_y_yyyyy_yyy_0,  \
                             ta1_y_yyyyy_yyy_1,  \
                             ta1_y_yyyyy_yyz_0,  \
                             ta1_y_yyyyy_yyz_1,  \
                             ta1_y_yyyyy_yz_0,   \
                             ta1_y_yyyyy_yz_1,   \
                             ta1_y_yyyyy_yzz_0,  \
                             ta1_y_yyyyy_yzz_1,  \
                             ta1_y_yyyyy_zz_0,   \
                             ta1_y_yyyyy_zz_1,   \
                             ta1_y_yyyyy_zzz_0,  \
                             ta1_y_yyyyy_zzz_1,  \
                             ta1_y_yyyyyy_xxx_0, \
                             ta1_y_yyyyyy_xxy_0, \
                             ta1_y_yyyyyy_xxz_0, \
                             ta1_y_yyyyyy_xyy_0, \
                             ta1_y_yyyyyy_xyz_0, \
                             ta1_y_yyyyyy_xzz_0, \
                             ta1_y_yyyyyy_yyy_0, \
                             ta1_y_yyyyyy_yyz_0, \
                             ta1_y_yyyyyy_yzz_0, \
                             ta1_y_yyyyyy_zzz_0, \
                             ta_yyyyy_xxx_1,     \
                             ta_yyyyy_xxy_1,     \
                             ta_yyyyy_xxz_1,     \
                             ta_yyyyy_xyy_1,     \
                             ta_yyyyy_xyz_1,     \
                             ta_yyyyy_xzz_1,     \
                             ta_yyyyy_yyy_1,     \
                             ta_yyyyy_yyz_1,     \
                             ta_yyyyy_yzz_1,     \
                             ta_yyyyy_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyy_xxx_0[i] = 5.0 * ta1_y_yyyy_xxx_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxx_1[i] * fe_0 + ta_yyyyy_xxx_1[i] +
                                ta1_y_yyyyy_xxx_0[i] * pa_y[i] - ta1_y_yyyyy_xxx_1[i] * pc_y[i];

        ta1_y_yyyyyy_xxy_0[i] = 5.0 * ta1_y_yyyy_xxy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxy_1[i] * fe_0 + ta1_y_yyyyy_xx_0[i] * fe_0 -
                                ta1_y_yyyyy_xx_1[i] * fe_0 + ta_yyyyy_xxy_1[i] + ta1_y_yyyyy_xxy_0[i] * pa_y[i] - ta1_y_yyyyy_xxy_1[i] * pc_y[i];

        ta1_y_yyyyyy_xxz_0[i] = 5.0 * ta1_y_yyyy_xxz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xxz_1[i] * fe_0 + ta_yyyyy_xxz_1[i] +
                                ta1_y_yyyyy_xxz_0[i] * pa_y[i] - ta1_y_yyyyy_xxz_1[i] * pc_y[i];

        ta1_y_yyyyyy_xyy_0[i] = 5.0 * ta1_y_yyyy_xyy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xyy_1[i] * fe_0 + 2.0 * ta1_y_yyyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_y_yyyyy_xy_1[i] * fe_0 + ta_yyyyy_xyy_1[i] + ta1_y_yyyyy_xyy_0[i] * pa_y[i] -
                                ta1_y_yyyyy_xyy_1[i] * pc_y[i];

        ta1_y_yyyyyy_xyz_0[i] = 5.0 * ta1_y_yyyy_xyz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xyz_1[i] * fe_0 + ta1_y_yyyyy_xz_0[i] * fe_0 -
                                ta1_y_yyyyy_xz_1[i] * fe_0 + ta_yyyyy_xyz_1[i] + ta1_y_yyyyy_xyz_0[i] * pa_y[i] - ta1_y_yyyyy_xyz_1[i] * pc_y[i];

        ta1_y_yyyyyy_xzz_0[i] = 5.0 * ta1_y_yyyy_xzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xzz_1[i] * fe_0 + ta_yyyyy_xzz_1[i] +
                                ta1_y_yyyyy_xzz_0[i] * pa_y[i] - ta1_y_yyyyy_xzz_1[i] * pc_y[i];

        ta1_y_yyyyyy_yyy_0[i] = 5.0 * ta1_y_yyyy_yyy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yyy_1[i] * fe_0 + 3.0 * ta1_y_yyyyy_yy_0[i] * fe_0 -
                                3.0 * ta1_y_yyyyy_yy_1[i] * fe_0 + ta_yyyyy_yyy_1[i] + ta1_y_yyyyy_yyy_0[i] * pa_y[i] -
                                ta1_y_yyyyy_yyy_1[i] * pc_y[i];

        ta1_y_yyyyyy_yyz_0[i] = 5.0 * ta1_y_yyyy_yyz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yyz_1[i] * fe_0 + 2.0 * ta1_y_yyyyy_yz_0[i] * fe_0 -
                                2.0 * ta1_y_yyyyy_yz_1[i] * fe_0 + ta_yyyyy_yyz_1[i] + ta1_y_yyyyy_yyz_0[i] * pa_y[i] -
                                ta1_y_yyyyy_yyz_1[i] * pc_y[i];

        ta1_y_yyyyyy_yzz_0[i] = 5.0 * ta1_y_yyyy_yzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yzz_1[i] * fe_0 + ta1_y_yyyyy_zz_0[i] * fe_0 -
                                ta1_y_yyyyy_zz_1[i] * fe_0 + ta_yyyyy_yzz_1[i] + ta1_y_yyyyy_yzz_0[i] * pa_y[i] - ta1_y_yyyyy_yzz_1[i] * pc_y[i];

        ta1_y_yyyyyy_zzz_0[i] = 5.0 * ta1_y_yyyy_zzz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_zzz_1[i] * fe_0 + ta_yyyyy_zzz_1[i] +
                                ta1_y_yyyyy_zzz_0[i] * pa_y[i] - ta1_y_yyyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 500-510 components of targeted buffer : IF

    auto ta1_y_yyyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 500);

    auto ta1_y_yyyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 501);

    auto ta1_y_yyyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 502);

    auto ta1_y_yyyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 503);

    auto ta1_y_yyyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 504);

    auto ta1_y_yyyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 505);

    auto ta1_y_yyyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 506);

    auto ta1_y_yyyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 507);

    auto ta1_y_yyyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 508);

    auto ta1_y_yyyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 509);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_yyyyy_xx_0,   \
                             ta1_y_yyyyy_xx_1,   \
                             ta1_y_yyyyy_xxx_0,  \
                             ta1_y_yyyyy_xxx_1,  \
                             ta1_y_yyyyy_xxy_0,  \
                             ta1_y_yyyyy_xxy_1,  \
                             ta1_y_yyyyy_xxz_0,  \
                             ta1_y_yyyyy_xxz_1,  \
                             ta1_y_yyyyy_xy_0,   \
                             ta1_y_yyyyy_xy_1,   \
                             ta1_y_yyyyy_xyy_0,  \
                             ta1_y_yyyyy_xyy_1,  \
                             ta1_y_yyyyy_xyz_0,  \
                             ta1_y_yyyyy_xyz_1,  \
                             ta1_y_yyyyy_xz_0,   \
                             ta1_y_yyyyy_xz_1,   \
                             ta1_y_yyyyy_xzz_0,  \
                             ta1_y_yyyyy_xzz_1,  \
                             ta1_y_yyyyy_yy_0,   \
                             ta1_y_yyyyy_yy_1,   \
                             ta1_y_yyyyy_yyy_0,  \
                             ta1_y_yyyyy_yyy_1,  \
                             ta1_y_yyyyy_yyz_0,  \
                             ta1_y_yyyyy_yyz_1,  \
                             ta1_y_yyyyy_yz_0,   \
                             ta1_y_yyyyy_yz_1,   \
                             ta1_y_yyyyy_yzz_0,  \
                             ta1_y_yyyyy_yzz_1,  \
                             ta1_y_yyyyy_zz_0,   \
                             ta1_y_yyyyy_zz_1,   \
                             ta1_y_yyyyy_zzz_0,  \
                             ta1_y_yyyyy_zzz_1,  \
                             ta1_y_yyyyyz_xxx_0, \
                             ta1_y_yyyyyz_xxy_0, \
                             ta1_y_yyyyyz_xxz_0, \
                             ta1_y_yyyyyz_xyy_0, \
                             ta1_y_yyyyyz_xyz_0, \
                             ta1_y_yyyyyz_xzz_0, \
                             ta1_y_yyyyyz_yyy_0, \
                             ta1_y_yyyyyz_yyz_0, \
                             ta1_y_yyyyyz_yzz_0, \
                             ta1_y_yyyyyz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyz_xxx_0[i] = ta1_y_yyyyy_xxx_0[i] * pa_z[i] - ta1_y_yyyyy_xxx_1[i] * pc_z[i];

        ta1_y_yyyyyz_xxy_0[i] = ta1_y_yyyyy_xxy_0[i] * pa_z[i] - ta1_y_yyyyy_xxy_1[i] * pc_z[i];

        ta1_y_yyyyyz_xxz_0[i] =
            ta1_y_yyyyy_xx_0[i] * fe_0 - ta1_y_yyyyy_xx_1[i] * fe_0 + ta1_y_yyyyy_xxz_0[i] * pa_z[i] - ta1_y_yyyyy_xxz_1[i] * pc_z[i];

        ta1_y_yyyyyz_xyy_0[i] = ta1_y_yyyyy_xyy_0[i] * pa_z[i] - ta1_y_yyyyy_xyy_1[i] * pc_z[i];

        ta1_y_yyyyyz_xyz_0[i] =
            ta1_y_yyyyy_xy_0[i] * fe_0 - ta1_y_yyyyy_xy_1[i] * fe_0 + ta1_y_yyyyy_xyz_0[i] * pa_z[i] - ta1_y_yyyyy_xyz_1[i] * pc_z[i];

        ta1_y_yyyyyz_xzz_0[i] =
            2.0 * ta1_y_yyyyy_xz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_xz_1[i] * fe_0 + ta1_y_yyyyy_xzz_0[i] * pa_z[i] - ta1_y_yyyyy_xzz_1[i] * pc_z[i];

        ta1_y_yyyyyz_yyy_0[i] = ta1_y_yyyyy_yyy_0[i] * pa_z[i] - ta1_y_yyyyy_yyy_1[i] * pc_z[i];

        ta1_y_yyyyyz_yyz_0[i] =
            ta1_y_yyyyy_yy_0[i] * fe_0 - ta1_y_yyyyy_yy_1[i] * fe_0 + ta1_y_yyyyy_yyz_0[i] * pa_z[i] - ta1_y_yyyyy_yyz_1[i] * pc_z[i];

        ta1_y_yyyyyz_yzz_0[i] =
            2.0 * ta1_y_yyyyy_yz_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_yz_1[i] * fe_0 + ta1_y_yyyyy_yzz_0[i] * pa_z[i] - ta1_y_yyyyy_yzz_1[i] * pc_z[i];

        ta1_y_yyyyyz_zzz_0[i] =
            3.0 * ta1_y_yyyyy_zz_0[i] * fe_0 - 3.0 * ta1_y_yyyyy_zz_1[i] * fe_0 + ta1_y_yyyyy_zzz_0[i] * pa_z[i] - ta1_y_yyyyy_zzz_1[i] * pc_z[i];
    }

    // Set up 510-520 components of targeted buffer : IF

    auto ta1_y_yyyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 510);

    auto ta1_y_yyyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 511);

    auto ta1_y_yyyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 512);

    auto ta1_y_yyyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 513);

    auto ta1_y_yyyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 514);

    auto ta1_y_yyyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 515);

    auto ta1_y_yyyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 516);

    auto ta1_y_yyyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 517);

    auto ta1_y_yyyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 518);

    auto ta1_y_yyyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 519);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yyyy_xxx_0,   \
                             ta1_y_yyyy_xxx_1,   \
                             ta1_y_yyyy_xxy_0,   \
                             ta1_y_yyyy_xxy_1,   \
                             ta1_y_yyyy_xyy_0,   \
                             ta1_y_yyyy_xyy_1,   \
                             ta1_y_yyyy_xyz_0,   \
                             ta1_y_yyyy_xyz_1,   \
                             ta1_y_yyyy_yyy_0,   \
                             ta1_y_yyyy_yyy_1,   \
                             ta1_y_yyyy_yyz_0,   \
                             ta1_y_yyyy_yyz_1,   \
                             ta1_y_yyyy_yzz_0,   \
                             ta1_y_yyyy_yzz_1,   \
                             ta1_y_yyyyz_xxx_0,  \
                             ta1_y_yyyyz_xxx_1,  \
                             ta1_y_yyyyz_xxy_0,  \
                             ta1_y_yyyyz_xxy_1,  \
                             ta1_y_yyyyz_xy_0,   \
                             ta1_y_yyyyz_xy_1,   \
                             ta1_y_yyyyz_xyy_0,  \
                             ta1_y_yyyyz_xyy_1,  \
                             ta1_y_yyyyz_xyz_0,  \
                             ta1_y_yyyyz_xyz_1,  \
                             ta1_y_yyyyz_yy_0,   \
                             ta1_y_yyyyz_yy_1,   \
                             ta1_y_yyyyz_yyy_0,  \
                             ta1_y_yyyyz_yyy_1,  \
                             ta1_y_yyyyz_yyz_0,  \
                             ta1_y_yyyyz_yyz_1,  \
                             ta1_y_yyyyz_yz_0,   \
                             ta1_y_yyyyz_yz_1,   \
                             ta1_y_yyyyz_yzz_0,  \
                             ta1_y_yyyyz_yzz_1,  \
                             ta1_y_yyyyzz_xxx_0, \
                             ta1_y_yyyyzz_xxy_0, \
                             ta1_y_yyyyzz_xxz_0, \
                             ta1_y_yyyyzz_xyy_0, \
                             ta1_y_yyyyzz_xyz_0, \
                             ta1_y_yyyyzz_xzz_0, \
                             ta1_y_yyyyzz_yyy_0, \
                             ta1_y_yyyyzz_yyz_0, \
                             ta1_y_yyyyzz_yzz_0, \
                             ta1_y_yyyyzz_zzz_0, \
                             ta1_y_yyyzz_xxz_0,  \
                             ta1_y_yyyzz_xxz_1,  \
                             ta1_y_yyyzz_xzz_0,  \
                             ta1_y_yyyzz_xzz_1,  \
                             ta1_y_yyyzz_zzz_0,  \
                             ta1_y_yyyzz_zzz_1,  \
                             ta1_y_yyzz_xxz_0,   \
                             ta1_y_yyzz_xxz_1,   \
                             ta1_y_yyzz_xzz_0,   \
                             ta1_y_yyzz_xzz_1,   \
                             ta1_y_yyzz_zzz_0,   \
                             ta1_y_yyzz_zzz_1,   \
                             ta_yyyzz_xxz_1,     \
                             ta_yyyzz_xzz_1,     \
                             ta_yyyzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyzz_xxx_0[i] =
            ta1_y_yyyy_xxx_0[i] * fe_0 - ta1_y_yyyy_xxx_1[i] * fe_0 + ta1_y_yyyyz_xxx_0[i] * pa_z[i] - ta1_y_yyyyz_xxx_1[i] * pc_z[i];

        ta1_y_yyyyzz_xxy_0[i] =
            ta1_y_yyyy_xxy_0[i] * fe_0 - ta1_y_yyyy_xxy_1[i] * fe_0 + ta1_y_yyyyz_xxy_0[i] * pa_z[i] - ta1_y_yyyyz_xxy_1[i] * pc_z[i];

        ta1_y_yyyyzz_xxz_0[i] = 3.0 * ta1_y_yyzz_xxz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxz_1[i] * fe_0 + ta_yyyzz_xxz_1[i] +
                                ta1_y_yyyzz_xxz_0[i] * pa_y[i] - ta1_y_yyyzz_xxz_1[i] * pc_y[i];

        ta1_y_yyyyzz_xyy_0[i] =
            ta1_y_yyyy_xyy_0[i] * fe_0 - ta1_y_yyyy_xyy_1[i] * fe_0 + ta1_y_yyyyz_xyy_0[i] * pa_z[i] - ta1_y_yyyyz_xyy_1[i] * pc_z[i];

        ta1_y_yyyyzz_xyz_0[i] = ta1_y_yyyy_xyz_0[i] * fe_0 - ta1_y_yyyy_xyz_1[i] * fe_0 + ta1_y_yyyyz_xy_0[i] * fe_0 - ta1_y_yyyyz_xy_1[i] * fe_0 +
                                ta1_y_yyyyz_xyz_0[i] * pa_z[i] - ta1_y_yyyyz_xyz_1[i] * pc_z[i];

        ta1_y_yyyyzz_xzz_0[i] = 3.0 * ta1_y_yyzz_xzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xzz_1[i] * fe_0 + ta_yyyzz_xzz_1[i] +
                                ta1_y_yyyzz_xzz_0[i] * pa_y[i] - ta1_y_yyyzz_xzz_1[i] * pc_y[i];

        ta1_y_yyyyzz_yyy_0[i] =
            ta1_y_yyyy_yyy_0[i] * fe_0 - ta1_y_yyyy_yyy_1[i] * fe_0 + ta1_y_yyyyz_yyy_0[i] * pa_z[i] - ta1_y_yyyyz_yyy_1[i] * pc_z[i];

        ta1_y_yyyyzz_yyz_0[i] = ta1_y_yyyy_yyz_0[i] * fe_0 - ta1_y_yyyy_yyz_1[i] * fe_0 + ta1_y_yyyyz_yy_0[i] * fe_0 - ta1_y_yyyyz_yy_1[i] * fe_0 +
                                ta1_y_yyyyz_yyz_0[i] * pa_z[i] - ta1_y_yyyyz_yyz_1[i] * pc_z[i];

        ta1_y_yyyyzz_yzz_0[i] = ta1_y_yyyy_yzz_0[i] * fe_0 - ta1_y_yyyy_yzz_1[i] * fe_0 + 2.0 * ta1_y_yyyyz_yz_0[i] * fe_0 -
                                2.0 * ta1_y_yyyyz_yz_1[i] * fe_0 + ta1_y_yyyyz_yzz_0[i] * pa_z[i] - ta1_y_yyyyz_yzz_1[i] * pc_z[i];

        ta1_y_yyyyzz_zzz_0[i] = 3.0 * ta1_y_yyzz_zzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_zzz_1[i] * fe_0 + ta_yyyzz_zzz_1[i] +
                                ta1_y_yyyzz_zzz_0[i] * pa_y[i] - ta1_y_yyyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 520-530 components of targeted buffer : IF

    auto ta1_y_yyyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 520);

    auto ta1_y_yyyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 521);

    auto ta1_y_yyyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 522);

    auto ta1_y_yyyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 523);

    auto ta1_y_yyyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 524);

    auto ta1_y_yyyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 525);

    auto ta1_y_yyyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 526);

    auto ta1_y_yyyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 527);

    auto ta1_y_yyyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 528);

    auto ta1_y_yyyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 529);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yyyz_xxx_0,   \
                             ta1_y_yyyz_xxx_1,   \
                             ta1_y_yyyz_xxy_0,   \
                             ta1_y_yyyz_xxy_1,   \
                             ta1_y_yyyz_xyy_0,   \
                             ta1_y_yyyz_xyy_1,   \
                             ta1_y_yyyz_xyz_0,   \
                             ta1_y_yyyz_xyz_1,   \
                             ta1_y_yyyz_yyy_0,   \
                             ta1_y_yyyz_yyy_1,   \
                             ta1_y_yyyz_yyz_0,   \
                             ta1_y_yyyz_yyz_1,   \
                             ta1_y_yyyz_yzz_0,   \
                             ta1_y_yyyz_yzz_1,   \
                             ta1_y_yyyzz_xxx_0,  \
                             ta1_y_yyyzz_xxx_1,  \
                             ta1_y_yyyzz_xxy_0,  \
                             ta1_y_yyyzz_xxy_1,  \
                             ta1_y_yyyzz_xy_0,   \
                             ta1_y_yyyzz_xy_1,   \
                             ta1_y_yyyzz_xyy_0,  \
                             ta1_y_yyyzz_xyy_1,  \
                             ta1_y_yyyzz_xyz_0,  \
                             ta1_y_yyyzz_xyz_1,  \
                             ta1_y_yyyzz_yy_0,   \
                             ta1_y_yyyzz_yy_1,   \
                             ta1_y_yyyzz_yyy_0,  \
                             ta1_y_yyyzz_yyy_1,  \
                             ta1_y_yyyzz_yyz_0,  \
                             ta1_y_yyyzz_yyz_1,  \
                             ta1_y_yyyzz_yz_0,   \
                             ta1_y_yyyzz_yz_1,   \
                             ta1_y_yyyzz_yzz_0,  \
                             ta1_y_yyyzz_yzz_1,  \
                             ta1_y_yyyzzz_xxx_0, \
                             ta1_y_yyyzzz_xxy_0, \
                             ta1_y_yyyzzz_xxz_0, \
                             ta1_y_yyyzzz_xyy_0, \
                             ta1_y_yyyzzz_xyz_0, \
                             ta1_y_yyyzzz_xzz_0, \
                             ta1_y_yyyzzz_yyy_0, \
                             ta1_y_yyyzzz_yyz_0, \
                             ta1_y_yyyzzz_yzz_0, \
                             ta1_y_yyyzzz_zzz_0, \
                             ta1_y_yyzzz_xxz_0,  \
                             ta1_y_yyzzz_xxz_1,  \
                             ta1_y_yyzzz_xzz_0,  \
                             ta1_y_yyzzz_xzz_1,  \
                             ta1_y_yyzzz_zzz_0,  \
                             ta1_y_yyzzz_zzz_1,  \
                             ta1_y_yzzz_xxz_0,   \
                             ta1_y_yzzz_xxz_1,   \
                             ta1_y_yzzz_xzz_0,   \
                             ta1_y_yzzz_xzz_1,   \
                             ta1_y_yzzz_zzz_0,   \
                             ta1_y_yzzz_zzz_1,   \
                             ta_yyzzz_xxz_1,     \
                             ta_yyzzz_xzz_1,     \
                             ta_yyzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzzz_xxx_0[i] =
            2.0 * ta1_y_yyyz_xxx_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xxx_1[i] * fe_0 + ta1_y_yyyzz_xxx_0[i] * pa_z[i] - ta1_y_yyyzz_xxx_1[i] * pc_z[i];

        ta1_y_yyyzzz_xxy_0[i] =
            2.0 * ta1_y_yyyz_xxy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xxy_1[i] * fe_0 + ta1_y_yyyzz_xxy_0[i] * pa_z[i] - ta1_y_yyyzz_xxy_1[i] * pc_z[i];

        ta1_y_yyyzzz_xxz_0[i] = 2.0 * ta1_y_yzzz_xxz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xxz_1[i] * fe_0 + ta_yyzzz_xxz_1[i] +
                                ta1_y_yyzzz_xxz_0[i] * pa_y[i] - ta1_y_yyzzz_xxz_1[i] * pc_y[i];

        ta1_y_yyyzzz_xyy_0[i] =
            2.0 * ta1_y_yyyz_xyy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xyy_1[i] * fe_0 + ta1_y_yyyzz_xyy_0[i] * pa_z[i] - ta1_y_yyyzz_xyy_1[i] * pc_z[i];

        ta1_y_yyyzzz_xyz_0[i] = 2.0 * ta1_y_yyyz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xyz_1[i] * fe_0 + ta1_y_yyyzz_xy_0[i] * fe_0 -
                                ta1_y_yyyzz_xy_1[i] * fe_0 + ta1_y_yyyzz_xyz_0[i] * pa_z[i] - ta1_y_yyyzz_xyz_1[i] * pc_z[i];

        ta1_y_yyyzzz_xzz_0[i] = 2.0 * ta1_y_yzzz_xzz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xzz_1[i] * fe_0 + ta_yyzzz_xzz_1[i] +
                                ta1_y_yyzzz_xzz_0[i] * pa_y[i] - ta1_y_yyzzz_xzz_1[i] * pc_y[i];

        ta1_y_yyyzzz_yyy_0[i] =
            2.0 * ta1_y_yyyz_yyy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yyy_1[i] * fe_0 + ta1_y_yyyzz_yyy_0[i] * pa_z[i] - ta1_y_yyyzz_yyy_1[i] * pc_z[i];

        ta1_y_yyyzzz_yyz_0[i] = 2.0 * ta1_y_yyyz_yyz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yyz_1[i] * fe_0 + ta1_y_yyyzz_yy_0[i] * fe_0 -
                                ta1_y_yyyzz_yy_1[i] * fe_0 + ta1_y_yyyzz_yyz_0[i] * pa_z[i] - ta1_y_yyyzz_yyz_1[i] * pc_z[i];

        ta1_y_yyyzzz_yzz_0[i] = 2.0 * ta1_y_yyyz_yzz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yzz_1[i] * fe_0 + 2.0 * ta1_y_yyyzz_yz_0[i] * fe_0 -
                                2.0 * ta1_y_yyyzz_yz_1[i] * fe_0 + ta1_y_yyyzz_yzz_0[i] * pa_z[i] - ta1_y_yyyzz_yzz_1[i] * pc_z[i];

        ta1_y_yyyzzz_zzz_0[i] = 2.0 * ta1_y_yzzz_zzz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_zzz_1[i] * fe_0 + ta_yyzzz_zzz_1[i] +
                                ta1_y_yyzzz_zzz_0[i] * pa_y[i] - ta1_y_yyzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 530-540 components of targeted buffer : IF

    auto ta1_y_yyzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 530);

    auto ta1_y_yyzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 531);

    auto ta1_y_yyzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 532);

    auto ta1_y_yyzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 533);

    auto ta1_y_yyzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 534);

    auto ta1_y_yyzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 535);

    auto ta1_y_yyzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 536);

    auto ta1_y_yyzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 537);

    auto ta1_y_yyzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 538);

    auto ta1_y_yyzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 539);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yyzz_xxx_0,   \
                             ta1_y_yyzz_xxx_1,   \
                             ta1_y_yyzz_xxy_0,   \
                             ta1_y_yyzz_xxy_1,   \
                             ta1_y_yyzz_xyy_0,   \
                             ta1_y_yyzz_xyy_1,   \
                             ta1_y_yyzz_xyz_0,   \
                             ta1_y_yyzz_xyz_1,   \
                             ta1_y_yyzz_yyy_0,   \
                             ta1_y_yyzz_yyy_1,   \
                             ta1_y_yyzz_yyz_0,   \
                             ta1_y_yyzz_yyz_1,   \
                             ta1_y_yyzz_yzz_0,   \
                             ta1_y_yyzz_yzz_1,   \
                             ta1_y_yyzzz_xxx_0,  \
                             ta1_y_yyzzz_xxx_1,  \
                             ta1_y_yyzzz_xxy_0,  \
                             ta1_y_yyzzz_xxy_1,  \
                             ta1_y_yyzzz_xy_0,   \
                             ta1_y_yyzzz_xy_1,   \
                             ta1_y_yyzzz_xyy_0,  \
                             ta1_y_yyzzz_xyy_1,  \
                             ta1_y_yyzzz_xyz_0,  \
                             ta1_y_yyzzz_xyz_1,  \
                             ta1_y_yyzzz_yy_0,   \
                             ta1_y_yyzzz_yy_1,   \
                             ta1_y_yyzzz_yyy_0,  \
                             ta1_y_yyzzz_yyy_1,  \
                             ta1_y_yyzzz_yyz_0,  \
                             ta1_y_yyzzz_yyz_1,  \
                             ta1_y_yyzzz_yz_0,   \
                             ta1_y_yyzzz_yz_1,   \
                             ta1_y_yyzzz_yzz_0,  \
                             ta1_y_yyzzz_yzz_1,  \
                             ta1_y_yyzzzz_xxx_0, \
                             ta1_y_yyzzzz_xxy_0, \
                             ta1_y_yyzzzz_xxz_0, \
                             ta1_y_yyzzzz_xyy_0, \
                             ta1_y_yyzzzz_xyz_0, \
                             ta1_y_yyzzzz_xzz_0, \
                             ta1_y_yyzzzz_yyy_0, \
                             ta1_y_yyzzzz_yyz_0, \
                             ta1_y_yyzzzz_yzz_0, \
                             ta1_y_yyzzzz_zzz_0, \
                             ta1_y_yzzzz_xxz_0,  \
                             ta1_y_yzzzz_xxz_1,  \
                             ta1_y_yzzzz_xzz_0,  \
                             ta1_y_yzzzz_xzz_1,  \
                             ta1_y_yzzzz_zzz_0,  \
                             ta1_y_yzzzz_zzz_1,  \
                             ta1_y_zzzz_xxz_0,   \
                             ta1_y_zzzz_xxz_1,   \
                             ta1_y_zzzz_xzz_0,   \
                             ta1_y_zzzz_xzz_1,   \
                             ta1_y_zzzz_zzz_0,   \
                             ta1_y_zzzz_zzz_1,   \
                             ta_yzzzz_xxz_1,     \
                             ta_yzzzz_xzz_1,     \
                             ta_yzzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzzz_xxx_0[i] =
            3.0 * ta1_y_yyzz_xxx_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxx_1[i] * fe_0 + ta1_y_yyzzz_xxx_0[i] * pa_z[i] - ta1_y_yyzzz_xxx_1[i] * pc_z[i];

        ta1_y_yyzzzz_xxy_0[i] =
            3.0 * ta1_y_yyzz_xxy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xxy_1[i] * fe_0 + ta1_y_yyzzz_xxy_0[i] * pa_z[i] - ta1_y_yyzzz_xxy_1[i] * pc_z[i];

        ta1_y_yyzzzz_xxz_0[i] = ta1_y_zzzz_xxz_0[i] * fe_0 - ta1_y_zzzz_xxz_1[i] * fe_0 + ta_yzzzz_xxz_1[i] + ta1_y_yzzzz_xxz_0[i] * pa_y[i] -
                                ta1_y_yzzzz_xxz_1[i] * pc_y[i];

        ta1_y_yyzzzz_xyy_0[i] =
            3.0 * ta1_y_yyzz_xyy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xyy_1[i] * fe_0 + ta1_y_yyzzz_xyy_0[i] * pa_z[i] - ta1_y_yyzzz_xyy_1[i] * pc_z[i];

        ta1_y_yyzzzz_xyz_0[i] = 3.0 * ta1_y_yyzz_xyz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xyz_1[i] * fe_0 + ta1_y_yyzzz_xy_0[i] * fe_0 -
                                ta1_y_yyzzz_xy_1[i] * fe_0 + ta1_y_yyzzz_xyz_0[i] * pa_z[i] - ta1_y_yyzzz_xyz_1[i] * pc_z[i];

        ta1_y_yyzzzz_xzz_0[i] = ta1_y_zzzz_xzz_0[i] * fe_0 - ta1_y_zzzz_xzz_1[i] * fe_0 + ta_yzzzz_xzz_1[i] + ta1_y_yzzzz_xzz_0[i] * pa_y[i] -
                                ta1_y_yzzzz_xzz_1[i] * pc_y[i];

        ta1_y_yyzzzz_yyy_0[i] =
            3.0 * ta1_y_yyzz_yyy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yyy_1[i] * fe_0 + ta1_y_yyzzz_yyy_0[i] * pa_z[i] - ta1_y_yyzzz_yyy_1[i] * pc_z[i];

        ta1_y_yyzzzz_yyz_0[i] = 3.0 * ta1_y_yyzz_yyz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yyz_1[i] * fe_0 + ta1_y_yyzzz_yy_0[i] * fe_0 -
                                ta1_y_yyzzz_yy_1[i] * fe_0 + ta1_y_yyzzz_yyz_0[i] * pa_z[i] - ta1_y_yyzzz_yyz_1[i] * pc_z[i];

        ta1_y_yyzzzz_yzz_0[i] = 3.0 * ta1_y_yyzz_yzz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yzz_1[i] * fe_0 + 2.0 * ta1_y_yyzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_y_yyzzz_yz_1[i] * fe_0 + ta1_y_yyzzz_yzz_0[i] * pa_z[i] - ta1_y_yyzzz_yzz_1[i] * pc_z[i];

        ta1_y_yyzzzz_zzz_0[i] = ta1_y_zzzz_zzz_0[i] * fe_0 - ta1_y_zzzz_zzz_1[i] * fe_0 + ta_yzzzz_zzz_1[i] + ta1_y_yzzzz_zzz_0[i] * pa_y[i] -
                                ta1_y_yzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 540-550 components of targeted buffer : IF

    auto ta1_y_yzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 540);

    auto ta1_y_yzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 541);

    auto ta1_y_yzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 542);

    auto ta1_y_yzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 543);

    auto ta1_y_yzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 544);

    auto ta1_y_yzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 545);

    auto ta1_y_yzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 546);

    auto ta1_y_yzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 547);

    auto ta1_y_yzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 548);

    auto ta1_y_yzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 549);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yzzz_xxy_0,   \
                             ta1_y_yzzz_xxy_1,   \
                             ta1_y_yzzz_xyy_0,   \
                             ta1_y_yzzz_xyy_1,   \
                             ta1_y_yzzz_yyy_0,   \
                             ta1_y_yzzz_yyy_1,   \
                             ta1_y_yzzzz_xxy_0,  \
                             ta1_y_yzzzz_xxy_1,  \
                             ta1_y_yzzzz_xyy_0,  \
                             ta1_y_yzzzz_xyy_1,  \
                             ta1_y_yzzzz_yyy_0,  \
                             ta1_y_yzzzz_yyy_1,  \
                             ta1_y_yzzzzz_xxx_0, \
                             ta1_y_yzzzzz_xxy_0, \
                             ta1_y_yzzzzz_xxz_0, \
                             ta1_y_yzzzzz_xyy_0, \
                             ta1_y_yzzzzz_xyz_0, \
                             ta1_y_yzzzzz_xzz_0, \
                             ta1_y_yzzzzz_yyy_0, \
                             ta1_y_yzzzzz_yyz_0, \
                             ta1_y_yzzzzz_yzz_0, \
                             ta1_y_yzzzzz_zzz_0, \
                             ta1_y_zzzzz_xxx_0,  \
                             ta1_y_zzzzz_xxx_1,  \
                             ta1_y_zzzzz_xxz_0,  \
                             ta1_y_zzzzz_xxz_1,  \
                             ta1_y_zzzzz_xyz_0,  \
                             ta1_y_zzzzz_xyz_1,  \
                             ta1_y_zzzzz_xz_0,   \
                             ta1_y_zzzzz_xz_1,   \
                             ta1_y_zzzzz_xzz_0,  \
                             ta1_y_zzzzz_xzz_1,  \
                             ta1_y_zzzzz_yyz_0,  \
                             ta1_y_zzzzz_yyz_1,  \
                             ta1_y_zzzzz_yz_0,   \
                             ta1_y_zzzzz_yz_1,   \
                             ta1_y_zzzzz_yzz_0,  \
                             ta1_y_zzzzz_yzz_1,  \
                             ta1_y_zzzzz_zz_0,   \
                             ta1_y_zzzzz_zz_1,   \
                             ta1_y_zzzzz_zzz_0,  \
                             ta1_y_zzzzz_zzz_1,  \
                             ta_zzzzz_xxx_1,     \
                             ta_zzzzz_xxz_1,     \
                             ta_zzzzz_xyz_1,     \
                             ta_zzzzz_xzz_1,     \
                             ta_zzzzz_yyz_1,     \
                             ta_zzzzz_yzz_1,     \
                             ta_zzzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzzz_xxx_0[i] = ta_zzzzz_xxx_1[i] + ta1_y_zzzzz_xxx_0[i] * pa_y[i] - ta1_y_zzzzz_xxx_1[i] * pc_y[i];

        ta1_y_yzzzzz_xxy_0[i] =
            4.0 * ta1_y_yzzz_xxy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_xxy_1[i] * fe_0 + ta1_y_yzzzz_xxy_0[i] * pa_z[i] - ta1_y_yzzzz_xxy_1[i] * pc_z[i];

        ta1_y_yzzzzz_xxz_0[i] = ta_zzzzz_xxz_1[i] + ta1_y_zzzzz_xxz_0[i] * pa_y[i] - ta1_y_zzzzz_xxz_1[i] * pc_y[i];

        ta1_y_yzzzzz_xyy_0[i] =
            4.0 * ta1_y_yzzz_xyy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_xyy_1[i] * fe_0 + ta1_y_yzzzz_xyy_0[i] * pa_z[i] - ta1_y_yzzzz_xyy_1[i] * pc_z[i];

        ta1_y_yzzzzz_xyz_0[i] = ta1_y_zzzzz_xz_0[i] * fe_0 - ta1_y_zzzzz_xz_1[i] * fe_0 + ta_zzzzz_xyz_1[i] + ta1_y_zzzzz_xyz_0[i] * pa_y[i] -
                                ta1_y_zzzzz_xyz_1[i] * pc_y[i];

        ta1_y_yzzzzz_xzz_0[i] = ta_zzzzz_xzz_1[i] + ta1_y_zzzzz_xzz_0[i] * pa_y[i] - ta1_y_zzzzz_xzz_1[i] * pc_y[i];

        ta1_y_yzzzzz_yyy_0[i] =
            4.0 * ta1_y_yzzz_yyy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_yyy_1[i] * fe_0 + ta1_y_yzzzz_yyy_0[i] * pa_z[i] - ta1_y_yzzzz_yyy_1[i] * pc_z[i];

        ta1_y_yzzzzz_yyz_0[i] = 2.0 * ta1_y_zzzzz_yz_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_yz_1[i] * fe_0 + ta_zzzzz_yyz_1[i] +
                                ta1_y_zzzzz_yyz_0[i] * pa_y[i] - ta1_y_zzzzz_yyz_1[i] * pc_y[i];

        ta1_y_yzzzzz_yzz_0[i] = ta1_y_zzzzz_zz_0[i] * fe_0 - ta1_y_zzzzz_zz_1[i] * fe_0 + ta_zzzzz_yzz_1[i] + ta1_y_zzzzz_yzz_0[i] * pa_y[i] -
                                ta1_y_zzzzz_yzz_1[i] * pc_y[i];

        ta1_y_yzzzzz_zzz_0[i] = ta_zzzzz_zzz_1[i] + ta1_y_zzzzz_zzz_0[i] * pa_y[i] - ta1_y_zzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 550-560 components of targeted buffer : IF

    auto ta1_y_zzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 550);

    auto ta1_y_zzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 551);

    auto ta1_y_zzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 552);

    auto ta1_y_zzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 553);

    auto ta1_y_zzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 554);

    auto ta1_y_zzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 555);

    auto ta1_y_zzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 556);

    auto ta1_y_zzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 557);

    auto ta1_y_zzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 558);

    auto ta1_y_zzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 559);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_zzzz_xxx_0,   \
                             ta1_y_zzzz_xxx_1,   \
                             ta1_y_zzzz_xxy_0,   \
                             ta1_y_zzzz_xxy_1,   \
                             ta1_y_zzzz_xxz_0,   \
                             ta1_y_zzzz_xxz_1,   \
                             ta1_y_zzzz_xyy_0,   \
                             ta1_y_zzzz_xyy_1,   \
                             ta1_y_zzzz_xyz_0,   \
                             ta1_y_zzzz_xyz_1,   \
                             ta1_y_zzzz_xzz_0,   \
                             ta1_y_zzzz_xzz_1,   \
                             ta1_y_zzzz_yyy_0,   \
                             ta1_y_zzzz_yyy_1,   \
                             ta1_y_zzzz_yyz_0,   \
                             ta1_y_zzzz_yyz_1,   \
                             ta1_y_zzzz_yzz_0,   \
                             ta1_y_zzzz_yzz_1,   \
                             ta1_y_zzzz_zzz_0,   \
                             ta1_y_zzzz_zzz_1,   \
                             ta1_y_zzzzz_xx_0,   \
                             ta1_y_zzzzz_xx_1,   \
                             ta1_y_zzzzz_xxx_0,  \
                             ta1_y_zzzzz_xxx_1,  \
                             ta1_y_zzzzz_xxy_0,  \
                             ta1_y_zzzzz_xxy_1,  \
                             ta1_y_zzzzz_xxz_0,  \
                             ta1_y_zzzzz_xxz_1,  \
                             ta1_y_zzzzz_xy_0,   \
                             ta1_y_zzzzz_xy_1,   \
                             ta1_y_zzzzz_xyy_0,  \
                             ta1_y_zzzzz_xyy_1,  \
                             ta1_y_zzzzz_xyz_0,  \
                             ta1_y_zzzzz_xyz_1,  \
                             ta1_y_zzzzz_xz_0,   \
                             ta1_y_zzzzz_xz_1,   \
                             ta1_y_zzzzz_xzz_0,  \
                             ta1_y_zzzzz_xzz_1,  \
                             ta1_y_zzzzz_yy_0,   \
                             ta1_y_zzzzz_yy_1,   \
                             ta1_y_zzzzz_yyy_0,  \
                             ta1_y_zzzzz_yyy_1,  \
                             ta1_y_zzzzz_yyz_0,  \
                             ta1_y_zzzzz_yyz_1,  \
                             ta1_y_zzzzz_yz_0,   \
                             ta1_y_zzzzz_yz_1,   \
                             ta1_y_zzzzz_yzz_0,  \
                             ta1_y_zzzzz_yzz_1,  \
                             ta1_y_zzzzz_zz_0,   \
                             ta1_y_zzzzz_zz_1,   \
                             ta1_y_zzzzz_zzz_0,  \
                             ta1_y_zzzzz_zzz_1,  \
                             ta1_y_zzzzzz_xxx_0, \
                             ta1_y_zzzzzz_xxy_0, \
                             ta1_y_zzzzzz_xxz_0, \
                             ta1_y_zzzzzz_xyy_0, \
                             ta1_y_zzzzzz_xyz_0, \
                             ta1_y_zzzzzz_xzz_0, \
                             ta1_y_zzzzzz_yyy_0, \
                             ta1_y_zzzzzz_yyz_0, \
                             ta1_y_zzzzzz_yzz_0, \
                             ta1_y_zzzzzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzzz_xxx_0[i] =
            5.0 * ta1_y_zzzz_xxx_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxx_1[i] * fe_0 + ta1_y_zzzzz_xxx_0[i] * pa_z[i] - ta1_y_zzzzz_xxx_1[i] * pc_z[i];

        ta1_y_zzzzzz_xxy_0[i] =
            5.0 * ta1_y_zzzz_xxy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxy_1[i] * fe_0 + ta1_y_zzzzz_xxy_0[i] * pa_z[i] - ta1_y_zzzzz_xxy_1[i] * pc_z[i];

        ta1_y_zzzzzz_xxz_0[i] = 5.0 * ta1_y_zzzz_xxz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xxz_1[i] * fe_0 + ta1_y_zzzzz_xx_0[i] * fe_0 -
                                ta1_y_zzzzz_xx_1[i] * fe_0 + ta1_y_zzzzz_xxz_0[i] * pa_z[i] - ta1_y_zzzzz_xxz_1[i] * pc_z[i];

        ta1_y_zzzzzz_xyy_0[i] =
            5.0 * ta1_y_zzzz_xyy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xyy_1[i] * fe_0 + ta1_y_zzzzz_xyy_0[i] * pa_z[i] - ta1_y_zzzzz_xyy_1[i] * pc_z[i];

        ta1_y_zzzzzz_xyz_0[i] = 5.0 * ta1_y_zzzz_xyz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xyz_1[i] * fe_0 + ta1_y_zzzzz_xy_0[i] * fe_0 -
                                ta1_y_zzzzz_xy_1[i] * fe_0 + ta1_y_zzzzz_xyz_0[i] * pa_z[i] - ta1_y_zzzzz_xyz_1[i] * pc_z[i];

        ta1_y_zzzzzz_xzz_0[i] = 5.0 * ta1_y_zzzz_xzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xzz_1[i] * fe_0 + 2.0 * ta1_y_zzzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_y_zzzzz_xz_1[i] * fe_0 + ta1_y_zzzzz_xzz_0[i] * pa_z[i] - ta1_y_zzzzz_xzz_1[i] * pc_z[i];

        ta1_y_zzzzzz_yyy_0[i] =
            5.0 * ta1_y_zzzz_yyy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yyy_1[i] * fe_0 + ta1_y_zzzzz_yyy_0[i] * pa_z[i] - ta1_y_zzzzz_yyy_1[i] * pc_z[i];

        ta1_y_zzzzzz_yyz_0[i] = 5.0 * ta1_y_zzzz_yyz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yyz_1[i] * fe_0 + ta1_y_zzzzz_yy_0[i] * fe_0 -
                                ta1_y_zzzzz_yy_1[i] * fe_0 + ta1_y_zzzzz_yyz_0[i] * pa_z[i] - ta1_y_zzzzz_yyz_1[i] * pc_z[i];

        ta1_y_zzzzzz_yzz_0[i] = 5.0 * ta1_y_zzzz_yzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yzz_1[i] * fe_0 + 2.0 * ta1_y_zzzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_y_zzzzz_yz_1[i] * fe_0 + ta1_y_zzzzz_yzz_0[i] * pa_z[i] - ta1_y_zzzzz_yzz_1[i] * pc_z[i];

        ta1_y_zzzzzz_zzz_0[i] = 5.0 * ta1_y_zzzz_zzz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_zzz_1[i] * fe_0 + 3.0 * ta1_y_zzzzz_zz_0[i] * fe_0 -
                                3.0 * ta1_y_zzzzz_zz_1[i] * fe_0 + ta1_y_zzzzz_zzz_0[i] * pa_z[i] - ta1_y_zzzzz_zzz_1[i] * pc_z[i];
    }

    // Set up 560-570 components of targeted buffer : IF

    auto ta1_z_xxxxxx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 560);

    auto ta1_z_xxxxxx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 561);

    auto ta1_z_xxxxxx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 562);

    auto ta1_z_xxxxxx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 563);

    auto ta1_z_xxxxxx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 564);

    auto ta1_z_xxxxxx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 565);

    auto ta1_z_xxxxxx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 566);

    auto ta1_z_xxxxxx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 567);

    auto ta1_z_xxxxxx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 568);

    auto ta1_z_xxxxxx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 569);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xxxx_xxx_0,   \
                             ta1_z_xxxx_xxx_1,   \
                             ta1_z_xxxx_xxy_0,   \
                             ta1_z_xxxx_xxy_1,   \
                             ta1_z_xxxx_xxz_0,   \
                             ta1_z_xxxx_xxz_1,   \
                             ta1_z_xxxx_xyy_0,   \
                             ta1_z_xxxx_xyy_1,   \
                             ta1_z_xxxx_xyz_0,   \
                             ta1_z_xxxx_xyz_1,   \
                             ta1_z_xxxx_xzz_0,   \
                             ta1_z_xxxx_xzz_1,   \
                             ta1_z_xxxx_yyy_0,   \
                             ta1_z_xxxx_yyy_1,   \
                             ta1_z_xxxx_yyz_0,   \
                             ta1_z_xxxx_yyz_1,   \
                             ta1_z_xxxx_yzz_0,   \
                             ta1_z_xxxx_yzz_1,   \
                             ta1_z_xxxx_zzz_0,   \
                             ta1_z_xxxx_zzz_1,   \
                             ta1_z_xxxxx_xx_0,   \
                             ta1_z_xxxxx_xx_1,   \
                             ta1_z_xxxxx_xxx_0,  \
                             ta1_z_xxxxx_xxx_1,  \
                             ta1_z_xxxxx_xxy_0,  \
                             ta1_z_xxxxx_xxy_1,  \
                             ta1_z_xxxxx_xxz_0,  \
                             ta1_z_xxxxx_xxz_1,  \
                             ta1_z_xxxxx_xy_0,   \
                             ta1_z_xxxxx_xy_1,   \
                             ta1_z_xxxxx_xyy_0,  \
                             ta1_z_xxxxx_xyy_1,  \
                             ta1_z_xxxxx_xyz_0,  \
                             ta1_z_xxxxx_xyz_1,  \
                             ta1_z_xxxxx_xz_0,   \
                             ta1_z_xxxxx_xz_1,   \
                             ta1_z_xxxxx_xzz_0,  \
                             ta1_z_xxxxx_xzz_1,  \
                             ta1_z_xxxxx_yy_0,   \
                             ta1_z_xxxxx_yy_1,   \
                             ta1_z_xxxxx_yyy_0,  \
                             ta1_z_xxxxx_yyy_1,  \
                             ta1_z_xxxxx_yyz_0,  \
                             ta1_z_xxxxx_yyz_1,  \
                             ta1_z_xxxxx_yz_0,   \
                             ta1_z_xxxxx_yz_1,   \
                             ta1_z_xxxxx_yzz_0,  \
                             ta1_z_xxxxx_yzz_1,  \
                             ta1_z_xxxxx_zz_0,   \
                             ta1_z_xxxxx_zz_1,   \
                             ta1_z_xxxxx_zzz_0,  \
                             ta1_z_xxxxx_zzz_1,  \
                             ta1_z_xxxxxx_xxx_0, \
                             ta1_z_xxxxxx_xxy_0, \
                             ta1_z_xxxxxx_xxz_0, \
                             ta1_z_xxxxxx_xyy_0, \
                             ta1_z_xxxxxx_xyz_0, \
                             ta1_z_xxxxxx_xzz_0, \
                             ta1_z_xxxxxx_yyy_0, \
                             ta1_z_xxxxxx_yyz_0, \
                             ta1_z_xxxxxx_yzz_0, \
                             ta1_z_xxxxxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxx_xxx_0[i] = 5.0 * ta1_z_xxxx_xxx_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxx_1[i] * fe_0 + 3.0 * ta1_z_xxxxx_xx_0[i] * fe_0 -
                                3.0 * ta1_z_xxxxx_xx_1[i] * fe_0 + ta1_z_xxxxx_xxx_0[i] * pa_x[i] - ta1_z_xxxxx_xxx_1[i] * pc_x[i];

        ta1_z_xxxxxx_xxy_0[i] = 5.0 * ta1_z_xxxx_xxy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxxxx_xy_0[i] * fe_0 -
                                2.0 * ta1_z_xxxxx_xy_1[i] * fe_0 + ta1_z_xxxxx_xxy_0[i] * pa_x[i] - ta1_z_xxxxx_xxy_1[i] * pc_x[i];

        ta1_z_xxxxxx_xxz_0[i] = 5.0 * ta1_z_xxxx_xxz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xxz_1[i] * fe_0 + 2.0 * ta1_z_xxxxx_xz_0[i] * fe_0 -
                                2.0 * ta1_z_xxxxx_xz_1[i] * fe_0 + ta1_z_xxxxx_xxz_0[i] * pa_x[i] - ta1_z_xxxxx_xxz_1[i] * pc_x[i];

        ta1_z_xxxxxx_xyy_0[i] = 5.0 * ta1_z_xxxx_xyy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xyy_1[i] * fe_0 + ta1_z_xxxxx_yy_0[i] * fe_0 -
                                ta1_z_xxxxx_yy_1[i] * fe_0 + ta1_z_xxxxx_xyy_0[i] * pa_x[i] - ta1_z_xxxxx_xyy_1[i] * pc_x[i];

        ta1_z_xxxxxx_xyz_0[i] = 5.0 * ta1_z_xxxx_xyz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xyz_1[i] * fe_0 + ta1_z_xxxxx_yz_0[i] * fe_0 -
                                ta1_z_xxxxx_yz_1[i] * fe_0 + ta1_z_xxxxx_xyz_0[i] * pa_x[i] - ta1_z_xxxxx_xyz_1[i] * pc_x[i];

        ta1_z_xxxxxx_xzz_0[i] = 5.0 * ta1_z_xxxx_xzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xzz_1[i] * fe_0 + ta1_z_xxxxx_zz_0[i] * fe_0 -
                                ta1_z_xxxxx_zz_1[i] * fe_0 + ta1_z_xxxxx_xzz_0[i] * pa_x[i] - ta1_z_xxxxx_xzz_1[i] * pc_x[i];

        ta1_z_xxxxxx_yyy_0[i] =
            5.0 * ta1_z_xxxx_yyy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yyy_1[i] * fe_0 + ta1_z_xxxxx_yyy_0[i] * pa_x[i] - ta1_z_xxxxx_yyy_1[i] * pc_x[i];

        ta1_z_xxxxxx_yyz_0[i] =
            5.0 * ta1_z_xxxx_yyz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yyz_1[i] * fe_0 + ta1_z_xxxxx_yyz_0[i] * pa_x[i] - ta1_z_xxxxx_yyz_1[i] * pc_x[i];

        ta1_z_xxxxxx_yzz_0[i] =
            5.0 * ta1_z_xxxx_yzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yzz_1[i] * fe_0 + ta1_z_xxxxx_yzz_0[i] * pa_x[i] - ta1_z_xxxxx_yzz_1[i] * pc_x[i];

        ta1_z_xxxxxx_zzz_0[i] =
            5.0 * ta1_z_xxxx_zzz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_zzz_1[i] * fe_0 + ta1_z_xxxxx_zzz_0[i] * pa_x[i] - ta1_z_xxxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 570-580 components of targeted buffer : IF

    auto ta1_z_xxxxxy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 570);

    auto ta1_z_xxxxxy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 571);

    auto ta1_z_xxxxxy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 572);

    auto ta1_z_xxxxxy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 573);

    auto ta1_z_xxxxxy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 574);

    auto ta1_z_xxxxxy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 575);

    auto ta1_z_xxxxxy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 576);

    auto ta1_z_xxxxxy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 577);

    auto ta1_z_xxxxxy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 578);

    auto ta1_z_xxxxxy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 579);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxxxx_xx_0,   \
                             ta1_z_xxxxx_xx_1,   \
                             ta1_z_xxxxx_xxx_0,  \
                             ta1_z_xxxxx_xxx_1,  \
                             ta1_z_xxxxx_xxy_0,  \
                             ta1_z_xxxxx_xxy_1,  \
                             ta1_z_xxxxx_xxz_0,  \
                             ta1_z_xxxxx_xxz_1,  \
                             ta1_z_xxxxx_xy_0,   \
                             ta1_z_xxxxx_xy_1,   \
                             ta1_z_xxxxx_xyy_0,  \
                             ta1_z_xxxxx_xyy_1,  \
                             ta1_z_xxxxx_xyz_0,  \
                             ta1_z_xxxxx_xyz_1,  \
                             ta1_z_xxxxx_xz_0,   \
                             ta1_z_xxxxx_xz_1,   \
                             ta1_z_xxxxx_xzz_0,  \
                             ta1_z_xxxxx_xzz_1,  \
                             ta1_z_xxxxx_zzz_0,  \
                             ta1_z_xxxxx_zzz_1,  \
                             ta1_z_xxxxxy_xxx_0, \
                             ta1_z_xxxxxy_xxy_0, \
                             ta1_z_xxxxxy_xxz_0, \
                             ta1_z_xxxxxy_xyy_0, \
                             ta1_z_xxxxxy_xyz_0, \
                             ta1_z_xxxxxy_xzz_0, \
                             ta1_z_xxxxxy_yyy_0, \
                             ta1_z_xxxxxy_yyz_0, \
                             ta1_z_xxxxxy_yzz_0, \
                             ta1_z_xxxxxy_zzz_0, \
                             ta1_z_xxxxy_yyy_0,  \
                             ta1_z_xxxxy_yyy_1,  \
                             ta1_z_xxxxy_yyz_0,  \
                             ta1_z_xxxxy_yyz_1,  \
                             ta1_z_xxxxy_yzz_0,  \
                             ta1_z_xxxxy_yzz_1,  \
                             ta1_z_xxxy_yyy_0,   \
                             ta1_z_xxxy_yyy_1,   \
                             ta1_z_xxxy_yyz_0,   \
                             ta1_z_xxxy_yyz_1,   \
                             ta1_z_xxxy_yzz_0,   \
                             ta1_z_xxxy_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxy_xxx_0[i] = ta1_z_xxxxx_xxx_0[i] * pa_y[i] - ta1_z_xxxxx_xxx_1[i] * pc_y[i];

        ta1_z_xxxxxy_xxy_0[i] =
            ta1_z_xxxxx_xx_0[i] * fe_0 - ta1_z_xxxxx_xx_1[i] * fe_0 + ta1_z_xxxxx_xxy_0[i] * pa_y[i] - ta1_z_xxxxx_xxy_1[i] * pc_y[i];

        ta1_z_xxxxxy_xxz_0[i] = ta1_z_xxxxx_xxz_0[i] * pa_y[i] - ta1_z_xxxxx_xxz_1[i] * pc_y[i];

        ta1_z_xxxxxy_xyy_0[i] =
            2.0 * ta1_z_xxxxx_xy_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xy_1[i] * fe_0 + ta1_z_xxxxx_xyy_0[i] * pa_y[i] - ta1_z_xxxxx_xyy_1[i] * pc_y[i];

        ta1_z_xxxxxy_xyz_0[i] =
            ta1_z_xxxxx_xz_0[i] * fe_0 - ta1_z_xxxxx_xz_1[i] * fe_0 + ta1_z_xxxxx_xyz_0[i] * pa_y[i] - ta1_z_xxxxx_xyz_1[i] * pc_y[i];

        ta1_z_xxxxxy_xzz_0[i] = ta1_z_xxxxx_xzz_0[i] * pa_y[i] - ta1_z_xxxxx_xzz_1[i] * pc_y[i];

        ta1_z_xxxxxy_yyy_0[i] =
            4.0 * ta1_z_xxxy_yyy_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yyy_1[i] * fe_0 + ta1_z_xxxxy_yyy_0[i] * pa_x[i] - ta1_z_xxxxy_yyy_1[i] * pc_x[i];

        ta1_z_xxxxxy_yyz_0[i] =
            4.0 * ta1_z_xxxy_yyz_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yyz_1[i] * fe_0 + ta1_z_xxxxy_yyz_0[i] * pa_x[i] - ta1_z_xxxxy_yyz_1[i] * pc_x[i];

        ta1_z_xxxxxy_yzz_0[i] =
            4.0 * ta1_z_xxxy_yzz_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yzz_1[i] * fe_0 + ta1_z_xxxxy_yzz_0[i] * pa_x[i] - ta1_z_xxxxy_yzz_1[i] * pc_x[i];

        ta1_z_xxxxxy_zzz_0[i] = ta1_z_xxxxx_zzz_0[i] * pa_y[i] - ta1_z_xxxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 580-590 components of targeted buffer : IF

    auto ta1_z_xxxxxz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 580);

    auto ta1_z_xxxxxz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 581);

    auto ta1_z_xxxxxz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 582);

    auto ta1_z_xxxxxz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 583);

    auto ta1_z_xxxxxz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 584);

    auto ta1_z_xxxxxz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 585);

    auto ta1_z_xxxxxz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 586);

    auto ta1_z_xxxxxz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 587);

    auto ta1_z_xxxxxz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 588);

    auto ta1_z_xxxxxz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 589);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxxxx_xx_0,   \
                             ta1_z_xxxxx_xx_1,   \
                             ta1_z_xxxxx_xxx_0,  \
                             ta1_z_xxxxx_xxx_1,  \
                             ta1_z_xxxxx_xxy_0,  \
                             ta1_z_xxxxx_xxy_1,  \
                             ta1_z_xxxxx_xxz_0,  \
                             ta1_z_xxxxx_xxz_1,  \
                             ta1_z_xxxxx_xy_0,   \
                             ta1_z_xxxxx_xy_1,   \
                             ta1_z_xxxxx_xyy_0,  \
                             ta1_z_xxxxx_xyy_1,  \
                             ta1_z_xxxxx_xyz_0,  \
                             ta1_z_xxxxx_xyz_1,  \
                             ta1_z_xxxxx_xz_0,   \
                             ta1_z_xxxxx_xz_1,   \
                             ta1_z_xxxxx_xzz_0,  \
                             ta1_z_xxxxx_xzz_1,  \
                             ta1_z_xxxxx_yyy_0,  \
                             ta1_z_xxxxx_yyy_1,  \
                             ta1_z_xxxxxz_xxx_0, \
                             ta1_z_xxxxxz_xxy_0, \
                             ta1_z_xxxxxz_xxz_0, \
                             ta1_z_xxxxxz_xyy_0, \
                             ta1_z_xxxxxz_xyz_0, \
                             ta1_z_xxxxxz_xzz_0, \
                             ta1_z_xxxxxz_yyy_0, \
                             ta1_z_xxxxxz_yyz_0, \
                             ta1_z_xxxxxz_yzz_0, \
                             ta1_z_xxxxxz_zzz_0, \
                             ta1_z_xxxxz_yyz_0,  \
                             ta1_z_xxxxz_yyz_1,  \
                             ta1_z_xxxxz_yzz_0,  \
                             ta1_z_xxxxz_yzz_1,  \
                             ta1_z_xxxxz_zzz_0,  \
                             ta1_z_xxxxz_zzz_1,  \
                             ta1_z_xxxz_yyz_0,   \
                             ta1_z_xxxz_yyz_1,   \
                             ta1_z_xxxz_yzz_0,   \
                             ta1_z_xxxz_yzz_1,   \
                             ta1_z_xxxz_zzz_0,   \
                             ta1_z_xxxz_zzz_1,   \
                             ta_xxxxx_xxx_1,     \
                             ta_xxxxx_xxy_1,     \
                             ta_xxxxx_xxz_1,     \
                             ta_xxxxx_xyy_1,     \
                             ta_xxxxx_xyz_1,     \
                             ta_xxxxx_xzz_1,     \
                             ta_xxxxx_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxz_xxx_0[i] = ta_xxxxx_xxx_1[i] + ta1_z_xxxxx_xxx_0[i] * pa_z[i] - ta1_z_xxxxx_xxx_1[i] * pc_z[i];

        ta1_z_xxxxxz_xxy_0[i] = ta_xxxxx_xxy_1[i] + ta1_z_xxxxx_xxy_0[i] * pa_z[i] - ta1_z_xxxxx_xxy_1[i] * pc_z[i];

        ta1_z_xxxxxz_xxz_0[i] = ta1_z_xxxxx_xx_0[i] * fe_0 - ta1_z_xxxxx_xx_1[i] * fe_0 + ta_xxxxx_xxz_1[i] + ta1_z_xxxxx_xxz_0[i] * pa_z[i] -
                                ta1_z_xxxxx_xxz_1[i] * pc_z[i];

        ta1_z_xxxxxz_xyy_0[i] = ta_xxxxx_xyy_1[i] + ta1_z_xxxxx_xyy_0[i] * pa_z[i] - ta1_z_xxxxx_xyy_1[i] * pc_z[i];

        ta1_z_xxxxxz_xyz_0[i] = ta1_z_xxxxx_xy_0[i] * fe_0 - ta1_z_xxxxx_xy_1[i] * fe_0 + ta_xxxxx_xyz_1[i] + ta1_z_xxxxx_xyz_0[i] * pa_z[i] -
                                ta1_z_xxxxx_xyz_1[i] * pc_z[i];

        ta1_z_xxxxxz_xzz_0[i] = 2.0 * ta1_z_xxxxx_xz_0[i] * fe_0 - 2.0 * ta1_z_xxxxx_xz_1[i] * fe_0 + ta_xxxxx_xzz_1[i] +
                                ta1_z_xxxxx_xzz_0[i] * pa_z[i] - ta1_z_xxxxx_xzz_1[i] * pc_z[i];

        ta1_z_xxxxxz_yyy_0[i] = ta_xxxxx_yyy_1[i] + ta1_z_xxxxx_yyy_0[i] * pa_z[i] - ta1_z_xxxxx_yyy_1[i] * pc_z[i];

        ta1_z_xxxxxz_yyz_0[i] =
            4.0 * ta1_z_xxxz_yyz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_yyz_1[i] * fe_0 + ta1_z_xxxxz_yyz_0[i] * pa_x[i] - ta1_z_xxxxz_yyz_1[i] * pc_x[i];

        ta1_z_xxxxxz_yzz_0[i] =
            4.0 * ta1_z_xxxz_yzz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_yzz_1[i] * fe_0 + ta1_z_xxxxz_yzz_0[i] * pa_x[i] - ta1_z_xxxxz_yzz_1[i] * pc_x[i];

        ta1_z_xxxxxz_zzz_0[i] =
            4.0 * ta1_z_xxxz_zzz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_zzz_1[i] * fe_0 + ta1_z_xxxxz_zzz_0[i] * pa_x[i] - ta1_z_xxxxz_zzz_1[i] * pc_x[i];
    }

    // Set up 590-600 components of targeted buffer : IF

    auto ta1_z_xxxxyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 590);

    auto ta1_z_xxxxyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 591);

    auto ta1_z_xxxxyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 592);

    auto ta1_z_xxxxyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 593);

    auto ta1_z_xxxxyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 594);

    auto ta1_z_xxxxyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 595);

    auto ta1_z_xxxxyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 596);

    auto ta1_z_xxxxyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 597);

    auto ta1_z_xxxxyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 598);

    auto ta1_z_xxxxyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 599);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxxx_xxx_0,   \
                             ta1_z_xxxx_xxx_1,   \
                             ta1_z_xxxx_xxz_0,   \
                             ta1_z_xxxx_xxz_1,   \
                             ta1_z_xxxx_xzz_0,   \
                             ta1_z_xxxx_xzz_1,   \
                             ta1_z_xxxxy_xxx_0,  \
                             ta1_z_xxxxy_xxx_1,  \
                             ta1_z_xxxxy_xxz_0,  \
                             ta1_z_xxxxy_xxz_1,  \
                             ta1_z_xxxxy_xzz_0,  \
                             ta1_z_xxxxy_xzz_1,  \
                             ta1_z_xxxxyy_xxx_0, \
                             ta1_z_xxxxyy_xxy_0, \
                             ta1_z_xxxxyy_xxz_0, \
                             ta1_z_xxxxyy_xyy_0, \
                             ta1_z_xxxxyy_xyz_0, \
                             ta1_z_xxxxyy_xzz_0, \
                             ta1_z_xxxxyy_yyy_0, \
                             ta1_z_xxxxyy_yyz_0, \
                             ta1_z_xxxxyy_yzz_0, \
                             ta1_z_xxxxyy_zzz_0, \
                             ta1_z_xxxyy_xxy_0,  \
                             ta1_z_xxxyy_xxy_1,  \
                             ta1_z_xxxyy_xy_0,   \
                             ta1_z_xxxyy_xy_1,   \
                             ta1_z_xxxyy_xyy_0,  \
                             ta1_z_xxxyy_xyy_1,  \
                             ta1_z_xxxyy_xyz_0,  \
                             ta1_z_xxxyy_xyz_1,  \
                             ta1_z_xxxyy_yy_0,   \
                             ta1_z_xxxyy_yy_1,   \
                             ta1_z_xxxyy_yyy_0,  \
                             ta1_z_xxxyy_yyy_1,  \
                             ta1_z_xxxyy_yyz_0,  \
                             ta1_z_xxxyy_yyz_1,  \
                             ta1_z_xxxyy_yz_0,   \
                             ta1_z_xxxyy_yz_1,   \
                             ta1_z_xxxyy_yzz_0,  \
                             ta1_z_xxxyy_yzz_1,  \
                             ta1_z_xxxyy_zzz_0,  \
                             ta1_z_xxxyy_zzz_1,  \
                             ta1_z_xxyy_xxy_0,   \
                             ta1_z_xxyy_xxy_1,   \
                             ta1_z_xxyy_xyy_0,   \
                             ta1_z_xxyy_xyy_1,   \
                             ta1_z_xxyy_xyz_0,   \
                             ta1_z_xxyy_xyz_1,   \
                             ta1_z_xxyy_yyy_0,   \
                             ta1_z_xxyy_yyy_1,   \
                             ta1_z_xxyy_yyz_0,   \
                             ta1_z_xxyy_yyz_1,   \
                             ta1_z_xxyy_yzz_0,   \
                             ta1_z_xxyy_yzz_1,   \
                             ta1_z_xxyy_zzz_0,   \
                             ta1_z_xxyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxyy_xxx_0[i] =
            ta1_z_xxxx_xxx_0[i] * fe_0 - ta1_z_xxxx_xxx_1[i] * fe_0 + ta1_z_xxxxy_xxx_0[i] * pa_y[i] - ta1_z_xxxxy_xxx_1[i] * pc_y[i];

        ta1_z_xxxxyy_xxy_0[i] = 3.0 * ta1_z_xxyy_xxy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxxyy_xy_0[i] * fe_0 -
                                2.0 * ta1_z_xxxyy_xy_1[i] * fe_0 + ta1_z_xxxyy_xxy_0[i] * pa_x[i] - ta1_z_xxxyy_xxy_1[i] * pc_x[i];

        ta1_z_xxxxyy_xxz_0[i] =
            ta1_z_xxxx_xxz_0[i] * fe_0 - ta1_z_xxxx_xxz_1[i] * fe_0 + ta1_z_xxxxy_xxz_0[i] * pa_y[i] - ta1_z_xxxxy_xxz_1[i] * pc_y[i];

        ta1_z_xxxxyy_xyy_0[i] = 3.0 * ta1_z_xxyy_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xyy_1[i] * fe_0 + ta1_z_xxxyy_yy_0[i] * fe_0 -
                                ta1_z_xxxyy_yy_1[i] * fe_0 + ta1_z_xxxyy_xyy_0[i] * pa_x[i] - ta1_z_xxxyy_xyy_1[i] * pc_x[i];

        ta1_z_xxxxyy_xyz_0[i] = 3.0 * ta1_z_xxyy_xyz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xyz_1[i] * fe_0 + ta1_z_xxxyy_yz_0[i] * fe_0 -
                                ta1_z_xxxyy_yz_1[i] * fe_0 + ta1_z_xxxyy_xyz_0[i] * pa_x[i] - ta1_z_xxxyy_xyz_1[i] * pc_x[i];

        ta1_z_xxxxyy_xzz_0[i] =
            ta1_z_xxxx_xzz_0[i] * fe_0 - ta1_z_xxxx_xzz_1[i] * fe_0 + ta1_z_xxxxy_xzz_0[i] * pa_y[i] - ta1_z_xxxxy_xzz_1[i] * pc_y[i];

        ta1_z_xxxxyy_yyy_0[i] =
            3.0 * ta1_z_xxyy_yyy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yyy_1[i] * fe_0 + ta1_z_xxxyy_yyy_0[i] * pa_x[i] - ta1_z_xxxyy_yyy_1[i] * pc_x[i];

        ta1_z_xxxxyy_yyz_0[i] =
            3.0 * ta1_z_xxyy_yyz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yyz_1[i] * fe_0 + ta1_z_xxxyy_yyz_0[i] * pa_x[i] - ta1_z_xxxyy_yyz_1[i] * pc_x[i];

        ta1_z_xxxxyy_yzz_0[i] =
            3.0 * ta1_z_xxyy_yzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yzz_1[i] * fe_0 + ta1_z_xxxyy_yzz_0[i] * pa_x[i] - ta1_z_xxxyy_yzz_1[i] * pc_x[i];

        ta1_z_xxxxyy_zzz_0[i] =
            3.0 * ta1_z_xxyy_zzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_zzz_1[i] * fe_0 + ta1_z_xxxyy_zzz_0[i] * pa_x[i] - ta1_z_xxxyy_zzz_1[i] * pc_x[i];
    }

    // Set up 600-610 components of targeted buffer : IF

    auto ta1_z_xxxxyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 600);

    auto ta1_z_xxxxyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 601);

    auto ta1_z_xxxxyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 602);

    auto ta1_z_xxxxyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 603);

    auto ta1_z_xxxxyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 604);

    auto ta1_z_xxxxyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 605);

    auto ta1_z_xxxxyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 606);

    auto ta1_z_xxxxyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 607);

    auto ta1_z_xxxxyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 608);

    auto ta1_z_xxxxyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 609);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxxxy_xxy_0,  \
                             ta1_z_xxxxy_xxy_1,  \
                             ta1_z_xxxxy_xyy_0,  \
                             ta1_z_xxxxy_xyy_1,  \
                             ta1_z_xxxxy_yyy_0,  \
                             ta1_z_xxxxy_yyy_1,  \
                             ta1_z_xxxxyz_xxx_0, \
                             ta1_z_xxxxyz_xxy_0, \
                             ta1_z_xxxxyz_xxz_0, \
                             ta1_z_xxxxyz_xyy_0, \
                             ta1_z_xxxxyz_xyz_0, \
                             ta1_z_xxxxyz_xzz_0, \
                             ta1_z_xxxxyz_yyy_0, \
                             ta1_z_xxxxyz_yyz_0, \
                             ta1_z_xxxxyz_yzz_0, \
                             ta1_z_xxxxyz_zzz_0, \
                             ta1_z_xxxxz_xxx_0,  \
                             ta1_z_xxxxz_xxx_1,  \
                             ta1_z_xxxxz_xxz_0,  \
                             ta1_z_xxxxz_xxz_1,  \
                             ta1_z_xxxxz_xyz_0,  \
                             ta1_z_xxxxz_xyz_1,  \
                             ta1_z_xxxxz_xz_0,   \
                             ta1_z_xxxxz_xz_1,   \
                             ta1_z_xxxxz_xzz_0,  \
                             ta1_z_xxxxz_xzz_1,  \
                             ta1_z_xxxxz_zzz_0,  \
                             ta1_z_xxxxz_zzz_1,  \
                             ta1_z_xxxyz_yyz_0,  \
                             ta1_z_xxxyz_yyz_1,  \
                             ta1_z_xxxyz_yzz_0,  \
                             ta1_z_xxxyz_yzz_1,  \
                             ta1_z_xxyz_yyz_0,   \
                             ta1_z_xxyz_yyz_1,   \
                             ta1_z_xxyz_yzz_0,   \
                             ta1_z_xxyz_yzz_1,   \
                             ta_xxxxy_xxy_1,     \
                             ta_xxxxy_xyy_1,     \
                             ta_xxxxy_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxyz_xxx_0[i] = ta1_z_xxxxz_xxx_0[i] * pa_y[i] - ta1_z_xxxxz_xxx_1[i] * pc_y[i];

        ta1_z_xxxxyz_xxy_0[i] = ta_xxxxy_xxy_1[i] + ta1_z_xxxxy_xxy_0[i] * pa_z[i] - ta1_z_xxxxy_xxy_1[i] * pc_z[i];

        ta1_z_xxxxyz_xxz_0[i] = ta1_z_xxxxz_xxz_0[i] * pa_y[i] - ta1_z_xxxxz_xxz_1[i] * pc_y[i];

        ta1_z_xxxxyz_xyy_0[i] = ta_xxxxy_xyy_1[i] + ta1_z_xxxxy_xyy_0[i] * pa_z[i] - ta1_z_xxxxy_xyy_1[i] * pc_z[i];

        ta1_z_xxxxyz_xyz_0[i] =
            ta1_z_xxxxz_xz_0[i] * fe_0 - ta1_z_xxxxz_xz_1[i] * fe_0 + ta1_z_xxxxz_xyz_0[i] * pa_y[i] - ta1_z_xxxxz_xyz_1[i] * pc_y[i];

        ta1_z_xxxxyz_xzz_0[i] = ta1_z_xxxxz_xzz_0[i] * pa_y[i] - ta1_z_xxxxz_xzz_1[i] * pc_y[i];

        ta1_z_xxxxyz_yyy_0[i] = ta_xxxxy_yyy_1[i] + ta1_z_xxxxy_yyy_0[i] * pa_z[i] - ta1_z_xxxxy_yyy_1[i] * pc_z[i];

        ta1_z_xxxxyz_yyz_0[i] =
            3.0 * ta1_z_xxyz_yyz_0[i] * fe_0 - 3.0 * ta1_z_xxyz_yyz_1[i] * fe_0 + ta1_z_xxxyz_yyz_0[i] * pa_x[i] - ta1_z_xxxyz_yyz_1[i] * pc_x[i];

        ta1_z_xxxxyz_yzz_0[i] =
            3.0 * ta1_z_xxyz_yzz_0[i] * fe_0 - 3.0 * ta1_z_xxyz_yzz_1[i] * fe_0 + ta1_z_xxxyz_yzz_0[i] * pa_x[i] - ta1_z_xxxyz_yzz_1[i] * pc_x[i];

        ta1_z_xxxxyz_zzz_0[i] = ta1_z_xxxxz_zzz_0[i] * pa_y[i] - ta1_z_xxxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 610-620 components of targeted buffer : IF

    auto ta1_z_xxxxzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 610);

    auto ta1_z_xxxxzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 611);

    auto ta1_z_xxxxzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 612);

    auto ta1_z_xxxxzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 613);

    auto ta1_z_xxxxzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 614);

    auto ta1_z_xxxxzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 615);

    auto ta1_z_xxxxzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 616);

    auto ta1_z_xxxxzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 617);

    auto ta1_z_xxxxzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 618);

    auto ta1_z_xxxxzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 619);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxxx_xxx_0,   \
                             ta1_z_xxxx_xxx_1,   \
                             ta1_z_xxxx_xxy_0,   \
                             ta1_z_xxxx_xxy_1,   \
                             ta1_z_xxxx_xyy_0,   \
                             ta1_z_xxxx_xyy_1,   \
                             ta1_z_xxxxz_xxx_0,  \
                             ta1_z_xxxxz_xxx_1,  \
                             ta1_z_xxxxz_xxy_0,  \
                             ta1_z_xxxxz_xxy_1,  \
                             ta1_z_xxxxz_xyy_0,  \
                             ta1_z_xxxxz_xyy_1,  \
                             ta1_z_xxxxzz_xxx_0, \
                             ta1_z_xxxxzz_xxy_0, \
                             ta1_z_xxxxzz_xxz_0, \
                             ta1_z_xxxxzz_xyy_0, \
                             ta1_z_xxxxzz_xyz_0, \
                             ta1_z_xxxxzz_xzz_0, \
                             ta1_z_xxxxzz_yyy_0, \
                             ta1_z_xxxxzz_yyz_0, \
                             ta1_z_xxxxzz_yzz_0, \
                             ta1_z_xxxxzz_zzz_0, \
                             ta1_z_xxxzz_xxz_0,  \
                             ta1_z_xxxzz_xxz_1,  \
                             ta1_z_xxxzz_xyz_0,  \
                             ta1_z_xxxzz_xyz_1,  \
                             ta1_z_xxxzz_xz_0,   \
                             ta1_z_xxxzz_xz_1,   \
                             ta1_z_xxxzz_xzz_0,  \
                             ta1_z_xxxzz_xzz_1,  \
                             ta1_z_xxxzz_yyy_0,  \
                             ta1_z_xxxzz_yyy_1,  \
                             ta1_z_xxxzz_yyz_0,  \
                             ta1_z_xxxzz_yyz_1,  \
                             ta1_z_xxxzz_yz_0,   \
                             ta1_z_xxxzz_yz_1,   \
                             ta1_z_xxxzz_yzz_0,  \
                             ta1_z_xxxzz_yzz_1,  \
                             ta1_z_xxxzz_zz_0,   \
                             ta1_z_xxxzz_zz_1,   \
                             ta1_z_xxxzz_zzz_0,  \
                             ta1_z_xxxzz_zzz_1,  \
                             ta1_z_xxzz_xxz_0,   \
                             ta1_z_xxzz_xxz_1,   \
                             ta1_z_xxzz_xyz_0,   \
                             ta1_z_xxzz_xyz_1,   \
                             ta1_z_xxzz_xzz_0,   \
                             ta1_z_xxzz_xzz_1,   \
                             ta1_z_xxzz_yyy_0,   \
                             ta1_z_xxzz_yyy_1,   \
                             ta1_z_xxzz_yyz_0,   \
                             ta1_z_xxzz_yyz_1,   \
                             ta1_z_xxzz_yzz_0,   \
                             ta1_z_xxzz_yzz_1,   \
                             ta1_z_xxzz_zzz_0,   \
                             ta1_z_xxzz_zzz_1,   \
                             ta_xxxxz_xxx_1,     \
                             ta_xxxxz_xxy_1,     \
                             ta_xxxxz_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxzz_xxx_0[i] = ta1_z_xxxx_xxx_0[i] * fe_0 - ta1_z_xxxx_xxx_1[i] * fe_0 + ta_xxxxz_xxx_1[i] + ta1_z_xxxxz_xxx_0[i] * pa_z[i] -
                                ta1_z_xxxxz_xxx_1[i] * pc_z[i];

        ta1_z_xxxxzz_xxy_0[i] = ta1_z_xxxx_xxy_0[i] * fe_0 - ta1_z_xxxx_xxy_1[i] * fe_0 + ta_xxxxz_xxy_1[i] + ta1_z_xxxxz_xxy_0[i] * pa_z[i] -
                                ta1_z_xxxxz_xxy_1[i] * pc_z[i];

        ta1_z_xxxxzz_xxz_0[i] = 3.0 * ta1_z_xxzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxz_1[i] * fe_0 + 2.0 * ta1_z_xxxzz_xz_0[i] * fe_0 -
                                2.0 * ta1_z_xxxzz_xz_1[i] * fe_0 + ta1_z_xxxzz_xxz_0[i] * pa_x[i] - ta1_z_xxxzz_xxz_1[i] * pc_x[i];

        ta1_z_xxxxzz_xyy_0[i] = ta1_z_xxxx_xyy_0[i] * fe_0 - ta1_z_xxxx_xyy_1[i] * fe_0 + ta_xxxxz_xyy_1[i] + ta1_z_xxxxz_xyy_0[i] * pa_z[i] -
                                ta1_z_xxxxz_xyy_1[i] * pc_z[i];

        ta1_z_xxxxzz_xyz_0[i] = 3.0 * ta1_z_xxzz_xyz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xyz_1[i] * fe_0 + ta1_z_xxxzz_yz_0[i] * fe_0 -
                                ta1_z_xxxzz_yz_1[i] * fe_0 + ta1_z_xxxzz_xyz_0[i] * pa_x[i] - ta1_z_xxxzz_xyz_1[i] * pc_x[i];

        ta1_z_xxxxzz_xzz_0[i] = 3.0 * ta1_z_xxzz_xzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xzz_1[i] * fe_0 + ta1_z_xxxzz_zz_0[i] * fe_0 -
                                ta1_z_xxxzz_zz_1[i] * fe_0 + ta1_z_xxxzz_xzz_0[i] * pa_x[i] - ta1_z_xxxzz_xzz_1[i] * pc_x[i];

        ta1_z_xxxxzz_yyy_0[i] =
            3.0 * ta1_z_xxzz_yyy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yyy_1[i] * fe_0 + ta1_z_xxxzz_yyy_0[i] * pa_x[i] - ta1_z_xxxzz_yyy_1[i] * pc_x[i];

        ta1_z_xxxxzz_yyz_0[i] =
            3.0 * ta1_z_xxzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yyz_1[i] * fe_0 + ta1_z_xxxzz_yyz_0[i] * pa_x[i] - ta1_z_xxxzz_yyz_1[i] * pc_x[i];

        ta1_z_xxxxzz_yzz_0[i] =
            3.0 * ta1_z_xxzz_yzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yzz_1[i] * fe_0 + ta1_z_xxxzz_yzz_0[i] * pa_x[i] - ta1_z_xxxzz_yzz_1[i] * pc_x[i];

        ta1_z_xxxxzz_zzz_0[i] =
            3.0 * ta1_z_xxzz_zzz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_zzz_1[i] * fe_0 + ta1_z_xxxzz_zzz_0[i] * pa_x[i] - ta1_z_xxxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 620-630 components of targeted buffer : IF

    auto ta1_z_xxxyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 620);

    auto ta1_z_xxxyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 621);

    auto ta1_z_xxxyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 622);

    auto ta1_z_xxxyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 623);

    auto ta1_z_xxxyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 624);

    auto ta1_z_xxxyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 625);

    auto ta1_z_xxxyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 626);

    auto ta1_z_xxxyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 627);

    auto ta1_z_xxxyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 628);

    auto ta1_z_xxxyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 629);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxxy_xxx_0,   \
                             ta1_z_xxxy_xxx_1,   \
                             ta1_z_xxxy_xxz_0,   \
                             ta1_z_xxxy_xxz_1,   \
                             ta1_z_xxxy_xzz_0,   \
                             ta1_z_xxxy_xzz_1,   \
                             ta1_z_xxxyy_xxx_0,  \
                             ta1_z_xxxyy_xxx_1,  \
                             ta1_z_xxxyy_xxz_0,  \
                             ta1_z_xxxyy_xxz_1,  \
                             ta1_z_xxxyy_xzz_0,  \
                             ta1_z_xxxyy_xzz_1,  \
                             ta1_z_xxxyyy_xxx_0, \
                             ta1_z_xxxyyy_xxy_0, \
                             ta1_z_xxxyyy_xxz_0, \
                             ta1_z_xxxyyy_xyy_0, \
                             ta1_z_xxxyyy_xyz_0, \
                             ta1_z_xxxyyy_xzz_0, \
                             ta1_z_xxxyyy_yyy_0, \
                             ta1_z_xxxyyy_yyz_0, \
                             ta1_z_xxxyyy_yzz_0, \
                             ta1_z_xxxyyy_zzz_0, \
                             ta1_z_xxyyy_xxy_0,  \
                             ta1_z_xxyyy_xxy_1,  \
                             ta1_z_xxyyy_xy_0,   \
                             ta1_z_xxyyy_xy_1,   \
                             ta1_z_xxyyy_xyy_0,  \
                             ta1_z_xxyyy_xyy_1,  \
                             ta1_z_xxyyy_xyz_0,  \
                             ta1_z_xxyyy_xyz_1,  \
                             ta1_z_xxyyy_yy_0,   \
                             ta1_z_xxyyy_yy_1,   \
                             ta1_z_xxyyy_yyy_0,  \
                             ta1_z_xxyyy_yyy_1,  \
                             ta1_z_xxyyy_yyz_0,  \
                             ta1_z_xxyyy_yyz_1,  \
                             ta1_z_xxyyy_yz_0,   \
                             ta1_z_xxyyy_yz_1,   \
                             ta1_z_xxyyy_yzz_0,  \
                             ta1_z_xxyyy_yzz_1,  \
                             ta1_z_xxyyy_zzz_0,  \
                             ta1_z_xxyyy_zzz_1,  \
                             ta1_z_xyyy_xxy_0,   \
                             ta1_z_xyyy_xxy_1,   \
                             ta1_z_xyyy_xyy_0,   \
                             ta1_z_xyyy_xyy_1,   \
                             ta1_z_xyyy_xyz_0,   \
                             ta1_z_xyyy_xyz_1,   \
                             ta1_z_xyyy_yyy_0,   \
                             ta1_z_xyyy_yyy_1,   \
                             ta1_z_xyyy_yyz_0,   \
                             ta1_z_xyyy_yyz_1,   \
                             ta1_z_xyyy_yzz_0,   \
                             ta1_z_xyyy_yzz_1,   \
                             ta1_z_xyyy_zzz_0,   \
                             ta1_z_xyyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyy_xxx_0[i] =
            2.0 * ta1_z_xxxy_xxx_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xxx_1[i] * fe_0 + ta1_z_xxxyy_xxx_0[i] * pa_y[i] - ta1_z_xxxyy_xxx_1[i] * pc_y[i];

        ta1_z_xxxyyy_xxy_0[i] = 2.0 * ta1_z_xyyy_xxy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_z_xxyyy_xy_1[i] * fe_0 + ta1_z_xxyyy_xxy_0[i] * pa_x[i] - ta1_z_xxyyy_xxy_1[i] * pc_x[i];

        ta1_z_xxxyyy_xxz_0[i] =
            2.0 * ta1_z_xxxy_xxz_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xxz_1[i] * fe_0 + ta1_z_xxxyy_xxz_0[i] * pa_y[i] - ta1_z_xxxyy_xxz_1[i] * pc_y[i];

        ta1_z_xxxyyy_xyy_0[i] = 2.0 * ta1_z_xyyy_xyy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xyy_1[i] * fe_0 + ta1_z_xxyyy_yy_0[i] * fe_0 -
                                ta1_z_xxyyy_yy_1[i] * fe_0 + ta1_z_xxyyy_xyy_0[i] * pa_x[i] - ta1_z_xxyyy_xyy_1[i] * pc_x[i];

        ta1_z_xxxyyy_xyz_0[i] = 2.0 * ta1_z_xyyy_xyz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xyz_1[i] * fe_0 + ta1_z_xxyyy_yz_0[i] * fe_0 -
                                ta1_z_xxyyy_yz_1[i] * fe_0 + ta1_z_xxyyy_xyz_0[i] * pa_x[i] - ta1_z_xxyyy_xyz_1[i] * pc_x[i];

        ta1_z_xxxyyy_xzz_0[i] =
            2.0 * ta1_z_xxxy_xzz_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xzz_1[i] * fe_0 + ta1_z_xxxyy_xzz_0[i] * pa_y[i] - ta1_z_xxxyy_xzz_1[i] * pc_y[i];

        ta1_z_xxxyyy_yyy_0[i] =
            2.0 * ta1_z_xyyy_yyy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yyy_1[i] * fe_0 + ta1_z_xxyyy_yyy_0[i] * pa_x[i] - ta1_z_xxyyy_yyy_1[i] * pc_x[i];

        ta1_z_xxxyyy_yyz_0[i] =
            2.0 * ta1_z_xyyy_yyz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yyz_1[i] * fe_0 + ta1_z_xxyyy_yyz_0[i] * pa_x[i] - ta1_z_xxyyy_yyz_1[i] * pc_x[i];

        ta1_z_xxxyyy_yzz_0[i] =
            2.0 * ta1_z_xyyy_yzz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yzz_1[i] * fe_0 + ta1_z_xxyyy_yzz_0[i] * pa_x[i] - ta1_z_xxyyy_yzz_1[i] * pc_x[i];

        ta1_z_xxxyyy_zzz_0[i] =
            2.0 * ta1_z_xyyy_zzz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_zzz_1[i] * fe_0 + ta1_z_xxyyy_zzz_0[i] * pa_x[i] - ta1_z_xxyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 630-640 components of targeted buffer : IF

    auto ta1_z_xxxyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 630);

    auto ta1_z_xxxyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 631);

    auto ta1_z_xxxyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 632);

    auto ta1_z_xxxyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 633);

    auto ta1_z_xxxyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 634);

    auto ta1_z_xxxyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 635);

    auto ta1_z_xxxyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 636);

    auto ta1_z_xxxyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 637);

    auto ta1_z_xxxyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 638);

    auto ta1_z_xxxyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 639);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxxyy_xxx_0,  \
                             ta1_z_xxxyy_xxx_1,  \
                             ta1_z_xxxyy_xxy_0,  \
                             ta1_z_xxxyy_xxy_1,  \
                             ta1_z_xxxyy_xy_0,   \
                             ta1_z_xxxyy_xy_1,   \
                             ta1_z_xxxyy_xyy_0,  \
                             ta1_z_xxxyy_xyy_1,  \
                             ta1_z_xxxyy_xyz_0,  \
                             ta1_z_xxxyy_xyz_1,  \
                             ta1_z_xxxyy_yyy_0,  \
                             ta1_z_xxxyy_yyy_1,  \
                             ta1_z_xxxyyz_xxx_0, \
                             ta1_z_xxxyyz_xxy_0, \
                             ta1_z_xxxyyz_xxz_0, \
                             ta1_z_xxxyyz_xyy_0, \
                             ta1_z_xxxyyz_xyz_0, \
                             ta1_z_xxxyyz_xzz_0, \
                             ta1_z_xxxyyz_yyy_0, \
                             ta1_z_xxxyyz_yyz_0, \
                             ta1_z_xxxyyz_yzz_0, \
                             ta1_z_xxxyyz_zzz_0, \
                             ta1_z_xxxyz_xxz_0,  \
                             ta1_z_xxxyz_xxz_1,  \
                             ta1_z_xxxyz_xzz_0,  \
                             ta1_z_xxxyz_xzz_1,  \
                             ta1_z_xxxz_xxz_0,   \
                             ta1_z_xxxz_xxz_1,   \
                             ta1_z_xxxz_xzz_0,   \
                             ta1_z_xxxz_xzz_1,   \
                             ta1_z_xxyyz_yyz_0,  \
                             ta1_z_xxyyz_yyz_1,  \
                             ta1_z_xxyyz_yzz_0,  \
                             ta1_z_xxyyz_yzz_1,  \
                             ta1_z_xxyyz_zzz_0,  \
                             ta1_z_xxyyz_zzz_1,  \
                             ta1_z_xyyz_yyz_0,   \
                             ta1_z_xyyz_yyz_1,   \
                             ta1_z_xyyz_yzz_0,   \
                             ta1_z_xyyz_yzz_1,   \
                             ta1_z_xyyz_zzz_0,   \
                             ta1_z_xyyz_zzz_1,   \
                             ta_xxxyy_xxx_1,     \
                             ta_xxxyy_xxy_1,     \
                             ta_xxxyy_xyy_1,     \
                             ta_xxxyy_xyz_1,     \
                             ta_xxxyy_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyz_xxx_0[i] = ta_xxxyy_xxx_1[i] + ta1_z_xxxyy_xxx_0[i] * pa_z[i] - ta1_z_xxxyy_xxx_1[i] * pc_z[i];

        ta1_z_xxxyyz_xxy_0[i] = ta_xxxyy_xxy_1[i] + ta1_z_xxxyy_xxy_0[i] * pa_z[i] - ta1_z_xxxyy_xxy_1[i] * pc_z[i];

        ta1_z_xxxyyz_xxz_0[i] =
            ta1_z_xxxz_xxz_0[i] * fe_0 - ta1_z_xxxz_xxz_1[i] * fe_0 + ta1_z_xxxyz_xxz_0[i] * pa_y[i] - ta1_z_xxxyz_xxz_1[i] * pc_y[i];

        ta1_z_xxxyyz_xyy_0[i] = ta_xxxyy_xyy_1[i] + ta1_z_xxxyy_xyy_0[i] * pa_z[i] - ta1_z_xxxyy_xyy_1[i] * pc_z[i];

        ta1_z_xxxyyz_xyz_0[i] = ta1_z_xxxyy_xy_0[i] * fe_0 - ta1_z_xxxyy_xy_1[i] * fe_0 + ta_xxxyy_xyz_1[i] + ta1_z_xxxyy_xyz_0[i] * pa_z[i] -
                                ta1_z_xxxyy_xyz_1[i] * pc_z[i];

        ta1_z_xxxyyz_xzz_0[i] =
            ta1_z_xxxz_xzz_0[i] * fe_0 - ta1_z_xxxz_xzz_1[i] * fe_0 + ta1_z_xxxyz_xzz_0[i] * pa_y[i] - ta1_z_xxxyz_xzz_1[i] * pc_y[i];

        ta1_z_xxxyyz_yyy_0[i] = ta_xxxyy_yyy_1[i] + ta1_z_xxxyy_yyy_0[i] * pa_z[i] - ta1_z_xxxyy_yyy_1[i] * pc_z[i];

        ta1_z_xxxyyz_yyz_0[i] =
            2.0 * ta1_z_xyyz_yyz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_yyz_1[i] * fe_0 + ta1_z_xxyyz_yyz_0[i] * pa_x[i] - ta1_z_xxyyz_yyz_1[i] * pc_x[i];

        ta1_z_xxxyyz_yzz_0[i] =
            2.0 * ta1_z_xyyz_yzz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_yzz_1[i] * fe_0 + ta1_z_xxyyz_yzz_0[i] * pa_x[i] - ta1_z_xxyyz_yzz_1[i] * pc_x[i];

        ta1_z_xxxyyz_zzz_0[i] =
            2.0 * ta1_z_xyyz_zzz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_zzz_1[i] * fe_0 + ta1_z_xxyyz_zzz_0[i] * pa_x[i] - ta1_z_xxyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 640-650 components of targeted buffer : IF

    auto ta1_z_xxxyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 640);

    auto ta1_z_xxxyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 641);

    auto ta1_z_xxxyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 642);

    auto ta1_z_xxxyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 643);

    auto ta1_z_xxxyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 644);

    auto ta1_z_xxxyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 645);

    auto ta1_z_xxxyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 646);

    auto ta1_z_xxxyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 647);

    auto ta1_z_xxxyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 648);

    auto ta1_z_xxxyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 649);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxxyzz_xxx_0, \
                             ta1_z_xxxyzz_xxy_0, \
                             ta1_z_xxxyzz_xxz_0, \
                             ta1_z_xxxyzz_xyy_0, \
                             ta1_z_xxxyzz_xyz_0, \
                             ta1_z_xxxyzz_xzz_0, \
                             ta1_z_xxxyzz_yyy_0, \
                             ta1_z_xxxyzz_yyz_0, \
                             ta1_z_xxxyzz_yzz_0, \
                             ta1_z_xxxyzz_zzz_0, \
                             ta1_z_xxxzz_xx_0,   \
                             ta1_z_xxxzz_xx_1,   \
                             ta1_z_xxxzz_xxx_0,  \
                             ta1_z_xxxzz_xxx_1,  \
                             ta1_z_xxxzz_xxy_0,  \
                             ta1_z_xxxzz_xxy_1,  \
                             ta1_z_xxxzz_xxz_0,  \
                             ta1_z_xxxzz_xxz_1,  \
                             ta1_z_xxxzz_xy_0,   \
                             ta1_z_xxxzz_xy_1,   \
                             ta1_z_xxxzz_xyy_0,  \
                             ta1_z_xxxzz_xyy_1,  \
                             ta1_z_xxxzz_xyz_0,  \
                             ta1_z_xxxzz_xyz_1,  \
                             ta1_z_xxxzz_xz_0,   \
                             ta1_z_xxxzz_xz_1,   \
                             ta1_z_xxxzz_xzz_0,  \
                             ta1_z_xxxzz_xzz_1,  \
                             ta1_z_xxxzz_zzz_0,  \
                             ta1_z_xxxzz_zzz_1,  \
                             ta1_z_xxyzz_yyy_0,  \
                             ta1_z_xxyzz_yyy_1,  \
                             ta1_z_xxyzz_yyz_0,  \
                             ta1_z_xxyzz_yyz_1,  \
                             ta1_z_xxyzz_yzz_0,  \
                             ta1_z_xxyzz_yzz_1,  \
                             ta1_z_xyzz_yyy_0,   \
                             ta1_z_xyzz_yyy_1,   \
                             ta1_z_xyzz_yyz_0,   \
                             ta1_z_xyzz_yyz_1,   \
                             ta1_z_xyzz_yzz_0,   \
                             ta1_z_xyzz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyzz_xxx_0[i] = ta1_z_xxxzz_xxx_0[i] * pa_y[i] - ta1_z_xxxzz_xxx_1[i] * pc_y[i];

        ta1_z_xxxyzz_xxy_0[i] =
            ta1_z_xxxzz_xx_0[i] * fe_0 - ta1_z_xxxzz_xx_1[i] * fe_0 + ta1_z_xxxzz_xxy_0[i] * pa_y[i] - ta1_z_xxxzz_xxy_1[i] * pc_y[i];

        ta1_z_xxxyzz_xxz_0[i] = ta1_z_xxxzz_xxz_0[i] * pa_y[i] - ta1_z_xxxzz_xxz_1[i] * pc_y[i];

        ta1_z_xxxyzz_xyy_0[i] =
            2.0 * ta1_z_xxxzz_xy_0[i] * fe_0 - 2.0 * ta1_z_xxxzz_xy_1[i] * fe_0 + ta1_z_xxxzz_xyy_0[i] * pa_y[i] - ta1_z_xxxzz_xyy_1[i] * pc_y[i];

        ta1_z_xxxyzz_xyz_0[i] =
            ta1_z_xxxzz_xz_0[i] * fe_0 - ta1_z_xxxzz_xz_1[i] * fe_0 + ta1_z_xxxzz_xyz_0[i] * pa_y[i] - ta1_z_xxxzz_xyz_1[i] * pc_y[i];

        ta1_z_xxxyzz_xzz_0[i] = ta1_z_xxxzz_xzz_0[i] * pa_y[i] - ta1_z_xxxzz_xzz_1[i] * pc_y[i];

        ta1_z_xxxyzz_yyy_0[i] =
            2.0 * ta1_z_xyzz_yyy_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yyy_1[i] * fe_0 + ta1_z_xxyzz_yyy_0[i] * pa_x[i] - ta1_z_xxyzz_yyy_1[i] * pc_x[i];

        ta1_z_xxxyzz_yyz_0[i] =
            2.0 * ta1_z_xyzz_yyz_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yyz_1[i] * fe_0 + ta1_z_xxyzz_yyz_0[i] * pa_x[i] - ta1_z_xxyzz_yyz_1[i] * pc_x[i];

        ta1_z_xxxyzz_yzz_0[i] =
            2.0 * ta1_z_xyzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yzz_1[i] * fe_0 + ta1_z_xxyzz_yzz_0[i] * pa_x[i] - ta1_z_xxyzz_yzz_1[i] * pc_x[i];

        ta1_z_xxxyzz_zzz_0[i] = ta1_z_xxxzz_zzz_0[i] * pa_y[i] - ta1_z_xxxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 650-660 components of targeted buffer : IF

    auto ta1_z_xxxzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 650);

    auto ta1_z_xxxzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 651);

    auto ta1_z_xxxzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 652);

    auto ta1_z_xxxzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 653);

    auto ta1_z_xxxzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 654);

    auto ta1_z_xxxzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 655);

    auto ta1_z_xxxzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 656);

    auto ta1_z_xxxzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 657);

    auto ta1_z_xxxzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 658);

    auto ta1_z_xxxzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 659);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxxz_xxx_0,   \
                             ta1_z_xxxz_xxx_1,   \
                             ta1_z_xxxz_xxy_0,   \
                             ta1_z_xxxz_xxy_1,   \
                             ta1_z_xxxz_xyy_0,   \
                             ta1_z_xxxz_xyy_1,   \
                             ta1_z_xxxzz_xxx_0,  \
                             ta1_z_xxxzz_xxx_1,  \
                             ta1_z_xxxzz_xxy_0,  \
                             ta1_z_xxxzz_xxy_1,  \
                             ta1_z_xxxzz_xyy_0,  \
                             ta1_z_xxxzz_xyy_1,  \
                             ta1_z_xxxzzz_xxx_0, \
                             ta1_z_xxxzzz_xxy_0, \
                             ta1_z_xxxzzz_xxz_0, \
                             ta1_z_xxxzzz_xyy_0, \
                             ta1_z_xxxzzz_xyz_0, \
                             ta1_z_xxxzzz_xzz_0, \
                             ta1_z_xxxzzz_yyy_0, \
                             ta1_z_xxxzzz_yyz_0, \
                             ta1_z_xxxzzz_yzz_0, \
                             ta1_z_xxxzzz_zzz_0, \
                             ta1_z_xxzzz_xxz_0,  \
                             ta1_z_xxzzz_xxz_1,  \
                             ta1_z_xxzzz_xyz_0,  \
                             ta1_z_xxzzz_xyz_1,  \
                             ta1_z_xxzzz_xz_0,   \
                             ta1_z_xxzzz_xz_1,   \
                             ta1_z_xxzzz_xzz_0,  \
                             ta1_z_xxzzz_xzz_1,  \
                             ta1_z_xxzzz_yyy_0,  \
                             ta1_z_xxzzz_yyy_1,  \
                             ta1_z_xxzzz_yyz_0,  \
                             ta1_z_xxzzz_yyz_1,  \
                             ta1_z_xxzzz_yz_0,   \
                             ta1_z_xxzzz_yz_1,   \
                             ta1_z_xxzzz_yzz_0,  \
                             ta1_z_xxzzz_yzz_1,  \
                             ta1_z_xxzzz_zz_0,   \
                             ta1_z_xxzzz_zz_1,   \
                             ta1_z_xxzzz_zzz_0,  \
                             ta1_z_xxzzz_zzz_1,  \
                             ta1_z_xzzz_xxz_0,   \
                             ta1_z_xzzz_xxz_1,   \
                             ta1_z_xzzz_xyz_0,   \
                             ta1_z_xzzz_xyz_1,   \
                             ta1_z_xzzz_xzz_0,   \
                             ta1_z_xzzz_xzz_1,   \
                             ta1_z_xzzz_yyy_0,   \
                             ta1_z_xzzz_yyy_1,   \
                             ta1_z_xzzz_yyz_0,   \
                             ta1_z_xzzz_yyz_1,   \
                             ta1_z_xzzz_yzz_0,   \
                             ta1_z_xzzz_yzz_1,   \
                             ta1_z_xzzz_zzz_0,   \
                             ta1_z_xzzz_zzz_1,   \
                             ta_xxxzz_xxx_1,     \
                             ta_xxxzz_xxy_1,     \
                             ta_xxxzz_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzzz_xxx_0[i] = 2.0 * ta1_z_xxxz_xxx_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xxx_1[i] * fe_0 + ta_xxxzz_xxx_1[i] +
                                ta1_z_xxxzz_xxx_0[i] * pa_z[i] - ta1_z_xxxzz_xxx_1[i] * pc_z[i];

        ta1_z_xxxzzz_xxy_0[i] = 2.0 * ta1_z_xxxz_xxy_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xxy_1[i] * fe_0 + ta_xxxzz_xxy_1[i] +
                                ta1_z_xxxzz_xxy_0[i] * pa_z[i] - ta1_z_xxxzz_xxy_1[i] * pc_z[i];

        ta1_z_xxxzzz_xxz_0[i] = 2.0 * ta1_z_xzzz_xxz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xxz_1[i] * fe_0 + 2.0 * ta1_z_xxzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_z_xxzzz_xz_1[i] * fe_0 + ta1_z_xxzzz_xxz_0[i] * pa_x[i] - ta1_z_xxzzz_xxz_1[i] * pc_x[i];

        ta1_z_xxxzzz_xyy_0[i] = 2.0 * ta1_z_xxxz_xyy_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xyy_1[i] * fe_0 + ta_xxxzz_xyy_1[i] +
                                ta1_z_xxxzz_xyy_0[i] * pa_z[i] - ta1_z_xxxzz_xyy_1[i] * pc_z[i];

        ta1_z_xxxzzz_xyz_0[i] = 2.0 * ta1_z_xzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xyz_1[i] * fe_0 + ta1_z_xxzzz_yz_0[i] * fe_0 -
                                ta1_z_xxzzz_yz_1[i] * fe_0 + ta1_z_xxzzz_xyz_0[i] * pa_x[i] - ta1_z_xxzzz_xyz_1[i] * pc_x[i];

        ta1_z_xxxzzz_xzz_0[i] = 2.0 * ta1_z_xzzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xzz_1[i] * fe_0 + ta1_z_xxzzz_zz_0[i] * fe_0 -
                                ta1_z_xxzzz_zz_1[i] * fe_0 + ta1_z_xxzzz_xzz_0[i] * pa_x[i] - ta1_z_xxzzz_xzz_1[i] * pc_x[i];

        ta1_z_xxxzzz_yyy_0[i] =
            2.0 * ta1_z_xzzz_yyy_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yyy_1[i] * fe_0 + ta1_z_xxzzz_yyy_0[i] * pa_x[i] - ta1_z_xxzzz_yyy_1[i] * pc_x[i];

        ta1_z_xxxzzz_yyz_0[i] =
            2.0 * ta1_z_xzzz_yyz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yyz_1[i] * fe_0 + ta1_z_xxzzz_yyz_0[i] * pa_x[i] - ta1_z_xxzzz_yyz_1[i] * pc_x[i];

        ta1_z_xxxzzz_yzz_0[i] =
            2.0 * ta1_z_xzzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yzz_1[i] * fe_0 + ta1_z_xxzzz_yzz_0[i] * pa_x[i] - ta1_z_xxzzz_yzz_1[i] * pc_x[i];

        ta1_z_xxxzzz_zzz_0[i] =
            2.0 * ta1_z_xzzz_zzz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_zzz_1[i] * fe_0 + ta1_z_xxzzz_zzz_0[i] * pa_x[i] - ta1_z_xxzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 660-670 components of targeted buffer : IF

    auto ta1_z_xxyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 660);

    auto ta1_z_xxyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 661);

    auto ta1_z_xxyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 662);

    auto ta1_z_xxyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 663);

    auto ta1_z_xxyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 664);

    auto ta1_z_xxyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 665);

    auto ta1_z_xxyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 666);

    auto ta1_z_xxyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 667);

    auto ta1_z_xxyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 668);

    auto ta1_z_xxyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 669);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxyy_xxx_0,   \
                             ta1_z_xxyy_xxx_1,   \
                             ta1_z_xxyy_xxz_0,   \
                             ta1_z_xxyy_xxz_1,   \
                             ta1_z_xxyy_xzz_0,   \
                             ta1_z_xxyy_xzz_1,   \
                             ta1_z_xxyyy_xxx_0,  \
                             ta1_z_xxyyy_xxx_1,  \
                             ta1_z_xxyyy_xxz_0,  \
                             ta1_z_xxyyy_xxz_1,  \
                             ta1_z_xxyyy_xzz_0,  \
                             ta1_z_xxyyy_xzz_1,  \
                             ta1_z_xxyyyy_xxx_0, \
                             ta1_z_xxyyyy_xxy_0, \
                             ta1_z_xxyyyy_xxz_0, \
                             ta1_z_xxyyyy_xyy_0, \
                             ta1_z_xxyyyy_xyz_0, \
                             ta1_z_xxyyyy_xzz_0, \
                             ta1_z_xxyyyy_yyy_0, \
                             ta1_z_xxyyyy_yyz_0, \
                             ta1_z_xxyyyy_yzz_0, \
                             ta1_z_xxyyyy_zzz_0, \
                             ta1_z_xyyyy_xxy_0,  \
                             ta1_z_xyyyy_xxy_1,  \
                             ta1_z_xyyyy_xy_0,   \
                             ta1_z_xyyyy_xy_1,   \
                             ta1_z_xyyyy_xyy_0,  \
                             ta1_z_xyyyy_xyy_1,  \
                             ta1_z_xyyyy_xyz_0,  \
                             ta1_z_xyyyy_xyz_1,  \
                             ta1_z_xyyyy_yy_0,   \
                             ta1_z_xyyyy_yy_1,   \
                             ta1_z_xyyyy_yyy_0,  \
                             ta1_z_xyyyy_yyy_1,  \
                             ta1_z_xyyyy_yyz_0,  \
                             ta1_z_xyyyy_yyz_1,  \
                             ta1_z_xyyyy_yz_0,   \
                             ta1_z_xyyyy_yz_1,   \
                             ta1_z_xyyyy_yzz_0,  \
                             ta1_z_xyyyy_yzz_1,  \
                             ta1_z_xyyyy_zzz_0,  \
                             ta1_z_xyyyy_zzz_1,  \
                             ta1_z_yyyy_xxy_0,   \
                             ta1_z_yyyy_xxy_1,   \
                             ta1_z_yyyy_xyy_0,   \
                             ta1_z_yyyy_xyy_1,   \
                             ta1_z_yyyy_xyz_0,   \
                             ta1_z_yyyy_xyz_1,   \
                             ta1_z_yyyy_yyy_0,   \
                             ta1_z_yyyy_yyy_1,   \
                             ta1_z_yyyy_yyz_0,   \
                             ta1_z_yyyy_yyz_1,   \
                             ta1_z_yyyy_yzz_0,   \
                             ta1_z_yyyy_yzz_1,   \
                             ta1_z_yyyy_zzz_0,   \
                             ta1_z_yyyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyy_xxx_0[i] =
            3.0 * ta1_z_xxyy_xxx_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxx_1[i] * fe_0 + ta1_z_xxyyy_xxx_0[i] * pa_y[i] - ta1_z_xxyyy_xxx_1[i] * pc_y[i];

        ta1_z_xxyyyy_xxy_0[i] = ta1_z_yyyy_xxy_0[i] * fe_0 - ta1_z_yyyy_xxy_1[i] * fe_0 + 2.0 * ta1_z_xyyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_z_xyyyy_xy_1[i] * fe_0 + ta1_z_xyyyy_xxy_0[i] * pa_x[i] - ta1_z_xyyyy_xxy_1[i] * pc_x[i];

        ta1_z_xxyyyy_xxz_0[i] =
            3.0 * ta1_z_xxyy_xxz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xxz_1[i] * fe_0 + ta1_z_xxyyy_xxz_0[i] * pa_y[i] - ta1_z_xxyyy_xxz_1[i] * pc_y[i];

        ta1_z_xxyyyy_xyy_0[i] = ta1_z_yyyy_xyy_0[i] * fe_0 - ta1_z_yyyy_xyy_1[i] * fe_0 + ta1_z_xyyyy_yy_0[i] * fe_0 - ta1_z_xyyyy_yy_1[i] * fe_0 +
                                ta1_z_xyyyy_xyy_0[i] * pa_x[i] - ta1_z_xyyyy_xyy_1[i] * pc_x[i];

        ta1_z_xxyyyy_xyz_0[i] = ta1_z_yyyy_xyz_0[i] * fe_0 - ta1_z_yyyy_xyz_1[i] * fe_0 + ta1_z_xyyyy_yz_0[i] * fe_0 - ta1_z_xyyyy_yz_1[i] * fe_0 +
                                ta1_z_xyyyy_xyz_0[i] * pa_x[i] - ta1_z_xyyyy_xyz_1[i] * pc_x[i];

        ta1_z_xxyyyy_xzz_0[i] =
            3.0 * ta1_z_xxyy_xzz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xzz_1[i] * fe_0 + ta1_z_xxyyy_xzz_0[i] * pa_y[i] - ta1_z_xxyyy_xzz_1[i] * pc_y[i];

        ta1_z_xxyyyy_yyy_0[i] =
            ta1_z_yyyy_yyy_0[i] * fe_0 - ta1_z_yyyy_yyy_1[i] * fe_0 + ta1_z_xyyyy_yyy_0[i] * pa_x[i] - ta1_z_xyyyy_yyy_1[i] * pc_x[i];

        ta1_z_xxyyyy_yyz_0[i] =
            ta1_z_yyyy_yyz_0[i] * fe_0 - ta1_z_yyyy_yyz_1[i] * fe_0 + ta1_z_xyyyy_yyz_0[i] * pa_x[i] - ta1_z_xyyyy_yyz_1[i] * pc_x[i];

        ta1_z_xxyyyy_yzz_0[i] =
            ta1_z_yyyy_yzz_0[i] * fe_0 - ta1_z_yyyy_yzz_1[i] * fe_0 + ta1_z_xyyyy_yzz_0[i] * pa_x[i] - ta1_z_xyyyy_yzz_1[i] * pc_x[i];

        ta1_z_xxyyyy_zzz_0[i] =
            ta1_z_yyyy_zzz_0[i] * fe_0 - ta1_z_yyyy_zzz_1[i] * fe_0 + ta1_z_xyyyy_zzz_0[i] * pa_x[i] - ta1_z_xyyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 670-680 components of targeted buffer : IF

    auto ta1_z_xxyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 670);

    auto ta1_z_xxyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 671);

    auto ta1_z_xxyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 672);

    auto ta1_z_xxyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 673);

    auto ta1_z_xxyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 674);

    auto ta1_z_xxyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 675);

    auto ta1_z_xxyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 676);

    auto ta1_z_xxyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 677);

    auto ta1_z_xxyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 678);

    auto ta1_z_xxyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 679);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxyyy_xxx_0,  \
                             ta1_z_xxyyy_xxx_1,  \
                             ta1_z_xxyyy_xxy_0,  \
                             ta1_z_xxyyy_xxy_1,  \
                             ta1_z_xxyyy_xy_0,   \
                             ta1_z_xxyyy_xy_1,   \
                             ta1_z_xxyyy_xyy_0,  \
                             ta1_z_xxyyy_xyy_1,  \
                             ta1_z_xxyyy_xyz_0,  \
                             ta1_z_xxyyy_xyz_1,  \
                             ta1_z_xxyyy_yyy_0,  \
                             ta1_z_xxyyy_yyy_1,  \
                             ta1_z_xxyyyz_xxx_0, \
                             ta1_z_xxyyyz_xxy_0, \
                             ta1_z_xxyyyz_xxz_0, \
                             ta1_z_xxyyyz_xyy_0, \
                             ta1_z_xxyyyz_xyz_0, \
                             ta1_z_xxyyyz_xzz_0, \
                             ta1_z_xxyyyz_yyy_0, \
                             ta1_z_xxyyyz_yyz_0, \
                             ta1_z_xxyyyz_yzz_0, \
                             ta1_z_xxyyyz_zzz_0, \
                             ta1_z_xxyyz_xxz_0,  \
                             ta1_z_xxyyz_xxz_1,  \
                             ta1_z_xxyyz_xzz_0,  \
                             ta1_z_xxyyz_xzz_1,  \
                             ta1_z_xxyz_xxz_0,   \
                             ta1_z_xxyz_xxz_1,   \
                             ta1_z_xxyz_xzz_0,   \
                             ta1_z_xxyz_xzz_1,   \
                             ta1_z_xyyyz_yyz_0,  \
                             ta1_z_xyyyz_yyz_1,  \
                             ta1_z_xyyyz_yzz_0,  \
                             ta1_z_xyyyz_yzz_1,  \
                             ta1_z_xyyyz_zzz_0,  \
                             ta1_z_xyyyz_zzz_1,  \
                             ta1_z_yyyz_yyz_0,   \
                             ta1_z_yyyz_yyz_1,   \
                             ta1_z_yyyz_yzz_0,   \
                             ta1_z_yyyz_yzz_1,   \
                             ta1_z_yyyz_zzz_0,   \
                             ta1_z_yyyz_zzz_1,   \
                             ta_xxyyy_xxx_1,     \
                             ta_xxyyy_xxy_1,     \
                             ta_xxyyy_xyy_1,     \
                             ta_xxyyy_xyz_1,     \
                             ta_xxyyy_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyz_xxx_0[i] = ta_xxyyy_xxx_1[i] + ta1_z_xxyyy_xxx_0[i] * pa_z[i] - ta1_z_xxyyy_xxx_1[i] * pc_z[i];

        ta1_z_xxyyyz_xxy_0[i] = ta_xxyyy_xxy_1[i] + ta1_z_xxyyy_xxy_0[i] * pa_z[i] - ta1_z_xxyyy_xxy_1[i] * pc_z[i];

        ta1_z_xxyyyz_xxz_0[i] =
            2.0 * ta1_z_xxyz_xxz_0[i] * fe_0 - 2.0 * ta1_z_xxyz_xxz_1[i] * fe_0 + ta1_z_xxyyz_xxz_0[i] * pa_y[i] - ta1_z_xxyyz_xxz_1[i] * pc_y[i];

        ta1_z_xxyyyz_xyy_0[i] = ta_xxyyy_xyy_1[i] + ta1_z_xxyyy_xyy_0[i] * pa_z[i] - ta1_z_xxyyy_xyy_1[i] * pc_z[i];

        ta1_z_xxyyyz_xyz_0[i] = ta1_z_xxyyy_xy_0[i] * fe_0 - ta1_z_xxyyy_xy_1[i] * fe_0 + ta_xxyyy_xyz_1[i] + ta1_z_xxyyy_xyz_0[i] * pa_z[i] -
                                ta1_z_xxyyy_xyz_1[i] * pc_z[i];

        ta1_z_xxyyyz_xzz_0[i] =
            2.0 * ta1_z_xxyz_xzz_0[i] * fe_0 - 2.0 * ta1_z_xxyz_xzz_1[i] * fe_0 + ta1_z_xxyyz_xzz_0[i] * pa_y[i] - ta1_z_xxyyz_xzz_1[i] * pc_y[i];

        ta1_z_xxyyyz_yyy_0[i] = ta_xxyyy_yyy_1[i] + ta1_z_xxyyy_yyy_0[i] * pa_z[i] - ta1_z_xxyyy_yyy_1[i] * pc_z[i];

        ta1_z_xxyyyz_yyz_0[i] =
            ta1_z_yyyz_yyz_0[i] * fe_0 - ta1_z_yyyz_yyz_1[i] * fe_0 + ta1_z_xyyyz_yyz_0[i] * pa_x[i] - ta1_z_xyyyz_yyz_1[i] * pc_x[i];

        ta1_z_xxyyyz_yzz_0[i] =
            ta1_z_yyyz_yzz_0[i] * fe_0 - ta1_z_yyyz_yzz_1[i] * fe_0 + ta1_z_xyyyz_yzz_0[i] * pa_x[i] - ta1_z_xyyyz_yzz_1[i] * pc_x[i];

        ta1_z_xxyyyz_zzz_0[i] =
            ta1_z_yyyz_zzz_0[i] * fe_0 - ta1_z_yyyz_zzz_1[i] * fe_0 + ta1_z_xyyyz_zzz_0[i] * pa_x[i] - ta1_z_xyyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 680-690 components of targeted buffer : IF

    auto ta1_z_xxyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 680);

    auto ta1_z_xxyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 681);

    auto ta1_z_xxyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 682);

    auto ta1_z_xxyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 683);

    auto ta1_z_xxyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 684);

    auto ta1_z_xxyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 685);

    auto ta1_z_xxyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 686);

    auto ta1_z_xxyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 687);

    auto ta1_z_xxyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 688);

    auto ta1_z_xxyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 689);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxyy_xxy_0,   \
                             ta1_z_xxyy_xxy_1,   \
                             ta1_z_xxyy_xyy_0,   \
                             ta1_z_xxyy_xyy_1,   \
                             ta1_z_xxyyz_xxy_0,  \
                             ta1_z_xxyyz_xxy_1,  \
                             ta1_z_xxyyz_xyy_0,  \
                             ta1_z_xxyyz_xyy_1,  \
                             ta1_z_xxyyzz_xxx_0, \
                             ta1_z_xxyyzz_xxy_0, \
                             ta1_z_xxyyzz_xxz_0, \
                             ta1_z_xxyyzz_xyy_0, \
                             ta1_z_xxyyzz_xyz_0, \
                             ta1_z_xxyyzz_xzz_0, \
                             ta1_z_xxyyzz_yyy_0, \
                             ta1_z_xxyyzz_yyz_0, \
                             ta1_z_xxyyzz_yzz_0, \
                             ta1_z_xxyyzz_zzz_0, \
                             ta1_z_xxyzz_xxx_0,  \
                             ta1_z_xxyzz_xxx_1,  \
                             ta1_z_xxyzz_xxz_0,  \
                             ta1_z_xxyzz_xxz_1,  \
                             ta1_z_xxyzz_xzz_0,  \
                             ta1_z_xxyzz_xzz_1,  \
                             ta1_z_xxzz_xxx_0,   \
                             ta1_z_xxzz_xxx_1,   \
                             ta1_z_xxzz_xxz_0,   \
                             ta1_z_xxzz_xxz_1,   \
                             ta1_z_xxzz_xzz_0,   \
                             ta1_z_xxzz_xzz_1,   \
                             ta1_z_xyyzz_xyz_0,  \
                             ta1_z_xyyzz_xyz_1,  \
                             ta1_z_xyyzz_yyy_0,  \
                             ta1_z_xyyzz_yyy_1,  \
                             ta1_z_xyyzz_yyz_0,  \
                             ta1_z_xyyzz_yyz_1,  \
                             ta1_z_xyyzz_yz_0,   \
                             ta1_z_xyyzz_yz_1,   \
                             ta1_z_xyyzz_yzz_0,  \
                             ta1_z_xyyzz_yzz_1,  \
                             ta1_z_xyyzz_zzz_0,  \
                             ta1_z_xyyzz_zzz_1,  \
                             ta1_z_yyzz_xyz_0,   \
                             ta1_z_yyzz_xyz_1,   \
                             ta1_z_yyzz_yyy_0,   \
                             ta1_z_yyzz_yyy_1,   \
                             ta1_z_yyzz_yyz_0,   \
                             ta1_z_yyzz_yyz_1,   \
                             ta1_z_yyzz_yzz_0,   \
                             ta1_z_yyzz_yzz_1,   \
                             ta1_z_yyzz_zzz_0,   \
                             ta1_z_yyzz_zzz_1,   \
                             ta_xxyyz_xxy_1,     \
                             ta_xxyyz_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyzz_xxx_0[i] =
            ta1_z_xxzz_xxx_0[i] * fe_0 - ta1_z_xxzz_xxx_1[i] * fe_0 + ta1_z_xxyzz_xxx_0[i] * pa_y[i] - ta1_z_xxyzz_xxx_1[i] * pc_y[i];

        ta1_z_xxyyzz_xxy_0[i] = ta1_z_xxyy_xxy_0[i] * fe_0 - ta1_z_xxyy_xxy_1[i] * fe_0 + ta_xxyyz_xxy_1[i] + ta1_z_xxyyz_xxy_0[i] * pa_z[i] -
                                ta1_z_xxyyz_xxy_1[i] * pc_z[i];

        ta1_z_xxyyzz_xxz_0[i] =
            ta1_z_xxzz_xxz_0[i] * fe_0 - ta1_z_xxzz_xxz_1[i] * fe_0 + ta1_z_xxyzz_xxz_0[i] * pa_y[i] - ta1_z_xxyzz_xxz_1[i] * pc_y[i];

        ta1_z_xxyyzz_xyy_0[i] = ta1_z_xxyy_xyy_0[i] * fe_0 - ta1_z_xxyy_xyy_1[i] * fe_0 + ta_xxyyz_xyy_1[i] + ta1_z_xxyyz_xyy_0[i] * pa_z[i] -
                                ta1_z_xxyyz_xyy_1[i] * pc_z[i];

        ta1_z_xxyyzz_xyz_0[i] = ta1_z_yyzz_xyz_0[i] * fe_0 - ta1_z_yyzz_xyz_1[i] * fe_0 + ta1_z_xyyzz_yz_0[i] * fe_0 - ta1_z_xyyzz_yz_1[i] * fe_0 +
                                ta1_z_xyyzz_xyz_0[i] * pa_x[i] - ta1_z_xyyzz_xyz_1[i] * pc_x[i];

        ta1_z_xxyyzz_xzz_0[i] =
            ta1_z_xxzz_xzz_0[i] * fe_0 - ta1_z_xxzz_xzz_1[i] * fe_0 + ta1_z_xxyzz_xzz_0[i] * pa_y[i] - ta1_z_xxyzz_xzz_1[i] * pc_y[i];

        ta1_z_xxyyzz_yyy_0[i] =
            ta1_z_yyzz_yyy_0[i] * fe_0 - ta1_z_yyzz_yyy_1[i] * fe_0 + ta1_z_xyyzz_yyy_0[i] * pa_x[i] - ta1_z_xyyzz_yyy_1[i] * pc_x[i];

        ta1_z_xxyyzz_yyz_0[i] =
            ta1_z_yyzz_yyz_0[i] * fe_0 - ta1_z_yyzz_yyz_1[i] * fe_0 + ta1_z_xyyzz_yyz_0[i] * pa_x[i] - ta1_z_xyyzz_yyz_1[i] * pc_x[i];

        ta1_z_xxyyzz_yzz_0[i] =
            ta1_z_yyzz_yzz_0[i] * fe_0 - ta1_z_yyzz_yzz_1[i] * fe_0 + ta1_z_xyyzz_yzz_0[i] * pa_x[i] - ta1_z_xyyzz_yzz_1[i] * pc_x[i];

        ta1_z_xxyyzz_zzz_0[i] =
            ta1_z_yyzz_zzz_0[i] * fe_0 - ta1_z_yyzz_zzz_1[i] * fe_0 + ta1_z_xyyzz_zzz_0[i] * pa_x[i] - ta1_z_xyyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 690-700 components of targeted buffer : IF

    auto ta1_z_xxyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 690);

    auto ta1_z_xxyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 691);

    auto ta1_z_xxyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 692);

    auto ta1_z_xxyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 693);

    auto ta1_z_xxyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 694);

    auto ta1_z_xxyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 695);

    auto ta1_z_xxyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 696);

    auto ta1_z_xxyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 697);

    auto ta1_z_xxyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 698);

    auto ta1_z_xxyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 699);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxyzzz_xxx_0, \
                             ta1_z_xxyzzz_xxy_0, \
                             ta1_z_xxyzzz_xxz_0, \
                             ta1_z_xxyzzz_xyy_0, \
                             ta1_z_xxyzzz_xyz_0, \
                             ta1_z_xxyzzz_xzz_0, \
                             ta1_z_xxyzzz_yyy_0, \
                             ta1_z_xxyzzz_yyz_0, \
                             ta1_z_xxyzzz_yzz_0, \
                             ta1_z_xxyzzz_zzz_0, \
                             ta1_z_xxzzz_xx_0,   \
                             ta1_z_xxzzz_xx_1,   \
                             ta1_z_xxzzz_xxx_0,  \
                             ta1_z_xxzzz_xxx_1,  \
                             ta1_z_xxzzz_xxy_0,  \
                             ta1_z_xxzzz_xxy_1,  \
                             ta1_z_xxzzz_xxz_0,  \
                             ta1_z_xxzzz_xxz_1,  \
                             ta1_z_xxzzz_xy_0,   \
                             ta1_z_xxzzz_xy_1,   \
                             ta1_z_xxzzz_xyy_0,  \
                             ta1_z_xxzzz_xyy_1,  \
                             ta1_z_xxzzz_xyz_0,  \
                             ta1_z_xxzzz_xyz_1,  \
                             ta1_z_xxzzz_xz_0,   \
                             ta1_z_xxzzz_xz_1,   \
                             ta1_z_xxzzz_xzz_0,  \
                             ta1_z_xxzzz_xzz_1,  \
                             ta1_z_xxzzz_zzz_0,  \
                             ta1_z_xxzzz_zzz_1,  \
                             ta1_z_xyzzz_yyy_0,  \
                             ta1_z_xyzzz_yyy_1,  \
                             ta1_z_xyzzz_yyz_0,  \
                             ta1_z_xyzzz_yyz_1,  \
                             ta1_z_xyzzz_yzz_0,  \
                             ta1_z_xyzzz_yzz_1,  \
                             ta1_z_yzzz_yyy_0,   \
                             ta1_z_yzzz_yyy_1,   \
                             ta1_z_yzzz_yyz_0,   \
                             ta1_z_yzzz_yyz_1,   \
                             ta1_z_yzzz_yzz_0,   \
                             ta1_z_yzzz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzzz_xxx_0[i] = ta1_z_xxzzz_xxx_0[i] * pa_y[i] - ta1_z_xxzzz_xxx_1[i] * pc_y[i];

        ta1_z_xxyzzz_xxy_0[i] =
            ta1_z_xxzzz_xx_0[i] * fe_0 - ta1_z_xxzzz_xx_1[i] * fe_0 + ta1_z_xxzzz_xxy_0[i] * pa_y[i] - ta1_z_xxzzz_xxy_1[i] * pc_y[i];

        ta1_z_xxyzzz_xxz_0[i] = ta1_z_xxzzz_xxz_0[i] * pa_y[i] - ta1_z_xxzzz_xxz_1[i] * pc_y[i];

        ta1_z_xxyzzz_xyy_0[i] =
            2.0 * ta1_z_xxzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_xxzzz_xy_1[i] * fe_0 + ta1_z_xxzzz_xyy_0[i] * pa_y[i] - ta1_z_xxzzz_xyy_1[i] * pc_y[i];

        ta1_z_xxyzzz_xyz_0[i] =
            ta1_z_xxzzz_xz_0[i] * fe_0 - ta1_z_xxzzz_xz_1[i] * fe_0 + ta1_z_xxzzz_xyz_0[i] * pa_y[i] - ta1_z_xxzzz_xyz_1[i] * pc_y[i];

        ta1_z_xxyzzz_xzz_0[i] = ta1_z_xxzzz_xzz_0[i] * pa_y[i] - ta1_z_xxzzz_xzz_1[i] * pc_y[i];

        ta1_z_xxyzzz_yyy_0[i] =
            ta1_z_yzzz_yyy_0[i] * fe_0 - ta1_z_yzzz_yyy_1[i] * fe_0 + ta1_z_xyzzz_yyy_0[i] * pa_x[i] - ta1_z_xyzzz_yyy_1[i] * pc_x[i];

        ta1_z_xxyzzz_yyz_0[i] =
            ta1_z_yzzz_yyz_0[i] * fe_0 - ta1_z_yzzz_yyz_1[i] * fe_0 + ta1_z_xyzzz_yyz_0[i] * pa_x[i] - ta1_z_xyzzz_yyz_1[i] * pc_x[i];

        ta1_z_xxyzzz_yzz_0[i] =
            ta1_z_yzzz_yzz_0[i] * fe_0 - ta1_z_yzzz_yzz_1[i] * fe_0 + ta1_z_xyzzz_yzz_0[i] * pa_x[i] - ta1_z_xyzzz_yzz_1[i] * pc_x[i];

        ta1_z_xxyzzz_zzz_0[i] = ta1_z_xxzzz_zzz_0[i] * pa_y[i] - ta1_z_xxzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 700-710 components of targeted buffer : IF

    auto ta1_z_xxzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 700);

    auto ta1_z_xxzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 701);

    auto ta1_z_xxzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 702);

    auto ta1_z_xxzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 703);

    auto ta1_z_xxzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 704);

    auto ta1_z_xxzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 705);

    auto ta1_z_xxzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 706);

    auto ta1_z_xxzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 707);

    auto ta1_z_xxzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 708);

    auto ta1_z_xxzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 709);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxzz_xxx_0,   \
                             ta1_z_xxzz_xxx_1,   \
                             ta1_z_xxzz_xxy_0,   \
                             ta1_z_xxzz_xxy_1,   \
                             ta1_z_xxzz_xyy_0,   \
                             ta1_z_xxzz_xyy_1,   \
                             ta1_z_xxzzz_xxx_0,  \
                             ta1_z_xxzzz_xxx_1,  \
                             ta1_z_xxzzz_xxy_0,  \
                             ta1_z_xxzzz_xxy_1,  \
                             ta1_z_xxzzz_xyy_0,  \
                             ta1_z_xxzzz_xyy_1,  \
                             ta1_z_xxzzzz_xxx_0, \
                             ta1_z_xxzzzz_xxy_0, \
                             ta1_z_xxzzzz_xxz_0, \
                             ta1_z_xxzzzz_xyy_0, \
                             ta1_z_xxzzzz_xyz_0, \
                             ta1_z_xxzzzz_xzz_0, \
                             ta1_z_xxzzzz_yyy_0, \
                             ta1_z_xxzzzz_yyz_0, \
                             ta1_z_xxzzzz_yzz_0, \
                             ta1_z_xxzzzz_zzz_0, \
                             ta1_z_xzzzz_xxz_0,  \
                             ta1_z_xzzzz_xxz_1,  \
                             ta1_z_xzzzz_xyz_0,  \
                             ta1_z_xzzzz_xyz_1,  \
                             ta1_z_xzzzz_xz_0,   \
                             ta1_z_xzzzz_xz_1,   \
                             ta1_z_xzzzz_xzz_0,  \
                             ta1_z_xzzzz_xzz_1,  \
                             ta1_z_xzzzz_yyy_0,  \
                             ta1_z_xzzzz_yyy_1,  \
                             ta1_z_xzzzz_yyz_0,  \
                             ta1_z_xzzzz_yyz_1,  \
                             ta1_z_xzzzz_yz_0,   \
                             ta1_z_xzzzz_yz_1,   \
                             ta1_z_xzzzz_yzz_0,  \
                             ta1_z_xzzzz_yzz_1,  \
                             ta1_z_xzzzz_zz_0,   \
                             ta1_z_xzzzz_zz_1,   \
                             ta1_z_xzzzz_zzz_0,  \
                             ta1_z_xzzzz_zzz_1,  \
                             ta1_z_zzzz_xxz_0,   \
                             ta1_z_zzzz_xxz_1,   \
                             ta1_z_zzzz_xyz_0,   \
                             ta1_z_zzzz_xyz_1,   \
                             ta1_z_zzzz_xzz_0,   \
                             ta1_z_zzzz_xzz_1,   \
                             ta1_z_zzzz_yyy_0,   \
                             ta1_z_zzzz_yyy_1,   \
                             ta1_z_zzzz_yyz_0,   \
                             ta1_z_zzzz_yyz_1,   \
                             ta1_z_zzzz_yzz_0,   \
                             ta1_z_zzzz_yzz_1,   \
                             ta1_z_zzzz_zzz_0,   \
                             ta1_z_zzzz_zzz_1,   \
                             ta_xxzzz_xxx_1,     \
                             ta_xxzzz_xxy_1,     \
                             ta_xxzzz_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzzz_xxx_0[i] = 3.0 * ta1_z_xxzz_xxx_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxx_1[i] * fe_0 + ta_xxzzz_xxx_1[i] +
                                ta1_z_xxzzz_xxx_0[i] * pa_z[i] - ta1_z_xxzzz_xxx_1[i] * pc_z[i];

        ta1_z_xxzzzz_xxy_0[i] = 3.0 * ta1_z_xxzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xxy_1[i] * fe_0 + ta_xxzzz_xxy_1[i] +
                                ta1_z_xxzzz_xxy_0[i] * pa_z[i] - ta1_z_xxzzz_xxy_1[i] * pc_z[i];

        ta1_z_xxzzzz_xxz_0[i] = ta1_z_zzzz_xxz_0[i] * fe_0 - ta1_z_zzzz_xxz_1[i] * fe_0 + 2.0 * ta1_z_xzzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_z_xzzzz_xz_1[i] * fe_0 + ta1_z_xzzzz_xxz_0[i] * pa_x[i] - ta1_z_xzzzz_xxz_1[i] * pc_x[i];

        ta1_z_xxzzzz_xyy_0[i] = 3.0 * ta1_z_xxzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xyy_1[i] * fe_0 + ta_xxzzz_xyy_1[i] +
                                ta1_z_xxzzz_xyy_0[i] * pa_z[i] - ta1_z_xxzzz_xyy_1[i] * pc_z[i];

        ta1_z_xxzzzz_xyz_0[i] = ta1_z_zzzz_xyz_0[i] * fe_0 - ta1_z_zzzz_xyz_1[i] * fe_0 + ta1_z_xzzzz_yz_0[i] * fe_0 - ta1_z_xzzzz_yz_1[i] * fe_0 +
                                ta1_z_xzzzz_xyz_0[i] * pa_x[i] - ta1_z_xzzzz_xyz_1[i] * pc_x[i];

        ta1_z_xxzzzz_xzz_0[i] = ta1_z_zzzz_xzz_0[i] * fe_0 - ta1_z_zzzz_xzz_1[i] * fe_0 + ta1_z_xzzzz_zz_0[i] * fe_0 - ta1_z_xzzzz_zz_1[i] * fe_0 +
                                ta1_z_xzzzz_xzz_0[i] * pa_x[i] - ta1_z_xzzzz_xzz_1[i] * pc_x[i];

        ta1_z_xxzzzz_yyy_0[i] =
            ta1_z_zzzz_yyy_0[i] * fe_0 - ta1_z_zzzz_yyy_1[i] * fe_0 + ta1_z_xzzzz_yyy_0[i] * pa_x[i] - ta1_z_xzzzz_yyy_1[i] * pc_x[i];

        ta1_z_xxzzzz_yyz_0[i] =
            ta1_z_zzzz_yyz_0[i] * fe_0 - ta1_z_zzzz_yyz_1[i] * fe_0 + ta1_z_xzzzz_yyz_0[i] * pa_x[i] - ta1_z_xzzzz_yyz_1[i] * pc_x[i];

        ta1_z_xxzzzz_yzz_0[i] =
            ta1_z_zzzz_yzz_0[i] * fe_0 - ta1_z_zzzz_yzz_1[i] * fe_0 + ta1_z_xzzzz_yzz_0[i] * pa_x[i] - ta1_z_xzzzz_yzz_1[i] * pc_x[i];

        ta1_z_xxzzzz_zzz_0[i] =
            ta1_z_zzzz_zzz_0[i] * fe_0 - ta1_z_zzzz_zzz_1[i] * fe_0 + ta1_z_xzzzz_zzz_0[i] * pa_x[i] - ta1_z_xzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 710-720 components of targeted buffer : IF

    auto ta1_z_xyyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 710);

    auto ta1_z_xyyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 711);

    auto ta1_z_xyyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 712);

    auto ta1_z_xyyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 713);

    auto ta1_z_xyyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 714);

    auto ta1_z_xyyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 715);

    auto ta1_z_xyyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 716);

    auto ta1_z_xyyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 717);

    auto ta1_z_xyyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 718);

    auto ta1_z_xyyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 719);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xyyyyy_xxx_0, \
                             ta1_z_xyyyyy_xxy_0, \
                             ta1_z_xyyyyy_xxz_0, \
                             ta1_z_xyyyyy_xyy_0, \
                             ta1_z_xyyyyy_xyz_0, \
                             ta1_z_xyyyyy_xzz_0, \
                             ta1_z_xyyyyy_yyy_0, \
                             ta1_z_xyyyyy_yyz_0, \
                             ta1_z_xyyyyy_yzz_0, \
                             ta1_z_xyyyyy_zzz_0, \
                             ta1_z_yyyyy_xx_0,   \
                             ta1_z_yyyyy_xx_1,   \
                             ta1_z_yyyyy_xxx_0,  \
                             ta1_z_yyyyy_xxx_1,  \
                             ta1_z_yyyyy_xxy_0,  \
                             ta1_z_yyyyy_xxy_1,  \
                             ta1_z_yyyyy_xxz_0,  \
                             ta1_z_yyyyy_xxz_1,  \
                             ta1_z_yyyyy_xy_0,   \
                             ta1_z_yyyyy_xy_1,   \
                             ta1_z_yyyyy_xyy_0,  \
                             ta1_z_yyyyy_xyy_1,  \
                             ta1_z_yyyyy_xyz_0,  \
                             ta1_z_yyyyy_xyz_1,  \
                             ta1_z_yyyyy_xz_0,   \
                             ta1_z_yyyyy_xz_1,   \
                             ta1_z_yyyyy_xzz_0,  \
                             ta1_z_yyyyy_xzz_1,  \
                             ta1_z_yyyyy_yy_0,   \
                             ta1_z_yyyyy_yy_1,   \
                             ta1_z_yyyyy_yyy_0,  \
                             ta1_z_yyyyy_yyy_1,  \
                             ta1_z_yyyyy_yyz_0,  \
                             ta1_z_yyyyy_yyz_1,  \
                             ta1_z_yyyyy_yz_0,   \
                             ta1_z_yyyyy_yz_1,   \
                             ta1_z_yyyyy_yzz_0,  \
                             ta1_z_yyyyy_yzz_1,  \
                             ta1_z_yyyyy_zz_0,   \
                             ta1_z_yyyyy_zz_1,   \
                             ta1_z_yyyyy_zzz_0,  \
                             ta1_z_yyyyy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyyy_xxx_0[i] =
            3.0 * ta1_z_yyyyy_xx_0[i] * fe_0 - 3.0 * ta1_z_yyyyy_xx_1[i] * fe_0 + ta1_z_yyyyy_xxx_0[i] * pa_x[i] - ta1_z_yyyyy_xxx_1[i] * pc_x[i];

        ta1_z_xyyyyy_xxy_0[i] =
            2.0 * ta1_z_yyyyy_xy_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xy_1[i] * fe_0 + ta1_z_yyyyy_xxy_0[i] * pa_x[i] - ta1_z_yyyyy_xxy_1[i] * pc_x[i];

        ta1_z_xyyyyy_xxz_0[i] =
            2.0 * ta1_z_yyyyy_xz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_xz_1[i] * fe_0 + ta1_z_yyyyy_xxz_0[i] * pa_x[i] - ta1_z_yyyyy_xxz_1[i] * pc_x[i];

        ta1_z_xyyyyy_xyy_0[i] =
            ta1_z_yyyyy_yy_0[i] * fe_0 - ta1_z_yyyyy_yy_1[i] * fe_0 + ta1_z_yyyyy_xyy_0[i] * pa_x[i] - ta1_z_yyyyy_xyy_1[i] * pc_x[i];

        ta1_z_xyyyyy_xyz_0[i] =
            ta1_z_yyyyy_yz_0[i] * fe_0 - ta1_z_yyyyy_yz_1[i] * fe_0 + ta1_z_yyyyy_xyz_0[i] * pa_x[i] - ta1_z_yyyyy_xyz_1[i] * pc_x[i];

        ta1_z_xyyyyy_xzz_0[i] =
            ta1_z_yyyyy_zz_0[i] * fe_0 - ta1_z_yyyyy_zz_1[i] * fe_0 + ta1_z_yyyyy_xzz_0[i] * pa_x[i] - ta1_z_yyyyy_xzz_1[i] * pc_x[i];

        ta1_z_xyyyyy_yyy_0[i] = ta1_z_yyyyy_yyy_0[i] * pa_x[i] - ta1_z_yyyyy_yyy_1[i] * pc_x[i];

        ta1_z_xyyyyy_yyz_0[i] = ta1_z_yyyyy_yyz_0[i] * pa_x[i] - ta1_z_yyyyy_yyz_1[i] * pc_x[i];

        ta1_z_xyyyyy_yzz_0[i] = ta1_z_yyyyy_yzz_0[i] * pa_x[i] - ta1_z_yyyyy_yzz_1[i] * pc_x[i];

        ta1_z_xyyyyy_zzz_0[i] = ta1_z_yyyyy_zzz_0[i] * pa_x[i] - ta1_z_yyyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 720-730 components of targeted buffer : IF

    auto ta1_z_xyyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 720);

    auto ta1_z_xyyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 721);

    auto ta1_z_xyyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 722);

    auto ta1_z_xyyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 723);

    auto ta1_z_xyyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 724);

    auto ta1_z_xyyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 725);

    auto ta1_z_xyyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 726);

    auto ta1_z_xyyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 727);

    auto ta1_z_xyyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 728);

    auto ta1_z_xyyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 729);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xyyyy_xxx_0,  \
                             ta1_z_xyyyy_xxx_1,  \
                             ta1_z_xyyyy_xxy_0,  \
                             ta1_z_xyyyy_xxy_1,  \
                             ta1_z_xyyyy_xyy_0,  \
                             ta1_z_xyyyy_xyy_1,  \
                             ta1_z_xyyyyz_xxx_0, \
                             ta1_z_xyyyyz_xxy_0, \
                             ta1_z_xyyyyz_xxz_0, \
                             ta1_z_xyyyyz_xyy_0, \
                             ta1_z_xyyyyz_xyz_0, \
                             ta1_z_xyyyyz_xzz_0, \
                             ta1_z_xyyyyz_yyy_0, \
                             ta1_z_xyyyyz_yyz_0, \
                             ta1_z_xyyyyz_yzz_0, \
                             ta1_z_xyyyyz_zzz_0, \
                             ta1_z_yyyyz_xxz_0,  \
                             ta1_z_yyyyz_xxz_1,  \
                             ta1_z_yyyyz_xyz_0,  \
                             ta1_z_yyyyz_xyz_1,  \
                             ta1_z_yyyyz_xz_0,   \
                             ta1_z_yyyyz_xz_1,   \
                             ta1_z_yyyyz_xzz_0,  \
                             ta1_z_yyyyz_xzz_1,  \
                             ta1_z_yyyyz_yyy_0,  \
                             ta1_z_yyyyz_yyy_1,  \
                             ta1_z_yyyyz_yyz_0,  \
                             ta1_z_yyyyz_yyz_1,  \
                             ta1_z_yyyyz_yz_0,   \
                             ta1_z_yyyyz_yz_1,   \
                             ta1_z_yyyyz_yzz_0,  \
                             ta1_z_yyyyz_yzz_1,  \
                             ta1_z_yyyyz_zz_0,   \
                             ta1_z_yyyyz_zz_1,   \
                             ta1_z_yyyyz_zzz_0,  \
                             ta1_z_yyyyz_zzz_1,  \
                             ta_xyyyy_xxx_1,     \
                             ta_xyyyy_xxy_1,     \
                             ta_xyyyy_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyyz_xxx_0[i] = ta_xyyyy_xxx_1[i] + ta1_z_xyyyy_xxx_0[i] * pa_z[i] - ta1_z_xyyyy_xxx_1[i] * pc_z[i];

        ta1_z_xyyyyz_xxy_0[i] = ta_xyyyy_xxy_1[i] + ta1_z_xyyyy_xxy_0[i] * pa_z[i] - ta1_z_xyyyy_xxy_1[i] * pc_z[i];

        ta1_z_xyyyyz_xxz_0[i] =
            2.0 * ta1_z_yyyyz_xz_0[i] * fe_0 - 2.0 * ta1_z_yyyyz_xz_1[i] * fe_0 + ta1_z_yyyyz_xxz_0[i] * pa_x[i] - ta1_z_yyyyz_xxz_1[i] * pc_x[i];

        ta1_z_xyyyyz_xyy_0[i] = ta_xyyyy_xyy_1[i] + ta1_z_xyyyy_xyy_0[i] * pa_z[i] - ta1_z_xyyyy_xyy_1[i] * pc_z[i];

        ta1_z_xyyyyz_xyz_0[i] =
            ta1_z_yyyyz_yz_0[i] * fe_0 - ta1_z_yyyyz_yz_1[i] * fe_0 + ta1_z_yyyyz_xyz_0[i] * pa_x[i] - ta1_z_yyyyz_xyz_1[i] * pc_x[i];

        ta1_z_xyyyyz_xzz_0[i] =
            ta1_z_yyyyz_zz_0[i] * fe_0 - ta1_z_yyyyz_zz_1[i] * fe_0 + ta1_z_yyyyz_xzz_0[i] * pa_x[i] - ta1_z_yyyyz_xzz_1[i] * pc_x[i];

        ta1_z_xyyyyz_yyy_0[i] = ta1_z_yyyyz_yyy_0[i] * pa_x[i] - ta1_z_yyyyz_yyy_1[i] * pc_x[i];

        ta1_z_xyyyyz_yyz_0[i] = ta1_z_yyyyz_yyz_0[i] * pa_x[i] - ta1_z_yyyyz_yyz_1[i] * pc_x[i];

        ta1_z_xyyyyz_yzz_0[i] = ta1_z_yyyyz_yzz_0[i] * pa_x[i] - ta1_z_yyyyz_yzz_1[i] * pc_x[i];

        ta1_z_xyyyyz_zzz_0[i] = ta1_z_yyyyz_zzz_0[i] * pa_x[i] - ta1_z_yyyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 730-740 components of targeted buffer : IF

    auto ta1_z_xyyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 730);

    auto ta1_z_xyyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 731);

    auto ta1_z_xyyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 732);

    auto ta1_z_xyyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 733);

    auto ta1_z_xyyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 734);

    auto ta1_z_xyyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 735);

    auto ta1_z_xyyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 736);

    auto ta1_z_xyyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 737);

    auto ta1_z_xyyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 738);

    auto ta1_z_xyyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 739);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xyyyzz_xxx_0, \
                             ta1_z_xyyyzz_xxy_0, \
                             ta1_z_xyyyzz_xxz_0, \
                             ta1_z_xyyyzz_xyy_0, \
                             ta1_z_xyyyzz_xyz_0, \
                             ta1_z_xyyyzz_xzz_0, \
                             ta1_z_xyyyzz_yyy_0, \
                             ta1_z_xyyyzz_yyz_0, \
                             ta1_z_xyyyzz_yzz_0, \
                             ta1_z_xyyyzz_zzz_0, \
                             ta1_z_yyyzz_xx_0,   \
                             ta1_z_yyyzz_xx_1,   \
                             ta1_z_yyyzz_xxx_0,  \
                             ta1_z_yyyzz_xxx_1,  \
                             ta1_z_yyyzz_xxy_0,  \
                             ta1_z_yyyzz_xxy_1,  \
                             ta1_z_yyyzz_xxz_0,  \
                             ta1_z_yyyzz_xxz_1,  \
                             ta1_z_yyyzz_xy_0,   \
                             ta1_z_yyyzz_xy_1,   \
                             ta1_z_yyyzz_xyy_0,  \
                             ta1_z_yyyzz_xyy_1,  \
                             ta1_z_yyyzz_xyz_0,  \
                             ta1_z_yyyzz_xyz_1,  \
                             ta1_z_yyyzz_xz_0,   \
                             ta1_z_yyyzz_xz_1,   \
                             ta1_z_yyyzz_xzz_0,  \
                             ta1_z_yyyzz_xzz_1,  \
                             ta1_z_yyyzz_yy_0,   \
                             ta1_z_yyyzz_yy_1,   \
                             ta1_z_yyyzz_yyy_0,  \
                             ta1_z_yyyzz_yyy_1,  \
                             ta1_z_yyyzz_yyz_0,  \
                             ta1_z_yyyzz_yyz_1,  \
                             ta1_z_yyyzz_yz_0,   \
                             ta1_z_yyyzz_yz_1,   \
                             ta1_z_yyyzz_yzz_0,  \
                             ta1_z_yyyzz_yzz_1,  \
                             ta1_z_yyyzz_zz_0,   \
                             ta1_z_yyyzz_zz_1,   \
                             ta1_z_yyyzz_zzz_0,  \
                             ta1_z_yyyzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyzz_xxx_0[i] =
            3.0 * ta1_z_yyyzz_xx_0[i] * fe_0 - 3.0 * ta1_z_yyyzz_xx_1[i] * fe_0 + ta1_z_yyyzz_xxx_0[i] * pa_x[i] - ta1_z_yyyzz_xxx_1[i] * pc_x[i];

        ta1_z_xyyyzz_xxy_0[i] =
            2.0 * ta1_z_yyyzz_xy_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_xy_1[i] * fe_0 + ta1_z_yyyzz_xxy_0[i] * pa_x[i] - ta1_z_yyyzz_xxy_1[i] * pc_x[i];

        ta1_z_xyyyzz_xxz_0[i] =
            2.0 * ta1_z_yyyzz_xz_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_xz_1[i] * fe_0 + ta1_z_yyyzz_xxz_0[i] * pa_x[i] - ta1_z_yyyzz_xxz_1[i] * pc_x[i];

        ta1_z_xyyyzz_xyy_0[i] =
            ta1_z_yyyzz_yy_0[i] * fe_0 - ta1_z_yyyzz_yy_1[i] * fe_0 + ta1_z_yyyzz_xyy_0[i] * pa_x[i] - ta1_z_yyyzz_xyy_1[i] * pc_x[i];

        ta1_z_xyyyzz_xyz_0[i] =
            ta1_z_yyyzz_yz_0[i] * fe_0 - ta1_z_yyyzz_yz_1[i] * fe_0 + ta1_z_yyyzz_xyz_0[i] * pa_x[i] - ta1_z_yyyzz_xyz_1[i] * pc_x[i];

        ta1_z_xyyyzz_xzz_0[i] =
            ta1_z_yyyzz_zz_0[i] * fe_0 - ta1_z_yyyzz_zz_1[i] * fe_0 + ta1_z_yyyzz_xzz_0[i] * pa_x[i] - ta1_z_yyyzz_xzz_1[i] * pc_x[i];

        ta1_z_xyyyzz_yyy_0[i] = ta1_z_yyyzz_yyy_0[i] * pa_x[i] - ta1_z_yyyzz_yyy_1[i] * pc_x[i];

        ta1_z_xyyyzz_yyz_0[i] = ta1_z_yyyzz_yyz_0[i] * pa_x[i] - ta1_z_yyyzz_yyz_1[i] * pc_x[i];

        ta1_z_xyyyzz_yzz_0[i] = ta1_z_yyyzz_yzz_0[i] * pa_x[i] - ta1_z_yyyzz_yzz_1[i] * pc_x[i];

        ta1_z_xyyyzz_zzz_0[i] = ta1_z_yyyzz_zzz_0[i] * pa_x[i] - ta1_z_yyyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 740-750 components of targeted buffer : IF

    auto ta1_z_xyyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 740);

    auto ta1_z_xyyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 741);

    auto ta1_z_xyyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 742);

    auto ta1_z_xyyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 743);

    auto ta1_z_xyyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 744);

    auto ta1_z_xyyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 745);

    auto ta1_z_xyyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 746);

    auto ta1_z_xyyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 747);

    auto ta1_z_xyyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 748);

    auto ta1_z_xyyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 749);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xyyzzz_xxx_0, \
                             ta1_z_xyyzzz_xxy_0, \
                             ta1_z_xyyzzz_xxz_0, \
                             ta1_z_xyyzzz_xyy_0, \
                             ta1_z_xyyzzz_xyz_0, \
                             ta1_z_xyyzzz_xzz_0, \
                             ta1_z_xyyzzz_yyy_0, \
                             ta1_z_xyyzzz_yyz_0, \
                             ta1_z_xyyzzz_yzz_0, \
                             ta1_z_xyyzzz_zzz_0, \
                             ta1_z_yyzzz_xx_0,   \
                             ta1_z_yyzzz_xx_1,   \
                             ta1_z_yyzzz_xxx_0,  \
                             ta1_z_yyzzz_xxx_1,  \
                             ta1_z_yyzzz_xxy_0,  \
                             ta1_z_yyzzz_xxy_1,  \
                             ta1_z_yyzzz_xxz_0,  \
                             ta1_z_yyzzz_xxz_1,  \
                             ta1_z_yyzzz_xy_0,   \
                             ta1_z_yyzzz_xy_1,   \
                             ta1_z_yyzzz_xyy_0,  \
                             ta1_z_yyzzz_xyy_1,  \
                             ta1_z_yyzzz_xyz_0,  \
                             ta1_z_yyzzz_xyz_1,  \
                             ta1_z_yyzzz_xz_0,   \
                             ta1_z_yyzzz_xz_1,   \
                             ta1_z_yyzzz_xzz_0,  \
                             ta1_z_yyzzz_xzz_1,  \
                             ta1_z_yyzzz_yy_0,   \
                             ta1_z_yyzzz_yy_1,   \
                             ta1_z_yyzzz_yyy_0,  \
                             ta1_z_yyzzz_yyy_1,  \
                             ta1_z_yyzzz_yyz_0,  \
                             ta1_z_yyzzz_yyz_1,  \
                             ta1_z_yyzzz_yz_0,   \
                             ta1_z_yyzzz_yz_1,   \
                             ta1_z_yyzzz_yzz_0,  \
                             ta1_z_yyzzz_yzz_1,  \
                             ta1_z_yyzzz_zz_0,   \
                             ta1_z_yyzzz_zz_1,   \
                             ta1_z_yyzzz_zzz_0,  \
                             ta1_z_yyzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzzz_xxx_0[i] =
            3.0 * ta1_z_yyzzz_xx_0[i] * fe_0 - 3.0 * ta1_z_yyzzz_xx_1[i] * fe_0 + ta1_z_yyzzz_xxx_0[i] * pa_x[i] - ta1_z_yyzzz_xxx_1[i] * pc_x[i];

        ta1_z_xyyzzz_xxy_0[i] =
            2.0 * ta1_z_yyzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_xy_1[i] * fe_0 + ta1_z_yyzzz_xxy_0[i] * pa_x[i] - ta1_z_yyzzz_xxy_1[i] * pc_x[i];

        ta1_z_xyyzzz_xxz_0[i] =
            2.0 * ta1_z_yyzzz_xz_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_xz_1[i] * fe_0 + ta1_z_yyzzz_xxz_0[i] * pa_x[i] - ta1_z_yyzzz_xxz_1[i] * pc_x[i];

        ta1_z_xyyzzz_xyy_0[i] =
            ta1_z_yyzzz_yy_0[i] * fe_0 - ta1_z_yyzzz_yy_1[i] * fe_0 + ta1_z_yyzzz_xyy_0[i] * pa_x[i] - ta1_z_yyzzz_xyy_1[i] * pc_x[i];

        ta1_z_xyyzzz_xyz_0[i] =
            ta1_z_yyzzz_yz_0[i] * fe_0 - ta1_z_yyzzz_yz_1[i] * fe_0 + ta1_z_yyzzz_xyz_0[i] * pa_x[i] - ta1_z_yyzzz_xyz_1[i] * pc_x[i];

        ta1_z_xyyzzz_xzz_0[i] =
            ta1_z_yyzzz_zz_0[i] * fe_0 - ta1_z_yyzzz_zz_1[i] * fe_0 + ta1_z_yyzzz_xzz_0[i] * pa_x[i] - ta1_z_yyzzz_xzz_1[i] * pc_x[i];

        ta1_z_xyyzzz_yyy_0[i] = ta1_z_yyzzz_yyy_0[i] * pa_x[i] - ta1_z_yyzzz_yyy_1[i] * pc_x[i];

        ta1_z_xyyzzz_yyz_0[i] = ta1_z_yyzzz_yyz_0[i] * pa_x[i] - ta1_z_yyzzz_yyz_1[i] * pc_x[i];

        ta1_z_xyyzzz_yzz_0[i] = ta1_z_yyzzz_yzz_0[i] * pa_x[i] - ta1_z_yyzzz_yzz_1[i] * pc_x[i];

        ta1_z_xyyzzz_zzz_0[i] = ta1_z_yyzzz_zzz_0[i] * pa_x[i] - ta1_z_yyzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 750-760 components of targeted buffer : IF

    auto ta1_z_xyzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 750);

    auto ta1_z_xyzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 751);

    auto ta1_z_xyzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 752);

    auto ta1_z_xyzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 753);

    auto ta1_z_xyzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 754);

    auto ta1_z_xyzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 755);

    auto ta1_z_xyzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 756);

    auto ta1_z_xyzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 757);

    auto ta1_z_xyzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 758);

    auto ta1_z_xyzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 759);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xyzzzz_xxx_0, \
                             ta1_z_xyzzzz_xxy_0, \
                             ta1_z_xyzzzz_xxz_0, \
                             ta1_z_xyzzzz_xyy_0, \
                             ta1_z_xyzzzz_xyz_0, \
                             ta1_z_xyzzzz_xzz_0, \
                             ta1_z_xyzzzz_yyy_0, \
                             ta1_z_xyzzzz_yyz_0, \
                             ta1_z_xyzzzz_yzz_0, \
                             ta1_z_xyzzzz_zzz_0, \
                             ta1_z_xzzzz_xxx_0,  \
                             ta1_z_xzzzz_xxx_1,  \
                             ta1_z_xzzzz_xxz_0,  \
                             ta1_z_xzzzz_xxz_1,  \
                             ta1_z_xzzzz_xzz_0,  \
                             ta1_z_xzzzz_xzz_1,  \
                             ta1_z_yzzzz_xxy_0,  \
                             ta1_z_yzzzz_xxy_1,  \
                             ta1_z_yzzzz_xy_0,   \
                             ta1_z_yzzzz_xy_1,   \
                             ta1_z_yzzzz_xyy_0,  \
                             ta1_z_yzzzz_xyy_1,  \
                             ta1_z_yzzzz_xyz_0,  \
                             ta1_z_yzzzz_xyz_1,  \
                             ta1_z_yzzzz_yy_0,   \
                             ta1_z_yzzzz_yy_1,   \
                             ta1_z_yzzzz_yyy_0,  \
                             ta1_z_yzzzz_yyy_1,  \
                             ta1_z_yzzzz_yyz_0,  \
                             ta1_z_yzzzz_yyz_1,  \
                             ta1_z_yzzzz_yz_0,   \
                             ta1_z_yzzzz_yz_1,   \
                             ta1_z_yzzzz_yzz_0,  \
                             ta1_z_yzzzz_yzz_1,  \
                             ta1_z_yzzzz_zzz_0,  \
                             ta1_z_yzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzzzz_xxx_0[i] = ta1_z_xzzzz_xxx_0[i] * pa_y[i] - ta1_z_xzzzz_xxx_1[i] * pc_y[i];

        ta1_z_xyzzzz_xxy_0[i] =
            2.0 * ta1_z_yzzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_yzzzz_xy_1[i] * fe_0 + ta1_z_yzzzz_xxy_0[i] * pa_x[i] - ta1_z_yzzzz_xxy_1[i] * pc_x[i];

        ta1_z_xyzzzz_xxz_0[i] = ta1_z_xzzzz_xxz_0[i] * pa_y[i] - ta1_z_xzzzz_xxz_1[i] * pc_y[i];

        ta1_z_xyzzzz_xyy_0[i] =
            ta1_z_yzzzz_yy_0[i] * fe_0 - ta1_z_yzzzz_yy_1[i] * fe_0 + ta1_z_yzzzz_xyy_0[i] * pa_x[i] - ta1_z_yzzzz_xyy_1[i] * pc_x[i];

        ta1_z_xyzzzz_xyz_0[i] =
            ta1_z_yzzzz_yz_0[i] * fe_0 - ta1_z_yzzzz_yz_1[i] * fe_0 + ta1_z_yzzzz_xyz_0[i] * pa_x[i] - ta1_z_yzzzz_xyz_1[i] * pc_x[i];

        ta1_z_xyzzzz_xzz_0[i] = ta1_z_xzzzz_xzz_0[i] * pa_y[i] - ta1_z_xzzzz_xzz_1[i] * pc_y[i];

        ta1_z_xyzzzz_yyy_0[i] = ta1_z_yzzzz_yyy_0[i] * pa_x[i] - ta1_z_yzzzz_yyy_1[i] * pc_x[i];

        ta1_z_xyzzzz_yyz_0[i] = ta1_z_yzzzz_yyz_0[i] * pa_x[i] - ta1_z_yzzzz_yyz_1[i] * pc_x[i];

        ta1_z_xyzzzz_yzz_0[i] = ta1_z_yzzzz_yzz_0[i] * pa_x[i] - ta1_z_yzzzz_yzz_1[i] * pc_x[i];

        ta1_z_xyzzzz_zzz_0[i] = ta1_z_yzzzz_zzz_0[i] * pa_x[i] - ta1_z_yzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 760-770 components of targeted buffer : IF

    auto ta1_z_xzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 760);

    auto ta1_z_xzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 761);

    auto ta1_z_xzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 762);

    auto ta1_z_xzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 763);

    auto ta1_z_xzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 764);

    auto ta1_z_xzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 765);

    auto ta1_z_xzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 766);

    auto ta1_z_xzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 767);

    auto ta1_z_xzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 768);

    auto ta1_z_xzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 769);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xzzzzz_xxx_0, \
                             ta1_z_xzzzzz_xxy_0, \
                             ta1_z_xzzzzz_xxz_0, \
                             ta1_z_xzzzzz_xyy_0, \
                             ta1_z_xzzzzz_xyz_0, \
                             ta1_z_xzzzzz_xzz_0, \
                             ta1_z_xzzzzz_yyy_0, \
                             ta1_z_xzzzzz_yyz_0, \
                             ta1_z_xzzzzz_yzz_0, \
                             ta1_z_xzzzzz_zzz_0, \
                             ta1_z_zzzzz_xx_0,   \
                             ta1_z_zzzzz_xx_1,   \
                             ta1_z_zzzzz_xxx_0,  \
                             ta1_z_zzzzz_xxx_1,  \
                             ta1_z_zzzzz_xxy_0,  \
                             ta1_z_zzzzz_xxy_1,  \
                             ta1_z_zzzzz_xxz_0,  \
                             ta1_z_zzzzz_xxz_1,  \
                             ta1_z_zzzzz_xy_0,   \
                             ta1_z_zzzzz_xy_1,   \
                             ta1_z_zzzzz_xyy_0,  \
                             ta1_z_zzzzz_xyy_1,  \
                             ta1_z_zzzzz_xyz_0,  \
                             ta1_z_zzzzz_xyz_1,  \
                             ta1_z_zzzzz_xz_0,   \
                             ta1_z_zzzzz_xz_1,   \
                             ta1_z_zzzzz_xzz_0,  \
                             ta1_z_zzzzz_xzz_1,  \
                             ta1_z_zzzzz_yy_0,   \
                             ta1_z_zzzzz_yy_1,   \
                             ta1_z_zzzzz_yyy_0,  \
                             ta1_z_zzzzz_yyy_1,  \
                             ta1_z_zzzzz_yyz_0,  \
                             ta1_z_zzzzz_yyz_1,  \
                             ta1_z_zzzzz_yz_0,   \
                             ta1_z_zzzzz_yz_1,   \
                             ta1_z_zzzzz_yzz_0,  \
                             ta1_z_zzzzz_yzz_1,  \
                             ta1_z_zzzzz_zz_0,   \
                             ta1_z_zzzzz_zz_1,   \
                             ta1_z_zzzzz_zzz_0,  \
                             ta1_z_zzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzzz_xxx_0[i] =
            3.0 * ta1_z_zzzzz_xx_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_xx_1[i] * fe_0 + ta1_z_zzzzz_xxx_0[i] * pa_x[i] - ta1_z_zzzzz_xxx_1[i] * pc_x[i];

        ta1_z_xzzzzz_xxy_0[i] =
            2.0 * ta1_z_zzzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xy_1[i] * fe_0 + ta1_z_zzzzz_xxy_0[i] * pa_x[i] - ta1_z_zzzzz_xxy_1[i] * pc_x[i];

        ta1_z_xzzzzz_xxz_0[i] =
            2.0 * ta1_z_zzzzz_xz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xz_1[i] * fe_0 + ta1_z_zzzzz_xxz_0[i] * pa_x[i] - ta1_z_zzzzz_xxz_1[i] * pc_x[i];

        ta1_z_xzzzzz_xyy_0[i] =
            ta1_z_zzzzz_yy_0[i] * fe_0 - ta1_z_zzzzz_yy_1[i] * fe_0 + ta1_z_zzzzz_xyy_0[i] * pa_x[i] - ta1_z_zzzzz_xyy_1[i] * pc_x[i];

        ta1_z_xzzzzz_xyz_0[i] =
            ta1_z_zzzzz_yz_0[i] * fe_0 - ta1_z_zzzzz_yz_1[i] * fe_0 + ta1_z_zzzzz_xyz_0[i] * pa_x[i] - ta1_z_zzzzz_xyz_1[i] * pc_x[i];

        ta1_z_xzzzzz_xzz_0[i] =
            ta1_z_zzzzz_zz_0[i] * fe_0 - ta1_z_zzzzz_zz_1[i] * fe_0 + ta1_z_zzzzz_xzz_0[i] * pa_x[i] - ta1_z_zzzzz_xzz_1[i] * pc_x[i];

        ta1_z_xzzzzz_yyy_0[i] = ta1_z_zzzzz_yyy_0[i] * pa_x[i] - ta1_z_zzzzz_yyy_1[i] * pc_x[i];

        ta1_z_xzzzzz_yyz_0[i] = ta1_z_zzzzz_yyz_0[i] * pa_x[i] - ta1_z_zzzzz_yyz_1[i] * pc_x[i];

        ta1_z_xzzzzz_yzz_0[i] = ta1_z_zzzzz_yzz_0[i] * pa_x[i] - ta1_z_zzzzz_yzz_1[i] * pc_x[i];

        ta1_z_xzzzzz_zzz_0[i] = ta1_z_zzzzz_zzz_0[i] * pa_x[i] - ta1_z_zzzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 770-780 components of targeted buffer : IF

    auto ta1_z_yyyyyy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 770);

    auto ta1_z_yyyyyy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 771);

    auto ta1_z_yyyyyy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 772);

    auto ta1_z_yyyyyy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 773);

    auto ta1_z_yyyyyy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 774);

    auto ta1_z_yyyyyy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 775);

    auto ta1_z_yyyyyy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 776);

    auto ta1_z_yyyyyy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 777);

    auto ta1_z_yyyyyy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 778);

    auto ta1_z_yyyyyy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 779);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yyyy_xxx_0,   \
                             ta1_z_yyyy_xxx_1,   \
                             ta1_z_yyyy_xxy_0,   \
                             ta1_z_yyyy_xxy_1,   \
                             ta1_z_yyyy_xxz_0,   \
                             ta1_z_yyyy_xxz_1,   \
                             ta1_z_yyyy_xyy_0,   \
                             ta1_z_yyyy_xyy_1,   \
                             ta1_z_yyyy_xyz_0,   \
                             ta1_z_yyyy_xyz_1,   \
                             ta1_z_yyyy_xzz_0,   \
                             ta1_z_yyyy_xzz_1,   \
                             ta1_z_yyyy_yyy_0,   \
                             ta1_z_yyyy_yyy_1,   \
                             ta1_z_yyyy_yyz_0,   \
                             ta1_z_yyyy_yyz_1,   \
                             ta1_z_yyyy_yzz_0,   \
                             ta1_z_yyyy_yzz_1,   \
                             ta1_z_yyyy_zzz_0,   \
                             ta1_z_yyyy_zzz_1,   \
                             ta1_z_yyyyy_xx_0,   \
                             ta1_z_yyyyy_xx_1,   \
                             ta1_z_yyyyy_xxx_0,  \
                             ta1_z_yyyyy_xxx_1,  \
                             ta1_z_yyyyy_xxy_0,  \
                             ta1_z_yyyyy_xxy_1,  \
                             ta1_z_yyyyy_xxz_0,  \
                             ta1_z_yyyyy_xxz_1,  \
                             ta1_z_yyyyy_xy_0,   \
                             ta1_z_yyyyy_xy_1,   \
                             ta1_z_yyyyy_xyy_0,  \
                             ta1_z_yyyyy_xyy_1,  \
                             ta1_z_yyyyy_xyz_0,  \
                             ta1_z_yyyyy_xyz_1,  \
                             ta1_z_yyyyy_xz_0,   \
                             ta1_z_yyyyy_xz_1,   \
                             ta1_z_yyyyy_xzz_0,  \
                             ta1_z_yyyyy_xzz_1,  \
                             ta1_z_yyyyy_yy_0,   \
                             ta1_z_yyyyy_yy_1,   \
                             ta1_z_yyyyy_yyy_0,  \
                             ta1_z_yyyyy_yyy_1,  \
                             ta1_z_yyyyy_yyz_0,  \
                             ta1_z_yyyyy_yyz_1,  \
                             ta1_z_yyyyy_yz_0,   \
                             ta1_z_yyyyy_yz_1,   \
                             ta1_z_yyyyy_yzz_0,  \
                             ta1_z_yyyyy_yzz_1,  \
                             ta1_z_yyyyy_zz_0,   \
                             ta1_z_yyyyy_zz_1,   \
                             ta1_z_yyyyy_zzz_0,  \
                             ta1_z_yyyyy_zzz_1,  \
                             ta1_z_yyyyyy_xxx_0, \
                             ta1_z_yyyyyy_xxy_0, \
                             ta1_z_yyyyyy_xxz_0, \
                             ta1_z_yyyyyy_xyy_0, \
                             ta1_z_yyyyyy_xyz_0, \
                             ta1_z_yyyyyy_xzz_0, \
                             ta1_z_yyyyyy_yyy_0, \
                             ta1_z_yyyyyy_yyz_0, \
                             ta1_z_yyyyyy_yzz_0, \
                             ta1_z_yyyyyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyy_xxx_0[i] =
            5.0 * ta1_z_yyyy_xxx_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxx_1[i] * fe_0 + ta1_z_yyyyy_xxx_0[i] * pa_y[i] - ta1_z_yyyyy_xxx_1[i] * pc_y[i];

        ta1_z_yyyyyy_xxy_0[i] = 5.0 * ta1_z_yyyy_xxy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxy_1[i] * fe_0 + ta1_z_yyyyy_xx_0[i] * fe_0 -
                                ta1_z_yyyyy_xx_1[i] * fe_0 + ta1_z_yyyyy_xxy_0[i] * pa_y[i] - ta1_z_yyyyy_xxy_1[i] * pc_y[i];

        ta1_z_yyyyyy_xxz_0[i] =
            5.0 * ta1_z_yyyy_xxz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xxz_1[i] * fe_0 + ta1_z_yyyyy_xxz_0[i] * pa_y[i] - ta1_z_yyyyy_xxz_1[i] * pc_y[i];

        ta1_z_yyyyyy_xyy_0[i] = 5.0 * ta1_z_yyyy_xyy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xyy_1[i] * fe_0 + 2.0 * ta1_z_yyyyy_xy_0[i] * fe_0 -
                                2.0 * ta1_z_yyyyy_xy_1[i] * fe_0 + ta1_z_yyyyy_xyy_0[i] * pa_y[i] - ta1_z_yyyyy_xyy_1[i] * pc_y[i];

        ta1_z_yyyyyy_xyz_0[i] = 5.0 * ta1_z_yyyy_xyz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xyz_1[i] * fe_0 + ta1_z_yyyyy_xz_0[i] * fe_0 -
                                ta1_z_yyyyy_xz_1[i] * fe_0 + ta1_z_yyyyy_xyz_0[i] * pa_y[i] - ta1_z_yyyyy_xyz_1[i] * pc_y[i];

        ta1_z_yyyyyy_xzz_0[i] =
            5.0 * ta1_z_yyyy_xzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xzz_1[i] * fe_0 + ta1_z_yyyyy_xzz_0[i] * pa_y[i] - ta1_z_yyyyy_xzz_1[i] * pc_y[i];

        ta1_z_yyyyyy_yyy_0[i] = 5.0 * ta1_z_yyyy_yyy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yyy_1[i] * fe_0 + 3.0 * ta1_z_yyyyy_yy_0[i] * fe_0 -
                                3.0 * ta1_z_yyyyy_yy_1[i] * fe_0 + ta1_z_yyyyy_yyy_0[i] * pa_y[i] - ta1_z_yyyyy_yyy_1[i] * pc_y[i];

        ta1_z_yyyyyy_yyz_0[i] = 5.0 * ta1_z_yyyy_yyz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yyz_1[i] * fe_0 + 2.0 * ta1_z_yyyyy_yz_0[i] * fe_0 -
                                2.0 * ta1_z_yyyyy_yz_1[i] * fe_0 + ta1_z_yyyyy_yyz_0[i] * pa_y[i] - ta1_z_yyyyy_yyz_1[i] * pc_y[i];

        ta1_z_yyyyyy_yzz_0[i] = 5.0 * ta1_z_yyyy_yzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yzz_1[i] * fe_0 + ta1_z_yyyyy_zz_0[i] * fe_0 -
                                ta1_z_yyyyy_zz_1[i] * fe_0 + ta1_z_yyyyy_yzz_0[i] * pa_y[i] - ta1_z_yyyyy_yzz_1[i] * pc_y[i];

        ta1_z_yyyyyy_zzz_0[i] =
            5.0 * ta1_z_yyyy_zzz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_zzz_1[i] * fe_0 + ta1_z_yyyyy_zzz_0[i] * pa_y[i] - ta1_z_yyyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 780-790 components of targeted buffer : IF

    auto ta1_z_yyyyyz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 780);

    auto ta1_z_yyyyyz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 781);

    auto ta1_z_yyyyyz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 782);

    auto ta1_z_yyyyyz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 783);

    auto ta1_z_yyyyyz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 784);

    auto ta1_z_yyyyyz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 785);

    auto ta1_z_yyyyyz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 786);

    auto ta1_z_yyyyyz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 787);

    auto ta1_z_yyyyyz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 788);

    auto ta1_z_yyyyyz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 789);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyyyy_xxx_0,  \
                             ta1_z_yyyyy_xxx_1,  \
                             ta1_z_yyyyy_xxy_0,  \
                             ta1_z_yyyyy_xxy_1,  \
                             ta1_z_yyyyy_xy_0,   \
                             ta1_z_yyyyy_xy_1,   \
                             ta1_z_yyyyy_xyy_0,  \
                             ta1_z_yyyyy_xyy_1,  \
                             ta1_z_yyyyy_xyz_0,  \
                             ta1_z_yyyyy_xyz_1,  \
                             ta1_z_yyyyy_yy_0,   \
                             ta1_z_yyyyy_yy_1,   \
                             ta1_z_yyyyy_yyy_0,  \
                             ta1_z_yyyyy_yyy_1,  \
                             ta1_z_yyyyy_yyz_0,  \
                             ta1_z_yyyyy_yyz_1,  \
                             ta1_z_yyyyy_yz_0,   \
                             ta1_z_yyyyy_yz_1,   \
                             ta1_z_yyyyy_yzz_0,  \
                             ta1_z_yyyyy_yzz_1,  \
                             ta1_z_yyyyyz_xxx_0, \
                             ta1_z_yyyyyz_xxy_0, \
                             ta1_z_yyyyyz_xxz_0, \
                             ta1_z_yyyyyz_xyy_0, \
                             ta1_z_yyyyyz_xyz_0, \
                             ta1_z_yyyyyz_xzz_0, \
                             ta1_z_yyyyyz_yyy_0, \
                             ta1_z_yyyyyz_yyz_0, \
                             ta1_z_yyyyyz_yzz_0, \
                             ta1_z_yyyyyz_zzz_0, \
                             ta1_z_yyyyz_xxz_0,  \
                             ta1_z_yyyyz_xxz_1,  \
                             ta1_z_yyyyz_xzz_0,  \
                             ta1_z_yyyyz_xzz_1,  \
                             ta1_z_yyyyz_zzz_0,  \
                             ta1_z_yyyyz_zzz_1,  \
                             ta1_z_yyyz_xxz_0,   \
                             ta1_z_yyyz_xxz_1,   \
                             ta1_z_yyyz_xzz_0,   \
                             ta1_z_yyyz_xzz_1,   \
                             ta1_z_yyyz_zzz_0,   \
                             ta1_z_yyyz_zzz_1,   \
                             ta_yyyyy_xxx_1,     \
                             ta_yyyyy_xxy_1,     \
                             ta_yyyyy_xyy_1,     \
                             ta_yyyyy_xyz_1,     \
                             ta_yyyyy_yyy_1,     \
                             ta_yyyyy_yyz_1,     \
                             ta_yyyyy_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyz_xxx_0[i] = ta_yyyyy_xxx_1[i] + ta1_z_yyyyy_xxx_0[i] * pa_z[i] - ta1_z_yyyyy_xxx_1[i] * pc_z[i];

        ta1_z_yyyyyz_xxy_0[i] = ta_yyyyy_xxy_1[i] + ta1_z_yyyyy_xxy_0[i] * pa_z[i] - ta1_z_yyyyy_xxy_1[i] * pc_z[i];

        ta1_z_yyyyyz_xxz_0[i] =
            4.0 * ta1_z_yyyz_xxz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_xxz_1[i] * fe_0 + ta1_z_yyyyz_xxz_0[i] * pa_y[i] - ta1_z_yyyyz_xxz_1[i] * pc_y[i];

        ta1_z_yyyyyz_xyy_0[i] = ta_yyyyy_xyy_1[i] + ta1_z_yyyyy_xyy_0[i] * pa_z[i] - ta1_z_yyyyy_xyy_1[i] * pc_z[i];

        ta1_z_yyyyyz_xyz_0[i] = ta1_z_yyyyy_xy_0[i] * fe_0 - ta1_z_yyyyy_xy_1[i] * fe_0 + ta_yyyyy_xyz_1[i] + ta1_z_yyyyy_xyz_0[i] * pa_z[i] -
                                ta1_z_yyyyy_xyz_1[i] * pc_z[i];

        ta1_z_yyyyyz_xzz_0[i] =
            4.0 * ta1_z_yyyz_xzz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_xzz_1[i] * fe_0 + ta1_z_yyyyz_xzz_0[i] * pa_y[i] - ta1_z_yyyyz_xzz_1[i] * pc_y[i];

        ta1_z_yyyyyz_yyy_0[i] = ta_yyyyy_yyy_1[i] + ta1_z_yyyyy_yyy_0[i] * pa_z[i] - ta1_z_yyyyy_yyy_1[i] * pc_z[i];

        ta1_z_yyyyyz_yyz_0[i] = ta1_z_yyyyy_yy_0[i] * fe_0 - ta1_z_yyyyy_yy_1[i] * fe_0 + ta_yyyyy_yyz_1[i] + ta1_z_yyyyy_yyz_0[i] * pa_z[i] -
                                ta1_z_yyyyy_yyz_1[i] * pc_z[i];

        ta1_z_yyyyyz_yzz_0[i] = 2.0 * ta1_z_yyyyy_yz_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_yz_1[i] * fe_0 + ta_yyyyy_yzz_1[i] +
                                ta1_z_yyyyy_yzz_0[i] * pa_z[i] - ta1_z_yyyyy_yzz_1[i] * pc_z[i];

        ta1_z_yyyyyz_zzz_0[i] =
            4.0 * ta1_z_yyyz_zzz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_zzz_1[i] * fe_0 + ta1_z_yyyyz_zzz_0[i] * pa_y[i] - ta1_z_yyyyz_zzz_1[i] * pc_y[i];
    }

    // Set up 790-800 components of targeted buffer : IF

    auto ta1_z_yyyyzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 790);

    auto ta1_z_yyyyzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 791);

    auto ta1_z_yyyyzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 792);

    auto ta1_z_yyyyzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 793);

    auto ta1_z_yyyyzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 794);

    auto ta1_z_yyyyzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 795);

    auto ta1_z_yyyyzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 796);

    auto ta1_z_yyyyzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 797);

    auto ta1_z_yyyyzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 798);

    auto ta1_z_yyyyzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 799);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyyy_xxy_0,   \
                             ta1_z_yyyy_xxy_1,   \
                             ta1_z_yyyy_xyy_0,   \
                             ta1_z_yyyy_xyy_1,   \
                             ta1_z_yyyy_yyy_0,   \
                             ta1_z_yyyy_yyy_1,   \
                             ta1_z_yyyyz_xxy_0,  \
                             ta1_z_yyyyz_xxy_1,  \
                             ta1_z_yyyyz_xyy_0,  \
                             ta1_z_yyyyz_xyy_1,  \
                             ta1_z_yyyyz_yyy_0,  \
                             ta1_z_yyyyz_yyy_1,  \
                             ta1_z_yyyyzz_xxx_0, \
                             ta1_z_yyyyzz_xxy_0, \
                             ta1_z_yyyyzz_xxz_0, \
                             ta1_z_yyyyzz_xyy_0, \
                             ta1_z_yyyyzz_xyz_0, \
                             ta1_z_yyyyzz_xzz_0, \
                             ta1_z_yyyyzz_yyy_0, \
                             ta1_z_yyyyzz_yyz_0, \
                             ta1_z_yyyyzz_yzz_0, \
                             ta1_z_yyyyzz_zzz_0, \
                             ta1_z_yyyzz_xxx_0,  \
                             ta1_z_yyyzz_xxx_1,  \
                             ta1_z_yyyzz_xxz_0,  \
                             ta1_z_yyyzz_xxz_1,  \
                             ta1_z_yyyzz_xyz_0,  \
                             ta1_z_yyyzz_xyz_1,  \
                             ta1_z_yyyzz_xz_0,   \
                             ta1_z_yyyzz_xz_1,   \
                             ta1_z_yyyzz_xzz_0,  \
                             ta1_z_yyyzz_xzz_1,  \
                             ta1_z_yyyzz_yyz_0,  \
                             ta1_z_yyyzz_yyz_1,  \
                             ta1_z_yyyzz_yz_0,   \
                             ta1_z_yyyzz_yz_1,   \
                             ta1_z_yyyzz_yzz_0,  \
                             ta1_z_yyyzz_yzz_1,  \
                             ta1_z_yyyzz_zz_0,   \
                             ta1_z_yyyzz_zz_1,   \
                             ta1_z_yyyzz_zzz_0,  \
                             ta1_z_yyyzz_zzz_1,  \
                             ta1_z_yyzz_xxx_0,   \
                             ta1_z_yyzz_xxx_1,   \
                             ta1_z_yyzz_xxz_0,   \
                             ta1_z_yyzz_xxz_1,   \
                             ta1_z_yyzz_xyz_0,   \
                             ta1_z_yyzz_xyz_1,   \
                             ta1_z_yyzz_xzz_0,   \
                             ta1_z_yyzz_xzz_1,   \
                             ta1_z_yyzz_yyz_0,   \
                             ta1_z_yyzz_yyz_1,   \
                             ta1_z_yyzz_yzz_0,   \
                             ta1_z_yyzz_yzz_1,   \
                             ta1_z_yyzz_zzz_0,   \
                             ta1_z_yyzz_zzz_1,   \
                             ta_yyyyz_xxy_1,     \
                             ta_yyyyz_xyy_1,     \
                             ta_yyyyz_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyzz_xxx_0[i] =
            3.0 * ta1_z_yyzz_xxx_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxx_1[i] * fe_0 + ta1_z_yyyzz_xxx_0[i] * pa_y[i] - ta1_z_yyyzz_xxx_1[i] * pc_y[i];

        ta1_z_yyyyzz_xxy_0[i] = ta1_z_yyyy_xxy_0[i] * fe_0 - ta1_z_yyyy_xxy_1[i] * fe_0 + ta_yyyyz_xxy_1[i] + ta1_z_yyyyz_xxy_0[i] * pa_z[i] -
                                ta1_z_yyyyz_xxy_1[i] * pc_z[i];

        ta1_z_yyyyzz_xxz_0[i] =
            3.0 * ta1_z_yyzz_xxz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxz_1[i] * fe_0 + ta1_z_yyyzz_xxz_0[i] * pa_y[i] - ta1_z_yyyzz_xxz_1[i] * pc_y[i];

        ta1_z_yyyyzz_xyy_0[i] = ta1_z_yyyy_xyy_0[i] * fe_0 - ta1_z_yyyy_xyy_1[i] * fe_0 + ta_yyyyz_xyy_1[i] + ta1_z_yyyyz_xyy_0[i] * pa_z[i] -
                                ta1_z_yyyyz_xyy_1[i] * pc_z[i];

        ta1_z_yyyyzz_xyz_0[i] = 3.0 * ta1_z_yyzz_xyz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xyz_1[i] * fe_0 + ta1_z_yyyzz_xz_0[i] * fe_0 -
                                ta1_z_yyyzz_xz_1[i] * fe_0 + ta1_z_yyyzz_xyz_0[i] * pa_y[i] - ta1_z_yyyzz_xyz_1[i] * pc_y[i];

        ta1_z_yyyyzz_xzz_0[i] =
            3.0 * ta1_z_yyzz_xzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xzz_1[i] * fe_0 + ta1_z_yyyzz_xzz_0[i] * pa_y[i] - ta1_z_yyyzz_xzz_1[i] * pc_y[i];

        ta1_z_yyyyzz_yyy_0[i] = ta1_z_yyyy_yyy_0[i] * fe_0 - ta1_z_yyyy_yyy_1[i] * fe_0 + ta_yyyyz_yyy_1[i] + ta1_z_yyyyz_yyy_0[i] * pa_z[i] -
                                ta1_z_yyyyz_yyy_1[i] * pc_z[i];

        ta1_z_yyyyzz_yyz_0[i] = 3.0 * ta1_z_yyzz_yyz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yyz_1[i] * fe_0 + 2.0 * ta1_z_yyyzz_yz_0[i] * fe_0 -
                                2.0 * ta1_z_yyyzz_yz_1[i] * fe_0 + ta1_z_yyyzz_yyz_0[i] * pa_y[i] - ta1_z_yyyzz_yyz_1[i] * pc_y[i];

        ta1_z_yyyyzz_yzz_0[i] = 3.0 * ta1_z_yyzz_yzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yzz_1[i] * fe_0 + ta1_z_yyyzz_zz_0[i] * fe_0 -
                                ta1_z_yyyzz_zz_1[i] * fe_0 + ta1_z_yyyzz_yzz_0[i] * pa_y[i] - ta1_z_yyyzz_yzz_1[i] * pc_y[i];

        ta1_z_yyyyzz_zzz_0[i] =
            3.0 * ta1_z_yyzz_zzz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_zzz_1[i] * fe_0 + ta1_z_yyyzz_zzz_0[i] * pa_y[i] - ta1_z_yyyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 800-810 components of targeted buffer : IF

    auto ta1_z_yyyzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 800);

    auto ta1_z_yyyzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 801);

    auto ta1_z_yyyzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 802);

    auto ta1_z_yyyzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 803);

    auto ta1_z_yyyzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 804);

    auto ta1_z_yyyzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 805);

    auto ta1_z_yyyzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 806);

    auto ta1_z_yyyzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 807);

    auto ta1_z_yyyzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 808);

    auto ta1_z_yyyzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 809);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyyz_xxy_0,   \
                             ta1_z_yyyz_xxy_1,   \
                             ta1_z_yyyz_xyy_0,   \
                             ta1_z_yyyz_xyy_1,   \
                             ta1_z_yyyz_yyy_0,   \
                             ta1_z_yyyz_yyy_1,   \
                             ta1_z_yyyzz_xxy_0,  \
                             ta1_z_yyyzz_xxy_1,  \
                             ta1_z_yyyzz_xyy_0,  \
                             ta1_z_yyyzz_xyy_1,  \
                             ta1_z_yyyzz_yyy_0,  \
                             ta1_z_yyyzz_yyy_1,  \
                             ta1_z_yyyzzz_xxx_0, \
                             ta1_z_yyyzzz_xxy_0, \
                             ta1_z_yyyzzz_xxz_0, \
                             ta1_z_yyyzzz_xyy_0, \
                             ta1_z_yyyzzz_xyz_0, \
                             ta1_z_yyyzzz_xzz_0, \
                             ta1_z_yyyzzz_yyy_0, \
                             ta1_z_yyyzzz_yyz_0, \
                             ta1_z_yyyzzz_yzz_0, \
                             ta1_z_yyyzzz_zzz_0, \
                             ta1_z_yyzzz_xxx_0,  \
                             ta1_z_yyzzz_xxx_1,  \
                             ta1_z_yyzzz_xxz_0,  \
                             ta1_z_yyzzz_xxz_1,  \
                             ta1_z_yyzzz_xyz_0,  \
                             ta1_z_yyzzz_xyz_1,  \
                             ta1_z_yyzzz_xz_0,   \
                             ta1_z_yyzzz_xz_1,   \
                             ta1_z_yyzzz_xzz_0,  \
                             ta1_z_yyzzz_xzz_1,  \
                             ta1_z_yyzzz_yyz_0,  \
                             ta1_z_yyzzz_yyz_1,  \
                             ta1_z_yyzzz_yz_0,   \
                             ta1_z_yyzzz_yz_1,   \
                             ta1_z_yyzzz_yzz_0,  \
                             ta1_z_yyzzz_yzz_1,  \
                             ta1_z_yyzzz_zz_0,   \
                             ta1_z_yyzzz_zz_1,   \
                             ta1_z_yyzzz_zzz_0,  \
                             ta1_z_yyzzz_zzz_1,  \
                             ta1_z_yzzz_xxx_0,   \
                             ta1_z_yzzz_xxx_1,   \
                             ta1_z_yzzz_xxz_0,   \
                             ta1_z_yzzz_xxz_1,   \
                             ta1_z_yzzz_xyz_0,   \
                             ta1_z_yzzz_xyz_1,   \
                             ta1_z_yzzz_xzz_0,   \
                             ta1_z_yzzz_xzz_1,   \
                             ta1_z_yzzz_yyz_0,   \
                             ta1_z_yzzz_yyz_1,   \
                             ta1_z_yzzz_yzz_0,   \
                             ta1_z_yzzz_yzz_1,   \
                             ta1_z_yzzz_zzz_0,   \
                             ta1_z_yzzz_zzz_1,   \
                             ta_yyyzz_xxy_1,     \
                             ta_yyyzz_xyy_1,     \
                             ta_yyyzz_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzzz_xxx_0[i] =
            2.0 * ta1_z_yzzz_xxx_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xxx_1[i] * fe_0 + ta1_z_yyzzz_xxx_0[i] * pa_y[i] - ta1_z_yyzzz_xxx_1[i] * pc_y[i];

        ta1_z_yyyzzz_xxy_0[i] = 2.0 * ta1_z_yyyz_xxy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xxy_1[i] * fe_0 + ta_yyyzz_xxy_1[i] +
                                ta1_z_yyyzz_xxy_0[i] * pa_z[i] - ta1_z_yyyzz_xxy_1[i] * pc_z[i];

        ta1_z_yyyzzz_xxz_0[i] =
            2.0 * ta1_z_yzzz_xxz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xxz_1[i] * fe_0 + ta1_z_yyzzz_xxz_0[i] * pa_y[i] - ta1_z_yyzzz_xxz_1[i] * pc_y[i];

        ta1_z_yyyzzz_xyy_0[i] = 2.0 * ta1_z_yyyz_xyy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xyy_1[i] * fe_0 + ta_yyyzz_xyy_1[i] +
                                ta1_z_yyyzz_xyy_0[i] * pa_z[i] - ta1_z_yyyzz_xyy_1[i] * pc_z[i];

        ta1_z_yyyzzz_xyz_0[i] = 2.0 * ta1_z_yzzz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xyz_1[i] * fe_0 + ta1_z_yyzzz_xz_0[i] * fe_0 -
                                ta1_z_yyzzz_xz_1[i] * fe_0 + ta1_z_yyzzz_xyz_0[i] * pa_y[i] - ta1_z_yyzzz_xyz_1[i] * pc_y[i];

        ta1_z_yyyzzz_xzz_0[i] =
            2.0 * ta1_z_yzzz_xzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xzz_1[i] * fe_0 + ta1_z_yyzzz_xzz_0[i] * pa_y[i] - ta1_z_yyzzz_xzz_1[i] * pc_y[i];

        ta1_z_yyyzzz_yyy_0[i] = 2.0 * ta1_z_yyyz_yyy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_yyy_1[i] * fe_0 + ta_yyyzz_yyy_1[i] +
                                ta1_z_yyyzz_yyy_0[i] * pa_z[i] - ta1_z_yyyzz_yyy_1[i] * pc_z[i];

        ta1_z_yyyzzz_yyz_0[i] = 2.0 * ta1_z_yzzz_yyz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_yyz_1[i] * fe_0 + 2.0 * ta1_z_yyzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_z_yyzzz_yz_1[i] * fe_0 + ta1_z_yyzzz_yyz_0[i] * pa_y[i] - ta1_z_yyzzz_yyz_1[i] * pc_y[i];

        ta1_z_yyyzzz_yzz_0[i] = 2.0 * ta1_z_yzzz_yzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_yzz_1[i] * fe_0 + ta1_z_yyzzz_zz_0[i] * fe_0 -
                                ta1_z_yyzzz_zz_1[i] * fe_0 + ta1_z_yyzzz_yzz_0[i] * pa_y[i] - ta1_z_yyzzz_yzz_1[i] * pc_y[i];

        ta1_z_yyyzzz_zzz_0[i] =
            2.0 * ta1_z_yzzz_zzz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_zzz_1[i] * fe_0 + ta1_z_yyzzz_zzz_0[i] * pa_y[i] - ta1_z_yyzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 810-820 components of targeted buffer : IF

    auto ta1_z_yyzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 810);

    auto ta1_z_yyzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 811);

    auto ta1_z_yyzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 812);

    auto ta1_z_yyzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 813);

    auto ta1_z_yyzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 814);

    auto ta1_z_yyzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 815);

    auto ta1_z_yyzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 816);

    auto ta1_z_yyzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 817);

    auto ta1_z_yyzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 818);

    auto ta1_z_yyzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 819);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyzz_xxy_0,   \
                             ta1_z_yyzz_xxy_1,   \
                             ta1_z_yyzz_xyy_0,   \
                             ta1_z_yyzz_xyy_1,   \
                             ta1_z_yyzz_yyy_0,   \
                             ta1_z_yyzz_yyy_1,   \
                             ta1_z_yyzzz_xxy_0,  \
                             ta1_z_yyzzz_xxy_1,  \
                             ta1_z_yyzzz_xyy_0,  \
                             ta1_z_yyzzz_xyy_1,  \
                             ta1_z_yyzzz_yyy_0,  \
                             ta1_z_yyzzz_yyy_1,  \
                             ta1_z_yyzzzz_xxx_0, \
                             ta1_z_yyzzzz_xxy_0, \
                             ta1_z_yyzzzz_xxz_0, \
                             ta1_z_yyzzzz_xyy_0, \
                             ta1_z_yyzzzz_xyz_0, \
                             ta1_z_yyzzzz_xzz_0, \
                             ta1_z_yyzzzz_yyy_0, \
                             ta1_z_yyzzzz_yyz_0, \
                             ta1_z_yyzzzz_yzz_0, \
                             ta1_z_yyzzzz_zzz_0, \
                             ta1_z_yzzzz_xxx_0,  \
                             ta1_z_yzzzz_xxx_1,  \
                             ta1_z_yzzzz_xxz_0,  \
                             ta1_z_yzzzz_xxz_1,  \
                             ta1_z_yzzzz_xyz_0,  \
                             ta1_z_yzzzz_xyz_1,  \
                             ta1_z_yzzzz_xz_0,   \
                             ta1_z_yzzzz_xz_1,   \
                             ta1_z_yzzzz_xzz_0,  \
                             ta1_z_yzzzz_xzz_1,  \
                             ta1_z_yzzzz_yyz_0,  \
                             ta1_z_yzzzz_yyz_1,  \
                             ta1_z_yzzzz_yz_0,   \
                             ta1_z_yzzzz_yz_1,   \
                             ta1_z_yzzzz_yzz_0,  \
                             ta1_z_yzzzz_yzz_1,  \
                             ta1_z_yzzzz_zz_0,   \
                             ta1_z_yzzzz_zz_1,   \
                             ta1_z_yzzzz_zzz_0,  \
                             ta1_z_yzzzz_zzz_1,  \
                             ta1_z_zzzz_xxx_0,   \
                             ta1_z_zzzz_xxx_1,   \
                             ta1_z_zzzz_xxz_0,   \
                             ta1_z_zzzz_xxz_1,   \
                             ta1_z_zzzz_xyz_0,   \
                             ta1_z_zzzz_xyz_1,   \
                             ta1_z_zzzz_xzz_0,   \
                             ta1_z_zzzz_xzz_1,   \
                             ta1_z_zzzz_yyz_0,   \
                             ta1_z_zzzz_yyz_1,   \
                             ta1_z_zzzz_yzz_0,   \
                             ta1_z_zzzz_yzz_1,   \
                             ta1_z_zzzz_zzz_0,   \
                             ta1_z_zzzz_zzz_1,   \
                             ta_yyzzz_xxy_1,     \
                             ta_yyzzz_xyy_1,     \
                             ta_yyzzz_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzzz_xxx_0[i] =
            ta1_z_zzzz_xxx_0[i] * fe_0 - ta1_z_zzzz_xxx_1[i] * fe_0 + ta1_z_yzzzz_xxx_0[i] * pa_y[i] - ta1_z_yzzzz_xxx_1[i] * pc_y[i];

        ta1_z_yyzzzz_xxy_0[i] = 3.0 * ta1_z_yyzz_xxy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xxy_1[i] * fe_0 + ta_yyzzz_xxy_1[i] +
                                ta1_z_yyzzz_xxy_0[i] * pa_z[i] - ta1_z_yyzzz_xxy_1[i] * pc_z[i];

        ta1_z_yyzzzz_xxz_0[i] =
            ta1_z_zzzz_xxz_0[i] * fe_0 - ta1_z_zzzz_xxz_1[i] * fe_0 + ta1_z_yzzzz_xxz_0[i] * pa_y[i] - ta1_z_yzzzz_xxz_1[i] * pc_y[i];

        ta1_z_yyzzzz_xyy_0[i] = 3.0 * ta1_z_yyzz_xyy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xyy_1[i] * fe_0 + ta_yyzzz_xyy_1[i] +
                                ta1_z_yyzzz_xyy_0[i] * pa_z[i] - ta1_z_yyzzz_xyy_1[i] * pc_z[i];

        ta1_z_yyzzzz_xyz_0[i] = ta1_z_zzzz_xyz_0[i] * fe_0 - ta1_z_zzzz_xyz_1[i] * fe_0 + ta1_z_yzzzz_xz_0[i] * fe_0 - ta1_z_yzzzz_xz_1[i] * fe_0 +
                                ta1_z_yzzzz_xyz_0[i] * pa_y[i] - ta1_z_yzzzz_xyz_1[i] * pc_y[i];

        ta1_z_yyzzzz_xzz_0[i] =
            ta1_z_zzzz_xzz_0[i] * fe_0 - ta1_z_zzzz_xzz_1[i] * fe_0 + ta1_z_yzzzz_xzz_0[i] * pa_y[i] - ta1_z_yzzzz_xzz_1[i] * pc_y[i];

        ta1_z_yyzzzz_yyy_0[i] = 3.0 * ta1_z_yyzz_yyy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yyy_1[i] * fe_0 + ta_yyzzz_yyy_1[i] +
                                ta1_z_yyzzz_yyy_0[i] * pa_z[i] - ta1_z_yyzzz_yyy_1[i] * pc_z[i];

        ta1_z_yyzzzz_yyz_0[i] = ta1_z_zzzz_yyz_0[i] * fe_0 - ta1_z_zzzz_yyz_1[i] * fe_0 + 2.0 * ta1_z_yzzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_z_yzzzz_yz_1[i] * fe_0 + ta1_z_yzzzz_yyz_0[i] * pa_y[i] - ta1_z_yzzzz_yyz_1[i] * pc_y[i];

        ta1_z_yyzzzz_yzz_0[i] = ta1_z_zzzz_yzz_0[i] * fe_0 - ta1_z_zzzz_yzz_1[i] * fe_0 + ta1_z_yzzzz_zz_0[i] * fe_0 - ta1_z_yzzzz_zz_1[i] * fe_0 +
                                ta1_z_yzzzz_yzz_0[i] * pa_y[i] - ta1_z_yzzzz_yzz_1[i] * pc_y[i];

        ta1_z_yyzzzz_zzz_0[i] =
            ta1_z_zzzz_zzz_0[i] * fe_0 - ta1_z_zzzz_zzz_1[i] * fe_0 + ta1_z_yzzzz_zzz_0[i] * pa_y[i] - ta1_z_yzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 820-830 components of targeted buffer : IF

    auto ta1_z_yzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 820);

    auto ta1_z_yzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 821);

    auto ta1_z_yzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 822);

    auto ta1_z_yzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 823);

    auto ta1_z_yzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 824);

    auto ta1_z_yzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 825);

    auto ta1_z_yzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 826);

    auto ta1_z_yzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 827);

    auto ta1_z_yzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 828);

    auto ta1_z_yzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 829);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yzzzzz_xxx_0, \
                             ta1_z_yzzzzz_xxy_0, \
                             ta1_z_yzzzzz_xxz_0, \
                             ta1_z_yzzzzz_xyy_0, \
                             ta1_z_yzzzzz_xyz_0, \
                             ta1_z_yzzzzz_xzz_0, \
                             ta1_z_yzzzzz_yyy_0, \
                             ta1_z_yzzzzz_yyz_0, \
                             ta1_z_yzzzzz_yzz_0, \
                             ta1_z_yzzzzz_zzz_0, \
                             ta1_z_zzzzz_xx_0,   \
                             ta1_z_zzzzz_xx_1,   \
                             ta1_z_zzzzz_xxx_0,  \
                             ta1_z_zzzzz_xxx_1,  \
                             ta1_z_zzzzz_xxy_0,  \
                             ta1_z_zzzzz_xxy_1,  \
                             ta1_z_zzzzz_xxz_0,  \
                             ta1_z_zzzzz_xxz_1,  \
                             ta1_z_zzzzz_xy_0,   \
                             ta1_z_zzzzz_xy_1,   \
                             ta1_z_zzzzz_xyy_0,  \
                             ta1_z_zzzzz_xyy_1,  \
                             ta1_z_zzzzz_xyz_0,  \
                             ta1_z_zzzzz_xyz_1,  \
                             ta1_z_zzzzz_xz_0,   \
                             ta1_z_zzzzz_xz_1,   \
                             ta1_z_zzzzz_xzz_0,  \
                             ta1_z_zzzzz_xzz_1,  \
                             ta1_z_zzzzz_yy_0,   \
                             ta1_z_zzzzz_yy_1,   \
                             ta1_z_zzzzz_yyy_0,  \
                             ta1_z_zzzzz_yyy_1,  \
                             ta1_z_zzzzz_yyz_0,  \
                             ta1_z_zzzzz_yyz_1,  \
                             ta1_z_zzzzz_yz_0,   \
                             ta1_z_zzzzz_yz_1,   \
                             ta1_z_zzzzz_yzz_0,  \
                             ta1_z_zzzzz_yzz_1,  \
                             ta1_z_zzzzz_zz_0,   \
                             ta1_z_zzzzz_zz_1,   \
                             ta1_z_zzzzz_zzz_0,  \
                             ta1_z_zzzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzzz_xxx_0[i] = ta1_z_zzzzz_xxx_0[i] * pa_y[i] - ta1_z_zzzzz_xxx_1[i] * pc_y[i];

        ta1_z_yzzzzz_xxy_0[i] =
            ta1_z_zzzzz_xx_0[i] * fe_0 - ta1_z_zzzzz_xx_1[i] * fe_0 + ta1_z_zzzzz_xxy_0[i] * pa_y[i] - ta1_z_zzzzz_xxy_1[i] * pc_y[i];

        ta1_z_yzzzzz_xxz_0[i] = ta1_z_zzzzz_xxz_0[i] * pa_y[i] - ta1_z_zzzzz_xxz_1[i] * pc_y[i];

        ta1_z_yzzzzz_xyy_0[i] =
            2.0 * ta1_z_zzzzz_xy_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_xy_1[i] * fe_0 + ta1_z_zzzzz_xyy_0[i] * pa_y[i] - ta1_z_zzzzz_xyy_1[i] * pc_y[i];

        ta1_z_yzzzzz_xyz_0[i] =
            ta1_z_zzzzz_xz_0[i] * fe_0 - ta1_z_zzzzz_xz_1[i] * fe_0 + ta1_z_zzzzz_xyz_0[i] * pa_y[i] - ta1_z_zzzzz_xyz_1[i] * pc_y[i];

        ta1_z_yzzzzz_xzz_0[i] = ta1_z_zzzzz_xzz_0[i] * pa_y[i] - ta1_z_zzzzz_xzz_1[i] * pc_y[i];

        ta1_z_yzzzzz_yyy_0[i] =
            3.0 * ta1_z_zzzzz_yy_0[i] * fe_0 - 3.0 * ta1_z_zzzzz_yy_1[i] * fe_0 + ta1_z_zzzzz_yyy_0[i] * pa_y[i] - ta1_z_zzzzz_yyy_1[i] * pc_y[i];

        ta1_z_yzzzzz_yyz_0[i] =
            2.0 * ta1_z_zzzzz_yz_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_yz_1[i] * fe_0 + ta1_z_zzzzz_yyz_0[i] * pa_y[i] - ta1_z_zzzzz_yyz_1[i] * pc_y[i];

        ta1_z_yzzzzz_yzz_0[i] =
            ta1_z_zzzzz_zz_0[i] * fe_0 - ta1_z_zzzzz_zz_1[i] * fe_0 + ta1_z_zzzzz_yzz_0[i] * pa_y[i] - ta1_z_zzzzz_yzz_1[i] * pc_y[i];

        ta1_z_yzzzzz_zzz_0[i] = ta1_z_zzzzz_zzz_0[i] * pa_y[i] - ta1_z_zzzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 830-840 components of targeted buffer : IF

    auto ta1_z_zzzzzz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_if + 830);

    auto ta1_z_zzzzzz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_if + 831);

    auto ta1_z_zzzzzz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_if + 832);

    auto ta1_z_zzzzzz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 833);

    auto ta1_z_zzzzzz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 834);

    auto ta1_z_zzzzzz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 835);

    auto ta1_z_zzzzzz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_if + 836);

    auto ta1_z_zzzzzz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_if + 837);

    auto ta1_z_zzzzzz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 838);

    auto ta1_z_zzzzzz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_if + 839);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_z_zzzz_xxx_0,   \
                             ta1_z_zzzz_xxx_1,   \
                             ta1_z_zzzz_xxy_0,   \
                             ta1_z_zzzz_xxy_1,   \
                             ta1_z_zzzz_xxz_0,   \
                             ta1_z_zzzz_xxz_1,   \
                             ta1_z_zzzz_xyy_0,   \
                             ta1_z_zzzz_xyy_1,   \
                             ta1_z_zzzz_xyz_0,   \
                             ta1_z_zzzz_xyz_1,   \
                             ta1_z_zzzz_xzz_0,   \
                             ta1_z_zzzz_xzz_1,   \
                             ta1_z_zzzz_yyy_0,   \
                             ta1_z_zzzz_yyy_1,   \
                             ta1_z_zzzz_yyz_0,   \
                             ta1_z_zzzz_yyz_1,   \
                             ta1_z_zzzz_yzz_0,   \
                             ta1_z_zzzz_yzz_1,   \
                             ta1_z_zzzz_zzz_0,   \
                             ta1_z_zzzz_zzz_1,   \
                             ta1_z_zzzzz_xx_0,   \
                             ta1_z_zzzzz_xx_1,   \
                             ta1_z_zzzzz_xxx_0,  \
                             ta1_z_zzzzz_xxx_1,  \
                             ta1_z_zzzzz_xxy_0,  \
                             ta1_z_zzzzz_xxy_1,  \
                             ta1_z_zzzzz_xxz_0,  \
                             ta1_z_zzzzz_xxz_1,  \
                             ta1_z_zzzzz_xy_0,   \
                             ta1_z_zzzzz_xy_1,   \
                             ta1_z_zzzzz_xyy_0,  \
                             ta1_z_zzzzz_xyy_1,  \
                             ta1_z_zzzzz_xyz_0,  \
                             ta1_z_zzzzz_xyz_1,  \
                             ta1_z_zzzzz_xz_0,   \
                             ta1_z_zzzzz_xz_1,   \
                             ta1_z_zzzzz_xzz_0,  \
                             ta1_z_zzzzz_xzz_1,  \
                             ta1_z_zzzzz_yy_0,   \
                             ta1_z_zzzzz_yy_1,   \
                             ta1_z_zzzzz_yyy_0,  \
                             ta1_z_zzzzz_yyy_1,  \
                             ta1_z_zzzzz_yyz_0,  \
                             ta1_z_zzzzz_yyz_1,  \
                             ta1_z_zzzzz_yz_0,   \
                             ta1_z_zzzzz_yz_1,   \
                             ta1_z_zzzzz_yzz_0,  \
                             ta1_z_zzzzz_yzz_1,  \
                             ta1_z_zzzzz_zz_0,   \
                             ta1_z_zzzzz_zz_1,   \
                             ta1_z_zzzzz_zzz_0,  \
                             ta1_z_zzzzz_zzz_1,  \
                             ta1_z_zzzzzz_xxx_0, \
                             ta1_z_zzzzzz_xxy_0, \
                             ta1_z_zzzzzz_xxz_0, \
                             ta1_z_zzzzzz_xyy_0, \
                             ta1_z_zzzzzz_xyz_0, \
                             ta1_z_zzzzzz_xzz_0, \
                             ta1_z_zzzzzz_yyy_0, \
                             ta1_z_zzzzzz_yyz_0, \
                             ta1_z_zzzzzz_yzz_0, \
                             ta1_z_zzzzzz_zzz_0, \
                             ta_zzzzz_xxx_1,     \
                             ta_zzzzz_xxy_1,     \
                             ta_zzzzz_xxz_1,     \
                             ta_zzzzz_xyy_1,     \
                             ta_zzzzz_xyz_1,     \
                             ta_zzzzz_xzz_1,     \
                             ta_zzzzz_yyy_1,     \
                             ta_zzzzz_yyz_1,     \
                             ta_zzzzz_yzz_1,     \
                             ta_zzzzz_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzzz_xxx_0[i] = 5.0 * ta1_z_zzzz_xxx_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxx_1[i] * fe_0 + ta_zzzzz_xxx_1[i] +
                                ta1_z_zzzzz_xxx_0[i] * pa_z[i] - ta1_z_zzzzz_xxx_1[i] * pc_z[i];

        ta1_z_zzzzzz_xxy_0[i] = 5.0 * ta1_z_zzzz_xxy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxy_1[i] * fe_0 + ta_zzzzz_xxy_1[i] +
                                ta1_z_zzzzz_xxy_0[i] * pa_z[i] - ta1_z_zzzzz_xxy_1[i] * pc_z[i];

        ta1_z_zzzzzz_xxz_0[i] = 5.0 * ta1_z_zzzz_xxz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xxz_1[i] * fe_0 + ta1_z_zzzzz_xx_0[i] * fe_0 -
                                ta1_z_zzzzz_xx_1[i] * fe_0 + ta_zzzzz_xxz_1[i] + ta1_z_zzzzz_xxz_0[i] * pa_z[i] - ta1_z_zzzzz_xxz_1[i] * pc_z[i];

        ta1_z_zzzzzz_xyy_0[i] = 5.0 * ta1_z_zzzz_xyy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xyy_1[i] * fe_0 + ta_zzzzz_xyy_1[i] +
                                ta1_z_zzzzz_xyy_0[i] * pa_z[i] - ta1_z_zzzzz_xyy_1[i] * pc_z[i];

        ta1_z_zzzzzz_xyz_0[i] = 5.0 * ta1_z_zzzz_xyz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xyz_1[i] * fe_0 + ta1_z_zzzzz_xy_0[i] * fe_0 -
                                ta1_z_zzzzz_xy_1[i] * fe_0 + ta_zzzzz_xyz_1[i] + ta1_z_zzzzz_xyz_0[i] * pa_z[i] - ta1_z_zzzzz_xyz_1[i] * pc_z[i];

        ta1_z_zzzzzz_xzz_0[i] = 5.0 * ta1_z_zzzz_xzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xzz_1[i] * fe_0 + 2.0 * ta1_z_zzzzz_xz_0[i] * fe_0 -
                                2.0 * ta1_z_zzzzz_xz_1[i] * fe_0 + ta_zzzzz_xzz_1[i] + ta1_z_zzzzz_xzz_0[i] * pa_z[i] -
                                ta1_z_zzzzz_xzz_1[i] * pc_z[i];

        ta1_z_zzzzzz_yyy_0[i] = 5.0 * ta1_z_zzzz_yyy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yyy_1[i] * fe_0 + ta_zzzzz_yyy_1[i] +
                                ta1_z_zzzzz_yyy_0[i] * pa_z[i] - ta1_z_zzzzz_yyy_1[i] * pc_z[i];

        ta1_z_zzzzzz_yyz_0[i] = 5.0 * ta1_z_zzzz_yyz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yyz_1[i] * fe_0 + ta1_z_zzzzz_yy_0[i] * fe_0 -
                                ta1_z_zzzzz_yy_1[i] * fe_0 + ta_zzzzz_yyz_1[i] + ta1_z_zzzzz_yyz_0[i] * pa_z[i] - ta1_z_zzzzz_yyz_1[i] * pc_z[i];

        ta1_z_zzzzzz_yzz_0[i] = 5.0 * ta1_z_zzzz_yzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yzz_1[i] * fe_0 + 2.0 * ta1_z_zzzzz_yz_0[i] * fe_0 -
                                2.0 * ta1_z_zzzzz_yz_1[i] * fe_0 + ta_zzzzz_yzz_1[i] + ta1_z_zzzzz_yzz_0[i] * pa_z[i] -
                                ta1_z_zzzzz_yzz_1[i] * pc_z[i];

        ta1_z_zzzzzz_zzz_0[i] = 5.0 * ta1_z_zzzz_zzz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_zzz_1[i] * fe_0 + 3.0 * ta1_z_zzzzz_zz_0[i] * fe_0 -
                                3.0 * ta1_z_zzzzz_zz_1[i] * fe_0 + ta_zzzzz_zzz_1[i] + ta1_z_zzzzz_zzz_0[i] * pa_z[i] -
                                ta1_z_zzzzz_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
