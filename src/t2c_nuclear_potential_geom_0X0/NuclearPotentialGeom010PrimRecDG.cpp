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

#include "NuclearPotentialGeom010PrimRecDG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_dg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_dg,
                                        const size_t              idx_npot_geom_010_0_sg,
                                        const size_t              idx_npot_geom_010_1_sg,
                                        const size_t              idx_npot_geom_010_0_pf,
                                        const size_t              idx_npot_geom_010_1_pf,
                                        const size_t              idx_npot_1_pg,
                                        const size_t              idx_npot_geom_010_0_pg,
                                        const size_t              idx_npot_geom_010_1_pg,
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

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg);

    auto ta1_x_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 1);

    auto ta1_x_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 2);

    auto ta1_x_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 3);

    auto ta1_x_0_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 4);

    auto ta1_x_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 5);

    auto ta1_x_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 6);

    auto ta1_x_0_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 7);

    auto ta1_x_0_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 8);

    auto ta1_x_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 9);

    auto ta1_x_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 10);

    auto ta1_x_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 11);

    auto ta1_x_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 12);

    auto ta1_x_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 13);

    auto ta1_x_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 14);

    auto ta1_y_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 15);

    auto ta1_y_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 16);

    auto ta1_y_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 17);

    auto ta1_y_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 18);

    auto ta1_y_0_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 19);

    auto ta1_y_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 20);

    auto ta1_y_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 21);

    auto ta1_y_0_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 22);

    auto ta1_y_0_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 23);

    auto ta1_y_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 24);

    auto ta1_y_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 25);

    auto ta1_y_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 26);

    auto ta1_y_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 27);

    auto ta1_y_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 28);

    auto ta1_y_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 29);

    auto ta1_z_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 30);

    auto ta1_z_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 31);

    auto ta1_z_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 32);

    auto ta1_z_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 33);

    auto ta1_z_0_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 34);

    auto ta1_z_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 35);

    auto ta1_z_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 36);

    auto ta1_z_0_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 37);

    auto ta1_z_0_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 38);

    auto ta1_z_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 39);

    auto ta1_z_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 40);

    auto ta1_z_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 41);

    auto ta1_z_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 42);

    auto ta1_z_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 43);

    auto ta1_z_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 44);

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg);

    auto ta1_x_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 1);

    auto ta1_x_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 2);

    auto ta1_x_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 3);

    auto ta1_x_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 4);

    auto ta1_x_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 5);

    auto ta1_x_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 6);

    auto ta1_x_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 7);

    auto ta1_x_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 8);

    auto ta1_x_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 9);

    auto ta1_x_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 10);

    auto ta1_x_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 11);

    auto ta1_x_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 12);

    auto ta1_x_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 13);

    auto ta1_x_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 14);

    auto ta1_y_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 15);

    auto ta1_y_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 16);

    auto ta1_y_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 17);

    auto ta1_y_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 18);

    auto ta1_y_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 19);

    auto ta1_y_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 20);

    auto ta1_y_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 21);

    auto ta1_y_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 22);

    auto ta1_y_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 23);

    auto ta1_y_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 24);

    auto ta1_y_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 25);

    auto ta1_y_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 26);

    auto ta1_y_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 27);

    auto ta1_y_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 28);

    auto ta1_y_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 29);

    auto ta1_z_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 30);

    auto ta1_z_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 31);

    auto ta1_z_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 32);

    auto ta1_z_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 33);

    auto ta1_z_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 34);

    auto ta1_z_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 35);

    auto ta1_z_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 36);

    auto ta1_z_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 37);

    auto ta1_z_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 38);

    auto ta1_z_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 39);

    auto ta1_z_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 40);

    auto ta1_z_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 41);

    auto ta1_z_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 42);

    auto ta1_z_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 43);

    auto ta1_z_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 44);

    // Set up components of auxiliary buffer : PF

    auto ta1_x_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf);

    auto ta1_x_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 1);

    auto ta1_x_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 2);

    auto ta1_x_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 3);

    auto ta1_x_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 4);

    auto ta1_x_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 5);

    auto ta1_x_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 6);

    auto ta1_x_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 7);

    auto ta1_x_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 8);

    auto ta1_x_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 9);

    auto ta1_x_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 10);

    auto ta1_x_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 11);

    auto ta1_x_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 12);

    auto ta1_x_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 13);

    auto ta1_x_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 14);

    auto ta1_x_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 15);

    auto ta1_x_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 16);

    auto ta1_x_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 17);

    auto ta1_x_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 18);

    auto ta1_x_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 19);

    auto ta1_x_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 20);

    auto ta1_x_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 21);

    auto ta1_x_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 22);

    auto ta1_x_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 23);

    auto ta1_x_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 24);

    auto ta1_x_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 25);

    auto ta1_x_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 26);

    auto ta1_x_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 27);

    auto ta1_x_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 28);

    auto ta1_x_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 29);

    auto ta1_y_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 30);

    auto ta1_y_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 31);

    auto ta1_y_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 32);

    auto ta1_y_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 33);

    auto ta1_y_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 34);

    auto ta1_y_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 35);

    auto ta1_y_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 36);

    auto ta1_y_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 37);

    auto ta1_y_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 38);

    auto ta1_y_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 39);

    auto ta1_y_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 40);

    auto ta1_y_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 41);

    auto ta1_y_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 42);

    auto ta1_y_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 43);

    auto ta1_y_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 44);

    auto ta1_y_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 45);

    auto ta1_y_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 46);

    auto ta1_y_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 47);

    auto ta1_y_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 48);

    auto ta1_y_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 49);

    auto ta1_y_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 50);

    auto ta1_y_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 51);

    auto ta1_y_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 52);

    auto ta1_y_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 53);

    auto ta1_y_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 54);

    auto ta1_y_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 55);

    auto ta1_y_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 56);

    auto ta1_y_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 57);

    auto ta1_y_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 58);

    auto ta1_y_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 59);

    auto ta1_z_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 60);

    auto ta1_z_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 61);

    auto ta1_z_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 62);

    auto ta1_z_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 63);

    auto ta1_z_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 64);

    auto ta1_z_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 65);

    auto ta1_z_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 66);

    auto ta1_z_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 67);

    auto ta1_z_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 68);

    auto ta1_z_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 69);

    auto ta1_z_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 70);

    auto ta1_z_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 71);

    auto ta1_z_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 72);

    auto ta1_z_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 73);

    auto ta1_z_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 74);

    auto ta1_z_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 75);

    auto ta1_z_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 76);

    auto ta1_z_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 77);

    auto ta1_z_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 78);

    auto ta1_z_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 79);

    auto ta1_z_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 80);

    auto ta1_z_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 81);

    auto ta1_z_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 82);

    auto ta1_z_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 83);

    auto ta1_z_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 84);

    auto ta1_z_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 85);

    auto ta1_z_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 86);

    auto ta1_z_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 87);

    auto ta1_z_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 88);

    auto ta1_z_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 89);

    // Set up components of auxiliary buffer : PF

    auto ta1_x_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf);

    auto ta1_x_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 1);

    auto ta1_x_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 2);

    auto ta1_x_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 3);

    auto ta1_x_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 4);

    auto ta1_x_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 5);

    auto ta1_x_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 6);

    auto ta1_x_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 7);

    auto ta1_x_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 8);

    auto ta1_x_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 9);

    auto ta1_x_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 10);

    auto ta1_x_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 11);

    auto ta1_x_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 12);

    auto ta1_x_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 13);

    auto ta1_x_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 14);

    auto ta1_x_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 15);

    auto ta1_x_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 16);

    auto ta1_x_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 17);

    auto ta1_x_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 18);

    auto ta1_x_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 19);

    auto ta1_x_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 20);

    auto ta1_x_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 21);

    auto ta1_x_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 22);

    auto ta1_x_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 23);

    auto ta1_x_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 24);

    auto ta1_x_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 25);

    auto ta1_x_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 26);

    auto ta1_x_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 27);

    auto ta1_x_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 28);

    auto ta1_x_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 29);

    auto ta1_y_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 30);

    auto ta1_y_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 31);

    auto ta1_y_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 32);

    auto ta1_y_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 33);

    auto ta1_y_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 34);

    auto ta1_y_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 35);

    auto ta1_y_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 36);

    auto ta1_y_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 37);

    auto ta1_y_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 38);

    auto ta1_y_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 39);

    auto ta1_y_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 40);

    auto ta1_y_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 41);

    auto ta1_y_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 42);

    auto ta1_y_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 43);

    auto ta1_y_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 44);

    auto ta1_y_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 45);

    auto ta1_y_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 46);

    auto ta1_y_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 47);

    auto ta1_y_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 48);

    auto ta1_y_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 49);

    auto ta1_y_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 50);

    auto ta1_y_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 51);

    auto ta1_y_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 52);

    auto ta1_y_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 53);

    auto ta1_y_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 54);

    auto ta1_y_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 55);

    auto ta1_y_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 56);

    auto ta1_y_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 57);

    auto ta1_y_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 58);

    auto ta1_y_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 59);

    auto ta1_z_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 60);

    auto ta1_z_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 61);

    auto ta1_z_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 62);

    auto ta1_z_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 63);

    auto ta1_z_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 64);

    auto ta1_z_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 65);

    auto ta1_z_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 66);

    auto ta1_z_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 67);

    auto ta1_z_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 68);

    auto ta1_z_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 69);

    auto ta1_z_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 70);

    auto ta1_z_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 71);

    auto ta1_z_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 72);

    auto ta1_z_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 73);

    auto ta1_z_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 74);

    auto ta1_z_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 75);

    auto ta1_z_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 76);

    auto ta1_z_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 77);

    auto ta1_z_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 78);

    auto ta1_z_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 79);

    auto ta1_z_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 80);

    auto ta1_z_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 81);

    auto ta1_z_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 82);

    auto ta1_z_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 83);

    auto ta1_z_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 84);

    auto ta1_z_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 85);

    auto ta1_z_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 86);

    auto ta1_z_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 87);

    auto ta1_z_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 88);

    auto ta1_z_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 89);

    // Set up components of auxiliary buffer : PG

    auto ta_x_xxxx_1 = pbuffer.data(idx_npot_1_pg);

    auto ta_x_xxxy_1 = pbuffer.data(idx_npot_1_pg + 1);

    auto ta_x_xxxz_1 = pbuffer.data(idx_npot_1_pg + 2);

    auto ta_x_xxyy_1 = pbuffer.data(idx_npot_1_pg + 3);

    auto ta_x_xxyz_1 = pbuffer.data(idx_npot_1_pg + 4);

    auto ta_x_xxzz_1 = pbuffer.data(idx_npot_1_pg + 5);

    auto ta_x_xyyy_1 = pbuffer.data(idx_npot_1_pg + 6);

    auto ta_x_xyyz_1 = pbuffer.data(idx_npot_1_pg + 7);

    auto ta_x_xyzz_1 = pbuffer.data(idx_npot_1_pg + 8);

    auto ta_x_xzzz_1 = pbuffer.data(idx_npot_1_pg + 9);

    auto ta_x_yyyy_1 = pbuffer.data(idx_npot_1_pg + 10);

    auto ta_x_yyyz_1 = pbuffer.data(idx_npot_1_pg + 11);

    auto ta_x_yyzz_1 = pbuffer.data(idx_npot_1_pg + 12);

    auto ta_x_yzzz_1 = pbuffer.data(idx_npot_1_pg + 13);

    auto ta_x_zzzz_1 = pbuffer.data(idx_npot_1_pg + 14);

    auto ta_y_xxxx_1 = pbuffer.data(idx_npot_1_pg + 15);

    auto ta_y_xxxy_1 = pbuffer.data(idx_npot_1_pg + 16);

    auto ta_y_xxxz_1 = pbuffer.data(idx_npot_1_pg + 17);

    auto ta_y_xxyy_1 = pbuffer.data(idx_npot_1_pg + 18);

    auto ta_y_xxyz_1 = pbuffer.data(idx_npot_1_pg + 19);

    auto ta_y_xxzz_1 = pbuffer.data(idx_npot_1_pg + 20);

    auto ta_y_xyyy_1 = pbuffer.data(idx_npot_1_pg + 21);

    auto ta_y_xyyz_1 = pbuffer.data(idx_npot_1_pg + 22);

    auto ta_y_xyzz_1 = pbuffer.data(idx_npot_1_pg + 23);

    auto ta_y_xzzz_1 = pbuffer.data(idx_npot_1_pg + 24);

    auto ta_y_yyyy_1 = pbuffer.data(idx_npot_1_pg + 25);

    auto ta_y_yyyz_1 = pbuffer.data(idx_npot_1_pg + 26);

    auto ta_y_yyzz_1 = pbuffer.data(idx_npot_1_pg + 27);

    auto ta_y_yzzz_1 = pbuffer.data(idx_npot_1_pg + 28);

    auto ta_y_zzzz_1 = pbuffer.data(idx_npot_1_pg + 29);

    auto ta_z_xxxx_1 = pbuffer.data(idx_npot_1_pg + 30);

    auto ta_z_xxxy_1 = pbuffer.data(idx_npot_1_pg + 31);

    auto ta_z_xxxz_1 = pbuffer.data(idx_npot_1_pg + 32);

    auto ta_z_xxyy_1 = pbuffer.data(idx_npot_1_pg + 33);

    auto ta_z_xxyz_1 = pbuffer.data(idx_npot_1_pg + 34);

    auto ta_z_xxzz_1 = pbuffer.data(idx_npot_1_pg + 35);

    auto ta_z_xyyy_1 = pbuffer.data(idx_npot_1_pg + 36);

    auto ta_z_xyyz_1 = pbuffer.data(idx_npot_1_pg + 37);

    auto ta_z_xyzz_1 = pbuffer.data(idx_npot_1_pg + 38);

    auto ta_z_xzzz_1 = pbuffer.data(idx_npot_1_pg + 39);

    auto ta_z_yyyy_1 = pbuffer.data(idx_npot_1_pg + 40);

    auto ta_z_yyyz_1 = pbuffer.data(idx_npot_1_pg + 41);

    auto ta_z_yyzz_1 = pbuffer.data(idx_npot_1_pg + 42);

    auto ta_z_yzzz_1 = pbuffer.data(idx_npot_1_pg + 43);

    auto ta_z_zzzz_1 = pbuffer.data(idx_npot_1_pg + 44);

    // Set up components of auxiliary buffer : PG

    auto ta1_x_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg);

    auto ta1_x_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 1);

    auto ta1_x_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 2);

    auto ta1_x_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 3);

    auto ta1_x_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 4);

    auto ta1_x_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 5);

    auto ta1_x_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 6);

    auto ta1_x_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 7);

    auto ta1_x_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 8);

    auto ta1_x_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 9);

    auto ta1_x_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 10);

    auto ta1_x_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 11);

    auto ta1_x_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 12);

    auto ta1_x_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 13);

    auto ta1_x_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 14);

    auto ta1_x_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 15);

    auto ta1_x_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 16);

    auto ta1_x_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 17);

    auto ta1_x_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 18);

    auto ta1_x_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 19);

    auto ta1_x_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 20);

    auto ta1_x_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 21);

    auto ta1_x_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 22);

    auto ta1_x_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 23);

    auto ta1_x_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 24);

    auto ta1_x_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 25);

    auto ta1_x_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 26);

    auto ta1_x_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 27);

    auto ta1_x_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 28);

    auto ta1_x_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 29);

    auto ta1_x_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 30);

    auto ta1_x_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 31);

    auto ta1_x_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 32);

    auto ta1_x_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 33);

    auto ta1_x_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 34);

    auto ta1_x_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 35);

    auto ta1_x_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 36);

    auto ta1_x_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 37);

    auto ta1_x_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 38);

    auto ta1_x_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 39);

    auto ta1_x_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 40);

    auto ta1_x_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 41);

    auto ta1_x_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 42);

    auto ta1_x_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 43);

    auto ta1_x_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 44);

    auto ta1_y_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 45);

    auto ta1_y_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 46);

    auto ta1_y_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 47);

    auto ta1_y_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 48);

    auto ta1_y_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 49);

    auto ta1_y_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 50);

    auto ta1_y_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 51);

    auto ta1_y_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 52);

    auto ta1_y_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 53);

    auto ta1_y_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 54);

    auto ta1_y_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 55);

    auto ta1_y_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 56);

    auto ta1_y_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 57);

    auto ta1_y_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 58);

    auto ta1_y_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 59);

    auto ta1_y_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 60);

    auto ta1_y_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 61);

    auto ta1_y_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 62);

    auto ta1_y_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 63);

    auto ta1_y_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 64);

    auto ta1_y_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 65);

    auto ta1_y_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 66);

    auto ta1_y_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 67);

    auto ta1_y_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 68);

    auto ta1_y_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 69);

    auto ta1_y_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 70);

    auto ta1_y_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 71);

    auto ta1_y_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 72);

    auto ta1_y_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 73);

    auto ta1_y_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 74);

    auto ta1_y_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 75);

    auto ta1_y_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 76);

    auto ta1_y_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 77);

    auto ta1_y_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 78);

    auto ta1_y_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 79);

    auto ta1_y_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 80);

    auto ta1_y_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 81);

    auto ta1_y_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 82);

    auto ta1_y_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 83);

    auto ta1_y_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 84);

    auto ta1_y_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 85);

    auto ta1_y_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 86);

    auto ta1_y_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 87);

    auto ta1_y_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 88);

    auto ta1_y_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 89);

    auto ta1_z_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 90);

    auto ta1_z_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 91);

    auto ta1_z_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 92);

    auto ta1_z_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 93);

    auto ta1_z_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 94);

    auto ta1_z_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 95);

    auto ta1_z_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 96);

    auto ta1_z_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 97);

    auto ta1_z_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 98);

    auto ta1_z_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 99);

    auto ta1_z_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 100);

    auto ta1_z_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 101);

    auto ta1_z_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 102);

    auto ta1_z_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 103);

    auto ta1_z_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 104);

    auto ta1_z_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 105);

    auto ta1_z_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 106);

    auto ta1_z_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 107);

    auto ta1_z_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 108);

    auto ta1_z_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 109);

    auto ta1_z_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 110);

    auto ta1_z_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 111);

    auto ta1_z_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 112);

    auto ta1_z_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 113);

    auto ta1_z_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 114);

    auto ta1_z_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 115);

    auto ta1_z_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 116);

    auto ta1_z_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 117);

    auto ta1_z_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 118);

    auto ta1_z_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 119);

    auto ta1_z_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 120);

    auto ta1_z_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 121);

    auto ta1_z_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 122);

    auto ta1_z_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 123);

    auto ta1_z_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 124);

    auto ta1_z_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 125);

    auto ta1_z_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 126);

    auto ta1_z_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 127);

    auto ta1_z_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 128);

    auto ta1_z_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 129);

    auto ta1_z_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 130);

    auto ta1_z_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 131);

    auto ta1_z_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 132);

    auto ta1_z_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 133);

    auto ta1_z_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 134);

    // Set up components of auxiliary buffer : PG

    auto ta1_x_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg);

    auto ta1_x_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 1);

    auto ta1_x_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 2);

    auto ta1_x_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 3);

    auto ta1_x_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 4);

    auto ta1_x_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 5);

    auto ta1_x_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 6);

    auto ta1_x_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 7);

    auto ta1_x_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 8);

    auto ta1_x_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 9);

    auto ta1_x_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 10);

    auto ta1_x_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 11);

    auto ta1_x_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 12);

    auto ta1_x_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 13);

    auto ta1_x_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 14);

    auto ta1_x_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 15);

    auto ta1_x_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 16);

    auto ta1_x_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 17);

    auto ta1_x_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 18);

    auto ta1_x_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 19);

    auto ta1_x_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 20);

    auto ta1_x_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 21);

    auto ta1_x_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 22);

    auto ta1_x_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 23);

    auto ta1_x_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 24);

    auto ta1_x_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 25);

    auto ta1_x_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 26);

    auto ta1_x_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 27);

    auto ta1_x_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 28);

    auto ta1_x_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 29);

    auto ta1_x_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 30);

    auto ta1_x_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 31);

    auto ta1_x_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 32);

    auto ta1_x_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 33);

    auto ta1_x_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 34);

    auto ta1_x_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 35);

    auto ta1_x_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 36);

    auto ta1_x_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 37);

    auto ta1_x_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 38);

    auto ta1_x_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 39);

    auto ta1_x_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 40);

    auto ta1_x_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 41);

    auto ta1_x_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 42);

    auto ta1_x_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 43);

    auto ta1_x_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 44);

    auto ta1_y_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 45);

    auto ta1_y_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 46);

    auto ta1_y_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 47);

    auto ta1_y_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 48);

    auto ta1_y_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 49);

    auto ta1_y_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 50);

    auto ta1_y_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 51);

    auto ta1_y_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 52);

    auto ta1_y_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 53);

    auto ta1_y_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 54);

    auto ta1_y_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 55);

    auto ta1_y_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 56);

    auto ta1_y_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 57);

    auto ta1_y_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 58);

    auto ta1_y_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 59);

    auto ta1_y_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 60);

    auto ta1_y_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 61);

    auto ta1_y_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 62);

    auto ta1_y_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 63);

    auto ta1_y_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 64);

    auto ta1_y_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 65);

    auto ta1_y_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 66);

    auto ta1_y_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 67);

    auto ta1_y_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 68);

    auto ta1_y_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 69);

    auto ta1_y_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 70);

    auto ta1_y_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 71);

    auto ta1_y_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 72);

    auto ta1_y_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 73);

    auto ta1_y_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 74);

    auto ta1_y_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 75);

    auto ta1_y_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 76);

    auto ta1_y_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 77);

    auto ta1_y_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 78);

    auto ta1_y_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 79);

    auto ta1_y_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 80);

    auto ta1_y_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 81);

    auto ta1_y_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 82);

    auto ta1_y_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 83);

    auto ta1_y_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 84);

    auto ta1_y_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 85);

    auto ta1_y_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 86);

    auto ta1_y_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 87);

    auto ta1_y_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 88);

    auto ta1_y_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 89);

    auto ta1_z_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 90);

    auto ta1_z_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 91);

    auto ta1_z_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 92);

    auto ta1_z_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 93);

    auto ta1_z_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 94);

    auto ta1_z_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 95);

    auto ta1_z_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 96);

    auto ta1_z_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 97);

    auto ta1_z_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 98);

    auto ta1_z_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 99);

    auto ta1_z_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 100);

    auto ta1_z_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 101);

    auto ta1_z_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 102);

    auto ta1_z_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 103);

    auto ta1_z_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 104);

    auto ta1_z_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 105);

    auto ta1_z_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 106);

    auto ta1_z_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 107);

    auto ta1_z_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 108);

    auto ta1_z_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 109);

    auto ta1_z_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 110);

    auto ta1_z_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 111);

    auto ta1_z_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 112);

    auto ta1_z_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 113);

    auto ta1_z_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 114);

    auto ta1_z_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 115);

    auto ta1_z_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 116);

    auto ta1_z_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 117);

    auto ta1_z_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 118);

    auto ta1_z_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 119);

    auto ta1_z_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 120);

    auto ta1_z_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 121);

    auto ta1_z_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 122);

    auto ta1_z_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 123);

    auto ta1_z_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 124);

    auto ta1_z_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 125);

    auto ta1_z_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 126);

    auto ta1_z_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 127);

    auto ta1_z_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 128);

    auto ta1_z_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 129);

    auto ta1_z_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 130);

    auto ta1_z_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 131);

    auto ta1_z_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 132);

    auto ta1_z_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 133);

    auto ta1_z_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 134);

    // Set up 0-15 components of targeted buffer : DG

    auto ta1_x_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg);

    auto ta1_x_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 1);

    auto ta1_x_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 2);

    auto ta1_x_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 3);

    auto ta1_x_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 4);

    auto ta1_x_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 5);

    auto ta1_x_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 6);

    auto ta1_x_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 7);

    auto ta1_x_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 8);

    auto ta1_x_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 9);

    auto ta1_x_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 10);

    auto ta1_x_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 11);

    auto ta1_x_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 12);

    auto ta1_x_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 13);

    auto ta1_x_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 14);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_x_0_xxxx_0,  \
                             ta1_x_0_xxxx_1,  \
                             ta1_x_0_xxxy_0,  \
                             ta1_x_0_xxxy_1,  \
                             ta1_x_0_xxxz_0,  \
                             ta1_x_0_xxxz_1,  \
                             ta1_x_0_xxyy_0,  \
                             ta1_x_0_xxyy_1,  \
                             ta1_x_0_xxyz_0,  \
                             ta1_x_0_xxyz_1,  \
                             ta1_x_0_xxzz_0,  \
                             ta1_x_0_xxzz_1,  \
                             ta1_x_0_xyyy_0,  \
                             ta1_x_0_xyyy_1,  \
                             ta1_x_0_xyyz_0,  \
                             ta1_x_0_xyyz_1,  \
                             ta1_x_0_xyzz_0,  \
                             ta1_x_0_xyzz_1,  \
                             ta1_x_0_xzzz_0,  \
                             ta1_x_0_xzzz_1,  \
                             ta1_x_0_yyyy_0,  \
                             ta1_x_0_yyyy_1,  \
                             ta1_x_0_yyyz_0,  \
                             ta1_x_0_yyyz_1,  \
                             ta1_x_0_yyzz_0,  \
                             ta1_x_0_yyzz_1,  \
                             ta1_x_0_yzzz_0,  \
                             ta1_x_0_yzzz_1,  \
                             ta1_x_0_zzzz_0,  \
                             ta1_x_0_zzzz_1,  \
                             ta1_x_x_xxx_0,   \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxxx_0,  \
                             ta1_x_x_xxxx_1,  \
                             ta1_x_x_xxxy_0,  \
                             ta1_x_x_xxxy_1,  \
                             ta1_x_x_xxxz_0,  \
                             ta1_x_x_xxxz_1,  \
                             ta1_x_x_xxy_0,   \
                             ta1_x_x_xxy_1,   \
                             ta1_x_x_xxyy_0,  \
                             ta1_x_x_xxyy_1,  \
                             ta1_x_x_xxyz_0,  \
                             ta1_x_x_xxyz_1,  \
                             ta1_x_x_xxz_0,   \
                             ta1_x_x_xxz_1,   \
                             ta1_x_x_xxzz_0,  \
                             ta1_x_x_xxzz_1,  \
                             ta1_x_x_xyy_0,   \
                             ta1_x_x_xyy_1,   \
                             ta1_x_x_xyyy_0,  \
                             ta1_x_x_xyyy_1,  \
                             ta1_x_x_xyyz_0,  \
                             ta1_x_x_xyyz_1,  \
                             ta1_x_x_xyz_0,   \
                             ta1_x_x_xyz_1,   \
                             ta1_x_x_xyzz_0,  \
                             ta1_x_x_xyzz_1,  \
                             ta1_x_x_xzz_0,   \
                             ta1_x_x_xzz_1,   \
                             ta1_x_x_xzzz_0,  \
                             ta1_x_x_xzzz_1,  \
                             ta1_x_x_yyy_0,   \
                             ta1_x_x_yyy_1,   \
                             ta1_x_x_yyyy_0,  \
                             ta1_x_x_yyyy_1,  \
                             ta1_x_x_yyyz_0,  \
                             ta1_x_x_yyyz_1,  \
                             ta1_x_x_yyz_0,   \
                             ta1_x_x_yyz_1,   \
                             ta1_x_x_yyzz_0,  \
                             ta1_x_x_yyzz_1,  \
                             ta1_x_x_yzz_0,   \
                             ta1_x_x_yzz_1,   \
                             ta1_x_x_yzzz_0,  \
                             ta1_x_x_yzzz_1,  \
                             ta1_x_x_zzz_0,   \
                             ta1_x_x_zzz_1,   \
                             ta1_x_x_zzzz_0,  \
                             ta1_x_x_zzzz_1,  \
                             ta1_x_xx_xxxx_0, \
                             ta1_x_xx_xxxy_0, \
                             ta1_x_xx_xxxz_0, \
                             ta1_x_xx_xxyy_0, \
                             ta1_x_xx_xxyz_0, \
                             ta1_x_xx_xxzz_0, \
                             ta1_x_xx_xyyy_0, \
                             ta1_x_xx_xyyz_0, \
                             ta1_x_xx_xyzz_0, \
                             ta1_x_xx_xzzz_0, \
                             ta1_x_xx_yyyy_0, \
                             ta1_x_xx_yyyz_0, \
                             ta1_x_xx_yyzz_0, \
                             ta1_x_xx_yzzz_0, \
                             ta1_x_xx_zzzz_0, \
                             ta_x_xxxx_1,     \
                             ta_x_xxxy_1,     \
                             ta_x_xxxz_1,     \
                             ta_x_xxyy_1,     \
                             ta_x_xxyz_1,     \
                             ta_x_xxzz_1,     \
                             ta_x_xyyy_1,     \
                             ta_x_xyyz_1,     \
                             ta_x_xyzz_1,     \
                             ta_x_xzzz_1,     \
                             ta_x_yyyy_1,     \
                             ta_x_yyyz_1,     \
                             ta_x_yyzz_1,     \
                             ta_x_yzzz_1,     \
                             ta_x_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xx_xxxx_0[i] = ta1_x_0_xxxx_0[i] * fe_0 - ta1_x_0_xxxx_1[i] * fe_0 + 4.0 * ta1_x_x_xxx_0[i] * fe_0 - 4.0 * ta1_x_x_xxx_1[i] * fe_0 +
                             ta_x_xxxx_1[i] + ta1_x_x_xxxx_0[i] * pa_x[i] - ta1_x_x_xxxx_1[i] * pc_x[i];

        ta1_x_xx_xxxy_0[i] = ta1_x_0_xxxy_0[i] * fe_0 - ta1_x_0_xxxy_1[i] * fe_0 + 3.0 * ta1_x_x_xxy_0[i] * fe_0 - 3.0 * ta1_x_x_xxy_1[i] * fe_0 +
                             ta_x_xxxy_1[i] + ta1_x_x_xxxy_0[i] * pa_x[i] - ta1_x_x_xxxy_1[i] * pc_x[i];

        ta1_x_xx_xxxz_0[i] = ta1_x_0_xxxz_0[i] * fe_0 - ta1_x_0_xxxz_1[i] * fe_0 + 3.0 * ta1_x_x_xxz_0[i] * fe_0 - 3.0 * ta1_x_x_xxz_1[i] * fe_0 +
                             ta_x_xxxz_1[i] + ta1_x_x_xxxz_0[i] * pa_x[i] - ta1_x_x_xxxz_1[i] * pc_x[i];

        ta1_x_xx_xxyy_0[i] = ta1_x_0_xxyy_0[i] * fe_0 - ta1_x_0_xxyy_1[i] * fe_0 + 2.0 * ta1_x_x_xyy_0[i] * fe_0 - 2.0 * ta1_x_x_xyy_1[i] * fe_0 +
                             ta_x_xxyy_1[i] + ta1_x_x_xxyy_0[i] * pa_x[i] - ta1_x_x_xxyy_1[i] * pc_x[i];

        ta1_x_xx_xxyz_0[i] = ta1_x_0_xxyz_0[i] * fe_0 - ta1_x_0_xxyz_1[i] * fe_0 + 2.0 * ta1_x_x_xyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyz_1[i] * fe_0 +
                             ta_x_xxyz_1[i] + ta1_x_x_xxyz_0[i] * pa_x[i] - ta1_x_x_xxyz_1[i] * pc_x[i];

        ta1_x_xx_xxzz_0[i] = ta1_x_0_xxzz_0[i] * fe_0 - ta1_x_0_xxzz_1[i] * fe_0 + 2.0 * ta1_x_x_xzz_0[i] * fe_0 - 2.0 * ta1_x_x_xzz_1[i] * fe_0 +
                             ta_x_xxzz_1[i] + ta1_x_x_xxzz_0[i] * pa_x[i] - ta1_x_x_xxzz_1[i] * pc_x[i];

        ta1_x_xx_xyyy_0[i] = ta1_x_0_xyyy_0[i] * fe_0 - ta1_x_0_xyyy_1[i] * fe_0 + ta1_x_x_yyy_0[i] * fe_0 - ta1_x_x_yyy_1[i] * fe_0 +
                             ta_x_xyyy_1[i] + ta1_x_x_xyyy_0[i] * pa_x[i] - ta1_x_x_xyyy_1[i] * pc_x[i];

        ta1_x_xx_xyyz_0[i] = ta1_x_0_xyyz_0[i] * fe_0 - ta1_x_0_xyyz_1[i] * fe_0 + ta1_x_x_yyz_0[i] * fe_0 - ta1_x_x_yyz_1[i] * fe_0 +
                             ta_x_xyyz_1[i] + ta1_x_x_xyyz_0[i] * pa_x[i] - ta1_x_x_xyyz_1[i] * pc_x[i];

        ta1_x_xx_xyzz_0[i] = ta1_x_0_xyzz_0[i] * fe_0 - ta1_x_0_xyzz_1[i] * fe_0 + ta1_x_x_yzz_0[i] * fe_0 - ta1_x_x_yzz_1[i] * fe_0 +
                             ta_x_xyzz_1[i] + ta1_x_x_xyzz_0[i] * pa_x[i] - ta1_x_x_xyzz_1[i] * pc_x[i];

        ta1_x_xx_xzzz_0[i] = ta1_x_0_xzzz_0[i] * fe_0 - ta1_x_0_xzzz_1[i] * fe_0 + ta1_x_x_zzz_0[i] * fe_0 - ta1_x_x_zzz_1[i] * fe_0 +
                             ta_x_xzzz_1[i] + ta1_x_x_xzzz_0[i] * pa_x[i] - ta1_x_x_xzzz_1[i] * pc_x[i];

        ta1_x_xx_yyyy_0[i] =
            ta1_x_0_yyyy_0[i] * fe_0 - ta1_x_0_yyyy_1[i] * fe_0 + ta_x_yyyy_1[i] + ta1_x_x_yyyy_0[i] * pa_x[i] - ta1_x_x_yyyy_1[i] * pc_x[i];

        ta1_x_xx_yyyz_0[i] =
            ta1_x_0_yyyz_0[i] * fe_0 - ta1_x_0_yyyz_1[i] * fe_0 + ta_x_yyyz_1[i] + ta1_x_x_yyyz_0[i] * pa_x[i] - ta1_x_x_yyyz_1[i] * pc_x[i];

        ta1_x_xx_yyzz_0[i] =
            ta1_x_0_yyzz_0[i] * fe_0 - ta1_x_0_yyzz_1[i] * fe_0 + ta_x_yyzz_1[i] + ta1_x_x_yyzz_0[i] * pa_x[i] - ta1_x_x_yyzz_1[i] * pc_x[i];

        ta1_x_xx_yzzz_0[i] =
            ta1_x_0_yzzz_0[i] * fe_0 - ta1_x_0_yzzz_1[i] * fe_0 + ta_x_yzzz_1[i] + ta1_x_x_yzzz_0[i] * pa_x[i] - ta1_x_x_yzzz_1[i] * pc_x[i];

        ta1_x_xx_zzzz_0[i] =
            ta1_x_0_zzzz_0[i] * fe_0 - ta1_x_0_zzzz_1[i] * fe_0 + ta_x_zzzz_1[i] + ta1_x_x_zzzz_0[i] * pa_x[i] - ta1_x_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : DG

    auto ta1_x_xy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 15);

    auto ta1_x_xy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 16);

    auto ta1_x_xy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 17);

    auto ta1_x_xy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 18);

    auto ta1_x_xy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 19);

    auto ta1_x_xy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 20);

    auto ta1_x_xy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 21);

    auto ta1_x_xy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 22);

    auto ta1_x_xy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 23);

    auto ta1_x_xy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 24);

    auto ta1_x_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 25);

    auto ta1_x_xy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 26);

    auto ta1_x_xy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 27);

    auto ta1_x_xy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 28);

    auto ta1_x_xy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 29);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_x_xxx_0,   \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxxx_0,  \
                             ta1_x_x_xxxx_1,  \
                             ta1_x_x_xxxy_0,  \
                             ta1_x_x_xxxy_1,  \
                             ta1_x_x_xxxz_0,  \
                             ta1_x_x_xxxz_1,  \
                             ta1_x_x_xxy_0,   \
                             ta1_x_x_xxy_1,   \
                             ta1_x_x_xxyy_0,  \
                             ta1_x_x_xxyy_1,  \
                             ta1_x_x_xxyz_0,  \
                             ta1_x_x_xxyz_1,  \
                             ta1_x_x_xxz_0,   \
                             ta1_x_x_xxz_1,   \
                             ta1_x_x_xxzz_0,  \
                             ta1_x_x_xxzz_1,  \
                             ta1_x_x_xyy_0,   \
                             ta1_x_x_xyy_1,   \
                             ta1_x_x_xyyy_0,  \
                             ta1_x_x_xyyy_1,  \
                             ta1_x_x_xyyz_0,  \
                             ta1_x_x_xyyz_1,  \
                             ta1_x_x_xyz_0,   \
                             ta1_x_x_xyz_1,   \
                             ta1_x_x_xyzz_0,  \
                             ta1_x_x_xyzz_1,  \
                             ta1_x_x_xzz_0,   \
                             ta1_x_x_xzz_1,   \
                             ta1_x_x_xzzz_0,  \
                             ta1_x_x_xzzz_1,  \
                             ta1_x_x_zzzz_0,  \
                             ta1_x_x_zzzz_1,  \
                             ta1_x_xy_xxxx_0, \
                             ta1_x_xy_xxxy_0, \
                             ta1_x_xy_xxxz_0, \
                             ta1_x_xy_xxyy_0, \
                             ta1_x_xy_xxyz_0, \
                             ta1_x_xy_xxzz_0, \
                             ta1_x_xy_xyyy_0, \
                             ta1_x_xy_xyyz_0, \
                             ta1_x_xy_xyzz_0, \
                             ta1_x_xy_xzzz_0, \
                             ta1_x_xy_yyyy_0, \
                             ta1_x_xy_yyyz_0, \
                             ta1_x_xy_yyzz_0, \
                             ta1_x_xy_yzzz_0, \
                             ta1_x_xy_zzzz_0, \
                             ta1_x_y_yyyy_0,  \
                             ta1_x_y_yyyy_1,  \
                             ta1_x_y_yyyz_0,  \
                             ta1_x_y_yyyz_1,  \
                             ta1_x_y_yyzz_0,  \
                             ta1_x_y_yyzz_1,  \
                             ta1_x_y_yzzz_0,  \
                             ta1_x_y_yzzz_1,  \
                             ta_y_yyyy_1,     \
                             ta_y_yyyz_1,     \
                             ta_y_yyzz_1,     \
                             ta_y_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xy_xxxx_0[i] = ta1_x_x_xxxx_0[i] * pa_y[i] - ta1_x_x_xxxx_1[i] * pc_y[i];

        ta1_x_xy_xxxy_0[i] = ta1_x_x_xxx_0[i] * fe_0 - ta1_x_x_xxx_1[i] * fe_0 + ta1_x_x_xxxy_0[i] * pa_y[i] - ta1_x_x_xxxy_1[i] * pc_y[i];

        ta1_x_xy_xxxz_0[i] = ta1_x_x_xxxz_0[i] * pa_y[i] - ta1_x_x_xxxz_1[i] * pc_y[i];

        ta1_x_xy_xxyy_0[i] =
            2.0 * ta1_x_x_xxy_0[i] * fe_0 - 2.0 * ta1_x_x_xxy_1[i] * fe_0 + ta1_x_x_xxyy_0[i] * pa_y[i] - ta1_x_x_xxyy_1[i] * pc_y[i];

        ta1_x_xy_xxyz_0[i] = ta1_x_x_xxz_0[i] * fe_0 - ta1_x_x_xxz_1[i] * fe_0 + ta1_x_x_xxyz_0[i] * pa_y[i] - ta1_x_x_xxyz_1[i] * pc_y[i];

        ta1_x_xy_xxzz_0[i] = ta1_x_x_xxzz_0[i] * pa_y[i] - ta1_x_x_xxzz_1[i] * pc_y[i];

        ta1_x_xy_xyyy_0[i] =
            3.0 * ta1_x_x_xyy_0[i] * fe_0 - 3.0 * ta1_x_x_xyy_1[i] * fe_0 + ta1_x_x_xyyy_0[i] * pa_y[i] - ta1_x_x_xyyy_1[i] * pc_y[i];

        ta1_x_xy_xyyz_0[i] =
            2.0 * ta1_x_x_xyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyz_1[i] * fe_0 + ta1_x_x_xyyz_0[i] * pa_y[i] - ta1_x_x_xyyz_1[i] * pc_y[i];

        ta1_x_xy_xyzz_0[i] = ta1_x_x_xzz_0[i] * fe_0 - ta1_x_x_xzz_1[i] * fe_0 + ta1_x_x_xyzz_0[i] * pa_y[i] - ta1_x_x_xyzz_1[i] * pc_y[i];

        ta1_x_xy_xzzz_0[i] = ta1_x_x_xzzz_0[i] * pa_y[i] - ta1_x_x_xzzz_1[i] * pc_y[i];

        ta1_x_xy_yyyy_0[i] = ta_y_yyyy_1[i] + ta1_x_y_yyyy_0[i] * pa_x[i] - ta1_x_y_yyyy_1[i] * pc_x[i];

        ta1_x_xy_yyyz_0[i] = ta_y_yyyz_1[i] + ta1_x_y_yyyz_0[i] * pa_x[i] - ta1_x_y_yyyz_1[i] * pc_x[i];

        ta1_x_xy_yyzz_0[i] = ta_y_yyzz_1[i] + ta1_x_y_yyzz_0[i] * pa_x[i] - ta1_x_y_yyzz_1[i] * pc_x[i];

        ta1_x_xy_yzzz_0[i] = ta_y_yzzz_1[i] + ta1_x_y_yzzz_0[i] * pa_x[i] - ta1_x_y_yzzz_1[i] * pc_x[i];

        ta1_x_xy_zzzz_0[i] = ta1_x_x_zzzz_0[i] * pa_y[i] - ta1_x_x_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : DG

    auto ta1_x_xz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 30);

    auto ta1_x_xz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 31);

    auto ta1_x_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 32);

    auto ta1_x_xz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 33);

    auto ta1_x_xz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 34);

    auto ta1_x_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 35);

    auto ta1_x_xz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 36);

    auto ta1_x_xz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 37);

    auto ta1_x_xz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 38);

    auto ta1_x_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 39);

    auto ta1_x_xz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 40);

    auto ta1_x_xz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 41);

    auto ta1_x_xz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 42);

    auto ta1_x_xz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 43);

    auto ta1_x_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 44);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_x_xxx_0,   \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxxx_0,  \
                             ta1_x_x_xxxx_1,  \
                             ta1_x_x_xxxy_0,  \
                             ta1_x_x_xxxy_1,  \
                             ta1_x_x_xxxz_0,  \
                             ta1_x_x_xxxz_1,  \
                             ta1_x_x_xxy_0,   \
                             ta1_x_x_xxy_1,   \
                             ta1_x_x_xxyy_0,  \
                             ta1_x_x_xxyy_1,  \
                             ta1_x_x_xxyz_0,  \
                             ta1_x_x_xxyz_1,  \
                             ta1_x_x_xxz_0,   \
                             ta1_x_x_xxz_1,   \
                             ta1_x_x_xxzz_0,  \
                             ta1_x_x_xxzz_1,  \
                             ta1_x_x_xyy_0,   \
                             ta1_x_x_xyy_1,   \
                             ta1_x_x_xyyy_0,  \
                             ta1_x_x_xyyy_1,  \
                             ta1_x_x_xyyz_0,  \
                             ta1_x_x_xyyz_1,  \
                             ta1_x_x_xyz_0,   \
                             ta1_x_x_xyz_1,   \
                             ta1_x_x_xyzz_0,  \
                             ta1_x_x_xyzz_1,  \
                             ta1_x_x_xzz_0,   \
                             ta1_x_x_xzz_1,   \
                             ta1_x_x_xzzz_0,  \
                             ta1_x_x_xzzz_1,  \
                             ta1_x_x_yyyy_0,  \
                             ta1_x_x_yyyy_1,  \
                             ta1_x_xz_xxxx_0, \
                             ta1_x_xz_xxxy_0, \
                             ta1_x_xz_xxxz_0, \
                             ta1_x_xz_xxyy_0, \
                             ta1_x_xz_xxyz_0, \
                             ta1_x_xz_xxzz_0, \
                             ta1_x_xz_xyyy_0, \
                             ta1_x_xz_xyyz_0, \
                             ta1_x_xz_xyzz_0, \
                             ta1_x_xz_xzzz_0, \
                             ta1_x_xz_yyyy_0, \
                             ta1_x_xz_yyyz_0, \
                             ta1_x_xz_yyzz_0, \
                             ta1_x_xz_yzzz_0, \
                             ta1_x_xz_zzzz_0, \
                             ta1_x_z_yyyz_0,  \
                             ta1_x_z_yyyz_1,  \
                             ta1_x_z_yyzz_0,  \
                             ta1_x_z_yyzz_1,  \
                             ta1_x_z_yzzz_0,  \
                             ta1_x_z_yzzz_1,  \
                             ta1_x_z_zzzz_0,  \
                             ta1_x_z_zzzz_1,  \
                             ta_z_yyyz_1,     \
                             ta_z_yyzz_1,     \
                             ta_z_yzzz_1,     \
                             ta_z_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xz_xxxx_0[i] = ta1_x_x_xxxx_0[i] * pa_z[i] - ta1_x_x_xxxx_1[i] * pc_z[i];

        ta1_x_xz_xxxy_0[i] = ta1_x_x_xxxy_0[i] * pa_z[i] - ta1_x_x_xxxy_1[i] * pc_z[i];

        ta1_x_xz_xxxz_0[i] = ta1_x_x_xxx_0[i] * fe_0 - ta1_x_x_xxx_1[i] * fe_0 + ta1_x_x_xxxz_0[i] * pa_z[i] - ta1_x_x_xxxz_1[i] * pc_z[i];

        ta1_x_xz_xxyy_0[i] = ta1_x_x_xxyy_0[i] * pa_z[i] - ta1_x_x_xxyy_1[i] * pc_z[i];

        ta1_x_xz_xxyz_0[i] = ta1_x_x_xxy_0[i] * fe_0 - ta1_x_x_xxy_1[i] * fe_0 + ta1_x_x_xxyz_0[i] * pa_z[i] - ta1_x_x_xxyz_1[i] * pc_z[i];

        ta1_x_xz_xxzz_0[i] =
            2.0 * ta1_x_x_xxz_0[i] * fe_0 - 2.0 * ta1_x_x_xxz_1[i] * fe_0 + ta1_x_x_xxzz_0[i] * pa_z[i] - ta1_x_x_xxzz_1[i] * pc_z[i];

        ta1_x_xz_xyyy_0[i] = ta1_x_x_xyyy_0[i] * pa_z[i] - ta1_x_x_xyyy_1[i] * pc_z[i];

        ta1_x_xz_xyyz_0[i] = ta1_x_x_xyy_0[i] * fe_0 - ta1_x_x_xyy_1[i] * fe_0 + ta1_x_x_xyyz_0[i] * pa_z[i] - ta1_x_x_xyyz_1[i] * pc_z[i];

        ta1_x_xz_xyzz_0[i] =
            2.0 * ta1_x_x_xyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyz_1[i] * fe_0 + ta1_x_x_xyzz_0[i] * pa_z[i] - ta1_x_x_xyzz_1[i] * pc_z[i];

        ta1_x_xz_xzzz_0[i] =
            3.0 * ta1_x_x_xzz_0[i] * fe_0 - 3.0 * ta1_x_x_xzz_1[i] * fe_0 + ta1_x_x_xzzz_0[i] * pa_z[i] - ta1_x_x_xzzz_1[i] * pc_z[i];

        ta1_x_xz_yyyy_0[i] = ta1_x_x_yyyy_0[i] * pa_z[i] - ta1_x_x_yyyy_1[i] * pc_z[i];

        ta1_x_xz_yyyz_0[i] = ta_z_yyyz_1[i] + ta1_x_z_yyyz_0[i] * pa_x[i] - ta1_x_z_yyyz_1[i] * pc_x[i];

        ta1_x_xz_yyzz_0[i] = ta_z_yyzz_1[i] + ta1_x_z_yyzz_0[i] * pa_x[i] - ta1_x_z_yyzz_1[i] * pc_x[i];

        ta1_x_xz_yzzz_0[i] = ta_z_yzzz_1[i] + ta1_x_z_yzzz_0[i] * pa_x[i] - ta1_x_z_yzzz_1[i] * pc_x[i];

        ta1_x_xz_zzzz_0[i] = ta_z_zzzz_1[i] + ta1_x_z_zzzz_0[i] * pa_x[i] - ta1_x_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 45-60 components of targeted buffer : DG

    auto ta1_x_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 45);

    auto ta1_x_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 46);

    auto ta1_x_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 47);

    auto ta1_x_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 48);

    auto ta1_x_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 49);

    auto ta1_x_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 50);

    auto ta1_x_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 51);

    auto ta1_x_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 52);

    auto ta1_x_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 53);

    auto ta1_x_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 54);

    auto ta1_x_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 55);

    auto ta1_x_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 56);

    auto ta1_x_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 57);

    auto ta1_x_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 58);

    auto ta1_x_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 59);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_0_xxxx_0,  \
                             ta1_x_0_xxxx_1,  \
                             ta1_x_0_xxxy_0,  \
                             ta1_x_0_xxxy_1,  \
                             ta1_x_0_xxxz_0,  \
                             ta1_x_0_xxxz_1,  \
                             ta1_x_0_xxyy_0,  \
                             ta1_x_0_xxyy_1,  \
                             ta1_x_0_xxyz_0,  \
                             ta1_x_0_xxyz_1,  \
                             ta1_x_0_xxzz_0,  \
                             ta1_x_0_xxzz_1,  \
                             ta1_x_0_xyyy_0,  \
                             ta1_x_0_xyyy_1,  \
                             ta1_x_0_xyyz_0,  \
                             ta1_x_0_xyyz_1,  \
                             ta1_x_0_xyzz_0,  \
                             ta1_x_0_xyzz_1,  \
                             ta1_x_0_xzzz_0,  \
                             ta1_x_0_xzzz_1,  \
                             ta1_x_0_yyyy_0,  \
                             ta1_x_0_yyyy_1,  \
                             ta1_x_0_yyyz_0,  \
                             ta1_x_0_yyyz_1,  \
                             ta1_x_0_yyzz_0,  \
                             ta1_x_0_yyzz_1,  \
                             ta1_x_0_yzzz_0,  \
                             ta1_x_0_yzzz_1,  \
                             ta1_x_0_zzzz_0,  \
                             ta1_x_0_zzzz_1,  \
                             ta1_x_y_xxx_0,   \
                             ta1_x_y_xxx_1,   \
                             ta1_x_y_xxxx_0,  \
                             ta1_x_y_xxxx_1,  \
                             ta1_x_y_xxxy_0,  \
                             ta1_x_y_xxxy_1,  \
                             ta1_x_y_xxxz_0,  \
                             ta1_x_y_xxxz_1,  \
                             ta1_x_y_xxy_0,   \
                             ta1_x_y_xxy_1,   \
                             ta1_x_y_xxyy_0,  \
                             ta1_x_y_xxyy_1,  \
                             ta1_x_y_xxyz_0,  \
                             ta1_x_y_xxyz_1,  \
                             ta1_x_y_xxz_0,   \
                             ta1_x_y_xxz_1,   \
                             ta1_x_y_xxzz_0,  \
                             ta1_x_y_xxzz_1,  \
                             ta1_x_y_xyy_0,   \
                             ta1_x_y_xyy_1,   \
                             ta1_x_y_xyyy_0,  \
                             ta1_x_y_xyyy_1,  \
                             ta1_x_y_xyyz_0,  \
                             ta1_x_y_xyyz_1,  \
                             ta1_x_y_xyz_0,   \
                             ta1_x_y_xyz_1,   \
                             ta1_x_y_xyzz_0,  \
                             ta1_x_y_xyzz_1,  \
                             ta1_x_y_xzz_0,   \
                             ta1_x_y_xzz_1,   \
                             ta1_x_y_xzzz_0,  \
                             ta1_x_y_xzzz_1,  \
                             ta1_x_y_yyy_0,   \
                             ta1_x_y_yyy_1,   \
                             ta1_x_y_yyyy_0,  \
                             ta1_x_y_yyyy_1,  \
                             ta1_x_y_yyyz_0,  \
                             ta1_x_y_yyyz_1,  \
                             ta1_x_y_yyz_0,   \
                             ta1_x_y_yyz_1,   \
                             ta1_x_y_yyzz_0,  \
                             ta1_x_y_yyzz_1,  \
                             ta1_x_y_yzz_0,   \
                             ta1_x_y_yzz_1,   \
                             ta1_x_y_yzzz_0,  \
                             ta1_x_y_yzzz_1,  \
                             ta1_x_y_zzz_0,   \
                             ta1_x_y_zzz_1,   \
                             ta1_x_y_zzzz_0,  \
                             ta1_x_y_zzzz_1,  \
                             ta1_x_yy_xxxx_0, \
                             ta1_x_yy_xxxy_0, \
                             ta1_x_yy_xxxz_0, \
                             ta1_x_yy_xxyy_0, \
                             ta1_x_yy_xxyz_0, \
                             ta1_x_yy_xxzz_0, \
                             ta1_x_yy_xyyy_0, \
                             ta1_x_yy_xyyz_0, \
                             ta1_x_yy_xyzz_0, \
                             ta1_x_yy_xzzz_0, \
                             ta1_x_yy_yyyy_0, \
                             ta1_x_yy_yyyz_0, \
                             ta1_x_yy_yyzz_0, \
                             ta1_x_yy_yzzz_0, \
                             ta1_x_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yy_xxxx_0[i] = ta1_x_0_xxxx_0[i] * fe_0 - ta1_x_0_xxxx_1[i] * fe_0 + ta1_x_y_xxxx_0[i] * pa_y[i] - ta1_x_y_xxxx_1[i] * pc_y[i];

        ta1_x_yy_xxxy_0[i] = ta1_x_0_xxxy_0[i] * fe_0 - ta1_x_0_xxxy_1[i] * fe_0 + ta1_x_y_xxx_0[i] * fe_0 - ta1_x_y_xxx_1[i] * fe_0 +
                             ta1_x_y_xxxy_0[i] * pa_y[i] - ta1_x_y_xxxy_1[i] * pc_y[i];

        ta1_x_yy_xxxz_0[i] = ta1_x_0_xxxz_0[i] * fe_0 - ta1_x_0_xxxz_1[i] * fe_0 + ta1_x_y_xxxz_0[i] * pa_y[i] - ta1_x_y_xxxz_1[i] * pc_y[i];

        ta1_x_yy_xxyy_0[i] = ta1_x_0_xxyy_0[i] * fe_0 - ta1_x_0_xxyy_1[i] * fe_0 + 2.0 * ta1_x_y_xxy_0[i] * fe_0 - 2.0 * ta1_x_y_xxy_1[i] * fe_0 +
                             ta1_x_y_xxyy_0[i] * pa_y[i] - ta1_x_y_xxyy_1[i] * pc_y[i];

        ta1_x_yy_xxyz_0[i] = ta1_x_0_xxyz_0[i] * fe_0 - ta1_x_0_xxyz_1[i] * fe_0 + ta1_x_y_xxz_0[i] * fe_0 - ta1_x_y_xxz_1[i] * fe_0 +
                             ta1_x_y_xxyz_0[i] * pa_y[i] - ta1_x_y_xxyz_1[i] * pc_y[i];

        ta1_x_yy_xxzz_0[i] = ta1_x_0_xxzz_0[i] * fe_0 - ta1_x_0_xxzz_1[i] * fe_0 + ta1_x_y_xxzz_0[i] * pa_y[i] - ta1_x_y_xxzz_1[i] * pc_y[i];

        ta1_x_yy_xyyy_0[i] = ta1_x_0_xyyy_0[i] * fe_0 - ta1_x_0_xyyy_1[i] * fe_0 + 3.0 * ta1_x_y_xyy_0[i] * fe_0 - 3.0 * ta1_x_y_xyy_1[i] * fe_0 +
                             ta1_x_y_xyyy_0[i] * pa_y[i] - ta1_x_y_xyyy_1[i] * pc_y[i];

        ta1_x_yy_xyyz_0[i] = ta1_x_0_xyyz_0[i] * fe_0 - ta1_x_0_xyyz_1[i] * fe_0 + 2.0 * ta1_x_y_xyz_0[i] * fe_0 - 2.0 * ta1_x_y_xyz_1[i] * fe_0 +
                             ta1_x_y_xyyz_0[i] * pa_y[i] - ta1_x_y_xyyz_1[i] * pc_y[i];

        ta1_x_yy_xyzz_0[i] = ta1_x_0_xyzz_0[i] * fe_0 - ta1_x_0_xyzz_1[i] * fe_0 + ta1_x_y_xzz_0[i] * fe_0 - ta1_x_y_xzz_1[i] * fe_0 +
                             ta1_x_y_xyzz_0[i] * pa_y[i] - ta1_x_y_xyzz_1[i] * pc_y[i];

        ta1_x_yy_xzzz_0[i] = ta1_x_0_xzzz_0[i] * fe_0 - ta1_x_0_xzzz_1[i] * fe_0 + ta1_x_y_xzzz_0[i] * pa_y[i] - ta1_x_y_xzzz_1[i] * pc_y[i];

        ta1_x_yy_yyyy_0[i] = ta1_x_0_yyyy_0[i] * fe_0 - ta1_x_0_yyyy_1[i] * fe_0 + 4.0 * ta1_x_y_yyy_0[i] * fe_0 - 4.0 * ta1_x_y_yyy_1[i] * fe_0 +
                             ta1_x_y_yyyy_0[i] * pa_y[i] - ta1_x_y_yyyy_1[i] * pc_y[i];

        ta1_x_yy_yyyz_0[i] = ta1_x_0_yyyz_0[i] * fe_0 - ta1_x_0_yyyz_1[i] * fe_0 + 3.0 * ta1_x_y_yyz_0[i] * fe_0 - 3.0 * ta1_x_y_yyz_1[i] * fe_0 +
                             ta1_x_y_yyyz_0[i] * pa_y[i] - ta1_x_y_yyyz_1[i] * pc_y[i];

        ta1_x_yy_yyzz_0[i] = ta1_x_0_yyzz_0[i] * fe_0 - ta1_x_0_yyzz_1[i] * fe_0 + 2.0 * ta1_x_y_yzz_0[i] * fe_0 - 2.0 * ta1_x_y_yzz_1[i] * fe_0 +
                             ta1_x_y_yyzz_0[i] * pa_y[i] - ta1_x_y_yyzz_1[i] * pc_y[i];

        ta1_x_yy_yzzz_0[i] = ta1_x_0_yzzz_0[i] * fe_0 - ta1_x_0_yzzz_1[i] * fe_0 + ta1_x_y_zzz_0[i] * fe_0 - ta1_x_y_zzz_1[i] * fe_0 +
                             ta1_x_y_yzzz_0[i] * pa_y[i] - ta1_x_y_yzzz_1[i] * pc_y[i];

        ta1_x_yy_zzzz_0[i] = ta1_x_0_zzzz_0[i] * fe_0 - ta1_x_0_zzzz_1[i] * fe_0 + ta1_x_y_zzzz_0[i] * pa_y[i] - ta1_x_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 60-75 components of targeted buffer : DG

    auto ta1_x_yz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 60);

    auto ta1_x_yz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 61);

    auto ta1_x_yz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 62);

    auto ta1_x_yz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 63);

    auto ta1_x_yz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 64);

    auto ta1_x_yz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 65);

    auto ta1_x_yz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 66);

    auto ta1_x_yz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 67);

    auto ta1_x_yz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 68);

    auto ta1_x_yz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 69);

    auto ta1_x_yz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 70);

    auto ta1_x_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 71);

    auto ta1_x_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 72);

    auto ta1_x_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 73);

    auto ta1_x_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 74);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_y_xxxy_0,  \
                             ta1_x_y_xxxy_1,  \
                             ta1_x_y_xxyy_0,  \
                             ta1_x_y_xxyy_1,  \
                             ta1_x_y_xyyy_0,  \
                             ta1_x_y_xyyy_1,  \
                             ta1_x_y_yyyy_0,  \
                             ta1_x_y_yyyy_1,  \
                             ta1_x_yz_xxxx_0, \
                             ta1_x_yz_xxxy_0, \
                             ta1_x_yz_xxxz_0, \
                             ta1_x_yz_xxyy_0, \
                             ta1_x_yz_xxyz_0, \
                             ta1_x_yz_xxzz_0, \
                             ta1_x_yz_xyyy_0, \
                             ta1_x_yz_xyyz_0, \
                             ta1_x_yz_xyzz_0, \
                             ta1_x_yz_xzzz_0, \
                             ta1_x_yz_yyyy_0, \
                             ta1_x_yz_yyyz_0, \
                             ta1_x_yz_yyzz_0, \
                             ta1_x_yz_yzzz_0, \
                             ta1_x_yz_zzzz_0, \
                             ta1_x_z_xxxx_0,  \
                             ta1_x_z_xxxx_1,  \
                             ta1_x_z_xxxz_0,  \
                             ta1_x_z_xxxz_1,  \
                             ta1_x_z_xxyz_0,  \
                             ta1_x_z_xxyz_1,  \
                             ta1_x_z_xxz_0,   \
                             ta1_x_z_xxz_1,   \
                             ta1_x_z_xxzz_0,  \
                             ta1_x_z_xxzz_1,  \
                             ta1_x_z_xyyz_0,  \
                             ta1_x_z_xyyz_1,  \
                             ta1_x_z_xyz_0,   \
                             ta1_x_z_xyz_1,   \
                             ta1_x_z_xyzz_0,  \
                             ta1_x_z_xyzz_1,  \
                             ta1_x_z_xzz_0,   \
                             ta1_x_z_xzz_1,   \
                             ta1_x_z_xzzz_0,  \
                             ta1_x_z_xzzz_1,  \
                             ta1_x_z_yyyz_0,  \
                             ta1_x_z_yyyz_1,  \
                             ta1_x_z_yyz_0,   \
                             ta1_x_z_yyz_1,   \
                             ta1_x_z_yyzz_0,  \
                             ta1_x_z_yyzz_1,  \
                             ta1_x_z_yzz_0,   \
                             ta1_x_z_yzz_1,   \
                             ta1_x_z_yzzz_0,  \
                             ta1_x_z_yzzz_1,  \
                             ta1_x_z_zzz_0,   \
                             ta1_x_z_zzz_1,   \
                             ta1_x_z_zzzz_0,  \
                             ta1_x_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yz_xxxx_0[i] = ta1_x_z_xxxx_0[i] * pa_y[i] - ta1_x_z_xxxx_1[i] * pc_y[i];

        ta1_x_yz_xxxy_0[i] = ta1_x_y_xxxy_0[i] * pa_z[i] - ta1_x_y_xxxy_1[i] * pc_z[i];

        ta1_x_yz_xxxz_0[i] = ta1_x_z_xxxz_0[i] * pa_y[i] - ta1_x_z_xxxz_1[i] * pc_y[i];

        ta1_x_yz_xxyy_0[i] = ta1_x_y_xxyy_0[i] * pa_z[i] - ta1_x_y_xxyy_1[i] * pc_z[i];

        ta1_x_yz_xxyz_0[i] = ta1_x_z_xxz_0[i] * fe_0 - ta1_x_z_xxz_1[i] * fe_0 + ta1_x_z_xxyz_0[i] * pa_y[i] - ta1_x_z_xxyz_1[i] * pc_y[i];

        ta1_x_yz_xxzz_0[i] = ta1_x_z_xxzz_0[i] * pa_y[i] - ta1_x_z_xxzz_1[i] * pc_y[i];

        ta1_x_yz_xyyy_0[i] = ta1_x_y_xyyy_0[i] * pa_z[i] - ta1_x_y_xyyy_1[i] * pc_z[i];

        ta1_x_yz_xyyz_0[i] =
            2.0 * ta1_x_z_xyz_0[i] * fe_0 - 2.0 * ta1_x_z_xyz_1[i] * fe_0 + ta1_x_z_xyyz_0[i] * pa_y[i] - ta1_x_z_xyyz_1[i] * pc_y[i];

        ta1_x_yz_xyzz_0[i] = ta1_x_z_xzz_0[i] * fe_0 - ta1_x_z_xzz_1[i] * fe_0 + ta1_x_z_xyzz_0[i] * pa_y[i] - ta1_x_z_xyzz_1[i] * pc_y[i];

        ta1_x_yz_xzzz_0[i] = ta1_x_z_xzzz_0[i] * pa_y[i] - ta1_x_z_xzzz_1[i] * pc_y[i];

        ta1_x_yz_yyyy_0[i] = ta1_x_y_yyyy_0[i] * pa_z[i] - ta1_x_y_yyyy_1[i] * pc_z[i];

        ta1_x_yz_yyyz_0[i] =
            3.0 * ta1_x_z_yyz_0[i] * fe_0 - 3.0 * ta1_x_z_yyz_1[i] * fe_0 + ta1_x_z_yyyz_0[i] * pa_y[i] - ta1_x_z_yyyz_1[i] * pc_y[i];

        ta1_x_yz_yyzz_0[i] =
            2.0 * ta1_x_z_yzz_0[i] * fe_0 - 2.0 * ta1_x_z_yzz_1[i] * fe_0 + ta1_x_z_yyzz_0[i] * pa_y[i] - ta1_x_z_yyzz_1[i] * pc_y[i];

        ta1_x_yz_yzzz_0[i] = ta1_x_z_zzz_0[i] * fe_0 - ta1_x_z_zzz_1[i] * fe_0 + ta1_x_z_yzzz_0[i] * pa_y[i] - ta1_x_z_yzzz_1[i] * pc_y[i];

        ta1_x_yz_zzzz_0[i] = ta1_x_z_zzzz_0[i] * pa_y[i] - ta1_x_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : DG

    auto ta1_x_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 75);

    auto ta1_x_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 76);

    auto ta1_x_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 77);

    auto ta1_x_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 78);

    auto ta1_x_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 79);

    auto ta1_x_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 80);

    auto ta1_x_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 81);

    auto ta1_x_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 82);

    auto ta1_x_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 83);

    auto ta1_x_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 84);

    auto ta1_x_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 85);

    auto ta1_x_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 86);

    auto ta1_x_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 87);

    auto ta1_x_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 88);

    auto ta1_x_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 89);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_0_xxxx_0,  \
                             ta1_x_0_xxxx_1,  \
                             ta1_x_0_xxxy_0,  \
                             ta1_x_0_xxxy_1,  \
                             ta1_x_0_xxxz_0,  \
                             ta1_x_0_xxxz_1,  \
                             ta1_x_0_xxyy_0,  \
                             ta1_x_0_xxyy_1,  \
                             ta1_x_0_xxyz_0,  \
                             ta1_x_0_xxyz_1,  \
                             ta1_x_0_xxzz_0,  \
                             ta1_x_0_xxzz_1,  \
                             ta1_x_0_xyyy_0,  \
                             ta1_x_0_xyyy_1,  \
                             ta1_x_0_xyyz_0,  \
                             ta1_x_0_xyyz_1,  \
                             ta1_x_0_xyzz_0,  \
                             ta1_x_0_xyzz_1,  \
                             ta1_x_0_xzzz_0,  \
                             ta1_x_0_xzzz_1,  \
                             ta1_x_0_yyyy_0,  \
                             ta1_x_0_yyyy_1,  \
                             ta1_x_0_yyyz_0,  \
                             ta1_x_0_yyyz_1,  \
                             ta1_x_0_yyzz_0,  \
                             ta1_x_0_yyzz_1,  \
                             ta1_x_0_yzzz_0,  \
                             ta1_x_0_yzzz_1,  \
                             ta1_x_0_zzzz_0,  \
                             ta1_x_0_zzzz_1,  \
                             ta1_x_z_xxx_0,   \
                             ta1_x_z_xxx_1,   \
                             ta1_x_z_xxxx_0,  \
                             ta1_x_z_xxxx_1,  \
                             ta1_x_z_xxxy_0,  \
                             ta1_x_z_xxxy_1,  \
                             ta1_x_z_xxxz_0,  \
                             ta1_x_z_xxxz_1,  \
                             ta1_x_z_xxy_0,   \
                             ta1_x_z_xxy_1,   \
                             ta1_x_z_xxyy_0,  \
                             ta1_x_z_xxyy_1,  \
                             ta1_x_z_xxyz_0,  \
                             ta1_x_z_xxyz_1,  \
                             ta1_x_z_xxz_0,   \
                             ta1_x_z_xxz_1,   \
                             ta1_x_z_xxzz_0,  \
                             ta1_x_z_xxzz_1,  \
                             ta1_x_z_xyy_0,   \
                             ta1_x_z_xyy_1,   \
                             ta1_x_z_xyyy_0,  \
                             ta1_x_z_xyyy_1,  \
                             ta1_x_z_xyyz_0,  \
                             ta1_x_z_xyyz_1,  \
                             ta1_x_z_xyz_0,   \
                             ta1_x_z_xyz_1,   \
                             ta1_x_z_xyzz_0,  \
                             ta1_x_z_xyzz_1,  \
                             ta1_x_z_xzz_0,   \
                             ta1_x_z_xzz_1,   \
                             ta1_x_z_xzzz_0,  \
                             ta1_x_z_xzzz_1,  \
                             ta1_x_z_yyy_0,   \
                             ta1_x_z_yyy_1,   \
                             ta1_x_z_yyyy_0,  \
                             ta1_x_z_yyyy_1,  \
                             ta1_x_z_yyyz_0,  \
                             ta1_x_z_yyyz_1,  \
                             ta1_x_z_yyz_0,   \
                             ta1_x_z_yyz_1,   \
                             ta1_x_z_yyzz_0,  \
                             ta1_x_z_yyzz_1,  \
                             ta1_x_z_yzz_0,   \
                             ta1_x_z_yzz_1,   \
                             ta1_x_z_yzzz_0,  \
                             ta1_x_z_yzzz_1,  \
                             ta1_x_z_zzz_0,   \
                             ta1_x_z_zzz_1,   \
                             ta1_x_z_zzzz_0,  \
                             ta1_x_z_zzzz_1,  \
                             ta1_x_zz_xxxx_0, \
                             ta1_x_zz_xxxy_0, \
                             ta1_x_zz_xxxz_0, \
                             ta1_x_zz_xxyy_0, \
                             ta1_x_zz_xxyz_0, \
                             ta1_x_zz_xxzz_0, \
                             ta1_x_zz_xyyy_0, \
                             ta1_x_zz_xyyz_0, \
                             ta1_x_zz_xyzz_0, \
                             ta1_x_zz_xzzz_0, \
                             ta1_x_zz_yyyy_0, \
                             ta1_x_zz_yyyz_0, \
                             ta1_x_zz_yyzz_0, \
                             ta1_x_zz_yzzz_0, \
                             ta1_x_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zz_xxxx_0[i] = ta1_x_0_xxxx_0[i] * fe_0 - ta1_x_0_xxxx_1[i] * fe_0 + ta1_x_z_xxxx_0[i] * pa_z[i] - ta1_x_z_xxxx_1[i] * pc_z[i];

        ta1_x_zz_xxxy_0[i] = ta1_x_0_xxxy_0[i] * fe_0 - ta1_x_0_xxxy_1[i] * fe_0 + ta1_x_z_xxxy_0[i] * pa_z[i] - ta1_x_z_xxxy_1[i] * pc_z[i];

        ta1_x_zz_xxxz_0[i] = ta1_x_0_xxxz_0[i] * fe_0 - ta1_x_0_xxxz_1[i] * fe_0 + ta1_x_z_xxx_0[i] * fe_0 - ta1_x_z_xxx_1[i] * fe_0 +
                             ta1_x_z_xxxz_0[i] * pa_z[i] - ta1_x_z_xxxz_1[i] * pc_z[i];

        ta1_x_zz_xxyy_0[i] = ta1_x_0_xxyy_0[i] * fe_0 - ta1_x_0_xxyy_1[i] * fe_0 + ta1_x_z_xxyy_0[i] * pa_z[i] - ta1_x_z_xxyy_1[i] * pc_z[i];

        ta1_x_zz_xxyz_0[i] = ta1_x_0_xxyz_0[i] * fe_0 - ta1_x_0_xxyz_1[i] * fe_0 + ta1_x_z_xxy_0[i] * fe_0 - ta1_x_z_xxy_1[i] * fe_0 +
                             ta1_x_z_xxyz_0[i] * pa_z[i] - ta1_x_z_xxyz_1[i] * pc_z[i];

        ta1_x_zz_xxzz_0[i] = ta1_x_0_xxzz_0[i] * fe_0 - ta1_x_0_xxzz_1[i] * fe_0 + 2.0 * ta1_x_z_xxz_0[i] * fe_0 - 2.0 * ta1_x_z_xxz_1[i] * fe_0 +
                             ta1_x_z_xxzz_0[i] * pa_z[i] - ta1_x_z_xxzz_1[i] * pc_z[i];

        ta1_x_zz_xyyy_0[i] = ta1_x_0_xyyy_0[i] * fe_0 - ta1_x_0_xyyy_1[i] * fe_0 + ta1_x_z_xyyy_0[i] * pa_z[i] - ta1_x_z_xyyy_1[i] * pc_z[i];

        ta1_x_zz_xyyz_0[i] = ta1_x_0_xyyz_0[i] * fe_0 - ta1_x_0_xyyz_1[i] * fe_0 + ta1_x_z_xyy_0[i] * fe_0 - ta1_x_z_xyy_1[i] * fe_0 +
                             ta1_x_z_xyyz_0[i] * pa_z[i] - ta1_x_z_xyyz_1[i] * pc_z[i];

        ta1_x_zz_xyzz_0[i] = ta1_x_0_xyzz_0[i] * fe_0 - ta1_x_0_xyzz_1[i] * fe_0 + 2.0 * ta1_x_z_xyz_0[i] * fe_0 - 2.0 * ta1_x_z_xyz_1[i] * fe_0 +
                             ta1_x_z_xyzz_0[i] * pa_z[i] - ta1_x_z_xyzz_1[i] * pc_z[i];

        ta1_x_zz_xzzz_0[i] = ta1_x_0_xzzz_0[i] * fe_0 - ta1_x_0_xzzz_1[i] * fe_0 + 3.0 * ta1_x_z_xzz_0[i] * fe_0 - 3.0 * ta1_x_z_xzz_1[i] * fe_0 +
                             ta1_x_z_xzzz_0[i] * pa_z[i] - ta1_x_z_xzzz_1[i] * pc_z[i];

        ta1_x_zz_yyyy_0[i] = ta1_x_0_yyyy_0[i] * fe_0 - ta1_x_0_yyyy_1[i] * fe_0 + ta1_x_z_yyyy_0[i] * pa_z[i] - ta1_x_z_yyyy_1[i] * pc_z[i];

        ta1_x_zz_yyyz_0[i] = ta1_x_0_yyyz_0[i] * fe_0 - ta1_x_0_yyyz_1[i] * fe_0 + ta1_x_z_yyy_0[i] * fe_0 - ta1_x_z_yyy_1[i] * fe_0 +
                             ta1_x_z_yyyz_0[i] * pa_z[i] - ta1_x_z_yyyz_1[i] * pc_z[i];

        ta1_x_zz_yyzz_0[i] = ta1_x_0_yyzz_0[i] * fe_0 - ta1_x_0_yyzz_1[i] * fe_0 + 2.0 * ta1_x_z_yyz_0[i] * fe_0 - 2.0 * ta1_x_z_yyz_1[i] * fe_0 +
                             ta1_x_z_yyzz_0[i] * pa_z[i] - ta1_x_z_yyzz_1[i] * pc_z[i];

        ta1_x_zz_yzzz_0[i] = ta1_x_0_yzzz_0[i] * fe_0 - ta1_x_0_yzzz_1[i] * fe_0 + 3.0 * ta1_x_z_yzz_0[i] * fe_0 - 3.0 * ta1_x_z_yzz_1[i] * fe_0 +
                             ta1_x_z_yzzz_0[i] * pa_z[i] - ta1_x_z_yzzz_1[i] * pc_z[i];

        ta1_x_zz_zzzz_0[i] = ta1_x_0_zzzz_0[i] * fe_0 - ta1_x_0_zzzz_1[i] * fe_0 + 4.0 * ta1_x_z_zzz_0[i] * fe_0 - 4.0 * ta1_x_z_zzz_1[i] * fe_0 +
                             ta1_x_z_zzzz_0[i] * pa_z[i] - ta1_x_z_zzzz_1[i] * pc_z[i];
    }

    // Set up 90-105 components of targeted buffer : DG

    auto ta1_y_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 90);

    auto ta1_y_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 91);

    auto ta1_y_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 92);

    auto ta1_y_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 93);

    auto ta1_y_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 94);

    auto ta1_y_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 95);

    auto ta1_y_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 96);

    auto ta1_y_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 97);

    auto ta1_y_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 98);

    auto ta1_y_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 99);

    auto ta1_y_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 100);

    auto ta1_y_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 101);

    auto ta1_y_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 102);

    auto ta1_y_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 103);

    auto ta1_y_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 104);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_0_xxxx_0,  \
                             ta1_y_0_xxxx_1,  \
                             ta1_y_0_xxxy_0,  \
                             ta1_y_0_xxxy_1,  \
                             ta1_y_0_xxxz_0,  \
                             ta1_y_0_xxxz_1,  \
                             ta1_y_0_xxyy_0,  \
                             ta1_y_0_xxyy_1,  \
                             ta1_y_0_xxyz_0,  \
                             ta1_y_0_xxyz_1,  \
                             ta1_y_0_xxzz_0,  \
                             ta1_y_0_xxzz_1,  \
                             ta1_y_0_xyyy_0,  \
                             ta1_y_0_xyyy_1,  \
                             ta1_y_0_xyyz_0,  \
                             ta1_y_0_xyyz_1,  \
                             ta1_y_0_xyzz_0,  \
                             ta1_y_0_xyzz_1,  \
                             ta1_y_0_xzzz_0,  \
                             ta1_y_0_xzzz_1,  \
                             ta1_y_0_yyyy_0,  \
                             ta1_y_0_yyyy_1,  \
                             ta1_y_0_yyyz_0,  \
                             ta1_y_0_yyyz_1,  \
                             ta1_y_0_yyzz_0,  \
                             ta1_y_0_yyzz_1,  \
                             ta1_y_0_yzzz_0,  \
                             ta1_y_0_yzzz_1,  \
                             ta1_y_0_zzzz_0,  \
                             ta1_y_0_zzzz_1,  \
                             ta1_y_x_xxx_0,   \
                             ta1_y_x_xxx_1,   \
                             ta1_y_x_xxxx_0,  \
                             ta1_y_x_xxxx_1,  \
                             ta1_y_x_xxxy_0,  \
                             ta1_y_x_xxxy_1,  \
                             ta1_y_x_xxxz_0,  \
                             ta1_y_x_xxxz_1,  \
                             ta1_y_x_xxy_0,   \
                             ta1_y_x_xxy_1,   \
                             ta1_y_x_xxyy_0,  \
                             ta1_y_x_xxyy_1,  \
                             ta1_y_x_xxyz_0,  \
                             ta1_y_x_xxyz_1,  \
                             ta1_y_x_xxz_0,   \
                             ta1_y_x_xxz_1,   \
                             ta1_y_x_xxzz_0,  \
                             ta1_y_x_xxzz_1,  \
                             ta1_y_x_xyy_0,   \
                             ta1_y_x_xyy_1,   \
                             ta1_y_x_xyyy_0,  \
                             ta1_y_x_xyyy_1,  \
                             ta1_y_x_xyyz_0,  \
                             ta1_y_x_xyyz_1,  \
                             ta1_y_x_xyz_0,   \
                             ta1_y_x_xyz_1,   \
                             ta1_y_x_xyzz_0,  \
                             ta1_y_x_xyzz_1,  \
                             ta1_y_x_xzz_0,   \
                             ta1_y_x_xzz_1,   \
                             ta1_y_x_xzzz_0,  \
                             ta1_y_x_xzzz_1,  \
                             ta1_y_x_yyy_0,   \
                             ta1_y_x_yyy_1,   \
                             ta1_y_x_yyyy_0,  \
                             ta1_y_x_yyyy_1,  \
                             ta1_y_x_yyyz_0,  \
                             ta1_y_x_yyyz_1,  \
                             ta1_y_x_yyz_0,   \
                             ta1_y_x_yyz_1,   \
                             ta1_y_x_yyzz_0,  \
                             ta1_y_x_yyzz_1,  \
                             ta1_y_x_yzz_0,   \
                             ta1_y_x_yzz_1,   \
                             ta1_y_x_yzzz_0,  \
                             ta1_y_x_yzzz_1,  \
                             ta1_y_x_zzz_0,   \
                             ta1_y_x_zzz_1,   \
                             ta1_y_x_zzzz_0,  \
                             ta1_y_x_zzzz_1,  \
                             ta1_y_xx_xxxx_0, \
                             ta1_y_xx_xxxy_0, \
                             ta1_y_xx_xxxz_0, \
                             ta1_y_xx_xxyy_0, \
                             ta1_y_xx_xxyz_0, \
                             ta1_y_xx_xxzz_0, \
                             ta1_y_xx_xyyy_0, \
                             ta1_y_xx_xyyz_0, \
                             ta1_y_xx_xyzz_0, \
                             ta1_y_xx_xzzz_0, \
                             ta1_y_xx_yyyy_0, \
                             ta1_y_xx_yyyz_0, \
                             ta1_y_xx_yyzz_0, \
                             ta1_y_xx_yzzz_0, \
                             ta1_y_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xx_xxxx_0[i] = ta1_y_0_xxxx_0[i] * fe_0 - ta1_y_0_xxxx_1[i] * fe_0 + 4.0 * ta1_y_x_xxx_0[i] * fe_0 - 4.0 * ta1_y_x_xxx_1[i] * fe_0 +
                             ta1_y_x_xxxx_0[i] * pa_x[i] - ta1_y_x_xxxx_1[i] * pc_x[i];

        ta1_y_xx_xxxy_0[i] = ta1_y_0_xxxy_0[i] * fe_0 - ta1_y_0_xxxy_1[i] * fe_0 + 3.0 * ta1_y_x_xxy_0[i] * fe_0 - 3.0 * ta1_y_x_xxy_1[i] * fe_0 +
                             ta1_y_x_xxxy_0[i] * pa_x[i] - ta1_y_x_xxxy_1[i] * pc_x[i];

        ta1_y_xx_xxxz_0[i] = ta1_y_0_xxxz_0[i] * fe_0 - ta1_y_0_xxxz_1[i] * fe_0 + 3.0 * ta1_y_x_xxz_0[i] * fe_0 - 3.0 * ta1_y_x_xxz_1[i] * fe_0 +
                             ta1_y_x_xxxz_0[i] * pa_x[i] - ta1_y_x_xxxz_1[i] * pc_x[i];

        ta1_y_xx_xxyy_0[i] = ta1_y_0_xxyy_0[i] * fe_0 - ta1_y_0_xxyy_1[i] * fe_0 + 2.0 * ta1_y_x_xyy_0[i] * fe_0 - 2.0 * ta1_y_x_xyy_1[i] * fe_0 +
                             ta1_y_x_xxyy_0[i] * pa_x[i] - ta1_y_x_xxyy_1[i] * pc_x[i];

        ta1_y_xx_xxyz_0[i] = ta1_y_0_xxyz_0[i] * fe_0 - ta1_y_0_xxyz_1[i] * fe_0 + 2.0 * ta1_y_x_xyz_0[i] * fe_0 - 2.0 * ta1_y_x_xyz_1[i] * fe_0 +
                             ta1_y_x_xxyz_0[i] * pa_x[i] - ta1_y_x_xxyz_1[i] * pc_x[i];

        ta1_y_xx_xxzz_0[i] = ta1_y_0_xxzz_0[i] * fe_0 - ta1_y_0_xxzz_1[i] * fe_0 + 2.0 * ta1_y_x_xzz_0[i] * fe_0 - 2.0 * ta1_y_x_xzz_1[i] * fe_0 +
                             ta1_y_x_xxzz_0[i] * pa_x[i] - ta1_y_x_xxzz_1[i] * pc_x[i];

        ta1_y_xx_xyyy_0[i] = ta1_y_0_xyyy_0[i] * fe_0 - ta1_y_0_xyyy_1[i] * fe_0 + ta1_y_x_yyy_0[i] * fe_0 - ta1_y_x_yyy_1[i] * fe_0 +
                             ta1_y_x_xyyy_0[i] * pa_x[i] - ta1_y_x_xyyy_1[i] * pc_x[i];

        ta1_y_xx_xyyz_0[i] = ta1_y_0_xyyz_0[i] * fe_0 - ta1_y_0_xyyz_1[i] * fe_0 + ta1_y_x_yyz_0[i] * fe_0 - ta1_y_x_yyz_1[i] * fe_0 +
                             ta1_y_x_xyyz_0[i] * pa_x[i] - ta1_y_x_xyyz_1[i] * pc_x[i];

        ta1_y_xx_xyzz_0[i] = ta1_y_0_xyzz_0[i] * fe_0 - ta1_y_0_xyzz_1[i] * fe_0 + ta1_y_x_yzz_0[i] * fe_0 - ta1_y_x_yzz_1[i] * fe_0 +
                             ta1_y_x_xyzz_0[i] * pa_x[i] - ta1_y_x_xyzz_1[i] * pc_x[i];

        ta1_y_xx_xzzz_0[i] = ta1_y_0_xzzz_0[i] * fe_0 - ta1_y_0_xzzz_1[i] * fe_0 + ta1_y_x_zzz_0[i] * fe_0 - ta1_y_x_zzz_1[i] * fe_0 +
                             ta1_y_x_xzzz_0[i] * pa_x[i] - ta1_y_x_xzzz_1[i] * pc_x[i];

        ta1_y_xx_yyyy_0[i] = ta1_y_0_yyyy_0[i] * fe_0 - ta1_y_0_yyyy_1[i] * fe_0 + ta1_y_x_yyyy_0[i] * pa_x[i] - ta1_y_x_yyyy_1[i] * pc_x[i];

        ta1_y_xx_yyyz_0[i] = ta1_y_0_yyyz_0[i] * fe_0 - ta1_y_0_yyyz_1[i] * fe_0 + ta1_y_x_yyyz_0[i] * pa_x[i] - ta1_y_x_yyyz_1[i] * pc_x[i];

        ta1_y_xx_yyzz_0[i] = ta1_y_0_yyzz_0[i] * fe_0 - ta1_y_0_yyzz_1[i] * fe_0 + ta1_y_x_yyzz_0[i] * pa_x[i] - ta1_y_x_yyzz_1[i] * pc_x[i];

        ta1_y_xx_yzzz_0[i] = ta1_y_0_yzzz_0[i] * fe_0 - ta1_y_0_yzzz_1[i] * fe_0 + ta1_y_x_yzzz_0[i] * pa_x[i] - ta1_y_x_yzzz_1[i] * pc_x[i];

        ta1_y_xx_zzzz_0[i] = ta1_y_0_zzzz_0[i] * fe_0 - ta1_y_0_zzzz_1[i] * fe_0 + ta1_y_x_zzzz_0[i] * pa_x[i] - ta1_y_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : DG

    auto ta1_y_xy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 105);

    auto ta1_y_xy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 106);

    auto ta1_y_xy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 107);

    auto ta1_y_xy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 108);

    auto ta1_y_xy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 109);

    auto ta1_y_xy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 110);

    auto ta1_y_xy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 111);

    auto ta1_y_xy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 112);

    auto ta1_y_xy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 113);

    auto ta1_y_xy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 114);

    auto ta1_y_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 115);

    auto ta1_y_xy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 116);

    auto ta1_y_xy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 117);

    auto ta1_y_xy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 118);

    auto ta1_y_xy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 119);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_x_xxxx_0,  \
                             ta1_y_x_xxxx_1,  \
                             ta1_y_x_xxxz_0,  \
                             ta1_y_x_xxxz_1,  \
                             ta1_y_x_xxzz_0,  \
                             ta1_y_x_xxzz_1,  \
                             ta1_y_x_xzzz_0,  \
                             ta1_y_x_xzzz_1,  \
                             ta1_y_xy_xxxx_0, \
                             ta1_y_xy_xxxy_0, \
                             ta1_y_xy_xxxz_0, \
                             ta1_y_xy_xxyy_0, \
                             ta1_y_xy_xxyz_0, \
                             ta1_y_xy_xxzz_0, \
                             ta1_y_xy_xyyy_0, \
                             ta1_y_xy_xyyz_0, \
                             ta1_y_xy_xyzz_0, \
                             ta1_y_xy_xzzz_0, \
                             ta1_y_xy_yyyy_0, \
                             ta1_y_xy_yyyz_0, \
                             ta1_y_xy_yyzz_0, \
                             ta1_y_xy_yzzz_0, \
                             ta1_y_xy_zzzz_0, \
                             ta1_y_y_xxxy_0,  \
                             ta1_y_y_xxxy_1,  \
                             ta1_y_y_xxy_0,   \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xxyy_0,  \
                             ta1_y_y_xxyy_1,  \
                             ta1_y_y_xxyz_0,  \
                             ta1_y_y_xxyz_1,  \
                             ta1_y_y_xyy_0,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_xyyy_0,  \
                             ta1_y_y_xyyy_1,  \
                             ta1_y_y_xyyz_0,  \
                             ta1_y_y_xyyz_1,  \
                             ta1_y_y_xyz_0,   \
                             ta1_y_y_xyz_1,   \
                             ta1_y_y_xyzz_0,  \
                             ta1_y_y_xyzz_1,  \
                             ta1_y_y_yyy_0,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_y_yyyy_0,  \
                             ta1_y_y_yyyy_1,  \
                             ta1_y_y_yyyz_0,  \
                             ta1_y_y_yyyz_1,  \
                             ta1_y_y_yyz_0,   \
                             ta1_y_y_yyz_1,   \
                             ta1_y_y_yyzz_0,  \
                             ta1_y_y_yyzz_1,  \
                             ta1_y_y_yzz_0,   \
                             ta1_y_y_yzz_1,   \
                             ta1_y_y_yzzz_0,  \
                             ta1_y_y_yzzz_1,  \
                             ta1_y_y_zzzz_0,  \
                             ta1_y_y_zzzz_1,  \
                             ta_x_xxxx_1,     \
                             ta_x_xxxz_1,     \
                             ta_x_xxzz_1,     \
                             ta_x_xzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xy_xxxx_0[i] = ta_x_xxxx_1[i] + ta1_y_x_xxxx_0[i] * pa_y[i] - ta1_y_x_xxxx_1[i] * pc_y[i];

        ta1_y_xy_xxxy_0[i] =
            3.0 * ta1_y_y_xxy_0[i] * fe_0 - 3.0 * ta1_y_y_xxy_1[i] * fe_0 + ta1_y_y_xxxy_0[i] * pa_x[i] - ta1_y_y_xxxy_1[i] * pc_x[i];

        ta1_y_xy_xxxz_0[i] = ta_x_xxxz_1[i] + ta1_y_x_xxxz_0[i] * pa_y[i] - ta1_y_x_xxxz_1[i] * pc_y[i];

        ta1_y_xy_xxyy_0[i] =
            2.0 * ta1_y_y_xyy_0[i] * fe_0 - 2.0 * ta1_y_y_xyy_1[i] * fe_0 + ta1_y_y_xxyy_0[i] * pa_x[i] - ta1_y_y_xxyy_1[i] * pc_x[i];

        ta1_y_xy_xxyz_0[i] =
            2.0 * ta1_y_y_xyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyz_1[i] * fe_0 + ta1_y_y_xxyz_0[i] * pa_x[i] - ta1_y_y_xxyz_1[i] * pc_x[i];

        ta1_y_xy_xxzz_0[i] = ta_x_xxzz_1[i] + ta1_y_x_xxzz_0[i] * pa_y[i] - ta1_y_x_xxzz_1[i] * pc_y[i];

        ta1_y_xy_xyyy_0[i] = ta1_y_y_yyy_0[i] * fe_0 - ta1_y_y_yyy_1[i] * fe_0 + ta1_y_y_xyyy_0[i] * pa_x[i] - ta1_y_y_xyyy_1[i] * pc_x[i];

        ta1_y_xy_xyyz_0[i] = ta1_y_y_yyz_0[i] * fe_0 - ta1_y_y_yyz_1[i] * fe_0 + ta1_y_y_xyyz_0[i] * pa_x[i] - ta1_y_y_xyyz_1[i] * pc_x[i];

        ta1_y_xy_xyzz_0[i] = ta1_y_y_yzz_0[i] * fe_0 - ta1_y_y_yzz_1[i] * fe_0 + ta1_y_y_xyzz_0[i] * pa_x[i] - ta1_y_y_xyzz_1[i] * pc_x[i];

        ta1_y_xy_xzzz_0[i] = ta_x_xzzz_1[i] + ta1_y_x_xzzz_0[i] * pa_y[i] - ta1_y_x_xzzz_1[i] * pc_y[i];

        ta1_y_xy_yyyy_0[i] = ta1_y_y_yyyy_0[i] * pa_x[i] - ta1_y_y_yyyy_1[i] * pc_x[i];

        ta1_y_xy_yyyz_0[i] = ta1_y_y_yyyz_0[i] * pa_x[i] - ta1_y_y_yyyz_1[i] * pc_x[i];

        ta1_y_xy_yyzz_0[i] = ta1_y_y_yyzz_0[i] * pa_x[i] - ta1_y_y_yyzz_1[i] * pc_x[i];

        ta1_y_xy_yzzz_0[i] = ta1_y_y_yzzz_0[i] * pa_x[i] - ta1_y_y_yzzz_1[i] * pc_x[i];

        ta1_y_xy_zzzz_0[i] = ta1_y_y_zzzz_0[i] * pa_x[i] - ta1_y_y_zzzz_1[i] * pc_x[i];
    }

    // Set up 120-135 components of targeted buffer : DG

    auto ta1_y_xz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 120);

    auto ta1_y_xz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 121);

    auto ta1_y_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 122);

    auto ta1_y_xz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 123);

    auto ta1_y_xz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 124);

    auto ta1_y_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 125);

    auto ta1_y_xz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 126);

    auto ta1_y_xz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 127);

    auto ta1_y_xz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 128);

    auto ta1_y_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 129);

    auto ta1_y_xz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 130);

    auto ta1_y_xz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 131);

    auto ta1_y_xz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 132);

    auto ta1_y_xz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 133);

    auto ta1_y_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 134);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_x_xxxx_0,  \
                             ta1_y_x_xxxx_1,  \
                             ta1_y_x_xxxy_0,  \
                             ta1_y_x_xxxy_1,  \
                             ta1_y_x_xxyy_0,  \
                             ta1_y_x_xxyy_1,  \
                             ta1_y_x_xyyy_0,  \
                             ta1_y_x_xyyy_1,  \
                             ta1_y_xz_xxxx_0, \
                             ta1_y_xz_xxxy_0, \
                             ta1_y_xz_xxxz_0, \
                             ta1_y_xz_xxyy_0, \
                             ta1_y_xz_xxyz_0, \
                             ta1_y_xz_xxzz_0, \
                             ta1_y_xz_xyyy_0, \
                             ta1_y_xz_xyyz_0, \
                             ta1_y_xz_xyzz_0, \
                             ta1_y_xz_xzzz_0, \
                             ta1_y_xz_yyyy_0, \
                             ta1_y_xz_yyyz_0, \
                             ta1_y_xz_yyzz_0, \
                             ta1_y_xz_yzzz_0, \
                             ta1_y_xz_zzzz_0, \
                             ta1_y_z_xxxz_0,  \
                             ta1_y_z_xxxz_1,  \
                             ta1_y_z_xxyz_0,  \
                             ta1_y_z_xxyz_1,  \
                             ta1_y_z_xxz_0,   \
                             ta1_y_z_xxz_1,   \
                             ta1_y_z_xxzz_0,  \
                             ta1_y_z_xxzz_1,  \
                             ta1_y_z_xyyz_0,  \
                             ta1_y_z_xyyz_1,  \
                             ta1_y_z_xyz_0,   \
                             ta1_y_z_xyz_1,   \
                             ta1_y_z_xyzz_0,  \
                             ta1_y_z_xyzz_1,  \
                             ta1_y_z_xzz_0,   \
                             ta1_y_z_xzz_1,   \
                             ta1_y_z_xzzz_0,  \
                             ta1_y_z_xzzz_1,  \
                             ta1_y_z_yyyy_0,  \
                             ta1_y_z_yyyy_1,  \
                             ta1_y_z_yyyz_0,  \
                             ta1_y_z_yyyz_1,  \
                             ta1_y_z_yyz_0,   \
                             ta1_y_z_yyz_1,   \
                             ta1_y_z_yyzz_0,  \
                             ta1_y_z_yyzz_1,  \
                             ta1_y_z_yzz_0,   \
                             ta1_y_z_yzz_1,   \
                             ta1_y_z_yzzz_0,  \
                             ta1_y_z_yzzz_1,  \
                             ta1_y_z_zzz_0,   \
                             ta1_y_z_zzz_1,   \
                             ta1_y_z_zzzz_0,  \
                             ta1_y_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xz_xxxx_0[i] = ta1_y_x_xxxx_0[i] * pa_z[i] - ta1_y_x_xxxx_1[i] * pc_z[i];

        ta1_y_xz_xxxy_0[i] = ta1_y_x_xxxy_0[i] * pa_z[i] - ta1_y_x_xxxy_1[i] * pc_z[i];

        ta1_y_xz_xxxz_0[i] =
            3.0 * ta1_y_z_xxz_0[i] * fe_0 - 3.0 * ta1_y_z_xxz_1[i] * fe_0 + ta1_y_z_xxxz_0[i] * pa_x[i] - ta1_y_z_xxxz_1[i] * pc_x[i];

        ta1_y_xz_xxyy_0[i] = ta1_y_x_xxyy_0[i] * pa_z[i] - ta1_y_x_xxyy_1[i] * pc_z[i];

        ta1_y_xz_xxyz_0[i] =
            2.0 * ta1_y_z_xyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyz_1[i] * fe_0 + ta1_y_z_xxyz_0[i] * pa_x[i] - ta1_y_z_xxyz_1[i] * pc_x[i];

        ta1_y_xz_xxzz_0[i] =
            2.0 * ta1_y_z_xzz_0[i] * fe_0 - 2.0 * ta1_y_z_xzz_1[i] * fe_0 + ta1_y_z_xxzz_0[i] * pa_x[i] - ta1_y_z_xxzz_1[i] * pc_x[i];

        ta1_y_xz_xyyy_0[i] = ta1_y_x_xyyy_0[i] * pa_z[i] - ta1_y_x_xyyy_1[i] * pc_z[i];

        ta1_y_xz_xyyz_0[i] = ta1_y_z_yyz_0[i] * fe_0 - ta1_y_z_yyz_1[i] * fe_0 + ta1_y_z_xyyz_0[i] * pa_x[i] - ta1_y_z_xyyz_1[i] * pc_x[i];

        ta1_y_xz_xyzz_0[i] = ta1_y_z_yzz_0[i] * fe_0 - ta1_y_z_yzz_1[i] * fe_0 + ta1_y_z_xyzz_0[i] * pa_x[i] - ta1_y_z_xyzz_1[i] * pc_x[i];

        ta1_y_xz_xzzz_0[i] = ta1_y_z_zzz_0[i] * fe_0 - ta1_y_z_zzz_1[i] * fe_0 + ta1_y_z_xzzz_0[i] * pa_x[i] - ta1_y_z_xzzz_1[i] * pc_x[i];

        ta1_y_xz_yyyy_0[i] = ta1_y_z_yyyy_0[i] * pa_x[i] - ta1_y_z_yyyy_1[i] * pc_x[i];

        ta1_y_xz_yyyz_0[i] = ta1_y_z_yyyz_0[i] * pa_x[i] - ta1_y_z_yyyz_1[i] * pc_x[i];

        ta1_y_xz_yyzz_0[i] = ta1_y_z_yyzz_0[i] * pa_x[i] - ta1_y_z_yyzz_1[i] * pc_x[i];

        ta1_y_xz_yzzz_0[i] = ta1_y_z_yzzz_0[i] * pa_x[i] - ta1_y_z_yzzz_1[i] * pc_x[i];

        ta1_y_xz_zzzz_0[i] = ta1_y_z_zzzz_0[i] * pa_x[i] - ta1_y_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 135-150 components of targeted buffer : DG

    auto ta1_y_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 135);

    auto ta1_y_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 136);

    auto ta1_y_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 137);

    auto ta1_y_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 138);

    auto ta1_y_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 139);

    auto ta1_y_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 140);

    auto ta1_y_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 141);

    auto ta1_y_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 142);

    auto ta1_y_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 143);

    auto ta1_y_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 144);

    auto ta1_y_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 145);

    auto ta1_y_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 146);

    auto ta1_y_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 147);

    auto ta1_y_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 148);

    auto ta1_y_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 149);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_y_0_xxxx_0,  \
                             ta1_y_0_xxxx_1,  \
                             ta1_y_0_xxxy_0,  \
                             ta1_y_0_xxxy_1,  \
                             ta1_y_0_xxxz_0,  \
                             ta1_y_0_xxxz_1,  \
                             ta1_y_0_xxyy_0,  \
                             ta1_y_0_xxyy_1,  \
                             ta1_y_0_xxyz_0,  \
                             ta1_y_0_xxyz_1,  \
                             ta1_y_0_xxzz_0,  \
                             ta1_y_0_xxzz_1,  \
                             ta1_y_0_xyyy_0,  \
                             ta1_y_0_xyyy_1,  \
                             ta1_y_0_xyyz_0,  \
                             ta1_y_0_xyyz_1,  \
                             ta1_y_0_xyzz_0,  \
                             ta1_y_0_xyzz_1,  \
                             ta1_y_0_xzzz_0,  \
                             ta1_y_0_xzzz_1,  \
                             ta1_y_0_yyyy_0,  \
                             ta1_y_0_yyyy_1,  \
                             ta1_y_0_yyyz_0,  \
                             ta1_y_0_yyyz_1,  \
                             ta1_y_0_yyzz_0,  \
                             ta1_y_0_yyzz_1,  \
                             ta1_y_0_yzzz_0,  \
                             ta1_y_0_yzzz_1,  \
                             ta1_y_0_zzzz_0,  \
                             ta1_y_0_zzzz_1,  \
                             ta1_y_y_xxx_0,   \
                             ta1_y_y_xxx_1,   \
                             ta1_y_y_xxxx_0,  \
                             ta1_y_y_xxxx_1,  \
                             ta1_y_y_xxxy_0,  \
                             ta1_y_y_xxxy_1,  \
                             ta1_y_y_xxxz_0,  \
                             ta1_y_y_xxxz_1,  \
                             ta1_y_y_xxy_0,   \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xxyy_0,  \
                             ta1_y_y_xxyy_1,  \
                             ta1_y_y_xxyz_0,  \
                             ta1_y_y_xxyz_1,  \
                             ta1_y_y_xxz_0,   \
                             ta1_y_y_xxz_1,   \
                             ta1_y_y_xxzz_0,  \
                             ta1_y_y_xxzz_1,  \
                             ta1_y_y_xyy_0,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_xyyy_0,  \
                             ta1_y_y_xyyy_1,  \
                             ta1_y_y_xyyz_0,  \
                             ta1_y_y_xyyz_1,  \
                             ta1_y_y_xyz_0,   \
                             ta1_y_y_xyz_1,   \
                             ta1_y_y_xyzz_0,  \
                             ta1_y_y_xyzz_1,  \
                             ta1_y_y_xzz_0,   \
                             ta1_y_y_xzz_1,   \
                             ta1_y_y_xzzz_0,  \
                             ta1_y_y_xzzz_1,  \
                             ta1_y_y_yyy_0,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_y_yyyy_0,  \
                             ta1_y_y_yyyy_1,  \
                             ta1_y_y_yyyz_0,  \
                             ta1_y_y_yyyz_1,  \
                             ta1_y_y_yyz_0,   \
                             ta1_y_y_yyz_1,   \
                             ta1_y_y_yyzz_0,  \
                             ta1_y_y_yyzz_1,  \
                             ta1_y_y_yzz_0,   \
                             ta1_y_y_yzz_1,   \
                             ta1_y_y_yzzz_0,  \
                             ta1_y_y_yzzz_1,  \
                             ta1_y_y_zzz_0,   \
                             ta1_y_y_zzz_1,   \
                             ta1_y_y_zzzz_0,  \
                             ta1_y_y_zzzz_1,  \
                             ta1_y_yy_xxxx_0, \
                             ta1_y_yy_xxxy_0, \
                             ta1_y_yy_xxxz_0, \
                             ta1_y_yy_xxyy_0, \
                             ta1_y_yy_xxyz_0, \
                             ta1_y_yy_xxzz_0, \
                             ta1_y_yy_xyyy_0, \
                             ta1_y_yy_xyyz_0, \
                             ta1_y_yy_xyzz_0, \
                             ta1_y_yy_xzzz_0, \
                             ta1_y_yy_yyyy_0, \
                             ta1_y_yy_yyyz_0, \
                             ta1_y_yy_yyzz_0, \
                             ta1_y_yy_yzzz_0, \
                             ta1_y_yy_zzzz_0, \
                             ta_y_xxxx_1,     \
                             ta_y_xxxy_1,     \
                             ta_y_xxxz_1,     \
                             ta_y_xxyy_1,     \
                             ta_y_xxyz_1,     \
                             ta_y_xxzz_1,     \
                             ta_y_xyyy_1,     \
                             ta_y_xyyz_1,     \
                             ta_y_xyzz_1,     \
                             ta_y_xzzz_1,     \
                             ta_y_yyyy_1,     \
                             ta_y_yyyz_1,     \
                             ta_y_yyzz_1,     \
                             ta_y_yzzz_1,     \
                             ta_y_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yy_xxxx_0[i] =
            ta1_y_0_xxxx_0[i] * fe_0 - ta1_y_0_xxxx_1[i] * fe_0 + ta_y_xxxx_1[i] + ta1_y_y_xxxx_0[i] * pa_y[i] - ta1_y_y_xxxx_1[i] * pc_y[i];

        ta1_y_yy_xxxy_0[i] = ta1_y_0_xxxy_0[i] * fe_0 - ta1_y_0_xxxy_1[i] * fe_0 + ta1_y_y_xxx_0[i] * fe_0 - ta1_y_y_xxx_1[i] * fe_0 +
                             ta_y_xxxy_1[i] + ta1_y_y_xxxy_0[i] * pa_y[i] - ta1_y_y_xxxy_1[i] * pc_y[i];

        ta1_y_yy_xxxz_0[i] =
            ta1_y_0_xxxz_0[i] * fe_0 - ta1_y_0_xxxz_1[i] * fe_0 + ta_y_xxxz_1[i] + ta1_y_y_xxxz_0[i] * pa_y[i] - ta1_y_y_xxxz_1[i] * pc_y[i];

        ta1_y_yy_xxyy_0[i] = ta1_y_0_xxyy_0[i] * fe_0 - ta1_y_0_xxyy_1[i] * fe_0 + 2.0 * ta1_y_y_xxy_0[i] * fe_0 - 2.0 * ta1_y_y_xxy_1[i] * fe_0 +
                             ta_y_xxyy_1[i] + ta1_y_y_xxyy_0[i] * pa_y[i] - ta1_y_y_xxyy_1[i] * pc_y[i];

        ta1_y_yy_xxyz_0[i] = ta1_y_0_xxyz_0[i] * fe_0 - ta1_y_0_xxyz_1[i] * fe_0 + ta1_y_y_xxz_0[i] * fe_0 - ta1_y_y_xxz_1[i] * fe_0 +
                             ta_y_xxyz_1[i] + ta1_y_y_xxyz_0[i] * pa_y[i] - ta1_y_y_xxyz_1[i] * pc_y[i];

        ta1_y_yy_xxzz_0[i] =
            ta1_y_0_xxzz_0[i] * fe_0 - ta1_y_0_xxzz_1[i] * fe_0 + ta_y_xxzz_1[i] + ta1_y_y_xxzz_0[i] * pa_y[i] - ta1_y_y_xxzz_1[i] * pc_y[i];

        ta1_y_yy_xyyy_0[i] = ta1_y_0_xyyy_0[i] * fe_0 - ta1_y_0_xyyy_1[i] * fe_0 + 3.0 * ta1_y_y_xyy_0[i] * fe_0 - 3.0 * ta1_y_y_xyy_1[i] * fe_0 +
                             ta_y_xyyy_1[i] + ta1_y_y_xyyy_0[i] * pa_y[i] - ta1_y_y_xyyy_1[i] * pc_y[i];

        ta1_y_yy_xyyz_0[i] = ta1_y_0_xyyz_0[i] * fe_0 - ta1_y_0_xyyz_1[i] * fe_0 + 2.0 * ta1_y_y_xyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyz_1[i] * fe_0 +
                             ta_y_xyyz_1[i] + ta1_y_y_xyyz_0[i] * pa_y[i] - ta1_y_y_xyyz_1[i] * pc_y[i];

        ta1_y_yy_xyzz_0[i] = ta1_y_0_xyzz_0[i] * fe_0 - ta1_y_0_xyzz_1[i] * fe_0 + ta1_y_y_xzz_0[i] * fe_0 - ta1_y_y_xzz_1[i] * fe_0 +
                             ta_y_xyzz_1[i] + ta1_y_y_xyzz_0[i] * pa_y[i] - ta1_y_y_xyzz_1[i] * pc_y[i];

        ta1_y_yy_xzzz_0[i] =
            ta1_y_0_xzzz_0[i] * fe_0 - ta1_y_0_xzzz_1[i] * fe_0 + ta_y_xzzz_1[i] + ta1_y_y_xzzz_0[i] * pa_y[i] - ta1_y_y_xzzz_1[i] * pc_y[i];

        ta1_y_yy_yyyy_0[i] = ta1_y_0_yyyy_0[i] * fe_0 - ta1_y_0_yyyy_1[i] * fe_0 + 4.0 * ta1_y_y_yyy_0[i] * fe_0 - 4.0 * ta1_y_y_yyy_1[i] * fe_0 +
                             ta_y_yyyy_1[i] + ta1_y_y_yyyy_0[i] * pa_y[i] - ta1_y_y_yyyy_1[i] * pc_y[i];

        ta1_y_yy_yyyz_0[i] = ta1_y_0_yyyz_0[i] * fe_0 - ta1_y_0_yyyz_1[i] * fe_0 + 3.0 * ta1_y_y_yyz_0[i] * fe_0 - 3.0 * ta1_y_y_yyz_1[i] * fe_0 +
                             ta_y_yyyz_1[i] + ta1_y_y_yyyz_0[i] * pa_y[i] - ta1_y_y_yyyz_1[i] * pc_y[i];

        ta1_y_yy_yyzz_0[i] = ta1_y_0_yyzz_0[i] * fe_0 - ta1_y_0_yyzz_1[i] * fe_0 + 2.0 * ta1_y_y_yzz_0[i] * fe_0 - 2.0 * ta1_y_y_yzz_1[i] * fe_0 +
                             ta_y_yyzz_1[i] + ta1_y_y_yyzz_0[i] * pa_y[i] - ta1_y_y_yyzz_1[i] * pc_y[i];

        ta1_y_yy_yzzz_0[i] = ta1_y_0_yzzz_0[i] * fe_0 - ta1_y_0_yzzz_1[i] * fe_0 + ta1_y_y_zzz_0[i] * fe_0 - ta1_y_y_zzz_1[i] * fe_0 +
                             ta_y_yzzz_1[i] + ta1_y_y_yzzz_0[i] * pa_y[i] - ta1_y_y_yzzz_1[i] * pc_y[i];

        ta1_y_yy_zzzz_0[i] =
            ta1_y_0_zzzz_0[i] * fe_0 - ta1_y_0_zzzz_1[i] * fe_0 + ta_y_zzzz_1[i] + ta1_y_y_zzzz_0[i] * pa_y[i] - ta1_y_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 150-165 components of targeted buffer : DG

    auto ta1_y_yz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 150);

    auto ta1_y_yz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 151);

    auto ta1_y_yz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 152);

    auto ta1_y_yz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 153);

    auto ta1_y_yz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 154);

    auto ta1_y_yz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 155);

    auto ta1_y_yz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 156);

    auto ta1_y_yz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 157);

    auto ta1_y_yz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 158);

    auto ta1_y_yz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 159);

    auto ta1_y_yz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 160);

    auto ta1_y_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 161);

    auto ta1_y_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 162);

    auto ta1_y_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 163);

    auto ta1_y_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 164);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_y_xxxx_0,  \
                             ta1_y_y_xxxx_1,  \
                             ta1_y_y_xxxy_0,  \
                             ta1_y_y_xxxy_1,  \
                             ta1_y_y_xxy_0,   \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xxyy_0,  \
                             ta1_y_y_xxyy_1,  \
                             ta1_y_y_xxyz_0,  \
                             ta1_y_y_xxyz_1,  \
                             ta1_y_y_xyy_0,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_xyyy_0,  \
                             ta1_y_y_xyyy_1,  \
                             ta1_y_y_xyyz_0,  \
                             ta1_y_y_xyyz_1,  \
                             ta1_y_y_xyz_0,   \
                             ta1_y_y_xyz_1,   \
                             ta1_y_y_xyzz_0,  \
                             ta1_y_y_xyzz_1,  \
                             ta1_y_y_yyy_0,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_y_yyyy_0,  \
                             ta1_y_y_yyyy_1,  \
                             ta1_y_y_yyyz_0,  \
                             ta1_y_y_yyyz_1,  \
                             ta1_y_y_yyz_0,   \
                             ta1_y_y_yyz_1,   \
                             ta1_y_y_yyzz_0,  \
                             ta1_y_y_yyzz_1,  \
                             ta1_y_y_yzz_0,   \
                             ta1_y_y_yzz_1,   \
                             ta1_y_y_yzzz_0,  \
                             ta1_y_y_yzzz_1,  \
                             ta1_y_yz_xxxx_0, \
                             ta1_y_yz_xxxy_0, \
                             ta1_y_yz_xxxz_0, \
                             ta1_y_yz_xxyy_0, \
                             ta1_y_yz_xxyz_0, \
                             ta1_y_yz_xxzz_0, \
                             ta1_y_yz_xyyy_0, \
                             ta1_y_yz_xyyz_0, \
                             ta1_y_yz_xyzz_0, \
                             ta1_y_yz_xzzz_0, \
                             ta1_y_yz_yyyy_0, \
                             ta1_y_yz_yyyz_0, \
                             ta1_y_yz_yyzz_0, \
                             ta1_y_yz_yzzz_0, \
                             ta1_y_yz_zzzz_0, \
                             ta1_y_z_xxxz_0,  \
                             ta1_y_z_xxxz_1,  \
                             ta1_y_z_xxzz_0,  \
                             ta1_y_z_xxzz_1,  \
                             ta1_y_z_xzzz_0,  \
                             ta1_y_z_xzzz_1,  \
                             ta1_y_z_zzzz_0,  \
                             ta1_y_z_zzzz_1,  \
                             ta_z_xxxz_1,     \
                             ta_z_xxzz_1,     \
                             ta_z_xzzz_1,     \
                             ta_z_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yz_xxxx_0[i] = ta1_y_y_xxxx_0[i] * pa_z[i] - ta1_y_y_xxxx_1[i] * pc_z[i];

        ta1_y_yz_xxxy_0[i] = ta1_y_y_xxxy_0[i] * pa_z[i] - ta1_y_y_xxxy_1[i] * pc_z[i];

        ta1_y_yz_xxxz_0[i] = ta_z_xxxz_1[i] + ta1_y_z_xxxz_0[i] * pa_y[i] - ta1_y_z_xxxz_1[i] * pc_y[i];

        ta1_y_yz_xxyy_0[i] = ta1_y_y_xxyy_0[i] * pa_z[i] - ta1_y_y_xxyy_1[i] * pc_z[i];

        ta1_y_yz_xxyz_0[i] = ta1_y_y_xxy_0[i] * fe_0 - ta1_y_y_xxy_1[i] * fe_0 + ta1_y_y_xxyz_0[i] * pa_z[i] - ta1_y_y_xxyz_1[i] * pc_z[i];

        ta1_y_yz_xxzz_0[i] = ta_z_xxzz_1[i] + ta1_y_z_xxzz_0[i] * pa_y[i] - ta1_y_z_xxzz_1[i] * pc_y[i];

        ta1_y_yz_xyyy_0[i] = ta1_y_y_xyyy_0[i] * pa_z[i] - ta1_y_y_xyyy_1[i] * pc_z[i];

        ta1_y_yz_xyyz_0[i] = ta1_y_y_xyy_0[i] * fe_0 - ta1_y_y_xyy_1[i] * fe_0 + ta1_y_y_xyyz_0[i] * pa_z[i] - ta1_y_y_xyyz_1[i] * pc_z[i];

        ta1_y_yz_xyzz_0[i] =
            2.0 * ta1_y_y_xyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyz_1[i] * fe_0 + ta1_y_y_xyzz_0[i] * pa_z[i] - ta1_y_y_xyzz_1[i] * pc_z[i];

        ta1_y_yz_xzzz_0[i] = ta_z_xzzz_1[i] + ta1_y_z_xzzz_0[i] * pa_y[i] - ta1_y_z_xzzz_1[i] * pc_y[i];

        ta1_y_yz_yyyy_0[i] = ta1_y_y_yyyy_0[i] * pa_z[i] - ta1_y_y_yyyy_1[i] * pc_z[i];

        ta1_y_yz_yyyz_0[i] = ta1_y_y_yyy_0[i] * fe_0 - ta1_y_y_yyy_1[i] * fe_0 + ta1_y_y_yyyz_0[i] * pa_z[i] - ta1_y_y_yyyz_1[i] * pc_z[i];

        ta1_y_yz_yyzz_0[i] =
            2.0 * ta1_y_y_yyz_0[i] * fe_0 - 2.0 * ta1_y_y_yyz_1[i] * fe_0 + ta1_y_y_yyzz_0[i] * pa_z[i] - ta1_y_y_yyzz_1[i] * pc_z[i];

        ta1_y_yz_yzzz_0[i] =
            3.0 * ta1_y_y_yzz_0[i] * fe_0 - 3.0 * ta1_y_y_yzz_1[i] * fe_0 + ta1_y_y_yzzz_0[i] * pa_z[i] - ta1_y_y_yzzz_1[i] * pc_z[i];

        ta1_y_yz_zzzz_0[i] = ta_z_zzzz_1[i] + ta1_y_z_zzzz_0[i] * pa_y[i] - ta1_y_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 165-180 components of targeted buffer : DG

    auto ta1_y_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 165);

    auto ta1_y_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 166);

    auto ta1_y_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 167);

    auto ta1_y_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 168);

    auto ta1_y_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 169);

    auto ta1_y_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 170);

    auto ta1_y_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 171);

    auto ta1_y_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 172);

    auto ta1_y_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 173);

    auto ta1_y_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 174);

    auto ta1_y_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 175);

    auto ta1_y_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 176);

    auto ta1_y_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 177);

    auto ta1_y_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 178);

    auto ta1_y_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 179);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_0_xxxx_0,  \
                             ta1_y_0_xxxx_1,  \
                             ta1_y_0_xxxy_0,  \
                             ta1_y_0_xxxy_1,  \
                             ta1_y_0_xxxz_0,  \
                             ta1_y_0_xxxz_1,  \
                             ta1_y_0_xxyy_0,  \
                             ta1_y_0_xxyy_1,  \
                             ta1_y_0_xxyz_0,  \
                             ta1_y_0_xxyz_1,  \
                             ta1_y_0_xxzz_0,  \
                             ta1_y_0_xxzz_1,  \
                             ta1_y_0_xyyy_0,  \
                             ta1_y_0_xyyy_1,  \
                             ta1_y_0_xyyz_0,  \
                             ta1_y_0_xyyz_1,  \
                             ta1_y_0_xyzz_0,  \
                             ta1_y_0_xyzz_1,  \
                             ta1_y_0_xzzz_0,  \
                             ta1_y_0_xzzz_1,  \
                             ta1_y_0_yyyy_0,  \
                             ta1_y_0_yyyy_1,  \
                             ta1_y_0_yyyz_0,  \
                             ta1_y_0_yyyz_1,  \
                             ta1_y_0_yyzz_0,  \
                             ta1_y_0_yyzz_1,  \
                             ta1_y_0_yzzz_0,  \
                             ta1_y_0_yzzz_1,  \
                             ta1_y_0_zzzz_0,  \
                             ta1_y_0_zzzz_1,  \
                             ta1_y_z_xxx_0,   \
                             ta1_y_z_xxx_1,   \
                             ta1_y_z_xxxx_0,  \
                             ta1_y_z_xxxx_1,  \
                             ta1_y_z_xxxy_0,  \
                             ta1_y_z_xxxy_1,  \
                             ta1_y_z_xxxz_0,  \
                             ta1_y_z_xxxz_1,  \
                             ta1_y_z_xxy_0,   \
                             ta1_y_z_xxy_1,   \
                             ta1_y_z_xxyy_0,  \
                             ta1_y_z_xxyy_1,  \
                             ta1_y_z_xxyz_0,  \
                             ta1_y_z_xxyz_1,  \
                             ta1_y_z_xxz_0,   \
                             ta1_y_z_xxz_1,   \
                             ta1_y_z_xxzz_0,  \
                             ta1_y_z_xxzz_1,  \
                             ta1_y_z_xyy_0,   \
                             ta1_y_z_xyy_1,   \
                             ta1_y_z_xyyy_0,  \
                             ta1_y_z_xyyy_1,  \
                             ta1_y_z_xyyz_0,  \
                             ta1_y_z_xyyz_1,  \
                             ta1_y_z_xyz_0,   \
                             ta1_y_z_xyz_1,   \
                             ta1_y_z_xyzz_0,  \
                             ta1_y_z_xyzz_1,  \
                             ta1_y_z_xzz_0,   \
                             ta1_y_z_xzz_1,   \
                             ta1_y_z_xzzz_0,  \
                             ta1_y_z_xzzz_1,  \
                             ta1_y_z_yyy_0,   \
                             ta1_y_z_yyy_1,   \
                             ta1_y_z_yyyy_0,  \
                             ta1_y_z_yyyy_1,  \
                             ta1_y_z_yyyz_0,  \
                             ta1_y_z_yyyz_1,  \
                             ta1_y_z_yyz_0,   \
                             ta1_y_z_yyz_1,   \
                             ta1_y_z_yyzz_0,  \
                             ta1_y_z_yyzz_1,  \
                             ta1_y_z_yzz_0,   \
                             ta1_y_z_yzz_1,   \
                             ta1_y_z_yzzz_0,  \
                             ta1_y_z_yzzz_1,  \
                             ta1_y_z_zzz_0,   \
                             ta1_y_z_zzz_1,   \
                             ta1_y_z_zzzz_0,  \
                             ta1_y_z_zzzz_1,  \
                             ta1_y_zz_xxxx_0, \
                             ta1_y_zz_xxxy_0, \
                             ta1_y_zz_xxxz_0, \
                             ta1_y_zz_xxyy_0, \
                             ta1_y_zz_xxyz_0, \
                             ta1_y_zz_xxzz_0, \
                             ta1_y_zz_xyyy_0, \
                             ta1_y_zz_xyyz_0, \
                             ta1_y_zz_xyzz_0, \
                             ta1_y_zz_xzzz_0, \
                             ta1_y_zz_yyyy_0, \
                             ta1_y_zz_yyyz_0, \
                             ta1_y_zz_yyzz_0, \
                             ta1_y_zz_yzzz_0, \
                             ta1_y_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zz_xxxx_0[i] = ta1_y_0_xxxx_0[i] * fe_0 - ta1_y_0_xxxx_1[i] * fe_0 + ta1_y_z_xxxx_0[i] * pa_z[i] - ta1_y_z_xxxx_1[i] * pc_z[i];

        ta1_y_zz_xxxy_0[i] = ta1_y_0_xxxy_0[i] * fe_0 - ta1_y_0_xxxy_1[i] * fe_0 + ta1_y_z_xxxy_0[i] * pa_z[i] - ta1_y_z_xxxy_1[i] * pc_z[i];

        ta1_y_zz_xxxz_0[i] = ta1_y_0_xxxz_0[i] * fe_0 - ta1_y_0_xxxz_1[i] * fe_0 + ta1_y_z_xxx_0[i] * fe_0 - ta1_y_z_xxx_1[i] * fe_0 +
                             ta1_y_z_xxxz_0[i] * pa_z[i] - ta1_y_z_xxxz_1[i] * pc_z[i];

        ta1_y_zz_xxyy_0[i] = ta1_y_0_xxyy_0[i] * fe_0 - ta1_y_0_xxyy_1[i] * fe_0 + ta1_y_z_xxyy_0[i] * pa_z[i] - ta1_y_z_xxyy_1[i] * pc_z[i];

        ta1_y_zz_xxyz_0[i] = ta1_y_0_xxyz_0[i] * fe_0 - ta1_y_0_xxyz_1[i] * fe_0 + ta1_y_z_xxy_0[i] * fe_0 - ta1_y_z_xxy_1[i] * fe_0 +
                             ta1_y_z_xxyz_0[i] * pa_z[i] - ta1_y_z_xxyz_1[i] * pc_z[i];

        ta1_y_zz_xxzz_0[i] = ta1_y_0_xxzz_0[i] * fe_0 - ta1_y_0_xxzz_1[i] * fe_0 + 2.0 * ta1_y_z_xxz_0[i] * fe_0 - 2.0 * ta1_y_z_xxz_1[i] * fe_0 +
                             ta1_y_z_xxzz_0[i] * pa_z[i] - ta1_y_z_xxzz_1[i] * pc_z[i];

        ta1_y_zz_xyyy_0[i] = ta1_y_0_xyyy_0[i] * fe_0 - ta1_y_0_xyyy_1[i] * fe_0 + ta1_y_z_xyyy_0[i] * pa_z[i] - ta1_y_z_xyyy_1[i] * pc_z[i];

        ta1_y_zz_xyyz_0[i] = ta1_y_0_xyyz_0[i] * fe_0 - ta1_y_0_xyyz_1[i] * fe_0 + ta1_y_z_xyy_0[i] * fe_0 - ta1_y_z_xyy_1[i] * fe_0 +
                             ta1_y_z_xyyz_0[i] * pa_z[i] - ta1_y_z_xyyz_1[i] * pc_z[i];

        ta1_y_zz_xyzz_0[i] = ta1_y_0_xyzz_0[i] * fe_0 - ta1_y_0_xyzz_1[i] * fe_0 + 2.0 * ta1_y_z_xyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyz_1[i] * fe_0 +
                             ta1_y_z_xyzz_0[i] * pa_z[i] - ta1_y_z_xyzz_1[i] * pc_z[i];

        ta1_y_zz_xzzz_0[i] = ta1_y_0_xzzz_0[i] * fe_0 - ta1_y_0_xzzz_1[i] * fe_0 + 3.0 * ta1_y_z_xzz_0[i] * fe_0 - 3.0 * ta1_y_z_xzz_1[i] * fe_0 +
                             ta1_y_z_xzzz_0[i] * pa_z[i] - ta1_y_z_xzzz_1[i] * pc_z[i];

        ta1_y_zz_yyyy_0[i] = ta1_y_0_yyyy_0[i] * fe_0 - ta1_y_0_yyyy_1[i] * fe_0 + ta1_y_z_yyyy_0[i] * pa_z[i] - ta1_y_z_yyyy_1[i] * pc_z[i];

        ta1_y_zz_yyyz_0[i] = ta1_y_0_yyyz_0[i] * fe_0 - ta1_y_0_yyyz_1[i] * fe_0 + ta1_y_z_yyy_0[i] * fe_0 - ta1_y_z_yyy_1[i] * fe_0 +
                             ta1_y_z_yyyz_0[i] * pa_z[i] - ta1_y_z_yyyz_1[i] * pc_z[i];

        ta1_y_zz_yyzz_0[i] = ta1_y_0_yyzz_0[i] * fe_0 - ta1_y_0_yyzz_1[i] * fe_0 + 2.0 * ta1_y_z_yyz_0[i] * fe_0 - 2.0 * ta1_y_z_yyz_1[i] * fe_0 +
                             ta1_y_z_yyzz_0[i] * pa_z[i] - ta1_y_z_yyzz_1[i] * pc_z[i];

        ta1_y_zz_yzzz_0[i] = ta1_y_0_yzzz_0[i] * fe_0 - ta1_y_0_yzzz_1[i] * fe_0 + 3.0 * ta1_y_z_yzz_0[i] * fe_0 - 3.0 * ta1_y_z_yzz_1[i] * fe_0 +
                             ta1_y_z_yzzz_0[i] * pa_z[i] - ta1_y_z_yzzz_1[i] * pc_z[i];

        ta1_y_zz_zzzz_0[i] = ta1_y_0_zzzz_0[i] * fe_0 - ta1_y_0_zzzz_1[i] * fe_0 + 4.0 * ta1_y_z_zzz_0[i] * fe_0 - 4.0 * ta1_y_z_zzz_1[i] * fe_0 +
                             ta1_y_z_zzzz_0[i] * pa_z[i] - ta1_y_z_zzzz_1[i] * pc_z[i];
    }

    // Set up 180-195 components of targeted buffer : DG

    auto ta1_z_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 180);

    auto ta1_z_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 181);

    auto ta1_z_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 182);

    auto ta1_z_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 183);

    auto ta1_z_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 184);

    auto ta1_z_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 185);

    auto ta1_z_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 186);

    auto ta1_z_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 187);

    auto ta1_z_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 188);

    auto ta1_z_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 189);

    auto ta1_z_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 190);

    auto ta1_z_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 191);

    auto ta1_z_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 192);

    auto ta1_z_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 193);

    auto ta1_z_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 194);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_0_xxxx_0,  \
                             ta1_z_0_xxxx_1,  \
                             ta1_z_0_xxxy_0,  \
                             ta1_z_0_xxxy_1,  \
                             ta1_z_0_xxxz_0,  \
                             ta1_z_0_xxxz_1,  \
                             ta1_z_0_xxyy_0,  \
                             ta1_z_0_xxyy_1,  \
                             ta1_z_0_xxyz_0,  \
                             ta1_z_0_xxyz_1,  \
                             ta1_z_0_xxzz_0,  \
                             ta1_z_0_xxzz_1,  \
                             ta1_z_0_xyyy_0,  \
                             ta1_z_0_xyyy_1,  \
                             ta1_z_0_xyyz_0,  \
                             ta1_z_0_xyyz_1,  \
                             ta1_z_0_xyzz_0,  \
                             ta1_z_0_xyzz_1,  \
                             ta1_z_0_xzzz_0,  \
                             ta1_z_0_xzzz_1,  \
                             ta1_z_0_yyyy_0,  \
                             ta1_z_0_yyyy_1,  \
                             ta1_z_0_yyyz_0,  \
                             ta1_z_0_yyyz_1,  \
                             ta1_z_0_yyzz_0,  \
                             ta1_z_0_yyzz_1,  \
                             ta1_z_0_yzzz_0,  \
                             ta1_z_0_yzzz_1,  \
                             ta1_z_0_zzzz_0,  \
                             ta1_z_0_zzzz_1,  \
                             ta1_z_x_xxx_0,   \
                             ta1_z_x_xxx_1,   \
                             ta1_z_x_xxxx_0,  \
                             ta1_z_x_xxxx_1,  \
                             ta1_z_x_xxxy_0,  \
                             ta1_z_x_xxxy_1,  \
                             ta1_z_x_xxxz_0,  \
                             ta1_z_x_xxxz_1,  \
                             ta1_z_x_xxy_0,   \
                             ta1_z_x_xxy_1,   \
                             ta1_z_x_xxyy_0,  \
                             ta1_z_x_xxyy_1,  \
                             ta1_z_x_xxyz_0,  \
                             ta1_z_x_xxyz_1,  \
                             ta1_z_x_xxz_0,   \
                             ta1_z_x_xxz_1,   \
                             ta1_z_x_xxzz_0,  \
                             ta1_z_x_xxzz_1,  \
                             ta1_z_x_xyy_0,   \
                             ta1_z_x_xyy_1,   \
                             ta1_z_x_xyyy_0,  \
                             ta1_z_x_xyyy_1,  \
                             ta1_z_x_xyyz_0,  \
                             ta1_z_x_xyyz_1,  \
                             ta1_z_x_xyz_0,   \
                             ta1_z_x_xyz_1,   \
                             ta1_z_x_xyzz_0,  \
                             ta1_z_x_xyzz_1,  \
                             ta1_z_x_xzz_0,   \
                             ta1_z_x_xzz_1,   \
                             ta1_z_x_xzzz_0,  \
                             ta1_z_x_xzzz_1,  \
                             ta1_z_x_yyy_0,   \
                             ta1_z_x_yyy_1,   \
                             ta1_z_x_yyyy_0,  \
                             ta1_z_x_yyyy_1,  \
                             ta1_z_x_yyyz_0,  \
                             ta1_z_x_yyyz_1,  \
                             ta1_z_x_yyz_0,   \
                             ta1_z_x_yyz_1,   \
                             ta1_z_x_yyzz_0,  \
                             ta1_z_x_yyzz_1,  \
                             ta1_z_x_yzz_0,   \
                             ta1_z_x_yzz_1,   \
                             ta1_z_x_yzzz_0,  \
                             ta1_z_x_yzzz_1,  \
                             ta1_z_x_zzz_0,   \
                             ta1_z_x_zzz_1,   \
                             ta1_z_x_zzzz_0,  \
                             ta1_z_x_zzzz_1,  \
                             ta1_z_xx_xxxx_0, \
                             ta1_z_xx_xxxy_0, \
                             ta1_z_xx_xxxz_0, \
                             ta1_z_xx_xxyy_0, \
                             ta1_z_xx_xxyz_0, \
                             ta1_z_xx_xxzz_0, \
                             ta1_z_xx_xyyy_0, \
                             ta1_z_xx_xyyz_0, \
                             ta1_z_xx_xyzz_0, \
                             ta1_z_xx_xzzz_0, \
                             ta1_z_xx_yyyy_0, \
                             ta1_z_xx_yyyz_0, \
                             ta1_z_xx_yyzz_0, \
                             ta1_z_xx_yzzz_0, \
                             ta1_z_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xx_xxxx_0[i] = ta1_z_0_xxxx_0[i] * fe_0 - ta1_z_0_xxxx_1[i] * fe_0 + 4.0 * ta1_z_x_xxx_0[i] * fe_0 - 4.0 * ta1_z_x_xxx_1[i] * fe_0 +
                             ta1_z_x_xxxx_0[i] * pa_x[i] - ta1_z_x_xxxx_1[i] * pc_x[i];

        ta1_z_xx_xxxy_0[i] = ta1_z_0_xxxy_0[i] * fe_0 - ta1_z_0_xxxy_1[i] * fe_0 + 3.0 * ta1_z_x_xxy_0[i] * fe_0 - 3.0 * ta1_z_x_xxy_1[i] * fe_0 +
                             ta1_z_x_xxxy_0[i] * pa_x[i] - ta1_z_x_xxxy_1[i] * pc_x[i];

        ta1_z_xx_xxxz_0[i] = ta1_z_0_xxxz_0[i] * fe_0 - ta1_z_0_xxxz_1[i] * fe_0 + 3.0 * ta1_z_x_xxz_0[i] * fe_0 - 3.0 * ta1_z_x_xxz_1[i] * fe_0 +
                             ta1_z_x_xxxz_0[i] * pa_x[i] - ta1_z_x_xxxz_1[i] * pc_x[i];

        ta1_z_xx_xxyy_0[i] = ta1_z_0_xxyy_0[i] * fe_0 - ta1_z_0_xxyy_1[i] * fe_0 + 2.0 * ta1_z_x_xyy_0[i] * fe_0 - 2.0 * ta1_z_x_xyy_1[i] * fe_0 +
                             ta1_z_x_xxyy_0[i] * pa_x[i] - ta1_z_x_xxyy_1[i] * pc_x[i];

        ta1_z_xx_xxyz_0[i] = ta1_z_0_xxyz_0[i] * fe_0 - ta1_z_0_xxyz_1[i] * fe_0 + 2.0 * ta1_z_x_xyz_0[i] * fe_0 - 2.0 * ta1_z_x_xyz_1[i] * fe_0 +
                             ta1_z_x_xxyz_0[i] * pa_x[i] - ta1_z_x_xxyz_1[i] * pc_x[i];

        ta1_z_xx_xxzz_0[i] = ta1_z_0_xxzz_0[i] * fe_0 - ta1_z_0_xxzz_1[i] * fe_0 + 2.0 * ta1_z_x_xzz_0[i] * fe_0 - 2.0 * ta1_z_x_xzz_1[i] * fe_0 +
                             ta1_z_x_xxzz_0[i] * pa_x[i] - ta1_z_x_xxzz_1[i] * pc_x[i];

        ta1_z_xx_xyyy_0[i] = ta1_z_0_xyyy_0[i] * fe_0 - ta1_z_0_xyyy_1[i] * fe_0 + ta1_z_x_yyy_0[i] * fe_0 - ta1_z_x_yyy_1[i] * fe_0 +
                             ta1_z_x_xyyy_0[i] * pa_x[i] - ta1_z_x_xyyy_1[i] * pc_x[i];

        ta1_z_xx_xyyz_0[i] = ta1_z_0_xyyz_0[i] * fe_0 - ta1_z_0_xyyz_1[i] * fe_0 + ta1_z_x_yyz_0[i] * fe_0 - ta1_z_x_yyz_1[i] * fe_0 +
                             ta1_z_x_xyyz_0[i] * pa_x[i] - ta1_z_x_xyyz_1[i] * pc_x[i];

        ta1_z_xx_xyzz_0[i] = ta1_z_0_xyzz_0[i] * fe_0 - ta1_z_0_xyzz_1[i] * fe_0 + ta1_z_x_yzz_0[i] * fe_0 - ta1_z_x_yzz_1[i] * fe_0 +
                             ta1_z_x_xyzz_0[i] * pa_x[i] - ta1_z_x_xyzz_1[i] * pc_x[i];

        ta1_z_xx_xzzz_0[i] = ta1_z_0_xzzz_0[i] * fe_0 - ta1_z_0_xzzz_1[i] * fe_0 + ta1_z_x_zzz_0[i] * fe_0 - ta1_z_x_zzz_1[i] * fe_0 +
                             ta1_z_x_xzzz_0[i] * pa_x[i] - ta1_z_x_xzzz_1[i] * pc_x[i];

        ta1_z_xx_yyyy_0[i] = ta1_z_0_yyyy_0[i] * fe_0 - ta1_z_0_yyyy_1[i] * fe_0 + ta1_z_x_yyyy_0[i] * pa_x[i] - ta1_z_x_yyyy_1[i] * pc_x[i];

        ta1_z_xx_yyyz_0[i] = ta1_z_0_yyyz_0[i] * fe_0 - ta1_z_0_yyyz_1[i] * fe_0 + ta1_z_x_yyyz_0[i] * pa_x[i] - ta1_z_x_yyyz_1[i] * pc_x[i];

        ta1_z_xx_yyzz_0[i] = ta1_z_0_yyzz_0[i] * fe_0 - ta1_z_0_yyzz_1[i] * fe_0 + ta1_z_x_yyzz_0[i] * pa_x[i] - ta1_z_x_yyzz_1[i] * pc_x[i];

        ta1_z_xx_yzzz_0[i] = ta1_z_0_yzzz_0[i] * fe_0 - ta1_z_0_yzzz_1[i] * fe_0 + ta1_z_x_yzzz_0[i] * pa_x[i] - ta1_z_x_yzzz_1[i] * pc_x[i];

        ta1_z_xx_zzzz_0[i] = ta1_z_0_zzzz_0[i] * fe_0 - ta1_z_0_zzzz_1[i] * fe_0 + ta1_z_x_zzzz_0[i] * pa_x[i] - ta1_z_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 195-210 components of targeted buffer : DG

    auto ta1_z_xy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 195);

    auto ta1_z_xy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 196);

    auto ta1_z_xy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 197);

    auto ta1_z_xy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 198);

    auto ta1_z_xy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 199);

    auto ta1_z_xy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 200);

    auto ta1_z_xy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 201);

    auto ta1_z_xy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 202);

    auto ta1_z_xy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 203);

    auto ta1_z_xy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 204);

    auto ta1_z_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 205);

    auto ta1_z_xy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 206);

    auto ta1_z_xy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 207);

    auto ta1_z_xy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 208);

    auto ta1_z_xy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 209);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_x_xxxx_0,  \
                             ta1_z_x_xxxx_1,  \
                             ta1_z_x_xxxz_0,  \
                             ta1_z_x_xxxz_1,  \
                             ta1_z_x_xxzz_0,  \
                             ta1_z_x_xxzz_1,  \
                             ta1_z_x_xzzz_0,  \
                             ta1_z_x_xzzz_1,  \
                             ta1_z_xy_xxxx_0, \
                             ta1_z_xy_xxxy_0, \
                             ta1_z_xy_xxxz_0, \
                             ta1_z_xy_xxyy_0, \
                             ta1_z_xy_xxyz_0, \
                             ta1_z_xy_xxzz_0, \
                             ta1_z_xy_xyyy_0, \
                             ta1_z_xy_xyyz_0, \
                             ta1_z_xy_xyzz_0, \
                             ta1_z_xy_xzzz_0, \
                             ta1_z_xy_yyyy_0, \
                             ta1_z_xy_yyyz_0, \
                             ta1_z_xy_yyzz_0, \
                             ta1_z_xy_yzzz_0, \
                             ta1_z_xy_zzzz_0, \
                             ta1_z_y_xxxy_0,  \
                             ta1_z_y_xxxy_1,  \
                             ta1_z_y_xxy_0,   \
                             ta1_z_y_xxy_1,   \
                             ta1_z_y_xxyy_0,  \
                             ta1_z_y_xxyy_1,  \
                             ta1_z_y_xxyz_0,  \
                             ta1_z_y_xxyz_1,  \
                             ta1_z_y_xyy_0,   \
                             ta1_z_y_xyy_1,   \
                             ta1_z_y_xyyy_0,  \
                             ta1_z_y_xyyy_1,  \
                             ta1_z_y_xyyz_0,  \
                             ta1_z_y_xyyz_1,  \
                             ta1_z_y_xyz_0,   \
                             ta1_z_y_xyz_1,   \
                             ta1_z_y_xyzz_0,  \
                             ta1_z_y_xyzz_1,  \
                             ta1_z_y_yyy_0,   \
                             ta1_z_y_yyy_1,   \
                             ta1_z_y_yyyy_0,  \
                             ta1_z_y_yyyy_1,  \
                             ta1_z_y_yyyz_0,  \
                             ta1_z_y_yyyz_1,  \
                             ta1_z_y_yyz_0,   \
                             ta1_z_y_yyz_1,   \
                             ta1_z_y_yyzz_0,  \
                             ta1_z_y_yyzz_1,  \
                             ta1_z_y_yzz_0,   \
                             ta1_z_y_yzz_1,   \
                             ta1_z_y_yzzz_0,  \
                             ta1_z_y_yzzz_1,  \
                             ta1_z_y_zzzz_0,  \
                             ta1_z_y_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xy_xxxx_0[i] = ta1_z_x_xxxx_0[i] * pa_y[i] - ta1_z_x_xxxx_1[i] * pc_y[i];

        ta1_z_xy_xxxy_0[i] =
            3.0 * ta1_z_y_xxy_0[i] * fe_0 - 3.0 * ta1_z_y_xxy_1[i] * fe_0 + ta1_z_y_xxxy_0[i] * pa_x[i] - ta1_z_y_xxxy_1[i] * pc_x[i];

        ta1_z_xy_xxxz_0[i] = ta1_z_x_xxxz_0[i] * pa_y[i] - ta1_z_x_xxxz_1[i] * pc_y[i];

        ta1_z_xy_xxyy_0[i] =
            2.0 * ta1_z_y_xyy_0[i] * fe_0 - 2.0 * ta1_z_y_xyy_1[i] * fe_0 + ta1_z_y_xxyy_0[i] * pa_x[i] - ta1_z_y_xxyy_1[i] * pc_x[i];

        ta1_z_xy_xxyz_0[i] =
            2.0 * ta1_z_y_xyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyz_1[i] * fe_0 + ta1_z_y_xxyz_0[i] * pa_x[i] - ta1_z_y_xxyz_1[i] * pc_x[i];

        ta1_z_xy_xxzz_0[i] = ta1_z_x_xxzz_0[i] * pa_y[i] - ta1_z_x_xxzz_1[i] * pc_y[i];

        ta1_z_xy_xyyy_0[i] = ta1_z_y_yyy_0[i] * fe_0 - ta1_z_y_yyy_1[i] * fe_0 + ta1_z_y_xyyy_0[i] * pa_x[i] - ta1_z_y_xyyy_1[i] * pc_x[i];

        ta1_z_xy_xyyz_0[i] = ta1_z_y_yyz_0[i] * fe_0 - ta1_z_y_yyz_1[i] * fe_0 + ta1_z_y_xyyz_0[i] * pa_x[i] - ta1_z_y_xyyz_1[i] * pc_x[i];

        ta1_z_xy_xyzz_0[i] = ta1_z_y_yzz_0[i] * fe_0 - ta1_z_y_yzz_1[i] * fe_0 + ta1_z_y_xyzz_0[i] * pa_x[i] - ta1_z_y_xyzz_1[i] * pc_x[i];

        ta1_z_xy_xzzz_0[i] = ta1_z_x_xzzz_0[i] * pa_y[i] - ta1_z_x_xzzz_1[i] * pc_y[i];

        ta1_z_xy_yyyy_0[i] = ta1_z_y_yyyy_0[i] * pa_x[i] - ta1_z_y_yyyy_1[i] * pc_x[i];

        ta1_z_xy_yyyz_0[i] = ta1_z_y_yyyz_0[i] * pa_x[i] - ta1_z_y_yyyz_1[i] * pc_x[i];

        ta1_z_xy_yyzz_0[i] = ta1_z_y_yyzz_0[i] * pa_x[i] - ta1_z_y_yyzz_1[i] * pc_x[i];

        ta1_z_xy_yzzz_0[i] = ta1_z_y_yzzz_0[i] * pa_x[i] - ta1_z_y_yzzz_1[i] * pc_x[i];

        ta1_z_xy_zzzz_0[i] = ta1_z_y_zzzz_0[i] * pa_x[i] - ta1_z_y_zzzz_1[i] * pc_x[i];
    }

    // Set up 210-225 components of targeted buffer : DG

    auto ta1_z_xz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 210);

    auto ta1_z_xz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 211);

    auto ta1_z_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 212);

    auto ta1_z_xz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 213);

    auto ta1_z_xz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 214);

    auto ta1_z_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 215);

    auto ta1_z_xz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 216);

    auto ta1_z_xz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 217);

    auto ta1_z_xz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 218);

    auto ta1_z_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 219);

    auto ta1_z_xz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 220);

    auto ta1_z_xz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 221);

    auto ta1_z_xz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 222);

    auto ta1_z_xz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 223);

    auto ta1_z_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 224);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_x_xxxx_0,  \
                             ta1_z_x_xxxx_1,  \
                             ta1_z_x_xxxy_0,  \
                             ta1_z_x_xxxy_1,  \
                             ta1_z_x_xxyy_0,  \
                             ta1_z_x_xxyy_1,  \
                             ta1_z_x_xyyy_0,  \
                             ta1_z_x_xyyy_1,  \
                             ta1_z_xz_xxxx_0, \
                             ta1_z_xz_xxxy_0, \
                             ta1_z_xz_xxxz_0, \
                             ta1_z_xz_xxyy_0, \
                             ta1_z_xz_xxyz_0, \
                             ta1_z_xz_xxzz_0, \
                             ta1_z_xz_xyyy_0, \
                             ta1_z_xz_xyyz_0, \
                             ta1_z_xz_xyzz_0, \
                             ta1_z_xz_xzzz_0, \
                             ta1_z_xz_yyyy_0, \
                             ta1_z_xz_yyyz_0, \
                             ta1_z_xz_yyzz_0, \
                             ta1_z_xz_yzzz_0, \
                             ta1_z_xz_zzzz_0, \
                             ta1_z_z_xxxz_0,  \
                             ta1_z_z_xxxz_1,  \
                             ta1_z_z_xxyz_0,  \
                             ta1_z_z_xxyz_1,  \
                             ta1_z_z_xxz_0,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xxzz_0,  \
                             ta1_z_z_xxzz_1,  \
                             ta1_z_z_xyyz_0,  \
                             ta1_z_z_xyyz_1,  \
                             ta1_z_z_xyz_0,   \
                             ta1_z_z_xyz_1,   \
                             ta1_z_z_xyzz_0,  \
                             ta1_z_z_xyzz_1,  \
                             ta1_z_z_xzz_0,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_xzzz_0,  \
                             ta1_z_z_xzzz_1,  \
                             ta1_z_z_yyyy_0,  \
                             ta1_z_z_yyyy_1,  \
                             ta1_z_z_yyyz_0,  \
                             ta1_z_z_yyyz_1,  \
                             ta1_z_z_yyz_0,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yyzz_0,  \
                             ta1_z_z_yyzz_1,  \
                             ta1_z_z_yzz_0,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_yzzz_0,  \
                             ta1_z_z_yzzz_1,  \
                             ta1_z_z_zzz_0,   \
                             ta1_z_z_zzz_1,   \
                             ta1_z_z_zzzz_0,  \
                             ta1_z_z_zzzz_1,  \
                             ta_x_xxxx_1,     \
                             ta_x_xxxy_1,     \
                             ta_x_xxyy_1,     \
                             ta_x_xyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xz_xxxx_0[i] = ta_x_xxxx_1[i] + ta1_z_x_xxxx_0[i] * pa_z[i] - ta1_z_x_xxxx_1[i] * pc_z[i];

        ta1_z_xz_xxxy_0[i] = ta_x_xxxy_1[i] + ta1_z_x_xxxy_0[i] * pa_z[i] - ta1_z_x_xxxy_1[i] * pc_z[i];

        ta1_z_xz_xxxz_0[i] =
            3.0 * ta1_z_z_xxz_0[i] * fe_0 - 3.0 * ta1_z_z_xxz_1[i] * fe_0 + ta1_z_z_xxxz_0[i] * pa_x[i] - ta1_z_z_xxxz_1[i] * pc_x[i];

        ta1_z_xz_xxyy_0[i] = ta_x_xxyy_1[i] + ta1_z_x_xxyy_0[i] * pa_z[i] - ta1_z_x_xxyy_1[i] * pc_z[i];

        ta1_z_xz_xxyz_0[i] =
            2.0 * ta1_z_z_xyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyz_1[i] * fe_0 + ta1_z_z_xxyz_0[i] * pa_x[i] - ta1_z_z_xxyz_1[i] * pc_x[i];

        ta1_z_xz_xxzz_0[i] =
            2.0 * ta1_z_z_xzz_0[i] * fe_0 - 2.0 * ta1_z_z_xzz_1[i] * fe_0 + ta1_z_z_xxzz_0[i] * pa_x[i] - ta1_z_z_xxzz_1[i] * pc_x[i];

        ta1_z_xz_xyyy_0[i] = ta_x_xyyy_1[i] + ta1_z_x_xyyy_0[i] * pa_z[i] - ta1_z_x_xyyy_1[i] * pc_z[i];

        ta1_z_xz_xyyz_0[i] = ta1_z_z_yyz_0[i] * fe_0 - ta1_z_z_yyz_1[i] * fe_0 + ta1_z_z_xyyz_0[i] * pa_x[i] - ta1_z_z_xyyz_1[i] * pc_x[i];

        ta1_z_xz_xyzz_0[i] = ta1_z_z_yzz_0[i] * fe_0 - ta1_z_z_yzz_1[i] * fe_0 + ta1_z_z_xyzz_0[i] * pa_x[i] - ta1_z_z_xyzz_1[i] * pc_x[i];

        ta1_z_xz_xzzz_0[i] = ta1_z_z_zzz_0[i] * fe_0 - ta1_z_z_zzz_1[i] * fe_0 + ta1_z_z_xzzz_0[i] * pa_x[i] - ta1_z_z_xzzz_1[i] * pc_x[i];

        ta1_z_xz_yyyy_0[i] = ta1_z_z_yyyy_0[i] * pa_x[i] - ta1_z_z_yyyy_1[i] * pc_x[i];

        ta1_z_xz_yyyz_0[i] = ta1_z_z_yyyz_0[i] * pa_x[i] - ta1_z_z_yyyz_1[i] * pc_x[i];

        ta1_z_xz_yyzz_0[i] = ta1_z_z_yyzz_0[i] * pa_x[i] - ta1_z_z_yyzz_1[i] * pc_x[i];

        ta1_z_xz_yzzz_0[i] = ta1_z_z_yzzz_0[i] * pa_x[i] - ta1_z_z_yzzz_1[i] * pc_x[i];

        ta1_z_xz_zzzz_0[i] = ta1_z_z_zzzz_0[i] * pa_x[i] - ta1_z_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : DG

    auto ta1_z_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 225);

    auto ta1_z_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 226);

    auto ta1_z_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 227);

    auto ta1_z_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 228);

    auto ta1_z_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 229);

    auto ta1_z_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 230);

    auto ta1_z_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 231);

    auto ta1_z_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 232);

    auto ta1_z_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 233);

    auto ta1_z_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 234);

    auto ta1_z_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 235);

    auto ta1_z_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 236);

    auto ta1_z_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 237);

    auto ta1_z_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 238);

    auto ta1_z_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 239);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_0_xxxx_0,  \
                             ta1_z_0_xxxx_1,  \
                             ta1_z_0_xxxy_0,  \
                             ta1_z_0_xxxy_1,  \
                             ta1_z_0_xxxz_0,  \
                             ta1_z_0_xxxz_1,  \
                             ta1_z_0_xxyy_0,  \
                             ta1_z_0_xxyy_1,  \
                             ta1_z_0_xxyz_0,  \
                             ta1_z_0_xxyz_1,  \
                             ta1_z_0_xxzz_0,  \
                             ta1_z_0_xxzz_1,  \
                             ta1_z_0_xyyy_0,  \
                             ta1_z_0_xyyy_1,  \
                             ta1_z_0_xyyz_0,  \
                             ta1_z_0_xyyz_1,  \
                             ta1_z_0_xyzz_0,  \
                             ta1_z_0_xyzz_1,  \
                             ta1_z_0_xzzz_0,  \
                             ta1_z_0_xzzz_1,  \
                             ta1_z_0_yyyy_0,  \
                             ta1_z_0_yyyy_1,  \
                             ta1_z_0_yyyz_0,  \
                             ta1_z_0_yyyz_1,  \
                             ta1_z_0_yyzz_0,  \
                             ta1_z_0_yyzz_1,  \
                             ta1_z_0_yzzz_0,  \
                             ta1_z_0_yzzz_1,  \
                             ta1_z_0_zzzz_0,  \
                             ta1_z_0_zzzz_1,  \
                             ta1_z_y_xxx_0,   \
                             ta1_z_y_xxx_1,   \
                             ta1_z_y_xxxx_0,  \
                             ta1_z_y_xxxx_1,  \
                             ta1_z_y_xxxy_0,  \
                             ta1_z_y_xxxy_1,  \
                             ta1_z_y_xxxz_0,  \
                             ta1_z_y_xxxz_1,  \
                             ta1_z_y_xxy_0,   \
                             ta1_z_y_xxy_1,   \
                             ta1_z_y_xxyy_0,  \
                             ta1_z_y_xxyy_1,  \
                             ta1_z_y_xxyz_0,  \
                             ta1_z_y_xxyz_1,  \
                             ta1_z_y_xxz_0,   \
                             ta1_z_y_xxz_1,   \
                             ta1_z_y_xxzz_0,  \
                             ta1_z_y_xxzz_1,  \
                             ta1_z_y_xyy_0,   \
                             ta1_z_y_xyy_1,   \
                             ta1_z_y_xyyy_0,  \
                             ta1_z_y_xyyy_1,  \
                             ta1_z_y_xyyz_0,  \
                             ta1_z_y_xyyz_1,  \
                             ta1_z_y_xyz_0,   \
                             ta1_z_y_xyz_1,   \
                             ta1_z_y_xyzz_0,  \
                             ta1_z_y_xyzz_1,  \
                             ta1_z_y_xzz_0,   \
                             ta1_z_y_xzz_1,   \
                             ta1_z_y_xzzz_0,  \
                             ta1_z_y_xzzz_1,  \
                             ta1_z_y_yyy_0,   \
                             ta1_z_y_yyy_1,   \
                             ta1_z_y_yyyy_0,  \
                             ta1_z_y_yyyy_1,  \
                             ta1_z_y_yyyz_0,  \
                             ta1_z_y_yyyz_1,  \
                             ta1_z_y_yyz_0,   \
                             ta1_z_y_yyz_1,   \
                             ta1_z_y_yyzz_0,  \
                             ta1_z_y_yyzz_1,  \
                             ta1_z_y_yzz_0,   \
                             ta1_z_y_yzz_1,   \
                             ta1_z_y_yzzz_0,  \
                             ta1_z_y_yzzz_1,  \
                             ta1_z_y_zzz_0,   \
                             ta1_z_y_zzz_1,   \
                             ta1_z_y_zzzz_0,  \
                             ta1_z_y_zzzz_1,  \
                             ta1_z_yy_xxxx_0, \
                             ta1_z_yy_xxxy_0, \
                             ta1_z_yy_xxxz_0, \
                             ta1_z_yy_xxyy_0, \
                             ta1_z_yy_xxyz_0, \
                             ta1_z_yy_xxzz_0, \
                             ta1_z_yy_xyyy_0, \
                             ta1_z_yy_xyyz_0, \
                             ta1_z_yy_xyzz_0, \
                             ta1_z_yy_xzzz_0, \
                             ta1_z_yy_yyyy_0, \
                             ta1_z_yy_yyyz_0, \
                             ta1_z_yy_yyzz_0, \
                             ta1_z_yy_yzzz_0, \
                             ta1_z_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yy_xxxx_0[i] = ta1_z_0_xxxx_0[i] * fe_0 - ta1_z_0_xxxx_1[i] * fe_0 + ta1_z_y_xxxx_0[i] * pa_y[i] - ta1_z_y_xxxx_1[i] * pc_y[i];

        ta1_z_yy_xxxy_0[i] = ta1_z_0_xxxy_0[i] * fe_0 - ta1_z_0_xxxy_1[i] * fe_0 + ta1_z_y_xxx_0[i] * fe_0 - ta1_z_y_xxx_1[i] * fe_0 +
                             ta1_z_y_xxxy_0[i] * pa_y[i] - ta1_z_y_xxxy_1[i] * pc_y[i];

        ta1_z_yy_xxxz_0[i] = ta1_z_0_xxxz_0[i] * fe_0 - ta1_z_0_xxxz_1[i] * fe_0 + ta1_z_y_xxxz_0[i] * pa_y[i] - ta1_z_y_xxxz_1[i] * pc_y[i];

        ta1_z_yy_xxyy_0[i] = ta1_z_0_xxyy_0[i] * fe_0 - ta1_z_0_xxyy_1[i] * fe_0 + 2.0 * ta1_z_y_xxy_0[i] * fe_0 - 2.0 * ta1_z_y_xxy_1[i] * fe_0 +
                             ta1_z_y_xxyy_0[i] * pa_y[i] - ta1_z_y_xxyy_1[i] * pc_y[i];

        ta1_z_yy_xxyz_0[i] = ta1_z_0_xxyz_0[i] * fe_0 - ta1_z_0_xxyz_1[i] * fe_0 + ta1_z_y_xxz_0[i] * fe_0 - ta1_z_y_xxz_1[i] * fe_0 +
                             ta1_z_y_xxyz_0[i] * pa_y[i] - ta1_z_y_xxyz_1[i] * pc_y[i];

        ta1_z_yy_xxzz_0[i] = ta1_z_0_xxzz_0[i] * fe_0 - ta1_z_0_xxzz_1[i] * fe_0 + ta1_z_y_xxzz_0[i] * pa_y[i] - ta1_z_y_xxzz_1[i] * pc_y[i];

        ta1_z_yy_xyyy_0[i] = ta1_z_0_xyyy_0[i] * fe_0 - ta1_z_0_xyyy_1[i] * fe_0 + 3.0 * ta1_z_y_xyy_0[i] * fe_0 - 3.0 * ta1_z_y_xyy_1[i] * fe_0 +
                             ta1_z_y_xyyy_0[i] * pa_y[i] - ta1_z_y_xyyy_1[i] * pc_y[i];

        ta1_z_yy_xyyz_0[i] = ta1_z_0_xyyz_0[i] * fe_0 - ta1_z_0_xyyz_1[i] * fe_0 + 2.0 * ta1_z_y_xyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyz_1[i] * fe_0 +
                             ta1_z_y_xyyz_0[i] * pa_y[i] - ta1_z_y_xyyz_1[i] * pc_y[i];

        ta1_z_yy_xyzz_0[i] = ta1_z_0_xyzz_0[i] * fe_0 - ta1_z_0_xyzz_1[i] * fe_0 + ta1_z_y_xzz_0[i] * fe_0 - ta1_z_y_xzz_1[i] * fe_0 +
                             ta1_z_y_xyzz_0[i] * pa_y[i] - ta1_z_y_xyzz_1[i] * pc_y[i];

        ta1_z_yy_xzzz_0[i] = ta1_z_0_xzzz_0[i] * fe_0 - ta1_z_0_xzzz_1[i] * fe_0 + ta1_z_y_xzzz_0[i] * pa_y[i] - ta1_z_y_xzzz_1[i] * pc_y[i];

        ta1_z_yy_yyyy_0[i] = ta1_z_0_yyyy_0[i] * fe_0 - ta1_z_0_yyyy_1[i] * fe_0 + 4.0 * ta1_z_y_yyy_0[i] * fe_0 - 4.0 * ta1_z_y_yyy_1[i] * fe_0 +
                             ta1_z_y_yyyy_0[i] * pa_y[i] - ta1_z_y_yyyy_1[i] * pc_y[i];

        ta1_z_yy_yyyz_0[i] = ta1_z_0_yyyz_0[i] * fe_0 - ta1_z_0_yyyz_1[i] * fe_0 + 3.0 * ta1_z_y_yyz_0[i] * fe_0 - 3.0 * ta1_z_y_yyz_1[i] * fe_0 +
                             ta1_z_y_yyyz_0[i] * pa_y[i] - ta1_z_y_yyyz_1[i] * pc_y[i];

        ta1_z_yy_yyzz_0[i] = ta1_z_0_yyzz_0[i] * fe_0 - ta1_z_0_yyzz_1[i] * fe_0 + 2.0 * ta1_z_y_yzz_0[i] * fe_0 - 2.0 * ta1_z_y_yzz_1[i] * fe_0 +
                             ta1_z_y_yyzz_0[i] * pa_y[i] - ta1_z_y_yyzz_1[i] * pc_y[i];

        ta1_z_yy_yzzz_0[i] = ta1_z_0_yzzz_0[i] * fe_0 - ta1_z_0_yzzz_1[i] * fe_0 + ta1_z_y_zzz_0[i] * fe_0 - ta1_z_y_zzz_1[i] * fe_0 +
                             ta1_z_y_yzzz_0[i] * pa_y[i] - ta1_z_y_yzzz_1[i] * pc_y[i];

        ta1_z_yy_zzzz_0[i] = ta1_z_0_zzzz_0[i] * fe_0 - ta1_z_0_zzzz_1[i] * fe_0 + ta1_z_y_zzzz_0[i] * pa_y[i] - ta1_z_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 240-255 components of targeted buffer : DG

    auto ta1_z_yz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 240);

    auto ta1_z_yz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 241);

    auto ta1_z_yz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 242);

    auto ta1_z_yz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 243);

    auto ta1_z_yz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 244);

    auto ta1_z_yz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 245);

    auto ta1_z_yz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 246);

    auto ta1_z_yz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 247);

    auto ta1_z_yz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 248);

    auto ta1_z_yz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 249);

    auto ta1_z_yz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 250);

    auto ta1_z_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 251);

    auto ta1_z_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 252);

    auto ta1_z_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 253);

    auto ta1_z_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 254);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_y_xxxy_0,  \
                             ta1_z_y_xxxy_1,  \
                             ta1_z_y_xxyy_0,  \
                             ta1_z_y_xxyy_1,  \
                             ta1_z_y_xyyy_0,  \
                             ta1_z_y_xyyy_1,  \
                             ta1_z_y_yyyy_0,  \
                             ta1_z_y_yyyy_1,  \
                             ta1_z_yz_xxxx_0, \
                             ta1_z_yz_xxxy_0, \
                             ta1_z_yz_xxxz_0, \
                             ta1_z_yz_xxyy_0, \
                             ta1_z_yz_xxyz_0, \
                             ta1_z_yz_xxzz_0, \
                             ta1_z_yz_xyyy_0, \
                             ta1_z_yz_xyyz_0, \
                             ta1_z_yz_xyzz_0, \
                             ta1_z_yz_xzzz_0, \
                             ta1_z_yz_yyyy_0, \
                             ta1_z_yz_yyyz_0, \
                             ta1_z_yz_yyzz_0, \
                             ta1_z_yz_yzzz_0, \
                             ta1_z_yz_zzzz_0, \
                             ta1_z_z_xxxx_0,  \
                             ta1_z_z_xxxx_1,  \
                             ta1_z_z_xxxz_0,  \
                             ta1_z_z_xxxz_1,  \
                             ta1_z_z_xxyz_0,  \
                             ta1_z_z_xxyz_1,  \
                             ta1_z_z_xxz_0,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xxzz_0,  \
                             ta1_z_z_xxzz_1,  \
                             ta1_z_z_xyyz_0,  \
                             ta1_z_z_xyyz_1,  \
                             ta1_z_z_xyz_0,   \
                             ta1_z_z_xyz_1,   \
                             ta1_z_z_xyzz_0,  \
                             ta1_z_z_xyzz_1,  \
                             ta1_z_z_xzz_0,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_xzzz_0,  \
                             ta1_z_z_xzzz_1,  \
                             ta1_z_z_yyyz_0,  \
                             ta1_z_z_yyyz_1,  \
                             ta1_z_z_yyz_0,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yyzz_0,  \
                             ta1_z_z_yyzz_1,  \
                             ta1_z_z_yzz_0,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_yzzz_0,  \
                             ta1_z_z_yzzz_1,  \
                             ta1_z_z_zzz_0,   \
                             ta1_z_z_zzz_1,   \
                             ta1_z_z_zzzz_0,  \
                             ta1_z_z_zzzz_1,  \
                             ta_y_xxxy_1,     \
                             ta_y_xxyy_1,     \
                             ta_y_xyyy_1,     \
                             ta_y_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yz_xxxx_0[i] = ta1_z_z_xxxx_0[i] * pa_y[i] - ta1_z_z_xxxx_1[i] * pc_y[i];

        ta1_z_yz_xxxy_0[i] = ta_y_xxxy_1[i] + ta1_z_y_xxxy_0[i] * pa_z[i] - ta1_z_y_xxxy_1[i] * pc_z[i];

        ta1_z_yz_xxxz_0[i] = ta1_z_z_xxxz_0[i] * pa_y[i] - ta1_z_z_xxxz_1[i] * pc_y[i];

        ta1_z_yz_xxyy_0[i] = ta_y_xxyy_1[i] + ta1_z_y_xxyy_0[i] * pa_z[i] - ta1_z_y_xxyy_1[i] * pc_z[i];

        ta1_z_yz_xxyz_0[i] = ta1_z_z_xxz_0[i] * fe_0 - ta1_z_z_xxz_1[i] * fe_0 + ta1_z_z_xxyz_0[i] * pa_y[i] - ta1_z_z_xxyz_1[i] * pc_y[i];

        ta1_z_yz_xxzz_0[i] = ta1_z_z_xxzz_0[i] * pa_y[i] - ta1_z_z_xxzz_1[i] * pc_y[i];

        ta1_z_yz_xyyy_0[i] = ta_y_xyyy_1[i] + ta1_z_y_xyyy_0[i] * pa_z[i] - ta1_z_y_xyyy_1[i] * pc_z[i];

        ta1_z_yz_xyyz_0[i] =
            2.0 * ta1_z_z_xyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyz_1[i] * fe_0 + ta1_z_z_xyyz_0[i] * pa_y[i] - ta1_z_z_xyyz_1[i] * pc_y[i];

        ta1_z_yz_xyzz_0[i] = ta1_z_z_xzz_0[i] * fe_0 - ta1_z_z_xzz_1[i] * fe_0 + ta1_z_z_xyzz_0[i] * pa_y[i] - ta1_z_z_xyzz_1[i] * pc_y[i];

        ta1_z_yz_xzzz_0[i] = ta1_z_z_xzzz_0[i] * pa_y[i] - ta1_z_z_xzzz_1[i] * pc_y[i];

        ta1_z_yz_yyyy_0[i] = ta_y_yyyy_1[i] + ta1_z_y_yyyy_0[i] * pa_z[i] - ta1_z_y_yyyy_1[i] * pc_z[i];

        ta1_z_yz_yyyz_0[i] =
            3.0 * ta1_z_z_yyz_0[i] * fe_0 - 3.0 * ta1_z_z_yyz_1[i] * fe_0 + ta1_z_z_yyyz_0[i] * pa_y[i] - ta1_z_z_yyyz_1[i] * pc_y[i];

        ta1_z_yz_yyzz_0[i] =
            2.0 * ta1_z_z_yzz_0[i] * fe_0 - 2.0 * ta1_z_z_yzz_1[i] * fe_0 + ta1_z_z_yyzz_0[i] * pa_y[i] - ta1_z_z_yyzz_1[i] * pc_y[i];

        ta1_z_yz_yzzz_0[i] = ta1_z_z_zzz_0[i] * fe_0 - ta1_z_z_zzz_1[i] * fe_0 + ta1_z_z_yzzz_0[i] * pa_y[i] - ta1_z_z_yzzz_1[i] * pc_y[i];

        ta1_z_yz_zzzz_0[i] = ta1_z_z_zzzz_0[i] * pa_y[i] - ta1_z_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : DG

    auto ta1_z_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 255);

    auto ta1_z_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 256);

    auto ta1_z_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 257);

    auto ta1_z_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 258);

    auto ta1_z_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 259);

    auto ta1_z_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 260);

    auto ta1_z_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 261);

    auto ta1_z_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 262);

    auto ta1_z_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 263);

    auto ta1_z_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 264);

    auto ta1_z_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 265);

    auto ta1_z_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 266);

    auto ta1_z_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 267);

    auto ta1_z_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 268);

    auto ta1_z_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 269);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_z_0_xxxx_0,  \
                             ta1_z_0_xxxx_1,  \
                             ta1_z_0_xxxy_0,  \
                             ta1_z_0_xxxy_1,  \
                             ta1_z_0_xxxz_0,  \
                             ta1_z_0_xxxz_1,  \
                             ta1_z_0_xxyy_0,  \
                             ta1_z_0_xxyy_1,  \
                             ta1_z_0_xxyz_0,  \
                             ta1_z_0_xxyz_1,  \
                             ta1_z_0_xxzz_0,  \
                             ta1_z_0_xxzz_1,  \
                             ta1_z_0_xyyy_0,  \
                             ta1_z_0_xyyy_1,  \
                             ta1_z_0_xyyz_0,  \
                             ta1_z_0_xyyz_1,  \
                             ta1_z_0_xyzz_0,  \
                             ta1_z_0_xyzz_1,  \
                             ta1_z_0_xzzz_0,  \
                             ta1_z_0_xzzz_1,  \
                             ta1_z_0_yyyy_0,  \
                             ta1_z_0_yyyy_1,  \
                             ta1_z_0_yyyz_0,  \
                             ta1_z_0_yyyz_1,  \
                             ta1_z_0_yyzz_0,  \
                             ta1_z_0_yyzz_1,  \
                             ta1_z_0_yzzz_0,  \
                             ta1_z_0_yzzz_1,  \
                             ta1_z_0_zzzz_0,  \
                             ta1_z_0_zzzz_1,  \
                             ta1_z_z_xxx_0,   \
                             ta1_z_z_xxx_1,   \
                             ta1_z_z_xxxx_0,  \
                             ta1_z_z_xxxx_1,  \
                             ta1_z_z_xxxy_0,  \
                             ta1_z_z_xxxy_1,  \
                             ta1_z_z_xxxz_0,  \
                             ta1_z_z_xxxz_1,  \
                             ta1_z_z_xxy_0,   \
                             ta1_z_z_xxy_1,   \
                             ta1_z_z_xxyy_0,  \
                             ta1_z_z_xxyy_1,  \
                             ta1_z_z_xxyz_0,  \
                             ta1_z_z_xxyz_1,  \
                             ta1_z_z_xxz_0,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xxzz_0,  \
                             ta1_z_z_xxzz_1,  \
                             ta1_z_z_xyy_0,   \
                             ta1_z_z_xyy_1,   \
                             ta1_z_z_xyyy_0,  \
                             ta1_z_z_xyyy_1,  \
                             ta1_z_z_xyyz_0,  \
                             ta1_z_z_xyyz_1,  \
                             ta1_z_z_xyz_0,   \
                             ta1_z_z_xyz_1,   \
                             ta1_z_z_xyzz_0,  \
                             ta1_z_z_xyzz_1,  \
                             ta1_z_z_xzz_0,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_xzzz_0,  \
                             ta1_z_z_xzzz_1,  \
                             ta1_z_z_yyy_0,   \
                             ta1_z_z_yyy_1,   \
                             ta1_z_z_yyyy_0,  \
                             ta1_z_z_yyyy_1,  \
                             ta1_z_z_yyyz_0,  \
                             ta1_z_z_yyyz_1,  \
                             ta1_z_z_yyz_0,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yyzz_0,  \
                             ta1_z_z_yyzz_1,  \
                             ta1_z_z_yzz_0,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_yzzz_0,  \
                             ta1_z_z_yzzz_1,  \
                             ta1_z_z_zzz_0,   \
                             ta1_z_z_zzz_1,   \
                             ta1_z_z_zzzz_0,  \
                             ta1_z_z_zzzz_1,  \
                             ta1_z_zz_xxxx_0, \
                             ta1_z_zz_xxxy_0, \
                             ta1_z_zz_xxxz_0, \
                             ta1_z_zz_xxyy_0, \
                             ta1_z_zz_xxyz_0, \
                             ta1_z_zz_xxzz_0, \
                             ta1_z_zz_xyyy_0, \
                             ta1_z_zz_xyyz_0, \
                             ta1_z_zz_xyzz_0, \
                             ta1_z_zz_xzzz_0, \
                             ta1_z_zz_yyyy_0, \
                             ta1_z_zz_yyyz_0, \
                             ta1_z_zz_yyzz_0, \
                             ta1_z_zz_yzzz_0, \
                             ta1_z_zz_zzzz_0, \
                             ta_z_xxxx_1,     \
                             ta_z_xxxy_1,     \
                             ta_z_xxxz_1,     \
                             ta_z_xxyy_1,     \
                             ta_z_xxyz_1,     \
                             ta_z_xxzz_1,     \
                             ta_z_xyyy_1,     \
                             ta_z_xyyz_1,     \
                             ta_z_xyzz_1,     \
                             ta_z_xzzz_1,     \
                             ta_z_yyyy_1,     \
                             ta_z_yyyz_1,     \
                             ta_z_yyzz_1,     \
                             ta_z_yzzz_1,     \
                             ta_z_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zz_xxxx_0[i] =
            ta1_z_0_xxxx_0[i] * fe_0 - ta1_z_0_xxxx_1[i] * fe_0 + ta_z_xxxx_1[i] + ta1_z_z_xxxx_0[i] * pa_z[i] - ta1_z_z_xxxx_1[i] * pc_z[i];

        ta1_z_zz_xxxy_0[i] =
            ta1_z_0_xxxy_0[i] * fe_0 - ta1_z_0_xxxy_1[i] * fe_0 + ta_z_xxxy_1[i] + ta1_z_z_xxxy_0[i] * pa_z[i] - ta1_z_z_xxxy_1[i] * pc_z[i];

        ta1_z_zz_xxxz_0[i] = ta1_z_0_xxxz_0[i] * fe_0 - ta1_z_0_xxxz_1[i] * fe_0 + ta1_z_z_xxx_0[i] * fe_0 - ta1_z_z_xxx_1[i] * fe_0 +
                             ta_z_xxxz_1[i] + ta1_z_z_xxxz_0[i] * pa_z[i] - ta1_z_z_xxxz_1[i] * pc_z[i];

        ta1_z_zz_xxyy_0[i] =
            ta1_z_0_xxyy_0[i] * fe_0 - ta1_z_0_xxyy_1[i] * fe_0 + ta_z_xxyy_1[i] + ta1_z_z_xxyy_0[i] * pa_z[i] - ta1_z_z_xxyy_1[i] * pc_z[i];

        ta1_z_zz_xxyz_0[i] = ta1_z_0_xxyz_0[i] * fe_0 - ta1_z_0_xxyz_1[i] * fe_0 + ta1_z_z_xxy_0[i] * fe_0 - ta1_z_z_xxy_1[i] * fe_0 +
                             ta_z_xxyz_1[i] + ta1_z_z_xxyz_0[i] * pa_z[i] - ta1_z_z_xxyz_1[i] * pc_z[i];

        ta1_z_zz_xxzz_0[i] = ta1_z_0_xxzz_0[i] * fe_0 - ta1_z_0_xxzz_1[i] * fe_0 + 2.0 * ta1_z_z_xxz_0[i] * fe_0 - 2.0 * ta1_z_z_xxz_1[i] * fe_0 +
                             ta_z_xxzz_1[i] + ta1_z_z_xxzz_0[i] * pa_z[i] - ta1_z_z_xxzz_1[i] * pc_z[i];

        ta1_z_zz_xyyy_0[i] =
            ta1_z_0_xyyy_0[i] * fe_0 - ta1_z_0_xyyy_1[i] * fe_0 + ta_z_xyyy_1[i] + ta1_z_z_xyyy_0[i] * pa_z[i] - ta1_z_z_xyyy_1[i] * pc_z[i];

        ta1_z_zz_xyyz_0[i] = ta1_z_0_xyyz_0[i] * fe_0 - ta1_z_0_xyyz_1[i] * fe_0 + ta1_z_z_xyy_0[i] * fe_0 - ta1_z_z_xyy_1[i] * fe_0 +
                             ta_z_xyyz_1[i] + ta1_z_z_xyyz_0[i] * pa_z[i] - ta1_z_z_xyyz_1[i] * pc_z[i];

        ta1_z_zz_xyzz_0[i] = ta1_z_0_xyzz_0[i] * fe_0 - ta1_z_0_xyzz_1[i] * fe_0 + 2.0 * ta1_z_z_xyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyz_1[i] * fe_0 +
                             ta_z_xyzz_1[i] + ta1_z_z_xyzz_0[i] * pa_z[i] - ta1_z_z_xyzz_1[i] * pc_z[i];

        ta1_z_zz_xzzz_0[i] = ta1_z_0_xzzz_0[i] * fe_0 - ta1_z_0_xzzz_1[i] * fe_0 + 3.0 * ta1_z_z_xzz_0[i] * fe_0 - 3.0 * ta1_z_z_xzz_1[i] * fe_0 +
                             ta_z_xzzz_1[i] + ta1_z_z_xzzz_0[i] * pa_z[i] - ta1_z_z_xzzz_1[i] * pc_z[i];

        ta1_z_zz_yyyy_0[i] =
            ta1_z_0_yyyy_0[i] * fe_0 - ta1_z_0_yyyy_1[i] * fe_0 + ta_z_yyyy_1[i] + ta1_z_z_yyyy_0[i] * pa_z[i] - ta1_z_z_yyyy_1[i] * pc_z[i];

        ta1_z_zz_yyyz_0[i] = ta1_z_0_yyyz_0[i] * fe_0 - ta1_z_0_yyyz_1[i] * fe_0 + ta1_z_z_yyy_0[i] * fe_0 - ta1_z_z_yyy_1[i] * fe_0 +
                             ta_z_yyyz_1[i] + ta1_z_z_yyyz_0[i] * pa_z[i] - ta1_z_z_yyyz_1[i] * pc_z[i];

        ta1_z_zz_yyzz_0[i] = ta1_z_0_yyzz_0[i] * fe_0 - ta1_z_0_yyzz_1[i] * fe_0 + 2.0 * ta1_z_z_yyz_0[i] * fe_0 - 2.0 * ta1_z_z_yyz_1[i] * fe_0 +
                             ta_z_yyzz_1[i] + ta1_z_z_yyzz_0[i] * pa_z[i] - ta1_z_z_yyzz_1[i] * pc_z[i];

        ta1_z_zz_yzzz_0[i] = ta1_z_0_yzzz_0[i] * fe_0 - ta1_z_0_yzzz_1[i] * fe_0 + 3.0 * ta1_z_z_yzz_0[i] * fe_0 - 3.0 * ta1_z_z_yzz_1[i] * fe_0 +
                             ta_z_yzzz_1[i] + ta1_z_z_yzzz_0[i] * pa_z[i] - ta1_z_z_yzzz_1[i] * pc_z[i];

        ta1_z_zz_zzzz_0[i] = ta1_z_0_zzzz_0[i] * fe_0 - ta1_z_0_zzzz_1[i] * fe_0 + 4.0 * ta1_z_z_zzz_0[i] * fe_0 - 4.0 * ta1_z_z_zzz_1[i] * fe_0 +
                             ta_z_zzzz_1[i] + ta1_z_z_zzzz_0[i] * pa_z[i] - ta1_z_z_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
