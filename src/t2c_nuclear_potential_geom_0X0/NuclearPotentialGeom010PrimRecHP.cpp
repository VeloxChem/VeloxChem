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

#include "NuclearPotentialGeom010PrimRecHP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_hp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_hp,
                                        const size_t              idx_npot_geom_010_0_fp,
                                        const size_t              idx_npot_geom_010_1_fp,
                                        const size_t              idx_npot_geom_010_0_gs,
                                        const size_t              idx_npot_geom_010_1_gs,
                                        const size_t              idx_npot_1_gp,
                                        const size_t              idx_npot_geom_010_0_gp,
                                        const size_t              idx_npot_geom_010_1_gp,
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

    // Set up components of auxiliary buffer : FP

    auto ta1_x_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp);

    auto ta1_x_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 1);

    auto ta1_x_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 2);

    auto ta1_x_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 3);

    auto ta1_x_xxy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 5);

    auto ta1_x_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 6);

    auto ta1_x_xxz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 7);

    auto ta1_x_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 8);

    auto ta1_x_xyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 9);

    auto ta1_x_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 10);

    auto ta1_x_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 15);

    auto ta1_x_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 17);

    auto ta1_x_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 18);

    auto ta1_x_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 19);

    auto ta1_x_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 20);

    auto ta1_x_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 22);

    auto ta1_x_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 23);

    auto ta1_x_yzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 24);

    auto ta1_x_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 26);

    auto ta1_x_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 27);

    auto ta1_x_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 28);

    auto ta1_x_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 29);

    auto ta1_y_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 30);

    auto ta1_y_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 31);

    auto ta1_y_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 32);

    auto ta1_y_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 33);

    auto ta1_y_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 34);

    auto ta1_y_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 36);

    auto ta1_y_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 38);

    auto ta1_y_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 40);

    auto ta1_y_xyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 41);

    auto ta1_y_xzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 46);

    auto ta1_y_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 47);

    auto ta1_y_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 48);

    auto ta1_y_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 49);

    auto ta1_y_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 50);

    auto ta1_y_yyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 51);

    auto ta1_y_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 52);

    auto ta1_y_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 53);

    auto ta1_y_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 55);

    auto ta1_y_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 56);

    auto ta1_y_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 57);

    auto ta1_y_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 58);

    auto ta1_y_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 59);

    auto ta1_z_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 60);

    auto ta1_z_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 61);

    auto ta1_z_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 62);

    auto ta1_z_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 63);

    auto ta1_z_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 64);

    auto ta1_z_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 66);

    auto ta1_z_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 68);

    auto ta1_z_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 70);

    auto ta1_z_xyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 71);

    auto ta1_z_xzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 76);

    auto ta1_z_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 77);

    auto ta1_z_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 78);

    auto ta1_z_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 79);

    auto ta1_z_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 80);

    auto ta1_z_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 82);

    auto ta1_z_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 83);

    auto ta1_z_yzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 84);

    auto ta1_z_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 85);

    auto ta1_z_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 86);

    auto ta1_z_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 87);

    auto ta1_z_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 88);

    auto ta1_z_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 89);

    // Set up components of auxiliary buffer : FP

    auto ta1_x_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp);

    auto ta1_x_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 1);

    auto ta1_x_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 2);

    auto ta1_x_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 3);

    auto ta1_x_xxy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 5);

    auto ta1_x_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 6);

    auto ta1_x_xxz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 7);

    auto ta1_x_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 8);

    auto ta1_x_xyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 9);

    auto ta1_x_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 10);

    auto ta1_x_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 15);

    auto ta1_x_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 17);

    auto ta1_x_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 18);

    auto ta1_x_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 19);

    auto ta1_x_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 20);

    auto ta1_x_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 22);

    auto ta1_x_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 23);

    auto ta1_x_yzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 24);

    auto ta1_x_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 26);

    auto ta1_x_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 27);

    auto ta1_x_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 28);

    auto ta1_x_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 29);

    auto ta1_y_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 30);

    auto ta1_y_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 31);

    auto ta1_y_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 32);

    auto ta1_y_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 33);

    auto ta1_y_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 34);

    auto ta1_y_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 36);

    auto ta1_y_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 38);

    auto ta1_y_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 40);

    auto ta1_y_xyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 41);

    auto ta1_y_xzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 46);

    auto ta1_y_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 47);

    auto ta1_y_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 48);

    auto ta1_y_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 49);

    auto ta1_y_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 50);

    auto ta1_y_yyz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 51);

    auto ta1_y_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 52);

    auto ta1_y_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 53);

    auto ta1_y_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 55);

    auto ta1_y_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 56);

    auto ta1_y_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 57);

    auto ta1_y_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 58);

    auto ta1_y_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 59);

    auto ta1_z_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 60);

    auto ta1_z_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 61);

    auto ta1_z_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 62);

    auto ta1_z_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 63);

    auto ta1_z_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 64);

    auto ta1_z_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 66);

    auto ta1_z_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 68);

    auto ta1_z_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 70);

    auto ta1_z_xyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 71);

    auto ta1_z_xzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 76);

    auto ta1_z_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 77);

    auto ta1_z_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 78);

    auto ta1_z_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 79);

    auto ta1_z_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 80);

    auto ta1_z_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 82);

    auto ta1_z_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 83);

    auto ta1_z_yzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 84);

    auto ta1_z_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 85);

    auto ta1_z_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 86);

    auto ta1_z_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 87);

    auto ta1_z_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 88);

    auto ta1_z_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 89);

    // Set up components of auxiliary buffer : GS

    auto ta1_x_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs);

    auto ta1_x_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 5);

    auto ta1_x_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 10);

    auto ta1_x_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 14);

    auto ta1_y_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 15);

    auto ta1_y_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 25);

    auto ta1_y_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 27);

    auto ta1_y_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 29);

    auto ta1_z_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 30);

    auto ta1_z_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 40);

    auto ta1_z_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 42);

    auto ta1_z_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 44);

    // Set up components of auxiliary buffer : GS

    auto ta1_x_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs);

    auto ta1_x_xxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 5);

    auto ta1_x_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 10);

    auto ta1_x_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 14);

    auto ta1_y_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 15);

    auto ta1_y_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 25);

    auto ta1_y_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 27);

    auto ta1_y_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 29);

    auto ta1_z_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 30);

    auto ta1_z_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 40);

    auto ta1_z_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 42);

    auto ta1_z_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 44);

    // Set up components of auxiliary buffer : GP

    auto ta_xxxx_x_1 = pbuffer.data(idx_npot_1_gp);

    auto ta_xxxx_y_1 = pbuffer.data(idx_npot_1_gp + 1);

    auto ta_xxxx_z_1 = pbuffer.data(idx_npot_1_gp + 2);

    auto ta_xxxy_x_1 = pbuffer.data(idx_npot_1_gp + 3);

    auto ta_xxxy_y_1 = pbuffer.data(idx_npot_1_gp + 4);

    auto ta_xxxz_x_1 = pbuffer.data(idx_npot_1_gp + 6);

    auto ta_xxxz_z_1 = pbuffer.data(idx_npot_1_gp + 8);

    auto ta_xxyy_x_1 = pbuffer.data(idx_npot_1_gp + 9);

    auto ta_xxyy_y_1 = pbuffer.data(idx_npot_1_gp + 10);

    auto ta_xxzz_x_1 = pbuffer.data(idx_npot_1_gp + 15);

    auto ta_xxzz_z_1 = pbuffer.data(idx_npot_1_gp + 17);

    auto ta_xyyy_x_1 = pbuffer.data(idx_npot_1_gp + 18);

    auto ta_xyyy_y_1 = pbuffer.data(idx_npot_1_gp + 19);

    auto ta_xzzz_x_1 = pbuffer.data(idx_npot_1_gp + 27);

    auto ta_xzzz_z_1 = pbuffer.data(idx_npot_1_gp + 29);

    auto ta_yyyy_x_1 = pbuffer.data(idx_npot_1_gp + 30);

    auto ta_yyyy_y_1 = pbuffer.data(idx_npot_1_gp + 31);

    auto ta_yyyy_z_1 = pbuffer.data(idx_npot_1_gp + 32);

    auto ta_yyyz_y_1 = pbuffer.data(idx_npot_1_gp + 34);

    auto ta_yyyz_z_1 = pbuffer.data(idx_npot_1_gp + 35);

    auto ta_yyzz_y_1 = pbuffer.data(idx_npot_1_gp + 37);

    auto ta_yyzz_z_1 = pbuffer.data(idx_npot_1_gp + 38);

    auto ta_yzzz_y_1 = pbuffer.data(idx_npot_1_gp + 40);

    auto ta_yzzz_z_1 = pbuffer.data(idx_npot_1_gp + 41);

    auto ta_zzzz_x_1 = pbuffer.data(idx_npot_1_gp + 42);

    auto ta_zzzz_y_1 = pbuffer.data(idx_npot_1_gp + 43);

    auto ta_zzzz_z_1 = pbuffer.data(idx_npot_1_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto ta1_x_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp);

    auto ta1_x_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 1);

    auto ta1_x_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 2);

    auto ta1_x_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 3);

    auto ta1_x_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 4);

    auto ta1_x_xxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 5);

    auto ta1_x_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 6);

    auto ta1_x_xxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 7);

    auto ta1_x_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 8);

    auto ta1_x_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 9);

    auto ta1_x_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 10);

    auto ta1_x_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 11);

    auto ta1_x_xxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 14);

    auto ta1_x_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 15);

    auto ta1_x_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 16);

    auto ta1_x_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 17);

    auto ta1_x_xyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 18);

    auto ta1_x_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 19);

    auto ta1_x_xyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 24);

    auto ta1_x_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 27);

    auto ta1_x_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 29);

    auto ta1_x_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 30);

    auto ta1_x_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 31);

    auto ta1_x_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 32);

    auto ta1_x_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 34);

    auto ta1_x_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 35);

    auto ta1_x_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 36);

    auto ta1_x_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 37);

    auto ta1_x_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 38);

    auto ta1_x_yzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 39);

    auto ta1_x_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 40);

    auto ta1_x_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 41);

    auto ta1_x_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 42);

    auto ta1_x_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 43);

    auto ta1_x_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 44);

    auto ta1_y_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 45);

    auto ta1_y_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 46);

    auto ta1_y_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 47);

    auto ta1_y_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 48);

    auto ta1_y_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 49);

    auto ta1_y_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 51);

    auto ta1_y_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 53);

    auto ta1_y_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 54);

    auto ta1_y_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 55);

    auto ta1_y_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 56);

    auto ta1_y_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 60);

    auto ta1_y_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 61);

    auto ta1_y_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 62);

    auto ta1_y_xyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 63);

    auto ta1_y_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 64);

    auto ta1_y_xyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 65);

    auto ta1_y_xyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 68);

    auto ta1_y_xyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 70);

    auto ta1_y_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 72);

    auto ta1_y_xzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 73);

    auto ta1_y_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 74);

    auto ta1_y_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 75);

    auto ta1_y_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 76);

    auto ta1_y_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 77);

    auto ta1_y_yyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 78);

    auto ta1_y_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 79);

    auto ta1_y_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 80);

    auto ta1_y_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 81);

    auto ta1_y_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 82);

    auto ta1_y_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 83);

    auto ta1_y_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 85);

    auto ta1_y_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 86);

    auto ta1_y_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 87);

    auto ta1_y_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 88);

    auto ta1_y_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 89);

    auto ta1_z_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 90);

    auto ta1_z_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 91);

    auto ta1_z_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 92);

    auto ta1_z_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 93);

    auto ta1_z_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 94);

    auto ta1_z_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 96);

    auto ta1_z_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 98);

    auto ta1_z_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 99);

    auto ta1_z_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 100);

    auto ta1_z_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 101);

    auto ta1_z_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 105);

    auto ta1_z_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 106);

    auto ta1_z_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 107);

    auto ta1_z_xyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 108);

    auto ta1_z_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 109);

    auto ta1_z_xyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 110);

    auto ta1_z_xyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 113);

    auto ta1_z_xyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 115);

    auto ta1_z_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 117);

    auto ta1_z_xzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 118);

    auto ta1_z_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 119);

    auto ta1_z_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 120);

    auto ta1_z_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 121);

    auto ta1_z_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 122);

    auto ta1_z_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 124);

    auto ta1_z_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 125);

    auto ta1_z_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 126);

    auto ta1_z_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 127);

    auto ta1_z_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 128);

    auto ta1_z_yzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 129);

    auto ta1_z_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 130);

    auto ta1_z_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 131);

    auto ta1_z_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 132);

    auto ta1_z_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 133);

    auto ta1_z_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 134);

    // Set up components of auxiliary buffer : GP

    auto ta1_x_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp);

    auto ta1_x_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 1);

    auto ta1_x_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 2);

    auto ta1_x_xxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 3);

    auto ta1_x_xxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 4);

    auto ta1_x_xxxy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 5);

    auto ta1_x_xxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 6);

    auto ta1_x_xxxz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 7);

    auto ta1_x_xxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 8);

    auto ta1_x_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 9);

    auto ta1_x_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 10);

    auto ta1_x_xxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 11);

    auto ta1_x_xxyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 14);

    auto ta1_x_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 15);

    auto ta1_x_xxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 16);

    auto ta1_x_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 17);

    auto ta1_x_xyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 18);

    auto ta1_x_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 19);

    auto ta1_x_xyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 24);

    auto ta1_x_xzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 27);

    auto ta1_x_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 29);

    auto ta1_x_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 30);

    auto ta1_x_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 31);

    auto ta1_x_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 32);

    auto ta1_x_yyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 34);

    auto ta1_x_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 35);

    auto ta1_x_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 36);

    auto ta1_x_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 37);

    auto ta1_x_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 38);

    auto ta1_x_yzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 39);

    auto ta1_x_yzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 40);

    auto ta1_x_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 41);

    auto ta1_x_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 42);

    auto ta1_x_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 43);

    auto ta1_x_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 44);

    auto ta1_y_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 45);

    auto ta1_y_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 46);

    auto ta1_y_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 47);

    auto ta1_y_xxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 48);

    auto ta1_y_xxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 49);

    auto ta1_y_xxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 51);

    auto ta1_y_xxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 53);

    auto ta1_y_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 54);

    auto ta1_y_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 55);

    auto ta1_y_xxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 56);

    auto ta1_y_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 60);

    auto ta1_y_xxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 61);

    auto ta1_y_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 62);

    auto ta1_y_xyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 63);

    auto ta1_y_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 64);

    auto ta1_y_xyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 65);

    auto ta1_y_xyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 68);

    auto ta1_y_xyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 70);

    auto ta1_y_xzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 72);

    auto ta1_y_xzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 73);

    auto ta1_y_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 74);

    auto ta1_y_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 75);

    auto ta1_y_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 76);

    auto ta1_y_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 77);

    auto ta1_y_yyyz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 78);

    auto ta1_y_yyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 79);

    auto ta1_y_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 80);

    auto ta1_y_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 81);

    auto ta1_y_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 82);

    auto ta1_y_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 83);

    auto ta1_y_yzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 85);

    auto ta1_y_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 86);

    auto ta1_y_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 87);

    auto ta1_y_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 88);

    auto ta1_y_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 89);

    auto ta1_z_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 90);

    auto ta1_z_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 91);

    auto ta1_z_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 92);

    auto ta1_z_xxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 93);

    auto ta1_z_xxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 94);

    auto ta1_z_xxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 96);

    auto ta1_z_xxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 98);

    auto ta1_z_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 99);

    auto ta1_z_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 100);

    auto ta1_z_xxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 101);

    auto ta1_z_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 105);

    auto ta1_z_xxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 106);

    auto ta1_z_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 107);

    auto ta1_z_xyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 108);

    auto ta1_z_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 109);

    auto ta1_z_xyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 110);

    auto ta1_z_xyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 113);

    auto ta1_z_xyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 115);

    auto ta1_z_xzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 117);

    auto ta1_z_xzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 118);

    auto ta1_z_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 119);

    auto ta1_z_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 120);

    auto ta1_z_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 121);

    auto ta1_z_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 122);

    auto ta1_z_yyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 124);

    auto ta1_z_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 125);

    auto ta1_z_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 126);

    auto ta1_z_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 127);

    auto ta1_z_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 128);

    auto ta1_z_yzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 129);

    auto ta1_z_yzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 130);

    auto ta1_z_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 131);

    auto ta1_z_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 132);

    auto ta1_z_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 133);

    auto ta1_z_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 134);

    // Set up 0-3 components of targeted buffer : HP

    auto ta1_x_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp);

    auto ta1_x_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 1);

    auto ta1_x_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 2);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_x_xxx_x_0,   \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_y_0,   \
                             ta1_x_xxx_y_1,   \
                             ta1_x_xxx_z_0,   \
                             ta1_x_xxx_z_1,   \
                             ta1_x_xxxx_0_0,  \
                             ta1_x_xxxx_0_1,  \
                             ta1_x_xxxx_x_0,  \
                             ta1_x_xxxx_x_1,  \
                             ta1_x_xxxx_y_0,  \
                             ta1_x_xxxx_y_1,  \
                             ta1_x_xxxx_z_0,  \
                             ta1_x_xxxx_z_1,  \
                             ta1_x_xxxxx_x_0, \
                             ta1_x_xxxxx_y_0, \
                             ta1_x_xxxxx_z_0, \
                             ta_xxxx_x_1,     \
                             ta_xxxx_y_1,     \
                             ta_xxxx_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxx_x_0[i] = 4.0 * ta1_x_xxx_x_0[i] * fe_0 - 4.0 * ta1_x_xxx_x_1[i] * fe_0 + ta1_x_xxxx_0_0[i] * fe_0 - ta1_x_xxxx_0_1[i] * fe_0 +
                             ta_xxxx_x_1[i] + ta1_x_xxxx_x_0[i] * pa_x[i] - ta1_x_xxxx_x_1[i] * pc_x[i];

        ta1_x_xxxxx_y_0[i] = 4.0 * ta1_x_xxx_y_0[i] * fe_0 - 4.0 * ta1_x_xxx_y_1[i] * fe_0 + ta_xxxx_y_1[i] + ta1_x_xxxx_y_0[i] * pa_x[i] -
                             ta1_x_xxxx_y_1[i] * pc_x[i];

        ta1_x_xxxxx_z_0[i] = 4.0 * ta1_x_xxx_z_0[i] * fe_0 - 4.0 * ta1_x_xxx_z_1[i] * fe_0 + ta_xxxx_z_1[i] + ta1_x_xxxx_z_0[i] * pa_x[i] -
                             ta1_x_xxxx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : HP

    auto ta1_x_xxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 3);

    auto ta1_x_xxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 4);

    auto ta1_x_xxxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 5);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_xxxx_0_0,  \
                             ta1_x_xxxx_0_1,  \
                             ta1_x_xxxx_x_0,  \
                             ta1_x_xxxx_x_1,  \
                             ta1_x_xxxx_y_0,  \
                             ta1_x_xxxx_y_1,  \
                             ta1_x_xxxx_z_0,  \
                             ta1_x_xxxx_z_1,  \
                             ta1_x_xxxxy_x_0, \
                             ta1_x_xxxxy_y_0, \
                             ta1_x_xxxxy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxy_x_0[i] = ta1_x_xxxx_x_0[i] * pa_y[i] - ta1_x_xxxx_x_1[i] * pc_y[i];

        ta1_x_xxxxy_y_0[i] = ta1_x_xxxx_0_0[i] * fe_0 - ta1_x_xxxx_0_1[i] * fe_0 + ta1_x_xxxx_y_0[i] * pa_y[i] - ta1_x_xxxx_y_1[i] * pc_y[i];

        ta1_x_xxxxy_z_0[i] = ta1_x_xxxx_z_0[i] * pa_y[i] - ta1_x_xxxx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : HP

    auto ta1_x_xxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 6);

    auto ta1_x_xxxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 7);

    auto ta1_x_xxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 8);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_xxxx_0_0,  \
                             ta1_x_xxxx_0_1,  \
                             ta1_x_xxxx_x_0,  \
                             ta1_x_xxxx_x_1,  \
                             ta1_x_xxxx_y_0,  \
                             ta1_x_xxxx_y_1,  \
                             ta1_x_xxxx_z_0,  \
                             ta1_x_xxxx_z_1,  \
                             ta1_x_xxxxz_x_0, \
                             ta1_x_xxxxz_y_0, \
                             ta1_x_xxxxz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxz_x_0[i] = ta1_x_xxxx_x_0[i] * pa_z[i] - ta1_x_xxxx_x_1[i] * pc_z[i];

        ta1_x_xxxxz_y_0[i] = ta1_x_xxxx_y_0[i] * pa_z[i] - ta1_x_xxxx_y_1[i] * pc_z[i];

        ta1_x_xxxxz_z_0[i] = ta1_x_xxxx_0_0[i] * fe_0 - ta1_x_xxxx_0_1[i] * fe_0 + ta1_x_xxxx_z_0[i] * pa_z[i] - ta1_x_xxxx_z_1[i] * pc_z[i];
    }

    // Set up 9-12 components of targeted buffer : HP

    auto ta1_x_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 9);

    auto ta1_x_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 10);

    auto ta1_x_xxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 11);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xxx_x_0,   \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_z_0,   \
                             ta1_x_xxx_z_1,   \
                             ta1_x_xxxy_x_0,  \
                             ta1_x_xxxy_x_1,  \
                             ta1_x_xxxy_z_0,  \
                             ta1_x_xxxy_z_1,  \
                             ta1_x_xxxyy_x_0, \
                             ta1_x_xxxyy_y_0, \
                             ta1_x_xxxyy_z_0, \
                             ta1_x_xxyy_y_0,  \
                             ta1_x_xxyy_y_1,  \
                             ta1_x_xyy_y_0,   \
                             ta1_x_xyy_y_1,   \
                             ta_xxyy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyy_x_0[i] = ta1_x_xxx_x_0[i] * fe_0 - ta1_x_xxx_x_1[i] * fe_0 + ta1_x_xxxy_x_0[i] * pa_y[i] - ta1_x_xxxy_x_1[i] * pc_y[i];

        ta1_x_xxxyy_y_0[i] = 2.0 * ta1_x_xyy_y_0[i] * fe_0 - 2.0 * ta1_x_xyy_y_1[i] * fe_0 + ta_xxyy_y_1[i] + ta1_x_xxyy_y_0[i] * pa_x[i] -
                             ta1_x_xxyy_y_1[i] * pc_x[i];

        ta1_x_xxxyy_z_0[i] = ta1_x_xxx_z_0[i] * fe_0 - ta1_x_xxx_z_1[i] * fe_0 + ta1_x_xxxy_z_0[i] * pa_y[i] - ta1_x_xxxy_z_1[i] * pc_y[i];
    }

    // Set up 12-15 components of targeted buffer : HP

    auto ta1_x_xxxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 12);

    auto ta1_x_xxxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 13);

    auto ta1_x_xxxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 14);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xxxy_y_0,  \
                             ta1_x_xxxy_y_1,  \
                             ta1_x_xxxyz_x_0, \
                             ta1_x_xxxyz_y_0, \
                             ta1_x_xxxyz_z_0, \
                             ta1_x_xxxz_x_0,  \
                             ta1_x_xxxz_x_1,  \
                             ta1_x_xxxz_z_0,  \
                             ta1_x_xxxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xxxyz_x_0[i] = ta1_x_xxxz_x_0[i] * pa_y[i] - ta1_x_xxxz_x_1[i] * pc_y[i];

        ta1_x_xxxyz_y_0[i] = ta1_x_xxxy_y_0[i] * pa_z[i] - ta1_x_xxxy_y_1[i] * pc_z[i];

        ta1_x_xxxyz_z_0[i] = ta1_x_xxxz_z_0[i] * pa_y[i] - ta1_x_xxxz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : HP

    auto ta1_x_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 15);

    auto ta1_x_xxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 16);

    auto ta1_x_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 17);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xxx_x_0,   \
                             ta1_x_xxx_x_1,   \
                             ta1_x_xxx_y_0,   \
                             ta1_x_xxx_y_1,   \
                             ta1_x_xxxz_x_0,  \
                             ta1_x_xxxz_x_1,  \
                             ta1_x_xxxz_y_0,  \
                             ta1_x_xxxz_y_1,  \
                             ta1_x_xxxzz_x_0, \
                             ta1_x_xxxzz_y_0, \
                             ta1_x_xxxzz_z_0, \
                             ta1_x_xxzz_z_0,  \
                             ta1_x_xxzz_z_1,  \
                             ta1_x_xzz_z_0,   \
                             ta1_x_xzz_z_1,   \
                             ta_xxzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzz_x_0[i] = ta1_x_xxx_x_0[i] * fe_0 - ta1_x_xxx_x_1[i] * fe_0 + ta1_x_xxxz_x_0[i] * pa_z[i] - ta1_x_xxxz_x_1[i] * pc_z[i];

        ta1_x_xxxzz_y_0[i] = ta1_x_xxx_y_0[i] * fe_0 - ta1_x_xxx_y_1[i] * fe_0 + ta1_x_xxxz_y_0[i] * pa_z[i] - ta1_x_xxxz_y_1[i] * pc_z[i];

        ta1_x_xxxzz_z_0[i] = 2.0 * ta1_x_xzz_z_0[i] * fe_0 - 2.0 * ta1_x_xzz_z_1[i] * fe_0 + ta_xxzz_z_1[i] + ta1_x_xxzz_z_0[i] * pa_x[i] -
                             ta1_x_xxzz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : HP

    auto ta1_x_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 18);

    auto ta1_x_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 19);

    auto ta1_x_xxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 20);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xxy_x_0,   \
                             ta1_x_xxy_x_1,   \
                             ta1_x_xxy_z_0,   \
                             ta1_x_xxy_z_1,   \
                             ta1_x_xxyy_x_0,  \
                             ta1_x_xxyy_x_1,  \
                             ta1_x_xxyy_z_0,  \
                             ta1_x_xxyy_z_1,  \
                             ta1_x_xxyyy_x_0, \
                             ta1_x_xxyyy_y_0, \
                             ta1_x_xxyyy_z_0, \
                             ta1_x_xyyy_y_0,  \
                             ta1_x_xyyy_y_1,  \
                             ta1_x_yyy_y_0,   \
                             ta1_x_yyy_y_1,   \
                             ta_xyyy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyy_x_0[i] =
            2.0 * ta1_x_xxy_x_0[i] * fe_0 - 2.0 * ta1_x_xxy_x_1[i] * fe_0 + ta1_x_xxyy_x_0[i] * pa_y[i] - ta1_x_xxyy_x_1[i] * pc_y[i];

        ta1_x_xxyyy_y_0[i] =
            ta1_x_yyy_y_0[i] * fe_0 - ta1_x_yyy_y_1[i] * fe_0 + ta_xyyy_y_1[i] + ta1_x_xyyy_y_0[i] * pa_x[i] - ta1_x_xyyy_y_1[i] * pc_x[i];

        ta1_x_xxyyy_z_0[i] =
            2.0 * ta1_x_xxy_z_0[i] * fe_0 - 2.0 * ta1_x_xxy_z_1[i] * fe_0 + ta1_x_xxyy_z_0[i] * pa_y[i] - ta1_x_xxyy_z_1[i] * pc_y[i];
    }

    // Set up 21-24 components of targeted buffer : HP

    auto ta1_x_xxyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 21);

    auto ta1_x_xxyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 22);

    auto ta1_x_xxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 23);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xxyy_x_0,  \
                             ta1_x_xxyy_x_1,  \
                             ta1_x_xxyy_y_0,  \
                             ta1_x_xxyy_y_1,  \
                             ta1_x_xxyyz_x_0, \
                             ta1_x_xxyyz_y_0, \
                             ta1_x_xxyyz_z_0, \
                             ta1_x_xxyz_z_0,  \
                             ta1_x_xxyz_z_1,  \
                             ta1_x_xxz_z_0,   \
                             ta1_x_xxz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyz_x_0[i] = ta1_x_xxyy_x_0[i] * pa_z[i] - ta1_x_xxyy_x_1[i] * pc_z[i];

        ta1_x_xxyyz_y_0[i] = ta1_x_xxyy_y_0[i] * pa_z[i] - ta1_x_xxyy_y_1[i] * pc_z[i];

        ta1_x_xxyyz_z_0[i] = ta1_x_xxz_z_0[i] * fe_0 - ta1_x_xxz_z_1[i] * fe_0 + ta1_x_xxyz_z_0[i] * pa_y[i] - ta1_x_xxyz_z_1[i] * pc_y[i];
    }

    // Set up 24-27 components of targeted buffer : HP

    auto ta1_x_xxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 24);

    auto ta1_x_xxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 25);

    auto ta1_x_xxyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 26);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_xxyzz_x_0, \
                             ta1_x_xxyzz_y_0, \
                             ta1_x_xxyzz_z_0, \
                             ta1_x_xxzz_0_0,  \
                             ta1_x_xxzz_0_1,  \
                             ta1_x_xxzz_x_0,  \
                             ta1_x_xxzz_x_1,  \
                             ta1_x_xxzz_y_0,  \
                             ta1_x_xxzz_y_1,  \
                             ta1_x_xxzz_z_0,  \
                             ta1_x_xxzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzz_x_0[i] = ta1_x_xxzz_x_0[i] * pa_y[i] - ta1_x_xxzz_x_1[i] * pc_y[i];

        ta1_x_xxyzz_y_0[i] = ta1_x_xxzz_0_0[i] * fe_0 - ta1_x_xxzz_0_1[i] * fe_0 + ta1_x_xxzz_y_0[i] * pa_y[i] - ta1_x_xxzz_y_1[i] * pc_y[i];

        ta1_x_xxyzz_z_0[i] = ta1_x_xxzz_z_0[i] * pa_y[i] - ta1_x_xxzz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : HP

    auto ta1_x_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 27);

    auto ta1_x_xxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 28);

    auto ta1_x_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 29);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xxz_x_0,   \
                             ta1_x_xxz_x_1,   \
                             ta1_x_xxz_y_0,   \
                             ta1_x_xxz_y_1,   \
                             ta1_x_xxzz_x_0,  \
                             ta1_x_xxzz_x_1,  \
                             ta1_x_xxzz_y_0,  \
                             ta1_x_xxzz_y_1,  \
                             ta1_x_xxzzz_x_0, \
                             ta1_x_xxzzz_y_0, \
                             ta1_x_xxzzz_z_0, \
                             ta1_x_xzzz_z_0,  \
                             ta1_x_xzzz_z_1,  \
                             ta1_x_zzz_z_0,   \
                             ta1_x_zzz_z_1,   \
                             ta_xzzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzz_x_0[i] =
            2.0 * ta1_x_xxz_x_0[i] * fe_0 - 2.0 * ta1_x_xxz_x_1[i] * fe_0 + ta1_x_xxzz_x_0[i] * pa_z[i] - ta1_x_xxzz_x_1[i] * pc_z[i];

        ta1_x_xxzzz_y_0[i] =
            2.0 * ta1_x_xxz_y_0[i] * fe_0 - 2.0 * ta1_x_xxz_y_1[i] * fe_0 + ta1_x_xxzz_y_0[i] * pa_z[i] - ta1_x_xxzz_y_1[i] * pc_z[i];

        ta1_x_xxzzz_z_0[i] =
            ta1_x_zzz_z_0[i] * fe_0 - ta1_x_zzz_z_1[i] * fe_0 + ta_xzzz_z_1[i] + ta1_x_xzzz_z_0[i] * pa_x[i] - ta1_x_xzzz_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : HP

    auto ta1_x_xyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 30);

    auto ta1_x_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 31);

    auto ta1_x_xyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 32);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xyy_x_0,   \
                             ta1_x_xyy_x_1,   \
                             ta1_x_xyyy_x_0,  \
                             ta1_x_xyyy_x_1,  \
                             ta1_x_xyyyy_x_0, \
                             ta1_x_xyyyy_y_0, \
                             ta1_x_xyyyy_z_0, \
                             ta1_x_yyyy_y_0,  \
                             ta1_x_yyyy_y_1,  \
                             ta1_x_yyyy_z_0,  \
                             ta1_x_yyyy_z_1,  \
                             ta_yyyy_y_1,     \
                             ta_yyyy_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyy_x_0[i] =
            3.0 * ta1_x_xyy_x_0[i] * fe_0 - 3.0 * ta1_x_xyy_x_1[i] * fe_0 + ta1_x_xyyy_x_0[i] * pa_y[i] - ta1_x_xyyy_x_1[i] * pc_y[i];

        ta1_x_xyyyy_y_0[i] = ta_yyyy_y_1[i] + ta1_x_yyyy_y_0[i] * pa_x[i] - ta1_x_yyyy_y_1[i] * pc_x[i];

        ta1_x_xyyyy_z_0[i] = ta_yyyy_z_1[i] + ta1_x_yyyy_z_0[i] * pa_x[i] - ta1_x_yyyy_z_1[i] * pc_x[i];
    }

    // Set up 33-36 components of targeted buffer : HP

    auto ta1_x_xyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 33);

    auto ta1_x_xyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 34);

    auto ta1_x_xyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 35);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xyyy_x_0,  \
                             ta1_x_xyyy_x_1,  \
                             ta1_x_xyyy_y_0,  \
                             ta1_x_xyyy_y_1,  \
                             ta1_x_xyyyz_x_0, \
                             ta1_x_xyyyz_y_0, \
                             ta1_x_xyyyz_z_0, \
                             ta1_x_yyyz_z_0,  \
                             ta1_x_yyyz_z_1,  \
                             ta_yyyz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyyyz_x_0[i] = ta1_x_xyyy_x_0[i] * pa_z[i] - ta1_x_xyyy_x_1[i] * pc_z[i];

        ta1_x_xyyyz_y_0[i] = ta1_x_xyyy_y_0[i] * pa_z[i] - ta1_x_xyyy_y_1[i] * pc_z[i];

        ta1_x_xyyyz_z_0[i] = ta_yyyz_z_1[i] + ta1_x_yyyz_z_0[i] * pa_x[i] - ta1_x_yyyz_z_1[i] * pc_x[i];
    }

    // Set up 36-39 components of targeted buffer : HP

    auto ta1_x_xyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 36);

    auto ta1_x_xyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 37);

    auto ta1_x_xyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 38);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xyyzz_x_0, \
                             ta1_x_xyyzz_y_0, \
                             ta1_x_xyyzz_z_0, \
                             ta1_x_xyzz_x_0,  \
                             ta1_x_xyzz_x_1,  \
                             ta1_x_xzz_x_0,   \
                             ta1_x_xzz_x_1,   \
                             ta1_x_yyzz_y_0,  \
                             ta1_x_yyzz_y_1,  \
                             ta1_x_yyzz_z_0,  \
                             ta1_x_yyzz_z_1,  \
                             ta_yyzz_y_1,     \
                             ta_yyzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzz_x_0[i] = ta1_x_xzz_x_0[i] * fe_0 - ta1_x_xzz_x_1[i] * fe_0 + ta1_x_xyzz_x_0[i] * pa_y[i] - ta1_x_xyzz_x_1[i] * pc_y[i];

        ta1_x_xyyzz_y_0[i] = ta_yyzz_y_1[i] + ta1_x_yyzz_y_0[i] * pa_x[i] - ta1_x_yyzz_y_1[i] * pc_x[i];

        ta1_x_xyyzz_z_0[i] = ta_yyzz_z_1[i] + ta1_x_yyzz_z_0[i] * pa_x[i] - ta1_x_yyzz_z_1[i] * pc_x[i];
    }

    // Set up 39-42 components of targeted buffer : HP

    auto ta1_x_xyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 39);

    auto ta1_x_xyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 40);

    auto ta1_x_xyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 41);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_xyzzz_x_0, \
                             ta1_x_xyzzz_y_0, \
                             ta1_x_xyzzz_z_0, \
                             ta1_x_xzzz_x_0,  \
                             ta1_x_xzzz_x_1,  \
                             ta1_x_xzzz_z_0,  \
                             ta1_x_xzzz_z_1,  \
                             ta1_x_yzzz_y_0,  \
                             ta1_x_yzzz_y_1,  \
                             ta_yzzz_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyzzz_x_0[i] = ta1_x_xzzz_x_0[i] * pa_y[i] - ta1_x_xzzz_x_1[i] * pc_y[i];

        ta1_x_xyzzz_y_0[i] = ta_yzzz_y_1[i] + ta1_x_yzzz_y_0[i] * pa_x[i] - ta1_x_yzzz_y_1[i] * pc_x[i];

        ta1_x_xyzzz_z_0[i] = ta1_x_xzzz_z_0[i] * pa_y[i] - ta1_x_xzzz_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : HP

    auto ta1_x_xzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 42);

    auto ta1_x_xzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 43);

    auto ta1_x_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 44);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_xzz_x_0,   \
                             ta1_x_xzz_x_1,   \
                             ta1_x_xzzz_x_0,  \
                             ta1_x_xzzz_x_1,  \
                             ta1_x_xzzzz_x_0, \
                             ta1_x_xzzzz_y_0, \
                             ta1_x_xzzzz_z_0, \
                             ta1_x_zzzz_y_0,  \
                             ta1_x_zzzz_y_1,  \
                             ta1_x_zzzz_z_0,  \
                             ta1_x_zzzz_z_1,  \
                             ta_zzzz_y_1,     \
                             ta_zzzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzz_x_0[i] =
            3.0 * ta1_x_xzz_x_0[i] * fe_0 - 3.0 * ta1_x_xzz_x_1[i] * fe_0 + ta1_x_xzzz_x_0[i] * pa_z[i] - ta1_x_xzzz_x_1[i] * pc_z[i];

        ta1_x_xzzzz_y_0[i] = ta_zzzz_y_1[i] + ta1_x_zzzz_y_0[i] * pa_x[i] - ta1_x_zzzz_y_1[i] * pc_x[i];

        ta1_x_xzzzz_z_0[i] = ta_zzzz_z_1[i] + ta1_x_zzzz_z_0[i] * pa_x[i] - ta1_x_zzzz_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : HP

    auto ta1_x_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 45);

    auto ta1_x_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 46);

    auto ta1_x_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 47);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_yyy_x_0,   \
                             ta1_x_yyy_x_1,   \
                             ta1_x_yyy_y_0,   \
                             ta1_x_yyy_y_1,   \
                             ta1_x_yyy_z_0,   \
                             ta1_x_yyy_z_1,   \
                             ta1_x_yyyy_0_0,  \
                             ta1_x_yyyy_0_1,  \
                             ta1_x_yyyy_x_0,  \
                             ta1_x_yyyy_x_1,  \
                             ta1_x_yyyy_y_0,  \
                             ta1_x_yyyy_y_1,  \
                             ta1_x_yyyy_z_0,  \
                             ta1_x_yyyy_z_1,  \
                             ta1_x_yyyyy_x_0, \
                             ta1_x_yyyyy_y_0, \
                             ta1_x_yyyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyy_x_0[i] =
            4.0 * ta1_x_yyy_x_0[i] * fe_0 - 4.0 * ta1_x_yyy_x_1[i] * fe_0 + ta1_x_yyyy_x_0[i] * pa_y[i] - ta1_x_yyyy_x_1[i] * pc_y[i];

        ta1_x_yyyyy_y_0[i] = 4.0 * ta1_x_yyy_y_0[i] * fe_0 - 4.0 * ta1_x_yyy_y_1[i] * fe_0 + ta1_x_yyyy_0_0[i] * fe_0 - ta1_x_yyyy_0_1[i] * fe_0 +
                             ta1_x_yyyy_y_0[i] * pa_y[i] - ta1_x_yyyy_y_1[i] * pc_y[i];

        ta1_x_yyyyy_z_0[i] =
            4.0 * ta1_x_yyy_z_0[i] * fe_0 - 4.0 * ta1_x_yyy_z_1[i] * fe_0 + ta1_x_yyyy_z_0[i] * pa_y[i] - ta1_x_yyyy_z_1[i] * pc_y[i];
    }

    // Set up 48-51 components of targeted buffer : HP

    auto ta1_x_yyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 48);

    auto ta1_x_yyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 49);

    auto ta1_x_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 50);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yyyy_x_0,  \
                             ta1_x_yyyy_x_1,  \
                             ta1_x_yyyy_y_0,  \
                             ta1_x_yyyy_y_1,  \
                             ta1_x_yyyyz_x_0, \
                             ta1_x_yyyyz_y_0, \
                             ta1_x_yyyyz_z_0, \
                             ta1_x_yyyz_z_0,  \
                             ta1_x_yyyz_z_1,  \
                             ta1_x_yyz_z_0,   \
                             ta1_x_yyz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyz_x_0[i] = ta1_x_yyyy_x_0[i] * pa_z[i] - ta1_x_yyyy_x_1[i] * pc_z[i];

        ta1_x_yyyyz_y_0[i] = ta1_x_yyyy_y_0[i] * pa_z[i] - ta1_x_yyyy_y_1[i] * pc_z[i];

        ta1_x_yyyyz_z_0[i] =
            3.0 * ta1_x_yyz_z_0[i] * fe_0 - 3.0 * ta1_x_yyz_z_1[i] * fe_0 + ta1_x_yyyz_z_0[i] * pa_y[i] - ta1_x_yyyz_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : HP

    auto ta1_x_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 51);

    auto ta1_x_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 52);

    auto ta1_x_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 53);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yyy_y_0,   \
                             ta1_x_yyy_y_1,   \
                             ta1_x_yyyz_y_0,  \
                             ta1_x_yyyz_y_1,  \
                             ta1_x_yyyzz_x_0, \
                             ta1_x_yyyzz_y_0, \
                             ta1_x_yyyzz_z_0, \
                             ta1_x_yyzz_x_0,  \
                             ta1_x_yyzz_x_1,  \
                             ta1_x_yyzz_z_0,  \
                             ta1_x_yyzz_z_1,  \
                             ta1_x_yzz_x_0,   \
                             ta1_x_yzz_x_1,   \
                             ta1_x_yzz_z_0,   \
                             ta1_x_yzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzz_x_0[i] =
            2.0 * ta1_x_yzz_x_0[i] * fe_0 - 2.0 * ta1_x_yzz_x_1[i] * fe_0 + ta1_x_yyzz_x_0[i] * pa_y[i] - ta1_x_yyzz_x_1[i] * pc_y[i];

        ta1_x_yyyzz_y_0[i] = ta1_x_yyy_y_0[i] * fe_0 - ta1_x_yyy_y_1[i] * fe_0 + ta1_x_yyyz_y_0[i] * pa_z[i] - ta1_x_yyyz_y_1[i] * pc_z[i];

        ta1_x_yyyzz_z_0[i] =
            2.0 * ta1_x_yzz_z_0[i] * fe_0 - 2.0 * ta1_x_yzz_z_1[i] * fe_0 + ta1_x_yyzz_z_0[i] * pa_y[i] - ta1_x_yyzz_z_1[i] * pc_y[i];
    }

    // Set up 54-57 components of targeted buffer : HP

    auto ta1_x_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 54);

    auto ta1_x_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 55);

    auto ta1_x_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 56);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_yyz_y_0,   \
                             ta1_x_yyz_y_1,   \
                             ta1_x_yyzz_y_0,  \
                             ta1_x_yyzz_y_1,  \
                             ta1_x_yyzzz_x_0, \
                             ta1_x_yyzzz_y_0, \
                             ta1_x_yyzzz_z_0, \
                             ta1_x_yzzz_x_0,  \
                             ta1_x_yzzz_x_1,  \
                             ta1_x_yzzz_z_0,  \
                             ta1_x_yzzz_z_1,  \
                             ta1_x_zzz_x_0,   \
                             ta1_x_zzz_x_1,   \
                             ta1_x_zzz_z_0,   \
                             ta1_x_zzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzz_x_0[i] = ta1_x_zzz_x_0[i] * fe_0 - ta1_x_zzz_x_1[i] * fe_0 + ta1_x_yzzz_x_0[i] * pa_y[i] - ta1_x_yzzz_x_1[i] * pc_y[i];

        ta1_x_yyzzz_y_0[i] =
            2.0 * ta1_x_yyz_y_0[i] * fe_0 - 2.0 * ta1_x_yyz_y_1[i] * fe_0 + ta1_x_yyzz_y_0[i] * pa_z[i] - ta1_x_yyzz_y_1[i] * pc_z[i];

        ta1_x_yyzzz_z_0[i] = ta1_x_zzz_z_0[i] * fe_0 - ta1_x_zzz_z_1[i] * fe_0 + ta1_x_yzzz_z_0[i] * pa_y[i] - ta1_x_yzzz_z_1[i] * pc_y[i];
    }

    // Set up 57-60 components of targeted buffer : HP

    auto ta1_x_yzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 57);

    auto ta1_x_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 58);

    auto ta1_x_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 59);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_yzzzz_x_0, \
                             ta1_x_yzzzz_y_0, \
                             ta1_x_yzzzz_z_0, \
                             ta1_x_zzzz_0_0,  \
                             ta1_x_zzzz_0_1,  \
                             ta1_x_zzzz_x_0,  \
                             ta1_x_zzzz_x_1,  \
                             ta1_x_zzzz_y_0,  \
                             ta1_x_zzzz_y_1,  \
                             ta1_x_zzzz_z_0,  \
                             ta1_x_zzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzz_x_0[i] = ta1_x_zzzz_x_0[i] * pa_y[i] - ta1_x_zzzz_x_1[i] * pc_y[i];

        ta1_x_yzzzz_y_0[i] = ta1_x_zzzz_0_0[i] * fe_0 - ta1_x_zzzz_0_1[i] * fe_0 + ta1_x_zzzz_y_0[i] * pa_y[i] - ta1_x_zzzz_y_1[i] * pc_y[i];

        ta1_x_yzzzz_z_0[i] = ta1_x_zzzz_z_0[i] * pa_y[i] - ta1_x_zzzz_z_1[i] * pc_y[i];
    }

    // Set up 60-63 components of targeted buffer : HP

    auto ta1_x_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 60);

    auto ta1_x_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 61);

    auto ta1_x_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 62);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_zzz_x_0,   \
                             ta1_x_zzz_x_1,   \
                             ta1_x_zzz_y_0,   \
                             ta1_x_zzz_y_1,   \
                             ta1_x_zzz_z_0,   \
                             ta1_x_zzz_z_1,   \
                             ta1_x_zzzz_0_0,  \
                             ta1_x_zzzz_0_1,  \
                             ta1_x_zzzz_x_0,  \
                             ta1_x_zzzz_x_1,  \
                             ta1_x_zzzz_y_0,  \
                             ta1_x_zzzz_y_1,  \
                             ta1_x_zzzz_z_0,  \
                             ta1_x_zzzz_z_1,  \
                             ta1_x_zzzzz_x_0, \
                             ta1_x_zzzzz_y_0, \
                             ta1_x_zzzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzz_x_0[i] =
            4.0 * ta1_x_zzz_x_0[i] * fe_0 - 4.0 * ta1_x_zzz_x_1[i] * fe_0 + ta1_x_zzzz_x_0[i] * pa_z[i] - ta1_x_zzzz_x_1[i] * pc_z[i];

        ta1_x_zzzzz_y_0[i] =
            4.0 * ta1_x_zzz_y_0[i] * fe_0 - 4.0 * ta1_x_zzz_y_1[i] * fe_0 + ta1_x_zzzz_y_0[i] * pa_z[i] - ta1_x_zzzz_y_1[i] * pc_z[i];

        ta1_x_zzzzz_z_0[i] = 4.0 * ta1_x_zzz_z_0[i] * fe_0 - 4.0 * ta1_x_zzz_z_1[i] * fe_0 + ta1_x_zzzz_0_0[i] * fe_0 - ta1_x_zzzz_0_1[i] * fe_0 +
                             ta1_x_zzzz_z_0[i] * pa_z[i] - ta1_x_zzzz_z_1[i] * pc_z[i];
    }

    // Set up 63-66 components of targeted buffer : HP

    auto ta1_y_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 63);

    auto ta1_y_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 64);

    auto ta1_y_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 65);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xxx_x_0,   \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxx_y_0,   \
                             ta1_y_xxx_y_1,   \
                             ta1_y_xxx_z_0,   \
                             ta1_y_xxx_z_1,   \
                             ta1_y_xxxx_0_0,  \
                             ta1_y_xxxx_0_1,  \
                             ta1_y_xxxx_x_0,  \
                             ta1_y_xxxx_x_1,  \
                             ta1_y_xxxx_y_0,  \
                             ta1_y_xxxx_y_1,  \
                             ta1_y_xxxx_z_0,  \
                             ta1_y_xxxx_z_1,  \
                             ta1_y_xxxxx_x_0, \
                             ta1_y_xxxxx_y_0, \
                             ta1_y_xxxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxx_x_0[i] = 4.0 * ta1_y_xxx_x_0[i] * fe_0 - 4.0 * ta1_y_xxx_x_1[i] * fe_0 + ta1_y_xxxx_0_0[i] * fe_0 - ta1_y_xxxx_0_1[i] * fe_0 +
                             ta1_y_xxxx_x_0[i] * pa_x[i] - ta1_y_xxxx_x_1[i] * pc_x[i];

        ta1_y_xxxxx_y_0[i] =
            4.0 * ta1_y_xxx_y_0[i] * fe_0 - 4.0 * ta1_y_xxx_y_1[i] * fe_0 + ta1_y_xxxx_y_0[i] * pa_x[i] - ta1_y_xxxx_y_1[i] * pc_x[i];

        ta1_y_xxxxx_z_0[i] =
            4.0 * ta1_y_xxx_z_0[i] * fe_0 - 4.0 * ta1_y_xxx_z_1[i] * fe_0 + ta1_y_xxxx_z_0[i] * pa_x[i] - ta1_y_xxxx_z_1[i] * pc_x[i];
    }

    // Set up 66-69 components of targeted buffer : HP

    auto ta1_y_xxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 66);

    auto ta1_y_xxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 67);

    auto ta1_y_xxxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 68);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xxxx_x_0,  \
                             ta1_y_xxxx_x_1,  \
                             ta1_y_xxxx_z_0,  \
                             ta1_y_xxxx_z_1,  \
                             ta1_y_xxxxy_x_0, \
                             ta1_y_xxxxy_y_0, \
                             ta1_y_xxxxy_z_0, \
                             ta1_y_xxxy_y_0,  \
                             ta1_y_xxxy_y_1,  \
                             ta1_y_xxy_y_0,   \
                             ta1_y_xxy_y_1,   \
                             ta_xxxx_x_1,     \
                             ta_xxxx_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxy_x_0[i] = ta_xxxx_x_1[i] + ta1_y_xxxx_x_0[i] * pa_y[i] - ta1_y_xxxx_x_1[i] * pc_y[i];

        ta1_y_xxxxy_y_0[i] =
            3.0 * ta1_y_xxy_y_0[i] * fe_0 - 3.0 * ta1_y_xxy_y_1[i] * fe_0 + ta1_y_xxxy_y_0[i] * pa_x[i] - ta1_y_xxxy_y_1[i] * pc_x[i];

        ta1_y_xxxxy_z_0[i] = ta_xxxx_z_1[i] + ta1_y_xxxx_z_0[i] * pa_y[i] - ta1_y_xxxx_z_1[i] * pc_y[i];
    }

    // Set up 69-72 components of targeted buffer : HP

    auto ta1_y_xxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 69);

    auto ta1_y_xxxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 70);

    auto ta1_y_xxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 71);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xxxx_x_0,  \
                             ta1_y_xxxx_x_1,  \
                             ta1_y_xxxx_y_0,  \
                             ta1_y_xxxx_y_1,  \
                             ta1_y_xxxxz_x_0, \
                             ta1_y_xxxxz_y_0, \
                             ta1_y_xxxxz_z_0, \
                             ta1_y_xxxz_z_0,  \
                             ta1_y_xxxz_z_1,  \
                             ta1_y_xxz_z_0,   \
                             ta1_y_xxz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxz_x_0[i] = ta1_y_xxxx_x_0[i] * pa_z[i] - ta1_y_xxxx_x_1[i] * pc_z[i];

        ta1_y_xxxxz_y_0[i] = ta1_y_xxxx_y_0[i] * pa_z[i] - ta1_y_xxxx_y_1[i] * pc_z[i];

        ta1_y_xxxxz_z_0[i] =
            3.0 * ta1_y_xxz_z_0[i] * fe_0 - 3.0 * ta1_y_xxz_z_1[i] * fe_0 + ta1_y_xxxz_z_0[i] * pa_x[i] - ta1_y_xxxz_z_1[i] * pc_x[i];
    }

    // Set up 72-75 components of targeted buffer : HP

    auto ta1_y_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 72);

    auto ta1_y_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 73);

    auto ta1_y_xxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 74);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xxx_x_0,   \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxxy_x_0,  \
                             ta1_y_xxxy_x_1,  \
                             ta1_y_xxxyy_x_0, \
                             ta1_y_xxxyy_y_0, \
                             ta1_y_xxxyy_z_0, \
                             ta1_y_xxyy_y_0,  \
                             ta1_y_xxyy_y_1,  \
                             ta1_y_xxyy_z_0,  \
                             ta1_y_xxyy_z_1,  \
                             ta1_y_xyy_y_0,   \
                             ta1_y_xyy_y_1,   \
                             ta1_y_xyy_z_0,   \
                             ta1_y_xyy_z_1,   \
                             ta_xxxy_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyy_x_0[i] =
            ta1_y_xxx_x_0[i] * fe_0 - ta1_y_xxx_x_1[i] * fe_0 + ta_xxxy_x_1[i] + ta1_y_xxxy_x_0[i] * pa_y[i] - ta1_y_xxxy_x_1[i] * pc_y[i];

        ta1_y_xxxyy_y_0[i] =
            2.0 * ta1_y_xyy_y_0[i] * fe_0 - 2.0 * ta1_y_xyy_y_1[i] * fe_0 + ta1_y_xxyy_y_0[i] * pa_x[i] - ta1_y_xxyy_y_1[i] * pc_x[i];

        ta1_y_xxxyy_z_0[i] =
            2.0 * ta1_y_xyy_z_0[i] * fe_0 - 2.0 * ta1_y_xyy_z_1[i] * fe_0 + ta1_y_xxyy_z_0[i] * pa_x[i] - ta1_y_xxyy_z_1[i] * pc_x[i];
    }

    // Set up 75-78 components of targeted buffer : HP

    auto ta1_y_xxxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 75);

    auto ta1_y_xxxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 76);

    auto ta1_y_xxxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 77);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_xxxy_x_0,  \
                             ta1_y_xxxy_x_1,  \
                             ta1_y_xxxy_y_0,  \
                             ta1_y_xxxy_y_1,  \
                             ta1_y_xxxyz_x_0, \
                             ta1_y_xxxyz_y_0, \
                             ta1_y_xxxyz_z_0, \
                             ta1_y_xxxz_z_0,  \
                             ta1_y_xxxz_z_1,  \
                             ta_xxxz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xxxyz_x_0[i] = ta1_y_xxxy_x_0[i] * pa_z[i] - ta1_y_xxxy_x_1[i] * pc_z[i];

        ta1_y_xxxyz_y_0[i] = ta1_y_xxxy_y_0[i] * pa_z[i] - ta1_y_xxxy_y_1[i] * pc_z[i];

        ta1_y_xxxyz_z_0[i] = ta_xxxz_z_1[i] + ta1_y_xxxz_z_0[i] * pa_y[i] - ta1_y_xxxz_z_1[i] * pc_y[i];
    }

    // Set up 78-81 components of targeted buffer : HP

    auto ta1_y_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 78);

    auto ta1_y_xxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 79);

    auto ta1_y_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 80);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xxx_x_0,   \
                             ta1_y_xxx_x_1,   \
                             ta1_y_xxxz_x_0,  \
                             ta1_y_xxxz_x_1,  \
                             ta1_y_xxxzz_x_0, \
                             ta1_y_xxxzz_y_0, \
                             ta1_y_xxxzz_z_0, \
                             ta1_y_xxzz_y_0,  \
                             ta1_y_xxzz_y_1,  \
                             ta1_y_xxzz_z_0,  \
                             ta1_y_xxzz_z_1,  \
                             ta1_y_xzz_y_0,   \
                             ta1_y_xzz_y_1,   \
                             ta1_y_xzz_z_0,   \
                             ta1_y_xzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzz_x_0[i] = ta1_y_xxx_x_0[i] * fe_0 - ta1_y_xxx_x_1[i] * fe_0 + ta1_y_xxxz_x_0[i] * pa_z[i] - ta1_y_xxxz_x_1[i] * pc_z[i];

        ta1_y_xxxzz_y_0[i] =
            2.0 * ta1_y_xzz_y_0[i] * fe_0 - 2.0 * ta1_y_xzz_y_1[i] * fe_0 + ta1_y_xxzz_y_0[i] * pa_x[i] - ta1_y_xxzz_y_1[i] * pc_x[i];

        ta1_y_xxxzz_z_0[i] =
            2.0 * ta1_y_xzz_z_0[i] * fe_0 - 2.0 * ta1_y_xzz_z_1[i] * fe_0 + ta1_y_xxzz_z_0[i] * pa_x[i] - ta1_y_xxzz_z_1[i] * pc_x[i];
    }

    // Set up 81-84 components of targeted buffer : HP

    auto ta1_y_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 81);

    auto ta1_y_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 82);

    auto ta1_y_xxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 83);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xxy_x_0,   \
                             ta1_y_xxy_x_1,   \
                             ta1_y_xxyy_x_0,  \
                             ta1_y_xxyy_x_1,  \
                             ta1_y_xxyyy_x_0, \
                             ta1_y_xxyyy_y_0, \
                             ta1_y_xxyyy_z_0, \
                             ta1_y_xyyy_y_0,  \
                             ta1_y_xyyy_y_1,  \
                             ta1_y_xyyy_z_0,  \
                             ta1_y_xyyy_z_1,  \
                             ta1_y_yyy_y_0,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyy_z_0,   \
                             ta1_y_yyy_z_1,   \
                             ta_xxyy_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyy_x_0[i] = 2.0 * ta1_y_xxy_x_0[i] * fe_0 - 2.0 * ta1_y_xxy_x_1[i] * fe_0 + ta_xxyy_x_1[i] + ta1_y_xxyy_x_0[i] * pa_y[i] -
                             ta1_y_xxyy_x_1[i] * pc_y[i];

        ta1_y_xxyyy_y_0[i] = ta1_y_yyy_y_0[i] * fe_0 - ta1_y_yyy_y_1[i] * fe_0 + ta1_y_xyyy_y_0[i] * pa_x[i] - ta1_y_xyyy_y_1[i] * pc_x[i];

        ta1_y_xxyyy_z_0[i] = ta1_y_yyy_z_0[i] * fe_0 - ta1_y_yyy_z_1[i] * fe_0 + ta1_y_xyyy_z_0[i] * pa_x[i] - ta1_y_xyyy_z_1[i] * pc_x[i];
    }

    // Set up 84-87 components of targeted buffer : HP

    auto ta1_y_xxyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 84);

    auto ta1_y_xxyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 85);

    auto ta1_y_xxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 86);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xxyy_x_0,  \
                             ta1_y_xxyy_x_1,  \
                             ta1_y_xxyy_y_0,  \
                             ta1_y_xxyy_y_1,  \
                             ta1_y_xxyyz_x_0, \
                             ta1_y_xxyyz_y_0, \
                             ta1_y_xxyyz_z_0, \
                             ta1_y_xyyz_z_0,  \
                             ta1_y_xyyz_z_1,  \
                             ta1_y_yyz_z_0,   \
                             ta1_y_yyz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyz_x_0[i] = ta1_y_xxyy_x_0[i] * pa_z[i] - ta1_y_xxyy_x_1[i] * pc_z[i];

        ta1_y_xxyyz_y_0[i] = ta1_y_xxyy_y_0[i] * pa_z[i] - ta1_y_xxyy_y_1[i] * pc_z[i];

        ta1_y_xxyyz_z_0[i] = ta1_y_yyz_z_0[i] * fe_0 - ta1_y_yyz_z_1[i] * fe_0 + ta1_y_xyyz_z_0[i] * pa_x[i] - ta1_y_xyyz_z_1[i] * pc_x[i];
    }

    // Set up 87-90 components of targeted buffer : HP

    auto ta1_y_xxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 87);

    auto ta1_y_xxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 88);

    auto ta1_y_xxyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 89);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xxyzz_x_0, \
                             ta1_y_xxyzz_y_0, \
                             ta1_y_xxyzz_z_0, \
                             ta1_y_xxzz_x_0,  \
                             ta1_y_xxzz_x_1,  \
                             ta1_y_xxzz_z_0,  \
                             ta1_y_xxzz_z_1,  \
                             ta1_y_xyzz_y_0,  \
                             ta1_y_xyzz_y_1,  \
                             ta1_y_yzz_y_0,   \
                             ta1_y_yzz_y_1,   \
                             ta_xxzz_x_1,     \
                             ta_xxzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzz_x_0[i] = ta_xxzz_x_1[i] + ta1_y_xxzz_x_0[i] * pa_y[i] - ta1_y_xxzz_x_1[i] * pc_y[i];

        ta1_y_xxyzz_y_0[i] = ta1_y_yzz_y_0[i] * fe_0 - ta1_y_yzz_y_1[i] * fe_0 + ta1_y_xyzz_y_0[i] * pa_x[i] - ta1_y_xyzz_y_1[i] * pc_x[i];

        ta1_y_xxyzz_z_0[i] = ta_xxzz_z_1[i] + ta1_y_xxzz_z_0[i] * pa_y[i] - ta1_y_xxzz_z_1[i] * pc_y[i];
    }

    // Set up 90-93 components of targeted buffer : HP

    auto ta1_y_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 90);

    auto ta1_y_xxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 91);

    auto ta1_y_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 92);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xxz_x_0,   \
                             ta1_y_xxz_x_1,   \
                             ta1_y_xxzz_x_0,  \
                             ta1_y_xxzz_x_1,  \
                             ta1_y_xxzzz_x_0, \
                             ta1_y_xxzzz_y_0, \
                             ta1_y_xxzzz_z_0, \
                             ta1_y_xzzz_y_0,  \
                             ta1_y_xzzz_y_1,  \
                             ta1_y_xzzz_z_0,  \
                             ta1_y_xzzz_z_1,  \
                             ta1_y_zzz_y_0,   \
                             ta1_y_zzz_y_1,   \
                             ta1_y_zzz_z_0,   \
                             ta1_y_zzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzz_x_0[i] =
            2.0 * ta1_y_xxz_x_0[i] * fe_0 - 2.0 * ta1_y_xxz_x_1[i] * fe_0 + ta1_y_xxzz_x_0[i] * pa_z[i] - ta1_y_xxzz_x_1[i] * pc_z[i];

        ta1_y_xxzzz_y_0[i] = ta1_y_zzz_y_0[i] * fe_0 - ta1_y_zzz_y_1[i] * fe_0 + ta1_y_xzzz_y_0[i] * pa_x[i] - ta1_y_xzzz_y_1[i] * pc_x[i];

        ta1_y_xxzzz_z_0[i] = ta1_y_zzz_z_0[i] * fe_0 - ta1_y_zzz_z_1[i] * fe_0 + ta1_y_xzzz_z_0[i] * pa_x[i] - ta1_y_xzzz_z_1[i] * pc_x[i];
    }

    // Set up 93-96 components of targeted buffer : HP

    auto ta1_y_xyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 93);

    auto ta1_y_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 94);

    auto ta1_y_xyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 95);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xyyyy_x_0, \
                             ta1_y_xyyyy_y_0, \
                             ta1_y_xyyyy_z_0, \
                             ta1_y_yyyy_0_0,  \
                             ta1_y_yyyy_0_1,  \
                             ta1_y_yyyy_x_0,  \
                             ta1_y_yyyy_x_1,  \
                             ta1_y_yyyy_y_0,  \
                             ta1_y_yyyy_y_1,  \
                             ta1_y_yyyy_z_0,  \
                             ta1_y_yyyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyy_x_0[i] = ta1_y_yyyy_0_0[i] * fe_0 - ta1_y_yyyy_0_1[i] * fe_0 + ta1_y_yyyy_x_0[i] * pa_x[i] - ta1_y_yyyy_x_1[i] * pc_x[i];

        ta1_y_xyyyy_y_0[i] = ta1_y_yyyy_y_0[i] * pa_x[i] - ta1_y_yyyy_y_1[i] * pc_x[i];

        ta1_y_xyyyy_z_0[i] = ta1_y_yyyy_z_0[i] * pa_x[i] - ta1_y_yyyy_z_1[i] * pc_x[i];
    }

    // Set up 96-99 components of targeted buffer : HP

    auto ta1_y_xyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 96);

    auto ta1_y_xyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 97);

    auto ta1_y_xyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 98);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_xyyy_x_0,  \
                             ta1_y_xyyy_x_1,  \
                             ta1_y_xyyyz_x_0, \
                             ta1_y_xyyyz_y_0, \
                             ta1_y_xyyyz_z_0, \
                             ta1_y_yyyz_y_0,  \
                             ta1_y_yyyz_y_1,  \
                             ta1_y_yyyz_z_0,  \
                             ta1_y_yyyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyyyz_x_0[i] = ta1_y_xyyy_x_0[i] * pa_z[i] - ta1_y_xyyy_x_1[i] * pc_z[i];

        ta1_y_xyyyz_y_0[i] = ta1_y_yyyz_y_0[i] * pa_x[i] - ta1_y_yyyz_y_1[i] * pc_x[i];

        ta1_y_xyyyz_z_0[i] = ta1_y_yyyz_z_0[i] * pa_x[i] - ta1_y_yyyz_z_1[i] * pc_x[i];
    }

    // Set up 99-102 components of targeted buffer : HP

    auto ta1_y_xyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 99);

    auto ta1_y_xyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 100);

    auto ta1_y_xyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 101);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xyyzz_x_0, \
                             ta1_y_xyyzz_y_0, \
                             ta1_y_xyyzz_z_0, \
                             ta1_y_yyzz_0_0,  \
                             ta1_y_yyzz_0_1,  \
                             ta1_y_yyzz_x_0,  \
                             ta1_y_yyzz_x_1,  \
                             ta1_y_yyzz_y_0,  \
                             ta1_y_yyzz_y_1,  \
                             ta1_y_yyzz_z_0,  \
                             ta1_y_yyzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzz_x_0[i] = ta1_y_yyzz_0_0[i] * fe_0 - ta1_y_yyzz_0_1[i] * fe_0 + ta1_y_yyzz_x_0[i] * pa_x[i] - ta1_y_yyzz_x_1[i] * pc_x[i];

        ta1_y_xyyzz_y_0[i] = ta1_y_yyzz_y_0[i] * pa_x[i] - ta1_y_yyzz_y_1[i] * pc_x[i];

        ta1_y_xyyzz_z_0[i] = ta1_y_yyzz_z_0[i] * pa_x[i] - ta1_y_yyzz_z_1[i] * pc_x[i];
    }

    // Set up 102-105 components of targeted buffer : HP

    auto ta1_y_xyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 102);

    auto ta1_y_xyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 103);

    auto ta1_y_xyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 104);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_xyzzz_x_0, \
                             ta1_y_xyzzz_y_0, \
                             ta1_y_xyzzz_z_0, \
                             ta1_y_xzzz_x_0,  \
                             ta1_y_xzzz_x_1,  \
                             ta1_y_yzzz_y_0,  \
                             ta1_y_yzzz_y_1,  \
                             ta1_y_yzzz_z_0,  \
                             ta1_y_yzzz_z_1,  \
                             ta_xzzz_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyzzz_x_0[i] = ta_xzzz_x_1[i] + ta1_y_xzzz_x_0[i] * pa_y[i] - ta1_y_xzzz_x_1[i] * pc_y[i];

        ta1_y_xyzzz_y_0[i] = ta1_y_yzzz_y_0[i] * pa_x[i] - ta1_y_yzzz_y_1[i] * pc_x[i];

        ta1_y_xyzzz_z_0[i] = ta1_y_yzzz_z_0[i] * pa_x[i] - ta1_y_yzzz_z_1[i] * pc_x[i];
    }

    // Set up 105-108 components of targeted buffer : HP

    auto ta1_y_xzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 105);

    auto ta1_y_xzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 106);

    auto ta1_y_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 107);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_xzzzz_x_0, \
                             ta1_y_xzzzz_y_0, \
                             ta1_y_xzzzz_z_0, \
                             ta1_y_zzzz_0_0,  \
                             ta1_y_zzzz_0_1,  \
                             ta1_y_zzzz_x_0,  \
                             ta1_y_zzzz_x_1,  \
                             ta1_y_zzzz_y_0,  \
                             ta1_y_zzzz_y_1,  \
                             ta1_y_zzzz_z_0,  \
                             ta1_y_zzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzz_x_0[i] = ta1_y_zzzz_0_0[i] * fe_0 - ta1_y_zzzz_0_1[i] * fe_0 + ta1_y_zzzz_x_0[i] * pa_x[i] - ta1_y_zzzz_x_1[i] * pc_x[i];

        ta1_y_xzzzz_y_0[i] = ta1_y_zzzz_y_0[i] * pa_x[i] - ta1_y_zzzz_y_1[i] * pc_x[i];

        ta1_y_xzzzz_z_0[i] = ta1_y_zzzz_z_0[i] * pa_x[i] - ta1_y_zzzz_z_1[i] * pc_x[i];
    }

    // Set up 108-111 components of targeted buffer : HP

    auto ta1_y_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 108);

    auto ta1_y_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 109);

    auto ta1_y_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 110);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_y_yyy_x_0,   \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_y_0,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyy_z_0,   \
                             ta1_y_yyy_z_1,   \
                             ta1_y_yyyy_0_0,  \
                             ta1_y_yyyy_0_1,  \
                             ta1_y_yyyy_x_0,  \
                             ta1_y_yyyy_x_1,  \
                             ta1_y_yyyy_y_0,  \
                             ta1_y_yyyy_y_1,  \
                             ta1_y_yyyy_z_0,  \
                             ta1_y_yyyy_z_1,  \
                             ta1_y_yyyyy_x_0, \
                             ta1_y_yyyyy_y_0, \
                             ta1_y_yyyyy_z_0, \
                             ta_yyyy_x_1,     \
                             ta_yyyy_y_1,     \
                             ta_yyyy_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyy_x_0[i] = 4.0 * ta1_y_yyy_x_0[i] * fe_0 - 4.0 * ta1_y_yyy_x_1[i] * fe_0 + ta_yyyy_x_1[i] + ta1_y_yyyy_x_0[i] * pa_y[i] -
                             ta1_y_yyyy_x_1[i] * pc_y[i];

        ta1_y_yyyyy_y_0[i] = 4.0 * ta1_y_yyy_y_0[i] * fe_0 - 4.0 * ta1_y_yyy_y_1[i] * fe_0 + ta1_y_yyyy_0_0[i] * fe_0 - ta1_y_yyyy_0_1[i] * fe_0 +
                             ta_yyyy_y_1[i] + ta1_y_yyyy_y_0[i] * pa_y[i] - ta1_y_yyyy_y_1[i] * pc_y[i];

        ta1_y_yyyyy_z_0[i] = 4.0 * ta1_y_yyy_z_0[i] * fe_0 - 4.0 * ta1_y_yyy_z_1[i] * fe_0 + ta_yyyy_z_1[i] + ta1_y_yyyy_z_0[i] * pa_y[i] -
                             ta1_y_yyyy_z_1[i] * pc_y[i];
    }

    // Set up 111-114 components of targeted buffer : HP

    auto ta1_y_yyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 111);

    auto ta1_y_yyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 112);

    auto ta1_y_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 113);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_yyyy_0_0,  \
                             ta1_y_yyyy_0_1,  \
                             ta1_y_yyyy_x_0,  \
                             ta1_y_yyyy_x_1,  \
                             ta1_y_yyyy_y_0,  \
                             ta1_y_yyyy_y_1,  \
                             ta1_y_yyyy_z_0,  \
                             ta1_y_yyyy_z_1,  \
                             ta1_y_yyyyz_x_0, \
                             ta1_y_yyyyz_y_0, \
                             ta1_y_yyyyz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyz_x_0[i] = ta1_y_yyyy_x_0[i] * pa_z[i] - ta1_y_yyyy_x_1[i] * pc_z[i];

        ta1_y_yyyyz_y_0[i] = ta1_y_yyyy_y_0[i] * pa_z[i] - ta1_y_yyyy_y_1[i] * pc_z[i];

        ta1_y_yyyyz_z_0[i] = ta1_y_yyyy_0_0[i] * fe_0 - ta1_y_yyyy_0_1[i] * fe_0 + ta1_y_yyyy_z_0[i] * pa_z[i] - ta1_y_yyyy_z_1[i] * pc_z[i];
    }

    // Set up 114-117 components of targeted buffer : HP

    auto ta1_y_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 114);

    auto ta1_y_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 115);

    auto ta1_y_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 116);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yyy_x_0,   \
                             ta1_y_yyy_x_1,   \
                             ta1_y_yyy_y_0,   \
                             ta1_y_yyy_y_1,   \
                             ta1_y_yyyz_x_0,  \
                             ta1_y_yyyz_x_1,  \
                             ta1_y_yyyz_y_0,  \
                             ta1_y_yyyz_y_1,  \
                             ta1_y_yyyzz_x_0, \
                             ta1_y_yyyzz_y_0, \
                             ta1_y_yyyzz_z_0, \
                             ta1_y_yyzz_z_0,  \
                             ta1_y_yyzz_z_1,  \
                             ta1_y_yzz_z_0,   \
                             ta1_y_yzz_z_1,   \
                             ta_yyzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzz_x_0[i] = ta1_y_yyy_x_0[i] * fe_0 - ta1_y_yyy_x_1[i] * fe_0 + ta1_y_yyyz_x_0[i] * pa_z[i] - ta1_y_yyyz_x_1[i] * pc_z[i];

        ta1_y_yyyzz_y_0[i] = ta1_y_yyy_y_0[i] * fe_0 - ta1_y_yyy_y_1[i] * fe_0 + ta1_y_yyyz_y_0[i] * pa_z[i] - ta1_y_yyyz_y_1[i] * pc_z[i];

        ta1_y_yyyzz_z_0[i] = 2.0 * ta1_y_yzz_z_0[i] * fe_0 - 2.0 * ta1_y_yzz_z_1[i] * fe_0 + ta_yyzz_z_1[i] + ta1_y_yyzz_z_0[i] * pa_y[i] -
                             ta1_y_yyzz_z_1[i] * pc_y[i];
    }

    // Set up 117-120 components of targeted buffer : HP

    auto ta1_y_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 117);

    auto ta1_y_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 118);

    auto ta1_y_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 119);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yyz_x_0,   \
                             ta1_y_yyz_x_1,   \
                             ta1_y_yyz_y_0,   \
                             ta1_y_yyz_y_1,   \
                             ta1_y_yyzz_x_0,  \
                             ta1_y_yyzz_x_1,  \
                             ta1_y_yyzz_y_0,  \
                             ta1_y_yyzz_y_1,  \
                             ta1_y_yyzzz_x_0, \
                             ta1_y_yyzzz_y_0, \
                             ta1_y_yyzzz_z_0, \
                             ta1_y_yzzz_z_0,  \
                             ta1_y_yzzz_z_1,  \
                             ta1_y_zzz_z_0,   \
                             ta1_y_zzz_z_1,   \
                             ta_yzzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzz_x_0[i] =
            2.0 * ta1_y_yyz_x_0[i] * fe_0 - 2.0 * ta1_y_yyz_x_1[i] * fe_0 + ta1_y_yyzz_x_0[i] * pa_z[i] - ta1_y_yyzz_x_1[i] * pc_z[i];

        ta1_y_yyzzz_y_0[i] =
            2.0 * ta1_y_yyz_y_0[i] * fe_0 - 2.0 * ta1_y_yyz_y_1[i] * fe_0 + ta1_y_yyzz_y_0[i] * pa_z[i] - ta1_y_yyzz_y_1[i] * pc_z[i];

        ta1_y_yyzzz_z_0[i] =
            ta1_y_zzz_z_0[i] * fe_0 - ta1_y_zzz_z_1[i] * fe_0 + ta_yzzz_z_1[i] + ta1_y_yzzz_z_0[i] * pa_y[i] - ta1_y_yzzz_z_1[i] * pc_y[i];
    }

    // Set up 120-123 components of targeted buffer : HP

    auto ta1_y_yzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 120);

    auto ta1_y_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 121);

    auto ta1_y_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 122);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_yzz_y_0,   \
                             ta1_y_yzz_y_1,   \
                             ta1_y_yzzz_y_0,  \
                             ta1_y_yzzz_y_1,  \
                             ta1_y_yzzzz_x_0, \
                             ta1_y_yzzzz_y_0, \
                             ta1_y_yzzzz_z_0, \
                             ta1_y_zzzz_x_0,  \
                             ta1_y_zzzz_x_1,  \
                             ta1_y_zzzz_z_0,  \
                             ta1_y_zzzz_z_1,  \
                             ta_zzzz_x_1,     \
                             ta_zzzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzz_x_0[i] = ta_zzzz_x_1[i] + ta1_y_zzzz_x_0[i] * pa_y[i] - ta1_y_zzzz_x_1[i] * pc_y[i];

        ta1_y_yzzzz_y_0[i] =
            3.0 * ta1_y_yzz_y_0[i] * fe_0 - 3.0 * ta1_y_yzz_y_1[i] * fe_0 + ta1_y_yzzz_y_0[i] * pa_z[i] - ta1_y_yzzz_y_1[i] * pc_z[i];

        ta1_y_yzzzz_z_0[i] = ta_zzzz_z_1[i] + ta1_y_zzzz_z_0[i] * pa_y[i] - ta1_y_zzzz_z_1[i] * pc_y[i];
    }

    // Set up 123-126 components of targeted buffer : HP

    auto ta1_y_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 123);

    auto ta1_y_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 124);

    auto ta1_y_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 125);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_zzz_x_0,   \
                             ta1_y_zzz_x_1,   \
                             ta1_y_zzz_y_0,   \
                             ta1_y_zzz_y_1,   \
                             ta1_y_zzz_z_0,   \
                             ta1_y_zzz_z_1,   \
                             ta1_y_zzzz_0_0,  \
                             ta1_y_zzzz_0_1,  \
                             ta1_y_zzzz_x_0,  \
                             ta1_y_zzzz_x_1,  \
                             ta1_y_zzzz_y_0,  \
                             ta1_y_zzzz_y_1,  \
                             ta1_y_zzzz_z_0,  \
                             ta1_y_zzzz_z_1,  \
                             ta1_y_zzzzz_x_0, \
                             ta1_y_zzzzz_y_0, \
                             ta1_y_zzzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzz_x_0[i] =
            4.0 * ta1_y_zzz_x_0[i] * fe_0 - 4.0 * ta1_y_zzz_x_1[i] * fe_0 + ta1_y_zzzz_x_0[i] * pa_z[i] - ta1_y_zzzz_x_1[i] * pc_z[i];

        ta1_y_zzzzz_y_0[i] =
            4.0 * ta1_y_zzz_y_0[i] * fe_0 - 4.0 * ta1_y_zzz_y_1[i] * fe_0 + ta1_y_zzzz_y_0[i] * pa_z[i] - ta1_y_zzzz_y_1[i] * pc_z[i];

        ta1_y_zzzzz_z_0[i] = 4.0 * ta1_y_zzz_z_0[i] * fe_0 - 4.0 * ta1_y_zzz_z_1[i] * fe_0 + ta1_y_zzzz_0_0[i] * fe_0 - ta1_y_zzzz_0_1[i] * fe_0 +
                             ta1_y_zzzz_z_0[i] * pa_z[i] - ta1_y_zzzz_z_1[i] * pc_z[i];
    }

    // Set up 126-129 components of targeted buffer : HP

    auto ta1_z_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 126);

    auto ta1_z_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 127);

    auto ta1_z_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 128);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xxx_x_0,   \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxx_y_0,   \
                             ta1_z_xxx_y_1,   \
                             ta1_z_xxx_z_0,   \
                             ta1_z_xxx_z_1,   \
                             ta1_z_xxxx_0_0,  \
                             ta1_z_xxxx_0_1,  \
                             ta1_z_xxxx_x_0,  \
                             ta1_z_xxxx_x_1,  \
                             ta1_z_xxxx_y_0,  \
                             ta1_z_xxxx_y_1,  \
                             ta1_z_xxxx_z_0,  \
                             ta1_z_xxxx_z_1,  \
                             ta1_z_xxxxx_x_0, \
                             ta1_z_xxxxx_y_0, \
                             ta1_z_xxxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxx_x_0[i] = 4.0 * ta1_z_xxx_x_0[i] * fe_0 - 4.0 * ta1_z_xxx_x_1[i] * fe_0 + ta1_z_xxxx_0_0[i] * fe_0 - ta1_z_xxxx_0_1[i] * fe_0 +
                             ta1_z_xxxx_x_0[i] * pa_x[i] - ta1_z_xxxx_x_1[i] * pc_x[i];

        ta1_z_xxxxx_y_0[i] =
            4.0 * ta1_z_xxx_y_0[i] * fe_0 - 4.0 * ta1_z_xxx_y_1[i] * fe_0 + ta1_z_xxxx_y_0[i] * pa_x[i] - ta1_z_xxxx_y_1[i] * pc_x[i];

        ta1_z_xxxxx_z_0[i] =
            4.0 * ta1_z_xxx_z_0[i] * fe_0 - 4.0 * ta1_z_xxx_z_1[i] * fe_0 + ta1_z_xxxx_z_0[i] * pa_x[i] - ta1_z_xxxx_z_1[i] * pc_x[i];
    }

    // Set up 129-132 components of targeted buffer : HP

    auto ta1_z_xxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 129);

    auto ta1_z_xxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 130);

    auto ta1_z_xxxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 131);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xxxx_x_0,  \
                             ta1_z_xxxx_x_1,  \
                             ta1_z_xxxx_z_0,  \
                             ta1_z_xxxx_z_1,  \
                             ta1_z_xxxxy_x_0, \
                             ta1_z_xxxxy_y_0, \
                             ta1_z_xxxxy_z_0, \
                             ta1_z_xxxy_y_0,  \
                             ta1_z_xxxy_y_1,  \
                             ta1_z_xxy_y_0,   \
                             ta1_z_xxy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxy_x_0[i] = ta1_z_xxxx_x_0[i] * pa_y[i] - ta1_z_xxxx_x_1[i] * pc_y[i];

        ta1_z_xxxxy_y_0[i] =
            3.0 * ta1_z_xxy_y_0[i] * fe_0 - 3.0 * ta1_z_xxy_y_1[i] * fe_0 + ta1_z_xxxy_y_0[i] * pa_x[i] - ta1_z_xxxy_y_1[i] * pc_x[i];

        ta1_z_xxxxy_z_0[i] = ta1_z_xxxx_z_0[i] * pa_y[i] - ta1_z_xxxx_z_1[i] * pc_y[i];
    }

    // Set up 132-135 components of targeted buffer : HP

    auto ta1_z_xxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 132);

    auto ta1_z_xxxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 133);

    auto ta1_z_xxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 134);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xxxx_x_0,  \
                             ta1_z_xxxx_x_1,  \
                             ta1_z_xxxx_y_0,  \
                             ta1_z_xxxx_y_1,  \
                             ta1_z_xxxxz_x_0, \
                             ta1_z_xxxxz_y_0, \
                             ta1_z_xxxxz_z_0, \
                             ta1_z_xxxz_z_0,  \
                             ta1_z_xxxz_z_1,  \
                             ta1_z_xxz_z_0,   \
                             ta1_z_xxz_z_1,   \
                             ta_xxxx_x_1,     \
                             ta_xxxx_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxz_x_0[i] = ta_xxxx_x_1[i] + ta1_z_xxxx_x_0[i] * pa_z[i] - ta1_z_xxxx_x_1[i] * pc_z[i];

        ta1_z_xxxxz_y_0[i] = ta_xxxx_y_1[i] + ta1_z_xxxx_y_0[i] * pa_z[i] - ta1_z_xxxx_y_1[i] * pc_z[i];

        ta1_z_xxxxz_z_0[i] =
            3.0 * ta1_z_xxz_z_0[i] * fe_0 - 3.0 * ta1_z_xxz_z_1[i] * fe_0 + ta1_z_xxxz_z_0[i] * pa_x[i] - ta1_z_xxxz_z_1[i] * pc_x[i];
    }

    // Set up 135-138 components of targeted buffer : HP

    auto ta1_z_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 135);

    auto ta1_z_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 136);

    auto ta1_z_xxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 137);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xxx_x_0,   \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxxy_x_0,  \
                             ta1_z_xxxy_x_1,  \
                             ta1_z_xxxyy_x_0, \
                             ta1_z_xxxyy_y_0, \
                             ta1_z_xxxyy_z_0, \
                             ta1_z_xxyy_y_0,  \
                             ta1_z_xxyy_y_1,  \
                             ta1_z_xxyy_z_0,  \
                             ta1_z_xxyy_z_1,  \
                             ta1_z_xyy_y_0,   \
                             ta1_z_xyy_y_1,   \
                             ta1_z_xyy_z_0,   \
                             ta1_z_xyy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyy_x_0[i] = ta1_z_xxx_x_0[i] * fe_0 - ta1_z_xxx_x_1[i] * fe_0 + ta1_z_xxxy_x_0[i] * pa_y[i] - ta1_z_xxxy_x_1[i] * pc_y[i];

        ta1_z_xxxyy_y_0[i] =
            2.0 * ta1_z_xyy_y_0[i] * fe_0 - 2.0 * ta1_z_xyy_y_1[i] * fe_0 + ta1_z_xxyy_y_0[i] * pa_x[i] - ta1_z_xxyy_y_1[i] * pc_x[i];

        ta1_z_xxxyy_z_0[i] =
            2.0 * ta1_z_xyy_z_0[i] * fe_0 - 2.0 * ta1_z_xyy_z_1[i] * fe_0 + ta1_z_xxyy_z_0[i] * pa_x[i] - ta1_z_xxyy_z_1[i] * pc_x[i];
    }

    // Set up 138-141 components of targeted buffer : HP

    auto ta1_z_xxxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 138);

    auto ta1_z_xxxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 139);

    auto ta1_z_xxxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 140);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_xxxy_y_0,  \
                             ta1_z_xxxy_y_1,  \
                             ta1_z_xxxyz_x_0, \
                             ta1_z_xxxyz_y_0, \
                             ta1_z_xxxyz_z_0, \
                             ta1_z_xxxz_x_0,  \
                             ta1_z_xxxz_x_1,  \
                             ta1_z_xxxz_z_0,  \
                             ta1_z_xxxz_z_1,  \
                             ta_xxxy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xxxyz_x_0[i] = ta1_z_xxxz_x_0[i] * pa_y[i] - ta1_z_xxxz_x_1[i] * pc_y[i];

        ta1_z_xxxyz_y_0[i] = ta_xxxy_y_1[i] + ta1_z_xxxy_y_0[i] * pa_z[i] - ta1_z_xxxy_y_1[i] * pc_z[i];

        ta1_z_xxxyz_z_0[i] = ta1_z_xxxz_z_0[i] * pa_y[i] - ta1_z_xxxz_z_1[i] * pc_y[i];
    }

    // Set up 141-144 components of targeted buffer : HP

    auto ta1_z_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 141);

    auto ta1_z_xxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 142);

    auto ta1_z_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 143);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xxx_x_0,   \
                             ta1_z_xxx_x_1,   \
                             ta1_z_xxxz_x_0,  \
                             ta1_z_xxxz_x_1,  \
                             ta1_z_xxxzz_x_0, \
                             ta1_z_xxxzz_y_0, \
                             ta1_z_xxxzz_z_0, \
                             ta1_z_xxzz_y_0,  \
                             ta1_z_xxzz_y_1,  \
                             ta1_z_xxzz_z_0,  \
                             ta1_z_xxzz_z_1,  \
                             ta1_z_xzz_y_0,   \
                             ta1_z_xzz_y_1,   \
                             ta1_z_xzz_z_0,   \
                             ta1_z_xzz_z_1,   \
                             ta_xxxz_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzz_x_0[i] =
            ta1_z_xxx_x_0[i] * fe_0 - ta1_z_xxx_x_1[i] * fe_0 + ta_xxxz_x_1[i] + ta1_z_xxxz_x_0[i] * pa_z[i] - ta1_z_xxxz_x_1[i] * pc_z[i];

        ta1_z_xxxzz_y_0[i] =
            2.0 * ta1_z_xzz_y_0[i] * fe_0 - 2.0 * ta1_z_xzz_y_1[i] * fe_0 + ta1_z_xxzz_y_0[i] * pa_x[i] - ta1_z_xxzz_y_1[i] * pc_x[i];

        ta1_z_xxxzz_z_0[i] =
            2.0 * ta1_z_xzz_z_0[i] * fe_0 - 2.0 * ta1_z_xzz_z_1[i] * fe_0 + ta1_z_xxzz_z_0[i] * pa_x[i] - ta1_z_xxzz_z_1[i] * pc_x[i];
    }

    // Set up 144-147 components of targeted buffer : HP

    auto ta1_z_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 144);

    auto ta1_z_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 145);

    auto ta1_z_xxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 146);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xxy_x_0,   \
                             ta1_z_xxy_x_1,   \
                             ta1_z_xxyy_x_0,  \
                             ta1_z_xxyy_x_1,  \
                             ta1_z_xxyyy_x_0, \
                             ta1_z_xxyyy_y_0, \
                             ta1_z_xxyyy_z_0, \
                             ta1_z_xyyy_y_0,  \
                             ta1_z_xyyy_y_1,  \
                             ta1_z_xyyy_z_0,  \
                             ta1_z_xyyy_z_1,  \
                             ta1_z_yyy_y_0,   \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyy_z_0,   \
                             ta1_z_yyy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyy_x_0[i] =
            2.0 * ta1_z_xxy_x_0[i] * fe_0 - 2.0 * ta1_z_xxy_x_1[i] * fe_0 + ta1_z_xxyy_x_0[i] * pa_y[i] - ta1_z_xxyy_x_1[i] * pc_y[i];

        ta1_z_xxyyy_y_0[i] = ta1_z_yyy_y_0[i] * fe_0 - ta1_z_yyy_y_1[i] * fe_0 + ta1_z_xyyy_y_0[i] * pa_x[i] - ta1_z_xyyy_y_1[i] * pc_x[i];

        ta1_z_xxyyy_z_0[i] = ta1_z_yyy_z_0[i] * fe_0 - ta1_z_yyy_z_1[i] * fe_0 + ta1_z_xyyy_z_0[i] * pa_x[i] - ta1_z_xyyy_z_1[i] * pc_x[i];
    }

    // Set up 147-150 components of targeted buffer : HP

    auto ta1_z_xxyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 147);

    auto ta1_z_xxyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 148);

    auto ta1_z_xxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 149);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xxyy_x_0,  \
                             ta1_z_xxyy_x_1,  \
                             ta1_z_xxyy_y_0,  \
                             ta1_z_xxyy_y_1,  \
                             ta1_z_xxyyz_x_0, \
                             ta1_z_xxyyz_y_0, \
                             ta1_z_xxyyz_z_0, \
                             ta1_z_xyyz_z_0,  \
                             ta1_z_xyyz_z_1,  \
                             ta1_z_yyz_z_0,   \
                             ta1_z_yyz_z_1,   \
                             ta_xxyy_x_1,     \
                             ta_xxyy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyz_x_0[i] = ta_xxyy_x_1[i] + ta1_z_xxyy_x_0[i] * pa_z[i] - ta1_z_xxyy_x_1[i] * pc_z[i];

        ta1_z_xxyyz_y_0[i] = ta_xxyy_y_1[i] + ta1_z_xxyy_y_0[i] * pa_z[i] - ta1_z_xxyy_y_1[i] * pc_z[i];

        ta1_z_xxyyz_z_0[i] = ta1_z_yyz_z_0[i] * fe_0 - ta1_z_yyz_z_1[i] * fe_0 + ta1_z_xyyz_z_0[i] * pa_x[i] - ta1_z_xyyz_z_1[i] * pc_x[i];
    }

    // Set up 150-153 components of targeted buffer : HP

    auto ta1_z_xxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 150);

    auto ta1_z_xxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 151);

    auto ta1_z_xxyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 152);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xxyzz_x_0, \
                             ta1_z_xxyzz_y_0, \
                             ta1_z_xxyzz_z_0, \
                             ta1_z_xxzz_x_0,  \
                             ta1_z_xxzz_x_1,  \
                             ta1_z_xxzz_z_0,  \
                             ta1_z_xxzz_z_1,  \
                             ta1_z_xyzz_y_0,  \
                             ta1_z_xyzz_y_1,  \
                             ta1_z_yzz_y_0,   \
                             ta1_z_yzz_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzz_x_0[i] = ta1_z_xxzz_x_0[i] * pa_y[i] - ta1_z_xxzz_x_1[i] * pc_y[i];

        ta1_z_xxyzz_y_0[i] = ta1_z_yzz_y_0[i] * fe_0 - ta1_z_yzz_y_1[i] * fe_0 + ta1_z_xyzz_y_0[i] * pa_x[i] - ta1_z_xyzz_y_1[i] * pc_x[i];

        ta1_z_xxyzz_z_0[i] = ta1_z_xxzz_z_0[i] * pa_y[i] - ta1_z_xxzz_z_1[i] * pc_y[i];
    }

    // Set up 153-156 components of targeted buffer : HP

    auto ta1_z_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 153);

    auto ta1_z_xxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 154);

    auto ta1_z_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 155);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xxz_x_0,   \
                             ta1_z_xxz_x_1,   \
                             ta1_z_xxzz_x_0,  \
                             ta1_z_xxzz_x_1,  \
                             ta1_z_xxzzz_x_0, \
                             ta1_z_xxzzz_y_0, \
                             ta1_z_xxzzz_z_0, \
                             ta1_z_xzzz_y_0,  \
                             ta1_z_xzzz_y_1,  \
                             ta1_z_xzzz_z_0,  \
                             ta1_z_xzzz_z_1,  \
                             ta1_z_zzz_y_0,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_z_0,   \
                             ta1_z_zzz_z_1,   \
                             ta_xxzz_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzz_x_0[i] = 2.0 * ta1_z_xxz_x_0[i] * fe_0 - 2.0 * ta1_z_xxz_x_1[i] * fe_0 + ta_xxzz_x_1[i] + ta1_z_xxzz_x_0[i] * pa_z[i] -
                             ta1_z_xxzz_x_1[i] * pc_z[i];

        ta1_z_xxzzz_y_0[i] = ta1_z_zzz_y_0[i] * fe_0 - ta1_z_zzz_y_1[i] * fe_0 + ta1_z_xzzz_y_0[i] * pa_x[i] - ta1_z_xzzz_y_1[i] * pc_x[i];

        ta1_z_xxzzz_z_0[i] = ta1_z_zzz_z_0[i] * fe_0 - ta1_z_zzz_z_1[i] * fe_0 + ta1_z_xzzz_z_0[i] * pa_x[i] - ta1_z_xzzz_z_1[i] * pc_x[i];
    }

    // Set up 156-159 components of targeted buffer : HP

    auto ta1_z_xyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 156);

    auto ta1_z_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 157);

    auto ta1_z_xyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 158);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xyyyy_x_0, \
                             ta1_z_xyyyy_y_0, \
                             ta1_z_xyyyy_z_0, \
                             ta1_z_yyyy_0_0,  \
                             ta1_z_yyyy_0_1,  \
                             ta1_z_yyyy_x_0,  \
                             ta1_z_yyyy_x_1,  \
                             ta1_z_yyyy_y_0,  \
                             ta1_z_yyyy_y_1,  \
                             ta1_z_yyyy_z_0,  \
                             ta1_z_yyyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyy_x_0[i] = ta1_z_yyyy_0_0[i] * fe_0 - ta1_z_yyyy_0_1[i] * fe_0 + ta1_z_yyyy_x_0[i] * pa_x[i] - ta1_z_yyyy_x_1[i] * pc_x[i];

        ta1_z_xyyyy_y_0[i] = ta1_z_yyyy_y_0[i] * pa_x[i] - ta1_z_yyyy_y_1[i] * pc_x[i];

        ta1_z_xyyyy_z_0[i] = ta1_z_yyyy_z_0[i] * pa_x[i] - ta1_z_yyyy_z_1[i] * pc_x[i];
    }

    // Set up 159-162 components of targeted buffer : HP

    auto ta1_z_xyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 159);

    auto ta1_z_xyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 160);

    auto ta1_z_xyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 161);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_xyyy_x_0,  \
                             ta1_z_xyyy_x_1,  \
                             ta1_z_xyyyz_x_0, \
                             ta1_z_xyyyz_y_0, \
                             ta1_z_xyyyz_z_0, \
                             ta1_z_yyyz_y_0,  \
                             ta1_z_yyyz_y_1,  \
                             ta1_z_yyyz_z_0,  \
                             ta1_z_yyyz_z_1,  \
                             ta_xyyy_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyyyz_x_0[i] = ta_xyyy_x_1[i] + ta1_z_xyyy_x_0[i] * pa_z[i] - ta1_z_xyyy_x_1[i] * pc_z[i];

        ta1_z_xyyyz_y_0[i] = ta1_z_yyyz_y_0[i] * pa_x[i] - ta1_z_yyyz_y_1[i] * pc_x[i];

        ta1_z_xyyyz_z_0[i] = ta1_z_yyyz_z_0[i] * pa_x[i] - ta1_z_yyyz_z_1[i] * pc_x[i];
    }

    // Set up 162-165 components of targeted buffer : HP

    auto ta1_z_xyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 162);

    auto ta1_z_xyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 163);

    auto ta1_z_xyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 164);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xyyzz_x_0, \
                             ta1_z_xyyzz_y_0, \
                             ta1_z_xyyzz_z_0, \
                             ta1_z_yyzz_0_0,  \
                             ta1_z_yyzz_0_1,  \
                             ta1_z_yyzz_x_0,  \
                             ta1_z_yyzz_x_1,  \
                             ta1_z_yyzz_y_0,  \
                             ta1_z_yyzz_y_1,  \
                             ta1_z_yyzz_z_0,  \
                             ta1_z_yyzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzz_x_0[i] = ta1_z_yyzz_0_0[i] * fe_0 - ta1_z_yyzz_0_1[i] * fe_0 + ta1_z_yyzz_x_0[i] * pa_x[i] - ta1_z_yyzz_x_1[i] * pc_x[i];

        ta1_z_xyyzz_y_0[i] = ta1_z_yyzz_y_0[i] * pa_x[i] - ta1_z_yyzz_y_1[i] * pc_x[i];

        ta1_z_xyyzz_z_0[i] = ta1_z_yyzz_z_0[i] * pa_x[i] - ta1_z_yyzz_z_1[i] * pc_x[i];
    }

    // Set up 165-168 components of targeted buffer : HP

    auto ta1_z_xyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 165);

    auto ta1_z_xyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 166);

    auto ta1_z_xyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 167);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_xyzzz_x_0, \
                             ta1_z_xyzzz_y_0, \
                             ta1_z_xyzzz_z_0, \
                             ta1_z_xzzz_x_0,  \
                             ta1_z_xzzz_x_1,  \
                             ta1_z_yzzz_y_0,  \
                             ta1_z_yzzz_y_1,  \
                             ta1_z_yzzz_z_0,  \
                             ta1_z_yzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyzzz_x_0[i] = ta1_z_xzzz_x_0[i] * pa_y[i] - ta1_z_xzzz_x_1[i] * pc_y[i];

        ta1_z_xyzzz_y_0[i] = ta1_z_yzzz_y_0[i] * pa_x[i] - ta1_z_yzzz_y_1[i] * pc_x[i];

        ta1_z_xyzzz_z_0[i] = ta1_z_yzzz_z_0[i] * pa_x[i] - ta1_z_yzzz_z_1[i] * pc_x[i];
    }

    // Set up 168-171 components of targeted buffer : HP

    auto ta1_z_xzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 168);

    auto ta1_z_xzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 169);

    auto ta1_z_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 170);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_xzzzz_x_0, \
                             ta1_z_xzzzz_y_0, \
                             ta1_z_xzzzz_z_0, \
                             ta1_z_zzzz_0_0,  \
                             ta1_z_zzzz_0_1,  \
                             ta1_z_zzzz_x_0,  \
                             ta1_z_zzzz_x_1,  \
                             ta1_z_zzzz_y_0,  \
                             ta1_z_zzzz_y_1,  \
                             ta1_z_zzzz_z_0,  \
                             ta1_z_zzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzz_x_0[i] = ta1_z_zzzz_0_0[i] * fe_0 - ta1_z_zzzz_0_1[i] * fe_0 + ta1_z_zzzz_x_0[i] * pa_x[i] - ta1_z_zzzz_x_1[i] * pc_x[i];

        ta1_z_xzzzz_y_0[i] = ta1_z_zzzz_y_0[i] * pa_x[i] - ta1_z_zzzz_y_1[i] * pc_x[i];

        ta1_z_xzzzz_z_0[i] = ta1_z_zzzz_z_0[i] * pa_x[i] - ta1_z_zzzz_z_1[i] * pc_x[i];
    }

    // Set up 171-174 components of targeted buffer : HP

    auto ta1_z_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 171);

    auto ta1_z_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 172);

    auto ta1_z_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 173);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_yyy_x_0,   \
                             ta1_z_yyy_x_1,   \
                             ta1_z_yyy_y_0,   \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyy_z_0,   \
                             ta1_z_yyy_z_1,   \
                             ta1_z_yyyy_0_0,  \
                             ta1_z_yyyy_0_1,  \
                             ta1_z_yyyy_x_0,  \
                             ta1_z_yyyy_x_1,  \
                             ta1_z_yyyy_y_0,  \
                             ta1_z_yyyy_y_1,  \
                             ta1_z_yyyy_z_0,  \
                             ta1_z_yyyy_z_1,  \
                             ta1_z_yyyyy_x_0, \
                             ta1_z_yyyyy_y_0, \
                             ta1_z_yyyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyy_x_0[i] =
            4.0 * ta1_z_yyy_x_0[i] * fe_0 - 4.0 * ta1_z_yyy_x_1[i] * fe_0 + ta1_z_yyyy_x_0[i] * pa_y[i] - ta1_z_yyyy_x_1[i] * pc_y[i];

        ta1_z_yyyyy_y_0[i] = 4.0 * ta1_z_yyy_y_0[i] * fe_0 - 4.0 * ta1_z_yyy_y_1[i] * fe_0 + ta1_z_yyyy_0_0[i] * fe_0 - ta1_z_yyyy_0_1[i] * fe_0 +
                             ta1_z_yyyy_y_0[i] * pa_y[i] - ta1_z_yyyy_y_1[i] * pc_y[i];

        ta1_z_yyyyy_z_0[i] =
            4.0 * ta1_z_yyy_z_0[i] * fe_0 - 4.0 * ta1_z_yyy_z_1[i] * fe_0 + ta1_z_yyyy_z_0[i] * pa_y[i] - ta1_z_yyyy_z_1[i] * pc_y[i];
    }

    // Set up 174-177 components of targeted buffer : HP

    auto ta1_z_yyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 174);

    auto ta1_z_yyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 175);

    auto ta1_z_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 176);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yyyy_x_0,  \
                             ta1_z_yyyy_x_1,  \
                             ta1_z_yyyy_y_0,  \
                             ta1_z_yyyy_y_1,  \
                             ta1_z_yyyyz_x_0, \
                             ta1_z_yyyyz_y_0, \
                             ta1_z_yyyyz_z_0, \
                             ta1_z_yyyz_z_0,  \
                             ta1_z_yyyz_z_1,  \
                             ta1_z_yyz_z_0,   \
                             ta1_z_yyz_z_1,   \
                             ta_yyyy_x_1,     \
                             ta_yyyy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyz_x_0[i] = ta_yyyy_x_1[i] + ta1_z_yyyy_x_0[i] * pa_z[i] - ta1_z_yyyy_x_1[i] * pc_z[i];

        ta1_z_yyyyz_y_0[i] = ta_yyyy_y_1[i] + ta1_z_yyyy_y_0[i] * pa_z[i] - ta1_z_yyyy_y_1[i] * pc_z[i];

        ta1_z_yyyyz_z_0[i] =
            3.0 * ta1_z_yyz_z_0[i] * fe_0 - 3.0 * ta1_z_yyz_z_1[i] * fe_0 + ta1_z_yyyz_z_0[i] * pa_y[i] - ta1_z_yyyz_z_1[i] * pc_y[i];
    }

    // Set up 177-180 components of targeted buffer : HP

    auto ta1_z_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 177);

    auto ta1_z_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 178);

    auto ta1_z_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 179);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yyy_y_0,   \
                             ta1_z_yyy_y_1,   \
                             ta1_z_yyyz_y_0,  \
                             ta1_z_yyyz_y_1,  \
                             ta1_z_yyyzz_x_0, \
                             ta1_z_yyyzz_y_0, \
                             ta1_z_yyyzz_z_0, \
                             ta1_z_yyzz_x_0,  \
                             ta1_z_yyzz_x_1,  \
                             ta1_z_yyzz_z_0,  \
                             ta1_z_yyzz_z_1,  \
                             ta1_z_yzz_x_0,   \
                             ta1_z_yzz_x_1,   \
                             ta1_z_yzz_z_0,   \
                             ta1_z_yzz_z_1,   \
                             ta_yyyz_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzz_x_0[i] =
            2.0 * ta1_z_yzz_x_0[i] * fe_0 - 2.0 * ta1_z_yzz_x_1[i] * fe_0 + ta1_z_yyzz_x_0[i] * pa_y[i] - ta1_z_yyzz_x_1[i] * pc_y[i];

        ta1_z_yyyzz_y_0[i] =
            ta1_z_yyy_y_0[i] * fe_0 - ta1_z_yyy_y_1[i] * fe_0 + ta_yyyz_y_1[i] + ta1_z_yyyz_y_0[i] * pa_z[i] - ta1_z_yyyz_y_1[i] * pc_z[i];

        ta1_z_yyyzz_z_0[i] =
            2.0 * ta1_z_yzz_z_0[i] * fe_0 - 2.0 * ta1_z_yzz_z_1[i] * fe_0 + ta1_z_yyzz_z_0[i] * pa_y[i] - ta1_z_yyzz_z_1[i] * pc_y[i];
    }

    // Set up 180-183 components of targeted buffer : HP

    auto ta1_z_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 180);

    auto ta1_z_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 181);

    auto ta1_z_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 182);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_yyz_y_0,   \
                             ta1_z_yyz_y_1,   \
                             ta1_z_yyzz_y_0,  \
                             ta1_z_yyzz_y_1,  \
                             ta1_z_yyzzz_x_0, \
                             ta1_z_yyzzz_y_0, \
                             ta1_z_yyzzz_z_0, \
                             ta1_z_yzzz_x_0,  \
                             ta1_z_yzzz_x_1,  \
                             ta1_z_yzzz_z_0,  \
                             ta1_z_yzzz_z_1,  \
                             ta1_z_zzz_x_0,   \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_z_0,   \
                             ta1_z_zzz_z_1,   \
                             ta_yyzz_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzz_x_0[i] = ta1_z_zzz_x_0[i] * fe_0 - ta1_z_zzz_x_1[i] * fe_0 + ta1_z_yzzz_x_0[i] * pa_y[i] - ta1_z_yzzz_x_1[i] * pc_y[i];

        ta1_z_yyzzz_y_0[i] = 2.0 * ta1_z_yyz_y_0[i] * fe_0 - 2.0 * ta1_z_yyz_y_1[i] * fe_0 + ta_yyzz_y_1[i] + ta1_z_yyzz_y_0[i] * pa_z[i] -
                             ta1_z_yyzz_y_1[i] * pc_z[i];

        ta1_z_yyzzz_z_0[i] = ta1_z_zzz_z_0[i] * fe_0 - ta1_z_zzz_z_1[i] * fe_0 + ta1_z_yzzz_z_0[i] * pa_y[i] - ta1_z_yzzz_z_1[i] * pc_y[i];
    }

    // Set up 183-186 components of targeted buffer : HP

    auto ta1_z_yzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 183);

    auto ta1_z_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 184);

    auto ta1_z_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 185);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_yzzzz_x_0, \
                             ta1_z_yzzzz_y_0, \
                             ta1_z_yzzzz_z_0, \
                             ta1_z_zzzz_0_0,  \
                             ta1_z_zzzz_0_1,  \
                             ta1_z_zzzz_x_0,  \
                             ta1_z_zzzz_x_1,  \
                             ta1_z_zzzz_y_0,  \
                             ta1_z_zzzz_y_1,  \
                             ta1_z_zzzz_z_0,  \
                             ta1_z_zzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzz_x_0[i] = ta1_z_zzzz_x_0[i] * pa_y[i] - ta1_z_zzzz_x_1[i] * pc_y[i];

        ta1_z_yzzzz_y_0[i] = ta1_z_zzzz_0_0[i] * fe_0 - ta1_z_zzzz_0_1[i] * fe_0 + ta1_z_zzzz_y_0[i] * pa_y[i] - ta1_z_zzzz_y_1[i] * pc_y[i];

        ta1_z_yzzzz_z_0[i] = ta1_z_zzzz_z_0[i] * pa_y[i] - ta1_z_zzzz_z_1[i] * pc_y[i];
    }

    // Set up 186-189 components of targeted buffer : HP

    auto ta1_z_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 186);

    auto ta1_z_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 187);

    auto ta1_z_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 188);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_z_zzz_x_0,   \
                             ta1_z_zzz_x_1,   \
                             ta1_z_zzz_y_0,   \
                             ta1_z_zzz_y_1,   \
                             ta1_z_zzz_z_0,   \
                             ta1_z_zzz_z_1,   \
                             ta1_z_zzzz_0_0,  \
                             ta1_z_zzzz_0_1,  \
                             ta1_z_zzzz_x_0,  \
                             ta1_z_zzzz_x_1,  \
                             ta1_z_zzzz_y_0,  \
                             ta1_z_zzzz_y_1,  \
                             ta1_z_zzzz_z_0,  \
                             ta1_z_zzzz_z_1,  \
                             ta1_z_zzzzz_x_0, \
                             ta1_z_zzzzz_y_0, \
                             ta1_z_zzzzz_z_0, \
                             ta_zzzz_x_1,     \
                             ta_zzzz_y_1,     \
                             ta_zzzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzz_x_0[i] = 4.0 * ta1_z_zzz_x_0[i] * fe_0 - 4.0 * ta1_z_zzz_x_1[i] * fe_0 + ta_zzzz_x_1[i] + ta1_z_zzzz_x_0[i] * pa_z[i] -
                             ta1_z_zzzz_x_1[i] * pc_z[i];

        ta1_z_zzzzz_y_0[i] = 4.0 * ta1_z_zzz_y_0[i] * fe_0 - 4.0 * ta1_z_zzz_y_1[i] * fe_0 + ta_zzzz_y_1[i] + ta1_z_zzzz_y_0[i] * pa_z[i] -
                             ta1_z_zzzz_y_1[i] * pc_z[i];

        ta1_z_zzzzz_z_0[i] = 4.0 * ta1_z_zzz_z_0[i] * fe_0 - 4.0 * ta1_z_zzz_z_1[i] * fe_0 + ta1_z_zzzz_0_0[i] * fe_0 - ta1_z_zzzz_0_1[i] * fe_0 +
                             ta_zzzz_z_1[i] + ta1_z_zzzz_z_0[i] * pa_z[i] - ta1_z_zzzz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
