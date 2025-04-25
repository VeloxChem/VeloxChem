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

#include "NuclearPotentialGeom010PrimRecPG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_pg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_pg,
                                        const size_t              idx_npot_geom_010_0_sf,
                                        const size_t              idx_npot_geom_010_1_sf,
                                        const size_t              idx_npot_1_sg,
                                        const size_t              idx_npot_geom_010_0_sg,
                                        const size_t              idx_npot_geom_010_1_sg,
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

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf);

    auto ta1_x_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 1);

    auto ta1_x_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 2);

    auto ta1_x_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 3);

    auto ta1_x_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 4);

    auto ta1_x_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 5);

    auto ta1_x_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 6);

    auto ta1_x_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 7);

    auto ta1_x_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 8);

    auto ta1_x_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 9);

    auto ta1_y_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 10);

    auto ta1_y_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 11);

    auto ta1_y_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 12);

    auto ta1_y_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 13);

    auto ta1_y_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 14);

    auto ta1_y_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 15);

    auto ta1_y_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 16);

    auto ta1_y_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 17);

    auto ta1_y_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 18);

    auto ta1_y_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 19);

    auto ta1_z_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 20);

    auto ta1_z_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 21);

    auto ta1_z_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 22);

    auto ta1_z_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 23);

    auto ta1_z_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 24);

    auto ta1_z_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 25);

    auto ta1_z_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 26);

    auto ta1_z_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 27);

    auto ta1_z_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 28);

    auto ta1_z_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 29);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf);

    auto ta1_x_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 1);

    auto ta1_x_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 2);

    auto ta1_x_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 3);

    auto ta1_x_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 4);

    auto ta1_x_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 5);

    auto ta1_x_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 6);

    auto ta1_x_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 7);

    auto ta1_x_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 8);

    auto ta1_x_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 9);

    auto ta1_y_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 10);

    auto ta1_y_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 11);

    auto ta1_y_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 12);

    auto ta1_y_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 13);

    auto ta1_y_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 14);

    auto ta1_y_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 15);

    auto ta1_y_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 16);

    auto ta1_y_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 17);

    auto ta1_y_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 18);

    auto ta1_y_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 19);

    auto ta1_z_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 20);

    auto ta1_z_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 21);

    auto ta1_z_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 22);

    auto ta1_z_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 23);

    auto ta1_z_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 24);

    auto ta1_z_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 25);

    auto ta1_z_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 26);

    auto ta1_z_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 27);

    auto ta1_z_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 28);

    auto ta1_z_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 29);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_1 = pbuffer.data(idx_npot_1_sg);

    auto ta_0_xxxy_1 = pbuffer.data(idx_npot_1_sg + 1);

    auto ta_0_xxxz_1 = pbuffer.data(idx_npot_1_sg + 2);

    auto ta_0_xxyy_1 = pbuffer.data(idx_npot_1_sg + 3);

    auto ta_0_xxyz_1 = pbuffer.data(idx_npot_1_sg + 4);

    auto ta_0_xxzz_1 = pbuffer.data(idx_npot_1_sg + 5);

    auto ta_0_xyyy_1 = pbuffer.data(idx_npot_1_sg + 6);

    auto ta_0_xyyz_1 = pbuffer.data(idx_npot_1_sg + 7);

    auto ta_0_xyzz_1 = pbuffer.data(idx_npot_1_sg + 8);

    auto ta_0_xzzz_1 = pbuffer.data(idx_npot_1_sg + 9);

    auto ta_0_yyyy_1 = pbuffer.data(idx_npot_1_sg + 10);

    auto ta_0_yyyz_1 = pbuffer.data(idx_npot_1_sg + 11);

    auto ta_0_yyzz_1 = pbuffer.data(idx_npot_1_sg + 12);

    auto ta_0_yzzz_1 = pbuffer.data(idx_npot_1_sg + 13);

    auto ta_0_zzzz_1 = pbuffer.data(idx_npot_1_sg + 14);

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

    // Set up 0-15 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_x_0_xxx_0,  \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxxx_0, \
                             ta1_x_0_xxxx_1, \
                             ta1_x_0_xxxy_0, \
                             ta1_x_0_xxxy_1, \
                             ta1_x_0_xxxz_0, \
                             ta1_x_0_xxxz_1, \
                             ta1_x_0_xxy_0,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxyy_0, \
                             ta1_x_0_xxyy_1, \
                             ta1_x_0_xxyz_0, \
                             ta1_x_0_xxyz_1, \
                             ta1_x_0_xxz_0,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xxzz_0, \
                             ta1_x_0_xxzz_1, \
                             ta1_x_0_xyy_0,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyyy_0, \
                             ta1_x_0_xyyy_1, \
                             ta1_x_0_xyyz_0, \
                             ta1_x_0_xyyz_1, \
                             ta1_x_0_xyz_0,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xyzz_0, \
                             ta1_x_0_xyzz_1, \
                             ta1_x_0_xzz_0,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_xzzz_0, \
                             ta1_x_0_xzzz_1, \
                             ta1_x_0_yyy_0,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyyy_0, \
                             ta1_x_0_yyyy_1, \
                             ta1_x_0_yyyz_0, \
                             ta1_x_0_yyyz_1, \
                             ta1_x_0_yyz_0,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yyzz_0, \
                             ta1_x_0_yyzz_1, \
                             ta1_x_0_yzz_0,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_yzzz_0, \
                             ta1_x_0_yzzz_1, \
                             ta1_x_0_zzz_0,  \
                             ta1_x_0_zzz_1,  \
                             ta1_x_0_zzzz_0, \
                             ta1_x_0_zzzz_1, \
                             ta1_x_x_xxxx_0, \
                             ta1_x_x_xxxy_0, \
                             ta1_x_x_xxxz_0, \
                             ta1_x_x_xxyy_0, \
                             ta1_x_x_xxyz_0, \
                             ta1_x_x_xxzz_0, \
                             ta1_x_x_xyyy_0, \
                             ta1_x_x_xyyz_0, \
                             ta1_x_x_xyzz_0, \
                             ta1_x_x_xzzz_0, \
                             ta1_x_x_yyyy_0, \
                             ta1_x_x_yyyz_0, \
                             ta1_x_x_yyzz_0, \
                             ta1_x_x_yzzz_0, \
                             ta1_x_x_zzzz_0, \
                             ta_0_xxxx_1,    \
                             ta_0_xxxy_1,    \
                             ta_0_xxxz_1,    \
                             ta_0_xxyy_1,    \
                             ta_0_xxyz_1,    \
                             ta_0_xxzz_1,    \
                             ta_0_xyyy_1,    \
                             ta_0_xyyz_1,    \
                             ta_0_xyzz_1,    \
                             ta_0_xzzz_1,    \
                             ta_0_yyyy_1,    \
                             ta_0_yyyz_1,    \
                             ta_0_yyzz_1,    \
                             ta_0_yzzz_1,    \
                             ta_0_zzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_x_xxxx_0[i] = 4.0 * ta1_x_0_xxx_0[i] * fe_0 - 4.0 * ta1_x_0_xxx_1[i] * fe_0 + ta_0_xxxx_1[i] + ta1_x_0_xxxx_0[i] * pa_x[i] -
                            ta1_x_0_xxxx_1[i] * pc_x[i];

        ta1_x_x_xxxy_0[i] = 3.0 * ta1_x_0_xxy_0[i] * fe_0 - 3.0 * ta1_x_0_xxy_1[i] * fe_0 + ta_0_xxxy_1[i] + ta1_x_0_xxxy_0[i] * pa_x[i] -
                            ta1_x_0_xxxy_1[i] * pc_x[i];

        ta1_x_x_xxxz_0[i] = 3.0 * ta1_x_0_xxz_0[i] * fe_0 - 3.0 * ta1_x_0_xxz_1[i] * fe_0 + ta_0_xxxz_1[i] + ta1_x_0_xxxz_0[i] * pa_x[i] -
                            ta1_x_0_xxxz_1[i] * pc_x[i];

        ta1_x_x_xxyy_0[i] = 2.0 * ta1_x_0_xyy_0[i] * fe_0 - 2.0 * ta1_x_0_xyy_1[i] * fe_0 + ta_0_xxyy_1[i] + ta1_x_0_xxyy_0[i] * pa_x[i] -
                            ta1_x_0_xxyy_1[i] * pc_x[i];

        ta1_x_x_xxyz_0[i] = 2.0 * ta1_x_0_xyz_0[i] * fe_0 - 2.0 * ta1_x_0_xyz_1[i] * fe_0 + ta_0_xxyz_1[i] + ta1_x_0_xxyz_0[i] * pa_x[i] -
                            ta1_x_0_xxyz_1[i] * pc_x[i];

        ta1_x_x_xxzz_0[i] = 2.0 * ta1_x_0_xzz_0[i] * fe_0 - 2.0 * ta1_x_0_xzz_1[i] * fe_0 + ta_0_xxzz_1[i] + ta1_x_0_xxzz_0[i] * pa_x[i] -
                            ta1_x_0_xxzz_1[i] * pc_x[i];

        ta1_x_x_xyyy_0[i] =
            ta1_x_0_yyy_0[i] * fe_0 - ta1_x_0_yyy_1[i] * fe_0 + ta_0_xyyy_1[i] + ta1_x_0_xyyy_0[i] * pa_x[i] - ta1_x_0_xyyy_1[i] * pc_x[i];

        ta1_x_x_xyyz_0[i] =
            ta1_x_0_yyz_0[i] * fe_0 - ta1_x_0_yyz_1[i] * fe_0 + ta_0_xyyz_1[i] + ta1_x_0_xyyz_0[i] * pa_x[i] - ta1_x_0_xyyz_1[i] * pc_x[i];

        ta1_x_x_xyzz_0[i] =
            ta1_x_0_yzz_0[i] * fe_0 - ta1_x_0_yzz_1[i] * fe_0 + ta_0_xyzz_1[i] + ta1_x_0_xyzz_0[i] * pa_x[i] - ta1_x_0_xyzz_1[i] * pc_x[i];

        ta1_x_x_xzzz_0[i] =
            ta1_x_0_zzz_0[i] * fe_0 - ta1_x_0_zzz_1[i] * fe_0 + ta_0_xzzz_1[i] + ta1_x_0_xzzz_0[i] * pa_x[i] - ta1_x_0_xzzz_1[i] * pc_x[i];

        ta1_x_x_yyyy_0[i] = ta_0_yyyy_1[i] + ta1_x_0_yyyy_0[i] * pa_x[i] - ta1_x_0_yyyy_1[i] * pc_x[i];

        ta1_x_x_yyyz_0[i] = ta_0_yyyz_1[i] + ta1_x_0_yyyz_0[i] * pa_x[i] - ta1_x_0_yyyz_1[i] * pc_x[i];

        ta1_x_x_yyzz_0[i] = ta_0_yyzz_1[i] + ta1_x_0_yyzz_0[i] * pa_x[i] - ta1_x_0_yyzz_1[i] * pc_x[i];

        ta1_x_x_yzzz_0[i] = ta_0_yzzz_1[i] + ta1_x_0_yzzz_0[i] * pa_x[i] - ta1_x_0_yzzz_1[i] * pc_x[i];

        ta1_x_x_zzzz_0[i] = ta_0_zzzz_1[i] + ta1_x_0_zzzz_0[i] * pa_x[i] - ta1_x_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_0_xxx_0,  \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxxx_0, \
                             ta1_x_0_xxxx_1, \
                             ta1_x_0_xxxy_0, \
                             ta1_x_0_xxxy_1, \
                             ta1_x_0_xxxz_0, \
                             ta1_x_0_xxxz_1, \
                             ta1_x_0_xxy_0,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxyy_0, \
                             ta1_x_0_xxyy_1, \
                             ta1_x_0_xxyz_0, \
                             ta1_x_0_xxyz_1, \
                             ta1_x_0_xxz_0,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xxzz_0, \
                             ta1_x_0_xxzz_1, \
                             ta1_x_0_xyy_0,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyyy_0, \
                             ta1_x_0_xyyy_1, \
                             ta1_x_0_xyyz_0, \
                             ta1_x_0_xyyz_1, \
                             ta1_x_0_xyz_0,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xyzz_0, \
                             ta1_x_0_xyzz_1, \
                             ta1_x_0_xzz_0,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_xzzz_0, \
                             ta1_x_0_xzzz_1, \
                             ta1_x_0_yyy_0,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyyy_0, \
                             ta1_x_0_yyyy_1, \
                             ta1_x_0_yyyz_0, \
                             ta1_x_0_yyyz_1, \
                             ta1_x_0_yyz_0,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yyzz_0, \
                             ta1_x_0_yyzz_1, \
                             ta1_x_0_yzz_0,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_yzzz_0, \
                             ta1_x_0_yzzz_1, \
                             ta1_x_0_zzz_0,  \
                             ta1_x_0_zzz_1,  \
                             ta1_x_0_zzzz_0, \
                             ta1_x_0_zzzz_1, \
                             ta1_x_y_xxxx_0, \
                             ta1_x_y_xxxy_0, \
                             ta1_x_y_xxxz_0, \
                             ta1_x_y_xxyy_0, \
                             ta1_x_y_xxyz_0, \
                             ta1_x_y_xxzz_0, \
                             ta1_x_y_xyyy_0, \
                             ta1_x_y_xyyz_0, \
                             ta1_x_y_xyzz_0, \
                             ta1_x_y_xzzz_0, \
                             ta1_x_y_yyyy_0, \
                             ta1_x_y_yyyz_0, \
                             ta1_x_y_yyzz_0, \
                             ta1_x_y_yzzz_0, \
                             ta1_x_y_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_y_xxxx_0[i] = ta1_x_0_xxxx_0[i] * pa_y[i] - ta1_x_0_xxxx_1[i] * pc_y[i];

        ta1_x_y_xxxy_0[i] = ta1_x_0_xxx_0[i] * fe_0 - ta1_x_0_xxx_1[i] * fe_0 + ta1_x_0_xxxy_0[i] * pa_y[i] - ta1_x_0_xxxy_1[i] * pc_y[i];

        ta1_x_y_xxxz_0[i] = ta1_x_0_xxxz_0[i] * pa_y[i] - ta1_x_0_xxxz_1[i] * pc_y[i];

        ta1_x_y_xxyy_0[i] = 2.0 * ta1_x_0_xxy_0[i] * fe_0 - 2.0 * ta1_x_0_xxy_1[i] * fe_0 + ta1_x_0_xxyy_0[i] * pa_y[i] - ta1_x_0_xxyy_1[i] * pc_y[i];

        ta1_x_y_xxyz_0[i] = ta1_x_0_xxz_0[i] * fe_0 - ta1_x_0_xxz_1[i] * fe_0 + ta1_x_0_xxyz_0[i] * pa_y[i] - ta1_x_0_xxyz_1[i] * pc_y[i];

        ta1_x_y_xxzz_0[i] = ta1_x_0_xxzz_0[i] * pa_y[i] - ta1_x_0_xxzz_1[i] * pc_y[i];

        ta1_x_y_xyyy_0[i] = 3.0 * ta1_x_0_xyy_0[i] * fe_0 - 3.0 * ta1_x_0_xyy_1[i] * fe_0 + ta1_x_0_xyyy_0[i] * pa_y[i] - ta1_x_0_xyyy_1[i] * pc_y[i];

        ta1_x_y_xyyz_0[i] = 2.0 * ta1_x_0_xyz_0[i] * fe_0 - 2.0 * ta1_x_0_xyz_1[i] * fe_0 + ta1_x_0_xyyz_0[i] * pa_y[i] - ta1_x_0_xyyz_1[i] * pc_y[i];

        ta1_x_y_xyzz_0[i] = ta1_x_0_xzz_0[i] * fe_0 - ta1_x_0_xzz_1[i] * fe_0 + ta1_x_0_xyzz_0[i] * pa_y[i] - ta1_x_0_xyzz_1[i] * pc_y[i];

        ta1_x_y_xzzz_0[i] = ta1_x_0_xzzz_0[i] * pa_y[i] - ta1_x_0_xzzz_1[i] * pc_y[i];

        ta1_x_y_yyyy_0[i] = 4.0 * ta1_x_0_yyy_0[i] * fe_0 - 4.0 * ta1_x_0_yyy_1[i] * fe_0 + ta1_x_0_yyyy_0[i] * pa_y[i] - ta1_x_0_yyyy_1[i] * pc_y[i];

        ta1_x_y_yyyz_0[i] = 3.0 * ta1_x_0_yyz_0[i] * fe_0 - 3.0 * ta1_x_0_yyz_1[i] * fe_0 + ta1_x_0_yyyz_0[i] * pa_y[i] - ta1_x_0_yyyz_1[i] * pc_y[i];

        ta1_x_y_yyzz_0[i] = 2.0 * ta1_x_0_yzz_0[i] * fe_0 - 2.0 * ta1_x_0_yzz_1[i] * fe_0 + ta1_x_0_yyzz_0[i] * pa_y[i] - ta1_x_0_yyzz_1[i] * pc_y[i];

        ta1_x_y_yzzz_0[i] = ta1_x_0_zzz_0[i] * fe_0 - ta1_x_0_zzz_1[i] * fe_0 + ta1_x_0_yzzz_0[i] * pa_y[i] - ta1_x_0_yzzz_1[i] * pc_y[i];

        ta1_x_y_zzzz_0[i] = ta1_x_0_zzzz_0[i] * pa_y[i] - ta1_x_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_0_xxx_0,  \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxxx_0, \
                             ta1_x_0_xxxx_1, \
                             ta1_x_0_xxxy_0, \
                             ta1_x_0_xxxy_1, \
                             ta1_x_0_xxxz_0, \
                             ta1_x_0_xxxz_1, \
                             ta1_x_0_xxy_0,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxyy_0, \
                             ta1_x_0_xxyy_1, \
                             ta1_x_0_xxyz_0, \
                             ta1_x_0_xxyz_1, \
                             ta1_x_0_xxz_0,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xxzz_0, \
                             ta1_x_0_xxzz_1, \
                             ta1_x_0_xyy_0,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyyy_0, \
                             ta1_x_0_xyyy_1, \
                             ta1_x_0_xyyz_0, \
                             ta1_x_0_xyyz_1, \
                             ta1_x_0_xyz_0,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xyzz_0, \
                             ta1_x_0_xyzz_1, \
                             ta1_x_0_xzz_0,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_xzzz_0, \
                             ta1_x_0_xzzz_1, \
                             ta1_x_0_yyy_0,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyyy_0, \
                             ta1_x_0_yyyy_1, \
                             ta1_x_0_yyyz_0, \
                             ta1_x_0_yyyz_1, \
                             ta1_x_0_yyz_0,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yyzz_0, \
                             ta1_x_0_yyzz_1, \
                             ta1_x_0_yzz_0,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_yzzz_0, \
                             ta1_x_0_yzzz_1, \
                             ta1_x_0_zzz_0,  \
                             ta1_x_0_zzz_1,  \
                             ta1_x_0_zzzz_0, \
                             ta1_x_0_zzzz_1, \
                             ta1_x_z_xxxx_0, \
                             ta1_x_z_xxxy_0, \
                             ta1_x_z_xxxz_0, \
                             ta1_x_z_xxyy_0, \
                             ta1_x_z_xxyz_0, \
                             ta1_x_z_xxzz_0, \
                             ta1_x_z_xyyy_0, \
                             ta1_x_z_xyyz_0, \
                             ta1_x_z_xyzz_0, \
                             ta1_x_z_xzzz_0, \
                             ta1_x_z_yyyy_0, \
                             ta1_x_z_yyyz_0, \
                             ta1_x_z_yyzz_0, \
                             ta1_x_z_yzzz_0, \
                             ta1_x_z_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_z_xxxx_0[i] = ta1_x_0_xxxx_0[i] * pa_z[i] - ta1_x_0_xxxx_1[i] * pc_z[i];

        ta1_x_z_xxxy_0[i] = ta1_x_0_xxxy_0[i] * pa_z[i] - ta1_x_0_xxxy_1[i] * pc_z[i];

        ta1_x_z_xxxz_0[i] = ta1_x_0_xxx_0[i] * fe_0 - ta1_x_0_xxx_1[i] * fe_0 + ta1_x_0_xxxz_0[i] * pa_z[i] - ta1_x_0_xxxz_1[i] * pc_z[i];

        ta1_x_z_xxyy_0[i] = ta1_x_0_xxyy_0[i] * pa_z[i] - ta1_x_0_xxyy_1[i] * pc_z[i];

        ta1_x_z_xxyz_0[i] = ta1_x_0_xxy_0[i] * fe_0 - ta1_x_0_xxy_1[i] * fe_0 + ta1_x_0_xxyz_0[i] * pa_z[i] - ta1_x_0_xxyz_1[i] * pc_z[i];

        ta1_x_z_xxzz_0[i] = 2.0 * ta1_x_0_xxz_0[i] * fe_0 - 2.0 * ta1_x_0_xxz_1[i] * fe_0 + ta1_x_0_xxzz_0[i] * pa_z[i] - ta1_x_0_xxzz_1[i] * pc_z[i];

        ta1_x_z_xyyy_0[i] = ta1_x_0_xyyy_0[i] * pa_z[i] - ta1_x_0_xyyy_1[i] * pc_z[i];

        ta1_x_z_xyyz_0[i] = ta1_x_0_xyy_0[i] * fe_0 - ta1_x_0_xyy_1[i] * fe_0 + ta1_x_0_xyyz_0[i] * pa_z[i] - ta1_x_0_xyyz_1[i] * pc_z[i];

        ta1_x_z_xyzz_0[i] = 2.0 * ta1_x_0_xyz_0[i] * fe_0 - 2.0 * ta1_x_0_xyz_1[i] * fe_0 + ta1_x_0_xyzz_0[i] * pa_z[i] - ta1_x_0_xyzz_1[i] * pc_z[i];

        ta1_x_z_xzzz_0[i] = 3.0 * ta1_x_0_xzz_0[i] * fe_0 - 3.0 * ta1_x_0_xzz_1[i] * fe_0 + ta1_x_0_xzzz_0[i] * pa_z[i] - ta1_x_0_xzzz_1[i] * pc_z[i];

        ta1_x_z_yyyy_0[i] = ta1_x_0_yyyy_0[i] * pa_z[i] - ta1_x_0_yyyy_1[i] * pc_z[i];

        ta1_x_z_yyyz_0[i] = ta1_x_0_yyy_0[i] * fe_0 - ta1_x_0_yyy_1[i] * fe_0 + ta1_x_0_yyyz_0[i] * pa_z[i] - ta1_x_0_yyyz_1[i] * pc_z[i];

        ta1_x_z_yyzz_0[i] = 2.0 * ta1_x_0_yyz_0[i] * fe_0 - 2.0 * ta1_x_0_yyz_1[i] * fe_0 + ta1_x_0_yyzz_0[i] * pa_z[i] - ta1_x_0_yyzz_1[i] * pc_z[i];

        ta1_x_z_yzzz_0[i] = 3.0 * ta1_x_0_yzz_0[i] * fe_0 - 3.0 * ta1_x_0_yzz_1[i] * fe_0 + ta1_x_0_yzzz_0[i] * pa_z[i] - ta1_x_0_yzzz_1[i] * pc_z[i];

        ta1_x_z_zzzz_0[i] = 4.0 * ta1_x_0_zzz_0[i] * fe_0 - 4.0 * ta1_x_0_zzz_1[i] * fe_0 + ta1_x_0_zzzz_0[i] * pa_z[i] - ta1_x_0_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_0_xxx_0,  \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxxx_0, \
                             ta1_y_0_xxxx_1, \
                             ta1_y_0_xxxy_0, \
                             ta1_y_0_xxxy_1, \
                             ta1_y_0_xxxz_0, \
                             ta1_y_0_xxxz_1, \
                             ta1_y_0_xxy_0,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxyy_0, \
                             ta1_y_0_xxyy_1, \
                             ta1_y_0_xxyz_0, \
                             ta1_y_0_xxyz_1, \
                             ta1_y_0_xxz_0,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xxzz_0, \
                             ta1_y_0_xxzz_1, \
                             ta1_y_0_xyy_0,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyyy_0, \
                             ta1_y_0_xyyy_1, \
                             ta1_y_0_xyyz_0, \
                             ta1_y_0_xyyz_1, \
                             ta1_y_0_xyz_0,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xyzz_0, \
                             ta1_y_0_xyzz_1, \
                             ta1_y_0_xzz_0,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_xzzz_0, \
                             ta1_y_0_xzzz_1, \
                             ta1_y_0_yyy_0,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyyy_0, \
                             ta1_y_0_yyyy_1, \
                             ta1_y_0_yyyz_0, \
                             ta1_y_0_yyyz_1, \
                             ta1_y_0_yyz_0,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yyzz_0, \
                             ta1_y_0_yyzz_1, \
                             ta1_y_0_yzz_0,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_yzzz_0, \
                             ta1_y_0_yzzz_1, \
                             ta1_y_0_zzz_0,  \
                             ta1_y_0_zzz_1,  \
                             ta1_y_0_zzzz_0, \
                             ta1_y_0_zzzz_1, \
                             ta1_y_x_xxxx_0, \
                             ta1_y_x_xxxy_0, \
                             ta1_y_x_xxxz_0, \
                             ta1_y_x_xxyy_0, \
                             ta1_y_x_xxyz_0, \
                             ta1_y_x_xxzz_0, \
                             ta1_y_x_xyyy_0, \
                             ta1_y_x_xyyz_0, \
                             ta1_y_x_xyzz_0, \
                             ta1_y_x_xzzz_0, \
                             ta1_y_x_yyyy_0, \
                             ta1_y_x_yyyz_0, \
                             ta1_y_x_yyzz_0, \
                             ta1_y_x_yzzz_0, \
                             ta1_y_x_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_x_xxxx_0[i] = 4.0 * ta1_y_0_xxx_0[i] * fe_0 - 4.0 * ta1_y_0_xxx_1[i] * fe_0 + ta1_y_0_xxxx_0[i] * pa_x[i] - ta1_y_0_xxxx_1[i] * pc_x[i];

        ta1_y_x_xxxy_0[i] = 3.0 * ta1_y_0_xxy_0[i] * fe_0 - 3.0 * ta1_y_0_xxy_1[i] * fe_0 + ta1_y_0_xxxy_0[i] * pa_x[i] - ta1_y_0_xxxy_1[i] * pc_x[i];

        ta1_y_x_xxxz_0[i] = 3.0 * ta1_y_0_xxz_0[i] * fe_0 - 3.0 * ta1_y_0_xxz_1[i] * fe_0 + ta1_y_0_xxxz_0[i] * pa_x[i] - ta1_y_0_xxxz_1[i] * pc_x[i];

        ta1_y_x_xxyy_0[i] = 2.0 * ta1_y_0_xyy_0[i] * fe_0 - 2.0 * ta1_y_0_xyy_1[i] * fe_0 + ta1_y_0_xxyy_0[i] * pa_x[i] - ta1_y_0_xxyy_1[i] * pc_x[i];

        ta1_y_x_xxyz_0[i] = 2.0 * ta1_y_0_xyz_0[i] * fe_0 - 2.0 * ta1_y_0_xyz_1[i] * fe_0 + ta1_y_0_xxyz_0[i] * pa_x[i] - ta1_y_0_xxyz_1[i] * pc_x[i];

        ta1_y_x_xxzz_0[i] = 2.0 * ta1_y_0_xzz_0[i] * fe_0 - 2.0 * ta1_y_0_xzz_1[i] * fe_0 + ta1_y_0_xxzz_0[i] * pa_x[i] - ta1_y_0_xxzz_1[i] * pc_x[i];

        ta1_y_x_xyyy_0[i] = ta1_y_0_yyy_0[i] * fe_0 - ta1_y_0_yyy_1[i] * fe_0 + ta1_y_0_xyyy_0[i] * pa_x[i] - ta1_y_0_xyyy_1[i] * pc_x[i];

        ta1_y_x_xyyz_0[i] = ta1_y_0_yyz_0[i] * fe_0 - ta1_y_0_yyz_1[i] * fe_0 + ta1_y_0_xyyz_0[i] * pa_x[i] - ta1_y_0_xyyz_1[i] * pc_x[i];

        ta1_y_x_xyzz_0[i] = ta1_y_0_yzz_0[i] * fe_0 - ta1_y_0_yzz_1[i] * fe_0 + ta1_y_0_xyzz_0[i] * pa_x[i] - ta1_y_0_xyzz_1[i] * pc_x[i];

        ta1_y_x_xzzz_0[i] = ta1_y_0_zzz_0[i] * fe_0 - ta1_y_0_zzz_1[i] * fe_0 + ta1_y_0_xzzz_0[i] * pa_x[i] - ta1_y_0_xzzz_1[i] * pc_x[i];

        ta1_y_x_yyyy_0[i] = ta1_y_0_yyyy_0[i] * pa_x[i] - ta1_y_0_yyyy_1[i] * pc_x[i];

        ta1_y_x_yyyz_0[i] = ta1_y_0_yyyz_0[i] * pa_x[i] - ta1_y_0_yyyz_1[i] * pc_x[i];

        ta1_y_x_yyzz_0[i] = ta1_y_0_yyzz_0[i] * pa_x[i] - ta1_y_0_yyzz_1[i] * pc_x[i];

        ta1_y_x_yzzz_0[i] = ta1_y_0_yzzz_0[i] * pa_x[i] - ta1_y_0_yzzz_1[i] * pc_x[i];

        ta1_y_x_zzzz_0[i] = ta1_y_0_zzzz_0[i] * pa_x[i] - ta1_y_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_y_0_xxx_0,  \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxxx_0, \
                             ta1_y_0_xxxx_1, \
                             ta1_y_0_xxxy_0, \
                             ta1_y_0_xxxy_1, \
                             ta1_y_0_xxxz_0, \
                             ta1_y_0_xxxz_1, \
                             ta1_y_0_xxy_0,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxyy_0, \
                             ta1_y_0_xxyy_1, \
                             ta1_y_0_xxyz_0, \
                             ta1_y_0_xxyz_1, \
                             ta1_y_0_xxz_0,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xxzz_0, \
                             ta1_y_0_xxzz_1, \
                             ta1_y_0_xyy_0,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyyy_0, \
                             ta1_y_0_xyyy_1, \
                             ta1_y_0_xyyz_0, \
                             ta1_y_0_xyyz_1, \
                             ta1_y_0_xyz_0,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xyzz_0, \
                             ta1_y_0_xyzz_1, \
                             ta1_y_0_xzz_0,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_xzzz_0, \
                             ta1_y_0_xzzz_1, \
                             ta1_y_0_yyy_0,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyyy_0, \
                             ta1_y_0_yyyy_1, \
                             ta1_y_0_yyyz_0, \
                             ta1_y_0_yyyz_1, \
                             ta1_y_0_yyz_0,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yyzz_0, \
                             ta1_y_0_yyzz_1, \
                             ta1_y_0_yzz_0,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_yzzz_0, \
                             ta1_y_0_yzzz_1, \
                             ta1_y_0_zzz_0,  \
                             ta1_y_0_zzz_1,  \
                             ta1_y_0_zzzz_0, \
                             ta1_y_0_zzzz_1, \
                             ta1_y_y_xxxx_0, \
                             ta1_y_y_xxxy_0, \
                             ta1_y_y_xxxz_0, \
                             ta1_y_y_xxyy_0, \
                             ta1_y_y_xxyz_0, \
                             ta1_y_y_xxzz_0, \
                             ta1_y_y_xyyy_0, \
                             ta1_y_y_xyyz_0, \
                             ta1_y_y_xyzz_0, \
                             ta1_y_y_xzzz_0, \
                             ta1_y_y_yyyy_0, \
                             ta1_y_y_yyyz_0, \
                             ta1_y_y_yyzz_0, \
                             ta1_y_y_yzzz_0, \
                             ta1_y_y_zzzz_0, \
                             ta_0_xxxx_1,    \
                             ta_0_xxxy_1,    \
                             ta_0_xxxz_1,    \
                             ta_0_xxyy_1,    \
                             ta_0_xxyz_1,    \
                             ta_0_xxzz_1,    \
                             ta_0_xyyy_1,    \
                             ta_0_xyyz_1,    \
                             ta_0_xyzz_1,    \
                             ta_0_xzzz_1,    \
                             ta_0_yyyy_1,    \
                             ta_0_yyyz_1,    \
                             ta_0_yyzz_1,    \
                             ta_0_yzzz_1,    \
                             ta_0_zzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_y_xxxx_0[i] = ta_0_xxxx_1[i] + ta1_y_0_xxxx_0[i] * pa_y[i] - ta1_y_0_xxxx_1[i] * pc_y[i];

        ta1_y_y_xxxy_0[i] =
            ta1_y_0_xxx_0[i] * fe_0 - ta1_y_0_xxx_1[i] * fe_0 + ta_0_xxxy_1[i] + ta1_y_0_xxxy_0[i] * pa_y[i] - ta1_y_0_xxxy_1[i] * pc_y[i];

        ta1_y_y_xxxz_0[i] = ta_0_xxxz_1[i] + ta1_y_0_xxxz_0[i] * pa_y[i] - ta1_y_0_xxxz_1[i] * pc_y[i];

        ta1_y_y_xxyy_0[i] = 2.0 * ta1_y_0_xxy_0[i] * fe_0 - 2.0 * ta1_y_0_xxy_1[i] * fe_0 + ta_0_xxyy_1[i] + ta1_y_0_xxyy_0[i] * pa_y[i] -
                            ta1_y_0_xxyy_1[i] * pc_y[i];

        ta1_y_y_xxyz_0[i] =
            ta1_y_0_xxz_0[i] * fe_0 - ta1_y_0_xxz_1[i] * fe_0 + ta_0_xxyz_1[i] + ta1_y_0_xxyz_0[i] * pa_y[i] - ta1_y_0_xxyz_1[i] * pc_y[i];

        ta1_y_y_xxzz_0[i] = ta_0_xxzz_1[i] + ta1_y_0_xxzz_0[i] * pa_y[i] - ta1_y_0_xxzz_1[i] * pc_y[i];

        ta1_y_y_xyyy_0[i] = 3.0 * ta1_y_0_xyy_0[i] * fe_0 - 3.0 * ta1_y_0_xyy_1[i] * fe_0 + ta_0_xyyy_1[i] + ta1_y_0_xyyy_0[i] * pa_y[i] -
                            ta1_y_0_xyyy_1[i] * pc_y[i];

        ta1_y_y_xyyz_0[i] = 2.0 * ta1_y_0_xyz_0[i] * fe_0 - 2.0 * ta1_y_0_xyz_1[i] * fe_0 + ta_0_xyyz_1[i] + ta1_y_0_xyyz_0[i] * pa_y[i] -
                            ta1_y_0_xyyz_1[i] * pc_y[i];

        ta1_y_y_xyzz_0[i] =
            ta1_y_0_xzz_0[i] * fe_0 - ta1_y_0_xzz_1[i] * fe_0 + ta_0_xyzz_1[i] + ta1_y_0_xyzz_0[i] * pa_y[i] - ta1_y_0_xyzz_1[i] * pc_y[i];

        ta1_y_y_xzzz_0[i] = ta_0_xzzz_1[i] + ta1_y_0_xzzz_0[i] * pa_y[i] - ta1_y_0_xzzz_1[i] * pc_y[i];

        ta1_y_y_yyyy_0[i] = 4.0 * ta1_y_0_yyy_0[i] * fe_0 - 4.0 * ta1_y_0_yyy_1[i] * fe_0 + ta_0_yyyy_1[i] + ta1_y_0_yyyy_0[i] * pa_y[i] -
                            ta1_y_0_yyyy_1[i] * pc_y[i];

        ta1_y_y_yyyz_0[i] = 3.0 * ta1_y_0_yyz_0[i] * fe_0 - 3.0 * ta1_y_0_yyz_1[i] * fe_0 + ta_0_yyyz_1[i] + ta1_y_0_yyyz_0[i] * pa_y[i] -
                            ta1_y_0_yyyz_1[i] * pc_y[i];

        ta1_y_y_yyzz_0[i] = 2.0 * ta1_y_0_yzz_0[i] * fe_0 - 2.0 * ta1_y_0_yzz_1[i] * fe_0 + ta_0_yyzz_1[i] + ta1_y_0_yyzz_0[i] * pa_y[i] -
                            ta1_y_0_yyzz_1[i] * pc_y[i];

        ta1_y_y_yzzz_0[i] =
            ta1_y_0_zzz_0[i] * fe_0 - ta1_y_0_zzz_1[i] * fe_0 + ta_0_yzzz_1[i] + ta1_y_0_yzzz_0[i] * pa_y[i] - ta1_y_0_yzzz_1[i] * pc_y[i];

        ta1_y_y_zzzz_0[i] = ta_0_zzzz_1[i] + ta1_y_0_zzzz_0[i] * pa_y[i] - ta1_y_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_0_xxx_0,  \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxxx_0, \
                             ta1_y_0_xxxx_1, \
                             ta1_y_0_xxxy_0, \
                             ta1_y_0_xxxy_1, \
                             ta1_y_0_xxxz_0, \
                             ta1_y_0_xxxz_1, \
                             ta1_y_0_xxy_0,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxyy_0, \
                             ta1_y_0_xxyy_1, \
                             ta1_y_0_xxyz_0, \
                             ta1_y_0_xxyz_1, \
                             ta1_y_0_xxz_0,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xxzz_0, \
                             ta1_y_0_xxzz_1, \
                             ta1_y_0_xyy_0,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyyy_0, \
                             ta1_y_0_xyyy_1, \
                             ta1_y_0_xyyz_0, \
                             ta1_y_0_xyyz_1, \
                             ta1_y_0_xyz_0,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xyzz_0, \
                             ta1_y_0_xyzz_1, \
                             ta1_y_0_xzz_0,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_xzzz_0, \
                             ta1_y_0_xzzz_1, \
                             ta1_y_0_yyy_0,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyyy_0, \
                             ta1_y_0_yyyy_1, \
                             ta1_y_0_yyyz_0, \
                             ta1_y_0_yyyz_1, \
                             ta1_y_0_yyz_0,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yyzz_0, \
                             ta1_y_0_yyzz_1, \
                             ta1_y_0_yzz_0,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_yzzz_0, \
                             ta1_y_0_yzzz_1, \
                             ta1_y_0_zzz_0,  \
                             ta1_y_0_zzz_1,  \
                             ta1_y_0_zzzz_0, \
                             ta1_y_0_zzzz_1, \
                             ta1_y_z_xxxx_0, \
                             ta1_y_z_xxxy_0, \
                             ta1_y_z_xxxz_0, \
                             ta1_y_z_xxyy_0, \
                             ta1_y_z_xxyz_0, \
                             ta1_y_z_xxzz_0, \
                             ta1_y_z_xyyy_0, \
                             ta1_y_z_xyyz_0, \
                             ta1_y_z_xyzz_0, \
                             ta1_y_z_xzzz_0, \
                             ta1_y_z_yyyy_0, \
                             ta1_y_z_yyyz_0, \
                             ta1_y_z_yyzz_0, \
                             ta1_y_z_yzzz_0, \
                             ta1_y_z_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_z_xxxx_0[i] = ta1_y_0_xxxx_0[i] * pa_z[i] - ta1_y_0_xxxx_1[i] * pc_z[i];

        ta1_y_z_xxxy_0[i] = ta1_y_0_xxxy_0[i] * pa_z[i] - ta1_y_0_xxxy_1[i] * pc_z[i];

        ta1_y_z_xxxz_0[i] = ta1_y_0_xxx_0[i] * fe_0 - ta1_y_0_xxx_1[i] * fe_0 + ta1_y_0_xxxz_0[i] * pa_z[i] - ta1_y_0_xxxz_1[i] * pc_z[i];

        ta1_y_z_xxyy_0[i] = ta1_y_0_xxyy_0[i] * pa_z[i] - ta1_y_0_xxyy_1[i] * pc_z[i];

        ta1_y_z_xxyz_0[i] = ta1_y_0_xxy_0[i] * fe_0 - ta1_y_0_xxy_1[i] * fe_0 + ta1_y_0_xxyz_0[i] * pa_z[i] - ta1_y_0_xxyz_1[i] * pc_z[i];

        ta1_y_z_xxzz_0[i] = 2.0 * ta1_y_0_xxz_0[i] * fe_0 - 2.0 * ta1_y_0_xxz_1[i] * fe_0 + ta1_y_0_xxzz_0[i] * pa_z[i] - ta1_y_0_xxzz_1[i] * pc_z[i];

        ta1_y_z_xyyy_0[i] = ta1_y_0_xyyy_0[i] * pa_z[i] - ta1_y_0_xyyy_1[i] * pc_z[i];

        ta1_y_z_xyyz_0[i] = ta1_y_0_xyy_0[i] * fe_0 - ta1_y_0_xyy_1[i] * fe_0 + ta1_y_0_xyyz_0[i] * pa_z[i] - ta1_y_0_xyyz_1[i] * pc_z[i];

        ta1_y_z_xyzz_0[i] = 2.0 * ta1_y_0_xyz_0[i] * fe_0 - 2.0 * ta1_y_0_xyz_1[i] * fe_0 + ta1_y_0_xyzz_0[i] * pa_z[i] - ta1_y_0_xyzz_1[i] * pc_z[i];

        ta1_y_z_xzzz_0[i] = 3.0 * ta1_y_0_xzz_0[i] * fe_0 - 3.0 * ta1_y_0_xzz_1[i] * fe_0 + ta1_y_0_xzzz_0[i] * pa_z[i] - ta1_y_0_xzzz_1[i] * pc_z[i];

        ta1_y_z_yyyy_0[i] = ta1_y_0_yyyy_0[i] * pa_z[i] - ta1_y_0_yyyy_1[i] * pc_z[i];

        ta1_y_z_yyyz_0[i] = ta1_y_0_yyy_0[i] * fe_0 - ta1_y_0_yyy_1[i] * fe_0 + ta1_y_0_yyyz_0[i] * pa_z[i] - ta1_y_0_yyyz_1[i] * pc_z[i];

        ta1_y_z_yyzz_0[i] = 2.0 * ta1_y_0_yyz_0[i] * fe_0 - 2.0 * ta1_y_0_yyz_1[i] * fe_0 + ta1_y_0_yyzz_0[i] * pa_z[i] - ta1_y_0_yyzz_1[i] * pc_z[i];

        ta1_y_z_yzzz_0[i] = 3.0 * ta1_y_0_yzz_0[i] * fe_0 - 3.0 * ta1_y_0_yzz_1[i] * fe_0 + ta1_y_0_yzzz_0[i] * pa_z[i] - ta1_y_0_yzzz_1[i] * pc_z[i];

        ta1_y_z_zzzz_0[i] = 4.0 * ta1_y_0_zzz_0[i] * fe_0 - 4.0 * ta1_y_0_zzz_1[i] * fe_0 + ta1_y_0_zzzz_0[i] * pa_z[i] - ta1_y_0_zzzz_1[i] * pc_z[i];
    }

    // Set up 90-105 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_0_xxx_0,  \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxxx_0, \
                             ta1_z_0_xxxx_1, \
                             ta1_z_0_xxxy_0, \
                             ta1_z_0_xxxy_1, \
                             ta1_z_0_xxxz_0, \
                             ta1_z_0_xxxz_1, \
                             ta1_z_0_xxy_0,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxyy_0, \
                             ta1_z_0_xxyy_1, \
                             ta1_z_0_xxyz_0, \
                             ta1_z_0_xxyz_1, \
                             ta1_z_0_xxz_0,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xxzz_0, \
                             ta1_z_0_xxzz_1, \
                             ta1_z_0_xyy_0,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyyy_0, \
                             ta1_z_0_xyyy_1, \
                             ta1_z_0_xyyz_0, \
                             ta1_z_0_xyyz_1, \
                             ta1_z_0_xyz_0,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xyzz_0, \
                             ta1_z_0_xyzz_1, \
                             ta1_z_0_xzz_0,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_xzzz_0, \
                             ta1_z_0_xzzz_1, \
                             ta1_z_0_yyy_0,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyyy_0, \
                             ta1_z_0_yyyy_1, \
                             ta1_z_0_yyyz_0, \
                             ta1_z_0_yyyz_1, \
                             ta1_z_0_yyz_0,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yyzz_0, \
                             ta1_z_0_yyzz_1, \
                             ta1_z_0_yzz_0,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_yzzz_0, \
                             ta1_z_0_yzzz_1, \
                             ta1_z_0_zzz_0,  \
                             ta1_z_0_zzz_1,  \
                             ta1_z_0_zzzz_0, \
                             ta1_z_0_zzzz_1, \
                             ta1_z_x_xxxx_0, \
                             ta1_z_x_xxxy_0, \
                             ta1_z_x_xxxz_0, \
                             ta1_z_x_xxyy_0, \
                             ta1_z_x_xxyz_0, \
                             ta1_z_x_xxzz_0, \
                             ta1_z_x_xyyy_0, \
                             ta1_z_x_xyyz_0, \
                             ta1_z_x_xyzz_0, \
                             ta1_z_x_xzzz_0, \
                             ta1_z_x_yyyy_0, \
                             ta1_z_x_yyyz_0, \
                             ta1_z_x_yyzz_0, \
                             ta1_z_x_yzzz_0, \
                             ta1_z_x_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_x_xxxx_0[i] = 4.0 * ta1_z_0_xxx_0[i] * fe_0 - 4.0 * ta1_z_0_xxx_1[i] * fe_0 + ta1_z_0_xxxx_0[i] * pa_x[i] - ta1_z_0_xxxx_1[i] * pc_x[i];

        ta1_z_x_xxxy_0[i] = 3.0 * ta1_z_0_xxy_0[i] * fe_0 - 3.0 * ta1_z_0_xxy_1[i] * fe_0 + ta1_z_0_xxxy_0[i] * pa_x[i] - ta1_z_0_xxxy_1[i] * pc_x[i];

        ta1_z_x_xxxz_0[i] = 3.0 * ta1_z_0_xxz_0[i] * fe_0 - 3.0 * ta1_z_0_xxz_1[i] * fe_0 + ta1_z_0_xxxz_0[i] * pa_x[i] - ta1_z_0_xxxz_1[i] * pc_x[i];

        ta1_z_x_xxyy_0[i] = 2.0 * ta1_z_0_xyy_0[i] * fe_0 - 2.0 * ta1_z_0_xyy_1[i] * fe_0 + ta1_z_0_xxyy_0[i] * pa_x[i] - ta1_z_0_xxyy_1[i] * pc_x[i];

        ta1_z_x_xxyz_0[i] = 2.0 * ta1_z_0_xyz_0[i] * fe_0 - 2.0 * ta1_z_0_xyz_1[i] * fe_0 + ta1_z_0_xxyz_0[i] * pa_x[i] - ta1_z_0_xxyz_1[i] * pc_x[i];

        ta1_z_x_xxzz_0[i] = 2.0 * ta1_z_0_xzz_0[i] * fe_0 - 2.0 * ta1_z_0_xzz_1[i] * fe_0 + ta1_z_0_xxzz_0[i] * pa_x[i] - ta1_z_0_xxzz_1[i] * pc_x[i];

        ta1_z_x_xyyy_0[i] = ta1_z_0_yyy_0[i] * fe_0 - ta1_z_0_yyy_1[i] * fe_0 + ta1_z_0_xyyy_0[i] * pa_x[i] - ta1_z_0_xyyy_1[i] * pc_x[i];

        ta1_z_x_xyyz_0[i] = ta1_z_0_yyz_0[i] * fe_0 - ta1_z_0_yyz_1[i] * fe_0 + ta1_z_0_xyyz_0[i] * pa_x[i] - ta1_z_0_xyyz_1[i] * pc_x[i];

        ta1_z_x_xyzz_0[i] = ta1_z_0_yzz_0[i] * fe_0 - ta1_z_0_yzz_1[i] * fe_0 + ta1_z_0_xyzz_0[i] * pa_x[i] - ta1_z_0_xyzz_1[i] * pc_x[i];

        ta1_z_x_xzzz_0[i] = ta1_z_0_zzz_0[i] * fe_0 - ta1_z_0_zzz_1[i] * fe_0 + ta1_z_0_xzzz_0[i] * pa_x[i] - ta1_z_0_xzzz_1[i] * pc_x[i];

        ta1_z_x_yyyy_0[i] = ta1_z_0_yyyy_0[i] * pa_x[i] - ta1_z_0_yyyy_1[i] * pc_x[i];

        ta1_z_x_yyyz_0[i] = ta1_z_0_yyyz_0[i] * pa_x[i] - ta1_z_0_yyyz_1[i] * pc_x[i];

        ta1_z_x_yyzz_0[i] = ta1_z_0_yyzz_0[i] * pa_x[i] - ta1_z_0_yyzz_1[i] * pc_x[i];

        ta1_z_x_yzzz_0[i] = ta1_z_0_yzzz_0[i] * pa_x[i] - ta1_z_0_yzzz_1[i] * pc_x[i];

        ta1_z_x_zzzz_0[i] = ta1_z_0_zzzz_0[i] * pa_x[i] - ta1_z_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_0_xxx_0,  \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxxx_0, \
                             ta1_z_0_xxxx_1, \
                             ta1_z_0_xxxy_0, \
                             ta1_z_0_xxxy_1, \
                             ta1_z_0_xxxz_0, \
                             ta1_z_0_xxxz_1, \
                             ta1_z_0_xxy_0,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxyy_0, \
                             ta1_z_0_xxyy_1, \
                             ta1_z_0_xxyz_0, \
                             ta1_z_0_xxyz_1, \
                             ta1_z_0_xxz_0,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xxzz_0, \
                             ta1_z_0_xxzz_1, \
                             ta1_z_0_xyy_0,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyyy_0, \
                             ta1_z_0_xyyy_1, \
                             ta1_z_0_xyyz_0, \
                             ta1_z_0_xyyz_1, \
                             ta1_z_0_xyz_0,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xyzz_0, \
                             ta1_z_0_xyzz_1, \
                             ta1_z_0_xzz_0,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_xzzz_0, \
                             ta1_z_0_xzzz_1, \
                             ta1_z_0_yyy_0,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyyy_0, \
                             ta1_z_0_yyyy_1, \
                             ta1_z_0_yyyz_0, \
                             ta1_z_0_yyyz_1, \
                             ta1_z_0_yyz_0,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yyzz_0, \
                             ta1_z_0_yyzz_1, \
                             ta1_z_0_yzz_0,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_yzzz_0, \
                             ta1_z_0_yzzz_1, \
                             ta1_z_0_zzz_0,  \
                             ta1_z_0_zzz_1,  \
                             ta1_z_0_zzzz_0, \
                             ta1_z_0_zzzz_1, \
                             ta1_z_y_xxxx_0, \
                             ta1_z_y_xxxy_0, \
                             ta1_z_y_xxxz_0, \
                             ta1_z_y_xxyy_0, \
                             ta1_z_y_xxyz_0, \
                             ta1_z_y_xxzz_0, \
                             ta1_z_y_xyyy_0, \
                             ta1_z_y_xyyz_0, \
                             ta1_z_y_xyzz_0, \
                             ta1_z_y_xzzz_0, \
                             ta1_z_y_yyyy_0, \
                             ta1_z_y_yyyz_0, \
                             ta1_z_y_yyzz_0, \
                             ta1_z_y_yzzz_0, \
                             ta1_z_y_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_y_xxxx_0[i] = ta1_z_0_xxxx_0[i] * pa_y[i] - ta1_z_0_xxxx_1[i] * pc_y[i];

        ta1_z_y_xxxy_0[i] = ta1_z_0_xxx_0[i] * fe_0 - ta1_z_0_xxx_1[i] * fe_0 + ta1_z_0_xxxy_0[i] * pa_y[i] - ta1_z_0_xxxy_1[i] * pc_y[i];

        ta1_z_y_xxxz_0[i] = ta1_z_0_xxxz_0[i] * pa_y[i] - ta1_z_0_xxxz_1[i] * pc_y[i];

        ta1_z_y_xxyy_0[i] = 2.0 * ta1_z_0_xxy_0[i] * fe_0 - 2.0 * ta1_z_0_xxy_1[i] * fe_0 + ta1_z_0_xxyy_0[i] * pa_y[i] - ta1_z_0_xxyy_1[i] * pc_y[i];

        ta1_z_y_xxyz_0[i] = ta1_z_0_xxz_0[i] * fe_0 - ta1_z_0_xxz_1[i] * fe_0 + ta1_z_0_xxyz_0[i] * pa_y[i] - ta1_z_0_xxyz_1[i] * pc_y[i];

        ta1_z_y_xxzz_0[i] = ta1_z_0_xxzz_0[i] * pa_y[i] - ta1_z_0_xxzz_1[i] * pc_y[i];

        ta1_z_y_xyyy_0[i] = 3.0 * ta1_z_0_xyy_0[i] * fe_0 - 3.0 * ta1_z_0_xyy_1[i] * fe_0 + ta1_z_0_xyyy_0[i] * pa_y[i] - ta1_z_0_xyyy_1[i] * pc_y[i];

        ta1_z_y_xyyz_0[i] = 2.0 * ta1_z_0_xyz_0[i] * fe_0 - 2.0 * ta1_z_0_xyz_1[i] * fe_0 + ta1_z_0_xyyz_0[i] * pa_y[i] - ta1_z_0_xyyz_1[i] * pc_y[i];

        ta1_z_y_xyzz_0[i] = ta1_z_0_xzz_0[i] * fe_0 - ta1_z_0_xzz_1[i] * fe_0 + ta1_z_0_xyzz_0[i] * pa_y[i] - ta1_z_0_xyzz_1[i] * pc_y[i];

        ta1_z_y_xzzz_0[i] = ta1_z_0_xzzz_0[i] * pa_y[i] - ta1_z_0_xzzz_1[i] * pc_y[i];

        ta1_z_y_yyyy_0[i] = 4.0 * ta1_z_0_yyy_0[i] * fe_0 - 4.0 * ta1_z_0_yyy_1[i] * fe_0 + ta1_z_0_yyyy_0[i] * pa_y[i] - ta1_z_0_yyyy_1[i] * pc_y[i];

        ta1_z_y_yyyz_0[i] = 3.0 * ta1_z_0_yyz_0[i] * fe_0 - 3.0 * ta1_z_0_yyz_1[i] * fe_0 + ta1_z_0_yyyz_0[i] * pa_y[i] - ta1_z_0_yyyz_1[i] * pc_y[i];

        ta1_z_y_yyzz_0[i] = 2.0 * ta1_z_0_yzz_0[i] * fe_0 - 2.0 * ta1_z_0_yzz_1[i] * fe_0 + ta1_z_0_yyzz_0[i] * pa_y[i] - ta1_z_0_yyzz_1[i] * pc_y[i];

        ta1_z_y_yzzz_0[i] = ta1_z_0_zzz_0[i] * fe_0 - ta1_z_0_zzz_1[i] * fe_0 + ta1_z_0_yzzz_0[i] * pa_y[i] - ta1_z_0_yzzz_1[i] * pc_y[i];

        ta1_z_y_zzzz_0[i] = ta1_z_0_zzzz_0[i] * pa_y[i] - ta1_z_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 120-135 components of targeted buffer : PG

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

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_z_0_xxx_0,  \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxxx_0, \
                             ta1_z_0_xxxx_1, \
                             ta1_z_0_xxxy_0, \
                             ta1_z_0_xxxy_1, \
                             ta1_z_0_xxxz_0, \
                             ta1_z_0_xxxz_1, \
                             ta1_z_0_xxy_0,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxyy_0, \
                             ta1_z_0_xxyy_1, \
                             ta1_z_0_xxyz_0, \
                             ta1_z_0_xxyz_1, \
                             ta1_z_0_xxz_0,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xxzz_0, \
                             ta1_z_0_xxzz_1, \
                             ta1_z_0_xyy_0,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyyy_0, \
                             ta1_z_0_xyyy_1, \
                             ta1_z_0_xyyz_0, \
                             ta1_z_0_xyyz_1, \
                             ta1_z_0_xyz_0,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xyzz_0, \
                             ta1_z_0_xyzz_1, \
                             ta1_z_0_xzz_0,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_xzzz_0, \
                             ta1_z_0_xzzz_1, \
                             ta1_z_0_yyy_0,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyyy_0, \
                             ta1_z_0_yyyy_1, \
                             ta1_z_0_yyyz_0, \
                             ta1_z_0_yyyz_1, \
                             ta1_z_0_yyz_0,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yyzz_0, \
                             ta1_z_0_yyzz_1, \
                             ta1_z_0_yzz_0,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_yzzz_0, \
                             ta1_z_0_yzzz_1, \
                             ta1_z_0_zzz_0,  \
                             ta1_z_0_zzz_1,  \
                             ta1_z_0_zzzz_0, \
                             ta1_z_0_zzzz_1, \
                             ta1_z_z_xxxx_0, \
                             ta1_z_z_xxxy_0, \
                             ta1_z_z_xxxz_0, \
                             ta1_z_z_xxyy_0, \
                             ta1_z_z_xxyz_0, \
                             ta1_z_z_xxzz_0, \
                             ta1_z_z_xyyy_0, \
                             ta1_z_z_xyyz_0, \
                             ta1_z_z_xyzz_0, \
                             ta1_z_z_xzzz_0, \
                             ta1_z_z_yyyy_0, \
                             ta1_z_z_yyyz_0, \
                             ta1_z_z_yyzz_0, \
                             ta1_z_z_yzzz_0, \
                             ta1_z_z_zzzz_0, \
                             ta_0_xxxx_1,    \
                             ta_0_xxxy_1,    \
                             ta_0_xxxz_1,    \
                             ta_0_xxyy_1,    \
                             ta_0_xxyz_1,    \
                             ta_0_xxzz_1,    \
                             ta_0_xyyy_1,    \
                             ta_0_xyyz_1,    \
                             ta_0_xyzz_1,    \
                             ta_0_xzzz_1,    \
                             ta_0_yyyy_1,    \
                             ta_0_yyyz_1,    \
                             ta_0_yyzz_1,    \
                             ta_0_yzzz_1,    \
                             ta_0_zzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_z_xxxx_0[i] = ta_0_xxxx_1[i] + ta1_z_0_xxxx_0[i] * pa_z[i] - ta1_z_0_xxxx_1[i] * pc_z[i];

        ta1_z_z_xxxy_0[i] = ta_0_xxxy_1[i] + ta1_z_0_xxxy_0[i] * pa_z[i] - ta1_z_0_xxxy_1[i] * pc_z[i];

        ta1_z_z_xxxz_0[i] =
            ta1_z_0_xxx_0[i] * fe_0 - ta1_z_0_xxx_1[i] * fe_0 + ta_0_xxxz_1[i] + ta1_z_0_xxxz_0[i] * pa_z[i] - ta1_z_0_xxxz_1[i] * pc_z[i];

        ta1_z_z_xxyy_0[i] = ta_0_xxyy_1[i] + ta1_z_0_xxyy_0[i] * pa_z[i] - ta1_z_0_xxyy_1[i] * pc_z[i];

        ta1_z_z_xxyz_0[i] =
            ta1_z_0_xxy_0[i] * fe_0 - ta1_z_0_xxy_1[i] * fe_0 + ta_0_xxyz_1[i] + ta1_z_0_xxyz_0[i] * pa_z[i] - ta1_z_0_xxyz_1[i] * pc_z[i];

        ta1_z_z_xxzz_0[i] = 2.0 * ta1_z_0_xxz_0[i] * fe_0 - 2.0 * ta1_z_0_xxz_1[i] * fe_0 + ta_0_xxzz_1[i] + ta1_z_0_xxzz_0[i] * pa_z[i] -
                            ta1_z_0_xxzz_1[i] * pc_z[i];

        ta1_z_z_xyyy_0[i] = ta_0_xyyy_1[i] + ta1_z_0_xyyy_0[i] * pa_z[i] - ta1_z_0_xyyy_1[i] * pc_z[i];

        ta1_z_z_xyyz_0[i] =
            ta1_z_0_xyy_0[i] * fe_0 - ta1_z_0_xyy_1[i] * fe_0 + ta_0_xyyz_1[i] + ta1_z_0_xyyz_0[i] * pa_z[i] - ta1_z_0_xyyz_1[i] * pc_z[i];

        ta1_z_z_xyzz_0[i] = 2.0 * ta1_z_0_xyz_0[i] * fe_0 - 2.0 * ta1_z_0_xyz_1[i] * fe_0 + ta_0_xyzz_1[i] + ta1_z_0_xyzz_0[i] * pa_z[i] -
                            ta1_z_0_xyzz_1[i] * pc_z[i];

        ta1_z_z_xzzz_0[i] = 3.0 * ta1_z_0_xzz_0[i] * fe_0 - 3.0 * ta1_z_0_xzz_1[i] * fe_0 + ta_0_xzzz_1[i] + ta1_z_0_xzzz_0[i] * pa_z[i] -
                            ta1_z_0_xzzz_1[i] * pc_z[i];

        ta1_z_z_yyyy_0[i] = ta_0_yyyy_1[i] + ta1_z_0_yyyy_0[i] * pa_z[i] - ta1_z_0_yyyy_1[i] * pc_z[i];

        ta1_z_z_yyyz_0[i] =
            ta1_z_0_yyy_0[i] * fe_0 - ta1_z_0_yyy_1[i] * fe_0 + ta_0_yyyz_1[i] + ta1_z_0_yyyz_0[i] * pa_z[i] - ta1_z_0_yyyz_1[i] * pc_z[i];

        ta1_z_z_yyzz_0[i] = 2.0 * ta1_z_0_yyz_0[i] * fe_0 - 2.0 * ta1_z_0_yyz_1[i] * fe_0 + ta_0_yyzz_1[i] + ta1_z_0_yyzz_0[i] * pa_z[i] -
                            ta1_z_0_yyzz_1[i] * pc_z[i];

        ta1_z_z_yzzz_0[i] = 3.0 * ta1_z_0_yzz_0[i] * fe_0 - 3.0 * ta1_z_0_yzz_1[i] * fe_0 + ta_0_yzzz_1[i] + ta1_z_0_yzzz_0[i] * pa_z[i] -
                            ta1_z_0_yzzz_1[i] * pc_z[i];

        ta1_z_z_zzzz_0[i] = 4.0 * ta1_z_0_zzz_0[i] * fe_0 - 4.0 * ta1_z_0_zzz_1[i] * fe_0 + ta_0_zzzz_1[i] + ta1_z_0_zzzz_0[i] * pa_z[i] -
                            ta1_z_0_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
