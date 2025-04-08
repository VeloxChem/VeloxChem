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

#include "NuclearPotentialGeom010PrimRecSG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_sg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_sg,
                                        const size_t              idx_npot_geom_010_0_sd,
                                        const size_t              idx_npot_geom_010_1_sd,
                                        const size_t              idx_npot_1_sf,
                                        const size_t              idx_npot_geom_010_0_sf,
                                        const size_t              idx_npot_geom_010_1_sf,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpb,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd);

    auto ta1_x_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 3);

    auto ta1_x_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 5);

    auto ta1_y_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 6);

    auto ta1_y_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 9);

    auto ta1_y_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 11);

    auto ta1_z_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 12);

    auto ta1_z_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 15);

    auto ta1_z_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 17);

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd);

    auto ta1_x_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 3);

    auto ta1_x_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 5);

    auto ta1_y_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 6);

    auto ta1_y_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 9);

    auto ta1_y_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 11);

    auto ta1_z_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 12);

    auto ta1_z_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 15);

    auto ta1_z_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 17);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_1 = pbuffer.data(idx_npot_1_sf);

    auto ta_0_yyy_1 = pbuffer.data(idx_npot_1_sf + 6);

    auto ta_0_zzz_1 = pbuffer.data(idx_npot_1_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf);

    auto ta1_x_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 1);

    auto ta1_x_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 2);

    auto ta1_x_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 3);

    auto ta1_x_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 5);

    auto ta1_x_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 6);

    auto ta1_x_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 8);

    auto ta1_x_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 9);

    auto ta1_y_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 10);

    auto ta1_y_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 11);

    auto ta1_y_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 13);

    auto ta1_y_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 15);

    auto ta1_y_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 16);

    auto ta1_y_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 17);

    auto ta1_y_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 18);

    auto ta1_y_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 19);

    auto ta1_z_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 20);

    auto ta1_z_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 22);

    auto ta1_z_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 23);

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

    auto ta1_x_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 5);

    auto ta1_x_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 6);

    auto ta1_x_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 8);

    auto ta1_x_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 9);

    auto ta1_y_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 10);

    auto ta1_y_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 11);

    auto ta1_y_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 13);

    auto ta1_y_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 15);

    auto ta1_y_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 16);

    auto ta1_y_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 17);

    auto ta1_y_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 18);

    auto ta1_y_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 19);

    auto ta1_z_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 20);

    auto ta1_z_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 22);

    auto ta1_z_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 23);

    auto ta1_z_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 25);

    auto ta1_z_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 26);

    auto ta1_z_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 27);

    auto ta1_z_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 28);

    auto ta1_z_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 29);

    // Set up components of targeted buffer : SG

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

#pragma omp simd aligned(pb_x,               \
                             pb_y,           \
                             pb_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_0_xx_0,   \
                             ta1_x_0_xx_1,   \
                             ta1_x_0_xxx_0,  \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxxx_0, \
                             ta1_x_0_xxxy_0, \
                             ta1_x_0_xxxz_0, \
                             ta1_x_0_xxy_0,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxyy_0, \
                             ta1_x_0_xxyz_0, \
                             ta1_x_0_xxz_0,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xxzz_0, \
                             ta1_x_0_xyy_0,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyyy_0, \
                             ta1_x_0_xyyz_0, \
                             ta1_x_0_xyzz_0, \
                             ta1_x_0_xzz_0,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_xzzz_0, \
                             ta1_x_0_yy_0,   \
                             ta1_x_0_yy_1,   \
                             ta1_x_0_yyy_0,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyyy_0, \
                             ta1_x_0_yyyz_0, \
                             ta1_x_0_yyzz_0, \
                             ta1_x_0_yzz_0,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_yzzz_0, \
                             ta1_x_0_zz_0,   \
                             ta1_x_0_zz_1,   \
                             ta1_x_0_zzz_0,  \
                             ta1_x_0_zzz_1,  \
                             ta1_x_0_zzzz_0, \
                             ta1_y_0_xx_0,   \
                             ta1_y_0_xx_1,   \
                             ta1_y_0_xxx_0,  \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxxx_0, \
                             ta1_y_0_xxxy_0, \
                             ta1_y_0_xxxz_0, \
                             ta1_y_0_xxy_0,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxyy_0, \
                             ta1_y_0_xxyz_0, \
                             ta1_y_0_xxzz_0, \
                             ta1_y_0_xyy_0,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyyy_0, \
                             ta1_y_0_xyyz_0, \
                             ta1_y_0_xyzz_0, \
                             ta1_y_0_xzz_0,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_xzzz_0, \
                             ta1_y_0_yy_0,   \
                             ta1_y_0_yy_1,   \
                             ta1_y_0_yyy_0,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyyy_0, \
                             ta1_y_0_yyyz_0, \
                             ta1_y_0_yyz_0,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yyzz_0, \
                             ta1_y_0_yzz_0,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_yzzz_0, \
                             ta1_y_0_zz_0,   \
                             ta1_y_0_zz_1,   \
                             ta1_y_0_zzz_0,  \
                             ta1_y_0_zzz_1,  \
                             ta1_y_0_zzzz_0, \
                             ta1_z_0_xx_0,   \
                             ta1_z_0_xx_1,   \
                             ta1_z_0_xxx_0,  \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxxx_0, \
                             ta1_z_0_xxxy_0, \
                             ta1_z_0_xxxz_0, \
                             ta1_z_0_xxyy_0, \
                             ta1_z_0_xxyz_0, \
                             ta1_z_0_xxz_0,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xxzz_0, \
                             ta1_z_0_xyy_0,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyyy_0, \
                             ta1_z_0_xyyz_0, \
                             ta1_z_0_xyzz_0, \
                             ta1_z_0_xzz_0,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_xzzz_0, \
                             ta1_z_0_yy_0,   \
                             ta1_z_0_yy_1,   \
                             ta1_z_0_yyy_0,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyyy_0, \
                             ta1_z_0_yyyz_0, \
                             ta1_z_0_yyz_0,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yyzz_0, \
                             ta1_z_0_yzz_0,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_yzzz_0, \
                             ta1_z_0_zz_0,   \
                             ta1_z_0_zz_1,   \
                             ta1_z_0_zzz_0,  \
                             ta1_z_0_zzz_1,  \
                             ta1_z_0_zzzz_0, \
                             ta_0_xxx_1,     \
                             ta_0_yyy_1,     \
                             ta_0_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_0_xxxx_0[i] =
            3.0 * ta1_x_0_xx_0[i] * fe_0 - 3.0 * ta1_x_0_xx_1[i] * fe_0 + ta_0_xxx_1[i] + ta1_x_0_xxx_0[i] * pb_x[i] - ta1_x_0_xxx_1[i] * pc_x[i];

        ta1_x_0_xxxy_0[i] = ta1_x_0_xxx_0[i] * pb_y[i] - ta1_x_0_xxx_1[i] * pc_y[i];

        ta1_x_0_xxxz_0[i] = ta1_x_0_xxx_0[i] * pb_z[i] - ta1_x_0_xxx_1[i] * pc_z[i];

        ta1_x_0_xxyy_0[i] = ta1_x_0_xx_0[i] * fe_0 - ta1_x_0_xx_1[i] * fe_0 + ta1_x_0_xxy_0[i] * pb_y[i] - ta1_x_0_xxy_1[i] * pc_y[i];

        ta1_x_0_xxyz_0[i] = ta1_x_0_xxz_0[i] * pb_y[i] - ta1_x_0_xxz_1[i] * pc_y[i];

        ta1_x_0_xxzz_0[i] = ta1_x_0_xx_0[i] * fe_0 - ta1_x_0_xx_1[i] * fe_0 + ta1_x_0_xxz_0[i] * pb_z[i] - ta1_x_0_xxz_1[i] * pc_z[i];

        ta1_x_0_xyyy_0[i] = ta_0_yyy_1[i] + ta1_x_0_yyy_0[i] * pb_x[i] - ta1_x_0_yyy_1[i] * pc_x[i];

        ta1_x_0_xyyz_0[i] = ta1_x_0_xyy_0[i] * pb_z[i] - ta1_x_0_xyy_1[i] * pc_z[i];

        ta1_x_0_xyzz_0[i] = ta1_x_0_xzz_0[i] * pb_y[i] - ta1_x_0_xzz_1[i] * pc_y[i];

        ta1_x_0_xzzz_0[i] = ta_0_zzz_1[i] + ta1_x_0_zzz_0[i] * pb_x[i] - ta1_x_0_zzz_1[i] * pc_x[i];

        ta1_x_0_yyyy_0[i] = 3.0 * ta1_x_0_yy_0[i] * fe_0 - 3.0 * ta1_x_0_yy_1[i] * fe_0 + ta1_x_0_yyy_0[i] * pb_y[i] - ta1_x_0_yyy_1[i] * pc_y[i];

        ta1_x_0_yyyz_0[i] = ta1_x_0_yyy_0[i] * pb_z[i] - ta1_x_0_yyy_1[i] * pc_z[i];

        ta1_x_0_yyzz_0[i] = ta1_x_0_zz_0[i] * fe_0 - ta1_x_0_zz_1[i] * fe_0 + ta1_x_0_yzz_0[i] * pb_y[i] - ta1_x_0_yzz_1[i] * pc_y[i];

        ta1_x_0_yzzz_0[i] = ta1_x_0_zzz_0[i] * pb_y[i] - ta1_x_0_zzz_1[i] * pc_y[i];

        ta1_x_0_zzzz_0[i] = 3.0 * ta1_x_0_zz_0[i] * fe_0 - 3.0 * ta1_x_0_zz_1[i] * fe_0 + ta1_x_0_zzz_0[i] * pb_z[i] - ta1_x_0_zzz_1[i] * pc_z[i];

        ta1_y_0_xxxx_0[i] = 3.0 * ta1_y_0_xx_0[i] * fe_0 - 3.0 * ta1_y_0_xx_1[i] * fe_0 + ta1_y_0_xxx_0[i] * pb_x[i] - ta1_y_0_xxx_1[i] * pc_x[i];

        ta1_y_0_xxxy_0[i] = ta_0_xxx_1[i] + ta1_y_0_xxx_0[i] * pb_y[i] - ta1_y_0_xxx_1[i] * pc_y[i];

        ta1_y_0_xxxz_0[i] = ta1_y_0_xxx_0[i] * pb_z[i] - ta1_y_0_xxx_1[i] * pc_z[i];

        ta1_y_0_xxyy_0[i] = ta1_y_0_yy_0[i] * fe_0 - ta1_y_0_yy_1[i] * fe_0 + ta1_y_0_xyy_0[i] * pb_x[i] - ta1_y_0_xyy_1[i] * pc_x[i];

        ta1_y_0_xxyz_0[i] = ta1_y_0_xxy_0[i] * pb_z[i] - ta1_y_0_xxy_1[i] * pc_z[i];

        ta1_y_0_xxzz_0[i] = ta1_y_0_zz_0[i] * fe_0 - ta1_y_0_zz_1[i] * fe_0 + ta1_y_0_xzz_0[i] * pb_x[i] - ta1_y_0_xzz_1[i] * pc_x[i];

        ta1_y_0_xyyy_0[i] = ta1_y_0_yyy_0[i] * pb_x[i] - ta1_y_0_yyy_1[i] * pc_x[i];

        ta1_y_0_xyyz_0[i] = ta1_y_0_yyz_0[i] * pb_x[i] - ta1_y_0_yyz_1[i] * pc_x[i];

        ta1_y_0_xyzz_0[i] = ta1_y_0_yzz_0[i] * pb_x[i] - ta1_y_0_yzz_1[i] * pc_x[i];

        ta1_y_0_xzzz_0[i] = ta1_y_0_zzz_0[i] * pb_x[i] - ta1_y_0_zzz_1[i] * pc_x[i];

        ta1_y_0_yyyy_0[i] =
            3.0 * ta1_y_0_yy_0[i] * fe_0 - 3.0 * ta1_y_0_yy_1[i] * fe_0 + ta_0_yyy_1[i] + ta1_y_0_yyy_0[i] * pb_y[i] - ta1_y_0_yyy_1[i] * pc_y[i];

        ta1_y_0_yyyz_0[i] = ta1_y_0_yyy_0[i] * pb_z[i] - ta1_y_0_yyy_1[i] * pc_z[i];

        ta1_y_0_yyzz_0[i] = ta1_y_0_yy_0[i] * fe_0 - ta1_y_0_yy_1[i] * fe_0 + ta1_y_0_yyz_0[i] * pb_z[i] - ta1_y_0_yyz_1[i] * pc_z[i];

        ta1_y_0_yzzz_0[i] = ta_0_zzz_1[i] + ta1_y_0_zzz_0[i] * pb_y[i] - ta1_y_0_zzz_1[i] * pc_y[i];

        ta1_y_0_zzzz_0[i] = 3.0 * ta1_y_0_zz_0[i] * fe_0 - 3.0 * ta1_y_0_zz_1[i] * fe_0 + ta1_y_0_zzz_0[i] * pb_z[i] - ta1_y_0_zzz_1[i] * pc_z[i];

        ta1_z_0_xxxx_0[i] = 3.0 * ta1_z_0_xx_0[i] * fe_0 - 3.0 * ta1_z_0_xx_1[i] * fe_0 + ta1_z_0_xxx_0[i] * pb_x[i] - ta1_z_0_xxx_1[i] * pc_x[i];

        ta1_z_0_xxxy_0[i] = ta1_z_0_xxx_0[i] * pb_y[i] - ta1_z_0_xxx_1[i] * pc_y[i];

        ta1_z_0_xxxz_0[i] = ta_0_xxx_1[i] + ta1_z_0_xxx_0[i] * pb_z[i] - ta1_z_0_xxx_1[i] * pc_z[i];

        ta1_z_0_xxyy_0[i] = ta1_z_0_yy_0[i] * fe_0 - ta1_z_0_yy_1[i] * fe_0 + ta1_z_0_xyy_0[i] * pb_x[i] - ta1_z_0_xyy_1[i] * pc_x[i];

        ta1_z_0_xxyz_0[i] = ta1_z_0_xxz_0[i] * pb_y[i] - ta1_z_0_xxz_1[i] * pc_y[i];

        ta1_z_0_xxzz_0[i] = ta1_z_0_zz_0[i] * fe_0 - ta1_z_0_zz_1[i] * fe_0 + ta1_z_0_xzz_0[i] * pb_x[i] - ta1_z_0_xzz_1[i] * pc_x[i];

        ta1_z_0_xyyy_0[i] = ta1_z_0_yyy_0[i] * pb_x[i] - ta1_z_0_yyy_1[i] * pc_x[i];

        ta1_z_0_xyyz_0[i] = ta1_z_0_yyz_0[i] * pb_x[i] - ta1_z_0_yyz_1[i] * pc_x[i];

        ta1_z_0_xyzz_0[i] = ta1_z_0_yzz_0[i] * pb_x[i] - ta1_z_0_yzz_1[i] * pc_x[i];

        ta1_z_0_xzzz_0[i] = ta1_z_0_zzz_0[i] * pb_x[i] - ta1_z_0_zzz_1[i] * pc_x[i];

        ta1_z_0_yyyy_0[i] = 3.0 * ta1_z_0_yy_0[i] * fe_0 - 3.0 * ta1_z_0_yy_1[i] * fe_0 + ta1_z_0_yyy_0[i] * pb_y[i] - ta1_z_0_yyy_1[i] * pc_y[i];

        ta1_z_0_yyyz_0[i] = ta_0_yyy_1[i] + ta1_z_0_yyy_0[i] * pb_z[i] - ta1_z_0_yyy_1[i] * pc_z[i];

        ta1_z_0_yyzz_0[i] = ta1_z_0_zz_0[i] * fe_0 - ta1_z_0_zz_1[i] * fe_0 + ta1_z_0_yzz_0[i] * pb_y[i] - ta1_z_0_yzz_1[i] * pc_y[i];

        ta1_z_0_yzzz_0[i] = ta1_z_0_zzz_0[i] * pb_y[i] - ta1_z_0_zzz_1[i] * pc_y[i];

        ta1_z_0_zzzz_0[i] =
            3.0 * ta1_z_0_zz_0[i] * fe_0 - 3.0 * ta1_z_0_zz_1[i] * fe_0 + ta_0_zzz_1[i] + ta1_z_0_zzz_0[i] * pb_z[i] - ta1_z_0_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
