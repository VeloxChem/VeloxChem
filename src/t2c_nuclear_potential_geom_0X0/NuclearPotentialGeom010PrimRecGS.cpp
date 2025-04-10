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

#include "NuclearPotentialGeom010PrimRecGS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_gs(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_gs,
                                        const size_t              idx_npot_geom_010_0_ds,
                                        const size_t              idx_npot_geom_010_1_ds,
                                        const size_t              idx_npot_1_fs,
                                        const size_t              idx_npot_geom_010_0_fs,
                                        const size_t              idx_npot_geom_010_1_fs,
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

    // Set up components of auxiliary buffer : DS

    auto ta1_x_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds);

    auto ta1_x_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 3);

    auto ta1_x_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 5);

    auto ta1_y_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 6);

    auto ta1_y_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 9);

    auto ta1_y_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 11);

    auto ta1_z_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 12);

    auto ta1_z_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 15);

    auto ta1_z_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 17);

    // Set up components of auxiliary buffer : DS

    auto ta1_x_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds);

    auto ta1_x_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 3);

    auto ta1_x_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 5);

    auto ta1_y_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 6);

    auto ta1_y_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 9);

    auto ta1_y_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 11);

    auto ta1_z_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 12);

    auto ta1_z_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 15);

    auto ta1_z_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 17);

    // Set up components of auxiliary buffer : FS

    auto ta_xxx_0_1 = pbuffer.data(idx_npot_1_fs);

    auto ta_yyy_0_1 = pbuffer.data(idx_npot_1_fs + 6);

    auto ta_zzz_0_1 = pbuffer.data(idx_npot_1_fs + 9);

    // Set up components of auxiliary buffer : FS

    auto ta1_x_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs);

    auto ta1_x_xxy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 1);

    auto ta1_x_xxz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 2);

    auto ta1_x_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 3);

    auto ta1_x_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 5);

    auto ta1_x_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 6);

    auto ta1_x_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 8);

    auto ta1_x_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 9);

    auto ta1_y_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 10);

    auto ta1_y_xxy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 11);

    auto ta1_y_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 13);

    auto ta1_y_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 15);

    auto ta1_y_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 16);

    auto ta1_y_yyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 17);

    auto ta1_y_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 18);

    auto ta1_y_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 19);

    auto ta1_z_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 20);

    auto ta1_z_xxz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 22);

    auto ta1_z_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 23);

    auto ta1_z_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 25);

    auto ta1_z_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 26);

    auto ta1_z_yyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 27);

    auto ta1_z_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 28);

    auto ta1_z_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 29);

    // Set up components of auxiliary buffer : FS

    auto ta1_x_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs);

    auto ta1_x_xxy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 1);

    auto ta1_x_xxz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 2);

    auto ta1_x_xyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 3);

    auto ta1_x_xzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 5);

    auto ta1_x_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 6);

    auto ta1_x_yzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 8);

    auto ta1_x_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 9);

    auto ta1_y_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 10);

    auto ta1_y_xxy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 11);

    auto ta1_y_xyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 13);

    auto ta1_y_xzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 15);

    auto ta1_y_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 16);

    auto ta1_y_yyz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 17);

    auto ta1_y_yzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 18);

    auto ta1_y_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 19);

    auto ta1_z_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 20);

    auto ta1_z_xxz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 22);

    auto ta1_z_xyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 23);

    auto ta1_z_xzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 25);

    auto ta1_z_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 26);

    auto ta1_z_yyz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 27);

    auto ta1_z_yzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 28);

    auto ta1_z_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 29);

    // Set up components of targeted buffer : GS

    auto ta1_x_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs);

    auto ta1_x_xxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 1);

    auto ta1_x_xxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 2);

    auto ta1_x_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 3);

    auto ta1_x_xxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 4);

    auto ta1_x_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 5);

    auto ta1_x_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 6);

    auto ta1_x_xyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 7);

    auto ta1_x_xyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 8);

    auto ta1_x_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 9);

    auto ta1_x_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 10);

    auto ta1_x_yyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 11);

    auto ta1_x_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 12);

    auto ta1_x_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 13);

    auto ta1_x_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 14);

    auto ta1_y_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 15);

    auto ta1_y_xxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 16);

    auto ta1_y_xxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 17);

    auto ta1_y_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 18);

    auto ta1_y_xxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 19);

    auto ta1_y_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 20);

    auto ta1_y_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 21);

    auto ta1_y_xyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 22);

    auto ta1_y_xyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 23);

    auto ta1_y_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 24);

    auto ta1_y_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 25);

    auto ta1_y_yyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 26);

    auto ta1_y_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 27);

    auto ta1_y_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 28);

    auto ta1_y_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 29);

    auto ta1_z_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 30);

    auto ta1_z_xxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 31);

    auto ta1_z_xxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 32);

    auto ta1_z_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 33);

    auto ta1_z_xxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 34);

    auto ta1_z_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 35);

    auto ta1_z_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 36);

    auto ta1_z_xyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 37);

    auto ta1_z_xyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 38);

    auto ta1_z_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 39);

    auto ta1_z_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 40);

    auto ta1_z_yyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 41);

    auto ta1_z_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 42);

    auto ta1_z_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 43);

    auto ta1_z_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 44);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_xx_0_0,   \
                             ta1_x_xx_0_1,   \
                             ta1_x_xxx_0_0,  \
                             ta1_x_xxx_0_1,  \
                             ta1_x_xxxx_0_0, \
                             ta1_x_xxxy_0_0, \
                             ta1_x_xxxz_0_0, \
                             ta1_x_xxy_0_0,  \
                             ta1_x_xxy_0_1,  \
                             ta1_x_xxyy_0_0, \
                             ta1_x_xxyz_0_0, \
                             ta1_x_xxz_0_0,  \
                             ta1_x_xxz_0_1,  \
                             ta1_x_xxzz_0_0, \
                             ta1_x_xyy_0_0,  \
                             ta1_x_xyy_0_1,  \
                             ta1_x_xyyy_0_0, \
                             ta1_x_xyyz_0_0, \
                             ta1_x_xyzz_0_0, \
                             ta1_x_xzz_0_0,  \
                             ta1_x_xzz_0_1,  \
                             ta1_x_xzzz_0_0, \
                             ta1_x_yy_0_0,   \
                             ta1_x_yy_0_1,   \
                             ta1_x_yyy_0_0,  \
                             ta1_x_yyy_0_1,  \
                             ta1_x_yyyy_0_0, \
                             ta1_x_yyyz_0_0, \
                             ta1_x_yyzz_0_0, \
                             ta1_x_yzz_0_0,  \
                             ta1_x_yzz_0_1,  \
                             ta1_x_yzzz_0_0, \
                             ta1_x_zz_0_0,   \
                             ta1_x_zz_0_1,   \
                             ta1_x_zzz_0_0,  \
                             ta1_x_zzz_0_1,  \
                             ta1_x_zzzz_0_0, \
                             ta1_y_xx_0_0,   \
                             ta1_y_xx_0_1,   \
                             ta1_y_xxx_0_0,  \
                             ta1_y_xxx_0_1,  \
                             ta1_y_xxxx_0_0, \
                             ta1_y_xxxy_0_0, \
                             ta1_y_xxxz_0_0, \
                             ta1_y_xxy_0_0,  \
                             ta1_y_xxy_0_1,  \
                             ta1_y_xxyy_0_0, \
                             ta1_y_xxyz_0_0, \
                             ta1_y_xxzz_0_0, \
                             ta1_y_xyy_0_0,  \
                             ta1_y_xyy_0_1,  \
                             ta1_y_xyyy_0_0, \
                             ta1_y_xyyz_0_0, \
                             ta1_y_xyzz_0_0, \
                             ta1_y_xzz_0_0,  \
                             ta1_y_xzz_0_1,  \
                             ta1_y_xzzz_0_0, \
                             ta1_y_yy_0_0,   \
                             ta1_y_yy_0_1,   \
                             ta1_y_yyy_0_0,  \
                             ta1_y_yyy_0_1,  \
                             ta1_y_yyyy_0_0, \
                             ta1_y_yyyz_0_0, \
                             ta1_y_yyz_0_0,  \
                             ta1_y_yyz_0_1,  \
                             ta1_y_yyzz_0_0, \
                             ta1_y_yzz_0_0,  \
                             ta1_y_yzz_0_1,  \
                             ta1_y_yzzz_0_0, \
                             ta1_y_zz_0_0,   \
                             ta1_y_zz_0_1,   \
                             ta1_y_zzz_0_0,  \
                             ta1_y_zzz_0_1,  \
                             ta1_y_zzzz_0_0, \
                             ta1_z_xx_0_0,   \
                             ta1_z_xx_0_1,   \
                             ta1_z_xxx_0_0,  \
                             ta1_z_xxx_0_1,  \
                             ta1_z_xxxx_0_0, \
                             ta1_z_xxxy_0_0, \
                             ta1_z_xxxz_0_0, \
                             ta1_z_xxyy_0_0, \
                             ta1_z_xxyz_0_0, \
                             ta1_z_xxz_0_0,  \
                             ta1_z_xxz_0_1,  \
                             ta1_z_xxzz_0_0, \
                             ta1_z_xyy_0_0,  \
                             ta1_z_xyy_0_1,  \
                             ta1_z_xyyy_0_0, \
                             ta1_z_xyyz_0_0, \
                             ta1_z_xyzz_0_0, \
                             ta1_z_xzz_0_0,  \
                             ta1_z_xzz_0_1,  \
                             ta1_z_xzzz_0_0, \
                             ta1_z_yy_0_0,   \
                             ta1_z_yy_0_1,   \
                             ta1_z_yyy_0_0,  \
                             ta1_z_yyy_0_1,  \
                             ta1_z_yyyy_0_0, \
                             ta1_z_yyyz_0_0, \
                             ta1_z_yyz_0_0,  \
                             ta1_z_yyz_0_1,  \
                             ta1_z_yyzz_0_0, \
                             ta1_z_yzz_0_0,  \
                             ta1_z_yzz_0_1,  \
                             ta1_z_yzzz_0_0, \
                             ta1_z_zz_0_0,   \
                             ta1_z_zz_0_1,   \
                             ta1_z_zzz_0_0,  \
                             ta1_z_zzz_0_1,  \
                             ta1_z_zzzz_0_0, \
                             ta_xxx_0_1,     \
                             ta_yyy_0_1,     \
                             ta_zzz_0_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxx_0_0[i] =
            3.0 * ta1_x_xx_0_0[i] * fe_0 - 3.0 * ta1_x_xx_0_1[i] * fe_0 + ta_xxx_0_1[i] + ta1_x_xxx_0_0[i] * pa_x[i] - ta1_x_xxx_0_1[i] * pc_x[i];

        ta1_x_xxxy_0_0[i] = ta1_x_xxx_0_0[i] * pa_y[i] - ta1_x_xxx_0_1[i] * pc_y[i];

        ta1_x_xxxz_0_0[i] = ta1_x_xxx_0_0[i] * pa_z[i] - ta1_x_xxx_0_1[i] * pc_z[i];

        ta1_x_xxyy_0_0[i] = ta1_x_xx_0_0[i] * fe_0 - ta1_x_xx_0_1[i] * fe_0 + ta1_x_xxy_0_0[i] * pa_y[i] - ta1_x_xxy_0_1[i] * pc_y[i];

        ta1_x_xxyz_0_0[i] = ta1_x_xxz_0_0[i] * pa_y[i] - ta1_x_xxz_0_1[i] * pc_y[i];

        ta1_x_xxzz_0_0[i] = ta1_x_xx_0_0[i] * fe_0 - ta1_x_xx_0_1[i] * fe_0 + ta1_x_xxz_0_0[i] * pa_z[i] - ta1_x_xxz_0_1[i] * pc_z[i];

        ta1_x_xyyy_0_0[i] = ta_yyy_0_1[i] + ta1_x_yyy_0_0[i] * pa_x[i] - ta1_x_yyy_0_1[i] * pc_x[i];

        ta1_x_xyyz_0_0[i] = ta1_x_xyy_0_0[i] * pa_z[i] - ta1_x_xyy_0_1[i] * pc_z[i];

        ta1_x_xyzz_0_0[i] = ta1_x_xzz_0_0[i] * pa_y[i] - ta1_x_xzz_0_1[i] * pc_y[i];

        ta1_x_xzzz_0_0[i] = ta_zzz_0_1[i] + ta1_x_zzz_0_0[i] * pa_x[i] - ta1_x_zzz_0_1[i] * pc_x[i];

        ta1_x_yyyy_0_0[i] = 3.0 * ta1_x_yy_0_0[i] * fe_0 - 3.0 * ta1_x_yy_0_1[i] * fe_0 + ta1_x_yyy_0_0[i] * pa_y[i] - ta1_x_yyy_0_1[i] * pc_y[i];

        ta1_x_yyyz_0_0[i] = ta1_x_yyy_0_0[i] * pa_z[i] - ta1_x_yyy_0_1[i] * pc_z[i];

        ta1_x_yyzz_0_0[i] = ta1_x_zz_0_0[i] * fe_0 - ta1_x_zz_0_1[i] * fe_0 + ta1_x_yzz_0_0[i] * pa_y[i] - ta1_x_yzz_0_1[i] * pc_y[i];

        ta1_x_yzzz_0_0[i] = ta1_x_zzz_0_0[i] * pa_y[i] - ta1_x_zzz_0_1[i] * pc_y[i];

        ta1_x_zzzz_0_0[i] = 3.0 * ta1_x_zz_0_0[i] * fe_0 - 3.0 * ta1_x_zz_0_1[i] * fe_0 + ta1_x_zzz_0_0[i] * pa_z[i] - ta1_x_zzz_0_1[i] * pc_z[i];

        ta1_y_xxxx_0_0[i] = 3.0 * ta1_y_xx_0_0[i] * fe_0 - 3.0 * ta1_y_xx_0_1[i] * fe_0 + ta1_y_xxx_0_0[i] * pa_x[i] - ta1_y_xxx_0_1[i] * pc_x[i];

        ta1_y_xxxy_0_0[i] = ta_xxx_0_1[i] + ta1_y_xxx_0_0[i] * pa_y[i] - ta1_y_xxx_0_1[i] * pc_y[i];

        ta1_y_xxxz_0_0[i] = ta1_y_xxx_0_0[i] * pa_z[i] - ta1_y_xxx_0_1[i] * pc_z[i];

        ta1_y_xxyy_0_0[i] = ta1_y_yy_0_0[i] * fe_0 - ta1_y_yy_0_1[i] * fe_0 + ta1_y_xyy_0_0[i] * pa_x[i] - ta1_y_xyy_0_1[i] * pc_x[i];

        ta1_y_xxyz_0_0[i] = ta1_y_xxy_0_0[i] * pa_z[i] - ta1_y_xxy_0_1[i] * pc_z[i];

        ta1_y_xxzz_0_0[i] = ta1_y_zz_0_0[i] * fe_0 - ta1_y_zz_0_1[i] * fe_0 + ta1_y_xzz_0_0[i] * pa_x[i] - ta1_y_xzz_0_1[i] * pc_x[i];

        ta1_y_xyyy_0_0[i] = ta1_y_yyy_0_0[i] * pa_x[i] - ta1_y_yyy_0_1[i] * pc_x[i];

        ta1_y_xyyz_0_0[i] = ta1_y_yyz_0_0[i] * pa_x[i] - ta1_y_yyz_0_1[i] * pc_x[i];

        ta1_y_xyzz_0_0[i] = ta1_y_yzz_0_0[i] * pa_x[i] - ta1_y_yzz_0_1[i] * pc_x[i];

        ta1_y_xzzz_0_0[i] = ta1_y_zzz_0_0[i] * pa_x[i] - ta1_y_zzz_0_1[i] * pc_x[i];

        ta1_y_yyyy_0_0[i] =
            3.0 * ta1_y_yy_0_0[i] * fe_0 - 3.0 * ta1_y_yy_0_1[i] * fe_0 + ta_yyy_0_1[i] + ta1_y_yyy_0_0[i] * pa_y[i] - ta1_y_yyy_0_1[i] * pc_y[i];

        ta1_y_yyyz_0_0[i] = ta1_y_yyy_0_0[i] * pa_z[i] - ta1_y_yyy_0_1[i] * pc_z[i];

        ta1_y_yyzz_0_0[i] = ta1_y_yy_0_0[i] * fe_0 - ta1_y_yy_0_1[i] * fe_0 + ta1_y_yyz_0_0[i] * pa_z[i] - ta1_y_yyz_0_1[i] * pc_z[i];

        ta1_y_yzzz_0_0[i] = ta_zzz_0_1[i] + ta1_y_zzz_0_0[i] * pa_y[i] - ta1_y_zzz_0_1[i] * pc_y[i];

        ta1_y_zzzz_0_0[i] = 3.0 * ta1_y_zz_0_0[i] * fe_0 - 3.0 * ta1_y_zz_0_1[i] * fe_0 + ta1_y_zzz_0_0[i] * pa_z[i] - ta1_y_zzz_0_1[i] * pc_z[i];

        ta1_z_xxxx_0_0[i] = 3.0 * ta1_z_xx_0_0[i] * fe_0 - 3.0 * ta1_z_xx_0_1[i] * fe_0 + ta1_z_xxx_0_0[i] * pa_x[i] - ta1_z_xxx_0_1[i] * pc_x[i];

        ta1_z_xxxy_0_0[i] = ta1_z_xxx_0_0[i] * pa_y[i] - ta1_z_xxx_0_1[i] * pc_y[i];

        ta1_z_xxxz_0_0[i] = ta_xxx_0_1[i] + ta1_z_xxx_0_0[i] * pa_z[i] - ta1_z_xxx_0_1[i] * pc_z[i];

        ta1_z_xxyy_0_0[i] = ta1_z_yy_0_0[i] * fe_0 - ta1_z_yy_0_1[i] * fe_0 + ta1_z_xyy_0_0[i] * pa_x[i] - ta1_z_xyy_0_1[i] * pc_x[i];

        ta1_z_xxyz_0_0[i] = ta1_z_xxz_0_0[i] * pa_y[i] - ta1_z_xxz_0_1[i] * pc_y[i];

        ta1_z_xxzz_0_0[i] = ta1_z_zz_0_0[i] * fe_0 - ta1_z_zz_0_1[i] * fe_0 + ta1_z_xzz_0_0[i] * pa_x[i] - ta1_z_xzz_0_1[i] * pc_x[i];

        ta1_z_xyyy_0_0[i] = ta1_z_yyy_0_0[i] * pa_x[i] - ta1_z_yyy_0_1[i] * pc_x[i];

        ta1_z_xyyz_0_0[i] = ta1_z_yyz_0_0[i] * pa_x[i] - ta1_z_yyz_0_1[i] * pc_x[i];

        ta1_z_xyzz_0_0[i] = ta1_z_yzz_0_0[i] * pa_x[i] - ta1_z_yzz_0_1[i] * pc_x[i];

        ta1_z_xzzz_0_0[i] = ta1_z_zzz_0_0[i] * pa_x[i] - ta1_z_zzz_0_1[i] * pc_x[i];

        ta1_z_yyyy_0_0[i] = 3.0 * ta1_z_yy_0_0[i] * fe_0 - 3.0 * ta1_z_yy_0_1[i] * fe_0 + ta1_z_yyy_0_0[i] * pa_y[i] - ta1_z_yyy_0_1[i] * pc_y[i];

        ta1_z_yyyz_0_0[i] = ta_yyy_0_1[i] + ta1_z_yyy_0_0[i] * pa_z[i] - ta1_z_yyy_0_1[i] * pc_z[i];

        ta1_z_yyzz_0_0[i] = ta1_z_zz_0_0[i] * fe_0 - ta1_z_zz_0_1[i] * fe_0 + ta1_z_yzz_0_0[i] * pa_y[i] - ta1_z_yzz_0_1[i] * pc_y[i];

        ta1_z_yzzz_0_0[i] = ta1_z_zzz_0_0[i] * pa_y[i] - ta1_z_zzz_0_1[i] * pc_y[i];

        ta1_z_zzzz_0_0[i] =
            3.0 * ta1_z_zz_0_0[i] * fe_0 - 3.0 * ta1_z_zz_0_1[i] * fe_0 + ta_zzz_0_1[i] + ta1_z_zzz_0_0[i] * pa_z[i] - ta1_z_zzz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
