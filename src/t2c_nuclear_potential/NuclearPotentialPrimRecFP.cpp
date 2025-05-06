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

#include "NuclearPotentialPrimRecFP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_fp(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_fp,
                               const size_t              idx_npot_0_pp,
                               const size_t              idx_npot_1_pp,
                               const size_t              idx_npot_0_ds,
                               const size_t              idx_npot_1_ds,
                               const size_t              idx_npot_0_dp,
                               const size_t              idx_npot_1_dp,
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

    // Set up components of auxiliary buffer : PP

    auto ta_x_x_0 = pbuffer.data(idx_npot_0_pp);

    auto ta_x_y_0 = pbuffer.data(idx_npot_0_pp + 1);

    auto ta_x_z_0 = pbuffer.data(idx_npot_0_pp + 2);

    auto ta_y_x_0 = pbuffer.data(idx_npot_0_pp + 3);

    auto ta_y_y_0 = pbuffer.data(idx_npot_0_pp + 4);

    auto ta_y_z_0 = pbuffer.data(idx_npot_0_pp + 5);

    auto ta_z_x_0 = pbuffer.data(idx_npot_0_pp + 6);

    auto ta_z_y_0 = pbuffer.data(idx_npot_0_pp + 7);

    auto ta_z_z_0 = pbuffer.data(idx_npot_0_pp + 8);

    // Set up components of auxiliary buffer : PP

    auto ta_x_x_1 = pbuffer.data(idx_npot_1_pp);

    auto ta_x_y_1 = pbuffer.data(idx_npot_1_pp + 1);

    auto ta_x_z_1 = pbuffer.data(idx_npot_1_pp + 2);

    auto ta_y_x_1 = pbuffer.data(idx_npot_1_pp + 3);

    auto ta_y_y_1 = pbuffer.data(idx_npot_1_pp + 4);

    auto ta_y_z_1 = pbuffer.data(idx_npot_1_pp + 5);

    auto ta_z_x_1 = pbuffer.data(idx_npot_1_pp + 6);

    auto ta_z_y_1 = pbuffer.data(idx_npot_1_pp + 7);

    auto ta_z_z_1 = pbuffer.data(idx_npot_1_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_0 = pbuffer.data(idx_npot_0_ds);

    auto ta_yy_0_0 = pbuffer.data(idx_npot_0_ds + 3);

    auto ta_zz_0_0 = pbuffer.data(idx_npot_0_ds + 5);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_1 = pbuffer.data(idx_npot_1_ds);

    auto ta_yy_0_1 = pbuffer.data(idx_npot_1_ds + 3);

    auto ta_zz_0_1 = pbuffer.data(idx_npot_1_ds + 5);

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_0 = pbuffer.data(idx_npot_0_dp);

    auto ta_xx_y_0 = pbuffer.data(idx_npot_0_dp + 1);

    auto ta_xx_z_0 = pbuffer.data(idx_npot_0_dp + 2);

    auto ta_xy_y_0 = pbuffer.data(idx_npot_0_dp + 4);

    auto ta_xz_x_0 = pbuffer.data(idx_npot_0_dp + 6);

    auto ta_xz_z_0 = pbuffer.data(idx_npot_0_dp + 8);

    auto ta_yy_x_0 = pbuffer.data(idx_npot_0_dp + 9);

    auto ta_yy_y_0 = pbuffer.data(idx_npot_0_dp + 10);

    auto ta_yy_z_0 = pbuffer.data(idx_npot_0_dp + 11);

    auto ta_yz_y_0 = pbuffer.data(idx_npot_0_dp + 13);

    auto ta_yz_z_0 = pbuffer.data(idx_npot_0_dp + 14);

    auto ta_zz_x_0 = pbuffer.data(idx_npot_0_dp + 15);

    auto ta_zz_y_0 = pbuffer.data(idx_npot_0_dp + 16);

    auto ta_zz_z_0 = pbuffer.data(idx_npot_0_dp + 17);

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_1 = pbuffer.data(idx_npot_1_dp);

    auto ta_xx_y_1 = pbuffer.data(idx_npot_1_dp + 1);

    auto ta_xx_z_1 = pbuffer.data(idx_npot_1_dp + 2);

    auto ta_xy_y_1 = pbuffer.data(idx_npot_1_dp + 4);

    auto ta_xz_x_1 = pbuffer.data(idx_npot_1_dp + 6);

    auto ta_xz_z_1 = pbuffer.data(idx_npot_1_dp + 8);

    auto ta_yy_x_1 = pbuffer.data(idx_npot_1_dp + 9);

    auto ta_yy_y_1 = pbuffer.data(idx_npot_1_dp + 10);

    auto ta_yy_z_1 = pbuffer.data(idx_npot_1_dp + 11);

    auto ta_yz_y_1 = pbuffer.data(idx_npot_1_dp + 13);

    auto ta_yz_z_1 = pbuffer.data(idx_npot_1_dp + 14);

    auto ta_zz_x_1 = pbuffer.data(idx_npot_1_dp + 15);

    auto ta_zz_y_1 = pbuffer.data(idx_npot_1_dp + 16);

    auto ta_zz_z_1 = pbuffer.data(idx_npot_1_dp + 17);

    // Set up 0-3 components of targeted buffer : FP

    auto ta_xxx_x_0 = pbuffer.data(idx_npot_0_fp);

    auto ta_xxx_y_0 = pbuffer.data(idx_npot_0_fp + 1);

    auto ta_xxx_z_0 = pbuffer.data(idx_npot_0_fp + 2);

#pragma omp simd aligned(pa_x,           \
                             pc_x,       \
                             ta_x_x_0,   \
                             ta_x_x_1,   \
                             ta_x_y_0,   \
                             ta_x_y_1,   \
                             ta_x_z_0,   \
                             ta_x_z_1,   \
                             ta_xx_0_0,  \
                             ta_xx_0_1,  \
                             ta_xx_x_0,  \
                             ta_xx_x_1,  \
                             ta_xx_y_0,  \
                             ta_xx_y_1,  \
                             ta_xx_z_0,  \
                             ta_xx_z_1,  \
                             ta_xxx_x_0, \
                             ta_xxx_y_0, \
                             ta_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxx_x_0[i] = 2.0 * ta_x_x_0[i] * fe_0 - 2.0 * ta_x_x_1[i] * fe_0 + ta_xx_0_0[i] * fe_0 - ta_xx_0_1[i] * fe_0 + ta_xx_x_0[i] * pa_x[i] -
                        ta_xx_x_1[i] * pc_x[i];

        ta_xxx_y_0[i] = 2.0 * ta_x_y_0[i] * fe_0 - 2.0 * ta_x_y_1[i] * fe_0 + ta_xx_y_0[i] * pa_x[i] - ta_xx_y_1[i] * pc_x[i];

        ta_xxx_z_0[i] = 2.0 * ta_x_z_0[i] * fe_0 - 2.0 * ta_x_z_1[i] * fe_0 + ta_xx_z_0[i] * pa_x[i] - ta_xx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto ta_xxy_x_0 = pbuffer.data(idx_npot_0_fp + 3);

    auto ta_xxy_y_0 = pbuffer.data(idx_npot_0_fp + 4);

    auto ta_xxy_z_0 = pbuffer.data(idx_npot_0_fp + 5);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             pc_x,       \
                             pc_y,       \
                             ta_xx_x_0,  \
                             ta_xx_x_1,  \
                             ta_xx_z_0,  \
                             ta_xx_z_1,  \
                             ta_xxy_x_0, \
                             ta_xxy_y_0, \
                             ta_xxy_z_0, \
                             ta_xy_y_0,  \
                             ta_xy_y_1,  \
                             ta_y_y_0,   \
                             ta_y_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxy_x_0[i] = ta_xx_x_0[i] * pa_y[i] - ta_xx_x_1[i] * pc_y[i];

        ta_xxy_y_0[i] = ta_y_y_0[i] * fe_0 - ta_y_y_1[i] * fe_0 + ta_xy_y_0[i] * pa_x[i] - ta_xy_y_1[i] * pc_x[i];

        ta_xxy_z_0[i] = ta_xx_z_0[i] * pa_y[i] - ta_xx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto ta_xxz_x_0 = pbuffer.data(idx_npot_0_fp + 6);

    auto ta_xxz_y_0 = pbuffer.data(idx_npot_0_fp + 7);

    auto ta_xxz_z_0 = pbuffer.data(idx_npot_0_fp + 8);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             pc_x,       \
                             pc_z,       \
                             ta_xx_x_0,  \
                             ta_xx_x_1,  \
                             ta_xx_y_0,  \
                             ta_xx_y_1,  \
                             ta_xxz_x_0, \
                             ta_xxz_y_0, \
                             ta_xxz_z_0, \
                             ta_xz_z_0,  \
                             ta_xz_z_1,  \
                             ta_z_z_0,   \
                             ta_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxz_x_0[i] = ta_xx_x_0[i] * pa_z[i] - ta_xx_x_1[i] * pc_z[i];

        ta_xxz_y_0[i] = ta_xx_y_0[i] * pa_z[i] - ta_xx_y_1[i] * pc_z[i];

        ta_xxz_z_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_xz_z_0[i] * pa_x[i] - ta_xz_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto ta_xyy_x_0 = pbuffer.data(idx_npot_0_fp + 9);

    auto ta_xyy_y_0 = pbuffer.data(idx_npot_0_fp + 10);

    auto ta_xyy_z_0 = pbuffer.data(idx_npot_0_fp + 11);

#pragma omp simd aligned(pa_x,           \
                             pc_x,       \
                             ta_xyy_x_0, \
                             ta_xyy_y_0, \
                             ta_xyy_z_0, \
                             ta_yy_0_0,  \
                             ta_yy_0_1,  \
                             ta_yy_x_0,  \
                             ta_yy_x_1,  \
                             ta_yy_y_0,  \
                             ta_yy_y_1,  \
                             ta_yy_z_0,  \
                             ta_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyy_x_0[i] = ta_yy_0_0[i] * fe_0 - ta_yy_0_1[i] * fe_0 + ta_yy_x_0[i] * pa_x[i] - ta_yy_x_1[i] * pc_x[i];

        ta_xyy_y_0[i] = ta_yy_y_0[i] * pa_x[i] - ta_yy_y_1[i] * pc_x[i];

        ta_xyy_z_0[i] = ta_yy_z_0[i] * pa_x[i] - ta_yy_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto ta_xyz_x_0 = pbuffer.data(idx_npot_0_fp + 12);

    auto ta_xyz_y_0 = pbuffer.data(idx_npot_0_fp + 13);

    auto ta_xyz_z_0 = pbuffer.data(idx_npot_0_fp + 14);

#pragma omp simd aligned( \
        pa_x, pa_y, pc_x, pc_y, ta_xyz_x_0, ta_xyz_y_0, ta_xyz_z_0, ta_xz_x_0, ta_xz_x_1, ta_yz_y_0, ta_yz_y_1, ta_yz_z_0, ta_yz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyz_x_0[i] = ta_xz_x_0[i] * pa_y[i] - ta_xz_x_1[i] * pc_y[i];

        ta_xyz_y_0[i] = ta_yz_y_0[i] * pa_x[i] - ta_yz_y_1[i] * pc_x[i];

        ta_xyz_z_0[i] = ta_yz_z_0[i] * pa_x[i] - ta_yz_z_1[i] * pc_x[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto ta_xzz_x_0 = pbuffer.data(idx_npot_0_fp + 15);

    auto ta_xzz_y_0 = pbuffer.data(idx_npot_0_fp + 16);

    auto ta_xzz_z_0 = pbuffer.data(idx_npot_0_fp + 17);

#pragma omp simd aligned(pa_x,           \
                             pc_x,       \
                             ta_xzz_x_0, \
                             ta_xzz_y_0, \
                             ta_xzz_z_0, \
                             ta_zz_0_0,  \
                             ta_zz_0_1,  \
                             ta_zz_x_0,  \
                             ta_zz_x_1,  \
                             ta_zz_y_0,  \
                             ta_zz_y_1,  \
                             ta_zz_z_0,  \
                             ta_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzz_x_0[i] = ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_zz_x_0[i] * pa_x[i] - ta_zz_x_1[i] * pc_x[i];

        ta_xzz_y_0[i] = ta_zz_y_0[i] * pa_x[i] - ta_zz_y_1[i] * pc_x[i];

        ta_xzz_z_0[i] = ta_zz_z_0[i] * pa_x[i] - ta_zz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto ta_yyy_x_0 = pbuffer.data(idx_npot_0_fp + 18);

    auto ta_yyy_y_0 = pbuffer.data(idx_npot_0_fp + 19);

    auto ta_yyy_z_0 = pbuffer.data(idx_npot_0_fp + 20);

#pragma omp simd aligned(pa_y,           \
                             pc_y,       \
                             ta_y_x_0,   \
                             ta_y_x_1,   \
                             ta_y_y_0,   \
                             ta_y_y_1,   \
                             ta_y_z_0,   \
                             ta_y_z_1,   \
                             ta_yy_0_0,  \
                             ta_yy_0_1,  \
                             ta_yy_x_0,  \
                             ta_yy_x_1,  \
                             ta_yy_y_0,  \
                             ta_yy_y_1,  \
                             ta_yy_z_0,  \
                             ta_yy_z_1,  \
                             ta_yyy_x_0, \
                             ta_yyy_y_0, \
                             ta_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyy_x_0[i] = 2.0 * ta_y_x_0[i] * fe_0 - 2.0 * ta_y_x_1[i] * fe_0 + ta_yy_x_0[i] * pa_y[i] - ta_yy_x_1[i] * pc_y[i];

        ta_yyy_y_0[i] = 2.0 * ta_y_y_0[i] * fe_0 - 2.0 * ta_y_y_1[i] * fe_0 + ta_yy_0_0[i] * fe_0 - ta_yy_0_1[i] * fe_0 + ta_yy_y_0[i] * pa_y[i] -
                        ta_yy_y_1[i] * pc_y[i];

        ta_yyy_z_0[i] = 2.0 * ta_y_z_0[i] * fe_0 - 2.0 * ta_y_z_1[i] * fe_0 + ta_yy_z_0[i] * pa_y[i] - ta_yy_z_1[i] * pc_y[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto ta_yyz_x_0 = pbuffer.data(idx_npot_0_fp + 21);

    auto ta_yyz_y_0 = pbuffer.data(idx_npot_0_fp + 22);

    auto ta_yyz_z_0 = pbuffer.data(idx_npot_0_fp + 23);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             pc_y,       \
                             pc_z,       \
                             ta_yy_x_0,  \
                             ta_yy_x_1,  \
                             ta_yy_y_0,  \
                             ta_yy_y_1,  \
                             ta_yyz_x_0, \
                             ta_yyz_y_0, \
                             ta_yyz_z_0, \
                             ta_yz_z_0,  \
                             ta_yz_z_1,  \
                             ta_z_z_0,   \
                             ta_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyz_x_0[i] = ta_yy_x_0[i] * pa_z[i] - ta_yy_x_1[i] * pc_z[i];

        ta_yyz_y_0[i] = ta_yy_y_0[i] * pa_z[i] - ta_yy_y_1[i] * pc_z[i];

        ta_yyz_z_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_yz_z_0[i] * pa_y[i] - ta_yz_z_1[i] * pc_y[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto ta_yzz_x_0 = pbuffer.data(idx_npot_0_fp + 24);

    auto ta_yzz_y_0 = pbuffer.data(idx_npot_0_fp + 25);

    auto ta_yzz_z_0 = pbuffer.data(idx_npot_0_fp + 26);

#pragma omp simd aligned(pa_y,           \
                             pc_y,       \
                             ta_yzz_x_0, \
                             ta_yzz_y_0, \
                             ta_yzz_z_0, \
                             ta_zz_0_0,  \
                             ta_zz_0_1,  \
                             ta_zz_x_0,  \
                             ta_zz_x_1,  \
                             ta_zz_y_0,  \
                             ta_zz_y_1,  \
                             ta_zz_z_0,  \
                             ta_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzz_x_0[i] = ta_zz_x_0[i] * pa_y[i] - ta_zz_x_1[i] * pc_y[i];

        ta_yzz_y_0[i] = ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_zz_y_0[i] * pa_y[i] - ta_zz_y_1[i] * pc_y[i];

        ta_yzz_z_0[i] = ta_zz_z_0[i] * pa_y[i] - ta_zz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto ta_zzz_x_0 = pbuffer.data(idx_npot_0_fp + 27);

    auto ta_zzz_y_0 = pbuffer.data(idx_npot_0_fp + 28);

    auto ta_zzz_z_0 = pbuffer.data(idx_npot_0_fp + 29);

#pragma omp simd aligned(pa_z,           \
                             pc_z,       \
                             ta_z_x_0,   \
                             ta_z_x_1,   \
                             ta_z_y_0,   \
                             ta_z_y_1,   \
                             ta_z_z_0,   \
                             ta_z_z_1,   \
                             ta_zz_0_0,  \
                             ta_zz_0_1,  \
                             ta_zz_x_0,  \
                             ta_zz_x_1,  \
                             ta_zz_y_0,  \
                             ta_zz_y_1,  \
                             ta_zz_z_0,  \
                             ta_zz_z_1,  \
                             ta_zzz_x_0, \
                             ta_zzz_y_0, \
                             ta_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzz_x_0[i] = 2.0 * ta_z_x_0[i] * fe_0 - 2.0 * ta_z_x_1[i] * fe_0 + ta_zz_x_0[i] * pa_z[i] - ta_zz_x_1[i] * pc_z[i];

        ta_zzz_y_0[i] = 2.0 * ta_z_y_0[i] * fe_0 - 2.0 * ta_z_y_1[i] * fe_0 + ta_zz_y_0[i] * pa_z[i] - ta_zz_y_1[i] * pc_z[i];

        ta_zzz_z_0[i] = 2.0 * ta_z_z_0[i] * fe_0 - 2.0 * ta_z_z_1[i] * fe_0 + ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_zz_z_0[i] * pa_z[i] -
                        ta_zz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
