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

#include "NuclearPotentialPrimRecPF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_pf(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_pf,
                               const size_t              idx_npot_0_sd,
                               const size_t              idx_npot_1_sd,
                               const size_t              idx_npot_0_sf,
                               const size_t              idx_npot_1_sf,
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

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_0 = pbuffer.data(idx_npot_0_sd);

    auto ta_0_xy_0 = pbuffer.data(idx_npot_0_sd + 1);

    auto ta_0_xz_0 = pbuffer.data(idx_npot_0_sd + 2);

    auto ta_0_yy_0 = pbuffer.data(idx_npot_0_sd + 3);

    auto ta_0_yz_0 = pbuffer.data(idx_npot_0_sd + 4);

    auto ta_0_zz_0 = pbuffer.data(idx_npot_0_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = pbuffer.data(idx_npot_1_sd);

    auto ta_0_xy_1 = pbuffer.data(idx_npot_1_sd + 1);

    auto ta_0_xz_1 = pbuffer.data(idx_npot_1_sd + 2);

    auto ta_0_yy_1 = pbuffer.data(idx_npot_1_sd + 3);

    auto ta_0_yz_1 = pbuffer.data(idx_npot_1_sd + 4);

    auto ta_0_zz_1 = pbuffer.data(idx_npot_1_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_0 = pbuffer.data(idx_npot_0_sf);

    auto ta_0_xxy_0 = pbuffer.data(idx_npot_0_sf + 1);

    auto ta_0_xxz_0 = pbuffer.data(idx_npot_0_sf + 2);

    auto ta_0_xyy_0 = pbuffer.data(idx_npot_0_sf + 3);

    auto ta_0_xyz_0 = pbuffer.data(idx_npot_0_sf + 4);

    auto ta_0_xzz_0 = pbuffer.data(idx_npot_0_sf + 5);

    auto ta_0_yyy_0 = pbuffer.data(idx_npot_0_sf + 6);

    auto ta_0_yyz_0 = pbuffer.data(idx_npot_0_sf + 7);

    auto ta_0_yzz_0 = pbuffer.data(idx_npot_0_sf + 8);

    auto ta_0_zzz_0 = pbuffer.data(idx_npot_0_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_1 = pbuffer.data(idx_npot_1_sf);

    auto ta_0_xxy_1 = pbuffer.data(idx_npot_1_sf + 1);

    auto ta_0_xxz_1 = pbuffer.data(idx_npot_1_sf + 2);

    auto ta_0_xyy_1 = pbuffer.data(idx_npot_1_sf + 3);

    auto ta_0_xyz_1 = pbuffer.data(idx_npot_1_sf + 4);

    auto ta_0_xzz_1 = pbuffer.data(idx_npot_1_sf + 5);

    auto ta_0_yyy_1 = pbuffer.data(idx_npot_1_sf + 6);

    auto ta_0_yyz_1 = pbuffer.data(idx_npot_1_sf + 7);

    auto ta_0_yzz_1 = pbuffer.data(idx_npot_1_sf + 8);

    auto ta_0_zzz_1 = pbuffer.data(idx_npot_1_sf + 9);

    // Set up 0-10 components of targeted buffer : PF

    auto ta_x_xxx_0 = pbuffer.data(idx_npot_0_pf);

    auto ta_x_xxy_0 = pbuffer.data(idx_npot_0_pf + 1);

    auto ta_x_xxz_0 = pbuffer.data(idx_npot_0_pf + 2);

    auto ta_x_xyy_0 = pbuffer.data(idx_npot_0_pf + 3);

    auto ta_x_xyz_0 = pbuffer.data(idx_npot_0_pf + 4);

    auto ta_x_xzz_0 = pbuffer.data(idx_npot_0_pf + 5);

    auto ta_x_yyy_0 = pbuffer.data(idx_npot_0_pf + 6);

    auto ta_x_yyz_0 = pbuffer.data(idx_npot_0_pf + 7);

    auto ta_x_yzz_0 = pbuffer.data(idx_npot_0_pf + 8);

    auto ta_x_zzz_0 = pbuffer.data(idx_npot_0_pf + 9);

#pragma omp simd aligned(pa_x,           \
                             pc_x,       \
                             ta_0_xx_0,  \
                             ta_0_xx_1,  \
                             ta_0_xxx_0, \
                             ta_0_xxx_1, \
                             ta_0_xxy_0, \
                             ta_0_xxy_1, \
                             ta_0_xxz_0, \
                             ta_0_xxz_1, \
                             ta_0_xy_0,  \
                             ta_0_xy_1,  \
                             ta_0_xyy_0, \
                             ta_0_xyy_1, \
                             ta_0_xyz_0, \
                             ta_0_xyz_1, \
                             ta_0_xz_0,  \
                             ta_0_xz_1,  \
                             ta_0_xzz_0, \
                             ta_0_xzz_1, \
                             ta_0_yy_0,  \
                             ta_0_yy_1,  \
                             ta_0_yyy_0, \
                             ta_0_yyy_1, \
                             ta_0_yyz_0, \
                             ta_0_yyz_1, \
                             ta_0_yz_0,  \
                             ta_0_yz_1,  \
                             ta_0_yzz_0, \
                             ta_0_yzz_1, \
                             ta_0_zz_0,  \
                             ta_0_zz_1,  \
                             ta_0_zzz_0, \
                             ta_0_zzz_1, \
                             ta_x_xxx_0, \
                             ta_x_xxy_0, \
                             ta_x_xxz_0, \
                             ta_x_xyy_0, \
                             ta_x_xyz_0, \
                             ta_x_xzz_0, \
                             ta_x_yyy_0, \
                             ta_x_yyz_0, \
                             ta_x_yzz_0, \
                             ta_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_x_xxx_0[i] = 3.0 * ta_0_xx_0[i] * fe_0 - 3.0 * ta_0_xx_1[i] * fe_0 + ta_0_xxx_0[i] * pa_x[i] - ta_0_xxx_1[i] * pc_x[i];

        ta_x_xxy_0[i] = 2.0 * ta_0_xy_0[i] * fe_0 - 2.0 * ta_0_xy_1[i] * fe_0 + ta_0_xxy_0[i] * pa_x[i] - ta_0_xxy_1[i] * pc_x[i];

        ta_x_xxz_0[i] = 2.0 * ta_0_xz_0[i] * fe_0 - 2.0 * ta_0_xz_1[i] * fe_0 + ta_0_xxz_0[i] * pa_x[i] - ta_0_xxz_1[i] * pc_x[i];

        ta_x_xyy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_0_xyy_0[i] * pa_x[i] - ta_0_xyy_1[i] * pc_x[i];

        ta_x_xyz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_0_xyz_0[i] * pa_x[i] - ta_0_xyz_1[i] * pc_x[i];

        ta_x_xzz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_0_xzz_0[i] * pa_x[i] - ta_0_xzz_1[i] * pc_x[i];

        ta_x_yyy_0[i] = ta_0_yyy_0[i] * pa_x[i] - ta_0_yyy_1[i] * pc_x[i];

        ta_x_yyz_0[i] = ta_0_yyz_0[i] * pa_x[i] - ta_0_yyz_1[i] * pc_x[i];

        ta_x_yzz_0[i] = ta_0_yzz_0[i] * pa_x[i] - ta_0_yzz_1[i] * pc_x[i];

        ta_x_zzz_0[i] = ta_0_zzz_0[i] * pa_x[i] - ta_0_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto ta_y_xxx_0 = pbuffer.data(idx_npot_0_pf + 10);

    auto ta_y_xxy_0 = pbuffer.data(idx_npot_0_pf + 11);

    auto ta_y_xxz_0 = pbuffer.data(idx_npot_0_pf + 12);

    auto ta_y_xyy_0 = pbuffer.data(idx_npot_0_pf + 13);

    auto ta_y_xyz_0 = pbuffer.data(idx_npot_0_pf + 14);

    auto ta_y_xzz_0 = pbuffer.data(idx_npot_0_pf + 15);

    auto ta_y_yyy_0 = pbuffer.data(idx_npot_0_pf + 16);

    auto ta_y_yyz_0 = pbuffer.data(idx_npot_0_pf + 17);

    auto ta_y_yzz_0 = pbuffer.data(idx_npot_0_pf + 18);

    auto ta_y_zzz_0 = pbuffer.data(idx_npot_0_pf + 19);

#pragma omp simd aligned(pa_y,           \
                             pc_y,       \
                             ta_0_xx_0,  \
                             ta_0_xx_1,  \
                             ta_0_xxx_0, \
                             ta_0_xxx_1, \
                             ta_0_xxy_0, \
                             ta_0_xxy_1, \
                             ta_0_xxz_0, \
                             ta_0_xxz_1, \
                             ta_0_xy_0,  \
                             ta_0_xy_1,  \
                             ta_0_xyy_0, \
                             ta_0_xyy_1, \
                             ta_0_xyz_0, \
                             ta_0_xyz_1, \
                             ta_0_xz_0,  \
                             ta_0_xz_1,  \
                             ta_0_xzz_0, \
                             ta_0_xzz_1, \
                             ta_0_yy_0,  \
                             ta_0_yy_1,  \
                             ta_0_yyy_0, \
                             ta_0_yyy_1, \
                             ta_0_yyz_0, \
                             ta_0_yyz_1, \
                             ta_0_yz_0,  \
                             ta_0_yz_1,  \
                             ta_0_yzz_0, \
                             ta_0_yzz_1, \
                             ta_0_zz_0,  \
                             ta_0_zz_1,  \
                             ta_0_zzz_0, \
                             ta_0_zzz_1, \
                             ta_y_xxx_0, \
                             ta_y_xxy_0, \
                             ta_y_xxz_0, \
                             ta_y_xyy_0, \
                             ta_y_xyz_0, \
                             ta_y_xzz_0, \
                             ta_y_yyy_0, \
                             ta_y_yyz_0, \
                             ta_y_yzz_0, \
                             ta_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_y_xxx_0[i] = ta_0_xxx_0[i] * pa_y[i] - ta_0_xxx_1[i] * pc_y[i];

        ta_y_xxy_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_0_xxy_0[i] * pa_y[i] - ta_0_xxy_1[i] * pc_y[i];

        ta_y_xxz_0[i] = ta_0_xxz_0[i] * pa_y[i] - ta_0_xxz_1[i] * pc_y[i];

        ta_y_xyy_0[i] = 2.0 * ta_0_xy_0[i] * fe_0 - 2.0 * ta_0_xy_1[i] * fe_0 + ta_0_xyy_0[i] * pa_y[i] - ta_0_xyy_1[i] * pc_y[i];

        ta_y_xyz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_0_xyz_0[i] * pa_y[i] - ta_0_xyz_1[i] * pc_y[i];

        ta_y_xzz_0[i] = ta_0_xzz_0[i] * pa_y[i] - ta_0_xzz_1[i] * pc_y[i];

        ta_y_yyy_0[i] = 3.0 * ta_0_yy_0[i] * fe_0 - 3.0 * ta_0_yy_1[i] * fe_0 + ta_0_yyy_0[i] * pa_y[i] - ta_0_yyy_1[i] * pc_y[i];

        ta_y_yyz_0[i] = 2.0 * ta_0_yz_0[i] * fe_0 - 2.0 * ta_0_yz_1[i] * fe_0 + ta_0_yyz_0[i] * pa_y[i] - ta_0_yyz_1[i] * pc_y[i];

        ta_y_yzz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_0_yzz_0[i] * pa_y[i] - ta_0_yzz_1[i] * pc_y[i];

        ta_y_zzz_0[i] = ta_0_zzz_0[i] * pa_y[i] - ta_0_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : PF

    auto ta_z_xxx_0 = pbuffer.data(idx_npot_0_pf + 20);

    auto ta_z_xxy_0 = pbuffer.data(idx_npot_0_pf + 21);

    auto ta_z_xxz_0 = pbuffer.data(idx_npot_0_pf + 22);

    auto ta_z_xyy_0 = pbuffer.data(idx_npot_0_pf + 23);

    auto ta_z_xyz_0 = pbuffer.data(idx_npot_0_pf + 24);

    auto ta_z_xzz_0 = pbuffer.data(idx_npot_0_pf + 25);

    auto ta_z_yyy_0 = pbuffer.data(idx_npot_0_pf + 26);

    auto ta_z_yyz_0 = pbuffer.data(idx_npot_0_pf + 27);

    auto ta_z_yzz_0 = pbuffer.data(idx_npot_0_pf + 28);

    auto ta_z_zzz_0 = pbuffer.data(idx_npot_0_pf + 29);

#pragma omp simd aligned(pa_z,           \
                             pc_z,       \
                             ta_0_xx_0,  \
                             ta_0_xx_1,  \
                             ta_0_xxx_0, \
                             ta_0_xxx_1, \
                             ta_0_xxy_0, \
                             ta_0_xxy_1, \
                             ta_0_xxz_0, \
                             ta_0_xxz_1, \
                             ta_0_xy_0,  \
                             ta_0_xy_1,  \
                             ta_0_xyy_0, \
                             ta_0_xyy_1, \
                             ta_0_xyz_0, \
                             ta_0_xyz_1, \
                             ta_0_xz_0,  \
                             ta_0_xz_1,  \
                             ta_0_xzz_0, \
                             ta_0_xzz_1, \
                             ta_0_yy_0,  \
                             ta_0_yy_1,  \
                             ta_0_yyy_0, \
                             ta_0_yyy_1, \
                             ta_0_yyz_0, \
                             ta_0_yyz_1, \
                             ta_0_yz_0,  \
                             ta_0_yz_1,  \
                             ta_0_yzz_0, \
                             ta_0_yzz_1, \
                             ta_0_zz_0,  \
                             ta_0_zz_1,  \
                             ta_0_zzz_0, \
                             ta_0_zzz_1, \
                             ta_z_xxx_0, \
                             ta_z_xxy_0, \
                             ta_z_xxz_0, \
                             ta_z_xyy_0, \
                             ta_z_xyz_0, \
                             ta_z_xzz_0, \
                             ta_z_yyy_0, \
                             ta_z_yyz_0, \
                             ta_z_yzz_0, \
                             ta_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_z_xxx_0[i] = ta_0_xxx_0[i] * pa_z[i] - ta_0_xxx_1[i] * pc_z[i];

        ta_z_xxy_0[i] = ta_0_xxy_0[i] * pa_z[i] - ta_0_xxy_1[i] * pc_z[i];

        ta_z_xxz_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_0_xxz_0[i] * pa_z[i] - ta_0_xxz_1[i] * pc_z[i];

        ta_z_xyy_0[i] = ta_0_xyy_0[i] * pa_z[i] - ta_0_xyy_1[i] * pc_z[i];

        ta_z_xyz_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_0_xyz_0[i] * pa_z[i] - ta_0_xyz_1[i] * pc_z[i];

        ta_z_xzz_0[i] = 2.0 * ta_0_xz_0[i] * fe_0 - 2.0 * ta_0_xz_1[i] * fe_0 + ta_0_xzz_0[i] * pa_z[i] - ta_0_xzz_1[i] * pc_z[i];

        ta_z_yyy_0[i] = ta_0_yyy_0[i] * pa_z[i] - ta_0_yyy_1[i] * pc_z[i];

        ta_z_yyz_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_0_yyz_0[i] * pa_z[i] - ta_0_yyz_1[i] * pc_z[i];

        ta_z_yzz_0[i] = 2.0 * ta_0_yz_0[i] * fe_0 - 2.0 * ta_0_yz_1[i] * fe_0 + ta_0_yzz_0[i] * pa_z[i] - ta_0_yzz_1[i] * pc_z[i];

        ta_z_zzz_0[i] = 3.0 * ta_0_zz_0[i] * fe_0 - 3.0 * ta_0_zz_1[i] * fe_0 + ta_0_zzz_0[i] * pa_z[i] - ta_0_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
