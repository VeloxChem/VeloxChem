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

#include "NuclearPotentialPrimRecFF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_ff(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_ff,
                               const size_t              idx_npot_0_pf,
                               const size_t              idx_npot_1_pf,
                               const size_t              idx_npot_0_dd,
                               const size_t              idx_npot_1_dd,
                               const size_t              idx_npot_0_df,
                               const size_t              idx_npot_1_df,
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

    // Set up components of auxiliary buffer : PF

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

    // Set up components of auxiliary buffer : PF

    auto ta_x_xxx_1 = pbuffer.data(idx_npot_1_pf);

    auto ta_x_xxy_1 = pbuffer.data(idx_npot_1_pf + 1);

    auto ta_x_xxz_1 = pbuffer.data(idx_npot_1_pf + 2);

    auto ta_x_xyy_1 = pbuffer.data(idx_npot_1_pf + 3);

    auto ta_x_xyz_1 = pbuffer.data(idx_npot_1_pf + 4);

    auto ta_x_xzz_1 = pbuffer.data(idx_npot_1_pf + 5);

    auto ta_x_yyy_1 = pbuffer.data(idx_npot_1_pf + 6);

    auto ta_x_yyz_1 = pbuffer.data(idx_npot_1_pf + 7);

    auto ta_x_yzz_1 = pbuffer.data(idx_npot_1_pf + 8);

    auto ta_x_zzz_1 = pbuffer.data(idx_npot_1_pf + 9);

    auto ta_y_xxx_1 = pbuffer.data(idx_npot_1_pf + 10);

    auto ta_y_xxy_1 = pbuffer.data(idx_npot_1_pf + 11);

    auto ta_y_xxz_1 = pbuffer.data(idx_npot_1_pf + 12);

    auto ta_y_xyy_1 = pbuffer.data(idx_npot_1_pf + 13);

    auto ta_y_xyz_1 = pbuffer.data(idx_npot_1_pf + 14);

    auto ta_y_xzz_1 = pbuffer.data(idx_npot_1_pf + 15);

    auto ta_y_yyy_1 = pbuffer.data(idx_npot_1_pf + 16);

    auto ta_y_yyz_1 = pbuffer.data(idx_npot_1_pf + 17);

    auto ta_y_yzz_1 = pbuffer.data(idx_npot_1_pf + 18);

    auto ta_y_zzz_1 = pbuffer.data(idx_npot_1_pf + 19);

    auto ta_z_xxx_1 = pbuffer.data(idx_npot_1_pf + 20);

    auto ta_z_xxy_1 = pbuffer.data(idx_npot_1_pf + 21);

    auto ta_z_xxz_1 = pbuffer.data(idx_npot_1_pf + 22);

    auto ta_z_xyy_1 = pbuffer.data(idx_npot_1_pf + 23);

    auto ta_z_xyz_1 = pbuffer.data(idx_npot_1_pf + 24);

    auto ta_z_xzz_1 = pbuffer.data(idx_npot_1_pf + 25);

    auto ta_z_yyy_1 = pbuffer.data(idx_npot_1_pf + 26);

    auto ta_z_yyz_1 = pbuffer.data(idx_npot_1_pf + 27);

    auto ta_z_yzz_1 = pbuffer.data(idx_npot_1_pf + 28);

    auto ta_z_zzz_1 = pbuffer.data(idx_npot_1_pf + 29);

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_0 = pbuffer.data(idx_npot_0_dd);

    auto ta_xx_xy_0 = pbuffer.data(idx_npot_0_dd + 1);

    auto ta_xx_xz_0 = pbuffer.data(idx_npot_0_dd + 2);

    auto ta_xx_yy_0 = pbuffer.data(idx_npot_0_dd + 3);

    auto ta_xx_yz_0 = pbuffer.data(idx_npot_0_dd + 4);

    auto ta_xx_zz_0 = pbuffer.data(idx_npot_0_dd + 5);

    auto ta_yy_xx_0 = pbuffer.data(idx_npot_0_dd + 18);

    auto ta_yy_xy_0 = pbuffer.data(idx_npot_0_dd + 19);

    auto ta_yy_xz_0 = pbuffer.data(idx_npot_0_dd + 20);

    auto ta_yy_yy_0 = pbuffer.data(idx_npot_0_dd + 21);

    auto ta_yy_yz_0 = pbuffer.data(idx_npot_0_dd + 22);

    auto ta_yy_zz_0 = pbuffer.data(idx_npot_0_dd + 23);

    auto ta_yz_yz_0 = pbuffer.data(idx_npot_0_dd + 28);

    auto ta_zz_xx_0 = pbuffer.data(idx_npot_0_dd + 30);

    auto ta_zz_xy_0 = pbuffer.data(idx_npot_0_dd + 31);

    auto ta_zz_xz_0 = pbuffer.data(idx_npot_0_dd + 32);

    auto ta_zz_yy_0 = pbuffer.data(idx_npot_0_dd + 33);

    auto ta_zz_yz_0 = pbuffer.data(idx_npot_0_dd + 34);

    auto ta_zz_zz_0 = pbuffer.data(idx_npot_0_dd + 35);

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_1 = pbuffer.data(idx_npot_1_dd);

    auto ta_xx_xy_1 = pbuffer.data(idx_npot_1_dd + 1);

    auto ta_xx_xz_1 = pbuffer.data(idx_npot_1_dd + 2);

    auto ta_xx_yy_1 = pbuffer.data(idx_npot_1_dd + 3);

    auto ta_xx_yz_1 = pbuffer.data(idx_npot_1_dd + 4);

    auto ta_xx_zz_1 = pbuffer.data(idx_npot_1_dd + 5);

    auto ta_yy_xx_1 = pbuffer.data(idx_npot_1_dd + 18);

    auto ta_yy_xy_1 = pbuffer.data(idx_npot_1_dd + 19);

    auto ta_yy_xz_1 = pbuffer.data(idx_npot_1_dd + 20);

    auto ta_yy_yy_1 = pbuffer.data(idx_npot_1_dd + 21);

    auto ta_yy_yz_1 = pbuffer.data(idx_npot_1_dd + 22);

    auto ta_yy_zz_1 = pbuffer.data(idx_npot_1_dd + 23);

    auto ta_yz_yz_1 = pbuffer.data(idx_npot_1_dd + 28);

    auto ta_zz_xx_1 = pbuffer.data(idx_npot_1_dd + 30);

    auto ta_zz_xy_1 = pbuffer.data(idx_npot_1_dd + 31);

    auto ta_zz_xz_1 = pbuffer.data(idx_npot_1_dd + 32);

    auto ta_zz_yy_1 = pbuffer.data(idx_npot_1_dd + 33);

    auto ta_zz_yz_1 = pbuffer.data(idx_npot_1_dd + 34);

    auto ta_zz_zz_1 = pbuffer.data(idx_npot_1_dd + 35);

    // Set up components of auxiliary buffer : DF

    auto ta_xx_xxx_0 = pbuffer.data(idx_npot_0_df);

    auto ta_xx_xxy_0 = pbuffer.data(idx_npot_0_df + 1);

    auto ta_xx_xxz_0 = pbuffer.data(idx_npot_0_df + 2);

    auto ta_xx_xyy_0 = pbuffer.data(idx_npot_0_df + 3);

    auto ta_xx_xyz_0 = pbuffer.data(idx_npot_0_df + 4);

    auto ta_xx_xzz_0 = pbuffer.data(idx_npot_0_df + 5);

    auto ta_xx_yyy_0 = pbuffer.data(idx_npot_0_df + 6);

    auto ta_xx_yyz_0 = pbuffer.data(idx_npot_0_df + 7);

    auto ta_xx_yzz_0 = pbuffer.data(idx_npot_0_df + 8);

    auto ta_xx_zzz_0 = pbuffer.data(idx_npot_0_df + 9);

    auto ta_xy_xxy_0 = pbuffer.data(idx_npot_0_df + 11);

    auto ta_xy_xyy_0 = pbuffer.data(idx_npot_0_df + 13);

    auto ta_xy_yyy_0 = pbuffer.data(idx_npot_0_df + 16);

    auto ta_xy_yyz_0 = pbuffer.data(idx_npot_0_df + 17);

    auto ta_xy_yzz_0 = pbuffer.data(idx_npot_0_df + 18);

    auto ta_xz_xxx_0 = pbuffer.data(idx_npot_0_df + 20);

    auto ta_xz_xxz_0 = pbuffer.data(idx_npot_0_df + 22);

    auto ta_xz_xzz_0 = pbuffer.data(idx_npot_0_df + 25);

    auto ta_xz_yyz_0 = pbuffer.data(idx_npot_0_df + 27);

    auto ta_xz_yzz_0 = pbuffer.data(idx_npot_0_df + 28);

    auto ta_xz_zzz_0 = pbuffer.data(idx_npot_0_df + 29);

    auto ta_yy_xxx_0 = pbuffer.data(idx_npot_0_df + 30);

    auto ta_yy_xxy_0 = pbuffer.data(idx_npot_0_df + 31);

    auto ta_yy_xxz_0 = pbuffer.data(idx_npot_0_df + 32);

    auto ta_yy_xyy_0 = pbuffer.data(idx_npot_0_df + 33);

    auto ta_yy_xyz_0 = pbuffer.data(idx_npot_0_df + 34);

    auto ta_yy_xzz_0 = pbuffer.data(idx_npot_0_df + 35);

    auto ta_yy_yyy_0 = pbuffer.data(idx_npot_0_df + 36);

    auto ta_yy_yyz_0 = pbuffer.data(idx_npot_0_df + 37);

    auto ta_yy_yzz_0 = pbuffer.data(idx_npot_0_df + 38);

    auto ta_yy_zzz_0 = pbuffer.data(idx_npot_0_df + 39);

    auto ta_yz_xxz_0 = pbuffer.data(idx_npot_0_df + 42);

    auto ta_yz_xyz_0 = pbuffer.data(idx_npot_0_df + 44);

    auto ta_yz_xzz_0 = pbuffer.data(idx_npot_0_df + 45);

    auto ta_yz_yyy_0 = pbuffer.data(idx_npot_0_df + 46);

    auto ta_yz_yyz_0 = pbuffer.data(idx_npot_0_df + 47);

    auto ta_yz_yzz_0 = pbuffer.data(idx_npot_0_df + 48);

    auto ta_yz_zzz_0 = pbuffer.data(idx_npot_0_df + 49);

    auto ta_zz_xxx_0 = pbuffer.data(idx_npot_0_df + 50);

    auto ta_zz_xxy_0 = pbuffer.data(idx_npot_0_df + 51);

    auto ta_zz_xxz_0 = pbuffer.data(idx_npot_0_df + 52);

    auto ta_zz_xyy_0 = pbuffer.data(idx_npot_0_df + 53);

    auto ta_zz_xyz_0 = pbuffer.data(idx_npot_0_df + 54);

    auto ta_zz_xzz_0 = pbuffer.data(idx_npot_0_df + 55);

    auto ta_zz_yyy_0 = pbuffer.data(idx_npot_0_df + 56);

    auto ta_zz_yyz_0 = pbuffer.data(idx_npot_0_df + 57);

    auto ta_zz_yzz_0 = pbuffer.data(idx_npot_0_df + 58);

    auto ta_zz_zzz_0 = pbuffer.data(idx_npot_0_df + 59);

    // Set up components of auxiliary buffer : DF

    auto ta_xx_xxx_1 = pbuffer.data(idx_npot_1_df);

    auto ta_xx_xxy_1 = pbuffer.data(idx_npot_1_df + 1);

    auto ta_xx_xxz_1 = pbuffer.data(idx_npot_1_df + 2);

    auto ta_xx_xyy_1 = pbuffer.data(idx_npot_1_df + 3);

    auto ta_xx_xyz_1 = pbuffer.data(idx_npot_1_df + 4);

    auto ta_xx_xzz_1 = pbuffer.data(idx_npot_1_df + 5);

    auto ta_xx_yyy_1 = pbuffer.data(idx_npot_1_df + 6);

    auto ta_xx_yyz_1 = pbuffer.data(idx_npot_1_df + 7);

    auto ta_xx_yzz_1 = pbuffer.data(idx_npot_1_df + 8);

    auto ta_xx_zzz_1 = pbuffer.data(idx_npot_1_df + 9);

    auto ta_xy_xxy_1 = pbuffer.data(idx_npot_1_df + 11);

    auto ta_xy_xyy_1 = pbuffer.data(idx_npot_1_df + 13);

    auto ta_xy_yyy_1 = pbuffer.data(idx_npot_1_df + 16);

    auto ta_xy_yyz_1 = pbuffer.data(idx_npot_1_df + 17);

    auto ta_xy_yzz_1 = pbuffer.data(idx_npot_1_df + 18);

    auto ta_xz_xxx_1 = pbuffer.data(idx_npot_1_df + 20);

    auto ta_xz_xxz_1 = pbuffer.data(idx_npot_1_df + 22);

    auto ta_xz_xzz_1 = pbuffer.data(idx_npot_1_df + 25);

    auto ta_xz_yyz_1 = pbuffer.data(idx_npot_1_df + 27);

    auto ta_xz_yzz_1 = pbuffer.data(idx_npot_1_df + 28);

    auto ta_xz_zzz_1 = pbuffer.data(idx_npot_1_df + 29);

    auto ta_yy_xxx_1 = pbuffer.data(idx_npot_1_df + 30);

    auto ta_yy_xxy_1 = pbuffer.data(idx_npot_1_df + 31);

    auto ta_yy_xxz_1 = pbuffer.data(idx_npot_1_df + 32);

    auto ta_yy_xyy_1 = pbuffer.data(idx_npot_1_df + 33);

    auto ta_yy_xyz_1 = pbuffer.data(idx_npot_1_df + 34);

    auto ta_yy_xzz_1 = pbuffer.data(idx_npot_1_df + 35);

    auto ta_yy_yyy_1 = pbuffer.data(idx_npot_1_df + 36);

    auto ta_yy_yyz_1 = pbuffer.data(idx_npot_1_df + 37);

    auto ta_yy_yzz_1 = pbuffer.data(idx_npot_1_df + 38);

    auto ta_yy_zzz_1 = pbuffer.data(idx_npot_1_df + 39);

    auto ta_yz_xxz_1 = pbuffer.data(idx_npot_1_df + 42);

    auto ta_yz_xyz_1 = pbuffer.data(idx_npot_1_df + 44);

    auto ta_yz_xzz_1 = pbuffer.data(idx_npot_1_df + 45);

    auto ta_yz_yyy_1 = pbuffer.data(idx_npot_1_df + 46);

    auto ta_yz_yyz_1 = pbuffer.data(idx_npot_1_df + 47);

    auto ta_yz_yzz_1 = pbuffer.data(idx_npot_1_df + 48);

    auto ta_yz_zzz_1 = pbuffer.data(idx_npot_1_df + 49);

    auto ta_zz_xxx_1 = pbuffer.data(idx_npot_1_df + 50);

    auto ta_zz_xxy_1 = pbuffer.data(idx_npot_1_df + 51);

    auto ta_zz_xxz_1 = pbuffer.data(idx_npot_1_df + 52);

    auto ta_zz_xyy_1 = pbuffer.data(idx_npot_1_df + 53);

    auto ta_zz_xyz_1 = pbuffer.data(idx_npot_1_df + 54);

    auto ta_zz_xzz_1 = pbuffer.data(idx_npot_1_df + 55);

    auto ta_zz_yyy_1 = pbuffer.data(idx_npot_1_df + 56);

    auto ta_zz_yyz_1 = pbuffer.data(idx_npot_1_df + 57);

    auto ta_zz_yzz_1 = pbuffer.data(idx_npot_1_df + 58);

    auto ta_zz_zzz_1 = pbuffer.data(idx_npot_1_df + 59);

    // Set up 0-10 components of targeted buffer : FF

    auto ta_xxx_xxx_0 = pbuffer.data(idx_npot_0_ff);

    auto ta_xxx_xxy_0 = pbuffer.data(idx_npot_0_ff + 1);

    auto ta_xxx_xxz_0 = pbuffer.data(idx_npot_0_ff + 2);

    auto ta_xxx_xyy_0 = pbuffer.data(idx_npot_0_ff + 3);

    auto ta_xxx_xyz_0 = pbuffer.data(idx_npot_0_ff + 4);

    auto ta_xxx_xzz_0 = pbuffer.data(idx_npot_0_ff + 5);

    auto ta_xxx_yyy_0 = pbuffer.data(idx_npot_0_ff + 6);

    auto ta_xxx_yyz_0 = pbuffer.data(idx_npot_0_ff + 7);

    auto ta_xxx_yzz_0 = pbuffer.data(idx_npot_0_ff + 8);

    auto ta_xxx_zzz_0 = pbuffer.data(idx_npot_0_ff + 9);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_x_xxx_0,   \
                             ta_x_xxx_1,   \
                             ta_x_xxy_0,   \
                             ta_x_xxy_1,   \
                             ta_x_xxz_0,   \
                             ta_x_xxz_1,   \
                             ta_x_xyy_0,   \
                             ta_x_xyy_1,   \
                             ta_x_xyz_0,   \
                             ta_x_xyz_1,   \
                             ta_x_xzz_0,   \
                             ta_x_xzz_1,   \
                             ta_x_yyy_0,   \
                             ta_x_yyy_1,   \
                             ta_x_yyz_0,   \
                             ta_x_yyz_1,   \
                             ta_x_yzz_0,   \
                             ta_x_yzz_1,   \
                             ta_x_zzz_0,   \
                             ta_x_zzz_1,   \
                             ta_xx_xx_0,   \
                             ta_xx_xx_1,   \
                             ta_xx_xxx_0,  \
                             ta_xx_xxx_1,  \
                             ta_xx_xxy_0,  \
                             ta_xx_xxy_1,  \
                             ta_xx_xxz_0,  \
                             ta_xx_xxz_1,  \
                             ta_xx_xy_0,   \
                             ta_xx_xy_1,   \
                             ta_xx_xyy_0,  \
                             ta_xx_xyy_1,  \
                             ta_xx_xyz_0,  \
                             ta_xx_xyz_1,  \
                             ta_xx_xz_0,   \
                             ta_xx_xz_1,   \
                             ta_xx_xzz_0,  \
                             ta_xx_xzz_1,  \
                             ta_xx_yy_0,   \
                             ta_xx_yy_1,   \
                             ta_xx_yyy_0,  \
                             ta_xx_yyy_1,  \
                             ta_xx_yyz_0,  \
                             ta_xx_yyz_1,  \
                             ta_xx_yz_0,   \
                             ta_xx_yz_1,   \
                             ta_xx_yzz_0,  \
                             ta_xx_yzz_1,  \
                             ta_xx_zz_0,   \
                             ta_xx_zz_1,   \
                             ta_xx_zzz_0,  \
                             ta_xx_zzz_1,  \
                             ta_xxx_xxx_0, \
                             ta_xxx_xxy_0, \
                             ta_xxx_xxz_0, \
                             ta_xxx_xyy_0, \
                             ta_xxx_xyz_0, \
                             ta_xxx_xzz_0, \
                             ta_xxx_yyy_0, \
                             ta_xxx_yyz_0, \
                             ta_xxx_yzz_0, \
                             ta_xxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxx_xxx_0[i] = 2.0 * ta_x_xxx_0[i] * fe_0 - 2.0 * ta_x_xxx_1[i] * fe_0 + 3.0 * ta_xx_xx_0[i] * fe_0 - 3.0 * ta_xx_xx_1[i] * fe_0 +
                          ta_xx_xxx_0[i] * pa_x[i] - ta_xx_xxx_1[i] * pc_x[i];

        ta_xxx_xxy_0[i] = 2.0 * ta_x_xxy_0[i] * fe_0 - 2.0 * ta_x_xxy_1[i] * fe_0 + 2.0 * ta_xx_xy_0[i] * fe_0 - 2.0 * ta_xx_xy_1[i] * fe_0 +
                          ta_xx_xxy_0[i] * pa_x[i] - ta_xx_xxy_1[i] * pc_x[i];

        ta_xxx_xxz_0[i] = 2.0 * ta_x_xxz_0[i] * fe_0 - 2.0 * ta_x_xxz_1[i] * fe_0 + 2.0 * ta_xx_xz_0[i] * fe_0 - 2.0 * ta_xx_xz_1[i] * fe_0 +
                          ta_xx_xxz_0[i] * pa_x[i] - ta_xx_xxz_1[i] * pc_x[i];

        ta_xxx_xyy_0[i] = 2.0 * ta_x_xyy_0[i] * fe_0 - 2.0 * ta_x_xyy_1[i] * fe_0 + ta_xx_yy_0[i] * fe_0 - ta_xx_yy_1[i] * fe_0 +
                          ta_xx_xyy_0[i] * pa_x[i] - ta_xx_xyy_1[i] * pc_x[i];

        ta_xxx_xyz_0[i] = 2.0 * ta_x_xyz_0[i] * fe_0 - 2.0 * ta_x_xyz_1[i] * fe_0 + ta_xx_yz_0[i] * fe_0 - ta_xx_yz_1[i] * fe_0 +
                          ta_xx_xyz_0[i] * pa_x[i] - ta_xx_xyz_1[i] * pc_x[i];

        ta_xxx_xzz_0[i] = 2.0 * ta_x_xzz_0[i] * fe_0 - 2.0 * ta_x_xzz_1[i] * fe_0 + ta_xx_zz_0[i] * fe_0 - ta_xx_zz_1[i] * fe_0 +
                          ta_xx_xzz_0[i] * pa_x[i] - ta_xx_xzz_1[i] * pc_x[i];

        ta_xxx_yyy_0[i] = 2.0 * ta_x_yyy_0[i] * fe_0 - 2.0 * ta_x_yyy_1[i] * fe_0 + ta_xx_yyy_0[i] * pa_x[i] - ta_xx_yyy_1[i] * pc_x[i];

        ta_xxx_yyz_0[i] = 2.0 * ta_x_yyz_0[i] * fe_0 - 2.0 * ta_x_yyz_1[i] * fe_0 + ta_xx_yyz_0[i] * pa_x[i] - ta_xx_yyz_1[i] * pc_x[i];

        ta_xxx_yzz_0[i] = 2.0 * ta_x_yzz_0[i] * fe_0 - 2.0 * ta_x_yzz_1[i] * fe_0 + ta_xx_yzz_0[i] * pa_x[i] - ta_xx_yzz_1[i] * pc_x[i];

        ta_xxx_zzz_0[i] = 2.0 * ta_x_zzz_0[i] * fe_0 - 2.0 * ta_x_zzz_1[i] * fe_0 + ta_xx_zzz_0[i] * pa_x[i] - ta_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto ta_xxy_xxx_0 = pbuffer.data(idx_npot_0_ff + 10);

    auto ta_xxy_xxy_0 = pbuffer.data(idx_npot_0_ff + 11);

    auto ta_xxy_xxz_0 = pbuffer.data(idx_npot_0_ff + 12);

    auto ta_xxy_xyy_0 = pbuffer.data(idx_npot_0_ff + 13);

    auto ta_xxy_xyz_0 = pbuffer.data(idx_npot_0_ff + 14);

    auto ta_xxy_xzz_0 = pbuffer.data(idx_npot_0_ff + 15);

    auto ta_xxy_yyy_0 = pbuffer.data(idx_npot_0_ff + 16);

    auto ta_xxy_yyz_0 = pbuffer.data(idx_npot_0_ff + 17);

    auto ta_xxy_yzz_0 = pbuffer.data(idx_npot_0_ff + 18);

    auto ta_xxy_zzz_0 = pbuffer.data(idx_npot_0_ff + 19);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xx_xx_0,   \
                             ta_xx_xx_1,   \
                             ta_xx_xxx_0,  \
                             ta_xx_xxx_1,  \
                             ta_xx_xxy_0,  \
                             ta_xx_xxy_1,  \
                             ta_xx_xxz_0,  \
                             ta_xx_xxz_1,  \
                             ta_xx_xy_0,   \
                             ta_xx_xy_1,   \
                             ta_xx_xyy_0,  \
                             ta_xx_xyy_1,  \
                             ta_xx_xyz_0,  \
                             ta_xx_xyz_1,  \
                             ta_xx_xz_0,   \
                             ta_xx_xz_1,   \
                             ta_xx_xzz_0,  \
                             ta_xx_xzz_1,  \
                             ta_xx_zzz_0,  \
                             ta_xx_zzz_1,  \
                             ta_xxy_xxx_0, \
                             ta_xxy_xxy_0, \
                             ta_xxy_xxz_0, \
                             ta_xxy_xyy_0, \
                             ta_xxy_xyz_0, \
                             ta_xxy_xzz_0, \
                             ta_xxy_yyy_0, \
                             ta_xxy_yyz_0, \
                             ta_xxy_yzz_0, \
                             ta_xxy_zzz_0, \
                             ta_xy_yyy_0,  \
                             ta_xy_yyy_1,  \
                             ta_xy_yyz_0,  \
                             ta_xy_yyz_1,  \
                             ta_xy_yzz_0,  \
                             ta_xy_yzz_1,  \
                             ta_y_yyy_0,   \
                             ta_y_yyy_1,   \
                             ta_y_yyz_0,   \
                             ta_y_yyz_1,   \
                             ta_y_yzz_0,   \
                             ta_y_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxy_xxx_0[i] = ta_xx_xxx_0[i] * pa_y[i] - ta_xx_xxx_1[i] * pc_y[i];

        ta_xxy_xxy_0[i] = ta_xx_xx_0[i] * fe_0 - ta_xx_xx_1[i] * fe_0 + ta_xx_xxy_0[i] * pa_y[i] - ta_xx_xxy_1[i] * pc_y[i];

        ta_xxy_xxz_0[i] = ta_xx_xxz_0[i] * pa_y[i] - ta_xx_xxz_1[i] * pc_y[i];

        ta_xxy_xyy_0[i] = 2.0 * ta_xx_xy_0[i] * fe_0 - 2.0 * ta_xx_xy_1[i] * fe_0 + ta_xx_xyy_0[i] * pa_y[i] - ta_xx_xyy_1[i] * pc_y[i];

        ta_xxy_xyz_0[i] = ta_xx_xz_0[i] * fe_0 - ta_xx_xz_1[i] * fe_0 + ta_xx_xyz_0[i] * pa_y[i] - ta_xx_xyz_1[i] * pc_y[i];

        ta_xxy_xzz_0[i] = ta_xx_xzz_0[i] * pa_y[i] - ta_xx_xzz_1[i] * pc_y[i];

        ta_xxy_yyy_0[i] = ta_y_yyy_0[i] * fe_0 - ta_y_yyy_1[i] * fe_0 + ta_xy_yyy_0[i] * pa_x[i] - ta_xy_yyy_1[i] * pc_x[i];

        ta_xxy_yyz_0[i] = ta_y_yyz_0[i] * fe_0 - ta_y_yyz_1[i] * fe_0 + ta_xy_yyz_0[i] * pa_x[i] - ta_xy_yyz_1[i] * pc_x[i];

        ta_xxy_yzz_0[i] = ta_y_yzz_0[i] * fe_0 - ta_y_yzz_1[i] * fe_0 + ta_xy_yzz_0[i] * pa_x[i] - ta_xy_yzz_1[i] * pc_x[i];

        ta_xxy_zzz_0[i] = ta_xx_zzz_0[i] * pa_y[i] - ta_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto ta_xxz_xxx_0 = pbuffer.data(idx_npot_0_ff + 20);

    auto ta_xxz_xxy_0 = pbuffer.data(idx_npot_0_ff + 21);

    auto ta_xxz_xxz_0 = pbuffer.data(idx_npot_0_ff + 22);

    auto ta_xxz_xyy_0 = pbuffer.data(idx_npot_0_ff + 23);

    auto ta_xxz_xyz_0 = pbuffer.data(idx_npot_0_ff + 24);

    auto ta_xxz_xzz_0 = pbuffer.data(idx_npot_0_ff + 25);

    auto ta_xxz_yyy_0 = pbuffer.data(idx_npot_0_ff + 26);

    auto ta_xxz_yyz_0 = pbuffer.data(idx_npot_0_ff + 27);

    auto ta_xxz_yzz_0 = pbuffer.data(idx_npot_0_ff + 28);

    auto ta_xxz_zzz_0 = pbuffer.data(idx_npot_0_ff + 29);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xx_xx_0,   \
                             ta_xx_xx_1,   \
                             ta_xx_xxx_0,  \
                             ta_xx_xxx_1,  \
                             ta_xx_xxy_0,  \
                             ta_xx_xxy_1,  \
                             ta_xx_xxz_0,  \
                             ta_xx_xxz_1,  \
                             ta_xx_xy_0,   \
                             ta_xx_xy_1,   \
                             ta_xx_xyy_0,  \
                             ta_xx_xyy_1,  \
                             ta_xx_xyz_0,  \
                             ta_xx_xyz_1,  \
                             ta_xx_xz_0,   \
                             ta_xx_xz_1,   \
                             ta_xx_xzz_0,  \
                             ta_xx_xzz_1,  \
                             ta_xx_yyy_0,  \
                             ta_xx_yyy_1,  \
                             ta_xxz_xxx_0, \
                             ta_xxz_xxy_0, \
                             ta_xxz_xxz_0, \
                             ta_xxz_xyy_0, \
                             ta_xxz_xyz_0, \
                             ta_xxz_xzz_0, \
                             ta_xxz_yyy_0, \
                             ta_xxz_yyz_0, \
                             ta_xxz_yzz_0, \
                             ta_xxz_zzz_0, \
                             ta_xz_yyz_0,  \
                             ta_xz_yyz_1,  \
                             ta_xz_yzz_0,  \
                             ta_xz_yzz_1,  \
                             ta_xz_zzz_0,  \
                             ta_xz_zzz_1,  \
                             ta_z_yyz_0,   \
                             ta_z_yyz_1,   \
                             ta_z_yzz_0,   \
                             ta_z_yzz_1,   \
                             ta_z_zzz_0,   \
                             ta_z_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxz_xxx_0[i] = ta_xx_xxx_0[i] * pa_z[i] - ta_xx_xxx_1[i] * pc_z[i];

        ta_xxz_xxy_0[i] = ta_xx_xxy_0[i] * pa_z[i] - ta_xx_xxy_1[i] * pc_z[i];

        ta_xxz_xxz_0[i] = ta_xx_xx_0[i] * fe_0 - ta_xx_xx_1[i] * fe_0 + ta_xx_xxz_0[i] * pa_z[i] - ta_xx_xxz_1[i] * pc_z[i];

        ta_xxz_xyy_0[i] = ta_xx_xyy_0[i] * pa_z[i] - ta_xx_xyy_1[i] * pc_z[i];

        ta_xxz_xyz_0[i] = ta_xx_xy_0[i] * fe_0 - ta_xx_xy_1[i] * fe_0 + ta_xx_xyz_0[i] * pa_z[i] - ta_xx_xyz_1[i] * pc_z[i];

        ta_xxz_xzz_0[i] = 2.0 * ta_xx_xz_0[i] * fe_0 - 2.0 * ta_xx_xz_1[i] * fe_0 + ta_xx_xzz_0[i] * pa_z[i] - ta_xx_xzz_1[i] * pc_z[i];

        ta_xxz_yyy_0[i] = ta_xx_yyy_0[i] * pa_z[i] - ta_xx_yyy_1[i] * pc_z[i];

        ta_xxz_yyz_0[i] = ta_z_yyz_0[i] * fe_0 - ta_z_yyz_1[i] * fe_0 + ta_xz_yyz_0[i] * pa_x[i] - ta_xz_yyz_1[i] * pc_x[i];

        ta_xxz_yzz_0[i] = ta_z_yzz_0[i] * fe_0 - ta_z_yzz_1[i] * fe_0 + ta_xz_yzz_0[i] * pa_x[i] - ta_xz_yzz_1[i] * pc_x[i];

        ta_xxz_zzz_0[i] = ta_z_zzz_0[i] * fe_0 - ta_z_zzz_1[i] * fe_0 + ta_xz_zzz_0[i] * pa_x[i] - ta_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto ta_xyy_xxx_0 = pbuffer.data(idx_npot_0_ff + 30);

    auto ta_xyy_xxy_0 = pbuffer.data(idx_npot_0_ff + 31);

    auto ta_xyy_xxz_0 = pbuffer.data(idx_npot_0_ff + 32);

    auto ta_xyy_xyy_0 = pbuffer.data(idx_npot_0_ff + 33);

    auto ta_xyy_xyz_0 = pbuffer.data(idx_npot_0_ff + 34);

    auto ta_xyy_xzz_0 = pbuffer.data(idx_npot_0_ff + 35);

    auto ta_xyy_yyy_0 = pbuffer.data(idx_npot_0_ff + 36);

    auto ta_xyy_yyz_0 = pbuffer.data(idx_npot_0_ff + 37);

    auto ta_xyy_yzz_0 = pbuffer.data(idx_npot_0_ff + 38);

    auto ta_xyy_zzz_0 = pbuffer.data(idx_npot_0_ff + 39);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xyy_xxx_0, \
                             ta_xyy_xxy_0, \
                             ta_xyy_xxz_0, \
                             ta_xyy_xyy_0, \
                             ta_xyy_xyz_0, \
                             ta_xyy_xzz_0, \
                             ta_xyy_yyy_0, \
                             ta_xyy_yyz_0, \
                             ta_xyy_yzz_0, \
                             ta_xyy_zzz_0, \
                             ta_yy_xx_0,   \
                             ta_yy_xx_1,   \
                             ta_yy_xxx_0,  \
                             ta_yy_xxx_1,  \
                             ta_yy_xxy_0,  \
                             ta_yy_xxy_1,  \
                             ta_yy_xxz_0,  \
                             ta_yy_xxz_1,  \
                             ta_yy_xy_0,   \
                             ta_yy_xy_1,   \
                             ta_yy_xyy_0,  \
                             ta_yy_xyy_1,  \
                             ta_yy_xyz_0,  \
                             ta_yy_xyz_1,  \
                             ta_yy_xz_0,   \
                             ta_yy_xz_1,   \
                             ta_yy_xzz_0,  \
                             ta_yy_xzz_1,  \
                             ta_yy_yy_0,   \
                             ta_yy_yy_1,   \
                             ta_yy_yyy_0,  \
                             ta_yy_yyy_1,  \
                             ta_yy_yyz_0,  \
                             ta_yy_yyz_1,  \
                             ta_yy_yz_0,   \
                             ta_yy_yz_1,   \
                             ta_yy_yzz_0,  \
                             ta_yy_yzz_1,  \
                             ta_yy_zz_0,   \
                             ta_yy_zz_1,   \
                             ta_yy_zzz_0,  \
                             ta_yy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyy_xxx_0[i] = 3.0 * ta_yy_xx_0[i] * fe_0 - 3.0 * ta_yy_xx_1[i] * fe_0 + ta_yy_xxx_0[i] * pa_x[i] - ta_yy_xxx_1[i] * pc_x[i];

        ta_xyy_xxy_0[i] = 2.0 * ta_yy_xy_0[i] * fe_0 - 2.0 * ta_yy_xy_1[i] * fe_0 + ta_yy_xxy_0[i] * pa_x[i] - ta_yy_xxy_1[i] * pc_x[i];

        ta_xyy_xxz_0[i] = 2.0 * ta_yy_xz_0[i] * fe_0 - 2.0 * ta_yy_xz_1[i] * fe_0 + ta_yy_xxz_0[i] * pa_x[i] - ta_yy_xxz_1[i] * pc_x[i];

        ta_xyy_xyy_0[i] = ta_yy_yy_0[i] * fe_0 - ta_yy_yy_1[i] * fe_0 + ta_yy_xyy_0[i] * pa_x[i] - ta_yy_xyy_1[i] * pc_x[i];

        ta_xyy_xyz_0[i] = ta_yy_yz_0[i] * fe_0 - ta_yy_yz_1[i] * fe_0 + ta_yy_xyz_0[i] * pa_x[i] - ta_yy_xyz_1[i] * pc_x[i];

        ta_xyy_xzz_0[i] = ta_yy_zz_0[i] * fe_0 - ta_yy_zz_1[i] * fe_0 + ta_yy_xzz_0[i] * pa_x[i] - ta_yy_xzz_1[i] * pc_x[i];

        ta_xyy_yyy_0[i] = ta_yy_yyy_0[i] * pa_x[i] - ta_yy_yyy_1[i] * pc_x[i];

        ta_xyy_yyz_0[i] = ta_yy_yyz_0[i] * pa_x[i] - ta_yy_yyz_1[i] * pc_x[i];

        ta_xyy_yzz_0[i] = ta_yy_yzz_0[i] * pa_x[i] - ta_yy_yzz_1[i] * pc_x[i];

        ta_xyy_zzz_0[i] = ta_yy_zzz_0[i] * pa_x[i] - ta_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto ta_xyz_xxx_0 = pbuffer.data(idx_npot_0_ff + 40);

    auto ta_xyz_xxy_0 = pbuffer.data(idx_npot_0_ff + 41);

    auto ta_xyz_xxz_0 = pbuffer.data(idx_npot_0_ff + 42);

    auto ta_xyz_xyy_0 = pbuffer.data(idx_npot_0_ff + 43);

    auto ta_xyz_xyz_0 = pbuffer.data(idx_npot_0_ff + 44);

    auto ta_xyz_xzz_0 = pbuffer.data(idx_npot_0_ff + 45);

    auto ta_xyz_yyy_0 = pbuffer.data(idx_npot_0_ff + 46);

    auto ta_xyz_yyz_0 = pbuffer.data(idx_npot_0_ff + 47);

    auto ta_xyz_yzz_0 = pbuffer.data(idx_npot_0_ff + 48);

    auto ta_xyz_zzz_0 = pbuffer.data(idx_npot_0_ff + 49);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             pc_x,         \
                             pc_y,         \
                             pc_z,         \
                             ta_xy_xxy_0,  \
                             ta_xy_xxy_1,  \
                             ta_xy_xyy_0,  \
                             ta_xy_xyy_1,  \
                             ta_xyz_xxx_0, \
                             ta_xyz_xxy_0, \
                             ta_xyz_xxz_0, \
                             ta_xyz_xyy_0, \
                             ta_xyz_xyz_0, \
                             ta_xyz_xzz_0, \
                             ta_xyz_yyy_0, \
                             ta_xyz_yyz_0, \
                             ta_xyz_yzz_0, \
                             ta_xyz_zzz_0, \
                             ta_xz_xxx_0,  \
                             ta_xz_xxx_1,  \
                             ta_xz_xxz_0,  \
                             ta_xz_xxz_1,  \
                             ta_xz_xzz_0,  \
                             ta_xz_xzz_1,  \
                             ta_yz_xyz_0,  \
                             ta_yz_xyz_1,  \
                             ta_yz_yyy_0,  \
                             ta_yz_yyy_1,  \
                             ta_yz_yyz_0,  \
                             ta_yz_yyz_1,  \
                             ta_yz_yz_0,   \
                             ta_yz_yz_1,   \
                             ta_yz_yzz_0,  \
                             ta_yz_yzz_1,  \
                             ta_yz_zzz_0,  \
                             ta_yz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyz_xxx_0[i] = ta_xz_xxx_0[i] * pa_y[i] - ta_xz_xxx_1[i] * pc_y[i];

        ta_xyz_xxy_0[i] = ta_xy_xxy_0[i] * pa_z[i] - ta_xy_xxy_1[i] * pc_z[i];

        ta_xyz_xxz_0[i] = ta_xz_xxz_0[i] * pa_y[i] - ta_xz_xxz_1[i] * pc_y[i];

        ta_xyz_xyy_0[i] = ta_xy_xyy_0[i] * pa_z[i] - ta_xy_xyy_1[i] * pc_z[i];

        ta_xyz_xyz_0[i] = ta_yz_yz_0[i] * fe_0 - ta_yz_yz_1[i] * fe_0 + ta_yz_xyz_0[i] * pa_x[i] - ta_yz_xyz_1[i] * pc_x[i];

        ta_xyz_xzz_0[i] = ta_xz_xzz_0[i] * pa_y[i] - ta_xz_xzz_1[i] * pc_y[i];

        ta_xyz_yyy_0[i] = ta_yz_yyy_0[i] * pa_x[i] - ta_yz_yyy_1[i] * pc_x[i];

        ta_xyz_yyz_0[i] = ta_yz_yyz_0[i] * pa_x[i] - ta_yz_yyz_1[i] * pc_x[i];

        ta_xyz_yzz_0[i] = ta_yz_yzz_0[i] * pa_x[i] - ta_yz_yzz_1[i] * pc_x[i];

        ta_xyz_zzz_0[i] = ta_yz_zzz_0[i] * pa_x[i] - ta_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto ta_xzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 50);

    auto ta_xzz_xxy_0 = pbuffer.data(idx_npot_0_ff + 51);

    auto ta_xzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 52);

    auto ta_xzz_xyy_0 = pbuffer.data(idx_npot_0_ff + 53);

    auto ta_xzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 54);

    auto ta_xzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 55);

    auto ta_xzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 56);

    auto ta_xzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 57);

    auto ta_xzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 58);

    auto ta_xzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 59);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xzz_xxx_0, \
                             ta_xzz_xxy_0, \
                             ta_xzz_xxz_0, \
                             ta_xzz_xyy_0, \
                             ta_xzz_xyz_0, \
                             ta_xzz_xzz_0, \
                             ta_xzz_yyy_0, \
                             ta_xzz_yyz_0, \
                             ta_xzz_yzz_0, \
                             ta_xzz_zzz_0, \
                             ta_zz_xx_0,   \
                             ta_zz_xx_1,   \
                             ta_zz_xxx_0,  \
                             ta_zz_xxx_1,  \
                             ta_zz_xxy_0,  \
                             ta_zz_xxy_1,  \
                             ta_zz_xxz_0,  \
                             ta_zz_xxz_1,  \
                             ta_zz_xy_0,   \
                             ta_zz_xy_1,   \
                             ta_zz_xyy_0,  \
                             ta_zz_xyy_1,  \
                             ta_zz_xyz_0,  \
                             ta_zz_xyz_1,  \
                             ta_zz_xz_0,   \
                             ta_zz_xz_1,   \
                             ta_zz_xzz_0,  \
                             ta_zz_xzz_1,  \
                             ta_zz_yy_0,   \
                             ta_zz_yy_1,   \
                             ta_zz_yyy_0,  \
                             ta_zz_yyy_1,  \
                             ta_zz_yyz_0,  \
                             ta_zz_yyz_1,  \
                             ta_zz_yz_0,   \
                             ta_zz_yz_1,   \
                             ta_zz_yzz_0,  \
                             ta_zz_yzz_1,  \
                             ta_zz_zz_0,   \
                             ta_zz_zz_1,   \
                             ta_zz_zzz_0,  \
                             ta_zz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzz_xxx_0[i] = 3.0 * ta_zz_xx_0[i] * fe_0 - 3.0 * ta_zz_xx_1[i] * fe_0 + ta_zz_xxx_0[i] * pa_x[i] - ta_zz_xxx_1[i] * pc_x[i];

        ta_xzz_xxy_0[i] = 2.0 * ta_zz_xy_0[i] * fe_0 - 2.0 * ta_zz_xy_1[i] * fe_0 + ta_zz_xxy_0[i] * pa_x[i] - ta_zz_xxy_1[i] * pc_x[i];

        ta_xzz_xxz_0[i] = 2.0 * ta_zz_xz_0[i] * fe_0 - 2.0 * ta_zz_xz_1[i] * fe_0 + ta_zz_xxz_0[i] * pa_x[i] - ta_zz_xxz_1[i] * pc_x[i];

        ta_xzz_xyy_0[i] = ta_zz_yy_0[i] * fe_0 - ta_zz_yy_1[i] * fe_0 + ta_zz_xyy_0[i] * pa_x[i] - ta_zz_xyy_1[i] * pc_x[i];

        ta_xzz_xyz_0[i] = ta_zz_yz_0[i] * fe_0 - ta_zz_yz_1[i] * fe_0 + ta_zz_xyz_0[i] * pa_x[i] - ta_zz_xyz_1[i] * pc_x[i];

        ta_xzz_xzz_0[i] = ta_zz_zz_0[i] * fe_0 - ta_zz_zz_1[i] * fe_0 + ta_zz_xzz_0[i] * pa_x[i] - ta_zz_xzz_1[i] * pc_x[i];

        ta_xzz_yyy_0[i] = ta_zz_yyy_0[i] * pa_x[i] - ta_zz_yyy_1[i] * pc_x[i];

        ta_xzz_yyz_0[i] = ta_zz_yyz_0[i] * pa_x[i] - ta_zz_yyz_1[i] * pc_x[i];

        ta_xzz_yzz_0[i] = ta_zz_yzz_0[i] * pa_x[i] - ta_zz_yzz_1[i] * pc_x[i];

        ta_xzz_zzz_0[i] = ta_zz_zzz_0[i] * pa_x[i] - ta_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto ta_yyy_xxx_0 = pbuffer.data(idx_npot_0_ff + 60);

    auto ta_yyy_xxy_0 = pbuffer.data(idx_npot_0_ff + 61);

    auto ta_yyy_xxz_0 = pbuffer.data(idx_npot_0_ff + 62);

    auto ta_yyy_xyy_0 = pbuffer.data(idx_npot_0_ff + 63);

    auto ta_yyy_xyz_0 = pbuffer.data(idx_npot_0_ff + 64);

    auto ta_yyy_xzz_0 = pbuffer.data(idx_npot_0_ff + 65);

    auto ta_yyy_yyy_0 = pbuffer.data(idx_npot_0_ff + 66);

    auto ta_yyy_yyz_0 = pbuffer.data(idx_npot_0_ff + 67);

    auto ta_yyy_yzz_0 = pbuffer.data(idx_npot_0_ff + 68);

    auto ta_yyy_zzz_0 = pbuffer.data(idx_npot_0_ff + 69);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta_y_xxx_0,   \
                             ta_y_xxx_1,   \
                             ta_y_xxy_0,   \
                             ta_y_xxy_1,   \
                             ta_y_xxz_0,   \
                             ta_y_xxz_1,   \
                             ta_y_xyy_0,   \
                             ta_y_xyy_1,   \
                             ta_y_xyz_0,   \
                             ta_y_xyz_1,   \
                             ta_y_xzz_0,   \
                             ta_y_xzz_1,   \
                             ta_y_yyy_0,   \
                             ta_y_yyy_1,   \
                             ta_y_yyz_0,   \
                             ta_y_yyz_1,   \
                             ta_y_yzz_0,   \
                             ta_y_yzz_1,   \
                             ta_y_zzz_0,   \
                             ta_y_zzz_1,   \
                             ta_yy_xx_0,   \
                             ta_yy_xx_1,   \
                             ta_yy_xxx_0,  \
                             ta_yy_xxx_1,  \
                             ta_yy_xxy_0,  \
                             ta_yy_xxy_1,  \
                             ta_yy_xxz_0,  \
                             ta_yy_xxz_1,  \
                             ta_yy_xy_0,   \
                             ta_yy_xy_1,   \
                             ta_yy_xyy_0,  \
                             ta_yy_xyy_1,  \
                             ta_yy_xyz_0,  \
                             ta_yy_xyz_1,  \
                             ta_yy_xz_0,   \
                             ta_yy_xz_1,   \
                             ta_yy_xzz_0,  \
                             ta_yy_xzz_1,  \
                             ta_yy_yy_0,   \
                             ta_yy_yy_1,   \
                             ta_yy_yyy_0,  \
                             ta_yy_yyy_1,  \
                             ta_yy_yyz_0,  \
                             ta_yy_yyz_1,  \
                             ta_yy_yz_0,   \
                             ta_yy_yz_1,   \
                             ta_yy_yzz_0,  \
                             ta_yy_yzz_1,  \
                             ta_yy_zz_0,   \
                             ta_yy_zz_1,   \
                             ta_yy_zzz_0,  \
                             ta_yy_zzz_1,  \
                             ta_yyy_xxx_0, \
                             ta_yyy_xxy_0, \
                             ta_yyy_xxz_0, \
                             ta_yyy_xyy_0, \
                             ta_yyy_xyz_0, \
                             ta_yyy_xzz_0, \
                             ta_yyy_yyy_0, \
                             ta_yyy_yyz_0, \
                             ta_yyy_yzz_0, \
                             ta_yyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyy_xxx_0[i] = 2.0 * ta_y_xxx_0[i] * fe_0 - 2.0 * ta_y_xxx_1[i] * fe_0 + ta_yy_xxx_0[i] * pa_y[i] - ta_yy_xxx_1[i] * pc_y[i];

        ta_yyy_xxy_0[i] = 2.0 * ta_y_xxy_0[i] * fe_0 - 2.0 * ta_y_xxy_1[i] * fe_0 + ta_yy_xx_0[i] * fe_0 - ta_yy_xx_1[i] * fe_0 +
                          ta_yy_xxy_0[i] * pa_y[i] - ta_yy_xxy_1[i] * pc_y[i];

        ta_yyy_xxz_0[i] = 2.0 * ta_y_xxz_0[i] * fe_0 - 2.0 * ta_y_xxz_1[i] * fe_0 + ta_yy_xxz_0[i] * pa_y[i] - ta_yy_xxz_1[i] * pc_y[i];

        ta_yyy_xyy_0[i] = 2.0 * ta_y_xyy_0[i] * fe_0 - 2.0 * ta_y_xyy_1[i] * fe_0 + 2.0 * ta_yy_xy_0[i] * fe_0 - 2.0 * ta_yy_xy_1[i] * fe_0 +
                          ta_yy_xyy_0[i] * pa_y[i] - ta_yy_xyy_1[i] * pc_y[i];

        ta_yyy_xyz_0[i] = 2.0 * ta_y_xyz_0[i] * fe_0 - 2.0 * ta_y_xyz_1[i] * fe_0 + ta_yy_xz_0[i] * fe_0 - ta_yy_xz_1[i] * fe_0 +
                          ta_yy_xyz_0[i] * pa_y[i] - ta_yy_xyz_1[i] * pc_y[i];

        ta_yyy_xzz_0[i] = 2.0 * ta_y_xzz_0[i] * fe_0 - 2.0 * ta_y_xzz_1[i] * fe_0 + ta_yy_xzz_0[i] * pa_y[i] - ta_yy_xzz_1[i] * pc_y[i];

        ta_yyy_yyy_0[i] = 2.0 * ta_y_yyy_0[i] * fe_0 - 2.0 * ta_y_yyy_1[i] * fe_0 + 3.0 * ta_yy_yy_0[i] * fe_0 - 3.0 * ta_yy_yy_1[i] * fe_0 +
                          ta_yy_yyy_0[i] * pa_y[i] - ta_yy_yyy_1[i] * pc_y[i];

        ta_yyy_yyz_0[i] = 2.0 * ta_y_yyz_0[i] * fe_0 - 2.0 * ta_y_yyz_1[i] * fe_0 + 2.0 * ta_yy_yz_0[i] * fe_0 - 2.0 * ta_yy_yz_1[i] * fe_0 +
                          ta_yy_yyz_0[i] * pa_y[i] - ta_yy_yyz_1[i] * pc_y[i];

        ta_yyy_yzz_0[i] = 2.0 * ta_y_yzz_0[i] * fe_0 - 2.0 * ta_y_yzz_1[i] * fe_0 + ta_yy_zz_0[i] * fe_0 - ta_yy_zz_1[i] * fe_0 +
                          ta_yy_yzz_0[i] * pa_y[i] - ta_yy_yzz_1[i] * pc_y[i];

        ta_yyy_zzz_0[i] = 2.0 * ta_y_zzz_0[i] * fe_0 - 2.0 * ta_y_zzz_1[i] * fe_0 + ta_yy_zzz_0[i] * pa_y[i] - ta_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto ta_yyz_xxx_0 = pbuffer.data(idx_npot_0_ff + 70);

    auto ta_yyz_xxy_0 = pbuffer.data(idx_npot_0_ff + 71);

    auto ta_yyz_xxz_0 = pbuffer.data(idx_npot_0_ff + 72);

    auto ta_yyz_xyy_0 = pbuffer.data(idx_npot_0_ff + 73);

    auto ta_yyz_xyz_0 = pbuffer.data(idx_npot_0_ff + 74);

    auto ta_yyz_xzz_0 = pbuffer.data(idx_npot_0_ff + 75);

    auto ta_yyz_yyy_0 = pbuffer.data(idx_npot_0_ff + 76);

    auto ta_yyz_yyz_0 = pbuffer.data(idx_npot_0_ff + 77);

    auto ta_yyz_yzz_0 = pbuffer.data(idx_npot_0_ff + 78);

    auto ta_yyz_zzz_0 = pbuffer.data(idx_npot_0_ff + 79);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_yy_xxx_0,  \
                             ta_yy_xxx_1,  \
                             ta_yy_xxy_0,  \
                             ta_yy_xxy_1,  \
                             ta_yy_xy_0,   \
                             ta_yy_xy_1,   \
                             ta_yy_xyy_0,  \
                             ta_yy_xyy_1,  \
                             ta_yy_xyz_0,  \
                             ta_yy_xyz_1,  \
                             ta_yy_yy_0,   \
                             ta_yy_yy_1,   \
                             ta_yy_yyy_0,  \
                             ta_yy_yyy_1,  \
                             ta_yy_yyz_0,  \
                             ta_yy_yyz_1,  \
                             ta_yy_yz_0,   \
                             ta_yy_yz_1,   \
                             ta_yy_yzz_0,  \
                             ta_yy_yzz_1,  \
                             ta_yyz_xxx_0, \
                             ta_yyz_xxy_0, \
                             ta_yyz_xxz_0, \
                             ta_yyz_xyy_0, \
                             ta_yyz_xyz_0, \
                             ta_yyz_xzz_0, \
                             ta_yyz_yyy_0, \
                             ta_yyz_yyz_0, \
                             ta_yyz_yzz_0, \
                             ta_yyz_zzz_0, \
                             ta_yz_xxz_0,  \
                             ta_yz_xxz_1,  \
                             ta_yz_xzz_0,  \
                             ta_yz_xzz_1,  \
                             ta_yz_zzz_0,  \
                             ta_yz_zzz_1,  \
                             ta_z_xxz_0,   \
                             ta_z_xxz_1,   \
                             ta_z_xzz_0,   \
                             ta_z_xzz_1,   \
                             ta_z_zzz_0,   \
                             ta_z_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyz_xxx_0[i] = ta_yy_xxx_0[i] * pa_z[i] - ta_yy_xxx_1[i] * pc_z[i];

        ta_yyz_xxy_0[i] = ta_yy_xxy_0[i] * pa_z[i] - ta_yy_xxy_1[i] * pc_z[i];

        ta_yyz_xxz_0[i] = ta_z_xxz_0[i] * fe_0 - ta_z_xxz_1[i] * fe_0 + ta_yz_xxz_0[i] * pa_y[i] - ta_yz_xxz_1[i] * pc_y[i];

        ta_yyz_xyy_0[i] = ta_yy_xyy_0[i] * pa_z[i] - ta_yy_xyy_1[i] * pc_z[i];

        ta_yyz_xyz_0[i] = ta_yy_xy_0[i] * fe_0 - ta_yy_xy_1[i] * fe_0 + ta_yy_xyz_0[i] * pa_z[i] - ta_yy_xyz_1[i] * pc_z[i];

        ta_yyz_xzz_0[i] = ta_z_xzz_0[i] * fe_0 - ta_z_xzz_1[i] * fe_0 + ta_yz_xzz_0[i] * pa_y[i] - ta_yz_xzz_1[i] * pc_y[i];

        ta_yyz_yyy_0[i] = ta_yy_yyy_0[i] * pa_z[i] - ta_yy_yyy_1[i] * pc_z[i];

        ta_yyz_yyz_0[i] = ta_yy_yy_0[i] * fe_0 - ta_yy_yy_1[i] * fe_0 + ta_yy_yyz_0[i] * pa_z[i] - ta_yy_yyz_1[i] * pc_z[i];

        ta_yyz_yzz_0[i] = 2.0 * ta_yy_yz_0[i] * fe_0 - 2.0 * ta_yy_yz_1[i] * fe_0 + ta_yy_yzz_0[i] * pa_z[i] - ta_yy_yzz_1[i] * pc_z[i];

        ta_yyz_zzz_0[i] = ta_z_zzz_0[i] * fe_0 - ta_z_zzz_1[i] * fe_0 + ta_yz_zzz_0[i] * pa_y[i] - ta_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto ta_yzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 80);

    auto ta_yzz_xxy_0 = pbuffer.data(idx_npot_0_ff + 81);

    auto ta_yzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 82);

    auto ta_yzz_xyy_0 = pbuffer.data(idx_npot_0_ff + 83);

    auto ta_yzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 84);

    auto ta_yzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 85);

    auto ta_yzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 86);

    auto ta_yzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 87);

    auto ta_yzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 88);

    auto ta_yzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 89);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta_yzz_xxx_0, \
                             ta_yzz_xxy_0, \
                             ta_yzz_xxz_0, \
                             ta_yzz_xyy_0, \
                             ta_yzz_xyz_0, \
                             ta_yzz_xzz_0, \
                             ta_yzz_yyy_0, \
                             ta_yzz_yyz_0, \
                             ta_yzz_yzz_0, \
                             ta_yzz_zzz_0, \
                             ta_zz_xx_0,   \
                             ta_zz_xx_1,   \
                             ta_zz_xxx_0,  \
                             ta_zz_xxx_1,  \
                             ta_zz_xxy_0,  \
                             ta_zz_xxy_1,  \
                             ta_zz_xxz_0,  \
                             ta_zz_xxz_1,  \
                             ta_zz_xy_0,   \
                             ta_zz_xy_1,   \
                             ta_zz_xyy_0,  \
                             ta_zz_xyy_1,  \
                             ta_zz_xyz_0,  \
                             ta_zz_xyz_1,  \
                             ta_zz_xz_0,   \
                             ta_zz_xz_1,   \
                             ta_zz_xzz_0,  \
                             ta_zz_xzz_1,  \
                             ta_zz_yy_0,   \
                             ta_zz_yy_1,   \
                             ta_zz_yyy_0,  \
                             ta_zz_yyy_1,  \
                             ta_zz_yyz_0,  \
                             ta_zz_yyz_1,  \
                             ta_zz_yz_0,   \
                             ta_zz_yz_1,   \
                             ta_zz_yzz_0,  \
                             ta_zz_yzz_1,  \
                             ta_zz_zz_0,   \
                             ta_zz_zz_1,   \
                             ta_zz_zzz_0,  \
                             ta_zz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzz_xxx_0[i] = ta_zz_xxx_0[i] * pa_y[i] - ta_zz_xxx_1[i] * pc_y[i];

        ta_yzz_xxy_0[i] = ta_zz_xx_0[i] * fe_0 - ta_zz_xx_1[i] * fe_0 + ta_zz_xxy_0[i] * pa_y[i] - ta_zz_xxy_1[i] * pc_y[i];

        ta_yzz_xxz_0[i] = ta_zz_xxz_0[i] * pa_y[i] - ta_zz_xxz_1[i] * pc_y[i];

        ta_yzz_xyy_0[i] = 2.0 * ta_zz_xy_0[i] * fe_0 - 2.0 * ta_zz_xy_1[i] * fe_0 + ta_zz_xyy_0[i] * pa_y[i] - ta_zz_xyy_1[i] * pc_y[i];

        ta_yzz_xyz_0[i] = ta_zz_xz_0[i] * fe_0 - ta_zz_xz_1[i] * fe_0 + ta_zz_xyz_0[i] * pa_y[i] - ta_zz_xyz_1[i] * pc_y[i];

        ta_yzz_xzz_0[i] = ta_zz_xzz_0[i] * pa_y[i] - ta_zz_xzz_1[i] * pc_y[i];

        ta_yzz_yyy_0[i] = 3.0 * ta_zz_yy_0[i] * fe_0 - 3.0 * ta_zz_yy_1[i] * fe_0 + ta_zz_yyy_0[i] * pa_y[i] - ta_zz_yyy_1[i] * pc_y[i];

        ta_yzz_yyz_0[i] = 2.0 * ta_zz_yz_0[i] * fe_0 - 2.0 * ta_zz_yz_1[i] * fe_0 + ta_zz_yyz_0[i] * pa_y[i] - ta_zz_yyz_1[i] * pc_y[i];

        ta_yzz_yzz_0[i] = ta_zz_zz_0[i] * fe_0 - ta_zz_zz_1[i] * fe_0 + ta_zz_yzz_0[i] * pa_y[i] - ta_zz_yzz_1[i] * pc_y[i];

        ta_yzz_zzz_0[i] = ta_zz_zzz_0[i] * pa_y[i] - ta_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto ta_zzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 90);

    auto ta_zzz_xxy_0 = pbuffer.data(idx_npot_0_ff + 91);

    auto ta_zzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 92);

    auto ta_zzz_xyy_0 = pbuffer.data(idx_npot_0_ff + 93);

    auto ta_zzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 94);

    auto ta_zzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 95);

    auto ta_zzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 96);

    auto ta_zzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 97);

    auto ta_zzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 98);

    auto ta_zzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 99);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta_z_xxx_0,   \
                             ta_z_xxx_1,   \
                             ta_z_xxy_0,   \
                             ta_z_xxy_1,   \
                             ta_z_xxz_0,   \
                             ta_z_xxz_1,   \
                             ta_z_xyy_0,   \
                             ta_z_xyy_1,   \
                             ta_z_xyz_0,   \
                             ta_z_xyz_1,   \
                             ta_z_xzz_0,   \
                             ta_z_xzz_1,   \
                             ta_z_yyy_0,   \
                             ta_z_yyy_1,   \
                             ta_z_yyz_0,   \
                             ta_z_yyz_1,   \
                             ta_z_yzz_0,   \
                             ta_z_yzz_1,   \
                             ta_z_zzz_0,   \
                             ta_z_zzz_1,   \
                             ta_zz_xx_0,   \
                             ta_zz_xx_1,   \
                             ta_zz_xxx_0,  \
                             ta_zz_xxx_1,  \
                             ta_zz_xxy_0,  \
                             ta_zz_xxy_1,  \
                             ta_zz_xxz_0,  \
                             ta_zz_xxz_1,  \
                             ta_zz_xy_0,   \
                             ta_zz_xy_1,   \
                             ta_zz_xyy_0,  \
                             ta_zz_xyy_1,  \
                             ta_zz_xyz_0,  \
                             ta_zz_xyz_1,  \
                             ta_zz_xz_0,   \
                             ta_zz_xz_1,   \
                             ta_zz_xzz_0,  \
                             ta_zz_xzz_1,  \
                             ta_zz_yy_0,   \
                             ta_zz_yy_1,   \
                             ta_zz_yyy_0,  \
                             ta_zz_yyy_1,  \
                             ta_zz_yyz_0,  \
                             ta_zz_yyz_1,  \
                             ta_zz_yz_0,   \
                             ta_zz_yz_1,   \
                             ta_zz_yzz_0,  \
                             ta_zz_yzz_1,  \
                             ta_zz_zz_0,   \
                             ta_zz_zz_1,   \
                             ta_zz_zzz_0,  \
                             ta_zz_zzz_1,  \
                             ta_zzz_xxx_0, \
                             ta_zzz_xxy_0, \
                             ta_zzz_xxz_0, \
                             ta_zzz_xyy_0, \
                             ta_zzz_xyz_0, \
                             ta_zzz_xzz_0, \
                             ta_zzz_yyy_0, \
                             ta_zzz_yyz_0, \
                             ta_zzz_yzz_0, \
                             ta_zzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzz_xxx_0[i] = 2.0 * ta_z_xxx_0[i] * fe_0 - 2.0 * ta_z_xxx_1[i] * fe_0 + ta_zz_xxx_0[i] * pa_z[i] - ta_zz_xxx_1[i] * pc_z[i];

        ta_zzz_xxy_0[i] = 2.0 * ta_z_xxy_0[i] * fe_0 - 2.0 * ta_z_xxy_1[i] * fe_0 + ta_zz_xxy_0[i] * pa_z[i] - ta_zz_xxy_1[i] * pc_z[i];

        ta_zzz_xxz_0[i] = 2.0 * ta_z_xxz_0[i] * fe_0 - 2.0 * ta_z_xxz_1[i] * fe_0 + ta_zz_xx_0[i] * fe_0 - ta_zz_xx_1[i] * fe_0 +
                          ta_zz_xxz_0[i] * pa_z[i] - ta_zz_xxz_1[i] * pc_z[i];

        ta_zzz_xyy_0[i] = 2.0 * ta_z_xyy_0[i] * fe_0 - 2.0 * ta_z_xyy_1[i] * fe_0 + ta_zz_xyy_0[i] * pa_z[i] - ta_zz_xyy_1[i] * pc_z[i];

        ta_zzz_xyz_0[i] = 2.0 * ta_z_xyz_0[i] * fe_0 - 2.0 * ta_z_xyz_1[i] * fe_0 + ta_zz_xy_0[i] * fe_0 - ta_zz_xy_1[i] * fe_0 +
                          ta_zz_xyz_0[i] * pa_z[i] - ta_zz_xyz_1[i] * pc_z[i];

        ta_zzz_xzz_0[i] = 2.0 * ta_z_xzz_0[i] * fe_0 - 2.0 * ta_z_xzz_1[i] * fe_0 + 2.0 * ta_zz_xz_0[i] * fe_0 - 2.0 * ta_zz_xz_1[i] * fe_0 +
                          ta_zz_xzz_0[i] * pa_z[i] - ta_zz_xzz_1[i] * pc_z[i];

        ta_zzz_yyy_0[i] = 2.0 * ta_z_yyy_0[i] * fe_0 - 2.0 * ta_z_yyy_1[i] * fe_0 + ta_zz_yyy_0[i] * pa_z[i] - ta_zz_yyy_1[i] * pc_z[i];

        ta_zzz_yyz_0[i] = 2.0 * ta_z_yyz_0[i] * fe_0 - 2.0 * ta_z_yyz_1[i] * fe_0 + ta_zz_yy_0[i] * fe_0 - ta_zz_yy_1[i] * fe_0 +
                          ta_zz_yyz_0[i] * pa_z[i] - ta_zz_yyz_1[i] * pc_z[i];

        ta_zzz_yzz_0[i] = 2.0 * ta_z_yzz_0[i] * fe_0 - 2.0 * ta_z_yzz_1[i] * fe_0 + 2.0 * ta_zz_yz_0[i] * fe_0 - 2.0 * ta_zz_yz_1[i] * fe_0 +
                          ta_zz_yzz_0[i] * pa_z[i] - ta_zz_yzz_1[i] * pc_z[i];

        ta_zzz_zzz_0[i] = 2.0 * ta_z_zzz_0[i] * fe_0 - 2.0 * ta_z_zzz_1[i] * fe_0 + 3.0 * ta_zz_zz_0[i] * fe_0 - 3.0 * ta_zz_zz_1[i] * fe_0 +
                          ta_zz_zzz_0[i] * pa_z[i] - ta_zz_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
