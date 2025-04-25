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

#include "NuclearPotentialPrimRecFS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_fs(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_fs,
                               const size_t              idx_npot_0_ps,
                               const size_t              idx_npot_1_ps,
                               const size_t              idx_npot_0_ds,
                               const size_t              idx_npot_1_ds,
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

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_0 = pbuffer.data(idx_npot_0_ps);

    auto ta_y_0_0 = pbuffer.data(idx_npot_0_ps + 1);

    auto ta_z_0_0 = pbuffer.data(idx_npot_0_ps + 2);

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_1 = pbuffer.data(idx_npot_1_ps);

    auto ta_y_0_1 = pbuffer.data(idx_npot_1_ps + 1);

    auto ta_z_0_1 = pbuffer.data(idx_npot_1_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_0 = pbuffer.data(idx_npot_0_ds);

    auto ta_yy_0_0 = pbuffer.data(idx_npot_0_ds + 3);

    auto ta_yz_0_0 = pbuffer.data(idx_npot_0_ds + 4);

    auto ta_zz_0_0 = pbuffer.data(idx_npot_0_ds + 5);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_1 = pbuffer.data(idx_npot_1_ds);

    auto ta_yy_0_1 = pbuffer.data(idx_npot_1_ds + 3);

    auto ta_yz_0_1 = pbuffer.data(idx_npot_1_ds + 4);

    auto ta_zz_0_1 = pbuffer.data(idx_npot_1_ds + 5);

    // Set up components of targeted buffer : FS

    auto ta_xxx_0_0 = pbuffer.data(idx_npot_0_fs);

    auto ta_xxy_0_0 = pbuffer.data(idx_npot_0_fs + 1);

    auto ta_xxz_0_0 = pbuffer.data(idx_npot_0_fs + 2);

    auto ta_xyy_0_0 = pbuffer.data(idx_npot_0_fs + 3);

    auto ta_xyz_0_0 = pbuffer.data(idx_npot_0_fs + 4);

    auto ta_xzz_0_0 = pbuffer.data(idx_npot_0_fs + 5);

    auto ta_yyy_0_0 = pbuffer.data(idx_npot_0_fs + 6);

    auto ta_yyz_0_0 = pbuffer.data(idx_npot_0_fs + 7);

    auto ta_yzz_0_0 = pbuffer.data(idx_npot_0_fs + 8);

    auto ta_zzz_0_0 = pbuffer.data(idx_npot_0_fs + 9);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             pa_z,       \
                             pc_x,       \
                             pc_y,       \
                             pc_z,       \
                             ta_x_0_0,   \
                             ta_x_0_1,   \
                             ta_xx_0_0,  \
                             ta_xx_0_1,  \
                             ta_xxx_0_0, \
                             ta_xxy_0_0, \
                             ta_xxz_0_0, \
                             ta_xyy_0_0, \
                             ta_xyz_0_0, \
                             ta_xzz_0_0, \
                             ta_y_0_0,   \
                             ta_y_0_1,   \
                             ta_yy_0_0,  \
                             ta_yy_0_1,  \
                             ta_yyy_0_0, \
                             ta_yyz_0_0, \
                             ta_yz_0_0,  \
                             ta_yz_0_1,  \
                             ta_yzz_0_0, \
                             ta_z_0_0,   \
                             ta_z_0_1,   \
                             ta_zz_0_0,  \
                             ta_zz_0_1,  \
                             ta_zzz_0_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxx_0_0[i] = 2.0 * ta_x_0_0[i] * fe_0 - 2.0 * ta_x_0_1[i] * fe_0 + ta_xx_0_0[i] * pa_x[i] - ta_xx_0_1[i] * pc_x[i];

        ta_xxy_0_0[i] = ta_xx_0_0[i] * pa_y[i] - ta_xx_0_1[i] * pc_y[i];

        ta_xxz_0_0[i] = ta_xx_0_0[i] * pa_z[i] - ta_xx_0_1[i] * pc_z[i];

        ta_xyy_0_0[i] = ta_yy_0_0[i] * pa_x[i] - ta_yy_0_1[i] * pc_x[i];

        ta_xyz_0_0[i] = ta_yz_0_0[i] * pa_x[i] - ta_yz_0_1[i] * pc_x[i];

        ta_xzz_0_0[i] = ta_zz_0_0[i] * pa_x[i] - ta_zz_0_1[i] * pc_x[i];

        ta_yyy_0_0[i] = 2.0 * ta_y_0_0[i] * fe_0 - 2.0 * ta_y_0_1[i] * fe_0 + ta_yy_0_0[i] * pa_y[i] - ta_yy_0_1[i] * pc_y[i];

        ta_yyz_0_0[i] = ta_yy_0_0[i] * pa_z[i] - ta_yy_0_1[i] * pc_z[i];

        ta_yzz_0_0[i] = ta_zz_0_0[i] * pa_y[i] - ta_zz_0_1[i] * pc_y[i];

        ta_zzz_0_0[i] = 2.0 * ta_z_0_0[i] * fe_0 - 2.0 * ta_z_0_1[i] * fe_0 + ta_zz_0_0[i] * pa_z[i] - ta_zz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
