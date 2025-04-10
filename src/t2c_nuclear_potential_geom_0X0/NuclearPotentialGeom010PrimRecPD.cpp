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

#include "NuclearPotentialGeom010PrimRecPD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_pd(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_pd,
                                        const size_t              idx_npot_geom_010_0_sp,
                                        const size_t              idx_npot_geom_010_1_sp,
                                        const size_t              idx_npot_1_sd,
                                        const size_t              idx_npot_geom_010_0_sd,
                                        const size_t              idx_npot_geom_010_1_sd,
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

    // Set up components of auxiliary buffer : SP

    auto ta1_x_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp);

    auto ta1_x_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 1);

    auto ta1_x_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 2);

    auto ta1_y_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp + 3);

    auto ta1_y_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 4);

    auto ta1_y_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 5);

    auto ta1_z_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp + 6);

    auto ta1_z_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 7);

    auto ta1_z_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 8);

    // Set up components of auxiliary buffer : SP

    auto ta1_x_0_x_1 = pbuffer.data(idx_npot_geom_010_1_sp);

    auto ta1_x_0_y_1 = pbuffer.data(idx_npot_geom_010_1_sp + 1);

    auto ta1_x_0_z_1 = pbuffer.data(idx_npot_geom_010_1_sp + 2);

    auto ta1_y_0_x_1 = pbuffer.data(idx_npot_geom_010_1_sp + 3);

    auto ta1_y_0_y_1 = pbuffer.data(idx_npot_geom_010_1_sp + 4);

    auto ta1_y_0_z_1 = pbuffer.data(idx_npot_geom_010_1_sp + 5);

    auto ta1_z_0_x_1 = pbuffer.data(idx_npot_geom_010_1_sp + 6);

    auto ta1_z_0_y_1 = pbuffer.data(idx_npot_geom_010_1_sp + 7);

    auto ta1_z_0_z_1 = pbuffer.data(idx_npot_geom_010_1_sp + 8);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = pbuffer.data(idx_npot_1_sd);

    auto ta_0_xy_1 = pbuffer.data(idx_npot_1_sd + 1);

    auto ta_0_xz_1 = pbuffer.data(idx_npot_1_sd + 2);

    auto ta_0_yy_1 = pbuffer.data(idx_npot_1_sd + 3);

    auto ta_0_yz_1 = pbuffer.data(idx_npot_1_sd + 4);

    auto ta_0_zz_1 = pbuffer.data(idx_npot_1_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd);

    auto ta1_x_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 1);

    auto ta1_x_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 2);

    auto ta1_x_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 3);

    auto ta1_x_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 4);

    auto ta1_x_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 5);

    auto ta1_y_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 6);

    auto ta1_y_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 7);

    auto ta1_y_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 8);

    auto ta1_y_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 9);

    auto ta1_y_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 10);

    auto ta1_y_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 11);

    auto ta1_z_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 12);

    auto ta1_z_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 13);

    auto ta1_z_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 14);

    auto ta1_z_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 15);

    auto ta1_z_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 16);

    auto ta1_z_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 17);

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd);

    auto ta1_x_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 1);

    auto ta1_x_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 2);

    auto ta1_x_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 3);

    auto ta1_x_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 4);

    auto ta1_x_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 5);

    auto ta1_y_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 6);

    auto ta1_y_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 7);

    auto ta1_y_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 8);

    auto ta1_y_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 9);

    auto ta1_y_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 10);

    auto ta1_y_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 11);

    auto ta1_z_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 12);

    auto ta1_z_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 13);

    auto ta1_z_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 14);

    auto ta1_z_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 15);

    auto ta1_z_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 16);

    auto ta1_z_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 17);

    // Set up 0-6 components of targeted buffer : PD

    auto ta1_x_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd);

    auto ta1_x_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 1);

    auto ta1_x_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 2);

    auto ta1_x_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 3);

    auto ta1_x_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 4);

    auto ta1_x_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 5);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_x_0_x_0,  \
                             ta1_x_0_x_1,  \
                             ta1_x_0_xx_0, \
                             ta1_x_0_xx_1, \
                             ta1_x_0_xy_0, \
                             ta1_x_0_xy_1, \
                             ta1_x_0_xz_0, \
                             ta1_x_0_xz_1, \
                             ta1_x_0_y_0,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_yy_0, \
                             ta1_x_0_yy_1, \
                             ta1_x_0_yz_0, \
                             ta1_x_0_yz_1, \
                             ta1_x_0_z_0,  \
                             ta1_x_0_z_1,  \
                             ta1_x_0_zz_0, \
                             ta1_x_0_zz_1, \
                             ta1_x_x_xx_0, \
                             ta1_x_x_xy_0, \
                             ta1_x_x_xz_0, \
                             ta1_x_x_yy_0, \
                             ta1_x_x_yz_0, \
                             ta1_x_x_zz_0, \
                             ta_0_xx_1,    \
                             ta_0_xy_1,    \
                             ta_0_xz_1,    \
                             ta_0_yy_1,    \
                             ta_0_yz_1,    \
                             ta_0_zz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_x_xx_0[i] =
            2.0 * ta1_x_0_x_0[i] * fe_0 - 2.0 * ta1_x_0_x_1[i] * fe_0 + ta_0_xx_1[i] + ta1_x_0_xx_0[i] * pa_x[i] - ta1_x_0_xx_1[i] * pc_x[i];

        ta1_x_x_xy_0[i] = ta1_x_0_y_0[i] * fe_0 - ta1_x_0_y_1[i] * fe_0 + ta_0_xy_1[i] + ta1_x_0_xy_0[i] * pa_x[i] - ta1_x_0_xy_1[i] * pc_x[i];

        ta1_x_x_xz_0[i] = ta1_x_0_z_0[i] * fe_0 - ta1_x_0_z_1[i] * fe_0 + ta_0_xz_1[i] + ta1_x_0_xz_0[i] * pa_x[i] - ta1_x_0_xz_1[i] * pc_x[i];

        ta1_x_x_yy_0[i] = ta_0_yy_1[i] + ta1_x_0_yy_0[i] * pa_x[i] - ta1_x_0_yy_1[i] * pc_x[i];

        ta1_x_x_yz_0[i] = ta_0_yz_1[i] + ta1_x_0_yz_0[i] * pa_x[i] - ta1_x_0_yz_1[i] * pc_x[i];

        ta1_x_x_zz_0[i] = ta_0_zz_1[i] + ta1_x_0_zz_0[i] * pa_x[i] - ta1_x_0_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto ta1_x_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 6);

    auto ta1_x_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 7);

    auto ta1_x_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 8);

    auto ta1_x_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 9);

    auto ta1_x_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 10);

    auto ta1_x_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 11);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_x_0_x_0,  \
                             ta1_x_0_x_1,  \
                             ta1_x_0_xx_0, \
                             ta1_x_0_xx_1, \
                             ta1_x_0_xy_0, \
                             ta1_x_0_xy_1, \
                             ta1_x_0_xz_0, \
                             ta1_x_0_xz_1, \
                             ta1_x_0_y_0,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_yy_0, \
                             ta1_x_0_yy_1, \
                             ta1_x_0_yz_0, \
                             ta1_x_0_yz_1, \
                             ta1_x_0_z_0,  \
                             ta1_x_0_z_1,  \
                             ta1_x_0_zz_0, \
                             ta1_x_0_zz_1, \
                             ta1_x_y_xx_0, \
                             ta1_x_y_xy_0, \
                             ta1_x_y_xz_0, \
                             ta1_x_y_yy_0, \
                             ta1_x_y_yz_0, \
                             ta1_x_y_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_y_xx_0[i] = ta1_x_0_xx_0[i] * pa_y[i] - ta1_x_0_xx_1[i] * pc_y[i];

        ta1_x_y_xy_0[i] = ta1_x_0_x_0[i] * fe_0 - ta1_x_0_x_1[i] * fe_0 + ta1_x_0_xy_0[i] * pa_y[i] - ta1_x_0_xy_1[i] * pc_y[i];

        ta1_x_y_xz_0[i] = ta1_x_0_xz_0[i] * pa_y[i] - ta1_x_0_xz_1[i] * pc_y[i];

        ta1_x_y_yy_0[i] = 2.0 * ta1_x_0_y_0[i] * fe_0 - 2.0 * ta1_x_0_y_1[i] * fe_0 + ta1_x_0_yy_0[i] * pa_y[i] - ta1_x_0_yy_1[i] * pc_y[i];

        ta1_x_y_yz_0[i] = ta1_x_0_z_0[i] * fe_0 - ta1_x_0_z_1[i] * fe_0 + ta1_x_0_yz_0[i] * pa_y[i] - ta1_x_0_yz_1[i] * pc_y[i];

        ta1_x_y_zz_0[i] = ta1_x_0_zz_0[i] * pa_y[i] - ta1_x_0_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto ta1_x_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 12);

    auto ta1_x_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 13);

    auto ta1_x_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 14);

    auto ta1_x_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 15);

    auto ta1_x_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 16);

    auto ta1_x_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 17);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_x_0_x_0,  \
                             ta1_x_0_x_1,  \
                             ta1_x_0_xx_0, \
                             ta1_x_0_xx_1, \
                             ta1_x_0_xy_0, \
                             ta1_x_0_xy_1, \
                             ta1_x_0_xz_0, \
                             ta1_x_0_xz_1, \
                             ta1_x_0_y_0,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_yy_0, \
                             ta1_x_0_yy_1, \
                             ta1_x_0_yz_0, \
                             ta1_x_0_yz_1, \
                             ta1_x_0_z_0,  \
                             ta1_x_0_z_1,  \
                             ta1_x_0_zz_0, \
                             ta1_x_0_zz_1, \
                             ta1_x_z_xx_0, \
                             ta1_x_z_xy_0, \
                             ta1_x_z_xz_0, \
                             ta1_x_z_yy_0, \
                             ta1_x_z_yz_0, \
                             ta1_x_z_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_z_xx_0[i] = ta1_x_0_xx_0[i] * pa_z[i] - ta1_x_0_xx_1[i] * pc_z[i];

        ta1_x_z_xy_0[i] = ta1_x_0_xy_0[i] * pa_z[i] - ta1_x_0_xy_1[i] * pc_z[i];

        ta1_x_z_xz_0[i] = ta1_x_0_x_0[i] * fe_0 - ta1_x_0_x_1[i] * fe_0 + ta1_x_0_xz_0[i] * pa_z[i] - ta1_x_0_xz_1[i] * pc_z[i];

        ta1_x_z_yy_0[i] = ta1_x_0_yy_0[i] * pa_z[i] - ta1_x_0_yy_1[i] * pc_z[i];

        ta1_x_z_yz_0[i] = ta1_x_0_y_0[i] * fe_0 - ta1_x_0_y_1[i] * fe_0 + ta1_x_0_yz_0[i] * pa_z[i] - ta1_x_0_yz_1[i] * pc_z[i];

        ta1_x_z_zz_0[i] = 2.0 * ta1_x_0_z_0[i] * fe_0 - 2.0 * ta1_x_0_z_1[i] * fe_0 + ta1_x_0_zz_0[i] * pa_z[i] - ta1_x_0_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : PD

    auto ta1_y_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 18);

    auto ta1_y_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 19);

    auto ta1_y_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 20);

    auto ta1_y_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 21);

    auto ta1_y_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 22);

    auto ta1_y_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 23);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_y_0_x_0,  \
                             ta1_y_0_x_1,  \
                             ta1_y_0_xx_0, \
                             ta1_y_0_xx_1, \
                             ta1_y_0_xy_0, \
                             ta1_y_0_xy_1, \
                             ta1_y_0_xz_0, \
                             ta1_y_0_xz_1, \
                             ta1_y_0_y_0,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_yy_0, \
                             ta1_y_0_yy_1, \
                             ta1_y_0_yz_0, \
                             ta1_y_0_yz_1, \
                             ta1_y_0_z_0,  \
                             ta1_y_0_z_1,  \
                             ta1_y_0_zz_0, \
                             ta1_y_0_zz_1, \
                             ta1_y_x_xx_0, \
                             ta1_y_x_xy_0, \
                             ta1_y_x_xz_0, \
                             ta1_y_x_yy_0, \
                             ta1_y_x_yz_0, \
                             ta1_y_x_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_x_xx_0[i] = 2.0 * ta1_y_0_x_0[i] * fe_0 - 2.0 * ta1_y_0_x_1[i] * fe_0 + ta1_y_0_xx_0[i] * pa_x[i] - ta1_y_0_xx_1[i] * pc_x[i];

        ta1_y_x_xy_0[i] = ta1_y_0_y_0[i] * fe_0 - ta1_y_0_y_1[i] * fe_0 + ta1_y_0_xy_0[i] * pa_x[i] - ta1_y_0_xy_1[i] * pc_x[i];

        ta1_y_x_xz_0[i] = ta1_y_0_z_0[i] * fe_0 - ta1_y_0_z_1[i] * fe_0 + ta1_y_0_xz_0[i] * pa_x[i] - ta1_y_0_xz_1[i] * pc_x[i];

        ta1_y_x_yy_0[i] = ta1_y_0_yy_0[i] * pa_x[i] - ta1_y_0_yy_1[i] * pc_x[i];

        ta1_y_x_yz_0[i] = ta1_y_0_yz_0[i] * pa_x[i] - ta1_y_0_yz_1[i] * pc_x[i];

        ta1_y_x_zz_0[i] = ta1_y_0_zz_0[i] * pa_x[i] - ta1_y_0_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : PD

    auto ta1_y_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 24);

    auto ta1_y_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 25);

    auto ta1_y_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 26);

    auto ta1_y_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 27);

    auto ta1_y_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 28);

    auto ta1_y_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 29);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_y_0_x_0,  \
                             ta1_y_0_x_1,  \
                             ta1_y_0_xx_0, \
                             ta1_y_0_xx_1, \
                             ta1_y_0_xy_0, \
                             ta1_y_0_xy_1, \
                             ta1_y_0_xz_0, \
                             ta1_y_0_xz_1, \
                             ta1_y_0_y_0,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_yy_0, \
                             ta1_y_0_yy_1, \
                             ta1_y_0_yz_0, \
                             ta1_y_0_yz_1, \
                             ta1_y_0_z_0,  \
                             ta1_y_0_z_1,  \
                             ta1_y_0_zz_0, \
                             ta1_y_0_zz_1, \
                             ta1_y_y_xx_0, \
                             ta1_y_y_xy_0, \
                             ta1_y_y_xz_0, \
                             ta1_y_y_yy_0, \
                             ta1_y_y_yz_0, \
                             ta1_y_y_zz_0, \
                             ta_0_xx_1,    \
                             ta_0_xy_1,    \
                             ta_0_xz_1,    \
                             ta_0_yy_1,    \
                             ta_0_yz_1,    \
                             ta_0_zz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_y_xx_0[i] = ta_0_xx_1[i] + ta1_y_0_xx_0[i] * pa_y[i] - ta1_y_0_xx_1[i] * pc_y[i];

        ta1_y_y_xy_0[i] = ta1_y_0_x_0[i] * fe_0 - ta1_y_0_x_1[i] * fe_0 + ta_0_xy_1[i] + ta1_y_0_xy_0[i] * pa_y[i] - ta1_y_0_xy_1[i] * pc_y[i];

        ta1_y_y_xz_0[i] = ta_0_xz_1[i] + ta1_y_0_xz_0[i] * pa_y[i] - ta1_y_0_xz_1[i] * pc_y[i];

        ta1_y_y_yy_0[i] =
            2.0 * ta1_y_0_y_0[i] * fe_0 - 2.0 * ta1_y_0_y_1[i] * fe_0 + ta_0_yy_1[i] + ta1_y_0_yy_0[i] * pa_y[i] - ta1_y_0_yy_1[i] * pc_y[i];

        ta1_y_y_yz_0[i] = ta1_y_0_z_0[i] * fe_0 - ta1_y_0_z_1[i] * fe_0 + ta_0_yz_1[i] + ta1_y_0_yz_0[i] * pa_y[i] - ta1_y_0_yz_1[i] * pc_y[i];

        ta1_y_y_zz_0[i] = ta_0_zz_1[i] + ta1_y_0_zz_0[i] * pa_y[i] - ta1_y_0_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : PD

    auto ta1_y_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 30);

    auto ta1_y_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 31);

    auto ta1_y_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 32);

    auto ta1_y_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 33);

    auto ta1_y_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 34);

    auto ta1_y_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 35);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_y_0_x_0,  \
                             ta1_y_0_x_1,  \
                             ta1_y_0_xx_0, \
                             ta1_y_0_xx_1, \
                             ta1_y_0_xy_0, \
                             ta1_y_0_xy_1, \
                             ta1_y_0_xz_0, \
                             ta1_y_0_xz_1, \
                             ta1_y_0_y_0,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_yy_0, \
                             ta1_y_0_yy_1, \
                             ta1_y_0_yz_0, \
                             ta1_y_0_yz_1, \
                             ta1_y_0_z_0,  \
                             ta1_y_0_z_1,  \
                             ta1_y_0_zz_0, \
                             ta1_y_0_zz_1, \
                             ta1_y_z_xx_0, \
                             ta1_y_z_xy_0, \
                             ta1_y_z_xz_0, \
                             ta1_y_z_yy_0, \
                             ta1_y_z_yz_0, \
                             ta1_y_z_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_z_xx_0[i] = ta1_y_0_xx_0[i] * pa_z[i] - ta1_y_0_xx_1[i] * pc_z[i];

        ta1_y_z_xy_0[i] = ta1_y_0_xy_0[i] * pa_z[i] - ta1_y_0_xy_1[i] * pc_z[i];

        ta1_y_z_xz_0[i] = ta1_y_0_x_0[i] * fe_0 - ta1_y_0_x_1[i] * fe_0 + ta1_y_0_xz_0[i] * pa_z[i] - ta1_y_0_xz_1[i] * pc_z[i];

        ta1_y_z_yy_0[i] = ta1_y_0_yy_0[i] * pa_z[i] - ta1_y_0_yy_1[i] * pc_z[i];

        ta1_y_z_yz_0[i] = ta1_y_0_y_0[i] * fe_0 - ta1_y_0_y_1[i] * fe_0 + ta1_y_0_yz_0[i] * pa_z[i] - ta1_y_0_yz_1[i] * pc_z[i];

        ta1_y_z_zz_0[i] = 2.0 * ta1_y_0_z_0[i] * fe_0 - 2.0 * ta1_y_0_z_1[i] * fe_0 + ta1_y_0_zz_0[i] * pa_z[i] - ta1_y_0_zz_1[i] * pc_z[i];
    }

    // Set up 36-42 components of targeted buffer : PD

    auto ta1_z_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 36);

    auto ta1_z_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 37);

    auto ta1_z_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 38);

    auto ta1_z_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 39);

    auto ta1_z_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 40);

    auto ta1_z_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 41);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_z_0_x_0,  \
                             ta1_z_0_x_1,  \
                             ta1_z_0_xx_0, \
                             ta1_z_0_xx_1, \
                             ta1_z_0_xy_0, \
                             ta1_z_0_xy_1, \
                             ta1_z_0_xz_0, \
                             ta1_z_0_xz_1, \
                             ta1_z_0_y_0,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_yy_0, \
                             ta1_z_0_yy_1, \
                             ta1_z_0_yz_0, \
                             ta1_z_0_yz_1, \
                             ta1_z_0_z_0,  \
                             ta1_z_0_z_1,  \
                             ta1_z_0_zz_0, \
                             ta1_z_0_zz_1, \
                             ta1_z_x_xx_0, \
                             ta1_z_x_xy_0, \
                             ta1_z_x_xz_0, \
                             ta1_z_x_yy_0, \
                             ta1_z_x_yz_0, \
                             ta1_z_x_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_x_xx_0[i] = 2.0 * ta1_z_0_x_0[i] * fe_0 - 2.0 * ta1_z_0_x_1[i] * fe_0 + ta1_z_0_xx_0[i] * pa_x[i] - ta1_z_0_xx_1[i] * pc_x[i];

        ta1_z_x_xy_0[i] = ta1_z_0_y_0[i] * fe_0 - ta1_z_0_y_1[i] * fe_0 + ta1_z_0_xy_0[i] * pa_x[i] - ta1_z_0_xy_1[i] * pc_x[i];

        ta1_z_x_xz_0[i] = ta1_z_0_z_0[i] * fe_0 - ta1_z_0_z_1[i] * fe_0 + ta1_z_0_xz_0[i] * pa_x[i] - ta1_z_0_xz_1[i] * pc_x[i];

        ta1_z_x_yy_0[i] = ta1_z_0_yy_0[i] * pa_x[i] - ta1_z_0_yy_1[i] * pc_x[i];

        ta1_z_x_yz_0[i] = ta1_z_0_yz_0[i] * pa_x[i] - ta1_z_0_yz_1[i] * pc_x[i];

        ta1_z_x_zz_0[i] = ta1_z_0_zz_0[i] * pa_x[i] - ta1_z_0_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : PD

    auto ta1_z_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 42);

    auto ta1_z_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 43);

    auto ta1_z_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 44);

    auto ta1_z_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 45);

    auto ta1_z_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 46);

    auto ta1_z_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 47);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_z_0_x_0,  \
                             ta1_z_0_x_1,  \
                             ta1_z_0_xx_0, \
                             ta1_z_0_xx_1, \
                             ta1_z_0_xy_0, \
                             ta1_z_0_xy_1, \
                             ta1_z_0_xz_0, \
                             ta1_z_0_xz_1, \
                             ta1_z_0_y_0,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_yy_0, \
                             ta1_z_0_yy_1, \
                             ta1_z_0_yz_0, \
                             ta1_z_0_yz_1, \
                             ta1_z_0_z_0,  \
                             ta1_z_0_z_1,  \
                             ta1_z_0_zz_0, \
                             ta1_z_0_zz_1, \
                             ta1_z_y_xx_0, \
                             ta1_z_y_xy_0, \
                             ta1_z_y_xz_0, \
                             ta1_z_y_yy_0, \
                             ta1_z_y_yz_0, \
                             ta1_z_y_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_y_xx_0[i] = ta1_z_0_xx_0[i] * pa_y[i] - ta1_z_0_xx_1[i] * pc_y[i];

        ta1_z_y_xy_0[i] = ta1_z_0_x_0[i] * fe_0 - ta1_z_0_x_1[i] * fe_0 + ta1_z_0_xy_0[i] * pa_y[i] - ta1_z_0_xy_1[i] * pc_y[i];

        ta1_z_y_xz_0[i] = ta1_z_0_xz_0[i] * pa_y[i] - ta1_z_0_xz_1[i] * pc_y[i];

        ta1_z_y_yy_0[i] = 2.0 * ta1_z_0_y_0[i] * fe_0 - 2.0 * ta1_z_0_y_1[i] * fe_0 + ta1_z_0_yy_0[i] * pa_y[i] - ta1_z_0_yy_1[i] * pc_y[i];

        ta1_z_y_yz_0[i] = ta1_z_0_z_0[i] * fe_0 - ta1_z_0_z_1[i] * fe_0 + ta1_z_0_yz_0[i] * pa_y[i] - ta1_z_0_yz_1[i] * pc_y[i];

        ta1_z_y_zz_0[i] = ta1_z_0_zz_0[i] * pa_y[i] - ta1_z_0_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : PD

    auto ta1_z_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 48);

    auto ta1_z_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 49);

    auto ta1_z_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 50);

    auto ta1_z_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 51);

    auto ta1_z_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 52);

    auto ta1_z_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 53);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_z_0_x_0,  \
                             ta1_z_0_x_1,  \
                             ta1_z_0_xx_0, \
                             ta1_z_0_xx_1, \
                             ta1_z_0_xy_0, \
                             ta1_z_0_xy_1, \
                             ta1_z_0_xz_0, \
                             ta1_z_0_xz_1, \
                             ta1_z_0_y_0,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_yy_0, \
                             ta1_z_0_yy_1, \
                             ta1_z_0_yz_0, \
                             ta1_z_0_yz_1, \
                             ta1_z_0_z_0,  \
                             ta1_z_0_z_1,  \
                             ta1_z_0_zz_0, \
                             ta1_z_0_zz_1, \
                             ta1_z_z_xx_0, \
                             ta1_z_z_xy_0, \
                             ta1_z_z_xz_0, \
                             ta1_z_z_yy_0, \
                             ta1_z_z_yz_0, \
                             ta1_z_z_zz_0, \
                             ta_0_xx_1,    \
                             ta_0_xy_1,    \
                             ta_0_xz_1,    \
                             ta_0_yy_1,    \
                             ta_0_yz_1,    \
                             ta_0_zz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_z_xx_0[i] = ta_0_xx_1[i] + ta1_z_0_xx_0[i] * pa_z[i] - ta1_z_0_xx_1[i] * pc_z[i];

        ta1_z_z_xy_0[i] = ta_0_xy_1[i] + ta1_z_0_xy_0[i] * pa_z[i] - ta1_z_0_xy_1[i] * pc_z[i];

        ta1_z_z_xz_0[i] = ta1_z_0_x_0[i] * fe_0 - ta1_z_0_x_1[i] * fe_0 + ta_0_xz_1[i] + ta1_z_0_xz_0[i] * pa_z[i] - ta1_z_0_xz_1[i] * pc_z[i];

        ta1_z_z_yy_0[i] = ta_0_yy_1[i] + ta1_z_0_yy_0[i] * pa_z[i] - ta1_z_0_yy_1[i] * pc_z[i];

        ta1_z_z_yz_0[i] = ta1_z_0_y_0[i] * fe_0 - ta1_z_0_y_1[i] * fe_0 + ta_0_yz_1[i] + ta1_z_0_yz_0[i] * pa_z[i] - ta1_z_0_yz_1[i] * pc_z[i];

        ta1_z_z_zz_0[i] =
            2.0 * ta1_z_0_z_0[i] * fe_0 - 2.0 * ta1_z_0_z_1[i] * fe_0 + ta_0_zz_1[i] + ta1_z_0_zz_0[i] * pa_z[i] - ta1_z_0_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
