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

#include "NuclearPotentialPrimRecSH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_sh(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_sh,
                               const size_t              idx_npot_0_sf,
                               const size_t              idx_npot_1_sf,
                               const size_t              idx_npot_0_sg,
                               const size_t              idx_npot_1_sg,
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

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_0 = pbuffer.data(idx_npot_0_sf);

    auto ta_0_xyy_0 = pbuffer.data(idx_npot_0_sf + 3);

    auto ta_0_xzz_0 = pbuffer.data(idx_npot_0_sf + 5);

    auto ta_0_yyy_0 = pbuffer.data(idx_npot_0_sf + 6);

    auto ta_0_yzz_0 = pbuffer.data(idx_npot_0_sf + 8);

    auto ta_0_zzz_0 = pbuffer.data(idx_npot_0_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_1 = pbuffer.data(idx_npot_1_sf);

    auto ta_0_xyy_1 = pbuffer.data(idx_npot_1_sf + 3);

    auto ta_0_xzz_1 = pbuffer.data(idx_npot_1_sf + 5);

    auto ta_0_yyy_1 = pbuffer.data(idx_npot_1_sf + 6);

    auto ta_0_yzz_1 = pbuffer.data(idx_npot_1_sf + 8);

    auto ta_0_zzz_1 = pbuffer.data(idx_npot_1_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_0 = pbuffer.data(idx_npot_0_sg);

    auto ta_0_xxxz_0 = pbuffer.data(idx_npot_0_sg + 2);

    auto ta_0_xxyy_0 = pbuffer.data(idx_npot_0_sg + 3);

    auto ta_0_xxzz_0 = pbuffer.data(idx_npot_0_sg + 5);

    auto ta_0_xyyy_0 = pbuffer.data(idx_npot_0_sg + 6);

    auto ta_0_xzzz_0 = pbuffer.data(idx_npot_0_sg + 9);

    auto ta_0_yyyy_0 = pbuffer.data(idx_npot_0_sg + 10);

    auto ta_0_yyyz_0 = pbuffer.data(idx_npot_0_sg + 11);

    auto ta_0_yyzz_0 = pbuffer.data(idx_npot_0_sg + 12);

    auto ta_0_yzzz_0 = pbuffer.data(idx_npot_0_sg + 13);

    auto ta_0_zzzz_0 = pbuffer.data(idx_npot_0_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_1 = pbuffer.data(idx_npot_1_sg);

    auto ta_0_xxxz_1 = pbuffer.data(idx_npot_1_sg + 2);

    auto ta_0_xxyy_1 = pbuffer.data(idx_npot_1_sg + 3);

    auto ta_0_xxzz_1 = pbuffer.data(idx_npot_1_sg + 5);

    auto ta_0_xyyy_1 = pbuffer.data(idx_npot_1_sg + 6);

    auto ta_0_xzzz_1 = pbuffer.data(idx_npot_1_sg + 9);

    auto ta_0_yyyy_1 = pbuffer.data(idx_npot_1_sg + 10);

    auto ta_0_yyyz_1 = pbuffer.data(idx_npot_1_sg + 11);

    auto ta_0_yyzz_1 = pbuffer.data(idx_npot_1_sg + 12);

    auto ta_0_yzzz_1 = pbuffer.data(idx_npot_1_sg + 13);

    auto ta_0_zzzz_1 = pbuffer.data(idx_npot_1_sg + 14);

    // Set up components of targeted buffer : SH

    auto ta_0_xxxxx_0 = pbuffer.data(idx_npot_0_sh);

    auto ta_0_xxxxy_0 = pbuffer.data(idx_npot_0_sh + 1);

    auto ta_0_xxxxz_0 = pbuffer.data(idx_npot_0_sh + 2);

    auto ta_0_xxxyy_0 = pbuffer.data(idx_npot_0_sh + 3);

    auto ta_0_xxxyz_0 = pbuffer.data(idx_npot_0_sh + 4);

    auto ta_0_xxxzz_0 = pbuffer.data(idx_npot_0_sh + 5);

    auto ta_0_xxyyy_0 = pbuffer.data(idx_npot_0_sh + 6);

    auto ta_0_xxyyz_0 = pbuffer.data(idx_npot_0_sh + 7);

    auto ta_0_xxyzz_0 = pbuffer.data(idx_npot_0_sh + 8);

    auto ta_0_xxzzz_0 = pbuffer.data(idx_npot_0_sh + 9);

    auto ta_0_xyyyy_0 = pbuffer.data(idx_npot_0_sh + 10);

    auto ta_0_xyyyz_0 = pbuffer.data(idx_npot_0_sh + 11);

    auto ta_0_xyyzz_0 = pbuffer.data(idx_npot_0_sh + 12);

    auto ta_0_xyzzz_0 = pbuffer.data(idx_npot_0_sh + 13);

    auto ta_0_xzzzz_0 = pbuffer.data(idx_npot_0_sh + 14);

    auto ta_0_yyyyy_0 = pbuffer.data(idx_npot_0_sh + 15);

    auto ta_0_yyyyz_0 = pbuffer.data(idx_npot_0_sh + 16);

    auto ta_0_yyyzz_0 = pbuffer.data(idx_npot_0_sh + 17);

    auto ta_0_yyzzz_0 = pbuffer.data(idx_npot_0_sh + 18);

    auto ta_0_yzzzz_0 = pbuffer.data(idx_npot_0_sh + 19);

    auto ta_0_zzzzz_0 = pbuffer.data(idx_npot_0_sh + 20);

#pragma omp simd aligned(pb_x,             \
                             pb_y,         \
                             pb_z,         \
                             pc_x,         \
                             pc_y,         \
                             pc_z,         \
                             ta_0_xxx_0,   \
                             ta_0_xxx_1,   \
                             ta_0_xxxx_0,  \
                             ta_0_xxxx_1,  \
                             ta_0_xxxxx_0, \
                             ta_0_xxxxy_0, \
                             ta_0_xxxxz_0, \
                             ta_0_xxxyy_0, \
                             ta_0_xxxyz_0, \
                             ta_0_xxxz_0,  \
                             ta_0_xxxz_1,  \
                             ta_0_xxxzz_0, \
                             ta_0_xxyy_0,  \
                             ta_0_xxyy_1,  \
                             ta_0_xxyyy_0, \
                             ta_0_xxyyz_0, \
                             ta_0_xxyzz_0, \
                             ta_0_xxzz_0,  \
                             ta_0_xxzz_1,  \
                             ta_0_xxzzz_0, \
                             ta_0_xyy_0,   \
                             ta_0_xyy_1,   \
                             ta_0_xyyy_0,  \
                             ta_0_xyyy_1,  \
                             ta_0_xyyyy_0, \
                             ta_0_xyyyz_0, \
                             ta_0_xyyzz_0, \
                             ta_0_xyzzz_0, \
                             ta_0_xzz_0,   \
                             ta_0_xzz_1,   \
                             ta_0_xzzz_0,  \
                             ta_0_xzzz_1,  \
                             ta_0_xzzzz_0, \
                             ta_0_yyy_0,   \
                             ta_0_yyy_1,   \
                             ta_0_yyyy_0,  \
                             ta_0_yyyy_1,  \
                             ta_0_yyyyy_0, \
                             ta_0_yyyyz_0, \
                             ta_0_yyyz_0,  \
                             ta_0_yyyz_1,  \
                             ta_0_yyyzz_0, \
                             ta_0_yyzz_0,  \
                             ta_0_yyzz_1,  \
                             ta_0_yyzzz_0, \
                             ta_0_yzz_0,   \
                             ta_0_yzz_1,   \
                             ta_0_yzzz_0,  \
                             ta_0_yzzz_1,  \
                             ta_0_yzzzz_0, \
                             ta_0_zzz_0,   \
                             ta_0_zzz_1,   \
                             ta_0_zzzz_0,  \
                             ta_0_zzzz_1,  \
                             ta_0_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_0_xxxxx_0[i] = 4.0 * ta_0_xxx_0[i] * fe_0 - 4.0 * ta_0_xxx_1[i] * fe_0 + ta_0_xxxx_0[i] * pb_x[i] - ta_0_xxxx_1[i] * pc_x[i];

        ta_0_xxxxy_0[i] = ta_0_xxxx_0[i] * pb_y[i] - ta_0_xxxx_1[i] * pc_y[i];

        ta_0_xxxxz_0[i] = ta_0_xxxx_0[i] * pb_z[i] - ta_0_xxxx_1[i] * pc_z[i];

        ta_0_xxxyy_0[i] = 2.0 * ta_0_xyy_0[i] * fe_0 - 2.0 * ta_0_xyy_1[i] * fe_0 + ta_0_xxyy_0[i] * pb_x[i] - ta_0_xxyy_1[i] * pc_x[i];

        ta_0_xxxyz_0[i] = ta_0_xxxz_0[i] * pb_y[i] - ta_0_xxxz_1[i] * pc_y[i];

        ta_0_xxxzz_0[i] = 2.0 * ta_0_xzz_0[i] * fe_0 - 2.0 * ta_0_xzz_1[i] * fe_0 + ta_0_xxzz_0[i] * pb_x[i] - ta_0_xxzz_1[i] * pc_x[i];

        ta_0_xxyyy_0[i] = ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + ta_0_xyyy_0[i] * pb_x[i] - ta_0_xyyy_1[i] * pc_x[i];

        ta_0_xxyyz_0[i] = ta_0_xxyy_0[i] * pb_z[i] - ta_0_xxyy_1[i] * pc_z[i];

        ta_0_xxyzz_0[i] = ta_0_xxzz_0[i] * pb_y[i] - ta_0_xxzz_1[i] * pc_y[i];

        ta_0_xxzzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_0_xzzz_0[i] * pb_x[i] - ta_0_xzzz_1[i] * pc_x[i];

        ta_0_xyyyy_0[i] = ta_0_yyyy_0[i] * pb_x[i] - ta_0_yyyy_1[i] * pc_x[i];

        ta_0_xyyyz_0[i] = ta_0_yyyz_0[i] * pb_x[i] - ta_0_yyyz_1[i] * pc_x[i];

        ta_0_xyyzz_0[i] = ta_0_yyzz_0[i] * pb_x[i] - ta_0_yyzz_1[i] * pc_x[i];

        ta_0_xyzzz_0[i] = ta_0_yzzz_0[i] * pb_x[i] - ta_0_yzzz_1[i] * pc_x[i];

        ta_0_xzzzz_0[i] = ta_0_zzzz_0[i] * pb_x[i] - ta_0_zzzz_1[i] * pc_x[i];

        ta_0_yyyyy_0[i] = 4.0 * ta_0_yyy_0[i] * fe_0 - 4.0 * ta_0_yyy_1[i] * fe_0 + ta_0_yyyy_0[i] * pb_y[i] - ta_0_yyyy_1[i] * pc_y[i];

        ta_0_yyyyz_0[i] = ta_0_yyyy_0[i] * pb_z[i] - ta_0_yyyy_1[i] * pc_z[i];

        ta_0_yyyzz_0[i] = 2.0 * ta_0_yzz_0[i] * fe_0 - 2.0 * ta_0_yzz_1[i] * fe_0 + ta_0_yyzz_0[i] * pb_y[i] - ta_0_yyzz_1[i] * pc_y[i];

        ta_0_yyzzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_0_yzzz_0[i] * pb_y[i] - ta_0_yzzz_1[i] * pc_y[i];

        ta_0_yzzzz_0[i] = ta_0_zzzz_0[i] * pb_y[i] - ta_0_zzzz_1[i] * pc_y[i];

        ta_0_zzzzz_0[i] = 4.0 * ta_0_zzz_0[i] * fe_0 - 4.0 * ta_0_zzz_1[i] * fe_0 + ta_0_zzzz_0[i] * pb_z[i] - ta_0_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
