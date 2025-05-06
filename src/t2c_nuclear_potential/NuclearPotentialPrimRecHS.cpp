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

#include "NuclearPotentialPrimRecHS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_hs(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_hs,
                               const size_t              idx_npot_0_fs,
                               const size_t              idx_npot_1_fs,
                               const size_t              idx_npot_0_gs,
                               const size_t              idx_npot_1_gs,
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

    // Set up components of auxiliary buffer : FS

    auto ta_xxx_0_0 = pbuffer.data(idx_npot_0_fs);

    auto ta_xyy_0_0 = pbuffer.data(idx_npot_0_fs + 3);

    auto ta_xzz_0_0 = pbuffer.data(idx_npot_0_fs + 5);

    auto ta_yyy_0_0 = pbuffer.data(idx_npot_0_fs + 6);

    auto ta_yzz_0_0 = pbuffer.data(idx_npot_0_fs + 8);

    auto ta_zzz_0_0 = pbuffer.data(idx_npot_0_fs + 9);

    // Set up components of auxiliary buffer : FS

    auto ta_xxx_0_1 = pbuffer.data(idx_npot_1_fs);

    auto ta_xyy_0_1 = pbuffer.data(idx_npot_1_fs + 3);

    auto ta_xzz_0_1 = pbuffer.data(idx_npot_1_fs + 5);

    auto ta_yyy_0_1 = pbuffer.data(idx_npot_1_fs + 6);

    auto ta_yzz_0_1 = pbuffer.data(idx_npot_1_fs + 8);

    auto ta_zzz_0_1 = pbuffer.data(idx_npot_1_fs + 9);

    // Set up components of auxiliary buffer : GS

    auto ta_xxxx_0_0 = pbuffer.data(idx_npot_0_gs);

    auto ta_xxxz_0_0 = pbuffer.data(idx_npot_0_gs + 2);

    auto ta_xxyy_0_0 = pbuffer.data(idx_npot_0_gs + 3);

    auto ta_xxzz_0_0 = pbuffer.data(idx_npot_0_gs + 5);

    auto ta_xyyy_0_0 = pbuffer.data(idx_npot_0_gs + 6);

    auto ta_xzzz_0_0 = pbuffer.data(idx_npot_0_gs + 9);

    auto ta_yyyy_0_0 = pbuffer.data(idx_npot_0_gs + 10);

    auto ta_yyyz_0_0 = pbuffer.data(idx_npot_0_gs + 11);

    auto ta_yyzz_0_0 = pbuffer.data(idx_npot_0_gs + 12);

    auto ta_yzzz_0_0 = pbuffer.data(idx_npot_0_gs + 13);

    auto ta_zzzz_0_0 = pbuffer.data(idx_npot_0_gs + 14);

    // Set up components of auxiliary buffer : GS

    auto ta_xxxx_0_1 = pbuffer.data(idx_npot_1_gs);

    auto ta_xxxz_0_1 = pbuffer.data(idx_npot_1_gs + 2);

    auto ta_xxyy_0_1 = pbuffer.data(idx_npot_1_gs + 3);

    auto ta_xxzz_0_1 = pbuffer.data(idx_npot_1_gs + 5);

    auto ta_xyyy_0_1 = pbuffer.data(idx_npot_1_gs + 6);

    auto ta_xzzz_0_1 = pbuffer.data(idx_npot_1_gs + 9);

    auto ta_yyyy_0_1 = pbuffer.data(idx_npot_1_gs + 10);

    auto ta_yyyz_0_1 = pbuffer.data(idx_npot_1_gs + 11);

    auto ta_yyzz_0_1 = pbuffer.data(idx_npot_1_gs + 12);

    auto ta_yzzz_0_1 = pbuffer.data(idx_npot_1_gs + 13);

    auto ta_zzzz_0_1 = pbuffer.data(idx_npot_1_gs + 14);

    // Set up components of targeted buffer : HS

    auto ta_xxxxx_0_0 = pbuffer.data(idx_npot_0_hs);

    auto ta_xxxxy_0_0 = pbuffer.data(idx_npot_0_hs + 1);

    auto ta_xxxxz_0_0 = pbuffer.data(idx_npot_0_hs + 2);

    auto ta_xxxyy_0_0 = pbuffer.data(idx_npot_0_hs + 3);

    auto ta_xxxyz_0_0 = pbuffer.data(idx_npot_0_hs + 4);

    auto ta_xxxzz_0_0 = pbuffer.data(idx_npot_0_hs + 5);

    auto ta_xxyyy_0_0 = pbuffer.data(idx_npot_0_hs + 6);

    auto ta_xxyyz_0_0 = pbuffer.data(idx_npot_0_hs + 7);

    auto ta_xxyzz_0_0 = pbuffer.data(idx_npot_0_hs + 8);

    auto ta_xxzzz_0_0 = pbuffer.data(idx_npot_0_hs + 9);

    auto ta_xyyyy_0_0 = pbuffer.data(idx_npot_0_hs + 10);

    auto ta_xyyyz_0_0 = pbuffer.data(idx_npot_0_hs + 11);

    auto ta_xyyzz_0_0 = pbuffer.data(idx_npot_0_hs + 12);

    auto ta_xyzzz_0_0 = pbuffer.data(idx_npot_0_hs + 13);

    auto ta_xzzzz_0_0 = pbuffer.data(idx_npot_0_hs + 14);

    auto ta_yyyyy_0_0 = pbuffer.data(idx_npot_0_hs + 15);

    auto ta_yyyyz_0_0 = pbuffer.data(idx_npot_0_hs + 16);

    auto ta_yyyzz_0_0 = pbuffer.data(idx_npot_0_hs + 17);

    auto ta_yyzzz_0_0 = pbuffer.data(idx_npot_0_hs + 18);

    auto ta_yzzzz_0_0 = pbuffer.data(idx_npot_0_hs + 19);

    auto ta_zzzzz_0_0 = pbuffer.data(idx_npot_0_hs + 20);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             pc_x,         \
                             pc_y,         \
                             pc_z,         \
                             ta_xxx_0_0,   \
                             ta_xxx_0_1,   \
                             ta_xxxx_0_0,  \
                             ta_xxxx_0_1,  \
                             ta_xxxxx_0_0, \
                             ta_xxxxy_0_0, \
                             ta_xxxxz_0_0, \
                             ta_xxxyy_0_0, \
                             ta_xxxyz_0_0, \
                             ta_xxxz_0_0,  \
                             ta_xxxz_0_1,  \
                             ta_xxxzz_0_0, \
                             ta_xxyy_0_0,  \
                             ta_xxyy_0_1,  \
                             ta_xxyyy_0_0, \
                             ta_xxyyz_0_0, \
                             ta_xxyzz_0_0, \
                             ta_xxzz_0_0,  \
                             ta_xxzz_0_1,  \
                             ta_xxzzz_0_0, \
                             ta_xyy_0_0,   \
                             ta_xyy_0_1,   \
                             ta_xyyy_0_0,  \
                             ta_xyyy_0_1,  \
                             ta_xyyyy_0_0, \
                             ta_xyyyz_0_0, \
                             ta_xyyzz_0_0, \
                             ta_xyzzz_0_0, \
                             ta_xzz_0_0,   \
                             ta_xzz_0_1,   \
                             ta_xzzz_0_0,  \
                             ta_xzzz_0_1,  \
                             ta_xzzzz_0_0, \
                             ta_yyy_0_0,   \
                             ta_yyy_0_1,   \
                             ta_yyyy_0_0,  \
                             ta_yyyy_0_1,  \
                             ta_yyyyy_0_0, \
                             ta_yyyyz_0_0, \
                             ta_yyyz_0_0,  \
                             ta_yyyz_0_1,  \
                             ta_yyyzz_0_0, \
                             ta_yyzz_0_0,  \
                             ta_yyzz_0_1,  \
                             ta_yyzzz_0_0, \
                             ta_yzz_0_0,   \
                             ta_yzz_0_1,   \
                             ta_yzzz_0_0,  \
                             ta_yzzz_0_1,  \
                             ta_yzzzz_0_0, \
                             ta_zzz_0_0,   \
                             ta_zzz_0_1,   \
                             ta_zzzz_0_0,  \
                             ta_zzzz_0_1,  \
                             ta_zzzzz_0_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxx_0_0[i] = 4.0 * ta_xxx_0_0[i] * fe_0 - 4.0 * ta_xxx_0_1[i] * fe_0 + ta_xxxx_0_0[i] * pa_x[i] - ta_xxxx_0_1[i] * pc_x[i];

        ta_xxxxy_0_0[i] = ta_xxxx_0_0[i] * pa_y[i] - ta_xxxx_0_1[i] * pc_y[i];

        ta_xxxxz_0_0[i] = ta_xxxx_0_0[i] * pa_z[i] - ta_xxxx_0_1[i] * pc_z[i];

        ta_xxxyy_0_0[i] = 2.0 * ta_xyy_0_0[i] * fe_0 - 2.0 * ta_xyy_0_1[i] * fe_0 + ta_xxyy_0_0[i] * pa_x[i] - ta_xxyy_0_1[i] * pc_x[i];

        ta_xxxyz_0_0[i] = ta_xxxz_0_0[i] * pa_y[i] - ta_xxxz_0_1[i] * pc_y[i];

        ta_xxxzz_0_0[i] = 2.0 * ta_xzz_0_0[i] * fe_0 - 2.0 * ta_xzz_0_1[i] * fe_0 + ta_xxzz_0_0[i] * pa_x[i] - ta_xxzz_0_1[i] * pc_x[i];

        ta_xxyyy_0_0[i] = ta_yyy_0_0[i] * fe_0 - ta_yyy_0_1[i] * fe_0 + ta_xyyy_0_0[i] * pa_x[i] - ta_xyyy_0_1[i] * pc_x[i];

        ta_xxyyz_0_0[i] = ta_xxyy_0_0[i] * pa_z[i] - ta_xxyy_0_1[i] * pc_z[i];

        ta_xxyzz_0_0[i] = ta_xxzz_0_0[i] * pa_y[i] - ta_xxzz_0_1[i] * pc_y[i];

        ta_xxzzz_0_0[i] = ta_zzz_0_0[i] * fe_0 - ta_zzz_0_1[i] * fe_0 + ta_xzzz_0_0[i] * pa_x[i] - ta_xzzz_0_1[i] * pc_x[i];

        ta_xyyyy_0_0[i] = ta_yyyy_0_0[i] * pa_x[i] - ta_yyyy_0_1[i] * pc_x[i];

        ta_xyyyz_0_0[i] = ta_yyyz_0_0[i] * pa_x[i] - ta_yyyz_0_1[i] * pc_x[i];

        ta_xyyzz_0_0[i] = ta_yyzz_0_0[i] * pa_x[i] - ta_yyzz_0_1[i] * pc_x[i];

        ta_xyzzz_0_0[i] = ta_yzzz_0_0[i] * pa_x[i] - ta_yzzz_0_1[i] * pc_x[i];

        ta_xzzzz_0_0[i] = ta_zzzz_0_0[i] * pa_x[i] - ta_zzzz_0_1[i] * pc_x[i];

        ta_yyyyy_0_0[i] = 4.0 * ta_yyy_0_0[i] * fe_0 - 4.0 * ta_yyy_0_1[i] * fe_0 + ta_yyyy_0_0[i] * pa_y[i] - ta_yyyy_0_1[i] * pc_y[i];

        ta_yyyyz_0_0[i] = ta_yyyy_0_0[i] * pa_z[i] - ta_yyyy_0_1[i] * pc_z[i];

        ta_yyyzz_0_0[i] = 2.0 * ta_yzz_0_0[i] * fe_0 - 2.0 * ta_yzz_0_1[i] * fe_0 + ta_yyzz_0_0[i] * pa_y[i] - ta_yyzz_0_1[i] * pc_y[i];

        ta_yyzzz_0_0[i] = ta_zzz_0_0[i] * fe_0 - ta_zzz_0_1[i] * fe_0 + ta_yzzz_0_0[i] * pa_y[i] - ta_yzzz_0_1[i] * pc_y[i];

        ta_yzzzz_0_0[i] = ta_zzzz_0_0[i] * pa_y[i] - ta_zzzz_0_1[i] * pc_y[i];

        ta_zzzzz_0_0[i] = 4.0 * ta_zzz_0_0[i] * fe_0 - 4.0 * ta_zzz_0_1[i] * fe_0 + ta_zzzz_0_0[i] * pa_z[i] - ta_zzzz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
