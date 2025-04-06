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

#include "NuclearPotentialPrimRecIP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_ip(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_ip,
                               const size_t              idx_npot_0_gp,
                               const size_t              idx_npot_1_gp,
                               const size_t              idx_npot_0_hs,
                               const size_t              idx_npot_1_hs,
                               const size_t              idx_npot_0_hp,
                               const size_t              idx_npot_1_hp,
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

    // Set up components of auxiliary buffer : GP

    auto ta_xxxx_x_0 = pbuffer.data(idx_npot_0_gp);

    auto ta_xxxx_y_0 = pbuffer.data(idx_npot_0_gp + 1);

    auto ta_xxxx_z_0 = pbuffer.data(idx_npot_0_gp + 2);

    auto ta_xxxy_x_0 = pbuffer.data(idx_npot_0_gp + 3);

    auto ta_xxxy_y_0 = pbuffer.data(idx_npot_0_gp + 4);

    auto ta_xxxz_x_0 = pbuffer.data(idx_npot_0_gp + 6);

    auto ta_xxxz_z_0 = pbuffer.data(idx_npot_0_gp + 8);

    auto ta_xxyy_x_0 = pbuffer.data(idx_npot_0_gp + 9);

    auto ta_xxyy_y_0 = pbuffer.data(idx_npot_0_gp + 10);

    auto ta_xxyy_z_0 = pbuffer.data(idx_npot_0_gp + 11);

    auto ta_xxzz_x_0 = pbuffer.data(idx_npot_0_gp + 15);

    auto ta_xxzz_y_0 = pbuffer.data(idx_npot_0_gp + 16);

    auto ta_xxzz_z_0 = pbuffer.data(idx_npot_0_gp + 17);

    auto ta_xyyy_y_0 = pbuffer.data(idx_npot_0_gp + 19);

    auto ta_xyyy_z_0 = pbuffer.data(idx_npot_0_gp + 20);

    auto ta_xyyz_z_0 = pbuffer.data(idx_npot_0_gp + 23);

    auto ta_xyzz_y_0 = pbuffer.data(idx_npot_0_gp + 25);

    auto ta_xzzz_y_0 = pbuffer.data(idx_npot_0_gp + 28);

    auto ta_xzzz_z_0 = pbuffer.data(idx_npot_0_gp + 29);

    auto ta_yyyy_x_0 = pbuffer.data(idx_npot_0_gp + 30);

    auto ta_yyyy_y_0 = pbuffer.data(idx_npot_0_gp + 31);

    auto ta_yyyy_z_0 = pbuffer.data(idx_npot_0_gp + 32);

    auto ta_yyyz_y_0 = pbuffer.data(idx_npot_0_gp + 34);

    auto ta_yyyz_z_0 = pbuffer.data(idx_npot_0_gp + 35);

    auto ta_yyzz_x_0 = pbuffer.data(idx_npot_0_gp + 36);

    auto ta_yyzz_y_0 = pbuffer.data(idx_npot_0_gp + 37);

    auto ta_yyzz_z_0 = pbuffer.data(idx_npot_0_gp + 38);

    auto ta_yzzz_x_0 = pbuffer.data(idx_npot_0_gp + 39);

    auto ta_yzzz_y_0 = pbuffer.data(idx_npot_0_gp + 40);

    auto ta_yzzz_z_0 = pbuffer.data(idx_npot_0_gp + 41);

    auto ta_zzzz_x_0 = pbuffer.data(idx_npot_0_gp + 42);

    auto ta_zzzz_y_0 = pbuffer.data(idx_npot_0_gp + 43);

    auto ta_zzzz_z_0 = pbuffer.data(idx_npot_0_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto ta_xxxx_x_1 = pbuffer.data(idx_npot_1_gp);

    auto ta_xxxx_y_1 = pbuffer.data(idx_npot_1_gp + 1);

    auto ta_xxxx_z_1 = pbuffer.data(idx_npot_1_gp + 2);

    auto ta_xxxy_x_1 = pbuffer.data(idx_npot_1_gp + 3);

    auto ta_xxxy_y_1 = pbuffer.data(idx_npot_1_gp + 4);

    auto ta_xxxz_x_1 = pbuffer.data(idx_npot_1_gp + 6);

    auto ta_xxxz_z_1 = pbuffer.data(idx_npot_1_gp + 8);

    auto ta_xxyy_x_1 = pbuffer.data(idx_npot_1_gp + 9);

    auto ta_xxyy_y_1 = pbuffer.data(idx_npot_1_gp + 10);

    auto ta_xxyy_z_1 = pbuffer.data(idx_npot_1_gp + 11);

    auto ta_xxzz_x_1 = pbuffer.data(idx_npot_1_gp + 15);

    auto ta_xxzz_y_1 = pbuffer.data(idx_npot_1_gp + 16);

    auto ta_xxzz_z_1 = pbuffer.data(idx_npot_1_gp + 17);

    auto ta_xyyy_y_1 = pbuffer.data(idx_npot_1_gp + 19);

    auto ta_xyyy_z_1 = pbuffer.data(idx_npot_1_gp + 20);

    auto ta_xyyz_z_1 = pbuffer.data(idx_npot_1_gp + 23);

    auto ta_xyzz_y_1 = pbuffer.data(idx_npot_1_gp + 25);

    auto ta_xzzz_y_1 = pbuffer.data(idx_npot_1_gp + 28);

    auto ta_xzzz_z_1 = pbuffer.data(idx_npot_1_gp + 29);

    auto ta_yyyy_x_1 = pbuffer.data(idx_npot_1_gp + 30);

    auto ta_yyyy_y_1 = pbuffer.data(idx_npot_1_gp + 31);

    auto ta_yyyy_z_1 = pbuffer.data(idx_npot_1_gp + 32);

    auto ta_yyyz_y_1 = pbuffer.data(idx_npot_1_gp + 34);

    auto ta_yyyz_z_1 = pbuffer.data(idx_npot_1_gp + 35);

    auto ta_yyzz_x_1 = pbuffer.data(idx_npot_1_gp + 36);

    auto ta_yyzz_y_1 = pbuffer.data(idx_npot_1_gp + 37);

    auto ta_yyzz_z_1 = pbuffer.data(idx_npot_1_gp + 38);

    auto ta_yzzz_x_1 = pbuffer.data(idx_npot_1_gp + 39);

    auto ta_yzzz_y_1 = pbuffer.data(idx_npot_1_gp + 40);

    auto ta_yzzz_z_1 = pbuffer.data(idx_npot_1_gp + 41);

    auto ta_zzzz_x_1 = pbuffer.data(idx_npot_1_gp + 42);

    auto ta_zzzz_y_1 = pbuffer.data(idx_npot_1_gp + 43);

    auto ta_zzzz_z_1 = pbuffer.data(idx_npot_1_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto ta_xxxxx_0_0 = pbuffer.data(idx_npot_0_hs);

    auto ta_yyyyy_0_0 = pbuffer.data(idx_npot_0_hs + 15);

    auto ta_yyyzz_0_0 = pbuffer.data(idx_npot_0_hs + 17);

    auto ta_yyzzz_0_0 = pbuffer.data(idx_npot_0_hs + 18);

    auto ta_zzzzz_0_0 = pbuffer.data(idx_npot_0_hs + 20);

    // Set up components of auxiliary buffer : HS

    auto ta_xxxxx_0_1 = pbuffer.data(idx_npot_1_hs);

    auto ta_yyyyy_0_1 = pbuffer.data(idx_npot_1_hs + 15);

    auto ta_yyyzz_0_1 = pbuffer.data(idx_npot_1_hs + 17);

    auto ta_yyzzz_0_1 = pbuffer.data(idx_npot_1_hs + 18);

    auto ta_zzzzz_0_1 = pbuffer.data(idx_npot_1_hs + 20);

    // Set up components of auxiliary buffer : HP

    auto ta_xxxxx_x_0 = pbuffer.data(idx_npot_0_hp);

    auto ta_xxxxx_y_0 = pbuffer.data(idx_npot_0_hp + 1);

    auto ta_xxxxx_z_0 = pbuffer.data(idx_npot_0_hp + 2);

    auto ta_xxxxy_x_0 = pbuffer.data(idx_npot_0_hp + 3);

    auto ta_xxxxy_y_0 = pbuffer.data(idx_npot_0_hp + 4);

    auto ta_xxxxz_x_0 = pbuffer.data(idx_npot_0_hp + 6);

    auto ta_xxxxz_z_0 = pbuffer.data(idx_npot_0_hp + 8);

    auto ta_xxxyy_x_0 = pbuffer.data(idx_npot_0_hp + 9);

    auto ta_xxxyy_y_0 = pbuffer.data(idx_npot_0_hp + 10);

    auto ta_xxxyy_z_0 = pbuffer.data(idx_npot_0_hp + 11);

    auto ta_xxxzz_x_0 = pbuffer.data(idx_npot_0_hp + 15);

    auto ta_xxxzz_y_0 = pbuffer.data(idx_npot_0_hp + 16);

    auto ta_xxxzz_z_0 = pbuffer.data(idx_npot_0_hp + 17);

    auto ta_xxyyy_x_0 = pbuffer.data(idx_npot_0_hp + 18);

    auto ta_xxyyy_y_0 = pbuffer.data(idx_npot_0_hp + 19);

    auto ta_xxyyy_z_0 = pbuffer.data(idx_npot_0_hp + 20);

    auto ta_xxyyz_z_0 = pbuffer.data(idx_npot_0_hp + 23);

    auto ta_xxyzz_x_0 = pbuffer.data(idx_npot_0_hp + 24);

    auto ta_xxyzz_y_0 = pbuffer.data(idx_npot_0_hp + 25);

    auto ta_xxzzz_x_0 = pbuffer.data(idx_npot_0_hp + 27);

    auto ta_xxzzz_y_0 = pbuffer.data(idx_npot_0_hp + 28);

    auto ta_xxzzz_z_0 = pbuffer.data(idx_npot_0_hp + 29);

    auto ta_xyyyy_x_0 = pbuffer.data(idx_npot_0_hp + 30);

    auto ta_xyyyy_y_0 = pbuffer.data(idx_npot_0_hp + 31);

    auto ta_xyyyy_z_0 = pbuffer.data(idx_npot_0_hp + 32);

    auto ta_xyyyz_z_0 = pbuffer.data(idx_npot_0_hp + 35);

    auto ta_xyyzz_y_0 = pbuffer.data(idx_npot_0_hp + 37);

    auto ta_xyyzz_z_0 = pbuffer.data(idx_npot_0_hp + 38);

    auto ta_xyzzz_y_0 = pbuffer.data(idx_npot_0_hp + 40);

    auto ta_xzzzz_x_0 = pbuffer.data(idx_npot_0_hp + 42);

    auto ta_xzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 43);

    auto ta_xzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 44);

    auto ta_yyyyy_x_0 = pbuffer.data(idx_npot_0_hp + 45);

    auto ta_yyyyy_y_0 = pbuffer.data(idx_npot_0_hp + 46);

    auto ta_yyyyy_z_0 = pbuffer.data(idx_npot_0_hp + 47);

    auto ta_yyyyz_y_0 = pbuffer.data(idx_npot_0_hp + 49);

    auto ta_yyyyz_z_0 = pbuffer.data(idx_npot_0_hp + 50);

    auto ta_yyyzz_x_0 = pbuffer.data(idx_npot_0_hp + 51);

    auto ta_yyyzz_y_0 = pbuffer.data(idx_npot_0_hp + 52);

    auto ta_yyyzz_z_0 = pbuffer.data(idx_npot_0_hp + 53);

    auto ta_yyzzz_x_0 = pbuffer.data(idx_npot_0_hp + 54);

    auto ta_yyzzz_y_0 = pbuffer.data(idx_npot_0_hp + 55);

    auto ta_yyzzz_z_0 = pbuffer.data(idx_npot_0_hp + 56);

    auto ta_yzzzz_x_0 = pbuffer.data(idx_npot_0_hp + 57);

    auto ta_yzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 58);

    auto ta_yzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 59);

    auto ta_zzzzz_x_0 = pbuffer.data(idx_npot_0_hp + 60);

    auto ta_zzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 61);

    auto ta_zzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 62);

    // Set up components of auxiliary buffer : HP

    auto ta_xxxxx_x_1 = pbuffer.data(idx_npot_1_hp);

    auto ta_xxxxx_y_1 = pbuffer.data(idx_npot_1_hp + 1);

    auto ta_xxxxx_z_1 = pbuffer.data(idx_npot_1_hp + 2);

    auto ta_xxxxy_x_1 = pbuffer.data(idx_npot_1_hp + 3);

    auto ta_xxxxy_y_1 = pbuffer.data(idx_npot_1_hp + 4);

    auto ta_xxxxz_x_1 = pbuffer.data(idx_npot_1_hp + 6);

    auto ta_xxxxz_z_1 = pbuffer.data(idx_npot_1_hp + 8);

    auto ta_xxxyy_x_1 = pbuffer.data(idx_npot_1_hp + 9);

    auto ta_xxxyy_y_1 = pbuffer.data(idx_npot_1_hp + 10);

    auto ta_xxxyy_z_1 = pbuffer.data(idx_npot_1_hp + 11);

    auto ta_xxxzz_x_1 = pbuffer.data(idx_npot_1_hp + 15);

    auto ta_xxxzz_y_1 = pbuffer.data(idx_npot_1_hp + 16);

    auto ta_xxxzz_z_1 = pbuffer.data(idx_npot_1_hp + 17);

    auto ta_xxyyy_x_1 = pbuffer.data(idx_npot_1_hp + 18);

    auto ta_xxyyy_y_1 = pbuffer.data(idx_npot_1_hp + 19);

    auto ta_xxyyy_z_1 = pbuffer.data(idx_npot_1_hp + 20);

    auto ta_xxyyz_z_1 = pbuffer.data(idx_npot_1_hp + 23);

    auto ta_xxyzz_x_1 = pbuffer.data(idx_npot_1_hp + 24);

    auto ta_xxyzz_y_1 = pbuffer.data(idx_npot_1_hp + 25);

    auto ta_xxzzz_x_1 = pbuffer.data(idx_npot_1_hp + 27);

    auto ta_xxzzz_y_1 = pbuffer.data(idx_npot_1_hp + 28);

    auto ta_xxzzz_z_1 = pbuffer.data(idx_npot_1_hp + 29);

    auto ta_xyyyy_x_1 = pbuffer.data(idx_npot_1_hp + 30);

    auto ta_xyyyy_y_1 = pbuffer.data(idx_npot_1_hp + 31);

    auto ta_xyyyy_z_1 = pbuffer.data(idx_npot_1_hp + 32);

    auto ta_xyyyz_z_1 = pbuffer.data(idx_npot_1_hp + 35);

    auto ta_xyyzz_y_1 = pbuffer.data(idx_npot_1_hp + 37);

    auto ta_xyyzz_z_1 = pbuffer.data(idx_npot_1_hp + 38);

    auto ta_xyzzz_y_1 = pbuffer.data(idx_npot_1_hp + 40);

    auto ta_xzzzz_x_1 = pbuffer.data(idx_npot_1_hp + 42);

    auto ta_xzzzz_y_1 = pbuffer.data(idx_npot_1_hp + 43);

    auto ta_xzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 44);

    auto ta_yyyyy_x_1 = pbuffer.data(idx_npot_1_hp + 45);

    auto ta_yyyyy_y_1 = pbuffer.data(idx_npot_1_hp + 46);

    auto ta_yyyyy_z_1 = pbuffer.data(idx_npot_1_hp + 47);

    auto ta_yyyyz_y_1 = pbuffer.data(idx_npot_1_hp + 49);

    auto ta_yyyyz_z_1 = pbuffer.data(idx_npot_1_hp + 50);

    auto ta_yyyzz_x_1 = pbuffer.data(idx_npot_1_hp + 51);

    auto ta_yyyzz_y_1 = pbuffer.data(idx_npot_1_hp + 52);

    auto ta_yyyzz_z_1 = pbuffer.data(idx_npot_1_hp + 53);

    auto ta_yyzzz_x_1 = pbuffer.data(idx_npot_1_hp + 54);

    auto ta_yyzzz_y_1 = pbuffer.data(idx_npot_1_hp + 55);

    auto ta_yyzzz_z_1 = pbuffer.data(idx_npot_1_hp + 56);

    auto ta_yzzzz_x_1 = pbuffer.data(idx_npot_1_hp + 57);

    auto ta_yzzzz_y_1 = pbuffer.data(idx_npot_1_hp + 58);

    auto ta_yzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 59);

    auto ta_zzzzz_x_1 = pbuffer.data(idx_npot_1_hp + 60);

    auto ta_zzzzz_y_1 = pbuffer.data(idx_npot_1_hp + 61);

    auto ta_zzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 62);

    // Set up 0-3 components of targeted buffer : IP

    auto ta_xxxxxx_x_0 = pbuffer.data(idx_npot_0_ip);

    auto ta_xxxxxx_y_0 = pbuffer.data(idx_npot_0_ip + 1);

    auto ta_xxxxxx_z_0 = pbuffer.data(idx_npot_0_ip + 2);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xxxx_x_0,   \
                             ta_xxxx_x_1,   \
                             ta_xxxx_y_0,   \
                             ta_xxxx_y_1,   \
                             ta_xxxx_z_0,   \
                             ta_xxxx_z_1,   \
                             ta_xxxxx_0_0,  \
                             ta_xxxxx_0_1,  \
                             ta_xxxxx_x_0,  \
                             ta_xxxxx_x_1,  \
                             ta_xxxxx_y_0,  \
                             ta_xxxxx_y_1,  \
                             ta_xxxxx_z_0,  \
                             ta_xxxxx_z_1,  \
                             ta_xxxxxx_x_0, \
                             ta_xxxxxx_y_0, \
                             ta_xxxxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxx_x_0[i] = 5.0 * ta_xxxx_x_0[i] * fe_0 - 5.0 * ta_xxxx_x_1[i] * fe_0 + ta_xxxxx_0_0[i] * fe_0 - ta_xxxxx_0_1[i] * fe_0 +
                           ta_xxxxx_x_0[i] * pa_x[i] - ta_xxxxx_x_1[i] * pc_x[i];

        ta_xxxxxx_y_0[i] = 5.0 * ta_xxxx_y_0[i] * fe_0 - 5.0 * ta_xxxx_y_1[i] * fe_0 + ta_xxxxx_y_0[i] * pa_x[i] - ta_xxxxx_y_1[i] * pc_x[i];

        ta_xxxxxx_z_0[i] = 5.0 * ta_xxxx_z_0[i] * fe_0 - 5.0 * ta_xxxx_z_1[i] * fe_0 + ta_xxxxx_z_0[i] * pa_x[i] - ta_xxxxx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : IP

    auto ta_xxxxxy_x_0 = pbuffer.data(idx_npot_0_ip + 3);

    auto ta_xxxxxy_y_0 = pbuffer.data(idx_npot_0_ip + 4);

    auto ta_xxxxxy_z_0 = pbuffer.data(idx_npot_0_ip + 5);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxxxx_x_0,  \
                             ta_xxxxx_x_1,  \
                             ta_xxxxx_z_0,  \
                             ta_xxxxx_z_1,  \
                             ta_xxxxxy_x_0, \
                             ta_xxxxxy_y_0, \
                             ta_xxxxxy_z_0, \
                             ta_xxxxy_y_0,  \
                             ta_xxxxy_y_1,  \
                             ta_xxxy_y_0,   \
                             ta_xxxy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxy_x_0[i] = ta_xxxxx_x_0[i] * pa_y[i] - ta_xxxxx_x_1[i] * pc_y[i];

        ta_xxxxxy_y_0[i] = 4.0 * ta_xxxy_y_0[i] * fe_0 - 4.0 * ta_xxxy_y_1[i] * fe_0 + ta_xxxxy_y_0[i] * pa_x[i] - ta_xxxxy_y_1[i] * pc_x[i];

        ta_xxxxxy_z_0[i] = ta_xxxxx_z_0[i] * pa_y[i] - ta_xxxxx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : IP

    auto ta_xxxxxz_x_0 = pbuffer.data(idx_npot_0_ip + 6);

    auto ta_xxxxxz_y_0 = pbuffer.data(idx_npot_0_ip + 7);

    auto ta_xxxxxz_z_0 = pbuffer.data(idx_npot_0_ip + 8);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxxxx_x_0,  \
                             ta_xxxxx_x_1,  \
                             ta_xxxxx_y_0,  \
                             ta_xxxxx_y_1,  \
                             ta_xxxxxz_x_0, \
                             ta_xxxxxz_y_0, \
                             ta_xxxxxz_z_0, \
                             ta_xxxxz_z_0,  \
                             ta_xxxxz_z_1,  \
                             ta_xxxz_z_0,   \
                             ta_xxxz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxz_x_0[i] = ta_xxxxx_x_0[i] * pa_z[i] - ta_xxxxx_x_1[i] * pc_z[i];

        ta_xxxxxz_y_0[i] = ta_xxxxx_y_0[i] * pa_z[i] - ta_xxxxx_y_1[i] * pc_z[i];

        ta_xxxxxz_z_0[i] = 4.0 * ta_xxxz_z_0[i] * fe_0 - 4.0 * ta_xxxz_z_1[i] * fe_0 + ta_xxxxz_z_0[i] * pa_x[i] - ta_xxxxz_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : IP

    auto ta_xxxxyy_x_0 = pbuffer.data(idx_npot_0_ip + 9);

    auto ta_xxxxyy_y_0 = pbuffer.data(idx_npot_0_ip + 10);

    auto ta_xxxxyy_z_0 = pbuffer.data(idx_npot_0_ip + 11);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxxx_x_0,   \
                             ta_xxxx_x_1,   \
                             ta_xxxxy_x_0,  \
                             ta_xxxxy_x_1,  \
                             ta_xxxxyy_x_0, \
                             ta_xxxxyy_y_0, \
                             ta_xxxxyy_z_0, \
                             ta_xxxyy_y_0,  \
                             ta_xxxyy_y_1,  \
                             ta_xxxyy_z_0,  \
                             ta_xxxyy_z_1,  \
                             ta_xxyy_y_0,   \
                             ta_xxyy_y_1,   \
                             ta_xxyy_z_0,   \
                             ta_xxyy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyy_x_0[i] = ta_xxxx_x_0[i] * fe_0 - ta_xxxx_x_1[i] * fe_0 + ta_xxxxy_x_0[i] * pa_y[i] - ta_xxxxy_x_1[i] * pc_y[i];

        ta_xxxxyy_y_0[i] = 3.0 * ta_xxyy_y_0[i] * fe_0 - 3.0 * ta_xxyy_y_1[i] * fe_0 + ta_xxxyy_y_0[i] * pa_x[i] - ta_xxxyy_y_1[i] * pc_x[i];

        ta_xxxxyy_z_0[i] = 3.0 * ta_xxyy_z_0[i] * fe_0 - 3.0 * ta_xxyy_z_1[i] * fe_0 + ta_xxxyy_z_0[i] * pa_x[i] - ta_xxxyy_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : IP

    auto ta_xxxxyz_x_0 = pbuffer.data(idx_npot_0_ip + 12);

    auto ta_xxxxyz_y_0 = pbuffer.data(idx_npot_0_ip + 13);

    auto ta_xxxxyz_z_0 = pbuffer.data(idx_npot_0_ip + 14);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_xxxxy_y_0,  \
                             ta_xxxxy_y_1,  \
                             ta_xxxxyz_x_0, \
                             ta_xxxxyz_y_0, \
                             ta_xxxxyz_z_0, \
                             ta_xxxxz_x_0,  \
                             ta_xxxxz_x_1,  \
                             ta_xxxxz_z_0,  \
                             ta_xxxxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xxxxyz_x_0[i] = ta_xxxxz_x_0[i] * pa_y[i] - ta_xxxxz_x_1[i] * pc_y[i];

        ta_xxxxyz_y_0[i] = ta_xxxxy_y_0[i] * pa_z[i] - ta_xxxxy_y_1[i] * pc_z[i];

        ta_xxxxyz_z_0[i] = ta_xxxxz_z_0[i] * pa_y[i] - ta_xxxxz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : IP

    auto ta_xxxxzz_x_0 = pbuffer.data(idx_npot_0_ip + 15);

    auto ta_xxxxzz_y_0 = pbuffer.data(idx_npot_0_ip + 16);

    auto ta_xxxxzz_z_0 = pbuffer.data(idx_npot_0_ip + 17);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxxx_x_0,   \
                             ta_xxxx_x_1,   \
                             ta_xxxxz_x_0,  \
                             ta_xxxxz_x_1,  \
                             ta_xxxxzz_x_0, \
                             ta_xxxxzz_y_0, \
                             ta_xxxxzz_z_0, \
                             ta_xxxzz_y_0,  \
                             ta_xxxzz_y_1,  \
                             ta_xxxzz_z_0,  \
                             ta_xxxzz_z_1,  \
                             ta_xxzz_y_0,   \
                             ta_xxzz_y_1,   \
                             ta_xxzz_z_0,   \
                             ta_xxzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxzz_x_0[i] = ta_xxxx_x_0[i] * fe_0 - ta_xxxx_x_1[i] * fe_0 + ta_xxxxz_x_0[i] * pa_z[i] - ta_xxxxz_x_1[i] * pc_z[i];

        ta_xxxxzz_y_0[i] = 3.0 * ta_xxzz_y_0[i] * fe_0 - 3.0 * ta_xxzz_y_1[i] * fe_0 + ta_xxxzz_y_0[i] * pa_x[i] - ta_xxxzz_y_1[i] * pc_x[i];

        ta_xxxxzz_z_0[i] = 3.0 * ta_xxzz_z_0[i] * fe_0 - 3.0 * ta_xxzz_z_1[i] * fe_0 + ta_xxxzz_z_0[i] * pa_x[i] - ta_xxxzz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : IP

    auto ta_xxxyyy_x_0 = pbuffer.data(idx_npot_0_ip + 18);

    auto ta_xxxyyy_y_0 = pbuffer.data(idx_npot_0_ip + 19);

    auto ta_xxxyyy_z_0 = pbuffer.data(idx_npot_0_ip + 20);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxxy_x_0,   \
                             ta_xxxy_x_1,   \
                             ta_xxxyy_x_0,  \
                             ta_xxxyy_x_1,  \
                             ta_xxxyyy_x_0, \
                             ta_xxxyyy_y_0, \
                             ta_xxxyyy_z_0, \
                             ta_xxyyy_y_0,  \
                             ta_xxyyy_y_1,  \
                             ta_xxyyy_z_0,  \
                             ta_xxyyy_z_1,  \
                             ta_xyyy_y_0,   \
                             ta_xyyy_y_1,   \
                             ta_xyyy_z_0,   \
                             ta_xyyy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyy_x_0[i] = 2.0 * ta_xxxy_x_0[i] * fe_0 - 2.0 * ta_xxxy_x_1[i] * fe_0 + ta_xxxyy_x_0[i] * pa_y[i] - ta_xxxyy_x_1[i] * pc_y[i];

        ta_xxxyyy_y_0[i] = 2.0 * ta_xyyy_y_0[i] * fe_0 - 2.0 * ta_xyyy_y_1[i] * fe_0 + ta_xxyyy_y_0[i] * pa_x[i] - ta_xxyyy_y_1[i] * pc_x[i];

        ta_xxxyyy_z_0[i] = 2.0 * ta_xyyy_z_0[i] * fe_0 - 2.0 * ta_xyyy_z_1[i] * fe_0 + ta_xxyyy_z_0[i] * pa_x[i] - ta_xxyyy_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : IP

    auto ta_xxxyyz_x_0 = pbuffer.data(idx_npot_0_ip + 21);

    auto ta_xxxyyz_y_0 = pbuffer.data(idx_npot_0_ip + 22);

    auto ta_xxxyyz_z_0 = pbuffer.data(idx_npot_0_ip + 23);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxxyy_x_0,  \
                             ta_xxxyy_x_1,  \
                             ta_xxxyy_y_0,  \
                             ta_xxxyy_y_1,  \
                             ta_xxxyyz_x_0, \
                             ta_xxxyyz_y_0, \
                             ta_xxxyyz_z_0, \
                             ta_xxyyz_z_0,  \
                             ta_xxyyz_z_1,  \
                             ta_xyyz_z_0,   \
                             ta_xyyz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyz_x_0[i] = ta_xxxyy_x_0[i] * pa_z[i] - ta_xxxyy_x_1[i] * pc_z[i];

        ta_xxxyyz_y_0[i] = ta_xxxyy_y_0[i] * pa_z[i] - ta_xxxyy_y_1[i] * pc_z[i];

        ta_xxxyyz_z_0[i] = 2.0 * ta_xyyz_z_0[i] * fe_0 - 2.0 * ta_xyyz_z_1[i] * fe_0 + ta_xxyyz_z_0[i] * pa_x[i] - ta_xxyyz_z_1[i] * pc_x[i];
    }

    // Set up 24-27 components of targeted buffer : IP

    auto ta_xxxyzz_x_0 = pbuffer.data(idx_npot_0_ip + 24);

    auto ta_xxxyzz_y_0 = pbuffer.data(idx_npot_0_ip + 25);

    auto ta_xxxyzz_z_0 = pbuffer.data(idx_npot_0_ip + 26);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxxyzz_x_0, \
                             ta_xxxyzz_y_0, \
                             ta_xxxyzz_z_0, \
                             ta_xxxzz_x_0,  \
                             ta_xxxzz_x_1,  \
                             ta_xxxzz_z_0,  \
                             ta_xxxzz_z_1,  \
                             ta_xxyzz_y_0,  \
                             ta_xxyzz_y_1,  \
                             ta_xyzz_y_0,   \
                             ta_xyzz_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyzz_x_0[i] = ta_xxxzz_x_0[i] * pa_y[i] - ta_xxxzz_x_1[i] * pc_y[i];

        ta_xxxyzz_y_0[i] = 2.0 * ta_xyzz_y_0[i] * fe_0 - 2.0 * ta_xyzz_y_1[i] * fe_0 + ta_xxyzz_y_0[i] * pa_x[i] - ta_xxyzz_y_1[i] * pc_x[i];

        ta_xxxyzz_z_0[i] = ta_xxxzz_z_0[i] * pa_y[i] - ta_xxxzz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : IP

    auto ta_xxxzzz_x_0 = pbuffer.data(idx_npot_0_ip + 27);

    auto ta_xxxzzz_y_0 = pbuffer.data(idx_npot_0_ip + 28);

    auto ta_xxxzzz_z_0 = pbuffer.data(idx_npot_0_ip + 29);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxxz_x_0,   \
                             ta_xxxz_x_1,   \
                             ta_xxxzz_x_0,  \
                             ta_xxxzz_x_1,  \
                             ta_xxxzzz_x_0, \
                             ta_xxxzzz_y_0, \
                             ta_xxxzzz_z_0, \
                             ta_xxzzz_y_0,  \
                             ta_xxzzz_y_1,  \
                             ta_xxzzz_z_0,  \
                             ta_xxzzz_z_1,  \
                             ta_xzzz_y_0,   \
                             ta_xzzz_y_1,   \
                             ta_xzzz_z_0,   \
                             ta_xzzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzzz_x_0[i] = 2.0 * ta_xxxz_x_0[i] * fe_0 - 2.0 * ta_xxxz_x_1[i] * fe_0 + ta_xxxzz_x_0[i] * pa_z[i] - ta_xxxzz_x_1[i] * pc_z[i];

        ta_xxxzzz_y_0[i] = 2.0 * ta_xzzz_y_0[i] * fe_0 - 2.0 * ta_xzzz_y_1[i] * fe_0 + ta_xxzzz_y_0[i] * pa_x[i] - ta_xxzzz_y_1[i] * pc_x[i];

        ta_xxxzzz_z_0[i] = 2.0 * ta_xzzz_z_0[i] * fe_0 - 2.0 * ta_xzzz_z_1[i] * fe_0 + ta_xxzzz_z_0[i] * pa_x[i] - ta_xxzzz_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : IP

    auto ta_xxyyyy_x_0 = pbuffer.data(idx_npot_0_ip + 30);

    auto ta_xxyyyy_y_0 = pbuffer.data(idx_npot_0_ip + 31);

    auto ta_xxyyyy_z_0 = pbuffer.data(idx_npot_0_ip + 32);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxyy_x_0,   \
                             ta_xxyy_x_1,   \
                             ta_xxyyy_x_0,  \
                             ta_xxyyy_x_1,  \
                             ta_xxyyyy_x_0, \
                             ta_xxyyyy_y_0, \
                             ta_xxyyyy_z_0, \
                             ta_xyyyy_y_0,  \
                             ta_xyyyy_y_1,  \
                             ta_xyyyy_z_0,  \
                             ta_xyyyy_z_1,  \
                             ta_yyyy_y_0,   \
                             ta_yyyy_y_1,   \
                             ta_yyyy_z_0,   \
                             ta_yyyy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyy_x_0[i] = 3.0 * ta_xxyy_x_0[i] * fe_0 - 3.0 * ta_xxyy_x_1[i] * fe_0 + ta_xxyyy_x_0[i] * pa_y[i] - ta_xxyyy_x_1[i] * pc_y[i];

        ta_xxyyyy_y_0[i] = ta_yyyy_y_0[i] * fe_0 - ta_yyyy_y_1[i] * fe_0 + ta_xyyyy_y_0[i] * pa_x[i] - ta_xyyyy_y_1[i] * pc_x[i];

        ta_xxyyyy_z_0[i] = ta_yyyy_z_0[i] * fe_0 - ta_yyyy_z_1[i] * fe_0 + ta_xyyyy_z_0[i] * pa_x[i] - ta_xyyyy_z_1[i] * pc_x[i];
    }

    // Set up 33-36 components of targeted buffer : IP

    auto ta_xxyyyz_x_0 = pbuffer.data(idx_npot_0_ip + 33);

    auto ta_xxyyyz_y_0 = pbuffer.data(idx_npot_0_ip + 34);

    auto ta_xxyyyz_z_0 = pbuffer.data(idx_npot_0_ip + 35);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxyyy_x_0,  \
                             ta_xxyyy_x_1,  \
                             ta_xxyyy_y_0,  \
                             ta_xxyyy_y_1,  \
                             ta_xxyyyz_x_0, \
                             ta_xxyyyz_y_0, \
                             ta_xxyyyz_z_0, \
                             ta_xyyyz_z_0,  \
                             ta_xyyyz_z_1,  \
                             ta_yyyz_z_0,   \
                             ta_yyyz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyz_x_0[i] = ta_xxyyy_x_0[i] * pa_z[i] - ta_xxyyy_x_1[i] * pc_z[i];

        ta_xxyyyz_y_0[i] = ta_xxyyy_y_0[i] * pa_z[i] - ta_xxyyy_y_1[i] * pc_z[i];

        ta_xxyyyz_z_0[i] = ta_yyyz_z_0[i] * fe_0 - ta_yyyz_z_1[i] * fe_0 + ta_xyyyz_z_0[i] * pa_x[i] - ta_xyyyz_z_1[i] * pc_x[i];
    }

    // Set up 36-39 components of targeted buffer : IP

    auto ta_xxyyzz_x_0 = pbuffer.data(idx_npot_0_ip + 36);

    auto ta_xxyyzz_y_0 = pbuffer.data(idx_npot_0_ip + 37);

    auto ta_xxyyzz_z_0 = pbuffer.data(idx_npot_0_ip + 38);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxyyzz_x_0, \
                             ta_xxyyzz_y_0, \
                             ta_xxyyzz_z_0, \
                             ta_xxyzz_x_0,  \
                             ta_xxyzz_x_1,  \
                             ta_xxzz_x_0,   \
                             ta_xxzz_x_1,   \
                             ta_xyyzz_y_0,  \
                             ta_xyyzz_y_1,  \
                             ta_xyyzz_z_0,  \
                             ta_xyyzz_z_1,  \
                             ta_yyzz_y_0,   \
                             ta_yyzz_y_1,   \
                             ta_yyzz_z_0,   \
                             ta_yyzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyzz_x_0[i] = ta_xxzz_x_0[i] * fe_0 - ta_xxzz_x_1[i] * fe_0 + ta_xxyzz_x_0[i] * pa_y[i] - ta_xxyzz_x_1[i] * pc_y[i];

        ta_xxyyzz_y_0[i] = ta_yyzz_y_0[i] * fe_0 - ta_yyzz_y_1[i] * fe_0 + ta_xyyzz_y_0[i] * pa_x[i] - ta_xyyzz_y_1[i] * pc_x[i];

        ta_xxyyzz_z_0[i] = ta_yyzz_z_0[i] * fe_0 - ta_yyzz_z_1[i] * fe_0 + ta_xyyzz_z_0[i] * pa_x[i] - ta_xyyzz_z_1[i] * pc_x[i];
    }

    // Set up 39-42 components of targeted buffer : IP

    auto ta_xxyzzz_x_0 = pbuffer.data(idx_npot_0_ip + 39);

    auto ta_xxyzzz_y_0 = pbuffer.data(idx_npot_0_ip + 40);

    auto ta_xxyzzz_z_0 = pbuffer.data(idx_npot_0_ip + 41);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxyzzz_x_0, \
                             ta_xxyzzz_y_0, \
                             ta_xxyzzz_z_0, \
                             ta_xxzzz_x_0,  \
                             ta_xxzzz_x_1,  \
                             ta_xxzzz_z_0,  \
                             ta_xxzzz_z_1,  \
                             ta_xyzzz_y_0,  \
                             ta_xyzzz_y_1,  \
                             ta_yzzz_y_0,   \
                             ta_yzzz_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzzz_x_0[i] = ta_xxzzz_x_0[i] * pa_y[i] - ta_xxzzz_x_1[i] * pc_y[i];

        ta_xxyzzz_y_0[i] = ta_yzzz_y_0[i] * fe_0 - ta_yzzz_y_1[i] * fe_0 + ta_xyzzz_y_0[i] * pa_x[i] - ta_xyzzz_y_1[i] * pc_x[i];

        ta_xxyzzz_z_0[i] = ta_xxzzz_z_0[i] * pa_y[i] - ta_xxzzz_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : IP

    auto ta_xxzzzz_x_0 = pbuffer.data(idx_npot_0_ip + 42);

    auto ta_xxzzzz_y_0 = pbuffer.data(idx_npot_0_ip + 43);

    auto ta_xxzzzz_z_0 = pbuffer.data(idx_npot_0_ip + 44);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxzz_x_0,   \
                             ta_xxzz_x_1,   \
                             ta_xxzzz_x_0,  \
                             ta_xxzzz_x_1,  \
                             ta_xxzzzz_x_0, \
                             ta_xxzzzz_y_0, \
                             ta_xxzzzz_z_0, \
                             ta_xzzzz_y_0,  \
                             ta_xzzzz_y_1,  \
                             ta_xzzzz_z_0,  \
                             ta_xzzzz_z_1,  \
                             ta_zzzz_y_0,   \
                             ta_zzzz_y_1,   \
                             ta_zzzz_z_0,   \
                             ta_zzzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzzz_x_0[i] = 3.0 * ta_xxzz_x_0[i] * fe_0 - 3.0 * ta_xxzz_x_1[i] * fe_0 + ta_xxzzz_x_0[i] * pa_z[i] - ta_xxzzz_x_1[i] * pc_z[i];

        ta_xxzzzz_y_0[i] = ta_zzzz_y_0[i] * fe_0 - ta_zzzz_y_1[i] * fe_0 + ta_xzzzz_y_0[i] * pa_x[i] - ta_xzzzz_y_1[i] * pc_x[i];

        ta_xxzzzz_z_0[i] = ta_zzzz_z_0[i] * fe_0 - ta_zzzz_z_1[i] * fe_0 + ta_xzzzz_z_0[i] * pa_x[i] - ta_xzzzz_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : IP

    auto ta_xyyyyy_x_0 = pbuffer.data(idx_npot_0_ip + 45);

    auto ta_xyyyyy_y_0 = pbuffer.data(idx_npot_0_ip + 46);

    auto ta_xyyyyy_z_0 = pbuffer.data(idx_npot_0_ip + 47);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xyyyyy_x_0, \
                             ta_xyyyyy_y_0, \
                             ta_xyyyyy_z_0, \
                             ta_yyyyy_0_0,  \
                             ta_yyyyy_0_1,  \
                             ta_yyyyy_x_0,  \
                             ta_yyyyy_x_1,  \
                             ta_yyyyy_y_0,  \
                             ta_yyyyy_y_1,  \
                             ta_yyyyy_z_0,  \
                             ta_yyyyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyy_x_0[i] = ta_yyyyy_0_0[i] * fe_0 - ta_yyyyy_0_1[i] * fe_0 + ta_yyyyy_x_0[i] * pa_x[i] - ta_yyyyy_x_1[i] * pc_x[i];

        ta_xyyyyy_y_0[i] = ta_yyyyy_y_0[i] * pa_x[i] - ta_yyyyy_y_1[i] * pc_x[i];

        ta_xyyyyy_z_0[i] = ta_yyyyy_z_0[i] * pa_x[i] - ta_yyyyy_z_1[i] * pc_x[i];
    }

    // Set up 48-51 components of targeted buffer : IP

    auto ta_xyyyyz_x_0 = pbuffer.data(idx_npot_0_ip + 48);

    auto ta_xyyyyz_y_0 = pbuffer.data(idx_npot_0_ip + 49);

    auto ta_xyyyyz_z_0 = pbuffer.data(idx_npot_0_ip + 50);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xyyyy_x_0,  \
                             ta_xyyyy_x_1,  \
                             ta_xyyyyz_x_0, \
                             ta_xyyyyz_y_0, \
                             ta_xyyyyz_z_0, \
                             ta_yyyyz_y_0,  \
                             ta_yyyyz_y_1,  \
                             ta_yyyyz_z_0,  \
                             ta_yyyyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyyyyz_x_0[i] = ta_xyyyy_x_0[i] * pa_z[i] - ta_xyyyy_x_1[i] * pc_z[i];

        ta_xyyyyz_y_0[i] = ta_yyyyz_y_0[i] * pa_x[i] - ta_yyyyz_y_1[i] * pc_x[i];

        ta_xyyyyz_z_0[i] = ta_yyyyz_z_0[i] * pa_x[i] - ta_yyyyz_z_1[i] * pc_x[i];
    }

    // Set up 51-54 components of targeted buffer : IP

    auto ta_xyyyzz_x_0 = pbuffer.data(idx_npot_0_ip + 51);

    auto ta_xyyyzz_y_0 = pbuffer.data(idx_npot_0_ip + 52);

    auto ta_xyyyzz_z_0 = pbuffer.data(idx_npot_0_ip + 53);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xyyyzz_x_0, \
                             ta_xyyyzz_y_0, \
                             ta_xyyyzz_z_0, \
                             ta_yyyzz_0_0,  \
                             ta_yyyzz_0_1,  \
                             ta_yyyzz_x_0,  \
                             ta_yyyzz_x_1,  \
                             ta_yyyzz_y_0,  \
                             ta_yyyzz_y_1,  \
                             ta_yyyzz_z_0,  \
                             ta_yyyzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyzz_x_0[i] = ta_yyyzz_0_0[i] * fe_0 - ta_yyyzz_0_1[i] * fe_0 + ta_yyyzz_x_0[i] * pa_x[i] - ta_yyyzz_x_1[i] * pc_x[i];

        ta_xyyyzz_y_0[i] = ta_yyyzz_y_0[i] * pa_x[i] - ta_yyyzz_y_1[i] * pc_x[i];

        ta_xyyyzz_z_0[i] = ta_yyyzz_z_0[i] * pa_x[i] - ta_yyyzz_z_1[i] * pc_x[i];
    }

    // Set up 54-57 components of targeted buffer : IP

    auto ta_xyyzzz_x_0 = pbuffer.data(idx_npot_0_ip + 54);

    auto ta_xyyzzz_y_0 = pbuffer.data(idx_npot_0_ip + 55);

    auto ta_xyyzzz_z_0 = pbuffer.data(idx_npot_0_ip + 56);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xyyzzz_x_0, \
                             ta_xyyzzz_y_0, \
                             ta_xyyzzz_z_0, \
                             ta_yyzzz_0_0,  \
                             ta_yyzzz_0_1,  \
                             ta_yyzzz_x_0,  \
                             ta_yyzzz_x_1,  \
                             ta_yyzzz_y_0,  \
                             ta_yyzzz_y_1,  \
                             ta_yyzzz_z_0,  \
                             ta_yyzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzzz_x_0[i] = ta_yyzzz_0_0[i] * fe_0 - ta_yyzzz_0_1[i] * fe_0 + ta_yyzzz_x_0[i] * pa_x[i] - ta_yyzzz_x_1[i] * pc_x[i];

        ta_xyyzzz_y_0[i] = ta_yyzzz_y_0[i] * pa_x[i] - ta_yyzzz_y_1[i] * pc_x[i];

        ta_xyyzzz_z_0[i] = ta_yyzzz_z_0[i] * pa_x[i] - ta_yyzzz_z_1[i] * pc_x[i];
    }

    // Set up 57-60 components of targeted buffer : IP

    auto ta_xyzzzz_x_0 = pbuffer.data(idx_npot_0_ip + 57);

    auto ta_xyzzzz_y_0 = pbuffer.data(idx_npot_0_ip + 58);

    auto ta_xyzzzz_z_0 = pbuffer.data(idx_npot_0_ip + 59);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xyzzzz_x_0, \
                             ta_xyzzzz_y_0, \
                             ta_xyzzzz_z_0, \
                             ta_xzzzz_x_0,  \
                             ta_xzzzz_x_1,  \
                             ta_yzzzz_y_0,  \
                             ta_yzzzz_y_1,  \
                             ta_yzzzz_z_0,  \
                             ta_yzzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyzzzz_x_0[i] = ta_xzzzz_x_0[i] * pa_y[i] - ta_xzzzz_x_1[i] * pc_y[i];

        ta_xyzzzz_y_0[i] = ta_yzzzz_y_0[i] * pa_x[i] - ta_yzzzz_y_1[i] * pc_x[i];

        ta_xyzzzz_z_0[i] = ta_yzzzz_z_0[i] * pa_x[i] - ta_yzzzz_z_1[i] * pc_x[i];
    }

    // Set up 60-63 components of targeted buffer : IP

    auto ta_xzzzzz_x_0 = pbuffer.data(idx_npot_0_ip + 60);

    auto ta_xzzzzz_y_0 = pbuffer.data(idx_npot_0_ip + 61);

    auto ta_xzzzzz_z_0 = pbuffer.data(idx_npot_0_ip + 62);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xzzzzz_x_0, \
                             ta_xzzzzz_y_0, \
                             ta_xzzzzz_z_0, \
                             ta_zzzzz_0_0,  \
                             ta_zzzzz_0_1,  \
                             ta_zzzzz_x_0,  \
                             ta_zzzzz_x_1,  \
                             ta_zzzzz_y_0,  \
                             ta_zzzzz_y_1,  \
                             ta_zzzzz_z_0,  \
                             ta_zzzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzzz_x_0[i] = ta_zzzzz_0_0[i] * fe_0 - ta_zzzzz_0_1[i] * fe_0 + ta_zzzzz_x_0[i] * pa_x[i] - ta_zzzzz_x_1[i] * pc_x[i];

        ta_xzzzzz_y_0[i] = ta_zzzzz_y_0[i] * pa_x[i] - ta_zzzzz_y_1[i] * pc_x[i];

        ta_xzzzzz_z_0[i] = ta_zzzzz_z_0[i] * pa_x[i] - ta_zzzzz_z_1[i] * pc_x[i];
    }

    // Set up 63-66 components of targeted buffer : IP

    auto ta_yyyyyy_x_0 = pbuffer.data(idx_npot_0_ip + 63);

    auto ta_yyyyyy_y_0 = pbuffer.data(idx_npot_0_ip + 64);

    auto ta_yyyyyy_z_0 = pbuffer.data(idx_npot_0_ip + 65);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_yyyy_x_0,   \
                             ta_yyyy_x_1,   \
                             ta_yyyy_y_0,   \
                             ta_yyyy_y_1,   \
                             ta_yyyy_z_0,   \
                             ta_yyyy_z_1,   \
                             ta_yyyyy_0_0,  \
                             ta_yyyyy_0_1,  \
                             ta_yyyyy_x_0,  \
                             ta_yyyyy_x_1,  \
                             ta_yyyyy_y_0,  \
                             ta_yyyyy_y_1,  \
                             ta_yyyyy_z_0,  \
                             ta_yyyyy_z_1,  \
                             ta_yyyyyy_x_0, \
                             ta_yyyyyy_y_0, \
                             ta_yyyyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyy_x_0[i] = 5.0 * ta_yyyy_x_0[i] * fe_0 - 5.0 * ta_yyyy_x_1[i] * fe_0 + ta_yyyyy_x_0[i] * pa_y[i] - ta_yyyyy_x_1[i] * pc_y[i];

        ta_yyyyyy_y_0[i] = 5.0 * ta_yyyy_y_0[i] * fe_0 - 5.0 * ta_yyyy_y_1[i] * fe_0 + ta_yyyyy_0_0[i] * fe_0 - ta_yyyyy_0_1[i] * fe_0 +
                           ta_yyyyy_y_0[i] * pa_y[i] - ta_yyyyy_y_1[i] * pc_y[i];

        ta_yyyyyy_z_0[i] = 5.0 * ta_yyyy_z_0[i] * fe_0 - 5.0 * ta_yyyy_z_1[i] * fe_0 + ta_yyyyy_z_0[i] * pa_y[i] - ta_yyyyy_z_1[i] * pc_y[i];
    }

    // Set up 66-69 components of targeted buffer : IP

    auto ta_yyyyyz_x_0 = pbuffer.data(idx_npot_0_ip + 66);

    auto ta_yyyyyz_y_0 = pbuffer.data(idx_npot_0_ip + 67);

    auto ta_yyyyyz_z_0 = pbuffer.data(idx_npot_0_ip + 68);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyyyy_x_0,  \
                             ta_yyyyy_x_1,  \
                             ta_yyyyy_y_0,  \
                             ta_yyyyy_y_1,  \
                             ta_yyyyyz_x_0, \
                             ta_yyyyyz_y_0, \
                             ta_yyyyyz_z_0, \
                             ta_yyyyz_z_0,  \
                             ta_yyyyz_z_1,  \
                             ta_yyyz_z_0,   \
                             ta_yyyz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyz_x_0[i] = ta_yyyyy_x_0[i] * pa_z[i] - ta_yyyyy_x_1[i] * pc_z[i];

        ta_yyyyyz_y_0[i] = ta_yyyyy_y_0[i] * pa_z[i] - ta_yyyyy_y_1[i] * pc_z[i];

        ta_yyyyyz_z_0[i] = 4.0 * ta_yyyz_z_0[i] * fe_0 - 4.0 * ta_yyyz_z_1[i] * fe_0 + ta_yyyyz_z_0[i] * pa_y[i] - ta_yyyyz_z_1[i] * pc_y[i];
    }

    // Set up 69-72 components of targeted buffer : IP

    auto ta_yyyyzz_x_0 = pbuffer.data(idx_npot_0_ip + 69);

    auto ta_yyyyzz_y_0 = pbuffer.data(idx_npot_0_ip + 70);

    auto ta_yyyyzz_z_0 = pbuffer.data(idx_npot_0_ip + 71);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyyy_y_0,   \
                             ta_yyyy_y_1,   \
                             ta_yyyyz_y_0,  \
                             ta_yyyyz_y_1,  \
                             ta_yyyyzz_x_0, \
                             ta_yyyyzz_y_0, \
                             ta_yyyyzz_z_0, \
                             ta_yyyzz_x_0,  \
                             ta_yyyzz_x_1,  \
                             ta_yyyzz_z_0,  \
                             ta_yyyzz_z_1,  \
                             ta_yyzz_x_0,   \
                             ta_yyzz_x_1,   \
                             ta_yyzz_z_0,   \
                             ta_yyzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyzz_x_0[i] = 3.0 * ta_yyzz_x_0[i] * fe_0 - 3.0 * ta_yyzz_x_1[i] * fe_0 + ta_yyyzz_x_0[i] * pa_y[i] - ta_yyyzz_x_1[i] * pc_y[i];

        ta_yyyyzz_y_0[i] = ta_yyyy_y_0[i] * fe_0 - ta_yyyy_y_1[i] * fe_0 + ta_yyyyz_y_0[i] * pa_z[i] - ta_yyyyz_y_1[i] * pc_z[i];

        ta_yyyyzz_z_0[i] = 3.0 * ta_yyzz_z_0[i] * fe_0 - 3.0 * ta_yyzz_z_1[i] * fe_0 + ta_yyyzz_z_0[i] * pa_y[i] - ta_yyyzz_z_1[i] * pc_y[i];
    }

    // Set up 72-75 components of targeted buffer : IP

    auto ta_yyyzzz_x_0 = pbuffer.data(idx_npot_0_ip + 72);

    auto ta_yyyzzz_y_0 = pbuffer.data(idx_npot_0_ip + 73);

    auto ta_yyyzzz_z_0 = pbuffer.data(idx_npot_0_ip + 74);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyyz_y_0,   \
                             ta_yyyz_y_1,   \
                             ta_yyyzz_y_0,  \
                             ta_yyyzz_y_1,  \
                             ta_yyyzzz_x_0, \
                             ta_yyyzzz_y_0, \
                             ta_yyyzzz_z_0, \
                             ta_yyzzz_x_0,  \
                             ta_yyzzz_x_1,  \
                             ta_yyzzz_z_0,  \
                             ta_yyzzz_z_1,  \
                             ta_yzzz_x_0,   \
                             ta_yzzz_x_1,   \
                             ta_yzzz_z_0,   \
                             ta_yzzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzzz_x_0[i] = 2.0 * ta_yzzz_x_0[i] * fe_0 - 2.0 * ta_yzzz_x_1[i] * fe_0 + ta_yyzzz_x_0[i] * pa_y[i] - ta_yyzzz_x_1[i] * pc_y[i];

        ta_yyyzzz_y_0[i] = 2.0 * ta_yyyz_y_0[i] * fe_0 - 2.0 * ta_yyyz_y_1[i] * fe_0 + ta_yyyzz_y_0[i] * pa_z[i] - ta_yyyzz_y_1[i] * pc_z[i];

        ta_yyyzzz_z_0[i] = 2.0 * ta_yzzz_z_0[i] * fe_0 - 2.0 * ta_yzzz_z_1[i] * fe_0 + ta_yyzzz_z_0[i] * pa_y[i] - ta_yyzzz_z_1[i] * pc_y[i];
    }

    // Set up 75-78 components of targeted buffer : IP

    auto ta_yyzzzz_x_0 = pbuffer.data(idx_npot_0_ip + 75);

    auto ta_yyzzzz_y_0 = pbuffer.data(idx_npot_0_ip + 76);

    auto ta_yyzzzz_z_0 = pbuffer.data(idx_npot_0_ip + 77);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyzz_y_0,   \
                             ta_yyzz_y_1,   \
                             ta_yyzzz_y_0,  \
                             ta_yyzzz_y_1,  \
                             ta_yyzzzz_x_0, \
                             ta_yyzzzz_y_0, \
                             ta_yyzzzz_z_0, \
                             ta_yzzzz_x_0,  \
                             ta_yzzzz_x_1,  \
                             ta_yzzzz_z_0,  \
                             ta_yzzzz_z_1,  \
                             ta_zzzz_x_0,   \
                             ta_zzzz_x_1,   \
                             ta_zzzz_z_0,   \
                             ta_zzzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzzz_x_0[i] = ta_zzzz_x_0[i] * fe_0 - ta_zzzz_x_1[i] * fe_0 + ta_yzzzz_x_0[i] * pa_y[i] - ta_yzzzz_x_1[i] * pc_y[i];

        ta_yyzzzz_y_0[i] = 3.0 * ta_yyzz_y_0[i] * fe_0 - 3.0 * ta_yyzz_y_1[i] * fe_0 + ta_yyzzz_y_0[i] * pa_z[i] - ta_yyzzz_y_1[i] * pc_z[i];

        ta_yyzzzz_z_0[i] = ta_zzzz_z_0[i] * fe_0 - ta_zzzz_z_1[i] * fe_0 + ta_yzzzz_z_0[i] * pa_y[i] - ta_yzzzz_z_1[i] * pc_y[i];
    }

    // Set up 78-81 components of targeted buffer : IP

    auto ta_yzzzzz_x_0 = pbuffer.data(idx_npot_0_ip + 78);

    auto ta_yzzzzz_y_0 = pbuffer.data(idx_npot_0_ip + 79);

    auto ta_yzzzzz_z_0 = pbuffer.data(idx_npot_0_ip + 80);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_yzzzzz_x_0, \
                             ta_yzzzzz_y_0, \
                             ta_yzzzzz_z_0, \
                             ta_zzzzz_0_0,  \
                             ta_zzzzz_0_1,  \
                             ta_zzzzz_x_0,  \
                             ta_zzzzz_x_1,  \
                             ta_zzzzz_y_0,  \
                             ta_zzzzz_y_1,  \
                             ta_zzzzz_z_0,  \
                             ta_zzzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzzz_x_0[i] = ta_zzzzz_x_0[i] * pa_y[i] - ta_zzzzz_x_1[i] * pc_y[i];

        ta_yzzzzz_y_0[i] = ta_zzzzz_0_0[i] * fe_0 - ta_zzzzz_0_1[i] * fe_0 + ta_zzzzz_y_0[i] * pa_y[i] - ta_zzzzz_y_1[i] * pc_y[i];

        ta_yzzzzz_z_0[i] = ta_zzzzz_z_0[i] * pa_y[i] - ta_zzzzz_z_1[i] * pc_y[i];
    }

    // Set up 81-84 components of targeted buffer : IP

    auto ta_zzzzzz_x_0 = pbuffer.data(idx_npot_0_ip + 81);

    auto ta_zzzzzz_y_0 = pbuffer.data(idx_npot_0_ip + 82);

    auto ta_zzzzzz_z_0 = pbuffer.data(idx_npot_0_ip + 83);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta_zzzz_x_0,   \
                             ta_zzzz_x_1,   \
                             ta_zzzz_y_0,   \
                             ta_zzzz_y_1,   \
                             ta_zzzz_z_0,   \
                             ta_zzzz_z_1,   \
                             ta_zzzzz_0_0,  \
                             ta_zzzzz_0_1,  \
                             ta_zzzzz_x_0,  \
                             ta_zzzzz_x_1,  \
                             ta_zzzzz_y_0,  \
                             ta_zzzzz_y_1,  \
                             ta_zzzzz_z_0,  \
                             ta_zzzzz_z_1,  \
                             ta_zzzzzz_x_0, \
                             ta_zzzzzz_y_0, \
                             ta_zzzzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzzz_x_0[i] = 5.0 * ta_zzzz_x_0[i] * fe_0 - 5.0 * ta_zzzz_x_1[i] * fe_0 + ta_zzzzz_x_0[i] * pa_z[i] - ta_zzzzz_x_1[i] * pc_z[i];

        ta_zzzzzz_y_0[i] = 5.0 * ta_zzzz_y_0[i] * fe_0 - 5.0 * ta_zzzz_y_1[i] * fe_0 + ta_zzzzz_y_0[i] * pa_z[i] - ta_zzzzz_y_1[i] * pc_z[i];

        ta_zzzzzz_z_0[i] = 5.0 * ta_zzzz_z_0[i] * fe_0 - 5.0 * ta_zzzz_z_1[i] * fe_0 + ta_zzzzz_0_0[i] * fe_0 - ta_zzzzz_0_1[i] * fe_0 +
                           ta_zzzzz_z_0[i] * pa_z[i] - ta_zzzzz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
