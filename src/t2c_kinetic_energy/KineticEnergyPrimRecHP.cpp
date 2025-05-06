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

#include "KineticEnergyPrimRecHP.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_hp(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_hp,
                            const size_t              idx_ovl_fp,
                            const size_t              idx_kin_fp,
                            const size_t              idx_kin_gs,
                            const size_t              idx_kin_gp,
                            const size_t              idx_ovl_hp,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_ovl_fp);

    auto ts_xxx_y = pbuffer.data(idx_ovl_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_ovl_fp + 2);

    auto ts_xxy_x = pbuffer.data(idx_ovl_fp + 3);

    auto ts_xxz_x = pbuffer.data(idx_ovl_fp + 6);

    auto ts_xyy_y = pbuffer.data(idx_ovl_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_ovl_fp + 11);

    auto ts_xzz_y = pbuffer.data(idx_ovl_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_ovl_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_ovl_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_ovl_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_ovl_fp + 20);

    auto ts_yyz_y = pbuffer.data(idx_ovl_fp + 22);

    auto ts_yzz_x = pbuffer.data(idx_ovl_fp + 24);

    auto ts_yzz_z = pbuffer.data(idx_ovl_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_ovl_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_ovl_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_ovl_fp + 29);

    // Set up components of auxiliary buffer : FP

    auto tk_xxx_x = pbuffer.data(idx_kin_fp);

    auto tk_xxx_y = pbuffer.data(idx_kin_fp + 1);

    auto tk_xxx_z = pbuffer.data(idx_kin_fp + 2);

    auto tk_xxy_x = pbuffer.data(idx_kin_fp + 3);

    auto tk_xxz_x = pbuffer.data(idx_kin_fp + 6);

    auto tk_xyy_y = pbuffer.data(idx_kin_fp + 10);

    auto tk_xyy_z = pbuffer.data(idx_kin_fp + 11);

    auto tk_xzz_y = pbuffer.data(idx_kin_fp + 16);

    auto tk_xzz_z = pbuffer.data(idx_kin_fp + 17);

    auto tk_yyy_x = pbuffer.data(idx_kin_fp + 18);

    auto tk_yyy_y = pbuffer.data(idx_kin_fp + 19);

    auto tk_yyy_z = pbuffer.data(idx_kin_fp + 20);

    auto tk_yyz_y = pbuffer.data(idx_kin_fp + 22);

    auto tk_yzz_x = pbuffer.data(idx_kin_fp + 24);

    auto tk_yzz_z = pbuffer.data(idx_kin_fp + 26);

    auto tk_zzz_x = pbuffer.data(idx_kin_fp + 27);

    auto tk_zzz_y = pbuffer.data(idx_kin_fp + 28);

    auto tk_zzz_z = pbuffer.data(idx_kin_fp + 29);

    // Set up components of auxiliary buffer : GS

    auto tk_xxxx_0 = pbuffer.data(idx_kin_gs);

    auto tk_xxyy_0 = pbuffer.data(idx_kin_gs + 3);

    auto tk_xxzz_0 = pbuffer.data(idx_kin_gs + 5);

    auto tk_yyyy_0 = pbuffer.data(idx_kin_gs + 10);

    auto tk_yyzz_0 = pbuffer.data(idx_kin_gs + 12);

    auto tk_zzzz_0 = pbuffer.data(idx_kin_gs + 14);

    // Set up components of auxiliary buffer : GP

    auto tk_xxxx_x = pbuffer.data(idx_kin_gp);

    auto tk_xxxx_y = pbuffer.data(idx_kin_gp + 1);

    auto tk_xxxx_z = pbuffer.data(idx_kin_gp + 2);

    auto tk_xxxy_x = pbuffer.data(idx_kin_gp + 3);

    auto tk_xxxy_y = pbuffer.data(idx_kin_gp + 4);

    auto tk_xxxz_x = pbuffer.data(idx_kin_gp + 6);

    auto tk_xxxz_z = pbuffer.data(idx_kin_gp + 8);

    auto tk_xxyy_x = pbuffer.data(idx_kin_gp + 9);

    auto tk_xxyy_y = pbuffer.data(idx_kin_gp + 10);

    auto tk_xxyy_z = pbuffer.data(idx_kin_gp + 11);

    auto tk_xxzz_x = pbuffer.data(idx_kin_gp + 15);

    auto tk_xxzz_y = pbuffer.data(idx_kin_gp + 16);

    auto tk_xxzz_z = pbuffer.data(idx_kin_gp + 17);

    auto tk_xyyy_x = pbuffer.data(idx_kin_gp + 18);

    auto tk_xyyy_y = pbuffer.data(idx_kin_gp + 19);

    auto tk_xyyy_z = pbuffer.data(idx_kin_gp + 20);

    auto tk_xzzz_x = pbuffer.data(idx_kin_gp + 27);

    auto tk_xzzz_y = pbuffer.data(idx_kin_gp + 28);

    auto tk_xzzz_z = pbuffer.data(idx_kin_gp + 29);

    auto tk_yyyy_x = pbuffer.data(idx_kin_gp + 30);

    auto tk_yyyy_y = pbuffer.data(idx_kin_gp + 31);

    auto tk_yyyy_z = pbuffer.data(idx_kin_gp + 32);

    auto tk_yyyz_y = pbuffer.data(idx_kin_gp + 34);

    auto tk_yyyz_z = pbuffer.data(idx_kin_gp + 35);

    auto tk_yyzz_x = pbuffer.data(idx_kin_gp + 36);

    auto tk_yyzz_y = pbuffer.data(idx_kin_gp + 37);

    auto tk_yyzz_z = pbuffer.data(idx_kin_gp + 38);

    auto tk_yzzz_x = pbuffer.data(idx_kin_gp + 39);

    auto tk_yzzz_y = pbuffer.data(idx_kin_gp + 40);

    auto tk_yzzz_z = pbuffer.data(idx_kin_gp + 41);

    auto tk_zzzz_x = pbuffer.data(idx_kin_gp + 42);

    auto tk_zzzz_y = pbuffer.data(idx_kin_gp + 43);

    auto tk_zzzz_z = pbuffer.data(idx_kin_gp + 44);

    // Set up components of auxiliary buffer : HP

    auto ts_xxxxx_x = pbuffer.data(idx_ovl_hp);

    auto ts_xxxxx_y = pbuffer.data(idx_ovl_hp + 1);

    auto ts_xxxxx_z = pbuffer.data(idx_ovl_hp + 2);

    auto ts_xxxxy_x = pbuffer.data(idx_ovl_hp + 3);

    auto ts_xxxxy_y = pbuffer.data(idx_ovl_hp + 4);

    auto ts_xxxxy_z = pbuffer.data(idx_ovl_hp + 5);

    auto ts_xxxxz_x = pbuffer.data(idx_ovl_hp + 6);

    auto ts_xxxxz_y = pbuffer.data(idx_ovl_hp + 7);

    auto ts_xxxxz_z = pbuffer.data(idx_ovl_hp + 8);

    auto ts_xxxyy_x = pbuffer.data(idx_ovl_hp + 9);

    auto ts_xxxyy_y = pbuffer.data(idx_ovl_hp + 10);

    auto ts_xxxyy_z = pbuffer.data(idx_ovl_hp + 11);

    auto ts_xxxyz_x = pbuffer.data(idx_ovl_hp + 12);

    auto ts_xxxyz_y = pbuffer.data(idx_ovl_hp + 13);

    auto ts_xxxyz_z = pbuffer.data(idx_ovl_hp + 14);

    auto ts_xxxzz_x = pbuffer.data(idx_ovl_hp + 15);

    auto ts_xxxzz_y = pbuffer.data(idx_ovl_hp + 16);

    auto ts_xxxzz_z = pbuffer.data(idx_ovl_hp + 17);

    auto ts_xxyyy_x = pbuffer.data(idx_ovl_hp + 18);

    auto ts_xxyyy_y = pbuffer.data(idx_ovl_hp + 19);

    auto ts_xxyyy_z = pbuffer.data(idx_ovl_hp + 20);

    auto ts_xxyyz_x = pbuffer.data(idx_ovl_hp + 21);

    auto ts_xxyyz_y = pbuffer.data(idx_ovl_hp + 22);

    auto ts_xxyyz_z = pbuffer.data(idx_ovl_hp + 23);

    auto ts_xxyzz_x = pbuffer.data(idx_ovl_hp + 24);

    auto ts_xxyzz_y = pbuffer.data(idx_ovl_hp + 25);

    auto ts_xxyzz_z = pbuffer.data(idx_ovl_hp + 26);

    auto ts_xxzzz_x = pbuffer.data(idx_ovl_hp + 27);

    auto ts_xxzzz_y = pbuffer.data(idx_ovl_hp + 28);

    auto ts_xxzzz_z = pbuffer.data(idx_ovl_hp + 29);

    auto ts_xyyyy_x = pbuffer.data(idx_ovl_hp + 30);

    auto ts_xyyyy_y = pbuffer.data(idx_ovl_hp + 31);

    auto ts_xyyyy_z = pbuffer.data(idx_ovl_hp + 32);

    auto ts_xyyyz_x = pbuffer.data(idx_ovl_hp + 33);

    auto ts_xyyyz_y = pbuffer.data(idx_ovl_hp + 34);

    auto ts_xyyyz_z = pbuffer.data(idx_ovl_hp + 35);

    auto ts_xyyzz_x = pbuffer.data(idx_ovl_hp + 36);

    auto ts_xyyzz_y = pbuffer.data(idx_ovl_hp + 37);

    auto ts_xyyzz_z = pbuffer.data(idx_ovl_hp + 38);

    auto ts_xyzzz_x = pbuffer.data(idx_ovl_hp + 39);

    auto ts_xyzzz_y = pbuffer.data(idx_ovl_hp + 40);

    auto ts_xyzzz_z = pbuffer.data(idx_ovl_hp + 41);

    auto ts_xzzzz_x = pbuffer.data(idx_ovl_hp + 42);

    auto ts_xzzzz_y = pbuffer.data(idx_ovl_hp + 43);

    auto ts_xzzzz_z = pbuffer.data(idx_ovl_hp + 44);

    auto ts_yyyyy_x = pbuffer.data(idx_ovl_hp + 45);

    auto ts_yyyyy_y = pbuffer.data(idx_ovl_hp + 46);

    auto ts_yyyyy_z = pbuffer.data(idx_ovl_hp + 47);

    auto ts_yyyyz_x = pbuffer.data(idx_ovl_hp + 48);

    auto ts_yyyyz_y = pbuffer.data(idx_ovl_hp + 49);

    auto ts_yyyyz_z = pbuffer.data(idx_ovl_hp + 50);

    auto ts_yyyzz_x = pbuffer.data(idx_ovl_hp + 51);

    auto ts_yyyzz_y = pbuffer.data(idx_ovl_hp + 52);

    auto ts_yyyzz_z = pbuffer.data(idx_ovl_hp + 53);

    auto ts_yyzzz_x = pbuffer.data(idx_ovl_hp + 54);

    auto ts_yyzzz_y = pbuffer.data(idx_ovl_hp + 55);

    auto ts_yyzzz_z = pbuffer.data(idx_ovl_hp + 56);

    auto ts_yzzzz_x = pbuffer.data(idx_ovl_hp + 57);

    auto ts_yzzzz_y = pbuffer.data(idx_ovl_hp + 58);

    auto ts_yzzzz_z = pbuffer.data(idx_ovl_hp + 59);

    auto ts_zzzzz_x = pbuffer.data(idx_ovl_hp + 60);

    auto ts_zzzzz_y = pbuffer.data(idx_ovl_hp + 61);

    auto ts_zzzzz_z = pbuffer.data(idx_ovl_hp + 62);

    // Set up 0-3 components of targeted buffer : HP

    auto tk_xxxxx_x = pbuffer.data(idx_kin_hp);

    auto tk_xxxxx_y = pbuffer.data(idx_kin_hp + 1);

    auto tk_xxxxx_z = pbuffer.data(idx_kin_hp + 2);

#pragma omp simd aligned(pa_x,           \
                             tk_xxx_x,   \
                             tk_xxx_y,   \
                             tk_xxx_z,   \
                             tk_xxxx_0,  \
                             tk_xxxx_x,  \
                             tk_xxxx_y,  \
                             tk_xxxx_z,  \
                             tk_xxxxx_x, \
                             tk_xxxxx_y, \
                             tk_xxxxx_z, \
                             ts_xxx_x,   \
                             ts_xxx_y,   \
                             ts_xxx_z,   \
                             ts_xxxxx_x, \
                             ts_xxxxx_y, \
                             ts_xxxxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxx_x[i] =
            -8.0 * ts_xxx_x[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_x[i] * fe_0 + tk_xxxx_0[i] * fe_0 + tk_xxxx_x[i] * pa_x[i] + 2.0 * ts_xxxxx_x[i] * fz_0;

        tk_xxxxx_y[i] = -8.0 * ts_xxx_y[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_y[i] * fe_0 + tk_xxxx_y[i] * pa_x[i] + 2.0 * ts_xxxxx_y[i] * fz_0;

        tk_xxxxx_z[i] = -8.0 * ts_xxx_z[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_z[i] * fe_0 + tk_xxxx_z[i] * pa_x[i] + 2.0 * ts_xxxxx_z[i] * fz_0;
    }

    // Set up 3-6 components of targeted buffer : HP

    auto tk_xxxxy_x = pbuffer.data(idx_kin_hp + 3);

    auto tk_xxxxy_y = pbuffer.data(idx_kin_hp + 4);

    auto tk_xxxxy_z = pbuffer.data(idx_kin_hp + 5);

#pragma omp simd aligned( \
        pa_y, tk_xxxx_0, tk_xxxx_x, tk_xxxx_y, tk_xxxx_z, tk_xxxxy_x, tk_xxxxy_y, tk_xxxxy_z, ts_xxxxy_x, ts_xxxxy_y, ts_xxxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxy_x[i] = tk_xxxx_x[i] * pa_y[i] + 2.0 * ts_xxxxy_x[i] * fz_0;

        tk_xxxxy_y[i] = tk_xxxx_0[i] * fe_0 + tk_xxxx_y[i] * pa_y[i] + 2.0 * ts_xxxxy_y[i] * fz_0;

        tk_xxxxy_z[i] = tk_xxxx_z[i] * pa_y[i] + 2.0 * ts_xxxxy_z[i] * fz_0;
    }

    // Set up 6-9 components of targeted buffer : HP

    auto tk_xxxxz_x = pbuffer.data(idx_kin_hp + 6);

    auto tk_xxxxz_y = pbuffer.data(idx_kin_hp + 7);

    auto tk_xxxxz_z = pbuffer.data(idx_kin_hp + 8);

#pragma omp simd aligned( \
        pa_z, tk_xxxx_0, tk_xxxx_x, tk_xxxx_y, tk_xxxx_z, tk_xxxxz_x, tk_xxxxz_y, tk_xxxxz_z, ts_xxxxz_x, ts_xxxxz_y, ts_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxz_x[i] = tk_xxxx_x[i] * pa_z[i] + 2.0 * ts_xxxxz_x[i] * fz_0;

        tk_xxxxz_y[i] = tk_xxxx_y[i] * pa_z[i] + 2.0 * ts_xxxxz_y[i] * fz_0;

        tk_xxxxz_z[i] = tk_xxxx_0[i] * fe_0 + tk_xxxx_z[i] * pa_z[i] + 2.0 * ts_xxxxz_z[i] * fz_0;
    }

    // Set up 9-12 components of targeted buffer : HP

    auto tk_xxxyy_x = pbuffer.data(idx_kin_hp + 9);

    auto tk_xxxyy_y = pbuffer.data(idx_kin_hp + 10);

    auto tk_xxxyy_z = pbuffer.data(idx_kin_hp + 11);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             tk_xxx_x,   \
                             tk_xxxy_x,  \
                             tk_xxxyy_x, \
                             tk_xxxyy_y, \
                             tk_xxxyy_z, \
                             tk_xxyy_y,  \
                             tk_xxyy_z,  \
                             tk_xyy_y,   \
                             tk_xyy_z,   \
                             ts_xxx_x,   \
                             ts_xxxyy_x, \
                             ts_xxxyy_y, \
                             ts_xxxyy_z, \
                             ts_xyy_y,   \
                             ts_xyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyy_x[i] = -2.0 * ts_xxx_x[i] * fbe_0 * fz_0 + tk_xxx_x[i] * fe_0 + tk_xxxy_x[i] * pa_y[i] + 2.0 * ts_xxxyy_x[i] * fz_0;

        tk_xxxyy_y[i] = -4.0 * ts_xyy_y[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_y[i] * fe_0 + tk_xxyy_y[i] * pa_x[i] + 2.0 * ts_xxxyy_y[i] * fz_0;

        tk_xxxyy_z[i] = -4.0 * ts_xyy_z[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_z[i] * fe_0 + tk_xxyy_z[i] * pa_x[i] + 2.0 * ts_xxxyy_z[i] * fz_0;
    }

    // Set up 12-15 components of targeted buffer : HP

    auto tk_xxxyz_x = pbuffer.data(idx_kin_hp + 12);

    auto tk_xxxyz_y = pbuffer.data(idx_kin_hp + 13);

    auto tk_xxxyz_z = pbuffer.data(idx_kin_hp + 14);

#pragma omp simd aligned( \
        pa_y, pa_z, tk_xxxy_y, tk_xxxyz_x, tk_xxxyz_y, tk_xxxyz_z, tk_xxxz_x, tk_xxxz_z, ts_xxxyz_x, ts_xxxyz_y, ts_xxxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xxxyz_x[i] = tk_xxxz_x[i] * pa_y[i] + 2.0 * ts_xxxyz_x[i] * fz_0;

        tk_xxxyz_y[i] = tk_xxxy_y[i] * pa_z[i] + 2.0 * ts_xxxyz_y[i] * fz_0;

        tk_xxxyz_z[i] = tk_xxxz_z[i] * pa_y[i] + 2.0 * ts_xxxyz_z[i] * fz_0;
    }

    // Set up 15-18 components of targeted buffer : HP

    auto tk_xxxzz_x = pbuffer.data(idx_kin_hp + 15);

    auto tk_xxxzz_y = pbuffer.data(idx_kin_hp + 16);

    auto tk_xxxzz_z = pbuffer.data(idx_kin_hp + 17);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             tk_xxx_x,   \
                             tk_xxxz_x,  \
                             tk_xxxzz_x, \
                             tk_xxxzz_y, \
                             tk_xxxzz_z, \
                             tk_xxzz_y,  \
                             tk_xxzz_z,  \
                             tk_xzz_y,   \
                             tk_xzz_z,   \
                             ts_xxx_x,   \
                             ts_xxxzz_x, \
                             ts_xxxzz_y, \
                             ts_xxxzz_z, \
                             ts_xzz_y,   \
                             ts_xzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzz_x[i] = -2.0 * ts_xxx_x[i] * fbe_0 * fz_0 + tk_xxx_x[i] * fe_0 + tk_xxxz_x[i] * pa_z[i] + 2.0 * ts_xxxzz_x[i] * fz_0;

        tk_xxxzz_y[i] = -4.0 * ts_xzz_y[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_y[i] * fe_0 + tk_xxzz_y[i] * pa_x[i] + 2.0 * ts_xxxzz_y[i] * fz_0;

        tk_xxxzz_z[i] = -4.0 * ts_xzz_z[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_z[i] * fe_0 + tk_xxzz_z[i] * pa_x[i] + 2.0 * ts_xxxzz_z[i] * fz_0;
    }

    // Set up 18-21 components of targeted buffer : HP

    auto tk_xxyyy_x = pbuffer.data(idx_kin_hp + 18);

    auto tk_xxyyy_y = pbuffer.data(idx_kin_hp + 19);

    auto tk_xxyyy_z = pbuffer.data(idx_kin_hp + 20);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             tk_xxy_x,   \
                             tk_xxyy_x,  \
                             tk_xxyyy_x, \
                             tk_xxyyy_y, \
                             tk_xxyyy_z, \
                             tk_xyyy_y,  \
                             tk_xyyy_z,  \
                             tk_yyy_y,   \
                             tk_yyy_z,   \
                             ts_xxy_x,   \
                             ts_xxyyy_x, \
                             ts_xxyyy_y, \
                             ts_xxyyy_z, \
                             ts_yyy_y,   \
                             ts_yyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyy_x[i] = -4.0 * ts_xxy_x[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_x[i] * fe_0 + tk_xxyy_x[i] * pa_y[i] + 2.0 * ts_xxyyy_x[i] * fz_0;

        tk_xxyyy_y[i] = -2.0 * ts_yyy_y[i] * fbe_0 * fz_0 + tk_yyy_y[i] * fe_0 + tk_xyyy_y[i] * pa_x[i] + 2.0 * ts_xxyyy_y[i] * fz_0;

        tk_xxyyy_z[i] = -2.0 * ts_yyy_z[i] * fbe_0 * fz_0 + tk_yyy_z[i] * fe_0 + tk_xyyy_z[i] * pa_x[i] + 2.0 * ts_xxyyy_z[i] * fz_0;
    }

    // Set up 21-24 components of targeted buffer : HP

    auto tk_xxyyz_x = pbuffer.data(idx_kin_hp + 21);

    auto tk_xxyyz_y = pbuffer.data(idx_kin_hp + 22);

    auto tk_xxyyz_z = pbuffer.data(idx_kin_hp + 23);

#pragma omp simd aligned( \
        pa_z, tk_xxyy_0, tk_xxyy_x, tk_xxyy_y, tk_xxyy_z, tk_xxyyz_x, tk_xxyyz_y, tk_xxyyz_z, ts_xxyyz_x, ts_xxyyz_y, ts_xxyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyz_x[i] = tk_xxyy_x[i] * pa_z[i] + 2.0 * ts_xxyyz_x[i] * fz_0;

        tk_xxyyz_y[i] = tk_xxyy_y[i] * pa_z[i] + 2.0 * ts_xxyyz_y[i] * fz_0;

        tk_xxyyz_z[i] = tk_xxyy_0[i] * fe_0 + tk_xxyy_z[i] * pa_z[i] + 2.0 * ts_xxyyz_z[i] * fz_0;
    }

    // Set up 24-27 components of targeted buffer : HP

    auto tk_xxyzz_x = pbuffer.data(idx_kin_hp + 24);

    auto tk_xxyzz_y = pbuffer.data(idx_kin_hp + 25);

    auto tk_xxyzz_z = pbuffer.data(idx_kin_hp + 26);

#pragma omp simd aligned( \
        pa_y, tk_xxyzz_x, tk_xxyzz_y, tk_xxyzz_z, tk_xxzz_0, tk_xxzz_x, tk_xxzz_y, tk_xxzz_z, ts_xxyzz_x, ts_xxyzz_y, ts_xxyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzz_x[i] = tk_xxzz_x[i] * pa_y[i] + 2.0 * ts_xxyzz_x[i] * fz_0;

        tk_xxyzz_y[i] = tk_xxzz_0[i] * fe_0 + tk_xxzz_y[i] * pa_y[i] + 2.0 * ts_xxyzz_y[i] * fz_0;

        tk_xxyzz_z[i] = tk_xxzz_z[i] * pa_y[i] + 2.0 * ts_xxyzz_z[i] * fz_0;
    }

    // Set up 27-30 components of targeted buffer : HP

    auto tk_xxzzz_x = pbuffer.data(idx_kin_hp + 27);

    auto tk_xxzzz_y = pbuffer.data(idx_kin_hp + 28);

    auto tk_xxzzz_z = pbuffer.data(idx_kin_hp + 29);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             tk_xxz_x,   \
                             tk_xxzz_x,  \
                             tk_xxzzz_x, \
                             tk_xxzzz_y, \
                             tk_xxzzz_z, \
                             tk_xzzz_y,  \
                             tk_xzzz_z,  \
                             tk_zzz_y,   \
                             tk_zzz_z,   \
                             ts_xxz_x,   \
                             ts_xxzzz_x, \
                             ts_xxzzz_y, \
                             ts_xxzzz_z, \
                             ts_zzz_y,   \
                             ts_zzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzz_x[i] = -4.0 * ts_xxz_x[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_x[i] * fe_0 + tk_xxzz_x[i] * pa_z[i] + 2.0 * ts_xxzzz_x[i] * fz_0;

        tk_xxzzz_y[i] = -2.0 * ts_zzz_y[i] * fbe_0 * fz_0 + tk_zzz_y[i] * fe_0 + tk_xzzz_y[i] * pa_x[i] + 2.0 * ts_xxzzz_y[i] * fz_0;

        tk_xxzzz_z[i] = -2.0 * ts_zzz_z[i] * fbe_0 * fz_0 + tk_zzz_z[i] * fe_0 + tk_xzzz_z[i] * pa_x[i] + 2.0 * ts_xxzzz_z[i] * fz_0;
    }

    // Set up 30-33 components of targeted buffer : HP

    auto tk_xyyyy_x = pbuffer.data(idx_kin_hp + 30);

    auto tk_xyyyy_y = pbuffer.data(idx_kin_hp + 31);

    auto tk_xyyyy_z = pbuffer.data(idx_kin_hp + 32);

#pragma omp simd aligned( \
        pa_x, tk_xyyyy_x, tk_xyyyy_y, tk_xyyyy_z, tk_yyyy_0, tk_yyyy_x, tk_yyyy_y, tk_yyyy_z, ts_xyyyy_x, ts_xyyyy_y, ts_xyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyy_x[i] = tk_yyyy_0[i] * fe_0 + tk_yyyy_x[i] * pa_x[i] + 2.0 * ts_xyyyy_x[i] * fz_0;

        tk_xyyyy_y[i] = tk_yyyy_y[i] * pa_x[i] + 2.0 * ts_xyyyy_y[i] * fz_0;

        tk_xyyyy_z[i] = tk_yyyy_z[i] * pa_x[i] + 2.0 * ts_xyyyy_z[i] * fz_0;
    }

    // Set up 33-36 components of targeted buffer : HP

    auto tk_xyyyz_x = pbuffer.data(idx_kin_hp + 33);

    auto tk_xyyyz_y = pbuffer.data(idx_kin_hp + 34);

    auto tk_xyyyz_z = pbuffer.data(idx_kin_hp + 35);

#pragma omp simd aligned( \
        pa_x, pa_z, tk_xyyy_x, tk_xyyyz_x, tk_xyyyz_y, tk_xyyyz_z, tk_yyyz_y, tk_yyyz_z, ts_xyyyz_x, ts_xyyyz_y, ts_xyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyyyz_x[i] = tk_xyyy_x[i] * pa_z[i] + 2.0 * ts_xyyyz_x[i] * fz_0;

        tk_xyyyz_y[i] = tk_yyyz_y[i] * pa_x[i] + 2.0 * ts_xyyyz_y[i] * fz_0;

        tk_xyyyz_z[i] = tk_yyyz_z[i] * pa_x[i] + 2.0 * ts_xyyyz_z[i] * fz_0;
    }

    // Set up 36-39 components of targeted buffer : HP

    auto tk_xyyzz_x = pbuffer.data(idx_kin_hp + 36);

    auto tk_xyyzz_y = pbuffer.data(idx_kin_hp + 37);

    auto tk_xyyzz_z = pbuffer.data(idx_kin_hp + 38);

#pragma omp simd aligned( \
        pa_x, tk_xyyzz_x, tk_xyyzz_y, tk_xyyzz_z, tk_yyzz_0, tk_yyzz_x, tk_yyzz_y, tk_yyzz_z, ts_xyyzz_x, ts_xyyzz_y, ts_xyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzz_x[i] = tk_yyzz_0[i] * fe_0 + tk_yyzz_x[i] * pa_x[i] + 2.0 * ts_xyyzz_x[i] * fz_0;

        tk_xyyzz_y[i] = tk_yyzz_y[i] * pa_x[i] + 2.0 * ts_xyyzz_y[i] * fz_0;

        tk_xyyzz_z[i] = tk_yyzz_z[i] * pa_x[i] + 2.0 * ts_xyyzz_z[i] * fz_0;
    }

    // Set up 39-42 components of targeted buffer : HP

    auto tk_xyzzz_x = pbuffer.data(idx_kin_hp + 39);

    auto tk_xyzzz_y = pbuffer.data(idx_kin_hp + 40);

    auto tk_xyzzz_z = pbuffer.data(idx_kin_hp + 41);

#pragma omp simd aligned( \
        pa_x, pa_y, tk_xyzzz_x, tk_xyzzz_y, tk_xyzzz_z, tk_xzzz_x, tk_yzzz_y, tk_yzzz_z, ts_xyzzz_x, ts_xyzzz_y, ts_xyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyzzz_x[i] = tk_xzzz_x[i] * pa_y[i] + 2.0 * ts_xyzzz_x[i] * fz_0;

        tk_xyzzz_y[i] = tk_yzzz_y[i] * pa_x[i] + 2.0 * ts_xyzzz_y[i] * fz_0;

        tk_xyzzz_z[i] = tk_yzzz_z[i] * pa_x[i] + 2.0 * ts_xyzzz_z[i] * fz_0;
    }

    // Set up 42-45 components of targeted buffer : HP

    auto tk_xzzzz_x = pbuffer.data(idx_kin_hp + 42);

    auto tk_xzzzz_y = pbuffer.data(idx_kin_hp + 43);

    auto tk_xzzzz_z = pbuffer.data(idx_kin_hp + 44);

#pragma omp simd aligned( \
        pa_x, tk_xzzzz_x, tk_xzzzz_y, tk_xzzzz_z, tk_zzzz_0, tk_zzzz_x, tk_zzzz_y, tk_zzzz_z, ts_xzzzz_x, ts_xzzzz_y, ts_xzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzz_x[i] = tk_zzzz_0[i] * fe_0 + tk_zzzz_x[i] * pa_x[i] + 2.0 * ts_xzzzz_x[i] * fz_0;

        tk_xzzzz_y[i] = tk_zzzz_y[i] * pa_x[i] + 2.0 * ts_xzzzz_y[i] * fz_0;

        tk_xzzzz_z[i] = tk_zzzz_z[i] * pa_x[i] + 2.0 * ts_xzzzz_z[i] * fz_0;
    }

    // Set up 45-48 components of targeted buffer : HP

    auto tk_yyyyy_x = pbuffer.data(idx_kin_hp + 45);

    auto tk_yyyyy_y = pbuffer.data(idx_kin_hp + 46);

    auto tk_yyyyy_z = pbuffer.data(idx_kin_hp + 47);

#pragma omp simd aligned(pa_y,           \
                             tk_yyy_x,   \
                             tk_yyy_y,   \
                             tk_yyy_z,   \
                             tk_yyyy_0,  \
                             tk_yyyy_x,  \
                             tk_yyyy_y,  \
                             tk_yyyy_z,  \
                             tk_yyyyy_x, \
                             tk_yyyyy_y, \
                             tk_yyyyy_z, \
                             ts_yyy_x,   \
                             ts_yyy_y,   \
                             ts_yyy_z,   \
                             ts_yyyyy_x, \
                             ts_yyyyy_y, \
                             ts_yyyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyy_x[i] = -8.0 * ts_yyy_x[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_x[i] * fe_0 + tk_yyyy_x[i] * pa_y[i] + 2.0 * ts_yyyyy_x[i] * fz_0;

        tk_yyyyy_y[i] =
            -8.0 * ts_yyy_y[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_y[i] * fe_0 + tk_yyyy_0[i] * fe_0 + tk_yyyy_y[i] * pa_y[i] + 2.0 * ts_yyyyy_y[i] * fz_0;

        tk_yyyyy_z[i] = -8.0 * ts_yyy_z[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_z[i] * fe_0 + tk_yyyy_z[i] * pa_y[i] + 2.0 * ts_yyyyy_z[i] * fz_0;
    }

    // Set up 48-51 components of targeted buffer : HP

    auto tk_yyyyz_x = pbuffer.data(idx_kin_hp + 48);

    auto tk_yyyyz_y = pbuffer.data(idx_kin_hp + 49);

    auto tk_yyyyz_z = pbuffer.data(idx_kin_hp + 50);

#pragma omp simd aligned( \
        pa_z, tk_yyyy_0, tk_yyyy_x, tk_yyyy_y, tk_yyyy_z, tk_yyyyz_x, tk_yyyyz_y, tk_yyyyz_z, ts_yyyyz_x, ts_yyyyz_y, ts_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyz_x[i] = tk_yyyy_x[i] * pa_z[i] + 2.0 * ts_yyyyz_x[i] * fz_0;

        tk_yyyyz_y[i] = tk_yyyy_y[i] * pa_z[i] + 2.0 * ts_yyyyz_y[i] * fz_0;

        tk_yyyyz_z[i] = tk_yyyy_0[i] * fe_0 + tk_yyyy_z[i] * pa_z[i] + 2.0 * ts_yyyyz_z[i] * fz_0;
    }

    // Set up 51-54 components of targeted buffer : HP

    auto tk_yyyzz_x = pbuffer.data(idx_kin_hp + 51);

    auto tk_yyyzz_y = pbuffer.data(idx_kin_hp + 52);

    auto tk_yyyzz_z = pbuffer.data(idx_kin_hp + 53);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             tk_yyy_y,   \
                             tk_yyyz_y,  \
                             tk_yyyzz_x, \
                             tk_yyyzz_y, \
                             tk_yyyzz_z, \
                             tk_yyzz_x,  \
                             tk_yyzz_z,  \
                             tk_yzz_x,   \
                             tk_yzz_z,   \
                             ts_yyy_y,   \
                             ts_yyyzz_x, \
                             ts_yyyzz_y, \
                             ts_yyyzz_z, \
                             ts_yzz_x,   \
                             ts_yzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzz_x[i] = -4.0 * ts_yzz_x[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_x[i] * fe_0 + tk_yyzz_x[i] * pa_y[i] + 2.0 * ts_yyyzz_x[i] * fz_0;

        tk_yyyzz_y[i] = -2.0 * ts_yyy_y[i] * fbe_0 * fz_0 + tk_yyy_y[i] * fe_0 + tk_yyyz_y[i] * pa_z[i] + 2.0 * ts_yyyzz_y[i] * fz_0;

        tk_yyyzz_z[i] = -4.0 * ts_yzz_z[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_z[i] * fe_0 + tk_yyzz_z[i] * pa_y[i] + 2.0 * ts_yyyzz_z[i] * fz_0;
    }

    // Set up 54-57 components of targeted buffer : HP

    auto tk_yyzzz_x = pbuffer.data(idx_kin_hp + 54);

    auto tk_yyzzz_y = pbuffer.data(idx_kin_hp + 55);

    auto tk_yyzzz_z = pbuffer.data(idx_kin_hp + 56);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             tk_yyz_y,   \
                             tk_yyzz_y,  \
                             tk_yyzzz_x, \
                             tk_yyzzz_y, \
                             tk_yyzzz_z, \
                             tk_yzzz_x,  \
                             tk_yzzz_z,  \
                             tk_zzz_x,   \
                             tk_zzz_z,   \
                             ts_yyz_y,   \
                             ts_yyzzz_x, \
                             ts_yyzzz_y, \
                             ts_yyzzz_z, \
                             ts_zzz_x,   \
                             ts_zzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzz_x[i] = -2.0 * ts_zzz_x[i] * fbe_0 * fz_0 + tk_zzz_x[i] * fe_0 + tk_yzzz_x[i] * pa_y[i] + 2.0 * ts_yyzzz_x[i] * fz_0;

        tk_yyzzz_y[i] = -4.0 * ts_yyz_y[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_y[i] * fe_0 + tk_yyzz_y[i] * pa_z[i] + 2.0 * ts_yyzzz_y[i] * fz_0;

        tk_yyzzz_z[i] = -2.0 * ts_zzz_z[i] * fbe_0 * fz_0 + tk_zzz_z[i] * fe_0 + tk_yzzz_z[i] * pa_y[i] + 2.0 * ts_yyzzz_z[i] * fz_0;
    }

    // Set up 57-60 components of targeted buffer : HP

    auto tk_yzzzz_x = pbuffer.data(idx_kin_hp + 57);

    auto tk_yzzzz_y = pbuffer.data(idx_kin_hp + 58);

    auto tk_yzzzz_z = pbuffer.data(idx_kin_hp + 59);

#pragma omp simd aligned( \
        pa_y, tk_yzzzz_x, tk_yzzzz_y, tk_yzzzz_z, tk_zzzz_0, tk_zzzz_x, tk_zzzz_y, tk_zzzz_z, ts_yzzzz_x, ts_yzzzz_y, ts_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzz_x[i] = tk_zzzz_x[i] * pa_y[i] + 2.0 * ts_yzzzz_x[i] * fz_0;

        tk_yzzzz_y[i] = tk_zzzz_0[i] * fe_0 + tk_zzzz_y[i] * pa_y[i] + 2.0 * ts_yzzzz_y[i] * fz_0;

        tk_yzzzz_z[i] = tk_zzzz_z[i] * pa_y[i] + 2.0 * ts_yzzzz_z[i] * fz_0;
    }

    // Set up 60-63 components of targeted buffer : HP

    auto tk_zzzzz_x = pbuffer.data(idx_kin_hp + 60);

    auto tk_zzzzz_y = pbuffer.data(idx_kin_hp + 61);

    auto tk_zzzzz_z = pbuffer.data(idx_kin_hp + 62);

#pragma omp simd aligned(pa_z,           \
                             tk_zzz_x,   \
                             tk_zzz_y,   \
                             tk_zzz_z,   \
                             tk_zzzz_0,  \
                             tk_zzzz_x,  \
                             tk_zzzz_y,  \
                             tk_zzzz_z,  \
                             tk_zzzzz_x, \
                             tk_zzzzz_y, \
                             tk_zzzzz_z, \
                             ts_zzz_x,   \
                             ts_zzz_y,   \
                             ts_zzz_z,   \
                             ts_zzzzz_x, \
                             ts_zzzzz_y, \
                             ts_zzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzz_x[i] = -8.0 * ts_zzz_x[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_x[i] * fe_0 + tk_zzzz_x[i] * pa_z[i] + 2.0 * ts_zzzzz_x[i] * fz_0;

        tk_zzzzz_y[i] = -8.0 * ts_zzz_y[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_y[i] * fe_0 + tk_zzzz_y[i] * pa_z[i] + 2.0 * ts_zzzzz_y[i] * fz_0;

        tk_zzzzz_z[i] =
            -8.0 * ts_zzz_z[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_z[i] * fe_0 + tk_zzzz_0[i] * fe_0 + tk_zzzz_z[i] * pa_z[i] + 2.0 * ts_zzzzz_z[i] * fz_0;
    }
}

}  // namespace kinrec
