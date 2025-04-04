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

#include "KineticEnergyPrimRecIP.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ip(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ip,
                            const size_t              idx_ovl_gp,
                            const size_t              idx_kin_gp,
                            const size_t              idx_kin_hs,
                            const size_t              idx_kin_hp,
                            const size_t              idx_ovl_ip,
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

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_ovl_gp);

    auto ts_xxxx_y = pbuffer.data(idx_ovl_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_ovl_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_ovl_gp + 3);

    auto ts_xxxz_x = pbuffer.data(idx_ovl_gp + 6);

    auto ts_xxyy_x = pbuffer.data(idx_ovl_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_ovl_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_ovl_gp + 11);

    auto ts_xxzz_x = pbuffer.data(idx_ovl_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_ovl_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_ovl_gp + 17);

    auto ts_xyyy_y = pbuffer.data(idx_ovl_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_ovl_gp + 20);

    auto ts_xzzz_y = pbuffer.data(idx_ovl_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_ovl_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_ovl_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_ovl_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_ovl_gp + 32);

    auto ts_yyyz_y = pbuffer.data(idx_ovl_gp + 34);

    auto ts_yyzz_x = pbuffer.data(idx_ovl_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_ovl_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_ovl_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_ovl_gp + 39);

    auto ts_yzzz_z = pbuffer.data(idx_ovl_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_ovl_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_ovl_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_ovl_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto tk_xxxx_x = pbuffer.data(idx_kin_gp);

    auto tk_xxxx_y = pbuffer.data(idx_kin_gp + 1);

    auto tk_xxxx_z = pbuffer.data(idx_kin_gp + 2);

    auto tk_xxxy_x = pbuffer.data(idx_kin_gp + 3);

    auto tk_xxxz_x = pbuffer.data(idx_kin_gp + 6);

    auto tk_xxyy_x = pbuffer.data(idx_kin_gp + 9);

    auto tk_xxyy_y = pbuffer.data(idx_kin_gp + 10);

    auto tk_xxyy_z = pbuffer.data(idx_kin_gp + 11);

    auto tk_xxzz_x = pbuffer.data(idx_kin_gp + 15);

    auto tk_xxzz_y = pbuffer.data(idx_kin_gp + 16);

    auto tk_xxzz_z = pbuffer.data(idx_kin_gp + 17);

    auto tk_xyyy_y = pbuffer.data(idx_kin_gp + 19);

    auto tk_xyyy_z = pbuffer.data(idx_kin_gp + 20);

    auto tk_xzzz_y = pbuffer.data(idx_kin_gp + 28);

    auto tk_xzzz_z = pbuffer.data(idx_kin_gp + 29);

    auto tk_yyyy_x = pbuffer.data(idx_kin_gp + 30);

    auto tk_yyyy_y = pbuffer.data(idx_kin_gp + 31);

    auto tk_yyyy_z = pbuffer.data(idx_kin_gp + 32);

    auto tk_yyyz_y = pbuffer.data(idx_kin_gp + 34);

    auto tk_yyzz_x = pbuffer.data(idx_kin_gp + 36);

    auto tk_yyzz_y = pbuffer.data(idx_kin_gp + 37);

    auto tk_yyzz_z = pbuffer.data(idx_kin_gp + 38);

    auto tk_yzzz_x = pbuffer.data(idx_kin_gp + 39);

    auto tk_yzzz_z = pbuffer.data(idx_kin_gp + 41);

    auto tk_zzzz_x = pbuffer.data(idx_kin_gp + 42);

    auto tk_zzzz_y = pbuffer.data(idx_kin_gp + 43);

    auto tk_zzzz_z = pbuffer.data(idx_kin_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto tk_xxxxx_0 = pbuffer.data(idx_kin_hs);

    auto tk_xxxyy_0 = pbuffer.data(idx_kin_hs + 3);

    auto tk_xxxzz_0 = pbuffer.data(idx_kin_hs + 5);

    auto tk_xxyyy_0 = pbuffer.data(idx_kin_hs + 6);

    auto tk_xxzzz_0 = pbuffer.data(idx_kin_hs + 9);

    auto tk_yyyyy_0 = pbuffer.data(idx_kin_hs + 15);

    auto tk_yyyzz_0 = pbuffer.data(idx_kin_hs + 17);

    auto tk_yyzzz_0 = pbuffer.data(idx_kin_hs + 18);

    auto tk_zzzzz_0 = pbuffer.data(idx_kin_hs + 20);

    // Set up components of auxiliary buffer : HP

    auto tk_xxxxx_x = pbuffer.data(idx_kin_hp);

    auto tk_xxxxx_y = pbuffer.data(idx_kin_hp + 1);

    auto tk_xxxxx_z = pbuffer.data(idx_kin_hp + 2);

    auto tk_xxxxy_x = pbuffer.data(idx_kin_hp + 3);

    auto tk_xxxxy_y = pbuffer.data(idx_kin_hp + 4);

    auto tk_xxxxz_x = pbuffer.data(idx_kin_hp + 6);

    auto tk_xxxxz_z = pbuffer.data(idx_kin_hp + 8);

    auto tk_xxxyy_x = pbuffer.data(idx_kin_hp + 9);

    auto tk_xxxyy_y = pbuffer.data(idx_kin_hp + 10);

    auto tk_xxxyy_z = pbuffer.data(idx_kin_hp + 11);

    auto tk_xxxzz_x = pbuffer.data(idx_kin_hp + 15);

    auto tk_xxxzz_y = pbuffer.data(idx_kin_hp + 16);

    auto tk_xxxzz_z = pbuffer.data(idx_kin_hp + 17);

    auto tk_xxyyy_x = pbuffer.data(idx_kin_hp + 18);

    auto tk_xxyyy_y = pbuffer.data(idx_kin_hp + 19);

    auto tk_xxyyy_z = pbuffer.data(idx_kin_hp + 20);

    auto tk_xxyzz_x = pbuffer.data(idx_kin_hp + 24);

    auto tk_xxzzz_x = pbuffer.data(idx_kin_hp + 27);

    auto tk_xxzzz_y = pbuffer.data(idx_kin_hp + 28);

    auto tk_xxzzz_z = pbuffer.data(idx_kin_hp + 29);

    auto tk_xyyyy_x = pbuffer.data(idx_kin_hp + 30);

    auto tk_xyyyy_y = pbuffer.data(idx_kin_hp + 31);

    auto tk_xyyyy_z = pbuffer.data(idx_kin_hp + 32);

    auto tk_xyyzz_y = pbuffer.data(idx_kin_hp + 37);

    auto tk_xyyzz_z = pbuffer.data(idx_kin_hp + 38);

    auto tk_xzzzz_x = pbuffer.data(idx_kin_hp + 42);

    auto tk_xzzzz_y = pbuffer.data(idx_kin_hp + 43);

    auto tk_xzzzz_z = pbuffer.data(idx_kin_hp + 44);

    auto tk_yyyyy_x = pbuffer.data(idx_kin_hp + 45);

    auto tk_yyyyy_y = pbuffer.data(idx_kin_hp + 46);

    auto tk_yyyyy_z = pbuffer.data(idx_kin_hp + 47);

    auto tk_yyyyz_y = pbuffer.data(idx_kin_hp + 49);

    auto tk_yyyyz_z = pbuffer.data(idx_kin_hp + 50);

    auto tk_yyyzz_x = pbuffer.data(idx_kin_hp + 51);

    auto tk_yyyzz_y = pbuffer.data(idx_kin_hp + 52);

    auto tk_yyyzz_z = pbuffer.data(idx_kin_hp + 53);

    auto tk_yyzzz_x = pbuffer.data(idx_kin_hp + 54);

    auto tk_yyzzz_y = pbuffer.data(idx_kin_hp + 55);

    auto tk_yyzzz_z = pbuffer.data(idx_kin_hp + 56);

    auto tk_yzzzz_x = pbuffer.data(idx_kin_hp + 57);

    auto tk_yzzzz_y = pbuffer.data(idx_kin_hp + 58);

    auto tk_yzzzz_z = pbuffer.data(idx_kin_hp + 59);

    auto tk_zzzzz_x = pbuffer.data(idx_kin_hp + 60);

    auto tk_zzzzz_y = pbuffer.data(idx_kin_hp + 61);

    auto tk_zzzzz_z = pbuffer.data(idx_kin_hp + 62);

    // Set up components of auxiliary buffer : IP

    auto ts_xxxxxx_x = pbuffer.data(idx_ovl_ip);

    auto ts_xxxxxx_y = pbuffer.data(idx_ovl_ip + 1);

    auto ts_xxxxxx_z = pbuffer.data(idx_ovl_ip + 2);

    auto ts_xxxxxy_x = pbuffer.data(idx_ovl_ip + 3);

    auto ts_xxxxxy_y = pbuffer.data(idx_ovl_ip + 4);

    auto ts_xxxxxy_z = pbuffer.data(idx_ovl_ip + 5);

    auto ts_xxxxxz_x = pbuffer.data(idx_ovl_ip + 6);

    auto ts_xxxxxz_y = pbuffer.data(idx_ovl_ip + 7);

    auto ts_xxxxxz_z = pbuffer.data(idx_ovl_ip + 8);

    auto ts_xxxxyy_x = pbuffer.data(idx_ovl_ip + 9);

    auto ts_xxxxyy_y = pbuffer.data(idx_ovl_ip + 10);

    auto ts_xxxxyy_z = pbuffer.data(idx_ovl_ip + 11);

    auto ts_xxxxyz_x = pbuffer.data(idx_ovl_ip + 12);

    auto ts_xxxxyz_y = pbuffer.data(idx_ovl_ip + 13);

    auto ts_xxxxyz_z = pbuffer.data(idx_ovl_ip + 14);

    auto ts_xxxxzz_x = pbuffer.data(idx_ovl_ip + 15);

    auto ts_xxxxzz_y = pbuffer.data(idx_ovl_ip + 16);

    auto ts_xxxxzz_z = pbuffer.data(idx_ovl_ip + 17);

    auto ts_xxxyyy_x = pbuffer.data(idx_ovl_ip + 18);

    auto ts_xxxyyy_y = pbuffer.data(idx_ovl_ip + 19);

    auto ts_xxxyyy_z = pbuffer.data(idx_ovl_ip + 20);

    auto ts_xxxyyz_x = pbuffer.data(idx_ovl_ip + 21);

    auto ts_xxxyyz_y = pbuffer.data(idx_ovl_ip + 22);

    auto ts_xxxyyz_z = pbuffer.data(idx_ovl_ip + 23);

    auto ts_xxxyzz_x = pbuffer.data(idx_ovl_ip + 24);

    auto ts_xxxyzz_y = pbuffer.data(idx_ovl_ip + 25);

    auto ts_xxxyzz_z = pbuffer.data(idx_ovl_ip + 26);

    auto ts_xxxzzz_x = pbuffer.data(idx_ovl_ip + 27);

    auto ts_xxxzzz_y = pbuffer.data(idx_ovl_ip + 28);

    auto ts_xxxzzz_z = pbuffer.data(idx_ovl_ip + 29);

    auto ts_xxyyyy_x = pbuffer.data(idx_ovl_ip + 30);

    auto ts_xxyyyy_y = pbuffer.data(idx_ovl_ip + 31);

    auto ts_xxyyyy_z = pbuffer.data(idx_ovl_ip + 32);

    auto ts_xxyyyz_x = pbuffer.data(idx_ovl_ip + 33);

    auto ts_xxyyyz_y = pbuffer.data(idx_ovl_ip + 34);

    auto ts_xxyyyz_z = pbuffer.data(idx_ovl_ip + 35);

    auto ts_xxyyzz_x = pbuffer.data(idx_ovl_ip + 36);

    auto ts_xxyyzz_y = pbuffer.data(idx_ovl_ip + 37);

    auto ts_xxyyzz_z = pbuffer.data(idx_ovl_ip + 38);

    auto ts_xxyzzz_x = pbuffer.data(idx_ovl_ip + 39);

    auto ts_xxyzzz_y = pbuffer.data(idx_ovl_ip + 40);

    auto ts_xxyzzz_z = pbuffer.data(idx_ovl_ip + 41);

    auto ts_xxzzzz_x = pbuffer.data(idx_ovl_ip + 42);

    auto ts_xxzzzz_y = pbuffer.data(idx_ovl_ip + 43);

    auto ts_xxzzzz_z = pbuffer.data(idx_ovl_ip + 44);

    auto ts_xyyyyy_x = pbuffer.data(idx_ovl_ip + 45);

    auto ts_xyyyyy_y = pbuffer.data(idx_ovl_ip + 46);

    auto ts_xyyyyy_z = pbuffer.data(idx_ovl_ip + 47);

    auto ts_xyyyyz_x = pbuffer.data(idx_ovl_ip + 48);

    auto ts_xyyyyz_y = pbuffer.data(idx_ovl_ip + 49);

    auto ts_xyyyyz_z = pbuffer.data(idx_ovl_ip + 50);

    auto ts_xyyyzz_x = pbuffer.data(idx_ovl_ip + 51);

    auto ts_xyyyzz_y = pbuffer.data(idx_ovl_ip + 52);

    auto ts_xyyyzz_z = pbuffer.data(idx_ovl_ip + 53);

    auto ts_xyyzzz_x = pbuffer.data(idx_ovl_ip + 54);

    auto ts_xyyzzz_y = pbuffer.data(idx_ovl_ip + 55);

    auto ts_xyyzzz_z = pbuffer.data(idx_ovl_ip + 56);

    auto ts_xyzzzz_x = pbuffer.data(idx_ovl_ip + 57);

    auto ts_xyzzzz_y = pbuffer.data(idx_ovl_ip + 58);

    auto ts_xyzzzz_z = pbuffer.data(idx_ovl_ip + 59);

    auto ts_xzzzzz_x = pbuffer.data(idx_ovl_ip + 60);

    auto ts_xzzzzz_y = pbuffer.data(idx_ovl_ip + 61);

    auto ts_xzzzzz_z = pbuffer.data(idx_ovl_ip + 62);

    auto ts_yyyyyy_x = pbuffer.data(idx_ovl_ip + 63);

    auto ts_yyyyyy_y = pbuffer.data(idx_ovl_ip + 64);

    auto ts_yyyyyy_z = pbuffer.data(idx_ovl_ip + 65);

    auto ts_yyyyyz_x = pbuffer.data(idx_ovl_ip + 66);

    auto ts_yyyyyz_y = pbuffer.data(idx_ovl_ip + 67);

    auto ts_yyyyyz_z = pbuffer.data(idx_ovl_ip + 68);

    auto ts_yyyyzz_x = pbuffer.data(idx_ovl_ip + 69);

    auto ts_yyyyzz_y = pbuffer.data(idx_ovl_ip + 70);

    auto ts_yyyyzz_z = pbuffer.data(idx_ovl_ip + 71);

    auto ts_yyyzzz_x = pbuffer.data(idx_ovl_ip + 72);

    auto ts_yyyzzz_y = pbuffer.data(idx_ovl_ip + 73);

    auto ts_yyyzzz_z = pbuffer.data(idx_ovl_ip + 74);

    auto ts_yyzzzz_x = pbuffer.data(idx_ovl_ip + 75);

    auto ts_yyzzzz_y = pbuffer.data(idx_ovl_ip + 76);

    auto ts_yyzzzz_z = pbuffer.data(idx_ovl_ip + 77);

    auto ts_yzzzzz_x = pbuffer.data(idx_ovl_ip + 78);

    auto ts_yzzzzz_y = pbuffer.data(idx_ovl_ip + 79);

    auto ts_yzzzzz_z = pbuffer.data(idx_ovl_ip + 80);

    auto ts_zzzzzz_x = pbuffer.data(idx_ovl_ip + 81);

    auto ts_zzzzzz_y = pbuffer.data(idx_ovl_ip + 82);

    auto ts_zzzzzz_z = pbuffer.data(idx_ovl_ip + 83);

    // Set up 0-3 components of targeted buffer : IP

    auto tk_xxxxxx_x = pbuffer.data(idx_kin_ip);

    auto tk_xxxxxx_y = pbuffer.data(idx_kin_ip + 1);

    auto tk_xxxxxx_z = pbuffer.data(idx_kin_ip + 2);

#pragma omp simd aligned(pa_x,            \
                             tk_xxxx_x,   \
                             tk_xxxx_y,   \
                             tk_xxxx_z,   \
                             tk_xxxxx_0,  \
                             tk_xxxxx_x,  \
                             tk_xxxxx_y,  \
                             tk_xxxxx_z,  \
                             tk_xxxxxx_x, \
                             tk_xxxxxx_y, \
                             tk_xxxxxx_z, \
                             ts_xxxx_x,   \
                             ts_xxxx_y,   \
                             ts_xxxx_z,   \
                             ts_xxxxxx_x, \
                             ts_xxxxxx_y, \
                             ts_xxxxxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxxx_x[i] = -10.0 * ts_xxxx_x[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_x[i] * fe_0 + tk_xxxxx_0[i] * fe_0 + tk_xxxxx_x[i] * pa_x[i] +
                         2.0 * ts_xxxxxx_x[i] * fz_0;

        tk_xxxxxx_y[i] = -10.0 * ts_xxxx_y[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_y[i] * fe_0 + tk_xxxxx_y[i] * pa_x[i] + 2.0 * ts_xxxxxx_y[i] * fz_0;

        tk_xxxxxx_z[i] = -10.0 * ts_xxxx_z[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_z[i] * fe_0 + tk_xxxxx_z[i] * pa_x[i] + 2.0 * ts_xxxxxx_z[i] * fz_0;
    }

    // Set up 3-6 components of targeted buffer : IP

    auto tk_xxxxxy_x = pbuffer.data(idx_kin_ip + 3);

    auto tk_xxxxxy_y = pbuffer.data(idx_kin_ip + 4);

    auto tk_xxxxxy_z = pbuffer.data(idx_kin_ip + 5);

#pragma omp simd aligned(pa_y,            \
                             tk_xxxxx_0,  \
                             tk_xxxxx_x,  \
                             tk_xxxxx_y,  \
                             tk_xxxxx_z,  \
                             tk_xxxxxy_x, \
                             tk_xxxxxy_y, \
                             tk_xxxxxy_z, \
                             ts_xxxxxy_x, \
                             ts_xxxxxy_y, \
                             ts_xxxxxy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxy_x[i] = tk_xxxxx_x[i] * pa_y[i] + 2.0 * ts_xxxxxy_x[i] * fz_0;

        tk_xxxxxy_y[i] = tk_xxxxx_0[i] * fe_0 + tk_xxxxx_y[i] * pa_y[i] + 2.0 * ts_xxxxxy_y[i] * fz_0;

        tk_xxxxxy_z[i] = tk_xxxxx_z[i] * pa_y[i] + 2.0 * ts_xxxxxy_z[i] * fz_0;
    }

    // Set up 6-9 components of targeted buffer : IP

    auto tk_xxxxxz_x = pbuffer.data(idx_kin_ip + 6);

    auto tk_xxxxxz_y = pbuffer.data(idx_kin_ip + 7);

    auto tk_xxxxxz_z = pbuffer.data(idx_kin_ip + 8);

#pragma omp simd aligned(pa_z,            \
                             tk_xxxxx_0,  \
                             tk_xxxxx_x,  \
                             tk_xxxxx_y,  \
                             tk_xxxxx_z,  \
                             tk_xxxxxz_x, \
                             tk_xxxxxz_y, \
                             tk_xxxxxz_z, \
                             ts_xxxxxz_x, \
                             ts_xxxxxz_y, \
                             ts_xxxxxz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxz_x[i] = tk_xxxxx_x[i] * pa_z[i] + 2.0 * ts_xxxxxz_x[i] * fz_0;

        tk_xxxxxz_y[i] = tk_xxxxx_y[i] * pa_z[i] + 2.0 * ts_xxxxxz_y[i] * fz_0;

        tk_xxxxxz_z[i] = tk_xxxxx_0[i] * fe_0 + tk_xxxxx_z[i] * pa_z[i] + 2.0 * ts_xxxxxz_z[i] * fz_0;
    }

    // Set up 9-12 components of targeted buffer : IP

    auto tk_xxxxyy_x = pbuffer.data(idx_kin_ip + 9);

    auto tk_xxxxyy_y = pbuffer.data(idx_kin_ip + 10);

    auto tk_xxxxyy_z = pbuffer.data(idx_kin_ip + 11);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tk_xxxx_x,   \
                             tk_xxxxy_x,  \
                             tk_xxxxyy_x, \
                             tk_xxxxyy_y, \
                             tk_xxxxyy_z, \
                             tk_xxxyy_y,  \
                             tk_xxxyy_z,  \
                             tk_xxyy_y,   \
                             tk_xxyy_z,   \
                             ts_xxxx_x,   \
                             ts_xxxxyy_x, \
                             ts_xxxxyy_y, \
                             ts_xxxxyy_z, \
                             ts_xxyy_y,   \
                             ts_xxyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxyy_x[i] = -2.0 * ts_xxxx_x[i] * fbe_0 * fz_0 + tk_xxxx_x[i] * fe_0 + tk_xxxxy_x[i] * pa_y[i] + 2.0 * ts_xxxxyy_x[i] * fz_0;

        tk_xxxxyy_y[i] = -6.0 * ts_xxyy_y[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_y[i] * fe_0 + tk_xxxyy_y[i] * pa_x[i] + 2.0 * ts_xxxxyy_y[i] * fz_0;

        tk_xxxxyy_z[i] = -6.0 * ts_xxyy_z[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_z[i] * fe_0 + tk_xxxyy_z[i] * pa_x[i] + 2.0 * ts_xxxxyy_z[i] * fz_0;
    }

    // Set up 12-15 components of targeted buffer : IP

    auto tk_xxxxyz_x = pbuffer.data(idx_kin_ip + 12);

    auto tk_xxxxyz_y = pbuffer.data(idx_kin_ip + 13);

    auto tk_xxxxyz_z = pbuffer.data(idx_kin_ip + 14);

#pragma omp simd aligned( \
        pa_y, pa_z, tk_xxxxy_y, tk_xxxxyz_x, tk_xxxxyz_y, tk_xxxxyz_z, tk_xxxxz_x, tk_xxxxz_z, ts_xxxxyz_x, ts_xxxxyz_y, ts_xxxxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xxxxyz_x[i] = tk_xxxxz_x[i] * pa_y[i] + 2.0 * ts_xxxxyz_x[i] * fz_0;

        tk_xxxxyz_y[i] = tk_xxxxy_y[i] * pa_z[i] + 2.0 * ts_xxxxyz_y[i] * fz_0;

        tk_xxxxyz_z[i] = tk_xxxxz_z[i] * pa_y[i] + 2.0 * ts_xxxxyz_z[i] * fz_0;
    }

    // Set up 15-18 components of targeted buffer : IP

    auto tk_xxxxzz_x = pbuffer.data(idx_kin_ip + 15);

    auto tk_xxxxzz_y = pbuffer.data(idx_kin_ip + 16);

    auto tk_xxxxzz_z = pbuffer.data(idx_kin_ip + 17);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tk_xxxx_x,   \
                             tk_xxxxz_x,  \
                             tk_xxxxzz_x, \
                             tk_xxxxzz_y, \
                             tk_xxxxzz_z, \
                             tk_xxxzz_y,  \
                             tk_xxxzz_z,  \
                             tk_xxzz_y,   \
                             tk_xxzz_z,   \
                             ts_xxxx_x,   \
                             ts_xxxxzz_x, \
                             ts_xxxxzz_y, \
                             ts_xxxxzz_z, \
                             ts_xxzz_y,   \
                             ts_xxzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxzz_x[i] = -2.0 * ts_xxxx_x[i] * fbe_0 * fz_0 + tk_xxxx_x[i] * fe_0 + tk_xxxxz_x[i] * pa_z[i] + 2.0 * ts_xxxxzz_x[i] * fz_0;

        tk_xxxxzz_y[i] = -6.0 * ts_xxzz_y[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_y[i] * fe_0 + tk_xxxzz_y[i] * pa_x[i] + 2.0 * ts_xxxxzz_y[i] * fz_0;

        tk_xxxxzz_z[i] = -6.0 * ts_xxzz_z[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_z[i] * fe_0 + tk_xxxzz_z[i] * pa_x[i] + 2.0 * ts_xxxxzz_z[i] * fz_0;
    }

    // Set up 18-21 components of targeted buffer : IP

    auto tk_xxxyyy_x = pbuffer.data(idx_kin_ip + 18);

    auto tk_xxxyyy_y = pbuffer.data(idx_kin_ip + 19);

    auto tk_xxxyyy_z = pbuffer.data(idx_kin_ip + 20);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tk_xxxy_x,   \
                             tk_xxxyy_x,  \
                             tk_xxxyyy_x, \
                             tk_xxxyyy_y, \
                             tk_xxxyyy_z, \
                             tk_xxyyy_y,  \
                             tk_xxyyy_z,  \
                             tk_xyyy_y,   \
                             tk_xyyy_z,   \
                             ts_xxxy_x,   \
                             ts_xxxyyy_x, \
                             ts_xxxyyy_y, \
                             ts_xxxyyy_z, \
                             ts_xyyy_y,   \
                             ts_xyyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyyy_x[i] = -4.0 * ts_xxxy_x[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_x[i] * fe_0 + tk_xxxyy_x[i] * pa_y[i] + 2.0 * ts_xxxyyy_x[i] * fz_0;

        tk_xxxyyy_y[i] = -4.0 * ts_xyyy_y[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_y[i] * fe_0 + tk_xxyyy_y[i] * pa_x[i] + 2.0 * ts_xxxyyy_y[i] * fz_0;

        tk_xxxyyy_z[i] = -4.0 * ts_xyyy_z[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_z[i] * fe_0 + tk_xxyyy_z[i] * pa_x[i] + 2.0 * ts_xxxyyy_z[i] * fz_0;
    }

    // Set up 21-24 components of targeted buffer : IP

    auto tk_xxxyyz_x = pbuffer.data(idx_kin_ip + 21);

    auto tk_xxxyyz_y = pbuffer.data(idx_kin_ip + 22);

    auto tk_xxxyyz_z = pbuffer.data(idx_kin_ip + 23);

#pragma omp simd aligned(pa_z,            \
                             tk_xxxyy_0,  \
                             tk_xxxyy_x,  \
                             tk_xxxyy_y,  \
                             tk_xxxyy_z,  \
                             tk_xxxyyz_x, \
                             tk_xxxyyz_y, \
                             tk_xxxyyz_z, \
                             ts_xxxyyz_x, \
                             ts_xxxyyz_y, \
                             ts_xxxyyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyyz_x[i] = tk_xxxyy_x[i] * pa_z[i] + 2.0 * ts_xxxyyz_x[i] * fz_0;

        tk_xxxyyz_y[i] = tk_xxxyy_y[i] * pa_z[i] + 2.0 * ts_xxxyyz_y[i] * fz_0;

        tk_xxxyyz_z[i] = tk_xxxyy_0[i] * fe_0 + tk_xxxyy_z[i] * pa_z[i] + 2.0 * ts_xxxyyz_z[i] * fz_0;
    }

    // Set up 24-27 components of targeted buffer : IP

    auto tk_xxxyzz_x = pbuffer.data(idx_kin_ip + 24);

    auto tk_xxxyzz_y = pbuffer.data(idx_kin_ip + 25);

    auto tk_xxxyzz_z = pbuffer.data(idx_kin_ip + 26);

#pragma omp simd aligned(pa_y,            \
                             tk_xxxyzz_x, \
                             tk_xxxyzz_y, \
                             tk_xxxyzz_z, \
                             tk_xxxzz_0,  \
                             tk_xxxzz_x,  \
                             tk_xxxzz_y,  \
                             tk_xxxzz_z,  \
                             ts_xxxyzz_x, \
                             ts_xxxyzz_y, \
                             ts_xxxyzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyzz_x[i] = tk_xxxzz_x[i] * pa_y[i] + 2.0 * ts_xxxyzz_x[i] * fz_0;

        tk_xxxyzz_y[i] = tk_xxxzz_0[i] * fe_0 + tk_xxxzz_y[i] * pa_y[i] + 2.0 * ts_xxxyzz_y[i] * fz_0;

        tk_xxxyzz_z[i] = tk_xxxzz_z[i] * pa_y[i] + 2.0 * ts_xxxyzz_z[i] * fz_0;
    }

    // Set up 27-30 components of targeted buffer : IP

    auto tk_xxxzzz_x = pbuffer.data(idx_kin_ip + 27);

    auto tk_xxxzzz_y = pbuffer.data(idx_kin_ip + 28);

    auto tk_xxxzzz_z = pbuffer.data(idx_kin_ip + 29);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tk_xxxz_x,   \
                             tk_xxxzz_x,  \
                             tk_xxxzzz_x, \
                             tk_xxxzzz_y, \
                             tk_xxxzzz_z, \
                             tk_xxzzz_y,  \
                             tk_xxzzz_z,  \
                             tk_xzzz_y,   \
                             tk_xzzz_z,   \
                             ts_xxxz_x,   \
                             ts_xxxzzz_x, \
                             ts_xxxzzz_y, \
                             ts_xxxzzz_z, \
                             ts_xzzz_y,   \
                             ts_xzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzzz_x[i] = -4.0 * ts_xxxz_x[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_x[i] * fe_0 + tk_xxxzz_x[i] * pa_z[i] + 2.0 * ts_xxxzzz_x[i] * fz_0;

        tk_xxxzzz_y[i] = -4.0 * ts_xzzz_y[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_y[i] * fe_0 + tk_xxzzz_y[i] * pa_x[i] + 2.0 * ts_xxxzzz_y[i] * fz_0;

        tk_xxxzzz_z[i] = -4.0 * ts_xzzz_z[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_z[i] * fe_0 + tk_xxzzz_z[i] * pa_x[i] + 2.0 * ts_xxxzzz_z[i] * fz_0;
    }

    // Set up 30-33 components of targeted buffer : IP

    auto tk_xxyyyy_x = pbuffer.data(idx_kin_ip + 30);

    auto tk_xxyyyy_y = pbuffer.data(idx_kin_ip + 31);

    auto tk_xxyyyy_z = pbuffer.data(idx_kin_ip + 32);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tk_xxyy_x,   \
                             tk_xxyyy_x,  \
                             tk_xxyyyy_x, \
                             tk_xxyyyy_y, \
                             tk_xxyyyy_z, \
                             tk_xyyyy_y,  \
                             tk_xyyyy_z,  \
                             tk_yyyy_y,   \
                             tk_yyyy_z,   \
                             ts_xxyy_x,   \
                             ts_xxyyyy_x, \
                             ts_xxyyyy_y, \
                             ts_xxyyyy_z, \
                             ts_yyyy_y,   \
                             ts_yyyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyyy_x[i] = -6.0 * ts_xxyy_x[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_x[i] * fe_0 + tk_xxyyy_x[i] * pa_y[i] + 2.0 * ts_xxyyyy_x[i] * fz_0;

        tk_xxyyyy_y[i] = -2.0 * ts_yyyy_y[i] * fbe_0 * fz_0 + tk_yyyy_y[i] * fe_0 + tk_xyyyy_y[i] * pa_x[i] + 2.0 * ts_xxyyyy_y[i] * fz_0;

        tk_xxyyyy_z[i] = -2.0 * ts_yyyy_z[i] * fbe_0 * fz_0 + tk_yyyy_z[i] * fe_0 + tk_xyyyy_z[i] * pa_x[i] + 2.0 * ts_xxyyyy_z[i] * fz_0;
    }

    // Set up 33-36 components of targeted buffer : IP

    auto tk_xxyyyz_x = pbuffer.data(idx_kin_ip + 33);

    auto tk_xxyyyz_y = pbuffer.data(idx_kin_ip + 34);

    auto tk_xxyyyz_z = pbuffer.data(idx_kin_ip + 35);

#pragma omp simd aligned(pa_z,            \
                             tk_xxyyy_0,  \
                             tk_xxyyy_x,  \
                             tk_xxyyy_y,  \
                             tk_xxyyy_z,  \
                             tk_xxyyyz_x, \
                             tk_xxyyyz_y, \
                             tk_xxyyyz_z, \
                             ts_xxyyyz_x, \
                             ts_xxyyyz_y, \
                             ts_xxyyyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyyz_x[i] = tk_xxyyy_x[i] * pa_z[i] + 2.0 * ts_xxyyyz_x[i] * fz_0;

        tk_xxyyyz_y[i] = tk_xxyyy_y[i] * pa_z[i] + 2.0 * ts_xxyyyz_y[i] * fz_0;

        tk_xxyyyz_z[i] = tk_xxyyy_0[i] * fe_0 + tk_xxyyy_z[i] * pa_z[i] + 2.0 * ts_xxyyyz_z[i] * fz_0;
    }

    // Set up 36-39 components of targeted buffer : IP

    auto tk_xxyyzz_x = pbuffer.data(idx_kin_ip + 36);

    auto tk_xxyyzz_y = pbuffer.data(idx_kin_ip + 37);

    auto tk_xxyyzz_z = pbuffer.data(idx_kin_ip + 38);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tk_xxyyzz_x, \
                             tk_xxyyzz_y, \
                             tk_xxyyzz_z, \
                             tk_xxyzz_x,  \
                             tk_xxzz_x,   \
                             tk_xyyzz_y,  \
                             tk_xyyzz_z,  \
                             tk_yyzz_y,   \
                             tk_yyzz_z,   \
                             ts_xxyyzz_x, \
                             ts_xxyyzz_y, \
                             ts_xxyyzz_z, \
                             ts_xxzz_x,   \
                             ts_yyzz_y,   \
                             ts_yyzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyzz_x[i] = -2.0 * ts_xxzz_x[i] * fbe_0 * fz_0 + tk_xxzz_x[i] * fe_0 + tk_xxyzz_x[i] * pa_y[i] + 2.0 * ts_xxyyzz_x[i] * fz_0;

        tk_xxyyzz_y[i] = -2.0 * ts_yyzz_y[i] * fbe_0 * fz_0 + tk_yyzz_y[i] * fe_0 + tk_xyyzz_y[i] * pa_x[i] + 2.0 * ts_xxyyzz_y[i] * fz_0;

        tk_xxyyzz_z[i] = -2.0 * ts_yyzz_z[i] * fbe_0 * fz_0 + tk_yyzz_z[i] * fe_0 + tk_xyyzz_z[i] * pa_x[i] + 2.0 * ts_xxyyzz_z[i] * fz_0;
    }

    // Set up 39-42 components of targeted buffer : IP

    auto tk_xxyzzz_x = pbuffer.data(idx_kin_ip + 39);

    auto tk_xxyzzz_y = pbuffer.data(idx_kin_ip + 40);

    auto tk_xxyzzz_z = pbuffer.data(idx_kin_ip + 41);

#pragma omp simd aligned(pa_y,            \
                             tk_xxyzzz_x, \
                             tk_xxyzzz_y, \
                             tk_xxyzzz_z, \
                             tk_xxzzz_0,  \
                             tk_xxzzz_x,  \
                             tk_xxzzz_y,  \
                             tk_xxzzz_z,  \
                             ts_xxyzzz_x, \
                             ts_xxyzzz_y, \
                             ts_xxyzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzzz_x[i] = tk_xxzzz_x[i] * pa_y[i] + 2.0 * ts_xxyzzz_x[i] * fz_0;

        tk_xxyzzz_y[i] = tk_xxzzz_0[i] * fe_0 + tk_xxzzz_y[i] * pa_y[i] + 2.0 * ts_xxyzzz_y[i] * fz_0;

        tk_xxyzzz_z[i] = tk_xxzzz_z[i] * pa_y[i] + 2.0 * ts_xxyzzz_z[i] * fz_0;
    }

    // Set up 42-45 components of targeted buffer : IP

    auto tk_xxzzzz_x = pbuffer.data(idx_kin_ip + 42);

    auto tk_xxzzzz_y = pbuffer.data(idx_kin_ip + 43);

    auto tk_xxzzzz_z = pbuffer.data(idx_kin_ip + 44);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tk_xxzz_x,   \
                             tk_xxzzz_x,  \
                             tk_xxzzzz_x, \
                             tk_xxzzzz_y, \
                             tk_xxzzzz_z, \
                             tk_xzzzz_y,  \
                             tk_xzzzz_z,  \
                             tk_zzzz_y,   \
                             tk_zzzz_z,   \
                             ts_xxzz_x,   \
                             ts_xxzzzz_x, \
                             ts_xxzzzz_y, \
                             ts_xxzzzz_z, \
                             ts_zzzz_y,   \
                             ts_zzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzzz_x[i] = -6.0 * ts_xxzz_x[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_x[i] * fe_0 + tk_xxzzz_x[i] * pa_z[i] + 2.0 * ts_xxzzzz_x[i] * fz_0;

        tk_xxzzzz_y[i] = -2.0 * ts_zzzz_y[i] * fbe_0 * fz_0 + tk_zzzz_y[i] * fe_0 + tk_xzzzz_y[i] * pa_x[i] + 2.0 * ts_xxzzzz_y[i] * fz_0;

        tk_xxzzzz_z[i] = -2.0 * ts_zzzz_z[i] * fbe_0 * fz_0 + tk_zzzz_z[i] * fe_0 + tk_xzzzz_z[i] * pa_x[i] + 2.0 * ts_xxzzzz_z[i] * fz_0;
    }

    // Set up 45-48 components of targeted buffer : IP

    auto tk_xyyyyy_x = pbuffer.data(idx_kin_ip + 45);

    auto tk_xyyyyy_y = pbuffer.data(idx_kin_ip + 46);

    auto tk_xyyyyy_z = pbuffer.data(idx_kin_ip + 47);

#pragma omp simd aligned(pa_x,            \
                             tk_xyyyyy_x, \
                             tk_xyyyyy_y, \
                             tk_xyyyyy_z, \
                             tk_yyyyy_0,  \
                             tk_yyyyy_x,  \
                             tk_yyyyy_y,  \
                             tk_yyyyy_z,  \
                             ts_xyyyyy_x, \
                             ts_xyyyyy_y, \
                             ts_xyyyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyy_x[i] = tk_yyyyy_0[i] * fe_0 + tk_yyyyy_x[i] * pa_x[i] + 2.0 * ts_xyyyyy_x[i] * fz_0;

        tk_xyyyyy_y[i] = tk_yyyyy_y[i] * pa_x[i] + 2.0 * ts_xyyyyy_y[i] * fz_0;

        tk_xyyyyy_z[i] = tk_yyyyy_z[i] * pa_x[i] + 2.0 * ts_xyyyyy_z[i] * fz_0;
    }

    // Set up 48-51 components of targeted buffer : IP

    auto tk_xyyyyz_x = pbuffer.data(idx_kin_ip + 48);

    auto tk_xyyyyz_y = pbuffer.data(idx_kin_ip + 49);

    auto tk_xyyyyz_z = pbuffer.data(idx_kin_ip + 50);

#pragma omp simd aligned( \
        pa_x, pa_z, tk_xyyyy_x, tk_xyyyyz_x, tk_xyyyyz_y, tk_xyyyyz_z, tk_yyyyz_y, tk_yyyyz_z, ts_xyyyyz_x, ts_xyyyyz_y, ts_xyyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyyyyz_x[i] = tk_xyyyy_x[i] * pa_z[i] + 2.0 * ts_xyyyyz_x[i] * fz_0;

        tk_xyyyyz_y[i] = tk_yyyyz_y[i] * pa_x[i] + 2.0 * ts_xyyyyz_y[i] * fz_0;

        tk_xyyyyz_z[i] = tk_yyyyz_z[i] * pa_x[i] + 2.0 * ts_xyyyyz_z[i] * fz_0;
    }

    // Set up 51-54 components of targeted buffer : IP

    auto tk_xyyyzz_x = pbuffer.data(idx_kin_ip + 51);

    auto tk_xyyyzz_y = pbuffer.data(idx_kin_ip + 52);

    auto tk_xyyyzz_z = pbuffer.data(idx_kin_ip + 53);

#pragma omp simd aligned(pa_x,            \
                             tk_xyyyzz_x, \
                             tk_xyyyzz_y, \
                             tk_xyyyzz_z, \
                             tk_yyyzz_0,  \
                             tk_yyyzz_x,  \
                             tk_yyyzz_y,  \
                             tk_yyyzz_z,  \
                             ts_xyyyzz_x, \
                             ts_xyyyzz_y, \
                             ts_xyyyzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyzz_x[i] = tk_yyyzz_0[i] * fe_0 + tk_yyyzz_x[i] * pa_x[i] + 2.0 * ts_xyyyzz_x[i] * fz_0;

        tk_xyyyzz_y[i] = tk_yyyzz_y[i] * pa_x[i] + 2.0 * ts_xyyyzz_y[i] * fz_0;

        tk_xyyyzz_z[i] = tk_yyyzz_z[i] * pa_x[i] + 2.0 * ts_xyyyzz_z[i] * fz_0;
    }

    // Set up 54-57 components of targeted buffer : IP

    auto tk_xyyzzz_x = pbuffer.data(idx_kin_ip + 54);

    auto tk_xyyzzz_y = pbuffer.data(idx_kin_ip + 55);

    auto tk_xyyzzz_z = pbuffer.data(idx_kin_ip + 56);

#pragma omp simd aligned(pa_x,            \
                             tk_xyyzzz_x, \
                             tk_xyyzzz_y, \
                             tk_xyyzzz_z, \
                             tk_yyzzz_0,  \
                             tk_yyzzz_x,  \
                             tk_yyzzz_y,  \
                             tk_yyzzz_z,  \
                             ts_xyyzzz_x, \
                             ts_xyyzzz_y, \
                             ts_xyyzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzzz_x[i] = tk_yyzzz_0[i] * fe_0 + tk_yyzzz_x[i] * pa_x[i] + 2.0 * ts_xyyzzz_x[i] * fz_0;

        tk_xyyzzz_y[i] = tk_yyzzz_y[i] * pa_x[i] + 2.0 * ts_xyyzzz_y[i] * fz_0;

        tk_xyyzzz_z[i] = tk_yyzzz_z[i] * pa_x[i] + 2.0 * ts_xyyzzz_z[i] * fz_0;
    }

    // Set up 57-60 components of targeted buffer : IP

    auto tk_xyzzzz_x = pbuffer.data(idx_kin_ip + 57);

    auto tk_xyzzzz_y = pbuffer.data(idx_kin_ip + 58);

    auto tk_xyzzzz_z = pbuffer.data(idx_kin_ip + 59);

#pragma omp simd aligned( \
        pa_x, pa_y, tk_xyzzzz_x, tk_xyzzzz_y, tk_xyzzzz_z, tk_xzzzz_x, tk_yzzzz_y, tk_yzzzz_z, ts_xyzzzz_x, ts_xyzzzz_y, ts_xyzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyzzzz_x[i] = tk_xzzzz_x[i] * pa_y[i] + 2.0 * ts_xyzzzz_x[i] * fz_0;

        tk_xyzzzz_y[i] = tk_yzzzz_y[i] * pa_x[i] + 2.0 * ts_xyzzzz_y[i] * fz_0;

        tk_xyzzzz_z[i] = tk_yzzzz_z[i] * pa_x[i] + 2.0 * ts_xyzzzz_z[i] * fz_0;
    }

    // Set up 60-63 components of targeted buffer : IP

    auto tk_xzzzzz_x = pbuffer.data(idx_kin_ip + 60);

    auto tk_xzzzzz_y = pbuffer.data(idx_kin_ip + 61);

    auto tk_xzzzzz_z = pbuffer.data(idx_kin_ip + 62);

#pragma omp simd aligned(pa_x,            \
                             tk_xzzzzz_x, \
                             tk_xzzzzz_y, \
                             tk_xzzzzz_z, \
                             tk_zzzzz_0,  \
                             tk_zzzzz_x,  \
                             tk_zzzzz_y,  \
                             tk_zzzzz_z,  \
                             ts_xzzzzz_x, \
                             ts_xzzzzz_y, \
                             ts_xzzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzzz_x[i] = tk_zzzzz_0[i] * fe_0 + tk_zzzzz_x[i] * pa_x[i] + 2.0 * ts_xzzzzz_x[i] * fz_0;

        tk_xzzzzz_y[i] = tk_zzzzz_y[i] * pa_x[i] + 2.0 * ts_xzzzzz_y[i] * fz_0;

        tk_xzzzzz_z[i] = tk_zzzzz_z[i] * pa_x[i] + 2.0 * ts_xzzzzz_z[i] * fz_0;
    }

    // Set up 63-66 components of targeted buffer : IP

    auto tk_yyyyyy_x = pbuffer.data(idx_kin_ip + 63);

    auto tk_yyyyyy_y = pbuffer.data(idx_kin_ip + 64);

    auto tk_yyyyyy_z = pbuffer.data(idx_kin_ip + 65);

#pragma omp simd aligned(pa_y,            \
                             tk_yyyy_x,   \
                             tk_yyyy_y,   \
                             tk_yyyy_z,   \
                             tk_yyyyy_0,  \
                             tk_yyyyy_x,  \
                             tk_yyyyy_y,  \
                             tk_yyyyy_z,  \
                             tk_yyyyyy_x, \
                             tk_yyyyyy_y, \
                             tk_yyyyyy_z, \
                             ts_yyyy_x,   \
                             ts_yyyy_y,   \
                             ts_yyyy_z,   \
                             ts_yyyyyy_x, \
                             ts_yyyyyy_y, \
                             ts_yyyyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyyy_x[i] = -10.0 * ts_yyyy_x[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_x[i] * fe_0 + tk_yyyyy_x[i] * pa_y[i] + 2.0 * ts_yyyyyy_x[i] * fz_0;

        tk_yyyyyy_y[i] = -10.0 * ts_yyyy_y[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_y[i] * fe_0 + tk_yyyyy_0[i] * fe_0 + tk_yyyyy_y[i] * pa_y[i] +
                         2.0 * ts_yyyyyy_y[i] * fz_0;

        tk_yyyyyy_z[i] = -10.0 * ts_yyyy_z[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_z[i] * fe_0 + tk_yyyyy_z[i] * pa_y[i] + 2.0 * ts_yyyyyy_z[i] * fz_0;
    }

    // Set up 66-69 components of targeted buffer : IP

    auto tk_yyyyyz_x = pbuffer.data(idx_kin_ip + 66);

    auto tk_yyyyyz_y = pbuffer.data(idx_kin_ip + 67);

    auto tk_yyyyyz_z = pbuffer.data(idx_kin_ip + 68);

#pragma omp simd aligned(pa_z,            \
                             tk_yyyyy_0,  \
                             tk_yyyyy_x,  \
                             tk_yyyyy_y,  \
                             tk_yyyyy_z,  \
                             tk_yyyyyz_x, \
                             tk_yyyyyz_y, \
                             tk_yyyyyz_z, \
                             ts_yyyyyz_x, \
                             ts_yyyyyz_y, \
                             ts_yyyyyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyyz_x[i] = tk_yyyyy_x[i] * pa_z[i] + 2.0 * ts_yyyyyz_x[i] * fz_0;

        tk_yyyyyz_y[i] = tk_yyyyy_y[i] * pa_z[i] + 2.0 * ts_yyyyyz_y[i] * fz_0;

        tk_yyyyyz_z[i] = tk_yyyyy_0[i] * fe_0 + tk_yyyyy_z[i] * pa_z[i] + 2.0 * ts_yyyyyz_z[i] * fz_0;
    }

    // Set up 69-72 components of targeted buffer : IP

    auto tk_yyyyzz_x = pbuffer.data(idx_kin_ip + 69);

    auto tk_yyyyzz_y = pbuffer.data(idx_kin_ip + 70);

    auto tk_yyyyzz_z = pbuffer.data(idx_kin_ip + 71);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tk_yyyy_y,   \
                             tk_yyyyz_y,  \
                             tk_yyyyzz_x, \
                             tk_yyyyzz_y, \
                             tk_yyyyzz_z, \
                             tk_yyyzz_x,  \
                             tk_yyyzz_z,  \
                             tk_yyzz_x,   \
                             tk_yyzz_z,   \
                             ts_yyyy_y,   \
                             ts_yyyyzz_x, \
                             ts_yyyyzz_y, \
                             ts_yyyyzz_z, \
                             ts_yyzz_x,   \
                             ts_yyzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyzz_x[i] = -6.0 * ts_yyzz_x[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_x[i] * fe_0 + tk_yyyzz_x[i] * pa_y[i] + 2.0 * ts_yyyyzz_x[i] * fz_0;

        tk_yyyyzz_y[i] = -2.0 * ts_yyyy_y[i] * fbe_0 * fz_0 + tk_yyyy_y[i] * fe_0 + tk_yyyyz_y[i] * pa_z[i] + 2.0 * ts_yyyyzz_y[i] * fz_0;

        tk_yyyyzz_z[i] = -6.0 * ts_yyzz_z[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_z[i] * fe_0 + tk_yyyzz_z[i] * pa_y[i] + 2.0 * ts_yyyyzz_z[i] * fz_0;
    }

    // Set up 72-75 components of targeted buffer : IP

    auto tk_yyyzzz_x = pbuffer.data(idx_kin_ip + 72);

    auto tk_yyyzzz_y = pbuffer.data(idx_kin_ip + 73);

    auto tk_yyyzzz_z = pbuffer.data(idx_kin_ip + 74);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tk_yyyz_y,   \
                             tk_yyyzz_y,  \
                             tk_yyyzzz_x, \
                             tk_yyyzzz_y, \
                             tk_yyyzzz_z, \
                             tk_yyzzz_x,  \
                             tk_yyzzz_z,  \
                             tk_yzzz_x,   \
                             tk_yzzz_z,   \
                             ts_yyyz_y,   \
                             ts_yyyzzz_x, \
                             ts_yyyzzz_y, \
                             ts_yyyzzz_z, \
                             ts_yzzz_x,   \
                             ts_yzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzzz_x[i] = -4.0 * ts_yzzz_x[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_x[i] * fe_0 + tk_yyzzz_x[i] * pa_y[i] + 2.0 * ts_yyyzzz_x[i] * fz_0;

        tk_yyyzzz_y[i] = -4.0 * ts_yyyz_y[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_y[i] * fe_0 + tk_yyyzz_y[i] * pa_z[i] + 2.0 * ts_yyyzzz_y[i] * fz_0;

        tk_yyyzzz_z[i] = -4.0 * ts_yzzz_z[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_z[i] * fe_0 + tk_yyzzz_z[i] * pa_y[i] + 2.0 * ts_yyyzzz_z[i] * fz_0;
    }

    // Set up 75-78 components of targeted buffer : IP

    auto tk_yyzzzz_x = pbuffer.data(idx_kin_ip + 75);

    auto tk_yyzzzz_y = pbuffer.data(idx_kin_ip + 76);

    auto tk_yyzzzz_z = pbuffer.data(idx_kin_ip + 77);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tk_yyzz_y,   \
                             tk_yyzzz_y,  \
                             tk_yyzzzz_x, \
                             tk_yyzzzz_y, \
                             tk_yyzzzz_z, \
                             tk_yzzzz_x,  \
                             tk_yzzzz_z,  \
                             tk_zzzz_x,   \
                             tk_zzzz_z,   \
                             ts_yyzz_y,   \
                             ts_yyzzzz_x, \
                             ts_yyzzzz_y, \
                             ts_yyzzzz_z, \
                             ts_zzzz_x,   \
                             ts_zzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzzz_x[i] = -2.0 * ts_zzzz_x[i] * fbe_0 * fz_0 + tk_zzzz_x[i] * fe_0 + tk_yzzzz_x[i] * pa_y[i] + 2.0 * ts_yyzzzz_x[i] * fz_0;

        tk_yyzzzz_y[i] = -6.0 * ts_yyzz_y[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_y[i] * fe_0 + tk_yyzzz_y[i] * pa_z[i] + 2.0 * ts_yyzzzz_y[i] * fz_0;

        tk_yyzzzz_z[i] = -2.0 * ts_zzzz_z[i] * fbe_0 * fz_0 + tk_zzzz_z[i] * fe_0 + tk_yzzzz_z[i] * pa_y[i] + 2.0 * ts_yyzzzz_z[i] * fz_0;
    }

    // Set up 78-81 components of targeted buffer : IP

    auto tk_yzzzzz_x = pbuffer.data(idx_kin_ip + 78);

    auto tk_yzzzzz_y = pbuffer.data(idx_kin_ip + 79);

    auto tk_yzzzzz_z = pbuffer.data(idx_kin_ip + 80);

#pragma omp simd aligned(pa_y,            \
                             tk_yzzzzz_x, \
                             tk_yzzzzz_y, \
                             tk_yzzzzz_z, \
                             tk_zzzzz_0,  \
                             tk_zzzzz_x,  \
                             tk_zzzzz_y,  \
                             tk_zzzzz_z,  \
                             ts_yzzzzz_x, \
                             ts_yzzzzz_y, \
                             ts_yzzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzzz_x[i] = tk_zzzzz_x[i] * pa_y[i] + 2.0 * ts_yzzzzz_x[i] * fz_0;

        tk_yzzzzz_y[i] = tk_zzzzz_0[i] * fe_0 + tk_zzzzz_y[i] * pa_y[i] + 2.0 * ts_yzzzzz_y[i] * fz_0;

        tk_yzzzzz_z[i] = tk_zzzzz_z[i] * pa_y[i] + 2.0 * ts_yzzzzz_z[i] * fz_0;
    }

    // Set up 81-84 components of targeted buffer : IP

    auto tk_zzzzzz_x = pbuffer.data(idx_kin_ip + 81);

    auto tk_zzzzzz_y = pbuffer.data(idx_kin_ip + 82);

    auto tk_zzzzzz_z = pbuffer.data(idx_kin_ip + 83);

#pragma omp simd aligned(pa_z,            \
                             tk_zzzz_x,   \
                             tk_zzzz_y,   \
                             tk_zzzz_z,   \
                             tk_zzzzz_0,  \
                             tk_zzzzz_x,  \
                             tk_zzzzz_y,  \
                             tk_zzzzz_z,  \
                             tk_zzzzzz_x, \
                             tk_zzzzzz_y, \
                             tk_zzzzzz_z, \
                             ts_zzzz_x,   \
                             ts_zzzz_y,   \
                             ts_zzzz_z,   \
                             ts_zzzzzz_x, \
                             ts_zzzzzz_y, \
                             ts_zzzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzzz_x[i] = -10.0 * ts_zzzz_x[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_x[i] * fe_0 + tk_zzzzz_x[i] * pa_z[i] + 2.0 * ts_zzzzzz_x[i] * fz_0;

        tk_zzzzzz_y[i] = -10.0 * ts_zzzz_y[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_y[i] * fe_0 + tk_zzzzz_y[i] * pa_z[i] + 2.0 * ts_zzzzzz_y[i] * fz_0;

        tk_zzzzzz_z[i] = -10.0 * ts_zzzz_z[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_z[i] * fe_0 + tk_zzzzz_0[i] * fe_0 + tk_zzzzz_z[i] * pa_z[i] +
                         2.0 * ts_zzzzzz_z[i] * fz_0;
    }
}

}  // namespace kinrec
