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

#include "KineticEnergyPrimRecSH.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_sh(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_sh,
                            const size_t              idx_ovl_sf,
                            const size_t              idx_kin_sf,
                            const size_t              idx_kin_sg,
                            const size_t              idx_ovl_sh,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpb,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_ovl_sf);

    auto ts_0_xyy = pbuffer.data(idx_ovl_sf + 3);

    auto ts_0_xzz = pbuffer.data(idx_ovl_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_ovl_sf + 6);

    auto ts_0_yzz = pbuffer.data(idx_ovl_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_ovl_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto tk_0_xxx = pbuffer.data(idx_kin_sf);

    auto tk_0_xyy = pbuffer.data(idx_kin_sf + 3);

    auto tk_0_xzz = pbuffer.data(idx_kin_sf + 5);

    auto tk_0_yyy = pbuffer.data(idx_kin_sf + 6);

    auto tk_0_yzz = pbuffer.data(idx_kin_sf + 8);

    auto tk_0_zzz = pbuffer.data(idx_kin_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto tk_0_xxxx = pbuffer.data(idx_kin_sg);

    auto tk_0_xxxz = pbuffer.data(idx_kin_sg + 2);

    auto tk_0_xxyy = pbuffer.data(idx_kin_sg + 3);

    auto tk_0_xxzz = pbuffer.data(idx_kin_sg + 5);

    auto tk_0_xyyy = pbuffer.data(idx_kin_sg + 6);

    auto tk_0_xzzz = pbuffer.data(idx_kin_sg + 9);

    auto tk_0_yyyy = pbuffer.data(idx_kin_sg + 10);

    auto tk_0_yyyz = pbuffer.data(idx_kin_sg + 11);

    auto tk_0_yyzz = pbuffer.data(idx_kin_sg + 12);

    auto tk_0_yzzz = pbuffer.data(idx_kin_sg + 13);

    auto tk_0_zzzz = pbuffer.data(idx_kin_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_ovl_sh);

    auto ts_0_xxxxy = pbuffer.data(idx_ovl_sh + 1);

    auto ts_0_xxxxz = pbuffer.data(idx_ovl_sh + 2);

    auto ts_0_xxxyy = pbuffer.data(idx_ovl_sh + 3);

    auto ts_0_xxxyz = pbuffer.data(idx_ovl_sh + 4);

    auto ts_0_xxxzz = pbuffer.data(idx_ovl_sh + 5);

    auto ts_0_xxyyy = pbuffer.data(idx_ovl_sh + 6);

    auto ts_0_xxyyz = pbuffer.data(idx_ovl_sh + 7);

    auto ts_0_xxyzz = pbuffer.data(idx_ovl_sh + 8);

    auto ts_0_xxzzz = pbuffer.data(idx_ovl_sh + 9);

    auto ts_0_xyyyy = pbuffer.data(idx_ovl_sh + 10);

    auto ts_0_xyyyz = pbuffer.data(idx_ovl_sh + 11);

    auto ts_0_xyyzz = pbuffer.data(idx_ovl_sh + 12);

    auto ts_0_xyzzz = pbuffer.data(idx_ovl_sh + 13);

    auto ts_0_xzzzz = pbuffer.data(idx_ovl_sh + 14);

    auto ts_0_yyyyy = pbuffer.data(idx_ovl_sh + 15);

    auto ts_0_yyyyz = pbuffer.data(idx_ovl_sh + 16);

    auto ts_0_yyyzz = pbuffer.data(idx_ovl_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_ovl_sh + 18);

    auto ts_0_yzzzz = pbuffer.data(idx_ovl_sh + 19);

    auto ts_0_zzzzz = pbuffer.data(idx_ovl_sh + 20);

    // Set up components of targeted buffer : SH

    auto tk_0_xxxxx = pbuffer.data(idx_kin_sh);

    auto tk_0_xxxxy = pbuffer.data(idx_kin_sh + 1);

    auto tk_0_xxxxz = pbuffer.data(idx_kin_sh + 2);

    auto tk_0_xxxyy = pbuffer.data(idx_kin_sh + 3);

    auto tk_0_xxxyz = pbuffer.data(idx_kin_sh + 4);

    auto tk_0_xxxzz = pbuffer.data(idx_kin_sh + 5);

    auto tk_0_xxyyy = pbuffer.data(idx_kin_sh + 6);

    auto tk_0_xxyyz = pbuffer.data(idx_kin_sh + 7);

    auto tk_0_xxyzz = pbuffer.data(idx_kin_sh + 8);

    auto tk_0_xxzzz = pbuffer.data(idx_kin_sh + 9);

    auto tk_0_xyyyy = pbuffer.data(idx_kin_sh + 10);

    auto tk_0_xyyyz = pbuffer.data(idx_kin_sh + 11);

    auto tk_0_xyyzz = pbuffer.data(idx_kin_sh + 12);

    auto tk_0_xyzzz = pbuffer.data(idx_kin_sh + 13);

    auto tk_0_xzzzz = pbuffer.data(idx_kin_sh + 14);

    auto tk_0_yyyyy = pbuffer.data(idx_kin_sh + 15);

    auto tk_0_yyyyz = pbuffer.data(idx_kin_sh + 16);

    auto tk_0_yyyzz = pbuffer.data(idx_kin_sh + 17);

    auto tk_0_yyzzz = pbuffer.data(idx_kin_sh + 18);

    auto tk_0_yzzzz = pbuffer.data(idx_kin_sh + 19);

    auto tk_0_zzzzz = pbuffer.data(idx_kin_sh + 20);

#pragma omp simd aligned(pb_x,           \
                             pb_y,       \
                             pb_z,       \
                             tk_0_xxx,   \
                             tk_0_xxxx,  \
                             tk_0_xxxxx, \
                             tk_0_xxxxy, \
                             tk_0_xxxxz, \
                             tk_0_xxxyy, \
                             tk_0_xxxyz, \
                             tk_0_xxxz,  \
                             tk_0_xxxzz, \
                             tk_0_xxyy,  \
                             tk_0_xxyyy, \
                             tk_0_xxyyz, \
                             tk_0_xxyzz, \
                             tk_0_xxzz,  \
                             tk_0_xxzzz, \
                             tk_0_xyy,   \
                             tk_0_xyyy,  \
                             tk_0_xyyyy, \
                             tk_0_xyyyz, \
                             tk_0_xyyzz, \
                             tk_0_xyzzz, \
                             tk_0_xzz,   \
                             tk_0_xzzz,  \
                             tk_0_xzzzz, \
                             tk_0_yyy,   \
                             tk_0_yyyy,  \
                             tk_0_yyyyy, \
                             tk_0_yyyyz, \
                             tk_0_yyyz,  \
                             tk_0_yyyzz, \
                             tk_0_yyzz,  \
                             tk_0_yyzzz, \
                             tk_0_yzz,   \
                             tk_0_yzzz,  \
                             tk_0_yzzzz, \
                             tk_0_zzz,   \
                             tk_0_zzzz,  \
                             tk_0_zzzzz, \
                             ts_0_xxx,   \
                             ts_0_xxxxx, \
                             ts_0_xxxxy, \
                             ts_0_xxxxz, \
                             ts_0_xxxyy, \
                             ts_0_xxxyz, \
                             ts_0_xxxzz, \
                             ts_0_xxyyy, \
                             ts_0_xxyyz, \
                             ts_0_xxyzz, \
                             ts_0_xxzzz, \
                             ts_0_xyy,   \
                             ts_0_xyyyy, \
                             ts_0_xyyyz, \
                             ts_0_xyyzz, \
                             ts_0_xyzzz, \
                             ts_0_xzz,   \
                             ts_0_xzzzz, \
                             ts_0_yyy,   \
                             ts_0_yyyyy, \
                             ts_0_yyyyz, \
                             ts_0_yyyzz, \
                             ts_0_yyzzz, \
                             ts_0_yzz,   \
                             ts_0_yzzzz, \
                             ts_0_zzz,   \
                             ts_0_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fke_0 = 0.5 / b_exps[i];

        tk_0_xxxxx[i] = -8.0 * ts_0_xxx[i] * fke_0 * fz_0 + 4.0 * tk_0_xxx[i] * fe_0 + tk_0_xxxx[i] * pb_x[i] + 2.0 * ts_0_xxxxx[i] * fz_0;

        tk_0_xxxxy[i] = tk_0_xxxx[i] * pb_y[i] + 2.0 * ts_0_xxxxy[i] * fz_0;

        tk_0_xxxxz[i] = tk_0_xxxx[i] * pb_z[i] + 2.0 * ts_0_xxxxz[i] * fz_0;

        tk_0_xxxyy[i] = -4.0 * ts_0_xyy[i] * fke_0 * fz_0 + 2.0 * tk_0_xyy[i] * fe_0 + tk_0_xxyy[i] * pb_x[i] + 2.0 * ts_0_xxxyy[i] * fz_0;

        tk_0_xxxyz[i] = tk_0_xxxz[i] * pb_y[i] + 2.0 * ts_0_xxxyz[i] * fz_0;

        tk_0_xxxzz[i] = -4.0 * ts_0_xzz[i] * fke_0 * fz_0 + 2.0 * tk_0_xzz[i] * fe_0 + tk_0_xxzz[i] * pb_x[i] + 2.0 * ts_0_xxxzz[i] * fz_0;

        tk_0_xxyyy[i] = -2.0 * ts_0_yyy[i] * fke_0 * fz_0 + tk_0_yyy[i] * fe_0 + tk_0_xyyy[i] * pb_x[i] + 2.0 * ts_0_xxyyy[i] * fz_0;

        tk_0_xxyyz[i] = tk_0_xxyy[i] * pb_z[i] + 2.0 * ts_0_xxyyz[i] * fz_0;

        tk_0_xxyzz[i] = tk_0_xxzz[i] * pb_y[i] + 2.0 * ts_0_xxyzz[i] * fz_0;

        tk_0_xxzzz[i] = -2.0 * ts_0_zzz[i] * fke_0 * fz_0 + tk_0_zzz[i] * fe_0 + tk_0_xzzz[i] * pb_x[i] + 2.0 * ts_0_xxzzz[i] * fz_0;

        tk_0_xyyyy[i] = tk_0_yyyy[i] * pb_x[i] + 2.0 * ts_0_xyyyy[i] * fz_0;

        tk_0_xyyyz[i] = tk_0_yyyz[i] * pb_x[i] + 2.0 * ts_0_xyyyz[i] * fz_0;

        tk_0_xyyzz[i] = tk_0_yyzz[i] * pb_x[i] + 2.0 * ts_0_xyyzz[i] * fz_0;

        tk_0_xyzzz[i] = tk_0_yzzz[i] * pb_x[i] + 2.0 * ts_0_xyzzz[i] * fz_0;

        tk_0_xzzzz[i] = tk_0_zzzz[i] * pb_x[i] + 2.0 * ts_0_xzzzz[i] * fz_0;

        tk_0_yyyyy[i] = -8.0 * ts_0_yyy[i] * fke_0 * fz_0 + 4.0 * tk_0_yyy[i] * fe_0 + tk_0_yyyy[i] * pb_y[i] + 2.0 * ts_0_yyyyy[i] * fz_0;

        tk_0_yyyyz[i] = tk_0_yyyy[i] * pb_z[i] + 2.0 * ts_0_yyyyz[i] * fz_0;

        tk_0_yyyzz[i] = -4.0 * ts_0_yzz[i] * fke_0 * fz_0 + 2.0 * tk_0_yzz[i] * fe_0 + tk_0_yyzz[i] * pb_y[i] + 2.0 * ts_0_yyyzz[i] * fz_0;

        tk_0_yyzzz[i] = -2.0 * ts_0_zzz[i] * fke_0 * fz_0 + tk_0_zzz[i] * fe_0 + tk_0_yzzz[i] * pb_y[i] + 2.0 * ts_0_yyzzz[i] * fz_0;

        tk_0_yzzzz[i] = tk_0_zzzz[i] * pb_y[i] + 2.0 * ts_0_yzzzz[i] * fz_0;

        tk_0_zzzzz[i] = -8.0 * ts_0_zzz[i] * fke_0 * fz_0 + 4.0 * tk_0_zzz[i] * fe_0 + tk_0_zzzz[i] * pb_z[i] + 2.0 * ts_0_zzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
