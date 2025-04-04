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

#include "KineticEnergyPrimRecSD.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_sd(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_sd,
                            const size_t              idx_ovl_ss,
                            const size_t              idx_kin_ss,
                            const size_t              idx_kin_sp,
                            const size_t              idx_ovl_sd,
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

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of auxiliary buffer : SS

    auto tk_0_0 = pbuffer.data(idx_kin_ss);

    // Set up components of auxiliary buffer : SP

    auto tk_0_x = pbuffer.data(idx_kin_sp);

    auto tk_0_y = pbuffer.data(idx_kin_sp + 1);

    auto tk_0_z = pbuffer.data(idx_kin_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_xy = pbuffer.data(idx_ovl_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_ovl_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up components of targeted buffer : SD

    auto tk_0_xx = pbuffer.data(idx_kin_sd);

    auto tk_0_xy = pbuffer.data(idx_kin_sd + 1);

    auto tk_0_xz = pbuffer.data(idx_kin_sd + 2);

    auto tk_0_yy = pbuffer.data(idx_kin_sd + 3);

    auto tk_0_yz = pbuffer.data(idx_kin_sd + 4);

    auto tk_0_zz = pbuffer.data(idx_kin_sd + 5);

#pragma omp simd aligned(pb_x,        \
                             pb_y,    \
                             pb_z,    \
                             tk_0_0,  \
                             tk_0_x,  \
                             tk_0_xx, \
                             tk_0_xy, \
                             tk_0_xz, \
                             tk_0_y,  \
                             tk_0_yy, \
                             tk_0_yz, \
                             tk_0_z,  \
                             tk_0_zz, \
                             ts_0_0,  \
                             ts_0_xx, \
                             ts_0_xy, \
                             ts_0_xz, \
                             ts_0_yy, \
                             ts_0_yz, \
                             ts_0_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fke_0 = 0.5 / b_exps[i];

        tk_0_xx[i] = -2.0 * ts_0_0[i] * fke_0 * fz_0 + tk_0_0[i] * fe_0 + tk_0_x[i] * pb_x[i] + 2.0 * ts_0_xx[i] * fz_0;

        tk_0_xy[i] = tk_0_y[i] * pb_x[i] + 2.0 * ts_0_xy[i] * fz_0;

        tk_0_xz[i] = tk_0_z[i] * pb_x[i] + 2.0 * ts_0_xz[i] * fz_0;

        tk_0_yy[i] = -2.0 * ts_0_0[i] * fke_0 * fz_0 + tk_0_0[i] * fe_0 + tk_0_y[i] * pb_y[i] + 2.0 * ts_0_yy[i] * fz_0;

        tk_0_yz[i] = tk_0_z[i] * pb_y[i] + 2.0 * ts_0_yz[i] * fz_0;

        tk_0_zz[i] = -2.0 * ts_0_0[i] * fke_0 * fz_0 + tk_0_0[i] * fe_0 + tk_0_z[i] * pb_z[i] + 2.0 * ts_0_zz[i] * fz_0;
    }
}

}  // namespace kinrec
