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

#include "OverlapPrimRecDS.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_ds(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_ds,
                     const size_t              idx_ovl_ss,
                     const size_t              idx_ovl_ps,
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

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ovl_ps);

    auto ts_y_0 = pbuffer.data(idx_ovl_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ovl_ps + 2);

    // Set up components of targeted buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ovl_ds);

    auto ts_xy_0 = pbuffer.data(idx_ovl_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ovl_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ovl_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ovl_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ovl_ds + 5);

#pragma omp simd aligned(pa_x, pa_y, pa_z, ts_0_0, ts_x_0, ts_xx_0, ts_xy_0, ts_xz_0, ts_y_0, ts_yy_0, ts_yz_0, ts_z_0, ts_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xx_0[i] = ts_0_0[i] * fe_0 + ts_x_0[i] * pa_x[i];

        ts_xy_0[i] = ts_y_0[i] * pa_x[i];

        ts_xz_0[i] = ts_z_0[i] * pa_x[i];

        ts_yy_0[i] = ts_0_0[i] * fe_0 + ts_y_0[i] * pa_y[i];

        ts_yz_0[i] = ts_z_0[i] * pa_y[i];

        ts_zz_0[i] = ts_0_0[i] * fe_0 + ts_z_0[i] * pa_z[i];
    }
}

}  // namespace ovlrec
