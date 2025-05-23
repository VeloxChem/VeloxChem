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

#include "OverlapPrimRecSD.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_sd(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_sd,
                     const size_t              idx_ovl_ss,
                     const size_t              idx_ovl_sp,
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

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of targeted buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_xy = pbuffer.data(idx_ovl_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_ovl_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

#pragma omp simd aligned(pb_x, pb_y, pb_z, ts_0_0, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_0_xx[i] = ts_0_0[i] * fe_0 + ts_0_x[i] * pb_x[i];

        ts_0_xy[i] = ts_0_y[i] * pb_x[i];

        ts_0_xz[i] = ts_0_z[i] * pb_x[i];

        ts_0_yy[i] = ts_0_0[i] * fe_0 + ts_0_y[i] * pb_y[i];

        ts_0_yz[i] = ts_0_z[i] * pb_y[i];

        ts_0_zz[i] = ts_0_0[i] * fe_0 + ts_0_z[i] * pb_z[i];
    }
}

}  // namespace ovlrec
