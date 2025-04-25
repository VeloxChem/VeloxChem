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

#include "OverlapPrimRecPD.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_pd(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_pd,
                     const size_t              idx_ovl_sp,
                     const size_t              idx_ovl_sd,
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

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_xy = pbuffer.data(idx_ovl_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_ovl_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up 0-6 components of targeted buffer : PD

    auto ts_x_xx = pbuffer.data(idx_ovl_pd);

    auto ts_x_xy = pbuffer.data(idx_ovl_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_ovl_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_ovl_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_ovl_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_ovl_pd + 5);

#pragma omp simd aligned(pa_x,        \
                             ts_0_x,  \
                             ts_0_xx, \
                             ts_0_xy, \
                             ts_0_xz, \
                             ts_0_y,  \
                             ts_0_yy, \
                             ts_0_yz, \
                             ts_0_z,  \
                             ts_0_zz, \
                             ts_x_xx, \
                             ts_x_xy, \
                             ts_x_xz, \
                             ts_x_yy, \
                             ts_x_yz, \
                             ts_x_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_x_xx[i] = 2.0 * ts_0_x[i] * fe_0 + ts_0_xx[i] * pa_x[i];

        ts_x_xy[i] = ts_0_y[i] * fe_0 + ts_0_xy[i] * pa_x[i];

        ts_x_xz[i] = ts_0_z[i] * fe_0 + ts_0_xz[i] * pa_x[i];

        ts_x_yy[i] = ts_0_yy[i] * pa_x[i];

        ts_x_yz[i] = ts_0_yz[i] * pa_x[i];

        ts_x_zz[i] = ts_0_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto ts_y_xx = pbuffer.data(idx_ovl_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_ovl_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_ovl_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_ovl_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_ovl_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_ovl_pd + 11);

#pragma omp simd aligned(pa_y,        \
                             ts_0_x,  \
                             ts_0_xx, \
                             ts_0_xy, \
                             ts_0_xz, \
                             ts_0_y,  \
                             ts_0_yy, \
                             ts_0_yz, \
                             ts_0_z,  \
                             ts_0_zz, \
                             ts_y_xx, \
                             ts_y_xy, \
                             ts_y_xz, \
                             ts_y_yy, \
                             ts_y_yz, \
                             ts_y_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_y_xx[i] = ts_0_xx[i] * pa_y[i];

        ts_y_xy[i] = ts_0_x[i] * fe_0 + ts_0_xy[i] * pa_y[i];

        ts_y_xz[i] = ts_0_xz[i] * pa_y[i];

        ts_y_yy[i] = 2.0 * ts_0_y[i] * fe_0 + ts_0_yy[i] * pa_y[i];

        ts_y_yz[i] = ts_0_z[i] * fe_0 + ts_0_yz[i] * pa_y[i];

        ts_y_zz[i] = ts_0_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto ts_z_xx = pbuffer.data(idx_ovl_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_ovl_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_ovl_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_ovl_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_ovl_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_ovl_pd + 17);

#pragma omp simd aligned(pa_z,        \
                             ts_0_x,  \
                             ts_0_xx, \
                             ts_0_xy, \
                             ts_0_xz, \
                             ts_0_y,  \
                             ts_0_yy, \
                             ts_0_yz, \
                             ts_0_z,  \
                             ts_0_zz, \
                             ts_z_xx, \
                             ts_z_xy, \
                             ts_z_xz, \
                             ts_z_yy, \
                             ts_z_yz, \
                             ts_z_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_z_xx[i] = ts_0_xx[i] * pa_z[i];

        ts_z_xy[i] = ts_0_xy[i] * pa_z[i];

        ts_z_xz[i] = ts_0_x[i] * fe_0 + ts_0_xz[i] * pa_z[i];

        ts_z_yy[i] = ts_0_yy[i] * pa_z[i];

        ts_z_yz[i] = ts_0_y[i] * fe_0 + ts_0_yz[i] * pa_z[i];

        ts_z_zz[i] = 2.0 * ts_0_z[i] * fe_0 + ts_0_zz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
