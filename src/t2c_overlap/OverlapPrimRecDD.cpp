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

#include "OverlapPrimRecDD.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_dd(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_dd,
                     const size_t              idx_ovl_sd,
                     const size_t              idx_ovl_pp,
                     const size_t              idx_ovl_pd,
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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_xy = pbuffer.data(idx_ovl_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_ovl_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_ovl_pp);

    auto ts_x_y = pbuffer.data(idx_ovl_pp + 1);

    auto ts_x_z = pbuffer.data(idx_ovl_pp + 2);

    auto ts_y_x = pbuffer.data(idx_ovl_pp + 3);

    auto ts_y_y = pbuffer.data(idx_ovl_pp + 4);

    auto ts_y_z = pbuffer.data(idx_ovl_pp + 5);

    auto ts_z_x = pbuffer.data(idx_ovl_pp + 6);

    auto ts_z_y = pbuffer.data(idx_ovl_pp + 7);

    auto ts_z_z = pbuffer.data(idx_ovl_pp + 8);

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_ovl_pd);

    auto ts_x_xy = pbuffer.data(idx_ovl_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_ovl_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_ovl_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_ovl_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_ovl_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_ovl_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_ovl_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_ovl_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_ovl_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_ovl_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_ovl_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_ovl_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_ovl_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_ovl_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_ovl_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_ovl_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_ovl_pd + 17);

    // Set up 0-6 components of targeted buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_ovl_dd);

    auto ts_xx_xy = pbuffer.data(idx_ovl_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_ovl_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_ovl_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_ovl_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_ovl_dd + 5);

#pragma omp simd aligned(pa_x,         \
                             ts_0_xx,  \
                             ts_0_xy,  \
                             ts_0_xz,  \
                             ts_0_yy,  \
                             ts_0_yz,  \
                             ts_0_zz,  \
                             ts_x_x,   \
                             ts_x_xx,  \
                             ts_x_xy,  \
                             ts_x_xz,  \
                             ts_x_y,   \
                             ts_x_yy,  \
                             ts_x_yz,  \
                             ts_x_z,   \
                             ts_x_zz,  \
                             ts_xx_xx, \
                             ts_xx_xy, \
                             ts_xx_xz, \
                             ts_xx_yy, \
                             ts_xx_yz, \
                             ts_xx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xx_xx[i] = ts_0_xx[i] * fe_0 + 2.0 * ts_x_x[i] * fe_0 + ts_x_xx[i] * pa_x[i];

        ts_xx_xy[i] = ts_0_xy[i] * fe_0 + ts_x_y[i] * fe_0 + ts_x_xy[i] * pa_x[i];

        ts_xx_xz[i] = ts_0_xz[i] * fe_0 + ts_x_z[i] * fe_0 + ts_x_xz[i] * pa_x[i];

        ts_xx_yy[i] = ts_0_yy[i] * fe_0 + ts_x_yy[i] * pa_x[i];

        ts_xx_yz[i] = ts_0_yz[i] * fe_0 + ts_x_yz[i] * pa_x[i];

        ts_xx_zz[i] = ts_0_zz[i] * fe_0 + ts_x_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto ts_xy_xx = pbuffer.data(idx_ovl_dd + 6);

    auto ts_xy_xy = pbuffer.data(idx_ovl_dd + 7);

    auto ts_xy_xz = pbuffer.data(idx_ovl_dd + 8);

    auto ts_xy_yy = pbuffer.data(idx_ovl_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_ovl_dd + 10);

    auto ts_xy_zz = pbuffer.data(idx_ovl_dd + 11);

#pragma omp simd aligned(pa_x,         \
                             pa_y,     \
                             ts_x_xx,  \
                             ts_x_xz,  \
                             ts_xy_xx, \
                             ts_xy_xy, \
                             ts_xy_xz, \
                             ts_xy_yy, \
                             ts_xy_yz, \
                             ts_xy_zz, \
                             ts_y_xy,  \
                             ts_y_y,   \
                             ts_y_yy,  \
                             ts_y_yz,  \
                             ts_y_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xy_xx[i] = ts_x_xx[i] * pa_y[i];

        ts_xy_xy[i] = ts_y_y[i] * fe_0 + ts_y_xy[i] * pa_x[i];

        ts_xy_xz[i] = ts_x_xz[i] * pa_y[i];

        ts_xy_yy[i] = ts_y_yy[i] * pa_x[i];

        ts_xy_yz[i] = ts_y_yz[i] * pa_x[i];

        ts_xy_zz[i] = ts_y_zz[i] * pa_x[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto ts_xz_xx = pbuffer.data(idx_ovl_dd + 12);

    auto ts_xz_xy = pbuffer.data(idx_ovl_dd + 13);

    auto ts_xz_xz = pbuffer.data(idx_ovl_dd + 14);

    auto ts_xz_yy = pbuffer.data(idx_ovl_dd + 15);

    auto ts_xz_yz = pbuffer.data(idx_ovl_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_ovl_dd + 17);

#pragma omp simd aligned(pa_x,         \
                             pa_z,     \
                             ts_x_xx,  \
                             ts_x_xy,  \
                             ts_xz_xx, \
                             ts_xz_xy, \
                             ts_xz_xz, \
                             ts_xz_yy, \
                             ts_xz_yz, \
                             ts_xz_zz, \
                             ts_z_xz,  \
                             ts_z_yy,  \
                             ts_z_yz,  \
                             ts_z_z,   \
                             ts_z_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xz_xx[i] = ts_x_xx[i] * pa_z[i];

        ts_xz_xy[i] = ts_x_xy[i] * pa_z[i];

        ts_xz_xz[i] = ts_z_z[i] * fe_0 + ts_z_xz[i] * pa_x[i];

        ts_xz_yy[i] = ts_z_yy[i] * pa_x[i];

        ts_xz_yz[i] = ts_z_yz[i] * pa_x[i];

        ts_xz_zz[i] = ts_z_zz[i] * pa_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto ts_yy_xx = pbuffer.data(idx_ovl_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_ovl_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_ovl_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_ovl_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_ovl_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_ovl_dd + 23);

#pragma omp simd aligned(pa_y,         \
                             ts_0_xx,  \
                             ts_0_xy,  \
                             ts_0_xz,  \
                             ts_0_yy,  \
                             ts_0_yz,  \
                             ts_0_zz,  \
                             ts_y_x,   \
                             ts_y_xx,  \
                             ts_y_xy,  \
                             ts_y_xz,  \
                             ts_y_y,   \
                             ts_y_yy,  \
                             ts_y_yz,  \
                             ts_y_z,   \
                             ts_y_zz,  \
                             ts_yy_xx, \
                             ts_yy_xy, \
                             ts_yy_xz, \
                             ts_yy_yy, \
                             ts_yy_yz, \
                             ts_yy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yy_xx[i] = ts_0_xx[i] * fe_0 + ts_y_xx[i] * pa_y[i];

        ts_yy_xy[i] = ts_0_xy[i] * fe_0 + ts_y_x[i] * fe_0 + ts_y_xy[i] * pa_y[i];

        ts_yy_xz[i] = ts_0_xz[i] * fe_0 + ts_y_xz[i] * pa_y[i];

        ts_yy_yy[i] = ts_0_yy[i] * fe_0 + 2.0 * ts_y_y[i] * fe_0 + ts_y_yy[i] * pa_y[i];

        ts_yy_yz[i] = ts_0_yz[i] * fe_0 + ts_y_z[i] * fe_0 + ts_y_yz[i] * pa_y[i];

        ts_yy_zz[i] = ts_0_zz[i] * fe_0 + ts_y_zz[i] * pa_y[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto ts_yz_xx = pbuffer.data(idx_ovl_dd + 24);

    auto ts_yz_xy = pbuffer.data(idx_ovl_dd + 25);

    auto ts_yz_xz = pbuffer.data(idx_ovl_dd + 26);

    auto ts_yz_yy = pbuffer.data(idx_ovl_dd + 27);

    auto ts_yz_yz = pbuffer.data(idx_ovl_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_ovl_dd + 29);

#pragma omp simd aligned(pa_y,         \
                             pa_z,     \
                             ts_y_xy,  \
                             ts_y_yy,  \
                             ts_yz_xx, \
                             ts_yz_xy, \
                             ts_yz_xz, \
                             ts_yz_yy, \
                             ts_yz_yz, \
                             ts_yz_zz, \
                             ts_z_xx,  \
                             ts_z_xz,  \
                             ts_z_yz,  \
                             ts_z_z,   \
                             ts_z_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yz_xx[i] = ts_z_xx[i] * pa_y[i];

        ts_yz_xy[i] = ts_y_xy[i] * pa_z[i];

        ts_yz_xz[i] = ts_z_xz[i] * pa_y[i];

        ts_yz_yy[i] = ts_y_yy[i] * pa_z[i];

        ts_yz_yz[i] = ts_z_z[i] * fe_0 + ts_z_yz[i] * pa_y[i];

        ts_yz_zz[i] = ts_z_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto ts_zz_xx = pbuffer.data(idx_ovl_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_ovl_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_ovl_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_ovl_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_ovl_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_ovl_dd + 35);

#pragma omp simd aligned(pa_z,         \
                             ts_0_xx,  \
                             ts_0_xy,  \
                             ts_0_xz,  \
                             ts_0_yy,  \
                             ts_0_yz,  \
                             ts_0_zz,  \
                             ts_z_x,   \
                             ts_z_xx,  \
                             ts_z_xy,  \
                             ts_z_xz,  \
                             ts_z_y,   \
                             ts_z_yy,  \
                             ts_z_yz,  \
                             ts_z_z,   \
                             ts_z_zz,  \
                             ts_zz_xx, \
                             ts_zz_xy, \
                             ts_zz_xz, \
                             ts_zz_yy, \
                             ts_zz_yz, \
                             ts_zz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zz_xx[i] = ts_0_xx[i] * fe_0 + ts_z_xx[i] * pa_z[i];

        ts_zz_xy[i] = ts_0_xy[i] * fe_0 + ts_z_xy[i] * pa_z[i];

        ts_zz_xz[i] = ts_0_xz[i] * fe_0 + ts_z_x[i] * fe_0 + ts_z_xz[i] * pa_z[i];

        ts_zz_yy[i] = ts_0_yy[i] * fe_0 + ts_z_yy[i] * pa_z[i];

        ts_zz_yz[i] = ts_0_yz[i] * fe_0 + ts_z_y[i] * fe_0 + ts_z_yz[i] * pa_z[i];

        ts_zz_zz[i] = ts_0_zz[i] * fe_0 + 2.0 * ts_z_z[i] * fe_0 + ts_z_zz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
