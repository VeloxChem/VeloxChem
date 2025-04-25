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

#include "OverlapPrimRecGS.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_gs(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_gs,
                     const size_t              idx_ovl_ds,
                     const size_t              idx_ovl_fs,
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

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ovl_ds);

    auto ts_yy_0 = pbuffer.data(idx_ovl_ds + 3);

    auto ts_zz_0 = pbuffer.data(idx_ovl_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_ovl_fs);

    auto ts_xxz_0 = pbuffer.data(idx_ovl_fs + 2);

    auto ts_xyy_0 = pbuffer.data(idx_ovl_fs + 3);

    auto ts_xzz_0 = pbuffer.data(idx_ovl_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_ovl_fs + 6);

    auto ts_yyz_0 = pbuffer.data(idx_ovl_fs + 7);

    auto ts_yzz_0 = pbuffer.data(idx_ovl_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_ovl_fs + 9);

    // Set up components of targeted buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_ovl_gs);

    auto ts_xxxy_0 = pbuffer.data(idx_ovl_gs + 1);

    auto ts_xxxz_0 = pbuffer.data(idx_ovl_gs + 2);

    auto ts_xxyy_0 = pbuffer.data(idx_ovl_gs + 3);

    auto ts_xxyz_0 = pbuffer.data(idx_ovl_gs + 4);

    auto ts_xxzz_0 = pbuffer.data(idx_ovl_gs + 5);

    auto ts_xyyy_0 = pbuffer.data(idx_ovl_gs + 6);

    auto ts_xyyz_0 = pbuffer.data(idx_ovl_gs + 7);

    auto ts_xyzz_0 = pbuffer.data(idx_ovl_gs + 8);

    auto ts_xzzz_0 = pbuffer.data(idx_ovl_gs + 9);

    auto ts_yyyy_0 = pbuffer.data(idx_ovl_gs + 10);

    auto ts_yyyz_0 = pbuffer.data(idx_ovl_gs + 11);

    auto ts_yyzz_0 = pbuffer.data(idx_ovl_gs + 12);

    auto ts_yzzz_0 = pbuffer.data(idx_ovl_gs + 13);

    auto ts_zzzz_0 = pbuffer.data(idx_ovl_gs + 14);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             pa_z,      \
                             ts_xx_0,   \
                             ts_xxx_0,  \
                             ts_xxxx_0, \
                             ts_xxxy_0, \
                             ts_xxxz_0, \
                             ts_xxyy_0, \
                             ts_xxyz_0, \
                             ts_xxz_0,  \
                             ts_xxzz_0, \
                             ts_xyy_0,  \
                             ts_xyyy_0, \
                             ts_xyyz_0, \
                             ts_xyzz_0, \
                             ts_xzz_0,  \
                             ts_xzzz_0, \
                             ts_yy_0,   \
                             ts_yyy_0,  \
                             ts_yyyy_0, \
                             ts_yyyz_0, \
                             ts_yyz_0,  \
                             ts_yyzz_0, \
                             ts_yzz_0,  \
                             ts_yzzz_0, \
                             ts_zz_0,   \
                             ts_zzz_0,  \
                             ts_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxx_0[i] = 3.0 * ts_xx_0[i] * fe_0 + ts_xxx_0[i] * pa_x[i];

        ts_xxxy_0[i] = ts_xxx_0[i] * pa_y[i];

        ts_xxxz_0[i] = ts_xxx_0[i] * pa_z[i];

        ts_xxyy_0[i] = ts_yy_0[i] * fe_0 + ts_xyy_0[i] * pa_x[i];

        ts_xxyz_0[i] = ts_xxz_0[i] * pa_y[i];

        ts_xxzz_0[i] = ts_zz_0[i] * fe_0 + ts_xzz_0[i] * pa_x[i];

        ts_xyyy_0[i] = ts_yyy_0[i] * pa_x[i];

        ts_xyyz_0[i] = ts_yyz_0[i] * pa_x[i];

        ts_xyzz_0[i] = ts_yzz_0[i] * pa_x[i];

        ts_xzzz_0[i] = ts_zzz_0[i] * pa_x[i];

        ts_yyyy_0[i] = 3.0 * ts_yy_0[i] * fe_0 + ts_yyy_0[i] * pa_y[i];

        ts_yyyz_0[i] = ts_yyy_0[i] * pa_z[i];

        ts_yyzz_0[i] = ts_zz_0[i] * fe_0 + ts_yzz_0[i] * pa_y[i];

        ts_yzzz_0[i] = ts_zzz_0[i] * pa_y[i];

        ts_zzzz_0[i] = 3.0 * ts_zz_0[i] * fe_0 + ts_zzz_0[i] * pa_z[i];
    }
}

}  // namespace ovlrec
