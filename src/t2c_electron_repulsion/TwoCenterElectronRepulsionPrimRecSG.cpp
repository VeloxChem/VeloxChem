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

#include "TwoCenterElectronRepulsionPrimRecSG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_sg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_sg,
                                const size_t idx_eri_0_sd,
                                const size_t idx_eri_1_sd,
                                const size_t idx_eri_1_sf,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpb,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SD

    auto g_0_xx_0 = pbuffer.data(idx_eri_0_sd);

    auto g_0_yy_0 = pbuffer.data(idx_eri_0_sd + 3);

    auto g_0_zz_0 = pbuffer.data(idx_eri_0_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto g_0_xx_1 = pbuffer.data(idx_eri_1_sd);

    auto g_0_yy_1 = pbuffer.data(idx_eri_1_sd + 3);

    auto g_0_zz_1 = pbuffer.data(idx_eri_1_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto g_0_xxx_1 = pbuffer.data(idx_eri_1_sf);

    auto g_0_xxz_1 = pbuffer.data(idx_eri_1_sf + 2);

    auto g_0_xyy_1 = pbuffer.data(idx_eri_1_sf + 3);

    auto g_0_xzz_1 = pbuffer.data(idx_eri_1_sf + 5);

    auto g_0_yyy_1 = pbuffer.data(idx_eri_1_sf + 6);

    auto g_0_yyz_1 = pbuffer.data(idx_eri_1_sf + 7);

    auto g_0_yzz_1 = pbuffer.data(idx_eri_1_sf + 8);

    auto g_0_zzz_1 = pbuffer.data(idx_eri_1_sf + 9);

    // Set up components of targeted buffer : SG

    auto g_0_xxxx_0 = pbuffer.data(idx_eri_0_sg);

    auto g_0_xxxy_0 = pbuffer.data(idx_eri_0_sg + 1);

    auto g_0_xxxz_0 = pbuffer.data(idx_eri_0_sg + 2);

    auto g_0_xxyy_0 = pbuffer.data(idx_eri_0_sg + 3);

    auto g_0_xxyz_0 = pbuffer.data(idx_eri_0_sg + 4);

    auto g_0_xxzz_0 = pbuffer.data(idx_eri_0_sg + 5);

    auto g_0_xyyy_0 = pbuffer.data(idx_eri_0_sg + 6);

    auto g_0_xyyz_0 = pbuffer.data(idx_eri_0_sg + 7);

    auto g_0_xyzz_0 = pbuffer.data(idx_eri_0_sg + 8);

    auto g_0_xzzz_0 = pbuffer.data(idx_eri_0_sg + 9);

    auto g_0_yyyy_0 = pbuffer.data(idx_eri_0_sg + 10);

    auto g_0_yyyz_0 = pbuffer.data(idx_eri_0_sg + 11);

    auto g_0_yyzz_0 = pbuffer.data(idx_eri_0_sg + 12);

    auto g_0_yzzz_0 = pbuffer.data(idx_eri_0_sg + 13);

    auto g_0_zzzz_0 = pbuffer.data(idx_eri_0_sg + 14);

    #pragma omp simd aligned(g_0_xx_0, g_0_xx_1, g_0_xxx_1, g_0_xxxx_0, g_0_xxxy_0, g_0_xxxz_0, g_0_xxyy_0, g_0_xxyz_0, g_0_xxz_1, g_0_xxzz_0, g_0_xyy_1, g_0_xyyy_0, g_0_xyyz_0, g_0_xyzz_0, g_0_xzz_1, g_0_xzzz_0, g_0_yy_0, g_0_yy_1, g_0_yyy_1, g_0_yyyy_0, g_0_yyyz_0, g_0_yyz_1, g_0_yyzz_0, g_0_yzz_1, g_0_yzzz_0, g_0_zz_0, g_0_zz_1, g_0_zzz_1, g_0_zzzz_0, pb_x, pb_y, pb_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fke_0 = 0.5 / b_exps[i];

        const double fz_ke_0 = a_exp * fke_0 / (a_exp + b_exps[i]);

        g_0_xxxx_0[i] = 3.0 * g_0_xx_0[i] * fke_0 - 3.0 * g_0_xx_1[i] * fz_ke_0 + g_0_xxx_1[i] * pb_x[i];

        g_0_xxxy_0[i] = g_0_xxx_1[i] * pb_y[i];

        g_0_xxxz_0[i] = g_0_xxx_1[i] * pb_z[i];

        g_0_xxyy_0[i] = g_0_yy_0[i] * fke_0 - g_0_yy_1[i] * fz_ke_0 + g_0_xyy_1[i] * pb_x[i];

        g_0_xxyz_0[i] = g_0_xxz_1[i] * pb_y[i];

        g_0_xxzz_0[i] = g_0_zz_0[i] * fke_0 - g_0_zz_1[i] * fz_ke_0 + g_0_xzz_1[i] * pb_x[i];

        g_0_xyyy_0[i] = g_0_yyy_1[i] * pb_x[i];

        g_0_xyyz_0[i] = g_0_yyz_1[i] * pb_x[i];

        g_0_xyzz_0[i] = g_0_yzz_1[i] * pb_x[i];

        g_0_xzzz_0[i] = g_0_zzz_1[i] * pb_x[i];

        g_0_yyyy_0[i] = 3.0 * g_0_yy_0[i] * fke_0 - 3.0 * g_0_yy_1[i] * fz_ke_0 + g_0_yyy_1[i] * pb_y[i];

        g_0_yyyz_0[i] = g_0_yyy_1[i] * pb_z[i];

        g_0_yyzz_0[i] = g_0_zz_0[i] * fke_0 - g_0_zz_1[i] * fz_ke_0 + g_0_yzz_1[i] * pb_y[i];

        g_0_yzzz_0[i] = g_0_zzz_1[i] * pb_y[i];

        g_0_zzzz_0[i] = 3.0 * g_0_zz_0[i] * fke_0 - 3.0 * g_0_zz_1[i] * fz_ke_0 + g_0_zzz_1[i] * pb_z[i];
    }
}

} // t2ceri namespace

