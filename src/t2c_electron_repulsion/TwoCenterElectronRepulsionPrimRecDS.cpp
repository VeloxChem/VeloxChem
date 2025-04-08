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

#include "TwoCenterElectronRepulsionPrimRecDS.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ds(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ds,
                                const size_t idx_eri_0_ss,
                                const size_t idx_eri_1_ss,
                                const size_t idx_eri_1_ps,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SS

    auto g_0_0_0 = pbuffer.data(idx_eri_0_ss);

    // Set up components of auxiliary buffer : SS

    auto g_0_0_1 = pbuffer.data(idx_eri_1_ss);

    // Set up components of auxiliary buffer : PS

    auto g_x_0_1 = pbuffer.data(idx_eri_1_ps);

    auto g_y_0_1 = pbuffer.data(idx_eri_1_ps + 1);

    auto g_z_0_1 = pbuffer.data(idx_eri_1_ps + 2);

    // Set up components of targeted buffer : DS

    auto g_xx_0_0 = pbuffer.data(idx_eri_0_ds);

    auto g_xy_0_0 = pbuffer.data(idx_eri_0_ds + 1);

    auto g_xz_0_0 = pbuffer.data(idx_eri_0_ds + 2);

    auto g_yy_0_0 = pbuffer.data(idx_eri_0_ds + 3);

    auto g_yz_0_0 = pbuffer.data(idx_eri_0_ds + 4);

    auto g_zz_0_0 = pbuffer.data(idx_eri_0_ds + 5);

    #pragma omp simd aligned(g_0_0_0, g_0_0_1, g_x_0_1, g_xx_0_0, g_xy_0_0, g_xz_0_0, g_y_0_1, g_yy_0_0, g_yz_0_0, g_z_0_1, g_zz_0_0, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xx_0_0[i] = g_0_0_0[i] * fbe_0 - g_0_0_1[i] * fz_be_0 + g_x_0_1[i] * pa_x[i];

        g_xy_0_0[i] = g_y_0_1[i] * pa_x[i];

        g_xz_0_0[i] = g_z_0_1[i] * pa_x[i];

        g_yy_0_0[i] = g_0_0_0[i] * fbe_0 - g_0_0_1[i] * fz_be_0 + g_y_0_1[i] * pa_y[i];

        g_yz_0_0[i] = g_z_0_1[i] * pa_y[i];

        g_zz_0_0[i] = g_0_0_0[i] * fbe_0 - g_0_0_1[i] * fz_be_0 + g_z_0_1[i] * pa_z[i];
    }
}

} // t2ceri namespace

