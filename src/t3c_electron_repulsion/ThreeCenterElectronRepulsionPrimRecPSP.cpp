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

#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psp,
                                 size_t idx_eri_1_sss,
                                 size_t idx_eri_1_ssp,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : SSS

    auto g_0_0_0_1 = pbuffer.data(idx_eri_1_sss);

    /// Set up components of auxilary buffer : SSP

    auto g_0_0_x_1 = pbuffer.data(idx_eri_1_ssp);

    auto g_0_0_y_1 = pbuffer.data(idx_eri_1_ssp + 1);

    auto g_0_0_z_1 = pbuffer.data(idx_eri_1_ssp + 2);

    /// Set up 0-3 components of targeted buffer : PSP

    auto g_x_0_x_0 = pbuffer.data(idx_eri_0_psp);

    auto g_x_0_y_0 = pbuffer.data(idx_eri_0_psp + 1);

    auto g_x_0_z_0 = pbuffer.data(idx_eri_0_psp + 2);

    #pragma omp simd aligned(g_0_0_0_1, g_0_0_x_1, g_0_0_y_1, g_0_0_z_1, g_x_0_x_0, g_x_0_y_0, g_x_0_z_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_x_0[i] = g_0_0_0_1[i] * fi_acd_0 + g_0_0_x_1[i] * wa_x[i];

        g_x_0_y_0[i] = g_0_0_y_1[i] * wa_x[i];

        g_x_0_z_0[i] = g_0_0_z_1[i] * wa_x[i];
    }

    /// Set up 3-6 components of targeted buffer : PSP

    auto g_y_0_x_0 = pbuffer.data(idx_eri_0_psp + 3);

    auto g_y_0_y_0 = pbuffer.data(idx_eri_0_psp + 4);

    auto g_y_0_z_0 = pbuffer.data(idx_eri_0_psp + 5);

    #pragma omp simd aligned(g_0_0_0_1, g_0_0_x_1, g_0_0_y_1, g_0_0_z_1, g_y_0_x_0, g_y_0_y_0, g_y_0_z_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_x_0[i] = g_0_0_x_1[i] * wa_y[i];

        g_y_0_y_0[i] = g_0_0_0_1[i] * fi_acd_0 + g_0_0_y_1[i] * wa_y[i];

        g_y_0_z_0[i] = g_0_0_z_1[i] * wa_y[i];
    }

    /// Set up 6-9 components of targeted buffer : PSP

    auto g_z_0_x_0 = pbuffer.data(idx_eri_0_psp + 6);

    auto g_z_0_y_0 = pbuffer.data(idx_eri_0_psp + 7);

    auto g_z_0_z_0 = pbuffer.data(idx_eri_0_psp + 8);

    #pragma omp simd aligned(g_0_0_0_1, g_0_0_x_1, g_0_0_y_1, g_0_0_z_1, g_z_0_x_0, g_z_0_y_0, g_z_0_z_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_x_0[i] = g_0_0_x_1[i] * wa_z[i];

        g_z_0_y_0[i] = g_0_0_y_1[i] * wa_z[i];

        g_z_0_z_0[i] = g_0_0_0_1[i] * fi_acd_0 + g_0_0_z_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

