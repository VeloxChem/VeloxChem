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

#include "TwoCenterElectronRepulsionPrimRecPS.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ps(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ps,
                                const size_t idx_eri_1_ss,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SS

    auto g_0_0_1 = pbuffer.data(idx_eri_1_ss);

    // Set up components of targeted buffer : PS

    auto g_x_0_0 = pbuffer.data(idx_eri_0_ps);

    auto g_y_0_0 = pbuffer.data(idx_eri_0_ps + 1);

    auto g_z_0_0 = pbuffer.data(idx_eri_0_ps + 2);

    #pragma omp simd aligned(g_0_0_1, g_x_0_0, g_y_0_0, g_z_0_0, pa_x, pa_y, pa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_x_0_0[i] = g_0_0_1[i] * pa_x[i];

        g_y_0_0[i] = g_0_0_1[i] * pa_y[i];

        g_z_0_0[i] = g_0_0_1[i] * pa_z[i];
    }
}

} // t2ceri namespace

