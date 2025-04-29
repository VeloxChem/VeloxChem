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

#include "ThreeCenterElectronRepulsionPrimRecSSP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssp,
                                 size_t idx_eri_0_sss,
                                 size_t idx_eri_1_sss,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(QD) distances

    auto qd_x = factors.data(idx_qd);

    auto qd_y = factors.data(idx_qd + 1);

    auto qd_z = factors.data(idx_qd + 2);

    // Set up R(WQ) distances

    auto wq_x = factors.data(idx_wq);

    auto wq_y = factors.data(idx_wq + 1);

    auto wq_z = factors.data(idx_wq + 2);

    /// Set up components of auxilary buffer : SSS

    auto g_0_0_0_0 = pbuffer.data(idx_eri_0_sss);

    /// Set up components of auxilary buffer : SSS

    auto g_0_0_0_1 = pbuffer.data(idx_eri_1_sss);

    /// Set up components of targeted buffer : SSP

    auto g_0_0_x_0 = pbuffer.data(idx_eri_0_ssp);

    auto g_0_0_y_0 = pbuffer.data(idx_eri_0_ssp + 1);

    auto g_0_0_z_0 = pbuffer.data(idx_eri_0_ssp + 2);

    #pragma omp simd aligned(g_0_0_0_0, g_0_0_0_1, g_0_0_x_0, g_0_0_y_0, g_0_0_z_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_0_x_0[i] = g_0_0_0_0[i] * qd_x[i] + g_0_0_0_1[i] * wq_x[i];

        g_0_0_y_0[i] = g_0_0_0_0[i] * qd_y[i] + g_0_0_0_1[i] * wq_y[i];

        g_0_0_z_0[i] = g_0_0_0_0[i] * qd_z[i] + g_0_0_0_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

