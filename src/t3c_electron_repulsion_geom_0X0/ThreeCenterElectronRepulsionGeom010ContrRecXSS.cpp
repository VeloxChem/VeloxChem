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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXSS.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xss(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xss,
                                        const size_t idx_xss,
                                        const size_t idx_xsp,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SSS

        const auto ss_off = idx_xss + i * 1;

        auto g_0_0 = cbuffer.data(ss_off + 0);

        /// Set up components of auxilary buffer : SSP

        const auto sp_off = idx_xsp + i * 3;

        auto g_0_x = cbuffer.data(sp_off + 0);

        auto g_0_y = cbuffer.data(sp_off + 1);

        auto g_0_z = cbuffer.data(sp_off + 2);

        /// set up bra offset for contr_buffer_xxss

        const auto ss_geom_10_off = idx_geom_10_xss + i * 1;

        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_x_0_0_0 = cbuffer.data(ss_geom_10_off + 0 * acomps  + 0);

        #pragma omp simd aligned(cd_x, g_0_0, g_0_x, g_x_0_0_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_0_0[k] = -g_0_0[k] * cd_x[k] + g_0_x[k];
        }
        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_y_0_0_0 = cbuffer.data(ss_geom_10_off + 1 * acomps  + 0);

        #pragma omp simd aligned(cd_y, g_0_0, g_0_y, g_y_0_0_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_0_0[k] = -g_0_0[k] * cd_y[k] + g_0_y[k];
        }
        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_z_0_0_0 = cbuffer.data(ss_geom_10_off + 2 * acomps  + 0);

        #pragma omp simd aligned(cd_z, g_0_0, g_0_z, g_z_0_0_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_0_0[k] = -g_0_0[k] * cd_z[k] + g_0_z[k];
        }
    }
}

} // t3ceri namespace

