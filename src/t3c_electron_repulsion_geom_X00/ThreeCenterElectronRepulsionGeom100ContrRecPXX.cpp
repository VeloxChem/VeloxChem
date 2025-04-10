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

#include "ThreeCenterElectronRepulsionGeom100ContrRecPXX.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_bra_geom1_electron_repulsion_pxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_pxx,
                                      const size_t idx_sxx,
                                      const size_t idx_dxx,
                                      const int c_angmom,
                                      const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_cartesian_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_cartesian_components(std::array<int, 1>{d_angmom,});

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SSS

            const auto s_off = idx_sxx + i * dcomps + j;

            auto g_0_0 = cbuffer.data(s_off + 0 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DSS

            const auto d_off = idx_dxx + i * dcomps + j;

            auto g_xx_0 = cbuffer.data(d_off + 0 * ccomps * dcomps);

            auto g_xy_0 = cbuffer.data(d_off + 1 * ccomps * dcomps);

            auto g_xz_0 = cbuffer.data(d_off + 2 * ccomps * dcomps);

            auto g_yy_0 = cbuffer.data(d_off + 3 * ccomps * dcomps);

            auto g_yz_0 = cbuffer.data(d_off + 4 * ccomps * dcomps);

            auto g_zz_0 = cbuffer.data(d_off + 5 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pxx

            const auto p_geom_100_off = idx_geom_100_pxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_x_0 = cbuffer.data(p_geom_100_off + 0 * ccomps * dcomps);

            auto g_x_0_0_y_0 = cbuffer.data(p_geom_100_off + 1 * ccomps * dcomps);

            auto g_x_0_0_z_0 = cbuffer.data(p_geom_100_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_x_0_0_x_0, g_x_0_0_y_0, g_x_0_0_z_0, g_xx_0, g_xy_0, g_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_x_0[k] = -g_0_0[k] + g_xx_0[k];

                g_x_0_0_y_0[k] = g_xy_0[k];

                g_x_0_0_z_0[k] = g_xz_0[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_x_0 = cbuffer.data(p_geom_100_off + 3 * ccomps * dcomps);

            auto g_y_0_0_y_0 = cbuffer.data(p_geom_100_off + 4 * ccomps * dcomps);

            auto g_y_0_0_z_0 = cbuffer.data(p_geom_100_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_xy_0, g_y_0_0_x_0, g_y_0_0_y_0, g_y_0_0_z_0, g_yy_0, g_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_x_0[k] = g_xy_0[k];

                g_y_0_0_y_0[k] = -g_0_0[k] + g_yy_0[k];

                g_y_0_0_z_0[k] = g_yz_0[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_x_0 = cbuffer.data(p_geom_100_off + 6 * ccomps * dcomps);

            auto g_z_0_0_y_0 = cbuffer.data(p_geom_100_off + 7 * ccomps * dcomps);

            auto g_z_0_0_z_0 = cbuffer.data(p_geom_100_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_xz_0, g_yz_0, g_z_0_0_x_0, g_z_0_0_y_0, g_z_0_0_z_0, g_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_x_0[k] = g_xz_0[k];

                g_z_0_0_y_0[k] = g_yz_0[k];

                g_z_0_0_z_0[k] = -g_0_0[k] + g_zz_0[k];
            }
        }
    }
}

} // t3ceri namespace

