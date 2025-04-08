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

#include "ElectronRepulsionGeom0100ContrRecPSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_psxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_psxx,
                                            const size_t idx_ssxx,
                                            const size_t idx_geom_01_ssxx,
                                            const size_t idx_geom_01_spxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSS

            const auto ss_off = idx_ssxx + i * dcomps + j;

            auto g_0_0 = cbuffer.data(ss_off + 0 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SSSS

            const auto ss_geom_01_off = idx_geom_01_ssxx + i * dcomps + j;

            auto g_0_x_0_0 = cbuffer.data(ss_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_y_0_0 = cbuffer.data(ss_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_z_0_0 = cbuffer.data(ss_geom_01_off + 2 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SPSS

            const auto sp_geom_01_off = idx_geom_01_spxx + i * dcomps + j;

            auto g_0_x_0_x = cbuffer.data(sp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_y = cbuffer.data(sp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_z = cbuffer.data(sp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_y_0_x = cbuffer.data(sp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_y_0_y = cbuffer.data(sp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_y_0_z = cbuffer.data(sp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_z_0_x = cbuffer.data(sp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_z_0_y = cbuffer.data(sp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_z_0_z = cbuffer.data(sp_geom_01_off + 8 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_psxx

            const auto ps_geom_01_off = idx_geom_01_psxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_0_x_x_0 = cbuffer.data(ps_geom_01_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_x_0_0, g_0_x_0_x, g_0_x_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_x_0[k] = g_0_0[k] - g_0_x_0_0[k] * ab_x + g_0_x_0_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_0_x_y_0 = cbuffer.data(ps_geom_01_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_0, g_0_x_0_y, g_0_x_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_y_0[k] = -g_0_x_0_0[k] * ab_y + g_0_x_0_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_z_0 = cbuffer.data(ps_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_0, g_0_x_0_z, g_0_x_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_z_0[k] = -g_0_x_0_0[k] * ab_z + g_0_x_0_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_0_y_x_0 = cbuffer.data(ps_geom_01_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_0, g_0_y_0_x, g_0_y_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_x_0[k] = -g_0_y_0_0[k] * ab_x + g_0_y_0_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_0_y_y_0 = cbuffer.data(ps_geom_01_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_y_0_0, g_0_y_0_y, g_0_y_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_y_0[k] = g_0_0[k] - g_0_y_0_0[k] * ab_y + g_0_y_0_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_0_y_z_0 = cbuffer.data(ps_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_0, g_0_y_0_z, g_0_y_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_z_0[k] = -g_0_y_0_0[k] * ab_z + g_0_y_0_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_0_z_x_0 = cbuffer.data(ps_geom_01_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_0, g_0_z_0_x, g_0_z_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_x_0[k] = -g_0_z_0_0[k] * ab_x + g_0_z_0_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_0_z_y_0 = cbuffer.data(ps_geom_01_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_0, g_0_z_0_y, g_0_z_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_y_0[k] = -g_0_z_0_0[k] * ab_y + g_0_z_0_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_0_z_z_0 = cbuffer.data(ps_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_z_0_0, g_0_z_0_z, g_0_z_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_z_0[k] = g_0_0[k] - g_0_z_0_0[k] * ab_z + g_0_z_0_z[k];
            }
        }
    }
}

} // erirec namespace

