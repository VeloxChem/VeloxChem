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

#include "ElectronRepulsionGeom0010ContrRecXXPS.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxps(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxps,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxss,
                                            const size_t idx_geom_10_xxss,
                                            const size_t idx_geom_10_xxsp,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSS

            const auto ss_off = idx_xxss + (i * bcomps + j) * 1;

            auto g_0_0 = pbuffer.data(ss_off + 0);

            /// Set up components of auxilary buffer : SSSS

            const auto ss_geom_10_off = idx_geom_10_xxss + (i * bcomps + j) * 1;

            auto g_x_0_0_0 = cbuffer.data(ss_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_y_0_0_0 = cbuffer.data(ss_geom_10_off + 1 * acomps * bcomps + 0);

            auto g_z_0_0_0 = cbuffer.data(ss_geom_10_off + 2 * acomps * bcomps + 0);

            /// Set up components of auxilary buffer : SSSP

            const auto sp_geom_10_off = idx_geom_10_xxsp + (i * bcomps + j) * 3;

            auto g_x_0_0_x = cbuffer.data(sp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_y = cbuffer.data(sp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_z = cbuffer.data(sp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_y_0_0_x = cbuffer.data(sp_geom_10_off + 3 * acomps * bcomps + 0);

            auto g_y_0_0_y = cbuffer.data(sp_geom_10_off + 3 * acomps * bcomps + 1);

            auto g_y_0_0_z = cbuffer.data(sp_geom_10_off + 3 * acomps * bcomps + 2);

            auto g_z_0_0_x = cbuffer.data(sp_geom_10_off + 6 * acomps * bcomps + 0);

            auto g_z_0_0_y = cbuffer.data(sp_geom_10_off + 6 * acomps * bcomps + 1);

            auto g_z_0_0_z = cbuffer.data(sp_geom_10_off + 6 * acomps * bcomps + 2);

            /// set up bra offset for contr_buffer_xxps

            const auto ps_geom_10_off = idx_geom_10_xxps + (i * bcomps + j) * 3;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0 = cbuffer.data(ps_geom_10_off + 0 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_0_0, g_x_0_0_0, g_x_0_0_x, g_x_0_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0[k] = -g_0_0[k] - g_x_0_0_0[k] * cd_x[k] + g_x_0_0_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0 = cbuffer.data(ps_geom_10_off + 0 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_x_0_0_0, g_x_0_0_y, g_x_0_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0[k] = -g_x_0_0_0[k] * cd_y[k] + g_x_0_0_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0 = cbuffer.data(ps_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_x_0_0_0, g_x_0_0_z, g_x_0_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0[k] = -g_x_0_0_0[k] * cd_z[k] + g_x_0_0_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0 = cbuffer.data(ps_geom_10_off + 3 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_y_0_0_0, g_y_0_0_x, g_y_0_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0[k] = -g_y_0_0_0[k] * cd_x[k] + g_y_0_0_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0 = cbuffer.data(ps_geom_10_off + 3 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_0_0, g_y_0_0_0, g_y_0_0_y, g_y_0_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0[k] = -g_0_0[k] - g_y_0_0_0[k] * cd_y[k] + g_y_0_0_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0 = cbuffer.data(ps_geom_10_off + 3 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_y_0_0_0, g_y_0_0_z, g_y_0_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0[k] = -g_y_0_0_0[k] * cd_z[k] + g_y_0_0_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0 = cbuffer.data(ps_geom_10_off + 6 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_z_0_0_0, g_z_0_0_x, g_z_0_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0[k] = -g_z_0_0_0[k] * cd_x[k] + g_z_0_0_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0 = cbuffer.data(ps_geom_10_off + 6 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_z_0_0_0, g_z_0_0_y, g_z_0_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0[k] = -g_z_0_0_0[k] * cd_y[k] + g_z_0_0_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0 = cbuffer.data(ps_geom_10_off + 6 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_0_0, g_z_0_0_0, g_z_0_0_z, g_z_0_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0[k] = -g_0_0[k] - g_z_0_0_0[k] * cd_z[k] + g_z_0_0_z[k];
            }
        }
    }
}

} // erirec namespace

