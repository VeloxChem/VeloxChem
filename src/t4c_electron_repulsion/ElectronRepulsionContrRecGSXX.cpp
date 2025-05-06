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

#include "ElectronRepulsionContrRecGSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_gsxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_gsxx,
                                     const size_t idx_fsxx,
                                     const size_t idx_fpxx,
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
            /// Set up components of auxilary buffer : FSSS

            const auto fs_off = idx_fsxx + i * dcomps + j;

            auto g_xxx_0 = cbuffer.data(fs_off + 0 * ccomps * dcomps);

            auto g_xxy_0 = cbuffer.data(fs_off + 1 * ccomps * dcomps);

            auto g_xxz_0 = cbuffer.data(fs_off + 2 * ccomps * dcomps);

            auto g_xyy_0 = cbuffer.data(fs_off + 3 * ccomps * dcomps);

            auto g_xyz_0 = cbuffer.data(fs_off + 4 * ccomps * dcomps);

            auto g_xzz_0 = cbuffer.data(fs_off + 5 * ccomps * dcomps);

            auto g_yyy_0 = cbuffer.data(fs_off + 6 * ccomps * dcomps);

            auto g_yyz_0 = cbuffer.data(fs_off + 7 * ccomps * dcomps);

            auto g_yzz_0 = cbuffer.data(fs_off + 8 * ccomps * dcomps);

            auto g_zzz_0 = cbuffer.data(fs_off + 9 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FPSS

            const auto fp_off = idx_fpxx + i * dcomps + j;

            auto g_xxx_x = cbuffer.data(fp_off + 0 * ccomps * dcomps);

            auto g_xxy_x = cbuffer.data(fp_off + 3 * ccomps * dcomps);

            auto g_xxz_x = cbuffer.data(fp_off + 6 * ccomps * dcomps);

            auto g_xyy_x = cbuffer.data(fp_off + 9 * ccomps * dcomps);

            auto g_xyz_x = cbuffer.data(fp_off + 12 * ccomps * dcomps);

            auto g_xzz_x = cbuffer.data(fp_off + 15 * ccomps * dcomps);

            auto g_yyy_x = cbuffer.data(fp_off + 18 * ccomps * dcomps);

            auto g_yyy_y = cbuffer.data(fp_off + 19 * ccomps * dcomps);

            auto g_yyz_x = cbuffer.data(fp_off + 21 * ccomps * dcomps);

            auto g_yyz_y = cbuffer.data(fp_off + 22 * ccomps * dcomps);

            auto g_yzz_x = cbuffer.data(fp_off + 24 * ccomps * dcomps);

            auto g_yzz_y = cbuffer.data(fp_off + 25 * ccomps * dcomps);

            auto g_zzz_x = cbuffer.data(fp_off + 27 * ccomps * dcomps);

            auto g_zzz_y = cbuffer.data(fp_off + 28 * ccomps * dcomps);

            auto g_zzz_z = cbuffer.data(fp_off + 29 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gsxx

            const auto gs_off = idx_gsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_xxxx_0 = cbuffer.data(gs_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxx_0, g_xxx_x, g_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxxx_0[k] = -g_xxx_0[k] * ab_x + g_xxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_xxxy_0 = cbuffer.data(gs_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxy_0, g_xxy_0, g_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxxy_0[k] = -g_xxy_0[k] * ab_x + g_xxy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_xxxz_0 = cbuffer.data(gs_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxz_0, g_xxz_0, g_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxxz_0[k] = -g_xxz_0[k] * ab_x + g_xxz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_xxyy_0 = cbuffer.data(gs_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxyy_0, g_xyy_0, g_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxyy_0[k] = -g_xyy_0[k] * ab_x + g_xyy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_xxyz_0 = cbuffer.data(gs_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxyz_0, g_xyz_0, g_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxyz_0[k] = -g_xyz_0[k] * ab_x + g_xyz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_xxzz_0 = cbuffer.data(gs_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxzz_0, g_xzz_0, g_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxzz_0[k] = -g_xzz_0[k] * ab_x + g_xzz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_xyyy_0 = cbuffer.data(gs_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyyy_0, g_yyy_0, g_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyyy_0[k] = -g_yyy_0[k] * ab_x + g_yyy_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_xyyz_0 = cbuffer.data(gs_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyyz_0, g_yyz_0, g_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyyz_0[k] = -g_yyz_0[k] * ab_x + g_yyz_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_xyzz_0 = cbuffer.data(gs_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyzz_0, g_yzz_0, g_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyzz_0[k] = -g_yzz_0[k] * ab_x + g_yzz_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_xzzz_0 = cbuffer.data(gs_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_xzzz_0, g_zzz_0, g_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xzzz_0[k] = -g_zzz_0[k] * ab_x + g_zzz_x[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_yyyy_0 = cbuffer.data(gs_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyy_0, g_yyy_y, g_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyyy_0[k] = -g_yyy_0[k] * ab_y + g_yyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_yyyz_0 = cbuffer.data(gs_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyyz_0, g_yyz_0, g_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyyz_0[k] = -g_yyz_0[k] * ab_y + g_yyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_yyzz_0 = cbuffer.data(gs_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyzz_0, g_yzz_0, g_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyzz_0[k] = -g_yzz_0[k] * ab_y + g_yzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_yzzz_0 = cbuffer.data(gs_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_yzzz_0, g_zzz_0, g_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yzzz_0[k] = -g_zzz_0[k] * ab_y + g_zzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_zzzz_0 = cbuffer.data(gs_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_zzz_0, g_zzz_z, g_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zzzz_0[k] = -g_zzz_0[k] * ab_z + g_zzz_z[k];
            }
        }
    }
}

} // erirec namespace

